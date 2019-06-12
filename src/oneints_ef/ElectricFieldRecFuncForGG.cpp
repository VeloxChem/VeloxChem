//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectricFieldRecFuncForGG.hpp"

namespace efieldrecfunc {  // efieldrecfunc namespace

void
compElectricFieldForGG(CMemBlock2D<double>&       primBuffer,
                       const CRecursionMap&       recursionMap,
                       const CMemBlock2D<double>& osFactors,
                       const CMemBlock2D<double>& paDistances,
                       const CMemBlock2D<double>& pcDistances,
                       const CGtoBlock&           braGtoBlock,
                       const CGtoBlock&           ketGtoBlock,
                       const int32_t              iContrGto)
{
    efieldrecfunc::compElectricFieldForGG_0_49(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForGG_49_98(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForGG_98_147(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForGG_147_195(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForGG_195_243(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForGG_243_291(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForGG_291_339(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForGG_339_387(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForGG_387_435(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForGG_435_483(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForGG_483_531(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForGG_531_579(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForGG_579_627(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForGG_627_675(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compElectricFieldForGG_0_49(CMemBlock2D<double>&       primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CMemBlock2D<double>& pcDistances,
                            const CGtoBlock&           braGtoBlock,
                            const CGtoBlock&           ketGtoBlock,
                            const int32_t              iContrGto)
{
    // Batch of Integrals (0,49)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_4_4_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_2_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tex_xxx_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx);

            auto tey_xxx_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx);

            auto tez_xxx_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx);

            auto tex_xxx_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 1);

            auto tey_xxx_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 1);

            auto tez_xxx_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 1);

            auto tex_xxx_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 2);

            auto tey_xxx_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 2);

            auto tez_xxx_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 2);

            auto tex_xxx_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 3);

            auto tey_xxx_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 3);

            auto tez_xxx_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 3);

            auto tex_xxx_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 4);

            auto tey_xxx_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 4);

            auto tez_xxx_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 4);

            auto tex_xxx_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 5);

            auto tey_xxx_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 5);

            auto tez_xxx_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 5);

            auto tex_xxx_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 6);

            auto tey_xxx_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 6);

            auto tez_xxx_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 6);

            auto tex_xxx_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 7);

            auto tey_xxx_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 7);

            auto tez_xxx_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 7);

            auto tex_xxx_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 8);

            auto tey_xxx_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 8);

            auto tez_xxx_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 8);

            auto tex_xxx_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 9);

            auto tey_xxx_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 9);

            auto tez_xxx_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 9);

            auto tex_xxx_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 10);

            auto tey_xxx_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 10);

            auto tez_xxx_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 10);

            auto tex_xxx_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 11);

            auto tey_xxx_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 11);

            auto tez_xxx_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 11);

            auto tex_xxx_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 12);

            auto tey_xxx_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 12);

            auto tez_xxx_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 12);

            auto tex_xxx_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 13);

            auto tey_xxx_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 13);

            auto tez_xxx_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 13);

            auto tex_xxx_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 14);

            auto tey_xxx_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 14);

            auto tez_xxx_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 14);

            auto tex_xxy_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 15);

            auto tey_xxy_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 15);

            auto tez_xxy_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 15);

            auto tex_xxy_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 16);

            auto tex_xxx_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx);

            auto tey_xxx_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx);

            auto tez_xxx_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx);

            auto tex_xxx_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 1);

            auto tey_xxx_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 1);

            auto tez_xxx_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 1);

            auto tex_xxx_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 2);

            auto tey_xxx_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 2);

            auto tez_xxx_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 2);

            auto tex_xxx_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 3);

            auto tey_xxx_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 3);

            auto tez_xxx_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 3);

            auto tex_xxx_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 4);

            auto tey_xxx_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 4);

            auto tez_xxx_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 4);

            auto tex_xxx_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 5);

            auto tey_xxx_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 5);

            auto tez_xxx_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 5);

            auto tex_xxx_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 6);

            auto tey_xxx_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 6);

            auto tez_xxx_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 6);

            auto tex_xxx_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 7);

            auto tey_xxx_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 7);

            auto tez_xxx_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 7);

            auto tex_xxx_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 8);

            auto tey_xxx_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 8);

            auto tez_xxx_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 8);

            auto tex_xxx_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 9);

            auto tey_xxx_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 9);

            auto tez_xxx_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 9);

            auto tex_xxx_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 10);

            auto tey_xxx_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 10);

            auto tez_xxx_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 10);

            auto tex_xxx_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 11);

            auto tey_xxx_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 11);

            auto tez_xxx_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 11);

            auto tex_xxx_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 12);

            auto tey_xxx_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 12);

            auto tez_xxx_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 12);

            auto tex_xxx_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 13);

            auto tey_xxx_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 13);

            auto tez_xxx_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 13);

            auto tex_xxx_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 14);

            auto tey_xxx_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 14);

            auto tez_xxx_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 14);

            auto tex_xxy_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 15);

            auto tey_xxy_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 15);

            auto tez_xxy_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 15);

            auto tex_xxy_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 16);

            auto tex_xx_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx);

            auto tey_xx_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx);

            auto tez_xx_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx);

            auto tex_xx_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 1);

            auto tey_xx_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 1);

            auto tez_xx_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 1);

            auto tex_xx_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 2);

            auto tey_xx_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 2);

            auto tez_xx_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 2);

            auto tex_xx_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 3);

            auto tey_xx_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 3);

            auto tez_xx_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 3);

            auto tex_xx_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 4);

            auto tey_xx_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 4);

            auto tez_xx_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 4);

            auto tex_xx_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 5);

            auto tey_xx_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 5);

            auto tez_xx_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 5);

            auto tex_xx_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 6);

            auto tey_xx_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 6);

            auto tez_xx_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 6);

            auto tex_xx_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 7);

            auto tey_xx_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 7);

            auto tez_xx_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 7);

            auto tex_xx_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 8);

            auto tey_xx_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 8);

            auto tez_xx_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 8);

            auto tex_xx_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 9);

            auto tey_xx_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 9);

            auto tez_xx_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 9);

            auto tex_xx_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 10);

            auto tey_xx_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 10);

            auto tez_xx_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 10);

            auto tex_xx_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 11);

            auto tey_xx_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 11);

            auto tez_xx_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 11);

            auto tex_xx_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 12);

            auto tey_xx_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 12);

            auto tez_xx_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 12);

            auto tex_xx_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 13);

            auto tey_xx_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 13);

            auto tez_xx_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 13);

            auto tex_xx_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 14);

            auto tey_xx_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 14);

            auto tez_xx_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 14);

            auto tex_xy_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 15);

            auto tey_xy_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 15);

            auto tez_xy_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 15);

            auto tex_xy_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 16);

            auto tex_xx_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx);

            auto tey_xx_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx);

            auto tez_xx_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx);

            auto tex_xx_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 1);

            auto tey_xx_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 1);

            auto tez_xx_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 1);

            auto tex_xx_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 2);

            auto tey_xx_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 2);

            auto tez_xx_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 2);

            auto tex_xx_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 3);

            auto tey_xx_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 3);

            auto tez_xx_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 3);

            auto tex_xx_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 4);

            auto tey_xx_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 4);

            auto tez_xx_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 4);

            auto tex_xx_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 5);

            auto tey_xx_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 5);

            auto tez_xx_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 5);

            auto tex_xx_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 6);

            auto tey_xx_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 6);

            auto tez_xx_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 6);

            auto tex_xx_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 7);

            auto tey_xx_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 7);

            auto tez_xx_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 7);

            auto tex_xx_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 8);

            auto tey_xx_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 8);

            auto tez_xx_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 8);

            auto tex_xx_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 9);

            auto tey_xx_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 9);

            auto tez_xx_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 9);

            auto tex_xx_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 10);

            auto tey_xx_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 10);

            auto tez_xx_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 10);

            auto tex_xx_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 11);

            auto tey_xx_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 11);

            auto tez_xx_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 11);

            auto tex_xx_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 12);

            auto tey_xx_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 12);

            auto tez_xx_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 12);

            auto tex_xx_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 13);

            auto tey_xx_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 13);

            auto tez_xx_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 13);

            auto tex_xx_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 14);

            auto tey_xx_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 14);

            auto tez_xx_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 14);

            auto tex_xy_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 15);

            auto tey_xy_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 15);

            auto tez_xy_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 15);

            auto tex_xy_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 16);

            auto tex_xxx_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx);

            auto tey_xxx_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx);

            auto tez_xxx_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx);

            auto tex_xxx_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 1);

            auto tey_xxx_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 1);

            auto tez_xxx_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 1);

            auto tex_xxx_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 2);

            auto tey_xxx_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 2);

            auto tez_xxx_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 2);

            auto tex_xxx_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 3);

            auto tey_xxx_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 3);

            auto tez_xxx_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 3);

            auto tex_xxx_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 4);

            auto tey_xxx_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 4);

            auto tez_xxx_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 4);

            auto tex_xxx_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 5);

            auto tey_xxx_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 5);

            auto tez_xxx_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 5);

            auto tex_xxx_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 6);

            auto tey_xxx_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 6);

            auto tez_xxx_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 6);

            auto tex_xxx_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 7);

            auto tey_xxx_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 7);

            auto tez_xxx_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 7);

            auto tex_xxx_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 8);

            auto tey_xxx_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 8);

            auto tez_xxx_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 8);

            auto tex_xxx_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 9);

            auto tey_xxx_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 9);

            auto tez_xxx_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 9);

            auto tex_xxy_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 10);

            auto tey_xxy_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 10);

            auto tez_xxy_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 10);

            auto tex_xxy_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 11);

            auto tex_xxx_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx);

            auto tey_xxx_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx);

            auto tez_xxx_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx);

            auto tex_xxx_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 1);

            auto tey_xxx_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 1);

            auto tez_xxx_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 1);

            auto tex_xxx_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 2);

            auto tey_xxx_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 2);

            auto tez_xxx_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 2);

            auto tex_xxx_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 3);

            auto tey_xxx_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 3);

            auto tez_xxx_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 3);

            auto tex_xxx_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 4);

            auto tey_xxx_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 4);

            auto tez_xxx_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 4);

            auto tex_xxx_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 5);

            auto tey_xxx_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 5);

            auto tez_xxx_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 5);

            auto tex_xxx_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 6);

            auto tey_xxx_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 6);

            auto tez_xxx_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 6);

            auto tex_xxx_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 7);

            auto tey_xxx_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 7);

            auto tez_xxx_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 7);

            auto tex_xxx_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 8);

            auto tey_xxx_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 8);

            auto tez_xxx_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 8);

            auto tex_xxx_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 9);

            auto tey_xxx_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 9);

            auto tez_xxx_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 9);

            auto tex_xxy_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 10);

            auto tey_xxy_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 10);

            auto tez_xxy_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 10);

            auto tex_xxy_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 11);

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

            // set up pointers to integrals

            auto tex_xxxx_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx);

            auto tey_xxxx_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx);

            auto tez_xxxx_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx);

            auto tex_xxxx_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 1);

            auto tey_xxxx_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 1);

            auto tez_xxxx_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 1);

            auto tex_xxxx_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 2);

            auto tey_xxxx_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 2);

            auto tez_xxxx_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 2);

            auto tex_xxxx_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 3);

            auto tey_xxxx_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 3);

            auto tez_xxxx_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 3);

            auto tex_xxxx_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 4);

            auto tey_xxxx_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 4);

            auto tez_xxxx_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 4);

            auto tex_xxxx_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 5);

            auto tey_xxxx_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 5);

            auto tez_xxxx_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 5);

            auto tex_xxxx_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 6);

            auto tey_xxxx_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 6);

            auto tez_xxxx_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 6);

            auto tex_xxxx_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 7);

            auto tey_xxxx_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 7);

            auto tez_xxxx_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 7);

            auto tex_xxxx_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 8);

            auto tey_xxxx_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 8);

            auto tez_xxxx_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 8);

            auto tex_xxxx_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 9);

            auto tey_xxxx_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 9);

            auto tez_xxxx_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 9);

            auto tex_xxxx_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 10);

            auto tey_xxxx_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 10);

            auto tez_xxxx_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 10);

            auto tex_xxxx_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 11);

            auto tey_xxxx_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 11);

            auto tez_xxxx_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 11);

            auto tex_xxxx_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 12);

            auto tey_xxxx_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 12);

            auto tez_xxxx_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 12);

            auto tex_xxxx_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 13);

            auto tey_xxxx_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 13);

            auto tez_xxxx_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 13);

            auto tex_xxxx_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 14);

            auto tey_xxxx_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 14);

            auto tez_xxxx_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 14);

            auto tex_xxxy_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 15);

            auto tey_xxxy_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 15);

            auto tez_xxxy_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 15);

            auto tex_xxxy_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 16);

            // Batch of Integrals (0,49)

            #pragma omp simd aligned(fx, pa_x, pc_x, ta_xxx_xxxx_1, ta_xxx_xxxy_1, ta_xxx_xxxz_1, \
                                         ta_xxx_xxyy_1, ta_xxx_xxyz_1, ta_xxx_xxzz_1, ta_xxx_xyyy_1, ta_xxx_xyyz_1, \
                                         ta_xxx_xyzz_1, ta_xxx_xzzz_1, ta_xxx_yyyy_1, ta_xxx_yyyz_1, ta_xxx_yyzz_1, \
                                         ta_xxx_yzzz_1, ta_xxx_zzzz_1, ta_xxy_xxxx_1, ta_xxy_xxxy_1, tex_xx_xxxx_0, \
                                         tex_xx_xxxx_1, tex_xx_xxxy_0, tex_xx_xxxy_1, tex_xx_xxxz_0, tex_xx_xxxz_1, \
                                         tex_xx_xxyy_0, tex_xx_xxyy_1, tex_xx_xxyz_0, tex_xx_xxyz_1, tex_xx_xxzz_0, \
                                         tex_xx_xxzz_1, tex_xx_xyyy_0, tex_xx_xyyy_1, tex_xx_xyyz_0, tex_xx_xyyz_1, \
                                         tex_xx_xyzz_0, tex_xx_xyzz_1, tex_xx_xzzz_0, tex_xx_xzzz_1, tex_xx_yyyy_0, \
                                         tex_xx_yyyy_1, tex_xx_yyyz_0, tex_xx_yyyz_1, tex_xx_yyzz_0, tex_xx_yyzz_1, \
                                         tex_xx_yzzz_0, tex_xx_yzzz_1, tex_xx_zzzz_0, tex_xx_zzzz_1, tex_xxx_xxx_0, \
                                         tex_xxx_xxx_1, tex_xxx_xxxx_0, tex_xxx_xxxx_1, tex_xxx_xxxy_0, tex_xxx_xxxy_1, \
                                         tex_xxx_xxxz_0, tex_xxx_xxxz_1, tex_xxx_xxy_0, tex_xxx_xxy_1, tex_xxx_xxyy_0, \
                                         tex_xxx_xxyy_1, tex_xxx_xxyz_0, tex_xxx_xxyz_1, tex_xxx_xxz_0, tex_xxx_xxz_1, \
                                         tex_xxx_xxzz_0, tex_xxx_xxzz_1, tex_xxx_xyy_0, tex_xxx_xyy_1, tex_xxx_xyyy_0, \
                                         tex_xxx_xyyy_1, tex_xxx_xyyz_0, tex_xxx_xyyz_1, tex_xxx_xyz_0, tex_xxx_xyz_1, \
                                         tex_xxx_xyzz_0, tex_xxx_xyzz_1, tex_xxx_xzz_0, tex_xxx_xzz_1, tex_xxx_xzzz_0, \
                                         tex_xxx_xzzz_1, tex_xxx_yyy_0, tex_xxx_yyy_1, tex_xxx_yyyy_0, tex_xxx_yyyy_1, \
                                         tex_xxx_yyyz_0, tex_xxx_yyyz_1, tex_xxx_yyz_0, tex_xxx_yyz_1, tex_xxx_yyzz_0, \
                                         tex_xxx_yyzz_1, tex_xxx_yzz_0, tex_xxx_yzz_1, tex_xxx_yzzz_0, tex_xxx_yzzz_1, \
                                         tex_xxx_zzz_0, tex_xxx_zzz_1, tex_xxx_zzzz_0, tex_xxx_zzzz_1, tex_xxxx_xxxx_0, \
                                         tex_xxxx_xxxy_0, tex_xxxx_xxxz_0, tex_xxxx_xxyy_0, tex_xxxx_xxyz_0, tex_xxxx_xxzz_0, \
                                         tex_xxxx_xyyy_0, tex_xxxx_xyyz_0, tex_xxxx_xyzz_0, tex_xxxx_xzzz_0, tex_xxxx_yyyy_0, \
                                         tex_xxxx_yyyz_0, tex_xxxx_yyzz_0, tex_xxxx_yzzz_0, tex_xxxx_zzzz_0, tex_xxxy_xxxx_0, \
                                         tex_xxxy_xxxy_0, tex_xxy_xxx_0, tex_xxy_xxx_1, tex_xxy_xxxx_0, tex_xxy_xxxx_1, \
                                         tex_xxy_xxxy_0, tex_xxy_xxxy_1, tex_xxy_xxy_0, tex_xxy_xxy_1, tex_xy_xxxx_0, \
                                         tex_xy_xxxx_1, tex_xy_xxxy_0, tex_xy_xxxy_1, tey_xx_xxxx_0, tey_xx_xxxx_1, \
                                         tey_xx_xxxy_0, tey_xx_xxxy_1, tey_xx_xxxz_0, tey_xx_xxxz_1, tey_xx_xxyy_0, \
                                         tey_xx_xxyy_1, tey_xx_xxyz_0, tey_xx_xxyz_1, tey_xx_xxzz_0, tey_xx_xxzz_1, \
                                         tey_xx_xyyy_0, tey_xx_xyyy_1, tey_xx_xyyz_0, tey_xx_xyyz_1, tey_xx_xyzz_0, \
                                         tey_xx_xyzz_1, tey_xx_xzzz_0, tey_xx_xzzz_1, tey_xx_yyyy_0, tey_xx_yyyy_1, \
                                         tey_xx_yyyz_0, tey_xx_yyyz_1, tey_xx_yyzz_0, tey_xx_yyzz_1, tey_xx_yzzz_0, \
                                         tey_xx_yzzz_1, tey_xx_zzzz_0, tey_xx_zzzz_1, tey_xxx_xxx_0, tey_xxx_xxx_1, \
                                         tey_xxx_xxxx_0, tey_xxx_xxxx_1, tey_xxx_xxxy_0, tey_xxx_xxxy_1, tey_xxx_xxxz_0, \
                                         tey_xxx_xxxz_1, tey_xxx_xxy_0, tey_xxx_xxy_1, tey_xxx_xxyy_0, tey_xxx_xxyy_1, \
                                         tey_xxx_xxyz_0, tey_xxx_xxyz_1, tey_xxx_xxz_0, tey_xxx_xxz_1, tey_xxx_xxzz_0, \
                                         tey_xxx_xxzz_1, tey_xxx_xyy_0, tey_xxx_xyy_1, tey_xxx_xyyy_0, tey_xxx_xyyy_1, \
                                         tey_xxx_xyyz_0, tey_xxx_xyyz_1, tey_xxx_xyz_0, tey_xxx_xyz_1, tey_xxx_xyzz_0, \
                                         tey_xxx_xyzz_1, tey_xxx_xzz_0, tey_xxx_xzz_1, tey_xxx_xzzz_0, tey_xxx_xzzz_1, \
                                         tey_xxx_yyy_0, tey_xxx_yyy_1, tey_xxx_yyyy_0, tey_xxx_yyyy_1, tey_xxx_yyyz_0, \
                                         tey_xxx_yyyz_1, tey_xxx_yyz_0, tey_xxx_yyz_1, tey_xxx_yyzz_0, tey_xxx_yyzz_1, \
                                         tey_xxx_yzz_0, tey_xxx_yzz_1, tey_xxx_yzzz_0, tey_xxx_yzzz_1, tey_xxx_zzz_0, \
                                         tey_xxx_zzz_1, tey_xxx_zzzz_0, tey_xxx_zzzz_1, tey_xxxx_xxxx_0, tey_xxxx_xxxy_0, \
                                         tey_xxxx_xxxz_0, tey_xxxx_xxyy_0, tey_xxxx_xxyz_0, tey_xxxx_xxzz_0, tey_xxxx_xyyy_0, \
                                         tey_xxxx_xyyz_0, tey_xxxx_xyzz_0, tey_xxxx_xzzz_0, tey_xxxx_yyyy_0, tey_xxxx_yyyz_0, \
                                         tey_xxxx_yyzz_0, tey_xxxx_yzzz_0, tey_xxxx_zzzz_0, tey_xxxy_xxxx_0, tey_xxy_xxx_0, \
                                         tey_xxy_xxx_1, tey_xxy_xxxx_0, tey_xxy_xxxx_1, tey_xy_xxxx_0, tey_xy_xxxx_1, \
                                         tez_xx_xxxx_0, tez_xx_xxxx_1, tez_xx_xxxy_0, tez_xx_xxxy_1, tez_xx_xxxz_0, \
                                         tez_xx_xxxz_1, tez_xx_xxyy_0, tez_xx_xxyy_1, tez_xx_xxyz_0, tez_xx_xxyz_1, \
                                         tez_xx_xxzz_0, tez_xx_xxzz_1, tez_xx_xyyy_0, tez_xx_xyyy_1, tez_xx_xyyz_0, \
                                         tez_xx_xyyz_1, tez_xx_xyzz_0, tez_xx_xyzz_1, tez_xx_xzzz_0, tez_xx_xzzz_1, \
                                         tez_xx_yyyy_0, tez_xx_yyyy_1, tez_xx_yyyz_0, tez_xx_yyyz_1, tez_xx_yyzz_0, \
                                         tez_xx_yyzz_1, tez_xx_yzzz_0, tez_xx_yzzz_1, tez_xx_zzzz_0, tez_xx_zzzz_1, \
                                         tez_xxx_xxx_0, tez_xxx_xxx_1, tez_xxx_xxxx_0, tez_xxx_xxxx_1, tez_xxx_xxxy_0, \
                                         tez_xxx_xxxy_1, tez_xxx_xxxz_0, tez_xxx_xxxz_1, tez_xxx_xxy_0, tez_xxx_xxy_1, \
                                         tez_xxx_xxyy_0, tez_xxx_xxyy_1, tez_xxx_xxyz_0, tez_xxx_xxyz_1, tez_xxx_xxz_0, \
                                         tez_xxx_xxz_1, tez_xxx_xxzz_0, tez_xxx_xxzz_1, tez_xxx_xyy_0, tez_xxx_xyy_1, \
                                         tez_xxx_xyyy_0, tez_xxx_xyyy_1, tez_xxx_xyyz_0, tez_xxx_xyyz_1, tez_xxx_xyz_0, \
                                         tez_xxx_xyz_1, tez_xxx_xyzz_0, tez_xxx_xyzz_1, tez_xxx_xzz_0, tez_xxx_xzz_1, \
                                         tez_xxx_xzzz_0, tez_xxx_xzzz_1, tez_xxx_yyy_0, tez_xxx_yyy_1, tez_xxx_yyyy_0, \
                                         tez_xxx_yyyy_1, tez_xxx_yyyz_0, tez_xxx_yyyz_1, tez_xxx_yyz_0, tez_xxx_yyz_1, \
                                         tez_xxx_yyzz_0, tez_xxx_yyzz_1, tez_xxx_yzz_0, tez_xxx_yzz_1, tez_xxx_yzzz_0, \
                                         tez_xxx_yzzz_1, tez_xxx_zzz_0, tez_xxx_zzz_1, tez_xxx_zzzz_0, tez_xxx_zzzz_1, \
                                         tez_xxxx_xxxx_0, tez_xxxx_xxxy_0, tez_xxxx_xxxz_0, tez_xxxx_xxyy_0, tez_xxxx_xxyz_0, \
                                         tez_xxxx_xxzz_0, tez_xxxx_xyyy_0, tez_xxxx_xyyz_0, tez_xxxx_xyzz_0, tez_xxxx_xzzz_0, \
                                         tez_xxxx_yyyy_0, tez_xxxx_yyyz_0, tez_xxxx_yyzz_0, tez_xxxx_yzzz_0, tez_xxxx_zzzz_0, \
                                         tez_xxxy_xxxx_0, tez_xxy_xxx_0, tez_xxy_xxx_1, tez_xxy_xxxx_0, tez_xxy_xxxx_1, \
                                         tez_xy_xxxx_0, tez_xy_xxxx_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_xxxx_xxxx_0[j] = pa_x[j] * tex_xxx_xxxx_0[j] - pc_x[j] * tex_xxx_xxxx_1[j] + 1.5 * fl1_fx * tex_xx_xxxx_0[j] -
                                     1.5 * fl1_fx * tex_xx_xxxx_1[j] + 2.0 * fl1_fx * tex_xxx_xxx_0[j] - 2.0 * fl1_fx * tex_xxx_xxx_1[j] +
                                     ta_xxx_xxxx_1[j];

                tey_xxxx_xxxx_0[j] = pa_x[j] * tey_xxx_xxxx_0[j] - pc_x[j] * tey_xxx_xxxx_1[j] + 1.5 * fl1_fx * tey_xx_xxxx_0[j] -
                                     1.5 * fl1_fx * tey_xx_xxxx_1[j] + 2.0 * fl1_fx * tey_xxx_xxx_0[j] - 2.0 * fl1_fx * tey_xxx_xxx_1[j];

                tez_xxxx_xxxx_0[j] = pa_x[j] * tez_xxx_xxxx_0[j] - pc_x[j] * tez_xxx_xxxx_1[j] + 1.5 * fl1_fx * tez_xx_xxxx_0[j] -
                                     1.5 * fl1_fx * tez_xx_xxxx_1[j] + 2.0 * fl1_fx * tez_xxx_xxx_0[j] - 2.0 * fl1_fx * tez_xxx_xxx_1[j];

                tex_xxxx_xxxy_0[j] = pa_x[j] * tex_xxx_xxxy_0[j] - pc_x[j] * tex_xxx_xxxy_1[j] + 1.5 * fl1_fx * tex_xx_xxxy_0[j] -
                                     1.5 * fl1_fx * tex_xx_xxxy_1[j] + 1.5 * fl1_fx * tex_xxx_xxy_0[j] - 1.5 * fl1_fx * tex_xxx_xxy_1[j] +
                                     ta_xxx_xxxy_1[j];

                tey_xxxx_xxxy_0[j] = pa_x[j] * tey_xxx_xxxy_0[j] - pc_x[j] * tey_xxx_xxxy_1[j] + 1.5 * fl1_fx * tey_xx_xxxy_0[j] -
                                     1.5 * fl1_fx * tey_xx_xxxy_1[j] + 1.5 * fl1_fx * tey_xxx_xxy_0[j] - 1.5 * fl1_fx * tey_xxx_xxy_1[j];

                tez_xxxx_xxxy_0[j] = pa_x[j] * tez_xxx_xxxy_0[j] - pc_x[j] * tez_xxx_xxxy_1[j] + 1.5 * fl1_fx * tez_xx_xxxy_0[j] -
                                     1.5 * fl1_fx * tez_xx_xxxy_1[j] + 1.5 * fl1_fx * tez_xxx_xxy_0[j] - 1.5 * fl1_fx * tez_xxx_xxy_1[j];

                tex_xxxx_xxxz_0[j] = pa_x[j] * tex_xxx_xxxz_0[j] - pc_x[j] * tex_xxx_xxxz_1[j] + 1.5 * fl1_fx * tex_xx_xxxz_0[j] -
                                     1.5 * fl1_fx * tex_xx_xxxz_1[j] + 1.5 * fl1_fx * tex_xxx_xxz_0[j] - 1.5 * fl1_fx * tex_xxx_xxz_1[j] +
                                     ta_xxx_xxxz_1[j];

                tey_xxxx_xxxz_0[j] = pa_x[j] * tey_xxx_xxxz_0[j] - pc_x[j] * tey_xxx_xxxz_1[j] + 1.5 * fl1_fx * tey_xx_xxxz_0[j] -
                                     1.5 * fl1_fx * tey_xx_xxxz_1[j] + 1.5 * fl1_fx * tey_xxx_xxz_0[j] - 1.5 * fl1_fx * tey_xxx_xxz_1[j];

                tez_xxxx_xxxz_0[j] = pa_x[j] * tez_xxx_xxxz_0[j] - pc_x[j] * tez_xxx_xxxz_1[j] + 1.5 * fl1_fx * tez_xx_xxxz_0[j] -
                                     1.5 * fl1_fx * tez_xx_xxxz_1[j] + 1.5 * fl1_fx * tez_xxx_xxz_0[j] - 1.5 * fl1_fx * tez_xxx_xxz_1[j];

                tex_xxxx_xxyy_0[j] = pa_x[j] * tex_xxx_xxyy_0[j] - pc_x[j] * tex_xxx_xxyy_1[j] + 1.5 * fl1_fx * tex_xx_xxyy_0[j] -
                                     1.5 * fl1_fx * tex_xx_xxyy_1[j] + fl1_fx * tex_xxx_xyy_0[j] - fl1_fx * tex_xxx_xyy_1[j] + ta_xxx_xxyy_1[j];

                tey_xxxx_xxyy_0[j] = pa_x[j] * tey_xxx_xxyy_0[j] - pc_x[j] * tey_xxx_xxyy_1[j] + 1.5 * fl1_fx * tey_xx_xxyy_0[j] -
                                     1.5 * fl1_fx * tey_xx_xxyy_1[j] + fl1_fx * tey_xxx_xyy_0[j] - fl1_fx * tey_xxx_xyy_1[j];

                tez_xxxx_xxyy_0[j] = pa_x[j] * tez_xxx_xxyy_0[j] - pc_x[j] * tez_xxx_xxyy_1[j] + 1.5 * fl1_fx * tez_xx_xxyy_0[j] -
                                     1.5 * fl1_fx * tez_xx_xxyy_1[j] + fl1_fx * tez_xxx_xyy_0[j] - fl1_fx * tez_xxx_xyy_1[j];

                tex_xxxx_xxyz_0[j] = pa_x[j] * tex_xxx_xxyz_0[j] - pc_x[j] * tex_xxx_xxyz_1[j] + 1.5 * fl1_fx * tex_xx_xxyz_0[j] -
                                     1.5 * fl1_fx * tex_xx_xxyz_1[j] + fl1_fx * tex_xxx_xyz_0[j] - fl1_fx * tex_xxx_xyz_1[j] + ta_xxx_xxyz_1[j];

                tey_xxxx_xxyz_0[j] = pa_x[j] * tey_xxx_xxyz_0[j] - pc_x[j] * tey_xxx_xxyz_1[j] + 1.5 * fl1_fx * tey_xx_xxyz_0[j] -
                                     1.5 * fl1_fx * tey_xx_xxyz_1[j] + fl1_fx * tey_xxx_xyz_0[j] - fl1_fx * tey_xxx_xyz_1[j];

                tez_xxxx_xxyz_0[j] = pa_x[j] * tez_xxx_xxyz_0[j] - pc_x[j] * tez_xxx_xxyz_1[j] + 1.5 * fl1_fx * tez_xx_xxyz_0[j] -
                                     1.5 * fl1_fx * tez_xx_xxyz_1[j] + fl1_fx * tez_xxx_xyz_0[j] - fl1_fx * tez_xxx_xyz_1[j];

                tex_xxxx_xxzz_0[j] = pa_x[j] * tex_xxx_xxzz_0[j] - pc_x[j] * tex_xxx_xxzz_1[j] + 1.5 * fl1_fx * tex_xx_xxzz_0[j] -
                                     1.5 * fl1_fx * tex_xx_xxzz_1[j] + fl1_fx * tex_xxx_xzz_0[j] - fl1_fx * tex_xxx_xzz_1[j] + ta_xxx_xxzz_1[j];

                tey_xxxx_xxzz_0[j] = pa_x[j] * tey_xxx_xxzz_0[j] - pc_x[j] * tey_xxx_xxzz_1[j] + 1.5 * fl1_fx * tey_xx_xxzz_0[j] -
                                     1.5 * fl1_fx * tey_xx_xxzz_1[j] + fl1_fx * tey_xxx_xzz_0[j] - fl1_fx * tey_xxx_xzz_1[j];

                tez_xxxx_xxzz_0[j] = pa_x[j] * tez_xxx_xxzz_0[j] - pc_x[j] * tez_xxx_xxzz_1[j] + 1.5 * fl1_fx * tez_xx_xxzz_0[j] -
                                     1.5 * fl1_fx * tez_xx_xxzz_1[j] + fl1_fx * tez_xxx_xzz_0[j] - fl1_fx * tez_xxx_xzz_1[j];

                tex_xxxx_xyyy_0[j] = pa_x[j] * tex_xxx_xyyy_0[j] - pc_x[j] * tex_xxx_xyyy_1[j] + 1.5 * fl1_fx * tex_xx_xyyy_0[j] -
                                     1.5 * fl1_fx * tex_xx_xyyy_1[j] + 0.5 * fl1_fx * tex_xxx_yyy_0[j] - 0.5 * fl1_fx * tex_xxx_yyy_1[j] +
                                     ta_xxx_xyyy_1[j];

                tey_xxxx_xyyy_0[j] = pa_x[j] * tey_xxx_xyyy_0[j] - pc_x[j] * tey_xxx_xyyy_1[j] + 1.5 * fl1_fx * tey_xx_xyyy_0[j] -
                                     1.5 * fl1_fx * tey_xx_xyyy_1[j] + 0.5 * fl1_fx * tey_xxx_yyy_0[j] - 0.5 * fl1_fx * tey_xxx_yyy_1[j];

                tez_xxxx_xyyy_0[j] = pa_x[j] * tez_xxx_xyyy_0[j] - pc_x[j] * tez_xxx_xyyy_1[j] + 1.5 * fl1_fx * tez_xx_xyyy_0[j] -
                                     1.5 * fl1_fx * tez_xx_xyyy_1[j] + 0.5 * fl1_fx * tez_xxx_yyy_0[j] - 0.5 * fl1_fx * tez_xxx_yyy_1[j];

                tex_xxxx_xyyz_0[j] = pa_x[j] * tex_xxx_xyyz_0[j] - pc_x[j] * tex_xxx_xyyz_1[j] + 1.5 * fl1_fx * tex_xx_xyyz_0[j] -
                                     1.5 * fl1_fx * tex_xx_xyyz_1[j] + 0.5 * fl1_fx * tex_xxx_yyz_0[j] - 0.5 * fl1_fx * tex_xxx_yyz_1[j] +
                                     ta_xxx_xyyz_1[j];

                tey_xxxx_xyyz_0[j] = pa_x[j] * tey_xxx_xyyz_0[j] - pc_x[j] * tey_xxx_xyyz_1[j] + 1.5 * fl1_fx * tey_xx_xyyz_0[j] -
                                     1.5 * fl1_fx * tey_xx_xyyz_1[j] + 0.5 * fl1_fx * tey_xxx_yyz_0[j] - 0.5 * fl1_fx * tey_xxx_yyz_1[j];

                tez_xxxx_xyyz_0[j] = pa_x[j] * tez_xxx_xyyz_0[j] - pc_x[j] * tez_xxx_xyyz_1[j] + 1.5 * fl1_fx * tez_xx_xyyz_0[j] -
                                     1.5 * fl1_fx * tez_xx_xyyz_1[j] + 0.5 * fl1_fx * tez_xxx_yyz_0[j] - 0.5 * fl1_fx * tez_xxx_yyz_1[j];

                tex_xxxx_xyzz_0[j] = pa_x[j] * tex_xxx_xyzz_0[j] - pc_x[j] * tex_xxx_xyzz_1[j] + 1.5 * fl1_fx * tex_xx_xyzz_0[j] -
                                     1.5 * fl1_fx * tex_xx_xyzz_1[j] + 0.5 * fl1_fx * tex_xxx_yzz_0[j] - 0.5 * fl1_fx * tex_xxx_yzz_1[j] +
                                     ta_xxx_xyzz_1[j];

                tey_xxxx_xyzz_0[j] = pa_x[j] * tey_xxx_xyzz_0[j] - pc_x[j] * tey_xxx_xyzz_1[j] + 1.5 * fl1_fx * tey_xx_xyzz_0[j] -
                                     1.5 * fl1_fx * tey_xx_xyzz_1[j] + 0.5 * fl1_fx * tey_xxx_yzz_0[j] - 0.5 * fl1_fx * tey_xxx_yzz_1[j];

                tez_xxxx_xyzz_0[j] = pa_x[j] * tez_xxx_xyzz_0[j] - pc_x[j] * tez_xxx_xyzz_1[j] + 1.5 * fl1_fx * tez_xx_xyzz_0[j] -
                                     1.5 * fl1_fx * tez_xx_xyzz_1[j] + 0.5 * fl1_fx * tez_xxx_yzz_0[j] - 0.5 * fl1_fx * tez_xxx_yzz_1[j];

                tex_xxxx_xzzz_0[j] = pa_x[j] * tex_xxx_xzzz_0[j] - pc_x[j] * tex_xxx_xzzz_1[j] + 1.5 * fl1_fx * tex_xx_xzzz_0[j] -
                                     1.5 * fl1_fx * tex_xx_xzzz_1[j] + 0.5 * fl1_fx * tex_xxx_zzz_0[j] - 0.5 * fl1_fx * tex_xxx_zzz_1[j] +
                                     ta_xxx_xzzz_1[j];

                tey_xxxx_xzzz_0[j] = pa_x[j] * tey_xxx_xzzz_0[j] - pc_x[j] * tey_xxx_xzzz_1[j] + 1.5 * fl1_fx * tey_xx_xzzz_0[j] -
                                     1.5 * fl1_fx * tey_xx_xzzz_1[j] + 0.5 * fl1_fx * tey_xxx_zzz_0[j] - 0.5 * fl1_fx * tey_xxx_zzz_1[j];

                tez_xxxx_xzzz_0[j] = pa_x[j] * tez_xxx_xzzz_0[j] - pc_x[j] * tez_xxx_xzzz_1[j] + 1.5 * fl1_fx * tez_xx_xzzz_0[j] -
                                     1.5 * fl1_fx * tez_xx_xzzz_1[j] + 0.5 * fl1_fx * tez_xxx_zzz_0[j] - 0.5 * fl1_fx * tez_xxx_zzz_1[j];

                tex_xxxx_yyyy_0[j] = pa_x[j] * tex_xxx_yyyy_0[j] - pc_x[j] * tex_xxx_yyyy_1[j] + 1.5 * fl1_fx * tex_xx_yyyy_0[j] -
                                     1.5 * fl1_fx * tex_xx_yyyy_1[j] + ta_xxx_yyyy_1[j];

                tey_xxxx_yyyy_0[j] =
                    pa_x[j] * tey_xxx_yyyy_0[j] - pc_x[j] * tey_xxx_yyyy_1[j] + 1.5 * fl1_fx * tey_xx_yyyy_0[j] - 1.5 * fl1_fx * tey_xx_yyyy_1[j];

                tez_xxxx_yyyy_0[j] =
                    pa_x[j] * tez_xxx_yyyy_0[j] - pc_x[j] * tez_xxx_yyyy_1[j] + 1.5 * fl1_fx * tez_xx_yyyy_0[j] - 1.5 * fl1_fx * tez_xx_yyyy_1[j];

                tex_xxxx_yyyz_0[j] = pa_x[j] * tex_xxx_yyyz_0[j] - pc_x[j] * tex_xxx_yyyz_1[j] + 1.5 * fl1_fx * tex_xx_yyyz_0[j] -
                                     1.5 * fl1_fx * tex_xx_yyyz_1[j] + ta_xxx_yyyz_1[j];

                tey_xxxx_yyyz_0[j] =
                    pa_x[j] * tey_xxx_yyyz_0[j] - pc_x[j] * tey_xxx_yyyz_1[j] + 1.5 * fl1_fx * tey_xx_yyyz_0[j] - 1.5 * fl1_fx * tey_xx_yyyz_1[j];

                tez_xxxx_yyyz_0[j] =
                    pa_x[j] * tez_xxx_yyyz_0[j] - pc_x[j] * tez_xxx_yyyz_1[j] + 1.5 * fl1_fx * tez_xx_yyyz_0[j] - 1.5 * fl1_fx * tez_xx_yyyz_1[j];

                tex_xxxx_yyzz_0[j] = pa_x[j] * tex_xxx_yyzz_0[j] - pc_x[j] * tex_xxx_yyzz_1[j] + 1.5 * fl1_fx * tex_xx_yyzz_0[j] -
                                     1.5 * fl1_fx * tex_xx_yyzz_1[j] + ta_xxx_yyzz_1[j];

                tey_xxxx_yyzz_0[j] =
                    pa_x[j] * tey_xxx_yyzz_0[j] - pc_x[j] * tey_xxx_yyzz_1[j] + 1.5 * fl1_fx * tey_xx_yyzz_0[j] - 1.5 * fl1_fx * tey_xx_yyzz_1[j];

                tez_xxxx_yyzz_0[j] =
                    pa_x[j] * tez_xxx_yyzz_0[j] - pc_x[j] * tez_xxx_yyzz_1[j] + 1.5 * fl1_fx * tez_xx_yyzz_0[j] - 1.5 * fl1_fx * tez_xx_yyzz_1[j];

                tex_xxxx_yzzz_0[j] = pa_x[j] * tex_xxx_yzzz_0[j] - pc_x[j] * tex_xxx_yzzz_1[j] + 1.5 * fl1_fx * tex_xx_yzzz_0[j] -
                                     1.5 * fl1_fx * tex_xx_yzzz_1[j] + ta_xxx_yzzz_1[j];

                tey_xxxx_yzzz_0[j] =
                    pa_x[j] * tey_xxx_yzzz_0[j] - pc_x[j] * tey_xxx_yzzz_1[j] + 1.5 * fl1_fx * tey_xx_yzzz_0[j] - 1.5 * fl1_fx * tey_xx_yzzz_1[j];

                tez_xxxx_yzzz_0[j] =
                    pa_x[j] * tez_xxx_yzzz_0[j] - pc_x[j] * tez_xxx_yzzz_1[j] + 1.5 * fl1_fx * tez_xx_yzzz_0[j] - 1.5 * fl1_fx * tez_xx_yzzz_1[j];

                tex_xxxx_zzzz_0[j] = pa_x[j] * tex_xxx_zzzz_0[j] - pc_x[j] * tex_xxx_zzzz_1[j] + 1.5 * fl1_fx * tex_xx_zzzz_0[j] -
                                     1.5 * fl1_fx * tex_xx_zzzz_1[j] + ta_xxx_zzzz_1[j];

                tey_xxxx_zzzz_0[j] =
                    pa_x[j] * tey_xxx_zzzz_0[j] - pc_x[j] * tey_xxx_zzzz_1[j] + 1.5 * fl1_fx * tey_xx_zzzz_0[j] - 1.5 * fl1_fx * tey_xx_zzzz_1[j];

                tez_xxxx_zzzz_0[j] =
                    pa_x[j] * tez_xxx_zzzz_0[j] - pc_x[j] * tez_xxx_zzzz_1[j] + 1.5 * fl1_fx * tez_xx_zzzz_0[j] - 1.5 * fl1_fx * tez_xx_zzzz_1[j];

                tex_xxxy_xxxx_0[j] = pa_x[j] * tex_xxy_xxxx_0[j] - pc_x[j] * tex_xxy_xxxx_1[j] + fl1_fx * tex_xy_xxxx_0[j] -
                                     fl1_fx * tex_xy_xxxx_1[j] + 2.0 * fl1_fx * tex_xxy_xxx_0[j] - 2.0 * fl1_fx * tex_xxy_xxx_1[j] + ta_xxy_xxxx_1[j];

                tey_xxxy_xxxx_0[j] = pa_x[j] * tey_xxy_xxxx_0[j] - pc_x[j] * tey_xxy_xxxx_1[j] + fl1_fx * tey_xy_xxxx_0[j] -
                                     fl1_fx * tey_xy_xxxx_1[j] + 2.0 * fl1_fx * tey_xxy_xxx_0[j] - 2.0 * fl1_fx * tey_xxy_xxx_1[j];

                tez_xxxy_xxxx_0[j] = pa_x[j] * tez_xxy_xxxx_0[j] - pc_x[j] * tez_xxy_xxxx_1[j] + fl1_fx * tez_xy_xxxx_0[j] -
                                     fl1_fx * tez_xy_xxxx_1[j] + 2.0 * fl1_fx * tez_xxy_xxx_0[j] - 2.0 * fl1_fx * tez_xxy_xxx_1[j];

                tex_xxxy_xxxy_0[j] = pa_x[j] * tex_xxy_xxxy_0[j] - pc_x[j] * tex_xxy_xxxy_1[j] + fl1_fx * tex_xy_xxxy_0[j] -
                                     fl1_fx * tex_xy_xxxy_1[j] + 1.5 * fl1_fx * tex_xxy_xxy_0[j] - 1.5 * fl1_fx * tex_xxy_xxy_1[j] + ta_xxy_xxxy_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForGG_49_98(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CMemBlock2D<double>& pcDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto)
{
    // Batch of Integrals (49,98)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_4_4_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_2_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tey_xxy_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 16);

            auto tez_xxy_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 16);

            auto tex_xxy_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 17);

            auto tey_xxy_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 17);

            auto tez_xxy_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 17);

            auto tex_xxy_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 18);

            auto tey_xxy_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 18);

            auto tez_xxy_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 18);

            auto tex_xxy_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 19);

            auto tey_xxy_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 19);

            auto tez_xxy_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 19);

            auto tex_xxy_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 20);

            auto tey_xxy_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 20);

            auto tez_xxy_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 20);

            auto tex_xxy_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 21);

            auto tey_xxy_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 21);

            auto tez_xxy_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 21);

            auto tex_xxy_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 22);

            auto tey_xxy_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 22);

            auto tez_xxy_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 22);

            auto tex_xxy_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 23);

            auto tey_xxy_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 23);

            auto tez_xxy_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 23);

            auto tex_xxy_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 24);

            auto tey_xxy_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 24);

            auto tez_xxy_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 24);

            auto tex_xxy_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 25);

            auto tey_xxy_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 25);

            auto tez_xxy_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 25);

            auto tex_xxy_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 26);

            auto tey_xxy_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 26);

            auto tez_xxy_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 26);

            auto tex_xxy_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 27);

            auto tey_xxy_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 27);

            auto tez_xxy_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 27);

            auto tex_xxy_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 28);

            auto tey_xxy_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 28);

            auto tez_xxy_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 28);

            auto tex_xxy_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 29);

            auto tey_xxy_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 29);

            auto tez_xxy_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 29);

            auto tex_xxz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 30);

            auto tey_xxz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 30);

            auto tez_xxz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 30);

            auto tex_xxz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 31);

            auto tey_xxz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 31);

            auto tez_xxz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 31);

            auto tex_xxz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 32);

            auto tey_xxz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 32);

            auto tey_xxy_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 16);

            auto tez_xxy_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 16);

            auto tex_xxy_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 17);

            auto tey_xxy_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 17);

            auto tez_xxy_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 17);

            auto tex_xxy_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 18);

            auto tey_xxy_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 18);

            auto tez_xxy_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 18);

            auto tex_xxy_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 19);

            auto tey_xxy_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 19);

            auto tez_xxy_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 19);

            auto tex_xxy_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 20);

            auto tey_xxy_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 20);

            auto tez_xxy_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 20);

            auto tex_xxy_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 21);

            auto tey_xxy_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 21);

            auto tez_xxy_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 21);

            auto tex_xxy_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 22);

            auto tey_xxy_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 22);

            auto tez_xxy_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 22);

            auto tex_xxy_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 23);

            auto tey_xxy_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 23);

            auto tez_xxy_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 23);

            auto tex_xxy_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 24);

            auto tey_xxy_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 24);

            auto tez_xxy_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 24);

            auto tex_xxy_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 25);

            auto tey_xxy_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 25);

            auto tez_xxy_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 25);

            auto tex_xxy_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 26);

            auto tey_xxy_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 26);

            auto tez_xxy_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 26);

            auto tex_xxy_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 27);

            auto tey_xxy_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 27);

            auto tez_xxy_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 27);

            auto tex_xxy_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 28);

            auto tey_xxy_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 28);

            auto tez_xxy_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 28);

            auto tex_xxy_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 29);

            auto tey_xxy_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 29);

            auto tez_xxy_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 29);

            auto tex_xxz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 30);

            auto tey_xxz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 30);

            auto tez_xxz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 30);

            auto tex_xxz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 31);

            auto tey_xxz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 31);

            auto tez_xxz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 31);

            auto tex_xxz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 32);

            auto tey_xxz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 32);

            auto tey_xy_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 16);

            auto tez_xy_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 16);

            auto tex_xy_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 17);

            auto tey_xy_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 17);

            auto tez_xy_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 17);

            auto tex_xy_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 18);

            auto tey_xy_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 18);

            auto tez_xy_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 18);

            auto tex_xy_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 19);

            auto tey_xy_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 19);

            auto tez_xy_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 19);

            auto tex_xy_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 20);

            auto tey_xy_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 20);

            auto tez_xy_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 20);

            auto tex_xy_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 21);

            auto tey_xy_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 21);

            auto tez_xy_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 21);

            auto tex_xy_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 22);

            auto tey_xy_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 22);

            auto tez_xy_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 22);

            auto tex_xy_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 23);

            auto tey_xy_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 23);

            auto tez_xy_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 23);

            auto tex_xy_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 24);

            auto tey_xy_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 24);

            auto tez_xy_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 24);

            auto tex_xy_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 25);

            auto tey_xy_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 25);

            auto tez_xy_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 25);

            auto tex_xy_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 26);

            auto tey_xy_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 26);

            auto tez_xy_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 26);

            auto tex_xy_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 27);

            auto tey_xy_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 27);

            auto tez_xy_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 27);

            auto tex_xy_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 28);

            auto tey_xy_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 28);

            auto tez_xy_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 28);

            auto tex_xy_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 29);

            auto tey_xy_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 29);

            auto tez_xy_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 29);

            auto tex_xz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 30);

            auto tey_xz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 30);

            auto tez_xz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 30);

            auto tex_xz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 31);

            auto tey_xz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 31);

            auto tez_xz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 31);

            auto tex_xz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 32);

            auto tey_xz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 32);

            auto tey_xy_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 16);

            auto tez_xy_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 16);

            auto tex_xy_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 17);

            auto tey_xy_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 17);

            auto tez_xy_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 17);

            auto tex_xy_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 18);

            auto tey_xy_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 18);

            auto tez_xy_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 18);

            auto tex_xy_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 19);

            auto tey_xy_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 19);

            auto tez_xy_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 19);

            auto tex_xy_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 20);

            auto tey_xy_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 20);

            auto tez_xy_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 20);

            auto tex_xy_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 21);

            auto tey_xy_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 21);

            auto tez_xy_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 21);

            auto tex_xy_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 22);

            auto tey_xy_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 22);

            auto tez_xy_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 22);

            auto tex_xy_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 23);

            auto tey_xy_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 23);

            auto tez_xy_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 23);

            auto tex_xy_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 24);

            auto tey_xy_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 24);

            auto tez_xy_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 24);

            auto tex_xy_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 25);

            auto tey_xy_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 25);

            auto tez_xy_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 25);

            auto tex_xy_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 26);

            auto tey_xy_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 26);

            auto tez_xy_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 26);

            auto tex_xy_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 27);

            auto tey_xy_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 27);

            auto tez_xy_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 27);

            auto tex_xy_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 28);

            auto tey_xy_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 28);

            auto tez_xy_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 28);

            auto tex_xy_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 29);

            auto tey_xy_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 29);

            auto tez_xy_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 29);

            auto tex_xz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 30);

            auto tey_xz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 30);

            auto tez_xz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 30);

            auto tex_xz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 31);

            auto tey_xz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 31);

            auto tez_xz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 31);

            auto tex_xz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 32);

            auto tey_xz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 32);

            auto tey_xxy_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 11);

            auto tez_xxy_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 11);

            auto tex_xxy_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 12);

            auto tey_xxy_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 12);

            auto tez_xxy_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 12);

            auto tex_xxy_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 13);

            auto tey_xxy_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 13);

            auto tez_xxy_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 13);

            auto tex_xxy_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 14);

            auto tey_xxy_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 14);

            auto tez_xxy_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 14);

            auto tex_xxy_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 15);

            auto tey_xxy_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 15);

            auto tez_xxy_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 15);

            auto tex_xxy_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 16);

            auto tey_xxy_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 16);

            auto tez_xxy_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 16);

            auto tex_xxy_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 17);

            auto tey_xxy_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 17);

            auto tez_xxy_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 17);

            auto tex_xxy_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 18);

            auto tey_xxy_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 18);

            auto tez_xxy_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 18);

            auto tex_xxy_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 19);

            auto tey_xxy_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 19);

            auto tez_xxy_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 19);

            auto tex_xxz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 20);

            auto tey_xxz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 20);

            auto tez_xxz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 20);

            auto tex_xxz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 21);

            auto tey_xxz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 21);

            auto tez_xxz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 21);

            auto tex_xxz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 22);

            auto tey_xxz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 22);

            auto tey_xxy_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 11);

            auto tez_xxy_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 11);

            auto tex_xxy_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 12);

            auto tey_xxy_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 12);

            auto tez_xxy_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 12);

            auto tex_xxy_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 13);

            auto tey_xxy_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 13);

            auto tez_xxy_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 13);

            auto tex_xxy_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 14);

            auto tey_xxy_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 14);

            auto tez_xxy_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 14);

            auto tex_xxy_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 15);

            auto tey_xxy_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 15);

            auto tez_xxy_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 15);

            auto tex_xxy_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 16);

            auto tey_xxy_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 16);

            auto tez_xxy_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 16);

            auto tex_xxy_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 17);

            auto tey_xxy_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 17);

            auto tez_xxy_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 17);

            auto tex_xxy_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 18);

            auto tey_xxy_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 18);

            auto tez_xxy_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 18);

            auto tex_xxy_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 19);

            auto tey_xxy_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 19);

            auto tez_xxy_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 19);

            auto tex_xxz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 20);

            auto tey_xxz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 20);

            auto tez_xxz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 20);

            auto tex_xxz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 21);

            auto tey_xxz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 21);

            auto tez_xxz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 21);

            auto tex_xxz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 22);

            auto tey_xxz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 22);

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

            // set up pointers to integrals

            auto tey_xxxy_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 16);

            auto tez_xxxy_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 16);

            auto tex_xxxy_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 17);

            auto tey_xxxy_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 17);

            auto tez_xxxy_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 17);

            auto tex_xxxy_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 18);

            auto tey_xxxy_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 18);

            auto tez_xxxy_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 18);

            auto tex_xxxy_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 19);

            auto tey_xxxy_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 19);

            auto tez_xxxy_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 19);

            auto tex_xxxy_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 20);

            auto tey_xxxy_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 20);

            auto tez_xxxy_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 20);

            auto tex_xxxy_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 21);

            auto tey_xxxy_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 21);

            auto tez_xxxy_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 21);

            auto tex_xxxy_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 22);

            auto tey_xxxy_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 22);

            auto tez_xxxy_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 22);

            auto tex_xxxy_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 23);

            auto tey_xxxy_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 23);

            auto tez_xxxy_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 23);

            auto tex_xxxy_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 24);

            auto tey_xxxy_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 24);

            auto tez_xxxy_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 24);

            auto tex_xxxy_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 25);

            auto tey_xxxy_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 25);

            auto tez_xxxy_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 25);

            auto tex_xxxy_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 26);

            auto tey_xxxy_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 26);

            auto tez_xxxy_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 26);

            auto tex_xxxy_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 27);

            auto tey_xxxy_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 27);

            auto tez_xxxy_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 27);

            auto tex_xxxy_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 28);

            auto tey_xxxy_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 28);

            auto tez_xxxy_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 28);

            auto tex_xxxy_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 29);

            auto tey_xxxy_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 29);

            auto tez_xxxy_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 29);

            auto tex_xxxz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 30);

            auto tey_xxxz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 30);

            auto tez_xxxz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 30);

            auto tex_xxxz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 31);

            auto tey_xxxz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 31);

            auto tez_xxxz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 31);

            auto tex_xxxz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 32);

            auto tey_xxxz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 32);

            // Batch of Integrals (49,98)

            #pragma omp simd aligned(fx, pa_x, pc_x, ta_xxy_xxxz_1, ta_xxy_xxyy_1, ta_xxy_xxyz_1, \
                                         ta_xxy_xxzz_1, ta_xxy_xyyy_1, ta_xxy_xyyz_1, ta_xxy_xyzz_1, ta_xxy_xzzz_1, \
                                         ta_xxy_yyyy_1, ta_xxy_yyyz_1, ta_xxy_yyzz_1, ta_xxy_yzzz_1, ta_xxy_zzzz_1, \
                                         ta_xxz_xxxx_1, ta_xxz_xxxy_1, ta_xxz_xxxz_1, tex_xxxy_xxxz_0, tex_xxxy_xxyy_0, \
                                         tex_xxxy_xxyz_0, tex_xxxy_xxzz_0, tex_xxxy_xyyy_0, tex_xxxy_xyyz_0, tex_xxxy_xyzz_0, \
                                         tex_xxxy_xzzz_0, tex_xxxy_yyyy_0, tex_xxxy_yyyz_0, tex_xxxy_yyzz_0, tex_xxxy_yzzz_0, \
                                         tex_xxxy_zzzz_0, tex_xxxz_xxxx_0, tex_xxxz_xxxy_0, tex_xxxz_xxxz_0, tex_xxy_xxxz_0, \
                                         tex_xxy_xxxz_1, tex_xxy_xxyy_0, tex_xxy_xxyy_1, tex_xxy_xxyz_0, tex_xxy_xxyz_1, \
                                         tex_xxy_xxz_0, tex_xxy_xxz_1, tex_xxy_xxzz_0, tex_xxy_xxzz_1, tex_xxy_xyy_0, \
                                         tex_xxy_xyy_1, tex_xxy_xyyy_0, tex_xxy_xyyy_1, tex_xxy_xyyz_0, tex_xxy_xyyz_1, \
                                         tex_xxy_xyz_0, tex_xxy_xyz_1, tex_xxy_xyzz_0, tex_xxy_xyzz_1, tex_xxy_xzz_0, \
                                         tex_xxy_xzz_1, tex_xxy_xzzz_0, tex_xxy_xzzz_1, tex_xxy_yyy_0, tex_xxy_yyy_1, \
                                         tex_xxy_yyyy_0, tex_xxy_yyyy_1, tex_xxy_yyyz_0, tex_xxy_yyyz_1, tex_xxy_yyz_0, \
                                         tex_xxy_yyz_1, tex_xxy_yyzz_0, tex_xxy_yyzz_1, tex_xxy_yzz_0, tex_xxy_yzz_1, \
                                         tex_xxy_yzzz_0, tex_xxy_yzzz_1, tex_xxy_zzz_0, tex_xxy_zzz_1, tex_xxy_zzzz_0, \
                                         tex_xxy_zzzz_1, tex_xxz_xxx_0, tex_xxz_xxx_1, tex_xxz_xxxx_0, tex_xxz_xxxx_1, \
                                         tex_xxz_xxxy_0, tex_xxz_xxxy_1, tex_xxz_xxxz_0, tex_xxz_xxxz_1, tex_xxz_xxy_0, \
                                         tex_xxz_xxy_1, tex_xxz_xxz_0, tex_xxz_xxz_1, tex_xy_xxxz_0, tex_xy_xxxz_1, \
                                         tex_xy_xxyy_0, tex_xy_xxyy_1, tex_xy_xxyz_0, tex_xy_xxyz_1, tex_xy_xxzz_0, \
                                         tex_xy_xxzz_1, tex_xy_xyyy_0, tex_xy_xyyy_1, tex_xy_xyyz_0, tex_xy_xyyz_1, \
                                         tex_xy_xyzz_0, tex_xy_xyzz_1, tex_xy_xzzz_0, tex_xy_xzzz_1, tex_xy_yyyy_0, \
                                         tex_xy_yyyy_1, tex_xy_yyyz_0, tex_xy_yyyz_1, tex_xy_yyzz_0, tex_xy_yyzz_1, \
                                         tex_xy_yzzz_0, tex_xy_yzzz_1, tex_xy_zzzz_0, tex_xy_zzzz_1, tex_xz_xxxx_0, \
                                         tex_xz_xxxx_1, tex_xz_xxxy_0, tex_xz_xxxy_1, tex_xz_xxxz_0, tex_xz_xxxz_1, \
                                         tey_xxxy_xxxy_0, tey_xxxy_xxxz_0, tey_xxxy_xxyy_0, tey_xxxy_xxyz_0, tey_xxxy_xxzz_0, \
                                         tey_xxxy_xyyy_0, tey_xxxy_xyyz_0, tey_xxxy_xyzz_0, tey_xxxy_xzzz_0, tey_xxxy_yyyy_0, \
                                         tey_xxxy_yyyz_0, tey_xxxy_yyzz_0, tey_xxxy_yzzz_0, tey_xxxy_zzzz_0, tey_xxxz_xxxx_0, \
                                         tey_xxxz_xxxy_0, tey_xxxz_xxxz_0, tey_xxy_xxxy_0, tey_xxy_xxxy_1, tey_xxy_xxxz_0, \
                                         tey_xxy_xxxz_1, tey_xxy_xxy_0, tey_xxy_xxy_1, tey_xxy_xxyy_0, tey_xxy_xxyy_1, \
                                         tey_xxy_xxyz_0, tey_xxy_xxyz_1, tey_xxy_xxz_0, tey_xxy_xxz_1, tey_xxy_xxzz_0, \
                                         tey_xxy_xxzz_1, tey_xxy_xyy_0, tey_xxy_xyy_1, tey_xxy_xyyy_0, tey_xxy_xyyy_1, \
                                         tey_xxy_xyyz_0, tey_xxy_xyyz_1, tey_xxy_xyz_0, tey_xxy_xyz_1, tey_xxy_xyzz_0, \
                                         tey_xxy_xyzz_1, tey_xxy_xzz_0, tey_xxy_xzz_1, tey_xxy_xzzz_0, tey_xxy_xzzz_1, \
                                         tey_xxy_yyy_0, tey_xxy_yyy_1, tey_xxy_yyyy_0, tey_xxy_yyyy_1, tey_xxy_yyyz_0, \
                                         tey_xxy_yyyz_1, tey_xxy_yyz_0, tey_xxy_yyz_1, tey_xxy_yyzz_0, tey_xxy_yyzz_1, \
                                         tey_xxy_yzz_0, tey_xxy_yzz_1, tey_xxy_yzzz_0, tey_xxy_yzzz_1, tey_xxy_zzz_0, \
                                         tey_xxy_zzz_1, tey_xxy_zzzz_0, tey_xxy_zzzz_1, tey_xxz_xxx_0, tey_xxz_xxx_1, \
                                         tey_xxz_xxxx_0, tey_xxz_xxxx_1, tey_xxz_xxxy_0, tey_xxz_xxxy_1, tey_xxz_xxxz_0, \
                                         tey_xxz_xxxz_1, tey_xxz_xxy_0, tey_xxz_xxy_1, tey_xxz_xxz_0, tey_xxz_xxz_1, \
                                         tey_xy_xxxy_0, tey_xy_xxxy_1, tey_xy_xxxz_0, tey_xy_xxxz_1, tey_xy_xxyy_0, \
                                         tey_xy_xxyy_1, tey_xy_xxyz_0, tey_xy_xxyz_1, tey_xy_xxzz_0, tey_xy_xxzz_1, \
                                         tey_xy_xyyy_0, tey_xy_xyyy_1, tey_xy_xyyz_0, tey_xy_xyyz_1, tey_xy_xyzz_0, \
                                         tey_xy_xyzz_1, tey_xy_xzzz_0, tey_xy_xzzz_1, tey_xy_yyyy_0, tey_xy_yyyy_1, \
                                         tey_xy_yyyz_0, tey_xy_yyyz_1, tey_xy_yyzz_0, tey_xy_yyzz_1, tey_xy_yzzz_0, \
                                         tey_xy_yzzz_1, tey_xy_zzzz_0, tey_xy_zzzz_1, tey_xz_xxxx_0, tey_xz_xxxx_1, \
                                         tey_xz_xxxy_0, tey_xz_xxxy_1, tey_xz_xxxz_0, tey_xz_xxxz_1, tez_xxxy_xxxy_0, \
                                         tez_xxxy_xxxz_0, tez_xxxy_xxyy_0, tez_xxxy_xxyz_0, tez_xxxy_xxzz_0, tez_xxxy_xyyy_0, \
                                         tez_xxxy_xyyz_0, tez_xxxy_xyzz_0, tez_xxxy_xzzz_0, tez_xxxy_yyyy_0, tez_xxxy_yyyz_0, \
                                         tez_xxxy_yyzz_0, tez_xxxy_yzzz_0, tez_xxxy_zzzz_0, tez_xxxz_xxxx_0, tez_xxxz_xxxy_0, \
                                         tez_xxy_xxxy_0, tez_xxy_xxxy_1, tez_xxy_xxxz_0, tez_xxy_xxxz_1, tez_xxy_xxy_0, \
                                         tez_xxy_xxy_1, tez_xxy_xxyy_0, tez_xxy_xxyy_1, tez_xxy_xxyz_0, tez_xxy_xxyz_1, \
                                         tez_xxy_xxz_0, tez_xxy_xxz_1, tez_xxy_xxzz_0, tez_xxy_xxzz_1, tez_xxy_xyy_0, \
                                         tez_xxy_xyy_1, tez_xxy_xyyy_0, tez_xxy_xyyy_1, tez_xxy_xyyz_0, tez_xxy_xyyz_1, \
                                         tez_xxy_xyz_0, tez_xxy_xyz_1, tez_xxy_xyzz_0, tez_xxy_xyzz_1, tez_xxy_xzz_0, \
                                         tez_xxy_xzz_1, tez_xxy_xzzz_0, tez_xxy_xzzz_1, tez_xxy_yyy_0, tez_xxy_yyy_1, \
                                         tez_xxy_yyyy_0, tez_xxy_yyyy_1, tez_xxy_yyyz_0, tez_xxy_yyyz_1, tez_xxy_yyz_0, \
                                         tez_xxy_yyz_1, tez_xxy_yyzz_0, tez_xxy_yyzz_1, tez_xxy_yzz_0, tez_xxy_yzz_1, \
                                         tez_xxy_yzzz_0, tez_xxy_yzzz_1, tez_xxy_zzz_0, tez_xxy_zzz_1, tez_xxy_zzzz_0, \
                                         tez_xxy_zzzz_1, tez_xxz_xxx_0, tez_xxz_xxx_1, tez_xxz_xxxx_0, tez_xxz_xxxx_1, \
                                         tez_xxz_xxxy_0, tez_xxz_xxxy_1, tez_xxz_xxy_0, tez_xxz_xxy_1, tez_xy_xxxy_0, \
                                         tez_xy_xxxy_1, tez_xy_xxxz_0, tez_xy_xxxz_1, tez_xy_xxyy_0, tez_xy_xxyy_1, \
                                         tez_xy_xxyz_0, tez_xy_xxyz_1, tez_xy_xxzz_0, tez_xy_xxzz_1, tez_xy_xyyy_0, \
                                         tez_xy_xyyy_1, tez_xy_xyyz_0, tez_xy_xyyz_1, tez_xy_xyzz_0, tez_xy_xyzz_1, \
                                         tez_xy_xzzz_0, tez_xy_xzzz_1, tez_xy_yyyy_0, tez_xy_yyyy_1, tez_xy_yyyz_0, \
                                         tez_xy_yyyz_1, tez_xy_yyzz_0, tez_xy_yyzz_1, tez_xy_yzzz_0, tez_xy_yzzz_1, \
                                         tez_xy_zzzz_0, tez_xy_zzzz_1, tez_xz_xxxx_0, tez_xz_xxxx_1, tez_xz_xxxy_0, \
                                         tez_xz_xxxy_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tey_xxxy_xxxy_0[j] = pa_x[j] * tey_xxy_xxxy_0[j] - pc_x[j] * tey_xxy_xxxy_1[j] + fl1_fx * tey_xy_xxxy_0[j] -
                                     fl1_fx * tey_xy_xxxy_1[j] + 1.5 * fl1_fx * tey_xxy_xxy_0[j] - 1.5 * fl1_fx * tey_xxy_xxy_1[j];

                tez_xxxy_xxxy_0[j] = pa_x[j] * tez_xxy_xxxy_0[j] - pc_x[j] * tez_xxy_xxxy_1[j] + fl1_fx * tez_xy_xxxy_0[j] -
                                     fl1_fx * tez_xy_xxxy_1[j] + 1.5 * fl1_fx * tez_xxy_xxy_0[j] - 1.5 * fl1_fx * tez_xxy_xxy_1[j];

                tex_xxxy_xxxz_0[j] = pa_x[j] * tex_xxy_xxxz_0[j] - pc_x[j] * tex_xxy_xxxz_1[j] + fl1_fx * tex_xy_xxxz_0[j] -
                                     fl1_fx * tex_xy_xxxz_1[j] + 1.5 * fl1_fx * tex_xxy_xxz_0[j] - 1.5 * fl1_fx * tex_xxy_xxz_1[j] + ta_xxy_xxxz_1[j];

                tey_xxxy_xxxz_0[j] = pa_x[j] * tey_xxy_xxxz_0[j] - pc_x[j] * tey_xxy_xxxz_1[j] + fl1_fx * tey_xy_xxxz_0[j] -
                                     fl1_fx * tey_xy_xxxz_1[j] + 1.5 * fl1_fx * tey_xxy_xxz_0[j] - 1.5 * fl1_fx * tey_xxy_xxz_1[j];

                tez_xxxy_xxxz_0[j] = pa_x[j] * tez_xxy_xxxz_0[j] - pc_x[j] * tez_xxy_xxxz_1[j] + fl1_fx * tez_xy_xxxz_0[j] -
                                     fl1_fx * tez_xy_xxxz_1[j] + 1.5 * fl1_fx * tez_xxy_xxz_0[j] - 1.5 * fl1_fx * tez_xxy_xxz_1[j];

                tex_xxxy_xxyy_0[j] = pa_x[j] * tex_xxy_xxyy_0[j] - pc_x[j] * tex_xxy_xxyy_1[j] + fl1_fx * tex_xy_xxyy_0[j] -
                                     fl1_fx * tex_xy_xxyy_1[j] + fl1_fx * tex_xxy_xyy_0[j] - fl1_fx * tex_xxy_xyy_1[j] + ta_xxy_xxyy_1[j];

                tey_xxxy_xxyy_0[j] = pa_x[j] * tey_xxy_xxyy_0[j] - pc_x[j] * tey_xxy_xxyy_1[j] + fl1_fx * tey_xy_xxyy_0[j] -
                                     fl1_fx * tey_xy_xxyy_1[j] + fl1_fx * tey_xxy_xyy_0[j] - fl1_fx * tey_xxy_xyy_1[j];

                tez_xxxy_xxyy_0[j] = pa_x[j] * tez_xxy_xxyy_0[j] - pc_x[j] * tez_xxy_xxyy_1[j] + fl1_fx * tez_xy_xxyy_0[j] -
                                     fl1_fx * tez_xy_xxyy_1[j] + fl1_fx * tez_xxy_xyy_0[j] - fl1_fx * tez_xxy_xyy_1[j];

                tex_xxxy_xxyz_0[j] = pa_x[j] * tex_xxy_xxyz_0[j] - pc_x[j] * tex_xxy_xxyz_1[j] + fl1_fx * tex_xy_xxyz_0[j] -
                                     fl1_fx * tex_xy_xxyz_1[j] + fl1_fx * tex_xxy_xyz_0[j] - fl1_fx * tex_xxy_xyz_1[j] + ta_xxy_xxyz_1[j];

                tey_xxxy_xxyz_0[j] = pa_x[j] * tey_xxy_xxyz_0[j] - pc_x[j] * tey_xxy_xxyz_1[j] + fl1_fx * tey_xy_xxyz_0[j] -
                                     fl1_fx * tey_xy_xxyz_1[j] + fl1_fx * tey_xxy_xyz_0[j] - fl1_fx * tey_xxy_xyz_1[j];

                tez_xxxy_xxyz_0[j] = pa_x[j] * tez_xxy_xxyz_0[j] - pc_x[j] * tez_xxy_xxyz_1[j] + fl1_fx * tez_xy_xxyz_0[j] -
                                     fl1_fx * tez_xy_xxyz_1[j] + fl1_fx * tez_xxy_xyz_0[j] - fl1_fx * tez_xxy_xyz_1[j];

                tex_xxxy_xxzz_0[j] = pa_x[j] * tex_xxy_xxzz_0[j] - pc_x[j] * tex_xxy_xxzz_1[j] + fl1_fx * tex_xy_xxzz_0[j] -
                                     fl1_fx * tex_xy_xxzz_1[j] + fl1_fx * tex_xxy_xzz_0[j] - fl1_fx * tex_xxy_xzz_1[j] + ta_xxy_xxzz_1[j];

                tey_xxxy_xxzz_0[j] = pa_x[j] * tey_xxy_xxzz_0[j] - pc_x[j] * tey_xxy_xxzz_1[j] + fl1_fx * tey_xy_xxzz_0[j] -
                                     fl1_fx * tey_xy_xxzz_1[j] + fl1_fx * tey_xxy_xzz_0[j] - fl1_fx * tey_xxy_xzz_1[j];

                tez_xxxy_xxzz_0[j] = pa_x[j] * tez_xxy_xxzz_0[j] - pc_x[j] * tez_xxy_xxzz_1[j] + fl1_fx * tez_xy_xxzz_0[j] -
                                     fl1_fx * tez_xy_xxzz_1[j] + fl1_fx * tez_xxy_xzz_0[j] - fl1_fx * tez_xxy_xzz_1[j];

                tex_xxxy_xyyy_0[j] = pa_x[j] * tex_xxy_xyyy_0[j] - pc_x[j] * tex_xxy_xyyy_1[j] + fl1_fx * tex_xy_xyyy_0[j] -
                                     fl1_fx * tex_xy_xyyy_1[j] + 0.5 * fl1_fx * tex_xxy_yyy_0[j] - 0.5 * fl1_fx * tex_xxy_yyy_1[j] + ta_xxy_xyyy_1[j];

                tey_xxxy_xyyy_0[j] = pa_x[j] * tey_xxy_xyyy_0[j] - pc_x[j] * tey_xxy_xyyy_1[j] + fl1_fx * tey_xy_xyyy_0[j] -
                                     fl1_fx * tey_xy_xyyy_1[j] + 0.5 * fl1_fx * tey_xxy_yyy_0[j] - 0.5 * fl1_fx * tey_xxy_yyy_1[j];

                tez_xxxy_xyyy_0[j] = pa_x[j] * tez_xxy_xyyy_0[j] - pc_x[j] * tez_xxy_xyyy_1[j] + fl1_fx * tez_xy_xyyy_0[j] -
                                     fl1_fx * tez_xy_xyyy_1[j] + 0.5 * fl1_fx * tez_xxy_yyy_0[j] - 0.5 * fl1_fx * tez_xxy_yyy_1[j];

                tex_xxxy_xyyz_0[j] = pa_x[j] * tex_xxy_xyyz_0[j] - pc_x[j] * tex_xxy_xyyz_1[j] + fl1_fx * tex_xy_xyyz_0[j] -
                                     fl1_fx * tex_xy_xyyz_1[j] + 0.5 * fl1_fx * tex_xxy_yyz_0[j] - 0.5 * fl1_fx * tex_xxy_yyz_1[j] + ta_xxy_xyyz_1[j];

                tey_xxxy_xyyz_0[j] = pa_x[j] * tey_xxy_xyyz_0[j] - pc_x[j] * tey_xxy_xyyz_1[j] + fl1_fx * tey_xy_xyyz_0[j] -
                                     fl1_fx * tey_xy_xyyz_1[j] + 0.5 * fl1_fx * tey_xxy_yyz_0[j] - 0.5 * fl1_fx * tey_xxy_yyz_1[j];

                tez_xxxy_xyyz_0[j] = pa_x[j] * tez_xxy_xyyz_0[j] - pc_x[j] * tez_xxy_xyyz_1[j] + fl1_fx * tez_xy_xyyz_0[j] -
                                     fl1_fx * tez_xy_xyyz_1[j] + 0.5 * fl1_fx * tez_xxy_yyz_0[j] - 0.5 * fl1_fx * tez_xxy_yyz_1[j];

                tex_xxxy_xyzz_0[j] = pa_x[j] * tex_xxy_xyzz_0[j] - pc_x[j] * tex_xxy_xyzz_1[j] + fl1_fx * tex_xy_xyzz_0[j] -
                                     fl1_fx * tex_xy_xyzz_1[j] + 0.5 * fl1_fx * tex_xxy_yzz_0[j] - 0.5 * fl1_fx * tex_xxy_yzz_1[j] + ta_xxy_xyzz_1[j];

                tey_xxxy_xyzz_0[j] = pa_x[j] * tey_xxy_xyzz_0[j] - pc_x[j] * tey_xxy_xyzz_1[j] + fl1_fx * tey_xy_xyzz_0[j] -
                                     fl1_fx * tey_xy_xyzz_1[j] + 0.5 * fl1_fx * tey_xxy_yzz_0[j] - 0.5 * fl1_fx * tey_xxy_yzz_1[j];

                tez_xxxy_xyzz_0[j] = pa_x[j] * tez_xxy_xyzz_0[j] - pc_x[j] * tez_xxy_xyzz_1[j] + fl1_fx * tez_xy_xyzz_0[j] -
                                     fl1_fx * tez_xy_xyzz_1[j] + 0.5 * fl1_fx * tez_xxy_yzz_0[j] - 0.5 * fl1_fx * tez_xxy_yzz_1[j];

                tex_xxxy_xzzz_0[j] = pa_x[j] * tex_xxy_xzzz_0[j] - pc_x[j] * tex_xxy_xzzz_1[j] + fl1_fx * tex_xy_xzzz_0[j] -
                                     fl1_fx * tex_xy_xzzz_1[j] + 0.5 * fl1_fx * tex_xxy_zzz_0[j] - 0.5 * fl1_fx * tex_xxy_zzz_1[j] + ta_xxy_xzzz_1[j];

                tey_xxxy_xzzz_0[j] = pa_x[j] * tey_xxy_xzzz_0[j] - pc_x[j] * tey_xxy_xzzz_1[j] + fl1_fx * tey_xy_xzzz_0[j] -
                                     fl1_fx * tey_xy_xzzz_1[j] + 0.5 * fl1_fx * tey_xxy_zzz_0[j] - 0.5 * fl1_fx * tey_xxy_zzz_1[j];

                tez_xxxy_xzzz_0[j] = pa_x[j] * tez_xxy_xzzz_0[j] - pc_x[j] * tez_xxy_xzzz_1[j] + fl1_fx * tez_xy_xzzz_0[j] -
                                     fl1_fx * tez_xy_xzzz_1[j] + 0.5 * fl1_fx * tez_xxy_zzz_0[j] - 0.5 * fl1_fx * tez_xxy_zzz_1[j];

                tex_xxxy_yyyy_0[j] = pa_x[j] * tex_xxy_yyyy_0[j] - pc_x[j] * tex_xxy_yyyy_1[j] + fl1_fx * tex_xy_yyyy_0[j] -
                                     fl1_fx * tex_xy_yyyy_1[j] + ta_xxy_yyyy_1[j];

                tey_xxxy_yyyy_0[j] =
                    pa_x[j] * tey_xxy_yyyy_0[j] - pc_x[j] * tey_xxy_yyyy_1[j] + fl1_fx * tey_xy_yyyy_0[j] - fl1_fx * tey_xy_yyyy_1[j];

                tez_xxxy_yyyy_0[j] =
                    pa_x[j] * tez_xxy_yyyy_0[j] - pc_x[j] * tez_xxy_yyyy_1[j] + fl1_fx * tez_xy_yyyy_0[j] - fl1_fx * tez_xy_yyyy_1[j];

                tex_xxxy_yyyz_0[j] = pa_x[j] * tex_xxy_yyyz_0[j] - pc_x[j] * tex_xxy_yyyz_1[j] + fl1_fx * tex_xy_yyyz_0[j] -
                                     fl1_fx * tex_xy_yyyz_1[j] + ta_xxy_yyyz_1[j];

                tey_xxxy_yyyz_0[j] =
                    pa_x[j] * tey_xxy_yyyz_0[j] - pc_x[j] * tey_xxy_yyyz_1[j] + fl1_fx * tey_xy_yyyz_0[j] - fl1_fx * tey_xy_yyyz_1[j];

                tez_xxxy_yyyz_0[j] =
                    pa_x[j] * tez_xxy_yyyz_0[j] - pc_x[j] * tez_xxy_yyyz_1[j] + fl1_fx * tez_xy_yyyz_0[j] - fl1_fx * tez_xy_yyyz_1[j];

                tex_xxxy_yyzz_0[j] = pa_x[j] * tex_xxy_yyzz_0[j] - pc_x[j] * tex_xxy_yyzz_1[j] + fl1_fx * tex_xy_yyzz_0[j] -
                                     fl1_fx * tex_xy_yyzz_1[j] + ta_xxy_yyzz_1[j];

                tey_xxxy_yyzz_0[j] =
                    pa_x[j] * tey_xxy_yyzz_0[j] - pc_x[j] * tey_xxy_yyzz_1[j] + fl1_fx * tey_xy_yyzz_0[j] - fl1_fx * tey_xy_yyzz_1[j];

                tez_xxxy_yyzz_0[j] =
                    pa_x[j] * tez_xxy_yyzz_0[j] - pc_x[j] * tez_xxy_yyzz_1[j] + fl1_fx * tez_xy_yyzz_0[j] - fl1_fx * tez_xy_yyzz_1[j];

                tex_xxxy_yzzz_0[j] = pa_x[j] * tex_xxy_yzzz_0[j] - pc_x[j] * tex_xxy_yzzz_1[j] + fl1_fx * tex_xy_yzzz_0[j] -
                                     fl1_fx * tex_xy_yzzz_1[j] + ta_xxy_yzzz_1[j];

                tey_xxxy_yzzz_0[j] =
                    pa_x[j] * tey_xxy_yzzz_0[j] - pc_x[j] * tey_xxy_yzzz_1[j] + fl1_fx * tey_xy_yzzz_0[j] - fl1_fx * tey_xy_yzzz_1[j];

                tez_xxxy_yzzz_0[j] =
                    pa_x[j] * tez_xxy_yzzz_0[j] - pc_x[j] * tez_xxy_yzzz_1[j] + fl1_fx * tez_xy_yzzz_0[j] - fl1_fx * tez_xy_yzzz_1[j];

                tex_xxxy_zzzz_0[j] = pa_x[j] * tex_xxy_zzzz_0[j] - pc_x[j] * tex_xxy_zzzz_1[j] + fl1_fx * tex_xy_zzzz_0[j] -
                                     fl1_fx * tex_xy_zzzz_1[j] + ta_xxy_zzzz_1[j];

                tey_xxxy_zzzz_0[j] =
                    pa_x[j] * tey_xxy_zzzz_0[j] - pc_x[j] * tey_xxy_zzzz_1[j] + fl1_fx * tey_xy_zzzz_0[j] - fl1_fx * tey_xy_zzzz_1[j];

                tez_xxxy_zzzz_0[j] =
                    pa_x[j] * tez_xxy_zzzz_0[j] - pc_x[j] * tez_xxy_zzzz_1[j] + fl1_fx * tez_xy_zzzz_0[j] - fl1_fx * tez_xy_zzzz_1[j];

                tex_xxxz_xxxx_0[j] = pa_x[j] * tex_xxz_xxxx_0[j] - pc_x[j] * tex_xxz_xxxx_1[j] + fl1_fx * tex_xz_xxxx_0[j] -
                                     fl1_fx * tex_xz_xxxx_1[j] + 2.0 * fl1_fx * tex_xxz_xxx_0[j] - 2.0 * fl1_fx * tex_xxz_xxx_1[j] + ta_xxz_xxxx_1[j];

                tey_xxxz_xxxx_0[j] = pa_x[j] * tey_xxz_xxxx_0[j] - pc_x[j] * tey_xxz_xxxx_1[j] + fl1_fx * tey_xz_xxxx_0[j] -
                                     fl1_fx * tey_xz_xxxx_1[j] + 2.0 * fl1_fx * tey_xxz_xxx_0[j] - 2.0 * fl1_fx * tey_xxz_xxx_1[j];

                tez_xxxz_xxxx_0[j] = pa_x[j] * tez_xxz_xxxx_0[j] - pc_x[j] * tez_xxz_xxxx_1[j] + fl1_fx * tez_xz_xxxx_0[j] -
                                     fl1_fx * tez_xz_xxxx_1[j] + 2.0 * fl1_fx * tez_xxz_xxx_0[j] - 2.0 * fl1_fx * tez_xxz_xxx_1[j];

                tex_xxxz_xxxy_0[j] = pa_x[j] * tex_xxz_xxxy_0[j] - pc_x[j] * tex_xxz_xxxy_1[j] + fl1_fx * tex_xz_xxxy_0[j] -
                                     fl1_fx * tex_xz_xxxy_1[j] + 1.5 * fl1_fx * tex_xxz_xxy_0[j] - 1.5 * fl1_fx * tex_xxz_xxy_1[j] + ta_xxz_xxxy_1[j];

                tey_xxxz_xxxy_0[j] = pa_x[j] * tey_xxz_xxxy_0[j] - pc_x[j] * tey_xxz_xxxy_1[j] + fl1_fx * tey_xz_xxxy_0[j] -
                                     fl1_fx * tey_xz_xxxy_1[j] + 1.5 * fl1_fx * tey_xxz_xxy_0[j] - 1.5 * fl1_fx * tey_xxz_xxy_1[j];

                tez_xxxz_xxxy_0[j] = pa_x[j] * tez_xxz_xxxy_0[j] - pc_x[j] * tez_xxz_xxxy_1[j] + fl1_fx * tez_xz_xxxy_0[j] -
                                     fl1_fx * tez_xz_xxxy_1[j] + 1.5 * fl1_fx * tez_xxz_xxy_0[j] - 1.5 * fl1_fx * tez_xxz_xxy_1[j];

                tex_xxxz_xxxz_0[j] = pa_x[j] * tex_xxz_xxxz_0[j] - pc_x[j] * tex_xxz_xxxz_1[j] + fl1_fx * tex_xz_xxxz_0[j] -
                                     fl1_fx * tex_xz_xxxz_1[j] + 1.5 * fl1_fx * tex_xxz_xxz_0[j] - 1.5 * fl1_fx * tex_xxz_xxz_1[j] + ta_xxz_xxxz_1[j];

                tey_xxxz_xxxz_0[j] = pa_x[j] * tey_xxz_xxxz_0[j] - pc_x[j] * tey_xxz_xxxz_1[j] + fl1_fx * tey_xz_xxxz_0[j] -
                                     fl1_fx * tey_xz_xxxz_1[j] + 1.5 * fl1_fx * tey_xxz_xxz_0[j] - 1.5 * fl1_fx * tey_xxz_xxz_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForGG_98_147(CMemBlock2D<double>&       primBuffer,
                              const CRecursionMap&       recursionMap,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CMemBlock2D<double>& pcDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto)
{
    // Batch of Integrals (98,147)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_4_4_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_2_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tez_xxz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 32);

            auto tex_xxz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 33);

            auto tey_xxz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 33);

            auto tez_xxz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 33);

            auto tex_xxz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 34);

            auto tey_xxz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 34);

            auto tez_xxz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 34);

            auto tex_xxz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 35);

            auto tey_xxz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 35);

            auto tez_xxz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 35);

            auto tex_xxz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 36);

            auto tey_xxz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 36);

            auto tez_xxz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 36);

            auto tex_xxz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 37);

            auto tey_xxz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 37);

            auto tez_xxz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 37);

            auto tex_xxz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 38);

            auto tey_xxz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 38);

            auto tez_xxz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 38);

            auto tex_xxz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 39);

            auto tey_xxz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 39);

            auto tez_xxz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 39);

            auto tex_xxz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 40);

            auto tey_xxz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 40);

            auto tez_xxz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 40);

            auto tex_xxz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 41);

            auto tey_xxz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 41);

            auto tez_xxz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 41);

            auto tex_xxz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 42);

            auto tey_xxz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 42);

            auto tez_xxz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 42);

            auto tex_xxz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 43);

            auto tey_xxz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 43);

            auto tez_xxz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 43);

            auto tex_xxz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 44);

            auto tey_xxz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 44);

            auto tez_xxz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 44);

            auto tex_xyy_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 45);

            auto tey_xyy_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 45);

            auto tez_xyy_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 45);

            auto tex_xyy_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 46);

            auto tey_xyy_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 46);

            auto tez_xyy_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 46);

            auto tex_xyy_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 47);

            auto tey_xyy_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 47);

            auto tez_xyy_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 47);

            auto tex_xyy_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 48);

            auto tey_xyy_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 48);

            auto tez_xyy_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 48);

            auto tez_xxz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 32);

            auto tex_xxz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 33);

            auto tey_xxz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 33);

            auto tez_xxz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 33);

            auto tex_xxz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 34);

            auto tey_xxz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 34);

            auto tez_xxz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 34);

            auto tex_xxz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 35);

            auto tey_xxz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 35);

            auto tez_xxz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 35);

            auto tex_xxz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 36);

            auto tey_xxz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 36);

            auto tez_xxz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 36);

            auto tex_xxz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 37);

            auto tey_xxz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 37);

            auto tez_xxz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 37);

            auto tex_xxz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 38);

            auto tey_xxz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 38);

            auto tez_xxz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 38);

            auto tex_xxz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 39);

            auto tey_xxz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 39);

            auto tez_xxz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 39);

            auto tex_xxz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 40);

            auto tey_xxz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 40);

            auto tez_xxz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 40);

            auto tex_xxz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 41);

            auto tey_xxz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 41);

            auto tez_xxz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 41);

            auto tex_xxz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 42);

            auto tey_xxz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 42);

            auto tez_xxz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 42);

            auto tex_xxz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 43);

            auto tey_xxz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 43);

            auto tez_xxz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 43);

            auto tex_xxz_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 44);

            auto tey_xxz_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 44);

            auto tez_xxz_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 44);

            auto tex_xyy_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 45);

            auto tey_xyy_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 45);

            auto tez_xyy_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 45);

            auto tex_xyy_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 46);

            auto tey_xyy_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 46);

            auto tez_xyy_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 46);

            auto tex_xyy_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 47);

            auto tey_xyy_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 47);

            auto tez_xyy_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 47);

            auto tex_xyy_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 48);

            auto tey_xyy_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 48);

            auto tez_xyy_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 48);

            auto tez_xz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 32);

            auto tex_xz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 33);

            auto tey_xz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 33);

            auto tez_xz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 33);

            auto tex_xz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 34);

            auto tey_xz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 34);

            auto tez_xz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 34);

            auto tex_xz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 35);

            auto tey_xz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 35);

            auto tez_xz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 35);

            auto tex_xz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 36);

            auto tey_xz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 36);

            auto tez_xz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 36);

            auto tex_xz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 37);

            auto tey_xz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 37);

            auto tez_xz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 37);

            auto tex_xz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 38);

            auto tey_xz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 38);

            auto tez_xz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 38);

            auto tex_xz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 39);

            auto tey_xz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 39);

            auto tez_xz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 39);

            auto tex_xz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 40);

            auto tey_xz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 40);

            auto tez_xz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 40);

            auto tex_xz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 41);

            auto tey_xz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 41);

            auto tez_xz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 41);

            auto tex_xz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 42);

            auto tey_xz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 42);

            auto tez_xz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 42);

            auto tex_xz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 43);

            auto tey_xz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 43);

            auto tez_xz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 43);

            auto tex_xz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 44);

            auto tey_xz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 44);

            auto tez_xz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 44);

            auto tex_yy_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 45);

            auto tey_yy_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 45);

            auto tez_yy_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 45);

            auto tex_yy_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 46);

            auto tey_yy_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 46);

            auto tez_yy_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 46);

            auto tex_yy_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 47);

            auto tey_yy_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 47);

            auto tez_yy_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 47);

            auto tex_yy_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 48);

            auto tey_yy_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 48);

            auto tez_yy_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 48);

            auto tez_xz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 32);

            auto tex_xz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 33);

            auto tey_xz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 33);

            auto tez_xz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 33);

            auto tex_xz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 34);

            auto tey_xz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 34);

            auto tez_xz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 34);

            auto tex_xz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 35);

            auto tey_xz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 35);

            auto tez_xz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 35);

            auto tex_xz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 36);

            auto tey_xz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 36);

            auto tez_xz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 36);

            auto tex_xz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 37);

            auto tey_xz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 37);

            auto tez_xz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 37);

            auto tex_xz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 38);

            auto tey_xz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 38);

            auto tez_xz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 38);

            auto tex_xz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 39);

            auto tey_xz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 39);

            auto tez_xz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 39);

            auto tex_xz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 40);

            auto tey_xz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 40);

            auto tez_xz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 40);

            auto tex_xz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 41);

            auto tey_xz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 41);

            auto tez_xz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 41);

            auto tex_xz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 42);

            auto tey_xz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 42);

            auto tez_xz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 42);

            auto tex_xz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 43);

            auto tey_xz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 43);

            auto tez_xz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 43);

            auto tex_xz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 44);

            auto tey_xz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 44);

            auto tez_xz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 44);

            auto tex_yy_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 45);

            auto tey_yy_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 45);

            auto tez_yy_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 45);

            auto tex_yy_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 46);

            auto tey_yy_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 46);

            auto tez_yy_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 46);

            auto tex_yy_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 47);

            auto tey_yy_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 47);

            auto tez_yy_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 47);

            auto tex_yy_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 48);

            auto tey_yy_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 48);

            auto tez_yy_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 48);

            auto tez_xxz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 22);

            auto tex_xxz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 23);

            auto tey_xxz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 23);

            auto tez_xxz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 23);

            auto tex_xxz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 24);

            auto tey_xxz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 24);

            auto tez_xxz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 24);

            auto tex_xxz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 25);

            auto tey_xxz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 25);

            auto tez_xxz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 25);

            auto tex_xxz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 26);

            auto tey_xxz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 26);

            auto tez_xxz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 26);

            auto tex_xxz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 27);

            auto tey_xxz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 27);

            auto tez_xxz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 27);

            auto tex_xxz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 28);

            auto tey_xxz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 28);

            auto tez_xxz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 28);

            auto tex_xxz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 29);

            auto tey_xxz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 29);

            auto tez_xxz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 29);

            auto tex_xyy_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 30);

            auto tey_xyy_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 30);

            auto tez_xyy_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 30);

            auto tex_xyy_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 31);

            auto tey_xyy_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 31);

            auto tez_xyy_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 31);

            auto tex_xyy_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 32);

            auto tey_xyy_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 32);

            auto tez_xyy_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 32);

            auto tex_xyy_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 33);

            auto tey_xyy_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 33);

            auto tez_xyy_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 33);

            auto tez_xxz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 22);

            auto tex_xxz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 23);

            auto tey_xxz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 23);

            auto tez_xxz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 23);

            auto tex_xxz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 24);

            auto tey_xxz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 24);

            auto tez_xxz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 24);

            auto tex_xxz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 25);

            auto tey_xxz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 25);

            auto tez_xxz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 25);

            auto tex_xxz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 26);

            auto tey_xxz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 26);

            auto tez_xxz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 26);

            auto tex_xxz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 27);

            auto tey_xxz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 27);

            auto tez_xxz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 27);

            auto tex_xxz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 28);

            auto tey_xxz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 28);

            auto tez_xxz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 28);

            auto tex_xxz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 29);

            auto tey_xxz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 29);

            auto tez_xxz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 29);

            auto tex_xyy_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 30);

            auto tey_xyy_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 30);

            auto tez_xyy_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 30);

            auto tex_xyy_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 31);

            auto tey_xyy_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 31);

            auto tez_xyy_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 31);

            auto tex_xyy_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 32);

            auto tey_xyy_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 32);

            auto tez_xyy_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 32);

            auto tex_xyy_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 33);

            auto tey_xyy_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 33);

            auto tez_xyy_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 33);

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

            auto ta_xyy_xxxx_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 45);

            auto ta_xyy_xxxy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 46);

            auto ta_xyy_xxxz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 47);

            auto ta_xyy_xxyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 48);

            // set up pointers to integrals

            auto tez_xxxz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 32);

            auto tex_xxxz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 33);

            auto tey_xxxz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 33);

            auto tez_xxxz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 33);

            auto tex_xxxz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 34);

            auto tey_xxxz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 34);

            auto tez_xxxz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 34);

            auto tex_xxxz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 35);

            auto tey_xxxz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 35);

            auto tez_xxxz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 35);

            auto tex_xxxz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 36);

            auto tey_xxxz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 36);

            auto tez_xxxz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 36);

            auto tex_xxxz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 37);

            auto tey_xxxz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 37);

            auto tez_xxxz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 37);

            auto tex_xxxz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 38);

            auto tey_xxxz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 38);

            auto tez_xxxz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 38);

            auto tex_xxxz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 39);

            auto tey_xxxz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 39);

            auto tez_xxxz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 39);

            auto tex_xxxz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 40);

            auto tey_xxxz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 40);

            auto tez_xxxz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 40);

            auto tex_xxxz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 41);

            auto tey_xxxz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 41);

            auto tez_xxxz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 41);

            auto tex_xxxz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 42);

            auto tey_xxxz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 42);

            auto tez_xxxz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 42);

            auto tex_xxxz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 43);

            auto tey_xxxz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 43);

            auto tez_xxxz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 43);

            auto tex_xxxz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 44);

            auto tey_xxxz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 44);

            auto tez_xxxz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 44);

            auto tex_xxyy_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 45);

            auto tey_xxyy_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 45);

            auto tez_xxyy_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 45);

            auto tex_xxyy_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 46);

            auto tey_xxyy_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 46);

            auto tez_xxyy_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 46);

            auto tex_xxyy_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 47);

            auto tey_xxyy_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 47);

            auto tez_xxyy_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 47);

            auto tex_xxyy_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 48);

            auto tey_xxyy_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 48);

            auto tez_xxyy_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 48);

            // Batch of Integrals (98,147)

            #pragma omp simd aligned(fx, pa_x, pc_x, ta_xxz_xxyy_1, ta_xxz_xxyz_1, ta_xxz_xxzz_1, \
                                         ta_xxz_xyyy_1, ta_xxz_xyyz_1, ta_xxz_xyzz_1, ta_xxz_xzzz_1, ta_xxz_yyyy_1, \
                                         ta_xxz_yyyz_1, ta_xxz_yyzz_1, ta_xxz_yzzz_1, ta_xxz_zzzz_1, ta_xyy_xxxx_1, \
                                         ta_xyy_xxxy_1, ta_xyy_xxxz_1, ta_xyy_xxyy_1, tex_xxxz_xxyy_0, tex_xxxz_xxyz_0, \
                                         tex_xxxz_xxzz_0, tex_xxxz_xyyy_0, tex_xxxz_xyyz_0, tex_xxxz_xyzz_0, tex_xxxz_xzzz_0, \
                                         tex_xxxz_yyyy_0, tex_xxxz_yyyz_0, tex_xxxz_yyzz_0, tex_xxxz_yzzz_0, tex_xxxz_zzzz_0, \
                                         tex_xxyy_xxxx_0, tex_xxyy_xxxy_0, tex_xxyy_xxxz_0, tex_xxyy_xxyy_0, tex_xxz_xxyy_0, \
                                         tex_xxz_xxyy_1, tex_xxz_xxyz_0, tex_xxz_xxyz_1, tex_xxz_xxzz_0, tex_xxz_xxzz_1, \
                                         tex_xxz_xyy_0, tex_xxz_xyy_1, tex_xxz_xyyy_0, tex_xxz_xyyy_1, tex_xxz_xyyz_0, \
                                         tex_xxz_xyyz_1, tex_xxz_xyz_0, tex_xxz_xyz_1, tex_xxz_xyzz_0, tex_xxz_xyzz_1, \
                                         tex_xxz_xzz_0, tex_xxz_xzz_1, tex_xxz_xzzz_0, tex_xxz_xzzz_1, tex_xxz_yyy_0, \
                                         tex_xxz_yyy_1, tex_xxz_yyyy_0, tex_xxz_yyyy_1, tex_xxz_yyyz_0, tex_xxz_yyyz_1, \
                                         tex_xxz_yyz_0, tex_xxz_yyz_1, tex_xxz_yyzz_0, tex_xxz_yyzz_1, tex_xxz_yzz_0, \
                                         tex_xxz_yzz_1, tex_xxz_yzzz_0, tex_xxz_yzzz_1, tex_xxz_zzz_0, tex_xxz_zzz_1, \
                                         tex_xxz_zzzz_0, tex_xxz_zzzz_1, tex_xyy_xxx_0, tex_xyy_xxx_1, tex_xyy_xxxx_0, \
                                         tex_xyy_xxxx_1, tex_xyy_xxxy_0, tex_xyy_xxxy_1, tex_xyy_xxxz_0, tex_xyy_xxxz_1, \
                                         tex_xyy_xxy_0, tex_xyy_xxy_1, tex_xyy_xxyy_0, tex_xyy_xxyy_1, tex_xyy_xxz_0, \
                                         tex_xyy_xxz_1, tex_xyy_xyy_0, tex_xyy_xyy_1, tex_xz_xxyy_0, tex_xz_xxyy_1, \
                                         tex_xz_xxyz_0, tex_xz_xxyz_1, tex_xz_xxzz_0, tex_xz_xxzz_1, tex_xz_xyyy_0, \
                                         tex_xz_xyyy_1, tex_xz_xyyz_0, tex_xz_xyyz_1, tex_xz_xyzz_0, tex_xz_xyzz_1, \
                                         tex_xz_xzzz_0, tex_xz_xzzz_1, tex_xz_yyyy_0, tex_xz_yyyy_1, tex_xz_yyyz_0, \
                                         tex_xz_yyyz_1, tex_xz_yyzz_0, tex_xz_yyzz_1, tex_xz_yzzz_0, tex_xz_yzzz_1, \
                                         tex_xz_zzzz_0, tex_xz_zzzz_1, tex_yy_xxxx_0, tex_yy_xxxx_1, tex_yy_xxxy_0, \
                                         tex_yy_xxxy_1, tex_yy_xxxz_0, tex_yy_xxxz_1, tex_yy_xxyy_0, tex_yy_xxyy_1, \
                                         tey_xxxz_xxyy_0, tey_xxxz_xxyz_0, tey_xxxz_xxzz_0, tey_xxxz_xyyy_0, tey_xxxz_xyyz_0, \
                                         tey_xxxz_xyzz_0, tey_xxxz_xzzz_0, tey_xxxz_yyyy_0, tey_xxxz_yyyz_0, tey_xxxz_yyzz_0, \
                                         tey_xxxz_yzzz_0, tey_xxxz_zzzz_0, tey_xxyy_xxxx_0, tey_xxyy_xxxy_0, tey_xxyy_xxxz_0, \
                                         tey_xxyy_xxyy_0, tey_xxz_xxyy_0, tey_xxz_xxyy_1, tey_xxz_xxyz_0, tey_xxz_xxyz_1, \
                                         tey_xxz_xxzz_0, tey_xxz_xxzz_1, tey_xxz_xyy_0, tey_xxz_xyy_1, tey_xxz_xyyy_0, \
                                         tey_xxz_xyyy_1, tey_xxz_xyyz_0, tey_xxz_xyyz_1, tey_xxz_xyz_0, tey_xxz_xyz_1, \
                                         tey_xxz_xyzz_0, tey_xxz_xyzz_1, tey_xxz_xzz_0, tey_xxz_xzz_1, tey_xxz_xzzz_0, \
                                         tey_xxz_xzzz_1, tey_xxz_yyy_0, tey_xxz_yyy_1, tey_xxz_yyyy_0, tey_xxz_yyyy_1, \
                                         tey_xxz_yyyz_0, tey_xxz_yyyz_1, tey_xxz_yyz_0, tey_xxz_yyz_1, tey_xxz_yyzz_0, \
                                         tey_xxz_yyzz_1, tey_xxz_yzz_0, tey_xxz_yzz_1, tey_xxz_yzzz_0, tey_xxz_yzzz_1, \
                                         tey_xxz_zzz_0, tey_xxz_zzz_1, tey_xxz_zzzz_0, tey_xxz_zzzz_1, tey_xyy_xxx_0, \
                                         tey_xyy_xxx_1, tey_xyy_xxxx_0, tey_xyy_xxxx_1, tey_xyy_xxxy_0, tey_xyy_xxxy_1, \
                                         tey_xyy_xxxz_0, tey_xyy_xxxz_1, tey_xyy_xxy_0, tey_xyy_xxy_1, tey_xyy_xxyy_0, \
                                         tey_xyy_xxyy_1, tey_xyy_xxz_0, tey_xyy_xxz_1, tey_xyy_xyy_0, tey_xyy_xyy_1, \
                                         tey_xz_xxyy_0, tey_xz_xxyy_1, tey_xz_xxyz_0, tey_xz_xxyz_1, tey_xz_xxzz_0, \
                                         tey_xz_xxzz_1, tey_xz_xyyy_0, tey_xz_xyyy_1, tey_xz_xyyz_0, tey_xz_xyyz_1, \
                                         tey_xz_xyzz_0, tey_xz_xyzz_1, tey_xz_xzzz_0, tey_xz_xzzz_1, tey_xz_yyyy_0, \
                                         tey_xz_yyyy_1, tey_xz_yyyz_0, tey_xz_yyyz_1, tey_xz_yyzz_0, tey_xz_yyzz_1, \
                                         tey_xz_yzzz_0, tey_xz_yzzz_1, tey_xz_zzzz_0, tey_xz_zzzz_1, tey_yy_xxxx_0, \
                                         tey_yy_xxxx_1, tey_yy_xxxy_0, tey_yy_xxxy_1, tey_yy_xxxz_0, tey_yy_xxxz_1, \
                                         tey_yy_xxyy_0, tey_yy_xxyy_1, tez_xxxz_xxxz_0, tez_xxxz_xxyy_0, tez_xxxz_xxyz_0, \
                                         tez_xxxz_xxzz_0, tez_xxxz_xyyy_0, tez_xxxz_xyyz_0, tez_xxxz_xyzz_0, tez_xxxz_xzzz_0, \
                                         tez_xxxz_yyyy_0, tez_xxxz_yyyz_0, tez_xxxz_yyzz_0, tez_xxxz_yzzz_0, tez_xxxz_zzzz_0, \
                                         tez_xxyy_xxxx_0, tez_xxyy_xxxy_0, tez_xxyy_xxxz_0, tez_xxyy_xxyy_0, tez_xxz_xxxz_0, \
                                         tez_xxz_xxxz_1, tez_xxz_xxyy_0, tez_xxz_xxyy_1, tez_xxz_xxyz_0, tez_xxz_xxyz_1, \
                                         tez_xxz_xxz_0, tez_xxz_xxz_1, tez_xxz_xxzz_0, tez_xxz_xxzz_1, tez_xxz_xyy_0, \
                                         tez_xxz_xyy_1, tez_xxz_xyyy_0, tez_xxz_xyyy_1, tez_xxz_xyyz_0, tez_xxz_xyyz_1, \
                                         tez_xxz_xyz_0, tez_xxz_xyz_1, tez_xxz_xyzz_0, tez_xxz_xyzz_1, tez_xxz_xzz_0, \
                                         tez_xxz_xzz_1, tez_xxz_xzzz_0, tez_xxz_xzzz_1, tez_xxz_yyy_0, tez_xxz_yyy_1, \
                                         tez_xxz_yyyy_0, tez_xxz_yyyy_1, tez_xxz_yyyz_0, tez_xxz_yyyz_1, tez_xxz_yyz_0, \
                                         tez_xxz_yyz_1, tez_xxz_yyzz_0, tez_xxz_yyzz_1, tez_xxz_yzz_0, tez_xxz_yzz_1, \
                                         tez_xxz_yzzz_0, tez_xxz_yzzz_1, tez_xxz_zzz_0, tez_xxz_zzz_1, tez_xxz_zzzz_0, \
                                         tez_xxz_zzzz_1, tez_xyy_xxx_0, tez_xyy_xxx_1, tez_xyy_xxxx_0, tez_xyy_xxxx_1, \
                                         tez_xyy_xxxy_0, tez_xyy_xxxy_1, tez_xyy_xxxz_0, tez_xyy_xxxz_1, tez_xyy_xxy_0, \
                                         tez_xyy_xxy_1, tez_xyy_xxyy_0, tez_xyy_xxyy_1, tez_xyy_xxz_0, tez_xyy_xxz_1, \
                                         tez_xyy_xyy_0, tez_xyy_xyy_1, tez_xz_xxxz_0, tez_xz_xxxz_1, tez_xz_xxyy_0, \
                                         tez_xz_xxyy_1, tez_xz_xxyz_0, tez_xz_xxyz_1, tez_xz_xxzz_0, tez_xz_xxzz_1, \
                                         tez_xz_xyyy_0, tez_xz_xyyy_1, tez_xz_xyyz_0, tez_xz_xyyz_1, tez_xz_xyzz_0, \
                                         tez_xz_xyzz_1, tez_xz_xzzz_0, tez_xz_xzzz_1, tez_xz_yyyy_0, tez_xz_yyyy_1, \
                                         tez_xz_yyyz_0, tez_xz_yyyz_1, tez_xz_yyzz_0, tez_xz_yyzz_1, tez_xz_yzzz_0, \
                                         tez_xz_yzzz_1, tez_xz_zzzz_0, tez_xz_zzzz_1, tez_yy_xxxx_0, tez_yy_xxxx_1, \
                                         tez_yy_xxxy_0, tez_yy_xxxy_1, tez_yy_xxxz_0, tez_yy_xxxz_1, tez_yy_xxyy_0, \
                                         tez_yy_xxyy_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tez_xxxz_xxxz_0[j] = pa_x[j] * tez_xxz_xxxz_0[j] - pc_x[j] * tez_xxz_xxxz_1[j] + fl1_fx * tez_xz_xxxz_0[j] -
                                     fl1_fx * tez_xz_xxxz_1[j] + 1.5 * fl1_fx * tez_xxz_xxz_0[j] - 1.5 * fl1_fx * tez_xxz_xxz_1[j];

                tex_xxxz_xxyy_0[j] = pa_x[j] * tex_xxz_xxyy_0[j] - pc_x[j] * tex_xxz_xxyy_1[j] + fl1_fx * tex_xz_xxyy_0[j] -
                                     fl1_fx * tex_xz_xxyy_1[j] + fl1_fx * tex_xxz_xyy_0[j] - fl1_fx * tex_xxz_xyy_1[j] + ta_xxz_xxyy_1[j];

                tey_xxxz_xxyy_0[j] = pa_x[j] * tey_xxz_xxyy_0[j] - pc_x[j] * tey_xxz_xxyy_1[j] + fl1_fx * tey_xz_xxyy_0[j] -
                                     fl1_fx * tey_xz_xxyy_1[j] + fl1_fx * tey_xxz_xyy_0[j] - fl1_fx * tey_xxz_xyy_1[j];

                tez_xxxz_xxyy_0[j] = pa_x[j] * tez_xxz_xxyy_0[j] - pc_x[j] * tez_xxz_xxyy_1[j] + fl1_fx * tez_xz_xxyy_0[j] -
                                     fl1_fx * tez_xz_xxyy_1[j] + fl1_fx * tez_xxz_xyy_0[j] - fl1_fx * tez_xxz_xyy_1[j];

                tex_xxxz_xxyz_0[j] = pa_x[j] * tex_xxz_xxyz_0[j] - pc_x[j] * tex_xxz_xxyz_1[j] + fl1_fx * tex_xz_xxyz_0[j] -
                                     fl1_fx * tex_xz_xxyz_1[j] + fl1_fx * tex_xxz_xyz_0[j] - fl1_fx * tex_xxz_xyz_1[j] + ta_xxz_xxyz_1[j];

                tey_xxxz_xxyz_0[j] = pa_x[j] * tey_xxz_xxyz_0[j] - pc_x[j] * tey_xxz_xxyz_1[j] + fl1_fx * tey_xz_xxyz_0[j] -
                                     fl1_fx * tey_xz_xxyz_1[j] + fl1_fx * tey_xxz_xyz_0[j] - fl1_fx * tey_xxz_xyz_1[j];

                tez_xxxz_xxyz_0[j] = pa_x[j] * tez_xxz_xxyz_0[j] - pc_x[j] * tez_xxz_xxyz_1[j] + fl1_fx * tez_xz_xxyz_0[j] -
                                     fl1_fx * tez_xz_xxyz_1[j] + fl1_fx * tez_xxz_xyz_0[j] - fl1_fx * tez_xxz_xyz_1[j];

                tex_xxxz_xxzz_0[j] = pa_x[j] * tex_xxz_xxzz_0[j] - pc_x[j] * tex_xxz_xxzz_1[j] + fl1_fx * tex_xz_xxzz_0[j] -
                                     fl1_fx * tex_xz_xxzz_1[j] + fl1_fx * tex_xxz_xzz_0[j] - fl1_fx * tex_xxz_xzz_1[j] + ta_xxz_xxzz_1[j];

                tey_xxxz_xxzz_0[j] = pa_x[j] * tey_xxz_xxzz_0[j] - pc_x[j] * tey_xxz_xxzz_1[j] + fl1_fx * tey_xz_xxzz_0[j] -
                                     fl1_fx * tey_xz_xxzz_1[j] + fl1_fx * tey_xxz_xzz_0[j] - fl1_fx * tey_xxz_xzz_1[j];

                tez_xxxz_xxzz_0[j] = pa_x[j] * tez_xxz_xxzz_0[j] - pc_x[j] * tez_xxz_xxzz_1[j] + fl1_fx * tez_xz_xxzz_0[j] -
                                     fl1_fx * tez_xz_xxzz_1[j] + fl1_fx * tez_xxz_xzz_0[j] - fl1_fx * tez_xxz_xzz_1[j];

                tex_xxxz_xyyy_0[j] = pa_x[j] * tex_xxz_xyyy_0[j] - pc_x[j] * tex_xxz_xyyy_1[j] + fl1_fx * tex_xz_xyyy_0[j] -
                                     fl1_fx * tex_xz_xyyy_1[j] + 0.5 * fl1_fx * tex_xxz_yyy_0[j] - 0.5 * fl1_fx * tex_xxz_yyy_1[j] + ta_xxz_xyyy_1[j];

                tey_xxxz_xyyy_0[j] = pa_x[j] * tey_xxz_xyyy_0[j] - pc_x[j] * tey_xxz_xyyy_1[j] + fl1_fx * tey_xz_xyyy_0[j] -
                                     fl1_fx * tey_xz_xyyy_1[j] + 0.5 * fl1_fx * tey_xxz_yyy_0[j] - 0.5 * fl1_fx * tey_xxz_yyy_1[j];

                tez_xxxz_xyyy_0[j] = pa_x[j] * tez_xxz_xyyy_0[j] - pc_x[j] * tez_xxz_xyyy_1[j] + fl1_fx * tez_xz_xyyy_0[j] -
                                     fl1_fx * tez_xz_xyyy_1[j] + 0.5 * fl1_fx * tez_xxz_yyy_0[j] - 0.5 * fl1_fx * tez_xxz_yyy_1[j];

                tex_xxxz_xyyz_0[j] = pa_x[j] * tex_xxz_xyyz_0[j] - pc_x[j] * tex_xxz_xyyz_1[j] + fl1_fx * tex_xz_xyyz_0[j] -
                                     fl1_fx * tex_xz_xyyz_1[j] + 0.5 * fl1_fx * tex_xxz_yyz_0[j] - 0.5 * fl1_fx * tex_xxz_yyz_1[j] + ta_xxz_xyyz_1[j];

                tey_xxxz_xyyz_0[j] = pa_x[j] * tey_xxz_xyyz_0[j] - pc_x[j] * tey_xxz_xyyz_1[j] + fl1_fx * tey_xz_xyyz_0[j] -
                                     fl1_fx * tey_xz_xyyz_1[j] + 0.5 * fl1_fx * tey_xxz_yyz_0[j] - 0.5 * fl1_fx * tey_xxz_yyz_1[j];

                tez_xxxz_xyyz_0[j] = pa_x[j] * tez_xxz_xyyz_0[j] - pc_x[j] * tez_xxz_xyyz_1[j] + fl1_fx * tez_xz_xyyz_0[j] -
                                     fl1_fx * tez_xz_xyyz_1[j] + 0.5 * fl1_fx * tez_xxz_yyz_0[j] - 0.5 * fl1_fx * tez_xxz_yyz_1[j];

                tex_xxxz_xyzz_0[j] = pa_x[j] * tex_xxz_xyzz_0[j] - pc_x[j] * tex_xxz_xyzz_1[j] + fl1_fx * tex_xz_xyzz_0[j] -
                                     fl1_fx * tex_xz_xyzz_1[j] + 0.5 * fl1_fx * tex_xxz_yzz_0[j] - 0.5 * fl1_fx * tex_xxz_yzz_1[j] + ta_xxz_xyzz_1[j];

                tey_xxxz_xyzz_0[j] = pa_x[j] * tey_xxz_xyzz_0[j] - pc_x[j] * tey_xxz_xyzz_1[j] + fl1_fx * tey_xz_xyzz_0[j] -
                                     fl1_fx * tey_xz_xyzz_1[j] + 0.5 * fl1_fx * tey_xxz_yzz_0[j] - 0.5 * fl1_fx * tey_xxz_yzz_1[j];

                tez_xxxz_xyzz_0[j] = pa_x[j] * tez_xxz_xyzz_0[j] - pc_x[j] * tez_xxz_xyzz_1[j] + fl1_fx * tez_xz_xyzz_0[j] -
                                     fl1_fx * tez_xz_xyzz_1[j] + 0.5 * fl1_fx * tez_xxz_yzz_0[j] - 0.5 * fl1_fx * tez_xxz_yzz_1[j];

                tex_xxxz_xzzz_0[j] = pa_x[j] * tex_xxz_xzzz_0[j] - pc_x[j] * tex_xxz_xzzz_1[j] + fl1_fx * tex_xz_xzzz_0[j] -
                                     fl1_fx * tex_xz_xzzz_1[j] + 0.5 * fl1_fx * tex_xxz_zzz_0[j] - 0.5 * fl1_fx * tex_xxz_zzz_1[j] + ta_xxz_xzzz_1[j];

                tey_xxxz_xzzz_0[j] = pa_x[j] * tey_xxz_xzzz_0[j] - pc_x[j] * tey_xxz_xzzz_1[j] + fl1_fx * tey_xz_xzzz_0[j] -
                                     fl1_fx * tey_xz_xzzz_1[j] + 0.5 * fl1_fx * tey_xxz_zzz_0[j] - 0.5 * fl1_fx * tey_xxz_zzz_1[j];

                tez_xxxz_xzzz_0[j] = pa_x[j] * tez_xxz_xzzz_0[j] - pc_x[j] * tez_xxz_xzzz_1[j] + fl1_fx * tez_xz_xzzz_0[j] -
                                     fl1_fx * tez_xz_xzzz_1[j] + 0.5 * fl1_fx * tez_xxz_zzz_0[j] - 0.5 * fl1_fx * tez_xxz_zzz_1[j];

                tex_xxxz_yyyy_0[j] = pa_x[j] * tex_xxz_yyyy_0[j] - pc_x[j] * tex_xxz_yyyy_1[j] + fl1_fx * tex_xz_yyyy_0[j] -
                                     fl1_fx * tex_xz_yyyy_1[j] + ta_xxz_yyyy_1[j];

                tey_xxxz_yyyy_0[j] =
                    pa_x[j] * tey_xxz_yyyy_0[j] - pc_x[j] * tey_xxz_yyyy_1[j] + fl1_fx * tey_xz_yyyy_0[j] - fl1_fx * tey_xz_yyyy_1[j];

                tez_xxxz_yyyy_0[j] =
                    pa_x[j] * tez_xxz_yyyy_0[j] - pc_x[j] * tez_xxz_yyyy_1[j] + fl1_fx * tez_xz_yyyy_0[j] - fl1_fx * tez_xz_yyyy_1[j];

                tex_xxxz_yyyz_0[j] = pa_x[j] * tex_xxz_yyyz_0[j] - pc_x[j] * tex_xxz_yyyz_1[j] + fl1_fx * tex_xz_yyyz_0[j] -
                                     fl1_fx * tex_xz_yyyz_1[j] + ta_xxz_yyyz_1[j];

                tey_xxxz_yyyz_0[j] =
                    pa_x[j] * tey_xxz_yyyz_0[j] - pc_x[j] * tey_xxz_yyyz_1[j] + fl1_fx * tey_xz_yyyz_0[j] - fl1_fx * tey_xz_yyyz_1[j];

                tez_xxxz_yyyz_0[j] =
                    pa_x[j] * tez_xxz_yyyz_0[j] - pc_x[j] * tez_xxz_yyyz_1[j] + fl1_fx * tez_xz_yyyz_0[j] - fl1_fx * tez_xz_yyyz_1[j];

                tex_xxxz_yyzz_0[j] = pa_x[j] * tex_xxz_yyzz_0[j] - pc_x[j] * tex_xxz_yyzz_1[j] + fl1_fx * tex_xz_yyzz_0[j] -
                                     fl1_fx * tex_xz_yyzz_1[j] + ta_xxz_yyzz_1[j];

                tey_xxxz_yyzz_0[j] =
                    pa_x[j] * tey_xxz_yyzz_0[j] - pc_x[j] * tey_xxz_yyzz_1[j] + fl1_fx * tey_xz_yyzz_0[j] - fl1_fx * tey_xz_yyzz_1[j];

                tez_xxxz_yyzz_0[j] =
                    pa_x[j] * tez_xxz_yyzz_0[j] - pc_x[j] * tez_xxz_yyzz_1[j] + fl1_fx * tez_xz_yyzz_0[j] - fl1_fx * tez_xz_yyzz_1[j];

                tex_xxxz_yzzz_0[j] = pa_x[j] * tex_xxz_yzzz_0[j] - pc_x[j] * tex_xxz_yzzz_1[j] + fl1_fx * tex_xz_yzzz_0[j] -
                                     fl1_fx * tex_xz_yzzz_1[j] + ta_xxz_yzzz_1[j];

                tey_xxxz_yzzz_0[j] =
                    pa_x[j] * tey_xxz_yzzz_0[j] - pc_x[j] * tey_xxz_yzzz_1[j] + fl1_fx * tey_xz_yzzz_0[j] - fl1_fx * tey_xz_yzzz_1[j];

                tez_xxxz_yzzz_0[j] =
                    pa_x[j] * tez_xxz_yzzz_0[j] - pc_x[j] * tez_xxz_yzzz_1[j] + fl1_fx * tez_xz_yzzz_0[j] - fl1_fx * tez_xz_yzzz_1[j];

                tex_xxxz_zzzz_0[j] = pa_x[j] * tex_xxz_zzzz_0[j] - pc_x[j] * tex_xxz_zzzz_1[j] + fl1_fx * tex_xz_zzzz_0[j] -
                                     fl1_fx * tex_xz_zzzz_1[j] + ta_xxz_zzzz_1[j];

                tey_xxxz_zzzz_0[j] =
                    pa_x[j] * tey_xxz_zzzz_0[j] - pc_x[j] * tey_xxz_zzzz_1[j] + fl1_fx * tey_xz_zzzz_0[j] - fl1_fx * tey_xz_zzzz_1[j];

                tez_xxxz_zzzz_0[j] =
                    pa_x[j] * tez_xxz_zzzz_0[j] - pc_x[j] * tez_xxz_zzzz_1[j] + fl1_fx * tez_xz_zzzz_0[j] - fl1_fx * tez_xz_zzzz_1[j];

                tex_xxyy_xxxx_0[j] = pa_x[j] * tex_xyy_xxxx_0[j] - pc_x[j] * tex_xyy_xxxx_1[j] + 0.5 * fl1_fx * tex_yy_xxxx_0[j] -
                                     0.5 * fl1_fx * tex_yy_xxxx_1[j] + 2.0 * fl1_fx * tex_xyy_xxx_0[j] - 2.0 * fl1_fx * tex_xyy_xxx_1[j] +
                                     ta_xyy_xxxx_1[j];

                tey_xxyy_xxxx_0[j] = pa_x[j] * tey_xyy_xxxx_0[j] - pc_x[j] * tey_xyy_xxxx_1[j] + 0.5 * fl1_fx * tey_yy_xxxx_0[j] -
                                     0.5 * fl1_fx * tey_yy_xxxx_1[j] + 2.0 * fl1_fx * tey_xyy_xxx_0[j] - 2.0 * fl1_fx * tey_xyy_xxx_1[j];

                tez_xxyy_xxxx_0[j] = pa_x[j] * tez_xyy_xxxx_0[j] - pc_x[j] * tez_xyy_xxxx_1[j] + 0.5 * fl1_fx * tez_yy_xxxx_0[j] -
                                     0.5 * fl1_fx * tez_yy_xxxx_1[j] + 2.0 * fl1_fx * tez_xyy_xxx_0[j] - 2.0 * fl1_fx * tez_xyy_xxx_1[j];

                tex_xxyy_xxxy_0[j] = pa_x[j] * tex_xyy_xxxy_0[j] - pc_x[j] * tex_xyy_xxxy_1[j] + 0.5 * fl1_fx * tex_yy_xxxy_0[j] -
                                     0.5 * fl1_fx * tex_yy_xxxy_1[j] + 1.5 * fl1_fx * tex_xyy_xxy_0[j] - 1.5 * fl1_fx * tex_xyy_xxy_1[j] +
                                     ta_xyy_xxxy_1[j];

                tey_xxyy_xxxy_0[j] = pa_x[j] * tey_xyy_xxxy_0[j] - pc_x[j] * tey_xyy_xxxy_1[j] + 0.5 * fl1_fx * tey_yy_xxxy_0[j] -
                                     0.5 * fl1_fx * tey_yy_xxxy_1[j] + 1.5 * fl1_fx * tey_xyy_xxy_0[j] - 1.5 * fl1_fx * tey_xyy_xxy_1[j];

                tez_xxyy_xxxy_0[j] = pa_x[j] * tez_xyy_xxxy_0[j] - pc_x[j] * tez_xyy_xxxy_1[j] + 0.5 * fl1_fx * tez_yy_xxxy_0[j] -
                                     0.5 * fl1_fx * tez_yy_xxxy_1[j] + 1.5 * fl1_fx * tez_xyy_xxy_0[j] - 1.5 * fl1_fx * tez_xyy_xxy_1[j];

                tex_xxyy_xxxz_0[j] = pa_x[j] * tex_xyy_xxxz_0[j] - pc_x[j] * tex_xyy_xxxz_1[j] + 0.5 * fl1_fx * tex_yy_xxxz_0[j] -
                                     0.5 * fl1_fx * tex_yy_xxxz_1[j] + 1.5 * fl1_fx * tex_xyy_xxz_0[j] - 1.5 * fl1_fx * tex_xyy_xxz_1[j] +
                                     ta_xyy_xxxz_1[j];

                tey_xxyy_xxxz_0[j] = pa_x[j] * tey_xyy_xxxz_0[j] - pc_x[j] * tey_xyy_xxxz_1[j] + 0.5 * fl1_fx * tey_yy_xxxz_0[j] -
                                     0.5 * fl1_fx * tey_yy_xxxz_1[j] + 1.5 * fl1_fx * tey_xyy_xxz_0[j] - 1.5 * fl1_fx * tey_xyy_xxz_1[j];

                tez_xxyy_xxxz_0[j] = pa_x[j] * tez_xyy_xxxz_0[j] - pc_x[j] * tez_xyy_xxxz_1[j] + 0.5 * fl1_fx * tez_yy_xxxz_0[j] -
                                     0.5 * fl1_fx * tez_yy_xxxz_1[j] + 1.5 * fl1_fx * tez_xyy_xxz_0[j] - 1.5 * fl1_fx * tez_xyy_xxz_1[j];

                tex_xxyy_xxyy_0[j] = pa_x[j] * tex_xyy_xxyy_0[j] - pc_x[j] * tex_xyy_xxyy_1[j] + 0.5 * fl1_fx * tex_yy_xxyy_0[j] -
                                     0.5 * fl1_fx * tex_yy_xxyy_1[j] + fl1_fx * tex_xyy_xyy_0[j] - fl1_fx * tex_xyy_xyy_1[j] + ta_xyy_xxyy_1[j];

                tey_xxyy_xxyy_0[j] = pa_x[j] * tey_xyy_xxyy_0[j] - pc_x[j] * tey_xyy_xxyy_1[j] + 0.5 * fl1_fx * tey_yy_xxyy_0[j] -
                                     0.5 * fl1_fx * tey_yy_xxyy_1[j] + fl1_fx * tey_xyy_xyy_0[j] - fl1_fx * tey_xyy_xyy_1[j];

                tez_xxyy_xxyy_0[j] = pa_x[j] * tez_xyy_xxyy_0[j] - pc_x[j] * tez_xyy_xxyy_1[j] + 0.5 * fl1_fx * tez_yy_xxyy_0[j] -
                                     0.5 * fl1_fx * tez_yy_xxyy_1[j] + fl1_fx * tez_xyy_xyy_0[j] - fl1_fx * tez_xyy_xyy_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForGG_147_195(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pcDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto)
{
    // Batch of Integrals (147,195)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_4_4_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_2_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tex_xyy_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 49);

            auto tey_xyy_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 49);

            auto tez_xyy_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 49);

            auto tex_xyy_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 50);

            auto tey_xyy_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 50);

            auto tez_xyy_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 50);

            auto tex_xyy_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 51);

            auto tey_xyy_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 51);

            auto tez_xyy_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 51);

            auto tex_xyy_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 52);

            auto tey_xyy_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 52);

            auto tez_xyy_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 52);

            auto tex_xyy_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 53);

            auto tey_xyy_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 53);

            auto tez_xyy_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 53);

            auto tex_xyy_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 54);

            auto tey_xyy_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 54);

            auto tez_xyy_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 54);

            auto tex_xyy_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 55);

            auto tey_xyy_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 55);

            auto tez_xyy_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 55);

            auto tex_xyy_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 56);

            auto tey_xyy_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 56);

            auto tez_xyy_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 56);

            auto tex_xyy_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 57);

            auto tey_xyy_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 57);

            auto tez_xyy_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 57);

            auto tex_xyy_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 58);

            auto tey_xyy_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 58);

            auto tez_xyy_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 58);

            auto tex_xyy_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 59);

            auto tey_xyy_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 59);

            auto tez_xyy_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 59);

            auto tex_xyz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 60);

            auto tey_xyz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 60);

            auto tez_xyz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 60);

            auto tex_xyz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 61);

            auto tey_xyz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 61);

            auto tez_xyz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 61);

            auto tex_xyz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 62);

            auto tey_xyz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 62);

            auto tez_xyz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 62);

            auto tex_xyz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 63);

            auto tey_xyz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 63);

            auto tez_xyz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 63);

            auto tex_xyz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 64);

            auto tey_xyz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 64);

            auto tez_xyz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 64);

            auto tex_xyy_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 49);

            auto tey_xyy_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 49);

            auto tez_xyy_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 49);

            auto tex_xyy_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 50);

            auto tey_xyy_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 50);

            auto tez_xyy_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 50);

            auto tex_xyy_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 51);

            auto tey_xyy_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 51);

            auto tez_xyy_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 51);

            auto tex_xyy_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 52);

            auto tey_xyy_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 52);

            auto tez_xyy_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 52);

            auto tex_xyy_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 53);

            auto tey_xyy_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 53);

            auto tez_xyy_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 53);

            auto tex_xyy_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 54);

            auto tey_xyy_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 54);

            auto tez_xyy_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 54);

            auto tex_xyy_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 55);

            auto tey_xyy_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 55);

            auto tez_xyy_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 55);

            auto tex_xyy_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 56);

            auto tey_xyy_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 56);

            auto tez_xyy_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 56);

            auto tex_xyy_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 57);

            auto tey_xyy_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 57);

            auto tez_xyy_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 57);

            auto tex_xyy_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 58);

            auto tey_xyy_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 58);

            auto tez_xyy_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 58);

            auto tex_xyy_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 59);

            auto tey_xyy_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 59);

            auto tez_xyy_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 59);

            auto tex_xyz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 60);

            auto tey_xyz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 60);

            auto tez_xyz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 60);

            auto tex_xyz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 61);

            auto tey_xyz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 61);

            auto tez_xyz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 61);

            auto tex_xyz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 62);

            auto tey_xyz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 62);

            auto tez_xyz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 62);

            auto tex_xyz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 63);

            auto tey_xyz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 63);

            auto tez_xyz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 63);

            auto tex_xyz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 64);

            auto tey_xyz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 64);

            auto tez_xyz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 64);

            auto tex_yy_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 49);

            auto tey_yy_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 49);

            auto tez_yy_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 49);

            auto tex_yy_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 50);

            auto tey_yy_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 50);

            auto tez_yy_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 50);

            auto tex_yy_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 51);

            auto tey_yy_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 51);

            auto tez_yy_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 51);

            auto tex_yy_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 52);

            auto tey_yy_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 52);

            auto tez_yy_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 52);

            auto tex_yy_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 53);

            auto tey_yy_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 53);

            auto tez_yy_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 53);

            auto tex_yy_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 54);

            auto tey_yy_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 54);

            auto tez_yy_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 54);

            auto tex_yy_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 55);

            auto tey_yy_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 55);

            auto tez_yy_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 55);

            auto tex_yy_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 56);

            auto tey_yy_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 56);

            auto tez_yy_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 56);

            auto tex_yy_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 57);

            auto tey_yy_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 57);

            auto tez_yy_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 57);

            auto tex_yy_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 58);

            auto tey_yy_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 58);

            auto tez_yy_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 58);

            auto tex_yy_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 59);

            auto tey_yy_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 59);

            auto tez_yy_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 59);

            auto tex_yz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 60);

            auto tey_yz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 60);

            auto tez_yz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 60);

            auto tex_yz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 61);

            auto tey_yz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 61);

            auto tez_yz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 61);

            auto tex_yz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 62);

            auto tey_yz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 62);

            auto tez_yz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 62);

            auto tex_yz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 63);

            auto tey_yz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 63);

            auto tez_yz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 63);

            auto tex_yz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 64);

            auto tey_yz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 64);

            auto tez_yz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 64);

            auto tex_yy_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 49);

            auto tey_yy_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 49);

            auto tez_yy_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 49);

            auto tex_yy_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 50);

            auto tey_yy_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 50);

            auto tez_yy_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 50);

            auto tex_yy_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 51);

            auto tey_yy_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 51);

            auto tez_yy_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 51);

            auto tex_yy_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 52);

            auto tey_yy_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 52);

            auto tez_yy_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 52);

            auto tex_yy_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 53);

            auto tey_yy_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 53);

            auto tez_yy_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 53);

            auto tex_yy_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 54);

            auto tey_yy_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 54);

            auto tez_yy_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 54);

            auto tex_yy_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 55);

            auto tey_yy_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 55);

            auto tez_yy_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 55);

            auto tex_yy_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 56);

            auto tey_yy_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 56);

            auto tez_yy_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 56);

            auto tex_yy_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 57);

            auto tey_yy_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 57);

            auto tez_yy_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 57);

            auto tex_yy_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 58);

            auto tey_yy_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 58);

            auto tez_yy_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 58);

            auto tex_yy_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 59);

            auto tey_yy_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 59);

            auto tez_yy_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 59);

            auto tex_yz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 60);

            auto tey_yz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 60);

            auto tez_yz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 60);

            auto tex_yz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 61);

            auto tey_yz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 61);

            auto tez_yz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 61);

            auto tex_yz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 62);

            auto tey_yz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 62);

            auto tez_yz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 62);

            auto tex_yz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 63);

            auto tey_yz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 63);

            auto tez_yz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 63);

            auto tex_yz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 64);

            auto tey_yz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 64);

            auto tez_yz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 64);

            auto tex_xyy_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 34);

            auto tey_xyy_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 34);

            auto tez_xyy_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 34);

            auto tex_xyy_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 35);

            auto tey_xyy_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 35);

            auto tez_xyy_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 35);

            auto tex_xyy_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 36);

            auto tey_xyy_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 36);

            auto tez_xyy_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 36);

            auto tex_xyy_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 37);

            auto tey_xyy_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 37);

            auto tez_xyy_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 37);

            auto tex_xyy_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 38);

            auto tey_xyy_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 38);

            auto tez_xyy_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 38);

            auto tex_xyy_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 39);

            auto tey_xyy_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 39);

            auto tez_xyy_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 39);

            auto tex_xyz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 40);

            auto tey_xyz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 40);

            auto tez_xyz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 40);

            auto tex_xyz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 41);

            auto tey_xyz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 41);

            auto tez_xyz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 41);

            auto tex_xyz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 42);

            auto tey_xyz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 42);

            auto tez_xyz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 42);

            auto tex_xyz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 43);

            auto tey_xyz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 43);

            auto tez_xyz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 43);

            auto tex_xyz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 44);

            auto tey_xyz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 44);

            auto tez_xyz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 44);

            auto tex_xyy_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 34);

            auto tey_xyy_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 34);

            auto tez_xyy_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 34);

            auto tex_xyy_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 35);

            auto tey_xyy_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 35);

            auto tez_xyy_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 35);

            auto tex_xyy_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 36);

            auto tey_xyy_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 36);

            auto tez_xyy_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 36);

            auto tex_xyy_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 37);

            auto tey_xyy_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 37);

            auto tez_xyy_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 37);

            auto tex_xyy_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 38);

            auto tey_xyy_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 38);

            auto tez_xyy_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 38);

            auto tex_xyy_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 39);

            auto tey_xyy_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 39);

            auto tez_xyy_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 39);

            auto tex_xyz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 40);

            auto tey_xyz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 40);

            auto tez_xyz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 40);

            auto tex_xyz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 41);

            auto tey_xyz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 41);

            auto tez_xyz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 41);

            auto tex_xyz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 42);

            auto tey_xyz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 42);

            auto tez_xyz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 42);

            auto tex_xyz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 43);

            auto tey_xyz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 43);

            auto tez_xyz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 43);

            auto tex_xyz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 44);

            auto tey_xyz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 44);

            auto tez_xyz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 44);

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

            // set up pointers to integrals

            auto tex_xxyy_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 49);

            auto tey_xxyy_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 49);

            auto tez_xxyy_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 49);

            auto tex_xxyy_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 50);

            auto tey_xxyy_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 50);

            auto tez_xxyy_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 50);

            auto tex_xxyy_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 51);

            auto tey_xxyy_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 51);

            auto tez_xxyy_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 51);

            auto tex_xxyy_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 52);

            auto tey_xxyy_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 52);

            auto tez_xxyy_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 52);

            auto tex_xxyy_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 53);

            auto tey_xxyy_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 53);

            auto tez_xxyy_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 53);

            auto tex_xxyy_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 54);

            auto tey_xxyy_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 54);

            auto tez_xxyy_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 54);

            auto tex_xxyy_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 55);

            auto tey_xxyy_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 55);

            auto tez_xxyy_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 55);

            auto tex_xxyy_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 56);

            auto tey_xxyy_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 56);

            auto tez_xxyy_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 56);

            auto tex_xxyy_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 57);

            auto tey_xxyy_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 57);

            auto tez_xxyy_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 57);

            auto tex_xxyy_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 58);

            auto tey_xxyy_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 58);

            auto tez_xxyy_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 58);

            auto tex_xxyy_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 59);

            auto tey_xxyy_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 59);

            auto tez_xxyy_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 59);

            auto tex_xxyz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 60);

            auto tey_xxyz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 60);

            auto tez_xxyz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 60);

            auto tex_xxyz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 61);

            auto tey_xxyz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 61);

            auto tez_xxyz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 61);

            auto tex_xxyz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 62);

            auto tey_xxyz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 62);

            auto tez_xxyz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 62);

            auto tex_xxyz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 63);

            auto tey_xxyz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 63);

            auto tez_xxyz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 63);

            auto tex_xxyz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 64);

            auto tey_xxyz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 64);

            auto tez_xxyz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 64);

            // Batch of Integrals (147,195)

            #pragma omp simd aligned(fx, pa_x, pc_x, ta_xyy_xxyz_1, ta_xyy_xxzz_1, ta_xyy_xyyy_1, \
                                         ta_xyy_xyyz_1, ta_xyy_xyzz_1, ta_xyy_xzzz_1, ta_xyy_yyyy_1, ta_xyy_yyyz_1, \
                                         ta_xyy_yyzz_1, ta_xyy_yzzz_1, ta_xyy_zzzz_1, ta_xyz_xxxx_1, ta_xyz_xxxy_1, \
                                         ta_xyz_xxxz_1, ta_xyz_xxyy_1, ta_xyz_xxyz_1, tex_xxyy_xxyz_0, tex_xxyy_xxzz_0, \
                                         tex_xxyy_xyyy_0, tex_xxyy_xyyz_0, tex_xxyy_xyzz_0, tex_xxyy_xzzz_0, tex_xxyy_yyyy_0, \
                                         tex_xxyy_yyyz_0, tex_xxyy_yyzz_0, tex_xxyy_yzzz_0, tex_xxyy_zzzz_0, tex_xxyz_xxxx_0, \
                                         tex_xxyz_xxxy_0, tex_xxyz_xxxz_0, tex_xxyz_xxyy_0, tex_xxyz_xxyz_0, tex_xyy_xxyz_0, \
                                         tex_xyy_xxyz_1, tex_xyy_xxzz_0, tex_xyy_xxzz_1, tex_xyy_xyyy_0, tex_xyy_xyyy_1, \
                                         tex_xyy_xyyz_0, tex_xyy_xyyz_1, tex_xyy_xyz_0, tex_xyy_xyz_1, tex_xyy_xyzz_0, \
                                         tex_xyy_xyzz_1, tex_xyy_xzz_0, tex_xyy_xzz_1, tex_xyy_xzzz_0, tex_xyy_xzzz_1, \
                                         tex_xyy_yyy_0, tex_xyy_yyy_1, tex_xyy_yyyy_0, tex_xyy_yyyy_1, tex_xyy_yyyz_0, \
                                         tex_xyy_yyyz_1, tex_xyy_yyz_0, tex_xyy_yyz_1, tex_xyy_yyzz_0, tex_xyy_yyzz_1, \
                                         tex_xyy_yzz_0, tex_xyy_yzz_1, tex_xyy_yzzz_0, tex_xyy_yzzz_1, tex_xyy_zzz_0, \
                                         tex_xyy_zzz_1, tex_xyy_zzzz_0, tex_xyy_zzzz_1, tex_xyz_xxx_0, tex_xyz_xxx_1, \
                                         tex_xyz_xxxx_0, tex_xyz_xxxx_1, tex_xyz_xxxy_0, tex_xyz_xxxy_1, tex_xyz_xxxz_0, \
                                         tex_xyz_xxxz_1, tex_xyz_xxy_0, tex_xyz_xxy_1, tex_xyz_xxyy_0, tex_xyz_xxyy_1, \
                                         tex_xyz_xxyz_0, tex_xyz_xxyz_1, tex_xyz_xxz_0, tex_xyz_xxz_1, tex_xyz_xyy_0, \
                                         tex_xyz_xyy_1, tex_xyz_xyz_0, tex_xyz_xyz_1, tex_yy_xxyz_0, tex_yy_xxyz_1, \
                                         tex_yy_xxzz_0, tex_yy_xxzz_1, tex_yy_xyyy_0, tex_yy_xyyy_1, tex_yy_xyyz_0, \
                                         tex_yy_xyyz_1, tex_yy_xyzz_0, tex_yy_xyzz_1, tex_yy_xzzz_0, tex_yy_xzzz_1, \
                                         tex_yy_yyyy_0, tex_yy_yyyy_1, tex_yy_yyyz_0, tex_yy_yyyz_1, tex_yy_yyzz_0, \
                                         tex_yy_yyzz_1, tex_yy_yzzz_0, tex_yy_yzzz_1, tex_yy_zzzz_0, tex_yy_zzzz_1, \
                                         tex_yz_xxxx_0, tex_yz_xxxx_1, tex_yz_xxxy_0, tex_yz_xxxy_1, tex_yz_xxxz_0, \
                                         tex_yz_xxxz_1, tex_yz_xxyy_0, tex_yz_xxyy_1, tex_yz_xxyz_0, tex_yz_xxyz_1, \
                                         tey_xxyy_xxyz_0, tey_xxyy_xxzz_0, tey_xxyy_xyyy_0, tey_xxyy_xyyz_0, tey_xxyy_xyzz_0, \
                                         tey_xxyy_xzzz_0, tey_xxyy_yyyy_0, tey_xxyy_yyyz_0, tey_xxyy_yyzz_0, tey_xxyy_yzzz_0, \
                                         tey_xxyy_zzzz_0, tey_xxyz_xxxx_0, tey_xxyz_xxxy_0, tey_xxyz_xxxz_0, tey_xxyz_xxyy_0, \
                                         tey_xxyz_xxyz_0, tey_xyy_xxyz_0, tey_xyy_xxyz_1, tey_xyy_xxzz_0, tey_xyy_xxzz_1, \
                                         tey_xyy_xyyy_0, tey_xyy_xyyy_1, tey_xyy_xyyz_0, tey_xyy_xyyz_1, tey_xyy_xyz_0, \
                                         tey_xyy_xyz_1, tey_xyy_xyzz_0, tey_xyy_xyzz_1, tey_xyy_xzz_0, tey_xyy_xzz_1, \
                                         tey_xyy_xzzz_0, tey_xyy_xzzz_1, tey_xyy_yyy_0, tey_xyy_yyy_1, tey_xyy_yyyy_0, \
                                         tey_xyy_yyyy_1, tey_xyy_yyyz_0, tey_xyy_yyyz_1, tey_xyy_yyz_0, tey_xyy_yyz_1, \
                                         tey_xyy_yyzz_0, tey_xyy_yyzz_1, tey_xyy_yzz_0, tey_xyy_yzz_1, tey_xyy_yzzz_0, \
                                         tey_xyy_yzzz_1, tey_xyy_zzz_0, tey_xyy_zzz_1, tey_xyy_zzzz_0, tey_xyy_zzzz_1, \
                                         tey_xyz_xxx_0, tey_xyz_xxx_1, tey_xyz_xxxx_0, tey_xyz_xxxx_1, tey_xyz_xxxy_0, \
                                         tey_xyz_xxxy_1, tey_xyz_xxxz_0, tey_xyz_xxxz_1, tey_xyz_xxy_0, tey_xyz_xxy_1, \
                                         tey_xyz_xxyy_0, tey_xyz_xxyy_1, tey_xyz_xxyz_0, tey_xyz_xxyz_1, tey_xyz_xxz_0, \
                                         tey_xyz_xxz_1, tey_xyz_xyy_0, tey_xyz_xyy_1, tey_xyz_xyz_0, tey_xyz_xyz_1, \
                                         tey_yy_xxyz_0, tey_yy_xxyz_1, tey_yy_xxzz_0, tey_yy_xxzz_1, tey_yy_xyyy_0, \
                                         tey_yy_xyyy_1, tey_yy_xyyz_0, tey_yy_xyyz_1, tey_yy_xyzz_0, tey_yy_xyzz_1, \
                                         tey_yy_xzzz_0, tey_yy_xzzz_1, tey_yy_yyyy_0, tey_yy_yyyy_1, tey_yy_yyyz_0, \
                                         tey_yy_yyyz_1, tey_yy_yyzz_0, tey_yy_yyzz_1, tey_yy_yzzz_0, tey_yy_yzzz_1, \
                                         tey_yy_zzzz_0, tey_yy_zzzz_1, tey_yz_xxxx_0, tey_yz_xxxx_1, tey_yz_xxxy_0, \
                                         tey_yz_xxxy_1, tey_yz_xxxz_0, tey_yz_xxxz_1, tey_yz_xxyy_0, tey_yz_xxyy_1, \
                                         tey_yz_xxyz_0, tey_yz_xxyz_1, tez_xxyy_xxyz_0, tez_xxyy_xxzz_0, tez_xxyy_xyyy_0, \
                                         tez_xxyy_xyyz_0, tez_xxyy_xyzz_0, tez_xxyy_xzzz_0, tez_xxyy_yyyy_0, tez_xxyy_yyyz_0, \
                                         tez_xxyy_yyzz_0, tez_xxyy_yzzz_0, tez_xxyy_zzzz_0, tez_xxyz_xxxx_0, tez_xxyz_xxxy_0, \
                                         tez_xxyz_xxxz_0, tez_xxyz_xxyy_0, tez_xxyz_xxyz_0, tez_xyy_xxyz_0, tez_xyy_xxyz_1, \
                                         tez_xyy_xxzz_0, tez_xyy_xxzz_1, tez_xyy_xyyy_0, tez_xyy_xyyy_1, tez_xyy_xyyz_0, \
                                         tez_xyy_xyyz_1, tez_xyy_xyz_0, tez_xyy_xyz_1, tez_xyy_xyzz_0, tez_xyy_xyzz_1, \
                                         tez_xyy_xzz_0, tez_xyy_xzz_1, tez_xyy_xzzz_0, tez_xyy_xzzz_1, tez_xyy_yyy_0, \
                                         tez_xyy_yyy_1, tez_xyy_yyyy_0, tez_xyy_yyyy_1, tez_xyy_yyyz_0, tez_xyy_yyyz_1, \
                                         tez_xyy_yyz_0, tez_xyy_yyz_1, tez_xyy_yyzz_0, tez_xyy_yyzz_1, tez_xyy_yzz_0, \
                                         tez_xyy_yzz_1, tez_xyy_yzzz_0, tez_xyy_yzzz_1, tez_xyy_zzz_0, tez_xyy_zzz_1, \
                                         tez_xyy_zzzz_0, tez_xyy_zzzz_1, tez_xyz_xxx_0, tez_xyz_xxx_1, tez_xyz_xxxx_0, \
                                         tez_xyz_xxxx_1, tez_xyz_xxxy_0, tez_xyz_xxxy_1, tez_xyz_xxxz_0, tez_xyz_xxxz_1, \
                                         tez_xyz_xxy_0, tez_xyz_xxy_1, tez_xyz_xxyy_0, tez_xyz_xxyy_1, tez_xyz_xxyz_0, \
                                         tez_xyz_xxyz_1, tez_xyz_xxz_0, tez_xyz_xxz_1, tez_xyz_xyy_0, tez_xyz_xyy_1, \
                                         tez_xyz_xyz_0, tez_xyz_xyz_1, tez_yy_xxyz_0, tez_yy_xxyz_1, tez_yy_xxzz_0, \
                                         tez_yy_xxzz_1, tez_yy_xyyy_0, tez_yy_xyyy_1, tez_yy_xyyz_0, tez_yy_xyyz_1, \
                                         tez_yy_xyzz_0, tez_yy_xyzz_1, tez_yy_xzzz_0, tez_yy_xzzz_1, tez_yy_yyyy_0, \
                                         tez_yy_yyyy_1, tez_yy_yyyz_0, tez_yy_yyyz_1, tez_yy_yyzz_0, tez_yy_yyzz_1, \
                                         tez_yy_yzzz_0, tez_yy_yzzz_1, tez_yy_zzzz_0, tez_yy_zzzz_1, tez_yz_xxxx_0, \
                                         tez_yz_xxxx_1, tez_yz_xxxy_0, tez_yz_xxxy_1, tez_yz_xxxz_0, tez_yz_xxxz_1, \
                                         tez_yz_xxyy_0, tez_yz_xxyy_1, tez_yz_xxyz_0, tez_yz_xxyz_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_xxyy_xxyz_0[j] = pa_x[j] * tex_xyy_xxyz_0[j] - pc_x[j] * tex_xyy_xxyz_1[j] + 0.5 * fl1_fx * tex_yy_xxyz_0[j] -
                                     0.5 * fl1_fx * tex_yy_xxyz_1[j] + fl1_fx * tex_xyy_xyz_0[j] - fl1_fx * tex_xyy_xyz_1[j] + ta_xyy_xxyz_1[j];

                tey_xxyy_xxyz_0[j] = pa_x[j] * tey_xyy_xxyz_0[j] - pc_x[j] * tey_xyy_xxyz_1[j] + 0.5 * fl1_fx * tey_yy_xxyz_0[j] -
                                     0.5 * fl1_fx * tey_yy_xxyz_1[j] + fl1_fx * tey_xyy_xyz_0[j] - fl1_fx * tey_xyy_xyz_1[j];

                tez_xxyy_xxyz_0[j] = pa_x[j] * tez_xyy_xxyz_0[j] - pc_x[j] * tez_xyy_xxyz_1[j] + 0.5 * fl1_fx * tez_yy_xxyz_0[j] -
                                     0.5 * fl1_fx * tez_yy_xxyz_1[j] + fl1_fx * tez_xyy_xyz_0[j] - fl1_fx * tez_xyy_xyz_1[j];

                tex_xxyy_xxzz_0[j] = pa_x[j] * tex_xyy_xxzz_0[j] - pc_x[j] * tex_xyy_xxzz_1[j] + 0.5 * fl1_fx * tex_yy_xxzz_0[j] -
                                     0.5 * fl1_fx * tex_yy_xxzz_1[j] + fl1_fx * tex_xyy_xzz_0[j] - fl1_fx * tex_xyy_xzz_1[j] + ta_xyy_xxzz_1[j];

                tey_xxyy_xxzz_0[j] = pa_x[j] * tey_xyy_xxzz_0[j] - pc_x[j] * tey_xyy_xxzz_1[j] + 0.5 * fl1_fx * tey_yy_xxzz_0[j] -
                                     0.5 * fl1_fx * tey_yy_xxzz_1[j] + fl1_fx * tey_xyy_xzz_0[j] - fl1_fx * tey_xyy_xzz_1[j];

                tez_xxyy_xxzz_0[j] = pa_x[j] * tez_xyy_xxzz_0[j] - pc_x[j] * tez_xyy_xxzz_1[j] + 0.5 * fl1_fx * tez_yy_xxzz_0[j] -
                                     0.5 * fl1_fx * tez_yy_xxzz_1[j] + fl1_fx * tez_xyy_xzz_0[j] - fl1_fx * tez_xyy_xzz_1[j];

                tex_xxyy_xyyy_0[j] = pa_x[j] * tex_xyy_xyyy_0[j] - pc_x[j] * tex_xyy_xyyy_1[j] + 0.5 * fl1_fx * tex_yy_xyyy_0[j] -
                                     0.5 * fl1_fx * tex_yy_xyyy_1[j] + 0.5 * fl1_fx * tex_xyy_yyy_0[j] - 0.5 * fl1_fx * tex_xyy_yyy_1[j] +
                                     ta_xyy_xyyy_1[j];

                tey_xxyy_xyyy_0[j] = pa_x[j] * tey_xyy_xyyy_0[j] - pc_x[j] * tey_xyy_xyyy_1[j] + 0.5 * fl1_fx * tey_yy_xyyy_0[j] -
                                     0.5 * fl1_fx * tey_yy_xyyy_1[j] + 0.5 * fl1_fx * tey_xyy_yyy_0[j] - 0.5 * fl1_fx * tey_xyy_yyy_1[j];

                tez_xxyy_xyyy_0[j] = pa_x[j] * tez_xyy_xyyy_0[j] - pc_x[j] * tez_xyy_xyyy_1[j] + 0.5 * fl1_fx * tez_yy_xyyy_0[j] -
                                     0.5 * fl1_fx * tez_yy_xyyy_1[j] + 0.5 * fl1_fx * tez_xyy_yyy_0[j] - 0.5 * fl1_fx * tez_xyy_yyy_1[j];

                tex_xxyy_xyyz_0[j] = pa_x[j] * tex_xyy_xyyz_0[j] - pc_x[j] * tex_xyy_xyyz_1[j] + 0.5 * fl1_fx * tex_yy_xyyz_0[j] -
                                     0.5 * fl1_fx * tex_yy_xyyz_1[j] + 0.5 * fl1_fx * tex_xyy_yyz_0[j] - 0.5 * fl1_fx * tex_xyy_yyz_1[j] +
                                     ta_xyy_xyyz_1[j];

                tey_xxyy_xyyz_0[j] = pa_x[j] * tey_xyy_xyyz_0[j] - pc_x[j] * tey_xyy_xyyz_1[j] + 0.5 * fl1_fx * tey_yy_xyyz_0[j] -
                                     0.5 * fl1_fx * tey_yy_xyyz_1[j] + 0.5 * fl1_fx * tey_xyy_yyz_0[j] - 0.5 * fl1_fx * tey_xyy_yyz_1[j];

                tez_xxyy_xyyz_0[j] = pa_x[j] * tez_xyy_xyyz_0[j] - pc_x[j] * tez_xyy_xyyz_1[j] + 0.5 * fl1_fx * tez_yy_xyyz_0[j] -
                                     0.5 * fl1_fx * tez_yy_xyyz_1[j] + 0.5 * fl1_fx * tez_xyy_yyz_0[j] - 0.5 * fl1_fx * tez_xyy_yyz_1[j];

                tex_xxyy_xyzz_0[j] = pa_x[j] * tex_xyy_xyzz_0[j] - pc_x[j] * tex_xyy_xyzz_1[j] + 0.5 * fl1_fx * tex_yy_xyzz_0[j] -
                                     0.5 * fl1_fx * tex_yy_xyzz_1[j] + 0.5 * fl1_fx * tex_xyy_yzz_0[j] - 0.5 * fl1_fx * tex_xyy_yzz_1[j] +
                                     ta_xyy_xyzz_1[j];

                tey_xxyy_xyzz_0[j] = pa_x[j] * tey_xyy_xyzz_0[j] - pc_x[j] * tey_xyy_xyzz_1[j] + 0.5 * fl1_fx * tey_yy_xyzz_0[j] -
                                     0.5 * fl1_fx * tey_yy_xyzz_1[j] + 0.5 * fl1_fx * tey_xyy_yzz_0[j] - 0.5 * fl1_fx * tey_xyy_yzz_1[j];

                tez_xxyy_xyzz_0[j] = pa_x[j] * tez_xyy_xyzz_0[j] - pc_x[j] * tez_xyy_xyzz_1[j] + 0.5 * fl1_fx * tez_yy_xyzz_0[j] -
                                     0.5 * fl1_fx * tez_yy_xyzz_1[j] + 0.5 * fl1_fx * tez_xyy_yzz_0[j] - 0.5 * fl1_fx * tez_xyy_yzz_1[j];

                tex_xxyy_xzzz_0[j] = pa_x[j] * tex_xyy_xzzz_0[j] - pc_x[j] * tex_xyy_xzzz_1[j] + 0.5 * fl1_fx * tex_yy_xzzz_0[j] -
                                     0.5 * fl1_fx * tex_yy_xzzz_1[j] + 0.5 * fl1_fx * tex_xyy_zzz_0[j] - 0.5 * fl1_fx * tex_xyy_zzz_1[j] +
                                     ta_xyy_xzzz_1[j];

                tey_xxyy_xzzz_0[j] = pa_x[j] * tey_xyy_xzzz_0[j] - pc_x[j] * tey_xyy_xzzz_1[j] + 0.5 * fl1_fx * tey_yy_xzzz_0[j] -
                                     0.5 * fl1_fx * tey_yy_xzzz_1[j] + 0.5 * fl1_fx * tey_xyy_zzz_0[j] - 0.5 * fl1_fx * tey_xyy_zzz_1[j];

                tez_xxyy_xzzz_0[j] = pa_x[j] * tez_xyy_xzzz_0[j] - pc_x[j] * tez_xyy_xzzz_1[j] + 0.5 * fl1_fx * tez_yy_xzzz_0[j] -
                                     0.5 * fl1_fx * tez_yy_xzzz_1[j] + 0.5 * fl1_fx * tez_xyy_zzz_0[j] - 0.5 * fl1_fx * tez_xyy_zzz_1[j];

                tex_xxyy_yyyy_0[j] = pa_x[j] * tex_xyy_yyyy_0[j] - pc_x[j] * tex_xyy_yyyy_1[j] + 0.5 * fl1_fx * tex_yy_yyyy_0[j] -
                                     0.5 * fl1_fx * tex_yy_yyyy_1[j] + ta_xyy_yyyy_1[j];

                tey_xxyy_yyyy_0[j] =
                    pa_x[j] * tey_xyy_yyyy_0[j] - pc_x[j] * tey_xyy_yyyy_1[j] + 0.5 * fl1_fx * tey_yy_yyyy_0[j] - 0.5 * fl1_fx * tey_yy_yyyy_1[j];

                tez_xxyy_yyyy_0[j] =
                    pa_x[j] * tez_xyy_yyyy_0[j] - pc_x[j] * tez_xyy_yyyy_1[j] + 0.5 * fl1_fx * tez_yy_yyyy_0[j] - 0.5 * fl1_fx * tez_yy_yyyy_1[j];

                tex_xxyy_yyyz_0[j] = pa_x[j] * tex_xyy_yyyz_0[j] - pc_x[j] * tex_xyy_yyyz_1[j] + 0.5 * fl1_fx * tex_yy_yyyz_0[j] -
                                     0.5 * fl1_fx * tex_yy_yyyz_1[j] + ta_xyy_yyyz_1[j];

                tey_xxyy_yyyz_0[j] =
                    pa_x[j] * tey_xyy_yyyz_0[j] - pc_x[j] * tey_xyy_yyyz_1[j] + 0.5 * fl1_fx * tey_yy_yyyz_0[j] - 0.5 * fl1_fx * tey_yy_yyyz_1[j];

                tez_xxyy_yyyz_0[j] =
                    pa_x[j] * tez_xyy_yyyz_0[j] - pc_x[j] * tez_xyy_yyyz_1[j] + 0.5 * fl1_fx * tez_yy_yyyz_0[j] - 0.5 * fl1_fx * tez_yy_yyyz_1[j];

                tex_xxyy_yyzz_0[j] = pa_x[j] * tex_xyy_yyzz_0[j] - pc_x[j] * tex_xyy_yyzz_1[j] + 0.5 * fl1_fx * tex_yy_yyzz_0[j] -
                                     0.5 * fl1_fx * tex_yy_yyzz_1[j] + ta_xyy_yyzz_1[j];

                tey_xxyy_yyzz_0[j] =
                    pa_x[j] * tey_xyy_yyzz_0[j] - pc_x[j] * tey_xyy_yyzz_1[j] + 0.5 * fl1_fx * tey_yy_yyzz_0[j] - 0.5 * fl1_fx * tey_yy_yyzz_1[j];

                tez_xxyy_yyzz_0[j] =
                    pa_x[j] * tez_xyy_yyzz_0[j] - pc_x[j] * tez_xyy_yyzz_1[j] + 0.5 * fl1_fx * tez_yy_yyzz_0[j] - 0.5 * fl1_fx * tez_yy_yyzz_1[j];

                tex_xxyy_yzzz_0[j] = pa_x[j] * tex_xyy_yzzz_0[j] - pc_x[j] * tex_xyy_yzzz_1[j] + 0.5 * fl1_fx * tex_yy_yzzz_0[j] -
                                     0.5 * fl1_fx * tex_yy_yzzz_1[j] + ta_xyy_yzzz_1[j];

                tey_xxyy_yzzz_0[j] =
                    pa_x[j] * tey_xyy_yzzz_0[j] - pc_x[j] * tey_xyy_yzzz_1[j] + 0.5 * fl1_fx * tey_yy_yzzz_0[j] - 0.5 * fl1_fx * tey_yy_yzzz_1[j];

                tez_xxyy_yzzz_0[j] =
                    pa_x[j] * tez_xyy_yzzz_0[j] - pc_x[j] * tez_xyy_yzzz_1[j] + 0.5 * fl1_fx * tez_yy_yzzz_0[j] - 0.5 * fl1_fx * tez_yy_yzzz_1[j];

                tex_xxyy_zzzz_0[j] = pa_x[j] * tex_xyy_zzzz_0[j] - pc_x[j] * tex_xyy_zzzz_1[j] + 0.5 * fl1_fx * tex_yy_zzzz_0[j] -
                                     0.5 * fl1_fx * tex_yy_zzzz_1[j] + ta_xyy_zzzz_1[j];

                tey_xxyy_zzzz_0[j] =
                    pa_x[j] * tey_xyy_zzzz_0[j] - pc_x[j] * tey_xyy_zzzz_1[j] + 0.5 * fl1_fx * tey_yy_zzzz_0[j] - 0.5 * fl1_fx * tey_yy_zzzz_1[j];

                tez_xxyy_zzzz_0[j] =
                    pa_x[j] * tez_xyy_zzzz_0[j] - pc_x[j] * tez_xyy_zzzz_1[j] + 0.5 * fl1_fx * tez_yy_zzzz_0[j] - 0.5 * fl1_fx * tez_yy_zzzz_1[j];

                tex_xxyz_xxxx_0[j] = pa_x[j] * tex_xyz_xxxx_0[j] - pc_x[j] * tex_xyz_xxxx_1[j] + 0.5 * fl1_fx * tex_yz_xxxx_0[j] -
                                     0.5 * fl1_fx * tex_yz_xxxx_1[j] + 2.0 * fl1_fx * tex_xyz_xxx_0[j] - 2.0 * fl1_fx * tex_xyz_xxx_1[j] +
                                     ta_xyz_xxxx_1[j];

                tey_xxyz_xxxx_0[j] = pa_x[j] * tey_xyz_xxxx_0[j] - pc_x[j] * tey_xyz_xxxx_1[j] + 0.5 * fl1_fx * tey_yz_xxxx_0[j] -
                                     0.5 * fl1_fx * tey_yz_xxxx_1[j] + 2.0 * fl1_fx * tey_xyz_xxx_0[j] - 2.0 * fl1_fx * tey_xyz_xxx_1[j];

                tez_xxyz_xxxx_0[j] = pa_x[j] * tez_xyz_xxxx_0[j] - pc_x[j] * tez_xyz_xxxx_1[j] + 0.5 * fl1_fx * tez_yz_xxxx_0[j] -
                                     0.5 * fl1_fx * tez_yz_xxxx_1[j] + 2.0 * fl1_fx * tez_xyz_xxx_0[j] - 2.0 * fl1_fx * tez_xyz_xxx_1[j];

                tex_xxyz_xxxy_0[j] = pa_x[j] * tex_xyz_xxxy_0[j] - pc_x[j] * tex_xyz_xxxy_1[j] + 0.5 * fl1_fx * tex_yz_xxxy_0[j] -
                                     0.5 * fl1_fx * tex_yz_xxxy_1[j] + 1.5 * fl1_fx * tex_xyz_xxy_0[j] - 1.5 * fl1_fx * tex_xyz_xxy_1[j] +
                                     ta_xyz_xxxy_1[j];

                tey_xxyz_xxxy_0[j] = pa_x[j] * tey_xyz_xxxy_0[j] - pc_x[j] * tey_xyz_xxxy_1[j] + 0.5 * fl1_fx * tey_yz_xxxy_0[j] -
                                     0.5 * fl1_fx * tey_yz_xxxy_1[j] + 1.5 * fl1_fx * tey_xyz_xxy_0[j] - 1.5 * fl1_fx * tey_xyz_xxy_1[j];

                tez_xxyz_xxxy_0[j] = pa_x[j] * tez_xyz_xxxy_0[j] - pc_x[j] * tez_xyz_xxxy_1[j] + 0.5 * fl1_fx * tez_yz_xxxy_0[j] -
                                     0.5 * fl1_fx * tez_yz_xxxy_1[j] + 1.5 * fl1_fx * tez_xyz_xxy_0[j] - 1.5 * fl1_fx * tez_xyz_xxy_1[j];

                tex_xxyz_xxxz_0[j] = pa_x[j] * tex_xyz_xxxz_0[j] - pc_x[j] * tex_xyz_xxxz_1[j] + 0.5 * fl1_fx * tex_yz_xxxz_0[j] -
                                     0.5 * fl1_fx * tex_yz_xxxz_1[j] + 1.5 * fl1_fx * tex_xyz_xxz_0[j] - 1.5 * fl1_fx * tex_xyz_xxz_1[j] +
                                     ta_xyz_xxxz_1[j];

                tey_xxyz_xxxz_0[j] = pa_x[j] * tey_xyz_xxxz_0[j] - pc_x[j] * tey_xyz_xxxz_1[j] + 0.5 * fl1_fx * tey_yz_xxxz_0[j] -
                                     0.5 * fl1_fx * tey_yz_xxxz_1[j] + 1.5 * fl1_fx * tey_xyz_xxz_0[j] - 1.5 * fl1_fx * tey_xyz_xxz_1[j];

                tez_xxyz_xxxz_0[j] = pa_x[j] * tez_xyz_xxxz_0[j] - pc_x[j] * tez_xyz_xxxz_1[j] + 0.5 * fl1_fx * tez_yz_xxxz_0[j] -
                                     0.5 * fl1_fx * tez_yz_xxxz_1[j] + 1.5 * fl1_fx * tez_xyz_xxz_0[j] - 1.5 * fl1_fx * tez_xyz_xxz_1[j];

                tex_xxyz_xxyy_0[j] = pa_x[j] * tex_xyz_xxyy_0[j] - pc_x[j] * tex_xyz_xxyy_1[j] + 0.5 * fl1_fx * tex_yz_xxyy_0[j] -
                                     0.5 * fl1_fx * tex_yz_xxyy_1[j] + fl1_fx * tex_xyz_xyy_0[j] - fl1_fx * tex_xyz_xyy_1[j] + ta_xyz_xxyy_1[j];

                tey_xxyz_xxyy_0[j] = pa_x[j] * tey_xyz_xxyy_0[j] - pc_x[j] * tey_xyz_xxyy_1[j] + 0.5 * fl1_fx * tey_yz_xxyy_0[j] -
                                     0.5 * fl1_fx * tey_yz_xxyy_1[j] + fl1_fx * tey_xyz_xyy_0[j] - fl1_fx * tey_xyz_xyy_1[j];

                tez_xxyz_xxyy_0[j] = pa_x[j] * tez_xyz_xxyy_0[j] - pc_x[j] * tez_xyz_xxyy_1[j] + 0.5 * fl1_fx * tez_yz_xxyy_0[j] -
                                     0.5 * fl1_fx * tez_yz_xxyy_1[j] + fl1_fx * tez_xyz_xyy_0[j] - fl1_fx * tez_xyz_xyy_1[j];

                tex_xxyz_xxyz_0[j] = pa_x[j] * tex_xyz_xxyz_0[j] - pc_x[j] * tex_xyz_xxyz_1[j] + 0.5 * fl1_fx * tex_yz_xxyz_0[j] -
                                     0.5 * fl1_fx * tex_yz_xxyz_1[j] + fl1_fx * tex_xyz_xyz_0[j] - fl1_fx * tex_xyz_xyz_1[j] + ta_xyz_xxyz_1[j];

                tey_xxyz_xxyz_0[j] = pa_x[j] * tey_xyz_xxyz_0[j] - pc_x[j] * tey_xyz_xxyz_1[j] + 0.5 * fl1_fx * tey_yz_xxyz_0[j] -
                                     0.5 * fl1_fx * tey_yz_xxyz_1[j] + fl1_fx * tey_xyz_xyz_0[j] - fl1_fx * tey_xyz_xyz_1[j];

                tez_xxyz_xxyz_0[j] = pa_x[j] * tez_xyz_xxyz_0[j] - pc_x[j] * tez_xyz_xxyz_1[j] + 0.5 * fl1_fx * tez_yz_xxyz_0[j] -
                                     0.5 * fl1_fx * tez_yz_xxyz_1[j] + fl1_fx * tez_xyz_xyz_0[j] - fl1_fx * tez_xyz_xyz_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForGG_195_243(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pcDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto)
{
    // Batch of Integrals (195,243)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_4_4_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_2_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tex_xyz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 65);

            auto tey_xyz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 65);

            auto tez_xyz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 65);

            auto tex_xyz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 66);

            auto tey_xyz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 66);

            auto tez_xyz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 66);

            auto tex_xyz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 67);

            auto tey_xyz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 67);

            auto tez_xyz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 67);

            auto tex_xyz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 68);

            auto tey_xyz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 68);

            auto tez_xyz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 68);

            auto tex_xyz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 69);

            auto tey_xyz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 69);

            auto tez_xyz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 69);

            auto tex_xyz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 70);

            auto tey_xyz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 70);

            auto tez_xyz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 70);

            auto tex_xyz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 71);

            auto tey_xyz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 71);

            auto tez_xyz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 71);

            auto tex_xyz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 72);

            auto tey_xyz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 72);

            auto tez_xyz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 72);

            auto tex_xyz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 73);

            auto tey_xyz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 73);

            auto tez_xyz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 73);

            auto tex_xyz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 74);

            auto tey_xyz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 74);

            auto tez_xyz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 74);

            auto tex_xzz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 75);

            auto tey_xzz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 75);

            auto tez_xzz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 75);

            auto tex_xzz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 76);

            auto tey_xzz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 76);

            auto tez_xzz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 76);

            auto tex_xzz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 77);

            auto tey_xzz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 77);

            auto tez_xzz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 77);

            auto tex_xzz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 78);

            auto tey_xzz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 78);

            auto tez_xzz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 78);

            auto tex_xzz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 79);

            auto tey_xzz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 79);

            auto tez_xzz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 79);

            auto tex_xzz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 80);

            auto tey_xzz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 80);

            auto tez_xzz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 80);

            auto tex_xyz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 65);

            auto tey_xyz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 65);

            auto tez_xyz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 65);

            auto tex_xyz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 66);

            auto tey_xyz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 66);

            auto tez_xyz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 66);

            auto tex_xyz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 67);

            auto tey_xyz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 67);

            auto tez_xyz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 67);

            auto tex_xyz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 68);

            auto tey_xyz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 68);

            auto tez_xyz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 68);

            auto tex_xyz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 69);

            auto tey_xyz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 69);

            auto tez_xyz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 69);

            auto tex_xyz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 70);

            auto tey_xyz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 70);

            auto tez_xyz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 70);

            auto tex_xyz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 71);

            auto tey_xyz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 71);

            auto tez_xyz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 71);

            auto tex_xyz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 72);

            auto tey_xyz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 72);

            auto tez_xyz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 72);

            auto tex_xyz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 73);

            auto tey_xyz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 73);

            auto tez_xyz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 73);

            auto tex_xyz_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 74);

            auto tey_xyz_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 74);

            auto tez_xyz_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 74);

            auto tex_xzz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 75);

            auto tey_xzz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 75);

            auto tez_xzz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 75);

            auto tex_xzz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 76);

            auto tey_xzz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 76);

            auto tez_xzz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 76);

            auto tex_xzz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 77);

            auto tey_xzz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 77);

            auto tez_xzz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 77);

            auto tex_xzz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 78);

            auto tey_xzz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 78);

            auto tez_xzz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 78);

            auto tex_xzz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 79);

            auto tey_xzz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 79);

            auto tez_xzz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 79);

            auto tex_xzz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 80);

            auto tey_xzz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 80);

            auto tez_xzz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 80);

            auto tex_yz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 65);

            auto tey_yz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 65);

            auto tez_yz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 65);

            auto tex_yz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 66);

            auto tey_yz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 66);

            auto tez_yz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 66);

            auto tex_yz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 67);

            auto tey_yz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 67);

            auto tez_yz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 67);

            auto tex_yz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 68);

            auto tey_yz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 68);

            auto tez_yz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 68);

            auto tex_yz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 69);

            auto tey_yz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 69);

            auto tez_yz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 69);

            auto tex_yz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 70);

            auto tey_yz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 70);

            auto tez_yz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 70);

            auto tex_yz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 71);

            auto tey_yz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 71);

            auto tez_yz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 71);

            auto tex_yz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 72);

            auto tey_yz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 72);

            auto tez_yz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 72);

            auto tex_yz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 73);

            auto tey_yz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 73);

            auto tez_yz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 73);

            auto tex_yz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 74);

            auto tey_yz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 74);

            auto tez_yz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 74);

            auto tex_zz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 75);

            auto tey_zz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 75);

            auto tez_zz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 75);

            auto tex_zz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 76);

            auto tey_zz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 76);

            auto tez_zz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 76);

            auto tex_zz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 77);

            auto tey_zz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 77);

            auto tez_zz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 77);

            auto tex_zz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 78);

            auto tey_zz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 78);

            auto tez_zz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 78);

            auto tex_zz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 79);

            auto tey_zz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 79);

            auto tez_zz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 79);

            auto tex_zz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 80);

            auto tey_zz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 80);

            auto tez_zz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 80);

            auto tex_yz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 65);

            auto tey_yz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 65);

            auto tez_yz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 65);

            auto tex_yz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 66);

            auto tey_yz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 66);

            auto tez_yz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 66);

            auto tex_yz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 67);

            auto tey_yz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 67);

            auto tez_yz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 67);

            auto tex_yz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 68);

            auto tey_yz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 68);

            auto tez_yz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 68);

            auto tex_yz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 69);

            auto tey_yz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 69);

            auto tez_yz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 69);

            auto tex_yz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 70);

            auto tey_yz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 70);

            auto tez_yz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 70);

            auto tex_yz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 71);

            auto tey_yz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 71);

            auto tez_yz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 71);

            auto tex_yz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 72);

            auto tey_yz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 72);

            auto tez_yz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 72);

            auto tex_yz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 73);

            auto tey_yz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 73);

            auto tez_yz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 73);

            auto tex_yz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 74);

            auto tey_yz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 74);

            auto tez_yz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 74);

            auto tex_zz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 75);

            auto tey_zz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 75);

            auto tez_zz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 75);

            auto tex_zz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 76);

            auto tey_zz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 76);

            auto tez_zz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 76);

            auto tex_zz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 77);

            auto tey_zz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 77);

            auto tez_zz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 77);

            auto tex_zz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 78);

            auto tey_zz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 78);

            auto tez_zz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 78);

            auto tex_zz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 79);

            auto tey_zz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 79);

            auto tez_zz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 79);

            auto tex_zz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 80);

            auto tey_zz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 80);

            auto tez_zz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 80);

            auto tex_xyz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 45);

            auto tey_xyz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 45);

            auto tez_xyz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 45);

            auto tex_xyz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 46);

            auto tey_xyz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 46);

            auto tez_xyz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 46);

            auto tex_xyz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 47);

            auto tey_xyz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 47);

            auto tez_xyz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 47);

            auto tex_xyz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 48);

            auto tey_xyz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 48);

            auto tez_xyz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 48);

            auto tex_xyz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 49);

            auto tey_xyz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 49);

            auto tez_xyz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 49);

            auto tex_xzz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 50);

            auto tey_xzz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 50);

            auto tez_xzz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 50);

            auto tex_xzz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 51);

            auto tey_xzz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 51);

            auto tez_xzz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 51);

            auto tex_xzz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 52);

            auto tey_xzz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 52);

            auto tez_xzz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 52);

            auto tex_xzz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 53);

            auto tey_xzz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 53);

            auto tez_xzz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 53);

            auto tex_xzz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 54);

            auto tey_xzz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 54);

            auto tez_xzz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 54);

            auto tex_xzz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 55);

            auto tey_xzz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 55);

            auto tez_xzz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 55);

            auto tex_xyz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 45);

            auto tey_xyz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 45);

            auto tez_xyz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 45);

            auto tex_xyz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 46);

            auto tey_xyz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 46);

            auto tez_xyz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 46);

            auto tex_xyz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 47);

            auto tey_xyz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 47);

            auto tez_xyz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 47);

            auto tex_xyz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 48);

            auto tey_xyz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 48);

            auto tez_xyz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 48);

            auto tex_xyz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 49);

            auto tey_xyz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 49);

            auto tez_xyz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 49);

            auto tex_xzz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 50);

            auto tey_xzz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 50);

            auto tez_xzz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 50);

            auto tex_xzz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 51);

            auto tey_xzz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 51);

            auto tez_xzz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 51);

            auto tex_xzz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 52);

            auto tey_xzz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 52);

            auto tez_xzz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 52);

            auto tex_xzz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 53);

            auto tey_xzz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 53);

            auto tez_xzz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 53);

            auto tex_xzz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 54);

            auto tey_xzz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 54);

            auto tez_xzz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 54);

            auto tex_xzz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 55);

            auto tey_xzz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 55);

            auto tez_xzz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 55);

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

            // set up pointers to integrals

            auto tex_xxyz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 65);

            auto tey_xxyz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 65);

            auto tez_xxyz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 65);

            auto tex_xxyz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 66);

            auto tey_xxyz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 66);

            auto tez_xxyz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 66);

            auto tex_xxyz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 67);

            auto tey_xxyz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 67);

            auto tez_xxyz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 67);

            auto tex_xxyz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 68);

            auto tey_xxyz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 68);

            auto tez_xxyz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 68);

            auto tex_xxyz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 69);

            auto tey_xxyz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 69);

            auto tez_xxyz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 69);

            auto tex_xxyz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 70);

            auto tey_xxyz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 70);

            auto tez_xxyz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 70);

            auto tex_xxyz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 71);

            auto tey_xxyz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 71);

            auto tez_xxyz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 71);

            auto tex_xxyz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 72);

            auto tey_xxyz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 72);

            auto tez_xxyz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 72);

            auto tex_xxyz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 73);

            auto tey_xxyz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 73);

            auto tez_xxyz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 73);

            auto tex_xxyz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 74);

            auto tey_xxyz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 74);

            auto tez_xxyz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 74);

            auto tex_xxzz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 75);

            auto tey_xxzz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 75);

            auto tez_xxzz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 75);

            auto tex_xxzz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 76);

            auto tey_xxzz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 76);

            auto tez_xxzz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 76);

            auto tex_xxzz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 77);

            auto tey_xxzz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 77);

            auto tez_xxzz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 77);

            auto tex_xxzz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 78);

            auto tey_xxzz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 78);

            auto tez_xxzz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 78);

            auto tex_xxzz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 79);

            auto tey_xxzz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 79);

            auto tez_xxzz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 79);

            auto tex_xxzz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 80);

            auto tey_xxzz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 80);

            auto tez_xxzz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 80);

            // Batch of Integrals (195,243)

            #pragma omp simd aligned(fx, pa_x, pc_x, ta_xyz_xxzz_1, ta_xyz_xyyy_1, ta_xyz_xyyz_1, \
                                         ta_xyz_xyzz_1, ta_xyz_xzzz_1, ta_xyz_yyyy_1, ta_xyz_yyyz_1, ta_xyz_yyzz_1, \
                                         ta_xyz_yzzz_1, ta_xyz_zzzz_1, ta_xzz_xxxx_1, ta_xzz_xxxy_1, ta_xzz_xxxz_1, \
                                         ta_xzz_xxyy_1, ta_xzz_xxyz_1, ta_xzz_xxzz_1, tex_xxyz_xxzz_0, tex_xxyz_xyyy_0, \
                                         tex_xxyz_xyyz_0, tex_xxyz_xyzz_0, tex_xxyz_xzzz_0, tex_xxyz_yyyy_0, tex_xxyz_yyyz_0, \
                                         tex_xxyz_yyzz_0, tex_xxyz_yzzz_0, tex_xxyz_zzzz_0, tex_xxzz_xxxx_0, tex_xxzz_xxxy_0, \
                                         tex_xxzz_xxxz_0, tex_xxzz_xxyy_0, tex_xxzz_xxyz_0, tex_xxzz_xxzz_0, tex_xyz_xxzz_0, \
                                         tex_xyz_xxzz_1, tex_xyz_xyyy_0, tex_xyz_xyyy_1, tex_xyz_xyyz_0, tex_xyz_xyyz_1, \
                                         tex_xyz_xyzz_0, tex_xyz_xyzz_1, tex_xyz_xzz_0, tex_xyz_xzz_1, tex_xyz_xzzz_0, \
                                         tex_xyz_xzzz_1, tex_xyz_yyy_0, tex_xyz_yyy_1, tex_xyz_yyyy_0, tex_xyz_yyyy_1, \
                                         tex_xyz_yyyz_0, tex_xyz_yyyz_1, tex_xyz_yyz_0, tex_xyz_yyz_1, tex_xyz_yyzz_0, \
                                         tex_xyz_yyzz_1, tex_xyz_yzz_0, tex_xyz_yzz_1, tex_xyz_yzzz_0, tex_xyz_yzzz_1, \
                                         tex_xyz_zzz_0, tex_xyz_zzz_1, tex_xyz_zzzz_0, tex_xyz_zzzz_1, tex_xzz_xxx_0, \
                                         tex_xzz_xxx_1, tex_xzz_xxxx_0, tex_xzz_xxxx_1, tex_xzz_xxxy_0, tex_xzz_xxxy_1, \
                                         tex_xzz_xxxz_0, tex_xzz_xxxz_1, tex_xzz_xxy_0, tex_xzz_xxy_1, tex_xzz_xxyy_0, \
                                         tex_xzz_xxyy_1, tex_xzz_xxyz_0, tex_xzz_xxyz_1, tex_xzz_xxz_0, tex_xzz_xxz_1, \
                                         tex_xzz_xxzz_0, tex_xzz_xxzz_1, tex_xzz_xyy_0, tex_xzz_xyy_1, tex_xzz_xyz_0, \
                                         tex_xzz_xyz_1, tex_xzz_xzz_0, tex_xzz_xzz_1, tex_yz_xxzz_0, tex_yz_xxzz_1, \
                                         tex_yz_xyyy_0, tex_yz_xyyy_1, tex_yz_xyyz_0, tex_yz_xyyz_1, tex_yz_xyzz_0, \
                                         tex_yz_xyzz_1, tex_yz_xzzz_0, tex_yz_xzzz_1, tex_yz_yyyy_0, tex_yz_yyyy_1, \
                                         tex_yz_yyyz_0, tex_yz_yyyz_1, tex_yz_yyzz_0, tex_yz_yyzz_1, tex_yz_yzzz_0, \
                                         tex_yz_yzzz_1, tex_yz_zzzz_0, tex_yz_zzzz_1, tex_zz_xxxx_0, tex_zz_xxxx_1, \
                                         tex_zz_xxxy_0, tex_zz_xxxy_1, tex_zz_xxxz_0, tex_zz_xxxz_1, tex_zz_xxyy_0, \
                                         tex_zz_xxyy_1, tex_zz_xxyz_0, tex_zz_xxyz_1, tex_zz_xxzz_0, tex_zz_xxzz_1, \
                                         tey_xxyz_xxzz_0, tey_xxyz_xyyy_0, tey_xxyz_xyyz_0, tey_xxyz_xyzz_0, tey_xxyz_xzzz_0, \
                                         tey_xxyz_yyyy_0, tey_xxyz_yyyz_0, tey_xxyz_yyzz_0, tey_xxyz_yzzz_0, tey_xxyz_zzzz_0, \
                                         tey_xxzz_xxxx_0, tey_xxzz_xxxy_0, tey_xxzz_xxxz_0, tey_xxzz_xxyy_0, tey_xxzz_xxyz_0, \
                                         tey_xxzz_xxzz_0, tey_xyz_xxzz_0, tey_xyz_xxzz_1, tey_xyz_xyyy_0, tey_xyz_xyyy_1, \
                                         tey_xyz_xyyz_0, tey_xyz_xyyz_1, tey_xyz_xyzz_0, tey_xyz_xyzz_1, tey_xyz_xzz_0, \
                                         tey_xyz_xzz_1, tey_xyz_xzzz_0, tey_xyz_xzzz_1, tey_xyz_yyy_0, tey_xyz_yyy_1, \
                                         tey_xyz_yyyy_0, tey_xyz_yyyy_1, tey_xyz_yyyz_0, tey_xyz_yyyz_1, tey_xyz_yyz_0, \
                                         tey_xyz_yyz_1, tey_xyz_yyzz_0, tey_xyz_yyzz_1, tey_xyz_yzz_0, tey_xyz_yzz_1, \
                                         tey_xyz_yzzz_0, tey_xyz_yzzz_1, tey_xyz_zzz_0, tey_xyz_zzz_1, tey_xyz_zzzz_0, \
                                         tey_xyz_zzzz_1, tey_xzz_xxx_0, tey_xzz_xxx_1, tey_xzz_xxxx_0, tey_xzz_xxxx_1, \
                                         tey_xzz_xxxy_0, tey_xzz_xxxy_1, tey_xzz_xxxz_0, tey_xzz_xxxz_1, tey_xzz_xxy_0, \
                                         tey_xzz_xxy_1, tey_xzz_xxyy_0, tey_xzz_xxyy_1, tey_xzz_xxyz_0, tey_xzz_xxyz_1, \
                                         tey_xzz_xxz_0, tey_xzz_xxz_1, tey_xzz_xxzz_0, tey_xzz_xxzz_1, tey_xzz_xyy_0, \
                                         tey_xzz_xyy_1, tey_xzz_xyz_0, tey_xzz_xyz_1, tey_xzz_xzz_0, tey_xzz_xzz_1, \
                                         tey_yz_xxzz_0, tey_yz_xxzz_1, tey_yz_xyyy_0, tey_yz_xyyy_1, tey_yz_xyyz_0, \
                                         tey_yz_xyyz_1, tey_yz_xyzz_0, tey_yz_xyzz_1, tey_yz_xzzz_0, tey_yz_xzzz_1, \
                                         tey_yz_yyyy_0, tey_yz_yyyy_1, tey_yz_yyyz_0, tey_yz_yyyz_1, tey_yz_yyzz_0, \
                                         tey_yz_yyzz_1, tey_yz_yzzz_0, tey_yz_yzzz_1, tey_yz_zzzz_0, tey_yz_zzzz_1, \
                                         tey_zz_xxxx_0, tey_zz_xxxx_1, tey_zz_xxxy_0, tey_zz_xxxy_1, tey_zz_xxxz_0, \
                                         tey_zz_xxxz_1, tey_zz_xxyy_0, tey_zz_xxyy_1, tey_zz_xxyz_0, tey_zz_xxyz_1, \
                                         tey_zz_xxzz_0, tey_zz_xxzz_1, tez_xxyz_xxzz_0, tez_xxyz_xyyy_0, tez_xxyz_xyyz_0, \
                                         tez_xxyz_xyzz_0, tez_xxyz_xzzz_0, tez_xxyz_yyyy_0, tez_xxyz_yyyz_0, tez_xxyz_yyzz_0, \
                                         tez_xxyz_yzzz_0, tez_xxyz_zzzz_0, tez_xxzz_xxxx_0, tez_xxzz_xxxy_0, tez_xxzz_xxxz_0, \
                                         tez_xxzz_xxyy_0, tez_xxzz_xxyz_0, tez_xxzz_xxzz_0, tez_xyz_xxzz_0, tez_xyz_xxzz_1, \
                                         tez_xyz_xyyy_0, tez_xyz_xyyy_1, tez_xyz_xyyz_0, tez_xyz_xyyz_1, tez_xyz_xyzz_0, \
                                         tez_xyz_xyzz_1, tez_xyz_xzz_0, tez_xyz_xzz_1, tez_xyz_xzzz_0, tez_xyz_xzzz_1, \
                                         tez_xyz_yyy_0, tez_xyz_yyy_1, tez_xyz_yyyy_0, tez_xyz_yyyy_1, tez_xyz_yyyz_0, \
                                         tez_xyz_yyyz_1, tez_xyz_yyz_0, tez_xyz_yyz_1, tez_xyz_yyzz_0, tez_xyz_yyzz_1, \
                                         tez_xyz_yzz_0, tez_xyz_yzz_1, tez_xyz_yzzz_0, tez_xyz_yzzz_1, tez_xyz_zzz_0, \
                                         tez_xyz_zzz_1, tez_xyz_zzzz_0, tez_xyz_zzzz_1, tez_xzz_xxx_0, tez_xzz_xxx_1, \
                                         tez_xzz_xxxx_0, tez_xzz_xxxx_1, tez_xzz_xxxy_0, tez_xzz_xxxy_1, tez_xzz_xxxz_0, \
                                         tez_xzz_xxxz_1, tez_xzz_xxy_0, tez_xzz_xxy_1, tez_xzz_xxyy_0, tez_xzz_xxyy_1, \
                                         tez_xzz_xxyz_0, tez_xzz_xxyz_1, tez_xzz_xxz_0, tez_xzz_xxz_1, tez_xzz_xxzz_0, \
                                         tez_xzz_xxzz_1, tez_xzz_xyy_0, tez_xzz_xyy_1, tez_xzz_xyz_0, tez_xzz_xyz_1, \
                                         tez_xzz_xzz_0, tez_xzz_xzz_1, tez_yz_xxzz_0, tez_yz_xxzz_1, tez_yz_xyyy_0, \
                                         tez_yz_xyyy_1, tez_yz_xyyz_0, tez_yz_xyyz_1, tez_yz_xyzz_0, tez_yz_xyzz_1, \
                                         tez_yz_xzzz_0, tez_yz_xzzz_1, tez_yz_yyyy_0, tez_yz_yyyy_1, tez_yz_yyyz_0, \
                                         tez_yz_yyyz_1, tez_yz_yyzz_0, tez_yz_yyzz_1, tez_yz_yzzz_0, tez_yz_yzzz_1, \
                                         tez_yz_zzzz_0, tez_yz_zzzz_1, tez_zz_xxxx_0, tez_zz_xxxx_1, tez_zz_xxxy_0, \
                                         tez_zz_xxxy_1, tez_zz_xxxz_0, tez_zz_xxxz_1, tez_zz_xxyy_0, tez_zz_xxyy_1, \
                                         tez_zz_xxyz_0, tez_zz_xxyz_1, tez_zz_xxzz_0, tez_zz_xxzz_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_xxyz_xxzz_0[j] = pa_x[j] * tex_xyz_xxzz_0[j] - pc_x[j] * tex_xyz_xxzz_1[j] + 0.5 * fl1_fx * tex_yz_xxzz_0[j] -
                                     0.5 * fl1_fx * tex_yz_xxzz_1[j] + fl1_fx * tex_xyz_xzz_0[j] - fl1_fx * tex_xyz_xzz_1[j] + ta_xyz_xxzz_1[j];

                tey_xxyz_xxzz_0[j] = pa_x[j] * tey_xyz_xxzz_0[j] - pc_x[j] * tey_xyz_xxzz_1[j] + 0.5 * fl1_fx * tey_yz_xxzz_0[j] -
                                     0.5 * fl1_fx * tey_yz_xxzz_1[j] + fl1_fx * tey_xyz_xzz_0[j] - fl1_fx * tey_xyz_xzz_1[j];

                tez_xxyz_xxzz_0[j] = pa_x[j] * tez_xyz_xxzz_0[j] - pc_x[j] * tez_xyz_xxzz_1[j] + 0.5 * fl1_fx * tez_yz_xxzz_0[j] -
                                     0.5 * fl1_fx * tez_yz_xxzz_1[j] + fl1_fx * tez_xyz_xzz_0[j] - fl1_fx * tez_xyz_xzz_1[j];

                tex_xxyz_xyyy_0[j] = pa_x[j] * tex_xyz_xyyy_0[j] - pc_x[j] * tex_xyz_xyyy_1[j] + 0.5 * fl1_fx * tex_yz_xyyy_0[j] -
                                     0.5 * fl1_fx * tex_yz_xyyy_1[j] + 0.5 * fl1_fx * tex_xyz_yyy_0[j] - 0.5 * fl1_fx * tex_xyz_yyy_1[j] +
                                     ta_xyz_xyyy_1[j];

                tey_xxyz_xyyy_0[j] = pa_x[j] * tey_xyz_xyyy_0[j] - pc_x[j] * tey_xyz_xyyy_1[j] + 0.5 * fl1_fx * tey_yz_xyyy_0[j] -
                                     0.5 * fl1_fx * tey_yz_xyyy_1[j] + 0.5 * fl1_fx * tey_xyz_yyy_0[j] - 0.5 * fl1_fx * tey_xyz_yyy_1[j];

                tez_xxyz_xyyy_0[j] = pa_x[j] * tez_xyz_xyyy_0[j] - pc_x[j] * tez_xyz_xyyy_1[j] + 0.5 * fl1_fx * tez_yz_xyyy_0[j] -
                                     0.5 * fl1_fx * tez_yz_xyyy_1[j] + 0.5 * fl1_fx * tez_xyz_yyy_0[j] - 0.5 * fl1_fx * tez_xyz_yyy_1[j];

                tex_xxyz_xyyz_0[j] = pa_x[j] * tex_xyz_xyyz_0[j] - pc_x[j] * tex_xyz_xyyz_1[j] + 0.5 * fl1_fx * tex_yz_xyyz_0[j] -
                                     0.5 * fl1_fx * tex_yz_xyyz_1[j] + 0.5 * fl1_fx * tex_xyz_yyz_0[j] - 0.5 * fl1_fx * tex_xyz_yyz_1[j] +
                                     ta_xyz_xyyz_1[j];

                tey_xxyz_xyyz_0[j] = pa_x[j] * tey_xyz_xyyz_0[j] - pc_x[j] * tey_xyz_xyyz_1[j] + 0.5 * fl1_fx * tey_yz_xyyz_0[j] -
                                     0.5 * fl1_fx * tey_yz_xyyz_1[j] + 0.5 * fl1_fx * tey_xyz_yyz_0[j] - 0.5 * fl1_fx * tey_xyz_yyz_1[j];

                tez_xxyz_xyyz_0[j] = pa_x[j] * tez_xyz_xyyz_0[j] - pc_x[j] * tez_xyz_xyyz_1[j] + 0.5 * fl1_fx * tez_yz_xyyz_0[j] -
                                     0.5 * fl1_fx * tez_yz_xyyz_1[j] + 0.5 * fl1_fx * tez_xyz_yyz_0[j] - 0.5 * fl1_fx * tez_xyz_yyz_1[j];

                tex_xxyz_xyzz_0[j] = pa_x[j] * tex_xyz_xyzz_0[j] - pc_x[j] * tex_xyz_xyzz_1[j] + 0.5 * fl1_fx * tex_yz_xyzz_0[j] -
                                     0.5 * fl1_fx * tex_yz_xyzz_1[j] + 0.5 * fl1_fx * tex_xyz_yzz_0[j] - 0.5 * fl1_fx * tex_xyz_yzz_1[j] +
                                     ta_xyz_xyzz_1[j];

                tey_xxyz_xyzz_0[j] = pa_x[j] * tey_xyz_xyzz_0[j] - pc_x[j] * tey_xyz_xyzz_1[j] + 0.5 * fl1_fx * tey_yz_xyzz_0[j] -
                                     0.5 * fl1_fx * tey_yz_xyzz_1[j] + 0.5 * fl1_fx * tey_xyz_yzz_0[j] - 0.5 * fl1_fx * tey_xyz_yzz_1[j];

                tez_xxyz_xyzz_0[j] = pa_x[j] * tez_xyz_xyzz_0[j] - pc_x[j] * tez_xyz_xyzz_1[j] + 0.5 * fl1_fx * tez_yz_xyzz_0[j] -
                                     0.5 * fl1_fx * tez_yz_xyzz_1[j] + 0.5 * fl1_fx * tez_xyz_yzz_0[j] - 0.5 * fl1_fx * tez_xyz_yzz_1[j];

                tex_xxyz_xzzz_0[j] = pa_x[j] * tex_xyz_xzzz_0[j] - pc_x[j] * tex_xyz_xzzz_1[j] + 0.5 * fl1_fx * tex_yz_xzzz_0[j] -
                                     0.5 * fl1_fx * tex_yz_xzzz_1[j] + 0.5 * fl1_fx * tex_xyz_zzz_0[j] - 0.5 * fl1_fx * tex_xyz_zzz_1[j] +
                                     ta_xyz_xzzz_1[j];

                tey_xxyz_xzzz_0[j] = pa_x[j] * tey_xyz_xzzz_0[j] - pc_x[j] * tey_xyz_xzzz_1[j] + 0.5 * fl1_fx * tey_yz_xzzz_0[j] -
                                     0.5 * fl1_fx * tey_yz_xzzz_1[j] + 0.5 * fl1_fx * tey_xyz_zzz_0[j] - 0.5 * fl1_fx * tey_xyz_zzz_1[j];

                tez_xxyz_xzzz_0[j] = pa_x[j] * tez_xyz_xzzz_0[j] - pc_x[j] * tez_xyz_xzzz_1[j] + 0.5 * fl1_fx * tez_yz_xzzz_0[j] -
                                     0.5 * fl1_fx * tez_yz_xzzz_1[j] + 0.5 * fl1_fx * tez_xyz_zzz_0[j] - 0.5 * fl1_fx * tez_xyz_zzz_1[j];

                tex_xxyz_yyyy_0[j] = pa_x[j] * tex_xyz_yyyy_0[j] - pc_x[j] * tex_xyz_yyyy_1[j] + 0.5 * fl1_fx * tex_yz_yyyy_0[j] -
                                     0.5 * fl1_fx * tex_yz_yyyy_1[j] + ta_xyz_yyyy_1[j];

                tey_xxyz_yyyy_0[j] =
                    pa_x[j] * tey_xyz_yyyy_0[j] - pc_x[j] * tey_xyz_yyyy_1[j] + 0.5 * fl1_fx * tey_yz_yyyy_0[j] - 0.5 * fl1_fx * tey_yz_yyyy_1[j];

                tez_xxyz_yyyy_0[j] =
                    pa_x[j] * tez_xyz_yyyy_0[j] - pc_x[j] * tez_xyz_yyyy_1[j] + 0.5 * fl1_fx * tez_yz_yyyy_0[j] - 0.5 * fl1_fx * tez_yz_yyyy_1[j];

                tex_xxyz_yyyz_0[j] = pa_x[j] * tex_xyz_yyyz_0[j] - pc_x[j] * tex_xyz_yyyz_1[j] + 0.5 * fl1_fx * tex_yz_yyyz_0[j] -
                                     0.5 * fl1_fx * tex_yz_yyyz_1[j] + ta_xyz_yyyz_1[j];

                tey_xxyz_yyyz_0[j] =
                    pa_x[j] * tey_xyz_yyyz_0[j] - pc_x[j] * tey_xyz_yyyz_1[j] + 0.5 * fl1_fx * tey_yz_yyyz_0[j] - 0.5 * fl1_fx * tey_yz_yyyz_1[j];

                tez_xxyz_yyyz_0[j] =
                    pa_x[j] * tez_xyz_yyyz_0[j] - pc_x[j] * tez_xyz_yyyz_1[j] + 0.5 * fl1_fx * tez_yz_yyyz_0[j] - 0.5 * fl1_fx * tez_yz_yyyz_1[j];

                tex_xxyz_yyzz_0[j] = pa_x[j] * tex_xyz_yyzz_0[j] - pc_x[j] * tex_xyz_yyzz_1[j] + 0.5 * fl1_fx * tex_yz_yyzz_0[j] -
                                     0.5 * fl1_fx * tex_yz_yyzz_1[j] + ta_xyz_yyzz_1[j];

                tey_xxyz_yyzz_0[j] =
                    pa_x[j] * tey_xyz_yyzz_0[j] - pc_x[j] * tey_xyz_yyzz_1[j] + 0.5 * fl1_fx * tey_yz_yyzz_0[j] - 0.5 * fl1_fx * tey_yz_yyzz_1[j];

                tez_xxyz_yyzz_0[j] =
                    pa_x[j] * tez_xyz_yyzz_0[j] - pc_x[j] * tez_xyz_yyzz_1[j] + 0.5 * fl1_fx * tez_yz_yyzz_0[j] - 0.5 * fl1_fx * tez_yz_yyzz_1[j];

                tex_xxyz_yzzz_0[j] = pa_x[j] * tex_xyz_yzzz_0[j] - pc_x[j] * tex_xyz_yzzz_1[j] + 0.5 * fl1_fx * tex_yz_yzzz_0[j] -
                                     0.5 * fl1_fx * tex_yz_yzzz_1[j] + ta_xyz_yzzz_1[j];

                tey_xxyz_yzzz_0[j] =
                    pa_x[j] * tey_xyz_yzzz_0[j] - pc_x[j] * tey_xyz_yzzz_1[j] + 0.5 * fl1_fx * tey_yz_yzzz_0[j] - 0.5 * fl1_fx * tey_yz_yzzz_1[j];

                tez_xxyz_yzzz_0[j] =
                    pa_x[j] * tez_xyz_yzzz_0[j] - pc_x[j] * tez_xyz_yzzz_1[j] + 0.5 * fl1_fx * tez_yz_yzzz_0[j] - 0.5 * fl1_fx * tez_yz_yzzz_1[j];

                tex_xxyz_zzzz_0[j] = pa_x[j] * tex_xyz_zzzz_0[j] - pc_x[j] * tex_xyz_zzzz_1[j] + 0.5 * fl1_fx * tex_yz_zzzz_0[j] -
                                     0.5 * fl1_fx * tex_yz_zzzz_1[j] + ta_xyz_zzzz_1[j];

                tey_xxyz_zzzz_0[j] =
                    pa_x[j] * tey_xyz_zzzz_0[j] - pc_x[j] * tey_xyz_zzzz_1[j] + 0.5 * fl1_fx * tey_yz_zzzz_0[j] - 0.5 * fl1_fx * tey_yz_zzzz_1[j];

                tez_xxyz_zzzz_0[j] =
                    pa_x[j] * tez_xyz_zzzz_0[j] - pc_x[j] * tez_xyz_zzzz_1[j] + 0.5 * fl1_fx * tez_yz_zzzz_0[j] - 0.5 * fl1_fx * tez_yz_zzzz_1[j];

                tex_xxzz_xxxx_0[j] = pa_x[j] * tex_xzz_xxxx_0[j] - pc_x[j] * tex_xzz_xxxx_1[j] + 0.5 * fl1_fx * tex_zz_xxxx_0[j] -
                                     0.5 * fl1_fx * tex_zz_xxxx_1[j] + 2.0 * fl1_fx * tex_xzz_xxx_0[j] - 2.0 * fl1_fx * tex_xzz_xxx_1[j] +
                                     ta_xzz_xxxx_1[j];

                tey_xxzz_xxxx_0[j] = pa_x[j] * tey_xzz_xxxx_0[j] - pc_x[j] * tey_xzz_xxxx_1[j] + 0.5 * fl1_fx * tey_zz_xxxx_0[j] -
                                     0.5 * fl1_fx * tey_zz_xxxx_1[j] + 2.0 * fl1_fx * tey_xzz_xxx_0[j] - 2.0 * fl1_fx * tey_xzz_xxx_1[j];

                tez_xxzz_xxxx_0[j] = pa_x[j] * tez_xzz_xxxx_0[j] - pc_x[j] * tez_xzz_xxxx_1[j] + 0.5 * fl1_fx * tez_zz_xxxx_0[j] -
                                     0.5 * fl1_fx * tez_zz_xxxx_1[j] + 2.0 * fl1_fx * tez_xzz_xxx_0[j] - 2.0 * fl1_fx * tez_xzz_xxx_1[j];

                tex_xxzz_xxxy_0[j] = pa_x[j] * tex_xzz_xxxy_0[j] - pc_x[j] * tex_xzz_xxxy_1[j] + 0.5 * fl1_fx * tex_zz_xxxy_0[j] -
                                     0.5 * fl1_fx * tex_zz_xxxy_1[j] + 1.5 * fl1_fx * tex_xzz_xxy_0[j] - 1.5 * fl1_fx * tex_xzz_xxy_1[j] +
                                     ta_xzz_xxxy_1[j];

                tey_xxzz_xxxy_0[j] = pa_x[j] * tey_xzz_xxxy_0[j] - pc_x[j] * tey_xzz_xxxy_1[j] + 0.5 * fl1_fx * tey_zz_xxxy_0[j] -
                                     0.5 * fl1_fx * tey_zz_xxxy_1[j] + 1.5 * fl1_fx * tey_xzz_xxy_0[j] - 1.5 * fl1_fx * tey_xzz_xxy_1[j];

                tez_xxzz_xxxy_0[j] = pa_x[j] * tez_xzz_xxxy_0[j] - pc_x[j] * tez_xzz_xxxy_1[j] + 0.5 * fl1_fx * tez_zz_xxxy_0[j] -
                                     0.5 * fl1_fx * tez_zz_xxxy_1[j] + 1.5 * fl1_fx * tez_xzz_xxy_0[j] - 1.5 * fl1_fx * tez_xzz_xxy_1[j];

                tex_xxzz_xxxz_0[j] = pa_x[j] * tex_xzz_xxxz_0[j] - pc_x[j] * tex_xzz_xxxz_1[j] + 0.5 * fl1_fx * tex_zz_xxxz_0[j] -
                                     0.5 * fl1_fx * tex_zz_xxxz_1[j] + 1.5 * fl1_fx * tex_xzz_xxz_0[j] - 1.5 * fl1_fx * tex_xzz_xxz_1[j] +
                                     ta_xzz_xxxz_1[j];

                tey_xxzz_xxxz_0[j] = pa_x[j] * tey_xzz_xxxz_0[j] - pc_x[j] * tey_xzz_xxxz_1[j] + 0.5 * fl1_fx * tey_zz_xxxz_0[j] -
                                     0.5 * fl1_fx * tey_zz_xxxz_1[j] + 1.5 * fl1_fx * tey_xzz_xxz_0[j] - 1.5 * fl1_fx * tey_xzz_xxz_1[j];

                tez_xxzz_xxxz_0[j] = pa_x[j] * tez_xzz_xxxz_0[j] - pc_x[j] * tez_xzz_xxxz_1[j] + 0.5 * fl1_fx * tez_zz_xxxz_0[j] -
                                     0.5 * fl1_fx * tez_zz_xxxz_1[j] + 1.5 * fl1_fx * tez_xzz_xxz_0[j] - 1.5 * fl1_fx * tez_xzz_xxz_1[j];

                tex_xxzz_xxyy_0[j] = pa_x[j] * tex_xzz_xxyy_0[j] - pc_x[j] * tex_xzz_xxyy_1[j] + 0.5 * fl1_fx * tex_zz_xxyy_0[j] -
                                     0.5 * fl1_fx * tex_zz_xxyy_1[j] + fl1_fx * tex_xzz_xyy_0[j] - fl1_fx * tex_xzz_xyy_1[j] + ta_xzz_xxyy_1[j];

                tey_xxzz_xxyy_0[j] = pa_x[j] * tey_xzz_xxyy_0[j] - pc_x[j] * tey_xzz_xxyy_1[j] + 0.5 * fl1_fx * tey_zz_xxyy_0[j] -
                                     0.5 * fl1_fx * tey_zz_xxyy_1[j] + fl1_fx * tey_xzz_xyy_0[j] - fl1_fx * tey_xzz_xyy_1[j];

                tez_xxzz_xxyy_0[j] = pa_x[j] * tez_xzz_xxyy_0[j] - pc_x[j] * tez_xzz_xxyy_1[j] + 0.5 * fl1_fx * tez_zz_xxyy_0[j] -
                                     0.5 * fl1_fx * tez_zz_xxyy_1[j] + fl1_fx * tez_xzz_xyy_0[j] - fl1_fx * tez_xzz_xyy_1[j];

                tex_xxzz_xxyz_0[j] = pa_x[j] * tex_xzz_xxyz_0[j] - pc_x[j] * tex_xzz_xxyz_1[j] + 0.5 * fl1_fx * tex_zz_xxyz_0[j] -
                                     0.5 * fl1_fx * tex_zz_xxyz_1[j] + fl1_fx * tex_xzz_xyz_0[j] - fl1_fx * tex_xzz_xyz_1[j] + ta_xzz_xxyz_1[j];

                tey_xxzz_xxyz_0[j] = pa_x[j] * tey_xzz_xxyz_0[j] - pc_x[j] * tey_xzz_xxyz_1[j] + 0.5 * fl1_fx * tey_zz_xxyz_0[j] -
                                     0.5 * fl1_fx * tey_zz_xxyz_1[j] + fl1_fx * tey_xzz_xyz_0[j] - fl1_fx * tey_xzz_xyz_1[j];

                tez_xxzz_xxyz_0[j] = pa_x[j] * tez_xzz_xxyz_0[j] - pc_x[j] * tez_xzz_xxyz_1[j] + 0.5 * fl1_fx * tez_zz_xxyz_0[j] -
                                     0.5 * fl1_fx * tez_zz_xxyz_1[j] + fl1_fx * tez_xzz_xyz_0[j] - fl1_fx * tez_xzz_xyz_1[j];

                tex_xxzz_xxzz_0[j] = pa_x[j] * tex_xzz_xxzz_0[j] - pc_x[j] * tex_xzz_xxzz_1[j] + 0.5 * fl1_fx * tex_zz_xxzz_0[j] -
                                     0.5 * fl1_fx * tex_zz_xxzz_1[j] + fl1_fx * tex_xzz_xzz_0[j] - fl1_fx * tex_xzz_xzz_1[j] + ta_xzz_xxzz_1[j];

                tey_xxzz_xxzz_0[j] = pa_x[j] * tey_xzz_xxzz_0[j] - pc_x[j] * tey_xzz_xxzz_1[j] + 0.5 * fl1_fx * tey_zz_xxzz_0[j] -
                                     0.5 * fl1_fx * tey_zz_xxzz_1[j] + fl1_fx * tey_xzz_xzz_0[j] - fl1_fx * tey_xzz_xzz_1[j];

                tez_xxzz_xxzz_0[j] = pa_x[j] * tez_xzz_xxzz_0[j] - pc_x[j] * tez_xzz_xxzz_1[j] + 0.5 * fl1_fx * tez_zz_xxzz_0[j] -
                                     0.5 * fl1_fx * tez_zz_xxzz_1[j] + fl1_fx * tez_xzz_xzz_0[j] - fl1_fx * tez_xzz_xzz_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForGG_243_291(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pcDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto)
{
    // Batch of Integrals (243,291)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_4_4_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_2_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tex_xzz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 81);

            auto tey_xzz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 81);

            auto tez_xzz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 81);

            auto tex_xzz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 82);

            auto tey_xzz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 82);

            auto tez_xzz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 82);

            auto tex_xzz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 83);

            auto tey_xzz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 83);

            auto tez_xzz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 83);

            auto tex_xzz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 84);

            auto tey_xzz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 84);

            auto tez_xzz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 84);

            auto tex_xzz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 85);

            auto tey_xzz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 85);

            auto tez_xzz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 85);

            auto tex_xzz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 86);

            auto tey_xzz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 86);

            auto tez_xzz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 86);

            auto tex_xzz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 87);

            auto tey_xzz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 87);

            auto tez_xzz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 87);

            auto tex_xzz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 88);

            auto tey_xzz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 88);

            auto tez_xzz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 88);

            auto tex_xzz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 89);

            auto tey_xzz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 89);

            auto tez_xzz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 89);

            auto tex_yyy_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 90);

            auto tey_yyy_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 90);

            auto tez_yyy_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 90);

            auto tex_yyy_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 91);

            auto tey_yyy_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 91);

            auto tez_yyy_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 91);

            auto tex_yyy_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 92);

            auto tey_yyy_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 92);

            auto tez_yyy_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 92);

            auto tex_yyy_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 93);

            auto tey_yyy_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 93);

            auto tez_yyy_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 93);

            auto tex_yyy_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 94);

            auto tey_yyy_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 94);

            auto tez_yyy_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 94);

            auto tex_yyy_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 95);

            auto tey_yyy_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 95);

            auto tez_yyy_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 95);

            auto tex_yyy_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 96);

            auto tey_yyy_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 96);

            auto tez_yyy_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 96);

            auto tex_xzz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 81);

            auto tey_xzz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 81);

            auto tez_xzz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 81);

            auto tex_xzz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 82);

            auto tey_xzz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 82);

            auto tez_xzz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 82);

            auto tex_xzz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 83);

            auto tey_xzz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 83);

            auto tez_xzz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 83);

            auto tex_xzz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 84);

            auto tey_xzz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 84);

            auto tez_xzz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 84);

            auto tex_xzz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 85);

            auto tey_xzz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 85);

            auto tez_xzz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 85);

            auto tex_xzz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 86);

            auto tey_xzz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 86);

            auto tez_xzz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 86);

            auto tex_xzz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 87);

            auto tey_xzz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 87);

            auto tez_xzz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 87);

            auto tex_xzz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 88);

            auto tey_xzz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 88);

            auto tez_xzz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 88);

            auto tex_xzz_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 89);

            auto tey_xzz_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 89);

            auto tez_xzz_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 89);

            auto tex_yyy_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 90);

            auto tey_yyy_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 90);

            auto tez_yyy_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 90);

            auto tex_yyy_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 91);

            auto tey_yyy_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 91);

            auto tez_yyy_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 91);

            auto tex_yyy_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 92);

            auto tey_yyy_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 92);

            auto tez_yyy_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 92);

            auto tex_yyy_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 93);

            auto tey_yyy_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 93);

            auto tez_yyy_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 93);

            auto tex_yyy_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 94);

            auto tey_yyy_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 94);

            auto tez_yyy_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 94);

            auto tex_yyy_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 95);

            auto tey_yyy_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 95);

            auto tez_yyy_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 95);

            auto tex_yyy_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 96);

            auto tey_yyy_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 96);

            auto tez_yyy_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 96);

            auto tex_zz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 81);

            auto tey_zz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 81);

            auto tez_zz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 81);

            auto tex_zz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 82);

            auto tey_zz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 82);

            auto tez_zz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 82);

            auto tex_zz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 83);

            auto tey_zz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 83);

            auto tez_zz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 83);

            auto tex_zz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 84);

            auto tey_zz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 84);

            auto tez_zz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 84);

            auto tex_zz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 85);

            auto tey_zz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 85);

            auto tez_zz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 85);

            auto tex_zz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 86);

            auto tey_zz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 86);

            auto tez_zz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 86);

            auto tex_zz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 87);

            auto tey_zz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 87);

            auto tez_zz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 87);

            auto tex_zz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 88);

            auto tey_zz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 88);

            auto tez_zz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 88);

            auto tex_zz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 89);

            auto tey_zz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 89);

            auto tez_zz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 89);

            auto tex_zz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 81);

            auto tey_zz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 81);

            auto tez_zz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 81);

            auto tex_zz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 82);

            auto tey_zz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 82);

            auto tez_zz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 82);

            auto tex_zz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 83);

            auto tey_zz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 83);

            auto tez_zz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 83);

            auto tex_zz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 84);

            auto tey_zz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 84);

            auto tez_zz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 84);

            auto tex_zz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 85);

            auto tey_zz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 85);

            auto tez_zz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 85);

            auto tex_zz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 86);

            auto tey_zz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 86);

            auto tez_zz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 86);

            auto tex_zz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 87);

            auto tey_zz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 87);

            auto tez_zz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 87);

            auto tex_zz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 88);

            auto tey_zz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 88);

            auto tez_zz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 88);

            auto tex_zz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 89);

            auto tey_zz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 89);

            auto tez_zz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 89);

            auto tex_xzz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 56);

            auto tey_xzz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 56);

            auto tez_xzz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 56);

            auto tex_xzz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 57);

            auto tey_xzz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 57);

            auto tez_xzz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 57);

            auto tex_xzz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 58);

            auto tey_xzz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 58);

            auto tez_xzz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 58);

            auto tex_xzz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 59);

            auto tey_xzz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 59);

            auto tez_xzz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 59);

            auto tex_yyy_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 60);

            auto tey_yyy_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 60);

            auto tez_yyy_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 60);

            auto tex_yyy_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 61);

            auto tey_yyy_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 61);

            auto tez_yyy_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 61);

            auto tex_yyy_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 62);

            auto tey_yyy_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 62);

            auto tez_yyy_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 62);

            auto tex_yyy_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 63);

            auto tey_yyy_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 63);

            auto tez_yyy_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 63);

            auto tex_yyy_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 64);

            auto tey_yyy_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 64);

            auto tez_yyy_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 64);

            auto tex_yyy_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 65);

            auto tey_yyy_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 65);

            auto tez_yyy_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 65);

            auto tex_yyy_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 66);

            auto tey_yyy_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 66);

            auto tez_yyy_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 66);

            auto tex_xzz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 56);

            auto tey_xzz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 56);

            auto tez_xzz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 56);

            auto tex_xzz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 57);

            auto tey_xzz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 57);

            auto tez_xzz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 57);

            auto tex_xzz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 58);

            auto tey_xzz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 58);

            auto tez_xzz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 58);

            auto tex_xzz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 59);

            auto tey_xzz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 59);

            auto tez_xzz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 59);

            auto tex_yyy_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 60);

            auto tey_yyy_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 60);

            auto tez_yyy_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 60);

            auto tex_yyy_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 61);

            auto tey_yyy_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 61);

            auto tez_yyy_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 61);

            auto tex_yyy_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 62);

            auto tey_yyy_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 62);

            auto tez_yyy_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 62);

            auto tex_yyy_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 63);

            auto tey_yyy_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 63);

            auto tez_yyy_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 63);

            auto tex_yyy_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 64);

            auto tey_yyy_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 64);

            auto tez_yyy_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 64);

            auto tex_yyy_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 65);

            auto tey_yyy_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 65);

            auto tez_yyy_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 65);

            auto tex_yyy_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 66);

            auto tey_yyy_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 66);

            auto tez_yyy_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 66);

            auto ta_xzz_xyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 81);

            auto ta_xzz_xyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 82);

            auto ta_xzz_xyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 83);

            auto ta_xzz_xzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 84);

            auto ta_xzz_yyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 85);

            auto ta_xzz_yyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 86);

            auto ta_xzz_yyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 87);

            auto ta_xzz_yzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 88);

            auto ta_xzz_zzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 89);

            auto ta_yyy_xxxx_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 90);

            auto ta_yyy_xxxy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 91);

            auto ta_yyy_xxxz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 92);

            auto ta_yyy_xxyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 93);

            auto ta_yyy_xxyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 94);

            auto ta_yyy_xxzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 95);

            auto ta_yyy_xyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 96);

            // set up pointers to integrals

            auto tex_xxzz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 81);

            auto tey_xxzz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 81);

            auto tez_xxzz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 81);

            auto tex_xxzz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 82);

            auto tey_xxzz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 82);

            auto tez_xxzz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 82);

            auto tex_xxzz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 83);

            auto tey_xxzz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 83);

            auto tez_xxzz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 83);

            auto tex_xxzz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 84);

            auto tey_xxzz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 84);

            auto tez_xxzz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 84);

            auto tex_xxzz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 85);

            auto tey_xxzz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 85);

            auto tez_xxzz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 85);

            auto tex_xxzz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 86);

            auto tey_xxzz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 86);

            auto tez_xxzz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 86);

            auto tex_xxzz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 87);

            auto tey_xxzz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 87);

            auto tez_xxzz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 87);

            auto tex_xxzz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 88);

            auto tey_xxzz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 88);

            auto tez_xxzz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 88);

            auto tex_xxzz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 89);

            auto tey_xxzz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 89);

            auto tez_xxzz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 89);

            auto tex_xyyy_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 90);

            auto tey_xyyy_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 90);

            auto tez_xyyy_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 90);

            auto tex_xyyy_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 91);

            auto tey_xyyy_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 91);

            auto tez_xyyy_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 91);

            auto tex_xyyy_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 92);

            auto tey_xyyy_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 92);

            auto tez_xyyy_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 92);

            auto tex_xyyy_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 93);

            auto tey_xyyy_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 93);

            auto tez_xyyy_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 93);

            auto tex_xyyy_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 94);

            auto tey_xyyy_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 94);

            auto tez_xyyy_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 94);

            auto tex_xyyy_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 95);

            auto tey_xyyy_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 95);

            auto tez_xyyy_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 95);

            auto tex_xyyy_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 96);

            auto tey_xyyy_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 96);

            auto tez_xyyy_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 96);

            // Batch of Integrals (243,291)

            #pragma omp simd aligned(fx, pa_x, pc_x, ta_xzz_xyyy_1, ta_xzz_xyyz_1, ta_xzz_xyzz_1, \
                                         ta_xzz_xzzz_1, ta_xzz_yyyy_1, ta_xzz_yyyz_1, ta_xzz_yyzz_1, ta_xzz_yzzz_1, \
                                         ta_xzz_zzzz_1, ta_yyy_xxxx_1, ta_yyy_xxxy_1, ta_yyy_xxxz_1, ta_yyy_xxyy_1, \
                                         ta_yyy_xxyz_1, ta_yyy_xxzz_1, ta_yyy_xyyy_1, tex_xxzz_xyyy_0, tex_xxzz_xyyz_0, \
                                         tex_xxzz_xyzz_0, tex_xxzz_xzzz_0, tex_xxzz_yyyy_0, tex_xxzz_yyyz_0, tex_xxzz_yyzz_0, \
                                         tex_xxzz_yzzz_0, tex_xxzz_zzzz_0, tex_xyyy_xxxx_0, tex_xyyy_xxxy_0, tex_xyyy_xxxz_0, \
                                         tex_xyyy_xxyy_0, tex_xyyy_xxyz_0, tex_xyyy_xxzz_0, tex_xyyy_xyyy_0, tex_xzz_xyyy_0, \
                                         tex_xzz_xyyy_1, tex_xzz_xyyz_0, tex_xzz_xyyz_1, tex_xzz_xyzz_0, tex_xzz_xyzz_1, \
                                         tex_xzz_xzzz_0, tex_xzz_xzzz_1, tex_xzz_yyy_0, tex_xzz_yyy_1, tex_xzz_yyyy_0, \
                                         tex_xzz_yyyy_1, tex_xzz_yyyz_0, tex_xzz_yyyz_1, tex_xzz_yyz_0, tex_xzz_yyz_1, \
                                         tex_xzz_yyzz_0, tex_xzz_yyzz_1, tex_xzz_yzz_0, tex_xzz_yzz_1, tex_xzz_yzzz_0, \
                                         tex_xzz_yzzz_1, tex_xzz_zzz_0, tex_xzz_zzz_1, tex_xzz_zzzz_0, tex_xzz_zzzz_1, \
                                         tex_yyy_xxx_0, tex_yyy_xxx_1, tex_yyy_xxxx_0, tex_yyy_xxxx_1, tex_yyy_xxxy_0, \
                                         tex_yyy_xxxy_1, tex_yyy_xxxz_0, tex_yyy_xxxz_1, tex_yyy_xxy_0, tex_yyy_xxy_1, \
                                         tex_yyy_xxyy_0, tex_yyy_xxyy_1, tex_yyy_xxyz_0, tex_yyy_xxyz_1, tex_yyy_xxz_0, \
                                         tex_yyy_xxz_1, tex_yyy_xxzz_0, tex_yyy_xxzz_1, tex_yyy_xyy_0, tex_yyy_xyy_1, \
                                         tex_yyy_xyyy_0, tex_yyy_xyyy_1, tex_yyy_xyz_0, tex_yyy_xyz_1, tex_yyy_xzz_0, \
                                         tex_yyy_xzz_1, tex_yyy_yyy_0, tex_yyy_yyy_1, tex_zz_xyyy_0, tex_zz_xyyy_1, \
                                         tex_zz_xyyz_0, tex_zz_xyyz_1, tex_zz_xyzz_0, tex_zz_xyzz_1, tex_zz_xzzz_0, \
                                         tex_zz_xzzz_1, tex_zz_yyyy_0, tex_zz_yyyy_1, tex_zz_yyyz_0, tex_zz_yyyz_1, \
                                         tex_zz_yyzz_0, tex_zz_yyzz_1, tex_zz_yzzz_0, tex_zz_yzzz_1, tex_zz_zzzz_0, \
                                         tex_zz_zzzz_1, tey_xxzz_xyyy_0, tey_xxzz_xyyz_0, tey_xxzz_xyzz_0, tey_xxzz_xzzz_0, \
                                         tey_xxzz_yyyy_0, tey_xxzz_yyyz_0, tey_xxzz_yyzz_0, tey_xxzz_yzzz_0, tey_xxzz_zzzz_0, \
                                         tey_xyyy_xxxx_0, tey_xyyy_xxxy_0, tey_xyyy_xxxz_0, tey_xyyy_xxyy_0, tey_xyyy_xxyz_0, \
                                         tey_xyyy_xxzz_0, tey_xyyy_xyyy_0, tey_xzz_xyyy_0, tey_xzz_xyyy_1, tey_xzz_xyyz_0, \
                                         tey_xzz_xyyz_1, tey_xzz_xyzz_0, tey_xzz_xyzz_1, tey_xzz_xzzz_0, tey_xzz_xzzz_1, \
                                         tey_xzz_yyy_0, tey_xzz_yyy_1, tey_xzz_yyyy_0, tey_xzz_yyyy_1, tey_xzz_yyyz_0, \
                                         tey_xzz_yyyz_1, tey_xzz_yyz_0, tey_xzz_yyz_1, tey_xzz_yyzz_0, tey_xzz_yyzz_1, \
                                         tey_xzz_yzz_0, tey_xzz_yzz_1, tey_xzz_yzzz_0, tey_xzz_yzzz_1, tey_xzz_zzz_0, \
                                         tey_xzz_zzz_1, tey_xzz_zzzz_0, tey_xzz_zzzz_1, tey_yyy_xxx_0, tey_yyy_xxx_1, \
                                         tey_yyy_xxxx_0, tey_yyy_xxxx_1, tey_yyy_xxxy_0, tey_yyy_xxxy_1, tey_yyy_xxxz_0, \
                                         tey_yyy_xxxz_1, tey_yyy_xxy_0, tey_yyy_xxy_1, tey_yyy_xxyy_0, tey_yyy_xxyy_1, \
                                         tey_yyy_xxyz_0, tey_yyy_xxyz_1, tey_yyy_xxz_0, tey_yyy_xxz_1, tey_yyy_xxzz_0, \
                                         tey_yyy_xxzz_1, tey_yyy_xyy_0, tey_yyy_xyy_1, tey_yyy_xyyy_0, tey_yyy_xyyy_1, \
                                         tey_yyy_xyz_0, tey_yyy_xyz_1, tey_yyy_xzz_0, tey_yyy_xzz_1, tey_yyy_yyy_0, \
                                         tey_yyy_yyy_1, tey_zz_xyyy_0, tey_zz_xyyy_1, tey_zz_xyyz_0, tey_zz_xyyz_1, \
                                         tey_zz_xyzz_0, tey_zz_xyzz_1, tey_zz_xzzz_0, tey_zz_xzzz_1, tey_zz_yyyy_0, \
                                         tey_zz_yyyy_1, tey_zz_yyyz_0, tey_zz_yyyz_1, tey_zz_yyzz_0, tey_zz_yyzz_1, \
                                         tey_zz_yzzz_0, tey_zz_yzzz_1, tey_zz_zzzz_0, tey_zz_zzzz_1, tez_xxzz_xyyy_0, \
                                         tez_xxzz_xyyz_0, tez_xxzz_xyzz_0, tez_xxzz_xzzz_0, tez_xxzz_yyyy_0, tez_xxzz_yyyz_0, \
                                         tez_xxzz_yyzz_0, tez_xxzz_yzzz_0, tez_xxzz_zzzz_0, tez_xyyy_xxxx_0, tez_xyyy_xxxy_0, \
                                         tez_xyyy_xxxz_0, tez_xyyy_xxyy_0, tez_xyyy_xxyz_0, tez_xyyy_xxzz_0, tez_xyyy_xyyy_0, \
                                         tez_xzz_xyyy_0, tez_xzz_xyyy_1, tez_xzz_xyyz_0, tez_xzz_xyyz_1, tez_xzz_xyzz_0, \
                                         tez_xzz_xyzz_1, tez_xzz_xzzz_0, tez_xzz_xzzz_1, tez_xzz_yyy_0, tez_xzz_yyy_1, \
                                         tez_xzz_yyyy_0, tez_xzz_yyyy_1, tez_xzz_yyyz_0, tez_xzz_yyyz_1, tez_xzz_yyz_0, \
                                         tez_xzz_yyz_1, tez_xzz_yyzz_0, tez_xzz_yyzz_1, tez_xzz_yzz_0, tez_xzz_yzz_1, \
                                         tez_xzz_yzzz_0, tez_xzz_yzzz_1, tez_xzz_zzz_0, tez_xzz_zzz_1, tez_xzz_zzzz_0, \
                                         tez_xzz_zzzz_1, tez_yyy_xxx_0, tez_yyy_xxx_1, tez_yyy_xxxx_0, tez_yyy_xxxx_1, \
                                         tez_yyy_xxxy_0, tez_yyy_xxxy_1, tez_yyy_xxxz_0, tez_yyy_xxxz_1, tez_yyy_xxy_0, \
                                         tez_yyy_xxy_1, tez_yyy_xxyy_0, tez_yyy_xxyy_1, tez_yyy_xxyz_0, tez_yyy_xxyz_1, \
                                         tez_yyy_xxz_0, tez_yyy_xxz_1, tez_yyy_xxzz_0, tez_yyy_xxzz_1, tez_yyy_xyy_0, \
                                         tez_yyy_xyy_1, tez_yyy_xyyy_0, tez_yyy_xyyy_1, tez_yyy_xyz_0, tez_yyy_xyz_1, \
                                         tez_yyy_xzz_0, tez_yyy_xzz_1, tez_yyy_yyy_0, tez_yyy_yyy_1, tez_zz_xyyy_0, \
                                         tez_zz_xyyy_1, tez_zz_xyyz_0, tez_zz_xyyz_1, tez_zz_xyzz_0, tez_zz_xyzz_1, \
                                         tez_zz_xzzz_0, tez_zz_xzzz_1, tez_zz_yyyy_0, tez_zz_yyyy_1, tez_zz_yyyz_0, \
                                         tez_zz_yyyz_1, tez_zz_yyzz_0, tez_zz_yyzz_1, tez_zz_yzzz_0, tez_zz_yzzz_1, \
                                         tez_zz_zzzz_0, tez_zz_zzzz_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_xxzz_xyyy_0[j] = pa_x[j] * tex_xzz_xyyy_0[j] - pc_x[j] * tex_xzz_xyyy_1[j] + 0.5 * fl1_fx * tex_zz_xyyy_0[j] -
                                     0.5 * fl1_fx * tex_zz_xyyy_1[j] + 0.5 * fl1_fx * tex_xzz_yyy_0[j] - 0.5 * fl1_fx * tex_xzz_yyy_1[j] +
                                     ta_xzz_xyyy_1[j];

                tey_xxzz_xyyy_0[j] = pa_x[j] * tey_xzz_xyyy_0[j] - pc_x[j] * tey_xzz_xyyy_1[j] + 0.5 * fl1_fx * tey_zz_xyyy_0[j] -
                                     0.5 * fl1_fx * tey_zz_xyyy_1[j] + 0.5 * fl1_fx * tey_xzz_yyy_0[j] - 0.5 * fl1_fx * tey_xzz_yyy_1[j];

                tez_xxzz_xyyy_0[j] = pa_x[j] * tez_xzz_xyyy_0[j] - pc_x[j] * tez_xzz_xyyy_1[j] + 0.5 * fl1_fx * tez_zz_xyyy_0[j] -
                                     0.5 * fl1_fx * tez_zz_xyyy_1[j] + 0.5 * fl1_fx * tez_xzz_yyy_0[j] - 0.5 * fl1_fx * tez_xzz_yyy_1[j];

                tex_xxzz_xyyz_0[j] = pa_x[j] * tex_xzz_xyyz_0[j] - pc_x[j] * tex_xzz_xyyz_1[j] + 0.5 * fl1_fx * tex_zz_xyyz_0[j] -
                                     0.5 * fl1_fx * tex_zz_xyyz_1[j] + 0.5 * fl1_fx * tex_xzz_yyz_0[j] - 0.5 * fl1_fx * tex_xzz_yyz_1[j] +
                                     ta_xzz_xyyz_1[j];

                tey_xxzz_xyyz_0[j] = pa_x[j] * tey_xzz_xyyz_0[j] - pc_x[j] * tey_xzz_xyyz_1[j] + 0.5 * fl1_fx * tey_zz_xyyz_0[j] -
                                     0.5 * fl1_fx * tey_zz_xyyz_1[j] + 0.5 * fl1_fx * tey_xzz_yyz_0[j] - 0.5 * fl1_fx * tey_xzz_yyz_1[j];

                tez_xxzz_xyyz_0[j] = pa_x[j] * tez_xzz_xyyz_0[j] - pc_x[j] * tez_xzz_xyyz_1[j] + 0.5 * fl1_fx * tez_zz_xyyz_0[j] -
                                     0.5 * fl1_fx * tez_zz_xyyz_1[j] + 0.5 * fl1_fx * tez_xzz_yyz_0[j] - 0.5 * fl1_fx * tez_xzz_yyz_1[j];

                tex_xxzz_xyzz_0[j] = pa_x[j] * tex_xzz_xyzz_0[j] - pc_x[j] * tex_xzz_xyzz_1[j] + 0.5 * fl1_fx * tex_zz_xyzz_0[j] -
                                     0.5 * fl1_fx * tex_zz_xyzz_1[j] + 0.5 * fl1_fx * tex_xzz_yzz_0[j] - 0.5 * fl1_fx * tex_xzz_yzz_1[j] +
                                     ta_xzz_xyzz_1[j];

                tey_xxzz_xyzz_0[j] = pa_x[j] * tey_xzz_xyzz_0[j] - pc_x[j] * tey_xzz_xyzz_1[j] + 0.5 * fl1_fx * tey_zz_xyzz_0[j] -
                                     0.5 * fl1_fx * tey_zz_xyzz_1[j] + 0.5 * fl1_fx * tey_xzz_yzz_0[j] - 0.5 * fl1_fx * tey_xzz_yzz_1[j];

                tez_xxzz_xyzz_0[j] = pa_x[j] * tez_xzz_xyzz_0[j] - pc_x[j] * tez_xzz_xyzz_1[j] + 0.5 * fl1_fx * tez_zz_xyzz_0[j] -
                                     0.5 * fl1_fx * tez_zz_xyzz_1[j] + 0.5 * fl1_fx * tez_xzz_yzz_0[j] - 0.5 * fl1_fx * tez_xzz_yzz_1[j];

                tex_xxzz_xzzz_0[j] = pa_x[j] * tex_xzz_xzzz_0[j] - pc_x[j] * tex_xzz_xzzz_1[j] + 0.5 * fl1_fx * tex_zz_xzzz_0[j] -
                                     0.5 * fl1_fx * tex_zz_xzzz_1[j] + 0.5 * fl1_fx * tex_xzz_zzz_0[j] - 0.5 * fl1_fx * tex_xzz_zzz_1[j] +
                                     ta_xzz_xzzz_1[j];

                tey_xxzz_xzzz_0[j] = pa_x[j] * tey_xzz_xzzz_0[j] - pc_x[j] * tey_xzz_xzzz_1[j] + 0.5 * fl1_fx * tey_zz_xzzz_0[j] -
                                     0.5 * fl1_fx * tey_zz_xzzz_1[j] + 0.5 * fl1_fx * tey_xzz_zzz_0[j] - 0.5 * fl1_fx * tey_xzz_zzz_1[j];

                tez_xxzz_xzzz_0[j] = pa_x[j] * tez_xzz_xzzz_0[j] - pc_x[j] * tez_xzz_xzzz_1[j] + 0.5 * fl1_fx * tez_zz_xzzz_0[j] -
                                     0.5 * fl1_fx * tez_zz_xzzz_1[j] + 0.5 * fl1_fx * tez_xzz_zzz_0[j] - 0.5 * fl1_fx * tez_xzz_zzz_1[j];

                tex_xxzz_yyyy_0[j] = pa_x[j] * tex_xzz_yyyy_0[j] - pc_x[j] * tex_xzz_yyyy_1[j] + 0.5 * fl1_fx * tex_zz_yyyy_0[j] -
                                     0.5 * fl1_fx * tex_zz_yyyy_1[j] + ta_xzz_yyyy_1[j];

                tey_xxzz_yyyy_0[j] =
                    pa_x[j] * tey_xzz_yyyy_0[j] - pc_x[j] * tey_xzz_yyyy_1[j] + 0.5 * fl1_fx * tey_zz_yyyy_0[j] - 0.5 * fl1_fx * tey_zz_yyyy_1[j];

                tez_xxzz_yyyy_0[j] =
                    pa_x[j] * tez_xzz_yyyy_0[j] - pc_x[j] * tez_xzz_yyyy_1[j] + 0.5 * fl1_fx * tez_zz_yyyy_0[j] - 0.5 * fl1_fx * tez_zz_yyyy_1[j];

                tex_xxzz_yyyz_0[j] = pa_x[j] * tex_xzz_yyyz_0[j] - pc_x[j] * tex_xzz_yyyz_1[j] + 0.5 * fl1_fx * tex_zz_yyyz_0[j] -
                                     0.5 * fl1_fx * tex_zz_yyyz_1[j] + ta_xzz_yyyz_1[j];

                tey_xxzz_yyyz_0[j] =
                    pa_x[j] * tey_xzz_yyyz_0[j] - pc_x[j] * tey_xzz_yyyz_1[j] + 0.5 * fl1_fx * tey_zz_yyyz_0[j] - 0.5 * fl1_fx * tey_zz_yyyz_1[j];

                tez_xxzz_yyyz_0[j] =
                    pa_x[j] * tez_xzz_yyyz_0[j] - pc_x[j] * tez_xzz_yyyz_1[j] + 0.5 * fl1_fx * tez_zz_yyyz_0[j] - 0.5 * fl1_fx * tez_zz_yyyz_1[j];

                tex_xxzz_yyzz_0[j] = pa_x[j] * tex_xzz_yyzz_0[j] - pc_x[j] * tex_xzz_yyzz_1[j] + 0.5 * fl1_fx * tex_zz_yyzz_0[j] -
                                     0.5 * fl1_fx * tex_zz_yyzz_1[j] + ta_xzz_yyzz_1[j];

                tey_xxzz_yyzz_0[j] =
                    pa_x[j] * tey_xzz_yyzz_0[j] - pc_x[j] * tey_xzz_yyzz_1[j] + 0.5 * fl1_fx * tey_zz_yyzz_0[j] - 0.5 * fl1_fx * tey_zz_yyzz_1[j];

                tez_xxzz_yyzz_0[j] =
                    pa_x[j] * tez_xzz_yyzz_0[j] - pc_x[j] * tez_xzz_yyzz_1[j] + 0.5 * fl1_fx * tez_zz_yyzz_0[j] - 0.5 * fl1_fx * tez_zz_yyzz_1[j];

                tex_xxzz_yzzz_0[j] = pa_x[j] * tex_xzz_yzzz_0[j] - pc_x[j] * tex_xzz_yzzz_1[j] + 0.5 * fl1_fx * tex_zz_yzzz_0[j] -
                                     0.5 * fl1_fx * tex_zz_yzzz_1[j] + ta_xzz_yzzz_1[j];

                tey_xxzz_yzzz_0[j] =
                    pa_x[j] * tey_xzz_yzzz_0[j] - pc_x[j] * tey_xzz_yzzz_1[j] + 0.5 * fl1_fx * tey_zz_yzzz_0[j] - 0.5 * fl1_fx * tey_zz_yzzz_1[j];

                tez_xxzz_yzzz_0[j] =
                    pa_x[j] * tez_xzz_yzzz_0[j] - pc_x[j] * tez_xzz_yzzz_1[j] + 0.5 * fl1_fx * tez_zz_yzzz_0[j] - 0.5 * fl1_fx * tez_zz_yzzz_1[j];

                tex_xxzz_zzzz_0[j] = pa_x[j] * tex_xzz_zzzz_0[j] - pc_x[j] * tex_xzz_zzzz_1[j] + 0.5 * fl1_fx * tex_zz_zzzz_0[j] -
                                     0.5 * fl1_fx * tex_zz_zzzz_1[j] + ta_xzz_zzzz_1[j];

                tey_xxzz_zzzz_0[j] =
                    pa_x[j] * tey_xzz_zzzz_0[j] - pc_x[j] * tey_xzz_zzzz_1[j] + 0.5 * fl1_fx * tey_zz_zzzz_0[j] - 0.5 * fl1_fx * tey_zz_zzzz_1[j];

                tez_xxzz_zzzz_0[j] =
                    pa_x[j] * tez_xzz_zzzz_0[j] - pc_x[j] * tez_xzz_zzzz_1[j] + 0.5 * fl1_fx * tez_zz_zzzz_0[j] - 0.5 * fl1_fx * tez_zz_zzzz_1[j];

                tex_xyyy_xxxx_0[j] = pa_x[j] * tex_yyy_xxxx_0[j] - pc_x[j] * tex_yyy_xxxx_1[j] + 2.0 * fl1_fx * tex_yyy_xxx_0[j] -
                                     2.0 * fl1_fx * tex_yyy_xxx_1[j] + ta_yyy_xxxx_1[j];

                tey_xyyy_xxxx_0[j] =
                    pa_x[j] * tey_yyy_xxxx_0[j] - pc_x[j] * tey_yyy_xxxx_1[j] + 2.0 * fl1_fx * tey_yyy_xxx_0[j] - 2.0 * fl1_fx * tey_yyy_xxx_1[j];

                tez_xyyy_xxxx_0[j] =
                    pa_x[j] * tez_yyy_xxxx_0[j] - pc_x[j] * tez_yyy_xxxx_1[j] + 2.0 * fl1_fx * tez_yyy_xxx_0[j] - 2.0 * fl1_fx * tez_yyy_xxx_1[j];

                tex_xyyy_xxxy_0[j] = pa_x[j] * tex_yyy_xxxy_0[j] - pc_x[j] * tex_yyy_xxxy_1[j] + 1.5 * fl1_fx * tex_yyy_xxy_0[j] -
                                     1.5 * fl1_fx * tex_yyy_xxy_1[j] + ta_yyy_xxxy_1[j];

                tey_xyyy_xxxy_0[j] =
                    pa_x[j] * tey_yyy_xxxy_0[j] - pc_x[j] * tey_yyy_xxxy_1[j] + 1.5 * fl1_fx * tey_yyy_xxy_0[j] - 1.5 * fl1_fx * tey_yyy_xxy_1[j];

                tez_xyyy_xxxy_0[j] =
                    pa_x[j] * tez_yyy_xxxy_0[j] - pc_x[j] * tez_yyy_xxxy_1[j] + 1.5 * fl1_fx * tez_yyy_xxy_0[j] - 1.5 * fl1_fx * tez_yyy_xxy_1[j];

                tex_xyyy_xxxz_0[j] = pa_x[j] * tex_yyy_xxxz_0[j] - pc_x[j] * tex_yyy_xxxz_1[j] + 1.5 * fl1_fx * tex_yyy_xxz_0[j] -
                                     1.5 * fl1_fx * tex_yyy_xxz_1[j] + ta_yyy_xxxz_1[j];

                tey_xyyy_xxxz_0[j] =
                    pa_x[j] * tey_yyy_xxxz_0[j] - pc_x[j] * tey_yyy_xxxz_1[j] + 1.5 * fl1_fx * tey_yyy_xxz_0[j] - 1.5 * fl1_fx * tey_yyy_xxz_1[j];

                tez_xyyy_xxxz_0[j] =
                    pa_x[j] * tez_yyy_xxxz_0[j] - pc_x[j] * tez_yyy_xxxz_1[j] + 1.5 * fl1_fx * tez_yyy_xxz_0[j] - 1.5 * fl1_fx * tez_yyy_xxz_1[j];

                tex_xyyy_xxyy_0[j] = pa_x[j] * tex_yyy_xxyy_0[j] - pc_x[j] * tex_yyy_xxyy_1[j] + fl1_fx * tex_yyy_xyy_0[j] -
                                     fl1_fx * tex_yyy_xyy_1[j] + ta_yyy_xxyy_1[j];

                tey_xyyy_xxyy_0[j] =
                    pa_x[j] * tey_yyy_xxyy_0[j] - pc_x[j] * tey_yyy_xxyy_1[j] + fl1_fx * tey_yyy_xyy_0[j] - fl1_fx * tey_yyy_xyy_1[j];

                tez_xyyy_xxyy_0[j] =
                    pa_x[j] * tez_yyy_xxyy_0[j] - pc_x[j] * tez_yyy_xxyy_1[j] + fl1_fx * tez_yyy_xyy_0[j] - fl1_fx * tez_yyy_xyy_1[j];

                tex_xyyy_xxyz_0[j] = pa_x[j] * tex_yyy_xxyz_0[j] - pc_x[j] * tex_yyy_xxyz_1[j] + fl1_fx * tex_yyy_xyz_0[j] -
                                     fl1_fx * tex_yyy_xyz_1[j] + ta_yyy_xxyz_1[j];

                tey_xyyy_xxyz_0[j] =
                    pa_x[j] * tey_yyy_xxyz_0[j] - pc_x[j] * tey_yyy_xxyz_1[j] + fl1_fx * tey_yyy_xyz_0[j] - fl1_fx * tey_yyy_xyz_1[j];

                tez_xyyy_xxyz_0[j] =
                    pa_x[j] * tez_yyy_xxyz_0[j] - pc_x[j] * tez_yyy_xxyz_1[j] + fl1_fx * tez_yyy_xyz_0[j] - fl1_fx * tez_yyy_xyz_1[j];

                tex_xyyy_xxzz_0[j] = pa_x[j] * tex_yyy_xxzz_0[j] - pc_x[j] * tex_yyy_xxzz_1[j] + fl1_fx * tex_yyy_xzz_0[j] -
                                     fl1_fx * tex_yyy_xzz_1[j] + ta_yyy_xxzz_1[j];

                tey_xyyy_xxzz_0[j] =
                    pa_x[j] * tey_yyy_xxzz_0[j] - pc_x[j] * tey_yyy_xxzz_1[j] + fl1_fx * tey_yyy_xzz_0[j] - fl1_fx * tey_yyy_xzz_1[j];

                tez_xyyy_xxzz_0[j] =
                    pa_x[j] * tez_yyy_xxzz_0[j] - pc_x[j] * tez_yyy_xxzz_1[j] + fl1_fx * tez_yyy_xzz_0[j] - fl1_fx * tez_yyy_xzz_1[j];

                tex_xyyy_xyyy_0[j] = pa_x[j] * tex_yyy_xyyy_0[j] - pc_x[j] * tex_yyy_xyyy_1[j] + 0.5 * fl1_fx * tex_yyy_yyy_0[j] -
                                     0.5 * fl1_fx * tex_yyy_yyy_1[j] + ta_yyy_xyyy_1[j];

                tey_xyyy_xyyy_0[j] =
                    pa_x[j] * tey_yyy_xyyy_0[j] - pc_x[j] * tey_yyy_xyyy_1[j] + 0.5 * fl1_fx * tey_yyy_yyy_0[j] - 0.5 * fl1_fx * tey_yyy_yyy_1[j];

                tez_xyyy_xyyy_0[j] =
                    pa_x[j] * tez_yyy_xyyy_0[j] - pc_x[j] * tez_yyy_xyyy_1[j] + 0.5 * fl1_fx * tez_yyy_yyy_0[j] - 0.5 * fl1_fx * tez_yyy_yyy_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForGG_291_339(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pcDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto)
{
    // Batch of Integrals (291,339)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_4_4_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tex_yyy_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 97);

            auto tey_yyy_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 97);

            auto tez_yyy_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 97);

            auto tex_yyy_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 98);

            auto tey_yyy_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 98);

            auto tez_yyy_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 98);

            auto tex_yyy_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 99);

            auto tey_yyy_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 99);

            auto tez_yyy_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 99);

            auto tex_yyy_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 100);

            auto tey_yyy_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 100);

            auto tez_yyy_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 100);

            auto tex_yyy_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 101);

            auto tey_yyy_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 101);

            auto tez_yyy_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 101);

            auto tex_yyy_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 102);

            auto tey_yyy_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 102);

            auto tez_yyy_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 102);

            auto tex_yyy_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 103);

            auto tey_yyy_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 103);

            auto tez_yyy_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 103);

            auto tex_yyy_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 104);

            auto tey_yyy_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 104);

            auto tez_yyy_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 104);

            auto tex_yyz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 105);

            auto tey_yyz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 105);

            auto tez_yyz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 105);

            auto tex_yyz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 106);

            auto tey_yyz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 106);

            auto tez_yyz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 106);

            auto tex_yyz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 107);

            auto tey_yyz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 107);

            auto tez_yyz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 107);

            auto tex_yyz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 108);

            auto tey_yyz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 108);

            auto tez_yyz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 108);

            auto tex_yyz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 109);

            auto tey_yyz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 109);

            auto tez_yyz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 109);

            auto tex_yyz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 110);

            auto tey_yyz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 110);

            auto tez_yyz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 110);

            auto tex_yyz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 111);

            auto tey_yyz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 111);

            auto tez_yyz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 111);

            auto tex_yyz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 112);

            auto tey_yyz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 112);

            auto tez_yyz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 112);

            auto tex_yyy_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 97);

            auto tey_yyy_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 97);

            auto tez_yyy_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 97);

            auto tex_yyy_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 98);

            auto tey_yyy_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 98);

            auto tez_yyy_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 98);

            auto tex_yyy_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 99);

            auto tey_yyy_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 99);

            auto tez_yyy_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 99);

            auto tex_yyy_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 100);

            auto tey_yyy_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 100);

            auto tez_yyy_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 100);

            auto tex_yyy_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 101);

            auto tey_yyy_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 101);

            auto tez_yyy_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 101);

            auto tex_yyy_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 102);

            auto tey_yyy_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 102);

            auto tez_yyy_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 102);

            auto tex_yyy_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 103);

            auto tey_yyy_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 103);

            auto tez_yyy_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 103);

            auto tex_yyy_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 104);

            auto tey_yyy_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 104);

            auto tez_yyy_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 104);

            auto tex_yyz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 105);

            auto tey_yyz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 105);

            auto tez_yyz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 105);

            auto tex_yyz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 106);

            auto tey_yyz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 106);

            auto tez_yyz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 106);

            auto tex_yyz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 107);

            auto tey_yyz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 107);

            auto tez_yyz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 107);

            auto tex_yyz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 108);

            auto tey_yyz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 108);

            auto tez_yyz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 108);

            auto tex_yyz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 109);

            auto tey_yyz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 109);

            auto tez_yyz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 109);

            auto tex_yyz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 110);

            auto tey_yyz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 110);

            auto tez_yyz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 110);

            auto tex_yyz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 111);

            auto tey_yyz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 111);

            auto tez_yyz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 111);

            auto tex_yyz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 112);

            auto tey_yyz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 112);

            auto tez_yyz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 112);

            auto tex_yyy_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 67);

            auto tey_yyy_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 67);

            auto tez_yyy_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 67);

            auto tex_yyy_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 68);

            auto tey_yyy_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 68);

            auto tez_yyy_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 68);

            auto tex_yyy_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 69);

            auto tey_yyy_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 69);

            auto tez_yyy_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 69);

            auto tex_yyz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 70);

            auto tey_yyz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 70);

            auto tez_yyz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 70);

            auto tex_yyz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 71);

            auto tey_yyz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 71);

            auto tez_yyz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 71);

            auto tex_yyz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 72);

            auto tey_yyz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 72);

            auto tez_yyz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 72);

            auto tex_yyz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 73);

            auto tey_yyz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 73);

            auto tez_yyz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 73);

            auto tex_yyz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 74);

            auto tey_yyz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 74);

            auto tez_yyz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 74);

            auto tex_yyz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 75);

            auto tey_yyz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 75);

            auto tez_yyz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 75);

            auto tex_yyz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 76);

            auto tey_yyz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 76);

            auto tez_yyz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 76);

            auto tex_yyz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 77);

            auto tey_yyz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 77);

            auto tez_yyz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 77);

            auto tex_yyy_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 67);

            auto tey_yyy_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 67);

            auto tez_yyy_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 67);

            auto tex_yyy_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 68);

            auto tey_yyy_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 68);

            auto tez_yyy_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 68);

            auto tex_yyy_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 69);

            auto tey_yyy_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 69);

            auto tez_yyy_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 69);

            auto tex_yyz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 70);

            auto tey_yyz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 70);

            auto tez_yyz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 70);

            auto tex_yyz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 71);

            auto tey_yyz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 71);

            auto tez_yyz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 71);

            auto tex_yyz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 72);

            auto tey_yyz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 72);

            auto tez_yyz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 72);

            auto tex_yyz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 73);

            auto tey_yyz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 73);

            auto tez_yyz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 73);

            auto tex_yyz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 74);

            auto tey_yyz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 74);

            auto tez_yyz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 74);

            auto tex_yyz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 75);

            auto tey_yyz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 75);

            auto tez_yyz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 75);

            auto tex_yyz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 76);

            auto tey_yyz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 76);

            auto tez_yyz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 76);

            auto tex_yyz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 77);

            auto tey_yyz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 77);

            auto tez_yyz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 77);

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

            // set up pointers to integrals

            auto tex_xyyy_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 97);

            auto tey_xyyy_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 97);

            auto tez_xyyy_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 97);

            auto tex_xyyy_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 98);

            auto tey_xyyy_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 98);

            auto tez_xyyy_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 98);

            auto tex_xyyy_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 99);

            auto tey_xyyy_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 99);

            auto tez_xyyy_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 99);

            auto tex_xyyy_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 100);

            auto tey_xyyy_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 100);

            auto tez_xyyy_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 100);

            auto tex_xyyy_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 101);

            auto tey_xyyy_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 101);

            auto tez_xyyy_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 101);

            auto tex_xyyy_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 102);

            auto tey_xyyy_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 102);

            auto tez_xyyy_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 102);

            auto tex_xyyy_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 103);

            auto tey_xyyy_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 103);

            auto tez_xyyy_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 103);

            auto tex_xyyy_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 104);

            auto tey_xyyy_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 104);

            auto tez_xyyy_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 104);

            auto tex_xyyz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 105);

            auto tey_xyyz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 105);

            auto tez_xyyz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 105);

            auto tex_xyyz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 106);

            auto tey_xyyz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 106);

            auto tez_xyyz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 106);

            auto tex_xyyz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 107);

            auto tey_xyyz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 107);

            auto tez_xyyz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 107);

            auto tex_xyyz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 108);

            auto tey_xyyz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 108);

            auto tez_xyyz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 108);

            auto tex_xyyz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 109);

            auto tey_xyyz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 109);

            auto tez_xyyz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 109);

            auto tex_xyyz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 110);

            auto tey_xyyz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 110);

            auto tez_xyyz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 110);

            auto tex_xyyz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 111);

            auto tey_xyyz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 111);

            auto tez_xyyz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 111);

            auto tex_xyyz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 112);

            auto tey_xyyz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 112);

            auto tez_xyyz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 112);

            // Batch of Integrals (291,339)

            #pragma omp simd aligned(fx, pa_x, pc_x, ta_yyy_xyyz_1, ta_yyy_xyzz_1, ta_yyy_xzzz_1, \
                                         ta_yyy_yyyy_1, ta_yyy_yyyz_1, ta_yyy_yyzz_1, ta_yyy_yzzz_1, ta_yyy_zzzz_1, \
                                         ta_yyz_xxxx_1, ta_yyz_xxxy_1, ta_yyz_xxxz_1, ta_yyz_xxyy_1, ta_yyz_xxyz_1, \
                                         ta_yyz_xxzz_1, ta_yyz_xyyy_1, ta_yyz_xyyz_1, tex_xyyy_xyyz_0, tex_xyyy_xyzz_0, \
                                         tex_xyyy_xzzz_0, tex_xyyy_yyyy_0, tex_xyyy_yyyz_0, tex_xyyy_yyzz_0, tex_xyyy_yzzz_0, \
                                         tex_xyyy_zzzz_0, tex_xyyz_xxxx_0, tex_xyyz_xxxy_0, tex_xyyz_xxxz_0, tex_xyyz_xxyy_0, \
                                         tex_xyyz_xxyz_0, tex_xyyz_xxzz_0, tex_xyyz_xyyy_0, tex_xyyz_xyyz_0, tex_yyy_xyyz_0, \
                                         tex_yyy_xyyz_1, tex_yyy_xyzz_0, tex_yyy_xyzz_1, tex_yyy_xzzz_0, tex_yyy_xzzz_1, \
                                         tex_yyy_yyyy_0, tex_yyy_yyyy_1, tex_yyy_yyyz_0, tex_yyy_yyyz_1, tex_yyy_yyz_0, \
                                         tex_yyy_yyz_1, tex_yyy_yyzz_0, tex_yyy_yyzz_1, tex_yyy_yzz_0, tex_yyy_yzz_1, \
                                         tex_yyy_yzzz_0, tex_yyy_yzzz_1, tex_yyy_zzz_0, tex_yyy_zzz_1, tex_yyy_zzzz_0, \
                                         tex_yyy_zzzz_1, tex_yyz_xxx_0, tex_yyz_xxx_1, tex_yyz_xxxx_0, tex_yyz_xxxx_1, \
                                         tex_yyz_xxxy_0, tex_yyz_xxxy_1, tex_yyz_xxxz_0, tex_yyz_xxxz_1, tex_yyz_xxy_0, \
                                         tex_yyz_xxy_1, tex_yyz_xxyy_0, tex_yyz_xxyy_1, tex_yyz_xxyz_0, tex_yyz_xxyz_1, \
                                         tex_yyz_xxz_0, tex_yyz_xxz_1, tex_yyz_xxzz_0, tex_yyz_xxzz_1, tex_yyz_xyy_0, \
                                         tex_yyz_xyy_1, tex_yyz_xyyy_0, tex_yyz_xyyy_1, tex_yyz_xyyz_0, tex_yyz_xyyz_1, \
                                         tex_yyz_xyz_0, tex_yyz_xyz_1, tex_yyz_xzz_0, tex_yyz_xzz_1, tex_yyz_yyy_0, \
                                         tex_yyz_yyy_1, tex_yyz_yyz_0, tex_yyz_yyz_1, tey_xyyy_xyyz_0, tey_xyyy_xyzz_0, \
                                         tey_xyyy_xzzz_0, tey_xyyy_yyyy_0, tey_xyyy_yyyz_0, tey_xyyy_yyzz_0, tey_xyyy_yzzz_0, \
                                         tey_xyyy_zzzz_0, tey_xyyz_xxxx_0, tey_xyyz_xxxy_0, tey_xyyz_xxxz_0, tey_xyyz_xxyy_0, \
                                         tey_xyyz_xxyz_0, tey_xyyz_xxzz_0, tey_xyyz_xyyy_0, tey_xyyz_xyyz_0, tey_yyy_xyyz_0, \
                                         tey_yyy_xyyz_1, tey_yyy_xyzz_0, tey_yyy_xyzz_1, tey_yyy_xzzz_0, tey_yyy_xzzz_1, \
                                         tey_yyy_yyyy_0, tey_yyy_yyyy_1, tey_yyy_yyyz_0, tey_yyy_yyyz_1, tey_yyy_yyz_0, \
                                         tey_yyy_yyz_1, tey_yyy_yyzz_0, tey_yyy_yyzz_1, tey_yyy_yzz_0, tey_yyy_yzz_1, \
                                         tey_yyy_yzzz_0, tey_yyy_yzzz_1, tey_yyy_zzz_0, tey_yyy_zzz_1, tey_yyy_zzzz_0, \
                                         tey_yyy_zzzz_1, tey_yyz_xxx_0, tey_yyz_xxx_1, tey_yyz_xxxx_0, tey_yyz_xxxx_1, \
                                         tey_yyz_xxxy_0, tey_yyz_xxxy_1, tey_yyz_xxxz_0, tey_yyz_xxxz_1, tey_yyz_xxy_0, \
                                         tey_yyz_xxy_1, tey_yyz_xxyy_0, tey_yyz_xxyy_1, tey_yyz_xxyz_0, tey_yyz_xxyz_1, \
                                         tey_yyz_xxz_0, tey_yyz_xxz_1, tey_yyz_xxzz_0, tey_yyz_xxzz_1, tey_yyz_xyy_0, \
                                         tey_yyz_xyy_1, tey_yyz_xyyy_0, tey_yyz_xyyy_1, tey_yyz_xyyz_0, tey_yyz_xyyz_1, \
                                         tey_yyz_xyz_0, tey_yyz_xyz_1, tey_yyz_xzz_0, tey_yyz_xzz_1, tey_yyz_yyy_0, \
                                         tey_yyz_yyy_1, tey_yyz_yyz_0, tey_yyz_yyz_1, tez_xyyy_xyyz_0, tez_xyyy_xyzz_0, \
                                         tez_xyyy_xzzz_0, tez_xyyy_yyyy_0, tez_xyyy_yyyz_0, tez_xyyy_yyzz_0, tez_xyyy_yzzz_0, \
                                         tez_xyyy_zzzz_0, tez_xyyz_xxxx_0, tez_xyyz_xxxy_0, tez_xyyz_xxxz_0, tez_xyyz_xxyy_0, \
                                         tez_xyyz_xxyz_0, tez_xyyz_xxzz_0, tez_xyyz_xyyy_0, tez_xyyz_xyyz_0, tez_yyy_xyyz_0, \
                                         tez_yyy_xyyz_1, tez_yyy_xyzz_0, tez_yyy_xyzz_1, tez_yyy_xzzz_0, tez_yyy_xzzz_1, \
                                         tez_yyy_yyyy_0, tez_yyy_yyyy_1, tez_yyy_yyyz_0, tez_yyy_yyyz_1, tez_yyy_yyz_0, \
                                         tez_yyy_yyz_1, tez_yyy_yyzz_0, tez_yyy_yyzz_1, tez_yyy_yzz_0, tez_yyy_yzz_1, \
                                         tez_yyy_yzzz_0, tez_yyy_yzzz_1, tez_yyy_zzz_0, tez_yyy_zzz_1, tez_yyy_zzzz_0, \
                                         tez_yyy_zzzz_1, tez_yyz_xxx_0, tez_yyz_xxx_1, tez_yyz_xxxx_0, tez_yyz_xxxx_1, \
                                         tez_yyz_xxxy_0, tez_yyz_xxxy_1, tez_yyz_xxxz_0, tez_yyz_xxxz_1, tez_yyz_xxy_0, \
                                         tez_yyz_xxy_1, tez_yyz_xxyy_0, tez_yyz_xxyy_1, tez_yyz_xxyz_0, tez_yyz_xxyz_1, \
                                         tez_yyz_xxz_0, tez_yyz_xxz_1, tez_yyz_xxzz_0, tez_yyz_xxzz_1, tez_yyz_xyy_0, \
                                         tez_yyz_xyy_1, tez_yyz_xyyy_0, tez_yyz_xyyy_1, tez_yyz_xyyz_0, tez_yyz_xyyz_1, \
                                         tez_yyz_xyz_0, tez_yyz_xyz_1, tez_yyz_xzz_0, tez_yyz_xzz_1, tez_yyz_yyy_0, \
                                         tez_yyz_yyy_1, tez_yyz_yyz_0, tez_yyz_yyz_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_xyyy_xyyz_0[j] = pa_x[j] * tex_yyy_xyyz_0[j] - pc_x[j] * tex_yyy_xyyz_1[j] + 0.5 * fl1_fx * tex_yyy_yyz_0[j] -
                                     0.5 * fl1_fx * tex_yyy_yyz_1[j] + ta_yyy_xyyz_1[j];

                tey_xyyy_xyyz_0[j] =
                    pa_x[j] * tey_yyy_xyyz_0[j] - pc_x[j] * tey_yyy_xyyz_1[j] + 0.5 * fl1_fx * tey_yyy_yyz_0[j] - 0.5 * fl1_fx * tey_yyy_yyz_1[j];

                tez_xyyy_xyyz_0[j] =
                    pa_x[j] * tez_yyy_xyyz_0[j] - pc_x[j] * tez_yyy_xyyz_1[j] + 0.5 * fl1_fx * tez_yyy_yyz_0[j] - 0.5 * fl1_fx * tez_yyy_yyz_1[j];

                tex_xyyy_xyzz_0[j] = pa_x[j] * tex_yyy_xyzz_0[j] - pc_x[j] * tex_yyy_xyzz_1[j] + 0.5 * fl1_fx * tex_yyy_yzz_0[j] -
                                     0.5 * fl1_fx * tex_yyy_yzz_1[j] + ta_yyy_xyzz_1[j];

                tey_xyyy_xyzz_0[j] =
                    pa_x[j] * tey_yyy_xyzz_0[j] - pc_x[j] * tey_yyy_xyzz_1[j] + 0.5 * fl1_fx * tey_yyy_yzz_0[j] - 0.5 * fl1_fx * tey_yyy_yzz_1[j];

                tez_xyyy_xyzz_0[j] =
                    pa_x[j] * tez_yyy_xyzz_0[j] - pc_x[j] * tez_yyy_xyzz_1[j] + 0.5 * fl1_fx * tez_yyy_yzz_0[j] - 0.5 * fl1_fx * tez_yyy_yzz_1[j];

                tex_xyyy_xzzz_0[j] = pa_x[j] * tex_yyy_xzzz_0[j] - pc_x[j] * tex_yyy_xzzz_1[j] + 0.5 * fl1_fx * tex_yyy_zzz_0[j] -
                                     0.5 * fl1_fx * tex_yyy_zzz_1[j] + ta_yyy_xzzz_1[j];

                tey_xyyy_xzzz_0[j] =
                    pa_x[j] * tey_yyy_xzzz_0[j] - pc_x[j] * tey_yyy_xzzz_1[j] + 0.5 * fl1_fx * tey_yyy_zzz_0[j] - 0.5 * fl1_fx * tey_yyy_zzz_1[j];

                tez_xyyy_xzzz_0[j] =
                    pa_x[j] * tez_yyy_xzzz_0[j] - pc_x[j] * tez_yyy_xzzz_1[j] + 0.5 * fl1_fx * tez_yyy_zzz_0[j] - 0.5 * fl1_fx * tez_yyy_zzz_1[j];

                tex_xyyy_yyyy_0[j] = pa_x[j] * tex_yyy_yyyy_0[j] - pc_x[j] * tex_yyy_yyyy_1[j] + ta_yyy_yyyy_1[j];

                tey_xyyy_yyyy_0[j] = pa_x[j] * tey_yyy_yyyy_0[j] - pc_x[j] * tey_yyy_yyyy_1[j];

                tez_xyyy_yyyy_0[j] = pa_x[j] * tez_yyy_yyyy_0[j] - pc_x[j] * tez_yyy_yyyy_1[j];

                tex_xyyy_yyyz_0[j] = pa_x[j] * tex_yyy_yyyz_0[j] - pc_x[j] * tex_yyy_yyyz_1[j] + ta_yyy_yyyz_1[j];

                tey_xyyy_yyyz_0[j] = pa_x[j] * tey_yyy_yyyz_0[j] - pc_x[j] * tey_yyy_yyyz_1[j];

                tez_xyyy_yyyz_0[j] = pa_x[j] * tez_yyy_yyyz_0[j] - pc_x[j] * tez_yyy_yyyz_1[j];

                tex_xyyy_yyzz_0[j] = pa_x[j] * tex_yyy_yyzz_0[j] - pc_x[j] * tex_yyy_yyzz_1[j] + ta_yyy_yyzz_1[j];

                tey_xyyy_yyzz_0[j] = pa_x[j] * tey_yyy_yyzz_0[j] - pc_x[j] * tey_yyy_yyzz_1[j];

                tez_xyyy_yyzz_0[j] = pa_x[j] * tez_yyy_yyzz_0[j] - pc_x[j] * tez_yyy_yyzz_1[j];

                tex_xyyy_yzzz_0[j] = pa_x[j] * tex_yyy_yzzz_0[j] - pc_x[j] * tex_yyy_yzzz_1[j] + ta_yyy_yzzz_1[j];

                tey_xyyy_yzzz_0[j] = pa_x[j] * tey_yyy_yzzz_0[j] - pc_x[j] * tey_yyy_yzzz_1[j];

                tez_xyyy_yzzz_0[j] = pa_x[j] * tez_yyy_yzzz_0[j] - pc_x[j] * tez_yyy_yzzz_1[j];

                tex_xyyy_zzzz_0[j] = pa_x[j] * tex_yyy_zzzz_0[j] - pc_x[j] * tex_yyy_zzzz_1[j] + ta_yyy_zzzz_1[j];

                tey_xyyy_zzzz_0[j] = pa_x[j] * tey_yyy_zzzz_0[j] - pc_x[j] * tey_yyy_zzzz_1[j];

                tez_xyyy_zzzz_0[j] = pa_x[j] * tez_yyy_zzzz_0[j] - pc_x[j] * tez_yyy_zzzz_1[j];

                tex_xyyz_xxxx_0[j] = pa_x[j] * tex_yyz_xxxx_0[j] - pc_x[j] * tex_yyz_xxxx_1[j] + 2.0 * fl1_fx * tex_yyz_xxx_0[j] -
                                     2.0 * fl1_fx * tex_yyz_xxx_1[j] + ta_yyz_xxxx_1[j];

                tey_xyyz_xxxx_0[j] =
                    pa_x[j] * tey_yyz_xxxx_0[j] - pc_x[j] * tey_yyz_xxxx_1[j] + 2.0 * fl1_fx * tey_yyz_xxx_0[j] - 2.0 * fl1_fx * tey_yyz_xxx_1[j];

                tez_xyyz_xxxx_0[j] =
                    pa_x[j] * tez_yyz_xxxx_0[j] - pc_x[j] * tez_yyz_xxxx_1[j] + 2.0 * fl1_fx * tez_yyz_xxx_0[j] - 2.0 * fl1_fx * tez_yyz_xxx_1[j];

                tex_xyyz_xxxy_0[j] = pa_x[j] * tex_yyz_xxxy_0[j] - pc_x[j] * tex_yyz_xxxy_1[j] + 1.5 * fl1_fx * tex_yyz_xxy_0[j] -
                                     1.5 * fl1_fx * tex_yyz_xxy_1[j] + ta_yyz_xxxy_1[j];

                tey_xyyz_xxxy_0[j] =
                    pa_x[j] * tey_yyz_xxxy_0[j] - pc_x[j] * tey_yyz_xxxy_1[j] + 1.5 * fl1_fx * tey_yyz_xxy_0[j] - 1.5 * fl1_fx * tey_yyz_xxy_1[j];

                tez_xyyz_xxxy_0[j] =
                    pa_x[j] * tez_yyz_xxxy_0[j] - pc_x[j] * tez_yyz_xxxy_1[j] + 1.5 * fl1_fx * tez_yyz_xxy_0[j] - 1.5 * fl1_fx * tez_yyz_xxy_1[j];

                tex_xyyz_xxxz_0[j] = pa_x[j] * tex_yyz_xxxz_0[j] - pc_x[j] * tex_yyz_xxxz_1[j] + 1.5 * fl1_fx * tex_yyz_xxz_0[j] -
                                     1.5 * fl1_fx * tex_yyz_xxz_1[j] + ta_yyz_xxxz_1[j];

                tey_xyyz_xxxz_0[j] =
                    pa_x[j] * tey_yyz_xxxz_0[j] - pc_x[j] * tey_yyz_xxxz_1[j] + 1.5 * fl1_fx * tey_yyz_xxz_0[j] - 1.5 * fl1_fx * tey_yyz_xxz_1[j];

                tez_xyyz_xxxz_0[j] =
                    pa_x[j] * tez_yyz_xxxz_0[j] - pc_x[j] * tez_yyz_xxxz_1[j] + 1.5 * fl1_fx * tez_yyz_xxz_0[j] - 1.5 * fl1_fx * tez_yyz_xxz_1[j];

                tex_xyyz_xxyy_0[j] = pa_x[j] * tex_yyz_xxyy_0[j] - pc_x[j] * tex_yyz_xxyy_1[j] + fl1_fx * tex_yyz_xyy_0[j] -
                                     fl1_fx * tex_yyz_xyy_1[j] + ta_yyz_xxyy_1[j];

                tey_xyyz_xxyy_0[j] =
                    pa_x[j] * tey_yyz_xxyy_0[j] - pc_x[j] * tey_yyz_xxyy_1[j] + fl1_fx * tey_yyz_xyy_0[j] - fl1_fx * tey_yyz_xyy_1[j];

                tez_xyyz_xxyy_0[j] =
                    pa_x[j] * tez_yyz_xxyy_0[j] - pc_x[j] * tez_yyz_xxyy_1[j] + fl1_fx * tez_yyz_xyy_0[j] - fl1_fx * tez_yyz_xyy_1[j];

                tex_xyyz_xxyz_0[j] = pa_x[j] * tex_yyz_xxyz_0[j] - pc_x[j] * tex_yyz_xxyz_1[j] + fl1_fx * tex_yyz_xyz_0[j] -
                                     fl1_fx * tex_yyz_xyz_1[j] + ta_yyz_xxyz_1[j];

                tey_xyyz_xxyz_0[j] =
                    pa_x[j] * tey_yyz_xxyz_0[j] - pc_x[j] * tey_yyz_xxyz_1[j] + fl1_fx * tey_yyz_xyz_0[j] - fl1_fx * tey_yyz_xyz_1[j];

                tez_xyyz_xxyz_0[j] =
                    pa_x[j] * tez_yyz_xxyz_0[j] - pc_x[j] * tez_yyz_xxyz_1[j] + fl1_fx * tez_yyz_xyz_0[j] - fl1_fx * tez_yyz_xyz_1[j];

                tex_xyyz_xxzz_0[j] = pa_x[j] * tex_yyz_xxzz_0[j] - pc_x[j] * tex_yyz_xxzz_1[j] + fl1_fx * tex_yyz_xzz_0[j] -
                                     fl1_fx * tex_yyz_xzz_1[j] + ta_yyz_xxzz_1[j];

                tey_xyyz_xxzz_0[j] =
                    pa_x[j] * tey_yyz_xxzz_0[j] - pc_x[j] * tey_yyz_xxzz_1[j] + fl1_fx * tey_yyz_xzz_0[j] - fl1_fx * tey_yyz_xzz_1[j];

                tez_xyyz_xxzz_0[j] =
                    pa_x[j] * tez_yyz_xxzz_0[j] - pc_x[j] * tez_yyz_xxzz_1[j] + fl1_fx * tez_yyz_xzz_0[j] - fl1_fx * tez_yyz_xzz_1[j];

                tex_xyyz_xyyy_0[j] = pa_x[j] * tex_yyz_xyyy_0[j] - pc_x[j] * tex_yyz_xyyy_1[j] + 0.5 * fl1_fx * tex_yyz_yyy_0[j] -
                                     0.5 * fl1_fx * tex_yyz_yyy_1[j] + ta_yyz_xyyy_1[j];

                tey_xyyz_xyyy_0[j] =
                    pa_x[j] * tey_yyz_xyyy_0[j] - pc_x[j] * tey_yyz_xyyy_1[j] + 0.5 * fl1_fx * tey_yyz_yyy_0[j] - 0.5 * fl1_fx * tey_yyz_yyy_1[j];

                tez_xyyz_xyyy_0[j] =
                    pa_x[j] * tez_yyz_xyyy_0[j] - pc_x[j] * tez_yyz_xyyy_1[j] + 0.5 * fl1_fx * tez_yyz_yyy_0[j] - 0.5 * fl1_fx * tez_yyz_yyy_1[j];

                tex_xyyz_xyyz_0[j] = pa_x[j] * tex_yyz_xyyz_0[j] - pc_x[j] * tex_yyz_xyyz_1[j] + 0.5 * fl1_fx * tex_yyz_yyz_0[j] -
                                     0.5 * fl1_fx * tex_yyz_yyz_1[j] + ta_yyz_xyyz_1[j];

                tey_xyyz_xyyz_0[j] =
                    pa_x[j] * tey_yyz_xyyz_0[j] - pc_x[j] * tey_yyz_xyyz_1[j] + 0.5 * fl1_fx * tey_yyz_yyz_0[j] - 0.5 * fl1_fx * tey_yyz_yyz_1[j];

                tez_xyyz_xyyz_0[j] =
                    pa_x[j] * tez_yyz_xyyz_0[j] - pc_x[j] * tez_yyz_xyyz_1[j] + 0.5 * fl1_fx * tez_yyz_yyz_0[j] - 0.5 * fl1_fx * tez_yyz_yyz_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForGG_339_387(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pcDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto)
{
    // Batch of Integrals (339,387)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_4_4_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tex_yyz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 113);

            auto tey_yyz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 113);

            auto tez_yyz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 113);

            auto tex_yyz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 114);

            auto tey_yyz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 114);

            auto tez_yyz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 114);

            auto tex_yyz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 115);

            auto tey_yyz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 115);

            auto tez_yyz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 115);

            auto tex_yyz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 116);

            auto tey_yyz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 116);

            auto tez_yyz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 116);

            auto tex_yyz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 117);

            auto tey_yyz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 117);

            auto tez_yyz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 117);

            auto tex_yyz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 118);

            auto tey_yyz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 118);

            auto tez_yyz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 118);

            auto tex_yyz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 119);

            auto tey_yyz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 119);

            auto tez_yyz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 119);

            auto tex_yzz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 120);

            auto tey_yzz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 120);

            auto tez_yzz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 120);

            auto tex_yzz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 121);

            auto tey_yzz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 121);

            auto tez_yzz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 121);

            auto tex_yzz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 122);

            auto tey_yzz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 122);

            auto tez_yzz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 122);

            auto tex_yzz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 123);

            auto tey_yzz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 123);

            auto tez_yzz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 123);

            auto tex_yzz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 124);

            auto tey_yzz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 124);

            auto tez_yzz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 124);

            auto tex_yzz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 125);

            auto tey_yzz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 125);

            auto tez_yzz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 125);

            auto tex_yzz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 126);

            auto tey_yzz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 126);

            auto tez_yzz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 126);

            auto tex_yzz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 127);

            auto tey_yzz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 127);

            auto tez_yzz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 127);

            auto tex_yzz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 128);

            auto tey_yzz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 128);

            auto tez_yzz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 128);

            auto tex_yyz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 113);

            auto tey_yyz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 113);

            auto tez_yyz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 113);

            auto tex_yyz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 114);

            auto tey_yyz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 114);

            auto tez_yyz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 114);

            auto tex_yyz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 115);

            auto tey_yyz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 115);

            auto tez_yyz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 115);

            auto tex_yyz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 116);

            auto tey_yyz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 116);

            auto tez_yyz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 116);

            auto tex_yyz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 117);

            auto tey_yyz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 117);

            auto tez_yyz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 117);

            auto tex_yyz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 118);

            auto tey_yyz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 118);

            auto tez_yyz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 118);

            auto tex_yyz_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 119);

            auto tey_yyz_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 119);

            auto tez_yyz_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 119);

            auto tex_yzz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 120);

            auto tey_yzz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 120);

            auto tez_yzz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 120);

            auto tex_yzz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 121);

            auto tey_yzz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 121);

            auto tez_yzz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 121);

            auto tex_yzz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 122);

            auto tey_yzz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 122);

            auto tez_yzz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 122);

            auto tex_yzz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 123);

            auto tey_yzz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 123);

            auto tez_yzz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 123);

            auto tex_yzz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 124);

            auto tey_yzz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 124);

            auto tez_yzz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 124);

            auto tex_yzz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 125);

            auto tey_yzz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 125);

            auto tez_yzz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 125);

            auto tex_yzz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 126);

            auto tey_yzz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 126);

            auto tez_yzz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 126);

            auto tex_yzz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 127);

            auto tey_yzz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 127);

            auto tez_yzz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 127);

            auto tex_yzz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 128);

            auto tey_yzz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 128);

            auto tez_yzz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 128);

            auto tex_yyz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 78);

            auto tey_yyz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 78);

            auto tez_yyz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 78);

            auto tex_yyz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 79);

            auto tey_yyz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 79);

            auto tez_yyz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 79);

            auto tex_yzz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 80);

            auto tey_yzz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 80);

            auto tez_yzz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 80);

            auto tex_yzz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 81);

            auto tey_yzz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 81);

            auto tez_yzz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 81);

            auto tex_yzz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 82);

            auto tey_yzz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 82);

            auto tez_yzz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 82);

            auto tex_yzz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 83);

            auto tey_yzz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 83);

            auto tez_yzz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 83);

            auto tex_yzz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 84);

            auto tey_yzz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 84);

            auto tez_yzz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 84);

            auto tex_yzz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 85);

            auto tey_yzz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 85);

            auto tez_yzz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 85);

            auto tex_yzz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 86);

            auto tey_yzz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 86);

            auto tez_yzz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 86);

            auto tex_yzz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 87);

            auto tey_yzz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 87);

            auto tez_yzz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 87);

            auto tex_yzz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 88);

            auto tey_yzz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 88);

            auto tez_yzz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 88);

            auto tex_yyz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 78);

            auto tey_yyz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 78);

            auto tez_yyz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 78);

            auto tex_yyz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 79);

            auto tey_yyz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 79);

            auto tez_yyz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 79);

            auto tex_yzz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 80);

            auto tey_yzz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 80);

            auto tez_yzz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 80);

            auto tex_yzz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 81);

            auto tey_yzz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 81);

            auto tez_yzz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 81);

            auto tex_yzz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 82);

            auto tey_yzz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 82);

            auto tez_yzz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 82);

            auto tex_yzz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 83);

            auto tey_yzz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 83);

            auto tez_yzz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 83);

            auto tex_yzz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 84);

            auto tey_yzz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 84);

            auto tez_yzz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 84);

            auto tex_yzz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 85);

            auto tey_yzz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 85);

            auto tez_yzz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 85);

            auto tex_yzz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 86);

            auto tey_yzz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 86);

            auto tez_yzz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 86);

            auto tex_yzz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 87);

            auto tey_yzz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 87);

            auto tez_yzz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 87);

            auto tex_yzz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 88);

            auto tey_yzz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 88);

            auto tez_yzz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 88);

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

            // set up pointers to integrals

            auto tex_xyyz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 113);

            auto tey_xyyz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 113);

            auto tez_xyyz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 113);

            auto tex_xyyz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 114);

            auto tey_xyyz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 114);

            auto tez_xyyz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 114);

            auto tex_xyyz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 115);

            auto tey_xyyz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 115);

            auto tez_xyyz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 115);

            auto tex_xyyz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 116);

            auto tey_xyyz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 116);

            auto tez_xyyz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 116);

            auto tex_xyyz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 117);

            auto tey_xyyz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 117);

            auto tez_xyyz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 117);

            auto tex_xyyz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 118);

            auto tey_xyyz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 118);

            auto tez_xyyz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 118);

            auto tex_xyyz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 119);

            auto tey_xyyz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 119);

            auto tez_xyyz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 119);

            auto tex_xyzz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 120);

            auto tey_xyzz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 120);

            auto tez_xyzz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 120);

            auto tex_xyzz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 121);

            auto tey_xyzz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 121);

            auto tez_xyzz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 121);

            auto tex_xyzz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 122);

            auto tey_xyzz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 122);

            auto tez_xyzz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 122);

            auto tex_xyzz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 123);

            auto tey_xyzz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 123);

            auto tez_xyzz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 123);

            auto tex_xyzz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 124);

            auto tey_xyzz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 124);

            auto tez_xyzz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 124);

            auto tex_xyzz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 125);

            auto tey_xyzz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 125);

            auto tez_xyzz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 125);

            auto tex_xyzz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 126);

            auto tey_xyzz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 126);

            auto tez_xyzz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 126);

            auto tex_xyzz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 127);

            auto tey_xyzz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 127);

            auto tez_xyzz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 127);

            auto tex_xyzz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 128);

            auto tey_xyzz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 128);

            auto tez_xyzz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 128);

            // Batch of Integrals (339,387)

            #pragma omp simd aligned(fx, pa_x, pc_x, ta_yyz_xyzz_1, ta_yyz_xzzz_1, ta_yyz_yyyy_1, \
                                         ta_yyz_yyyz_1, ta_yyz_yyzz_1, ta_yyz_yzzz_1, ta_yyz_zzzz_1, ta_yzz_xxxx_1, \
                                         ta_yzz_xxxy_1, ta_yzz_xxxz_1, ta_yzz_xxyy_1, ta_yzz_xxyz_1, ta_yzz_xxzz_1, \
                                         ta_yzz_xyyy_1, ta_yzz_xyyz_1, ta_yzz_xyzz_1, tex_xyyz_xyzz_0, tex_xyyz_xzzz_0, \
                                         tex_xyyz_yyyy_0, tex_xyyz_yyyz_0, tex_xyyz_yyzz_0, tex_xyyz_yzzz_0, tex_xyyz_zzzz_0, \
                                         tex_xyzz_xxxx_0, tex_xyzz_xxxy_0, tex_xyzz_xxxz_0, tex_xyzz_xxyy_0, tex_xyzz_xxyz_0, \
                                         tex_xyzz_xxzz_0, tex_xyzz_xyyy_0, tex_xyzz_xyyz_0, tex_xyzz_xyzz_0, tex_yyz_xyzz_0, \
                                         tex_yyz_xyzz_1, tex_yyz_xzzz_0, tex_yyz_xzzz_1, tex_yyz_yyyy_0, tex_yyz_yyyy_1, \
                                         tex_yyz_yyyz_0, tex_yyz_yyyz_1, tex_yyz_yyzz_0, tex_yyz_yyzz_1, tex_yyz_yzz_0, \
                                         tex_yyz_yzz_1, tex_yyz_yzzz_0, tex_yyz_yzzz_1, tex_yyz_zzz_0, tex_yyz_zzz_1, \
                                         tex_yyz_zzzz_0, tex_yyz_zzzz_1, tex_yzz_xxx_0, tex_yzz_xxx_1, tex_yzz_xxxx_0, \
                                         tex_yzz_xxxx_1, tex_yzz_xxxy_0, tex_yzz_xxxy_1, tex_yzz_xxxz_0, tex_yzz_xxxz_1, \
                                         tex_yzz_xxy_0, tex_yzz_xxy_1, tex_yzz_xxyy_0, tex_yzz_xxyy_1, tex_yzz_xxyz_0, \
                                         tex_yzz_xxyz_1, tex_yzz_xxz_0, tex_yzz_xxz_1, tex_yzz_xxzz_0, tex_yzz_xxzz_1, \
                                         tex_yzz_xyy_0, tex_yzz_xyy_1, tex_yzz_xyyy_0, tex_yzz_xyyy_1, tex_yzz_xyyz_0, \
                                         tex_yzz_xyyz_1, tex_yzz_xyz_0, tex_yzz_xyz_1, tex_yzz_xyzz_0, tex_yzz_xyzz_1, \
                                         tex_yzz_xzz_0, tex_yzz_xzz_1, tex_yzz_yyy_0, tex_yzz_yyy_1, tex_yzz_yyz_0, \
                                         tex_yzz_yyz_1, tex_yzz_yzz_0, tex_yzz_yzz_1, tey_xyyz_xyzz_0, tey_xyyz_xzzz_0, \
                                         tey_xyyz_yyyy_0, tey_xyyz_yyyz_0, tey_xyyz_yyzz_0, tey_xyyz_yzzz_0, tey_xyyz_zzzz_0, \
                                         tey_xyzz_xxxx_0, tey_xyzz_xxxy_0, tey_xyzz_xxxz_0, tey_xyzz_xxyy_0, tey_xyzz_xxyz_0, \
                                         tey_xyzz_xxzz_0, tey_xyzz_xyyy_0, tey_xyzz_xyyz_0, tey_xyzz_xyzz_0, tey_yyz_xyzz_0, \
                                         tey_yyz_xyzz_1, tey_yyz_xzzz_0, tey_yyz_xzzz_1, tey_yyz_yyyy_0, tey_yyz_yyyy_1, \
                                         tey_yyz_yyyz_0, tey_yyz_yyyz_1, tey_yyz_yyzz_0, tey_yyz_yyzz_1, tey_yyz_yzz_0, \
                                         tey_yyz_yzz_1, tey_yyz_yzzz_0, tey_yyz_yzzz_1, tey_yyz_zzz_0, tey_yyz_zzz_1, \
                                         tey_yyz_zzzz_0, tey_yyz_zzzz_1, tey_yzz_xxx_0, tey_yzz_xxx_1, tey_yzz_xxxx_0, \
                                         tey_yzz_xxxx_1, tey_yzz_xxxy_0, tey_yzz_xxxy_1, tey_yzz_xxxz_0, tey_yzz_xxxz_1, \
                                         tey_yzz_xxy_0, tey_yzz_xxy_1, tey_yzz_xxyy_0, tey_yzz_xxyy_1, tey_yzz_xxyz_0, \
                                         tey_yzz_xxyz_1, tey_yzz_xxz_0, tey_yzz_xxz_1, tey_yzz_xxzz_0, tey_yzz_xxzz_1, \
                                         tey_yzz_xyy_0, tey_yzz_xyy_1, tey_yzz_xyyy_0, tey_yzz_xyyy_1, tey_yzz_xyyz_0, \
                                         tey_yzz_xyyz_1, tey_yzz_xyz_0, tey_yzz_xyz_1, tey_yzz_xyzz_0, tey_yzz_xyzz_1, \
                                         tey_yzz_xzz_0, tey_yzz_xzz_1, tey_yzz_yyy_0, tey_yzz_yyy_1, tey_yzz_yyz_0, \
                                         tey_yzz_yyz_1, tey_yzz_yzz_0, tey_yzz_yzz_1, tez_xyyz_xyzz_0, tez_xyyz_xzzz_0, \
                                         tez_xyyz_yyyy_0, tez_xyyz_yyyz_0, tez_xyyz_yyzz_0, tez_xyyz_yzzz_0, tez_xyyz_zzzz_0, \
                                         tez_xyzz_xxxx_0, tez_xyzz_xxxy_0, tez_xyzz_xxxz_0, tez_xyzz_xxyy_0, tez_xyzz_xxyz_0, \
                                         tez_xyzz_xxzz_0, tez_xyzz_xyyy_0, tez_xyzz_xyyz_0, tez_xyzz_xyzz_0, tez_yyz_xyzz_0, \
                                         tez_yyz_xyzz_1, tez_yyz_xzzz_0, tez_yyz_xzzz_1, tez_yyz_yyyy_0, tez_yyz_yyyy_1, \
                                         tez_yyz_yyyz_0, tez_yyz_yyyz_1, tez_yyz_yyzz_0, tez_yyz_yyzz_1, tez_yyz_yzz_0, \
                                         tez_yyz_yzz_1, tez_yyz_yzzz_0, tez_yyz_yzzz_1, tez_yyz_zzz_0, tez_yyz_zzz_1, \
                                         tez_yyz_zzzz_0, tez_yyz_zzzz_1, tez_yzz_xxx_0, tez_yzz_xxx_1, tez_yzz_xxxx_0, \
                                         tez_yzz_xxxx_1, tez_yzz_xxxy_0, tez_yzz_xxxy_1, tez_yzz_xxxz_0, tez_yzz_xxxz_1, \
                                         tez_yzz_xxy_0, tez_yzz_xxy_1, tez_yzz_xxyy_0, tez_yzz_xxyy_1, tez_yzz_xxyz_0, \
                                         tez_yzz_xxyz_1, tez_yzz_xxz_0, tez_yzz_xxz_1, tez_yzz_xxzz_0, tez_yzz_xxzz_1, \
                                         tez_yzz_xyy_0, tez_yzz_xyy_1, tez_yzz_xyyy_0, tez_yzz_xyyy_1, tez_yzz_xyyz_0, \
                                         tez_yzz_xyyz_1, tez_yzz_xyz_0, tez_yzz_xyz_1, tez_yzz_xyzz_0, tez_yzz_xyzz_1, \
                                         tez_yzz_xzz_0, tez_yzz_xzz_1, tez_yzz_yyy_0, tez_yzz_yyy_1, tez_yzz_yyz_0, \
                                         tez_yzz_yyz_1, tez_yzz_yzz_0, tez_yzz_yzz_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_xyyz_xyzz_0[j] = pa_x[j] * tex_yyz_xyzz_0[j] - pc_x[j] * tex_yyz_xyzz_1[j] + 0.5 * fl1_fx * tex_yyz_yzz_0[j] -
                                     0.5 * fl1_fx * tex_yyz_yzz_1[j] + ta_yyz_xyzz_1[j];

                tey_xyyz_xyzz_0[j] =
                    pa_x[j] * tey_yyz_xyzz_0[j] - pc_x[j] * tey_yyz_xyzz_1[j] + 0.5 * fl1_fx * tey_yyz_yzz_0[j] - 0.5 * fl1_fx * tey_yyz_yzz_1[j];

                tez_xyyz_xyzz_0[j] =
                    pa_x[j] * tez_yyz_xyzz_0[j] - pc_x[j] * tez_yyz_xyzz_1[j] + 0.5 * fl1_fx * tez_yyz_yzz_0[j] - 0.5 * fl1_fx * tez_yyz_yzz_1[j];

                tex_xyyz_xzzz_0[j] = pa_x[j] * tex_yyz_xzzz_0[j] - pc_x[j] * tex_yyz_xzzz_1[j] + 0.5 * fl1_fx * tex_yyz_zzz_0[j] -
                                     0.5 * fl1_fx * tex_yyz_zzz_1[j] + ta_yyz_xzzz_1[j];

                tey_xyyz_xzzz_0[j] =
                    pa_x[j] * tey_yyz_xzzz_0[j] - pc_x[j] * tey_yyz_xzzz_1[j] + 0.5 * fl1_fx * tey_yyz_zzz_0[j] - 0.5 * fl1_fx * tey_yyz_zzz_1[j];

                tez_xyyz_xzzz_0[j] =
                    pa_x[j] * tez_yyz_xzzz_0[j] - pc_x[j] * tez_yyz_xzzz_1[j] + 0.5 * fl1_fx * tez_yyz_zzz_0[j] - 0.5 * fl1_fx * tez_yyz_zzz_1[j];

                tex_xyyz_yyyy_0[j] = pa_x[j] * tex_yyz_yyyy_0[j] - pc_x[j] * tex_yyz_yyyy_1[j] + ta_yyz_yyyy_1[j];

                tey_xyyz_yyyy_0[j] = pa_x[j] * tey_yyz_yyyy_0[j] - pc_x[j] * tey_yyz_yyyy_1[j];

                tez_xyyz_yyyy_0[j] = pa_x[j] * tez_yyz_yyyy_0[j] - pc_x[j] * tez_yyz_yyyy_1[j];

                tex_xyyz_yyyz_0[j] = pa_x[j] * tex_yyz_yyyz_0[j] - pc_x[j] * tex_yyz_yyyz_1[j] + ta_yyz_yyyz_1[j];

                tey_xyyz_yyyz_0[j] = pa_x[j] * tey_yyz_yyyz_0[j] - pc_x[j] * tey_yyz_yyyz_1[j];

                tez_xyyz_yyyz_0[j] = pa_x[j] * tez_yyz_yyyz_0[j] - pc_x[j] * tez_yyz_yyyz_1[j];

                tex_xyyz_yyzz_0[j] = pa_x[j] * tex_yyz_yyzz_0[j] - pc_x[j] * tex_yyz_yyzz_1[j] + ta_yyz_yyzz_1[j];

                tey_xyyz_yyzz_0[j] = pa_x[j] * tey_yyz_yyzz_0[j] - pc_x[j] * tey_yyz_yyzz_1[j];

                tez_xyyz_yyzz_0[j] = pa_x[j] * tez_yyz_yyzz_0[j] - pc_x[j] * tez_yyz_yyzz_1[j];

                tex_xyyz_yzzz_0[j] = pa_x[j] * tex_yyz_yzzz_0[j] - pc_x[j] * tex_yyz_yzzz_1[j] + ta_yyz_yzzz_1[j];

                tey_xyyz_yzzz_0[j] = pa_x[j] * tey_yyz_yzzz_0[j] - pc_x[j] * tey_yyz_yzzz_1[j];

                tez_xyyz_yzzz_0[j] = pa_x[j] * tez_yyz_yzzz_0[j] - pc_x[j] * tez_yyz_yzzz_1[j];

                tex_xyyz_zzzz_0[j] = pa_x[j] * tex_yyz_zzzz_0[j] - pc_x[j] * tex_yyz_zzzz_1[j] + ta_yyz_zzzz_1[j];

                tey_xyyz_zzzz_0[j] = pa_x[j] * tey_yyz_zzzz_0[j] - pc_x[j] * tey_yyz_zzzz_1[j];

                tez_xyyz_zzzz_0[j] = pa_x[j] * tez_yyz_zzzz_0[j] - pc_x[j] * tez_yyz_zzzz_1[j];

                tex_xyzz_xxxx_0[j] = pa_x[j] * tex_yzz_xxxx_0[j] - pc_x[j] * tex_yzz_xxxx_1[j] + 2.0 * fl1_fx * tex_yzz_xxx_0[j] -
                                     2.0 * fl1_fx * tex_yzz_xxx_1[j] + ta_yzz_xxxx_1[j];

                tey_xyzz_xxxx_0[j] =
                    pa_x[j] * tey_yzz_xxxx_0[j] - pc_x[j] * tey_yzz_xxxx_1[j] + 2.0 * fl1_fx * tey_yzz_xxx_0[j] - 2.0 * fl1_fx * tey_yzz_xxx_1[j];

                tez_xyzz_xxxx_0[j] =
                    pa_x[j] * tez_yzz_xxxx_0[j] - pc_x[j] * tez_yzz_xxxx_1[j] + 2.0 * fl1_fx * tez_yzz_xxx_0[j] - 2.0 * fl1_fx * tez_yzz_xxx_1[j];

                tex_xyzz_xxxy_0[j] = pa_x[j] * tex_yzz_xxxy_0[j] - pc_x[j] * tex_yzz_xxxy_1[j] + 1.5 * fl1_fx * tex_yzz_xxy_0[j] -
                                     1.5 * fl1_fx * tex_yzz_xxy_1[j] + ta_yzz_xxxy_1[j];

                tey_xyzz_xxxy_0[j] =
                    pa_x[j] * tey_yzz_xxxy_0[j] - pc_x[j] * tey_yzz_xxxy_1[j] + 1.5 * fl1_fx * tey_yzz_xxy_0[j] - 1.5 * fl1_fx * tey_yzz_xxy_1[j];

                tez_xyzz_xxxy_0[j] =
                    pa_x[j] * tez_yzz_xxxy_0[j] - pc_x[j] * tez_yzz_xxxy_1[j] + 1.5 * fl1_fx * tez_yzz_xxy_0[j] - 1.5 * fl1_fx * tez_yzz_xxy_1[j];

                tex_xyzz_xxxz_0[j] = pa_x[j] * tex_yzz_xxxz_0[j] - pc_x[j] * tex_yzz_xxxz_1[j] + 1.5 * fl1_fx * tex_yzz_xxz_0[j] -
                                     1.5 * fl1_fx * tex_yzz_xxz_1[j] + ta_yzz_xxxz_1[j];

                tey_xyzz_xxxz_0[j] =
                    pa_x[j] * tey_yzz_xxxz_0[j] - pc_x[j] * tey_yzz_xxxz_1[j] + 1.5 * fl1_fx * tey_yzz_xxz_0[j] - 1.5 * fl1_fx * tey_yzz_xxz_1[j];

                tez_xyzz_xxxz_0[j] =
                    pa_x[j] * tez_yzz_xxxz_0[j] - pc_x[j] * tez_yzz_xxxz_1[j] + 1.5 * fl1_fx * tez_yzz_xxz_0[j] - 1.5 * fl1_fx * tez_yzz_xxz_1[j];

                tex_xyzz_xxyy_0[j] = pa_x[j] * tex_yzz_xxyy_0[j] - pc_x[j] * tex_yzz_xxyy_1[j] + fl1_fx * tex_yzz_xyy_0[j] -
                                     fl1_fx * tex_yzz_xyy_1[j] + ta_yzz_xxyy_1[j];

                tey_xyzz_xxyy_0[j] =
                    pa_x[j] * tey_yzz_xxyy_0[j] - pc_x[j] * tey_yzz_xxyy_1[j] + fl1_fx * tey_yzz_xyy_0[j] - fl1_fx * tey_yzz_xyy_1[j];

                tez_xyzz_xxyy_0[j] =
                    pa_x[j] * tez_yzz_xxyy_0[j] - pc_x[j] * tez_yzz_xxyy_1[j] + fl1_fx * tez_yzz_xyy_0[j] - fl1_fx * tez_yzz_xyy_1[j];

                tex_xyzz_xxyz_0[j] = pa_x[j] * tex_yzz_xxyz_0[j] - pc_x[j] * tex_yzz_xxyz_1[j] + fl1_fx * tex_yzz_xyz_0[j] -
                                     fl1_fx * tex_yzz_xyz_1[j] + ta_yzz_xxyz_1[j];

                tey_xyzz_xxyz_0[j] =
                    pa_x[j] * tey_yzz_xxyz_0[j] - pc_x[j] * tey_yzz_xxyz_1[j] + fl1_fx * tey_yzz_xyz_0[j] - fl1_fx * tey_yzz_xyz_1[j];

                tez_xyzz_xxyz_0[j] =
                    pa_x[j] * tez_yzz_xxyz_0[j] - pc_x[j] * tez_yzz_xxyz_1[j] + fl1_fx * tez_yzz_xyz_0[j] - fl1_fx * tez_yzz_xyz_1[j];

                tex_xyzz_xxzz_0[j] = pa_x[j] * tex_yzz_xxzz_0[j] - pc_x[j] * tex_yzz_xxzz_1[j] + fl1_fx * tex_yzz_xzz_0[j] -
                                     fl1_fx * tex_yzz_xzz_1[j] + ta_yzz_xxzz_1[j];

                tey_xyzz_xxzz_0[j] =
                    pa_x[j] * tey_yzz_xxzz_0[j] - pc_x[j] * tey_yzz_xxzz_1[j] + fl1_fx * tey_yzz_xzz_0[j] - fl1_fx * tey_yzz_xzz_1[j];

                tez_xyzz_xxzz_0[j] =
                    pa_x[j] * tez_yzz_xxzz_0[j] - pc_x[j] * tez_yzz_xxzz_1[j] + fl1_fx * tez_yzz_xzz_0[j] - fl1_fx * tez_yzz_xzz_1[j];

                tex_xyzz_xyyy_0[j] = pa_x[j] * tex_yzz_xyyy_0[j] - pc_x[j] * tex_yzz_xyyy_1[j] + 0.5 * fl1_fx * tex_yzz_yyy_0[j] -
                                     0.5 * fl1_fx * tex_yzz_yyy_1[j] + ta_yzz_xyyy_1[j];

                tey_xyzz_xyyy_0[j] =
                    pa_x[j] * tey_yzz_xyyy_0[j] - pc_x[j] * tey_yzz_xyyy_1[j] + 0.5 * fl1_fx * tey_yzz_yyy_0[j] - 0.5 * fl1_fx * tey_yzz_yyy_1[j];

                tez_xyzz_xyyy_0[j] =
                    pa_x[j] * tez_yzz_xyyy_0[j] - pc_x[j] * tez_yzz_xyyy_1[j] + 0.5 * fl1_fx * tez_yzz_yyy_0[j] - 0.5 * fl1_fx * tez_yzz_yyy_1[j];

                tex_xyzz_xyyz_0[j] = pa_x[j] * tex_yzz_xyyz_0[j] - pc_x[j] * tex_yzz_xyyz_1[j] + 0.5 * fl1_fx * tex_yzz_yyz_0[j] -
                                     0.5 * fl1_fx * tex_yzz_yyz_1[j] + ta_yzz_xyyz_1[j];

                tey_xyzz_xyyz_0[j] =
                    pa_x[j] * tey_yzz_xyyz_0[j] - pc_x[j] * tey_yzz_xyyz_1[j] + 0.5 * fl1_fx * tey_yzz_yyz_0[j] - 0.5 * fl1_fx * tey_yzz_yyz_1[j];

                tez_xyzz_xyyz_0[j] =
                    pa_x[j] * tez_yzz_xyyz_0[j] - pc_x[j] * tez_yzz_xyyz_1[j] + 0.5 * fl1_fx * tez_yzz_yyz_0[j] - 0.5 * fl1_fx * tez_yzz_yyz_1[j];

                tex_xyzz_xyzz_0[j] = pa_x[j] * tex_yzz_xyzz_0[j] - pc_x[j] * tex_yzz_xyzz_1[j] + 0.5 * fl1_fx * tex_yzz_yzz_0[j] -
                                     0.5 * fl1_fx * tex_yzz_yzz_1[j] + ta_yzz_xyzz_1[j];

                tey_xyzz_xyzz_0[j] =
                    pa_x[j] * tey_yzz_xyzz_0[j] - pc_x[j] * tey_yzz_xyzz_1[j] + 0.5 * fl1_fx * tey_yzz_yzz_0[j] - 0.5 * fl1_fx * tey_yzz_yzz_1[j];

                tez_xyzz_xyzz_0[j] =
                    pa_x[j] * tez_yzz_xyzz_0[j] - pc_x[j] * tez_yzz_xyzz_1[j] + 0.5 * fl1_fx * tez_yzz_yzz_0[j] - 0.5 * fl1_fx * tez_yzz_yzz_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForGG_387_435(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pcDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto)
{
    // Batch of Integrals (387,435)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_4_4_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tex_yzz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 129);

            auto tey_yzz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 129);

            auto tez_yzz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 129);

            auto tex_yzz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 130);

            auto tey_yzz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 130);

            auto tez_yzz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 130);

            auto tex_yzz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 131);

            auto tey_yzz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 131);

            auto tez_yzz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 131);

            auto tex_yzz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 132);

            auto tey_yzz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 132);

            auto tez_yzz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 132);

            auto tex_yzz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 133);

            auto tey_yzz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 133);

            auto tez_yzz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 133);

            auto tex_yzz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 134);

            auto tey_yzz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 134);

            auto tez_yzz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 134);

            auto tex_zzz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 135);

            auto tey_zzz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 135);

            auto tez_zzz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 135);

            auto tex_zzz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 136);

            auto tey_zzz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 136);

            auto tez_zzz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 136);

            auto tex_zzz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 137);

            auto tey_zzz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 137);

            auto tez_zzz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 137);

            auto tex_zzz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 138);

            auto tey_zzz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 138);

            auto tez_zzz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 138);

            auto tex_zzz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 139);

            auto tey_zzz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 139);

            auto tez_zzz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 139);

            auto tex_zzz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 140);

            auto tey_zzz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 140);

            auto tez_zzz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 140);

            auto tex_zzz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 141);

            auto tey_zzz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 141);

            auto tez_zzz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 141);

            auto tex_zzz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 142);

            auto tey_zzz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 142);

            auto tez_zzz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 142);

            auto tex_zzz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 143);

            auto tey_zzz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 143);

            auto tez_zzz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 143);

            auto tex_zzz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 144);

            auto tey_zzz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 144);

            auto tez_zzz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 144);

            auto tex_yzz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 129);

            auto tey_yzz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 129);

            auto tez_yzz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 129);

            auto tex_yzz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 130);

            auto tey_yzz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 130);

            auto tez_yzz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 130);

            auto tex_yzz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 131);

            auto tey_yzz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 131);

            auto tez_yzz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 131);

            auto tex_yzz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 132);

            auto tey_yzz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 132);

            auto tez_yzz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 132);

            auto tex_yzz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 133);

            auto tey_yzz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 133);

            auto tez_yzz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 133);

            auto tex_yzz_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 134);

            auto tey_yzz_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 134);

            auto tez_yzz_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 134);

            auto tex_zzz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 135);

            auto tey_zzz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 135);

            auto tez_zzz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 135);

            auto tex_zzz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 136);

            auto tey_zzz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 136);

            auto tez_zzz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 136);

            auto tex_zzz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 137);

            auto tey_zzz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 137);

            auto tez_zzz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 137);

            auto tex_zzz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 138);

            auto tey_zzz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 138);

            auto tez_zzz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 138);

            auto tex_zzz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 139);

            auto tey_zzz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 139);

            auto tez_zzz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 139);

            auto tex_zzz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 140);

            auto tey_zzz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 140);

            auto tez_zzz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 140);

            auto tex_zzz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 141);

            auto tey_zzz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 141);

            auto tez_zzz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 141);

            auto tex_zzz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 142);

            auto tey_zzz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 142);

            auto tez_zzz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 142);

            auto tex_zzz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 143);

            auto tey_zzz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 143);

            auto tez_zzz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 143);

            auto tex_zzz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 144);

            auto tey_zzz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 144);

            auto tez_zzz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 144);

            auto tex_yzz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 89);

            auto tey_yzz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 89);

            auto tez_yzz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 89);

            auto tex_zzz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 90);

            auto tey_zzz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 90);

            auto tez_zzz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 90);

            auto tex_zzz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 91);

            auto tey_zzz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 91);

            auto tez_zzz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 91);

            auto tex_zzz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 92);

            auto tey_zzz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 92);

            auto tez_zzz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 92);

            auto tex_zzz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 93);

            auto tey_zzz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 93);

            auto tez_zzz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 93);

            auto tex_zzz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 94);

            auto tey_zzz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 94);

            auto tez_zzz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 94);

            auto tex_zzz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 95);

            auto tey_zzz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 95);

            auto tez_zzz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 95);

            auto tex_zzz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 96);

            auto tey_zzz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 96);

            auto tez_zzz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 96);

            auto tex_zzz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 97);

            auto tey_zzz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 97);

            auto tez_zzz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 97);

            auto tex_zzz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 98);

            auto tey_zzz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 98);

            auto tez_zzz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 98);

            auto tex_zzz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 99);

            auto tey_zzz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 99);

            auto tez_zzz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 99);

            auto tex_yzz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 89);

            auto tey_yzz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 89);

            auto tez_yzz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 89);

            auto tex_zzz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 90);

            auto tey_zzz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 90);

            auto tez_zzz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 90);

            auto tex_zzz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 91);

            auto tey_zzz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 91);

            auto tez_zzz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 91);

            auto tex_zzz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 92);

            auto tey_zzz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 92);

            auto tez_zzz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 92);

            auto tex_zzz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 93);

            auto tey_zzz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 93);

            auto tez_zzz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 93);

            auto tex_zzz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 94);

            auto tey_zzz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 94);

            auto tez_zzz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 94);

            auto tex_zzz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 95);

            auto tey_zzz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 95);

            auto tez_zzz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 95);

            auto tex_zzz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 96);

            auto tey_zzz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 96);

            auto tez_zzz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 96);

            auto tex_zzz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 97);

            auto tey_zzz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 97);

            auto tez_zzz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 97);

            auto tex_zzz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 98);

            auto tey_zzz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 98);

            auto tez_zzz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 98);

            auto tex_zzz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 99);

            auto tey_zzz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 99);

            auto tez_zzz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 99);

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

            // set up pointers to integrals

            auto tex_xyzz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 129);

            auto tey_xyzz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 129);

            auto tez_xyzz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 129);

            auto tex_xyzz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 130);

            auto tey_xyzz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 130);

            auto tez_xyzz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 130);

            auto tex_xyzz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 131);

            auto tey_xyzz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 131);

            auto tez_xyzz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 131);

            auto tex_xyzz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 132);

            auto tey_xyzz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 132);

            auto tez_xyzz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 132);

            auto tex_xyzz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 133);

            auto tey_xyzz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 133);

            auto tez_xyzz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 133);

            auto tex_xyzz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 134);

            auto tey_xyzz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 134);

            auto tez_xyzz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 134);

            auto tex_xzzz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 135);

            auto tey_xzzz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 135);

            auto tez_xzzz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 135);

            auto tex_xzzz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 136);

            auto tey_xzzz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 136);

            auto tez_xzzz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 136);

            auto tex_xzzz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 137);

            auto tey_xzzz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 137);

            auto tez_xzzz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 137);

            auto tex_xzzz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 138);

            auto tey_xzzz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 138);

            auto tez_xzzz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 138);

            auto tex_xzzz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 139);

            auto tey_xzzz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 139);

            auto tez_xzzz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 139);

            auto tex_xzzz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 140);

            auto tey_xzzz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 140);

            auto tez_xzzz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 140);

            auto tex_xzzz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 141);

            auto tey_xzzz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 141);

            auto tez_xzzz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 141);

            auto tex_xzzz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 142);

            auto tey_xzzz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 142);

            auto tez_xzzz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 142);

            auto tex_xzzz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 143);

            auto tey_xzzz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 143);

            auto tez_xzzz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 143);

            auto tex_xzzz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 144);

            auto tey_xzzz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 144);

            auto tez_xzzz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 144);

            // Batch of Integrals (387,435)

            #pragma omp simd aligned(fx, pa_x, pc_x, ta_yzz_xzzz_1, ta_yzz_yyyy_1, ta_yzz_yyyz_1, \
                                         ta_yzz_yyzz_1, ta_yzz_yzzz_1, ta_yzz_zzzz_1, ta_zzz_xxxx_1, ta_zzz_xxxy_1, \
                                         ta_zzz_xxxz_1, ta_zzz_xxyy_1, ta_zzz_xxyz_1, ta_zzz_xxzz_1, ta_zzz_xyyy_1, \
                                         ta_zzz_xyyz_1, ta_zzz_xyzz_1, ta_zzz_xzzz_1, tex_xyzz_xzzz_0, tex_xyzz_yyyy_0, \
                                         tex_xyzz_yyyz_0, tex_xyzz_yyzz_0, tex_xyzz_yzzz_0, tex_xyzz_zzzz_0, tex_xzzz_xxxx_0, \
                                         tex_xzzz_xxxy_0, tex_xzzz_xxxz_0, tex_xzzz_xxyy_0, tex_xzzz_xxyz_0, tex_xzzz_xxzz_0, \
                                         tex_xzzz_xyyy_0, tex_xzzz_xyyz_0, tex_xzzz_xyzz_0, tex_xzzz_xzzz_0, tex_yzz_xzzz_0, \
                                         tex_yzz_xzzz_1, tex_yzz_yyyy_0, tex_yzz_yyyy_1, tex_yzz_yyyz_0, tex_yzz_yyyz_1, \
                                         tex_yzz_yyzz_0, tex_yzz_yyzz_1, tex_yzz_yzzz_0, tex_yzz_yzzz_1, tex_yzz_zzz_0, \
                                         tex_yzz_zzz_1, tex_yzz_zzzz_0, tex_yzz_zzzz_1, tex_zzz_xxx_0, tex_zzz_xxx_1, \
                                         tex_zzz_xxxx_0, tex_zzz_xxxx_1, tex_zzz_xxxy_0, tex_zzz_xxxy_1, tex_zzz_xxxz_0, \
                                         tex_zzz_xxxz_1, tex_zzz_xxy_0, tex_zzz_xxy_1, tex_zzz_xxyy_0, tex_zzz_xxyy_1, \
                                         tex_zzz_xxyz_0, tex_zzz_xxyz_1, tex_zzz_xxz_0, tex_zzz_xxz_1, tex_zzz_xxzz_0, \
                                         tex_zzz_xxzz_1, tex_zzz_xyy_0, tex_zzz_xyy_1, tex_zzz_xyyy_0, tex_zzz_xyyy_1, \
                                         tex_zzz_xyyz_0, tex_zzz_xyyz_1, tex_zzz_xyz_0, tex_zzz_xyz_1, tex_zzz_xyzz_0, \
                                         tex_zzz_xyzz_1, tex_zzz_xzz_0, tex_zzz_xzz_1, tex_zzz_xzzz_0, tex_zzz_xzzz_1, \
                                         tex_zzz_yyy_0, tex_zzz_yyy_1, tex_zzz_yyz_0, tex_zzz_yyz_1, tex_zzz_yzz_0, \
                                         tex_zzz_yzz_1, tex_zzz_zzz_0, tex_zzz_zzz_1, tey_xyzz_xzzz_0, tey_xyzz_yyyy_0, \
                                         tey_xyzz_yyyz_0, tey_xyzz_yyzz_0, tey_xyzz_yzzz_0, tey_xyzz_zzzz_0, tey_xzzz_xxxx_0, \
                                         tey_xzzz_xxxy_0, tey_xzzz_xxxz_0, tey_xzzz_xxyy_0, tey_xzzz_xxyz_0, tey_xzzz_xxzz_0, \
                                         tey_xzzz_xyyy_0, tey_xzzz_xyyz_0, tey_xzzz_xyzz_0, tey_xzzz_xzzz_0, tey_yzz_xzzz_0, \
                                         tey_yzz_xzzz_1, tey_yzz_yyyy_0, tey_yzz_yyyy_1, tey_yzz_yyyz_0, tey_yzz_yyyz_1, \
                                         tey_yzz_yyzz_0, tey_yzz_yyzz_1, tey_yzz_yzzz_0, tey_yzz_yzzz_1, tey_yzz_zzz_0, \
                                         tey_yzz_zzz_1, tey_yzz_zzzz_0, tey_yzz_zzzz_1, tey_zzz_xxx_0, tey_zzz_xxx_1, \
                                         tey_zzz_xxxx_0, tey_zzz_xxxx_1, tey_zzz_xxxy_0, tey_zzz_xxxy_1, tey_zzz_xxxz_0, \
                                         tey_zzz_xxxz_1, tey_zzz_xxy_0, tey_zzz_xxy_1, tey_zzz_xxyy_0, tey_zzz_xxyy_1, \
                                         tey_zzz_xxyz_0, tey_zzz_xxyz_1, tey_zzz_xxz_0, tey_zzz_xxz_1, tey_zzz_xxzz_0, \
                                         tey_zzz_xxzz_1, tey_zzz_xyy_0, tey_zzz_xyy_1, tey_zzz_xyyy_0, tey_zzz_xyyy_1, \
                                         tey_zzz_xyyz_0, tey_zzz_xyyz_1, tey_zzz_xyz_0, tey_zzz_xyz_1, tey_zzz_xyzz_0, \
                                         tey_zzz_xyzz_1, tey_zzz_xzz_0, tey_zzz_xzz_1, tey_zzz_xzzz_0, tey_zzz_xzzz_1, \
                                         tey_zzz_yyy_0, tey_zzz_yyy_1, tey_zzz_yyz_0, tey_zzz_yyz_1, tey_zzz_yzz_0, \
                                         tey_zzz_yzz_1, tey_zzz_zzz_0, tey_zzz_zzz_1, tez_xyzz_xzzz_0, tez_xyzz_yyyy_0, \
                                         tez_xyzz_yyyz_0, tez_xyzz_yyzz_0, tez_xyzz_yzzz_0, tez_xyzz_zzzz_0, tez_xzzz_xxxx_0, \
                                         tez_xzzz_xxxy_0, tez_xzzz_xxxz_0, tez_xzzz_xxyy_0, tez_xzzz_xxyz_0, tez_xzzz_xxzz_0, \
                                         tez_xzzz_xyyy_0, tez_xzzz_xyyz_0, tez_xzzz_xyzz_0, tez_xzzz_xzzz_0, tez_yzz_xzzz_0, \
                                         tez_yzz_xzzz_1, tez_yzz_yyyy_0, tez_yzz_yyyy_1, tez_yzz_yyyz_0, tez_yzz_yyyz_1, \
                                         tez_yzz_yyzz_0, tez_yzz_yyzz_1, tez_yzz_yzzz_0, tez_yzz_yzzz_1, tez_yzz_zzz_0, \
                                         tez_yzz_zzz_1, tez_yzz_zzzz_0, tez_yzz_zzzz_1, tez_zzz_xxx_0, tez_zzz_xxx_1, \
                                         tez_zzz_xxxx_0, tez_zzz_xxxx_1, tez_zzz_xxxy_0, tez_zzz_xxxy_1, tez_zzz_xxxz_0, \
                                         tez_zzz_xxxz_1, tez_zzz_xxy_0, tez_zzz_xxy_1, tez_zzz_xxyy_0, tez_zzz_xxyy_1, \
                                         tez_zzz_xxyz_0, tez_zzz_xxyz_1, tez_zzz_xxz_0, tez_zzz_xxz_1, tez_zzz_xxzz_0, \
                                         tez_zzz_xxzz_1, tez_zzz_xyy_0, tez_zzz_xyy_1, tez_zzz_xyyy_0, tez_zzz_xyyy_1, \
                                         tez_zzz_xyyz_0, tez_zzz_xyyz_1, tez_zzz_xyz_0, tez_zzz_xyz_1, tez_zzz_xyzz_0, \
                                         tez_zzz_xyzz_1, tez_zzz_xzz_0, tez_zzz_xzz_1, tez_zzz_xzzz_0, tez_zzz_xzzz_1, \
                                         tez_zzz_yyy_0, tez_zzz_yyy_1, tez_zzz_yyz_0, tez_zzz_yyz_1, tez_zzz_yzz_0, \
                                         tez_zzz_yzz_1, tez_zzz_zzz_0, tez_zzz_zzz_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_xyzz_xzzz_0[j] = pa_x[j] * tex_yzz_xzzz_0[j] - pc_x[j] * tex_yzz_xzzz_1[j] + 0.5 * fl1_fx * tex_yzz_zzz_0[j] -
                                     0.5 * fl1_fx * tex_yzz_zzz_1[j] + ta_yzz_xzzz_1[j];

                tey_xyzz_xzzz_0[j] =
                    pa_x[j] * tey_yzz_xzzz_0[j] - pc_x[j] * tey_yzz_xzzz_1[j] + 0.5 * fl1_fx * tey_yzz_zzz_0[j] - 0.5 * fl1_fx * tey_yzz_zzz_1[j];

                tez_xyzz_xzzz_0[j] =
                    pa_x[j] * tez_yzz_xzzz_0[j] - pc_x[j] * tez_yzz_xzzz_1[j] + 0.5 * fl1_fx * tez_yzz_zzz_0[j] - 0.5 * fl1_fx * tez_yzz_zzz_1[j];

                tex_xyzz_yyyy_0[j] = pa_x[j] * tex_yzz_yyyy_0[j] - pc_x[j] * tex_yzz_yyyy_1[j] + ta_yzz_yyyy_1[j];

                tey_xyzz_yyyy_0[j] = pa_x[j] * tey_yzz_yyyy_0[j] - pc_x[j] * tey_yzz_yyyy_1[j];

                tez_xyzz_yyyy_0[j] = pa_x[j] * tez_yzz_yyyy_0[j] - pc_x[j] * tez_yzz_yyyy_1[j];

                tex_xyzz_yyyz_0[j] = pa_x[j] * tex_yzz_yyyz_0[j] - pc_x[j] * tex_yzz_yyyz_1[j] + ta_yzz_yyyz_1[j];

                tey_xyzz_yyyz_0[j] = pa_x[j] * tey_yzz_yyyz_0[j] - pc_x[j] * tey_yzz_yyyz_1[j];

                tez_xyzz_yyyz_0[j] = pa_x[j] * tez_yzz_yyyz_0[j] - pc_x[j] * tez_yzz_yyyz_1[j];

                tex_xyzz_yyzz_0[j] = pa_x[j] * tex_yzz_yyzz_0[j] - pc_x[j] * tex_yzz_yyzz_1[j] + ta_yzz_yyzz_1[j];

                tey_xyzz_yyzz_0[j] = pa_x[j] * tey_yzz_yyzz_0[j] - pc_x[j] * tey_yzz_yyzz_1[j];

                tez_xyzz_yyzz_0[j] = pa_x[j] * tez_yzz_yyzz_0[j] - pc_x[j] * tez_yzz_yyzz_1[j];

                tex_xyzz_yzzz_0[j] = pa_x[j] * tex_yzz_yzzz_0[j] - pc_x[j] * tex_yzz_yzzz_1[j] + ta_yzz_yzzz_1[j];

                tey_xyzz_yzzz_0[j] = pa_x[j] * tey_yzz_yzzz_0[j] - pc_x[j] * tey_yzz_yzzz_1[j];

                tez_xyzz_yzzz_0[j] = pa_x[j] * tez_yzz_yzzz_0[j] - pc_x[j] * tez_yzz_yzzz_1[j];

                tex_xyzz_zzzz_0[j] = pa_x[j] * tex_yzz_zzzz_0[j] - pc_x[j] * tex_yzz_zzzz_1[j] + ta_yzz_zzzz_1[j];

                tey_xyzz_zzzz_0[j] = pa_x[j] * tey_yzz_zzzz_0[j] - pc_x[j] * tey_yzz_zzzz_1[j];

                tez_xyzz_zzzz_0[j] = pa_x[j] * tez_yzz_zzzz_0[j] - pc_x[j] * tez_yzz_zzzz_1[j];

                tex_xzzz_xxxx_0[j] = pa_x[j] * tex_zzz_xxxx_0[j] - pc_x[j] * tex_zzz_xxxx_1[j] + 2.0 * fl1_fx * tex_zzz_xxx_0[j] -
                                     2.0 * fl1_fx * tex_zzz_xxx_1[j] + ta_zzz_xxxx_1[j];

                tey_xzzz_xxxx_0[j] =
                    pa_x[j] * tey_zzz_xxxx_0[j] - pc_x[j] * tey_zzz_xxxx_1[j] + 2.0 * fl1_fx * tey_zzz_xxx_0[j] - 2.0 * fl1_fx * tey_zzz_xxx_1[j];

                tez_xzzz_xxxx_0[j] =
                    pa_x[j] * tez_zzz_xxxx_0[j] - pc_x[j] * tez_zzz_xxxx_1[j] + 2.0 * fl1_fx * tez_zzz_xxx_0[j] - 2.0 * fl1_fx * tez_zzz_xxx_1[j];

                tex_xzzz_xxxy_0[j] = pa_x[j] * tex_zzz_xxxy_0[j] - pc_x[j] * tex_zzz_xxxy_1[j] + 1.5 * fl1_fx * tex_zzz_xxy_0[j] -
                                     1.5 * fl1_fx * tex_zzz_xxy_1[j] + ta_zzz_xxxy_1[j];

                tey_xzzz_xxxy_0[j] =
                    pa_x[j] * tey_zzz_xxxy_0[j] - pc_x[j] * tey_zzz_xxxy_1[j] + 1.5 * fl1_fx * tey_zzz_xxy_0[j] - 1.5 * fl1_fx * tey_zzz_xxy_1[j];

                tez_xzzz_xxxy_0[j] =
                    pa_x[j] * tez_zzz_xxxy_0[j] - pc_x[j] * tez_zzz_xxxy_1[j] + 1.5 * fl1_fx * tez_zzz_xxy_0[j] - 1.5 * fl1_fx * tez_zzz_xxy_1[j];

                tex_xzzz_xxxz_0[j] = pa_x[j] * tex_zzz_xxxz_0[j] - pc_x[j] * tex_zzz_xxxz_1[j] + 1.5 * fl1_fx * tex_zzz_xxz_0[j] -
                                     1.5 * fl1_fx * tex_zzz_xxz_1[j] + ta_zzz_xxxz_1[j];

                tey_xzzz_xxxz_0[j] =
                    pa_x[j] * tey_zzz_xxxz_0[j] - pc_x[j] * tey_zzz_xxxz_1[j] + 1.5 * fl1_fx * tey_zzz_xxz_0[j] - 1.5 * fl1_fx * tey_zzz_xxz_1[j];

                tez_xzzz_xxxz_0[j] =
                    pa_x[j] * tez_zzz_xxxz_0[j] - pc_x[j] * tez_zzz_xxxz_1[j] + 1.5 * fl1_fx * tez_zzz_xxz_0[j] - 1.5 * fl1_fx * tez_zzz_xxz_1[j];

                tex_xzzz_xxyy_0[j] = pa_x[j] * tex_zzz_xxyy_0[j] - pc_x[j] * tex_zzz_xxyy_1[j] + fl1_fx * tex_zzz_xyy_0[j] -
                                     fl1_fx * tex_zzz_xyy_1[j] + ta_zzz_xxyy_1[j];

                tey_xzzz_xxyy_0[j] =
                    pa_x[j] * tey_zzz_xxyy_0[j] - pc_x[j] * tey_zzz_xxyy_1[j] + fl1_fx * tey_zzz_xyy_0[j] - fl1_fx * tey_zzz_xyy_1[j];

                tez_xzzz_xxyy_0[j] =
                    pa_x[j] * tez_zzz_xxyy_0[j] - pc_x[j] * tez_zzz_xxyy_1[j] + fl1_fx * tez_zzz_xyy_0[j] - fl1_fx * tez_zzz_xyy_1[j];

                tex_xzzz_xxyz_0[j] = pa_x[j] * tex_zzz_xxyz_0[j] - pc_x[j] * tex_zzz_xxyz_1[j] + fl1_fx * tex_zzz_xyz_0[j] -
                                     fl1_fx * tex_zzz_xyz_1[j] + ta_zzz_xxyz_1[j];

                tey_xzzz_xxyz_0[j] =
                    pa_x[j] * tey_zzz_xxyz_0[j] - pc_x[j] * tey_zzz_xxyz_1[j] + fl1_fx * tey_zzz_xyz_0[j] - fl1_fx * tey_zzz_xyz_1[j];

                tez_xzzz_xxyz_0[j] =
                    pa_x[j] * tez_zzz_xxyz_0[j] - pc_x[j] * tez_zzz_xxyz_1[j] + fl1_fx * tez_zzz_xyz_0[j] - fl1_fx * tez_zzz_xyz_1[j];

                tex_xzzz_xxzz_0[j] = pa_x[j] * tex_zzz_xxzz_0[j] - pc_x[j] * tex_zzz_xxzz_1[j] + fl1_fx * tex_zzz_xzz_0[j] -
                                     fl1_fx * tex_zzz_xzz_1[j] + ta_zzz_xxzz_1[j];

                tey_xzzz_xxzz_0[j] =
                    pa_x[j] * tey_zzz_xxzz_0[j] - pc_x[j] * tey_zzz_xxzz_1[j] + fl1_fx * tey_zzz_xzz_0[j] - fl1_fx * tey_zzz_xzz_1[j];

                tez_xzzz_xxzz_0[j] =
                    pa_x[j] * tez_zzz_xxzz_0[j] - pc_x[j] * tez_zzz_xxzz_1[j] + fl1_fx * tez_zzz_xzz_0[j] - fl1_fx * tez_zzz_xzz_1[j];

                tex_xzzz_xyyy_0[j] = pa_x[j] * tex_zzz_xyyy_0[j] - pc_x[j] * tex_zzz_xyyy_1[j] + 0.5 * fl1_fx * tex_zzz_yyy_0[j] -
                                     0.5 * fl1_fx * tex_zzz_yyy_1[j] + ta_zzz_xyyy_1[j];

                tey_xzzz_xyyy_0[j] =
                    pa_x[j] * tey_zzz_xyyy_0[j] - pc_x[j] * tey_zzz_xyyy_1[j] + 0.5 * fl1_fx * tey_zzz_yyy_0[j] - 0.5 * fl1_fx * tey_zzz_yyy_1[j];

                tez_xzzz_xyyy_0[j] =
                    pa_x[j] * tez_zzz_xyyy_0[j] - pc_x[j] * tez_zzz_xyyy_1[j] + 0.5 * fl1_fx * tez_zzz_yyy_0[j] - 0.5 * fl1_fx * tez_zzz_yyy_1[j];

                tex_xzzz_xyyz_0[j] = pa_x[j] * tex_zzz_xyyz_0[j] - pc_x[j] * tex_zzz_xyyz_1[j] + 0.5 * fl1_fx * tex_zzz_yyz_0[j] -
                                     0.5 * fl1_fx * tex_zzz_yyz_1[j] + ta_zzz_xyyz_1[j];

                tey_xzzz_xyyz_0[j] =
                    pa_x[j] * tey_zzz_xyyz_0[j] - pc_x[j] * tey_zzz_xyyz_1[j] + 0.5 * fl1_fx * tey_zzz_yyz_0[j] - 0.5 * fl1_fx * tey_zzz_yyz_1[j];

                tez_xzzz_xyyz_0[j] =
                    pa_x[j] * tez_zzz_xyyz_0[j] - pc_x[j] * tez_zzz_xyyz_1[j] + 0.5 * fl1_fx * tez_zzz_yyz_0[j] - 0.5 * fl1_fx * tez_zzz_yyz_1[j];

                tex_xzzz_xyzz_0[j] = pa_x[j] * tex_zzz_xyzz_0[j] - pc_x[j] * tex_zzz_xyzz_1[j] + 0.5 * fl1_fx * tex_zzz_yzz_0[j] -
                                     0.5 * fl1_fx * tex_zzz_yzz_1[j] + ta_zzz_xyzz_1[j];

                tey_xzzz_xyzz_0[j] =
                    pa_x[j] * tey_zzz_xyzz_0[j] - pc_x[j] * tey_zzz_xyzz_1[j] + 0.5 * fl1_fx * tey_zzz_yzz_0[j] - 0.5 * fl1_fx * tey_zzz_yzz_1[j];

                tez_xzzz_xyzz_0[j] =
                    pa_x[j] * tez_zzz_xyzz_0[j] - pc_x[j] * tez_zzz_xyzz_1[j] + 0.5 * fl1_fx * tez_zzz_yzz_0[j] - 0.5 * fl1_fx * tez_zzz_yzz_1[j];

                tex_xzzz_xzzz_0[j] = pa_x[j] * tex_zzz_xzzz_0[j] - pc_x[j] * tex_zzz_xzzz_1[j] + 0.5 * fl1_fx * tex_zzz_zzz_0[j] -
                                     0.5 * fl1_fx * tex_zzz_zzz_1[j] + ta_zzz_xzzz_1[j];

                tey_xzzz_xzzz_0[j] =
                    pa_x[j] * tey_zzz_xzzz_0[j] - pc_x[j] * tey_zzz_xzzz_1[j] + 0.5 * fl1_fx * tey_zzz_zzz_0[j] - 0.5 * fl1_fx * tey_zzz_zzz_1[j];

                tez_xzzz_xzzz_0[j] =
                    pa_x[j] * tez_zzz_xzzz_0[j] - pc_x[j] * tez_zzz_xzzz_1[j] + 0.5 * fl1_fx * tez_zzz_zzz_0[j] - 0.5 * fl1_fx * tez_zzz_zzz_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForGG_435_483(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pcDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto)
{
    // Batch of Integrals (435,483)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_4_4_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_2_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tex_yyy_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 90);

            auto tey_yyy_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 90);

            auto tez_yyy_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 90);

            auto tex_yyy_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 91);

            auto tey_yyy_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 91);

            auto tez_yyy_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 91);

            auto tex_yyy_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 92);

            auto tey_yyy_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 92);

            auto tez_yyy_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 92);

            auto tex_yyy_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 93);

            auto tey_yyy_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 93);

            auto tez_yyy_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 93);

            auto tex_yyy_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 94);

            auto tey_yyy_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 94);

            auto tez_yyy_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 94);

            auto tex_yyy_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 95);

            auto tey_yyy_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 95);

            auto tez_yyy_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 95);

            auto tex_yyy_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 96);

            auto tey_yyy_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 96);

            auto tez_yyy_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 96);

            auto tex_yyy_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 97);

            auto tey_yyy_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 97);

            auto tez_yyy_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 97);

            auto tex_yyy_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 98);

            auto tey_yyy_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 98);

            auto tez_yyy_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 98);

            auto tex_yyy_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 99);

            auto tey_yyy_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 99);

            auto tez_yyy_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 99);

            auto tex_yyy_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 100);

            auto tey_yyy_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 100);

            auto tez_yyy_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 100);

            auto tex_zzz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 145);

            auto tey_zzz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 145);

            auto tez_zzz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 145);

            auto tex_zzz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 146);

            auto tey_zzz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 146);

            auto tez_zzz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 146);

            auto tex_zzz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 147);

            auto tey_zzz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 147);

            auto tez_zzz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 147);

            auto tex_zzz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 148);

            auto tey_zzz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 148);

            auto tez_zzz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 148);

            auto tex_zzz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 149);

            auto tey_zzz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 149);

            auto tez_zzz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 149);

            auto tex_yyy_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 90);

            auto tey_yyy_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 90);

            auto tez_yyy_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 90);

            auto tex_yyy_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 91);

            auto tey_yyy_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 91);

            auto tez_yyy_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 91);

            auto tex_yyy_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 92);

            auto tey_yyy_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 92);

            auto tez_yyy_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 92);

            auto tex_yyy_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 93);

            auto tey_yyy_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 93);

            auto tez_yyy_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 93);

            auto tex_yyy_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 94);

            auto tey_yyy_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 94);

            auto tez_yyy_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 94);

            auto tex_yyy_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 95);

            auto tey_yyy_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 95);

            auto tez_yyy_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 95);

            auto tex_yyy_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 96);

            auto tey_yyy_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 96);

            auto tez_yyy_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 96);

            auto tex_yyy_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 97);

            auto tey_yyy_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 97);

            auto tez_yyy_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 97);

            auto tex_yyy_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 98);

            auto tey_yyy_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 98);

            auto tez_yyy_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 98);

            auto tex_yyy_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 99);

            auto tey_yyy_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 99);

            auto tez_yyy_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 99);

            auto tex_yyy_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 100);

            auto tey_yyy_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 100);

            auto tez_yyy_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 100);

            auto tex_zzz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 145);

            auto tey_zzz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 145);

            auto tez_zzz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 145);

            auto tex_zzz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 146);

            auto tey_zzz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 146);

            auto tez_zzz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 146);

            auto tex_zzz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 147);

            auto tey_zzz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 147);

            auto tez_zzz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 147);

            auto tex_zzz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 148);

            auto tey_zzz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 148);

            auto tez_zzz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 148);

            auto tex_zzz_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 149);

            auto tey_zzz_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 149);

            auto tez_zzz_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 149);

            auto tex_yy_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 45);

            auto tey_yy_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 45);

            auto tez_yy_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 45);

            auto tex_yy_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 46);

            auto tey_yy_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 46);

            auto tez_yy_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 46);

            auto tex_yy_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 47);

            auto tey_yy_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 47);

            auto tez_yy_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 47);

            auto tex_yy_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 48);

            auto tey_yy_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 48);

            auto tez_yy_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 48);

            auto tex_yy_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 49);

            auto tey_yy_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 49);

            auto tez_yy_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 49);

            auto tex_yy_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 50);

            auto tey_yy_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 50);

            auto tez_yy_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 50);

            auto tex_yy_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 51);

            auto tey_yy_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 51);

            auto tez_yy_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 51);

            auto tex_yy_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 52);

            auto tey_yy_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 52);

            auto tez_yy_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 52);

            auto tex_yy_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 53);

            auto tey_yy_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 53);

            auto tez_yy_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 53);

            auto tex_yy_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 54);

            auto tey_yy_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 54);

            auto tez_yy_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 54);

            auto tex_yy_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 55);

            auto tey_yy_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 55);

            auto tez_yy_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 55);

            auto tex_yy_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 45);

            auto tey_yy_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 45);

            auto tez_yy_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 45);

            auto tex_yy_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 46);

            auto tey_yy_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 46);

            auto tez_yy_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 46);

            auto tex_yy_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 47);

            auto tey_yy_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 47);

            auto tez_yy_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 47);

            auto tex_yy_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 48);

            auto tey_yy_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 48);

            auto tez_yy_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 48);

            auto tex_yy_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 49);

            auto tey_yy_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 49);

            auto tez_yy_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 49);

            auto tex_yy_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 50);

            auto tey_yy_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 50);

            auto tez_yy_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 50);

            auto tex_yy_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 51);

            auto tey_yy_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 51);

            auto tez_yy_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 51);

            auto tex_yy_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 52);

            auto tey_yy_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 52);

            auto tez_yy_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 52);

            auto tex_yy_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 53);

            auto tey_yy_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 53);

            auto tez_yy_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 53);

            auto tex_yy_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 54);

            auto tey_yy_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 54);

            auto tez_yy_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 54);

            auto tex_yy_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 55);

            auto tey_yy_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 55);

            auto tez_yy_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 55);

            auto tex_yyy_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 60);

            auto tey_yyy_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 60);

            auto tez_yyy_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 60);

            auto tex_yyy_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 61);

            auto tey_yyy_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 61);

            auto tez_yyy_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 61);

            auto tex_yyy_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 62);

            auto tey_yyy_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 62);

            auto tez_yyy_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 62);

            auto tex_yyy_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 63);

            auto tey_yyy_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 63);

            auto tez_yyy_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 63);

            auto tex_yyy_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 64);

            auto tey_yyy_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 64);

            auto tez_yyy_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 64);

            auto tex_yyy_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 65);

            auto tey_yyy_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 65);

            auto tez_yyy_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 65);

            auto tex_yyy_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 66);

            auto tey_yyy_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 66);

            auto tez_yyy_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 66);

            auto tex_yyy_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 60);

            auto tey_yyy_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 60);

            auto tez_yyy_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 60);

            auto tex_yyy_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 61);

            auto tey_yyy_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 61);

            auto tez_yyy_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 61);

            auto tex_yyy_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 62);

            auto tey_yyy_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 62);

            auto tez_yyy_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 62);

            auto tex_yyy_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 63);

            auto tey_yyy_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 63);

            auto tez_yyy_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 63);

            auto tex_yyy_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 64);

            auto tey_yyy_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 64);

            auto tez_yyy_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 64);

            auto tex_yyy_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 65);

            auto tey_yyy_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 65);

            auto tez_yyy_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 65);

            auto tex_yyy_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 66);

            auto tey_yyy_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 66);

            auto tez_yyy_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 66);

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

            auto ta_zzz_yyyy_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 145);

            auto ta_zzz_yyyz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 146);

            auto ta_zzz_yyzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 147);

            auto ta_zzz_yzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 148);

            auto ta_zzz_zzzz_1 = primBuffer.data(pidx_a_3_4_m1 + 150 * idx + 149);

            // set up pointers to integrals

            auto tex_xzzz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 145);

            auto tey_xzzz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 145);

            auto tez_xzzz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 145);

            auto tex_xzzz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 146);

            auto tey_xzzz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 146);

            auto tez_xzzz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 146);

            auto tex_xzzz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 147);

            auto tey_xzzz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 147);

            auto tez_xzzz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 147);

            auto tex_xzzz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 148);

            auto tey_xzzz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 148);

            auto tez_xzzz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 148);

            auto tex_xzzz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 149);

            auto tey_xzzz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 149);

            auto tez_xzzz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 149);

            auto tex_yyyy_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 150);

            auto tey_yyyy_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 150);

            auto tez_yyyy_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 150);

            auto tex_yyyy_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 151);

            auto tey_yyyy_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 151);

            auto tez_yyyy_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 151);

            auto tex_yyyy_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 152);

            auto tey_yyyy_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 152);

            auto tez_yyyy_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 152);

            auto tex_yyyy_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 153);

            auto tey_yyyy_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 153);

            auto tez_yyyy_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 153);

            auto tex_yyyy_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 154);

            auto tey_yyyy_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 154);

            auto tez_yyyy_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 154);

            auto tex_yyyy_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 155);

            auto tey_yyyy_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 155);

            auto tez_yyyy_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 155);

            auto tex_yyyy_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 156);

            auto tey_yyyy_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 156);

            auto tez_yyyy_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 156);

            auto tex_yyyy_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 157);

            auto tey_yyyy_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 157);

            auto tez_yyyy_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 157);

            auto tex_yyyy_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 158);

            auto tey_yyyy_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 158);

            auto tez_yyyy_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 158);

            auto tex_yyyy_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 159);

            auto tey_yyyy_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 159);

            auto tez_yyyy_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 159);

            auto tex_yyyy_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 160);

            auto tey_yyyy_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 160);

            auto tez_yyyy_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 160);

            // Batch of Integrals (435,483)

            #pragma omp simd aligned(fx, pa_x, pa_y, pc_x, pc_y, ta_yyy_xxxx_1, ta_yyy_xxxy_1, ta_yyy_xxxz_1, \
                                         ta_yyy_xxyy_1, ta_yyy_xxyz_1, ta_yyy_xxzz_1, ta_yyy_xyyy_1, ta_yyy_xyyz_1, \
                                         ta_yyy_xyzz_1, ta_yyy_xzzz_1, ta_yyy_yyyy_1, ta_zzz_yyyy_1, ta_zzz_yyyz_1, \
                                         ta_zzz_yyzz_1, ta_zzz_yzzz_1, ta_zzz_zzzz_1, tex_xzzz_yyyy_0, tex_xzzz_yyyz_0, \
                                         tex_xzzz_yyzz_0, tex_xzzz_yzzz_0, tex_xzzz_zzzz_0, tex_yy_xxxx_0, tex_yy_xxxx_1, \
                                         tex_yy_xxxy_0, tex_yy_xxxy_1, tex_yy_xxxz_0, tex_yy_xxxz_1, tex_yy_xxyy_0, \
                                         tex_yy_xxyy_1, tex_yy_xxyz_0, tex_yy_xxyz_1, tex_yy_xxzz_0, tex_yy_xxzz_1, \
                                         tex_yy_xyyy_0, tex_yy_xyyy_1, tex_yy_xyyz_0, tex_yy_xyyz_1, tex_yy_xyzz_0, \
                                         tex_yy_xyzz_1, tex_yy_xzzz_0, tex_yy_xzzz_1, tex_yy_yyyy_0, tex_yy_yyyy_1, \
                                         tex_yyy_xxx_0, tex_yyy_xxx_1, tex_yyy_xxxx_0, tex_yyy_xxxx_1, tex_yyy_xxxy_0, \
                                         tex_yyy_xxxy_1, tex_yyy_xxxz_0, tex_yyy_xxxz_1, tex_yyy_xxy_0, tex_yyy_xxy_1, \
                                         tex_yyy_xxyy_0, tex_yyy_xxyy_1, tex_yyy_xxyz_0, tex_yyy_xxyz_1, tex_yyy_xxz_0, \
                                         tex_yyy_xxz_1, tex_yyy_xxzz_0, tex_yyy_xxzz_1, tex_yyy_xyy_0, tex_yyy_xyy_1, \
                                         tex_yyy_xyyy_0, tex_yyy_xyyy_1, tex_yyy_xyyz_0, tex_yyy_xyyz_1, tex_yyy_xyz_0, \
                                         tex_yyy_xyz_1, tex_yyy_xyzz_0, tex_yyy_xyzz_1, tex_yyy_xzz_0, tex_yyy_xzz_1, \
                                         tex_yyy_xzzz_0, tex_yyy_xzzz_1, tex_yyy_yyy_0, tex_yyy_yyy_1, tex_yyy_yyyy_0, \
                                         tex_yyy_yyyy_1, tex_yyyy_xxxx_0, tex_yyyy_xxxy_0, tex_yyyy_xxxz_0, tex_yyyy_xxyy_0, \
                                         tex_yyyy_xxyz_0, tex_yyyy_xxzz_0, tex_yyyy_xyyy_0, tex_yyyy_xyyz_0, tex_yyyy_xyzz_0, \
                                         tex_yyyy_xzzz_0, tex_yyyy_yyyy_0, tex_zzz_yyyy_0, tex_zzz_yyyy_1, tex_zzz_yyyz_0, \
                                         tex_zzz_yyyz_1, tex_zzz_yyzz_0, tex_zzz_yyzz_1, tex_zzz_yzzz_0, tex_zzz_yzzz_1, \
                                         tex_zzz_zzzz_0, tex_zzz_zzzz_1, tey_xzzz_yyyy_0, tey_xzzz_yyyz_0, tey_xzzz_yyzz_0, \
                                         tey_xzzz_yzzz_0, tey_xzzz_zzzz_0, tey_yy_xxxx_0, tey_yy_xxxx_1, tey_yy_xxxy_0, \
                                         tey_yy_xxxy_1, tey_yy_xxxz_0, tey_yy_xxxz_1, tey_yy_xxyy_0, tey_yy_xxyy_1, \
                                         tey_yy_xxyz_0, tey_yy_xxyz_1, tey_yy_xxzz_0, tey_yy_xxzz_1, tey_yy_xyyy_0, \
                                         tey_yy_xyyy_1, tey_yy_xyyz_0, tey_yy_xyyz_1, tey_yy_xyzz_0, tey_yy_xyzz_1, \
                                         tey_yy_xzzz_0, tey_yy_xzzz_1, tey_yy_yyyy_0, tey_yy_yyyy_1, tey_yyy_xxx_0, \
                                         tey_yyy_xxx_1, tey_yyy_xxxx_0, tey_yyy_xxxx_1, tey_yyy_xxxy_0, tey_yyy_xxxy_1, \
                                         tey_yyy_xxxz_0, tey_yyy_xxxz_1, tey_yyy_xxy_0, tey_yyy_xxy_1, tey_yyy_xxyy_0, \
                                         tey_yyy_xxyy_1, tey_yyy_xxyz_0, tey_yyy_xxyz_1, tey_yyy_xxz_0, tey_yyy_xxz_1, \
                                         tey_yyy_xxzz_0, tey_yyy_xxzz_1, tey_yyy_xyy_0, tey_yyy_xyy_1, tey_yyy_xyyy_0, \
                                         tey_yyy_xyyy_1, tey_yyy_xyyz_0, tey_yyy_xyyz_1, tey_yyy_xyz_0, tey_yyy_xyz_1, \
                                         tey_yyy_xyzz_0, tey_yyy_xyzz_1, tey_yyy_xzz_0, tey_yyy_xzz_1, tey_yyy_xzzz_0, \
                                         tey_yyy_xzzz_1, tey_yyy_yyy_0, tey_yyy_yyy_1, tey_yyy_yyyy_0, tey_yyy_yyyy_1, \
                                         tey_yyyy_xxxx_0, tey_yyyy_xxxy_0, tey_yyyy_xxxz_0, tey_yyyy_xxyy_0, tey_yyyy_xxyz_0, \
                                         tey_yyyy_xxzz_0, tey_yyyy_xyyy_0, tey_yyyy_xyyz_0, tey_yyyy_xyzz_0, tey_yyyy_xzzz_0, \
                                         tey_yyyy_yyyy_0, tey_zzz_yyyy_0, tey_zzz_yyyy_1, tey_zzz_yyyz_0, tey_zzz_yyyz_1, \
                                         tey_zzz_yyzz_0, tey_zzz_yyzz_1, tey_zzz_yzzz_0, tey_zzz_yzzz_1, tey_zzz_zzzz_0, \
                                         tey_zzz_zzzz_1, tez_xzzz_yyyy_0, tez_xzzz_yyyz_0, tez_xzzz_yyzz_0, tez_xzzz_yzzz_0, \
                                         tez_xzzz_zzzz_0, tez_yy_xxxx_0, tez_yy_xxxx_1, tez_yy_xxxy_0, tez_yy_xxxy_1, \
                                         tez_yy_xxxz_0, tez_yy_xxxz_1, tez_yy_xxyy_0, tez_yy_xxyy_1, tez_yy_xxyz_0, \
                                         tez_yy_xxyz_1, tez_yy_xxzz_0, tez_yy_xxzz_1, tez_yy_xyyy_0, tez_yy_xyyy_1, \
                                         tez_yy_xyyz_0, tez_yy_xyyz_1, tez_yy_xyzz_0, tez_yy_xyzz_1, tez_yy_xzzz_0, \
                                         tez_yy_xzzz_1, tez_yy_yyyy_0, tez_yy_yyyy_1, tez_yyy_xxx_0, tez_yyy_xxx_1, \
                                         tez_yyy_xxxx_0, tez_yyy_xxxx_1, tez_yyy_xxxy_0, tez_yyy_xxxy_1, tez_yyy_xxxz_0, \
                                         tez_yyy_xxxz_1, tez_yyy_xxy_0, tez_yyy_xxy_1, tez_yyy_xxyy_0, tez_yyy_xxyy_1, \
                                         tez_yyy_xxyz_0, tez_yyy_xxyz_1, tez_yyy_xxz_0, tez_yyy_xxz_1, tez_yyy_xxzz_0, \
                                         tez_yyy_xxzz_1, tez_yyy_xyy_0, tez_yyy_xyy_1, tez_yyy_xyyy_0, tez_yyy_xyyy_1, \
                                         tez_yyy_xyyz_0, tez_yyy_xyyz_1, tez_yyy_xyz_0, tez_yyy_xyz_1, tez_yyy_xyzz_0, \
                                         tez_yyy_xyzz_1, tez_yyy_xzz_0, tez_yyy_xzz_1, tez_yyy_xzzz_0, tez_yyy_xzzz_1, \
                                         tez_yyy_yyy_0, tez_yyy_yyy_1, tez_yyy_yyyy_0, tez_yyy_yyyy_1, tez_yyyy_xxxx_0, \
                                         tez_yyyy_xxxy_0, tez_yyyy_xxxz_0, tez_yyyy_xxyy_0, tez_yyyy_xxyz_0, tez_yyyy_xxzz_0, \
                                         tez_yyyy_xyyy_0, tez_yyyy_xyyz_0, tez_yyyy_xyzz_0, tez_yyyy_xzzz_0, tez_yyyy_yyyy_0, \
                                         tez_zzz_yyyy_0, tez_zzz_yyyy_1, tez_zzz_yyyz_0, tez_zzz_yyyz_1, tez_zzz_yyzz_0, \
                                         tez_zzz_yyzz_1, tez_zzz_yzzz_0, tez_zzz_yzzz_1, tez_zzz_zzzz_0, tez_zzz_zzzz_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_xzzz_yyyy_0[j] = pa_x[j] * tex_zzz_yyyy_0[j] - pc_x[j] * tex_zzz_yyyy_1[j] + ta_zzz_yyyy_1[j];

                tey_xzzz_yyyy_0[j] = pa_x[j] * tey_zzz_yyyy_0[j] - pc_x[j] * tey_zzz_yyyy_1[j];

                tez_xzzz_yyyy_0[j] = pa_x[j] * tez_zzz_yyyy_0[j] - pc_x[j] * tez_zzz_yyyy_1[j];

                tex_xzzz_yyyz_0[j] = pa_x[j] * tex_zzz_yyyz_0[j] - pc_x[j] * tex_zzz_yyyz_1[j] + ta_zzz_yyyz_1[j];

                tey_xzzz_yyyz_0[j] = pa_x[j] * tey_zzz_yyyz_0[j] - pc_x[j] * tey_zzz_yyyz_1[j];

                tez_xzzz_yyyz_0[j] = pa_x[j] * tez_zzz_yyyz_0[j] - pc_x[j] * tez_zzz_yyyz_1[j];

                tex_xzzz_yyzz_0[j] = pa_x[j] * tex_zzz_yyzz_0[j] - pc_x[j] * tex_zzz_yyzz_1[j] + ta_zzz_yyzz_1[j];

                tey_xzzz_yyzz_0[j] = pa_x[j] * tey_zzz_yyzz_0[j] - pc_x[j] * tey_zzz_yyzz_1[j];

                tez_xzzz_yyzz_0[j] = pa_x[j] * tez_zzz_yyzz_0[j] - pc_x[j] * tez_zzz_yyzz_1[j];

                tex_xzzz_yzzz_0[j] = pa_x[j] * tex_zzz_yzzz_0[j] - pc_x[j] * tex_zzz_yzzz_1[j] + ta_zzz_yzzz_1[j];

                tey_xzzz_yzzz_0[j] = pa_x[j] * tey_zzz_yzzz_0[j] - pc_x[j] * tey_zzz_yzzz_1[j];

                tez_xzzz_yzzz_0[j] = pa_x[j] * tez_zzz_yzzz_0[j] - pc_x[j] * tez_zzz_yzzz_1[j];

                tex_xzzz_zzzz_0[j] = pa_x[j] * tex_zzz_zzzz_0[j] - pc_x[j] * tex_zzz_zzzz_1[j] + ta_zzz_zzzz_1[j];

                tey_xzzz_zzzz_0[j] = pa_x[j] * tey_zzz_zzzz_0[j] - pc_x[j] * tey_zzz_zzzz_1[j];

                tez_xzzz_zzzz_0[j] = pa_x[j] * tez_zzz_zzzz_0[j] - pc_x[j] * tez_zzz_zzzz_1[j];

                tex_yyyy_xxxx_0[j] =
                    pa_y[j] * tex_yyy_xxxx_0[j] - pc_y[j] * tex_yyy_xxxx_1[j] + 1.5 * fl1_fx * tex_yy_xxxx_0[j] - 1.5 * fl1_fx * tex_yy_xxxx_1[j];

                tey_yyyy_xxxx_0[j] = pa_y[j] * tey_yyy_xxxx_0[j] - pc_y[j] * tey_yyy_xxxx_1[j] + 1.5 * fl1_fx * tey_yy_xxxx_0[j] -
                                     1.5 * fl1_fx * tey_yy_xxxx_1[j] + ta_yyy_xxxx_1[j];

                tez_yyyy_xxxx_0[j] =
                    pa_y[j] * tez_yyy_xxxx_0[j] - pc_y[j] * tez_yyy_xxxx_1[j] + 1.5 * fl1_fx * tez_yy_xxxx_0[j] - 1.5 * fl1_fx * tez_yy_xxxx_1[j];

                tex_yyyy_xxxy_0[j] = pa_y[j] * tex_yyy_xxxy_0[j] - pc_y[j] * tex_yyy_xxxy_1[j] + 1.5 * fl1_fx * tex_yy_xxxy_0[j] -
                                     1.5 * fl1_fx * tex_yy_xxxy_1[j] + 0.5 * fl1_fx * tex_yyy_xxx_0[j] - 0.5 * fl1_fx * tex_yyy_xxx_1[j];

                tey_yyyy_xxxy_0[j] = pa_y[j] * tey_yyy_xxxy_0[j] - pc_y[j] * tey_yyy_xxxy_1[j] + 1.5 * fl1_fx * tey_yy_xxxy_0[j] -
                                     1.5 * fl1_fx * tey_yy_xxxy_1[j] + 0.5 * fl1_fx * tey_yyy_xxx_0[j] - 0.5 * fl1_fx * tey_yyy_xxx_1[j] +
                                     ta_yyy_xxxy_1[j];

                tez_yyyy_xxxy_0[j] = pa_y[j] * tez_yyy_xxxy_0[j] - pc_y[j] * tez_yyy_xxxy_1[j] + 1.5 * fl1_fx * tez_yy_xxxy_0[j] -
                                     1.5 * fl1_fx * tez_yy_xxxy_1[j] + 0.5 * fl1_fx * tez_yyy_xxx_0[j] - 0.5 * fl1_fx * tez_yyy_xxx_1[j];

                tex_yyyy_xxxz_0[j] =
                    pa_y[j] * tex_yyy_xxxz_0[j] - pc_y[j] * tex_yyy_xxxz_1[j] + 1.5 * fl1_fx * tex_yy_xxxz_0[j] - 1.5 * fl1_fx * tex_yy_xxxz_1[j];

                tey_yyyy_xxxz_0[j] = pa_y[j] * tey_yyy_xxxz_0[j] - pc_y[j] * tey_yyy_xxxz_1[j] + 1.5 * fl1_fx * tey_yy_xxxz_0[j] -
                                     1.5 * fl1_fx * tey_yy_xxxz_1[j] + ta_yyy_xxxz_1[j];

                tez_yyyy_xxxz_0[j] =
                    pa_y[j] * tez_yyy_xxxz_0[j] - pc_y[j] * tez_yyy_xxxz_1[j] + 1.5 * fl1_fx * tez_yy_xxxz_0[j] - 1.5 * fl1_fx * tez_yy_xxxz_1[j];

                tex_yyyy_xxyy_0[j] = pa_y[j] * tex_yyy_xxyy_0[j] - pc_y[j] * tex_yyy_xxyy_1[j] + 1.5 * fl1_fx * tex_yy_xxyy_0[j] -
                                     1.5 * fl1_fx * tex_yy_xxyy_1[j] + fl1_fx * tex_yyy_xxy_0[j] - fl1_fx * tex_yyy_xxy_1[j];

                tey_yyyy_xxyy_0[j] = pa_y[j] * tey_yyy_xxyy_0[j] - pc_y[j] * tey_yyy_xxyy_1[j] + 1.5 * fl1_fx * tey_yy_xxyy_0[j] -
                                     1.5 * fl1_fx * tey_yy_xxyy_1[j] + fl1_fx * tey_yyy_xxy_0[j] - fl1_fx * tey_yyy_xxy_1[j] + ta_yyy_xxyy_1[j];

                tez_yyyy_xxyy_0[j] = pa_y[j] * tez_yyy_xxyy_0[j] - pc_y[j] * tez_yyy_xxyy_1[j] + 1.5 * fl1_fx * tez_yy_xxyy_0[j] -
                                     1.5 * fl1_fx * tez_yy_xxyy_1[j] + fl1_fx * tez_yyy_xxy_0[j] - fl1_fx * tez_yyy_xxy_1[j];

                tex_yyyy_xxyz_0[j] = pa_y[j] * tex_yyy_xxyz_0[j] - pc_y[j] * tex_yyy_xxyz_1[j] + 1.5 * fl1_fx * tex_yy_xxyz_0[j] -
                                     1.5 * fl1_fx * tex_yy_xxyz_1[j] + 0.5 * fl1_fx * tex_yyy_xxz_0[j] - 0.5 * fl1_fx * tex_yyy_xxz_1[j];

                tey_yyyy_xxyz_0[j] = pa_y[j] * tey_yyy_xxyz_0[j] - pc_y[j] * tey_yyy_xxyz_1[j] + 1.5 * fl1_fx * tey_yy_xxyz_0[j] -
                                     1.5 * fl1_fx * tey_yy_xxyz_1[j] + 0.5 * fl1_fx * tey_yyy_xxz_0[j] - 0.5 * fl1_fx * tey_yyy_xxz_1[j] +
                                     ta_yyy_xxyz_1[j];

                tez_yyyy_xxyz_0[j] = pa_y[j] * tez_yyy_xxyz_0[j] - pc_y[j] * tez_yyy_xxyz_1[j] + 1.5 * fl1_fx * tez_yy_xxyz_0[j] -
                                     1.5 * fl1_fx * tez_yy_xxyz_1[j] + 0.5 * fl1_fx * tez_yyy_xxz_0[j] - 0.5 * fl1_fx * tez_yyy_xxz_1[j];

                tex_yyyy_xxzz_0[j] =
                    pa_y[j] * tex_yyy_xxzz_0[j] - pc_y[j] * tex_yyy_xxzz_1[j] + 1.5 * fl1_fx * tex_yy_xxzz_0[j] - 1.5 * fl1_fx * tex_yy_xxzz_1[j];

                tey_yyyy_xxzz_0[j] = pa_y[j] * tey_yyy_xxzz_0[j] - pc_y[j] * tey_yyy_xxzz_1[j] + 1.5 * fl1_fx * tey_yy_xxzz_0[j] -
                                     1.5 * fl1_fx * tey_yy_xxzz_1[j] + ta_yyy_xxzz_1[j];

                tez_yyyy_xxzz_0[j] =
                    pa_y[j] * tez_yyy_xxzz_0[j] - pc_y[j] * tez_yyy_xxzz_1[j] + 1.5 * fl1_fx * tez_yy_xxzz_0[j] - 1.5 * fl1_fx * tez_yy_xxzz_1[j];

                tex_yyyy_xyyy_0[j] = pa_y[j] * tex_yyy_xyyy_0[j] - pc_y[j] * tex_yyy_xyyy_1[j] + 1.5 * fl1_fx * tex_yy_xyyy_0[j] -
                                     1.5 * fl1_fx * tex_yy_xyyy_1[j] + 1.5 * fl1_fx * tex_yyy_xyy_0[j] - 1.5 * fl1_fx * tex_yyy_xyy_1[j];

                tey_yyyy_xyyy_0[j] = pa_y[j] * tey_yyy_xyyy_0[j] - pc_y[j] * tey_yyy_xyyy_1[j] + 1.5 * fl1_fx * tey_yy_xyyy_0[j] -
                                     1.5 * fl1_fx * tey_yy_xyyy_1[j] + 1.5 * fl1_fx * tey_yyy_xyy_0[j] - 1.5 * fl1_fx * tey_yyy_xyy_1[j] +
                                     ta_yyy_xyyy_1[j];

                tez_yyyy_xyyy_0[j] = pa_y[j] * tez_yyy_xyyy_0[j] - pc_y[j] * tez_yyy_xyyy_1[j] + 1.5 * fl1_fx * tez_yy_xyyy_0[j] -
                                     1.5 * fl1_fx * tez_yy_xyyy_1[j] + 1.5 * fl1_fx * tez_yyy_xyy_0[j] - 1.5 * fl1_fx * tez_yyy_xyy_1[j];

                tex_yyyy_xyyz_0[j] = pa_y[j] * tex_yyy_xyyz_0[j] - pc_y[j] * tex_yyy_xyyz_1[j] + 1.5 * fl1_fx * tex_yy_xyyz_0[j] -
                                     1.5 * fl1_fx * tex_yy_xyyz_1[j] + fl1_fx * tex_yyy_xyz_0[j] - fl1_fx * tex_yyy_xyz_1[j];

                tey_yyyy_xyyz_0[j] = pa_y[j] * tey_yyy_xyyz_0[j] - pc_y[j] * tey_yyy_xyyz_1[j] + 1.5 * fl1_fx * tey_yy_xyyz_0[j] -
                                     1.5 * fl1_fx * tey_yy_xyyz_1[j] + fl1_fx * tey_yyy_xyz_0[j] - fl1_fx * tey_yyy_xyz_1[j] + ta_yyy_xyyz_1[j];

                tez_yyyy_xyyz_0[j] = pa_y[j] * tez_yyy_xyyz_0[j] - pc_y[j] * tez_yyy_xyyz_1[j] + 1.5 * fl1_fx * tez_yy_xyyz_0[j] -
                                     1.5 * fl1_fx * tez_yy_xyyz_1[j] + fl1_fx * tez_yyy_xyz_0[j] - fl1_fx * tez_yyy_xyz_1[j];

                tex_yyyy_xyzz_0[j] = pa_y[j] * tex_yyy_xyzz_0[j] - pc_y[j] * tex_yyy_xyzz_1[j] + 1.5 * fl1_fx * tex_yy_xyzz_0[j] -
                                     1.5 * fl1_fx * tex_yy_xyzz_1[j] + 0.5 * fl1_fx * tex_yyy_xzz_0[j] - 0.5 * fl1_fx * tex_yyy_xzz_1[j];

                tey_yyyy_xyzz_0[j] = pa_y[j] * tey_yyy_xyzz_0[j] - pc_y[j] * tey_yyy_xyzz_1[j] + 1.5 * fl1_fx * tey_yy_xyzz_0[j] -
                                     1.5 * fl1_fx * tey_yy_xyzz_1[j] + 0.5 * fl1_fx * tey_yyy_xzz_0[j] - 0.5 * fl1_fx * tey_yyy_xzz_1[j] +
                                     ta_yyy_xyzz_1[j];

                tez_yyyy_xyzz_0[j] = pa_y[j] * tez_yyy_xyzz_0[j] - pc_y[j] * tez_yyy_xyzz_1[j] + 1.5 * fl1_fx * tez_yy_xyzz_0[j] -
                                     1.5 * fl1_fx * tez_yy_xyzz_1[j] + 0.5 * fl1_fx * tez_yyy_xzz_0[j] - 0.5 * fl1_fx * tez_yyy_xzz_1[j];

                tex_yyyy_xzzz_0[j] =
                    pa_y[j] * tex_yyy_xzzz_0[j] - pc_y[j] * tex_yyy_xzzz_1[j] + 1.5 * fl1_fx * tex_yy_xzzz_0[j] - 1.5 * fl1_fx * tex_yy_xzzz_1[j];

                tey_yyyy_xzzz_0[j] = pa_y[j] * tey_yyy_xzzz_0[j] - pc_y[j] * tey_yyy_xzzz_1[j] + 1.5 * fl1_fx * tey_yy_xzzz_0[j] -
                                     1.5 * fl1_fx * tey_yy_xzzz_1[j] + ta_yyy_xzzz_1[j];

                tez_yyyy_xzzz_0[j] =
                    pa_y[j] * tez_yyy_xzzz_0[j] - pc_y[j] * tez_yyy_xzzz_1[j] + 1.5 * fl1_fx * tez_yy_xzzz_0[j] - 1.5 * fl1_fx * tez_yy_xzzz_1[j];

                tex_yyyy_yyyy_0[j] = pa_y[j] * tex_yyy_yyyy_0[j] - pc_y[j] * tex_yyy_yyyy_1[j] + 1.5 * fl1_fx * tex_yy_yyyy_0[j] -
                                     1.5 * fl1_fx * tex_yy_yyyy_1[j] + 2.0 * fl1_fx * tex_yyy_yyy_0[j] - 2.0 * fl1_fx * tex_yyy_yyy_1[j];

                tey_yyyy_yyyy_0[j] = pa_y[j] * tey_yyy_yyyy_0[j] - pc_y[j] * tey_yyy_yyyy_1[j] + 1.5 * fl1_fx * tey_yy_yyyy_0[j] -
                                     1.5 * fl1_fx * tey_yy_yyyy_1[j] + 2.0 * fl1_fx * tey_yyy_yyy_0[j] - 2.0 * fl1_fx * tey_yyy_yyy_1[j] +
                                     ta_yyy_yyyy_1[j];

                tez_yyyy_yyyy_0[j] = pa_y[j] * tez_yyy_yyyy_0[j] - pc_y[j] * tez_yyy_yyyy_1[j] + 1.5 * fl1_fx * tez_yy_yyyy_0[j] -
                                     1.5 * fl1_fx * tez_yy_yyyy_1[j] + 2.0 * fl1_fx * tez_yyy_yyy_0[j] - 2.0 * fl1_fx * tez_yyy_yyy_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForGG_483_531(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pcDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto)
{
    // Batch of Integrals (483,531)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_4_4_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_2_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(3 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            // set up pointers to tensors product of distances R(PC) = P - C

            auto pc_y = pcDistances.data(3 * idx + 1);

            // set up pointers to auxilary integrals

            auto tex_yyy_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 101);

            auto tey_yyy_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 101);

            auto tez_yyy_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 101);

            auto tex_yyy_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 102);

            auto tey_yyy_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 102);

            auto tez_yyy_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 102);

            auto tex_yyy_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 103);

            auto tey_yyy_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 103);

            auto tez_yyy_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 103);

            auto tex_yyy_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 104);

            auto tey_yyy_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 104);

            auto tez_yyy_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 104);

            auto tex_yyz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 105);

            auto tey_yyz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 105);

            auto tez_yyz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 105);

            auto tex_yyz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 106);

            auto tey_yyz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 106);

            auto tez_yyz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 106);

            auto tex_yyz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 107);

            auto tey_yyz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 107);

            auto tez_yyz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 107);

            auto tex_yyz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 108);

            auto tey_yyz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 108);

            auto tez_yyz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 108);

            auto tex_yyz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 109);

            auto tey_yyz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 109);

            auto tez_yyz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 109);

            auto tex_yyz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 110);

            auto tey_yyz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 110);

            auto tez_yyz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 110);

            auto tex_yyz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 111);

            auto tey_yyz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 111);

            auto tez_yyz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 111);

            auto tex_yyz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 112);

            auto tey_yyz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 112);

            auto tez_yyz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 112);

            auto tex_yyz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 113);

            auto tey_yyz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 113);

            auto tez_yyz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 113);

            auto tex_yyz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 114);

            auto tey_yyz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 114);

            auto tez_yyz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 114);

            auto tex_yyz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 115);

            auto tey_yyz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 115);

            auto tez_yyz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 115);

            auto tex_yyz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 116);

            auto tey_yyz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 116);

            auto tez_yyz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 116);

            auto tex_yyy_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 101);

            auto tey_yyy_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 101);

            auto tez_yyy_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 101);

            auto tex_yyy_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 102);

            auto tey_yyy_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 102);

            auto tez_yyy_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 102);

            auto tex_yyy_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 103);

            auto tey_yyy_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 103);

            auto tez_yyy_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 103);

            auto tex_yyy_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 104);

            auto tey_yyy_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 104);

            auto tez_yyy_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 104);

            auto tex_yyz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 105);

            auto tey_yyz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 105);

            auto tez_yyz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 105);

            auto tex_yyz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 106);

            auto tey_yyz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 106);

            auto tez_yyz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 106);

            auto tex_yyz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 107);

            auto tey_yyz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 107);

            auto tez_yyz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 107);

            auto tex_yyz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 108);

            auto tey_yyz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 108);

            auto tez_yyz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 108);

            auto tex_yyz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 109);

            auto tey_yyz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 109);

            auto tez_yyz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 109);

            auto tex_yyz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 110);

            auto tey_yyz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 110);

            auto tez_yyz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 110);

            auto tex_yyz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 111);

            auto tey_yyz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 111);

            auto tez_yyz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 111);

            auto tex_yyz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 112);

            auto tey_yyz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 112);

            auto tez_yyz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 112);

            auto tex_yyz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 113);

            auto tey_yyz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 113);

            auto tez_yyz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 113);

            auto tex_yyz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 114);

            auto tey_yyz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 114);

            auto tez_yyz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 114);

            auto tex_yyz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 115);

            auto tey_yyz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 115);

            auto tez_yyz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 115);

            auto tex_yyz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 116);

            auto tey_yyz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 116);

            auto tez_yyz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 116);

            auto tex_yy_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 56);

            auto tey_yy_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 56);

            auto tez_yy_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 56);

            auto tex_yy_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 57);

            auto tey_yy_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 57);

            auto tez_yy_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 57);

            auto tex_yy_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 58);

            auto tey_yy_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 58);

            auto tez_yy_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 58);

            auto tex_yy_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 59);

            auto tey_yy_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 59);

            auto tez_yy_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 59);

            auto tex_yz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 60);

            auto tey_yz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 60);

            auto tez_yz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 60);

            auto tex_yz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 61);

            auto tey_yz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 61);

            auto tez_yz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 61);

            auto tex_yz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 62);

            auto tey_yz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 62);

            auto tez_yz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 62);

            auto tex_yz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 63);

            auto tey_yz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 63);

            auto tez_yz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 63);

            auto tex_yz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 64);

            auto tey_yz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 64);

            auto tez_yz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 64);

            auto tex_yz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 65);

            auto tey_yz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 65);

            auto tez_yz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 65);

            auto tex_yz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 66);

            auto tey_yz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 66);

            auto tez_yz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 66);

            auto tex_yz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 67);

            auto tey_yz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 67);

            auto tez_yz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 67);

            auto tex_yz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 68);

            auto tey_yz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 68);

            auto tez_yz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 68);

            auto tex_yz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 69);

            auto tey_yz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 69);

            auto tez_yz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 69);

            auto tex_yz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 70);

            auto tey_yz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 70);

            auto tez_yz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 70);

            auto tex_yz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 71);

            auto tey_yz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 71);

            auto tez_yz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 71);

            auto tex_yy_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 56);

            auto tey_yy_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 56);

            auto tez_yy_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 56);

            auto tex_yy_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 57);

            auto tey_yy_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 57);

            auto tez_yy_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 57);

            auto tex_yy_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 58);

            auto tey_yy_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 58);

            auto tez_yy_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 58);

            auto tex_yy_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 59);

            auto tey_yy_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 59);

            auto tez_yy_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 59);

            auto tex_yz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 60);

            auto tey_yz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 60);

            auto tez_yz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 60);

            auto tex_yz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 61);

            auto tey_yz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 61);

            auto tez_yz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 61);

            auto tex_yz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 62);

            auto tey_yz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 62);

            auto tez_yz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 62);

            auto tex_yz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 63);

            auto tey_yz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 63);

            auto tez_yz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 63);

            auto tex_yz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 64);

            auto tey_yz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 64);

            auto tez_yz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 64);

            auto tex_yz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 65);

            auto tey_yz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 65);

            auto tez_yz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 65);

            auto tex_yz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 66);

            auto tey_yz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 66);

            auto tez_yz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 66);

            auto tex_yz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 67);

            auto tey_yz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 67);

            auto tez_yz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 67);

            auto tex_yz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 68);

            auto tey_yz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 68);

            auto tez_yz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 68);

            auto tex_yz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 69);

            auto tey_yz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 69);

            auto tez_yz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 69);

            auto tex_yz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 70);

            auto tey_yz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 70);

            auto tez_yz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 70);

            auto tex_yz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 71);

            auto tey_yz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 71);

            auto tez_yz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 71);

            auto tex_yyy_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 67);

            auto tey_yyy_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 67);

            auto tez_yyy_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 67);

            auto tex_yyy_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 68);

            auto tey_yyy_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 68);

            auto tez_yyy_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 68);

            auto tex_yyy_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 69);

            auto tey_yyy_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 69);

            auto tez_yyy_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 69);

            auto tex_yyz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 70);

            auto tey_yyz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 70);

            auto tez_yyz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 70);

            auto tex_yyz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 71);

            auto tey_yyz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 71);

            auto tez_yyz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 71);

            auto tex_yyz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 72);

            auto tey_yyz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 72);

            auto tez_yyz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 72);

            auto tex_yyz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 73);

            auto tey_yyz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 73);

            auto tez_yyz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 73);

            auto tex_yyz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 74);

            auto tey_yyz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 74);

            auto tez_yyz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 74);

            auto tex_yyz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 75);

            auto tey_yyz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 75);

            auto tez_yyz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 75);

            auto tex_yyz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 76);

            auto tey_yyz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 76);

            auto tez_yyz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 76);

            auto tex_yyz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 77);

            auto tey_yyz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 77);

            auto tez_yyz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 77);

            auto tex_yyy_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 67);

            auto tey_yyy_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 67);

            auto tez_yyy_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 67);

            auto tex_yyy_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 68);

            auto tey_yyy_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 68);

            auto tez_yyy_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 68);

            auto tex_yyy_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 69);

            auto tey_yyy_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 69);

            auto tez_yyy_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 69);

            auto tex_yyz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 70);

            auto tey_yyz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 70);

            auto tez_yyz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 70);

            auto tex_yyz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 71);

            auto tey_yyz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 71);

            auto tez_yyz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 71);

            auto tex_yyz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 72);

            auto tey_yyz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 72);

            auto tez_yyz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 72);

            auto tex_yyz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 73);

            auto tey_yyz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 73);

            auto tez_yyz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 73);

            auto tex_yyz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 74);

            auto tey_yyz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 74);

            auto tez_yyz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 74);

            auto tex_yyz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 75);

            auto tey_yyz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 75);

            auto tez_yyz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 75);

            auto tex_yyz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 76);

            auto tey_yyz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 76);

            auto tez_yyz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 76);

            auto tex_yyz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 77);

            auto tey_yyz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 77);

            auto tez_yyz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 77);

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

            // set up pointers to integrals

            auto tex_yyyy_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 161);

            auto tey_yyyy_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 161);

            auto tez_yyyy_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 161);

            auto tex_yyyy_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 162);

            auto tey_yyyy_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 162);

            auto tez_yyyy_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 162);

            auto tex_yyyy_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 163);

            auto tey_yyyy_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 163);

            auto tez_yyyy_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 163);

            auto tex_yyyy_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 164);

            auto tey_yyyy_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 164);

            auto tez_yyyy_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 164);

            auto tex_yyyz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 165);

            auto tey_yyyz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 165);

            auto tez_yyyz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 165);

            auto tex_yyyz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 166);

            auto tey_yyyz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 166);

            auto tez_yyyz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 166);

            auto tex_yyyz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 167);

            auto tey_yyyz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 167);

            auto tez_yyyz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 167);

            auto tex_yyyz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 168);

            auto tey_yyyz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 168);

            auto tez_yyyz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 168);

            auto tex_yyyz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 169);

            auto tey_yyyz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 169);

            auto tez_yyyz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 169);

            auto tex_yyyz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 170);

            auto tey_yyyz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 170);

            auto tez_yyyz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 170);

            auto tex_yyyz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 171);

            auto tey_yyyz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 171);

            auto tez_yyyz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 171);

            auto tex_yyyz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 172);

            auto tey_yyyz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 172);

            auto tez_yyyz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 172);

            auto tex_yyyz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 173);

            auto tey_yyyz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 173);

            auto tez_yyyz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 173);

            auto tex_yyyz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 174);

            auto tey_yyyz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 174);

            auto tez_yyyz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 174);

            auto tex_yyyz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 175);

            auto tey_yyyz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 175);

            auto tez_yyyz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 175);

            auto tex_yyyz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 176);

            auto tey_yyyz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 176);

            auto tez_yyyz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 176);

            // Batch of Integrals (483,531)

            #pragma omp simd aligned(fx, pa_y, pc_y, ta_yyy_yyyz_1, ta_yyy_yyzz_1, ta_yyy_yzzz_1, \
                                         ta_yyy_zzzz_1, ta_yyz_xxxx_1, ta_yyz_xxxy_1, ta_yyz_xxxz_1, ta_yyz_xxyy_1, \
                                         ta_yyz_xxyz_1, ta_yyz_xxzz_1, ta_yyz_xyyy_1, ta_yyz_xyyz_1, ta_yyz_xyzz_1, \
                                         ta_yyz_xzzz_1, ta_yyz_yyyy_1, ta_yyz_yyyz_1, tex_yy_yyyz_0, tex_yy_yyyz_1, \
                                         tex_yy_yyzz_0, tex_yy_yyzz_1, tex_yy_yzzz_0, tex_yy_yzzz_1, tex_yy_zzzz_0, \
                                         tex_yy_zzzz_1, tex_yyy_yyyz_0, tex_yyy_yyyz_1, tex_yyy_yyz_0, tex_yyy_yyz_1, \
                                         tex_yyy_yyzz_0, tex_yyy_yyzz_1, tex_yyy_yzz_0, tex_yyy_yzz_1, tex_yyy_yzzz_0, \
                                         tex_yyy_yzzz_1, tex_yyy_zzz_0, tex_yyy_zzz_1, tex_yyy_zzzz_0, tex_yyy_zzzz_1, \
                                         tex_yyyy_yyyz_0, tex_yyyy_yyzz_0, tex_yyyy_yzzz_0, tex_yyyy_zzzz_0, tex_yyyz_xxxx_0, \
                                         tex_yyyz_xxxy_0, tex_yyyz_xxxz_0, tex_yyyz_xxyy_0, tex_yyyz_xxyz_0, tex_yyyz_xxzz_0, \
                                         tex_yyyz_xyyy_0, tex_yyyz_xyyz_0, tex_yyyz_xyzz_0, tex_yyyz_xzzz_0, tex_yyyz_yyyy_0, \
                                         tex_yyyz_yyyz_0, tex_yyz_xxx_0, tex_yyz_xxx_1, tex_yyz_xxxx_0, tex_yyz_xxxx_1, \
                                         tex_yyz_xxxy_0, tex_yyz_xxxy_1, tex_yyz_xxxz_0, tex_yyz_xxxz_1, tex_yyz_xxy_0, \
                                         tex_yyz_xxy_1, tex_yyz_xxyy_0, tex_yyz_xxyy_1, tex_yyz_xxyz_0, tex_yyz_xxyz_1, \
                                         tex_yyz_xxz_0, tex_yyz_xxz_1, tex_yyz_xxzz_0, tex_yyz_xxzz_1, tex_yyz_xyy_0, \
                                         tex_yyz_xyy_1, tex_yyz_xyyy_0, tex_yyz_xyyy_1, tex_yyz_xyyz_0, tex_yyz_xyyz_1, \
                                         tex_yyz_xyz_0, tex_yyz_xyz_1, tex_yyz_xyzz_0, tex_yyz_xyzz_1, tex_yyz_xzz_0, \
                                         tex_yyz_xzz_1, tex_yyz_xzzz_0, tex_yyz_xzzz_1, tex_yyz_yyy_0, tex_yyz_yyy_1, \
                                         tex_yyz_yyyy_0, tex_yyz_yyyy_1, tex_yyz_yyyz_0, tex_yyz_yyyz_1, tex_yyz_yyz_0, \
                                         tex_yyz_yyz_1, tex_yz_xxxx_0, tex_yz_xxxx_1, tex_yz_xxxy_0, tex_yz_xxxy_1, \
                                         tex_yz_xxxz_0, tex_yz_xxxz_1, tex_yz_xxyy_0, tex_yz_xxyy_1, tex_yz_xxyz_0, \
                                         tex_yz_xxyz_1, tex_yz_xxzz_0, tex_yz_xxzz_1, tex_yz_xyyy_0, tex_yz_xyyy_1, \
                                         tex_yz_xyyz_0, tex_yz_xyyz_1, tex_yz_xyzz_0, tex_yz_xyzz_1, tex_yz_xzzz_0, \
                                         tex_yz_xzzz_1, tex_yz_yyyy_0, tex_yz_yyyy_1, tex_yz_yyyz_0, tex_yz_yyyz_1, \
                                         tey_yy_yyyz_0, tey_yy_yyyz_1, tey_yy_yyzz_0, tey_yy_yyzz_1, tey_yy_yzzz_0, \
                                         tey_yy_yzzz_1, tey_yy_zzzz_0, tey_yy_zzzz_1, tey_yyy_yyyz_0, tey_yyy_yyyz_1, \
                                         tey_yyy_yyz_0, tey_yyy_yyz_1, tey_yyy_yyzz_0, tey_yyy_yyzz_1, tey_yyy_yzz_0, \
                                         tey_yyy_yzz_1, tey_yyy_yzzz_0, tey_yyy_yzzz_1, tey_yyy_zzz_0, tey_yyy_zzz_1, \
                                         tey_yyy_zzzz_0, tey_yyy_zzzz_1, tey_yyyy_yyyz_0, tey_yyyy_yyzz_0, tey_yyyy_yzzz_0, \
                                         tey_yyyy_zzzz_0, tey_yyyz_xxxx_0, tey_yyyz_xxxy_0, tey_yyyz_xxxz_0, tey_yyyz_xxyy_0, \
                                         tey_yyyz_xxyz_0, tey_yyyz_xxzz_0, tey_yyyz_xyyy_0, tey_yyyz_xyyz_0, tey_yyyz_xyzz_0, \
                                         tey_yyyz_xzzz_0, tey_yyyz_yyyy_0, tey_yyyz_yyyz_0, tey_yyz_xxx_0, tey_yyz_xxx_1, \
                                         tey_yyz_xxxx_0, tey_yyz_xxxx_1, tey_yyz_xxxy_0, tey_yyz_xxxy_1, tey_yyz_xxxz_0, \
                                         tey_yyz_xxxz_1, tey_yyz_xxy_0, tey_yyz_xxy_1, tey_yyz_xxyy_0, tey_yyz_xxyy_1, \
                                         tey_yyz_xxyz_0, tey_yyz_xxyz_1, tey_yyz_xxz_0, tey_yyz_xxz_1, tey_yyz_xxzz_0, \
                                         tey_yyz_xxzz_1, tey_yyz_xyy_0, tey_yyz_xyy_1, tey_yyz_xyyy_0, tey_yyz_xyyy_1, \
                                         tey_yyz_xyyz_0, tey_yyz_xyyz_1, tey_yyz_xyz_0, tey_yyz_xyz_1, tey_yyz_xyzz_0, \
                                         tey_yyz_xyzz_1, tey_yyz_xzz_0, tey_yyz_xzz_1, tey_yyz_xzzz_0, tey_yyz_xzzz_1, \
                                         tey_yyz_yyy_0, tey_yyz_yyy_1, tey_yyz_yyyy_0, tey_yyz_yyyy_1, tey_yyz_yyyz_0, \
                                         tey_yyz_yyyz_1, tey_yyz_yyz_0, tey_yyz_yyz_1, tey_yz_xxxx_0, tey_yz_xxxx_1, \
                                         tey_yz_xxxy_0, tey_yz_xxxy_1, tey_yz_xxxz_0, tey_yz_xxxz_1, tey_yz_xxyy_0, \
                                         tey_yz_xxyy_1, tey_yz_xxyz_0, tey_yz_xxyz_1, tey_yz_xxzz_0, tey_yz_xxzz_1, \
                                         tey_yz_xyyy_0, tey_yz_xyyy_1, tey_yz_xyyz_0, tey_yz_xyyz_1, tey_yz_xyzz_0, \
                                         tey_yz_xyzz_1, tey_yz_xzzz_0, tey_yz_xzzz_1, tey_yz_yyyy_0, tey_yz_yyyy_1, \
                                         tey_yz_yyyz_0, tey_yz_yyyz_1, tez_yy_yyyz_0, tez_yy_yyyz_1, tez_yy_yyzz_0, \
                                         tez_yy_yyzz_1, tez_yy_yzzz_0, tez_yy_yzzz_1, tez_yy_zzzz_0, tez_yy_zzzz_1, \
                                         tez_yyy_yyyz_0, tez_yyy_yyyz_1, tez_yyy_yyz_0, tez_yyy_yyz_1, tez_yyy_yyzz_0, \
                                         tez_yyy_yyzz_1, tez_yyy_yzz_0, tez_yyy_yzz_1, tez_yyy_yzzz_0, tez_yyy_yzzz_1, \
                                         tez_yyy_zzz_0, tez_yyy_zzz_1, tez_yyy_zzzz_0, tez_yyy_zzzz_1, tez_yyyy_yyyz_0, \
                                         tez_yyyy_yyzz_0, tez_yyyy_yzzz_0, tez_yyyy_zzzz_0, tez_yyyz_xxxx_0, tez_yyyz_xxxy_0, \
                                         tez_yyyz_xxxz_0, tez_yyyz_xxyy_0, tez_yyyz_xxyz_0, tez_yyyz_xxzz_0, tez_yyyz_xyyy_0, \
                                         tez_yyyz_xyyz_0, tez_yyyz_xyzz_0, tez_yyyz_xzzz_0, tez_yyyz_yyyy_0, tez_yyyz_yyyz_0, \
                                         tez_yyz_xxx_0, tez_yyz_xxx_1, tez_yyz_xxxx_0, tez_yyz_xxxx_1, tez_yyz_xxxy_0, \
                                         tez_yyz_xxxy_1, tez_yyz_xxxz_0, tez_yyz_xxxz_1, tez_yyz_xxy_0, tez_yyz_xxy_1, \
                                         tez_yyz_xxyy_0, tez_yyz_xxyy_1, tez_yyz_xxyz_0, tez_yyz_xxyz_1, tez_yyz_xxz_0, \
                                         tez_yyz_xxz_1, tez_yyz_xxzz_0, tez_yyz_xxzz_1, tez_yyz_xyy_0, tez_yyz_xyy_1, \
                                         tez_yyz_xyyy_0, tez_yyz_xyyy_1, tez_yyz_xyyz_0, tez_yyz_xyyz_1, tez_yyz_xyz_0, \
                                         tez_yyz_xyz_1, tez_yyz_xyzz_0, tez_yyz_xyzz_1, tez_yyz_xzz_0, tez_yyz_xzz_1, \
                                         tez_yyz_xzzz_0, tez_yyz_xzzz_1, tez_yyz_yyy_0, tez_yyz_yyy_1, tez_yyz_yyyy_0, \
                                         tez_yyz_yyyy_1, tez_yyz_yyyz_0, tez_yyz_yyyz_1, tez_yyz_yyz_0, tez_yyz_yyz_1, \
                                         tez_yz_xxxx_0, tez_yz_xxxx_1, tez_yz_xxxy_0, tez_yz_xxxy_1, tez_yz_xxxz_0, \
                                         tez_yz_xxxz_1, tez_yz_xxyy_0, tez_yz_xxyy_1, tez_yz_xxyz_0, tez_yz_xxyz_1, \
                                         tez_yz_xxzz_0, tez_yz_xxzz_1, tez_yz_xyyy_0, tez_yz_xyyy_1, tez_yz_xyyz_0, \
                                         tez_yz_xyyz_1, tez_yz_xyzz_0, tez_yz_xyzz_1, tez_yz_xzzz_0, tez_yz_xzzz_1, \
                                         tez_yz_yyyy_0, tez_yz_yyyy_1, tez_yz_yyyz_0, tez_yz_yyyz_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_yyyy_yyyz_0[j] = pa_y[j] * tex_yyy_yyyz_0[j] - pc_y[j] * tex_yyy_yyyz_1[j] + 1.5 * fl1_fx * tex_yy_yyyz_0[j] -
                                     1.5 * fl1_fx * tex_yy_yyyz_1[j] + 1.5 * fl1_fx * tex_yyy_yyz_0[j] - 1.5 * fl1_fx * tex_yyy_yyz_1[j];

                tey_yyyy_yyyz_0[j] = pa_y[j] * tey_yyy_yyyz_0[j] - pc_y[j] * tey_yyy_yyyz_1[j] + 1.5 * fl1_fx * tey_yy_yyyz_0[j] -
                                     1.5 * fl1_fx * tey_yy_yyyz_1[j] + 1.5 * fl1_fx * tey_yyy_yyz_0[j] - 1.5 * fl1_fx * tey_yyy_yyz_1[j] +
                                     ta_yyy_yyyz_1[j];

                tez_yyyy_yyyz_0[j] = pa_y[j] * tez_yyy_yyyz_0[j] - pc_y[j] * tez_yyy_yyyz_1[j] + 1.5 * fl1_fx * tez_yy_yyyz_0[j] -
                                     1.5 * fl1_fx * tez_yy_yyyz_1[j] + 1.5 * fl1_fx * tez_yyy_yyz_0[j] - 1.5 * fl1_fx * tez_yyy_yyz_1[j];

                tex_yyyy_yyzz_0[j] = pa_y[j] * tex_yyy_yyzz_0[j] - pc_y[j] * tex_yyy_yyzz_1[j] + 1.5 * fl1_fx * tex_yy_yyzz_0[j] -
                                     1.5 * fl1_fx * tex_yy_yyzz_1[j] + fl1_fx * tex_yyy_yzz_0[j] - fl1_fx * tex_yyy_yzz_1[j];

                tey_yyyy_yyzz_0[j] = pa_y[j] * tey_yyy_yyzz_0[j] - pc_y[j] * tey_yyy_yyzz_1[j] + 1.5 * fl1_fx * tey_yy_yyzz_0[j] -
                                     1.5 * fl1_fx * tey_yy_yyzz_1[j] + fl1_fx * tey_yyy_yzz_0[j] - fl1_fx * tey_yyy_yzz_1[j] + ta_yyy_yyzz_1[j];

                tez_yyyy_yyzz_0[j] = pa_y[j] * tez_yyy_yyzz_0[j] - pc_y[j] * tez_yyy_yyzz_1[j] + 1.5 * fl1_fx * tez_yy_yyzz_0[j] -
                                     1.5 * fl1_fx * tez_yy_yyzz_1[j] + fl1_fx * tez_yyy_yzz_0[j] - fl1_fx * tez_yyy_yzz_1[j];

                tex_yyyy_yzzz_0[j] = pa_y[j] * tex_yyy_yzzz_0[j] - pc_y[j] * tex_yyy_yzzz_1[j] + 1.5 * fl1_fx * tex_yy_yzzz_0[j] -
                                     1.5 * fl1_fx * tex_yy_yzzz_1[j] + 0.5 * fl1_fx * tex_yyy_zzz_0[j] - 0.5 * fl1_fx * tex_yyy_zzz_1[j];

                tey_yyyy_yzzz_0[j] = pa_y[j] * tey_yyy_yzzz_0[j] - pc_y[j] * tey_yyy_yzzz_1[j] + 1.5 * fl1_fx * tey_yy_yzzz_0[j] -
                                     1.5 * fl1_fx * tey_yy_yzzz_1[j] + 0.5 * fl1_fx * tey_yyy_zzz_0[j] - 0.5 * fl1_fx * tey_yyy_zzz_1[j] +
                                     ta_yyy_yzzz_1[j];

                tez_yyyy_yzzz_0[j] = pa_y[j] * tez_yyy_yzzz_0[j] - pc_y[j] * tez_yyy_yzzz_1[j] + 1.5 * fl1_fx * tez_yy_yzzz_0[j] -
                                     1.5 * fl1_fx * tez_yy_yzzz_1[j] + 0.5 * fl1_fx * tez_yyy_zzz_0[j] - 0.5 * fl1_fx * tez_yyy_zzz_1[j];

                tex_yyyy_zzzz_0[j] =
                    pa_y[j] * tex_yyy_zzzz_0[j] - pc_y[j] * tex_yyy_zzzz_1[j] + 1.5 * fl1_fx * tex_yy_zzzz_0[j] - 1.5 * fl1_fx * tex_yy_zzzz_1[j];

                tey_yyyy_zzzz_0[j] = pa_y[j] * tey_yyy_zzzz_0[j] - pc_y[j] * tey_yyy_zzzz_1[j] + 1.5 * fl1_fx * tey_yy_zzzz_0[j] -
                                     1.5 * fl1_fx * tey_yy_zzzz_1[j] + ta_yyy_zzzz_1[j];

                tez_yyyy_zzzz_0[j] =
                    pa_y[j] * tez_yyy_zzzz_0[j] - pc_y[j] * tez_yyy_zzzz_1[j] + 1.5 * fl1_fx * tez_yy_zzzz_0[j] - 1.5 * fl1_fx * tez_yy_zzzz_1[j];

                tex_yyyz_xxxx_0[j] =
                    pa_y[j] * tex_yyz_xxxx_0[j] - pc_y[j] * tex_yyz_xxxx_1[j] + fl1_fx * tex_yz_xxxx_0[j] - fl1_fx * tex_yz_xxxx_1[j];

                tey_yyyz_xxxx_0[j] = pa_y[j] * tey_yyz_xxxx_0[j] - pc_y[j] * tey_yyz_xxxx_1[j] + fl1_fx * tey_yz_xxxx_0[j] -
                                     fl1_fx * tey_yz_xxxx_1[j] + ta_yyz_xxxx_1[j];

                tez_yyyz_xxxx_0[j] =
                    pa_y[j] * tez_yyz_xxxx_0[j] - pc_y[j] * tez_yyz_xxxx_1[j] + fl1_fx * tez_yz_xxxx_0[j] - fl1_fx * tez_yz_xxxx_1[j];

                tex_yyyz_xxxy_0[j] = pa_y[j] * tex_yyz_xxxy_0[j] - pc_y[j] * tex_yyz_xxxy_1[j] + fl1_fx * tex_yz_xxxy_0[j] -
                                     fl1_fx * tex_yz_xxxy_1[j] + 0.5 * fl1_fx * tex_yyz_xxx_0[j] - 0.5 * fl1_fx * tex_yyz_xxx_1[j];

                tey_yyyz_xxxy_0[j] = pa_y[j] * tey_yyz_xxxy_0[j] - pc_y[j] * tey_yyz_xxxy_1[j] + fl1_fx * tey_yz_xxxy_0[j] -
                                     fl1_fx * tey_yz_xxxy_1[j] + 0.5 * fl1_fx * tey_yyz_xxx_0[j] - 0.5 * fl1_fx * tey_yyz_xxx_1[j] + ta_yyz_xxxy_1[j];

                tez_yyyz_xxxy_0[j] = pa_y[j] * tez_yyz_xxxy_0[j] - pc_y[j] * tez_yyz_xxxy_1[j] + fl1_fx * tez_yz_xxxy_0[j] -
                                     fl1_fx * tez_yz_xxxy_1[j] + 0.5 * fl1_fx * tez_yyz_xxx_0[j] - 0.5 * fl1_fx * tez_yyz_xxx_1[j];

                tex_yyyz_xxxz_0[j] =
                    pa_y[j] * tex_yyz_xxxz_0[j] - pc_y[j] * tex_yyz_xxxz_1[j] + fl1_fx * tex_yz_xxxz_0[j] - fl1_fx * tex_yz_xxxz_1[j];

                tey_yyyz_xxxz_0[j] = pa_y[j] * tey_yyz_xxxz_0[j] - pc_y[j] * tey_yyz_xxxz_1[j] + fl1_fx * tey_yz_xxxz_0[j] -
                                     fl1_fx * tey_yz_xxxz_1[j] + ta_yyz_xxxz_1[j];

                tez_yyyz_xxxz_0[j] =
                    pa_y[j] * tez_yyz_xxxz_0[j] - pc_y[j] * tez_yyz_xxxz_1[j] + fl1_fx * tez_yz_xxxz_0[j] - fl1_fx * tez_yz_xxxz_1[j];

                tex_yyyz_xxyy_0[j] = pa_y[j] * tex_yyz_xxyy_0[j] - pc_y[j] * tex_yyz_xxyy_1[j] + fl1_fx * tex_yz_xxyy_0[j] -
                                     fl1_fx * tex_yz_xxyy_1[j] + fl1_fx * tex_yyz_xxy_0[j] - fl1_fx * tex_yyz_xxy_1[j];

                tey_yyyz_xxyy_0[j] = pa_y[j] * tey_yyz_xxyy_0[j] - pc_y[j] * tey_yyz_xxyy_1[j] + fl1_fx * tey_yz_xxyy_0[j] -
                                     fl1_fx * tey_yz_xxyy_1[j] + fl1_fx * tey_yyz_xxy_0[j] - fl1_fx * tey_yyz_xxy_1[j] + ta_yyz_xxyy_1[j];

                tez_yyyz_xxyy_0[j] = pa_y[j] * tez_yyz_xxyy_0[j] - pc_y[j] * tez_yyz_xxyy_1[j] + fl1_fx * tez_yz_xxyy_0[j] -
                                     fl1_fx * tez_yz_xxyy_1[j] + fl1_fx * tez_yyz_xxy_0[j] - fl1_fx * tez_yyz_xxy_1[j];

                tex_yyyz_xxyz_0[j] = pa_y[j] * tex_yyz_xxyz_0[j] - pc_y[j] * tex_yyz_xxyz_1[j] + fl1_fx * tex_yz_xxyz_0[j] -
                                     fl1_fx * tex_yz_xxyz_1[j] + 0.5 * fl1_fx * tex_yyz_xxz_0[j] - 0.5 * fl1_fx * tex_yyz_xxz_1[j];

                tey_yyyz_xxyz_0[j] = pa_y[j] * tey_yyz_xxyz_0[j] - pc_y[j] * tey_yyz_xxyz_1[j] + fl1_fx * tey_yz_xxyz_0[j] -
                                     fl1_fx * tey_yz_xxyz_1[j] + 0.5 * fl1_fx * tey_yyz_xxz_0[j] - 0.5 * fl1_fx * tey_yyz_xxz_1[j] + ta_yyz_xxyz_1[j];

                tez_yyyz_xxyz_0[j] = pa_y[j] * tez_yyz_xxyz_0[j] - pc_y[j] * tez_yyz_xxyz_1[j] + fl1_fx * tez_yz_xxyz_0[j] -
                                     fl1_fx * tez_yz_xxyz_1[j] + 0.5 * fl1_fx * tez_yyz_xxz_0[j] - 0.5 * fl1_fx * tez_yyz_xxz_1[j];

                tex_yyyz_xxzz_0[j] =
                    pa_y[j] * tex_yyz_xxzz_0[j] - pc_y[j] * tex_yyz_xxzz_1[j] + fl1_fx * tex_yz_xxzz_0[j] - fl1_fx * tex_yz_xxzz_1[j];

                tey_yyyz_xxzz_0[j] = pa_y[j] * tey_yyz_xxzz_0[j] - pc_y[j] * tey_yyz_xxzz_1[j] + fl1_fx * tey_yz_xxzz_0[j] -
                                     fl1_fx * tey_yz_xxzz_1[j] + ta_yyz_xxzz_1[j];

                tez_yyyz_xxzz_0[j] =
                    pa_y[j] * tez_yyz_xxzz_0[j] - pc_y[j] * tez_yyz_xxzz_1[j] + fl1_fx * tez_yz_xxzz_0[j] - fl1_fx * tez_yz_xxzz_1[j];

                tex_yyyz_xyyy_0[j] = pa_y[j] * tex_yyz_xyyy_0[j] - pc_y[j] * tex_yyz_xyyy_1[j] + fl1_fx * tex_yz_xyyy_0[j] -
                                     fl1_fx * tex_yz_xyyy_1[j] + 1.5 * fl1_fx * tex_yyz_xyy_0[j] - 1.5 * fl1_fx * tex_yyz_xyy_1[j];

                tey_yyyz_xyyy_0[j] = pa_y[j] * tey_yyz_xyyy_0[j] - pc_y[j] * tey_yyz_xyyy_1[j] + fl1_fx * tey_yz_xyyy_0[j] -
                                     fl1_fx * tey_yz_xyyy_1[j] + 1.5 * fl1_fx * tey_yyz_xyy_0[j] - 1.5 * fl1_fx * tey_yyz_xyy_1[j] + ta_yyz_xyyy_1[j];

                tez_yyyz_xyyy_0[j] = pa_y[j] * tez_yyz_xyyy_0[j] - pc_y[j] * tez_yyz_xyyy_1[j] + fl1_fx * tez_yz_xyyy_0[j] -
                                     fl1_fx * tez_yz_xyyy_1[j] + 1.5 * fl1_fx * tez_yyz_xyy_0[j] - 1.5 * fl1_fx * tez_yyz_xyy_1[j];

                tex_yyyz_xyyz_0[j] = pa_y[j] * tex_yyz_xyyz_0[j] - pc_y[j] * tex_yyz_xyyz_1[j] + fl1_fx * tex_yz_xyyz_0[j] -
                                     fl1_fx * tex_yz_xyyz_1[j] + fl1_fx * tex_yyz_xyz_0[j] - fl1_fx * tex_yyz_xyz_1[j];

                tey_yyyz_xyyz_0[j] = pa_y[j] * tey_yyz_xyyz_0[j] - pc_y[j] * tey_yyz_xyyz_1[j] + fl1_fx * tey_yz_xyyz_0[j] -
                                     fl1_fx * tey_yz_xyyz_1[j] + fl1_fx * tey_yyz_xyz_0[j] - fl1_fx * tey_yyz_xyz_1[j] + ta_yyz_xyyz_1[j];

                tez_yyyz_xyyz_0[j] = pa_y[j] * tez_yyz_xyyz_0[j] - pc_y[j] * tez_yyz_xyyz_1[j] + fl1_fx * tez_yz_xyyz_0[j] -
                                     fl1_fx * tez_yz_xyyz_1[j] + fl1_fx * tez_yyz_xyz_0[j] - fl1_fx * tez_yyz_xyz_1[j];

                tex_yyyz_xyzz_0[j] = pa_y[j] * tex_yyz_xyzz_0[j] - pc_y[j] * tex_yyz_xyzz_1[j] + fl1_fx * tex_yz_xyzz_0[j] -
                                     fl1_fx * tex_yz_xyzz_1[j] + 0.5 * fl1_fx * tex_yyz_xzz_0[j] - 0.5 * fl1_fx * tex_yyz_xzz_1[j];

                tey_yyyz_xyzz_0[j] = pa_y[j] * tey_yyz_xyzz_0[j] - pc_y[j] * tey_yyz_xyzz_1[j] + fl1_fx * tey_yz_xyzz_0[j] -
                                     fl1_fx * tey_yz_xyzz_1[j] + 0.5 * fl1_fx * tey_yyz_xzz_0[j] - 0.5 * fl1_fx * tey_yyz_xzz_1[j] + ta_yyz_xyzz_1[j];

                tez_yyyz_xyzz_0[j] = pa_y[j] * tez_yyz_xyzz_0[j] - pc_y[j] * tez_yyz_xyzz_1[j] + fl1_fx * tez_yz_xyzz_0[j] -
                                     fl1_fx * tez_yz_xyzz_1[j] + 0.5 * fl1_fx * tez_yyz_xzz_0[j] - 0.5 * fl1_fx * tez_yyz_xzz_1[j];

                tex_yyyz_xzzz_0[j] =
                    pa_y[j] * tex_yyz_xzzz_0[j] - pc_y[j] * tex_yyz_xzzz_1[j] + fl1_fx * tex_yz_xzzz_0[j] - fl1_fx * tex_yz_xzzz_1[j];

                tey_yyyz_xzzz_0[j] = pa_y[j] * tey_yyz_xzzz_0[j] - pc_y[j] * tey_yyz_xzzz_1[j] + fl1_fx * tey_yz_xzzz_0[j] -
                                     fl1_fx * tey_yz_xzzz_1[j] + ta_yyz_xzzz_1[j];

                tez_yyyz_xzzz_0[j] =
                    pa_y[j] * tez_yyz_xzzz_0[j] - pc_y[j] * tez_yyz_xzzz_1[j] + fl1_fx * tez_yz_xzzz_0[j] - fl1_fx * tez_yz_xzzz_1[j];

                tex_yyyz_yyyy_0[j] = pa_y[j] * tex_yyz_yyyy_0[j] - pc_y[j] * tex_yyz_yyyy_1[j] + fl1_fx * tex_yz_yyyy_0[j] -
                                     fl1_fx * tex_yz_yyyy_1[j] + 2.0 * fl1_fx * tex_yyz_yyy_0[j] - 2.0 * fl1_fx * tex_yyz_yyy_1[j];

                tey_yyyz_yyyy_0[j] = pa_y[j] * tey_yyz_yyyy_0[j] - pc_y[j] * tey_yyz_yyyy_1[j] + fl1_fx * tey_yz_yyyy_0[j] -
                                     fl1_fx * tey_yz_yyyy_1[j] + 2.0 * fl1_fx * tey_yyz_yyy_0[j] - 2.0 * fl1_fx * tey_yyz_yyy_1[j] + ta_yyz_yyyy_1[j];

                tez_yyyz_yyyy_0[j] = pa_y[j] * tez_yyz_yyyy_0[j] - pc_y[j] * tez_yyz_yyyy_1[j] + fl1_fx * tez_yz_yyyy_0[j] -
                                     fl1_fx * tez_yz_yyyy_1[j] + 2.0 * fl1_fx * tez_yyz_yyy_0[j] - 2.0 * fl1_fx * tez_yyz_yyy_1[j];

                tex_yyyz_yyyz_0[j] = pa_y[j] * tex_yyz_yyyz_0[j] - pc_y[j] * tex_yyz_yyyz_1[j] + fl1_fx * tex_yz_yyyz_0[j] -
                                     fl1_fx * tex_yz_yyyz_1[j] + 1.5 * fl1_fx * tex_yyz_yyz_0[j] - 1.5 * fl1_fx * tex_yyz_yyz_1[j];

                tey_yyyz_yyyz_0[j] = pa_y[j] * tey_yyz_yyyz_0[j] - pc_y[j] * tey_yyz_yyyz_1[j] + fl1_fx * tey_yz_yyyz_0[j] -
                                     fl1_fx * tey_yz_yyyz_1[j] + 1.5 * fl1_fx * tey_yyz_yyz_0[j] - 1.5 * fl1_fx * tey_yyz_yyz_1[j] + ta_yyz_yyyz_1[j];

                tez_yyyz_yyyz_0[j] = pa_y[j] * tez_yyz_yyyz_0[j] - pc_y[j] * tez_yyz_yyyz_1[j] + fl1_fx * tez_yz_yyyz_0[j] -
                                     fl1_fx * tez_yz_yyyz_1[j] + 1.5 * fl1_fx * tez_yyz_yyz_0[j] - 1.5 * fl1_fx * tez_yyz_yyz_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForGG_531_579(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pcDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto)
{
    // Batch of Integrals (531,579)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_4_4_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_2_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(3 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            // set up pointers to tensors product of distances R(PC) = P - C

            auto pc_y = pcDistances.data(3 * idx + 1);

            // set up pointers to auxilary integrals

            auto tex_yyz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 117);

            auto tey_yyz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 117);

            auto tez_yyz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 117);

            auto tex_yyz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 118);

            auto tey_yyz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 118);

            auto tez_yyz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 118);

            auto tex_yyz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 119);

            auto tey_yyz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 119);

            auto tez_yyz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 119);

            auto tex_yzz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 120);

            auto tey_yzz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 120);

            auto tez_yzz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 120);

            auto tex_yzz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 121);

            auto tey_yzz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 121);

            auto tez_yzz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 121);

            auto tex_yzz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 122);

            auto tey_yzz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 122);

            auto tez_yzz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 122);

            auto tex_yzz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 123);

            auto tey_yzz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 123);

            auto tez_yzz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 123);

            auto tex_yzz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 124);

            auto tey_yzz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 124);

            auto tez_yzz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 124);

            auto tex_yzz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 125);

            auto tey_yzz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 125);

            auto tez_yzz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 125);

            auto tex_yzz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 126);

            auto tey_yzz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 126);

            auto tez_yzz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 126);

            auto tex_yzz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 127);

            auto tey_yzz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 127);

            auto tez_yzz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 127);

            auto tex_yzz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 128);

            auto tey_yzz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 128);

            auto tez_yzz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 128);

            auto tex_yzz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 129);

            auto tey_yzz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 129);

            auto tez_yzz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 129);

            auto tex_yzz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 130);

            auto tey_yzz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 130);

            auto tez_yzz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 130);

            auto tex_yzz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 131);

            auto tey_yzz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 131);

            auto tez_yzz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 131);

            auto tex_yzz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 132);

            auto tey_yzz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 132);

            auto tez_yzz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 132);

            auto tex_yyz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 117);

            auto tey_yyz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 117);

            auto tez_yyz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 117);

            auto tex_yyz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 118);

            auto tey_yyz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 118);

            auto tez_yyz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 118);

            auto tex_yyz_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 119);

            auto tey_yyz_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 119);

            auto tez_yyz_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 119);

            auto tex_yzz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 120);

            auto tey_yzz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 120);

            auto tez_yzz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 120);

            auto tex_yzz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 121);

            auto tey_yzz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 121);

            auto tez_yzz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 121);

            auto tex_yzz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 122);

            auto tey_yzz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 122);

            auto tez_yzz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 122);

            auto tex_yzz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 123);

            auto tey_yzz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 123);

            auto tez_yzz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 123);

            auto tex_yzz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 124);

            auto tey_yzz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 124);

            auto tez_yzz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 124);

            auto tex_yzz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 125);

            auto tey_yzz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 125);

            auto tez_yzz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 125);

            auto tex_yzz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 126);

            auto tey_yzz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 126);

            auto tez_yzz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 126);

            auto tex_yzz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 127);

            auto tey_yzz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 127);

            auto tez_yzz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 127);

            auto tex_yzz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 128);

            auto tey_yzz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 128);

            auto tez_yzz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 128);

            auto tex_yzz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 129);

            auto tey_yzz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 129);

            auto tez_yzz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 129);

            auto tex_yzz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 130);

            auto tey_yzz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 130);

            auto tez_yzz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 130);

            auto tex_yzz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 131);

            auto tey_yzz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 131);

            auto tez_yzz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 131);

            auto tex_yzz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 132);

            auto tey_yzz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 132);

            auto tez_yzz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 132);

            auto tex_yz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 72);

            auto tey_yz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 72);

            auto tez_yz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 72);

            auto tex_yz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 73);

            auto tey_yz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 73);

            auto tez_yz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 73);

            auto tex_yz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 74);

            auto tey_yz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 74);

            auto tez_yz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 74);

            auto tex_zz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 75);

            auto tey_zz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 75);

            auto tez_zz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 75);

            auto tex_zz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 76);

            auto tey_zz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 76);

            auto tez_zz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 76);

            auto tex_zz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 77);

            auto tey_zz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 77);

            auto tez_zz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 77);

            auto tex_zz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 78);

            auto tey_zz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 78);

            auto tez_zz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 78);

            auto tex_zz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 79);

            auto tey_zz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 79);

            auto tez_zz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 79);

            auto tex_zz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 80);

            auto tey_zz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 80);

            auto tez_zz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 80);

            auto tex_zz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 81);

            auto tey_zz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 81);

            auto tez_zz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 81);

            auto tex_zz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 82);

            auto tey_zz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 82);

            auto tez_zz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 82);

            auto tex_zz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 83);

            auto tey_zz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 83);

            auto tez_zz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 83);

            auto tex_zz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 84);

            auto tey_zz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 84);

            auto tez_zz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 84);

            auto tex_zz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 85);

            auto tey_zz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 85);

            auto tez_zz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 85);

            auto tex_zz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 86);

            auto tey_zz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 86);

            auto tez_zz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 86);

            auto tex_zz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 87);

            auto tey_zz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 87);

            auto tez_zz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 87);

            auto tex_yz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 72);

            auto tey_yz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 72);

            auto tez_yz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 72);

            auto tex_yz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 73);

            auto tey_yz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 73);

            auto tez_yz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 73);

            auto tex_yz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 74);

            auto tey_yz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 74);

            auto tez_yz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 74);

            auto tex_zz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 75);

            auto tey_zz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 75);

            auto tez_zz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 75);

            auto tex_zz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 76);

            auto tey_zz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 76);

            auto tez_zz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 76);

            auto tex_zz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 77);

            auto tey_zz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 77);

            auto tez_zz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 77);

            auto tex_zz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 78);

            auto tey_zz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 78);

            auto tez_zz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 78);

            auto tex_zz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 79);

            auto tey_zz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 79);

            auto tez_zz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 79);

            auto tex_zz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 80);

            auto tey_zz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 80);

            auto tez_zz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 80);

            auto tex_zz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 81);

            auto tey_zz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 81);

            auto tez_zz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 81);

            auto tex_zz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 82);

            auto tey_zz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 82);

            auto tez_zz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 82);

            auto tex_zz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 83);

            auto tey_zz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 83);

            auto tez_zz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 83);

            auto tex_zz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 84);

            auto tey_zz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 84);

            auto tez_zz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 84);

            auto tex_zz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 85);

            auto tey_zz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 85);

            auto tez_zz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 85);

            auto tex_zz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 86);

            auto tey_zz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 86);

            auto tez_zz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 86);

            auto tex_zz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 87);

            auto tey_zz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 87);

            auto tez_zz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 87);

            auto tex_yyz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 78);

            auto tey_yyz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 78);

            auto tez_yyz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 78);

            auto tex_yyz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 79);

            auto tey_yyz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 79);

            auto tez_yyz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 79);

            auto tex_yzz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 80);

            auto tey_yzz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 80);

            auto tez_yzz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 80);

            auto tex_yzz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 81);

            auto tey_yzz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 81);

            auto tez_yzz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 81);

            auto tex_yzz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 82);

            auto tey_yzz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 82);

            auto tez_yzz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 82);

            auto tex_yzz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 83);

            auto tey_yzz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 83);

            auto tez_yzz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 83);

            auto tex_yzz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 84);

            auto tey_yzz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 84);

            auto tez_yzz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 84);

            auto tex_yzz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 85);

            auto tey_yzz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 85);

            auto tez_yzz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 85);

            auto tex_yzz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 86);

            auto tey_yzz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 86);

            auto tez_yzz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 86);

            auto tex_yzz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 87);

            auto tey_yzz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 87);

            auto tez_yzz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 87);

            auto tex_yzz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 88);

            auto tey_yzz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 88);

            auto tez_yzz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 88);

            auto tex_yyz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 78);

            auto tey_yyz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 78);

            auto tez_yyz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 78);

            auto tex_yyz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 79);

            auto tey_yyz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 79);

            auto tez_yyz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 79);

            auto tex_yzz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 80);

            auto tey_yzz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 80);

            auto tez_yzz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 80);

            auto tex_yzz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 81);

            auto tey_yzz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 81);

            auto tez_yzz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 81);

            auto tex_yzz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 82);

            auto tey_yzz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 82);

            auto tez_yzz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 82);

            auto tex_yzz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 83);

            auto tey_yzz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 83);

            auto tez_yzz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 83);

            auto tex_yzz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 84);

            auto tey_yzz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 84);

            auto tez_yzz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 84);

            auto tex_yzz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 85);

            auto tey_yzz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 85);

            auto tez_yzz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 85);

            auto tex_yzz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 86);

            auto tey_yzz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 86);

            auto tez_yzz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 86);

            auto tex_yzz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 87);

            auto tey_yzz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 87);

            auto tez_yzz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 87);

            auto tex_yzz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 88);

            auto tey_yzz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 88);

            auto tez_yzz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 88);

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

            // set up pointers to integrals

            auto tex_yyyz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 177);

            auto tey_yyyz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 177);

            auto tez_yyyz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 177);

            auto tex_yyyz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 178);

            auto tey_yyyz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 178);

            auto tez_yyyz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 178);

            auto tex_yyyz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 179);

            auto tey_yyyz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 179);

            auto tez_yyyz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 179);

            auto tex_yyzz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 180);

            auto tey_yyzz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 180);

            auto tez_yyzz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 180);

            auto tex_yyzz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 181);

            auto tey_yyzz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 181);

            auto tez_yyzz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 181);

            auto tex_yyzz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 182);

            auto tey_yyzz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 182);

            auto tez_yyzz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 182);

            auto tex_yyzz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 183);

            auto tey_yyzz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 183);

            auto tez_yyzz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 183);

            auto tex_yyzz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 184);

            auto tey_yyzz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 184);

            auto tez_yyzz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 184);

            auto tex_yyzz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 185);

            auto tey_yyzz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 185);

            auto tez_yyzz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 185);

            auto tex_yyzz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 186);

            auto tey_yyzz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 186);

            auto tez_yyzz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 186);

            auto tex_yyzz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 187);

            auto tey_yyzz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 187);

            auto tez_yyzz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 187);

            auto tex_yyzz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 188);

            auto tey_yyzz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 188);

            auto tez_yyzz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 188);

            auto tex_yyzz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 189);

            auto tey_yyzz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 189);

            auto tez_yyzz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 189);

            auto tex_yyzz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 190);

            auto tey_yyzz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 190);

            auto tez_yyzz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 190);

            auto tex_yyzz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 191);

            auto tey_yyzz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 191);

            auto tez_yyzz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 191);

            auto tex_yyzz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 192);

            auto tey_yyzz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 192);

            auto tez_yyzz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 192);

            // Batch of Integrals (531,579)

            #pragma omp simd aligned(fx, pa_y, pc_y, ta_yyz_yyzz_1, ta_yyz_yzzz_1, ta_yyz_zzzz_1, \
                                         ta_yzz_xxxx_1, ta_yzz_xxxy_1, ta_yzz_xxxz_1, ta_yzz_xxyy_1, ta_yzz_xxyz_1, \
                                         ta_yzz_xxzz_1, ta_yzz_xyyy_1, ta_yzz_xyyz_1, ta_yzz_xyzz_1, ta_yzz_xzzz_1, \
                                         ta_yzz_yyyy_1, ta_yzz_yyyz_1, ta_yzz_yyzz_1, tex_yyyz_yyzz_0, tex_yyyz_yzzz_0, \
                                         tex_yyyz_zzzz_0, tex_yyz_yyzz_0, tex_yyz_yyzz_1, tex_yyz_yzz_0, tex_yyz_yzz_1, \
                                         tex_yyz_yzzz_0, tex_yyz_yzzz_1, tex_yyz_zzz_0, tex_yyz_zzz_1, tex_yyz_zzzz_0, \
                                         tex_yyz_zzzz_1, tex_yyzz_xxxx_0, tex_yyzz_xxxy_0, tex_yyzz_xxxz_0, tex_yyzz_xxyy_0, \
                                         tex_yyzz_xxyz_0, tex_yyzz_xxzz_0, tex_yyzz_xyyy_0, tex_yyzz_xyyz_0, tex_yyzz_xyzz_0, \
                                         tex_yyzz_xzzz_0, tex_yyzz_yyyy_0, tex_yyzz_yyyz_0, tex_yyzz_yyzz_0, tex_yz_yyzz_0, \
                                         tex_yz_yyzz_1, tex_yz_yzzz_0, tex_yz_yzzz_1, tex_yz_zzzz_0, tex_yz_zzzz_1, \
                                         tex_yzz_xxx_0, tex_yzz_xxx_1, tex_yzz_xxxx_0, tex_yzz_xxxx_1, tex_yzz_xxxy_0, \
                                         tex_yzz_xxxy_1, tex_yzz_xxxz_0, tex_yzz_xxxz_1, tex_yzz_xxy_0, tex_yzz_xxy_1, \
                                         tex_yzz_xxyy_0, tex_yzz_xxyy_1, tex_yzz_xxyz_0, tex_yzz_xxyz_1, tex_yzz_xxz_0, \
                                         tex_yzz_xxz_1, tex_yzz_xxzz_0, tex_yzz_xxzz_1, tex_yzz_xyy_0, tex_yzz_xyy_1, \
                                         tex_yzz_xyyy_0, tex_yzz_xyyy_1, tex_yzz_xyyz_0, tex_yzz_xyyz_1, tex_yzz_xyz_0, \
                                         tex_yzz_xyz_1, tex_yzz_xyzz_0, tex_yzz_xyzz_1, tex_yzz_xzz_0, tex_yzz_xzz_1, \
                                         tex_yzz_xzzz_0, tex_yzz_xzzz_1, tex_yzz_yyy_0, tex_yzz_yyy_1, tex_yzz_yyyy_0, \
                                         tex_yzz_yyyy_1, tex_yzz_yyyz_0, tex_yzz_yyyz_1, tex_yzz_yyz_0, tex_yzz_yyz_1, \
                                         tex_yzz_yyzz_0, tex_yzz_yyzz_1, tex_yzz_yzz_0, tex_yzz_yzz_1, tex_zz_xxxx_0, \
                                         tex_zz_xxxx_1, tex_zz_xxxy_0, tex_zz_xxxy_1, tex_zz_xxxz_0, tex_zz_xxxz_1, \
                                         tex_zz_xxyy_0, tex_zz_xxyy_1, tex_zz_xxyz_0, tex_zz_xxyz_1, tex_zz_xxzz_0, \
                                         tex_zz_xxzz_1, tex_zz_xyyy_0, tex_zz_xyyy_1, tex_zz_xyyz_0, tex_zz_xyyz_1, \
                                         tex_zz_xyzz_0, tex_zz_xyzz_1, tex_zz_xzzz_0, tex_zz_xzzz_1, tex_zz_yyyy_0, \
                                         tex_zz_yyyy_1, tex_zz_yyyz_0, tex_zz_yyyz_1, tex_zz_yyzz_0, tex_zz_yyzz_1, \
                                         tey_yyyz_yyzz_0, tey_yyyz_yzzz_0, tey_yyyz_zzzz_0, tey_yyz_yyzz_0, tey_yyz_yyzz_1, \
                                         tey_yyz_yzz_0, tey_yyz_yzz_1, tey_yyz_yzzz_0, tey_yyz_yzzz_1, tey_yyz_zzz_0, \
                                         tey_yyz_zzz_1, tey_yyz_zzzz_0, tey_yyz_zzzz_1, tey_yyzz_xxxx_0, tey_yyzz_xxxy_0, \
                                         tey_yyzz_xxxz_0, tey_yyzz_xxyy_0, tey_yyzz_xxyz_0, tey_yyzz_xxzz_0, tey_yyzz_xyyy_0, \
                                         tey_yyzz_xyyz_0, tey_yyzz_xyzz_0, tey_yyzz_xzzz_0, tey_yyzz_yyyy_0, tey_yyzz_yyyz_0, \
                                         tey_yyzz_yyzz_0, tey_yz_yyzz_0, tey_yz_yyzz_1, tey_yz_yzzz_0, tey_yz_yzzz_1, \
                                         tey_yz_zzzz_0, tey_yz_zzzz_1, tey_yzz_xxx_0, tey_yzz_xxx_1, tey_yzz_xxxx_0, \
                                         tey_yzz_xxxx_1, tey_yzz_xxxy_0, tey_yzz_xxxy_1, tey_yzz_xxxz_0, tey_yzz_xxxz_1, \
                                         tey_yzz_xxy_0, tey_yzz_xxy_1, tey_yzz_xxyy_0, tey_yzz_xxyy_1, tey_yzz_xxyz_0, \
                                         tey_yzz_xxyz_1, tey_yzz_xxz_0, tey_yzz_xxz_1, tey_yzz_xxzz_0, tey_yzz_xxzz_1, \
                                         tey_yzz_xyy_0, tey_yzz_xyy_1, tey_yzz_xyyy_0, tey_yzz_xyyy_1, tey_yzz_xyyz_0, \
                                         tey_yzz_xyyz_1, tey_yzz_xyz_0, tey_yzz_xyz_1, tey_yzz_xyzz_0, tey_yzz_xyzz_1, \
                                         tey_yzz_xzz_0, tey_yzz_xzz_1, tey_yzz_xzzz_0, tey_yzz_xzzz_1, tey_yzz_yyy_0, \
                                         tey_yzz_yyy_1, tey_yzz_yyyy_0, tey_yzz_yyyy_1, tey_yzz_yyyz_0, tey_yzz_yyyz_1, \
                                         tey_yzz_yyz_0, tey_yzz_yyz_1, tey_yzz_yyzz_0, tey_yzz_yyzz_1, tey_yzz_yzz_0, \
                                         tey_yzz_yzz_1, tey_zz_xxxx_0, tey_zz_xxxx_1, tey_zz_xxxy_0, tey_zz_xxxy_1, \
                                         tey_zz_xxxz_0, tey_zz_xxxz_1, tey_zz_xxyy_0, tey_zz_xxyy_1, tey_zz_xxyz_0, \
                                         tey_zz_xxyz_1, tey_zz_xxzz_0, tey_zz_xxzz_1, tey_zz_xyyy_0, tey_zz_xyyy_1, \
                                         tey_zz_xyyz_0, tey_zz_xyyz_1, tey_zz_xyzz_0, tey_zz_xyzz_1, tey_zz_xzzz_0, \
                                         tey_zz_xzzz_1, tey_zz_yyyy_0, tey_zz_yyyy_1, tey_zz_yyyz_0, tey_zz_yyyz_1, \
                                         tey_zz_yyzz_0, tey_zz_yyzz_1, tez_yyyz_yyzz_0, tez_yyyz_yzzz_0, tez_yyyz_zzzz_0, \
                                         tez_yyz_yyzz_0, tez_yyz_yyzz_1, tez_yyz_yzz_0, tez_yyz_yzz_1, tez_yyz_yzzz_0, \
                                         tez_yyz_yzzz_1, tez_yyz_zzz_0, tez_yyz_zzz_1, tez_yyz_zzzz_0, tez_yyz_zzzz_1, \
                                         tez_yyzz_xxxx_0, tez_yyzz_xxxy_0, tez_yyzz_xxxz_0, tez_yyzz_xxyy_0, tez_yyzz_xxyz_0, \
                                         tez_yyzz_xxzz_0, tez_yyzz_xyyy_0, tez_yyzz_xyyz_0, tez_yyzz_xyzz_0, tez_yyzz_xzzz_0, \
                                         tez_yyzz_yyyy_0, tez_yyzz_yyyz_0, tez_yyzz_yyzz_0, tez_yz_yyzz_0, tez_yz_yyzz_1, \
                                         tez_yz_yzzz_0, tez_yz_yzzz_1, tez_yz_zzzz_0, tez_yz_zzzz_1, tez_yzz_xxx_0, \
                                         tez_yzz_xxx_1, tez_yzz_xxxx_0, tez_yzz_xxxx_1, tez_yzz_xxxy_0, tez_yzz_xxxy_1, \
                                         tez_yzz_xxxz_0, tez_yzz_xxxz_1, tez_yzz_xxy_0, tez_yzz_xxy_1, tez_yzz_xxyy_0, \
                                         tez_yzz_xxyy_1, tez_yzz_xxyz_0, tez_yzz_xxyz_1, tez_yzz_xxz_0, tez_yzz_xxz_1, \
                                         tez_yzz_xxzz_0, tez_yzz_xxzz_1, tez_yzz_xyy_0, tez_yzz_xyy_1, tez_yzz_xyyy_0, \
                                         tez_yzz_xyyy_1, tez_yzz_xyyz_0, tez_yzz_xyyz_1, tez_yzz_xyz_0, tez_yzz_xyz_1, \
                                         tez_yzz_xyzz_0, tez_yzz_xyzz_1, tez_yzz_xzz_0, tez_yzz_xzz_1, tez_yzz_xzzz_0, \
                                         tez_yzz_xzzz_1, tez_yzz_yyy_0, tez_yzz_yyy_1, tez_yzz_yyyy_0, tez_yzz_yyyy_1, \
                                         tez_yzz_yyyz_0, tez_yzz_yyyz_1, tez_yzz_yyz_0, tez_yzz_yyz_1, tez_yzz_yyzz_0, \
                                         tez_yzz_yyzz_1, tez_yzz_yzz_0, tez_yzz_yzz_1, tez_zz_xxxx_0, tez_zz_xxxx_1, \
                                         tez_zz_xxxy_0, tez_zz_xxxy_1, tez_zz_xxxz_0, tez_zz_xxxz_1, tez_zz_xxyy_0, \
                                         tez_zz_xxyy_1, tez_zz_xxyz_0, tez_zz_xxyz_1, tez_zz_xxzz_0, tez_zz_xxzz_1, \
                                         tez_zz_xyyy_0, tez_zz_xyyy_1, tez_zz_xyyz_0, tez_zz_xyyz_1, tez_zz_xyzz_0, \
                                         tez_zz_xyzz_1, tez_zz_xzzz_0, tez_zz_xzzz_1, tez_zz_yyyy_0, tez_zz_yyyy_1, \
                                         tez_zz_yyyz_0, tez_zz_yyyz_1, tez_zz_yyzz_0, tez_zz_yyzz_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_yyyz_yyzz_0[j] = pa_y[j] * tex_yyz_yyzz_0[j] - pc_y[j] * tex_yyz_yyzz_1[j] + fl1_fx * tex_yz_yyzz_0[j] -
                                     fl1_fx * tex_yz_yyzz_1[j] + fl1_fx * tex_yyz_yzz_0[j] - fl1_fx * tex_yyz_yzz_1[j];

                tey_yyyz_yyzz_0[j] = pa_y[j] * tey_yyz_yyzz_0[j] - pc_y[j] * tey_yyz_yyzz_1[j] + fl1_fx * tey_yz_yyzz_0[j] -
                                     fl1_fx * tey_yz_yyzz_1[j] + fl1_fx * tey_yyz_yzz_0[j] - fl1_fx * tey_yyz_yzz_1[j] + ta_yyz_yyzz_1[j];

                tez_yyyz_yyzz_0[j] = pa_y[j] * tez_yyz_yyzz_0[j] - pc_y[j] * tez_yyz_yyzz_1[j] + fl1_fx * tez_yz_yyzz_0[j] -
                                     fl1_fx * tez_yz_yyzz_1[j] + fl1_fx * tez_yyz_yzz_0[j] - fl1_fx * tez_yyz_yzz_1[j];

                tex_yyyz_yzzz_0[j] = pa_y[j] * tex_yyz_yzzz_0[j] - pc_y[j] * tex_yyz_yzzz_1[j] + fl1_fx * tex_yz_yzzz_0[j] -
                                     fl1_fx * tex_yz_yzzz_1[j] + 0.5 * fl1_fx * tex_yyz_zzz_0[j] - 0.5 * fl1_fx * tex_yyz_zzz_1[j];

                tey_yyyz_yzzz_0[j] = pa_y[j] * tey_yyz_yzzz_0[j] - pc_y[j] * tey_yyz_yzzz_1[j] + fl1_fx * tey_yz_yzzz_0[j] -
                                     fl1_fx * tey_yz_yzzz_1[j] + 0.5 * fl1_fx * tey_yyz_zzz_0[j] - 0.5 * fl1_fx * tey_yyz_zzz_1[j] + ta_yyz_yzzz_1[j];

                tez_yyyz_yzzz_0[j] = pa_y[j] * tez_yyz_yzzz_0[j] - pc_y[j] * tez_yyz_yzzz_1[j] + fl1_fx * tez_yz_yzzz_0[j] -
                                     fl1_fx * tez_yz_yzzz_1[j] + 0.5 * fl1_fx * tez_yyz_zzz_0[j] - 0.5 * fl1_fx * tez_yyz_zzz_1[j];

                tex_yyyz_zzzz_0[j] =
                    pa_y[j] * tex_yyz_zzzz_0[j] - pc_y[j] * tex_yyz_zzzz_1[j] + fl1_fx * tex_yz_zzzz_0[j] - fl1_fx * tex_yz_zzzz_1[j];

                tey_yyyz_zzzz_0[j] = pa_y[j] * tey_yyz_zzzz_0[j] - pc_y[j] * tey_yyz_zzzz_1[j] + fl1_fx * tey_yz_zzzz_0[j] -
                                     fl1_fx * tey_yz_zzzz_1[j] + ta_yyz_zzzz_1[j];

                tez_yyyz_zzzz_0[j] =
                    pa_y[j] * tez_yyz_zzzz_0[j] - pc_y[j] * tez_yyz_zzzz_1[j] + fl1_fx * tez_yz_zzzz_0[j] - fl1_fx * tez_yz_zzzz_1[j];

                tex_yyzz_xxxx_0[j] =
                    pa_y[j] * tex_yzz_xxxx_0[j] - pc_y[j] * tex_yzz_xxxx_1[j] + 0.5 * fl1_fx * tex_zz_xxxx_0[j] - 0.5 * fl1_fx * tex_zz_xxxx_1[j];

                tey_yyzz_xxxx_0[j] = pa_y[j] * tey_yzz_xxxx_0[j] - pc_y[j] * tey_yzz_xxxx_1[j] + 0.5 * fl1_fx * tey_zz_xxxx_0[j] -
                                     0.5 * fl1_fx * tey_zz_xxxx_1[j] + ta_yzz_xxxx_1[j];

                tez_yyzz_xxxx_0[j] =
                    pa_y[j] * tez_yzz_xxxx_0[j] - pc_y[j] * tez_yzz_xxxx_1[j] + 0.5 * fl1_fx * tez_zz_xxxx_0[j] - 0.5 * fl1_fx * tez_zz_xxxx_1[j];

                tex_yyzz_xxxy_0[j] = pa_y[j] * tex_yzz_xxxy_0[j] - pc_y[j] * tex_yzz_xxxy_1[j] + 0.5 * fl1_fx * tex_zz_xxxy_0[j] -
                                     0.5 * fl1_fx * tex_zz_xxxy_1[j] + 0.5 * fl1_fx * tex_yzz_xxx_0[j] - 0.5 * fl1_fx * tex_yzz_xxx_1[j];

                tey_yyzz_xxxy_0[j] = pa_y[j] * tey_yzz_xxxy_0[j] - pc_y[j] * tey_yzz_xxxy_1[j] + 0.5 * fl1_fx * tey_zz_xxxy_0[j] -
                                     0.5 * fl1_fx * tey_zz_xxxy_1[j] + 0.5 * fl1_fx * tey_yzz_xxx_0[j] - 0.5 * fl1_fx * tey_yzz_xxx_1[j] +
                                     ta_yzz_xxxy_1[j];

                tez_yyzz_xxxy_0[j] = pa_y[j] * tez_yzz_xxxy_0[j] - pc_y[j] * tez_yzz_xxxy_1[j] + 0.5 * fl1_fx * tez_zz_xxxy_0[j] -
                                     0.5 * fl1_fx * tez_zz_xxxy_1[j] + 0.5 * fl1_fx * tez_yzz_xxx_0[j] - 0.5 * fl1_fx * tez_yzz_xxx_1[j];

                tex_yyzz_xxxz_0[j] =
                    pa_y[j] * tex_yzz_xxxz_0[j] - pc_y[j] * tex_yzz_xxxz_1[j] + 0.5 * fl1_fx * tex_zz_xxxz_0[j] - 0.5 * fl1_fx * tex_zz_xxxz_1[j];

                tey_yyzz_xxxz_0[j] = pa_y[j] * tey_yzz_xxxz_0[j] - pc_y[j] * tey_yzz_xxxz_1[j] + 0.5 * fl1_fx * tey_zz_xxxz_0[j] -
                                     0.5 * fl1_fx * tey_zz_xxxz_1[j] + ta_yzz_xxxz_1[j];

                tez_yyzz_xxxz_0[j] =
                    pa_y[j] * tez_yzz_xxxz_0[j] - pc_y[j] * tez_yzz_xxxz_1[j] + 0.5 * fl1_fx * tez_zz_xxxz_0[j] - 0.5 * fl1_fx * tez_zz_xxxz_1[j];

                tex_yyzz_xxyy_0[j] = pa_y[j] * tex_yzz_xxyy_0[j] - pc_y[j] * tex_yzz_xxyy_1[j] + 0.5 * fl1_fx * tex_zz_xxyy_0[j] -
                                     0.5 * fl1_fx * tex_zz_xxyy_1[j] + fl1_fx * tex_yzz_xxy_0[j] - fl1_fx * tex_yzz_xxy_1[j];

                tey_yyzz_xxyy_0[j] = pa_y[j] * tey_yzz_xxyy_0[j] - pc_y[j] * tey_yzz_xxyy_1[j] + 0.5 * fl1_fx * tey_zz_xxyy_0[j] -
                                     0.5 * fl1_fx * tey_zz_xxyy_1[j] + fl1_fx * tey_yzz_xxy_0[j] - fl1_fx * tey_yzz_xxy_1[j] + ta_yzz_xxyy_1[j];

                tez_yyzz_xxyy_0[j] = pa_y[j] * tez_yzz_xxyy_0[j] - pc_y[j] * tez_yzz_xxyy_1[j] + 0.5 * fl1_fx * tez_zz_xxyy_0[j] -
                                     0.5 * fl1_fx * tez_zz_xxyy_1[j] + fl1_fx * tez_yzz_xxy_0[j] - fl1_fx * tez_yzz_xxy_1[j];

                tex_yyzz_xxyz_0[j] = pa_y[j] * tex_yzz_xxyz_0[j] - pc_y[j] * tex_yzz_xxyz_1[j] + 0.5 * fl1_fx * tex_zz_xxyz_0[j] -
                                     0.5 * fl1_fx * tex_zz_xxyz_1[j] + 0.5 * fl1_fx * tex_yzz_xxz_0[j] - 0.5 * fl1_fx * tex_yzz_xxz_1[j];

                tey_yyzz_xxyz_0[j] = pa_y[j] * tey_yzz_xxyz_0[j] - pc_y[j] * tey_yzz_xxyz_1[j] + 0.5 * fl1_fx * tey_zz_xxyz_0[j] -
                                     0.5 * fl1_fx * tey_zz_xxyz_1[j] + 0.5 * fl1_fx * tey_yzz_xxz_0[j] - 0.5 * fl1_fx * tey_yzz_xxz_1[j] +
                                     ta_yzz_xxyz_1[j];

                tez_yyzz_xxyz_0[j] = pa_y[j] * tez_yzz_xxyz_0[j] - pc_y[j] * tez_yzz_xxyz_1[j] + 0.5 * fl1_fx * tez_zz_xxyz_0[j] -
                                     0.5 * fl1_fx * tez_zz_xxyz_1[j] + 0.5 * fl1_fx * tez_yzz_xxz_0[j] - 0.5 * fl1_fx * tez_yzz_xxz_1[j];

                tex_yyzz_xxzz_0[j] =
                    pa_y[j] * tex_yzz_xxzz_0[j] - pc_y[j] * tex_yzz_xxzz_1[j] + 0.5 * fl1_fx * tex_zz_xxzz_0[j] - 0.5 * fl1_fx * tex_zz_xxzz_1[j];

                tey_yyzz_xxzz_0[j] = pa_y[j] * tey_yzz_xxzz_0[j] - pc_y[j] * tey_yzz_xxzz_1[j] + 0.5 * fl1_fx * tey_zz_xxzz_0[j] -
                                     0.5 * fl1_fx * tey_zz_xxzz_1[j] + ta_yzz_xxzz_1[j];

                tez_yyzz_xxzz_0[j] =
                    pa_y[j] * tez_yzz_xxzz_0[j] - pc_y[j] * tez_yzz_xxzz_1[j] + 0.5 * fl1_fx * tez_zz_xxzz_0[j] - 0.5 * fl1_fx * tez_zz_xxzz_1[j];

                tex_yyzz_xyyy_0[j] = pa_y[j] * tex_yzz_xyyy_0[j] - pc_y[j] * tex_yzz_xyyy_1[j] + 0.5 * fl1_fx * tex_zz_xyyy_0[j] -
                                     0.5 * fl1_fx * tex_zz_xyyy_1[j] + 1.5 * fl1_fx * tex_yzz_xyy_0[j] - 1.5 * fl1_fx * tex_yzz_xyy_1[j];

                tey_yyzz_xyyy_0[j] = pa_y[j] * tey_yzz_xyyy_0[j] - pc_y[j] * tey_yzz_xyyy_1[j] + 0.5 * fl1_fx * tey_zz_xyyy_0[j] -
                                     0.5 * fl1_fx * tey_zz_xyyy_1[j] + 1.5 * fl1_fx * tey_yzz_xyy_0[j] - 1.5 * fl1_fx * tey_yzz_xyy_1[j] +
                                     ta_yzz_xyyy_1[j];

                tez_yyzz_xyyy_0[j] = pa_y[j] * tez_yzz_xyyy_0[j] - pc_y[j] * tez_yzz_xyyy_1[j] + 0.5 * fl1_fx * tez_zz_xyyy_0[j] -
                                     0.5 * fl1_fx * tez_zz_xyyy_1[j] + 1.5 * fl1_fx * tez_yzz_xyy_0[j] - 1.5 * fl1_fx * tez_yzz_xyy_1[j];

                tex_yyzz_xyyz_0[j] = pa_y[j] * tex_yzz_xyyz_0[j] - pc_y[j] * tex_yzz_xyyz_1[j] + 0.5 * fl1_fx * tex_zz_xyyz_0[j] -
                                     0.5 * fl1_fx * tex_zz_xyyz_1[j] + fl1_fx * tex_yzz_xyz_0[j] - fl1_fx * tex_yzz_xyz_1[j];

                tey_yyzz_xyyz_0[j] = pa_y[j] * tey_yzz_xyyz_0[j] - pc_y[j] * tey_yzz_xyyz_1[j] + 0.5 * fl1_fx * tey_zz_xyyz_0[j] -
                                     0.5 * fl1_fx * tey_zz_xyyz_1[j] + fl1_fx * tey_yzz_xyz_0[j] - fl1_fx * tey_yzz_xyz_1[j] + ta_yzz_xyyz_1[j];

                tez_yyzz_xyyz_0[j] = pa_y[j] * tez_yzz_xyyz_0[j] - pc_y[j] * tez_yzz_xyyz_1[j] + 0.5 * fl1_fx * tez_zz_xyyz_0[j] -
                                     0.5 * fl1_fx * tez_zz_xyyz_1[j] + fl1_fx * tez_yzz_xyz_0[j] - fl1_fx * tez_yzz_xyz_1[j];

                tex_yyzz_xyzz_0[j] = pa_y[j] * tex_yzz_xyzz_0[j] - pc_y[j] * tex_yzz_xyzz_1[j] + 0.5 * fl1_fx * tex_zz_xyzz_0[j] -
                                     0.5 * fl1_fx * tex_zz_xyzz_1[j] + 0.5 * fl1_fx * tex_yzz_xzz_0[j] - 0.5 * fl1_fx * tex_yzz_xzz_1[j];

                tey_yyzz_xyzz_0[j] = pa_y[j] * tey_yzz_xyzz_0[j] - pc_y[j] * tey_yzz_xyzz_1[j] + 0.5 * fl1_fx * tey_zz_xyzz_0[j] -
                                     0.5 * fl1_fx * tey_zz_xyzz_1[j] + 0.5 * fl1_fx * tey_yzz_xzz_0[j] - 0.5 * fl1_fx * tey_yzz_xzz_1[j] +
                                     ta_yzz_xyzz_1[j];

                tez_yyzz_xyzz_0[j] = pa_y[j] * tez_yzz_xyzz_0[j] - pc_y[j] * tez_yzz_xyzz_1[j] + 0.5 * fl1_fx * tez_zz_xyzz_0[j] -
                                     0.5 * fl1_fx * tez_zz_xyzz_1[j] + 0.5 * fl1_fx * tez_yzz_xzz_0[j] - 0.5 * fl1_fx * tez_yzz_xzz_1[j];

                tex_yyzz_xzzz_0[j] =
                    pa_y[j] * tex_yzz_xzzz_0[j] - pc_y[j] * tex_yzz_xzzz_1[j] + 0.5 * fl1_fx * tex_zz_xzzz_0[j] - 0.5 * fl1_fx * tex_zz_xzzz_1[j];

                tey_yyzz_xzzz_0[j] = pa_y[j] * tey_yzz_xzzz_0[j] - pc_y[j] * tey_yzz_xzzz_1[j] + 0.5 * fl1_fx * tey_zz_xzzz_0[j] -
                                     0.5 * fl1_fx * tey_zz_xzzz_1[j] + ta_yzz_xzzz_1[j];

                tez_yyzz_xzzz_0[j] =
                    pa_y[j] * tez_yzz_xzzz_0[j] - pc_y[j] * tez_yzz_xzzz_1[j] + 0.5 * fl1_fx * tez_zz_xzzz_0[j] - 0.5 * fl1_fx * tez_zz_xzzz_1[j];

                tex_yyzz_yyyy_0[j] = pa_y[j] * tex_yzz_yyyy_0[j] - pc_y[j] * tex_yzz_yyyy_1[j] + 0.5 * fl1_fx * tex_zz_yyyy_0[j] -
                                     0.5 * fl1_fx * tex_zz_yyyy_1[j] + 2.0 * fl1_fx * tex_yzz_yyy_0[j] - 2.0 * fl1_fx * tex_yzz_yyy_1[j];

                tey_yyzz_yyyy_0[j] = pa_y[j] * tey_yzz_yyyy_0[j] - pc_y[j] * tey_yzz_yyyy_1[j] + 0.5 * fl1_fx * tey_zz_yyyy_0[j] -
                                     0.5 * fl1_fx * tey_zz_yyyy_1[j] + 2.0 * fl1_fx * tey_yzz_yyy_0[j] - 2.0 * fl1_fx * tey_yzz_yyy_1[j] +
                                     ta_yzz_yyyy_1[j];

                tez_yyzz_yyyy_0[j] = pa_y[j] * tez_yzz_yyyy_0[j] - pc_y[j] * tez_yzz_yyyy_1[j] + 0.5 * fl1_fx * tez_zz_yyyy_0[j] -
                                     0.5 * fl1_fx * tez_zz_yyyy_1[j] + 2.0 * fl1_fx * tez_yzz_yyy_0[j] - 2.0 * fl1_fx * tez_yzz_yyy_1[j];

                tex_yyzz_yyyz_0[j] = pa_y[j] * tex_yzz_yyyz_0[j] - pc_y[j] * tex_yzz_yyyz_1[j] + 0.5 * fl1_fx * tex_zz_yyyz_0[j] -
                                     0.5 * fl1_fx * tex_zz_yyyz_1[j] + 1.5 * fl1_fx * tex_yzz_yyz_0[j] - 1.5 * fl1_fx * tex_yzz_yyz_1[j];

                tey_yyzz_yyyz_0[j] = pa_y[j] * tey_yzz_yyyz_0[j] - pc_y[j] * tey_yzz_yyyz_1[j] + 0.5 * fl1_fx * tey_zz_yyyz_0[j] -
                                     0.5 * fl1_fx * tey_zz_yyyz_1[j] + 1.5 * fl1_fx * tey_yzz_yyz_0[j] - 1.5 * fl1_fx * tey_yzz_yyz_1[j] +
                                     ta_yzz_yyyz_1[j];

                tez_yyzz_yyyz_0[j] = pa_y[j] * tez_yzz_yyyz_0[j] - pc_y[j] * tez_yzz_yyyz_1[j] + 0.5 * fl1_fx * tez_zz_yyyz_0[j] -
                                     0.5 * fl1_fx * tez_zz_yyyz_1[j] + 1.5 * fl1_fx * tez_yzz_yyz_0[j] - 1.5 * fl1_fx * tez_yzz_yyz_1[j];

                tex_yyzz_yyzz_0[j] = pa_y[j] * tex_yzz_yyzz_0[j] - pc_y[j] * tex_yzz_yyzz_1[j] + 0.5 * fl1_fx * tex_zz_yyzz_0[j] -
                                     0.5 * fl1_fx * tex_zz_yyzz_1[j] + fl1_fx * tex_yzz_yzz_0[j] - fl1_fx * tex_yzz_yzz_1[j];

                tey_yyzz_yyzz_0[j] = pa_y[j] * tey_yzz_yyzz_0[j] - pc_y[j] * tey_yzz_yyzz_1[j] + 0.5 * fl1_fx * tey_zz_yyzz_0[j] -
                                     0.5 * fl1_fx * tey_zz_yyzz_1[j] + fl1_fx * tey_yzz_yzz_0[j] - fl1_fx * tey_yzz_yzz_1[j] + ta_yzz_yyzz_1[j];

                tez_yyzz_yyzz_0[j] = pa_y[j] * tez_yzz_yyzz_0[j] - pc_y[j] * tez_yzz_yyzz_1[j] + 0.5 * fl1_fx * tez_zz_yyzz_0[j] -
                                     0.5 * fl1_fx * tez_zz_yyzz_1[j] + fl1_fx * tez_yzz_yzz_0[j] - fl1_fx * tez_yzz_yzz_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForGG_579_627(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pcDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto)
{
    // Batch of Integrals (579,627)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_4_4_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_2_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(3 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            // set up pointers to tensors product of distances R(PC) = P - C

            auto pc_y = pcDistances.data(3 * idx + 1);

            // set up pointers to auxilary integrals

            auto tex_yzz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 133);

            auto tey_yzz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 133);

            auto tez_yzz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 133);

            auto tex_yzz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 134);

            auto tey_yzz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 134);

            auto tez_yzz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 134);

            auto tex_zzz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 135);

            auto tey_zzz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 135);

            auto tez_zzz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 135);

            auto tex_zzz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 136);

            auto tey_zzz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 136);

            auto tez_zzz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 136);

            auto tex_zzz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 137);

            auto tey_zzz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 137);

            auto tez_zzz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 137);

            auto tex_zzz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 138);

            auto tey_zzz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 138);

            auto tez_zzz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 138);

            auto tex_zzz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 139);

            auto tey_zzz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 139);

            auto tez_zzz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 139);

            auto tex_zzz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 140);

            auto tey_zzz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 140);

            auto tez_zzz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 140);

            auto tex_zzz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 141);

            auto tey_zzz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 141);

            auto tez_zzz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 141);

            auto tex_zzz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 142);

            auto tey_zzz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 142);

            auto tez_zzz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 142);

            auto tex_zzz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 143);

            auto tey_zzz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 143);

            auto tez_zzz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 143);

            auto tex_zzz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 144);

            auto tey_zzz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 144);

            auto tez_zzz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 144);

            auto tex_zzz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 145);

            auto tey_zzz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 145);

            auto tez_zzz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 145);

            auto tex_zzz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 146);

            auto tey_zzz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 146);

            auto tez_zzz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 146);

            auto tex_zzz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 147);

            auto tey_zzz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 147);

            auto tez_zzz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 147);

            auto tex_zzz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 148);

            auto tey_zzz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 148);

            auto tez_zzz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 148);

            auto tex_yzz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 133);

            auto tey_yzz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 133);

            auto tez_yzz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 133);

            auto tex_yzz_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 134);

            auto tey_yzz_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 134);

            auto tez_yzz_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 134);

            auto tex_zzz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 135);

            auto tey_zzz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 135);

            auto tez_zzz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 135);

            auto tex_zzz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 136);

            auto tey_zzz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 136);

            auto tez_zzz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 136);

            auto tex_zzz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 137);

            auto tey_zzz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 137);

            auto tez_zzz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 137);

            auto tex_zzz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 138);

            auto tey_zzz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 138);

            auto tez_zzz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 138);

            auto tex_zzz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 139);

            auto tey_zzz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 139);

            auto tez_zzz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 139);

            auto tex_zzz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 140);

            auto tey_zzz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 140);

            auto tez_zzz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 140);

            auto tex_zzz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 141);

            auto tey_zzz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 141);

            auto tez_zzz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 141);

            auto tex_zzz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 142);

            auto tey_zzz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 142);

            auto tez_zzz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 142);

            auto tex_zzz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 143);

            auto tey_zzz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 143);

            auto tez_zzz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 143);

            auto tex_zzz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 144);

            auto tey_zzz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 144);

            auto tez_zzz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 144);

            auto tex_zzz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 145);

            auto tey_zzz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 145);

            auto tez_zzz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 145);

            auto tex_zzz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 146);

            auto tey_zzz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 146);

            auto tez_zzz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 146);

            auto tex_zzz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 147);

            auto tey_zzz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 147);

            auto tez_zzz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 147);

            auto tex_zzz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 148);

            auto tey_zzz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 148);

            auto tez_zzz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 148);

            auto tex_zz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 88);

            auto tey_zz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 88);

            auto tez_zz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 88);

            auto tex_zz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 89);

            auto tey_zz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 89);

            auto tez_zz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 89);

            auto tex_zz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 88);

            auto tey_zz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 88);

            auto tez_zz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 88);

            auto tex_zz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 89);

            auto tey_zz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 89);

            auto tez_zz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 89);

            auto tex_yzz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 89);

            auto tey_yzz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 89);

            auto tez_yzz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 89);

            auto tex_zzz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 90);

            auto tey_zzz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 90);

            auto tez_zzz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 90);

            auto tex_zzz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 91);

            auto tey_zzz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 91);

            auto tez_zzz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 91);

            auto tex_zzz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 92);

            auto tey_zzz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 92);

            auto tez_zzz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 92);

            auto tex_zzz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 93);

            auto tey_zzz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 93);

            auto tez_zzz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 93);

            auto tex_zzz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 94);

            auto tey_zzz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 94);

            auto tez_zzz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 94);

            auto tex_zzz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 95);

            auto tey_zzz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 95);

            auto tez_zzz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 95);

            auto tex_zzz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 96);

            auto tey_zzz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 96);

            auto tez_zzz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 96);

            auto tex_zzz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 97);

            auto tey_zzz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 97);

            auto tez_zzz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 97);

            auto tex_zzz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 98);

            auto tey_zzz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 98);

            auto tez_zzz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 98);

            auto tex_zzz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 99);

            auto tey_zzz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 99);

            auto tez_zzz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 99);

            auto tex_yzz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 89);

            auto tey_yzz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 89);

            auto tez_yzz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 89);

            auto tex_zzz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 90);

            auto tey_zzz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 90);

            auto tez_zzz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 90);

            auto tex_zzz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 91);

            auto tey_zzz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 91);

            auto tez_zzz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 91);

            auto tex_zzz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 92);

            auto tey_zzz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 92);

            auto tez_zzz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 92);

            auto tex_zzz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 93);

            auto tey_zzz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 93);

            auto tez_zzz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 93);

            auto tex_zzz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 94);

            auto tey_zzz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 94);

            auto tez_zzz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 94);

            auto tex_zzz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 95);

            auto tey_zzz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 95);

            auto tez_zzz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 95);

            auto tex_zzz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 96);

            auto tey_zzz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 96);

            auto tez_zzz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 96);

            auto tex_zzz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 97);

            auto tey_zzz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 97);

            auto tez_zzz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 97);

            auto tex_zzz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 98);

            auto tey_zzz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 98);

            auto tez_zzz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 98);

            auto tex_zzz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 99);

            auto tey_zzz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 99);

            auto tez_zzz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 99);

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

            // set up pointers to integrals

            auto tex_yyzz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 193);

            auto tey_yyzz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 193);

            auto tez_yyzz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 193);

            auto tex_yyzz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 194);

            auto tey_yyzz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 194);

            auto tez_yyzz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 194);

            auto tex_yzzz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 195);

            auto tey_yzzz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 195);

            auto tez_yzzz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 195);

            auto tex_yzzz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 196);

            auto tey_yzzz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 196);

            auto tez_yzzz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 196);

            auto tex_yzzz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 197);

            auto tey_yzzz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 197);

            auto tez_yzzz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 197);

            auto tex_yzzz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 198);

            auto tey_yzzz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 198);

            auto tez_yzzz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 198);

            auto tex_yzzz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 199);

            auto tey_yzzz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 199);

            auto tez_yzzz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 199);

            auto tex_yzzz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 200);

            auto tey_yzzz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 200);

            auto tez_yzzz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 200);

            auto tex_yzzz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 201);

            auto tey_yzzz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 201);

            auto tez_yzzz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 201);

            auto tex_yzzz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 202);

            auto tey_yzzz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 202);

            auto tez_yzzz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 202);

            auto tex_yzzz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 203);

            auto tey_yzzz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 203);

            auto tez_yzzz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 203);

            auto tex_yzzz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 204);

            auto tey_yzzz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 204);

            auto tez_yzzz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 204);

            auto tex_yzzz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 205);

            auto tey_yzzz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 205);

            auto tez_yzzz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 205);

            auto tex_yzzz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 206);

            auto tey_yzzz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 206);

            auto tez_yzzz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 206);

            auto tex_yzzz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 207);

            auto tey_yzzz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 207);

            auto tez_yzzz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 207);

            auto tex_yzzz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 208);

            auto tey_yzzz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 208);

            auto tez_yzzz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 208);

            // Batch of Integrals (579,627)

            #pragma omp simd aligned(fx, pa_y, pc_y, ta_yzz_yzzz_1, ta_yzz_zzzz_1, ta_zzz_xxxx_1, \
                                         ta_zzz_xxxy_1, ta_zzz_xxxz_1, ta_zzz_xxyy_1, ta_zzz_xxyz_1, ta_zzz_xxzz_1, \
                                         ta_zzz_xyyy_1, ta_zzz_xyyz_1, ta_zzz_xyzz_1, ta_zzz_xzzz_1, ta_zzz_yyyy_1, \
                                         ta_zzz_yyyz_1, ta_zzz_yyzz_1, ta_zzz_yzzz_1, tex_yyzz_yzzz_0, tex_yyzz_zzzz_0, \
                                         tex_yzz_yzzz_0, tex_yzz_yzzz_1, tex_yzz_zzz_0, tex_yzz_zzz_1, tex_yzz_zzzz_0, \
                                         tex_yzz_zzzz_1, tex_yzzz_xxxx_0, tex_yzzz_xxxy_0, tex_yzzz_xxxz_0, tex_yzzz_xxyy_0, \
                                         tex_yzzz_xxyz_0, tex_yzzz_xxzz_0, tex_yzzz_xyyy_0, tex_yzzz_xyyz_0, tex_yzzz_xyzz_0, \
                                         tex_yzzz_xzzz_0, tex_yzzz_yyyy_0, tex_yzzz_yyyz_0, tex_yzzz_yyzz_0, tex_yzzz_yzzz_0, \
                                         tex_zz_yzzz_0, tex_zz_yzzz_1, tex_zz_zzzz_0, tex_zz_zzzz_1, tex_zzz_xxx_0, \
                                         tex_zzz_xxx_1, tex_zzz_xxxx_0, tex_zzz_xxxx_1, tex_zzz_xxxy_0, tex_zzz_xxxy_1, \
                                         tex_zzz_xxxz_0, tex_zzz_xxxz_1, tex_zzz_xxy_0, tex_zzz_xxy_1, tex_zzz_xxyy_0, \
                                         tex_zzz_xxyy_1, tex_zzz_xxyz_0, tex_zzz_xxyz_1, tex_zzz_xxz_0, tex_zzz_xxz_1, \
                                         tex_zzz_xxzz_0, tex_zzz_xxzz_1, tex_zzz_xyy_0, tex_zzz_xyy_1, tex_zzz_xyyy_0, \
                                         tex_zzz_xyyy_1, tex_zzz_xyyz_0, tex_zzz_xyyz_1, tex_zzz_xyz_0, tex_zzz_xyz_1, \
                                         tex_zzz_xyzz_0, tex_zzz_xyzz_1, tex_zzz_xzz_0, tex_zzz_xzz_1, tex_zzz_xzzz_0, \
                                         tex_zzz_xzzz_1, tex_zzz_yyy_0, tex_zzz_yyy_1, tex_zzz_yyyy_0, tex_zzz_yyyy_1, \
                                         tex_zzz_yyyz_0, tex_zzz_yyyz_1, tex_zzz_yyz_0, tex_zzz_yyz_1, tex_zzz_yyzz_0, \
                                         tex_zzz_yyzz_1, tex_zzz_yzz_0, tex_zzz_yzz_1, tex_zzz_yzzz_0, tex_zzz_yzzz_1, \
                                         tex_zzz_zzz_0, tex_zzz_zzz_1, tey_yyzz_yzzz_0, tey_yyzz_zzzz_0, tey_yzz_yzzz_0, \
                                         tey_yzz_yzzz_1, tey_yzz_zzz_0, tey_yzz_zzz_1, tey_yzz_zzzz_0, tey_yzz_zzzz_1, \
                                         tey_yzzz_xxxx_0, tey_yzzz_xxxy_0, tey_yzzz_xxxz_0, tey_yzzz_xxyy_0, tey_yzzz_xxyz_0, \
                                         tey_yzzz_xxzz_0, tey_yzzz_xyyy_0, tey_yzzz_xyyz_0, tey_yzzz_xyzz_0, tey_yzzz_xzzz_0, \
                                         tey_yzzz_yyyy_0, tey_yzzz_yyyz_0, tey_yzzz_yyzz_0, tey_yzzz_yzzz_0, tey_zz_yzzz_0, \
                                         tey_zz_yzzz_1, tey_zz_zzzz_0, tey_zz_zzzz_1, tey_zzz_xxx_0, tey_zzz_xxx_1, \
                                         tey_zzz_xxxx_0, tey_zzz_xxxx_1, tey_zzz_xxxy_0, tey_zzz_xxxy_1, tey_zzz_xxxz_0, \
                                         tey_zzz_xxxz_1, tey_zzz_xxy_0, tey_zzz_xxy_1, tey_zzz_xxyy_0, tey_zzz_xxyy_1, \
                                         tey_zzz_xxyz_0, tey_zzz_xxyz_1, tey_zzz_xxz_0, tey_zzz_xxz_1, tey_zzz_xxzz_0, \
                                         tey_zzz_xxzz_1, tey_zzz_xyy_0, tey_zzz_xyy_1, tey_zzz_xyyy_0, tey_zzz_xyyy_1, \
                                         tey_zzz_xyyz_0, tey_zzz_xyyz_1, tey_zzz_xyz_0, tey_zzz_xyz_1, tey_zzz_xyzz_0, \
                                         tey_zzz_xyzz_1, tey_zzz_xzz_0, tey_zzz_xzz_1, tey_zzz_xzzz_0, tey_zzz_xzzz_1, \
                                         tey_zzz_yyy_0, tey_zzz_yyy_1, tey_zzz_yyyy_0, tey_zzz_yyyy_1, tey_zzz_yyyz_0, \
                                         tey_zzz_yyyz_1, tey_zzz_yyz_0, tey_zzz_yyz_1, tey_zzz_yyzz_0, tey_zzz_yyzz_1, \
                                         tey_zzz_yzz_0, tey_zzz_yzz_1, tey_zzz_yzzz_0, tey_zzz_yzzz_1, tey_zzz_zzz_0, \
                                         tey_zzz_zzz_1, tez_yyzz_yzzz_0, tez_yyzz_zzzz_0, tez_yzz_yzzz_0, tez_yzz_yzzz_1, \
                                         tez_yzz_zzz_0, tez_yzz_zzz_1, tez_yzz_zzzz_0, tez_yzz_zzzz_1, tez_yzzz_xxxx_0, \
                                         tez_yzzz_xxxy_0, tez_yzzz_xxxz_0, tez_yzzz_xxyy_0, tez_yzzz_xxyz_0, tez_yzzz_xxzz_0, \
                                         tez_yzzz_xyyy_0, tez_yzzz_xyyz_0, tez_yzzz_xyzz_0, tez_yzzz_xzzz_0, tez_yzzz_yyyy_0, \
                                         tez_yzzz_yyyz_0, tez_yzzz_yyzz_0, tez_yzzz_yzzz_0, tez_zz_yzzz_0, tez_zz_yzzz_1, \
                                         tez_zz_zzzz_0, tez_zz_zzzz_1, tez_zzz_xxx_0, tez_zzz_xxx_1, tez_zzz_xxxx_0, \
                                         tez_zzz_xxxx_1, tez_zzz_xxxy_0, tez_zzz_xxxy_1, tez_zzz_xxxz_0, tez_zzz_xxxz_1, \
                                         tez_zzz_xxy_0, tez_zzz_xxy_1, tez_zzz_xxyy_0, tez_zzz_xxyy_1, tez_zzz_xxyz_0, \
                                         tez_zzz_xxyz_1, tez_zzz_xxz_0, tez_zzz_xxz_1, tez_zzz_xxzz_0, tez_zzz_xxzz_1, \
                                         tez_zzz_xyy_0, tez_zzz_xyy_1, tez_zzz_xyyy_0, tez_zzz_xyyy_1, tez_zzz_xyyz_0, \
                                         tez_zzz_xyyz_1, tez_zzz_xyz_0, tez_zzz_xyz_1, tez_zzz_xyzz_0, tez_zzz_xyzz_1, \
                                         tez_zzz_xzz_0, tez_zzz_xzz_1, tez_zzz_xzzz_0, tez_zzz_xzzz_1, tez_zzz_yyy_0, \
                                         tez_zzz_yyy_1, tez_zzz_yyyy_0, tez_zzz_yyyy_1, tez_zzz_yyyz_0, tez_zzz_yyyz_1, \
                                         tez_zzz_yyz_0, tez_zzz_yyz_1, tez_zzz_yyzz_0, tez_zzz_yyzz_1, tez_zzz_yzz_0, \
                                         tez_zzz_yzz_1, tez_zzz_yzzz_0, tez_zzz_yzzz_1, tez_zzz_zzz_0, tez_zzz_zzz_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_yyzz_yzzz_0[j] = pa_y[j] * tex_yzz_yzzz_0[j] - pc_y[j] * tex_yzz_yzzz_1[j] + 0.5 * fl1_fx * tex_zz_yzzz_0[j] -
                                     0.5 * fl1_fx * tex_zz_yzzz_1[j] + 0.5 * fl1_fx * tex_yzz_zzz_0[j] - 0.5 * fl1_fx * tex_yzz_zzz_1[j];

                tey_yyzz_yzzz_0[j] = pa_y[j] * tey_yzz_yzzz_0[j] - pc_y[j] * tey_yzz_yzzz_1[j] + 0.5 * fl1_fx * tey_zz_yzzz_0[j] -
                                     0.5 * fl1_fx * tey_zz_yzzz_1[j] + 0.5 * fl1_fx * tey_yzz_zzz_0[j] - 0.5 * fl1_fx * tey_yzz_zzz_1[j] +
                                     ta_yzz_yzzz_1[j];

                tez_yyzz_yzzz_0[j] = pa_y[j] * tez_yzz_yzzz_0[j] - pc_y[j] * tez_yzz_yzzz_1[j] + 0.5 * fl1_fx * tez_zz_yzzz_0[j] -
                                     0.5 * fl1_fx * tez_zz_yzzz_1[j] + 0.5 * fl1_fx * tez_yzz_zzz_0[j] - 0.5 * fl1_fx * tez_yzz_zzz_1[j];

                tex_yyzz_zzzz_0[j] =
                    pa_y[j] * tex_yzz_zzzz_0[j] - pc_y[j] * tex_yzz_zzzz_1[j] + 0.5 * fl1_fx * tex_zz_zzzz_0[j] - 0.5 * fl1_fx * tex_zz_zzzz_1[j];

                tey_yyzz_zzzz_0[j] = pa_y[j] * tey_yzz_zzzz_0[j] - pc_y[j] * tey_yzz_zzzz_1[j] + 0.5 * fl1_fx * tey_zz_zzzz_0[j] -
                                     0.5 * fl1_fx * tey_zz_zzzz_1[j] + ta_yzz_zzzz_1[j];

                tez_yyzz_zzzz_0[j] =
                    pa_y[j] * tez_yzz_zzzz_0[j] - pc_y[j] * tez_yzz_zzzz_1[j] + 0.5 * fl1_fx * tez_zz_zzzz_0[j] - 0.5 * fl1_fx * tez_zz_zzzz_1[j];

                tex_yzzz_xxxx_0[j] = pa_y[j] * tex_zzz_xxxx_0[j] - pc_y[j] * tex_zzz_xxxx_1[j];

                tey_yzzz_xxxx_0[j] = pa_y[j] * tey_zzz_xxxx_0[j] - pc_y[j] * tey_zzz_xxxx_1[j] + ta_zzz_xxxx_1[j];

                tez_yzzz_xxxx_0[j] = pa_y[j] * tez_zzz_xxxx_0[j] - pc_y[j] * tez_zzz_xxxx_1[j];

                tex_yzzz_xxxy_0[j] =
                    pa_y[j] * tex_zzz_xxxy_0[j] - pc_y[j] * tex_zzz_xxxy_1[j] + 0.5 * fl1_fx * tex_zzz_xxx_0[j] - 0.5 * fl1_fx * tex_zzz_xxx_1[j];

                tey_yzzz_xxxy_0[j] = pa_y[j] * tey_zzz_xxxy_0[j] - pc_y[j] * tey_zzz_xxxy_1[j] + 0.5 * fl1_fx * tey_zzz_xxx_0[j] -
                                     0.5 * fl1_fx * tey_zzz_xxx_1[j] + ta_zzz_xxxy_1[j];

                tez_yzzz_xxxy_0[j] =
                    pa_y[j] * tez_zzz_xxxy_0[j] - pc_y[j] * tez_zzz_xxxy_1[j] + 0.5 * fl1_fx * tez_zzz_xxx_0[j] - 0.5 * fl1_fx * tez_zzz_xxx_1[j];

                tex_yzzz_xxxz_0[j] = pa_y[j] * tex_zzz_xxxz_0[j] - pc_y[j] * tex_zzz_xxxz_1[j];

                tey_yzzz_xxxz_0[j] = pa_y[j] * tey_zzz_xxxz_0[j] - pc_y[j] * tey_zzz_xxxz_1[j] + ta_zzz_xxxz_1[j];

                tez_yzzz_xxxz_0[j] = pa_y[j] * tez_zzz_xxxz_0[j] - pc_y[j] * tez_zzz_xxxz_1[j];

                tex_yzzz_xxyy_0[j] =
                    pa_y[j] * tex_zzz_xxyy_0[j] - pc_y[j] * tex_zzz_xxyy_1[j] + fl1_fx * tex_zzz_xxy_0[j] - fl1_fx * tex_zzz_xxy_1[j];

                tey_yzzz_xxyy_0[j] = pa_y[j] * tey_zzz_xxyy_0[j] - pc_y[j] * tey_zzz_xxyy_1[j] + fl1_fx * tey_zzz_xxy_0[j] -
                                     fl1_fx * tey_zzz_xxy_1[j] + ta_zzz_xxyy_1[j];

                tez_yzzz_xxyy_0[j] =
                    pa_y[j] * tez_zzz_xxyy_0[j] - pc_y[j] * tez_zzz_xxyy_1[j] + fl1_fx * tez_zzz_xxy_0[j] - fl1_fx * tez_zzz_xxy_1[j];

                tex_yzzz_xxyz_0[j] =
                    pa_y[j] * tex_zzz_xxyz_0[j] - pc_y[j] * tex_zzz_xxyz_1[j] + 0.5 * fl1_fx * tex_zzz_xxz_0[j] - 0.5 * fl1_fx * tex_zzz_xxz_1[j];

                tey_yzzz_xxyz_0[j] = pa_y[j] * tey_zzz_xxyz_0[j] - pc_y[j] * tey_zzz_xxyz_1[j] + 0.5 * fl1_fx * tey_zzz_xxz_0[j] -
                                     0.5 * fl1_fx * tey_zzz_xxz_1[j] + ta_zzz_xxyz_1[j];

                tez_yzzz_xxyz_0[j] =
                    pa_y[j] * tez_zzz_xxyz_0[j] - pc_y[j] * tez_zzz_xxyz_1[j] + 0.5 * fl1_fx * tez_zzz_xxz_0[j] - 0.5 * fl1_fx * tez_zzz_xxz_1[j];

                tex_yzzz_xxzz_0[j] = pa_y[j] * tex_zzz_xxzz_0[j] - pc_y[j] * tex_zzz_xxzz_1[j];

                tey_yzzz_xxzz_0[j] = pa_y[j] * tey_zzz_xxzz_0[j] - pc_y[j] * tey_zzz_xxzz_1[j] + ta_zzz_xxzz_1[j];

                tez_yzzz_xxzz_0[j] = pa_y[j] * tez_zzz_xxzz_0[j] - pc_y[j] * tez_zzz_xxzz_1[j];

                tex_yzzz_xyyy_0[j] =
                    pa_y[j] * tex_zzz_xyyy_0[j] - pc_y[j] * tex_zzz_xyyy_1[j] + 1.5 * fl1_fx * tex_zzz_xyy_0[j] - 1.5 * fl1_fx * tex_zzz_xyy_1[j];

                tey_yzzz_xyyy_0[j] = pa_y[j] * tey_zzz_xyyy_0[j] - pc_y[j] * tey_zzz_xyyy_1[j] + 1.5 * fl1_fx * tey_zzz_xyy_0[j] -
                                     1.5 * fl1_fx * tey_zzz_xyy_1[j] + ta_zzz_xyyy_1[j];

                tez_yzzz_xyyy_0[j] =
                    pa_y[j] * tez_zzz_xyyy_0[j] - pc_y[j] * tez_zzz_xyyy_1[j] + 1.5 * fl1_fx * tez_zzz_xyy_0[j] - 1.5 * fl1_fx * tez_zzz_xyy_1[j];

                tex_yzzz_xyyz_0[j] =
                    pa_y[j] * tex_zzz_xyyz_0[j] - pc_y[j] * tex_zzz_xyyz_1[j] + fl1_fx * tex_zzz_xyz_0[j] - fl1_fx * tex_zzz_xyz_1[j];

                tey_yzzz_xyyz_0[j] = pa_y[j] * tey_zzz_xyyz_0[j] - pc_y[j] * tey_zzz_xyyz_1[j] + fl1_fx * tey_zzz_xyz_0[j] -
                                     fl1_fx * tey_zzz_xyz_1[j] + ta_zzz_xyyz_1[j];

                tez_yzzz_xyyz_0[j] =
                    pa_y[j] * tez_zzz_xyyz_0[j] - pc_y[j] * tez_zzz_xyyz_1[j] + fl1_fx * tez_zzz_xyz_0[j] - fl1_fx * tez_zzz_xyz_1[j];

                tex_yzzz_xyzz_0[j] =
                    pa_y[j] * tex_zzz_xyzz_0[j] - pc_y[j] * tex_zzz_xyzz_1[j] + 0.5 * fl1_fx * tex_zzz_xzz_0[j] - 0.5 * fl1_fx * tex_zzz_xzz_1[j];

                tey_yzzz_xyzz_0[j] = pa_y[j] * tey_zzz_xyzz_0[j] - pc_y[j] * tey_zzz_xyzz_1[j] + 0.5 * fl1_fx * tey_zzz_xzz_0[j] -
                                     0.5 * fl1_fx * tey_zzz_xzz_1[j] + ta_zzz_xyzz_1[j];

                tez_yzzz_xyzz_0[j] =
                    pa_y[j] * tez_zzz_xyzz_0[j] - pc_y[j] * tez_zzz_xyzz_1[j] + 0.5 * fl1_fx * tez_zzz_xzz_0[j] - 0.5 * fl1_fx * tez_zzz_xzz_1[j];

                tex_yzzz_xzzz_0[j] = pa_y[j] * tex_zzz_xzzz_0[j] - pc_y[j] * tex_zzz_xzzz_1[j];

                tey_yzzz_xzzz_0[j] = pa_y[j] * tey_zzz_xzzz_0[j] - pc_y[j] * tey_zzz_xzzz_1[j] + ta_zzz_xzzz_1[j];

                tez_yzzz_xzzz_0[j] = pa_y[j] * tez_zzz_xzzz_0[j] - pc_y[j] * tez_zzz_xzzz_1[j];

                tex_yzzz_yyyy_0[j] =
                    pa_y[j] * tex_zzz_yyyy_0[j] - pc_y[j] * tex_zzz_yyyy_1[j] + 2.0 * fl1_fx * tex_zzz_yyy_0[j] - 2.0 * fl1_fx * tex_zzz_yyy_1[j];

                tey_yzzz_yyyy_0[j] = pa_y[j] * tey_zzz_yyyy_0[j] - pc_y[j] * tey_zzz_yyyy_1[j] + 2.0 * fl1_fx * tey_zzz_yyy_0[j] -
                                     2.0 * fl1_fx * tey_zzz_yyy_1[j] + ta_zzz_yyyy_1[j];

                tez_yzzz_yyyy_0[j] =
                    pa_y[j] * tez_zzz_yyyy_0[j] - pc_y[j] * tez_zzz_yyyy_1[j] + 2.0 * fl1_fx * tez_zzz_yyy_0[j] - 2.0 * fl1_fx * tez_zzz_yyy_1[j];

                tex_yzzz_yyyz_0[j] =
                    pa_y[j] * tex_zzz_yyyz_0[j] - pc_y[j] * tex_zzz_yyyz_1[j] + 1.5 * fl1_fx * tex_zzz_yyz_0[j] - 1.5 * fl1_fx * tex_zzz_yyz_1[j];

                tey_yzzz_yyyz_0[j] = pa_y[j] * tey_zzz_yyyz_0[j] - pc_y[j] * tey_zzz_yyyz_1[j] + 1.5 * fl1_fx * tey_zzz_yyz_0[j] -
                                     1.5 * fl1_fx * tey_zzz_yyz_1[j] + ta_zzz_yyyz_1[j];

                tez_yzzz_yyyz_0[j] =
                    pa_y[j] * tez_zzz_yyyz_0[j] - pc_y[j] * tez_zzz_yyyz_1[j] + 1.5 * fl1_fx * tez_zzz_yyz_0[j] - 1.5 * fl1_fx * tez_zzz_yyz_1[j];

                tex_yzzz_yyzz_0[j] =
                    pa_y[j] * tex_zzz_yyzz_0[j] - pc_y[j] * tex_zzz_yyzz_1[j] + fl1_fx * tex_zzz_yzz_0[j] - fl1_fx * tex_zzz_yzz_1[j];

                tey_yzzz_yyzz_0[j] = pa_y[j] * tey_zzz_yyzz_0[j] - pc_y[j] * tey_zzz_yyzz_1[j] + fl1_fx * tey_zzz_yzz_0[j] -
                                     fl1_fx * tey_zzz_yzz_1[j] + ta_zzz_yyzz_1[j];

                tez_yzzz_yyzz_0[j] =
                    pa_y[j] * tez_zzz_yyzz_0[j] - pc_y[j] * tez_zzz_yyzz_1[j] + fl1_fx * tez_zzz_yzz_0[j] - fl1_fx * tez_zzz_yzz_1[j];

                tex_yzzz_yzzz_0[j] =
                    pa_y[j] * tex_zzz_yzzz_0[j] - pc_y[j] * tex_zzz_yzzz_1[j] + 0.5 * fl1_fx * tex_zzz_zzz_0[j] - 0.5 * fl1_fx * tex_zzz_zzz_1[j];

                tey_yzzz_yzzz_0[j] = pa_y[j] * tey_zzz_yzzz_0[j] - pc_y[j] * tey_zzz_yzzz_1[j] + 0.5 * fl1_fx * tey_zzz_zzz_0[j] -
                                     0.5 * fl1_fx * tey_zzz_zzz_1[j] + ta_zzz_yzzz_1[j];

                tez_yzzz_yzzz_0[j] =
                    pa_y[j] * tez_zzz_yzzz_0[j] - pc_y[j] * tez_zzz_yzzz_1[j] + 0.5 * fl1_fx * tez_zzz_zzz_0[j] - 0.5 * fl1_fx * tez_zzz_zzz_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForGG_627_675(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pcDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto)
{
    // Batch of Integrals (627,675)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_4_4_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_3_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_2_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tex_zzz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 135);

            auto tey_zzz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 135);

            auto tez_zzz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 135);

            auto tex_zzz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 136);

            auto tey_zzz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 136);

            auto tez_zzz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 136);

            auto tex_zzz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 137);

            auto tey_zzz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 137);

            auto tez_zzz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 137);

            auto tex_zzz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 138);

            auto tey_zzz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 138);

            auto tez_zzz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 138);

            auto tex_zzz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 139);

            auto tey_zzz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 139);

            auto tez_zzz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 139);

            auto tex_zzz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 140);

            auto tey_zzz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 140);

            auto tez_zzz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 140);

            auto tex_zzz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 141);

            auto tey_zzz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 141);

            auto tez_zzz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 141);

            auto tex_zzz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 142);

            auto tey_zzz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 142);

            auto tez_zzz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 142);

            auto tex_zzz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 143);

            auto tey_zzz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 143);

            auto tez_zzz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 143);

            auto tex_zzz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 144);

            auto tey_zzz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 144);

            auto tez_zzz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 144);

            auto tex_zzz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 145);

            auto tey_zzz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 145);

            auto tez_zzz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 145);

            auto tex_zzz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 146);

            auto tey_zzz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 146);

            auto tez_zzz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 146);

            auto tex_zzz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 147);

            auto tey_zzz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 147);

            auto tez_zzz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 147);

            auto tex_zzz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 148);

            auto tey_zzz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 148);

            auto tez_zzz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 148);

            auto tex_zzz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 149);

            auto tey_zzz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 149);

            auto tez_zzz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 149);

            auto tex_zzz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 135);

            auto tey_zzz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 135);

            auto tez_zzz_xxxx_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 135);

            auto tex_zzz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 136);

            auto tey_zzz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 136);

            auto tez_zzz_xxxy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 136);

            auto tex_zzz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 137);

            auto tey_zzz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 137);

            auto tez_zzz_xxxz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 137);

            auto tex_zzz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 138);

            auto tey_zzz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 138);

            auto tez_zzz_xxyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 138);

            auto tex_zzz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 139);

            auto tey_zzz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 139);

            auto tez_zzz_xxyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 139);

            auto tex_zzz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 140);

            auto tey_zzz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 140);

            auto tez_zzz_xxzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 140);

            auto tex_zzz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 141);

            auto tey_zzz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 141);

            auto tez_zzz_xyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 141);

            auto tex_zzz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 142);

            auto tey_zzz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 142);

            auto tez_zzz_xyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 142);

            auto tex_zzz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 143);

            auto tey_zzz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 143);

            auto tez_zzz_xyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 143);

            auto tex_zzz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 144);

            auto tey_zzz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 144);

            auto tez_zzz_xzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 144);

            auto tex_zzz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 145);

            auto tey_zzz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 145);

            auto tez_zzz_yyyy_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 145);

            auto tex_zzz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 146);

            auto tey_zzz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 146);

            auto tez_zzz_yyyz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 146);

            auto tex_zzz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 147);

            auto tey_zzz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 147);

            auto tez_zzz_yyzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 147);

            auto tex_zzz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 148);

            auto tey_zzz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 148);

            auto tez_zzz_yzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 148);

            auto tex_zzz_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * idx + 149);

            auto tey_zzz_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 150 * bdim + 150 * idx + 149);

            auto tez_zzz_zzzz_1 = primBuffer.data(pidx_e_3_4_m1 + 300 * bdim + 150 * idx + 149);

            auto tex_zz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 75);

            auto tey_zz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 75);

            auto tez_zz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 75);

            auto tex_zz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 76);

            auto tey_zz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 76);

            auto tez_zz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 76);

            auto tex_zz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 77);

            auto tey_zz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 77);

            auto tez_zz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 77);

            auto tex_zz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 78);

            auto tey_zz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 78);

            auto tez_zz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 78);

            auto tex_zz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 79);

            auto tey_zz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 79);

            auto tez_zz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 79);

            auto tex_zz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 80);

            auto tey_zz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 80);

            auto tez_zz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 80);

            auto tex_zz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 81);

            auto tey_zz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 81);

            auto tez_zz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 81);

            auto tex_zz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 82);

            auto tey_zz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 82);

            auto tez_zz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 82);

            auto tex_zz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 83);

            auto tey_zz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 83);

            auto tez_zz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 83);

            auto tex_zz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 84);

            auto tey_zz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 84);

            auto tez_zz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 84);

            auto tex_zz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 85);

            auto tey_zz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 85);

            auto tez_zz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 85);

            auto tex_zz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 86);

            auto tey_zz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 86);

            auto tez_zz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 86);

            auto tex_zz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 87);

            auto tey_zz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 87);

            auto tez_zz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 87);

            auto tex_zz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 88);

            auto tey_zz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 88);

            auto tez_zz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 88);

            auto tex_zz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 89);

            auto tey_zz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 89);

            auto tez_zz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 89);

            auto tex_zz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 75);

            auto tey_zz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 75);

            auto tez_zz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 75);

            auto tex_zz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 76);

            auto tey_zz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 76);

            auto tez_zz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 76);

            auto tex_zz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 77);

            auto tey_zz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 77);

            auto tez_zz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 77);

            auto tex_zz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 78);

            auto tey_zz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 78);

            auto tez_zz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 78);

            auto tex_zz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 79);

            auto tey_zz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 79);

            auto tez_zz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 79);

            auto tex_zz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 80);

            auto tey_zz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 80);

            auto tez_zz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 80);

            auto tex_zz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 81);

            auto tey_zz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 81);

            auto tez_zz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 81);

            auto tex_zz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 82);

            auto tey_zz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 82);

            auto tez_zz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 82);

            auto tex_zz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 83);

            auto tey_zz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 83);

            auto tez_zz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 83);

            auto tex_zz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 84);

            auto tey_zz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 84);

            auto tez_zz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 84);

            auto tex_zz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 85);

            auto tey_zz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 85);

            auto tez_zz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 85);

            auto tex_zz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 86);

            auto tey_zz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 86);

            auto tez_zz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 86);

            auto tex_zz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 87);

            auto tey_zz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 87);

            auto tez_zz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 87);

            auto tex_zz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 88);

            auto tey_zz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 88);

            auto tez_zz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 88);

            auto tex_zz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 89);

            auto tey_zz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 89);

            auto tez_zz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 89);

            auto tex_zzz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 90);

            auto tey_zzz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 90);

            auto tez_zzz_xxx_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 90);

            auto tex_zzz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 91);

            auto tey_zzz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 91);

            auto tez_zzz_xxy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 91);

            auto tex_zzz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 92);

            auto tey_zzz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 92);

            auto tez_zzz_xxz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 92);

            auto tex_zzz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 93);

            auto tey_zzz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 93);

            auto tez_zzz_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 93);

            auto tex_zzz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 94);

            auto tey_zzz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 94);

            auto tez_zzz_xyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 94);

            auto tex_zzz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 95);

            auto tey_zzz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 95);

            auto tez_zzz_xzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 95);

            auto tex_zzz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 96);

            auto tey_zzz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 96);

            auto tez_zzz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 96);

            auto tex_zzz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 97);

            auto tey_zzz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 97);

            auto tez_zzz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 97);

            auto tex_zzz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 98);

            auto tey_zzz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 98);

            auto tez_zzz_yzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 98);

            auto tex_zzz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 99);

            auto tey_zzz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 99);

            auto tez_zzz_zzz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 99);

            auto tex_zzz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 90);

            auto tey_zzz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 90);

            auto tez_zzz_xxx_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 90);

            auto tex_zzz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 91);

            auto tey_zzz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 91);

            auto tez_zzz_xxy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 91);

            auto tex_zzz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 92);

            auto tey_zzz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 92);

            auto tez_zzz_xxz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 92);

            auto tex_zzz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 93);

            auto tey_zzz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 93);

            auto tez_zzz_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 93);

            auto tex_zzz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 94);

            auto tey_zzz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 94);

            auto tez_zzz_xyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 94);

            auto tex_zzz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 95);

            auto tey_zzz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 95);

            auto tez_zzz_xzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 95);

            auto tex_zzz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 96);

            auto tey_zzz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 96);

            auto tez_zzz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 96);

            auto tex_zzz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 97);

            auto tey_zzz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 97);

            auto tez_zzz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 97);

            auto tex_zzz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 98);

            auto tey_zzz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 98);

            auto tez_zzz_yzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 98);

            auto tex_zzz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 99);

            auto tey_zzz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 99);

            auto tez_zzz_zzz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 99);

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

            // set up pointers to integrals

            auto tex_yzzz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 209);

            auto tey_yzzz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 209);

            auto tez_yzzz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 209);

            auto tex_zzzz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 210);

            auto tey_zzzz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 210);

            auto tez_zzzz_xxxx_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 210);

            auto tex_zzzz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 211);

            auto tey_zzzz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 211);

            auto tez_zzzz_xxxy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 211);

            auto tex_zzzz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 212);

            auto tey_zzzz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 212);

            auto tez_zzzz_xxxz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 212);

            auto tex_zzzz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 213);

            auto tey_zzzz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 213);

            auto tez_zzzz_xxyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 213);

            auto tex_zzzz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 214);

            auto tey_zzzz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 214);

            auto tez_zzzz_xxyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 214);

            auto tex_zzzz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 215);

            auto tey_zzzz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 215);

            auto tez_zzzz_xxzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 215);

            auto tex_zzzz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 216);

            auto tey_zzzz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 216);

            auto tez_zzzz_xyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 216);

            auto tex_zzzz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 217);

            auto tey_zzzz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 217);

            auto tez_zzzz_xyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 217);

            auto tex_zzzz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 218);

            auto tey_zzzz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 218);

            auto tez_zzzz_xyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 218);

            auto tex_zzzz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 219);

            auto tey_zzzz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 219);

            auto tez_zzzz_xzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 219);

            auto tex_zzzz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 220);

            auto tey_zzzz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 220);

            auto tez_zzzz_yyyy_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 220);

            auto tex_zzzz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 221);

            auto tey_zzzz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 221);

            auto tez_zzzz_yyyz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 221);

            auto tex_zzzz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 222);

            auto tey_zzzz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 222);

            auto tez_zzzz_yyzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 222);

            auto tex_zzzz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 223);

            auto tey_zzzz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 223);

            auto tez_zzzz_yzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 223);

            auto tex_zzzz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * idx + 224);

            auto tey_zzzz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 225 * bdim + 225 * idx + 224);

            auto tez_zzzz_zzzz_0 = primBuffer.data(pidx_e_4_4_m0 + 450 * bdim + 225 * idx + 224);

            // Batch of Integrals (627,675)

            #pragma omp simd aligned(fx, pa_y, pa_z, pc_y, pc_z, ta_zzz_xxxx_1, ta_zzz_xxxy_1, ta_zzz_xxxz_1, \
                                         ta_zzz_xxyy_1, ta_zzz_xxyz_1, ta_zzz_xxzz_1, ta_zzz_xyyy_1, ta_zzz_xyyz_1, \
                                         ta_zzz_xyzz_1, ta_zzz_xzzz_1, ta_zzz_yyyy_1, ta_zzz_yyyz_1, ta_zzz_yyzz_1, \
                                         ta_zzz_yzzz_1, ta_zzz_zzzz_1, tex_yzzz_zzzz_0, tex_zz_xxxx_0, tex_zz_xxxx_1, \
                                         tex_zz_xxxy_0, tex_zz_xxxy_1, tex_zz_xxxz_0, tex_zz_xxxz_1, tex_zz_xxyy_0, \
                                         tex_zz_xxyy_1, tex_zz_xxyz_0, tex_zz_xxyz_1, tex_zz_xxzz_0, tex_zz_xxzz_1, \
                                         tex_zz_xyyy_0, tex_zz_xyyy_1, tex_zz_xyyz_0, tex_zz_xyyz_1, tex_zz_xyzz_0, \
                                         tex_zz_xyzz_1, tex_zz_xzzz_0, tex_zz_xzzz_1, tex_zz_yyyy_0, tex_zz_yyyy_1, \
                                         tex_zz_yyyz_0, tex_zz_yyyz_1, tex_zz_yyzz_0, tex_zz_yyzz_1, tex_zz_yzzz_0, \
                                         tex_zz_yzzz_1, tex_zz_zzzz_0, tex_zz_zzzz_1, tex_zzz_xxx_0, tex_zzz_xxx_1, \
                                         tex_zzz_xxxx_0, tex_zzz_xxxx_1, tex_zzz_xxxy_0, tex_zzz_xxxy_1, tex_zzz_xxxz_0, \
                                         tex_zzz_xxxz_1, tex_zzz_xxy_0, tex_zzz_xxy_1, tex_zzz_xxyy_0, tex_zzz_xxyy_1, \
                                         tex_zzz_xxyz_0, tex_zzz_xxyz_1, tex_zzz_xxz_0, tex_zzz_xxz_1, tex_zzz_xxzz_0, \
                                         tex_zzz_xxzz_1, tex_zzz_xyy_0, tex_zzz_xyy_1, tex_zzz_xyyy_0, tex_zzz_xyyy_1, \
                                         tex_zzz_xyyz_0, tex_zzz_xyyz_1, tex_zzz_xyz_0, tex_zzz_xyz_1, tex_zzz_xyzz_0, \
                                         tex_zzz_xyzz_1, tex_zzz_xzz_0, tex_zzz_xzz_1, tex_zzz_xzzz_0, tex_zzz_xzzz_1, \
                                         tex_zzz_yyy_0, tex_zzz_yyy_1, tex_zzz_yyyy_0, tex_zzz_yyyy_1, tex_zzz_yyyz_0, \
                                         tex_zzz_yyyz_1, tex_zzz_yyz_0, tex_zzz_yyz_1, tex_zzz_yyzz_0, tex_zzz_yyzz_1, \
                                         tex_zzz_yzz_0, tex_zzz_yzz_1, tex_zzz_yzzz_0, tex_zzz_yzzz_1, tex_zzz_zzz_0, \
                                         tex_zzz_zzz_1, tex_zzz_zzzz_0, tex_zzz_zzzz_1, tex_zzzz_xxxx_0, tex_zzzz_xxxy_0, \
                                         tex_zzzz_xxxz_0, tex_zzzz_xxyy_0, tex_zzzz_xxyz_0, tex_zzzz_xxzz_0, tex_zzzz_xyyy_0, \
                                         tex_zzzz_xyyz_0, tex_zzzz_xyzz_0, tex_zzzz_xzzz_0, tex_zzzz_yyyy_0, tex_zzzz_yyyz_0, \
                                         tex_zzzz_yyzz_0, tex_zzzz_yzzz_0, tex_zzzz_zzzz_0, tey_yzzz_zzzz_0, tey_zz_xxxx_0, \
                                         tey_zz_xxxx_1, tey_zz_xxxy_0, tey_zz_xxxy_1, tey_zz_xxxz_0, tey_zz_xxxz_1, \
                                         tey_zz_xxyy_0, tey_zz_xxyy_1, tey_zz_xxyz_0, tey_zz_xxyz_1, tey_zz_xxzz_0, \
                                         tey_zz_xxzz_1, tey_zz_xyyy_0, tey_zz_xyyy_1, tey_zz_xyyz_0, tey_zz_xyyz_1, \
                                         tey_zz_xyzz_0, tey_zz_xyzz_1, tey_zz_xzzz_0, tey_zz_xzzz_1, tey_zz_yyyy_0, \
                                         tey_zz_yyyy_1, tey_zz_yyyz_0, tey_zz_yyyz_1, tey_zz_yyzz_0, tey_zz_yyzz_1, \
                                         tey_zz_yzzz_0, tey_zz_yzzz_1, tey_zz_zzzz_0, tey_zz_zzzz_1, tey_zzz_xxx_0, \
                                         tey_zzz_xxx_1, tey_zzz_xxxx_0, tey_zzz_xxxx_1, tey_zzz_xxxy_0, tey_zzz_xxxy_1, \
                                         tey_zzz_xxxz_0, tey_zzz_xxxz_1, tey_zzz_xxy_0, tey_zzz_xxy_1, tey_zzz_xxyy_0, \
                                         tey_zzz_xxyy_1, tey_zzz_xxyz_0, tey_zzz_xxyz_1, tey_zzz_xxz_0, tey_zzz_xxz_1, \
                                         tey_zzz_xxzz_0, tey_zzz_xxzz_1, tey_zzz_xyy_0, tey_zzz_xyy_1, tey_zzz_xyyy_0, \
                                         tey_zzz_xyyy_1, tey_zzz_xyyz_0, tey_zzz_xyyz_1, tey_zzz_xyz_0, tey_zzz_xyz_1, \
                                         tey_zzz_xyzz_0, tey_zzz_xyzz_1, tey_zzz_xzz_0, tey_zzz_xzz_1, tey_zzz_xzzz_0, \
                                         tey_zzz_xzzz_1, tey_zzz_yyy_0, tey_zzz_yyy_1, tey_zzz_yyyy_0, tey_zzz_yyyy_1, \
                                         tey_zzz_yyyz_0, tey_zzz_yyyz_1, tey_zzz_yyz_0, tey_zzz_yyz_1, tey_zzz_yyzz_0, \
                                         tey_zzz_yyzz_1, tey_zzz_yzz_0, tey_zzz_yzz_1, tey_zzz_yzzz_0, tey_zzz_yzzz_1, \
                                         tey_zzz_zzz_0, tey_zzz_zzz_1, tey_zzz_zzzz_0, tey_zzz_zzzz_1, tey_zzzz_xxxx_0, \
                                         tey_zzzz_xxxy_0, tey_zzzz_xxxz_0, tey_zzzz_xxyy_0, tey_zzzz_xxyz_0, tey_zzzz_xxzz_0, \
                                         tey_zzzz_xyyy_0, tey_zzzz_xyyz_0, tey_zzzz_xyzz_0, tey_zzzz_xzzz_0, tey_zzzz_yyyy_0, \
                                         tey_zzzz_yyyz_0, tey_zzzz_yyzz_0, tey_zzzz_yzzz_0, tey_zzzz_zzzz_0, tez_yzzz_zzzz_0, \
                                         tez_zz_xxxx_0, tez_zz_xxxx_1, tez_zz_xxxy_0, tez_zz_xxxy_1, tez_zz_xxxz_0, \
                                         tez_zz_xxxz_1, tez_zz_xxyy_0, tez_zz_xxyy_1, tez_zz_xxyz_0, tez_zz_xxyz_1, \
                                         tez_zz_xxzz_0, tez_zz_xxzz_1, tez_zz_xyyy_0, tez_zz_xyyy_1, tez_zz_xyyz_0, \
                                         tez_zz_xyyz_1, tez_zz_xyzz_0, tez_zz_xyzz_1, tez_zz_xzzz_0, tez_zz_xzzz_1, \
                                         tez_zz_yyyy_0, tez_zz_yyyy_1, tez_zz_yyyz_0, tez_zz_yyyz_1, tez_zz_yyzz_0, \
                                         tez_zz_yyzz_1, tez_zz_yzzz_0, tez_zz_yzzz_1, tez_zz_zzzz_0, tez_zz_zzzz_1, \
                                         tez_zzz_xxx_0, tez_zzz_xxx_1, tez_zzz_xxxx_0, tez_zzz_xxxx_1, tez_zzz_xxxy_0, \
                                         tez_zzz_xxxy_1, tez_zzz_xxxz_0, tez_zzz_xxxz_1, tez_zzz_xxy_0, tez_zzz_xxy_1, \
                                         tez_zzz_xxyy_0, tez_zzz_xxyy_1, tez_zzz_xxyz_0, tez_zzz_xxyz_1, tez_zzz_xxz_0, \
                                         tez_zzz_xxz_1, tez_zzz_xxzz_0, tez_zzz_xxzz_1, tez_zzz_xyy_0, tez_zzz_xyy_1, \
                                         tez_zzz_xyyy_0, tez_zzz_xyyy_1, tez_zzz_xyyz_0, tez_zzz_xyyz_1, tez_zzz_xyz_0, \
                                         tez_zzz_xyz_1, tez_zzz_xyzz_0, tez_zzz_xyzz_1, tez_zzz_xzz_0, tez_zzz_xzz_1, \
                                         tez_zzz_xzzz_0, tez_zzz_xzzz_1, tez_zzz_yyy_0, tez_zzz_yyy_1, tez_zzz_yyyy_0, \
                                         tez_zzz_yyyy_1, tez_zzz_yyyz_0, tez_zzz_yyyz_1, tez_zzz_yyz_0, tez_zzz_yyz_1, \
                                         tez_zzz_yyzz_0, tez_zzz_yyzz_1, tez_zzz_yzz_0, tez_zzz_yzz_1, tez_zzz_yzzz_0, \
                                         tez_zzz_yzzz_1, tez_zzz_zzz_0, tez_zzz_zzz_1, tez_zzz_zzzz_0, tez_zzz_zzzz_1, \
                                         tez_zzzz_xxxx_0, tez_zzzz_xxxy_0, tez_zzzz_xxxz_0, tez_zzzz_xxyy_0, tez_zzzz_xxyz_0, \
                                         tez_zzzz_xxzz_0, tez_zzzz_xyyy_0, tez_zzzz_xyyz_0, tez_zzzz_xyzz_0, tez_zzzz_xzzz_0, \
                                         tez_zzzz_yyyy_0, tez_zzzz_yyyz_0, tez_zzzz_yyzz_0, tez_zzzz_yzzz_0, tez_zzzz_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_yzzz_zzzz_0[j] = pa_y[j] * tex_zzz_zzzz_0[j] - pc_y[j] * tex_zzz_zzzz_1[j];

                tey_yzzz_zzzz_0[j] = pa_y[j] * tey_zzz_zzzz_0[j] - pc_y[j] * tey_zzz_zzzz_1[j] + ta_zzz_zzzz_1[j];

                tez_yzzz_zzzz_0[j] = pa_y[j] * tez_zzz_zzzz_0[j] - pc_y[j] * tez_zzz_zzzz_1[j];

                tex_zzzz_xxxx_0[j] =
                    pa_z[j] * tex_zzz_xxxx_0[j] - pc_z[j] * tex_zzz_xxxx_1[j] + 1.5 * fl1_fx * tex_zz_xxxx_0[j] - 1.5 * fl1_fx * tex_zz_xxxx_1[j];

                tey_zzzz_xxxx_0[j] =
                    pa_z[j] * tey_zzz_xxxx_0[j] - pc_z[j] * tey_zzz_xxxx_1[j] + 1.5 * fl1_fx * tey_zz_xxxx_0[j] - 1.5 * fl1_fx * tey_zz_xxxx_1[j];

                tez_zzzz_xxxx_0[j] = pa_z[j] * tez_zzz_xxxx_0[j] - pc_z[j] * tez_zzz_xxxx_1[j] + 1.5 * fl1_fx * tez_zz_xxxx_0[j] -
                                     1.5 * fl1_fx * tez_zz_xxxx_1[j] + ta_zzz_xxxx_1[j];

                tex_zzzz_xxxy_0[j] =
                    pa_z[j] * tex_zzz_xxxy_0[j] - pc_z[j] * tex_zzz_xxxy_1[j] + 1.5 * fl1_fx * tex_zz_xxxy_0[j] - 1.5 * fl1_fx * tex_zz_xxxy_1[j];

                tey_zzzz_xxxy_0[j] =
                    pa_z[j] * tey_zzz_xxxy_0[j] - pc_z[j] * tey_zzz_xxxy_1[j] + 1.5 * fl1_fx * tey_zz_xxxy_0[j] - 1.5 * fl1_fx * tey_zz_xxxy_1[j];

                tez_zzzz_xxxy_0[j] = pa_z[j] * tez_zzz_xxxy_0[j] - pc_z[j] * tez_zzz_xxxy_1[j] + 1.5 * fl1_fx * tez_zz_xxxy_0[j] -
                                     1.5 * fl1_fx * tez_zz_xxxy_1[j] + ta_zzz_xxxy_1[j];

                tex_zzzz_xxxz_0[j] = pa_z[j] * tex_zzz_xxxz_0[j] - pc_z[j] * tex_zzz_xxxz_1[j] + 1.5 * fl1_fx * tex_zz_xxxz_0[j] -
                                     1.5 * fl1_fx * tex_zz_xxxz_1[j] + 0.5 * fl1_fx * tex_zzz_xxx_0[j] - 0.5 * fl1_fx * tex_zzz_xxx_1[j];

                tey_zzzz_xxxz_0[j] = pa_z[j] * tey_zzz_xxxz_0[j] - pc_z[j] * tey_zzz_xxxz_1[j] + 1.5 * fl1_fx * tey_zz_xxxz_0[j] -
                                     1.5 * fl1_fx * tey_zz_xxxz_1[j] + 0.5 * fl1_fx * tey_zzz_xxx_0[j] - 0.5 * fl1_fx * tey_zzz_xxx_1[j];

                tez_zzzz_xxxz_0[j] = pa_z[j] * tez_zzz_xxxz_0[j] - pc_z[j] * tez_zzz_xxxz_1[j] + 1.5 * fl1_fx * tez_zz_xxxz_0[j] -
                                     1.5 * fl1_fx * tez_zz_xxxz_1[j] + 0.5 * fl1_fx * tez_zzz_xxx_0[j] - 0.5 * fl1_fx * tez_zzz_xxx_1[j] +
                                     ta_zzz_xxxz_1[j];

                tex_zzzz_xxyy_0[j] =
                    pa_z[j] * tex_zzz_xxyy_0[j] - pc_z[j] * tex_zzz_xxyy_1[j] + 1.5 * fl1_fx * tex_zz_xxyy_0[j] - 1.5 * fl1_fx * tex_zz_xxyy_1[j];

                tey_zzzz_xxyy_0[j] =
                    pa_z[j] * tey_zzz_xxyy_0[j] - pc_z[j] * tey_zzz_xxyy_1[j] + 1.5 * fl1_fx * tey_zz_xxyy_0[j] - 1.5 * fl1_fx * tey_zz_xxyy_1[j];

                tez_zzzz_xxyy_0[j] = pa_z[j] * tez_zzz_xxyy_0[j] - pc_z[j] * tez_zzz_xxyy_1[j] + 1.5 * fl1_fx * tez_zz_xxyy_0[j] -
                                     1.5 * fl1_fx * tez_zz_xxyy_1[j] + ta_zzz_xxyy_1[j];

                tex_zzzz_xxyz_0[j] = pa_z[j] * tex_zzz_xxyz_0[j] - pc_z[j] * tex_zzz_xxyz_1[j] + 1.5 * fl1_fx * tex_zz_xxyz_0[j] -
                                     1.5 * fl1_fx * tex_zz_xxyz_1[j] + 0.5 * fl1_fx * tex_zzz_xxy_0[j] - 0.5 * fl1_fx * tex_zzz_xxy_1[j];

                tey_zzzz_xxyz_0[j] = pa_z[j] * tey_zzz_xxyz_0[j] - pc_z[j] * tey_zzz_xxyz_1[j] + 1.5 * fl1_fx * tey_zz_xxyz_0[j] -
                                     1.5 * fl1_fx * tey_zz_xxyz_1[j] + 0.5 * fl1_fx * tey_zzz_xxy_0[j] - 0.5 * fl1_fx * tey_zzz_xxy_1[j];

                tez_zzzz_xxyz_0[j] = pa_z[j] * tez_zzz_xxyz_0[j] - pc_z[j] * tez_zzz_xxyz_1[j] + 1.5 * fl1_fx * tez_zz_xxyz_0[j] -
                                     1.5 * fl1_fx * tez_zz_xxyz_1[j] + 0.5 * fl1_fx * tez_zzz_xxy_0[j] - 0.5 * fl1_fx * tez_zzz_xxy_1[j] +
                                     ta_zzz_xxyz_1[j];

                tex_zzzz_xxzz_0[j] = pa_z[j] * tex_zzz_xxzz_0[j] - pc_z[j] * tex_zzz_xxzz_1[j] + 1.5 * fl1_fx * tex_zz_xxzz_0[j] -
                                     1.5 * fl1_fx * tex_zz_xxzz_1[j] + fl1_fx * tex_zzz_xxz_0[j] - fl1_fx * tex_zzz_xxz_1[j];

                tey_zzzz_xxzz_0[j] = pa_z[j] * tey_zzz_xxzz_0[j] - pc_z[j] * tey_zzz_xxzz_1[j] + 1.5 * fl1_fx * tey_zz_xxzz_0[j] -
                                     1.5 * fl1_fx * tey_zz_xxzz_1[j] + fl1_fx * tey_zzz_xxz_0[j] - fl1_fx * tey_zzz_xxz_1[j];

                tez_zzzz_xxzz_0[j] = pa_z[j] * tez_zzz_xxzz_0[j] - pc_z[j] * tez_zzz_xxzz_1[j] + 1.5 * fl1_fx * tez_zz_xxzz_0[j] -
                                     1.5 * fl1_fx * tez_zz_xxzz_1[j] + fl1_fx * tez_zzz_xxz_0[j] - fl1_fx * tez_zzz_xxz_1[j] + ta_zzz_xxzz_1[j];

                tex_zzzz_xyyy_0[j] =
                    pa_z[j] * tex_zzz_xyyy_0[j] - pc_z[j] * tex_zzz_xyyy_1[j] + 1.5 * fl1_fx * tex_zz_xyyy_0[j] - 1.5 * fl1_fx * tex_zz_xyyy_1[j];

                tey_zzzz_xyyy_0[j] =
                    pa_z[j] * tey_zzz_xyyy_0[j] - pc_z[j] * tey_zzz_xyyy_1[j] + 1.5 * fl1_fx * tey_zz_xyyy_0[j] - 1.5 * fl1_fx * tey_zz_xyyy_1[j];

                tez_zzzz_xyyy_0[j] = pa_z[j] * tez_zzz_xyyy_0[j] - pc_z[j] * tez_zzz_xyyy_1[j] + 1.5 * fl1_fx * tez_zz_xyyy_0[j] -
                                     1.5 * fl1_fx * tez_zz_xyyy_1[j] + ta_zzz_xyyy_1[j];

                tex_zzzz_xyyz_0[j] = pa_z[j] * tex_zzz_xyyz_0[j] - pc_z[j] * tex_zzz_xyyz_1[j] + 1.5 * fl1_fx * tex_zz_xyyz_0[j] -
                                     1.5 * fl1_fx * tex_zz_xyyz_1[j] + 0.5 * fl1_fx * tex_zzz_xyy_0[j] - 0.5 * fl1_fx * tex_zzz_xyy_1[j];

                tey_zzzz_xyyz_0[j] = pa_z[j] * tey_zzz_xyyz_0[j] - pc_z[j] * tey_zzz_xyyz_1[j] + 1.5 * fl1_fx * tey_zz_xyyz_0[j] -
                                     1.5 * fl1_fx * tey_zz_xyyz_1[j] + 0.5 * fl1_fx * tey_zzz_xyy_0[j] - 0.5 * fl1_fx * tey_zzz_xyy_1[j];

                tez_zzzz_xyyz_0[j] = pa_z[j] * tez_zzz_xyyz_0[j] - pc_z[j] * tez_zzz_xyyz_1[j] + 1.5 * fl1_fx * tez_zz_xyyz_0[j] -
                                     1.5 * fl1_fx * tez_zz_xyyz_1[j] + 0.5 * fl1_fx * tez_zzz_xyy_0[j] - 0.5 * fl1_fx * tez_zzz_xyy_1[j] +
                                     ta_zzz_xyyz_1[j];

                tex_zzzz_xyzz_0[j] = pa_z[j] * tex_zzz_xyzz_0[j] - pc_z[j] * tex_zzz_xyzz_1[j] + 1.5 * fl1_fx * tex_zz_xyzz_0[j] -
                                     1.5 * fl1_fx * tex_zz_xyzz_1[j] + fl1_fx * tex_zzz_xyz_0[j] - fl1_fx * tex_zzz_xyz_1[j];

                tey_zzzz_xyzz_0[j] = pa_z[j] * tey_zzz_xyzz_0[j] - pc_z[j] * tey_zzz_xyzz_1[j] + 1.5 * fl1_fx * tey_zz_xyzz_0[j] -
                                     1.5 * fl1_fx * tey_zz_xyzz_1[j] + fl1_fx * tey_zzz_xyz_0[j] - fl1_fx * tey_zzz_xyz_1[j];

                tez_zzzz_xyzz_0[j] = pa_z[j] * tez_zzz_xyzz_0[j] - pc_z[j] * tez_zzz_xyzz_1[j] + 1.5 * fl1_fx * tez_zz_xyzz_0[j] -
                                     1.5 * fl1_fx * tez_zz_xyzz_1[j] + fl1_fx * tez_zzz_xyz_0[j] - fl1_fx * tez_zzz_xyz_1[j] + ta_zzz_xyzz_1[j];

                tex_zzzz_xzzz_0[j] = pa_z[j] * tex_zzz_xzzz_0[j] - pc_z[j] * tex_zzz_xzzz_1[j] + 1.5 * fl1_fx * tex_zz_xzzz_0[j] -
                                     1.5 * fl1_fx * tex_zz_xzzz_1[j] + 1.5 * fl1_fx * tex_zzz_xzz_0[j] - 1.5 * fl1_fx * tex_zzz_xzz_1[j];

                tey_zzzz_xzzz_0[j] = pa_z[j] * tey_zzz_xzzz_0[j] - pc_z[j] * tey_zzz_xzzz_1[j] + 1.5 * fl1_fx * tey_zz_xzzz_0[j] -
                                     1.5 * fl1_fx * tey_zz_xzzz_1[j] + 1.5 * fl1_fx * tey_zzz_xzz_0[j] - 1.5 * fl1_fx * tey_zzz_xzz_1[j];

                tez_zzzz_xzzz_0[j] = pa_z[j] * tez_zzz_xzzz_0[j] - pc_z[j] * tez_zzz_xzzz_1[j] + 1.5 * fl1_fx * tez_zz_xzzz_0[j] -
                                     1.5 * fl1_fx * tez_zz_xzzz_1[j] + 1.5 * fl1_fx * tez_zzz_xzz_0[j] - 1.5 * fl1_fx * tez_zzz_xzz_1[j] +
                                     ta_zzz_xzzz_1[j];

                tex_zzzz_yyyy_0[j] =
                    pa_z[j] * tex_zzz_yyyy_0[j] - pc_z[j] * tex_zzz_yyyy_1[j] + 1.5 * fl1_fx * tex_zz_yyyy_0[j] - 1.5 * fl1_fx * tex_zz_yyyy_1[j];

                tey_zzzz_yyyy_0[j] =
                    pa_z[j] * tey_zzz_yyyy_0[j] - pc_z[j] * tey_zzz_yyyy_1[j] + 1.5 * fl1_fx * tey_zz_yyyy_0[j] - 1.5 * fl1_fx * tey_zz_yyyy_1[j];

                tez_zzzz_yyyy_0[j] = pa_z[j] * tez_zzz_yyyy_0[j] - pc_z[j] * tez_zzz_yyyy_1[j] + 1.5 * fl1_fx * tez_zz_yyyy_0[j] -
                                     1.5 * fl1_fx * tez_zz_yyyy_1[j] + ta_zzz_yyyy_1[j];

                tex_zzzz_yyyz_0[j] = pa_z[j] * tex_zzz_yyyz_0[j] - pc_z[j] * tex_zzz_yyyz_1[j] + 1.5 * fl1_fx * tex_zz_yyyz_0[j] -
                                     1.5 * fl1_fx * tex_zz_yyyz_1[j] + 0.5 * fl1_fx * tex_zzz_yyy_0[j] - 0.5 * fl1_fx * tex_zzz_yyy_1[j];

                tey_zzzz_yyyz_0[j] = pa_z[j] * tey_zzz_yyyz_0[j] - pc_z[j] * tey_zzz_yyyz_1[j] + 1.5 * fl1_fx * tey_zz_yyyz_0[j] -
                                     1.5 * fl1_fx * tey_zz_yyyz_1[j] + 0.5 * fl1_fx * tey_zzz_yyy_0[j] - 0.5 * fl1_fx * tey_zzz_yyy_1[j];

                tez_zzzz_yyyz_0[j] = pa_z[j] * tez_zzz_yyyz_0[j] - pc_z[j] * tez_zzz_yyyz_1[j] + 1.5 * fl1_fx * tez_zz_yyyz_0[j] -
                                     1.5 * fl1_fx * tez_zz_yyyz_1[j] + 0.5 * fl1_fx * tez_zzz_yyy_0[j] - 0.5 * fl1_fx * tez_zzz_yyy_1[j] +
                                     ta_zzz_yyyz_1[j];

                tex_zzzz_yyzz_0[j] = pa_z[j] * tex_zzz_yyzz_0[j] - pc_z[j] * tex_zzz_yyzz_1[j] + 1.5 * fl1_fx * tex_zz_yyzz_0[j] -
                                     1.5 * fl1_fx * tex_zz_yyzz_1[j] + fl1_fx * tex_zzz_yyz_0[j] - fl1_fx * tex_zzz_yyz_1[j];

                tey_zzzz_yyzz_0[j] = pa_z[j] * tey_zzz_yyzz_0[j] - pc_z[j] * tey_zzz_yyzz_1[j] + 1.5 * fl1_fx * tey_zz_yyzz_0[j] -
                                     1.5 * fl1_fx * tey_zz_yyzz_1[j] + fl1_fx * tey_zzz_yyz_0[j] - fl1_fx * tey_zzz_yyz_1[j];

                tez_zzzz_yyzz_0[j] = pa_z[j] * tez_zzz_yyzz_0[j] - pc_z[j] * tez_zzz_yyzz_1[j] + 1.5 * fl1_fx * tez_zz_yyzz_0[j] -
                                     1.5 * fl1_fx * tez_zz_yyzz_1[j] + fl1_fx * tez_zzz_yyz_0[j] - fl1_fx * tez_zzz_yyz_1[j] + ta_zzz_yyzz_1[j];

                tex_zzzz_yzzz_0[j] = pa_z[j] * tex_zzz_yzzz_0[j] - pc_z[j] * tex_zzz_yzzz_1[j] + 1.5 * fl1_fx * tex_zz_yzzz_0[j] -
                                     1.5 * fl1_fx * tex_zz_yzzz_1[j] + 1.5 * fl1_fx * tex_zzz_yzz_0[j] - 1.5 * fl1_fx * tex_zzz_yzz_1[j];

                tey_zzzz_yzzz_0[j] = pa_z[j] * tey_zzz_yzzz_0[j] - pc_z[j] * tey_zzz_yzzz_1[j] + 1.5 * fl1_fx * tey_zz_yzzz_0[j] -
                                     1.5 * fl1_fx * tey_zz_yzzz_1[j] + 1.5 * fl1_fx * tey_zzz_yzz_0[j] - 1.5 * fl1_fx * tey_zzz_yzz_1[j];

                tez_zzzz_yzzz_0[j] = pa_z[j] * tez_zzz_yzzz_0[j] - pc_z[j] * tez_zzz_yzzz_1[j] + 1.5 * fl1_fx * tez_zz_yzzz_0[j] -
                                     1.5 * fl1_fx * tez_zz_yzzz_1[j] + 1.5 * fl1_fx * tez_zzz_yzz_0[j] - 1.5 * fl1_fx * tez_zzz_yzz_1[j] +
                                     ta_zzz_yzzz_1[j];

                tex_zzzz_zzzz_0[j] = pa_z[j] * tex_zzz_zzzz_0[j] - pc_z[j] * tex_zzz_zzzz_1[j] + 1.5 * fl1_fx * tex_zz_zzzz_0[j] -
                                     1.5 * fl1_fx * tex_zz_zzzz_1[j] + 2.0 * fl1_fx * tex_zzz_zzz_0[j] - 2.0 * fl1_fx * tex_zzz_zzz_1[j];

                tey_zzzz_zzzz_0[j] = pa_z[j] * tey_zzz_zzzz_0[j] - pc_z[j] * tey_zzz_zzzz_1[j] + 1.5 * fl1_fx * tey_zz_zzzz_0[j] -
                                     1.5 * fl1_fx * tey_zz_zzzz_1[j] + 2.0 * fl1_fx * tey_zzz_zzz_0[j] - 2.0 * fl1_fx * tey_zzz_zzz_1[j];

                tez_zzzz_zzzz_0[j] = pa_z[j] * tez_zzz_zzzz_0[j] - pc_z[j] * tez_zzz_zzzz_1[j] + 1.5 * fl1_fx * tez_zz_zzzz_0[j] -
                                     1.5 * fl1_fx * tez_zz_zzzz_1[j] + 2.0 * fl1_fx * tez_zzz_zzz_0[j] - 2.0 * fl1_fx * tez_zzz_zzz_1[j] +
                                     ta_zzz_zzzz_1[j];
            }

            idx++;
        }
    }
}

}  // namespace efieldrecfunc
