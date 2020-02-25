//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ElectricFieldRecFuncForGF.hpp"

namespace efieldrecfunc {  // efieldrecfunc namespace

void
compElectricFieldForGF(CMemBlock2D<double>&       primBuffer,
                       const CRecursionMap&       recursionMap,
                       const CMemBlock2D<double>& osFactors,
                       const CMemBlock2D<double>& paDistances,
                       const CMemBlock2D<double>& pcDistances,
                       const CGtoBlock&           braGtoBlock,
                       const CGtoBlock&           ketGtoBlock,
                       const int32_t              iContrGto)
{
    efieldrecfunc::compElectricFieldForGF_0_50(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForGF_50_100(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForGF_100_150(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForGF_150_200(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForGF_200_250(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForGF_250_300(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForGF_300_350(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForGF_350_400(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForGF_400_450(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compElectricFieldForGF_0_50(CMemBlock2D<double>&       primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CMemBlock2D<double>& pcDistances,
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

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_4_3_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_2_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tex_xx_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx);

            auto tey_xx_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx);

            auto tez_xx_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx);

            auto tex_xx_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 1);

            auto tey_xx_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 1);

            auto tez_xx_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 1);

            auto tex_xx_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 2);

            auto tey_xx_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 2);

            auto tez_xx_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 2);

            auto tex_xx_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 3);

            auto tey_xx_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 3);

            auto tez_xx_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 3);

            auto tex_xx_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 4);

            auto tey_xx_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 4);

            auto tez_xx_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 4);

            auto tex_xx_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 5);

            auto tey_xx_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 5);

            auto tez_xx_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 5);

            auto tex_xx_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 6);

            auto tey_xx_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 6);

            auto tez_xx_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 6);

            auto tex_xx_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 7);

            auto tey_xx_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 7);

            auto tez_xx_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 7);

            auto tex_xx_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 8);

            auto tey_xx_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 8);

            auto tez_xx_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 8);

            auto tex_xx_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 9);

            auto tey_xx_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 9);

            auto tez_xx_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 9);

            auto tex_xy_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 10);

            auto tey_xy_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 10);

            auto tez_xy_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 10);

            auto tex_xy_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 11);

            auto tey_xy_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 11);

            auto tez_xy_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 11);

            auto tex_xy_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 12);

            auto tey_xy_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 12);

            auto tez_xy_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 12);

            auto tex_xy_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 13);

            auto tey_xy_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 13);

            auto tez_xy_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 13);

            auto tex_xy_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 14);

            auto tey_xy_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 14);

            auto tez_xy_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 14);

            auto tex_xy_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 15);

            auto tey_xy_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 15);

            auto tez_xy_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 15);

            auto tex_xy_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 16);

            auto tey_xy_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 16);

            auto tex_xx_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx);

            auto tey_xx_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx);

            auto tez_xx_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx);

            auto tex_xx_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 1);

            auto tey_xx_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 1);

            auto tez_xx_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 1);

            auto tex_xx_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 2);

            auto tey_xx_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 2);

            auto tez_xx_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 2);

            auto tex_xx_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 3);

            auto tey_xx_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 3);

            auto tez_xx_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 3);

            auto tex_xx_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 4);

            auto tey_xx_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 4);

            auto tez_xx_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 4);

            auto tex_xx_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 5);

            auto tey_xx_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 5);

            auto tez_xx_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 5);

            auto tex_xx_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 6);

            auto tey_xx_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 6);

            auto tez_xx_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 6);

            auto tex_xx_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 7);

            auto tey_xx_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 7);

            auto tez_xx_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 7);

            auto tex_xx_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 8);

            auto tey_xx_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 8);

            auto tez_xx_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 8);

            auto tex_xx_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 9);

            auto tey_xx_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 9);

            auto tez_xx_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 9);

            auto tex_xy_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 10);

            auto tey_xy_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 10);

            auto tez_xy_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 10);

            auto tex_xy_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 11);

            auto tey_xy_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 11);

            auto tez_xy_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 11);

            auto tex_xy_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 12);

            auto tey_xy_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 12);

            auto tez_xy_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 12);

            auto tex_xy_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 13);

            auto tey_xy_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 13);

            auto tez_xy_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 13);

            auto tex_xy_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 14);

            auto tey_xy_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 14);

            auto tez_xy_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 14);

            auto tex_xy_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 15);

            auto tey_xy_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 15);

            auto tez_xy_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 15);

            auto tex_xy_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 16);

            auto tey_xy_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 16);

            auto tex_xxx_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx);

            auto tey_xxx_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx);

            auto tez_xxx_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx);

            auto tex_xxx_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 1);

            auto tey_xxx_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 1);

            auto tez_xxx_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 1);

            auto tex_xxx_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 2);

            auto tey_xxx_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 2);

            auto tez_xxx_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 2);

            auto tex_xxx_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 3);

            auto tey_xxx_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 3);

            auto tez_xxx_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 3);

            auto tex_xxx_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 4);

            auto tey_xxx_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 4);

            auto tez_xxx_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 4);

            auto tex_xxx_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 5);

            auto tey_xxx_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 5);

            auto tez_xxx_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 5);

            auto tex_xxy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 6);

            auto tey_xxy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 6);

            auto tez_xxy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 6);

            auto tex_xxy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 7);

            auto tey_xxy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 7);

            auto tez_xxy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 7);

            auto tex_xxy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 8);

            auto tey_xxy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 8);

            auto tez_xxy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 8);

            auto tex_xxy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 9);

            auto tey_xxy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 9);

            auto tez_xxy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 9);

            auto tex_xxy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 10);

            auto tey_xxy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 10);

            auto tez_xxy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 10);

            auto tex_xxy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 11);

            auto tey_xxy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 11);

            auto tez_xxy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 11);

            auto tex_xxx_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx);

            auto tey_xxx_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx);

            auto tez_xxx_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx);

            auto tex_xxx_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 1);

            auto tey_xxx_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 1);

            auto tez_xxx_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 1);

            auto tex_xxx_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 2);

            auto tey_xxx_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 2);

            auto tez_xxx_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 2);

            auto tex_xxx_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 3);

            auto tey_xxx_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 3);

            auto tez_xxx_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 3);

            auto tex_xxx_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 4);

            auto tey_xxx_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 4);

            auto tez_xxx_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 4);

            auto tex_xxx_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 5);

            auto tey_xxx_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 5);

            auto tez_xxx_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 5);

            auto tex_xxy_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 6);

            auto tey_xxy_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 6);

            auto tez_xxy_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 6);

            auto tex_xxy_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 7);

            auto tey_xxy_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 7);

            auto tez_xxy_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 7);

            auto tex_xxy_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 8);

            auto tey_xxy_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 8);

            auto tez_xxy_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 8);

            auto tex_xxy_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 9);

            auto tey_xxy_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 9);

            auto tez_xxy_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 9);

            auto tex_xxy_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 10);

            auto tey_xxy_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 10);

            auto tez_xxy_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 10);

            auto tex_xxy_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 11);

            auto tey_xxy_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 11);

            auto tez_xxy_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 11);

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

            // set up pointers to integrals

            auto tex_xxxx_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx);

            auto tey_xxxx_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx);

            auto tez_xxxx_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx);

            auto tex_xxxx_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 1);

            auto tey_xxxx_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 1);

            auto tez_xxxx_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 1);

            auto tex_xxxx_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 2);

            auto tey_xxxx_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 2);

            auto tez_xxxx_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 2);

            auto tex_xxxx_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 3);

            auto tey_xxxx_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 3);

            auto tez_xxxx_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 3);

            auto tex_xxxx_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 4);

            auto tey_xxxx_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 4);

            auto tez_xxxx_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 4);

            auto tex_xxxx_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 5);

            auto tey_xxxx_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 5);

            auto tez_xxxx_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 5);

            auto tex_xxxx_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 6);

            auto tey_xxxx_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 6);

            auto tez_xxxx_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 6);

            auto tex_xxxx_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 7);

            auto tey_xxxx_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 7);

            auto tez_xxxx_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 7);

            auto tex_xxxx_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 8);

            auto tey_xxxx_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 8);

            auto tez_xxxx_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 8);

            auto tex_xxxx_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 9);

            auto tey_xxxx_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 9);

            auto tez_xxxx_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 9);

            auto tex_xxxy_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 10);

            auto tey_xxxy_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 10);

            auto tez_xxxy_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 10);

            auto tex_xxxy_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 11);

            auto tey_xxxy_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 11);

            auto tez_xxxy_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 11);

            auto tex_xxxy_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 12);

            auto tey_xxxy_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 12);

            auto tez_xxxy_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 12);

            auto tex_xxxy_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 13);

            auto tey_xxxy_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 13);

            auto tez_xxxy_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 13);

            auto tex_xxxy_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 14);

            auto tey_xxxy_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 14);

            auto tez_xxxy_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 14);

            auto tex_xxxy_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 15);

            auto tey_xxxy_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 15);

            auto tez_xxxy_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 15);

            auto tex_xxxy_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 16);

            auto tey_xxxy_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 16);

            // Batch of Integrals (0,50)

            #pragma omp simd aligned(fx, pa_x, pc_x, ta_xxx_xxx_1, ta_xxx_xxy_1, ta_xxx_xxz_1, ta_xxx_xyy_1, \
                                         ta_xxx_xyz_1, ta_xxx_xzz_1, ta_xxx_yyy_1, ta_xxx_yyz_1, ta_xxx_yzz_1, ta_xxx_zzz_1, \
                                         ta_xxy_xxx_1, ta_xxy_xxy_1, ta_xxy_xxz_1, ta_xxy_xyy_1, ta_xxy_xyz_1, ta_xxy_xzz_1, \
                                         ta_xxy_yyy_1, tex_xx_xxx_0, tex_xx_xxx_1, tex_xx_xxy_0, tex_xx_xxy_1, tex_xx_xxz_0, \
                                         tex_xx_xxz_1, tex_xx_xyy_0, tex_xx_xyy_1, tex_xx_xyz_0, tex_xx_xyz_1, tex_xx_xzz_0, \
                                         tex_xx_xzz_1, tex_xx_yyy_0, tex_xx_yyy_1, tex_xx_yyz_0, tex_xx_yyz_1, tex_xx_yzz_0, \
                                         tex_xx_yzz_1, tex_xx_zzz_0, tex_xx_zzz_1, tex_xxx_xx_0, tex_xxx_xx_1, \
                                         tex_xxx_xxx_0, tex_xxx_xxx_1, tex_xxx_xxy_0, tex_xxx_xxy_1, tex_xxx_xxz_0, \
                                         tex_xxx_xxz_1, tex_xxx_xy_0, tex_xxx_xy_1, tex_xxx_xyy_0, tex_xxx_xyy_1, \
                                         tex_xxx_xyz_0, tex_xxx_xyz_1, tex_xxx_xz_0, tex_xxx_xz_1, tex_xxx_xzz_0, \
                                         tex_xxx_xzz_1, tex_xxx_yy_0, tex_xxx_yy_1, tex_xxx_yyy_0, tex_xxx_yyy_1, \
                                         tex_xxx_yyz_0, tex_xxx_yyz_1, tex_xxx_yz_0, tex_xxx_yz_1, tex_xxx_yzz_0, \
                                         tex_xxx_yzz_1, tex_xxx_zz_0, tex_xxx_zz_1, tex_xxx_zzz_0, tex_xxx_zzz_1, \
                                         tex_xxxx_xxx_0, tex_xxxx_xxy_0, tex_xxxx_xxz_0, tex_xxxx_xyy_0, tex_xxxx_xyz_0, \
                                         tex_xxxx_xzz_0, tex_xxxx_yyy_0, tex_xxxx_yyz_0, tex_xxxx_yzz_0, tex_xxxx_zzz_0, \
                                         tex_xxxy_xxx_0, tex_xxxy_xxy_0, tex_xxxy_xxz_0, tex_xxxy_xyy_0, tex_xxxy_xyz_0, \
                                         tex_xxxy_xzz_0, tex_xxxy_yyy_0, tex_xxy_xx_0, tex_xxy_xx_1, tex_xxy_xxx_0, \
                                         tex_xxy_xxx_1, tex_xxy_xxy_0, tex_xxy_xxy_1, tex_xxy_xxz_0, tex_xxy_xxz_1, \
                                         tex_xxy_xy_0, tex_xxy_xy_1, tex_xxy_xyy_0, tex_xxy_xyy_1, tex_xxy_xyz_0, \
                                         tex_xxy_xyz_1, tex_xxy_xz_0, tex_xxy_xz_1, tex_xxy_xzz_0, tex_xxy_xzz_1, \
                                         tex_xxy_yy_0, tex_xxy_yy_1, tex_xxy_yyy_0, tex_xxy_yyy_1, tex_xxy_yz_0, \
                                         tex_xxy_yz_1, tex_xxy_zz_0, tex_xxy_zz_1, tex_xy_xxx_0, tex_xy_xxx_1, tex_xy_xxy_0, \
                                         tex_xy_xxy_1, tex_xy_xxz_0, tex_xy_xxz_1, tex_xy_xyy_0, tex_xy_xyy_1, tex_xy_xyz_0, \
                                         tex_xy_xyz_1, tex_xy_xzz_0, tex_xy_xzz_1, tex_xy_yyy_0, tex_xy_yyy_1, tey_xx_xxx_0, \
                                         tey_xx_xxx_1, tey_xx_xxy_0, tey_xx_xxy_1, tey_xx_xxz_0, tey_xx_xxz_1, tey_xx_xyy_0, \
                                         tey_xx_xyy_1, tey_xx_xyz_0, tey_xx_xyz_1, tey_xx_xzz_0, tey_xx_xzz_1, tey_xx_yyy_0, \
                                         tey_xx_yyy_1, tey_xx_yyz_0, tey_xx_yyz_1, tey_xx_yzz_0, tey_xx_yzz_1, tey_xx_zzz_0, \
                                         tey_xx_zzz_1, tey_xxx_xx_0, tey_xxx_xx_1, tey_xxx_xxx_0, tey_xxx_xxx_1, \
                                         tey_xxx_xxy_0, tey_xxx_xxy_1, tey_xxx_xxz_0, tey_xxx_xxz_1, tey_xxx_xy_0, \
                                         tey_xxx_xy_1, tey_xxx_xyy_0, tey_xxx_xyy_1, tey_xxx_xyz_0, tey_xxx_xyz_1, \
                                         tey_xxx_xz_0, tey_xxx_xz_1, tey_xxx_xzz_0, tey_xxx_xzz_1, tey_xxx_yy_0, \
                                         tey_xxx_yy_1, tey_xxx_yyy_0, tey_xxx_yyy_1, tey_xxx_yyz_0, tey_xxx_yyz_1, \
                                         tey_xxx_yz_0, tey_xxx_yz_1, tey_xxx_yzz_0, tey_xxx_yzz_1, tey_xxx_zz_0, \
                                         tey_xxx_zz_1, tey_xxx_zzz_0, tey_xxx_zzz_1, tey_xxxx_xxx_0, tey_xxxx_xxy_0, \
                                         tey_xxxx_xxz_0, tey_xxxx_xyy_0, tey_xxxx_xyz_0, tey_xxxx_xzz_0, tey_xxxx_yyy_0, \
                                         tey_xxxx_yyz_0, tey_xxxx_yzz_0, tey_xxxx_zzz_0, tey_xxxy_xxx_0, tey_xxxy_xxy_0, \
                                         tey_xxxy_xxz_0, tey_xxxy_xyy_0, tey_xxxy_xyz_0, tey_xxxy_xzz_0, tey_xxxy_yyy_0, \
                                         tey_xxy_xx_0, tey_xxy_xx_1, tey_xxy_xxx_0, tey_xxy_xxx_1, tey_xxy_xxy_0, \
                                         tey_xxy_xxy_1, tey_xxy_xxz_0, tey_xxy_xxz_1, tey_xxy_xy_0, tey_xxy_xy_1, \
                                         tey_xxy_xyy_0, tey_xxy_xyy_1, tey_xxy_xyz_0, tey_xxy_xyz_1, tey_xxy_xz_0, \
                                         tey_xxy_xz_1, tey_xxy_xzz_0, tey_xxy_xzz_1, tey_xxy_yy_0, tey_xxy_yy_1, \
                                         tey_xxy_yyy_0, tey_xxy_yyy_1, tey_xxy_yz_0, tey_xxy_yz_1, tey_xxy_zz_0, \
                                         tey_xxy_zz_1, tey_xy_xxx_0, tey_xy_xxx_1, tey_xy_xxy_0, tey_xy_xxy_1, tey_xy_xxz_0, \
                                         tey_xy_xxz_1, tey_xy_xyy_0, tey_xy_xyy_1, tey_xy_xyz_0, tey_xy_xyz_1, tey_xy_xzz_0, \
                                         tey_xy_xzz_1, tey_xy_yyy_0, tey_xy_yyy_1, tez_xx_xxx_0, tez_xx_xxx_1, tez_xx_xxy_0, \
                                         tez_xx_xxy_1, tez_xx_xxz_0, tez_xx_xxz_1, tez_xx_xyy_0, tez_xx_xyy_1, tez_xx_xyz_0, \
                                         tez_xx_xyz_1, tez_xx_xzz_0, tez_xx_xzz_1, tez_xx_yyy_0, tez_xx_yyy_1, tez_xx_yyz_0, \
                                         tez_xx_yyz_1, tez_xx_yzz_0, tez_xx_yzz_1, tez_xx_zzz_0, tez_xx_zzz_1, tez_xxx_xx_0, \
                                         tez_xxx_xx_1, tez_xxx_xxx_0, tez_xxx_xxx_1, tez_xxx_xxy_0, tez_xxx_xxy_1, \
                                         tez_xxx_xxz_0, tez_xxx_xxz_1, tez_xxx_xy_0, tez_xxx_xy_1, tez_xxx_xyy_0, \
                                         tez_xxx_xyy_1, tez_xxx_xyz_0, tez_xxx_xyz_1, tez_xxx_xz_0, tez_xxx_xz_1, \
                                         tez_xxx_xzz_0, tez_xxx_xzz_1, tez_xxx_yy_0, tez_xxx_yy_1, tez_xxx_yyy_0, \
                                         tez_xxx_yyy_1, tez_xxx_yyz_0, tez_xxx_yyz_1, tez_xxx_yz_0, tez_xxx_yz_1, \
                                         tez_xxx_yzz_0, tez_xxx_yzz_1, tez_xxx_zz_0, tez_xxx_zz_1, tez_xxx_zzz_0, \
                                         tez_xxx_zzz_1, tez_xxxx_xxx_0, tez_xxxx_xxy_0, tez_xxxx_xxz_0, tez_xxxx_xyy_0, \
                                         tez_xxxx_xyz_0, tez_xxxx_xzz_0, tez_xxxx_yyy_0, tez_xxxx_yyz_0, tez_xxxx_yzz_0, \
                                         tez_xxxx_zzz_0, tez_xxxy_xxx_0, tez_xxxy_xxy_0, tez_xxxy_xxz_0, tez_xxxy_xyy_0, \
                                         tez_xxxy_xyz_0, tez_xxxy_xzz_0, tez_xxy_xx_0, tez_xxy_xx_1, tez_xxy_xxx_0, \
                                         tez_xxy_xxx_1, tez_xxy_xxy_0, tez_xxy_xxy_1, tez_xxy_xxz_0, tez_xxy_xxz_1, \
                                         tez_xxy_xy_0, tez_xxy_xy_1, tez_xxy_xyy_0, tez_xxy_xyy_1, tez_xxy_xyz_0, \
                                         tez_xxy_xyz_1, tez_xxy_xz_0, tez_xxy_xz_1, tez_xxy_xzz_0, tez_xxy_xzz_1, \
                                         tez_xxy_yy_0, tez_xxy_yy_1, tez_xxy_yz_0, tez_xxy_yz_1, tez_xxy_zz_0, tez_xxy_zz_1, \
                                         tez_xy_xxx_0, tez_xy_xxx_1, tez_xy_xxy_0, tez_xy_xxy_1, tez_xy_xxz_0, tez_xy_xxz_1, \
                                         tez_xy_xyy_0, tez_xy_xyy_1, tez_xy_xyz_0, tez_xy_xyz_1, tez_xy_xzz_0, tez_xy_xzz_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_xxxx_xxx_0[j] = pa_x[j] * tex_xxx_xxx_0[j] - pc_x[j] * tex_xxx_xxx_1[j] + 1.5 * fl1_fx * tex_xx_xxx_0[j] -
                                    1.5 * fl1_fx * tex_xx_xxx_1[j] + 1.5 * fl1_fx * tex_xxx_xx_0[j] - 1.5 * fl1_fx * tex_xxx_xx_1[j] +
                                    ta_xxx_xxx_1[j];

                tey_xxxx_xxx_0[j] = pa_x[j] * tey_xxx_xxx_0[j] - pc_x[j] * tey_xxx_xxx_1[j] + 1.5 * fl1_fx * tey_xx_xxx_0[j] -
                                    1.5 * fl1_fx * tey_xx_xxx_1[j] + 1.5 * fl1_fx * tey_xxx_xx_0[j] - 1.5 * fl1_fx * tey_xxx_xx_1[j];

                tez_xxxx_xxx_0[j] = pa_x[j] * tez_xxx_xxx_0[j] - pc_x[j] * tez_xxx_xxx_1[j] + 1.5 * fl1_fx * tez_xx_xxx_0[j] -
                                    1.5 * fl1_fx * tez_xx_xxx_1[j] + 1.5 * fl1_fx * tez_xxx_xx_0[j] - 1.5 * fl1_fx * tez_xxx_xx_1[j];

                tex_xxxx_xxy_0[j] = pa_x[j] * tex_xxx_xxy_0[j] - pc_x[j] * tex_xxx_xxy_1[j] + 1.5 * fl1_fx * tex_xx_xxy_0[j] -
                                    1.5 * fl1_fx * tex_xx_xxy_1[j] + fl1_fx * tex_xxx_xy_0[j] - fl1_fx * tex_xxx_xy_1[j] + ta_xxx_xxy_1[j];

                tey_xxxx_xxy_0[j] = pa_x[j] * tey_xxx_xxy_0[j] - pc_x[j] * tey_xxx_xxy_1[j] + 1.5 * fl1_fx * tey_xx_xxy_0[j] -
                                    1.5 * fl1_fx * tey_xx_xxy_1[j] + fl1_fx * tey_xxx_xy_0[j] - fl1_fx * tey_xxx_xy_1[j];

                tez_xxxx_xxy_0[j] = pa_x[j] * tez_xxx_xxy_0[j] - pc_x[j] * tez_xxx_xxy_1[j] + 1.5 * fl1_fx * tez_xx_xxy_0[j] -
                                    1.5 * fl1_fx * tez_xx_xxy_1[j] + fl1_fx * tez_xxx_xy_0[j] - fl1_fx * tez_xxx_xy_1[j];

                tex_xxxx_xxz_0[j] = pa_x[j] * tex_xxx_xxz_0[j] - pc_x[j] * tex_xxx_xxz_1[j] + 1.5 * fl1_fx * tex_xx_xxz_0[j] -
                                    1.5 * fl1_fx * tex_xx_xxz_1[j] + fl1_fx * tex_xxx_xz_0[j] - fl1_fx * tex_xxx_xz_1[j] + ta_xxx_xxz_1[j];

                tey_xxxx_xxz_0[j] = pa_x[j] * tey_xxx_xxz_0[j] - pc_x[j] * tey_xxx_xxz_1[j] + 1.5 * fl1_fx * tey_xx_xxz_0[j] -
                                    1.5 * fl1_fx * tey_xx_xxz_1[j] + fl1_fx * tey_xxx_xz_0[j] - fl1_fx * tey_xxx_xz_1[j];

                tez_xxxx_xxz_0[j] = pa_x[j] * tez_xxx_xxz_0[j] - pc_x[j] * tez_xxx_xxz_1[j] + 1.5 * fl1_fx * tez_xx_xxz_0[j] -
                                    1.5 * fl1_fx * tez_xx_xxz_1[j] + fl1_fx * tez_xxx_xz_0[j] - fl1_fx * tez_xxx_xz_1[j];

                tex_xxxx_xyy_0[j] = pa_x[j] * tex_xxx_xyy_0[j] - pc_x[j] * tex_xxx_xyy_1[j] + 1.5 * fl1_fx * tex_xx_xyy_0[j] -
                                    1.5 * fl1_fx * tex_xx_xyy_1[j] + 0.5 * fl1_fx * tex_xxx_yy_0[j] - 0.5 * fl1_fx * tex_xxx_yy_1[j] +
                                    ta_xxx_xyy_1[j];

                tey_xxxx_xyy_0[j] = pa_x[j] * tey_xxx_xyy_0[j] - pc_x[j] * tey_xxx_xyy_1[j] + 1.5 * fl1_fx * tey_xx_xyy_0[j] -
                                    1.5 * fl1_fx * tey_xx_xyy_1[j] + 0.5 * fl1_fx * tey_xxx_yy_0[j] - 0.5 * fl1_fx * tey_xxx_yy_1[j];

                tez_xxxx_xyy_0[j] = pa_x[j] * tez_xxx_xyy_0[j] - pc_x[j] * tez_xxx_xyy_1[j] + 1.5 * fl1_fx * tez_xx_xyy_0[j] -
                                    1.5 * fl1_fx * tez_xx_xyy_1[j] + 0.5 * fl1_fx * tez_xxx_yy_0[j] - 0.5 * fl1_fx * tez_xxx_yy_1[j];

                tex_xxxx_xyz_0[j] = pa_x[j] * tex_xxx_xyz_0[j] - pc_x[j] * tex_xxx_xyz_1[j] + 1.5 * fl1_fx * tex_xx_xyz_0[j] -
                                    1.5 * fl1_fx * tex_xx_xyz_1[j] + 0.5 * fl1_fx * tex_xxx_yz_0[j] - 0.5 * fl1_fx * tex_xxx_yz_1[j] +
                                    ta_xxx_xyz_1[j];

                tey_xxxx_xyz_0[j] = pa_x[j] * tey_xxx_xyz_0[j] - pc_x[j] * tey_xxx_xyz_1[j] + 1.5 * fl1_fx * tey_xx_xyz_0[j] -
                                    1.5 * fl1_fx * tey_xx_xyz_1[j] + 0.5 * fl1_fx * tey_xxx_yz_0[j] - 0.5 * fl1_fx * tey_xxx_yz_1[j];

                tez_xxxx_xyz_0[j] = pa_x[j] * tez_xxx_xyz_0[j] - pc_x[j] * tez_xxx_xyz_1[j] + 1.5 * fl1_fx * tez_xx_xyz_0[j] -
                                    1.5 * fl1_fx * tez_xx_xyz_1[j] + 0.5 * fl1_fx * tez_xxx_yz_0[j] - 0.5 * fl1_fx * tez_xxx_yz_1[j];

                tex_xxxx_xzz_0[j] = pa_x[j] * tex_xxx_xzz_0[j] - pc_x[j] * tex_xxx_xzz_1[j] + 1.5 * fl1_fx * tex_xx_xzz_0[j] -
                                    1.5 * fl1_fx * tex_xx_xzz_1[j] + 0.5 * fl1_fx * tex_xxx_zz_0[j] - 0.5 * fl1_fx * tex_xxx_zz_1[j] +
                                    ta_xxx_xzz_1[j];

                tey_xxxx_xzz_0[j] = pa_x[j] * tey_xxx_xzz_0[j] - pc_x[j] * tey_xxx_xzz_1[j] + 1.5 * fl1_fx * tey_xx_xzz_0[j] -
                                    1.5 * fl1_fx * tey_xx_xzz_1[j] + 0.5 * fl1_fx * tey_xxx_zz_0[j] - 0.5 * fl1_fx * tey_xxx_zz_1[j];

                tez_xxxx_xzz_0[j] = pa_x[j] * tez_xxx_xzz_0[j] - pc_x[j] * tez_xxx_xzz_1[j] + 1.5 * fl1_fx * tez_xx_xzz_0[j] -
                                    1.5 * fl1_fx * tez_xx_xzz_1[j] + 0.5 * fl1_fx * tez_xxx_zz_0[j] - 0.5 * fl1_fx * tez_xxx_zz_1[j];

                tex_xxxx_yyy_0[j] = pa_x[j] * tex_xxx_yyy_0[j] - pc_x[j] * tex_xxx_yyy_1[j] + 1.5 * fl1_fx * tex_xx_yyy_0[j] -
                                    1.5 * fl1_fx * tex_xx_yyy_1[j] + ta_xxx_yyy_1[j];

                tey_xxxx_yyy_0[j] =
                    pa_x[j] * tey_xxx_yyy_0[j] - pc_x[j] * tey_xxx_yyy_1[j] + 1.5 * fl1_fx * tey_xx_yyy_0[j] - 1.5 * fl1_fx * tey_xx_yyy_1[j];

                tez_xxxx_yyy_0[j] =
                    pa_x[j] * tez_xxx_yyy_0[j] - pc_x[j] * tez_xxx_yyy_1[j] + 1.5 * fl1_fx * tez_xx_yyy_0[j] - 1.5 * fl1_fx * tez_xx_yyy_1[j];

                tex_xxxx_yyz_0[j] = pa_x[j] * tex_xxx_yyz_0[j] - pc_x[j] * tex_xxx_yyz_1[j] + 1.5 * fl1_fx * tex_xx_yyz_0[j] -
                                    1.5 * fl1_fx * tex_xx_yyz_1[j] + ta_xxx_yyz_1[j];

                tey_xxxx_yyz_0[j] =
                    pa_x[j] * tey_xxx_yyz_0[j] - pc_x[j] * tey_xxx_yyz_1[j] + 1.5 * fl1_fx * tey_xx_yyz_0[j] - 1.5 * fl1_fx * tey_xx_yyz_1[j];

                tez_xxxx_yyz_0[j] =
                    pa_x[j] * tez_xxx_yyz_0[j] - pc_x[j] * tez_xxx_yyz_1[j] + 1.5 * fl1_fx * tez_xx_yyz_0[j] - 1.5 * fl1_fx * tez_xx_yyz_1[j];

                tex_xxxx_yzz_0[j] = pa_x[j] * tex_xxx_yzz_0[j] - pc_x[j] * tex_xxx_yzz_1[j] + 1.5 * fl1_fx * tex_xx_yzz_0[j] -
                                    1.5 * fl1_fx * tex_xx_yzz_1[j] + ta_xxx_yzz_1[j];

                tey_xxxx_yzz_0[j] =
                    pa_x[j] * tey_xxx_yzz_0[j] - pc_x[j] * tey_xxx_yzz_1[j] + 1.5 * fl1_fx * tey_xx_yzz_0[j] - 1.5 * fl1_fx * tey_xx_yzz_1[j];

                tez_xxxx_yzz_0[j] =
                    pa_x[j] * tez_xxx_yzz_0[j] - pc_x[j] * tez_xxx_yzz_1[j] + 1.5 * fl1_fx * tez_xx_yzz_0[j] - 1.5 * fl1_fx * tez_xx_yzz_1[j];

                tex_xxxx_zzz_0[j] = pa_x[j] * tex_xxx_zzz_0[j] - pc_x[j] * tex_xxx_zzz_1[j] + 1.5 * fl1_fx * tex_xx_zzz_0[j] -
                                    1.5 * fl1_fx * tex_xx_zzz_1[j] + ta_xxx_zzz_1[j];

                tey_xxxx_zzz_0[j] =
                    pa_x[j] * tey_xxx_zzz_0[j] - pc_x[j] * tey_xxx_zzz_1[j] + 1.5 * fl1_fx * tey_xx_zzz_0[j] - 1.5 * fl1_fx * tey_xx_zzz_1[j];

                tez_xxxx_zzz_0[j] =
                    pa_x[j] * tez_xxx_zzz_0[j] - pc_x[j] * tez_xxx_zzz_1[j] + 1.5 * fl1_fx * tez_xx_zzz_0[j] - 1.5 * fl1_fx * tez_xx_zzz_1[j];

                tex_xxxy_xxx_0[j] = pa_x[j] * tex_xxy_xxx_0[j] - pc_x[j] * tex_xxy_xxx_1[j] + fl1_fx * tex_xy_xxx_0[j] - fl1_fx * tex_xy_xxx_1[j] +
                                    1.5 * fl1_fx * tex_xxy_xx_0[j] - 1.5 * fl1_fx * tex_xxy_xx_1[j] + ta_xxy_xxx_1[j];

                tey_xxxy_xxx_0[j] = pa_x[j] * tey_xxy_xxx_0[j] - pc_x[j] * tey_xxy_xxx_1[j] + fl1_fx * tey_xy_xxx_0[j] - fl1_fx * tey_xy_xxx_1[j] +
                                    1.5 * fl1_fx * tey_xxy_xx_0[j] - 1.5 * fl1_fx * tey_xxy_xx_1[j];

                tez_xxxy_xxx_0[j] = pa_x[j] * tez_xxy_xxx_0[j] - pc_x[j] * tez_xxy_xxx_1[j] + fl1_fx * tez_xy_xxx_0[j] - fl1_fx * tez_xy_xxx_1[j] +
                                    1.5 * fl1_fx * tez_xxy_xx_0[j] - 1.5 * fl1_fx * tez_xxy_xx_1[j];

                tex_xxxy_xxy_0[j] = pa_x[j] * tex_xxy_xxy_0[j] - pc_x[j] * tex_xxy_xxy_1[j] + fl1_fx * tex_xy_xxy_0[j] - fl1_fx * tex_xy_xxy_1[j] +
                                    fl1_fx * tex_xxy_xy_0[j] - fl1_fx * tex_xxy_xy_1[j] + ta_xxy_xxy_1[j];

                tey_xxxy_xxy_0[j] = pa_x[j] * tey_xxy_xxy_0[j] - pc_x[j] * tey_xxy_xxy_1[j] + fl1_fx * tey_xy_xxy_0[j] - fl1_fx * tey_xy_xxy_1[j] +
                                    fl1_fx * tey_xxy_xy_0[j] - fl1_fx * tey_xxy_xy_1[j];

                tez_xxxy_xxy_0[j] = pa_x[j] * tez_xxy_xxy_0[j] - pc_x[j] * tez_xxy_xxy_1[j] + fl1_fx * tez_xy_xxy_0[j] - fl1_fx * tez_xy_xxy_1[j] +
                                    fl1_fx * tez_xxy_xy_0[j] - fl1_fx * tez_xxy_xy_1[j];

                tex_xxxy_xxz_0[j] = pa_x[j] * tex_xxy_xxz_0[j] - pc_x[j] * tex_xxy_xxz_1[j] + fl1_fx * tex_xy_xxz_0[j] - fl1_fx * tex_xy_xxz_1[j] +
                                    fl1_fx * tex_xxy_xz_0[j] - fl1_fx * tex_xxy_xz_1[j] + ta_xxy_xxz_1[j];

                tey_xxxy_xxz_0[j] = pa_x[j] * tey_xxy_xxz_0[j] - pc_x[j] * tey_xxy_xxz_1[j] + fl1_fx * tey_xy_xxz_0[j] - fl1_fx * tey_xy_xxz_1[j] +
                                    fl1_fx * tey_xxy_xz_0[j] - fl1_fx * tey_xxy_xz_1[j];

                tez_xxxy_xxz_0[j] = pa_x[j] * tez_xxy_xxz_0[j] - pc_x[j] * tez_xxy_xxz_1[j] + fl1_fx * tez_xy_xxz_0[j] - fl1_fx * tez_xy_xxz_1[j] +
                                    fl1_fx * tez_xxy_xz_0[j] - fl1_fx * tez_xxy_xz_1[j];

                tex_xxxy_xyy_0[j] = pa_x[j] * tex_xxy_xyy_0[j] - pc_x[j] * tex_xxy_xyy_1[j] + fl1_fx * tex_xy_xyy_0[j] - fl1_fx * tex_xy_xyy_1[j] +
                                    0.5 * fl1_fx * tex_xxy_yy_0[j] - 0.5 * fl1_fx * tex_xxy_yy_1[j] + ta_xxy_xyy_1[j];

                tey_xxxy_xyy_0[j] = pa_x[j] * tey_xxy_xyy_0[j] - pc_x[j] * tey_xxy_xyy_1[j] + fl1_fx * tey_xy_xyy_0[j] - fl1_fx * tey_xy_xyy_1[j] +
                                    0.5 * fl1_fx * tey_xxy_yy_0[j] - 0.5 * fl1_fx * tey_xxy_yy_1[j];

                tez_xxxy_xyy_0[j] = pa_x[j] * tez_xxy_xyy_0[j] - pc_x[j] * tez_xxy_xyy_1[j] + fl1_fx * tez_xy_xyy_0[j] - fl1_fx * tez_xy_xyy_1[j] +
                                    0.5 * fl1_fx * tez_xxy_yy_0[j] - 0.5 * fl1_fx * tez_xxy_yy_1[j];

                tex_xxxy_xyz_0[j] = pa_x[j] * tex_xxy_xyz_0[j] - pc_x[j] * tex_xxy_xyz_1[j] + fl1_fx * tex_xy_xyz_0[j] - fl1_fx * tex_xy_xyz_1[j] +
                                    0.5 * fl1_fx * tex_xxy_yz_0[j] - 0.5 * fl1_fx * tex_xxy_yz_1[j] + ta_xxy_xyz_1[j];

                tey_xxxy_xyz_0[j] = pa_x[j] * tey_xxy_xyz_0[j] - pc_x[j] * tey_xxy_xyz_1[j] + fl1_fx * tey_xy_xyz_0[j] - fl1_fx * tey_xy_xyz_1[j] +
                                    0.5 * fl1_fx * tey_xxy_yz_0[j] - 0.5 * fl1_fx * tey_xxy_yz_1[j];

                tez_xxxy_xyz_0[j] = pa_x[j] * tez_xxy_xyz_0[j] - pc_x[j] * tez_xxy_xyz_1[j] + fl1_fx * tez_xy_xyz_0[j] - fl1_fx * tez_xy_xyz_1[j] +
                                    0.5 * fl1_fx * tez_xxy_yz_0[j] - 0.5 * fl1_fx * tez_xxy_yz_1[j];

                tex_xxxy_xzz_0[j] = pa_x[j] * tex_xxy_xzz_0[j] - pc_x[j] * tex_xxy_xzz_1[j] + fl1_fx * tex_xy_xzz_0[j] - fl1_fx * tex_xy_xzz_1[j] +
                                    0.5 * fl1_fx * tex_xxy_zz_0[j] - 0.5 * fl1_fx * tex_xxy_zz_1[j] + ta_xxy_xzz_1[j];

                tey_xxxy_xzz_0[j] = pa_x[j] * tey_xxy_xzz_0[j] - pc_x[j] * tey_xxy_xzz_1[j] + fl1_fx * tey_xy_xzz_0[j] - fl1_fx * tey_xy_xzz_1[j] +
                                    0.5 * fl1_fx * tey_xxy_zz_0[j] - 0.5 * fl1_fx * tey_xxy_zz_1[j];

                tez_xxxy_xzz_0[j] = pa_x[j] * tez_xxy_xzz_0[j] - pc_x[j] * tez_xxy_xzz_1[j] + fl1_fx * tez_xy_xzz_0[j] - fl1_fx * tez_xy_xzz_1[j] +
                                    0.5 * fl1_fx * tez_xxy_zz_0[j] - 0.5 * fl1_fx * tez_xxy_zz_1[j];

                tex_xxxy_yyy_0[j] =
                    pa_x[j] * tex_xxy_yyy_0[j] - pc_x[j] * tex_xxy_yyy_1[j] + fl1_fx * tex_xy_yyy_0[j] - fl1_fx * tex_xy_yyy_1[j] + ta_xxy_yyy_1[j];

                tey_xxxy_yyy_0[j] = pa_x[j] * tey_xxy_yyy_0[j] - pc_x[j] * tey_xxy_yyy_1[j] + fl1_fx * tey_xy_yyy_0[j] - fl1_fx * tey_xy_yyy_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForGF_50_100(CMemBlock2D<double>&       primBuffer,
                              const CRecursionMap&       recursionMap,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CMemBlock2D<double>& pcDistances,
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

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_4_3_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_2_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tez_xy_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 16);

            auto tex_xy_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 17);

            auto tey_xy_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 17);

            auto tez_xy_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 17);

            auto tex_xy_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 18);

            auto tey_xy_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 18);

            auto tez_xy_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 18);

            auto tex_xy_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 19);

            auto tey_xy_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 19);

            auto tez_xy_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 19);

            auto tex_xz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 20);

            auto tey_xz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 20);

            auto tez_xz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 20);

            auto tex_xz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 21);

            auto tey_xz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 21);

            auto tez_xz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 21);

            auto tex_xz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 22);

            auto tey_xz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 22);

            auto tez_xz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 22);

            auto tex_xz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 23);

            auto tey_xz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 23);

            auto tez_xz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 23);

            auto tex_xz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 24);

            auto tey_xz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 24);

            auto tez_xz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 24);

            auto tex_xz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 25);

            auto tey_xz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 25);

            auto tez_xz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 25);

            auto tex_xz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 26);

            auto tey_xz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 26);

            auto tez_xz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 26);

            auto tex_xz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 27);

            auto tey_xz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 27);

            auto tez_xz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 27);

            auto tex_xz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 28);

            auto tey_xz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 28);

            auto tez_xz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 28);

            auto tex_xz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 29);

            auto tey_xz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 29);

            auto tez_xz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 29);

            auto tex_yy_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 30);

            auto tey_yy_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 30);

            auto tez_yy_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 30);

            auto tex_yy_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 31);

            auto tey_yy_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 31);

            auto tez_yy_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 31);

            auto tex_yy_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 32);

            auto tey_yy_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 32);

            auto tez_yy_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 32);

            auto tex_yy_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 33);

            auto tez_xy_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 16);

            auto tex_xy_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 17);

            auto tey_xy_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 17);

            auto tez_xy_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 17);

            auto tex_xy_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 18);

            auto tey_xy_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 18);

            auto tez_xy_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 18);

            auto tex_xy_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 19);

            auto tey_xy_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 19);

            auto tez_xy_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 19);

            auto tex_xz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 20);

            auto tey_xz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 20);

            auto tez_xz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 20);

            auto tex_xz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 21);

            auto tey_xz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 21);

            auto tez_xz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 21);

            auto tex_xz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 22);

            auto tey_xz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 22);

            auto tez_xz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 22);

            auto tex_xz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 23);

            auto tey_xz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 23);

            auto tez_xz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 23);

            auto tex_xz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 24);

            auto tey_xz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 24);

            auto tez_xz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 24);

            auto tex_xz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 25);

            auto tey_xz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 25);

            auto tez_xz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 25);

            auto tex_xz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 26);

            auto tey_xz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 26);

            auto tez_xz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 26);

            auto tex_xz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 27);

            auto tey_xz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 27);

            auto tez_xz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 27);

            auto tex_xz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 28);

            auto tey_xz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 28);

            auto tez_xz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 28);

            auto tex_xz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 29);

            auto tey_xz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 29);

            auto tez_xz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 29);

            auto tex_yy_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 30);

            auto tey_yy_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 30);

            auto tez_yy_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 30);

            auto tex_yy_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 31);

            auto tey_yy_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 31);

            auto tez_yy_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 31);

            auto tex_yy_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 32);

            auto tey_yy_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 32);

            auto tez_yy_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 32);

            auto tex_yy_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 33);

            auto tex_xxz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 12);

            auto tey_xxz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 12);

            auto tez_xxz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 12);

            auto tex_xxz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 13);

            auto tey_xxz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 13);

            auto tez_xxz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 13);

            auto tex_xxz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 14);

            auto tey_xxz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 14);

            auto tez_xxz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 14);

            auto tex_xxz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 15);

            auto tey_xxz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 15);

            auto tez_xxz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 15);

            auto tex_xxz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 16);

            auto tey_xxz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 16);

            auto tez_xxz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 16);

            auto tex_xxz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 17);

            auto tey_xxz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 17);

            auto tez_xxz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 17);

            auto tex_xyy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 18);

            auto tey_xyy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 18);

            auto tez_xyy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 18);

            auto tex_xyy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 19);

            auto tey_xyy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 19);

            auto tez_xyy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 19);

            auto tex_xyy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 20);

            auto tey_xyy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 20);

            auto tez_xyy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 20);

            auto tex_xyy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 21);

            auto tex_xxz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 12);

            auto tey_xxz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 12);

            auto tez_xxz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 12);

            auto tex_xxz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 13);

            auto tey_xxz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 13);

            auto tez_xxz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 13);

            auto tex_xxz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 14);

            auto tey_xxz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 14);

            auto tez_xxz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 14);

            auto tex_xxz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 15);

            auto tey_xxz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 15);

            auto tez_xxz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 15);

            auto tex_xxz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 16);

            auto tey_xxz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 16);

            auto tez_xxz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 16);

            auto tex_xxz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 17);

            auto tey_xxz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 17);

            auto tez_xxz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 17);

            auto tex_xyy_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 18);

            auto tey_xyy_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 18);

            auto tez_xyy_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 18);

            auto tex_xyy_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 19);

            auto tey_xyy_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 19);

            auto tez_xyy_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 19);

            auto tex_xyy_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 20);

            auto tey_xyy_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 20);

            auto tez_xyy_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 20);

            auto tex_xyy_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 21);

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

            auto ta_xyy_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 30);

            auto ta_xyy_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 31);

            auto ta_xyy_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 32);

            auto ta_xyy_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 33);

            // set up pointers to integrals

            auto tez_xxxy_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 16);

            auto tex_xxxy_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 17);

            auto tey_xxxy_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 17);

            auto tez_xxxy_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 17);

            auto tex_xxxy_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 18);

            auto tey_xxxy_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 18);

            auto tez_xxxy_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 18);

            auto tex_xxxy_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 19);

            auto tey_xxxy_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 19);

            auto tez_xxxy_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 19);

            auto tex_xxxz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 20);

            auto tey_xxxz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 20);

            auto tez_xxxz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 20);

            auto tex_xxxz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 21);

            auto tey_xxxz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 21);

            auto tez_xxxz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 21);

            auto tex_xxxz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 22);

            auto tey_xxxz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 22);

            auto tez_xxxz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 22);

            auto tex_xxxz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 23);

            auto tey_xxxz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 23);

            auto tez_xxxz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 23);

            auto tex_xxxz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 24);

            auto tey_xxxz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 24);

            auto tez_xxxz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 24);

            auto tex_xxxz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 25);

            auto tey_xxxz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 25);

            auto tez_xxxz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 25);

            auto tex_xxxz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 26);

            auto tey_xxxz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 26);

            auto tez_xxxz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 26);

            auto tex_xxxz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 27);

            auto tey_xxxz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 27);

            auto tez_xxxz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 27);

            auto tex_xxxz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 28);

            auto tey_xxxz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 28);

            auto tez_xxxz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 28);

            auto tex_xxxz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 29);

            auto tey_xxxz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 29);

            auto tez_xxxz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 29);

            auto tex_xxyy_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 30);

            auto tey_xxyy_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 30);

            auto tez_xxyy_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 30);

            auto tex_xxyy_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 31);

            auto tey_xxyy_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 31);

            auto tez_xxyy_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 31);

            auto tex_xxyy_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 32);

            auto tey_xxyy_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 32);

            auto tez_xxyy_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 32);

            auto tex_xxyy_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 33);

            // Batch of Integrals (50,100)

            #pragma omp simd aligned(fx, pa_x, pc_x, ta_xxy_yyz_1, ta_xxy_yzz_1, ta_xxy_zzz_1, ta_xxz_xxx_1, \
                                         ta_xxz_xxy_1, ta_xxz_xxz_1, ta_xxz_xyy_1, ta_xxz_xyz_1, ta_xxz_xzz_1, ta_xxz_yyy_1, \
                                         ta_xxz_yyz_1, ta_xxz_yzz_1, ta_xxz_zzz_1, ta_xyy_xxx_1, ta_xyy_xxy_1, ta_xyy_xxz_1, \
                                         ta_xyy_xyy_1, tex_xxxy_yyz_0, tex_xxxy_yzz_0, tex_xxxy_zzz_0, tex_xxxz_xxx_0, \
                                         tex_xxxz_xxy_0, tex_xxxz_xxz_0, tex_xxxz_xyy_0, tex_xxxz_xyz_0, tex_xxxz_xzz_0, \
                                         tex_xxxz_yyy_0, tex_xxxz_yyz_0, tex_xxxz_yzz_0, tex_xxxz_zzz_0, tex_xxy_yyz_0, \
                                         tex_xxy_yyz_1, tex_xxy_yzz_0, tex_xxy_yzz_1, tex_xxy_zzz_0, tex_xxy_zzz_1, \
                                         tex_xxyy_xxx_0, tex_xxyy_xxy_0, tex_xxyy_xxz_0, tex_xxyy_xyy_0, tex_xxz_xx_0, \
                                         tex_xxz_xx_1, tex_xxz_xxx_0, tex_xxz_xxx_1, tex_xxz_xxy_0, tex_xxz_xxy_1, \
                                         tex_xxz_xxz_0, tex_xxz_xxz_1, tex_xxz_xy_0, tex_xxz_xy_1, tex_xxz_xyy_0, \
                                         tex_xxz_xyy_1, tex_xxz_xyz_0, tex_xxz_xyz_1, tex_xxz_xz_0, tex_xxz_xz_1, \
                                         tex_xxz_xzz_0, tex_xxz_xzz_1, tex_xxz_yy_0, tex_xxz_yy_1, tex_xxz_yyy_0, \
                                         tex_xxz_yyy_1, tex_xxz_yyz_0, tex_xxz_yyz_1, tex_xxz_yz_0, tex_xxz_yz_1, \
                                         tex_xxz_yzz_0, tex_xxz_yzz_1, tex_xxz_zz_0, tex_xxz_zz_1, tex_xxz_zzz_0, \
                                         tex_xxz_zzz_1, tex_xy_yyz_0, tex_xy_yyz_1, tex_xy_yzz_0, tex_xy_yzz_1, tex_xy_zzz_0, \
                                         tex_xy_zzz_1, tex_xyy_xx_0, tex_xyy_xx_1, tex_xyy_xxx_0, tex_xyy_xxx_1, \
                                         tex_xyy_xxy_0, tex_xyy_xxy_1, tex_xyy_xxz_0, tex_xyy_xxz_1, tex_xyy_xy_0, \
                                         tex_xyy_xy_1, tex_xyy_xyy_0, tex_xyy_xyy_1, tex_xyy_xz_0, tex_xyy_xz_1, \
                                         tex_xyy_yy_0, tex_xyy_yy_1, tex_xz_xxx_0, tex_xz_xxx_1, tex_xz_xxy_0, tex_xz_xxy_1, \
                                         tex_xz_xxz_0, tex_xz_xxz_1, tex_xz_xyy_0, tex_xz_xyy_1, tex_xz_xyz_0, tex_xz_xyz_1, \
                                         tex_xz_xzz_0, tex_xz_xzz_1, tex_xz_yyy_0, tex_xz_yyy_1, tex_xz_yyz_0, tex_xz_yyz_1, \
                                         tex_xz_yzz_0, tex_xz_yzz_1, tex_xz_zzz_0, tex_xz_zzz_1, tex_yy_xxx_0, tex_yy_xxx_1, \
                                         tex_yy_xxy_0, tex_yy_xxy_1, tex_yy_xxz_0, tex_yy_xxz_1, tex_yy_xyy_0, tex_yy_xyy_1, \
                                         tey_xxxy_yyz_0, tey_xxxy_yzz_0, tey_xxxy_zzz_0, tey_xxxz_xxx_0, tey_xxxz_xxy_0, \
                                         tey_xxxz_xxz_0, tey_xxxz_xyy_0, tey_xxxz_xyz_0, tey_xxxz_xzz_0, tey_xxxz_yyy_0, \
                                         tey_xxxz_yyz_0, tey_xxxz_yzz_0, tey_xxxz_zzz_0, tey_xxy_yyz_0, tey_xxy_yyz_1, \
                                         tey_xxy_yzz_0, tey_xxy_yzz_1, tey_xxy_zzz_0, tey_xxy_zzz_1, tey_xxyy_xxx_0, \
                                         tey_xxyy_xxy_0, tey_xxyy_xxz_0, tey_xxz_xx_0, tey_xxz_xx_1, tey_xxz_xxx_0, \
                                         tey_xxz_xxx_1, tey_xxz_xxy_0, tey_xxz_xxy_1, tey_xxz_xxz_0, tey_xxz_xxz_1, \
                                         tey_xxz_xy_0, tey_xxz_xy_1, tey_xxz_xyy_0, tey_xxz_xyy_1, tey_xxz_xyz_0, \
                                         tey_xxz_xyz_1, tey_xxz_xz_0, tey_xxz_xz_1, tey_xxz_xzz_0, tey_xxz_xzz_1, \
                                         tey_xxz_yy_0, tey_xxz_yy_1, tey_xxz_yyy_0, tey_xxz_yyy_1, tey_xxz_yyz_0, \
                                         tey_xxz_yyz_1, tey_xxz_yz_0, tey_xxz_yz_1, tey_xxz_yzz_0, tey_xxz_yzz_1, \
                                         tey_xxz_zz_0, tey_xxz_zz_1, tey_xxz_zzz_0, tey_xxz_zzz_1, tey_xy_yyz_0, \
                                         tey_xy_yyz_1, tey_xy_yzz_0, tey_xy_yzz_1, tey_xy_zzz_0, tey_xy_zzz_1, tey_xyy_xx_0, \
                                         tey_xyy_xx_1, tey_xyy_xxx_0, tey_xyy_xxx_1, tey_xyy_xxy_0, tey_xyy_xxy_1, \
                                         tey_xyy_xxz_0, tey_xyy_xxz_1, tey_xyy_xy_0, tey_xyy_xy_1, tey_xyy_xz_0, \
                                         tey_xyy_xz_1, tey_xz_xxx_0, tey_xz_xxx_1, tey_xz_xxy_0, tey_xz_xxy_1, tey_xz_xxz_0, \
                                         tey_xz_xxz_1, tey_xz_xyy_0, tey_xz_xyy_1, tey_xz_xyz_0, tey_xz_xyz_1, tey_xz_xzz_0, \
                                         tey_xz_xzz_1, tey_xz_yyy_0, tey_xz_yyy_1, tey_xz_yyz_0, tey_xz_yyz_1, tey_xz_yzz_0, \
                                         tey_xz_yzz_1, tey_xz_zzz_0, tey_xz_zzz_1, tey_yy_xxx_0, tey_yy_xxx_1, tey_yy_xxy_0, \
                                         tey_yy_xxy_1, tey_yy_xxz_0, tey_yy_xxz_1, tez_xxxy_yyy_0, tez_xxxy_yyz_0, \
                                         tez_xxxy_yzz_0, tez_xxxy_zzz_0, tez_xxxz_xxx_0, tez_xxxz_xxy_0, tez_xxxz_xxz_0, \
                                         tez_xxxz_xyy_0, tez_xxxz_xyz_0, tez_xxxz_xzz_0, tez_xxxz_yyy_0, tez_xxxz_yyz_0, \
                                         tez_xxxz_yzz_0, tez_xxxz_zzz_0, tez_xxy_yyy_0, tez_xxy_yyy_1, tez_xxy_yyz_0, \
                                         tez_xxy_yyz_1, tez_xxy_yzz_0, tez_xxy_yzz_1, tez_xxy_zzz_0, tez_xxy_zzz_1, \
                                         tez_xxyy_xxx_0, tez_xxyy_xxy_0, tez_xxyy_xxz_0, tez_xxz_xx_0, tez_xxz_xx_1, \
                                         tez_xxz_xxx_0, tez_xxz_xxx_1, tez_xxz_xxy_0, tez_xxz_xxy_1, tez_xxz_xxz_0, \
                                         tez_xxz_xxz_1, tez_xxz_xy_0, tez_xxz_xy_1, tez_xxz_xyy_0, tez_xxz_xyy_1, \
                                         tez_xxz_xyz_0, tez_xxz_xyz_1, tez_xxz_xz_0, tez_xxz_xz_1, tez_xxz_xzz_0, \
                                         tez_xxz_xzz_1, tez_xxz_yy_0, tez_xxz_yy_1, tez_xxz_yyy_0, tez_xxz_yyy_1, \
                                         tez_xxz_yyz_0, tez_xxz_yyz_1, tez_xxz_yz_0, tez_xxz_yz_1, tez_xxz_yzz_0, \
                                         tez_xxz_yzz_1, tez_xxz_zz_0, tez_xxz_zz_1, tez_xxz_zzz_0, tez_xxz_zzz_1, \
                                         tez_xy_yyy_0, tez_xy_yyy_1, tez_xy_yyz_0, tez_xy_yyz_1, tez_xy_yzz_0, tez_xy_yzz_1, \
                                         tez_xy_zzz_0, tez_xy_zzz_1, tez_xyy_xx_0, tez_xyy_xx_1, tez_xyy_xxx_0, \
                                         tez_xyy_xxx_1, tez_xyy_xxy_0, tez_xyy_xxy_1, tez_xyy_xxz_0, tez_xyy_xxz_1, \
                                         tez_xyy_xy_0, tez_xyy_xy_1, tez_xyy_xz_0, tez_xyy_xz_1, tez_xz_xxx_0, tez_xz_xxx_1, \
                                         tez_xz_xxy_0, tez_xz_xxy_1, tez_xz_xxz_0, tez_xz_xxz_1, tez_xz_xyy_0, tez_xz_xyy_1, \
                                         tez_xz_xyz_0, tez_xz_xyz_1, tez_xz_xzz_0, tez_xz_xzz_1, tez_xz_yyy_0, tez_xz_yyy_1, \
                                         tez_xz_yyz_0, tez_xz_yyz_1, tez_xz_yzz_0, tez_xz_yzz_1, tez_xz_zzz_0, tez_xz_zzz_1, \
                                         tez_yy_xxx_0, tez_yy_xxx_1, tez_yy_xxy_0, tez_yy_xxy_1, tez_yy_xxz_0, tez_yy_xxz_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tez_xxxy_yyy_0[j] = pa_x[j] * tez_xxy_yyy_0[j] - pc_x[j] * tez_xxy_yyy_1[j] + fl1_fx * tez_xy_yyy_0[j] - fl1_fx * tez_xy_yyy_1[j];

                tex_xxxy_yyz_0[j] =
                    pa_x[j] * tex_xxy_yyz_0[j] - pc_x[j] * tex_xxy_yyz_1[j] + fl1_fx * tex_xy_yyz_0[j] - fl1_fx * tex_xy_yyz_1[j] + ta_xxy_yyz_1[j];

                tey_xxxy_yyz_0[j] = pa_x[j] * tey_xxy_yyz_0[j] - pc_x[j] * tey_xxy_yyz_1[j] + fl1_fx * tey_xy_yyz_0[j] - fl1_fx * tey_xy_yyz_1[j];

                tez_xxxy_yyz_0[j] = pa_x[j] * tez_xxy_yyz_0[j] - pc_x[j] * tez_xxy_yyz_1[j] + fl1_fx * tez_xy_yyz_0[j] - fl1_fx * tez_xy_yyz_1[j];

                tex_xxxy_yzz_0[j] =
                    pa_x[j] * tex_xxy_yzz_0[j] - pc_x[j] * tex_xxy_yzz_1[j] + fl1_fx * tex_xy_yzz_0[j] - fl1_fx * tex_xy_yzz_1[j] + ta_xxy_yzz_1[j];

                tey_xxxy_yzz_0[j] = pa_x[j] * tey_xxy_yzz_0[j] - pc_x[j] * tey_xxy_yzz_1[j] + fl1_fx * tey_xy_yzz_0[j] - fl1_fx * tey_xy_yzz_1[j];

                tez_xxxy_yzz_0[j] = pa_x[j] * tez_xxy_yzz_0[j] - pc_x[j] * tez_xxy_yzz_1[j] + fl1_fx * tez_xy_yzz_0[j] - fl1_fx * tez_xy_yzz_1[j];

                tex_xxxy_zzz_0[j] =
                    pa_x[j] * tex_xxy_zzz_0[j] - pc_x[j] * tex_xxy_zzz_1[j] + fl1_fx * tex_xy_zzz_0[j] - fl1_fx * tex_xy_zzz_1[j] + ta_xxy_zzz_1[j];

                tey_xxxy_zzz_0[j] = pa_x[j] * tey_xxy_zzz_0[j] - pc_x[j] * tey_xxy_zzz_1[j] + fl1_fx * tey_xy_zzz_0[j] - fl1_fx * tey_xy_zzz_1[j];

                tez_xxxy_zzz_0[j] = pa_x[j] * tez_xxy_zzz_0[j] - pc_x[j] * tez_xxy_zzz_1[j] + fl1_fx * tez_xy_zzz_0[j] - fl1_fx * tez_xy_zzz_1[j];

                tex_xxxz_xxx_0[j] = pa_x[j] * tex_xxz_xxx_0[j] - pc_x[j] * tex_xxz_xxx_1[j] + fl1_fx * tex_xz_xxx_0[j] - fl1_fx * tex_xz_xxx_1[j] +
                                    1.5 * fl1_fx * tex_xxz_xx_0[j] - 1.5 * fl1_fx * tex_xxz_xx_1[j] + ta_xxz_xxx_1[j];

                tey_xxxz_xxx_0[j] = pa_x[j] * tey_xxz_xxx_0[j] - pc_x[j] * tey_xxz_xxx_1[j] + fl1_fx * tey_xz_xxx_0[j] - fl1_fx * tey_xz_xxx_1[j] +
                                    1.5 * fl1_fx * tey_xxz_xx_0[j] - 1.5 * fl1_fx * tey_xxz_xx_1[j];

                tez_xxxz_xxx_0[j] = pa_x[j] * tez_xxz_xxx_0[j] - pc_x[j] * tez_xxz_xxx_1[j] + fl1_fx * tez_xz_xxx_0[j] - fl1_fx * tez_xz_xxx_1[j] +
                                    1.5 * fl1_fx * tez_xxz_xx_0[j] - 1.5 * fl1_fx * tez_xxz_xx_1[j];

                tex_xxxz_xxy_0[j] = pa_x[j] * tex_xxz_xxy_0[j] - pc_x[j] * tex_xxz_xxy_1[j] + fl1_fx * tex_xz_xxy_0[j] - fl1_fx * tex_xz_xxy_1[j] +
                                    fl1_fx * tex_xxz_xy_0[j] - fl1_fx * tex_xxz_xy_1[j] + ta_xxz_xxy_1[j];

                tey_xxxz_xxy_0[j] = pa_x[j] * tey_xxz_xxy_0[j] - pc_x[j] * tey_xxz_xxy_1[j] + fl1_fx * tey_xz_xxy_0[j] - fl1_fx * tey_xz_xxy_1[j] +
                                    fl1_fx * tey_xxz_xy_0[j] - fl1_fx * tey_xxz_xy_1[j];

                tez_xxxz_xxy_0[j] = pa_x[j] * tez_xxz_xxy_0[j] - pc_x[j] * tez_xxz_xxy_1[j] + fl1_fx * tez_xz_xxy_0[j] - fl1_fx * tez_xz_xxy_1[j] +
                                    fl1_fx * tez_xxz_xy_0[j] - fl1_fx * tez_xxz_xy_1[j];

                tex_xxxz_xxz_0[j] = pa_x[j] * tex_xxz_xxz_0[j] - pc_x[j] * tex_xxz_xxz_1[j] + fl1_fx * tex_xz_xxz_0[j] - fl1_fx * tex_xz_xxz_1[j] +
                                    fl1_fx * tex_xxz_xz_0[j] - fl1_fx * tex_xxz_xz_1[j] + ta_xxz_xxz_1[j];

                tey_xxxz_xxz_0[j] = pa_x[j] * tey_xxz_xxz_0[j] - pc_x[j] * tey_xxz_xxz_1[j] + fl1_fx * tey_xz_xxz_0[j] - fl1_fx * tey_xz_xxz_1[j] +
                                    fl1_fx * tey_xxz_xz_0[j] - fl1_fx * tey_xxz_xz_1[j];

                tez_xxxz_xxz_0[j] = pa_x[j] * tez_xxz_xxz_0[j] - pc_x[j] * tez_xxz_xxz_1[j] + fl1_fx * tez_xz_xxz_0[j] - fl1_fx * tez_xz_xxz_1[j] +
                                    fl1_fx * tez_xxz_xz_0[j] - fl1_fx * tez_xxz_xz_1[j];

                tex_xxxz_xyy_0[j] = pa_x[j] * tex_xxz_xyy_0[j] - pc_x[j] * tex_xxz_xyy_1[j] + fl1_fx * tex_xz_xyy_0[j] - fl1_fx * tex_xz_xyy_1[j] +
                                    0.5 * fl1_fx * tex_xxz_yy_0[j] - 0.5 * fl1_fx * tex_xxz_yy_1[j] + ta_xxz_xyy_1[j];

                tey_xxxz_xyy_0[j] = pa_x[j] * tey_xxz_xyy_0[j] - pc_x[j] * tey_xxz_xyy_1[j] + fl1_fx * tey_xz_xyy_0[j] - fl1_fx * tey_xz_xyy_1[j] +
                                    0.5 * fl1_fx * tey_xxz_yy_0[j] - 0.5 * fl1_fx * tey_xxz_yy_1[j];

                tez_xxxz_xyy_0[j] = pa_x[j] * tez_xxz_xyy_0[j] - pc_x[j] * tez_xxz_xyy_1[j] + fl1_fx * tez_xz_xyy_0[j] - fl1_fx * tez_xz_xyy_1[j] +
                                    0.5 * fl1_fx * tez_xxz_yy_0[j] - 0.5 * fl1_fx * tez_xxz_yy_1[j];

                tex_xxxz_xyz_0[j] = pa_x[j] * tex_xxz_xyz_0[j] - pc_x[j] * tex_xxz_xyz_1[j] + fl1_fx * tex_xz_xyz_0[j] - fl1_fx * tex_xz_xyz_1[j] +
                                    0.5 * fl1_fx * tex_xxz_yz_0[j] - 0.5 * fl1_fx * tex_xxz_yz_1[j] + ta_xxz_xyz_1[j];

                tey_xxxz_xyz_0[j] = pa_x[j] * tey_xxz_xyz_0[j] - pc_x[j] * tey_xxz_xyz_1[j] + fl1_fx * tey_xz_xyz_0[j] - fl1_fx * tey_xz_xyz_1[j] +
                                    0.5 * fl1_fx * tey_xxz_yz_0[j] - 0.5 * fl1_fx * tey_xxz_yz_1[j];

                tez_xxxz_xyz_0[j] = pa_x[j] * tez_xxz_xyz_0[j] - pc_x[j] * tez_xxz_xyz_1[j] + fl1_fx * tez_xz_xyz_0[j] - fl1_fx * tez_xz_xyz_1[j] +
                                    0.5 * fl1_fx * tez_xxz_yz_0[j] - 0.5 * fl1_fx * tez_xxz_yz_1[j];

                tex_xxxz_xzz_0[j] = pa_x[j] * tex_xxz_xzz_0[j] - pc_x[j] * tex_xxz_xzz_1[j] + fl1_fx * tex_xz_xzz_0[j] - fl1_fx * tex_xz_xzz_1[j] +
                                    0.5 * fl1_fx * tex_xxz_zz_0[j] - 0.5 * fl1_fx * tex_xxz_zz_1[j] + ta_xxz_xzz_1[j];

                tey_xxxz_xzz_0[j] = pa_x[j] * tey_xxz_xzz_0[j] - pc_x[j] * tey_xxz_xzz_1[j] + fl1_fx * tey_xz_xzz_0[j] - fl1_fx * tey_xz_xzz_1[j] +
                                    0.5 * fl1_fx * tey_xxz_zz_0[j] - 0.5 * fl1_fx * tey_xxz_zz_1[j];

                tez_xxxz_xzz_0[j] = pa_x[j] * tez_xxz_xzz_0[j] - pc_x[j] * tez_xxz_xzz_1[j] + fl1_fx * tez_xz_xzz_0[j] - fl1_fx * tez_xz_xzz_1[j] +
                                    0.5 * fl1_fx * tez_xxz_zz_0[j] - 0.5 * fl1_fx * tez_xxz_zz_1[j];

                tex_xxxz_yyy_0[j] =
                    pa_x[j] * tex_xxz_yyy_0[j] - pc_x[j] * tex_xxz_yyy_1[j] + fl1_fx * tex_xz_yyy_0[j] - fl1_fx * tex_xz_yyy_1[j] + ta_xxz_yyy_1[j];

                tey_xxxz_yyy_0[j] = pa_x[j] * tey_xxz_yyy_0[j] - pc_x[j] * tey_xxz_yyy_1[j] + fl1_fx * tey_xz_yyy_0[j] - fl1_fx * tey_xz_yyy_1[j];

                tez_xxxz_yyy_0[j] = pa_x[j] * tez_xxz_yyy_0[j] - pc_x[j] * tez_xxz_yyy_1[j] + fl1_fx * tez_xz_yyy_0[j] - fl1_fx * tez_xz_yyy_1[j];

                tex_xxxz_yyz_0[j] =
                    pa_x[j] * tex_xxz_yyz_0[j] - pc_x[j] * tex_xxz_yyz_1[j] + fl1_fx * tex_xz_yyz_0[j] - fl1_fx * tex_xz_yyz_1[j] + ta_xxz_yyz_1[j];

                tey_xxxz_yyz_0[j] = pa_x[j] * tey_xxz_yyz_0[j] - pc_x[j] * tey_xxz_yyz_1[j] + fl1_fx * tey_xz_yyz_0[j] - fl1_fx * tey_xz_yyz_1[j];

                tez_xxxz_yyz_0[j] = pa_x[j] * tez_xxz_yyz_0[j] - pc_x[j] * tez_xxz_yyz_1[j] + fl1_fx * tez_xz_yyz_0[j] - fl1_fx * tez_xz_yyz_1[j];

                tex_xxxz_yzz_0[j] =
                    pa_x[j] * tex_xxz_yzz_0[j] - pc_x[j] * tex_xxz_yzz_1[j] + fl1_fx * tex_xz_yzz_0[j] - fl1_fx * tex_xz_yzz_1[j] + ta_xxz_yzz_1[j];

                tey_xxxz_yzz_0[j] = pa_x[j] * tey_xxz_yzz_0[j] - pc_x[j] * tey_xxz_yzz_1[j] + fl1_fx * tey_xz_yzz_0[j] - fl1_fx * tey_xz_yzz_1[j];

                tez_xxxz_yzz_0[j] = pa_x[j] * tez_xxz_yzz_0[j] - pc_x[j] * tez_xxz_yzz_1[j] + fl1_fx * tez_xz_yzz_0[j] - fl1_fx * tez_xz_yzz_1[j];

                tex_xxxz_zzz_0[j] =
                    pa_x[j] * tex_xxz_zzz_0[j] - pc_x[j] * tex_xxz_zzz_1[j] + fl1_fx * tex_xz_zzz_0[j] - fl1_fx * tex_xz_zzz_1[j] + ta_xxz_zzz_1[j];

                tey_xxxz_zzz_0[j] = pa_x[j] * tey_xxz_zzz_0[j] - pc_x[j] * tey_xxz_zzz_1[j] + fl1_fx * tey_xz_zzz_0[j] - fl1_fx * tey_xz_zzz_1[j];

                tez_xxxz_zzz_0[j] = pa_x[j] * tez_xxz_zzz_0[j] - pc_x[j] * tez_xxz_zzz_1[j] + fl1_fx * tez_xz_zzz_0[j] - fl1_fx * tez_xz_zzz_1[j];

                tex_xxyy_xxx_0[j] = pa_x[j] * tex_xyy_xxx_0[j] - pc_x[j] * tex_xyy_xxx_1[j] + 0.5 * fl1_fx * tex_yy_xxx_0[j] -
                                    0.5 * fl1_fx * tex_yy_xxx_1[j] + 1.5 * fl1_fx * tex_xyy_xx_0[j] - 1.5 * fl1_fx * tex_xyy_xx_1[j] +
                                    ta_xyy_xxx_1[j];

                tey_xxyy_xxx_0[j] = pa_x[j] * tey_xyy_xxx_0[j] - pc_x[j] * tey_xyy_xxx_1[j] + 0.5 * fl1_fx * tey_yy_xxx_0[j] -
                                    0.5 * fl1_fx * tey_yy_xxx_1[j] + 1.5 * fl1_fx * tey_xyy_xx_0[j] - 1.5 * fl1_fx * tey_xyy_xx_1[j];

                tez_xxyy_xxx_0[j] = pa_x[j] * tez_xyy_xxx_0[j] - pc_x[j] * tez_xyy_xxx_1[j] + 0.5 * fl1_fx * tez_yy_xxx_0[j] -
                                    0.5 * fl1_fx * tez_yy_xxx_1[j] + 1.5 * fl1_fx * tez_xyy_xx_0[j] - 1.5 * fl1_fx * tez_xyy_xx_1[j];

                tex_xxyy_xxy_0[j] = pa_x[j] * tex_xyy_xxy_0[j] - pc_x[j] * tex_xyy_xxy_1[j] + 0.5 * fl1_fx * tex_yy_xxy_0[j] -
                                    0.5 * fl1_fx * tex_yy_xxy_1[j] + fl1_fx * tex_xyy_xy_0[j] - fl1_fx * tex_xyy_xy_1[j] + ta_xyy_xxy_1[j];

                tey_xxyy_xxy_0[j] = pa_x[j] * tey_xyy_xxy_0[j] - pc_x[j] * tey_xyy_xxy_1[j] + 0.5 * fl1_fx * tey_yy_xxy_0[j] -
                                    0.5 * fl1_fx * tey_yy_xxy_1[j] + fl1_fx * tey_xyy_xy_0[j] - fl1_fx * tey_xyy_xy_1[j];

                tez_xxyy_xxy_0[j] = pa_x[j] * tez_xyy_xxy_0[j] - pc_x[j] * tez_xyy_xxy_1[j] + 0.5 * fl1_fx * tez_yy_xxy_0[j] -
                                    0.5 * fl1_fx * tez_yy_xxy_1[j] + fl1_fx * tez_xyy_xy_0[j] - fl1_fx * tez_xyy_xy_1[j];

                tex_xxyy_xxz_0[j] = pa_x[j] * tex_xyy_xxz_0[j] - pc_x[j] * tex_xyy_xxz_1[j] + 0.5 * fl1_fx * tex_yy_xxz_0[j] -
                                    0.5 * fl1_fx * tex_yy_xxz_1[j] + fl1_fx * tex_xyy_xz_0[j] - fl1_fx * tex_xyy_xz_1[j] + ta_xyy_xxz_1[j];

                tey_xxyy_xxz_0[j] = pa_x[j] * tey_xyy_xxz_0[j] - pc_x[j] * tey_xyy_xxz_1[j] + 0.5 * fl1_fx * tey_yy_xxz_0[j] -
                                    0.5 * fl1_fx * tey_yy_xxz_1[j] + fl1_fx * tey_xyy_xz_0[j] - fl1_fx * tey_xyy_xz_1[j];

                tez_xxyy_xxz_0[j] = pa_x[j] * tez_xyy_xxz_0[j] - pc_x[j] * tez_xyy_xxz_1[j] + 0.5 * fl1_fx * tez_yy_xxz_0[j] -
                                    0.5 * fl1_fx * tez_yy_xxz_1[j] + fl1_fx * tez_xyy_xz_0[j] - fl1_fx * tez_xyy_xz_1[j];

                tex_xxyy_xyy_0[j] = pa_x[j] * tex_xyy_xyy_0[j] - pc_x[j] * tex_xyy_xyy_1[j] + 0.5 * fl1_fx * tex_yy_xyy_0[j] -
                                    0.5 * fl1_fx * tex_yy_xyy_1[j] + 0.5 * fl1_fx * tex_xyy_yy_0[j] - 0.5 * fl1_fx * tex_xyy_yy_1[j] +
                                    ta_xyy_xyy_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForGF_100_150(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pcDistances,
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

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_4_3_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_2_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tey_xyy_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 33);

            auto tez_xyy_xyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 33);

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

            auto tey_xyy_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 33);

            auto tez_xyy_xyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 33);

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

            auto tey_yy_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 33);

            auto tez_yy_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 33);

            auto tex_yy_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 34);

            auto tey_yy_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 34);

            auto tez_yy_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 34);

            auto tex_yy_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 35);

            auto tey_yy_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 35);

            auto tez_yy_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 35);

            auto tex_yy_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 36);

            auto tey_yy_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 36);

            auto tez_yy_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 36);

            auto tex_yy_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 37);

            auto tey_yy_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 37);

            auto tez_yy_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 37);

            auto tex_yy_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 38);

            auto tey_yy_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 38);

            auto tez_yy_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 38);

            auto tex_yy_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 39);

            auto tey_yy_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 39);

            auto tez_yy_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 39);

            auto tex_yz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 40);

            auto tey_yz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 40);

            auto tez_yz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 40);

            auto tex_yz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 41);

            auto tey_yz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 41);

            auto tez_yz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 41);

            auto tex_yz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 42);

            auto tey_yz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 42);

            auto tez_yz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 42);

            auto tex_yz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 43);

            auto tey_yz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 43);

            auto tez_yz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 43);

            auto tex_yz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 44);

            auto tey_yz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 44);

            auto tez_yz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 44);

            auto tex_yz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 45);

            auto tey_yz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 45);

            auto tez_yz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 45);

            auto tex_yz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 46);

            auto tey_yz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 46);

            auto tez_yz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 46);

            auto tex_yz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 47);

            auto tey_yz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 47);

            auto tez_yz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 47);

            auto tex_yz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 48);

            auto tey_yz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 48);

            auto tez_yz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 48);

            auto tex_yz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 49);

            auto tey_yz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 49);

            auto tez_yz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 49);

            auto tey_yy_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 33);

            auto tez_yy_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 33);

            auto tex_yy_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 34);

            auto tey_yy_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 34);

            auto tez_yy_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 34);

            auto tex_yy_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 35);

            auto tey_yy_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 35);

            auto tez_yy_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 35);

            auto tex_yy_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 36);

            auto tey_yy_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 36);

            auto tez_yy_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 36);

            auto tex_yy_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 37);

            auto tey_yy_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 37);

            auto tez_yy_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 37);

            auto tex_yy_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 38);

            auto tey_yy_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 38);

            auto tez_yy_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 38);

            auto tex_yy_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 39);

            auto tey_yy_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 39);

            auto tez_yy_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 39);

            auto tex_yz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 40);

            auto tey_yz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 40);

            auto tez_yz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 40);

            auto tex_yz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 41);

            auto tey_yz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 41);

            auto tez_yz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 41);

            auto tex_yz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 42);

            auto tey_yz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 42);

            auto tez_yz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 42);

            auto tex_yz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 43);

            auto tey_yz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 43);

            auto tez_yz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 43);

            auto tex_yz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 44);

            auto tey_yz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 44);

            auto tez_yz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 44);

            auto tex_yz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 45);

            auto tey_yz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 45);

            auto tez_yz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 45);

            auto tex_yz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 46);

            auto tey_yz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 46);

            auto tez_yz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 46);

            auto tex_yz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 47);

            auto tey_yz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 47);

            auto tez_yz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 47);

            auto tex_yz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 48);

            auto tey_yz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 48);

            auto tez_yz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 48);

            auto tex_yz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 49);

            auto tey_yz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 49);

            auto tez_yz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 49);

            auto tey_xyy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 21);

            auto tez_xyy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 21);

            auto tex_xyy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 22);

            auto tey_xyy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 22);

            auto tez_xyy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 22);

            auto tex_xyy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 23);

            auto tey_xyy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 23);

            auto tez_xyy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 23);

            auto tex_xyz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 24);

            auto tey_xyz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 24);

            auto tez_xyz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 24);

            auto tex_xyz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 25);

            auto tey_xyz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 25);

            auto tez_xyz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 25);

            auto tex_xyz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 26);

            auto tey_xyz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 26);

            auto tez_xyz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 26);

            auto tex_xyz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 27);

            auto tey_xyz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 27);

            auto tez_xyz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 27);

            auto tex_xyz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 28);

            auto tey_xyz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 28);

            auto tez_xyz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 28);

            auto tex_xyz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 29);

            auto tey_xyz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 29);

            auto tez_xyz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 29);

            auto tey_xyy_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 21);

            auto tez_xyy_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 21);

            auto tex_xyy_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 22);

            auto tey_xyy_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 22);

            auto tez_xyy_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 22);

            auto tex_xyy_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 23);

            auto tey_xyy_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 23);

            auto tez_xyy_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 23);

            auto tex_xyz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 24);

            auto tey_xyz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 24);

            auto tez_xyz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 24);

            auto tex_xyz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 25);

            auto tey_xyz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 25);

            auto tez_xyz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 25);

            auto tex_xyz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 26);

            auto tey_xyz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 26);

            auto tez_xyz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 26);

            auto tex_xyz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 27);

            auto tey_xyz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 27);

            auto tez_xyz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 27);

            auto tex_xyz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 28);

            auto tey_xyz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 28);

            auto tez_xyz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 28);

            auto tex_xyz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 29);

            auto tey_xyz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 29);

            auto tez_xyz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 29);

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

            // set up pointers to integrals

            auto tey_xxyy_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 33);

            auto tez_xxyy_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 33);

            auto tex_xxyy_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 34);

            auto tey_xxyy_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 34);

            auto tez_xxyy_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 34);

            auto tex_xxyy_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 35);

            auto tey_xxyy_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 35);

            auto tez_xxyy_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 35);

            auto tex_xxyy_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 36);

            auto tey_xxyy_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 36);

            auto tez_xxyy_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 36);

            auto tex_xxyy_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 37);

            auto tey_xxyy_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 37);

            auto tez_xxyy_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 37);

            auto tex_xxyy_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 38);

            auto tey_xxyy_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 38);

            auto tez_xxyy_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 38);

            auto tex_xxyy_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 39);

            auto tey_xxyy_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 39);

            auto tez_xxyy_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 39);

            auto tex_xxyz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 40);

            auto tey_xxyz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 40);

            auto tez_xxyz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 40);

            auto tex_xxyz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 41);

            auto tey_xxyz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 41);

            auto tez_xxyz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 41);

            auto tex_xxyz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 42);

            auto tey_xxyz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 42);

            auto tez_xxyz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 42);

            auto tex_xxyz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 43);

            auto tey_xxyz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 43);

            auto tez_xxyz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 43);

            auto tex_xxyz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 44);

            auto tey_xxyz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 44);

            auto tez_xxyz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 44);

            auto tex_xxyz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 45);

            auto tey_xxyz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 45);

            auto tez_xxyz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 45);

            auto tex_xxyz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 46);

            auto tey_xxyz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 46);

            auto tez_xxyz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 46);

            auto tex_xxyz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 47);

            auto tey_xxyz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 47);

            auto tez_xxyz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 47);

            auto tex_xxyz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 48);

            auto tey_xxyz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 48);

            auto tez_xxyz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 48);

            auto tex_xxyz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 49);

            auto tey_xxyz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 49);

            auto tez_xxyz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 49);

            // Batch of Integrals (100,150)

            #pragma omp simd aligned(fx, pa_x, pc_x, ta_xyy_xyz_1, ta_xyy_xzz_1, ta_xyy_yyy_1, ta_xyy_yyz_1, \
                                         ta_xyy_yzz_1, ta_xyy_zzz_1, ta_xyz_xxx_1, ta_xyz_xxy_1, ta_xyz_xxz_1, ta_xyz_xyy_1, \
                                         ta_xyz_xyz_1, ta_xyz_xzz_1, ta_xyz_yyy_1, ta_xyz_yyz_1, ta_xyz_yzz_1, ta_xyz_zzz_1, \
                                         tex_xxyy_xyz_0, tex_xxyy_xzz_0, tex_xxyy_yyy_0, tex_xxyy_yyz_0, tex_xxyy_yzz_0, \
                                         tex_xxyy_zzz_0, tex_xxyz_xxx_0, tex_xxyz_xxy_0, tex_xxyz_xxz_0, tex_xxyz_xyy_0, \
                                         tex_xxyz_xyz_0, tex_xxyz_xzz_0, tex_xxyz_yyy_0, tex_xxyz_yyz_0, tex_xxyz_yzz_0, \
                                         tex_xxyz_zzz_0, tex_xyy_xyz_0, tex_xyy_xyz_1, tex_xyy_xzz_0, tex_xyy_xzz_1, \
                                         tex_xyy_yyy_0, tex_xyy_yyy_1, tex_xyy_yyz_0, tex_xyy_yyz_1, tex_xyy_yz_0, \
                                         tex_xyy_yz_1, tex_xyy_yzz_0, tex_xyy_yzz_1, tex_xyy_zz_0, tex_xyy_zz_1, \
                                         tex_xyy_zzz_0, tex_xyy_zzz_1, tex_xyz_xx_0, tex_xyz_xx_1, tex_xyz_xxx_0, \
                                         tex_xyz_xxx_1, tex_xyz_xxy_0, tex_xyz_xxy_1, tex_xyz_xxz_0, tex_xyz_xxz_1, \
                                         tex_xyz_xy_0, tex_xyz_xy_1, tex_xyz_xyy_0, tex_xyz_xyy_1, tex_xyz_xyz_0, \
                                         tex_xyz_xyz_1, tex_xyz_xz_0, tex_xyz_xz_1, tex_xyz_xzz_0, tex_xyz_xzz_1, \
                                         tex_xyz_yy_0, tex_xyz_yy_1, tex_xyz_yyy_0, tex_xyz_yyy_1, tex_xyz_yyz_0, \
                                         tex_xyz_yyz_1, tex_xyz_yz_0, tex_xyz_yz_1, tex_xyz_yzz_0, tex_xyz_yzz_1, \
                                         tex_xyz_zz_0, tex_xyz_zz_1, tex_xyz_zzz_0, tex_xyz_zzz_1, tex_yy_xyz_0, \
                                         tex_yy_xyz_1, tex_yy_xzz_0, tex_yy_xzz_1, tex_yy_yyy_0, tex_yy_yyy_1, tex_yy_yyz_0, \
                                         tex_yy_yyz_1, tex_yy_yzz_0, tex_yy_yzz_1, tex_yy_zzz_0, tex_yy_zzz_1, tex_yz_xxx_0, \
                                         tex_yz_xxx_1, tex_yz_xxy_0, tex_yz_xxy_1, tex_yz_xxz_0, tex_yz_xxz_1, tex_yz_xyy_0, \
                                         tex_yz_xyy_1, tex_yz_xyz_0, tex_yz_xyz_1, tex_yz_xzz_0, tex_yz_xzz_1, tex_yz_yyy_0, \
                                         tex_yz_yyy_1, tex_yz_yyz_0, tex_yz_yyz_1, tex_yz_yzz_0, tex_yz_yzz_1, tex_yz_zzz_0, \
                                         tex_yz_zzz_1, tey_xxyy_xyy_0, tey_xxyy_xyz_0, tey_xxyy_xzz_0, tey_xxyy_yyy_0, \
                                         tey_xxyy_yyz_0, tey_xxyy_yzz_0, tey_xxyy_zzz_0, tey_xxyz_xxx_0, tey_xxyz_xxy_0, \
                                         tey_xxyz_xxz_0, tey_xxyz_xyy_0, tey_xxyz_xyz_0, tey_xxyz_xzz_0, tey_xxyz_yyy_0, \
                                         tey_xxyz_yyz_0, tey_xxyz_yzz_0, tey_xxyz_zzz_0, tey_xyy_xyy_0, tey_xyy_xyy_1, \
                                         tey_xyy_xyz_0, tey_xyy_xyz_1, tey_xyy_xzz_0, tey_xyy_xzz_1, tey_xyy_yy_0, \
                                         tey_xyy_yy_1, tey_xyy_yyy_0, tey_xyy_yyy_1, tey_xyy_yyz_0, tey_xyy_yyz_1, \
                                         tey_xyy_yz_0, tey_xyy_yz_1, tey_xyy_yzz_0, tey_xyy_yzz_1, tey_xyy_zz_0, \
                                         tey_xyy_zz_1, tey_xyy_zzz_0, tey_xyy_zzz_1, tey_xyz_xx_0, tey_xyz_xx_1, \
                                         tey_xyz_xxx_0, tey_xyz_xxx_1, tey_xyz_xxy_0, tey_xyz_xxy_1, tey_xyz_xxz_0, \
                                         tey_xyz_xxz_1, tey_xyz_xy_0, tey_xyz_xy_1, tey_xyz_xyy_0, tey_xyz_xyy_1, \
                                         tey_xyz_xyz_0, tey_xyz_xyz_1, tey_xyz_xz_0, tey_xyz_xz_1, tey_xyz_xzz_0, \
                                         tey_xyz_xzz_1, tey_xyz_yy_0, tey_xyz_yy_1, tey_xyz_yyy_0, tey_xyz_yyy_1, \
                                         tey_xyz_yyz_0, tey_xyz_yyz_1, tey_xyz_yz_0, tey_xyz_yz_1, tey_xyz_yzz_0, \
                                         tey_xyz_yzz_1, tey_xyz_zz_0, tey_xyz_zz_1, tey_xyz_zzz_0, tey_xyz_zzz_1, \
                                         tey_yy_xyy_0, tey_yy_xyy_1, tey_yy_xyz_0, tey_yy_xyz_1, tey_yy_xzz_0, tey_yy_xzz_1, \
                                         tey_yy_yyy_0, tey_yy_yyy_1, tey_yy_yyz_0, tey_yy_yyz_1, tey_yy_yzz_0, tey_yy_yzz_1, \
                                         tey_yy_zzz_0, tey_yy_zzz_1, tey_yz_xxx_0, tey_yz_xxx_1, tey_yz_xxy_0, tey_yz_xxy_1, \
                                         tey_yz_xxz_0, tey_yz_xxz_1, tey_yz_xyy_0, tey_yz_xyy_1, tey_yz_xyz_0, tey_yz_xyz_1, \
                                         tey_yz_xzz_0, tey_yz_xzz_1, tey_yz_yyy_0, tey_yz_yyy_1, tey_yz_yyz_0, tey_yz_yyz_1, \
                                         tey_yz_yzz_0, tey_yz_yzz_1, tey_yz_zzz_0, tey_yz_zzz_1, tez_xxyy_xyy_0, \
                                         tez_xxyy_xyz_0, tez_xxyy_xzz_0, tez_xxyy_yyy_0, tez_xxyy_yyz_0, tez_xxyy_yzz_0, \
                                         tez_xxyy_zzz_0, tez_xxyz_xxx_0, tez_xxyz_xxy_0, tez_xxyz_xxz_0, tez_xxyz_xyy_0, \
                                         tez_xxyz_xyz_0, tez_xxyz_xzz_0, tez_xxyz_yyy_0, tez_xxyz_yyz_0, tez_xxyz_yzz_0, \
                                         tez_xxyz_zzz_0, tez_xyy_xyy_0, tez_xyy_xyy_1, tez_xyy_xyz_0, tez_xyy_xyz_1, \
                                         tez_xyy_xzz_0, tez_xyy_xzz_1, tez_xyy_yy_0, tez_xyy_yy_1, tez_xyy_yyy_0, \
                                         tez_xyy_yyy_1, tez_xyy_yyz_0, tez_xyy_yyz_1, tez_xyy_yz_0, tez_xyy_yz_1, \
                                         tez_xyy_yzz_0, tez_xyy_yzz_1, tez_xyy_zz_0, tez_xyy_zz_1, tez_xyy_zzz_0, \
                                         tez_xyy_zzz_1, tez_xyz_xx_0, tez_xyz_xx_1, tez_xyz_xxx_0, tez_xyz_xxx_1, \
                                         tez_xyz_xxy_0, tez_xyz_xxy_1, tez_xyz_xxz_0, tez_xyz_xxz_1, tez_xyz_xy_0, \
                                         tez_xyz_xy_1, tez_xyz_xyy_0, tez_xyz_xyy_1, tez_xyz_xyz_0, tez_xyz_xyz_1, \
                                         tez_xyz_xz_0, tez_xyz_xz_1, tez_xyz_xzz_0, tez_xyz_xzz_1, tez_xyz_yy_0, \
                                         tez_xyz_yy_1, tez_xyz_yyy_0, tez_xyz_yyy_1, tez_xyz_yyz_0, tez_xyz_yyz_1, \
                                         tez_xyz_yz_0, tez_xyz_yz_1, tez_xyz_yzz_0, tez_xyz_yzz_1, tez_xyz_zz_0, \
                                         tez_xyz_zz_1, tez_xyz_zzz_0, tez_xyz_zzz_1, tez_yy_xyy_0, tez_yy_xyy_1, \
                                         tez_yy_xyz_0, tez_yy_xyz_1, tez_yy_xzz_0, tez_yy_xzz_1, tez_yy_yyy_0, tez_yy_yyy_1, \
                                         tez_yy_yyz_0, tez_yy_yyz_1, tez_yy_yzz_0, tez_yy_yzz_1, tez_yy_zzz_0, tez_yy_zzz_1, \
                                         tez_yz_xxx_0, tez_yz_xxx_1, tez_yz_xxy_0, tez_yz_xxy_1, tez_yz_xxz_0, tez_yz_xxz_1, \
                                         tez_yz_xyy_0, tez_yz_xyy_1, tez_yz_xyz_0, tez_yz_xyz_1, tez_yz_xzz_0, tez_yz_xzz_1, \
                                         tez_yz_yyy_0, tez_yz_yyy_1, tez_yz_yyz_0, tez_yz_yyz_1, tez_yz_yzz_0, tez_yz_yzz_1, \
                                         tez_yz_zzz_0, tez_yz_zzz_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tey_xxyy_xyy_0[j] = pa_x[j] * tey_xyy_xyy_0[j] - pc_x[j] * tey_xyy_xyy_1[j] + 0.5 * fl1_fx * tey_yy_xyy_0[j] -
                                    0.5 * fl1_fx * tey_yy_xyy_1[j] + 0.5 * fl1_fx * tey_xyy_yy_0[j] - 0.5 * fl1_fx * tey_xyy_yy_1[j];

                tez_xxyy_xyy_0[j] = pa_x[j] * tez_xyy_xyy_0[j] - pc_x[j] * tez_xyy_xyy_1[j] + 0.5 * fl1_fx * tez_yy_xyy_0[j] -
                                    0.5 * fl1_fx * tez_yy_xyy_1[j] + 0.5 * fl1_fx * tez_xyy_yy_0[j] - 0.5 * fl1_fx * tez_xyy_yy_1[j];

                tex_xxyy_xyz_0[j] = pa_x[j] * tex_xyy_xyz_0[j] - pc_x[j] * tex_xyy_xyz_1[j] + 0.5 * fl1_fx * tex_yy_xyz_0[j] -
                                    0.5 * fl1_fx * tex_yy_xyz_1[j] + 0.5 * fl1_fx * tex_xyy_yz_0[j] - 0.5 * fl1_fx * tex_xyy_yz_1[j] +
                                    ta_xyy_xyz_1[j];

                tey_xxyy_xyz_0[j] = pa_x[j] * tey_xyy_xyz_0[j] - pc_x[j] * tey_xyy_xyz_1[j] + 0.5 * fl1_fx * tey_yy_xyz_0[j] -
                                    0.5 * fl1_fx * tey_yy_xyz_1[j] + 0.5 * fl1_fx * tey_xyy_yz_0[j] - 0.5 * fl1_fx * tey_xyy_yz_1[j];

                tez_xxyy_xyz_0[j] = pa_x[j] * tez_xyy_xyz_0[j] - pc_x[j] * tez_xyy_xyz_1[j] + 0.5 * fl1_fx * tez_yy_xyz_0[j] -
                                    0.5 * fl1_fx * tez_yy_xyz_1[j] + 0.5 * fl1_fx * tez_xyy_yz_0[j] - 0.5 * fl1_fx * tez_xyy_yz_1[j];

                tex_xxyy_xzz_0[j] = pa_x[j] * tex_xyy_xzz_0[j] - pc_x[j] * tex_xyy_xzz_1[j] + 0.5 * fl1_fx * tex_yy_xzz_0[j] -
                                    0.5 * fl1_fx * tex_yy_xzz_1[j] + 0.5 * fl1_fx * tex_xyy_zz_0[j] - 0.5 * fl1_fx * tex_xyy_zz_1[j] +
                                    ta_xyy_xzz_1[j];

                tey_xxyy_xzz_0[j] = pa_x[j] * tey_xyy_xzz_0[j] - pc_x[j] * tey_xyy_xzz_1[j] + 0.5 * fl1_fx * tey_yy_xzz_0[j] -
                                    0.5 * fl1_fx * tey_yy_xzz_1[j] + 0.5 * fl1_fx * tey_xyy_zz_0[j] - 0.5 * fl1_fx * tey_xyy_zz_1[j];

                tez_xxyy_xzz_0[j] = pa_x[j] * tez_xyy_xzz_0[j] - pc_x[j] * tez_xyy_xzz_1[j] + 0.5 * fl1_fx * tez_yy_xzz_0[j] -
                                    0.5 * fl1_fx * tez_yy_xzz_1[j] + 0.5 * fl1_fx * tez_xyy_zz_0[j] - 0.5 * fl1_fx * tez_xyy_zz_1[j];

                tex_xxyy_yyy_0[j] = pa_x[j] * tex_xyy_yyy_0[j] - pc_x[j] * tex_xyy_yyy_1[j] + 0.5 * fl1_fx * tex_yy_yyy_0[j] -
                                    0.5 * fl1_fx * tex_yy_yyy_1[j] + ta_xyy_yyy_1[j];

                tey_xxyy_yyy_0[j] =
                    pa_x[j] * tey_xyy_yyy_0[j] - pc_x[j] * tey_xyy_yyy_1[j] + 0.5 * fl1_fx * tey_yy_yyy_0[j] - 0.5 * fl1_fx * tey_yy_yyy_1[j];

                tez_xxyy_yyy_0[j] =
                    pa_x[j] * tez_xyy_yyy_0[j] - pc_x[j] * tez_xyy_yyy_1[j] + 0.5 * fl1_fx * tez_yy_yyy_0[j] - 0.5 * fl1_fx * tez_yy_yyy_1[j];

                tex_xxyy_yyz_0[j] = pa_x[j] * tex_xyy_yyz_0[j] - pc_x[j] * tex_xyy_yyz_1[j] + 0.5 * fl1_fx * tex_yy_yyz_0[j] -
                                    0.5 * fl1_fx * tex_yy_yyz_1[j] + ta_xyy_yyz_1[j];

                tey_xxyy_yyz_0[j] =
                    pa_x[j] * tey_xyy_yyz_0[j] - pc_x[j] * tey_xyy_yyz_1[j] + 0.5 * fl1_fx * tey_yy_yyz_0[j] - 0.5 * fl1_fx * tey_yy_yyz_1[j];

                tez_xxyy_yyz_0[j] =
                    pa_x[j] * tez_xyy_yyz_0[j] - pc_x[j] * tez_xyy_yyz_1[j] + 0.5 * fl1_fx * tez_yy_yyz_0[j] - 0.5 * fl1_fx * tez_yy_yyz_1[j];

                tex_xxyy_yzz_0[j] = pa_x[j] * tex_xyy_yzz_0[j] - pc_x[j] * tex_xyy_yzz_1[j] + 0.5 * fl1_fx * tex_yy_yzz_0[j] -
                                    0.5 * fl1_fx * tex_yy_yzz_1[j] + ta_xyy_yzz_1[j];

                tey_xxyy_yzz_0[j] =
                    pa_x[j] * tey_xyy_yzz_0[j] - pc_x[j] * tey_xyy_yzz_1[j] + 0.5 * fl1_fx * tey_yy_yzz_0[j] - 0.5 * fl1_fx * tey_yy_yzz_1[j];

                tez_xxyy_yzz_0[j] =
                    pa_x[j] * tez_xyy_yzz_0[j] - pc_x[j] * tez_xyy_yzz_1[j] + 0.5 * fl1_fx * tez_yy_yzz_0[j] - 0.5 * fl1_fx * tez_yy_yzz_1[j];

                tex_xxyy_zzz_0[j] = pa_x[j] * tex_xyy_zzz_0[j] - pc_x[j] * tex_xyy_zzz_1[j] + 0.5 * fl1_fx * tex_yy_zzz_0[j] -
                                    0.5 * fl1_fx * tex_yy_zzz_1[j] + ta_xyy_zzz_1[j];

                tey_xxyy_zzz_0[j] =
                    pa_x[j] * tey_xyy_zzz_0[j] - pc_x[j] * tey_xyy_zzz_1[j] + 0.5 * fl1_fx * tey_yy_zzz_0[j] - 0.5 * fl1_fx * tey_yy_zzz_1[j];

                tez_xxyy_zzz_0[j] =
                    pa_x[j] * tez_xyy_zzz_0[j] - pc_x[j] * tez_xyy_zzz_1[j] + 0.5 * fl1_fx * tez_yy_zzz_0[j] - 0.5 * fl1_fx * tez_yy_zzz_1[j];

                tex_xxyz_xxx_0[j] = pa_x[j] * tex_xyz_xxx_0[j] - pc_x[j] * tex_xyz_xxx_1[j] + 0.5 * fl1_fx * tex_yz_xxx_0[j] -
                                    0.5 * fl1_fx * tex_yz_xxx_1[j] + 1.5 * fl1_fx * tex_xyz_xx_0[j] - 1.5 * fl1_fx * tex_xyz_xx_1[j] +
                                    ta_xyz_xxx_1[j];

                tey_xxyz_xxx_0[j] = pa_x[j] * tey_xyz_xxx_0[j] - pc_x[j] * tey_xyz_xxx_1[j] + 0.5 * fl1_fx * tey_yz_xxx_0[j] -
                                    0.5 * fl1_fx * tey_yz_xxx_1[j] + 1.5 * fl1_fx * tey_xyz_xx_0[j] - 1.5 * fl1_fx * tey_xyz_xx_1[j];

                tez_xxyz_xxx_0[j] = pa_x[j] * tez_xyz_xxx_0[j] - pc_x[j] * tez_xyz_xxx_1[j] + 0.5 * fl1_fx * tez_yz_xxx_0[j] -
                                    0.5 * fl1_fx * tez_yz_xxx_1[j] + 1.5 * fl1_fx * tez_xyz_xx_0[j] - 1.5 * fl1_fx * tez_xyz_xx_1[j];

                tex_xxyz_xxy_0[j] = pa_x[j] * tex_xyz_xxy_0[j] - pc_x[j] * tex_xyz_xxy_1[j] + 0.5 * fl1_fx * tex_yz_xxy_0[j] -
                                    0.5 * fl1_fx * tex_yz_xxy_1[j] + fl1_fx * tex_xyz_xy_0[j] - fl1_fx * tex_xyz_xy_1[j] + ta_xyz_xxy_1[j];

                tey_xxyz_xxy_0[j] = pa_x[j] * tey_xyz_xxy_0[j] - pc_x[j] * tey_xyz_xxy_1[j] + 0.5 * fl1_fx * tey_yz_xxy_0[j] -
                                    0.5 * fl1_fx * tey_yz_xxy_1[j] + fl1_fx * tey_xyz_xy_0[j] - fl1_fx * tey_xyz_xy_1[j];

                tez_xxyz_xxy_0[j] = pa_x[j] * tez_xyz_xxy_0[j] - pc_x[j] * tez_xyz_xxy_1[j] + 0.5 * fl1_fx * tez_yz_xxy_0[j] -
                                    0.5 * fl1_fx * tez_yz_xxy_1[j] + fl1_fx * tez_xyz_xy_0[j] - fl1_fx * tez_xyz_xy_1[j];

                tex_xxyz_xxz_0[j] = pa_x[j] * tex_xyz_xxz_0[j] - pc_x[j] * tex_xyz_xxz_1[j] + 0.5 * fl1_fx * tex_yz_xxz_0[j] -
                                    0.5 * fl1_fx * tex_yz_xxz_1[j] + fl1_fx * tex_xyz_xz_0[j] - fl1_fx * tex_xyz_xz_1[j] + ta_xyz_xxz_1[j];

                tey_xxyz_xxz_0[j] = pa_x[j] * tey_xyz_xxz_0[j] - pc_x[j] * tey_xyz_xxz_1[j] + 0.5 * fl1_fx * tey_yz_xxz_0[j] -
                                    0.5 * fl1_fx * tey_yz_xxz_1[j] + fl1_fx * tey_xyz_xz_0[j] - fl1_fx * tey_xyz_xz_1[j];

                tez_xxyz_xxz_0[j] = pa_x[j] * tez_xyz_xxz_0[j] - pc_x[j] * tez_xyz_xxz_1[j] + 0.5 * fl1_fx * tez_yz_xxz_0[j] -
                                    0.5 * fl1_fx * tez_yz_xxz_1[j] + fl1_fx * tez_xyz_xz_0[j] - fl1_fx * tez_xyz_xz_1[j];

                tex_xxyz_xyy_0[j] = pa_x[j] * tex_xyz_xyy_0[j] - pc_x[j] * tex_xyz_xyy_1[j] + 0.5 * fl1_fx * tex_yz_xyy_0[j] -
                                    0.5 * fl1_fx * tex_yz_xyy_1[j] + 0.5 * fl1_fx * tex_xyz_yy_0[j] - 0.5 * fl1_fx * tex_xyz_yy_1[j] +
                                    ta_xyz_xyy_1[j];

                tey_xxyz_xyy_0[j] = pa_x[j] * tey_xyz_xyy_0[j] - pc_x[j] * tey_xyz_xyy_1[j] + 0.5 * fl1_fx * tey_yz_xyy_0[j] -
                                    0.5 * fl1_fx * tey_yz_xyy_1[j] + 0.5 * fl1_fx * tey_xyz_yy_0[j] - 0.5 * fl1_fx * tey_xyz_yy_1[j];

                tez_xxyz_xyy_0[j] = pa_x[j] * tez_xyz_xyy_0[j] - pc_x[j] * tez_xyz_xyy_1[j] + 0.5 * fl1_fx * tez_yz_xyy_0[j] -
                                    0.5 * fl1_fx * tez_yz_xyy_1[j] + 0.5 * fl1_fx * tez_xyz_yy_0[j] - 0.5 * fl1_fx * tez_xyz_yy_1[j];

                tex_xxyz_xyz_0[j] = pa_x[j] * tex_xyz_xyz_0[j] - pc_x[j] * tex_xyz_xyz_1[j] + 0.5 * fl1_fx * tex_yz_xyz_0[j] -
                                    0.5 * fl1_fx * tex_yz_xyz_1[j] + 0.5 * fl1_fx * tex_xyz_yz_0[j] - 0.5 * fl1_fx * tex_xyz_yz_1[j] +
                                    ta_xyz_xyz_1[j];

                tey_xxyz_xyz_0[j] = pa_x[j] * tey_xyz_xyz_0[j] - pc_x[j] * tey_xyz_xyz_1[j] + 0.5 * fl1_fx * tey_yz_xyz_0[j] -
                                    0.5 * fl1_fx * tey_yz_xyz_1[j] + 0.5 * fl1_fx * tey_xyz_yz_0[j] - 0.5 * fl1_fx * tey_xyz_yz_1[j];

                tez_xxyz_xyz_0[j] = pa_x[j] * tez_xyz_xyz_0[j] - pc_x[j] * tez_xyz_xyz_1[j] + 0.5 * fl1_fx * tez_yz_xyz_0[j] -
                                    0.5 * fl1_fx * tez_yz_xyz_1[j] + 0.5 * fl1_fx * tez_xyz_yz_0[j] - 0.5 * fl1_fx * tez_xyz_yz_1[j];

                tex_xxyz_xzz_0[j] = pa_x[j] * tex_xyz_xzz_0[j] - pc_x[j] * tex_xyz_xzz_1[j] + 0.5 * fl1_fx * tex_yz_xzz_0[j] -
                                    0.5 * fl1_fx * tex_yz_xzz_1[j] + 0.5 * fl1_fx * tex_xyz_zz_0[j] - 0.5 * fl1_fx * tex_xyz_zz_1[j] +
                                    ta_xyz_xzz_1[j];

                tey_xxyz_xzz_0[j] = pa_x[j] * tey_xyz_xzz_0[j] - pc_x[j] * tey_xyz_xzz_1[j] + 0.5 * fl1_fx * tey_yz_xzz_0[j] -
                                    0.5 * fl1_fx * tey_yz_xzz_1[j] + 0.5 * fl1_fx * tey_xyz_zz_0[j] - 0.5 * fl1_fx * tey_xyz_zz_1[j];

                tez_xxyz_xzz_0[j] = pa_x[j] * tez_xyz_xzz_0[j] - pc_x[j] * tez_xyz_xzz_1[j] + 0.5 * fl1_fx * tez_yz_xzz_0[j] -
                                    0.5 * fl1_fx * tez_yz_xzz_1[j] + 0.5 * fl1_fx * tez_xyz_zz_0[j] - 0.5 * fl1_fx * tez_xyz_zz_1[j];

                tex_xxyz_yyy_0[j] = pa_x[j] * tex_xyz_yyy_0[j] - pc_x[j] * tex_xyz_yyy_1[j] + 0.5 * fl1_fx * tex_yz_yyy_0[j] -
                                    0.5 * fl1_fx * tex_yz_yyy_1[j] + ta_xyz_yyy_1[j];

                tey_xxyz_yyy_0[j] =
                    pa_x[j] * tey_xyz_yyy_0[j] - pc_x[j] * tey_xyz_yyy_1[j] + 0.5 * fl1_fx * tey_yz_yyy_0[j] - 0.5 * fl1_fx * tey_yz_yyy_1[j];

                tez_xxyz_yyy_0[j] =
                    pa_x[j] * tez_xyz_yyy_0[j] - pc_x[j] * tez_xyz_yyy_1[j] + 0.5 * fl1_fx * tez_yz_yyy_0[j] - 0.5 * fl1_fx * tez_yz_yyy_1[j];

                tex_xxyz_yyz_0[j] = pa_x[j] * tex_xyz_yyz_0[j] - pc_x[j] * tex_xyz_yyz_1[j] + 0.5 * fl1_fx * tex_yz_yyz_0[j] -
                                    0.5 * fl1_fx * tex_yz_yyz_1[j] + ta_xyz_yyz_1[j];

                tey_xxyz_yyz_0[j] =
                    pa_x[j] * tey_xyz_yyz_0[j] - pc_x[j] * tey_xyz_yyz_1[j] + 0.5 * fl1_fx * tey_yz_yyz_0[j] - 0.5 * fl1_fx * tey_yz_yyz_1[j];

                tez_xxyz_yyz_0[j] =
                    pa_x[j] * tez_xyz_yyz_0[j] - pc_x[j] * tez_xyz_yyz_1[j] + 0.5 * fl1_fx * tez_yz_yyz_0[j] - 0.5 * fl1_fx * tez_yz_yyz_1[j];

                tex_xxyz_yzz_0[j] = pa_x[j] * tex_xyz_yzz_0[j] - pc_x[j] * tex_xyz_yzz_1[j] + 0.5 * fl1_fx * tex_yz_yzz_0[j] -
                                    0.5 * fl1_fx * tex_yz_yzz_1[j] + ta_xyz_yzz_1[j];

                tey_xxyz_yzz_0[j] =
                    pa_x[j] * tey_xyz_yzz_0[j] - pc_x[j] * tey_xyz_yzz_1[j] + 0.5 * fl1_fx * tey_yz_yzz_0[j] - 0.5 * fl1_fx * tey_yz_yzz_1[j];

                tez_xxyz_yzz_0[j] =
                    pa_x[j] * tez_xyz_yzz_0[j] - pc_x[j] * tez_xyz_yzz_1[j] + 0.5 * fl1_fx * tez_yz_yzz_0[j] - 0.5 * fl1_fx * tez_yz_yzz_1[j];

                tex_xxyz_zzz_0[j] = pa_x[j] * tex_xyz_zzz_0[j] - pc_x[j] * tex_xyz_zzz_1[j] + 0.5 * fl1_fx * tex_yz_zzz_0[j] -
                                    0.5 * fl1_fx * tex_yz_zzz_1[j] + ta_xyz_zzz_1[j];

                tey_xxyz_zzz_0[j] =
                    pa_x[j] * tey_xyz_zzz_0[j] - pc_x[j] * tey_xyz_zzz_1[j] + 0.5 * fl1_fx * tey_yz_zzz_0[j] - 0.5 * fl1_fx * tey_yz_zzz_1[j];

                tez_xxyz_zzz_0[j] =
                    pa_x[j] * tez_xyz_zzz_0[j] - pc_x[j] * tez_xyz_zzz_1[j] + 0.5 * fl1_fx * tez_yz_zzz_0[j] - 0.5 * fl1_fx * tez_yz_zzz_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForGF_150_200(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pcDistances,
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

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_4_3_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_2_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tex_zz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 50);

            auto tey_zz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 50);

            auto tez_zz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 50);

            auto tex_zz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 51);

            auto tey_zz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 51);

            auto tez_zz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 51);

            auto tex_zz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 52);

            auto tey_zz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 52);

            auto tez_zz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 52);

            auto tex_zz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 53);

            auto tey_zz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 53);

            auto tez_zz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 53);

            auto tex_zz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 54);

            auto tey_zz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 54);

            auto tez_zz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 54);

            auto tex_zz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 55);

            auto tey_zz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 55);

            auto tez_zz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 55);

            auto tex_zz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 56);

            auto tey_zz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 56);

            auto tez_zz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 56);

            auto tex_zz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 57);

            auto tey_zz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 57);

            auto tez_zz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 57);

            auto tex_zz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 58);

            auto tey_zz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 58);

            auto tez_zz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 58);

            auto tex_zz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 59);

            auto tey_zz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 59);

            auto tez_zz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 59);

            auto tex_zz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 50);

            auto tey_zz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 50);

            auto tez_zz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 50);

            auto tex_zz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 51);

            auto tey_zz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 51);

            auto tez_zz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 51);

            auto tex_zz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 52);

            auto tey_zz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 52);

            auto tez_zz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 52);

            auto tex_zz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 53);

            auto tey_zz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 53);

            auto tez_zz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 53);

            auto tex_zz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 54);

            auto tey_zz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 54);

            auto tez_zz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 54);

            auto tex_zz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 55);

            auto tey_zz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 55);

            auto tez_zz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 55);

            auto tex_zz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 56);

            auto tey_zz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 56);

            auto tez_zz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 56);

            auto tex_zz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 57);

            auto tey_zz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 57);

            auto tez_zz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 57);

            auto tex_zz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 58);

            auto tey_zz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 58);

            auto tez_zz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 58);

            auto tex_zz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 59);

            auto tey_zz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 59);

            auto tez_zz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 59);

            auto tex_xzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 30);

            auto tey_xzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 30);

            auto tez_xzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 30);

            auto tex_xzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 31);

            auto tey_xzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 31);

            auto tez_xzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 31);

            auto tex_xzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 32);

            auto tey_xzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 32);

            auto tez_xzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 32);

            auto tex_xzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 33);

            auto tey_xzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 33);

            auto tez_xzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 33);

            auto tex_xzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 34);

            auto tey_xzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 34);

            auto tez_xzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 34);

            auto tex_xzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 35);

            auto tey_xzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 35);

            auto tez_xzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 35);

            auto tex_yyy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 36);

            auto tey_yyy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 36);

            auto tez_yyy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 36);

            auto tex_yyy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 37);

            auto tey_yyy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 37);

            auto tez_yyy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 37);

            auto tex_yyy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 38);

            auto tey_yyy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 38);

            auto tez_yyy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 38);

            auto tex_yyy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 39);

            auto tey_yyy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 39);

            auto tez_yyy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 39);

            auto tex_yyy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 40);

            auto tey_yyy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 40);

            auto tez_yyy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 40);

            auto tex_yyy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 41);

            auto tey_yyy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 41);

            auto tez_yyy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 41);

            auto tex_xzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 30);

            auto tey_xzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 30);

            auto tez_xzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 30);

            auto tex_xzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 31);

            auto tey_xzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 31);

            auto tez_xzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 31);

            auto tex_xzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 32);

            auto tey_xzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 32);

            auto tez_xzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 32);

            auto tex_xzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 33);

            auto tey_xzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 33);

            auto tez_xzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 33);

            auto tex_xzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 34);

            auto tey_xzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 34);

            auto tez_xzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 34);

            auto tex_xzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 35);

            auto tey_xzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 35);

            auto tez_xzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 35);

            auto tex_yyy_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 36);

            auto tey_yyy_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 36);

            auto tez_yyy_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 36);

            auto tex_yyy_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 37);

            auto tey_yyy_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 37);

            auto tez_yyy_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 37);

            auto tex_yyy_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 38);

            auto tey_yyy_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 38);

            auto tez_yyy_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 38);

            auto tex_yyy_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 39);

            auto tey_yyy_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 39);

            auto tez_yyy_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 39);

            auto tex_yyy_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 40);

            auto tey_yyy_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 40);

            auto tez_yyy_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 40);

            auto tex_yyy_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 41);

            auto tey_yyy_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 41);

            auto tez_yyy_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 41);

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

            auto ta_yyy_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 60);

            auto ta_yyy_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 61);

            auto ta_yyy_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 62);

            auto ta_yyy_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 63);

            auto ta_yyy_xyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 64);

            auto ta_yyy_xzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 65);

            auto ta_yyy_yyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 66);

            // set up pointers to integrals

            auto tex_xxzz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 50);

            auto tey_xxzz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 50);

            auto tez_xxzz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 50);

            auto tex_xxzz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 51);

            auto tey_xxzz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 51);

            auto tez_xxzz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 51);

            auto tex_xxzz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 52);

            auto tey_xxzz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 52);

            auto tez_xxzz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 52);

            auto tex_xxzz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 53);

            auto tey_xxzz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 53);

            auto tez_xxzz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 53);

            auto tex_xxzz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 54);

            auto tey_xxzz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 54);

            auto tez_xxzz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 54);

            auto tex_xxzz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 55);

            auto tey_xxzz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 55);

            auto tez_xxzz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 55);

            auto tex_xxzz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 56);

            auto tey_xxzz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 56);

            auto tez_xxzz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 56);

            auto tex_xxzz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 57);

            auto tey_xxzz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 57);

            auto tez_xxzz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 57);

            auto tex_xxzz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 58);

            auto tey_xxzz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 58);

            auto tez_xxzz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 58);

            auto tex_xxzz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 59);

            auto tey_xxzz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 59);

            auto tez_xxzz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 59);

            auto tex_xyyy_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 60);

            auto tey_xyyy_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 60);

            auto tez_xyyy_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 60);

            auto tex_xyyy_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 61);

            auto tey_xyyy_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 61);

            auto tez_xyyy_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 61);

            auto tex_xyyy_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 62);

            auto tey_xyyy_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 62);

            auto tez_xyyy_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 62);

            auto tex_xyyy_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 63);

            auto tey_xyyy_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 63);

            auto tez_xyyy_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 63);

            auto tex_xyyy_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 64);

            auto tey_xyyy_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 64);

            auto tez_xyyy_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 64);

            auto tex_xyyy_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 65);

            auto tey_xyyy_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 65);

            auto tez_xyyy_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 65);

            auto tex_xyyy_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 66);

            auto tey_xyyy_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 66);

            // Batch of Integrals (150,200)

            #pragma omp simd aligned(fx, pa_x, pc_x, ta_xzz_xxx_1, ta_xzz_xxy_1, ta_xzz_xxz_1, ta_xzz_xyy_1, \
                                         ta_xzz_xyz_1, ta_xzz_xzz_1, ta_xzz_yyy_1, ta_xzz_yyz_1, ta_xzz_yzz_1, ta_xzz_zzz_1, \
                                         ta_yyy_xxx_1, ta_yyy_xxy_1, ta_yyy_xxz_1, ta_yyy_xyy_1, ta_yyy_xyz_1, ta_yyy_xzz_1, \
                                         ta_yyy_yyy_1, tex_xxzz_xxx_0, tex_xxzz_xxy_0, tex_xxzz_xxz_0, tex_xxzz_xyy_0, \
                                         tex_xxzz_xyz_0, tex_xxzz_xzz_0, tex_xxzz_yyy_0, tex_xxzz_yyz_0, tex_xxzz_yzz_0, \
                                         tex_xxzz_zzz_0, tex_xyyy_xxx_0, tex_xyyy_xxy_0, tex_xyyy_xxz_0, tex_xyyy_xyy_0, \
                                         tex_xyyy_xyz_0, tex_xyyy_xzz_0, tex_xyyy_yyy_0, tex_xzz_xx_0, tex_xzz_xx_1, \
                                         tex_xzz_xxx_0, tex_xzz_xxx_1, tex_xzz_xxy_0, tex_xzz_xxy_1, tex_xzz_xxz_0, \
                                         tex_xzz_xxz_1, tex_xzz_xy_0, tex_xzz_xy_1, tex_xzz_xyy_0, tex_xzz_xyy_1, \
                                         tex_xzz_xyz_0, tex_xzz_xyz_1, tex_xzz_xz_0, tex_xzz_xz_1, tex_xzz_xzz_0, \
                                         tex_xzz_xzz_1, tex_xzz_yy_0, tex_xzz_yy_1, tex_xzz_yyy_0, tex_xzz_yyy_1, \
                                         tex_xzz_yyz_0, tex_xzz_yyz_1, tex_xzz_yz_0, tex_xzz_yz_1, tex_xzz_yzz_0, \
                                         tex_xzz_yzz_1, tex_xzz_zz_0, tex_xzz_zz_1, tex_xzz_zzz_0, tex_xzz_zzz_1, \
                                         tex_yyy_xx_0, tex_yyy_xx_1, tex_yyy_xxx_0, tex_yyy_xxx_1, tex_yyy_xxy_0, \
                                         tex_yyy_xxy_1, tex_yyy_xxz_0, tex_yyy_xxz_1, tex_yyy_xy_0, tex_yyy_xy_1, \
                                         tex_yyy_xyy_0, tex_yyy_xyy_1, tex_yyy_xyz_0, tex_yyy_xyz_1, tex_yyy_xz_0, \
                                         tex_yyy_xz_1, tex_yyy_xzz_0, tex_yyy_xzz_1, tex_yyy_yy_0, tex_yyy_yy_1, \
                                         tex_yyy_yyy_0, tex_yyy_yyy_1, tex_yyy_yz_0, tex_yyy_yz_1, tex_yyy_zz_0, \
                                         tex_yyy_zz_1, tex_zz_xxx_0, tex_zz_xxx_1, tex_zz_xxy_0, tex_zz_xxy_1, tex_zz_xxz_0, \
                                         tex_zz_xxz_1, tex_zz_xyy_0, tex_zz_xyy_1, tex_zz_xyz_0, tex_zz_xyz_1, tex_zz_xzz_0, \
                                         tex_zz_xzz_1, tex_zz_yyy_0, tex_zz_yyy_1, tex_zz_yyz_0, tex_zz_yyz_1, tex_zz_yzz_0, \
                                         tex_zz_yzz_1, tex_zz_zzz_0, tex_zz_zzz_1, tey_xxzz_xxx_0, tey_xxzz_xxy_0, \
                                         tey_xxzz_xxz_0, tey_xxzz_xyy_0, tey_xxzz_xyz_0, tey_xxzz_xzz_0, tey_xxzz_yyy_0, \
                                         tey_xxzz_yyz_0, tey_xxzz_yzz_0, tey_xxzz_zzz_0, tey_xyyy_xxx_0, tey_xyyy_xxy_0, \
                                         tey_xyyy_xxz_0, tey_xyyy_xyy_0, tey_xyyy_xyz_0, tey_xyyy_xzz_0, tey_xyyy_yyy_0, \
                                         tey_xzz_xx_0, tey_xzz_xx_1, tey_xzz_xxx_0, tey_xzz_xxx_1, tey_xzz_xxy_0, \
                                         tey_xzz_xxy_1, tey_xzz_xxz_0, tey_xzz_xxz_1, tey_xzz_xy_0, tey_xzz_xy_1, \
                                         tey_xzz_xyy_0, tey_xzz_xyy_1, tey_xzz_xyz_0, tey_xzz_xyz_1, tey_xzz_xz_0, \
                                         tey_xzz_xz_1, tey_xzz_xzz_0, tey_xzz_xzz_1, tey_xzz_yy_0, tey_xzz_yy_1, \
                                         tey_xzz_yyy_0, tey_xzz_yyy_1, tey_xzz_yyz_0, tey_xzz_yyz_1, tey_xzz_yz_0, \
                                         tey_xzz_yz_1, tey_xzz_yzz_0, tey_xzz_yzz_1, tey_xzz_zz_0, tey_xzz_zz_1, \
                                         tey_xzz_zzz_0, tey_xzz_zzz_1, tey_yyy_xx_0, tey_yyy_xx_1, tey_yyy_xxx_0, \
                                         tey_yyy_xxx_1, tey_yyy_xxy_0, tey_yyy_xxy_1, tey_yyy_xxz_0, tey_yyy_xxz_1, \
                                         tey_yyy_xy_0, tey_yyy_xy_1, tey_yyy_xyy_0, tey_yyy_xyy_1, tey_yyy_xyz_0, \
                                         tey_yyy_xyz_1, tey_yyy_xz_0, tey_yyy_xz_1, tey_yyy_xzz_0, tey_yyy_xzz_1, \
                                         tey_yyy_yy_0, tey_yyy_yy_1, tey_yyy_yyy_0, tey_yyy_yyy_1, tey_yyy_yz_0, \
                                         tey_yyy_yz_1, tey_yyy_zz_0, tey_yyy_zz_1, tey_zz_xxx_0, tey_zz_xxx_1, tey_zz_xxy_0, \
                                         tey_zz_xxy_1, tey_zz_xxz_0, tey_zz_xxz_1, tey_zz_xyy_0, tey_zz_xyy_1, tey_zz_xyz_0, \
                                         tey_zz_xyz_1, tey_zz_xzz_0, tey_zz_xzz_1, tey_zz_yyy_0, tey_zz_yyy_1, tey_zz_yyz_0, \
                                         tey_zz_yyz_1, tey_zz_yzz_0, tey_zz_yzz_1, tey_zz_zzz_0, tey_zz_zzz_1, \
                                         tez_xxzz_xxx_0, tez_xxzz_xxy_0, tez_xxzz_xxz_0, tez_xxzz_xyy_0, tez_xxzz_xyz_0, \
                                         tez_xxzz_xzz_0, tez_xxzz_yyy_0, tez_xxzz_yyz_0, tez_xxzz_yzz_0, tez_xxzz_zzz_0, \
                                         tez_xyyy_xxx_0, tez_xyyy_xxy_0, tez_xyyy_xxz_0, tez_xyyy_xyy_0, tez_xyyy_xyz_0, \
                                         tez_xyyy_xzz_0, tez_xzz_xx_0, tez_xzz_xx_1, tez_xzz_xxx_0, tez_xzz_xxx_1, \
                                         tez_xzz_xxy_0, tez_xzz_xxy_1, tez_xzz_xxz_0, tez_xzz_xxz_1, tez_xzz_xy_0, \
                                         tez_xzz_xy_1, tez_xzz_xyy_0, tez_xzz_xyy_1, tez_xzz_xyz_0, tez_xzz_xyz_1, \
                                         tez_xzz_xz_0, tez_xzz_xz_1, tez_xzz_xzz_0, tez_xzz_xzz_1, tez_xzz_yy_0, \
                                         tez_xzz_yy_1, tez_xzz_yyy_0, tez_xzz_yyy_1, tez_xzz_yyz_0, tez_xzz_yyz_1, \
                                         tez_xzz_yz_0, tez_xzz_yz_1, tez_xzz_yzz_0, tez_xzz_yzz_1, tez_xzz_zz_0, \
                                         tez_xzz_zz_1, tez_xzz_zzz_0, tez_xzz_zzz_1, tez_yyy_xx_0, tez_yyy_xx_1, \
                                         tez_yyy_xxx_0, tez_yyy_xxx_1, tez_yyy_xxy_0, tez_yyy_xxy_1, tez_yyy_xxz_0, \
                                         tez_yyy_xxz_1, tez_yyy_xy_0, tez_yyy_xy_1, tez_yyy_xyy_0, tez_yyy_xyy_1, \
                                         tez_yyy_xyz_0, tez_yyy_xyz_1, tez_yyy_xz_0, tez_yyy_xz_1, tez_yyy_xzz_0, \
                                         tez_yyy_xzz_1, tez_yyy_yy_0, tez_yyy_yy_1, tez_yyy_yz_0, tez_yyy_yz_1, tez_yyy_zz_0, \
                                         tez_yyy_zz_1, tez_zz_xxx_0, tez_zz_xxx_1, tez_zz_xxy_0, tez_zz_xxy_1, tez_zz_xxz_0, \
                                         tez_zz_xxz_1, tez_zz_xyy_0, tez_zz_xyy_1, tez_zz_xyz_0, tez_zz_xyz_1, tez_zz_xzz_0, \
                                         tez_zz_xzz_1, tez_zz_yyy_0, tez_zz_yyy_1, tez_zz_yyz_0, tez_zz_yyz_1, tez_zz_yzz_0, \
                                         tez_zz_yzz_1, tez_zz_zzz_0, tez_zz_zzz_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_xxzz_xxx_0[j] = pa_x[j] * tex_xzz_xxx_0[j] - pc_x[j] * tex_xzz_xxx_1[j] + 0.5 * fl1_fx * tex_zz_xxx_0[j] -
                                    0.5 * fl1_fx * tex_zz_xxx_1[j] + 1.5 * fl1_fx * tex_xzz_xx_0[j] - 1.5 * fl1_fx * tex_xzz_xx_1[j] +
                                    ta_xzz_xxx_1[j];

                tey_xxzz_xxx_0[j] = pa_x[j] * tey_xzz_xxx_0[j] - pc_x[j] * tey_xzz_xxx_1[j] + 0.5 * fl1_fx * tey_zz_xxx_0[j] -
                                    0.5 * fl1_fx * tey_zz_xxx_1[j] + 1.5 * fl1_fx * tey_xzz_xx_0[j] - 1.5 * fl1_fx * tey_xzz_xx_1[j];

                tez_xxzz_xxx_0[j] = pa_x[j] * tez_xzz_xxx_0[j] - pc_x[j] * tez_xzz_xxx_1[j] + 0.5 * fl1_fx * tez_zz_xxx_0[j] -
                                    0.5 * fl1_fx * tez_zz_xxx_1[j] + 1.5 * fl1_fx * tez_xzz_xx_0[j] - 1.5 * fl1_fx * tez_xzz_xx_1[j];

                tex_xxzz_xxy_0[j] = pa_x[j] * tex_xzz_xxy_0[j] - pc_x[j] * tex_xzz_xxy_1[j] + 0.5 * fl1_fx * tex_zz_xxy_0[j] -
                                    0.5 * fl1_fx * tex_zz_xxy_1[j] + fl1_fx * tex_xzz_xy_0[j] - fl1_fx * tex_xzz_xy_1[j] + ta_xzz_xxy_1[j];

                tey_xxzz_xxy_0[j] = pa_x[j] * tey_xzz_xxy_0[j] - pc_x[j] * tey_xzz_xxy_1[j] + 0.5 * fl1_fx * tey_zz_xxy_0[j] -
                                    0.5 * fl1_fx * tey_zz_xxy_1[j] + fl1_fx * tey_xzz_xy_0[j] - fl1_fx * tey_xzz_xy_1[j];

                tez_xxzz_xxy_0[j] = pa_x[j] * tez_xzz_xxy_0[j] - pc_x[j] * tez_xzz_xxy_1[j] + 0.5 * fl1_fx * tez_zz_xxy_0[j] -
                                    0.5 * fl1_fx * tez_zz_xxy_1[j] + fl1_fx * tez_xzz_xy_0[j] - fl1_fx * tez_xzz_xy_1[j];

                tex_xxzz_xxz_0[j] = pa_x[j] * tex_xzz_xxz_0[j] - pc_x[j] * tex_xzz_xxz_1[j] + 0.5 * fl1_fx * tex_zz_xxz_0[j] -
                                    0.5 * fl1_fx * tex_zz_xxz_1[j] + fl1_fx * tex_xzz_xz_0[j] - fl1_fx * tex_xzz_xz_1[j] + ta_xzz_xxz_1[j];

                tey_xxzz_xxz_0[j] = pa_x[j] * tey_xzz_xxz_0[j] - pc_x[j] * tey_xzz_xxz_1[j] + 0.5 * fl1_fx * tey_zz_xxz_0[j] -
                                    0.5 * fl1_fx * tey_zz_xxz_1[j] + fl1_fx * tey_xzz_xz_0[j] - fl1_fx * tey_xzz_xz_1[j];

                tez_xxzz_xxz_0[j] = pa_x[j] * tez_xzz_xxz_0[j] - pc_x[j] * tez_xzz_xxz_1[j] + 0.5 * fl1_fx * tez_zz_xxz_0[j] -
                                    0.5 * fl1_fx * tez_zz_xxz_1[j] + fl1_fx * tez_xzz_xz_0[j] - fl1_fx * tez_xzz_xz_1[j];

                tex_xxzz_xyy_0[j] = pa_x[j] * tex_xzz_xyy_0[j] - pc_x[j] * tex_xzz_xyy_1[j] + 0.5 * fl1_fx * tex_zz_xyy_0[j] -
                                    0.5 * fl1_fx * tex_zz_xyy_1[j] + 0.5 * fl1_fx * tex_xzz_yy_0[j] - 0.5 * fl1_fx * tex_xzz_yy_1[j] +
                                    ta_xzz_xyy_1[j];

                tey_xxzz_xyy_0[j] = pa_x[j] * tey_xzz_xyy_0[j] - pc_x[j] * tey_xzz_xyy_1[j] + 0.5 * fl1_fx * tey_zz_xyy_0[j] -
                                    0.5 * fl1_fx * tey_zz_xyy_1[j] + 0.5 * fl1_fx * tey_xzz_yy_0[j] - 0.5 * fl1_fx * tey_xzz_yy_1[j];

                tez_xxzz_xyy_0[j] = pa_x[j] * tez_xzz_xyy_0[j] - pc_x[j] * tez_xzz_xyy_1[j] + 0.5 * fl1_fx * tez_zz_xyy_0[j] -
                                    0.5 * fl1_fx * tez_zz_xyy_1[j] + 0.5 * fl1_fx * tez_xzz_yy_0[j] - 0.5 * fl1_fx * tez_xzz_yy_1[j];

                tex_xxzz_xyz_0[j] = pa_x[j] * tex_xzz_xyz_0[j] - pc_x[j] * tex_xzz_xyz_1[j] + 0.5 * fl1_fx * tex_zz_xyz_0[j] -
                                    0.5 * fl1_fx * tex_zz_xyz_1[j] + 0.5 * fl1_fx * tex_xzz_yz_0[j] - 0.5 * fl1_fx * tex_xzz_yz_1[j] +
                                    ta_xzz_xyz_1[j];

                tey_xxzz_xyz_0[j] = pa_x[j] * tey_xzz_xyz_0[j] - pc_x[j] * tey_xzz_xyz_1[j] + 0.5 * fl1_fx * tey_zz_xyz_0[j] -
                                    0.5 * fl1_fx * tey_zz_xyz_1[j] + 0.5 * fl1_fx * tey_xzz_yz_0[j] - 0.5 * fl1_fx * tey_xzz_yz_1[j];

                tez_xxzz_xyz_0[j] = pa_x[j] * tez_xzz_xyz_0[j] - pc_x[j] * tez_xzz_xyz_1[j] + 0.5 * fl1_fx * tez_zz_xyz_0[j] -
                                    0.5 * fl1_fx * tez_zz_xyz_1[j] + 0.5 * fl1_fx * tez_xzz_yz_0[j] - 0.5 * fl1_fx * tez_xzz_yz_1[j];

                tex_xxzz_xzz_0[j] = pa_x[j] * tex_xzz_xzz_0[j] - pc_x[j] * tex_xzz_xzz_1[j] + 0.5 * fl1_fx * tex_zz_xzz_0[j] -
                                    0.5 * fl1_fx * tex_zz_xzz_1[j] + 0.5 * fl1_fx * tex_xzz_zz_0[j] - 0.5 * fl1_fx * tex_xzz_zz_1[j] +
                                    ta_xzz_xzz_1[j];

                tey_xxzz_xzz_0[j] = pa_x[j] * tey_xzz_xzz_0[j] - pc_x[j] * tey_xzz_xzz_1[j] + 0.5 * fl1_fx * tey_zz_xzz_0[j] -
                                    0.5 * fl1_fx * tey_zz_xzz_1[j] + 0.5 * fl1_fx * tey_xzz_zz_0[j] - 0.5 * fl1_fx * tey_xzz_zz_1[j];

                tez_xxzz_xzz_0[j] = pa_x[j] * tez_xzz_xzz_0[j] - pc_x[j] * tez_xzz_xzz_1[j] + 0.5 * fl1_fx * tez_zz_xzz_0[j] -
                                    0.5 * fl1_fx * tez_zz_xzz_1[j] + 0.5 * fl1_fx * tez_xzz_zz_0[j] - 0.5 * fl1_fx * tez_xzz_zz_1[j];

                tex_xxzz_yyy_0[j] = pa_x[j] * tex_xzz_yyy_0[j] - pc_x[j] * tex_xzz_yyy_1[j] + 0.5 * fl1_fx * tex_zz_yyy_0[j] -
                                    0.5 * fl1_fx * tex_zz_yyy_1[j] + ta_xzz_yyy_1[j];

                tey_xxzz_yyy_0[j] =
                    pa_x[j] * tey_xzz_yyy_0[j] - pc_x[j] * tey_xzz_yyy_1[j] + 0.5 * fl1_fx * tey_zz_yyy_0[j] - 0.5 * fl1_fx * tey_zz_yyy_1[j];

                tez_xxzz_yyy_0[j] =
                    pa_x[j] * tez_xzz_yyy_0[j] - pc_x[j] * tez_xzz_yyy_1[j] + 0.5 * fl1_fx * tez_zz_yyy_0[j] - 0.5 * fl1_fx * tez_zz_yyy_1[j];

                tex_xxzz_yyz_0[j] = pa_x[j] * tex_xzz_yyz_0[j] - pc_x[j] * tex_xzz_yyz_1[j] + 0.5 * fl1_fx * tex_zz_yyz_0[j] -
                                    0.5 * fl1_fx * tex_zz_yyz_1[j] + ta_xzz_yyz_1[j];

                tey_xxzz_yyz_0[j] =
                    pa_x[j] * tey_xzz_yyz_0[j] - pc_x[j] * tey_xzz_yyz_1[j] + 0.5 * fl1_fx * tey_zz_yyz_0[j] - 0.5 * fl1_fx * tey_zz_yyz_1[j];

                tez_xxzz_yyz_0[j] =
                    pa_x[j] * tez_xzz_yyz_0[j] - pc_x[j] * tez_xzz_yyz_1[j] + 0.5 * fl1_fx * tez_zz_yyz_0[j] - 0.5 * fl1_fx * tez_zz_yyz_1[j];

                tex_xxzz_yzz_0[j] = pa_x[j] * tex_xzz_yzz_0[j] - pc_x[j] * tex_xzz_yzz_1[j] + 0.5 * fl1_fx * tex_zz_yzz_0[j] -
                                    0.5 * fl1_fx * tex_zz_yzz_1[j] + ta_xzz_yzz_1[j];

                tey_xxzz_yzz_0[j] =
                    pa_x[j] * tey_xzz_yzz_0[j] - pc_x[j] * tey_xzz_yzz_1[j] + 0.5 * fl1_fx * tey_zz_yzz_0[j] - 0.5 * fl1_fx * tey_zz_yzz_1[j];

                tez_xxzz_yzz_0[j] =
                    pa_x[j] * tez_xzz_yzz_0[j] - pc_x[j] * tez_xzz_yzz_1[j] + 0.5 * fl1_fx * tez_zz_yzz_0[j] - 0.5 * fl1_fx * tez_zz_yzz_1[j];

                tex_xxzz_zzz_0[j] = pa_x[j] * tex_xzz_zzz_0[j] - pc_x[j] * tex_xzz_zzz_1[j] + 0.5 * fl1_fx * tex_zz_zzz_0[j] -
                                    0.5 * fl1_fx * tex_zz_zzz_1[j] + ta_xzz_zzz_1[j];

                tey_xxzz_zzz_0[j] =
                    pa_x[j] * tey_xzz_zzz_0[j] - pc_x[j] * tey_xzz_zzz_1[j] + 0.5 * fl1_fx * tey_zz_zzz_0[j] - 0.5 * fl1_fx * tey_zz_zzz_1[j];

                tez_xxzz_zzz_0[j] =
                    pa_x[j] * tez_xzz_zzz_0[j] - pc_x[j] * tez_xzz_zzz_1[j] + 0.5 * fl1_fx * tez_zz_zzz_0[j] - 0.5 * fl1_fx * tez_zz_zzz_1[j];

                tex_xyyy_xxx_0[j] = pa_x[j] * tex_yyy_xxx_0[j] - pc_x[j] * tex_yyy_xxx_1[j] + 1.5 * fl1_fx * tex_yyy_xx_0[j] -
                                    1.5 * fl1_fx * tex_yyy_xx_1[j] + ta_yyy_xxx_1[j];

                tey_xyyy_xxx_0[j] =
                    pa_x[j] * tey_yyy_xxx_0[j] - pc_x[j] * tey_yyy_xxx_1[j] + 1.5 * fl1_fx * tey_yyy_xx_0[j] - 1.5 * fl1_fx * tey_yyy_xx_1[j];

                tez_xyyy_xxx_0[j] =
                    pa_x[j] * tez_yyy_xxx_0[j] - pc_x[j] * tez_yyy_xxx_1[j] + 1.5 * fl1_fx * tez_yyy_xx_0[j] - 1.5 * fl1_fx * tez_yyy_xx_1[j];

                tex_xyyy_xxy_0[j] =
                    pa_x[j] * tex_yyy_xxy_0[j] - pc_x[j] * tex_yyy_xxy_1[j] + fl1_fx * tex_yyy_xy_0[j] - fl1_fx * tex_yyy_xy_1[j] + ta_yyy_xxy_1[j];

                tey_xyyy_xxy_0[j] = pa_x[j] * tey_yyy_xxy_0[j] - pc_x[j] * tey_yyy_xxy_1[j] + fl1_fx * tey_yyy_xy_0[j] - fl1_fx * tey_yyy_xy_1[j];

                tez_xyyy_xxy_0[j] = pa_x[j] * tez_yyy_xxy_0[j] - pc_x[j] * tez_yyy_xxy_1[j] + fl1_fx * tez_yyy_xy_0[j] - fl1_fx * tez_yyy_xy_1[j];

                tex_xyyy_xxz_0[j] =
                    pa_x[j] * tex_yyy_xxz_0[j] - pc_x[j] * tex_yyy_xxz_1[j] + fl1_fx * tex_yyy_xz_0[j] - fl1_fx * tex_yyy_xz_1[j] + ta_yyy_xxz_1[j];

                tey_xyyy_xxz_0[j] = pa_x[j] * tey_yyy_xxz_0[j] - pc_x[j] * tey_yyy_xxz_1[j] + fl1_fx * tey_yyy_xz_0[j] - fl1_fx * tey_yyy_xz_1[j];

                tez_xyyy_xxz_0[j] = pa_x[j] * tez_yyy_xxz_0[j] - pc_x[j] * tez_yyy_xxz_1[j] + fl1_fx * tez_yyy_xz_0[j] - fl1_fx * tez_yyy_xz_1[j];

                tex_xyyy_xyy_0[j] = pa_x[j] * tex_yyy_xyy_0[j] - pc_x[j] * tex_yyy_xyy_1[j] + 0.5 * fl1_fx * tex_yyy_yy_0[j] -
                                    0.5 * fl1_fx * tex_yyy_yy_1[j] + ta_yyy_xyy_1[j];

                tey_xyyy_xyy_0[j] =
                    pa_x[j] * tey_yyy_xyy_0[j] - pc_x[j] * tey_yyy_xyy_1[j] + 0.5 * fl1_fx * tey_yyy_yy_0[j] - 0.5 * fl1_fx * tey_yyy_yy_1[j];

                tez_xyyy_xyy_0[j] =
                    pa_x[j] * tez_yyy_xyy_0[j] - pc_x[j] * tez_yyy_xyy_1[j] + 0.5 * fl1_fx * tez_yyy_yy_0[j] - 0.5 * fl1_fx * tez_yyy_yy_1[j];

                tex_xyyy_xyz_0[j] = pa_x[j] * tex_yyy_xyz_0[j] - pc_x[j] * tex_yyy_xyz_1[j] + 0.5 * fl1_fx * tex_yyy_yz_0[j] -
                                    0.5 * fl1_fx * tex_yyy_yz_1[j] + ta_yyy_xyz_1[j];

                tey_xyyy_xyz_0[j] =
                    pa_x[j] * tey_yyy_xyz_0[j] - pc_x[j] * tey_yyy_xyz_1[j] + 0.5 * fl1_fx * tey_yyy_yz_0[j] - 0.5 * fl1_fx * tey_yyy_yz_1[j];

                tez_xyyy_xyz_0[j] =
                    pa_x[j] * tez_yyy_xyz_0[j] - pc_x[j] * tez_yyy_xyz_1[j] + 0.5 * fl1_fx * tez_yyy_yz_0[j] - 0.5 * fl1_fx * tez_yyy_yz_1[j];

                tex_xyyy_xzz_0[j] = pa_x[j] * tex_yyy_xzz_0[j] - pc_x[j] * tex_yyy_xzz_1[j] + 0.5 * fl1_fx * tex_yyy_zz_0[j] -
                                    0.5 * fl1_fx * tex_yyy_zz_1[j] + ta_yyy_xzz_1[j];

                tey_xyyy_xzz_0[j] =
                    pa_x[j] * tey_yyy_xzz_0[j] - pc_x[j] * tey_yyy_xzz_1[j] + 0.5 * fl1_fx * tey_yyy_zz_0[j] - 0.5 * fl1_fx * tey_yyy_zz_1[j];

                tez_xyyy_xzz_0[j] =
                    pa_x[j] * tez_yyy_xzz_0[j] - pc_x[j] * tez_yyy_xzz_1[j] + 0.5 * fl1_fx * tez_yyy_zz_0[j] - 0.5 * fl1_fx * tez_yyy_zz_1[j];

                tex_xyyy_yyy_0[j] = pa_x[j] * tex_yyy_yyy_0[j] - pc_x[j] * tex_yyy_yyy_1[j] + ta_yyy_yyy_1[j];

                tey_xyyy_yyy_0[j] = pa_x[j] * tey_yyy_yyy_0[j] - pc_x[j] * tey_yyy_yyy_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForGF_200_250(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pcDistances,
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

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_4_3_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tez_yyy_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 66);

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

            auto tez_yyy_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 66);

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

            auto tex_yyz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 42);

            auto tey_yyz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 42);

            auto tez_yyz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 42);

            auto tex_yyz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 43);

            auto tey_yyz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 43);

            auto tez_yyz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 43);

            auto tex_yyz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 44);

            auto tey_yyz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 44);

            auto tez_yyz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 44);

            auto tex_yyz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 45);

            auto tey_yyz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 45);

            auto tez_yyz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 45);

            auto tex_yyz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 46);

            auto tey_yyz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 46);

            auto tez_yyz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 46);

            auto tex_yyz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 47);

            auto tey_yyz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 47);

            auto tez_yyz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 47);

            auto tex_yzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 48);

            auto tey_yzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 48);

            auto tez_yzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 48);

            auto tex_yzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 49);

            auto tey_yzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 49);

            auto tez_yzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 49);

            auto tex_yzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 50);

            auto tey_yzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 50);

            auto tez_yzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 50);

            auto tex_yzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 51);

            auto tex_yyz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 42);

            auto tey_yyz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 42);

            auto tez_yyz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 42);

            auto tex_yyz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 43);

            auto tey_yyz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 43);

            auto tez_yyz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 43);

            auto tex_yyz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 44);

            auto tey_yyz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 44);

            auto tez_yyz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 44);

            auto tex_yyz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 45);

            auto tey_yyz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 45);

            auto tez_yyz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 45);

            auto tex_yyz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 46);

            auto tey_yyz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 46);

            auto tez_yyz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 46);

            auto tex_yyz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 47);

            auto tey_yyz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 47);

            auto tez_yyz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 47);

            auto tex_yzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 48);

            auto tey_yzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 48);

            auto tez_yzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 48);

            auto tex_yzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 49);

            auto tey_yzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 49);

            auto tez_yzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 49);

            auto tex_yzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 50);

            auto tey_yzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 50);

            auto tez_yzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 50);

            auto tex_yzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 51);

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

            // set up pointers to integrals

            auto tez_xyyy_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 66);

            auto tex_xyyy_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 67);

            auto tey_xyyy_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 67);

            auto tez_xyyy_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 67);

            auto tex_xyyy_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 68);

            auto tey_xyyy_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 68);

            auto tez_xyyy_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 68);

            auto tex_xyyy_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 69);

            auto tey_xyyy_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 69);

            auto tez_xyyy_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 69);

            auto tex_xyyz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 70);

            auto tey_xyyz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 70);

            auto tez_xyyz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 70);

            auto tex_xyyz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 71);

            auto tey_xyyz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 71);

            auto tez_xyyz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 71);

            auto tex_xyyz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 72);

            auto tey_xyyz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 72);

            auto tez_xyyz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 72);

            auto tex_xyyz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 73);

            auto tey_xyyz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 73);

            auto tez_xyyz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 73);

            auto tex_xyyz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 74);

            auto tey_xyyz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 74);

            auto tez_xyyz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 74);

            auto tex_xyyz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 75);

            auto tey_xyyz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 75);

            auto tez_xyyz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 75);

            auto tex_xyyz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 76);

            auto tey_xyyz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 76);

            auto tez_xyyz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 76);

            auto tex_xyyz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 77);

            auto tey_xyyz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 77);

            auto tez_xyyz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 77);

            auto tex_xyyz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 78);

            auto tey_xyyz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 78);

            auto tez_xyyz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 78);

            auto tex_xyyz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 79);

            auto tey_xyyz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 79);

            auto tez_xyyz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 79);

            auto tex_xyzz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 80);

            auto tey_xyzz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 80);

            auto tez_xyzz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 80);

            auto tex_xyzz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 81);

            auto tey_xyzz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 81);

            auto tez_xyzz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 81);

            auto tex_xyzz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 82);

            auto tey_xyzz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 82);

            auto tez_xyzz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 82);

            auto tex_xyzz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 83);

            // Batch of Integrals (200,250)

            #pragma omp simd aligned(fx, pa_x, pc_x, ta_yyy_yyz_1, ta_yyy_yzz_1, ta_yyy_zzz_1, ta_yyz_xxx_1, \
                                         ta_yyz_xxy_1, ta_yyz_xxz_1, ta_yyz_xyy_1, ta_yyz_xyz_1, ta_yyz_xzz_1, ta_yyz_yyy_1, \
                                         ta_yyz_yyz_1, ta_yyz_yzz_1, ta_yyz_zzz_1, ta_yzz_xxx_1, ta_yzz_xxy_1, ta_yzz_xxz_1, \
                                         ta_yzz_xyy_1, tex_xyyy_yyz_0, tex_xyyy_yzz_0, tex_xyyy_zzz_0, tex_xyyz_xxx_0, \
                                         tex_xyyz_xxy_0, tex_xyyz_xxz_0, tex_xyyz_xyy_0, tex_xyyz_xyz_0, tex_xyyz_xzz_0, \
                                         tex_xyyz_yyy_0, tex_xyyz_yyz_0, tex_xyyz_yzz_0, tex_xyyz_zzz_0, tex_xyzz_xxx_0, \
                                         tex_xyzz_xxy_0, tex_xyzz_xxz_0, tex_xyzz_xyy_0, tex_yyy_yyz_0, tex_yyy_yyz_1, \
                                         tex_yyy_yzz_0, tex_yyy_yzz_1, tex_yyy_zzz_0, tex_yyy_zzz_1, tex_yyz_xx_0, \
                                         tex_yyz_xx_1, tex_yyz_xxx_0, tex_yyz_xxx_1, tex_yyz_xxy_0, tex_yyz_xxy_1, \
                                         tex_yyz_xxz_0, tex_yyz_xxz_1, tex_yyz_xy_0, tex_yyz_xy_1, tex_yyz_xyy_0, \
                                         tex_yyz_xyy_1, tex_yyz_xyz_0, tex_yyz_xyz_1, tex_yyz_xz_0, tex_yyz_xz_1, \
                                         tex_yyz_xzz_0, tex_yyz_xzz_1, tex_yyz_yy_0, tex_yyz_yy_1, tex_yyz_yyy_0, \
                                         tex_yyz_yyy_1, tex_yyz_yyz_0, tex_yyz_yyz_1, tex_yyz_yz_0, tex_yyz_yz_1, \
                                         tex_yyz_yzz_0, tex_yyz_yzz_1, tex_yyz_zz_0, tex_yyz_zz_1, tex_yyz_zzz_0, \
                                         tex_yyz_zzz_1, tex_yzz_xx_0, tex_yzz_xx_1, tex_yzz_xxx_0, tex_yzz_xxx_1, \
                                         tex_yzz_xxy_0, tex_yzz_xxy_1, tex_yzz_xxz_0, tex_yzz_xxz_1, tex_yzz_xy_0, \
                                         tex_yzz_xy_1, tex_yzz_xyy_0, tex_yzz_xyy_1, tex_yzz_xz_0, tex_yzz_xz_1, \
                                         tex_yzz_yy_0, tex_yzz_yy_1, tey_xyyy_yyz_0, tey_xyyy_yzz_0, tey_xyyy_zzz_0, \
                                         tey_xyyz_xxx_0, tey_xyyz_xxy_0, tey_xyyz_xxz_0, tey_xyyz_xyy_0, tey_xyyz_xyz_0, \
                                         tey_xyyz_xzz_0, tey_xyyz_yyy_0, tey_xyyz_yyz_0, tey_xyyz_yzz_0, tey_xyyz_zzz_0, \
                                         tey_xyzz_xxx_0, tey_xyzz_xxy_0, tey_xyzz_xxz_0, tey_yyy_yyz_0, tey_yyy_yyz_1, \
                                         tey_yyy_yzz_0, tey_yyy_yzz_1, tey_yyy_zzz_0, tey_yyy_zzz_1, tey_yyz_xx_0, \
                                         tey_yyz_xx_1, tey_yyz_xxx_0, tey_yyz_xxx_1, tey_yyz_xxy_0, tey_yyz_xxy_1, \
                                         tey_yyz_xxz_0, tey_yyz_xxz_1, tey_yyz_xy_0, tey_yyz_xy_1, tey_yyz_xyy_0, \
                                         tey_yyz_xyy_1, tey_yyz_xyz_0, tey_yyz_xyz_1, tey_yyz_xz_0, tey_yyz_xz_1, \
                                         tey_yyz_xzz_0, tey_yyz_xzz_1, tey_yyz_yy_0, tey_yyz_yy_1, tey_yyz_yyy_0, \
                                         tey_yyz_yyy_1, tey_yyz_yyz_0, tey_yyz_yyz_1, tey_yyz_yz_0, tey_yyz_yz_1, \
                                         tey_yyz_yzz_0, tey_yyz_yzz_1, tey_yyz_zz_0, tey_yyz_zz_1, tey_yyz_zzz_0, \
                                         tey_yyz_zzz_1, tey_yzz_xx_0, tey_yzz_xx_1, tey_yzz_xxx_0, tey_yzz_xxx_1, \
                                         tey_yzz_xxy_0, tey_yzz_xxy_1, tey_yzz_xxz_0, tey_yzz_xxz_1, tey_yzz_xy_0, \
                                         tey_yzz_xy_1, tey_yzz_xz_0, tey_yzz_xz_1, tez_xyyy_yyy_0, tez_xyyy_yyz_0, \
                                         tez_xyyy_yzz_0, tez_xyyy_zzz_0, tez_xyyz_xxx_0, tez_xyyz_xxy_0, tez_xyyz_xxz_0, \
                                         tez_xyyz_xyy_0, tez_xyyz_xyz_0, tez_xyyz_xzz_0, tez_xyyz_yyy_0, tez_xyyz_yyz_0, \
                                         tez_xyyz_yzz_0, tez_xyyz_zzz_0, tez_xyzz_xxx_0, tez_xyzz_xxy_0, tez_xyzz_xxz_0, \
                                         tez_yyy_yyy_0, tez_yyy_yyy_1, tez_yyy_yyz_0, tez_yyy_yyz_1, tez_yyy_yzz_0, \
                                         tez_yyy_yzz_1, tez_yyy_zzz_0, tez_yyy_zzz_1, tez_yyz_xx_0, tez_yyz_xx_1, \
                                         tez_yyz_xxx_0, tez_yyz_xxx_1, tez_yyz_xxy_0, tez_yyz_xxy_1, tez_yyz_xxz_0, \
                                         tez_yyz_xxz_1, tez_yyz_xy_0, tez_yyz_xy_1, tez_yyz_xyy_0, tez_yyz_xyy_1, \
                                         tez_yyz_xyz_0, tez_yyz_xyz_1, tez_yyz_xz_0, tez_yyz_xz_1, tez_yyz_xzz_0, \
                                         tez_yyz_xzz_1, tez_yyz_yy_0, tez_yyz_yy_1, tez_yyz_yyy_0, tez_yyz_yyy_1, \
                                         tez_yyz_yyz_0, tez_yyz_yyz_1, tez_yyz_yz_0, tez_yyz_yz_1, tez_yyz_yzz_0, \
                                         tez_yyz_yzz_1, tez_yyz_zz_0, tez_yyz_zz_1, tez_yyz_zzz_0, tez_yyz_zzz_1, \
                                         tez_yzz_xx_0, tez_yzz_xx_1, tez_yzz_xxx_0, tez_yzz_xxx_1, tez_yzz_xxy_0, \
                                         tez_yzz_xxy_1, tez_yzz_xxz_0, tez_yzz_xxz_1, tez_yzz_xy_0, tez_yzz_xy_1, \
                                         tez_yzz_xz_0, tez_yzz_xz_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tez_xyyy_yyy_0[j] = pa_x[j] * tez_yyy_yyy_0[j] - pc_x[j] * tez_yyy_yyy_1[j];

                tex_xyyy_yyz_0[j] = pa_x[j] * tex_yyy_yyz_0[j] - pc_x[j] * tex_yyy_yyz_1[j] + ta_yyy_yyz_1[j];

                tey_xyyy_yyz_0[j] = pa_x[j] * tey_yyy_yyz_0[j] - pc_x[j] * tey_yyy_yyz_1[j];

                tez_xyyy_yyz_0[j] = pa_x[j] * tez_yyy_yyz_0[j] - pc_x[j] * tez_yyy_yyz_1[j];

                tex_xyyy_yzz_0[j] = pa_x[j] * tex_yyy_yzz_0[j] - pc_x[j] * tex_yyy_yzz_1[j] + ta_yyy_yzz_1[j];

                tey_xyyy_yzz_0[j] = pa_x[j] * tey_yyy_yzz_0[j] - pc_x[j] * tey_yyy_yzz_1[j];

                tez_xyyy_yzz_0[j] = pa_x[j] * tez_yyy_yzz_0[j] - pc_x[j] * tez_yyy_yzz_1[j];

                tex_xyyy_zzz_0[j] = pa_x[j] * tex_yyy_zzz_0[j] - pc_x[j] * tex_yyy_zzz_1[j] + ta_yyy_zzz_1[j];

                tey_xyyy_zzz_0[j] = pa_x[j] * tey_yyy_zzz_0[j] - pc_x[j] * tey_yyy_zzz_1[j];

                tez_xyyy_zzz_0[j] = pa_x[j] * tez_yyy_zzz_0[j] - pc_x[j] * tez_yyy_zzz_1[j];

                tex_xyyz_xxx_0[j] = pa_x[j] * tex_yyz_xxx_0[j] - pc_x[j] * tex_yyz_xxx_1[j] + 1.5 * fl1_fx * tex_yyz_xx_0[j] -
                                    1.5 * fl1_fx * tex_yyz_xx_1[j] + ta_yyz_xxx_1[j];

                tey_xyyz_xxx_0[j] =
                    pa_x[j] * tey_yyz_xxx_0[j] - pc_x[j] * tey_yyz_xxx_1[j] + 1.5 * fl1_fx * tey_yyz_xx_0[j] - 1.5 * fl1_fx * tey_yyz_xx_1[j];

                tez_xyyz_xxx_0[j] =
                    pa_x[j] * tez_yyz_xxx_0[j] - pc_x[j] * tez_yyz_xxx_1[j] + 1.5 * fl1_fx * tez_yyz_xx_0[j] - 1.5 * fl1_fx * tez_yyz_xx_1[j];

                tex_xyyz_xxy_0[j] =
                    pa_x[j] * tex_yyz_xxy_0[j] - pc_x[j] * tex_yyz_xxy_1[j] + fl1_fx * tex_yyz_xy_0[j] - fl1_fx * tex_yyz_xy_1[j] + ta_yyz_xxy_1[j];

                tey_xyyz_xxy_0[j] = pa_x[j] * tey_yyz_xxy_0[j] - pc_x[j] * tey_yyz_xxy_1[j] + fl1_fx * tey_yyz_xy_0[j] - fl1_fx * tey_yyz_xy_1[j];

                tez_xyyz_xxy_0[j] = pa_x[j] * tez_yyz_xxy_0[j] - pc_x[j] * tez_yyz_xxy_1[j] + fl1_fx * tez_yyz_xy_0[j] - fl1_fx * tez_yyz_xy_1[j];

                tex_xyyz_xxz_0[j] =
                    pa_x[j] * tex_yyz_xxz_0[j] - pc_x[j] * tex_yyz_xxz_1[j] + fl1_fx * tex_yyz_xz_0[j] - fl1_fx * tex_yyz_xz_1[j] + ta_yyz_xxz_1[j];

                tey_xyyz_xxz_0[j] = pa_x[j] * tey_yyz_xxz_0[j] - pc_x[j] * tey_yyz_xxz_1[j] + fl1_fx * tey_yyz_xz_0[j] - fl1_fx * tey_yyz_xz_1[j];

                tez_xyyz_xxz_0[j] = pa_x[j] * tez_yyz_xxz_0[j] - pc_x[j] * tez_yyz_xxz_1[j] + fl1_fx * tez_yyz_xz_0[j] - fl1_fx * tez_yyz_xz_1[j];

                tex_xyyz_xyy_0[j] = pa_x[j] * tex_yyz_xyy_0[j] - pc_x[j] * tex_yyz_xyy_1[j] + 0.5 * fl1_fx * tex_yyz_yy_0[j] -
                                    0.5 * fl1_fx * tex_yyz_yy_1[j] + ta_yyz_xyy_1[j];

                tey_xyyz_xyy_0[j] =
                    pa_x[j] * tey_yyz_xyy_0[j] - pc_x[j] * tey_yyz_xyy_1[j] + 0.5 * fl1_fx * tey_yyz_yy_0[j] - 0.5 * fl1_fx * tey_yyz_yy_1[j];

                tez_xyyz_xyy_0[j] =
                    pa_x[j] * tez_yyz_xyy_0[j] - pc_x[j] * tez_yyz_xyy_1[j] + 0.5 * fl1_fx * tez_yyz_yy_0[j] - 0.5 * fl1_fx * tez_yyz_yy_1[j];

                tex_xyyz_xyz_0[j] = pa_x[j] * tex_yyz_xyz_0[j] - pc_x[j] * tex_yyz_xyz_1[j] + 0.5 * fl1_fx * tex_yyz_yz_0[j] -
                                    0.5 * fl1_fx * tex_yyz_yz_1[j] + ta_yyz_xyz_1[j];

                tey_xyyz_xyz_0[j] =
                    pa_x[j] * tey_yyz_xyz_0[j] - pc_x[j] * tey_yyz_xyz_1[j] + 0.5 * fl1_fx * tey_yyz_yz_0[j] - 0.5 * fl1_fx * tey_yyz_yz_1[j];

                tez_xyyz_xyz_0[j] =
                    pa_x[j] * tez_yyz_xyz_0[j] - pc_x[j] * tez_yyz_xyz_1[j] + 0.5 * fl1_fx * tez_yyz_yz_0[j] - 0.5 * fl1_fx * tez_yyz_yz_1[j];

                tex_xyyz_xzz_0[j] = pa_x[j] * tex_yyz_xzz_0[j] - pc_x[j] * tex_yyz_xzz_1[j] + 0.5 * fl1_fx * tex_yyz_zz_0[j] -
                                    0.5 * fl1_fx * tex_yyz_zz_1[j] + ta_yyz_xzz_1[j];

                tey_xyyz_xzz_0[j] =
                    pa_x[j] * tey_yyz_xzz_0[j] - pc_x[j] * tey_yyz_xzz_1[j] + 0.5 * fl1_fx * tey_yyz_zz_0[j] - 0.5 * fl1_fx * tey_yyz_zz_1[j];

                tez_xyyz_xzz_0[j] =
                    pa_x[j] * tez_yyz_xzz_0[j] - pc_x[j] * tez_yyz_xzz_1[j] + 0.5 * fl1_fx * tez_yyz_zz_0[j] - 0.5 * fl1_fx * tez_yyz_zz_1[j];

                tex_xyyz_yyy_0[j] = pa_x[j] * tex_yyz_yyy_0[j] - pc_x[j] * tex_yyz_yyy_1[j] + ta_yyz_yyy_1[j];

                tey_xyyz_yyy_0[j] = pa_x[j] * tey_yyz_yyy_0[j] - pc_x[j] * tey_yyz_yyy_1[j];

                tez_xyyz_yyy_0[j] = pa_x[j] * tez_yyz_yyy_0[j] - pc_x[j] * tez_yyz_yyy_1[j];

                tex_xyyz_yyz_0[j] = pa_x[j] * tex_yyz_yyz_0[j] - pc_x[j] * tex_yyz_yyz_1[j] + ta_yyz_yyz_1[j];

                tey_xyyz_yyz_0[j] = pa_x[j] * tey_yyz_yyz_0[j] - pc_x[j] * tey_yyz_yyz_1[j];

                tez_xyyz_yyz_0[j] = pa_x[j] * tez_yyz_yyz_0[j] - pc_x[j] * tez_yyz_yyz_1[j];

                tex_xyyz_yzz_0[j] = pa_x[j] * tex_yyz_yzz_0[j] - pc_x[j] * tex_yyz_yzz_1[j] + ta_yyz_yzz_1[j];

                tey_xyyz_yzz_0[j] = pa_x[j] * tey_yyz_yzz_0[j] - pc_x[j] * tey_yyz_yzz_1[j];

                tez_xyyz_yzz_0[j] = pa_x[j] * tez_yyz_yzz_0[j] - pc_x[j] * tez_yyz_yzz_1[j];

                tex_xyyz_zzz_0[j] = pa_x[j] * tex_yyz_zzz_0[j] - pc_x[j] * tex_yyz_zzz_1[j] + ta_yyz_zzz_1[j];

                tey_xyyz_zzz_0[j] = pa_x[j] * tey_yyz_zzz_0[j] - pc_x[j] * tey_yyz_zzz_1[j];

                tez_xyyz_zzz_0[j] = pa_x[j] * tez_yyz_zzz_0[j] - pc_x[j] * tez_yyz_zzz_1[j];

                tex_xyzz_xxx_0[j] = pa_x[j] * tex_yzz_xxx_0[j] - pc_x[j] * tex_yzz_xxx_1[j] + 1.5 * fl1_fx * tex_yzz_xx_0[j] -
                                    1.5 * fl1_fx * tex_yzz_xx_1[j] + ta_yzz_xxx_1[j];

                tey_xyzz_xxx_0[j] =
                    pa_x[j] * tey_yzz_xxx_0[j] - pc_x[j] * tey_yzz_xxx_1[j] + 1.5 * fl1_fx * tey_yzz_xx_0[j] - 1.5 * fl1_fx * tey_yzz_xx_1[j];

                tez_xyzz_xxx_0[j] =
                    pa_x[j] * tez_yzz_xxx_0[j] - pc_x[j] * tez_yzz_xxx_1[j] + 1.5 * fl1_fx * tez_yzz_xx_0[j] - 1.5 * fl1_fx * tez_yzz_xx_1[j];

                tex_xyzz_xxy_0[j] =
                    pa_x[j] * tex_yzz_xxy_0[j] - pc_x[j] * tex_yzz_xxy_1[j] + fl1_fx * tex_yzz_xy_0[j] - fl1_fx * tex_yzz_xy_1[j] + ta_yzz_xxy_1[j];

                tey_xyzz_xxy_0[j] = pa_x[j] * tey_yzz_xxy_0[j] - pc_x[j] * tey_yzz_xxy_1[j] + fl1_fx * tey_yzz_xy_0[j] - fl1_fx * tey_yzz_xy_1[j];

                tez_xyzz_xxy_0[j] = pa_x[j] * tez_yzz_xxy_0[j] - pc_x[j] * tez_yzz_xxy_1[j] + fl1_fx * tez_yzz_xy_0[j] - fl1_fx * tez_yzz_xy_1[j];

                tex_xyzz_xxz_0[j] =
                    pa_x[j] * tex_yzz_xxz_0[j] - pc_x[j] * tex_yzz_xxz_1[j] + fl1_fx * tex_yzz_xz_0[j] - fl1_fx * tex_yzz_xz_1[j] + ta_yzz_xxz_1[j];

                tey_xyzz_xxz_0[j] = pa_x[j] * tey_yzz_xxz_0[j] - pc_x[j] * tey_yzz_xxz_1[j] + fl1_fx * tey_yzz_xz_0[j] - fl1_fx * tey_yzz_xz_1[j];

                tez_xyzz_xxz_0[j] = pa_x[j] * tez_yzz_xxz_0[j] - pc_x[j] * tez_yzz_xxz_1[j] + fl1_fx * tez_yzz_xz_0[j] - fl1_fx * tez_yzz_xz_1[j];

                tex_xyzz_xyy_0[j] = pa_x[j] * tex_yzz_xyy_0[j] - pc_x[j] * tex_yzz_xyy_1[j] + 0.5 * fl1_fx * tex_yzz_yy_0[j] -
                                    0.5 * fl1_fx * tex_yzz_yy_1[j] + ta_yzz_xyy_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForGF_250_300(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pcDistances,
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

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_4_3_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tey_yzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 51);

            auto tez_yzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 51);

            auto tex_yzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 52);

            auto tey_yzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 52);

            auto tez_yzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 52);

            auto tex_yzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 53);

            auto tey_yzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 53);

            auto tez_yzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 53);

            auto tex_zzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 54);

            auto tey_zzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 54);

            auto tez_zzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 54);

            auto tex_zzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 55);

            auto tey_zzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 55);

            auto tez_zzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 55);

            auto tex_zzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 56);

            auto tey_zzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 56);

            auto tez_zzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 56);

            auto tex_zzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 57);

            auto tey_zzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 57);

            auto tez_zzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 57);

            auto tex_zzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 58);

            auto tey_zzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 58);

            auto tez_zzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 58);

            auto tex_zzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 59);

            auto tey_zzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 59);

            auto tez_zzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 59);

            auto tey_yzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 51);

            auto tez_yzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 51);

            auto tex_yzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 52);

            auto tey_yzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 52);

            auto tez_yzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 52);

            auto tex_yzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 53);

            auto tey_yzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 53);

            auto tez_yzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 53);

            auto tex_zzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 54);

            auto tey_zzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 54);

            auto tez_zzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 54);

            auto tex_zzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 55);

            auto tey_zzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 55);

            auto tez_zzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 55);

            auto tex_zzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 56);

            auto tey_zzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 56);

            auto tez_zzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 56);

            auto tex_zzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 57);

            auto tey_zzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 57);

            auto tez_zzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 57);

            auto tex_zzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 58);

            auto tey_zzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 58);

            auto tez_zzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 58);

            auto tex_zzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 59);

            auto tey_zzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 59);

            auto tez_zzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 59);

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

            auto tey_xyzz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 83);

            auto tez_xyzz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 83);

            auto tex_xyzz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 84);

            auto tey_xyzz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 84);

            auto tez_xyzz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 84);

            auto tex_xyzz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 85);

            auto tey_xyzz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 85);

            auto tez_xyzz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 85);

            auto tex_xyzz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 86);

            auto tey_xyzz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 86);

            auto tez_xyzz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 86);

            auto tex_xyzz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 87);

            auto tey_xyzz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 87);

            auto tez_xyzz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 87);

            auto tex_xyzz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 88);

            auto tey_xyzz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 88);

            auto tez_xyzz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 88);

            auto tex_xyzz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 89);

            auto tey_xyzz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 89);

            auto tez_xyzz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 89);

            auto tex_xzzz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 90);

            auto tey_xzzz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 90);

            auto tez_xzzz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 90);

            auto tex_xzzz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 91);

            auto tey_xzzz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 91);

            auto tez_xzzz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 91);

            auto tex_xzzz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 92);

            auto tey_xzzz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 92);

            auto tez_xzzz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 92);

            auto tex_xzzz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 93);

            auto tey_xzzz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 93);

            auto tez_xzzz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 93);

            auto tex_xzzz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 94);

            auto tey_xzzz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 94);

            auto tez_xzzz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 94);

            auto tex_xzzz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 95);

            auto tey_xzzz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 95);

            auto tez_xzzz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 95);

            auto tex_xzzz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 96);

            auto tey_xzzz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 96);

            auto tez_xzzz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 96);

            auto tex_xzzz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 97);

            auto tey_xzzz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 97);

            auto tez_xzzz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 97);

            auto tex_xzzz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 98);

            auto tey_xzzz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 98);

            auto tez_xzzz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 98);

            auto tex_xzzz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 99);

            auto tey_xzzz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 99);

            auto tez_xzzz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 99);

            // Batch of Integrals (250,300)

            #pragma omp simd aligned(fx, pa_x, pc_x, ta_yzz_xyz_1, ta_yzz_xzz_1, ta_yzz_yyy_1, ta_yzz_yyz_1, \
                                         ta_yzz_yzz_1, ta_yzz_zzz_1, ta_zzz_xxx_1, ta_zzz_xxy_1, ta_zzz_xxz_1, ta_zzz_xyy_1, \
                                         ta_zzz_xyz_1, ta_zzz_xzz_1, ta_zzz_yyy_1, ta_zzz_yyz_1, ta_zzz_yzz_1, ta_zzz_zzz_1, \
                                         tex_xyzz_xyz_0, tex_xyzz_xzz_0, tex_xyzz_yyy_0, tex_xyzz_yyz_0, tex_xyzz_yzz_0, \
                                         tex_xyzz_zzz_0, tex_xzzz_xxx_0, tex_xzzz_xxy_0, tex_xzzz_xxz_0, tex_xzzz_xyy_0, \
                                         tex_xzzz_xyz_0, tex_xzzz_xzz_0, tex_xzzz_yyy_0, tex_xzzz_yyz_0, tex_xzzz_yzz_0, \
                                         tex_xzzz_zzz_0, tex_yzz_xyz_0, tex_yzz_xyz_1, tex_yzz_xzz_0, tex_yzz_xzz_1, \
                                         tex_yzz_yyy_0, tex_yzz_yyy_1, tex_yzz_yyz_0, tex_yzz_yyz_1, tex_yzz_yz_0, \
                                         tex_yzz_yz_1, tex_yzz_yzz_0, tex_yzz_yzz_1, tex_yzz_zz_0, tex_yzz_zz_1, \
                                         tex_yzz_zzz_0, tex_yzz_zzz_1, tex_zzz_xx_0, tex_zzz_xx_1, tex_zzz_xxx_0, \
                                         tex_zzz_xxx_1, tex_zzz_xxy_0, tex_zzz_xxy_1, tex_zzz_xxz_0, tex_zzz_xxz_1, \
                                         tex_zzz_xy_0, tex_zzz_xy_1, tex_zzz_xyy_0, tex_zzz_xyy_1, tex_zzz_xyz_0, \
                                         tex_zzz_xyz_1, tex_zzz_xz_0, tex_zzz_xz_1, tex_zzz_xzz_0, tex_zzz_xzz_1, \
                                         tex_zzz_yy_0, tex_zzz_yy_1, tex_zzz_yyy_0, tex_zzz_yyy_1, tex_zzz_yyz_0, \
                                         tex_zzz_yyz_1, tex_zzz_yz_0, tex_zzz_yz_1, tex_zzz_yzz_0, tex_zzz_yzz_1, \
                                         tex_zzz_zz_0, tex_zzz_zz_1, tex_zzz_zzz_0, tex_zzz_zzz_1, tey_xyzz_xyy_0, \
                                         tey_xyzz_xyz_0, tey_xyzz_xzz_0, tey_xyzz_yyy_0, tey_xyzz_yyz_0, tey_xyzz_yzz_0, \
                                         tey_xyzz_zzz_0, tey_xzzz_xxx_0, tey_xzzz_xxy_0, tey_xzzz_xxz_0, tey_xzzz_xyy_0, \
                                         tey_xzzz_xyz_0, tey_xzzz_xzz_0, tey_xzzz_yyy_0, tey_xzzz_yyz_0, tey_xzzz_yzz_0, \
                                         tey_xzzz_zzz_0, tey_yzz_xyy_0, tey_yzz_xyy_1, tey_yzz_xyz_0, tey_yzz_xyz_1, \
                                         tey_yzz_xzz_0, tey_yzz_xzz_1, tey_yzz_yy_0, tey_yzz_yy_1, tey_yzz_yyy_0, \
                                         tey_yzz_yyy_1, tey_yzz_yyz_0, tey_yzz_yyz_1, tey_yzz_yz_0, tey_yzz_yz_1, \
                                         tey_yzz_yzz_0, tey_yzz_yzz_1, tey_yzz_zz_0, tey_yzz_zz_1, tey_yzz_zzz_0, \
                                         tey_yzz_zzz_1, tey_zzz_xx_0, tey_zzz_xx_1, tey_zzz_xxx_0, tey_zzz_xxx_1, \
                                         tey_zzz_xxy_0, tey_zzz_xxy_1, tey_zzz_xxz_0, tey_zzz_xxz_1, tey_zzz_xy_0, \
                                         tey_zzz_xy_1, tey_zzz_xyy_0, tey_zzz_xyy_1, tey_zzz_xyz_0, tey_zzz_xyz_1, \
                                         tey_zzz_xz_0, tey_zzz_xz_1, tey_zzz_xzz_0, tey_zzz_xzz_1, tey_zzz_yy_0, \
                                         tey_zzz_yy_1, tey_zzz_yyy_0, tey_zzz_yyy_1, tey_zzz_yyz_0, tey_zzz_yyz_1, \
                                         tey_zzz_yz_0, tey_zzz_yz_1, tey_zzz_yzz_0, tey_zzz_yzz_1, tey_zzz_zz_0, \
                                         tey_zzz_zz_1, tey_zzz_zzz_0, tey_zzz_zzz_1, tez_xyzz_xyy_0, tez_xyzz_xyz_0, \
                                         tez_xyzz_xzz_0, tez_xyzz_yyy_0, tez_xyzz_yyz_0, tez_xyzz_yzz_0, tez_xyzz_zzz_0, \
                                         tez_xzzz_xxx_0, tez_xzzz_xxy_0, tez_xzzz_xxz_0, tez_xzzz_xyy_0, tez_xzzz_xyz_0, \
                                         tez_xzzz_xzz_0, tez_xzzz_yyy_0, tez_xzzz_yyz_0, tez_xzzz_yzz_0, tez_xzzz_zzz_0, \
                                         tez_yzz_xyy_0, tez_yzz_xyy_1, tez_yzz_xyz_0, tez_yzz_xyz_1, tez_yzz_xzz_0, \
                                         tez_yzz_xzz_1, tez_yzz_yy_0, tez_yzz_yy_1, tez_yzz_yyy_0, tez_yzz_yyy_1, \
                                         tez_yzz_yyz_0, tez_yzz_yyz_1, tez_yzz_yz_0, tez_yzz_yz_1, tez_yzz_yzz_0, \
                                         tez_yzz_yzz_1, tez_yzz_zz_0, tez_yzz_zz_1, tez_yzz_zzz_0, tez_yzz_zzz_1, \
                                         tez_zzz_xx_0, tez_zzz_xx_1, tez_zzz_xxx_0, tez_zzz_xxx_1, tez_zzz_xxy_0, \
                                         tez_zzz_xxy_1, tez_zzz_xxz_0, tez_zzz_xxz_1, tez_zzz_xy_0, tez_zzz_xy_1, \
                                         tez_zzz_xyy_0, tez_zzz_xyy_1, tez_zzz_xyz_0, tez_zzz_xyz_1, tez_zzz_xz_0, \
                                         tez_zzz_xz_1, tez_zzz_xzz_0, tez_zzz_xzz_1, tez_zzz_yy_0, tez_zzz_yy_1, \
                                         tez_zzz_yyy_0, tez_zzz_yyy_1, tez_zzz_yyz_0, tez_zzz_yyz_1, tez_zzz_yz_0, \
                                         tez_zzz_yz_1, tez_zzz_yzz_0, tez_zzz_yzz_1, tez_zzz_zz_0, tez_zzz_zz_1, \
                                         tez_zzz_zzz_0, tez_zzz_zzz_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tey_xyzz_xyy_0[j] =
                    pa_x[j] * tey_yzz_xyy_0[j] - pc_x[j] * tey_yzz_xyy_1[j] + 0.5 * fl1_fx * tey_yzz_yy_0[j] - 0.5 * fl1_fx * tey_yzz_yy_1[j];

                tez_xyzz_xyy_0[j] =
                    pa_x[j] * tez_yzz_xyy_0[j] - pc_x[j] * tez_yzz_xyy_1[j] + 0.5 * fl1_fx * tez_yzz_yy_0[j] - 0.5 * fl1_fx * tez_yzz_yy_1[j];

                tex_xyzz_xyz_0[j] = pa_x[j] * tex_yzz_xyz_0[j] - pc_x[j] * tex_yzz_xyz_1[j] + 0.5 * fl1_fx * tex_yzz_yz_0[j] -
                                    0.5 * fl1_fx * tex_yzz_yz_1[j] + ta_yzz_xyz_1[j];

                tey_xyzz_xyz_0[j] =
                    pa_x[j] * tey_yzz_xyz_0[j] - pc_x[j] * tey_yzz_xyz_1[j] + 0.5 * fl1_fx * tey_yzz_yz_0[j] - 0.5 * fl1_fx * tey_yzz_yz_1[j];

                tez_xyzz_xyz_0[j] =
                    pa_x[j] * tez_yzz_xyz_0[j] - pc_x[j] * tez_yzz_xyz_1[j] + 0.5 * fl1_fx * tez_yzz_yz_0[j] - 0.5 * fl1_fx * tez_yzz_yz_1[j];

                tex_xyzz_xzz_0[j] = pa_x[j] * tex_yzz_xzz_0[j] - pc_x[j] * tex_yzz_xzz_1[j] + 0.5 * fl1_fx * tex_yzz_zz_0[j] -
                                    0.5 * fl1_fx * tex_yzz_zz_1[j] + ta_yzz_xzz_1[j];

                tey_xyzz_xzz_0[j] =
                    pa_x[j] * tey_yzz_xzz_0[j] - pc_x[j] * tey_yzz_xzz_1[j] + 0.5 * fl1_fx * tey_yzz_zz_0[j] - 0.5 * fl1_fx * tey_yzz_zz_1[j];

                tez_xyzz_xzz_0[j] =
                    pa_x[j] * tez_yzz_xzz_0[j] - pc_x[j] * tez_yzz_xzz_1[j] + 0.5 * fl1_fx * tez_yzz_zz_0[j] - 0.5 * fl1_fx * tez_yzz_zz_1[j];

                tex_xyzz_yyy_0[j] = pa_x[j] * tex_yzz_yyy_0[j] - pc_x[j] * tex_yzz_yyy_1[j] + ta_yzz_yyy_1[j];

                tey_xyzz_yyy_0[j] = pa_x[j] * tey_yzz_yyy_0[j] - pc_x[j] * tey_yzz_yyy_1[j];

                tez_xyzz_yyy_0[j] = pa_x[j] * tez_yzz_yyy_0[j] - pc_x[j] * tez_yzz_yyy_1[j];

                tex_xyzz_yyz_0[j] = pa_x[j] * tex_yzz_yyz_0[j] - pc_x[j] * tex_yzz_yyz_1[j] + ta_yzz_yyz_1[j];

                tey_xyzz_yyz_0[j] = pa_x[j] * tey_yzz_yyz_0[j] - pc_x[j] * tey_yzz_yyz_1[j];

                tez_xyzz_yyz_0[j] = pa_x[j] * tez_yzz_yyz_0[j] - pc_x[j] * tez_yzz_yyz_1[j];

                tex_xyzz_yzz_0[j] = pa_x[j] * tex_yzz_yzz_0[j] - pc_x[j] * tex_yzz_yzz_1[j] + ta_yzz_yzz_1[j];

                tey_xyzz_yzz_0[j] = pa_x[j] * tey_yzz_yzz_0[j] - pc_x[j] * tey_yzz_yzz_1[j];

                tez_xyzz_yzz_0[j] = pa_x[j] * tez_yzz_yzz_0[j] - pc_x[j] * tez_yzz_yzz_1[j];

                tex_xyzz_zzz_0[j] = pa_x[j] * tex_yzz_zzz_0[j] - pc_x[j] * tex_yzz_zzz_1[j] + ta_yzz_zzz_1[j];

                tey_xyzz_zzz_0[j] = pa_x[j] * tey_yzz_zzz_0[j] - pc_x[j] * tey_yzz_zzz_1[j];

                tez_xyzz_zzz_0[j] = pa_x[j] * tez_yzz_zzz_0[j] - pc_x[j] * tez_yzz_zzz_1[j];

                tex_xzzz_xxx_0[j] = pa_x[j] * tex_zzz_xxx_0[j] - pc_x[j] * tex_zzz_xxx_1[j] + 1.5 * fl1_fx * tex_zzz_xx_0[j] -
                                    1.5 * fl1_fx * tex_zzz_xx_1[j] + ta_zzz_xxx_1[j];

                tey_xzzz_xxx_0[j] =
                    pa_x[j] * tey_zzz_xxx_0[j] - pc_x[j] * tey_zzz_xxx_1[j] + 1.5 * fl1_fx * tey_zzz_xx_0[j] - 1.5 * fl1_fx * tey_zzz_xx_1[j];

                tez_xzzz_xxx_0[j] =
                    pa_x[j] * tez_zzz_xxx_0[j] - pc_x[j] * tez_zzz_xxx_1[j] + 1.5 * fl1_fx * tez_zzz_xx_0[j] - 1.5 * fl1_fx * tez_zzz_xx_1[j];

                tex_xzzz_xxy_0[j] =
                    pa_x[j] * tex_zzz_xxy_0[j] - pc_x[j] * tex_zzz_xxy_1[j] + fl1_fx * tex_zzz_xy_0[j] - fl1_fx * tex_zzz_xy_1[j] + ta_zzz_xxy_1[j];

                tey_xzzz_xxy_0[j] = pa_x[j] * tey_zzz_xxy_0[j] - pc_x[j] * tey_zzz_xxy_1[j] + fl1_fx * tey_zzz_xy_0[j] - fl1_fx * tey_zzz_xy_1[j];

                tez_xzzz_xxy_0[j] = pa_x[j] * tez_zzz_xxy_0[j] - pc_x[j] * tez_zzz_xxy_1[j] + fl1_fx * tez_zzz_xy_0[j] - fl1_fx * tez_zzz_xy_1[j];

                tex_xzzz_xxz_0[j] =
                    pa_x[j] * tex_zzz_xxz_0[j] - pc_x[j] * tex_zzz_xxz_1[j] + fl1_fx * tex_zzz_xz_0[j] - fl1_fx * tex_zzz_xz_1[j] + ta_zzz_xxz_1[j];

                tey_xzzz_xxz_0[j] = pa_x[j] * tey_zzz_xxz_0[j] - pc_x[j] * tey_zzz_xxz_1[j] + fl1_fx * tey_zzz_xz_0[j] - fl1_fx * tey_zzz_xz_1[j];

                tez_xzzz_xxz_0[j] = pa_x[j] * tez_zzz_xxz_0[j] - pc_x[j] * tez_zzz_xxz_1[j] + fl1_fx * tez_zzz_xz_0[j] - fl1_fx * tez_zzz_xz_1[j];

                tex_xzzz_xyy_0[j] = pa_x[j] * tex_zzz_xyy_0[j] - pc_x[j] * tex_zzz_xyy_1[j] + 0.5 * fl1_fx * tex_zzz_yy_0[j] -
                                    0.5 * fl1_fx * tex_zzz_yy_1[j] + ta_zzz_xyy_1[j];

                tey_xzzz_xyy_0[j] =
                    pa_x[j] * tey_zzz_xyy_0[j] - pc_x[j] * tey_zzz_xyy_1[j] + 0.5 * fl1_fx * tey_zzz_yy_0[j] - 0.5 * fl1_fx * tey_zzz_yy_1[j];

                tez_xzzz_xyy_0[j] =
                    pa_x[j] * tez_zzz_xyy_0[j] - pc_x[j] * tez_zzz_xyy_1[j] + 0.5 * fl1_fx * tez_zzz_yy_0[j] - 0.5 * fl1_fx * tez_zzz_yy_1[j];

                tex_xzzz_xyz_0[j] = pa_x[j] * tex_zzz_xyz_0[j] - pc_x[j] * tex_zzz_xyz_1[j] + 0.5 * fl1_fx * tex_zzz_yz_0[j] -
                                    0.5 * fl1_fx * tex_zzz_yz_1[j] + ta_zzz_xyz_1[j];

                tey_xzzz_xyz_0[j] =
                    pa_x[j] * tey_zzz_xyz_0[j] - pc_x[j] * tey_zzz_xyz_1[j] + 0.5 * fl1_fx * tey_zzz_yz_0[j] - 0.5 * fl1_fx * tey_zzz_yz_1[j];

                tez_xzzz_xyz_0[j] =
                    pa_x[j] * tez_zzz_xyz_0[j] - pc_x[j] * tez_zzz_xyz_1[j] + 0.5 * fl1_fx * tez_zzz_yz_0[j] - 0.5 * fl1_fx * tez_zzz_yz_1[j];

                tex_xzzz_xzz_0[j] = pa_x[j] * tex_zzz_xzz_0[j] - pc_x[j] * tex_zzz_xzz_1[j] + 0.5 * fl1_fx * tex_zzz_zz_0[j] -
                                    0.5 * fl1_fx * tex_zzz_zz_1[j] + ta_zzz_xzz_1[j];

                tey_xzzz_xzz_0[j] =
                    pa_x[j] * tey_zzz_xzz_0[j] - pc_x[j] * tey_zzz_xzz_1[j] + 0.5 * fl1_fx * tey_zzz_zz_0[j] - 0.5 * fl1_fx * tey_zzz_zz_1[j];

                tez_xzzz_xzz_0[j] =
                    pa_x[j] * tez_zzz_xzz_0[j] - pc_x[j] * tez_zzz_xzz_1[j] + 0.5 * fl1_fx * tez_zzz_zz_0[j] - 0.5 * fl1_fx * tez_zzz_zz_1[j];

                tex_xzzz_yyy_0[j] = pa_x[j] * tex_zzz_yyy_0[j] - pc_x[j] * tex_zzz_yyy_1[j] + ta_zzz_yyy_1[j];

                tey_xzzz_yyy_0[j] = pa_x[j] * tey_zzz_yyy_0[j] - pc_x[j] * tey_zzz_yyy_1[j];

                tez_xzzz_yyy_0[j] = pa_x[j] * tez_zzz_yyy_0[j] - pc_x[j] * tez_zzz_yyy_1[j];

                tex_xzzz_yyz_0[j] = pa_x[j] * tex_zzz_yyz_0[j] - pc_x[j] * tex_zzz_yyz_1[j] + ta_zzz_yyz_1[j];

                tey_xzzz_yyz_0[j] = pa_x[j] * tey_zzz_yyz_0[j] - pc_x[j] * tey_zzz_yyz_1[j];

                tez_xzzz_yyz_0[j] = pa_x[j] * tez_zzz_yyz_0[j] - pc_x[j] * tez_zzz_yyz_1[j];

                tex_xzzz_yzz_0[j] = pa_x[j] * tex_zzz_yzz_0[j] - pc_x[j] * tex_zzz_yzz_1[j] + ta_zzz_yzz_1[j];

                tey_xzzz_yzz_0[j] = pa_x[j] * tey_zzz_yzz_0[j] - pc_x[j] * tey_zzz_yzz_1[j];

                tez_xzzz_yzz_0[j] = pa_x[j] * tez_zzz_yzz_0[j] - pc_x[j] * tez_zzz_yzz_1[j];

                tex_xzzz_zzz_0[j] = pa_x[j] * tex_zzz_zzz_0[j] - pc_x[j] * tex_zzz_zzz_1[j] + ta_zzz_zzz_1[j];

                tey_xzzz_zzz_0[j] = pa_x[j] * tey_zzz_zzz_0[j] - pc_x[j] * tey_zzz_zzz_1[j];

                tez_xzzz_zzz_0[j] = pa_x[j] * tez_zzz_zzz_0[j] - pc_x[j] * tez_zzz_zzz_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForGF_300_350(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pcDistances,
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

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_4_3_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_2_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tex_yy_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 30);

            auto tey_yy_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 30);

            auto tez_yy_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 30);

            auto tex_yy_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 31);

            auto tey_yy_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 31);

            auto tez_yy_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 31);

            auto tex_yy_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 32);

            auto tey_yy_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 32);

            auto tez_yy_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 32);

            auto tex_yy_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 33);

            auto tey_yy_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 33);

            auto tez_yy_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 33);

            auto tex_yy_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 34);

            auto tey_yy_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 34);

            auto tez_yy_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 34);

            auto tex_yy_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 35);

            auto tey_yy_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 35);

            auto tez_yy_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 35);

            auto tex_yy_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 36);

            auto tey_yy_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 36);

            auto tez_yy_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 36);

            auto tex_yy_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 37);

            auto tey_yy_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 37);

            auto tez_yy_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 37);

            auto tex_yy_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 38);

            auto tey_yy_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 38);

            auto tez_yy_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 38);

            auto tex_yy_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 39);

            auto tey_yy_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 39);

            auto tez_yy_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 39);

            auto tex_yz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 40);

            auto tey_yz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 40);

            auto tez_yz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 40);

            auto tex_yz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 41);

            auto tey_yz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 41);

            auto tez_yz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 41);

            auto tex_yz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 42);

            auto tey_yz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 42);

            auto tez_yz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 42);

            auto tex_yz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 43);

            auto tey_yz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 43);

            auto tez_yz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 43);

            auto tex_yz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 44);

            auto tey_yz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 44);

            auto tez_yz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 44);

            auto tex_yz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 45);

            auto tey_yz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 45);

            auto tez_yz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 45);

            auto tex_yz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 46);

            auto tey_yz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 46);

            auto tex_yy_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 30);

            auto tey_yy_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 30);

            auto tez_yy_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 30);

            auto tex_yy_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 31);

            auto tey_yy_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 31);

            auto tez_yy_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 31);

            auto tex_yy_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 32);

            auto tey_yy_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 32);

            auto tez_yy_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 32);

            auto tex_yy_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 33);

            auto tey_yy_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 33);

            auto tez_yy_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 33);

            auto tex_yy_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 34);

            auto tey_yy_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 34);

            auto tez_yy_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 34);

            auto tex_yy_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 35);

            auto tey_yy_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 35);

            auto tez_yy_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 35);

            auto tex_yy_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 36);

            auto tey_yy_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 36);

            auto tez_yy_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 36);

            auto tex_yy_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 37);

            auto tey_yy_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 37);

            auto tez_yy_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 37);

            auto tex_yy_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 38);

            auto tey_yy_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 38);

            auto tez_yy_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 38);

            auto tex_yy_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 39);

            auto tey_yy_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 39);

            auto tez_yy_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 39);

            auto tex_yz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 40);

            auto tey_yz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 40);

            auto tez_yz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 40);

            auto tex_yz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 41);

            auto tey_yz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 41);

            auto tez_yz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 41);

            auto tex_yz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 42);

            auto tey_yz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 42);

            auto tez_yz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 42);

            auto tex_yz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 43);

            auto tey_yz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 43);

            auto tez_yz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 43);

            auto tex_yz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 44);

            auto tey_yz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 44);

            auto tez_yz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 44);

            auto tex_yz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 45);

            auto tey_yz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 45);

            auto tez_yz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 45);

            auto tex_yz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 46);

            auto tey_yz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 46);

            auto tex_yyy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 36);

            auto tey_yyy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 36);

            auto tez_yyy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 36);

            auto tex_yyy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 37);

            auto tey_yyy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 37);

            auto tez_yyy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 37);

            auto tex_yyy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 38);

            auto tey_yyy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 38);

            auto tez_yyy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 38);

            auto tex_yyy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 39);

            auto tey_yyy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 39);

            auto tez_yyy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 39);

            auto tex_yyy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 40);

            auto tey_yyy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 40);

            auto tez_yyy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 40);

            auto tex_yyy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 41);

            auto tey_yyy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 41);

            auto tez_yyy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 41);

            auto tex_yyz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 42);

            auto tey_yyz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 42);

            auto tez_yyz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 42);

            auto tex_yyz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 43);

            auto tey_yyz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 43);

            auto tez_yyz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 43);

            auto tex_yyz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 44);

            auto tey_yyz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 44);

            auto tez_yyz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 44);

            auto tex_yyz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 45);

            auto tey_yyz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 45);

            auto tex_yyy_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 36);

            auto tey_yyy_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 36);

            auto tez_yyy_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 36);

            auto tex_yyy_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 37);

            auto tey_yyy_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 37);

            auto tez_yyy_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 37);

            auto tex_yyy_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 38);

            auto tey_yyy_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 38);

            auto tez_yyy_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 38);

            auto tex_yyy_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 39);

            auto tey_yyy_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 39);

            auto tez_yyy_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 39);

            auto tex_yyy_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 40);

            auto tey_yyy_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 40);

            auto tez_yyy_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 40);

            auto tex_yyy_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 41);

            auto tey_yyy_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 41);

            auto tez_yyy_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 41);

            auto tex_yyz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 42);

            auto tey_yyz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 42);

            auto tez_yyz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 42);

            auto tex_yyz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 43);

            auto tey_yyz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 43);

            auto tez_yyz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 43);

            auto tex_yyz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 44);

            auto tey_yyz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 44);

            auto tez_yyz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 44);

            auto tex_yyz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 45);

            auto tey_yyz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 45);

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

            // set up pointers to integrals

            auto tex_yyyy_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 100);

            auto tey_yyyy_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 100);

            auto tez_yyyy_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 100);

            auto tex_yyyy_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 101);

            auto tey_yyyy_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 101);

            auto tez_yyyy_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 101);

            auto tex_yyyy_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 102);

            auto tey_yyyy_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 102);

            auto tez_yyyy_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 102);

            auto tex_yyyy_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 103);

            auto tey_yyyy_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 103);

            auto tez_yyyy_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 103);

            auto tex_yyyy_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 104);

            auto tey_yyyy_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 104);

            auto tez_yyyy_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 104);

            auto tex_yyyy_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 105);

            auto tey_yyyy_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 105);

            auto tez_yyyy_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 105);

            auto tex_yyyy_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 106);

            auto tey_yyyy_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 106);

            auto tez_yyyy_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 106);

            auto tex_yyyy_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 107);

            auto tey_yyyy_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 107);

            auto tez_yyyy_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 107);

            auto tex_yyyy_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 108);

            auto tey_yyyy_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 108);

            auto tez_yyyy_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 108);

            auto tex_yyyy_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 109);

            auto tey_yyyy_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 109);

            auto tez_yyyy_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 109);

            auto tex_yyyz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 110);

            auto tey_yyyz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 110);

            auto tez_yyyz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 110);

            auto tex_yyyz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 111);

            auto tey_yyyz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 111);

            auto tez_yyyz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 111);

            auto tex_yyyz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 112);

            auto tey_yyyz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 112);

            auto tez_yyyz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 112);

            auto tex_yyyz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 113);

            auto tey_yyyz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 113);

            auto tez_yyyz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 113);

            auto tex_yyyz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 114);

            auto tey_yyyz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 114);

            auto tez_yyyz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 114);

            auto tex_yyyz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 115);

            auto tey_yyyz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 115);

            auto tez_yyyz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 115);

            auto tex_yyyz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 116);

            auto tey_yyyz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 116);

            // Batch of Integrals (300,350)

            #pragma omp simd aligned(fx, pa_y, pc_y, ta_yyy_xxx_1, ta_yyy_xxy_1, ta_yyy_xxz_1, ta_yyy_xyy_1, \
                                         ta_yyy_xyz_1, ta_yyy_xzz_1, ta_yyy_yyy_1, ta_yyy_yyz_1, ta_yyy_yzz_1, ta_yyy_zzz_1, \
                                         ta_yyz_xxx_1, ta_yyz_xxy_1, ta_yyz_xxz_1, ta_yyz_xyy_1, ta_yyz_xyz_1, ta_yyz_xzz_1, \
                                         ta_yyz_yyy_1, tex_yy_xxx_0, tex_yy_xxx_1, tex_yy_xxy_0, tex_yy_xxy_1, tex_yy_xxz_0, \
                                         tex_yy_xxz_1, tex_yy_xyy_0, tex_yy_xyy_1, tex_yy_xyz_0, tex_yy_xyz_1, tex_yy_xzz_0, \
                                         tex_yy_xzz_1, tex_yy_yyy_0, tex_yy_yyy_1, tex_yy_yyz_0, tex_yy_yyz_1, tex_yy_yzz_0, \
                                         tex_yy_yzz_1, tex_yy_zzz_0, tex_yy_zzz_1, tex_yyy_xx_0, tex_yyy_xx_1, \
                                         tex_yyy_xxx_0, tex_yyy_xxx_1, tex_yyy_xxy_0, tex_yyy_xxy_1, tex_yyy_xxz_0, \
                                         tex_yyy_xxz_1, tex_yyy_xy_0, tex_yyy_xy_1, tex_yyy_xyy_0, tex_yyy_xyy_1, \
                                         tex_yyy_xyz_0, tex_yyy_xyz_1, tex_yyy_xz_0, tex_yyy_xz_1, tex_yyy_xzz_0, \
                                         tex_yyy_xzz_1, tex_yyy_yy_0, tex_yyy_yy_1, tex_yyy_yyy_0, tex_yyy_yyy_1, \
                                         tex_yyy_yyz_0, tex_yyy_yyz_1, tex_yyy_yz_0, tex_yyy_yz_1, tex_yyy_yzz_0, \
                                         tex_yyy_yzz_1, tex_yyy_zz_0, tex_yyy_zz_1, tex_yyy_zzz_0, tex_yyy_zzz_1, \
                                         tex_yyyy_xxx_0, tex_yyyy_xxy_0, tex_yyyy_xxz_0, tex_yyyy_xyy_0, tex_yyyy_xyz_0, \
                                         tex_yyyy_xzz_0, tex_yyyy_yyy_0, tex_yyyy_yyz_0, tex_yyyy_yzz_0, tex_yyyy_zzz_0, \
                                         tex_yyyz_xxx_0, tex_yyyz_xxy_0, tex_yyyz_xxz_0, tex_yyyz_xyy_0, tex_yyyz_xyz_0, \
                                         tex_yyyz_xzz_0, tex_yyyz_yyy_0, tex_yyz_xx_0, tex_yyz_xx_1, tex_yyz_xxx_0, \
                                         tex_yyz_xxx_1, tex_yyz_xxy_0, tex_yyz_xxy_1, tex_yyz_xxz_0, tex_yyz_xxz_1, \
                                         tex_yyz_xy_0, tex_yyz_xy_1, tex_yyz_xyy_0, tex_yyz_xyy_1, tex_yyz_xyz_0, \
                                         tex_yyz_xyz_1, tex_yyz_xz_0, tex_yyz_xz_1, tex_yyz_xzz_0, tex_yyz_xzz_1, \
                                         tex_yyz_yy_0, tex_yyz_yy_1, tex_yyz_yyy_0, tex_yyz_yyy_1, tex_yz_xxx_0, \
                                         tex_yz_xxx_1, tex_yz_xxy_0, tex_yz_xxy_1, tex_yz_xxz_0, tex_yz_xxz_1, tex_yz_xyy_0, \
                                         tex_yz_xyy_1, tex_yz_xyz_0, tex_yz_xyz_1, tex_yz_xzz_0, tex_yz_xzz_1, tex_yz_yyy_0, \
                                         tex_yz_yyy_1, tey_yy_xxx_0, tey_yy_xxx_1, tey_yy_xxy_0, tey_yy_xxy_1, tey_yy_xxz_0, \
                                         tey_yy_xxz_1, tey_yy_xyy_0, tey_yy_xyy_1, tey_yy_xyz_0, tey_yy_xyz_1, tey_yy_xzz_0, \
                                         tey_yy_xzz_1, tey_yy_yyy_0, tey_yy_yyy_1, tey_yy_yyz_0, tey_yy_yyz_1, tey_yy_yzz_0, \
                                         tey_yy_yzz_1, tey_yy_zzz_0, tey_yy_zzz_1, tey_yyy_xx_0, tey_yyy_xx_1, \
                                         tey_yyy_xxx_0, tey_yyy_xxx_1, tey_yyy_xxy_0, tey_yyy_xxy_1, tey_yyy_xxz_0, \
                                         tey_yyy_xxz_1, tey_yyy_xy_0, tey_yyy_xy_1, tey_yyy_xyy_0, tey_yyy_xyy_1, \
                                         tey_yyy_xyz_0, tey_yyy_xyz_1, tey_yyy_xz_0, tey_yyy_xz_1, tey_yyy_xzz_0, \
                                         tey_yyy_xzz_1, tey_yyy_yy_0, tey_yyy_yy_1, tey_yyy_yyy_0, tey_yyy_yyy_1, \
                                         tey_yyy_yyz_0, tey_yyy_yyz_1, tey_yyy_yz_0, tey_yyy_yz_1, tey_yyy_yzz_0, \
                                         tey_yyy_yzz_1, tey_yyy_zz_0, tey_yyy_zz_1, tey_yyy_zzz_0, tey_yyy_zzz_1, \
                                         tey_yyyy_xxx_0, tey_yyyy_xxy_0, tey_yyyy_xxz_0, tey_yyyy_xyy_0, tey_yyyy_xyz_0, \
                                         tey_yyyy_xzz_0, tey_yyyy_yyy_0, tey_yyyy_yyz_0, tey_yyyy_yzz_0, tey_yyyy_zzz_0, \
                                         tey_yyyz_xxx_0, tey_yyyz_xxy_0, tey_yyyz_xxz_0, tey_yyyz_xyy_0, tey_yyyz_xyz_0, \
                                         tey_yyyz_xzz_0, tey_yyyz_yyy_0, tey_yyz_xx_0, tey_yyz_xx_1, tey_yyz_xxx_0, \
                                         tey_yyz_xxx_1, tey_yyz_xxy_0, tey_yyz_xxy_1, tey_yyz_xxz_0, tey_yyz_xxz_1, \
                                         tey_yyz_xy_0, tey_yyz_xy_1, tey_yyz_xyy_0, tey_yyz_xyy_1, tey_yyz_xyz_0, \
                                         tey_yyz_xyz_1, tey_yyz_xz_0, tey_yyz_xz_1, tey_yyz_xzz_0, tey_yyz_xzz_1, \
                                         tey_yyz_yy_0, tey_yyz_yy_1, tey_yyz_yyy_0, tey_yyz_yyy_1, tey_yz_xxx_0, \
                                         tey_yz_xxx_1, tey_yz_xxy_0, tey_yz_xxy_1, tey_yz_xxz_0, tey_yz_xxz_1, tey_yz_xyy_0, \
                                         tey_yz_xyy_1, tey_yz_xyz_0, tey_yz_xyz_1, tey_yz_xzz_0, tey_yz_xzz_1, tey_yz_yyy_0, \
                                         tey_yz_yyy_1, tez_yy_xxx_0, tez_yy_xxx_1, tez_yy_xxy_0, tez_yy_xxy_1, tez_yy_xxz_0, \
                                         tez_yy_xxz_1, tez_yy_xyy_0, tez_yy_xyy_1, tez_yy_xyz_0, tez_yy_xyz_1, tez_yy_xzz_0, \
                                         tez_yy_xzz_1, tez_yy_yyy_0, tez_yy_yyy_1, tez_yy_yyz_0, tez_yy_yyz_1, tez_yy_yzz_0, \
                                         tez_yy_yzz_1, tez_yy_zzz_0, tez_yy_zzz_1, tez_yyy_xx_0, tez_yyy_xx_1, \
                                         tez_yyy_xxx_0, tez_yyy_xxx_1, tez_yyy_xxy_0, tez_yyy_xxy_1, tez_yyy_xxz_0, \
                                         tez_yyy_xxz_1, tez_yyy_xy_0, tez_yyy_xy_1, tez_yyy_xyy_0, tez_yyy_xyy_1, \
                                         tez_yyy_xyz_0, tez_yyy_xyz_1, tez_yyy_xz_0, tez_yyy_xz_1, tez_yyy_xzz_0, \
                                         tez_yyy_xzz_1, tez_yyy_yy_0, tez_yyy_yy_1, tez_yyy_yyy_0, tez_yyy_yyy_1, \
                                         tez_yyy_yyz_0, tez_yyy_yyz_1, tez_yyy_yz_0, tez_yyy_yz_1, tez_yyy_yzz_0, \
                                         tez_yyy_yzz_1, tez_yyy_zz_0, tez_yyy_zz_1, tez_yyy_zzz_0, tez_yyy_zzz_1, \
                                         tez_yyyy_xxx_0, tez_yyyy_xxy_0, tez_yyyy_xxz_0, tez_yyyy_xyy_0, tez_yyyy_xyz_0, \
                                         tez_yyyy_xzz_0, tez_yyyy_yyy_0, tez_yyyy_yyz_0, tez_yyyy_yzz_0, tez_yyyy_zzz_0, \
                                         tez_yyyz_xxx_0, tez_yyyz_xxy_0, tez_yyyz_xxz_0, tez_yyyz_xyy_0, tez_yyyz_xyz_0, \
                                         tez_yyyz_xzz_0, tez_yyz_xx_0, tez_yyz_xx_1, tez_yyz_xxx_0, tez_yyz_xxx_1, \
                                         tez_yyz_xxy_0, tez_yyz_xxy_1, tez_yyz_xxz_0, tez_yyz_xxz_1, tez_yyz_xy_0, \
                                         tez_yyz_xy_1, tez_yyz_xyy_0, tez_yyz_xyy_1, tez_yyz_xyz_0, tez_yyz_xyz_1, \
                                         tez_yyz_xz_0, tez_yyz_xz_1, tez_yyz_xzz_0, tez_yyz_xzz_1, tez_yz_xxx_0, \
                                         tez_yz_xxx_1, tez_yz_xxy_0, tez_yz_xxy_1, tez_yz_xxz_0, tez_yz_xxz_1, tez_yz_xyy_0, \
                                         tez_yz_xyy_1, tez_yz_xyz_0, tez_yz_xyz_1, tez_yz_xzz_0, tez_yz_xzz_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_yyyy_xxx_0[j] =
                    pa_y[j] * tex_yyy_xxx_0[j] - pc_y[j] * tex_yyy_xxx_1[j] + 1.5 * fl1_fx * tex_yy_xxx_0[j] - 1.5 * fl1_fx * tex_yy_xxx_1[j];

                tey_yyyy_xxx_0[j] = pa_y[j] * tey_yyy_xxx_0[j] - pc_y[j] * tey_yyy_xxx_1[j] + 1.5 * fl1_fx * tey_yy_xxx_0[j] -
                                    1.5 * fl1_fx * tey_yy_xxx_1[j] + ta_yyy_xxx_1[j];

                tez_yyyy_xxx_0[j] =
                    pa_y[j] * tez_yyy_xxx_0[j] - pc_y[j] * tez_yyy_xxx_1[j] + 1.5 * fl1_fx * tez_yy_xxx_0[j] - 1.5 * fl1_fx * tez_yy_xxx_1[j];

                tex_yyyy_xxy_0[j] = pa_y[j] * tex_yyy_xxy_0[j] - pc_y[j] * tex_yyy_xxy_1[j] + 1.5 * fl1_fx * tex_yy_xxy_0[j] -
                                    1.5 * fl1_fx * tex_yy_xxy_1[j] + 0.5 * fl1_fx * tex_yyy_xx_0[j] - 0.5 * fl1_fx * tex_yyy_xx_1[j];

                tey_yyyy_xxy_0[j] = pa_y[j] * tey_yyy_xxy_0[j] - pc_y[j] * tey_yyy_xxy_1[j] + 1.5 * fl1_fx * tey_yy_xxy_0[j] -
                                    1.5 * fl1_fx * tey_yy_xxy_1[j] + 0.5 * fl1_fx * tey_yyy_xx_0[j] - 0.5 * fl1_fx * tey_yyy_xx_1[j] +
                                    ta_yyy_xxy_1[j];

                tez_yyyy_xxy_0[j] = pa_y[j] * tez_yyy_xxy_0[j] - pc_y[j] * tez_yyy_xxy_1[j] + 1.5 * fl1_fx * tez_yy_xxy_0[j] -
                                    1.5 * fl1_fx * tez_yy_xxy_1[j] + 0.5 * fl1_fx * tez_yyy_xx_0[j] - 0.5 * fl1_fx * tez_yyy_xx_1[j];

                tex_yyyy_xxz_0[j] =
                    pa_y[j] * tex_yyy_xxz_0[j] - pc_y[j] * tex_yyy_xxz_1[j] + 1.5 * fl1_fx * tex_yy_xxz_0[j] - 1.5 * fl1_fx * tex_yy_xxz_1[j];

                tey_yyyy_xxz_0[j] = pa_y[j] * tey_yyy_xxz_0[j] - pc_y[j] * tey_yyy_xxz_1[j] + 1.5 * fl1_fx * tey_yy_xxz_0[j] -
                                    1.5 * fl1_fx * tey_yy_xxz_1[j] + ta_yyy_xxz_1[j];

                tez_yyyy_xxz_0[j] =
                    pa_y[j] * tez_yyy_xxz_0[j] - pc_y[j] * tez_yyy_xxz_1[j] + 1.5 * fl1_fx * tez_yy_xxz_0[j] - 1.5 * fl1_fx * tez_yy_xxz_1[j];

                tex_yyyy_xyy_0[j] = pa_y[j] * tex_yyy_xyy_0[j] - pc_y[j] * tex_yyy_xyy_1[j] + 1.5 * fl1_fx * tex_yy_xyy_0[j] -
                                    1.5 * fl1_fx * tex_yy_xyy_1[j] + fl1_fx * tex_yyy_xy_0[j] - fl1_fx * tex_yyy_xy_1[j];

                tey_yyyy_xyy_0[j] = pa_y[j] * tey_yyy_xyy_0[j] - pc_y[j] * tey_yyy_xyy_1[j] + 1.5 * fl1_fx * tey_yy_xyy_0[j] -
                                    1.5 * fl1_fx * tey_yy_xyy_1[j] + fl1_fx * tey_yyy_xy_0[j] - fl1_fx * tey_yyy_xy_1[j] + ta_yyy_xyy_1[j];

                tez_yyyy_xyy_0[j] = pa_y[j] * tez_yyy_xyy_0[j] - pc_y[j] * tez_yyy_xyy_1[j] + 1.5 * fl1_fx * tez_yy_xyy_0[j] -
                                    1.5 * fl1_fx * tez_yy_xyy_1[j] + fl1_fx * tez_yyy_xy_0[j] - fl1_fx * tez_yyy_xy_1[j];

                tex_yyyy_xyz_0[j] = pa_y[j] * tex_yyy_xyz_0[j] - pc_y[j] * tex_yyy_xyz_1[j] + 1.5 * fl1_fx * tex_yy_xyz_0[j] -
                                    1.5 * fl1_fx * tex_yy_xyz_1[j] + 0.5 * fl1_fx * tex_yyy_xz_0[j] - 0.5 * fl1_fx * tex_yyy_xz_1[j];

                tey_yyyy_xyz_0[j] = pa_y[j] * tey_yyy_xyz_0[j] - pc_y[j] * tey_yyy_xyz_1[j] + 1.5 * fl1_fx * tey_yy_xyz_0[j] -
                                    1.5 * fl1_fx * tey_yy_xyz_1[j] + 0.5 * fl1_fx * tey_yyy_xz_0[j] - 0.5 * fl1_fx * tey_yyy_xz_1[j] +
                                    ta_yyy_xyz_1[j];

                tez_yyyy_xyz_0[j] = pa_y[j] * tez_yyy_xyz_0[j] - pc_y[j] * tez_yyy_xyz_1[j] + 1.5 * fl1_fx * tez_yy_xyz_0[j] -
                                    1.5 * fl1_fx * tez_yy_xyz_1[j] + 0.5 * fl1_fx * tez_yyy_xz_0[j] - 0.5 * fl1_fx * tez_yyy_xz_1[j];

                tex_yyyy_xzz_0[j] =
                    pa_y[j] * tex_yyy_xzz_0[j] - pc_y[j] * tex_yyy_xzz_1[j] + 1.5 * fl1_fx * tex_yy_xzz_0[j] - 1.5 * fl1_fx * tex_yy_xzz_1[j];

                tey_yyyy_xzz_0[j] = pa_y[j] * tey_yyy_xzz_0[j] - pc_y[j] * tey_yyy_xzz_1[j] + 1.5 * fl1_fx * tey_yy_xzz_0[j] -
                                    1.5 * fl1_fx * tey_yy_xzz_1[j] + ta_yyy_xzz_1[j];

                tez_yyyy_xzz_0[j] =
                    pa_y[j] * tez_yyy_xzz_0[j] - pc_y[j] * tez_yyy_xzz_1[j] + 1.5 * fl1_fx * tez_yy_xzz_0[j] - 1.5 * fl1_fx * tez_yy_xzz_1[j];

                tex_yyyy_yyy_0[j] = pa_y[j] * tex_yyy_yyy_0[j] - pc_y[j] * tex_yyy_yyy_1[j] + 1.5 * fl1_fx * tex_yy_yyy_0[j] -
                                    1.5 * fl1_fx * tex_yy_yyy_1[j] + 1.5 * fl1_fx * tex_yyy_yy_0[j] - 1.5 * fl1_fx * tex_yyy_yy_1[j];

                tey_yyyy_yyy_0[j] = pa_y[j] * tey_yyy_yyy_0[j] - pc_y[j] * tey_yyy_yyy_1[j] + 1.5 * fl1_fx * tey_yy_yyy_0[j] -
                                    1.5 * fl1_fx * tey_yy_yyy_1[j] + 1.5 * fl1_fx * tey_yyy_yy_0[j] - 1.5 * fl1_fx * tey_yyy_yy_1[j] +
                                    ta_yyy_yyy_1[j];

                tez_yyyy_yyy_0[j] = pa_y[j] * tez_yyy_yyy_0[j] - pc_y[j] * tez_yyy_yyy_1[j] + 1.5 * fl1_fx * tez_yy_yyy_0[j] -
                                    1.5 * fl1_fx * tez_yy_yyy_1[j] + 1.5 * fl1_fx * tez_yyy_yy_0[j] - 1.5 * fl1_fx * tez_yyy_yy_1[j];

                tex_yyyy_yyz_0[j] = pa_y[j] * tex_yyy_yyz_0[j] - pc_y[j] * tex_yyy_yyz_1[j] + 1.5 * fl1_fx * tex_yy_yyz_0[j] -
                                    1.5 * fl1_fx * tex_yy_yyz_1[j] + fl1_fx * tex_yyy_yz_0[j] - fl1_fx * tex_yyy_yz_1[j];

                tey_yyyy_yyz_0[j] = pa_y[j] * tey_yyy_yyz_0[j] - pc_y[j] * tey_yyy_yyz_1[j] + 1.5 * fl1_fx * tey_yy_yyz_0[j] -
                                    1.5 * fl1_fx * tey_yy_yyz_1[j] + fl1_fx * tey_yyy_yz_0[j] - fl1_fx * tey_yyy_yz_1[j] + ta_yyy_yyz_1[j];

                tez_yyyy_yyz_0[j] = pa_y[j] * tez_yyy_yyz_0[j] - pc_y[j] * tez_yyy_yyz_1[j] + 1.5 * fl1_fx * tez_yy_yyz_0[j] -
                                    1.5 * fl1_fx * tez_yy_yyz_1[j] + fl1_fx * tez_yyy_yz_0[j] - fl1_fx * tez_yyy_yz_1[j];

                tex_yyyy_yzz_0[j] = pa_y[j] * tex_yyy_yzz_0[j] - pc_y[j] * tex_yyy_yzz_1[j] + 1.5 * fl1_fx * tex_yy_yzz_0[j] -
                                    1.5 * fl1_fx * tex_yy_yzz_1[j] + 0.5 * fl1_fx * tex_yyy_zz_0[j] - 0.5 * fl1_fx * tex_yyy_zz_1[j];

                tey_yyyy_yzz_0[j] = pa_y[j] * tey_yyy_yzz_0[j] - pc_y[j] * tey_yyy_yzz_1[j] + 1.5 * fl1_fx * tey_yy_yzz_0[j] -
                                    1.5 * fl1_fx * tey_yy_yzz_1[j] + 0.5 * fl1_fx * tey_yyy_zz_0[j] - 0.5 * fl1_fx * tey_yyy_zz_1[j] +
                                    ta_yyy_yzz_1[j];

                tez_yyyy_yzz_0[j] = pa_y[j] * tez_yyy_yzz_0[j] - pc_y[j] * tez_yyy_yzz_1[j] + 1.5 * fl1_fx * tez_yy_yzz_0[j] -
                                    1.5 * fl1_fx * tez_yy_yzz_1[j] + 0.5 * fl1_fx * tez_yyy_zz_0[j] - 0.5 * fl1_fx * tez_yyy_zz_1[j];

                tex_yyyy_zzz_0[j] =
                    pa_y[j] * tex_yyy_zzz_0[j] - pc_y[j] * tex_yyy_zzz_1[j] + 1.5 * fl1_fx * tex_yy_zzz_0[j] - 1.5 * fl1_fx * tex_yy_zzz_1[j];

                tey_yyyy_zzz_0[j] = pa_y[j] * tey_yyy_zzz_0[j] - pc_y[j] * tey_yyy_zzz_1[j] + 1.5 * fl1_fx * tey_yy_zzz_0[j] -
                                    1.5 * fl1_fx * tey_yy_zzz_1[j] + ta_yyy_zzz_1[j];

                tez_yyyy_zzz_0[j] =
                    pa_y[j] * tez_yyy_zzz_0[j] - pc_y[j] * tez_yyy_zzz_1[j] + 1.5 * fl1_fx * tez_yy_zzz_0[j] - 1.5 * fl1_fx * tez_yy_zzz_1[j];

                tex_yyyz_xxx_0[j] = pa_y[j] * tex_yyz_xxx_0[j] - pc_y[j] * tex_yyz_xxx_1[j] + fl1_fx * tex_yz_xxx_0[j] - fl1_fx * tex_yz_xxx_1[j];

                tey_yyyz_xxx_0[j] =
                    pa_y[j] * tey_yyz_xxx_0[j] - pc_y[j] * tey_yyz_xxx_1[j] + fl1_fx * tey_yz_xxx_0[j] - fl1_fx * tey_yz_xxx_1[j] + ta_yyz_xxx_1[j];

                tez_yyyz_xxx_0[j] = pa_y[j] * tez_yyz_xxx_0[j] - pc_y[j] * tez_yyz_xxx_1[j] + fl1_fx * tez_yz_xxx_0[j] - fl1_fx * tez_yz_xxx_1[j];

                tex_yyyz_xxy_0[j] = pa_y[j] * tex_yyz_xxy_0[j] - pc_y[j] * tex_yyz_xxy_1[j] + fl1_fx * tex_yz_xxy_0[j] - fl1_fx * tex_yz_xxy_1[j] +
                                    0.5 * fl1_fx * tex_yyz_xx_0[j] - 0.5 * fl1_fx * tex_yyz_xx_1[j];

                tey_yyyz_xxy_0[j] = pa_y[j] * tey_yyz_xxy_0[j] - pc_y[j] * tey_yyz_xxy_1[j] + fl1_fx * tey_yz_xxy_0[j] - fl1_fx * tey_yz_xxy_1[j] +
                                    0.5 * fl1_fx * tey_yyz_xx_0[j] - 0.5 * fl1_fx * tey_yyz_xx_1[j] + ta_yyz_xxy_1[j];

                tez_yyyz_xxy_0[j] = pa_y[j] * tez_yyz_xxy_0[j] - pc_y[j] * tez_yyz_xxy_1[j] + fl1_fx * tez_yz_xxy_0[j] - fl1_fx * tez_yz_xxy_1[j] +
                                    0.5 * fl1_fx * tez_yyz_xx_0[j] - 0.5 * fl1_fx * tez_yyz_xx_1[j];

                tex_yyyz_xxz_0[j] = pa_y[j] * tex_yyz_xxz_0[j] - pc_y[j] * tex_yyz_xxz_1[j] + fl1_fx * tex_yz_xxz_0[j] - fl1_fx * tex_yz_xxz_1[j];

                tey_yyyz_xxz_0[j] =
                    pa_y[j] * tey_yyz_xxz_0[j] - pc_y[j] * tey_yyz_xxz_1[j] + fl1_fx * tey_yz_xxz_0[j] - fl1_fx * tey_yz_xxz_1[j] + ta_yyz_xxz_1[j];

                tez_yyyz_xxz_0[j] = pa_y[j] * tez_yyz_xxz_0[j] - pc_y[j] * tez_yyz_xxz_1[j] + fl1_fx * tez_yz_xxz_0[j] - fl1_fx * tez_yz_xxz_1[j];

                tex_yyyz_xyy_0[j] = pa_y[j] * tex_yyz_xyy_0[j] - pc_y[j] * tex_yyz_xyy_1[j] + fl1_fx * tex_yz_xyy_0[j] - fl1_fx * tex_yz_xyy_1[j] +
                                    fl1_fx * tex_yyz_xy_0[j] - fl1_fx * tex_yyz_xy_1[j];

                tey_yyyz_xyy_0[j] = pa_y[j] * tey_yyz_xyy_0[j] - pc_y[j] * tey_yyz_xyy_1[j] + fl1_fx * tey_yz_xyy_0[j] - fl1_fx * tey_yz_xyy_1[j] +
                                    fl1_fx * tey_yyz_xy_0[j] - fl1_fx * tey_yyz_xy_1[j] + ta_yyz_xyy_1[j];

                tez_yyyz_xyy_0[j] = pa_y[j] * tez_yyz_xyy_0[j] - pc_y[j] * tez_yyz_xyy_1[j] + fl1_fx * tez_yz_xyy_0[j] - fl1_fx * tez_yz_xyy_1[j] +
                                    fl1_fx * tez_yyz_xy_0[j] - fl1_fx * tez_yyz_xy_1[j];

                tex_yyyz_xyz_0[j] = pa_y[j] * tex_yyz_xyz_0[j] - pc_y[j] * tex_yyz_xyz_1[j] + fl1_fx * tex_yz_xyz_0[j] - fl1_fx * tex_yz_xyz_1[j] +
                                    0.5 * fl1_fx * tex_yyz_xz_0[j] - 0.5 * fl1_fx * tex_yyz_xz_1[j];

                tey_yyyz_xyz_0[j] = pa_y[j] * tey_yyz_xyz_0[j] - pc_y[j] * tey_yyz_xyz_1[j] + fl1_fx * tey_yz_xyz_0[j] - fl1_fx * tey_yz_xyz_1[j] +
                                    0.5 * fl1_fx * tey_yyz_xz_0[j] - 0.5 * fl1_fx * tey_yyz_xz_1[j] + ta_yyz_xyz_1[j];

                tez_yyyz_xyz_0[j] = pa_y[j] * tez_yyz_xyz_0[j] - pc_y[j] * tez_yyz_xyz_1[j] + fl1_fx * tez_yz_xyz_0[j] - fl1_fx * tez_yz_xyz_1[j] +
                                    0.5 * fl1_fx * tez_yyz_xz_0[j] - 0.5 * fl1_fx * tez_yyz_xz_1[j];

                tex_yyyz_xzz_0[j] = pa_y[j] * tex_yyz_xzz_0[j] - pc_y[j] * tex_yyz_xzz_1[j] + fl1_fx * tex_yz_xzz_0[j] - fl1_fx * tex_yz_xzz_1[j];

                tey_yyyz_xzz_0[j] =
                    pa_y[j] * tey_yyz_xzz_0[j] - pc_y[j] * tey_yyz_xzz_1[j] + fl1_fx * tey_yz_xzz_0[j] - fl1_fx * tey_yz_xzz_1[j] + ta_yyz_xzz_1[j];

                tez_yyyz_xzz_0[j] = pa_y[j] * tez_yyz_xzz_0[j] - pc_y[j] * tez_yyz_xzz_1[j] + fl1_fx * tez_yz_xzz_0[j] - fl1_fx * tez_yz_xzz_1[j];

                tex_yyyz_yyy_0[j] = pa_y[j] * tex_yyz_yyy_0[j] - pc_y[j] * tex_yyz_yyy_1[j] + fl1_fx * tex_yz_yyy_0[j] - fl1_fx * tex_yz_yyy_1[j] +
                                    1.5 * fl1_fx * tex_yyz_yy_0[j] - 1.5 * fl1_fx * tex_yyz_yy_1[j];

                tey_yyyz_yyy_0[j] = pa_y[j] * tey_yyz_yyy_0[j] - pc_y[j] * tey_yyz_yyy_1[j] + fl1_fx * tey_yz_yyy_0[j] - fl1_fx * tey_yz_yyy_1[j] +
                                    1.5 * fl1_fx * tey_yyz_yy_0[j] - 1.5 * fl1_fx * tey_yyz_yy_1[j] + ta_yyz_yyy_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForGF_350_400(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pcDistances,
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

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_4_3_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_2_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

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

            // set up pointers to tensors product of distances R(PC) = P - C

            auto pc_y = pcDistances.data(3 * idx + 1);

            // set up pointers to auxilary integrals

            auto tez_yyz_yyy_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 76);

            auto tex_yyz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * idx + 77);

            auto tey_yyz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 100 * bdim + 100 * idx + 77);

            auto tez_yyz_yyz_0 = primBuffer.data(pidx_e_3_3_m0 + 200 * bdim + 100 * idx + 77);

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

            auto tez_yyz_yyy_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 76);

            auto tex_yyz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * idx + 77);

            auto tey_yyz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 100 * bdim + 100 * idx + 77);

            auto tez_yyz_yyz_1 = primBuffer.data(pidx_e_3_3_m1 + 200 * bdim + 100 * idx + 77);

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

            auto tez_yz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 46);

            auto tex_yz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 47);

            auto tey_yz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 47);

            auto tez_yz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 47);

            auto tex_yz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 48);

            auto tey_yz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 48);

            auto tez_yz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 48);

            auto tex_yz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 49);

            auto tey_yz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 49);

            auto tez_yz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 49);

            auto tex_zz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 50);

            auto tey_zz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 50);

            auto tez_zz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 50);

            auto tex_zz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 51);

            auto tey_zz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 51);

            auto tez_zz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 51);

            auto tex_zz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 52);

            auto tey_zz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 52);

            auto tez_zz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 52);

            auto tex_zz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 53);

            auto tey_zz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 53);

            auto tez_zz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 53);

            auto tex_zz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 54);

            auto tey_zz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 54);

            auto tez_zz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 54);

            auto tex_zz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 55);

            auto tey_zz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 55);

            auto tez_zz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 55);

            auto tex_zz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 56);

            auto tey_zz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 56);

            auto tez_zz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 56);

            auto tex_zz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 57);

            auto tey_zz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 57);

            auto tez_zz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 57);

            auto tex_zz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 58);

            auto tey_zz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 58);

            auto tez_zz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 58);

            auto tex_zz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 59);

            auto tey_zz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 59);

            auto tez_zz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 59);

            auto tez_yz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 46);

            auto tex_yz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 47);

            auto tey_yz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 47);

            auto tez_yz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 47);

            auto tex_yz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 48);

            auto tey_yz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 48);

            auto tez_yz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 48);

            auto tex_yz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 49);

            auto tey_yz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 49);

            auto tez_yz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 49);

            auto tex_zz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 50);

            auto tey_zz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 50);

            auto tez_zz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 50);

            auto tex_zz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 51);

            auto tey_zz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 51);

            auto tez_zz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 51);

            auto tex_zz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 52);

            auto tey_zz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 52);

            auto tez_zz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 52);

            auto tex_zz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 53);

            auto tey_zz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 53);

            auto tez_zz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 53);

            auto tex_zz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 54);

            auto tey_zz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 54);

            auto tez_zz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 54);

            auto tex_zz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 55);

            auto tey_zz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 55);

            auto tez_zz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 55);

            auto tex_zz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 56);

            auto tey_zz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 56);

            auto tez_zz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 56);

            auto tex_zz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 57);

            auto tey_zz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 57);

            auto tez_zz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 57);

            auto tex_zz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 58);

            auto tey_zz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 58);

            auto tez_zz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 58);

            auto tex_zz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 59);

            auto tey_zz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 59);

            auto tez_zz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 59);

            auto tez_yyz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 45);

            auto tex_yyz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 46);

            auto tey_yyz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 46);

            auto tez_yyz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 46);

            auto tex_yyz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 47);

            auto tey_yyz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 47);

            auto tez_yyz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 47);

            auto tex_yzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 48);

            auto tey_yzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 48);

            auto tez_yzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 48);

            auto tex_yzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 49);

            auto tey_yzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 49);

            auto tez_yzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 49);

            auto tex_yzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 50);

            auto tey_yzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 50);

            auto tez_yzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 50);

            auto tex_yzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 51);

            auto tey_yzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 51);

            auto tez_yzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 51);

            auto tex_yzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 52);

            auto tey_yzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 52);

            auto tez_yzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 52);

            auto tex_yzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 53);

            auto tey_yzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 53);

            auto tez_yzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 53);

            auto tex_zzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 54);

            auto tey_zzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 54);

            auto tez_zzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 54);

            auto tex_zzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 55);

            auto tez_yyz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 45);

            auto tex_yyz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 46);

            auto tey_yyz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 46);

            auto tez_yyz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 46);

            auto tex_yyz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 47);

            auto tey_yyz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 47);

            auto tez_yyz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 47);

            auto tex_yzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 48);

            auto tey_yzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 48);

            auto tez_yzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 48);

            auto tex_yzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 49);

            auto tey_yzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 49);

            auto tez_yzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 49);

            auto tex_yzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 50);

            auto tey_yzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 50);

            auto tez_yzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 50);

            auto tex_yzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 51);

            auto tey_yzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 51);

            auto tez_yzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 51);

            auto tex_yzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 52);

            auto tey_yzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 52);

            auto tez_yzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 52);

            auto tex_yzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 53);

            auto tey_yzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 53);

            auto tez_yzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 53);

            auto tex_zzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 54);

            auto tey_zzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 54);

            auto tez_zzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 54);

            auto tex_zzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 55);

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

            auto ta_zzz_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 90);

            auto ta_zzz_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 91);

            auto ta_zzz_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 92);

            // set up pointers to integrals

            auto tez_yyyz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 116);

            auto tex_yyyz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 117);

            auto tey_yyyz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 117);

            auto tez_yyyz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 117);

            auto tex_yyyz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 118);

            auto tey_yyyz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 118);

            auto tez_yyyz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 118);

            auto tex_yyyz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 119);

            auto tey_yyyz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 119);

            auto tez_yyyz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 119);

            auto tex_yyzz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 120);

            auto tey_yyzz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 120);

            auto tez_yyzz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 120);

            auto tex_yyzz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 121);

            auto tey_yyzz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 121);

            auto tez_yyzz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 121);

            auto tex_yyzz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 122);

            auto tey_yyzz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 122);

            auto tez_yyzz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 122);

            auto tex_yyzz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 123);

            auto tey_yyzz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 123);

            auto tez_yyzz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 123);

            auto tex_yyzz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 124);

            auto tey_yyzz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 124);

            auto tez_yyzz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 124);

            auto tex_yyzz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 125);

            auto tey_yyzz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 125);

            auto tez_yyzz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 125);

            auto tex_yyzz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 126);

            auto tey_yyzz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 126);

            auto tez_yyzz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 126);

            auto tex_yyzz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 127);

            auto tey_yyzz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 127);

            auto tez_yyzz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 127);

            auto tex_yyzz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 128);

            auto tey_yyzz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 128);

            auto tez_yyzz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 128);

            auto tex_yyzz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 129);

            auto tey_yyzz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 129);

            auto tez_yyzz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 129);

            auto tex_yzzz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 130);

            auto tey_yzzz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 130);

            auto tez_yzzz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 130);

            auto tex_yzzz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 131);

            auto tey_yzzz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 131);

            auto tez_yzzz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 131);

            auto tex_yzzz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 132);

            auto tey_yzzz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 132);

            auto tez_yzzz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 132);

            auto tex_yzzz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 133);

            // Batch of Integrals (350,400)

            #pragma omp simd aligned(fx, pa_y, pc_y, ta_yyz_yyz_1, ta_yyz_yzz_1, ta_yyz_zzz_1, ta_yzz_xxx_1, \
                                         ta_yzz_xxy_1, ta_yzz_xxz_1, ta_yzz_xyy_1, ta_yzz_xyz_1, ta_yzz_xzz_1, ta_yzz_yyy_1, \
                                         ta_yzz_yyz_1, ta_yzz_yzz_1, ta_yzz_zzz_1, ta_zzz_xxx_1, ta_zzz_xxy_1, ta_zzz_xxz_1, \
                                         tex_yyyz_yyz_0, tex_yyyz_yzz_0, tex_yyyz_zzz_0, tex_yyz_yyz_0, tex_yyz_yyz_1, \
                                         tex_yyz_yz_0, tex_yyz_yz_1, tex_yyz_yzz_0, tex_yyz_yzz_1, tex_yyz_zz_0, \
                                         tex_yyz_zz_1, tex_yyz_zzz_0, tex_yyz_zzz_1, tex_yyzz_xxx_0, tex_yyzz_xxy_0, \
                                         tex_yyzz_xxz_0, tex_yyzz_xyy_0, tex_yyzz_xyz_0, tex_yyzz_xzz_0, tex_yyzz_yyy_0, \
                                         tex_yyzz_yyz_0, tex_yyzz_yzz_0, tex_yyzz_zzz_0, tex_yz_yyz_0, tex_yz_yyz_1, \
                                         tex_yz_yzz_0, tex_yz_yzz_1, tex_yz_zzz_0, tex_yz_zzz_1, tex_yzz_xx_0, tex_yzz_xx_1, \
                                         tex_yzz_xxx_0, tex_yzz_xxx_1, tex_yzz_xxy_0, tex_yzz_xxy_1, tex_yzz_xxz_0, \
                                         tex_yzz_xxz_1, tex_yzz_xy_0, tex_yzz_xy_1, tex_yzz_xyy_0, tex_yzz_xyy_1, \
                                         tex_yzz_xyz_0, tex_yzz_xyz_1, tex_yzz_xz_0, tex_yzz_xz_1, tex_yzz_xzz_0, \
                                         tex_yzz_xzz_1, tex_yzz_yy_0, tex_yzz_yy_1, tex_yzz_yyy_0, tex_yzz_yyy_1, \
                                         tex_yzz_yyz_0, tex_yzz_yyz_1, tex_yzz_yz_0, tex_yzz_yz_1, tex_yzz_yzz_0, \
                                         tex_yzz_yzz_1, tex_yzz_zz_0, tex_yzz_zz_1, tex_yzz_zzz_0, tex_yzz_zzz_1, \
                                         tex_yzzz_xxx_0, tex_yzzz_xxy_0, tex_yzzz_xxz_0, tex_yzzz_xyy_0, tex_zz_xxx_0, \
                                         tex_zz_xxx_1, tex_zz_xxy_0, tex_zz_xxy_1, tex_zz_xxz_0, tex_zz_xxz_1, tex_zz_xyy_0, \
                                         tex_zz_xyy_1, tex_zz_xyz_0, tex_zz_xyz_1, tex_zz_xzz_0, tex_zz_xzz_1, tex_zz_yyy_0, \
                                         tex_zz_yyy_1, tex_zz_yyz_0, tex_zz_yyz_1, tex_zz_yzz_0, tex_zz_yzz_1, tex_zz_zzz_0, \
                                         tex_zz_zzz_1, tex_zzz_xx_0, tex_zzz_xx_1, tex_zzz_xxx_0, tex_zzz_xxx_1, \
                                         tex_zzz_xxy_0, tex_zzz_xxy_1, tex_zzz_xxz_0, tex_zzz_xxz_1, tex_zzz_xy_0, \
                                         tex_zzz_xy_1, tex_zzz_xyy_0, tex_zzz_xyy_1, tey_yyyz_yyz_0, tey_yyyz_yzz_0, \
                                         tey_yyyz_zzz_0, tey_yyz_yyz_0, tey_yyz_yyz_1, tey_yyz_yz_0, tey_yyz_yz_1, \
                                         tey_yyz_yzz_0, tey_yyz_yzz_1, tey_yyz_zz_0, tey_yyz_zz_1, tey_yyz_zzz_0, \
                                         tey_yyz_zzz_1, tey_yyzz_xxx_0, tey_yyzz_xxy_0, tey_yyzz_xxz_0, tey_yyzz_xyy_0, \
                                         tey_yyzz_xyz_0, tey_yyzz_xzz_0, tey_yyzz_yyy_0, tey_yyzz_yyz_0, tey_yyzz_yzz_0, \
                                         tey_yyzz_zzz_0, tey_yz_yyz_0, tey_yz_yyz_1, tey_yz_yzz_0, tey_yz_yzz_1, tey_yz_zzz_0, \
                                         tey_yz_zzz_1, tey_yzz_xx_0, tey_yzz_xx_1, tey_yzz_xxx_0, tey_yzz_xxx_1, \
                                         tey_yzz_xxy_0, tey_yzz_xxy_1, tey_yzz_xxz_0, tey_yzz_xxz_1, tey_yzz_xy_0, \
                                         tey_yzz_xy_1, tey_yzz_xyy_0, tey_yzz_xyy_1, tey_yzz_xyz_0, tey_yzz_xyz_1, \
                                         tey_yzz_xz_0, tey_yzz_xz_1, tey_yzz_xzz_0, tey_yzz_xzz_1, tey_yzz_yy_0, \
                                         tey_yzz_yy_1, tey_yzz_yyy_0, tey_yzz_yyy_1, tey_yzz_yyz_0, tey_yzz_yyz_1, \
                                         tey_yzz_yz_0, tey_yzz_yz_1, tey_yzz_yzz_0, tey_yzz_yzz_1, tey_yzz_zz_0, \
                                         tey_yzz_zz_1, tey_yzz_zzz_0, tey_yzz_zzz_1, tey_yzzz_xxx_0, tey_yzzz_xxy_0, \
                                         tey_yzzz_xxz_0, tey_zz_xxx_0, tey_zz_xxx_1, tey_zz_xxy_0, tey_zz_xxy_1, tey_zz_xxz_0, \
                                         tey_zz_xxz_1, tey_zz_xyy_0, tey_zz_xyy_1, tey_zz_xyz_0, tey_zz_xyz_1, tey_zz_xzz_0, \
                                         tey_zz_xzz_1, tey_zz_yyy_0, tey_zz_yyy_1, tey_zz_yyz_0, tey_zz_yyz_1, tey_zz_yzz_0, \
                                         tey_zz_yzz_1, tey_zz_zzz_0, tey_zz_zzz_1, tey_zzz_xx_0, tey_zzz_xx_1, \
                                         tey_zzz_xxx_0, tey_zzz_xxx_1, tey_zzz_xxy_0, tey_zzz_xxy_1, tey_zzz_xxz_0, \
                                         tey_zzz_xxz_1, tez_yyyz_yyy_0, tez_yyyz_yyz_0, tez_yyyz_yzz_0, tez_yyyz_zzz_0, \
                                         tez_yyz_yy_0, tez_yyz_yy_1, tez_yyz_yyy_0, tez_yyz_yyy_1, tez_yyz_yyz_0, \
                                         tez_yyz_yyz_1, tez_yyz_yz_0, tez_yyz_yz_1, tez_yyz_yzz_0, tez_yyz_yzz_1, \
                                         tez_yyz_zz_0, tez_yyz_zz_1, tez_yyz_zzz_0, tez_yyz_zzz_1, tez_yyzz_xxx_0, \
                                         tez_yyzz_xxy_0, tez_yyzz_xxz_0, tez_yyzz_xyy_0, tez_yyzz_xyz_0, tez_yyzz_xzz_0, \
                                         tez_yyzz_yyy_0, tez_yyzz_yyz_0, tez_yyzz_yzz_0, tez_yyzz_zzz_0, tez_yz_yyy_0, \
                                         tez_yz_yyy_1, tez_yz_yyz_0, tez_yz_yyz_1, tez_yz_yzz_0, tez_yz_yzz_1, tez_yz_zzz_0, \
                                         tez_yz_zzz_1, tez_yzz_xx_0, tez_yzz_xx_1, tez_yzz_xxx_0, tez_yzz_xxx_1, \
                                         tez_yzz_xxy_0, tez_yzz_xxy_1, tez_yzz_xxz_0, tez_yzz_xxz_1, tez_yzz_xy_0, \
                                         tez_yzz_xy_1, tez_yzz_xyy_0, tez_yzz_xyy_1, tez_yzz_xyz_0, tez_yzz_xyz_1, \
                                         tez_yzz_xz_0, tez_yzz_xz_1, tez_yzz_xzz_0, tez_yzz_xzz_1, tez_yzz_yy_0, \
                                         tez_yzz_yy_1, tez_yzz_yyy_0, tez_yzz_yyy_1, tez_yzz_yyz_0, tez_yzz_yyz_1, \
                                         tez_yzz_yz_0, tez_yzz_yz_1, tez_yzz_yzz_0, tez_yzz_yzz_1, tez_yzz_zz_0, \
                                         tez_yzz_zz_1, tez_yzz_zzz_0, tez_yzz_zzz_1, tez_yzzz_xxx_0, tez_yzzz_xxy_0, \
                                         tez_yzzz_xxz_0, tez_zz_xxx_0, tez_zz_xxx_1, tez_zz_xxy_0, tez_zz_xxy_1, tez_zz_xxz_0, \
                                         tez_zz_xxz_1, tez_zz_xyy_0, tez_zz_xyy_1, tez_zz_xyz_0, tez_zz_xyz_1, tez_zz_xzz_0, \
                                         tez_zz_xzz_1, tez_zz_yyy_0, tez_zz_yyy_1, tez_zz_yyz_0, tez_zz_yyz_1, tez_zz_yzz_0, \
                                         tez_zz_yzz_1, tez_zz_zzz_0, tez_zz_zzz_1, tez_zzz_xx_0, tez_zzz_xx_1, \
                                         tez_zzz_xxx_0, tez_zzz_xxx_1, tez_zzz_xxy_0, tez_zzz_xxy_1, tez_zzz_xxz_0, \
                                         tez_zzz_xxz_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tez_yyyz_yyy_0[j] = pa_y[j] * tez_yyz_yyy_0[j] - pc_y[j] * tez_yyz_yyy_1[j] + fl1_fx * tez_yz_yyy_0[j] - fl1_fx * tez_yz_yyy_1[j] +
                                    1.5 * fl1_fx * tez_yyz_yy_0[j] - 1.5 * fl1_fx * tez_yyz_yy_1[j];

                tex_yyyz_yyz_0[j] = pa_y[j] * tex_yyz_yyz_0[j] - pc_y[j] * tex_yyz_yyz_1[j] + fl1_fx * tex_yz_yyz_0[j] - fl1_fx * tex_yz_yyz_1[j] +
                                    fl1_fx * tex_yyz_yz_0[j] - fl1_fx * tex_yyz_yz_1[j];

                tey_yyyz_yyz_0[j] = pa_y[j] * tey_yyz_yyz_0[j] - pc_y[j] * tey_yyz_yyz_1[j] + fl1_fx * tey_yz_yyz_0[j] - fl1_fx * tey_yz_yyz_1[j] +
                                    fl1_fx * tey_yyz_yz_0[j] - fl1_fx * tey_yyz_yz_1[j] + ta_yyz_yyz_1[j];

                tez_yyyz_yyz_0[j] = pa_y[j] * tez_yyz_yyz_0[j] - pc_y[j] * tez_yyz_yyz_1[j] + fl1_fx * tez_yz_yyz_0[j] - fl1_fx * tez_yz_yyz_1[j] +
                                    fl1_fx * tez_yyz_yz_0[j] - fl1_fx * tez_yyz_yz_1[j];

                tex_yyyz_yzz_0[j] = pa_y[j] * tex_yyz_yzz_0[j] - pc_y[j] * tex_yyz_yzz_1[j] + fl1_fx * tex_yz_yzz_0[j] - fl1_fx * tex_yz_yzz_1[j] +
                                    0.5 * fl1_fx * tex_yyz_zz_0[j] - 0.5 * fl1_fx * tex_yyz_zz_1[j];

                tey_yyyz_yzz_0[j] = pa_y[j] * tey_yyz_yzz_0[j] - pc_y[j] * tey_yyz_yzz_1[j] + fl1_fx * tey_yz_yzz_0[j] - fl1_fx * tey_yz_yzz_1[j] +
                                    0.5 * fl1_fx * tey_yyz_zz_0[j] - 0.5 * fl1_fx * tey_yyz_zz_1[j] + ta_yyz_yzz_1[j];

                tez_yyyz_yzz_0[j] = pa_y[j] * tez_yyz_yzz_0[j] - pc_y[j] * tez_yyz_yzz_1[j] + fl1_fx * tez_yz_yzz_0[j] - fl1_fx * tez_yz_yzz_1[j] +
                                    0.5 * fl1_fx * tez_yyz_zz_0[j] - 0.5 * fl1_fx * tez_yyz_zz_1[j];

                tex_yyyz_zzz_0[j] = pa_y[j] * tex_yyz_zzz_0[j] - pc_y[j] * tex_yyz_zzz_1[j] + fl1_fx * tex_yz_zzz_0[j] - fl1_fx * tex_yz_zzz_1[j];

                tey_yyyz_zzz_0[j] =
                    pa_y[j] * tey_yyz_zzz_0[j] - pc_y[j] * tey_yyz_zzz_1[j] + fl1_fx * tey_yz_zzz_0[j] - fl1_fx * tey_yz_zzz_1[j] + ta_yyz_zzz_1[j];

                tez_yyyz_zzz_0[j] = pa_y[j] * tez_yyz_zzz_0[j] - pc_y[j] * tez_yyz_zzz_1[j] + fl1_fx * tez_yz_zzz_0[j] - fl1_fx * tez_yz_zzz_1[j];

                tex_yyzz_xxx_0[j] =
                    pa_y[j] * tex_yzz_xxx_0[j] - pc_y[j] * tex_yzz_xxx_1[j] + 0.5 * fl1_fx * tex_zz_xxx_0[j] - 0.5 * fl1_fx * tex_zz_xxx_1[j];

                tey_yyzz_xxx_0[j] = pa_y[j] * tey_yzz_xxx_0[j] - pc_y[j] * tey_yzz_xxx_1[j] + 0.5 * fl1_fx * tey_zz_xxx_0[j] -
                                    0.5 * fl1_fx * tey_zz_xxx_1[j] + ta_yzz_xxx_1[j];

                tez_yyzz_xxx_0[j] =
                    pa_y[j] * tez_yzz_xxx_0[j] - pc_y[j] * tez_yzz_xxx_1[j] + 0.5 * fl1_fx * tez_zz_xxx_0[j] - 0.5 * fl1_fx * tez_zz_xxx_1[j];

                tex_yyzz_xxy_0[j] = pa_y[j] * tex_yzz_xxy_0[j] - pc_y[j] * tex_yzz_xxy_1[j] + 0.5 * fl1_fx * tex_zz_xxy_0[j] -
                                    0.5 * fl1_fx * tex_zz_xxy_1[j] + 0.5 * fl1_fx * tex_yzz_xx_0[j] - 0.5 * fl1_fx * tex_yzz_xx_1[j];

                tey_yyzz_xxy_0[j] = pa_y[j] * tey_yzz_xxy_0[j] - pc_y[j] * tey_yzz_xxy_1[j] + 0.5 * fl1_fx * tey_zz_xxy_0[j] -
                                    0.5 * fl1_fx * tey_zz_xxy_1[j] + 0.5 * fl1_fx * tey_yzz_xx_0[j] - 0.5 * fl1_fx * tey_yzz_xx_1[j] +
                                    ta_yzz_xxy_1[j];

                tez_yyzz_xxy_0[j] = pa_y[j] * tez_yzz_xxy_0[j] - pc_y[j] * tez_yzz_xxy_1[j] + 0.5 * fl1_fx * tez_zz_xxy_0[j] -
                                    0.5 * fl1_fx * tez_zz_xxy_1[j] + 0.5 * fl1_fx * tez_yzz_xx_0[j] - 0.5 * fl1_fx * tez_yzz_xx_1[j];

                tex_yyzz_xxz_0[j] =
                    pa_y[j] * tex_yzz_xxz_0[j] - pc_y[j] * tex_yzz_xxz_1[j] + 0.5 * fl1_fx * tex_zz_xxz_0[j] - 0.5 * fl1_fx * tex_zz_xxz_1[j];

                tey_yyzz_xxz_0[j] = pa_y[j] * tey_yzz_xxz_0[j] - pc_y[j] * tey_yzz_xxz_1[j] + 0.5 * fl1_fx * tey_zz_xxz_0[j] -
                                    0.5 * fl1_fx * tey_zz_xxz_1[j] + ta_yzz_xxz_1[j];

                tez_yyzz_xxz_0[j] =
                    pa_y[j] * tez_yzz_xxz_0[j] - pc_y[j] * tez_yzz_xxz_1[j] + 0.5 * fl1_fx * tez_zz_xxz_0[j] - 0.5 * fl1_fx * tez_zz_xxz_1[j];

                tex_yyzz_xyy_0[j] = pa_y[j] * tex_yzz_xyy_0[j] - pc_y[j] * tex_yzz_xyy_1[j] + 0.5 * fl1_fx * tex_zz_xyy_0[j] -
                                    0.5 * fl1_fx * tex_zz_xyy_1[j] + fl1_fx * tex_yzz_xy_0[j] - fl1_fx * tex_yzz_xy_1[j];

                tey_yyzz_xyy_0[j] = pa_y[j] * tey_yzz_xyy_0[j] - pc_y[j] * tey_yzz_xyy_1[j] + 0.5 * fl1_fx * tey_zz_xyy_0[j] -
                                    0.5 * fl1_fx * tey_zz_xyy_1[j] + fl1_fx * tey_yzz_xy_0[j] - fl1_fx * tey_yzz_xy_1[j] + ta_yzz_xyy_1[j];

                tez_yyzz_xyy_0[j] = pa_y[j] * tez_yzz_xyy_0[j] - pc_y[j] * tez_yzz_xyy_1[j] + 0.5 * fl1_fx * tez_zz_xyy_0[j] -
                                    0.5 * fl1_fx * tez_zz_xyy_1[j] + fl1_fx * tez_yzz_xy_0[j] - fl1_fx * tez_yzz_xy_1[j];

                tex_yyzz_xyz_0[j] = pa_y[j] * tex_yzz_xyz_0[j] - pc_y[j] * tex_yzz_xyz_1[j] + 0.5 * fl1_fx * tex_zz_xyz_0[j] -
                                    0.5 * fl1_fx * tex_zz_xyz_1[j] + 0.5 * fl1_fx * tex_yzz_xz_0[j] - 0.5 * fl1_fx * tex_yzz_xz_1[j];

                tey_yyzz_xyz_0[j] = pa_y[j] * tey_yzz_xyz_0[j] - pc_y[j] * tey_yzz_xyz_1[j] + 0.5 * fl1_fx * tey_zz_xyz_0[j] -
                                    0.5 * fl1_fx * tey_zz_xyz_1[j] + 0.5 * fl1_fx * tey_yzz_xz_0[j] - 0.5 * fl1_fx * tey_yzz_xz_1[j] +
                                    ta_yzz_xyz_1[j];

                tez_yyzz_xyz_0[j] = pa_y[j] * tez_yzz_xyz_0[j] - pc_y[j] * tez_yzz_xyz_1[j] + 0.5 * fl1_fx * tez_zz_xyz_0[j] -
                                    0.5 * fl1_fx * tez_zz_xyz_1[j] + 0.5 * fl1_fx * tez_yzz_xz_0[j] - 0.5 * fl1_fx * tez_yzz_xz_1[j];

                tex_yyzz_xzz_0[j] =
                    pa_y[j] * tex_yzz_xzz_0[j] - pc_y[j] * tex_yzz_xzz_1[j] + 0.5 * fl1_fx * tex_zz_xzz_0[j] - 0.5 * fl1_fx * tex_zz_xzz_1[j];

                tey_yyzz_xzz_0[j] = pa_y[j] * tey_yzz_xzz_0[j] - pc_y[j] * tey_yzz_xzz_1[j] + 0.5 * fl1_fx * tey_zz_xzz_0[j] -
                                    0.5 * fl1_fx * tey_zz_xzz_1[j] + ta_yzz_xzz_1[j];

                tez_yyzz_xzz_0[j] =
                    pa_y[j] * tez_yzz_xzz_0[j] - pc_y[j] * tez_yzz_xzz_1[j] + 0.5 * fl1_fx * tez_zz_xzz_0[j] - 0.5 * fl1_fx * tez_zz_xzz_1[j];

                tex_yyzz_yyy_0[j] = pa_y[j] * tex_yzz_yyy_0[j] - pc_y[j] * tex_yzz_yyy_1[j] + 0.5 * fl1_fx * tex_zz_yyy_0[j] -
                                    0.5 * fl1_fx * tex_zz_yyy_1[j] + 1.5 * fl1_fx * tex_yzz_yy_0[j] - 1.5 * fl1_fx * tex_yzz_yy_1[j];

                tey_yyzz_yyy_0[j] = pa_y[j] * tey_yzz_yyy_0[j] - pc_y[j] * tey_yzz_yyy_1[j] + 0.5 * fl1_fx * tey_zz_yyy_0[j] -
                                    0.5 * fl1_fx * tey_zz_yyy_1[j] + 1.5 * fl1_fx * tey_yzz_yy_0[j] - 1.5 * fl1_fx * tey_yzz_yy_1[j] +
                                    ta_yzz_yyy_1[j];

                tez_yyzz_yyy_0[j] = pa_y[j] * tez_yzz_yyy_0[j] - pc_y[j] * tez_yzz_yyy_1[j] + 0.5 * fl1_fx * tez_zz_yyy_0[j] -
                                    0.5 * fl1_fx * tez_zz_yyy_1[j] + 1.5 * fl1_fx * tez_yzz_yy_0[j] - 1.5 * fl1_fx * tez_yzz_yy_1[j];

                tex_yyzz_yyz_0[j] = pa_y[j] * tex_yzz_yyz_0[j] - pc_y[j] * tex_yzz_yyz_1[j] + 0.5 * fl1_fx * tex_zz_yyz_0[j] -
                                    0.5 * fl1_fx * tex_zz_yyz_1[j] + fl1_fx * tex_yzz_yz_0[j] - fl1_fx * tex_yzz_yz_1[j];

                tey_yyzz_yyz_0[j] = pa_y[j] * tey_yzz_yyz_0[j] - pc_y[j] * tey_yzz_yyz_1[j] + 0.5 * fl1_fx * tey_zz_yyz_0[j] -
                                    0.5 * fl1_fx * tey_zz_yyz_1[j] + fl1_fx * tey_yzz_yz_0[j] - fl1_fx * tey_yzz_yz_1[j] + ta_yzz_yyz_1[j];

                tez_yyzz_yyz_0[j] = pa_y[j] * tez_yzz_yyz_0[j] - pc_y[j] * tez_yzz_yyz_1[j] + 0.5 * fl1_fx * tez_zz_yyz_0[j] -
                                    0.5 * fl1_fx * tez_zz_yyz_1[j] + fl1_fx * tez_yzz_yz_0[j] - fl1_fx * tez_yzz_yz_1[j];

                tex_yyzz_yzz_0[j] = pa_y[j] * tex_yzz_yzz_0[j] - pc_y[j] * tex_yzz_yzz_1[j] + 0.5 * fl1_fx * tex_zz_yzz_0[j] -
                                    0.5 * fl1_fx * tex_zz_yzz_1[j] + 0.5 * fl1_fx * tex_yzz_zz_0[j] - 0.5 * fl1_fx * tex_yzz_zz_1[j];

                tey_yyzz_yzz_0[j] = pa_y[j] * tey_yzz_yzz_0[j] - pc_y[j] * tey_yzz_yzz_1[j] + 0.5 * fl1_fx * tey_zz_yzz_0[j] -
                                    0.5 * fl1_fx * tey_zz_yzz_1[j] + 0.5 * fl1_fx * tey_yzz_zz_0[j] - 0.5 * fl1_fx * tey_yzz_zz_1[j] +
                                    ta_yzz_yzz_1[j];

                tez_yyzz_yzz_0[j] = pa_y[j] * tez_yzz_yzz_0[j] - pc_y[j] * tez_yzz_yzz_1[j] + 0.5 * fl1_fx * tez_zz_yzz_0[j] -
                                    0.5 * fl1_fx * tez_zz_yzz_1[j] + 0.5 * fl1_fx * tez_yzz_zz_0[j] - 0.5 * fl1_fx * tez_yzz_zz_1[j];

                tex_yyzz_zzz_0[j] =
                    pa_y[j] * tex_yzz_zzz_0[j] - pc_y[j] * tex_yzz_zzz_1[j] + 0.5 * fl1_fx * tex_zz_zzz_0[j] - 0.5 * fl1_fx * tex_zz_zzz_1[j];

                tey_yyzz_zzz_0[j] = pa_y[j] * tey_yzz_zzz_0[j] - pc_y[j] * tey_yzz_zzz_1[j] + 0.5 * fl1_fx * tey_zz_zzz_0[j] -
                                    0.5 * fl1_fx * tey_zz_zzz_1[j] + ta_yzz_zzz_1[j];

                tez_yyzz_zzz_0[j] =
                    pa_y[j] * tez_yzz_zzz_0[j] - pc_y[j] * tez_yzz_zzz_1[j] + 0.5 * fl1_fx * tez_zz_zzz_0[j] - 0.5 * fl1_fx * tez_zz_zzz_1[j];

                tex_yzzz_xxx_0[j] = pa_y[j] * tex_zzz_xxx_0[j] - pc_y[j] * tex_zzz_xxx_1[j];

                tey_yzzz_xxx_0[j] = pa_y[j] * tey_zzz_xxx_0[j] - pc_y[j] * tey_zzz_xxx_1[j] + ta_zzz_xxx_1[j];

                tez_yzzz_xxx_0[j] = pa_y[j] * tez_zzz_xxx_0[j] - pc_y[j] * tez_zzz_xxx_1[j];

                tex_yzzz_xxy_0[j] =
                    pa_y[j] * tex_zzz_xxy_0[j] - pc_y[j] * tex_zzz_xxy_1[j] + 0.5 * fl1_fx * tex_zzz_xx_0[j] - 0.5 * fl1_fx * tex_zzz_xx_1[j];

                tey_yzzz_xxy_0[j] = pa_y[j] * tey_zzz_xxy_0[j] - pc_y[j] * tey_zzz_xxy_1[j] + 0.5 * fl1_fx * tey_zzz_xx_0[j] -
                                    0.5 * fl1_fx * tey_zzz_xx_1[j] + ta_zzz_xxy_1[j];

                tez_yzzz_xxy_0[j] =
                    pa_y[j] * tez_zzz_xxy_0[j] - pc_y[j] * tez_zzz_xxy_1[j] + 0.5 * fl1_fx * tez_zzz_xx_0[j] - 0.5 * fl1_fx * tez_zzz_xx_1[j];

                tex_yzzz_xxz_0[j] = pa_y[j] * tex_zzz_xxz_0[j] - pc_y[j] * tex_zzz_xxz_1[j];

                tey_yzzz_xxz_0[j] = pa_y[j] * tey_zzz_xxz_0[j] - pc_y[j] * tey_zzz_xxz_1[j] + ta_zzz_xxz_1[j];

                tez_yzzz_xxz_0[j] = pa_y[j] * tez_zzz_xxz_0[j] - pc_y[j] * tez_zzz_xxz_1[j];

                tex_yzzz_xyy_0[j] = pa_y[j] * tex_zzz_xyy_0[j] - pc_y[j] * tex_zzz_xyy_1[j] + fl1_fx * tex_zzz_xy_0[j] - fl1_fx * tex_zzz_xy_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForGF_400_450(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pcDistances,
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

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_4_3_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_3_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_2_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tex_zz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 50);

            auto tey_zz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 50);

            auto tez_zz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 50);

            auto tex_zz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 51);

            auto tey_zz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 51);

            auto tez_zz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 51);

            auto tex_zz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 52);

            auto tey_zz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 52);

            auto tez_zz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 52);

            auto tex_zz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 53);

            auto tey_zz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 53);

            auto tez_zz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 53);

            auto tex_zz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 54);

            auto tey_zz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 54);

            auto tez_zz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 54);

            auto tex_zz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 55);

            auto tey_zz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 55);

            auto tez_zz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 55);

            auto tex_zz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 56);

            auto tey_zz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 56);

            auto tez_zz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 56);

            auto tex_zz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 57);

            auto tey_zz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 57);

            auto tez_zz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 57);

            auto tex_zz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 58);

            auto tey_zz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 58);

            auto tez_zz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 58);

            auto tex_zz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 59);

            auto tey_zz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 59);

            auto tez_zz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 59);

            auto tex_zz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 50);

            auto tey_zz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 50);

            auto tez_zz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 50);

            auto tex_zz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 51);

            auto tey_zz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 51);

            auto tez_zz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 51);

            auto tex_zz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 52);

            auto tey_zz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 52);

            auto tez_zz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 52);

            auto tex_zz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 53);

            auto tey_zz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 53);

            auto tez_zz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 53);

            auto tex_zz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 54);

            auto tey_zz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 54);

            auto tez_zz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 54);

            auto tex_zz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 55);

            auto tey_zz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 55);

            auto tez_zz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 55);

            auto tex_zz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 56);

            auto tey_zz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 56);

            auto tez_zz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 56);

            auto tex_zz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 57);

            auto tey_zz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 57);

            auto tez_zz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 57);

            auto tex_zz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 58);

            auto tey_zz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 58);

            auto tez_zz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 58);

            auto tex_zz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 59);

            auto tey_zz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 59);

            auto tez_zz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 59);

            auto tex_zzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 54);

            auto tey_zzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 54);

            auto tez_zzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 54);

            auto tex_zzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 55);

            auto tey_zzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 55);

            auto tez_zzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 55);

            auto tex_zzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 56);

            auto tey_zzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 56);

            auto tez_zzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 56);

            auto tex_zzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 57);

            auto tey_zzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 57);

            auto tez_zzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 57);

            auto tex_zzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 58);

            auto tey_zzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 58);

            auto tez_zzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 58);

            auto tex_zzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 59);

            auto tey_zzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 59);

            auto tez_zzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 59);

            auto tex_zzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 54);

            auto tey_zzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 54);

            auto tez_zzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 54);

            auto tex_zzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 55);

            auto tey_zzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 55);

            auto tez_zzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 55);

            auto tex_zzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 56);

            auto tey_zzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 56);

            auto tez_zzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 56);

            auto tex_zzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 57);

            auto tey_zzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 57);

            auto tez_zzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 57);

            auto tex_zzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 58);

            auto tey_zzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 58);

            auto tez_zzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 58);

            auto tex_zzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 59);

            auto tey_zzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 59);

            auto tez_zzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 59);

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

            auto tey_yzzz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 133);

            auto tez_yzzz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 133);

            auto tex_yzzz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 134);

            auto tey_yzzz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 134);

            auto tez_yzzz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 134);

            auto tex_yzzz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 135);

            auto tey_yzzz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 135);

            auto tez_yzzz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 135);

            auto tex_yzzz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 136);

            auto tey_yzzz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 136);

            auto tez_yzzz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 136);

            auto tex_yzzz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 137);

            auto tey_yzzz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 137);

            auto tez_yzzz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 137);

            auto tex_yzzz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 138);

            auto tey_yzzz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 138);

            auto tez_yzzz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 138);

            auto tex_yzzz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 139);

            auto tey_yzzz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 139);

            auto tez_yzzz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 139);

            auto tex_zzzz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 140);

            auto tey_zzzz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 140);

            auto tez_zzzz_xxx_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 140);

            auto tex_zzzz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 141);

            auto tey_zzzz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 141);

            auto tez_zzzz_xxy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 141);

            auto tex_zzzz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 142);

            auto tey_zzzz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 142);

            auto tez_zzzz_xxz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 142);

            auto tex_zzzz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 143);

            auto tey_zzzz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 143);

            auto tez_zzzz_xyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 143);

            auto tex_zzzz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 144);

            auto tey_zzzz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 144);

            auto tez_zzzz_xyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 144);

            auto tex_zzzz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 145);

            auto tey_zzzz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 145);

            auto tez_zzzz_xzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 145);

            auto tex_zzzz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 146);

            auto tey_zzzz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 146);

            auto tez_zzzz_yyy_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 146);

            auto tex_zzzz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 147);

            auto tey_zzzz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 147);

            auto tez_zzzz_yyz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 147);

            auto tex_zzzz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 148);

            auto tey_zzzz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 148);

            auto tez_zzzz_yzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 148);

            auto tex_zzzz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * idx + 149);

            auto tey_zzzz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 150 * bdim + 150 * idx + 149);

            auto tez_zzzz_zzz_0 = primBuffer.data(pidx_e_4_3_m0 + 300 * bdim + 150 * idx + 149);

            // Batch of Integrals (400,450)

            #pragma omp simd aligned(fx, pa_y, pa_z, pc_y, pc_z, ta_zzz_xxx_1, ta_zzz_xxy_1, ta_zzz_xxz_1, \
                                         ta_zzz_xyy_1, ta_zzz_xyz_1, ta_zzz_xzz_1, ta_zzz_yyy_1, ta_zzz_yyz_1, ta_zzz_yzz_1, \
                                         ta_zzz_zzz_1, tex_yzzz_xyz_0, tex_yzzz_xzz_0, tex_yzzz_yyy_0, tex_yzzz_yyz_0, \
                                         tex_yzzz_yzz_0, tex_yzzz_zzz_0, tex_zz_xxx_0, tex_zz_xxx_1, tex_zz_xxy_0, \
                                         tex_zz_xxy_1, tex_zz_xxz_0, tex_zz_xxz_1, tex_zz_xyy_0, tex_zz_xyy_1, tex_zz_xyz_0, \
                                         tex_zz_xyz_1, tex_zz_xzz_0, tex_zz_xzz_1, tex_zz_yyy_0, tex_zz_yyy_1, tex_zz_yyz_0, \
                                         tex_zz_yyz_1, tex_zz_yzz_0, tex_zz_yzz_1, tex_zz_zzz_0, tex_zz_zzz_1, tex_zzz_xx_0, \
                                         tex_zzz_xx_1, tex_zzz_xxx_0, tex_zzz_xxx_1, tex_zzz_xxy_0, tex_zzz_xxy_1, \
                                         tex_zzz_xxz_0, tex_zzz_xxz_1, tex_zzz_xy_0, tex_zzz_xy_1, tex_zzz_xyy_0, \
                                         tex_zzz_xyy_1, tex_zzz_xyz_0, tex_zzz_xyz_1, tex_zzz_xz_0, tex_zzz_xz_1, \
                                         tex_zzz_xzz_0, tex_zzz_xzz_1, tex_zzz_yy_0, tex_zzz_yy_1, tex_zzz_yyy_0, \
                                         tex_zzz_yyy_1, tex_zzz_yyz_0, tex_zzz_yyz_1, tex_zzz_yz_0, tex_zzz_yz_1, \
                                         tex_zzz_yzz_0, tex_zzz_yzz_1, tex_zzz_zz_0, tex_zzz_zz_1, tex_zzz_zzz_0, \
                                         tex_zzz_zzz_1, tex_zzzz_xxx_0, tex_zzzz_xxy_0, tex_zzzz_xxz_0, tex_zzzz_xyy_0, \
                                         tex_zzzz_xyz_0, tex_zzzz_xzz_0, tex_zzzz_yyy_0, tex_zzzz_yyz_0, tex_zzzz_yzz_0, \
                                         tex_zzzz_zzz_0, tey_yzzz_xyy_0, tey_yzzz_xyz_0, tey_yzzz_xzz_0, tey_yzzz_yyy_0, \
                                         tey_yzzz_yyz_0, tey_yzzz_yzz_0, tey_yzzz_zzz_0, tey_zz_xxx_0, tey_zz_xxx_1, \
                                         tey_zz_xxy_0, tey_zz_xxy_1, tey_zz_xxz_0, tey_zz_xxz_1, tey_zz_xyy_0, tey_zz_xyy_1, \
                                         tey_zz_xyz_0, tey_zz_xyz_1, tey_zz_xzz_0, tey_zz_xzz_1, tey_zz_yyy_0, tey_zz_yyy_1, \
                                         tey_zz_yyz_0, tey_zz_yyz_1, tey_zz_yzz_0, tey_zz_yzz_1, tey_zz_zzz_0, tey_zz_zzz_1, \
                                         tey_zzz_xx_0, tey_zzz_xx_1, tey_zzz_xxx_0, tey_zzz_xxx_1, tey_zzz_xxy_0, \
                                         tey_zzz_xxy_1, tey_zzz_xxz_0, tey_zzz_xxz_1, tey_zzz_xy_0, tey_zzz_xy_1, \
                                         tey_zzz_xyy_0, tey_zzz_xyy_1, tey_zzz_xyz_0, tey_zzz_xyz_1, tey_zzz_xz_0, \
                                         tey_zzz_xz_1, tey_zzz_xzz_0, tey_zzz_xzz_1, tey_zzz_yy_0, tey_zzz_yy_1, \
                                         tey_zzz_yyy_0, tey_zzz_yyy_1, tey_zzz_yyz_0, tey_zzz_yyz_1, tey_zzz_yz_0, \
                                         tey_zzz_yz_1, tey_zzz_yzz_0, tey_zzz_yzz_1, tey_zzz_zz_0, tey_zzz_zz_1, \
                                         tey_zzz_zzz_0, tey_zzz_zzz_1, tey_zzzz_xxx_0, tey_zzzz_xxy_0, tey_zzzz_xxz_0, \
                                         tey_zzzz_xyy_0, tey_zzzz_xyz_0, tey_zzzz_xzz_0, tey_zzzz_yyy_0, tey_zzzz_yyz_0, \
                                         tey_zzzz_yzz_0, tey_zzzz_zzz_0, tez_yzzz_xyy_0, tez_yzzz_xyz_0, tez_yzzz_xzz_0, \
                                         tez_yzzz_yyy_0, tez_yzzz_yyz_0, tez_yzzz_yzz_0, tez_yzzz_zzz_0, tez_zz_xxx_0, \
                                         tez_zz_xxx_1, tez_zz_xxy_0, tez_zz_xxy_1, tez_zz_xxz_0, tez_zz_xxz_1, tez_zz_xyy_0, \
                                         tez_zz_xyy_1, tez_zz_xyz_0, tez_zz_xyz_1, tez_zz_xzz_0, tez_zz_xzz_1, tez_zz_yyy_0, \
                                         tez_zz_yyy_1, tez_zz_yyz_0, tez_zz_yyz_1, tez_zz_yzz_0, tez_zz_yzz_1, tez_zz_zzz_0, \
                                         tez_zz_zzz_1, tez_zzz_xx_0, tez_zzz_xx_1, tez_zzz_xxx_0, tez_zzz_xxx_1, \
                                         tez_zzz_xxy_0, tez_zzz_xxy_1, tez_zzz_xxz_0, tez_zzz_xxz_1, tez_zzz_xy_0, \
                                         tez_zzz_xy_1, tez_zzz_xyy_0, tez_zzz_xyy_1, tez_zzz_xyz_0, tez_zzz_xyz_1, \
                                         tez_zzz_xz_0, tez_zzz_xz_1, tez_zzz_xzz_0, tez_zzz_xzz_1, tez_zzz_yy_0, \
                                         tez_zzz_yy_1, tez_zzz_yyy_0, tez_zzz_yyy_1, tez_zzz_yyz_0, tez_zzz_yyz_1, \
                                         tez_zzz_yz_0, tez_zzz_yz_1, tez_zzz_yzz_0, tez_zzz_yzz_1, tez_zzz_zz_0, \
                                         tez_zzz_zz_1, tez_zzz_zzz_0, tez_zzz_zzz_1, tez_zzzz_xxx_0, tez_zzzz_xxy_0, \
                                         tez_zzzz_xxz_0, tez_zzzz_xyy_0, tez_zzzz_xyz_0, tez_zzzz_xzz_0, tez_zzzz_yyy_0, \
                                         tez_zzzz_yyz_0, tez_zzzz_yzz_0, tez_zzzz_zzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tey_yzzz_xyy_0[j] =
                    pa_y[j] * tey_zzz_xyy_0[j] - pc_y[j] * tey_zzz_xyy_1[j] + fl1_fx * tey_zzz_xy_0[j] - fl1_fx * tey_zzz_xy_1[j] + ta_zzz_xyy_1[j];

                tez_yzzz_xyy_0[j] = pa_y[j] * tez_zzz_xyy_0[j] - pc_y[j] * tez_zzz_xyy_1[j] + fl1_fx * tez_zzz_xy_0[j] - fl1_fx * tez_zzz_xy_1[j];

                tex_yzzz_xyz_0[j] =
                    pa_y[j] * tex_zzz_xyz_0[j] - pc_y[j] * tex_zzz_xyz_1[j] + 0.5 * fl1_fx * tex_zzz_xz_0[j] - 0.5 * fl1_fx * tex_zzz_xz_1[j];

                tey_yzzz_xyz_0[j] = pa_y[j] * tey_zzz_xyz_0[j] - pc_y[j] * tey_zzz_xyz_1[j] + 0.5 * fl1_fx * tey_zzz_xz_0[j] -
                                    0.5 * fl1_fx * tey_zzz_xz_1[j] + ta_zzz_xyz_1[j];

                tez_yzzz_xyz_0[j] =
                    pa_y[j] * tez_zzz_xyz_0[j] - pc_y[j] * tez_zzz_xyz_1[j] + 0.5 * fl1_fx * tez_zzz_xz_0[j] - 0.5 * fl1_fx * tez_zzz_xz_1[j];

                tex_yzzz_xzz_0[j] = pa_y[j] * tex_zzz_xzz_0[j] - pc_y[j] * tex_zzz_xzz_1[j];

                tey_yzzz_xzz_0[j] = pa_y[j] * tey_zzz_xzz_0[j] - pc_y[j] * tey_zzz_xzz_1[j] + ta_zzz_xzz_1[j];

                tez_yzzz_xzz_0[j] = pa_y[j] * tez_zzz_xzz_0[j] - pc_y[j] * tez_zzz_xzz_1[j];

                tex_yzzz_yyy_0[j] =
                    pa_y[j] * tex_zzz_yyy_0[j] - pc_y[j] * tex_zzz_yyy_1[j] + 1.5 * fl1_fx * tex_zzz_yy_0[j] - 1.5 * fl1_fx * tex_zzz_yy_1[j];

                tey_yzzz_yyy_0[j] = pa_y[j] * tey_zzz_yyy_0[j] - pc_y[j] * tey_zzz_yyy_1[j] + 1.5 * fl1_fx * tey_zzz_yy_0[j] -
                                    1.5 * fl1_fx * tey_zzz_yy_1[j] + ta_zzz_yyy_1[j];

                tez_yzzz_yyy_0[j] =
                    pa_y[j] * tez_zzz_yyy_0[j] - pc_y[j] * tez_zzz_yyy_1[j] + 1.5 * fl1_fx * tez_zzz_yy_0[j] - 1.5 * fl1_fx * tez_zzz_yy_1[j];

                tex_yzzz_yyz_0[j] = pa_y[j] * tex_zzz_yyz_0[j] - pc_y[j] * tex_zzz_yyz_1[j] + fl1_fx * tex_zzz_yz_0[j] - fl1_fx * tex_zzz_yz_1[j];

                tey_yzzz_yyz_0[j] =
                    pa_y[j] * tey_zzz_yyz_0[j] - pc_y[j] * tey_zzz_yyz_1[j] + fl1_fx * tey_zzz_yz_0[j] - fl1_fx * tey_zzz_yz_1[j] + ta_zzz_yyz_1[j];

                tez_yzzz_yyz_0[j] = pa_y[j] * tez_zzz_yyz_0[j] - pc_y[j] * tez_zzz_yyz_1[j] + fl1_fx * tez_zzz_yz_0[j] - fl1_fx * tez_zzz_yz_1[j];

                tex_yzzz_yzz_0[j] =
                    pa_y[j] * tex_zzz_yzz_0[j] - pc_y[j] * tex_zzz_yzz_1[j] + 0.5 * fl1_fx * tex_zzz_zz_0[j] - 0.5 * fl1_fx * tex_zzz_zz_1[j];

                tey_yzzz_yzz_0[j] = pa_y[j] * tey_zzz_yzz_0[j] - pc_y[j] * tey_zzz_yzz_1[j] + 0.5 * fl1_fx * tey_zzz_zz_0[j] -
                                    0.5 * fl1_fx * tey_zzz_zz_1[j] + ta_zzz_yzz_1[j];

                tez_yzzz_yzz_0[j] =
                    pa_y[j] * tez_zzz_yzz_0[j] - pc_y[j] * tez_zzz_yzz_1[j] + 0.5 * fl1_fx * tez_zzz_zz_0[j] - 0.5 * fl1_fx * tez_zzz_zz_1[j];

                tex_yzzz_zzz_0[j] = pa_y[j] * tex_zzz_zzz_0[j] - pc_y[j] * tex_zzz_zzz_1[j];

                tey_yzzz_zzz_0[j] = pa_y[j] * tey_zzz_zzz_0[j] - pc_y[j] * tey_zzz_zzz_1[j] + ta_zzz_zzz_1[j];

                tez_yzzz_zzz_0[j] = pa_y[j] * tez_zzz_zzz_0[j] - pc_y[j] * tez_zzz_zzz_1[j];

                tex_zzzz_xxx_0[j] =
                    pa_z[j] * tex_zzz_xxx_0[j] - pc_z[j] * tex_zzz_xxx_1[j] + 1.5 * fl1_fx * tex_zz_xxx_0[j] - 1.5 * fl1_fx * tex_zz_xxx_1[j];

                tey_zzzz_xxx_0[j] =
                    pa_z[j] * tey_zzz_xxx_0[j] - pc_z[j] * tey_zzz_xxx_1[j] + 1.5 * fl1_fx * tey_zz_xxx_0[j] - 1.5 * fl1_fx * tey_zz_xxx_1[j];

                tez_zzzz_xxx_0[j] = pa_z[j] * tez_zzz_xxx_0[j] - pc_z[j] * tez_zzz_xxx_1[j] + 1.5 * fl1_fx * tez_zz_xxx_0[j] -
                                    1.5 * fl1_fx * tez_zz_xxx_1[j] + ta_zzz_xxx_1[j];

                tex_zzzz_xxy_0[j] =
                    pa_z[j] * tex_zzz_xxy_0[j] - pc_z[j] * tex_zzz_xxy_1[j] + 1.5 * fl1_fx * tex_zz_xxy_0[j] - 1.5 * fl1_fx * tex_zz_xxy_1[j];

                tey_zzzz_xxy_0[j] =
                    pa_z[j] * tey_zzz_xxy_0[j] - pc_z[j] * tey_zzz_xxy_1[j] + 1.5 * fl1_fx * tey_zz_xxy_0[j] - 1.5 * fl1_fx * tey_zz_xxy_1[j];

                tez_zzzz_xxy_0[j] = pa_z[j] * tez_zzz_xxy_0[j] - pc_z[j] * tez_zzz_xxy_1[j] + 1.5 * fl1_fx * tez_zz_xxy_0[j] -
                                    1.5 * fl1_fx * tez_zz_xxy_1[j] + ta_zzz_xxy_1[j];

                tex_zzzz_xxz_0[j] = pa_z[j] * tex_zzz_xxz_0[j] - pc_z[j] * tex_zzz_xxz_1[j] + 1.5 * fl1_fx * tex_zz_xxz_0[j] -
                                    1.5 * fl1_fx * tex_zz_xxz_1[j] + 0.5 * fl1_fx * tex_zzz_xx_0[j] - 0.5 * fl1_fx * tex_zzz_xx_1[j];

                tey_zzzz_xxz_0[j] = pa_z[j] * tey_zzz_xxz_0[j] - pc_z[j] * tey_zzz_xxz_1[j] + 1.5 * fl1_fx * tey_zz_xxz_0[j] -
                                    1.5 * fl1_fx * tey_zz_xxz_1[j] + 0.5 * fl1_fx * tey_zzz_xx_0[j] - 0.5 * fl1_fx * tey_zzz_xx_1[j];

                tez_zzzz_xxz_0[j] = pa_z[j] * tez_zzz_xxz_0[j] - pc_z[j] * tez_zzz_xxz_1[j] + 1.5 * fl1_fx * tez_zz_xxz_0[j] -
                                    1.5 * fl1_fx * tez_zz_xxz_1[j] + 0.5 * fl1_fx * tez_zzz_xx_0[j] - 0.5 * fl1_fx * tez_zzz_xx_1[j] +
                                    ta_zzz_xxz_1[j];

                tex_zzzz_xyy_0[j] =
                    pa_z[j] * tex_zzz_xyy_0[j] - pc_z[j] * tex_zzz_xyy_1[j] + 1.5 * fl1_fx * tex_zz_xyy_0[j] - 1.5 * fl1_fx * tex_zz_xyy_1[j];

                tey_zzzz_xyy_0[j] =
                    pa_z[j] * tey_zzz_xyy_0[j] - pc_z[j] * tey_zzz_xyy_1[j] + 1.5 * fl1_fx * tey_zz_xyy_0[j] - 1.5 * fl1_fx * tey_zz_xyy_1[j];

                tez_zzzz_xyy_0[j] = pa_z[j] * tez_zzz_xyy_0[j] - pc_z[j] * tez_zzz_xyy_1[j] + 1.5 * fl1_fx * tez_zz_xyy_0[j] -
                                    1.5 * fl1_fx * tez_zz_xyy_1[j] + ta_zzz_xyy_1[j];

                tex_zzzz_xyz_0[j] = pa_z[j] * tex_zzz_xyz_0[j] - pc_z[j] * tex_zzz_xyz_1[j] + 1.5 * fl1_fx * tex_zz_xyz_0[j] -
                                    1.5 * fl1_fx * tex_zz_xyz_1[j] + 0.5 * fl1_fx * tex_zzz_xy_0[j] - 0.5 * fl1_fx * tex_zzz_xy_1[j];

                tey_zzzz_xyz_0[j] = pa_z[j] * tey_zzz_xyz_0[j] - pc_z[j] * tey_zzz_xyz_1[j] + 1.5 * fl1_fx * tey_zz_xyz_0[j] -
                                    1.5 * fl1_fx * tey_zz_xyz_1[j] + 0.5 * fl1_fx * tey_zzz_xy_0[j] - 0.5 * fl1_fx * tey_zzz_xy_1[j];

                tez_zzzz_xyz_0[j] = pa_z[j] * tez_zzz_xyz_0[j] - pc_z[j] * tez_zzz_xyz_1[j] + 1.5 * fl1_fx * tez_zz_xyz_0[j] -
                                    1.5 * fl1_fx * tez_zz_xyz_1[j] + 0.5 * fl1_fx * tez_zzz_xy_0[j] - 0.5 * fl1_fx * tez_zzz_xy_1[j] +
                                    ta_zzz_xyz_1[j];

                tex_zzzz_xzz_0[j] = pa_z[j] * tex_zzz_xzz_0[j] - pc_z[j] * tex_zzz_xzz_1[j] + 1.5 * fl1_fx * tex_zz_xzz_0[j] -
                                    1.5 * fl1_fx * tex_zz_xzz_1[j] + fl1_fx * tex_zzz_xz_0[j] - fl1_fx * tex_zzz_xz_1[j];

                tey_zzzz_xzz_0[j] = pa_z[j] * tey_zzz_xzz_0[j] - pc_z[j] * tey_zzz_xzz_1[j] + 1.5 * fl1_fx * tey_zz_xzz_0[j] -
                                    1.5 * fl1_fx * tey_zz_xzz_1[j] + fl1_fx * tey_zzz_xz_0[j] - fl1_fx * tey_zzz_xz_1[j];

                tez_zzzz_xzz_0[j] = pa_z[j] * tez_zzz_xzz_0[j] - pc_z[j] * tez_zzz_xzz_1[j] + 1.5 * fl1_fx * tez_zz_xzz_0[j] -
                                    1.5 * fl1_fx * tez_zz_xzz_1[j] + fl1_fx * tez_zzz_xz_0[j] - fl1_fx * tez_zzz_xz_1[j] + ta_zzz_xzz_1[j];

                tex_zzzz_yyy_0[j] =
                    pa_z[j] * tex_zzz_yyy_0[j] - pc_z[j] * tex_zzz_yyy_1[j] + 1.5 * fl1_fx * tex_zz_yyy_0[j] - 1.5 * fl1_fx * tex_zz_yyy_1[j];

                tey_zzzz_yyy_0[j] =
                    pa_z[j] * tey_zzz_yyy_0[j] - pc_z[j] * tey_zzz_yyy_1[j] + 1.5 * fl1_fx * tey_zz_yyy_0[j] - 1.5 * fl1_fx * tey_zz_yyy_1[j];

                tez_zzzz_yyy_0[j] = pa_z[j] * tez_zzz_yyy_0[j] - pc_z[j] * tez_zzz_yyy_1[j] + 1.5 * fl1_fx * tez_zz_yyy_0[j] -
                                    1.5 * fl1_fx * tez_zz_yyy_1[j] + ta_zzz_yyy_1[j];

                tex_zzzz_yyz_0[j] = pa_z[j] * tex_zzz_yyz_0[j] - pc_z[j] * tex_zzz_yyz_1[j] + 1.5 * fl1_fx * tex_zz_yyz_0[j] -
                                    1.5 * fl1_fx * tex_zz_yyz_1[j] + 0.5 * fl1_fx * tex_zzz_yy_0[j] - 0.5 * fl1_fx * tex_zzz_yy_1[j];

                tey_zzzz_yyz_0[j] = pa_z[j] * tey_zzz_yyz_0[j] - pc_z[j] * tey_zzz_yyz_1[j] + 1.5 * fl1_fx * tey_zz_yyz_0[j] -
                                    1.5 * fl1_fx * tey_zz_yyz_1[j] + 0.5 * fl1_fx * tey_zzz_yy_0[j] - 0.5 * fl1_fx * tey_zzz_yy_1[j];

                tez_zzzz_yyz_0[j] = pa_z[j] * tez_zzz_yyz_0[j] - pc_z[j] * tez_zzz_yyz_1[j] + 1.5 * fl1_fx * tez_zz_yyz_0[j] -
                                    1.5 * fl1_fx * tez_zz_yyz_1[j] + 0.5 * fl1_fx * tez_zzz_yy_0[j] - 0.5 * fl1_fx * tez_zzz_yy_1[j] +
                                    ta_zzz_yyz_1[j];

                tex_zzzz_yzz_0[j] = pa_z[j] * tex_zzz_yzz_0[j] - pc_z[j] * tex_zzz_yzz_1[j] + 1.5 * fl1_fx * tex_zz_yzz_0[j] -
                                    1.5 * fl1_fx * tex_zz_yzz_1[j] + fl1_fx * tex_zzz_yz_0[j] - fl1_fx * tex_zzz_yz_1[j];

                tey_zzzz_yzz_0[j] = pa_z[j] * tey_zzz_yzz_0[j] - pc_z[j] * tey_zzz_yzz_1[j] + 1.5 * fl1_fx * tey_zz_yzz_0[j] -
                                    1.5 * fl1_fx * tey_zz_yzz_1[j] + fl1_fx * tey_zzz_yz_0[j] - fl1_fx * tey_zzz_yz_1[j];

                tez_zzzz_yzz_0[j] = pa_z[j] * tez_zzz_yzz_0[j] - pc_z[j] * tez_zzz_yzz_1[j] + 1.5 * fl1_fx * tez_zz_yzz_0[j] -
                                    1.5 * fl1_fx * tez_zz_yzz_1[j] + fl1_fx * tez_zzz_yz_0[j] - fl1_fx * tez_zzz_yz_1[j] + ta_zzz_yzz_1[j];

                tex_zzzz_zzz_0[j] = pa_z[j] * tex_zzz_zzz_0[j] - pc_z[j] * tex_zzz_zzz_1[j] + 1.5 * fl1_fx * tex_zz_zzz_0[j] -
                                    1.5 * fl1_fx * tex_zz_zzz_1[j] + 1.5 * fl1_fx * tex_zzz_zz_0[j] - 1.5 * fl1_fx * tex_zzz_zz_1[j];

                tey_zzzz_zzz_0[j] = pa_z[j] * tey_zzz_zzz_0[j] - pc_z[j] * tey_zzz_zzz_1[j] + 1.5 * fl1_fx * tey_zz_zzz_0[j] -
                                    1.5 * fl1_fx * tey_zz_zzz_1[j] + 1.5 * fl1_fx * tey_zzz_zz_0[j] - 1.5 * fl1_fx * tey_zzz_zz_1[j];

                tez_zzzz_zzz_0[j] = pa_z[j] * tez_zzz_zzz_0[j] - pc_z[j] * tez_zzz_zzz_1[j] + 1.5 * fl1_fx * tez_zz_zzz_0[j] -
                                    1.5 * fl1_fx * tez_zz_zzz_1[j] + 1.5 * fl1_fx * tez_zzz_zz_0[j] - 1.5 * fl1_fx * tez_zzz_zz_1[j] +
                                    ta_zzz_zzz_1[j];
            }

            idx++;
        }
    }
}

}  // namespace efieldrecfunc
