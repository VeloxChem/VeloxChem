//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "NuclearPotentialRecFuncForFF.hpp"

namespace npotrecfunc {  // npotrecfunc namespace

void
compNuclearPotentialForFF(CMemBlock2D<double>&       primBuffer,
                          const CRecursionMap&       recursionMap,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pcDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto)
{
    npotrecfunc::compNuclearPotentialForFF_0_50(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForFF_50_100(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compNuclearPotentialForFF_0_50(CMemBlock2D<double>&       primBuffer,
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

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_3_3_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_3_3_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_2_3_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_2_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_1_3_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_1_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_2_2_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_2_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

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

            auto ta_xx_xxx_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx);

            auto ta_xx_xxy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 1);

            auto ta_xx_xxz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 2);

            auto ta_xx_xyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 3);

            auto ta_xx_xyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 4);

            auto ta_xx_xzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 5);

            auto ta_xx_yyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 6);

            auto ta_xx_yyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 7);

            auto ta_xx_yzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 8);

            auto ta_xx_zzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 9);

            auto ta_xy_xxx_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 10);

            auto ta_xy_xxy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 11);

            auto ta_xy_xxz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 12);

            auto ta_xy_xyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 13);

            auto ta_xy_xyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 14);

            auto ta_xy_xzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 15);

            auto ta_xy_yyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 16);

            auto ta_xy_yyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 17);

            auto ta_xy_yzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 18);

            auto ta_xy_zzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 19);

            auto ta_xz_xxx_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 20);

            auto ta_xz_xxy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 21);

            auto ta_xz_xxz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 22);

            auto ta_xz_xyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 23);

            auto ta_xz_xyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 24);

            auto ta_xz_xzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 25);

            auto ta_xz_yyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 26);

            auto ta_xz_yyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 27);

            auto ta_xz_yzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 28);

            auto ta_xz_zzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 29);

            auto ta_yy_xxx_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 30);

            auto ta_yy_xxy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 31);

            auto ta_yy_xxz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 32);

            auto ta_yy_xyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 33);

            auto ta_yy_xyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 34);

            auto ta_yy_xzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 35);

            auto ta_yy_yyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 36);

            auto ta_yy_yyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 37);

            auto ta_yy_yzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 38);

            auto ta_yy_zzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 39);

            auto ta_yz_xxx_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 40);

            auto ta_yz_xxy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 41);

            auto ta_yz_xxz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 42);

            auto ta_yz_xyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 43);

            auto ta_yz_xyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 44);

            auto ta_yz_xzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 45);

            auto ta_yz_yyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 46);

            auto ta_yz_yyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 47);

            auto ta_yz_yzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 48);

            auto ta_yz_zzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 49);

            auto ta_xx_xxx_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx);

            auto ta_xx_xxy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 1);

            auto ta_xx_xxz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 2);

            auto ta_xx_xyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 3);

            auto ta_xx_xyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 4);

            auto ta_xx_xzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 5);

            auto ta_xx_yyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 6);

            auto ta_xx_yyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 7);

            auto ta_xx_yzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 8);

            auto ta_xx_zzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 9);

            auto ta_xy_xxx_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 10);

            auto ta_xy_xxy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 11);

            auto ta_xy_xxz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 12);

            auto ta_xy_xyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 13);

            auto ta_xy_xyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 14);

            auto ta_xy_xzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 15);

            auto ta_xy_yyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 16);

            auto ta_xy_yyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 17);

            auto ta_xy_yzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 18);

            auto ta_xy_zzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 19);

            auto ta_xz_xxx_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 20);

            auto ta_xz_xxy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 21);

            auto ta_xz_xxz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 22);

            auto ta_xz_xyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 23);

            auto ta_xz_xyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 24);

            auto ta_xz_xzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 25);

            auto ta_xz_yyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 26);

            auto ta_xz_yyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 27);

            auto ta_xz_yzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 28);

            auto ta_xz_zzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 29);

            auto ta_yy_xxx_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 30);

            auto ta_yy_xxy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 31);

            auto ta_yy_xxz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 32);

            auto ta_yy_xyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 33);

            auto ta_yy_xyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 34);

            auto ta_yy_xzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 35);

            auto ta_yy_yyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 36);

            auto ta_yy_yyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 37);

            auto ta_yy_yzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 38);

            auto ta_yy_zzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 39);

            auto ta_yz_xxx_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 40);

            auto ta_yz_xxy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 41);

            auto ta_yz_xxz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 42);

            auto ta_yz_xyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 43);

            auto ta_yz_xyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 44);

            auto ta_yz_xzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 45);

            auto ta_yz_yyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 46);

            auto ta_yz_yyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 47);

            auto ta_yz_yzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 48);

            auto ta_yz_zzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 49);

            auto ta_x_xxx_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx);

            auto ta_x_xxy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 1);

            auto ta_x_xxz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 2);

            auto ta_x_xyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 3);

            auto ta_x_xyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 4);

            auto ta_x_xzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 5);

            auto ta_x_yyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 6);

            auto ta_x_yyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 7);

            auto ta_x_yzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 8);

            auto ta_x_zzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 9);

            auto ta_y_xxx_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 10);

            auto ta_y_xxy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 11);

            auto ta_y_xxz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 12);

            auto ta_y_xyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 13);

            auto ta_y_xyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 14);

            auto ta_y_xzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 15);

            auto ta_y_yyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 16);

            auto ta_y_yyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 17);

            auto ta_y_yzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 18);

            auto ta_y_zzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 19);

            auto ta_z_xxx_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 20);

            auto ta_z_xxy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 21);

            auto ta_z_xxz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 22);

            auto ta_z_xyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 23);

            auto ta_z_xyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 24);

            auto ta_z_xzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 25);

            auto ta_z_yyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 26);

            auto ta_z_yyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 27);

            auto ta_z_yzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 28);

            auto ta_z_zzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 29);

            auto ta_x_xxx_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx);

            auto ta_x_xxy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 1);

            auto ta_x_xxz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 2);

            auto ta_x_xyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 3);

            auto ta_x_xyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 4);

            auto ta_x_xzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 5);

            auto ta_x_yyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 6);

            auto ta_x_yyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 7);

            auto ta_x_yzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 8);

            auto ta_x_zzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 9);

            auto ta_y_xxx_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 10);

            auto ta_y_xxy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 11);

            auto ta_y_xxz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 12);

            auto ta_y_xyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 13);

            auto ta_y_xyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 14);

            auto ta_y_xzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 15);

            auto ta_y_yyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 16);

            auto ta_y_yyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 17);

            auto ta_y_yzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 18);

            auto ta_y_zzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 19);

            auto ta_z_xxx_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 20);

            auto ta_z_xxy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 21);

            auto ta_z_xxz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 22);

            auto ta_z_xyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 23);

            auto ta_z_xyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 24);

            auto ta_z_xzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 25);

            auto ta_z_yyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 26);

            auto ta_z_yyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 27);

            auto ta_z_yzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 28);

            auto ta_z_zzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 29);

            auto ta_xx_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx);

            auto ta_xx_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 1);

            auto ta_xx_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 2);

            auto ta_xx_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 3);

            auto ta_xx_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 4);

            auto ta_xx_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 5);

            auto ta_xy_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 6);

            auto ta_xy_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 7);

            auto ta_xy_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 8);

            auto ta_xy_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 9);

            auto ta_xy_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 10);

            auto ta_xy_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 11);

            auto ta_xz_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 12);

            auto ta_xz_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 13);

            auto ta_xz_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 14);

            auto ta_xz_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 15);

            auto ta_xz_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 16);

            auto ta_xz_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 17);

            auto ta_yy_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 18);

            auto ta_yy_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 19);

            auto ta_yy_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 20);

            auto ta_yy_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 21);

            auto ta_yy_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 22);

            auto ta_yy_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 23);

            auto ta_yz_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 24);

            auto ta_yz_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 25);

            auto ta_yz_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 26);

            auto ta_yz_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 27);

            auto ta_yz_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 28);

            auto ta_yz_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 29);

            auto ta_xx_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx);

            auto ta_xx_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 1);

            auto ta_xx_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 2);

            auto ta_xx_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 3);

            auto ta_xx_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 4);

            auto ta_xx_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 5);

            auto ta_xy_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 6);

            auto ta_xy_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 7);

            auto ta_xy_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 8);

            auto ta_xy_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 9);

            auto ta_xy_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 10);

            auto ta_xy_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 11);

            auto ta_xz_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 12);

            auto ta_xz_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 13);

            auto ta_xz_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 14);

            auto ta_xz_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 15);

            auto ta_xz_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 16);

            auto ta_xz_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 17);

            auto ta_yy_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 18);

            auto ta_yy_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 19);

            auto ta_yy_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 20);

            auto ta_yy_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 21);

            auto ta_yy_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 22);

            auto ta_yy_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 23);

            auto ta_yz_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 24);

            auto ta_yz_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 25);

            auto ta_yz_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 26);

            auto ta_yz_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 27);

            auto ta_yz_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 28);

            auto ta_yz_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 29);

            // set up pointers to integrals

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

            // Batch of Integrals (0,50)

            #pragma omp simd aligned(fx, pa_x, pc_x, ta_x_xxx_0, ta_x_xxx_1, ta_x_xxy_0, ta_x_xxy_1, ta_x_xxz_0, \
                                         ta_x_xxz_1, ta_x_xyy_0, ta_x_xyy_1, ta_x_xyz_0, ta_x_xyz_1, ta_x_xzz_0, ta_x_xzz_1, \
                                         ta_x_yyy_0, ta_x_yyy_1, ta_x_yyz_0, ta_x_yyz_1, ta_x_yzz_0, ta_x_yzz_1, ta_x_zzz_0, \
                                         ta_x_zzz_1, ta_xx_xx_0, ta_xx_xx_1, ta_xx_xxx_0, ta_xx_xxx_1, ta_xx_xxy_0, \
                                         ta_xx_xxy_1, ta_xx_xxz_0, ta_xx_xxz_1, ta_xx_xy_0, ta_xx_xy_1, ta_xx_xyy_0, \
                                         ta_xx_xyy_1, ta_xx_xyz_0, ta_xx_xyz_1, ta_xx_xz_0, ta_xx_xz_1, ta_xx_xzz_0, \
                                         ta_xx_xzz_1, ta_xx_yy_0, ta_xx_yy_1, ta_xx_yyy_0, ta_xx_yyy_1, ta_xx_yyz_0, \
                                         ta_xx_yyz_1, ta_xx_yz_0, ta_xx_yz_1, ta_xx_yzz_0, ta_xx_yzz_1, ta_xx_zz_0, \
                                         ta_xx_zz_1, ta_xx_zzz_0, ta_xx_zzz_1, ta_xxx_xxx_0, ta_xxx_xxy_0, ta_xxx_xxz_0, \
                                         ta_xxx_xyy_0, ta_xxx_xyz_0, ta_xxx_xzz_0, ta_xxx_yyy_0, ta_xxx_yyz_0, ta_xxx_yzz_0, \
                                         ta_xxx_zzz_0, ta_xxy_xxx_0, ta_xxy_xxy_0, ta_xxy_xxz_0, ta_xxy_xyy_0, ta_xxy_xyz_0, \
                                         ta_xxy_xzz_0, ta_xxy_yyy_0, ta_xxy_yyz_0, ta_xxy_yzz_0, ta_xxy_zzz_0, ta_xxz_xxx_0, \
                                         ta_xxz_xxy_0, ta_xxz_xxz_0, ta_xxz_xyy_0, ta_xxz_xyz_0, ta_xxz_xzz_0, ta_xxz_yyy_0, \
                                         ta_xxz_yyz_0, ta_xxz_yzz_0, ta_xxz_zzz_0, ta_xy_xx_0, ta_xy_xx_1, ta_xy_xxx_0, \
                                         ta_xy_xxx_1, ta_xy_xxy_0, ta_xy_xxy_1, ta_xy_xxz_0, ta_xy_xxz_1, ta_xy_xy_0, \
                                         ta_xy_xy_1, ta_xy_xyy_0, ta_xy_xyy_1, ta_xy_xyz_0, ta_xy_xyz_1, ta_xy_xz_0, \
                                         ta_xy_xz_1, ta_xy_xzz_0, ta_xy_xzz_1, ta_xy_yy_0, ta_xy_yy_1, ta_xy_yyy_0, \
                                         ta_xy_yyy_1, ta_xy_yyz_0, ta_xy_yyz_1, ta_xy_yz_0, ta_xy_yz_1, ta_xy_yzz_0, \
                                         ta_xy_yzz_1, ta_xy_zz_0, ta_xy_zz_1, ta_xy_zzz_0, ta_xy_zzz_1, ta_xyy_xxx_0, \
                                         ta_xyy_xxy_0, ta_xyy_xxz_0, ta_xyy_xyy_0, ta_xyy_xyz_0, ta_xyy_xzz_0, ta_xyy_yyy_0, \
                                         ta_xyy_yyz_0, ta_xyy_yzz_0, ta_xyy_zzz_0, ta_xyz_xxx_0, ta_xyz_xxy_0, ta_xyz_xxz_0, \
                                         ta_xyz_xyy_0, ta_xyz_xyz_0, ta_xyz_xzz_0, ta_xyz_yyy_0, ta_xyz_yyz_0, ta_xyz_yzz_0, \
                                         ta_xyz_zzz_0, ta_xz_xx_0, ta_xz_xx_1, ta_xz_xxx_0, ta_xz_xxx_1, ta_xz_xxy_0, \
                                         ta_xz_xxy_1, ta_xz_xxz_0, ta_xz_xxz_1, ta_xz_xy_0, ta_xz_xy_1, ta_xz_xyy_0, \
                                         ta_xz_xyy_1, ta_xz_xyz_0, ta_xz_xyz_1, ta_xz_xz_0, ta_xz_xz_1, ta_xz_xzz_0, \
                                         ta_xz_xzz_1, ta_xz_yy_0, ta_xz_yy_1, ta_xz_yyy_0, ta_xz_yyy_1, ta_xz_yyz_0, \
                                         ta_xz_yyz_1, ta_xz_yz_0, ta_xz_yz_1, ta_xz_yzz_0, ta_xz_yzz_1, ta_xz_zz_0, \
                                         ta_xz_zz_1, ta_xz_zzz_0, ta_xz_zzz_1, ta_y_xxx_0, ta_y_xxx_1, ta_y_xxy_0, \
                                         ta_y_xxy_1, ta_y_xxz_0, ta_y_xxz_1, ta_y_xyy_0, ta_y_xyy_1, ta_y_xyz_0, ta_y_xyz_1, \
                                         ta_y_xzz_0, ta_y_xzz_1, ta_y_yyy_0, ta_y_yyy_1, ta_y_yyz_0, ta_y_yyz_1, ta_y_yzz_0, \
                                         ta_y_yzz_1, ta_y_zzz_0, ta_y_zzz_1, ta_yy_xx_0, ta_yy_xx_1, ta_yy_xxx_0, \
                                         ta_yy_xxx_1, ta_yy_xxy_0, ta_yy_xxy_1, ta_yy_xxz_0, ta_yy_xxz_1, ta_yy_xy_0, \
                                         ta_yy_xy_1, ta_yy_xyy_0, ta_yy_xyy_1, ta_yy_xyz_0, ta_yy_xyz_1, ta_yy_xz_0, \
                                         ta_yy_xz_1, ta_yy_xzz_0, ta_yy_xzz_1, ta_yy_yy_0, ta_yy_yy_1, ta_yy_yyy_0, \
                                         ta_yy_yyy_1, ta_yy_yyz_0, ta_yy_yyz_1, ta_yy_yz_0, ta_yy_yz_1, ta_yy_yzz_0, \
                                         ta_yy_yzz_1, ta_yy_zz_0, ta_yy_zz_1, ta_yy_zzz_0, ta_yy_zzz_1, ta_yz_xx_0, \
                                         ta_yz_xx_1, ta_yz_xxx_0, ta_yz_xxx_1, ta_yz_xxy_0, ta_yz_xxy_1, ta_yz_xxz_0, \
                                         ta_yz_xxz_1, ta_yz_xy_0, ta_yz_xy_1, ta_yz_xyy_0, ta_yz_xyy_1, ta_yz_xyz_0, \
                                         ta_yz_xyz_1, ta_yz_xz_0, ta_yz_xz_1, ta_yz_xzz_0, ta_yz_xzz_1, ta_yz_yy_0, \
                                         ta_yz_yy_1, ta_yz_yyy_0, ta_yz_yyy_1, ta_yz_yyz_0, ta_yz_yyz_1, ta_yz_yz_0, \
                                         ta_yz_yz_1, ta_yz_yzz_0, ta_yz_yzz_1, ta_yz_zz_0, ta_yz_zz_1, ta_yz_zzz_0, \
                                         ta_yz_zzz_1, ta_z_xxx_0, ta_z_xxx_1, ta_z_xxy_0, ta_z_xxy_1, ta_z_xxz_0, ta_z_xxz_1, \
                                         ta_z_xyy_0, ta_z_xyy_1, ta_z_xyz_0, ta_z_xyz_1, ta_z_xzz_0, ta_z_xzz_1, ta_z_yyy_0, \
                                         ta_z_yyy_1, ta_z_yyz_0, ta_z_yyz_1, ta_z_yzz_0, ta_z_yzz_1, ta_z_zzz_0, ta_z_zzz_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_xxx_xxx_0[j] = pa_x[j] * ta_xx_xxx_0[j] - pc_x[j] * ta_xx_xxx_1[j] + fl1_fx * ta_x_xxx_0[j] - fl1_fx * ta_x_xxx_1[j] +
                                  1.5 * fl1_fx * ta_xx_xx_0[j] - 1.5 * fl1_fx * ta_xx_xx_1[j];

                ta_xxx_xxy_0[j] = pa_x[j] * ta_xx_xxy_0[j] - pc_x[j] * ta_xx_xxy_1[j] + fl1_fx * ta_x_xxy_0[j] - fl1_fx * ta_x_xxy_1[j] +
                                  fl1_fx * ta_xx_xy_0[j] - fl1_fx * ta_xx_xy_1[j];

                ta_xxx_xxz_0[j] = pa_x[j] * ta_xx_xxz_0[j] - pc_x[j] * ta_xx_xxz_1[j] + fl1_fx * ta_x_xxz_0[j] - fl1_fx * ta_x_xxz_1[j] +
                                  fl1_fx * ta_xx_xz_0[j] - fl1_fx * ta_xx_xz_1[j];

                ta_xxx_xyy_0[j] = pa_x[j] * ta_xx_xyy_0[j] - pc_x[j] * ta_xx_xyy_1[j] + fl1_fx * ta_x_xyy_0[j] - fl1_fx * ta_x_xyy_1[j] +
                                  0.5 * fl1_fx * ta_xx_yy_0[j] - 0.5 * fl1_fx * ta_xx_yy_1[j];

                ta_xxx_xyz_0[j] = pa_x[j] * ta_xx_xyz_0[j] - pc_x[j] * ta_xx_xyz_1[j] + fl1_fx * ta_x_xyz_0[j] - fl1_fx * ta_x_xyz_1[j] +
                                  0.5 * fl1_fx * ta_xx_yz_0[j] - 0.5 * fl1_fx * ta_xx_yz_1[j];

                ta_xxx_xzz_0[j] = pa_x[j] * ta_xx_xzz_0[j] - pc_x[j] * ta_xx_xzz_1[j] + fl1_fx * ta_x_xzz_0[j] - fl1_fx * ta_x_xzz_1[j] +
                                  0.5 * fl1_fx * ta_xx_zz_0[j] - 0.5 * fl1_fx * ta_xx_zz_1[j];

                ta_xxx_yyy_0[j] = pa_x[j] * ta_xx_yyy_0[j] - pc_x[j] * ta_xx_yyy_1[j] + fl1_fx * ta_x_yyy_0[j] - fl1_fx * ta_x_yyy_1[j];

                ta_xxx_yyz_0[j] = pa_x[j] * ta_xx_yyz_0[j] - pc_x[j] * ta_xx_yyz_1[j] + fl1_fx * ta_x_yyz_0[j] - fl1_fx * ta_x_yyz_1[j];

                ta_xxx_yzz_0[j] = pa_x[j] * ta_xx_yzz_0[j] - pc_x[j] * ta_xx_yzz_1[j] + fl1_fx * ta_x_yzz_0[j] - fl1_fx * ta_x_yzz_1[j];

                ta_xxx_zzz_0[j] = pa_x[j] * ta_xx_zzz_0[j] - pc_x[j] * ta_xx_zzz_1[j] + fl1_fx * ta_x_zzz_0[j] - fl1_fx * ta_x_zzz_1[j];

                ta_xxy_xxx_0[j] = pa_x[j] * ta_xy_xxx_0[j] - pc_x[j] * ta_xy_xxx_1[j] + 0.5 * fl1_fx * ta_y_xxx_0[j] - 0.5 * fl1_fx * ta_y_xxx_1[j] +
                                  1.5 * fl1_fx * ta_xy_xx_0[j] - 1.5 * fl1_fx * ta_xy_xx_1[j];

                ta_xxy_xxy_0[j] = pa_x[j] * ta_xy_xxy_0[j] - pc_x[j] * ta_xy_xxy_1[j] + 0.5 * fl1_fx * ta_y_xxy_0[j] - 0.5 * fl1_fx * ta_y_xxy_1[j] +
                                  fl1_fx * ta_xy_xy_0[j] - fl1_fx * ta_xy_xy_1[j];

                ta_xxy_xxz_0[j] = pa_x[j] * ta_xy_xxz_0[j] - pc_x[j] * ta_xy_xxz_1[j] + 0.5 * fl1_fx * ta_y_xxz_0[j] - 0.5 * fl1_fx * ta_y_xxz_1[j] +
                                  fl1_fx * ta_xy_xz_0[j] - fl1_fx * ta_xy_xz_1[j];

                ta_xxy_xyy_0[j] = pa_x[j] * ta_xy_xyy_0[j] - pc_x[j] * ta_xy_xyy_1[j] + 0.5 * fl1_fx * ta_y_xyy_0[j] - 0.5 * fl1_fx * ta_y_xyy_1[j] +
                                  0.5 * fl1_fx * ta_xy_yy_0[j] - 0.5 * fl1_fx * ta_xy_yy_1[j];

                ta_xxy_xyz_0[j] = pa_x[j] * ta_xy_xyz_0[j] - pc_x[j] * ta_xy_xyz_1[j] + 0.5 * fl1_fx * ta_y_xyz_0[j] - 0.5 * fl1_fx * ta_y_xyz_1[j] +
                                  0.5 * fl1_fx * ta_xy_yz_0[j] - 0.5 * fl1_fx * ta_xy_yz_1[j];

                ta_xxy_xzz_0[j] = pa_x[j] * ta_xy_xzz_0[j] - pc_x[j] * ta_xy_xzz_1[j] + 0.5 * fl1_fx * ta_y_xzz_0[j] - 0.5 * fl1_fx * ta_y_xzz_1[j] +
                                  0.5 * fl1_fx * ta_xy_zz_0[j] - 0.5 * fl1_fx * ta_xy_zz_1[j];

                ta_xxy_yyy_0[j] = pa_x[j] * ta_xy_yyy_0[j] - pc_x[j] * ta_xy_yyy_1[j] + 0.5 * fl1_fx * ta_y_yyy_0[j] - 0.5 * fl1_fx * ta_y_yyy_1[j];

                ta_xxy_yyz_0[j] = pa_x[j] * ta_xy_yyz_0[j] - pc_x[j] * ta_xy_yyz_1[j] + 0.5 * fl1_fx * ta_y_yyz_0[j] - 0.5 * fl1_fx * ta_y_yyz_1[j];

                ta_xxy_yzz_0[j] = pa_x[j] * ta_xy_yzz_0[j] - pc_x[j] * ta_xy_yzz_1[j] + 0.5 * fl1_fx * ta_y_yzz_0[j] - 0.5 * fl1_fx * ta_y_yzz_1[j];

                ta_xxy_zzz_0[j] = pa_x[j] * ta_xy_zzz_0[j] - pc_x[j] * ta_xy_zzz_1[j] + 0.5 * fl1_fx * ta_y_zzz_0[j] - 0.5 * fl1_fx * ta_y_zzz_1[j];

                ta_xxz_xxx_0[j] = pa_x[j] * ta_xz_xxx_0[j] - pc_x[j] * ta_xz_xxx_1[j] + 0.5 * fl1_fx * ta_z_xxx_0[j] - 0.5 * fl1_fx * ta_z_xxx_1[j] +
                                  1.5 * fl1_fx * ta_xz_xx_0[j] - 1.5 * fl1_fx * ta_xz_xx_1[j];

                ta_xxz_xxy_0[j] = pa_x[j] * ta_xz_xxy_0[j] - pc_x[j] * ta_xz_xxy_1[j] + 0.5 * fl1_fx * ta_z_xxy_0[j] - 0.5 * fl1_fx * ta_z_xxy_1[j] +
                                  fl1_fx * ta_xz_xy_0[j] - fl1_fx * ta_xz_xy_1[j];

                ta_xxz_xxz_0[j] = pa_x[j] * ta_xz_xxz_0[j] - pc_x[j] * ta_xz_xxz_1[j] + 0.5 * fl1_fx * ta_z_xxz_0[j] - 0.5 * fl1_fx * ta_z_xxz_1[j] +
                                  fl1_fx * ta_xz_xz_0[j] - fl1_fx * ta_xz_xz_1[j];

                ta_xxz_xyy_0[j] = pa_x[j] * ta_xz_xyy_0[j] - pc_x[j] * ta_xz_xyy_1[j] + 0.5 * fl1_fx * ta_z_xyy_0[j] - 0.5 * fl1_fx * ta_z_xyy_1[j] +
                                  0.5 * fl1_fx * ta_xz_yy_0[j] - 0.5 * fl1_fx * ta_xz_yy_1[j];

                ta_xxz_xyz_0[j] = pa_x[j] * ta_xz_xyz_0[j] - pc_x[j] * ta_xz_xyz_1[j] + 0.5 * fl1_fx * ta_z_xyz_0[j] - 0.5 * fl1_fx * ta_z_xyz_1[j] +
                                  0.5 * fl1_fx * ta_xz_yz_0[j] - 0.5 * fl1_fx * ta_xz_yz_1[j];

                ta_xxz_xzz_0[j] = pa_x[j] * ta_xz_xzz_0[j] - pc_x[j] * ta_xz_xzz_1[j] + 0.5 * fl1_fx * ta_z_xzz_0[j] - 0.5 * fl1_fx * ta_z_xzz_1[j] +
                                  0.5 * fl1_fx * ta_xz_zz_0[j] - 0.5 * fl1_fx * ta_xz_zz_1[j];

                ta_xxz_yyy_0[j] = pa_x[j] * ta_xz_yyy_0[j] - pc_x[j] * ta_xz_yyy_1[j] + 0.5 * fl1_fx * ta_z_yyy_0[j] - 0.5 * fl1_fx * ta_z_yyy_1[j];

                ta_xxz_yyz_0[j] = pa_x[j] * ta_xz_yyz_0[j] - pc_x[j] * ta_xz_yyz_1[j] + 0.5 * fl1_fx * ta_z_yyz_0[j] - 0.5 * fl1_fx * ta_z_yyz_1[j];

                ta_xxz_yzz_0[j] = pa_x[j] * ta_xz_yzz_0[j] - pc_x[j] * ta_xz_yzz_1[j] + 0.5 * fl1_fx * ta_z_yzz_0[j] - 0.5 * fl1_fx * ta_z_yzz_1[j];

                ta_xxz_zzz_0[j] = pa_x[j] * ta_xz_zzz_0[j] - pc_x[j] * ta_xz_zzz_1[j] + 0.5 * fl1_fx * ta_z_zzz_0[j] - 0.5 * fl1_fx * ta_z_zzz_1[j];

                ta_xyy_xxx_0[j] = pa_x[j] * ta_yy_xxx_0[j] - pc_x[j] * ta_yy_xxx_1[j] + 1.5 * fl1_fx * ta_yy_xx_0[j] - 1.5 * fl1_fx * ta_yy_xx_1[j];

                ta_xyy_xxy_0[j] = pa_x[j] * ta_yy_xxy_0[j] - pc_x[j] * ta_yy_xxy_1[j] + fl1_fx * ta_yy_xy_0[j] - fl1_fx * ta_yy_xy_1[j];

                ta_xyy_xxz_0[j] = pa_x[j] * ta_yy_xxz_0[j] - pc_x[j] * ta_yy_xxz_1[j] + fl1_fx * ta_yy_xz_0[j] - fl1_fx * ta_yy_xz_1[j];

                ta_xyy_xyy_0[j] = pa_x[j] * ta_yy_xyy_0[j] - pc_x[j] * ta_yy_xyy_1[j] + 0.5 * fl1_fx * ta_yy_yy_0[j] - 0.5 * fl1_fx * ta_yy_yy_1[j];

                ta_xyy_xyz_0[j] = pa_x[j] * ta_yy_xyz_0[j] - pc_x[j] * ta_yy_xyz_1[j] + 0.5 * fl1_fx * ta_yy_yz_0[j] - 0.5 * fl1_fx * ta_yy_yz_1[j];

                ta_xyy_xzz_0[j] = pa_x[j] * ta_yy_xzz_0[j] - pc_x[j] * ta_yy_xzz_1[j] + 0.5 * fl1_fx * ta_yy_zz_0[j] - 0.5 * fl1_fx * ta_yy_zz_1[j];

                ta_xyy_yyy_0[j] = pa_x[j] * ta_yy_yyy_0[j] - pc_x[j] * ta_yy_yyy_1[j];

                ta_xyy_yyz_0[j] = pa_x[j] * ta_yy_yyz_0[j] - pc_x[j] * ta_yy_yyz_1[j];

                ta_xyy_yzz_0[j] = pa_x[j] * ta_yy_yzz_0[j] - pc_x[j] * ta_yy_yzz_1[j];

                ta_xyy_zzz_0[j] = pa_x[j] * ta_yy_zzz_0[j] - pc_x[j] * ta_yy_zzz_1[j];

                ta_xyz_xxx_0[j] = pa_x[j] * ta_yz_xxx_0[j] - pc_x[j] * ta_yz_xxx_1[j] + 1.5 * fl1_fx * ta_yz_xx_0[j] - 1.5 * fl1_fx * ta_yz_xx_1[j];

                ta_xyz_xxy_0[j] = pa_x[j] * ta_yz_xxy_0[j] - pc_x[j] * ta_yz_xxy_1[j] + fl1_fx * ta_yz_xy_0[j] - fl1_fx * ta_yz_xy_1[j];

                ta_xyz_xxz_0[j] = pa_x[j] * ta_yz_xxz_0[j] - pc_x[j] * ta_yz_xxz_1[j] + fl1_fx * ta_yz_xz_0[j] - fl1_fx * ta_yz_xz_1[j];

                ta_xyz_xyy_0[j] = pa_x[j] * ta_yz_xyy_0[j] - pc_x[j] * ta_yz_xyy_1[j] + 0.5 * fl1_fx * ta_yz_yy_0[j] - 0.5 * fl1_fx * ta_yz_yy_1[j];

                ta_xyz_xyz_0[j] = pa_x[j] * ta_yz_xyz_0[j] - pc_x[j] * ta_yz_xyz_1[j] + 0.5 * fl1_fx * ta_yz_yz_0[j] - 0.5 * fl1_fx * ta_yz_yz_1[j];

                ta_xyz_xzz_0[j] = pa_x[j] * ta_yz_xzz_0[j] - pc_x[j] * ta_yz_xzz_1[j] + 0.5 * fl1_fx * ta_yz_zz_0[j] - 0.5 * fl1_fx * ta_yz_zz_1[j];

                ta_xyz_yyy_0[j] = pa_x[j] * ta_yz_yyy_0[j] - pc_x[j] * ta_yz_yyy_1[j];

                ta_xyz_yyz_0[j] = pa_x[j] * ta_yz_yyz_0[j] - pc_x[j] * ta_yz_yyz_1[j];

                ta_xyz_yzz_0[j] = pa_x[j] * ta_yz_yzz_0[j] - pc_x[j] * ta_yz_yzz_1[j];

                ta_xyz_zzz_0[j] = pa_x[j] * ta_yz_zzz_0[j] - pc_x[j] * ta_yz_zzz_1[j];
            }

            idx++;
        }
    }
}

void
compNuclearPotentialForFF_50_100(CMemBlock2D<double>&       primBuffer,
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

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_3_3_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_3_3_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_2_3_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_2_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_2_2_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_2_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_1_3_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_1_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(3 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to tensors product of distances R(PC) = P - C

            auto pc_x = pcDistances.data(3 * idx);

            auto pc_y = pcDistances.data(3 * idx + 1);

            auto pc_z = pcDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto ta_yy_xxx_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 30);

            auto ta_yy_xxy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 31);

            auto ta_yy_xxz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 32);

            auto ta_yy_xyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 33);

            auto ta_yy_xyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 34);

            auto ta_yy_xzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 35);

            auto ta_yy_yyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 36);

            auto ta_yy_yyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 37);

            auto ta_yy_yzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 38);

            auto ta_yy_zzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 39);

            auto ta_yz_xxx_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 40);

            auto ta_yz_xxy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 41);

            auto ta_yz_xxz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 42);

            auto ta_yz_xyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 43);

            auto ta_yz_xyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 44);

            auto ta_yz_xzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 45);

            auto ta_yz_yyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 46);

            auto ta_yz_yyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 47);

            auto ta_yz_yzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 48);

            auto ta_yz_zzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 49);

            auto ta_zz_xxx_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 50);

            auto ta_zz_xxy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 51);

            auto ta_zz_xxz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 52);

            auto ta_zz_xyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 53);

            auto ta_zz_xyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 54);

            auto ta_zz_xzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 55);

            auto ta_zz_yyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 56);

            auto ta_zz_yyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 57);

            auto ta_zz_yzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 58);

            auto ta_zz_zzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 59);

            auto ta_yy_xxx_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 30);

            auto ta_yy_xxy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 31);

            auto ta_yy_xxz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 32);

            auto ta_yy_xyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 33);

            auto ta_yy_xyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 34);

            auto ta_yy_xzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 35);

            auto ta_yy_yyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 36);

            auto ta_yy_yyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 37);

            auto ta_yy_yzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 38);

            auto ta_yy_zzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 39);

            auto ta_yz_xxx_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 40);

            auto ta_yz_xxy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 41);

            auto ta_yz_xxz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 42);

            auto ta_yz_xyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 43);

            auto ta_yz_xyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 44);

            auto ta_yz_xzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 45);

            auto ta_yz_yyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 46);

            auto ta_yz_yyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 47);

            auto ta_yz_yzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 48);

            auto ta_yz_zzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 49);

            auto ta_zz_xxx_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 50);

            auto ta_zz_xxy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 51);

            auto ta_zz_xxz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 52);

            auto ta_zz_xyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 53);

            auto ta_zz_xyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 54);

            auto ta_zz_xzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 55);

            auto ta_zz_yyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 56);

            auto ta_zz_yyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 57);

            auto ta_zz_yzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 58);

            auto ta_zz_zzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 59);

            auto ta_y_xxx_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 10);

            auto ta_y_xxy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 11);

            auto ta_y_xxz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 12);

            auto ta_y_xyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 13);

            auto ta_y_xyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 14);

            auto ta_y_xzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 15);

            auto ta_y_yyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 16);

            auto ta_y_yyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 17);

            auto ta_y_yzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 18);

            auto ta_y_zzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 19);

            auto ta_z_xxx_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 20);

            auto ta_z_xxy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 21);

            auto ta_z_xxz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 22);

            auto ta_z_xyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 23);

            auto ta_z_xyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 24);

            auto ta_z_xzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 25);

            auto ta_z_yyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 26);

            auto ta_z_yyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 27);

            auto ta_z_yzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 28);

            auto ta_z_zzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 29);

            auto ta_y_xxx_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 10);

            auto ta_y_xxy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 11);

            auto ta_y_xxz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 12);

            auto ta_y_xyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 13);

            auto ta_y_xyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 14);

            auto ta_y_xzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 15);

            auto ta_y_yyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 16);

            auto ta_y_yyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 17);

            auto ta_y_yzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 18);

            auto ta_y_zzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 19);

            auto ta_z_xxx_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 20);

            auto ta_z_xxy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 21);

            auto ta_z_xxz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 22);

            auto ta_z_xyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 23);

            auto ta_z_xyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 24);

            auto ta_z_xzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 25);

            auto ta_z_yyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 26);

            auto ta_z_yyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 27);

            auto ta_z_yzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 28);

            auto ta_z_zzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 29);

            auto ta_yy_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 18);

            auto ta_yy_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 19);

            auto ta_yy_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 20);

            auto ta_yy_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 21);

            auto ta_yy_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 22);

            auto ta_yy_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 23);

            auto ta_yz_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 24);

            auto ta_yz_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 25);

            auto ta_yz_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 26);

            auto ta_yz_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 27);

            auto ta_yz_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 28);

            auto ta_yz_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 29);

            auto ta_zz_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 30);

            auto ta_zz_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 31);

            auto ta_zz_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 32);

            auto ta_zz_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 33);

            auto ta_zz_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 34);

            auto ta_zz_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 35);

            auto ta_yy_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 18);

            auto ta_yy_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 19);

            auto ta_yy_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 20);

            auto ta_yy_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 21);

            auto ta_yy_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 22);

            auto ta_yy_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 23);

            auto ta_yz_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 24);

            auto ta_yz_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 25);

            auto ta_yz_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 26);

            auto ta_yz_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 27);

            auto ta_yz_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 28);

            auto ta_yz_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 29);

            auto ta_zz_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 30);

            auto ta_zz_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 31);

            auto ta_zz_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 32);

            auto ta_zz_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 33);

            auto ta_zz_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 34);

            auto ta_zz_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 35);

            // set up pointers to integrals

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

            // Batch of Integrals (50,100)

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_xzz_xxx_0, ta_xzz_xxy_0, \
                                         ta_xzz_xxz_0, ta_xzz_xyy_0, ta_xzz_xyz_0, ta_xzz_xzz_0, ta_xzz_yyy_0, ta_xzz_yyz_0, \
                                         ta_xzz_yzz_0, ta_xzz_zzz_0, ta_y_xxx_0, ta_y_xxx_1, ta_y_xxy_0, ta_y_xxy_1, \
                                         ta_y_xxz_0, ta_y_xxz_1, ta_y_xyy_0, ta_y_xyy_1, ta_y_xyz_0, ta_y_xyz_1, ta_y_xzz_0, \
                                         ta_y_xzz_1, ta_y_yyy_0, ta_y_yyy_1, ta_y_yyz_0, ta_y_yyz_1, ta_y_yzz_0, ta_y_yzz_1, \
                                         ta_y_zzz_0, ta_y_zzz_1, ta_yy_xx_0, ta_yy_xx_1, ta_yy_xxx_0, ta_yy_xxx_1, \
                                         ta_yy_xxy_0, ta_yy_xxy_1, ta_yy_xxz_0, ta_yy_xxz_1, ta_yy_xy_0, ta_yy_xy_1, \
                                         ta_yy_xyy_0, ta_yy_xyy_1, ta_yy_xyz_0, ta_yy_xyz_1, ta_yy_xz_0, ta_yy_xz_1, \
                                         ta_yy_xzz_0, ta_yy_xzz_1, ta_yy_yy_0, ta_yy_yy_1, ta_yy_yyy_0, ta_yy_yyy_1, \
                                         ta_yy_yyz_0, ta_yy_yyz_1, ta_yy_yz_0, ta_yy_yz_1, ta_yy_yzz_0, ta_yy_yzz_1, \
                                         ta_yy_zz_0, ta_yy_zz_1, ta_yy_zzz_0, ta_yy_zzz_1, ta_yyy_xxx_0, ta_yyy_xxy_0, \
                                         ta_yyy_xxz_0, ta_yyy_xyy_0, ta_yyy_xyz_0, ta_yyy_xzz_0, ta_yyy_yyy_0, ta_yyy_yyz_0, \
                                         ta_yyy_yzz_0, ta_yyy_zzz_0, ta_yyz_xxx_0, ta_yyz_xxy_0, ta_yyz_xxz_0, ta_yyz_xyy_0, \
                                         ta_yyz_xyz_0, ta_yyz_xzz_0, ta_yyz_yyy_0, ta_yyz_yyz_0, ta_yyz_yzz_0, ta_yyz_zzz_0, \
                                         ta_yz_xx_0, ta_yz_xx_1, ta_yz_xxx_0, ta_yz_xxx_1, ta_yz_xxy_0, ta_yz_xxy_1, \
                                         ta_yz_xxz_0, ta_yz_xxz_1, ta_yz_xy_0, ta_yz_xy_1, ta_yz_xyy_0, ta_yz_xyy_1, \
                                         ta_yz_xyz_0, ta_yz_xyz_1, ta_yz_xz_0, ta_yz_xz_1, ta_yz_xzz_0, ta_yz_xzz_1, \
                                         ta_yz_yy_0, ta_yz_yy_1, ta_yz_yyy_0, ta_yz_yyy_1, ta_yz_yyz_0, ta_yz_yyz_1, \
                                         ta_yz_yz_0, ta_yz_yz_1, ta_yz_yzz_0, ta_yz_yzz_1, ta_yz_zz_0, ta_yz_zz_1, \
                                         ta_yz_zzz_0, ta_yz_zzz_1, ta_yzz_xxx_0, ta_yzz_xxy_0, ta_yzz_xxz_0, ta_yzz_xyy_0, \
                                         ta_yzz_xyz_0, ta_yzz_xzz_0, ta_yzz_yyy_0, ta_yzz_yyz_0, ta_yzz_yzz_0, ta_yzz_zzz_0, \
                                         ta_z_xxx_0, ta_z_xxx_1, ta_z_xxy_0, ta_z_xxy_1, ta_z_xxz_0, ta_z_xxz_1, ta_z_xyy_0, \
                                         ta_z_xyy_1, ta_z_xyz_0, ta_z_xyz_1, ta_z_xzz_0, ta_z_xzz_1, ta_z_yyy_0, ta_z_yyy_1, \
                                         ta_z_yyz_0, ta_z_yyz_1, ta_z_yzz_0, ta_z_yzz_1, ta_z_zzz_0, ta_z_zzz_1, ta_zz_xx_0, \
                                         ta_zz_xx_1, ta_zz_xxx_0, ta_zz_xxx_1, ta_zz_xxy_0, ta_zz_xxy_1, ta_zz_xxz_0, \
                                         ta_zz_xxz_1, ta_zz_xy_0, ta_zz_xy_1, ta_zz_xyy_0, ta_zz_xyy_1, ta_zz_xyz_0, \
                                         ta_zz_xyz_1, ta_zz_xz_0, ta_zz_xz_1, ta_zz_xzz_0, ta_zz_xzz_1, ta_zz_yy_0, \
                                         ta_zz_yy_1, ta_zz_yyy_0, ta_zz_yyy_1, ta_zz_yyz_0, ta_zz_yyz_1, ta_zz_yz_0, \
                                         ta_zz_yz_1, ta_zz_yzz_0, ta_zz_yzz_1, ta_zz_zz_0, ta_zz_zz_1, ta_zz_zzz_0, \
                                         ta_zz_zzz_1, ta_zzz_xxx_0, ta_zzz_xxy_0, ta_zzz_xxz_0, ta_zzz_xyy_0, ta_zzz_xyz_0, \
                                         ta_zzz_xzz_0, ta_zzz_yyy_0, ta_zzz_yyz_0, ta_zzz_yzz_0, ta_zzz_zzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_xzz_xxx_0[j] = pa_x[j] * ta_zz_xxx_0[j] - pc_x[j] * ta_zz_xxx_1[j] + 1.5 * fl1_fx * ta_zz_xx_0[j] - 1.5 * fl1_fx * ta_zz_xx_1[j];

                ta_xzz_xxy_0[j] = pa_x[j] * ta_zz_xxy_0[j] - pc_x[j] * ta_zz_xxy_1[j] + fl1_fx * ta_zz_xy_0[j] - fl1_fx * ta_zz_xy_1[j];

                ta_xzz_xxz_0[j] = pa_x[j] * ta_zz_xxz_0[j] - pc_x[j] * ta_zz_xxz_1[j] + fl1_fx * ta_zz_xz_0[j] - fl1_fx * ta_zz_xz_1[j];

                ta_xzz_xyy_0[j] = pa_x[j] * ta_zz_xyy_0[j] - pc_x[j] * ta_zz_xyy_1[j] + 0.5 * fl1_fx * ta_zz_yy_0[j] - 0.5 * fl1_fx * ta_zz_yy_1[j];

                ta_xzz_xyz_0[j] = pa_x[j] * ta_zz_xyz_0[j] - pc_x[j] * ta_zz_xyz_1[j] + 0.5 * fl1_fx * ta_zz_yz_0[j] - 0.5 * fl1_fx * ta_zz_yz_1[j];

                ta_xzz_xzz_0[j] = pa_x[j] * ta_zz_xzz_0[j] - pc_x[j] * ta_zz_xzz_1[j] + 0.5 * fl1_fx * ta_zz_zz_0[j] - 0.5 * fl1_fx * ta_zz_zz_1[j];

                ta_xzz_yyy_0[j] = pa_x[j] * ta_zz_yyy_0[j] - pc_x[j] * ta_zz_yyy_1[j];

                ta_xzz_yyz_0[j] = pa_x[j] * ta_zz_yyz_0[j] - pc_x[j] * ta_zz_yyz_1[j];

                ta_xzz_yzz_0[j] = pa_x[j] * ta_zz_yzz_0[j] - pc_x[j] * ta_zz_yzz_1[j];

                ta_xzz_zzz_0[j] = pa_x[j] * ta_zz_zzz_0[j] - pc_x[j] * ta_zz_zzz_1[j];

                ta_yyy_xxx_0[j] = pa_y[j] * ta_yy_xxx_0[j] - pc_y[j] * ta_yy_xxx_1[j] + fl1_fx * ta_y_xxx_0[j] - fl1_fx * ta_y_xxx_1[j];

                ta_yyy_xxy_0[j] = pa_y[j] * ta_yy_xxy_0[j] - pc_y[j] * ta_yy_xxy_1[j] + fl1_fx * ta_y_xxy_0[j] - fl1_fx * ta_y_xxy_1[j] +
                                  0.5 * fl1_fx * ta_yy_xx_0[j] - 0.5 * fl1_fx * ta_yy_xx_1[j];

                ta_yyy_xxz_0[j] = pa_y[j] * ta_yy_xxz_0[j] - pc_y[j] * ta_yy_xxz_1[j] + fl1_fx * ta_y_xxz_0[j] - fl1_fx * ta_y_xxz_1[j];

                ta_yyy_xyy_0[j] = pa_y[j] * ta_yy_xyy_0[j] - pc_y[j] * ta_yy_xyy_1[j] + fl1_fx * ta_y_xyy_0[j] - fl1_fx * ta_y_xyy_1[j] +
                                  fl1_fx * ta_yy_xy_0[j] - fl1_fx * ta_yy_xy_1[j];

                ta_yyy_xyz_0[j] = pa_y[j] * ta_yy_xyz_0[j] - pc_y[j] * ta_yy_xyz_1[j] + fl1_fx * ta_y_xyz_0[j] - fl1_fx * ta_y_xyz_1[j] +
                                  0.5 * fl1_fx * ta_yy_xz_0[j] - 0.5 * fl1_fx * ta_yy_xz_1[j];

                ta_yyy_xzz_0[j] = pa_y[j] * ta_yy_xzz_0[j] - pc_y[j] * ta_yy_xzz_1[j] + fl1_fx * ta_y_xzz_0[j] - fl1_fx * ta_y_xzz_1[j];

                ta_yyy_yyy_0[j] = pa_y[j] * ta_yy_yyy_0[j] - pc_y[j] * ta_yy_yyy_1[j] + fl1_fx * ta_y_yyy_0[j] - fl1_fx * ta_y_yyy_1[j] +
                                  1.5 * fl1_fx * ta_yy_yy_0[j] - 1.5 * fl1_fx * ta_yy_yy_1[j];

                ta_yyy_yyz_0[j] = pa_y[j] * ta_yy_yyz_0[j] - pc_y[j] * ta_yy_yyz_1[j] + fl1_fx * ta_y_yyz_0[j] - fl1_fx * ta_y_yyz_1[j] +
                                  fl1_fx * ta_yy_yz_0[j] - fl1_fx * ta_yy_yz_1[j];

                ta_yyy_yzz_0[j] = pa_y[j] * ta_yy_yzz_0[j] - pc_y[j] * ta_yy_yzz_1[j] + fl1_fx * ta_y_yzz_0[j] - fl1_fx * ta_y_yzz_1[j] +
                                  0.5 * fl1_fx * ta_yy_zz_0[j] - 0.5 * fl1_fx * ta_yy_zz_1[j];

                ta_yyy_zzz_0[j] = pa_y[j] * ta_yy_zzz_0[j] - pc_y[j] * ta_yy_zzz_1[j] + fl1_fx * ta_y_zzz_0[j] - fl1_fx * ta_y_zzz_1[j];

                ta_yyz_xxx_0[j] = pa_y[j] * ta_yz_xxx_0[j] - pc_y[j] * ta_yz_xxx_1[j] + 0.5 * fl1_fx * ta_z_xxx_0[j] - 0.5 * fl1_fx * ta_z_xxx_1[j];

                ta_yyz_xxy_0[j] = pa_y[j] * ta_yz_xxy_0[j] - pc_y[j] * ta_yz_xxy_1[j] + 0.5 * fl1_fx * ta_z_xxy_0[j] - 0.5 * fl1_fx * ta_z_xxy_1[j] +
                                  0.5 * fl1_fx * ta_yz_xx_0[j] - 0.5 * fl1_fx * ta_yz_xx_1[j];

                ta_yyz_xxz_0[j] = pa_y[j] * ta_yz_xxz_0[j] - pc_y[j] * ta_yz_xxz_1[j] + 0.5 * fl1_fx * ta_z_xxz_0[j] - 0.5 * fl1_fx * ta_z_xxz_1[j];

                ta_yyz_xyy_0[j] = pa_y[j] * ta_yz_xyy_0[j] - pc_y[j] * ta_yz_xyy_1[j] + 0.5 * fl1_fx * ta_z_xyy_0[j] - 0.5 * fl1_fx * ta_z_xyy_1[j] +
                                  fl1_fx * ta_yz_xy_0[j] - fl1_fx * ta_yz_xy_1[j];

                ta_yyz_xyz_0[j] = pa_y[j] * ta_yz_xyz_0[j] - pc_y[j] * ta_yz_xyz_1[j] + 0.5 * fl1_fx * ta_z_xyz_0[j] - 0.5 * fl1_fx * ta_z_xyz_1[j] +
                                  0.5 * fl1_fx * ta_yz_xz_0[j] - 0.5 * fl1_fx * ta_yz_xz_1[j];

                ta_yyz_xzz_0[j] = pa_y[j] * ta_yz_xzz_0[j] - pc_y[j] * ta_yz_xzz_1[j] + 0.5 * fl1_fx * ta_z_xzz_0[j] - 0.5 * fl1_fx * ta_z_xzz_1[j];

                ta_yyz_yyy_0[j] = pa_y[j] * ta_yz_yyy_0[j] - pc_y[j] * ta_yz_yyy_1[j] + 0.5 * fl1_fx * ta_z_yyy_0[j] - 0.5 * fl1_fx * ta_z_yyy_1[j] +
                                  1.5 * fl1_fx * ta_yz_yy_0[j] - 1.5 * fl1_fx * ta_yz_yy_1[j];

                ta_yyz_yyz_0[j] = pa_y[j] * ta_yz_yyz_0[j] - pc_y[j] * ta_yz_yyz_1[j] + 0.5 * fl1_fx * ta_z_yyz_0[j] - 0.5 * fl1_fx * ta_z_yyz_1[j] +
                                  fl1_fx * ta_yz_yz_0[j] - fl1_fx * ta_yz_yz_1[j];

                ta_yyz_yzz_0[j] = pa_y[j] * ta_yz_yzz_0[j] - pc_y[j] * ta_yz_yzz_1[j] + 0.5 * fl1_fx * ta_z_yzz_0[j] - 0.5 * fl1_fx * ta_z_yzz_1[j] +
                                  0.5 * fl1_fx * ta_yz_zz_0[j] - 0.5 * fl1_fx * ta_yz_zz_1[j];

                ta_yyz_zzz_0[j] = pa_y[j] * ta_yz_zzz_0[j] - pc_y[j] * ta_yz_zzz_1[j] + 0.5 * fl1_fx * ta_z_zzz_0[j] - 0.5 * fl1_fx * ta_z_zzz_1[j];

                ta_yzz_xxx_0[j] = pa_y[j] * ta_zz_xxx_0[j] - pc_y[j] * ta_zz_xxx_1[j];

                ta_yzz_xxy_0[j] = pa_y[j] * ta_zz_xxy_0[j] - pc_y[j] * ta_zz_xxy_1[j] + 0.5 * fl1_fx * ta_zz_xx_0[j] - 0.5 * fl1_fx * ta_zz_xx_1[j];

                ta_yzz_xxz_0[j] = pa_y[j] * ta_zz_xxz_0[j] - pc_y[j] * ta_zz_xxz_1[j];

                ta_yzz_xyy_0[j] = pa_y[j] * ta_zz_xyy_0[j] - pc_y[j] * ta_zz_xyy_1[j] + fl1_fx * ta_zz_xy_0[j] - fl1_fx * ta_zz_xy_1[j];

                ta_yzz_xyz_0[j] = pa_y[j] * ta_zz_xyz_0[j] - pc_y[j] * ta_zz_xyz_1[j] + 0.5 * fl1_fx * ta_zz_xz_0[j] - 0.5 * fl1_fx * ta_zz_xz_1[j];

                ta_yzz_xzz_0[j] = pa_y[j] * ta_zz_xzz_0[j] - pc_y[j] * ta_zz_xzz_1[j];

                ta_yzz_yyy_0[j] = pa_y[j] * ta_zz_yyy_0[j] - pc_y[j] * ta_zz_yyy_1[j] + 1.5 * fl1_fx * ta_zz_yy_0[j] - 1.5 * fl1_fx * ta_zz_yy_1[j];

                ta_yzz_yyz_0[j] = pa_y[j] * ta_zz_yyz_0[j] - pc_y[j] * ta_zz_yyz_1[j] + fl1_fx * ta_zz_yz_0[j] - fl1_fx * ta_zz_yz_1[j];

                ta_yzz_yzz_0[j] = pa_y[j] * ta_zz_yzz_0[j] - pc_y[j] * ta_zz_yzz_1[j] + 0.5 * fl1_fx * ta_zz_zz_0[j] - 0.5 * fl1_fx * ta_zz_zz_1[j];

                ta_yzz_zzz_0[j] = pa_y[j] * ta_zz_zzz_0[j] - pc_y[j] * ta_zz_zzz_1[j];

                ta_zzz_xxx_0[j] = pa_z[j] * ta_zz_xxx_0[j] - pc_z[j] * ta_zz_xxx_1[j] + fl1_fx * ta_z_xxx_0[j] - fl1_fx * ta_z_xxx_1[j];

                ta_zzz_xxy_0[j] = pa_z[j] * ta_zz_xxy_0[j] - pc_z[j] * ta_zz_xxy_1[j] + fl1_fx * ta_z_xxy_0[j] - fl1_fx * ta_z_xxy_1[j];

                ta_zzz_xxz_0[j] = pa_z[j] * ta_zz_xxz_0[j] - pc_z[j] * ta_zz_xxz_1[j] + fl1_fx * ta_z_xxz_0[j] - fl1_fx * ta_z_xxz_1[j] +
                                  0.5 * fl1_fx * ta_zz_xx_0[j] - 0.5 * fl1_fx * ta_zz_xx_1[j];

                ta_zzz_xyy_0[j] = pa_z[j] * ta_zz_xyy_0[j] - pc_z[j] * ta_zz_xyy_1[j] + fl1_fx * ta_z_xyy_0[j] - fl1_fx * ta_z_xyy_1[j];

                ta_zzz_xyz_0[j] = pa_z[j] * ta_zz_xyz_0[j] - pc_z[j] * ta_zz_xyz_1[j] + fl1_fx * ta_z_xyz_0[j] - fl1_fx * ta_z_xyz_1[j] +
                                  0.5 * fl1_fx * ta_zz_xy_0[j] - 0.5 * fl1_fx * ta_zz_xy_1[j];

                ta_zzz_xzz_0[j] = pa_z[j] * ta_zz_xzz_0[j] - pc_z[j] * ta_zz_xzz_1[j] + fl1_fx * ta_z_xzz_0[j] - fl1_fx * ta_z_xzz_1[j] +
                                  fl1_fx * ta_zz_xz_0[j] - fl1_fx * ta_zz_xz_1[j];

                ta_zzz_yyy_0[j] = pa_z[j] * ta_zz_yyy_0[j] - pc_z[j] * ta_zz_yyy_1[j] + fl1_fx * ta_z_yyy_0[j] - fl1_fx * ta_z_yyy_1[j];

                ta_zzz_yyz_0[j] = pa_z[j] * ta_zz_yyz_0[j] - pc_z[j] * ta_zz_yyz_1[j] + fl1_fx * ta_z_yyz_0[j] - fl1_fx * ta_z_yyz_1[j] +
                                  0.5 * fl1_fx * ta_zz_yy_0[j] - 0.5 * fl1_fx * ta_zz_yy_1[j];

                ta_zzz_yzz_0[j] = pa_z[j] * ta_zz_yzz_0[j] - pc_z[j] * ta_zz_yzz_1[j] + fl1_fx * ta_z_yzz_0[j] - fl1_fx * ta_z_yzz_1[j] +
                                  fl1_fx * ta_zz_yz_0[j] - fl1_fx * ta_zz_yz_1[j];

                ta_zzz_zzz_0[j] = pa_z[j] * ta_zz_zzz_0[j] - pc_z[j] * ta_zz_zzz_1[j] + fl1_fx * ta_z_zzz_0[j] - fl1_fx * ta_z_zzz_1[j] +
                                  1.5 * fl1_fx * ta_zz_zz_0[j] - 1.5 * fl1_fx * ta_zz_zz_1[j];
            }

            idx++;
        }
    }
}

}  // namespace npotrecfunc
