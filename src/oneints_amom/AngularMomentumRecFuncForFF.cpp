//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "AngularMomentumRecFuncForFF.hpp"

namespace amomrecfunc {  // amomrecfunc namespace

void
compAngularMomentumForFF(CMemBlock2D<double>&       primBuffer,
                         const CRecursionMap&       recursionMap,
                         const CMemBlock2D<double>& osFactors,
                         const CMemBlock2D<double>& paDistances,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const int32_t              iContrGto)
{
    amomrecfunc::compAngularMomentumForFF_0_50(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForFF_50_100(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForFF_100_150(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForFF_150_200(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForFF_200_250(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForFF_250_300(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compAngularMomentumForFF_0_50(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_l_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_3_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tlx_y_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 15);

        auto tly_y_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 15);

        auto tlz_y_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 15);

        auto tlx_y_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 16);

        auto tly_y_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 16);

        auto tlx_xx_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx);

        auto tly_xx_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx);

        auto tlz_xx_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx);

        auto tlx_xx_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 1);

        auto tly_xx_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 1);

        auto tlz_xx_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 1);

        auto tlx_xx_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 2);

        auto tly_xx_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 2);

        auto tlz_xx_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 2);

        auto tlx_xx_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 3);

        auto tly_xx_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 3);

        auto tlz_xx_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 3);

        auto tlx_xx_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 4);

        auto tly_xx_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 4);

        auto tlz_xx_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 4);

        auto tlx_xx_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 5);

        auto tly_xx_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 5);

        auto tlz_xx_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 5);

        auto tlx_xy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 6);

        auto tly_xy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 6);

        auto tlz_xy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 6);

        auto tlx_xy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 7);

        auto tly_xy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 7);

        auto tlz_xy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 7);

        auto tlx_xy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 8);

        auto tly_xy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 8);

        auto tlz_xy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 8);

        auto tlx_xy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 9);

        auto tly_xy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 9);

        auto tlz_xy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 9);

        auto tlx_xy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 10);

        auto tly_xy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 10);

        auto tlz_xy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 10);

        auto tlx_xy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 11);

        auto tly_xy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 11);

        auto tlz_xy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 11);

        auto tpy_xx_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx);

        auto tpz_xx_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx);

        auto tpy_xx_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 1);

        auto tpz_xx_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 1);

        auto tpy_xx_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 2);

        auto tpz_xx_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 2);

        auto tpy_xx_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 3);

        auto tpz_xx_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 3);

        auto tpy_xx_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 4);

        auto tpz_xx_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 4);

        auto tpy_xx_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 5);

        auto tpz_xx_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 5);

        auto tpy_xx_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 6);

        auto tpz_xx_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 6);

        auto tpy_xx_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 7);

        auto tpz_xx_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 7);

        auto tpy_xx_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 8);

        auto tpz_xx_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 8);

        auto tpy_xx_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 9);

        auto tpz_xx_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 9);

        auto tpy_xy_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 10);

        auto tpz_xy_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 10);

        auto tpy_xy_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 11);

        auto tpz_xy_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 11);

        auto tpy_xy_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 12);

        auto tpz_xy_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 12);

        auto tpy_xy_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 13);

        auto tpz_xy_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 13);

        auto tpy_xy_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 14);

        auto tpz_xy_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 14);

        auto tpy_xy_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 15);

        auto tpz_xy_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 15);

        auto tpz_xy_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 16);

        auto tdy_xx_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx);

        auto tdz_xx_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx);

        auto tdy_xx_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 1);

        auto tdz_xx_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 1);

        auto tdy_xx_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 2);

        auto tdz_xx_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 2);

        auto tdy_xx_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 3);

        auto tdz_xx_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 3);

        auto tdy_xx_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 4);

        auto tdz_xx_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 4);

        auto tdy_xx_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 5);

        auto tdz_xx_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 5);

        auto tdy_xx_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 6);

        auto tdz_xx_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 6);

        auto tdy_xx_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 7);

        auto tdz_xx_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 7);

        auto tdy_xx_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 8);

        auto tdz_xx_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 8);

        auto tdy_xx_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 9);

        auto tdz_xx_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 9);

        auto tdy_xy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 10);

        auto tdz_xy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 10);

        auto tdy_xy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 11);

        auto tdz_xy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 11);

        auto tdy_xy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 12);

        auto tdz_xy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 12);

        auto tdy_xy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 13);

        auto tdz_xy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 13);

        auto tdy_xy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 14);

        auto tdz_xy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 14);

        auto tdy_xy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 15);

        auto tdz_xy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 15);

        auto tdz_xy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 16);

        // set up pointers to integrals

        auto tlx_xxx_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx);

        auto tly_xxx_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx);

        auto tlz_xxx_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx);

        auto tlx_xxx_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 1);

        auto tly_xxx_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 1);

        auto tlz_xxx_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 1);

        auto tlx_xxx_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 2);

        auto tly_xxx_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 2);

        auto tlz_xxx_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 2);

        auto tlx_xxx_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 3);

        auto tly_xxx_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 3);

        auto tlz_xxx_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 3);

        auto tlx_xxx_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 4);

        auto tly_xxx_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 4);

        auto tlz_xxx_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 4);

        auto tlx_xxx_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 5);

        auto tly_xxx_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 5);

        auto tlz_xxx_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 5);

        auto tlx_xxx_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 6);

        auto tly_xxx_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 6);

        auto tlz_xxx_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 6);

        auto tlx_xxx_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 7);

        auto tly_xxx_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 7);

        auto tlz_xxx_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 7);

        auto tlx_xxx_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 8);

        auto tly_xxx_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 8);

        auto tlz_xxx_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 8);

        auto tlx_xxx_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 9);

        auto tly_xxx_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 9);

        auto tlz_xxx_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 9);

        auto tlx_xxy_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 10);

        auto tly_xxy_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 10);

        auto tlz_xxy_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 10);

        auto tlx_xxy_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 11);

        auto tly_xxy_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 11);

        auto tlz_xxy_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 11);

        auto tlx_xxy_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 12);

        auto tly_xxy_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 12);

        auto tlz_xxy_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 12);

        auto tlx_xxy_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 13);

        auto tly_xxy_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 13);

        auto tlz_xxy_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 13);

        auto tlx_xxy_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 14);

        auto tly_xxy_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 14);

        auto tlz_xxy_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 14);

        auto tlx_xxy_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 15);

        auto tly_xxy_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 15);

        auto tlz_xxy_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 15);

        auto tlx_xxy_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 16);

        auto tly_xxy_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 16);

        // Batch of Integrals (0,50)

        #pragma omp simd aligned(fgb, fx, pa_x, tdy_xx_xxx_0, tdy_xx_xxy_0, tdy_xx_xxz_0, tdy_xx_xyy_0, \
                                     tdy_xx_xyz_0, tdy_xx_xzz_0, tdy_xx_yyy_0, tdy_xx_yyz_0, tdy_xx_yzz_0, tdy_xx_zzz_0, \
                                     tdy_xy_xxx_0, tdy_xy_xxy_0, tdy_xy_xxz_0, tdy_xy_xyy_0, tdy_xy_xyz_0, tdy_xy_xzz_0, \
                                     tdz_xx_xxx_0, tdz_xx_xxy_0, tdz_xx_xxz_0, tdz_xx_xyy_0, tdz_xx_xyz_0, tdz_xx_xzz_0, \
                                     tdz_xx_yyy_0, tdz_xx_yyz_0, tdz_xx_yzz_0, tdz_xx_zzz_0, tdz_xy_xxx_0, tdz_xy_xxy_0, \
                                     tdz_xy_xxz_0, tdz_xy_xyy_0, tdz_xy_xyz_0, tdz_xy_xzz_0, tdz_xy_yyy_0, tlx_x_xxx_0, \
                                     tlx_x_xxy_0, tlx_x_xxz_0, tlx_x_xyy_0, tlx_x_xyz_0, tlx_x_xzz_0, tlx_x_yyy_0, \
                                     tlx_x_yyz_0, tlx_x_yzz_0, tlx_x_zzz_0, tlx_xx_xx_0, tlx_xx_xxx_0, tlx_xx_xxy_0, \
                                     tlx_xx_xxz_0, tlx_xx_xy_0, tlx_xx_xyy_0, tlx_xx_xyz_0, tlx_xx_xz_0, tlx_xx_xzz_0, \
                                     tlx_xx_yy_0, tlx_xx_yyy_0, tlx_xx_yyz_0, tlx_xx_yz_0, tlx_xx_yzz_0, tlx_xx_zz_0, \
                                     tlx_xx_zzz_0, tlx_xxx_xxx_0, tlx_xxx_xxy_0, tlx_xxx_xxz_0, tlx_xxx_xyy_0, \
                                     tlx_xxx_xyz_0, tlx_xxx_xzz_0, tlx_xxx_yyy_0, tlx_xxx_yyz_0, tlx_xxx_yzz_0, \
                                     tlx_xxx_zzz_0, tlx_xxy_xxx_0, tlx_xxy_xxy_0, tlx_xxy_xxz_0, tlx_xxy_xyy_0, \
                                     tlx_xxy_xyz_0, tlx_xxy_xzz_0, tlx_xxy_yyy_0, tlx_xy_xx_0, tlx_xy_xxx_0, \
                                     tlx_xy_xxy_0, tlx_xy_xxz_0, tlx_xy_xy_0, tlx_xy_xyy_0, tlx_xy_xyz_0, tlx_xy_xz_0, \
                                     tlx_xy_xzz_0, tlx_xy_yy_0, tlx_xy_yyy_0, tlx_xy_yz_0, tlx_xy_zz_0, tlx_y_xxx_0, \
                                     tlx_y_xxy_0, tlx_y_xxz_0, tlx_y_xyy_0, tlx_y_xyz_0, tlx_y_xzz_0, tlx_y_yyy_0, \
                                     tly_x_xxx_0, tly_x_xxy_0, tly_x_xxz_0, tly_x_xyy_0, tly_x_xyz_0, tly_x_xzz_0, \
                                     tly_x_yyy_0, tly_x_yyz_0, tly_x_yzz_0, tly_x_zzz_0, tly_xx_xx_0, tly_xx_xxx_0, \
                                     tly_xx_xxy_0, tly_xx_xxz_0, tly_xx_xy_0, tly_xx_xyy_0, tly_xx_xyz_0, tly_xx_xz_0, \
                                     tly_xx_xzz_0, tly_xx_yy_0, tly_xx_yyy_0, tly_xx_yyz_0, tly_xx_yz_0, tly_xx_yzz_0, \
                                     tly_xx_zz_0, tly_xx_zzz_0, tly_xxx_xxx_0, tly_xxx_xxy_0, tly_xxx_xxz_0, \
                                     tly_xxx_xyy_0, tly_xxx_xyz_0, tly_xxx_xzz_0, tly_xxx_yyy_0, tly_xxx_yyz_0, \
                                     tly_xxx_yzz_0, tly_xxx_zzz_0, tly_xxy_xxx_0, tly_xxy_xxy_0, tly_xxy_xxz_0, \
                                     tly_xxy_xyy_0, tly_xxy_xyz_0, tly_xxy_xzz_0, tly_xxy_yyy_0, tly_xy_xx_0, \
                                     tly_xy_xxx_0, tly_xy_xxy_0, tly_xy_xxz_0, tly_xy_xy_0, tly_xy_xyy_0, tly_xy_xyz_0, \
                                     tly_xy_xz_0, tly_xy_xzz_0, tly_xy_yy_0, tly_xy_yyy_0, tly_xy_yz_0, tly_xy_zz_0, \
                                     tly_y_xxx_0, tly_y_xxy_0, tly_y_xxz_0, tly_y_xyy_0, tly_y_xyz_0, tly_y_xzz_0, \
                                     tly_y_yyy_0, tlz_x_xxx_0, tlz_x_xxy_0, tlz_x_xxz_0, tlz_x_xyy_0, tlz_x_xyz_0, \
                                     tlz_x_xzz_0, tlz_x_yyy_0, tlz_x_yyz_0, tlz_x_yzz_0, tlz_x_zzz_0, tlz_xx_xx_0, \
                                     tlz_xx_xxx_0, tlz_xx_xxy_0, tlz_xx_xxz_0, tlz_xx_xy_0, tlz_xx_xyy_0, tlz_xx_xyz_0, \
                                     tlz_xx_xz_0, tlz_xx_xzz_0, tlz_xx_yy_0, tlz_xx_yyy_0, tlz_xx_yyz_0, tlz_xx_yz_0, \
                                     tlz_xx_yzz_0, tlz_xx_zz_0, tlz_xx_zzz_0, tlz_xxx_xxx_0, tlz_xxx_xxy_0, \
                                     tlz_xxx_xxz_0, tlz_xxx_xyy_0, tlz_xxx_xyz_0, tlz_xxx_xzz_0, tlz_xxx_yyy_0, \
                                     tlz_xxx_yyz_0, tlz_xxx_yzz_0, tlz_xxx_zzz_0, tlz_xxy_xxx_0, tlz_xxy_xxy_0, \
                                     tlz_xxy_xxz_0, tlz_xxy_xyy_0, tlz_xxy_xyz_0, tlz_xxy_xzz_0, tlz_xy_xx_0, \
                                     tlz_xy_xxx_0, tlz_xy_xxy_0, tlz_xy_xxz_0, tlz_xy_xy_0, tlz_xy_xyy_0, tlz_xy_xyz_0, \
                                     tlz_xy_xz_0, tlz_xy_xzz_0, tlz_xy_yy_0, tlz_xy_yz_0, tlz_xy_zz_0, tlz_y_xxx_0, \
                                     tlz_y_xxy_0, tlz_y_xxz_0, tlz_y_xyy_0, tlz_y_xyz_0, tlz_y_xzz_0, tpy_xx_xxx_0, \
                                     tpy_xx_xxy_0, tpy_xx_xxz_0, tpy_xx_xyy_0, tpy_xx_xyz_0, tpy_xx_xzz_0, tpy_xx_yyy_0, \
                                     tpy_xx_yyz_0, tpy_xx_yzz_0, tpy_xx_zzz_0, tpy_xy_xxx_0, tpy_xy_xxy_0, tpy_xy_xxz_0, \
                                     tpy_xy_xyy_0, tpy_xy_xyz_0, tpy_xy_xzz_0, tpz_xx_xxx_0, tpz_xx_xxy_0, tpz_xx_xxz_0, \
                                     tpz_xx_xyy_0, tpz_xx_xyz_0, tpz_xx_xzz_0, tpz_xx_yyy_0, tpz_xx_yyz_0, tpz_xx_yzz_0, \
                                     tpz_xx_zzz_0, tpz_xy_xxx_0, tpz_xy_xxy_0, tpz_xy_xxz_0, tpz_xy_xyy_0, tpz_xy_xyz_0, \
                                     tpz_xy_xzz_0, tpz_xy_yyy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_xxx_xxx_0[j] = pa_x[j] * tlx_xx_xxx_0[j] + fl1_fx * tlx_x_xxx_0[j] + 1.5 * fl1_fx * tlx_xx_xx_0[j];

            tly_xxx_xxx_0[j] = pa_x[j] * tly_xx_xxx_0[j] + fl1_fx * tly_x_xxx_0[j] + 1.5 * fl1_fx * tly_xx_xx_0[j] + 0.5 * fl1_fx * tpz_xx_xxx_0[j] +
                               fl1_fx * fl1_fgb * tdz_xx_xxx_0[j];

            tlz_xxx_xxx_0[j] = pa_x[j] * tlz_xx_xxx_0[j] + fl1_fx * tlz_x_xxx_0[j] + 1.5 * fl1_fx * tlz_xx_xx_0[j] - 0.5 * fl1_fx * tpy_xx_xxx_0[j] -
                               fl1_fx * fl1_fgb * tdy_xx_xxx_0[j];

            tlx_xxx_xxy_0[j] = pa_x[j] * tlx_xx_xxy_0[j] + fl1_fx * tlx_x_xxy_0[j] + fl1_fx * tlx_xx_xy_0[j];

            tly_xxx_xxy_0[j] = pa_x[j] * tly_xx_xxy_0[j] + fl1_fx * tly_x_xxy_0[j] + fl1_fx * tly_xx_xy_0[j] + 0.5 * fl1_fx * tpz_xx_xxy_0[j] +
                               fl1_fx * fl1_fgb * tdz_xx_xxy_0[j];

            tlz_xxx_xxy_0[j] = pa_x[j] * tlz_xx_xxy_0[j] + fl1_fx * tlz_x_xxy_0[j] + fl1_fx * tlz_xx_xy_0[j] - 0.5 * fl1_fx * tpy_xx_xxy_0[j] -
                               fl1_fx * fl1_fgb * tdy_xx_xxy_0[j];

            tlx_xxx_xxz_0[j] = pa_x[j] * tlx_xx_xxz_0[j] + fl1_fx * tlx_x_xxz_0[j] + fl1_fx * tlx_xx_xz_0[j];

            tly_xxx_xxz_0[j] = pa_x[j] * tly_xx_xxz_0[j] + fl1_fx * tly_x_xxz_0[j] + fl1_fx * tly_xx_xz_0[j] + 0.5 * fl1_fx * tpz_xx_xxz_0[j] +
                               fl1_fx * fl1_fgb * tdz_xx_xxz_0[j];

            tlz_xxx_xxz_0[j] = pa_x[j] * tlz_xx_xxz_0[j] + fl1_fx * tlz_x_xxz_0[j] + fl1_fx * tlz_xx_xz_0[j] - 0.5 * fl1_fx * tpy_xx_xxz_0[j] -
                               fl1_fx * fl1_fgb * tdy_xx_xxz_0[j];

            tlx_xxx_xyy_0[j] = pa_x[j] * tlx_xx_xyy_0[j] + fl1_fx * tlx_x_xyy_0[j] + 0.5 * fl1_fx * tlx_xx_yy_0[j];

            tly_xxx_xyy_0[j] = pa_x[j] * tly_xx_xyy_0[j] + fl1_fx * tly_x_xyy_0[j] + 0.5 * fl1_fx * tly_xx_yy_0[j] + 0.5 * fl1_fx * tpz_xx_xyy_0[j] +
                               fl1_fx * fl1_fgb * tdz_xx_xyy_0[j];

            tlz_xxx_xyy_0[j] = pa_x[j] * tlz_xx_xyy_0[j] + fl1_fx * tlz_x_xyy_0[j] + 0.5 * fl1_fx * tlz_xx_yy_0[j] - 0.5 * fl1_fx * tpy_xx_xyy_0[j] -
                               fl1_fx * fl1_fgb * tdy_xx_xyy_0[j];

            tlx_xxx_xyz_0[j] = pa_x[j] * tlx_xx_xyz_0[j] + fl1_fx * tlx_x_xyz_0[j] + 0.5 * fl1_fx * tlx_xx_yz_0[j];

            tly_xxx_xyz_0[j] = pa_x[j] * tly_xx_xyz_0[j] + fl1_fx * tly_x_xyz_0[j] + 0.5 * fl1_fx * tly_xx_yz_0[j] + 0.5 * fl1_fx * tpz_xx_xyz_0[j] +
                               fl1_fx * fl1_fgb * tdz_xx_xyz_0[j];

            tlz_xxx_xyz_0[j] = pa_x[j] * tlz_xx_xyz_0[j] + fl1_fx * tlz_x_xyz_0[j] + 0.5 * fl1_fx * tlz_xx_yz_0[j] - 0.5 * fl1_fx * tpy_xx_xyz_0[j] -
                               fl1_fx * fl1_fgb * tdy_xx_xyz_0[j];

            tlx_xxx_xzz_0[j] = pa_x[j] * tlx_xx_xzz_0[j] + fl1_fx * tlx_x_xzz_0[j] + 0.5 * fl1_fx * tlx_xx_zz_0[j];

            tly_xxx_xzz_0[j] = pa_x[j] * tly_xx_xzz_0[j] + fl1_fx * tly_x_xzz_0[j] + 0.5 * fl1_fx * tly_xx_zz_0[j] + 0.5 * fl1_fx * tpz_xx_xzz_0[j] +
                               fl1_fx * fl1_fgb * tdz_xx_xzz_0[j];

            tlz_xxx_xzz_0[j] = pa_x[j] * tlz_xx_xzz_0[j] + fl1_fx * tlz_x_xzz_0[j] + 0.5 * fl1_fx * tlz_xx_zz_0[j] - 0.5 * fl1_fx * tpy_xx_xzz_0[j] -
                               fl1_fx * fl1_fgb * tdy_xx_xzz_0[j];

            tlx_xxx_yyy_0[j] = pa_x[j] * tlx_xx_yyy_0[j] + fl1_fx * tlx_x_yyy_0[j];

            tly_xxx_yyy_0[j] =
                pa_x[j] * tly_xx_yyy_0[j] + fl1_fx * tly_x_yyy_0[j] + 0.5 * fl1_fx * tpz_xx_yyy_0[j] + fl1_fx * fl1_fgb * tdz_xx_yyy_0[j];

            tlz_xxx_yyy_0[j] =
                pa_x[j] * tlz_xx_yyy_0[j] + fl1_fx * tlz_x_yyy_0[j] - 0.5 * fl1_fx * tpy_xx_yyy_0[j] - fl1_fx * fl1_fgb * tdy_xx_yyy_0[j];

            tlx_xxx_yyz_0[j] = pa_x[j] * tlx_xx_yyz_0[j] + fl1_fx * tlx_x_yyz_0[j];

            tly_xxx_yyz_0[j] =
                pa_x[j] * tly_xx_yyz_0[j] + fl1_fx * tly_x_yyz_0[j] + 0.5 * fl1_fx * tpz_xx_yyz_0[j] + fl1_fx * fl1_fgb * tdz_xx_yyz_0[j];

            tlz_xxx_yyz_0[j] =
                pa_x[j] * tlz_xx_yyz_0[j] + fl1_fx * tlz_x_yyz_0[j] - 0.5 * fl1_fx * tpy_xx_yyz_0[j] - fl1_fx * fl1_fgb * tdy_xx_yyz_0[j];

            tlx_xxx_yzz_0[j] = pa_x[j] * tlx_xx_yzz_0[j] + fl1_fx * tlx_x_yzz_0[j];

            tly_xxx_yzz_0[j] =
                pa_x[j] * tly_xx_yzz_0[j] + fl1_fx * tly_x_yzz_0[j] + 0.5 * fl1_fx * tpz_xx_yzz_0[j] + fl1_fx * fl1_fgb * tdz_xx_yzz_0[j];

            tlz_xxx_yzz_0[j] =
                pa_x[j] * tlz_xx_yzz_0[j] + fl1_fx * tlz_x_yzz_0[j] - 0.5 * fl1_fx * tpy_xx_yzz_0[j] - fl1_fx * fl1_fgb * tdy_xx_yzz_0[j];

            tlx_xxx_zzz_0[j] = pa_x[j] * tlx_xx_zzz_0[j] + fl1_fx * tlx_x_zzz_0[j];

            tly_xxx_zzz_0[j] =
                pa_x[j] * tly_xx_zzz_0[j] + fl1_fx * tly_x_zzz_0[j] + 0.5 * fl1_fx * tpz_xx_zzz_0[j] + fl1_fx * fl1_fgb * tdz_xx_zzz_0[j];

            tlz_xxx_zzz_0[j] =
                pa_x[j] * tlz_xx_zzz_0[j] + fl1_fx * tlz_x_zzz_0[j] - 0.5 * fl1_fx * tpy_xx_zzz_0[j] - fl1_fx * fl1_fgb * tdy_xx_zzz_0[j];

            tlx_xxy_xxx_0[j] = pa_x[j] * tlx_xy_xxx_0[j] + 0.5 * fl1_fx * tlx_y_xxx_0[j] + 1.5 * fl1_fx * tlx_xy_xx_0[j];

            tly_xxy_xxx_0[j] = pa_x[j] * tly_xy_xxx_0[j] + 0.5 * fl1_fx * tly_y_xxx_0[j] + 1.5 * fl1_fx * tly_xy_xx_0[j] +
                               0.5 * fl1_fx * tpz_xy_xxx_0[j] + fl1_fx * fl1_fgb * tdz_xy_xxx_0[j];

            tlz_xxy_xxx_0[j] = pa_x[j] * tlz_xy_xxx_0[j] + 0.5 * fl1_fx * tlz_y_xxx_0[j] + 1.5 * fl1_fx * tlz_xy_xx_0[j] -
                               0.5 * fl1_fx * tpy_xy_xxx_0[j] - fl1_fx * fl1_fgb * tdy_xy_xxx_0[j];

            tlx_xxy_xxy_0[j] = pa_x[j] * tlx_xy_xxy_0[j] + 0.5 * fl1_fx * tlx_y_xxy_0[j] + fl1_fx * tlx_xy_xy_0[j];

            tly_xxy_xxy_0[j] = pa_x[j] * tly_xy_xxy_0[j] + 0.5 * fl1_fx * tly_y_xxy_0[j] + fl1_fx * tly_xy_xy_0[j] + 0.5 * fl1_fx * tpz_xy_xxy_0[j] +
                               fl1_fx * fl1_fgb * tdz_xy_xxy_0[j];

            tlz_xxy_xxy_0[j] = pa_x[j] * tlz_xy_xxy_0[j] + 0.5 * fl1_fx * tlz_y_xxy_0[j] + fl1_fx * tlz_xy_xy_0[j] - 0.5 * fl1_fx * tpy_xy_xxy_0[j] -
                               fl1_fx * fl1_fgb * tdy_xy_xxy_0[j];

            tlx_xxy_xxz_0[j] = pa_x[j] * tlx_xy_xxz_0[j] + 0.5 * fl1_fx * tlx_y_xxz_0[j] + fl1_fx * tlx_xy_xz_0[j];

            tly_xxy_xxz_0[j] = pa_x[j] * tly_xy_xxz_0[j] + 0.5 * fl1_fx * tly_y_xxz_0[j] + fl1_fx * tly_xy_xz_0[j] + 0.5 * fl1_fx * tpz_xy_xxz_0[j] +
                               fl1_fx * fl1_fgb * tdz_xy_xxz_0[j];

            tlz_xxy_xxz_0[j] = pa_x[j] * tlz_xy_xxz_0[j] + 0.5 * fl1_fx * tlz_y_xxz_0[j] + fl1_fx * tlz_xy_xz_0[j] - 0.5 * fl1_fx * tpy_xy_xxz_0[j] -
                               fl1_fx * fl1_fgb * tdy_xy_xxz_0[j];

            tlx_xxy_xyy_0[j] = pa_x[j] * tlx_xy_xyy_0[j] + 0.5 * fl1_fx * tlx_y_xyy_0[j] + 0.5 * fl1_fx * tlx_xy_yy_0[j];

            tly_xxy_xyy_0[j] = pa_x[j] * tly_xy_xyy_0[j] + 0.5 * fl1_fx * tly_y_xyy_0[j] + 0.5 * fl1_fx * tly_xy_yy_0[j] +
                               0.5 * fl1_fx * tpz_xy_xyy_0[j] + fl1_fx * fl1_fgb * tdz_xy_xyy_0[j];

            tlz_xxy_xyy_0[j] = pa_x[j] * tlz_xy_xyy_0[j] + 0.5 * fl1_fx * tlz_y_xyy_0[j] + 0.5 * fl1_fx * tlz_xy_yy_0[j] -
                               0.5 * fl1_fx * tpy_xy_xyy_0[j] - fl1_fx * fl1_fgb * tdy_xy_xyy_0[j];

            tlx_xxy_xyz_0[j] = pa_x[j] * tlx_xy_xyz_0[j] + 0.5 * fl1_fx * tlx_y_xyz_0[j] + 0.5 * fl1_fx * tlx_xy_yz_0[j];

            tly_xxy_xyz_0[j] = pa_x[j] * tly_xy_xyz_0[j] + 0.5 * fl1_fx * tly_y_xyz_0[j] + 0.5 * fl1_fx * tly_xy_yz_0[j] +
                               0.5 * fl1_fx * tpz_xy_xyz_0[j] + fl1_fx * fl1_fgb * tdz_xy_xyz_0[j];

            tlz_xxy_xyz_0[j] = pa_x[j] * tlz_xy_xyz_0[j] + 0.5 * fl1_fx * tlz_y_xyz_0[j] + 0.5 * fl1_fx * tlz_xy_yz_0[j] -
                               0.5 * fl1_fx * tpy_xy_xyz_0[j] - fl1_fx * fl1_fgb * tdy_xy_xyz_0[j];

            tlx_xxy_xzz_0[j] = pa_x[j] * tlx_xy_xzz_0[j] + 0.5 * fl1_fx * tlx_y_xzz_0[j] + 0.5 * fl1_fx * tlx_xy_zz_0[j];

            tly_xxy_xzz_0[j] = pa_x[j] * tly_xy_xzz_0[j] + 0.5 * fl1_fx * tly_y_xzz_0[j] + 0.5 * fl1_fx * tly_xy_zz_0[j] +
                               0.5 * fl1_fx * tpz_xy_xzz_0[j] + fl1_fx * fl1_fgb * tdz_xy_xzz_0[j];

            tlz_xxy_xzz_0[j] = pa_x[j] * tlz_xy_xzz_0[j] + 0.5 * fl1_fx * tlz_y_xzz_0[j] + 0.5 * fl1_fx * tlz_xy_zz_0[j] -
                               0.5 * fl1_fx * tpy_xy_xzz_0[j] - fl1_fx * fl1_fgb * tdy_xy_xzz_0[j];

            tlx_xxy_yyy_0[j] = pa_x[j] * tlx_xy_yyy_0[j] + 0.5 * fl1_fx * tlx_y_yyy_0[j];

            tly_xxy_yyy_0[j] =
                pa_x[j] * tly_xy_yyy_0[j] + 0.5 * fl1_fx * tly_y_yyy_0[j] + 0.5 * fl1_fx * tpz_xy_yyy_0[j] + fl1_fx * fl1_fgb * tdz_xy_yyy_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForFF_50_100(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_l_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_3_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

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

        auto tlx_xz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 12);

        auto tly_xz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 12);

        auto tlz_xz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 12);

        auto tlx_xz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 13);

        auto tly_xz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 13);

        auto tlz_xz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 13);

        auto tlx_xz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 14);

        auto tly_xz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 14);

        auto tlz_xz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 14);

        auto tlx_xz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 15);

        auto tly_xz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 15);

        auto tlz_xz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 15);

        auto tlx_xz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 16);

        auto tly_xz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 16);

        auto tlz_xz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 16);

        auto tlx_xz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 17);

        auto tly_xz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 17);

        auto tlz_xz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 17);

        auto tlx_yy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 18);

        auto tly_yy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 18);

        auto tlz_yy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 18);

        auto tlx_yy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 19);

        auto tly_yy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 19);

        auto tlz_yy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 19);

        auto tlx_yy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 20);

        auto tly_yy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 20);

        auto tlz_yy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 20);

        auto tlx_yy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 21);

        auto tpy_xy_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 16);

        auto tpy_xy_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 17);

        auto tpz_xy_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 17);

        auto tpy_xy_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 18);

        auto tpz_xy_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 18);

        auto tpy_xy_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 19);

        auto tpz_xy_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 19);

        auto tpy_xz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 20);

        auto tpz_xz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 20);

        auto tpy_xz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 21);

        auto tpz_xz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 21);

        auto tpy_xz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 22);

        auto tpz_xz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 22);

        auto tpy_xz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 23);

        auto tpz_xz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 23);

        auto tpy_xz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 24);

        auto tpz_xz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 24);

        auto tpy_xz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 25);

        auto tpz_xz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 25);

        auto tpy_xz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 26);

        auto tpz_xz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 26);

        auto tpy_xz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 27);

        auto tpz_xz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 27);

        auto tpy_xz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 28);

        auto tpz_xz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 28);

        auto tpy_xz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 29);

        auto tpz_xz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 29);

        auto tpy_yy_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 30);

        auto tpz_yy_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 30);

        auto tpy_yy_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 31);

        auto tpz_yy_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 31);

        auto tpy_yy_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 32);

        auto tpz_yy_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 32);

        auto tdy_xy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 16);

        auto tdy_xy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 17);

        auto tdz_xy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 17);

        auto tdy_xy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 18);

        auto tdz_xy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 18);

        auto tdy_xy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 19);

        auto tdz_xy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 19);

        auto tdy_xz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 20);

        auto tdz_xz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 20);

        auto tdy_xz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 21);

        auto tdz_xz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 21);

        auto tdy_xz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 22);

        auto tdz_xz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 22);

        auto tdy_xz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 23);

        auto tdz_xz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 23);

        auto tdy_xz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 24);

        auto tdz_xz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 24);

        auto tdy_xz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 25);

        auto tdz_xz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 25);

        auto tdy_xz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 26);

        auto tdz_xz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 26);

        auto tdy_xz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 27);

        auto tdz_xz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 27);

        auto tdy_xz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 28);

        auto tdz_xz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 28);

        auto tdy_xz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 29);

        auto tdz_xz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 29);

        auto tdy_yy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 30);

        auto tdz_yy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 30);

        auto tdy_yy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 31);

        auto tdz_yy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 31);

        auto tdy_yy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 32);

        auto tdz_yy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 32);

        // set up pointers to integrals

        auto tlz_xxy_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 16);

        auto tlx_xxy_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 17);

        auto tly_xxy_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 17);

        auto tlz_xxy_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 17);

        auto tlx_xxy_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 18);

        auto tly_xxy_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 18);

        auto tlz_xxy_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 18);

        auto tlx_xxy_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 19);

        auto tly_xxy_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 19);

        auto tlz_xxy_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 19);

        auto tlx_xxz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 20);

        auto tly_xxz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 20);

        auto tlz_xxz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 20);

        auto tlx_xxz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 21);

        auto tly_xxz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 21);

        auto tlz_xxz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 21);

        auto tlx_xxz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 22);

        auto tly_xxz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 22);

        auto tlz_xxz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 22);

        auto tlx_xxz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 23);

        auto tly_xxz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 23);

        auto tlz_xxz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 23);

        auto tlx_xxz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 24);

        auto tly_xxz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 24);

        auto tlz_xxz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 24);

        auto tlx_xxz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 25);

        auto tly_xxz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 25);

        auto tlz_xxz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 25);

        auto tlx_xxz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 26);

        auto tly_xxz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 26);

        auto tlz_xxz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 26);

        auto tlx_xxz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 27);

        auto tly_xxz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 27);

        auto tlz_xxz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 27);

        auto tlx_xxz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 28);

        auto tly_xxz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 28);

        auto tlz_xxz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 28);

        auto tlx_xxz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 29);

        auto tly_xxz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 29);

        auto tlz_xxz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 29);

        auto tlx_xyy_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 30);

        auto tly_xyy_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 30);

        auto tlz_xyy_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 30);

        auto tlx_xyy_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 31);

        auto tly_xyy_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 31);

        auto tlz_xyy_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 31);

        auto tlx_xyy_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 32);

        auto tly_xyy_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 32);

        auto tlz_xyy_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 32);

        auto tlx_xyy_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 33);

        // Batch of Integrals (50,100)

        #pragma omp simd aligned(fgb, fx, pa_x, tdy_xy_yyy_0, tdy_xy_yyz_0, tdy_xy_yzz_0, tdy_xy_zzz_0, \
                                     tdy_xz_xxx_0, tdy_xz_xxy_0, tdy_xz_xxz_0, tdy_xz_xyy_0, tdy_xz_xyz_0, tdy_xz_xzz_0, \
                                     tdy_xz_yyy_0, tdy_xz_yyz_0, tdy_xz_yzz_0, tdy_xz_zzz_0, tdy_yy_xxx_0, tdy_yy_xxy_0, \
                                     tdy_yy_xxz_0, tdz_xy_yyz_0, tdz_xy_yzz_0, tdz_xy_zzz_0, tdz_xz_xxx_0, tdz_xz_xxy_0, \
                                     tdz_xz_xxz_0, tdz_xz_xyy_0, tdz_xz_xyz_0, tdz_xz_xzz_0, tdz_xz_yyy_0, tdz_xz_yyz_0, \
                                     tdz_xz_yzz_0, tdz_xz_zzz_0, tdz_yy_xxx_0, tdz_yy_xxy_0, tdz_yy_xxz_0, \
                                     tlx_xxy_yyz_0, tlx_xxy_yzz_0, tlx_xxy_zzz_0, tlx_xxz_xxx_0, tlx_xxz_xxy_0, \
                                     tlx_xxz_xxz_0, tlx_xxz_xyy_0, tlx_xxz_xyz_0, tlx_xxz_xzz_0, tlx_xxz_yyy_0, \
                                     tlx_xxz_yyz_0, tlx_xxz_yzz_0, tlx_xxz_zzz_0, tlx_xy_yyz_0, tlx_xy_yzz_0, \
                                     tlx_xy_zzz_0, tlx_xyy_xxx_0, tlx_xyy_xxy_0, tlx_xyy_xxz_0, tlx_xyy_xyy_0, \
                                     tlx_xz_xx_0, tlx_xz_xxx_0, tlx_xz_xxy_0, tlx_xz_xxz_0, tlx_xz_xy_0, tlx_xz_xyy_0, \
                                     tlx_xz_xyz_0, tlx_xz_xz_0, tlx_xz_xzz_0, tlx_xz_yy_0, tlx_xz_yyy_0, tlx_xz_yyz_0, \
                                     tlx_xz_yz_0, tlx_xz_yzz_0, tlx_xz_zz_0, tlx_xz_zzz_0, tlx_y_yyz_0, tlx_y_yzz_0, \
                                     tlx_y_zzz_0, tlx_yy_xx_0, tlx_yy_xxx_0, tlx_yy_xxy_0, tlx_yy_xxz_0, tlx_yy_xy_0, \
                                     tlx_yy_xyy_0, tlx_yy_xz_0, tlx_yy_yy_0, tlx_z_xxx_0, tlx_z_xxy_0, tlx_z_xxz_0, \
                                     tlx_z_xyy_0, tlx_z_xyz_0, tlx_z_xzz_0, tlx_z_yyy_0, tlx_z_yyz_0, tlx_z_yzz_0, \
                                     tlx_z_zzz_0, tly_xxy_yyz_0, tly_xxy_yzz_0, tly_xxy_zzz_0, tly_xxz_xxx_0, \
                                     tly_xxz_xxy_0, tly_xxz_xxz_0, tly_xxz_xyy_0, tly_xxz_xyz_0, tly_xxz_xzz_0, \
                                     tly_xxz_yyy_0, tly_xxz_yyz_0, tly_xxz_yzz_0, tly_xxz_zzz_0, tly_xy_yyz_0, \
                                     tly_xy_yzz_0, tly_xy_zzz_0, tly_xyy_xxx_0, tly_xyy_xxy_0, tly_xyy_xxz_0, \
                                     tly_xz_xx_0, tly_xz_xxx_0, tly_xz_xxy_0, tly_xz_xxz_0, tly_xz_xy_0, tly_xz_xyy_0, \
                                     tly_xz_xyz_0, tly_xz_xz_0, tly_xz_xzz_0, tly_xz_yy_0, tly_xz_yyy_0, tly_xz_yyz_0, \
                                     tly_xz_yz_0, tly_xz_yzz_0, tly_xz_zz_0, tly_xz_zzz_0, tly_y_yyz_0, tly_y_yzz_0, \
                                     tly_y_zzz_0, tly_yy_xx_0, tly_yy_xxx_0, tly_yy_xxy_0, tly_yy_xxz_0, tly_yy_xy_0, \
                                     tly_yy_xz_0, tly_z_xxx_0, tly_z_xxy_0, tly_z_xxz_0, tly_z_xyy_0, tly_z_xyz_0, \
                                     tly_z_xzz_0, tly_z_yyy_0, tly_z_yyz_0, tly_z_yzz_0, tly_z_zzz_0, tlz_xxy_yyy_0, \
                                     tlz_xxy_yyz_0, tlz_xxy_yzz_0, tlz_xxy_zzz_0, tlz_xxz_xxx_0, tlz_xxz_xxy_0, \
                                     tlz_xxz_xxz_0, tlz_xxz_xyy_0, tlz_xxz_xyz_0, tlz_xxz_xzz_0, tlz_xxz_yyy_0, \
                                     tlz_xxz_yyz_0, tlz_xxz_yzz_0, tlz_xxz_zzz_0, tlz_xy_yyy_0, tlz_xy_yyz_0, \
                                     tlz_xy_yzz_0, tlz_xy_zzz_0, tlz_xyy_xxx_0, tlz_xyy_xxy_0, tlz_xyy_xxz_0, \
                                     tlz_xz_xx_0, tlz_xz_xxx_0, tlz_xz_xxy_0, tlz_xz_xxz_0, tlz_xz_xy_0, tlz_xz_xyy_0, \
                                     tlz_xz_xyz_0, tlz_xz_xz_0, tlz_xz_xzz_0, tlz_xz_yy_0, tlz_xz_yyy_0, tlz_xz_yyz_0, \
                                     tlz_xz_yz_0, tlz_xz_yzz_0, tlz_xz_zz_0, tlz_xz_zzz_0, tlz_y_yyy_0, tlz_y_yyz_0, \
                                     tlz_y_yzz_0, tlz_y_zzz_0, tlz_yy_xx_0, tlz_yy_xxx_0, tlz_yy_xxy_0, tlz_yy_xxz_0, \
                                     tlz_yy_xy_0, tlz_yy_xz_0, tlz_z_xxx_0, tlz_z_xxy_0, tlz_z_xxz_0, tlz_z_xyy_0, \
                                     tlz_z_xyz_0, tlz_z_xzz_0, tlz_z_yyy_0, tlz_z_yyz_0, tlz_z_yzz_0, tlz_z_zzz_0, \
                                     tpy_xy_yyy_0, tpy_xy_yyz_0, tpy_xy_yzz_0, tpy_xy_zzz_0, tpy_xz_xxx_0, tpy_xz_xxy_0, \
                                     tpy_xz_xxz_0, tpy_xz_xyy_0, tpy_xz_xyz_0, tpy_xz_xzz_0, tpy_xz_yyy_0, tpy_xz_yyz_0, \
                                     tpy_xz_yzz_0, tpy_xz_zzz_0, tpy_yy_xxx_0, tpy_yy_xxy_0, tpy_yy_xxz_0, tpz_xy_yyz_0, \
                                     tpz_xy_yzz_0, tpz_xy_zzz_0, tpz_xz_xxx_0, tpz_xz_xxy_0, tpz_xz_xxz_0, tpz_xz_xyy_0, \
                                     tpz_xz_xyz_0, tpz_xz_xzz_0, tpz_xz_yyy_0, tpz_xz_yyz_0, tpz_xz_yzz_0, tpz_xz_zzz_0, \
                                     tpz_yy_xxx_0, tpz_yy_xxy_0, tpz_yy_xxz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlz_xxy_yyy_0[j] =
                pa_x[j] * tlz_xy_yyy_0[j] + 0.5 * fl1_fx * tlz_y_yyy_0[j] - 0.5 * fl1_fx * tpy_xy_yyy_0[j] - fl1_fx * fl1_fgb * tdy_xy_yyy_0[j];

            tlx_xxy_yyz_0[j] = pa_x[j] * tlx_xy_yyz_0[j] + 0.5 * fl1_fx * tlx_y_yyz_0[j];

            tly_xxy_yyz_0[j] =
                pa_x[j] * tly_xy_yyz_0[j] + 0.5 * fl1_fx * tly_y_yyz_0[j] + 0.5 * fl1_fx * tpz_xy_yyz_0[j] + fl1_fx * fl1_fgb * tdz_xy_yyz_0[j];

            tlz_xxy_yyz_0[j] =
                pa_x[j] * tlz_xy_yyz_0[j] + 0.5 * fl1_fx * tlz_y_yyz_0[j] - 0.5 * fl1_fx * tpy_xy_yyz_0[j] - fl1_fx * fl1_fgb * tdy_xy_yyz_0[j];

            tlx_xxy_yzz_0[j] = pa_x[j] * tlx_xy_yzz_0[j] + 0.5 * fl1_fx * tlx_y_yzz_0[j];

            tly_xxy_yzz_0[j] =
                pa_x[j] * tly_xy_yzz_0[j] + 0.5 * fl1_fx * tly_y_yzz_0[j] + 0.5 * fl1_fx * tpz_xy_yzz_0[j] + fl1_fx * fl1_fgb * tdz_xy_yzz_0[j];

            tlz_xxy_yzz_0[j] =
                pa_x[j] * tlz_xy_yzz_0[j] + 0.5 * fl1_fx * tlz_y_yzz_0[j] - 0.5 * fl1_fx * tpy_xy_yzz_0[j] - fl1_fx * fl1_fgb * tdy_xy_yzz_0[j];

            tlx_xxy_zzz_0[j] = pa_x[j] * tlx_xy_zzz_0[j] + 0.5 * fl1_fx * tlx_y_zzz_0[j];

            tly_xxy_zzz_0[j] =
                pa_x[j] * tly_xy_zzz_0[j] + 0.5 * fl1_fx * tly_y_zzz_0[j] + 0.5 * fl1_fx * tpz_xy_zzz_0[j] + fl1_fx * fl1_fgb * tdz_xy_zzz_0[j];

            tlz_xxy_zzz_0[j] =
                pa_x[j] * tlz_xy_zzz_0[j] + 0.5 * fl1_fx * tlz_y_zzz_0[j] - 0.5 * fl1_fx * tpy_xy_zzz_0[j] - fl1_fx * fl1_fgb * tdy_xy_zzz_0[j];

            tlx_xxz_xxx_0[j] = pa_x[j] * tlx_xz_xxx_0[j] + 0.5 * fl1_fx * tlx_z_xxx_0[j] + 1.5 * fl1_fx * tlx_xz_xx_0[j];

            tly_xxz_xxx_0[j] = pa_x[j] * tly_xz_xxx_0[j] + 0.5 * fl1_fx * tly_z_xxx_0[j] + 1.5 * fl1_fx * tly_xz_xx_0[j] +
                               0.5 * fl1_fx * tpz_xz_xxx_0[j] + fl1_fx * fl1_fgb * tdz_xz_xxx_0[j];

            tlz_xxz_xxx_0[j] = pa_x[j] * tlz_xz_xxx_0[j] + 0.5 * fl1_fx * tlz_z_xxx_0[j] + 1.5 * fl1_fx * tlz_xz_xx_0[j] -
                               0.5 * fl1_fx * tpy_xz_xxx_0[j] - fl1_fx * fl1_fgb * tdy_xz_xxx_0[j];

            tlx_xxz_xxy_0[j] = pa_x[j] * tlx_xz_xxy_0[j] + 0.5 * fl1_fx * tlx_z_xxy_0[j] + fl1_fx * tlx_xz_xy_0[j];

            tly_xxz_xxy_0[j] = pa_x[j] * tly_xz_xxy_0[j] + 0.5 * fl1_fx * tly_z_xxy_0[j] + fl1_fx * tly_xz_xy_0[j] + 0.5 * fl1_fx * tpz_xz_xxy_0[j] +
                               fl1_fx * fl1_fgb * tdz_xz_xxy_0[j];

            tlz_xxz_xxy_0[j] = pa_x[j] * tlz_xz_xxy_0[j] + 0.5 * fl1_fx * tlz_z_xxy_0[j] + fl1_fx * tlz_xz_xy_0[j] - 0.5 * fl1_fx * tpy_xz_xxy_0[j] -
                               fl1_fx * fl1_fgb * tdy_xz_xxy_0[j];

            tlx_xxz_xxz_0[j] = pa_x[j] * tlx_xz_xxz_0[j] + 0.5 * fl1_fx * tlx_z_xxz_0[j] + fl1_fx * tlx_xz_xz_0[j];

            tly_xxz_xxz_0[j] = pa_x[j] * tly_xz_xxz_0[j] + 0.5 * fl1_fx * tly_z_xxz_0[j] + fl1_fx * tly_xz_xz_0[j] + 0.5 * fl1_fx * tpz_xz_xxz_0[j] +
                               fl1_fx * fl1_fgb * tdz_xz_xxz_0[j];

            tlz_xxz_xxz_0[j] = pa_x[j] * tlz_xz_xxz_0[j] + 0.5 * fl1_fx * tlz_z_xxz_0[j] + fl1_fx * tlz_xz_xz_0[j] - 0.5 * fl1_fx * tpy_xz_xxz_0[j] -
                               fl1_fx * fl1_fgb * tdy_xz_xxz_0[j];

            tlx_xxz_xyy_0[j] = pa_x[j] * tlx_xz_xyy_0[j] + 0.5 * fl1_fx * tlx_z_xyy_0[j] + 0.5 * fl1_fx * tlx_xz_yy_0[j];

            tly_xxz_xyy_0[j] = pa_x[j] * tly_xz_xyy_0[j] + 0.5 * fl1_fx * tly_z_xyy_0[j] + 0.5 * fl1_fx * tly_xz_yy_0[j] +
                               0.5 * fl1_fx * tpz_xz_xyy_0[j] + fl1_fx * fl1_fgb * tdz_xz_xyy_0[j];

            tlz_xxz_xyy_0[j] = pa_x[j] * tlz_xz_xyy_0[j] + 0.5 * fl1_fx * tlz_z_xyy_0[j] + 0.5 * fl1_fx * tlz_xz_yy_0[j] -
                               0.5 * fl1_fx * tpy_xz_xyy_0[j] - fl1_fx * fl1_fgb * tdy_xz_xyy_0[j];

            tlx_xxz_xyz_0[j] = pa_x[j] * tlx_xz_xyz_0[j] + 0.5 * fl1_fx * tlx_z_xyz_0[j] + 0.5 * fl1_fx * tlx_xz_yz_0[j];

            tly_xxz_xyz_0[j] = pa_x[j] * tly_xz_xyz_0[j] + 0.5 * fl1_fx * tly_z_xyz_0[j] + 0.5 * fl1_fx * tly_xz_yz_0[j] +
                               0.5 * fl1_fx * tpz_xz_xyz_0[j] + fl1_fx * fl1_fgb * tdz_xz_xyz_0[j];

            tlz_xxz_xyz_0[j] = pa_x[j] * tlz_xz_xyz_0[j] + 0.5 * fl1_fx * tlz_z_xyz_0[j] + 0.5 * fl1_fx * tlz_xz_yz_0[j] -
                               0.5 * fl1_fx * tpy_xz_xyz_0[j] - fl1_fx * fl1_fgb * tdy_xz_xyz_0[j];

            tlx_xxz_xzz_0[j] = pa_x[j] * tlx_xz_xzz_0[j] + 0.5 * fl1_fx * tlx_z_xzz_0[j] + 0.5 * fl1_fx * tlx_xz_zz_0[j];

            tly_xxz_xzz_0[j] = pa_x[j] * tly_xz_xzz_0[j] + 0.5 * fl1_fx * tly_z_xzz_0[j] + 0.5 * fl1_fx * tly_xz_zz_0[j] +
                               0.5 * fl1_fx * tpz_xz_xzz_0[j] + fl1_fx * fl1_fgb * tdz_xz_xzz_0[j];

            tlz_xxz_xzz_0[j] = pa_x[j] * tlz_xz_xzz_0[j] + 0.5 * fl1_fx * tlz_z_xzz_0[j] + 0.5 * fl1_fx * tlz_xz_zz_0[j] -
                               0.5 * fl1_fx * tpy_xz_xzz_0[j] - fl1_fx * fl1_fgb * tdy_xz_xzz_0[j];

            tlx_xxz_yyy_0[j] = pa_x[j] * tlx_xz_yyy_0[j] + 0.5 * fl1_fx * tlx_z_yyy_0[j];

            tly_xxz_yyy_0[j] =
                pa_x[j] * tly_xz_yyy_0[j] + 0.5 * fl1_fx * tly_z_yyy_0[j] + 0.5 * fl1_fx * tpz_xz_yyy_0[j] + fl1_fx * fl1_fgb * tdz_xz_yyy_0[j];

            tlz_xxz_yyy_0[j] =
                pa_x[j] * tlz_xz_yyy_0[j] + 0.5 * fl1_fx * tlz_z_yyy_0[j] - 0.5 * fl1_fx * tpy_xz_yyy_0[j] - fl1_fx * fl1_fgb * tdy_xz_yyy_0[j];

            tlx_xxz_yyz_0[j] = pa_x[j] * tlx_xz_yyz_0[j] + 0.5 * fl1_fx * tlx_z_yyz_0[j];

            tly_xxz_yyz_0[j] =
                pa_x[j] * tly_xz_yyz_0[j] + 0.5 * fl1_fx * tly_z_yyz_0[j] + 0.5 * fl1_fx * tpz_xz_yyz_0[j] + fl1_fx * fl1_fgb * tdz_xz_yyz_0[j];

            tlz_xxz_yyz_0[j] =
                pa_x[j] * tlz_xz_yyz_0[j] + 0.5 * fl1_fx * tlz_z_yyz_0[j] - 0.5 * fl1_fx * tpy_xz_yyz_0[j] - fl1_fx * fl1_fgb * tdy_xz_yyz_0[j];

            tlx_xxz_yzz_0[j] = pa_x[j] * tlx_xz_yzz_0[j] + 0.5 * fl1_fx * tlx_z_yzz_0[j];

            tly_xxz_yzz_0[j] =
                pa_x[j] * tly_xz_yzz_0[j] + 0.5 * fl1_fx * tly_z_yzz_0[j] + 0.5 * fl1_fx * tpz_xz_yzz_0[j] + fl1_fx * fl1_fgb * tdz_xz_yzz_0[j];

            tlz_xxz_yzz_0[j] =
                pa_x[j] * tlz_xz_yzz_0[j] + 0.5 * fl1_fx * tlz_z_yzz_0[j] - 0.5 * fl1_fx * tpy_xz_yzz_0[j] - fl1_fx * fl1_fgb * tdy_xz_yzz_0[j];

            tlx_xxz_zzz_0[j] = pa_x[j] * tlx_xz_zzz_0[j] + 0.5 * fl1_fx * tlx_z_zzz_0[j];

            tly_xxz_zzz_0[j] =
                pa_x[j] * tly_xz_zzz_0[j] + 0.5 * fl1_fx * tly_z_zzz_0[j] + 0.5 * fl1_fx * tpz_xz_zzz_0[j] + fl1_fx * fl1_fgb * tdz_xz_zzz_0[j];

            tlz_xxz_zzz_0[j] =
                pa_x[j] * tlz_xz_zzz_0[j] + 0.5 * fl1_fx * tlz_z_zzz_0[j] - 0.5 * fl1_fx * tpy_xz_zzz_0[j] - fl1_fx * fl1_fgb * tdy_xz_zzz_0[j];

            tlx_xyy_xxx_0[j] = pa_x[j] * tlx_yy_xxx_0[j] + 1.5 * fl1_fx * tlx_yy_xx_0[j];

            tly_xyy_xxx_0[j] =
                pa_x[j] * tly_yy_xxx_0[j] + 1.5 * fl1_fx * tly_yy_xx_0[j] + 0.5 * fl1_fx * tpz_yy_xxx_0[j] + fl1_fx * fl1_fgb * tdz_yy_xxx_0[j];

            tlz_xyy_xxx_0[j] =
                pa_x[j] * tlz_yy_xxx_0[j] + 1.5 * fl1_fx * tlz_yy_xx_0[j] - 0.5 * fl1_fx * tpy_yy_xxx_0[j] - fl1_fx * fl1_fgb * tdy_yy_xxx_0[j];

            tlx_xyy_xxy_0[j] = pa_x[j] * tlx_yy_xxy_0[j] + fl1_fx * tlx_yy_xy_0[j];

            tly_xyy_xxy_0[j] =
                pa_x[j] * tly_yy_xxy_0[j] + fl1_fx * tly_yy_xy_0[j] + 0.5 * fl1_fx * tpz_yy_xxy_0[j] + fl1_fx * fl1_fgb * tdz_yy_xxy_0[j];

            tlz_xyy_xxy_0[j] =
                pa_x[j] * tlz_yy_xxy_0[j] + fl1_fx * tlz_yy_xy_0[j] - 0.5 * fl1_fx * tpy_yy_xxy_0[j] - fl1_fx * fl1_fgb * tdy_yy_xxy_0[j];

            tlx_xyy_xxz_0[j] = pa_x[j] * tlx_yy_xxz_0[j] + fl1_fx * tlx_yy_xz_0[j];

            tly_xyy_xxz_0[j] =
                pa_x[j] * tly_yy_xxz_0[j] + fl1_fx * tly_yy_xz_0[j] + 0.5 * fl1_fx * tpz_yy_xxz_0[j] + fl1_fx * fl1_fgb * tdz_yy_xxz_0[j];

            tlz_xyy_xxz_0[j] =
                pa_x[j] * tlz_yy_xxz_0[j] + fl1_fx * tlz_yy_xz_0[j] - 0.5 * fl1_fx * tpy_yy_xxz_0[j] - fl1_fx * fl1_fgb * tdy_yy_xxz_0[j];

            tlx_xyy_xyy_0[j] = pa_x[j] * tlx_yy_xyy_0[j] + 0.5 * fl1_fx * tlx_yy_yy_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForFF_100_150(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_l_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_3_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tly_yy_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 33);

        auto tlz_yy_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 33);

        auto tlx_yy_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 34);

        auto tly_yy_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 34);

        auto tlz_yy_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 34);

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

        auto tly_yy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 21);

        auto tlz_yy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 21);

        auto tlx_yy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 22);

        auto tly_yy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 22);

        auto tlz_yy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 22);

        auto tlx_yy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 23);

        auto tly_yy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 23);

        auto tlz_yy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 23);

        auto tlx_yz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 24);

        auto tly_yz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 24);

        auto tlz_yz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 24);

        auto tlx_yz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 25);

        auto tly_yz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 25);

        auto tlz_yz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 25);

        auto tlx_yz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 26);

        auto tly_yz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 26);

        auto tlz_yz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 26);

        auto tlx_yz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 27);

        auto tly_yz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 27);

        auto tlz_yz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 27);

        auto tlx_yz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 28);

        auto tly_yz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 28);

        auto tlz_yz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 28);

        auto tlx_yz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 29);

        auto tly_yz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 29);

        auto tlz_yz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 29);

        auto tpy_yy_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 33);

        auto tpz_yy_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 33);

        auto tpy_yy_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 34);

        auto tpz_yy_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 34);

        auto tpy_yy_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 35);

        auto tpz_yy_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 35);

        auto tpy_yy_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 36);

        auto tpz_yy_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 36);

        auto tpy_yy_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 37);

        auto tpz_yy_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 37);

        auto tpy_yy_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 38);

        auto tpz_yy_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 38);

        auto tpy_yy_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 39);

        auto tpz_yy_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 39);

        auto tpy_yz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 40);

        auto tpz_yz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 40);

        auto tpy_yz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 41);

        auto tpz_yz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 41);

        auto tpy_yz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 42);

        auto tpz_yz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 42);

        auto tpy_yz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 43);

        auto tpz_yz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 43);

        auto tpy_yz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 44);

        auto tpz_yz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 44);

        auto tpy_yz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 45);

        auto tpz_yz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 45);

        auto tpy_yz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 46);

        auto tpz_yz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 46);

        auto tpy_yz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 47);

        auto tpz_yz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 47);

        auto tpy_yz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 48);

        auto tpz_yz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 48);

        auto tpy_yz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 49);

        auto tpz_yz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 49);

        auto tdy_yy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 33);

        auto tdz_yy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 33);

        auto tdy_yy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 34);

        auto tdz_yy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 34);

        auto tdy_yy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 35);

        auto tdz_yy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 35);

        auto tdy_yy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 36);

        auto tdz_yy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 36);

        auto tdy_yy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 37);

        auto tdz_yy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 37);

        auto tdy_yy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 38);

        auto tdz_yy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 38);

        auto tdy_yy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 39);

        auto tdz_yy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 39);

        auto tdy_yz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 40);

        auto tdz_yz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 40);

        auto tdy_yz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 41);

        auto tdz_yz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 41);

        auto tdy_yz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 42);

        auto tdz_yz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 42);

        auto tdy_yz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 43);

        auto tdz_yz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 43);

        auto tdy_yz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 44);

        auto tdz_yz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 44);

        auto tdy_yz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 45);

        auto tdz_yz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 45);

        auto tdy_yz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 46);

        auto tdz_yz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 46);

        auto tdy_yz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 47);

        auto tdz_yz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 47);

        auto tdy_yz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 48);

        auto tdz_yz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 48);

        auto tdy_yz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 49);

        auto tdz_yz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 49);

        // set up pointers to integrals

        auto tly_xyy_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 33);

        auto tlz_xyy_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 33);

        auto tlx_xyy_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 34);

        auto tly_xyy_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 34);

        auto tlz_xyy_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 34);

        auto tlx_xyy_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 35);

        auto tly_xyy_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 35);

        auto tlz_xyy_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 35);

        auto tlx_xyy_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 36);

        auto tly_xyy_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 36);

        auto tlz_xyy_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 36);

        auto tlx_xyy_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 37);

        auto tly_xyy_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 37);

        auto tlz_xyy_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 37);

        auto tlx_xyy_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 38);

        auto tly_xyy_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 38);

        auto tlz_xyy_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 38);

        auto tlx_xyy_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 39);

        auto tly_xyy_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 39);

        auto tlz_xyy_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 39);

        auto tlx_xyz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 40);

        auto tly_xyz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 40);

        auto tlz_xyz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 40);

        auto tlx_xyz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 41);

        auto tly_xyz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 41);

        auto tlz_xyz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 41);

        auto tlx_xyz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 42);

        auto tly_xyz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 42);

        auto tlz_xyz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 42);

        auto tlx_xyz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 43);

        auto tly_xyz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 43);

        auto tlz_xyz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 43);

        auto tlx_xyz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 44);

        auto tly_xyz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 44);

        auto tlz_xyz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 44);

        auto tlx_xyz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 45);

        auto tly_xyz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 45);

        auto tlz_xyz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 45);

        auto tlx_xyz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 46);

        auto tly_xyz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 46);

        auto tlz_xyz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 46);

        auto tlx_xyz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 47);

        auto tly_xyz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 47);

        auto tlz_xyz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 47);

        auto tlx_xyz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 48);

        auto tly_xyz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 48);

        auto tlz_xyz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 48);

        auto tlx_xyz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 49);

        auto tly_xyz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 49);

        auto tlz_xyz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 49);

        // Batch of Integrals (100,150)

        #pragma omp simd aligned(fgb, fx, pa_x, tdy_yy_xyy_0, tdy_yy_xyz_0, tdy_yy_xzz_0, tdy_yy_yyy_0, \
                                     tdy_yy_yyz_0, tdy_yy_yzz_0, tdy_yy_zzz_0, tdy_yz_xxx_0, tdy_yz_xxy_0, tdy_yz_xxz_0, \
                                     tdy_yz_xyy_0, tdy_yz_xyz_0, tdy_yz_xzz_0, tdy_yz_yyy_0, tdy_yz_yyz_0, tdy_yz_yzz_0, \
                                     tdy_yz_zzz_0, tdz_yy_xyy_0, tdz_yy_xyz_0, tdz_yy_xzz_0, tdz_yy_yyy_0, tdz_yy_yyz_0, \
                                     tdz_yy_yzz_0, tdz_yy_zzz_0, tdz_yz_xxx_0, tdz_yz_xxy_0, tdz_yz_xxz_0, tdz_yz_xyy_0, \
                                     tdz_yz_xyz_0, tdz_yz_xzz_0, tdz_yz_yyy_0, tdz_yz_yyz_0, tdz_yz_yzz_0, tdz_yz_zzz_0, \
                                     tlx_xyy_xyz_0, tlx_xyy_xzz_0, tlx_xyy_yyy_0, tlx_xyy_yyz_0, tlx_xyy_yzz_0, \
                                     tlx_xyy_zzz_0, tlx_xyz_xxx_0, tlx_xyz_xxy_0, tlx_xyz_xxz_0, tlx_xyz_xyy_0, \
                                     tlx_xyz_xyz_0, tlx_xyz_xzz_0, tlx_xyz_yyy_0, tlx_xyz_yyz_0, tlx_xyz_yzz_0, \
                                     tlx_xyz_zzz_0, tlx_yy_xyz_0, tlx_yy_xzz_0, tlx_yy_yyy_0, tlx_yy_yyz_0, tlx_yy_yz_0, \
                                     tlx_yy_yzz_0, tlx_yy_zz_0, tlx_yy_zzz_0, tlx_yz_xx_0, tlx_yz_xxx_0, tlx_yz_xxy_0, \
                                     tlx_yz_xxz_0, tlx_yz_xy_0, tlx_yz_xyy_0, tlx_yz_xyz_0, tlx_yz_xz_0, tlx_yz_xzz_0, \
                                     tlx_yz_yy_0, tlx_yz_yyy_0, tlx_yz_yyz_0, tlx_yz_yz_0, tlx_yz_yzz_0, tlx_yz_zz_0, \
                                     tlx_yz_zzz_0, tly_xyy_xyy_0, tly_xyy_xyz_0, tly_xyy_xzz_0, tly_xyy_yyy_0, \
                                     tly_xyy_yyz_0, tly_xyy_yzz_0, tly_xyy_zzz_0, tly_xyz_xxx_0, tly_xyz_xxy_0, \
                                     tly_xyz_xxz_0, tly_xyz_xyy_0, tly_xyz_xyz_0, tly_xyz_xzz_0, tly_xyz_yyy_0, \
                                     tly_xyz_yyz_0, tly_xyz_yzz_0, tly_xyz_zzz_0, tly_yy_xyy_0, tly_yy_xyz_0, \
                                     tly_yy_xzz_0, tly_yy_yy_0, tly_yy_yyy_0, tly_yy_yyz_0, tly_yy_yz_0, tly_yy_yzz_0, \
                                     tly_yy_zz_0, tly_yy_zzz_0, tly_yz_xx_0, tly_yz_xxx_0, tly_yz_xxy_0, tly_yz_xxz_0, \
                                     tly_yz_xy_0, tly_yz_xyy_0, tly_yz_xyz_0, tly_yz_xz_0, tly_yz_xzz_0, tly_yz_yy_0, \
                                     tly_yz_yyy_0, tly_yz_yyz_0, tly_yz_yz_0, tly_yz_yzz_0, tly_yz_zz_0, tly_yz_zzz_0, \
                                     tlz_xyy_xyy_0, tlz_xyy_xyz_0, tlz_xyy_xzz_0, tlz_xyy_yyy_0, tlz_xyy_yyz_0, \
                                     tlz_xyy_yzz_0, tlz_xyy_zzz_0, tlz_xyz_xxx_0, tlz_xyz_xxy_0, tlz_xyz_xxz_0, \
                                     tlz_xyz_xyy_0, tlz_xyz_xyz_0, tlz_xyz_xzz_0, tlz_xyz_yyy_0, tlz_xyz_yyz_0, \
                                     tlz_xyz_yzz_0, tlz_xyz_zzz_0, tlz_yy_xyy_0, tlz_yy_xyz_0, tlz_yy_xzz_0, tlz_yy_yy_0, \
                                     tlz_yy_yyy_0, tlz_yy_yyz_0, tlz_yy_yz_0, tlz_yy_yzz_0, tlz_yy_zz_0, tlz_yy_zzz_0, \
                                     tlz_yz_xx_0, tlz_yz_xxx_0, tlz_yz_xxy_0, tlz_yz_xxz_0, tlz_yz_xy_0, tlz_yz_xyy_0, \
                                     tlz_yz_xyz_0, tlz_yz_xz_0, tlz_yz_xzz_0, tlz_yz_yy_0, tlz_yz_yyy_0, tlz_yz_yyz_0, \
                                     tlz_yz_yz_0, tlz_yz_yzz_0, tlz_yz_zz_0, tlz_yz_zzz_0, tpy_yy_xyy_0, tpy_yy_xyz_0, \
                                     tpy_yy_xzz_0, tpy_yy_yyy_0, tpy_yy_yyz_0, tpy_yy_yzz_0, tpy_yy_zzz_0, tpy_yz_xxx_0, \
                                     tpy_yz_xxy_0, tpy_yz_xxz_0, tpy_yz_xyy_0, tpy_yz_xyz_0, tpy_yz_xzz_0, tpy_yz_yyy_0, \
                                     tpy_yz_yyz_0, tpy_yz_yzz_0, tpy_yz_zzz_0, tpz_yy_xyy_0, tpz_yy_xyz_0, tpz_yy_xzz_0, \
                                     tpz_yy_yyy_0, tpz_yy_yyz_0, tpz_yy_yzz_0, tpz_yy_zzz_0, tpz_yz_xxx_0, tpz_yz_xxy_0, \
                                     tpz_yz_xxz_0, tpz_yz_xyy_0, tpz_yz_xyz_0, tpz_yz_xzz_0, tpz_yz_yyy_0, tpz_yz_yyz_0, \
                                     tpz_yz_yzz_0, tpz_yz_zzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tly_xyy_xyy_0[j] =
                pa_x[j] * tly_yy_xyy_0[j] + 0.5 * fl1_fx * tly_yy_yy_0[j] + 0.5 * fl1_fx * tpz_yy_xyy_0[j] + fl1_fx * fl1_fgb * tdz_yy_xyy_0[j];

            tlz_xyy_xyy_0[j] =
                pa_x[j] * tlz_yy_xyy_0[j] + 0.5 * fl1_fx * tlz_yy_yy_0[j] - 0.5 * fl1_fx * tpy_yy_xyy_0[j] - fl1_fx * fl1_fgb * tdy_yy_xyy_0[j];

            tlx_xyy_xyz_0[j] = pa_x[j] * tlx_yy_xyz_0[j] + 0.5 * fl1_fx * tlx_yy_yz_0[j];

            tly_xyy_xyz_0[j] =
                pa_x[j] * tly_yy_xyz_0[j] + 0.5 * fl1_fx * tly_yy_yz_0[j] + 0.5 * fl1_fx * tpz_yy_xyz_0[j] + fl1_fx * fl1_fgb * tdz_yy_xyz_0[j];

            tlz_xyy_xyz_0[j] =
                pa_x[j] * tlz_yy_xyz_0[j] + 0.5 * fl1_fx * tlz_yy_yz_0[j] - 0.5 * fl1_fx * tpy_yy_xyz_0[j] - fl1_fx * fl1_fgb * tdy_yy_xyz_0[j];

            tlx_xyy_xzz_0[j] = pa_x[j] * tlx_yy_xzz_0[j] + 0.5 * fl1_fx * tlx_yy_zz_0[j];

            tly_xyy_xzz_0[j] =
                pa_x[j] * tly_yy_xzz_0[j] + 0.5 * fl1_fx * tly_yy_zz_0[j] + 0.5 * fl1_fx * tpz_yy_xzz_0[j] + fl1_fx * fl1_fgb * tdz_yy_xzz_0[j];

            tlz_xyy_xzz_0[j] =
                pa_x[j] * tlz_yy_xzz_0[j] + 0.5 * fl1_fx * tlz_yy_zz_0[j] - 0.5 * fl1_fx * tpy_yy_xzz_0[j] - fl1_fx * fl1_fgb * tdy_yy_xzz_0[j];

            tlx_xyy_yyy_0[j] = pa_x[j] * tlx_yy_yyy_0[j];

            tly_xyy_yyy_0[j] = pa_x[j] * tly_yy_yyy_0[j] + 0.5 * fl1_fx * tpz_yy_yyy_0[j] + fl1_fx * fl1_fgb * tdz_yy_yyy_0[j];

            tlz_xyy_yyy_0[j] = pa_x[j] * tlz_yy_yyy_0[j] - 0.5 * fl1_fx * tpy_yy_yyy_0[j] - fl1_fx * fl1_fgb * tdy_yy_yyy_0[j];

            tlx_xyy_yyz_0[j] = pa_x[j] * tlx_yy_yyz_0[j];

            tly_xyy_yyz_0[j] = pa_x[j] * tly_yy_yyz_0[j] + 0.5 * fl1_fx * tpz_yy_yyz_0[j] + fl1_fx * fl1_fgb * tdz_yy_yyz_0[j];

            tlz_xyy_yyz_0[j] = pa_x[j] * tlz_yy_yyz_0[j] - 0.5 * fl1_fx * tpy_yy_yyz_0[j] - fl1_fx * fl1_fgb * tdy_yy_yyz_0[j];

            tlx_xyy_yzz_0[j] = pa_x[j] * tlx_yy_yzz_0[j];

            tly_xyy_yzz_0[j] = pa_x[j] * tly_yy_yzz_0[j] + 0.5 * fl1_fx * tpz_yy_yzz_0[j] + fl1_fx * fl1_fgb * tdz_yy_yzz_0[j];

            tlz_xyy_yzz_0[j] = pa_x[j] * tlz_yy_yzz_0[j] - 0.5 * fl1_fx * tpy_yy_yzz_0[j] - fl1_fx * fl1_fgb * tdy_yy_yzz_0[j];

            tlx_xyy_zzz_0[j] = pa_x[j] * tlx_yy_zzz_0[j];

            tly_xyy_zzz_0[j] = pa_x[j] * tly_yy_zzz_0[j] + 0.5 * fl1_fx * tpz_yy_zzz_0[j] + fl1_fx * fl1_fgb * tdz_yy_zzz_0[j];

            tlz_xyy_zzz_0[j] = pa_x[j] * tlz_yy_zzz_0[j] - 0.5 * fl1_fx * tpy_yy_zzz_0[j] - fl1_fx * fl1_fgb * tdy_yy_zzz_0[j];

            tlx_xyz_xxx_0[j] = pa_x[j] * tlx_yz_xxx_0[j] + 1.5 * fl1_fx * tlx_yz_xx_0[j];

            tly_xyz_xxx_0[j] =
                pa_x[j] * tly_yz_xxx_0[j] + 1.5 * fl1_fx * tly_yz_xx_0[j] + 0.5 * fl1_fx * tpz_yz_xxx_0[j] + fl1_fx * fl1_fgb * tdz_yz_xxx_0[j];

            tlz_xyz_xxx_0[j] =
                pa_x[j] * tlz_yz_xxx_0[j] + 1.5 * fl1_fx * tlz_yz_xx_0[j] - 0.5 * fl1_fx * tpy_yz_xxx_0[j] - fl1_fx * fl1_fgb * tdy_yz_xxx_0[j];

            tlx_xyz_xxy_0[j] = pa_x[j] * tlx_yz_xxy_0[j] + fl1_fx * tlx_yz_xy_0[j];

            tly_xyz_xxy_0[j] =
                pa_x[j] * tly_yz_xxy_0[j] + fl1_fx * tly_yz_xy_0[j] + 0.5 * fl1_fx * tpz_yz_xxy_0[j] + fl1_fx * fl1_fgb * tdz_yz_xxy_0[j];

            tlz_xyz_xxy_0[j] =
                pa_x[j] * tlz_yz_xxy_0[j] + fl1_fx * tlz_yz_xy_0[j] - 0.5 * fl1_fx * tpy_yz_xxy_0[j] - fl1_fx * fl1_fgb * tdy_yz_xxy_0[j];

            tlx_xyz_xxz_0[j] = pa_x[j] * tlx_yz_xxz_0[j] + fl1_fx * tlx_yz_xz_0[j];

            tly_xyz_xxz_0[j] =
                pa_x[j] * tly_yz_xxz_0[j] + fl1_fx * tly_yz_xz_0[j] + 0.5 * fl1_fx * tpz_yz_xxz_0[j] + fl1_fx * fl1_fgb * tdz_yz_xxz_0[j];

            tlz_xyz_xxz_0[j] =
                pa_x[j] * tlz_yz_xxz_0[j] + fl1_fx * tlz_yz_xz_0[j] - 0.5 * fl1_fx * tpy_yz_xxz_0[j] - fl1_fx * fl1_fgb * tdy_yz_xxz_0[j];

            tlx_xyz_xyy_0[j] = pa_x[j] * tlx_yz_xyy_0[j] + 0.5 * fl1_fx * tlx_yz_yy_0[j];

            tly_xyz_xyy_0[j] =
                pa_x[j] * tly_yz_xyy_0[j] + 0.5 * fl1_fx * tly_yz_yy_0[j] + 0.5 * fl1_fx * tpz_yz_xyy_0[j] + fl1_fx * fl1_fgb * tdz_yz_xyy_0[j];

            tlz_xyz_xyy_0[j] =
                pa_x[j] * tlz_yz_xyy_0[j] + 0.5 * fl1_fx * tlz_yz_yy_0[j] - 0.5 * fl1_fx * tpy_yz_xyy_0[j] - fl1_fx * fl1_fgb * tdy_yz_xyy_0[j];

            tlx_xyz_xyz_0[j] = pa_x[j] * tlx_yz_xyz_0[j] + 0.5 * fl1_fx * tlx_yz_yz_0[j];

            tly_xyz_xyz_0[j] =
                pa_x[j] * tly_yz_xyz_0[j] + 0.5 * fl1_fx * tly_yz_yz_0[j] + 0.5 * fl1_fx * tpz_yz_xyz_0[j] + fl1_fx * fl1_fgb * tdz_yz_xyz_0[j];

            tlz_xyz_xyz_0[j] =
                pa_x[j] * tlz_yz_xyz_0[j] + 0.5 * fl1_fx * tlz_yz_yz_0[j] - 0.5 * fl1_fx * tpy_yz_xyz_0[j] - fl1_fx * fl1_fgb * tdy_yz_xyz_0[j];

            tlx_xyz_xzz_0[j] = pa_x[j] * tlx_yz_xzz_0[j] + 0.5 * fl1_fx * tlx_yz_zz_0[j];

            tly_xyz_xzz_0[j] =
                pa_x[j] * tly_yz_xzz_0[j] + 0.5 * fl1_fx * tly_yz_zz_0[j] + 0.5 * fl1_fx * tpz_yz_xzz_0[j] + fl1_fx * fl1_fgb * tdz_yz_xzz_0[j];

            tlz_xyz_xzz_0[j] =
                pa_x[j] * tlz_yz_xzz_0[j] + 0.5 * fl1_fx * tlz_yz_zz_0[j] - 0.5 * fl1_fx * tpy_yz_xzz_0[j] - fl1_fx * fl1_fgb * tdy_yz_xzz_0[j];

            tlx_xyz_yyy_0[j] = pa_x[j] * tlx_yz_yyy_0[j];

            tly_xyz_yyy_0[j] = pa_x[j] * tly_yz_yyy_0[j] + 0.5 * fl1_fx * tpz_yz_yyy_0[j] + fl1_fx * fl1_fgb * tdz_yz_yyy_0[j];

            tlz_xyz_yyy_0[j] = pa_x[j] * tlz_yz_yyy_0[j] - 0.5 * fl1_fx * tpy_yz_yyy_0[j] - fl1_fx * fl1_fgb * tdy_yz_yyy_0[j];

            tlx_xyz_yyz_0[j] = pa_x[j] * tlx_yz_yyz_0[j];

            tly_xyz_yyz_0[j] = pa_x[j] * tly_yz_yyz_0[j] + 0.5 * fl1_fx * tpz_yz_yyz_0[j] + fl1_fx * fl1_fgb * tdz_yz_yyz_0[j];

            tlz_xyz_yyz_0[j] = pa_x[j] * tlz_yz_yyz_0[j] - 0.5 * fl1_fx * tpy_yz_yyz_0[j] - fl1_fx * fl1_fgb * tdy_yz_yyz_0[j];

            tlx_xyz_yzz_0[j] = pa_x[j] * tlx_yz_yzz_0[j];

            tly_xyz_yzz_0[j] = pa_x[j] * tly_yz_yzz_0[j] + 0.5 * fl1_fx * tpz_yz_yzz_0[j] + fl1_fx * fl1_fgb * tdz_yz_yzz_0[j];

            tlz_xyz_yzz_0[j] = pa_x[j] * tlz_yz_yzz_0[j] - 0.5 * fl1_fx * tpy_yz_yzz_0[j] - fl1_fx * fl1_fgb * tdy_yz_yzz_0[j];

            tlx_xyz_zzz_0[j] = pa_x[j] * tlx_yz_zzz_0[j];

            tly_xyz_zzz_0[j] = pa_x[j] * tly_yz_zzz_0[j] + 0.5 * fl1_fx * tpz_yz_zzz_0[j] + fl1_fx * fl1_fgb * tdz_yz_zzz_0[j];

            tlz_xyz_zzz_0[j] = pa_x[j] * tlz_yz_zzz_0[j] - 0.5 * fl1_fx * tpy_yz_zzz_0[j] - fl1_fx * fl1_fgb * tdy_yz_zzz_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForFF_150_200(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_l_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_3_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tlx_yy_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 36);

        auto tly_yy_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 36);

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

        auto tlx_y_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 15);

        auto tly_y_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 15);

        auto tlz_y_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 15);

        auto tlx_y_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 16);

        auto tly_y_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 16);

        auto tlx_yy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 18);

        auto tly_yy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 18);

        auto tlz_yy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 18);

        auto tlx_yy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 19);

        auto tly_yy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 19);

        auto tlz_yy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 19);

        auto tlx_yy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 20);

        auto tly_yy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 20);

        auto tlz_yy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 20);

        auto tlx_yy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 21);

        auto tly_yy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 21);

        auto tlx_zz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 30);

        auto tly_zz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 30);

        auto tlz_zz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 30);

        auto tlx_zz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 31);

        auto tly_zz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 31);

        auto tlz_zz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 31);

        auto tlx_zz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 32);

        auto tly_zz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 32);

        auto tlz_zz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 32);

        auto tlx_zz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 33);

        auto tly_zz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 33);

        auto tlz_zz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 33);

        auto tlx_zz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 34);

        auto tly_zz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 34);

        auto tlz_zz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 34);

        auto tlx_zz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 35);

        auto tly_zz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 35);

        auto tlz_zz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 35);

        auto tpx_yy_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 30);

        auto tpz_yy_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 30);

        auto tpx_yy_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 31);

        auto tpz_yy_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 31);

        auto tpx_yy_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 32);

        auto tpz_yy_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 32);

        auto tpx_yy_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 33);

        auto tpz_yy_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 33);

        auto tpx_yy_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 34);

        auto tpz_yy_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 34);

        auto tpx_yy_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 35);

        auto tpz_yy_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 35);

        auto tpz_yy_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 36);

        auto tpy_zz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 50);

        auto tpz_zz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 50);

        auto tpy_zz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 51);

        auto tpz_zz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 51);

        auto tpy_zz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 52);

        auto tpz_zz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 52);

        auto tpy_zz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 53);

        auto tpz_zz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 53);

        auto tpy_zz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 54);

        auto tpz_zz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 54);

        auto tpy_zz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 55);

        auto tpz_zz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 55);

        auto tpy_zz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 56);

        auto tpz_zz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 56);

        auto tpy_zz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 57);

        auto tpz_zz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 57);

        auto tpy_zz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 58);

        auto tpz_zz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 58);

        auto tpy_zz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 59);

        auto tpz_zz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 59);

        auto tdx_yy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 30);

        auto tdz_yy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 30);

        auto tdx_yy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 31);

        auto tdz_yy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 31);

        auto tdx_yy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 32);

        auto tdz_yy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 32);

        auto tdx_yy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 33);

        auto tdz_yy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 33);

        auto tdx_yy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 34);

        auto tdz_yy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 34);

        auto tdx_yy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 35);

        auto tdz_yy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 35);

        auto tdz_yy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 36);

        auto tdy_zz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 50);

        auto tdz_zz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 50);

        auto tdy_zz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 51);

        auto tdz_zz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 51);

        auto tdy_zz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 52);

        auto tdz_zz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 52);

        auto tdy_zz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 53);

        auto tdz_zz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 53);

        auto tdy_zz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 54);

        auto tdz_zz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 54);

        auto tdy_zz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 55);

        auto tdz_zz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 55);

        auto tdy_zz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 56);

        auto tdz_zz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 56);

        auto tdy_zz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 57);

        auto tdz_zz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 57);

        auto tdy_zz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 58);

        auto tdz_zz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 58);

        auto tdy_zz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 59);

        auto tdz_zz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 59);

        // set up pointers to integrals

        auto tlx_xzz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 50);

        auto tly_xzz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 50);

        auto tlz_xzz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 50);

        auto tlx_xzz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 51);

        auto tly_xzz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 51);

        auto tlz_xzz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 51);

        auto tlx_xzz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 52);

        auto tly_xzz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 52);

        auto tlz_xzz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 52);

        auto tlx_xzz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 53);

        auto tly_xzz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 53);

        auto tlz_xzz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 53);

        auto tlx_xzz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 54);

        auto tly_xzz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 54);

        auto tlz_xzz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 54);

        auto tlx_xzz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 55);

        auto tly_xzz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 55);

        auto tlz_xzz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 55);

        auto tlx_xzz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 56);

        auto tly_xzz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 56);

        auto tlz_xzz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 56);

        auto tlx_xzz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 57);

        auto tly_xzz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 57);

        auto tlz_xzz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 57);

        auto tlx_xzz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 58);

        auto tly_xzz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 58);

        auto tlz_xzz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 58);

        auto tlx_xzz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 59);

        auto tly_xzz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 59);

        auto tlz_xzz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 59);

        auto tlx_yyy_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 60);

        auto tly_yyy_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 60);

        auto tlz_yyy_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 60);

        auto tlx_yyy_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 61);

        auto tly_yyy_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 61);

        auto tlz_yyy_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 61);

        auto tlx_yyy_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 62);

        auto tly_yyy_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 62);

        auto tlz_yyy_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 62);

        auto tlx_yyy_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 63);

        auto tly_yyy_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 63);

        auto tlz_yyy_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 63);

        auto tlx_yyy_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 64);

        auto tly_yyy_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 64);

        auto tlz_yyy_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 64);

        auto tlx_yyy_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 65);

        auto tly_yyy_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 65);

        auto tlz_yyy_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 65);

        auto tlx_yyy_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 66);

        auto tly_yyy_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 66);

        // Batch of Integrals (150,200)

        #pragma omp simd aligned(fgb, fx, pa_x, pa_y, tdx_yy_xxx_0, tdx_yy_xxy_0, tdx_yy_xxz_0, \
                                     tdx_yy_xyy_0, tdx_yy_xyz_0, tdx_yy_xzz_0, tdy_zz_xxx_0, tdy_zz_xxy_0, tdy_zz_xxz_0, \
                                     tdy_zz_xyy_0, tdy_zz_xyz_0, tdy_zz_xzz_0, tdy_zz_yyy_0, tdy_zz_yyz_0, tdy_zz_yzz_0, \
                                     tdy_zz_zzz_0, tdz_yy_xxx_0, tdz_yy_xxy_0, tdz_yy_xxz_0, tdz_yy_xyy_0, tdz_yy_xyz_0, \
                                     tdz_yy_xzz_0, tdz_yy_yyy_0, tdz_zz_xxx_0, tdz_zz_xxy_0, tdz_zz_xxz_0, tdz_zz_xyy_0, \
                                     tdz_zz_xyz_0, tdz_zz_xzz_0, tdz_zz_yyy_0, tdz_zz_yyz_0, tdz_zz_yzz_0, tdz_zz_zzz_0, \
                                     tlx_xzz_xxx_0, tlx_xzz_xxy_0, tlx_xzz_xxz_0, tlx_xzz_xyy_0, tlx_xzz_xyz_0, \
                                     tlx_xzz_xzz_0, tlx_xzz_yyy_0, tlx_xzz_yyz_0, tlx_xzz_yzz_0, tlx_xzz_zzz_0, \
                                     tlx_y_xxx_0, tlx_y_xxy_0, tlx_y_xxz_0, tlx_y_xyy_0, tlx_y_xyz_0, tlx_y_xzz_0, \
                                     tlx_y_yyy_0, tlx_yy_xx_0, tlx_yy_xxx_0, tlx_yy_xxy_0, tlx_yy_xxz_0, tlx_yy_xy_0, \
                                     tlx_yy_xyy_0, tlx_yy_xyz_0, tlx_yy_xz_0, tlx_yy_xzz_0, tlx_yy_yy_0, tlx_yy_yyy_0, \
                                     tlx_yyy_xxx_0, tlx_yyy_xxy_0, tlx_yyy_xxz_0, tlx_yyy_xyy_0, tlx_yyy_xyz_0, \
                                     tlx_yyy_xzz_0, tlx_yyy_yyy_0, tlx_zz_xx_0, tlx_zz_xxx_0, tlx_zz_xxy_0, tlx_zz_xxz_0, \
                                     tlx_zz_xy_0, tlx_zz_xyy_0, tlx_zz_xyz_0, tlx_zz_xz_0, tlx_zz_xzz_0, tlx_zz_yy_0, \
                                     tlx_zz_yyy_0, tlx_zz_yyz_0, tlx_zz_yz_0, tlx_zz_yzz_0, tlx_zz_zz_0, tlx_zz_zzz_0, \
                                     tly_xzz_xxx_0, tly_xzz_xxy_0, tly_xzz_xxz_0, tly_xzz_xyy_0, tly_xzz_xyz_0, \
                                     tly_xzz_xzz_0, tly_xzz_yyy_0, tly_xzz_yyz_0, tly_xzz_yzz_0, tly_xzz_zzz_0, \
                                     tly_y_xxx_0, tly_y_xxy_0, tly_y_xxz_0, tly_y_xyy_0, tly_y_xyz_0, tly_y_xzz_0, \
                                     tly_y_yyy_0, tly_yy_xx_0, tly_yy_xxx_0, tly_yy_xxy_0, tly_yy_xxz_0, tly_yy_xy_0, \
                                     tly_yy_xyy_0, tly_yy_xyz_0, tly_yy_xz_0, tly_yy_xzz_0, tly_yy_yy_0, tly_yy_yyy_0, \
                                     tly_yyy_xxx_0, tly_yyy_xxy_0, tly_yyy_xxz_0, tly_yyy_xyy_0, tly_yyy_xyz_0, \
                                     tly_yyy_xzz_0, tly_yyy_yyy_0, tly_zz_xx_0, tly_zz_xxx_0, tly_zz_xxy_0, tly_zz_xxz_0, \
                                     tly_zz_xy_0, tly_zz_xyy_0, tly_zz_xyz_0, tly_zz_xz_0, tly_zz_xzz_0, tly_zz_yy_0, \
                                     tly_zz_yyy_0, tly_zz_yyz_0, tly_zz_yz_0, tly_zz_yzz_0, tly_zz_zz_0, tly_zz_zzz_0, \
                                     tlz_xzz_xxx_0, tlz_xzz_xxy_0, tlz_xzz_xxz_0, tlz_xzz_xyy_0, tlz_xzz_xyz_0, \
                                     tlz_xzz_xzz_0, tlz_xzz_yyy_0, tlz_xzz_yyz_0, tlz_xzz_yzz_0, tlz_xzz_zzz_0, \
                                     tlz_y_xxx_0, tlz_y_xxy_0, tlz_y_xxz_0, tlz_y_xyy_0, tlz_y_xyz_0, tlz_y_xzz_0, \
                                     tlz_yy_xx_0, tlz_yy_xxx_0, tlz_yy_xxy_0, tlz_yy_xxz_0, tlz_yy_xy_0, tlz_yy_xyy_0, \
                                     tlz_yy_xyz_0, tlz_yy_xz_0, tlz_yy_xzz_0, tlz_yyy_xxx_0, tlz_yyy_xxy_0, \
                                     tlz_yyy_xxz_0, tlz_yyy_xyy_0, tlz_yyy_xyz_0, tlz_yyy_xzz_0, tlz_zz_xx_0, \
                                     tlz_zz_xxx_0, tlz_zz_xxy_0, tlz_zz_xxz_0, tlz_zz_xy_0, tlz_zz_xyy_0, tlz_zz_xyz_0, \
                                     tlz_zz_xz_0, tlz_zz_xzz_0, tlz_zz_yy_0, tlz_zz_yyy_0, tlz_zz_yyz_0, tlz_zz_yz_0, \
                                     tlz_zz_yzz_0, tlz_zz_zz_0, tlz_zz_zzz_0, tpx_yy_xxx_0, tpx_yy_xxy_0, tpx_yy_xxz_0, \
                                     tpx_yy_xyy_0, tpx_yy_xyz_0, tpx_yy_xzz_0, tpy_zz_xxx_0, tpy_zz_xxy_0, tpy_zz_xxz_0, \
                                     tpy_zz_xyy_0, tpy_zz_xyz_0, tpy_zz_xzz_0, tpy_zz_yyy_0, tpy_zz_yyz_0, tpy_zz_yzz_0, \
                                     tpy_zz_zzz_0, tpz_yy_xxx_0, tpz_yy_xxy_0, tpz_yy_xxz_0, tpz_yy_xyy_0, tpz_yy_xyz_0, \
                                     tpz_yy_xzz_0, tpz_yy_yyy_0, tpz_zz_xxx_0, tpz_zz_xxy_0, tpz_zz_xxz_0, tpz_zz_xyy_0, \
                                     tpz_zz_xyz_0, tpz_zz_xzz_0, tpz_zz_yyy_0, tpz_zz_yyz_0, tpz_zz_yzz_0, tpz_zz_zzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_xzz_xxx_0[j] = pa_x[j] * tlx_zz_xxx_0[j] + 1.5 * fl1_fx * tlx_zz_xx_0[j];

            tly_xzz_xxx_0[j] =
                pa_x[j] * tly_zz_xxx_0[j] + 1.5 * fl1_fx * tly_zz_xx_0[j] + 0.5 * fl1_fx * tpz_zz_xxx_0[j] + fl1_fx * fl1_fgb * tdz_zz_xxx_0[j];

            tlz_xzz_xxx_0[j] =
                pa_x[j] * tlz_zz_xxx_0[j] + 1.5 * fl1_fx * tlz_zz_xx_0[j] - 0.5 * fl1_fx * tpy_zz_xxx_0[j] - fl1_fx * fl1_fgb * tdy_zz_xxx_0[j];

            tlx_xzz_xxy_0[j] = pa_x[j] * tlx_zz_xxy_0[j] + fl1_fx * tlx_zz_xy_0[j];

            tly_xzz_xxy_0[j] =
                pa_x[j] * tly_zz_xxy_0[j] + fl1_fx * tly_zz_xy_0[j] + 0.5 * fl1_fx * tpz_zz_xxy_0[j] + fl1_fx * fl1_fgb * tdz_zz_xxy_0[j];

            tlz_xzz_xxy_0[j] =
                pa_x[j] * tlz_zz_xxy_0[j] + fl1_fx * tlz_zz_xy_0[j] - 0.5 * fl1_fx * tpy_zz_xxy_0[j] - fl1_fx * fl1_fgb * tdy_zz_xxy_0[j];

            tlx_xzz_xxz_0[j] = pa_x[j] * tlx_zz_xxz_0[j] + fl1_fx * tlx_zz_xz_0[j];

            tly_xzz_xxz_0[j] =
                pa_x[j] * tly_zz_xxz_0[j] + fl1_fx * tly_zz_xz_0[j] + 0.5 * fl1_fx * tpz_zz_xxz_0[j] + fl1_fx * fl1_fgb * tdz_zz_xxz_0[j];

            tlz_xzz_xxz_0[j] =
                pa_x[j] * tlz_zz_xxz_0[j] + fl1_fx * tlz_zz_xz_0[j] - 0.5 * fl1_fx * tpy_zz_xxz_0[j] - fl1_fx * fl1_fgb * tdy_zz_xxz_0[j];

            tlx_xzz_xyy_0[j] = pa_x[j] * tlx_zz_xyy_0[j] + 0.5 * fl1_fx * tlx_zz_yy_0[j];

            tly_xzz_xyy_0[j] =
                pa_x[j] * tly_zz_xyy_0[j] + 0.5 * fl1_fx * tly_zz_yy_0[j] + 0.5 * fl1_fx * tpz_zz_xyy_0[j] + fl1_fx * fl1_fgb * tdz_zz_xyy_0[j];

            tlz_xzz_xyy_0[j] =
                pa_x[j] * tlz_zz_xyy_0[j] + 0.5 * fl1_fx * tlz_zz_yy_0[j] - 0.5 * fl1_fx * tpy_zz_xyy_0[j] - fl1_fx * fl1_fgb * tdy_zz_xyy_0[j];

            tlx_xzz_xyz_0[j] = pa_x[j] * tlx_zz_xyz_0[j] + 0.5 * fl1_fx * tlx_zz_yz_0[j];

            tly_xzz_xyz_0[j] =
                pa_x[j] * tly_zz_xyz_0[j] + 0.5 * fl1_fx * tly_zz_yz_0[j] + 0.5 * fl1_fx * tpz_zz_xyz_0[j] + fl1_fx * fl1_fgb * tdz_zz_xyz_0[j];

            tlz_xzz_xyz_0[j] =
                pa_x[j] * tlz_zz_xyz_0[j] + 0.5 * fl1_fx * tlz_zz_yz_0[j] - 0.5 * fl1_fx * tpy_zz_xyz_0[j] - fl1_fx * fl1_fgb * tdy_zz_xyz_0[j];

            tlx_xzz_xzz_0[j] = pa_x[j] * tlx_zz_xzz_0[j] + 0.5 * fl1_fx * tlx_zz_zz_0[j];

            tly_xzz_xzz_0[j] =
                pa_x[j] * tly_zz_xzz_0[j] + 0.5 * fl1_fx * tly_zz_zz_0[j] + 0.5 * fl1_fx * tpz_zz_xzz_0[j] + fl1_fx * fl1_fgb * tdz_zz_xzz_0[j];

            tlz_xzz_xzz_0[j] =
                pa_x[j] * tlz_zz_xzz_0[j] + 0.5 * fl1_fx * tlz_zz_zz_0[j] - 0.5 * fl1_fx * tpy_zz_xzz_0[j] - fl1_fx * fl1_fgb * tdy_zz_xzz_0[j];

            tlx_xzz_yyy_0[j] = pa_x[j] * tlx_zz_yyy_0[j];

            tly_xzz_yyy_0[j] = pa_x[j] * tly_zz_yyy_0[j] + 0.5 * fl1_fx * tpz_zz_yyy_0[j] + fl1_fx * fl1_fgb * tdz_zz_yyy_0[j];

            tlz_xzz_yyy_0[j] = pa_x[j] * tlz_zz_yyy_0[j] - 0.5 * fl1_fx * tpy_zz_yyy_0[j] - fl1_fx * fl1_fgb * tdy_zz_yyy_0[j];

            tlx_xzz_yyz_0[j] = pa_x[j] * tlx_zz_yyz_0[j];

            tly_xzz_yyz_0[j] = pa_x[j] * tly_zz_yyz_0[j] + 0.5 * fl1_fx * tpz_zz_yyz_0[j] + fl1_fx * fl1_fgb * tdz_zz_yyz_0[j];

            tlz_xzz_yyz_0[j] = pa_x[j] * tlz_zz_yyz_0[j] - 0.5 * fl1_fx * tpy_zz_yyz_0[j] - fl1_fx * fl1_fgb * tdy_zz_yyz_0[j];

            tlx_xzz_yzz_0[j] = pa_x[j] * tlx_zz_yzz_0[j];

            tly_xzz_yzz_0[j] = pa_x[j] * tly_zz_yzz_0[j] + 0.5 * fl1_fx * tpz_zz_yzz_0[j] + fl1_fx * fl1_fgb * tdz_zz_yzz_0[j];

            tlz_xzz_yzz_0[j] = pa_x[j] * tlz_zz_yzz_0[j] - 0.5 * fl1_fx * tpy_zz_yzz_0[j] - fl1_fx * fl1_fgb * tdy_zz_yzz_0[j];

            tlx_xzz_zzz_0[j] = pa_x[j] * tlx_zz_zzz_0[j];

            tly_xzz_zzz_0[j] = pa_x[j] * tly_zz_zzz_0[j] + 0.5 * fl1_fx * tpz_zz_zzz_0[j] + fl1_fx * fl1_fgb * tdz_zz_zzz_0[j];

            tlz_xzz_zzz_0[j] = pa_x[j] * tlz_zz_zzz_0[j] - 0.5 * fl1_fx * tpy_zz_zzz_0[j] - fl1_fx * fl1_fgb * tdy_zz_zzz_0[j];

            tlx_yyy_xxx_0[j] =
                pa_y[j] * tlx_yy_xxx_0[j] + fl1_fx * tlx_y_xxx_0[j] - 0.5 * fl1_fx * tpz_yy_xxx_0[j] - fl1_fx * fl1_fgb * tdz_yy_xxx_0[j];

            tly_yyy_xxx_0[j] = pa_y[j] * tly_yy_xxx_0[j] + fl1_fx * tly_y_xxx_0[j];

            tlz_yyy_xxx_0[j] =
                pa_y[j] * tlz_yy_xxx_0[j] + fl1_fx * tlz_y_xxx_0[j] + 0.5 * fl1_fx * tpx_yy_xxx_0[j] + fl1_fx * fl1_fgb * tdx_yy_xxx_0[j];

            tlx_yyy_xxy_0[j] = pa_y[j] * tlx_yy_xxy_0[j] + fl1_fx * tlx_y_xxy_0[j] + 0.5 * fl1_fx * tlx_yy_xx_0[j] - 0.5 * fl1_fx * tpz_yy_xxy_0[j] -
                               fl1_fx * fl1_fgb * tdz_yy_xxy_0[j];

            tly_yyy_xxy_0[j] = pa_y[j] * tly_yy_xxy_0[j] + fl1_fx * tly_y_xxy_0[j] + 0.5 * fl1_fx * tly_yy_xx_0[j];

            tlz_yyy_xxy_0[j] = pa_y[j] * tlz_yy_xxy_0[j] + fl1_fx * tlz_y_xxy_0[j] + 0.5 * fl1_fx * tlz_yy_xx_0[j] + 0.5 * fl1_fx * tpx_yy_xxy_0[j] +
                               fl1_fx * fl1_fgb * tdx_yy_xxy_0[j];

            tlx_yyy_xxz_0[j] =
                pa_y[j] * tlx_yy_xxz_0[j] + fl1_fx * tlx_y_xxz_0[j] - 0.5 * fl1_fx * tpz_yy_xxz_0[j] - fl1_fx * fl1_fgb * tdz_yy_xxz_0[j];

            tly_yyy_xxz_0[j] = pa_y[j] * tly_yy_xxz_0[j] + fl1_fx * tly_y_xxz_0[j];

            tlz_yyy_xxz_0[j] =
                pa_y[j] * tlz_yy_xxz_0[j] + fl1_fx * tlz_y_xxz_0[j] + 0.5 * fl1_fx * tpx_yy_xxz_0[j] + fl1_fx * fl1_fgb * tdx_yy_xxz_0[j];

            tlx_yyy_xyy_0[j] = pa_y[j] * tlx_yy_xyy_0[j] + fl1_fx * tlx_y_xyy_0[j] + fl1_fx * tlx_yy_xy_0[j] - 0.5 * fl1_fx * tpz_yy_xyy_0[j] -
                               fl1_fx * fl1_fgb * tdz_yy_xyy_0[j];

            tly_yyy_xyy_0[j] = pa_y[j] * tly_yy_xyy_0[j] + fl1_fx * tly_y_xyy_0[j] + fl1_fx * tly_yy_xy_0[j];

            tlz_yyy_xyy_0[j] = pa_y[j] * tlz_yy_xyy_0[j] + fl1_fx * tlz_y_xyy_0[j] + fl1_fx * tlz_yy_xy_0[j] + 0.5 * fl1_fx * tpx_yy_xyy_0[j] +
                               fl1_fx * fl1_fgb * tdx_yy_xyy_0[j];

            tlx_yyy_xyz_0[j] = pa_y[j] * tlx_yy_xyz_0[j] + fl1_fx * tlx_y_xyz_0[j] + 0.5 * fl1_fx * tlx_yy_xz_0[j] - 0.5 * fl1_fx * tpz_yy_xyz_0[j] -
                               fl1_fx * fl1_fgb * tdz_yy_xyz_0[j];

            tly_yyy_xyz_0[j] = pa_y[j] * tly_yy_xyz_0[j] + fl1_fx * tly_y_xyz_0[j] + 0.5 * fl1_fx * tly_yy_xz_0[j];

            tlz_yyy_xyz_0[j] = pa_y[j] * tlz_yy_xyz_0[j] + fl1_fx * tlz_y_xyz_0[j] + 0.5 * fl1_fx * tlz_yy_xz_0[j] + 0.5 * fl1_fx * tpx_yy_xyz_0[j] +
                               fl1_fx * fl1_fgb * tdx_yy_xyz_0[j];

            tlx_yyy_xzz_0[j] =
                pa_y[j] * tlx_yy_xzz_0[j] + fl1_fx * tlx_y_xzz_0[j] - 0.5 * fl1_fx * tpz_yy_xzz_0[j] - fl1_fx * fl1_fgb * tdz_yy_xzz_0[j];

            tly_yyy_xzz_0[j] = pa_y[j] * tly_yy_xzz_0[j] + fl1_fx * tly_y_xzz_0[j];

            tlz_yyy_xzz_0[j] =
                pa_y[j] * tlz_yy_xzz_0[j] + fl1_fx * tlz_y_xzz_0[j] + 0.5 * fl1_fx * tpx_yy_xzz_0[j] + fl1_fx * fl1_fgb * tdx_yy_xzz_0[j];

            tlx_yyy_yyy_0[j] = pa_y[j] * tlx_yy_yyy_0[j] + fl1_fx * tlx_y_yyy_0[j] + 1.5 * fl1_fx * tlx_yy_yy_0[j] - 0.5 * fl1_fx * tpz_yy_yyy_0[j] -
                               fl1_fx * fl1_fgb * tdz_yy_yyy_0[j];

            tly_yyy_yyy_0[j] = pa_y[j] * tly_yy_yyy_0[j] + fl1_fx * tly_y_yyy_0[j] + 1.5 * fl1_fx * tly_yy_yy_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForFF_200_250(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_l_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_3_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tlz_yy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 21);

        auto tlx_yy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 22);

        auto tly_yy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 22);

        auto tlz_yy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 22);

        auto tlx_yy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 23);

        auto tly_yy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 23);

        auto tlz_yy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 23);

        auto tlx_yz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 24);

        auto tly_yz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 24);

        auto tlz_yz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 24);

        auto tlx_yz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 25);

        auto tly_yz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 25);

        auto tlz_yz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 25);

        auto tlx_yz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 26);

        auto tly_yz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 26);

        auto tlz_yz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 26);

        auto tlx_yz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 27);

        auto tly_yz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 27);

        auto tlz_yz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 27);

        auto tlx_yz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 28);

        auto tly_yz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 28);

        auto tlz_yz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 28);

        auto tlx_yz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 29);

        auto tly_yz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 29);

        auto tlz_yz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 29);

        auto tlx_zz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 30);

        auto tly_zz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 30);

        auto tlz_zz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 30);

        auto tlx_zz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 31);

        auto tpx_yy_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 36);

        auto tpx_yy_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 37);

        auto tpz_yy_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 37);

        auto tpx_yy_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 38);

        auto tpz_yy_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 38);

        auto tpx_yy_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 39);

        auto tpz_yy_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 39);

        auto tpx_yz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 40);

        auto tpz_yz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 40);

        auto tpx_yz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 41);

        auto tpz_yz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 41);

        auto tpx_yz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 42);

        auto tpz_yz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 42);

        auto tpx_yz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 43);

        auto tpz_yz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 43);

        auto tpx_yz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 44);

        auto tpz_yz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 44);

        auto tpx_yz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 45);

        auto tpz_yz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 45);

        auto tpx_yz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 46);

        auto tpz_yz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 46);

        auto tpx_yz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 47);

        auto tpz_yz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 47);

        auto tpx_yz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 48);

        auto tpz_yz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 48);

        auto tpx_yz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 49);

        auto tpz_yz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 49);

        auto tpx_zz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 50);

        auto tpz_zz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 50);

        auto tpx_zz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 51);

        auto tpz_zz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 51);

        auto tpx_zz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 52);

        auto tpz_zz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 52);

        auto tpz_zz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 53);

        auto tdx_yy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 36);

        auto tdx_yy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 37);

        auto tdz_yy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 37);

        auto tdx_yy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 38);

        auto tdz_yy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 38);

        auto tdx_yy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 39);

        auto tdz_yy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 39);

        auto tdx_yz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 40);

        auto tdz_yz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 40);

        auto tdx_yz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 41);

        auto tdz_yz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 41);

        auto tdx_yz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 42);

        auto tdz_yz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 42);

        auto tdx_yz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 43);

        auto tdz_yz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 43);

        auto tdx_yz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 44);

        auto tdz_yz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 44);

        auto tdx_yz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 45);

        auto tdz_yz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 45);

        auto tdx_yz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 46);

        auto tdz_yz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 46);

        auto tdx_yz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 47);

        auto tdz_yz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 47);

        auto tdx_yz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 48);

        auto tdz_yz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 48);

        auto tdx_yz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 49);

        auto tdz_yz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 49);

        auto tdx_zz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 50);

        auto tdz_zz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 50);

        auto tdx_zz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 51);

        auto tdz_zz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 51);

        auto tdx_zz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 52);

        auto tdz_zz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 52);

        auto tdz_zz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 53);

        // set up pointers to integrals

        auto tlz_yyy_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 66);

        auto tlx_yyy_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 67);

        auto tly_yyy_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 67);

        auto tlz_yyy_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 67);

        auto tlx_yyy_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 68);

        auto tly_yyy_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 68);

        auto tlz_yyy_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 68);

        auto tlx_yyy_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 69);

        auto tly_yyy_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 69);

        auto tlz_yyy_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 69);

        auto tlx_yyz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 70);

        auto tly_yyz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 70);

        auto tlz_yyz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 70);

        auto tlx_yyz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 71);

        auto tly_yyz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 71);

        auto tlz_yyz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 71);

        auto tlx_yyz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 72);

        auto tly_yyz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 72);

        auto tlz_yyz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 72);

        auto tlx_yyz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 73);

        auto tly_yyz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 73);

        auto tlz_yyz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 73);

        auto tlx_yyz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 74);

        auto tly_yyz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 74);

        auto tlz_yyz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 74);

        auto tlx_yyz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 75);

        auto tly_yyz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 75);

        auto tlz_yyz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 75);

        auto tlx_yyz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 76);

        auto tly_yyz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 76);

        auto tlz_yyz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 76);

        auto tlx_yyz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 77);

        auto tly_yyz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 77);

        auto tlz_yyz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 77);

        auto tlx_yyz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 78);

        auto tly_yyz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 78);

        auto tlz_yyz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 78);

        auto tlx_yyz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 79);

        auto tly_yyz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 79);

        auto tlz_yyz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 79);

        auto tlx_yzz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 80);

        auto tly_yzz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 80);

        auto tlz_yzz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 80);

        auto tlx_yzz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 81);

        auto tly_yzz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 81);

        auto tlz_yzz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 81);

        auto tlx_yzz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 82);

        auto tly_yzz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 82);

        auto tlz_yzz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 82);

        auto tlx_yzz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 83);

        // Batch of Integrals (200,250)

        #pragma omp simd aligned(fgb, fx, pa_y, tdx_yy_yyy_0, tdx_yy_yyz_0, tdx_yy_yzz_0, tdx_yy_zzz_0, \
                                     tdx_yz_xxx_0, tdx_yz_xxy_0, tdx_yz_xxz_0, tdx_yz_xyy_0, tdx_yz_xyz_0, tdx_yz_xzz_0, \
                                     tdx_yz_yyy_0, tdx_yz_yyz_0, tdx_yz_yzz_0, tdx_yz_zzz_0, tdx_zz_xxx_0, tdx_zz_xxy_0, \
                                     tdx_zz_xxz_0, tdz_yy_yyz_0, tdz_yy_yzz_0, tdz_yy_zzz_0, tdz_yz_xxx_0, tdz_yz_xxy_0, \
                                     tdz_yz_xxz_0, tdz_yz_xyy_0, tdz_yz_xyz_0, tdz_yz_xzz_0, tdz_yz_yyy_0, tdz_yz_yyz_0, \
                                     tdz_yz_yzz_0, tdz_yz_zzz_0, tdz_zz_xxx_0, tdz_zz_xxy_0, tdz_zz_xxz_0, tdz_zz_xyy_0, \
                                     tlx_y_yyz_0, tlx_y_yzz_0, tlx_y_zzz_0, tlx_yy_yyz_0, tlx_yy_yz_0, tlx_yy_yzz_0, \
                                     tlx_yy_zz_0, tlx_yy_zzz_0, tlx_yyy_yyz_0, tlx_yyy_yzz_0, tlx_yyy_zzz_0, \
                                     tlx_yyz_xxx_0, tlx_yyz_xxy_0, tlx_yyz_xxz_0, tlx_yyz_xyy_0, tlx_yyz_xyz_0, \
                                     tlx_yyz_xzz_0, tlx_yyz_yyy_0, tlx_yyz_yyz_0, tlx_yyz_yzz_0, tlx_yyz_zzz_0, \
                                     tlx_yz_xx_0, tlx_yz_xxx_0, tlx_yz_xxy_0, tlx_yz_xxz_0, tlx_yz_xy_0, tlx_yz_xyy_0, \
                                     tlx_yz_xyz_0, tlx_yz_xz_0, tlx_yz_xzz_0, tlx_yz_yy_0, tlx_yz_yyy_0, tlx_yz_yyz_0, \
                                     tlx_yz_yz_0, tlx_yz_yzz_0, tlx_yz_zz_0, tlx_yz_zzz_0, tlx_yzz_xxx_0, \
                                     tlx_yzz_xxy_0, tlx_yzz_xxz_0, tlx_yzz_xyy_0, tlx_z_xxx_0, tlx_z_xxy_0, tlx_z_xxz_0, \
                                     tlx_z_xyy_0, tlx_z_xyz_0, tlx_z_xzz_0, tlx_z_yyy_0, tlx_z_yyz_0, tlx_z_yzz_0, \
                                     tlx_z_zzz_0, tlx_zz_xx_0, tlx_zz_xxx_0, tlx_zz_xxy_0, tlx_zz_xxz_0, tlx_zz_xy_0, \
                                     tlx_zz_xyy_0, tly_y_yyz_0, tly_y_yzz_0, tly_y_zzz_0, tly_yy_yyz_0, tly_yy_yz_0, \
                                     tly_yy_yzz_0, tly_yy_zz_0, tly_yy_zzz_0, tly_yyy_yyz_0, tly_yyy_yzz_0, \
                                     tly_yyy_zzz_0, tly_yyz_xxx_0, tly_yyz_xxy_0, tly_yyz_xxz_0, tly_yyz_xyy_0, \
                                     tly_yyz_xyz_0, tly_yyz_xzz_0, tly_yyz_yyy_0, tly_yyz_yyz_0, tly_yyz_yzz_0, \
                                     tly_yyz_zzz_0, tly_yz_xx_0, tly_yz_xxx_0, tly_yz_xxy_0, tly_yz_xxz_0, tly_yz_xy_0, \
                                     tly_yz_xyy_0, tly_yz_xyz_0, tly_yz_xz_0, tly_yz_xzz_0, tly_yz_yy_0, tly_yz_yyy_0, \
                                     tly_yz_yyz_0, tly_yz_yz_0, tly_yz_yzz_0, tly_yz_zz_0, tly_yz_zzz_0, tly_yzz_xxx_0, \
                                     tly_yzz_xxy_0, tly_yzz_xxz_0, tly_z_xxx_0, tly_z_xxy_0, tly_z_xxz_0, tly_z_xyy_0, \
                                     tly_z_xyz_0, tly_z_xzz_0, tly_z_yyy_0, tly_z_yyz_0, tly_z_yzz_0, tly_z_zzz_0, \
                                     tly_zz_xx_0, tly_zz_xxx_0, tly_zz_xxy_0, tly_zz_xxz_0, tlz_y_yyy_0, tlz_y_yyz_0, \
                                     tlz_y_yzz_0, tlz_y_zzz_0, tlz_yy_yy_0, tlz_yy_yyy_0, tlz_yy_yyz_0, tlz_yy_yz_0, \
                                     tlz_yy_yzz_0, tlz_yy_zz_0, tlz_yy_zzz_0, tlz_yyy_yyy_0, tlz_yyy_yyz_0, \
                                     tlz_yyy_yzz_0, tlz_yyy_zzz_0, tlz_yyz_xxx_0, tlz_yyz_xxy_0, tlz_yyz_xxz_0, \
                                     tlz_yyz_xyy_0, tlz_yyz_xyz_0, tlz_yyz_xzz_0, tlz_yyz_yyy_0, tlz_yyz_yyz_0, \
                                     tlz_yyz_yzz_0, tlz_yyz_zzz_0, tlz_yz_xx_0, tlz_yz_xxx_0, tlz_yz_xxy_0, tlz_yz_xxz_0, \
                                     tlz_yz_xy_0, tlz_yz_xyy_0, tlz_yz_xyz_0, tlz_yz_xz_0, tlz_yz_xzz_0, tlz_yz_yy_0, \
                                     tlz_yz_yyy_0, tlz_yz_yyz_0, tlz_yz_yz_0, tlz_yz_yzz_0, tlz_yz_zz_0, tlz_yz_zzz_0, \
                                     tlz_yzz_xxx_0, tlz_yzz_xxy_0, tlz_yzz_xxz_0, tlz_z_xxx_0, tlz_z_xxy_0, tlz_z_xxz_0, \
                                     tlz_z_xyy_0, tlz_z_xyz_0, tlz_z_xzz_0, tlz_z_yyy_0, tlz_z_yyz_0, tlz_z_yzz_0, \
                                     tlz_z_zzz_0, tlz_zz_xx_0, tlz_zz_xxx_0, tlz_zz_xxy_0, tlz_zz_xxz_0, tpx_yy_yyy_0, \
                                     tpx_yy_yyz_0, tpx_yy_yzz_0, tpx_yy_zzz_0, tpx_yz_xxx_0, tpx_yz_xxy_0, tpx_yz_xxz_0, \
                                     tpx_yz_xyy_0, tpx_yz_xyz_0, tpx_yz_xzz_0, tpx_yz_yyy_0, tpx_yz_yyz_0, tpx_yz_yzz_0, \
                                     tpx_yz_zzz_0, tpx_zz_xxx_0, tpx_zz_xxy_0, tpx_zz_xxz_0, tpz_yy_yyz_0, tpz_yy_yzz_0, \
                                     tpz_yy_zzz_0, tpz_yz_xxx_0, tpz_yz_xxy_0, tpz_yz_xxz_0, tpz_yz_xyy_0, tpz_yz_xyz_0, \
                                     tpz_yz_xzz_0, tpz_yz_yyy_0, tpz_yz_yyz_0, tpz_yz_yzz_0, tpz_yz_zzz_0, tpz_zz_xxx_0, \
                                     tpz_zz_xxy_0, tpz_zz_xxz_0, tpz_zz_xyy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlz_yyy_yyy_0[j] = pa_y[j] * tlz_yy_yyy_0[j] + fl1_fx * tlz_y_yyy_0[j] + 1.5 * fl1_fx * tlz_yy_yy_0[j] + 0.5 * fl1_fx * tpx_yy_yyy_0[j] +
                               fl1_fx * fl1_fgb * tdx_yy_yyy_0[j];

            tlx_yyy_yyz_0[j] = pa_y[j] * tlx_yy_yyz_0[j] + fl1_fx * tlx_y_yyz_0[j] + fl1_fx * tlx_yy_yz_0[j] - 0.5 * fl1_fx * tpz_yy_yyz_0[j] -
                               fl1_fx * fl1_fgb * tdz_yy_yyz_0[j];

            tly_yyy_yyz_0[j] = pa_y[j] * tly_yy_yyz_0[j] + fl1_fx * tly_y_yyz_0[j] + fl1_fx * tly_yy_yz_0[j];

            tlz_yyy_yyz_0[j] = pa_y[j] * tlz_yy_yyz_0[j] + fl1_fx * tlz_y_yyz_0[j] + fl1_fx * tlz_yy_yz_0[j] + 0.5 * fl1_fx * tpx_yy_yyz_0[j] +
                               fl1_fx * fl1_fgb * tdx_yy_yyz_0[j];

            tlx_yyy_yzz_0[j] = pa_y[j] * tlx_yy_yzz_0[j] + fl1_fx * tlx_y_yzz_0[j] + 0.5 * fl1_fx * tlx_yy_zz_0[j] - 0.5 * fl1_fx * tpz_yy_yzz_0[j] -
                               fl1_fx * fl1_fgb * tdz_yy_yzz_0[j];

            tly_yyy_yzz_0[j] = pa_y[j] * tly_yy_yzz_0[j] + fl1_fx * tly_y_yzz_0[j] + 0.5 * fl1_fx * tly_yy_zz_0[j];

            tlz_yyy_yzz_0[j] = pa_y[j] * tlz_yy_yzz_0[j] + fl1_fx * tlz_y_yzz_0[j] + 0.5 * fl1_fx * tlz_yy_zz_0[j] + 0.5 * fl1_fx * tpx_yy_yzz_0[j] +
                               fl1_fx * fl1_fgb * tdx_yy_yzz_0[j];

            tlx_yyy_zzz_0[j] =
                pa_y[j] * tlx_yy_zzz_0[j] + fl1_fx * tlx_y_zzz_0[j] - 0.5 * fl1_fx * tpz_yy_zzz_0[j] - fl1_fx * fl1_fgb * tdz_yy_zzz_0[j];

            tly_yyy_zzz_0[j] = pa_y[j] * tly_yy_zzz_0[j] + fl1_fx * tly_y_zzz_0[j];

            tlz_yyy_zzz_0[j] =
                pa_y[j] * tlz_yy_zzz_0[j] + fl1_fx * tlz_y_zzz_0[j] + 0.5 * fl1_fx * tpx_yy_zzz_0[j] + fl1_fx * fl1_fgb * tdx_yy_zzz_0[j];

            tlx_yyz_xxx_0[j] =
                pa_y[j] * tlx_yz_xxx_0[j] + 0.5 * fl1_fx * tlx_z_xxx_0[j] - 0.5 * fl1_fx * tpz_yz_xxx_0[j] - fl1_fx * fl1_fgb * tdz_yz_xxx_0[j];

            tly_yyz_xxx_0[j] = pa_y[j] * tly_yz_xxx_0[j] + 0.5 * fl1_fx * tly_z_xxx_0[j];

            tlz_yyz_xxx_0[j] =
                pa_y[j] * tlz_yz_xxx_0[j] + 0.5 * fl1_fx * tlz_z_xxx_0[j] + 0.5 * fl1_fx * tpx_yz_xxx_0[j] + fl1_fx * fl1_fgb * tdx_yz_xxx_0[j];

            tlx_yyz_xxy_0[j] = pa_y[j] * tlx_yz_xxy_0[j] + 0.5 * fl1_fx * tlx_z_xxy_0[j] + 0.5 * fl1_fx * tlx_yz_xx_0[j] -
                               0.5 * fl1_fx * tpz_yz_xxy_0[j] - fl1_fx * fl1_fgb * tdz_yz_xxy_0[j];

            tly_yyz_xxy_0[j] = pa_y[j] * tly_yz_xxy_0[j] + 0.5 * fl1_fx * tly_z_xxy_0[j] + 0.5 * fl1_fx * tly_yz_xx_0[j];

            tlz_yyz_xxy_0[j] = pa_y[j] * tlz_yz_xxy_0[j] + 0.5 * fl1_fx * tlz_z_xxy_0[j] + 0.5 * fl1_fx * tlz_yz_xx_0[j] +
                               0.5 * fl1_fx * tpx_yz_xxy_0[j] + fl1_fx * fl1_fgb * tdx_yz_xxy_0[j];

            tlx_yyz_xxz_0[j] =
                pa_y[j] * tlx_yz_xxz_0[j] + 0.5 * fl1_fx * tlx_z_xxz_0[j] - 0.5 * fl1_fx * tpz_yz_xxz_0[j] - fl1_fx * fl1_fgb * tdz_yz_xxz_0[j];

            tly_yyz_xxz_0[j] = pa_y[j] * tly_yz_xxz_0[j] + 0.5 * fl1_fx * tly_z_xxz_0[j];

            tlz_yyz_xxz_0[j] =
                pa_y[j] * tlz_yz_xxz_0[j] + 0.5 * fl1_fx * tlz_z_xxz_0[j] + 0.5 * fl1_fx * tpx_yz_xxz_0[j] + fl1_fx * fl1_fgb * tdx_yz_xxz_0[j];

            tlx_yyz_xyy_0[j] = pa_y[j] * tlx_yz_xyy_0[j] + 0.5 * fl1_fx * tlx_z_xyy_0[j] + fl1_fx * tlx_yz_xy_0[j] - 0.5 * fl1_fx * tpz_yz_xyy_0[j] -
                               fl1_fx * fl1_fgb * tdz_yz_xyy_0[j];

            tly_yyz_xyy_0[j] = pa_y[j] * tly_yz_xyy_0[j] + 0.5 * fl1_fx * tly_z_xyy_0[j] + fl1_fx * tly_yz_xy_0[j];

            tlz_yyz_xyy_0[j] = pa_y[j] * tlz_yz_xyy_0[j] + 0.5 * fl1_fx * tlz_z_xyy_0[j] + fl1_fx * tlz_yz_xy_0[j] + 0.5 * fl1_fx * tpx_yz_xyy_0[j] +
                               fl1_fx * fl1_fgb * tdx_yz_xyy_0[j];

            tlx_yyz_xyz_0[j] = pa_y[j] * tlx_yz_xyz_0[j] + 0.5 * fl1_fx * tlx_z_xyz_0[j] + 0.5 * fl1_fx * tlx_yz_xz_0[j] -
                               0.5 * fl1_fx * tpz_yz_xyz_0[j] - fl1_fx * fl1_fgb * tdz_yz_xyz_0[j];

            tly_yyz_xyz_0[j] = pa_y[j] * tly_yz_xyz_0[j] + 0.5 * fl1_fx * tly_z_xyz_0[j] + 0.5 * fl1_fx * tly_yz_xz_0[j];

            tlz_yyz_xyz_0[j] = pa_y[j] * tlz_yz_xyz_0[j] + 0.5 * fl1_fx * tlz_z_xyz_0[j] + 0.5 * fl1_fx * tlz_yz_xz_0[j] +
                               0.5 * fl1_fx * tpx_yz_xyz_0[j] + fl1_fx * fl1_fgb * tdx_yz_xyz_0[j];

            tlx_yyz_xzz_0[j] =
                pa_y[j] * tlx_yz_xzz_0[j] + 0.5 * fl1_fx * tlx_z_xzz_0[j] - 0.5 * fl1_fx * tpz_yz_xzz_0[j] - fl1_fx * fl1_fgb * tdz_yz_xzz_0[j];

            tly_yyz_xzz_0[j] = pa_y[j] * tly_yz_xzz_0[j] + 0.5 * fl1_fx * tly_z_xzz_0[j];

            tlz_yyz_xzz_0[j] =
                pa_y[j] * tlz_yz_xzz_0[j] + 0.5 * fl1_fx * tlz_z_xzz_0[j] + 0.5 * fl1_fx * tpx_yz_xzz_0[j] + fl1_fx * fl1_fgb * tdx_yz_xzz_0[j];

            tlx_yyz_yyy_0[j] = pa_y[j] * tlx_yz_yyy_0[j] + 0.5 * fl1_fx * tlx_z_yyy_0[j] + 1.5 * fl1_fx * tlx_yz_yy_0[j] -
                               0.5 * fl1_fx * tpz_yz_yyy_0[j] - fl1_fx * fl1_fgb * tdz_yz_yyy_0[j];

            tly_yyz_yyy_0[j] = pa_y[j] * tly_yz_yyy_0[j] + 0.5 * fl1_fx * tly_z_yyy_0[j] + 1.5 * fl1_fx * tly_yz_yy_0[j];

            tlz_yyz_yyy_0[j] = pa_y[j] * tlz_yz_yyy_0[j] + 0.5 * fl1_fx * tlz_z_yyy_0[j] + 1.5 * fl1_fx * tlz_yz_yy_0[j] +
                               0.5 * fl1_fx * tpx_yz_yyy_0[j] + fl1_fx * fl1_fgb * tdx_yz_yyy_0[j];

            tlx_yyz_yyz_0[j] = pa_y[j] * tlx_yz_yyz_0[j] + 0.5 * fl1_fx * tlx_z_yyz_0[j] + fl1_fx * tlx_yz_yz_0[j] - 0.5 * fl1_fx * tpz_yz_yyz_0[j] -
                               fl1_fx * fl1_fgb * tdz_yz_yyz_0[j];

            tly_yyz_yyz_0[j] = pa_y[j] * tly_yz_yyz_0[j] + 0.5 * fl1_fx * tly_z_yyz_0[j] + fl1_fx * tly_yz_yz_0[j];

            tlz_yyz_yyz_0[j] = pa_y[j] * tlz_yz_yyz_0[j] + 0.5 * fl1_fx * tlz_z_yyz_0[j] + fl1_fx * tlz_yz_yz_0[j] + 0.5 * fl1_fx * tpx_yz_yyz_0[j] +
                               fl1_fx * fl1_fgb * tdx_yz_yyz_0[j];

            tlx_yyz_yzz_0[j] = pa_y[j] * tlx_yz_yzz_0[j] + 0.5 * fl1_fx * tlx_z_yzz_0[j] + 0.5 * fl1_fx * tlx_yz_zz_0[j] -
                               0.5 * fl1_fx * tpz_yz_yzz_0[j] - fl1_fx * fl1_fgb * tdz_yz_yzz_0[j];

            tly_yyz_yzz_0[j] = pa_y[j] * tly_yz_yzz_0[j] + 0.5 * fl1_fx * tly_z_yzz_0[j] + 0.5 * fl1_fx * tly_yz_zz_0[j];

            tlz_yyz_yzz_0[j] = pa_y[j] * tlz_yz_yzz_0[j] + 0.5 * fl1_fx * tlz_z_yzz_0[j] + 0.5 * fl1_fx * tlz_yz_zz_0[j] +
                               0.5 * fl1_fx * tpx_yz_yzz_0[j] + fl1_fx * fl1_fgb * tdx_yz_yzz_0[j];

            tlx_yyz_zzz_0[j] =
                pa_y[j] * tlx_yz_zzz_0[j] + 0.5 * fl1_fx * tlx_z_zzz_0[j] - 0.5 * fl1_fx * tpz_yz_zzz_0[j] - fl1_fx * fl1_fgb * tdz_yz_zzz_0[j];

            tly_yyz_zzz_0[j] = pa_y[j] * tly_yz_zzz_0[j] + 0.5 * fl1_fx * tly_z_zzz_0[j];

            tlz_yyz_zzz_0[j] =
                pa_y[j] * tlz_yz_zzz_0[j] + 0.5 * fl1_fx * tlz_z_zzz_0[j] + 0.5 * fl1_fx * tpx_yz_zzz_0[j] + fl1_fx * fl1_fgb * tdx_yz_zzz_0[j];

            tlx_yzz_xxx_0[j] = pa_y[j] * tlx_zz_xxx_0[j] - 0.5 * fl1_fx * tpz_zz_xxx_0[j] - fl1_fx * fl1_fgb * tdz_zz_xxx_0[j];

            tly_yzz_xxx_0[j] = pa_y[j] * tly_zz_xxx_0[j];

            tlz_yzz_xxx_0[j] = pa_y[j] * tlz_zz_xxx_0[j] + 0.5 * fl1_fx * tpx_zz_xxx_0[j] + fl1_fx * fl1_fgb * tdx_zz_xxx_0[j];

            tlx_yzz_xxy_0[j] =
                pa_y[j] * tlx_zz_xxy_0[j] + 0.5 * fl1_fx * tlx_zz_xx_0[j] - 0.5 * fl1_fx * tpz_zz_xxy_0[j] - fl1_fx * fl1_fgb * tdz_zz_xxy_0[j];

            tly_yzz_xxy_0[j] = pa_y[j] * tly_zz_xxy_0[j] + 0.5 * fl1_fx * tly_zz_xx_0[j];

            tlz_yzz_xxy_0[j] =
                pa_y[j] * tlz_zz_xxy_0[j] + 0.5 * fl1_fx * tlz_zz_xx_0[j] + 0.5 * fl1_fx * tpx_zz_xxy_0[j] + fl1_fx * fl1_fgb * tdx_zz_xxy_0[j];

            tlx_yzz_xxz_0[j] = pa_y[j] * tlx_zz_xxz_0[j] - 0.5 * fl1_fx * tpz_zz_xxz_0[j] - fl1_fx * fl1_fgb * tdz_zz_xxz_0[j];

            tly_yzz_xxz_0[j] = pa_y[j] * tly_zz_xxz_0[j];

            tlz_yzz_xxz_0[j] = pa_y[j] * tlz_zz_xxz_0[j] + 0.5 * fl1_fx * tpx_zz_xxz_0[j] + fl1_fx * fl1_fgb * tdx_zz_xxz_0[j];

            tlx_yzz_xyy_0[j] =
                pa_y[j] * tlx_zz_xyy_0[j] + fl1_fx * tlx_zz_xy_0[j] - 0.5 * fl1_fx * tpz_zz_xyy_0[j] - fl1_fx * fl1_fgb * tdz_zz_xyy_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForFF_250_300(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_l_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_3_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tlx_zz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 30);

        auto tly_zz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 30);

        auto tlz_zz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 30);

        auto tlx_zz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 31);

        auto tly_zz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 31);

        auto tlz_zz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 31);

        auto tlx_zz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 32);

        auto tly_zz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 32);

        auto tlz_zz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 32);

        auto tlx_zz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 33);

        auto tly_zz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 33);

        auto tlz_zz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 33);

        auto tlx_zz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 34);

        auto tly_zz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 34);

        auto tlz_zz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 34);

        auto tlx_zz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 35);

        auto tly_zz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 35);

        auto tlz_zz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 35);

        auto tpx_zz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 50);

        auto tpy_zz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 50);

        auto tpx_zz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 51);

        auto tpy_zz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 51);

        auto tpx_zz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 52);

        auto tpy_zz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 52);

        auto tpx_zz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 53);

        auto tpy_zz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 53);

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

        auto tdx_zz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 50);

        auto tdy_zz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 50);

        auto tdx_zz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 51);

        auto tdy_zz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 51);

        auto tdx_zz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 52);

        auto tdy_zz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 52);

        auto tdx_zz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 53);

        auto tdy_zz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 53);

        auto tdx_zz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 54);

        auto tdy_zz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 54);

        auto tdz_zz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 54);

        auto tdx_zz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 55);

        auto tdy_zz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 55);

        auto tdz_zz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 55);

        auto tdx_zz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 56);

        auto tdy_zz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 56);

        auto tdz_zz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 56);

        auto tdx_zz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 57);

        auto tdy_zz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 57);

        auto tdz_zz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 57);

        auto tdx_zz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 58);

        auto tdy_zz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 58);

        auto tdz_zz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 58);

        auto tdx_zz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 59);

        auto tdy_zz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 59);

        auto tdz_zz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 59);

        // set up pointers to integrals

        auto tly_yzz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 83);

        auto tlz_yzz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 83);

        auto tlx_yzz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 84);

        auto tly_yzz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 84);

        auto tlz_yzz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 84);

        auto tlx_yzz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 85);

        auto tly_yzz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 85);

        auto tlz_yzz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 85);

        auto tlx_yzz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 86);

        auto tly_yzz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 86);

        auto tlz_yzz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 86);

        auto tlx_yzz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 87);

        auto tly_yzz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 87);

        auto tlz_yzz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 87);

        auto tlx_yzz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 88);

        auto tly_yzz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 88);

        auto tlz_yzz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 88);

        auto tlx_yzz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 89);

        auto tly_yzz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 89);

        auto tlz_yzz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 89);

        auto tlx_zzz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 90);

        auto tly_zzz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 90);

        auto tlz_zzz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 90);

        auto tlx_zzz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 91);

        auto tly_zzz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 91);

        auto tlz_zzz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 91);

        auto tlx_zzz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 92);

        auto tly_zzz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 92);

        auto tlz_zzz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 92);

        auto tlx_zzz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 93);

        auto tly_zzz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 93);

        auto tlz_zzz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 93);

        auto tlx_zzz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 94);

        auto tly_zzz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 94);

        auto tlz_zzz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 94);

        auto tlx_zzz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 95);

        auto tly_zzz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 95);

        auto tlz_zzz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 95);

        auto tlx_zzz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 96);

        auto tly_zzz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 96);

        auto tlz_zzz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 96);

        auto tlx_zzz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 97);

        auto tly_zzz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 97);

        auto tlz_zzz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 97);

        auto tlx_zzz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 98);

        auto tly_zzz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 98);

        auto tlz_zzz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 98);

        auto tlx_zzz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 99);

        auto tly_zzz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 99);

        auto tlz_zzz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 99);

        // Batch of Integrals (250,300)

        #pragma omp simd aligned(fgb, fx, pa_y, pa_z, tdx_zz_xxx_0, tdx_zz_xxy_0, tdx_zz_xxz_0, \
                                     tdx_zz_xyy_0, tdx_zz_xyz_0, tdx_zz_xzz_0, tdx_zz_yyy_0, tdx_zz_yyz_0, tdx_zz_yzz_0, \
                                     tdx_zz_zzz_0, tdy_zz_xxx_0, tdy_zz_xxy_0, tdy_zz_xxz_0, tdy_zz_xyy_0, tdy_zz_xyz_0, \
                                     tdy_zz_xzz_0, tdy_zz_yyy_0, tdy_zz_yyz_0, tdy_zz_yzz_0, tdy_zz_zzz_0, tdz_zz_xyz_0, \
                                     tdz_zz_xzz_0, tdz_zz_yyy_0, tdz_zz_yyz_0, tdz_zz_yzz_0, tdz_zz_zzz_0, \
                                     tlx_yzz_xyz_0, tlx_yzz_xzz_0, tlx_yzz_yyy_0, tlx_yzz_yyz_0, tlx_yzz_yzz_0, \
                                     tlx_yzz_zzz_0, tlx_z_xxx_0, tlx_z_xxy_0, tlx_z_xxz_0, tlx_z_xyy_0, tlx_z_xyz_0, \
                                     tlx_z_xzz_0, tlx_z_yyy_0, tlx_z_yyz_0, tlx_z_yzz_0, tlx_z_zzz_0, tlx_zz_xx_0, \
                                     tlx_zz_xxx_0, tlx_zz_xxy_0, tlx_zz_xxz_0, tlx_zz_xy_0, tlx_zz_xyy_0, tlx_zz_xyz_0, \
                                     tlx_zz_xz_0, tlx_zz_xzz_0, tlx_zz_yy_0, tlx_zz_yyy_0, tlx_zz_yyz_0, tlx_zz_yz_0, \
                                     tlx_zz_yzz_0, tlx_zz_zz_0, tlx_zz_zzz_0, tlx_zzz_xxx_0, tlx_zzz_xxy_0, \
                                     tlx_zzz_xxz_0, tlx_zzz_xyy_0, tlx_zzz_xyz_0, tlx_zzz_xzz_0, tlx_zzz_yyy_0, \
                                     tlx_zzz_yyz_0, tlx_zzz_yzz_0, tlx_zzz_zzz_0, tly_yzz_xyy_0, tly_yzz_xyz_0, \
                                     tly_yzz_xzz_0, tly_yzz_yyy_0, tly_yzz_yyz_0, tly_yzz_yzz_0, tly_yzz_zzz_0, \
                                     tly_z_xxx_0, tly_z_xxy_0, tly_z_xxz_0, tly_z_xyy_0, tly_z_xyz_0, tly_z_xzz_0, \
                                     tly_z_yyy_0, tly_z_yyz_0, tly_z_yzz_0, tly_z_zzz_0, tly_zz_xx_0, tly_zz_xxx_0, \
                                     tly_zz_xxy_0, tly_zz_xxz_0, tly_zz_xy_0, tly_zz_xyy_0, tly_zz_xyz_0, tly_zz_xz_0, \
                                     tly_zz_xzz_0, tly_zz_yy_0, tly_zz_yyy_0, tly_zz_yyz_0, tly_zz_yz_0, tly_zz_yzz_0, \
                                     tly_zz_zz_0, tly_zz_zzz_0, tly_zzz_xxx_0, tly_zzz_xxy_0, tly_zzz_xxz_0, \
                                     tly_zzz_xyy_0, tly_zzz_xyz_0, tly_zzz_xzz_0, tly_zzz_yyy_0, tly_zzz_yyz_0, \
                                     tly_zzz_yzz_0, tly_zzz_zzz_0, tlz_yzz_xyy_0, tlz_yzz_xyz_0, tlz_yzz_xzz_0, \
                                     tlz_yzz_yyy_0, tlz_yzz_yyz_0, tlz_yzz_yzz_0, tlz_yzz_zzz_0, tlz_z_xxx_0, \
                                     tlz_z_xxy_0, tlz_z_xxz_0, tlz_z_xyy_0, tlz_z_xyz_0, tlz_z_xzz_0, tlz_z_yyy_0, \
                                     tlz_z_yyz_0, tlz_z_yzz_0, tlz_z_zzz_0, tlz_zz_xx_0, tlz_zz_xxx_0, tlz_zz_xxy_0, \
                                     tlz_zz_xxz_0, tlz_zz_xy_0, tlz_zz_xyy_0, tlz_zz_xyz_0, tlz_zz_xz_0, tlz_zz_xzz_0, \
                                     tlz_zz_yy_0, tlz_zz_yyy_0, tlz_zz_yyz_0, tlz_zz_yz_0, tlz_zz_yzz_0, tlz_zz_zz_0, \
                                     tlz_zz_zzz_0, tlz_zzz_xxx_0, tlz_zzz_xxy_0, tlz_zzz_xxz_0, tlz_zzz_xyy_0, \
                                     tlz_zzz_xyz_0, tlz_zzz_xzz_0, tlz_zzz_yyy_0, tlz_zzz_yyz_0, tlz_zzz_yzz_0, \
                                     tlz_zzz_zzz_0, tpx_zz_xxx_0, tpx_zz_xxy_0, tpx_zz_xxz_0, tpx_zz_xyy_0, tpx_zz_xyz_0, \
                                     tpx_zz_xzz_0, tpx_zz_yyy_0, tpx_zz_yyz_0, tpx_zz_yzz_0, tpx_zz_zzz_0, tpy_zz_xxx_0, \
                                     tpy_zz_xxy_0, tpy_zz_xxz_0, tpy_zz_xyy_0, tpy_zz_xyz_0, tpy_zz_xzz_0, tpy_zz_yyy_0, \
                                     tpy_zz_yyz_0, tpy_zz_yzz_0, tpy_zz_zzz_0, tpz_zz_xyz_0, tpz_zz_xzz_0, tpz_zz_yyy_0, \
                                     tpz_zz_yyz_0, tpz_zz_yzz_0, tpz_zz_zzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tly_yzz_xyy_0[j] = pa_y[j] * tly_zz_xyy_0[j] + fl1_fx * tly_zz_xy_0[j];

            tlz_yzz_xyy_0[j] =
                pa_y[j] * tlz_zz_xyy_0[j] + fl1_fx * tlz_zz_xy_0[j] + 0.5 * fl1_fx * tpx_zz_xyy_0[j] + fl1_fx * fl1_fgb * tdx_zz_xyy_0[j];

            tlx_yzz_xyz_0[j] =
                pa_y[j] * tlx_zz_xyz_0[j] + 0.5 * fl1_fx * tlx_zz_xz_0[j] - 0.5 * fl1_fx * tpz_zz_xyz_0[j] - fl1_fx * fl1_fgb * tdz_zz_xyz_0[j];

            tly_yzz_xyz_0[j] = pa_y[j] * tly_zz_xyz_0[j] + 0.5 * fl1_fx * tly_zz_xz_0[j];

            tlz_yzz_xyz_0[j] =
                pa_y[j] * tlz_zz_xyz_0[j] + 0.5 * fl1_fx * tlz_zz_xz_0[j] + 0.5 * fl1_fx * tpx_zz_xyz_0[j] + fl1_fx * fl1_fgb * tdx_zz_xyz_0[j];

            tlx_yzz_xzz_0[j] = pa_y[j] * tlx_zz_xzz_0[j] - 0.5 * fl1_fx * tpz_zz_xzz_0[j] - fl1_fx * fl1_fgb * tdz_zz_xzz_0[j];

            tly_yzz_xzz_0[j] = pa_y[j] * tly_zz_xzz_0[j];

            tlz_yzz_xzz_0[j] = pa_y[j] * tlz_zz_xzz_0[j] + 0.5 * fl1_fx * tpx_zz_xzz_0[j] + fl1_fx * fl1_fgb * tdx_zz_xzz_0[j];

            tlx_yzz_yyy_0[j] =
                pa_y[j] * tlx_zz_yyy_0[j] + 1.5 * fl1_fx * tlx_zz_yy_0[j] - 0.5 * fl1_fx * tpz_zz_yyy_0[j] - fl1_fx * fl1_fgb * tdz_zz_yyy_0[j];

            tly_yzz_yyy_0[j] = pa_y[j] * tly_zz_yyy_0[j] + 1.5 * fl1_fx * tly_zz_yy_0[j];

            tlz_yzz_yyy_0[j] =
                pa_y[j] * tlz_zz_yyy_0[j] + 1.5 * fl1_fx * tlz_zz_yy_0[j] + 0.5 * fl1_fx * tpx_zz_yyy_0[j] + fl1_fx * fl1_fgb * tdx_zz_yyy_0[j];

            tlx_yzz_yyz_0[j] =
                pa_y[j] * tlx_zz_yyz_0[j] + fl1_fx * tlx_zz_yz_0[j] - 0.5 * fl1_fx * tpz_zz_yyz_0[j] - fl1_fx * fl1_fgb * tdz_zz_yyz_0[j];

            tly_yzz_yyz_0[j] = pa_y[j] * tly_zz_yyz_0[j] + fl1_fx * tly_zz_yz_0[j];

            tlz_yzz_yyz_0[j] =
                pa_y[j] * tlz_zz_yyz_0[j] + fl1_fx * tlz_zz_yz_0[j] + 0.5 * fl1_fx * tpx_zz_yyz_0[j] + fl1_fx * fl1_fgb * tdx_zz_yyz_0[j];

            tlx_yzz_yzz_0[j] =
                pa_y[j] * tlx_zz_yzz_0[j] + 0.5 * fl1_fx * tlx_zz_zz_0[j] - 0.5 * fl1_fx * tpz_zz_yzz_0[j] - fl1_fx * fl1_fgb * tdz_zz_yzz_0[j];

            tly_yzz_yzz_0[j] = pa_y[j] * tly_zz_yzz_0[j] + 0.5 * fl1_fx * tly_zz_zz_0[j];

            tlz_yzz_yzz_0[j] =
                pa_y[j] * tlz_zz_yzz_0[j] + 0.5 * fl1_fx * tlz_zz_zz_0[j] + 0.5 * fl1_fx * tpx_zz_yzz_0[j] + fl1_fx * fl1_fgb * tdx_zz_yzz_0[j];

            tlx_yzz_zzz_0[j] = pa_y[j] * tlx_zz_zzz_0[j] - 0.5 * fl1_fx * tpz_zz_zzz_0[j] - fl1_fx * fl1_fgb * tdz_zz_zzz_0[j];

            tly_yzz_zzz_0[j] = pa_y[j] * tly_zz_zzz_0[j];

            tlz_yzz_zzz_0[j] = pa_y[j] * tlz_zz_zzz_0[j] + 0.5 * fl1_fx * tpx_zz_zzz_0[j] + fl1_fx * fl1_fgb * tdx_zz_zzz_0[j];

            tlx_zzz_xxx_0[j] =
                pa_z[j] * tlx_zz_xxx_0[j] + fl1_fx * tlx_z_xxx_0[j] + 0.5 * fl1_fx * tpy_zz_xxx_0[j] + fl1_fx * fl1_fgb * tdy_zz_xxx_0[j];

            tly_zzz_xxx_0[j] =
                pa_z[j] * tly_zz_xxx_0[j] + fl1_fx * tly_z_xxx_0[j] - 0.5 * fl1_fx * tpx_zz_xxx_0[j] - fl1_fx * fl1_fgb * tdx_zz_xxx_0[j];

            tlz_zzz_xxx_0[j] = pa_z[j] * tlz_zz_xxx_0[j] + fl1_fx * tlz_z_xxx_0[j];

            tlx_zzz_xxy_0[j] =
                pa_z[j] * tlx_zz_xxy_0[j] + fl1_fx * tlx_z_xxy_0[j] + 0.5 * fl1_fx * tpy_zz_xxy_0[j] + fl1_fx * fl1_fgb * tdy_zz_xxy_0[j];

            tly_zzz_xxy_0[j] =
                pa_z[j] * tly_zz_xxy_0[j] + fl1_fx * tly_z_xxy_0[j] - 0.5 * fl1_fx * tpx_zz_xxy_0[j] - fl1_fx * fl1_fgb * tdx_zz_xxy_0[j];

            tlz_zzz_xxy_0[j] = pa_z[j] * tlz_zz_xxy_0[j] + fl1_fx * tlz_z_xxy_0[j];

            tlx_zzz_xxz_0[j] = pa_z[j] * tlx_zz_xxz_0[j] + fl1_fx * tlx_z_xxz_0[j] + 0.5 * fl1_fx * tlx_zz_xx_0[j] + 0.5 * fl1_fx * tpy_zz_xxz_0[j] +
                               fl1_fx * fl1_fgb * tdy_zz_xxz_0[j];

            tly_zzz_xxz_0[j] = pa_z[j] * tly_zz_xxz_0[j] + fl1_fx * tly_z_xxz_0[j] + 0.5 * fl1_fx * tly_zz_xx_0[j] - 0.5 * fl1_fx * tpx_zz_xxz_0[j] -
                               fl1_fx * fl1_fgb * tdx_zz_xxz_0[j];

            tlz_zzz_xxz_0[j] = pa_z[j] * tlz_zz_xxz_0[j] + fl1_fx * tlz_z_xxz_0[j] + 0.5 * fl1_fx * tlz_zz_xx_0[j];

            tlx_zzz_xyy_0[j] =
                pa_z[j] * tlx_zz_xyy_0[j] + fl1_fx * tlx_z_xyy_0[j] + 0.5 * fl1_fx * tpy_zz_xyy_0[j] + fl1_fx * fl1_fgb * tdy_zz_xyy_0[j];

            tly_zzz_xyy_0[j] =
                pa_z[j] * tly_zz_xyy_0[j] + fl1_fx * tly_z_xyy_0[j] - 0.5 * fl1_fx * tpx_zz_xyy_0[j] - fl1_fx * fl1_fgb * tdx_zz_xyy_0[j];

            tlz_zzz_xyy_0[j] = pa_z[j] * tlz_zz_xyy_0[j] + fl1_fx * tlz_z_xyy_0[j];

            tlx_zzz_xyz_0[j] = pa_z[j] * tlx_zz_xyz_0[j] + fl1_fx * tlx_z_xyz_0[j] + 0.5 * fl1_fx * tlx_zz_xy_0[j] + 0.5 * fl1_fx * tpy_zz_xyz_0[j] +
                               fl1_fx * fl1_fgb * tdy_zz_xyz_0[j];

            tly_zzz_xyz_0[j] = pa_z[j] * tly_zz_xyz_0[j] + fl1_fx * tly_z_xyz_0[j] + 0.5 * fl1_fx * tly_zz_xy_0[j] - 0.5 * fl1_fx * tpx_zz_xyz_0[j] -
                               fl1_fx * fl1_fgb * tdx_zz_xyz_0[j];

            tlz_zzz_xyz_0[j] = pa_z[j] * tlz_zz_xyz_0[j] + fl1_fx * tlz_z_xyz_0[j] + 0.5 * fl1_fx * tlz_zz_xy_0[j];

            tlx_zzz_xzz_0[j] = pa_z[j] * tlx_zz_xzz_0[j] + fl1_fx * tlx_z_xzz_0[j] + fl1_fx * tlx_zz_xz_0[j] + 0.5 * fl1_fx * tpy_zz_xzz_0[j] +
                               fl1_fx * fl1_fgb * tdy_zz_xzz_0[j];

            tly_zzz_xzz_0[j] = pa_z[j] * tly_zz_xzz_0[j] + fl1_fx * tly_z_xzz_0[j] + fl1_fx * tly_zz_xz_0[j] - 0.5 * fl1_fx * tpx_zz_xzz_0[j] -
                               fl1_fx * fl1_fgb * tdx_zz_xzz_0[j];

            tlz_zzz_xzz_0[j] = pa_z[j] * tlz_zz_xzz_0[j] + fl1_fx * tlz_z_xzz_0[j] + fl1_fx * tlz_zz_xz_0[j];

            tlx_zzz_yyy_0[j] =
                pa_z[j] * tlx_zz_yyy_0[j] + fl1_fx * tlx_z_yyy_0[j] + 0.5 * fl1_fx * tpy_zz_yyy_0[j] + fl1_fx * fl1_fgb * tdy_zz_yyy_0[j];

            tly_zzz_yyy_0[j] =
                pa_z[j] * tly_zz_yyy_0[j] + fl1_fx * tly_z_yyy_0[j] - 0.5 * fl1_fx * tpx_zz_yyy_0[j] - fl1_fx * fl1_fgb * tdx_zz_yyy_0[j];

            tlz_zzz_yyy_0[j] = pa_z[j] * tlz_zz_yyy_0[j] + fl1_fx * tlz_z_yyy_0[j];

            tlx_zzz_yyz_0[j] = pa_z[j] * tlx_zz_yyz_0[j] + fl1_fx * tlx_z_yyz_0[j] + 0.5 * fl1_fx * tlx_zz_yy_0[j] + 0.5 * fl1_fx * tpy_zz_yyz_0[j] +
                               fl1_fx * fl1_fgb * tdy_zz_yyz_0[j];

            tly_zzz_yyz_0[j] = pa_z[j] * tly_zz_yyz_0[j] + fl1_fx * tly_z_yyz_0[j] + 0.5 * fl1_fx * tly_zz_yy_0[j] - 0.5 * fl1_fx * tpx_zz_yyz_0[j] -
                               fl1_fx * fl1_fgb * tdx_zz_yyz_0[j];

            tlz_zzz_yyz_0[j] = pa_z[j] * tlz_zz_yyz_0[j] + fl1_fx * tlz_z_yyz_0[j] + 0.5 * fl1_fx * tlz_zz_yy_0[j];

            tlx_zzz_yzz_0[j] = pa_z[j] * tlx_zz_yzz_0[j] + fl1_fx * tlx_z_yzz_0[j] + fl1_fx * tlx_zz_yz_0[j] + 0.5 * fl1_fx * tpy_zz_yzz_0[j] +
                               fl1_fx * fl1_fgb * tdy_zz_yzz_0[j];

            tly_zzz_yzz_0[j] = pa_z[j] * tly_zz_yzz_0[j] + fl1_fx * tly_z_yzz_0[j] + fl1_fx * tly_zz_yz_0[j] - 0.5 * fl1_fx * tpx_zz_yzz_0[j] -
                               fl1_fx * fl1_fgb * tdx_zz_yzz_0[j];

            tlz_zzz_yzz_0[j] = pa_z[j] * tlz_zz_yzz_0[j] + fl1_fx * tlz_z_yzz_0[j] + fl1_fx * tlz_zz_yz_0[j];

            tlx_zzz_zzz_0[j] = pa_z[j] * tlx_zz_zzz_0[j] + fl1_fx * tlx_z_zzz_0[j] + 1.5 * fl1_fx * tlx_zz_zz_0[j] + 0.5 * fl1_fx * tpy_zz_zzz_0[j] +
                               fl1_fx * fl1_fgb * tdy_zz_zzz_0[j];

            tly_zzz_zzz_0[j] = pa_z[j] * tly_zz_zzz_0[j] + fl1_fx * tly_z_zzz_0[j] + 1.5 * fl1_fx * tly_zz_zz_0[j] - 0.5 * fl1_fx * tpx_zz_zzz_0[j] -
                               fl1_fx * fl1_fgb * tdx_zz_zzz_0[j];

            tlz_zzz_zzz_0[j] = pa_z[j] * tlz_zz_zzz_0[j] + fl1_fx * tlz_z_zzz_0[j] + 1.5 * fl1_fx * tlz_zz_zz_0[j];
        }

        idx++;
    }
}

}  // namespace amomrecfunc
