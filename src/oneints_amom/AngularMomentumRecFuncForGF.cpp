//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "AngularMomentumRecFuncForGF.hpp"

namespace amomrecfunc {  // amomrecfunc namespace

void
compAngularMomentumForGF(CMemBlock2D<double>&       primBuffer,
                         const CRecursionMap&       recursionMap,
                         const CMemBlock2D<double>& osFactors,
                         const CMemBlock2D<double>& paDistances,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const int32_t              iContrGto)
{
    amomrecfunc::compAngularMomentumForGF_0_50(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForGF_50_100(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForGF_100_150(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForGF_150_200(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForGF_200_250(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForGF_250_300(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForGF_300_350(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForGF_350_400(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForGF_400_450(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compAngularMomentumForGF_0_50(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_l_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_4_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tlx_xxx_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx);

        auto tly_xxx_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx);

        auto tlz_xxx_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx);

        auto tlx_xxx_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 1);

        auto tly_xxx_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 1);

        auto tlz_xxx_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 1);

        auto tlx_xxx_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 2);

        auto tly_xxx_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 2);

        auto tlz_xxx_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 2);

        auto tlx_xxx_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 3);

        auto tly_xxx_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 3);

        auto tlz_xxx_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 3);

        auto tlx_xxx_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 4);

        auto tly_xxx_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 4);

        auto tlz_xxx_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 4);

        auto tlx_xxx_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 5);

        auto tly_xxx_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 5);

        auto tlz_xxx_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 5);

        auto tlx_xxy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 6);

        auto tly_xxy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 6);

        auto tlz_xxy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 6);

        auto tlx_xxy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 7);

        auto tly_xxy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 7);

        auto tlz_xxy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 7);

        auto tlx_xxy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 8);

        auto tly_xxy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 8);

        auto tlz_xxy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 8);

        auto tlx_xxy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 9);

        auto tly_xxy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 9);

        auto tlz_xxy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 9);

        auto tlx_xxy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 10);

        auto tly_xxy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 10);

        auto tlz_xxy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 10);

        auto tlx_xxy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 11);

        auto tly_xxy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 11);

        auto tlz_xxy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 11);

        auto tpy_xxx_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx);

        auto tpz_xxx_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx);

        auto tpy_xxx_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 1);

        auto tpz_xxx_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 1);

        auto tpy_xxx_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 2);

        auto tpz_xxx_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 2);

        auto tpy_xxx_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 3);

        auto tpz_xxx_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 3);

        auto tpy_xxx_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 4);

        auto tpz_xxx_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 4);

        auto tpy_xxx_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 5);

        auto tpz_xxx_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 5);

        auto tpy_xxx_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 6);

        auto tpz_xxx_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 6);

        auto tpy_xxx_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 7);

        auto tpz_xxx_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 7);

        auto tpy_xxx_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 8);

        auto tpz_xxx_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 8);

        auto tpy_xxx_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 9);

        auto tpz_xxx_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 9);

        auto tpy_xxy_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 10);

        auto tpz_xxy_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 10);

        auto tpy_xxy_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 11);

        auto tpz_xxy_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 11);

        auto tpy_xxy_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 12);

        auto tpz_xxy_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 12);

        auto tpy_xxy_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 13);

        auto tpz_xxy_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 13);

        auto tpy_xxy_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 14);

        auto tpz_xxy_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 14);

        auto tpy_xxy_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 15);

        auto tpz_xxy_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 15);

        auto tpz_xxy_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 16);

        auto tdy_xxx_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx);

        auto tdz_xxx_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx);

        auto tdy_xxx_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 1);

        auto tdz_xxx_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 1);

        auto tdy_xxx_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 2);

        auto tdz_xxx_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 2);

        auto tdy_xxx_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 3);

        auto tdz_xxx_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 3);

        auto tdy_xxx_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 4);

        auto tdz_xxx_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 4);

        auto tdy_xxx_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 5);

        auto tdz_xxx_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 5);

        auto tdy_xxx_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 6);

        auto tdz_xxx_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 6);

        auto tdy_xxx_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 7);

        auto tdz_xxx_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 7);

        auto tdy_xxx_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 8);

        auto tdz_xxx_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 8);

        auto tdy_xxx_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 9);

        auto tdz_xxx_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 9);

        auto tdy_xxy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 10);

        auto tdz_xxy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 10);

        auto tdy_xxy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 11);

        auto tdz_xxy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 11);

        auto tdy_xxy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 12);

        auto tdz_xxy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 12);

        auto tdy_xxy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 13);

        auto tdz_xxy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 13);

        auto tdy_xxy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 14);

        auto tdz_xxy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 14);

        auto tdy_xxy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 15);

        auto tdz_xxy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 15);

        auto tdz_xxy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 16);

        // set up pointers to integrals

        auto tlx_xxxx_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx);

        auto tly_xxxx_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx);

        auto tlz_xxxx_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx);

        auto tlx_xxxx_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 1);

        auto tly_xxxx_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 1);

        auto tlz_xxxx_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 1);

        auto tlx_xxxx_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 2);

        auto tly_xxxx_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 2);

        auto tlz_xxxx_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 2);

        auto tlx_xxxx_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 3);

        auto tly_xxxx_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 3);

        auto tlz_xxxx_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 3);

        auto tlx_xxxx_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 4);

        auto tly_xxxx_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 4);

        auto tlz_xxxx_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 4);

        auto tlx_xxxx_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 5);

        auto tly_xxxx_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 5);

        auto tlz_xxxx_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 5);

        auto tlx_xxxx_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 6);

        auto tly_xxxx_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 6);

        auto tlz_xxxx_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 6);

        auto tlx_xxxx_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 7);

        auto tly_xxxx_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 7);

        auto tlz_xxxx_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 7);

        auto tlx_xxxx_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 8);

        auto tly_xxxx_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 8);

        auto tlz_xxxx_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 8);

        auto tlx_xxxx_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 9);

        auto tly_xxxx_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 9);

        auto tlz_xxxx_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 9);

        auto tlx_xxxy_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 10);

        auto tly_xxxy_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 10);

        auto tlz_xxxy_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 10);

        auto tlx_xxxy_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 11);

        auto tly_xxxy_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 11);

        auto tlz_xxxy_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 11);

        auto tlx_xxxy_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 12);

        auto tly_xxxy_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 12);

        auto tlz_xxxy_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 12);

        auto tlx_xxxy_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 13);

        auto tly_xxxy_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 13);

        auto tlz_xxxy_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 13);

        auto tlx_xxxy_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 14);

        auto tly_xxxy_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 14);

        auto tlz_xxxy_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 14);

        auto tlx_xxxy_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 15);

        auto tly_xxxy_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 15);

        auto tlz_xxxy_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 15);

        auto tlx_xxxy_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 16);

        auto tly_xxxy_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 16);

        // Batch of Integrals (0,50)

        #pragma omp simd aligned(fgb, fx, pa_x, tdy_xxx_xxx_0, tdy_xxx_xxy_0, tdy_xxx_xxz_0, \
                                     tdy_xxx_xyy_0, tdy_xxx_xyz_0, tdy_xxx_xzz_0, tdy_xxx_yyy_0, tdy_xxx_yyz_0, \
                                     tdy_xxx_yzz_0, tdy_xxx_zzz_0, tdy_xxy_xxx_0, tdy_xxy_xxy_0, tdy_xxy_xxz_0, \
                                     tdy_xxy_xyy_0, tdy_xxy_xyz_0, tdy_xxy_xzz_0, tdz_xxx_xxx_0, tdz_xxx_xxy_0, \
                                     tdz_xxx_xxz_0, tdz_xxx_xyy_0, tdz_xxx_xyz_0, tdz_xxx_xzz_0, tdz_xxx_yyy_0, \
                                     tdz_xxx_yyz_0, tdz_xxx_yzz_0, tdz_xxx_zzz_0, tdz_xxy_xxx_0, tdz_xxy_xxy_0, \
                                     tdz_xxy_xxz_0, tdz_xxy_xyy_0, tdz_xxy_xyz_0, tdz_xxy_xzz_0, tdz_xxy_yyy_0, \
                                     tlx_xx_xxx_0, tlx_xx_xxy_0, tlx_xx_xxz_0, tlx_xx_xyy_0, tlx_xx_xyz_0, tlx_xx_xzz_0, \
                                     tlx_xx_yyy_0, tlx_xx_yyz_0, tlx_xx_yzz_0, tlx_xx_zzz_0, tlx_xxx_xx_0, \
                                     tlx_xxx_xxx_0, tlx_xxx_xxy_0, tlx_xxx_xxz_0, tlx_xxx_xy_0, tlx_xxx_xyy_0, \
                                     tlx_xxx_xyz_0, tlx_xxx_xz_0, tlx_xxx_xzz_0, tlx_xxx_yy_0, tlx_xxx_yyy_0, \
                                     tlx_xxx_yyz_0, tlx_xxx_yz_0, tlx_xxx_yzz_0, tlx_xxx_zz_0, tlx_xxx_zzz_0, \
                                     tlx_xxxx_xxx_0, tlx_xxxx_xxy_0, tlx_xxxx_xxz_0, tlx_xxxx_xyy_0, tlx_xxxx_xyz_0, \
                                     tlx_xxxx_xzz_0, tlx_xxxx_yyy_0, tlx_xxxx_yyz_0, tlx_xxxx_yzz_0, tlx_xxxx_zzz_0, \
                                     tlx_xxxy_xxx_0, tlx_xxxy_xxy_0, tlx_xxxy_xxz_0, tlx_xxxy_xyy_0, tlx_xxxy_xyz_0, \
                                     tlx_xxxy_xzz_0, tlx_xxxy_yyy_0, tlx_xxy_xx_0, tlx_xxy_xxx_0, tlx_xxy_xxy_0, \
                                     tlx_xxy_xxz_0, tlx_xxy_xy_0, tlx_xxy_xyy_0, tlx_xxy_xyz_0, tlx_xxy_xz_0, \
                                     tlx_xxy_xzz_0, tlx_xxy_yy_0, tlx_xxy_yyy_0, tlx_xxy_yz_0, tlx_xxy_zz_0, \
                                     tlx_xy_xxx_0, tlx_xy_xxy_0, tlx_xy_xxz_0, tlx_xy_xyy_0, tlx_xy_xyz_0, tlx_xy_xzz_0, \
                                     tlx_xy_yyy_0, tly_xx_xxx_0, tly_xx_xxy_0, tly_xx_xxz_0, tly_xx_xyy_0, tly_xx_xyz_0, \
                                     tly_xx_xzz_0, tly_xx_yyy_0, tly_xx_yyz_0, tly_xx_yzz_0, tly_xx_zzz_0, tly_xxx_xx_0, \
                                     tly_xxx_xxx_0, tly_xxx_xxy_0, tly_xxx_xxz_0, tly_xxx_xy_0, tly_xxx_xyy_0, \
                                     tly_xxx_xyz_0, tly_xxx_xz_0, tly_xxx_xzz_0, tly_xxx_yy_0, tly_xxx_yyy_0, \
                                     tly_xxx_yyz_0, tly_xxx_yz_0, tly_xxx_yzz_0, tly_xxx_zz_0, tly_xxx_zzz_0, \
                                     tly_xxxx_xxx_0, tly_xxxx_xxy_0, tly_xxxx_xxz_0, tly_xxxx_xyy_0, tly_xxxx_xyz_0, \
                                     tly_xxxx_xzz_0, tly_xxxx_yyy_0, tly_xxxx_yyz_0, tly_xxxx_yzz_0, tly_xxxx_zzz_0, \
                                     tly_xxxy_xxx_0, tly_xxxy_xxy_0, tly_xxxy_xxz_0, tly_xxxy_xyy_0, tly_xxxy_xyz_0, \
                                     tly_xxxy_xzz_0, tly_xxxy_yyy_0, tly_xxy_xx_0, tly_xxy_xxx_0, tly_xxy_xxy_0, \
                                     tly_xxy_xxz_0, tly_xxy_xy_0, tly_xxy_xyy_0, tly_xxy_xyz_0, tly_xxy_xz_0, \
                                     tly_xxy_xzz_0, tly_xxy_yy_0, tly_xxy_yyy_0, tly_xxy_yz_0, tly_xxy_zz_0, \
                                     tly_xy_xxx_0, tly_xy_xxy_0, tly_xy_xxz_0, tly_xy_xyy_0, tly_xy_xyz_0, tly_xy_xzz_0, \
                                     tly_xy_yyy_0, tlz_xx_xxx_0, tlz_xx_xxy_0, tlz_xx_xxz_0, tlz_xx_xyy_0, tlz_xx_xyz_0, \
                                     tlz_xx_xzz_0, tlz_xx_yyy_0, tlz_xx_yyz_0, tlz_xx_yzz_0, tlz_xx_zzz_0, tlz_xxx_xx_0, \
                                     tlz_xxx_xxx_0, tlz_xxx_xxy_0, tlz_xxx_xxz_0, tlz_xxx_xy_0, tlz_xxx_xyy_0, \
                                     tlz_xxx_xyz_0, tlz_xxx_xz_0, tlz_xxx_xzz_0, tlz_xxx_yy_0, tlz_xxx_yyy_0, \
                                     tlz_xxx_yyz_0, tlz_xxx_yz_0, tlz_xxx_yzz_0, tlz_xxx_zz_0, tlz_xxx_zzz_0, \
                                     tlz_xxxx_xxx_0, tlz_xxxx_xxy_0, tlz_xxxx_xxz_0, tlz_xxxx_xyy_0, tlz_xxxx_xyz_0, \
                                     tlz_xxxx_xzz_0, tlz_xxxx_yyy_0, tlz_xxxx_yyz_0, tlz_xxxx_yzz_0, tlz_xxxx_zzz_0, \
                                     tlz_xxxy_xxx_0, tlz_xxxy_xxy_0, tlz_xxxy_xxz_0, tlz_xxxy_xyy_0, tlz_xxxy_xyz_0, \
                                     tlz_xxxy_xzz_0, tlz_xxy_xx_0, tlz_xxy_xxx_0, tlz_xxy_xxy_0, tlz_xxy_xxz_0, \
                                     tlz_xxy_xy_0, tlz_xxy_xyy_0, tlz_xxy_xyz_0, tlz_xxy_xz_0, tlz_xxy_xzz_0, \
                                     tlz_xxy_yy_0, tlz_xxy_yz_0, tlz_xxy_zz_0, tlz_xy_xxx_0, tlz_xy_xxy_0, tlz_xy_xxz_0, \
                                     tlz_xy_xyy_0, tlz_xy_xyz_0, tlz_xy_xzz_0, tpy_xxx_xxx_0, tpy_xxx_xxy_0, \
                                     tpy_xxx_xxz_0, tpy_xxx_xyy_0, tpy_xxx_xyz_0, tpy_xxx_xzz_0, tpy_xxx_yyy_0, \
                                     tpy_xxx_yyz_0, tpy_xxx_yzz_0, tpy_xxx_zzz_0, tpy_xxy_xxx_0, tpy_xxy_xxy_0, \
                                     tpy_xxy_xxz_0, tpy_xxy_xyy_0, tpy_xxy_xyz_0, tpy_xxy_xzz_0, tpz_xxx_xxx_0, \
                                     tpz_xxx_xxy_0, tpz_xxx_xxz_0, tpz_xxx_xyy_0, tpz_xxx_xyz_0, tpz_xxx_xzz_0, \
                                     tpz_xxx_yyy_0, tpz_xxx_yyz_0, tpz_xxx_yzz_0, tpz_xxx_zzz_0, tpz_xxy_xxx_0, \
                                     tpz_xxy_xxy_0, tpz_xxy_xxz_0, tpz_xxy_xyy_0, tpz_xxy_xyz_0, tpz_xxy_xzz_0, \
                                     tpz_xxy_yyy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_xxxx_xxx_0[j] = pa_x[j] * tlx_xxx_xxx_0[j] + 1.5 * fl1_fx * tlx_xx_xxx_0[j] + 1.5 * fl1_fx * tlx_xxx_xx_0[j];

            tly_xxxx_xxx_0[j] = pa_x[j] * tly_xxx_xxx_0[j] + 1.5 * fl1_fx * tly_xx_xxx_0[j] + 1.5 * fl1_fx * tly_xxx_xx_0[j] +
                                0.5 * fl1_fx * tpz_xxx_xxx_0[j] + fl1_fx * fl1_fgb * tdz_xxx_xxx_0[j];

            tlz_xxxx_xxx_0[j] = pa_x[j] * tlz_xxx_xxx_0[j] + 1.5 * fl1_fx * tlz_xx_xxx_0[j] + 1.5 * fl1_fx * tlz_xxx_xx_0[j] -
                                0.5 * fl1_fx * tpy_xxx_xxx_0[j] - fl1_fx * fl1_fgb * tdy_xxx_xxx_0[j];

            tlx_xxxx_xxy_0[j] = pa_x[j] * tlx_xxx_xxy_0[j] + 1.5 * fl1_fx * tlx_xx_xxy_0[j] + fl1_fx * tlx_xxx_xy_0[j];

            tly_xxxx_xxy_0[j] = pa_x[j] * tly_xxx_xxy_0[j] + 1.5 * fl1_fx * tly_xx_xxy_0[j] + fl1_fx * tly_xxx_xy_0[j] +
                                0.5 * fl1_fx * tpz_xxx_xxy_0[j] + fl1_fx * fl1_fgb * tdz_xxx_xxy_0[j];

            tlz_xxxx_xxy_0[j] = pa_x[j] * tlz_xxx_xxy_0[j] + 1.5 * fl1_fx * tlz_xx_xxy_0[j] + fl1_fx * tlz_xxx_xy_0[j] -
                                0.5 * fl1_fx * tpy_xxx_xxy_0[j] - fl1_fx * fl1_fgb * tdy_xxx_xxy_0[j];

            tlx_xxxx_xxz_0[j] = pa_x[j] * tlx_xxx_xxz_0[j] + 1.5 * fl1_fx * tlx_xx_xxz_0[j] + fl1_fx * tlx_xxx_xz_0[j];

            tly_xxxx_xxz_0[j] = pa_x[j] * tly_xxx_xxz_0[j] + 1.5 * fl1_fx * tly_xx_xxz_0[j] + fl1_fx * tly_xxx_xz_0[j] +
                                0.5 * fl1_fx * tpz_xxx_xxz_0[j] + fl1_fx * fl1_fgb * tdz_xxx_xxz_0[j];

            tlz_xxxx_xxz_0[j] = pa_x[j] * tlz_xxx_xxz_0[j] + 1.5 * fl1_fx * tlz_xx_xxz_0[j] + fl1_fx * tlz_xxx_xz_0[j] -
                                0.5 * fl1_fx * tpy_xxx_xxz_0[j] - fl1_fx * fl1_fgb * tdy_xxx_xxz_0[j];

            tlx_xxxx_xyy_0[j] = pa_x[j] * tlx_xxx_xyy_0[j] + 1.5 * fl1_fx * tlx_xx_xyy_0[j] + 0.5 * fl1_fx * tlx_xxx_yy_0[j];

            tly_xxxx_xyy_0[j] = pa_x[j] * tly_xxx_xyy_0[j] + 1.5 * fl1_fx * tly_xx_xyy_0[j] + 0.5 * fl1_fx * tly_xxx_yy_0[j] +
                                0.5 * fl1_fx * tpz_xxx_xyy_0[j] + fl1_fx * fl1_fgb * tdz_xxx_xyy_0[j];

            tlz_xxxx_xyy_0[j] = pa_x[j] * tlz_xxx_xyy_0[j] + 1.5 * fl1_fx * tlz_xx_xyy_0[j] + 0.5 * fl1_fx * tlz_xxx_yy_0[j] -
                                0.5 * fl1_fx * tpy_xxx_xyy_0[j] - fl1_fx * fl1_fgb * tdy_xxx_xyy_0[j];

            tlx_xxxx_xyz_0[j] = pa_x[j] * tlx_xxx_xyz_0[j] + 1.5 * fl1_fx * tlx_xx_xyz_0[j] + 0.5 * fl1_fx * tlx_xxx_yz_0[j];

            tly_xxxx_xyz_0[j] = pa_x[j] * tly_xxx_xyz_0[j] + 1.5 * fl1_fx * tly_xx_xyz_0[j] + 0.5 * fl1_fx * tly_xxx_yz_0[j] +
                                0.5 * fl1_fx * tpz_xxx_xyz_0[j] + fl1_fx * fl1_fgb * tdz_xxx_xyz_0[j];

            tlz_xxxx_xyz_0[j] = pa_x[j] * tlz_xxx_xyz_0[j] + 1.5 * fl1_fx * tlz_xx_xyz_0[j] + 0.5 * fl1_fx * tlz_xxx_yz_0[j] -
                                0.5 * fl1_fx * tpy_xxx_xyz_0[j] - fl1_fx * fl1_fgb * tdy_xxx_xyz_0[j];

            tlx_xxxx_xzz_0[j] = pa_x[j] * tlx_xxx_xzz_0[j] + 1.5 * fl1_fx * tlx_xx_xzz_0[j] + 0.5 * fl1_fx * tlx_xxx_zz_0[j];

            tly_xxxx_xzz_0[j] = pa_x[j] * tly_xxx_xzz_0[j] + 1.5 * fl1_fx * tly_xx_xzz_0[j] + 0.5 * fl1_fx * tly_xxx_zz_0[j] +
                                0.5 * fl1_fx * tpz_xxx_xzz_0[j] + fl1_fx * fl1_fgb * tdz_xxx_xzz_0[j];

            tlz_xxxx_xzz_0[j] = pa_x[j] * tlz_xxx_xzz_0[j] + 1.5 * fl1_fx * tlz_xx_xzz_0[j] + 0.5 * fl1_fx * tlz_xxx_zz_0[j] -
                                0.5 * fl1_fx * tpy_xxx_xzz_0[j] - fl1_fx * fl1_fgb * tdy_xxx_xzz_0[j];

            tlx_xxxx_yyy_0[j] = pa_x[j] * tlx_xxx_yyy_0[j] + 1.5 * fl1_fx * tlx_xx_yyy_0[j];

            tly_xxxx_yyy_0[j] =
                pa_x[j] * tly_xxx_yyy_0[j] + 1.5 * fl1_fx * tly_xx_yyy_0[j] + 0.5 * fl1_fx * tpz_xxx_yyy_0[j] + fl1_fx * fl1_fgb * tdz_xxx_yyy_0[j];

            tlz_xxxx_yyy_0[j] =
                pa_x[j] * tlz_xxx_yyy_0[j] + 1.5 * fl1_fx * tlz_xx_yyy_0[j] - 0.5 * fl1_fx * tpy_xxx_yyy_0[j] - fl1_fx * fl1_fgb * tdy_xxx_yyy_0[j];

            tlx_xxxx_yyz_0[j] = pa_x[j] * tlx_xxx_yyz_0[j] + 1.5 * fl1_fx * tlx_xx_yyz_0[j];

            tly_xxxx_yyz_0[j] =
                pa_x[j] * tly_xxx_yyz_0[j] + 1.5 * fl1_fx * tly_xx_yyz_0[j] + 0.5 * fl1_fx * tpz_xxx_yyz_0[j] + fl1_fx * fl1_fgb * tdz_xxx_yyz_0[j];

            tlz_xxxx_yyz_0[j] =
                pa_x[j] * tlz_xxx_yyz_0[j] + 1.5 * fl1_fx * tlz_xx_yyz_0[j] - 0.5 * fl1_fx * tpy_xxx_yyz_0[j] - fl1_fx * fl1_fgb * tdy_xxx_yyz_0[j];

            tlx_xxxx_yzz_0[j] = pa_x[j] * tlx_xxx_yzz_0[j] + 1.5 * fl1_fx * tlx_xx_yzz_0[j];

            tly_xxxx_yzz_0[j] =
                pa_x[j] * tly_xxx_yzz_0[j] + 1.5 * fl1_fx * tly_xx_yzz_0[j] + 0.5 * fl1_fx * tpz_xxx_yzz_0[j] + fl1_fx * fl1_fgb * tdz_xxx_yzz_0[j];

            tlz_xxxx_yzz_0[j] =
                pa_x[j] * tlz_xxx_yzz_0[j] + 1.5 * fl1_fx * tlz_xx_yzz_0[j] - 0.5 * fl1_fx * tpy_xxx_yzz_0[j] - fl1_fx * fl1_fgb * tdy_xxx_yzz_0[j];

            tlx_xxxx_zzz_0[j] = pa_x[j] * tlx_xxx_zzz_0[j] + 1.5 * fl1_fx * tlx_xx_zzz_0[j];

            tly_xxxx_zzz_0[j] =
                pa_x[j] * tly_xxx_zzz_0[j] + 1.5 * fl1_fx * tly_xx_zzz_0[j] + 0.5 * fl1_fx * tpz_xxx_zzz_0[j] + fl1_fx * fl1_fgb * tdz_xxx_zzz_0[j];

            tlz_xxxx_zzz_0[j] =
                pa_x[j] * tlz_xxx_zzz_0[j] + 1.5 * fl1_fx * tlz_xx_zzz_0[j] - 0.5 * fl1_fx * tpy_xxx_zzz_0[j] - fl1_fx * fl1_fgb * tdy_xxx_zzz_0[j];

            tlx_xxxy_xxx_0[j] = pa_x[j] * tlx_xxy_xxx_0[j] + fl1_fx * tlx_xy_xxx_0[j] + 1.5 * fl1_fx * tlx_xxy_xx_0[j];

            tly_xxxy_xxx_0[j] = pa_x[j] * tly_xxy_xxx_0[j] + fl1_fx * tly_xy_xxx_0[j] + 1.5 * fl1_fx * tly_xxy_xx_0[j] +
                                0.5 * fl1_fx * tpz_xxy_xxx_0[j] + fl1_fx * fl1_fgb * tdz_xxy_xxx_0[j];

            tlz_xxxy_xxx_0[j] = pa_x[j] * tlz_xxy_xxx_0[j] + fl1_fx * tlz_xy_xxx_0[j] + 1.5 * fl1_fx * tlz_xxy_xx_0[j] -
                                0.5 * fl1_fx * tpy_xxy_xxx_0[j] - fl1_fx * fl1_fgb * tdy_xxy_xxx_0[j];

            tlx_xxxy_xxy_0[j] = pa_x[j] * tlx_xxy_xxy_0[j] + fl1_fx * tlx_xy_xxy_0[j] + fl1_fx * tlx_xxy_xy_0[j];

            tly_xxxy_xxy_0[j] = pa_x[j] * tly_xxy_xxy_0[j] + fl1_fx * tly_xy_xxy_0[j] + fl1_fx * tly_xxy_xy_0[j] + 0.5 * fl1_fx * tpz_xxy_xxy_0[j] +
                                fl1_fx * fl1_fgb * tdz_xxy_xxy_0[j];

            tlz_xxxy_xxy_0[j] = pa_x[j] * tlz_xxy_xxy_0[j] + fl1_fx * tlz_xy_xxy_0[j] + fl1_fx * tlz_xxy_xy_0[j] - 0.5 * fl1_fx * tpy_xxy_xxy_0[j] -
                                fl1_fx * fl1_fgb * tdy_xxy_xxy_0[j];

            tlx_xxxy_xxz_0[j] = pa_x[j] * tlx_xxy_xxz_0[j] + fl1_fx * tlx_xy_xxz_0[j] + fl1_fx * tlx_xxy_xz_0[j];

            tly_xxxy_xxz_0[j] = pa_x[j] * tly_xxy_xxz_0[j] + fl1_fx * tly_xy_xxz_0[j] + fl1_fx * tly_xxy_xz_0[j] + 0.5 * fl1_fx * tpz_xxy_xxz_0[j] +
                                fl1_fx * fl1_fgb * tdz_xxy_xxz_0[j];

            tlz_xxxy_xxz_0[j] = pa_x[j] * tlz_xxy_xxz_0[j] + fl1_fx * tlz_xy_xxz_0[j] + fl1_fx * tlz_xxy_xz_0[j] - 0.5 * fl1_fx * tpy_xxy_xxz_0[j] -
                                fl1_fx * fl1_fgb * tdy_xxy_xxz_0[j];

            tlx_xxxy_xyy_0[j] = pa_x[j] * tlx_xxy_xyy_0[j] + fl1_fx * tlx_xy_xyy_0[j] + 0.5 * fl1_fx * tlx_xxy_yy_0[j];

            tly_xxxy_xyy_0[j] = pa_x[j] * tly_xxy_xyy_0[j] + fl1_fx * tly_xy_xyy_0[j] + 0.5 * fl1_fx * tly_xxy_yy_0[j] +
                                0.5 * fl1_fx * tpz_xxy_xyy_0[j] + fl1_fx * fl1_fgb * tdz_xxy_xyy_0[j];

            tlz_xxxy_xyy_0[j] = pa_x[j] * tlz_xxy_xyy_0[j] + fl1_fx * tlz_xy_xyy_0[j] + 0.5 * fl1_fx * tlz_xxy_yy_0[j] -
                                0.5 * fl1_fx * tpy_xxy_xyy_0[j] - fl1_fx * fl1_fgb * tdy_xxy_xyy_0[j];

            tlx_xxxy_xyz_0[j] = pa_x[j] * tlx_xxy_xyz_0[j] + fl1_fx * tlx_xy_xyz_0[j] + 0.5 * fl1_fx * tlx_xxy_yz_0[j];

            tly_xxxy_xyz_0[j] = pa_x[j] * tly_xxy_xyz_0[j] + fl1_fx * tly_xy_xyz_0[j] + 0.5 * fl1_fx * tly_xxy_yz_0[j] +
                                0.5 * fl1_fx * tpz_xxy_xyz_0[j] + fl1_fx * fl1_fgb * tdz_xxy_xyz_0[j];

            tlz_xxxy_xyz_0[j] = pa_x[j] * tlz_xxy_xyz_0[j] + fl1_fx * tlz_xy_xyz_0[j] + 0.5 * fl1_fx * tlz_xxy_yz_0[j] -
                                0.5 * fl1_fx * tpy_xxy_xyz_0[j] - fl1_fx * fl1_fgb * tdy_xxy_xyz_0[j];

            tlx_xxxy_xzz_0[j] = pa_x[j] * tlx_xxy_xzz_0[j] + fl1_fx * tlx_xy_xzz_0[j] + 0.5 * fl1_fx * tlx_xxy_zz_0[j];

            tly_xxxy_xzz_0[j] = pa_x[j] * tly_xxy_xzz_0[j] + fl1_fx * tly_xy_xzz_0[j] + 0.5 * fl1_fx * tly_xxy_zz_0[j] +
                                0.5 * fl1_fx * tpz_xxy_xzz_0[j] + fl1_fx * fl1_fgb * tdz_xxy_xzz_0[j];

            tlz_xxxy_xzz_0[j] = pa_x[j] * tlz_xxy_xzz_0[j] + fl1_fx * tlz_xy_xzz_0[j] + 0.5 * fl1_fx * tlz_xxy_zz_0[j] -
                                0.5 * fl1_fx * tpy_xxy_xzz_0[j] - fl1_fx * fl1_fgb * tdy_xxy_xzz_0[j];

            tlx_xxxy_yyy_0[j] = pa_x[j] * tlx_xxy_yyy_0[j] + fl1_fx * tlx_xy_yyy_0[j];

            tly_xxxy_yyy_0[j] =
                pa_x[j] * tly_xxy_yyy_0[j] + fl1_fx * tly_xy_yyy_0[j] + 0.5 * fl1_fx * tpz_xxy_yyy_0[j] + fl1_fx * fl1_fgb * tdz_xxy_yyy_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForGF_50_100(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_l_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_4_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

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

        auto tlx_xxz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 12);

        auto tly_xxz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 12);

        auto tlz_xxz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 12);

        auto tlx_xxz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 13);

        auto tly_xxz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 13);

        auto tlz_xxz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 13);

        auto tlx_xxz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 14);

        auto tly_xxz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 14);

        auto tlz_xxz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 14);

        auto tlx_xxz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 15);

        auto tly_xxz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 15);

        auto tlz_xxz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 15);

        auto tlx_xxz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 16);

        auto tly_xxz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 16);

        auto tlz_xxz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 16);

        auto tlx_xxz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 17);

        auto tly_xxz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 17);

        auto tlz_xxz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 17);

        auto tlx_xyy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 18);

        auto tly_xyy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 18);

        auto tlz_xyy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 18);

        auto tlx_xyy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 19);

        auto tly_xyy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 19);

        auto tlz_xyy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 19);

        auto tlx_xyy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 20);

        auto tly_xyy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 20);

        auto tlz_xyy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 20);

        auto tlx_xyy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 21);

        auto tpy_xxy_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 16);

        auto tpy_xxy_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 17);

        auto tpz_xxy_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 17);

        auto tpy_xxy_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 18);

        auto tpz_xxy_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 18);

        auto tpy_xxy_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 19);

        auto tpz_xxy_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 19);

        auto tpy_xxz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 20);

        auto tpz_xxz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 20);

        auto tpy_xxz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 21);

        auto tpz_xxz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 21);

        auto tpy_xxz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 22);

        auto tpz_xxz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 22);

        auto tpy_xxz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 23);

        auto tpz_xxz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 23);

        auto tpy_xxz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 24);

        auto tpz_xxz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 24);

        auto tpy_xxz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 25);

        auto tpz_xxz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 25);

        auto tpy_xxz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 26);

        auto tpz_xxz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 26);

        auto tpy_xxz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 27);

        auto tpz_xxz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 27);

        auto tpy_xxz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 28);

        auto tpz_xxz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 28);

        auto tpy_xxz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 29);

        auto tpz_xxz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 29);

        auto tpy_xyy_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 30);

        auto tpz_xyy_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 30);

        auto tpy_xyy_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 31);

        auto tpz_xyy_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 31);

        auto tpy_xyy_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 32);

        auto tpz_xyy_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 32);

        auto tdy_xxy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 16);

        auto tdy_xxy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 17);

        auto tdz_xxy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 17);

        auto tdy_xxy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 18);

        auto tdz_xxy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 18);

        auto tdy_xxy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 19);

        auto tdz_xxy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 19);

        auto tdy_xxz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 20);

        auto tdz_xxz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 20);

        auto tdy_xxz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 21);

        auto tdz_xxz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 21);

        auto tdy_xxz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 22);

        auto tdz_xxz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 22);

        auto tdy_xxz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 23);

        auto tdz_xxz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 23);

        auto tdy_xxz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 24);

        auto tdz_xxz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 24);

        auto tdy_xxz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 25);

        auto tdz_xxz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 25);

        auto tdy_xxz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 26);

        auto tdz_xxz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 26);

        auto tdy_xxz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 27);

        auto tdz_xxz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 27);

        auto tdy_xxz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 28);

        auto tdz_xxz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 28);

        auto tdy_xxz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 29);

        auto tdz_xxz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 29);

        auto tdy_xyy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 30);

        auto tdz_xyy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 30);

        auto tdy_xyy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 31);

        auto tdz_xyy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 31);

        auto tdy_xyy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 32);

        auto tdz_xyy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 32);

        // set up pointers to integrals

        auto tlz_xxxy_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 16);

        auto tlx_xxxy_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 17);

        auto tly_xxxy_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 17);

        auto tlz_xxxy_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 17);

        auto tlx_xxxy_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 18);

        auto tly_xxxy_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 18);

        auto tlz_xxxy_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 18);

        auto tlx_xxxy_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 19);

        auto tly_xxxy_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 19);

        auto tlz_xxxy_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 19);

        auto tlx_xxxz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 20);

        auto tly_xxxz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 20);

        auto tlz_xxxz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 20);

        auto tlx_xxxz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 21);

        auto tly_xxxz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 21);

        auto tlz_xxxz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 21);

        auto tlx_xxxz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 22);

        auto tly_xxxz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 22);

        auto tlz_xxxz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 22);

        auto tlx_xxxz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 23);

        auto tly_xxxz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 23);

        auto tlz_xxxz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 23);

        auto tlx_xxxz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 24);

        auto tly_xxxz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 24);

        auto tlz_xxxz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 24);

        auto tlx_xxxz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 25);

        auto tly_xxxz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 25);

        auto tlz_xxxz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 25);

        auto tlx_xxxz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 26);

        auto tly_xxxz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 26);

        auto tlz_xxxz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 26);

        auto tlx_xxxz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 27);

        auto tly_xxxz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 27);

        auto tlz_xxxz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 27);

        auto tlx_xxxz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 28);

        auto tly_xxxz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 28);

        auto tlz_xxxz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 28);

        auto tlx_xxxz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 29);

        auto tly_xxxz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 29);

        auto tlz_xxxz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 29);

        auto tlx_xxyy_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 30);

        auto tly_xxyy_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 30);

        auto tlz_xxyy_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 30);

        auto tlx_xxyy_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 31);

        auto tly_xxyy_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 31);

        auto tlz_xxyy_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 31);

        auto tlx_xxyy_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 32);

        auto tly_xxyy_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 32);

        auto tlz_xxyy_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 32);

        auto tlx_xxyy_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 33);

        // Batch of Integrals (50,100)

        #pragma omp simd aligned(fgb, fx, pa_x, tdy_xxy_yyy_0, tdy_xxy_yyz_0, tdy_xxy_yzz_0, \
                                     tdy_xxy_zzz_0, tdy_xxz_xxx_0, tdy_xxz_xxy_0, tdy_xxz_xxz_0, tdy_xxz_xyy_0, \
                                     tdy_xxz_xyz_0, tdy_xxz_xzz_0, tdy_xxz_yyy_0, tdy_xxz_yyz_0, tdy_xxz_yzz_0, \
                                     tdy_xxz_zzz_0, tdy_xyy_xxx_0, tdy_xyy_xxy_0, tdy_xyy_xxz_0, tdz_xxy_yyz_0, \
                                     tdz_xxy_yzz_0, tdz_xxy_zzz_0, tdz_xxz_xxx_0, tdz_xxz_xxy_0, tdz_xxz_xxz_0, \
                                     tdz_xxz_xyy_0, tdz_xxz_xyz_0, tdz_xxz_xzz_0, tdz_xxz_yyy_0, tdz_xxz_yyz_0, \
                                     tdz_xxz_yzz_0, tdz_xxz_zzz_0, tdz_xyy_xxx_0, tdz_xyy_xxy_0, tdz_xyy_xxz_0, \
                                     tlx_xxxy_yyz_0, tlx_xxxy_yzz_0, tlx_xxxy_zzz_0, tlx_xxxz_xxx_0, tlx_xxxz_xxy_0, \
                                     tlx_xxxz_xxz_0, tlx_xxxz_xyy_0, tlx_xxxz_xyz_0, tlx_xxxz_xzz_0, tlx_xxxz_yyy_0, \
                                     tlx_xxxz_yyz_0, tlx_xxxz_yzz_0, tlx_xxxz_zzz_0, tlx_xxy_yyz_0, tlx_xxy_yzz_0, \
                                     tlx_xxy_zzz_0, tlx_xxyy_xxx_0, tlx_xxyy_xxy_0, tlx_xxyy_xxz_0, tlx_xxyy_xyy_0, \
                                     tlx_xxz_xx_0, tlx_xxz_xxx_0, tlx_xxz_xxy_0, tlx_xxz_xxz_0, tlx_xxz_xy_0, \
                                     tlx_xxz_xyy_0, tlx_xxz_xyz_0, tlx_xxz_xz_0, tlx_xxz_xzz_0, tlx_xxz_yy_0, \
                                     tlx_xxz_yyy_0, tlx_xxz_yyz_0, tlx_xxz_yz_0, tlx_xxz_yzz_0, tlx_xxz_zz_0, \
                                     tlx_xxz_zzz_0, tlx_xy_yyz_0, tlx_xy_yzz_0, tlx_xy_zzz_0, tlx_xyy_xx_0, \
                                     tlx_xyy_xxx_0, tlx_xyy_xxy_0, tlx_xyy_xxz_0, tlx_xyy_xy_0, tlx_xyy_xyy_0, \
                                     tlx_xyy_xz_0, tlx_xyy_yy_0, tlx_xz_xxx_0, tlx_xz_xxy_0, tlx_xz_xxz_0, tlx_xz_xyy_0, \
                                     tlx_xz_xyz_0, tlx_xz_xzz_0, tlx_xz_yyy_0, tlx_xz_yyz_0, tlx_xz_yzz_0, tlx_xz_zzz_0, \
                                     tlx_yy_xxx_0, tlx_yy_xxy_0, tlx_yy_xxz_0, tlx_yy_xyy_0, tly_xxxy_yyz_0, \
                                     tly_xxxy_yzz_0, tly_xxxy_zzz_0, tly_xxxz_xxx_0, tly_xxxz_xxy_0, tly_xxxz_xxz_0, \
                                     tly_xxxz_xyy_0, tly_xxxz_xyz_0, tly_xxxz_xzz_0, tly_xxxz_yyy_0, tly_xxxz_yyz_0, \
                                     tly_xxxz_yzz_0, tly_xxxz_zzz_0, tly_xxy_yyz_0, tly_xxy_yzz_0, tly_xxy_zzz_0, \
                                     tly_xxyy_xxx_0, tly_xxyy_xxy_0, tly_xxyy_xxz_0, tly_xxz_xx_0, tly_xxz_xxx_0, \
                                     tly_xxz_xxy_0, tly_xxz_xxz_0, tly_xxz_xy_0, tly_xxz_xyy_0, tly_xxz_xyz_0, \
                                     tly_xxz_xz_0, tly_xxz_xzz_0, tly_xxz_yy_0, tly_xxz_yyy_0, tly_xxz_yyz_0, \
                                     tly_xxz_yz_0, tly_xxz_yzz_0, tly_xxz_zz_0, tly_xxz_zzz_0, tly_xy_yyz_0, \
                                     tly_xy_yzz_0, tly_xy_zzz_0, tly_xyy_xx_0, tly_xyy_xxx_0, tly_xyy_xxy_0, \
                                     tly_xyy_xxz_0, tly_xyy_xy_0, tly_xyy_xz_0, tly_xz_xxx_0, tly_xz_xxy_0, tly_xz_xxz_0, \
                                     tly_xz_xyy_0, tly_xz_xyz_0, tly_xz_xzz_0, tly_xz_yyy_0, tly_xz_yyz_0, tly_xz_yzz_0, \
                                     tly_xz_zzz_0, tly_yy_xxx_0, tly_yy_xxy_0, tly_yy_xxz_0, tlz_xxxy_yyy_0, \
                                     tlz_xxxy_yyz_0, tlz_xxxy_yzz_0, tlz_xxxy_zzz_0, tlz_xxxz_xxx_0, tlz_xxxz_xxy_0, \
                                     tlz_xxxz_xxz_0, tlz_xxxz_xyy_0, tlz_xxxz_xyz_0, tlz_xxxz_xzz_0, tlz_xxxz_yyy_0, \
                                     tlz_xxxz_yyz_0, tlz_xxxz_yzz_0, tlz_xxxz_zzz_0, tlz_xxy_yyy_0, tlz_xxy_yyz_0, \
                                     tlz_xxy_yzz_0, tlz_xxy_zzz_0, tlz_xxyy_xxx_0, tlz_xxyy_xxy_0, tlz_xxyy_xxz_0, \
                                     tlz_xxz_xx_0, tlz_xxz_xxx_0, tlz_xxz_xxy_0, tlz_xxz_xxz_0, tlz_xxz_xy_0, \
                                     tlz_xxz_xyy_0, tlz_xxz_xyz_0, tlz_xxz_xz_0, tlz_xxz_xzz_0, tlz_xxz_yy_0, \
                                     tlz_xxz_yyy_0, tlz_xxz_yyz_0, tlz_xxz_yz_0, tlz_xxz_yzz_0, tlz_xxz_zz_0, \
                                     tlz_xxz_zzz_0, tlz_xy_yyy_0, tlz_xy_yyz_0, tlz_xy_yzz_0, tlz_xy_zzz_0, tlz_xyy_xx_0, \
                                     tlz_xyy_xxx_0, tlz_xyy_xxy_0, tlz_xyy_xxz_0, tlz_xyy_xy_0, tlz_xyy_xz_0, \
                                     tlz_xz_xxx_0, tlz_xz_xxy_0, tlz_xz_xxz_0, tlz_xz_xyy_0, tlz_xz_xyz_0, tlz_xz_xzz_0, \
                                     tlz_xz_yyy_0, tlz_xz_yyz_0, tlz_xz_yzz_0, tlz_xz_zzz_0, tlz_yy_xxx_0, tlz_yy_xxy_0, \
                                     tlz_yy_xxz_0, tpy_xxy_yyy_0, tpy_xxy_yyz_0, tpy_xxy_yzz_0, tpy_xxy_zzz_0, \
                                     tpy_xxz_xxx_0, tpy_xxz_xxy_0, tpy_xxz_xxz_0, tpy_xxz_xyy_0, tpy_xxz_xyz_0, \
                                     tpy_xxz_xzz_0, tpy_xxz_yyy_0, tpy_xxz_yyz_0, tpy_xxz_yzz_0, tpy_xxz_zzz_0, \
                                     tpy_xyy_xxx_0, tpy_xyy_xxy_0, tpy_xyy_xxz_0, tpz_xxy_yyz_0, tpz_xxy_yzz_0, \
                                     tpz_xxy_zzz_0, tpz_xxz_xxx_0, tpz_xxz_xxy_0, tpz_xxz_xxz_0, tpz_xxz_xyy_0, \
                                     tpz_xxz_xyz_0, tpz_xxz_xzz_0, tpz_xxz_yyy_0, tpz_xxz_yyz_0, tpz_xxz_yzz_0, \
                                     tpz_xxz_zzz_0, tpz_xyy_xxx_0, tpz_xyy_xxy_0, tpz_xyy_xxz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlz_xxxy_yyy_0[j] =
                pa_x[j] * tlz_xxy_yyy_0[j] + fl1_fx * tlz_xy_yyy_0[j] - 0.5 * fl1_fx * tpy_xxy_yyy_0[j] - fl1_fx * fl1_fgb * tdy_xxy_yyy_0[j];

            tlx_xxxy_yyz_0[j] = pa_x[j] * tlx_xxy_yyz_0[j] + fl1_fx * tlx_xy_yyz_0[j];

            tly_xxxy_yyz_0[j] =
                pa_x[j] * tly_xxy_yyz_0[j] + fl1_fx * tly_xy_yyz_0[j] + 0.5 * fl1_fx * tpz_xxy_yyz_0[j] + fl1_fx * fl1_fgb * tdz_xxy_yyz_0[j];

            tlz_xxxy_yyz_0[j] =
                pa_x[j] * tlz_xxy_yyz_0[j] + fl1_fx * tlz_xy_yyz_0[j] - 0.5 * fl1_fx * tpy_xxy_yyz_0[j] - fl1_fx * fl1_fgb * tdy_xxy_yyz_0[j];

            tlx_xxxy_yzz_0[j] = pa_x[j] * tlx_xxy_yzz_0[j] + fl1_fx * tlx_xy_yzz_0[j];

            tly_xxxy_yzz_0[j] =
                pa_x[j] * tly_xxy_yzz_0[j] + fl1_fx * tly_xy_yzz_0[j] + 0.5 * fl1_fx * tpz_xxy_yzz_0[j] + fl1_fx * fl1_fgb * tdz_xxy_yzz_0[j];

            tlz_xxxy_yzz_0[j] =
                pa_x[j] * tlz_xxy_yzz_0[j] + fl1_fx * tlz_xy_yzz_0[j] - 0.5 * fl1_fx * tpy_xxy_yzz_0[j] - fl1_fx * fl1_fgb * tdy_xxy_yzz_0[j];

            tlx_xxxy_zzz_0[j] = pa_x[j] * tlx_xxy_zzz_0[j] + fl1_fx * tlx_xy_zzz_0[j];

            tly_xxxy_zzz_0[j] =
                pa_x[j] * tly_xxy_zzz_0[j] + fl1_fx * tly_xy_zzz_0[j] + 0.5 * fl1_fx * tpz_xxy_zzz_0[j] + fl1_fx * fl1_fgb * tdz_xxy_zzz_0[j];

            tlz_xxxy_zzz_0[j] =
                pa_x[j] * tlz_xxy_zzz_0[j] + fl1_fx * tlz_xy_zzz_0[j] - 0.5 * fl1_fx * tpy_xxy_zzz_0[j] - fl1_fx * fl1_fgb * tdy_xxy_zzz_0[j];

            tlx_xxxz_xxx_0[j] = pa_x[j] * tlx_xxz_xxx_0[j] + fl1_fx * tlx_xz_xxx_0[j] + 1.5 * fl1_fx * tlx_xxz_xx_0[j];

            tly_xxxz_xxx_0[j] = pa_x[j] * tly_xxz_xxx_0[j] + fl1_fx * tly_xz_xxx_0[j] + 1.5 * fl1_fx * tly_xxz_xx_0[j] +
                                0.5 * fl1_fx * tpz_xxz_xxx_0[j] + fl1_fx * fl1_fgb * tdz_xxz_xxx_0[j];

            tlz_xxxz_xxx_0[j] = pa_x[j] * tlz_xxz_xxx_0[j] + fl1_fx * tlz_xz_xxx_0[j] + 1.5 * fl1_fx * tlz_xxz_xx_0[j] -
                                0.5 * fl1_fx * tpy_xxz_xxx_0[j] - fl1_fx * fl1_fgb * tdy_xxz_xxx_0[j];

            tlx_xxxz_xxy_0[j] = pa_x[j] * tlx_xxz_xxy_0[j] + fl1_fx * tlx_xz_xxy_0[j] + fl1_fx * tlx_xxz_xy_0[j];

            tly_xxxz_xxy_0[j] = pa_x[j] * tly_xxz_xxy_0[j] + fl1_fx * tly_xz_xxy_0[j] + fl1_fx * tly_xxz_xy_0[j] + 0.5 * fl1_fx * tpz_xxz_xxy_0[j] +
                                fl1_fx * fl1_fgb * tdz_xxz_xxy_0[j];

            tlz_xxxz_xxy_0[j] = pa_x[j] * tlz_xxz_xxy_0[j] + fl1_fx * tlz_xz_xxy_0[j] + fl1_fx * tlz_xxz_xy_0[j] - 0.5 * fl1_fx * tpy_xxz_xxy_0[j] -
                                fl1_fx * fl1_fgb * tdy_xxz_xxy_0[j];

            tlx_xxxz_xxz_0[j] = pa_x[j] * tlx_xxz_xxz_0[j] + fl1_fx * tlx_xz_xxz_0[j] + fl1_fx * tlx_xxz_xz_0[j];

            tly_xxxz_xxz_0[j] = pa_x[j] * tly_xxz_xxz_0[j] + fl1_fx * tly_xz_xxz_0[j] + fl1_fx * tly_xxz_xz_0[j] + 0.5 * fl1_fx * tpz_xxz_xxz_0[j] +
                                fl1_fx * fl1_fgb * tdz_xxz_xxz_0[j];

            tlz_xxxz_xxz_0[j] = pa_x[j] * tlz_xxz_xxz_0[j] + fl1_fx * tlz_xz_xxz_0[j] + fl1_fx * tlz_xxz_xz_0[j] - 0.5 * fl1_fx * tpy_xxz_xxz_0[j] -
                                fl1_fx * fl1_fgb * tdy_xxz_xxz_0[j];

            tlx_xxxz_xyy_0[j] = pa_x[j] * tlx_xxz_xyy_0[j] + fl1_fx * tlx_xz_xyy_0[j] + 0.5 * fl1_fx * tlx_xxz_yy_0[j];

            tly_xxxz_xyy_0[j] = pa_x[j] * tly_xxz_xyy_0[j] + fl1_fx * tly_xz_xyy_0[j] + 0.5 * fl1_fx * tly_xxz_yy_0[j] +
                                0.5 * fl1_fx * tpz_xxz_xyy_0[j] + fl1_fx * fl1_fgb * tdz_xxz_xyy_0[j];

            tlz_xxxz_xyy_0[j] = pa_x[j] * tlz_xxz_xyy_0[j] + fl1_fx * tlz_xz_xyy_0[j] + 0.5 * fl1_fx * tlz_xxz_yy_0[j] -
                                0.5 * fl1_fx * tpy_xxz_xyy_0[j] - fl1_fx * fl1_fgb * tdy_xxz_xyy_0[j];

            tlx_xxxz_xyz_0[j] = pa_x[j] * tlx_xxz_xyz_0[j] + fl1_fx * tlx_xz_xyz_0[j] + 0.5 * fl1_fx * tlx_xxz_yz_0[j];

            tly_xxxz_xyz_0[j] = pa_x[j] * tly_xxz_xyz_0[j] + fl1_fx * tly_xz_xyz_0[j] + 0.5 * fl1_fx * tly_xxz_yz_0[j] +
                                0.5 * fl1_fx * tpz_xxz_xyz_0[j] + fl1_fx * fl1_fgb * tdz_xxz_xyz_0[j];

            tlz_xxxz_xyz_0[j] = pa_x[j] * tlz_xxz_xyz_0[j] + fl1_fx * tlz_xz_xyz_0[j] + 0.5 * fl1_fx * tlz_xxz_yz_0[j] -
                                0.5 * fl1_fx * tpy_xxz_xyz_0[j] - fl1_fx * fl1_fgb * tdy_xxz_xyz_0[j];

            tlx_xxxz_xzz_0[j] = pa_x[j] * tlx_xxz_xzz_0[j] + fl1_fx * tlx_xz_xzz_0[j] + 0.5 * fl1_fx * tlx_xxz_zz_0[j];

            tly_xxxz_xzz_0[j] = pa_x[j] * tly_xxz_xzz_0[j] + fl1_fx * tly_xz_xzz_0[j] + 0.5 * fl1_fx * tly_xxz_zz_0[j] +
                                0.5 * fl1_fx * tpz_xxz_xzz_0[j] + fl1_fx * fl1_fgb * tdz_xxz_xzz_0[j];

            tlz_xxxz_xzz_0[j] = pa_x[j] * tlz_xxz_xzz_0[j] + fl1_fx * tlz_xz_xzz_0[j] + 0.5 * fl1_fx * tlz_xxz_zz_0[j] -
                                0.5 * fl1_fx * tpy_xxz_xzz_0[j] - fl1_fx * fl1_fgb * tdy_xxz_xzz_0[j];

            tlx_xxxz_yyy_0[j] = pa_x[j] * tlx_xxz_yyy_0[j] + fl1_fx * tlx_xz_yyy_0[j];

            tly_xxxz_yyy_0[j] =
                pa_x[j] * tly_xxz_yyy_0[j] + fl1_fx * tly_xz_yyy_0[j] + 0.5 * fl1_fx * tpz_xxz_yyy_0[j] + fl1_fx * fl1_fgb * tdz_xxz_yyy_0[j];

            tlz_xxxz_yyy_0[j] =
                pa_x[j] * tlz_xxz_yyy_0[j] + fl1_fx * tlz_xz_yyy_0[j] - 0.5 * fl1_fx * tpy_xxz_yyy_0[j] - fl1_fx * fl1_fgb * tdy_xxz_yyy_0[j];

            tlx_xxxz_yyz_0[j] = pa_x[j] * tlx_xxz_yyz_0[j] + fl1_fx * tlx_xz_yyz_0[j];

            tly_xxxz_yyz_0[j] =
                pa_x[j] * tly_xxz_yyz_0[j] + fl1_fx * tly_xz_yyz_0[j] + 0.5 * fl1_fx * tpz_xxz_yyz_0[j] + fl1_fx * fl1_fgb * tdz_xxz_yyz_0[j];

            tlz_xxxz_yyz_0[j] =
                pa_x[j] * tlz_xxz_yyz_0[j] + fl1_fx * tlz_xz_yyz_0[j] - 0.5 * fl1_fx * tpy_xxz_yyz_0[j] - fl1_fx * fl1_fgb * tdy_xxz_yyz_0[j];

            tlx_xxxz_yzz_0[j] = pa_x[j] * tlx_xxz_yzz_0[j] + fl1_fx * tlx_xz_yzz_0[j];

            tly_xxxz_yzz_0[j] =
                pa_x[j] * tly_xxz_yzz_0[j] + fl1_fx * tly_xz_yzz_0[j] + 0.5 * fl1_fx * tpz_xxz_yzz_0[j] + fl1_fx * fl1_fgb * tdz_xxz_yzz_0[j];

            tlz_xxxz_yzz_0[j] =
                pa_x[j] * tlz_xxz_yzz_0[j] + fl1_fx * tlz_xz_yzz_0[j] - 0.5 * fl1_fx * tpy_xxz_yzz_0[j] - fl1_fx * fl1_fgb * tdy_xxz_yzz_0[j];

            tlx_xxxz_zzz_0[j] = pa_x[j] * tlx_xxz_zzz_0[j] + fl1_fx * tlx_xz_zzz_0[j];

            tly_xxxz_zzz_0[j] =
                pa_x[j] * tly_xxz_zzz_0[j] + fl1_fx * tly_xz_zzz_0[j] + 0.5 * fl1_fx * tpz_xxz_zzz_0[j] + fl1_fx * fl1_fgb * tdz_xxz_zzz_0[j];

            tlz_xxxz_zzz_0[j] =
                pa_x[j] * tlz_xxz_zzz_0[j] + fl1_fx * tlz_xz_zzz_0[j] - 0.5 * fl1_fx * tpy_xxz_zzz_0[j] - fl1_fx * fl1_fgb * tdy_xxz_zzz_0[j];

            tlx_xxyy_xxx_0[j] = pa_x[j] * tlx_xyy_xxx_0[j] + 0.5 * fl1_fx * tlx_yy_xxx_0[j] + 1.5 * fl1_fx * tlx_xyy_xx_0[j];

            tly_xxyy_xxx_0[j] = pa_x[j] * tly_xyy_xxx_0[j] + 0.5 * fl1_fx * tly_yy_xxx_0[j] + 1.5 * fl1_fx * tly_xyy_xx_0[j] +
                                0.5 * fl1_fx * tpz_xyy_xxx_0[j] + fl1_fx * fl1_fgb * tdz_xyy_xxx_0[j];

            tlz_xxyy_xxx_0[j] = pa_x[j] * tlz_xyy_xxx_0[j] + 0.5 * fl1_fx * tlz_yy_xxx_0[j] + 1.5 * fl1_fx * tlz_xyy_xx_0[j] -
                                0.5 * fl1_fx * tpy_xyy_xxx_0[j] - fl1_fx * fl1_fgb * tdy_xyy_xxx_0[j];

            tlx_xxyy_xxy_0[j] = pa_x[j] * tlx_xyy_xxy_0[j] + 0.5 * fl1_fx * tlx_yy_xxy_0[j] + fl1_fx * tlx_xyy_xy_0[j];

            tly_xxyy_xxy_0[j] = pa_x[j] * tly_xyy_xxy_0[j] + 0.5 * fl1_fx * tly_yy_xxy_0[j] + fl1_fx * tly_xyy_xy_0[j] +
                                0.5 * fl1_fx * tpz_xyy_xxy_0[j] + fl1_fx * fl1_fgb * tdz_xyy_xxy_0[j];

            tlz_xxyy_xxy_0[j] = pa_x[j] * tlz_xyy_xxy_0[j] + 0.5 * fl1_fx * tlz_yy_xxy_0[j] + fl1_fx * tlz_xyy_xy_0[j] -
                                0.5 * fl1_fx * tpy_xyy_xxy_0[j] - fl1_fx * fl1_fgb * tdy_xyy_xxy_0[j];

            tlx_xxyy_xxz_0[j] = pa_x[j] * tlx_xyy_xxz_0[j] + 0.5 * fl1_fx * tlx_yy_xxz_0[j] + fl1_fx * tlx_xyy_xz_0[j];

            tly_xxyy_xxz_0[j] = pa_x[j] * tly_xyy_xxz_0[j] + 0.5 * fl1_fx * tly_yy_xxz_0[j] + fl1_fx * tly_xyy_xz_0[j] +
                                0.5 * fl1_fx * tpz_xyy_xxz_0[j] + fl1_fx * fl1_fgb * tdz_xyy_xxz_0[j];

            tlz_xxyy_xxz_0[j] = pa_x[j] * tlz_xyy_xxz_0[j] + 0.5 * fl1_fx * tlz_yy_xxz_0[j] + fl1_fx * tlz_xyy_xz_0[j] -
                                0.5 * fl1_fx * tpy_xyy_xxz_0[j] - fl1_fx * fl1_fgb * tdy_xyy_xxz_0[j];

            tlx_xxyy_xyy_0[j] = pa_x[j] * tlx_xyy_xyy_0[j] + 0.5 * fl1_fx * tlx_yy_xyy_0[j] + 0.5 * fl1_fx * tlx_xyy_yy_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForGF_100_150(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_l_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_4_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tly_xyy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 21);

        auto tlz_xyy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 21);

        auto tlx_xyy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 22);

        auto tly_xyy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 22);

        auto tlz_xyy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 22);

        auto tlx_xyy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 23);

        auto tly_xyy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 23);

        auto tlz_xyy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 23);

        auto tlx_xyz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 24);

        auto tly_xyz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 24);

        auto tlz_xyz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 24);

        auto tlx_xyz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 25);

        auto tly_xyz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 25);

        auto tlz_xyz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 25);

        auto tlx_xyz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 26);

        auto tly_xyz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 26);

        auto tlz_xyz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 26);

        auto tlx_xyz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 27);

        auto tly_xyz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 27);

        auto tlz_xyz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 27);

        auto tlx_xyz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 28);

        auto tly_xyz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 28);

        auto tlz_xyz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 28);

        auto tlx_xyz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 29);

        auto tly_xyz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 29);

        auto tlz_xyz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 29);

        auto tpy_xyy_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 33);

        auto tpz_xyy_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 33);

        auto tpy_xyy_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 34);

        auto tpz_xyy_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 34);

        auto tpy_xyy_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 35);

        auto tpz_xyy_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 35);

        auto tpy_xyy_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 36);

        auto tpz_xyy_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 36);

        auto tpy_xyy_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 37);

        auto tpz_xyy_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 37);

        auto tpy_xyy_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 38);

        auto tpz_xyy_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 38);

        auto tpy_xyy_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 39);

        auto tpz_xyy_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 39);

        auto tpy_xyz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 40);

        auto tpz_xyz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 40);

        auto tpy_xyz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 41);

        auto tpz_xyz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 41);

        auto tpy_xyz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 42);

        auto tpz_xyz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 42);

        auto tpy_xyz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 43);

        auto tpz_xyz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 43);

        auto tpy_xyz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 44);

        auto tpz_xyz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 44);

        auto tpy_xyz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 45);

        auto tpz_xyz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 45);

        auto tpy_xyz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 46);

        auto tpz_xyz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 46);

        auto tpy_xyz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 47);

        auto tpz_xyz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 47);

        auto tpy_xyz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 48);

        auto tpz_xyz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 48);

        auto tpy_xyz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 49);

        auto tpz_xyz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 49);

        auto tdy_xyy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 33);

        auto tdz_xyy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 33);

        auto tdy_xyy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 34);

        auto tdz_xyy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 34);

        auto tdy_xyy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 35);

        auto tdz_xyy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 35);

        auto tdy_xyy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 36);

        auto tdz_xyy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 36);

        auto tdy_xyy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 37);

        auto tdz_xyy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 37);

        auto tdy_xyy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 38);

        auto tdz_xyy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 38);

        auto tdy_xyy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 39);

        auto tdz_xyy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 39);

        auto tdy_xyz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 40);

        auto tdz_xyz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 40);

        auto tdy_xyz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 41);

        auto tdz_xyz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 41);

        auto tdy_xyz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 42);

        auto tdz_xyz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 42);

        auto tdy_xyz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 43);

        auto tdz_xyz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 43);

        auto tdy_xyz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 44);

        auto tdz_xyz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 44);

        auto tdy_xyz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 45);

        auto tdz_xyz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 45);

        auto tdy_xyz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 46);

        auto tdz_xyz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 46);

        auto tdy_xyz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 47);

        auto tdz_xyz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 47);

        auto tdy_xyz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 48);

        auto tdz_xyz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 48);

        auto tdy_xyz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 49);

        auto tdz_xyz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 49);

        // set up pointers to integrals

        auto tly_xxyy_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 33);

        auto tlz_xxyy_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 33);

        auto tlx_xxyy_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 34);

        auto tly_xxyy_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 34);

        auto tlz_xxyy_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 34);

        auto tlx_xxyy_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 35);

        auto tly_xxyy_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 35);

        auto tlz_xxyy_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 35);

        auto tlx_xxyy_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 36);

        auto tly_xxyy_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 36);

        auto tlz_xxyy_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 36);

        auto tlx_xxyy_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 37);

        auto tly_xxyy_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 37);

        auto tlz_xxyy_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 37);

        auto tlx_xxyy_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 38);

        auto tly_xxyy_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 38);

        auto tlz_xxyy_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 38);

        auto tlx_xxyy_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 39);

        auto tly_xxyy_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 39);

        auto tlz_xxyy_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 39);

        auto tlx_xxyz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 40);

        auto tly_xxyz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 40);

        auto tlz_xxyz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 40);

        auto tlx_xxyz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 41);

        auto tly_xxyz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 41);

        auto tlz_xxyz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 41);

        auto tlx_xxyz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 42);

        auto tly_xxyz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 42);

        auto tlz_xxyz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 42);

        auto tlx_xxyz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 43);

        auto tly_xxyz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 43);

        auto tlz_xxyz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 43);

        auto tlx_xxyz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 44);

        auto tly_xxyz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 44);

        auto tlz_xxyz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 44);

        auto tlx_xxyz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 45);

        auto tly_xxyz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 45);

        auto tlz_xxyz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 45);

        auto tlx_xxyz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 46);

        auto tly_xxyz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 46);

        auto tlz_xxyz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 46);

        auto tlx_xxyz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 47);

        auto tly_xxyz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 47);

        auto tlz_xxyz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 47);

        auto tlx_xxyz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 48);

        auto tly_xxyz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 48);

        auto tlz_xxyz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 48);

        auto tlx_xxyz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 49);

        auto tly_xxyz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 49);

        auto tlz_xxyz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 49);

        // Batch of Integrals (100,150)

        #pragma omp simd aligned(fgb, fx, pa_x, tdy_xyy_xyy_0, tdy_xyy_xyz_0, tdy_xyy_xzz_0, \
                                     tdy_xyy_yyy_0, tdy_xyy_yyz_0, tdy_xyy_yzz_0, tdy_xyy_zzz_0, tdy_xyz_xxx_0, \
                                     tdy_xyz_xxy_0, tdy_xyz_xxz_0, tdy_xyz_xyy_0, tdy_xyz_xyz_0, tdy_xyz_xzz_0, \
                                     tdy_xyz_yyy_0, tdy_xyz_yyz_0, tdy_xyz_yzz_0, tdy_xyz_zzz_0, tdz_xyy_xyy_0, \
                                     tdz_xyy_xyz_0, tdz_xyy_xzz_0, tdz_xyy_yyy_0, tdz_xyy_yyz_0, tdz_xyy_yzz_0, \
                                     tdz_xyy_zzz_0, tdz_xyz_xxx_0, tdz_xyz_xxy_0, tdz_xyz_xxz_0, tdz_xyz_xyy_0, \
                                     tdz_xyz_xyz_0, tdz_xyz_xzz_0, tdz_xyz_yyy_0, tdz_xyz_yyz_0, tdz_xyz_yzz_0, \
                                     tdz_xyz_zzz_0, tlx_xxyy_xyz_0, tlx_xxyy_xzz_0, tlx_xxyy_yyy_0, tlx_xxyy_yyz_0, \
                                     tlx_xxyy_yzz_0, tlx_xxyy_zzz_0, tlx_xxyz_xxx_0, tlx_xxyz_xxy_0, tlx_xxyz_xxz_0, \
                                     tlx_xxyz_xyy_0, tlx_xxyz_xyz_0, tlx_xxyz_xzz_0, tlx_xxyz_yyy_0, tlx_xxyz_yyz_0, \
                                     tlx_xxyz_yzz_0, tlx_xxyz_zzz_0, tlx_xyy_xyz_0, tlx_xyy_xzz_0, tlx_xyy_yyy_0, \
                                     tlx_xyy_yyz_0, tlx_xyy_yz_0, tlx_xyy_yzz_0, tlx_xyy_zz_0, tlx_xyy_zzz_0, \
                                     tlx_xyz_xx_0, tlx_xyz_xxx_0, tlx_xyz_xxy_0, tlx_xyz_xxz_0, tlx_xyz_xy_0, \
                                     tlx_xyz_xyy_0, tlx_xyz_xyz_0, tlx_xyz_xz_0, tlx_xyz_xzz_0, tlx_xyz_yy_0, \
                                     tlx_xyz_yyy_0, tlx_xyz_yyz_0, tlx_xyz_yz_0, tlx_xyz_yzz_0, tlx_xyz_zz_0, \
                                     tlx_xyz_zzz_0, tlx_yy_xyz_0, tlx_yy_xzz_0, tlx_yy_yyy_0, tlx_yy_yyz_0, tlx_yy_yzz_0, \
                                     tlx_yy_zzz_0, tlx_yz_xxx_0, tlx_yz_xxy_0, tlx_yz_xxz_0, tlx_yz_xyy_0, tlx_yz_xyz_0, \
                                     tlx_yz_xzz_0, tlx_yz_yyy_0, tlx_yz_yyz_0, tlx_yz_yzz_0, tlx_yz_zzz_0, \
                                     tly_xxyy_xyy_0, tly_xxyy_xyz_0, tly_xxyy_xzz_0, tly_xxyy_yyy_0, tly_xxyy_yyz_0, \
                                     tly_xxyy_yzz_0, tly_xxyy_zzz_0, tly_xxyz_xxx_0, tly_xxyz_xxy_0, tly_xxyz_xxz_0, \
                                     tly_xxyz_xyy_0, tly_xxyz_xyz_0, tly_xxyz_xzz_0, tly_xxyz_yyy_0, tly_xxyz_yyz_0, \
                                     tly_xxyz_yzz_0, tly_xxyz_zzz_0, tly_xyy_xyy_0, tly_xyy_xyz_0, tly_xyy_xzz_0, \
                                     tly_xyy_yy_0, tly_xyy_yyy_0, tly_xyy_yyz_0, tly_xyy_yz_0, tly_xyy_yzz_0, \
                                     tly_xyy_zz_0, tly_xyy_zzz_0, tly_xyz_xx_0, tly_xyz_xxx_0, tly_xyz_xxy_0, \
                                     tly_xyz_xxz_0, tly_xyz_xy_0, tly_xyz_xyy_0, tly_xyz_xyz_0, tly_xyz_xz_0, \
                                     tly_xyz_xzz_0, tly_xyz_yy_0, tly_xyz_yyy_0, tly_xyz_yyz_0, tly_xyz_yz_0, \
                                     tly_xyz_yzz_0, tly_xyz_zz_0, tly_xyz_zzz_0, tly_yy_xyy_0, tly_yy_xyz_0, \
                                     tly_yy_xzz_0, tly_yy_yyy_0, tly_yy_yyz_0, tly_yy_yzz_0, tly_yy_zzz_0, tly_yz_xxx_0, \
                                     tly_yz_xxy_0, tly_yz_xxz_0, tly_yz_xyy_0, tly_yz_xyz_0, tly_yz_xzz_0, tly_yz_yyy_0, \
                                     tly_yz_yyz_0, tly_yz_yzz_0, tly_yz_zzz_0, tlz_xxyy_xyy_0, tlz_xxyy_xyz_0, \
                                     tlz_xxyy_xzz_0, tlz_xxyy_yyy_0, tlz_xxyy_yyz_0, tlz_xxyy_yzz_0, tlz_xxyy_zzz_0, \
                                     tlz_xxyz_xxx_0, tlz_xxyz_xxy_0, tlz_xxyz_xxz_0, tlz_xxyz_xyy_0, tlz_xxyz_xyz_0, \
                                     tlz_xxyz_xzz_0, tlz_xxyz_yyy_0, tlz_xxyz_yyz_0, tlz_xxyz_yzz_0, tlz_xxyz_zzz_0, \
                                     tlz_xyy_xyy_0, tlz_xyy_xyz_0, tlz_xyy_xzz_0, tlz_xyy_yy_0, tlz_xyy_yyy_0, \
                                     tlz_xyy_yyz_0, tlz_xyy_yz_0, tlz_xyy_yzz_0, tlz_xyy_zz_0, tlz_xyy_zzz_0, \
                                     tlz_xyz_xx_0, tlz_xyz_xxx_0, tlz_xyz_xxy_0, tlz_xyz_xxz_0, tlz_xyz_xy_0, \
                                     tlz_xyz_xyy_0, tlz_xyz_xyz_0, tlz_xyz_xz_0, tlz_xyz_xzz_0, tlz_xyz_yy_0, \
                                     tlz_xyz_yyy_0, tlz_xyz_yyz_0, tlz_xyz_yz_0, tlz_xyz_yzz_0, tlz_xyz_zz_0, \
                                     tlz_xyz_zzz_0, tlz_yy_xyy_0, tlz_yy_xyz_0, tlz_yy_xzz_0, tlz_yy_yyy_0, tlz_yy_yyz_0, \
                                     tlz_yy_yzz_0, tlz_yy_zzz_0, tlz_yz_xxx_0, tlz_yz_xxy_0, tlz_yz_xxz_0, tlz_yz_xyy_0, \
                                     tlz_yz_xyz_0, tlz_yz_xzz_0, tlz_yz_yyy_0, tlz_yz_yyz_0, tlz_yz_yzz_0, tlz_yz_zzz_0, \
                                     tpy_xyy_xyy_0, tpy_xyy_xyz_0, tpy_xyy_xzz_0, tpy_xyy_yyy_0, tpy_xyy_yyz_0, \
                                     tpy_xyy_yzz_0, tpy_xyy_zzz_0, tpy_xyz_xxx_0, tpy_xyz_xxy_0, tpy_xyz_xxz_0, \
                                     tpy_xyz_xyy_0, tpy_xyz_xyz_0, tpy_xyz_xzz_0, tpy_xyz_yyy_0, tpy_xyz_yyz_0, \
                                     tpy_xyz_yzz_0, tpy_xyz_zzz_0, tpz_xyy_xyy_0, tpz_xyy_xyz_0, tpz_xyy_xzz_0, \
                                     tpz_xyy_yyy_0, tpz_xyy_yyz_0, tpz_xyy_yzz_0, tpz_xyy_zzz_0, tpz_xyz_xxx_0, \
                                     tpz_xyz_xxy_0, tpz_xyz_xxz_0, tpz_xyz_xyy_0, tpz_xyz_xyz_0, tpz_xyz_xzz_0, \
                                     tpz_xyz_yyy_0, tpz_xyz_yyz_0, tpz_xyz_yzz_0, tpz_xyz_zzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tly_xxyy_xyy_0[j] = pa_x[j] * tly_xyy_xyy_0[j] + 0.5 * fl1_fx * tly_yy_xyy_0[j] + 0.5 * fl1_fx * tly_xyy_yy_0[j] +
                                0.5 * fl1_fx * tpz_xyy_xyy_0[j] + fl1_fx * fl1_fgb * tdz_xyy_xyy_0[j];

            tlz_xxyy_xyy_0[j] = pa_x[j] * tlz_xyy_xyy_0[j] + 0.5 * fl1_fx * tlz_yy_xyy_0[j] + 0.5 * fl1_fx * tlz_xyy_yy_0[j] -
                                0.5 * fl1_fx * tpy_xyy_xyy_0[j] - fl1_fx * fl1_fgb * tdy_xyy_xyy_0[j];

            tlx_xxyy_xyz_0[j] = pa_x[j] * tlx_xyy_xyz_0[j] + 0.5 * fl1_fx * tlx_yy_xyz_0[j] + 0.5 * fl1_fx * tlx_xyy_yz_0[j];

            tly_xxyy_xyz_0[j] = pa_x[j] * tly_xyy_xyz_0[j] + 0.5 * fl1_fx * tly_yy_xyz_0[j] + 0.5 * fl1_fx * tly_xyy_yz_0[j] +
                                0.5 * fl1_fx * tpz_xyy_xyz_0[j] + fl1_fx * fl1_fgb * tdz_xyy_xyz_0[j];

            tlz_xxyy_xyz_0[j] = pa_x[j] * tlz_xyy_xyz_0[j] + 0.5 * fl1_fx * tlz_yy_xyz_0[j] + 0.5 * fl1_fx * tlz_xyy_yz_0[j] -
                                0.5 * fl1_fx * tpy_xyy_xyz_0[j] - fl1_fx * fl1_fgb * tdy_xyy_xyz_0[j];

            tlx_xxyy_xzz_0[j] = pa_x[j] * tlx_xyy_xzz_0[j] + 0.5 * fl1_fx * tlx_yy_xzz_0[j] + 0.5 * fl1_fx * tlx_xyy_zz_0[j];

            tly_xxyy_xzz_0[j] = pa_x[j] * tly_xyy_xzz_0[j] + 0.5 * fl1_fx * tly_yy_xzz_0[j] + 0.5 * fl1_fx * tly_xyy_zz_0[j] +
                                0.5 * fl1_fx * tpz_xyy_xzz_0[j] + fl1_fx * fl1_fgb * tdz_xyy_xzz_0[j];

            tlz_xxyy_xzz_0[j] = pa_x[j] * tlz_xyy_xzz_0[j] + 0.5 * fl1_fx * tlz_yy_xzz_0[j] + 0.5 * fl1_fx * tlz_xyy_zz_0[j] -
                                0.5 * fl1_fx * tpy_xyy_xzz_0[j] - fl1_fx * fl1_fgb * tdy_xyy_xzz_0[j];

            tlx_xxyy_yyy_0[j] = pa_x[j] * tlx_xyy_yyy_0[j] + 0.5 * fl1_fx * tlx_yy_yyy_0[j];

            tly_xxyy_yyy_0[j] =
                pa_x[j] * tly_xyy_yyy_0[j] + 0.5 * fl1_fx * tly_yy_yyy_0[j] + 0.5 * fl1_fx * tpz_xyy_yyy_0[j] + fl1_fx * fl1_fgb * tdz_xyy_yyy_0[j];

            tlz_xxyy_yyy_0[j] =
                pa_x[j] * tlz_xyy_yyy_0[j] + 0.5 * fl1_fx * tlz_yy_yyy_0[j] - 0.5 * fl1_fx * tpy_xyy_yyy_0[j] - fl1_fx * fl1_fgb * tdy_xyy_yyy_0[j];

            tlx_xxyy_yyz_0[j] = pa_x[j] * tlx_xyy_yyz_0[j] + 0.5 * fl1_fx * tlx_yy_yyz_0[j];

            tly_xxyy_yyz_0[j] =
                pa_x[j] * tly_xyy_yyz_0[j] + 0.5 * fl1_fx * tly_yy_yyz_0[j] + 0.5 * fl1_fx * tpz_xyy_yyz_0[j] + fl1_fx * fl1_fgb * tdz_xyy_yyz_0[j];

            tlz_xxyy_yyz_0[j] =
                pa_x[j] * tlz_xyy_yyz_0[j] + 0.5 * fl1_fx * tlz_yy_yyz_0[j] - 0.5 * fl1_fx * tpy_xyy_yyz_0[j] - fl1_fx * fl1_fgb * tdy_xyy_yyz_0[j];

            tlx_xxyy_yzz_0[j] = pa_x[j] * tlx_xyy_yzz_0[j] + 0.5 * fl1_fx * tlx_yy_yzz_0[j];

            tly_xxyy_yzz_0[j] =
                pa_x[j] * tly_xyy_yzz_0[j] + 0.5 * fl1_fx * tly_yy_yzz_0[j] + 0.5 * fl1_fx * tpz_xyy_yzz_0[j] + fl1_fx * fl1_fgb * tdz_xyy_yzz_0[j];

            tlz_xxyy_yzz_0[j] =
                pa_x[j] * tlz_xyy_yzz_0[j] + 0.5 * fl1_fx * tlz_yy_yzz_0[j] - 0.5 * fl1_fx * tpy_xyy_yzz_0[j] - fl1_fx * fl1_fgb * tdy_xyy_yzz_0[j];

            tlx_xxyy_zzz_0[j] = pa_x[j] * tlx_xyy_zzz_0[j] + 0.5 * fl1_fx * tlx_yy_zzz_0[j];

            tly_xxyy_zzz_0[j] =
                pa_x[j] * tly_xyy_zzz_0[j] + 0.5 * fl1_fx * tly_yy_zzz_0[j] + 0.5 * fl1_fx * tpz_xyy_zzz_0[j] + fl1_fx * fl1_fgb * tdz_xyy_zzz_0[j];

            tlz_xxyy_zzz_0[j] =
                pa_x[j] * tlz_xyy_zzz_0[j] + 0.5 * fl1_fx * tlz_yy_zzz_0[j] - 0.5 * fl1_fx * tpy_xyy_zzz_0[j] - fl1_fx * fl1_fgb * tdy_xyy_zzz_0[j];

            tlx_xxyz_xxx_0[j] = pa_x[j] * tlx_xyz_xxx_0[j] + 0.5 * fl1_fx * tlx_yz_xxx_0[j] + 1.5 * fl1_fx * tlx_xyz_xx_0[j];

            tly_xxyz_xxx_0[j] = pa_x[j] * tly_xyz_xxx_0[j] + 0.5 * fl1_fx * tly_yz_xxx_0[j] + 1.5 * fl1_fx * tly_xyz_xx_0[j] +
                                0.5 * fl1_fx * tpz_xyz_xxx_0[j] + fl1_fx * fl1_fgb * tdz_xyz_xxx_0[j];

            tlz_xxyz_xxx_0[j] = pa_x[j] * tlz_xyz_xxx_0[j] + 0.5 * fl1_fx * tlz_yz_xxx_0[j] + 1.5 * fl1_fx * tlz_xyz_xx_0[j] -
                                0.5 * fl1_fx * tpy_xyz_xxx_0[j] - fl1_fx * fl1_fgb * tdy_xyz_xxx_0[j];

            tlx_xxyz_xxy_0[j] = pa_x[j] * tlx_xyz_xxy_0[j] + 0.5 * fl1_fx * tlx_yz_xxy_0[j] + fl1_fx * tlx_xyz_xy_0[j];

            tly_xxyz_xxy_0[j] = pa_x[j] * tly_xyz_xxy_0[j] + 0.5 * fl1_fx * tly_yz_xxy_0[j] + fl1_fx * tly_xyz_xy_0[j] +
                                0.5 * fl1_fx * tpz_xyz_xxy_0[j] + fl1_fx * fl1_fgb * tdz_xyz_xxy_0[j];

            tlz_xxyz_xxy_0[j] = pa_x[j] * tlz_xyz_xxy_0[j] + 0.5 * fl1_fx * tlz_yz_xxy_0[j] + fl1_fx * tlz_xyz_xy_0[j] -
                                0.5 * fl1_fx * tpy_xyz_xxy_0[j] - fl1_fx * fl1_fgb * tdy_xyz_xxy_0[j];

            tlx_xxyz_xxz_0[j] = pa_x[j] * tlx_xyz_xxz_0[j] + 0.5 * fl1_fx * tlx_yz_xxz_0[j] + fl1_fx * tlx_xyz_xz_0[j];

            tly_xxyz_xxz_0[j] = pa_x[j] * tly_xyz_xxz_0[j] + 0.5 * fl1_fx * tly_yz_xxz_0[j] + fl1_fx * tly_xyz_xz_0[j] +
                                0.5 * fl1_fx * tpz_xyz_xxz_0[j] + fl1_fx * fl1_fgb * tdz_xyz_xxz_0[j];

            tlz_xxyz_xxz_0[j] = pa_x[j] * tlz_xyz_xxz_0[j] + 0.5 * fl1_fx * tlz_yz_xxz_0[j] + fl1_fx * tlz_xyz_xz_0[j] -
                                0.5 * fl1_fx * tpy_xyz_xxz_0[j] - fl1_fx * fl1_fgb * tdy_xyz_xxz_0[j];

            tlx_xxyz_xyy_0[j] = pa_x[j] * tlx_xyz_xyy_0[j] + 0.5 * fl1_fx * tlx_yz_xyy_0[j] + 0.5 * fl1_fx * tlx_xyz_yy_0[j];

            tly_xxyz_xyy_0[j] = pa_x[j] * tly_xyz_xyy_0[j] + 0.5 * fl1_fx * tly_yz_xyy_0[j] + 0.5 * fl1_fx * tly_xyz_yy_0[j] +
                                0.5 * fl1_fx * tpz_xyz_xyy_0[j] + fl1_fx * fl1_fgb * tdz_xyz_xyy_0[j];

            tlz_xxyz_xyy_0[j] = pa_x[j] * tlz_xyz_xyy_0[j] + 0.5 * fl1_fx * tlz_yz_xyy_0[j] + 0.5 * fl1_fx * tlz_xyz_yy_0[j] -
                                0.5 * fl1_fx * tpy_xyz_xyy_0[j] - fl1_fx * fl1_fgb * tdy_xyz_xyy_0[j];

            tlx_xxyz_xyz_0[j] = pa_x[j] * tlx_xyz_xyz_0[j] + 0.5 * fl1_fx * tlx_yz_xyz_0[j] + 0.5 * fl1_fx * tlx_xyz_yz_0[j];

            tly_xxyz_xyz_0[j] = pa_x[j] * tly_xyz_xyz_0[j] + 0.5 * fl1_fx * tly_yz_xyz_0[j] + 0.5 * fl1_fx * tly_xyz_yz_0[j] +
                                0.5 * fl1_fx * tpz_xyz_xyz_0[j] + fl1_fx * fl1_fgb * tdz_xyz_xyz_0[j];

            tlz_xxyz_xyz_0[j] = pa_x[j] * tlz_xyz_xyz_0[j] + 0.5 * fl1_fx * tlz_yz_xyz_0[j] + 0.5 * fl1_fx * tlz_xyz_yz_0[j] -
                                0.5 * fl1_fx * tpy_xyz_xyz_0[j] - fl1_fx * fl1_fgb * tdy_xyz_xyz_0[j];

            tlx_xxyz_xzz_0[j] = pa_x[j] * tlx_xyz_xzz_0[j] + 0.5 * fl1_fx * tlx_yz_xzz_0[j] + 0.5 * fl1_fx * tlx_xyz_zz_0[j];

            tly_xxyz_xzz_0[j] = pa_x[j] * tly_xyz_xzz_0[j] + 0.5 * fl1_fx * tly_yz_xzz_0[j] + 0.5 * fl1_fx * tly_xyz_zz_0[j] +
                                0.5 * fl1_fx * tpz_xyz_xzz_0[j] + fl1_fx * fl1_fgb * tdz_xyz_xzz_0[j];

            tlz_xxyz_xzz_0[j] = pa_x[j] * tlz_xyz_xzz_0[j] + 0.5 * fl1_fx * tlz_yz_xzz_0[j] + 0.5 * fl1_fx * tlz_xyz_zz_0[j] -
                                0.5 * fl1_fx * tpy_xyz_xzz_0[j] - fl1_fx * fl1_fgb * tdy_xyz_xzz_0[j];

            tlx_xxyz_yyy_0[j] = pa_x[j] * tlx_xyz_yyy_0[j] + 0.5 * fl1_fx * tlx_yz_yyy_0[j];

            tly_xxyz_yyy_0[j] =
                pa_x[j] * tly_xyz_yyy_0[j] + 0.5 * fl1_fx * tly_yz_yyy_0[j] + 0.5 * fl1_fx * tpz_xyz_yyy_0[j] + fl1_fx * fl1_fgb * tdz_xyz_yyy_0[j];

            tlz_xxyz_yyy_0[j] =
                pa_x[j] * tlz_xyz_yyy_0[j] + 0.5 * fl1_fx * tlz_yz_yyy_0[j] - 0.5 * fl1_fx * tpy_xyz_yyy_0[j] - fl1_fx * fl1_fgb * tdy_xyz_yyy_0[j];

            tlx_xxyz_yyz_0[j] = pa_x[j] * tlx_xyz_yyz_0[j] + 0.5 * fl1_fx * tlx_yz_yyz_0[j];

            tly_xxyz_yyz_0[j] =
                pa_x[j] * tly_xyz_yyz_0[j] + 0.5 * fl1_fx * tly_yz_yyz_0[j] + 0.5 * fl1_fx * tpz_xyz_yyz_0[j] + fl1_fx * fl1_fgb * tdz_xyz_yyz_0[j];

            tlz_xxyz_yyz_0[j] =
                pa_x[j] * tlz_xyz_yyz_0[j] + 0.5 * fl1_fx * tlz_yz_yyz_0[j] - 0.5 * fl1_fx * tpy_xyz_yyz_0[j] - fl1_fx * fl1_fgb * tdy_xyz_yyz_0[j];

            tlx_xxyz_yzz_0[j] = pa_x[j] * tlx_xyz_yzz_0[j] + 0.5 * fl1_fx * tlx_yz_yzz_0[j];

            tly_xxyz_yzz_0[j] =
                pa_x[j] * tly_xyz_yzz_0[j] + 0.5 * fl1_fx * tly_yz_yzz_0[j] + 0.5 * fl1_fx * tpz_xyz_yzz_0[j] + fl1_fx * fl1_fgb * tdz_xyz_yzz_0[j];

            tlz_xxyz_yzz_0[j] =
                pa_x[j] * tlz_xyz_yzz_0[j] + 0.5 * fl1_fx * tlz_yz_yzz_0[j] - 0.5 * fl1_fx * tpy_xyz_yzz_0[j] - fl1_fx * fl1_fgb * tdy_xyz_yzz_0[j];

            tlx_xxyz_zzz_0[j] = pa_x[j] * tlx_xyz_zzz_0[j] + 0.5 * fl1_fx * tlx_yz_zzz_0[j];

            tly_xxyz_zzz_0[j] =
                pa_x[j] * tly_xyz_zzz_0[j] + 0.5 * fl1_fx * tly_yz_zzz_0[j] + 0.5 * fl1_fx * tpz_xyz_zzz_0[j] + fl1_fx * fl1_fgb * tdz_xyz_zzz_0[j];

            tlz_xxyz_zzz_0[j] =
                pa_x[j] * tlz_xyz_zzz_0[j] + 0.5 * fl1_fx * tlz_yz_zzz_0[j] - 0.5 * fl1_fx * tpy_xyz_zzz_0[j] - fl1_fx * fl1_fgb * tdy_xyz_zzz_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForGF_150_200(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_l_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_4_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tlx_xzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 30);

        auto tly_xzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 30);

        auto tlz_xzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 30);

        auto tlx_xzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 31);

        auto tly_xzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 31);

        auto tlz_xzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 31);

        auto tlx_xzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 32);

        auto tly_xzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 32);

        auto tlz_xzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 32);

        auto tlx_xzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 33);

        auto tly_xzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 33);

        auto tlz_xzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 33);

        auto tlx_xzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 34);

        auto tly_xzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 34);

        auto tlz_xzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 34);

        auto tlx_xzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 35);

        auto tly_xzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 35);

        auto tlz_xzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 35);

        auto tlx_yyy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 36);

        auto tly_yyy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 36);

        auto tlz_yyy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 36);

        auto tlx_yyy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 37);

        auto tly_yyy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 37);

        auto tlz_yyy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 37);

        auto tlx_yyy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 38);

        auto tly_yyy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 38);

        auto tlz_yyy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 38);

        auto tlx_yyy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 39);

        auto tly_yyy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 39);

        auto tlz_yyy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 39);

        auto tlx_yyy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 40);

        auto tly_yyy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 40);

        auto tlz_yyy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 40);

        auto tlx_yyy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 41);

        auto tly_yyy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 41);

        auto tlz_yyy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 41);

        auto tpy_xzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 50);

        auto tpz_xzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 50);

        auto tpy_xzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 51);

        auto tpz_xzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 51);

        auto tpy_xzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 52);

        auto tpz_xzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 52);

        auto tpy_xzz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 53);

        auto tpz_xzz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 53);

        auto tpy_xzz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 54);

        auto tpz_xzz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 54);

        auto tpy_xzz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 55);

        auto tpz_xzz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 55);

        auto tpy_xzz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 56);

        auto tpz_xzz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 56);

        auto tpy_xzz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 57);

        auto tpz_xzz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 57);

        auto tpy_xzz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 58);

        auto tpz_xzz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 58);

        auto tpy_xzz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 59);

        auto tpz_xzz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 59);

        auto tpy_yyy_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 60);

        auto tpz_yyy_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 60);

        auto tpy_yyy_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 61);

        auto tpz_yyy_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 61);

        auto tpy_yyy_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 62);

        auto tpz_yyy_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 62);

        auto tpy_yyy_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 63);

        auto tpz_yyy_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 63);

        auto tpy_yyy_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 64);

        auto tpz_yyy_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 64);

        auto tpy_yyy_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 65);

        auto tpz_yyy_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 65);

        auto tpz_yyy_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 66);

        auto tdy_xzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 50);

        auto tdz_xzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 50);

        auto tdy_xzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 51);

        auto tdz_xzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 51);

        auto tdy_xzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 52);

        auto tdz_xzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 52);

        auto tdy_xzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 53);

        auto tdz_xzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 53);

        auto tdy_xzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 54);

        auto tdz_xzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 54);

        auto tdy_xzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 55);

        auto tdz_xzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 55);

        auto tdy_xzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 56);

        auto tdz_xzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 56);

        auto tdy_xzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 57);

        auto tdz_xzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 57);

        auto tdy_xzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 58);

        auto tdz_xzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 58);

        auto tdy_xzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 59);

        auto tdz_xzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 59);

        auto tdy_yyy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 60);

        auto tdz_yyy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 60);

        auto tdy_yyy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 61);

        auto tdz_yyy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 61);

        auto tdy_yyy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 62);

        auto tdz_yyy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 62);

        auto tdy_yyy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 63);

        auto tdz_yyy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 63);

        auto tdy_yyy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 64);

        auto tdz_yyy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 64);

        auto tdy_yyy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 65);

        auto tdz_yyy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 65);

        auto tdz_yyy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 66);

        // set up pointers to integrals

        auto tlx_xxzz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 50);

        auto tly_xxzz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 50);

        auto tlz_xxzz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 50);

        auto tlx_xxzz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 51);

        auto tly_xxzz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 51);

        auto tlz_xxzz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 51);

        auto tlx_xxzz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 52);

        auto tly_xxzz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 52);

        auto tlz_xxzz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 52);

        auto tlx_xxzz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 53);

        auto tly_xxzz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 53);

        auto tlz_xxzz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 53);

        auto tlx_xxzz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 54);

        auto tly_xxzz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 54);

        auto tlz_xxzz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 54);

        auto tlx_xxzz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 55);

        auto tly_xxzz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 55);

        auto tlz_xxzz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 55);

        auto tlx_xxzz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 56);

        auto tly_xxzz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 56);

        auto tlz_xxzz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 56);

        auto tlx_xxzz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 57);

        auto tly_xxzz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 57);

        auto tlz_xxzz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 57);

        auto tlx_xxzz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 58);

        auto tly_xxzz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 58);

        auto tlz_xxzz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 58);

        auto tlx_xxzz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 59);

        auto tly_xxzz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 59);

        auto tlz_xxzz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 59);

        auto tlx_xyyy_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 60);

        auto tly_xyyy_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 60);

        auto tlz_xyyy_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 60);

        auto tlx_xyyy_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 61);

        auto tly_xyyy_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 61);

        auto tlz_xyyy_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 61);

        auto tlx_xyyy_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 62);

        auto tly_xyyy_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 62);

        auto tlz_xyyy_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 62);

        auto tlx_xyyy_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 63);

        auto tly_xyyy_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 63);

        auto tlz_xyyy_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 63);

        auto tlx_xyyy_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 64);

        auto tly_xyyy_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 64);

        auto tlz_xyyy_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 64);

        auto tlx_xyyy_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 65);

        auto tly_xyyy_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 65);

        auto tlz_xyyy_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 65);

        auto tlx_xyyy_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 66);

        auto tly_xyyy_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 66);

        // Batch of Integrals (150,200)

        #pragma omp simd aligned(fgb, fx, pa_x, tdy_xzz_xxx_0, tdy_xzz_xxy_0, tdy_xzz_xxz_0, \
                                     tdy_xzz_xyy_0, tdy_xzz_xyz_0, tdy_xzz_xzz_0, tdy_xzz_yyy_0, tdy_xzz_yyz_0, \
                                     tdy_xzz_yzz_0, tdy_xzz_zzz_0, tdy_yyy_xxx_0, tdy_yyy_xxy_0, tdy_yyy_xxz_0, \
                                     tdy_yyy_xyy_0, tdy_yyy_xyz_0, tdy_yyy_xzz_0, tdz_xzz_xxx_0, tdz_xzz_xxy_0, \
                                     tdz_xzz_xxz_0, tdz_xzz_xyy_0, tdz_xzz_xyz_0, tdz_xzz_xzz_0, tdz_xzz_yyy_0, \
                                     tdz_xzz_yyz_0, tdz_xzz_yzz_0, tdz_xzz_zzz_0, tdz_yyy_xxx_0, tdz_yyy_xxy_0, \
                                     tdz_yyy_xxz_0, tdz_yyy_xyy_0, tdz_yyy_xyz_0, tdz_yyy_xzz_0, tdz_yyy_yyy_0, \
                                     tlx_xxzz_xxx_0, tlx_xxzz_xxy_0, tlx_xxzz_xxz_0, tlx_xxzz_xyy_0, tlx_xxzz_xyz_0, \
                                     tlx_xxzz_xzz_0, tlx_xxzz_yyy_0, tlx_xxzz_yyz_0, tlx_xxzz_yzz_0, tlx_xxzz_zzz_0, \
                                     tlx_xyyy_xxx_0, tlx_xyyy_xxy_0, tlx_xyyy_xxz_0, tlx_xyyy_xyy_0, tlx_xyyy_xyz_0, \
                                     tlx_xyyy_xzz_0, tlx_xyyy_yyy_0, tlx_xzz_xx_0, tlx_xzz_xxx_0, tlx_xzz_xxy_0, \
                                     tlx_xzz_xxz_0, tlx_xzz_xy_0, tlx_xzz_xyy_0, tlx_xzz_xyz_0, tlx_xzz_xz_0, \
                                     tlx_xzz_xzz_0, tlx_xzz_yy_0, tlx_xzz_yyy_0, tlx_xzz_yyz_0, tlx_xzz_yz_0, \
                                     tlx_xzz_yzz_0, tlx_xzz_zz_0, tlx_xzz_zzz_0, tlx_yyy_xx_0, tlx_yyy_xxx_0, \
                                     tlx_yyy_xxy_0, tlx_yyy_xxz_0, tlx_yyy_xy_0, tlx_yyy_xyy_0, tlx_yyy_xyz_0, \
                                     tlx_yyy_xz_0, tlx_yyy_xzz_0, tlx_yyy_yy_0, tlx_yyy_yyy_0, tlx_yyy_yz_0, \
                                     tlx_yyy_zz_0, tlx_zz_xxx_0, tlx_zz_xxy_0, tlx_zz_xxz_0, tlx_zz_xyy_0, tlx_zz_xyz_0, \
                                     tlx_zz_xzz_0, tlx_zz_yyy_0, tlx_zz_yyz_0, tlx_zz_yzz_0, tlx_zz_zzz_0, \
                                     tly_xxzz_xxx_0, tly_xxzz_xxy_0, tly_xxzz_xxz_0, tly_xxzz_xyy_0, tly_xxzz_xyz_0, \
                                     tly_xxzz_xzz_0, tly_xxzz_yyy_0, tly_xxzz_yyz_0, tly_xxzz_yzz_0, tly_xxzz_zzz_0, \
                                     tly_xyyy_xxx_0, tly_xyyy_xxy_0, tly_xyyy_xxz_0, tly_xyyy_xyy_0, tly_xyyy_xyz_0, \
                                     tly_xyyy_xzz_0, tly_xyyy_yyy_0, tly_xzz_xx_0, tly_xzz_xxx_0, tly_xzz_xxy_0, \
                                     tly_xzz_xxz_0, tly_xzz_xy_0, tly_xzz_xyy_0, tly_xzz_xyz_0, tly_xzz_xz_0, \
                                     tly_xzz_xzz_0, tly_xzz_yy_0, tly_xzz_yyy_0, tly_xzz_yyz_0, tly_xzz_yz_0, \
                                     tly_xzz_yzz_0, tly_xzz_zz_0, tly_xzz_zzz_0, tly_yyy_xx_0, tly_yyy_xxx_0, \
                                     tly_yyy_xxy_0, tly_yyy_xxz_0, tly_yyy_xy_0, tly_yyy_xyy_0, tly_yyy_xyz_0, \
                                     tly_yyy_xz_0, tly_yyy_xzz_0, tly_yyy_yy_0, tly_yyy_yyy_0, tly_yyy_yz_0, \
                                     tly_yyy_zz_0, tly_zz_xxx_0, tly_zz_xxy_0, tly_zz_xxz_0, tly_zz_xyy_0, tly_zz_xyz_0, \
                                     tly_zz_xzz_0, tly_zz_yyy_0, tly_zz_yyz_0, tly_zz_yzz_0, tly_zz_zzz_0, \
                                     tlz_xxzz_xxx_0, tlz_xxzz_xxy_0, tlz_xxzz_xxz_0, tlz_xxzz_xyy_0, tlz_xxzz_xyz_0, \
                                     tlz_xxzz_xzz_0, tlz_xxzz_yyy_0, tlz_xxzz_yyz_0, tlz_xxzz_yzz_0, tlz_xxzz_zzz_0, \
                                     tlz_xyyy_xxx_0, tlz_xyyy_xxy_0, tlz_xyyy_xxz_0, tlz_xyyy_xyy_0, tlz_xyyy_xyz_0, \
                                     tlz_xyyy_xzz_0, tlz_xzz_xx_0, tlz_xzz_xxx_0, tlz_xzz_xxy_0, tlz_xzz_xxz_0, \
                                     tlz_xzz_xy_0, tlz_xzz_xyy_0, tlz_xzz_xyz_0, tlz_xzz_xz_0, tlz_xzz_xzz_0, \
                                     tlz_xzz_yy_0, tlz_xzz_yyy_0, tlz_xzz_yyz_0, tlz_xzz_yz_0, tlz_xzz_yzz_0, \
                                     tlz_xzz_zz_0, tlz_xzz_zzz_0, tlz_yyy_xx_0, tlz_yyy_xxx_0, tlz_yyy_xxy_0, \
                                     tlz_yyy_xxz_0, tlz_yyy_xy_0, tlz_yyy_xyy_0, tlz_yyy_xyz_0, tlz_yyy_xz_0, \
                                     tlz_yyy_xzz_0, tlz_yyy_yy_0, tlz_yyy_yz_0, tlz_yyy_zz_0, tlz_zz_xxx_0, tlz_zz_xxy_0, \
                                     tlz_zz_xxz_0, tlz_zz_xyy_0, tlz_zz_xyz_0, tlz_zz_xzz_0, tlz_zz_yyy_0, tlz_zz_yyz_0, \
                                     tlz_zz_yzz_0, tlz_zz_zzz_0, tpy_xzz_xxx_0, tpy_xzz_xxy_0, tpy_xzz_xxz_0, \
                                     tpy_xzz_xyy_0, tpy_xzz_xyz_0, tpy_xzz_xzz_0, tpy_xzz_yyy_0, tpy_xzz_yyz_0, \
                                     tpy_xzz_yzz_0, tpy_xzz_zzz_0, tpy_yyy_xxx_0, tpy_yyy_xxy_0, tpy_yyy_xxz_0, \
                                     tpy_yyy_xyy_0, tpy_yyy_xyz_0, tpy_yyy_xzz_0, tpz_xzz_xxx_0, tpz_xzz_xxy_0, \
                                     tpz_xzz_xxz_0, tpz_xzz_xyy_0, tpz_xzz_xyz_0, tpz_xzz_xzz_0, tpz_xzz_yyy_0, \
                                     tpz_xzz_yyz_0, tpz_xzz_yzz_0, tpz_xzz_zzz_0, tpz_yyy_xxx_0, tpz_yyy_xxy_0, \
                                     tpz_yyy_xxz_0, tpz_yyy_xyy_0, tpz_yyy_xyz_0, tpz_yyy_xzz_0, tpz_yyy_yyy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_xxzz_xxx_0[j] = pa_x[j] * tlx_xzz_xxx_0[j] + 0.5 * fl1_fx * tlx_zz_xxx_0[j] + 1.5 * fl1_fx * tlx_xzz_xx_0[j];

            tly_xxzz_xxx_0[j] = pa_x[j] * tly_xzz_xxx_0[j] + 0.5 * fl1_fx * tly_zz_xxx_0[j] + 1.5 * fl1_fx * tly_xzz_xx_0[j] +
                                0.5 * fl1_fx * tpz_xzz_xxx_0[j] + fl1_fx * fl1_fgb * tdz_xzz_xxx_0[j];

            tlz_xxzz_xxx_0[j] = pa_x[j] * tlz_xzz_xxx_0[j] + 0.5 * fl1_fx * tlz_zz_xxx_0[j] + 1.5 * fl1_fx * tlz_xzz_xx_0[j] -
                                0.5 * fl1_fx * tpy_xzz_xxx_0[j] - fl1_fx * fl1_fgb * tdy_xzz_xxx_0[j];

            tlx_xxzz_xxy_0[j] = pa_x[j] * tlx_xzz_xxy_0[j] + 0.5 * fl1_fx * tlx_zz_xxy_0[j] + fl1_fx * tlx_xzz_xy_0[j];

            tly_xxzz_xxy_0[j] = pa_x[j] * tly_xzz_xxy_0[j] + 0.5 * fl1_fx * tly_zz_xxy_0[j] + fl1_fx * tly_xzz_xy_0[j] +
                                0.5 * fl1_fx * tpz_xzz_xxy_0[j] + fl1_fx * fl1_fgb * tdz_xzz_xxy_0[j];

            tlz_xxzz_xxy_0[j] = pa_x[j] * tlz_xzz_xxy_0[j] + 0.5 * fl1_fx * tlz_zz_xxy_0[j] + fl1_fx * tlz_xzz_xy_0[j] -
                                0.5 * fl1_fx * tpy_xzz_xxy_0[j] - fl1_fx * fl1_fgb * tdy_xzz_xxy_0[j];

            tlx_xxzz_xxz_0[j] = pa_x[j] * tlx_xzz_xxz_0[j] + 0.5 * fl1_fx * tlx_zz_xxz_0[j] + fl1_fx * tlx_xzz_xz_0[j];

            tly_xxzz_xxz_0[j] = pa_x[j] * tly_xzz_xxz_0[j] + 0.5 * fl1_fx * tly_zz_xxz_0[j] + fl1_fx * tly_xzz_xz_0[j] +
                                0.5 * fl1_fx * tpz_xzz_xxz_0[j] + fl1_fx * fl1_fgb * tdz_xzz_xxz_0[j];

            tlz_xxzz_xxz_0[j] = pa_x[j] * tlz_xzz_xxz_0[j] + 0.5 * fl1_fx * tlz_zz_xxz_0[j] + fl1_fx * tlz_xzz_xz_0[j] -
                                0.5 * fl1_fx * tpy_xzz_xxz_0[j] - fl1_fx * fl1_fgb * tdy_xzz_xxz_0[j];

            tlx_xxzz_xyy_0[j] = pa_x[j] * tlx_xzz_xyy_0[j] + 0.5 * fl1_fx * tlx_zz_xyy_0[j] + 0.5 * fl1_fx * tlx_xzz_yy_0[j];

            tly_xxzz_xyy_0[j] = pa_x[j] * tly_xzz_xyy_0[j] + 0.5 * fl1_fx * tly_zz_xyy_0[j] + 0.5 * fl1_fx * tly_xzz_yy_0[j] +
                                0.5 * fl1_fx * tpz_xzz_xyy_0[j] + fl1_fx * fl1_fgb * tdz_xzz_xyy_0[j];

            tlz_xxzz_xyy_0[j] = pa_x[j] * tlz_xzz_xyy_0[j] + 0.5 * fl1_fx * tlz_zz_xyy_0[j] + 0.5 * fl1_fx * tlz_xzz_yy_0[j] -
                                0.5 * fl1_fx * tpy_xzz_xyy_0[j] - fl1_fx * fl1_fgb * tdy_xzz_xyy_0[j];

            tlx_xxzz_xyz_0[j] = pa_x[j] * tlx_xzz_xyz_0[j] + 0.5 * fl1_fx * tlx_zz_xyz_0[j] + 0.5 * fl1_fx * tlx_xzz_yz_0[j];

            tly_xxzz_xyz_0[j] = pa_x[j] * tly_xzz_xyz_0[j] + 0.5 * fl1_fx * tly_zz_xyz_0[j] + 0.5 * fl1_fx * tly_xzz_yz_0[j] +
                                0.5 * fl1_fx * tpz_xzz_xyz_0[j] + fl1_fx * fl1_fgb * tdz_xzz_xyz_0[j];

            tlz_xxzz_xyz_0[j] = pa_x[j] * tlz_xzz_xyz_0[j] + 0.5 * fl1_fx * tlz_zz_xyz_0[j] + 0.5 * fl1_fx * tlz_xzz_yz_0[j] -
                                0.5 * fl1_fx * tpy_xzz_xyz_0[j] - fl1_fx * fl1_fgb * tdy_xzz_xyz_0[j];

            tlx_xxzz_xzz_0[j] = pa_x[j] * tlx_xzz_xzz_0[j] + 0.5 * fl1_fx * tlx_zz_xzz_0[j] + 0.5 * fl1_fx * tlx_xzz_zz_0[j];

            tly_xxzz_xzz_0[j] = pa_x[j] * tly_xzz_xzz_0[j] + 0.5 * fl1_fx * tly_zz_xzz_0[j] + 0.5 * fl1_fx * tly_xzz_zz_0[j] +
                                0.5 * fl1_fx * tpz_xzz_xzz_0[j] + fl1_fx * fl1_fgb * tdz_xzz_xzz_0[j];

            tlz_xxzz_xzz_0[j] = pa_x[j] * tlz_xzz_xzz_0[j] + 0.5 * fl1_fx * tlz_zz_xzz_0[j] + 0.5 * fl1_fx * tlz_xzz_zz_0[j] -
                                0.5 * fl1_fx * tpy_xzz_xzz_0[j] - fl1_fx * fl1_fgb * tdy_xzz_xzz_0[j];

            tlx_xxzz_yyy_0[j] = pa_x[j] * tlx_xzz_yyy_0[j] + 0.5 * fl1_fx * tlx_zz_yyy_0[j];

            tly_xxzz_yyy_0[j] =
                pa_x[j] * tly_xzz_yyy_0[j] + 0.5 * fl1_fx * tly_zz_yyy_0[j] + 0.5 * fl1_fx * tpz_xzz_yyy_0[j] + fl1_fx * fl1_fgb * tdz_xzz_yyy_0[j];

            tlz_xxzz_yyy_0[j] =
                pa_x[j] * tlz_xzz_yyy_0[j] + 0.5 * fl1_fx * tlz_zz_yyy_0[j] - 0.5 * fl1_fx * tpy_xzz_yyy_0[j] - fl1_fx * fl1_fgb * tdy_xzz_yyy_0[j];

            tlx_xxzz_yyz_0[j] = pa_x[j] * tlx_xzz_yyz_0[j] + 0.5 * fl1_fx * tlx_zz_yyz_0[j];

            tly_xxzz_yyz_0[j] =
                pa_x[j] * tly_xzz_yyz_0[j] + 0.5 * fl1_fx * tly_zz_yyz_0[j] + 0.5 * fl1_fx * tpz_xzz_yyz_0[j] + fl1_fx * fl1_fgb * tdz_xzz_yyz_0[j];

            tlz_xxzz_yyz_0[j] =
                pa_x[j] * tlz_xzz_yyz_0[j] + 0.5 * fl1_fx * tlz_zz_yyz_0[j] - 0.5 * fl1_fx * tpy_xzz_yyz_0[j] - fl1_fx * fl1_fgb * tdy_xzz_yyz_0[j];

            tlx_xxzz_yzz_0[j] = pa_x[j] * tlx_xzz_yzz_0[j] + 0.5 * fl1_fx * tlx_zz_yzz_0[j];

            tly_xxzz_yzz_0[j] =
                pa_x[j] * tly_xzz_yzz_0[j] + 0.5 * fl1_fx * tly_zz_yzz_0[j] + 0.5 * fl1_fx * tpz_xzz_yzz_0[j] + fl1_fx * fl1_fgb * tdz_xzz_yzz_0[j];

            tlz_xxzz_yzz_0[j] =
                pa_x[j] * tlz_xzz_yzz_0[j] + 0.5 * fl1_fx * tlz_zz_yzz_0[j] - 0.5 * fl1_fx * tpy_xzz_yzz_0[j] - fl1_fx * fl1_fgb * tdy_xzz_yzz_0[j];

            tlx_xxzz_zzz_0[j] = pa_x[j] * tlx_xzz_zzz_0[j] + 0.5 * fl1_fx * tlx_zz_zzz_0[j];

            tly_xxzz_zzz_0[j] =
                pa_x[j] * tly_xzz_zzz_0[j] + 0.5 * fl1_fx * tly_zz_zzz_0[j] + 0.5 * fl1_fx * tpz_xzz_zzz_0[j] + fl1_fx * fl1_fgb * tdz_xzz_zzz_0[j];

            tlz_xxzz_zzz_0[j] =
                pa_x[j] * tlz_xzz_zzz_0[j] + 0.5 * fl1_fx * tlz_zz_zzz_0[j] - 0.5 * fl1_fx * tpy_xzz_zzz_0[j] - fl1_fx * fl1_fgb * tdy_xzz_zzz_0[j];

            tlx_xyyy_xxx_0[j] = pa_x[j] * tlx_yyy_xxx_0[j] + 1.5 * fl1_fx * tlx_yyy_xx_0[j];

            tly_xyyy_xxx_0[j] =
                pa_x[j] * tly_yyy_xxx_0[j] + 1.5 * fl1_fx * tly_yyy_xx_0[j] + 0.5 * fl1_fx * tpz_yyy_xxx_0[j] + fl1_fx * fl1_fgb * tdz_yyy_xxx_0[j];

            tlz_xyyy_xxx_0[j] =
                pa_x[j] * tlz_yyy_xxx_0[j] + 1.5 * fl1_fx * tlz_yyy_xx_0[j] - 0.5 * fl1_fx * tpy_yyy_xxx_0[j] - fl1_fx * fl1_fgb * tdy_yyy_xxx_0[j];

            tlx_xyyy_xxy_0[j] = pa_x[j] * tlx_yyy_xxy_0[j] + fl1_fx * tlx_yyy_xy_0[j];

            tly_xyyy_xxy_0[j] =
                pa_x[j] * tly_yyy_xxy_0[j] + fl1_fx * tly_yyy_xy_0[j] + 0.5 * fl1_fx * tpz_yyy_xxy_0[j] + fl1_fx * fl1_fgb * tdz_yyy_xxy_0[j];

            tlz_xyyy_xxy_0[j] =
                pa_x[j] * tlz_yyy_xxy_0[j] + fl1_fx * tlz_yyy_xy_0[j] - 0.5 * fl1_fx * tpy_yyy_xxy_0[j] - fl1_fx * fl1_fgb * tdy_yyy_xxy_0[j];

            tlx_xyyy_xxz_0[j] = pa_x[j] * tlx_yyy_xxz_0[j] + fl1_fx * tlx_yyy_xz_0[j];

            tly_xyyy_xxz_0[j] =
                pa_x[j] * tly_yyy_xxz_0[j] + fl1_fx * tly_yyy_xz_0[j] + 0.5 * fl1_fx * tpz_yyy_xxz_0[j] + fl1_fx * fl1_fgb * tdz_yyy_xxz_0[j];

            tlz_xyyy_xxz_0[j] =
                pa_x[j] * tlz_yyy_xxz_0[j] + fl1_fx * tlz_yyy_xz_0[j] - 0.5 * fl1_fx * tpy_yyy_xxz_0[j] - fl1_fx * fl1_fgb * tdy_yyy_xxz_0[j];

            tlx_xyyy_xyy_0[j] = pa_x[j] * tlx_yyy_xyy_0[j] + 0.5 * fl1_fx * tlx_yyy_yy_0[j];

            tly_xyyy_xyy_0[j] =
                pa_x[j] * tly_yyy_xyy_0[j] + 0.5 * fl1_fx * tly_yyy_yy_0[j] + 0.5 * fl1_fx * tpz_yyy_xyy_0[j] + fl1_fx * fl1_fgb * tdz_yyy_xyy_0[j];

            tlz_xyyy_xyy_0[j] =
                pa_x[j] * tlz_yyy_xyy_0[j] + 0.5 * fl1_fx * tlz_yyy_yy_0[j] - 0.5 * fl1_fx * tpy_yyy_xyy_0[j] - fl1_fx * fl1_fgb * tdy_yyy_xyy_0[j];

            tlx_xyyy_xyz_0[j] = pa_x[j] * tlx_yyy_xyz_0[j] + 0.5 * fl1_fx * tlx_yyy_yz_0[j];

            tly_xyyy_xyz_0[j] =
                pa_x[j] * tly_yyy_xyz_0[j] + 0.5 * fl1_fx * tly_yyy_yz_0[j] + 0.5 * fl1_fx * tpz_yyy_xyz_0[j] + fl1_fx * fl1_fgb * tdz_yyy_xyz_0[j];

            tlz_xyyy_xyz_0[j] =
                pa_x[j] * tlz_yyy_xyz_0[j] + 0.5 * fl1_fx * tlz_yyy_yz_0[j] - 0.5 * fl1_fx * tpy_yyy_xyz_0[j] - fl1_fx * fl1_fgb * tdy_yyy_xyz_0[j];

            tlx_xyyy_xzz_0[j] = pa_x[j] * tlx_yyy_xzz_0[j] + 0.5 * fl1_fx * tlx_yyy_zz_0[j];

            tly_xyyy_xzz_0[j] =
                pa_x[j] * tly_yyy_xzz_0[j] + 0.5 * fl1_fx * tly_yyy_zz_0[j] + 0.5 * fl1_fx * tpz_yyy_xzz_0[j] + fl1_fx * fl1_fgb * tdz_yyy_xzz_0[j];

            tlz_xyyy_xzz_0[j] =
                pa_x[j] * tlz_yyy_xzz_0[j] + 0.5 * fl1_fx * tlz_yyy_zz_0[j] - 0.5 * fl1_fx * tpy_yyy_xzz_0[j] - fl1_fx * fl1_fgb * tdy_yyy_xzz_0[j];

            tlx_xyyy_yyy_0[j] = pa_x[j] * tlx_yyy_yyy_0[j];

            tly_xyyy_yyy_0[j] = pa_x[j] * tly_yyy_yyy_0[j] + 0.5 * fl1_fx * tpz_yyy_yyy_0[j] + fl1_fx * fl1_fgb * tdz_yyy_yyy_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForGF_200_250(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_l_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_4_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

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

        auto tlx_yyz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 42);

        auto tly_yyz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 42);

        auto tlz_yyz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 42);

        auto tlx_yyz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 43);

        auto tly_yyz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 43);

        auto tlz_yyz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 43);

        auto tlx_yyz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 44);

        auto tly_yyz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 44);

        auto tlz_yyz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 44);

        auto tlx_yyz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 45);

        auto tly_yyz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 45);

        auto tlz_yyz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 45);

        auto tlx_yyz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 46);

        auto tly_yyz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 46);

        auto tlz_yyz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 46);

        auto tlx_yyz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 47);

        auto tly_yyz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 47);

        auto tlz_yyz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 47);

        auto tlx_yzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 48);

        auto tly_yzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 48);

        auto tlz_yzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 48);

        auto tlx_yzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 49);

        auto tly_yzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 49);

        auto tlz_yzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 49);

        auto tlx_yzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 50);

        auto tly_yzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 50);

        auto tlz_yzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 50);

        auto tlx_yzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 51);

        auto tpy_yyy_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 66);

        auto tpy_yyy_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 67);

        auto tpz_yyy_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 67);

        auto tpy_yyy_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 68);

        auto tpz_yyy_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 68);

        auto tpy_yyy_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 69);

        auto tpz_yyy_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 69);

        auto tpy_yyz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 70);

        auto tpz_yyz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 70);

        auto tpy_yyz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 71);

        auto tpz_yyz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 71);

        auto tpy_yyz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 72);

        auto tpz_yyz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 72);

        auto tpy_yyz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 73);

        auto tpz_yyz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 73);

        auto tpy_yyz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 74);

        auto tpz_yyz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 74);

        auto tpy_yyz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 75);

        auto tpz_yyz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 75);

        auto tpy_yyz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 76);

        auto tpz_yyz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 76);

        auto tpy_yyz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 77);

        auto tpz_yyz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 77);

        auto tpy_yyz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 78);

        auto tpz_yyz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 78);

        auto tpy_yyz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 79);

        auto tpz_yyz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 79);

        auto tpy_yzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 80);

        auto tpz_yzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 80);

        auto tpy_yzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 81);

        auto tpz_yzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 81);

        auto tpy_yzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 82);

        auto tpz_yzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 82);

        auto tdy_yyy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 66);

        auto tdy_yyy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 67);

        auto tdz_yyy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 67);

        auto tdy_yyy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 68);

        auto tdz_yyy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 68);

        auto tdy_yyy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 69);

        auto tdz_yyy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 69);

        auto tdy_yyz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 70);

        auto tdz_yyz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 70);

        auto tdy_yyz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 71);

        auto tdz_yyz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 71);

        auto tdy_yyz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 72);

        auto tdz_yyz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 72);

        auto tdy_yyz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 73);

        auto tdz_yyz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 73);

        auto tdy_yyz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 74);

        auto tdz_yyz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 74);

        auto tdy_yyz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 75);

        auto tdz_yyz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 75);

        auto tdy_yyz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 76);

        auto tdz_yyz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 76);

        auto tdy_yyz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 77);

        auto tdz_yyz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 77);

        auto tdy_yyz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 78);

        auto tdz_yyz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 78);

        auto tdy_yyz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 79);

        auto tdz_yyz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 79);

        auto tdy_yzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 80);

        auto tdz_yzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 80);

        auto tdy_yzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 81);

        auto tdz_yzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 81);

        auto tdy_yzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 82);

        auto tdz_yzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 82);

        // set up pointers to integrals

        auto tlz_xyyy_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 66);

        auto tlx_xyyy_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 67);

        auto tly_xyyy_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 67);

        auto tlz_xyyy_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 67);

        auto tlx_xyyy_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 68);

        auto tly_xyyy_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 68);

        auto tlz_xyyy_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 68);

        auto tlx_xyyy_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 69);

        auto tly_xyyy_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 69);

        auto tlz_xyyy_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 69);

        auto tlx_xyyz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 70);

        auto tly_xyyz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 70);

        auto tlz_xyyz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 70);

        auto tlx_xyyz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 71);

        auto tly_xyyz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 71);

        auto tlz_xyyz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 71);

        auto tlx_xyyz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 72);

        auto tly_xyyz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 72);

        auto tlz_xyyz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 72);

        auto tlx_xyyz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 73);

        auto tly_xyyz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 73);

        auto tlz_xyyz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 73);

        auto tlx_xyyz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 74);

        auto tly_xyyz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 74);

        auto tlz_xyyz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 74);

        auto tlx_xyyz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 75);

        auto tly_xyyz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 75);

        auto tlz_xyyz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 75);

        auto tlx_xyyz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 76);

        auto tly_xyyz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 76);

        auto tlz_xyyz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 76);

        auto tlx_xyyz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 77);

        auto tly_xyyz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 77);

        auto tlz_xyyz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 77);

        auto tlx_xyyz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 78);

        auto tly_xyyz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 78);

        auto tlz_xyyz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 78);

        auto tlx_xyyz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 79);

        auto tly_xyyz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 79);

        auto tlz_xyyz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 79);

        auto tlx_xyzz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 80);

        auto tly_xyzz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 80);

        auto tlz_xyzz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 80);

        auto tlx_xyzz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 81);

        auto tly_xyzz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 81);

        auto tlz_xyzz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 81);

        auto tlx_xyzz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 82);

        auto tly_xyzz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 82);

        auto tlz_xyzz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 82);

        auto tlx_xyzz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 83);

        // Batch of Integrals (200,250)

        #pragma omp simd aligned(fgb, fx, pa_x, tdy_yyy_yyy_0, tdy_yyy_yyz_0, tdy_yyy_yzz_0, \
                                     tdy_yyy_zzz_0, tdy_yyz_xxx_0, tdy_yyz_xxy_0, tdy_yyz_xxz_0, tdy_yyz_xyy_0, \
                                     tdy_yyz_xyz_0, tdy_yyz_xzz_0, tdy_yyz_yyy_0, tdy_yyz_yyz_0, tdy_yyz_yzz_0, \
                                     tdy_yyz_zzz_0, tdy_yzz_xxx_0, tdy_yzz_xxy_0, tdy_yzz_xxz_0, tdz_yyy_yyz_0, \
                                     tdz_yyy_yzz_0, tdz_yyy_zzz_0, tdz_yyz_xxx_0, tdz_yyz_xxy_0, tdz_yyz_xxz_0, \
                                     tdz_yyz_xyy_0, tdz_yyz_xyz_0, tdz_yyz_xzz_0, tdz_yyz_yyy_0, tdz_yyz_yyz_0, \
                                     tdz_yyz_yzz_0, tdz_yyz_zzz_0, tdz_yzz_xxx_0, tdz_yzz_xxy_0, tdz_yzz_xxz_0, \
                                     tlx_xyyy_yyz_0, tlx_xyyy_yzz_0, tlx_xyyy_zzz_0, tlx_xyyz_xxx_0, tlx_xyyz_xxy_0, \
                                     tlx_xyyz_xxz_0, tlx_xyyz_xyy_0, tlx_xyyz_xyz_0, tlx_xyyz_xzz_0, tlx_xyyz_yyy_0, \
                                     tlx_xyyz_yyz_0, tlx_xyyz_yzz_0, tlx_xyyz_zzz_0, tlx_xyzz_xxx_0, tlx_xyzz_xxy_0, \
                                     tlx_xyzz_xxz_0, tlx_xyzz_xyy_0, tlx_yyy_yyz_0, tlx_yyy_yzz_0, tlx_yyy_zzz_0, \
                                     tlx_yyz_xx_0, tlx_yyz_xxx_0, tlx_yyz_xxy_0, tlx_yyz_xxz_0, tlx_yyz_xy_0, \
                                     tlx_yyz_xyy_0, tlx_yyz_xyz_0, tlx_yyz_xz_0, tlx_yyz_xzz_0, tlx_yyz_yy_0, \
                                     tlx_yyz_yyy_0, tlx_yyz_yyz_0, tlx_yyz_yz_0, tlx_yyz_yzz_0, tlx_yyz_zz_0, \
                                     tlx_yyz_zzz_0, tlx_yzz_xx_0, tlx_yzz_xxx_0, tlx_yzz_xxy_0, tlx_yzz_xxz_0, \
                                     tlx_yzz_xy_0, tlx_yzz_xyy_0, tlx_yzz_xz_0, tlx_yzz_yy_0, tly_xyyy_yyz_0, \
                                     tly_xyyy_yzz_0, tly_xyyy_zzz_0, tly_xyyz_xxx_0, tly_xyyz_xxy_0, tly_xyyz_xxz_0, \
                                     tly_xyyz_xyy_0, tly_xyyz_xyz_0, tly_xyyz_xzz_0, tly_xyyz_yyy_0, tly_xyyz_yyz_0, \
                                     tly_xyyz_yzz_0, tly_xyyz_zzz_0, tly_xyzz_xxx_0, tly_xyzz_xxy_0, tly_xyzz_xxz_0, \
                                     tly_yyy_yyz_0, tly_yyy_yzz_0, tly_yyy_zzz_0, tly_yyz_xx_0, tly_yyz_xxx_0, \
                                     tly_yyz_xxy_0, tly_yyz_xxz_0, tly_yyz_xy_0, tly_yyz_xyy_0, tly_yyz_xyz_0, \
                                     tly_yyz_xz_0, tly_yyz_xzz_0, tly_yyz_yy_0, tly_yyz_yyy_0, tly_yyz_yyz_0, \
                                     tly_yyz_yz_0, tly_yyz_yzz_0, tly_yyz_zz_0, tly_yyz_zzz_0, tly_yzz_xx_0, \
                                     tly_yzz_xxx_0, tly_yzz_xxy_0, tly_yzz_xxz_0, tly_yzz_xy_0, tly_yzz_xz_0, \
                                     tlz_xyyy_yyy_0, tlz_xyyy_yyz_0, tlz_xyyy_yzz_0, tlz_xyyy_zzz_0, tlz_xyyz_xxx_0, \
                                     tlz_xyyz_xxy_0, tlz_xyyz_xxz_0, tlz_xyyz_xyy_0, tlz_xyyz_xyz_0, tlz_xyyz_xzz_0, \
                                     tlz_xyyz_yyy_0, tlz_xyyz_yyz_0, tlz_xyyz_yzz_0, tlz_xyyz_zzz_0, tlz_xyzz_xxx_0, \
                                     tlz_xyzz_xxy_0, tlz_xyzz_xxz_0, tlz_yyy_yyy_0, tlz_yyy_yyz_0, tlz_yyy_yzz_0, \
                                     tlz_yyy_zzz_0, tlz_yyz_xx_0, tlz_yyz_xxx_0, tlz_yyz_xxy_0, tlz_yyz_xxz_0, \
                                     tlz_yyz_xy_0, tlz_yyz_xyy_0, tlz_yyz_xyz_0, tlz_yyz_xz_0, tlz_yyz_xzz_0, \
                                     tlz_yyz_yy_0, tlz_yyz_yyy_0, tlz_yyz_yyz_0, tlz_yyz_yz_0, tlz_yyz_yzz_0, \
                                     tlz_yyz_zz_0, tlz_yyz_zzz_0, tlz_yzz_xx_0, tlz_yzz_xxx_0, tlz_yzz_xxy_0, \
                                     tlz_yzz_xxz_0, tlz_yzz_xy_0, tlz_yzz_xz_0, tpy_yyy_yyy_0, tpy_yyy_yyz_0, \
                                     tpy_yyy_yzz_0, tpy_yyy_zzz_0, tpy_yyz_xxx_0, tpy_yyz_xxy_0, tpy_yyz_xxz_0, \
                                     tpy_yyz_xyy_0, tpy_yyz_xyz_0, tpy_yyz_xzz_0, tpy_yyz_yyy_0, tpy_yyz_yyz_0, \
                                     tpy_yyz_yzz_0, tpy_yyz_zzz_0, tpy_yzz_xxx_0, tpy_yzz_xxy_0, tpy_yzz_xxz_0, \
                                     tpz_yyy_yyz_0, tpz_yyy_yzz_0, tpz_yyy_zzz_0, tpz_yyz_xxx_0, tpz_yyz_xxy_0, \
                                     tpz_yyz_xxz_0, tpz_yyz_xyy_0, tpz_yyz_xyz_0, tpz_yyz_xzz_0, tpz_yyz_yyy_0, \
                                     tpz_yyz_yyz_0, tpz_yyz_yzz_0, tpz_yyz_zzz_0, tpz_yzz_xxx_0, tpz_yzz_xxy_0, \
                                     tpz_yzz_xxz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlz_xyyy_yyy_0[j] = pa_x[j] * tlz_yyy_yyy_0[j] - 0.5 * fl1_fx * tpy_yyy_yyy_0[j] - fl1_fx * fl1_fgb * tdy_yyy_yyy_0[j];

            tlx_xyyy_yyz_0[j] = pa_x[j] * tlx_yyy_yyz_0[j];

            tly_xyyy_yyz_0[j] = pa_x[j] * tly_yyy_yyz_0[j] + 0.5 * fl1_fx * tpz_yyy_yyz_0[j] + fl1_fx * fl1_fgb * tdz_yyy_yyz_0[j];

            tlz_xyyy_yyz_0[j] = pa_x[j] * tlz_yyy_yyz_0[j] - 0.5 * fl1_fx * tpy_yyy_yyz_0[j] - fl1_fx * fl1_fgb * tdy_yyy_yyz_0[j];

            tlx_xyyy_yzz_0[j] = pa_x[j] * tlx_yyy_yzz_0[j];

            tly_xyyy_yzz_0[j] = pa_x[j] * tly_yyy_yzz_0[j] + 0.5 * fl1_fx * tpz_yyy_yzz_0[j] + fl1_fx * fl1_fgb * tdz_yyy_yzz_0[j];

            tlz_xyyy_yzz_0[j] = pa_x[j] * tlz_yyy_yzz_0[j] - 0.5 * fl1_fx * tpy_yyy_yzz_0[j] - fl1_fx * fl1_fgb * tdy_yyy_yzz_0[j];

            tlx_xyyy_zzz_0[j] = pa_x[j] * tlx_yyy_zzz_0[j];

            tly_xyyy_zzz_0[j] = pa_x[j] * tly_yyy_zzz_0[j] + 0.5 * fl1_fx * tpz_yyy_zzz_0[j] + fl1_fx * fl1_fgb * tdz_yyy_zzz_0[j];

            tlz_xyyy_zzz_0[j] = pa_x[j] * tlz_yyy_zzz_0[j] - 0.5 * fl1_fx * tpy_yyy_zzz_0[j] - fl1_fx * fl1_fgb * tdy_yyy_zzz_0[j];

            tlx_xyyz_xxx_0[j] = pa_x[j] * tlx_yyz_xxx_0[j] + 1.5 * fl1_fx * tlx_yyz_xx_0[j];

            tly_xyyz_xxx_0[j] =
                pa_x[j] * tly_yyz_xxx_0[j] + 1.5 * fl1_fx * tly_yyz_xx_0[j] + 0.5 * fl1_fx * tpz_yyz_xxx_0[j] + fl1_fx * fl1_fgb * tdz_yyz_xxx_0[j];

            tlz_xyyz_xxx_0[j] =
                pa_x[j] * tlz_yyz_xxx_0[j] + 1.5 * fl1_fx * tlz_yyz_xx_0[j] - 0.5 * fl1_fx * tpy_yyz_xxx_0[j] - fl1_fx * fl1_fgb * tdy_yyz_xxx_0[j];

            tlx_xyyz_xxy_0[j] = pa_x[j] * tlx_yyz_xxy_0[j] + fl1_fx * tlx_yyz_xy_0[j];

            tly_xyyz_xxy_0[j] =
                pa_x[j] * tly_yyz_xxy_0[j] + fl1_fx * tly_yyz_xy_0[j] + 0.5 * fl1_fx * tpz_yyz_xxy_0[j] + fl1_fx * fl1_fgb * tdz_yyz_xxy_0[j];

            tlz_xyyz_xxy_0[j] =
                pa_x[j] * tlz_yyz_xxy_0[j] + fl1_fx * tlz_yyz_xy_0[j] - 0.5 * fl1_fx * tpy_yyz_xxy_0[j] - fl1_fx * fl1_fgb * tdy_yyz_xxy_0[j];

            tlx_xyyz_xxz_0[j] = pa_x[j] * tlx_yyz_xxz_0[j] + fl1_fx * tlx_yyz_xz_0[j];

            tly_xyyz_xxz_0[j] =
                pa_x[j] * tly_yyz_xxz_0[j] + fl1_fx * tly_yyz_xz_0[j] + 0.5 * fl1_fx * tpz_yyz_xxz_0[j] + fl1_fx * fl1_fgb * tdz_yyz_xxz_0[j];

            tlz_xyyz_xxz_0[j] =
                pa_x[j] * tlz_yyz_xxz_0[j] + fl1_fx * tlz_yyz_xz_0[j] - 0.5 * fl1_fx * tpy_yyz_xxz_0[j] - fl1_fx * fl1_fgb * tdy_yyz_xxz_0[j];

            tlx_xyyz_xyy_0[j] = pa_x[j] * tlx_yyz_xyy_0[j] + 0.5 * fl1_fx * tlx_yyz_yy_0[j];

            tly_xyyz_xyy_0[j] =
                pa_x[j] * tly_yyz_xyy_0[j] + 0.5 * fl1_fx * tly_yyz_yy_0[j] + 0.5 * fl1_fx * tpz_yyz_xyy_0[j] + fl1_fx * fl1_fgb * tdz_yyz_xyy_0[j];

            tlz_xyyz_xyy_0[j] =
                pa_x[j] * tlz_yyz_xyy_0[j] + 0.5 * fl1_fx * tlz_yyz_yy_0[j] - 0.5 * fl1_fx * tpy_yyz_xyy_0[j] - fl1_fx * fl1_fgb * tdy_yyz_xyy_0[j];

            tlx_xyyz_xyz_0[j] = pa_x[j] * tlx_yyz_xyz_0[j] + 0.5 * fl1_fx * tlx_yyz_yz_0[j];

            tly_xyyz_xyz_0[j] =
                pa_x[j] * tly_yyz_xyz_0[j] + 0.5 * fl1_fx * tly_yyz_yz_0[j] + 0.5 * fl1_fx * tpz_yyz_xyz_0[j] + fl1_fx * fl1_fgb * tdz_yyz_xyz_0[j];

            tlz_xyyz_xyz_0[j] =
                pa_x[j] * tlz_yyz_xyz_0[j] + 0.5 * fl1_fx * tlz_yyz_yz_0[j] - 0.5 * fl1_fx * tpy_yyz_xyz_0[j] - fl1_fx * fl1_fgb * tdy_yyz_xyz_0[j];

            tlx_xyyz_xzz_0[j] = pa_x[j] * tlx_yyz_xzz_0[j] + 0.5 * fl1_fx * tlx_yyz_zz_0[j];

            tly_xyyz_xzz_0[j] =
                pa_x[j] * tly_yyz_xzz_0[j] + 0.5 * fl1_fx * tly_yyz_zz_0[j] + 0.5 * fl1_fx * tpz_yyz_xzz_0[j] + fl1_fx * fl1_fgb * tdz_yyz_xzz_0[j];

            tlz_xyyz_xzz_0[j] =
                pa_x[j] * tlz_yyz_xzz_0[j] + 0.5 * fl1_fx * tlz_yyz_zz_0[j] - 0.5 * fl1_fx * tpy_yyz_xzz_0[j] - fl1_fx * fl1_fgb * tdy_yyz_xzz_0[j];

            tlx_xyyz_yyy_0[j] = pa_x[j] * tlx_yyz_yyy_0[j];

            tly_xyyz_yyy_0[j] = pa_x[j] * tly_yyz_yyy_0[j] + 0.5 * fl1_fx * tpz_yyz_yyy_0[j] + fl1_fx * fl1_fgb * tdz_yyz_yyy_0[j];

            tlz_xyyz_yyy_0[j] = pa_x[j] * tlz_yyz_yyy_0[j] - 0.5 * fl1_fx * tpy_yyz_yyy_0[j] - fl1_fx * fl1_fgb * tdy_yyz_yyy_0[j];

            tlx_xyyz_yyz_0[j] = pa_x[j] * tlx_yyz_yyz_0[j];

            tly_xyyz_yyz_0[j] = pa_x[j] * tly_yyz_yyz_0[j] + 0.5 * fl1_fx * tpz_yyz_yyz_0[j] + fl1_fx * fl1_fgb * tdz_yyz_yyz_0[j];

            tlz_xyyz_yyz_0[j] = pa_x[j] * tlz_yyz_yyz_0[j] - 0.5 * fl1_fx * tpy_yyz_yyz_0[j] - fl1_fx * fl1_fgb * tdy_yyz_yyz_0[j];

            tlx_xyyz_yzz_0[j] = pa_x[j] * tlx_yyz_yzz_0[j];

            tly_xyyz_yzz_0[j] = pa_x[j] * tly_yyz_yzz_0[j] + 0.5 * fl1_fx * tpz_yyz_yzz_0[j] + fl1_fx * fl1_fgb * tdz_yyz_yzz_0[j];

            tlz_xyyz_yzz_0[j] = pa_x[j] * tlz_yyz_yzz_0[j] - 0.5 * fl1_fx * tpy_yyz_yzz_0[j] - fl1_fx * fl1_fgb * tdy_yyz_yzz_0[j];

            tlx_xyyz_zzz_0[j] = pa_x[j] * tlx_yyz_zzz_0[j];

            tly_xyyz_zzz_0[j] = pa_x[j] * tly_yyz_zzz_0[j] + 0.5 * fl1_fx * tpz_yyz_zzz_0[j] + fl1_fx * fl1_fgb * tdz_yyz_zzz_0[j];

            tlz_xyyz_zzz_0[j] = pa_x[j] * tlz_yyz_zzz_0[j] - 0.5 * fl1_fx * tpy_yyz_zzz_0[j] - fl1_fx * fl1_fgb * tdy_yyz_zzz_0[j];

            tlx_xyzz_xxx_0[j] = pa_x[j] * tlx_yzz_xxx_0[j] + 1.5 * fl1_fx * tlx_yzz_xx_0[j];

            tly_xyzz_xxx_0[j] =
                pa_x[j] * tly_yzz_xxx_0[j] + 1.5 * fl1_fx * tly_yzz_xx_0[j] + 0.5 * fl1_fx * tpz_yzz_xxx_0[j] + fl1_fx * fl1_fgb * tdz_yzz_xxx_0[j];

            tlz_xyzz_xxx_0[j] =
                pa_x[j] * tlz_yzz_xxx_0[j] + 1.5 * fl1_fx * tlz_yzz_xx_0[j] - 0.5 * fl1_fx * tpy_yzz_xxx_0[j] - fl1_fx * fl1_fgb * tdy_yzz_xxx_0[j];

            tlx_xyzz_xxy_0[j] = pa_x[j] * tlx_yzz_xxy_0[j] + fl1_fx * tlx_yzz_xy_0[j];

            tly_xyzz_xxy_0[j] =
                pa_x[j] * tly_yzz_xxy_0[j] + fl1_fx * tly_yzz_xy_0[j] + 0.5 * fl1_fx * tpz_yzz_xxy_0[j] + fl1_fx * fl1_fgb * tdz_yzz_xxy_0[j];

            tlz_xyzz_xxy_0[j] =
                pa_x[j] * tlz_yzz_xxy_0[j] + fl1_fx * tlz_yzz_xy_0[j] - 0.5 * fl1_fx * tpy_yzz_xxy_0[j] - fl1_fx * fl1_fgb * tdy_yzz_xxy_0[j];

            tlx_xyzz_xxz_0[j] = pa_x[j] * tlx_yzz_xxz_0[j] + fl1_fx * tlx_yzz_xz_0[j];

            tly_xyzz_xxz_0[j] =
                pa_x[j] * tly_yzz_xxz_0[j] + fl1_fx * tly_yzz_xz_0[j] + 0.5 * fl1_fx * tpz_yzz_xxz_0[j] + fl1_fx * fl1_fgb * tdz_yzz_xxz_0[j];

            tlz_xyzz_xxz_0[j] =
                pa_x[j] * tlz_yzz_xxz_0[j] + fl1_fx * tlz_yzz_xz_0[j] - 0.5 * fl1_fx * tpy_yzz_xxz_0[j] - fl1_fx * fl1_fgb * tdy_yzz_xxz_0[j];

            tlx_xyzz_xyy_0[j] = pa_x[j] * tlx_yzz_xyy_0[j] + 0.5 * fl1_fx * tlx_yzz_yy_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForGF_250_300(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_l_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_4_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tly_yzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 51);

        auto tlz_yzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 51);

        auto tlx_yzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 52);

        auto tly_yzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 52);

        auto tlz_yzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 52);

        auto tlx_yzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 53);

        auto tly_yzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 53);

        auto tlz_yzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 53);

        auto tlx_zzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 54);

        auto tly_zzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 54);

        auto tlz_zzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 54);

        auto tlx_zzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 55);

        auto tly_zzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 55);

        auto tlz_zzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 55);

        auto tlx_zzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 56);

        auto tly_zzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 56);

        auto tlz_zzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 56);

        auto tlx_zzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 57);

        auto tly_zzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 57);

        auto tlz_zzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 57);

        auto tlx_zzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 58);

        auto tly_zzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 58);

        auto tlz_zzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 58);

        auto tlx_zzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 59);

        auto tly_zzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 59);

        auto tlz_zzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 59);

        auto tpy_yzz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 83);

        auto tpz_yzz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 83);

        auto tpy_yzz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 84);

        auto tpz_yzz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 84);

        auto tpy_yzz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 85);

        auto tpz_yzz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 85);

        auto tpy_yzz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 86);

        auto tpz_yzz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 86);

        auto tpy_yzz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 87);

        auto tpz_yzz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 87);

        auto tpy_yzz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 88);

        auto tpz_yzz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 88);

        auto tpy_yzz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 89);

        auto tpz_yzz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 89);

        auto tpy_zzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 90);

        auto tpz_zzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 90);

        auto tpy_zzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 91);

        auto tpz_zzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 91);

        auto tpy_zzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 92);

        auto tpz_zzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 92);

        auto tpy_zzz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 93);

        auto tpz_zzz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 93);

        auto tpy_zzz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 94);

        auto tpz_zzz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 94);

        auto tpy_zzz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 95);

        auto tpz_zzz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 95);

        auto tpy_zzz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 96);

        auto tpz_zzz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 96);

        auto tpy_zzz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 97);

        auto tpz_zzz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 97);

        auto tpy_zzz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 98);

        auto tpz_zzz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 98);

        auto tpy_zzz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 99);

        auto tpz_zzz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 99);

        auto tdy_yzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 83);

        auto tdz_yzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 83);

        auto tdy_yzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 84);

        auto tdz_yzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 84);

        auto tdy_yzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 85);

        auto tdz_yzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 85);

        auto tdy_yzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 86);

        auto tdz_yzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 86);

        auto tdy_yzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 87);

        auto tdz_yzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 87);

        auto tdy_yzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 88);

        auto tdz_yzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 88);

        auto tdy_yzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 89);

        auto tdz_yzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 89);

        auto tdy_zzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 90);

        auto tdz_zzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 90);

        auto tdy_zzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 91);

        auto tdz_zzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 91);

        auto tdy_zzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 92);

        auto tdz_zzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 92);

        auto tdy_zzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 93);

        auto tdz_zzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 93);

        auto tdy_zzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 94);

        auto tdz_zzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 94);

        auto tdy_zzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 95);

        auto tdz_zzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 95);

        auto tdy_zzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 96);

        auto tdz_zzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 96);

        auto tdy_zzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 97);

        auto tdz_zzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 97);

        auto tdy_zzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 98);

        auto tdz_zzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 98);

        auto tdy_zzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 99);

        auto tdz_zzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 99);

        // set up pointers to integrals

        auto tly_xyzz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 83);

        auto tlz_xyzz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 83);

        auto tlx_xyzz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 84);

        auto tly_xyzz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 84);

        auto tlz_xyzz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 84);

        auto tlx_xyzz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 85);

        auto tly_xyzz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 85);

        auto tlz_xyzz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 85);

        auto tlx_xyzz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 86);

        auto tly_xyzz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 86);

        auto tlz_xyzz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 86);

        auto tlx_xyzz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 87);

        auto tly_xyzz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 87);

        auto tlz_xyzz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 87);

        auto tlx_xyzz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 88);

        auto tly_xyzz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 88);

        auto tlz_xyzz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 88);

        auto tlx_xyzz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 89);

        auto tly_xyzz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 89);

        auto tlz_xyzz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 89);

        auto tlx_xzzz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 90);

        auto tly_xzzz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 90);

        auto tlz_xzzz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 90);

        auto tlx_xzzz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 91);

        auto tly_xzzz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 91);

        auto tlz_xzzz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 91);

        auto tlx_xzzz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 92);

        auto tly_xzzz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 92);

        auto tlz_xzzz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 92);

        auto tlx_xzzz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 93);

        auto tly_xzzz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 93);

        auto tlz_xzzz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 93);

        auto tlx_xzzz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 94);

        auto tly_xzzz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 94);

        auto tlz_xzzz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 94);

        auto tlx_xzzz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 95);

        auto tly_xzzz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 95);

        auto tlz_xzzz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 95);

        auto tlx_xzzz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 96);

        auto tly_xzzz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 96);

        auto tlz_xzzz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 96);

        auto tlx_xzzz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 97);

        auto tly_xzzz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 97);

        auto tlz_xzzz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 97);

        auto tlx_xzzz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 98);

        auto tly_xzzz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 98);

        auto tlz_xzzz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 98);

        auto tlx_xzzz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 99);

        auto tly_xzzz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 99);

        auto tlz_xzzz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 99);

        // Batch of Integrals (250,300)

        #pragma omp simd aligned(fgb, fx, pa_x, tdy_yzz_xyy_0, tdy_yzz_xyz_0, tdy_yzz_xzz_0, \
                                     tdy_yzz_yyy_0, tdy_yzz_yyz_0, tdy_yzz_yzz_0, tdy_yzz_zzz_0, tdy_zzz_xxx_0, \
                                     tdy_zzz_xxy_0, tdy_zzz_xxz_0, tdy_zzz_xyy_0, tdy_zzz_xyz_0, tdy_zzz_xzz_0, \
                                     tdy_zzz_yyy_0, tdy_zzz_yyz_0, tdy_zzz_yzz_0, tdy_zzz_zzz_0, tdz_yzz_xyy_0, \
                                     tdz_yzz_xyz_0, tdz_yzz_xzz_0, tdz_yzz_yyy_0, tdz_yzz_yyz_0, tdz_yzz_yzz_0, \
                                     tdz_yzz_zzz_0, tdz_zzz_xxx_0, tdz_zzz_xxy_0, tdz_zzz_xxz_0, tdz_zzz_xyy_0, \
                                     tdz_zzz_xyz_0, tdz_zzz_xzz_0, tdz_zzz_yyy_0, tdz_zzz_yyz_0, tdz_zzz_yzz_0, \
                                     tdz_zzz_zzz_0, tlx_xyzz_xyz_0, tlx_xyzz_xzz_0, tlx_xyzz_yyy_0, tlx_xyzz_yyz_0, \
                                     tlx_xyzz_yzz_0, tlx_xyzz_zzz_0, tlx_xzzz_xxx_0, tlx_xzzz_xxy_0, tlx_xzzz_xxz_0, \
                                     tlx_xzzz_xyy_0, tlx_xzzz_xyz_0, tlx_xzzz_xzz_0, tlx_xzzz_yyy_0, tlx_xzzz_yyz_0, \
                                     tlx_xzzz_yzz_0, tlx_xzzz_zzz_0, tlx_yzz_xyz_0, tlx_yzz_xzz_0, tlx_yzz_yyy_0, \
                                     tlx_yzz_yyz_0, tlx_yzz_yz_0, tlx_yzz_yzz_0, tlx_yzz_zz_0, tlx_yzz_zzz_0, \
                                     tlx_zzz_xx_0, tlx_zzz_xxx_0, tlx_zzz_xxy_0, tlx_zzz_xxz_0, tlx_zzz_xy_0, \
                                     tlx_zzz_xyy_0, tlx_zzz_xyz_0, tlx_zzz_xz_0, tlx_zzz_xzz_0, tlx_zzz_yy_0, \
                                     tlx_zzz_yyy_0, tlx_zzz_yyz_0, tlx_zzz_yz_0, tlx_zzz_yzz_0, tlx_zzz_zz_0, \
                                     tlx_zzz_zzz_0, tly_xyzz_xyy_0, tly_xyzz_xyz_0, tly_xyzz_xzz_0, tly_xyzz_yyy_0, \
                                     tly_xyzz_yyz_0, tly_xyzz_yzz_0, tly_xyzz_zzz_0, tly_xzzz_xxx_0, tly_xzzz_xxy_0, \
                                     tly_xzzz_xxz_0, tly_xzzz_xyy_0, tly_xzzz_xyz_0, tly_xzzz_xzz_0, tly_xzzz_yyy_0, \
                                     tly_xzzz_yyz_0, tly_xzzz_yzz_0, tly_xzzz_zzz_0, tly_yzz_xyy_0, tly_yzz_xyz_0, \
                                     tly_yzz_xzz_0, tly_yzz_yy_0, tly_yzz_yyy_0, tly_yzz_yyz_0, tly_yzz_yz_0, \
                                     tly_yzz_yzz_0, tly_yzz_zz_0, tly_yzz_zzz_0, tly_zzz_xx_0, tly_zzz_xxx_0, \
                                     tly_zzz_xxy_0, tly_zzz_xxz_0, tly_zzz_xy_0, tly_zzz_xyy_0, tly_zzz_xyz_0, \
                                     tly_zzz_xz_0, tly_zzz_xzz_0, tly_zzz_yy_0, tly_zzz_yyy_0, tly_zzz_yyz_0, \
                                     tly_zzz_yz_0, tly_zzz_yzz_0, tly_zzz_zz_0, tly_zzz_zzz_0, tlz_xyzz_xyy_0, \
                                     tlz_xyzz_xyz_0, tlz_xyzz_xzz_0, tlz_xyzz_yyy_0, tlz_xyzz_yyz_0, tlz_xyzz_yzz_0, \
                                     tlz_xyzz_zzz_0, tlz_xzzz_xxx_0, tlz_xzzz_xxy_0, tlz_xzzz_xxz_0, tlz_xzzz_xyy_0, \
                                     tlz_xzzz_xyz_0, tlz_xzzz_xzz_0, tlz_xzzz_yyy_0, tlz_xzzz_yyz_0, tlz_xzzz_yzz_0, \
                                     tlz_xzzz_zzz_0, tlz_yzz_xyy_0, tlz_yzz_xyz_0, tlz_yzz_xzz_0, tlz_yzz_yy_0, \
                                     tlz_yzz_yyy_0, tlz_yzz_yyz_0, tlz_yzz_yz_0, tlz_yzz_yzz_0, tlz_yzz_zz_0, \
                                     tlz_yzz_zzz_0, tlz_zzz_xx_0, tlz_zzz_xxx_0, tlz_zzz_xxy_0, tlz_zzz_xxz_0, \
                                     tlz_zzz_xy_0, tlz_zzz_xyy_0, tlz_zzz_xyz_0, tlz_zzz_xz_0, tlz_zzz_xzz_0, \
                                     tlz_zzz_yy_0, tlz_zzz_yyy_0, tlz_zzz_yyz_0, tlz_zzz_yz_0, tlz_zzz_yzz_0, \
                                     tlz_zzz_zz_0, tlz_zzz_zzz_0, tpy_yzz_xyy_0, tpy_yzz_xyz_0, tpy_yzz_xzz_0, \
                                     tpy_yzz_yyy_0, tpy_yzz_yyz_0, tpy_yzz_yzz_0, tpy_yzz_zzz_0, tpy_zzz_xxx_0, \
                                     tpy_zzz_xxy_0, tpy_zzz_xxz_0, tpy_zzz_xyy_0, tpy_zzz_xyz_0, tpy_zzz_xzz_0, \
                                     tpy_zzz_yyy_0, tpy_zzz_yyz_0, tpy_zzz_yzz_0, tpy_zzz_zzz_0, tpz_yzz_xyy_0, \
                                     tpz_yzz_xyz_0, tpz_yzz_xzz_0, tpz_yzz_yyy_0, tpz_yzz_yyz_0, tpz_yzz_yzz_0, \
                                     tpz_yzz_zzz_0, tpz_zzz_xxx_0, tpz_zzz_xxy_0, tpz_zzz_xxz_0, tpz_zzz_xyy_0, \
                                     tpz_zzz_xyz_0, tpz_zzz_xzz_0, tpz_zzz_yyy_0, tpz_zzz_yyz_0, tpz_zzz_yzz_0, \
                                     tpz_zzz_zzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tly_xyzz_xyy_0[j] =
                pa_x[j] * tly_yzz_xyy_0[j] + 0.5 * fl1_fx * tly_yzz_yy_0[j] + 0.5 * fl1_fx * tpz_yzz_xyy_0[j] + fl1_fx * fl1_fgb * tdz_yzz_xyy_0[j];

            tlz_xyzz_xyy_0[j] =
                pa_x[j] * tlz_yzz_xyy_0[j] + 0.5 * fl1_fx * tlz_yzz_yy_0[j] - 0.5 * fl1_fx * tpy_yzz_xyy_0[j] - fl1_fx * fl1_fgb * tdy_yzz_xyy_0[j];

            tlx_xyzz_xyz_0[j] = pa_x[j] * tlx_yzz_xyz_0[j] + 0.5 * fl1_fx * tlx_yzz_yz_0[j];

            tly_xyzz_xyz_0[j] =
                pa_x[j] * tly_yzz_xyz_0[j] + 0.5 * fl1_fx * tly_yzz_yz_0[j] + 0.5 * fl1_fx * tpz_yzz_xyz_0[j] + fl1_fx * fl1_fgb * tdz_yzz_xyz_0[j];

            tlz_xyzz_xyz_0[j] =
                pa_x[j] * tlz_yzz_xyz_0[j] + 0.5 * fl1_fx * tlz_yzz_yz_0[j] - 0.5 * fl1_fx * tpy_yzz_xyz_0[j] - fl1_fx * fl1_fgb * tdy_yzz_xyz_0[j];

            tlx_xyzz_xzz_0[j] = pa_x[j] * tlx_yzz_xzz_0[j] + 0.5 * fl1_fx * tlx_yzz_zz_0[j];

            tly_xyzz_xzz_0[j] =
                pa_x[j] * tly_yzz_xzz_0[j] + 0.5 * fl1_fx * tly_yzz_zz_0[j] + 0.5 * fl1_fx * tpz_yzz_xzz_0[j] + fl1_fx * fl1_fgb * tdz_yzz_xzz_0[j];

            tlz_xyzz_xzz_0[j] =
                pa_x[j] * tlz_yzz_xzz_0[j] + 0.5 * fl1_fx * tlz_yzz_zz_0[j] - 0.5 * fl1_fx * tpy_yzz_xzz_0[j] - fl1_fx * fl1_fgb * tdy_yzz_xzz_0[j];

            tlx_xyzz_yyy_0[j] = pa_x[j] * tlx_yzz_yyy_0[j];

            tly_xyzz_yyy_0[j] = pa_x[j] * tly_yzz_yyy_0[j] + 0.5 * fl1_fx * tpz_yzz_yyy_0[j] + fl1_fx * fl1_fgb * tdz_yzz_yyy_0[j];

            tlz_xyzz_yyy_0[j] = pa_x[j] * tlz_yzz_yyy_0[j] - 0.5 * fl1_fx * tpy_yzz_yyy_0[j] - fl1_fx * fl1_fgb * tdy_yzz_yyy_0[j];

            tlx_xyzz_yyz_0[j] = pa_x[j] * tlx_yzz_yyz_0[j];

            tly_xyzz_yyz_0[j] = pa_x[j] * tly_yzz_yyz_0[j] + 0.5 * fl1_fx * tpz_yzz_yyz_0[j] + fl1_fx * fl1_fgb * tdz_yzz_yyz_0[j];

            tlz_xyzz_yyz_0[j] = pa_x[j] * tlz_yzz_yyz_0[j] - 0.5 * fl1_fx * tpy_yzz_yyz_0[j] - fl1_fx * fl1_fgb * tdy_yzz_yyz_0[j];

            tlx_xyzz_yzz_0[j] = pa_x[j] * tlx_yzz_yzz_0[j];

            tly_xyzz_yzz_0[j] = pa_x[j] * tly_yzz_yzz_0[j] + 0.5 * fl1_fx * tpz_yzz_yzz_0[j] + fl1_fx * fl1_fgb * tdz_yzz_yzz_0[j];

            tlz_xyzz_yzz_0[j] = pa_x[j] * tlz_yzz_yzz_0[j] - 0.5 * fl1_fx * tpy_yzz_yzz_0[j] - fl1_fx * fl1_fgb * tdy_yzz_yzz_0[j];

            tlx_xyzz_zzz_0[j] = pa_x[j] * tlx_yzz_zzz_0[j];

            tly_xyzz_zzz_0[j] = pa_x[j] * tly_yzz_zzz_0[j] + 0.5 * fl1_fx * tpz_yzz_zzz_0[j] + fl1_fx * fl1_fgb * tdz_yzz_zzz_0[j];

            tlz_xyzz_zzz_0[j] = pa_x[j] * tlz_yzz_zzz_0[j] - 0.5 * fl1_fx * tpy_yzz_zzz_0[j] - fl1_fx * fl1_fgb * tdy_yzz_zzz_0[j];

            tlx_xzzz_xxx_0[j] = pa_x[j] * tlx_zzz_xxx_0[j] + 1.5 * fl1_fx * tlx_zzz_xx_0[j];

            tly_xzzz_xxx_0[j] =
                pa_x[j] * tly_zzz_xxx_0[j] + 1.5 * fl1_fx * tly_zzz_xx_0[j] + 0.5 * fl1_fx * tpz_zzz_xxx_0[j] + fl1_fx * fl1_fgb * tdz_zzz_xxx_0[j];

            tlz_xzzz_xxx_0[j] =
                pa_x[j] * tlz_zzz_xxx_0[j] + 1.5 * fl1_fx * tlz_zzz_xx_0[j] - 0.5 * fl1_fx * tpy_zzz_xxx_0[j] - fl1_fx * fl1_fgb * tdy_zzz_xxx_0[j];

            tlx_xzzz_xxy_0[j] = pa_x[j] * tlx_zzz_xxy_0[j] + fl1_fx * tlx_zzz_xy_0[j];

            tly_xzzz_xxy_0[j] =
                pa_x[j] * tly_zzz_xxy_0[j] + fl1_fx * tly_zzz_xy_0[j] + 0.5 * fl1_fx * tpz_zzz_xxy_0[j] + fl1_fx * fl1_fgb * tdz_zzz_xxy_0[j];

            tlz_xzzz_xxy_0[j] =
                pa_x[j] * tlz_zzz_xxy_0[j] + fl1_fx * tlz_zzz_xy_0[j] - 0.5 * fl1_fx * tpy_zzz_xxy_0[j] - fl1_fx * fl1_fgb * tdy_zzz_xxy_0[j];

            tlx_xzzz_xxz_0[j] = pa_x[j] * tlx_zzz_xxz_0[j] + fl1_fx * tlx_zzz_xz_0[j];

            tly_xzzz_xxz_0[j] =
                pa_x[j] * tly_zzz_xxz_0[j] + fl1_fx * tly_zzz_xz_0[j] + 0.5 * fl1_fx * tpz_zzz_xxz_0[j] + fl1_fx * fl1_fgb * tdz_zzz_xxz_0[j];

            tlz_xzzz_xxz_0[j] =
                pa_x[j] * tlz_zzz_xxz_0[j] + fl1_fx * tlz_zzz_xz_0[j] - 0.5 * fl1_fx * tpy_zzz_xxz_0[j] - fl1_fx * fl1_fgb * tdy_zzz_xxz_0[j];

            tlx_xzzz_xyy_0[j] = pa_x[j] * tlx_zzz_xyy_0[j] + 0.5 * fl1_fx * tlx_zzz_yy_0[j];

            tly_xzzz_xyy_0[j] =
                pa_x[j] * tly_zzz_xyy_0[j] + 0.5 * fl1_fx * tly_zzz_yy_0[j] + 0.5 * fl1_fx * tpz_zzz_xyy_0[j] + fl1_fx * fl1_fgb * tdz_zzz_xyy_0[j];

            tlz_xzzz_xyy_0[j] =
                pa_x[j] * tlz_zzz_xyy_0[j] + 0.5 * fl1_fx * tlz_zzz_yy_0[j] - 0.5 * fl1_fx * tpy_zzz_xyy_0[j] - fl1_fx * fl1_fgb * tdy_zzz_xyy_0[j];

            tlx_xzzz_xyz_0[j] = pa_x[j] * tlx_zzz_xyz_0[j] + 0.5 * fl1_fx * tlx_zzz_yz_0[j];

            tly_xzzz_xyz_0[j] =
                pa_x[j] * tly_zzz_xyz_0[j] + 0.5 * fl1_fx * tly_zzz_yz_0[j] + 0.5 * fl1_fx * tpz_zzz_xyz_0[j] + fl1_fx * fl1_fgb * tdz_zzz_xyz_0[j];

            tlz_xzzz_xyz_0[j] =
                pa_x[j] * tlz_zzz_xyz_0[j] + 0.5 * fl1_fx * tlz_zzz_yz_0[j] - 0.5 * fl1_fx * tpy_zzz_xyz_0[j] - fl1_fx * fl1_fgb * tdy_zzz_xyz_0[j];

            tlx_xzzz_xzz_0[j] = pa_x[j] * tlx_zzz_xzz_0[j] + 0.5 * fl1_fx * tlx_zzz_zz_0[j];

            tly_xzzz_xzz_0[j] =
                pa_x[j] * tly_zzz_xzz_0[j] + 0.5 * fl1_fx * tly_zzz_zz_0[j] + 0.5 * fl1_fx * tpz_zzz_xzz_0[j] + fl1_fx * fl1_fgb * tdz_zzz_xzz_0[j];

            tlz_xzzz_xzz_0[j] =
                pa_x[j] * tlz_zzz_xzz_0[j] + 0.5 * fl1_fx * tlz_zzz_zz_0[j] - 0.5 * fl1_fx * tpy_zzz_xzz_0[j] - fl1_fx * fl1_fgb * tdy_zzz_xzz_0[j];

            tlx_xzzz_yyy_0[j] = pa_x[j] * tlx_zzz_yyy_0[j];

            tly_xzzz_yyy_0[j] = pa_x[j] * tly_zzz_yyy_0[j] + 0.5 * fl1_fx * tpz_zzz_yyy_0[j] + fl1_fx * fl1_fgb * tdz_zzz_yyy_0[j];

            tlz_xzzz_yyy_0[j] = pa_x[j] * tlz_zzz_yyy_0[j] - 0.5 * fl1_fx * tpy_zzz_yyy_0[j] - fl1_fx * fl1_fgb * tdy_zzz_yyy_0[j];

            tlx_xzzz_yyz_0[j] = pa_x[j] * tlx_zzz_yyz_0[j];

            tly_xzzz_yyz_0[j] = pa_x[j] * tly_zzz_yyz_0[j] + 0.5 * fl1_fx * tpz_zzz_yyz_0[j] + fl1_fx * fl1_fgb * tdz_zzz_yyz_0[j];

            tlz_xzzz_yyz_0[j] = pa_x[j] * tlz_zzz_yyz_0[j] - 0.5 * fl1_fx * tpy_zzz_yyz_0[j] - fl1_fx * fl1_fgb * tdy_zzz_yyz_0[j];

            tlx_xzzz_yzz_0[j] = pa_x[j] * tlx_zzz_yzz_0[j];

            tly_xzzz_yzz_0[j] = pa_x[j] * tly_zzz_yzz_0[j] + 0.5 * fl1_fx * tpz_zzz_yzz_0[j] + fl1_fx * fl1_fgb * tdz_zzz_yzz_0[j];

            tlz_xzzz_yzz_0[j] = pa_x[j] * tlz_zzz_yzz_0[j] - 0.5 * fl1_fx * tpy_zzz_yzz_0[j] - fl1_fx * fl1_fgb * tdy_zzz_yzz_0[j];

            tlx_xzzz_zzz_0[j] = pa_x[j] * tlx_zzz_zzz_0[j];

            tly_xzzz_zzz_0[j] = pa_x[j] * tly_zzz_zzz_0[j] + 0.5 * fl1_fx * tpz_zzz_zzz_0[j] + fl1_fx * fl1_fgb * tdz_zzz_zzz_0[j];

            tlz_xzzz_zzz_0[j] = pa_x[j] * tlz_zzz_zzz_0[j] - 0.5 * fl1_fx * tpy_zzz_zzz_0[j] - fl1_fx * fl1_fgb * tdy_zzz_zzz_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForGF_300_350(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_l_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_4_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

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

        auto tlx_yyy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 36);

        auto tly_yyy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 36);

        auto tlz_yyy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 36);

        auto tlx_yyy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 37);

        auto tly_yyy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 37);

        auto tlz_yyy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 37);

        auto tlx_yyy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 38);

        auto tly_yyy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 38);

        auto tlz_yyy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 38);

        auto tlx_yyy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 39);

        auto tly_yyy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 39);

        auto tlz_yyy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 39);

        auto tlx_yyy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 40);

        auto tly_yyy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 40);

        auto tlz_yyy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 40);

        auto tlx_yyy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 41);

        auto tly_yyy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 41);

        auto tlz_yyy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 41);

        auto tlx_yyz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 42);

        auto tly_yyz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 42);

        auto tlz_yyz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 42);

        auto tlx_yyz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 43);

        auto tly_yyz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 43);

        auto tlz_yyz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 43);

        auto tlx_yyz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 44);

        auto tly_yyz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 44);

        auto tlz_yyz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 44);

        auto tlx_yyz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 45);

        auto tly_yyz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 45);

        auto tpx_yyy_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 60);

        auto tpz_yyy_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 60);

        auto tpx_yyy_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 61);

        auto tpz_yyy_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 61);

        auto tpx_yyy_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 62);

        auto tpz_yyy_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 62);

        auto tpx_yyy_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 63);

        auto tpz_yyy_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 63);

        auto tpx_yyy_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 64);

        auto tpz_yyy_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 64);

        auto tpx_yyy_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 65);

        auto tpz_yyy_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 65);

        auto tpx_yyy_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 66);

        auto tpz_yyy_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 66);

        auto tpx_yyy_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 67);

        auto tpz_yyy_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 67);

        auto tpx_yyy_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 68);

        auto tpz_yyy_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 68);

        auto tpx_yyy_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 69);

        auto tpz_yyy_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 69);

        auto tpx_yyz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 70);

        auto tpz_yyz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 70);

        auto tpx_yyz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 71);

        auto tpz_yyz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 71);

        auto tpx_yyz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 72);

        auto tpz_yyz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 72);

        auto tpx_yyz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 73);

        auto tpz_yyz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 73);

        auto tpx_yyz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 74);

        auto tpz_yyz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 74);

        auto tpx_yyz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 75);

        auto tpz_yyz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 75);

        auto tpz_yyz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 76);

        auto tdx_yyy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 60);

        auto tdz_yyy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 60);

        auto tdx_yyy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 61);

        auto tdz_yyy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 61);

        auto tdx_yyy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 62);

        auto tdz_yyy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 62);

        auto tdx_yyy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 63);

        auto tdz_yyy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 63);

        auto tdx_yyy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 64);

        auto tdz_yyy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 64);

        auto tdx_yyy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 65);

        auto tdz_yyy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 65);

        auto tdx_yyy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 66);

        auto tdz_yyy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 66);

        auto tdx_yyy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 67);

        auto tdz_yyy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 67);

        auto tdx_yyy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 68);

        auto tdz_yyy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 68);

        auto tdx_yyy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 69);

        auto tdz_yyy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 69);

        auto tdx_yyz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 70);

        auto tdz_yyz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 70);

        auto tdx_yyz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 71);

        auto tdz_yyz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 71);

        auto tdx_yyz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 72);

        auto tdz_yyz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 72);

        auto tdx_yyz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 73);

        auto tdz_yyz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 73);

        auto tdx_yyz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 74);

        auto tdz_yyz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 74);

        auto tdx_yyz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 75);

        auto tdz_yyz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 75);

        auto tdz_yyz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 76);

        // set up pointers to integrals

        auto tlx_yyyy_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 100);

        auto tly_yyyy_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 100);

        auto tlz_yyyy_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 100);

        auto tlx_yyyy_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 101);

        auto tly_yyyy_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 101);

        auto tlz_yyyy_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 101);

        auto tlx_yyyy_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 102);

        auto tly_yyyy_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 102);

        auto tlz_yyyy_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 102);

        auto tlx_yyyy_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 103);

        auto tly_yyyy_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 103);

        auto tlz_yyyy_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 103);

        auto tlx_yyyy_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 104);

        auto tly_yyyy_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 104);

        auto tlz_yyyy_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 104);

        auto tlx_yyyy_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 105);

        auto tly_yyyy_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 105);

        auto tlz_yyyy_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 105);

        auto tlx_yyyy_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 106);

        auto tly_yyyy_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 106);

        auto tlz_yyyy_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 106);

        auto tlx_yyyy_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 107);

        auto tly_yyyy_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 107);

        auto tlz_yyyy_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 107);

        auto tlx_yyyy_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 108);

        auto tly_yyyy_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 108);

        auto tlz_yyyy_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 108);

        auto tlx_yyyy_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 109);

        auto tly_yyyy_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 109);

        auto tlz_yyyy_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 109);

        auto tlx_yyyz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 110);

        auto tly_yyyz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 110);

        auto tlz_yyyz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 110);

        auto tlx_yyyz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 111);

        auto tly_yyyz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 111);

        auto tlz_yyyz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 111);

        auto tlx_yyyz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 112);

        auto tly_yyyz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 112);

        auto tlz_yyyz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 112);

        auto tlx_yyyz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 113);

        auto tly_yyyz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 113);

        auto tlz_yyyz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 113);

        auto tlx_yyyz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 114);

        auto tly_yyyz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 114);

        auto tlz_yyyz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 114);

        auto tlx_yyyz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 115);

        auto tly_yyyz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 115);

        auto tlz_yyyz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 115);

        auto tlx_yyyz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 116);

        auto tly_yyyz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 116);

        // Batch of Integrals (300,350)

        #pragma omp simd aligned(fgb, fx, pa_y, tdx_yyy_xxx_0, tdx_yyy_xxy_0, tdx_yyy_xxz_0, \
                                     tdx_yyy_xyy_0, tdx_yyy_xyz_0, tdx_yyy_xzz_0, tdx_yyy_yyy_0, tdx_yyy_yyz_0, \
                                     tdx_yyy_yzz_0, tdx_yyy_zzz_0, tdx_yyz_xxx_0, tdx_yyz_xxy_0, tdx_yyz_xxz_0, \
                                     tdx_yyz_xyy_0, tdx_yyz_xyz_0, tdx_yyz_xzz_0, tdz_yyy_xxx_0, tdz_yyy_xxy_0, \
                                     tdz_yyy_xxz_0, tdz_yyy_xyy_0, tdz_yyy_xyz_0, tdz_yyy_xzz_0, tdz_yyy_yyy_0, \
                                     tdz_yyy_yyz_0, tdz_yyy_yzz_0, tdz_yyy_zzz_0, tdz_yyz_xxx_0, tdz_yyz_xxy_0, \
                                     tdz_yyz_xxz_0, tdz_yyz_xyy_0, tdz_yyz_xyz_0, tdz_yyz_xzz_0, tdz_yyz_yyy_0, \
                                     tlx_yy_xxx_0, tlx_yy_xxy_0, tlx_yy_xxz_0, tlx_yy_xyy_0, tlx_yy_xyz_0, tlx_yy_xzz_0, \
                                     tlx_yy_yyy_0, tlx_yy_yyz_0, tlx_yy_yzz_0, tlx_yy_zzz_0, tlx_yyy_xx_0, \
                                     tlx_yyy_xxx_0, tlx_yyy_xxy_0, tlx_yyy_xxz_0, tlx_yyy_xy_0, tlx_yyy_xyy_0, \
                                     tlx_yyy_xyz_0, tlx_yyy_xz_0, tlx_yyy_xzz_0, tlx_yyy_yy_0, tlx_yyy_yyy_0, \
                                     tlx_yyy_yyz_0, tlx_yyy_yz_0, tlx_yyy_yzz_0, tlx_yyy_zz_0, tlx_yyy_zzz_0, \
                                     tlx_yyyy_xxx_0, tlx_yyyy_xxy_0, tlx_yyyy_xxz_0, tlx_yyyy_xyy_0, tlx_yyyy_xyz_0, \
                                     tlx_yyyy_xzz_0, tlx_yyyy_yyy_0, tlx_yyyy_yyz_0, tlx_yyyy_yzz_0, tlx_yyyy_zzz_0, \
                                     tlx_yyyz_xxx_0, tlx_yyyz_xxy_0, tlx_yyyz_xxz_0, tlx_yyyz_xyy_0, tlx_yyyz_xyz_0, \
                                     tlx_yyyz_xzz_0, tlx_yyyz_yyy_0, tlx_yyz_xx_0, tlx_yyz_xxx_0, tlx_yyz_xxy_0, \
                                     tlx_yyz_xxz_0, tlx_yyz_xy_0, tlx_yyz_xyy_0, tlx_yyz_xyz_0, tlx_yyz_xz_0, \
                                     tlx_yyz_xzz_0, tlx_yyz_yy_0, tlx_yyz_yyy_0, tlx_yz_xxx_0, tlx_yz_xxy_0, \
                                     tlx_yz_xxz_0, tlx_yz_xyy_0, tlx_yz_xyz_0, tlx_yz_xzz_0, tlx_yz_yyy_0, tly_yy_xxx_0, \
                                     tly_yy_xxy_0, tly_yy_xxz_0, tly_yy_xyy_0, tly_yy_xyz_0, tly_yy_xzz_0, tly_yy_yyy_0, \
                                     tly_yy_yyz_0, tly_yy_yzz_0, tly_yy_zzz_0, tly_yyy_xx_0, tly_yyy_xxx_0, \
                                     tly_yyy_xxy_0, tly_yyy_xxz_0, tly_yyy_xy_0, tly_yyy_xyy_0, tly_yyy_xyz_0, \
                                     tly_yyy_xz_0, tly_yyy_xzz_0, tly_yyy_yy_0, tly_yyy_yyy_0, tly_yyy_yyz_0, \
                                     tly_yyy_yz_0, tly_yyy_yzz_0, tly_yyy_zz_0, tly_yyy_zzz_0, tly_yyyy_xxx_0, \
                                     tly_yyyy_xxy_0, tly_yyyy_xxz_0, tly_yyyy_xyy_0, tly_yyyy_xyz_0, tly_yyyy_xzz_0, \
                                     tly_yyyy_yyy_0, tly_yyyy_yyz_0, tly_yyyy_yzz_0, tly_yyyy_zzz_0, tly_yyyz_xxx_0, \
                                     tly_yyyz_xxy_0, tly_yyyz_xxz_0, tly_yyyz_xyy_0, tly_yyyz_xyz_0, tly_yyyz_xzz_0, \
                                     tly_yyyz_yyy_0, tly_yyz_xx_0, tly_yyz_xxx_0, tly_yyz_xxy_0, tly_yyz_xxz_0, \
                                     tly_yyz_xy_0, tly_yyz_xyy_0, tly_yyz_xyz_0, tly_yyz_xz_0, tly_yyz_xzz_0, \
                                     tly_yyz_yy_0, tly_yyz_yyy_0, tly_yz_xxx_0, tly_yz_xxy_0, tly_yz_xxz_0, \
                                     tly_yz_xyy_0, tly_yz_xyz_0, tly_yz_xzz_0, tly_yz_yyy_0, tlz_yy_xxx_0, tlz_yy_xxy_0, \
                                     tlz_yy_xxz_0, tlz_yy_xyy_0, tlz_yy_xyz_0, tlz_yy_xzz_0, tlz_yy_yyy_0, tlz_yy_yyz_0, \
                                     tlz_yy_yzz_0, tlz_yy_zzz_0, tlz_yyy_xx_0, tlz_yyy_xxx_0, tlz_yyy_xxy_0, \
                                     tlz_yyy_xxz_0, tlz_yyy_xy_0, tlz_yyy_xyy_0, tlz_yyy_xyz_0, tlz_yyy_xz_0, \
                                     tlz_yyy_xzz_0, tlz_yyy_yy_0, tlz_yyy_yyy_0, tlz_yyy_yyz_0, tlz_yyy_yz_0, \
                                     tlz_yyy_yzz_0, tlz_yyy_zz_0, tlz_yyy_zzz_0, tlz_yyyy_xxx_0, tlz_yyyy_xxy_0, \
                                     tlz_yyyy_xxz_0, tlz_yyyy_xyy_0, tlz_yyyy_xyz_0, tlz_yyyy_xzz_0, tlz_yyyy_yyy_0, \
                                     tlz_yyyy_yyz_0, tlz_yyyy_yzz_0, tlz_yyyy_zzz_0, tlz_yyyz_xxx_0, tlz_yyyz_xxy_0, \
                                     tlz_yyyz_xxz_0, tlz_yyyz_xyy_0, tlz_yyyz_xyz_0, tlz_yyyz_xzz_0, tlz_yyz_xx_0, \
                                     tlz_yyz_xxx_0, tlz_yyz_xxy_0, tlz_yyz_xxz_0, tlz_yyz_xy_0, tlz_yyz_xyy_0, \
                                     tlz_yyz_xyz_0, tlz_yyz_xz_0, tlz_yyz_xzz_0, tlz_yz_xxx_0, tlz_yz_xxy_0, \
                                     tlz_yz_xxz_0, tlz_yz_xyy_0, tlz_yz_xyz_0, tlz_yz_xzz_0, tpx_yyy_xxx_0, \
                                     tpx_yyy_xxy_0, tpx_yyy_xxz_0, tpx_yyy_xyy_0, tpx_yyy_xyz_0, tpx_yyy_xzz_0, \
                                     tpx_yyy_yyy_0, tpx_yyy_yyz_0, tpx_yyy_yzz_0, tpx_yyy_zzz_0, tpx_yyz_xxx_0, \
                                     tpx_yyz_xxy_0, tpx_yyz_xxz_0, tpx_yyz_xyy_0, tpx_yyz_xyz_0, tpx_yyz_xzz_0, \
                                     tpz_yyy_xxx_0, tpz_yyy_xxy_0, tpz_yyy_xxz_0, tpz_yyy_xyy_0, tpz_yyy_xyz_0, \
                                     tpz_yyy_xzz_0, tpz_yyy_yyy_0, tpz_yyy_yyz_0, tpz_yyy_yzz_0, tpz_yyy_zzz_0, \
                                     tpz_yyz_xxx_0, tpz_yyz_xxy_0, tpz_yyz_xxz_0, tpz_yyz_xyy_0, tpz_yyz_xyz_0, \
                                     tpz_yyz_xzz_0, tpz_yyz_yyy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_yyyy_xxx_0[j] =
                pa_y[j] * tlx_yyy_xxx_0[j] + 1.5 * fl1_fx * tlx_yy_xxx_0[j] - 0.5 * fl1_fx * tpz_yyy_xxx_0[j] - fl1_fx * fl1_fgb * tdz_yyy_xxx_0[j];

            tly_yyyy_xxx_0[j] = pa_y[j] * tly_yyy_xxx_0[j] + 1.5 * fl1_fx * tly_yy_xxx_0[j];

            tlz_yyyy_xxx_0[j] =
                pa_y[j] * tlz_yyy_xxx_0[j] + 1.5 * fl1_fx * tlz_yy_xxx_0[j] + 0.5 * fl1_fx * tpx_yyy_xxx_0[j] + fl1_fx * fl1_fgb * tdx_yyy_xxx_0[j];

            tlx_yyyy_xxy_0[j] = pa_y[j] * tlx_yyy_xxy_0[j] + 1.5 * fl1_fx * tlx_yy_xxy_0[j] + 0.5 * fl1_fx * tlx_yyy_xx_0[j] -
                                0.5 * fl1_fx * tpz_yyy_xxy_0[j] - fl1_fx * fl1_fgb * tdz_yyy_xxy_0[j];

            tly_yyyy_xxy_0[j] = pa_y[j] * tly_yyy_xxy_0[j] + 1.5 * fl1_fx * tly_yy_xxy_0[j] + 0.5 * fl1_fx * tly_yyy_xx_0[j];

            tlz_yyyy_xxy_0[j] = pa_y[j] * tlz_yyy_xxy_0[j] + 1.5 * fl1_fx * tlz_yy_xxy_0[j] + 0.5 * fl1_fx * tlz_yyy_xx_0[j] +
                                0.5 * fl1_fx * tpx_yyy_xxy_0[j] + fl1_fx * fl1_fgb * tdx_yyy_xxy_0[j];

            tlx_yyyy_xxz_0[j] =
                pa_y[j] * tlx_yyy_xxz_0[j] + 1.5 * fl1_fx * tlx_yy_xxz_0[j] - 0.5 * fl1_fx * tpz_yyy_xxz_0[j] - fl1_fx * fl1_fgb * tdz_yyy_xxz_0[j];

            tly_yyyy_xxz_0[j] = pa_y[j] * tly_yyy_xxz_0[j] + 1.5 * fl1_fx * tly_yy_xxz_0[j];

            tlz_yyyy_xxz_0[j] =
                pa_y[j] * tlz_yyy_xxz_0[j] + 1.5 * fl1_fx * tlz_yy_xxz_0[j] + 0.5 * fl1_fx * tpx_yyy_xxz_0[j] + fl1_fx * fl1_fgb * tdx_yyy_xxz_0[j];

            tlx_yyyy_xyy_0[j] = pa_y[j] * tlx_yyy_xyy_0[j] + 1.5 * fl1_fx * tlx_yy_xyy_0[j] + fl1_fx * tlx_yyy_xy_0[j] -
                                0.5 * fl1_fx * tpz_yyy_xyy_0[j] - fl1_fx * fl1_fgb * tdz_yyy_xyy_0[j];

            tly_yyyy_xyy_0[j] = pa_y[j] * tly_yyy_xyy_0[j] + 1.5 * fl1_fx * tly_yy_xyy_0[j] + fl1_fx * tly_yyy_xy_0[j];

            tlz_yyyy_xyy_0[j] = pa_y[j] * tlz_yyy_xyy_0[j] + 1.5 * fl1_fx * tlz_yy_xyy_0[j] + fl1_fx * tlz_yyy_xy_0[j] +
                                0.5 * fl1_fx * tpx_yyy_xyy_0[j] + fl1_fx * fl1_fgb * tdx_yyy_xyy_0[j];

            tlx_yyyy_xyz_0[j] = pa_y[j] * tlx_yyy_xyz_0[j] + 1.5 * fl1_fx * tlx_yy_xyz_0[j] + 0.5 * fl1_fx * tlx_yyy_xz_0[j] -
                                0.5 * fl1_fx * tpz_yyy_xyz_0[j] - fl1_fx * fl1_fgb * tdz_yyy_xyz_0[j];

            tly_yyyy_xyz_0[j] = pa_y[j] * tly_yyy_xyz_0[j] + 1.5 * fl1_fx * tly_yy_xyz_0[j] + 0.5 * fl1_fx * tly_yyy_xz_0[j];

            tlz_yyyy_xyz_0[j] = pa_y[j] * tlz_yyy_xyz_0[j] + 1.5 * fl1_fx * tlz_yy_xyz_0[j] + 0.5 * fl1_fx * tlz_yyy_xz_0[j] +
                                0.5 * fl1_fx * tpx_yyy_xyz_0[j] + fl1_fx * fl1_fgb * tdx_yyy_xyz_0[j];

            tlx_yyyy_xzz_0[j] =
                pa_y[j] * tlx_yyy_xzz_0[j] + 1.5 * fl1_fx * tlx_yy_xzz_0[j] - 0.5 * fl1_fx * tpz_yyy_xzz_0[j] - fl1_fx * fl1_fgb * tdz_yyy_xzz_0[j];

            tly_yyyy_xzz_0[j] = pa_y[j] * tly_yyy_xzz_0[j] + 1.5 * fl1_fx * tly_yy_xzz_0[j];

            tlz_yyyy_xzz_0[j] =
                pa_y[j] * tlz_yyy_xzz_0[j] + 1.5 * fl1_fx * tlz_yy_xzz_0[j] + 0.5 * fl1_fx * tpx_yyy_xzz_0[j] + fl1_fx * fl1_fgb * tdx_yyy_xzz_0[j];

            tlx_yyyy_yyy_0[j] = pa_y[j] * tlx_yyy_yyy_0[j] + 1.5 * fl1_fx * tlx_yy_yyy_0[j] + 1.5 * fl1_fx * tlx_yyy_yy_0[j] -
                                0.5 * fl1_fx * tpz_yyy_yyy_0[j] - fl1_fx * fl1_fgb * tdz_yyy_yyy_0[j];

            tly_yyyy_yyy_0[j] = pa_y[j] * tly_yyy_yyy_0[j] + 1.5 * fl1_fx * tly_yy_yyy_0[j] + 1.5 * fl1_fx * tly_yyy_yy_0[j];

            tlz_yyyy_yyy_0[j] = pa_y[j] * tlz_yyy_yyy_0[j] + 1.5 * fl1_fx * tlz_yy_yyy_0[j] + 1.5 * fl1_fx * tlz_yyy_yy_0[j] +
                                0.5 * fl1_fx * tpx_yyy_yyy_0[j] + fl1_fx * fl1_fgb * tdx_yyy_yyy_0[j];

            tlx_yyyy_yyz_0[j] = pa_y[j] * tlx_yyy_yyz_0[j] + 1.5 * fl1_fx * tlx_yy_yyz_0[j] + fl1_fx * tlx_yyy_yz_0[j] -
                                0.5 * fl1_fx * tpz_yyy_yyz_0[j] - fl1_fx * fl1_fgb * tdz_yyy_yyz_0[j];

            tly_yyyy_yyz_0[j] = pa_y[j] * tly_yyy_yyz_0[j] + 1.5 * fl1_fx * tly_yy_yyz_0[j] + fl1_fx * tly_yyy_yz_0[j];

            tlz_yyyy_yyz_0[j] = pa_y[j] * tlz_yyy_yyz_0[j] + 1.5 * fl1_fx * tlz_yy_yyz_0[j] + fl1_fx * tlz_yyy_yz_0[j] +
                                0.5 * fl1_fx * tpx_yyy_yyz_0[j] + fl1_fx * fl1_fgb * tdx_yyy_yyz_0[j];

            tlx_yyyy_yzz_0[j] = pa_y[j] * tlx_yyy_yzz_0[j] + 1.5 * fl1_fx * tlx_yy_yzz_0[j] + 0.5 * fl1_fx * tlx_yyy_zz_0[j] -
                                0.5 * fl1_fx * tpz_yyy_yzz_0[j] - fl1_fx * fl1_fgb * tdz_yyy_yzz_0[j];

            tly_yyyy_yzz_0[j] = pa_y[j] * tly_yyy_yzz_0[j] + 1.5 * fl1_fx * tly_yy_yzz_0[j] + 0.5 * fl1_fx * tly_yyy_zz_0[j];

            tlz_yyyy_yzz_0[j] = pa_y[j] * tlz_yyy_yzz_0[j] + 1.5 * fl1_fx * tlz_yy_yzz_0[j] + 0.5 * fl1_fx * tlz_yyy_zz_0[j] +
                                0.5 * fl1_fx * tpx_yyy_yzz_0[j] + fl1_fx * fl1_fgb * tdx_yyy_yzz_0[j];

            tlx_yyyy_zzz_0[j] =
                pa_y[j] * tlx_yyy_zzz_0[j] + 1.5 * fl1_fx * tlx_yy_zzz_0[j] - 0.5 * fl1_fx * tpz_yyy_zzz_0[j] - fl1_fx * fl1_fgb * tdz_yyy_zzz_0[j];

            tly_yyyy_zzz_0[j] = pa_y[j] * tly_yyy_zzz_0[j] + 1.5 * fl1_fx * tly_yy_zzz_0[j];

            tlz_yyyy_zzz_0[j] =
                pa_y[j] * tlz_yyy_zzz_0[j] + 1.5 * fl1_fx * tlz_yy_zzz_0[j] + 0.5 * fl1_fx * tpx_yyy_zzz_0[j] + fl1_fx * fl1_fgb * tdx_yyy_zzz_0[j];

            tlx_yyyz_xxx_0[j] =
                pa_y[j] * tlx_yyz_xxx_0[j] + fl1_fx * tlx_yz_xxx_0[j] - 0.5 * fl1_fx * tpz_yyz_xxx_0[j] - fl1_fx * fl1_fgb * tdz_yyz_xxx_0[j];

            tly_yyyz_xxx_0[j] = pa_y[j] * tly_yyz_xxx_0[j] + fl1_fx * tly_yz_xxx_0[j];

            tlz_yyyz_xxx_0[j] =
                pa_y[j] * tlz_yyz_xxx_0[j] + fl1_fx * tlz_yz_xxx_0[j] + 0.5 * fl1_fx * tpx_yyz_xxx_0[j] + fl1_fx * fl1_fgb * tdx_yyz_xxx_0[j];

            tlx_yyyz_xxy_0[j] = pa_y[j] * tlx_yyz_xxy_0[j] + fl1_fx * tlx_yz_xxy_0[j] + 0.5 * fl1_fx * tlx_yyz_xx_0[j] -
                                0.5 * fl1_fx * tpz_yyz_xxy_0[j] - fl1_fx * fl1_fgb * tdz_yyz_xxy_0[j];

            tly_yyyz_xxy_0[j] = pa_y[j] * tly_yyz_xxy_0[j] + fl1_fx * tly_yz_xxy_0[j] + 0.5 * fl1_fx * tly_yyz_xx_0[j];

            tlz_yyyz_xxy_0[j] = pa_y[j] * tlz_yyz_xxy_0[j] + fl1_fx * tlz_yz_xxy_0[j] + 0.5 * fl1_fx * tlz_yyz_xx_0[j] +
                                0.5 * fl1_fx * tpx_yyz_xxy_0[j] + fl1_fx * fl1_fgb * tdx_yyz_xxy_0[j];

            tlx_yyyz_xxz_0[j] =
                pa_y[j] * tlx_yyz_xxz_0[j] + fl1_fx * tlx_yz_xxz_0[j] - 0.5 * fl1_fx * tpz_yyz_xxz_0[j] - fl1_fx * fl1_fgb * tdz_yyz_xxz_0[j];

            tly_yyyz_xxz_0[j] = pa_y[j] * tly_yyz_xxz_0[j] + fl1_fx * tly_yz_xxz_0[j];

            tlz_yyyz_xxz_0[j] =
                pa_y[j] * tlz_yyz_xxz_0[j] + fl1_fx * tlz_yz_xxz_0[j] + 0.5 * fl1_fx * tpx_yyz_xxz_0[j] + fl1_fx * fl1_fgb * tdx_yyz_xxz_0[j];

            tlx_yyyz_xyy_0[j] = pa_y[j] * tlx_yyz_xyy_0[j] + fl1_fx * tlx_yz_xyy_0[j] + fl1_fx * tlx_yyz_xy_0[j] - 0.5 * fl1_fx * tpz_yyz_xyy_0[j] -
                                fl1_fx * fl1_fgb * tdz_yyz_xyy_0[j];

            tly_yyyz_xyy_0[j] = pa_y[j] * tly_yyz_xyy_0[j] + fl1_fx * tly_yz_xyy_0[j] + fl1_fx * tly_yyz_xy_0[j];

            tlz_yyyz_xyy_0[j] = pa_y[j] * tlz_yyz_xyy_0[j] + fl1_fx * tlz_yz_xyy_0[j] + fl1_fx * tlz_yyz_xy_0[j] + 0.5 * fl1_fx * tpx_yyz_xyy_0[j] +
                                fl1_fx * fl1_fgb * tdx_yyz_xyy_0[j];

            tlx_yyyz_xyz_0[j] = pa_y[j] * tlx_yyz_xyz_0[j] + fl1_fx * tlx_yz_xyz_0[j] + 0.5 * fl1_fx * tlx_yyz_xz_0[j] -
                                0.5 * fl1_fx * tpz_yyz_xyz_0[j] - fl1_fx * fl1_fgb * tdz_yyz_xyz_0[j];

            tly_yyyz_xyz_0[j] = pa_y[j] * tly_yyz_xyz_0[j] + fl1_fx * tly_yz_xyz_0[j] + 0.5 * fl1_fx * tly_yyz_xz_0[j];

            tlz_yyyz_xyz_0[j] = pa_y[j] * tlz_yyz_xyz_0[j] + fl1_fx * tlz_yz_xyz_0[j] + 0.5 * fl1_fx * tlz_yyz_xz_0[j] +
                                0.5 * fl1_fx * tpx_yyz_xyz_0[j] + fl1_fx * fl1_fgb * tdx_yyz_xyz_0[j];

            tlx_yyyz_xzz_0[j] =
                pa_y[j] * tlx_yyz_xzz_0[j] + fl1_fx * tlx_yz_xzz_0[j] - 0.5 * fl1_fx * tpz_yyz_xzz_0[j] - fl1_fx * fl1_fgb * tdz_yyz_xzz_0[j];

            tly_yyyz_xzz_0[j] = pa_y[j] * tly_yyz_xzz_0[j] + fl1_fx * tly_yz_xzz_0[j];

            tlz_yyyz_xzz_0[j] =
                pa_y[j] * tlz_yyz_xzz_0[j] + fl1_fx * tlz_yz_xzz_0[j] + 0.5 * fl1_fx * tpx_yyz_xzz_0[j] + fl1_fx * fl1_fgb * tdx_yyz_xzz_0[j];

            tlx_yyyz_yyy_0[j] = pa_y[j] * tlx_yyz_yyy_0[j] + fl1_fx * tlx_yz_yyy_0[j] + 1.5 * fl1_fx * tlx_yyz_yy_0[j] -
                                0.5 * fl1_fx * tpz_yyz_yyy_0[j] - fl1_fx * fl1_fgb * tdz_yyz_yyy_0[j];

            tly_yyyz_yyy_0[j] = pa_y[j] * tly_yyz_yyy_0[j] + fl1_fx * tly_yz_yyy_0[j] + 1.5 * fl1_fx * tly_yyz_yy_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForGF_350_400(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_l_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_4_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tly_zz_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 58);

        auto tlz_zz_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 58);

        auto tlx_zz_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 59);

        auto tly_zz_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 59);

        auto tlz_zz_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 59);

        auto tlz_yyz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 45);

        auto tlx_yyz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 46);

        auto tly_yyz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 46);

        auto tlz_yyz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 46);

        auto tlx_yyz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 47);

        auto tly_yyz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 47);

        auto tlz_yyz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 47);

        auto tlx_yzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 48);

        auto tly_yzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 48);

        auto tlz_yzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 48);

        auto tlx_yzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 49);

        auto tly_yzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 49);

        auto tlz_yzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 49);

        auto tlx_yzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 50);

        auto tly_yzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 50);

        auto tlz_yzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 50);

        auto tlx_yzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 51);

        auto tly_yzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 51);

        auto tlz_yzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 51);

        auto tlx_yzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 52);

        auto tly_yzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 52);

        auto tlz_yzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 52);

        auto tlx_yzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 53);

        auto tly_yzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 53);

        auto tlz_yzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 53);

        auto tlx_zzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 54);

        auto tly_zzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 54);

        auto tlz_zzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 54);

        auto tlx_zzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 55);

        auto tpx_yyz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 76);

        auto tpx_yyz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 77);

        auto tpz_yyz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 77);

        auto tpx_yyz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 78);

        auto tpz_yyz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 78);

        auto tpx_yyz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 79);

        auto tpz_yyz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 79);

        auto tpx_yzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 80);

        auto tpz_yzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 80);

        auto tpx_yzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 81);

        auto tpz_yzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 81);

        auto tpx_yzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 82);

        auto tpz_yzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 82);

        auto tpx_yzz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 83);

        auto tpz_yzz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 83);

        auto tpx_yzz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 84);

        auto tpz_yzz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 84);

        auto tpx_yzz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 85);

        auto tpz_yzz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 85);

        auto tpx_yzz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 86);

        auto tpz_yzz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 86);

        auto tpx_yzz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 87);

        auto tpz_yzz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 87);

        auto tpx_yzz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 88);

        auto tpz_yzz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 88);

        auto tpx_yzz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 89);

        auto tpz_yzz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 89);

        auto tpx_zzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 90);

        auto tpz_zzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 90);

        auto tpx_zzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 91);

        auto tpz_zzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 91);

        auto tpx_zzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 92);

        auto tpz_zzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 92);

        auto tpz_zzz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 93);

        auto tdx_yyz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 76);

        auto tdx_yyz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 77);

        auto tdz_yyz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 77);

        auto tdx_yyz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 78);

        auto tdz_yyz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 78);

        auto tdx_yyz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 79);

        auto tdz_yyz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 79);

        auto tdx_yzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 80);

        auto tdz_yzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 80);

        auto tdx_yzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 81);

        auto tdz_yzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 81);

        auto tdx_yzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 82);

        auto tdz_yzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 82);

        auto tdx_yzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 83);

        auto tdz_yzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 83);

        auto tdx_yzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 84);

        auto tdz_yzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 84);

        auto tdx_yzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 85);

        auto tdz_yzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 85);

        auto tdx_yzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 86);

        auto tdz_yzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 86);

        auto tdx_yzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 87);

        auto tdz_yzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 87);

        auto tdx_yzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 88);

        auto tdz_yzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 88);

        auto tdx_yzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 89);

        auto tdz_yzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 89);

        auto tdx_zzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 90);

        auto tdz_zzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 90);

        auto tdx_zzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 91);

        auto tdz_zzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 91);

        auto tdx_zzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 92);

        auto tdz_zzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 92);

        auto tdz_zzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 93);

        // set up pointers to integrals

        auto tlz_yyyz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 116);

        auto tlx_yyyz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 117);

        auto tly_yyyz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 117);

        auto tlz_yyyz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 117);

        auto tlx_yyyz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 118);

        auto tly_yyyz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 118);

        auto tlz_yyyz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 118);

        auto tlx_yyyz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 119);

        auto tly_yyyz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 119);

        auto tlz_yyyz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 119);

        auto tlx_yyzz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 120);

        auto tly_yyzz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 120);

        auto tlz_yyzz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 120);

        auto tlx_yyzz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 121);

        auto tly_yyzz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 121);

        auto tlz_yyzz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 121);

        auto tlx_yyzz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 122);

        auto tly_yyzz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 122);

        auto tlz_yyzz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 122);

        auto tlx_yyzz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 123);

        auto tly_yyzz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 123);

        auto tlz_yyzz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 123);

        auto tlx_yyzz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 124);

        auto tly_yyzz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 124);

        auto tlz_yyzz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 124);

        auto tlx_yyzz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 125);

        auto tly_yyzz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 125);

        auto tlz_yyzz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 125);

        auto tlx_yyzz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 126);

        auto tly_yyzz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 126);

        auto tlz_yyzz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 126);

        auto tlx_yyzz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 127);

        auto tly_yyzz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 127);

        auto tlz_yyzz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 127);

        auto tlx_yyzz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 128);

        auto tly_yyzz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 128);

        auto tlz_yyzz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 128);

        auto tlx_yyzz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 129);

        auto tly_yyzz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 129);

        auto tlz_yyzz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 129);

        auto tlx_yzzz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 130);

        auto tly_yzzz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 130);

        auto tlz_yzzz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 130);

        auto tlx_yzzz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 131);

        auto tly_yzzz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 131);

        auto tlz_yzzz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 131);

        auto tlx_yzzz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 132);

        auto tly_yzzz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 132);

        auto tlz_yzzz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 132);

        auto tlx_yzzz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 133);

        // Batch of Integrals (350,400)

        #pragma omp simd aligned(fgb, fx, pa_y, tdx_yyz_yyy_0, tdx_yyz_yyz_0, tdx_yyz_yzz_0, \
                                     tdx_yyz_zzz_0, tdx_yzz_xxx_0, tdx_yzz_xxy_0, tdx_yzz_xxz_0, tdx_yzz_xyy_0, \
                                     tdx_yzz_xyz_0, tdx_yzz_xzz_0, tdx_yzz_yyy_0, tdx_yzz_yyz_0, tdx_yzz_yzz_0, \
                                     tdx_yzz_zzz_0, tdx_zzz_xxx_0, tdx_zzz_xxy_0, tdx_zzz_xxz_0, tdz_yyz_yyz_0, \
                                     tdz_yyz_yzz_0, tdz_yyz_zzz_0, tdz_yzz_xxx_0, tdz_yzz_xxy_0, tdz_yzz_xxz_0, \
                                     tdz_yzz_xyy_0, tdz_yzz_xyz_0, tdz_yzz_xzz_0, tdz_yzz_yyy_0, tdz_yzz_yyz_0, \
                                     tdz_yzz_yzz_0, tdz_yzz_zzz_0, tdz_zzz_xxx_0, tdz_zzz_xxy_0, tdz_zzz_xxz_0, \
                                     tdz_zzz_xyy_0, tlx_yyyz_yyz_0, tlx_yyyz_yzz_0, tlx_yyyz_zzz_0, tlx_yyz_yyz_0, \
                                     tlx_yyz_yz_0, tlx_yyz_yzz_0, tlx_yyz_zz_0, tlx_yyz_zzz_0, tlx_yyzz_xxx_0, \
                                     tlx_yyzz_xxy_0, tlx_yyzz_xxz_0, tlx_yyzz_xyy_0, tlx_yyzz_xyz_0, tlx_yyzz_xzz_0, \
                                     tlx_yyzz_yyy_0, tlx_yyzz_yyz_0, tlx_yyzz_yzz_0, tlx_yyzz_zzz_0, tlx_yz_yyz_0, \
                                     tlx_yz_yzz_0, tlx_yz_zzz_0, tlx_yzz_xx_0, tlx_yzz_xxx_0, tlx_yzz_xxy_0, \
                                     tlx_yzz_xxz_0, tlx_yzz_xy_0, tlx_yzz_xyy_0, tlx_yzz_xyz_0, tlx_yzz_xz_0, \
                                     tlx_yzz_xzz_0, tlx_yzz_yy_0, tlx_yzz_yyy_0, tlx_yzz_yyz_0, tlx_yzz_yz_0, \
                                     tlx_yzz_yzz_0, tlx_yzz_zz_0, tlx_yzz_zzz_0, tlx_yzzz_xxx_0, tlx_yzzz_xxy_0, \
                                     tlx_yzzz_xxz_0, tlx_yzzz_xyy_0, tlx_zz_xxx_0, tlx_zz_xxy_0, tlx_zz_xxz_0, \
                                     tlx_zz_xyy_0, tlx_zz_xyz_0, tlx_zz_xzz_0, tlx_zz_yyy_0, tlx_zz_yyz_0, tlx_zz_yzz_0, \
                                     tlx_zz_zzz_0, tlx_zzz_xx_0, tlx_zzz_xxx_0, tlx_zzz_xxy_0, tlx_zzz_xxz_0, \
                                     tlx_zzz_xy_0, tlx_zzz_xyy_0, tly_yyyz_yyz_0, tly_yyyz_yzz_0, tly_yyyz_zzz_0, \
                                     tly_yyz_yyz_0, tly_yyz_yz_0, tly_yyz_yzz_0, tly_yyz_zz_0, tly_yyz_zzz_0, \
                                     tly_yyzz_xxx_0, tly_yyzz_xxy_0, tly_yyzz_xxz_0, tly_yyzz_xyy_0, tly_yyzz_xyz_0, \
                                     tly_yyzz_xzz_0, tly_yyzz_yyy_0, tly_yyzz_yyz_0, tly_yyzz_yzz_0, tly_yyzz_zzz_0, \
                                     tly_yz_yyz_0, tly_yz_yzz_0, tly_yz_zzz_0, tly_yzz_xx_0, tly_yzz_xxx_0, \
                                     tly_yzz_xxy_0, tly_yzz_xxz_0, tly_yzz_xy_0, tly_yzz_xyy_0, tly_yzz_xyz_0, \
                                     tly_yzz_xz_0, tly_yzz_xzz_0, tly_yzz_yy_0, tly_yzz_yyy_0, tly_yzz_yyz_0, \
                                     tly_yzz_yz_0, tly_yzz_yzz_0, tly_yzz_zz_0, tly_yzz_zzz_0, tly_yzzz_xxx_0, \
                                     tly_yzzz_xxy_0, tly_yzzz_xxz_0, tly_zz_xxx_0, tly_zz_xxy_0, tly_zz_xxz_0, \
                                     tly_zz_xyy_0, tly_zz_xyz_0, tly_zz_xzz_0, tly_zz_yyy_0, tly_zz_yyz_0, tly_zz_yzz_0, \
                                     tly_zz_zzz_0, tly_zzz_xx_0, tly_zzz_xxx_0, tly_zzz_xxy_0, tly_zzz_xxz_0, \
                                     tlz_yyyz_yyy_0, tlz_yyyz_yyz_0, tlz_yyyz_yzz_0, tlz_yyyz_zzz_0, tlz_yyz_yy_0, \
                                     tlz_yyz_yyy_0, tlz_yyz_yyz_0, tlz_yyz_yz_0, tlz_yyz_yzz_0, tlz_yyz_zz_0, \
                                     tlz_yyz_zzz_0, tlz_yyzz_xxx_0, tlz_yyzz_xxy_0, tlz_yyzz_xxz_0, tlz_yyzz_xyy_0, \
                                     tlz_yyzz_xyz_0, tlz_yyzz_xzz_0, tlz_yyzz_yyy_0, tlz_yyzz_yyz_0, tlz_yyzz_yzz_0, \
                                     tlz_yyzz_zzz_0, tlz_yz_yyy_0, tlz_yz_yyz_0, tlz_yz_yzz_0, tlz_yz_zzz_0, tlz_yzz_xx_0, \
                                     tlz_yzz_xxx_0, tlz_yzz_xxy_0, tlz_yzz_xxz_0, tlz_yzz_xy_0, tlz_yzz_xyy_0, \
                                     tlz_yzz_xyz_0, tlz_yzz_xz_0, tlz_yzz_xzz_0, tlz_yzz_yy_0, tlz_yzz_yyy_0, \
                                     tlz_yzz_yyz_0, tlz_yzz_yz_0, tlz_yzz_yzz_0, tlz_yzz_zz_0, tlz_yzz_zzz_0, \
                                     tlz_yzzz_xxx_0, tlz_yzzz_xxy_0, tlz_yzzz_xxz_0, tlz_zz_xxx_0, tlz_zz_xxy_0, \
                                     tlz_zz_xxz_0, tlz_zz_xyy_0, tlz_zz_xyz_0, tlz_zz_xzz_0, tlz_zz_yyy_0, tlz_zz_yyz_0, \
                                     tlz_zz_yzz_0, tlz_zz_zzz_0, tlz_zzz_xx_0, tlz_zzz_xxx_0, tlz_zzz_xxy_0, \
                                     tlz_zzz_xxz_0, tpx_yyz_yyy_0, tpx_yyz_yyz_0, tpx_yyz_yzz_0, tpx_yyz_zzz_0, \
                                     tpx_yzz_xxx_0, tpx_yzz_xxy_0, tpx_yzz_xxz_0, tpx_yzz_xyy_0, tpx_yzz_xyz_0, \
                                     tpx_yzz_xzz_0, tpx_yzz_yyy_0, tpx_yzz_yyz_0, tpx_yzz_yzz_0, tpx_yzz_zzz_0, \
                                     tpx_zzz_xxx_0, tpx_zzz_xxy_0, tpx_zzz_xxz_0, tpz_yyz_yyz_0, tpz_yyz_yzz_0, \
                                     tpz_yyz_zzz_0, tpz_yzz_xxx_0, tpz_yzz_xxy_0, tpz_yzz_xxz_0, tpz_yzz_xyy_0, \
                                     tpz_yzz_xyz_0, tpz_yzz_xzz_0, tpz_yzz_yyy_0, tpz_yzz_yyz_0, tpz_yzz_yzz_0, \
                                     tpz_yzz_zzz_0, tpz_zzz_xxx_0, tpz_zzz_xxy_0, tpz_zzz_xxz_0, tpz_zzz_xyy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlz_yyyz_yyy_0[j] = pa_y[j] * tlz_yyz_yyy_0[j] + fl1_fx * tlz_yz_yyy_0[j] + 1.5 * fl1_fx * tlz_yyz_yy_0[j] +
                                0.5 * fl1_fx * tpx_yyz_yyy_0[j] + fl1_fx * fl1_fgb * tdx_yyz_yyy_0[j];

            tlx_yyyz_yyz_0[j] = pa_y[j] * tlx_yyz_yyz_0[j] + fl1_fx * tlx_yz_yyz_0[j] + fl1_fx * tlx_yyz_yz_0[j] - 0.5 * fl1_fx * tpz_yyz_yyz_0[j] -
                                fl1_fx * fl1_fgb * tdz_yyz_yyz_0[j];

            tly_yyyz_yyz_0[j] = pa_y[j] * tly_yyz_yyz_0[j] + fl1_fx * tly_yz_yyz_0[j] + fl1_fx * tly_yyz_yz_0[j];

            tlz_yyyz_yyz_0[j] = pa_y[j] * tlz_yyz_yyz_0[j] + fl1_fx * tlz_yz_yyz_0[j] + fl1_fx * tlz_yyz_yz_0[j] + 0.5 * fl1_fx * tpx_yyz_yyz_0[j] +
                                fl1_fx * fl1_fgb * tdx_yyz_yyz_0[j];

            tlx_yyyz_yzz_0[j] = pa_y[j] * tlx_yyz_yzz_0[j] + fl1_fx * tlx_yz_yzz_0[j] + 0.5 * fl1_fx * tlx_yyz_zz_0[j] -
                                0.5 * fl1_fx * tpz_yyz_yzz_0[j] - fl1_fx * fl1_fgb * tdz_yyz_yzz_0[j];

            tly_yyyz_yzz_0[j] = pa_y[j] * tly_yyz_yzz_0[j] + fl1_fx * tly_yz_yzz_0[j] + 0.5 * fl1_fx * tly_yyz_zz_0[j];

            tlz_yyyz_yzz_0[j] = pa_y[j] * tlz_yyz_yzz_0[j] + fl1_fx * tlz_yz_yzz_0[j] + 0.5 * fl1_fx * tlz_yyz_zz_0[j] +
                                0.5 * fl1_fx * tpx_yyz_yzz_0[j] + fl1_fx * fl1_fgb * tdx_yyz_yzz_0[j];

            tlx_yyyz_zzz_0[j] =
                pa_y[j] * tlx_yyz_zzz_0[j] + fl1_fx * tlx_yz_zzz_0[j] - 0.5 * fl1_fx * tpz_yyz_zzz_0[j] - fl1_fx * fl1_fgb * tdz_yyz_zzz_0[j];

            tly_yyyz_zzz_0[j] = pa_y[j] * tly_yyz_zzz_0[j] + fl1_fx * tly_yz_zzz_0[j];

            tlz_yyyz_zzz_0[j] =
                pa_y[j] * tlz_yyz_zzz_0[j] + fl1_fx * tlz_yz_zzz_0[j] + 0.5 * fl1_fx * tpx_yyz_zzz_0[j] + fl1_fx * fl1_fgb * tdx_yyz_zzz_0[j];

            tlx_yyzz_xxx_0[j] =
                pa_y[j] * tlx_yzz_xxx_0[j] + 0.5 * fl1_fx * tlx_zz_xxx_0[j] - 0.5 * fl1_fx * tpz_yzz_xxx_0[j] - fl1_fx * fl1_fgb * tdz_yzz_xxx_0[j];

            tly_yyzz_xxx_0[j] = pa_y[j] * tly_yzz_xxx_0[j] + 0.5 * fl1_fx * tly_zz_xxx_0[j];

            tlz_yyzz_xxx_0[j] =
                pa_y[j] * tlz_yzz_xxx_0[j] + 0.5 * fl1_fx * tlz_zz_xxx_0[j] + 0.5 * fl1_fx * tpx_yzz_xxx_0[j] + fl1_fx * fl1_fgb * tdx_yzz_xxx_0[j];

            tlx_yyzz_xxy_0[j] = pa_y[j] * tlx_yzz_xxy_0[j] + 0.5 * fl1_fx * tlx_zz_xxy_0[j] + 0.5 * fl1_fx * tlx_yzz_xx_0[j] -
                                0.5 * fl1_fx * tpz_yzz_xxy_0[j] - fl1_fx * fl1_fgb * tdz_yzz_xxy_0[j];

            tly_yyzz_xxy_0[j] = pa_y[j] * tly_yzz_xxy_0[j] + 0.5 * fl1_fx * tly_zz_xxy_0[j] + 0.5 * fl1_fx * tly_yzz_xx_0[j];

            tlz_yyzz_xxy_0[j] = pa_y[j] * tlz_yzz_xxy_0[j] + 0.5 * fl1_fx * tlz_zz_xxy_0[j] + 0.5 * fl1_fx * tlz_yzz_xx_0[j] +
                                0.5 * fl1_fx * tpx_yzz_xxy_0[j] + fl1_fx * fl1_fgb * tdx_yzz_xxy_0[j];

            tlx_yyzz_xxz_0[j] =
                pa_y[j] * tlx_yzz_xxz_0[j] + 0.5 * fl1_fx * tlx_zz_xxz_0[j] - 0.5 * fl1_fx * tpz_yzz_xxz_0[j] - fl1_fx * fl1_fgb * tdz_yzz_xxz_0[j];

            tly_yyzz_xxz_0[j] = pa_y[j] * tly_yzz_xxz_0[j] + 0.5 * fl1_fx * tly_zz_xxz_0[j];

            tlz_yyzz_xxz_0[j] =
                pa_y[j] * tlz_yzz_xxz_0[j] + 0.5 * fl1_fx * tlz_zz_xxz_0[j] + 0.5 * fl1_fx * tpx_yzz_xxz_0[j] + fl1_fx * fl1_fgb * tdx_yzz_xxz_0[j];

            tlx_yyzz_xyy_0[j] = pa_y[j] * tlx_yzz_xyy_0[j] + 0.5 * fl1_fx * tlx_zz_xyy_0[j] + fl1_fx * tlx_yzz_xy_0[j] -
                                0.5 * fl1_fx * tpz_yzz_xyy_0[j] - fl1_fx * fl1_fgb * tdz_yzz_xyy_0[j];

            tly_yyzz_xyy_0[j] = pa_y[j] * tly_yzz_xyy_0[j] + 0.5 * fl1_fx * tly_zz_xyy_0[j] + fl1_fx * tly_yzz_xy_0[j];

            tlz_yyzz_xyy_0[j] = pa_y[j] * tlz_yzz_xyy_0[j] + 0.5 * fl1_fx * tlz_zz_xyy_0[j] + fl1_fx * tlz_yzz_xy_0[j] +
                                0.5 * fl1_fx * tpx_yzz_xyy_0[j] + fl1_fx * fl1_fgb * tdx_yzz_xyy_0[j];

            tlx_yyzz_xyz_0[j] = pa_y[j] * tlx_yzz_xyz_0[j] + 0.5 * fl1_fx * tlx_zz_xyz_0[j] + 0.5 * fl1_fx * tlx_yzz_xz_0[j] -
                                0.5 * fl1_fx * tpz_yzz_xyz_0[j] - fl1_fx * fl1_fgb * tdz_yzz_xyz_0[j];

            tly_yyzz_xyz_0[j] = pa_y[j] * tly_yzz_xyz_0[j] + 0.5 * fl1_fx * tly_zz_xyz_0[j] + 0.5 * fl1_fx * tly_yzz_xz_0[j];

            tlz_yyzz_xyz_0[j] = pa_y[j] * tlz_yzz_xyz_0[j] + 0.5 * fl1_fx * tlz_zz_xyz_0[j] + 0.5 * fl1_fx * tlz_yzz_xz_0[j] +
                                0.5 * fl1_fx * tpx_yzz_xyz_0[j] + fl1_fx * fl1_fgb * tdx_yzz_xyz_0[j];

            tlx_yyzz_xzz_0[j] =
                pa_y[j] * tlx_yzz_xzz_0[j] + 0.5 * fl1_fx * tlx_zz_xzz_0[j] - 0.5 * fl1_fx * tpz_yzz_xzz_0[j] - fl1_fx * fl1_fgb * tdz_yzz_xzz_0[j];

            tly_yyzz_xzz_0[j] = pa_y[j] * tly_yzz_xzz_0[j] + 0.5 * fl1_fx * tly_zz_xzz_0[j];

            tlz_yyzz_xzz_0[j] =
                pa_y[j] * tlz_yzz_xzz_0[j] + 0.5 * fl1_fx * tlz_zz_xzz_0[j] + 0.5 * fl1_fx * tpx_yzz_xzz_0[j] + fl1_fx * fl1_fgb * tdx_yzz_xzz_0[j];

            tlx_yyzz_yyy_0[j] = pa_y[j] * tlx_yzz_yyy_0[j] + 0.5 * fl1_fx * tlx_zz_yyy_0[j] + 1.5 * fl1_fx * tlx_yzz_yy_0[j] -
                                0.5 * fl1_fx * tpz_yzz_yyy_0[j] - fl1_fx * fl1_fgb * tdz_yzz_yyy_0[j];

            tly_yyzz_yyy_0[j] = pa_y[j] * tly_yzz_yyy_0[j] + 0.5 * fl1_fx * tly_zz_yyy_0[j] + 1.5 * fl1_fx * tly_yzz_yy_0[j];

            tlz_yyzz_yyy_0[j] = pa_y[j] * tlz_yzz_yyy_0[j] + 0.5 * fl1_fx * tlz_zz_yyy_0[j] + 1.5 * fl1_fx * tlz_yzz_yy_0[j] +
                                0.5 * fl1_fx * tpx_yzz_yyy_0[j] + fl1_fx * fl1_fgb * tdx_yzz_yyy_0[j];

            tlx_yyzz_yyz_0[j] = pa_y[j] * tlx_yzz_yyz_0[j] + 0.5 * fl1_fx * tlx_zz_yyz_0[j] + fl1_fx * tlx_yzz_yz_0[j] -
                                0.5 * fl1_fx * tpz_yzz_yyz_0[j] - fl1_fx * fl1_fgb * tdz_yzz_yyz_0[j];

            tly_yyzz_yyz_0[j] = pa_y[j] * tly_yzz_yyz_0[j] + 0.5 * fl1_fx * tly_zz_yyz_0[j] + fl1_fx * tly_yzz_yz_0[j];

            tlz_yyzz_yyz_0[j] = pa_y[j] * tlz_yzz_yyz_0[j] + 0.5 * fl1_fx * tlz_zz_yyz_0[j] + fl1_fx * tlz_yzz_yz_0[j] +
                                0.5 * fl1_fx * tpx_yzz_yyz_0[j] + fl1_fx * fl1_fgb * tdx_yzz_yyz_0[j];

            tlx_yyzz_yzz_0[j] = pa_y[j] * tlx_yzz_yzz_0[j] + 0.5 * fl1_fx * tlx_zz_yzz_0[j] + 0.5 * fl1_fx * tlx_yzz_zz_0[j] -
                                0.5 * fl1_fx * tpz_yzz_yzz_0[j] - fl1_fx * fl1_fgb * tdz_yzz_yzz_0[j];

            tly_yyzz_yzz_0[j] = pa_y[j] * tly_yzz_yzz_0[j] + 0.5 * fl1_fx * tly_zz_yzz_0[j] + 0.5 * fl1_fx * tly_yzz_zz_0[j];

            tlz_yyzz_yzz_0[j] = pa_y[j] * tlz_yzz_yzz_0[j] + 0.5 * fl1_fx * tlz_zz_yzz_0[j] + 0.5 * fl1_fx * tlz_yzz_zz_0[j] +
                                0.5 * fl1_fx * tpx_yzz_yzz_0[j] + fl1_fx * fl1_fgb * tdx_yzz_yzz_0[j];

            tlx_yyzz_zzz_0[j] =
                pa_y[j] * tlx_yzz_zzz_0[j] + 0.5 * fl1_fx * tlx_zz_zzz_0[j] - 0.5 * fl1_fx * tpz_yzz_zzz_0[j] - fl1_fx * fl1_fgb * tdz_yzz_zzz_0[j];

            tly_yyzz_zzz_0[j] = pa_y[j] * tly_yzz_zzz_0[j] + 0.5 * fl1_fx * tly_zz_zzz_0[j];

            tlz_yyzz_zzz_0[j] =
                pa_y[j] * tlz_yzz_zzz_0[j] + 0.5 * fl1_fx * tlz_zz_zzz_0[j] + 0.5 * fl1_fx * tpx_yzz_zzz_0[j] + fl1_fx * fl1_fgb * tdx_yzz_zzz_0[j];

            tlx_yzzz_xxx_0[j] = pa_y[j] * tlx_zzz_xxx_0[j] - 0.5 * fl1_fx * tpz_zzz_xxx_0[j] - fl1_fx * fl1_fgb * tdz_zzz_xxx_0[j];

            tly_yzzz_xxx_0[j] = pa_y[j] * tly_zzz_xxx_0[j];

            tlz_yzzz_xxx_0[j] = pa_y[j] * tlz_zzz_xxx_0[j] + 0.5 * fl1_fx * tpx_zzz_xxx_0[j] + fl1_fx * fl1_fgb * tdx_zzz_xxx_0[j];

            tlx_yzzz_xxy_0[j] =
                pa_y[j] * tlx_zzz_xxy_0[j] + 0.5 * fl1_fx * tlx_zzz_xx_0[j] - 0.5 * fl1_fx * tpz_zzz_xxy_0[j] - fl1_fx * fl1_fgb * tdz_zzz_xxy_0[j];

            tly_yzzz_xxy_0[j] = pa_y[j] * tly_zzz_xxy_0[j] + 0.5 * fl1_fx * tly_zzz_xx_0[j];

            tlz_yzzz_xxy_0[j] =
                pa_y[j] * tlz_zzz_xxy_0[j] + 0.5 * fl1_fx * tlz_zzz_xx_0[j] + 0.5 * fl1_fx * tpx_zzz_xxy_0[j] + fl1_fx * fl1_fgb * tdx_zzz_xxy_0[j];

            tlx_yzzz_xxz_0[j] = pa_y[j] * tlx_zzz_xxz_0[j] - 0.5 * fl1_fx * tpz_zzz_xxz_0[j] - fl1_fx * fl1_fgb * tdz_zzz_xxz_0[j];

            tly_yzzz_xxz_0[j] = pa_y[j] * tly_zzz_xxz_0[j];

            tlz_yzzz_xxz_0[j] = pa_y[j] * tlz_zzz_xxz_0[j] + 0.5 * fl1_fx * tpx_zzz_xxz_0[j] + fl1_fx * fl1_fgb * tdx_zzz_xxz_0[j];

            tlx_yzzz_xyy_0[j] =
                pa_y[j] * tlx_zzz_xyy_0[j] + fl1_fx * tlx_zzz_xy_0[j] - 0.5 * fl1_fx * tpz_zzz_xyy_0[j] - fl1_fx * fl1_fgb * tdz_zzz_xyy_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForGF_400_450(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_l_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_4_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tlx_zzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 54);

        auto tly_zzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 54);

        auto tlz_zzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 54);

        auto tlx_zzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 55);

        auto tly_zzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 55);

        auto tlz_zzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 55);

        auto tlx_zzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 56);

        auto tly_zzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 56);

        auto tlz_zzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 56);

        auto tlx_zzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 57);

        auto tly_zzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 57);

        auto tlz_zzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 57);

        auto tlx_zzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 58);

        auto tly_zzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 58);

        auto tlz_zzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 58);

        auto tlx_zzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 59);

        auto tly_zzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 59);

        auto tlz_zzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 59);

        auto tpx_zzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 90);

        auto tpy_zzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 90);

        auto tpx_zzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 91);

        auto tpy_zzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 91);

        auto tpx_zzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 92);

        auto tpy_zzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 92);

        auto tpx_zzz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 93);

        auto tpy_zzz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 93);

        auto tpx_zzz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 94);

        auto tpy_zzz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 94);

        auto tpz_zzz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 94);

        auto tpx_zzz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 95);

        auto tpy_zzz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 95);

        auto tpz_zzz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 95);

        auto tpx_zzz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 96);

        auto tpy_zzz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 96);

        auto tpz_zzz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 96);

        auto tpx_zzz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 97);

        auto tpy_zzz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 97);

        auto tpz_zzz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 97);

        auto tpx_zzz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 98);

        auto tpy_zzz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 98);

        auto tpz_zzz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 98);

        auto tpx_zzz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 99);

        auto tpy_zzz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 99);

        auto tpz_zzz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 99);

        auto tdx_zzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 90);

        auto tdy_zzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 90);

        auto tdx_zzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 91);

        auto tdy_zzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 91);

        auto tdx_zzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 92);

        auto tdy_zzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 92);

        auto tdx_zzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 93);

        auto tdy_zzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 93);

        auto tdx_zzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 94);

        auto tdy_zzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 94);

        auto tdz_zzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 94);

        auto tdx_zzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 95);

        auto tdy_zzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 95);

        auto tdz_zzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 95);

        auto tdx_zzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 96);

        auto tdy_zzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 96);

        auto tdz_zzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 96);

        auto tdx_zzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 97);

        auto tdy_zzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 97);

        auto tdz_zzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 97);

        auto tdx_zzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 98);

        auto tdy_zzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 98);

        auto tdz_zzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 98);

        auto tdx_zzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 99);

        auto tdy_zzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 99);

        auto tdz_zzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 99);

        // set up pointers to integrals

        auto tly_yzzz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 133);

        auto tlz_yzzz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 133);

        auto tlx_yzzz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 134);

        auto tly_yzzz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 134);

        auto tlz_yzzz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 134);

        auto tlx_yzzz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 135);

        auto tly_yzzz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 135);

        auto tlz_yzzz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 135);

        auto tlx_yzzz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 136);

        auto tly_yzzz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 136);

        auto tlz_yzzz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 136);

        auto tlx_yzzz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 137);

        auto tly_yzzz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 137);

        auto tlz_yzzz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 137);

        auto tlx_yzzz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 138);

        auto tly_yzzz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 138);

        auto tlz_yzzz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 138);

        auto tlx_yzzz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 139);

        auto tly_yzzz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 139);

        auto tlz_yzzz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 139);

        auto tlx_zzzz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 140);

        auto tly_zzzz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 140);

        auto tlz_zzzz_xxx_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 140);

        auto tlx_zzzz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 141);

        auto tly_zzzz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 141);

        auto tlz_zzzz_xxy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 141);

        auto tlx_zzzz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 142);

        auto tly_zzzz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 142);

        auto tlz_zzzz_xxz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 142);

        auto tlx_zzzz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 143);

        auto tly_zzzz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 143);

        auto tlz_zzzz_xyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 143);

        auto tlx_zzzz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 144);

        auto tly_zzzz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 144);

        auto tlz_zzzz_xyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 144);

        auto tlx_zzzz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 145);

        auto tly_zzzz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 145);

        auto tlz_zzzz_xzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 145);

        auto tlx_zzzz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 146);

        auto tly_zzzz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 146);

        auto tlz_zzzz_yyy_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 146);

        auto tlx_zzzz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 147);

        auto tly_zzzz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 147);

        auto tlz_zzzz_yyz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 147);

        auto tlx_zzzz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 148);

        auto tly_zzzz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 148);

        auto tlz_zzzz_yzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 148);

        auto tlx_zzzz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * idx + 149);

        auto tly_zzzz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 150 * bdim + 150 * idx + 149);

        auto tlz_zzzz_zzz_0 = primBuffer.data(pidx_l_4_3_m0 + 300 * bdim + 150 * idx + 149);

        // Batch of Integrals (400,450)

        #pragma omp simd aligned(fgb, fx, pa_y, pa_z, tdx_zzz_xxx_0, tdx_zzz_xxy_0, tdx_zzz_xxz_0, \
                                     tdx_zzz_xyy_0, tdx_zzz_xyz_0, tdx_zzz_xzz_0, tdx_zzz_yyy_0, tdx_zzz_yyz_0, \
                                     tdx_zzz_yzz_0, tdx_zzz_zzz_0, tdy_zzz_xxx_0, tdy_zzz_xxy_0, tdy_zzz_xxz_0, \
                                     tdy_zzz_xyy_0, tdy_zzz_xyz_0, tdy_zzz_xzz_0, tdy_zzz_yyy_0, tdy_zzz_yyz_0, \
                                     tdy_zzz_yzz_0, tdy_zzz_zzz_0, tdz_zzz_xyz_0, tdz_zzz_xzz_0, tdz_zzz_yyy_0, \
                                     tdz_zzz_yyz_0, tdz_zzz_yzz_0, tdz_zzz_zzz_0, tlx_yzzz_xyz_0, tlx_yzzz_xzz_0, \
                                     tlx_yzzz_yyy_0, tlx_yzzz_yyz_0, tlx_yzzz_yzz_0, tlx_yzzz_zzz_0, tlx_zz_xxx_0, \
                                     tlx_zz_xxy_0, tlx_zz_xxz_0, tlx_zz_xyy_0, tlx_zz_xyz_0, tlx_zz_xzz_0, tlx_zz_yyy_0, \
                                     tlx_zz_yyz_0, tlx_zz_yzz_0, tlx_zz_zzz_0, tlx_zzz_xx_0, tlx_zzz_xxx_0, \
                                     tlx_zzz_xxy_0, tlx_zzz_xxz_0, tlx_zzz_xy_0, tlx_zzz_xyy_0, tlx_zzz_xyz_0, \
                                     tlx_zzz_xz_0, tlx_zzz_xzz_0, tlx_zzz_yy_0, tlx_zzz_yyy_0, tlx_zzz_yyz_0, \
                                     tlx_zzz_yz_0, tlx_zzz_yzz_0, tlx_zzz_zz_0, tlx_zzz_zzz_0, tlx_zzzz_xxx_0, \
                                     tlx_zzzz_xxy_0, tlx_zzzz_xxz_0, tlx_zzzz_xyy_0, tlx_zzzz_xyz_0, tlx_zzzz_xzz_0, \
                                     tlx_zzzz_yyy_0, tlx_zzzz_yyz_0, tlx_zzzz_yzz_0, tlx_zzzz_zzz_0, tly_yzzz_xyy_0, \
                                     tly_yzzz_xyz_0, tly_yzzz_xzz_0, tly_yzzz_yyy_0, tly_yzzz_yyz_0, tly_yzzz_yzz_0, \
                                     tly_yzzz_zzz_0, tly_zz_xxx_0, tly_zz_xxy_0, tly_zz_xxz_0, tly_zz_xyy_0, tly_zz_xyz_0, \
                                     tly_zz_xzz_0, tly_zz_yyy_0, tly_zz_yyz_0, tly_zz_yzz_0, tly_zz_zzz_0, tly_zzz_xx_0, \
                                     tly_zzz_xxx_0, tly_zzz_xxy_0, tly_zzz_xxz_0, tly_zzz_xy_0, tly_zzz_xyy_0, \
                                     tly_zzz_xyz_0, tly_zzz_xz_0, tly_zzz_xzz_0, tly_zzz_yy_0, tly_zzz_yyy_0, \
                                     tly_zzz_yyz_0, tly_zzz_yz_0, tly_zzz_yzz_0, tly_zzz_zz_0, tly_zzz_zzz_0, \
                                     tly_zzzz_xxx_0, tly_zzzz_xxy_0, tly_zzzz_xxz_0, tly_zzzz_xyy_0, tly_zzzz_xyz_0, \
                                     tly_zzzz_xzz_0, tly_zzzz_yyy_0, tly_zzzz_yyz_0, tly_zzzz_yzz_0, tly_zzzz_zzz_0, \
                                     tlz_yzzz_xyy_0, tlz_yzzz_xyz_0, tlz_yzzz_xzz_0, tlz_yzzz_yyy_0, tlz_yzzz_yyz_0, \
                                     tlz_yzzz_yzz_0, tlz_yzzz_zzz_0, tlz_zz_xxx_0, tlz_zz_xxy_0, tlz_zz_xxz_0, \
                                     tlz_zz_xyy_0, tlz_zz_xyz_0, tlz_zz_xzz_0, tlz_zz_yyy_0, tlz_zz_yyz_0, tlz_zz_yzz_0, \
                                     tlz_zz_zzz_0, tlz_zzz_xx_0, tlz_zzz_xxx_0, tlz_zzz_xxy_0, tlz_zzz_xxz_0, \
                                     tlz_zzz_xy_0, tlz_zzz_xyy_0, tlz_zzz_xyz_0, tlz_zzz_xz_0, tlz_zzz_xzz_0, \
                                     tlz_zzz_yy_0, tlz_zzz_yyy_0, tlz_zzz_yyz_0, tlz_zzz_yz_0, tlz_zzz_yzz_0, \
                                     tlz_zzz_zz_0, tlz_zzz_zzz_0, tlz_zzzz_xxx_0, tlz_zzzz_xxy_0, tlz_zzzz_xxz_0, \
                                     tlz_zzzz_xyy_0, tlz_zzzz_xyz_0, tlz_zzzz_xzz_0, tlz_zzzz_yyy_0, tlz_zzzz_yyz_0, \
                                     tlz_zzzz_yzz_0, tlz_zzzz_zzz_0, tpx_zzz_xxx_0, tpx_zzz_xxy_0, tpx_zzz_xxz_0, \
                                     tpx_zzz_xyy_0, tpx_zzz_xyz_0, tpx_zzz_xzz_0, tpx_zzz_yyy_0, tpx_zzz_yyz_0, \
                                     tpx_zzz_yzz_0, tpx_zzz_zzz_0, tpy_zzz_xxx_0, tpy_zzz_xxy_0, tpy_zzz_xxz_0, \
                                     tpy_zzz_xyy_0, tpy_zzz_xyz_0, tpy_zzz_xzz_0, tpy_zzz_yyy_0, tpy_zzz_yyz_0, \
                                     tpy_zzz_yzz_0, tpy_zzz_zzz_0, tpz_zzz_xyz_0, tpz_zzz_xzz_0, tpz_zzz_yyy_0, \
                                     tpz_zzz_yyz_0, tpz_zzz_yzz_0, tpz_zzz_zzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tly_yzzz_xyy_0[j] = pa_y[j] * tly_zzz_xyy_0[j] + fl1_fx * tly_zzz_xy_0[j];

            tlz_yzzz_xyy_0[j] =
                pa_y[j] * tlz_zzz_xyy_0[j] + fl1_fx * tlz_zzz_xy_0[j] + 0.5 * fl1_fx * tpx_zzz_xyy_0[j] + fl1_fx * fl1_fgb * tdx_zzz_xyy_0[j];

            tlx_yzzz_xyz_0[j] =
                pa_y[j] * tlx_zzz_xyz_0[j] + 0.5 * fl1_fx * tlx_zzz_xz_0[j] - 0.5 * fl1_fx * tpz_zzz_xyz_0[j] - fl1_fx * fl1_fgb * tdz_zzz_xyz_0[j];

            tly_yzzz_xyz_0[j] = pa_y[j] * tly_zzz_xyz_0[j] + 0.5 * fl1_fx * tly_zzz_xz_0[j];

            tlz_yzzz_xyz_0[j] =
                pa_y[j] * tlz_zzz_xyz_0[j] + 0.5 * fl1_fx * tlz_zzz_xz_0[j] + 0.5 * fl1_fx * tpx_zzz_xyz_0[j] + fl1_fx * fl1_fgb * tdx_zzz_xyz_0[j];

            tlx_yzzz_xzz_0[j] = pa_y[j] * tlx_zzz_xzz_0[j] - 0.5 * fl1_fx * tpz_zzz_xzz_0[j] - fl1_fx * fl1_fgb * tdz_zzz_xzz_0[j];

            tly_yzzz_xzz_0[j] = pa_y[j] * tly_zzz_xzz_0[j];

            tlz_yzzz_xzz_0[j] = pa_y[j] * tlz_zzz_xzz_0[j] + 0.5 * fl1_fx * tpx_zzz_xzz_0[j] + fl1_fx * fl1_fgb * tdx_zzz_xzz_0[j];

            tlx_yzzz_yyy_0[j] =
                pa_y[j] * tlx_zzz_yyy_0[j] + 1.5 * fl1_fx * tlx_zzz_yy_0[j] - 0.5 * fl1_fx * tpz_zzz_yyy_0[j] - fl1_fx * fl1_fgb * tdz_zzz_yyy_0[j];

            tly_yzzz_yyy_0[j] = pa_y[j] * tly_zzz_yyy_0[j] + 1.5 * fl1_fx * tly_zzz_yy_0[j];

            tlz_yzzz_yyy_0[j] =
                pa_y[j] * tlz_zzz_yyy_0[j] + 1.5 * fl1_fx * tlz_zzz_yy_0[j] + 0.5 * fl1_fx * tpx_zzz_yyy_0[j] + fl1_fx * fl1_fgb * tdx_zzz_yyy_0[j];

            tlx_yzzz_yyz_0[j] =
                pa_y[j] * tlx_zzz_yyz_0[j] + fl1_fx * tlx_zzz_yz_0[j] - 0.5 * fl1_fx * tpz_zzz_yyz_0[j] - fl1_fx * fl1_fgb * tdz_zzz_yyz_0[j];

            tly_yzzz_yyz_0[j] = pa_y[j] * tly_zzz_yyz_0[j] + fl1_fx * tly_zzz_yz_0[j];

            tlz_yzzz_yyz_0[j] =
                pa_y[j] * tlz_zzz_yyz_0[j] + fl1_fx * tlz_zzz_yz_0[j] + 0.5 * fl1_fx * tpx_zzz_yyz_0[j] + fl1_fx * fl1_fgb * tdx_zzz_yyz_0[j];

            tlx_yzzz_yzz_0[j] =
                pa_y[j] * tlx_zzz_yzz_0[j] + 0.5 * fl1_fx * tlx_zzz_zz_0[j] - 0.5 * fl1_fx * tpz_zzz_yzz_0[j] - fl1_fx * fl1_fgb * tdz_zzz_yzz_0[j];

            tly_yzzz_yzz_0[j] = pa_y[j] * tly_zzz_yzz_0[j] + 0.5 * fl1_fx * tly_zzz_zz_0[j];

            tlz_yzzz_yzz_0[j] =
                pa_y[j] * tlz_zzz_yzz_0[j] + 0.5 * fl1_fx * tlz_zzz_zz_0[j] + 0.5 * fl1_fx * tpx_zzz_yzz_0[j] + fl1_fx * fl1_fgb * tdx_zzz_yzz_0[j];

            tlx_yzzz_zzz_0[j] = pa_y[j] * tlx_zzz_zzz_0[j] - 0.5 * fl1_fx * tpz_zzz_zzz_0[j] - fl1_fx * fl1_fgb * tdz_zzz_zzz_0[j];

            tly_yzzz_zzz_0[j] = pa_y[j] * tly_zzz_zzz_0[j];

            tlz_yzzz_zzz_0[j] = pa_y[j] * tlz_zzz_zzz_0[j] + 0.5 * fl1_fx * tpx_zzz_zzz_0[j] + fl1_fx * fl1_fgb * tdx_zzz_zzz_0[j];

            tlx_zzzz_xxx_0[j] =
                pa_z[j] * tlx_zzz_xxx_0[j] + 1.5 * fl1_fx * tlx_zz_xxx_0[j] + 0.5 * fl1_fx * tpy_zzz_xxx_0[j] + fl1_fx * fl1_fgb * tdy_zzz_xxx_0[j];

            tly_zzzz_xxx_0[j] =
                pa_z[j] * tly_zzz_xxx_0[j] + 1.5 * fl1_fx * tly_zz_xxx_0[j] - 0.5 * fl1_fx * tpx_zzz_xxx_0[j] - fl1_fx * fl1_fgb * tdx_zzz_xxx_0[j];

            tlz_zzzz_xxx_0[j] = pa_z[j] * tlz_zzz_xxx_0[j] + 1.5 * fl1_fx * tlz_zz_xxx_0[j];

            tlx_zzzz_xxy_0[j] =
                pa_z[j] * tlx_zzz_xxy_0[j] + 1.5 * fl1_fx * tlx_zz_xxy_0[j] + 0.5 * fl1_fx * tpy_zzz_xxy_0[j] + fl1_fx * fl1_fgb * tdy_zzz_xxy_0[j];

            tly_zzzz_xxy_0[j] =
                pa_z[j] * tly_zzz_xxy_0[j] + 1.5 * fl1_fx * tly_zz_xxy_0[j] - 0.5 * fl1_fx * tpx_zzz_xxy_0[j] - fl1_fx * fl1_fgb * tdx_zzz_xxy_0[j];

            tlz_zzzz_xxy_0[j] = pa_z[j] * tlz_zzz_xxy_0[j] + 1.5 * fl1_fx * tlz_zz_xxy_0[j];

            tlx_zzzz_xxz_0[j] = pa_z[j] * tlx_zzz_xxz_0[j] + 1.5 * fl1_fx * tlx_zz_xxz_0[j] + 0.5 * fl1_fx * tlx_zzz_xx_0[j] +
                                0.5 * fl1_fx * tpy_zzz_xxz_0[j] + fl1_fx * fl1_fgb * tdy_zzz_xxz_0[j];

            tly_zzzz_xxz_0[j] = pa_z[j] * tly_zzz_xxz_0[j] + 1.5 * fl1_fx * tly_zz_xxz_0[j] + 0.5 * fl1_fx * tly_zzz_xx_0[j] -
                                0.5 * fl1_fx * tpx_zzz_xxz_0[j] - fl1_fx * fl1_fgb * tdx_zzz_xxz_0[j];

            tlz_zzzz_xxz_0[j] = pa_z[j] * tlz_zzz_xxz_0[j] + 1.5 * fl1_fx * tlz_zz_xxz_0[j] + 0.5 * fl1_fx * tlz_zzz_xx_0[j];

            tlx_zzzz_xyy_0[j] =
                pa_z[j] * tlx_zzz_xyy_0[j] + 1.5 * fl1_fx * tlx_zz_xyy_0[j] + 0.5 * fl1_fx * tpy_zzz_xyy_0[j] + fl1_fx * fl1_fgb * tdy_zzz_xyy_0[j];

            tly_zzzz_xyy_0[j] =
                pa_z[j] * tly_zzz_xyy_0[j] + 1.5 * fl1_fx * tly_zz_xyy_0[j] - 0.5 * fl1_fx * tpx_zzz_xyy_0[j] - fl1_fx * fl1_fgb * tdx_zzz_xyy_0[j];

            tlz_zzzz_xyy_0[j] = pa_z[j] * tlz_zzz_xyy_0[j] + 1.5 * fl1_fx * tlz_zz_xyy_0[j];

            tlx_zzzz_xyz_0[j] = pa_z[j] * tlx_zzz_xyz_0[j] + 1.5 * fl1_fx * tlx_zz_xyz_0[j] + 0.5 * fl1_fx * tlx_zzz_xy_0[j] +
                                0.5 * fl1_fx * tpy_zzz_xyz_0[j] + fl1_fx * fl1_fgb * tdy_zzz_xyz_0[j];

            tly_zzzz_xyz_0[j] = pa_z[j] * tly_zzz_xyz_0[j] + 1.5 * fl1_fx * tly_zz_xyz_0[j] + 0.5 * fl1_fx * tly_zzz_xy_0[j] -
                                0.5 * fl1_fx * tpx_zzz_xyz_0[j] - fl1_fx * fl1_fgb * tdx_zzz_xyz_0[j];

            tlz_zzzz_xyz_0[j] = pa_z[j] * tlz_zzz_xyz_0[j] + 1.5 * fl1_fx * tlz_zz_xyz_0[j] + 0.5 * fl1_fx * tlz_zzz_xy_0[j];

            tlx_zzzz_xzz_0[j] = pa_z[j] * tlx_zzz_xzz_0[j] + 1.5 * fl1_fx * tlx_zz_xzz_0[j] + fl1_fx * tlx_zzz_xz_0[j] +
                                0.5 * fl1_fx * tpy_zzz_xzz_0[j] + fl1_fx * fl1_fgb * tdy_zzz_xzz_0[j];

            tly_zzzz_xzz_0[j] = pa_z[j] * tly_zzz_xzz_0[j] + 1.5 * fl1_fx * tly_zz_xzz_0[j] + fl1_fx * tly_zzz_xz_0[j] -
                                0.5 * fl1_fx * tpx_zzz_xzz_0[j] - fl1_fx * fl1_fgb * tdx_zzz_xzz_0[j];

            tlz_zzzz_xzz_0[j] = pa_z[j] * tlz_zzz_xzz_0[j] + 1.5 * fl1_fx * tlz_zz_xzz_0[j] + fl1_fx * tlz_zzz_xz_0[j];

            tlx_zzzz_yyy_0[j] =
                pa_z[j] * tlx_zzz_yyy_0[j] + 1.5 * fl1_fx * tlx_zz_yyy_0[j] + 0.5 * fl1_fx * tpy_zzz_yyy_0[j] + fl1_fx * fl1_fgb * tdy_zzz_yyy_0[j];

            tly_zzzz_yyy_0[j] =
                pa_z[j] * tly_zzz_yyy_0[j] + 1.5 * fl1_fx * tly_zz_yyy_0[j] - 0.5 * fl1_fx * tpx_zzz_yyy_0[j] - fl1_fx * fl1_fgb * tdx_zzz_yyy_0[j];

            tlz_zzzz_yyy_0[j] = pa_z[j] * tlz_zzz_yyy_0[j] + 1.5 * fl1_fx * tlz_zz_yyy_0[j];

            tlx_zzzz_yyz_0[j] = pa_z[j] * tlx_zzz_yyz_0[j] + 1.5 * fl1_fx * tlx_zz_yyz_0[j] + 0.5 * fl1_fx * tlx_zzz_yy_0[j] +
                                0.5 * fl1_fx * tpy_zzz_yyz_0[j] + fl1_fx * fl1_fgb * tdy_zzz_yyz_0[j];

            tly_zzzz_yyz_0[j] = pa_z[j] * tly_zzz_yyz_0[j] + 1.5 * fl1_fx * tly_zz_yyz_0[j] + 0.5 * fl1_fx * tly_zzz_yy_0[j] -
                                0.5 * fl1_fx * tpx_zzz_yyz_0[j] - fl1_fx * fl1_fgb * tdx_zzz_yyz_0[j];

            tlz_zzzz_yyz_0[j] = pa_z[j] * tlz_zzz_yyz_0[j] + 1.5 * fl1_fx * tlz_zz_yyz_0[j] + 0.5 * fl1_fx * tlz_zzz_yy_0[j];

            tlx_zzzz_yzz_0[j] = pa_z[j] * tlx_zzz_yzz_0[j] + 1.5 * fl1_fx * tlx_zz_yzz_0[j] + fl1_fx * tlx_zzz_yz_0[j] +
                                0.5 * fl1_fx * tpy_zzz_yzz_0[j] + fl1_fx * fl1_fgb * tdy_zzz_yzz_0[j];

            tly_zzzz_yzz_0[j] = pa_z[j] * tly_zzz_yzz_0[j] + 1.5 * fl1_fx * tly_zz_yzz_0[j] + fl1_fx * tly_zzz_yz_0[j] -
                                0.5 * fl1_fx * tpx_zzz_yzz_0[j] - fl1_fx * fl1_fgb * tdx_zzz_yzz_0[j];

            tlz_zzzz_yzz_0[j] = pa_z[j] * tlz_zzz_yzz_0[j] + 1.5 * fl1_fx * tlz_zz_yzz_0[j] + fl1_fx * tlz_zzz_yz_0[j];

            tlx_zzzz_zzz_0[j] = pa_z[j] * tlx_zzz_zzz_0[j] + 1.5 * fl1_fx * tlx_zz_zzz_0[j] + 1.5 * fl1_fx * tlx_zzz_zz_0[j] +
                                0.5 * fl1_fx * tpy_zzz_zzz_0[j] + fl1_fx * fl1_fgb * tdy_zzz_zzz_0[j];

            tly_zzzz_zzz_0[j] = pa_z[j] * tly_zzz_zzz_0[j] + 1.5 * fl1_fx * tly_zz_zzz_0[j] + 1.5 * fl1_fx * tly_zzz_zz_0[j] -
                                0.5 * fl1_fx * tpx_zzz_zzz_0[j] - fl1_fx * fl1_fgb * tdx_zzz_zzz_0[j];

            tlz_zzzz_zzz_0[j] = pa_z[j] * tlz_zzz_zzz_0[j] + 1.5 * fl1_fx * tlz_zz_zzz_0[j] + 1.5 * fl1_fx * tlz_zzz_zz_0[j];
        }

        idx++;
    }
}

}  // namespace amomrecfunc
