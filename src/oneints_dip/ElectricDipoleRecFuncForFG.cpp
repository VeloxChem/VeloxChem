//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ElectricDipoleRecFuncForFG.hpp"

namespace ediprecfunc {  // ediprecfunc namespace

void
compElectricDipoleForFG(CMemBlock2D<double>&       primBuffer,
                        const CRecursionMap&       recursionMap,
                        const CMemBlock2D<double>& osFactors,
                        const int32_t              nOSFactors,
                        const CMemBlock2D<double>& paDistances,
                        const CGtoBlock&           braGtoBlock,
                        const CGtoBlock&           ketGtoBlock,
                        const int32_t              iContrGto)
{
    ediprecfunc::compElectricDipoleForFG_0_50(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForFG_50_100(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForFG_100_150(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForFG_150_200(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForFG_200_250(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForFG_250_300(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForFG_300_350(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForFG_350_400(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForFG_400_450(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compElectricDipoleForFG_0_50(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const int32_t              nOSFactors,
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

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_d_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

        auto tdx_xx_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx);

        auto tdy_xx_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx);

        auto tdz_xx_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx);

        auto tdx_xx_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 1);

        auto tdy_xx_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 1);

        auto tdz_xx_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 1);

        auto tdx_xx_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 2);

        auto tdy_xx_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 2);

        auto tdz_xx_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 2);

        auto tdx_xx_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 3);

        auto tdy_xx_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 3);

        auto tdz_xx_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 3);

        auto tdx_xx_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 4);

        auto tdy_xx_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 4);

        auto tdz_xx_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 4);

        auto tdx_xx_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 5);

        auto tdy_xx_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 5);

        auto tdz_xx_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 5);

        auto tdx_xx_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 6);

        auto tdy_xx_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 6);

        auto tdz_xx_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 6);

        auto tdx_xx_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 7);

        auto tdy_xx_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 7);

        auto tdz_xx_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 7);

        auto tdx_xx_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 8);

        auto tdy_xx_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 8);

        auto tdz_xx_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 8);

        auto tdx_xx_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 9);

        auto tdy_xx_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 9);

        auto tdz_xx_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 9);

        auto tdx_xx_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 10);

        auto tdy_xx_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 10);

        auto tdz_xx_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 10);

        auto tdx_xx_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 11);

        auto tdy_xx_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 11);

        auto tdz_xx_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 11);

        auto tdx_xx_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 12);

        auto tdy_xx_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 12);

        auto tdz_xx_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 12);

        auto tdx_xx_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 13);

        auto tdy_xx_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 13);

        auto tdz_xx_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 13);

        auto tdx_xx_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 14);

        auto tdy_xx_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 14);

        auto tdz_xx_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 14);

        auto tdx_xy_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 15);

        auto tdy_xy_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 15);

        auto tdz_xy_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 15);

        auto tdx_xy_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 16);

        auto tdy_xy_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 16);

        auto tdx_x_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx);

        auto tdy_x_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx);

        auto tdz_x_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx);

        auto tdx_x_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 1);

        auto tdy_x_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 1);

        auto tdz_x_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 1);

        auto tdx_x_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 2);

        auto tdy_x_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 2);

        auto tdz_x_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 2);

        auto tdx_x_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 3);

        auto tdy_x_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 3);

        auto tdz_x_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 3);

        auto tdx_x_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 4);

        auto tdy_x_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 4);

        auto tdz_x_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 4);

        auto tdx_x_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 5);

        auto tdy_x_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 5);

        auto tdz_x_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 5);

        auto tdx_x_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 6);

        auto tdy_x_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 6);

        auto tdz_x_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 6);

        auto tdx_x_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 7);

        auto tdy_x_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 7);

        auto tdz_x_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 7);

        auto tdx_x_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 8);

        auto tdy_x_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 8);

        auto tdz_x_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 8);

        auto tdx_x_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 9);

        auto tdy_x_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 9);

        auto tdz_x_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 9);

        auto tdx_x_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 10);

        auto tdy_x_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 10);

        auto tdz_x_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 10);

        auto tdx_x_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 11);

        auto tdy_x_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 11);

        auto tdz_x_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 11);

        auto tdx_x_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 12);

        auto tdy_x_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 12);

        auto tdz_x_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 12);

        auto tdx_x_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 13);

        auto tdy_x_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 13);

        auto tdz_x_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 13);

        auto tdx_x_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 14);

        auto tdy_x_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 14);

        auto tdz_x_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 14);

        auto tdx_y_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 15);

        auto tdy_y_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 15);

        auto tdz_y_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 15);

        auto tdx_y_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 16);

        auto tdy_y_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 16);

        auto tdx_xx_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx);

        auto tdy_xx_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx);

        auto tdz_xx_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx);

        auto tdx_xx_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 1);

        auto tdy_xx_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 1);

        auto tdz_xx_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 1);

        auto tdx_xx_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 2);

        auto tdy_xx_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 2);

        auto tdz_xx_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 2);

        auto tdx_xx_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 3);

        auto tdy_xx_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 3);

        auto tdz_xx_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 3);

        auto tdx_xx_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 4);

        auto tdy_xx_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 4);

        auto tdz_xx_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 4);

        auto tdx_xx_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 5);

        auto tdy_xx_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 5);

        auto tdz_xx_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 5);

        auto tdx_xx_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 6);

        auto tdy_xx_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 6);

        auto tdz_xx_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 6);

        auto tdx_xx_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 7);

        auto tdy_xx_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 7);

        auto tdz_xx_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 7);

        auto tdx_xx_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 8);

        auto tdy_xx_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 8);

        auto tdz_xx_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 8);

        auto tdx_xx_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 9);

        auto tdy_xx_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 9);

        auto tdz_xx_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 9);

        auto tdx_xy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 10);

        auto tdy_xy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 10);

        auto tdz_xy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 10);

        auto tdx_xy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 11);

        auto tdy_xy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 11);

        auto ts_xx_xxxx_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx);

        auto ts_xx_xxxy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 1);

        auto ts_xx_xxxz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 2);

        auto ts_xx_xxyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 3);

        auto ts_xx_xxyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 4);

        auto ts_xx_xxzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 5);

        auto ts_xx_xyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 6);

        auto ts_xx_xyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 7);

        auto ts_xx_xyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 8);

        auto ts_xx_xzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 9);

        auto ts_xx_yyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 10);

        auto ts_xx_yyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 11);

        auto ts_xx_yyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 12);

        auto ts_xx_yzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 13);

        auto ts_xx_zzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 14);

        auto ts_xy_xxxx_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 15);

        auto ts_xy_xxxy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 16);

        // set up pointers to integrals

        auto tdx_xxx_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx);

        auto tdy_xxx_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx);

        auto tdz_xxx_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx);

        auto tdx_xxx_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 1);

        auto tdy_xxx_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 1);

        auto tdz_xxx_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 1);

        auto tdx_xxx_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 2);

        auto tdy_xxx_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 2);

        auto tdz_xxx_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 2);

        auto tdx_xxx_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 3);

        auto tdy_xxx_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 3);

        auto tdz_xxx_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 3);

        auto tdx_xxx_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 4);

        auto tdy_xxx_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 4);

        auto tdz_xxx_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 4);

        auto tdx_xxx_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 5);

        auto tdy_xxx_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 5);

        auto tdz_xxx_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 5);

        auto tdx_xxx_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 6);

        auto tdy_xxx_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 6);

        auto tdz_xxx_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 6);

        auto tdx_xxx_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 7);

        auto tdy_xxx_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 7);

        auto tdz_xxx_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 7);

        auto tdx_xxx_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 8);

        auto tdy_xxx_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 8);

        auto tdz_xxx_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 8);

        auto tdx_xxx_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 9);

        auto tdy_xxx_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 9);

        auto tdz_xxx_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 9);

        auto tdx_xxx_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 10);

        auto tdy_xxx_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 10);

        auto tdz_xxx_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 10);

        auto tdx_xxx_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 11);

        auto tdy_xxx_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 11);

        auto tdz_xxx_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 11);

        auto tdx_xxx_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 12);

        auto tdy_xxx_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 12);

        auto tdz_xxx_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 12);

        auto tdx_xxx_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 13);

        auto tdy_xxx_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 13);

        auto tdz_xxx_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 13);

        auto tdx_xxx_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 14);

        auto tdy_xxx_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 14);

        auto tdz_xxx_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 14);

        auto tdx_xxy_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 15);

        auto tdy_xxy_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 15);

        auto tdz_xxy_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 15);

        auto tdx_xxy_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 16);

        auto tdy_xxy_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 16);

        // Batch of Integrals (0,50)

        #pragma omp simd aligned(fx, pa_x, tdx_x_xxxx_0, tdx_x_xxxy_0, tdx_x_xxxz_0, tdx_x_xxyy_0, \
                                     tdx_x_xxyz_0, tdx_x_xxzz_0, tdx_x_xyyy_0, tdx_x_xyyz_0, tdx_x_xyzz_0, tdx_x_xzzz_0, \
                                     tdx_x_yyyy_0, tdx_x_yyyz_0, tdx_x_yyzz_0, tdx_x_yzzz_0, tdx_x_zzzz_0, tdx_xx_xxx_0, \
                                     tdx_xx_xxxx_0, tdx_xx_xxxy_0, tdx_xx_xxxz_0, tdx_xx_xxy_0, tdx_xx_xxyy_0, \
                                     tdx_xx_xxyz_0, tdx_xx_xxz_0, tdx_xx_xxzz_0, tdx_xx_xyy_0, tdx_xx_xyyy_0, \
                                     tdx_xx_xyyz_0, tdx_xx_xyz_0, tdx_xx_xyzz_0, tdx_xx_xzz_0, tdx_xx_xzzz_0, \
                                     tdx_xx_yyy_0, tdx_xx_yyyy_0, tdx_xx_yyyz_0, tdx_xx_yyz_0, tdx_xx_yyzz_0, \
                                     tdx_xx_yzz_0, tdx_xx_yzzz_0, tdx_xx_zzz_0, tdx_xx_zzzz_0, tdx_xxx_xxxx_0, \
                                     tdx_xxx_xxxy_0, tdx_xxx_xxxz_0, tdx_xxx_xxyy_0, tdx_xxx_xxyz_0, tdx_xxx_xxzz_0, \
                                     tdx_xxx_xyyy_0, tdx_xxx_xyyz_0, tdx_xxx_xyzz_0, tdx_xxx_xzzz_0, tdx_xxx_yyyy_0, \
                                     tdx_xxx_yyyz_0, tdx_xxx_yyzz_0, tdx_xxx_yzzz_0, tdx_xxx_zzzz_0, tdx_xxy_xxxx_0, \
                                     tdx_xxy_xxxy_0, tdx_xy_xxx_0, tdx_xy_xxxx_0, tdx_xy_xxxy_0, tdx_xy_xxy_0, \
                                     tdx_y_xxxx_0, tdx_y_xxxy_0, tdy_x_xxxx_0, tdy_x_xxxy_0, tdy_x_xxxz_0, tdy_x_xxyy_0, \
                                     tdy_x_xxyz_0, tdy_x_xxzz_0, tdy_x_xyyy_0, tdy_x_xyyz_0, tdy_x_xyzz_0, tdy_x_xzzz_0, \
                                     tdy_x_yyyy_0, tdy_x_yyyz_0, tdy_x_yyzz_0, tdy_x_yzzz_0, tdy_x_zzzz_0, tdy_xx_xxx_0, \
                                     tdy_xx_xxxx_0, tdy_xx_xxxy_0, tdy_xx_xxxz_0, tdy_xx_xxy_0, tdy_xx_xxyy_0, \
                                     tdy_xx_xxyz_0, tdy_xx_xxz_0, tdy_xx_xxzz_0, tdy_xx_xyy_0, tdy_xx_xyyy_0, \
                                     tdy_xx_xyyz_0, tdy_xx_xyz_0, tdy_xx_xyzz_0, tdy_xx_xzz_0, tdy_xx_xzzz_0, \
                                     tdy_xx_yyy_0, tdy_xx_yyyy_0, tdy_xx_yyyz_0, tdy_xx_yyz_0, tdy_xx_yyzz_0, \
                                     tdy_xx_yzz_0, tdy_xx_yzzz_0, tdy_xx_zzz_0, tdy_xx_zzzz_0, tdy_xxx_xxxx_0, \
                                     tdy_xxx_xxxy_0, tdy_xxx_xxxz_0, tdy_xxx_xxyy_0, tdy_xxx_xxyz_0, tdy_xxx_xxzz_0, \
                                     tdy_xxx_xyyy_0, tdy_xxx_xyyz_0, tdy_xxx_xyzz_0, tdy_xxx_xzzz_0, tdy_xxx_yyyy_0, \
                                     tdy_xxx_yyyz_0, tdy_xxx_yyzz_0, tdy_xxx_yzzz_0, tdy_xxx_zzzz_0, tdy_xxy_xxxx_0, \
                                     tdy_xxy_xxxy_0, tdy_xy_xxx_0, tdy_xy_xxxx_0, tdy_xy_xxxy_0, tdy_xy_xxy_0, \
                                     tdy_y_xxxx_0, tdy_y_xxxy_0, tdz_x_xxxx_0, tdz_x_xxxy_0, tdz_x_xxxz_0, tdz_x_xxyy_0, \
                                     tdz_x_xxyz_0, tdz_x_xxzz_0, tdz_x_xyyy_0, tdz_x_xyyz_0, tdz_x_xyzz_0, tdz_x_xzzz_0, \
                                     tdz_x_yyyy_0, tdz_x_yyyz_0, tdz_x_yyzz_0, tdz_x_yzzz_0, tdz_x_zzzz_0, tdz_xx_xxx_0, \
                                     tdz_xx_xxxx_0, tdz_xx_xxxy_0, tdz_xx_xxxz_0, tdz_xx_xxy_0, tdz_xx_xxyy_0, \
                                     tdz_xx_xxyz_0, tdz_xx_xxz_0, tdz_xx_xxzz_0, tdz_xx_xyy_0, tdz_xx_xyyy_0, \
                                     tdz_xx_xyyz_0, tdz_xx_xyz_0, tdz_xx_xyzz_0, tdz_xx_xzz_0, tdz_xx_xzzz_0, \
                                     tdz_xx_yyy_0, tdz_xx_yyyy_0, tdz_xx_yyyz_0, tdz_xx_yyz_0, tdz_xx_yyzz_0, \
                                     tdz_xx_yzz_0, tdz_xx_yzzz_0, tdz_xx_zzz_0, tdz_xx_zzzz_0, tdz_xxx_xxxx_0, \
                                     tdz_xxx_xxxy_0, tdz_xxx_xxxz_0, tdz_xxx_xxyy_0, tdz_xxx_xxyz_0, tdz_xxx_xxzz_0, \
                                     tdz_xxx_xyyy_0, tdz_xxx_xyyz_0, tdz_xxx_xyzz_0, tdz_xxx_xzzz_0, tdz_xxx_yyyy_0, \
                                     tdz_xxx_yyyz_0, tdz_xxx_yyzz_0, tdz_xxx_yzzz_0, tdz_xxx_zzzz_0, tdz_xxy_xxxx_0, \
                                     tdz_xy_xxx_0, tdz_xy_xxxx_0, tdz_y_xxxx_0, ts_xx_xxxx_0, ts_xx_xxxy_0, \
                                     ts_xx_xxxz_0, ts_xx_xxyy_0, ts_xx_xxyz_0, ts_xx_xxzz_0, ts_xx_xyyy_0, ts_xx_xyyz_0, \
                                     ts_xx_xyzz_0, ts_xx_xzzz_0, ts_xx_yyyy_0, ts_xx_yyyz_0, ts_xx_yyzz_0, ts_xx_yzzz_0, \
                                     ts_xx_zzzz_0, ts_xy_xxxx_0, ts_xy_xxxy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            tdx_xxx_xxxx_0[j] =
                pa_x[j] * tdx_xx_xxxx_0[j] + fl1_fx * tdx_x_xxxx_0[j] + 2.0 * fl1_fx * tdx_xx_xxx_0[j] + 0.5 * fl1_fx * ts_xx_xxxx_0[j];

            tdy_xxx_xxxx_0[j] = pa_x[j] * tdy_xx_xxxx_0[j] + fl1_fx * tdy_x_xxxx_0[j] + 2.0 * fl1_fx * tdy_xx_xxx_0[j];

            tdz_xxx_xxxx_0[j] = pa_x[j] * tdz_xx_xxxx_0[j] + fl1_fx * tdz_x_xxxx_0[j] + 2.0 * fl1_fx * tdz_xx_xxx_0[j];

            tdx_xxx_xxxy_0[j] =
                pa_x[j] * tdx_xx_xxxy_0[j] + fl1_fx * tdx_x_xxxy_0[j] + 1.5 * fl1_fx * tdx_xx_xxy_0[j] + 0.5 * fl1_fx * ts_xx_xxxy_0[j];

            tdy_xxx_xxxy_0[j] = pa_x[j] * tdy_xx_xxxy_0[j] + fl1_fx * tdy_x_xxxy_0[j] + 1.5 * fl1_fx * tdy_xx_xxy_0[j];

            tdz_xxx_xxxy_0[j] = pa_x[j] * tdz_xx_xxxy_0[j] + fl1_fx * tdz_x_xxxy_0[j] + 1.5 * fl1_fx * tdz_xx_xxy_0[j];

            tdx_xxx_xxxz_0[j] =
                pa_x[j] * tdx_xx_xxxz_0[j] + fl1_fx * tdx_x_xxxz_0[j] + 1.5 * fl1_fx * tdx_xx_xxz_0[j] + 0.5 * fl1_fx * ts_xx_xxxz_0[j];

            tdy_xxx_xxxz_0[j] = pa_x[j] * tdy_xx_xxxz_0[j] + fl1_fx * tdy_x_xxxz_0[j] + 1.5 * fl1_fx * tdy_xx_xxz_0[j];

            tdz_xxx_xxxz_0[j] = pa_x[j] * tdz_xx_xxxz_0[j] + fl1_fx * tdz_x_xxxz_0[j] + 1.5 * fl1_fx * tdz_xx_xxz_0[j];

            tdx_xxx_xxyy_0[j] = pa_x[j] * tdx_xx_xxyy_0[j] + fl1_fx * tdx_x_xxyy_0[j] + fl1_fx * tdx_xx_xyy_0[j] + 0.5 * fl1_fx * ts_xx_xxyy_0[j];

            tdy_xxx_xxyy_0[j] = pa_x[j] * tdy_xx_xxyy_0[j] + fl1_fx * tdy_x_xxyy_0[j] + fl1_fx * tdy_xx_xyy_0[j];

            tdz_xxx_xxyy_0[j] = pa_x[j] * tdz_xx_xxyy_0[j] + fl1_fx * tdz_x_xxyy_0[j] + fl1_fx * tdz_xx_xyy_0[j];

            tdx_xxx_xxyz_0[j] = pa_x[j] * tdx_xx_xxyz_0[j] + fl1_fx * tdx_x_xxyz_0[j] + fl1_fx * tdx_xx_xyz_0[j] + 0.5 * fl1_fx * ts_xx_xxyz_0[j];

            tdy_xxx_xxyz_0[j] = pa_x[j] * tdy_xx_xxyz_0[j] + fl1_fx * tdy_x_xxyz_0[j] + fl1_fx * tdy_xx_xyz_0[j];

            tdz_xxx_xxyz_0[j] = pa_x[j] * tdz_xx_xxyz_0[j] + fl1_fx * tdz_x_xxyz_0[j] + fl1_fx * tdz_xx_xyz_0[j];

            tdx_xxx_xxzz_0[j] = pa_x[j] * tdx_xx_xxzz_0[j] + fl1_fx * tdx_x_xxzz_0[j] + fl1_fx * tdx_xx_xzz_0[j] + 0.5 * fl1_fx * ts_xx_xxzz_0[j];

            tdy_xxx_xxzz_0[j] = pa_x[j] * tdy_xx_xxzz_0[j] + fl1_fx * tdy_x_xxzz_0[j] + fl1_fx * tdy_xx_xzz_0[j];

            tdz_xxx_xxzz_0[j] = pa_x[j] * tdz_xx_xxzz_0[j] + fl1_fx * tdz_x_xxzz_0[j] + fl1_fx * tdz_xx_xzz_0[j];

            tdx_xxx_xyyy_0[j] =
                pa_x[j] * tdx_xx_xyyy_0[j] + fl1_fx * tdx_x_xyyy_0[j] + 0.5 * fl1_fx * tdx_xx_yyy_0[j] + 0.5 * fl1_fx * ts_xx_xyyy_0[j];

            tdy_xxx_xyyy_0[j] = pa_x[j] * tdy_xx_xyyy_0[j] + fl1_fx * tdy_x_xyyy_0[j] + 0.5 * fl1_fx * tdy_xx_yyy_0[j];

            tdz_xxx_xyyy_0[j] = pa_x[j] * tdz_xx_xyyy_0[j] + fl1_fx * tdz_x_xyyy_0[j] + 0.5 * fl1_fx * tdz_xx_yyy_0[j];

            tdx_xxx_xyyz_0[j] =
                pa_x[j] * tdx_xx_xyyz_0[j] + fl1_fx * tdx_x_xyyz_0[j] + 0.5 * fl1_fx * tdx_xx_yyz_0[j] + 0.5 * fl1_fx * ts_xx_xyyz_0[j];

            tdy_xxx_xyyz_0[j] = pa_x[j] * tdy_xx_xyyz_0[j] + fl1_fx * tdy_x_xyyz_0[j] + 0.5 * fl1_fx * tdy_xx_yyz_0[j];

            tdz_xxx_xyyz_0[j] = pa_x[j] * tdz_xx_xyyz_0[j] + fl1_fx * tdz_x_xyyz_0[j] + 0.5 * fl1_fx * tdz_xx_yyz_0[j];

            tdx_xxx_xyzz_0[j] =
                pa_x[j] * tdx_xx_xyzz_0[j] + fl1_fx * tdx_x_xyzz_0[j] + 0.5 * fl1_fx * tdx_xx_yzz_0[j] + 0.5 * fl1_fx * ts_xx_xyzz_0[j];

            tdy_xxx_xyzz_0[j] = pa_x[j] * tdy_xx_xyzz_0[j] + fl1_fx * tdy_x_xyzz_0[j] + 0.5 * fl1_fx * tdy_xx_yzz_0[j];

            tdz_xxx_xyzz_0[j] = pa_x[j] * tdz_xx_xyzz_0[j] + fl1_fx * tdz_x_xyzz_0[j] + 0.5 * fl1_fx * tdz_xx_yzz_0[j];

            tdx_xxx_xzzz_0[j] =
                pa_x[j] * tdx_xx_xzzz_0[j] + fl1_fx * tdx_x_xzzz_0[j] + 0.5 * fl1_fx * tdx_xx_zzz_0[j] + 0.5 * fl1_fx * ts_xx_xzzz_0[j];

            tdy_xxx_xzzz_0[j] = pa_x[j] * tdy_xx_xzzz_0[j] + fl1_fx * tdy_x_xzzz_0[j] + 0.5 * fl1_fx * tdy_xx_zzz_0[j];

            tdz_xxx_xzzz_0[j] = pa_x[j] * tdz_xx_xzzz_0[j] + fl1_fx * tdz_x_xzzz_0[j] + 0.5 * fl1_fx * tdz_xx_zzz_0[j];

            tdx_xxx_yyyy_0[j] = pa_x[j] * tdx_xx_yyyy_0[j] + fl1_fx * tdx_x_yyyy_0[j] + 0.5 * fl1_fx * ts_xx_yyyy_0[j];

            tdy_xxx_yyyy_0[j] = pa_x[j] * tdy_xx_yyyy_0[j] + fl1_fx * tdy_x_yyyy_0[j];

            tdz_xxx_yyyy_0[j] = pa_x[j] * tdz_xx_yyyy_0[j] + fl1_fx * tdz_x_yyyy_0[j];

            tdx_xxx_yyyz_0[j] = pa_x[j] * tdx_xx_yyyz_0[j] + fl1_fx * tdx_x_yyyz_0[j] + 0.5 * fl1_fx * ts_xx_yyyz_0[j];

            tdy_xxx_yyyz_0[j] = pa_x[j] * tdy_xx_yyyz_0[j] + fl1_fx * tdy_x_yyyz_0[j];

            tdz_xxx_yyyz_0[j] = pa_x[j] * tdz_xx_yyyz_0[j] + fl1_fx * tdz_x_yyyz_0[j];

            tdx_xxx_yyzz_0[j] = pa_x[j] * tdx_xx_yyzz_0[j] + fl1_fx * tdx_x_yyzz_0[j] + 0.5 * fl1_fx * ts_xx_yyzz_0[j];

            tdy_xxx_yyzz_0[j] = pa_x[j] * tdy_xx_yyzz_0[j] + fl1_fx * tdy_x_yyzz_0[j];

            tdz_xxx_yyzz_0[j] = pa_x[j] * tdz_xx_yyzz_0[j] + fl1_fx * tdz_x_yyzz_0[j];

            tdx_xxx_yzzz_0[j] = pa_x[j] * tdx_xx_yzzz_0[j] + fl1_fx * tdx_x_yzzz_0[j] + 0.5 * fl1_fx * ts_xx_yzzz_0[j];

            tdy_xxx_yzzz_0[j] = pa_x[j] * tdy_xx_yzzz_0[j] + fl1_fx * tdy_x_yzzz_0[j];

            tdz_xxx_yzzz_0[j] = pa_x[j] * tdz_xx_yzzz_0[j] + fl1_fx * tdz_x_yzzz_0[j];

            tdx_xxx_zzzz_0[j] = pa_x[j] * tdx_xx_zzzz_0[j] + fl1_fx * tdx_x_zzzz_0[j] + 0.5 * fl1_fx * ts_xx_zzzz_0[j];

            tdy_xxx_zzzz_0[j] = pa_x[j] * tdy_xx_zzzz_0[j] + fl1_fx * tdy_x_zzzz_0[j];

            tdz_xxx_zzzz_0[j] = pa_x[j] * tdz_xx_zzzz_0[j] + fl1_fx * tdz_x_zzzz_0[j];

            tdx_xxy_xxxx_0[j] =
                pa_x[j] * tdx_xy_xxxx_0[j] + 0.5 * fl1_fx * tdx_y_xxxx_0[j] + 2.0 * fl1_fx * tdx_xy_xxx_0[j] + 0.5 * fl1_fx * ts_xy_xxxx_0[j];

            tdy_xxy_xxxx_0[j] = pa_x[j] * tdy_xy_xxxx_0[j] + 0.5 * fl1_fx * tdy_y_xxxx_0[j] + 2.0 * fl1_fx * tdy_xy_xxx_0[j];

            tdz_xxy_xxxx_0[j] = pa_x[j] * tdz_xy_xxxx_0[j] + 0.5 * fl1_fx * tdz_y_xxxx_0[j] + 2.0 * fl1_fx * tdz_xy_xxx_0[j];

            tdx_xxy_xxxy_0[j] =
                pa_x[j] * tdx_xy_xxxy_0[j] + 0.5 * fl1_fx * tdx_y_xxxy_0[j] + 1.5 * fl1_fx * tdx_xy_xxy_0[j] + 0.5 * fl1_fx * ts_xy_xxxy_0[j];

            tdy_xxy_xxxy_0[j] = pa_x[j] * tdy_xy_xxxy_0[j] + 0.5 * fl1_fx * tdy_y_xxxy_0[j] + 1.5 * fl1_fx * tdy_xy_xxy_0[j];
        }

        idx++;
    }
}

void
compElectricDipoleForFG_50_100(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const int32_t              nOSFactors,
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

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_d_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

        auto tdz_xy_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 16);

        auto tdx_xy_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 17);

        auto tdy_xy_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 17);

        auto tdz_xy_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 17);

        auto tdx_xy_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 18);

        auto tdy_xy_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 18);

        auto tdz_xy_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 18);

        auto tdx_xy_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 19);

        auto tdy_xy_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 19);

        auto tdz_xy_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 19);

        auto tdx_xy_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 20);

        auto tdy_xy_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 20);

        auto tdz_xy_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 20);

        auto tdx_xy_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 21);

        auto tdy_xy_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 21);

        auto tdz_xy_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 21);

        auto tdx_xy_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 22);

        auto tdy_xy_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 22);

        auto tdz_xy_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 22);

        auto tdx_xy_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 23);

        auto tdy_xy_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 23);

        auto tdz_xy_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 23);

        auto tdx_xy_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 24);

        auto tdy_xy_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 24);

        auto tdz_xy_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 24);

        auto tdx_xy_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 25);

        auto tdy_xy_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 25);

        auto tdz_xy_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 25);

        auto tdx_xy_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 26);

        auto tdy_xy_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 26);

        auto tdz_xy_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 26);

        auto tdx_xy_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 27);

        auto tdy_xy_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 27);

        auto tdz_xy_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 27);

        auto tdx_xy_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 28);

        auto tdy_xy_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 28);

        auto tdz_xy_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 28);

        auto tdx_xy_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 29);

        auto tdy_xy_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 29);

        auto tdz_xy_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 29);

        auto tdx_xz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 30);

        auto tdy_xz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 30);

        auto tdz_xz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 30);

        auto tdx_xz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 31);

        auto tdy_xz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 31);

        auto tdz_xz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 31);

        auto tdx_xz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 32);

        auto tdy_xz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 32);

        auto tdz_xz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 32);

        auto tdx_xz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 33);

        auto tdz_y_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 16);

        auto tdx_y_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 17);

        auto tdy_y_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 17);

        auto tdz_y_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 17);

        auto tdx_y_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 18);

        auto tdy_y_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 18);

        auto tdz_y_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 18);

        auto tdx_y_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 19);

        auto tdy_y_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 19);

        auto tdz_y_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 19);

        auto tdx_y_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 20);

        auto tdy_y_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 20);

        auto tdz_y_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 20);

        auto tdx_y_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 21);

        auto tdy_y_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 21);

        auto tdz_y_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 21);

        auto tdx_y_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 22);

        auto tdy_y_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 22);

        auto tdz_y_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 22);

        auto tdx_y_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 23);

        auto tdy_y_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 23);

        auto tdz_y_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 23);

        auto tdx_y_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 24);

        auto tdy_y_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 24);

        auto tdz_y_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 24);

        auto tdx_y_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 25);

        auto tdy_y_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 25);

        auto tdz_y_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 25);

        auto tdx_y_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 26);

        auto tdy_y_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 26);

        auto tdz_y_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 26);

        auto tdx_y_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 27);

        auto tdy_y_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 27);

        auto tdz_y_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 27);

        auto tdx_y_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 28);

        auto tdy_y_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 28);

        auto tdz_y_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 28);

        auto tdx_y_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 29);

        auto tdy_y_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 29);

        auto tdz_y_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 29);

        auto tdx_z_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 30);

        auto tdy_z_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 30);

        auto tdz_z_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 30);

        auto tdx_z_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 31);

        auto tdy_z_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 31);

        auto tdz_z_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 31);

        auto tdx_z_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 32);

        auto tdy_z_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 32);

        auto tdz_z_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 32);

        auto tdx_z_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 33);

        auto tdz_xy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 11);

        auto tdx_xy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 12);

        auto tdy_xy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 12);

        auto tdz_xy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 12);

        auto tdx_xy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 13);

        auto tdy_xy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 13);

        auto tdz_xy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 13);

        auto tdx_xy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 14);

        auto tdy_xy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 14);

        auto tdz_xy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 14);

        auto tdx_xy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 15);

        auto tdy_xy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 15);

        auto tdz_xy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 15);

        auto tdx_xy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 16);

        auto tdy_xy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 16);

        auto tdz_xy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 16);

        auto tdx_xy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 17);

        auto tdy_xy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 17);

        auto tdz_xy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 17);

        auto tdx_xy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 18);

        auto tdy_xy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 18);

        auto tdz_xy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 18);

        auto tdx_xy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 19);

        auto tdy_xy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 19);

        auto tdz_xy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 19);

        auto tdx_xz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 20);

        auto tdy_xz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 20);

        auto tdz_xz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 20);

        auto tdx_xz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 21);

        auto tdy_xz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 21);

        auto tdz_xz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 21);

        auto tdx_xz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 22);

        auto tdy_xz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 22);

        auto tdz_xz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 22);

        auto tdx_xz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 23);

        auto ts_xy_xxxz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 17);

        auto ts_xy_xxyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 18);

        auto ts_xy_xxyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 19);

        auto ts_xy_xxzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 20);

        auto ts_xy_xyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 21);

        auto ts_xy_xyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 22);

        auto ts_xy_xyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 23);

        auto ts_xy_xzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 24);

        auto ts_xy_yyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 25);

        auto ts_xy_yyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 26);

        auto ts_xy_yyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 27);

        auto ts_xy_yzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 28);

        auto ts_xy_zzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 29);

        auto ts_xz_xxxx_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 30);

        auto ts_xz_xxxy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 31);

        auto ts_xz_xxxz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 32);

        auto ts_xz_xxyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 33);

        // set up pointers to integrals

        auto tdz_xxy_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 16);

        auto tdx_xxy_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 17);

        auto tdy_xxy_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 17);

        auto tdz_xxy_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 17);

        auto tdx_xxy_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 18);

        auto tdy_xxy_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 18);

        auto tdz_xxy_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 18);

        auto tdx_xxy_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 19);

        auto tdy_xxy_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 19);

        auto tdz_xxy_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 19);

        auto tdx_xxy_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 20);

        auto tdy_xxy_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 20);

        auto tdz_xxy_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 20);

        auto tdx_xxy_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 21);

        auto tdy_xxy_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 21);

        auto tdz_xxy_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 21);

        auto tdx_xxy_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 22);

        auto tdy_xxy_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 22);

        auto tdz_xxy_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 22);

        auto tdx_xxy_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 23);

        auto tdy_xxy_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 23);

        auto tdz_xxy_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 23);

        auto tdx_xxy_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 24);

        auto tdy_xxy_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 24);

        auto tdz_xxy_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 24);

        auto tdx_xxy_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 25);

        auto tdy_xxy_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 25);

        auto tdz_xxy_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 25);

        auto tdx_xxy_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 26);

        auto tdy_xxy_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 26);

        auto tdz_xxy_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 26);

        auto tdx_xxy_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 27);

        auto tdy_xxy_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 27);

        auto tdz_xxy_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 27);

        auto tdx_xxy_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 28);

        auto tdy_xxy_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 28);

        auto tdz_xxy_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 28);

        auto tdx_xxy_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 29);

        auto tdy_xxy_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 29);

        auto tdz_xxy_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 29);

        auto tdx_xxz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 30);

        auto tdy_xxz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 30);

        auto tdz_xxz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 30);

        auto tdx_xxz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 31);

        auto tdy_xxz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 31);

        auto tdz_xxz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 31);

        auto tdx_xxz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 32);

        auto tdy_xxz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 32);

        auto tdz_xxz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 32);

        auto tdx_xxz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 33);

        // Batch of Integrals (50,100)

        #pragma omp simd aligned(fx, pa_x, tdx_xxy_xxxz_0, tdx_xxy_xxyy_0, tdx_xxy_xxyz_0, \
                                     tdx_xxy_xxzz_0, tdx_xxy_xyyy_0, tdx_xxy_xyyz_0, tdx_xxy_xyzz_0, tdx_xxy_xzzz_0, \
                                     tdx_xxy_yyyy_0, tdx_xxy_yyyz_0, tdx_xxy_yyzz_0, tdx_xxy_yzzz_0, tdx_xxy_zzzz_0, \
                                     tdx_xxz_xxxx_0, tdx_xxz_xxxy_0, tdx_xxz_xxxz_0, tdx_xxz_xxyy_0, tdx_xy_xxxz_0, \
                                     tdx_xy_xxyy_0, tdx_xy_xxyz_0, tdx_xy_xxz_0, tdx_xy_xxzz_0, tdx_xy_xyy_0, \
                                     tdx_xy_xyyy_0, tdx_xy_xyyz_0, tdx_xy_xyz_0, tdx_xy_xyzz_0, tdx_xy_xzz_0, \
                                     tdx_xy_xzzz_0, tdx_xy_yyy_0, tdx_xy_yyyy_0, tdx_xy_yyyz_0, tdx_xy_yyz_0, \
                                     tdx_xy_yyzz_0, tdx_xy_yzz_0, tdx_xy_yzzz_0, tdx_xy_zzz_0, tdx_xy_zzzz_0, \
                                     tdx_xz_xxx_0, tdx_xz_xxxx_0, tdx_xz_xxxy_0, tdx_xz_xxxz_0, tdx_xz_xxy_0, \
                                     tdx_xz_xxyy_0, tdx_xz_xxz_0, tdx_xz_xyy_0, tdx_y_xxxz_0, tdx_y_xxyy_0, tdx_y_xxyz_0, \
                                     tdx_y_xxzz_0, tdx_y_xyyy_0, tdx_y_xyyz_0, tdx_y_xyzz_0, tdx_y_xzzz_0, tdx_y_yyyy_0, \
                                     tdx_y_yyyz_0, tdx_y_yyzz_0, tdx_y_yzzz_0, tdx_y_zzzz_0, tdx_z_xxxx_0, tdx_z_xxxy_0, \
                                     tdx_z_xxxz_0, tdx_z_xxyy_0, tdy_xxy_xxxz_0, tdy_xxy_xxyy_0, tdy_xxy_xxyz_0, \
                                     tdy_xxy_xxzz_0, tdy_xxy_xyyy_0, tdy_xxy_xyyz_0, tdy_xxy_xyzz_0, tdy_xxy_xzzz_0, \
                                     tdy_xxy_yyyy_0, tdy_xxy_yyyz_0, tdy_xxy_yyzz_0, tdy_xxy_yzzz_0, tdy_xxy_zzzz_0, \
                                     tdy_xxz_xxxx_0, tdy_xxz_xxxy_0, tdy_xxz_xxxz_0, tdy_xy_xxxz_0, tdy_xy_xxyy_0, \
                                     tdy_xy_xxyz_0, tdy_xy_xxz_0, tdy_xy_xxzz_0, tdy_xy_xyy_0, tdy_xy_xyyy_0, \
                                     tdy_xy_xyyz_0, tdy_xy_xyz_0, tdy_xy_xyzz_0, tdy_xy_xzz_0, tdy_xy_xzzz_0, \
                                     tdy_xy_yyy_0, tdy_xy_yyyy_0, tdy_xy_yyyz_0, tdy_xy_yyz_0, tdy_xy_yyzz_0, \
                                     tdy_xy_yzz_0, tdy_xy_yzzz_0, tdy_xy_zzz_0, tdy_xy_zzzz_0, tdy_xz_xxx_0, \
                                     tdy_xz_xxxx_0, tdy_xz_xxxy_0, tdy_xz_xxxz_0, tdy_xz_xxy_0, tdy_xz_xxz_0, \
                                     tdy_y_xxxz_0, tdy_y_xxyy_0, tdy_y_xxyz_0, tdy_y_xxzz_0, tdy_y_xyyy_0, tdy_y_xyyz_0, \
                                     tdy_y_xyzz_0, tdy_y_xzzz_0, tdy_y_yyyy_0, tdy_y_yyyz_0, tdy_y_yyzz_0, tdy_y_yzzz_0, \
                                     tdy_y_zzzz_0, tdy_z_xxxx_0, tdy_z_xxxy_0, tdy_z_xxxz_0, tdz_xxy_xxxy_0, \
                                     tdz_xxy_xxxz_0, tdz_xxy_xxyy_0, tdz_xxy_xxyz_0, tdz_xxy_xxzz_0, tdz_xxy_xyyy_0, \
                                     tdz_xxy_xyyz_0, tdz_xxy_xyzz_0, tdz_xxy_xzzz_0, tdz_xxy_yyyy_0, tdz_xxy_yyyz_0, \
                                     tdz_xxy_yyzz_0, tdz_xxy_yzzz_0, tdz_xxy_zzzz_0, tdz_xxz_xxxx_0, tdz_xxz_xxxy_0, \
                                     tdz_xxz_xxxz_0, tdz_xy_xxxy_0, tdz_xy_xxxz_0, tdz_xy_xxy_0, tdz_xy_xxyy_0, \
                                     tdz_xy_xxyz_0, tdz_xy_xxz_0, tdz_xy_xxzz_0, tdz_xy_xyy_0, tdz_xy_xyyy_0, \
                                     tdz_xy_xyyz_0, tdz_xy_xyz_0, tdz_xy_xyzz_0, tdz_xy_xzz_0, tdz_xy_xzzz_0, \
                                     tdz_xy_yyy_0, tdz_xy_yyyy_0, tdz_xy_yyyz_0, tdz_xy_yyz_0, tdz_xy_yyzz_0, \
                                     tdz_xy_yzz_0, tdz_xy_yzzz_0, tdz_xy_zzz_0, tdz_xy_zzzz_0, tdz_xz_xxx_0, \
                                     tdz_xz_xxxx_0, tdz_xz_xxxy_0, tdz_xz_xxxz_0, tdz_xz_xxy_0, tdz_xz_xxz_0, \
                                     tdz_y_xxxy_0, tdz_y_xxxz_0, tdz_y_xxyy_0, tdz_y_xxyz_0, tdz_y_xxzz_0, tdz_y_xyyy_0, \
                                     tdz_y_xyyz_0, tdz_y_xyzz_0, tdz_y_xzzz_0, tdz_y_yyyy_0, tdz_y_yyyz_0, tdz_y_yyzz_0, \
                                     tdz_y_yzzz_0, tdz_y_zzzz_0, tdz_z_xxxx_0, tdz_z_xxxy_0, tdz_z_xxxz_0, ts_xy_xxxz_0, \
                                     ts_xy_xxyy_0, ts_xy_xxyz_0, ts_xy_xxzz_0, ts_xy_xyyy_0, ts_xy_xyyz_0, ts_xy_xyzz_0, \
                                     ts_xy_xzzz_0, ts_xy_yyyy_0, ts_xy_yyyz_0, ts_xy_yyzz_0, ts_xy_yzzz_0, ts_xy_zzzz_0, \
                                     ts_xz_xxxx_0, ts_xz_xxxy_0, ts_xz_xxxz_0, ts_xz_xxyy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            tdz_xxy_xxxy_0[j] = pa_x[j] * tdz_xy_xxxy_0[j] + 0.5 * fl1_fx * tdz_y_xxxy_0[j] + 1.5 * fl1_fx * tdz_xy_xxy_0[j];

            tdx_xxy_xxxz_0[j] =
                pa_x[j] * tdx_xy_xxxz_0[j] + 0.5 * fl1_fx * tdx_y_xxxz_0[j] + 1.5 * fl1_fx * tdx_xy_xxz_0[j] + 0.5 * fl1_fx * ts_xy_xxxz_0[j];

            tdy_xxy_xxxz_0[j] = pa_x[j] * tdy_xy_xxxz_0[j] + 0.5 * fl1_fx * tdy_y_xxxz_0[j] + 1.5 * fl1_fx * tdy_xy_xxz_0[j];

            tdz_xxy_xxxz_0[j] = pa_x[j] * tdz_xy_xxxz_0[j] + 0.5 * fl1_fx * tdz_y_xxxz_0[j] + 1.5 * fl1_fx * tdz_xy_xxz_0[j];

            tdx_xxy_xxyy_0[j] =
                pa_x[j] * tdx_xy_xxyy_0[j] + 0.5 * fl1_fx * tdx_y_xxyy_0[j] + fl1_fx * tdx_xy_xyy_0[j] + 0.5 * fl1_fx * ts_xy_xxyy_0[j];

            tdy_xxy_xxyy_0[j] = pa_x[j] * tdy_xy_xxyy_0[j] + 0.5 * fl1_fx * tdy_y_xxyy_0[j] + fl1_fx * tdy_xy_xyy_0[j];

            tdz_xxy_xxyy_0[j] = pa_x[j] * tdz_xy_xxyy_0[j] + 0.5 * fl1_fx * tdz_y_xxyy_0[j] + fl1_fx * tdz_xy_xyy_0[j];

            tdx_xxy_xxyz_0[j] =
                pa_x[j] * tdx_xy_xxyz_0[j] + 0.5 * fl1_fx * tdx_y_xxyz_0[j] + fl1_fx * tdx_xy_xyz_0[j] + 0.5 * fl1_fx * ts_xy_xxyz_0[j];

            tdy_xxy_xxyz_0[j] = pa_x[j] * tdy_xy_xxyz_0[j] + 0.5 * fl1_fx * tdy_y_xxyz_0[j] + fl1_fx * tdy_xy_xyz_0[j];

            tdz_xxy_xxyz_0[j] = pa_x[j] * tdz_xy_xxyz_0[j] + 0.5 * fl1_fx * tdz_y_xxyz_0[j] + fl1_fx * tdz_xy_xyz_0[j];

            tdx_xxy_xxzz_0[j] =
                pa_x[j] * tdx_xy_xxzz_0[j] + 0.5 * fl1_fx * tdx_y_xxzz_0[j] + fl1_fx * tdx_xy_xzz_0[j] + 0.5 * fl1_fx * ts_xy_xxzz_0[j];

            tdy_xxy_xxzz_0[j] = pa_x[j] * tdy_xy_xxzz_0[j] + 0.5 * fl1_fx * tdy_y_xxzz_0[j] + fl1_fx * tdy_xy_xzz_0[j];

            tdz_xxy_xxzz_0[j] = pa_x[j] * tdz_xy_xxzz_0[j] + 0.5 * fl1_fx * tdz_y_xxzz_0[j] + fl1_fx * tdz_xy_xzz_0[j];

            tdx_xxy_xyyy_0[j] =
                pa_x[j] * tdx_xy_xyyy_0[j] + 0.5 * fl1_fx * tdx_y_xyyy_0[j] + 0.5 * fl1_fx * tdx_xy_yyy_0[j] + 0.5 * fl1_fx * ts_xy_xyyy_0[j];

            tdy_xxy_xyyy_0[j] = pa_x[j] * tdy_xy_xyyy_0[j] + 0.5 * fl1_fx * tdy_y_xyyy_0[j] + 0.5 * fl1_fx * tdy_xy_yyy_0[j];

            tdz_xxy_xyyy_0[j] = pa_x[j] * tdz_xy_xyyy_0[j] + 0.5 * fl1_fx * tdz_y_xyyy_0[j] + 0.5 * fl1_fx * tdz_xy_yyy_0[j];

            tdx_xxy_xyyz_0[j] =
                pa_x[j] * tdx_xy_xyyz_0[j] + 0.5 * fl1_fx * tdx_y_xyyz_0[j] + 0.5 * fl1_fx * tdx_xy_yyz_0[j] + 0.5 * fl1_fx * ts_xy_xyyz_0[j];

            tdy_xxy_xyyz_0[j] = pa_x[j] * tdy_xy_xyyz_0[j] + 0.5 * fl1_fx * tdy_y_xyyz_0[j] + 0.5 * fl1_fx * tdy_xy_yyz_0[j];

            tdz_xxy_xyyz_0[j] = pa_x[j] * tdz_xy_xyyz_0[j] + 0.5 * fl1_fx * tdz_y_xyyz_0[j] + 0.5 * fl1_fx * tdz_xy_yyz_0[j];

            tdx_xxy_xyzz_0[j] =
                pa_x[j] * tdx_xy_xyzz_0[j] + 0.5 * fl1_fx * tdx_y_xyzz_0[j] + 0.5 * fl1_fx * tdx_xy_yzz_0[j] + 0.5 * fl1_fx * ts_xy_xyzz_0[j];

            tdy_xxy_xyzz_0[j] = pa_x[j] * tdy_xy_xyzz_0[j] + 0.5 * fl1_fx * tdy_y_xyzz_0[j] + 0.5 * fl1_fx * tdy_xy_yzz_0[j];

            tdz_xxy_xyzz_0[j] = pa_x[j] * tdz_xy_xyzz_0[j] + 0.5 * fl1_fx * tdz_y_xyzz_0[j] + 0.5 * fl1_fx * tdz_xy_yzz_0[j];

            tdx_xxy_xzzz_0[j] =
                pa_x[j] * tdx_xy_xzzz_0[j] + 0.5 * fl1_fx * tdx_y_xzzz_0[j] + 0.5 * fl1_fx * tdx_xy_zzz_0[j] + 0.5 * fl1_fx * ts_xy_xzzz_0[j];

            tdy_xxy_xzzz_0[j] = pa_x[j] * tdy_xy_xzzz_0[j] + 0.5 * fl1_fx * tdy_y_xzzz_0[j] + 0.5 * fl1_fx * tdy_xy_zzz_0[j];

            tdz_xxy_xzzz_0[j] = pa_x[j] * tdz_xy_xzzz_0[j] + 0.5 * fl1_fx * tdz_y_xzzz_0[j] + 0.5 * fl1_fx * tdz_xy_zzz_0[j];

            tdx_xxy_yyyy_0[j] = pa_x[j] * tdx_xy_yyyy_0[j] + 0.5 * fl1_fx * tdx_y_yyyy_0[j] + 0.5 * fl1_fx * ts_xy_yyyy_0[j];

            tdy_xxy_yyyy_0[j] = pa_x[j] * tdy_xy_yyyy_0[j] + 0.5 * fl1_fx * tdy_y_yyyy_0[j];

            tdz_xxy_yyyy_0[j] = pa_x[j] * tdz_xy_yyyy_0[j] + 0.5 * fl1_fx * tdz_y_yyyy_0[j];

            tdx_xxy_yyyz_0[j] = pa_x[j] * tdx_xy_yyyz_0[j] + 0.5 * fl1_fx * tdx_y_yyyz_0[j] + 0.5 * fl1_fx * ts_xy_yyyz_0[j];

            tdy_xxy_yyyz_0[j] = pa_x[j] * tdy_xy_yyyz_0[j] + 0.5 * fl1_fx * tdy_y_yyyz_0[j];

            tdz_xxy_yyyz_0[j] = pa_x[j] * tdz_xy_yyyz_0[j] + 0.5 * fl1_fx * tdz_y_yyyz_0[j];

            tdx_xxy_yyzz_0[j] = pa_x[j] * tdx_xy_yyzz_0[j] + 0.5 * fl1_fx * tdx_y_yyzz_0[j] + 0.5 * fl1_fx * ts_xy_yyzz_0[j];

            tdy_xxy_yyzz_0[j] = pa_x[j] * tdy_xy_yyzz_0[j] + 0.5 * fl1_fx * tdy_y_yyzz_0[j];

            tdz_xxy_yyzz_0[j] = pa_x[j] * tdz_xy_yyzz_0[j] + 0.5 * fl1_fx * tdz_y_yyzz_0[j];

            tdx_xxy_yzzz_0[j] = pa_x[j] * tdx_xy_yzzz_0[j] + 0.5 * fl1_fx * tdx_y_yzzz_0[j] + 0.5 * fl1_fx * ts_xy_yzzz_0[j];

            tdy_xxy_yzzz_0[j] = pa_x[j] * tdy_xy_yzzz_0[j] + 0.5 * fl1_fx * tdy_y_yzzz_0[j];

            tdz_xxy_yzzz_0[j] = pa_x[j] * tdz_xy_yzzz_0[j] + 0.5 * fl1_fx * tdz_y_yzzz_0[j];

            tdx_xxy_zzzz_0[j] = pa_x[j] * tdx_xy_zzzz_0[j] + 0.5 * fl1_fx * tdx_y_zzzz_0[j] + 0.5 * fl1_fx * ts_xy_zzzz_0[j];

            tdy_xxy_zzzz_0[j] = pa_x[j] * tdy_xy_zzzz_0[j] + 0.5 * fl1_fx * tdy_y_zzzz_0[j];

            tdz_xxy_zzzz_0[j] = pa_x[j] * tdz_xy_zzzz_0[j] + 0.5 * fl1_fx * tdz_y_zzzz_0[j];

            tdx_xxz_xxxx_0[j] =
                pa_x[j] * tdx_xz_xxxx_0[j] + 0.5 * fl1_fx * tdx_z_xxxx_0[j] + 2.0 * fl1_fx * tdx_xz_xxx_0[j] + 0.5 * fl1_fx * ts_xz_xxxx_0[j];

            tdy_xxz_xxxx_0[j] = pa_x[j] * tdy_xz_xxxx_0[j] + 0.5 * fl1_fx * tdy_z_xxxx_0[j] + 2.0 * fl1_fx * tdy_xz_xxx_0[j];

            tdz_xxz_xxxx_0[j] = pa_x[j] * tdz_xz_xxxx_0[j] + 0.5 * fl1_fx * tdz_z_xxxx_0[j] + 2.0 * fl1_fx * tdz_xz_xxx_0[j];

            tdx_xxz_xxxy_0[j] =
                pa_x[j] * tdx_xz_xxxy_0[j] + 0.5 * fl1_fx * tdx_z_xxxy_0[j] + 1.5 * fl1_fx * tdx_xz_xxy_0[j] + 0.5 * fl1_fx * ts_xz_xxxy_0[j];

            tdy_xxz_xxxy_0[j] = pa_x[j] * tdy_xz_xxxy_0[j] + 0.5 * fl1_fx * tdy_z_xxxy_0[j] + 1.5 * fl1_fx * tdy_xz_xxy_0[j];

            tdz_xxz_xxxy_0[j] = pa_x[j] * tdz_xz_xxxy_0[j] + 0.5 * fl1_fx * tdz_z_xxxy_0[j] + 1.5 * fl1_fx * tdz_xz_xxy_0[j];

            tdx_xxz_xxxz_0[j] =
                pa_x[j] * tdx_xz_xxxz_0[j] + 0.5 * fl1_fx * tdx_z_xxxz_0[j] + 1.5 * fl1_fx * tdx_xz_xxz_0[j] + 0.5 * fl1_fx * ts_xz_xxxz_0[j];

            tdy_xxz_xxxz_0[j] = pa_x[j] * tdy_xz_xxxz_0[j] + 0.5 * fl1_fx * tdy_z_xxxz_0[j] + 1.5 * fl1_fx * tdy_xz_xxz_0[j];

            tdz_xxz_xxxz_0[j] = pa_x[j] * tdz_xz_xxxz_0[j] + 0.5 * fl1_fx * tdz_z_xxxz_0[j] + 1.5 * fl1_fx * tdz_xz_xxz_0[j];

            tdx_xxz_xxyy_0[j] =
                pa_x[j] * tdx_xz_xxyy_0[j] + 0.5 * fl1_fx * tdx_z_xxyy_0[j] + fl1_fx * tdx_xz_xyy_0[j] + 0.5 * fl1_fx * ts_xz_xxyy_0[j];
        }

        idx++;
    }
}

void
compElectricDipoleForFG_100_150(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const int32_t              nOSFactors,
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

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_d_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

        auto tdy_xz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 33);

        auto tdz_xz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 33);

        auto tdx_xz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 34);

        auto tdy_xz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 34);

        auto tdz_xz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 34);

        auto tdx_xz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 35);

        auto tdy_xz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 35);

        auto tdz_xz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 35);

        auto tdx_xz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 36);

        auto tdy_xz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 36);

        auto tdz_xz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 36);

        auto tdx_xz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 37);

        auto tdy_xz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 37);

        auto tdz_xz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 37);

        auto tdx_xz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 38);

        auto tdy_xz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 38);

        auto tdz_xz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 38);

        auto tdx_xz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 39);

        auto tdy_xz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 39);

        auto tdz_xz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 39);

        auto tdx_xz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 40);

        auto tdy_xz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 40);

        auto tdz_xz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 40);

        auto tdx_xz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 41);

        auto tdy_xz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 41);

        auto tdz_xz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 41);

        auto tdx_xz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 42);

        auto tdy_xz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 42);

        auto tdz_xz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 42);

        auto tdx_xz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 43);

        auto tdy_xz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 43);

        auto tdz_xz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 43);

        auto tdx_xz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 44);

        auto tdy_xz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 44);

        auto tdz_xz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 44);

        auto tdx_yy_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 45);

        auto tdy_yy_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 45);

        auto tdz_yy_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 45);

        auto tdx_yy_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 46);

        auto tdy_yy_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 46);

        auto tdz_yy_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 46);

        auto tdx_yy_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 47);

        auto tdy_yy_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 47);

        auto tdz_yy_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 47);

        auto tdx_yy_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 48);

        auto tdy_yy_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 48);

        auto tdz_yy_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 48);

        auto tdx_yy_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 49);

        auto tdy_yy_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 49);

        auto tdz_yy_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 49);

        auto tdy_z_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 33);

        auto tdz_z_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 33);

        auto tdx_z_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 34);

        auto tdy_z_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 34);

        auto tdz_z_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 34);

        auto tdx_z_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 35);

        auto tdy_z_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 35);

        auto tdz_z_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 35);

        auto tdx_z_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 36);

        auto tdy_z_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 36);

        auto tdz_z_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 36);

        auto tdx_z_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 37);

        auto tdy_z_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 37);

        auto tdz_z_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 37);

        auto tdx_z_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 38);

        auto tdy_z_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 38);

        auto tdz_z_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 38);

        auto tdx_z_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 39);

        auto tdy_z_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 39);

        auto tdz_z_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 39);

        auto tdx_z_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 40);

        auto tdy_z_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 40);

        auto tdz_z_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 40);

        auto tdx_z_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 41);

        auto tdy_z_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 41);

        auto tdz_z_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 41);

        auto tdx_z_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 42);

        auto tdy_z_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 42);

        auto tdz_z_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 42);

        auto tdx_z_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 43);

        auto tdy_z_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 43);

        auto tdz_z_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 43);

        auto tdx_z_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 44);

        auto tdy_z_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 44);

        auto tdz_z_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 44);

        auto tdy_xz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 23);

        auto tdz_xz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 23);

        auto tdx_xz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 24);

        auto tdy_xz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 24);

        auto tdz_xz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 24);

        auto tdx_xz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 25);

        auto tdy_xz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 25);

        auto tdz_xz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 25);

        auto tdx_xz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 26);

        auto tdy_xz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 26);

        auto tdz_xz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 26);

        auto tdx_xz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 27);

        auto tdy_xz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 27);

        auto tdz_xz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 27);

        auto tdx_xz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 28);

        auto tdy_xz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 28);

        auto tdz_xz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 28);

        auto tdx_xz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 29);

        auto tdy_xz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 29);

        auto tdz_xz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 29);

        auto tdx_yy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 30);

        auto tdy_yy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 30);

        auto tdz_yy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 30);

        auto tdx_yy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 31);

        auto tdy_yy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 31);

        auto tdz_yy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 31);

        auto tdx_yy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 32);

        auto tdy_yy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 32);

        auto tdz_yy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 32);

        auto tdx_yy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 33);

        auto tdy_yy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 33);

        auto tdz_yy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 33);

        auto tdx_yy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 34);

        auto tdy_yy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 34);

        auto tdz_yy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 34);

        auto ts_xz_xxyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 34);

        auto ts_xz_xxzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 35);

        auto ts_xz_xyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 36);

        auto ts_xz_xyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 37);

        auto ts_xz_xyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 38);

        auto ts_xz_xzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 39);

        auto ts_xz_yyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 40);

        auto ts_xz_yyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 41);

        auto ts_xz_yyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 42);

        auto ts_xz_yzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 43);

        auto ts_xz_zzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 44);

        auto ts_yy_xxxx_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 45);

        auto ts_yy_xxxy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 46);

        auto ts_yy_xxxz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 47);

        auto ts_yy_xxyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 48);

        auto ts_yy_xxyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 49);

        // set up pointers to integrals

        auto tdy_xxz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 33);

        auto tdz_xxz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 33);

        auto tdx_xxz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 34);

        auto tdy_xxz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 34);

        auto tdz_xxz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 34);

        auto tdx_xxz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 35);

        auto tdy_xxz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 35);

        auto tdz_xxz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 35);

        auto tdx_xxz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 36);

        auto tdy_xxz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 36);

        auto tdz_xxz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 36);

        auto tdx_xxz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 37);

        auto tdy_xxz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 37);

        auto tdz_xxz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 37);

        auto tdx_xxz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 38);

        auto tdy_xxz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 38);

        auto tdz_xxz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 38);

        auto tdx_xxz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 39);

        auto tdy_xxz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 39);

        auto tdz_xxz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 39);

        auto tdx_xxz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 40);

        auto tdy_xxz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 40);

        auto tdz_xxz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 40);

        auto tdx_xxz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 41);

        auto tdy_xxz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 41);

        auto tdz_xxz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 41);

        auto tdx_xxz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 42);

        auto tdy_xxz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 42);

        auto tdz_xxz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 42);

        auto tdx_xxz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 43);

        auto tdy_xxz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 43);

        auto tdz_xxz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 43);

        auto tdx_xxz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 44);

        auto tdy_xxz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 44);

        auto tdz_xxz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 44);

        auto tdx_xyy_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 45);

        auto tdy_xyy_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 45);

        auto tdz_xyy_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 45);

        auto tdx_xyy_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 46);

        auto tdy_xyy_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 46);

        auto tdz_xyy_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 46);

        auto tdx_xyy_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 47);

        auto tdy_xyy_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 47);

        auto tdz_xyy_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 47);

        auto tdx_xyy_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 48);

        auto tdy_xyy_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 48);

        auto tdz_xyy_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 48);

        auto tdx_xyy_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 49);

        auto tdy_xyy_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 49);

        auto tdz_xyy_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 49);

        // Batch of Integrals (100,150)

        #pragma omp simd aligned(fx, pa_x, tdx_xxz_xxyz_0, tdx_xxz_xxzz_0, tdx_xxz_xyyy_0, \
                                     tdx_xxz_xyyz_0, tdx_xxz_xyzz_0, tdx_xxz_xzzz_0, tdx_xxz_yyyy_0, tdx_xxz_yyyz_0, \
                                     tdx_xxz_yyzz_0, tdx_xxz_yzzz_0, tdx_xxz_zzzz_0, tdx_xyy_xxxx_0, tdx_xyy_xxxy_0, \
                                     tdx_xyy_xxxz_0, tdx_xyy_xxyy_0, tdx_xyy_xxyz_0, tdx_xz_xxyz_0, tdx_xz_xxzz_0, \
                                     tdx_xz_xyyy_0, tdx_xz_xyyz_0, tdx_xz_xyz_0, tdx_xz_xyzz_0, tdx_xz_xzz_0, \
                                     tdx_xz_xzzz_0, tdx_xz_yyy_0, tdx_xz_yyyy_0, tdx_xz_yyyz_0, tdx_xz_yyz_0, \
                                     tdx_xz_yyzz_0, tdx_xz_yzz_0, tdx_xz_yzzz_0, tdx_xz_zzz_0, tdx_xz_zzzz_0, \
                                     tdx_yy_xxx_0, tdx_yy_xxxx_0, tdx_yy_xxxy_0, tdx_yy_xxxz_0, tdx_yy_xxy_0, \
                                     tdx_yy_xxyy_0, tdx_yy_xxyz_0, tdx_yy_xxz_0, tdx_yy_xyy_0, tdx_yy_xyz_0, \
                                     tdx_z_xxyz_0, tdx_z_xxzz_0, tdx_z_xyyy_0, tdx_z_xyyz_0, tdx_z_xyzz_0, tdx_z_xzzz_0, \
                                     tdx_z_yyyy_0, tdx_z_yyyz_0, tdx_z_yyzz_0, tdx_z_yzzz_0, tdx_z_zzzz_0, \
                                     tdy_xxz_xxyy_0, tdy_xxz_xxyz_0, tdy_xxz_xxzz_0, tdy_xxz_xyyy_0, tdy_xxz_xyyz_0, \
                                     tdy_xxz_xyzz_0, tdy_xxz_xzzz_0, tdy_xxz_yyyy_0, tdy_xxz_yyyz_0, tdy_xxz_yyzz_0, \
                                     tdy_xxz_yzzz_0, tdy_xxz_zzzz_0, tdy_xyy_xxxx_0, tdy_xyy_xxxy_0, tdy_xyy_xxxz_0, \
                                     tdy_xyy_xxyy_0, tdy_xyy_xxyz_0, tdy_xz_xxyy_0, tdy_xz_xxyz_0, tdy_xz_xxzz_0, \
                                     tdy_xz_xyy_0, tdy_xz_xyyy_0, tdy_xz_xyyz_0, tdy_xz_xyz_0, tdy_xz_xyzz_0, \
                                     tdy_xz_xzz_0, tdy_xz_xzzz_0, tdy_xz_yyy_0, tdy_xz_yyyy_0, tdy_xz_yyyz_0, \
                                     tdy_xz_yyz_0, tdy_xz_yyzz_0, tdy_xz_yzz_0, tdy_xz_yzzz_0, tdy_xz_zzz_0, \
                                     tdy_xz_zzzz_0, tdy_yy_xxx_0, tdy_yy_xxxx_0, tdy_yy_xxxy_0, tdy_yy_xxxz_0, \
                                     tdy_yy_xxy_0, tdy_yy_xxyy_0, tdy_yy_xxyz_0, tdy_yy_xxz_0, tdy_yy_xyy_0, \
                                     tdy_yy_xyz_0, tdy_z_xxyy_0, tdy_z_xxyz_0, tdy_z_xxzz_0, tdy_z_xyyy_0, tdy_z_xyyz_0, \
                                     tdy_z_xyzz_0, tdy_z_xzzz_0, tdy_z_yyyy_0, tdy_z_yyyz_0, tdy_z_yyzz_0, tdy_z_yzzz_0, \
                                     tdy_z_zzzz_0, tdz_xxz_xxyy_0, tdz_xxz_xxyz_0, tdz_xxz_xxzz_0, tdz_xxz_xyyy_0, \
                                     tdz_xxz_xyyz_0, tdz_xxz_xyzz_0, tdz_xxz_xzzz_0, tdz_xxz_yyyy_0, tdz_xxz_yyyz_0, \
                                     tdz_xxz_yyzz_0, tdz_xxz_yzzz_0, tdz_xxz_zzzz_0, tdz_xyy_xxxx_0, tdz_xyy_xxxy_0, \
                                     tdz_xyy_xxxz_0, tdz_xyy_xxyy_0, tdz_xyy_xxyz_0, tdz_xz_xxyy_0, tdz_xz_xxyz_0, \
                                     tdz_xz_xxzz_0, tdz_xz_xyy_0, tdz_xz_xyyy_0, tdz_xz_xyyz_0, tdz_xz_xyz_0, \
                                     tdz_xz_xyzz_0, tdz_xz_xzz_0, tdz_xz_xzzz_0, tdz_xz_yyy_0, tdz_xz_yyyy_0, \
                                     tdz_xz_yyyz_0, tdz_xz_yyz_0, tdz_xz_yyzz_0, tdz_xz_yzz_0, tdz_xz_yzzz_0, \
                                     tdz_xz_zzz_0, tdz_xz_zzzz_0, tdz_yy_xxx_0, tdz_yy_xxxx_0, tdz_yy_xxxy_0, \
                                     tdz_yy_xxxz_0, tdz_yy_xxy_0, tdz_yy_xxyy_0, tdz_yy_xxyz_0, tdz_yy_xxz_0, \
                                     tdz_yy_xyy_0, tdz_yy_xyz_0, tdz_z_xxyy_0, tdz_z_xxyz_0, tdz_z_xxzz_0, tdz_z_xyyy_0, \
                                     tdz_z_xyyz_0, tdz_z_xyzz_0, tdz_z_xzzz_0, tdz_z_yyyy_0, tdz_z_yyyz_0, tdz_z_yyzz_0, \
                                     tdz_z_yzzz_0, tdz_z_zzzz_0, ts_xz_xxyz_0, ts_xz_xxzz_0, ts_xz_xyyy_0, ts_xz_xyyz_0, \
                                     ts_xz_xyzz_0, ts_xz_xzzz_0, ts_xz_yyyy_0, ts_xz_yyyz_0, ts_xz_yyzz_0, ts_xz_yzzz_0, \
                                     ts_xz_zzzz_0, ts_yy_xxxx_0, ts_yy_xxxy_0, ts_yy_xxxz_0, ts_yy_xxyy_0, ts_yy_xxyz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            tdy_xxz_xxyy_0[j] = pa_x[j] * tdy_xz_xxyy_0[j] + 0.5 * fl1_fx * tdy_z_xxyy_0[j] + fl1_fx * tdy_xz_xyy_0[j];

            tdz_xxz_xxyy_0[j] = pa_x[j] * tdz_xz_xxyy_0[j] + 0.5 * fl1_fx * tdz_z_xxyy_0[j] + fl1_fx * tdz_xz_xyy_0[j];

            tdx_xxz_xxyz_0[j] =
                pa_x[j] * tdx_xz_xxyz_0[j] + 0.5 * fl1_fx * tdx_z_xxyz_0[j] + fl1_fx * tdx_xz_xyz_0[j] + 0.5 * fl1_fx * ts_xz_xxyz_0[j];

            tdy_xxz_xxyz_0[j] = pa_x[j] * tdy_xz_xxyz_0[j] + 0.5 * fl1_fx * tdy_z_xxyz_0[j] + fl1_fx * tdy_xz_xyz_0[j];

            tdz_xxz_xxyz_0[j] = pa_x[j] * tdz_xz_xxyz_0[j] + 0.5 * fl1_fx * tdz_z_xxyz_0[j] + fl1_fx * tdz_xz_xyz_0[j];

            tdx_xxz_xxzz_0[j] =
                pa_x[j] * tdx_xz_xxzz_0[j] + 0.5 * fl1_fx * tdx_z_xxzz_0[j] + fl1_fx * tdx_xz_xzz_0[j] + 0.5 * fl1_fx * ts_xz_xxzz_0[j];

            tdy_xxz_xxzz_0[j] = pa_x[j] * tdy_xz_xxzz_0[j] + 0.5 * fl1_fx * tdy_z_xxzz_0[j] + fl1_fx * tdy_xz_xzz_0[j];

            tdz_xxz_xxzz_0[j] = pa_x[j] * tdz_xz_xxzz_0[j] + 0.5 * fl1_fx * tdz_z_xxzz_0[j] + fl1_fx * tdz_xz_xzz_0[j];

            tdx_xxz_xyyy_0[j] =
                pa_x[j] * tdx_xz_xyyy_0[j] + 0.5 * fl1_fx * tdx_z_xyyy_0[j] + 0.5 * fl1_fx * tdx_xz_yyy_0[j] + 0.5 * fl1_fx * ts_xz_xyyy_0[j];

            tdy_xxz_xyyy_0[j] = pa_x[j] * tdy_xz_xyyy_0[j] + 0.5 * fl1_fx * tdy_z_xyyy_0[j] + 0.5 * fl1_fx * tdy_xz_yyy_0[j];

            tdz_xxz_xyyy_0[j] = pa_x[j] * tdz_xz_xyyy_0[j] + 0.5 * fl1_fx * tdz_z_xyyy_0[j] + 0.5 * fl1_fx * tdz_xz_yyy_0[j];

            tdx_xxz_xyyz_0[j] =
                pa_x[j] * tdx_xz_xyyz_0[j] + 0.5 * fl1_fx * tdx_z_xyyz_0[j] + 0.5 * fl1_fx * tdx_xz_yyz_0[j] + 0.5 * fl1_fx * ts_xz_xyyz_0[j];

            tdy_xxz_xyyz_0[j] = pa_x[j] * tdy_xz_xyyz_0[j] + 0.5 * fl1_fx * tdy_z_xyyz_0[j] + 0.5 * fl1_fx * tdy_xz_yyz_0[j];

            tdz_xxz_xyyz_0[j] = pa_x[j] * tdz_xz_xyyz_0[j] + 0.5 * fl1_fx * tdz_z_xyyz_0[j] + 0.5 * fl1_fx * tdz_xz_yyz_0[j];

            tdx_xxz_xyzz_0[j] =
                pa_x[j] * tdx_xz_xyzz_0[j] + 0.5 * fl1_fx * tdx_z_xyzz_0[j] + 0.5 * fl1_fx * tdx_xz_yzz_0[j] + 0.5 * fl1_fx * ts_xz_xyzz_0[j];

            tdy_xxz_xyzz_0[j] = pa_x[j] * tdy_xz_xyzz_0[j] + 0.5 * fl1_fx * tdy_z_xyzz_0[j] + 0.5 * fl1_fx * tdy_xz_yzz_0[j];

            tdz_xxz_xyzz_0[j] = pa_x[j] * tdz_xz_xyzz_0[j] + 0.5 * fl1_fx * tdz_z_xyzz_0[j] + 0.5 * fl1_fx * tdz_xz_yzz_0[j];

            tdx_xxz_xzzz_0[j] =
                pa_x[j] * tdx_xz_xzzz_0[j] + 0.5 * fl1_fx * tdx_z_xzzz_0[j] + 0.5 * fl1_fx * tdx_xz_zzz_0[j] + 0.5 * fl1_fx * ts_xz_xzzz_0[j];

            tdy_xxz_xzzz_0[j] = pa_x[j] * tdy_xz_xzzz_0[j] + 0.5 * fl1_fx * tdy_z_xzzz_0[j] + 0.5 * fl1_fx * tdy_xz_zzz_0[j];

            tdz_xxz_xzzz_0[j] = pa_x[j] * tdz_xz_xzzz_0[j] + 0.5 * fl1_fx * tdz_z_xzzz_0[j] + 0.5 * fl1_fx * tdz_xz_zzz_0[j];

            tdx_xxz_yyyy_0[j] = pa_x[j] * tdx_xz_yyyy_0[j] + 0.5 * fl1_fx * tdx_z_yyyy_0[j] + 0.5 * fl1_fx * ts_xz_yyyy_0[j];

            tdy_xxz_yyyy_0[j] = pa_x[j] * tdy_xz_yyyy_0[j] + 0.5 * fl1_fx * tdy_z_yyyy_0[j];

            tdz_xxz_yyyy_0[j] = pa_x[j] * tdz_xz_yyyy_0[j] + 0.5 * fl1_fx * tdz_z_yyyy_0[j];

            tdx_xxz_yyyz_0[j] = pa_x[j] * tdx_xz_yyyz_0[j] + 0.5 * fl1_fx * tdx_z_yyyz_0[j] + 0.5 * fl1_fx * ts_xz_yyyz_0[j];

            tdy_xxz_yyyz_0[j] = pa_x[j] * tdy_xz_yyyz_0[j] + 0.5 * fl1_fx * tdy_z_yyyz_0[j];

            tdz_xxz_yyyz_0[j] = pa_x[j] * tdz_xz_yyyz_0[j] + 0.5 * fl1_fx * tdz_z_yyyz_0[j];

            tdx_xxz_yyzz_0[j] = pa_x[j] * tdx_xz_yyzz_0[j] + 0.5 * fl1_fx * tdx_z_yyzz_0[j] + 0.5 * fl1_fx * ts_xz_yyzz_0[j];

            tdy_xxz_yyzz_0[j] = pa_x[j] * tdy_xz_yyzz_0[j] + 0.5 * fl1_fx * tdy_z_yyzz_0[j];

            tdz_xxz_yyzz_0[j] = pa_x[j] * tdz_xz_yyzz_0[j] + 0.5 * fl1_fx * tdz_z_yyzz_0[j];

            tdx_xxz_yzzz_0[j] = pa_x[j] * tdx_xz_yzzz_0[j] + 0.5 * fl1_fx * tdx_z_yzzz_0[j] + 0.5 * fl1_fx * ts_xz_yzzz_0[j];

            tdy_xxz_yzzz_0[j] = pa_x[j] * tdy_xz_yzzz_0[j] + 0.5 * fl1_fx * tdy_z_yzzz_0[j];

            tdz_xxz_yzzz_0[j] = pa_x[j] * tdz_xz_yzzz_0[j] + 0.5 * fl1_fx * tdz_z_yzzz_0[j];

            tdx_xxz_zzzz_0[j] = pa_x[j] * tdx_xz_zzzz_0[j] + 0.5 * fl1_fx * tdx_z_zzzz_0[j] + 0.5 * fl1_fx * ts_xz_zzzz_0[j];

            tdy_xxz_zzzz_0[j] = pa_x[j] * tdy_xz_zzzz_0[j] + 0.5 * fl1_fx * tdy_z_zzzz_0[j];

            tdz_xxz_zzzz_0[j] = pa_x[j] * tdz_xz_zzzz_0[j] + 0.5 * fl1_fx * tdz_z_zzzz_0[j];

            tdx_xyy_xxxx_0[j] = pa_x[j] * tdx_yy_xxxx_0[j] + 2.0 * fl1_fx * tdx_yy_xxx_0[j] + 0.5 * fl1_fx * ts_yy_xxxx_0[j];

            tdy_xyy_xxxx_0[j] = pa_x[j] * tdy_yy_xxxx_0[j] + 2.0 * fl1_fx * tdy_yy_xxx_0[j];

            tdz_xyy_xxxx_0[j] = pa_x[j] * tdz_yy_xxxx_0[j] + 2.0 * fl1_fx * tdz_yy_xxx_0[j];

            tdx_xyy_xxxy_0[j] = pa_x[j] * tdx_yy_xxxy_0[j] + 1.5 * fl1_fx * tdx_yy_xxy_0[j] + 0.5 * fl1_fx * ts_yy_xxxy_0[j];

            tdy_xyy_xxxy_0[j] = pa_x[j] * tdy_yy_xxxy_0[j] + 1.5 * fl1_fx * tdy_yy_xxy_0[j];

            tdz_xyy_xxxy_0[j] = pa_x[j] * tdz_yy_xxxy_0[j] + 1.5 * fl1_fx * tdz_yy_xxy_0[j];

            tdx_xyy_xxxz_0[j] = pa_x[j] * tdx_yy_xxxz_0[j] + 1.5 * fl1_fx * tdx_yy_xxz_0[j] + 0.5 * fl1_fx * ts_yy_xxxz_0[j];

            tdy_xyy_xxxz_0[j] = pa_x[j] * tdy_yy_xxxz_0[j] + 1.5 * fl1_fx * tdy_yy_xxz_0[j];

            tdz_xyy_xxxz_0[j] = pa_x[j] * tdz_yy_xxxz_0[j] + 1.5 * fl1_fx * tdz_yy_xxz_0[j];

            tdx_xyy_xxyy_0[j] = pa_x[j] * tdx_yy_xxyy_0[j] + fl1_fx * tdx_yy_xyy_0[j] + 0.5 * fl1_fx * ts_yy_xxyy_0[j];

            tdy_xyy_xxyy_0[j] = pa_x[j] * tdy_yy_xxyy_0[j] + fl1_fx * tdy_yy_xyy_0[j];

            tdz_xyy_xxyy_0[j] = pa_x[j] * tdz_yy_xxyy_0[j] + fl1_fx * tdz_yy_xyy_0[j];

            tdx_xyy_xxyz_0[j] = pa_x[j] * tdx_yy_xxyz_0[j] + fl1_fx * tdx_yy_xyz_0[j] + 0.5 * fl1_fx * ts_yy_xxyz_0[j];

            tdy_xyy_xxyz_0[j] = pa_x[j] * tdy_yy_xxyz_0[j] + fl1_fx * tdy_yy_xyz_0[j];

            tdz_xyy_xxyz_0[j] = pa_x[j] * tdz_yy_xxyz_0[j] + fl1_fx * tdz_yy_xyz_0[j];
        }

        idx++;
    }
}

void
compElectricDipoleForFG_150_200(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const int32_t              nOSFactors,
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

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_d_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

        auto tdx_yy_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 50);

        auto tdy_yy_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 50);

        auto tdz_yy_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 50);

        auto tdx_yy_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 51);

        auto tdy_yy_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 51);

        auto tdz_yy_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 51);

        auto tdx_yy_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 52);

        auto tdy_yy_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 52);

        auto tdz_yy_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 52);

        auto tdx_yy_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 53);

        auto tdy_yy_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 53);

        auto tdz_yy_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 53);

        auto tdx_yy_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 54);

        auto tdy_yy_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 54);

        auto tdz_yy_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 54);

        auto tdx_yy_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 55);

        auto tdy_yy_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 55);

        auto tdz_yy_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 55);

        auto tdx_yy_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 56);

        auto tdy_yy_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 56);

        auto tdz_yy_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 56);

        auto tdx_yy_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 57);

        auto tdy_yy_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 57);

        auto tdz_yy_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 57);

        auto tdx_yy_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 58);

        auto tdy_yy_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 58);

        auto tdz_yy_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 58);

        auto tdx_yy_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 59);

        auto tdy_yy_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 59);

        auto tdz_yy_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 59);

        auto tdx_yz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 60);

        auto tdy_yz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 60);

        auto tdz_yz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 60);

        auto tdx_yz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 61);

        auto tdy_yz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 61);

        auto tdz_yz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 61);

        auto tdx_yz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 62);

        auto tdy_yz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 62);

        auto tdz_yz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 62);

        auto tdx_yz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 63);

        auto tdy_yz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 63);

        auto tdz_yz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 63);

        auto tdx_yz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 64);

        auto tdy_yz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 64);

        auto tdz_yz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 64);

        auto tdx_yz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 65);

        auto tdy_yz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 65);

        auto tdz_yz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 65);

        auto tdx_yz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 66);

        auto tdy_yz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 66);

        auto tdx_yy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 35);

        auto tdy_yy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 35);

        auto tdz_yy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 35);

        auto tdx_yy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 36);

        auto tdy_yy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 36);

        auto tdz_yy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 36);

        auto tdx_yy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 37);

        auto tdy_yy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 37);

        auto tdz_yy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 37);

        auto tdx_yy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 38);

        auto tdy_yy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 38);

        auto tdz_yy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 38);

        auto tdx_yy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 39);

        auto tdy_yy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 39);

        auto tdz_yy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 39);

        auto tdx_yz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 40);

        auto tdy_yz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 40);

        auto tdz_yz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 40);

        auto tdx_yz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 41);

        auto tdy_yz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 41);

        auto tdz_yz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 41);

        auto tdx_yz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 42);

        auto tdy_yz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 42);

        auto tdz_yz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 42);

        auto tdx_yz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 43);

        auto tdy_yz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 43);

        auto tdz_yz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 43);

        auto tdx_yz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 44);

        auto tdy_yz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 44);

        auto tdz_yz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 44);

        auto tdx_yz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 45);

        auto tdy_yz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 45);

        auto tdz_yz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 45);

        auto tdx_yz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 46);

        auto tdy_yz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 46);

        auto ts_yy_xxzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 50);

        auto ts_yy_xyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 51);

        auto ts_yy_xyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 52);

        auto ts_yy_xyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 53);

        auto ts_yy_xzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 54);

        auto ts_yy_yyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 55);

        auto ts_yy_yyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 56);

        auto ts_yy_yyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 57);

        auto ts_yy_yzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 58);

        auto ts_yy_zzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 59);

        auto ts_yz_xxxx_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 60);

        auto ts_yz_xxxy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 61);

        auto ts_yz_xxxz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 62);

        auto ts_yz_xxyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 63);

        auto ts_yz_xxyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 64);

        auto ts_yz_xxzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 65);

        auto ts_yz_xyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 66);

        // set up pointers to integrals

        auto tdx_xyy_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 50);

        auto tdy_xyy_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 50);

        auto tdz_xyy_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 50);

        auto tdx_xyy_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 51);

        auto tdy_xyy_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 51);

        auto tdz_xyy_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 51);

        auto tdx_xyy_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 52);

        auto tdy_xyy_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 52);

        auto tdz_xyy_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 52);

        auto tdx_xyy_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 53);

        auto tdy_xyy_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 53);

        auto tdz_xyy_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 53);

        auto tdx_xyy_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 54);

        auto tdy_xyy_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 54);

        auto tdz_xyy_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 54);

        auto tdx_xyy_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 55);

        auto tdy_xyy_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 55);

        auto tdz_xyy_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 55);

        auto tdx_xyy_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 56);

        auto tdy_xyy_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 56);

        auto tdz_xyy_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 56);

        auto tdx_xyy_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 57);

        auto tdy_xyy_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 57);

        auto tdz_xyy_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 57);

        auto tdx_xyy_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 58);

        auto tdy_xyy_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 58);

        auto tdz_xyy_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 58);

        auto tdx_xyy_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 59);

        auto tdy_xyy_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 59);

        auto tdz_xyy_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 59);

        auto tdx_xyz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 60);

        auto tdy_xyz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 60);

        auto tdz_xyz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 60);

        auto tdx_xyz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 61);

        auto tdy_xyz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 61);

        auto tdz_xyz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 61);

        auto tdx_xyz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 62);

        auto tdy_xyz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 62);

        auto tdz_xyz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 62);

        auto tdx_xyz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 63);

        auto tdy_xyz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 63);

        auto tdz_xyz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 63);

        auto tdx_xyz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 64);

        auto tdy_xyz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 64);

        auto tdz_xyz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 64);

        auto tdx_xyz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 65);

        auto tdy_xyz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 65);

        auto tdz_xyz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 65);

        auto tdx_xyz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 66);

        auto tdy_xyz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 66);

        // Batch of Integrals (150,200)

        #pragma omp simd aligned(fx, pa_x, tdx_xyy_xxzz_0, tdx_xyy_xyyy_0, tdx_xyy_xyyz_0, \
                                     tdx_xyy_xyzz_0, tdx_xyy_xzzz_0, tdx_xyy_yyyy_0, tdx_xyy_yyyz_0, tdx_xyy_yyzz_0, \
                                     tdx_xyy_yzzz_0, tdx_xyy_zzzz_0, tdx_xyz_xxxx_0, tdx_xyz_xxxy_0, tdx_xyz_xxxz_0, \
                                     tdx_xyz_xxyy_0, tdx_xyz_xxyz_0, tdx_xyz_xxzz_0, tdx_xyz_xyyy_0, tdx_yy_xxzz_0, \
                                     tdx_yy_xyyy_0, tdx_yy_xyyz_0, tdx_yy_xyzz_0, tdx_yy_xzz_0, tdx_yy_xzzz_0, \
                                     tdx_yy_yyy_0, tdx_yy_yyyy_0, tdx_yy_yyyz_0, tdx_yy_yyz_0, tdx_yy_yyzz_0, \
                                     tdx_yy_yzz_0, tdx_yy_yzzz_0, tdx_yy_zzz_0, tdx_yy_zzzz_0, tdx_yz_xxx_0, \
                                     tdx_yz_xxxx_0, tdx_yz_xxxy_0, tdx_yz_xxxz_0, tdx_yz_xxy_0, tdx_yz_xxyy_0, \
                                     tdx_yz_xxyz_0, tdx_yz_xxz_0, tdx_yz_xxzz_0, tdx_yz_xyy_0, tdx_yz_xyyy_0, \
                                     tdx_yz_xyz_0, tdx_yz_xzz_0, tdx_yz_yyy_0, tdy_xyy_xxzz_0, tdy_xyy_xyyy_0, \
                                     tdy_xyy_xyyz_0, tdy_xyy_xyzz_0, tdy_xyy_xzzz_0, tdy_xyy_yyyy_0, tdy_xyy_yyyz_0, \
                                     tdy_xyy_yyzz_0, tdy_xyy_yzzz_0, tdy_xyy_zzzz_0, tdy_xyz_xxxx_0, tdy_xyz_xxxy_0, \
                                     tdy_xyz_xxxz_0, tdy_xyz_xxyy_0, tdy_xyz_xxyz_0, tdy_xyz_xxzz_0, tdy_xyz_xyyy_0, \
                                     tdy_yy_xxzz_0, tdy_yy_xyyy_0, tdy_yy_xyyz_0, tdy_yy_xyzz_0, tdy_yy_xzz_0, \
                                     tdy_yy_xzzz_0, tdy_yy_yyy_0, tdy_yy_yyyy_0, tdy_yy_yyyz_0, tdy_yy_yyz_0, \
                                     tdy_yy_yyzz_0, tdy_yy_yzz_0, tdy_yy_yzzz_0, tdy_yy_zzz_0, tdy_yy_zzzz_0, \
                                     tdy_yz_xxx_0, tdy_yz_xxxx_0, tdy_yz_xxxy_0, tdy_yz_xxxz_0, tdy_yz_xxy_0, \
                                     tdy_yz_xxyy_0, tdy_yz_xxyz_0, tdy_yz_xxz_0, tdy_yz_xxzz_0, tdy_yz_xyy_0, \
                                     tdy_yz_xyyy_0, tdy_yz_xyz_0, tdy_yz_xzz_0, tdy_yz_yyy_0, tdz_xyy_xxzz_0, \
                                     tdz_xyy_xyyy_0, tdz_xyy_xyyz_0, tdz_xyy_xyzz_0, tdz_xyy_xzzz_0, tdz_xyy_yyyy_0, \
                                     tdz_xyy_yyyz_0, tdz_xyy_yyzz_0, tdz_xyy_yzzz_0, tdz_xyy_zzzz_0, tdz_xyz_xxxx_0, \
                                     tdz_xyz_xxxy_0, tdz_xyz_xxxz_0, tdz_xyz_xxyy_0, tdz_xyz_xxyz_0, tdz_xyz_xxzz_0, \
                                     tdz_yy_xxzz_0, tdz_yy_xyyy_0, tdz_yy_xyyz_0, tdz_yy_xyzz_0, tdz_yy_xzz_0, \
                                     tdz_yy_xzzz_0, tdz_yy_yyy_0, tdz_yy_yyyy_0, tdz_yy_yyyz_0, tdz_yy_yyz_0, \
                                     tdz_yy_yyzz_0, tdz_yy_yzz_0, tdz_yy_yzzz_0, tdz_yy_zzz_0, tdz_yy_zzzz_0, \
                                     tdz_yz_xxx_0, tdz_yz_xxxx_0, tdz_yz_xxxy_0, tdz_yz_xxxz_0, tdz_yz_xxy_0, \
                                     tdz_yz_xxyy_0, tdz_yz_xxyz_0, tdz_yz_xxz_0, tdz_yz_xxzz_0, tdz_yz_xyy_0, \
                                     tdz_yz_xyz_0, tdz_yz_xzz_0, ts_yy_xxzz_0, ts_yy_xyyy_0, ts_yy_xyyz_0, ts_yy_xyzz_0, \
                                     ts_yy_xzzz_0, ts_yy_yyyy_0, ts_yy_yyyz_0, ts_yy_yyzz_0, ts_yy_yzzz_0, ts_yy_zzzz_0, \
                                     ts_yz_xxxx_0, ts_yz_xxxy_0, ts_yz_xxxz_0, ts_yz_xxyy_0, ts_yz_xxyz_0, ts_yz_xxzz_0, \
                                     ts_yz_xyyy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            tdx_xyy_xxzz_0[j] = pa_x[j] * tdx_yy_xxzz_0[j] + fl1_fx * tdx_yy_xzz_0[j] + 0.5 * fl1_fx * ts_yy_xxzz_0[j];

            tdy_xyy_xxzz_0[j] = pa_x[j] * tdy_yy_xxzz_0[j] + fl1_fx * tdy_yy_xzz_0[j];

            tdz_xyy_xxzz_0[j] = pa_x[j] * tdz_yy_xxzz_0[j] + fl1_fx * tdz_yy_xzz_0[j];

            tdx_xyy_xyyy_0[j] = pa_x[j] * tdx_yy_xyyy_0[j] + 0.5 * fl1_fx * tdx_yy_yyy_0[j] + 0.5 * fl1_fx * ts_yy_xyyy_0[j];

            tdy_xyy_xyyy_0[j] = pa_x[j] * tdy_yy_xyyy_0[j] + 0.5 * fl1_fx * tdy_yy_yyy_0[j];

            tdz_xyy_xyyy_0[j] = pa_x[j] * tdz_yy_xyyy_0[j] + 0.5 * fl1_fx * tdz_yy_yyy_0[j];

            tdx_xyy_xyyz_0[j] = pa_x[j] * tdx_yy_xyyz_0[j] + 0.5 * fl1_fx * tdx_yy_yyz_0[j] + 0.5 * fl1_fx * ts_yy_xyyz_0[j];

            tdy_xyy_xyyz_0[j] = pa_x[j] * tdy_yy_xyyz_0[j] + 0.5 * fl1_fx * tdy_yy_yyz_0[j];

            tdz_xyy_xyyz_0[j] = pa_x[j] * tdz_yy_xyyz_0[j] + 0.5 * fl1_fx * tdz_yy_yyz_0[j];

            tdx_xyy_xyzz_0[j] = pa_x[j] * tdx_yy_xyzz_0[j] + 0.5 * fl1_fx * tdx_yy_yzz_0[j] + 0.5 * fl1_fx * ts_yy_xyzz_0[j];

            tdy_xyy_xyzz_0[j] = pa_x[j] * tdy_yy_xyzz_0[j] + 0.5 * fl1_fx * tdy_yy_yzz_0[j];

            tdz_xyy_xyzz_0[j] = pa_x[j] * tdz_yy_xyzz_0[j] + 0.5 * fl1_fx * tdz_yy_yzz_0[j];

            tdx_xyy_xzzz_0[j] = pa_x[j] * tdx_yy_xzzz_0[j] + 0.5 * fl1_fx * tdx_yy_zzz_0[j] + 0.5 * fl1_fx * ts_yy_xzzz_0[j];

            tdy_xyy_xzzz_0[j] = pa_x[j] * tdy_yy_xzzz_0[j] + 0.5 * fl1_fx * tdy_yy_zzz_0[j];

            tdz_xyy_xzzz_0[j] = pa_x[j] * tdz_yy_xzzz_0[j] + 0.5 * fl1_fx * tdz_yy_zzz_0[j];

            tdx_xyy_yyyy_0[j] = pa_x[j] * tdx_yy_yyyy_0[j] + 0.5 * fl1_fx * ts_yy_yyyy_0[j];

            tdy_xyy_yyyy_0[j] = pa_x[j] * tdy_yy_yyyy_0[j];

            tdz_xyy_yyyy_0[j] = pa_x[j] * tdz_yy_yyyy_0[j];

            tdx_xyy_yyyz_0[j] = pa_x[j] * tdx_yy_yyyz_0[j] + 0.5 * fl1_fx * ts_yy_yyyz_0[j];

            tdy_xyy_yyyz_0[j] = pa_x[j] * tdy_yy_yyyz_0[j];

            tdz_xyy_yyyz_0[j] = pa_x[j] * tdz_yy_yyyz_0[j];

            tdx_xyy_yyzz_0[j] = pa_x[j] * tdx_yy_yyzz_0[j] + 0.5 * fl1_fx * ts_yy_yyzz_0[j];

            tdy_xyy_yyzz_0[j] = pa_x[j] * tdy_yy_yyzz_0[j];

            tdz_xyy_yyzz_0[j] = pa_x[j] * tdz_yy_yyzz_0[j];

            tdx_xyy_yzzz_0[j] = pa_x[j] * tdx_yy_yzzz_0[j] + 0.5 * fl1_fx * ts_yy_yzzz_0[j];

            tdy_xyy_yzzz_0[j] = pa_x[j] * tdy_yy_yzzz_0[j];

            tdz_xyy_yzzz_0[j] = pa_x[j] * tdz_yy_yzzz_0[j];

            tdx_xyy_zzzz_0[j] = pa_x[j] * tdx_yy_zzzz_0[j] + 0.5 * fl1_fx * ts_yy_zzzz_0[j];

            tdy_xyy_zzzz_0[j] = pa_x[j] * tdy_yy_zzzz_0[j];

            tdz_xyy_zzzz_0[j] = pa_x[j] * tdz_yy_zzzz_0[j];

            tdx_xyz_xxxx_0[j] = pa_x[j] * tdx_yz_xxxx_0[j] + 2.0 * fl1_fx * tdx_yz_xxx_0[j] + 0.5 * fl1_fx * ts_yz_xxxx_0[j];

            tdy_xyz_xxxx_0[j] = pa_x[j] * tdy_yz_xxxx_0[j] + 2.0 * fl1_fx * tdy_yz_xxx_0[j];

            tdz_xyz_xxxx_0[j] = pa_x[j] * tdz_yz_xxxx_0[j] + 2.0 * fl1_fx * tdz_yz_xxx_0[j];

            tdx_xyz_xxxy_0[j] = pa_x[j] * tdx_yz_xxxy_0[j] + 1.5 * fl1_fx * tdx_yz_xxy_0[j] + 0.5 * fl1_fx * ts_yz_xxxy_0[j];

            tdy_xyz_xxxy_0[j] = pa_x[j] * tdy_yz_xxxy_0[j] + 1.5 * fl1_fx * tdy_yz_xxy_0[j];

            tdz_xyz_xxxy_0[j] = pa_x[j] * tdz_yz_xxxy_0[j] + 1.5 * fl1_fx * tdz_yz_xxy_0[j];

            tdx_xyz_xxxz_0[j] = pa_x[j] * tdx_yz_xxxz_0[j] + 1.5 * fl1_fx * tdx_yz_xxz_0[j] + 0.5 * fl1_fx * ts_yz_xxxz_0[j];

            tdy_xyz_xxxz_0[j] = pa_x[j] * tdy_yz_xxxz_0[j] + 1.5 * fl1_fx * tdy_yz_xxz_0[j];

            tdz_xyz_xxxz_0[j] = pa_x[j] * tdz_yz_xxxz_0[j] + 1.5 * fl1_fx * tdz_yz_xxz_0[j];

            tdx_xyz_xxyy_0[j] = pa_x[j] * tdx_yz_xxyy_0[j] + fl1_fx * tdx_yz_xyy_0[j] + 0.5 * fl1_fx * ts_yz_xxyy_0[j];

            tdy_xyz_xxyy_0[j] = pa_x[j] * tdy_yz_xxyy_0[j] + fl1_fx * tdy_yz_xyy_0[j];

            tdz_xyz_xxyy_0[j] = pa_x[j] * tdz_yz_xxyy_0[j] + fl1_fx * tdz_yz_xyy_0[j];

            tdx_xyz_xxyz_0[j] = pa_x[j] * tdx_yz_xxyz_0[j] + fl1_fx * tdx_yz_xyz_0[j] + 0.5 * fl1_fx * ts_yz_xxyz_0[j];

            tdy_xyz_xxyz_0[j] = pa_x[j] * tdy_yz_xxyz_0[j] + fl1_fx * tdy_yz_xyz_0[j];

            tdz_xyz_xxyz_0[j] = pa_x[j] * tdz_yz_xxyz_0[j] + fl1_fx * tdz_yz_xyz_0[j];

            tdx_xyz_xxzz_0[j] = pa_x[j] * tdx_yz_xxzz_0[j] + fl1_fx * tdx_yz_xzz_0[j] + 0.5 * fl1_fx * ts_yz_xxzz_0[j];

            tdy_xyz_xxzz_0[j] = pa_x[j] * tdy_yz_xxzz_0[j] + fl1_fx * tdy_yz_xzz_0[j];

            tdz_xyz_xxzz_0[j] = pa_x[j] * tdz_yz_xxzz_0[j] + fl1_fx * tdz_yz_xzz_0[j];

            tdx_xyz_xyyy_0[j] = pa_x[j] * tdx_yz_xyyy_0[j] + 0.5 * fl1_fx * tdx_yz_yyy_0[j] + 0.5 * fl1_fx * ts_yz_xyyy_0[j];

            tdy_xyz_xyyy_0[j] = pa_x[j] * tdy_yz_xyyy_0[j] + 0.5 * fl1_fx * tdy_yz_yyy_0[j];
        }

        idx++;
    }
}

void
compElectricDipoleForFG_200_250(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const int32_t              nOSFactors,
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

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_d_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

        auto tdz_yz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 66);

        auto tdx_yz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 67);

        auto tdy_yz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 67);

        auto tdz_yz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 67);

        auto tdx_yz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 68);

        auto tdy_yz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 68);

        auto tdz_yz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 68);

        auto tdx_yz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 69);

        auto tdy_yz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 69);

        auto tdz_yz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 69);

        auto tdx_yz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 70);

        auto tdy_yz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 70);

        auto tdz_yz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 70);

        auto tdx_yz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 71);

        auto tdy_yz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 71);

        auto tdz_yz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 71);

        auto tdx_yz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 72);

        auto tdy_yz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 72);

        auto tdz_yz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 72);

        auto tdx_yz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 73);

        auto tdy_yz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 73);

        auto tdz_yz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 73);

        auto tdx_yz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 74);

        auto tdy_yz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 74);

        auto tdz_yz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 74);

        auto tdx_zz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 75);

        auto tdy_zz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 75);

        auto tdz_zz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 75);

        auto tdx_zz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 76);

        auto tdy_zz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 76);

        auto tdz_zz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 76);

        auto tdx_zz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 77);

        auto tdy_zz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 77);

        auto tdz_zz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 77);

        auto tdx_zz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 78);

        auto tdy_zz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 78);

        auto tdz_zz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 78);

        auto tdx_zz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 79);

        auto tdy_zz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 79);

        auto tdz_zz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 79);

        auto tdx_zz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 80);

        auto tdy_zz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 80);

        auto tdz_zz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 80);

        auto tdx_zz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 81);

        auto tdy_zz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 81);

        auto tdz_zz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 81);

        auto tdx_zz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 82);

        auto tdy_zz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 82);

        auto tdz_zz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 82);

        auto tdx_zz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 83);

        auto tdz_yz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 46);

        auto tdx_yz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 47);

        auto tdy_yz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 47);

        auto tdz_yz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 47);

        auto tdx_yz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 48);

        auto tdy_yz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 48);

        auto tdz_yz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 48);

        auto tdx_yz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 49);

        auto tdy_yz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 49);

        auto tdz_yz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 49);

        auto tdx_zz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 50);

        auto tdy_zz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 50);

        auto tdz_zz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 50);

        auto tdx_zz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 51);

        auto tdy_zz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 51);

        auto tdz_zz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 51);

        auto tdx_zz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 52);

        auto tdy_zz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 52);

        auto tdz_zz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 52);

        auto tdx_zz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 53);

        auto tdy_zz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 53);

        auto tdz_zz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 53);

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

        auto ts_yz_xyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 67);

        auto ts_yz_xyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 68);

        auto ts_yz_xzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 69);

        auto ts_yz_yyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 70);

        auto ts_yz_yyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 71);

        auto ts_yz_yyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 72);

        auto ts_yz_yzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 73);

        auto ts_yz_zzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 74);

        auto ts_zz_xxxx_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 75);

        auto ts_zz_xxxy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 76);

        auto ts_zz_xxxz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 77);

        auto ts_zz_xxyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 78);

        auto ts_zz_xxyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 79);

        auto ts_zz_xxzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 80);

        auto ts_zz_xyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 81);

        auto ts_zz_xyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 82);

        auto ts_zz_xyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 83);

        // set up pointers to integrals

        auto tdz_xyz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 66);

        auto tdx_xyz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 67);

        auto tdy_xyz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 67);

        auto tdz_xyz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 67);

        auto tdx_xyz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 68);

        auto tdy_xyz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 68);

        auto tdz_xyz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 68);

        auto tdx_xyz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 69);

        auto tdy_xyz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 69);

        auto tdz_xyz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 69);

        auto tdx_xyz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 70);

        auto tdy_xyz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 70);

        auto tdz_xyz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 70);

        auto tdx_xyz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 71);

        auto tdy_xyz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 71);

        auto tdz_xyz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 71);

        auto tdx_xyz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 72);

        auto tdy_xyz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 72);

        auto tdz_xyz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 72);

        auto tdx_xyz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 73);

        auto tdy_xyz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 73);

        auto tdz_xyz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 73);

        auto tdx_xyz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 74);

        auto tdy_xyz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 74);

        auto tdz_xyz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 74);

        auto tdx_xzz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 75);

        auto tdy_xzz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 75);

        auto tdz_xzz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 75);

        auto tdx_xzz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 76);

        auto tdy_xzz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 76);

        auto tdz_xzz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 76);

        auto tdx_xzz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 77);

        auto tdy_xzz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 77);

        auto tdz_xzz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 77);

        auto tdx_xzz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 78);

        auto tdy_xzz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 78);

        auto tdz_xzz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 78);

        auto tdx_xzz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 79);

        auto tdy_xzz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 79);

        auto tdz_xzz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 79);

        auto tdx_xzz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 80);

        auto tdy_xzz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 80);

        auto tdz_xzz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 80);

        auto tdx_xzz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 81);

        auto tdy_xzz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 81);

        auto tdz_xzz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 81);

        auto tdx_xzz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 82);

        auto tdy_xzz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 82);

        auto tdz_xzz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 82);

        auto tdx_xzz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 83);

        // Batch of Integrals (200,250)

        #pragma omp simd aligned(fx, pa_x, tdx_xyz_xyyz_0, tdx_xyz_xyzz_0, tdx_xyz_xzzz_0, \
                                     tdx_xyz_yyyy_0, tdx_xyz_yyyz_0, tdx_xyz_yyzz_0, tdx_xyz_yzzz_0, tdx_xyz_zzzz_0, \
                                     tdx_xzz_xxxx_0, tdx_xzz_xxxy_0, tdx_xzz_xxxz_0, tdx_xzz_xxyy_0, tdx_xzz_xxyz_0, \
                                     tdx_xzz_xxzz_0, tdx_xzz_xyyy_0, tdx_xzz_xyyz_0, tdx_xzz_xyzz_0, tdx_yz_xyyz_0, \
                                     tdx_yz_xyzz_0, tdx_yz_xzzz_0, tdx_yz_yyyy_0, tdx_yz_yyyz_0, tdx_yz_yyz_0, \
                                     tdx_yz_yyzz_0, tdx_yz_yzz_0, tdx_yz_yzzz_0, tdx_yz_zzz_0, tdx_yz_zzzz_0, \
                                     tdx_zz_xxx_0, tdx_zz_xxxx_0, tdx_zz_xxxy_0, tdx_zz_xxxz_0, tdx_zz_xxy_0, \
                                     tdx_zz_xxyy_0, tdx_zz_xxyz_0, tdx_zz_xxz_0, tdx_zz_xxzz_0, tdx_zz_xyy_0, \
                                     tdx_zz_xyyy_0, tdx_zz_xyyz_0, tdx_zz_xyz_0, tdx_zz_xyzz_0, tdx_zz_xzz_0, \
                                     tdx_zz_yyy_0, tdx_zz_yyz_0, tdx_zz_yzz_0, tdy_xyz_xyyz_0, tdy_xyz_xyzz_0, \
                                     tdy_xyz_xzzz_0, tdy_xyz_yyyy_0, tdy_xyz_yyyz_0, tdy_xyz_yyzz_0, tdy_xyz_yzzz_0, \
                                     tdy_xyz_zzzz_0, tdy_xzz_xxxx_0, tdy_xzz_xxxy_0, tdy_xzz_xxxz_0, tdy_xzz_xxyy_0, \
                                     tdy_xzz_xxyz_0, tdy_xzz_xxzz_0, tdy_xzz_xyyy_0, tdy_xzz_xyyz_0, tdy_yz_xyyz_0, \
                                     tdy_yz_xyzz_0, tdy_yz_xzzz_0, tdy_yz_yyyy_0, tdy_yz_yyyz_0, tdy_yz_yyz_0, \
                                     tdy_yz_yyzz_0, tdy_yz_yzz_0, tdy_yz_yzzz_0, tdy_yz_zzz_0, tdy_yz_zzzz_0, \
                                     tdy_zz_xxx_0, tdy_zz_xxxx_0, tdy_zz_xxxy_0, tdy_zz_xxxz_0, tdy_zz_xxy_0, \
                                     tdy_zz_xxyy_0, tdy_zz_xxyz_0, tdy_zz_xxz_0, tdy_zz_xxzz_0, tdy_zz_xyy_0, \
                                     tdy_zz_xyyy_0, tdy_zz_xyyz_0, tdy_zz_xyz_0, tdy_zz_xzz_0, tdy_zz_yyy_0, \
                                     tdy_zz_yyz_0, tdz_xyz_xyyy_0, tdz_xyz_xyyz_0, tdz_xyz_xyzz_0, tdz_xyz_xzzz_0, \
                                     tdz_xyz_yyyy_0, tdz_xyz_yyyz_0, tdz_xyz_yyzz_0, tdz_xyz_yzzz_0, tdz_xyz_zzzz_0, \
                                     tdz_xzz_xxxx_0, tdz_xzz_xxxy_0, tdz_xzz_xxxz_0, tdz_xzz_xxyy_0, tdz_xzz_xxyz_0, \
                                     tdz_xzz_xxzz_0, tdz_xzz_xyyy_0, tdz_xzz_xyyz_0, tdz_yz_xyyy_0, tdz_yz_xyyz_0, \
                                     tdz_yz_xyzz_0, tdz_yz_xzzz_0, tdz_yz_yyy_0, tdz_yz_yyyy_0, tdz_yz_yyyz_0, \
                                     tdz_yz_yyz_0, tdz_yz_yyzz_0, tdz_yz_yzz_0, tdz_yz_yzzz_0, tdz_yz_zzz_0, \
                                     tdz_yz_zzzz_0, tdz_zz_xxx_0, tdz_zz_xxxx_0, tdz_zz_xxxy_0, tdz_zz_xxxz_0, \
                                     tdz_zz_xxy_0, tdz_zz_xxyy_0, tdz_zz_xxyz_0, tdz_zz_xxz_0, tdz_zz_xxzz_0, \
                                     tdz_zz_xyy_0, tdz_zz_xyyy_0, tdz_zz_xyyz_0, tdz_zz_xyz_0, tdz_zz_xzz_0, \
                                     tdz_zz_yyy_0, tdz_zz_yyz_0, ts_yz_xyyz_0, ts_yz_xyzz_0, ts_yz_xzzz_0, ts_yz_yyyy_0, \
                                     ts_yz_yyyz_0, ts_yz_yyzz_0, ts_yz_yzzz_0, ts_yz_zzzz_0, ts_zz_xxxx_0, ts_zz_xxxy_0, \
                                     ts_zz_xxxz_0, ts_zz_xxyy_0, ts_zz_xxyz_0, ts_zz_xxzz_0, ts_zz_xyyy_0, ts_zz_xyyz_0, \
                                     ts_zz_xyzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            tdz_xyz_xyyy_0[j] = pa_x[j] * tdz_yz_xyyy_0[j] + 0.5 * fl1_fx * tdz_yz_yyy_0[j];

            tdx_xyz_xyyz_0[j] = pa_x[j] * tdx_yz_xyyz_0[j] + 0.5 * fl1_fx * tdx_yz_yyz_0[j] + 0.5 * fl1_fx * ts_yz_xyyz_0[j];

            tdy_xyz_xyyz_0[j] = pa_x[j] * tdy_yz_xyyz_0[j] + 0.5 * fl1_fx * tdy_yz_yyz_0[j];

            tdz_xyz_xyyz_0[j] = pa_x[j] * tdz_yz_xyyz_0[j] + 0.5 * fl1_fx * tdz_yz_yyz_0[j];

            tdx_xyz_xyzz_0[j] = pa_x[j] * tdx_yz_xyzz_0[j] + 0.5 * fl1_fx * tdx_yz_yzz_0[j] + 0.5 * fl1_fx * ts_yz_xyzz_0[j];

            tdy_xyz_xyzz_0[j] = pa_x[j] * tdy_yz_xyzz_0[j] + 0.5 * fl1_fx * tdy_yz_yzz_0[j];

            tdz_xyz_xyzz_0[j] = pa_x[j] * tdz_yz_xyzz_0[j] + 0.5 * fl1_fx * tdz_yz_yzz_0[j];

            tdx_xyz_xzzz_0[j] = pa_x[j] * tdx_yz_xzzz_0[j] + 0.5 * fl1_fx * tdx_yz_zzz_0[j] + 0.5 * fl1_fx * ts_yz_xzzz_0[j];

            tdy_xyz_xzzz_0[j] = pa_x[j] * tdy_yz_xzzz_0[j] + 0.5 * fl1_fx * tdy_yz_zzz_0[j];

            tdz_xyz_xzzz_0[j] = pa_x[j] * tdz_yz_xzzz_0[j] + 0.5 * fl1_fx * tdz_yz_zzz_0[j];

            tdx_xyz_yyyy_0[j] = pa_x[j] * tdx_yz_yyyy_0[j] + 0.5 * fl1_fx * ts_yz_yyyy_0[j];

            tdy_xyz_yyyy_0[j] = pa_x[j] * tdy_yz_yyyy_0[j];

            tdz_xyz_yyyy_0[j] = pa_x[j] * tdz_yz_yyyy_0[j];

            tdx_xyz_yyyz_0[j] = pa_x[j] * tdx_yz_yyyz_0[j] + 0.5 * fl1_fx * ts_yz_yyyz_0[j];

            tdy_xyz_yyyz_0[j] = pa_x[j] * tdy_yz_yyyz_0[j];

            tdz_xyz_yyyz_0[j] = pa_x[j] * tdz_yz_yyyz_0[j];

            tdx_xyz_yyzz_0[j] = pa_x[j] * tdx_yz_yyzz_0[j] + 0.5 * fl1_fx * ts_yz_yyzz_0[j];

            tdy_xyz_yyzz_0[j] = pa_x[j] * tdy_yz_yyzz_0[j];

            tdz_xyz_yyzz_0[j] = pa_x[j] * tdz_yz_yyzz_0[j];

            tdx_xyz_yzzz_0[j] = pa_x[j] * tdx_yz_yzzz_0[j] + 0.5 * fl1_fx * ts_yz_yzzz_0[j];

            tdy_xyz_yzzz_0[j] = pa_x[j] * tdy_yz_yzzz_0[j];

            tdz_xyz_yzzz_0[j] = pa_x[j] * tdz_yz_yzzz_0[j];

            tdx_xyz_zzzz_0[j] = pa_x[j] * tdx_yz_zzzz_0[j] + 0.5 * fl1_fx * ts_yz_zzzz_0[j];

            tdy_xyz_zzzz_0[j] = pa_x[j] * tdy_yz_zzzz_0[j];

            tdz_xyz_zzzz_0[j] = pa_x[j] * tdz_yz_zzzz_0[j];

            tdx_xzz_xxxx_0[j] = pa_x[j] * tdx_zz_xxxx_0[j] + 2.0 * fl1_fx * tdx_zz_xxx_0[j] + 0.5 * fl1_fx * ts_zz_xxxx_0[j];

            tdy_xzz_xxxx_0[j] = pa_x[j] * tdy_zz_xxxx_0[j] + 2.0 * fl1_fx * tdy_zz_xxx_0[j];

            tdz_xzz_xxxx_0[j] = pa_x[j] * tdz_zz_xxxx_0[j] + 2.0 * fl1_fx * tdz_zz_xxx_0[j];

            tdx_xzz_xxxy_0[j] = pa_x[j] * tdx_zz_xxxy_0[j] + 1.5 * fl1_fx * tdx_zz_xxy_0[j] + 0.5 * fl1_fx * ts_zz_xxxy_0[j];

            tdy_xzz_xxxy_0[j] = pa_x[j] * tdy_zz_xxxy_0[j] + 1.5 * fl1_fx * tdy_zz_xxy_0[j];

            tdz_xzz_xxxy_0[j] = pa_x[j] * tdz_zz_xxxy_0[j] + 1.5 * fl1_fx * tdz_zz_xxy_0[j];

            tdx_xzz_xxxz_0[j] = pa_x[j] * tdx_zz_xxxz_0[j] + 1.5 * fl1_fx * tdx_zz_xxz_0[j] + 0.5 * fl1_fx * ts_zz_xxxz_0[j];

            tdy_xzz_xxxz_0[j] = pa_x[j] * tdy_zz_xxxz_0[j] + 1.5 * fl1_fx * tdy_zz_xxz_0[j];

            tdz_xzz_xxxz_0[j] = pa_x[j] * tdz_zz_xxxz_0[j] + 1.5 * fl1_fx * tdz_zz_xxz_0[j];

            tdx_xzz_xxyy_0[j] = pa_x[j] * tdx_zz_xxyy_0[j] + fl1_fx * tdx_zz_xyy_0[j] + 0.5 * fl1_fx * ts_zz_xxyy_0[j];

            tdy_xzz_xxyy_0[j] = pa_x[j] * tdy_zz_xxyy_0[j] + fl1_fx * tdy_zz_xyy_0[j];

            tdz_xzz_xxyy_0[j] = pa_x[j] * tdz_zz_xxyy_0[j] + fl1_fx * tdz_zz_xyy_0[j];

            tdx_xzz_xxyz_0[j] = pa_x[j] * tdx_zz_xxyz_0[j] + fl1_fx * tdx_zz_xyz_0[j] + 0.5 * fl1_fx * ts_zz_xxyz_0[j];

            tdy_xzz_xxyz_0[j] = pa_x[j] * tdy_zz_xxyz_0[j] + fl1_fx * tdy_zz_xyz_0[j];

            tdz_xzz_xxyz_0[j] = pa_x[j] * tdz_zz_xxyz_0[j] + fl1_fx * tdz_zz_xyz_0[j];

            tdx_xzz_xxzz_0[j] = pa_x[j] * tdx_zz_xxzz_0[j] + fl1_fx * tdx_zz_xzz_0[j] + 0.5 * fl1_fx * ts_zz_xxzz_0[j];

            tdy_xzz_xxzz_0[j] = pa_x[j] * tdy_zz_xxzz_0[j] + fl1_fx * tdy_zz_xzz_0[j];

            tdz_xzz_xxzz_0[j] = pa_x[j] * tdz_zz_xxzz_0[j] + fl1_fx * tdz_zz_xzz_0[j];

            tdx_xzz_xyyy_0[j] = pa_x[j] * tdx_zz_xyyy_0[j] + 0.5 * fl1_fx * tdx_zz_yyy_0[j] + 0.5 * fl1_fx * ts_zz_xyyy_0[j];

            tdy_xzz_xyyy_0[j] = pa_x[j] * tdy_zz_xyyy_0[j] + 0.5 * fl1_fx * tdy_zz_yyy_0[j];

            tdz_xzz_xyyy_0[j] = pa_x[j] * tdz_zz_xyyy_0[j] + 0.5 * fl1_fx * tdz_zz_yyy_0[j];

            tdx_xzz_xyyz_0[j] = pa_x[j] * tdx_zz_xyyz_0[j] + 0.5 * fl1_fx * tdx_zz_yyz_0[j] + 0.5 * fl1_fx * ts_zz_xyyz_0[j];

            tdy_xzz_xyyz_0[j] = pa_x[j] * tdy_zz_xyyz_0[j] + 0.5 * fl1_fx * tdy_zz_yyz_0[j];

            tdz_xzz_xyyz_0[j] = pa_x[j] * tdz_zz_xyyz_0[j] + 0.5 * fl1_fx * tdz_zz_yyz_0[j];

            tdx_xzz_xyzz_0[j] = pa_x[j] * tdx_zz_xyzz_0[j] + 0.5 * fl1_fx * tdx_zz_yzz_0[j] + 0.5 * fl1_fx * ts_zz_xyzz_0[j];
        }

        idx++;
    }
}

void
compElectricDipoleForFG_250_300(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const int32_t              nOSFactors,
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

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_d_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        auto pa_y = paDistances.data(3 * idx + 1);

        // set up pointers to auxilary integrals

        auto tdx_yy_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 45);

        auto tdy_yy_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 45);

        auto tdz_yy_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 45);

        auto tdx_yy_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 46);

        auto tdy_yy_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 46);

        auto tdz_yy_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 46);

        auto tdx_yy_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 47);

        auto tdy_yy_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 47);

        auto tdz_yy_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 47);

        auto tdx_yy_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 48);

        auto tdy_yy_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 48);

        auto tdz_yy_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 48);

        auto tdx_yy_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 49);

        auto tdy_yy_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 49);

        auto tdz_yy_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 49);

        auto tdx_yy_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 50);

        auto tdy_yy_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 50);

        auto tdz_yy_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 50);

        auto tdx_yy_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 51);

        auto tdy_yy_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 51);

        auto tdz_yy_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 51);

        auto tdx_yy_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 52);

        auto tdy_yy_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 52);

        auto tdz_yy_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 52);

        auto tdx_yy_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 53);

        auto tdy_yy_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 53);

        auto tdz_yy_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 53);

        auto tdx_yy_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 54);

        auto tdy_yy_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 54);

        auto tdz_yy_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 54);

        auto tdy_zz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 83);

        auto tdz_zz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 83);

        auto tdx_zz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 84);

        auto tdy_zz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 84);

        auto tdz_zz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 84);

        auto tdx_zz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 85);

        auto tdy_zz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 85);

        auto tdz_zz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 85);

        auto tdx_zz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 86);

        auto tdy_zz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 86);

        auto tdz_zz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 86);

        auto tdx_zz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 87);

        auto tdy_zz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 87);

        auto tdz_zz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 87);

        auto tdx_zz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 88);

        auto tdy_zz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 88);

        auto tdz_zz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 88);

        auto tdx_zz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 89);

        auto tdy_zz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 89);

        auto tdz_zz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 89);

        auto tdx_y_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 15);

        auto tdy_y_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 15);

        auto tdz_y_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 15);

        auto tdx_y_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 16);

        auto tdy_y_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 16);

        auto tdz_y_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 16);

        auto tdx_y_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 17);

        auto tdy_y_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 17);

        auto tdz_y_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 17);

        auto tdx_y_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 18);

        auto tdy_y_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 18);

        auto tdz_y_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 18);

        auto tdx_y_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 19);

        auto tdy_y_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 19);

        auto tdz_y_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 19);

        auto tdx_y_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 20);

        auto tdy_y_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 20);

        auto tdz_y_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 20);

        auto tdx_y_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 21);

        auto tdy_y_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 21);

        auto tdz_y_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 21);

        auto tdx_y_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 22);

        auto tdy_y_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 22);

        auto tdz_y_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 22);

        auto tdx_y_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 23);

        auto tdy_y_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 23);

        auto tdz_y_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 23);

        auto tdx_y_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 24);

        auto tdy_y_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 24);

        auto tdz_y_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 24);

        auto tdx_yy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 30);

        auto tdy_yy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 30);

        auto tdz_yy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 30);

        auto tdx_yy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 31);

        auto tdy_yy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 31);

        auto tdz_yy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 31);

        auto tdx_yy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 32);

        auto tdy_yy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 32);

        auto tdz_yy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 32);

        auto tdx_yy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 33);

        auto tdy_yy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 33);

        auto tdz_yy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 33);

        auto tdx_yy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 34);

        auto tdy_yy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 34);

        auto tdz_yy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 34);

        auto tdx_yy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 35);

        auto tdy_yy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 35);

        auto tdz_yy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 35);

        auto tdy_zz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 58);

        auto tdz_zz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 58);

        auto tdx_zz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 59);

        auto tdy_zz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 59);

        auto tdz_zz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 59);

        auto ts_yy_xxxx_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 45);

        auto ts_yy_xxxy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 46);

        auto ts_yy_xxxz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 47);

        auto ts_yy_xxyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 48);

        auto ts_yy_xxyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 49);

        auto ts_yy_xxzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 50);

        auto ts_yy_xyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 51);

        auto ts_yy_xyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 52);

        auto ts_yy_xyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 53);

        auto ts_yy_xzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 54);

        auto ts_zz_xzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 84);

        auto ts_zz_yyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 85);

        auto ts_zz_yyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 86);

        auto ts_zz_yyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 87);

        auto ts_zz_yzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 88);

        auto ts_zz_zzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 89);

        // set up pointers to integrals

        auto tdy_xzz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 83);

        auto tdz_xzz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 83);

        auto tdx_xzz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 84);

        auto tdy_xzz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 84);

        auto tdz_xzz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 84);

        auto tdx_xzz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 85);

        auto tdy_xzz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 85);

        auto tdz_xzz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 85);

        auto tdx_xzz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 86);

        auto tdy_xzz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 86);

        auto tdz_xzz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 86);

        auto tdx_xzz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 87);

        auto tdy_xzz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 87);

        auto tdz_xzz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 87);

        auto tdx_xzz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 88);

        auto tdy_xzz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 88);

        auto tdz_xzz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 88);

        auto tdx_xzz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 89);

        auto tdy_xzz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 89);

        auto tdz_xzz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 89);

        auto tdx_yyy_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 90);

        auto tdy_yyy_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 90);

        auto tdz_yyy_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 90);

        auto tdx_yyy_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 91);

        auto tdy_yyy_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 91);

        auto tdz_yyy_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 91);

        auto tdx_yyy_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 92);

        auto tdy_yyy_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 92);

        auto tdz_yyy_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 92);

        auto tdx_yyy_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 93);

        auto tdy_yyy_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 93);

        auto tdz_yyy_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 93);

        auto tdx_yyy_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 94);

        auto tdy_yyy_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 94);

        auto tdz_yyy_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 94);

        auto tdx_yyy_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 95);

        auto tdy_yyy_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 95);

        auto tdz_yyy_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 95);

        auto tdx_yyy_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 96);

        auto tdy_yyy_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 96);

        auto tdz_yyy_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 96);

        auto tdx_yyy_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 97);

        auto tdy_yyy_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 97);

        auto tdz_yyy_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 97);

        auto tdx_yyy_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 98);

        auto tdy_yyy_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 98);

        auto tdz_yyy_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 98);

        auto tdx_yyy_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 99);

        auto tdy_yyy_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 99);

        auto tdz_yyy_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 99);

        // Batch of Integrals (250,300)

        #pragma omp simd aligned(fx, pa_x, pa_y, tdx_xzz_xzzz_0, tdx_xzz_yyyy_0, tdx_xzz_yyyz_0, \
                                     tdx_xzz_yyzz_0, tdx_xzz_yzzz_0, tdx_xzz_zzzz_0, tdx_y_xxxx_0, tdx_y_xxxy_0, \
                                     tdx_y_xxxz_0, tdx_y_xxyy_0, tdx_y_xxyz_0, tdx_y_xxzz_0, tdx_y_xyyy_0, tdx_y_xyyz_0, \
                                     tdx_y_xyzz_0, tdx_y_xzzz_0, tdx_yy_xxx_0, tdx_yy_xxxx_0, tdx_yy_xxxy_0, \
                                     tdx_yy_xxxz_0, tdx_yy_xxy_0, tdx_yy_xxyy_0, tdx_yy_xxyz_0, tdx_yy_xxz_0, \
                                     tdx_yy_xxzz_0, tdx_yy_xyy_0, tdx_yy_xyyy_0, tdx_yy_xyyz_0, tdx_yy_xyz_0, \
                                     tdx_yy_xyzz_0, tdx_yy_xzz_0, tdx_yy_xzzz_0, tdx_yyy_xxxx_0, tdx_yyy_xxxy_0, \
                                     tdx_yyy_xxxz_0, tdx_yyy_xxyy_0, tdx_yyy_xxyz_0, tdx_yyy_xxzz_0, tdx_yyy_xyyy_0, \
                                     tdx_yyy_xyyz_0, tdx_yyy_xyzz_0, tdx_yyy_xzzz_0, tdx_zz_xzzz_0, tdx_zz_yyyy_0, \
                                     tdx_zz_yyyz_0, tdx_zz_yyzz_0, tdx_zz_yzzz_0, tdx_zz_zzz_0, tdx_zz_zzzz_0, \
                                     tdy_xzz_xyzz_0, tdy_xzz_xzzz_0, tdy_xzz_yyyy_0, tdy_xzz_yyyz_0, tdy_xzz_yyzz_0, \
                                     tdy_xzz_yzzz_0, tdy_xzz_zzzz_0, tdy_y_xxxx_0, tdy_y_xxxy_0, tdy_y_xxxz_0, \
                                     tdy_y_xxyy_0, tdy_y_xxyz_0, tdy_y_xxzz_0, tdy_y_xyyy_0, tdy_y_xyyz_0, tdy_y_xyzz_0, \
                                     tdy_y_xzzz_0, tdy_yy_xxx_0, tdy_yy_xxxx_0, tdy_yy_xxxy_0, tdy_yy_xxxz_0, \
                                     tdy_yy_xxy_0, tdy_yy_xxyy_0, tdy_yy_xxyz_0, tdy_yy_xxz_0, tdy_yy_xxzz_0, \
                                     tdy_yy_xyy_0, tdy_yy_xyyy_0, tdy_yy_xyyz_0, tdy_yy_xyz_0, tdy_yy_xyzz_0, \
                                     tdy_yy_xzz_0, tdy_yy_xzzz_0, tdy_yyy_xxxx_0, tdy_yyy_xxxy_0, tdy_yyy_xxxz_0, \
                                     tdy_yyy_xxyy_0, tdy_yyy_xxyz_0, tdy_yyy_xxzz_0, tdy_yyy_xyyy_0, tdy_yyy_xyyz_0, \
                                     tdy_yyy_xyzz_0, tdy_yyy_xzzz_0, tdy_zz_xyzz_0, tdy_zz_xzzz_0, tdy_zz_yyyy_0, \
                                     tdy_zz_yyyz_0, tdy_zz_yyzz_0, tdy_zz_yzz_0, tdy_zz_yzzz_0, tdy_zz_zzz_0, \
                                     tdy_zz_zzzz_0, tdz_xzz_xyzz_0, tdz_xzz_xzzz_0, tdz_xzz_yyyy_0, tdz_xzz_yyyz_0, \
                                     tdz_xzz_yyzz_0, tdz_xzz_yzzz_0, tdz_xzz_zzzz_0, tdz_y_xxxx_0, tdz_y_xxxy_0, \
                                     tdz_y_xxxz_0, tdz_y_xxyy_0, tdz_y_xxyz_0, tdz_y_xxzz_0, tdz_y_xyyy_0, tdz_y_xyyz_0, \
                                     tdz_y_xyzz_0, tdz_y_xzzz_0, tdz_yy_xxx_0, tdz_yy_xxxx_0, tdz_yy_xxxy_0, \
                                     tdz_yy_xxxz_0, tdz_yy_xxy_0, tdz_yy_xxyy_0, tdz_yy_xxyz_0, tdz_yy_xxz_0, \
                                     tdz_yy_xxzz_0, tdz_yy_xyy_0, tdz_yy_xyyy_0, tdz_yy_xyyz_0, tdz_yy_xyz_0, \
                                     tdz_yy_xyzz_0, tdz_yy_xzz_0, tdz_yy_xzzz_0, tdz_yyy_xxxx_0, tdz_yyy_xxxy_0, \
                                     tdz_yyy_xxxz_0, tdz_yyy_xxyy_0, tdz_yyy_xxyz_0, tdz_yyy_xxzz_0, tdz_yyy_xyyy_0, \
                                     tdz_yyy_xyyz_0, tdz_yyy_xyzz_0, tdz_yyy_xzzz_0, tdz_zz_xyzz_0, tdz_zz_xzzz_0, \
                                     tdz_zz_yyyy_0, tdz_zz_yyyz_0, tdz_zz_yyzz_0, tdz_zz_yzz_0, tdz_zz_yzzz_0, \
                                     tdz_zz_zzz_0, tdz_zz_zzzz_0, ts_yy_xxxx_0, ts_yy_xxxy_0, ts_yy_xxxz_0, \
                                     ts_yy_xxyy_0, ts_yy_xxyz_0, ts_yy_xxzz_0, ts_yy_xyyy_0, ts_yy_xyyz_0, ts_yy_xyzz_0, \
                                     ts_yy_xzzz_0, ts_zz_xzzz_0, ts_zz_yyyy_0, ts_zz_yyyz_0, ts_zz_yyzz_0, ts_zz_yzzz_0, \
                                     ts_zz_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            tdy_xzz_xyzz_0[j] = pa_x[j] * tdy_zz_xyzz_0[j] + 0.5 * fl1_fx * tdy_zz_yzz_0[j];

            tdz_xzz_xyzz_0[j] = pa_x[j] * tdz_zz_xyzz_0[j] + 0.5 * fl1_fx * tdz_zz_yzz_0[j];

            tdx_xzz_xzzz_0[j] = pa_x[j] * tdx_zz_xzzz_0[j] + 0.5 * fl1_fx * tdx_zz_zzz_0[j] + 0.5 * fl1_fx * ts_zz_xzzz_0[j];

            tdy_xzz_xzzz_0[j] = pa_x[j] * tdy_zz_xzzz_0[j] + 0.5 * fl1_fx * tdy_zz_zzz_0[j];

            tdz_xzz_xzzz_0[j] = pa_x[j] * tdz_zz_xzzz_0[j] + 0.5 * fl1_fx * tdz_zz_zzz_0[j];

            tdx_xzz_yyyy_0[j] = pa_x[j] * tdx_zz_yyyy_0[j] + 0.5 * fl1_fx * ts_zz_yyyy_0[j];

            tdy_xzz_yyyy_0[j] = pa_x[j] * tdy_zz_yyyy_0[j];

            tdz_xzz_yyyy_0[j] = pa_x[j] * tdz_zz_yyyy_0[j];

            tdx_xzz_yyyz_0[j] = pa_x[j] * tdx_zz_yyyz_0[j] + 0.5 * fl1_fx * ts_zz_yyyz_0[j];

            tdy_xzz_yyyz_0[j] = pa_x[j] * tdy_zz_yyyz_0[j];

            tdz_xzz_yyyz_0[j] = pa_x[j] * tdz_zz_yyyz_0[j];

            tdx_xzz_yyzz_0[j] = pa_x[j] * tdx_zz_yyzz_0[j] + 0.5 * fl1_fx * ts_zz_yyzz_0[j];

            tdy_xzz_yyzz_0[j] = pa_x[j] * tdy_zz_yyzz_0[j];

            tdz_xzz_yyzz_0[j] = pa_x[j] * tdz_zz_yyzz_0[j];

            tdx_xzz_yzzz_0[j] = pa_x[j] * tdx_zz_yzzz_0[j] + 0.5 * fl1_fx * ts_zz_yzzz_0[j];

            tdy_xzz_yzzz_0[j] = pa_x[j] * tdy_zz_yzzz_0[j];

            tdz_xzz_yzzz_0[j] = pa_x[j] * tdz_zz_yzzz_0[j];

            tdx_xzz_zzzz_0[j] = pa_x[j] * tdx_zz_zzzz_0[j] + 0.5 * fl1_fx * ts_zz_zzzz_0[j];

            tdy_xzz_zzzz_0[j] = pa_x[j] * tdy_zz_zzzz_0[j];

            tdz_xzz_zzzz_0[j] = pa_x[j] * tdz_zz_zzzz_0[j];

            tdx_yyy_xxxx_0[j] = pa_y[j] * tdx_yy_xxxx_0[j] + fl1_fx * tdx_y_xxxx_0[j];

            tdy_yyy_xxxx_0[j] = pa_y[j] * tdy_yy_xxxx_0[j] + fl1_fx * tdy_y_xxxx_0[j] + 0.5 * fl1_fx * ts_yy_xxxx_0[j];

            tdz_yyy_xxxx_0[j] = pa_y[j] * tdz_yy_xxxx_0[j] + fl1_fx * tdz_y_xxxx_0[j];

            tdx_yyy_xxxy_0[j] = pa_y[j] * tdx_yy_xxxy_0[j] + fl1_fx * tdx_y_xxxy_0[j] + 0.5 * fl1_fx * tdx_yy_xxx_0[j];

            tdy_yyy_xxxy_0[j] =
                pa_y[j] * tdy_yy_xxxy_0[j] + fl1_fx * tdy_y_xxxy_0[j] + 0.5 * fl1_fx * tdy_yy_xxx_0[j] + 0.5 * fl1_fx * ts_yy_xxxy_0[j];

            tdz_yyy_xxxy_0[j] = pa_y[j] * tdz_yy_xxxy_0[j] + fl1_fx * tdz_y_xxxy_0[j] + 0.5 * fl1_fx * tdz_yy_xxx_0[j];

            tdx_yyy_xxxz_0[j] = pa_y[j] * tdx_yy_xxxz_0[j] + fl1_fx * tdx_y_xxxz_0[j];

            tdy_yyy_xxxz_0[j] = pa_y[j] * tdy_yy_xxxz_0[j] + fl1_fx * tdy_y_xxxz_0[j] + 0.5 * fl1_fx * ts_yy_xxxz_0[j];

            tdz_yyy_xxxz_0[j] = pa_y[j] * tdz_yy_xxxz_0[j] + fl1_fx * tdz_y_xxxz_0[j];

            tdx_yyy_xxyy_0[j] = pa_y[j] * tdx_yy_xxyy_0[j] + fl1_fx * tdx_y_xxyy_0[j] + fl1_fx * tdx_yy_xxy_0[j];

            tdy_yyy_xxyy_0[j] = pa_y[j] * tdy_yy_xxyy_0[j] + fl1_fx * tdy_y_xxyy_0[j] + fl1_fx * tdy_yy_xxy_0[j] + 0.5 * fl1_fx * ts_yy_xxyy_0[j];

            tdz_yyy_xxyy_0[j] = pa_y[j] * tdz_yy_xxyy_0[j] + fl1_fx * tdz_y_xxyy_0[j] + fl1_fx * tdz_yy_xxy_0[j];

            tdx_yyy_xxyz_0[j] = pa_y[j] * tdx_yy_xxyz_0[j] + fl1_fx * tdx_y_xxyz_0[j] + 0.5 * fl1_fx * tdx_yy_xxz_0[j];

            tdy_yyy_xxyz_0[j] =
                pa_y[j] * tdy_yy_xxyz_0[j] + fl1_fx * tdy_y_xxyz_0[j] + 0.5 * fl1_fx * tdy_yy_xxz_0[j] + 0.5 * fl1_fx * ts_yy_xxyz_0[j];

            tdz_yyy_xxyz_0[j] = pa_y[j] * tdz_yy_xxyz_0[j] + fl1_fx * tdz_y_xxyz_0[j] + 0.5 * fl1_fx * tdz_yy_xxz_0[j];

            tdx_yyy_xxzz_0[j] = pa_y[j] * tdx_yy_xxzz_0[j] + fl1_fx * tdx_y_xxzz_0[j];

            tdy_yyy_xxzz_0[j] = pa_y[j] * tdy_yy_xxzz_0[j] + fl1_fx * tdy_y_xxzz_0[j] + 0.5 * fl1_fx * ts_yy_xxzz_0[j];

            tdz_yyy_xxzz_0[j] = pa_y[j] * tdz_yy_xxzz_0[j] + fl1_fx * tdz_y_xxzz_0[j];

            tdx_yyy_xyyy_0[j] = pa_y[j] * tdx_yy_xyyy_0[j] + fl1_fx * tdx_y_xyyy_0[j] + 1.5 * fl1_fx * tdx_yy_xyy_0[j];

            tdy_yyy_xyyy_0[j] =
                pa_y[j] * tdy_yy_xyyy_0[j] + fl1_fx * tdy_y_xyyy_0[j] + 1.5 * fl1_fx * tdy_yy_xyy_0[j] + 0.5 * fl1_fx * ts_yy_xyyy_0[j];

            tdz_yyy_xyyy_0[j] = pa_y[j] * tdz_yy_xyyy_0[j] + fl1_fx * tdz_y_xyyy_0[j] + 1.5 * fl1_fx * tdz_yy_xyy_0[j];

            tdx_yyy_xyyz_0[j] = pa_y[j] * tdx_yy_xyyz_0[j] + fl1_fx * tdx_y_xyyz_0[j] + fl1_fx * tdx_yy_xyz_0[j];

            tdy_yyy_xyyz_0[j] = pa_y[j] * tdy_yy_xyyz_0[j] + fl1_fx * tdy_y_xyyz_0[j] + fl1_fx * tdy_yy_xyz_0[j] + 0.5 * fl1_fx * ts_yy_xyyz_0[j];

            tdz_yyy_xyyz_0[j] = pa_y[j] * tdz_yy_xyyz_0[j] + fl1_fx * tdz_y_xyyz_0[j] + fl1_fx * tdz_yy_xyz_0[j];

            tdx_yyy_xyzz_0[j] = pa_y[j] * tdx_yy_xyzz_0[j] + fl1_fx * tdx_y_xyzz_0[j] + 0.5 * fl1_fx * tdx_yy_xzz_0[j];

            tdy_yyy_xyzz_0[j] =
                pa_y[j] * tdy_yy_xyzz_0[j] + fl1_fx * tdy_y_xyzz_0[j] + 0.5 * fl1_fx * tdy_yy_xzz_0[j] + 0.5 * fl1_fx * ts_yy_xyzz_0[j];

            tdz_yyy_xyzz_0[j] = pa_y[j] * tdz_yy_xyzz_0[j] + fl1_fx * tdz_y_xyzz_0[j] + 0.5 * fl1_fx * tdz_yy_xzz_0[j];

            tdx_yyy_xzzz_0[j] = pa_y[j] * tdx_yy_xzzz_0[j] + fl1_fx * tdx_y_xzzz_0[j];

            tdy_yyy_xzzz_0[j] = pa_y[j] * tdy_yy_xzzz_0[j] + fl1_fx * tdy_y_xzzz_0[j] + 0.5 * fl1_fx * ts_yy_xzzz_0[j];

            tdz_yyy_xzzz_0[j] = pa_y[j] * tdz_yy_xzzz_0[j] + fl1_fx * tdz_y_xzzz_0[j];
        }

        idx++;
    }
}

void
compElectricDipoleForFG_300_350(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const int32_t              nOSFactors,
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

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_d_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_y = paDistances.data(3 * idx + 1);

        // set up pointers to auxilary integrals

        auto tdx_yy_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 55);

        auto tdy_yy_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 55);

        auto tdz_yy_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 55);

        auto tdx_yy_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 56);

        auto tdy_yy_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 56);

        auto tdz_yy_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 56);

        auto tdx_yy_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 57);

        auto tdy_yy_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 57);

        auto tdz_yy_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 57);

        auto tdx_yy_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 58);

        auto tdy_yy_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 58);

        auto tdz_yy_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 58);

        auto tdx_yy_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 59);

        auto tdy_yy_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 59);

        auto tdz_yy_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 59);

        auto tdx_yz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 60);

        auto tdy_yz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 60);

        auto tdz_yz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 60);

        auto tdx_yz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 61);

        auto tdy_yz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 61);

        auto tdz_yz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 61);

        auto tdx_yz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 62);

        auto tdy_yz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 62);

        auto tdz_yz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 62);

        auto tdx_yz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 63);

        auto tdy_yz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 63);

        auto tdz_yz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 63);

        auto tdx_yz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 64);

        auto tdy_yz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 64);

        auto tdz_yz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 64);

        auto tdx_yz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 65);

        auto tdy_yz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 65);

        auto tdz_yz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 65);

        auto tdx_yz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 66);

        auto tdy_yz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 66);

        auto tdz_yz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 66);

        auto tdx_yz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 67);

        auto tdy_yz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 67);

        auto tdz_yz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 67);

        auto tdx_yz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 68);

        auto tdy_yz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 68);

        auto tdz_yz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 68);

        auto tdx_yz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 69);

        auto tdy_yz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 69);

        auto tdz_yz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 69);

        auto tdx_yz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 70);

        auto tdy_yz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 70);

        auto tdz_yz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 70);

        auto tdx_yz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 71);

        auto tdy_yz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 71);

        auto tdx_y_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 25);

        auto tdy_y_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 25);

        auto tdz_y_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 25);

        auto tdx_y_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 26);

        auto tdy_y_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 26);

        auto tdz_y_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 26);

        auto tdx_y_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 27);

        auto tdy_y_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 27);

        auto tdz_y_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 27);

        auto tdx_y_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 28);

        auto tdy_y_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 28);

        auto tdz_y_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 28);

        auto tdx_y_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 29);

        auto tdy_y_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 29);

        auto tdz_y_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 29);

        auto tdx_z_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 30);

        auto tdy_z_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 30);

        auto tdz_z_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 30);

        auto tdx_z_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 31);

        auto tdy_z_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 31);

        auto tdz_z_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 31);

        auto tdx_z_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 32);

        auto tdy_z_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 32);

        auto tdz_z_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 32);

        auto tdx_z_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 33);

        auto tdy_z_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 33);

        auto tdz_z_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 33);

        auto tdx_z_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 34);

        auto tdy_z_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 34);

        auto tdz_z_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 34);

        auto tdx_z_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 35);

        auto tdy_z_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 35);

        auto tdz_z_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 35);

        auto tdx_z_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 36);

        auto tdy_z_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 36);

        auto tdz_z_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 36);

        auto tdx_z_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 37);

        auto tdy_z_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 37);

        auto tdz_z_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 37);

        auto tdx_z_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 38);

        auto tdy_z_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 38);

        auto tdz_z_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 38);

        auto tdx_z_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 39);

        auto tdy_z_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 39);

        auto tdz_z_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 39);

        auto tdx_z_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 40);

        auto tdy_z_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 40);

        auto tdz_z_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 40);

        auto tdx_z_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 41);

        auto tdy_z_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 41);

        auto tdx_yy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 36);

        auto tdy_yy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 36);

        auto tdz_yy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 36);

        auto tdx_yy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 37);

        auto tdy_yy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 37);

        auto tdz_yy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 37);

        auto tdx_yy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 38);

        auto tdy_yy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 38);

        auto tdz_yy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 38);

        auto tdx_yy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 39);

        auto tdy_yy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 39);

        auto tdz_yy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 39);

        auto tdx_yz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 40);

        auto tdy_yz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 40);

        auto tdz_yz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 40);

        auto tdx_yz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 41);

        auto tdy_yz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 41);

        auto tdz_yz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 41);

        auto tdx_yz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 42);

        auto tdy_yz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 42);

        auto tdz_yz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 42);

        auto tdx_yz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 43);

        auto tdy_yz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 43);

        auto tdz_yz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 43);

        auto tdx_yz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 44);

        auto tdy_yz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 44);

        auto tdz_yz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 44);

        auto tdx_yz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 45);

        auto tdy_yz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 45);

        auto tdz_yz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 45);

        auto tdx_yz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 46);

        auto tdy_yz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 46);

        auto tdz_yz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 46);

        auto tdx_yz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 47);

        auto tdy_yz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 47);

        auto ts_yy_yyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 55);

        auto ts_yy_yyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 56);

        auto ts_yy_yyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 57);

        auto ts_yy_yzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 58);

        auto ts_yy_zzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 59);

        auto ts_yz_xxxx_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 60);

        auto ts_yz_xxxy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 61);

        auto ts_yz_xxxz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 62);

        auto ts_yz_xxyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 63);

        auto ts_yz_xxyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 64);

        auto ts_yz_xxzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 65);

        auto ts_yz_xyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 66);

        auto ts_yz_xyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 67);

        auto ts_yz_xyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 68);

        auto ts_yz_xzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 69);

        auto ts_yz_yyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 70);

        auto ts_yz_yyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 71);

        // set up pointers to integrals

        auto tdx_yyy_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 100);

        auto tdy_yyy_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 100);

        auto tdz_yyy_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 100);

        auto tdx_yyy_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 101);

        auto tdy_yyy_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 101);

        auto tdz_yyy_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 101);

        auto tdx_yyy_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 102);

        auto tdy_yyy_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 102);

        auto tdz_yyy_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 102);

        auto tdx_yyy_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 103);

        auto tdy_yyy_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 103);

        auto tdz_yyy_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 103);

        auto tdx_yyy_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 104);

        auto tdy_yyy_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 104);

        auto tdz_yyy_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 104);

        auto tdx_yyz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 105);

        auto tdy_yyz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 105);

        auto tdz_yyz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 105);

        auto tdx_yyz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 106);

        auto tdy_yyz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 106);

        auto tdz_yyz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 106);

        auto tdx_yyz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 107);

        auto tdy_yyz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 107);

        auto tdz_yyz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 107);

        auto tdx_yyz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 108);

        auto tdy_yyz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 108);

        auto tdz_yyz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 108);

        auto tdx_yyz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 109);

        auto tdy_yyz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 109);

        auto tdz_yyz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 109);

        auto tdx_yyz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 110);

        auto tdy_yyz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 110);

        auto tdz_yyz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 110);

        auto tdx_yyz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 111);

        auto tdy_yyz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 111);

        auto tdz_yyz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 111);

        auto tdx_yyz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 112);

        auto tdy_yyz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 112);

        auto tdz_yyz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 112);

        auto tdx_yyz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 113);

        auto tdy_yyz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 113);

        auto tdz_yyz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 113);

        auto tdx_yyz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 114);

        auto tdy_yyz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 114);

        auto tdz_yyz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 114);

        auto tdx_yyz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 115);

        auto tdy_yyz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 115);

        auto tdz_yyz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 115);

        auto tdx_yyz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 116);

        auto tdy_yyz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 116);

        // Batch of Integrals (300,350)

        #pragma omp simd aligned(fx, pa_y, tdx_y_yyyy_0, tdx_y_yyyz_0, tdx_y_yyzz_0, tdx_y_yzzz_0, \
                                     tdx_y_zzzz_0, tdx_yy_yyy_0, tdx_yy_yyyy_0, tdx_yy_yyyz_0, tdx_yy_yyz_0, \
                                     tdx_yy_yyzz_0, tdx_yy_yzz_0, tdx_yy_yzzz_0, tdx_yy_zzz_0, tdx_yy_zzzz_0, \
                                     tdx_yyy_yyyy_0, tdx_yyy_yyyz_0, tdx_yyy_yyzz_0, tdx_yyy_yzzz_0, tdx_yyy_zzzz_0, \
                                     tdx_yyz_xxxx_0, tdx_yyz_xxxy_0, tdx_yyz_xxxz_0, tdx_yyz_xxyy_0, tdx_yyz_xxyz_0, \
                                     tdx_yyz_xxzz_0, tdx_yyz_xyyy_0, tdx_yyz_xyyz_0, tdx_yyz_xyzz_0, tdx_yyz_xzzz_0, \
                                     tdx_yyz_yyyy_0, tdx_yyz_yyyz_0, tdx_yz_xxx_0, tdx_yz_xxxx_0, tdx_yz_xxxy_0, \
                                     tdx_yz_xxxz_0, tdx_yz_xxy_0, tdx_yz_xxyy_0, tdx_yz_xxyz_0, tdx_yz_xxz_0, \
                                     tdx_yz_xxzz_0, tdx_yz_xyy_0, tdx_yz_xyyy_0, tdx_yz_xyyz_0, tdx_yz_xyz_0, \
                                     tdx_yz_xyzz_0, tdx_yz_xzz_0, tdx_yz_xzzz_0, tdx_yz_yyy_0, tdx_yz_yyyy_0, \
                                     tdx_yz_yyyz_0, tdx_yz_yyz_0, tdx_z_xxxx_0, tdx_z_xxxy_0, tdx_z_xxxz_0, tdx_z_xxyy_0, \
                                     tdx_z_xxyz_0, tdx_z_xxzz_0, tdx_z_xyyy_0, tdx_z_xyyz_0, tdx_z_xyzz_0, tdx_z_xzzz_0, \
                                     tdx_z_yyyy_0, tdx_z_yyyz_0, tdy_y_yyyy_0, tdy_y_yyyz_0, tdy_y_yyzz_0, tdy_y_yzzz_0, \
                                     tdy_y_zzzz_0, tdy_yy_yyy_0, tdy_yy_yyyy_0, tdy_yy_yyyz_0, tdy_yy_yyz_0, \
                                     tdy_yy_yyzz_0, tdy_yy_yzz_0, tdy_yy_yzzz_0, tdy_yy_zzz_0, tdy_yy_zzzz_0, \
                                     tdy_yyy_yyyy_0, tdy_yyy_yyyz_0, tdy_yyy_yyzz_0, tdy_yyy_yzzz_0, tdy_yyy_zzzz_0, \
                                     tdy_yyz_xxxx_0, tdy_yyz_xxxy_0, tdy_yyz_xxxz_0, tdy_yyz_xxyy_0, tdy_yyz_xxyz_0, \
                                     tdy_yyz_xxzz_0, tdy_yyz_xyyy_0, tdy_yyz_xyyz_0, tdy_yyz_xyzz_0, tdy_yyz_xzzz_0, \
                                     tdy_yyz_yyyy_0, tdy_yyz_yyyz_0, tdy_yz_xxx_0, tdy_yz_xxxx_0, tdy_yz_xxxy_0, \
                                     tdy_yz_xxxz_0, tdy_yz_xxy_0, tdy_yz_xxyy_0, tdy_yz_xxyz_0, tdy_yz_xxz_0, \
                                     tdy_yz_xxzz_0, tdy_yz_xyy_0, tdy_yz_xyyy_0, tdy_yz_xyyz_0, tdy_yz_xyz_0, \
                                     tdy_yz_xyzz_0, tdy_yz_xzz_0, tdy_yz_xzzz_0, tdy_yz_yyy_0, tdy_yz_yyyy_0, \
                                     tdy_yz_yyyz_0, tdy_yz_yyz_0, tdy_z_xxxx_0, tdy_z_xxxy_0, tdy_z_xxxz_0, tdy_z_xxyy_0, \
                                     tdy_z_xxyz_0, tdy_z_xxzz_0, tdy_z_xyyy_0, tdy_z_xyyz_0, tdy_z_xyzz_0, tdy_z_xzzz_0, \
                                     tdy_z_yyyy_0, tdy_z_yyyz_0, tdz_y_yyyy_0, tdz_y_yyyz_0, tdz_y_yyzz_0, tdz_y_yzzz_0, \
                                     tdz_y_zzzz_0, tdz_yy_yyy_0, tdz_yy_yyyy_0, tdz_yy_yyyz_0, tdz_yy_yyz_0, \
                                     tdz_yy_yyzz_0, tdz_yy_yzz_0, tdz_yy_yzzz_0, tdz_yy_zzz_0, tdz_yy_zzzz_0, \
                                     tdz_yyy_yyyy_0, tdz_yyy_yyyz_0, tdz_yyy_yyzz_0, tdz_yyy_yzzz_0, tdz_yyy_zzzz_0, \
                                     tdz_yyz_xxxx_0, tdz_yyz_xxxy_0, tdz_yyz_xxxz_0, tdz_yyz_xxyy_0, tdz_yyz_xxyz_0, \
                                     tdz_yyz_xxzz_0, tdz_yyz_xyyy_0, tdz_yyz_xyyz_0, tdz_yyz_xyzz_0, tdz_yyz_xzzz_0, \
                                     tdz_yyz_yyyy_0, tdz_yz_xxx_0, tdz_yz_xxxx_0, tdz_yz_xxxy_0, tdz_yz_xxxz_0, \
                                     tdz_yz_xxy_0, tdz_yz_xxyy_0, tdz_yz_xxyz_0, tdz_yz_xxz_0, tdz_yz_xxzz_0, \
                                     tdz_yz_xyy_0, tdz_yz_xyyy_0, tdz_yz_xyyz_0, tdz_yz_xyz_0, tdz_yz_xyzz_0, \
                                     tdz_yz_xzz_0, tdz_yz_xzzz_0, tdz_yz_yyy_0, tdz_yz_yyyy_0, tdz_z_xxxx_0, \
                                     tdz_z_xxxy_0, tdz_z_xxxz_0, tdz_z_xxyy_0, tdz_z_xxyz_0, tdz_z_xxzz_0, tdz_z_xyyy_0, \
                                     tdz_z_xyyz_0, tdz_z_xyzz_0, tdz_z_xzzz_0, tdz_z_yyyy_0, ts_yy_yyyy_0, ts_yy_yyyz_0, \
                                     ts_yy_yyzz_0, ts_yy_yzzz_0, ts_yy_zzzz_0, ts_yz_xxxx_0, ts_yz_xxxy_0, ts_yz_xxxz_0, \
                                     ts_yz_xxyy_0, ts_yz_xxyz_0, ts_yz_xxzz_0, ts_yz_xyyy_0, ts_yz_xyyz_0, ts_yz_xyzz_0, \
                                     ts_yz_xzzz_0, ts_yz_yyyy_0, ts_yz_yyyz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            tdx_yyy_yyyy_0[j] = pa_y[j] * tdx_yy_yyyy_0[j] + fl1_fx * tdx_y_yyyy_0[j] + 2.0 * fl1_fx * tdx_yy_yyy_0[j];

            tdy_yyy_yyyy_0[j] =
                pa_y[j] * tdy_yy_yyyy_0[j] + fl1_fx * tdy_y_yyyy_0[j] + 2.0 * fl1_fx * tdy_yy_yyy_0[j] + 0.5 * fl1_fx * ts_yy_yyyy_0[j];

            tdz_yyy_yyyy_0[j] = pa_y[j] * tdz_yy_yyyy_0[j] + fl1_fx * tdz_y_yyyy_0[j] + 2.0 * fl1_fx * tdz_yy_yyy_0[j];

            tdx_yyy_yyyz_0[j] = pa_y[j] * tdx_yy_yyyz_0[j] + fl1_fx * tdx_y_yyyz_0[j] + 1.5 * fl1_fx * tdx_yy_yyz_0[j];

            tdy_yyy_yyyz_0[j] =
                pa_y[j] * tdy_yy_yyyz_0[j] + fl1_fx * tdy_y_yyyz_0[j] + 1.5 * fl1_fx * tdy_yy_yyz_0[j] + 0.5 * fl1_fx * ts_yy_yyyz_0[j];

            tdz_yyy_yyyz_0[j] = pa_y[j] * tdz_yy_yyyz_0[j] + fl1_fx * tdz_y_yyyz_0[j] + 1.5 * fl1_fx * tdz_yy_yyz_0[j];

            tdx_yyy_yyzz_0[j] = pa_y[j] * tdx_yy_yyzz_0[j] + fl1_fx * tdx_y_yyzz_0[j] + fl1_fx * tdx_yy_yzz_0[j];

            tdy_yyy_yyzz_0[j] = pa_y[j] * tdy_yy_yyzz_0[j] + fl1_fx * tdy_y_yyzz_0[j] + fl1_fx * tdy_yy_yzz_0[j] + 0.5 * fl1_fx * ts_yy_yyzz_0[j];

            tdz_yyy_yyzz_0[j] = pa_y[j] * tdz_yy_yyzz_0[j] + fl1_fx * tdz_y_yyzz_0[j] + fl1_fx * tdz_yy_yzz_0[j];

            tdx_yyy_yzzz_0[j] = pa_y[j] * tdx_yy_yzzz_0[j] + fl1_fx * tdx_y_yzzz_0[j] + 0.5 * fl1_fx * tdx_yy_zzz_0[j];

            tdy_yyy_yzzz_0[j] =
                pa_y[j] * tdy_yy_yzzz_0[j] + fl1_fx * tdy_y_yzzz_0[j] + 0.5 * fl1_fx * tdy_yy_zzz_0[j] + 0.5 * fl1_fx * ts_yy_yzzz_0[j];

            tdz_yyy_yzzz_0[j] = pa_y[j] * tdz_yy_yzzz_0[j] + fl1_fx * tdz_y_yzzz_0[j] + 0.5 * fl1_fx * tdz_yy_zzz_0[j];

            tdx_yyy_zzzz_0[j] = pa_y[j] * tdx_yy_zzzz_0[j] + fl1_fx * tdx_y_zzzz_0[j];

            tdy_yyy_zzzz_0[j] = pa_y[j] * tdy_yy_zzzz_0[j] + fl1_fx * tdy_y_zzzz_0[j] + 0.5 * fl1_fx * ts_yy_zzzz_0[j];

            tdz_yyy_zzzz_0[j] = pa_y[j] * tdz_yy_zzzz_0[j] + fl1_fx * tdz_y_zzzz_0[j];

            tdx_yyz_xxxx_0[j] = pa_y[j] * tdx_yz_xxxx_0[j] + 0.5 * fl1_fx * tdx_z_xxxx_0[j];

            tdy_yyz_xxxx_0[j] = pa_y[j] * tdy_yz_xxxx_0[j] + 0.5 * fl1_fx * tdy_z_xxxx_0[j] + 0.5 * fl1_fx * ts_yz_xxxx_0[j];

            tdz_yyz_xxxx_0[j] = pa_y[j] * tdz_yz_xxxx_0[j] + 0.5 * fl1_fx * tdz_z_xxxx_0[j];

            tdx_yyz_xxxy_0[j] = pa_y[j] * tdx_yz_xxxy_0[j] + 0.5 * fl1_fx * tdx_z_xxxy_0[j] + 0.5 * fl1_fx * tdx_yz_xxx_0[j];

            tdy_yyz_xxxy_0[j] =
                pa_y[j] * tdy_yz_xxxy_0[j] + 0.5 * fl1_fx * tdy_z_xxxy_0[j] + 0.5 * fl1_fx * tdy_yz_xxx_0[j] + 0.5 * fl1_fx * ts_yz_xxxy_0[j];

            tdz_yyz_xxxy_0[j] = pa_y[j] * tdz_yz_xxxy_0[j] + 0.5 * fl1_fx * tdz_z_xxxy_0[j] + 0.5 * fl1_fx * tdz_yz_xxx_0[j];

            tdx_yyz_xxxz_0[j] = pa_y[j] * tdx_yz_xxxz_0[j] + 0.5 * fl1_fx * tdx_z_xxxz_0[j];

            tdy_yyz_xxxz_0[j] = pa_y[j] * tdy_yz_xxxz_0[j] + 0.5 * fl1_fx * tdy_z_xxxz_0[j] + 0.5 * fl1_fx * ts_yz_xxxz_0[j];

            tdz_yyz_xxxz_0[j] = pa_y[j] * tdz_yz_xxxz_0[j] + 0.5 * fl1_fx * tdz_z_xxxz_0[j];

            tdx_yyz_xxyy_0[j] = pa_y[j] * tdx_yz_xxyy_0[j] + 0.5 * fl1_fx * tdx_z_xxyy_0[j] + fl1_fx * tdx_yz_xxy_0[j];

            tdy_yyz_xxyy_0[j] =
                pa_y[j] * tdy_yz_xxyy_0[j] + 0.5 * fl1_fx * tdy_z_xxyy_0[j] + fl1_fx * tdy_yz_xxy_0[j] + 0.5 * fl1_fx * ts_yz_xxyy_0[j];

            tdz_yyz_xxyy_0[j] = pa_y[j] * tdz_yz_xxyy_0[j] + 0.5 * fl1_fx * tdz_z_xxyy_0[j] + fl1_fx * tdz_yz_xxy_0[j];

            tdx_yyz_xxyz_0[j] = pa_y[j] * tdx_yz_xxyz_0[j] + 0.5 * fl1_fx * tdx_z_xxyz_0[j] + 0.5 * fl1_fx * tdx_yz_xxz_0[j];

            tdy_yyz_xxyz_0[j] =
                pa_y[j] * tdy_yz_xxyz_0[j] + 0.5 * fl1_fx * tdy_z_xxyz_0[j] + 0.5 * fl1_fx * tdy_yz_xxz_0[j] + 0.5 * fl1_fx * ts_yz_xxyz_0[j];

            tdz_yyz_xxyz_0[j] = pa_y[j] * tdz_yz_xxyz_0[j] + 0.5 * fl1_fx * tdz_z_xxyz_0[j] + 0.5 * fl1_fx * tdz_yz_xxz_0[j];

            tdx_yyz_xxzz_0[j] = pa_y[j] * tdx_yz_xxzz_0[j] + 0.5 * fl1_fx * tdx_z_xxzz_0[j];

            tdy_yyz_xxzz_0[j] = pa_y[j] * tdy_yz_xxzz_0[j] + 0.5 * fl1_fx * tdy_z_xxzz_0[j] + 0.5 * fl1_fx * ts_yz_xxzz_0[j];

            tdz_yyz_xxzz_0[j] = pa_y[j] * tdz_yz_xxzz_0[j] + 0.5 * fl1_fx * tdz_z_xxzz_0[j];

            tdx_yyz_xyyy_0[j] = pa_y[j] * tdx_yz_xyyy_0[j] + 0.5 * fl1_fx * tdx_z_xyyy_0[j] + 1.5 * fl1_fx * tdx_yz_xyy_0[j];

            tdy_yyz_xyyy_0[j] =
                pa_y[j] * tdy_yz_xyyy_0[j] + 0.5 * fl1_fx * tdy_z_xyyy_0[j] + 1.5 * fl1_fx * tdy_yz_xyy_0[j] + 0.5 * fl1_fx * ts_yz_xyyy_0[j];

            tdz_yyz_xyyy_0[j] = pa_y[j] * tdz_yz_xyyy_0[j] + 0.5 * fl1_fx * tdz_z_xyyy_0[j] + 1.5 * fl1_fx * tdz_yz_xyy_0[j];

            tdx_yyz_xyyz_0[j] = pa_y[j] * tdx_yz_xyyz_0[j] + 0.5 * fl1_fx * tdx_z_xyyz_0[j] + fl1_fx * tdx_yz_xyz_0[j];

            tdy_yyz_xyyz_0[j] =
                pa_y[j] * tdy_yz_xyyz_0[j] + 0.5 * fl1_fx * tdy_z_xyyz_0[j] + fl1_fx * tdy_yz_xyz_0[j] + 0.5 * fl1_fx * ts_yz_xyyz_0[j];

            tdz_yyz_xyyz_0[j] = pa_y[j] * tdz_yz_xyyz_0[j] + 0.5 * fl1_fx * tdz_z_xyyz_0[j] + fl1_fx * tdz_yz_xyz_0[j];

            tdx_yyz_xyzz_0[j] = pa_y[j] * tdx_yz_xyzz_0[j] + 0.5 * fl1_fx * tdx_z_xyzz_0[j] + 0.5 * fl1_fx * tdx_yz_xzz_0[j];

            tdy_yyz_xyzz_0[j] =
                pa_y[j] * tdy_yz_xyzz_0[j] + 0.5 * fl1_fx * tdy_z_xyzz_0[j] + 0.5 * fl1_fx * tdy_yz_xzz_0[j] + 0.5 * fl1_fx * ts_yz_xyzz_0[j];

            tdz_yyz_xyzz_0[j] = pa_y[j] * tdz_yz_xyzz_0[j] + 0.5 * fl1_fx * tdz_z_xyzz_0[j] + 0.5 * fl1_fx * tdz_yz_xzz_0[j];

            tdx_yyz_xzzz_0[j] = pa_y[j] * tdx_yz_xzzz_0[j] + 0.5 * fl1_fx * tdx_z_xzzz_0[j];

            tdy_yyz_xzzz_0[j] = pa_y[j] * tdy_yz_xzzz_0[j] + 0.5 * fl1_fx * tdy_z_xzzz_0[j] + 0.5 * fl1_fx * ts_yz_xzzz_0[j];

            tdz_yyz_xzzz_0[j] = pa_y[j] * tdz_yz_xzzz_0[j] + 0.5 * fl1_fx * tdz_z_xzzz_0[j];

            tdx_yyz_yyyy_0[j] = pa_y[j] * tdx_yz_yyyy_0[j] + 0.5 * fl1_fx * tdx_z_yyyy_0[j] + 2.0 * fl1_fx * tdx_yz_yyy_0[j];

            tdy_yyz_yyyy_0[j] =
                pa_y[j] * tdy_yz_yyyy_0[j] + 0.5 * fl1_fx * tdy_z_yyyy_0[j] + 2.0 * fl1_fx * tdy_yz_yyy_0[j] + 0.5 * fl1_fx * ts_yz_yyyy_0[j];

            tdz_yyz_yyyy_0[j] = pa_y[j] * tdz_yz_yyyy_0[j] + 0.5 * fl1_fx * tdz_z_yyyy_0[j] + 2.0 * fl1_fx * tdz_yz_yyy_0[j];

            tdx_yyz_yyyz_0[j] = pa_y[j] * tdx_yz_yyyz_0[j] + 0.5 * fl1_fx * tdx_z_yyyz_0[j] + 1.5 * fl1_fx * tdx_yz_yyz_0[j];

            tdy_yyz_yyyz_0[j] =
                pa_y[j] * tdy_yz_yyyz_0[j] + 0.5 * fl1_fx * tdy_z_yyyz_0[j] + 1.5 * fl1_fx * tdy_yz_yyz_0[j] + 0.5 * fl1_fx * ts_yz_yyyz_0[j];
        }

        idx++;
    }
}

void
compElectricDipoleForFG_350_400(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const int32_t              nOSFactors,
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

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_d_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_y = paDistances.data(3 * idx + 1);

        // set up pointers to auxilary integrals

        auto tdz_yz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 71);

        auto tdx_yz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 72);

        auto tdy_yz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 72);

        auto tdz_yz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 72);

        auto tdx_yz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 73);

        auto tdy_yz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 73);

        auto tdz_yz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 73);

        auto tdx_yz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 74);

        auto tdy_yz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 74);

        auto tdz_yz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 74);

        auto tdx_zz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 75);

        auto tdy_zz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 75);

        auto tdz_zz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 75);

        auto tdx_zz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 76);

        auto tdy_zz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 76);

        auto tdz_zz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 76);

        auto tdx_zz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 77);

        auto tdy_zz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 77);

        auto tdz_zz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 77);

        auto tdx_zz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 78);

        auto tdy_zz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 78);

        auto tdz_zz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 78);

        auto tdx_zz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 79);

        auto tdy_zz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 79);

        auto tdz_zz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 79);

        auto tdx_zz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 80);

        auto tdy_zz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 80);

        auto tdz_zz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 80);

        auto tdx_zz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 81);

        auto tdy_zz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 81);

        auto tdz_zz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 81);

        auto tdx_zz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 82);

        auto tdy_zz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 82);

        auto tdz_zz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 82);

        auto tdx_zz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 83);

        auto tdy_zz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 83);

        auto tdz_zz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 83);

        auto tdx_zz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 84);

        auto tdy_zz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 84);

        auto tdz_zz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 84);

        auto tdx_zz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 85);

        auto tdy_zz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 85);

        auto tdz_zz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 85);

        auto tdx_zz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 86);

        auto tdy_zz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 86);

        auto tdz_zz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 86);

        auto tdx_zz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 87);

        auto tdy_zz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 87);

        auto tdz_zz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 87);

        auto tdx_zz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 88);

        auto tdz_z_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 41);

        auto tdx_z_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 42);

        auto tdy_z_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 42);

        auto tdz_z_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 42);

        auto tdx_z_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 43);

        auto tdy_z_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 43);

        auto tdz_z_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 43);

        auto tdx_z_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 44);

        auto tdy_z_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 44);

        auto tdz_z_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 44);

        auto tdz_yz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 47);

        auto tdx_yz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 48);

        auto tdy_yz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 48);

        auto tdz_yz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 48);

        auto tdx_yz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 49);

        auto tdy_yz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 49);

        auto tdz_yz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 49);

        auto tdx_zz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 50);

        auto tdy_zz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 50);

        auto tdz_zz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 50);

        auto tdx_zz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 51);

        auto tdy_zz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 51);

        auto tdz_zz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 51);

        auto tdx_zz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 52);

        auto tdy_zz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 52);

        auto tdz_zz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 52);

        auto tdx_zz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 53);

        auto tdy_zz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 53);

        auto tdz_zz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 53);

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

        auto ts_yz_yyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 72);

        auto ts_yz_yzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 73);

        auto ts_yz_zzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 74);

        auto ts_zz_xxxx_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 75);

        auto ts_zz_xxxy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 76);

        auto ts_zz_xxxz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 77);

        auto ts_zz_xxyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 78);

        auto ts_zz_xxyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 79);

        auto ts_zz_xxzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 80);

        auto ts_zz_xyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 81);

        auto ts_zz_xyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 82);

        auto ts_zz_xyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 83);

        auto ts_zz_xzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 84);

        auto ts_zz_yyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 85);

        auto ts_zz_yyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 86);

        auto ts_zz_yyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 87);

        // set up pointers to integrals

        auto tdz_yyz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 116);

        auto tdx_yyz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 117);

        auto tdy_yyz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 117);

        auto tdz_yyz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 117);

        auto tdx_yyz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 118);

        auto tdy_yyz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 118);

        auto tdz_yyz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 118);

        auto tdx_yyz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 119);

        auto tdy_yyz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 119);

        auto tdz_yyz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 119);

        auto tdx_yzz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 120);

        auto tdy_yzz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 120);

        auto tdz_yzz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 120);

        auto tdx_yzz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 121);

        auto tdy_yzz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 121);

        auto tdz_yzz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 121);

        auto tdx_yzz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 122);

        auto tdy_yzz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 122);

        auto tdz_yzz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 122);

        auto tdx_yzz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 123);

        auto tdy_yzz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 123);

        auto tdz_yzz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 123);

        auto tdx_yzz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 124);

        auto tdy_yzz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 124);

        auto tdz_yzz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 124);

        auto tdx_yzz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 125);

        auto tdy_yzz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 125);

        auto tdz_yzz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 125);

        auto tdx_yzz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 126);

        auto tdy_yzz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 126);

        auto tdz_yzz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 126);

        auto tdx_yzz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 127);

        auto tdy_yzz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 127);

        auto tdz_yzz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 127);

        auto tdx_yzz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 128);

        auto tdy_yzz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 128);

        auto tdz_yzz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 128);

        auto tdx_yzz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 129);

        auto tdy_yzz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 129);

        auto tdz_yzz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 129);

        auto tdx_yzz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 130);

        auto tdy_yzz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 130);

        auto tdz_yzz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 130);

        auto tdx_yzz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 131);

        auto tdy_yzz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 131);

        auto tdz_yzz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 131);

        auto tdx_yzz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 132);

        auto tdy_yzz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 132);

        auto tdz_yzz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 132);

        auto tdx_yzz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 133);

        // Batch of Integrals (350,400)

        #pragma omp simd aligned(fx, pa_y, tdx_yyz_yyzz_0, tdx_yyz_yzzz_0, tdx_yyz_zzzz_0, \
                                     tdx_yz_yyzz_0, tdx_yz_yzz_0, tdx_yz_yzzz_0, tdx_yz_zzz_0, tdx_yz_zzzz_0, \
                                     tdx_yzz_xxxx_0, tdx_yzz_xxxy_0, tdx_yzz_xxxz_0, tdx_yzz_xxyy_0, tdx_yzz_xxyz_0, \
                                     tdx_yzz_xxzz_0, tdx_yzz_xyyy_0, tdx_yzz_xyyz_0, tdx_yzz_xyzz_0, tdx_yzz_xzzz_0, \
                                     tdx_yzz_yyyy_0, tdx_yzz_yyyz_0, tdx_yzz_yyzz_0, tdx_yzz_yzzz_0, tdx_z_yyzz_0, \
                                     tdx_z_yzzz_0, tdx_z_zzzz_0, tdx_zz_xxx_0, tdx_zz_xxxx_0, tdx_zz_xxxy_0, \
                                     tdx_zz_xxxz_0, tdx_zz_xxy_0, tdx_zz_xxyy_0, tdx_zz_xxyz_0, tdx_zz_xxz_0, \
                                     tdx_zz_xxzz_0, tdx_zz_xyy_0, tdx_zz_xyyy_0, tdx_zz_xyyz_0, tdx_zz_xyz_0, \
                                     tdx_zz_xyzz_0, tdx_zz_xzz_0, tdx_zz_xzzz_0, tdx_zz_yyy_0, tdx_zz_yyyy_0, \
                                     tdx_zz_yyyz_0, tdx_zz_yyz_0, tdx_zz_yyzz_0, tdx_zz_yzz_0, tdx_zz_yzzz_0, \
                                     tdx_zz_zzz_0, tdy_yyz_yyzz_0, tdy_yyz_yzzz_0, tdy_yyz_zzzz_0, tdy_yz_yyzz_0, \
                                     tdy_yz_yzz_0, tdy_yz_yzzz_0, tdy_yz_zzz_0, tdy_yz_zzzz_0, tdy_yzz_xxxx_0, \
                                     tdy_yzz_xxxy_0, tdy_yzz_xxxz_0, tdy_yzz_xxyy_0, tdy_yzz_xxyz_0, tdy_yzz_xxzz_0, \
                                     tdy_yzz_xyyy_0, tdy_yzz_xyyz_0, tdy_yzz_xyzz_0, tdy_yzz_xzzz_0, tdy_yzz_yyyy_0, \
                                     tdy_yzz_yyyz_0, tdy_yzz_yyzz_0, tdy_z_yyzz_0, tdy_z_yzzz_0, tdy_z_zzzz_0, \
                                     tdy_zz_xxx_0, tdy_zz_xxxx_0, tdy_zz_xxxy_0, tdy_zz_xxxz_0, tdy_zz_xxy_0, \
                                     tdy_zz_xxyy_0, tdy_zz_xxyz_0, tdy_zz_xxz_0, tdy_zz_xxzz_0, tdy_zz_xyy_0, \
                                     tdy_zz_xyyy_0, tdy_zz_xyyz_0, tdy_zz_xyz_0, tdy_zz_xyzz_0, tdy_zz_xzz_0, \
                                     tdy_zz_xzzz_0, tdy_zz_yyy_0, tdy_zz_yyyy_0, tdy_zz_yyyz_0, tdy_zz_yyz_0, \
                                     tdy_zz_yyzz_0, tdy_zz_yzz_0, tdz_yyz_yyyz_0, tdz_yyz_yyzz_0, tdz_yyz_yzzz_0, \
                                     tdz_yyz_zzzz_0, tdz_yz_yyyz_0, tdz_yz_yyz_0, tdz_yz_yyzz_0, tdz_yz_yzz_0, \
                                     tdz_yz_yzzz_0, tdz_yz_zzz_0, tdz_yz_zzzz_0, tdz_yzz_xxxx_0, tdz_yzz_xxxy_0, \
                                     tdz_yzz_xxxz_0, tdz_yzz_xxyy_0, tdz_yzz_xxyz_0, tdz_yzz_xxzz_0, tdz_yzz_xyyy_0, \
                                     tdz_yzz_xyyz_0, tdz_yzz_xyzz_0, tdz_yzz_xzzz_0, tdz_yzz_yyyy_0, tdz_yzz_yyyz_0, \
                                     tdz_yzz_yyzz_0, tdz_z_yyyz_0, tdz_z_yyzz_0, tdz_z_yzzz_0, tdz_z_zzzz_0, tdz_zz_xxx_0, \
                                     tdz_zz_xxxx_0, tdz_zz_xxxy_0, tdz_zz_xxxz_0, tdz_zz_xxy_0, tdz_zz_xxyy_0, \
                                     tdz_zz_xxyz_0, tdz_zz_xxz_0, tdz_zz_xxzz_0, tdz_zz_xyy_0, tdz_zz_xyyy_0, \
                                     tdz_zz_xyyz_0, tdz_zz_xyz_0, tdz_zz_xyzz_0, tdz_zz_xzz_0, tdz_zz_xzzz_0, \
                                     tdz_zz_yyy_0, tdz_zz_yyyy_0, tdz_zz_yyyz_0, tdz_zz_yyz_0, tdz_zz_yyzz_0, \
                                     tdz_zz_yzz_0, ts_yz_yyzz_0, ts_yz_yzzz_0, ts_yz_zzzz_0, ts_zz_xxxx_0, ts_zz_xxxy_0, \
                                     ts_zz_xxxz_0, ts_zz_xxyy_0, ts_zz_xxyz_0, ts_zz_xxzz_0, ts_zz_xyyy_0, ts_zz_xyyz_0, \
                                     ts_zz_xyzz_0, ts_zz_xzzz_0, ts_zz_yyyy_0, ts_zz_yyyz_0, ts_zz_yyzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            tdz_yyz_yyyz_0[j] = pa_y[j] * tdz_yz_yyyz_0[j] + 0.5 * fl1_fx * tdz_z_yyyz_0[j] + 1.5 * fl1_fx * tdz_yz_yyz_0[j];

            tdx_yyz_yyzz_0[j] = pa_y[j] * tdx_yz_yyzz_0[j] + 0.5 * fl1_fx * tdx_z_yyzz_0[j] + fl1_fx * tdx_yz_yzz_0[j];

            tdy_yyz_yyzz_0[j] =
                pa_y[j] * tdy_yz_yyzz_0[j] + 0.5 * fl1_fx * tdy_z_yyzz_0[j] + fl1_fx * tdy_yz_yzz_0[j] + 0.5 * fl1_fx * ts_yz_yyzz_0[j];

            tdz_yyz_yyzz_0[j] = pa_y[j] * tdz_yz_yyzz_0[j] + 0.5 * fl1_fx * tdz_z_yyzz_0[j] + fl1_fx * tdz_yz_yzz_0[j];

            tdx_yyz_yzzz_0[j] = pa_y[j] * tdx_yz_yzzz_0[j] + 0.5 * fl1_fx * tdx_z_yzzz_0[j] + 0.5 * fl1_fx * tdx_yz_zzz_0[j];

            tdy_yyz_yzzz_0[j] =
                pa_y[j] * tdy_yz_yzzz_0[j] + 0.5 * fl1_fx * tdy_z_yzzz_0[j] + 0.5 * fl1_fx * tdy_yz_zzz_0[j] + 0.5 * fl1_fx * ts_yz_yzzz_0[j];

            tdz_yyz_yzzz_0[j] = pa_y[j] * tdz_yz_yzzz_0[j] + 0.5 * fl1_fx * tdz_z_yzzz_0[j] + 0.5 * fl1_fx * tdz_yz_zzz_0[j];

            tdx_yyz_zzzz_0[j] = pa_y[j] * tdx_yz_zzzz_0[j] + 0.5 * fl1_fx * tdx_z_zzzz_0[j];

            tdy_yyz_zzzz_0[j] = pa_y[j] * tdy_yz_zzzz_0[j] + 0.5 * fl1_fx * tdy_z_zzzz_0[j] + 0.5 * fl1_fx * ts_yz_zzzz_0[j];

            tdz_yyz_zzzz_0[j] = pa_y[j] * tdz_yz_zzzz_0[j] + 0.5 * fl1_fx * tdz_z_zzzz_0[j];

            tdx_yzz_xxxx_0[j] = pa_y[j] * tdx_zz_xxxx_0[j];

            tdy_yzz_xxxx_0[j] = pa_y[j] * tdy_zz_xxxx_0[j] + 0.5 * fl1_fx * ts_zz_xxxx_0[j];

            tdz_yzz_xxxx_0[j] = pa_y[j] * tdz_zz_xxxx_0[j];

            tdx_yzz_xxxy_0[j] = pa_y[j] * tdx_zz_xxxy_0[j] + 0.5 * fl1_fx * tdx_zz_xxx_0[j];

            tdy_yzz_xxxy_0[j] = pa_y[j] * tdy_zz_xxxy_0[j] + 0.5 * fl1_fx * tdy_zz_xxx_0[j] + 0.5 * fl1_fx * ts_zz_xxxy_0[j];

            tdz_yzz_xxxy_0[j] = pa_y[j] * tdz_zz_xxxy_0[j] + 0.5 * fl1_fx * tdz_zz_xxx_0[j];

            tdx_yzz_xxxz_0[j] = pa_y[j] * tdx_zz_xxxz_0[j];

            tdy_yzz_xxxz_0[j] = pa_y[j] * tdy_zz_xxxz_0[j] + 0.5 * fl1_fx * ts_zz_xxxz_0[j];

            tdz_yzz_xxxz_0[j] = pa_y[j] * tdz_zz_xxxz_0[j];

            tdx_yzz_xxyy_0[j] = pa_y[j] * tdx_zz_xxyy_0[j] + fl1_fx * tdx_zz_xxy_0[j];

            tdy_yzz_xxyy_0[j] = pa_y[j] * tdy_zz_xxyy_0[j] + fl1_fx * tdy_zz_xxy_0[j] + 0.5 * fl1_fx * ts_zz_xxyy_0[j];

            tdz_yzz_xxyy_0[j] = pa_y[j] * tdz_zz_xxyy_0[j] + fl1_fx * tdz_zz_xxy_0[j];

            tdx_yzz_xxyz_0[j] = pa_y[j] * tdx_zz_xxyz_0[j] + 0.5 * fl1_fx * tdx_zz_xxz_0[j];

            tdy_yzz_xxyz_0[j] = pa_y[j] * tdy_zz_xxyz_0[j] + 0.5 * fl1_fx * tdy_zz_xxz_0[j] + 0.5 * fl1_fx * ts_zz_xxyz_0[j];

            tdz_yzz_xxyz_0[j] = pa_y[j] * tdz_zz_xxyz_0[j] + 0.5 * fl1_fx * tdz_zz_xxz_0[j];

            tdx_yzz_xxzz_0[j] = pa_y[j] * tdx_zz_xxzz_0[j];

            tdy_yzz_xxzz_0[j] = pa_y[j] * tdy_zz_xxzz_0[j] + 0.5 * fl1_fx * ts_zz_xxzz_0[j];

            tdz_yzz_xxzz_0[j] = pa_y[j] * tdz_zz_xxzz_0[j];

            tdx_yzz_xyyy_0[j] = pa_y[j] * tdx_zz_xyyy_0[j] + 1.5 * fl1_fx * tdx_zz_xyy_0[j];

            tdy_yzz_xyyy_0[j] = pa_y[j] * tdy_zz_xyyy_0[j] + 1.5 * fl1_fx * tdy_zz_xyy_0[j] + 0.5 * fl1_fx * ts_zz_xyyy_0[j];

            tdz_yzz_xyyy_0[j] = pa_y[j] * tdz_zz_xyyy_0[j] + 1.5 * fl1_fx * tdz_zz_xyy_0[j];

            tdx_yzz_xyyz_0[j] = pa_y[j] * tdx_zz_xyyz_0[j] + fl1_fx * tdx_zz_xyz_0[j];

            tdy_yzz_xyyz_0[j] = pa_y[j] * tdy_zz_xyyz_0[j] + fl1_fx * tdy_zz_xyz_0[j] + 0.5 * fl1_fx * ts_zz_xyyz_0[j];

            tdz_yzz_xyyz_0[j] = pa_y[j] * tdz_zz_xyyz_0[j] + fl1_fx * tdz_zz_xyz_0[j];

            tdx_yzz_xyzz_0[j] = pa_y[j] * tdx_zz_xyzz_0[j] + 0.5 * fl1_fx * tdx_zz_xzz_0[j];

            tdy_yzz_xyzz_0[j] = pa_y[j] * tdy_zz_xyzz_0[j] + 0.5 * fl1_fx * tdy_zz_xzz_0[j] + 0.5 * fl1_fx * ts_zz_xyzz_0[j];

            tdz_yzz_xyzz_0[j] = pa_y[j] * tdz_zz_xyzz_0[j] + 0.5 * fl1_fx * tdz_zz_xzz_0[j];

            tdx_yzz_xzzz_0[j] = pa_y[j] * tdx_zz_xzzz_0[j];

            tdy_yzz_xzzz_0[j] = pa_y[j] * tdy_zz_xzzz_0[j] + 0.5 * fl1_fx * ts_zz_xzzz_0[j];

            tdz_yzz_xzzz_0[j] = pa_y[j] * tdz_zz_xzzz_0[j];

            tdx_yzz_yyyy_0[j] = pa_y[j] * tdx_zz_yyyy_0[j] + 2.0 * fl1_fx * tdx_zz_yyy_0[j];

            tdy_yzz_yyyy_0[j] = pa_y[j] * tdy_zz_yyyy_0[j] + 2.0 * fl1_fx * tdy_zz_yyy_0[j] + 0.5 * fl1_fx * ts_zz_yyyy_0[j];

            tdz_yzz_yyyy_0[j] = pa_y[j] * tdz_zz_yyyy_0[j] + 2.0 * fl1_fx * tdz_zz_yyy_0[j];

            tdx_yzz_yyyz_0[j] = pa_y[j] * tdx_zz_yyyz_0[j] + 1.5 * fl1_fx * tdx_zz_yyz_0[j];

            tdy_yzz_yyyz_0[j] = pa_y[j] * tdy_zz_yyyz_0[j] + 1.5 * fl1_fx * tdy_zz_yyz_0[j] + 0.5 * fl1_fx * ts_zz_yyyz_0[j];

            tdz_yzz_yyyz_0[j] = pa_y[j] * tdz_zz_yyyz_0[j] + 1.5 * fl1_fx * tdz_zz_yyz_0[j];

            tdx_yzz_yyzz_0[j] = pa_y[j] * tdx_zz_yyzz_0[j] + fl1_fx * tdx_zz_yzz_0[j];

            tdy_yzz_yyzz_0[j] = pa_y[j] * tdy_zz_yyzz_0[j] + fl1_fx * tdy_zz_yzz_0[j] + 0.5 * fl1_fx * ts_zz_yyzz_0[j];

            tdz_yzz_yyzz_0[j] = pa_y[j] * tdz_zz_yyzz_0[j] + fl1_fx * tdz_zz_yzz_0[j];

            tdx_yzz_yzzz_0[j] = pa_y[j] * tdx_zz_yzzz_0[j] + 0.5 * fl1_fx * tdx_zz_zzz_0[j];
        }

        idx++;
    }
}

void
compElectricDipoleForFG_400_450(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const int32_t              nOSFactors,
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

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_d_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_y = paDistances.data(3 * idx + 1);

        auto pa_z = paDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto tdx_zz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 75);

        auto tdy_zz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 75);

        auto tdz_zz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 75);

        auto tdx_zz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 76);

        auto tdy_zz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 76);

        auto tdz_zz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 76);

        auto tdx_zz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 77);

        auto tdy_zz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 77);

        auto tdz_zz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 77);

        auto tdx_zz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 78);

        auto tdy_zz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 78);

        auto tdz_zz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 78);

        auto tdx_zz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 79);

        auto tdy_zz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 79);

        auto tdz_zz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 79);

        auto tdx_zz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 80);

        auto tdy_zz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 80);

        auto tdz_zz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 80);

        auto tdx_zz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 81);

        auto tdy_zz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 81);

        auto tdz_zz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 81);

        auto tdx_zz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 82);

        auto tdy_zz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 82);

        auto tdz_zz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 82);

        auto tdx_zz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 83);

        auto tdy_zz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 83);

        auto tdz_zz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 83);

        auto tdx_zz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 84);

        auto tdy_zz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 84);

        auto tdz_zz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 84);

        auto tdx_zz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 85);

        auto tdy_zz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 85);

        auto tdz_zz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 85);

        auto tdx_zz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 86);

        auto tdy_zz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 86);

        auto tdz_zz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 86);

        auto tdx_zz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 87);

        auto tdy_zz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 87);

        auto tdz_zz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 87);

        auto tdx_zz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 88);

        auto tdy_zz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 88);

        auto tdz_zz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 88);

        auto tdx_zz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 89);

        auto tdy_zz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 89);

        auto tdz_zz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 89);

        auto tdx_z_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 30);

        auto tdy_z_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 30);

        auto tdz_z_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 30);

        auto tdx_z_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 31);

        auto tdy_z_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 31);

        auto tdz_z_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 31);

        auto tdx_z_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 32);

        auto tdy_z_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 32);

        auto tdz_z_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 32);

        auto tdx_z_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 33);

        auto tdy_z_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 33);

        auto tdz_z_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 33);

        auto tdx_z_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 34);

        auto tdy_z_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 34);

        auto tdz_z_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 34);

        auto tdx_z_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 35);

        auto tdy_z_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 35);

        auto tdz_z_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 35);

        auto tdx_z_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 36);

        auto tdy_z_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 36);

        auto tdz_z_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 36);

        auto tdx_z_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 37);

        auto tdy_z_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 37);

        auto tdz_z_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 37);

        auto tdx_z_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 38);

        auto tdy_z_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 38);

        auto tdz_z_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 38);

        auto tdx_z_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 39);

        auto tdy_z_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 39);

        auto tdz_z_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 39);

        auto tdx_z_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 40);

        auto tdy_z_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 40);

        auto tdz_z_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 40);

        auto tdx_z_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 41);

        auto tdy_z_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 41);

        auto tdz_z_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 41);

        auto tdx_z_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 42);

        auto tdy_z_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 42);

        auto tdz_z_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 42);

        auto tdx_z_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 43);

        auto tdy_z_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 43);

        auto tdz_z_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 43);

        auto tdx_z_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 44);

        auto tdy_z_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 44);

        auto tdz_z_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 44);

        auto tdx_zz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 50);

        auto tdy_zz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 50);

        auto tdz_zz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 50);

        auto tdx_zz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 51);

        auto tdy_zz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 51);

        auto tdz_zz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 51);

        auto tdx_zz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 52);

        auto tdy_zz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 52);

        auto tdz_zz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 52);

        auto tdx_zz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 53);

        auto tdy_zz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 53);

        auto tdz_zz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 53);

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

        auto ts_zz_xxxx_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 75);

        auto ts_zz_xxxy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 76);

        auto ts_zz_xxxz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 77);

        auto ts_zz_xxyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 78);

        auto ts_zz_xxyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 79);

        auto ts_zz_xxzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 80);

        auto ts_zz_xyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 81);

        auto ts_zz_xyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 82);

        auto ts_zz_xyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 83);

        auto ts_zz_xzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 84);

        auto ts_zz_yyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 85);

        auto ts_zz_yyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 86);

        auto ts_zz_yyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 87);

        auto ts_zz_yzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 88);

        auto ts_zz_zzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 89);

        // set up pointers to integrals

        auto tdy_yzz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 133);

        auto tdz_yzz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 133);

        auto tdx_yzz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 134);

        auto tdy_yzz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 134);

        auto tdz_yzz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 134);

        auto tdx_zzz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 135);

        auto tdy_zzz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 135);

        auto tdz_zzz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 135);

        auto tdx_zzz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 136);

        auto tdy_zzz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 136);

        auto tdz_zzz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 136);

        auto tdx_zzz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 137);

        auto tdy_zzz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 137);

        auto tdz_zzz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 137);

        auto tdx_zzz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 138);

        auto tdy_zzz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 138);

        auto tdz_zzz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 138);

        auto tdx_zzz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 139);

        auto tdy_zzz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 139);

        auto tdz_zzz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 139);

        auto tdx_zzz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 140);

        auto tdy_zzz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 140);

        auto tdz_zzz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 140);

        auto tdx_zzz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 141);

        auto tdy_zzz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 141);

        auto tdz_zzz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 141);

        auto tdx_zzz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 142);

        auto tdy_zzz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 142);

        auto tdz_zzz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 142);

        auto tdx_zzz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 143);

        auto tdy_zzz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 143);

        auto tdz_zzz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 143);

        auto tdx_zzz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 144);

        auto tdy_zzz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 144);

        auto tdz_zzz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 144);

        auto tdx_zzz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 145);

        auto tdy_zzz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 145);

        auto tdz_zzz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 145);

        auto tdx_zzz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 146);

        auto tdy_zzz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 146);

        auto tdz_zzz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 146);

        auto tdx_zzz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 147);

        auto tdy_zzz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 147);

        auto tdz_zzz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 147);

        auto tdx_zzz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 148);

        auto tdy_zzz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 148);

        auto tdz_zzz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 148);

        auto tdx_zzz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 149);

        auto tdy_zzz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 149);

        auto tdz_zzz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 149);

        // Batch of Integrals (400,450)

        #pragma omp simd aligned(fx, pa_y, pa_z, tdx_yzz_zzzz_0, tdx_z_xxxx_0, tdx_z_xxxy_0, tdx_z_xxxz_0, \
                                     tdx_z_xxyy_0, tdx_z_xxyz_0, tdx_z_xxzz_0, tdx_z_xyyy_0, tdx_z_xyyz_0, tdx_z_xyzz_0, \
                                     tdx_z_xzzz_0, tdx_z_yyyy_0, tdx_z_yyyz_0, tdx_z_yyzz_0, tdx_z_yzzz_0, tdx_z_zzzz_0, \
                                     tdx_zz_xxx_0, tdx_zz_xxxx_0, tdx_zz_xxxy_0, tdx_zz_xxxz_0, tdx_zz_xxy_0, \
                                     tdx_zz_xxyy_0, tdx_zz_xxyz_0, tdx_zz_xxz_0, tdx_zz_xxzz_0, tdx_zz_xyy_0, \
                                     tdx_zz_xyyy_0, tdx_zz_xyyz_0, tdx_zz_xyz_0, tdx_zz_xyzz_0, tdx_zz_xzz_0, \
                                     tdx_zz_xzzz_0, tdx_zz_yyy_0, tdx_zz_yyyy_0, tdx_zz_yyyz_0, tdx_zz_yyz_0, \
                                     tdx_zz_yyzz_0, tdx_zz_yzz_0, tdx_zz_yzzz_0, tdx_zz_zzz_0, tdx_zz_zzzz_0, \
                                     tdx_zzz_xxxx_0, tdx_zzz_xxxy_0, tdx_zzz_xxxz_0, tdx_zzz_xxyy_0, tdx_zzz_xxyz_0, \
                                     tdx_zzz_xxzz_0, tdx_zzz_xyyy_0, tdx_zzz_xyyz_0, tdx_zzz_xyzz_0, tdx_zzz_xzzz_0, \
                                     tdx_zzz_yyyy_0, tdx_zzz_yyyz_0, tdx_zzz_yyzz_0, tdx_zzz_yzzz_0, tdx_zzz_zzzz_0, \
                                     tdy_yzz_yzzz_0, tdy_yzz_zzzz_0, tdy_z_xxxx_0, tdy_z_xxxy_0, tdy_z_xxxz_0, \
                                     tdy_z_xxyy_0, tdy_z_xxyz_0, tdy_z_xxzz_0, tdy_z_xyyy_0, tdy_z_xyyz_0, tdy_z_xyzz_0, \
                                     tdy_z_xzzz_0, tdy_z_yyyy_0, tdy_z_yyyz_0, tdy_z_yyzz_0, tdy_z_yzzz_0, tdy_z_zzzz_0, \
                                     tdy_zz_xxx_0, tdy_zz_xxxx_0, tdy_zz_xxxy_0, tdy_zz_xxxz_0, tdy_zz_xxy_0, \
                                     tdy_zz_xxyy_0, tdy_zz_xxyz_0, tdy_zz_xxz_0, tdy_zz_xxzz_0, tdy_zz_xyy_0, \
                                     tdy_zz_xyyy_0, tdy_zz_xyyz_0, tdy_zz_xyz_0, tdy_zz_xyzz_0, tdy_zz_xzz_0, \
                                     tdy_zz_xzzz_0, tdy_zz_yyy_0, tdy_zz_yyyy_0, tdy_zz_yyyz_0, tdy_zz_yyz_0, \
                                     tdy_zz_yyzz_0, tdy_zz_yzz_0, tdy_zz_yzzz_0, tdy_zz_zzz_0, tdy_zz_zzzz_0, \
                                     tdy_zzz_xxxx_0, tdy_zzz_xxxy_0, tdy_zzz_xxxz_0, tdy_zzz_xxyy_0, tdy_zzz_xxyz_0, \
                                     tdy_zzz_xxzz_0, tdy_zzz_xyyy_0, tdy_zzz_xyyz_0, tdy_zzz_xyzz_0, tdy_zzz_xzzz_0, \
                                     tdy_zzz_yyyy_0, tdy_zzz_yyyz_0, tdy_zzz_yyzz_0, tdy_zzz_yzzz_0, tdy_zzz_zzzz_0, \
                                     tdz_yzz_yzzz_0, tdz_yzz_zzzz_0, tdz_z_xxxx_0, tdz_z_xxxy_0, tdz_z_xxxz_0, \
                                     tdz_z_xxyy_0, tdz_z_xxyz_0, tdz_z_xxzz_0, tdz_z_xyyy_0, tdz_z_xyyz_0, tdz_z_xyzz_0, \
                                     tdz_z_xzzz_0, tdz_z_yyyy_0, tdz_z_yyyz_0, tdz_z_yyzz_0, tdz_z_yzzz_0, tdz_z_zzzz_0, \
                                     tdz_zz_xxx_0, tdz_zz_xxxx_0, tdz_zz_xxxy_0, tdz_zz_xxxz_0, tdz_zz_xxy_0, \
                                     tdz_zz_xxyy_0, tdz_zz_xxyz_0, tdz_zz_xxz_0, tdz_zz_xxzz_0, tdz_zz_xyy_0, \
                                     tdz_zz_xyyy_0, tdz_zz_xyyz_0, tdz_zz_xyz_0, tdz_zz_xyzz_0, tdz_zz_xzz_0, \
                                     tdz_zz_xzzz_0, tdz_zz_yyy_0, tdz_zz_yyyy_0, tdz_zz_yyyz_0, tdz_zz_yyz_0, \
                                     tdz_zz_yyzz_0, tdz_zz_yzz_0, tdz_zz_yzzz_0, tdz_zz_zzz_0, tdz_zz_zzzz_0, \
                                     tdz_zzz_xxxx_0, tdz_zzz_xxxy_0, tdz_zzz_xxxz_0, tdz_zzz_xxyy_0, tdz_zzz_xxyz_0, \
                                     tdz_zzz_xxzz_0, tdz_zzz_xyyy_0, tdz_zzz_xyyz_0, tdz_zzz_xyzz_0, tdz_zzz_xzzz_0, \
                                     tdz_zzz_yyyy_0, tdz_zzz_yyyz_0, tdz_zzz_yyzz_0, tdz_zzz_yzzz_0, tdz_zzz_zzzz_0, \
                                     ts_zz_xxxx_0, ts_zz_xxxy_0, ts_zz_xxxz_0, ts_zz_xxyy_0, ts_zz_xxyz_0, ts_zz_xxzz_0, \
                                     ts_zz_xyyy_0, ts_zz_xyyz_0, ts_zz_xyzz_0, ts_zz_xzzz_0, ts_zz_yyyy_0, ts_zz_yyyz_0, \
                                     ts_zz_yyzz_0, ts_zz_yzzz_0, ts_zz_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            tdy_yzz_yzzz_0[j] = pa_y[j] * tdy_zz_yzzz_0[j] + 0.5 * fl1_fx * tdy_zz_zzz_0[j] + 0.5 * fl1_fx * ts_zz_yzzz_0[j];

            tdz_yzz_yzzz_0[j] = pa_y[j] * tdz_zz_yzzz_0[j] + 0.5 * fl1_fx * tdz_zz_zzz_0[j];

            tdx_yzz_zzzz_0[j] = pa_y[j] * tdx_zz_zzzz_0[j];

            tdy_yzz_zzzz_0[j] = pa_y[j] * tdy_zz_zzzz_0[j] + 0.5 * fl1_fx * ts_zz_zzzz_0[j];

            tdz_yzz_zzzz_0[j] = pa_y[j] * tdz_zz_zzzz_0[j];

            tdx_zzz_xxxx_0[j] = pa_z[j] * tdx_zz_xxxx_0[j] + fl1_fx * tdx_z_xxxx_0[j];

            tdy_zzz_xxxx_0[j] = pa_z[j] * tdy_zz_xxxx_0[j] + fl1_fx * tdy_z_xxxx_0[j];

            tdz_zzz_xxxx_0[j] = pa_z[j] * tdz_zz_xxxx_0[j] + fl1_fx * tdz_z_xxxx_0[j] + 0.5 * fl1_fx * ts_zz_xxxx_0[j];

            tdx_zzz_xxxy_0[j] = pa_z[j] * tdx_zz_xxxy_0[j] + fl1_fx * tdx_z_xxxy_0[j];

            tdy_zzz_xxxy_0[j] = pa_z[j] * tdy_zz_xxxy_0[j] + fl1_fx * tdy_z_xxxy_0[j];

            tdz_zzz_xxxy_0[j] = pa_z[j] * tdz_zz_xxxy_0[j] + fl1_fx * tdz_z_xxxy_0[j] + 0.5 * fl1_fx * ts_zz_xxxy_0[j];

            tdx_zzz_xxxz_0[j] = pa_z[j] * tdx_zz_xxxz_0[j] + fl1_fx * tdx_z_xxxz_0[j] + 0.5 * fl1_fx * tdx_zz_xxx_0[j];

            tdy_zzz_xxxz_0[j] = pa_z[j] * tdy_zz_xxxz_0[j] + fl1_fx * tdy_z_xxxz_0[j] + 0.5 * fl1_fx * tdy_zz_xxx_0[j];

            tdz_zzz_xxxz_0[j] =
                pa_z[j] * tdz_zz_xxxz_0[j] + fl1_fx * tdz_z_xxxz_0[j] + 0.5 * fl1_fx * tdz_zz_xxx_0[j] + 0.5 * fl1_fx * ts_zz_xxxz_0[j];

            tdx_zzz_xxyy_0[j] = pa_z[j] * tdx_zz_xxyy_0[j] + fl1_fx * tdx_z_xxyy_0[j];

            tdy_zzz_xxyy_0[j] = pa_z[j] * tdy_zz_xxyy_0[j] + fl1_fx * tdy_z_xxyy_0[j];

            tdz_zzz_xxyy_0[j] = pa_z[j] * tdz_zz_xxyy_0[j] + fl1_fx * tdz_z_xxyy_0[j] + 0.5 * fl1_fx * ts_zz_xxyy_0[j];

            tdx_zzz_xxyz_0[j] = pa_z[j] * tdx_zz_xxyz_0[j] + fl1_fx * tdx_z_xxyz_0[j] + 0.5 * fl1_fx * tdx_zz_xxy_0[j];

            tdy_zzz_xxyz_0[j] = pa_z[j] * tdy_zz_xxyz_0[j] + fl1_fx * tdy_z_xxyz_0[j] + 0.5 * fl1_fx * tdy_zz_xxy_0[j];

            tdz_zzz_xxyz_0[j] =
                pa_z[j] * tdz_zz_xxyz_0[j] + fl1_fx * tdz_z_xxyz_0[j] + 0.5 * fl1_fx * tdz_zz_xxy_0[j] + 0.5 * fl1_fx * ts_zz_xxyz_0[j];

            tdx_zzz_xxzz_0[j] = pa_z[j] * tdx_zz_xxzz_0[j] + fl1_fx * tdx_z_xxzz_0[j] + fl1_fx * tdx_zz_xxz_0[j];

            tdy_zzz_xxzz_0[j] = pa_z[j] * tdy_zz_xxzz_0[j] + fl1_fx * tdy_z_xxzz_0[j] + fl1_fx * tdy_zz_xxz_0[j];

            tdz_zzz_xxzz_0[j] = pa_z[j] * tdz_zz_xxzz_0[j] + fl1_fx * tdz_z_xxzz_0[j] + fl1_fx * tdz_zz_xxz_0[j] + 0.5 * fl1_fx * ts_zz_xxzz_0[j];

            tdx_zzz_xyyy_0[j] = pa_z[j] * tdx_zz_xyyy_0[j] + fl1_fx * tdx_z_xyyy_0[j];

            tdy_zzz_xyyy_0[j] = pa_z[j] * tdy_zz_xyyy_0[j] + fl1_fx * tdy_z_xyyy_0[j];

            tdz_zzz_xyyy_0[j] = pa_z[j] * tdz_zz_xyyy_0[j] + fl1_fx * tdz_z_xyyy_0[j] + 0.5 * fl1_fx * ts_zz_xyyy_0[j];

            tdx_zzz_xyyz_0[j] = pa_z[j] * tdx_zz_xyyz_0[j] + fl1_fx * tdx_z_xyyz_0[j] + 0.5 * fl1_fx * tdx_zz_xyy_0[j];

            tdy_zzz_xyyz_0[j] = pa_z[j] * tdy_zz_xyyz_0[j] + fl1_fx * tdy_z_xyyz_0[j] + 0.5 * fl1_fx * tdy_zz_xyy_0[j];

            tdz_zzz_xyyz_0[j] =
                pa_z[j] * tdz_zz_xyyz_0[j] + fl1_fx * tdz_z_xyyz_0[j] + 0.5 * fl1_fx * tdz_zz_xyy_0[j] + 0.5 * fl1_fx * ts_zz_xyyz_0[j];

            tdx_zzz_xyzz_0[j] = pa_z[j] * tdx_zz_xyzz_0[j] + fl1_fx * tdx_z_xyzz_0[j] + fl1_fx * tdx_zz_xyz_0[j];

            tdy_zzz_xyzz_0[j] = pa_z[j] * tdy_zz_xyzz_0[j] + fl1_fx * tdy_z_xyzz_0[j] + fl1_fx * tdy_zz_xyz_0[j];

            tdz_zzz_xyzz_0[j] = pa_z[j] * tdz_zz_xyzz_0[j] + fl1_fx * tdz_z_xyzz_0[j] + fl1_fx * tdz_zz_xyz_0[j] + 0.5 * fl1_fx * ts_zz_xyzz_0[j];

            tdx_zzz_xzzz_0[j] = pa_z[j] * tdx_zz_xzzz_0[j] + fl1_fx * tdx_z_xzzz_0[j] + 1.5 * fl1_fx * tdx_zz_xzz_0[j];

            tdy_zzz_xzzz_0[j] = pa_z[j] * tdy_zz_xzzz_0[j] + fl1_fx * tdy_z_xzzz_0[j] + 1.5 * fl1_fx * tdy_zz_xzz_0[j];

            tdz_zzz_xzzz_0[j] =
                pa_z[j] * tdz_zz_xzzz_0[j] + fl1_fx * tdz_z_xzzz_0[j] + 1.5 * fl1_fx * tdz_zz_xzz_0[j] + 0.5 * fl1_fx * ts_zz_xzzz_0[j];

            tdx_zzz_yyyy_0[j] = pa_z[j] * tdx_zz_yyyy_0[j] + fl1_fx * tdx_z_yyyy_0[j];

            tdy_zzz_yyyy_0[j] = pa_z[j] * tdy_zz_yyyy_0[j] + fl1_fx * tdy_z_yyyy_0[j];

            tdz_zzz_yyyy_0[j] = pa_z[j] * tdz_zz_yyyy_0[j] + fl1_fx * tdz_z_yyyy_0[j] + 0.5 * fl1_fx * ts_zz_yyyy_0[j];

            tdx_zzz_yyyz_0[j] = pa_z[j] * tdx_zz_yyyz_0[j] + fl1_fx * tdx_z_yyyz_0[j] + 0.5 * fl1_fx * tdx_zz_yyy_0[j];

            tdy_zzz_yyyz_0[j] = pa_z[j] * tdy_zz_yyyz_0[j] + fl1_fx * tdy_z_yyyz_0[j] + 0.5 * fl1_fx * tdy_zz_yyy_0[j];

            tdz_zzz_yyyz_0[j] =
                pa_z[j] * tdz_zz_yyyz_0[j] + fl1_fx * tdz_z_yyyz_0[j] + 0.5 * fl1_fx * tdz_zz_yyy_0[j] + 0.5 * fl1_fx * ts_zz_yyyz_0[j];

            tdx_zzz_yyzz_0[j] = pa_z[j] * tdx_zz_yyzz_0[j] + fl1_fx * tdx_z_yyzz_0[j] + fl1_fx * tdx_zz_yyz_0[j];

            tdy_zzz_yyzz_0[j] = pa_z[j] * tdy_zz_yyzz_0[j] + fl1_fx * tdy_z_yyzz_0[j] + fl1_fx * tdy_zz_yyz_0[j];

            tdz_zzz_yyzz_0[j] = pa_z[j] * tdz_zz_yyzz_0[j] + fl1_fx * tdz_z_yyzz_0[j] + fl1_fx * tdz_zz_yyz_0[j] + 0.5 * fl1_fx * ts_zz_yyzz_0[j];

            tdx_zzz_yzzz_0[j] = pa_z[j] * tdx_zz_yzzz_0[j] + fl1_fx * tdx_z_yzzz_0[j] + 1.5 * fl1_fx * tdx_zz_yzz_0[j];

            tdy_zzz_yzzz_0[j] = pa_z[j] * tdy_zz_yzzz_0[j] + fl1_fx * tdy_z_yzzz_0[j] + 1.5 * fl1_fx * tdy_zz_yzz_0[j];

            tdz_zzz_yzzz_0[j] =
                pa_z[j] * tdz_zz_yzzz_0[j] + fl1_fx * tdz_z_yzzz_0[j] + 1.5 * fl1_fx * tdz_zz_yzz_0[j] + 0.5 * fl1_fx * ts_zz_yzzz_0[j];

            tdx_zzz_zzzz_0[j] = pa_z[j] * tdx_zz_zzzz_0[j] + fl1_fx * tdx_z_zzzz_0[j] + 2.0 * fl1_fx * tdx_zz_zzz_0[j];

            tdy_zzz_zzzz_0[j] = pa_z[j] * tdy_zz_zzzz_0[j] + fl1_fx * tdy_z_zzzz_0[j] + 2.0 * fl1_fx * tdy_zz_zzz_0[j];

            tdz_zzz_zzzz_0[j] =
                pa_z[j] * tdz_zz_zzzz_0[j] + fl1_fx * tdz_z_zzzz_0[j] + 2.0 * fl1_fx * tdz_zz_zzz_0[j] + 0.5 * fl1_fx * ts_zz_zzzz_0[j];
        }

        idx++;
    }
}

}  // namespace ediprecfunc
