//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ElectricDipoleRecFuncForGG.hpp"

namespace ediprecfunc {  // ediprecfunc namespace

void
compElectricDipoleForGG(CMemBlock2D<double>&       primBuffer,
                        const CRecursionMap&       recursionMap,
                        const CMemBlock2D<double>& osFactors,
                        const int32_t              nOSFactors,
                        const CMemBlock2D<double>& paDistances,
                        const CGtoBlock&           braGtoBlock,
                        const CGtoBlock&           ketGtoBlock,
                        const int32_t              iContrGto)
{
    ediprecfunc::compElectricDipoleForGG_0_49(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForGG_49_98(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForGG_98_147(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForGG_147_195(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForGG_195_243(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForGG_243_291(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForGG_291_339(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForGG_339_387(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForGG_387_435(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForGG_435_483(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForGG_483_531(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForGG_531_579(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForGG_579_627(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForGG_627_675(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compElectricDipoleForGG_0_49(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const int32_t              nOSFactors,
                             const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_d_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_d_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

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

        auto tdx_xxx_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx);

        auto tdy_xxx_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx);

        auto tdz_xxx_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx);

        auto tdx_xxx_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 1);

        auto tdy_xxx_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 1);

        auto tdz_xxx_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 1);

        auto tdx_xxx_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 2);

        auto tdy_xxx_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 2);

        auto tdz_xxx_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 2);

        auto tdx_xxx_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 3);

        auto tdy_xxx_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 3);

        auto tdz_xxx_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 3);

        auto tdx_xxx_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 4);

        auto tdy_xxx_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 4);

        auto tdz_xxx_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 4);

        auto tdx_xxx_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 5);

        auto tdy_xxx_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 5);

        auto tdz_xxx_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 5);

        auto tdx_xxx_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 6);

        auto tdy_xxx_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 6);

        auto tdz_xxx_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 6);

        auto tdx_xxx_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 7);

        auto tdy_xxx_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 7);

        auto tdz_xxx_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 7);

        auto tdx_xxx_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 8);

        auto tdy_xxx_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 8);

        auto tdz_xxx_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 8);

        auto tdx_xxx_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 9);

        auto tdy_xxx_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 9);

        auto tdz_xxx_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 9);

        auto tdx_xxy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 10);

        auto tdy_xxy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 10);

        auto tdz_xxy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 10);

        auto tdx_xxy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 11);

        auto ts_xxx_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx);

        auto ts_xxx_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 1);

        auto ts_xxx_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 2);

        auto ts_xxx_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 3);

        auto ts_xxx_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 4);

        auto ts_xxx_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 5);

        auto ts_xxx_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 6);

        auto ts_xxx_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 7);

        auto ts_xxx_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 8);

        auto ts_xxx_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 9);

        auto ts_xxx_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 10);

        auto ts_xxx_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 11);

        auto ts_xxx_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 12);

        auto ts_xxx_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 13);

        auto ts_xxx_zzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 14);

        auto ts_xxy_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 15);

        auto ts_xxy_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 16);

        // set up pointers to integrals

        auto tdx_xxxx_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx);

        auto tdy_xxxx_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx);

        auto tdz_xxxx_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx);

        auto tdx_xxxx_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 1);

        auto tdy_xxxx_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 1);

        auto tdz_xxxx_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 1);

        auto tdx_xxxx_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 2);

        auto tdy_xxxx_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 2);

        auto tdz_xxxx_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 2);

        auto tdx_xxxx_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 3);

        auto tdy_xxxx_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 3);

        auto tdz_xxxx_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 3);

        auto tdx_xxxx_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 4);

        auto tdy_xxxx_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 4);

        auto tdz_xxxx_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 4);

        auto tdx_xxxx_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 5);

        auto tdy_xxxx_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 5);

        auto tdz_xxxx_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 5);

        auto tdx_xxxx_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 6);

        auto tdy_xxxx_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 6);

        auto tdz_xxxx_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 6);

        auto tdx_xxxx_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 7);

        auto tdy_xxxx_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 7);

        auto tdz_xxxx_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 7);

        auto tdx_xxxx_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 8);

        auto tdy_xxxx_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 8);

        auto tdz_xxxx_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 8);

        auto tdx_xxxx_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 9);

        auto tdy_xxxx_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 9);

        auto tdz_xxxx_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 9);

        auto tdx_xxxx_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 10);

        auto tdy_xxxx_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 10);

        auto tdz_xxxx_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 10);

        auto tdx_xxxx_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 11);

        auto tdy_xxxx_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 11);

        auto tdz_xxxx_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 11);

        auto tdx_xxxx_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 12);

        auto tdy_xxxx_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 12);

        auto tdz_xxxx_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 12);

        auto tdx_xxxx_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 13);

        auto tdy_xxxx_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 13);

        auto tdz_xxxx_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 13);

        auto tdx_xxxx_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 14);

        auto tdy_xxxx_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 14);

        auto tdz_xxxx_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 14);

        auto tdx_xxxy_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 15);

        auto tdy_xxxy_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 15);

        auto tdz_xxxy_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 15);

        auto tdx_xxxy_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 16);

        // Batch of Integrals (0,49)

        #pragma omp simd aligned(fx, pa_x, tdx_xx_xxxx_0, tdx_xx_xxxy_0, tdx_xx_xxxz_0, tdx_xx_xxyy_0, \
                                     tdx_xx_xxyz_0, tdx_xx_xxzz_0, tdx_xx_xyyy_0, tdx_xx_xyyz_0, tdx_xx_xyzz_0, \
                                     tdx_xx_xzzz_0, tdx_xx_yyyy_0, tdx_xx_yyyz_0, tdx_xx_yyzz_0, tdx_xx_yzzz_0, \
                                     tdx_xx_zzzz_0, tdx_xxx_xxx_0, tdx_xxx_xxxx_0, tdx_xxx_xxxy_0, tdx_xxx_xxxz_0, \
                                     tdx_xxx_xxy_0, tdx_xxx_xxyy_0, tdx_xxx_xxyz_0, tdx_xxx_xxz_0, tdx_xxx_xxzz_0, \
                                     tdx_xxx_xyy_0, tdx_xxx_xyyy_0, tdx_xxx_xyyz_0, tdx_xxx_xyz_0, tdx_xxx_xyzz_0, \
                                     tdx_xxx_xzz_0, tdx_xxx_xzzz_0, tdx_xxx_yyy_0, tdx_xxx_yyyy_0, tdx_xxx_yyyz_0, \
                                     tdx_xxx_yyz_0, tdx_xxx_yyzz_0, tdx_xxx_yzz_0, tdx_xxx_yzzz_0, tdx_xxx_zzz_0, \
                                     tdx_xxx_zzzz_0, tdx_xxxx_xxxx_0, tdx_xxxx_xxxy_0, tdx_xxxx_xxxz_0, tdx_xxxx_xxyy_0, \
                                     tdx_xxxx_xxyz_0, tdx_xxxx_xxzz_0, tdx_xxxx_xyyy_0, tdx_xxxx_xyyz_0, tdx_xxxx_xyzz_0, \
                                     tdx_xxxx_xzzz_0, tdx_xxxx_yyyy_0, tdx_xxxx_yyyz_0, tdx_xxxx_yyzz_0, tdx_xxxx_yzzz_0, \
                                     tdx_xxxx_zzzz_0, tdx_xxxy_xxxx_0, tdx_xxxy_xxxy_0, tdx_xxy_xxx_0, tdx_xxy_xxxx_0, \
                                     tdx_xxy_xxxy_0, tdx_xxy_xxy_0, tdx_xy_xxxx_0, tdx_xy_xxxy_0, tdy_xx_xxxx_0, \
                                     tdy_xx_xxxy_0, tdy_xx_xxxz_0, tdy_xx_xxyy_0, tdy_xx_xxyz_0, tdy_xx_xxzz_0, \
                                     tdy_xx_xyyy_0, tdy_xx_xyyz_0, tdy_xx_xyzz_0, tdy_xx_xzzz_0, tdy_xx_yyyy_0, \
                                     tdy_xx_yyyz_0, tdy_xx_yyzz_0, tdy_xx_yzzz_0, tdy_xx_zzzz_0, tdy_xxx_xxx_0, \
                                     tdy_xxx_xxxx_0, tdy_xxx_xxxy_0, tdy_xxx_xxxz_0, tdy_xxx_xxy_0, tdy_xxx_xxyy_0, \
                                     tdy_xxx_xxyz_0, tdy_xxx_xxz_0, tdy_xxx_xxzz_0, tdy_xxx_xyy_0, tdy_xxx_xyyy_0, \
                                     tdy_xxx_xyyz_0, tdy_xxx_xyz_0, tdy_xxx_xyzz_0, tdy_xxx_xzz_0, tdy_xxx_xzzz_0, \
                                     tdy_xxx_yyy_0, tdy_xxx_yyyy_0, tdy_xxx_yyyz_0, tdy_xxx_yyz_0, tdy_xxx_yyzz_0, \
                                     tdy_xxx_yzz_0, tdy_xxx_yzzz_0, tdy_xxx_zzz_0, tdy_xxx_zzzz_0, tdy_xxxx_xxxx_0, \
                                     tdy_xxxx_xxxy_0, tdy_xxxx_xxxz_0, tdy_xxxx_xxyy_0, tdy_xxxx_xxyz_0, tdy_xxxx_xxzz_0, \
                                     tdy_xxxx_xyyy_0, tdy_xxxx_xyyz_0, tdy_xxxx_xyzz_0, tdy_xxxx_xzzz_0, tdy_xxxx_yyyy_0, \
                                     tdy_xxxx_yyyz_0, tdy_xxxx_yyzz_0, tdy_xxxx_yzzz_0, tdy_xxxx_zzzz_0, tdy_xxxy_xxxx_0, \
                                     tdy_xxy_xxx_0, tdy_xxy_xxxx_0, tdy_xy_xxxx_0, tdz_xx_xxxx_0, tdz_xx_xxxy_0, \
                                     tdz_xx_xxxz_0, tdz_xx_xxyy_0, tdz_xx_xxyz_0, tdz_xx_xxzz_0, tdz_xx_xyyy_0, \
                                     tdz_xx_xyyz_0, tdz_xx_xyzz_0, tdz_xx_xzzz_0, tdz_xx_yyyy_0, tdz_xx_yyyz_0, \
                                     tdz_xx_yyzz_0, tdz_xx_yzzz_0, tdz_xx_zzzz_0, tdz_xxx_xxx_0, tdz_xxx_xxxx_0, \
                                     tdz_xxx_xxxy_0, tdz_xxx_xxxz_0, tdz_xxx_xxy_0, tdz_xxx_xxyy_0, tdz_xxx_xxyz_0, \
                                     tdz_xxx_xxz_0, tdz_xxx_xxzz_0, tdz_xxx_xyy_0, tdz_xxx_xyyy_0, tdz_xxx_xyyz_0, \
                                     tdz_xxx_xyz_0, tdz_xxx_xyzz_0, tdz_xxx_xzz_0, tdz_xxx_xzzz_0, tdz_xxx_yyy_0, \
                                     tdz_xxx_yyyy_0, tdz_xxx_yyyz_0, tdz_xxx_yyz_0, tdz_xxx_yyzz_0, tdz_xxx_yzz_0, \
                                     tdz_xxx_yzzz_0, tdz_xxx_zzz_0, tdz_xxx_zzzz_0, tdz_xxxx_xxxx_0, tdz_xxxx_xxxy_0, \
                                     tdz_xxxx_xxxz_0, tdz_xxxx_xxyy_0, tdz_xxxx_xxyz_0, tdz_xxxx_xxzz_0, tdz_xxxx_xyyy_0, \
                                     tdz_xxxx_xyyz_0, tdz_xxxx_xyzz_0, tdz_xxxx_xzzz_0, tdz_xxxx_yyyy_0, tdz_xxxx_yyyz_0, \
                                     tdz_xxxx_yyzz_0, tdz_xxxx_yzzz_0, tdz_xxxx_zzzz_0, tdz_xxxy_xxxx_0, tdz_xxy_xxx_0, \
                                     tdz_xxy_xxxx_0, tdz_xy_xxxx_0, ts_xxx_xxxx_0, ts_xxx_xxxy_0, ts_xxx_xxxz_0, \
                                     ts_xxx_xxyy_0, ts_xxx_xxyz_0, ts_xxx_xxzz_0, ts_xxx_xyyy_0, ts_xxx_xyyz_0, \
                                     ts_xxx_xyzz_0, ts_xxx_xzzz_0, ts_xxx_yyyy_0, ts_xxx_yyyz_0, ts_xxx_yyzz_0, \
                                     ts_xxx_yzzz_0, ts_xxx_zzzz_0, ts_xxy_xxxx_0, ts_xxy_xxxy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            tdx_xxxx_xxxx_0[j] =
                pa_x[j] * tdx_xxx_xxxx_0[j] + 1.5 * fl1_fx * tdx_xx_xxxx_0[j] + 2.0 * fl1_fx * tdx_xxx_xxx_0[j] + 0.5 * fl1_fx * ts_xxx_xxxx_0[j];

            tdy_xxxx_xxxx_0[j] = pa_x[j] * tdy_xxx_xxxx_0[j] + 1.5 * fl1_fx * tdy_xx_xxxx_0[j] + 2.0 * fl1_fx * tdy_xxx_xxx_0[j];

            tdz_xxxx_xxxx_0[j] = pa_x[j] * tdz_xxx_xxxx_0[j] + 1.5 * fl1_fx * tdz_xx_xxxx_0[j] + 2.0 * fl1_fx * tdz_xxx_xxx_0[j];

            tdx_xxxx_xxxy_0[j] =
                pa_x[j] * tdx_xxx_xxxy_0[j] + 1.5 * fl1_fx * tdx_xx_xxxy_0[j] + 1.5 * fl1_fx * tdx_xxx_xxy_0[j] + 0.5 * fl1_fx * ts_xxx_xxxy_0[j];

            tdy_xxxx_xxxy_0[j] = pa_x[j] * tdy_xxx_xxxy_0[j] + 1.5 * fl1_fx * tdy_xx_xxxy_0[j] + 1.5 * fl1_fx * tdy_xxx_xxy_0[j];

            tdz_xxxx_xxxy_0[j] = pa_x[j] * tdz_xxx_xxxy_0[j] + 1.5 * fl1_fx * tdz_xx_xxxy_0[j] + 1.5 * fl1_fx * tdz_xxx_xxy_0[j];

            tdx_xxxx_xxxz_0[j] =
                pa_x[j] * tdx_xxx_xxxz_0[j] + 1.5 * fl1_fx * tdx_xx_xxxz_0[j] + 1.5 * fl1_fx * tdx_xxx_xxz_0[j] + 0.5 * fl1_fx * ts_xxx_xxxz_0[j];

            tdy_xxxx_xxxz_0[j] = pa_x[j] * tdy_xxx_xxxz_0[j] + 1.5 * fl1_fx * tdy_xx_xxxz_0[j] + 1.5 * fl1_fx * tdy_xxx_xxz_0[j];

            tdz_xxxx_xxxz_0[j] = pa_x[j] * tdz_xxx_xxxz_0[j] + 1.5 * fl1_fx * tdz_xx_xxxz_0[j] + 1.5 * fl1_fx * tdz_xxx_xxz_0[j];

            tdx_xxxx_xxyy_0[j] =
                pa_x[j] * tdx_xxx_xxyy_0[j] + 1.5 * fl1_fx * tdx_xx_xxyy_0[j] + fl1_fx * tdx_xxx_xyy_0[j] + 0.5 * fl1_fx * ts_xxx_xxyy_0[j];

            tdy_xxxx_xxyy_0[j] = pa_x[j] * tdy_xxx_xxyy_0[j] + 1.5 * fl1_fx * tdy_xx_xxyy_0[j] + fl1_fx * tdy_xxx_xyy_0[j];

            tdz_xxxx_xxyy_0[j] = pa_x[j] * tdz_xxx_xxyy_0[j] + 1.5 * fl1_fx * tdz_xx_xxyy_0[j] + fl1_fx * tdz_xxx_xyy_0[j];

            tdx_xxxx_xxyz_0[j] =
                pa_x[j] * tdx_xxx_xxyz_0[j] + 1.5 * fl1_fx * tdx_xx_xxyz_0[j] + fl1_fx * tdx_xxx_xyz_0[j] + 0.5 * fl1_fx * ts_xxx_xxyz_0[j];

            tdy_xxxx_xxyz_0[j] = pa_x[j] * tdy_xxx_xxyz_0[j] + 1.5 * fl1_fx * tdy_xx_xxyz_0[j] + fl1_fx * tdy_xxx_xyz_0[j];

            tdz_xxxx_xxyz_0[j] = pa_x[j] * tdz_xxx_xxyz_0[j] + 1.5 * fl1_fx * tdz_xx_xxyz_0[j] + fl1_fx * tdz_xxx_xyz_0[j];

            tdx_xxxx_xxzz_0[j] =
                pa_x[j] * tdx_xxx_xxzz_0[j] + 1.5 * fl1_fx * tdx_xx_xxzz_0[j] + fl1_fx * tdx_xxx_xzz_0[j] + 0.5 * fl1_fx * ts_xxx_xxzz_0[j];

            tdy_xxxx_xxzz_0[j] = pa_x[j] * tdy_xxx_xxzz_0[j] + 1.5 * fl1_fx * tdy_xx_xxzz_0[j] + fl1_fx * tdy_xxx_xzz_0[j];

            tdz_xxxx_xxzz_0[j] = pa_x[j] * tdz_xxx_xxzz_0[j] + 1.5 * fl1_fx * tdz_xx_xxzz_0[j] + fl1_fx * tdz_xxx_xzz_0[j];

            tdx_xxxx_xyyy_0[j] =
                pa_x[j] * tdx_xxx_xyyy_0[j] + 1.5 * fl1_fx * tdx_xx_xyyy_0[j] + 0.5 * fl1_fx * tdx_xxx_yyy_0[j] + 0.5 * fl1_fx * ts_xxx_xyyy_0[j];

            tdy_xxxx_xyyy_0[j] = pa_x[j] * tdy_xxx_xyyy_0[j] + 1.5 * fl1_fx * tdy_xx_xyyy_0[j] + 0.5 * fl1_fx * tdy_xxx_yyy_0[j];

            tdz_xxxx_xyyy_0[j] = pa_x[j] * tdz_xxx_xyyy_0[j] + 1.5 * fl1_fx * tdz_xx_xyyy_0[j] + 0.5 * fl1_fx * tdz_xxx_yyy_0[j];

            tdx_xxxx_xyyz_0[j] =
                pa_x[j] * tdx_xxx_xyyz_0[j] + 1.5 * fl1_fx * tdx_xx_xyyz_0[j] + 0.5 * fl1_fx * tdx_xxx_yyz_0[j] + 0.5 * fl1_fx * ts_xxx_xyyz_0[j];

            tdy_xxxx_xyyz_0[j] = pa_x[j] * tdy_xxx_xyyz_0[j] + 1.5 * fl1_fx * tdy_xx_xyyz_0[j] + 0.5 * fl1_fx * tdy_xxx_yyz_0[j];

            tdz_xxxx_xyyz_0[j] = pa_x[j] * tdz_xxx_xyyz_0[j] + 1.5 * fl1_fx * tdz_xx_xyyz_0[j] + 0.5 * fl1_fx * tdz_xxx_yyz_0[j];

            tdx_xxxx_xyzz_0[j] =
                pa_x[j] * tdx_xxx_xyzz_0[j] + 1.5 * fl1_fx * tdx_xx_xyzz_0[j] + 0.5 * fl1_fx * tdx_xxx_yzz_0[j] + 0.5 * fl1_fx * ts_xxx_xyzz_0[j];

            tdy_xxxx_xyzz_0[j] = pa_x[j] * tdy_xxx_xyzz_0[j] + 1.5 * fl1_fx * tdy_xx_xyzz_0[j] + 0.5 * fl1_fx * tdy_xxx_yzz_0[j];

            tdz_xxxx_xyzz_0[j] = pa_x[j] * tdz_xxx_xyzz_0[j] + 1.5 * fl1_fx * tdz_xx_xyzz_0[j] + 0.5 * fl1_fx * tdz_xxx_yzz_0[j];

            tdx_xxxx_xzzz_0[j] =
                pa_x[j] * tdx_xxx_xzzz_0[j] + 1.5 * fl1_fx * tdx_xx_xzzz_0[j] + 0.5 * fl1_fx * tdx_xxx_zzz_0[j] + 0.5 * fl1_fx * ts_xxx_xzzz_0[j];

            tdy_xxxx_xzzz_0[j] = pa_x[j] * tdy_xxx_xzzz_0[j] + 1.5 * fl1_fx * tdy_xx_xzzz_0[j] + 0.5 * fl1_fx * tdy_xxx_zzz_0[j];

            tdz_xxxx_xzzz_0[j] = pa_x[j] * tdz_xxx_xzzz_0[j] + 1.5 * fl1_fx * tdz_xx_xzzz_0[j] + 0.5 * fl1_fx * tdz_xxx_zzz_0[j];

            tdx_xxxx_yyyy_0[j] = pa_x[j] * tdx_xxx_yyyy_0[j] + 1.5 * fl1_fx * tdx_xx_yyyy_0[j] + 0.5 * fl1_fx * ts_xxx_yyyy_0[j];

            tdy_xxxx_yyyy_0[j] = pa_x[j] * tdy_xxx_yyyy_0[j] + 1.5 * fl1_fx * tdy_xx_yyyy_0[j];

            tdz_xxxx_yyyy_0[j] = pa_x[j] * tdz_xxx_yyyy_0[j] + 1.5 * fl1_fx * tdz_xx_yyyy_0[j];

            tdx_xxxx_yyyz_0[j] = pa_x[j] * tdx_xxx_yyyz_0[j] + 1.5 * fl1_fx * tdx_xx_yyyz_0[j] + 0.5 * fl1_fx * ts_xxx_yyyz_0[j];

            tdy_xxxx_yyyz_0[j] = pa_x[j] * tdy_xxx_yyyz_0[j] + 1.5 * fl1_fx * tdy_xx_yyyz_0[j];

            tdz_xxxx_yyyz_0[j] = pa_x[j] * tdz_xxx_yyyz_0[j] + 1.5 * fl1_fx * tdz_xx_yyyz_0[j];

            tdx_xxxx_yyzz_0[j] = pa_x[j] * tdx_xxx_yyzz_0[j] + 1.5 * fl1_fx * tdx_xx_yyzz_0[j] + 0.5 * fl1_fx * ts_xxx_yyzz_0[j];

            tdy_xxxx_yyzz_0[j] = pa_x[j] * tdy_xxx_yyzz_0[j] + 1.5 * fl1_fx * tdy_xx_yyzz_0[j];

            tdz_xxxx_yyzz_0[j] = pa_x[j] * tdz_xxx_yyzz_0[j] + 1.5 * fl1_fx * tdz_xx_yyzz_0[j];

            tdx_xxxx_yzzz_0[j] = pa_x[j] * tdx_xxx_yzzz_0[j] + 1.5 * fl1_fx * tdx_xx_yzzz_0[j] + 0.5 * fl1_fx * ts_xxx_yzzz_0[j];

            tdy_xxxx_yzzz_0[j] = pa_x[j] * tdy_xxx_yzzz_0[j] + 1.5 * fl1_fx * tdy_xx_yzzz_0[j];

            tdz_xxxx_yzzz_0[j] = pa_x[j] * tdz_xxx_yzzz_0[j] + 1.5 * fl1_fx * tdz_xx_yzzz_0[j];

            tdx_xxxx_zzzz_0[j] = pa_x[j] * tdx_xxx_zzzz_0[j] + 1.5 * fl1_fx * tdx_xx_zzzz_0[j] + 0.5 * fl1_fx * ts_xxx_zzzz_0[j];

            tdy_xxxx_zzzz_0[j] = pa_x[j] * tdy_xxx_zzzz_0[j] + 1.5 * fl1_fx * tdy_xx_zzzz_0[j];

            tdz_xxxx_zzzz_0[j] = pa_x[j] * tdz_xxx_zzzz_0[j] + 1.5 * fl1_fx * tdz_xx_zzzz_0[j];

            tdx_xxxy_xxxx_0[j] =
                pa_x[j] * tdx_xxy_xxxx_0[j] + fl1_fx * tdx_xy_xxxx_0[j] + 2.0 * fl1_fx * tdx_xxy_xxx_0[j] + 0.5 * fl1_fx * ts_xxy_xxxx_0[j];

            tdy_xxxy_xxxx_0[j] = pa_x[j] * tdy_xxy_xxxx_0[j] + fl1_fx * tdy_xy_xxxx_0[j] + 2.0 * fl1_fx * tdy_xxy_xxx_0[j];

            tdz_xxxy_xxxx_0[j] = pa_x[j] * tdz_xxy_xxxx_0[j] + fl1_fx * tdz_xy_xxxx_0[j] + 2.0 * fl1_fx * tdz_xxy_xxx_0[j];

            tdx_xxxy_xxxy_0[j] =
                pa_x[j] * tdx_xxy_xxxy_0[j] + fl1_fx * tdx_xy_xxxy_0[j] + 1.5 * fl1_fx * tdx_xxy_xxy_0[j] + 0.5 * fl1_fx * ts_xxy_xxxy_0[j];
        }

        idx++;
    }
}

void
compElectricDipoleForGG_49_98(CMemBlock2D<double>&       primBuffer,
                              const CRecursionMap&       recursionMap,
                              const CMemBlock2D<double>& osFactors,
                              const int32_t              nOSFactors,
                              const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_d_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_d_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

        auto tdy_xxy_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 16);

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

        auto tdy_xy_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 16);

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

        auto tdy_xxy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 11);

        auto tdz_xxy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 11);

        auto tdx_xxy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 12);

        auto tdy_xxy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 12);

        auto tdz_xxy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 12);

        auto tdx_xxy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 13);

        auto tdy_xxy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 13);

        auto tdz_xxy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 13);

        auto tdx_xxy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 14);

        auto tdy_xxy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 14);

        auto tdz_xxy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 14);

        auto tdx_xxy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 15);

        auto tdy_xxy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 15);

        auto tdz_xxy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 15);

        auto tdx_xxy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 16);

        auto tdy_xxy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 16);

        auto tdz_xxy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 16);

        auto tdx_xxy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 17);

        auto tdy_xxy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 17);

        auto tdz_xxy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 17);

        auto tdx_xxy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 18);

        auto tdy_xxy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 18);

        auto tdz_xxy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 18);

        auto tdx_xxy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 19);

        auto tdy_xxy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 19);

        auto tdz_xxy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 19);

        auto tdx_xxz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 20);

        auto tdy_xxz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 20);

        auto tdz_xxz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 20);

        auto tdx_xxz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 21);

        auto tdy_xxz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 21);

        auto tdz_xxz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 21);

        auto tdx_xxz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 22);

        auto tdy_xxz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 22);

        auto ts_xxy_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 17);

        auto ts_xxy_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 18);

        auto ts_xxy_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 19);

        auto ts_xxy_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 20);

        auto ts_xxy_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 21);

        auto ts_xxy_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 22);

        auto ts_xxy_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 23);

        auto ts_xxy_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 24);

        auto ts_xxy_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 25);

        auto ts_xxy_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 26);

        auto ts_xxy_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 27);

        auto ts_xxy_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 28);

        auto ts_xxy_zzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 29);

        auto ts_xxz_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 30);

        auto ts_xxz_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 31);

        auto ts_xxz_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 32);

        // set up pointers to integrals

        auto tdy_xxxy_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 16);

        auto tdz_xxxy_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 16);

        auto tdx_xxxy_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 17);

        auto tdy_xxxy_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 17);

        auto tdz_xxxy_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 17);

        auto tdx_xxxy_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 18);

        auto tdy_xxxy_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 18);

        auto tdz_xxxy_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 18);

        auto tdx_xxxy_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 19);

        auto tdy_xxxy_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 19);

        auto tdz_xxxy_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 19);

        auto tdx_xxxy_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 20);

        auto tdy_xxxy_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 20);

        auto tdz_xxxy_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 20);

        auto tdx_xxxy_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 21);

        auto tdy_xxxy_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 21);

        auto tdz_xxxy_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 21);

        auto tdx_xxxy_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 22);

        auto tdy_xxxy_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 22);

        auto tdz_xxxy_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 22);

        auto tdx_xxxy_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 23);

        auto tdy_xxxy_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 23);

        auto tdz_xxxy_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 23);

        auto tdx_xxxy_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 24);

        auto tdy_xxxy_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 24);

        auto tdz_xxxy_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 24);

        auto tdx_xxxy_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 25);

        auto tdy_xxxy_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 25);

        auto tdz_xxxy_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 25);

        auto tdx_xxxy_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 26);

        auto tdy_xxxy_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 26);

        auto tdz_xxxy_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 26);

        auto tdx_xxxy_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 27);

        auto tdy_xxxy_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 27);

        auto tdz_xxxy_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 27);

        auto tdx_xxxy_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 28);

        auto tdy_xxxy_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 28);

        auto tdz_xxxy_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 28);

        auto tdx_xxxy_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 29);

        auto tdy_xxxy_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 29);

        auto tdz_xxxy_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 29);

        auto tdx_xxxz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 30);

        auto tdy_xxxz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 30);

        auto tdz_xxxz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 30);

        auto tdx_xxxz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 31);

        auto tdy_xxxz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 31);

        auto tdz_xxxz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 31);

        auto tdx_xxxz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 32);

        auto tdy_xxxz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 32);

        // Batch of Integrals (49,98)

        #pragma omp simd aligned(fx, pa_x, tdx_xxxy_xxxz_0, tdx_xxxy_xxyy_0, tdx_xxxy_xxyz_0, \
                                     tdx_xxxy_xxzz_0, tdx_xxxy_xyyy_0, tdx_xxxy_xyyz_0, tdx_xxxy_xyzz_0, tdx_xxxy_xzzz_0, \
                                     tdx_xxxy_yyyy_0, tdx_xxxy_yyyz_0, tdx_xxxy_yyzz_0, tdx_xxxy_yzzz_0, tdx_xxxy_zzzz_0, \
                                     tdx_xxxz_xxxx_0, tdx_xxxz_xxxy_0, tdx_xxxz_xxxz_0, tdx_xxy_xxxz_0, tdx_xxy_xxyy_0, \
                                     tdx_xxy_xxyz_0, tdx_xxy_xxz_0, tdx_xxy_xxzz_0, tdx_xxy_xyy_0, tdx_xxy_xyyy_0, \
                                     tdx_xxy_xyyz_0, tdx_xxy_xyz_0, tdx_xxy_xyzz_0, tdx_xxy_xzz_0, tdx_xxy_xzzz_0, \
                                     tdx_xxy_yyy_0, tdx_xxy_yyyy_0, tdx_xxy_yyyz_0, tdx_xxy_yyz_0, tdx_xxy_yyzz_0, \
                                     tdx_xxy_yzz_0, tdx_xxy_yzzz_0, tdx_xxy_zzz_0, tdx_xxy_zzzz_0, tdx_xxz_xxx_0, \
                                     tdx_xxz_xxxx_0, tdx_xxz_xxxy_0, tdx_xxz_xxxz_0, tdx_xxz_xxy_0, tdx_xxz_xxz_0, \
                                     tdx_xy_xxxz_0, tdx_xy_xxyy_0, tdx_xy_xxyz_0, tdx_xy_xxzz_0, tdx_xy_xyyy_0, \
                                     tdx_xy_xyyz_0, tdx_xy_xyzz_0, tdx_xy_xzzz_0, tdx_xy_yyyy_0, tdx_xy_yyyz_0, \
                                     tdx_xy_yyzz_0, tdx_xy_yzzz_0, tdx_xy_zzzz_0, tdx_xz_xxxx_0, tdx_xz_xxxy_0, \
                                     tdx_xz_xxxz_0, tdy_xxxy_xxxy_0, tdy_xxxy_xxxz_0, tdy_xxxy_xxyy_0, tdy_xxxy_xxyz_0, \
                                     tdy_xxxy_xxzz_0, tdy_xxxy_xyyy_0, tdy_xxxy_xyyz_0, tdy_xxxy_xyzz_0, tdy_xxxy_xzzz_0, \
                                     tdy_xxxy_yyyy_0, tdy_xxxy_yyyz_0, tdy_xxxy_yyzz_0, tdy_xxxy_yzzz_0, tdy_xxxy_zzzz_0, \
                                     tdy_xxxz_xxxx_0, tdy_xxxz_xxxy_0, tdy_xxxz_xxxz_0, tdy_xxy_xxxy_0, tdy_xxy_xxxz_0, \
                                     tdy_xxy_xxy_0, tdy_xxy_xxyy_0, tdy_xxy_xxyz_0, tdy_xxy_xxz_0, tdy_xxy_xxzz_0, \
                                     tdy_xxy_xyy_0, tdy_xxy_xyyy_0, tdy_xxy_xyyz_0, tdy_xxy_xyz_0, tdy_xxy_xyzz_0, \
                                     tdy_xxy_xzz_0, tdy_xxy_xzzz_0, tdy_xxy_yyy_0, tdy_xxy_yyyy_0, tdy_xxy_yyyz_0, \
                                     tdy_xxy_yyz_0, tdy_xxy_yyzz_0, tdy_xxy_yzz_0, tdy_xxy_yzzz_0, tdy_xxy_zzz_0, \
                                     tdy_xxy_zzzz_0, tdy_xxz_xxx_0, tdy_xxz_xxxx_0, tdy_xxz_xxxy_0, tdy_xxz_xxxz_0, \
                                     tdy_xxz_xxy_0, tdy_xxz_xxz_0, tdy_xy_xxxy_0, tdy_xy_xxxz_0, tdy_xy_xxyy_0, \
                                     tdy_xy_xxyz_0, tdy_xy_xxzz_0, tdy_xy_xyyy_0, tdy_xy_xyyz_0, tdy_xy_xyzz_0, \
                                     tdy_xy_xzzz_0, tdy_xy_yyyy_0, tdy_xy_yyyz_0, tdy_xy_yyzz_0, tdy_xy_yzzz_0, \
                                     tdy_xy_zzzz_0, tdy_xz_xxxx_0, tdy_xz_xxxy_0, tdy_xz_xxxz_0, tdz_xxxy_xxxy_0, \
                                     tdz_xxxy_xxxz_0, tdz_xxxy_xxyy_0, tdz_xxxy_xxyz_0, tdz_xxxy_xxzz_0, tdz_xxxy_xyyy_0, \
                                     tdz_xxxy_xyyz_0, tdz_xxxy_xyzz_0, tdz_xxxy_xzzz_0, tdz_xxxy_yyyy_0, tdz_xxxy_yyyz_0, \
                                     tdz_xxxy_yyzz_0, tdz_xxxy_yzzz_0, tdz_xxxy_zzzz_0, tdz_xxxz_xxxx_0, tdz_xxxz_xxxy_0, \
                                     tdz_xxy_xxxy_0, tdz_xxy_xxxz_0, tdz_xxy_xxy_0, tdz_xxy_xxyy_0, tdz_xxy_xxyz_0, \
                                     tdz_xxy_xxz_0, tdz_xxy_xxzz_0, tdz_xxy_xyy_0, tdz_xxy_xyyy_0, tdz_xxy_xyyz_0, \
                                     tdz_xxy_xyz_0, tdz_xxy_xyzz_0, tdz_xxy_xzz_0, tdz_xxy_xzzz_0, tdz_xxy_yyy_0, \
                                     tdz_xxy_yyyy_0, tdz_xxy_yyyz_0, tdz_xxy_yyz_0, tdz_xxy_yyzz_0, tdz_xxy_yzz_0, \
                                     tdz_xxy_yzzz_0, tdz_xxy_zzz_0, tdz_xxy_zzzz_0, tdz_xxz_xxx_0, tdz_xxz_xxxx_0, \
                                     tdz_xxz_xxxy_0, tdz_xxz_xxy_0, tdz_xy_xxxy_0, tdz_xy_xxxz_0, tdz_xy_xxyy_0, \
                                     tdz_xy_xxyz_0, tdz_xy_xxzz_0, tdz_xy_xyyy_0, tdz_xy_xyyz_0, tdz_xy_xyzz_0, \
                                     tdz_xy_xzzz_0, tdz_xy_yyyy_0, tdz_xy_yyyz_0, tdz_xy_yyzz_0, tdz_xy_yzzz_0, \
                                     tdz_xy_zzzz_0, tdz_xz_xxxx_0, tdz_xz_xxxy_0, ts_xxy_xxxz_0, ts_xxy_xxyy_0, \
                                     ts_xxy_xxyz_0, ts_xxy_xxzz_0, ts_xxy_xyyy_0, ts_xxy_xyyz_0, ts_xxy_xyzz_0, \
                                     ts_xxy_xzzz_0, ts_xxy_yyyy_0, ts_xxy_yyyz_0, ts_xxy_yyzz_0, ts_xxy_yzzz_0, \
                                     ts_xxy_zzzz_0, ts_xxz_xxxx_0, ts_xxz_xxxy_0, ts_xxz_xxxz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            tdy_xxxy_xxxy_0[j] = pa_x[j] * tdy_xxy_xxxy_0[j] + fl1_fx * tdy_xy_xxxy_0[j] + 1.5 * fl1_fx * tdy_xxy_xxy_0[j];

            tdz_xxxy_xxxy_0[j] = pa_x[j] * tdz_xxy_xxxy_0[j] + fl1_fx * tdz_xy_xxxy_0[j] + 1.5 * fl1_fx * tdz_xxy_xxy_0[j];

            tdx_xxxy_xxxz_0[j] =
                pa_x[j] * tdx_xxy_xxxz_0[j] + fl1_fx * tdx_xy_xxxz_0[j] + 1.5 * fl1_fx * tdx_xxy_xxz_0[j] + 0.5 * fl1_fx * ts_xxy_xxxz_0[j];

            tdy_xxxy_xxxz_0[j] = pa_x[j] * tdy_xxy_xxxz_0[j] + fl1_fx * tdy_xy_xxxz_0[j] + 1.5 * fl1_fx * tdy_xxy_xxz_0[j];

            tdz_xxxy_xxxz_0[j] = pa_x[j] * tdz_xxy_xxxz_0[j] + fl1_fx * tdz_xy_xxxz_0[j] + 1.5 * fl1_fx * tdz_xxy_xxz_0[j];

            tdx_xxxy_xxyy_0[j] =
                pa_x[j] * tdx_xxy_xxyy_0[j] + fl1_fx * tdx_xy_xxyy_0[j] + fl1_fx * tdx_xxy_xyy_0[j] + 0.5 * fl1_fx * ts_xxy_xxyy_0[j];

            tdy_xxxy_xxyy_0[j] = pa_x[j] * tdy_xxy_xxyy_0[j] + fl1_fx * tdy_xy_xxyy_0[j] + fl1_fx * tdy_xxy_xyy_0[j];

            tdz_xxxy_xxyy_0[j] = pa_x[j] * tdz_xxy_xxyy_0[j] + fl1_fx * tdz_xy_xxyy_0[j] + fl1_fx * tdz_xxy_xyy_0[j];

            tdx_xxxy_xxyz_0[j] =
                pa_x[j] * tdx_xxy_xxyz_0[j] + fl1_fx * tdx_xy_xxyz_0[j] + fl1_fx * tdx_xxy_xyz_0[j] + 0.5 * fl1_fx * ts_xxy_xxyz_0[j];

            tdy_xxxy_xxyz_0[j] = pa_x[j] * tdy_xxy_xxyz_0[j] + fl1_fx * tdy_xy_xxyz_0[j] + fl1_fx * tdy_xxy_xyz_0[j];

            tdz_xxxy_xxyz_0[j] = pa_x[j] * tdz_xxy_xxyz_0[j] + fl1_fx * tdz_xy_xxyz_0[j] + fl1_fx * tdz_xxy_xyz_0[j];

            tdx_xxxy_xxzz_0[j] =
                pa_x[j] * tdx_xxy_xxzz_0[j] + fl1_fx * tdx_xy_xxzz_0[j] + fl1_fx * tdx_xxy_xzz_0[j] + 0.5 * fl1_fx * ts_xxy_xxzz_0[j];

            tdy_xxxy_xxzz_0[j] = pa_x[j] * tdy_xxy_xxzz_0[j] + fl1_fx * tdy_xy_xxzz_0[j] + fl1_fx * tdy_xxy_xzz_0[j];

            tdz_xxxy_xxzz_0[j] = pa_x[j] * tdz_xxy_xxzz_0[j] + fl1_fx * tdz_xy_xxzz_0[j] + fl1_fx * tdz_xxy_xzz_0[j];

            tdx_xxxy_xyyy_0[j] =
                pa_x[j] * tdx_xxy_xyyy_0[j] + fl1_fx * tdx_xy_xyyy_0[j] + 0.5 * fl1_fx * tdx_xxy_yyy_0[j] + 0.5 * fl1_fx * ts_xxy_xyyy_0[j];

            tdy_xxxy_xyyy_0[j] = pa_x[j] * tdy_xxy_xyyy_0[j] + fl1_fx * tdy_xy_xyyy_0[j] + 0.5 * fl1_fx * tdy_xxy_yyy_0[j];

            tdz_xxxy_xyyy_0[j] = pa_x[j] * tdz_xxy_xyyy_0[j] + fl1_fx * tdz_xy_xyyy_0[j] + 0.5 * fl1_fx * tdz_xxy_yyy_0[j];

            tdx_xxxy_xyyz_0[j] =
                pa_x[j] * tdx_xxy_xyyz_0[j] + fl1_fx * tdx_xy_xyyz_0[j] + 0.5 * fl1_fx * tdx_xxy_yyz_0[j] + 0.5 * fl1_fx * ts_xxy_xyyz_0[j];

            tdy_xxxy_xyyz_0[j] = pa_x[j] * tdy_xxy_xyyz_0[j] + fl1_fx * tdy_xy_xyyz_0[j] + 0.5 * fl1_fx * tdy_xxy_yyz_0[j];

            tdz_xxxy_xyyz_0[j] = pa_x[j] * tdz_xxy_xyyz_0[j] + fl1_fx * tdz_xy_xyyz_0[j] + 0.5 * fl1_fx * tdz_xxy_yyz_0[j];

            tdx_xxxy_xyzz_0[j] =
                pa_x[j] * tdx_xxy_xyzz_0[j] + fl1_fx * tdx_xy_xyzz_0[j] + 0.5 * fl1_fx * tdx_xxy_yzz_0[j] + 0.5 * fl1_fx * ts_xxy_xyzz_0[j];

            tdy_xxxy_xyzz_0[j] = pa_x[j] * tdy_xxy_xyzz_0[j] + fl1_fx * tdy_xy_xyzz_0[j] + 0.5 * fl1_fx * tdy_xxy_yzz_0[j];

            tdz_xxxy_xyzz_0[j] = pa_x[j] * tdz_xxy_xyzz_0[j] + fl1_fx * tdz_xy_xyzz_0[j] + 0.5 * fl1_fx * tdz_xxy_yzz_0[j];

            tdx_xxxy_xzzz_0[j] =
                pa_x[j] * tdx_xxy_xzzz_0[j] + fl1_fx * tdx_xy_xzzz_0[j] + 0.5 * fl1_fx * tdx_xxy_zzz_0[j] + 0.5 * fl1_fx * ts_xxy_xzzz_0[j];

            tdy_xxxy_xzzz_0[j] = pa_x[j] * tdy_xxy_xzzz_0[j] + fl1_fx * tdy_xy_xzzz_0[j] + 0.5 * fl1_fx * tdy_xxy_zzz_0[j];

            tdz_xxxy_xzzz_0[j] = pa_x[j] * tdz_xxy_xzzz_0[j] + fl1_fx * tdz_xy_xzzz_0[j] + 0.5 * fl1_fx * tdz_xxy_zzz_0[j];

            tdx_xxxy_yyyy_0[j] = pa_x[j] * tdx_xxy_yyyy_0[j] + fl1_fx * tdx_xy_yyyy_0[j] + 0.5 * fl1_fx * ts_xxy_yyyy_0[j];

            tdy_xxxy_yyyy_0[j] = pa_x[j] * tdy_xxy_yyyy_0[j] + fl1_fx * tdy_xy_yyyy_0[j];

            tdz_xxxy_yyyy_0[j] = pa_x[j] * tdz_xxy_yyyy_0[j] + fl1_fx * tdz_xy_yyyy_0[j];

            tdx_xxxy_yyyz_0[j] = pa_x[j] * tdx_xxy_yyyz_0[j] + fl1_fx * tdx_xy_yyyz_0[j] + 0.5 * fl1_fx * ts_xxy_yyyz_0[j];

            tdy_xxxy_yyyz_0[j] = pa_x[j] * tdy_xxy_yyyz_0[j] + fl1_fx * tdy_xy_yyyz_0[j];

            tdz_xxxy_yyyz_0[j] = pa_x[j] * tdz_xxy_yyyz_0[j] + fl1_fx * tdz_xy_yyyz_0[j];

            tdx_xxxy_yyzz_0[j] = pa_x[j] * tdx_xxy_yyzz_0[j] + fl1_fx * tdx_xy_yyzz_0[j] + 0.5 * fl1_fx * ts_xxy_yyzz_0[j];

            tdy_xxxy_yyzz_0[j] = pa_x[j] * tdy_xxy_yyzz_0[j] + fl1_fx * tdy_xy_yyzz_0[j];

            tdz_xxxy_yyzz_0[j] = pa_x[j] * tdz_xxy_yyzz_0[j] + fl1_fx * tdz_xy_yyzz_0[j];

            tdx_xxxy_yzzz_0[j] = pa_x[j] * tdx_xxy_yzzz_0[j] + fl1_fx * tdx_xy_yzzz_0[j] + 0.5 * fl1_fx * ts_xxy_yzzz_0[j];

            tdy_xxxy_yzzz_0[j] = pa_x[j] * tdy_xxy_yzzz_0[j] + fl1_fx * tdy_xy_yzzz_0[j];

            tdz_xxxy_yzzz_0[j] = pa_x[j] * tdz_xxy_yzzz_0[j] + fl1_fx * tdz_xy_yzzz_0[j];

            tdx_xxxy_zzzz_0[j] = pa_x[j] * tdx_xxy_zzzz_0[j] + fl1_fx * tdx_xy_zzzz_0[j] + 0.5 * fl1_fx * ts_xxy_zzzz_0[j];

            tdy_xxxy_zzzz_0[j] = pa_x[j] * tdy_xxy_zzzz_0[j] + fl1_fx * tdy_xy_zzzz_0[j];

            tdz_xxxy_zzzz_0[j] = pa_x[j] * tdz_xxy_zzzz_0[j] + fl1_fx * tdz_xy_zzzz_0[j];

            tdx_xxxz_xxxx_0[j] =
                pa_x[j] * tdx_xxz_xxxx_0[j] + fl1_fx * tdx_xz_xxxx_0[j] + 2.0 * fl1_fx * tdx_xxz_xxx_0[j] + 0.5 * fl1_fx * ts_xxz_xxxx_0[j];

            tdy_xxxz_xxxx_0[j] = pa_x[j] * tdy_xxz_xxxx_0[j] + fl1_fx * tdy_xz_xxxx_0[j] + 2.0 * fl1_fx * tdy_xxz_xxx_0[j];

            tdz_xxxz_xxxx_0[j] = pa_x[j] * tdz_xxz_xxxx_0[j] + fl1_fx * tdz_xz_xxxx_0[j] + 2.0 * fl1_fx * tdz_xxz_xxx_0[j];

            tdx_xxxz_xxxy_0[j] =
                pa_x[j] * tdx_xxz_xxxy_0[j] + fl1_fx * tdx_xz_xxxy_0[j] + 1.5 * fl1_fx * tdx_xxz_xxy_0[j] + 0.5 * fl1_fx * ts_xxz_xxxy_0[j];

            tdy_xxxz_xxxy_0[j] = pa_x[j] * tdy_xxz_xxxy_0[j] + fl1_fx * tdy_xz_xxxy_0[j] + 1.5 * fl1_fx * tdy_xxz_xxy_0[j];

            tdz_xxxz_xxxy_0[j] = pa_x[j] * tdz_xxz_xxxy_0[j] + fl1_fx * tdz_xz_xxxy_0[j] + 1.5 * fl1_fx * tdz_xxz_xxy_0[j];

            tdx_xxxz_xxxz_0[j] =
                pa_x[j] * tdx_xxz_xxxz_0[j] + fl1_fx * tdx_xz_xxxz_0[j] + 1.5 * fl1_fx * tdx_xxz_xxz_0[j] + 0.5 * fl1_fx * ts_xxz_xxxz_0[j];

            tdy_xxxz_xxxz_0[j] = pa_x[j] * tdy_xxz_xxxz_0[j] + fl1_fx * tdy_xz_xxxz_0[j] + 1.5 * fl1_fx * tdy_xxz_xxz_0[j];
        }

        idx++;
    }
}

void
compElectricDipoleForGG_98_147(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const int32_t              nOSFactors,
                               const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_d_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_d_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

        auto tdz_xxz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 32);

        auto tdx_xxz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 33);

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

        auto tdz_xz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 32);

        auto tdx_xz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 33);

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

        auto tdz_xxz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 22);

        auto tdx_xxz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 23);

        auto tdy_xxz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 23);

        auto tdz_xxz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 23);

        auto tdx_xxz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 24);

        auto tdy_xxz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 24);

        auto tdz_xxz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 24);

        auto tdx_xxz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 25);

        auto tdy_xxz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 25);

        auto tdz_xxz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 25);

        auto tdx_xxz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 26);

        auto tdy_xxz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 26);

        auto tdz_xxz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 26);

        auto tdx_xxz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 27);

        auto tdy_xxz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 27);

        auto tdz_xxz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 27);

        auto tdx_xxz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 28);

        auto tdy_xxz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 28);

        auto tdz_xxz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 28);

        auto tdx_xxz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 29);

        auto tdy_xxz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 29);

        auto tdz_xxz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 29);

        auto tdx_xyy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 30);

        auto tdy_xyy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 30);

        auto tdz_xyy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 30);

        auto tdx_xyy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 31);

        auto tdy_xyy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 31);

        auto tdz_xyy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 31);

        auto tdx_xyy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 32);

        auto tdy_xyy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 32);

        auto tdz_xyy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 32);

        auto tdx_xyy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 33);

        auto tdy_xyy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 33);

        auto tdz_xyy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 33);

        auto ts_xxz_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 33);

        auto ts_xxz_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 34);

        auto ts_xxz_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 35);

        auto ts_xxz_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 36);

        auto ts_xxz_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 37);

        auto ts_xxz_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 38);

        auto ts_xxz_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 39);

        auto ts_xxz_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 40);

        auto ts_xxz_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 41);

        auto ts_xxz_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 42);

        auto ts_xxz_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 43);

        auto ts_xxz_zzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 44);

        auto ts_xyy_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 45);

        auto ts_xyy_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 46);

        auto ts_xyy_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 47);

        auto ts_xyy_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 48);

        // set up pointers to integrals

        auto tdz_xxxz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 32);

        auto tdx_xxxz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 33);

        auto tdy_xxxz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 33);

        auto tdz_xxxz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 33);

        auto tdx_xxxz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 34);

        auto tdy_xxxz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 34);

        auto tdz_xxxz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 34);

        auto tdx_xxxz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 35);

        auto tdy_xxxz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 35);

        auto tdz_xxxz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 35);

        auto tdx_xxxz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 36);

        auto tdy_xxxz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 36);

        auto tdz_xxxz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 36);

        auto tdx_xxxz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 37);

        auto tdy_xxxz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 37);

        auto tdz_xxxz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 37);

        auto tdx_xxxz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 38);

        auto tdy_xxxz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 38);

        auto tdz_xxxz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 38);

        auto tdx_xxxz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 39);

        auto tdy_xxxz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 39);

        auto tdz_xxxz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 39);

        auto tdx_xxxz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 40);

        auto tdy_xxxz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 40);

        auto tdz_xxxz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 40);

        auto tdx_xxxz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 41);

        auto tdy_xxxz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 41);

        auto tdz_xxxz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 41);

        auto tdx_xxxz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 42);

        auto tdy_xxxz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 42);

        auto tdz_xxxz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 42);

        auto tdx_xxxz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 43);

        auto tdy_xxxz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 43);

        auto tdz_xxxz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 43);

        auto tdx_xxxz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 44);

        auto tdy_xxxz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 44);

        auto tdz_xxxz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 44);

        auto tdx_xxyy_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 45);

        auto tdy_xxyy_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 45);

        auto tdz_xxyy_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 45);

        auto tdx_xxyy_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 46);

        auto tdy_xxyy_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 46);

        auto tdz_xxyy_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 46);

        auto tdx_xxyy_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 47);

        auto tdy_xxyy_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 47);

        auto tdz_xxyy_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 47);

        auto tdx_xxyy_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 48);

        auto tdy_xxyy_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 48);

        auto tdz_xxyy_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 48);

        // Batch of Integrals (98,147)

        #pragma omp simd aligned(fx, pa_x, tdx_xxxz_xxyy_0, tdx_xxxz_xxyz_0, tdx_xxxz_xxzz_0, \
                                     tdx_xxxz_xyyy_0, tdx_xxxz_xyyz_0, tdx_xxxz_xyzz_0, tdx_xxxz_xzzz_0, tdx_xxxz_yyyy_0, \
                                     tdx_xxxz_yyyz_0, tdx_xxxz_yyzz_0, tdx_xxxz_yzzz_0, tdx_xxxz_zzzz_0, tdx_xxyy_xxxx_0, \
                                     tdx_xxyy_xxxy_0, tdx_xxyy_xxxz_0, tdx_xxyy_xxyy_0, tdx_xxz_xxyy_0, tdx_xxz_xxyz_0, \
                                     tdx_xxz_xxzz_0, tdx_xxz_xyy_0, tdx_xxz_xyyy_0, tdx_xxz_xyyz_0, tdx_xxz_xyz_0, \
                                     tdx_xxz_xyzz_0, tdx_xxz_xzz_0, tdx_xxz_xzzz_0, tdx_xxz_yyy_0, tdx_xxz_yyyy_0, \
                                     tdx_xxz_yyyz_0, tdx_xxz_yyz_0, tdx_xxz_yyzz_0, tdx_xxz_yzz_0, tdx_xxz_yzzz_0, \
                                     tdx_xxz_zzz_0, tdx_xxz_zzzz_0, tdx_xyy_xxx_0, tdx_xyy_xxxx_0, tdx_xyy_xxxy_0, \
                                     tdx_xyy_xxxz_0, tdx_xyy_xxy_0, tdx_xyy_xxyy_0, tdx_xyy_xxz_0, tdx_xyy_xyy_0, \
                                     tdx_xz_xxyy_0, tdx_xz_xxyz_0, tdx_xz_xxzz_0, tdx_xz_xyyy_0, tdx_xz_xyyz_0, \
                                     tdx_xz_xyzz_0, tdx_xz_xzzz_0, tdx_xz_yyyy_0, tdx_xz_yyyz_0, tdx_xz_yyzz_0, \
                                     tdx_xz_yzzz_0, tdx_xz_zzzz_0, tdx_yy_xxxx_0, tdx_yy_xxxy_0, tdx_yy_xxxz_0, \
                                     tdx_yy_xxyy_0, tdy_xxxz_xxyy_0, tdy_xxxz_xxyz_0, tdy_xxxz_xxzz_0, tdy_xxxz_xyyy_0, \
                                     tdy_xxxz_xyyz_0, tdy_xxxz_xyzz_0, tdy_xxxz_xzzz_0, tdy_xxxz_yyyy_0, tdy_xxxz_yyyz_0, \
                                     tdy_xxxz_yyzz_0, tdy_xxxz_yzzz_0, tdy_xxxz_zzzz_0, tdy_xxyy_xxxx_0, tdy_xxyy_xxxy_0, \
                                     tdy_xxyy_xxxz_0, tdy_xxyy_xxyy_0, tdy_xxz_xxyy_0, tdy_xxz_xxyz_0, tdy_xxz_xxzz_0, \
                                     tdy_xxz_xyy_0, tdy_xxz_xyyy_0, tdy_xxz_xyyz_0, tdy_xxz_xyz_0, tdy_xxz_xyzz_0, \
                                     tdy_xxz_xzz_0, tdy_xxz_xzzz_0, tdy_xxz_yyy_0, tdy_xxz_yyyy_0, tdy_xxz_yyyz_0, \
                                     tdy_xxz_yyz_0, tdy_xxz_yyzz_0, tdy_xxz_yzz_0, tdy_xxz_yzzz_0, tdy_xxz_zzz_0, \
                                     tdy_xxz_zzzz_0, tdy_xyy_xxx_0, tdy_xyy_xxxx_0, tdy_xyy_xxxy_0, tdy_xyy_xxxz_0, \
                                     tdy_xyy_xxy_0, tdy_xyy_xxyy_0, tdy_xyy_xxz_0, tdy_xyy_xyy_0, tdy_xz_xxyy_0, \
                                     tdy_xz_xxyz_0, tdy_xz_xxzz_0, tdy_xz_xyyy_0, tdy_xz_xyyz_0, tdy_xz_xyzz_0, \
                                     tdy_xz_xzzz_0, tdy_xz_yyyy_0, tdy_xz_yyyz_0, tdy_xz_yyzz_0, tdy_xz_yzzz_0, \
                                     tdy_xz_zzzz_0, tdy_yy_xxxx_0, tdy_yy_xxxy_0, tdy_yy_xxxz_0, tdy_yy_xxyy_0, \
                                     tdz_xxxz_xxxz_0, tdz_xxxz_xxyy_0, tdz_xxxz_xxyz_0, tdz_xxxz_xxzz_0, tdz_xxxz_xyyy_0, \
                                     tdz_xxxz_xyyz_0, tdz_xxxz_xyzz_0, tdz_xxxz_xzzz_0, tdz_xxxz_yyyy_0, tdz_xxxz_yyyz_0, \
                                     tdz_xxxz_yyzz_0, tdz_xxxz_yzzz_0, tdz_xxxz_zzzz_0, tdz_xxyy_xxxx_0, tdz_xxyy_xxxy_0, \
                                     tdz_xxyy_xxxz_0, tdz_xxyy_xxyy_0, tdz_xxz_xxxz_0, tdz_xxz_xxyy_0, tdz_xxz_xxyz_0, \
                                     tdz_xxz_xxz_0, tdz_xxz_xxzz_0, tdz_xxz_xyy_0, tdz_xxz_xyyy_0, tdz_xxz_xyyz_0, \
                                     tdz_xxz_xyz_0, tdz_xxz_xyzz_0, tdz_xxz_xzz_0, tdz_xxz_xzzz_0, tdz_xxz_yyy_0, \
                                     tdz_xxz_yyyy_0, tdz_xxz_yyyz_0, tdz_xxz_yyz_0, tdz_xxz_yyzz_0, tdz_xxz_yzz_0, \
                                     tdz_xxz_yzzz_0, tdz_xxz_zzz_0, tdz_xxz_zzzz_0, tdz_xyy_xxx_0, tdz_xyy_xxxx_0, \
                                     tdz_xyy_xxxy_0, tdz_xyy_xxxz_0, tdz_xyy_xxy_0, tdz_xyy_xxyy_0, tdz_xyy_xxz_0, \
                                     tdz_xyy_xyy_0, tdz_xz_xxxz_0, tdz_xz_xxyy_0, tdz_xz_xxyz_0, tdz_xz_xxzz_0, \
                                     tdz_xz_xyyy_0, tdz_xz_xyyz_0, tdz_xz_xyzz_0, tdz_xz_xzzz_0, tdz_xz_yyyy_0, \
                                     tdz_xz_yyyz_0, tdz_xz_yyzz_0, tdz_xz_yzzz_0, tdz_xz_zzzz_0, tdz_yy_xxxx_0, \
                                     tdz_yy_xxxy_0, tdz_yy_xxxz_0, tdz_yy_xxyy_0, ts_xxz_xxyy_0, ts_xxz_xxyz_0, \
                                     ts_xxz_xxzz_0, ts_xxz_xyyy_0, ts_xxz_xyyz_0, ts_xxz_xyzz_0, ts_xxz_xzzz_0, \
                                     ts_xxz_yyyy_0, ts_xxz_yyyz_0, ts_xxz_yyzz_0, ts_xxz_yzzz_0, ts_xxz_zzzz_0, \
                                     ts_xyy_xxxx_0, ts_xyy_xxxy_0, ts_xyy_xxxz_0, ts_xyy_xxyy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            tdz_xxxz_xxxz_0[j] = pa_x[j] * tdz_xxz_xxxz_0[j] + fl1_fx * tdz_xz_xxxz_0[j] + 1.5 * fl1_fx * tdz_xxz_xxz_0[j];

            tdx_xxxz_xxyy_0[j] =
                pa_x[j] * tdx_xxz_xxyy_0[j] + fl1_fx * tdx_xz_xxyy_0[j] + fl1_fx * tdx_xxz_xyy_0[j] + 0.5 * fl1_fx * ts_xxz_xxyy_0[j];

            tdy_xxxz_xxyy_0[j] = pa_x[j] * tdy_xxz_xxyy_0[j] + fl1_fx * tdy_xz_xxyy_0[j] + fl1_fx * tdy_xxz_xyy_0[j];

            tdz_xxxz_xxyy_0[j] = pa_x[j] * tdz_xxz_xxyy_0[j] + fl1_fx * tdz_xz_xxyy_0[j] + fl1_fx * tdz_xxz_xyy_0[j];

            tdx_xxxz_xxyz_0[j] =
                pa_x[j] * tdx_xxz_xxyz_0[j] + fl1_fx * tdx_xz_xxyz_0[j] + fl1_fx * tdx_xxz_xyz_0[j] + 0.5 * fl1_fx * ts_xxz_xxyz_0[j];

            tdy_xxxz_xxyz_0[j] = pa_x[j] * tdy_xxz_xxyz_0[j] + fl1_fx * tdy_xz_xxyz_0[j] + fl1_fx * tdy_xxz_xyz_0[j];

            tdz_xxxz_xxyz_0[j] = pa_x[j] * tdz_xxz_xxyz_0[j] + fl1_fx * tdz_xz_xxyz_0[j] + fl1_fx * tdz_xxz_xyz_0[j];

            tdx_xxxz_xxzz_0[j] =
                pa_x[j] * tdx_xxz_xxzz_0[j] + fl1_fx * tdx_xz_xxzz_0[j] + fl1_fx * tdx_xxz_xzz_0[j] + 0.5 * fl1_fx * ts_xxz_xxzz_0[j];

            tdy_xxxz_xxzz_0[j] = pa_x[j] * tdy_xxz_xxzz_0[j] + fl1_fx * tdy_xz_xxzz_0[j] + fl1_fx * tdy_xxz_xzz_0[j];

            tdz_xxxz_xxzz_0[j] = pa_x[j] * tdz_xxz_xxzz_0[j] + fl1_fx * tdz_xz_xxzz_0[j] + fl1_fx * tdz_xxz_xzz_0[j];

            tdx_xxxz_xyyy_0[j] =
                pa_x[j] * tdx_xxz_xyyy_0[j] + fl1_fx * tdx_xz_xyyy_0[j] + 0.5 * fl1_fx * tdx_xxz_yyy_0[j] + 0.5 * fl1_fx * ts_xxz_xyyy_0[j];

            tdy_xxxz_xyyy_0[j] = pa_x[j] * tdy_xxz_xyyy_0[j] + fl1_fx * tdy_xz_xyyy_0[j] + 0.5 * fl1_fx * tdy_xxz_yyy_0[j];

            tdz_xxxz_xyyy_0[j] = pa_x[j] * tdz_xxz_xyyy_0[j] + fl1_fx * tdz_xz_xyyy_0[j] + 0.5 * fl1_fx * tdz_xxz_yyy_0[j];

            tdx_xxxz_xyyz_0[j] =
                pa_x[j] * tdx_xxz_xyyz_0[j] + fl1_fx * tdx_xz_xyyz_0[j] + 0.5 * fl1_fx * tdx_xxz_yyz_0[j] + 0.5 * fl1_fx * ts_xxz_xyyz_0[j];

            tdy_xxxz_xyyz_0[j] = pa_x[j] * tdy_xxz_xyyz_0[j] + fl1_fx * tdy_xz_xyyz_0[j] + 0.5 * fl1_fx * tdy_xxz_yyz_0[j];

            tdz_xxxz_xyyz_0[j] = pa_x[j] * tdz_xxz_xyyz_0[j] + fl1_fx * tdz_xz_xyyz_0[j] + 0.5 * fl1_fx * tdz_xxz_yyz_0[j];

            tdx_xxxz_xyzz_0[j] =
                pa_x[j] * tdx_xxz_xyzz_0[j] + fl1_fx * tdx_xz_xyzz_0[j] + 0.5 * fl1_fx * tdx_xxz_yzz_0[j] + 0.5 * fl1_fx * ts_xxz_xyzz_0[j];

            tdy_xxxz_xyzz_0[j] = pa_x[j] * tdy_xxz_xyzz_0[j] + fl1_fx * tdy_xz_xyzz_0[j] + 0.5 * fl1_fx * tdy_xxz_yzz_0[j];

            tdz_xxxz_xyzz_0[j] = pa_x[j] * tdz_xxz_xyzz_0[j] + fl1_fx * tdz_xz_xyzz_0[j] + 0.5 * fl1_fx * tdz_xxz_yzz_0[j];

            tdx_xxxz_xzzz_0[j] =
                pa_x[j] * tdx_xxz_xzzz_0[j] + fl1_fx * tdx_xz_xzzz_0[j] + 0.5 * fl1_fx * tdx_xxz_zzz_0[j] + 0.5 * fl1_fx * ts_xxz_xzzz_0[j];

            tdy_xxxz_xzzz_0[j] = pa_x[j] * tdy_xxz_xzzz_0[j] + fl1_fx * tdy_xz_xzzz_0[j] + 0.5 * fl1_fx * tdy_xxz_zzz_0[j];

            tdz_xxxz_xzzz_0[j] = pa_x[j] * tdz_xxz_xzzz_0[j] + fl1_fx * tdz_xz_xzzz_0[j] + 0.5 * fl1_fx * tdz_xxz_zzz_0[j];

            tdx_xxxz_yyyy_0[j] = pa_x[j] * tdx_xxz_yyyy_0[j] + fl1_fx * tdx_xz_yyyy_0[j] + 0.5 * fl1_fx * ts_xxz_yyyy_0[j];

            tdy_xxxz_yyyy_0[j] = pa_x[j] * tdy_xxz_yyyy_0[j] + fl1_fx * tdy_xz_yyyy_0[j];

            tdz_xxxz_yyyy_0[j] = pa_x[j] * tdz_xxz_yyyy_0[j] + fl1_fx * tdz_xz_yyyy_0[j];

            tdx_xxxz_yyyz_0[j] = pa_x[j] * tdx_xxz_yyyz_0[j] + fl1_fx * tdx_xz_yyyz_0[j] + 0.5 * fl1_fx * ts_xxz_yyyz_0[j];

            tdy_xxxz_yyyz_0[j] = pa_x[j] * tdy_xxz_yyyz_0[j] + fl1_fx * tdy_xz_yyyz_0[j];

            tdz_xxxz_yyyz_0[j] = pa_x[j] * tdz_xxz_yyyz_0[j] + fl1_fx * tdz_xz_yyyz_0[j];

            tdx_xxxz_yyzz_0[j] = pa_x[j] * tdx_xxz_yyzz_0[j] + fl1_fx * tdx_xz_yyzz_0[j] + 0.5 * fl1_fx * ts_xxz_yyzz_0[j];

            tdy_xxxz_yyzz_0[j] = pa_x[j] * tdy_xxz_yyzz_0[j] + fl1_fx * tdy_xz_yyzz_0[j];

            tdz_xxxz_yyzz_0[j] = pa_x[j] * tdz_xxz_yyzz_0[j] + fl1_fx * tdz_xz_yyzz_0[j];

            tdx_xxxz_yzzz_0[j] = pa_x[j] * tdx_xxz_yzzz_0[j] + fl1_fx * tdx_xz_yzzz_0[j] + 0.5 * fl1_fx * ts_xxz_yzzz_0[j];

            tdy_xxxz_yzzz_0[j] = pa_x[j] * tdy_xxz_yzzz_0[j] + fl1_fx * tdy_xz_yzzz_0[j];

            tdz_xxxz_yzzz_0[j] = pa_x[j] * tdz_xxz_yzzz_0[j] + fl1_fx * tdz_xz_yzzz_0[j];

            tdx_xxxz_zzzz_0[j] = pa_x[j] * tdx_xxz_zzzz_0[j] + fl1_fx * tdx_xz_zzzz_0[j] + 0.5 * fl1_fx * ts_xxz_zzzz_0[j];

            tdy_xxxz_zzzz_0[j] = pa_x[j] * tdy_xxz_zzzz_0[j] + fl1_fx * tdy_xz_zzzz_0[j];

            tdz_xxxz_zzzz_0[j] = pa_x[j] * tdz_xxz_zzzz_0[j] + fl1_fx * tdz_xz_zzzz_0[j];

            tdx_xxyy_xxxx_0[j] =
                pa_x[j] * tdx_xyy_xxxx_0[j] + 0.5 * fl1_fx * tdx_yy_xxxx_0[j] + 2.0 * fl1_fx * tdx_xyy_xxx_0[j] + 0.5 * fl1_fx * ts_xyy_xxxx_0[j];

            tdy_xxyy_xxxx_0[j] = pa_x[j] * tdy_xyy_xxxx_0[j] + 0.5 * fl1_fx * tdy_yy_xxxx_0[j] + 2.0 * fl1_fx * tdy_xyy_xxx_0[j];

            tdz_xxyy_xxxx_0[j] = pa_x[j] * tdz_xyy_xxxx_0[j] + 0.5 * fl1_fx * tdz_yy_xxxx_0[j] + 2.0 * fl1_fx * tdz_xyy_xxx_0[j];

            tdx_xxyy_xxxy_0[j] =
                pa_x[j] * tdx_xyy_xxxy_0[j] + 0.5 * fl1_fx * tdx_yy_xxxy_0[j] + 1.5 * fl1_fx * tdx_xyy_xxy_0[j] + 0.5 * fl1_fx * ts_xyy_xxxy_0[j];

            tdy_xxyy_xxxy_0[j] = pa_x[j] * tdy_xyy_xxxy_0[j] + 0.5 * fl1_fx * tdy_yy_xxxy_0[j] + 1.5 * fl1_fx * tdy_xyy_xxy_0[j];

            tdz_xxyy_xxxy_0[j] = pa_x[j] * tdz_xyy_xxxy_0[j] + 0.5 * fl1_fx * tdz_yy_xxxy_0[j] + 1.5 * fl1_fx * tdz_xyy_xxy_0[j];

            tdx_xxyy_xxxz_0[j] =
                pa_x[j] * tdx_xyy_xxxz_0[j] + 0.5 * fl1_fx * tdx_yy_xxxz_0[j] + 1.5 * fl1_fx * tdx_xyy_xxz_0[j] + 0.5 * fl1_fx * ts_xyy_xxxz_0[j];

            tdy_xxyy_xxxz_0[j] = pa_x[j] * tdy_xyy_xxxz_0[j] + 0.5 * fl1_fx * tdy_yy_xxxz_0[j] + 1.5 * fl1_fx * tdy_xyy_xxz_0[j];

            tdz_xxyy_xxxz_0[j] = pa_x[j] * tdz_xyy_xxxz_0[j] + 0.5 * fl1_fx * tdz_yy_xxxz_0[j] + 1.5 * fl1_fx * tdz_xyy_xxz_0[j];

            tdx_xxyy_xxyy_0[j] =
                pa_x[j] * tdx_xyy_xxyy_0[j] + 0.5 * fl1_fx * tdx_yy_xxyy_0[j] + fl1_fx * tdx_xyy_xyy_0[j] + 0.5 * fl1_fx * ts_xyy_xxyy_0[j];

            tdy_xxyy_xxyy_0[j] = pa_x[j] * tdy_xyy_xxyy_0[j] + 0.5 * fl1_fx * tdy_yy_xxyy_0[j] + fl1_fx * tdy_xyy_xyy_0[j];

            tdz_xxyy_xxyy_0[j] = pa_x[j] * tdz_xyy_xxyy_0[j] + 0.5 * fl1_fx * tdz_yy_xxyy_0[j] + fl1_fx * tdz_xyy_xyy_0[j];
        }

        idx++;
    }
}

void
compElectricDipoleForGG_147_195(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const int32_t              nOSFactors,
                                const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_d_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_d_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

        auto tdx_xyy_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 49);

        auto tdy_xyy_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 49);

        auto tdz_xyy_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 49);

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

        auto tdx_xyy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 34);

        auto tdy_xyy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 34);

        auto tdz_xyy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 34);

        auto tdx_xyy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 35);

        auto tdy_xyy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 35);

        auto tdz_xyy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 35);

        auto tdx_xyy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 36);

        auto tdy_xyy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 36);

        auto tdz_xyy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 36);

        auto tdx_xyy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 37);

        auto tdy_xyy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 37);

        auto tdz_xyy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 37);

        auto tdx_xyy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 38);

        auto tdy_xyy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 38);

        auto tdz_xyy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 38);

        auto tdx_xyy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 39);

        auto tdy_xyy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 39);

        auto tdz_xyy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 39);

        auto tdx_xyz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 40);

        auto tdy_xyz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 40);

        auto tdz_xyz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 40);

        auto tdx_xyz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 41);

        auto tdy_xyz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 41);

        auto tdz_xyz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 41);

        auto tdx_xyz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 42);

        auto tdy_xyz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 42);

        auto tdz_xyz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 42);

        auto tdx_xyz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 43);

        auto tdy_xyz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 43);

        auto tdz_xyz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 43);

        auto tdx_xyz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 44);

        auto tdy_xyz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 44);

        auto tdz_xyz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 44);

        auto ts_xyy_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 49);

        auto ts_xyy_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 50);

        auto ts_xyy_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 51);

        auto ts_xyy_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 52);

        auto ts_xyy_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 53);

        auto ts_xyy_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 54);

        auto ts_xyy_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 55);

        auto ts_xyy_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 56);

        auto ts_xyy_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 57);

        auto ts_xyy_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 58);

        auto ts_xyy_zzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 59);

        auto ts_xyz_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 60);

        auto ts_xyz_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 61);

        auto ts_xyz_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 62);

        auto ts_xyz_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 63);

        auto ts_xyz_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 64);

        // set up pointers to integrals

        auto tdx_xxyy_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 49);

        auto tdy_xxyy_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 49);

        auto tdz_xxyy_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 49);

        auto tdx_xxyy_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 50);

        auto tdy_xxyy_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 50);

        auto tdz_xxyy_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 50);

        auto tdx_xxyy_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 51);

        auto tdy_xxyy_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 51);

        auto tdz_xxyy_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 51);

        auto tdx_xxyy_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 52);

        auto tdy_xxyy_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 52);

        auto tdz_xxyy_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 52);

        auto tdx_xxyy_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 53);

        auto tdy_xxyy_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 53);

        auto tdz_xxyy_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 53);

        auto tdx_xxyy_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 54);

        auto tdy_xxyy_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 54);

        auto tdz_xxyy_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 54);

        auto tdx_xxyy_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 55);

        auto tdy_xxyy_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 55);

        auto tdz_xxyy_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 55);

        auto tdx_xxyy_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 56);

        auto tdy_xxyy_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 56);

        auto tdz_xxyy_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 56);

        auto tdx_xxyy_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 57);

        auto tdy_xxyy_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 57);

        auto tdz_xxyy_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 57);

        auto tdx_xxyy_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 58);

        auto tdy_xxyy_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 58);

        auto tdz_xxyy_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 58);

        auto tdx_xxyy_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 59);

        auto tdy_xxyy_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 59);

        auto tdz_xxyy_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 59);

        auto tdx_xxyz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 60);

        auto tdy_xxyz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 60);

        auto tdz_xxyz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 60);

        auto tdx_xxyz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 61);

        auto tdy_xxyz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 61);

        auto tdz_xxyz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 61);

        auto tdx_xxyz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 62);

        auto tdy_xxyz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 62);

        auto tdz_xxyz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 62);

        auto tdx_xxyz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 63);

        auto tdy_xxyz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 63);

        auto tdz_xxyz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 63);

        auto tdx_xxyz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 64);

        auto tdy_xxyz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 64);

        auto tdz_xxyz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 64);

        // Batch of Integrals (147,195)

        #pragma omp simd aligned(fx, pa_x, tdx_xxyy_xxyz_0, tdx_xxyy_xxzz_0, tdx_xxyy_xyyy_0, \
                                     tdx_xxyy_xyyz_0, tdx_xxyy_xyzz_0, tdx_xxyy_xzzz_0, tdx_xxyy_yyyy_0, tdx_xxyy_yyyz_0, \
                                     tdx_xxyy_yyzz_0, tdx_xxyy_yzzz_0, tdx_xxyy_zzzz_0, tdx_xxyz_xxxx_0, tdx_xxyz_xxxy_0, \
                                     tdx_xxyz_xxxz_0, tdx_xxyz_xxyy_0, tdx_xxyz_xxyz_0, tdx_xyy_xxyz_0, tdx_xyy_xxzz_0, \
                                     tdx_xyy_xyyy_0, tdx_xyy_xyyz_0, tdx_xyy_xyz_0, tdx_xyy_xyzz_0, tdx_xyy_xzz_0, \
                                     tdx_xyy_xzzz_0, tdx_xyy_yyy_0, tdx_xyy_yyyy_0, tdx_xyy_yyyz_0, tdx_xyy_yyz_0, \
                                     tdx_xyy_yyzz_0, tdx_xyy_yzz_0, tdx_xyy_yzzz_0, tdx_xyy_zzz_0, tdx_xyy_zzzz_0, \
                                     tdx_xyz_xxx_0, tdx_xyz_xxxx_0, tdx_xyz_xxxy_0, tdx_xyz_xxxz_0, tdx_xyz_xxy_0, \
                                     tdx_xyz_xxyy_0, tdx_xyz_xxyz_0, tdx_xyz_xxz_0, tdx_xyz_xyy_0, tdx_xyz_xyz_0, \
                                     tdx_yy_xxyz_0, tdx_yy_xxzz_0, tdx_yy_xyyy_0, tdx_yy_xyyz_0, tdx_yy_xyzz_0, \
                                     tdx_yy_xzzz_0, tdx_yy_yyyy_0, tdx_yy_yyyz_0, tdx_yy_yyzz_0, tdx_yy_yzzz_0, \
                                     tdx_yy_zzzz_0, tdx_yz_xxxx_0, tdx_yz_xxxy_0, tdx_yz_xxxz_0, tdx_yz_xxyy_0, \
                                     tdx_yz_xxyz_0, tdy_xxyy_xxyz_0, tdy_xxyy_xxzz_0, tdy_xxyy_xyyy_0, tdy_xxyy_xyyz_0, \
                                     tdy_xxyy_xyzz_0, tdy_xxyy_xzzz_0, tdy_xxyy_yyyy_0, tdy_xxyy_yyyz_0, tdy_xxyy_yyzz_0, \
                                     tdy_xxyy_yzzz_0, tdy_xxyy_zzzz_0, tdy_xxyz_xxxx_0, tdy_xxyz_xxxy_0, tdy_xxyz_xxxz_0, \
                                     tdy_xxyz_xxyy_0, tdy_xxyz_xxyz_0, tdy_xyy_xxyz_0, tdy_xyy_xxzz_0, tdy_xyy_xyyy_0, \
                                     tdy_xyy_xyyz_0, tdy_xyy_xyz_0, tdy_xyy_xyzz_0, tdy_xyy_xzz_0, tdy_xyy_xzzz_0, \
                                     tdy_xyy_yyy_0, tdy_xyy_yyyy_0, tdy_xyy_yyyz_0, tdy_xyy_yyz_0, tdy_xyy_yyzz_0, \
                                     tdy_xyy_yzz_0, tdy_xyy_yzzz_0, tdy_xyy_zzz_0, tdy_xyy_zzzz_0, tdy_xyz_xxx_0, \
                                     tdy_xyz_xxxx_0, tdy_xyz_xxxy_0, tdy_xyz_xxxz_0, tdy_xyz_xxy_0, tdy_xyz_xxyy_0, \
                                     tdy_xyz_xxyz_0, tdy_xyz_xxz_0, tdy_xyz_xyy_0, tdy_xyz_xyz_0, tdy_yy_xxyz_0, \
                                     tdy_yy_xxzz_0, tdy_yy_xyyy_0, tdy_yy_xyyz_0, tdy_yy_xyzz_0, tdy_yy_xzzz_0, \
                                     tdy_yy_yyyy_0, tdy_yy_yyyz_0, tdy_yy_yyzz_0, tdy_yy_yzzz_0, tdy_yy_zzzz_0, \
                                     tdy_yz_xxxx_0, tdy_yz_xxxy_0, tdy_yz_xxxz_0, tdy_yz_xxyy_0, tdy_yz_xxyz_0, \
                                     tdz_xxyy_xxyz_0, tdz_xxyy_xxzz_0, tdz_xxyy_xyyy_0, tdz_xxyy_xyyz_0, tdz_xxyy_xyzz_0, \
                                     tdz_xxyy_xzzz_0, tdz_xxyy_yyyy_0, tdz_xxyy_yyyz_0, tdz_xxyy_yyzz_0, tdz_xxyy_yzzz_0, \
                                     tdz_xxyy_zzzz_0, tdz_xxyz_xxxx_0, tdz_xxyz_xxxy_0, tdz_xxyz_xxxz_0, tdz_xxyz_xxyy_0, \
                                     tdz_xxyz_xxyz_0, tdz_xyy_xxyz_0, tdz_xyy_xxzz_0, tdz_xyy_xyyy_0, tdz_xyy_xyyz_0, \
                                     tdz_xyy_xyz_0, tdz_xyy_xyzz_0, tdz_xyy_xzz_0, tdz_xyy_xzzz_0, tdz_xyy_yyy_0, \
                                     tdz_xyy_yyyy_0, tdz_xyy_yyyz_0, tdz_xyy_yyz_0, tdz_xyy_yyzz_0, tdz_xyy_yzz_0, \
                                     tdz_xyy_yzzz_0, tdz_xyy_zzz_0, tdz_xyy_zzzz_0, tdz_xyz_xxx_0, tdz_xyz_xxxx_0, \
                                     tdz_xyz_xxxy_0, tdz_xyz_xxxz_0, tdz_xyz_xxy_0, tdz_xyz_xxyy_0, tdz_xyz_xxyz_0, \
                                     tdz_xyz_xxz_0, tdz_xyz_xyy_0, tdz_xyz_xyz_0, tdz_yy_xxyz_0, tdz_yy_xxzz_0, \
                                     tdz_yy_xyyy_0, tdz_yy_xyyz_0, tdz_yy_xyzz_0, tdz_yy_xzzz_0, tdz_yy_yyyy_0, \
                                     tdz_yy_yyyz_0, tdz_yy_yyzz_0, tdz_yy_yzzz_0, tdz_yy_zzzz_0, tdz_yz_xxxx_0, \
                                     tdz_yz_xxxy_0, tdz_yz_xxxz_0, tdz_yz_xxyy_0, tdz_yz_xxyz_0, ts_xyy_xxyz_0, \
                                     ts_xyy_xxzz_0, ts_xyy_xyyy_0, ts_xyy_xyyz_0, ts_xyy_xyzz_0, ts_xyy_xzzz_0, \
                                     ts_xyy_yyyy_0, ts_xyy_yyyz_0, ts_xyy_yyzz_0, ts_xyy_yzzz_0, ts_xyy_zzzz_0, \
                                     ts_xyz_xxxx_0, ts_xyz_xxxy_0, ts_xyz_xxxz_0, ts_xyz_xxyy_0, ts_xyz_xxyz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            tdx_xxyy_xxyz_0[j] =
                pa_x[j] * tdx_xyy_xxyz_0[j] + 0.5 * fl1_fx * tdx_yy_xxyz_0[j] + fl1_fx * tdx_xyy_xyz_0[j] + 0.5 * fl1_fx * ts_xyy_xxyz_0[j];

            tdy_xxyy_xxyz_0[j] = pa_x[j] * tdy_xyy_xxyz_0[j] + 0.5 * fl1_fx * tdy_yy_xxyz_0[j] + fl1_fx * tdy_xyy_xyz_0[j];

            tdz_xxyy_xxyz_0[j] = pa_x[j] * tdz_xyy_xxyz_0[j] + 0.5 * fl1_fx * tdz_yy_xxyz_0[j] + fl1_fx * tdz_xyy_xyz_0[j];

            tdx_xxyy_xxzz_0[j] =
                pa_x[j] * tdx_xyy_xxzz_0[j] + 0.5 * fl1_fx * tdx_yy_xxzz_0[j] + fl1_fx * tdx_xyy_xzz_0[j] + 0.5 * fl1_fx * ts_xyy_xxzz_0[j];

            tdy_xxyy_xxzz_0[j] = pa_x[j] * tdy_xyy_xxzz_0[j] + 0.5 * fl1_fx * tdy_yy_xxzz_0[j] + fl1_fx * tdy_xyy_xzz_0[j];

            tdz_xxyy_xxzz_0[j] = pa_x[j] * tdz_xyy_xxzz_0[j] + 0.5 * fl1_fx * tdz_yy_xxzz_0[j] + fl1_fx * tdz_xyy_xzz_0[j];

            tdx_xxyy_xyyy_0[j] =
                pa_x[j] * tdx_xyy_xyyy_0[j] + 0.5 * fl1_fx * tdx_yy_xyyy_0[j] + 0.5 * fl1_fx * tdx_xyy_yyy_0[j] + 0.5 * fl1_fx * ts_xyy_xyyy_0[j];

            tdy_xxyy_xyyy_0[j] = pa_x[j] * tdy_xyy_xyyy_0[j] + 0.5 * fl1_fx * tdy_yy_xyyy_0[j] + 0.5 * fl1_fx * tdy_xyy_yyy_0[j];

            tdz_xxyy_xyyy_0[j] = pa_x[j] * tdz_xyy_xyyy_0[j] + 0.5 * fl1_fx * tdz_yy_xyyy_0[j] + 0.5 * fl1_fx * tdz_xyy_yyy_0[j];

            tdx_xxyy_xyyz_0[j] =
                pa_x[j] * tdx_xyy_xyyz_0[j] + 0.5 * fl1_fx * tdx_yy_xyyz_0[j] + 0.5 * fl1_fx * tdx_xyy_yyz_0[j] + 0.5 * fl1_fx * ts_xyy_xyyz_0[j];

            tdy_xxyy_xyyz_0[j] = pa_x[j] * tdy_xyy_xyyz_0[j] + 0.5 * fl1_fx * tdy_yy_xyyz_0[j] + 0.5 * fl1_fx * tdy_xyy_yyz_0[j];

            tdz_xxyy_xyyz_0[j] = pa_x[j] * tdz_xyy_xyyz_0[j] + 0.5 * fl1_fx * tdz_yy_xyyz_0[j] + 0.5 * fl1_fx * tdz_xyy_yyz_0[j];

            tdx_xxyy_xyzz_0[j] =
                pa_x[j] * tdx_xyy_xyzz_0[j] + 0.5 * fl1_fx * tdx_yy_xyzz_0[j] + 0.5 * fl1_fx * tdx_xyy_yzz_0[j] + 0.5 * fl1_fx * ts_xyy_xyzz_0[j];

            tdy_xxyy_xyzz_0[j] = pa_x[j] * tdy_xyy_xyzz_0[j] + 0.5 * fl1_fx * tdy_yy_xyzz_0[j] + 0.5 * fl1_fx * tdy_xyy_yzz_0[j];

            tdz_xxyy_xyzz_0[j] = pa_x[j] * tdz_xyy_xyzz_0[j] + 0.5 * fl1_fx * tdz_yy_xyzz_0[j] + 0.5 * fl1_fx * tdz_xyy_yzz_0[j];

            tdx_xxyy_xzzz_0[j] =
                pa_x[j] * tdx_xyy_xzzz_0[j] + 0.5 * fl1_fx * tdx_yy_xzzz_0[j] + 0.5 * fl1_fx * tdx_xyy_zzz_0[j] + 0.5 * fl1_fx * ts_xyy_xzzz_0[j];

            tdy_xxyy_xzzz_0[j] = pa_x[j] * tdy_xyy_xzzz_0[j] + 0.5 * fl1_fx * tdy_yy_xzzz_0[j] + 0.5 * fl1_fx * tdy_xyy_zzz_0[j];

            tdz_xxyy_xzzz_0[j] = pa_x[j] * tdz_xyy_xzzz_0[j] + 0.5 * fl1_fx * tdz_yy_xzzz_0[j] + 0.5 * fl1_fx * tdz_xyy_zzz_0[j];

            tdx_xxyy_yyyy_0[j] = pa_x[j] * tdx_xyy_yyyy_0[j] + 0.5 * fl1_fx * tdx_yy_yyyy_0[j] + 0.5 * fl1_fx * ts_xyy_yyyy_0[j];

            tdy_xxyy_yyyy_0[j] = pa_x[j] * tdy_xyy_yyyy_0[j] + 0.5 * fl1_fx * tdy_yy_yyyy_0[j];

            tdz_xxyy_yyyy_0[j] = pa_x[j] * tdz_xyy_yyyy_0[j] + 0.5 * fl1_fx * tdz_yy_yyyy_0[j];

            tdx_xxyy_yyyz_0[j] = pa_x[j] * tdx_xyy_yyyz_0[j] + 0.5 * fl1_fx * tdx_yy_yyyz_0[j] + 0.5 * fl1_fx * ts_xyy_yyyz_0[j];

            tdy_xxyy_yyyz_0[j] = pa_x[j] * tdy_xyy_yyyz_0[j] + 0.5 * fl1_fx * tdy_yy_yyyz_0[j];

            tdz_xxyy_yyyz_0[j] = pa_x[j] * tdz_xyy_yyyz_0[j] + 0.5 * fl1_fx * tdz_yy_yyyz_0[j];

            tdx_xxyy_yyzz_0[j] = pa_x[j] * tdx_xyy_yyzz_0[j] + 0.5 * fl1_fx * tdx_yy_yyzz_0[j] + 0.5 * fl1_fx * ts_xyy_yyzz_0[j];

            tdy_xxyy_yyzz_0[j] = pa_x[j] * tdy_xyy_yyzz_0[j] + 0.5 * fl1_fx * tdy_yy_yyzz_0[j];

            tdz_xxyy_yyzz_0[j] = pa_x[j] * tdz_xyy_yyzz_0[j] + 0.5 * fl1_fx * tdz_yy_yyzz_0[j];

            tdx_xxyy_yzzz_0[j] = pa_x[j] * tdx_xyy_yzzz_0[j] + 0.5 * fl1_fx * tdx_yy_yzzz_0[j] + 0.5 * fl1_fx * ts_xyy_yzzz_0[j];

            tdy_xxyy_yzzz_0[j] = pa_x[j] * tdy_xyy_yzzz_0[j] + 0.5 * fl1_fx * tdy_yy_yzzz_0[j];

            tdz_xxyy_yzzz_0[j] = pa_x[j] * tdz_xyy_yzzz_0[j] + 0.5 * fl1_fx * tdz_yy_yzzz_0[j];

            tdx_xxyy_zzzz_0[j] = pa_x[j] * tdx_xyy_zzzz_0[j] + 0.5 * fl1_fx * tdx_yy_zzzz_0[j] + 0.5 * fl1_fx * ts_xyy_zzzz_0[j];

            tdy_xxyy_zzzz_0[j] = pa_x[j] * tdy_xyy_zzzz_0[j] + 0.5 * fl1_fx * tdy_yy_zzzz_0[j];

            tdz_xxyy_zzzz_0[j] = pa_x[j] * tdz_xyy_zzzz_0[j] + 0.5 * fl1_fx * tdz_yy_zzzz_0[j];

            tdx_xxyz_xxxx_0[j] =
                pa_x[j] * tdx_xyz_xxxx_0[j] + 0.5 * fl1_fx * tdx_yz_xxxx_0[j] + 2.0 * fl1_fx * tdx_xyz_xxx_0[j] + 0.5 * fl1_fx * ts_xyz_xxxx_0[j];

            tdy_xxyz_xxxx_0[j] = pa_x[j] * tdy_xyz_xxxx_0[j] + 0.5 * fl1_fx * tdy_yz_xxxx_0[j] + 2.0 * fl1_fx * tdy_xyz_xxx_0[j];

            tdz_xxyz_xxxx_0[j] = pa_x[j] * tdz_xyz_xxxx_0[j] + 0.5 * fl1_fx * tdz_yz_xxxx_0[j] + 2.0 * fl1_fx * tdz_xyz_xxx_0[j];

            tdx_xxyz_xxxy_0[j] =
                pa_x[j] * tdx_xyz_xxxy_0[j] + 0.5 * fl1_fx * tdx_yz_xxxy_0[j] + 1.5 * fl1_fx * tdx_xyz_xxy_0[j] + 0.5 * fl1_fx * ts_xyz_xxxy_0[j];

            tdy_xxyz_xxxy_0[j] = pa_x[j] * tdy_xyz_xxxy_0[j] + 0.5 * fl1_fx * tdy_yz_xxxy_0[j] + 1.5 * fl1_fx * tdy_xyz_xxy_0[j];

            tdz_xxyz_xxxy_0[j] = pa_x[j] * tdz_xyz_xxxy_0[j] + 0.5 * fl1_fx * tdz_yz_xxxy_0[j] + 1.5 * fl1_fx * tdz_xyz_xxy_0[j];

            tdx_xxyz_xxxz_0[j] =
                pa_x[j] * tdx_xyz_xxxz_0[j] + 0.5 * fl1_fx * tdx_yz_xxxz_0[j] + 1.5 * fl1_fx * tdx_xyz_xxz_0[j] + 0.5 * fl1_fx * ts_xyz_xxxz_0[j];

            tdy_xxyz_xxxz_0[j] = pa_x[j] * tdy_xyz_xxxz_0[j] + 0.5 * fl1_fx * tdy_yz_xxxz_0[j] + 1.5 * fl1_fx * tdy_xyz_xxz_0[j];

            tdz_xxyz_xxxz_0[j] = pa_x[j] * tdz_xyz_xxxz_0[j] + 0.5 * fl1_fx * tdz_yz_xxxz_0[j] + 1.5 * fl1_fx * tdz_xyz_xxz_0[j];

            tdx_xxyz_xxyy_0[j] =
                pa_x[j] * tdx_xyz_xxyy_0[j] + 0.5 * fl1_fx * tdx_yz_xxyy_0[j] + fl1_fx * tdx_xyz_xyy_0[j] + 0.5 * fl1_fx * ts_xyz_xxyy_0[j];

            tdy_xxyz_xxyy_0[j] = pa_x[j] * tdy_xyz_xxyy_0[j] + 0.5 * fl1_fx * tdy_yz_xxyy_0[j] + fl1_fx * tdy_xyz_xyy_0[j];

            tdz_xxyz_xxyy_0[j] = pa_x[j] * tdz_xyz_xxyy_0[j] + 0.5 * fl1_fx * tdz_yz_xxyy_0[j] + fl1_fx * tdz_xyz_xyy_0[j];

            tdx_xxyz_xxyz_0[j] =
                pa_x[j] * tdx_xyz_xxyz_0[j] + 0.5 * fl1_fx * tdx_yz_xxyz_0[j] + fl1_fx * tdx_xyz_xyz_0[j] + 0.5 * fl1_fx * ts_xyz_xxyz_0[j];

            tdy_xxyz_xxyz_0[j] = pa_x[j] * tdy_xyz_xxyz_0[j] + 0.5 * fl1_fx * tdy_yz_xxyz_0[j] + fl1_fx * tdy_xyz_xyz_0[j];

            tdz_xxyz_xxyz_0[j] = pa_x[j] * tdz_xyz_xxyz_0[j] + 0.5 * fl1_fx * tdz_yz_xxyz_0[j] + fl1_fx * tdz_xyz_xyz_0[j];
        }

        idx++;
    }
}

void
compElectricDipoleForGG_195_243(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const int32_t              nOSFactors,
                                const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_d_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_d_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

        auto tdx_xyz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 65);

        auto tdy_xyz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 65);

        auto tdz_xyz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 65);

        auto tdx_xyz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 66);

        auto tdy_xyz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 66);

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

        auto tdx_xyz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 45);

        auto tdy_xyz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 45);

        auto tdz_xyz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 45);

        auto tdx_xyz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 46);

        auto tdy_xyz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 46);

        auto tdz_xyz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 46);

        auto tdx_xyz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 47);

        auto tdy_xyz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 47);

        auto tdz_xyz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 47);

        auto tdx_xyz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 48);

        auto tdy_xyz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 48);

        auto tdz_xyz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 48);

        auto tdx_xyz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 49);

        auto tdy_xyz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 49);

        auto tdz_xyz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 49);

        auto tdx_xzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 50);

        auto tdy_xzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 50);

        auto tdz_xzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 50);

        auto tdx_xzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 51);

        auto tdy_xzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 51);

        auto tdz_xzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 51);

        auto tdx_xzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 52);

        auto tdy_xzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 52);

        auto tdz_xzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 52);

        auto tdx_xzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 53);

        auto tdy_xzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 53);

        auto tdz_xzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 53);

        auto tdx_xzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 54);

        auto tdy_xzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 54);

        auto tdz_xzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 54);

        auto tdx_xzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 55);

        auto tdy_xzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 55);

        auto tdz_xzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 55);

        auto ts_xyz_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 65);

        auto ts_xyz_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 66);

        auto ts_xyz_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 67);

        auto ts_xyz_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 68);

        auto ts_xyz_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 69);

        auto ts_xyz_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 70);

        auto ts_xyz_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 71);

        auto ts_xyz_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 72);

        auto ts_xyz_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 73);

        auto ts_xyz_zzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 74);

        auto ts_xzz_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 75);

        auto ts_xzz_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 76);

        auto ts_xzz_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 77);

        auto ts_xzz_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 78);

        auto ts_xzz_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 79);

        auto ts_xzz_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 80);

        // set up pointers to integrals

        auto tdx_xxyz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 65);

        auto tdy_xxyz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 65);

        auto tdz_xxyz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 65);

        auto tdx_xxyz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 66);

        auto tdy_xxyz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 66);

        auto tdz_xxyz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 66);

        auto tdx_xxyz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 67);

        auto tdy_xxyz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 67);

        auto tdz_xxyz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 67);

        auto tdx_xxyz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 68);

        auto tdy_xxyz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 68);

        auto tdz_xxyz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 68);

        auto tdx_xxyz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 69);

        auto tdy_xxyz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 69);

        auto tdz_xxyz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 69);

        auto tdx_xxyz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 70);

        auto tdy_xxyz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 70);

        auto tdz_xxyz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 70);

        auto tdx_xxyz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 71);

        auto tdy_xxyz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 71);

        auto tdz_xxyz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 71);

        auto tdx_xxyz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 72);

        auto tdy_xxyz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 72);

        auto tdz_xxyz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 72);

        auto tdx_xxyz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 73);

        auto tdy_xxyz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 73);

        auto tdz_xxyz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 73);

        auto tdx_xxyz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 74);

        auto tdy_xxyz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 74);

        auto tdz_xxyz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 74);

        auto tdx_xxzz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 75);

        auto tdy_xxzz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 75);

        auto tdz_xxzz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 75);

        auto tdx_xxzz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 76);

        auto tdy_xxzz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 76);

        auto tdz_xxzz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 76);

        auto tdx_xxzz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 77);

        auto tdy_xxzz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 77);

        auto tdz_xxzz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 77);

        auto tdx_xxzz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 78);

        auto tdy_xxzz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 78);

        auto tdz_xxzz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 78);

        auto tdx_xxzz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 79);

        auto tdy_xxzz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 79);

        auto tdz_xxzz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 79);

        auto tdx_xxzz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 80);

        auto tdy_xxzz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 80);

        auto tdz_xxzz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 80);

        // Batch of Integrals (195,243)

        #pragma omp simd aligned(fx, pa_x, tdx_xxyz_xxzz_0, tdx_xxyz_xyyy_0, tdx_xxyz_xyyz_0, \
                                     tdx_xxyz_xyzz_0, tdx_xxyz_xzzz_0, tdx_xxyz_yyyy_0, tdx_xxyz_yyyz_0, tdx_xxyz_yyzz_0, \
                                     tdx_xxyz_yzzz_0, tdx_xxyz_zzzz_0, tdx_xxzz_xxxx_0, tdx_xxzz_xxxy_0, tdx_xxzz_xxxz_0, \
                                     tdx_xxzz_xxyy_0, tdx_xxzz_xxyz_0, tdx_xxzz_xxzz_0, tdx_xyz_xxzz_0, tdx_xyz_xyyy_0, \
                                     tdx_xyz_xyyz_0, tdx_xyz_xyzz_0, tdx_xyz_xzz_0, tdx_xyz_xzzz_0, tdx_xyz_yyy_0, \
                                     tdx_xyz_yyyy_0, tdx_xyz_yyyz_0, tdx_xyz_yyz_0, tdx_xyz_yyzz_0, tdx_xyz_yzz_0, \
                                     tdx_xyz_yzzz_0, tdx_xyz_zzz_0, tdx_xyz_zzzz_0, tdx_xzz_xxx_0, tdx_xzz_xxxx_0, \
                                     tdx_xzz_xxxy_0, tdx_xzz_xxxz_0, tdx_xzz_xxy_0, tdx_xzz_xxyy_0, tdx_xzz_xxyz_0, \
                                     tdx_xzz_xxz_0, tdx_xzz_xxzz_0, tdx_xzz_xyy_0, tdx_xzz_xyz_0, tdx_xzz_xzz_0, \
                                     tdx_yz_xxzz_0, tdx_yz_xyyy_0, tdx_yz_xyyz_0, tdx_yz_xyzz_0, tdx_yz_xzzz_0, \
                                     tdx_yz_yyyy_0, tdx_yz_yyyz_0, tdx_yz_yyzz_0, tdx_yz_yzzz_0, tdx_yz_zzzz_0, \
                                     tdx_zz_xxxx_0, tdx_zz_xxxy_0, tdx_zz_xxxz_0, tdx_zz_xxyy_0, tdx_zz_xxyz_0, \
                                     tdx_zz_xxzz_0, tdy_xxyz_xxzz_0, tdy_xxyz_xyyy_0, tdy_xxyz_xyyz_0, tdy_xxyz_xyzz_0, \
                                     tdy_xxyz_xzzz_0, tdy_xxyz_yyyy_0, tdy_xxyz_yyyz_0, tdy_xxyz_yyzz_0, tdy_xxyz_yzzz_0, \
                                     tdy_xxyz_zzzz_0, tdy_xxzz_xxxx_0, tdy_xxzz_xxxy_0, tdy_xxzz_xxxz_0, tdy_xxzz_xxyy_0, \
                                     tdy_xxzz_xxyz_0, tdy_xxzz_xxzz_0, tdy_xyz_xxzz_0, tdy_xyz_xyyy_0, tdy_xyz_xyyz_0, \
                                     tdy_xyz_xyzz_0, tdy_xyz_xzz_0, tdy_xyz_xzzz_0, tdy_xyz_yyy_0, tdy_xyz_yyyy_0, \
                                     tdy_xyz_yyyz_0, tdy_xyz_yyz_0, tdy_xyz_yyzz_0, tdy_xyz_yzz_0, tdy_xyz_yzzz_0, \
                                     tdy_xyz_zzz_0, tdy_xyz_zzzz_0, tdy_xzz_xxx_0, tdy_xzz_xxxx_0, tdy_xzz_xxxy_0, \
                                     tdy_xzz_xxxz_0, tdy_xzz_xxy_0, tdy_xzz_xxyy_0, tdy_xzz_xxyz_0, tdy_xzz_xxz_0, \
                                     tdy_xzz_xxzz_0, tdy_xzz_xyy_0, tdy_xzz_xyz_0, tdy_xzz_xzz_0, tdy_yz_xxzz_0, \
                                     tdy_yz_xyyy_0, tdy_yz_xyyz_0, tdy_yz_xyzz_0, tdy_yz_xzzz_0, tdy_yz_yyyy_0, \
                                     tdy_yz_yyyz_0, tdy_yz_yyzz_0, tdy_yz_yzzz_0, tdy_yz_zzzz_0, tdy_zz_xxxx_0, \
                                     tdy_zz_xxxy_0, tdy_zz_xxxz_0, tdy_zz_xxyy_0, tdy_zz_xxyz_0, tdy_zz_xxzz_0, \
                                     tdz_xxyz_xxzz_0, tdz_xxyz_xyyy_0, tdz_xxyz_xyyz_0, tdz_xxyz_xyzz_0, tdz_xxyz_xzzz_0, \
                                     tdz_xxyz_yyyy_0, tdz_xxyz_yyyz_0, tdz_xxyz_yyzz_0, tdz_xxyz_yzzz_0, tdz_xxyz_zzzz_0, \
                                     tdz_xxzz_xxxx_0, tdz_xxzz_xxxy_0, tdz_xxzz_xxxz_0, tdz_xxzz_xxyy_0, tdz_xxzz_xxyz_0, \
                                     tdz_xxzz_xxzz_0, tdz_xyz_xxzz_0, tdz_xyz_xyyy_0, tdz_xyz_xyyz_0, tdz_xyz_xyzz_0, \
                                     tdz_xyz_xzz_0, tdz_xyz_xzzz_0, tdz_xyz_yyy_0, tdz_xyz_yyyy_0, tdz_xyz_yyyz_0, \
                                     tdz_xyz_yyz_0, tdz_xyz_yyzz_0, tdz_xyz_yzz_0, tdz_xyz_yzzz_0, tdz_xyz_zzz_0, \
                                     tdz_xyz_zzzz_0, tdz_xzz_xxx_0, tdz_xzz_xxxx_0, tdz_xzz_xxxy_0, tdz_xzz_xxxz_0, \
                                     tdz_xzz_xxy_0, tdz_xzz_xxyy_0, tdz_xzz_xxyz_0, tdz_xzz_xxz_0, tdz_xzz_xxzz_0, \
                                     tdz_xzz_xyy_0, tdz_xzz_xyz_0, tdz_xzz_xzz_0, tdz_yz_xxzz_0, tdz_yz_xyyy_0, \
                                     tdz_yz_xyyz_0, tdz_yz_xyzz_0, tdz_yz_xzzz_0, tdz_yz_yyyy_0, tdz_yz_yyyz_0, \
                                     tdz_yz_yyzz_0, tdz_yz_yzzz_0, tdz_yz_zzzz_0, tdz_zz_xxxx_0, tdz_zz_xxxy_0, \
                                     tdz_zz_xxxz_0, tdz_zz_xxyy_0, tdz_zz_xxyz_0, tdz_zz_xxzz_0, ts_xyz_xxzz_0, \
                                     ts_xyz_xyyy_0, ts_xyz_xyyz_0, ts_xyz_xyzz_0, ts_xyz_xzzz_0, ts_xyz_yyyy_0, \
                                     ts_xyz_yyyz_0, ts_xyz_yyzz_0, ts_xyz_yzzz_0, ts_xyz_zzzz_0, ts_xzz_xxxx_0, \
                                     ts_xzz_xxxy_0, ts_xzz_xxxz_0, ts_xzz_xxyy_0, ts_xzz_xxyz_0, ts_xzz_xxzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            tdx_xxyz_xxzz_0[j] =
                pa_x[j] * tdx_xyz_xxzz_0[j] + 0.5 * fl1_fx * tdx_yz_xxzz_0[j] + fl1_fx * tdx_xyz_xzz_0[j] + 0.5 * fl1_fx * ts_xyz_xxzz_0[j];

            tdy_xxyz_xxzz_0[j] = pa_x[j] * tdy_xyz_xxzz_0[j] + 0.5 * fl1_fx * tdy_yz_xxzz_0[j] + fl1_fx * tdy_xyz_xzz_0[j];

            tdz_xxyz_xxzz_0[j] = pa_x[j] * tdz_xyz_xxzz_0[j] + 0.5 * fl1_fx * tdz_yz_xxzz_0[j] + fl1_fx * tdz_xyz_xzz_0[j];

            tdx_xxyz_xyyy_0[j] =
                pa_x[j] * tdx_xyz_xyyy_0[j] + 0.5 * fl1_fx * tdx_yz_xyyy_0[j] + 0.5 * fl1_fx * tdx_xyz_yyy_0[j] + 0.5 * fl1_fx * ts_xyz_xyyy_0[j];

            tdy_xxyz_xyyy_0[j] = pa_x[j] * tdy_xyz_xyyy_0[j] + 0.5 * fl1_fx * tdy_yz_xyyy_0[j] + 0.5 * fl1_fx * tdy_xyz_yyy_0[j];

            tdz_xxyz_xyyy_0[j] = pa_x[j] * tdz_xyz_xyyy_0[j] + 0.5 * fl1_fx * tdz_yz_xyyy_0[j] + 0.5 * fl1_fx * tdz_xyz_yyy_0[j];

            tdx_xxyz_xyyz_0[j] =
                pa_x[j] * tdx_xyz_xyyz_0[j] + 0.5 * fl1_fx * tdx_yz_xyyz_0[j] + 0.5 * fl1_fx * tdx_xyz_yyz_0[j] + 0.5 * fl1_fx * ts_xyz_xyyz_0[j];

            tdy_xxyz_xyyz_0[j] = pa_x[j] * tdy_xyz_xyyz_0[j] + 0.5 * fl1_fx * tdy_yz_xyyz_0[j] + 0.5 * fl1_fx * tdy_xyz_yyz_0[j];

            tdz_xxyz_xyyz_0[j] = pa_x[j] * tdz_xyz_xyyz_0[j] + 0.5 * fl1_fx * tdz_yz_xyyz_0[j] + 0.5 * fl1_fx * tdz_xyz_yyz_0[j];

            tdx_xxyz_xyzz_0[j] =
                pa_x[j] * tdx_xyz_xyzz_0[j] + 0.5 * fl1_fx * tdx_yz_xyzz_0[j] + 0.5 * fl1_fx * tdx_xyz_yzz_0[j] + 0.5 * fl1_fx * ts_xyz_xyzz_0[j];

            tdy_xxyz_xyzz_0[j] = pa_x[j] * tdy_xyz_xyzz_0[j] + 0.5 * fl1_fx * tdy_yz_xyzz_0[j] + 0.5 * fl1_fx * tdy_xyz_yzz_0[j];

            tdz_xxyz_xyzz_0[j] = pa_x[j] * tdz_xyz_xyzz_0[j] + 0.5 * fl1_fx * tdz_yz_xyzz_0[j] + 0.5 * fl1_fx * tdz_xyz_yzz_0[j];

            tdx_xxyz_xzzz_0[j] =
                pa_x[j] * tdx_xyz_xzzz_0[j] + 0.5 * fl1_fx * tdx_yz_xzzz_0[j] + 0.5 * fl1_fx * tdx_xyz_zzz_0[j] + 0.5 * fl1_fx * ts_xyz_xzzz_0[j];

            tdy_xxyz_xzzz_0[j] = pa_x[j] * tdy_xyz_xzzz_0[j] + 0.5 * fl1_fx * tdy_yz_xzzz_0[j] + 0.5 * fl1_fx * tdy_xyz_zzz_0[j];

            tdz_xxyz_xzzz_0[j] = pa_x[j] * tdz_xyz_xzzz_0[j] + 0.5 * fl1_fx * tdz_yz_xzzz_0[j] + 0.5 * fl1_fx * tdz_xyz_zzz_0[j];

            tdx_xxyz_yyyy_0[j] = pa_x[j] * tdx_xyz_yyyy_0[j] + 0.5 * fl1_fx * tdx_yz_yyyy_0[j] + 0.5 * fl1_fx * ts_xyz_yyyy_0[j];

            tdy_xxyz_yyyy_0[j] = pa_x[j] * tdy_xyz_yyyy_0[j] + 0.5 * fl1_fx * tdy_yz_yyyy_0[j];

            tdz_xxyz_yyyy_0[j] = pa_x[j] * tdz_xyz_yyyy_0[j] + 0.5 * fl1_fx * tdz_yz_yyyy_0[j];

            tdx_xxyz_yyyz_0[j] = pa_x[j] * tdx_xyz_yyyz_0[j] + 0.5 * fl1_fx * tdx_yz_yyyz_0[j] + 0.5 * fl1_fx * ts_xyz_yyyz_0[j];

            tdy_xxyz_yyyz_0[j] = pa_x[j] * tdy_xyz_yyyz_0[j] + 0.5 * fl1_fx * tdy_yz_yyyz_0[j];

            tdz_xxyz_yyyz_0[j] = pa_x[j] * tdz_xyz_yyyz_0[j] + 0.5 * fl1_fx * tdz_yz_yyyz_0[j];

            tdx_xxyz_yyzz_0[j] = pa_x[j] * tdx_xyz_yyzz_0[j] + 0.5 * fl1_fx * tdx_yz_yyzz_0[j] + 0.5 * fl1_fx * ts_xyz_yyzz_0[j];

            tdy_xxyz_yyzz_0[j] = pa_x[j] * tdy_xyz_yyzz_0[j] + 0.5 * fl1_fx * tdy_yz_yyzz_0[j];

            tdz_xxyz_yyzz_0[j] = pa_x[j] * tdz_xyz_yyzz_0[j] + 0.5 * fl1_fx * tdz_yz_yyzz_0[j];

            tdx_xxyz_yzzz_0[j] = pa_x[j] * tdx_xyz_yzzz_0[j] + 0.5 * fl1_fx * tdx_yz_yzzz_0[j] + 0.5 * fl1_fx * ts_xyz_yzzz_0[j];

            tdy_xxyz_yzzz_0[j] = pa_x[j] * tdy_xyz_yzzz_0[j] + 0.5 * fl1_fx * tdy_yz_yzzz_0[j];

            tdz_xxyz_yzzz_0[j] = pa_x[j] * tdz_xyz_yzzz_0[j] + 0.5 * fl1_fx * tdz_yz_yzzz_0[j];

            tdx_xxyz_zzzz_0[j] = pa_x[j] * tdx_xyz_zzzz_0[j] + 0.5 * fl1_fx * tdx_yz_zzzz_0[j] + 0.5 * fl1_fx * ts_xyz_zzzz_0[j];

            tdy_xxyz_zzzz_0[j] = pa_x[j] * tdy_xyz_zzzz_0[j] + 0.5 * fl1_fx * tdy_yz_zzzz_0[j];

            tdz_xxyz_zzzz_0[j] = pa_x[j] * tdz_xyz_zzzz_0[j] + 0.5 * fl1_fx * tdz_yz_zzzz_0[j];

            tdx_xxzz_xxxx_0[j] =
                pa_x[j] * tdx_xzz_xxxx_0[j] + 0.5 * fl1_fx * tdx_zz_xxxx_0[j] + 2.0 * fl1_fx * tdx_xzz_xxx_0[j] + 0.5 * fl1_fx * ts_xzz_xxxx_0[j];

            tdy_xxzz_xxxx_0[j] = pa_x[j] * tdy_xzz_xxxx_0[j] + 0.5 * fl1_fx * tdy_zz_xxxx_0[j] + 2.0 * fl1_fx * tdy_xzz_xxx_0[j];

            tdz_xxzz_xxxx_0[j] = pa_x[j] * tdz_xzz_xxxx_0[j] + 0.5 * fl1_fx * tdz_zz_xxxx_0[j] + 2.0 * fl1_fx * tdz_xzz_xxx_0[j];

            tdx_xxzz_xxxy_0[j] =
                pa_x[j] * tdx_xzz_xxxy_0[j] + 0.5 * fl1_fx * tdx_zz_xxxy_0[j] + 1.5 * fl1_fx * tdx_xzz_xxy_0[j] + 0.5 * fl1_fx * ts_xzz_xxxy_0[j];

            tdy_xxzz_xxxy_0[j] = pa_x[j] * tdy_xzz_xxxy_0[j] + 0.5 * fl1_fx * tdy_zz_xxxy_0[j] + 1.5 * fl1_fx * tdy_xzz_xxy_0[j];

            tdz_xxzz_xxxy_0[j] = pa_x[j] * tdz_xzz_xxxy_0[j] + 0.5 * fl1_fx * tdz_zz_xxxy_0[j] + 1.5 * fl1_fx * tdz_xzz_xxy_0[j];

            tdx_xxzz_xxxz_0[j] =
                pa_x[j] * tdx_xzz_xxxz_0[j] + 0.5 * fl1_fx * tdx_zz_xxxz_0[j] + 1.5 * fl1_fx * tdx_xzz_xxz_0[j] + 0.5 * fl1_fx * ts_xzz_xxxz_0[j];

            tdy_xxzz_xxxz_0[j] = pa_x[j] * tdy_xzz_xxxz_0[j] + 0.5 * fl1_fx * tdy_zz_xxxz_0[j] + 1.5 * fl1_fx * tdy_xzz_xxz_0[j];

            tdz_xxzz_xxxz_0[j] = pa_x[j] * tdz_xzz_xxxz_0[j] + 0.5 * fl1_fx * tdz_zz_xxxz_0[j] + 1.5 * fl1_fx * tdz_xzz_xxz_0[j];

            tdx_xxzz_xxyy_0[j] =
                pa_x[j] * tdx_xzz_xxyy_0[j] + 0.5 * fl1_fx * tdx_zz_xxyy_0[j] + fl1_fx * tdx_xzz_xyy_0[j] + 0.5 * fl1_fx * ts_xzz_xxyy_0[j];

            tdy_xxzz_xxyy_0[j] = pa_x[j] * tdy_xzz_xxyy_0[j] + 0.5 * fl1_fx * tdy_zz_xxyy_0[j] + fl1_fx * tdy_xzz_xyy_0[j];

            tdz_xxzz_xxyy_0[j] = pa_x[j] * tdz_xzz_xxyy_0[j] + 0.5 * fl1_fx * tdz_zz_xxyy_0[j] + fl1_fx * tdz_xzz_xyy_0[j];

            tdx_xxzz_xxyz_0[j] =
                pa_x[j] * tdx_xzz_xxyz_0[j] + 0.5 * fl1_fx * tdx_zz_xxyz_0[j] + fl1_fx * tdx_xzz_xyz_0[j] + 0.5 * fl1_fx * ts_xzz_xxyz_0[j];

            tdy_xxzz_xxyz_0[j] = pa_x[j] * tdy_xzz_xxyz_0[j] + 0.5 * fl1_fx * tdy_zz_xxyz_0[j] + fl1_fx * tdy_xzz_xyz_0[j];

            tdz_xxzz_xxyz_0[j] = pa_x[j] * tdz_xzz_xxyz_0[j] + 0.5 * fl1_fx * tdz_zz_xxyz_0[j] + fl1_fx * tdz_xzz_xyz_0[j];

            tdx_xxzz_xxzz_0[j] =
                pa_x[j] * tdx_xzz_xxzz_0[j] + 0.5 * fl1_fx * tdx_zz_xxzz_0[j] + fl1_fx * tdx_xzz_xzz_0[j] + 0.5 * fl1_fx * ts_xzz_xxzz_0[j];

            tdy_xxzz_xxzz_0[j] = pa_x[j] * tdy_xzz_xxzz_0[j] + 0.5 * fl1_fx * tdy_zz_xxzz_0[j] + fl1_fx * tdy_xzz_xzz_0[j];

            tdz_xxzz_xxzz_0[j] = pa_x[j] * tdz_xzz_xxzz_0[j] + 0.5 * fl1_fx * tdz_zz_xxzz_0[j] + fl1_fx * tdz_xzz_xzz_0[j];
        }

        idx++;
    }
}

void
compElectricDipoleForGG_243_291(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const int32_t              nOSFactors,
                                const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_d_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_d_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

        auto tdx_xzz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 81);

        auto tdy_xzz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 81);

        auto tdz_xzz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 81);

        auto tdx_xzz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 82);

        auto tdy_xzz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 82);

        auto tdz_xzz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 82);

        auto tdx_xzz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 83);

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

        auto tdx_xzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 56);

        auto tdy_xzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 56);

        auto tdz_xzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 56);

        auto tdx_xzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 57);

        auto tdy_xzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 57);

        auto tdz_xzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 57);

        auto tdx_xzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 58);

        auto tdy_xzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 58);

        auto tdz_xzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 58);

        auto tdx_xzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 59);

        auto tdy_xzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 59);

        auto tdz_xzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 59);

        auto tdx_yyy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 60);

        auto tdy_yyy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 60);

        auto tdz_yyy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 60);

        auto tdx_yyy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 61);

        auto tdy_yyy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 61);

        auto tdz_yyy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 61);

        auto tdx_yyy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 62);

        auto tdy_yyy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 62);

        auto tdz_yyy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 62);

        auto tdx_yyy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 63);

        auto tdy_yyy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 63);

        auto tdz_yyy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 63);

        auto tdx_yyy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 64);

        auto tdy_yyy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 64);

        auto tdz_yyy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 64);

        auto tdx_yyy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 65);

        auto tdy_yyy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 65);

        auto tdz_yyy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 65);

        auto tdx_yyy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 66);

        auto tdy_yyy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 66);

        auto tdz_yyy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 66);

        auto ts_xzz_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 81);

        auto ts_xzz_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 82);

        auto ts_xzz_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 83);

        auto ts_xzz_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 84);

        auto ts_xzz_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 85);

        auto ts_xzz_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 86);

        auto ts_xzz_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 87);

        auto ts_xzz_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 88);

        auto ts_xzz_zzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 89);

        auto ts_yyy_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 90);

        auto ts_yyy_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 91);

        auto ts_yyy_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 92);

        auto ts_yyy_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 93);

        auto ts_yyy_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 94);

        auto ts_yyy_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 95);

        auto ts_yyy_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 96);

        // set up pointers to integrals

        auto tdx_xxzz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 81);

        auto tdy_xxzz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 81);

        auto tdz_xxzz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 81);

        auto tdx_xxzz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 82);

        auto tdy_xxzz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 82);

        auto tdz_xxzz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 82);

        auto tdx_xxzz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 83);

        auto tdy_xxzz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 83);

        auto tdz_xxzz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 83);

        auto tdx_xxzz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 84);

        auto tdy_xxzz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 84);

        auto tdz_xxzz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 84);

        auto tdx_xxzz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 85);

        auto tdy_xxzz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 85);

        auto tdz_xxzz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 85);

        auto tdx_xxzz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 86);

        auto tdy_xxzz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 86);

        auto tdz_xxzz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 86);

        auto tdx_xxzz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 87);

        auto tdy_xxzz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 87);

        auto tdz_xxzz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 87);

        auto tdx_xxzz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 88);

        auto tdy_xxzz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 88);

        auto tdz_xxzz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 88);

        auto tdx_xxzz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 89);

        auto tdy_xxzz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 89);

        auto tdz_xxzz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 89);

        auto tdx_xyyy_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 90);

        auto tdy_xyyy_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 90);

        auto tdz_xyyy_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 90);

        auto tdx_xyyy_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 91);

        auto tdy_xyyy_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 91);

        auto tdz_xyyy_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 91);

        auto tdx_xyyy_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 92);

        auto tdy_xyyy_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 92);

        auto tdz_xyyy_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 92);

        auto tdx_xyyy_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 93);

        auto tdy_xyyy_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 93);

        auto tdz_xyyy_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 93);

        auto tdx_xyyy_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 94);

        auto tdy_xyyy_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 94);

        auto tdz_xyyy_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 94);

        auto tdx_xyyy_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 95);

        auto tdy_xyyy_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 95);

        auto tdz_xyyy_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 95);

        auto tdx_xyyy_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 96);

        auto tdy_xyyy_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 96);

        auto tdz_xyyy_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 96);

        // Batch of Integrals (243,291)

        #pragma omp simd aligned(fx, pa_x, tdx_xxzz_xyyy_0, tdx_xxzz_xyyz_0, tdx_xxzz_xyzz_0, \
                                     tdx_xxzz_xzzz_0, tdx_xxzz_yyyy_0, tdx_xxzz_yyyz_0, tdx_xxzz_yyzz_0, tdx_xxzz_yzzz_0, \
                                     tdx_xxzz_zzzz_0, tdx_xyyy_xxxx_0, tdx_xyyy_xxxy_0, tdx_xyyy_xxxz_0, tdx_xyyy_xxyy_0, \
                                     tdx_xyyy_xxyz_0, tdx_xyyy_xxzz_0, tdx_xyyy_xyyy_0, tdx_xzz_xyyy_0, tdx_xzz_xyyz_0, \
                                     tdx_xzz_xyzz_0, tdx_xzz_xzzz_0, tdx_xzz_yyy_0, tdx_xzz_yyyy_0, tdx_xzz_yyyz_0, \
                                     tdx_xzz_yyz_0, tdx_xzz_yyzz_0, tdx_xzz_yzz_0, tdx_xzz_yzzz_0, tdx_xzz_zzz_0, \
                                     tdx_xzz_zzzz_0, tdx_yyy_xxx_0, tdx_yyy_xxxx_0, tdx_yyy_xxxy_0, tdx_yyy_xxxz_0, \
                                     tdx_yyy_xxy_0, tdx_yyy_xxyy_0, tdx_yyy_xxyz_0, tdx_yyy_xxz_0, tdx_yyy_xxzz_0, \
                                     tdx_yyy_xyy_0, tdx_yyy_xyyy_0, tdx_yyy_xyz_0, tdx_yyy_xzz_0, tdx_yyy_yyy_0, \
                                     tdx_zz_xyyy_0, tdx_zz_xyyz_0, tdx_zz_xyzz_0, tdx_zz_xzzz_0, tdx_zz_yyyy_0, \
                                     tdx_zz_yyyz_0, tdx_zz_yyzz_0, tdx_zz_yzzz_0, tdx_zz_zzzz_0, tdy_xxzz_xyyy_0, \
                                     tdy_xxzz_xyyz_0, tdy_xxzz_xyzz_0, tdy_xxzz_xzzz_0, tdy_xxzz_yyyy_0, tdy_xxzz_yyyz_0, \
                                     tdy_xxzz_yyzz_0, tdy_xxzz_yzzz_0, tdy_xxzz_zzzz_0, tdy_xyyy_xxxx_0, tdy_xyyy_xxxy_0, \
                                     tdy_xyyy_xxxz_0, tdy_xyyy_xxyy_0, tdy_xyyy_xxyz_0, tdy_xyyy_xxzz_0, tdy_xyyy_xyyy_0, \
                                     tdy_xzz_xyyy_0, tdy_xzz_xyyz_0, tdy_xzz_xyzz_0, tdy_xzz_xzzz_0, tdy_xzz_yyy_0, \
                                     tdy_xzz_yyyy_0, tdy_xzz_yyyz_0, tdy_xzz_yyz_0, tdy_xzz_yyzz_0, tdy_xzz_yzz_0, \
                                     tdy_xzz_yzzz_0, tdy_xzz_zzz_0, tdy_xzz_zzzz_0, tdy_yyy_xxx_0, tdy_yyy_xxxx_0, \
                                     tdy_yyy_xxxy_0, tdy_yyy_xxxz_0, tdy_yyy_xxy_0, tdy_yyy_xxyy_0, tdy_yyy_xxyz_0, \
                                     tdy_yyy_xxz_0, tdy_yyy_xxzz_0, tdy_yyy_xyy_0, tdy_yyy_xyyy_0, tdy_yyy_xyz_0, \
                                     tdy_yyy_xzz_0, tdy_yyy_yyy_0, tdy_zz_xyyy_0, tdy_zz_xyyz_0, tdy_zz_xyzz_0, \
                                     tdy_zz_xzzz_0, tdy_zz_yyyy_0, tdy_zz_yyyz_0, tdy_zz_yyzz_0, tdy_zz_yzzz_0, \
                                     tdy_zz_zzzz_0, tdz_xxzz_xyyy_0, tdz_xxzz_xyyz_0, tdz_xxzz_xyzz_0, tdz_xxzz_xzzz_0, \
                                     tdz_xxzz_yyyy_0, tdz_xxzz_yyyz_0, tdz_xxzz_yyzz_0, tdz_xxzz_yzzz_0, tdz_xxzz_zzzz_0, \
                                     tdz_xyyy_xxxx_0, tdz_xyyy_xxxy_0, tdz_xyyy_xxxz_0, tdz_xyyy_xxyy_0, tdz_xyyy_xxyz_0, \
                                     tdz_xyyy_xxzz_0, tdz_xyyy_xyyy_0, tdz_xzz_xyyy_0, tdz_xzz_xyyz_0, tdz_xzz_xyzz_0, \
                                     tdz_xzz_xzzz_0, tdz_xzz_yyy_0, tdz_xzz_yyyy_0, tdz_xzz_yyyz_0, tdz_xzz_yyz_0, \
                                     tdz_xzz_yyzz_0, tdz_xzz_yzz_0, tdz_xzz_yzzz_0, tdz_xzz_zzz_0, tdz_xzz_zzzz_0, \
                                     tdz_yyy_xxx_0, tdz_yyy_xxxx_0, tdz_yyy_xxxy_0, tdz_yyy_xxxz_0, tdz_yyy_xxy_0, \
                                     tdz_yyy_xxyy_0, tdz_yyy_xxyz_0, tdz_yyy_xxz_0, tdz_yyy_xxzz_0, tdz_yyy_xyy_0, \
                                     tdz_yyy_xyyy_0, tdz_yyy_xyz_0, tdz_yyy_xzz_0, tdz_yyy_yyy_0, tdz_zz_xyyy_0, \
                                     tdz_zz_xyyz_0, tdz_zz_xyzz_0, tdz_zz_xzzz_0, tdz_zz_yyyy_0, tdz_zz_yyyz_0, \
                                     tdz_zz_yyzz_0, tdz_zz_yzzz_0, tdz_zz_zzzz_0, ts_xzz_xyyy_0, ts_xzz_xyyz_0, \
                                     ts_xzz_xyzz_0, ts_xzz_xzzz_0, ts_xzz_yyyy_0, ts_xzz_yyyz_0, ts_xzz_yyzz_0, \
                                     ts_xzz_yzzz_0, ts_xzz_zzzz_0, ts_yyy_xxxx_0, ts_yyy_xxxy_0, ts_yyy_xxxz_0, \
                                     ts_yyy_xxyy_0, ts_yyy_xxyz_0, ts_yyy_xxzz_0, ts_yyy_xyyy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            tdx_xxzz_xyyy_0[j] =
                pa_x[j] * tdx_xzz_xyyy_0[j] + 0.5 * fl1_fx * tdx_zz_xyyy_0[j] + 0.5 * fl1_fx * tdx_xzz_yyy_0[j] + 0.5 * fl1_fx * ts_xzz_xyyy_0[j];

            tdy_xxzz_xyyy_0[j] = pa_x[j] * tdy_xzz_xyyy_0[j] + 0.5 * fl1_fx * tdy_zz_xyyy_0[j] + 0.5 * fl1_fx * tdy_xzz_yyy_0[j];

            tdz_xxzz_xyyy_0[j] = pa_x[j] * tdz_xzz_xyyy_0[j] + 0.5 * fl1_fx * tdz_zz_xyyy_0[j] + 0.5 * fl1_fx * tdz_xzz_yyy_0[j];

            tdx_xxzz_xyyz_0[j] =
                pa_x[j] * tdx_xzz_xyyz_0[j] + 0.5 * fl1_fx * tdx_zz_xyyz_0[j] + 0.5 * fl1_fx * tdx_xzz_yyz_0[j] + 0.5 * fl1_fx * ts_xzz_xyyz_0[j];

            tdy_xxzz_xyyz_0[j] = pa_x[j] * tdy_xzz_xyyz_0[j] + 0.5 * fl1_fx * tdy_zz_xyyz_0[j] + 0.5 * fl1_fx * tdy_xzz_yyz_0[j];

            tdz_xxzz_xyyz_0[j] = pa_x[j] * tdz_xzz_xyyz_0[j] + 0.5 * fl1_fx * tdz_zz_xyyz_0[j] + 0.5 * fl1_fx * tdz_xzz_yyz_0[j];

            tdx_xxzz_xyzz_0[j] =
                pa_x[j] * tdx_xzz_xyzz_0[j] + 0.5 * fl1_fx * tdx_zz_xyzz_0[j] + 0.5 * fl1_fx * tdx_xzz_yzz_0[j] + 0.5 * fl1_fx * ts_xzz_xyzz_0[j];

            tdy_xxzz_xyzz_0[j] = pa_x[j] * tdy_xzz_xyzz_0[j] + 0.5 * fl1_fx * tdy_zz_xyzz_0[j] + 0.5 * fl1_fx * tdy_xzz_yzz_0[j];

            tdz_xxzz_xyzz_0[j] = pa_x[j] * tdz_xzz_xyzz_0[j] + 0.5 * fl1_fx * tdz_zz_xyzz_0[j] + 0.5 * fl1_fx * tdz_xzz_yzz_0[j];

            tdx_xxzz_xzzz_0[j] =
                pa_x[j] * tdx_xzz_xzzz_0[j] + 0.5 * fl1_fx * tdx_zz_xzzz_0[j] + 0.5 * fl1_fx * tdx_xzz_zzz_0[j] + 0.5 * fl1_fx * ts_xzz_xzzz_0[j];

            tdy_xxzz_xzzz_0[j] = pa_x[j] * tdy_xzz_xzzz_0[j] + 0.5 * fl1_fx * tdy_zz_xzzz_0[j] + 0.5 * fl1_fx * tdy_xzz_zzz_0[j];

            tdz_xxzz_xzzz_0[j] = pa_x[j] * tdz_xzz_xzzz_0[j] + 0.5 * fl1_fx * tdz_zz_xzzz_0[j] + 0.5 * fl1_fx * tdz_xzz_zzz_0[j];

            tdx_xxzz_yyyy_0[j] = pa_x[j] * tdx_xzz_yyyy_0[j] + 0.5 * fl1_fx * tdx_zz_yyyy_0[j] + 0.5 * fl1_fx * ts_xzz_yyyy_0[j];

            tdy_xxzz_yyyy_0[j] = pa_x[j] * tdy_xzz_yyyy_0[j] + 0.5 * fl1_fx * tdy_zz_yyyy_0[j];

            tdz_xxzz_yyyy_0[j] = pa_x[j] * tdz_xzz_yyyy_0[j] + 0.5 * fl1_fx * tdz_zz_yyyy_0[j];

            tdx_xxzz_yyyz_0[j] = pa_x[j] * tdx_xzz_yyyz_0[j] + 0.5 * fl1_fx * tdx_zz_yyyz_0[j] + 0.5 * fl1_fx * ts_xzz_yyyz_0[j];

            tdy_xxzz_yyyz_0[j] = pa_x[j] * tdy_xzz_yyyz_0[j] + 0.5 * fl1_fx * tdy_zz_yyyz_0[j];

            tdz_xxzz_yyyz_0[j] = pa_x[j] * tdz_xzz_yyyz_0[j] + 0.5 * fl1_fx * tdz_zz_yyyz_0[j];

            tdx_xxzz_yyzz_0[j] = pa_x[j] * tdx_xzz_yyzz_0[j] + 0.5 * fl1_fx * tdx_zz_yyzz_0[j] + 0.5 * fl1_fx * ts_xzz_yyzz_0[j];

            tdy_xxzz_yyzz_0[j] = pa_x[j] * tdy_xzz_yyzz_0[j] + 0.5 * fl1_fx * tdy_zz_yyzz_0[j];

            tdz_xxzz_yyzz_0[j] = pa_x[j] * tdz_xzz_yyzz_0[j] + 0.5 * fl1_fx * tdz_zz_yyzz_0[j];

            tdx_xxzz_yzzz_0[j] = pa_x[j] * tdx_xzz_yzzz_0[j] + 0.5 * fl1_fx * tdx_zz_yzzz_0[j] + 0.5 * fl1_fx * ts_xzz_yzzz_0[j];

            tdy_xxzz_yzzz_0[j] = pa_x[j] * tdy_xzz_yzzz_0[j] + 0.5 * fl1_fx * tdy_zz_yzzz_0[j];

            tdz_xxzz_yzzz_0[j] = pa_x[j] * tdz_xzz_yzzz_0[j] + 0.5 * fl1_fx * tdz_zz_yzzz_0[j];

            tdx_xxzz_zzzz_0[j] = pa_x[j] * tdx_xzz_zzzz_0[j] + 0.5 * fl1_fx * tdx_zz_zzzz_0[j] + 0.5 * fl1_fx * ts_xzz_zzzz_0[j];

            tdy_xxzz_zzzz_0[j] = pa_x[j] * tdy_xzz_zzzz_0[j] + 0.5 * fl1_fx * tdy_zz_zzzz_0[j];

            tdz_xxzz_zzzz_0[j] = pa_x[j] * tdz_xzz_zzzz_0[j] + 0.5 * fl1_fx * tdz_zz_zzzz_0[j];

            tdx_xyyy_xxxx_0[j] = pa_x[j] * tdx_yyy_xxxx_0[j] + 2.0 * fl1_fx * tdx_yyy_xxx_0[j] + 0.5 * fl1_fx * ts_yyy_xxxx_0[j];

            tdy_xyyy_xxxx_0[j] = pa_x[j] * tdy_yyy_xxxx_0[j] + 2.0 * fl1_fx * tdy_yyy_xxx_0[j];

            tdz_xyyy_xxxx_0[j] = pa_x[j] * tdz_yyy_xxxx_0[j] + 2.0 * fl1_fx * tdz_yyy_xxx_0[j];

            tdx_xyyy_xxxy_0[j] = pa_x[j] * tdx_yyy_xxxy_0[j] + 1.5 * fl1_fx * tdx_yyy_xxy_0[j] + 0.5 * fl1_fx * ts_yyy_xxxy_0[j];

            tdy_xyyy_xxxy_0[j] = pa_x[j] * tdy_yyy_xxxy_0[j] + 1.5 * fl1_fx * tdy_yyy_xxy_0[j];

            tdz_xyyy_xxxy_0[j] = pa_x[j] * tdz_yyy_xxxy_0[j] + 1.5 * fl1_fx * tdz_yyy_xxy_0[j];

            tdx_xyyy_xxxz_0[j] = pa_x[j] * tdx_yyy_xxxz_0[j] + 1.5 * fl1_fx * tdx_yyy_xxz_0[j] + 0.5 * fl1_fx * ts_yyy_xxxz_0[j];

            tdy_xyyy_xxxz_0[j] = pa_x[j] * tdy_yyy_xxxz_0[j] + 1.5 * fl1_fx * tdy_yyy_xxz_0[j];

            tdz_xyyy_xxxz_0[j] = pa_x[j] * tdz_yyy_xxxz_0[j] + 1.5 * fl1_fx * tdz_yyy_xxz_0[j];

            tdx_xyyy_xxyy_0[j] = pa_x[j] * tdx_yyy_xxyy_0[j] + fl1_fx * tdx_yyy_xyy_0[j] + 0.5 * fl1_fx * ts_yyy_xxyy_0[j];

            tdy_xyyy_xxyy_0[j] = pa_x[j] * tdy_yyy_xxyy_0[j] + fl1_fx * tdy_yyy_xyy_0[j];

            tdz_xyyy_xxyy_0[j] = pa_x[j] * tdz_yyy_xxyy_0[j] + fl1_fx * tdz_yyy_xyy_0[j];

            tdx_xyyy_xxyz_0[j] = pa_x[j] * tdx_yyy_xxyz_0[j] + fl1_fx * tdx_yyy_xyz_0[j] + 0.5 * fl1_fx * ts_yyy_xxyz_0[j];

            tdy_xyyy_xxyz_0[j] = pa_x[j] * tdy_yyy_xxyz_0[j] + fl1_fx * tdy_yyy_xyz_0[j];

            tdz_xyyy_xxyz_0[j] = pa_x[j] * tdz_yyy_xxyz_0[j] + fl1_fx * tdz_yyy_xyz_0[j];

            tdx_xyyy_xxzz_0[j] = pa_x[j] * tdx_yyy_xxzz_0[j] + fl1_fx * tdx_yyy_xzz_0[j] + 0.5 * fl1_fx * ts_yyy_xxzz_0[j];

            tdy_xyyy_xxzz_0[j] = pa_x[j] * tdy_yyy_xxzz_0[j] + fl1_fx * tdy_yyy_xzz_0[j];

            tdz_xyyy_xxzz_0[j] = pa_x[j] * tdz_yyy_xxzz_0[j] + fl1_fx * tdz_yyy_xzz_0[j];

            tdx_xyyy_xyyy_0[j] = pa_x[j] * tdx_yyy_xyyy_0[j] + 0.5 * fl1_fx * tdx_yyy_yyy_0[j] + 0.5 * fl1_fx * ts_yyy_xyyy_0[j];

            tdy_xyyy_xyyy_0[j] = pa_x[j] * tdy_yyy_xyyy_0[j] + 0.5 * fl1_fx * tdy_yyy_yyy_0[j];

            tdz_xyyy_xyyy_0[j] = pa_x[j] * tdz_yyy_xyyy_0[j] + 0.5 * fl1_fx * tdz_yyy_yyy_0[j];
        }

        idx++;
    }
}

void
compElectricDipoleForGG_291_339(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const int32_t              nOSFactors,
                                const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_d_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_d_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

        auto tdx_yyy_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 97);

        auto tdy_yyy_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 97);

        auto tdz_yyy_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 97);

        auto tdx_yyy_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 98);

        auto tdy_yyy_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 98);

        auto tdz_yyy_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 98);

        auto tdx_yyy_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 99);

        auto tdy_yyy_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 99);

        auto tdz_yyy_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 99);

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

        auto tdx_yyy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 67);

        auto tdy_yyy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 67);

        auto tdz_yyy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 67);

        auto tdx_yyy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 68);

        auto tdy_yyy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 68);

        auto tdz_yyy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 68);

        auto tdx_yyy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 69);

        auto tdy_yyy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 69);

        auto tdz_yyy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 69);

        auto tdx_yyz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 70);

        auto tdy_yyz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 70);

        auto tdz_yyz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 70);

        auto tdx_yyz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 71);

        auto tdy_yyz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 71);

        auto tdz_yyz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 71);

        auto tdx_yyz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 72);

        auto tdy_yyz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 72);

        auto tdz_yyz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 72);

        auto tdx_yyz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 73);

        auto tdy_yyz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 73);

        auto tdz_yyz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 73);

        auto tdx_yyz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 74);

        auto tdy_yyz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 74);

        auto tdz_yyz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 74);

        auto tdx_yyz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 75);

        auto tdy_yyz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 75);

        auto tdz_yyz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 75);

        auto tdx_yyz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 76);

        auto tdy_yyz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 76);

        auto tdz_yyz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 76);

        auto tdx_yyz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 77);

        auto tdy_yyz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 77);

        auto tdz_yyz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 77);

        auto ts_yyy_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 97);

        auto ts_yyy_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 98);

        auto ts_yyy_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 99);

        auto ts_yyy_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 100);

        auto ts_yyy_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 101);

        auto ts_yyy_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 102);

        auto ts_yyy_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 103);

        auto ts_yyy_zzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 104);

        auto ts_yyz_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 105);

        auto ts_yyz_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 106);

        auto ts_yyz_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 107);

        auto ts_yyz_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 108);

        auto ts_yyz_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 109);

        auto ts_yyz_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 110);

        auto ts_yyz_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 111);

        auto ts_yyz_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 112);

        // set up pointers to integrals

        auto tdx_xyyy_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 97);

        auto tdy_xyyy_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 97);

        auto tdz_xyyy_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 97);

        auto tdx_xyyy_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 98);

        auto tdy_xyyy_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 98);

        auto tdz_xyyy_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 98);

        auto tdx_xyyy_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 99);

        auto tdy_xyyy_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 99);

        auto tdz_xyyy_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 99);

        auto tdx_xyyy_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 100);

        auto tdy_xyyy_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 100);

        auto tdz_xyyy_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 100);

        auto tdx_xyyy_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 101);

        auto tdy_xyyy_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 101);

        auto tdz_xyyy_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 101);

        auto tdx_xyyy_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 102);

        auto tdy_xyyy_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 102);

        auto tdz_xyyy_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 102);

        auto tdx_xyyy_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 103);

        auto tdy_xyyy_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 103);

        auto tdz_xyyy_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 103);

        auto tdx_xyyy_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 104);

        auto tdy_xyyy_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 104);

        auto tdz_xyyy_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 104);

        auto tdx_xyyz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 105);

        auto tdy_xyyz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 105);

        auto tdz_xyyz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 105);

        auto tdx_xyyz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 106);

        auto tdy_xyyz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 106);

        auto tdz_xyyz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 106);

        auto tdx_xyyz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 107);

        auto tdy_xyyz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 107);

        auto tdz_xyyz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 107);

        auto tdx_xyyz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 108);

        auto tdy_xyyz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 108);

        auto tdz_xyyz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 108);

        auto tdx_xyyz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 109);

        auto tdy_xyyz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 109);

        auto tdz_xyyz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 109);

        auto tdx_xyyz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 110);

        auto tdy_xyyz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 110);

        auto tdz_xyyz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 110);

        auto tdx_xyyz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 111);

        auto tdy_xyyz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 111);

        auto tdz_xyyz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 111);

        auto tdx_xyyz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 112);

        auto tdy_xyyz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 112);

        auto tdz_xyyz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 112);

        // Batch of Integrals (291,339)

        #pragma omp simd aligned(fx, pa_x, tdx_xyyy_xyyz_0, tdx_xyyy_xyzz_0, tdx_xyyy_xzzz_0, \
                                     tdx_xyyy_yyyy_0, tdx_xyyy_yyyz_0, tdx_xyyy_yyzz_0, tdx_xyyy_yzzz_0, tdx_xyyy_zzzz_0, \
                                     tdx_xyyz_xxxx_0, tdx_xyyz_xxxy_0, tdx_xyyz_xxxz_0, tdx_xyyz_xxyy_0, tdx_xyyz_xxyz_0, \
                                     tdx_xyyz_xxzz_0, tdx_xyyz_xyyy_0, tdx_xyyz_xyyz_0, tdx_yyy_xyyz_0, tdx_yyy_xyzz_0, \
                                     tdx_yyy_xzzz_0, tdx_yyy_yyyy_0, tdx_yyy_yyyz_0, tdx_yyy_yyz_0, tdx_yyy_yyzz_0, \
                                     tdx_yyy_yzz_0, tdx_yyy_yzzz_0, tdx_yyy_zzz_0, tdx_yyy_zzzz_0, tdx_yyz_xxx_0, \
                                     tdx_yyz_xxxx_0, tdx_yyz_xxxy_0, tdx_yyz_xxxz_0, tdx_yyz_xxy_0, tdx_yyz_xxyy_0, \
                                     tdx_yyz_xxyz_0, tdx_yyz_xxz_0, tdx_yyz_xxzz_0, tdx_yyz_xyy_0, tdx_yyz_xyyy_0, \
                                     tdx_yyz_xyyz_0, tdx_yyz_xyz_0, tdx_yyz_xzz_0, tdx_yyz_yyy_0, tdx_yyz_yyz_0, \
                                     tdy_xyyy_xyyz_0, tdy_xyyy_xyzz_0, tdy_xyyy_xzzz_0, tdy_xyyy_yyyy_0, tdy_xyyy_yyyz_0, \
                                     tdy_xyyy_yyzz_0, tdy_xyyy_yzzz_0, tdy_xyyy_zzzz_0, tdy_xyyz_xxxx_0, tdy_xyyz_xxxy_0, \
                                     tdy_xyyz_xxxz_0, tdy_xyyz_xxyy_0, tdy_xyyz_xxyz_0, tdy_xyyz_xxzz_0, tdy_xyyz_xyyy_0, \
                                     tdy_xyyz_xyyz_0, tdy_yyy_xyyz_0, tdy_yyy_xyzz_0, tdy_yyy_xzzz_0, tdy_yyy_yyyy_0, \
                                     tdy_yyy_yyyz_0, tdy_yyy_yyz_0, tdy_yyy_yyzz_0, tdy_yyy_yzz_0, tdy_yyy_yzzz_0, \
                                     tdy_yyy_zzz_0, tdy_yyy_zzzz_0, tdy_yyz_xxx_0, tdy_yyz_xxxx_0, tdy_yyz_xxxy_0, \
                                     tdy_yyz_xxxz_0, tdy_yyz_xxy_0, tdy_yyz_xxyy_0, tdy_yyz_xxyz_0, tdy_yyz_xxz_0, \
                                     tdy_yyz_xxzz_0, tdy_yyz_xyy_0, tdy_yyz_xyyy_0, tdy_yyz_xyyz_0, tdy_yyz_xyz_0, \
                                     tdy_yyz_xzz_0, tdy_yyz_yyy_0, tdy_yyz_yyz_0, tdz_xyyy_xyyz_0, tdz_xyyy_xyzz_0, \
                                     tdz_xyyy_xzzz_0, tdz_xyyy_yyyy_0, tdz_xyyy_yyyz_0, tdz_xyyy_yyzz_0, tdz_xyyy_yzzz_0, \
                                     tdz_xyyy_zzzz_0, tdz_xyyz_xxxx_0, tdz_xyyz_xxxy_0, tdz_xyyz_xxxz_0, tdz_xyyz_xxyy_0, \
                                     tdz_xyyz_xxyz_0, tdz_xyyz_xxzz_0, tdz_xyyz_xyyy_0, tdz_xyyz_xyyz_0, tdz_yyy_xyyz_0, \
                                     tdz_yyy_xyzz_0, tdz_yyy_xzzz_0, tdz_yyy_yyyy_0, tdz_yyy_yyyz_0, tdz_yyy_yyz_0, \
                                     tdz_yyy_yyzz_0, tdz_yyy_yzz_0, tdz_yyy_yzzz_0, tdz_yyy_zzz_0, tdz_yyy_zzzz_0, \
                                     tdz_yyz_xxx_0, tdz_yyz_xxxx_0, tdz_yyz_xxxy_0, tdz_yyz_xxxz_0, tdz_yyz_xxy_0, \
                                     tdz_yyz_xxyy_0, tdz_yyz_xxyz_0, tdz_yyz_xxz_0, tdz_yyz_xxzz_0, tdz_yyz_xyy_0, \
                                     tdz_yyz_xyyy_0, tdz_yyz_xyyz_0, tdz_yyz_xyz_0, tdz_yyz_xzz_0, tdz_yyz_yyy_0, \
                                     tdz_yyz_yyz_0, ts_yyy_xyyz_0, ts_yyy_xyzz_0, ts_yyy_xzzz_0, ts_yyy_yyyy_0, \
                                     ts_yyy_yyyz_0, ts_yyy_yyzz_0, ts_yyy_yzzz_0, ts_yyy_zzzz_0, ts_yyz_xxxx_0, \
                                     ts_yyz_xxxy_0, ts_yyz_xxxz_0, ts_yyz_xxyy_0, ts_yyz_xxyz_0, ts_yyz_xxzz_0, \
                                     ts_yyz_xyyy_0, ts_yyz_xyyz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            tdx_xyyy_xyyz_0[j] = pa_x[j] * tdx_yyy_xyyz_0[j] + 0.5 * fl1_fx * tdx_yyy_yyz_0[j] + 0.5 * fl1_fx * ts_yyy_xyyz_0[j];

            tdy_xyyy_xyyz_0[j] = pa_x[j] * tdy_yyy_xyyz_0[j] + 0.5 * fl1_fx * tdy_yyy_yyz_0[j];

            tdz_xyyy_xyyz_0[j] = pa_x[j] * tdz_yyy_xyyz_0[j] + 0.5 * fl1_fx * tdz_yyy_yyz_0[j];

            tdx_xyyy_xyzz_0[j] = pa_x[j] * tdx_yyy_xyzz_0[j] + 0.5 * fl1_fx * tdx_yyy_yzz_0[j] + 0.5 * fl1_fx * ts_yyy_xyzz_0[j];

            tdy_xyyy_xyzz_0[j] = pa_x[j] * tdy_yyy_xyzz_0[j] + 0.5 * fl1_fx * tdy_yyy_yzz_0[j];

            tdz_xyyy_xyzz_0[j] = pa_x[j] * tdz_yyy_xyzz_0[j] + 0.5 * fl1_fx * tdz_yyy_yzz_0[j];

            tdx_xyyy_xzzz_0[j] = pa_x[j] * tdx_yyy_xzzz_0[j] + 0.5 * fl1_fx * tdx_yyy_zzz_0[j] + 0.5 * fl1_fx * ts_yyy_xzzz_0[j];

            tdy_xyyy_xzzz_0[j] = pa_x[j] * tdy_yyy_xzzz_0[j] + 0.5 * fl1_fx * tdy_yyy_zzz_0[j];

            tdz_xyyy_xzzz_0[j] = pa_x[j] * tdz_yyy_xzzz_0[j] + 0.5 * fl1_fx * tdz_yyy_zzz_0[j];

            tdx_xyyy_yyyy_0[j] = pa_x[j] * tdx_yyy_yyyy_0[j] + 0.5 * fl1_fx * ts_yyy_yyyy_0[j];

            tdy_xyyy_yyyy_0[j] = pa_x[j] * tdy_yyy_yyyy_0[j];

            tdz_xyyy_yyyy_0[j] = pa_x[j] * tdz_yyy_yyyy_0[j];

            tdx_xyyy_yyyz_0[j] = pa_x[j] * tdx_yyy_yyyz_0[j] + 0.5 * fl1_fx * ts_yyy_yyyz_0[j];

            tdy_xyyy_yyyz_0[j] = pa_x[j] * tdy_yyy_yyyz_0[j];

            tdz_xyyy_yyyz_0[j] = pa_x[j] * tdz_yyy_yyyz_0[j];

            tdx_xyyy_yyzz_0[j] = pa_x[j] * tdx_yyy_yyzz_0[j] + 0.5 * fl1_fx * ts_yyy_yyzz_0[j];

            tdy_xyyy_yyzz_0[j] = pa_x[j] * tdy_yyy_yyzz_0[j];

            tdz_xyyy_yyzz_0[j] = pa_x[j] * tdz_yyy_yyzz_0[j];

            tdx_xyyy_yzzz_0[j] = pa_x[j] * tdx_yyy_yzzz_0[j] + 0.5 * fl1_fx * ts_yyy_yzzz_0[j];

            tdy_xyyy_yzzz_0[j] = pa_x[j] * tdy_yyy_yzzz_0[j];

            tdz_xyyy_yzzz_0[j] = pa_x[j] * tdz_yyy_yzzz_0[j];

            tdx_xyyy_zzzz_0[j] = pa_x[j] * tdx_yyy_zzzz_0[j] + 0.5 * fl1_fx * ts_yyy_zzzz_0[j];

            tdy_xyyy_zzzz_0[j] = pa_x[j] * tdy_yyy_zzzz_0[j];

            tdz_xyyy_zzzz_0[j] = pa_x[j] * tdz_yyy_zzzz_0[j];

            tdx_xyyz_xxxx_0[j] = pa_x[j] * tdx_yyz_xxxx_0[j] + 2.0 * fl1_fx * tdx_yyz_xxx_0[j] + 0.5 * fl1_fx * ts_yyz_xxxx_0[j];

            tdy_xyyz_xxxx_0[j] = pa_x[j] * tdy_yyz_xxxx_0[j] + 2.0 * fl1_fx * tdy_yyz_xxx_0[j];

            tdz_xyyz_xxxx_0[j] = pa_x[j] * tdz_yyz_xxxx_0[j] + 2.0 * fl1_fx * tdz_yyz_xxx_0[j];

            tdx_xyyz_xxxy_0[j] = pa_x[j] * tdx_yyz_xxxy_0[j] + 1.5 * fl1_fx * tdx_yyz_xxy_0[j] + 0.5 * fl1_fx * ts_yyz_xxxy_0[j];

            tdy_xyyz_xxxy_0[j] = pa_x[j] * tdy_yyz_xxxy_0[j] + 1.5 * fl1_fx * tdy_yyz_xxy_0[j];

            tdz_xyyz_xxxy_0[j] = pa_x[j] * tdz_yyz_xxxy_0[j] + 1.5 * fl1_fx * tdz_yyz_xxy_0[j];

            tdx_xyyz_xxxz_0[j] = pa_x[j] * tdx_yyz_xxxz_0[j] + 1.5 * fl1_fx * tdx_yyz_xxz_0[j] + 0.5 * fl1_fx * ts_yyz_xxxz_0[j];

            tdy_xyyz_xxxz_0[j] = pa_x[j] * tdy_yyz_xxxz_0[j] + 1.5 * fl1_fx * tdy_yyz_xxz_0[j];

            tdz_xyyz_xxxz_0[j] = pa_x[j] * tdz_yyz_xxxz_0[j] + 1.5 * fl1_fx * tdz_yyz_xxz_0[j];

            tdx_xyyz_xxyy_0[j] = pa_x[j] * tdx_yyz_xxyy_0[j] + fl1_fx * tdx_yyz_xyy_0[j] + 0.5 * fl1_fx * ts_yyz_xxyy_0[j];

            tdy_xyyz_xxyy_0[j] = pa_x[j] * tdy_yyz_xxyy_0[j] + fl1_fx * tdy_yyz_xyy_0[j];

            tdz_xyyz_xxyy_0[j] = pa_x[j] * tdz_yyz_xxyy_0[j] + fl1_fx * tdz_yyz_xyy_0[j];

            tdx_xyyz_xxyz_0[j] = pa_x[j] * tdx_yyz_xxyz_0[j] + fl1_fx * tdx_yyz_xyz_0[j] + 0.5 * fl1_fx * ts_yyz_xxyz_0[j];

            tdy_xyyz_xxyz_0[j] = pa_x[j] * tdy_yyz_xxyz_0[j] + fl1_fx * tdy_yyz_xyz_0[j];

            tdz_xyyz_xxyz_0[j] = pa_x[j] * tdz_yyz_xxyz_0[j] + fl1_fx * tdz_yyz_xyz_0[j];

            tdx_xyyz_xxzz_0[j] = pa_x[j] * tdx_yyz_xxzz_0[j] + fl1_fx * tdx_yyz_xzz_0[j] + 0.5 * fl1_fx * ts_yyz_xxzz_0[j];

            tdy_xyyz_xxzz_0[j] = pa_x[j] * tdy_yyz_xxzz_0[j] + fl1_fx * tdy_yyz_xzz_0[j];

            tdz_xyyz_xxzz_0[j] = pa_x[j] * tdz_yyz_xxzz_0[j] + fl1_fx * tdz_yyz_xzz_0[j];

            tdx_xyyz_xyyy_0[j] = pa_x[j] * tdx_yyz_xyyy_0[j] + 0.5 * fl1_fx * tdx_yyz_yyy_0[j] + 0.5 * fl1_fx * ts_yyz_xyyy_0[j];

            tdy_xyyz_xyyy_0[j] = pa_x[j] * tdy_yyz_xyyy_0[j] + 0.5 * fl1_fx * tdy_yyz_yyy_0[j];

            tdz_xyyz_xyyy_0[j] = pa_x[j] * tdz_yyz_xyyy_0[j] + 0.5 * fl1_fx * tdz_yyz_yyy_0[j];

            tdx_xyyz_xyyz_0[j] = pa_x[j] * tdx_yyz_xyyz_0[j] + 0.5 * fl1_fx * tdx_yyz_yyz_0[j] + 0.5 * fl1_fx * ts_yyz_xyyz_0[j];

            tdy_xyyz_xyyz_0[j] = pa_x[j] * tdy_yyz_xyyz_0[j] + 0.5 * fl1_fx * tdy_yyz_yyz_0[j];

            tdz_xyyz_xyyz_0[j] = pa_x[j] * tdz_yyz_xyyz_0[j] + 0.5 * fl1_fx * tdz_yyz_yyz_0[j];
        }

        idx++;
    }
}

void
compElectricDipoleForGG_339_387(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const int32_t              nOSFactors,
                                const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_d_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_d_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

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

        auto tdx_yyz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 78);

        auto tdy_yyz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 78);

        auto tdz_yyz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 78);

        auto tdx_yyz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 79);

        auto tdy_yyz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 79);

        auto tdz_yyz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 79);

        auto tdx_yzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 80);

        auto tdy_yzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 80);

        auto tdz_yzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 80);

        auto tdx_yzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 81);

        auto tdy_yzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 81);

        auto tdz_yzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 81);

        auto tdx_yzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 82);

        auto tdy_yzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 82);

        auto tdz_yzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 82);

        auto tdx_yzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 83);

        auto tdy_yzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 83);

        auto tdz_yzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 83);

        auto tdx_yzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 84);

        auto tdy_yzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 84);

        auto tdz_yzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 84);

        auto tdx_yzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 85);

        auto tdy_yzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 85);

        auto tdz_yzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 85);

        auto tdx_yzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 86);

        auto tdy_yzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 86);

        auto tdz_yzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 86);

        auto tdx_yzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 87);

        auto tdy_yzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 87);

        auto tdz_yzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 87);

        auto tdx_yzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 88);

        auto tdy_yzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 88);

        auto tdz_yzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 88);

        auto ts_yyz_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 113);

        auto ts_yyz_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 114);

        auto ts_yyz_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 115);

        auto ts_yyz_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 116);

        auto ts_yyz_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 117);

        auto ts_yyz_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 118);

        auto ts_yyz_zzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 119);

        auto ts_yzz_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 120);

        auto ts_yzz_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 121);

        auto ts_yzz_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 122);

        auto ts_yzz_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 123);

        auto ts_yzz_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 124);

        auto ts_yzz_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 125);

        auto ts_yzz_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 126);

        auto ts_yzz_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 127);

        auto ts_yzz_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 128);

        // set up pointers to integrals

        auto tdx_xyyz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 113);

        auto tdy_xyyz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 113);

        auto tdz_xyyz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 113);

        auto tdx_xyyz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 114);

        auto tdy_xyyz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 114);

        auto tdz_xyyz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 114);

        auto tdx_xyyz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 115);

        auto tdy_xyyz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 115);

        auto tdz_xyyz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 115);

        auto tdx_xyyz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 116);

        auto tdy_xyyz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 116);

        auto tdz_xyyz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 116);

        auto tdx_xyyz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 117);

        auto tdy_xyyz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 117);

        auto tdz_xyyz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 117);

        auto tdx_xyyz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 118);

        auto tdy_xyyz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 118);

        auto tdz_xyyz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 118);

        auto tdx_xyyz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 119);

        auto tdy_xyyz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 119);

        auto tdz_xyyz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 119);

        auto tdx_xyzz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 120);

        auto tdy_xyzz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 120);

        auto tdz_xyzz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 120);

        auto tdx_xyzz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 121);

        auto tdy_xyzz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 121);

        auto tdz_xyzz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 121);

        auto tdx_xyzz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 122);

        auto tdy_xyzz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 122);

        auto tdz_xyzz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 122);

        auto tdx_xyzz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 123);

        auto tdy_xyzz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 123);

        auto tdz_xyzz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 123);

        auto tdx_xyzz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 124);

        auto tdy_xyzz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 124);

        auto tdz_xyzz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 124);

        auto tdx_xyzz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 125);

        auto tdy_xyzz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 125);

        auto tdz_xyzz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 125);

        auto tdx_xyzz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 126);

        auto tdy_xyzz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 126);

        auto tdz_xyzz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 126);

        auto tdx_xyzz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 127);

        auto tdy_xyzz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 127);

        auto tdz_xyzz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 127);

        auto tdx_xyzz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 128);

        auto tdy_xyzz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 128);

        auto tdz_xyzz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 128);

        // Batch of Integrals (339,387)

        #pragma omp simd aligned(fx, pa_x, tdx_xyyz_xyzz_0, tdx_xyyz_xzzz_0, tdx_xyyz_yyyy_0, \
                                     tdx_xyyz_yyyz_0, tdx_xyyz_yyzz_0, tdx_xyyz_yzzz_0, tdx_xyyz_zzzz_0, tdx_xyzz_xxxx_0, \
                                     tdx_xyzz_xxxy_0, tdx_xyzz_xxxz_0, tdx_xyzz_xxyy_0, tdx_xyzz_xxyz_0, tdx_xyzz_xxzz_0, \
                                     tdx_xyzz_xyyy_0, tdx_xyzz_xyyz_0, tdx_xyzz_xyzz_0, tdx_yyz_xyzz_0, tdx_yyz_xzzz_0, \
                                     tdx_yyz_yyyy_0, tdx_yyz_yyyz_0, tdx_yyz_yyzz_0, tdx_yyz_yzz_0, tdx_yyz_yzzz_0, \
                                     tdx_yyz_zzz_0, tdx_yyz_zzzz_0, tdx_yzz_xxx_0, tdx_yzz_xxxx_0, tdx_yzz_xxxy_0, \
                                     tdx_yzz_xxxz_0, tdx_yzz_xxy_0, tdx_yzz_xxyy_0, tdx_yzz_xxyz_0, tdx_yzz_xxz_0, \
                                     tdx_yzz_xxzz_0, tdx_yzz_xyy_0, tdx_yzz_xyyy_0, tdx_yzz_xyyz_0, tdx_yzz_xyz_0, \
                                     tdx_yzz_xyzz_0, tdx_yzz_xzz_0, tdx_yzz_yyy_0, tdx_yzz_yyz_0, tdx_yzz_yzz_0, \
                                     tdy_xyyz_xyzz_0, tdy_xyyz_xzzz_0, tdy_xyyz_yyyy_0, tdy_xyyz_yyyz_0, tdy_xyyz_yyzz_0, \
                                     tdy_xyyz_yzzz_0, tdy_xyyz_zzzz_0, tdy_xyzz_xxxx_0, tdy_xyzz_xxxy_0, tdy_xyzz_xxxz_0, \
                                     tdy_xyzz_xxyy_0, tdy_xyzz_xxyz_0, tdy_xyzz_xxzz_0, tdy_xyzz_xyyy_0, tdy_xyzz_xyyz_0, \
                                     tdy_xyzz_xyzz_0, tdy_yyz_xyzz_0, tdy_yyz_xzzz_0, tdy_yyz_yyyy_0, tdy_yyz_yyyz_0, \
                                     tdy_yyz_yyzz_0, tdy_yyz_yzz_0, tdy_yyz_yzzz_0, tdy_yyz_zzz_0, tdy_yyz_zzzz_0, \
                                     tdy_yzz_xxx_0, tdy_yzz_xxxx_0, tdy_yzz_xxxy_0, tdy_yzz_xxxz_0, tdy_yzz_xxy_0, \
                                     tdy_yzz_xxyy_0, tdy_yzz_xxyz_0, tdy_yzz_xxz_0, tdy_yzz_xxzz_0, tdy_yzz_xyy_0, \
                                     tdy_yzz_xyyy_0, tdy_yzz_xyyz_0, tdy_yzz_xyz_0, tdy_yzz_xyzz_0, tdy_yzz_xzz_0, \
                                     tdy_yzz_yyy_0, tdy_yzz_yyz_0, tdy_yzz_yzz_0, tdz_xyyz_xyzz_0, tdz_xyyz_xzzz_0, \
                                     tdz_xyyz_yyyy_0, tdz_xyyz_yyyz_0, tdz_xyyz_yyzz_0, tdz_xyyz_yzzz_0, tdz_xyyz_zzzz_0, \
                                     tdz_xyzz_xxxx_0, tdz_xyzz_xxxy_0, tdz_xyzz_xxxz_0, tdz_xyzz_xxyy_0, tdz_xyzz_xxyz_0, \
                                     tdz_xyzz_xxzz_0, tdz_xyzz_xyyy_0, tdz_xyzz_xyyz_0, tdz_xyzz_xyzz_0, tdz_yyz_xyzz_0, \
                                     tdz_yyz_xzzz_0, tdz_yyz_yyyy_0, tdz_yyz_yyyz_0, tdz_yyz_yyzz_0, tdz_yyz_yzz_0, \
                                     tdz_yyz_yzzz_0, tdz_yyz_zzz_0, tdz_yyz_zzzz_0, tdz_yzz_xxx_0, tdz_yzz_xxxx_0, \
                                     tdz_yzz_xxxy_0, tdz_yzz_xxxz_0, tdz_yzz_xxy_0, tdz_yzz_xxyy_0, tdz_yzz_xxyz_0, \
                                     tdz_yzz_xxz_0, tdz_yzz_xxzz_0, tdz_yzz_xyy_0, tdz_yzz_xyyy_0, tdz_yzz_xyyz_0, \
                                     tdz_yzz_xyz_0, tdz_yzz_xyzz_0, tdz_yzz_xzz_0, tdz_yzz_yyy_0, tdz_yzz_yyz_0, \
                                     tdz_yzz_yzz_0, ts_yyz_xyzz_0, ts_yyz_xzzz_0, ts_yyz_yyyy_0, ts_yyz_yyyz_0, \
                                     ts_yyz_yyzz_0, ts_yyz_yzzz_0, ts_yyz_zzzz_0, ts_yzz_xxxx_0, ts_yzz_xxxy_0, \
                                     ts_yzz_xxxz_0, ts_yzz_xxyy_0, ts_yzz_xxyz_0, ts_yzz_xxzz_0, ts_yzz_xyyy_0, \
                                     ts_yzz_xyyz_0, ts_yzz_xyzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            tdx_xyyz_xyzz_0[j] = pa_x[j] * tdx_yyz_xyzz_0[j] + 0.5 * fl1_fx * tdx_yyz_yzz_0[j] + 0.5 * fl1_fx * ts_yyz_xyzz_0[j];

            tdy_xyyz_xyzz_0[j] = pa_x[j] * tdy_yyz_xyzz_0[j] + 0.5 * fl1_fx * tdy_yyz_yzz_0[j];

            tdz_xyyz_xyzz_0[j] = pa_x[j] * tdz_yyz_xyzz_0[j] + 0.5 * fl1_fx * tdz_yyz_yzz_0[j];

            tdx_xyyz_xzzz_0[j] = pa_x[j] * tdx_yyz_xzzz_0[j] + 0.5 * fl1_fx * tdx_yyz_zzz_0[j] + 0.5 * fl1_fx * ts_yyz_xzzz_0[j];

            tdy_xyyz_xzzz_0[j] = pa_x[j] * tdy_yyz_xzzz_0[j] + 0.5 * fl1_fx * tdy_yyz_zzz_0[j];

            tdz_xyyz_xzzz_0[j] = pa_x[j] * tdz_yyz_xzzz_0[j] + 0.5 * fl1_fx * tdz_yyz_zzz_0[j];

            tdx_xyyz_yyyy_0[j] = pa_x[j] * tdx_yyz_yyyy_0[j] + 0.5 * fl1_fx * ts_yyz_yyyy_0[j];

            tdy_xyyz_yyyy_0[j] = pa_x[j] * tdy_yyz_yyyy_0[j];

            tdz_xyyz_yyyy_0[j] = pa_x[j] * tdz_yyz_yyyy_0[j];

            tdx_xyyz_yyyz_0[j] = pa_x[j] * tdx_yyz_yyyz_0[j] + 0.5 * fl1_fx * ts_yyz_yyyz_0[j];

            tdy_xyyz_yyyz_0[j] = pa_x[j] * tdy_yyz_yyyz_0[j];

            tdz_xyyz_yyyz_0[j] = pa_x[j] * tdz_yyz_yyyz_0[j];

            tdx_xyyz_yyzz_0[j] = pa_x[j] * tdx_yyz_yyzz_0[j] + 0.5 * fl1_fx * ts_yyz_yyzz_0[j];

            tdy_xyyz_yyzz_0[j] = pa_x[j] * tdy_yyz_yyzz_0[j];

            tdz_xyyz_yyzz_0[j] = pa_x[j] * tdz_yyz_yyzz_0[j];

            tdx_xyyz_yzzz_0[j] = pa_x[j] * tdx_yyz_yzzz_0[j] + 0.5 * fl1_fx * ts_yyz_yzzz_0[j];

            tdy_xyyz_yzzz_0[j] = pa_x[j] * tdy_yyz_yzzz_0[j];

            tdz_xyyz_yzzz_0[j] = pa_x[j] * tdz_yyz_yzzz_0[j];

            tdx_xyyz_zzzz_0[j] = pa_x[j] * tdx_yyz_zzzz_0[j] + 0.5 * fl1_fx * ts_yyz_zzzz_0[j];

            tdy_xyyz_zzzz_0[j] = pa_x[j] * tdy_yyz_zzzz_0[j];

            tdz_xyyz_zzzz_0[j] = pa_x[j] * tdz_yyz_zzzz_0[j];

            tdx_xyzz_xxxx_0[j] = pa_x[j] * tdx_yzz_xxxx_0[j] + 2.0 * fl1_fx * tdx_yzz_xxx_0[j] + 0.5 * fl1_fx * ts_yzz_xxxx_0[j];

            tdy_xyzz_xxxx_0[j] = pa_x[j] * tdy_yzz_xxxx_0[j] + 2.0 * fl1_fx * tdy_yzz_xxx_0[j];

            tdz_xyzz_xxxx_0[j] = pa_x[j] * tdz_yzz_xxxx_0[j] + 2.0 * fl1_fx * tdz_yzz_xxx_0[j];

            tdx_xyzz_xxxy_0[j] = pa_x[j] * tdx_yzz_xxxy_0[j] + 1.5 * fl1_fx * tdx_yzz_xxy_0[j] + 0.5 * fl1_fx * ts_yzz_xxxy_0[j];

            tdy_xyzz_xxxy_0[j] = pa_x[j] * tdy_yzz_xxxy_0[j] + 1.5 * fl1_fx * tdy_yzz_xxy_0[j];

            tdz_xyzz_xxxy_0[j] = pa_x[j] * tdz_yzz_xxxy_0[j] + 1.5 * fl1_fx * tdz_yzz_xxy_0[j];

            tdx_xyzz_xxxz_0[j] = pa_x[j] * tdx_yzz_xxxz_0[j] + 1.5 * fl1_fx * tdx_yzz_xxz_0[j] + 0.5 * fl1_fx * ts_yzz_xxxz_0[j];

            tdy_xyzz_xxxz_0[j] = pa_x[j] * tdy_yzz_xxxz_0[j] + 1.5 * fl1_fx * tdy_yzz_xxz_0[j];

            tdz_xyzz_xxxz_0[j] = pa_x[j] * tdz_yzz_xxxz_0[j] + 1.5 * fl1_fx * tdz_yzz_xxz_0[j];

            tdx_xyzz_xxyy_0[j] = pa_x[j] * tdx_yzz_xxyy_0[j] + fl1_fx * tdx_yzz_xyy_0[j] + 0.5 * fl1_fx * ts_yzz_xxyy_0[j];

            tdy_xyzz_xxyy_0[j] = pa_x[j] * tdy_yzz_xxyy_0[j] + fl1_fx * tdy_yzz_xyy_0[j];

            tdz_xyzz_xxyy_0[j] = pa_x[j] * tdz_yzz_xxyy_0[j] + fl1_fx * tdz_yzz_xyy_0[j];

            tdx_xyzz_xxyz_0[j] = pa_x[j] * tdx_yzz_xxyz_0[j] + fl1_fx * tdx_yzz_xyz_0[j] + 0.5 * fl1_fx * ts_yzz_xxyz_0[j];

            tdy_xyzz_xxyz_0[j] = pa_x[j] * tdy_yzz_xxyz_0[j] + fl1_fx * tdy_yzz_xyz_0[j];

            tdz_xyzz_xxyz_0[j] = pa_x[j] * tdz_yzz_xxyz_0[j] + fl1_fx * tdz_yzz_xyz_0[j];

            tdx_xyzz_xxzz_0[j] = pa_x[j] * tdx_yzz_xxzz_0[j] + fl1_fx * tdx_yzz_xzz_0[j] + 0.5 * fl1_fx * ts_yzz_xxzz_0[j];

            tdy_xyzz_xxzz_0[j] = pa_x[j] * tdy_yzz_xxzz_0[j] + fl1_fx * tdy_yzz_xzz_0[j];

            tdz_xyzz_xxzz_0[j] = pa_x[j] * tdz_yzz_xxzz_0[j] + fl1_fx * tdz_yzz_xzz_0[j];

            tdx_xyzz_xyyy_0[j] = pa_x[j] * tdx_yzz_xyyy_0[j] + 0.5 * fl1_fx * tdx_yzz_yyy_0[j] + 0.5 * fl1_fx * ts_yzz_xyyy_0[j];

            tdy_xyzz_xyyy_0[j] = pa_x[j] * tdy_yzz_xyyy_0[j] + 0.5 * fl1_fx * tdy_yzz_yyy_0[j];

            tdz_xyzz_xyyy_0[j] = pa_x[j] * tdz_yzz_xyyy_0[j] + 0.5 * fl1_fx * tdz_yzz_yyy_0[j];

            tdx_xyzz_xyyz_0[j] = pa_x[j] * tdx_yzz_xyyz_0[j] + 0.5 * fl1_fx * tdx_yzz_yyz_0[j] + 0.5 * fl1_fx * ts_yzz_xyyz_0[j];

            tdy_xyzz_xyyz_0[j] = pa_x[j] * tdy_yzz_xyyz_0[j] + 0.5 * fl1_fx * tdy_yzz_yyz_0[j];

            tdz_xyzz_xyyz_0[j] = pa_x[j] * tdz_yzz_xyyz_0[j] + 0.5 * fl1_fx * tdz_yzz_yyz_0[j];

            tdx_xyzz_xyzz_0[j] = pa_x[j] * tdx_yzz_xyzz_0[j] + 0.5 * fl1_fx * tdx_yzz_yzz_0[j] + 0.5 * fl1_fx * ts_yzz_xyzz_0[j];

            tdy_xyzz_xyzz_0[j] = pa_x[j] * tdy_yzz_xyzz_0[j] + 0.5 * fl1_fx * tdy_yzz_yzz_0[j];

            tdz_xyzz_xyzz_0[j] = pa_x[j] * tdz_yzz_xyzz_0[j] + 0.5 * fl1_fx * tdz_yzz_yzz_0[j];
        }

        idx++;
    }
}

void
compElectricDipoleForGG_387_435(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const int32_t              nOSFactors,
                                const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_d_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_d_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

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

        auto tdx_yzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 89);

        auto tdy_yzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 89);

        auto tdz_yzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 89);

        auto tdx_zzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 90);

        auto tdy_zzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 90);

        auto tdz_zzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 90);

        auto tdx_zzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 91);

        auto tdy_zzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 91);

        auto tdz_zzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 91);

        auto tdx_zzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 92);

        auto tdy_zzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 92);

        auto tdz_zzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 92);

        auto tdx_zzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 93);

        auto tdy_zzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 93);

        auto tdz_zzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 93);

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

        auto ts_yzz_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 129);

        auto ts_yzz_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 130);

        auto ts_yzz_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 131);

        auto ts_yzz_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 132);

        auto ts_yzz_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 133);

        auto ts_yzz_zzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 134);

        auto ts_zzz_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 135);

        auto ts_zzz_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 136);

        auto ts_zzz_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 137);

        auto ts_zzz_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 138);

        auto ts_zzz_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 139);

        auto ts_zzz_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 140);

        auto ts_zzz_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 141);

        auto ts_zzz_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 142);

        auto ts_zzz_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 143);

        auto ts_zzz_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 144);

        // set up pointers to integrals

        auto tdx_xyzz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 129);

        auto tdy_xyzz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 129);

        auto tdz_xyzz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 129);

        auto tdx_xyzz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 130);

        auto tdy_xyzz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 130);

        auto tdz_xyzz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 130);

        auto tdx_xyzz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 131);

        auto tdy_xyzz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 131);

        auto tdz_xyzz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 131);

        auto tdx_xyzz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 132);

        auto tdy_xyzz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 132);

        auto tdz_xyzz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 132);

        auto tdx_xyzz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 133);

        auto tdy_xyzz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 133);

        auto tdz_xyzz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 133);

        auto tdx_xyzz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 134);

        auto tdy_xyzz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 134);

        auto tdz_xyzz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 134);

        auto tdx_xzzz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 135);

        auto tdy_xzzz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 135);

        auto tdz_xzzz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 135);

        auto tdx_xzzz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 136);

        auto tdy_xzzz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 136);

        auto tdz_xzzz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 136);

        auto tdx_xzzz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 137);

        auto tdy_xzzz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 137);

        auto tdz_xzzz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 137);

        auto tdx_xzzz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 138);

        auto tdy_xzzz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 138);

        auto tdz_xzzz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 138);

        auto tdx_xzzz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 139);

        auto tdy_xzzz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 139);

        auto tdz_xzzz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 139);

        auto tdx_xzzz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 140);

        auto tdy_xzzz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 140);

        auto tdz_xzzz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 140);

        auto tdx_xzzz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 141);

        auto tdy_xzzz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 141);

        auto tdz_xzzz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 141);

        auto tdx_xzzz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 142);

        auto tdy_xzzz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 142);

        auto tdz_xzzz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 142);

        auto tdx_xzzz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 143);

        auto tdy_xzzz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 143);

        auto tdz_xzzz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 143);

        auto tdx_xzzz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 144);

        auto tdy_xzzz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 144);

        auto tdz_xzzz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 144);

        // Batch of Integrals (387,435)

        #pragma omp simd aligned(fx, pa_x, tdx_xyzz_xzzz_0, tdx_xyzz_yyyy_0, tdx_xyzz_yyyz_0, \
                                     tdx_xyzz_yyzz_0, tdx_xyzz_yzzz_0, tdx_xyzz_zzzz_0, tdx_xzzz_xxxx_0, tdx_xzzz_xxxy_0, \
                                     tdx_xzzz_xxxz_0, tdx_xzzz_xxyy_0, tdx_xzzz_xxyz_0, tdx_xzzz_xxzz_0, tdx_xzzz_xyyy_0, \
                                     tdx_xzzz_xyyz_0, tdx_xzzz_xyzz_0, tdx_xzzz_xzzz_0, tdx_yzz_xzzz_0, tdx_yzz_yyyy_0, \
                                     tdx_yzz_yyyz_0, tdx_yzz_yyzz_0, tdx_yzz_yzzz_0, tdx_yzz_zzz_0, tdx_yzz_zzzz_0, \
                                     tdx_zzz_xxx_0, tdx_zzz_xxxx_0, tdx_zzz_xxxy_0, tdx_zzz_xxxz_0, tdx_zzz_xxy_0, \
                                     tdx_zzz_xxyy_0, tdx_zzz_xxyz_0, tdx_zzz_xxz_0, tdx_zzz_xxzz_0, tdx_zzz_xyy_0, \
                                     tdx_zzz_xyyy_0, tdx_zzz_xyyz_0, tdx_zzz_xyz_0, tdx_zzz_xyzz_0, tdx_zzz_xzz_0, \
                                     tdx_zzz_xzzz_0, tdx_zzz_yyy_0, tdx_zzz_yyz_0, tdx_zzz_yzz_0, tdx_zzz_zzz_0, \
                                     tdy_xyzz_xzzz_0, tdy_xyzz_yyyy_0, tdy_xyzz_yyyz_0, tdy_xyzz_yyzz_0, tdy_xyzz_yzzz_0, \
                                     tdy_xyzz_zzzz_0, tdy_xzzz_xxxx_0, tdy_xzzz_xxxy_0, tdy_xzzz_xxxz_0, tdy_xzzz_xxyy_0, \
                                     tdy_xzzz_xxyz_0, tdy_xzzz_xxzz_0, tdy_xzzz_xyyy_0, tdy_xzzz_xyyz_0, tdy_xzzz_xyzz_0, \
                                     tdy_xzzz_xzzz_0, tdy_yzz_xzzz_0, tdy_yzz_yyyy_0, tdy_yzz_yyyz_0, tdy_yzz_yyzz_0, \
                                     tdy_yzz_yzzz_0, tdy_yzz_zzz_0, tdy_yzz_zzzz_0, tdy_zzz_xxx_0, tdy_zzz_xxxx_0, \
                                     tdy_zzz_xxxy_0, tdy_zzz_xxxz_0, tdy_zzz_xxy_0, tdy_zzz_xxyy_0, tdy_zzz_xxyz_0, \
                                     tdy_zzz_xxz_0, tdy_zzz_xxzz_0, tdy_zzz_xyy_0, tdy_zzz_xyyy_0, tdy_zzz_xyyz_0, \
                                     tdy_zzz_xyz_0, tdy_zzz_xyzz_0, tdy_zzz_xzz_0, tdy_zzz_xzzz_0, tdy_zzz_yyy_0, \
                                     tdy_zzz_yyz_0, tdy_zzz_yzz_0, tdy_zzz_zzz_0, tdz_xyzz_xzzz_0, tdz_xyzz_yyyy_0, \
                                     tdz_xyzz_yyyz_0, tdz_xyzz_yyzz_0, tdz_xyzz_yzzz_0, tdz_xyzz_zzzz_0, tdz_xzzz_xxxx_0, \
                                     tdz_xzzz_xxxy_0, tdz_xzzz_xxxz_0, tdz_xzzz_xxyy_0, tdz_xzzz_xxyz_0, tdz_xzzz_xxzz_0, \
                                     tdz_xzzz_xyyy_0, tdz_xzzz_xyyz_0, tdz_xzzz_xyzz_0, tdz_xzzz_xzzz_0, tdz_yzz_xzzz_0, \
                                     tdz_yzz_yyyy_0, tdz_yzz_yyyz_0, tdz_yzz_yyzz_0, tdz_yzz_yzzz_0, tdz_yzz_zzz_0, \
                                     tdz_yzz_zzzz_0, tdz_zzz_xxx_0, tdz_zzz_xxxx_0, tdz_zzz_xxxy_0, tdz_zzz_xxxz_0, \
                                     tdz_zzz_xxy_0, tdz_zzz_xxyy_0, tdz_zzz_xxyz_0, tdz_zzz_xxz_0, tdz_zzz_xxzz_0, \
                                     tdz_zzz_xyy_0, tdz_zzz_xyyy_0, tdz_zzz_xyyz_0, tdz_zzz_xyz_0, tdz_zzz_xyzz_0, \
                                     tdz_zzz_xzz_0, tdz_zzz_xzzz_0, tdz_zzz_yyy_0, tdz_zzz_yyz_0, tdz_zzz_yzz_0, \
                                     tdz_zzz_zzz_0, ts_yzz_xzzz_0, ts_yzz_yyyy_0, ts_yzz_yyyz_0, ts_yzz_yyzz_0, \
                                     ts_yzz_yzzz_0, ts_yzz_zzzz_0, ts_zzz_xxxx_0, ts_zzz_xxxy_0, ts_zzz_xxxz_0, \
                                     ts_zzz_xxyy_0, ts_zzz_xxyz_0, ts_zzz_xxzz_0, ts_zzz_xyyy_0, ts_zzz_xyyz_0, \
                                     ts_zzz_xyzz_0, ts_zzz_xzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            tdx_xyzz_xzzz_0[j] = pa_x[j] * tdx_yzz_xzzz_0[j] + 0.5 * fl1_fx * tdx_yzz_zzz_0[j] + 0.5 * fl1_fx * ts_yzz_xzzz_0[j];

            tdy_xyzz_xzzz_0[j] = pa_x[j] * tdy_yzz_xzzz_0[j] + 0.5 * fl1_fx * tdy_yzz_zzz_0[j];

            tdz_xyzz_xzzz_0[j] = pa_x[j] * tdz_yzz_xzzz_0[j] + 0.5 * fl1_fx * tdz_yzz_zzz_0[j];

            tdx_xyzz_yyyy_0[j] = pa_x[j] * tdx_yzz_yyyy_0[j] + 0.5 * fl1_fx * ts_yzz_yyyy_0[j];

            tdy_xyzz_yyyy_0[j] = pa_x[j] * tdy_yzz_yyyy_0[j];

            tdz_xyzz_yyyy_0[j] = pa_x[j] * tdz_yzz_yyyy_0[j];

            tdx_xyzz_yyyz_0[j] = pa_x[j] * tdx_yzz_yyyz_0[j] + 0.5 * fl1_fx * ts_yzz_yyyz_0[j];

            tdy_xyzz_yyyz_0[j] = pa_x[j] * tdy_yzz_yyyz_0[j];

            tdz_xyzz_yyyz_0[j] = pa_x[j] * tdz_yzz_yyyz_0[j];

            tdx_xyzz_yyzz_0[j] = pa_x[j] * tdx_yzz_yyzz_0[j] + 0.5 * fl1_fx * ts_yzz_yyzz_0[j];

            tdy_xyzz_yyzz_0[j] = pa_x[j] * tdy_yzz_yyzz_0[j];

            tdz_xyzz_yyzz_0[j] = pa_x[j] * tdz_yzz_yyzz_0[j];

            tdx_xyzz_yzzz_0[j] = pa_x[j] * tdx_yzz_yzzz_0[j] + 0.5 * fl1_fx * ts_yzz_yzzz_0[j];

            tdy_xyzz_yzzz_0[j] = pa_x[j] * tdy_yzz_yzzz_0[j];

            tdz_xyzz_yzzz_0[j] = pa_x[j] * tdz_yzz_yzzz_0[j];

            tdx_xyzz_zzzz_0[j] = pa_x[j] * tdx_yzz_zzzz_0[j] + 0.5 * fl1_fx * ts_yzz_zzzz_0[j];

            tdy_xyzz_zzzz_0[j] = pa_x[j] * tdy_yzz_zzzz_0[j];

            tdz_xyzz_zzzz_0[j] = pa_x[j] * tdz_yzz_zzzz_0[j];

            tdx_xzzz_xxxx_0[j] = pa_x[j] * tdx_zzz_xxxx_0[j] + 2.0 * fl1_fx * tdx_zzz_xxx_0[j] + 0.5 * fl1_fx * ts_zzz_xxxx_0[j];

            tdy_xzzz_xxxx_0[j] = pa_x[j] * tdy_zzz_xxxx_0[j] + 2.0 * fl1_fx * tdy_zzz_xxx_0[j];

            tdz_xzzz_xxxx_0[j] = pa_x[j] * tdz_zzz_xxxx_0[j] + 2.0 * fl1_fx * tdz_zzz_xxx_0[j];

            tdx_xzzz_xxxy_0[j] = pa_x[j] * tdx_zzz_xxxy_0[j] + 1.5 * fl1_fx * tdx_zzz_xxy_0[j] + 0.5 * fl1_fx * ts_zzz_xxxy_0[j];

            tdy_xzzz_xxxy_0[j] = pa_x[j] * tdy_zzz_xxxy_0[j] + 1.5 * fl1_fx * tdy_zzz_xxy_0[j];

            tdz_xzzz_xxxy_0[j] = pa_x[j] * tdz_zzz_xxxy_0[j] + 1.5 * fl1_fx * tdz_zzz_xxy_0[j];

            tdx_xzzz_xxxz_0[j] = pa_x[j] * tdx_zzz_xxxz_0[j] + 1.5 * fl1_fx * tdx_zzz_xxz_0[j] + 0.5 * fl1_fx * ts_zzz_xxxz_0[j];

            tdy_xzzz_xxxz_0[j] = pa_x[j] * tdy_zzz_xxxz_0[j] + 1.5 * fl1_fx * tdy_zzz_xxz_0[j];

            tdz_xzzz_xxxz_0[j] = pa_x[j] * tdz_zzz_xxxz_0[j] + 1.5 * fl1_fx * tdz_zzz_xxz_0[j];

            tdx_xzzz_xxyy_0[j] = pa_x[j] * tdx_zzz_xxyy_0[j] + fl1_fx * tdx_zzz_xyy_0[j] + 0.5 * fl1_fx * ts_zzz_xxyy_0[j];

            tdy_xzzz_xxyy_0[j] = pa_x[j] * tdy_zzz_xxyy_0[j] + fl1_fx * tdy_zzz_xyy_0[j];

            tdz_xzzz_xxyy_0[j] = pa_x[j] * tdz_zzz_xxyy_0[j] + fl1_fx * tdz_zzz_xyy_0[j];

            tdx_xzzz_xxyz_0[j] = pa_x[j] * tdx_zzz_xxyz_0[j] + fl1_fx * tdx_zzz_xyz_0[j] + 0.5 * fl1_fx * ts_zzz_xxyz_0[j];

            tdy_xzzz_xxyz_0[j] = pa_x[j] * tdy_zzz_xxyz_0[j] + fl1_fx * tdy_zzz_xyz_0[j];

            tdz_xzzz_xxyz_0[j] = pa_x[j] * tdz_zzz_xxyz_0[j] + fl1_fx * tdz_zzz_xyz_0[j];

            tdx_xzzz_xxzz_0[j] = pa_x[j] * tdx_zzz_xxzz_0[j] + fl1_fx * tdx_zzz_xzz_0[j] + 0.5 * fl1_fx * ts_zzz_xxzz_0[j];

            tdy_xzzz_xxzz_0[j] = pa_x[j] * tdy_zzz_xxzz_0[j] + fl1_fx * tdy_zzz_xzz_0[j];

            tdz_xzzz_xxzz_0[j] = pa_x[j] * tdz_zzz_xxzz_0[j] + fl1_fx * tdz_zzz_xzz_0[j];

            tdx_xzzz_xyyy_0[j] = pa_x[j] * tdx_zzz_xyyy_0[j] + 0.5 * fl1_fx * tdx_zzz_yyy_0[j] + 0.5 * fl1_fx * ts_zzz_xyyy_0[j];

            tdy_xzzz_xyyy_0[j] = pa_x[j] * tdy_zzz_xyyy_0[j] + 0.5 * fl1_fx * tdy_zzz_yyy_0[j];

            tdz_xzzz_xyyy_0[j] = pa_x[j] * tdz_zzz_xyyy_0[j] + 0.5 * fl1_fx * tdz_zzz_yyy_0[j];

            tdx_xzzz_xyyz_0[j] = pa_x[j] * tdx_zzz_xyyz_0[j] + 0.5 * fl1_fx * tdx_zzz_yyz_0[j] + 0.5 * fl1_fx * ts_zzz_xyyz_0[j];

            tdy_xzzz_xyyz_0[j] = pa_x[j] * tdy_zzz_xyyz_0[j] + 0.5 * fl1_fx * tdy_zzz_yyz_0[j];

            tdz_xzzz_xyyz_0[j] = pa_x[j] * tdz_zzz_xyyz_0[j] + 0.5 * fl1_fx * tdz_zzz_yyz_0[j];

            tdx_xzzz_xyzz_0[j] = pa_x[j] * tdx_zzz_xyzz_0[j] + 0.5 * fl1_fx * tdx_zzz_yzz_0[j] + 0.5 * fl1_fx * ts_zzz_xyzz_0[j];

            tdy_xzzz_xyzz_0[j] = pa_x[j] * tdy_zzz_xyzz_0[j] + 0.5 * fl1_fx * tdy_zzz_yzz_0[j];

            tdz_xzzz_xyzz_0[j] = pa_x[j] * tdz_zzz_xyzz_0[j] + 0.5 * fl1_fx * tdz_zzz_yzz_0[j];

            tdx_xzzz_xzzz_0[j] = pa_x[j] * tdx_zzz_xzzz_0[j] + 0.5 * fl1_fx * tdx_zzz_zzz_0[j] + 0.5 * fl1_fx * ts_zzz_xzzz_0[j];

            tdy_xzzz_xzzz_0[j] = pa_x[j] * tdy_zzz_xzzz_0[j] + 0.5 * fl1_fx * tdy_zzz_zzz_0[j];

            tdz_xzzz_xzzz_0[j] = pa_x[j] * tdz_zzz_xzzz_0[j] + 0.5 * fl1_fx * tdz_zzz_zzz_0[j];
        }

        idx++;
    }
}

void
compElectricDipoleForGG_435_483(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const int32_t              nOSFactors,
                                const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_d_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_d_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tdx_yyy_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 100);

        auto tdy_yyy_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 100);

        auto tdz_yyy_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 100);

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

        auto tdx_yy_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 55);

        auto tdy_yy_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 55);

        auto tdz_yy_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 55);

        auto tdx_yyy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 60);

        auto tdy_yyy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 60);

        auto tdz_yyy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 60);

        auto tdx_yyy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 61);

        auto tdy_yyy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 61);

        auto tdz_yyy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 61);

        auto tdx_yyy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 62);

        auto tdy_yyy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 62);

        auto tdz_yyy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 62);

        auto tdx_yyy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 63);

        auto tdy_yyy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 63);

        auto tdz_yyy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 63);

        auto tdx_yyy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 64);

        auto tdy_yyy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 64);

        auto tdz_yyy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 64);

        auto tdx_yyy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 65);

        auto tdy_yyy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 65);

        auto tdz_yyy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 65);

        auto tdx_yyy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 66);

        auto tdy_yyy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 66);

        auto tdz_yyy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 66);

        auto ts_yyy_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 90);

        auto ts_yyy_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 91);

        auto ts_yyy_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 92);

        auto ts_yyy_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 93);

        auto ts_yyy_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 94);

        auto ts_yyy_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 95);

        auto ts_yyy_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 96);

        auto ts_yyy_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 97);

        auto ts_yyy_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 98);

        auto ts_yyy_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 99);

        auto ts_yyy_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 100);

        auto ts_zzz_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 145);

        auto ts_zzz_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 146);

        auto ts_zzz_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 147);

        auto ts_zzz_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 148);

        auto ts_zzz_zzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 149);

        // set up pointers to integrals

        auto tdx_xzzz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 145);

        auto tdy_xzzz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 145);

        auto tdz_xzzz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 145);

        auto tdx_xzzz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 146);

        auto tdy_xzzz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 146);

        auto tdz_xzzz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 146);

        auto tdx_xzzz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 147);

        auto tdy_xzzz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 147);

        auto tdz_xzzz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 147);

        auto tdx_xzzz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 148);

        auto tdy_xzzz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 148);

        auto tdz_xzzz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 148);

        auto tdx_xzzz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 149);

        auto tdy_xzzz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 149);

        auto tdz_xzzz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 149);

        auto tdx_yyyy_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 150);

        auto tdy_yyyy_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 150);

        auto tdz_yyyy_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 150);

        auto tdx_yyyy_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 151);

        auto tdy_yyyy_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 151);

        auto tdz_yyyy_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 151);

        auto tdx_yyyy_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 152);

        auto tdy_yyyy_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 152);

        auto tdz_yyyy_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 152);

        auto tdx_yyyy_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 153);

        auto tdy_yyyy_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 153);

        auto tdz_yyyy_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 153);

        auto tdx_yyyy_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 154);

        auto tdy_yyyy_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 154);

        auto tdz_yyyy_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 154);

        auto tdx_yyyy_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 155);

        auto tdy_yyyy_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 155);

        auto tdz_yyyy_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 155);

        auto tdx_yyyy_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 156);

        auto tdy_yyyy_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 156);

        auto tdz_yyyy_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 156);

        auto tdx_yyyy_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 157);

        auto tdy_yyyy_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 157);

        auto tdz_yyyy_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 157);

        auto tdx_yyyy_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 158);

        auto tdy_yyyy_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 158);

        auto tdz_yyyy_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 158);

        auto tdx_yyyy_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 159);

        auto tdy_yyyy_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 159);

        auto tdz_yyyy_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 159);

        auto tdx_yyyy_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 160);

        auto tdy_yyyy_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 160);

        auto tdz_yyyy_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 160);

        // Batch of Integrals (435,483)

        #pragma omp simd aligned(fx, pa_x, pa_y, tdx_xzzz_yyyy_0, tdx_xzzz_yyyz_0, tdx_xzzz_yyzz_0, \
                                     tdx_xzzz_yzzz_0, tdx_xzzz_zzzz_0, tdx_yy_xxxx_0, tdx_yy_xxxy_0, tdx_yy_xxxz_0, \
                                     tdx_yy_xxyy_0, tdx_yy_xxyz_0, tdx_yy_xxzz_0, tdx_yy_xyyy_0, tdx_yy_xyyz_0, \
                                     tdx_yy_xyzz_0, tdx_yy_xzzz_0, tdx_yy_yyyy_0, tdx_yyy_xxx_0, tdx_yyy_xxxx_0, \
                                     tdx_yyy_xxxy_0, tdx_yyy_xxxz_0, tdx_yyy_xxy_0, tdx_yyy_xxyy_0, tdx_yyy_xxyz_0, \
                                     tdx_yyy_xxz_0, tdx_yyy_xxzz_0, tdx_yyy_xyy_0, tdx_yyy_xyyy_0, tdx_yyy_xyyz_0, \
                                     tdx_yyy_xyz_0, tdx_yyy_xyzz_0, tdx_yyy_xzz_0, tdx_yyy_xzzz_0, tdx_yyy_yyy_0, \
                                     tdx_yyy_yyyy_0, tdx_yyyy_xxxx_0, tdx_yyyy_xxxy_0, tdx_yyyy_xxxz_0, tdx_yyyy_xxyy_0, \
                                     tdx_yyyy_xxyz_0, tdx_yyyy_xxzz_0, tdx_yyyy_xyyy_0, tdx_yyyy_xyyz_0, tdx_yyyy_xyzz_0, \
                                     tdx_yyyy_xzzz_0, tdx_yyyy_yyyy_0, tdx_zzz_yyyy_0, tdx_zzz_yyyz_0, tdx_zzz_yyzz_0, \
                                     tdx_zzz_yzzz_0, tdx_zzz_zzzz_0, tdy_xzzz_yyyy_0, tdy_xzzz_yyyz_0, tdy_xzzz_yyzz_0, \
                                     tdy_xzzz_yzzz_0, tdy_xzzz_zzzz_0, tdy_yy_xxxx_0, tdy_yy_xxxy_0, tdy_yy_xxxz_0, \
                                     tdy_yy_xxyy_0, tdy_yy_xxyz_0, tdy_yy_xxzz_0, tdy_yy_xyyy_0, tdy_yy_xyyz_0, \
                                     tdy_yy_xyzz_0, tdy_yy_xzzz_0, tdy_yy_yyyy_0, tdy_yyy_xxx_0, tdy_yyy_xxxx_0, \
                                     tdy_yyy_xxxy_0, tdy_yyy_xxxz_0, tdy_yyy_xxy_0, tdy_yyy_xxyy_0, tdy_yyy_xxyz_0, \
                                     tdy_yyy_xxz_0, tdy_yyy_xxzz_0, tdy_yyy_xyy_0, tdy_yyy_xyyy_0, tdy_yyy_xyyz_0, \
                                     tdy_yyy_xyz_0, tdy_yyy_xyzz_0, tdy_yyy_xzz_0, tdy_yyy_xzzz_0, tdy_yyy_yyy_0, \
                                     tdy_yyy_yyyy_0, tdy_yyyy_xxxx_0, tdy_yyyy_xxxy_0, tdy_yyyy_xxxz_0, tdy_yyyy_xxyy_0, \
                                     tdy_yyyy_xxyz_0, tdy_yyyy_xxzz_0, tdy_yyyy_xyyy_0, tdy_yyyy_xyyz_0, tdy_yyyy_xyzz_0, \
                                     tdy_yyyy_xzzz_0, tdy_yyyy_yyyy_0, tdy_zzz_yyyy_0, tdy_zzz_yyyz_0, tdy_zzz_yyzz_0, \
                                     tdy_zzz_yzzz_0, tdy_zzz_zzzz_0, tdz_xzzz_yyyy_0, tdz_xzzz_yyyz_0, tdz_xzzz_yyzz_0, \
                                     tdz_xzzz_yzzz_0, tdz_xzzz_zzzz_0, tdz_yy_xxxx_0, tdz_yy_xxxy_0, tdz_yy_xxxz_0, \
                                     tdz_yy_xxyy_0, tdz_yy_xxyz_0, tdz_yy_xxzz_0, tdz_yy_xyyy_0, tdz_yy_xyyz_0, \
                                     tdz_yy_xyzz_0, tdz_yy_xzzz_0, tdz_yy_yyyy_0, tdz_yyy_xxx_0, tdz_yyy_xxxx_0, \
                                     tdz_yyy_xxxy_0, tdz_yyy_xxxz_0, tdz_yyy_xxy_0, tdz_yyy_xxyy_0, tdz_yyy_xxyz_0, \
                                     tdz_yyy_xxz_0, tdz_yyy_xxzz_0, tdz_yyy_xyy_0, tdz_yyy_xyyy_0, tdz_yyy_xyyz_0, \
                                     tdz_yyy_xyz_0, tdz_yyy_xyzz_0, tdz_yyy_xzz_0, tdz_yyy_xzzz_0, tdz_yyy_yyy_0, \
                                     tdz_yyy_yyyy_0, tdz_yyyy_xxxx_0, tdz_yyyy_xxxy_0, tdz_yyyy_xxxz_0, tdz_yyyy_xxyy_0, \
                                     tdz_yyyy_xxyz_0, tdz_yyyy_xxzz_0, tdz_yyyy_xyyy_0, tdz_yyyy_xyyz_0, tdz_yyyy_xyzz_0, \
                                     tdz_yyyy_xzzz_0, tdz_yyyy_yyyy_0, tdz_zzz_yyyy_0, tdz_zzz_yyyz_0, tdz_zzz_yyzz_0, \
                                     tdz_zzz_yzzz_0, tdz_zzz_zzzz_0, ts_yyy_xxxx_0, ts_yyy_xxxy_0, ts_yyy_xxxz_0, \
                                     ts_yyy_xxyy_0, ts_yyy_xxyz_0, ts_yyy_xxzz_0, ts_yyy_xyyy_0, ts_yyy_xyyz_0, \
                                     ts_yyy_xyzz_0, ts_yyy_xzzz_0, ts_yyy_yyyy_0, ts_zzz_yyyy_0, ts_zzz_yyyz_0, \
                                     ts_zzz_yyzz_0, ts_zzz_yzzz_0, ts_zzz_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            tdx_xzzz_yyyy_0[j] = pa_x[j] * tdx_zzz_yyyy_0[j] + 0.5 * fl1_fx * ts_zzz_yyyy_0[j];

            tdy_xzzz_yyyy_0[j] = pa_x[j] * tdy_zzz_yyyy_0[j];

            tdz_xzzz_yyyy_0[j] = pa_x[j] * tdz_zzz_yyyy_0[j];

            tdx_xzzz_yyyz_0[j] = pa_x[j] * tdx_zzz_yyyz_0[j] + 0.5 * fl1_fx * ts_zzz_yyyz_0[j];

            tdy_xzzz_yyyz_0[j] = pa_x[j] * tdy_zzz_yyyz_0[j];

            tdz_xzzz_yyyz_0[j] = pa_x[j] * tdz_zzz_yyyz_0[j];

            tdx_xzzz_yyzz_0[j] = pa_x[j] * tdx_zzz_yyzz_0[j] + 0.5 * fl1_fx * ts_zzz_yyzz_0[j];

            tdy_xzzz_yyzz_0[j] = pa_x[j] * tdy_zzz_yyzz_0[j];

            tdz_xzzz_yyzz_0[j] = pa_x[j] * tdz_zzz_yyzz_0[j];

            tdx_xzzz_yzzz_0[j] = pa_x[j] * tdx_zzz_yzzz_0[j] + 0.5 * fl1_fx * ts_zzz_yzzz_0[j];

            tdy_xzzz_yzzz_0[j] = pa_x[j] * tdy_zzz_yzzz_0[j];

            tdz_xzzz_yzzz_0[j] = pa_x[j] * tdz_zzz_yzzz_0[j];

            tdx_xzzz_zzzz_0[j] = pa_x[j] * tdx_zzz_zzzz_0[j] + 0.5 * fl1_fx * ts_zzz_zzzz_0[j];

            tdy_xzzz_zzzz_0[j] = pa_x[j] * tdy_zzz_zzzz_0[j];

            tdz_xzzz_zzzz_0[j] = pa_x[j] * tdz_zzz_zzzz_0[j];

            tdx_yyyy_xxxx_0[j] = pa_y[j] * tdx_yyy_xxxx_0[j] + 1.5 * fl1_fx * tdx_yy_xxxx_0[j];

            tdy_yyyy_xxxx_0[j] = pa_y[j] * tdy_yyy_xxxx_0[j] + 1.5 * fl1_fx * tdy_yy_xxxx_0[j] + 0.5 * fl1_fx * ts_yyy_xxxx_0[j];

            tdz_yyyy_xxxx_0[j] = pa_y[j] * tdz_yyy_xxxx_0[j] + 1.5 * fl1_fx * tdz_yy_xxxx_0[j];

            tdx_yyyy_xxxy_0[j] = pa_y[j] * tdx_yyy_xxxy_0[j] + 1.5 * fl1_fx * tdx_yy_xxxy_0[j] + 0.5 * fl1_fx * tdx_yyy_xxx_0[j];

            tdy_yyyy_xxxy_0[j] =
                pa_y[j] * tdy_yyy_xxxy_0[j] + 1.5 * fl1_fx * tdy_yy_xxxy_0[j] + 0.5 * fl1_fx * tdy_yyy_xxx_0[j] + 0.5 * fl1_fx * ts_yyy_xxxy_0[j];

            tdz_yyyy_xxxy_0[j] = pa_y[j] * tdz_yyy_xxxy_0[j] + 1.5 * fl1_fx * tdz_yy_xxxy_0[j] + 0.5 * fl1_fx * tdz_yyy_xxx_0[j];

            tdx_yyyy_xxxz_0[j] = pa_y[j] * tdx_yyy_xxxz_0[j] + 1.5 * fl1_fx * tdx_yy_xxxz_0[j];

            tdy_yyyy_xxxz_0[j] = pa_y[j] * tdy_yyy_xxxz_0[j] + 1.5 * fl1_fx * tdy_yy_xxxz_0[j] + 0.5 * fl1_fx * ts_yyy_xxxz_0[j];

            tdz_yyyy_xxxz_0[j] = pa_y[j] * tdz_yyy_xxxz_0[j] + 1.5 * fl1_fx * tdz_yy_xxxz_0[j];

            tdx_yyyy_xxyy_0[j] = pa_y[j] * tdx_yyy_xxyy_0[j] + 1.5 * fl1_fx * tdx_yy_xxyy_0[j] + fl1_fx * tdx_yyy_xxy_0[j];

            tdy_yyyy_xxyy_0[j] =
                pa_y[j] * tdy_yyy_xxyy_0[j] + 1.5 * fl1_fx * tdy_yy_xxyy_0[j] + fl1_fx * tdy_yyy_xxy_0[j] + 0.5 * fl1_fx * ts_yyy_xxyy_0[j];

            tdz_yyyy_xxyy_0[j] = pa_y[j] * tdz_yyy_xxyy_0[j] + 1.5 * fl1_fx * tdz_yy_xxyy_0[j] + fl1_fx * tdz_yyy_xxy_0[j];

            tdx_yyyy_xxyz_0[j] = pa_y[j] * tdx_yyy_xxyz_0[j] + 1.5 * fl1_fx * tdx_yy_xxyz_0[j] + 0.5 * fl1_fx * tdx_yyy_xxz_0[j];

            tdy_yyyy_xxyz_0[j] =
                pa_y[j] * tdy_yyy_xxyz_0[j] + 1.5 * fl1_fx * tdy_yy_xxyz_0[j] + 0.5 * fl1_fx * tdy_yyy_xxz_0[j] + 0.5 * fl1_fx * ts_yyy_xxyz_0[j];

            tdz_yyyy_xxyz_0[j] = pa_y[j] * tdz_yyy_xxyz_0[j] + 1.5 * fl1_fx * tdz_yy_xxyz_0[j] + 0.5 * fl1_fx * tdz_yyy_xxz_0[j];

            tdx_yyyy_xxzz_0[j] = pa_y[j] * tdx_yyy_xxzz_0[j] + 1.5 * fl1_fx * tdx_yy_xxzz_0[j];

            tdy_yyyy_xxzz_0[j] = pa_y[j] * tdy_yyy_xxzz_0[j] + 1.5 * fl1_fx * tdy_yy_xxzz_0[j] + 0.5 * fl1_fx * ts_yyy_xxzz_0[j];

            tdz_yyyy_xxzz_0[j] = pa_y[j] * tdz_yyy_xxzz_0[j] + 1.5 * fl1_fx * tdz_yy_xxzz_0[j];

            tdx_yyyy_xyyy_0[j] = pa_y[j] * tdx_yyy_xyyy_0[j] + 1.5 * fl1_fx * tdx_yy_xyyy_0[j] + 1.5 * fl1_fx * tdx_yyy_xyy_0[j];

            tdy_yyyy_xyyy_0[j] =
                pa_y[j] * tdy_yyy_xyyy_0[j] + 1.5 * fl1_fx * tdy_yy_xyyy_0[j] + 1.5 * fl1_fx * tdy_yyy_xyy_0[j] + 0.5 * fl1_fx * ts_yyy_xyyy_0[j];

            tdz_yyyy_xyyy_0[j] = pa_y[j] * tdz_yyy_xyyy_0[j] + 1.5 * fl1_fx * tdz_yy_xyyy_0[j] + 1.5 * fl1_fx * tdz_yyy_xyy_0[j];

            tdx_yyyy_xyyz_0[j] = pa_y[j] * tdx_yyy_xyyz_0[j] + 1.5 * fl1_fx * tdx_yy_xyyz_0[j] + fl1_fx * tdx_yyy_xyz_0[j];

            tdy_yyyy_xyyz_0[j] =
                pa_y[j] * tdy_yyy_xyyz_0[j] + 1.5 * fl1_fx * tdy_yy_xyyz_0[j] + fl1_fx * tdy_yyy_xyz_0[j] + 0.5 * fl1_fx * ts_yyy_xyyz_0[j];

            tdz_yyyy_xyyz_0[j] = pa_y[j] * tdz_yyy_xyyz_0[j] + 1.5 * fl1_fx * tdz_yy_xyyz_0[j] + fl1_fx * tdz_yyy_xyz_0[j];

            tdx_yyyy_xyzz_0[j] = pa_y[j] * tdx_yyy_xyzz_0[j] + 1.5 * fl1_fx * tdx_yy_xyzz_0[j] + 0.5 * fl1_fx * tdx_yyy_xzz_0[j];

            tdy_yyyy_xyzz_0[j] =
                pa_y[j] * tdy_yyy_xyzz_0[j] + 1.5 * fl1_fx * tdy_yy_xyzz_0[j] + 0.5 * fl1_fx * tdy_yyy_xzz_0[j] + 0.5 * fl1_fx * ts_yyy_xyzz_0[j];

            tdz_yyyy_xyzz_0[j] = pa_y[j] * tdz_yyy_xyzz_0[j] + 1.5 * fl1_fx * tdz_yy_xyzz_0[j] + 0.5 * fl1_fx * tdz_yyy_xzz_0[j];

            tdx_yyyy_xzzz_0[j] = pa_y[j] * tdx_yyy_xzzz_0[j] + 1.5 * fl1_fx * tdx_yy_xzzz_0[j];

            tdy_yyyy_xzzz_0[j] = pa_y[j] * tdy_yyy_xzzz_0[j] + 1.5 * fl1_fx * tdy_yy_xzzz_0[j] + 0.5 * fl1_fx * ts_yyy_xzzz_0[j];

            tdz_yyyy_xzzz_0[j] = pa_y[j] * tdz_yyy_xzzz_0[j] + 1.5 * fl1_fx * tdz_yy_xzzz_0[j];

            tdx_yyyy_yyyy_0[j] = pa_y[j] * tdx_yyy_yyyy_0[j] + 1.5 * fl1_fx * tdx_yy_yyyy_0[j] + 2.0 * fl1_fx * tdx_yyy_yyy_0[j];

            tdy_yyyy_yyyy_0[j] =
                pa_y[j] * tdy_yyy_yyyy_0[j] + 1.5 * fl1_fx * tdy_yy_yyyy_0[j] + 2.0 * fl1_fx * tdy_yyy_yyy_0[j] + 0.5 * fl1_fx * ts_yyy_yyyy_0[j];

            tdz_yyyy_yyyy_0[j] = pa_y[j] * tdz_yyy_yyyy_0[j] + 1.5 * fl1_fx * tdz_yy_yyyy_0[j] + 2.0 * fl1_fx * tdz_yyy_yyy_0[j];
        }

        idx++;
    }
}

void
compElectricDipoleForGG_483_531(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const int32_t              nOSFactors,
                                const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_d_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_d_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_y = paDistances.data(3 * idx + 1);

        // set up pointers to auxilary integrals

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

        auto tdz_yyz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 116);

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

        auto tdz_yz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 71);

        auto tdx_yyy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 67);

        auto tdy_yyy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 67);

        auto tdz_yyy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 67);

        auto tdx_yyy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 68);

        auto tdy_yyy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 68);

        auto tdz_yyy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 68);

        auto tdx_yyy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 69);

        auto tdy_yyy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 69);

        auto tdz_yyy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 69);

        auto tdx_yyz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 70);

        auto tdy_yyz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 70);

        auto tdz_yyz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 70);

        auto tdx_yyz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 71);

        auto tdy_yyz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 71);

        auto tdz_yyz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 71);

        auto tdx_yyz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 72);

        auto tdy_yyz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 72);

        auto tdz_yyz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 72);

        auto tdx_yyz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 73);

        auto tdy_yyz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 73);

        auto tdz_yyz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 73);

        auto tdx_yyz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 74);

        auto tdy_yyz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 74);

        auto tdz_yyz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 74);

        auto tdx_yyz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 75);

        auto tdy_yyz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 75);

        auto tdz_yyz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 75);

        auto tdx_yyz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 76);

        auto tdy_yyz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 76);

        auto tdz_yyz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 76);

        auto tdx_yyz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 77);

        auto tdy_yyz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 77);

        auto tdz_yyz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 77);

        auto ts_yyy_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 101);

        auto ts_yyy_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 102);

        auto ts_yyy_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 103);

        auto ts_yyy_zzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 104);

        auto ts_yyz_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 105);

        auto ts_yyz_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 106);

        auto ts_yyz_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 107);

        auto ts_yyz_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 108);

        auto ts_yyz_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 109);

        auto ts_yyz_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 110);

        auto ts_yyz_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 111);

        auto ts_yyz_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 112);

        auto ts_yyz_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 113);

        auto ts_yyz_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 114);

        auto ts_yyz_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 115);

        auto ts_yyz_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 116);

        // set up pointers to integrals

        auto tdx_yyyy_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 161);

        auto tdy_yyyy_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 161);

        auto tdz_yyyy_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 161);

        auto tdx_yyyy_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 162);

        auto tdy_yyyy_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 162);

        auto tdz_yyyy_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 162);

        auto tdx_yyyy_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 163);

        auto tdy_yyyy_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 163);

        auto tdz_yyyy_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 163);

        auto tdx_yyyy_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 164);

        auto tdy_yyyy_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 164);

        auto tdz_yyyy_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 164);

        auto tdx_yyyz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 165);

        auto tdy_yyyz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 165);

        auto tdz_yyyz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 165);

        auto tdx_yyyz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 166);

        auto tdy_yyyz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 166);

        auto tdz_yyyz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 166);

        auto tdx_yyyz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 167);

        auto tdy_yyyz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 167);

        auto tdz_yyyz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 167);

        auto tdx_yyyz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 168);

        auto tdy_yyyz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 168);

        auto tdz_yyyz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 168);

        auto tdx_yyyz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 169);

        auto tdy_yyyz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 169);

        auto tdz_yyyz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 169);

        auto tdx_yyyz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 170);

        auto tdy_yyyz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 170);

        auto tdz_yyyz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 170);

        auto tdx_yyyz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 171);

        auto tdy_yyyz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 171);

        auto tdz_yyyz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 171);

        auto tdx_yyyz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 172);

        auto tdy_yyyz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 172);

        auto tdz_yyyz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 172);

        auto tdx_yyyz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 173);

        auto tdy_yyyz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 173);

        auto tdz_yyyz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 173);

        auto tdx_yyyz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 174);

        auto tdy_yyyz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 174);

        auto tdz_yyyz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 174);

        auto tdx_yyyz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 175);

        auto tdy_yyyz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 175);

        auto tdz_yyyz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 175);

        auto tdx_yyyz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 176);

        auto tdy_yyyz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 176);

        auto tdz_yyyz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 176);

        // Batch of Integrals (483,531)

        #pragma omp simd aligned(fx, pa_y, tdx_yy_yyyz_0, tdx_yy_yyzz_0, tdx_yy_yzzz_0, tdx_yy_zzzz_0, \
                                     tdx_yyy_yyyz_0, tdx_yyy_yyz_0, tdx_yyy_yyzz_0, tdx_yyy_yzz_0, tdx_yyy_yzzz_0, \
                                     tdx_yyy_zzz_0, tdx_yyy_zzzz_0, tdx_yyyy_yyyz_0, tdx_yyyy_yyzz_0, tdx_yyyy_yzzz_0, \
                                     tdx_yyyy_zzzz_0, tdx_yyyz_xxxx_0, tdx_yyyz_xxxy_0, tdx_yyyz_xxxz_0, tdx_yyyz_xxyy_0, \
                                     tdx_yyyz_xxyz_0, tdx_yyyz_xxzz_0, tdx_yyyz_xyyy_0, tdx_yyyz_xyyz_0, tdx_yyyz_xyzz_0, \
                                     tdx_yyyz_xzzz_0, tdx_yyyz_yyyy_0, tdx_yyyz_yyyz_0, tdx_yyz_xxx_0, tdx_yyz_xxxx_0, \
                                     tdx_yyz_xxxy_0, tdx_yyz_xxxz_0, tdx_yyz_xxy_0, tdx_yyz_xxyy_0, tdx_yyz_xxyz_0, \
                                     tdx_yyz_xxz_0, tdx_yyz_xxzz_0, tdx_yyz_xyy_0, tdx_yyz_xyyy_0, tdx_yyz_xyyz_0, \
                                     tdx_yyz_xyz_0, tdx_yyz_xyzz_0, tdx_yyz_xzz_0, tdx_yyz_xzzz_0, tdx_yyz_yyy_0, \
                                     tdx_yyz_yyyy_0, tdx_yyz_yyyz_0, tdx_yyz_yyz_0, tdx_yz_xxxx_0, tdx_yz_xxxy_0, \
                                     tdx_yz_xxxz_0, tdx_yz_xxyy_0, tdx_yz_xxyz_0, tdx_yz_xxzz_0, tdx_yz_xyyy_0, \
                                     tdx_yz_xyyz_0, tdx_yz_xyzz_0, tdx_yz_xzzz_0, tdx_yz_yyyy_0, tdx_yz_yyyz_0, \
                                     tdy_yy_yyyz_0, tdy_yy_yyzz_0, tdy_yy_yzzz_0, tdy_yy_zzzz_0, tdy_yyy_yyyz_0, \
                                     tdy_yyy_yyz_0, tdy_yyy_yyzz_0, tdy_yyy_yzz_0, tdy_yyy_yzzz_0, tdy_yyy_zzz_0, \
                                     tdy_yyy_zzzz_0, tdy_yyyy_yyyz_0, tdy_yyyy_yyzz_0, tdy_yyyy_yzzz_0, tdy_yyyy_zzzz_0, \
                                     tdy_yyyz_xxxx_0, tdy_yyyz_xxxy_0, tdy_yyyz_xxxz_0, tdy_yyyz_xxyy_0, tdy_yyyz_xxyz_0, \
                                     tdy_yyyz_xxzz_0, tdy_yyyz_xyyy_0, tdy_yyyz_xyyz_0, tdy_yyyz_xyzz_0, tdy_yyyz_xzzz_0, \
                                     tdy_yyyz_yyyy_0, tdy_yyyz_yyyz_0, tdy_yyz_xxx_0, tdy_yyz_xxxx_0, tdy_yyz_xxxy_0, \
                                     tdy_yyz_xxxz_0, tdy_yyz_xxy_0, tdy_yyz_xxyy_0, tdy_yyz_xxyz_0, tdy_yyz_xxz_0, \
                                     tdy_yyz_xxzz_0, tdy_yyz_xyy_0, tdy_yyz_xyyy_0, tdy_yyz_xyyz_0, tdy_yyz_xyz_0, \
                                     tdy_yyz_xyzz_0, tdy_yyz_xzz_0, tdy_yyz_xzzz_0, tdy_yyz_yyy_0, tdy_yyz_yyyy_0, \
                                     tdy_yyz_yyyz_0, tdy_yyz_yyz_0, tdy_yz_xxxx_0, tdy_yz_xxxy_0, tdy_yz_xxxz_0, \
                                     tdy_yz_xxyy_0, tdy_yz_xxyz_0, tdy_yz_xxzz_0, tdy_yz_xyyy_0, tdy_yz_xyyz_0, \
                                     tdy_yz_xyzz_0, tdy_yz_xzzz_0, tdy_yz_yyyy_0, tdy_yz_yyyz_0, tdz_yy_yyyz_0, \
                                     tdz_yy_yyzz_0, tdz_yy_yzzz_0, tdz_yy_zzzz_0, tdz_yyy_yyyz_0, tdz_yyy_yyz_0, \
                                     tdz_yyy_yyzz_0, tdz_yyy_yzz_0, tdz_yyy_yzzz_0, tdz_yyy_zzz_0, tdz_yyy_zzzz_0, \
                                     tdz_yyyy_yyyz_0, tdz_yyyy_yyzz_0, tdz_yyyy_yzzz_0, tdz_yyyy_zzzz_0, tdz_yyyz_xxxx_0, \
                                     tdz_yyyz_xxxy_0, tdz_yyyz_xxxz_0, tdz_yyyz_xxyy_0, tdz_yyyz_xxyz_0, tdz_yyyz_xxzz_0, \
                                     tdz_yyyz_xyyy_0, tdz_yyyz_xyyz_0, tdz_yyyz_xyzz_0, tdz_yyyz_xzzz_0, tdz_yyyz_yyyy_0, \
                                     tdz_yyyz_yyyz_0, tdz_yyz_xxx_0, tdz_yyz_xxxx_0, tdz_yyz_xxxy_0, tdz_yyz_xxxz_0, \
                                     tdz_yyz_xxy_0, tdz_yyz_xxyy_0, tdz_yyz_xxyz_0, tdz_yyz_xxz_0, tdz_yyz_xxzz_0, \
                                     tdz_yyz_xyy_0, tdz_yyz_xyyy_0, tdz_yyz_xyyz_0, tdz_yyz_xyz_0, tdz_yyz_xyzz_0, \
                                     tdz_yyz_xzz_0, tdz_yyz_xzzz_0, tdz_yyz_yyy_0, tdz_yyz_yyyy_0, tdz_yyz_yyyz_0, \
                                     tdz_yyz_yyz_0, tdz_yz_xxxx_0, tdz_yz_xxxy_0, tdz_yz_xxxz_0, tdz_yz_xxyy_0, \
                                     tdz_yz_xxyz_0, tdz_yz_xxzz_0, tdz_yz_xyyy_0, tdz_yz_xyyz_0, tdz_yz_xyzz_0, \
                                     tdz_yz_xzzz_0, tdz_yz_yyyy_0, tdz_yz_yyyz_0, ts_yyy_yyyz_0, ts_yyy_yyzz_0, \
                                     ts_yyy_yzzz_0, ts_yyy_zzzz_0, ts_yyz_xxxx_0, ts_yyz_xxxy_0, ts_yyz_xxxz_0, \
                                     ts_yyz_xxyy_0, ts_yyz_xxyz_0, ts_yyz_xxzz_0, ts_yyz_xyyy_0, ts_yyz_xyyz_0, \
                                     ts_yyz_xyzz_0, ts_yyz_xzzz_0, ts_yyz_yyyy_0, ts_yyz_yyyz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            tdx_yyyy_yyyz_0[j] = pa_y[j] * tdx_yyy_yyyz_0[j] + 1.5 * fl1_fx * tdx_yy_yyyz_0[j] + 1.5 * fl1_fx * tdx_yyy_yyz_0[j];

            tdy_yyyy_yyyz_0[j] =
                pa_y[j] * tdy_yyy_yyyz_0[j] + 1.5 * fl1_fx * tdy_yy_yyyz_0[j] + 1.5 * fl1_fx * tdy_yyy_yyz_0[j] + 0.5 * fl1_fx * ts_yyy_yyyz_0[j];

            tdz_yyyy_yyyz_0[j] = pa_y[j] * tdz_yyy_yyyz_0[j] + 1.5 * fl1_fx * tdz_yy_yyyz_0[j] + 1.5 * fl1_fx * tdz_yyy_yyz_0[j];

            tdx_yyyy_yyzz_0[j] = pa_y[j] * tdx_yyy_yyzz_0[j] + 1.5 * fl1_fx * tdx_yy_yyzz_0[j] + fl1_fx * tdx_yyy_yzz_0[j];

            tdy_yyyy_yyzz_0[j] =
                pa_y[j] * tdy_yyy_yyzz_0[j] + 1.5 * fl1_fx * tdy_yy_yyzz_0[j] + fl1_fx * tdy_yyy_yzz_0[j] + 0.5 * fl1_fx * ts_yyy_yyzz_0[j];

            tdz_yyyy_yyzz_0[j] = pa_y[j] * tdz_yyy_yyzz_0[j] + 1.5 * fl1_fx * tdz_yy_yyzz_0[j] + fl1_fx * tdz_yyy_yzz_0[j];

            tdx_yyyy_yzzz_0[j] = pa_y[j] * tdx_yyy_yzzz_0[j] + 1.5 * fl1_fx * tdx_yy_yzzz_0[j] + 0.5 * fl1_fx * tdx_yyy_zzz_0[j];

            tdy_yyyy_yzzz_0[j] =
                pa_y[j] * tdy_yyy_yzzz_0[j] + 1.5 * fl1_fx * tdy_yy_yzzz_0[j] + 0.5 * fl1_fx * tdy_yyy_zzz_0[j] + 0.5 * fl1_fx * ts_yyy_yzzz_0[j];

            tdz_yyyy_yzzz_0[j] = pa_y[j] * tdz_yyy_yzzz_0[j] + 1.5 * fl1_fx * tdz_yy_yzzz_0[j] + 0.5 * fl1_fx * tdz_yyy_zzz_0[j];

            tdx_yyyy_zzzz_0[j] = pa_y[j] * tdx_yyy_zzzz_0[j] + 1.5 * fl1_fx * tdx_yy_zzzz_0[j];

            tdy_yyyy_zzzz_0[j] = pa_y[j] * tdy_yyy_zzzz_0[j] + 1.5 * fl1_fx * tdy_yy_zzzz_0[j] + 0.5 * fl1_fx * ts_yyy_zzzz_0[j];

            tdz_yyyy_zzzz_0[j] = pa_y[j] * tdz_yyy_zzzz_0[j] + 1.5 * fl1_fx * tdz_yy_zzzz_0[j];

            tdx_yyyz_xxxx_0[j] = pa_y[j] * tdx_yyz_xxxx_0[j] + fl1_fx * tdx_yz_xxxx_0[j];

            tdy_yyyz_xxxx_0[j] = pa_y[j] * tdy_yyz_xxxx_0[j] + fl1_fx * tdy_yz_xxxx_0[j] + 0.5 * fl1_fx * ts_yyz_xxxx_0[j];

            tdz_yyyz_xxxx_0[j] = pa_y[j] * tdz_yyz_xxxx_0[j] + fl1_fx * tdz_yz_xxxx_0[j];

            tdx_yyyz_xxxy_0[j] = pa_y[j] * tdx_yyz_xxxy_0[j] + fl1_fx * tdx_yz_xxxy_0[j] + 0.5 * fl1_fx * tdx_yyz_xxx_0[j];

            tdy_yyyz_xxxy_0[j] =
                pa_y[j] * tdy_yyz_xxxy_0[j] + fl1_fx * tdy_yz_xxxy_0[j] + 0.5 * fl1_fx * tdy_yyz_xxx_0[j] + 0.5 * fl1_fx * ts_yyz_xxxy_0[j];

            tdz_yyyz_xxxy_0[j] = pa_y[j] * tdz_yyz_xxxy_0[j] + fl1_fx * tdz_yz_xxxy_0[j] + 0.5 * fl1_fx * tdz_yyz_xxx_0[j];

            tdx_yyyz_xxxz_0[j] = pa_y[j] * tdx_yyz_xxxz_0[j] + fl1_fx * tdx_yz_xxxz_0[j];

            tdy_yyyz_xxxz_0[j] = pa_y[j] * tdy_yyz_xxxz_0[j] + fl1_fx * tdy_yz_xxxz_0[j] + 0.5 * fl1_fx * ts_yyz_xxxz_0[j];

            tdz_yyyz_xxxz_0[j] = pa_y[j] * tdz_yyz_xxxz_0[j] + fl1_fx * tdz_yz_xxxz_0[j];

            tdx_yyyz_xxyy_0[j] = pa_y[j] * tdx_yyz_xxyy_0[j] + fl1_fx * tdx_yz_xxyy_0[j] + fl1_fx * tdx_yyz_xxy_0[j];

            tdy_yyyz_xxyy_0[j] =
                pa_y[j] * tdy_yyz_xxyy_0[j] + fl1_fx * tdy_yz_xxyy_0[j] + fl1_fx * tdy_yyz_xxy_0[j] + 0.5 * fl1_fx * ts_yyz_xxyy_0[j];

            tdz_yyyz_xxyy_0[j] = pa_y[j] * tdz_yyz_xxyy_0[j] + fl1_fx * tdz_yz_xxyy_0[j] + fl1_fx * tdz_yyz_xxy_0[j];

            tdx_yyyz_xxyz_0[j] = pa_y[j] * tdx_yyz_xxyz_0[j] + fl1_fx * tdx_yz_xxyz_0[j] + 0.5 * fl1_fx * tdx_yyz_xxz_0[j];

            tdy_yyyz_xxyz_0[j] =
                pa_y[j] * tdy_yyz_xxyz_0[j] + fl1_fx * tdy_yz_xxyz_0[j] + 0.5 * fl1_fx * tdy_yyz_xxz_0[j] + 0.5 * fl1_fx * ts_yyz_xxyz_0[j];

            tdz_yyyz_xxyz_0[j] = pa_y[j] * tdz_yyz_xxyz_0[j] + fl1_fx * tdz_yz_xxyz_0[j] + 0.5 * fl1_fx * tdz_yyz_xxz_0[j];

            tdx_yyyz_xxzz_0[j] = pa_y[j] * tdx_yyz_xxzz_0[j] + fl1_fx * tdx_yz_xxzz_0[j];

            tdy_yyyz_xxzz_0[j] = pa_y[j] * tdy_yyz_xxzz_0[j] + fl1_fx * tdy_yz_xxzz_0[j] + 0.5 * fl1_fx * ts_yyz_xxzz_0[j];

            tdz_yyyz_xxzz_0[j] = pa_y[j] * tdz_yyz_xxzz_0[j] + fl1_fx * tdz_yz_xxzz_0[j];

            tdx_yyyz_xyyy_0[j] = pa_y[j] * tdx_yyz_xyyy_0[j] + fl1_fx * tdx_yz_xyyy_0[j] + 1.5 * fl1_fx * tdx_yyz_xyy_0[j];

            tdy_yyyz_xyyy_0[j] =
                pa_y[j] * tdy_yyz_xyyy_0[j] + fl1_fx * tdy_yz_xyyy_0[j] + 1.5 * fl1_fx * tdy_yyz_xyy_0[j] + 0.5 * fl1_fx * ts_yyz_xyyy_0[j];

            tdz_yyyz_xyyy_0[j] = pa_y[j] * tdz_yyz_xyyy_0[j] + fl1_fx * tdz_yz_xyyy_0[j] + 1.5 * fl1_fx * tdz_yyz_xyy_0[j];

            tdx_yyyz_xyyz_0[j] = pa_y[j] * tdx_yyz_xyyz_0[j] + fl1_fx * tdx_yz_xyyz_0[j] + fl1_fx * tdx_yyz_xyz_0[j];

            tdy_yyyz_xyyz_0[j] =
                pa_y[j] * tdy_yyz_xyyz_0[j] + fl1_fx * tdy_yz_xyyz_0[j] + fl1_fx * tdy_yyz_xyz_0[j] + 0.5 * fl1_fx * ts_yyz_xyyz_0[j];

            tdz_yyyz_xyyz_0[j] = pa_y[j] * tdz_yyz_xyyz_0[j] + fl1_fx * tdz_yz_xyyz_0[j] + fl1_fx * tdz_yyz_xyz_0[j];

            tdx_yyyz_xyzz_0[j] = pa_y[j] * tdx_yyz_xyzz_0[j] + fl1_fx * tdx_yz_xyzz_0[j] + 0.5 * fl1_fx * tdx_yyz_xzz_0[j];

            tdy_yyyz_xyzz_0[j] =
                pa_y[j] * tdy_yyz_xyzz_0[j] + fl1_fx * tdy_yz_xyzz_0[j] + 0.5 * fl1_fx * tdy_yyz_xzz_0[j] + 0.5 * fl1_fx * ts_yyz_xyzz_0[j];

            tdz_yyyz_xyzz_0[j] = pa_y[j] * tdz_yyz_xyzz_0[j] + fl1_fx * tdz_yz_xyzz_0[j] + 0.5 * fl1_fx * tdz_yyz_xzz_0[j];

            tdx_yyyz_xzzz_0[j] = pa_y[j] * tdx_yyz_xzzz_0[j] + fl1_fx * tdx_yz_xzzz_0[j];

            tdy_yyyz_xzzz_0[j] = pa_y[j] * tdy_yyz_xzzz_0[j] + fl1_fx * tdy_yz_xzzz_0[j] + 0.5 * fl1_fx * ts_yyz_xzzz_0[j];

            tdz_yyyz_xzzz_0[j] = pa_y[j] * tdz_yyz_xzzz_0[j] + fl1_fx * tdz_yz_xzzz_0[j];

            tdx_yyyz_yyyy_0[j] = pa_y[j] * tdx_yyz_yyyy_0[j] + fl1_fx * tdx_yz_yyyy_0[j] + 2.0 * fl1_fx * tdx_yyz_yyy_0[j];

            tdy_yyyz_yyyy_0[j] =
                pa_y[j] * tdy_yyz_yyyy_0[j] + fl1_fx * tdy_yz_yyyy_0[j] + 2.0 * fl1_fx * tdy_yyz_yyy_0[j] + 0.5 * fl1_fx * ts_yyz_yyyy_0[j];

            tdz_yyyz_yyyy_0[j] = pa_y[j] * tdz_yyz_yyyy_0[j] + fl1_fx * tdz_yz_yyyy_0[j] + 2.0 * fl1_fx * tdz_yyz_yyy_0[j];

            tdx_yyyz_yyyz_0[j] = pa_y[j] * tdx_yyz_yyyz_0[j] + fl1_fx * tdx_yz_yyyz_0[j] + 1.5 * fl1_fx * tdx_yyz_yyz_0[j];

            tdy_yyyz_yyyz_0[j] =
                pa_y[j] * tdy_yyz_yyyz_0[j] + fl1_fx * tdy_yz_yyyz_0[j] + 1.5 * fl1_fx * tdy_yyz_yyz_0[j] + 0.5 * fl1_fx * ts_yyz_yyyz_0[j];

            tdz_yyyz_yyyz_0[j] = pa_y[j] * tdz_yyz_yyyz_0[j] + fl1_fx * tdz_yz_yyyz_0[j] + 1.5 * fl1_fx * tdz_yyz_yyz_0[j];
        }

        idx++;
    }
}

void
compElectricDipoleForGG_531_579(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const int32_t              nOSFactors,
                                const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_d_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_d_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_y = paDistances.data(3 * idx + 1);

        // set up pointers to auxilary integrals

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

        auto tdx_yyz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 78);

        auto tdy_yyz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 78);

        auto tdz_yyz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 78);

        auto tdx_yyz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 79);

        auto tdy_yyz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 79);

        auto tdz_yyz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 79);

        auto tdx_yzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 80);

        auto tdy_yzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 80);

        auto tdz_yzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 80);

        auto tdx_yzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 81);

        auto tdy_yzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 81);

        auto tdz_yzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 81);

        auto tdx_yzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 82);

        auto tdy_yzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 82);

        auto tdz_yzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 82);

        auto tdx_yzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 83);

        auto tdy_yzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 83);

        auto tdz_yzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 83);

        auto tdx_yzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 84);

        auto tdy_yzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 84);

        auto tdz_yzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 84);

        auto tdx_yzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 85);

        auto tdy_yzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 85);

        auto tdz_yzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 85);

        auto tdx_yzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 86);

        auto tdy_yzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 86);

        auto tdz_yzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 86);

        auto tdx_yzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 87);

        auto tdy_yzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 87);

        auto tdz_yzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 87);

        auto tdx_yzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 88);

        auto tdy_yzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 88);

        auto tdz_yzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 88);

        auto ts_yyz_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 117);

        auto ts_yyz_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 118);

        auto ts_yyz_zzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 119);

        auto ts_yzz_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 120);

        auto ts_yzz_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 121);

        auto ts_yzz_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 122);

        auto ts_yzz_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 123);

        auto ts_yzz_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 124);

        auto ts_yzz_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 125);

        auto ts_yzz_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 126);

        auto ts_yzz_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 127);

        auto ts_yzz_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 128);

        auto ts_yzz_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 129);

        auto ts_yzz_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 130);

        auto ts_yzz_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 131);

        auto ts_yzz_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 132);

        // set up pointers to integrals

        auto tdx_yyyz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 177);

        auto tdy_yyyz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 177);

        auto tdz_yyyz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 177);

        auto tdx_yyyz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 178);

        auto tdy_yyyz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 178);

        auto tdz_yyyz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 178);

        auto tdx_yyyz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 179);

        auto tdy_yyyz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 179);

        auto tdz_yyyz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 179);

        auto tdx_yyzz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 180);

        auto tdy_yyzz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 180);

        auto tdz_yyzz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 180);

        auto tdx_yyzz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 181);

        auto tdy_yyzz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 181);

        auto tdz_yyzz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 181);

        auto tdx_yyzz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 182);

        auto tdy_yyzz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 182);

        auto tdz_yyzz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 182);

        auto tdx_yyzz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 183);

        auto tdy_yyzz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 183);

        auto tdz_yyzz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 183);

        auto tdx_yyzz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 184);

        auto tdy_yyzz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 184);

        auto tdz_yyzz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 184);

        auto tdx_yyzz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 185);

        auto tdy_yyzz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 185);

        auto tdz_yyzz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 185);

        auto tdx_yyzz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 186);

        auto tdy_yyzz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 186);

        auto tdz_yyzz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 186);

        auto tdx_yyzz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 187);

        auto tdy_yyzz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 187);

        auto tdz_yyzz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 187);

        auto tdx_yyzz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 188);

        auto tdy_yyzz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 188);

        auto tdz_yyzz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 188);

        auto tdx_yyzz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 189);

        auto tdy_yyzz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 189);

        auto tdz_yyzz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 189);

        auto tdx_yyzz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 190);

        auto tdy_yyzz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 190);

        auto tdz_yyzz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 190);

        auto tdx_yyzz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 191);

        auto tdy_yyzz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 191);

        auto tdz_yyzz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 191);

        auto tdx_yyzz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 192);

        auto tdy_yyzz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 192);

        auto tdz_yyzz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 192);

        // Batch of Integrals (531,579)

        #pragma omp simd aligned(fx, pa_y, tdx_yyyz_yyzz_0, tdx_yyyz_yzzz_0, tdx_yyyz_zzzz_0, \
                                     tdx_yyz_yyzz_0, tdx_yyz_yzz_0, tdx_yyz_yzzz_0, tdx_yyz_zzz_0, tdx_yyz_zzzz_0, \
                                     tdx_yyzz_xxxx_0, tdx_yyzz_xxxy_0, tdx_yyzz_xxxz_0, tdx_yyzz_xxyy_0, tdx_yyzz_xxyz_0, \
                                     tdx_yyzz_xxzz_0, tdx_yyzz_xyyy_0, tdx_yyzz_xyyz_0, tdx_yyzz_xyzz_0, tdx_yyzz_xzzz_0, \
                                     tdx_yyzz_yyyy_0, tdx_yyzz_yyyz_0, tdx_yyzz_yyzz_0, tdx_yz_yyzz_0, tdx_yz_yzzz_0, \
                                     tdx_yz_zzzz_0, tdx_yzz_xxx_0, tdx_yzz_xxxx_0, tdx_yzz_xxxy_0, tdx_yzz_xxxz_0, \
                                     tdx_yzz_xxy_0, tdx_yzz_xxyy_0, tdx_yzz_xxyz_0, tdx_yzz_xxz_0, tdx_yzz_xxzz_0, \
                                     tdx_yzz_xyy_0, tdx_yzz_xyyy_0, tdx_yzz_xyyz_0, tdx_yzz_xyz_0, tdx_yzz_xyzz_0, \
                                     tdx_yzz_xzz_0, tdx_yzz_xzzz_0, tdx_yzz_yyy_0, tdx_yzz_yyyy_0, tdx_yzz_yyyz_0, \
                                     tdx_yzz_yyz_0, tdx_yzz_yyzz_0, tdx_yzz_yzz_0, tdx_zz_xxxx_0, tdx_zz_xxxy_0, \
                                     tdx_zz_xxxz_0, tdx_zz_xxyy_0, tdx_zz_xxyz_0, tdx_zz_xxzz_0, tdx_zz_xyyy_0, \
                                     tdx_zz_xyyz_0, tdx_zz_xyzz_0, tdx_zz_xzzz_0, tdx_zz_yyyy_0, tdx_zz_yyyz_0, \
                                     tdx_zz_yyzz_0, tdy_yyyz_yyzz_0, tdy_yyyz_yzzz_0, tdy_yyyz_zzzz_0, tdy_yyz_yyzz_0, \
                                     tdy_yyz_yzz_0, tdy_yyz_yzzz_0, tdy_yyz_zzz_0, tdy_yyz_zzzz_0, tdy_yyzz_xxxx_0, \
                                     tdy_yyzz_xxxy_0, tdy_yyzz_xxxz_0, tdy_yyzz_xxyy_0, tdy_yyzz_xxyz_0, tdy_yyzz_xxzz_0, \
                                     tdy_yyzz_xyyy_0, tdy_yyzz_xyyz_0, tdy_yyzz_xyzz_0, tdy_yyzz_xzzz_0, tdy_yyzz_yyyy_0, \
                                     tdy_yyzz_yyyz_0, tdy_yyzz_yyzz_0, tdy_yz_yyzz_0, tdy_yz_yzzz_0, tdy_yz_zzzz_0, \
                                     tdy_yzz_xxx_0, tdy_yzz_xxxx_0, tdy_yzz_xxxy_0, tdy_yzz_xxxz_0, tdy_yzz_xxy_0, \
                                     tdy_yzz_xxyy_0, tdy_yzz_xxyz_0, tdy_yzz_xxz_0, tdy_yzz_xxzz_0, tdy_yzz_xyy_0, \
                                     tdy_yzz_xyyy_0, tdy_yzz_xyyz_0, tdy_yzz_xyz_0, tdy_yzz_xyzz_0, tdy_yzz_xzz_0, \
                                     tdy_yzz_xzzz_0, tdy_yzz_yyy_0, tdy_yzz_yyyy_0, tdy_yzz_yyyz_0, tdy_yzz_yyz_0, \
                                     tdy_yzz_yyzz_0, tdy_yzz_yzz_0, tdy_zz_xxxx_0, tdy_zz_xxxy_0, tdy_zz_xxxz_0, \
                                     tdy_zz_xxyy_0, tdy_zz_xxyz_0, tdy_zz_xxzz_0, tdy_zz_xyyy_0, tdy_zz_xyyz_0, \
                                     tdy_zz_xyzz_0, tdy_zz_xzzz_0, tdy_zz_yyyy_0, tdy_zz_yyyz_0, tdy_zz_yyzz_0, \
                                     tdz_yyyz_yyzz_0, tdz_yyyz_yzzz_0, tdz_yyyz_zzzz_0, tdz_yyz_yyzz_0, tdz_yyz_yzz_0, \
                                     tdz_yyz_yzzz_0, tdz_yyz_zzz_0, tdz_yyz_zzzz_0, tdz_yyzz_xxxx_0, tdz_yyzz_xxxy_0, \
                                     tdz_yyzz_xxxz_0, tdz_yyzz_xxyy_0, tdz_yyzz_xxyz_0, tdz_yyzz_xxzz_0, tdz_yyzz_xyyy_0, \
                                     tdz_yyzz_xyyz_0, tdz_yyzz_xyzz_0, tdz_yyzz_xzzz_0, tdz_yyzz_yyyy_0, tdz_yyzz_yyyz_0, \
                                     tdz_yyzz_yyzz_0, tdz_yz_yyzz_0, tdz_yz_yzzz_0, tdz_yz_zzzz_0, tdz_yzz_xxx_0, \
                                     tdz_yzz_xxxx_0, tdz_yzz_xxxy_0, tdz_yzz_xxxz_0, tdz_yzz_xxy_0, tdz_yzz_xxyy_0, \
                                     tdz_yzz_xxyz_0, tdz_yzz_xxz_0, tdz_yzz_xxzz_0, tdz_yzz_xyy_0, tdz_yzz_xyyy_0, \
                                     tdz_yzz_xyyz_0, tdz_yzz_xyz_0, tdz_yzz_xyzz_0, tdz_yzz_xzz_0, tdz_yzz_xzzz_0, \
                                     tdz_yzz_yyy_0, tdz_yzz_yyyy_0, tdz_yzz_yyyz_0, tdz_yzz_yyz_0, tdz_yzz_yyzz_0, \
                                     tdz_yzz_yzz_0, tdz_zz_xxxx_0, tdz_zz_xxxy_0, tdz_zz_xxxz_0, tdz_zz_xxyy_0, \
                                     tdz_zz_xxyz_0, tdz_zz_xxzz_0, tdz_zz_xyyy_0, tdz_zz_xyyz_0, tdz_zz_xyzz_0, \
                                     tdz_zz_xzzz_0, tdz_zz_yyyy_0, tdz_zz_yyyz_0, tdz_zz_yyzz_0, ts_yyz_yyzz_0, \
                                     ts_yyz_yzzz_0, ts_yyz_zzzz_0, ts_yzz_xxxx_0, ts_yzz_xxxy_0, ts_yzz_xxxz_0, \
                                     ts_yzz_xxyy_0, ts_yzz_xxyz_0, ts_yzz_xxzz_0, ts_yzz_xyyy_0, ts_yzz_xyyz_0, \
                                     ts_yzz_xyzz_0, ts_yzz_xzzz_0, ts_yzz_yyyy_0, ts_yzz_yyyz_0, ts_yzz_yyzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            tdx_yyyz_yyzz_0[j] = pa_y[j] * tdx_yyz_yyzz_0[j] + fl1_fx * tdx_yz_yyzz_0[j] + fl1_fx * tdx_yyz_yzz_0[j];

            tdy_yyyz_yyzz_0[j] =
                pa_y[j] * tdy_yyz_yyzz_0[j] + fl1_fx * tdy_yz_yyzz_0[j] + fl1_fx * tdy_yyz_yzz_0[j] + 0.5 * fl1_fx * ts_yyz_yyzz_0[j];

            tdz_yyyz_yyzz_0[j] = pa_y[j] * tdz_yyz_yyzz_0[j] + fl1_fx * tdz_yz_yyzz_0[j] + fl1_fx * tdz_yyz_yzz_0[j];

            tdx_yyyz_yzzz_0[j] = pa_y[j] * tdx_yyz_yzzz_0[j] + fl1_fx * tdx_yz_yzzz_0[j] + 0.5 * fl1_fx * tdx_yyz_zzz_0[j];

            tdy_yyyz_yzzz_0[j] =
                pa_y[j] * tdy_yyz_yzzz_0[j] + fl1_fx * tdy_yz_yzzz_0[j] + 0.5 * fl1_fx * tdy_yyz_zzz_0[j] + 0.5 * fl1_fx * ts_yyz_yzzz_0[j];

            tdz_yyyz_yzzz_0[j] = pa_y[j] * tdz_yyz_yzzz_0[j] + fl1_fx * tdz_yz_yzzz_0[j] + 0.5 * fl1_fx * tdz_yyz_zzz_0[j];

            tdx_yyyz_zzzz_0[j] = pa_y[j] * tdx_yyz_zzzz_0[j] + fl1_fx * tdx_yz_zzzz_0[j];

            tdy_yyyz_zzzz_0[j] = pa_y[j] * tdy_yyz_zzzz_0[j] + fl1_fx * tdy_yz_zzzz_0[j] + 0.5 * fl1_fx * ts_yyz_zzzz_0[j];

            tdz_yyyz_zzzz_0[j] = pa_y[j] * tdz_yyz_zzzz_0[j] + fl1_fx * tdz_yz_zzzz_0[j];

            tdx_yyzz_xxxx_0[j] = pa_y[j] * tdx_yzz_xxxx_0[j] + 0.5 * fl1_fx * tdx_zz_xxxx_0[j];

            tdy_yyzz_xxxx_0[j] = pa_y[j] * tdy_yzz_xxxx_0[j] + 0.5 * fl1_fx * tdy_zz_xxxx_0[j] + 0.5 * fl1_fx * ts_yzz_xxxx_0[j];

            tdz_yyzz_xxxx_0[j] = pa_y[j] * tdz_yzz_xxxx_0[j] + 0.5 * fl1_fx * tdz_zz_xxxx_0[j];

            tdx_yyzz_xxxy_0[j] = pa_y[j] * tdx_yzz_xxxy_0[j] + 0.5 * fl1_fx * tdx_zz_xxxy_0[j] + 0.5 * fl1_fx * tdx_yzz_xxx_0[j];

            tdy_yyzz_xxxy_0[j] =
                pa_y[j] * tdy_yzz_xxxy_0[j] + 0.5 * fl1_fx * tdy_zz_xxxy_0[j] + 0.5 * fl1_fx * tdy_yzz_xxx_0[j] + 0.5 * fl1_fx * ts_yzz_xxxy_0[j];

            tdz_yyzz_xxxy_0[j] = pa_y[j] * tdz_yzz_xxxy_0[j] + 0.5 * fl1_fx * tdz_zz_xxxy_0[j] + 0.5 * fl1_fx * tdz_yzz_xxx_0[j];

            tdx_yyzz_xxxz_0[j] = pa_y[j] * tdx_yzz_xxxz_0[j] + 0.5 * fl1_fx * tdx_zz_xxxz_0[j];

            tdy_yyzz_xxxz_0[j] = pa_y[j] * tdy_yzz_xxxz_0[j] + 0.5 * fl1_fx * tdy_zz_xxxz_0[j] + 0.5 * fl1_fx * ts_yzz_xxxz_0[j];

            tdz_yyzz_xxxz_0[j] = pa_y[j] * tdz_yzz_xxxz_0[j] + 0.5 * fl1_fx * tdz_zz_xxxz_0[j];

            tdx_yyzz_xxyy_0[j] = pa_y[j] * tdx_yzz_xxyy_0[j] + 0.5 * fl1_fx * tdx_zz_xxyy_0[j] + fl1_fx * tdx_yzz_xxy_0[j];

            tdy_yyzz_xxyy_0[j] =
                pa_y[j] * tdy_yzz_xxyy_0[j] + 0.5 * fl1_fx * tdy_zz_xxyy_0[j] + fl1_fx * tdy_yzz_xxy_0[j] + 0.5 * fl1_fx * ts_yzz_xxyy_0[j];

            tdz_yyzz_xxyy_0[j] = pa_y[j] * tdz_yzz_xxyy_0[j] + 0.5 * fl1_fx * tdz_zz_xxyy_0[j] + fl1_fx * tdz_yzz_xxy_0[j];

            tdx_yyzz_xxyz_0[j] = pa_y[j] * tdx_yzz_xxyz_0[j] + 0.5 * fl1_fx * tdx_zz_xxyz_0[j] + 0.5 * fl1_fx * tdx_yzz_xxz_0[j];

            tdy_yyzz_xxyz_0[j] =
                pa_y[j] * tdy_yzz_xxyz_0[j] + 0.5 * fl1_fx * tdy_zz_xxyz_0[j] + 0.5 * fl1_fx * tdy_yzz_xxz_0[j] + 0.5 * fl1_fx * ts_yzz_xxyz_0[j];

            tdz_yyzz_xxyz_0[j] = pa_y[j] * tdz_yzz_xxyz_0[j] + 0.5 * fl1_fx * tdz_zz_xxyz_0[j] + 0.5 * fl1_fx * tdz_yzz_xxz_0[j];

            tdx_yyzz_xxzz_0[j] = pa_y[j] * tdx_yzz_xxzz_0[j] + 0.5 * fl1_fx * tdx_zz_xxzz_0[j];

            tdy_yyzz_xxzz_0[j] = pa_y[j] * tdy_yzz_xxzz_0[j] + 0.5 * fl1_fx * tdy_zz_xxzz_0[j] + 0.5 * fl1_fx * ts_yzz_xxzz_0[j];

            tdz_yyzz_xxzz_0[j] = pa_y[j] * tdz_yzz_xxzz_0[j] + 0.5 * fl1_fx * tdz_zz_xxzz_0[j];

            tdx_yyzz_xyyy_0[j] = pa_y[j] * tdx_yzz_xyyy_0[j] + 0.5 * fl1_fx * tdx_zz_xyyy_0[j] + 1.5 * fl1_fx * tdx_yzz_xyy_0[j];

            tdy_yyzz_xyyy_0[j] =
                pa_y[j] * tdy_yzz_xyyy_0[j] + 0.5 * fl1_fx * tdy_zz_xyyy_0[j] + 1.5 * fl1_fx * tdy_yzz_xyy_0[j] + 0.5 * fl1_fx * ts_yzz_xyyy_0[j];

            tdz_yyzz_xyyy_0[j] = pa_y[j] * tdz_yzz_xyyy_0[j] + 0.5 * fl1_fx * tdz_zz_xyyy_0[j] + 1.5 * fl1_fx * tdz_yzz_xyy_0[j];

            tdx_yyzz_xyyz_0[j] = pa_y[j] * tdx_yzz_xyyz_0[j] + 0.5 * fl1_fx * tdx_zz_xyyz_0[j] + fl1_fx * tdx_yzz_xyz_0[j];

            tdy_yyzz_xyyz_0[j] =
                pa_y[j] * tdy_yzz_xyyz_0[j] + 0.5 * fl1_fx * tdy_zz_xyyz_0[j] + fl1_fx * tdy_yzz_xyz_0[j] + 0.5 * fl1_fx * ts_yzz_xyyz_0[j];

            tdz_yyzz_xyyz_0[j] = pa_y[j] * tdz_yzz_xyyz_0[j] + 0.5 * fl1_fx * tdz_zz_xyyz_0[j] + fl1_fx * tdz_yzz_xyz_0[j];

            tdx_yyzz_xyzz_0[j] = pa_y[j] * tdx_yzz_xyzz_0[j] + 0.5 * fl1_fx * tdx_zz_xyzz_0[j] + 0.5 * fl1_fx * tdx_yzz_xzz_0[j];

            tdy_yyzz_xyzz_0[j] =
                pa_y[j] * tdy_yzz_xyzz_0[j] + 0.5 * fl1_fx * tdy_zz_xyzz_0[j] + 0.5 * fl1_fx * tdy_yzz_xzz_0[j] + 0.5 * fl1_fx * ts_yzz_xyzz_0[j];

            tdz_yyzz_xyzz_0[j] = pa_y[j] * tdz_yzz_xyzz_0[j] + 0.5 * fl1_fx * tdz_zz_xyzz_0[j] + 0.5 * fl1_fx * tdz_yzz_xzz_0[j];

            tdx_yyzz_xzzz_0[j] = pa_y[j] * tdx_yzz_xzzz_0[j] + 0.5 * fl1_fx * tdx_zz_xzzz_0[j];

            tdy_yyzz_xzzz_0[j] = pa_y[j] * tdy_yzz_xzzz_0[j] + 0.5 * fl1_fx * tdy_zz_xzzz_0[j] + 0.5 * fl1_fx * ts_yzz_xzzz_0[j];

            tdz_yyzz_xzzz_0[j] = pa_y[j] * tdz_yzz_xzzz_0[j] + 0.5 * fl1_fx * tdz_zz_xzzz_0[j];

            tdx_yyzz_yyyy_0[j] = pa_y[j] * tdx_yzz_yyyy_0[j] + 0.5 * fl1_fx * tdx_zz_yyyy_0[j] + 2.0 * fl1_fx * tdx_yzz_yyy_0[j];

            tdy_yyzz_yyyy_0[j] =
                pa_y[j] * tdy_yzz_yyyy_0[j] + 0.5 * fl1_fx * tdy_zz_yyyy_0[j] + 2.0 * fl1_fx * tdy_yzz_yyy_0[j] + 0.5 * fl1_fx * ts_yzz_yyyy_0[j];

            tdz_yyzz_yyyy_0[j] = pa_y[j] * tdz_yzz_yyyy_0[j] + 0.5 * fl1_fx * tdz_zz_yyyy_0[j] + 2.0 * fl1_fx * tdz_yzz_yyy_0[j];

            tdx_yyzz_yyyz_0[j] = pa_y[j] * tdx_yzz_yyyz_0[j] + 0.5 * fl1_fx * tdx_zz_yyyz_0[j] + 1.5 * fl1_fx * tdx_yzz_yyz_0[j];

            tdy_yyzz_yyyz_0[j] =
                pa_y[j] * tdy_yzz_yyyz_0[j] + 0.5 * fl1_fx * tdy_zz_yyyz_0[j] + 1.5 * fl1_fx * tdy_yzz_yyz_0[j] + 0.5 * fl1_fx * ts_yzz_yyyz_0[j];

            tdz_yyzz_yyyz_0[j] = pa_y[j] * tdz_yzz_yyyz_0[j] + 0.5 * fl1_fx * tdz_zz_yyyz_0[j] + 1.5 * fl1_fx * tdz_yzz_yyz_0[j];

            tdx_yyzz_yyzz_0[j] = pa_y[j] * tdx_yzz_yyzz_0[j] + 0.5 * fl1_fx * tdx_zz_yyzz_0[j] + fl1_fx * tdx_yzz_yzz_0[j];

            tdy_yyzz_yyzz_0[j] =
                pa_y[j] * tdy_yzz_yyzz_0[j] + 0.5 * fl1_fx * tdy_zz_yyzz_0[j] + fl1_fx * tdy_yzz_yzz_0[j] + 0.5 * fl1_fx * ts_yzz_yyzz_0[j];

            tdz_yyzz_yyzz_0[j] = pa_y[j] * tdz_yzz_yyzz_0[j] + 0.5 * fl1_fx * tdz_zz_yyzz_0[j] + fl1_fx * tdz_yzz_yzz_0[j];
        }

        idx++;
    }
}

void
compElectricDipoleForGG_579_627(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const int32_t              nOSFactors,
                                const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_d_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_d_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_y = paDistances.data(3 * idx + 1);

        // set up pointers to auxilary integrals

        auto tdx_yzz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 133);

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

        auto tdx_zz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 88);

        auto tdy_zz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 88);

        auto tdz_zz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 88);

        auto tdx_zz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 89);

        auto tdy_zz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 89);

        auto tdz_zz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 89);

        auto tdx_yzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 89);

        auto tdy_yzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 89);

        auto tdz_yzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 89);

        auto tdx_zzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 90);

        auto tdy_zzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 90);

        auto tdz_zzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 90);

        auto tdx_zzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 91);

        auto tdy_zzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 91);

        auto tdz_zzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 91);

        auto tdx_zzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 92);

        auto tdy_zzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 92);

        auto tdz_zzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 92);

        auto tdx_zzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 93);

        auto tdy_zzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 93);

        auto tdz_zzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 93);

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

        auto ts_yzz_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 133);

        auto ts_yzz_zzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 134);

        auto ts_zzz_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 135);

        auto ts_zzz_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 136);

        auto ts_zzz_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 137);

        auto ts_zzz_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 138);

        auto ts_zzz_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 139);

        auto ts_zzz_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 140);

        auto ts_zzz_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 141);

        auto ts_zzz_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 142);

        auto ts_zzz_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 143);

        auto ts_zzz_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 144);

        auto ts_zzz_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 145);

        auto ts_zzz_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 146);

        auto ts_zzz_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 147);

        auto ts_zzz_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 148);

        // set up pointers to integrals

        auto tdx_yyzz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 193);

        auto tdy_yyzz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 193);

        auto tdz_yyzz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 193);

        auto tdx_yyzz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 194);

        auto tdy_yyzz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 194);

        auto tdz_yyzz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 194);

        auto tdx_yzzz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 195);

        auto tdy_yzzz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 195);

        auto tdz_yzzz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 195);

        auto tdx_yzzz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 196);

        auto tdy_yzzz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 196);

        auto tdz_yzzz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 196);

        auto tdx_yzzz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 197);

        auto tdy_yzzz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 197);

        auto tdz_yzzz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 197);

        auto tdx_yzzz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 198);

        auto tdy_yzzz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 198);

        auto tdz_yzzz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 198);

        auto tdx_yzzz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 199);

        auto tdy_yzzz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 199);

        auto tdz_yzzz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 199);

        auto tdx_yzzz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 200);

        auto tdy_yzzz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 200);

        auto tdz_yzzz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 200);

        auto tdx_yzzz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 201);

        auto tdy_yzzz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 201);

        auto tdz_yzzz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 201);

        auto tdx_yzzz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 202);

        auto tdy_yzzz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 202);

        auto tdz_yzzz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 202);

        auto tdx_yzzz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 203);

        auto tdy_yzzz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 203);

        auto tdz_yzzz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 203);

        auto tdx_yzzz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 204);

        auto tdy_yzzz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 204);

        auto tdz_yzzz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 204);

        auto tdx_yzzz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 205);

        auto tdy_yzzz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 205);

        auto tdz_yzzz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 205);

        auto tdx_yzzz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 206);

        auto tdy_yzzz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 206);

        auto tdz_yzzz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 206);

        auto tdx_yzzz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 207);

        auto tdy_yzzz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 207);

        auto tdz_yzzz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 207);

        auto tdx_yzzz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 208);

        auto tdy_yzzz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 208);

        auto tdz_yzzz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 208);

        // Batch of Integrals (579,627)

        #pragma omp simd aligned(fx, pa_y, tdx_yyzz_yzzz_0, tdx_yyzz_zzzz_0, tdx_yzz_yzzz_0, \
                                     tdx_yzz_zzz_0, tdx_yzz_zzzz_0, tdx_yzzz_xxxx_0, tdx_yzzz_xxxy_0, tdx_yzzz_xxxz_0, \
                                     tdx_yzzz_xxyy_0, tdx_yzzz_xxyz_0, tdx_yzzz_xxzz_0, tdx_yzzz_xyyy_0, tdx_yzzz_xyyz_0, \
                                     tdx_yzzz_xyzz_0, tdx_yzzz_xzzz_0, tdx_yzzz_yyyy_0, tdx_yzzz_yyyz_0, tdx_yzzz_yyzz_0, \
                                     tdx_yzzz_yzzz_0, tdx_zz_yzzz_0, tdx_zz_zzzz_0, tdx_zzz_xxx_0, tdx_zzz_xxxx_0, \
                                     tdx_zzz_xxxy_0, tdx_zzz_xxxz_0, tdx_zzz_xxy_0, tdx_zzz_xxyy_0, tdx_zzz_xxyz_0, \
                                     tdx_zzz_xxz_0, tdx_zzz_xxzz_0, tdx_zzz_xyy_0, tdx_zzz_xyyy_0, tdx_zzz_xyyz_0, \
                                     tdx_zzz_xyz_0, tdx_zzz_xyzz_0, tdx_zzz_xzz_0, tdx_zzz_xzzz_0, tdx_zzz_yyy_0, \
                                     tdx_zzz_yyyy_0, tdx_zzz_yyyz_0, tdx_zzz_yyz_0, tdx_zzz_yyzz_0, tdx_zzz_yzz_0, \
                                     tdx_zzz_yzzz_0, tdx_zzz_zzz_0, tdy_yyzz_yzzz_0, tdy_yyzz_zzzz_0, tdy_yzz_yzzz_0, \
                                     tdy_yzz_zzz_0, tdy_yzz_zzzz_0, tdy_yzzz_xxxx_0, tdy_yzzz_xxxy_0, tdy_yzzz_xxxz_0, \
                                     tdy_yzzz_xxyy_0, tdy_yzzz_xxyz_0, tdy_yzzz_xxzz_0, tdy_yzzz_xyyy_0, tdy_yzzz_xyyz_0, \
                                     tdy_yzzz_xyzz_0, tdy_yzzz_xzzz_0, tdy_yzzz_yyyy_0, tdy_yzzz_yyyz_0, tdy_yzzz_yyzz_0, \
                                     tdy_yzzz_yzzz_0, tdy_zz_yzzz_0, tdy_zz_zzzz_0, tdy_zzz_xxx_0, tdy_zzz_xxxx_0, \
                                     tdy_zzz_xxxy_0, tdy_zzz_xxxz_0, tdy_zzz_xxy_0, tdy_zzz_xxyy_0, tdy_zzz_xxyz_0, \
                                     tdy_zzz_xxz_0, tdy_zzz_xxzz_0, tdy_zzz_xyy_0, tdy_zzz_xyyy_0, tdy_zzz_xyyz_0, \
                                     tdy_zzz_xyz_0, tdy_zzz_xyzz_0, tdy_zzz_xzz_0, tdy_zzz_xzzz_0, tdy_zzz_yyy_0, \
                                     tdy_zzz_yyyy_0, tdy_zzz_yyyz_0, tdy_zzz_yyz_0, tdy_zzz_yyzz_0, tdy_zzz_yzz_0, \
                                     tdy_zzz_yzzz_0, tdy_zzz_zzz_0, tdz_yyzz_yzzz_0, tdz_yyzz_zzzz_0, tdz_yzz_yzzz_0, \
                                     tdz_yzz_zzz_0, tdz_yzz_zzzz_0, tdz_yzzz_xxxx_0, tdz_yzzz_xxxy_0, tdz_yzzz_xxxz_0, \
                                     tdz_yzzz_xxyy_0, tdz_yzzz_xxyz_0, tdz_yzzz_xxzz_0, tdz_yzzz_xyyy_0, tdz_yzzz_xyyz_0, \
                                     tdz_yzzz_xyzz_0, tdz_yzzz_xzzz_0, tdz_yzzz_yyyy_0, tdz_yzzz_yyyz_0, tdz_yzzz_yyzz_0, \
                                     tdz_yzzz_yzzz_0, tdz_zz_yzzz_0, tdz_zz_zzzz_0, tdz_zzz_xxx_0, tdz_zzz_xxxx_0, \
                                     tdz_zzz_xxxy_0, tdz_zzz_xxxz_0, tdz_zzz_xxy_0, tdz_zzz_xxyy_0, tdz_zzz_xxyz_0, \
                                     tdz_zzz_xxz_0, tdz_zzz_xxzz_0, tdz_zzz_xyy_0, tdz_zzz_xyyy_0, tdz_zzz_xyyz_0, \
                                     tdz_zzz_xyz_0, tdz_zzz_xyzz_0, tdz_zzz_xzz_0, tdz_zzz_xzzz_0, tdz_zzz_yyy_0, \
                                     tdz_zzz_yyyy_0, tdz_zzz_yyyz_0, tdz_zzz_yyz_0, tdz_zzz_yyzz_0, tdz_zzz_yzz_0, \
                                     tdz_zzz_yzzz_0, tdz_zzz_zzz_0, ts_yzz_yzzz_0, ts_yzz_zzzz_0, ts_zzz_xxxx_0, \
                                     ts_zzz_xxxy_0, ts_zzz_xxxz_0, ts_zzz_xxyy_0, ts_zzz_xxyz_0, ts_zzz_xxzz_0, \
                                     ts_zzz_xyyy_0, ts_zzz_xyyz_0, ts_zzz_xyzz_0, ts_zzz_xzzz_0, ts_zzz_yyyy_0, \
                                     ts_zzz_yyyz_0, ts_zzz_yyzz_0, ts_zzz_yzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            tdx_yyzz_yzzz_0[j] = pa_y[j] * tdx_yzz_yzzz_0[j] + 0.5 * fl1_fx * tdx_zz_yzzz_0[j] + 0.5 * fl1_fx * tdx_yzz_zzz_0[j];

            tdy_yyzz_yzzz_0[j] =
                pa_y[j] * tdy_yzz_yzzz_0[j] + 0.5 * fl1_fx * tdy_zz_yzzz_0[j] + 0.5 * fl1_fx * tdy_yzz_zzz_0[j] + 0.5 * fl1_fx * ts_yzz_yzzz_0[j];

            tdz_yyzz_yzzz_0[j] = pa_y[j] * tdz_yzz_yzzz_0[j] + 0.5 * fl1_fx * tdz_zz_yzzz_0[j] + 0.5 * fl1_fx * tdz_yzz_zzz_0[j];

            tdx_yyzz_zzzz_0[j] = pa_y[j] * tdx_yzz_zzzz_0[j] + 0.5 * fl1_fx * tdx_zz_zzzz_0[j];

            tdy_yyzz_zzzz_0[j] = pa_y[j] * tdy_yzz_zzzz_0[j] + 0.5 * fl1_fx * tdy_zz_zzzz_0[j] + 0.5 * fl1_fx * ts_yzz_zzzz_0[j];

            tdz_yyzz_zzzz_0[j] = pa_y[j] * tdz_yzz_zzzz_0[j] + 0.5 * fl1_fx * tdz_zz_zzzz_0[j];

            tdx_yzzz_xxxx_0[j] = pa_y[j] * tdx_zzz_xxxx_0[j];

            tdy_yzzz_xxxx_0[j] = pa_y[j] * tdy_zzz_xxxx_0[j] + 0.5 * fl1_fx * ts_zzz_xxxx_0[j];

            tdz_yzzz_xxxx_0[j] = pa_y[j] * tdz_zzz_xxxx_0[j];

            tdx_yzzz_xxxy_0[j] = pa_y[j] * tdx_zzz_xxxy_0[j] + 0.5 * fl1_fx * tdx_zzz_xxx_0[j];

            tdy_yzzz_xxxy_0[j] = pa_y[j] * tdy_zzz_xxxy_0[j] + 0.5 * fl1_fx * tdy_zzz_xxx_0[j] + 0.5 * fl1_fx * ts_zzz_xxxy_0[j];

            tdz_yzzz_xxxy_0[j] = pa_y[j] * tdz_zzz_xxxy_0[j] + 0.5 * fl1_fx * tdz_zzz_xxx_0[j];

            tdx_yzzz_xxxz_0[j] = pa_y[j] * tdx_zzz_xxxz_0[j];

            tdy_yzzz_xxxz_0[j] = pa_y[j] * tdy_zzz_xxxz_0[j] + 0.5 * fl1_fx * ts_zzz_xxxz_0[j];

            tdz_yzzz_xxxz_0[j] = pa_y[j] * tdz_zzz_xxxz_0[j];

            tdx_yzzz_xxyy_0[j] = pa_y[j] * tdx_zzz_xxyy_0[j] + fl1_fx * tdx_zzz_xxy_0[j];

            tdy_yzzz_xxyy_0[j] = pa_y[j] * tdy_zzz_xxyy_0[j] + fl1_fx * tdy_zzz_xxy_0[j] + 0.5 * fl1_fx * ts_zzz_xxyy_0[j];

            tdz_yzzz_xxyy_0[j] = pa_y[j] * tdz_zzz_xxyy_0[j] + fl1_fx * tdz_zzz_xxy_0[j];

            tdx_yzzz_xxyz_0[j] = pa_y[j] * tdx_zzz_xxyz_0[j] + 0.5 * fl1_fx * tdx_zzz_xxz_0[j];

            tdy_yzzz_xxyz_0[j] = pa_y[j] * tdy_zzz_xxyz_0[j] + 0.5 * fl1_fx * tdy_zzz_xxz_0[j] + 0.5 * fl1_fx * ts_zzz_xxyz_0[j];

            tdz_yzzz_xxyz_0[j] = pa_y[j] * tdz_zzz_xxyz_0[j] + 0.5 * fl1_fx * tdz_zzz_xxz_0[j];

            tdx_yzzz_xxzz_0[j] = pa_y[j] * tdx_zzz_xxzz_0[j];

            tdy_yzzz_xxzz_0[j] = pa_y[j] * tdy_zzz_xxzz_0[j] + 0.5 * fl1_fx * ts_zzz_xxzz_0[j];

            tdz_yzzz_xxzz_0[j] = pa_y[j] * tdz_zzz_xxzz_0[j];

            tdx_yzzz_xyyy_0[j] = pa_y[j] * tdx_zzz_xyyy_0[j] + 1.5 * fl1_fx * tdx_zzz_xyy_0[j];

            tdy_yzzz_xyyy_0[j] = pa_y[j] * tdy_zzz_xyyy_0[j] + 1.5 * fl1_fx * tdy_zzz_xyy_0[j] + 0.5 * fl1_fx * ts_zzz_xyyy_0[j];

            tdz_yzzz_xyyy_0[j] = pa_y[j] * tdz_zzz_xyyy_0[j] + 1.5 * fl1_fx * tdz_zzz_xyy_0[j];

            tdx_yzzz_xyyz_0[j] = pa_y[j] * tdx_zzz_xyyz_0[j] + fl1_fx * tdx_zzz_xyz_0[j];

            tdy_yzzz_xyyz_0[j] = pa_y[j] * tdy_zzz_xyyz_0[j] + fl1_fx * tdy_zzz_xyz_0[j] + 0.5 * fl1_fx * ts_zzz_xyyz_0[j];

            tdz_yzzz_xyyz_0[j] = pa_y[j] * tdz_zzz_xyyz_0[j] + fl1_fx * tdz_zzz_xyz_0[j];

            tdx_yzzz_xyzz_0[j] = pa_y[j] * tdx_zzz_xyzz_0[j] + 0.5 * fl1_fx * tdx_zzz_xzz_0[j];

            tdy_yzzz_xyzz_0[j] = pa_y[j] * tdy_zzz_xyzz_0[j] + 0.5 * fl1_fx * tdy_zzz_xzz_0[j] + 0.5 * fl1_fx * ts_zzz_xyzz_0[j];

            tdz_yzzz_xyzz_0[j] = pa_y[j] * tdz_zzz_xyzz_0[j] + 0.5 * fl1_fx * tdz_zzz_xzz_0[j];

            tdx_yzzz_xzzz_0[j] = pa_y[j] * tdx_zzz_xzzz_0[j];

            tdy_yzzz_xzzz_0[j] = pa_y[j] * tdy_zzz_xzzz_0[j] + 0.5 * fl1_fx * ts_zzz_xzzz_0[j];

            tdz_yzzz_xzzz_0[j] = pa_y[j] * tdz_zzz_xzzz_0[j];

            tdx_yzzz_yyyy_0[j] = pa_y[j] * tdx_zzz_yyyy_0[j] + 2.0 * fl1_fx * tdx_zzz_yyy_0[j];

            tdy_yzzz_yyyy_0[j] = pa_y[j] * tdy_zzz_yyyy_0[j] + 2.0 * fl1_fx * tdy_zzz_yyy_0[j] + 0.5 * fl1_fx * ts_zzz_yyyy_0[j];

            tdz_yzzz_yyyy_0[j] = pa_y[j] * tdz_zzz_yyyy_0[j] + 2.0 * fl1_fx * tdz_zzz_yyy_0[j];

            tdx_yzzz_yyyz_0[j] = pa_y[j] * tdx_zzz_yyyz_0[j] + 1.5 * fl1_fx * tdx_zzz_yyz_0[j];

            tdy_yzzz_yyyz_0[j] = pa_y[j] * tdy_zzz_yyyz_0[j] + 1.5 * fl1_fx * tdy_zzz_yyz_0[j] + 0.5 * fl1_fx * ts_zzz_yyyz_0[j];

            tdz_yzzz_yyyz_0[j] = pa_y[j] * tdz_zzz_yyyz_0[j] + 1.5 * fl1_fx * tdz_zzz_yyz_0[j];

            tdx_yzzz_yyzz_0[j] = pa_y[j] * tdx_zzz_yyzz_0[j] + fl1_fx * tdx_zzz_yzz_0[j];

            tdy_yzzz_yyzz_0[j] = pa_y[j] * tdy_zzz_yyzz_0[j] + fl1_fx * tdy_zzz_yzz_0[j] + 0.5 * fl1_fx * ts_zzz_yyzz_0[j];

            tdz_yzzz_yyzz_0[j] = pa_y[j] * tdz_zzz_yyzz_0[j] + fl1_fx * tdz_zzz_yzz_0[j];

            tdx_yzzz_yzzz_0[j] = pa_y[j] * tdx_zzz_yzzz_0[j] + 0.5 * fl1_fx * tdx_zzz_zzz_0[j];

            tdy_yzzz_yzzz_0[j] = pa_y[j] * tdy_zzz_yzzz_0[j] + 0.5 * fl1_fx * tdy_zzz_zzz_0[j] + 0.5 * fl1_fx * ts_zzz_yzzz_0[j];

            tdz_yzzz_yzzz_0[j] = pa_y[j] * tdz_zzz_yzzz_0[j] + 0.5 * fl1_fx * tdz_zzz_zzz_0[j];
        }

        idx++;
    }
}

void
compElectricDipoleForGG_627_675(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const int32_t              nOSFactors,
                                const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_d_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_d_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tdx_zzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 90);

        auto tdy_zzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 90);

        auto tdz_zzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 90);

        auto tdx_zzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 91);

        auto tdy_zzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 91);

        auto tdz_zzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 91);

        auto tdx_zzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 92);

        auto tdy_zzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 92);

        auto tdz_zzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 92);

        auto tdx_zzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 93);

        auto tdy_zzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 93);

        auto tdz_zzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 93);

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

        auto ts_zzz_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 135);

        auto ts_zzz_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 136);

        auto ts_zzz_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 137);

        auto ts_zzz_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 138);

        auto ts_zzz_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 139);

        auto ts_zzz_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 140);

        auto ts_zzz_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 141);

        auto ts_zzz_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 142);

        auto ts_zzz_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 143);

        auto ts_zzz_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 144);

        auto ts_zzz_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 145);

        auto ts_zzz_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 146);

        auto ts_zzz_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 147);

        auto ts_zzz_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 148);

        auto ts_zzz_zzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 149);

        // set up pointers to integrals

        auto tdx_yzzz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 209);

        auto tdy_yzzz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 209);

        auto tdz_yzzz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 209);

        auto tdx_zzzz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 210);

        auto tdy_zzzz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 210);

        auto tdz_zzzz_xxxx_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 210);

        auto tdx_zzzz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 211);

        auto tdy_zzzz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 211);

        auto tdz_zzzz_xxxy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 211);

        auto tdx_zzzz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 212);

        auto tdy_zzzz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 212);

        auto tdz_zzzz_xxxz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 212);

        auto tdx_zzzz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 213);

        auto tdy_zzzz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 213);

        auto tdz_zzzz_xxyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 213);

        auto tdx_zzzz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 214);

        auto tdy_zzzz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 214);

        auto tdz_zzzz_xxyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 214);

        auto tdx_zzzz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 215);

        auto tdy_zzzz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 215);

        auto tdz_zzzz_xxzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 215);

        auto tdx_zzzz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 216);

        auto tdy_zzzz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 216);

        auto tdz_zzzz_xyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 216);

        auto tdx_zzzz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 217);

        auto tdy_zzzz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 217);

        auto tdz_zzzz_xyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 217);

        auto tdx_zzzz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 218);

        auto tdy_zzzz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 218);

        auto tdz_zzzz_xyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 218);

        auto tdx_zzzz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 219);

        auto tdy_zzzz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 219);

        auto tdz_zzzz_xzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 219);

        auto tdx_zzzz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 220);

        auto tdy_zzzz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 220);

        auto tdz_zzzz_yyyy_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 220);

        auto tdx_zzzz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 221);

        auto tdy_zzzz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 221);

        auto tdz_zzzz_yyyz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 221);

        auto tdx_zzzz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 222);

        auto tdy_zzzz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 222);

        auto tdz_zzzz_yyzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 222);

        auto tdx_zzzz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 223);

        auto tdy_zzzz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 223);

        auto tdz_zzzz_yzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 223);

        auto tdx_zzzz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * idx + 224);

        auto tdy_zzzz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 225 * bdim + 225 * idx + 224);

        auto tdz_zzzz_zzzz_0 = primBuffer.data(pidx_d_4_4_m0 + 450 * bdim + 225 * idx + 224);

        // Batch of Integrals (627,675)

        #pragma omp simd aligned(fx, pa_y, pa_z, tdx_yzzz_zzzz_0, tdx_zz_xxxx_0, tdx_zz_xxxy_0, \
                                     tdx_zz_xxxz_0, tdx_zz_xxyy_0, tdx_zz_xxyz_0, tdx_zz_xxzz_0, tdx_zz_xyyy_0, \
                                     tdx_zz_xyyz_0, tdx_zz_xyzz_0, tdx_zz_xzzz_0, tdx_zz_yyyy_0, tdx_zz_yyyz_0, \
                                     tdx_zz_yyzz_0, tdx_zz_yzzz_0, tdx_zz_zzzz_0, tdx_zzz_xxx_0, tdx_zzz_xxxx_0, \
                                     tdx_zzz_xxxy_0, tdx_zzz_xxxz_0, tdx_zzz_xxy_0, tdx_zzz_xxyy_0, tdx_zzz_xxyz_0, \
                                     tdx_zzz_xxz_0, tdx_zzz_xxzz_0, tdx_zzz_xyy_0, tdx_zzz_xyyy_0, tdx_zzz_xyyz_0, \
                                     tdx_zzz_xyz_0, tdx_zzz_xyzz_0, tdx_zzz_xzz_0, tdx_zzz_xzzz_0, tdx_zzz_yyy_0, \
                                     tdx_zzz_yyyy_0, tdx_zzz_yyyz_0, tdx_zzz_yyz_0, tdx_zzz_yyzz_0, tdx_zzz_yzz_0, \
                                     tdx_zzz_yzzz_0, tdx_zzz_zzz_0, tdx_zzz_zzzz_0, tdx_zzzz_xxxx_0, tdx_zzzz_xxxy_0, \
                                     tdx_zzzz_xxxz_0, tdx_zzzz_xxyy_0, tdx_zzzz_xxyz_0, tdx_zzzz_xxzz_0, tdx_zzzz_xyyy_0, \
                                     tdx_zzzz_xyyz_0, tdx_zzzz_xyzz_0, tdx_zzzz_xzzz_0, tdx_zzzz_yyyy_0, tdx_zzzz_yyyz_0, \
                                     tdx_zzzz_yyzz_0, tdx_zzzz_yzzz_0, tdx_zzzz_zzzz_0, tdy_yzzz_zzzz_0, tdy_zz_xxxx_0, \
                                     tdy_zz_xxxy_0, tdy_zz_xxxz_0, tdy_zz_xxyy_0, tdy_zz_xxyz_0, tdy_zz_xxzz_0, \
                                     tdy_zz_xyyy_0, tdy_zz_xyyz_0, tdy_zz_xyzz_0, tdy_zz_xzzz_0, tdy_zz_yyyy_0, \
                                     tdy_zz_yyyz_0, tdy_zz_yyzz_0, tdy_zz_yzzz_0, tdy_zz_zzzz_0, tdy_zzz_xxx_0, \
                                     tdy_zzz_xxxx_0, tdy_zzz_xxxy_0, tdy_zzz_xxxz_0, tdy_zzz_xxy_0, tdy_zzz_xxyy_0, \
                                     tdy_zzz_xxyz_0, tdy_zzz_xxz_0, tdy_zzz_xxzz_0, tdy_zzz_xyy_0, tdy_zzz_xyyy_0, \
                                     tdy_zzz_xyyz_0, tdy_zzz_xyz_0, tdy_zzz_xyzz_0, tdy_zzz_xzz_0, tdy_zzz_xzzz_0, \
                                     tdy_zzz_yyy_0, tdy_zzz_yyyy_0, tdy_zzz_yyyz_0, tdy_zzz_yyz_0, tdy_zzz_yyzz_0, \
                                     tdy_zzz_yzz_0, tdy_zzz_yzzz_0, tdy_zzz_zzz_0, tdy_zzz_zzzz_0, tdy_zzzz_xxxx_0, \
                                     tdy_zzzz_xxxy_0, tdy_zzzz_xxxz_0, tdy_zzzz_xxyy_0, tdy_zzzz_xxyz_0, tdy_zzzz_xxzz_0, \
                                     tdy_zzzz_xyyy_0, tdy_zzzz_xyyz_0, tdy_zzzz_xyzz_0, tdy_zzzz_xzzz_0, tdy_zzzz_yyyy_0, \
                                     tdy_zzzz_yyyz_0, tdy_zzzz_yyzz_0, tdy_zzzz_yzzz_0, tdy_zzzz_zzzz_0, tdz_yzzz_zzzz_0, \
                                     tdz_zz_xxxx_0, tdz_zz_xxxy_0, tdz_zz_xxxz_0, tdz_zz_xxyy_0, tdz_zz_xxyz_0, \
                                     tdz_zz_xxzz_0, tdz_zz_xyyy_0, tdz_zz_xyyz_0, tdz_zz_xyzz_0, tdz_zz_xzzz_0, \
                                     tdz_zz_yyyy_0, tdz_zz_yyyz_0, tdz_zz_yyzz_0, tdz_zz_yzzz_0, tdz_zz_zzzz_0, \
                                     tdz_zzz_xxx_0, tdz_zzz_xxxx_0, tdz_zzz_xxxy_0, tdz_zzz_xxxz_0, tdz_zzz_xxy_0, \
                                     tdz_zzz_xxyy_0, tdz_zzz_xxyz_0, tdz_zzz_xxz_0, tdz_zzz_xxzz_0, tdz_zzz_xyy_0, \
                                     tdz_zzz_xyyy_0, tdz_zzz_xyyz_0, tdz_zzz_xyz_0, tdz_zzz_xyzz_0, tdz_zzz_xzz_0, \
                                     tdz_zzz_xzzz_0, tdz_zzz_yyy_0, tdz_zzz_yyyy_0, tdz_zzz_yyyz_0, tdz_zzz_yyz_0, \
                                     tdz_zzz_yyzz_0, tdz_zzz_yzz_0, tdz_zzz_yzzz_0, tdz_zzz_zzz_0, tdz_zzz_zzzz_0, \
                                     tdz_zzzz_xxxx_0, tdz_zzzz_xxxy_0, tdz_zzzz_xxxz_0, tdz_zzzz_xxyy_0, tdz_zzzz_xxyz_0, \
                                     tdz_zzzz_xxzz_0, tdz_zzzz_xyyy_0, tdz_zzzz_xyyz_0, tdz_zzzz_xyzz_0, tdz_zzzz_xzzz_0, \
                                     tdz_zzzz_yyyy_0, tdz_zzzz_yyyz_0, tdz_zzzz_yyzz_0, tdz_zzzz_yzzz_0, tdz_zzzz_zzzz_0, \
                                     ts_zzz_xxxx_0, ts_zzz_xxxy_0, ts_zzz_xxxz_0, ts_zzz_xxyy_0, ts_zzz_xxyz_0, \
                                     ts_zzz_xxzz_0, ts_zzz_xyyy_0, ts_zzz_xyyz_0, ts_zzz_xyzz_0, ts_zzz_xzzz_0, \
                                     ts_zzz_yyyy_0, ts_zzz_yyyz_0, ts_zzz_yyzz_0, ts_zzz_yzzz_0, ts_zzz_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            tdx_yzzz_zzzz_0[j] = pa_y[j] * tdx_zzz_zzzz_0[j];

            tdy_yzzz_zzzz_0[j] = pa_y[j] * tdy_zzz_zzzz_0[j] + 0.5 * fl1_fx * ts_zzz_zzzz_0[j];

            tdz_yzzz_zzzz_0[j] = pa_y[j] * tdz_zzz_zzzz_0[j];

            tdx_zzzz_xxxx_0[j] = pa_z[j] * tdx_zzz_xxxx_0[j] + 1.5 * fl1_fx * tdx_zz_xxxx_0[j];

            tdy_zzzz_xxxx_0[j] = pa_z[j] * tdy_zzz_xxxx_0[j] + 1.5 * fl1_fx * tdy_zz_xxxx_0[j];

            tdz_zzzz_xxxx_0[j] = pa_z[j] * tdz_zzz_xxxx_0[j] + 1.5 * fl1_fx * tdz_zz_xxxx_0[j] + 0.5 * fl1_fx * ts_zzz_xxxx_0[j];

            tdx_zzzz_xxxy_0[j] = pa_z[j] * tdx_zzz_xxxy_0[j] + 1.5 * fl1_fx * tdx_zz_xxxy_0[j];

            tdy_zzzz_xxxy_0[j] = pa_z[j] * tdy_zzz_xxxy_0[j] + 1.5 * fl1_fx * tdy_zz_xxxy_0[j];

            tdz_zzzz_xxxy_0[j] = pa_z[j] * tdz_zzz_xxxy_0[j] + 1.5 * fl1_fx * tdz_zz_xxxy_0[j] + 0.5 * fl1_fx * ts_zzz_xxxy_0[j];

            tdx_zzzz_xxxz_0[j] = pa_z[j] * tdx_zzz_xxxz_0[j] + 1.5 * fl1_fx * tdx_zz_xxxz_0[j] + 0.5 * fl1_fx * tdx_zzz_xxx_0[j];

            tdy_zzzz_xxxz_0[j] = pa_z[j] * tdy_zzz_xxxz_0[j] + 1.5 * fl1_fx * tdy_zz_xxxz_0[j] + 0.5 * fl1_fx * tdy_zzz_xxx_0[j];

            tdz_zzzz_xxxz_0[j] =
                pa_z[j] * tdz_zzz_xxxz_0[j] + 1.5 * fl1_fx * tdz_zz_xxxz_0[j] + 0.5 * fl1_fx * tdz_zzz_xxx_0[j] + 0.5 * fl1_fx * ts_zzz_xxxz_0[j];

            tdx_zzzz_xxyy_0[j] = pa_z[j] * tdx_zzz_xxyy_0[j] + 1.5 * fl1_fx * tdx_zz_xxyy_0[j];

            tdy_zzzz_xxyy_0[j] = pa_z[j] * tdy_zzz_xxyy_0[j] + 1.5 * fl1_fx * tdy_zz_xxyy_0[j];

            tdz_zzzz_xxyy_0[j] = pa_z[j] * tdz_zzz_xxyy_0[j] + 1.5 * fl1_fx * tdz_zz_xxyy_0[j] + 0.5 * fl1_fx * ts_zzz_xxyy_0[j];

            tdx_zzzz_xxyz_0[j] = pa_z[j] * tdx_zzz_xxyz_0[j] + 1.5 * fl1_fx * tdx_zz_xxyz_0[j] + 0.5 * fl1_fx * tdx_zzz_xxy_0[j];

            tdy_zzzz_xxyz_0[j] = pa_z[j] * tdy_zzz_xxyz_0[j] + 1.5 * fl1_fx * tdy_zz_xxyz_0[j] + 0.5 * fl1_fx * tdy_zzz_xxy_0[j];

            tdz_zzzz_xxyz_0[j] =
                pa_z[j] * tdz_zzz_xxyz_0[j] + 1.5 * fl1_fx * tdz_zz_xxyz_0[j] + 0.5 * fl1_fx * tdz_zzz_xxy_0[j] + 0.5 * fl1_fx * ts_zzz_xxyz_0[j];

            tdx_zzzz_xxzz_0[j] = pa_z[j] * tdx_zzz_xxzz_0[j] + 1.5 * fl1_fx * tdx_zz_xxzz_0[j] + fl1_fx * tdx_zzz_xxz_0[j];

            tdy_zzzz_xxzz_0[j] = pa_z[j] * tdy_zzz_xxzz_0[j] + 1.5 * fl1_fx * tdy_zz_xxzz_0[j] + fl1_fx * tdy_zzz_xxz_0[j];

            tdz_zzzz_xxzz_0[j] =
                pa_z[j] * tdz_zzz_xxzz_0[j] + 1.5 * fl1_fx * tdz_zz_xxzz_0[j] + fl1_fx * tdz_zzz_xxz_0[j] + 0.5 * fl1_fx * ts_zzz_xxzz_0[j];

            tdx_zzzz_xyyy_0[j] = pa_z[j] * tdx_zzz_xyyy_0[j] + 1.5 * fl1_fx * tdx_zz_xyyy_0[j];

            tdy_zzzz_xyyy_0[j] = pa_z[j] * tdy_zzz_xyyy_0[j] + 1.5 * fl1_fx * tdy_zz_xyyy_0[j];

            tdz_zzzz_xyyy_0[j] = pa_z[j] * tdz_zzz_xyyy_0[j] + 1.5 * fl1_fx * tdz_zz_xyyy_0[j] + 0.5 * fl1_fx * ts_zzz_xyyy_0[j];

            tdx_zzzz_xyyz_0[j] = pa_z[j] * tdx_zzz_xyyz_0[j] + 1.5 * fl1_fx * tdx_zz_xyyz_0[j] + 0.5 * fl1_fx * tdx_zzz_xyy_0[j];

            tdy_zzzz_xyyz_0[j] = pa_z[j] * tdy_zzz_xyyz_0[j] + 1.5 * fl1_fx * tdy_zz_xyyz_0[j] + 0.5 * fl1_fx * tdy_zzz_xyy_0[j];

            tdz_zzzz_xyyz_0[j] =
                pa_z[j] * tdz_zzz_xyyz_0[j] + 1.5 * fl1_fx * tdz_zz_xyyz_0[j] + 0.5 * fl1_fx * tdz_zzz_xyy_0[j] + 0.5 * fl1_fx * ts_zzz_xyyz_0[j];

            tdx_zzzz_xyzz_0[j] = pa_z[j] * tdx_zzz_xyzz_0[j] + 1.5 * fl1_fx * tdx_zz_xyzz_0[j] + fl1_fx * tdx_zzz_xyz_0[j];

            tdy_zzzz_xyzz_0[j] = pa_z[j] * tdy_zzz_xyzz_0[j] + 1.5 * fl1_fx * tdy_zz_xyzz_0[j] + fl1_fx * tdy_zzz_xyz_0[j];

            tdz_zzzz_xyzz_0[j] =
                pa_z[j] * tdz_zzz_xyzz_0[j] + 1.5 * fl1_fx * tdz_zz_xyzz_0[j] + fl1_fx * tdz_zzz_xyz_0[j] + 0.5 * fl1_fx * ts_zzz_xyzz_0[j];

            tdx_zzzz_xzzz_0[j] = pa_z[j] * tdx_zzz_xzzz_0[j] + 1.5 * fl1_fx * tdx_zz_xzzz_0[j] + 1.5 * fl1_fx * tdx_zzz_xzz_0[j];

            tdy_zzzz_xzzz_0[j] = pa_z[j] * tdy_zzz_xzzz_0[j] + 1.5 * fl1_fx * tdy_zz_xzzz_0[j] + 1.5 * fl1_fx * tdy_zzz_xzz_0[j];

            tdz_zzzz_xzzz_0[j] =
                pa_z[j] * tdz_zzz_xzzz_0[j] + 1.5 * fl1_fx * tdz_zz_xzzz_0[j] + 1.5 * fl1_fx * tdz_zzz_xzz_0[j] + 0.5 * fl1_fx * ts_zzz_xzzz_0[j];

            tdx_zzzz_yyyy_0[j] = pa_z[j] * tdx_zzz_yyyy_0[j] + 1.5 * fl1_fx * tdx_zz_yyyy_0[j];

            tdy_zzzz_yyyy_0[j] = pa_z[j] * tdy_zzz_yyyy_0[j] + 1.5 * fl1_fx * tdy_zz_yyyy_0[j];

            tdz_zzzz_yyyy_0[j] = pa_z[j] * tdz_zzz_yyyy_0[j] + 1.5 * fl1_fx * tdz_zz_yyyy_0[j] + 0.5 * fl1_fx * ts_zzz_yyyy_0[j];

            tdx_zzzz_yyyz_0[j] = pa_z[j] * tdx_zzz_yyyz_0[j] + 1.5 * fl1_fx * tdx_zz_yyyz_0[j] + 0.5 * fl1_fx * tdx_zzz_yyy_0[j];

            tdy_zzzz_yyyz_0[j] = pa_z[j] * tdy_zzz_yyyz_0[j] + 1.5 * fl1_fx * tdy_zz_yyyz_0[j] + 0.5 * fl1_fx * tdy_zzz_yyy_0[j];

            tdz_zzzz_yyyz_0[j] =
                pa_z[j] * tdz_zzz_yyyz_0[j] + 1.5 * fl1_fx * tdz_zz_yyyz_0[j] + 0.5 * fl1_fx * tdz_zzz_yyy_0[j] + 0.5 * fl1_fx * ts_zzz_yyyz_0[j];

            tdx_zzzz_yyzz_0[j] = pa_z[j] * tdx_zzz_yyzz_0[j] + 1.5 * fl1_fx * tdx_zz_yyzz_0[j] + fl1_fx * tdx_zzz_yyz_0[j];

            tdy_zzzz_yyzz_0[j] = pa_z[j] * tdy_zzz_yyzz_0[j] + 1.5 * fl1_fx * tdy_zz_yyzz_0[j] + fl1_fx * tdy_zzz_yyz_0[j];

            tdz_zzzz_yyzz_0[j] =
                pa_z[j] * tdz_zzz_yyzz_0[j] + 1.5 * fl1_fx * tdz_zz_yyzz_0[j] + fl1_fx * tdz_zzz_yyz_0[j] + 0.5 * fl1_fx * ts_zzz_yyzz_0[j];

            tdx_zzzz_yzzz_0[j] = pa_z[j] * tdx_zzz_yzzz_0[j] + 1.5 * fl1_fx * tdx_zz_yzzz_0[j] + 1.5 * fl1_fx * tdx_zzz_yzz_0[j];

            tdy_zzzz_yzzz_0[j] = pa_z[j] * tdy_zzz_yzzz_0[j] + 1.5 * fl1_fx * tdy_zz_yzzz_0[j] + 1.5 * fl1_fx * tdy_zzz_yzz_0[j];

            tdz_zzzz_yzzz_0[j] =
                pa_z[j] * tdz_zzz_yzzz_0[j] + 1.5 * fl1_fx * tdz_zz_yzzz_0[j] + 1.5 * fl1_fx * tdz_zzz_yzz_0[j] + 0.5 * fl1_fx * ts_zzz_yzzz_0[j];

            tdx_zzzz_zzzz_0[j] = pa_z[j] * tdx_zzz_zzzz_0[j] + 1.5 * fl1_fx * tdx_zz_zzzz_0[j] + 2.0 * fl1_fx * tdx_zzz_zzz_0[j];

            tdy_zzzz_zzzz_0[j] = pa_z[j] * tdy_zzz_zzzz_0[j] + 1.5 * fl1_fx * tdy_zz_zzzz_0[j] + 2.0 * fl1_fx * tdy_zzz_zzz_0[j];

            tdz_zzzz_zzzz_0[j] =
                pa_z[j] * tdz_zzz_zzzz_0[j] + 1.5 * fl1_fx * tdz_zz_zzzz_0[j] + 2.0 * fl1_fx * tdz_zzz_zzz_0[j] + 0.5 * fl1_fx * ts_zzz_zzzz_0[j];
        }

        idx++;
    }
}

}  // namespace ediprecfunc
