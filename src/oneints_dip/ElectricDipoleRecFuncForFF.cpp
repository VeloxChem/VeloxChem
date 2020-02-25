//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ElectricDipoleRecFuncForFF.hpp"

namespace ediprecfunc {  // ediprecfunc namespace

void
compElectricDipoleForFF(CMemBlock2D<double>&       primBuffer,
                        const CRecursionMap&       recursionMap,
                        const CMemBlock2D<double>& osFactors,
                        const int32_t              nOSFactors,
                        const CMemBlock2D<double>& paDistances,
                        const CGtoBlock&           braGtoBlock,
                        const CGtoBlock&           ketGtoBlock,
                        const int32_t              iContrGto)
{
    ediprecfunc::compElectricDipoleForFF_0_50(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForFF_50_100(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForFF_100_150(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForFF_150_200(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForFF_200_250(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForFF_250_300(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compElectricDipoleForFF_0_50(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_d_3_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

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

        auto tdx_x_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx);

        auto tdy_x_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx);

        auto tdz_x_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx);

        auto tdx_x_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 1);

        auto tdy_x_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 1);

        auto tdz_x_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 1);

        auto tdx_x_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 2);

        auto tdy_x_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 2);

        auto tdz_x_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 2);

        auto tdx_x_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 3);

        auto tdy_x_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 3);

        auto tdz_x_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 3);

        auto tdx_x_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 4);

        auto tdy_x_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 4);

        auto tdz_x_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 4);

        auto tdx_x_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 5);

        auto tdy_x_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 5);

        auto tdz_x_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 5);

        auto tdx_x_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 6);

        auto tdy_x_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 6);

        auto tdz_x_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 6);

        auto tdx_x_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 7);

        auto tdy_x_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 7);

        auto tdz_x_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 7);

        auto tdx_x_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 8);

        auto tdy_x_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 8);

        auto tdz_x_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 8);

        auto tdx_x_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 9);

        auto tdy_x_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 9);

        auto tdz_x_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 9);

        auto tdx_y_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 10);

        auto tdy_y_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 10);

        auto tdz_y_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 10);

        auto tdx_y_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 11);

        auto tdy_y_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 11);

        auto tdz_y_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 11);

        auto tdx_y_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 12);

        auto tdy_y_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 12);

        auto tdz_y_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 12);

        auto tdx_y_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 13);

        auto tdy_y_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 13);

        auto tdz_y_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 13);

        auto tdx_y_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 14);

        auto tdy_y_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 14);

        auto tdz_y_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 14);

        auto tdx_y_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 15);

        auto tdy_y_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 15);

        auto tdz_y_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 15);

        auto tdx_y_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 16);

        auto tdy_y_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 16);

        auto tdx_xx_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx);

        auto tdy_xx_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx);

        auto tdz_xx_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx);

        auto tdx_xx_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 1);

        auto tdy_xx_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 1);

        auto tdz_xx_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 1);

        auto tdx_xx_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 2);

        auto tdy_xx_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 2);

        auto tdz_xx_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 2);

        auto tdx_xx_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 3);

        auto tdy_xx_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 3);

        auto tdz_xx_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 3);

        auto tdx_xx_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 4);

        auto tdy_xx_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 4);

        auto tdz_xx_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 4);

        auto tdx_xx_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 5);

        auto tdy_xx_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 5);

        auto tdz_xx_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 5);

        auto tdx_xy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 6);

        auto tdy_xy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 6);

        auto tdz_xy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 6);

        auto tdx_xy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 7);

        auto tdy_xy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 7);

        auto tdz_xy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 7);

        auto tdx_xy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 8);

        auto tdy_xy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 8);

        auto tdz_xy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 8);

        auto tdx_xy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 9);

        auto tdy_xy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 9);

        auto tdz_xy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 9);

        auto tdx_xy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 10);

        auto tdy_xy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 10);

        auto tdz_xy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 10);

        auto tdx_xy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 11);

        auto tdy_xy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 11);

        auto tdz_xy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 11);

        auto ts_xx_xxx_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx);

        auto ts_xx_xxy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 1);

        auto ts_xx_xxz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 2);

        auto ts_xx_xyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 3);

        auto ts_xx_xyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 4);

        auto ts_xx_xzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 5);

        auto ts_xx_yyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 6);

        auto ts_xx_yyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 7);

        auto ts_xx_yzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 8);

        auto ts_xx_zzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 9);

        auto ts_xy_xxx_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 10);

        auto ts_xy_xxy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 11);

        auto ts_xy_xxz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 12);

        auto ts_xy_xyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 13);

        auto ts_xy_xyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 14);

        auto ts_xy_xzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 15);

        auto ts_xy_yyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 16);

        // set up pointers to integrals

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

        // Batch of Integrals (0,50)

        #pragma omp simd aligned(fx, pa_x, tdx_x_xxx_0, tdx_x_xxy_0, tdx_x_xxz_0, tdx_x_xyy_0, \
                                     tdx_x_xyz_0, tdx_x_xzz_0, tdx_x_yyy_0, tdx_x_yyz_0, tdx_x_yzz_0, tdx_x_zzz_0, \
                                     tdx_xx_xx_0, tdx_xx_xxx_0, tdx_xx_xxy_0, tdx_xx_xxz_0, tdx_xx_xy_0, tdx_xx_xyy_0, \
                                     tdx_xx_xyz_0, tdx_xx_xz_0, tdx_xx_xzz_0, tdx_xx_yy_0, tdx_xx_yyy_0, tdx_xx_yyz_0, \
                                     tdx_xx_yz_0, tdx_xx_yzz_0, tdx_xx_zz_0, tdx_xx_zzz_0, tdx_xxx_xxx_0, \
                                     tdx_xxx_xxy_0, tdx_xxx_xxz_0, tdx_xxx_xyy_0, tdx_xxx_xyz_0, tdx_xxx_xzz_0, \
                                     tdx_xxx_yyy_0, tdx_xxx_yyz_0, tdx_xxx_yzz_0, tdx_xxx_zzz_0, tdx_xxy_xxx_0, \
                                     tdx_xxy_xxy_0, tdx_xxy_xxz_0, tdx_xxy_xyy_0, tdx_xxy_xyz_0, tdx_xxy_xzz_0, \
                                     tdx_xxy_yyy_0, tdx_xy_xx_0, tdx_xy_xxx_0, tdx_xy_xxy_0, tdx_xy_xxz_0, tdx_xy_xy_0, \
                                     tdx_xy_xyy_0, tdx_xy_xyz_0, tdx_xy_xz_0, tdx_xy_xzz_0, tdx_xy_yy_0, tdx_xy_yyy_0, \
                                     tdx_xy_yz_0, tdx_xy_zz_0, tdx_y_xxx_0, tdx_y_xxy_0, tdx_y_xxz_0, tdx_y_xyy_0, \
                                     tdx_y_xyz_0, tdx_y_xzz_0, tdx_y_yyy_0, tdy_x_xxx_0, tdy_x_xxy_0, tdy_x_xxz_0, \
                                     tdy_x_xyy_0, tdy_x_xyz_0, tdy_x_xzz_0, tdy_x_yyy_0, tdy_x_yyz_0, tdy_x_yzz_0, \
                                     tdy_x_zzz_0, tdy_xx_xx_0, tdy_xx_xxx_0, tdy_xx_xxy_0, tdy_xx_xxz_0, tdy_xx_xy_0, \
                                     tdy_xx_xyy_0, tdy_xx_xyz_0, tdy_xx_xz_0, tdy_xx_xzz_0, tdy_xx_yy_0, tdy_xx_yyy_0, \
                                     tdy_xx_yyz_0, tdy_xx_yz_0, tdy_xx_yzz_0, tdy_xx_zz_0, tdy_xx_zzz_0, tdy_xxx_xxx_0, \
                                     tdy_xxx_xxy_0, tdy_xxx_xxz_0, tdy_xxx_xyy_0, tdy_xxx_xyz_0, tdy_xxx_xzz_0, \
                                     tdy_xxx_yyy_0, tdy_xxx_yyz_0, tdy_xxx_yzz_0, tdy_xxx_zzz_0, tdy_xxy_xxx_0, \
                                     tdy_xxy_xxy_0, tdy_xxy_xxz_0, tdy_xxy_xyy_0, tdy_xxy_xyz_0, tdy_xxy_xzz_0, \
                                     tdy_xxy_yyy_0, tdy_xy_xx_0, tdy_xy_xxx_0, tdy_xy_xxy_0, tdy_xy_xxz_0, tdy_xy_xy_0, \
                                     tdy_xy_xyy_0, tdy_xy_xyz_0, tdy_xy_xz_0, tdy_xy_xzz_0, tdy_xy_yy_0, tdy_xy_yyy_0, \
                                     tdy_xy_yz_0, tdy_xy_zz_0, tdy_y_xxx_0, tdy_y_xxy_0, tdy_y_xxz_0, tdy_y_xyy_0, \
                                     tdy_y_xyz_0, tdy_y_xzz_0, tdy_y_yyy_0, tdz_x_xxx_0, tdz_x_xxy_0, tdz_x_xxz_0, \
                                     tdz_x_xyy_0, tdz_x_xyz_0, tdz_x_xzz_0, tdz_x_yyy_0, tdz_x_yyz_0, tdz_x_yzz_0, \
                                     tdz_x_zzz_0, tdz_xx_xx_0, tdz_xx_xxx_0, tdz_xx_xxy_0, tdz_xx_xxz_0, tdz_xx_xy_0, \
                                     tdz_xx_xyy_0, tdz_xx_xyz_0, tdz_xx_xz_0, tdz_xx_xzz_0, tdz_xx_yy_0, tdz_xx_yyy_0, \
                                     tdz_xx_yyz_0, tdz_xx_yz_0, tdz_xx_yzz_0, tdz_xx_zz_0, tdz_xx_zzz_0, tdz_xxx_xxx_0, \
                                     tdz_xxx_xxy_0, tdz_xxx_xxz_0, tdz_xxx_xyy_0, tdz_xxx_xyz_0, tdz_xxx_xzz_0, \
                                     tdz_xxx_yyy_0, tdz_xxx_yyz_0, tdz_xxx_yzz_0, tdz_xxx_zzz_0, tdz_xxy_xxx_0, \
                                     tdz_xxy_xxy_0, tdz_xxy_xxz_0, tdz_xxy_xyy_0, tdz_xxy_xyz_0, tdz_xxy_xzz_0, \
                                     tdz_xy_xx_0, tdz_xy_xxx_0, tdz_xy_xxy_0, tdz_xy_xxz_0, tdz_xy_xy_0, tdz_xy_xyy_0, \
                                     tdz_xy_xyz_0, tdz_xy_xz_0, tdz_xy_xzz_0, tdz_xy_yy_0, tdz_xy_yz_0, tdz_xy_zz_0, \
                                     tdz_y_xxx_0, tdz_y_xxy_0, tdz_y_xxz_0, tdz_y_xyy_0, tdz_y_xyz_0, tdz_y_xzz_0, \
                                     ts_xx_xxx_0, ts_xx_xxy_0, ts_xx_xxz_0, ts_xx_xyy_0, ts_xx_xyz_0, ts_xx_xzz_0, \
                                     ts_xx_yyy_0, ts_xx_yyz_0, ts_xx_yzz_0, ts_xx_zzz_0, ts_xy_xxx_0, ts_xy_xxy_0, \
                                     ts_xy_xxz_0, ts_xy_xyy_0, ts_xy_xyz_0, ts_xy_xzz_0, ts_xy_yyy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            tdx_xxx_xxx_0[j] = pa_x[j] * tdx_xx_xxx_0[j] + fl1_fx * tdx_x_xxx_0[j] + 1.5 * fl1_fx * tdx_xx_xx_0[j] + 0.5 * fl1_fx * ts_xx_xxx_0[j];

            tdy_xxx_xxx_0[j] = pa_x[j] * tdy_xx_xxx_0[j] + fl1_fx * tdy_x_xxx_0[j] + 1.5 * fl1_fx * tdy_xx_xx_0[j];

            tdz_xxx_xxx_0[j] = pa_x[j] * tdz_xx_xxx_0[j] + fl1_fx * tdz_x_xxx_0[j] + 1.5 * fl1_fx * tdz_xx_xx_0[j];

            tdx_xxx_xxy_0[j] = pa_x[j] * tdx_xx_xxy_0[j] + fl1_fx * tdx_x_xxy_0[j] + fl1_fx * tdx_xx_xy_0[j] + 0.5 * fl1_fx * ts_xx_xxy_0[j];

            tdy_xxx_xxy_0[j] = pa_x[j] * tdy_xx_xxy_0[j] + fl1_fx * tdy_x_xxy_0[j] + fl1_fx * tdy_xx_xy_0[j];

            tdz_xxx_xxy_0[j] = pa_x[j] * tdz_xx_xxy_0[j] + fl1_fx * tdz_x_xxy_0[j] + fl1_fx * tdz_xx_xy_0[j];

            tdx_xxx_xxz_0[j] = pa_x[j] * tdx_xx_xxz_0[j] + fl1_fx * tdx_x_xxz_0[j] + fl1_fx * tdx_xx_xz_0[j] + 0.5 * fl1_fx * ts_xx_xxz_0[j];

            tdy_xxx_xxz_0[j] = pa_x[j] * tdy_xx_xxz_0[j] + fl1_fx * tdy_x_xxz_0[j] + fl1_fx * tdy_xx_xz_0[j];

            tdz_xxx_xxz_0[j] = pa_x[j] * tdz_xx_xxz_0[j] + fl1_fx * tdz_x_xxz_0[j] + fl1_fx * tdz_xx_xz_0[j];

            tdx_xxx_xyy_0[j] = pa_x[j] * tdx_xx_xyy_0[j] + fl1_fx * tdx_x_xyy_0[j] + 0.5 * fl1_fx * tdx_xx_yy_0[j] + 0.5 * fl1_fx * ts_xx_xyy_0[j];

            tdy_xxx_xyy_0[j] = pa_x[j] * tdy_xx_xyy_0[j] + fl1_fx * tdy_x_xyy_0[j] + 0.5 * fl1_fx * tdy_xx_yy_0[j];

            tdz_xxx_xyy_0[j] = pa_x[j] * tdz_xx_xyy_0[j] + fl1_fx * tdz_x_xyy_0[j] + 0.5 * fl1_fx * tdz_xx_yy_0[j];

            tdx_xxx_xyz_0[j] = pa_x[j] * tdx_xx_xyz_0[j] + fl1_fx * tdx_x_xyz_0[j] + 0.5 * fl1_fx * tdx_xx_yz_0[j] + 0.5 * fl1_fx * ts_xx_xyz_0[j];

            tdy_xxx_xyz_0[j] = pa_x[j] * tdy_xx_xyz_0[j] + fl1_fx * tdy_x_xyz_0[j] + 0.5 * fl1_fx * tdy_xx_yz_0[j];

            tdz_xxx_xyz_0[j] = pa_x[j] * tdz_xx_xyz_0[j] + fl1_fx * tdz_x_xyz_0[j] + 0.5 * fl1_fx * tdz_xx_yz_0[j];

            tdx_xxx_xzz_0[j] = pa_x[j] * tdx_xx_xzz_0[j] + fl1_fx * tdx_x_xzz_0[j] + 0.5 * fl1_fx * tdx_xx_zz_0[j] + 0.5 * fl1_fx * ts_xx_xzz_0[j];

            tdy_xxx_xzz_0[j] = pa_x[j] * tdy_xx_xzz_0[j] + fl1_fx * tdy_x_xzz_0[j] + 0.5 * fl1_fx * tdy_xx_zz_0[j];

            tdz_xxx_xzz_0[j] = pa_x[j] * tdz_xx_xzz_0[j] + fl1_fx * tdz_x_xzz_0[j] + 0.5 * fl1_fx * tdz_xx_zz_0[j];

            tdx_xxx_yyy_0[j] = pa_x[j] * tdx_xx_yyy_0[j] + fl1_fx * tdx_x_yyy_0[j] + 0.5 * fl1_fx * ts_xx_yyy_0[j];

            tdy_xxx_yyy_0[j] = pa_x[j] * tdy_xx_yyy_0[j] + fl1_fx * tdy_x_yyy_0[j];

            tdz_xxx_yyy_0[j] = pa_x[j] * tdz_xx_yyy_0[j] + fl1_fx * tdz_x_yyy_0[j];

            tdx_xxx_yyz_0[j] = pa_x[j] * tdx_xx_yyz_0[j] + fl1_fx * tdx_x_yyz_0[j] + 0.5 * fl1_fx * ts_xx_yyz_0[j];

            tdy_xxx_yyz_0[j] = pa_x[j] * tdy_xx_yyz_0[j] + fl1_fx * tdy_x_yyz_0[j];

            tdz_xxx_yyz_0[j] = pa_x[j] * tdz_xx_yyz_0[j] + fl1_fx * tdz_x_yyz_0[j];

            tdx_xxx_yzz_0[j] = pa_x[j] * tdx_xx_yzz_0[j] + fl1_fx * tdx_x_yzz_0[j] + 0.5 * fl1_fx * ts_xx_yzz_0[j];

            tdy_xxx_yzz_0[j] = pa_x[j] * tdy_xx_yzz_0[j] + fl1_fx * tdy_x_yzz_0[j];

            tdz_xxx_yzz_0[j] = pa_x[j] * tdz_xx_yzz_0[j] + fl1_fx * tdz_x_yzz_0[j];

            tdx_xxx_zzz_0[j] = pa_x[j] * tdx_xx_zzz_0[j] + fl1_fx * tdx_x_zzz_0[j] + 0.5 * fl1_fx * ts_xx_zzz_0[j];

            tdy_xxx_zzz_0[j] = pa_x[j] * tdy_xx_zzz_0[j] + fl1_fx * tdy_x_zzz_0[j];

            tdz_xxx_zzz_0[j] = pa_x[j] * tdz_xx_zzz_0[j] + fl1_fx * tdz_x_zzz_0[j];

            tdx_xxy_xxx_0[j] =
                pa_x[j] * tdx_xy_xxx_0[j] + 0.5 * fl1_fx * tdx_y_xxx_0[j] + 1.5 * fl1_fx * tdx_xy_xx_0[j] + 0.5 * fl1_fx * ts_xy_xxx_0[j];

            tdy_xxy_xxx_0[j] = pa_x[j] * tdy_xy_xxx_0[j] + 0.5 * fl1_fx * tdy_y_xxx_0[j] + 1.5 * fl1_fx * tdy_xy_xx_0[j];

            tdz_xxy_xxx_0[j] = pa_x[j] * tdz_xy_xxx_0[j] + 0.5 * fl1_fx * tdz_y_xxx_0[j] + 1.5 * fl1_fx * tdz_xy_xx_0[j];

            tdx_xxy_xxy_0[j] = pa_x[j] * tdx_xy_xxy_0[j] + 0.5 * fl1_fx * tdx_y_xxy_0[j] + fl1_fx * tdx_xy_xy_0[j] + 0.5 * fl1_fx * ts_xy_xxy_0[j];

            tdy_xxy_xxy_0[j] = pa_x[j] * tdy_xy_xxy_0[j] + 0.5 * fl1_fx * tdy_y_xxy_0[j] + fl1_fx * tdy_xy_xy_0[j];

            tdz_xxy_xxy_0[j] = pa_x[j] * tdz_xy_xxy_0[j] + 0.5 * fl1_fx * tdz_y_xxy_0[j] + fl1_fx * tdz_xy_xy_0[j];

            tdx_xxy_xxz_0[j] = pa_x[j] * tdx_xy_xxz_0[j] + 0.5 * fl1_fx * tdx_y_xxz_0[j] + fl1_fx * tdx_xy_xz_0[j] + 0.5 * fl1_fx * ts_xy_xxz_0[j];

            tdy_xxy_xxz_0[j] = pa_x[j] * tdy_xy_xxz_0[j] + 0.5 * fl1_fx * tdy_y_xxz_0[j] + fl1_fx * tdy_xy_xz_0[j];

            tdz_xxy_xxz_0[j] = pa_x[j] * tdz_xy_xxz_0[j] + 0.5 * fl1_fx * tdz_y_xxz_0[j] + fl1_fx * tdz_xy_xz_0[j];

            tdx_xxy_xyy_0[j] =
                pa_x[j] * tdx_xy_xyy_0[j] + 0.5 * fl1_fx * tdx_y_xyy_0[j] + 0.5 * fl1_fx * tdx_xy_yy_0[j] + 0.5 * fl1_fx * ts_xy_xyy_0[j];

            tdy_xxy_xyy_0[j] = pa_x[j] * tdy_xy_xyy_0[j] + 0.5 * fl1_fx * tdy_y_xyy_0[j] + 0.5 * fl1_fx * tdy_xy_yy_0[j];

            tdz_xxy_xyy_0[j] = pa_x[j] * tdz_xy_xyy_0[j] + 0.5 * fl1_fx * tdz_y_xyy_0[j] + 0.5 * fl1_fx * tdz_xy_yy_0[j];

            tdx_xxy_xyz_0[j] =
                pa_x[j] * tdx_xy_xyz_0[j] + 0.5 * fl1_fx * tdx_y_xyz_0[j] + 0.5 * fl1_fx * tdx_xy_yz_0[j] + 0.5 * fl1_fx * ts_xy_xyz_0[j];

            tdy_xxy_xyz_0[j] = pa_x[j] * tdy_xy_xyz_0[j] + 0.5 * fl1_fx * tdy_y_xyz_0[j] + 0.5 * fl1_fx * tdy_xy_yz_0[j];

            tdz_xxy_xyz_0[j] = pa_x[j] * tdz_xy_xyz_0[j] + 0.5 * fl1_fx * tdz_y_xyz_0[j] + 0.5 * fl1_fx * tdz_xy_yz_0[j];

            tdx_xxy_xzz_0[j] =
                pa_x[j] * tdx_xy_xzz_0[j] + 0.5 * fl1_fx * tdx_y_xzz_0[j] + 0.5 * fl1_fx * tdx_xy_zz_0[j] + 0.5 * fl1_fx * ts_xy_xzz_0[j];

            tdy_xxy_xzz_0[j] = pa_x[j] * tdy_xy_xzz_0[j] + 0.5 * fl1_fx * tdy_y_xzz_0[j] + 0.5 * fl1_fx * tdy_xy_zz_0[j];

            tdz_xxy_xzz_0[j] = pa_x[j] * tdz_xy_xzz_0[j] + 0.5 * fl1_fx * tdz_y_xzz_0[j] + 0.5 * fl1_fx * tdz_xy_zz_0[j];

            tdx_xxy_yyy_0[j] = pa_x[j] * tdx_xy_yyy_0[j] + 0.5 * fl1_fx * tdx_y_yyy_0[j] + 0.5 * fl1_fx * ts_xy_yyy_0[j];

            tdy_xxy_yyy_0[j] = pa_x[j] * tdy_xy_yyy_0[j] + 0.5 * fl1_fx * tdy_y_yyy_0[j];
        }

        idx++;
    }
}

void
compElectricDipoleForFF_50_100(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_d_3_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

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

        auto tdz_y_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 16);

        auto tdx_y_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 17);

        auto tdy_y_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 17);

        auto tdz_y_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 17);

        auto tdx_y_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 18);

        auto tdy_y_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 18);

        auto tdz_y_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 18);

        auto tdx_y_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 19);

        auto tdy_y_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 19);

        auto tdz_y_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 19);

        auto tdx_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 20);

        auto tdy_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 20);

        auto tdz_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 20);

        auto tdx_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 21);

        auto tdy_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 21);

        auto tdz_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 21);

        auto tdx_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 22);

        auto tdy_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 22);

        auto tdz_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 22);

        auto tdx_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 23);

        auto tdy_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 23);

        auto tdz_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 23);

        auto tdx_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 24);

        auto tdy_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 24);

        auto tdz_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 24);

        auto tdx_z_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 25);

        auto tdy_z_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 25);

        auto tdz_z_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 25);

        auto tdx_z_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 26);

        auto tdy_z_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 26);

        auto tdz_z_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 26);

        auto tdx_z_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 27);

        auto tdy_z_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 27);

        auto tdz_z_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 27);

        auto tdx_z_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 28);

        auto tdy_z_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 28);

        auto tdz_z_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 28);

        auto tdx_z_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 29);

        auto tdy_z_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 29);

        auto tdz_z_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 29);

        auto tdx_xz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 12);

        auto tdy_xz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 12);

        auto tdz_xz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 12);

        auto tdx_xz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 13);

        auto tdy_xz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 13);

        auto tdz_xz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 13);

        auto tdx_xz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 14);

        auto tdy_xz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 14);

        auto tdz_xz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 14);

        auto tdx_xz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 15);

        auto tdy_xz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 15);

        auto tdz_xz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 15);

        auto tdx_xz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 16);

        auto tdy_xz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 16);

        auto tdz_xz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 16);

        auto tdx_xz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 17);

        auto tdy_xz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 17);

        auto tdz_xz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 17);

        auto tdx_yy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 18);

        auto tdy_yy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 18);

        auto tdz_yy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 18);

        auto tdx_yy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 19);

        auto tdy_yy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 19);

        auto tdz_yy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 19);

        auto tdx_yy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 20);

        auto tdy_yy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 20);

        auto tdz_yy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 20);

        auto tdx_yy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 21);

        auto ts_xy_yyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 17);

        auto ts_xy_yzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 18);

        auto ts_xy_zzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 19);

        auto ts_xz_xxx_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 20);

        auto ts_xz_xxy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 21);

        auto ts_xz_xxz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 22);

        auto ts_xz_xyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 23);

        auto ts_xz_xyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 24);

        auto ts_xz_xzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 25);

        auto ts_xz_yyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 26);

        auto ts_xz_yyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 27);

        auto ts_xz_yzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 28);

        auto ts_xz_zzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 29);

        auto ts_yy_xxx_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 30);

        auto ts_yy_xxy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 31);

        auto ts_yy_xxz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 32);

        auto ts_yy_xyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 33);

        // set up pointers to integrals

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

        // Batch of Integrals (50,100)

        #pragma omp simd aligned(fx, pa_x, tdx_xxy_yyz_0, tdx_xxy_yzz_0, tdx_xxy_zzz_0, tdx_xxz_xxx_0, \
                                     tdx_xxz_xxy_0, tdx_xxz_xxz_0, tdx_xxz_xyy_0, tdx_xxz_xyz_0, tdx_xxz_xzz_0, \
                                     tdx_xxz_yyy_0, tdx_xxz_yyz_0, tdx_xxz_yzz_0, tdx_xxz_zzz_0, tdx_xy_yyz_0, \
                                     tdx_xy_yzz_0, tdx_xy_zzz_0, tdx_xyy_xxx_0, tdx_xyy_xxy_0, tdx_xyy_xxz_0, \
                                     tdx_xyy_xyy_0, tdx_xz_xx_0, tdx_xz_xxx_0, tdx_xz_xxy_0, tdx_xz_xxz_0, tdx_xz_xy_0, \
                                     tdx_xz_xyy_0, tdx_xz_xyz_0, tdx_xz_xz_0, tdx_xz_xzz_0, tdx_xz_yy_0, tdx_xz_yyy_0, \
                                     tdx_xz_yyz_0, tdx_xz_yz_0, tdx_xz_yzz_0, tdx_xz_zz_0, tdx_xz_zzz_0, tdx_y_yyz_0, \
                                     tdx_y_yzz_0, tdx_y_zzz_0, tdx_yy_xx_0, tdx_yy_xxx_0, tdx_yy_xxy_0, tdx_yy_xxz_0, \
                                     tdx_yy_xy_0, tdx_yy_xyy_0, tdx_yy_xz_0, tdx_yy_yy_0, tdx_z_xxx_0, tdx_z_xxy_0, \
                                     tdx_z_xxz_0, tdx_z_xyy_0, tdx_z_xyz_0, tdx_z_xzz_0, tdx_z_yyy_0, tdx_z_yyz_0, \
                                     tdx_z_yzz_0, tdx_z_zzz_0, tdy_xxy_yyz_0, tdy_xxy_yzz_0, tdy_xxy_zzz_0, \
                                     tdy_xxz_xxx_0, tdy_xxz_xxy_0, tdy_xxz_xxz_0, tdy_xxz_xyy_0, tdy_xxz_xyz_0, \
                                     tdy_xxz_xzz_0, tdy_xxz_yyy_0, tdy_xxz_yyz_0, tdy_xxz_yzz_0, tdy_xxz_zzz_0, \
                                     tdy_xy_yyz_0, tdy_xy_yzz_0, tdy_xy_zzz_0, tdy_xyy_xxx_0, tdy_xyy_xxy_0, \
                                     tdy_xyy_xxz_0, tdy_xz_xx_0, tdy_xz_xxx_0, tdy_xz_xxy_0, tdy_xz_xxz_0, tdy_xz_xy_0, \
                                     tdy_xz_xyy_0, tdy_xz_xyz_0, tdy_xz_xz_0, tdy_xz_xzz_0, tdy_xz_yy_0, tdy_xz_yyy_0, \
                                     tdy_xz_yyz_0, tdy_xz_yz_0, tdy_xz_yzz_0, tdy_xz_zz_0, tdy_xz_zzz_0, tdy_y_yyz_0, \
                                     tdy_y_yzz_0, tdy_y_zzz_0, tdy_yy_xx_0, tdy_yy_xxx_0, tdy_yy_xxy_0, tdy_yy_xxz_0, \
                                     tdy_yy_xy_0, tdy_yy_xz_0, tdy_z_xxx_0, tdy_z_xxy_0, tdy_z_xxz_0, tdy_z_xyy_0, \
                                     tdy_z_xyz_0, tdy_z_xzz_0, tdy_z_yyy_0, tdy_z_yyz_0, tdy_z_yzz_0, tdy_z_zzz_0, \
                                     tdz_xxy_yyy_0, tdz_xxy_yyz_0, tdz_xxy_yzz_0, tdz_xxy_zzz_0, tdz_xxz_xxx_0, \
                                     tdz_xxz_xxy_0, tdz_xxz_xxz_0, tdz_xxz_xyy_0, tdz_xxz_xyz_0, tdz_xxz_xzz_0, \
                                     tdz_xxz_yyy_0, tdz_xxz_yyz_0, tdz_xxz_yzz_0, tdz_xxz_zzz_0, tdz_xy_yyy_0, \
                                     tdz_xy_yyz_0, tdz_xy_yzz_0, tdz_xy_zzz_0, tdz_xyy_xxx_0, tdz_xyy_xxy_0, \
                                     tdz_xyy_xxz_0, tdz_xz_xx_0, tdz_xz_xxx_0, tdz_xz_xxy_0, tdz_xz_xxz_0, tdz_xz_xy_0, \
                                     tdz_xz_xyy_0, tdz_xz_xyz_0, tdz_xz_xz_0, tdz_xz_xzz_0, tdz_xz_yy_0, tdz_xz_yyy_0, \
                                     tdz_xz_yyz_0, tdz_xz_yz_0, tdz_xz_yzz_0, tdz_xz_zz_0, tdz_xz_zzz_0, tdz_y_yyy_0, \
                                     tdz_y_yyz_0, tdz_y_yzz_0, tdz_y_zzz_0, tdz_yy_xx_0, tdz_yy_xxx_0, tdz_yy_xxy_0, \
                                     tdz_yy_xxz_0, tdz_yy_xy_0, tdz_yy_xz_0, tdz_z_xxx_0, tdz_z_xxy_0, tdz_z_xxz_0, \
                                     tdz_z_xyy_0, tdz_z_xyz_0, tdz_z_xzz_0, tdz_z_yyy_0, tdz_z_yyz_0, tdz_z_yzz_0, \
                                     tdz_z_zzz_0, ts_xy_yyz_0, ts_xy_yzz_0, ts_xy_zzz_0, ts_xz_xxx_0, ts_xz_xxy_0, \
                                     ts_xz_xxz_0, ts_xz_xyy_0, ts_xz_xyz_0, ts_xz_xzz_0, ts_xz_yyy_0, ts_xz_yyz_0, \
                                     ts_xz_yzz_0, ts_xz_zzz_0, ts_yy_xxx_0, ts_yy_xxy_0, ts_yy_xxz_0, ts_yy_xyy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            tdz_xxy_yyy_0[j] = pa_x[j] * tdz_xy_yyy_0[j] + 0.5 * fl1_fx * tdz_y_yyy_0[j];

            tdx_xxy_yyz_0[j] = pa_x[j] * tdx_xy_yyz_0[j] + 0.5 * fl1_fx * tdx_y_yyz_0[j] + 0.5 * fl1_fx * ts_xy_yyz_0[j];

            tdy_xxy_yyz_0[j] = pa_x[j] * tdy_xy_yyz_0[j] + 0.5 * fl1_fx * tdy_y_yyz_0[j];

            tdz_xxy_yyz_0[j] = pa_x[j] * tdz_xy_yyz_0[j] + 0.5 * fl1_fx * tdz_y_yyz_0[j];

            tdx_xxy_yzz_0[j] = pa_x[j] * tdx_xy_yzz_0[j] + 0.5 * fl1_fx * tdx_y_yzz_0[j] + 0.5 * fl1_fx * ts_xy_yzz_0[j];

            tdy_xxy_yzz_0[j] = pa_x[j] * tdy_xy_yzz_0[j] + 0.5 * fl1_fx * tdy_y_yzz_0[j];

            tdz_xxy_yzz_0[j] = pa_x[j] * tdz_xy_yzz_0[j] + 0.5 * fl1_fx * tdz_y_yzz_0[j];

            tdx_xxy_zzz_0[j] = pa_x[j] * tdx_xy_zzz_0[j] + 0.5 * fl1_fx * tdx_y_zzz_0[j] + 0.5 * fl1_fx * ts_xy_zzz_0[j];

            tdy_xxy_zzz_0[j] = pa_x[j] * tdy_xy_zzz_0[j] + 0.5 * fl1_fx * tdy_y_zzz_0[j];

            tdz_xxy_zzz_0[j] = pa_x[j] * tdz_xy_zzz_0[j] + 0.5 * fl1_fx * tdz_y_zzz_0[j];

            tdx_xxz_xxx_0[j] =
                pa_x[j] * tdx_xz_xxx_0[j] + 0.5 * fl1_fx * tdx_z_xxx_0[j] + 1.5 * fl1_fx * tdx_xz_xx_0[j] + 0.5 * fl1_fx * ts_xz_xxx_0[j];

            tdy_xxz_xxx_0[j] = pa_x[j] * tdy_xz_xxx_0[j] + 0.5 * fl1_fx * tdy_z_xxx_0[j] + 1.5 * fl1_fx * tdy_xz_xx_0[j];

            tdz_xxz_xxx_0[j] = pa_x[j] * tdz_xz_xxx_0[j] + 0.5 * fl1_fx * tdz_z_xxx_0[j] + 1.5 * fl1_fx * tdz_xz_xx_0[j];

            tdx_xxz_xxy_0[j] = pa_x[j] * tdx_xz_xxy_0[j] + 0.5 * fl1_fx * tdx_z_xxy_0[j] + fl1_fx * tdx_xz_xy_0[j] + 0.5 * fl1_fx * ts_xz_xxy_0[j];

            tdy_xxz_xxy_0[j] = pa_x[j] * tdy_xz_xxy_0[j] + 0.5 * fl1_fx * tdy_z_xxy_0[j] + fl1_fx * tdy_xz_xy_0[j];

            tdz_xxz_xxy_0[j] = pa_x[j] * tdz_xz_xxy_0[j] + 0.5 * fl1_fx * tdz_z_xxy_0[j] + fl1_fx * tdz_xz_xy_0[j];

            tdx_xxz_xxz_0[j] = pa_x[j] * tdx_xz_xxz_0[j] + 0.5 * fl1_fx * tdx_z_xxz_0[j] + fl1_fx * tdx_xz_xz_0[j] + 0.5 * fl1_fx * ts_xz_xxz_0[j];

            tdy_xxz_xxz_0[j] = pa_x[j] * tdy_xz_xxz_0[j] + 0.5 * fl1_fx * tdy_z_xxz_0[j] + fl1_fx * tdy_xz_xz_0[j];

            tdz_xxz_xxz_0[j] = pa_x[j] * tdz_xz_xxz_0[j] + 0.5 * fl1_fx * tdz_z_xxz_0[j] + fl1_fx * tdz_xz_xz_0[j];

            tdx_xxz_xyy_0[j] =
                pa_x[j] * tdx_xz_xyy_0[j] + 0.5 * fl1_fx * tdx_z_xyy_0[j] + 0.5 * fl1_fx * tdx_xz_yy_0[j] + 0.5 * fl1_fx * ts_xz_xyy_0[j];

            tdy_xxz_xyy_0[j] = pa_x[j] * tdy_xz_xyy_0[j] + 0.5 * fl1_fx * tdy_z_xyy_0[j] + 0.5 * fl1_fx * tdy_xz_yy_0[j];

            tdz_xxz_xyy_0[j] = pa_x[j] * tdz_xz_xyy_0[j] + 0.5 * fl1_fx * tdz_z_xyy_0[j] + 0.5 * fl1_fx * tdz_xz_yy_0[j];

            tdx_xxz_xyz_0[j] =
                pa_x[j] * tdx_xz_xyz_0[j] + 0.5 * fl1_fx * tdx_z_xyz_0[j] + 0.5 * fl1_fx * tdx_xz_yz_0[j] + 0.5 * fl1_fx * ts_xz_xyz_0[j];

            tdy_xxz_xyz_0[j] = pa_x[j] * tdy_xz_xyz_0[j] + 0.5 * fl1_fx * tdy_z_xyz_0[j] + 0.5 * fl1_fx * tdy_xz_yz_0[j];

            tdz_xxz_xyz_0[j] = pa_x[j] * tdz_xz_xyz_0[j] + 0.5 * fl1_fx * tdz_z_xyz_0[j] + 0.5 * fl1_fx * tdz_xz_yz_0[j];

            tdx_xxz_xzz_0[j] =
                pa_x[j] * tdx_xz_xzz_0[j] + 0.5 * fl1_fx * tdx_z_xzz_0[j] + 0.5 * fl1_fx * tdx_xz_zz_0[j] + 0.5 * fl1_fx * ts_xz_xzz_0[j];

            tdy_xxz_xzz_0[j] = pa_x[j] * tdy_xz_xzz_0[j] + 0.5 * fl1_fx * tdy_z_xzz_0[j] + 0.5 * fl1_fx * tdy_xz_zz_0[j];

            tdz_xxz_xzz_0[j] = pa_x[j] * tdz_xz_xzz_0[j] + 0.5 * fl1_fx * tdz_z_xzz_0[j] + 0.5 * fl1_fx * tdz_xz_zz_0[j];

            tdx_xxz_yyy_0[j] = pa_x[j] * tdx_xz_yyy_0[j] + 0.5 * fl1_fx * tdx_z_yyy_0[j] + 0.5 * fl1_fx * ts_xz_yyy_0[j];

            tdy_xxz_yyy_0[j] = pa_x[j] * tdy_xz_yyy_0[j] + 0.5 * fl1_fx * tdy_z_yyy_0[j];

            tdz_xxz_yyy_0[j] = pa_x[j] * tdz_xz_yyy_0[j] + 0.5 * fl1_fx * tdz_z_yyy_0[j];

            tdx_xxz_yyz_0[j] = pa_x[j] * tdx_xz_yyz_0[j] + 0.5 * fl1_fx * tdx_z_yyz_0[j] + 0.5 * fl1_fx * ts_xz_yyz_0[j];

            tdy_xxz_yyz_0[j] = pa_x[j] * tdy_xz_yyz_0[j] + 0.5 * fl1_fx * tdy_z_yyz_0[j];

            tdz_xxz_yyz_0[j] = pa_x[j] * tdz_xz_yyz_0[j] + 0.5 * fl1_fx * tdz_z_yyz_0[j];

            tdx_xxz_yzz_0[j] = pa_x[j] * tdx_xz_yzz_0[j] + 0.5 * fl1_fx * tdx_z_yzz_0[j] + 0.5 * fl1_fx * ts_xz_yzz_0[j];

            tdy_xxz_yzz_0[j] = pa_x[j] * tdy_xz_yzz_0[j] + 0.5 * fl1_fx * tdy_z_yzz_0[j];

            tdz_xxz_yzz_0[j] = pa_x[j] * tdz_xz_yzz_0[j] + 0.5 * fl1_fx * tdz_z_yzz_0[j];

            tdx_xxz_zzz_0[j] = pa_x[j] * tdx_xz_zzz_0[j] + 0.5 * fl1_fx * tdx_z_zzz_0[j] + 0.5 * fl1_fx * ts_xz_zzz_0[j];

            tdy_xxz_zzz_0[j] = pa_x[j] * tdy_xz_zzz_0[j] + 0.5 * fl1_fx * tdy_z_zzz_0[j];

            tdz_xxz_zzz_0[j] = pa_x[j] * tdz_xz_zzz_0[j] + 0.5 * fl1_fx * tdz_z_zzz_0[j];

            tdx_xyy_xxx_0[j] = pa_x[j] * tdx_yy_xxx_0[j] + 1.5 * fl1_fx * tdx_yy_xx_0[j] + 0.5 * fl1_fx * ts_yy_xxx_0[j];

            tdy_xyy_xxx_0[j] = pa_x[j] * tdy_yy_xxx_0[j] + 1.5 * fl1_fx * tdy_yy_xx_0[j];

            tdz_xyy_xxx_0[j] = pa_x[j] * tdz_yy_xxx_0[j] + 1.5 * fl1_fx * tdz_yy_xx_0[j];

            tdx_xyy_xxy_0[j] = pa_x[j] * tdx_yy_xxy_0[j] + fl1_fx * tdx_yy_xy_0[j] + 0.5 * fl1_fx * ts_yy_xxy_0[j];

            tdy_xyy_xxy_0[j] = pa_x[j] * tdy_yy_xxy_0[j] + fl1_fx * tdy_yy_xy_0[j];

            tdz_xyy_xxy_0[j] = pa_x[j] * tdz_yy_xxy_0[j] + fl1_fx * tdz_yy_xy_0[j];

            tdx_xyy_xxz_0[j] = pa_x[j] * tdx_yy_xxz_0[j] + fl1_fx * tdx_yy_xz_0[j] + 0.5 * fl1_fx * ts_yy_xxz_0[j];

            tdy_xyy_xxz_0[j] = pa_x[j] * tdy_yy_xxz_0[j] + fl1_fx * tdy_yy_xz_0[j];

            tdz_xyy_xxz_0[j] = pa_x[j] * tdz_yy_xxz_0[j] + fl1_fx * tdz_yy_xz_0[j];

            tdx_xyy_xyy_0[j] = pa_x[j] * tdx_yy_xyy_0[j] + 0.5 * fl1_fx * tdx_yy_yy_0[j] + 0.5 * fl1_fx * ts_yy_xyy_0[j];
        }

        idx++;
    }
}

void
compElectricDipoleForFF_100_150(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_d_3_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

        auto tdy_yy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 33);

        auto tdz_yy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 33);

        auto tdx_yy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 34);

        auto tdy_yy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 34);

        auto tdz_yy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 34);

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

        auto tdy_yy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 21);

        auto tdz_yy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 21);

        auto tdx_yy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 22);

        auto tdy_yy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 22);

        auto tdz_yy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 22);

        auto tdx_yy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 23);

        auto tdy_yy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 23);

        auto tdz_yy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 23);

        auto tdx_yz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 24);

        auto tdy_yz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 24);

        auto tdz_yz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 24);

        auto tdx_yz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 25);

        auto tdy_yz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 25);

        auto tdz_yz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 25);

        auto tdx_yz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 26);

        auto tdy_yz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 26);

        auto tdz_yz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 26);

        auto tdx_yz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 27);

        auto tdy_yz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 27);

        auto tdz_yz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 27);

        auto tdx_yz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 28);

        auto tdy_yz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 28);

        auto tdz_yz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 28);

        auto tdx_yz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 29);

        auto tdy_yz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 29);

        auto tdz_yz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 29);

        auto ts_yy_xyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 34);

        auto ts_yy_xzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 35);

        auto ts_yy_yyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 36);

        auto ts_yy_yyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 37);

        auto ts_yy_yzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 38);

        auto ts_yy_zzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 39);

        auto ts_yz_xxx_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 40);

        auto ts_yz_xxy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 41);

        auto ts_yz_xxz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 42);

        auto ts_yz_xyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 43);

        auto ts_yz_xyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 44);

        auto ts_yz_xzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 45);

        auto ts_yz_yyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 46);

        auto ts_yz_yyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 47);

        auto ts_yz_yzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 48);

        auto ts_yz_zzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 49);

        // set up pointers to integrals

        auto tdy_xyy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 33);

        auto tdz_xyy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 33);

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

        // Batch of Integrals (100,150)

        #pragma omp simd aligned(fx, pa_x, tdx_xyy_xyz_0, tdx_xyy_xzz_0, tdx_xyy_yyy_0, tdx_xyy_yyz_0, \
                                     tdx_xyy_yzz_0, tdx_xyy_zzz_0, tdx_xyz_xxx_0, tdx_xyz_xxy_0, tdx_xyz_xxz_0, \
                                     tdx_xyz_xyy_0, tdx_xyz_xyz_0, tdx_xyz_xzz_0, tdx_xyz_yyy_0, tdx_xyz_yyz_0, \
                                     tdx_xyz_yzz_0, tdx_xyz_zzz_0, tdx_yy_xyz_0, tdx_yy_xzz_0, tdx_yy_yyy_0, \
                                     tdx_yy_yyz_0, tdx_yy_yz_0, tdx_yy_yzz_0, tdx_yy_zz_0, tdx_yy_zzz_0, tdx_yz_xx_0, \
                                     tdx_yz_xxx_0, tdx_yz_xxy_0, tdx_yz_xxz_0, tdx_yz_xy_0, tdx_yz_xyy_0, tdx_yz_xyz_0, \
                                     tdx_yz_xz_0, tdx_yz_xzz_0, tdx_yz_yy_0, tdx_yz_yyy_0, tdx_yz_yyz_0, tdx_yz_yz_0, \
                                     tdx_yz_yzz_0, tdx_yz_zz_0, tdx_yz_zzz_0, tdy_xyy_xyy_0, tdy_xyy_xyz_0, \
                                     tdy_xyy_xzz_0, tdy_xyy_yyy_0, tdy_xyy_yyz_0, tdy_xyy_yzz_0, tdy_xyy_zzz_0, \
                                     tdy_xyz_xxx_0, tdy_xyz_xxy_0, tdy_xyz_xxz_0, tdy_xyz_xyy_0, tdy_xyz_xyz_0, \
                                     tdy_xyz_xzz_0, tdy_xyz_yyy_0, tdy_xyz_yyz_0, tdy_xyz_yzz_0, tdy_xyz_zzz_0, \
                                     tdy_yy_xyy_0, tdy_yy_xyz_0, tdy_yy_xzz_0, tdy_yy_yy_0, tdy_yy_yyy_0, tdy_yy_yyz_0, \
                                     tdy_yy_yz_0, tdy_yy_yzz_0, tdy_yy_zz_0, tdy_yy_zzz_0, tdy_yz_xx_0, tdy_yz_xxx_0, \
                                     tdy_yz_xxy_0, tdy_yz_xxz_0, tdy_yz_xy_0, tdy_yz_xyy_0, tdy_yz_xyz_0, tdy_yz_xz_0, \
                                     tdy_yz_xzz_0, tdy_yz_yy_0, tdy_yz_yyy_0, tdy_yz_yyz_0, tdy_yz_yz_0, tdy_yz_yzz_0, \
                                     tdy_yz_zz_0, tdy_yz_zzz_0, tdz_xyy_xyy_0, tdz_xyy_xyz_0, tdz_xyy_xzz_0, \
                                     tdz_xyy_yyy_0, tdz_xyy_yyz_0, tdz_xyy_yzz_0, tdz_xyy_zzz_0, tdz_xyz_xxx_0, \
                                     tdz_xyz_xxy_0, tdz_xyz_xxz_0, tdz_xyz_xyy_0, tdz_xyz_xyz_0, tdz_xyz_xzz_0, \
                                     tdz_xyz_yyy_0, tdz_xyz_yyz_0, tdz_xyz_yzz_0, tdz_xyz_zzz_0, tdz_yy_xyy_0, \
                                     tdz_yy_xyz_0, tdz_yy_xzz_0, tdz_yy_yy_0, tdz_yy_yyy_0, tdz_yy_yyz_0, tdz_yy_yz_0, \
                                     tdz_yy_yzz_0, tdz_yy_zz_0, tdz_yy_zzz_0, tdz_yz_xx_0, tdz_yz_xxx_0, tdz_yz_xxy_0, \
                                     tdz_yz_xxz_0, tdz_yz_xy_0, tdz_yz_xyy_0, tdz_yz_xyz_0, tdz_yz_xz_0, tdz_yz_xzz_0, \
                                     tdz_yz_yy_0, tdz_yz_yyy_0, tdz_yz_yyz_0, tdz_yz_yz_0, tdz_yz_yzz_0, tdz_yz_zz_0, \
                                     tdz_yz_zzz_0, ts_yy_xyz_0, ts_yy_xzz_0, ts_yy_yyy_0, ts_yy_yyz_0, ts_yy_yzz_0, \
                                     ts_yy_zzz_0, ts_yz_xxx_0, ts_yz_xxy_0, ts_yz_xxz_0, ts_yz_xyy_0, ts_yz_xyz_0, \
                                     ts_yz_xzz_0, ts_yz_yyy_0, ts_yz_yyz_0, ts_yz_yzz_0, ts_yz_zzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            tdy_xyy_xyy_0[j] = pa_x[j] * tdy_yy_xyy_0[j] + 0.5 * fl1_fx * tdy_yy_yy_0[j];

            tdz_xyy_xyy_0[j] = pa_x[j] * tdz_yy_xyy_0[j] + 0.5 * fl1_fx * tdz_yy_yy_0[j];

            tdx_xyy_xyz_0[j] = pa_x[j] * tdx_yy_xyz_0[j] + 0.5 * fl1_fx * tdx_yy_yz_0[j] + 0.5 * fl1_fx * ts_yy_xyz_0[j];

            tdy_xyy_xyz_0[j] = pa_x[j] * tdy_yy_xyz_0[j] + 0.5 * fl1_fx * tdy_yy_yz_0[j];

            tdz_xyy_xyz_0[j] = pa_x[j] * tdz_yy_xyz_0[j] + 0.5 * fl1_fx * tdz_yy_yz_0[j];

            tdx_xyy_xzz_0[j] = pa_x[j] * tdx_yy_xzz_0[j] + 0.5 * fl1_fx * tdx_yy_zz_0[j] + 0.5 * fl1_fx * ts_yy_xzz_0[j];

            tdy_xyy_xzz_0[j] = pa_x[j] * tdy_yy_xzz_0[j] + 0.5 * fl1_fx * tdy_yy_zz_0[j];

            tdz_xyy_xzz_0[j] = pa_x[j] * tdz_yy_xzz_0[j] + 0.5 * fl1_fx * tdz_yy_zz_0[j];

            tdx_xyy_yyy_0[j] = pa_x[j] * tdx_yy_yyy_0[j] + 0.5 * fl1_fx * ts_yy_yyy_0[j];

            tdy_xyy_yyy_0[j] = pa_x[j] * tdy_yy_yyy_0[j];

            tdz_xyy_yyy_0[j] = pa_x[j] * tdz_yy_yyy_0[j];

            tdx_xyy_yyz_0[j] = pa_x[j] * tdx_yy_yyz_0[j] + 0.5 * fl1_fx * ts_yy_yyz_0[j];

            tdy_xyy_yyz_0[j] = pa_x[j] * tdy_yy_yyz_0[j];

            tdz_xyy_yyz_0[j] = pa_x[j] * tdz_yy_yyz_0[j];

            tdx_xyy_yzz_0[j] = pa_x[j] * tdx_yy_yzz_0[j] + 0.5 * fl1_fx * ts_yy_yzz_0[j];

            tdy_xyy_yzz_0[j] = pa_x[j] * tdy_yy_yzz_0[j];

            tdz_xyy_yzz_0[j] = pa_x[j] * tdz_yy_yzz_0[j];

            tdx_xyy_zzz_0[j] = pa_x[j] * tdx_yy_zzz_0[j] + 0.5 * fl1_fx * ts_yy_zzz_0[j];

            tdy_xyy_zzz_0[j] = pa_x[j] * tdy_yy_zzz_0[j];

            tdz_xyy_zzz_0[j] = pa_x[j] * tdz_yy_zzz_0[j];

            tdx_xyz_xxx_0[j] = pa_x[j] * tdx_yz_xxx_0[j] + 1.5 * fl1_fx * tdx_yz_xx_0[j] + 0.5 * fl1_fx * ts_yz_xxx_0[j];

            tdy_xyz_xxx_0[j] = pa_x[j] * tdy_yz_xxx_0[j] + 1.5 * fl1_fx * tdy_yz_xx_0[j];

            tdz_xyz_xxx_0[j] = pa_x[j] * tdz_yz_xxx_0[j] + 1.5 * fl1_fx * tdz_yz_xx_0[j];

            tdx_xyz_xxy_0[j] = pa_x[j] * tdx_yz_xxy_0[j] + fl1_fx * tdx_yz_xy_0[j] + 0.5 * fl1_fx * ts_yz_xxy_0[j];

            tdy_xyz_xxy_0[j] = pa_x[j] * tdy_yz_xxy_0[j] + fl1_fx * tdy_yz_xy_0[j];

            tdz_xyz_xxy_0[j] = pa_x[j] * tdz_yz_xxy_0[j] + fl1_fx * tdz_yz_xy_0[j];

            tdx_xyz_xxz_0[j] = pa_x[j] * tdx_yz_xxz_0[j] + fl1_fx * tdx_yz_xz_0[j] + 0.5 * fl1_fx * ts_yz_xxz_0[j];

            tdy_xyz_xxz_0[j] = pa_x[j] * tdy_yz_xxz_0[j] + fl1_fx * tdy_yz_xz_0[j];

            tdz_xyz_xxz_0[j] = pa_x[j] * tdz_yz_xxz_0[j] + fl1_fx * tdz_yz_xz_0[j];

            tdx_xyz_xyy_0[j] = pa_x[j] * tdx_yz_xyy_0[j] + 0.5 * fl1_fx * tdx_yz_yy_0[j] + 0.5 * fl1_fx * ts_yz_xyy_0[j];

            tdy_xyz_xyy_0[j] = pa_x[j] * tdy_yz_xyy_0[j] + 0.5 * fl1_fx * tdy_yz_yy_0[j];

            tdz_xyz_xyy_0[j] = pa_x[j] * tdz_yz_xyy_0[j] + 0.5 * fl1_fx * tdz_yz_yy_0[j];

            tdx_xyz_xyz_0[j] = pa_x[j] * tdx_yz_xyz_0[j] + 0.5 * fl1_fx * tdx_yz_yz_0[j] + 0.5 * fl1_fx * ts_yz_xyz_0[j];

            tdy_xyz_xyz_0[j] = pa_x[j] * tdy_yz_xyz_0[j] + 0.5 * fl1_fx * tdy_yz_yz_0[j];

            tdz_xyz_xyz_0[j] = pa_x[j] * tdz_yz_xyz_0[j] + 0.5 * fl1_fx * tdz_yz_yz_0[j];

            tdx_xyz_xzz_0[j] = pa_x[j] * tdx_yz_xzz_0[j] + 0.5 * fl1_fx * tdx_yz_zz_0[j] + 0.5 * fl1_fx * ts_yz_xzz_0[j];

            tdy_xyz_xzz_0[j] = pa_x[j] * tdy_yz_xzz_0[j] + 0.5 * fl1_fx * tdy_yz_zz_0[j];

            tdz_xyz_xzz_0[j] = pa_x[j] * tdz_yz_xzz_0[j] + 0.5 * fl1_fx * tdz_yz_zz_0[j];

            tdx_xyz_yyy_0[j] = pa_x[j] * tdx_yz_yyy_0[j] + 0.5 * fl1_fx * ts_yz_yyy_0[j];

            tdy_xyz_yyy_0[j] = pa_x[j] * tdy_yz_yyy_0[j];

            tdz_xyz_yyy_0[j] = pa_x[j] * tdz_yz_yyy_0[j];

            tdx_xyz_yyz_0[j] = pa_x[j] * tdx_yz_yyz_0[j] + 0.5 * fl1_fx * ts_yz_yyz_0[j];

            tdy_xyz_yyz_0[j] = pa_x[j] * tdy_yz_yyz_0[j];

            tdz_xyz_yyz_0[j] = pa_x[j] * tdz_yz_yyz_0[j];

            tdx_xyz_yzz_0[j] = pa_x[j] * tdx_yz_yzz_0[j] + 0.5 * fl1_fx * ts_yz_yzz_0[j];

            tdy_xyz_yzz_0[j] = pa_x[j] * tdy_yz_yzz_0[j];

            tdz_xyz_yzz_0[j] = pa_x[j] * tdz_yz_yzz_0[j];

            tdx_xyz_zzz_0[j] = pa_x[j] * tdx_yz_zzz_0[j] + 0.5 * fl1_fx * ts_yz_zzz_0[j];

            tdy_xyz_zzz_0[j] = pa_x[j] * tdy_yz_zzz_0[j];

            tdz_xyz_zzz_0[j] = pa_x[j] * tdz_yz_zzz_0[j];
        }

        idx++;
    }
}

void
compElectricDipoleForFF_150_200(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_d_3_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tdx_yy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 36);

        auto tdy_yy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 36);

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

        auto tdx_y_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 10);

        auto tdy_y_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 10);

        auto tdz_y_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 10);

        auto tdx_y_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 11);

        auto tdy_y_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 11);

        auto tdz_y_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 11);

        auto tdx_y_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 12);

        auto tdy_y_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 12);

        auto tdz_y_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 12);

        auto tdx_y_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 13);

        auto tdy_y_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 13);

        auto tdz_y_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 13);

        auto tdx_y_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 14);

        auto tdy_y_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 14);

        auto tdz_y_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 14);

        auto tdx_y_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 15);

        auto tdy_y_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 15);

        auto tdz_y_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 15);

        auto tdx_y_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 16);

        auto tdy_y_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 16);

        auto tdx_yy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 18);

        auto tdy_yy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 18);

        auto tdz_yy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 18);

        auto tdx_yy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 19);

        auto tdy_yy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 19);

        auto tdz_yy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 19);

        auto tdx_yy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 20);

        auto tdy_yy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 20);

        auto tdz_yy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 20);

        auto tdx_yy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 21);

        auto tdy_yy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 21);

        auto tdx_zz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 30);

        auto tdy_zz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 30);

        auto tdz_zz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 30);

        auto tdx_zz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 31);

        auto tdy_zz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 31);

        auto tdz_zz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 31);

        auto tdx_zz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 32);

        auto tdy_zz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 32);

        auto tdz_zz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 32);

        auto tdx_zz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 33);

        auto tdy_zz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 33);

        auto tdz_zz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 33);

        auto tdx_zz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 34);

        auto tdy_zz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 34);

        auto tdz_zz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 34);

        auto tdx_zz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 35);

        auto tdy_zz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 35);

        auto tdz_zz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 35);

        auto ts_yy_xxx_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 30);

        auto ts_yy_xxy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 31);

        auto ts_yy_xxz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 32);

        auto ts_yy_xyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 33);

        auto ts_yy_xyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 34);

        auto ts_yy_xzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 35);

        auto ts_yy_yyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 36);

        auto ts_zz_xxx_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 50);

        auto ts_zz_xxy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 51);

        auto ts_zz_xxz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 52);

        auto ts_zz_xyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 53);

        auto ts_zz_xyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 54);

        auto ts_zz_xzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 55);

        auto ts_zz_yyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 56);

        auto ts_zz_yyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 57);

        auto ts_zz_yzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 58);

        auto ts_zz_zzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 59);

        // set up pointers to integrals

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

        // Batch of Integrals (150,200)

        #pragma omp simd aligned(fx, pa_x, pa_y, tdx_xzz_xxx_0, tdx_xzz_xxy_0, tdx_xzz_xxz_0, \
                                     tdx_xzz_xyy_0, tdx_xzz_xyz_0, tdx_xzz_xzz_0, tdx_xzz_yyy_0, tdx_xzz_yyz_0, \
                                     tdx_xzz_yzz_0, tdx_xzz_zzz_0, tdx_y_xxx_0, tdx_y_xxy_0, tdx_y_xxz_0, tdx_y_xyy_0, \
                                     tdx_y_xyz_0, tdx_y_xzz_0, tdx_y_yyy_0, tdx_yy_xx_0, tdx_yy_xxx_0, tdx_yy_xxy_0, \
                                     tdx_yy_xxz_0, tdx_yy_xy_0, tdx_yy_xyy_0, tdx_yy_xyz_0, tdx_yy_xz_0, tdx_yy_xzz_0, \
                                     tdx_yy_yy_0, tdx_yy_yyy_0, tdx_yyy_xxx_0, tdx_yyy_xxy_0, tdx_yyy_xxz_0, \
                                     tdx_yyy_xyy_0, tdx_yyy_xyz_0, tdx_yyy_xzz_0, tdx_yyy_yyy_0, tdx_zz_xx_0, \
                                     tdx_zz_xxx_0, tdx_zz_xxy_0, tdx_zz_xxz_0, tdx_zz_xy_0, tdx_zz_xyy_0, tdx_zz_xyz_0, \
                                     tdx_zz_xz_0, tdx_zz_xzz_0, tdx_zz_yy_0, tdx_zz_yyy_0, tdx_zz_yyz_0, tdx_zz_yz_0, \
                                     tdx_zz_yzz_0, tdx_zz_zz_0, tdx_zz_zzz_0, tdy_xzz_xxx_0, tdy_xzz_xxy_0, \
                                     tdy_xzz_xxz_0, tdy_xzz_xyy_0, tdy_xzz_xyz_0, tdy_xzz_xzz_0, tdy_xzz_yyy_0, \
                                     tdy_xzz_yyz_0, tdy_xzz_yzz_0, tdy_xzz_zzz_0, tdy_y_xxx_0, tdy_y_xxy_0, tdy_y_xxz_0, \
                                     tdy_y_xyy_0, tdy_y_xyz_0, tdy_y_xzz_0, tdy_y_yyy_0, tdy_yy_xx_0, tdy_yy_xxx_0, \
                                     tdy_yy_xxy_0, tdy_yy_xxz_0, tdy_yy_xy_0, tdy_yy_xyy_0, tdy_yy_xyz_0, tdy_yy_xz_0, \
                                     tdy_yy_xzz_0, tdy_yy_yy_0, tdy_yy_yyy_0, tdy_yyy_xxx_0, tdy_yyy_xxy_0, \
                                     tdy_yyy_xxz_0, tdy_yyy_xyy_0, tdy_yyy_xyz_0, tdy_yyy_xzz_0, tdy_yyy_yyy_0, \
                                     tdy_zz_xx_0, tdy_zz_xxx_0, tdy_zz_xxy_0, tdy_zz_xxz_0, tdy_zz_xy_0, tdy_zz_xyy_0, \
                                     tdy_zz_xyz_0, tdy_zz_xz_0, tdy_zz_xzz_0, tdy_zz_yy_0, tdy_zz_yyy_0, tdy_zz_yyz_0, \
                                     tdy_zz_yz_0, tdy_zz_yzz_0, tdy_zz_zz_0, tdy_zz_zzz_0, tdz_xzz_xxx_0, \
                                     tdz_xzz_xxy_0, tdz_xzz_xxz_0, tdz_xzz_xyy_0, tdz_xzz_xyz_0, tdz_xzz_xzz_0, \
                                     tdz_xzz_yyy_0, tdz_xzz_yyz_0, tdz_xzz_yzz_0, tdz_xzz_zzz_0, tdz_y_xxx_0, \
                                     tdz_y_xxy_0, tdz_y_xxz_0, tdz_y_xyy_0, tdz_y_xyz_0, tdz_y_xzz_0, tdz_yy_xx_0, \
                                     tdz_yy_xxx_0, tdz_yy_xxy_0, tdz_yy_xxz_0, tdz_yy_xy_0, tdz_yy_xyy_0, tdz_yy_xyz_0, \
                                     tdz_yy_xz_0, tdz_yy_xzz_0, tdz_yyy_xxx_0, tdz_yyy_xxy_0, tdz_yyy_xxz_0, \
                                     tdz_yyy_xyy_0, tdz_yyy_xyz_0, tdz_yyy_xzz_0, tdz_zz_xx_0, tdz_zz_xxx_0, \
                                     tdz_zz_xxy_0, tdz_zz_xxz_0, tdz_zz_xy_0, tdz_zz_xyy_0, tdz_zz_xyz_0, tdz_zz_xz_0, \
                                     tdz_zz_xzz_0, tdz_zz_yy_0, tdz_zz_yyy_0, tdz_zz_yyz_0, tdz_zz_yz_0, tdz_zz_yzz_0, \
                                     tdz_zz_zz_0, tdz_zz_zzz_0, ts_yy_xxx_0, ts_yy_xxy_0, ts_yy_xxz_0, ts_yy_xyy_0, \
                                     ts_yy_xyz_0, ts_yy_xzz_0, ts_yy_yyy_0, ts_zz_xxx_0, ts_zz_xxy_0, ts_zz_xxz_0, \
                                     ts_zz_xyy_0, ts_zz_xyz_0, ts_zz_xzz_0, ts_zz_yyy_0, ts_zz_yyz_0, ts_zz_yzz_0, \
                                     ts_zz_zzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            tdx_xzz_xxx_0[j] = pa_x[j] * tdx_zz_xxx_0[j] + 1.5 * fl1_fx * tdx_zz_xx_0[j] + 0.5 * fl1_fx * ts_zz_xxx_0[j];

            tdy_xzz_xxx_0[j] = pa_x[j] * tdy_zz_xxx_0[j] + 1.5 * fl1_fx * tdy_zz_xx_0[j];

            tdz_xzz_xxx_0[j] = pa_x[j] * tdz_zz_xxx_0[j] + 1.5 * fl1_fx * tdz_zz_xx_0[j];

            tdx_xzz_xxy_0[j] = pa_x[j] * tdx_zz_xxy_0[j] + fl1_fx * tdx_zz_xy_0[j] + 0.5 * fl1_fx * ts_zz_xxy_0[j];

            tdy_xzz_xxy_0[j] = pa_x[j] * tdy_zz_xxy_0[j] + fl1_fx * tdy_zz_xy_0[j];

            tdz_xzz_xxy_0[j] = pa_x[j] * tdz_zz_xxy_0[j] + fl1_fx * tdz_zz_xy_0[j];

            tdx_xzz_xxz_0[j] = pa_x[j] * tdx_zz_xxz_0[j] + fl1_fx * tdx_zz_xz_0[j] + 0.5 * fl1_fx * ts_zz_xxz_0[j];

            tdy_xzz_xxz_0[j] = pa_x[j] * tdy_zz_xxz_0[j] + fl1_fx * tdy_zz_xz_0[j];

            tdz_xzz_xxz_0[j] = pa_x[j] * tdz_zz_xxz_0[j] + fl1_fx * tdz_zz_xz_0[j];

            tdx_xzz_xyy_0[j] = pa_x[j] * tdx_zz_xyy_0[j] + 0.5 * fl1_fx * tdx_zz_yy_0[j] + 0.5 * fl1_fx * ts_zz_xyy_0[j];

            tdy_xzz_xyy_0[j] = pa_x[j] * tdy_zz_xyy_0[j] + 0.5 * fl1_fx * tdy_zz_yy_0[j];

            tdz_xzz_xyy_0[j] = pa_x[j] * tdz_zz_xyy_0[j] + 0.5 * fl1_fx * tdz_zz_yy_0[j];

            tdx_xzz_xyz_0[j] = pa_x[j] * tdx_zz_xyz_0[j] + 0.5 * fl1_fx * tdx_zz_yz_0[j] + 0.5 * fl1_fx * ts_zz_xyz_0[j];

            tdy_xzz_xyz_0[j] = pa_x[j] * tdy_zz_xyz_0[j] + 0.5 * fl1_fx * tdy_zz_yz_0[j];

            tdz_xzz_xyz_0[j] = pa_x[j] * tdz_zz_xyz_0[j] + 0.5 * fl1_fx * tdz_zz_yz_0[j];

            tdx_xzz_xzz_0[j] = pa_x[j] * tdx_zz_xzz_0[j] + 0.5 * fl1_fx * tdx_zz_zz_0[j] + 0.5 * fl1_fx * ts_zz_xzz_0[j];

            tdy_xzz_xzz_0[j] = pa_x[j] * tdy_zz_xzz_0[j] + 0.5 * fl1_fx * tdy_zz_zz_0[j];

            tdz_xzz_xzz_0[j] = pa_x[j] * tdz_zz_xzz_0[j] + 0.5 * fl1_fx * tdz_zz_zz_0[j];

            tdx_xzz_yyy_0[j] = pa_x[j] * tdx_zz_yyy_0[j] + 0.5 * fl1_fx * ts_zz_yyy_0[j];

            tdy_xzz_yyy_0[j] = pa_x[j] * tdy_zz_yyy_0[j];

            tdz_xzz_yyy_0[j] = pa_x[j] * tdz_zz_yyy_0[j];

            tdx_xzz_yyz_0[j] = pa_x[j] * tdx_zz_yyz_0[j] + 0.5 * fl1_fx * ts_zz_yyz_0[j];

            tdy_xzz_yyz_0[j] = pa_x[j] * tdy_zz_yyz_0[j];

            tdz_xzz_yyz_0[j] = pa_x[j] * tdz_zz_yyz_0[j];

            tdx_xzz_yzz_0[j] = pa_x[j] * tdx_zz_yzz_0[j] + 0.5 * fl1_fx * ts_zz_yzz_0[j];

            tdy_xzz_yzz_0[j] = pa_x[j] * tdy_zz_yzz_0[j];

            tdz_xzz_yzz_0[j] = pa_x[j] * tdz_zz_yzz_0[j];

            tdx_xzz_zzz_0[j] = pa_x[j] * tdx_zz_zzz_0[j] + 0.5 * fl1_fx * ts_zz_zzz_0[j];

            tdy_xzz_zzz_0[j] = pa_x[j] * tdy_zz_zzz_0[j];

            tdz_xzz_zzz_0[j] = pa_x[j] * tdz_zz_zzz_0[j];

            tdx_yyy_xxx_0[j] = pa_y[j] * tdx_yy_xxx_0[j] + fl1_fx * tdx_y_xxx_0[j];

            tdy_yyy_xxx_0[j] = pa_y[j] * tdy_yy_xxx_0[j] + fl1_fx * tdy_y_xxx_0[j] + 0.5 * fl1_fx * ts_yy_xxx_0[j];

            tdz_yyy_xxx_0[j] = pa_y[j] * tdz_yy_xxx_0[j] + fl1_fx * tdz_y_xxx_0[j];

            tdx_yyy_xxy_0[j] = pa_y[j] * tdx_yy_xxy_0[j] + fl1_fx * tdx_y_xxy_0[j] + 0.5 * fl1_fx * tdx_yy_xx_0[j];

            tdy_yyy_xxy_0[j] = pa_y[j] * tdy_yy_xxy_0[j] + fl1_fx * tdy_y_xxy_0[j] + 0.5 * fl1_fx * tdy_yy_xx_0[j] + 0.5 * fl1_fx * ts_yy_xxy_0[j];

            tdz_yyy_xxy_0[j] = pa_y[j] * tdz_yy_xxy_0[j] + fl1_fx * tdz_y_xxy_0[j] + 0.5 * fl1_fx * tdz_yy_xx_0[j];

            tdx_yyy_xxz_0[j] = pa_y[j] * tdx_yy_xxz_0[j] + fl1_fx * tdx_y_xxz_0[j];

            tdy_yyy_xxz_0[j] = pa_y[j] * tdy_yy_xxz_0[j] + fl1_fx * tdy_y_xxz_0[j] + 0.5 * fl1_fx * ts_yy_xxz_0[j];

            tdz_yyy_xxz_0[j] = pa_y[j] * tdz_yy_xxz_0[j] + fl1_fx * tdz_y_xxz_0[j];

            tdx_yyy_xyy_0[j] = pa_y[j] * tdx_yy_xyy_0[j] + fl1_fx * tdx_y_xyy_0[j] + fl1_fx * tdx_yy_xy_0[j];

            tdy_yyy_xyy_0[j] = pa_y[j] * tdy_yy_xyy_0[j] + fl1_fx * tdy_y_xyy_0[j] + fl1_fx * tdy_yy_xy_0[j] + 0.5 * fl1_fx * ts_yy_xyy_0[j];

            tdz_yyy_xyy_0[j] = pa_y[j] * tdz_yy_xyy_0[j] + fl1_fx * tdz_y_xyy_0[j] + fl1_fx * tdz_yy_xy_0[j];

            tdx_yyy_xyz_0[j] = pa_y[j] * tdx_yy_xyz_0[j] + fl1_fx * tdx_y_xyz_0[j] + 0.5 * fl1_fx * tdx_yy_xz_0[j];

            tdy_yyy_xyz_0[j] = pa_y[j] * tdy_yy_xyz_0[j] + fl1_fx * tdy_y_xyz_0[j] + 0.5 * fl1_fx * tdy_yy_xz_0[j] + 0.5 * fl1_fx * ts_yy_xyz_0[j];

            tdz_yyy_xyz_0[j] = pa_y[j] * tdz_yy_xyz_0[j] + fl1_fx * tdz_y_xyz_0[j] + 0.5 * fl1_fx * tdz_yy_xz_0[j];

            tdx_yyy_xzz_0[j] = pa_y[j] * tdx_yy_xzz_0[j] + fl1_fx * tdx_y_xzz_0[j];

            tdy_yyy_xzz_0[j] = pa_y[j] * tdy_yy_xzz_0[j] + fl1_fx * tdy_y_xzz_0[j] + 0.5 * fl1_fx * ts_yy_xzz_0[j];

            tdz_yyy_xzz_0[j] = pa_y[j] * tdz_yy_xzz_0[j] + fl1_fx * tdz_y_xzz_0[j];

            tdx_yyy_yyy_0[j] = pa_y[j] * tdx_yy_yyy_0[j] + fl1_fx * tdx_y_yyy_0[j] + 1.5 * fl1_fx * tdx_yy_yy_0[j];

            tdy_yyy_yyy_0[j] = pa_y[j] * tdy_yy_yyy_0[j] + fl1_fx * tdy_y_yyy_0[j] + 1.5 * fl1_fx * tdy_yy_yy_0[j] + 0.5 * fl1_fx * ts_yy_yyy_0[j];
        }

        idx++;
    }
}

void
compElectricDipoleForFF_200_250(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_d_3_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_y = paDistances.data(3 * idx + 1);

        // set up pointers to auxilary integrals

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

        auto tdz_y_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 16);

        auto tdx_y_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 17);

        auto tdy_y_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 17);

        auto tdz_y_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 17);

        auto tdx_y_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 18);

        auto tdy_y_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 18);

        auto tdz_y_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 18);

        auto tdx_y_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 19);

        auto tdy_y_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 19);

        auto tdz_y_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 19);

        auto tdx_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 20);

        auto tdy_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 20);

        auto tdz_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 20);

        auto tdx_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 21);

        auto tdy_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 21);

        auto tdz_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 21);

        auto tdx_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 22);

        auto tdy_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 22);

        auto tdz_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 22);

        auto tdx_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 23);

        auto tdy_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 23);

        auto tdz_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 23);

        auto tdx_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 24);

        auto tdy_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 24);

        auto tdz_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 24);

        auto tdx_z_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 25);

        auto tdy_z_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 25);

        auto tdz_z_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 25);

        auto tdx_z_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 26);

        auto tdy_z_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 26);

        auto tdz_z_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 26);

        auto tdx_z_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 27);

        auto tdy_z_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 27);

        auto tdz_z_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 27);

        auto tdx_z_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 28);

        auto tdy_z_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 28);

        auto tdz_z_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 28);

        auto tdx_z_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 29);

        auto tdy_z_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 29);

        auto tdz_z_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 29);

        auto tdz_yy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 21);

        auto tdx_yy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 22);

        auto tdy_yy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 22);

        auto tdz_yy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 22);

        auto tdx_yy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 23);

        auto tdy_yy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 23);

        auto tdz_yy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 23);

        auto tdx_yz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 24);

        auto tdy_yz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 24);

        auto tdz_yz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 24);

        auto tdx_yz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 25);

        auto tdy_yz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 25);

        auto tdz_yz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 25);

        auto tdx_yz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 26);

        auto tdy_yz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 26);

        auto tdz_yz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 26);

        auto tdx_yz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 27);

        auto tdy_yz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 27);

        auto tdz_yz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 27);

        auto tdx_yz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 28);

        auto tdy_yz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 28);

        auto tdz_yz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 28);

        auto tdx_yz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 29);

        auto tdy_yz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 29);

        auto tdz_yz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 29);

        auto tdx_zz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 30);

        auto tdy_zz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 30);

        auto tdz_zz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 30);

        auto tdx_zz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 31);

        auto ts_yy_yyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 37);

        auto ts_yy_yzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 38);

        auto ts_yy_zzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 39);

        auto ts_yz_xxx_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 40);

        auto ts_yz_xxy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 41);

        auto ts_yz_xxz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 42);

        auto ts_yz_xyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 43);

        auto ts_yz_xyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 44);

        auto ts_yz_xzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 45);

        auto ts_yz_yyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 46);

        auto ts_yz_yyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 47);

        auto ts_yz_yzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 48);

        auto ts_yz_zzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 49);

        auto ts_zz_xxx_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 50);

        auto ts_zz_xxy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 51);

        auto ts_zz_xxz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 52);

        // set up pointers to integrals

        auto tdz_yyy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 66);

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

        // Batch of Integrals (200,250)

        #pragma omp simd aligned(fx, pa_y, tdx_y_yyz_0, tdx_y_yzz_0, tdx_y_zzz_0, tdx_yy_yyz_0, \
                                     tdx_yy_yz_0, tdx_yy_yzz_0, tdx_yy_zz_0, tdx_yy_zzz_0, tdx_yyy_yyz_0, \
                                     tdx_yyy_yzz_0, tdx_yyy_zzz_0, tdx_yyz_xxx_0, tdx_yyz_xxy_0, tdx_yyz_xxz_0, \
                                     tdx_yyz_xyy_0, tdx_yyz_xyz_0, tdx_yyz_xzz_0, tdx_yyz_yyy_0, tdx_yyz_yyz_0, \
                                     tdx_yyz_yzz_0, tdx_yyz_zzz_0, tdx_yz_xx_0, tdx_yz_xxx_0, tdx_yz_xxy_0, tdx_yz_xxz_0, \
                                     tdx_yz_xy_0, tdx_yz_xyy_0, tdx_yz_xyz_0, tdx_yz_xz_0, tdx_yz_xzz_0, tdx_yz_yy_0, \
                                     tdx_yz_yyy_0, tdx_yz_yyz_0, tdx_yz_yz_0, tdx_yz_yzz_0, tdx_yz_zz_0, tdx_yz_zzz_0, \
                                     tdx_yzz_xxx_0, tdx_yzz_xxy_0, tdx_yzz_xxz_0, tdx_yzz_xyy_0, tdx_z_xxx_0, \
                                     tdx_z_xxy_0, tdx_z_xxz_0, tdx_z_xyy_0, tdx_z_xyz_0, tdx_z_xzz_0, tdx_z_yyy_0, \
                                     tdx_z_yyz_0, tdx_z_yzz_0, tdx_z_zzz_0, tdx_zz_xx_0, tdx_zz_xxx_0, tdx_zz_xxy_0, \
                                     tdx_zz_xxz_0, tdx_zz_xy_0, tdx_zz_xyy_0, tdy_y_yyz_0, tdy_y_yzz_0, tdy_y_zzz_0, \
                                     tdy_yy_yyz_0, tdy_yy_yz_0, tdy_yy_yzz_0, tdy_yy_zz_0, tdy_yy_zzz_0, tdy_yyy_yyz_0, \
                                     tdy_yyy_yzz_0, tdy_yyy_zzz_0, tdy_yyz_xxx_0, tdy_yyz_xxy_0, tdy_yyz_xxz_0, \
                                     tdy_yyz_xyy_0, tdy_yyz_xyz_0, tdy_yyz_xzz_0, tdy_yyz_yyy_0, tdy_yyz_yyz_0, \
                                     tdy_yyz_yzz_0, tdy_yyz_zzz_0, tdy_yz_xx_0, tdy_yz_xxx_0, tdy_yz_xxy_0, tdy_yz_xxz_0, \
                                     tdy_yz_xy_0, tdy_yz_xyy_0, tdy_yz_xyz_0, tdy_yz_xz_0, tdy_yz_xzz_0, tdy_yz_yy_0, \
                                     tdy_yz_yyy_0, tdy_yz_yyz_0, tdy_yz_yz_0, tdy_yz_yzz_0, tdy_yz_zz_0, tdy_yz_zzz_0, \
                                     tdy_yzz_xxx_0, tdy_yzz_xxy_0, tdy_yzz_xxz_0, tdy_z_xxx_0, tdy_z_xxy_0, tdy_z_xxz_0, \
                                     tdy_z_xyy_0, tdy_z_xyz_0, tdy_z_xzz_0, tdy_z_yyy_0, tdy_z_yyz_0, tdy_z_yzz_0, \
                                     tdy_z_zzz_0, tdy_zz_xx_0, tdy_zz_xxx_0, tdy_zz_xxy_0, tdy_zz_xxz_0, tdz_y_yyy_0, \
                                     tdz_y_yyz_0, tdz_y_yzz_0, tdz_y_zzz_0, tdz_yy_yy_0, tdz_yy_yyy_0, tdz_yy_yyz_0, \
                                     tdz_yy_yz_0, tdz_yy_yzz_0, tdz_yy_zz_0, tdz_yy_zzz_0, tdz_yyy_yyy_0, \
                                     tdz_yyy_yyz_0, tdz_yyy_yzz_0, tdz_yyy_zzz_0, tdz_yyz_xxx_0, tdz_yyz_xxy_0, \
                                     tdz_yyz_xxz_0, tdz_yyz_xyy_0, tdz_yyz_xyz_0, tdz_yyz_xzz_0, tdz_yyz_yyy_0, \
                                     tdz_yyz_yyz_0, tdz_yyz_yzz_0, tdz_yyz_zzz_0, tdz_yz_xx_0, tdz_yz_xxx_0, \
                                     tdz_yz_xxy_0, tdz_yz_xxz_0, tdz_yz_xy_0, tdz_yz_xyy_0, tdz_yz_xyz_0, tdz_yz_xz_0, \
                                     tdz_yz_xzz_0, tdz_yz_yy_0, tdz_yz_yyy_0, tdz_yz_yyz_0, tdz_yz_yz_0, tdz_yz_yzz_0, \
                                     tdz_yz_zz_0, tdz_yz_zzz_0, tdz_yzz_xxx_0, tdz_yzz_xxy_0, tdz_yzz_xxz_0, \
                                     tdz_z_xxx_0, tdz_z_xxy_0, tdz_z_xxz_0, tdz_z_xyy_0, tdz_z_xyz_0, tdz_z_xzz_0, \
                                     tdz_z_yyy_0, tdz_z_yyz_0, tdz_z_yzz_0, tdz_z_zzz_0, tdz_zz_xx_0, tdz_zz_xxx_0, \
                                     tdz_zz_xxy_0, tdz_zz_xxz_0, ts_yy_yyz_0, ts_yy_yzz_0, ts_yy_zzz_0, ts_yz_xxx_0, \
                                     ts_yz_xxy_0, ts_yz_xxz_0, ts_yz_xyy_0, ts_yz_xyz_0, ts_yz_xzz_0, ts_yz_yyy_0, \
                                     ts_yz_yyz_0, ts_yz_yzz_0, ts_yz_zzz_0, ts_zz_xxx_0, ts_zz_xxy_0, ts_zz_xxz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            tdz_yyy_yyy_0[j] = pa_y[j] * tdz_yy_yyy_0[j] + fl1_fx * tdz_y_yyy_0[j] + 1.5 * fl1_fx * tdz_yy_yy_0[j];

            tdx_yyy_yyz_0[j] = pa_y[j] * tdx_yy_yyz_0[j] + fl1_fx * tdx_y_yyz_0[j] + fl1_fx * tdx_yy_yz_0[j];

            tdy_yyy_yyz_0[j] = pa_y[j] * tdy_yy_yyz_0[j] + fl1_fx * tdy_y_yyz_0[j] + fl1_fx * tdy_yy_yz_0[j] + 0.5 * fl1_fx * ts_yy_yyz_0[j];

            tdz_yyy_yyz_0[j] = pa_y[j] * tdz_yy_yyz_0[j] + fl1_fx * tdz_y_yyz_0[j] + fl1_fx * tdz_yy_yz_0[j];

            tdx_yyy_yzz_0[j] = pa_y[j] * tdx_yy_yzz_0[j] + fl1_fx * tdx_y_yzz_0[j] + 0.5 * fl1_fx * tdx_yy_zz_0[j];

            tdy_yyy_yzz_0[j] = pa_y[j] * tdy_yy_yzz_0[j] + fl1_fx * tdy_y_yzz_0[j] + 0.5 * fl1_fx * tdy_yy_zz_0[j] + 0.5 * fl1_fx * ts_yy_yzz_0[j];

            tdz_yyy_yzz_0[j] = pa_y[j] * tdz_yy_yzz_0[j] + fl1_fx * tdz_y_yzz_0[j] + 0.5 * fl1_fx * tdz_yy_zz_0[j];

            tdx_yyy_zzz_0[j] = pa_y[j] * tdx_yy_zzz_0[j] + fl1_fx * tdx_y_zzz_0[j];

            tdy_yyy_zzz_0[j] = pa_y[j] * tdy_yy_zzz_0[j] + fl1_fx * tdy_y_zzz_0[j] + 0.5 * fl1_fx * ts_yy_zzz_0[j];

            tdz_yyy_zzz_0[j] = pa_y[j] * tdz_yy_zzz_0[j] + fl1_fx * tdz_y_zzz_0[j];

            tdx_yyz_xxx_0[j] = pa_y[j] * tdx_yz_xxx_0[j] + 0.5 * fl1_fx * tdx_z_xxx_0[j];

            tdy_yyz_xxx_0[j] = pa_y[j] * tdy_yz_xxx_0[j] + 0.5 * fl1_fx * tdy_z_xxx_0[j] + 0.5 * fl1_fx * ts_yz_xxx_0[j];

            tdz_yyz_xxx_0[j] = pa_y[j] * tdz_yz_xxx_0[j] + 0.5 * fl1_fx * tdz_z_xxx_0[j];

            tdx_yyz_xxy_0[j] = pa_y[j] * tdx_yz_xxy_0[j] + 0.5 * fl1_fx * tdx_z_xxy_0[j] + 0.5 * fl1_fx * tdx_yz_xx_0[j];

            tdy_yyz_xxy_0[j] =
                pa_y[j] * tdy_yz_xxy_0[j] + 0.5 * fl1_fx * tdy_z_xxy_0[j] + 0.5 * fl1_fx * tdy_yz_xx_0[j] + 0.5 * fl1_fx * ts_yz_xxy_0[j];

            tdz_yyz_xxy_0[j] = pa_y[j] * tdz_yz_xxy_0[j] + 0.5 * fl1_fx * tdz_z_xxy_0[j] + 0.5 * fl1_fx * tdz_yz_xx_0[j];

            tdx_yyz_xxz_0[j] = pa_y[j] * tdx_yz_xxz_0[j] + 0.5 * fl1_fx * tdx_z_xxz_0[j];

            tdy_yyz_xxz_0[j] = pa_y[j] * tdy_yz_xxz_0[j] + 0.5 * fl1_fx * tdy_z_xxz_0[j] + 0.5 * fl1_fx * ts_yz_xxz_0[j];

            tdz_yyz_xxz_0[j] = pa_y[j] * tdz_yz_xxz_0[j] + 0.5 * fl1_fx * tdz_z_xxz_0[j];

            tdx_yyz_xyy_0[j] = pa_y[j] * tdx_yz_xyy_0[j] + 0.5 * fl1_fx * tdx_z_xyy_0[j] + fl1_fx * tdx_yz_xy_0[j];

            tdy_yyz_xyy_0[j] = pa_y[j] * tdy_yz_xyy_0[j] + 0.5 * fl1_fx * tdy_z_xyy_0[j] + fl1_fx * tdy_yz_xy_0[j] + 0.5 * fl1_fx * ts_yz_xyy_0[j];

            tdz_yyz_xyy_0[j] = pa_y[j] * tdz_yz_xyy_0[j] + 0.5 * fl1_fx * tdz_z_xyy_0[j] + fl1_fx * tdz_yz_xy_0[j];

            tdx_yyz_xyz_0[j] = pa_y[j] * tdx_yz_xyz_0[j] + 0.5 * fl1_fx * tdx_z_xyz_0[j] + 0.5 * fl1_fx * tdx_yz_xz_0[j];

            tdy_yyz_xyz_0[j] =
                pa_y[j] * tdy_yz_xyz_0[j] + 0.5 * fl1_fx * tdy_z_xyz_0[j] + 0.5 * fl1_fx * tdy_yz_xz_0[j] + 0.5 * fl1_fx * ts_yz_xyz_0[j];

            tdz_yyz_xyz_0[j] = pa_y[j] * tdz_yz_xyz_0[j] + 0.5 * fl1_fx * tdz_z_xyz_0[j] + 0.5 * fl1_fx * tdz_yz_xz_0[j];

            tdx_yyz_xzz_0[j] = pa_y[j] * tdx_yz_xzz_0[j] + 0.5 * fl1_fx * tdx_z_xzz_0[j];

            tdy_yyz_xzz_0[j] = pa_y[j] * tdy_yz_xzz_0[j] + 0.5 * fl1_fx * tdy_z_xzz_0[j] + 0.5 * fl1_fx * ts_yz_xzz_0[j];

            tdz_yyz_xzz_0[j] = pa_y[j] * tdz_yz_xzz_0[j] + 0.5 * fl1_fx * tdz_z_xzz_0[j];

            tdx_yyz_yyy_0[j] = pa_y[j] * tdx_yz_yyy_0[j] + 0.5 * fl1_fx * tdx_z_yyy_0[j] + 1.5 * fl1_fx * tdx_yz_yy_0[j];

            tdy_yyz_yyy_0[j] =
                pa_y[j] * tdy_yz_yyy_0[j] + 0.5 * fl1_fx * tdy_z_yyy_0[j] + 1.5 * fl1_fx * tdy_yz_yy_0[j] + 0.5 * fl1_fx * ts_yz_yyy_0[j];

            tdz_yyz_yyy_0[j] = pa_y[j] * tdz_yz_yyy_0[j] + 0.5 * fl1_fx * tdz_z_yyy_0[j] + 1.5 * fl1_fx * tdz_yz_yy_0[j];

            tdx_yyz_yyz_0[j] = pa_y[j] * tdx_yz_yyz_0[j] + 0.5 * fl1_fx * tdx_z_yyz_0[j] + fl1_fx * tdx_yz_yz_0[j];

            tdy_yyz_yyz_0[j] = pa_y[j] * tdy_yz_yyz_0[j] + 0.5 * fl1_fx * tdy_z_yyz_0[j] + fl1_fx * tdy_yz_yz_0[j] + 0.5 * fl1_fx * ts_yz_yyz_0[j];

            tdz_yyz_yyz_0[j] = pa_y[j] * tdz_yz_yyz_0[j] + 0.5 * fl1_fx * tdz_z_yyz_0[j] + fl1_fx * tdz_yz_yz_0[j];

            tdx_yyz_yzz_0[j] = pa_y[j] * tdx_yz_yzz_0[j] + 0.5 * fl1_fx * tdx_z_yzz_0[j] + 0.5 * fl1_fx * tdx_yz_zz_0[j];

            tdy_yyz_yzz_0[j] =
                pa_y[j] * tdy_yz_yzz_0[j] + 0.5 * fl1_fx * tdy_z_yzz_0[j] + 0.5 * fl1_fx * tdy_yz_zz_0[j] + 0.5 * fl1_fx * ts_yz_yzz_0[j];

            tdz_yyz_yzz_0[j] = pa_y[j] * tdz_yz_yzz_0[j] + 0.5 * fl1_fx * tdz_z_yzz_0[j] + 0.5 * fl1_fx * tdz_yz_zz_0[j];

            tdx_yyz_zzz_0[j] = pa_y[j] * tdx_yz_zzz_0[j] + 0.5 * fl1_fx * tdx_z_zzz_0[j];

            tdy_yyz_zzz_0[j] = pa_y[j] * tdy_yz_zzz_0[j] + 0.5 * fl1_fx * tdy_z_zzz_0[j] + 0.5 * fl1_fx * ts_yz_zzz_0[j];

            tdz_yyz_zzz_0[j] = pa_y[j] * tdz_yz_zzz_0[j] + 0.5 * fl1_fx * tdz_z_zzz_0[j];

            tdx_yzz_xxx_0[j] = pa_y[j] * tdx_zz_xxx_0[j];

            tdy_yzz_xxx_0[j] = pa_y[j] * tdy_zz_xxx_0[j] + 0.5 * fl1_fx * ts_zz_xxx_0[j];

            tdz_yzz_xxx_0[j] = pa_y[j] * tdz_zz_xxx_0[j];

            tdx_yzz_xxy_0[j] = pa_y[j] * tdx_zz_xxy_0[j] + 0.5 * fl1_fx * tdx_zz_xx_0[j];

            tdy_yzz_xxy_0[j] = pa_y[j] * tdy_zz_xxy_0[j] + 0.5 * fl1_fx * tdy_zz_xx_0[j] + 0.5 * fl1_fx * ts_zz_xxy_0[j];

            tdz_yzz_xxy_0[j] = pa_y[j] * tdz_zz_xxy_0[j] + 0.5 * fl1_fx * tdz_zz_xx_0[j];

            tdx_yzz_xxz_0[j] = pa_y[j] * tdx_zz_xxz_0[j];

            tdy_yzz_xxz_0[j] = pa_y[j] * tdy_zz_xxz_0[j] + 0.5 * fl1_fx * ts_zz_xxz_0[j];

            tdz_yzz_xxz_0[j] = pa_y[j] * tdz_zz_xxz_0[j];

            tdx_yzz_xyy_0[j] = pa_y[j] * tdx_zz_xyy_0[j] + fl1_fx * tdx_zz_xy_0[j];
        }

        idx++;
    }
}

void
compElectricDipoleForFF_250_300(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_d_3_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tdx_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 20);

        auto tdy_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 20);

        auto tdz_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 20);

        auto tdx_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 21);

        auto tdy_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 21);

        auto tdz_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 21);

        auto tdx_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 22);

        auto tdy_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 22);

        auto tdz_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 22);

        auto tdx_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 23);

        auto tdy_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 23);

        auto tdz_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 23);

        auto tdx_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 24);

        auto tdy_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 24);

        auto tdz_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 24);

        auto tdx_z_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 25);

        auto tdy_z_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 25);

        auto tdz_z_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 25);

        auto tdx_z_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 26);

        auto tdy_z_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 26);

        auto tdz_z_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 26);

        auto tdx_z_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 27);

        auto tdy_z_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 27);

        auto tdz_z_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 27);

        auto tdx_z_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 28);

        auto tdy_z_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 28);

        auto tdz_z_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 28);

        auto tdx_z_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 29);

        auto tdy_z_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 29);

        auto tdz_z_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 29);

        auto tdx_zz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 30);

        auto tdy_zz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 30);

        auto tdz_zz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 30);

        auto tdx_zz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 31);

        auto tdy_zz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 31);

        auto tdz_zz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 31);

        auto tdx_zz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 32);

        auto tdy_zz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 32);

        auto tdz_zz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 32);

        auto tdx_zz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 33);

        auto tdy_zz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 33);

        auto tdz_zz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 33);

        auto tdx_zz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 34);

        auto tdy_zz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 34);

        auto tdz_zz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 34);

        auto tdx_zz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 35);

        auto tdy_zz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 35);

        auto tdz_zz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 35);

        auto ts_zz_xxx_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 50);

        auto ts_zz_xxy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 51);

        auto ts_zz_xxz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 52);

        auto ts_zz_xyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 53);

        auto ts_zz_xyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 54);

        auto ts_zz_xzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 55);

        auto ts_zz_yyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 56);

        auto ts_zz_yyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 57);

        auto ts_zz_yzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 58);

        auto ts_zz_zzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 59);

        // set up pointers to integrals

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

        // Batch of Integrals (250,300)

        #pragma omp simd aligned(fx, pa_y, pa_z, tdx_yzz_xyz_0, tdx_yzz_xzz_0, tdx_yzz_yyy_0, \
                                     tdx_yzz_yyz_0, tdx_yzz_yzz_0, tdx_yzz_zzz_0, tdx_z_xxx_0, tdx_z_xxy_0, tdx_z_xxz_0, \
                                     tdx_z_xyy_0, tdx_z_xyz_0, tdx_z_xzz_0, tdx_z_yyy_0, tdx_z_yyz_0, tdx_z_yzz_0, \
                                     tdx_z_zzz_0, tdx_zz_xx_0, tdx_zz_xxx_0, tdx_zz_xxy_0, tdx_zz_xxz_0, tdx_zz_xy_0, \
                                     tdx_zz_xyy_0, tdx_zz_xyz_0, tdx_zz_xz_0, tdx_zz_xzz_0, tdx_zz_yy_0, tdx_zz_yyy_0, \
                                     tdx_zz_yyz_0, tdx_zz_yz_0, tdx_zz_yzz_0, tdx_zz_zz_0, tdx_zz_zzz_0, tdx_zzz_xxx_0, \
                                     tdx_zzz_xxy_0, tdx_zzz_xxz_0, tdx_zzz_xyy_0, tdx_zzz_xyz_0, tdx_zzz_xzz_0, \
                                     tdx_zzz_yyy_0, tdx_zzz_yyz_0, tdx_zzz_yzz_0, tdx_zzz_zzz_0, tdy_yzz_xyy_0, \
                                     tdy_yzz_xyz_0, tdy_yzz_xzz_0, tdy_yzz_yyy_0, tdy_yzz_yyz_0, tdy_yzz_yzz_0, \
                                     tdy_yzz_zzz_0, tdy_z_xxx_0, tdy_z_xxy_0, tdy_z_xxz_0, tdy_z_xyy_0, tdy_z_xyz_0, \
                                     tdy_z_xzz_0, tdy_z_yyy_0, tdy_z_yyz_0, tdy_z_yzz_0, tdy_z_zzz_0, tdy_zz_xx_0, \
                                     tdy_zz_xxx_0, tdy_zz_xxy_0, tdy_zz_xxz_0, tdy_zz_xy_0, tdy_zz_xyy_0, tdy_zz_xyz_0, \
                                     tdy_zz_xz_0, tdy_zz_xzz_0, tdy_zz_yy_0, tdy_zz_yyy_0, tdy_zz_yyz_0, tdy_zz_yz_0, \
                                     tdy_zz_yzz_0, tdy_zz_zz_0, tdy_zz_zzz_0, tdy_zzz_xxx_0, tdy_zzz_xxy_0, \
                                     tdy_zzz_xxz_0, tdy_zzz_xyy_0, tdy_zzz_xyz_0, tdy_zzz_xzz_0, tdy_zzz_yyy_0, \
                                     tdy_zzz_yyz_0, tdy_zzz_yzz_0, tdy_zzz_zzz_0, tdz_yzz_xyy_0, tdz_yzz_xyz_0, \
                                     tdz_yzz_xzz_0, tdz_yzz_yyy_0, tdz_yzz_yyz_0, tdz_yzz_yzz_0, tdz_yzz_zzz_0, \
                                     tdz_z_xxx_0, tdz_z_xxy_0, tdz_z_xxz_0, tdz_z_xyy_0, tdz_z_xyz_0, tdz_z_xzz_0, \
                                     tdz_z_yyy_0, tdz_z_yyz_0, tdz_z_yzz_0, tdz_z_zzz_0, tdz_zz_xx_0, tdz_zz_xxx_0, \
                                     tdz_zz_xxy_0, tdz_zz_xxz_0, tdz_zz_xy_0, tdz_zz_xyy_0, tdz_zz_xyz_0, tdz_zz_xz_0, \
                                     tdz_zz_xzz_0, tdz_zz_yy_0, tdz_zz_yyy_0, tdz_zz_yyz_0, tdz_zz_yz_0, tdz_zz_yzz_0, \
                                     tdz_zz_zz_0, tdz_zz_zzz_0, tdz_zzz_xxx_0, tdz_zzz_xxy_0, tdz_zzz_xxz_0, \
                                     tdz_zzz_xyy_0, tdz_zzz_xyz_0, tdz_zzz_xzz_0, tdz_zzz_yyy_0, tdz_zzz_yyz_0, \
                                     tdz_zzz_yzz_0, tdz_zzz_zzz_0, ts_zz_xxx_0, ts_zz_xxy_0, ts_zz_xxz_0, ts_zz_xyy_0, \
                                     ts_zz_xyz_0, ts_zz_xzz_0, ts_zz_yyy_0, ts_zz_yyz_0, ts_zz_yzz_0, ts_zz_zzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            tdy_yzz_xyy_0[j] = pa_y[j] * tdy_zz_xyy_0[j] + fl1_fx * tdy_zz_xy_0[j] + 0.5 * fl1_fx * ts_zz_xyy_0[j];

            tdz_yzz_xyy_0[j] = pa_y[j] * tdz_zz_xyy_0[j] + fl1_fx * tdz_zz_xy_0[j];

            tdx_yzz_xyz_0[j] = pa_y[j] * tdx_zz_xyz_0[j] + 0.5 * fl1_fx * tdx_zz_xz_0[j];

            tdy_yzz_xyz_0[j] = pa_y[j] * tdy_zz_xyz_0[j] + 0.5 * fl1_fx * tdy_zz_xz_0[j] + 0.5 * fl1_fx * ts_zz_xyz_0[j];

            tdz_yzz_xyz_0[j] = pa_y[j] * tdz_zz_xyz_0[j] + 0.5 * fl1_fx * tdz_zz_xz_0[j];

            tdx_yzz_xzz_0[j] = pa_y[j] * tdx_zz_xzz_0[j];

            tdy_yzz_xzz_0[j] = pa_y[j] * tdy_zz_xzz_0[j] + 0.5 * fl1_fx * ts_zz_xzz_0[j];

            tdz_yzz_xzz_0[j] = pa_y[j] * tdz_zz_xzz_0[j];

            tdx_yzz_yyy_0[j] = pa_y[j] * tdx_zz_yyy_0[j] + 1.5 * fl1_fx * tdx_zz_yy_0[j];

            tdy_yzz_yyy_0[j] = pa_y[j] * tdy_zz_yyy_0[j] + 1.5 * fl1_fx * tdy_zz_yy_0[j] + 0.5 * fl1_fx * ts_zz_yyy_0[j];

            tdz_yzz_yyy_0[j] = pa_y[j] * tdz_zz_yyy_0[j] + 1.5 * fl1_fx * tdz_zz_yy_0[j];

            tdx_yzz_yyz_0[j] = pa_y[j] * tdx_zz_yyz_0[j] + fl1_fx * tdx_zz_yz_0[j];

            tdy_yzz_yyz_0[j] = pa_y[j] * tdy_zz_yyz_0[j] + fl1_fx * tdy_zz_yz_0[j] + 0.5 * fl1_fx * ts_zz_yyz_0[j];

            tdz_yzz_yyz_0[j] = pa_y[j] * tdz_zz_yyz_0[j] + fl1_fx * tdz_zz_yz_0[j];

            tdx_yzz_yzz_0[j] = pa_y[j] * tdx_zz_yzz_0[j] + 0.5 * fl1_fx * tdx_zz_zz_0[j];

            tdy_yzz_yzz_0[j] = pa_y[j] * tdy_zz_yzz_0[j] + 0.5 * fl1_fx * tdy_zz_zz_0[j] + 0.5 * fl1_fx * ts_zz_yzz_0[j];

            tdz_yzz_yzz_0[j] = pa_y[j] * tdz_zz_yzz_0[j] + 0.5 * fl1_fx * tdz_zz_zz_0[j];

            tdx_yzz_zzz_0[j] = pa_y[j] * tdx_zz_zzz_0[j];

            tdy_yzz_zzz_0[j] = pa_y[j] * tdy_zz_zzz_0[j] + 0.5 * fl1_fx * ts_zz_zzz_0[j];

            tdz_yzz_zzz_0[j] = pa_y[j] * tdz_zz_zzz_0[j];

            tdx_zzz_xxx_0[j] = pa_z[j] * tdx_zz_xxx_0[j] + fl1_fx * tdx_z_xxx_0[j];

            tdy_zzz_xxx_0[j] = pa_z[j] * tdy_zz_xxx_0[j] + fl1_fx * tdy_z_xxx_0[j];

            tdz_zzz_xxx_0[j] = pa_z[j] * tdz_zz_xxx_0[j] + fl1_fx * tdz_z_xxx_0[j] + 0.5 * fl1_fx * ts_zz_xxx_0[j];

            tdx_zzz_xxy_0[j] = pa_z[j] * tdx_zz_xxy_0[j] + fl1_fx * tdx_z_xxy_0[j];

            tdy_zzz_xxy_0[j] = pa_z[j] * tdy_zz_xxy_0[j] + fl1_fx * tdy_z_xxy_0[j];

            tdz_zzz_xxy_0[j] = pa_z[j] * tdz_zz_xxy_0[j] + fl1_fx * tdz_z_xxy_0[j] + 0.5 * fl1_fx * ts_zz_xxy_0[j];

            tdx_zzz_xxz_0[j] = pa_z[j] * tdx_zz_xxz_0[j] + fl1_fx * tdx_z_xxz_0[j] + 0.5 * fl1_fx * tdx_zz_xx_0[j];

            tdy_zzz_xxz_0[j] = pa_z[j] * tdy_zz_xxz_0[j] + fl1_fx * tdy_z_xxz_0[j] + 0.5 * fl1_fx * tdy_zz_xx_0[j];

            tdz_zzz_xxz_0[j] = pa_z[j] * tdz_zz_xxz_0[j] + fl1_fx * tdz_z_xxz_0[j] + 0.5 * fl1_fx * tdz_zz_xx_0[j] + 0.5 * fl1_fx * ts_zz_xxz_0[j];

            tdx_zzz_xyy_0[j] = pa_z[j] * tdx_zz_xyy_0[j] + fl1_fx * tdx_z_xyy_0[j];

            tdy_zzz_xyy_0[j] = pa_z[j] * tdy_zz_xyy_0[j] + fl1_fx * tdy_z_xyy_0[j];

            tdz_zzz_xyy_0[j] = pa_z[j] * tdz_zz_xyy_0[j] + fl1_fx * tdz_z_xyy_0[j] + 0.5 * fl1_fx * ts_zz_xyy_0[j];

            tdx_zzz_xyz_0[j] = pa_z[j] * tdx_zz_xyz_0[j] + fl1_fx * tdx_z_xyz_0[j] + 0.5 * fl1_fx * tdx_zz_xy_0[j];

            tdy_zzz_xyz_0[j] = pa_z[j] * tdy_zz_xyz_0[j] + fl1_fx * tdy_z_xyz_0[j] + 0.5 * fl1_fx * tdy_zz_xy_0[j];

            tdz_zzz_xyz_0[j] = pa_z[j] * tdz_zz_xyz_0[j] + fl1_fx * tdz_z_xyz_0[j] + 0.5 * fl1_fx * tdz_zz_xy_0[j] + 0.5 * fl1_fx * ts_zz_xyz_0[j];

            tdx_zzz_xzz_0[j] = pa_z[j] * tdx_zz_xzz_0[j] + fl1_fx * tdx_z_xzz_0[j] + fl1_fx * tdx_zz_xz_0[j];

            tdy_zzz_xzz_0[j] = pa_z[j] * tdy_zz_xzz_0[j] + fl1_fx * tdy_z_xzz_0[j] + fl1_fx * tdy_zz_xz_0[j];

            tdz_zzz_xzz_0[j] = pa_z[j] * tdz_zz_xzz_0[j] + fl1_fx * tdz_z_xzz_0[j] + fl1_fx * tdz_zz_xz_0[j] + 0.5 * fl1_fx * ts_zz_xzz_0[j];

            tdx_zzz_yyy_0[j] = pa_z[j] * tdx_zz_yyy_0[j] + fl1_fx * tdx_z_yyy_0[j];

            tdy_zzz_yyy_0[j] = pa_z[j] * tdy_zz_yyy_0[j] + fl1_fx * tdy_z_yyy_0[j];

            tdz_zzz_yyy_0[j] = pa_z[j] * tdz_zz_yyy_0[j] + fl1_fx * tdz_z_yyy_0[j] + 0.5 * fl1_fx * ts_zz_yyy_0[j];

            tdx_zzz_yyz_0[j] = pa_z[j] * tdx_zz_yyz_0[j] + fl1_fx * tdx_z_yyz_0[j] + 0.5 * fl1_fx * tdx_zz_yy_0[j];

            tdy_zzz_yyz_0[j] = pa_z[j] * tdy_zz_yyz_0[j] + fl1_fx * tdy_z_yyz_0[j] + 0.5 * fl1_fx * tdy_zz_yy_0[j];

            tdz_zzz_yyz_0[j] = pa_z[j] * tdz_zz_yyz_0[j] + fl1_fx * tdz_z_yyz_0[j] + 0.5 * fl1_fx * tdz_zz_yy_0[j] + 0.5 * fl1_fx * ts_zz_yyz_0[j];

            tdx_zzz_yzz_0[j] = pa_z[j] * tdx_zz_yzz_0[j] + fl1_fx * tdx_z_yzz_0[j] + fl1_fx * tdx_zz_yz_0[j];

            tdy_zzz_yzz_0[j] = pa_z[j] * tdy_zz_yzz_0[j] + fl1_fx * tdy_z_yzz_0[j] + fl1_fx * tdy_zz_yz_0[j];

            tdz_zzz_yzz_0[j] = pa_z[j] * tdz_zz_yzz_0[j] + fl1_fx * tdz_z_yzz_0[j] + fl1_fx * tdz_zz_yz_0[j] + 0.5 * fl1_fx * ts_zz_yzz_0[j];

            tdx_zzz_zzz_0[j] = pa_z[j] * tdx_zz_zzz_0[j] + fl1_fx * tdx_z_zzz_0[j] + 1.5 * fl1_fx * tdx_zz_zz_0[j];

            tdy_zzz_zzz_0[j] = pa_z[j] * tdy_zz_zzz_0[j] + fl1_fx * tdy_z_zzz_0[j] + 1.5 * fl1_fx * tdy_zz_zz_0[j];

            tdz_zzz_zzz_0[j] = pa_z[j] * tdz_zz_zzz_0[j] + fl1_fx * tdz_z_zzz_0[j] + 1.5 * fl1_fx * tdz_zz_zz_0[j] + 0.5 * fl1_fx * ts_zz_zzz_0[j];
        }

        idx++;
    }
}

}  // namespace ediprecfunc
