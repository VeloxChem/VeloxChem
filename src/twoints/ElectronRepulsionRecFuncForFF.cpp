//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ElectronRepulsionRecFuncForFF.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSFSF(      CMemBlock2D<double>* primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        auto r_pb_z = braGtoPairsBlock.getDistancesPBZ();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {3, -1, -1, -1},
                                             {3, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_3_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_3_3_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_2_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_2_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_1_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_1_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_2_2_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
            {
                // set up pointers to Obara-Saika factors

                double fx = b_fx[i];

                auto fza = osFactors.data(4 * idx + 2);

                auto fxn = osFactors.data(4 * idx);

                // set up distances R(PB) = P - B

                auto pb_x = r_pb_x[i];

                auto pb_y = r_pb_y[i];

                auto pb_z = r_pb_z[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                auto wp_y = wpDistances.data(3 * idx + 1);

                auto wp_z = wpDistances.data(3 * idx + 2);

                // set up pointers to auxilary integrals

                auto tg_xx_xxx_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx); 

                auto tg_xx_xxy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 1); 

                auto tg_xx_xxz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 2); 

                auto tg_xx_xyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 3); 

                auto tg_xx_xyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 4); 

                auto tg_xx_xzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 5); 

                auto tg_xx_yyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 6); 

                auto tg_xx_yyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 7); 

                auto tg_xx_yzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 8); 

                auto tg_xx_zzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 9); 

                auto tg_xy_xxx_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 10); 

                auto tg_xy_xxy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 11); 

                auto tg_xy_xxz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 12); 

                auto tg_xy_xyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 13); 

                auto tg_xy_xyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 14); 

                auto tg_xy_xzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 15); 

                auto tg_xy_yyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 16); 

                auto tg_xy_yyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 17); 

                auto tg_xy_yzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 18); 

                auto tg_xy_zzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 19); 

                auto tg_xz_xxx_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 20); 

                auto tg_xz_xxy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 21); 

                auto tg_xz_xxz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 22); 

                auto tg_xz_xyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 23); 

                auto tg_xz_xyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 24); 

                auto tg_xz_xzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 25); 

                auto tg_xz_yyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 26); 

                auto tg_xz_yyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 27); 

                auto tg_xz_yzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 28); 

                auto tg_xz_zzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 29); 

                auto tg_yy_xxx_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 30); 

                auto tg_yy_xxy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 31); 

                auto tg_yy_xxz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 32); 

                auto tg_yy_xyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 33); 

                auto tg_yy_xyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 34); 

                auto tg_yy_xzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 35); 

                auto tg_yy_yyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 36); 

                auto tg_yy_yyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 37); 

                auto tg_yy_yzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 38); 

                auto tg_yy_zzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 39); 

                auto tg_yz_xxx_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 40); 

                auto tg_yz_xxy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 41); 

                auto tg_yz_xxz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 42); 

                auto tg_yz_xyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 43); 

                auto tg_yz_xyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 44); 

                auto tg_yz_xzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 45); 

                auto tg_yz_yyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 46); 

                auto tg_yz_yyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 47); 

                auto tg_yz_yzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 48); 

                auto tg_yz_zzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 49); 

                auto tg_zz_xxx_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 50); 

                auto tg_zz_xxy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 51); 

                auto tg_zz_xxz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 52); 

                auto tg_zz_xyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 53); 

                auto tg_zz_xyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 54); 

                auto tg_zz_xzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 55); 

                auto tg_zz_yyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 56); 

                auto tg_zz_yyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 57); 

                auto tg_zz_yzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 58); 

                auto tg_zz_zzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 59); 

                auto tg_xx_xxx_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx); 

                auto tg_xx_xxy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 1); 

                auto tg_xx_xxz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 2); 

                auto tg_xx_xyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 3); 

                auto tg_xx_xyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 4); 

                auto tg_xx_xzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 5); 

                auto tg_xx_yyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 6); 

                auto tg_xx_yyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 7); 

                auto tg_xx_yzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 8); 

                auto tg_xx_zzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 9); 

                auto tg_xy_xxx_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 10); 

                auto tg_xy_xxy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 11); 

                auto tg_xy_xxz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 12); 

                auto tg_xy_xyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 13); 

                auto tg_xy_xyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 14); 

                auto tg_xy_xzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 15); 

                auto tg_xy_yyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 16); 

                auto tg_xy_yyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 17); 

                auto tg_xy_yzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 18); 

                auto tg_xy_zzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 19); 

                auto tg_xz_xxx_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 20); 

                auto tg_xz_xxy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 21); 

                auto tg_xz_xxz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 22); 

                auto tg_xz_xyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 23); 

                auto tg_xz_xyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 24); 

                auto tg_xz_xzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 25); 

                auto tg_xz_yyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 26); 

                auto tg_xz_yyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 27); 

                auto tg_xz_yzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 28); 

                auto tg_xz_zzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 29); 

                auto tg_yy_xxx_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 30); 

                auto tg_yy_xxy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 31); 

                auto tg_yy_xxz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 32); 

                auto tg_yy_xyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 33); 

                auto tg_yy_xyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 34); 

                auto tg_yy_xzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 35); 

                auto tg_yy_yyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 36); 

                auto tg_yy_yyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 37); 

                auto tg_yy_yzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 38); 

                auto tg_yy_zzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 39); 

                auto tg_yz_xxx_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 40); 

                auto tg_yz_xxy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 41); 

                auto tg_yz_xxz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 42); 

                auto tg_yz_xyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 43); 

                auto tg_yz_xyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 44); 

                auto tg_yz_xzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 45); 

                auto tg_yz_yyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 46); 

                auto tg_yz_yyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 47); 

                auto tg_yz_yzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 48); 

                auto tg_yz_zzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 49); 

                auto tg_zz_xxx_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 50); 

                auto tg_zz_xxy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 51); 

                auto tg_zz_xxz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 52); 

                auto tg_zz_xyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 53); 

                auto tg_zz_xyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 54); 

                auto tg_zz_xzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 55); 

                auto tg_zz_yyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 56); 

                auto tg_zz_yyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 57); 

                auto tg_zz_yzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 58); 

                auto tg_zz_zzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 59); 

                auto tg_x_xxx_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx); 

                auto tg_x_xxy_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx + 1); 

                auto tg_x_xxz_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx + 2); 

                auto tg_x_xyy_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx + 3); 

                auto tg_x_xyz_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx + 4); 

                auto tg_x_xzz_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx + 5); 

                auto tg_x_yyy_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx + 6); 

                auto tg_x_yyz_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx + 7); 

                auto tg_x_yzz_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx + 8); 

                auto tg_x_zzz_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx + 9); 

                auto tg_y_xxx_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx + 10); 

                auto tg_y_xxy_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx + 11); 

                auto tg_y_xxz_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx + 12); 

                auto tg_y_xyy_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx + 13); 

                auto tg_y_xyz_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx + 14); 

                auto tg_y_xzz_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx + 15); 

                auto tg_y_yyy_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx + 16); 

                auto tg_y_yyz_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx + 17); 

                auto tg_y_yzz_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx + 18); 

                auto tg_y_zzz_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx + 19); 

                auto tg_z_xxx_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx + 20); 

                auto tg_z_xxy_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx + 21); 

                auto tg_z_xxz_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx + 22); 

                auto tg_z_xyy_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx + 23); 

                auto tg_z_xyz_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx + 24); 

                auto tg_z_xzz_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx + 25); 

                auto tg_z_yyy_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx + 26); 

                auto tg_z_yyz_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx + 27); 

                auto tg_z_yzz_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx + 28); 

                auto tg_z_zzz_0 = primBuffer[pidx_g_1_3_m0].data(30 * idx + 29); 

                auto tg_x_xxx_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx); 

                auto tg_x_xxy_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx + 1); 

                auto tg_x_xxz_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx + 2); 

                auto tg_x_xyy_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx + 3); 

                auto tg_x_xyz_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx + 4); 

                auto tg_x_xzz_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx + 5); 

                auto tg_x_yyy_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx + 6); 

                auto tg_x_yyz_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx + 7); 

                auto tg_x_yzz_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx + 8); 

                auto tg_x_zzz_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx + 9); 

                auto tg_y_xxx_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx + 10); 

                auto tg_y_xxy_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx + 11); 

                auto tg_y_xxz_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx + 12); 

                auto tg_y_xyy_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx + 13); 

                auto tg_y_xyz_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx + 14); 

                auto tg_y_xzz_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx + 15); 

                auto tg_y_yyy_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx + 16); 

                auto tg_y_yyz_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx + 17); 

                auto tg_y_yzz_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx + 18); 

                auto tg_y_zzz_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx + 19); 

                auto tg_z_xxx_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx + 20); 

                auto tg_z_xxy_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx + 21); 

                auto tg_z_xxz_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx + 22); 

                auto tg_z_xyy_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx + 23); 

                auto tg_z_xyz_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx + 24); 

                auto tg_z_xzz_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx + 25); 

                auto tg_z_yyy_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx + 26); 

                auto tg_z_yyz_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx + 27); 

                auto tg_z_yzz_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx + 28); 

                auto tg_z_zzz_1 = primBuffer[pidx_g_1_3_m1].data(30 * idx + 29); 

                auto tg_xx_xx_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx); 

                auto tg_xx_xy_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 1); 

                auto tg_xx_xz_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 2); 

                auto tg_xx_yy_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 3); 

                auto tg_xx_yz_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 4); 

                auto tg_xx_zz_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 5); 

                auto tg_xy_xx_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 6); 

                auto tg_xy_xy_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 7); 

                auto tg_xy_xz_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 8); 

                auto tg_xy_yy_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 9); 

                auto tg_xy_yz_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 10); 

                auto tg_xy_zz_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 11); 

                auto tg_xz_xx_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 12); 

                auto tg_xz_xy_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 13); 

                auto tg_xz_xz_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 14); 

                auto tg_xz_yy_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 15); 

                auto tg_xz_yz_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 16); 

                auto tg_xz_zz_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 17); 

                auto tg_yy_xx_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 18); 

                auto tg_yy_xy_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 19); 

                auto tg_yy_xz_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 20); 

                auto tg_yy_yy_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 21); 

                auto tg_yy_yz_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 22); 

                auto tg_yy_zz_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 23); 

                auto tg_yz_xx_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 24); 

                auto tg_yz_xy_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 25); 

                auto tg_yz_xz_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 26); 

                auto tg_yz_yy_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 27); 

                auto tg_yz_yz_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 28); 

                auto tg_yz_zz_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 29); 

                auto tg_zz_xx_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 30); 

                auto tg_zz_xy_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 31); 

                auto tg_zz_xz_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 32); 

                auto tg_zz_yy_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 33); 

                auto tg_zz_yz_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 34); 

                auto tg_zz_zz_1 = primBuffer[pidx_g_2_2_m1].data(36 * idx + 35); 

                // set up pointers to integrals

                auto tg_xxx_xxx_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx); 

                auto tg_xxx_xxy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 1); 

                auto tg_xxx_xxz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 2); 

                auto tg_xxx_xyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 3); 

                auto tg_xxx_xyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 4); 

                auto tg_xxx_xzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 5); 

                auto tg_xxx_yyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 6); 

                auto tg_xxx_yyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 7); 

                auto tg_xxx_yzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 8); 

                auto tg_xxx_zzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 9); 

                auto tg_xxy_xxx_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 10); 

                auto tg_xxy_xxy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 11); 

                auto tg_xxy_xxz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 12); 

                auto tg_xxy_xyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 13); 

                auto tg_xxy_xyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 14); 

                auto tg_xxy_xzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 15); 

                auto tg_xxy_yyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 16); 

                auto tg_xxy_yyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 17); 

                auto tg_xxy_yzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 18); 

                auto tg_xxy_zzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 19); 

                auto tg_xxz_xxx_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 20); 

                auto tg_xxz_xxy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 21); 

                auto tg_xxz_xxz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 22); 

                auto tg_xxz_xyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 23); 

                auto tg_xxz_xyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 24); 

                auto tg_xxz_xzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 25); 

                auto tg_xxz_yyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 26); 

                auto tg_xxz_yyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 27); 

                auto tg_xxz_yzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 28); 

                auto tg_xxz_zzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 29); 

                auto tg_xyy_xxx_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 30); 

                auto tg_xyy_xxy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 31); 

                auto tg_xyy_xxz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 32); 

                auto tg_xyy_xyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 33); 

                auto tg_xyy_xyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 34); 

                auto tg_xyy_xzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 35); 

                auto tg_xyy_yyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 36); 

                auto tg_xyy_yyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 37); 

                auto tg_xyy_yzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 38); 

                auto tg_xyy_zzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 39); 

                auto tg_xyz_xxx_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 40); 

                auto tg_xyz_xxy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 41); 

                auto tg_xyz_xxz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 42); 

                auto tg_xyz_xyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 43); 

                auto tg_xyz_xyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 44); 

                auto tg_xyz_xzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 45); 

                auto tg_xyz_yyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 46); 

                auto tg_xyz_yyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 47); 

                auto tg_xyz_yzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 48); 

                auto tg_xyz_zzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 49); 

                auto tg_xzz_xxx_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 50); 

                auto tg_xzz_xxy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 51); 

                auto tg_xzz_xxz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 52); 

                auto tg_xzz_xyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 53); 

                auto tg_xzz_xyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 54); 

                auto tg_xzz_xzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 55); 

                auto tg_xzz_yyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 56); 

                auto tg_xzz_yyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 57); 

                auto tg_xzz_yzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 58); 

                auto tg_xzz_zzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 59); 

                auto tg_yyy_xxx_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 60); 

                auto tg_yyy_xxy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 61); 

                auto tg_yyy_xxz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 62); 

                auto tg_yyy_xyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 63); 

                auto tg_yyy_xyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 64); 

                auto tg_yyy_xzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 65); 

                auto tg_yyy_yyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 66); 

                auto tg_yyy_yyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 67); 

                auto tg_yyy_yzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 68); 

                auto tg_yyy_zzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 69); 

                auto tg_yyz_xxx_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 70); 

                auto tg_yyz_xxy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 71); 

                auto tg_yyz_xxz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 72); 

                auto tg_yyz_xyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 73); 

                auto tg_yyz_xyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 74); 

                auto tg_yyz_xzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 75); 

                auto tg_yyz_yyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 76); 

                auto tg_yyz_yyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 77); 

                auto tg_yyz_yzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 78); 

                auto tg_yyz_zzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 79); 

                auto tg_yzz_xxx_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 80); 

                auto tg_yzz_xxy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 81); 

                auto tg_yzz_xxz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 82); 

                auto tg_yzz_xyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 83); 

                auto tg_yzz_xyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 84); 

                auto tg_yzz_xzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 85); 

                auto tg_yzz_yyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 86); 

                auto tg_yzz_yyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 87); 

                auto tg_yzz_yzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 88); 

                auto tg_yzz_zzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 89); 

                auto tg_zzz_xxx_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 90); 

                auto tg_zzz_xxy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 91); 

                auto tg_zzz_xxz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 92); 

                auto tg_zzz_xyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 93); 

                auto tg_zzz_xyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 94); 

                auto tg_zzz_xzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 95); 

                auto tg_zzz_yyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 96); 

                auto tg_zzz_yyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 97); 

                auto tg_zzz_yzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 98); 

                auto tg_zzz_zzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 99); 

                #pragma omp simd aligned(fxn, fza, tg_x_xxx_0, tg_x_xxx_1, tg_x_xxy_0, tg_x_xxy_1, tg_x_xxz_0, \
                                         tg_x_xxz_1, tg_x_xyy_0, tg_x_xyy_1, tg_x_xyz_0, tg_x_xyz_1, tg_x_xzz_0, tg_x_xzz_1, \
                                         tg_x_yyy_0, tg_x_yyy_1, tg_x_yyz_0, tg_x_yyz_1, tg_x_yzz_0, tg_x_yzz_1, tg_x_zzz_0, \
                                         tg_x_zzz_1, tg_xx_xx_1, tg_xx_xxx_0, tg_xx_xxx_1, tg_xx_xxy_0, tg_xx_xxy_1, \
                                         tg_xx_xxz_0, tg_xx_xxz_1, tg_xx_xy_1, tg_xx_xyy_0, tg_xx_xyy_1, tg_xx_xyz_0, \
                                         tg_xx_xyz_1, tg_xx_xz_1, tg_xx_xzz_0, tg_xx_xzz_1, tg_xx_yy_1, tg_xx_yyy_0, \
                                         tg_xx_yyy_1, tg_xx_yyz_0, tg_xx_yyz_1, tg_xx_yz_1, tg_xx_yzz_0, tg_xx_yzz_1, \
                                         tg_xx_zz_1, tg_xx_zzz_0, tg_xx_zzz_1, tg_xxx_xxx_0, tg_xxx_xxy_0, tg_xxx_xxz_0, \
                                         tg_xxx_xyy_0, tg_xxx_xyz_0, tg_xxx_xzz_0, tg_xxx_yyy_0, tg_xxx_yyz_0, tg_xxx_yzz_0, \
                                         tg_xxx_zzz_0, tg_xxy_xxx_0, tg_xxy_xxy_0, tg_xxy_xxz_0, tg_xxy_xyy_0, tg_xxy_xyz_0, \
                                         tg_xxy_xzz_0, tg_xxy_yyy_0, tg_xxy_yyz_0, tg_xxy_yzz_0, tg_xxy_zzz_0, tg_xxz_xxx_0, \
                                         tg_xxz_xxy_0, tg_xxz_xxz_0, tg_xxz_xyy_0, tg_xxz_xyz_0, tg_xxz_xzz_0, tg_xxz_yyy_0, \
                                         tg_xxz_yyz_0, tg_xxz_yzz_0, tg_xxz_zzz_0, tg_xy_xx_1, tg_xy_xxx_0, tg_xy_xxx_1, \
                                         tg_xy_xxy_0, tg_xy_xxy_1, tg_xy_xxz_0, tg_xy_xxz_1, tg_xy_xy_1, tg_xy_xyy_0, \
                                         tg_xy_xyy_1, tg_xy_xyz_0, tg_xy_xyz_1, tg_xy_xz_1, tg_xy_xzz_0, tg_xy_xzz_1, \
                                         tg_xy_yy_1, tg_xy_yyy_0, tg_xy_yyy_1, tg_xy_yyz_0, tg_xy_yyz_1, tg_xy_yz_1, \
                                         tg_xy_yzz_0, tg_xy_yzz_1, tg_xy_zz_1, tg_xy_zzz_0, tg_xy_zzz_1, tg_xyy_xxx_0, \
                                         tg_xyy_xxy_0, tg_xyy_xxz_0, tg_xyy_xyy_0, tg_xyy_xyz_0, tg_xyy_xzz_0, tg_xyy_yyy_0, \
                                         tg_xyy_yyz_0, tg_xyy_yzz_0, tg_xyy_zzz_0, tg_xyz_xxx_0, tg_xyz_xxy_0, tg_xyz_xxz_0, \
                                         tg_xyz_xyy_0, tg_xyz_xyz_0, tg_xyz_xzz_0, tg_xyz_yyy_0, tg_xyz_yyz_0, tg_xyz_yzz_0, \
                                         tg_xyz_zzz_0, tg_xz_xx_1, tg_xz_xxx_0, tg_xz_xxx_1, tg_xz_xxy_0, tg_xz_xxy_1, \
                                         tg_xz_xxz_0, tg_xz_xxz_1, tg_xz_xy_1, tg_xz_xyy_0, tg_xz_xyy_1, tg_xz_xyz_0, \
                                         tg_xz_xyz_1, tg_xz_xz_1, tg_xz_xzz_0, tg_xz_xzz_1, tg_xz_yy_1, tg_xz_yyy_0, \
                                         tg_xz_yyy_1, tg_xz_yyz_0, tg_xz_yyz_1, tg_xz_yz_1, tg_xz_yzz_0, tg_xz_yzz_1, \
                                         tg_xz_zz_1, tg_xz_zzz_0, tg_xz_zzz_1, tg_xzz_xxx_0, tg_xzz_xxy_0, tg_xzz_xxz_0, \
                                         tg_xzz_xyy_0, tg_xzz_xyz_0, tg_xzz_xzz_0, tg_xzz_yyy_0, tg_xzz_yyz_0, tg_xzz_yzz_0, \
                                         tg_xzz_zzz_0, tg_y_xxx_0, tg_y_xxx_1, tg_y_xxy_0, tg_y_xxy_1, tg_y_xxz_0, tg_y_xxz_1, \
                                         tg_y_xyy_0, tg_y_xyy_1, tg_y_xyz_0, tg_y_xyz_1, tg_y_xzz_0, tg_y_xzz_1, tg_y_yyy_0, \
                                         tg_y_yyy_1, tg_y_yyz_0, tg_y_yyz_1, tg_y_yzz_0, tg_y_yzz_1, tg_y_zzz_0, tg_y_zzz_1, \
                                         tg_yy_xx_1, tg_yy_xxx_0, tg_yy_xxx_1, tg_yy_xxy_0, tg_yy_xxy_1, tg_yy_xxz_0, \
                                         tg_yy_xxz_1, tg_yy_xy_1, tg_yy_xyy_0, tg_yy_xyy_1, tg_yy_xyz_0, tg_yy_xyz_1, \
                                         tg_yy_xz_1, tg_yy_xzz_0, tg_yy_xzz_1, tg_yy_yy_1, tg_yy_yyy_0, tg_yy_yyy_1, \
                                         tg_yy_yyz_0, tg_yy_yyz_1, tg_yy_yz_1, tg_yy_yzz_0, tg_yy_yzz_1, tg_yy_zz_1, \
                                         tg_yy_zzz_0, tg_yy_zzz_1, tg_yyy_xxx_0, tg_yyy_xxy_0, tg_yyy_xxz_0, tg_yyy_xyy_0, \
                                         tg_yyy_xyz_0, tg_yyy_xzz_0, tg_yyy_yyy_0, tg_yyy_yyz_0, tg_yyy_yzz_0, tg_yyy_zzz_0, \
                                         tg_yyz_xxx_0, tg_yyz_xxy_0, tg_yyz_xxz_0, tg_yyz_xyy_0, tg_yyz_xyz_0, tg_yyz_xzz_0, \
                                         tg_yyz_yyy_0, tg_yyz_yyz_0, tg_yyz_yzz_0, tg_yyz_zzz_0, tg_yz_xx_1, tg_yz_xxx_0, \
                                         tg_yz_xxx_1, tg_yz_xxy_0, tg_yz_xxy_1, tg_yz_xxz_0, tg_yz_xxz_1, tg_yz_xy_1, \
                                         tg_yz_xyy_0, tg_yz_xyy_1, tg_yz_xyz_0, tg_yz_xyz_1, tg_yz_xz_1, tg_yz_xzz_0, \
                                         tg_yz_xzz_1, tg_yz_yy_1, tg_yz_yyy_0, tg_yz_yyy_1, tg_yz_yyz_0, tg_yz_yyz_1, \
                                         tg_yz_yz_1, tg_yz_yzz_0, tg_yz_yzz_1, tg_yz_zz_1, tg_yz_zzz_0, tg_yz_zzz_1, \
                                         tg_yzz_xxx_0, tg_yzz_xxy_0, tg_yzz_xxz_0, tg_yzz_xyy_0, tg_yzz_xyz_0, tg_yzz_xzz_0, \
                                         tg_yzz_yyy_0, tg_yzz_yyz_0, tg_yzz_yzz_0, tg_yzz_zzz_0, tg_z_xxx_0, tg_z_xxx_1, \
                                         tg_z_xxy_0, tg_z_xxy_1, tg_z_xxz_0, tg_z_xxz_1, tg_z_xyy_0, tg_z_xyy_1, tg_z_xyz_0, \
                                         tg_z_xyz_1, tg_z_xzz_0, tg_z_xzz_1, tg_z_yyy_0, tg_z_yyy_1, tg_z_yyz_0, tg_z_yyz_1, \
                                         tg_z_yzz_0, tg_z_yzz_1, tg_z_zzz_0, tg_z_zzz_1, tg_zz_xx_1, tg_zz_xxx_0, \
                                         tg_zz_xxx_1, tg_zz_xxy_0, tg_zz_xxy_1, tg_zz_xxz_0, tg_zz_xxz_1, tg_zz_xy_1, \
                                         tg_zz_xyy_0, tg_zz_xyy_1, tg_zz_xyz_0, tg_zz_xyz_1, tg_zz_xz_1, tg_zz_xzz_0, \
                                         tg_zz_xzz_1, tg_zz_yy_1, tg_zz_yyy_0, tg_zz_yyy_1, tg_zz_yyz_0, tg_zz_yyz_1, \
                                         tg_zz_yz_1, tg_zz_yzz_0, tg_zz_yzz_1, tg_zz_zz_1, tg_zz_zzz_0, tg_zz_zzz_1, \
                                         tg_zzz_xxx_0, tg_zzz_xxy_0, tg_zzz_xxz_0, tg_zzz_xyy_0, tg_zzz_xyz_0, tg_zzz_xzz_0, \
                                         tg_zzz_yyy_0, tg_zzz_yyz_0, tg_zzz_yzz_0, tg_zzz_zzz_0, wp_x, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxx_xxx_0[j] = pb_x * tg_xx_xxx_0[j] + fr * tg_xx_xxx_1[j] + fl1_fx * (tg_x_xxx_0[j] - tg_x_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xx_xx_1[j];

                    tg_xxx_xxy_0[j] = pb_x * tg_xx_xxy_0[j] + fr * tg_xx_xxy_1[j] + fl1_fx * (tg_x_xxy_0[j] - tg_x_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xx_xy_1[j];

                    tg_xxx_xxz_0[j] = pb_x * tg_xx_xxz_0[j] + fr * tg_xx_xxz_1[j] + fl1_fx * (tg_x_xxz_0[j] - tg_x_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xx_xz_1[j];

                    tg_xxx_xyy_0[j] = pb_x * tg_xx_xyy_0[j] + fr * tg_xx_xyy_1[j] + fl1_fx * (tg_x_xyy_0[j] - tg_x_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xx_yy_1[j];

                    tg_xxx_xyz_0[j] = pb_x * tg_xx_xyz_0[j] + fr * tg_xx_xyz_1[j] + fl1_fx * (tg_x_xyz_0[j] - tg_x_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xx_yz_1[j];

                    tg_xxx_xzz_0[j] = pb_x * tg_xx_xzz_0[j] + fr * tg_xx_xzz_1[j] + fl1_fx * (tg_x_xzz_0[j] - tg_x_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xx_zz_1[j];

                    tg_xxx_yyy_0[j] = pb_x * tg_xx_yyy_0[j] + fr * tg_xx_yyy_1[j] + fl1_fx * (tg_x_yyy_0[j] - tg_x_yyy_1[j] * fl1_fza);

                    tg_xxx_yyz_0[j] = pb_x * tg_xx_yyz_0[j] + fr * tg_xx_yyz_1[j] + fl1_fx * (tg_x_yyz_0[j] - tg_x_yyz_1[j] * fl1_fza);

                    tg_xxx_yzz_0[j] = pb_x * tg_xx_yzz_0[j] + fr * tg_xx_yzz_1[j] + fl1_fx * (tg_x_yzz_0[j] - tg_x_yzz_1[j] * fl1_fza);

                    tg_xxx_zzz_0[j] = pb_x * tg_xx_zzz_0[j] + fr * tg_xx_zzz_1[j] + fl1_fx * (tg_x_zzz_0[j] - tg_x_zzz_1[j] * fl1_fza);

                    tg_xxy_xxx_0[j] = pb_x * tg_xy_xxx_0[j] + fr * tg_xy_xxx_1[j] + 0.5 * fl1_fx * (tg_y_xxx_0[j] - tg_y_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xy_xx_1[j];

                    tg_xxy_xxy_0[j] = pb_x * tg_xy_xxy_0[j] + fr * tg_xy_xxy_1[j] + 0.5 * fl1_fx * (tg_y_xxy_0[j] - tg_y_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xy_xy_1[j];

                    tg_xxy_xxz_0[j] = pb_x * tg_xy_xxz_0[j] + fr * tg_xy_xxz_1[j] + 0.5 * fl1_fx * (tg_y_xxz_0[j] - tg_y_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xy_xz_1[j];

                    tg_xxy_xyy_0[j] = pb_x * tg_xy_xyy_0[j] + fr * tg_xy_xyy_1[j] + 0.5 * fl1_fx * (tg_y_xyy_0[j] - tg_y_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xy_yy_1[j];

                    tg_xxy_xyz_0[j] = pb_x * tg_xy_xyz_0[j] + fr * tg_xy_xyz_1[j] + 0.5 * fl1_fx * (tg_y_xyz_0[j] - tg_y_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xy_yz_1[j];

                    tg_xxy_xzz_0[j] = pb_x * tg_xy_xzz_0[j] + fr * tg_xy_xzz_1[j] + 0.5 * fl1_fx * (tg_y_xzz_0[j] - tg_y_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xy_zz_1[j];

                    tg_xxy_yyy_0[j] = pb_x * tg_xy_yyy_0[j] + fr * tg_xy_yyy_1[j] + 0.5 * fl1_fx * (tg_y_yyy_0[j] - tg_y_yyy_1[j] * fl1_fza);

                    tg_xxy_yyz_0[j] = pb_x * tg_xy_yyz_0[j] + fr * tg_xy_yyz_1[j] + 0.5 * fl1_fx * (tg_y_yyz_0[j] - tg_y_yyz_1[j] * fl1_fza);

                    tg_xxy_yzz_0[j] = pb_x * tg_xy_yzz_0[j] + fr * tg_xy_yzz_1[j] + 0.5 * fl1_fx * (tg_y_yzz_0[j] - tg_y_yzz_1[j] * fl1_fza);

                    tg_xxy_zzz_0[j] = pb_x * tg_xy_zzz_0[j] + fr * tg_xy_zzz_1[j] + 0.5 * fl1_fx * (tg_y_zzz_0[j] - tg_y_zzz_1[j] * fl1_fza);

                    tg_xxz_xxx_0[j] = pb_x * tg_xz_xxx_0[j] + fr * tg_xz_xxx_1[j] + 0.5 * fl1_fx * (tg_z_xxx_0[j] - tg_z_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xz_xx_1[j];

                    tg_xxz_xxy_0[j] = pb_x * tg_xz_xxy_0[j] + fr * tg_xz_xxy_1[j] + 0.5 * fl1_fx * (tg_z_xxy_0[j] - tg_z_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xz_xy_1[j];

                    tg_xxz_xxz_0[j] = pb_x * tg_xz_xxz_0[j] + fr * tg_xz_xxz_1[j] + 0.5 * fl1_fx * (tg_z_xxz_0[j] - tg_z_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xz_xz_1[j];

                    tg_xxz_xyy_0[j] = pb_x * tg_xz_xyy_0[j] + fr * tg_xz_xyy_1[j] + 0.5 * fl1_fx * (tg_z_xyy_0[j] - tg_z_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xz_yy_1[j];

                    tg_xxz_xyz_0[j] = pb_x * tg_xz_xyz_0[j] + fr * tg_xz_xyz_1[j] + 0.5 * fl1_fx * (tg_z_xyz_0[j] - tg_z_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xz_yz_1[j];

                    tg_xxz_xzz_0[j] = pb_x * tg_xz_xzz_0[j] + fr * tg_xz_xzz_1[j] + 0.5 * fl1_fx * (tg_z_xzz_0[j] - tg_z_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xz_zz_1[j];

                    tg_xxz_yyy_0[j] = pb_x * tg_xz_yyy_0[j] + fr * tg_xz_yyy_1[j] + 0.5 * fl1_fx * (tg_z_yyy_0[j] - tg_z_yyy_1[j] * fl1_fza);

                    tg_xxz_yyz_0[j] = pb_x * tg_xz_yyz_0[j] + fr * tg_xz_yyz_1[j] + 0.5 * fl1_fx * (tg_z_yyz_0[j] - tg_z_yyz_1[j] * fl1_fza);

                    tg_xxz_yzz_0[j] = pb_x * tg_xz_yzz_0[j] + fr * tg_xz_yzz_1[j] + 0.5 * fl1_fx * (tg_z_yzz_0[j] - tg_z_yzz_1[j] * fl1_fza);

                    tg_xxz_zzz_0[j] = pb_x * tg_xz_zzz_0[j] + fr * tg_xz_zzz_1[j] + 0.5 * fl1_fx * (tg_z_zzz_0[j] - tg_z_zzz_1[j] * fl1_fza);

                    tg_xyy_xxx_0[j] = pb_x * tg_yy_xxx_0[j] + fr * tg_yy_xxx_1[j] + 1.5 * fl1_fxn * tg_yy_xx_1[j];

                    tg_xyy_xxy_0[j] = pb_x * tg_yy_xxy_0[j] + fr * tg_yy_xxy_1[j] + fl1_fxn * tg_yy_xy_1[j];

                    tg_xyy_xxz_0[j] = pb_x * tg_yy_xxz_0[j] + fr * tg_yy_xxz_1[j] + fl1_fxn * tg_yy_xz_1[j];

                    tg_xyy_xyy_0[j] = pb_x * tg_yy_xyy_0[j] + fr * tg_yy_xyy_1[j] + 0.5 * fl1_fxn * tg_yy_yy_1[j];

                    tg_xyy_xyz_0[j] = pb_x * tg_yy_xyz_0[j] + fr * tg_yy_xyz_1[j] + 0.5 * fl1_fxn * tg_yy_yz_1[j];

                    tg_xyy_xzz_0[j] = pb_x * tg_yy_xzz_0[j] + fr * tg_yy_xzz_1[j] + 0.5 * fl1_fxn * tg_yy_zz_1[j];

                    tg_xyy_yyy_0[j] = pb_x * tg_yy_yyy_0[j] + fr * tg_yy_yyy_1[j];

                    tg_xyy_yyz_0[j] = pb_x * tg_yy_yyz_0[j] + fr * tg_yy_yyz_1[j];

                    tg_xyy_yzz_0[j] = pb_x * tg_yy_yzz_0[j] + fr * tg_yy_yzz_1[j];

                    tg_xyy_zzz_0[j] = pb_x * tg_yy_zzz_0[j] + fr * tg_yy_zzz_1[j];

                    tg_xyz_xxx_0[j] = pb_x * tg_yz_xxx_0[j] + fr * tg_yz_xxx_1[j] + 1.5 * fl1_fxn * tg_yz_xx_1[j];

                    tg_xyz_xxy_0[j] = pb_x * tg_yz_xxy_0[j] + fr * tg_yz_xxy_1[j] + fl1_fxn * tg_yz_xy_1[j];

                    tg_xyz_xxz_0[j] = pb_x * tg_yz_xxz_0[j] + fr * tg_yz_xxz_1[j] + fl1_fxn * tg_yz_xz_1[j];

                    tg_xyz_xyy_0[j] = pb_x * tg_yz_xyy_0[j] + fr * tg_yz_xyy_1[j] + 0.5 * fl1_fxn * tg_yz_yy_1[j];

                    tg_xyz_xyz_0[j] = pb_x * tg_yz_xyz_0[j] + fr * tg_yz_xyz_1[j] + 0.5 * fl1_fxn * tg_yz_yz_1[j];

                    tg_xyz_xzz_0[j] = pb_x * tg_yz_xzz_0[j] + fr * tg_yz_xzz_1[j] + 0.5 * fl1_fxn * tg_yz_zz_1[j];

                    tg_xyz_yyy_0[j] = pb_x * tg_yz_yyy_0[j] + fr * tg_yz_yyy_1[j];

                    tg_xyz_yyz_0[j] = pb_x * tg_yz_yyz_0[j] + fr * tg_yz_yyz_1[j];

                    tg_xyz_yzz_0[j] = pb_x * tg_yz_yzz_0[j] + fr * tg_yz_yzz_1[j];

                    tg_xyz_zzz_0[j] = pb_x * tg_yz_zzz_0[j] + fr * tg_yz_zzz_1[j];

                    tg_xzz_xxx_0[j] = pb_x * tg_zz_xxx_0[j] + fr * tg_zz_xxx_1[j] + 1.5 * fl1_fxn * tg_zz_xx_1[j];

                    tg_xzz_xxy_0[j] = pb_x * tg_zz_xxy_0[j] + fr * tg_zz_xxy_1[j] + fl1_fxn * tg_zz_xy_1[j];

                    tg_xzz_xxz_0[j] = pb_x * tg_zz_xxz_0[j] + fr * tg_zz_xxz_1[j] + fl1_fxn * tg_zz_xz_1[j];

                    tg_xzz_xyy_0[j] = pb_x * tg_zz_xyy_0[j] + fr * tg_zz_xyy_1[j] + 0.5 * fl1_fxn * tg_zz_yy_1[j];

                    tg_xzz_xyz_0[j] = pb_x * tg_zz_xyz_0[j] + fr * tg_zz_xyz_1[j] + 0.5 * fl1_fxn * tg_zz_yz_1[j];

                    tg_xzz_xzz_0[j] = pb_x * tg_zz_xzz_0[j] + fr * tg_zz_xzz_1[j] + 0.5 * fl1_fxn * tg_zz_zz_1[j];

                    tg_xzz_yyy_0[j] = pb_x * tg_zz_yyy_0[j] + fr * tg_zz_yyy_1[j];

                    tg_xzz_yyz_0[j] = pb_x * tg_zz_yyz_0[j] + fr * tg_zz_yyz_1[j];

                    tg_xzz_yzz_0[j] = pb_x * tg_zz_yzz_0[j] + fr * tg_zz_yzz_1[j];

                    tg_xzz_zzz_0[j] = pb_x * tg_zz_zzz_0[j] + fr * tg_zz_zzz_1[j];

                    fr = wp_y[j]; 

                    tg_yyy_xxx_0[j] = pb_y * tg_yy_xxx_0[j] + fr * tg_yy_xxx_1[j] + fl1_fx * (tg_y_xxx_0[j] - tg_y_xxx_1[j] * fl1_fza);

                    tg_yyy_xxy_0[j] = pb_y * tg_yy_xxy_0[j] + fr * tg_yy_xxy_1[j] + fl1_fx * (tg_y_xxy_0[j] - tg_y_xxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yy_xx_1[j];

                    tg_yyy_xxz_0[j] = pb_y * tg_yy_xxz_0[j] + fr * tg_yy_xxz_1[j] + fl1_fx * (tg_y_xxz_0[j] - tg_y_xxz_1[j] * fl1_fza);

                    tg_yyy_xyy_0[j] = pb_y * tg_yy_xyy_0[j] + fr * tg_yy_xyy_1[j] + fl1_fx * (tg_y_xyy_0[j] - tg_y_xyy_1[j] * fl1_fza) + fl1_fxn * tg_yy_xy_1[j];

                    tg_yyy_xyz_0[j] = pb_y * tg_yy_xyz_0[j] + fr * tg_yy_xyz_1[j] + fl1_fx * (tg_y_xyz_0[j] - tg_y_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yy_xz_1[j];

                    tg_yyy_xzz_0[j] = pb_y * tg_yy_xzz_0[j] + fr * tg_yy_xzz_1[j] + fl1_fx * (tg_y_xzz_0[j] - tg_y_xzz_1[j] * fl1_fza);

                    tg_yyy_yyy_0[j] = pb_y * tg_yy_yyy_0[j] + fr * tg_yy_yyy_1[j] + fl1_fx * (tg_y_yyy_0[j] - tg_y_yyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yy_yy_1[j];

                    tg_yyy_yyz_0[j] = pb_y * tg_yy_yyz_0[j] + fr * tg_yy_yyz_1[j] + fl1_fx * (tg_y_yyz_0[j] - tg_y_yyz_1[j] * fl1_fza) + fl1_fxn * tg_yy_yz_1[j];

                    tg_yyy_yzz_0[j] = pb_y * tg_yy_yzz_0[j] + fr * tg_yy_yzz_1[j] + fl1_fx * (tg_y_yzz_0[j] - tg_y_yzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yy_zz_1[j];

                    tg_yyy_zzz_0[j] = pb_y * tg_yy_zzz_0[j] + fr * tg_yy_zzz_1[j] + fl1_fx * (tg_y_zzz_0[j] - tg_y_zzz_1[j] * fl1_fza);

                    tg_yyz_xxx_0[j] = pb_y * tg_yz_xxx_0[j] + fr * tg_yz_xxx_1[j] + 0.5 * fl1_fx * (tg_z_xxx_0[j] - tg_z_xxx_1[j] * fl1_fza);

                    tg_yyz_xxy_0[j] = pb_y * tg_yz_xxy_0[j] + fr * tg_yz_xxy_1[j] + 0.5 * fl1_fx * (tg_z_xxy_0[j] - tg_z_xxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yz_xx_1[j];

                    tg_yyz_xxz_0[j] = pb_y * tg_yz_xxz_0[j] + fr * tg_yz_xxz_1[j] + 0.5 * fl1_fx * (tg_z_xxz_0[j] - tg_z_xxz_1[j] * fl1_fza);

                    tg_yyz_xyy_0[j] = pb_y * tg_yz_xyy_0[j] + fr * tg_yz_xyy_1[j] + 0.5 * fl1_fx * (tg_z_xyy_0[j] - tg_z_xyy_1[j] * fl1_fza) + fl1_fxn * tg_yz_xy_1[j];

                    tg_yyz_xyz_0[j] = pb_y * tg_yz_xyz_0[j] + fr * tg_yz_xyz_1[j] + 0.5 * fl1_fx * (tg_z_xyz_0[j] - tg_z_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yz_xz_1[j];

                    tg_yyz_xzz_0[j] = pb_y * tg_yz_xzz_0[j] + fr * tg_yz_xzz_1[j] + 0.5 * fl1_fx * (tg_z_xzz_0[j] - tg_z_xzz_1[j] * fl1_fza);

                    tg_yyz_yyy_0[j] = pb_y * tg_yz_yyy_0[j] + fr * tg_yz_yyy_1[j] + 0.5 * fl1_fx * (tg_z_yyy_0[j] - tg_z_yyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yz_yy_1[j];

                    tg_yyz_yyz_0[j] = pb_y * tg_yz_yyz_0[j] + fr * tg_yz_yyz_1[j] + 0.5 * fl1_fx * (tg_z_yyz_0[j] - tg_z_yyz_1[j] * fl1_fza) + fl1_fxn * tg_yz_yz_1[j];

                    tg_yyz_yzz_0[j] = pb_y * tg_yz_yzz_0[j] + fr * tg_yz_yzz_1[j] + 0.5 * fl1_fx * (tg_z_yzz_0[j] - tg_z_yzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yz_zz_1[j];

                    tg_yyz_zzz_0[j] = pb_y * tg_yz_zzz_0[j] + fr * tg_yz_zzz_1[j] + 0.5 * fl1_fx * (tg_z_zzz_0[j] - tg_z_zzz_1[j] * fl1_fza);

                    tg_yzz_xxx_0[j] = pb_y * tg_zz_xxx_0[j] + fr * tg_zz_xxx_1[j];

                    tg_yzz_xxy_0[j] = pb_y * tg_zz_xxy_0[j] + fr * tg_zz_xxy_1[j] + 0.5 * fl1_fxn * tg_zz_xx_1[j];

                    tg_yzz_xxz_0[j] = pb_y * tg_zz_xxz_0[j] + fr * tg_zz_xxz_1[j];

                    tg_yzz_xyy_0[j] = pb_y * tg_zz_xyy_0[j] + fr * tg_zz_xyy_1[j] + fl1_fxn * tg_zz_xy_1[j];

                    tg_yzz_xyz_0[j] = pb_y * tg_zz_xyz_0[j] + fr * tg_zz_xyz_1[j] + 0.5 * fl1_fxn * tg_zz_xz_1[j];

                    tg_yzz_xzz_0[j] = pb_y * tg_zz_xzz_0[j] + fr * tg_zz_xzz_1[j];

                    tg_yzz_yyy_0[j] = pb_y * tg_zz_yyy_0[j] + fr * tg_zz_yyy_1[j] + 1.5 * fl1_fxn * tg_zz_yy_1[j];

                    tg_yzz_yyz_0[j] = pb_y * tg_zz_yyz_0[j] + fr * tg_zz_yyz_1[j] + fl1_fxn * tg_zz_yz_1[j];

                    tg_yzz_yzz_0[j] = pb_y * tg_zz_yzz_0[j] + fr * tg_zz_yzz_1[j] + 0.5 * fl1_fxn * tg_zz_zz_1[j];

                    tg_yzz_zzz_0[j] = pb_y * tg_zz_zzz_0[j] + fr * tg_zz_zzz_1[j];

                    fr = wp_z[j]; 

                    tg_zzz_xxx_0[j] = pb_z * tg_zz_xxx_0[j] + fr * tg_zz_xxx_1[j] + fl1_fx * (tg_z_xxx_0[j] - tg_z_xxx_1[j] * fl1_fza);

                    tg_zzz_xxy_0[j] = pb_z * tg_zz_xxy_0[j] + fr * tg_zz_xxy_1[j] + fl1_fx * (tg_z_xxy_0[j] - tg_z_xxy_1[j] * fl1_fza);

                    tg_zzz_xxz_0[j] = pb_z * tg_zz_xxz_0[j] + fr * tg_zz_xxz_1[j] + fl1_fx * (tg_z_xxz_0[j] - tg_z_xxz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zz_xx_1[j];

                    tg_zzz_xyy_0[j] = pb_z * tg_zz_xyy_0[j] + fr * tg_zz_xyy_1[j] + fl1_fx * (tg_z_xyy_0[j] - tg_z_xyy_1[j] * fl1_fza);

                    tg_zzz_xyz_0[j] = pb_z * tg_zz_xyz_0[j] + fr * tg_zz_xyz_1[j] + fl1_fx * (tg_z_xyz_0[j] - tg_z_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zz_xy_1[j];

                    tg_zzz_xzz_0[j] = pb_z * tg_zz_xzz_0[j] + fr * tg_zz_xzz_1[j] + fl1_fx * (tg_z_xzz_0[j] - tg_z_xzz_1[j] * fl1_fza) + fl1_fxn * tg_zz_xz_1[j];

                    tg_zzz_yyy_0[j] = pb_z * tg_zz_yyy_0[j] + fr * tg_zz_yyy_1[j] + fl1_fx * (tg_z_yyy_0[j] - tg_z_yyy_1[j] * fl1_fza);

                    tg_zzz_yyz_0[j] = pb_z * tg_zz_yyz_0[j] + fr * tg_zz_yyz_1[j] + fl1_fx * (tg_z_yyz_0[j] - tg_z_yyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zz_yy_1[j];

                    tg_zzz_yzz_0[j] = pb_z * tg_zz_yzz_0[j] + fr * tg_zz_yzz_1[j] + fl1_fx * (tg_z_yzz_0[j] - tg_z_yzz_1[j] * fl1_fza) + fl1_fxn * tg_zz_yz_1[j];

                    tg_zzz_zzz_0[j] = pb_z * tg_zz_zzz_0[j] + fr * tg_zz_zzz_1[j] + fl1_fx * (tg_z_zzz_0[j] - tg_z_zzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zz_zz_1[j];
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

