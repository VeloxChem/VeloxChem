//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ElectronRepulsionRecFuncForFG.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSFSG(      CMemBlock2D<double>* primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSFSG_0_75(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSFSG_75_150(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 
    }

    void
    compElectronRepulsionForSFSG_0_75(      CMemBlock2D<double>* primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& wpDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetPrimPairs,
                                      const int32_t              iContrPair)
    {
        // Batch of Integrals (0,75)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {3, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_3_4_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_3_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_2_4_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_2_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_1_4_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_1_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_2_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {3, -1, -1, -1}, 
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

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                // set up pointers to auxilary integrals

                auto tg_xx_xxxx_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx); 

                auto tg_xx_xxxy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 1); 

                auto tg_xx_xxxz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 2); 

                auto tg_xx_xxyy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 3); 

                auto tg_xx_xxyz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 4); 

                auto tg_xx_xxzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 5); 

                auto tg_xx_xyyy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 6); 

                auto tg_xx_xyyz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 7); 

                auto tg_xx_xyzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 8); 

                auto tg_xx_xzzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 9); 

                auto tg_xx_yyyy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 10); 

                auto tg_xx_yyyz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 11); 

                auto tg_xx_yyzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 12); 

                auto tg_xx_yzzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 13); 

                auto tg_xx_zzzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 14); 

                auto tg_xy_xxxx_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 15); 

                auto tg_xy_xxxy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 16); 

                auto tg_xy_xxxz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 17); 

                auto tg_xy_xxyy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 18); 

                auto tg_xy_xxyz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 19); 

                auto tg_xy_xxzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 20); 

                auto tg_xy_xyyy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 21); 

                auto tg_xy_xyyz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 22); 

                auto tg_xy_xyzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 23); 

                auto tg_xy_xzzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 24); 

                auto tg_xy_yyyy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 25); 

                auto tg_xy_yyyz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 26); 

                auto tg_xy_yyzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 27); 

                auto tg_xy_yzzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 28); 

                auto tg_xy_zzzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 29); 

                auto tg_xz_xxxx_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 30); 

                auto tg_xz_xxxy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 31); 

                auto tg_xz_xxxz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 32); 

                auto tg_xz_xxyy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 33); 

                auto tg_xz_xxyz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 34); 

                auto tg_xz_xxzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 35); 

                auto tg_xz_xyyy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 36); 

                auto tg_xz_xyyz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 37); 

                auto tg_xz_xyzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 38); 

                auto tg_xz_xzzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 39); 

                auto tg_xz_yyyy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 40); 

                auto tg_xz_yyyz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 41); 

                auto tg_xz_yyzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 42); 

                auto tg_xz_yzzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 43); 

                auto tg_xz_zzzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 44); 

                auto tg_yy_xxxx_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 45); 

                auto tg_yy_xxxy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 46); 

                auto tg_yy_xxxz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 47); 

                auto tg_yy_xxyy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 48); 

                auto tg_yy_xxyz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 49); 

                auto tg_yy_xxzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 50); 

                auto tg_yy_xyyy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 51); 

                auto tg_yy_xyyz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 52); 

                auto tg_yy_xyzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 53); 

                auto tg_yy_xzzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 54); 

                auto tg_yy_yyyy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 55); 

                auto tg_yy_yyyz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 56); 

                auto tg_yy_yyzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 57); 

                auto tg_yy_yzzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 58); 

                auto tg_yy_zzzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 59); 

                auto tg_yz_xxxx_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 60); 

                auto tg_yz_xxxy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 61); 

                auto tg_yz_xxxz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 62); 

                auto tg_yz_xxyy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 63); 

                auto tg_yz_xxyz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 64); 

                auto tg_yz_xxzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 65); 

                auto tg_yz_xyyy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 66); 

                auto tg_yz_xyyz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 67); 

                auto tg_yz_xyzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 68); 

                auto tg_yz_xzzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 69); 

                auto tg_yz_yyyy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 70); 

                auto tg_yz_yyyz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 71); 

                auto tg_yz_yyzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 72); 

                auto tg_yz_yzzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 73); 

                auto tg_yz_zzzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 74); 

                auto tg_xx_xxxx_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx); 

                auto tg_xx_xxxy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 1); 

                auto tg_xx_xxxz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 2); 

                auto tg_xx_xxyy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 3); 

                auto tg_xx_xxyz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 4); 

                auto tg_xx_xxzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 5); 

                auto tg_xx_xyyy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 6); 

                auto tg_xx_xyyz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 7); 

                auto tg_xx_xyzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 8); 

                auto tg_xx_xzzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 9); 

                auto tg_xx_yyyy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 10); 

                auto tg_xx_yyyz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 11); 

                auto tg_xx_yyzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 12); 

                auto tg_xx_yzzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 13); 

                auto tg_xx_zzzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 14); 

                auto tg_xy_xxxx_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 15); 

                auto tg_xy_xxxy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 16); 

                auto tg_xy_xxxz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 17); 

                auto tg_xy_xxyy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 18); 

                auto tg_xy_xxyz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 19); 

                auto tg_xy_xxzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 20); 

                auto tg_xy_xyyy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 21); 

                auto tg_xy_xyyz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 22); 

                auto tg_xy_xyzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 23); 

                auto tg_xy_xzzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 24); 

                auto tg_xy_yyyy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 25); 

                auto tg_xy_yyyz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 26); 

                auto tg_xy_yyzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 27); 

                auto tg_xy_yzzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 28); 

                auto tg_xy_zzzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 29); 

                auto tg_xz_xxxx_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 30); 

                auto tg_xz_xxxy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 31); 

                auto tg_xz_xxxz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 32); 

                auto tg_xz_xxyy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 33); 

                auto tg_xz_xxyz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 34); 

                auto tg_xz_xxzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 35); 

                auto tg_xz_xyyy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 36); 

                auto tg_xz_xyyz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 37); 

                auto tg_xz_xyzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 38); 

                auto tg_xz_xzzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 39); 

                auto tg_xz_yyyy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 40); 

                auto tg_xz_yyyz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 41); 

                auto tg_xz_yyzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 42); 

                auto tg_xz_yzzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 43); 

                auto tg_xz_zzzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 44); 

                auto tg_yy_xxxx_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 45); 

                auto tg_yy_xxxy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 46); 

                auto tg_yy_xxxz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 47); 

                auto tg_yy_xxyy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 48); 

                auto tg_yy_xxyz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 49); 

                auto tg_yy_xxzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 50); 

                auto tg_yy_xyyy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 51); 

                auto tg_yy_xyyz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 52); 

                auto tg_yy_xyzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 53); 

                auto tg_yy_xzzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 54); 

                auto tg_yy_yyyy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 55); 

                auto tg_yy_yyyz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 56); 

                auto tg_yy_yyzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 57); 

                auto tg_yy_yzzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 58); 

                auto tg_yy_zzzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 59); 

                auto tg_yz_xxxx_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 60); 

                auto tg_yz_xxxy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 61); 

                auto tg_yz_xxxz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 62); 

                auto tg_yz_xxyy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 63); 

                auto tg_yz_xxyz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 64); 

                auto tg_yz_xxzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 65); 

                auto tg_yz_xyyy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 66); 

                auto tg_yz_xyyz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 67); 

                auto tg_yz_xyzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 68); 

                auto tg_yz_xzzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 69); 

                auto tg_yz_yyyy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 70); 

                auto tg_yz_yyyz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 71); 

                auto tg_yz_yyzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 72); 

                auto tg_yz_yzzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 73); 

                auto tg_yz_zzzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 74); 

                auto tg_x_xxxx_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx); 

                auto tg_x_xxxy_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 1); 

                auto tg_x_xxxz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 2); 

                auto tg_x_xxyy_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 3); 

                auto tg_x_xxyz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 4); 

                auto tg_x_xxzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 5); 

                auto tg_x_xyyy_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 6); 

                auto tg_x_xyyz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 7); 

                auto tg_x_xyzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 8); 

                auto tg_x_xzzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 9); 

                auto tg_x_yyyy_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 10); 

                auto tg_x_yyyz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 11); 

                auto tg_x_yyzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 12); 

                auto tg_x_yzzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 13); 

                auto tg_x_zzzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 14); 

                auto tg_y_xxxx_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 15); 

                auto tg_y_xxxy_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 16); 

                auto tg_y_xxxz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 17); 

                auto tg_y_xxyy_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 18); 

                auto tg_y_xxyz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 19); 

                auto tg_y_xxzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 20); 

                auto tg_y_xyyy_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 21); 

                auto tg_y_xyyz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 22); 

                auto tg_y_xyzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 23); 

                auto tg_y_xzzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 24); 

                auto tg_y_yyyy_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 25); 

                auto tg_y_yyyz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 26); 

                auto tg_y_yyzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 27); 

                auto tg_y_yzzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 28); 

                auto tg_y_zzzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 29); 

                auto tg_z_xxxx_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 30); 

                auto tg_z_xxxy_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 31); 

                auto tg_z_xxxz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 32); 

                auto tg_z_xxyy_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 33); 

                auto tg_z_xxyz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 34); 

                auto tg_z_xxzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 35); 

                auto tg_z_xyyy_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 36); 

                auto tg_z_xyyz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 37); 

                auto tg_z_xyzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 38); 

                auto tg_z_xzzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 39); 

                auto tg_z_yyyy_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 40); 

                auto tg_z_yyyz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 41); 

                auto tg_z_yyzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 42); 

                auto tg_z_yzzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 43); 

                auto tg_z_zzzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 44); 

                auto tg_x_xxxx_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx); 

                auto tg_x_xxxy_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 1); 

                auto tg_x_xxxz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 2); 

                auto tg_x_xxyy_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 3); 

                auto tg_x_xxyz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 4); 

                auto tg_x_xxzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 5); 

                auto tg_x_xyyy_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 6); 

                auto tg_x_xyyz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 7); 

                auto tg_x_xyzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 8); 

                auto tg_x_xzzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 9); 

                auto tg_x_yyyy_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 10); 

                auto tg_x_yyyz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 11); 

                auto tg_x_yyzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 12); 

                auto tg_x_yzzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 13); 

                auto tg_x_zzzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 14); 

                auto tg_y_xxxx_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 15); 

                auto tg_y_xxxy_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 16); 

                auto tg_y_xxxz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 17); 

                auto tg_y_xxyy_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 18); 

                auto tg_y_xxyz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 19); 

                auto tg_y_xxzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 20); 

                auto tg_y_xyyy_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 21); 

                auto tg_y_xyyz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 22); 

                auto tg_y_xyzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 23); 

                auto tg_y_xzzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 24); 

                auto tg_y_yyyy_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 25); 

                auto tg_y_yyyz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 26); 

                auto tg_y_yyzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 27); 

                auto tg_y_yzzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 28); 

                auto tg_y_zzzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 29); 

                auto tg_z_xxxx_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 30); 

                auto tg_z_xxxy_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 31); 

                auto tg_z_xxxz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 32); 

                auto tg_z_xxyy_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 33); 

                auto tg_z_xxyz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 34); 

                auto tg_z_xxzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 35); 

                auto tg_z_xyyy_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 36); 

                auto tg_z_xyyz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 37); 

                auto tg_z_xyzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 38); 

                auto tg_z_xzzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 39); 

                auto tg_z_yyyy_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 40); 

                auto tg_z_yyyz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 41); 

                auto tg_z_yyzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 42); 

                auto tg_z_yzzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 43); 

                auto tg_z_zzzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 44); 

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

                // set up pointers to integrals

                auto tg_xxx_xxxx_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx); 

                auto tg_xxx_xxxy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 1); 

                auto tg_xxx_xxxz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 2); 

                auto tg_xxx_xxyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 3); 

                auto tg_xxx_xxyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 4); 

                auto tg_xxx_xxzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 5); 

                auto tg_xxx_xyyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 6); 

                auto tg_xxx_xyyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 7); 

                auto tg_xxx_xyzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 8); 

                auto tg_xxx_xzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 9); 

                auto tg_xxx_yyyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 10); 

                auto tg_xxx_yyyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 11); 

                auto tg_xxx_yyzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 12); 

                auto tg_xxx_yzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 13); 

                auto tg_xxx_zzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 14); 

                auto tg_xxy_xxxx_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 15); 

                auto tg_xxy_xxxy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 16); 

                auto tg_xxy_xxxz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 17); 

                auto tg_xxy_xxyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 18); 

                auto tg_xxy_xxyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 19); 

                auto tg_xxy_xxzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 20); 

                auto tg_xxy_xyyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 21); 

                auto tg_xxy_xyyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 22); 

                auto tg_xxy_xyzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 23); 

                auto tg_xxy_xzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 24); 

                auto tg_xxy_yyyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 25); 

                auto tg_xxy_yyyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 26); 

                auto tg_xxy_yyzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 27); 

                auto tg_xxy_yzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 28); 

                auto tg_xxy_zzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 29); 

                auto tg_xxz_xxxx_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 30); 

                auto tg_xxz_xxxy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 31); 

                auto tg_xxz_xxxz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 32); 

                auto tg_xxz_xxyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 33); 

                auto tg_xxz_xxyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 34); 

                auto tg_xxz_xxzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 35); 

                auto tg_xxz_xyyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 36); 

                auto tg_xxz_xyyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 37); 

                auto tg_xxz_xyzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 38); 

                auto tg_xxz_xzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 39); 

                auto tg_xxz_yyyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 40); 

                auto tg_xxz_yyyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 41); 

                auto tg_xxz_yyzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 42); 

                auto tg_xxz_yzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 43); 

                auto tg_xxz_zzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 44); 

                auto tg_xyy_xxxx_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 45); 

                auto tg_xyy_xxxy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 46); 

                auto tg_xyy_xxxz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 47); 

                auto tg_xyy_xxyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 48); 

                auto tg_xyy_xxyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 49); 

                auto tg_xyy_xxzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 50); 

                auto tg_xyy_xyyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 51); 

                auto tg_xyy_xyyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 52); 

                auto tg_xyy_xyzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 53); 

                auto tg_xyy_xzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 54); 

                auto tg_xyy_yyyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 55); 

                auto tg_xyy_yyyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 56); 

                auto tg_xyy_yyzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 57); 

                auto tg_xyy_yzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 58); 

                auto tg_xyy_zzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 59); 

                auto tg_xyz_xxxx_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 60); 

                auto tg_xyz_xxxy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 61); 

                auto tg_xyz_xxxz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 62); 

                auto tg_xyz_xxyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 63); 

                auto tg_xyz_xxyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 64); 

                auto tg_xyz_xxzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 65); 

                auto tg_xyz_xyyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 66); 

                auto tg_xyz_xyyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 67); 

                auto tg_xyz_xyzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 68); 

                auto tg_xyz_xzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 69); 

                auto tg_xyz_yyyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 70); 

                auto tg_xyz_yyyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 71); 

                auto tg_xyz_yyzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 72); 

                auto tg_xyz_yzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 73); 

                auto tg_xyz_zzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 74); 

                // Batch of Integrals (0,75)

                #pragma omp simd aligned(fxn, fza, tg_x_xxxx_0, tg_x_xxxx_1, tg_x_xxxy_0, tg_x_xxxy_1, \
                                         tg_x_xxxz_0, tg_x_xxxz_1, tg_x_xxyy_0, tg_x_xxyy_1, tg_x_xxyz_0, tg_x_xxyz_1, \
                                         tg_x_xxzz_0, tg_x_xxzz_1, tg_x_xyyy_0, tg_x_xyyy_1, tg_x_xyyz_0, tg_x_xyyz_1, \
                                         tg_x_xyzz_0, tg_x_xyzz_1, tg_x_xzzz_0, tg_x_xzzz_1, tg_x_yyyy_0, tg_x_yyyy_1, \
                                         tg_x_yyyz_0, tg_x_yyyz_1, tg_x_yyzz_0, tg_x_yyzz_1, tg_x_yzzz_0, tg_x_yzzz_1, \
                                         tg_x_zzzz_0, tg_x_zzzz_1, tg_xx_xxx_1, tg_xx_xxxx_0, tg_xx_xxxx_1, tg_xx_xxxy_0, \
                                         tg_xx_xxxy_1, tg_xx_xxxz_0, tg_xx_xxxz_1, tg_xx_xxy_1, tg_xx_xxyy_0, tg_xx_xxyy_1, \
                                         tg_xx_xxyz_0, tg_xx_xxyz_1, tg_xx_xxz_1, tg_xx_xxzz_0, tg_xx_xxzz_1, tg_xx_xyy_1, \
                                         tg_xx_xyyy_0, tg_xx_xyyy_1, tg_xx_xyyz_0, tg_xx_xyyz_1, tg_xx_xyz_1, tg_xx_xyzz_0, \
                                         tg_xx_xyzz_1, tg_xx_xzz_1, tg_xx_xzzz_0, tg_xx_xzzz_1, tg_xx_yyy_1, tg_xx_yyyy_0, \
                                         tg_xx_yyyy_1, tg_xx_yyyz_0, tg_xx_yyyz_1, tg_xx_yyz_1, tg_xx_yyzz_0, tg_xx_yyzz_1, \
                                         tg_xx_yzz_1, tg_xx_yzzz_0, tg_xx_yzzz_1, tg_xx_zzz_1, tg_xx_zzzz_0, tg_xx_zzzz_1, \
                                         tg_xxx_xxxx_0, tg_xxx_xxxy_0, tg_xxx_xxxz_0, tg_xxx_xxyy_0, tg_xxx_xxyz_0, \
                                         tg_xxx_xxzz_0, tg_xxx_xyyy_0, tg_xxx_xyyz_0, tg_xxx_xyzz_0, tg_xxx_xzzz_0, \
                                         tg_xxx_yyyy_0, tg_xxx_yyyz_0, tg_xxx_yyzz_0, tg_xxx_yzzz_0, tg_xxx_zzzz_0, \
                                         tg_xxy_xxxx_0, tg_xxy_xxxy_0, tg_xxy_xxxz_0, tg_xxy_xxyy_0, tg_xxy_xxyz_0, \
                                         tg_xxy_xxzz_0, tg_xxy_xyyy_0, tg_xxy_xyyz_0, tg_xxy_xyzz_0, tg_xxy_xzzz_0, \
                                         tg_xxy_yyyy_0, tg_xxy_yyyz_0, tg_xxy_yyzz_0, tg_xxy_yzzz_0, tg_xxy_zzzz_0, \
                                         tg_xxz_xxxx_0, tg_xxz_xxxy_0, tg_xxz_xxxz_0, tg_xxz_xxyy_0, tg_xxz_xxyz_0, \
                                         tg_xxz_xxzz_0, tg_xxz_xyyy_0, tg_xxz_xyyz_0, tg_xxz_xyzz_0, tg_xxz_xzzz_0, \
                                         tg_xxz_yyyy_0, tg_xxz_yyyz_0, tg_xxz_yyzz_0, tg_xxz_yzzz_0, tg_xxz_zzzz_0, \
                                         tg_xy_xxx_1, tg_xy_xxxx_0, tg_xy_xxxx_1, tg_xy_xxxy_0, tg_xy_xxxy_1, tg_xy_xxxz_0, \
                                         tg_xy_xxxz_1, tg_xy_xxy_1, tg_xy_xxyy_0, tg_xy_xxyy_1, tg_xy_xxyz_0, tg_xy_xxyz_1, \
                                         tg_xy_xxz_1, tg_xy_xxzz_0, tg_xy_xxzz_1, tg_xy_xyy_1, tg_xy_xyyy_0, tg_xy_xyyy_1, \
                                         tg_xy_xyyz_0, tg_xy_xyyz_1, tg_xy_xyz_1, tg_xy_xyzz_0, tg_xy_xyzz_1, tg_xy_xzz_1, \
                                         tg_xy_xzzz_0, tg_xy_xzzz_1, tg_xy_yyy_1, tg_xy_yyyy_0, tg_xy_yyyy_1, tg_xy_yyyz_0, \
                                         tg_xy_yyyz_1, tg_xy_yyz_1, tg_xy_yyzz_0, tg_xy_yyzz_1, tg_xy_yzz_1, tg_xy_yzzz_0, \
                                         tg_xy_yzzz_1, tg_xy_zzz_1, tg_xy_zzzz_0, tg_xy_zzzz_1, tg_xyy_xxxx_0, \
                                         tg_xyy_xxxy_0, tg_xyy_xxxz_0, tg_xyy_xxyy_0, tg_xyy_xxyz_0, tg_xyy_xxzz_0, \
                                         tg_xyy_xyyy_0, tg_xyy_xyyz_0, tg_xyy_xyzz_0, tg_xyy_xzzz_0, tg_xyy_yyyy_0, \
                                         tg_xyy_yyyz_0, tg_xyy_yyzz_0, tg_xyy_yzzz_0, tg_xyy_zzzz_0, tg_xyz_xxxx_0, \
                                         tg_xyz_xxxy_0, tg_xyz_xxxz_0, tg_xyz_xxyy_0, tg_xyz_xxyz_0, tg_xyz_xxzz_0, \
                                         tg_xyz_xyyy_0, tg_xyz_xyyz_0, tg_xyz_xyzz_0, tg_xyz_xzzz_0, tg_xyz_yyyy_0, \
                                         tg_xyz_yyyz_0, tg_xyz_yyzz_0, tg_xyz_yzzz_0, tg_xyz_zzzz_0, tg_xz_xxx_1, \
                                         tg_xz_xxxx_0, tg_xz_xxxx_1, tg_xz_xxxy_0, tg_xz_xxxy_1, tg_xz_xxxz_0, tg_xz_xxxz_1, \
                                         tg_xz_xxy_1, tg_xz_xxyy_0, tg_xz_xxyy_1, tg_xz_xxyz_0, tg_xz_xxyz_1, tg_xz_xxz_1, \
                                         tg_xz_xxzz_0, tg_xz_xxzz_1, tg_xz_xyy_1, tg_xz_xyyy_0, tg_xz_xyyy_1, tg_xz_xyyz_0, \
                                         tg_xz_xyyz_1, tg_xz_xyz_1, tg_xz_xyzz_0, tg_xz_xyzz_1, tg_xz_xzz_1, tg_xz_xzzz_0, \
                                         tg_xz_xzzz_1, tg_xz_yyy_1, tg_xz_yyyy_0, tg_xz_yyyy_1, tg_xz_yyyz_0, tg_xz_yyyz_1, \
                                         tg_xz_yyz_1, tg_xz_yyzz_0, tg_xz_yyzz_1, tg_xz_yzz_1, tg_xz_yzzz_0, tg_xz_yzzz_1, \
                                         tg_xz_zzz_1, tg_xz_zzzz_0, tg_xz_zzzz_1, tg_y_xxxx_0, tg_y_xxxx_1, tg_y_xxxy_0, \
                                         tg_y_xxxy_1, tg_y_xxxz_0, tg_y_xxxz_1, tg_y_xxyy_0, tg_y_xxyy_1, tg_y_xxyz_0, \
                                         tg_y_xxyz_1, tg_y_xxzz_0, tg_y_xxzz_1, tg_y_xyyy_0, tg_y_xyyy_1, tg_y_xyyz_0, \
                                         tg_y_xyyz_1, tg_y_xyzz_0, tg_y_xyzz_1, tg_y_xzzz_0, tg_y_xzzz_1, tg_y_yyyy_0, \
                                         tg_y_yyyy_1, tg_y_yyyz_0, tg_y_yyyz_1, tg_y_yyzz_0, tg_y_yyzz_1, tg_y_yzzz_0, \
                                         tg_y_yzzz_1, tg_y_zzzz_0, tg_y_zzzz_1, tg_yy_xxx_1, tg_yy_xxxx_0, tg_yy_xxxx_1, \
                                         tg_yy_xxxy_0, tg_yy_xxxy_1, tg_yy_xxxz_0, tg_yy_xxxz_1, tg_yy_xxy_1, tg_yy_xxyy_0, \
                                         tg_yy_xxyy_1, tg_yy_xxyz_0, tg_yy_xxyz_1, tg_yy_xxz_1, tg_yy_xxzz_0, tg_yy_xxzz_1, \
                                         tg_yy_xyy_1, tg_yy_xyyy_0, tg_yy_xyyy_1, tg_yy_xyyz_0, tg_yy_xyyz_1, tg_yy_xyz_1, \
                                         tg_yy_xyzz_0, tg_yy_xyzz_1, tg_yy_xzz_1, tg_yy_xzzz_0, tg_yy_xzzz_1, tg_yy_yyy_1, \
                                         tg_yy_yyyy_0, tg_yy_yyyy_1, tg_yy_yyyz_0, tg_yy_yyyz_1, tg_yy_yyz_1, tg_yy_yyzz_0, \
                                         tg_yy_yyzz_1, tg_yy_yzz_1, tg_yy_yzzz_0, tg_yy_yzzz_1, tg_yy_zzz_1, tg_yy_zzzz_0, \
                                         tg_yy_zzzz_1, tg_yz_xxx_1, tg_yz_xxxx_0, tg_yz_xxxx_1, tg_yz_xxxy_0, tg_yz_xxxy_1, \
                                         tg_yz_xxxz_0, tg_yz_xxxz_1, tg_yz_xxy_1, tg_yz_xxyy_0, tg_yz_xxyy_1, tg_yz_xxyz_0, \
                                         tg_yz_xxyz_1, tg_yz_xxz_1, tg_yz_xxzz_0, tg_yz_xxzz_1, tg_yz_xyy_1, tg_yz_xyyy_0, \
                                         tg_yz_xyyy_1, tg_yz_xyyz_0, tg_yz_xyyz_1, tg_yz_xyz_1, tg_yz_xyzz_0, tg_yz_xyzz_1, \
                                         tg_yz_xzz_1, tg_yz_xzzz_0, tg_yz_xzzz_1, tg_yz_yyy_1, tg_yz_yyyy_0, tg_yz_yyyy_1, \
                                         tg_yz_yyyz_0, tg_yz_yyyz_1, tg_yz_yyz_1, tg_yz_yyzz_0, tg_yz_yyzz_1, tg_yz_yzz_1, \
                                         tg_yz_yzzz_0, tg_yz_yzzz_1, tg_yz_zzz_1, tg_yz_zzzz_0, tg_yz_zzzz_1, tg_z_xxxx_0, \
                                         tg_z_xxxx_1, tg_z_xxxy_0, tg_z_xxxy_1, tg_z_xxxz_0, tg_z_xxxz_1, tg_z_xxyy_0, \
                                         tg_z_xxyy_1, tg_z_xxyz_0, tg_z_xxyz_1, tg_z_xxzz_0, tg_z_xxzz_1, tg_z_xyyy_0, \
                                         tg_z_xyyy_1, tg_z_xyyz_0, tg_z_xyyz_1, tg_z_xyzz_0, tg_z_xyzz_1, tg_z_xzzz_0, \
                                         tg_z_xzzz_1, tg_z_yyyy_0, tg_z_yyyy_1, tg_z_yyyz_0, tg_z_yyyz_1, tg_z_yyzz_0, \
                                         tg_z_yyzz_1, tg_z_yzzz_0, tg_z_yzzz_1, tg_z_zzzz_0, tg_z_zzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxx_xxxx_0[j] = pb_x * tg_xx_xxxx_0[j] + fr * tg_xx_xxxx_1[j] + fl1_fx * (tg_x_xxxx_0[j] - tg_x_xxxx_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xx_xxx_1[j];

                    tg_xxx_xxxy_0[j] = pb_x * tg_xx_xxxy_0[j] + fr * tg_xx_xxxy_1[j] + fl1_fx * (tg_x_xxxy_0[j] - tg_x_xxxy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xx_xxy_1[j];

                    tg_xxx_xxxz_0[j] = pb_x * tg_xx_xxxz_0[j] + fr * tg_xx_xxxz_1[j] + fl1_fx * (tg_x_xxxz_0[j] - tg_x_xxxz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xx_xxz_1[j];

                    tg_xxx_xxyy_0[j] = pb_x * tg_xx_xxyy_0[j] + fr * tg_xx_xxyy_1[j] + fl1_fx * (tg_x_xxyy_0[j] - tg_x_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_xx_xyy_1[j];

                    tg_xxx_xxyz_0[j] = pb_x * tg_xx_xxyz_0[j] + fr * tg_xx_xxyz_1[j] + fl1_fx * (tg_x_xxyz_0[j] - tg_x_xxyz_1[j] * fl1_fza) + fl1_fxn * tg_xx_xyz_1[j];

                    tg_xxx_xxzz_0[j] = pb_x * tg_xx_xxzz_0[j] + fr * tg_xx_xxzz_1[j] + fl1_fx * (tg_x_xxzz_0[j] - tg_x_xxzz_1[j] * fl1_fza) + fl1_fxn * tg_xx_xzz_1[j];

                    tg_xxx_xyyy_0[j] = pb_x * tg_xx_xyyy_0[j] + fr * tg_xx_xyyy_1[j] + fl1_fx * (tg_x_xyyy_0[j] - tg_x_xyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xx_yyy_1[j];

                    tg_xxx_xyyz_0[j] = pb_x * tg_xx_xyyz_0[j] + fr * tg_xx_xyyz_1[j] + fl1_fx * (tg_x_xyyz_0[j] - tg_x_xyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xx_yyz_1[j];

                    tg_xxx_xyzz_0[j] = pb_x * tg_xx_xyzz_0[j] + fr * tg_xx_xyzz_1[j] + fl1_fx * (tg_x_xyzz_0[j] - tg_x_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xx_yzz_1[j];

                    tg_xxx_xzzz_0[j] = pb_x * tg_xx_xzzz_0[j] + fr * tg_xx_xzzz_1[j] + fl1_fx * (tg_x_xzzz_0[j] - tg_x_xzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xx_zzz_1[j];

                    tg_xxx_yyyy_0[j] = pb_x * tg_xx_yyyy_0[j] + fr * tg_xx_yyyy_1[j] + fl1_fx * (tg_x_yyyy_0[j] - tg_x_yyyy_1[j] * fl1_fza);

                    tg_xxx_yyyz_0[j] = pb_x * tg_xx_yyyz_0[j] + fr * tg_xx_yyyz_1[j] + fl1_fx * (tg_x_yyyz_0[j] - tg_x_yyyz_1[j] * fl1_fza);

                    tg_xxx_yyzz_0[j] = pb_x * tg_xx_yyzz_0[j] + fr * tg_xx_yyzz_1[j] + fl1_fx * (tg_x_yyzz_0[j] - tg_x_yyzz_1[j] * fl1_fza);

                    tg_xxx_yzzz_0[j] = pb_x * tg_xx_yzzz_0[j] + fr * tg_xx_yzzz_1[j] + fl1_fx * (tg_x_yzzz_0[j] - tg_x_yzzz_1[j] * fl1_fza);

                    tg_xxx_zzzz_0[j] = pb_x * tg_xx_zzzz_0[j] + fr * tg_xx_zzzz_1[j] + fl1_fx * (tg_x_zzzz_0[j] - tg_x_zzzz_1[j] * fl1_fza);

                    tg_xxy_xxxx_0[j] = pb_x * tg_xy_xxxx_0[j] + fr * tg_xy_xxxx_1[j] + 0.5 * fl1_fx * (tg_y_xxxx_0[j] - tg_y_xxxx_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xy_xxx_1[j];

                    tg_xxy_xxxy_0[j] = pb_x * tg_xy_xxxy_0[j] + fr * tg_xy_xxxy_1[j] + 0.5 * fl1_fx * (tg_y_xxxy_0[j] - tg_y_xxxy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xy_xxy_1[j];

                    tg_xxy_xxxz_0[j] = pb_x * tg_xy_xxxz_0[j] + fr * tg_xy_xxxz_1[j] + 0.5 * fl1_fx * (tg_y_xxxz_0[j] - tg_y_xxxz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xy_xxz_1[j];

                    tg_xxy_xxyy_0[j] = pb_x * tg_xy_xxyy_0[j] + fr * tg_xy_xxyy_1[j] + 0.5 * fl1_fx * (tg_y_xxyy_0[j] - tg_y_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_xy_xyy_1[j];

                    tg_xxy_xxyz_0[j] = pb_x * tg_xy_xxyz_0[j] + fr * tg_xy_xxyz_1[j] + 0.5 * fl1_fx * (tg_y_xxyz_0[j] - tg_y_xxyz_1[j] * fl1_fza) + fl1_fxn * tg_xy_xyz_1[j];

                    tg_xxy_xxzz_0[j] = pb_x * tg_xy_xxzz_0[j] + fr * tg_xy_xxzz_1[j] + 0.5 * fl1_fx * (tg_y_xxzz_0[j] - tg_y_xxzz_1[j] * fl1_fza) + fl1_fxn * tg_xy_xzz_1[j];

                    tg_xxy_xyyy_0[j] = pb_x * tg_xy_xyyy_0[j] + fr * tg_xy_xyyy_1[j] + 0.5 * fl1_fx * (tg_y_xyyy_0[j] - tg_y_xyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xy_yyy_1[j];

                    tg_xxy_xyyz_0[j] = pb_x * tg_xy_xyyz_0[j] + fr * tg_xy_xyyz_1[j] + 0.5 * fl1_fx * (tg_y_xyyz_0[j] - tg_y_xyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xy_yyz_1[j];

                    tg_xxy_xyzz_0[j] = pb_x * tg_xy_xyzz_0[j] + fr * tg_xy_xyzz_1[j] + 0.5 * fl1_fx * (tg_y_xyzz_0[j] - tg_y_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xy_yzz_1[j];

                    tg_xxy_xzzz_0[j] = pb_x * tg_xy_xzzz_0[j] + fr * tg_xy_xzzz_1[j] + 0.5 * fl1_fx * (tg_y_xzzz_0[j] - tg_y_xzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xy_zzz_1[j];

                    tg_xxy_yyyy_0[j] = pb_x * tg_xy_yyyy_0[j] + fr * tg_xy_yyyy_1[j] + 0.5 * fl1_fx * (tg_y_yyyy_0[j] - tg_y_yyyy_1[j] * fl1_fza);

                    tg_xxy_yyyz_0[j] = pb_x * tg_xy_yyyz_0[j] + fr * tg_xy_yyyz_1[j] + 0.5 * fl1_fx * (tg_y_yyyz_0[j] - tg_y_yyyz_1[j] * fl1_fza);

                    tg_xxy_yyzz_0[j] = pb_x * tg_xy_yyzz_0[j] + fr * tg_xy_yyzz_1[j] + 0.5 * fl1_fx * (tg_y_yyzz_0[j] - tg_y_yyzz_1[j] * fl1_fza);

                    tg_xxy_yzzz_0[j] = pb_x * tg_xy_yzzz_0[j] + fr * tg_xy_yzzz_1[j] + 0.5 * fl1_fx * (tg_y_yzzz_0[j] - tg_y_yzzz_1[j] * fl1_fza);

                    tg_xxy_zzzz_0[j] = pb_x * tg_xy_zzzz_0[j] + fr * tg_xy_zzzz_1[j] + 0.5 * fl1_fx * (tg_y_zzzz_0[j] - tg_y_zzzz_1[j] * fl1_fza);

                    tg_xxz_xxxx_0[j] = pb_x * tg_xz_xxxx_0[j] + fr * tg_xz_xxxx_1[j] + 0.5 * fl1_fx * (tg_z_xxxx_0[j] - tg_z_xxxx_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xz_xxx_1[j];

                    tg_xxz_xxxy_0[j] = pb_x * tg_xz_xxxy_0[j] + fr * tg_xz_xxxy_1[j] + 0.5 * fl1_fx * (tg_z_xxxy_0[j] - tg_z_xxxy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xz_xxy_1[j];

                    tg_xxz_xxxz_0[j] = pb_x * tg_xz_xxxz_0[j] + fr * tg_xz_xxxz_1[j] + 0.5 * fl1_fx * (tg_z_xxxz_0[j] - tg_z_xxxz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xz_xxz_1[j];

                    tg_xxz_xxyy_0[j] = pb_x * tg_xz_xxyy_0[j] + fr * tg_xz_xxyy_1[j] + 0.5 * fl1_fx * (tg_z_xxyy_0[j] - tg_z_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_xz_xyy_1[j];

                    tg_xxz_xxyz_0[j] = pb_x * tg_xz_xxyz_0[j] + fr * tg_xz_xxyz_1[j] + 0.5 * fl1_fx * (tg_z_xxyz_0[j] - tg_z_xxyz_1[j] * fl1_fza) + fl1_fxn * tg_xz_xyz_1[j];

                    tg_xxz_xxzz_0[j] = pb_x * tg_xz_xxzz_0[j] + fr * tg_xz_xxzz_1[j] + 0.5 * fl1_fx * (tg_z_xxzz_0[j] - tg_z_xxzz_1[j] * fl1_fza) + fl1_fxn * tg_xz_xzz_1[j];

                    tg_xxz_xyyy_0[j] = pb_x * tg_xz_xyyy_0[j] + fr * tg_xz_xyyy_1[j] + 0.5 * fl1_fx * (tg_z_xyyy_0[j] - tg_z_xyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xz_yyy_1[j];

                    tg_xxz_xyyz_0[j] = pb_x * tg_xz_xyyz_0[j] + fr * tg_xz_xyyz_1[j] + 0.5 * fl1_fx * (tg_z_xyyz_0[j] - tg_z_xyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xz_yyz_1[j];

                    tg_xxz_xyzz_0[j] = pb_x * tg_xz_xyzz_0[j] + fr * tg_xz_xyzz_1[j] + 0.5 * fl1_fx * (tg_z_xyzz_0[j] - tg_z_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xz_yzz_1[j];

                    tg_xxz_xzzz_0[j] = pb_x * tg_xz_xzzz_0[j] + fr * tg_xz_xzzz_1[j] + 0.5 * fl1_fx * (tg_z_xzzz_0[j] - tg_z_xzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xz_zzz_1[j];

                    tg_xxz_yyyy_0[j] = pb_x * tg_xz_yyyy_0[j] + fr * tg_xz_yyyy_1[j] + 0.5 * fl1_fx * (tg_z_yyyy_0[j] - tg_z_yyyy_1[j] * fl1_fza);

                    tg_xxz_yyyz_0[j] = pb_x * tg_xz_yyyz_0[j] + fr * tg_xz_yyyz_1[j] + 0.5 * fl1_fx * (tg_z_yyyz_0[j] - tg_z_yyyz_1[j] * fl1_fza);

                    tg_xxz_yyzz_0[j] = pb_x * tg_xz_yyzz_0[j] + fr * tg_xz_yyzz_1[j] + 0.5 * fl1_fx * (tg_z_yyzz_0[j] - tg_z_yyzz_1[j] * fl1_fza);

                    tg_xxz_yzzz_0[j] = pb_x * tg_xz_yzzz_0[j] + fr * tg_xz_yzzz_1[j] + 0.5 * fl1_fx * (tg_z_yzzz_0[j] - tg_z_yzzz_1[j] * fl1_fza);

                    tg_xxz_zzzz_0[j] = pb_x * tg_xz_zzzz_0[j] + fr * tg_xz_zzzz_1[j] + 0.5 * fl1_fx * (tg_z_zzzz_0[j] - tg_z_zzzz_1[j] * fl1_fza);

                    tg_xyy_xxxx_0[j] = pb_x * tg_yy_xxxx_0[j] + fr * tg_yy_xxxx_1[j] + 2.0 * fl1_fxn * tg_yy_xxx_1[j];

                    tg_xyy_xxxy_0[j] = pb_x * tg_yy_xxxy_0[j] + fr * tg_yy_xxxy_1[j] + 1.5 * fl1_fxn * tg_yy_xxy_1[j];

                    tg_xyy_xxxz_0[j] = pb_x * tg_yy_xxxz_0[j] + fr * tg_yy_xxxz_1[j] + 1.5 * fl1_fxn * tg_yy_xxz_1[j];

                    tg_xyy_xxyy_0[j] = pb_x * tg_yy_xxyy_0[j] + fr * tg_yy_xxyy_1[j] + fl1_fxn * tg_yy_xyy_1[j];

                    tg_xyy_xxyz_0[j] = pb_x * tg_yy_xxyz_0[j] + fr * tg_yy_xxyz_1[j] + fl1_fxn * tg_yy_xyz_1[j];

                    tg_xyy_xxzz_0[j] = pb_x * tg_yy_xxzz_0[j] + fr * tg_yy_xxzz_1[j] + fl1_fxn * tg_yy_xzz_1[j];

                    tg_xyy_xyyy_0[j] = pb_x * tg_yy_xyyy_0[j] + fr * tg_yy_xyyy_1[j] + 0.5 * fl1_fxn * tg_yy_yyy_1[j];

                    tg_xyy_xyyz_0[j] = pb_x * tg_yy_xyyz_0[j] + fr * tg_yy_xyyz_1[j] + 0.5 * fl1_fxn * tg_yy_yyz_1[j];

                    tg_xyy_xyzz_0[j] = pb_x * tg_yy_xyzz_0[j] + fr * tg_yy_xyzz_1[j] + 0.5 * fl1_fxn * tg_yy_yzz_1[j];

                    tg_xyy_xzzz_0[j] = pb_x * tg_yy_xzzz_0[j] + fr * tg_yy_xzzz_1[j] + 0.5 * fl1_fxn * tg_yy_zzz_1[j];

                    tg_xyy_yyyy_0[j] = pb_x * tg_yy_yyyy_0[j] + fr * tg_yy_yyyy_1[j];

                    tg_xyy_yyyz_0[j] = pb_x * tg_yy_yyyz_0[j] + fr * tg_yy_yyyz_1[j];

                    tg_xyy_yyzz_0[j] = pb_x * tg_yy_yyzz_0[j] + fr * tg_yy_yyzz_1[j];

                    tg_xyy_yzzz_0[j] = pb_x * tg_yy_yzzz_0[j] + fr * tg_yy_yzzz_1[j];

                    tg_xyy_zzzz_0[j] = pb_x * tg_yy_zzzz_0[j] + fr * tg_yy_zzzz_1[j];

                    tg_xyz_xxxx_0[j] = pb_x * tg_yz_xxxx_0[j] + fr * tg_yz_xxxx_1[j] + 2.0 * fl1_fxn * tg_yz_xxx_1[j];

                    tg_xyz_xxxy_0[j] = pb_x * tg_yz_xxxy_0[j] + fr * tg_yz_xxxy_1[j] + 1.5 * fl1_fxn * tg_yz_xxy_1[j];

                    tg_xyz_xxxz_0[j] = pb_x * tg_yz_xxxz_0[j] + fr * tg_yz_xxxz_1[j] + 1.5 * fl1_fxn * tg_yz_xxz_1[j];

                    tg_xyz_xxyy_0[j] = pb_x * tg_yz_xxyy_0[j] + fr * tg_yz_xxyy_1[j] + fl1_fxn * tg_yz_xyy_1[j];

                    tg_xyz_xxyz_0[j] = pb_x * tg_yz_xxyz_0[j] + fr * tg_yz_xxyz_1[j] + fl1_fxn * tg_yz_xyz_1[j];

                    tg_xyz_xxzz_0[j] = pb_x * tg_yz_xxzz_0[j] + fr * tg_yz_xxzz_1[j] + fl1_fxn * tg_yz_xzz_1[j];

                    tg_xyz_xyyy_0[j] = pb_x * tg_yz_xyyy_0[j] + fr * tg_yz_xyyy_1[j] + 0.5 * fl1_fxn * tg_yz_yyy_1[j];

                    tg_xyz_xyyz_0[j] = pb_x * tg_yz_xyyz_0[j] + fr * tg_yz_xyyz_1[j] + 0.5 * fl1_fxn * tg_yz_yyz_1[j];

                    tg_xyz_xyzz_0[j] = pb_x * tg_yz_xyzz_0[j] + fr * tg_yz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yz_yzz_1[j];

                    tg_xyz_xzzz_0[j] = pb_x * tg_yz_xzzz_0[j] + fr * tg_yz_xzzz_1[j] + 0.5 * fl1_fxn * tg_yz_zzz_1[j];

                    tg_xyz_yyyy_0[j] = pb_x * tg_yz_yyyy_0[j] + fr * tg_yz_yyyy_1[j];

                    tg_xyz_yyyz_0[j] = pb_x * tg_yz_yyyz_0[j] + fr * tg_yz_yyyz_1[j];

                    tg_xyz_yyzz_0[j] = pb_x * tg_yz_yyzz_0[j] + fr * tg_yz_yyzz_1[j];

                    tg_xyz_yzzz_0[j] = pb_x * tg_yz_yzzz_0[j] + fr * tg_yz_yzzz_1[j];

                    tg_xyz_zzzz_0[j] = pb_x * tg_yz_zzzz_0[j] + fr * tg_yz_zzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSFSG_75_150(      CMemBlock2D<double>* primBuffer,
                                        const CRecursionMap&       recursionMap,
                                        const CMemBlock2D<double>& osFactors,
                                        const CMemBlock2D<double>& wpDistances,
                                        const CGtoPairsBlock&      braGtoPairsBlock,
                                        const CGtoPairsBlock&      ketGtoPairsBlock,
                                        const int32_t              nKetPrimPairs,
                                        const int32_t              iContrPair)
    {
        // Batch of Integrals (75,150)

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
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_3_4_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_3_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_2_4_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_2_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_2_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_1_4_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_1_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {1, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tg_yy_xxxx_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 45); 

                auto tg_yy_xxxy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 46); 

                auto tg_yy_xxxz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 47); 

                auto tg_yy_xxyy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 48); 

                auto tg_yy_xxyz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 49); 

                auto tg_yy_xxzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 50); 

                auto tg_yy_xyyy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 51); 

                auto tg_yy_xyyz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 52); 

                auto tg_yy_xyzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 53); 

                auto tg_yy_xzzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 54); 

                auto tg_yy_yyyy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 55); 

                auto tg_yy_yyyz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 56); 

                auto tg_yy_yyzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 57); 

                auto tg_yy_yzzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 58); 

                auto tg_yy_zzzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 59); 

                auto tg_yz_xxxx_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 60); 

                auto tg_yz_xxxy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 61); 

                auto tg_yz_xxxz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 62); 

                auto tg_yz_xxyy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 63); 

                auto tg_yz_xxyz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 64); 

                auto tg_yz_xxzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 65); 

                auto tg_yz_xyyy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 66); 

                auto tg_yz_xyyz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 67); 

                auto tg_yz_xyzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 68); 

                auto tg_yz_xzzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 69); 

                auto tg_yz_yyyy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 70); 

                auto tg_yz_yyyz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 71); 

                auto tg_yz_yyzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 72); 

                auto tg_yz_yzzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 73); 

                auto tg_yz_zzzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 74); 

                auto tg_zz_xxxx_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 75); 

                auto tg_zz_xxxy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 76); 

                auto tg_zz_xxxz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 77); 

                auto tg_zz_xxyy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 78); 

                auto tg_zz_xxyz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 79); 

                auto tg_zz_xxzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 80); 

                auto tg_zz_xyyy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 81); 

                auto tg_zz_xyyz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 82); 

                auto tg_zz_xyzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 83); 

                auto tg_zz_xzzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 84); 

                auto tg_zz_yyyy_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 85); 

                auto tg_zz_yyyz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 86); 

                auto tg_zz_yyzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 87); 

                auto tg_zz_yzzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 88); 

                auto tg_zz_zzzz_0 = primBuffer[pidx_g_2_4_m0].data(90 * idx + 89); 

                auto tg_yy_xxxx_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 45); 

                auto tg_yy_xxxy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 46); 

                auto tg_yy_xxxz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 47); 

                auto tg_yy_xxyy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 48); 

                auto tg_yy_xxyz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 49); 

                auto tg_yy_xxzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 50); 

                auto tg_yy_xyyy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 51); 

                auto tg_yy_xyyz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 52); 

                auto tg_yy_xyzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 53); 

                auto tg_yy_xzzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 54); 

                auto tg_yy_yyyy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 55); 

                auto tg_yy_yyyz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 56); 

                auto tg_yy_yyzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 57); 

                auto tg_yy_yzzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 58); 

                auto tg_yy_zzzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 59); 

                auto tg_yz_xxxx_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 60); 

                auto tg_yz_xxxy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 61); 

                auto tg_yz_xxxz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 62); 

                auto tg_yz_xxyy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 63); 

                auto tg_yz_xxyz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 64); 

                auto tg_yz_xxzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 65); 

                auto tg_yz_xyyy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 66); 

                auto tg_yz_xyyz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 67); 

                auto tg_yz_xyzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 68); 

                auto tg_yz_xzzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 69); 

                auto tg_yz_yyyy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 70); 

                auto tg_yz_yyyz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 71); 

                auto tg_yz_yyzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 72); 

                auto tg_yz_yzzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 73); 

                auto tg_yz_zzzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 74); 

                auto tg_zz_xxxx_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 75); 

                auto tg_zz_xxxy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 76); 

                auto tg_zz_xxxz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 77); 

                auto tg_zz_xxyy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 78); 

                auto tg_zz_xxyz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 79); 

                auto tg_zz_xxzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 80); 

                auto tg_zz_xyyy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 81); 

                auto tg_zz_xyyz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 82); 

                auto tg_zz_xyzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 83); 

                auto tg_zz_xzzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 84); 

                auto tg_zz_yyyy_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 85); 

                auto tg_zz_yyyz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 86); 

                auto tg_zz_yyzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 87); 

                auto tg_zz_yzzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 88); 

                auto tg_zz_zzzz_1 = primBuffer[pidx_g_2_4_m1].data(90 * idx + 89); 

                auto tg_y_xxxx_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 15); 

                auto tg_y_xxxy_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 16); 

                auto tg_y_xxxz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 17); 

                auto tg_y_xxyy_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 18); 

                auto tg_y_xxyz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 19); 

                auto tg_y_xxzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 20); 

                auto tg_y_xyyy_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 21); 

                auto tg_y_xyyz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 22); 

                auto tg_y_xyzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 23); 

                auto tg_y_xzzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 24); 

                auto tg_y_yyyy_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 25); 

                auto tg_y_yyyz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 26); 

                auto tg_y_yyzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 27); 

                auto tg_y_yzzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 28); 

                auto tg_y_zzzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 29); 

                auto tg_z_xxxx_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 30); 

                auto tg_z_xxxy_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 31); 

                auto tg_z_xxxz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 32); 

                auto tg_z_xxyy_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 33); 

                auto tg_z_xxyz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 34); 

                auto tg_z_xxzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 35); 

                auto tg_z_xyyy_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 36); 

                auto tg_z_xyyz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 37); 

                auto tg_z_xyzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 38); 

                auto tg_z_xzzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 39); 

                auto tg_z_yyyy_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 40); 

                auto tg_z_yyyz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 41); 

                auto tg_z_yyzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 42); 

                auto tg_z_yzzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 43); 

                auto tg_z_zzzz_0 = primBuffer[pidx_g_1_4_m0].data(45 * idx + 44); 

                auto tg_y_xxxx_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 15); 

                auto tg_y_xxxy_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 16); 

                auto tg_y_xxxz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 17); 

                auto tg_y_xxyy_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 18); 

                auto tg_y_xxyz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 19); 

                auto tg_y_xxzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 20); 

                auto tg_y_xyyy_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 21); 

                auto tg_y_xyyz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 22); 

                auto tg_y_xyzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 23); 

                auto tg_y_xzzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 24); 

                auto tg_y_yyyy_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 25); 

                auto tg_y_yyyz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 26); 

                auto tg_y_yyzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 27); 

                auto tg_y_yzzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 28); 

                auto tg_y_zzzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 29); 

                auto tg_z_xxxx_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 30); 

                auto tg_z_xxxy_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 31); 

                auto tg_z_xxxz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 32); 

                auto tg_z_xxyy_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 33); 

                auto tg_z_xxyz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 34); 

                auto tg_z_xxzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 35); 

                auto tg_z_xyyy_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 36); 

                auto tg_z_xyyz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 37); 

                auto tg_z_xyzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 38); 

                auto tg_z_xzzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 39); 

                auto tg_z_yyyy_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 40); 

                auto tg_z_yyyz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 41); 

                auto tg_z_yyzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 42); 

                auto tg_z_yzzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 43); 

                auto tg_z_zzzz_1 = primBuffer[pidx_g_1_4_m1].data(45 * idx + 44); 

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

                // set up pointers to integrals

                auto tg_xzz_xxxx_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 75); 

                auto tg_xzz_xxxy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 76); 

                auto tg_xzz_xxxz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 77); 

                auto tg_xzz_xxyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 78); 

                auto tg_xzz_xxyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 79); 

                auto tg_xzz_xxzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 80); 

                auto tg_xzz_xyyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 81); 

                auto tg_xzz_xyyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 82); 

                auto tg_xzz_xyzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 83); 

                auto tg_xzz_xzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 84); 

                auto tg_xzz_yyyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 85); 

                auto tg_xzz_yyyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 86); 

                auto tg_xzz_yyzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 87); 

                auto tg_xzz_yzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 88); 

                auto tg_xzz_zzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 89); 

                auto tg_yyy_xxxx_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 90); 

                auto tg_yyy_xxxy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 91); 

                auto tg_yyy_xxxz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 92); 

                auto tg_yyy_xxyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 93); 

                auto tg_yyy_xxyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 94); 

                auto tg_yyy_xxzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 95); 

                auto tg_yyy_xyyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 96); 

                auto tg_yyy_xyyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 97); 

                auto tg_yyy_xyzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 98); 

                auto tg_yyy_xzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 99); 

                auto tg_yyy_yyyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 100); 

                auto tg_yyy_yyyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 101); 

                auto tg_yyy_yyzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 102); 

                auto tg_yyy_yzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 103); 

                auto tg_yyy_zzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 104); 

                auto tg_yyz_xxxx_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 105); 

                auto tg_yyz_xxxy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 106); 

                auto tg_yyz_xxxz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 107); 

                auto tg_yyz_xxyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 108); 

                auto tg_yyz_xxyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 109); 

                auto tg_yyz_xxzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 110); 

                auto tg_yyz_xyyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 111); 

                auto tg_yyz_xyyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 112); 

                auto tg_yyz_xyzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 113); 

                auto tg_yyz_xzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 114); 

                auto tg_yyz_yyyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 115); 

                auto tg_yyz_yyyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 116); 

                auto tg_yyz_yyzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 117); 

                auto tg_yyz_yzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 118); 

                auto tg_yyz_zzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 119); 

                auto tg_yzz_xxxx_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 120); 

                auto tg_yzz_xxxy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 121); 

                auto tg_yzz_xxxz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 122); 

                auto tg_yzz_xxyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 123); 

                auto tg_yzz_xxyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 124); 

                auto tg_yzz_xxzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 125); 

                auto tg_yzz_xyyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 126); 

                auto tg_yzz_xyyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 127); 

                auto tg_yzz_xyzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 128); 

                auto tg_yzz_xzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 129); 

                auto tg_yzz_yyyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 130); 

                auto tg_yzz_yyyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 131); 

                auto tg_yzz_yyzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 132); 

                auto tg_yzz_yzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 133); 

                auto tg_yzz_zzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 134); 

                auto tg_zzz_xxxx_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 135); 

                auto tg_zzz_xxxy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 136); 

                auto tg_zzz_xxxz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 137); 

                auto tg_zzz_xxyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 138); 

                auto tg_zzz_xxyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 139); 

                auto tg_zzz_xxzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 140); 

                auto tg_zzz_xyyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 141); 

                auto tg_zzz_xyyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 142); 

                auto tg_zzz_xyzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 143); 

                auto tg_zzz_xzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 144); 

                auto tg_zzz_yyyy_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 145); 

                auto tg_zzz_yyyz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 146); 

                auto tg_zzz_yyzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 147); 

                auto tg_zzz_yzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 148); 

                auto tg_zzz_zzzz_0 = primBuffer[pidx_g_3_4_m0].data(150 * idx + 149); 

                // Batch of Integrals (75,150)

                #pragma omp simd aligned(fxn, fza, tg_xzz_xxxx_0, tg_xzz_xxxy_0, tg_xzz_xxxz_0, tg_xzz_xxyy_0, \
                                         tg_xzz_xxyz_0, tg_xzz_xxzz_0, tg_xzz_xyyy_0, tg_xzz_xyyz_0, tg_xzz_xyzz_0, \
                                         tg_xzz_xzzz_0, tg_xzz_yyyy_0, tg_xzz_yyyz_0, tg_xzz_yyzz_0, tg_xzz_yzzz_0, \
                                         tg_xzz_zzzz_0, tg_y_xxxx_0, tg_y_xxxx_1, tg_y_xxxy_0, tg_y_xxxy_1, tg_y_xxxz_0, \
                                         tg_y_xxxz_1, tg_y_xxyy_0, tg_y_xxyy_1, tg_y_xxyz_0, tg_y_xxyz_1, tg_y_xxzz_0, \
                                         tg_y_xxzz_1, tg_y_xyyy_0, tg_y_xyyy_1, tg_y_xyyz_0, tg_y_xyyz_1, tg_y_xyzz_0, \
                                         tg_y_xyzz_1, tg_y_xzzz_0, tg_y_xzzz_1, tg_y_yyyy_0, tg_y_yyyy_1, tg_y_yyyz_0, \
                                         tg_y_yyyz_1, tg_y_yyzz_0, tg_y_yyzz_1, tg_y_yzzz_0, tg_y_yzzz_1, tg_y_zzzz_0, \
                                         tg_y_zzzz_1, tg_yy_xxx_1, tg_yy_xxxx_0, tg_yy_xxxx_1, tg_yy_xxxy_0, tg_yy_xxxy_1, \
                                         tg_yy_xxxz_0, tg_yy_xxxz_1, tg_yy_xxy_1, tg_yy_xxyy_0, tg_yy_xxyy_1, tg_yy_xxyz_0, \
                                         tg_yy_xxyz_1, tg_yy_xxz_1, tg_yy_xxzz_0, tg_yy_xxzz_1, tg_yy_xyy_1, tg_yy_xyyy_0, \
                                         tg_yy_xyyy_1, tg_yy_xyyz_0, tg_yy_xyyz_1, tg_yy_xyz_1, tg_yy_xyzz_0, tg_yy_xyzz_1, \
                                         tg_yy_xzz_1, tg_yy_xzzz_0, tg_yy_xzzz_1, tg_yy_yyy_1, tg_yy_yyyy_0, tg_yy_yyyy_1, \
                                         tg_yy_yyyz_0, tg_yy_yyyz_1, tg_yy_yyz_1, tg_yy_yyzz_0, tg_yy_yyzz_1, tg_yy_yzz_1, \
                                         tg_yy_yzzz_0, tg_yy_yzzz_1, tg_yy_zzz_1, tg_yy_zzzz_0, tg_yy_zzzz_1, tg_yyy_xxxx_0, \
                                         tg_yyy_xxxy_0, tg_yyy_xxxz_0, tg_yyy_xxyy_0, tg_yyy_xxyz_0, tg_yyy_xxzz_0, \
                                         tg_yyy_xyyy_0, tg_yyy_xyyz_0, tg_yyy_xyzz_0, tg_yyy_xzzz_0, tg_yyy_yyyy_0, \
                                         tg_yyy_yyyz_0, tg_yyy_yyzz_0, tg_yyy_yzzz_0, tg_yyy_zzzz_0, tg_yyz_xxxx_0, \
                                         tg_yyz_xxxy_0, tg_yyz_xxxz_0, tg_yyz_xxyy_0, tg_yyz_xxyz_0, tg_yyz_xxzz_0, \
                                         tg_yyz_xyyy_0, tg_yyz_xyyz_0, tg_yyz_xyzz_0, tg_yyz_xzzz_0, tg_yyz_yyyy_0, \
                                         tg_yyz_yyyz_0, tg_yyz_yyzz_0, tg_yyz_yzzz_0, tg_yyz_zzzz_0, tg_yz_xxx_1, \
                                         tg_yz_xxxx_0, tg_yz_xxxx_1, tg_yz_xxxy_0, tg_yz_xxxy_1, tg_yz_xxxz_0, tg_yz_xxxz_1, \
                                         tg_yz_xxy_1, tg_yz_xxyy_0, tg_yz_xxyy_1, tg_yz_xxyz_0, tg_yz_xxyz_1, tg_yz_xxz_1, \
                                         tg_yz_xxzz_0, tg_yz_xxzz_1, tg_yz_xyy_1, tg_yz_xyyy_0, tg_yz_xyyy_1, tg_yz_xyyz_0, \
                                         tg_yz_xyyz_1, tg_yz_xyz_1, tg_yz_xyzz_0, tg_yz_xyzz_1, tg_yz_xzz_1, tg_yz_xzzz_0, \
                                         tg_yz_xzzz_1, tg_yz_yyy_1, tg_yz_yyyy_0, tg_yz_yyyy_1, tg_yz_yyyz_0, tg_yz_yyyz_1, \
                                         tg_yz_yyz_1, tg_yz_yyzz_0, tg_yz_yyzz_1, tg_yz_yzz_1, tg_yz_yzzz_0, tg_yz_yzzz_1, \
                                         tg_yz_zzz_1, tg_yz_zzzz_0, tg_yz_zzzz_1, tg_yzz_xxxx_0, tg_yzz_xxxy_0, \
                                         tg_yzz_xxxz_0, tg_yzz_xxyy_0, tg_yzz_xxyz_0, tg_yzz_xxzz_0, tg_yzz_xyyy_0, \
                                         tg_yzz_xyyz_0, tg_yzz_xyzz_0, tg_yzz_xzzz_0, tg_yzz_yyyy_0, tg_yzz_yyyz_0, \
                                         tg_yzz_yyzz_0, tg_yzz_yzzz_0, tg_yzz_zzzz_0, tg_z_xxxx_0, tg_z_xxxx_1, tg_z_xxxy_0, \
                                         tg_z_xxxy_1, tg_z_xxxz_0, tg_z_xxxz_1, tg_z_xxyy_0, tg_z_xxyy_1, tg_z_xxyz_0, \
                                         tg_z_xxyz_1, tg_z_xxzz_0, tg_z_xxzz_1, tg_z_xyyy_0, tg_z_xyyy_1, tg_z_xyyz_0, \
                                         tg_z_xyyz_1, tg_z_xyzz_0, tg_z_xyzz_1, tg_z_xzzz_0, tg_z_xzzz_1, tg_z_yyyy_0, \
                                         tg_z_yyyy_1, tg_z_yyyz_0, tg_z_yyyz_1, tg_z_yyzz_0, tg_z_yyzz_1, tg_z_yzzz_0, \
                                         tg_z_yzzz_1, tg_z_zzzz_0, tg_z_zzzz_1, tg_zz_xxx_1, tg_zz_xxxx_0, tg_zz_xxxx_1, \
                                         tg_zz_xxxy_0, tg_zz_xxxy_1, tg_zz_xxxz_0, tg_zz_xxxz_1, tg_zz_xxy_1, tg_zz_xxyy_0, \
                                         tg_zz_xxyy_1, tg_zz_xxyz_0, tg_zz_xxyz_1, tg_zz_xxz_1, tg_zz_xxzz_0, tg_zz_xxzz_1, \
                                         tg_zz_xyy_1, tg_zz_xyyy_0, tg_zz_xyyy_1, tg_zz_xyyz_0, tg_zz_xyyz_1, tg_zz_xyz_1, \
                                         tg_zz_xyzz_0, tg_zz_xyzz_1, tg_zz_xzz_1, tg_zz_xzzz_0, tg_zz_xzzz_1, tg_zz_yyy_1, \
                                         tg_zz_yyyy_0, tg_zz_yyyy_1, tg_zz_yyyz_0, tg_zz_yyyz_1, tg_zz_yyz_1, tg_zz_yyzz_0, \
                                         tg_zz_yyzz_1, tg_zz_yzz_1, tg_zz_yzzz_0, tg_zz_yzzz_1, tg_zz_zzz_1, tg_zz_zzzz_0, \
                                         tg_zz_zzzz_1, tg_zzz_xxxx_0, tg_zzz_xxxy_0, tg_zzz_xxxz_0, tg_zzz_xxyy_0, \
                                         tg_zzz_xxyz_0, tg_zzz_xxzz_0, tg_zzz_xyyy_0, tg_zzz_xyyz_0, tg_zzz_xyzz_0, \
                                         tg_zzz_xzzz_0, tg_zzz_yyyy_0, tg_zzz_yyyz_0, tg_zzz_yyzz_0, tg_zzz_yzzz_0, \
                                         tg_zzz_zzzz_0, wp_x, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xzz_xxxx_0[j] = pb_x * tg_zz_xxxx_0[j] + fr * tg_zz_xxxx_1[j] + 2.0 * fl1_fxn * tg_zz_xxx_1[j];

                    tg_xzz_xxxy_0[j] = pb_x * tg_zz_xxxy_0[j] + fr * tg_zz_xxxy_1[j] + 1.5 * fl1_fxn * tg_zz_xxy_1[j];

                    tg_xzz_xxxz_0[j] = pb_x * tg_zz_xxxz_0[j] + fr * tg_zz_xxxz_1[j] + 1.5 * fl1_fxn * tg_zz_xxz_1[j];

                    tg_xzz_xxyy_0[j] = pb_x * tg_zz_xxyy_0[j] + fr * tg_zz_xxyy_1[j] + fl1_fxn * tg_zz_xyy_1[j];

                    tg_xzz_xxyz_0[j] = pb_x * tg_zz_xxyz_0[j] + fr * tg_zz_xxyz_1[j] + fl1_fxn * tg_zz_xyz_1[j];

                    tg_xzz_xxzz_0[j] = pb_x * tg_zz_xxzz_0[j] + fr * tg_zz_xxzz_1[j] + fl1_fxn * tg_zz_xzz_1[j];

                    tg_xzz_xyyy_0[j] = pb_x * tg_zz_xyyy_0[j] + fr * tg_zz_xyyy_1[j] + 0.5 * fl1_fxn * tg_zz_yyy_1[j];

                    tg_xzz_xyyz_0[j] = pb_x * tg_zz_xyyz_0[j] + fr * tg_zz_xyyz_1[j] + 0.5 * fl1_fxn * tg_zz_yyz_1[j];

                    tg_xzz_xyzz_0[j] = pb_x * tg_zz_xyzz_0[j] + fr * tg_zz_xyzz_1[j] + 0.5 * fl1_fxn * tg_zz_yzz_1[j];

                    tg_xzz_xzzz_0[j] = pb_x * tg_zz_xzzz_0[j] + fr * tg_zz_xzzz_1[j] + 0.5 * fl1_fxn * tg_zz_zzz_1[j];

                    tg_xzz_yyyy_0[j] = pb_x * tg_zz_yyyy_0[j] + fr * tg_zz_yyyy_1[j];

                    tg_xzz_yyyz_0[j] = pb_x * tg_zz_yyyz_0[j] + fr * tg_zz_yyyz_1[j];

                    tg_xzz_yyzz_0[j] = pb_x * tg_zz_yyzz_0[j] + fr * tg_zz_yyzz_1[j];

                    tg_xzz_yzzz_0[j] = pb_x * tg_zz_yzzz_0[j] + fr * tg_zz_yzzz_1[j];

                    tg_xzz_zzzz_0[j] = pb_x * tg_zz_zzzz_0[j] + fr * tg_zz_zzzz_1[j];

                    fr = wp_y[j]; 

                    tg_yyy_xxxx_0[j] = pb_y * tg_yy_xxxx_0[j] + fr * tg_yy_xxxx_1[j] + fl1_fx * (tg_y_xxxx_0[j] - tg_y_xxxx_1[j] * fl1_fza);

                    tg_yyy_xxxy_0[j] = pb_y * tg_yy_xxxy_0[j] + fr * tg_yy_xxxy_1[j] + fl1_fx * (tg_y_xxxy_0[j] - tg_y_xxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yy_xxx_1[j];

                    tg_yyy_xxxz_0[j] = pb_y * tg_yy_xxxz_0[j] + fr * tg_yy_xxxz_1[j] + fl1_fx * (tg_y_xxxz_0[j] - tg_y_xxxz_1[j] * fl1_fza);

                    tg_yyy_xxyy_0[j] = pb_y * tg_yy_xxyy_0[j] + fr * tg_yy_xxyy_1[j] + fl1_fx * (tg_y_xxyy_0[j] - tg_y_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_yy_xxy_1[j];

                    tg_yyy_xxyz_0[j] = pb_y * tg_yy_xxyz_0[j] + fr * tg_yy_xxyz_1[j] + fl1_fx * (tg_y_xxyz_0[j] - tg_y_xxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yy_xxz_1[j];

                    tg_yyy_xxzz_0[j] = pb_y * tg_yy_xxzz_0[j] + fr * tg_yy_xxzz_1[j] + fl1_fx * (tg_y_xxzz_0[j] - tg_y_xxzz_1[j] * fl1_fza);

                    tg_yyy_xyyy_0[j] = pb_y * tg_yy_xyyy_0[j] + fr * tg_yy_xyyy_1[j] + fl1_fx * (tg_y_xyyy_0[j] - tg_y_xyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yy_xyy_1[j];

                    tg_yyy_xyyz_0[j] = pb_y * tg_yy_xyyz_0[j] + fr * tg_yy_xyyz_1[j] + fl1_fx * (tg_y_xyyz_0[j] - tg_y_xyyz_1[j] * fl1_fza) + fl1_fxn * tg_yy_xyz_1[j];

                    tg_yyy_xyzz_0[j] = pb_y * tg_yy_xyzz_0[j] + fr * tg_yy_xyzz_1[j] + fl1_fx * (tg_y_xyzz_0[j] - tg_y_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yy_xzz_1[j];

                    tg_yyy_xzzz_0[j] = pb_y * tg_yy_xzzz_0[j] + fr * tg_yy_xzzz_1[j] + fl1_fx * (tg_y_xzzz_0[j] - tg_y_xzzz_1[j] * fl1_fza);

                    tg_yyy_yyyy_0[j] = pb_y * tg_yy_yyyy_0[j] + fr * tg_yy_yyyy_1[j] + fl1_fx * (tg_y_yyyy_0[j] - tg_y_yyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yy_yyy_1[j];

                    tg_yyy_yyyz_0[j] = pb_y * tg_yy_yyyz_0[j] + fr * tg_yy_yyyz_1[j] + fl1_fx * (tg_y_yyyz_0[j] - tg_y_yyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yy_yyz_1[j];

                    tg_yyy_yyzz_0[j] = pb_y * tg_yy_yyzz_0[j] + fr * tg_yy_yyzz_1[j] + fl1_fx * (tg_y_yyzz_0[j] - tg_y_yyzz_1[j] * fl1_fza) + fl1_fxn * tg_yy_yzz_1[j];

                    tg_yyy_yzzz_0[j] = pb_y * tg_yy_yzzz_0[j] + fr * tg_yy_yzzz_1[j] + fl1_fx * (tg_y_yzzz_0[j] - tg_y_yzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yy_zzz_1[j];

                    tg_yyy_zzzz_0[j] = pb_y * tg_yy_zzzz_0[j] + fr * tg_yy_zzzz_1[j] + fl1_fx * (tg_y_zzzz_0[j] - tg_y_zzzz_1[j] * fl1_fza);

                    tg_yyz_xxxx_0[j] = pb_y * tg_yz_xxxx_0[j] + fr * tg_yz_xxxx_1[j] + 0.5 * fl1_fx * (tg_z_xxxx_0[j] - tg_z_xxxx_1[j] * fl1_fza);

                    tg_yyz_xxxy_0[j] = pb_y * tg_yz_xxxy_0[j] + fr * tg_yz_xxxy_1[j] + 0.5 * fl1_fx * (tg_z_xxxy_0[j] - tg_z_xxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yz_xxx_1[j];

                    tg_yyz_xxxz_0[j] = pb_y * tg_yz_xxxz_0[j] + fr * tg_yz_xxxz_1[j] + 0.5 * fl1_fx * (tg_z_xxxz_0[j] - tg_z_xxxz_1[j] * fl1_fza);

                    tg_yyz_xxyy_0[j] = pb_y * tg_yz_xxyy_0[j] + fr * tg_yz_xxyy_1[j] + 0.5 * fl1_fx * (tg_z_xxyy_0[j] - tg_z_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_yz_xxy_1[j];

                    tg_yyz_xxyz_0[j] = pb_y * tg_yz_xxyz_0[j] + fr * tg_yz_xxyz_1[j] + 0.5 * fl1_fx * (tg_z_xxyz_0[j] - tg_z_xxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yz_xxz_1[j];

                    tg_yyz_xxzz_0[j] = pb_y * tg_yz_xxzz_0[j] + fr * tg_yz_xxzz_1[j] + 0.5 * fl1_fx * (tg_z_xxzz_0[j] - tg_z_xxzz_1[j] * fl1_fza);

                    tg_yyz_xyyy_0[j] = pb_y * tg_yz_xyyy_0[j] + fr * tg_yz_xyyy_1[j] + 0.5 * fl1_fx * (tg_z_xyyy_0[j] - tg_z_xyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yz_xyy_1[j];

                    tg_yyz_xyyz_0[j] = pb_y * tg_yz_xyyz_0[j] + fr * tg_yz_xyyz_1[j] + 0.5 * fl1_fx * (tg_z_xyyz_0[j] - tg_z_xyyz_1[j] * fl1_fza) + fl1_fxn * tg_yz_xyz_1[j];

                    tg_yyz_xyzz_0[j] = pb_y * tg_yz_xyzz_0[j] + fr * tg_yz_xyzz_1[j] + 0.5 * fl1_fx * (tg_z_xyzz_0[j] - tg_z_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yz_xzz_1[j];

                    tg_yyz_xzzz_0[j] = pb_y * tg_yz_xzzz_0[j] + fr * tg_yz_xzzz_1[j] + 0.5 * fl1_fx * (tg_z_xzzz_0[j] - tg_z_xzzz_1[j] * fl1_fza);

                    tg_yyz_yyyy_0[j] = pb_y * tg_yz_yyyy_0[j] + fr * tg_yz_yyyy_1[j] + 0.5 * fl1_fx * (tg_z_yyyy_0[j] - tg_z_yyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yz_yyy_1[j];

                    tg_yyz_yyyz_0[j] = pb_y * tg_yz_yyyz_0[j] + fr * tg_yz_yyyz_1[j] + 0.5 * fl1_fx * (tg_z_yyyz_0[j] - tg_z_yyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yz_yyz_1[j];

                    tg_yyz_yyzz_0[j] = pb_y * tg_yz_yyzz_0[j] + fr * tg_yz_yyzz_1[j] + 0.5 * fl1_fx * (tg_z_yyzz_0[j] - tg_z_yyzz_1[j] * fl1_fza) + fl1_fxn * tg_yz_yzz_1[j];

                    tg_yyz_yzzz_0[j] = pb_y * tg_yz_yzzz_0[j] + fr * tg_yz_yzzz_1[j] + 0.5 * fl1_fx * (tg_z_yzzz_0[j] - tg_z_yzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yz_zzz_1[j];

                    tg_yyz_zzzz_0[j] = pb_y * tg_yz_zzzz_0[j] + fr * tg_yz_zzzz_1[j] + 0.5 * fl1_fx * (tg_z_zzzz_0[j] - tg_z_zzzz_1[j] * fl1_fza);

                    tg_yzz_xxxx_0[j] = pb_y * tg_zz_xxxx_0[j] + fr * tg_zz_xxxx_1[j];

                    tg_yzz_xxxy_0[j] = pb_y * tg_zz_xxxy_0[j] + fr * tg_zz_xxxy_1[j] + 0.5 * fl1_fxn * tg_zz_xxx_1[j];

                    tg_yzz_xxxz_0[j] = pb_y * tg_zz_xxxz_0[j] + fr * tg_zz_xxxz_1[j];

                    tg_yzz_xxyy_0[j] = pb_y * tg_zz_xxyy_0[j] + fr * tg_zz_xxyy_1[j] + fl1_fxn * tg_zz_xxy_1[j];

                    tg_yzz_xxyz_0[j] = pb_y * tg_zz_xxyz_0[j] + fr * tg_zz_xxyz_1[j] + 0.5 * fl1_fxn * tg_zz_xxz_1[j];

                    tg_yzz_xxzz_0[j] = pb_y * tg_zz_xxzz_0[j] + fr * tg_zz_xxzz_1[j];

                    tg_yzz_xyyy_0[j] = pb_y * tg_zz_xyyy_0[j] + fr * tg_zz_xyyy_1[j] + 1.5 * fl1_fxn * tg_zz_xyy_1[j];

                    tg_yzz_xyyz_0[j] = pb_y * tg_zz_xyyz_0[j] + fr * tg_zz_xyyz_1[j] + fl1_fxn * tg_zz_xyz_1[j];

                    tg_yzz_xyzz_0[j] = pb_y * tg_zz_xyzz_0[j] + fr * tg_zz_xyzz_1[j] + 0.5 * fl1_fxn * tg_zz_xzz_1[j];

                    tg_yzz_xzzz_0[j] = pb_y * tg_zz_xzzz_0[j] + fr * tg_zz_xzzz_1[j];

                    tg_yzz_yyyy_0[j] = pb_y * tg_zz_yyyy_0[j] + fr * tg_zz_yyyy_1[j] + 2.0 * fl1_fxn * tg_zz_yyy_1[j];

                    tg_yzz_yyyz_0[j] = pb_y * tg_zz_yyyz_0[j] + fr * tg_zz_yyyz_1[j] + 1.5 * fl1_fxn * tg_zz_yyz_1[j];

                    tg_yzz_yyzz_0[j] = pb_y * tg_zz_yyzz_0[j] + fr * tg_zz_yyzz_1[j] + fl1_fxn * tg_zz_yzz_1[j];

                    tg_yzz_yzzz_0[j] = pb_y * tg_zz_yzzz_0[j] + fr * tg_zz_yzzz_1[j] + 0.5 * fl1_fxn * tg_zz_zzz_1[j];

                    tg_yzz_zzzz_0[j] = pb_y * tg_zz_zzzz_0[j] + fr * tg_zz_zzzz_1[j];

                    fr = wp_z[j]; 

                    tg_zzz_xxxx_0[j] = pb_z * tg_zz_xxxx_0[j] + fr * tg_zz_xxxx_1[j] + fl1_fx * (tg_z_xxxx_0[j] - tg_z_xxxx_1[j] * fl1_fza);

                    tg_zzz_xxxy_0[j] = pb_z * tg_zz_xxxy_0[j] + fr * tg_zz_xxxy_1[j] + fl1_fx * (tg_z_xxxy_0[j] - tg_z_xxxy_1[j] * fl1_fza);

                    tg_zzz_xxxz_0[j] = pb_z * tg_zz_xxxz_0[j] + fr * tg_zz_xxxz_1[j] + fl1_fx * (tg_z_xxxz_0[j] - tg_z_xxxz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zz_xxx_1[j];

                    tg_zzz_xxyy_0[j] = pb_z * tg_zz_xxyy_0[j] + fr * tg_zz_xxyy_1[j] + fl1_fx * (tg_z_xxyy_0[j] - tg_z_xxyy_1[j] * fl1_fza);

                    tg_zzz_xxyz_0[j] = pb_z * tg_zz_xxyz_0[j] + fr * tg_zz_xxyz_1[j] + fl1_fx * (tg_z_xxyz_0[j] - tg_z_xxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zz_xxy_1[j];

                    tg_zzz_xxzz_0[j] = pb_z * tg_zz_xxzz_0[j] + fr * tg_zz_xxzz_1[j] + fl1_fx * (tg_z_xxzz_0[j] - tg_z_xxzz_1[j] * fl1_fza) + fl1_fxn * tg_zz_xxz_1[j];

                    tg_zzz_xyyy_0[j] = pb_z * tg_zz_xyyy_0[j] + fr * tg_zz_xyyy_1[j] + fl1_fx * (tg_z_xyyy_0[j] - tg_z_xyyy_1[j] * fl1_fza);

                    tg_zzz_xyyz_0[j] = pb_z * tg_zz_xyyz_0[j] + fr * tg_zz_xyyz_1[j] + fl1_fx * (tg_z_xyyz_0[j] - tg_z_xyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zz_xyy_1[j];

                    tg_zzz_xyzz_0[j] = pb_z * tg_zz_xyzz_0[j] + fr * tg_zz_xyzz_1[j] + fl1_fx * (tg_z_xyzz_0[j] - tg_z_xyzz_1[j] * fl1_fza) + fl1_fxn * tg_zz_xyz_1[j];

                    tg_zzz_xzzz_0[j] = pb_z * tg_zz_xzzz_0[j] + fr * tg_zz_xzzz_1[j] + fl1_fx * (tg_z_xzzz_0[j] - tg_z_xzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zz_xzz_1[j];

                    tg_zzz_yyyy_0[j] = pb_z * tg_zz_yyyy_0[j] + fr * tg_zz_yyyy_1[j] + fl1_fx * (tg_z_yyyy_0[j] - tg_z_yyyy_1[j] * fl1_fza);

                    tg_zzz_yyyz_0[j] = pb_z * tg_zz_yyyz_0[j] + fr * tg_zz_yyyz_1[j] + fl1_fx * (tg_z_yyyz_0[j] - tg_z_yyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zz_yyy_1[j];

                    tg_zzz_yyzz_0[j] = pb_z * tg_zz_yyzz_0[j] + fr * tg_zz_yyzz_1[j] + fl1_fx * (tg_z_yyzz_0[j] - tg_z_yyzz_1[j] * fl1_fza) + fl1_fxn * tg_zz_yyz_1[j];

                    tg_zzz_yzzz_0[j] = pb_z * tg_zz_yzzz_0[j] + fr * tg_zz_yzzz_1[j] + fl1_fx * (tg_z_yzzz_0[j] - tg_z_yzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zz_yzz_1[j];

                    tg_zzz_zzzz_0[j] = pb_z * tg_zz_zzzz_0[j] + fr * tg_zz_zzzz_1[j] + fl1_fx * (tg_z_zzzz_0[j] - tg_z_zzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zz_zzz_1[j];
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

