//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectronRepulsionRecFuncForLG.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSLSG(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSLSG_0_49(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSLSG_49_98(primBuffer,
                                                       recursionMap,
                                                       osFactors,
                                                       wpDistances, 
                                                       braGtoPairsBlock,
                                                       ketGtoPairsBlock,
                                                       nKetPrimPairs,
                                                       iContrPair); 

        erirecfunc::compElectronRepulsionForSLSG_98_147(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 

        erirecfunc::compElectronRepulsionForSLSG_147_195(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSLSG_195_243(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSLSG_243_291(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSLSG_291_339(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSLSG_339_387(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSLSG_387_435(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSLSG_435_483(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSLSG_483_531(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSLSG_531_579(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSLSG_579_627(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSLSG_627_675(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 
    }

    void
    compElectronRepulsionForSLSG_0_49(      CMemBlock2D<double>& primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& wpDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetPrimPairs,
                                      const int32_t              iContrPair)
    {
        // Batch of Integrals (0,49)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {8, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {8, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_7_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_6_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_6_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_7_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_xxxxxxx_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx); 

                auto tg_xxxxxxx_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 1); 

                auto tg_xxxxxxx_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 2); 

                auto tg_xxxxxxx_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 3); 

                auto tg_xxxxxxx_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 4); 

                auto tg_xxxxxxx_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 5); 

                auto tg_xxxxxxx_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 6); 

                auto tg_xxxxxxx_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 7); 

                auto tg_xxxxxxx_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 8); 

                auto tg_xxxxxxx_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 9); 

                auto tg_xxxxxxx_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 10); 

                auto tg_xxxxxxx_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 11); 

                auto tg_xxxxxxx_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 12); 

                auto tg_xxxxxxx_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 13); 

                auto tg_xxxxxxx_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 14); 

                auto tg_xxxxxxy_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 15); 

                auto tg_xxxxxxy_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 16); 

                auto tg_xxxxxxy_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 17); 

                auto tg_xxxxxxy_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 18); 

                auto tg_xxxxxxy_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 19); 

                auto tg_xxxxxxy_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 20); 

                auto tg_xxxxxxy_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 21); 

                auto tg_xxxxxxy_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 22); 

                auto tg_xxxxxxy_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 23); 

                auto tg_xxxxxxy_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 24); 

                auto tg_xxxxxxy_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 25); 

                auto tg_xxxxxxy_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 26); 

                auto tg_xxxxxxy_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 27); 

                auto tg_xxxxxxy_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 28); 

                auto tg_xxxxxxy_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 29); 

                auto tg_xxxxxxz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 30); 

                auto tg_xxxxxxz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 31); 

                auto tg_xxxxxxz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 32); 

                auto tg_xxxxxxz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 33); 

                auto tg_xxxxxxz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 34); 

                auto tg_xxxxxxz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 35); 

                auto tg_xxxxxxz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 36); 

                auto tg_xxxxxxz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 37); 

                auto tg_xxxxxxz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 38); 

                auto tg_xxxxxxz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 39); 

                auto tg_xxxxxxz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 40); 

                auto tg_xxxxxxz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 41); 

                auto tg_xxxxxxz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 42); 

                auto tg_xxxxxxz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 43); 

                auto tg_xxxxxxz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 44); 

                auto tg_xxxxxyy_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 45); 

                auto tg_xxxxxyy_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 46); 

                auto tg_xxxxxyy_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 47); 

                auto tg_xxxxxyy_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 48); 

                auto tg_xxxxxxx_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx); 

                auto tg_xxxxxxx_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 1); 

                auto tg_xxxxxxx_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 2); 

                auto tg_xxxxxxx_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 3); 

                auto tg_xxxxxxx_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 4); 

                auto tg_xxxxxxx_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 5); 

                auto tg_xxxxxxx_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 6); 

                auto tg_xxxxxxx_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 7); 

                auto tg_xxxxxxx_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 8); 

                auto tg_xxxxxxx_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 9); 

                auto tg_xxxxxxx_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 10); 

                auto tg_xxxxxxx_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 11); 

                auto tg_xxxxxxx_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 12); 

                auto tg_xxxxxxx_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 13); 

                auto tg_xxxxxxx_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 14); 

                auto tg_xxxxxxy_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 15); 

                auto tg_xxxxxxy_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 16); 

                auto tg_xxxxxxy_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 17); 

                auto tg_xxxxxxy_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 18); 

                auto tg_xxxxxxy_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 19); 

                auto tg_xxxxxxy_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 20); 

                auto tg_xxxxxxy_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 21); 

                auto tg_xxxxxxy_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 22); 

                auto tg_xxxxxxy_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 23); 

                auto tg_xxxxxxy_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 24); 

                auto tg_xxxxxxy_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 25); 

                auto tg_xxxxxxy_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 26); 

                auto tg_xxxxxxy_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 27); 

                auto tg_xxxxxxy_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 28); 

                auto tg_xxxxxxy_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 29); 

                auto tg_xxxxxxz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 30); 

                auto tg_xxxxxxz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 31); 

                auto tg_xxxxxxz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 32); 

                auto tg_xxxxxxz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 33); 

                auto tg_xxxxxxz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 34); 

                auto tg_xxxxxxz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 35); 

                auto tg_xxxxxxz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 36); 

                auto tg_xxxxxxz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 37); 

                auto tg_xxxxxxz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 38); 

                auto tg_xxxxxxz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 39); 

                auto tg_xxxxxxz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 40); 

                auto tg_xxxxxxz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 41); 

                auto tg_xxxxxxz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 42); 

                auto tg_xxxxxxz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 43); 

                auto tg_xxxxxxz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 44); 

                auto tg_xxxxxyy_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 45); 

                auto tg_xxxxxyy_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 46); 

                auto tg_xxxxxyy_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 47); 

                auto tg_xxxxxyy_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 48); 

                auto tg_xxxxxx_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx); 

                auto tg_xxxxxx_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 1); 

                auto tg_xxxxxx_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 2); 

                auto tg_xxxxxx_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 3); 

                auto tg_xxxxxx_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 4); 

                auto tg_xxxxxx_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 5); 

                auto tg_xxxxxx_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 6); 

                auto tg_xxxxxx_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 7); 

                auto tg_xxxxxx_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 8); 

                auto tg_xxxxxx_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 9); 

                auto tg_xxxxxx_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 10); 

                auto tg_xxxxxx_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 11); 

                auto tg_xxxxxx_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 12); 

                auto tg_xxxxxx_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 13); 

                auto tg_xxxxxx_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 14); 

                auto tg_xxxxxy_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 15); 

                auto tg_xxxxxy_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 16); 

                auto tg_xxxxxy_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 17); 

                auto tg_xxxxxy_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 18); 

                auto tg_xxxxxy_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 19); 

                auto tg_xxxxxy_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 20); 

                auto tg_xxxxxy_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 21); 

                auto tg_xxxxxy_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 22); 

                auto tg_xxxxxy_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 23); 

                auto tg_xxxxxy_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 24); 

                auto tg_xxxxxy_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 25); 

                auto tg_xxxxxy_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 26); 

                auto tg_xxxxxy_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 27); 

                auto tg_xxxxxy_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 28); 

                auto tg_xxxxxy_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 29); 

                auto tg_xxxxxz_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 30); 

                auto tg_xxxxxz_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 31); 

                auto tg_xxxxxz_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 32); 

                auto tg_xxxxxz_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 33); 

                auto tg_xxxxxz_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 34); 

                auto tg_xxxxxz_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 35); 

                auto tg_xxxxxz_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 36); 

                auto tg_xxxxxz_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 37); 

                auto tg_xxxxxz_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 38); 

                auto tg_xxxxxz_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 39); 

                auto tg_xxxxxz_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 40); 

                auto tg_xxxxxz_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 41); 

                auto tg_xxxxxz_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 42); 

                auto tg_xxxxxz_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 43); 

                auto tg_xxxxxz_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 44); 

                auto tg_xxxxyy_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 45); 

                auto tg_xxxxyy_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 46); 

                auto tg_xxxxyy_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 47); 

                auto tg_xxxxyy_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 48); 

                auto tg_xxxxxx_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx); 

                auto tg_xxxxxx_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 1); 

                auto tg_xxxxxx_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 2); 

                auto tg_xxxxxx_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 3); 

                auto tg_xxxxxx_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 4); 

                auto tg_xxxxxx_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 5); 

                auto tg_xxxxxx_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 6); 

                auto tg_xxxxxx_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 7); 

                auto tg_xxxxxx_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 8); 

                auto tg_xxxxxx_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 9); 

                auto tg_xxxxxx_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 10); 

                auto tg_xxxxxx_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 11); 

                auto tg_xxxxxx_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 12); 

                auto tg_xxxxxx_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 13); 

                auto tg_xxxxxx_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 14); 

                auto tg_xxxxxy_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 15); 

                auto tg_xxxxxy_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 16); 

                auto tg_xxxxxy_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 17); 

                auto tg_xxxxxy_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 18); 

                auto tg_xxxxxy_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 19); 

                auto tg_xxxxxy_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 20); 

                auto tg_xxxxxy_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 21); 

                auto tg_xxxxxy_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 22); 

                auto tg_xxxxxy_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 23); 

                auto tg_xxxxxy_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 24); 

                auto tg_xxxxxy_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 25); 

                auto tg_xxxxxy_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 26); 

                auto tg_xxxxxy_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 27); 

                auto tg_xxxxxy_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 28); 

                auto tg_xxxxxy_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 29); 

                auto tg_xxxxxz_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 30); 

                auto tg_xxxxxz_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 31); 

                auto tg_xxxxxz_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 32); 

                auto tg_xxxxxz_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 33); 

                auto tg_xxxxxz_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 34); 

                auto tg_xxxxxz_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 35); 

                auto tg_xxxxxz_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 36); 

                auto tg_xxxxxz_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 37); 

                auto tg_xxxxxz_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 38); 

                auto tg_xxxxxz_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 39); 

                auto tg_xxxxxz_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 40); 

                auto tg_xxxxxz_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 41); 

                auto tg_xxxxxz_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 42); 

                auto tg_xxxxxz_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 43); 

                auto tg_xxxxxz_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 44); 

                auto tg_xxxxyy_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 45); 

                auto tg_xxxxyy_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 46); 

                auto tg_xxxxyy_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 47); 

                auto tg_xxxxyy_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 48); 

                auto tg_xxxxxxx_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx); 

                auto tg_xxxxxxx_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 1); 

                auto tg_xxxxxxx_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 2); 

                auto tg_xxxxxxx_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 3); 

                auto tg_xxxxxxx_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 4); 

                auto tg_xxxxxxx_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 5); 

                auto tg_xxxxxxx_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 6); 

                auto tg_xxxxxxx_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 7); 

                auto tg_xxxxxxx_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 8); 

                auto tg_xxxxxxx_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 9); 

                auto tg_xxxxxxy_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 10); 

                auto tg_xxxxxxy_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 11); 

                auto tg_xxxxxxy_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 12); 

                auto tg_xxxxxxy_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 13); 

                auto tg_xxxxxxy_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 14); 

                auto tg_xxxxxxy_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 15); 

                auto tg_xxxxxxy_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 16); 

                auto tg_xxxxxxy_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 17); 

                auto tg_xxxxxxy_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 18); 

                auto tg_xxxxxxy_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 19); 

                auto tg_xxxxxxz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 20); 

                auto tg_xxxxxxz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 21); 

                auto tg_xxxxxxz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 22); 

                auto tg_xxxxxxz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 23); 

                auto tg_xxxxxxz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 24); 

                auto tg_xxxxxxz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 25); 

                auto tg_xxxxxxz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 26); 

                auto tg_xxxxxxz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 27); 

                auto tg_xxxxxxz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 28); 

                auto tg_xxxxxxz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 29); 

                auto tg_xxxxxyy_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 30); 

                auto tg_xxxxxyy_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 31); 

                auto tg_xxxxxyy_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 32); 

                auto tg_xxxxxyy_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 33); 

                // set up pointers to integrals

                auto tg_xxxxxxxx_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx); 

                auto tg_xxxxxxxx_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 1); 

                auto tg_xxxxxxxx_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 2); 

                auto tg_xxxxxxxx_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 3); 

                auto tg_xxxxxxxx_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 4); 

                auto tg_xxxxxxxx_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 5); 

                auto tg_xxxxxxxx_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 6); 

                auto tg_xxxxxxxx_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 7); 

                auto tg_xxxxxxxx_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 8); 

                auto tg_xxxxxxxx_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 9); 

                auto tg_xxxxxxxx_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 10); 

                auto tg_xxxxxxxx_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 11); 

                auto tg_xxxxxxxx_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 12); 

                auto tg_xxxxxxxx_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 13); 

                auto tg_xxxxxxxx_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 14); 

                auto tg_xxxxxxxy_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 15); 

                auto tg_xxxxxxxy_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 16); 

                auto tg_xxxxxxxy_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 17); 

                auto tg_xxxxxxxy_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 18); 

                auto tg_xxxxxxxy_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 19); 

                auto tg_xxxxxxxy_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 20); 

                auto tg_xxxxxxxy_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 21); 

                auto tg_xxxxxxxy_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 22); 

                auto tg_xxxxxxxy_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 23); 

                auto tg_xxxxxxxy_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 24); 

                auto tg_xxxxxxxy_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 25); 

                auto tg_xxxxxxxy_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 26); 

                auto tg_xxxxxxxy_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 27); 

                auto tg_xxxxxxxy_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 28); 

                auto tg_xxxxxxxy_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 29); 

                auto tg_xxxxxxxz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 30); 

                auto tg_xxxxxxxz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 31); 

                auto tg_xxxxxxxz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 32); 

                auto tg_xxxxxxxz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 33); 

                auto tg_xxxxxxxz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 34); 

                auto tg_xxxxxxxz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 35); 

                auto tg_xxxxxxxz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 36); 

                auto tg_xxxxxxxz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 37); 

                auto tg_xxxxxxxz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 38); 

                auto tg_xxxxxxxz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 39); 

                auto tg_xxxxxxxz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 40); 

                auto tg_xxxxxxxz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 41); 

                auto tg_xxxxxxxz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 42); 

                auto tg_xxxxxxxz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 43); 

                auto tg_xxxxxxxz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 44); 

                auto tg_xxxxxxyy_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 45); 

                auto tg_xxxxxxyy_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 46); 

                auto tg_xxxxxxyy_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 47); 

                auto tg_xxxxxxyy_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 48); 

                // Batch of Integrals (0,49)

                #pragma omp simd aligned(fxn, fza, tg_xxxxxx_xxxx_0, tg_xxxxxx_xxxx_1, tg_xxxxxx_xxxy_0, \
                                         tg_xxxxxx_xxxy_1, tg_xxxxxx_xxxz_0, tg_xxxxxx_xxxz_1, tg_xxxxxx_xxyy_0, \
                                         tg_xxxxxx_xxyy_1, tg_xxxxxx_xxyz_0, tg_xxxxxx_xxyz_1, tg_xxxxxx_xxzz_0, \
                                         tg_xxxxxx_xxzz_1, tg_xxxxxx_xyyy_0, tg_xxxxxx_xyyy_1, tg_xxxxxx_xyyz_0, \
                                         tg_xxxxxx_xyyz_1, tg_xxxxxx_xyzz_0, tg_xxxxxx_xyzz_1, tg_xxxxxx_xzzz_0, \
                                         tg_xxxxxx_xzzz_1, tg_xxxxxx_yyyy_0, tg_xxxxxx_yyyy_1, tg_xxxxxx_yyyz_0, \
                                         tg_xxxxxx_yyyz_1, tg_xxxxxx_yyzz_0, tg_xxxxxx_yyzz_1, tg_xxxxxx_yzzz_0, \
                                         tg_xxxxxx_yzzz_1, tg_xxxxxx_zzzz_0, tg_xxxxxx_zzzz_1, tg_xxxxxxx_xxx_1, \
                                         tg_xxxxxxx_xxxx_0, tg_xxxxxxx_xxxx_1, tg_xxxxxxx_xxxy_0, tg_xxxxxxx_xxxy_1, \
                                         tg_xxxxxxx_xxxz_0, tg_xxxxxxx_xxxz_1, tg_xxxxxxx_xxy_1, tg_xxxxxxx_xxyy_0, \
                                         tg_xxxxxxx_xxyy_1, tg_xxxxxxx_xxyz_0, tg_xxxxxxx_xxyz_1, tg_xxxxxxx_xxz_1, \
                                         tg_xxxxxxx_xxzz_0, tg_xxxxxxx_xxzz_1, tg_xxxxxxx_xyy_1, tg_xxxxxxx_xyyy_0, \
                                         tg_xxxxxxx_xyyy_1, tg_xxxxxxx_xyyz_0, tg_xxxxxxx_xyyz_1, tg_xxxxxxx_xyz_1, \
                                         tg_xxxxxxx_xyzz_0, tg_xxxxxxx_xyzz_1, tg_xxxxxxx_xzz_1, tg_xxxxxxx_xzzz_0, \
                                         tg_xxxxxxx_xzzz_1, tg_xxxxxxx_yyy_1, tg_xxxxxxx_yyyy_0, tg_xxxxxxx_yyyy_1, \
                                         tg_xxxxxxx_yyyz_0, tg_xxxxxxx_yyyz_1, tg_xxxxxxx_yyz_1, tg_xxxxxxx_yyzz_0, \
                                         tg_xxxxxxx_yyzz_1, tg_xxxxxxx_yzz_1, tg_xxxxxxx_yzzz_0, tg_xxxxxxx_yzzz_1, \
                                         tg_xxxxxxx_zzz_1, tg_xxxxxxx_zzzz_0, tg_xxxxxxx_zzzz_1, tg_xxxxxxxx_xxxx_0, \
                                         tg_xxxxxxxx_xxxy_0, tg_xxxxxxxx_xxxz_0, tg_xxxxxxxx_xxyy_0, tg_xxxxxxxx_xxyz_0, \
                                         tg_xxxxxxxx_xxzz_0, tg_xxxxxxxx_xyyy_0, tg_xxxxxxxx_xyyz_0, tg_xxxxxxxx_xyzz_0, \
                                         tg_xxxxxxxx_xzzz_0, tg_xxxxxxxx_yyyy_0, tg_xxxxxxxx_yyyz_0, tg_xxxxxxxx_yyzz_0, \
                                         tg_xxxxxxxx_yzzz_0, tg_xxxxxxxx_zzzz_0, tg_xxxxxxxy_xxxx_0, tg_xxxxxxxy_xxxy_0, \
                                         tg_xxxxxxxy_xxxz_0, tg_xxxxxxxy_xxyy_0, tg_xxxxxxxy_xxyz_0, tg_xxxxxxxy_xxzz_0, \
                                         tg_xxxxxxxy_xyyy_0, tg_xxxxxxxy_xyyz_0, tg_xxxxxxxy_xyzz_0, tg_xxxxxxxy_xzzz_0, \
                                         tg_xxxxxxxy_yyyy_0, tg_xxxxxxxy_yyyz_0, tg_xxxxxxxy_yyzz_0, tg_xxxxxxxy_yzzz_0, \
                                         tg_xxxxxxxy_zzzz_0, tg_xxxxxxxz_xxxx_0, tg_xxxxxxxz_xxxy_0, tg_xxxxxxxz_xxxz_0, \
                                         tg_xxxxxxxz_xxyy_0, tg_xxxxxxxz_xxyz_0, tg_xxxxxxxz_xxzz_0, tg_xxxxxxxz_xyyy_0, \
                                         tg_xxxxxxxz_xyyz_0, tg_xxxxxxxz_xyzz_0, tg_xxxxxxxz_xzzz_0, tg_xxxxxxxz_yyyy_0, \
                                         tg_xxxxxxxz_yyyz_0, tg_xxxxxxxz_yyzz_0, tg_xxxxxxxz_yzzz_0, tg_xxxxxxxz_zzzz_0, \
                                         tg_xxxxxxy_xxx_1, tg_xxxxxxy_xxxx_0, tg_xxxxxxy_xxxx_1, tg_xxxxxxy_xxxy_0, \
                                         tg_xxxxxxy_xxxy_1, tg_xxxxxxy_xxxz_0, tg_xxxxxxy_xxxz_1, tg_xxxxxxy_xxy_1, \
                                         tg_xxxxxxy_xxyy_0, tg_xxxxxxy_xxyy_1, tg_xxxxxxy_xxyz_0, tg_xxxxxxy_xxyz_1, \
                                         tg_xxxxxxy_xxz_1, tg_xxxxxxy_xxzz_0, tg_xxxxxxy_xxzz_1, tg_xxxxxxy_xyy_1, \
                                         tg_xxxxxxy_xyyy_0, tg_xxxxxxy_xyyy_1, tg_xxxxxxy_xyyz_0, tg_xxxxxxy_xyyz_1, \
                                         tg_xxxxxxy_xyz_1, tg_xxxxxxy_xyzz_0, tg_xxxxxxy_xyzz_1, tg_xxxxxxy_xzz_1, \
                                         tg_xxxxxxy_xzzz_0, tg_xxxxxxy_xzzz_1, tg_xxxxxxy_yyy_1, tg_xxxxxxy_yyyy_0, \
                                         tg_xxxxxxy_yyyy_1, tg_xxxxxxy_yyyz_0, tg_xxxxxxy_yyyz_1, tg_xxxxxxy_yyz_1, \
                                         tg_xxxxxxy_yyzz_0, tg_xxxxxxy_yyzz_1, tg_xxxxxxy_yzz_1, tg_xxxxxxy_yzzz_0, \
                                         tg_xxxxxxy_yzzz_1, tg_xxxxxxy_zzz_1, tg_xxxxxxy_zzzz_0, tg_xxxxxxy_zzzz_1, \
                                         tg_xxxxxxyy_xxxx_0, tg_xxxxxxyy_xxxy_0, tg_xxxxxxyy_xxxz_0, tg_xxxxxxyy_xxyy_0, \
                                         tg_xxxxxxz_xxx_1, tg_xxxxxxz_xxxx_0, tg_xxxxxxz_xxxx_1, tg_xxxxxxz_xxxy_0, \
                                         tg_xxxxxxz_xxxy_1, tg_xxxxxxz_xxxz_0, tg_xxxxxxz_xxxz_1, tg_xxxxxxz_xxy_1, \
                                         tg_xxxxxxz_xxyy_0, tg_xxxxxxz_xxyy_1, tg_xxxxxxz_xxyz_0, tg_xxxxxxz_xxyz_1, \
                                         tg_xxxxxxz_xxz_1, tg_xxxxxxz_xxzz_0, tg_xxxxxxz_xxzz_1, tg_xxxxxxz_xyy_1, \
                                         tg_xxxxxxz_xyyy_0, tg_xxxxxxz_xyyy_1, tg_xxxxxxz_xyyz_0, tg_xxxxxxz_xyyz_1, \
                                         tg_xxxxxxz_xyz_1, tg_xxxxxxz_xyzz_0, tg_xxxxxxz_xyzz_1, tg_xxxxxxz_xzz_1, \
                                         tg_xxxxxxz_xzzz_0, tg_xxxxxxz_xzzz_1, tg_xxxxxxz_yyy_1, tg_xxxxxxz_yyyy_0, \
                                         tg_xxxxxxz_yyyy_1, tg_xxxxxxz_yyyz_0, tg_xxxxxxz_yyyz_1, tg_xxxxxxz_yyz_1, \
                                         tg_xxxxxxz_yyzz_0, tg_xxxxxxz_yyzz_1, tg_xxxxxxz_yzz_1, tg_xxxxxxz_yzzz_0, \
                                         tg_xxxxxxz_yzzz_1, tg_xxxxxxz_zzz_1, tg_xxxxxxz_zzzz_0, tg_xxxxxxz_zzzz_1, \
                                         tg_xxxxxy_xxxx_0, tg_xxxxxy_xxxx_1, tg_xxxxxy_xxxy_0, tg_xxxxxy_xxxy_1, \
                                         tg_xxxxxy_xxxz_0, tg_xxxxxy_xxxz_1, tg_xxxxxy_xxyy_0, tg_xxxxxy_xxyy_1, \
                                         tg_xxxxxy_xxyz_0, tg_xxxxxy_xxyz_1, tg_xxxxxy_xxzz_0, tg_xxxxxy_xxzz_1, \
                                         tg_xxxxxy_xyyy_0, tg_xxxxxy_xyyy_1, tg_xxxxxy_xyyz_0, tg_xxxxxy_xyyz_1, \
                                         tg_xxxxxy_xyzz_0, tg_xxxxxy_xyzz_1, tg_xxxxxy_xzzz_0, tg_xxxxxy_xzzz_1, \
                                         tg_xxxxxy_yyyy_0, tg_xxxxxy_yyyy_1, tg_xxxxxy_yyyz_0, tg_xxxxxy_yyyz_1, \
                                         tg_xxxxxy_yyzz_0, tg_xxxxxy_yyzz_1, tg_xxxxxy_yzzz_0, tg_xxxxxy_yzzz_1, \
                                         tg_xxxxxy_zzzz_0, tg_xxxxxy_zzzz_1, tg_xxxxxyy_xxx_1, tg_xxxxxyy_xxxx_0, \
                                         tg_xxxxxyy_xxxx_1, tg_xxxxxyy_xxxy_0, tg_xxxxxyy_xxxy_1, tg_xxxxxyy_xxxz_0, \
                                         tg_xxxxxyy_xxxz_1, tg_xxxxxyy_xxy_1, tg_xxxxxyy_xxyy_0, tg_xxxxxyy_xxyy_1, \
                                         tg_xxxxxyy_xxz_1, tg_xxxxxyy_xyy_1, tg_xxxxxz_xxxx_0, tg_xxxxxz_xxxx_1, \
                                         tg_xxxxxz_xxxy_0, tg_xxxxxz_xxxy_1, tg_xxxxxz_xxxz_0, tg_xxxxxz_xxxz_1, \
                                         tg_xxxxxz_xxyy_0, tg_xxxxxz_xxyy_1, tg_xxxxxz_xxyz_0, tg_xxxxxz_xxyz_1, \
                                         tg_xxxxxz_xxzz_0, tg_xxxxxz_xxzz_1, tg_xxxxxz_xyyy_0, tg_xxxxxz_xyyy_1, \
                                         tg_xxxxxz_xyyz_0, tg_xxxxxz_xyyz_1, tg_xxxxxz_xyzz_0, tg_xxxxxz_xyzz_1, \
                                         tg_xxxxxz_xzzz_0, tg_xxxxxz_xzzz_1, tg_xxxxxz_yyyy_0, tg_xxxxxz_yyyy_1, \
                                         tg_xxxxxz_yyyz_0, tg_xxxxxz_yyyz_1, tg_xxxxxz_yyzz_0, tg_xxxxxz_yyzz_1, \
                                         tg_xxxxxz_yzzz_0, tg_xxxxxz_yzzz_1, tg_xxxxxz_zzzz_0, tg_xxxxxz_zzzz_1, \
                                         tg_xxxxyy_xxxx_0, tg_xxxxyy_xxxx_1, tg_xxxxyy_xxxy_0, tg_xxxxyy_xxxy_1, \
                                         tg_xxxxyy_xxxz_0, tg_xxxxyy_xxxz_1, tg_xxxxyy_xxyy_0, tg_xxxxyy_xxyy_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxxxxxxx_xxxx_0[j] = pb_x * tg_xxxxxxx_xxxx_0[j] + wp_x[j] * tg_xxxxxxx_xxxx_1[j] + 3.5 * fl1_fx * tg_xxxxxx_xxxx_0[j] - 3.5 * fl1_fx * fl1_fza * tg_xxxxxx_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxxxxxx_xxx_1[j];

                    tg_xxxxxxxx_xxxy_0[j] = pb_x * tg_xxxxxxx_xxxy_0[j] + wp_x[j] * tg_xxxxxxx_xxxy_1[j] + 3.5 * fl1_fx * tg_xxxxxx_xxxy_0[j] - 3.5 * fl1_fx * fl1_fza * tg_xxxxxx_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxxxxxx_xxy_1[j];

                    tg_xxxxxxxx_xxxz_0[j] = pb_x * tg_xxxxxxx_xxxz_0[j] + wp_x[j] * tg_xxxxxxx_xxxz_1[j] + 3.5 * fl1_fx * tg_xxxxxx_xxxz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_xxxxxx_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxxxxxx_xxz_1[j];

                    tg_xxxxxxxx_xxyy_0[j] = pb_x * tg_xxxxxxx_xxyy_0[j] + wp_x[j] * tg_xxxxxxx_xxyy_1[j] + 3.5 * fl1_fx * tg_xxxxxx_xxyy_0[j] - 3.5 * fl1_fx * fl1_fza * tg_xxxxxx_xxyy_1[j] + fl1_fxn * tg_xxxxxxx_xyy_1[j];

                    tg_xxxxxxxx_xxyz_0[j] = pb_x * tg_xxxxxxx_xxyz_0[j] + wp_x[j] * tg_xxxxxxx_xxyz_1[j] + 3.5 * fl1_fx * tg_xxxxxx_xxyz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_xxxxxx_xxyz_1[j] + fl1_fxn * tg_xxxxxxx_xyz_1[j];

                    tg_xxxxxxxx_xxzz_0[j] = pb_x * tg_xxxxxxx_xxzz_0[j] + wp_x[j] * tg_xxxxxxx_xxzz_1[j] + 3.5 * fl1_fx * tg_xxxxxx_xxzz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_xxxxxx_xxzz_1[j] + fl1_fxn * tg_xxxxxxx_xzz_1[j];

                    tg_xxxxxxxx_xyyy_0[j] = pb_x * tg_xxxxxxx_xyyy_0[j] + wp_x[j] * tg_xxxxxxx_xyyy_1[j] + 3.5 * fl1_fx * tg_xxxxxx_xyyy_0[j] - 3.5 * fl1_fx * fl1_fza * tg_xxxxxx_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxxxxxx_yyy_1[j];

                    tg_xxxxxxxx_xyyz_0[j] = pb_x * tg_xxxxxxx_xyyz_0[j] + wp_x[j] * tg_xxxxxxx_xyyz_1[j] + 3.5 * fl1_fx * tg_xxxxxx_xyyz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_xxxxxx_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxxxxxx_yyz_1[j];

                    tg_xxxxxxxx_xyzz_0[j] = pb_x * tg_xxxxxxx_xyzz_0[j] + wp_x[j] * tg_xxxxxxx_xyzz_1[j] + 3.5 * fl1_fx * tg_xxxxxx_xyzz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_xxxxxx_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxxxxxx_yzz_1[j];

                    tg_xxxxxxxx_xzzz_0[j] = pb_x * tg_xxxxxxx_xzzz_0[j] + wp_x[j] * tg_xxxxxxx_xzzz_1[j] + 3.5 * fl1_fx * tg_xxxxxx_xzzz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_xxxxxx_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxxxxxx_zzz_1[j];

                    tg_xxxxxxxx_yyyy_0[j] = pb_x * tg_xxxxxxx_yyyy_0[j] + wp_x[j] * tg_xxxxxxx_yyyy_1[j] + 3.5 * fl1_fx * tg_xxxxxx_yyyy_0[j] - 3.5 * fl1_fx * fl1_fza * tg_xxxxxx_yyyy_1[j];

                    tg_xxxxxxxx_yyyz_0[j] = pb_x * tg_xxxxxxx_yyyz_0[j] + wp_x[j] * tg_xxxxxxx_yyyz_1[j] + 3.5 * fl1_fx * tg_xxxxxx_yyyz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_xxxxxx_yyyz_1[j];

                    tg_xxxxxxxx_yyzz_0[j] = pb_x * tg_xxxxxxx_yyzz_0[j] + wp_x[j] * tg_xxxxxxx_yyzz_1[j] + 3.5 * fl1_fx * tg_xxxxxx_yyzz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_xxxxxx_yyzz_1[j];

                    tg_xxxxxxxx_yzzz_0[j] = pb_x * tg_xxxxxxx_yzzz_0[j] + wp_x[j] * tg_xxxxxxx_yzzz_1[j] + 3.5 * fl1_fx * tg_xxxxxx_yzzz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_xxxxxx_yzzz_1[j];

                    tg_xxxxxxxx_zzzz_0[j] = pb_x * tg_xxxxxxx_zzzz_0[j] + wp_x[j] * tg_xxxxxxx_zzzz_1[j] + 3.5 * fl1_fx * tg_xxxxxx_zzzz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_xxxxxx_zzzz_1[j];

                    tg_xxxxxxxy_xxxx_0[j] = pb_x * tg_xxxxxxy_xxxx_0[j] + wp_x[j] * tg_xxxxxxy_xxxx_1[j] + 3.0 * fl1_fx * tg_xxxxxy_xxxx_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxy_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxxxxxy_xxx_1[j];

                    tg_xxxxxxxy_xxxy_0[j] = pb_x * tg_xxxxxxy_xxxy_0[j] + wp_x[j] * tg_xxxxxxy_xxxy_1[j] + 3.0 * fl1_fx * tg_xxxxxy_xxxy_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxy_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxxxxxy_xxy_1[j];

                    tg_xxxxxxxy_xxxz_0[j] = pb_x * tg_xxxxxxy_xxxz_0[j] + wp_x[j] * tg_xxxxxxy_xxxz_1[j] + 3.0 * fl1_fx * tg_xxxxxy_xxxz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxy_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxxxxxy_xxz_1[j];

                    tg_xxxxxxxy_xxyy_0[j] = pb_x * tg_xxxxxxy_xxyy_0[j] + wp_x[j] * tg_xxxxxxy_xxyy_1[j] + 3.0 * fl1_fx * tg_xxxxxy_xxyy_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxy_xxyy_1[j] + fl1_fxn * tg_xxxxxxy_xyy_1[j];

                    tg_xxxxxxxy_xxyz_0[j] = pb_x * tg_xxxxxxy_xxyz_0[j] + wp_x[j] * tg_xxxxxxy_xxyz_1[j] + 3.0 * fl1_fx * tg_xxxxxy_xxyz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxy_xxyz_1[j] + fl1_fxn * tg_xxxxxxy_xyz_1[j];

                    tg_xxxxxxxy_xxzz_0[j] = pb_x * tg_xxxxxxy_xxzz_0[j] + wp_x[j] * tg_xxxxxxy_xxzz_1[j] + 3.0 * fl1_fx * tg_xxxxxy_xxzz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxy_xxzz_1[j] + fl1_fxn * tg_xxxxxxy_xzz_1[j];

                    tg_xxxxxxxy_xyyy_0[j] = pb_x * tg_xxxxxxy_xyyy_0[j] + wp_x[j] * tg_xxxxxxy_xyyy_1[j] + 3.0 * fl1_fx * tg_xxxxxy_xyyy_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxy_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxxxxxy_yyy_1[j];

                    tg_xxxxxxxy_xyyz_0[j] = pb_x * tg_xxxxxxy_xyyz_0[j] + wp_x[j] * tg_xxxxxxy_xyyz_1[j] + 3.0 * fl1_fx * tg_xxxxxy_xyyz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxy_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxxxxxy_yyz_1[j];

                    tg_xxxxxxxy_xyzz_0[j] = pb_x * tg_xxxxxxy_xyzz_0[j] + wp_x[j] * tg_xxxxxxy_xyzz_1[j] + 3.0 * fl1_fx * tg_xxxxxy_xyzz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxy_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxxxxxy_yzz_1[j];

                    tg_xxxxxxxy_xzzz_0[j] = pb_x * tg_xxxxxxy_xzzz_0[j] + wp_x[j] * tg_xxxxxxy_xzzz_1[j] + 3.0 * fl1_fx * tg_xxxxxy_xzzz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxy_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxxxxxy_zzz_1[j];

                    tg_xxxxxxxy_yyyy_0[j] = pb_x * tg_xxxxxxy_yyyy_0[j] + wp_x[j] * tg_xxxxxxy_yyyy_1[j] + 3.0 * fl1_fx * tg_xxxxxy_yyyy_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxy_yyyy_1[j];

                    tg_xxxxxxxy_yyyz_0[j] = pb_x * tg_xxxxxxy_yyyz_0[j] + wp_x[j] * tg_xxxxxxy_yyyz_1[j] + 3.0 * fl1_fx * tg_xxxxxy_yyyz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxy_yyyz_1[j];

                    tg_xxxxxxxy_yyzz_0[j] = pb_x * tg_xxxxxxy_yyzz_0[j] + wp_x[j] * tg_xxxxxxy_yyzz_1[j] + 3.0 * fl1_fx * tg_xxxxxy_yyzz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxy_yyzz_1[j];

                    tg_xxxxxxxy_yzzz_0[j] = pb_x * tg_xxxxxxy_yzzz_0[j] + wp_x[j] * tg_xxxxxxy_yzzz_1[j] + 3.0 * fl1_fx * tg_xxxxxy_yzzz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxy_yzzz_1[j];

                    tg_xxxxxxxy_zzzz_0[j] = pb_x * tg_xxxxxxy_zzzz_0[j] + wp_x[j] * tg_xxxxxxy_zzzz_1[j] + 3.0 * fl1_fx * tg_xxxxxy_zzzz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxy_zzzz_1[j];

                    tg_xxxxxxxz_xxxx_0[j] = pb_x * tg_xxxxxxz_xxxx_0[j] + wp_x[j] * tg_xxxxxxz_xxxx_1[j] + 3.0 * fl1_fx * tg_xxxxxz_xxxx_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxxxxxz_xxx_1[j];

                    tg_xxxxxxxz_xxxy_0[j] = pb_x * tg_xxxxxxz_xxxy_0[j] + wp_x[j] * tg_xxxxxxz_xxxy_1[j] + 3.0 * fl1_fx * tg_xxxxxz_xxxy_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxxxxxz_xxy_1[j];

                    tg_xxxxxxxz_xxxz_0[j] = pb_x * tg_xxxxxxz_xxxz_0[j] + wp_x[j] * tg_xxxxxxz_xxxz_1[j] + 3.0 * fl1_fx * tg_xxxxxz_xxxz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxxxxxz_xxz_1[j];

                    tg_xxxxxxxz_xxyy_0[j] = pb_x * tg_xxxxxxz_xxyy_0[j] + wp_x[j] * tg_xxxxxxz_xxyy_1[j] + 3.0 * fl1_fx * tg_xxxxxz_xxyy_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxz_xxyy_1[j] + fl1_fxn * tg_xxxxxxz_xyy_1[j];

                    tg_xxxxxxxz_xxyz_0[j] = pb_x * tg_xxxxxxz_xxyz_0[j] + wp_x[j] * tg_xxxxxxz_xxyz_1[j] + 3.0 * fl1_fx * tg_xxxxxz_xxyz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxz_xxyz_1[j] + fl1_fxn * tg_xxxxxxz_xyz_1[j];

                    tg_xxxxxxxz_xxzz_0[j] = pb_x * tg_xxxxxxz_xxzz_0[j] + wp_x[j] * tg_xxxxxxz_xxzz_1[j] + 3.0 * fl1_fx * tg_xxxxxz_xxzz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxz_xxzz_1[j] + fl1_fxn * tg_xxxxxxz_xzz_1[j];

                    tg_xxxxxxxz_xyyy_0[j] = pb_x * tg_xxxxxxz_xyyy_0[j] + wp_x[j] * tg_xxxxxxz_xyyy_1[j] + 3.0 * fl1_fx * tg_xxxxxz_xyyy_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxxxxxz_yyy_1[j];

                    tg_xxxxxxxz_xyyz_0[j] = pb_x * tg_xxxxxxz_xyyz_0[j] + wp_x[j] * tg_xxxxxxz_xyyz_1[j] + 3.0 * fl1_fx * tg_xxxxxz_xyyz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxxxxxz_yyz_1[j];

                    tg_xxxxxxxz_xyzz_0[j] = pb_x * tg_xxxxxxz_xyzz_0[j] + wp_x[j] * tg_xxxxxxz_xyzz_1[j] + 3.0 * fl1_fx * tg_xxxxxz_xyzz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxxxxxz_yzz_1[j];

                    tg_xxxxxxxz_xzzz_0[j] = pb_x * tg_xxxxxxz_xzzz_0[j] + wp_x[j] * tg_xxxxxxz_xzzz_1[j] + 3.0 * fl1_fx * tg_xxxxxz_xzzz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxxxxxz_zzz_1[j];

                    tg_xxxxxxxz_yyyy_0[j] = pb_x * tg_xxxxxxz_yyyy_0[j] + wp_x[j] * tg_xxxxxxz_yyyy_1[j] + 3.0 * fl1_fx * tg_xxxxxz_yyyy_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxz_yyyy_1[j];

                    tg_xxxxxxxz_yyyz_0[j] = pb_x * tg_xxxxxxz_yyyz_0[j] + wp_x[j] * tg_xxxxxxz_yyyz_1[j] + 3.0 * fl1_fx * tg_xxxxxz_yyyz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxz_yyyz_1[j];

                    tg_xxxxxxxz_yyzz_0[j] = pb_x * tg_xxxxxxz_yyzz_0[j] + wp_x[j] * tg_xxxxxxz_yyzz_1[j] + 3.0 * fl1_fx * tg_xxxxxz_yyzz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxz_yyzz_1[j];

                    tg_xxxxxxxz_yzzz_0[j] = pb_x * tg_xxxxxxz_yzzz_0[j] + wp_x[j] * tg_xxxxxxz_yzzz_1[j] + 3.0 * fl1_fx * tg_xxxxxz_yzzz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxz_yzzz_1[j];

                    tg_xxxxxxxz_zzzz_0[j] = pb_x * tg_xxxxxxz_zzzz_0[j] + wp_x[j] * tg_xxxxxxz_zzzz_1[j] + 3.0 * fl1_fx * tg_xxxxxz_zzzz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxz_zzzz_1[j];

                    tg_xxxxxxyy_xxxx_0[j] = pb_x * tg_xxxxxyy_xxxx_0[j] + wp_x[j] * tg_xxxxxyy_xxxx_1[j] + 2.5 * fl1_fx * tg_xxxxyy_xxxx_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyy_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxxxxyy_xxx_1[j];

                    tg_xxxxxxyy_xxxy_0[j] = pb_x * tg_xxxxxyy_xxxy_0[j] + wp_x[j] * tg_xxxxxyy_xxxy_1[j] + 2.5 * fl1_fx * tg_xxxxyy_xxxy_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyy_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxxxxyy_xxy_1[j];

                    tg_xxxxxxyy_xxxz_0[j] = pb_x * tg_xxxxxyy_xxxz_0[j] + wp_x[j] * tg_xxxxxyy_xxxz_1[j] + 2.5 * fl1_fx * tg_xxxxyy_xxxz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyy_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxxxxyy_xxz_1[j];

                    tg_xxxxxxyy_xxyy_0[j] = pb_x * tg_xxxxxyy_xxyy_0[j] + wp_x[j] * tg_xxxxxyy_xxyy_1[j] + 2.5 * fl1_fx * tg_xxxxyy_xxyy_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyy_xxyy_1[j] + fl1_fxn * tg_xxxxxyy_xyy_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSLSG_49_98(      CMemBlock2D<double>& primBuffer,
                                       const CRecursionMap&       recursionMap,
                                       const CMemBlock2D<double>& osFactors,
                                       const CMemBlock2D<double>& wpDistances,
                                       const CGtoPairsBlock&      braGtoPairsBlock,
                                       const CGtoPairsBlock&      ketGtoPairsBlock,
                                       const int32_t              nKetPrimPairs,
                                       const int32_t              iContrPair)
    {
        // Batch of Integrals (49,98)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {8, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {8, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_7_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_6_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_6_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_7_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_xxxxxyy_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 49); 

                auto tg_xxxxxyy_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 50); 

                auto tg_xxxxxyy_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 51); 

                auto tg_xxxxxyy_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 52); 

                auto tg_xxxxxyy_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 53); 

                auto tg_xxxxxyy_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 54); 

                auto tg_xxxxxyy_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 55); 

                auto tg_xxxxxyy_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 56); 

                auto tg_xxxxxyy_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 57); 

                auto tg_xxxxxyy_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 58); 

                auto tg_xxxxxyy_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 59); 

                auto tg_xxxxxyz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 60); 

                auto tg_xxxxxyz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 61); 

                auto tg_xxxxxyz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 62); 

                auto tg_xxxxxyz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 63); 

                auto tg_xxxxxyz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 64); 

                auto tg_xxxxxyz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 65); 

                auto tg_xxxxxyz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 66); 

                auto tg_xxxxxyz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 67); 

                auto tg_xxxxxyz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 68); 

                auto tg_xxxxxyz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 69); 

                auto tg_xxxxxyz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 70); 

                auto tg_xxxxxyz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 71); 

                auto tg_xxxxxyz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 72); 

                auto tg_xxxxxyz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 73); 

                auto tg_xxxxxyz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 74); 

                auto tg_xxxxxzz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 75); 

                auto tg_xxxxxzz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 76); 

                auto tg_xxxxxzz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 77); 

                auto tg_xxxxxzz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 78); 

                auto tg_xxxxxzz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 79); 

                auto tg_xxxxxzz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 80); 

                auto tg_xxxxxzz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 81); 

                auto tg_xxxxxzz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 82); 

                auto tg_xxxxxzz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 83); 

                auto tg_xxxxxzz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 84); 

                auto tg_xxxxxzz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 85); 

                auto tg_xxxxxzz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 86); 

                auto tg_xxxxxzz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 87); 

                auto tg_xxxxxzz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 88); 

                auto tg_xxxxxzz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 89); 

                auto tg_xxxxyyy_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 90); 

                auto tg_xxxxyyy_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 91); 

                auto tg_xxxxyyy_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 92); 

                auto tg_xxxxyyy_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 93); 

                auto tg_xxxxyyy_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 94); 

                auto tg_xxxxyyy_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 95); 

                auto tg_xxxxyyy_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 96); 

                auto tg_xxxxyyy_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 97); 

                auto tg_xxxxxyy_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 49); 

                auto tg_xxxxxyy_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 50); 

                auto tg_xxxxxyy_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 51); 

                auto tg_xxxxxyy_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 52); 

                auto tg_xxxxxyy_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 53); 

                auto tg_xxxxxyy_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 54); 

                auto tg_xxxxxyy_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 55); 

                auto tg_xxxxxyy_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 56); 

                auto tg_xxxxxyy_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 57); 

                auto tg_xxxxxyy_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 58); 

                auto tg_xxxxxyy_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 59); 

                auto tg_xxxxxyz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 60); 

                auto tg_xxxxxyz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 61); 

                auto tg_xxxxxyz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 62); 

                auto tg_xxxxxyz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 63); 

                auto tg_xxxxxyz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 64); 

                auto tg_xxxxxyz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 65); 

                auto tg_xxxxxyz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 66); 

                auto tg_xxxxxyz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 67); 

                auto tg_xxxxxyz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 68); 

                auto tg_xxxxxyz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 69); 

                auto tg_xxxxxyz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 70); 

                auto tg_xxxxxyz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 71); 

                auto tg_xxxxxyz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 72); 

                auto tg_xxxxxyz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 73); 

                auto tg_xxxxxyz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 74); 

                auto tg_xxxxxzz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 75); 

                auto tg_xxxxxzz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 76); 

                auto tg_xxxxxzz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 77); 

                auto tg_xxxxxzz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 78); 

                auto tg_xxxxxzz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 79); 

                auto tg_xxxxxzz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 80); 

                auto tg_xxxxxzz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 81); 

                auto tg_xxxxxzz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 82); 

                auto tg_xxxxxzz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 83); 

                auto tg_xxxxxzz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 84); 

                auto tg_xxxxxzz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 85); 

                auto tg_xxxxxzz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 86); 

                auto tg_xxxxxzz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 87); 

                auto tg_xxxxxzz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 88); 

                auto tg_xxxxxzz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 89); 

                auto tg_xxxxyyy_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 90); 

                auto tg_xxxxyyy_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 91); 

                auto tg_xxxxyyy_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 92); 

                auto tg_xxxxyyy_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 93); 

                auto tg_xxxxyyy_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 94); 

                auto tg_xxxxyyy_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 95); 

                auto tg_xxxxyyy_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 96); 

                auto tg_xxxxyyy_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 97); 

                auto tg_xxxxyy_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 49); 

                auto tg_xxxxyy_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 50); 

                auto tg_xxxxyy_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 51); 

                auto tg_xxxxyy_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 52); 

                auto tg_xxxxyy_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 53); 

                auto tg_xxxxyy_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 54); 

                auto tg_xxxxyy_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 55); 

                auto tg_xxxxyy_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 56); 

                auto tg_xxxxyy_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 57); 

                auto tg_xxxxyy_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 58); 

                auto tg_xxxxyy_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 59); 

                auto tg_xxxxyz_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 60); 

                auto tg_xxxxyz_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 61); 

                auto tg_xxxxyz_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 62); 

                auto tg_xxxxyz_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 63); 

                auto tg_xxxxyz_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 64); 

                auto tg_xxxxyz_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 65); 

                auto tg_xxxxyz_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 66); 

                auto tg_xxxxyz_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 67); 

                auto tg_xxxxyz_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 68); 

                auto tg_xxxxyz_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 69); 

                auto tg_xxxxyz_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 70); 

                auto tg_xxxxyz_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 71); 

                auto tg_xxxxyz_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 72); 

                auto tg_xxxxyz_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 73); 

                auto tg_xxxxyz_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 74); 

                auto tg_xxxxzz_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 75); 

                auto tg_xxxxzz_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 76); 

                auto tg_xxxxzz_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 77); 

                auto tg_xxxxzz_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 78); 

                auto tg_xxxxzz_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 79); 

                auto tg_xxxxzz_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 80); 

                auto tg_xxxxzz_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 81); 

                auto tg_xxxxzz_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 82); 

                auto tg_xxxxzz_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 83); 

                auto tg_xxxxzz_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 84); 

                auto tg_xxxxzz_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 85); 

                auto tg_xxxxzz_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 86); 

                auto tg_xxxxzz_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 87); 

                auto tg_xxxxzz_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 88); 

                auto tg_xxxxzz_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 89); 

                auto tg_xxxyyy_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 90); 

                auto tg_xxxyyy_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 91); 

                auto tg_xxxyyy_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 92); 

                auto tg_xxxyyy_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 93); 

                auto tg_xxxyyy_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 94); 

                auto tg_xxxyyy_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 95); 

                auto tg_xxxyyy_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 96); 

                auto tg_xxxyyy_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 97); 

                auto tg_xxxxyy_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 49); 

                auto tg_xxxxyy_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 50); 

                auto tg_xxxxyy_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 51); 

                auto tg_xxxxyy_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 52); 

                auto tg_xxxxyy_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 53); 

                auto tg_xxxxyy_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 54); 

                auto tg_xxxxyy_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 55); 

                auto tg_xxxxyy_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 56); 

                auto tg_xxxxyy_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 57); 

                auto tg_xxxxyy_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 58); 

                auto tg_xxxxyy_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 59); 

                auto tg_xxxxyz_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 60); 

                auto tg_xxxxyz_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 61); 

                auto tg_xxxxyz_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 62); 

                auto tg_xxxxyz_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 63); 

                auto tg_xxxxyz_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 64); 

                auto tg_xxxxyz_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 65); 

                auto tg_xxxxyz_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 66); 

                auto tg_xxxxyz_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 67); 

                auto tg_xxxxyz_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 68); 

                auto tg_xxxxyz_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 69); 

                auto tg_xxxxyz_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 70); 

                auto tg_xxxxyz_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 71); 

                auto tg_xxxxyz_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 72); 

                auto tg_xxxxyz_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 73); 

                auto tg_xxxxyz_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 74); 

                auto tg_xxxxzz_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 75); 

                auto tg_xxxxzz_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 76); 

                auto tg_xxxxzz_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 77); 

                auto tg_xxxxzz_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 78); 

                auto tg_xxxxzz_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 79); 

                auto tg_xxxxzz_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 80); 

                auto tg_xxxxzz_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 81); 

                auto tg_xxxxzz_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 82); 

                auto tg_xxxxzz_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 83); 

                auto tg_xxxxzz_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 84); 

                auto tg_xxxxzz_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 85); 

                auto tg_xxxxzz_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 86); 

                auto tg_xxxxzz_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 87); 

                auto tg_xxxxzz_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 88); 

                auto tg_xxxxzz_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 89); 

                auto tg_xxxyyy_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 90); 

                auto tg_xxxyyy_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 91); 

                auto tg_xxxyyy_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 92); 

                auto tg_xxxyyy_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 93); 

                auto tg_xxxyyy_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 94); 

                auto tg_xxxyyy_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 95); 

                auto tg_xxxyyy_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 96); 

                auto tg_xxxyyy_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 97); 

                auto tg_xxxxxyy_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 34); 

                auto tg_xxxxxyy_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 35); 

                auto tg_xxxxxyy_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 36); 

                auto tg_xxxxxyy_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 37); 

                auto tg_xxxxxyy_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 38); 

                auto tg_xxxxxyy_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 39); 

                auto tg_xxxxxyz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 40); 

                auto tg_xxxxxyz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 41); 

                auto tg_xxxxxyz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 42); 

                auto tg_xxxxxyz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 43); 

                auto tg_xxxxxyz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 44); 

                auto tg_xxxxxyz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 45); 

                auto tg_xxxxxyz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 46); 

                auto tg_xxxxxyz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 47); 

                auto tg_xxxxxyz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 48); 

                auto tg_xxxxxyz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 49); 

                auto tg_xxxxxzz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 50); 

                auto tg_xxxxxzz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 51); 

                auto tg_xxxxxzz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 52); 

                auto tg_xxxxxzz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 53); 

                auto tg_xxxxxzz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 54); 

                auto tg_xxxxxzz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 55); 

                auto tg_xxxxxzz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 56); 

                auto tg_xxxxxzz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 57); 

                auto tg_xxxxxzz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 58); 

                auto tg_xxxxxzz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 59); 

                auto tg_xxxxyyy_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 60); 

                auto tg_xxxxyyy_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 61); 

                auto tg_xxxxyyy_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 62); 

                auto tg_xxxxyyy_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 63); 

                auto tg_xxxxyyy_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 64); 

                auto tg_xxxxyyy_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 65); 

                auto tg_xxxxyyy_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 66); 

                auto tg_xxxxyyy_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 67); 

                // set up pointers to integrals

                auto tg_xxxxxxyy_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 49); 

                auto tg_xxxxxxyy_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 50); 

                auto tg_xxxxxxyy_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 51); 

                auto tg_xxxxxxyy_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 52); 

                auto tg_xxxxxxyy_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 53); 

                auto tg_xxxxxxyy_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 54); 

                auto tg_xxxxxxyy_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 55); 

                auto tg_xxxxxxyy_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 56); 

                auto tg_xxxxxxyy_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 57); 

                auto tg_xxxxxxyy_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 58); 

                auto tg_xxxxxxyy_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 59); 

                auto tg_xxxxxxyz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 60); 

                auto tg_xxxxxxyz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 61); 

                auto tg_xxxxxxyz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 62); 

                auto tg_xxxxxxyz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 63); 

                auto tg_xxxxxxyz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 64); 

                auto tg_xxxxxxyz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 65); 

                auto tg_xxxxxxyz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 66); 

                auto tg_xxxxxxyz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 67); 

                auto tg_xxxxxxyz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 68); 

                auto tg_xxxxxxyz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 69); 

                auto tg_xxxxxxyz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 70); 

                auto tg_xxxxxxyz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 71); 

                auto tg_xxxxxxyz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 72); 

                auto tg_xxxxxxyz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 73); 

                auto tg_xxxxxxyz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 74); 

                auto tg_xxxxxxzz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 75); 

                auto tg_xxxxxxzz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 76); 

                auto tg_xxxxxxzz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 77); 

                auto tg_xxxxxxzz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 78); 

                auto tg_xxxxxxzz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 79); 

                auto tg_xxxxxxzz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 80); 

                auto tg_xxxxxxzz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 81); 

                auto tg_xxxxxxzz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 82); 

                auto tg_xxxxxxzz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 83); 

                auto tg_xxxxxxzz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 84); 

                auto tg_xxxxxxzz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 85); 

                auto tg_xxxxxxzz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 86); 

                auto tg_xxxxxxzz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 87); 

                auto tg_xxxxxxzz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 88); 

                auto tg_xxxxxxzz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 89); 

                auto tg_xxxxxyyy_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 90); 

                auto tg_xxxxxyyy_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 91); 

                auto tg_xxxxxyyy_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 92); 

                auto tg_xxxxxyyy_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 93); 

                auto tg_xxxxxyyy_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 94); 

                auto tg_xxxxxyyy_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 95); 

                auto tg_xxxxxyyy_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 96); 

                auto tg_xxxxxyyy_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 97); 

                // Batch of Integrals (49,98)

                #pragma omp simd aligned(fxn, fza, tg_xxxxxxyy_xxyz_0, tg_xxxxxxyy_xxzz_0, tg_xxxxxxyy_xyyy_0, \
                                         tg_xxxxxxyy_xyyz_0, tg_xxxxxxyy_xyzz_0, tg_xxxxxxyy_xzzz_0, tg_xxxxxxyy_yyyy_0, \
                                         tg_xxxxxxyy_yyyz_0, tg_xxxxxxyy_yyzz_0, tg_xxxxxxyy_yzzz_0, tg_xxxxxxyy_zzzz_0, \
                                         tg_xxxxxxyz_xxxx_0, tg_xxxxxxyz_xxxy_0, tg_xxxxxxyz_xxxz_0, tg_xxxxxxyz_xxyy_0, \
                                         tg_xxxxxxyz_xxyz_0, tg_xxxxxxyz_xxzz_0, tg_xxxxxxyz_xyyy_0, tg_xxxxxxyz_xyyz_0, \
                                         tg_xxxxxxyz_xyzz_0, tg_xxxxxxyz_xzzz_0, tg_xxxxxxyz_yyyy_0, tg_xxxxxxyz_yyyz_0, \
                                         tg_xxxxxxyz_yyzz_0, tg_xxxxxxyz_yzzz_0, tg_xxxxxxyz_zzzz_0, tg_xxxxxxzz_xxxx_0, \
                                         tg_xxxxxxzz_xxxy_0, tg_xxxxxxzz_xxxz_0, tg_xxxxxxzz_xxyy_0, tg_xxxxxxzz_xxyz_0, \
                                         tg_xxxxxxzz_xxzz_0, tg_xxxxxxzz_xyyy_0, tg_xxxxxxzz_xyyz_0, tg_xxxxxxzz_xyzz_0, \
                                         tg_xxxxxxzz_xzzz_0, tg_xxxxxxzz_yyyy_0, tg_xxxxxxzz_yyyz_0, tg_xxxxxxzz_yyzz_0, \
                                         tg_xxxxxxzz_yzzz_0, tg_xxxxxxzz_zzzz_0, tg_xxxxxyy_xxyz_0, tg_xxxxxyy_xxyz_1, \
                                         tg_xxxxxyy_xxzz_0, tg_xxxxxyy_xxzz_1, tg_xxxxxyy_xyyy_0, tg_xxxxxyy_xyyy_1, \
                                         tg_xxxxxyy_xyyz_0, tg_xxxxxyy_xyyz_1, tg_xxxxxyy_xyz_1, tg_xxxxxyy_xyzz_0, \
                                         tg_xxxxxyy_xyzz_1, tg_xxxxxyy_xzz_1, tg_xxxxxyy_xzzz_0, tg_xxxxxyy_xzzz_1, \
                                         tg_xxxxxyy_yyy_1, tg_xxxxxyy_yyyy_0, tg_xxxxxyy_yyyy_1, tg_xxxxxyy_yyyz_0, \
                                         tg_xxxxxyy_yyyz_1, tg_xxxxxyy_yyz_1, tg_xxxxxyy_yyzz_0, tg_xxxxxyy_yyzz_1, \
                                         tg_xxxxxyy_yzz_1, tg_xxxxxyy_yzzz_0, tg_xxxxxyy_yzzz_1, tg_xxxxxyy_zzz_1, \
                                         tg_xxxxxyy_zzzz_0, tg_xxxxxyy_zzzz_1, tg_xxxxxyyy_xxxx_0, tg_xxxxxyyy_xxxy_0, \
                                         tg_xxxxxyyy_xxxz_0, tg_xxxxxyyy_xxyy_0, tg_xxxxxyyy_xxyz_0, tg_xxxxxyyy_xxzz_0, \
                                         tg_xxxxxyyy_xyyy_0, tg_xxxxxyyy_xyyz_0, tg_xxxxxyz_xxx_1, tg_xxxxxyz_xxxx_0, \
                                         tg_xxxxxyz_xxxx_1, tg_xxxxxyz_xxxy_0, tg_xxxxxyz_xxxy_1, tg_xxxxxyz_xxxz_0, \
                                         tg_xxxxxyz_xxxz_1, tg_xxxxxyz_xxy_1, tg_xxxxxyz_xxyy_0, tg_xxxxxyz_xxyy_1, \
                                         tg_xxxxxyz_xxyz_0, tg_xxxxxyz_xxyz_1, tg_xxxxxyz_xxz_1, tg_xxxxxyz_xxzz_0, \
                                         tg_xxxxxyz_xxzz_1, tg_xxxxxyz_xyy_1, tg_xxxxxyz_xyyy_0, tg_xxxxxyz_xyyy_1, \
                                         tg_xxxxxyz_xyyz_0, tg_xxxxxyz_xyyz_1, tg_xxxxxyz_xyz_1, tg_xxxxxyz_xyzz_0, \
                                         tg_xxxxxyz_xyzz_1, tg_xxxxxyz_xzz_1, tg_xxxxxyz_xzzz_0, tg_xxxxxyz_xzzz_1, \
                                         tg_xxxxxyz_yyy_1, tg_xxxxxyz_yyyy_0, tg_xxxxxyz_yyyy_1, tg_xxxxxyz_yyyz_0, \
                                         tg_xxxxxyz_yyyz_1, tg_xxxxxyz_yyz_1, tg_xxxxxyz_yyzz_0, tg_xxxxxyz_yyzz_1, \
                                         tg_xxxxxyz_yzz_1, tg_xxxxxyz_yzzz_0, tg_xxxxxyz_yzzz_1, tg_xxxxxyz_zzz_1, \
                                         tg_xxxxxyz_zzzz_0, tg_xxxxxyz_zzzz_1, tg_xxxxxzz_xxx_1, tg_xxxxxzz_xxxx_0, \
                                         tg_xxxxxzz_xxxx_1, tg_xxxxxzz_xxxy_0, tg_xxxxxzz_xxxy_1, tg_xxxxxzz_xxxz_0, \
                                         tg_xxxxxzz_xxxz_1, tg_xxxxxzz_xxy_1, tg_xxxxxzz_xxyy_0, tg_xxxxxzz_xxyy_1, \
                                         tg_xxxxxzz_xxyz_0, tg_xxxxxzz_xxyz_1, tg_xxxxxzz_xxz_1, tg_xxxxxzz_xxzz_0, \
                                         tg_xxxxxzz_xxzz_1, tg_xxxxxzz_xyy_1, tg_xxxxxzz_xyyy_0, tg_xxxxxzz_xyyy_1, \
                                         tg_xxxxxzz_xyyz_0, tg_xxxxxzz_xyyz_1, tg_xxxxxzz_xyz_1, tg_xxxxxzz_xyzz_0, \
                                         tg_xxxxxzz_xyzz_1, tg_xxxxxzz_xzz_1, tg_xxxxxzz_xzzz_0, tg_xxxxxzz_xzzz_1, \
                                         tg_xxxxxzz_yyy_1, tg_xxxxxzz_yyyy_0, tg_xxxxxzz_yyyy_1, tg_xxxxxzz_yyyz_0, \
                                         tg_xxxxxzz_yyyz_1, tg_xxxxxzz_yyz_1, tg_xxxxxzz_yyzz_0, tg_xxxxxzz_yyzz_1, \
                                         tg_xxxxxzz_yzz_1, tg_xxxxxzz_yzzz_0, tg_xxxxxzz_yzzz_1, tg_xxxxxzz_zzz_1, \
                                         tg_xxxxxzz_zzzz_0, tg_xxxxxzz_zzzz_1, tg_xxxxyy_xxyz_0, tg_xxxxyy_xxyz_1, \
                                         tg_xxxxyy_xxzz_0, tg_xxxxyy_xxzz_1, tg_xxxxyy_xyyy_0, tg_xxxxyy_xyyy_1, \
                                         tg_xxxxyy_xyyz_0, tg_xxxxyy_xyyz_1, tg_xxxxyy_xyzz_0, tg_xxxxyy_xyzz_1, \
                                         tg_xxxxyy_xzzz_0, tg_xxxxyy_xzzz_1, tg_xxxxyy_yyyy_0, tg_xxxxyy_yyyy_1, \
                                         tg_xxxxyy_yyyz_0, tg_xxxxyy_yyyz_1, tg_xxxxyy_yyzz_0, tg_xxxxyy_yyzz_1, \
                                         tg_xxxxyy_yzzz_0, tg_xxxxyy_yzzz_1, tg_xxxxyy_zzzz_0, tg_xxxxyy_zzzz_1, \
                                         tg_xxxxyyy_xxx_1, tg_xxxxyyy_xxxx_0, tg_xxxxyyy_xxxx_1, tg_xxxxyyy_xxxy_0, \
                                         tg_xxxxyyy_xxxy_1, tg_xxxxyyy_xxxz_0, tg_xxxxyyy_xxxz_1, tg_xxxxyyy_xxy_1, \
                                         tg_xxxxyyy_xxyy_0, tg_xxxxyyy_xxyy_1, tg_xxxxyyy_xxyz_0, tg_xxxxyyy_xxyz_1, \
                                         tg_xxxxyyy_xxz_1, tg_xxxxyyy_xxzz_0, tg_xxxxyyy_xxzz_1, tg_xxxxyyy_xyy_1, \
                                         tg_xxxxyyy_xyyy_0, tg_xxxxyyy_xyyy_1, tg_xxxxyyy_xyyz_0, tg_xxxxyyy_xyyz_1, \
                                         tg_xxxxyyy_xyz_1, tg_xxxxyyy_xzz_1, tg_xxxxyyy_yyy_1, tg_xxxxyyy_yyz_1, \
                                         tg_xxxxyz_xxxx_0, tg_xxxxyz_xxxx_1, tg_xxxxyz_xxxy_0, tg_xxxxyz_xxxy_1, \
                                         tg_xxxxyz_xxxz_0, tg_xxxxyz_xxxz_1, tg_xxxxyz_xxyy_0, tg_xxxxyz_xxyy_1, \
                                         tg_xxxxyz_xxyz_0, tg_xxxxyz_xxyz_1, tg_xxxxyz_xxzz_0, tg_xxxxyz_xxzz_1, \
                                         tg_xxxxyz_xyyy_0, tg_xxxxyz_xyyy_1, tg_xxxxyz_xyyz_0, tg_xxxxyz_xyyz_1, \
                                         tg_xxxxyz_xyzz_0, tg_xxxxyz_xyzz_1, tg_xxxxyz_xzzz_0, tg_xxxxyz_xzzz_1, \
                                         tg_xxxxyz_yyyy_0, tg_xxxxyz_yyyy_1, tg_xxxxyz_yyyz_0, tg_xxxxyz_yyyz_1, \
                                         tg_xxxxyz_yyzz_0, tg_xxxxyz_yyzz_1, tg_xxxxyz_yzzz_0, tg_xxxxyz_yzzz_1, \
                                         tg_xxxxyz_zzzz_0, tg_xxxxyz_zzzz_1, tg_xxxxzz_xxxx_0, tg_xxxxzz_xxxx_1, \
                                         tg_xxxxzz_xxxy_0, tg_xxxxzz_xxxy_1, tg_xxxxzz_xxxz_0, tg_xxxxzz_xxxz_1, \
                                         tg_xxxxzz_xxyy_0, tg_xxxxzz_xxyy_1, tg_xxxxzz_xxyz_0, tg_xxxxzz_xxyz_1, \
                                         tg_xxxxzz_xxzz_0, tg_xxxxzz_xxzz_1, tg_xxxxzz_xyyy_0, tg_xxxxzz_xyyy_1, \
                                         tg_xxxxzz_xyyz_0, tg_xxxxzz_xyyz_1, tg_xxxxzz_xyzz_0, tg_xxxxzz_xyzz_1, \
                                         tg_xxxxzz_xzzz_0, tg_xxxxzz_xzzz_1, tg_xxxxzz_yyyy_0, tg_xxxxzz_yyyy_1, \
                                         tg_xxxxzz_yyyz_0, tg_xxxxzz_yyyz_1, tg_xxxxzz_yyzz_0, tg_xxxxzz_yyzz_1, \
                                         tg_xxxxzz_yzzz_0, tg_xxxxzz_yzzz_1, tg_xxxxzz_zzzz_0, tg_xxxxzz_zzzz_1, \
                                         tg_xxxyyy_xxxx_0, tg_xxxyyy_xxxx_1, tg_xxxyyy_xxxy_0, tg_xxxyyy_xxxy_1, \
                                         tg_xxxyyy_xxxz_0, tg_xxxyyy_xxxz_1, tg_xxxyyy_xxyy_0, tg_xxxyyy_xxyy_1, \
                                         tg_xxxyyy_xxyz_0, tg_xxxyyy_xxyz_1, tg_xxxyyy_xxzz_0, tg_xxxyyy_xxzz_1, \
                                         tg_xxxyyy_xyyy_0, tg_xxxyyy_xyyy_1, tg_xxxyyy_xyyz_0, tg_xxxyyy_xyyz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxxxxxyy_xxyz_0[j] = pb_x * tg_xxxxxyy_xxyz_0[j] + wp_x[j] * tg_xxxxxyy_xxyz_1[j] + 2.5 * fl1_fx * tg_xxxxyy_xxyz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyy_xxyz_1[j] + fl1_fxn * tg_xxxxxyy_xyz_1[j];

                    tg_xxxxxxyy_xxzz_0[j] = pb_x * tg_xxxxxyy_xxzz_0[j] + wp_x[j] * tg_xxxxxyy_xxzz_1[j] + 2.5 * fl1_fx * tg_xxxxyy_xxzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyy_xxzz_1[j] + fl1_fxn * tg_xxxxxyy_xzz_1[j];

                    tg_xxxxxxyy_xyyy_0[j] = pb_x * tg_xxxxxyy_xyyy_0[j] + wp_x[j] * tg_xxxxxyy_xyyy_1[j] + 2.5 * fl1_fx * tg_xxxxyy_xyyy_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyy_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxxxxyy_yyy_1[j];

                    tg_xxxxxxyy_xyyz_0[j] = pb_x * tg_xxxxxyy_xyyz_0[j] + wp_x[j] * tg_xxxxxyy_xyyz_1[j] + 2.5 * fl1_fx * tg_xxxxyy_xyyz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyy_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxxxxyy_yyz_1[j];

                    tg_xxxxxxyy_xyzz_0[j] = pb_x * tg_xxxxxyy_xyzz_0[j] + wp_x[j] * tg_xxxxxyy_xyzz_1[j] + 2.5 * fl1_fx * tg_xxxxyy_xyzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyy_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxxxxyy_yzz_1[j];

                    tg_xxxxxxyy_xzzz_0[j] = pb_x * tg_xxxxxyy_xzzz_0[j] + wp_x[j] * tg_xxxxxyy_xzzz_1[j] + 2.5 * fl1_fx * tg_xxxxyy_xzzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyy_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxxxxyy_zzz_1[j];

                    tg_xxxxxxyy_yyyy_0[j] = pb_x * tg_xxxxxyy_yyyy_0[j] + wp_x[j] * tg_xxxxxyy_yyyy_1[j] + 2.5 * fl1_fx * tg_xxxxyy_yyyy_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyy_yyyy_1[j];

                    tg_xxxxxxyy_yyyz_0[j] = pb_x * tg_xxxxxyy_yyyz_0[j] + wp_x[j] * tg_xxxxxyy_yyyz_1[j] + 2.5 * fl1_fx * tg_xxxxyy_yyyz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyy_yyyz_1[j];

                    tg_xxxxxxyy_yyzz_0[j] = pb_x * tg_xxxxxyy_yyzz_0[j] + wp_x[j] * tg_xxxxxyy_yyzz_1[j] + 2.5 * fl1_fx * tg_xxxxyy_yyzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyy_yyzz_1[j];

                    tg_xxxxxxyy_yzzz_0[j] = pb_x * tg_xxxxxyy_yzzz_0[j] + wp_x[j] * tg_xxxxxyy_yzzz_1[j] + 2.5 * fl1_fx * tg_xxxxyy_yzzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyy_yzzz_1[j];

                    tg_xxxxxxyy_zzzz_0[j] = pb_x * tg_xxxxxyy_zzzz_0[j] + wp_x[j] * tg_xxxxxyy_zzzz_1[j] + 2.5 * fl1_fx * tg_xxxxyy_zzzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyy_zzzz_1[j];

                    tg_xxxxxxyz_xxxx_0[j] = pb_x * tg_xxxxxyz_xxxx_0[j] + wp_x[j] * tg_xxxxxyz_xxxx_1[j] + 2.5 * fl1_fx * tg_xxxxyz_xxxx_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxxxxyz_xxx_1[j];

                    tg_xxxxxxyz_xxxy_0[j] = pb_x * tg_xxxxxyz_xxxy_0[j] + wp_x[j] * tg_xxxxxyz_xxxy_1[j] + 2.5 * fl1_fx * tg_xxxxyz_xxxy_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxxxxyz_xxy_1[j];

                    tg_xxxxxxyz_xxxz_0[j] = pb_x * tg_xxxxxyz_xxxz_0[j] + wp_x[j] * tg_xxxxxyz_xxxz_1[j] + 2.5 * fl1_fx * tg_xxxxyz_xxxz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxxxxyz_xxz_1[j];

                    tg_xxxxxxyz_xxyy_0[j] = pb_x * tg_xxxxxyz_xxyy_0[j] + wp_x[j] * tg_xxxxxyz_xxyy_1[j] + 2.5 * fl1_fx * tg_xxxxyz_xxyy_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyz_xxyy_1[j] + fl1_fxn * tg_xxxxxyz_xyy_1[j];

                    tg_xxxxxxyz_xxyz_0[j] = pb_x * tg_xxxxxyz_xxyz_0[j] + wp_x[j] * tg_xxxxxyz_xxyz_1[j] + 2.5 * fl1_fx * tg_xxxxyz_xxyz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyz_xxyz_1[j] + fl1_fxn * tg_xxxxxyz_xyz_1[j];

                    tg_xxxxxxyz_xxzz_0[j] = pb_x * tg_xxxxxyz_xxzz_0[j] + wp_x[j] * tg_xxxxxyz_xxzz_1[j] + 2.5 * fl1_fx * tg_xxxxyz_xxzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyz_xxzz_1[j] + fl1_fxn * tg_xxxxxyz_xzz_1[j];

                    tg_xxxxxxyz_xyyy_0[j] = pb_x * tg_xxxxxyz_xyyy_0[j] + wp_x[j] * tg_xxxxxyz_xyyy_1[j] + 2.5 * fl1_fx * tg_xxxxyz_xyyy_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxxxxyz_yyy_1[j];

                    tg_xxxxxxyz_xyyz_0[j] = pb_x * tg_xxxxxyz_xyyz_0[j] + wp_x[j] * tg_xxxxxyz_xyyz_1[j] + 2.5 * fl1_fx * tg_xxxxyz_xyyz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxxxxyz_yyz_1[j];

                    tg_xxxxxxyz_xyzz_0[j] = pb_x * tg_xxxxxyz_xyzz_0[j] + wp_x[j] * tg_xxxxxyz_xyzz_1[j] + 2.5 * fl1_fx * tg_xxxxyz_xyzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxxxxyz_yzz_1[j];

                    tg_xxxxxxyz_xzzz_0[j] = pb_x * tg_xxxxxyz_xzzz_0[j] + wp_x[j] * tg_xxxxxyz_xzzz_1[j] + 2.5 * fl1_fx * tg_xxxxyz_xzzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxxxxyz_zzz_1[j];

                    tg_xxxxxxyz_yyyy_0[j] = pb_x * tg_xxxxxyz_yyyy_0[j] + wp_x[j] * tg_xxxxxyz_yyyy_1[j] + 2.5 * fl1_fx * tg_xxxxyz_yyyy_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyz_yyyy_1[j];

                    tg_xxxxxxyz_yyyz_0[j] = pb_x * tg_xxxxxyz_yyyz_0[j] + wp_x[j] * tg_xxxxxyz_yyyz_1[j] + 2.5 * fl1_fx * tg_xxxxyz_yyyz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyz_yyyz_1[j];

                    tg_xxxxxxyz_yyzz_0[j] = pb_x * tg_xxxxxyz_yyzz_0[j] + wp_x[j] * tg_xxxxxyz_yyzz_1[j] + 2.5 * fl1_fx * tg_xxxxyz_yyzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyz_yyzz_1[j];

                    tg_xxxxxxyz_yzzz_0[j] = pb_x * tg_xxxxxyz_yzzz_0[j] + wp_x[j] * tg_xxxxxyz_yzzz_1[j] + 2.5 * fl1_fx * tg_xxxxyz_yzzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyz_yzzz_1[j];

                    tg_xxxxxxyz_zzzz_0[j] = pb_x * tg_xxxxxyz_zzzz_0[j] + wp_x[j] * tg_xxxxxyz_zzzz_1[j] + 2.5 * fl1_fx * tg_xxxxyz_zzzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyz_zzzz_1[j];

                    tg_xxxxxxzz_xxxx_0[j] = pb_x * tg_xxxxxzz_xxxx_0[j] + wp_x[j] * tg_xxxxxzz_xxxx_1[j] + 2.5 * fl1_fx * tg_xxxxzz_xxxx_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxxxxzz_xxx_1[j];

                    tg_xxxxxxzz_xxxy_0[j] = pb_x * tg_xxxxxzz_xxxy_0[j] + wp_x[j] * tg_xxxxxzz_xxxy_1[j] + 2.5 * fl1_fx * tg_xxxxzz_xxxy_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxxxxzz_xxy_1[j];

                    tg_xxxxxxzz_xxxz_0[j] = pb_x * tg_xxxxxzz_xxxz_0[j] + wp_x[j] * tg_xxxxxzz_xxxz_1[j] + 2.5 * fl1_fx * tg_xxxxzz_xxxz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxxxxzz_xxz_1[j];

                    tg_xxxxxxzz_xxyy_0[j] = pb_x * tg_xxxxxzz_xxyy_0[j] + wp_x[j] * tg_xxxxxzz_xxyy_1[j] + 2.5 * fl1_fx * tg_xxxxzz_xxyy_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxzz_xxyy_1[j] + fl1_fxn * tg_xxxxxzz_xyy_1[j];

                    tg_xxxxxxzz_xxyz_0[j] = pb_x * tg_xxxxxzz_xxyz_0[j] + wp_x[j] * tg_xxxxxzz_xxyz_1[j] + 2.5 * fl1_fx * tg_xxxxzz_xxyz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxzz_xxyz_1[j] + fl1_fxn * tg_xxxxxzz_xyz_1[j];

                    tg_xxxxxxzz_xxzz_0[j] = pb_x * tg_xxxxxzz_xxzz_0[j] + wp_x[j] * tg_xxxxxzz_xxzz_1[j] + 2.5 * fl1_fx * tg_xxxxzz_xxzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxzz_xxzz_1[j] + fl1_fxn * tg_xxxxxzz_xzz_1[j];

                    tg_xxxxxxzz_xyyy_0[j] = pb_x * tg_xxxxxzz_xyyy_0[j] + wp_x[j] * tg_xxxxxzz_xyyy_1[j] + 2.5 * fl1_fx * tg_xxxxzz_xyyy_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxxxxzz_yyy_1[j];

                    tg_xxxxxxzz_xyyz_0[j] = pb_x * tg_xxxxxzz_xyyz_0[j] + wp_x[j] * tg_xxxxxzz_xyyz_1[j] + 2.5 * fl1_fx * tg_xxxxzz_xyyz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxxxxzz_yyz_1[j];

                    tg_xxxxxxzz_xyzz_0[j] = pb_x * tg_xxxxxzz_xyzz_0[j] + wp_x[j] * tg_xxxxxzz_xyzz_1[j] + 2.5 * fl1_fx * tg_xxxxzz_xyzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxxxxzz_yzz_1[j];

                    tg_xxxxxxzz_xzzz_0[j] = pb_x * tg_xxxxxzz_xzzz_0[j] + wp_x[j] * tg_xxxxxzz_xzzz_1[j] + 2.5 * fl1_fx * tg_xxxxzz_xzzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxxxxzz_zzz_1[j];

                    tg_xxxxxxzz_yyyy_0[j] = pb_x * tg_xxxxxzz_yyyy_0[j] + wp_x[j] * tg_xxxxxzz_yyyy_1[j] + 2.5 * fl1_fx * tg_xxxxzz_yyyy_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxzz_yyyy_1[j];

                    tg_xxxxxxzz_yyyz_0[j] = pb_x * tg_xxxxxzz_yyyz_0[j] + wp_x[j] * tg_xxxxxzz_yyyz_1[j] + 2.5 * fl1_fx * tg_xxxxzz_yyyz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxzz_yyyz_1[j];

                    tg_xxxxxxzz_yyzz_0[j] = pb_x * tg_xxxxxzz_yyzz_0[j] + wp_x[j] * tg_xxxxxzz_yyzz_1[j] + 2.5 * fl1_fx * tg_xxxxzz_yyzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxzz_yyzz_1[j];

                    tg_xxxxxxzz_yzzz_0[j] = pb_x * tg_xxxxxzz_yzzz_0[j] + wp_x[j] * tg_xxxxxzz_yzzz_1[j] + 2.5 * fl1_fx * tg_xxxxzz_yzzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxzz_yzzz_1[j];

                    tg_xxxxxxzz_zzzz_0[j] = pb_x * tg_xxxxxzz_zzzz_0[j] + wp_x[j] * tg_xxxxxzz_zzzz_1[j] + 2.5 * fl1_fx * tg_xxxxzz_zzzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxzz_zzzz_1[j];

                    tg_xxxxxyyy_xxxx_0[j] = pb_x * tg_xxxxyyy_xxxx_0[j] + wp_x[j] * tg_xxxxyyy_xxxx_1[j] + 2.0 * fl1_fx * tg_xxxyyy_xxxx_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyy_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxxxyyy_xxx_1[j];

                    tg_xxxxxyyy_xxxy_0[j] = pb_x * tg_xxxxyyy_xxxy_0[j] + wp_x[j] * tg_xxxxyyy_xxxy_1[j] + 2.0 * fl1_fx * tg_xxxyyy_xxxy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyy_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxxxyyy_xxy_1[j];

                    tg_xxxxxyyy_xxxz_0[j] = pb_x * tg_xxxxyyy_xxxz_0[j] + wp_x[j] * tg_xxxxyyy_xxxz_1[j] + 2.0 * fl1_fx * tg_xxxyyy_xxxz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyy_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxxxyyy_xxz_1[j];

                    tg_xxxxxyyy_xxyy_0[j] = pb_x * tg_xxxxyyy_xxyy_0[j] + wp_x[j] * tg_xxxxyyy_xxyy_1[j] + 2.0 * fl1_fx * tg_xxxyyy_xxyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyy_xxyy_1[j] + fl1_fxn * tg_xxxxyyy_xyy_1[j];

                    tg_xxxxxyyy_xxyz_0[j] = pb_x * tg_xxxxyyy_xxyz_0[j] + wp_x[j] * tg_xxxxyyy_xxyz_1[j] + 2.0 * fl1_fx * tg_xxxyyy_xxyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyy_xxyz_1[j] + fl1_fxn * tg_xxxxyyy_xyz_1[j];

                    tg_xxxxxyyy_xxzz_0[j] = pb_x * tg_xxxxyyy_xxzz_0[j] + wp_x[j] * tg_xxxxyyy_xxzz_1[j] + 2.0 * fl1_fx * tg_xxxyyy_xxzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyy_xxzz_1[j] + fl1_fxn * tg_xxxxyyy_xzz_1[j];

                    tg_xxxxxyyy_xyyy_0[j] = pb_x * tg_xxxxyyy_xyyy_0[j] + wp_x[j] * tg_xxxxyyy_xyyy_1[j] + 2.0 * fl1_fx * tg_xxxyyy_xyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyy_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxxxyyy_yyy_1[j];

                    tg_xxxxxyyy_xyyz_0[j] = pb_x * tg_xxxxyyy_xyyz_0[j] + wp_x[j] * tg_xxxxyyy_xyyz_1[j] + 2.0 * fl1_fx * tg_xxxyyy_xyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyy_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxxxyyy_yyz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSLSG_98_147(      CMemBlock2D<double>& primBuffer,
                                        const CRecursionMap&       recursionMap,
                                        const CMemBlock2D<double>& osFactors,
                                        const CMemBlock2D<double>& wpDistances,
                                        const CGtoPairsBlock&      braGtoPairsBlock,
                                        const CGtoPairsBlock&      ketGtoPairsBlock,
                                        const int32_t              nKetPrimPairs,
                                        const int32_t              iContrPair)
    {
        // Batch of Integrals (98,147)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {8, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {8, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_7_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_6_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_6_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_7_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_xxxxyyy_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 98); 

                auto tg_xxxxyyy_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 99); 

                auto tg_xxxxyyy_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 100); 

                auto tg_xxxxyyy_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 101); 

                auto tg_xxxxyyy_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 102); 

                auto tg_xxxxyyy_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 103); 

                auto tg_xxxxyyy_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 104); 

                auto tg_xxxxyyz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 105); 

                auto tg_xxxxyyz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 106); 

                auto tg_xxxxyyz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 107); 

                auto tg_xxxxyyz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 108); 

                auto tg_xxxxyyz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 109); 

                auto tg_xxxxyyz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 110); 

                auto tg_xxxxyyz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 111); 

                auto tg_xxxxyyz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 112); 

                auto tg_xxxxyyz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 113); 

                auto tg_xxxxyyz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 114); 

                auto tg_xxxxyyz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 115); 

                auto tg_xxxxyyz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 116); 

                auto tg_xxxxyyz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 117); 

                auto tg_xxxxyyz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 118); 

                auto tg_xxxxyyz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 119); 

                auto tg_xxxxyzz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 120); 

                auto tg_xxxxyzz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 121); 

                auto tg_xxxxyzz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 122); 

                auto tg_xxxxyzz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 123); 

                auto tg_xxxxyzz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 124); 

                auto tg_xxxxyzz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 125); 

                auto tg_xxxxyzz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 126); 

                auto tg_xxxxyzz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 127); 

                auto tg_xxxxyzz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 128); 

                auto tg_xxxxyzz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 129); 

                auto tg_xxxxyzz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 130); 

                auto tg_xxxxyzz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 131); 

                auto tg_xxxxyzz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 132); 

                auto tg_xxxxyzz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 133); 

                auto tg_xxxxyzz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 134); 

                auto tg_xxxxzzz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 135); 

                auto tg_xxxxzzz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 136); 

                auto tg_xxxxzzz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 137); 

                auto tg_xxxxzzz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 138); 

                auto tg_xxxxzzz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 139); 

                auto tg_xxxxzzz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 140); 

                auto tg_xxxxzzz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 141); 

                auto tg_xxxxzzz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 142); 

                auto tg_xxxxzzz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 143); 

                auto tg_xxxxzzz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 144); 

                auto tg_xxxxzzz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 145); 

                auto tg_xxxxzzz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 146); 

                auto tg_xxxxyyy_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 98); 

                auto tg_xxxxyyy_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 99); 

                auto tg_xxxxyyy_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 100); 

                auto tg_xxxxyyy_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 101); 

                auto tg_xxxxyyy_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 102); 

                auto tg_xxxxyyy_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 103); 

                auto tg_xxxxyyy_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 104); 

                auto tg_xxxxyyz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 105); 

                auto tg_xxxxyyz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 106); 

                auto tg_xxxxyyz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 107); 

                auto tg_xxxxyyz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 108); 

                auto tg_xxxxyyz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 109); 

                auto tg_xxxxyyz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 110); 

                auto tg_xxxxyyz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 111); 

                auto tg_xxxxyyz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 112); 

                auto tg_xxxxyyz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 113); 

                auto tg_xxxxyyz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 114); 

                auto tg_xxxxyyz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 115); 

                auto tg_xxxxyyz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 116); 

                auto tg_xxxxyyz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 117); 

                auto tg_xxxxyyz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 118); 

                auto tg_xxxxyyz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 119); 

                auto tg_xxxxyzz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 120); 

                auto tg_xxxxyzz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 121); 

                auto tg_xxxxyzz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 122); 

                auto tg_xxxxyzz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 123); 

                auto tg_xxxxyzz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 124); 

                auto tg_xxxxyzz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 125); 

                auto tg_xxxxyzz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 126); 

                auto tg_xxxxyzz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 127); 

                auto tg_xxxxyzz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 128); 

                auto tg_xxxxyzz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 129); 

                auto tg_xxxxyzz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 130); 

                auto tg_xxxxyzz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 131); 

                auto tg_xxxxyzz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 132); 

                auto tg_xxxxyzz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 133); 

                auto tg_xxxxyzz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 134); 

                auto tg_xxxxzzz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 135); 

                auto tg_xxxxzzz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 136); 

                auto tg_xxxxzzz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 137); 

                auto tg_xxxxzzz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 138); 

                auto tg_xxxxzzz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 139); 

                auto tg_xxxxzzz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 140); 

                auto tg_xxxxzzz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 141); 

                auto tg_xxxxzzz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 142); 

                auto tg_xxxxzzz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 143); 

                auto tg_xxxxzzz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 144); 

                auto tg_xxxxzzz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 145); 

                auto tg_xxxxzzz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 146); 

                auto tg_xxxyyy_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 98); 

                auto tg_xxxyyy_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 99); 

                auto tg_xxxyyy_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 100); 

                auto tg_xxxyyy_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 101); 

                auto tg_xxxyyy_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 102); 

                auto tg_xxxyyy_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 103); 

                auto tg_xxxyyy_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 104); 

                auto tg_xxxyyz_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 105); 

                auto tg_xxxyyz_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 106); 

                auto tg_xxxyyz_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 107); 

                auto tg_xxxyyz_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 108); 

                auto tg_xxxyyz_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 109); 

                auto tg_xxxyyz_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 110); 

                auto tg_xxxyyz_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 111); 

                auto tg_xxxyyz_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 112); 

                auto tg_xxxyyz_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 113); 

                auto tg_xxxyyz_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 114); 

                auto tg_xxxyyz_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 115); 

                auto tg_xxxyyz_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 116); 

                auto tg_xxxyyz_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 117); 

                auto tg_xxxyyz_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 118); 

                auto tg_xxxyyz_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 119); 

                auto tg_xxxyzz_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 120); 

                auto tg_xxxyzz_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 121); 

                auto tg_xxxyzz_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 122); 

                auto tg_xxxyzz_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 123); 

                auto tg_xxxyzz_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 124); 

                auto tg_xxxyzz_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 125); 

                auto tg_xxxyzz_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 126); 

                auto tg_xxxyzz_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 127); 

                auto tg_xxxyzz_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 128); 

                auto tg_xxxyzz_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 129); 

                auto tg_xxxyzz_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 130); 

                auto tg_xxxyzz_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 131); 

                auto tg_xxxyzz_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 132); 

                auto tg_xxxyzz_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 133); 

                auto tg_xxxyzz_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 134); 

                auto tg_xxxzzz_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 135); 

                auto tg_xxxzzz_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 136); 

                auto tg_xxxzzz_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 137); 

                auto tg_xxxzzz_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 138); 

                auto tg_xxxzzz_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 139); 

                auto tg_xxxzzz_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 140); 

                auto tg_xxxzzz_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 141); 

                auto tg_xxxzzz_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 142); 

                auto tg_xxxzzz_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 143); 

                auto tg_xxxzzz_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 144); 

                auto tg_xxxzzz_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 145); 

                auto tg_xxxzzz_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 146); 

                auto tg_xxxyyy_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 98); 

                auto tg_xxxyyy_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 99); 

                auto tg_xxxyyy_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 100); 

                auto tg_xxxyyy_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 101); 

                auto tg_xxxyyy_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 102); 

                auto tg_xxxyyy_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 103); 

                auto tg_xxxyyy_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 104); 

                auto tg_xxxyyz_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 105); 

                auto tg_xxxyyz_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 106); 

                auto tg_xxxyyz_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 107); 

                auto tg_xxxyyz_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 108); 

                auto tg_xxxyyz_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 109); 

                auto tg_xxxyyz_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 110); 

                auto tg_xxxyyz_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 111); 

                auto tg_xxxyyz_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 112); 

                auto tg_xxxyyz_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 113); 

                auto tg_xxxyyz_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 114); 

                auto tg_xxxyyz_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 115); 

                auto tg_xxxyyz_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 116); 

                auto tg_xxxyyz_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 117); 

                auto tg_xxxyyz_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 118); 

                auto tg_xxxyyz_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 119); 

                auto tg_xxxyzz_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 120); 

                auto tg_xxxyzz_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 121); 

                auto tg_xxxyzz_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 122); 

                auto tg_xxxyzz_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 123); 

                auto tg_xxxyzz_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 124); 

                auto tg_xxxyzz_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 125); 

                auto tg_xxxyzz_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 126); 

                auto tg_xxxyzz_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 127); 

                auto tg_xxxyzz_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 128); 

                auto tg_xxxyzz_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 129); 

                auto tg_xxxyzz_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 130); 

                auto tg_xxxyzz_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 131); 

                auto tg_xxxyzz_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 132); 

                auto tg_xxxyzz_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 133); 

                auto tg_xxxyzz_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 134); 

                auto tg_xxxzzz_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 135); 

                auto tg_xxxzzz_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 136); 

                auto tg_xxxzzz_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 137); 

                auto tg_xxxzzz_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 138); 

                auto tg_xxxzzz_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 139); 

                auto tg_xxxzzz_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 140); 

                auto tg_xxxzzz_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 141); 

                auto tg_xxxzzz_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 142); 

                auto tg_xxxzzz_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 143); 

                auto tg_xxxzzz_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 144); 

                auto tg_xxxzzz_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 145); 

                auto tg_xxxzzz_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 146); 

                auto tg_xxxxyyy_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 68); 

                auto tg_xxxxyyy_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 69); 

                auto tg_xxxxyyz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 70); 

                auto tg_xxxxyyz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 71); 

                auto tg_xxxxyyz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 72); 

                auto tg_xxxxyyz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 73); 

                auto tg_xxxxyyz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 74); 

                auto tg_xxxxyyz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 75); 

                auto tg_xxxxyyz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 76); 

                auto tg_xxxxyyz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 77); 

                auto tg_xxxxyyz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 78); 

                auto tg_xxxxyyz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 79); 

                auto tg_xxxxyzz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 80); 

                auto tg_xxxxyzz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 81); 

                auto tg_xxxxyzz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 82); 

                auto tg_xxxxyzz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 83); 

                auto tg_xxxxyzz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 84); 

                auto tg_xxxxyzz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 85); 

                auto tg_xxxxyzz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 86); 

                auto tg_xxxxyzz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 87); 

                auto tg_xxxxyzz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 88); 

                auto tg_xxxxyzz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 89); 

                auto tg_xxxxzzz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 90); 

                auto tg_xxxxzzz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 91); 

                auto tg_xxxxzzz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 92); 

                auto tg_xxxxzzz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 93); 

                auto tg_xxxxzzz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 94); 

                auto tg_xxxxzzz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 95); 

                auto tg_xxxxzzz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 96); 

                auto tg_xxxxzzz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 97); 

                auto tg_xxxxzzz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 98); 

                auto tg_xxxxzzz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 99); 

                // set up pointers to integrals

                auto tg_xxxxxyyy_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 98); 

                auto tg_xxxxxyyy_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 99); 

                auto tg_xxxxxyyy_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 100); 

                auto tg_xxxxxyyy_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 101); 

                auto tg_xxxxxyyy_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 102); 

                auto tg_xxxxxyyy_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 103); 

                auto tg_xxxxxyyy_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 104); 

                auto tg_xxxxxyyz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 105); 

                auto tg_xxxxxyyz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 106); 

                auto tg_xxxxxyyz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 107); 

                auto tg_xxxxxyyz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 108); 

                auto tg_xxxxxyyz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 109); 

                auto tg_xxxxxyyz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 110); 

                auto tg_xxxxxyyz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 111); 

                auto tg_xxxxxyyz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 112); 

                auto tg_xxxxxyyz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 113); 

                auto tg_xxxxxyyz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 114); 

                auto tg_xxxxxyyz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 115); 

                auto tg_xxxxxyyz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 116); 

                auto tg_xxxxxyyz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 117); 

                auto tg_xxxxxyyz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 118); 

                auto tg_xxxxxyyz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 119); 

                auto tg_xxxxxyzz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 120); 

                auto tg_xxxxxyzz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 121); 

                auto tg_xxxxxyzz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 122); 

                auto tg_xxxxxyzz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 123); 

                auto tg_xxxxxyzz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 124); 

                auto tg_xxxxxyzz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 125); 

                auto tg_xxxxxyzz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 126); 

                auto tg_xxxxxyzz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 127); 

                auto tg_xxxxxyzz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 128); 

                auto tg_xxxxxyzz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 129); 

                auto tg_xxxxxyzz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 130); 

                auto tg_xxxxxyzz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 131); 

                auto tg_xxxxxyzz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 132); 

                auto tg_xxxxxyzz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 133); 

                auto tg_xxxxxyzz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 134); 

                auto tg_xxxxxzzz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 135); 

                auto tg_xxxxxzzz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 136); 

                auto tg_xxxxxzzz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 137); 

                auto tg_xxxxxzzz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 138); 

                auto tg_xxxxxzzz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 139); 

                auto tg_xxxxxzzz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 140); 

                auto tg_xxxxxzzz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 141); 

                auto tg_xxxxxzzz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 142); 

                auto tg_xxxxxzzz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 143); 

                auto tg_xxxxxzzz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 144); 

                auto tg_xxxxxzzz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 145); 

                auto tg_xxxxxzzz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 146); 

                // Batch of Integrals (98,147)

                #pragma omp simd aligned(fxn, fza, tg_xxxxxyyy_xyzz_0, tg_xxxxxyyy_xzzz_0, tg_xxxxxyyy_yyyy_0, \
                                         tg_xxxxxyyy_yyyz_0, tg_xxxxxyyy_yyzz_0, tg_xxxxxyyy_yzzz_0, tg_xxxxxyyy_zzzz_0, \
                                         tg_xxxxxyyz_xxxx_0, tg_xxxxxyyz_xxxy_0, tg_xxxxxyyz_xxxz_0, tg_xxxxxyyz_xxyy_0, \
                                         tg_xxxxxyyz_xxyz_0, tg_xxxxxyyz_xxzz_0, tg_xxxxxyyz_xyyy_0, tg_xxxxxyyz_xyyz_0, \
                                         tg_xxxxxyyz_xyzz_0, tg_xxxxxyyz_xzzz_0, tg_xxxxxyyz_yyyy_0, tg_xxxxxyyz_yyyz_0, \
                                         tg_xxxxxyyz_yyzz_0, tg_xxxxxyyz_yzzz_0, tg_xxxxxyyz_zzzz_0, tg_xxxxxyzz_xxxx_0, \
                                         tg_xxxxxyzz_xxxy_0, tg_xxxxxyzz_xxxz_0, tg_xxxxxyzz_xxyy_0, tg_xxxxxyzz_xxyz_0, \
                                         tg_xxxxxyzz_xxzz_0, tg_xxxxxyzz_xyyy_0, tg_xxxxxyzz_xyyz_0, tg_xxxxxyzz_xyzz_0, \
                                         tg_xxxxxyzz_xzzz_0, tg_xxxxxyzz_yyyy_0, tg_xxxxxyzz_yyyz_0, tg_xxxxxyzz_yyzz_0, \
                                         tg_xxxxxyzz_yzzz_0, tg_xxxxxyzz_zzzz_0, tg_xxxxxzzz_xxxx_0, tg_xxxxxzzz_xxxy_0, \
                                         tg_xxxxxzzz_xxxz_0, tg_xxxxxzzz_xxyy_0, tg_xxxxxzzz_xxyz_0, tg_xxxxxzzz_xxzz_0, \
                                         tg_xxxxxzzz_xyyy_0, tg_xxxxxzzz_xyyz_0, tg_xxxxxzzz_xyzz_0, tg_xxxxxzzz_xzzz_0, \
                                         tg_xxxxxzzz_yyyy_0, tg_xxxxxzzz_yyyz_0, tg_xxxxyyy_xyzz_0, tg_xxxxyyy_xyzz_1, \
                                         tg_xxxxyyy_xzzz_0, tg_xxxxyyy_xzzz_1, tg_xxxxyyy_yyyy_0, tg_xxxxyyy_yyyy_1, \
                                         tg_xxxxyyy_yyyz_0, tg_xxxxyyy_yyyz_1, tg_xxxxyyy_yyzz_0, tg_xxxxyyy_yyzz_1, \
                                         tg_xxxxyyy_yzz_1, tg_xxxxyyy_yzzz_0, tg_xxxxyyy_yzzz_1, tg_xxxxyyy_zzz_1, \
                                         tg_xxxxyyy_zzzz_0, tg_xxxxyyy_zzzz_1, tg_xxxxyyz_xxx_1, tg_xxxxyyz_xxxx_0, \
                                         tg_xxxxyyz_xxxx_1, tg_xxxxyyz_xxxy_0, tg_xxxxyyz_xxxy_1, tg_xxxxyyz_xxxz_0, \
                                         tg_xxxxyyz_xxxz_1, tg_xxxxyyz_xxy_1, tg_xxxxyyz_xxyy_0, tg_xxxxyyz_xxyy_1, \
                                         tg_xxxxyyz_xxyz_0, tg_xxxxyyz_xxyz_1, tg_xxxxyyz_xxz_1, tg_xxxxyyz_xxzz_0, \
                                         tg_xxxxyyz_xxzz_1, tg_xxxxyyz_xyy_1, tg_xxxxyyz_xyyy_0, tg_xxxxyyz_xyyy_1, \
                                         tg_xxxxyyz_xyyz_0, tg_xxxxyyz_xyyz_1, tg_xxxxyyz_xyz_1, tg_xxxxyyz_xyzz_0, \
                                         tg_xxxxyyz_xyzz_1, tg_xxxxyyz_xzz_1, tg_xxxxyyz_xzzz_0, tg_xxxxyyz_xzzz_1, \
                                         tg_xxxxyyz_yyy_1, tg_xxxxyyz_yyyy_0, tg_xxxxyyz_yyyy_1, tg_xxxxyyz_yyyz_0, \
                                         tg_xxxxyyz_yyyz_1, tg_xxxxyyz_yyz_1, tg_xxxxyyz_yyzz_0, tg_xxxxyyz_yyzz_1, \
                                         tg_xxxxyyz_yzz_1, tg_xxxxyyz_yzzz_0, tg_xxxxyyz_yzzz_1, tg_xxxxyyz_zzz_1, \
                                         tg_xxxxyyz_zzzz_0, tg_xxxxyyz_zzzz_1, tg_xxxxyzz_xxx_1, tg_xxxxyzz_xxxx_0, \
                                         tg_xxxxyzz_xxxx_1, tg_xxxxyzz_xxxy_0, tg_xxxxyzz_xxxy_1, tg_xxxxyzz_xxxz_0, \
                                         tg_xxxxyzz_xxxz_1, tg_xxxxyzz_xxy_1, tg_xxxxyzz_xxyy_0, tg_xxxxyzz_xxyy_1, \
                                         tg_xxxxyzz_xxyz_0, tg_xxxxyzz_xxyz_1, tg_xxxxyzz_xxz_1, tg_xxxxyzz_xxzz_0, \
                                         tg_xxxxyzz_xxzz_1, tg_xxxxyzz_xyy_1, tg_xxxxyzz_xyyy_0, tg_xxxxyzz_xyyy_1, \
                                         tg_xxxxyzz_xyyz_0, tg_xxxxyzz_xyyz_1, tg_xxxxyzz_xyz_1, tg_xxxxyzz_xyzz_0, \
                                         tg_xxxxyzz_xyzz_1, tg_xxxxyzz_xzz_1, tg_xxxxyzz_xzzz_0, tg_xxxxyzz_xzzz_1, \
                                         tg_xxxxyzz_yyy_1, tg_xxxxyzz_yyyy_0, tg_xxxxyzz_yyyy_1, tg_xxxxyzz_yyyz_0, \
                                         tg_xxxxyzz_yyyz_1, tg_xxxxyzz_yyz_1, tg_xxxxyzz_yyzz_0, tg_xxxxyzz_yyzz_1, \
                                         tg_xxxxyzz_yzz_1, tg_xxxxyzz_yzzz_0, tg_xxxxyzz_yzzz_1, tg_xxxxyzz_zzz_1, \
                                         tg_xxxxyzz_zzzz_0, tg_xxxxyzz_zzzz_1, tg_xxxxzzz_xxx_1, tg_xxxxzzz_xxxx_0, \
                                         tg_xxxxzzz_xxxx_1, tg_xxxxzzz_xxxy_0, tg_xxxxzzz_xxxy_1, tg_xxxxzzz_xxxz_0, \
                                         tg_xxxxzzz_xxxz_1, tg_xxxxzzz_xxy_1, tg_xxxxzzz_xxyy_0, tg_xxxxzzz_xxyy_1, \
                                         tg_xxxxzzz_xxyz_0, tg_xxxxzzz_xxyz_1, tg_xxxxzzz_xxz_1, tg_xxxxzzz_xxzz_0, \
                                         tg_xxxxzzz_xxzz_1, tg_xxxxzzz_xyy_1, tg_xxxxzzz_xyyy_0, tg_xxxxzzz_xyyy_1, \
                                         tg_xxxxzzz_xyyz_0, tg_xxxxzzz_xyyz_1, tg_xxxxzzz_xyz_1, tg_xxxxzzz_xyzz_0, \
                                         tg_xxxxzzz_xyzz_1, tg_xxxxzzz_xzz_1, tg_xxxxzzz_xzzz_0, tg_xxxxzzz_xzzz_1, \
                                         tg_xxxxzzz_yyy_1, tg_xxxxzzz_yyyy_0, tg_xxxxzzz_yyyy_1, tg_xxxxzzz_yyyz_0, \
                                         tg_xxxxzzz_yyyz_1, tg_xxxxzzz_yyz_1, tg_xxxxzzz_yzz_1, tg_xxxxzzz_zzz_1, \
                                         tg_xxxyyy_xyzz_0, tg_xxxyyy_xyzz_1, tg_xxxyyy_xzzz_0, tg_xxxyyy_xzzz_1, \
                                         tg_xxxyyy_yyyy_0, tg_xxxyyy_yyyy_1, tg_xxxyyy_yyyz_0, tg_xxxyyy_yyyz_1, \
                                         tg_xxxyyy_yyzz_0, tg_xxxyyy_yyzz_1, tg_xxxyyy_yzzz_0, tg_xxxyyy_yzzz_1, \
                                         tg_xxxyyy_zzzz_0, tg_xxxyyy_zzzz_1, tg_xxxyyz_xxxx_0, tg_xxxyyz_xxxx_1, \
                                         tg_xxxyyz_xxxy_0, tg_xxxyyz_xxxy_1, tg_xxxyyz_xxxz_0, tg_xxxyyz_xxxz_1, \
                                         tg_xxxyyz_xxyy_0, tg_xxxyyz_xxyy_1, tg_xxxyyz_xxyz_0, tg_xxxyyz_xxyz_1, \
                                         tg_xxxyyz_xxzz_0, tg_xxxyyz_xxzz_1, tg_xxxyyz_xyyy_0, tg_xxxyyz_xyyy_1, \
                                         tg_xxxyyz_xyyz_0, tg_xxxyyz_xyyz_1, tg_xxxyyz_xyzz_0, tg_xxxyyz_xyzz_1, \
                                         tg_xxxyyz_xzzz_0, tg_xxxyyz_xzzz_1, tg_xxxyyz_yyyy_0, tg_xxxyyz_yyyy_1, \
                                         tg_xxxyyz_yyyz_0, tg_xxxyyz_yyyz_1, tg_xxxyyz_yyzz_0, tg_xxxyyz_yyzz_1, \
                                         tg_xxxyyz_yzzz_0, tg_xxxyyz_yzzz_1, tg_xxxyyz_zzzz_0, tg_xxxyyz_zzzz_1, \
                                         tg_xxxyzz_xxxx_0, tg_xxxyzz_xxxx_1, tg_xxxyzz_xxxy_0, tg_xxxyzz_xxxy_1, \
                                         tg_xxxyzz_xxxz_0, tg_xxxyzz_xxxz_1, tg_xxxyzz_xxyy_0, tg_xxxyzz_xxyy_1, \
                                         tg_xxxyzz_xxyz_0, tg_xxxyzz_xxyz_1, tg_xxxyzz_xxzz_0, tg_xxxyzz_xxzz_1, \
                                         tg_xxxyzz_xyyy_0, tg_xxxyzz_xyyy_1, tg_xxxyzz_xyyz_0, tg_xxxyzz_xyyz_1, \
                                         tg_xxxyzz_xyzz_0, tg_xxxyzz_xyzz_1, tg_xxxyzz_xzzz_0, tg_xxxyzz_xzzz_1, \
                                         tg_xxxyzz_yyyy_0, tg_xxxyzz_yyyy_1, tg_xxxyzz_yyyz_0, tg_xxxyzz_yyyz_1, \
                                         tg_xxxyzz_yyzz_0, tg_xxxyzz_yyzz_1, tg_xxxyzz_yzzz_0, tg_xxxyzz_yzzz_1, \
                                         tg_xxxyzz_zzzz_0, tg_xxxyzz_zzzz_1, tg_xxxzzz_xxxx_0, tg_xxxzzz_xxxx_1, \
                                         tg_xxxzzz_xxxy_0, tg_xxxzzz_xxxy_1, tg_xxxzzz_xxxz_0, tg_xxxzzz_xxxz_1, \
                                         tg_xxxzzz_xxyy_0, tg_xxxzzz_xxyy_1, tg_xxxzzz_xxyz_0, tg_xxxzzz_xxyz_1, \
                                         tg_xxxzzz_xxzz_0, tg_xxxzzz_xxzz_1, tg_xxxzzz_xyyy_0, tg_xxxzzz_xyyy_1, \
                                         tg_xxxzzz_xyyz_0, tg_xxxzzz_xyyz_1, tg_xxxzzz_xyzz_0, tg_xxxzzz_xyzz_1, \
                                         tg_xxxzzz_xzzz_0, tg_xxxzzz_xzzz_1, tg_xxxzzz_yyyy_0, tg_xxxzzz_yyyy_1, \
                                         tg_xxxzzz_yyyz_0, tg_xxxzzz_yyyz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxxxxyyy_xyzz_0[j] = pb_x * tg_xxxxyyy_xyzz_0[j] + wp_x[j] * tg_xxxxyyy_xyzz_1[j] + 2.0 * fl1_fx * tg_xxxyyy_xyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyy_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxxxyyy_yzz_1[j];

                    tg_xxxxxyyy_xzzz_0[j] = pb_x * tg_xxxxyyy_xzzz_0[j] + wp_x[j] * tg_xxxxyyy_xzzz_1[j] + 2.0 * fl1_fx * tg_xxxyyy_xzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyy_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxxxyyy_zzz_1[j];

                    tg_xxxxxyyy_yyyy_0[j] = pb_x * tg_xxxxyyy_yyyy_0[j] + wp_x[j] * tg_xxxxyyy_yyyy_1[j] + 2.0 * fl1_fx * tg_xxxyyy_yyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyy_yyyy_1[j];

                    tg_xxxxxyyy_yyyz_0[j] = pb_x * tg_xxxxyyy_yyyz_0[j] + wp_x[j] * tg_xxxxyyy_yyyz_1[j] + 2.0 * fl1_fx * tg_xxxyyy_yyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyy_yyyz_1[j];

                    tg_xxxxxyyy_yyzz_0[j] = pb_x * tg_xxxxyyy_yyzz_0[j] + wp_x[j] * tg_xxxxyyy_yyzz_1[j] + 2.0 * fl1_fx * tg_xxxyyy_yyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyy_yyzz_1[j];

                    tg_xxxxxyyy_yzzz_0[j] = pb_x * tg_xxxxyyy_yzzz_0[j] + wp_x[j] * tg_xxxxyyy_yzzz_1[j] + 2.0 * fl1_fx * tg_xxxyyy_yzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyy_yzzz_1[j];

                    tg_xxxxxyyy_zzzz_0[j] = pb_x * tg_xxxxyyy_zzzz_0[j] + wp_x[j] * tg_xxxxyyy_zzzz_1[j] + 2.0 * fl1_fx * tg_xxxyyy_zzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyy_zzzz_1[j];

                    tg_xxxxxyyz_xxxx_0[j] = pb_x * tg_xxxxyyz_xxxx_0[j] + wp_x[j] * tg_xxxxyyz_xxxx_1[j] + 2.0 * fl1_fx * tg_xxxyyz_xxxx_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxxxyyz_xxx_1[j];

                    tg_xxxxxyyz_xxxy_0[j] = pb_x * tg_xxxxyyz_xxxy_0[j] + wp_x[j] * tg_xxxxyyz_xxxy_1[j] + 2.0 * fl1_fx * tg_xxxyyz_xxxy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxxxyyz_xxy_1[j];

                    tg_xxxxxyyz_xxxz_0[j] = pb_x * tg_xxxxyyz_xxxz_0[j] + wp_x[j] * tg_xxxxyyz_xxxz_1[j] + 2.0 * fl1_fx * tg_xxxyyz_xxxz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxxxyyz_xxz_1[j];

                    tg_xxxxxyyz_xxyy_0[j] = pb_x * tg_xxxxyyz_xxyy_0[j] + wp_x[j] * tg_xxxxyyz_xxyy_1[j] + 2.0 * fl1_fx * tg_xxxyyz_xxyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyz_xxyy_1[j] + fl1_fxn * tg_xxxxyyz_xyy_1[j];

                    tg_xxxxxyyz_xxyz_0[j] = pb_x * tg_xxxxyyz_xxyz_0[j] + wp_x[j] * tg_xxxxyyz_xxyz_1[j] + 2.0 * fl1_fx * tg_xxxyyz_xxyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyz_xxyz_1[j] + fl1_fxn * tg_xxxxyyz_xyz_1[j];

                    tg_xxxxxyyz_xxzz_0[j] = pb_x * tg_xxxxyyz_xxzz_0[j] + wp_x[j] * tg_xxxxyyz_xxzz_1[j] + 2.0 * fl1_fx * tg_xxxyyz_xxzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyz_xxzz_1[j] + fl1_fxn * tg_xxxxyyz_xzz_1[j];

                    tg_xxxxxyyz_xyyy_0[j] = pb_x * tg_xxxxyyz_xyyy_0[j] + wp_x[j] * tg_xxxxyyz_xyyy_1[j] + 2.0 * fl1_fx * tg_xxxyyz_xyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxxxyyz_yyy_1[j];

                    tg_xxxxxyyz_xyyz_0[j] = pb_x * tg_xxxxyyz_xyyz_0[j] + wp_x[j] * tg_xxxxyyz_xyyz_1[j] + 2.0 * fl1_fx * tg_xxxyyz_xyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxxxyyz_yyz_1[j];

                    tg_xxxxxyyz_xyzz_0[j] = pb_x * tg_xxxxyyz_xyzz_0[j] + wp_x[j] * tg_xxxxyyz_xyzz_1[j] + 2.0 * fl1_fx * tg_xxxyyz_xyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxxxyyz_yzz_1[j];

                    tg_xxxxxyyz_xzzz_0[j] = pb_x * tg_xxxxyyz_xzzz_0[j] + wp_x[j] * tg_xxxxyyz_xzzz_1[j] + 2.0 * fl1_fx * tg_xxxyyz_xzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxxxyyz_zzz_1[j];

                    tg_xxxxxyyz_yyyy_0[j] = pb_x * tg_xxxxyyz_yyyy_0[j] + wp_x[j] * tg_xxxxyyz_yyyy_1[j] + 2.0 * fl1_fx * tg_xxxyyz_yyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyz_yyyy_1[j];

                    tg_xxxxxyyz_yyyz_0[j] = pb_x * tg_xxxxyyz_yyyz_0[j] + wp_x[j] * tg_xxxxyyz_yyyz_1[j] + 2.0 * fl1_fx * tg_xxxyyz_yyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyz_yyyz_1[j];

                    tg_xxxxxyyz_yyzz_0[j] = pb_x * tg_xxxxyyz_yyzz_0[j] + wp_x[j] * tg_xxxxyyz_yyzz_1[j] + 2.0 * fl1_fx * tg_xxxyyz_yyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyz_yyzz_1[j];

                    tg_xxxxxyyz_yzzz_0[j] = pb_x * tg_xxxxyyz_yzzz_0[j] + wp_x[j] * tg_xxxxyyz_yzzz_1[j] + 2.0 * fl1_fx * tg_xxxyyz_yzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyz_yzzz_1[j];

                    tg_xxxxxyyz_zzzz_0[j] = pb_x * tg_xxxxyyz_zzzz_0[j] + wp_x[j] * tg_xxxxyyz_zzzz_1[j] + 2.0 * fl1_fx * tg_xxxyyz_zzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyz_zzzz_1[j];

                    tg_xxxxxyzz_xxxx_0[j] = pb_x * tg_xxxxyzz_xxxx_0[j] + wp_x[j] * tg_xxxxyzz_xxxx_1[j] + 2.0 * fl1_fx * tg_xxxyzz_xxxx_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxxxyzz_xxx_1[j];

                    tg_xxxxxyzz_xxxy_0[j] = pb_x * tg_xxxxyzz_xxxy_0[j] + wp_x[j] * tg_xxxxyzz_xxxy_1[j] + 2.0 * fl1_fx * tg_xxxyzz_xxxy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxxxyzz_xxy_1[j];

                    tg_xxxxxyzz_xxxz_0[j] = pb_x * tg_xxxxyzz_xxxz_0[j] + wp_x[j] * tg_xxxxyzz_xxxz_1[j] + 2.0 * fl1_fx * tg_xxxyzz_xxxz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxxxyzz_xxz_1[j];

                    tg_xxxxxyzz_xxyy_0[j] = pb_x * tg_xxxxyzz_xxyy_0[j] + wp_x[j] * tg_xxxxyzz_xxyy_1[j] + 2.0 * fl1_fx * tg_xxxyzz_xxyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyzz_xxyy_1[j] + fl1_fxn * tg_xxxxyzz_xyy_1[j];

                    tg_xxxxxyzz_xxyz_0[j] = pb_x * tg_xxxxyzz_xxyz_0[j] + wp_x[j] * tg_xxxxyzz_xxyz_1[j] + 2.0 * fl1_fx * tg_xxxyzz_xxyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyzz_xxyz_1[j] + fl1_fxn * tg_xxxxyzz_xyz_1[j];

                    tg_xxxxxyzz_xxzz_0[j] = pb_x * tg_xxxxyzz_xxzz_0[j] + wp_x[j] * tg_xxxxyzz_xxzz_1[j] + 2.0 * fl1_fx * tg_xxxyzz_xxzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyzz_xxzz_1[j] + fl1_fxn * tg_xxxxyzz_xzz_1[j];

                    tg_xxxxxyzz_xyyy_0[j] = pb_x * tg_xxxxyzz_xyyy_0[j] + wp_x[j] * tg_xxxxyzz_xyyy_1[j] + 2.0 * fl1_fx * tg_xxxyzz_xyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxxxyzz_yyy_1[j];

                    tg_xxxxxyzz_xyyz_0[j] = pb_x * tg_xxxxyzz_xyyz_0[j] + wp_x[j] * tg_xxxxyzz_xyyz_1[j] + 2.0 * fl1_fx * tg_xxxyzz_xyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxxxyzz_yyz_1[j];

                    tg_xxxxxyzz_xyzz_0[j] = pb_x * tg_xxxxyzz_xyzz_0[j] + wp_x[j] * tg_xxxxyzz_xyzz_1[j] + 2.0 * fl1_fx * tg_xxxyzz_xyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxxxyzz_yzz_1[j];

                    tg_xxxxxyzz_xzzz_0[j] = pb_x * tg_xxxxyzz_xzzz_0[j] + wp_x[j] * tg_xxxxyzz_xzzz_1[j] + 2.0 * fl1_fx * tg_xxxyzz_xzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxxxyzz_zzz_1[j];

                    tg_xxxxxyzz_yyyy_0[j] = pb_x * tg_xxxxyzz_yyyy_0[j] + wp_x[j] * tg_xxxxyzz_yyyy_1[j] + 2.0 * fl1_fx * tg_xxxyzz_yyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyzz_yyyy_1[j];

                    tg_xxxxxyzz_yyyz_0[j] = pb_x * tg_xxxxyzz_yyyz_0[j] + wp_x[j] * tg_xxxxyzz_yyyz_1[j] + 2.0 * fl1_fx * tg_xxxyzz_yyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyzz_yyyz_1[j];

                    tg_xxxxxyzz_yyzz_0[j] = pb_x * tg_xxxxyzz_yyzz_0[j] + wp_x[j] * tg_xxxxyzz_yyzz_1[j] + 2.0 * fl1_fx * tg_xxxyzz_yyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyzz_yyzz_1[j];

                    tg_xxxxxyzz_yzzz_0[j] = pb_x * tg_xxxxyzz_yzzz_0[j] + wp_x[j] * tg_xxxxyzz_yzzz_1[j] + 2.0 * fl1_fx * tg_xxxyzz_yzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyzz_yzzz_1[j];

                    tg_xxxxxyzz_zzzz_0[j] = pb_x * tg_xxxxyzz_zzzz_0[j] + wp_x[j] * tg_xxxxyzz_zzzz_1[j] + 2.0 * fl1_fx * tg_xxxyzz_zzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyzz_zzzz_1[j];

                    tg_xxxxxzzz_xxxx_0[j] = pb_x * tg_xxxxzzz_xxxx_0[j] + wp_x[j] * tg_xxxxzzz_xxxx_1[j] + 2.0 * fl1_fx * tg_xxxzzz_xxxx_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxzzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxxxzzz_xxx_1[j];

                    tg_xxxxxzzz_xxxy_0[j] = pb_x * tg_xxxxzzz_xxxy_0[j] + wp_x[j] * tg_xxxxzzz_xxxy_1[j] + 2.0 * fl1_fx * tg_xxxzzz_xxxy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxzzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxxxzzz_xxy_1[j];

                    tg_xxxxxzzz_xxxz_0[j] = pb_x * tg_xxxxzzz_xxxz_0[j] + wp_x[j] * tg_xxxxzzz_xxxz_1[j] + 2.0 * fl1_fx * tg_xxxzzz_xxxz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxzzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxxxzzz_xxz_1[j];

                    tg_xxxxxzzz_xxyy_0[j] = pb_x * tg_xxxxzzz_xxyy_0[j] + wp_x[j] * tg_xxxxzzz_xxyy_1[j] + 2.0 * fl1_fx * tg_xxxzzz_xxyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxzzz_xxyy_1[j] + fl1_fxn * tg_xxxxzzz_xyy_1[j];

                    tg_xxxxxzzz_xxyz_0[j] = pb_x * tg_xxxxzzz_xxyz_0[j] + wp_x[j] * tg_xxxxzzz_xxyz_1[j] + 2.0 * fl1_fx * tg_xxxzzz_xxyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxzzz_xxyz_1[j] + fl1_fxn * tg_xxxxzzz_xyz_1[j];

                    tg_xxxxxzzz_xxzz_0[j] = pb_x * tg_xxxxzzz_xxzz_0[j] + wp_x[j] * tg_xxxxzzz_xxzz_1[j] + 2.0 * fl1_fx * tg_xxxzzz_xxzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxzzz_xxzz_1[j] + fl1_fxn * tg_xxxxzzz_xzz_1[j];

                    tg_xxxxxzzz_xyyy_0[j] = pb_x * tg_xxxxzzz_xyyy_0[j] + wp_x[j] * tg_xxxxzzz_xyyy_1[j] + 2.0 * fl1_fx * tg_xxxzzz_xyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxzzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxxxzzz_yyy_1[j];

                    tg_xxxxxzzz_xyyz_0[j] = pb_x * tg_xxxxzzz_xyyz_0[j] + wp_x[j] * tg_xxxxzzz_xyyz_1[j] + 2.0 * fl1_fx * tg_xxxzzz_xyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxzzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxxxzzz_yyz_1[j];

                    tg_xxxxxzzz_xyzz_0[j] = pb_x * tg_xxxxzzz_xyzz_0[j] + wp_x[j] * tg_xxxxzzz_xyzz_1[j] + 2.0 * fl1_fx * tg_xxxzzz_xyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxxxzzz_yzz_1[j];

                    tg_xxxxxzzz_xzzz_0[j] = pb_x * tg_xxxxzzz_xzzz_0[j] + wp_x[j] * tg_xxxxzzz_xzzz_1[j] + 2.0 * fl1_fx * tg_xxxzzz_xzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxzzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxxxzzz_zzz_1[j];

                    tg_xxxxxzzz_yyyy_0[j] = pb_x * tg_xxxxzzz_yyyy_0[j] + wp_x[j] * tg_xxxxzzz_yyyy_1[j] + 2.0 * fl1_fx * tg_xxxzzz_yyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxzzz_yyyy_1[j];

                    tg_xxxxxzzz_yyyz_0[j] = pb_x * tg_xxxxzzz_yyyz_0[j] + wp_x[j] * tg_xxxxzzz_yyyz_1[j] + 2.0 * fl1_fx * tg_xxxzzz_yyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxzzz_yyyz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSLSG_147_195(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (147,195)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {8, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {8, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_7_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_6_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_6_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_7_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_xxxxzzz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 147); 

                auto tg_xxxxzzz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 148); 

                auto tg_xxxxzzz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 149); 

                auto tg_xxxyyyy_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 150); 

                auto tg_xxxyyyy_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 151); 

                auto tg_xxxyyyy_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 152); 

                auto tg_xxxyyyy_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 153); 

                auto tg_xxxyyyy_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 154); 

                auto tg_xxxyyyy_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 155); 

                auto tg_xxxyyyy_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 156); 

                auto tg_xxxyyyy_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 157); 

                auto tg_xxxyyyy_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 158); 

                auto tg_xxxyyyy_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 159); 

                auto tg_xxxyyyy_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 160); 

                auto tg_xxxyyyy_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 161); 

                auto tg_xxxyyyy_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 162); 

                auto tg_xxxyyyy_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 163); 

                auto tg_xxxyyyy_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 164); 

                auto tg_xxxyyyz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 165); 

                auto tg_xxxyyyz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 166); 

                auto tg_xxxyyyz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 167); 

                auto tg_xxxyyyz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 168); 

                auto tg_xxxyyyz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 169); 

                auto tg_xxxyyyz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 170); 

                auto tg_xxxyyyz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 171); 

                auto tg_xxxyyyz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 172); 

                auto tg_xxxyyyz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 173); 

                auto tg_xxxyyyz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 174); 

                auto tg_xxxyyyz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 175); 

                auto tg_xxxyyyz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 176); 

                auto tg_xxxyyyz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 177); 

                auto tg_xxxyyyz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 178); 

                auto tg_xxxyyyz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 179); 

                auto tg_xxxyyzz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 180); 

                auto tg_xxxyyzz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 181); 

                auto tg_xxxyyzz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 182); 

                auto tg_xxxyyzz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 183); 

                auto tg_xxxyyzz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 184); 

                auto tg_xxxyyzz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 185); 

                auto tg_xxxyyzz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 186); 

                auto tg_xxxyyzz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 187); 

                auto tg_xxxyyzz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 188); 

                auto tg_xxxyyzz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 189); 

                auto tg_xxxyyzz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 190); 

                auto tg_xxxyyzz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 191); 

                auto tg_xxxyyzz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 192); 

                auto tg_xxxyyzz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 193); 

                auto tg_xxxyyzz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 194); 

                auto tg_xxxxzzz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 147); 

                auto tg_xxxxzzz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 148); 

                auto tg_xxxxzzz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 149); 

                auto tg_xxxyyyy_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 150); 

                auto tg_xxxyyyy_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 151); 

                auto tg_xxxyyyy_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 152); 

                auto tg_xxxyyyy_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 153); 

                auto tg_xxxyyyy_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 154); 

                auto tg_xxxyyyy_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 155); 

                auto tg_xxxyyyy_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 156); 

                auto tg_xxxyyyy_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 157); 

                auto tg_xxxyyyy_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 158); 

                auto tg_xxxyyyy_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 159); 

                auto tg_xxxyyyy_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 160); 

                auto tg_xxxyyyy_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 161); 

                auto tg_xxxyyyy_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 162); 

                auto tg_xxxyyyy_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 163); 

                auto tg_xxxyyyy_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 164); 

                auto tg_xxxyyyz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 165); 

                auto tg_xxxyyyz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 166); 

                auto tg_xxxyyyz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 167); 

                auto tg_xxxyyyz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 168); 

                auto tg_xxxyyyz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 169); 

                auto tg_xxxyyyz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 170); 

                auto tg_xxxyyyz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 171); 

                auto tg_xxxyyyz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 172); 

                auto tg_xxxyyyz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 173); 

                auto tg_xxxyyyz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 174); 

                auto tg_xxxyyyz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 175); 

                auto tg_xxxyyyz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 176); 

                auto tg_xxxyyyz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 177); 

                auto tg_xxxyyyz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 178); 

                auto tg_xxxyyyz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 179); 

                auto tg_xxxyyzz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 180); 

                auto tg_xxxyyzz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 181); 

                auto tg_xxxyyzz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 182); 

                auto tg_xxxyyzz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 183); 

                auto tg_xxxyyzz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 184); 

                auto tg_xxxyyzz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 185); 

                auto tg_xxxyyzz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 186); 

                auto tg_xxxyyzz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 187); 

                auto tg_xxxyyzz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 188); 

                auto tg_xxxyyzz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 189); 

                auto tg_xxxyyzz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 190); 

                auto tg_xxxyyzz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 191); 

                auto tg_xxxyyzz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 192); 

                auto tg_xxxyyzz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 193); 

                auto tg_xxxyyzz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 194); 

                auto tg_xxxzzz_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 147); 

                auto tg_xxxzzz_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 148); 

                auto tg_xxxzzz_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 149); 

                auto tg_xxyyyy_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 150); 

                auto tg_xxyyyy_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 151); 

                auto tg_xxyyyy_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 152); 

                auto tg_xxyyyy_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 153); 

                auto tg_xxyyyy_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 154); 

                auto tg_xxyyyy_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 155); 

                auto tg_xxyyyy_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 156); 

                auto tg_xxyyyy_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 157); 

                auto tg_xxyyyy_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 158); 

                auto tg_xxyyyy_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 159); 

                auto tg_xxyyyy_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 160); 

                auto tg_xxyyyy_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 161); 

                auto tg_xxyyyy_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 162); 

                auto tg_xxyyyy_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 163); 

                auto tg_xxyyyy_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 164); 

                auto tg_xxyyyz_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 165); 

                auto tg_xxyyyz_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 166); 

                auto tg_xxyyyz_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 167); 

                auto tg_xxyyyz_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 168); 

                auto tg_xxyyyz_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 169); 

                auto tg_xxyyyz_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 170); 

                auto tg_xxyyyz_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 171); 

                auto tg_xxyyyz_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 172); 

                auto tg_xxyyyz_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 173); 

                auto tg_xxyyyz_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 174); 

                auto tg_xxyyyz_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 175); 

                auto tg_xxyyyz_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 176); 

                auto tg_xxyyyz_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 177); 

                auto tg_xxyyyz_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 178); 

                auto tg_xxyyyz_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 179); 

                auto tg_xxyyzz_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 180); 

                auto tg_xxyyzz_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 181); 

                auto tg_xxyyzz_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 182); 

                auto tg_xxyyzz_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 183); 

                auto tg_xxyyzz_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 184); 

                auto tg_xxyyzz_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 185); 

                auto tg_xxyyzz_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 186); 

                auto tg_xxyyzz_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 187); 

                auto tg_xxyyzz_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 188); 

                auto tg_xxyyzz_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 189); 

                auto tg_xxyyzz_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 190); 

                auto tg_xxyyzz_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 191); 

                auto tg_xxyyzz_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 192); 

                auto tg_xxyyzz_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 193); 

                auto tg_xxyyzz_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 194); 

                auto tg_xxxzzz_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 147); 

                auto tg_xxxzzz_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 148); 

                auto tg_xxxzzz_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 149); 

                auto tg_xxyyyy_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 150); 

                auto tg_xxyyyy_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 151); 

                auto tg_xxyyyy_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 152); 

                auto tg_xxyyyy_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 153); 

                auto tg_xxyyyy_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 154); 

                auto tg_xxyyyy_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 155); 

                auto tg_xxyyyy_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 156); 

                auto tg_xxyyyy_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 157); 

                auto tg_xxyyyy_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 158); 

                auto tg_xxyyyy_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 159); 

                auto tg_xxyyyy_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 160); 

                auto tg_xxyyyy_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 161); 

                auto tg_xxyyyy_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 162); 

                auto tg_xxyyyy_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 163); 

                auto tg_xxyyyy_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 164); 

                auto tg_xxyyyz_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 165); 

                auto tg_xxyyyz_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 166); 

                auto tg_xxyyyz_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 167); 

                auto tg_xxyyyz_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 168); 

                auto tg_xxyyyz_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 169); 

                auto tg_xxyyyz_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 170); 

                auto tg_xxyyyz_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 171); 

                auto tg_xxyyyz_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 172); 

                auto tg_xxyyyz_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 173); 

                auto tg_xxyyyz_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 174); 

                auto tg_xxyyyz_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 175); 

                auto tg_xxyyyz_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 176); 

                auto tg_xxyyyz_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 177); 

                auto tg_xxyyyz_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 178); 

                auto tg_xxyyyz_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 179); 

                auto tg_xxyyzz_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 180); 

                auto tg_xxyyzz_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 181); 

                auto tg_xxyyzz_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 182); 

                auto tg_xxyyzz_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 183); 

                auto tg_xxyyzz_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 184); 

                auto tg_xxyyzz_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 185); 

                auto tg_xxyyzz_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 186); 

                auto tg_xxyyzz_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 187); 

                auto tg_xxyyzz_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 188); 

                auto tg_xxyyzz_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 189); 

                auto tg_xxyyzz_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 190); 

                auto tg_xxyyzz_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 191); 

                auto tg_xxyyzz_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 192); 

                auto tg_xxyyzz_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 193); 

                auto tg_xxyyzz_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 194); 

                auto tg_xxxyyyy_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 100); 

                auto tg_xxxyyyy_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 101); 

                auto tg_xxxyyyy_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 102); 

                auto tg_xxxyyyy_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 103); 

                auto tg_xxxyyyy_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 104); 

                auto tg_xxxyyyy_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 105); 

                auto tg_xxxyyyy_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 106); 

                auto tg_xxxyyyy_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 107); 

                auto tg_xxxyyyy_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 108); 

                auto tg_xxxyyyy_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 109); 

                auto tg_xxxyyyz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 110); 

                auto tg_xxxyyyz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 111); 

                auto tg_xxxyyyz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 112); 

                auto tg_xxxyyyz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 113); 

                auto tg_xxxyyyz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 114); 

                auto tg_xxxyyyz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 115); 

                auto tg_xxxyyyz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 116); 

                auto tg_xxxyyyz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 117); 

                auto tg_xxxyyyz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 118); 

                auto tg_xxxyyyz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 119); 

                auto tg_xxxyyzz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 120); 

                auto tg_xxxyyzz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 121); 

                auto tg_xxxyyzz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 122); 

                auto tg_xxxyyzz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 123); 

                auto tg_xxxyyzz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 124); 

                auto tg_xxxyyzz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 125); 

                auto tg_xxxyyzz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 126); 

                auto tg_xxxyyzz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 127); 

                auto tg_xxxyyzz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 128); 

                auto tg_xxxyyzz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 129); 

                // set up pointers to integrals

                auto tg_xxxxxzzz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 147); 

                auto tg_xxxxxzzz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 148); 

                auto tg_xxxxxzzz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 149); 

                auto tg_xxxxyyyy_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 150); 

                auto tg_xxxxyyyy_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 151); 

                auto tg_xxxxyyyy_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 152); 

                auto tg_xxxxyyyy_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 153); 

                auto tg_xxxxyyyy_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 154); 

                auto tg_xxxxyyyy_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 155); 

                auto tg_xxxxyyyy_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 156); 

                auto tg_xxxxyyyy_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 157); 

                auto tg_xxxxyyyy_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 158); 

                auto tg_xxxxyyyy_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 159); 

                auto tg_xxxxyyyy_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 160); 

                auto tg_xxxxyyyy_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 161); 

                auto tg_xxxxyyyy_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 162); 

                auto tg_xxxxyyyy_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 163); 

                auto tg_xxxxyyyy_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 164); 

                auto tg_xxxxyyyz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 165); 

                auto tg_xxxxyyyz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 166); 

                auto tg_xxxxyyyz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 167); 

                auto tg_xxxxyyyz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 168); 

                auto tg_xxxxyyyz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 169); 

                auto tg_xxxxyyyz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 170); 

                auto tg_xxxxyyyz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 171); 

                auto tg_xxxxyyyz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 172); 

                auto tg_xxxxyyyz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 173); 

                auto tg_xxxxyyyz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 174); 

                auto tg_xxxxyyyz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 175); 

                auto tg_xxxxyyyz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 176); 

                auto tg_xxxxyyyz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 177); 

                auto tg_xxxxyyyz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 178); 

                auto tg_xxxxyyyz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 179); 

                auto tg_xxxxyyzz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 180); 

                auto tg_xxxxyyzz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 181); 

                auto tg_xxxxyyzz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 182); 

                auto tg_xxxxyyzz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 183); 

                auto tg_xxxxyyzz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 184); 

                auto tg_xxxxyyzz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 185); 

                auto tg_xxxxyyzz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 186); 

                auto tg_xxxxyyzz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 187); 

                auto tg_xxxxyyzz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 188); 

                auto tg_xxxxyyzz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 189); 

                auto tg_xxxxyyzz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 190); 

                auto tg_xxxxyyzz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 191); 

                auto tg_xxxxyyzz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 192); 

                auto tg_xxxxyyzz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 193); 

                auto tg_xxxxyyzz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 194); 

                // Batch of Integrals (147,195)

                #pragma omp simd aligned(fxn, fza, tg_xxxxxzzz_yyzz_0, tg_xxxxxzzz_yzzz_0, tg_xxxxxzzz_zzzz_0, \
                                         tg_xxxxyyyy_xxxx_0, tg_xxxxyyyy_xxxy_0, tg_xxxxyyyy_xxxz_0, tg_xxxxyyyy_xxyy_0, \
                                         tg_xxxxyyyy_xxyz_0, tg_xxxxyyyy_xxzz_0, tg_xxxxyyyy_xyyy_0, tg_xxxxyyyy_xyyz_0, \
                                         tg_xxxxyyyy_xyzz_0, tg_xxxxyyyy_xzzz_0, tg_xxxxyyyy_yyyy_0, tg_xxxxyyyy_yyyz_0, \
                                         tg_xxxxyyyy_yyzz_0, tg_xxxxyyyy_yzzz_0, tg_xxxxyyyy_zzzz_0, tg_xxxxyyyz_xxxx_0, \
                                         tg_xxxxyyyz_xxxy_0, tg_xxxxyyyz_xxxz_0, tg_xxxxyyyz_xxyy_0, tg_xxxxyyyz_xxyz_0, \
                                         tg_xxxxyyyz_xxzz_0, tg_xxxxyyyz_xyyy_0, tg_xxxxyyyz_xyyz_0, tg_xxxxyyyz_xyzz_0, \
                                         tg_xxxxyyyz_xzzz_0, tg_xxxxyyyz_yyyy_0, tg_xxxxyyyz_yyyz_0, tg_xxxxyyyz_yyzz_0, \
                                         tg_xxxxyyyz_yzzz_0, tg_xxxxyyyz_zzzz_0, tg_xxxxyyzz_xxxx_0, tg_xxxxyyzz_xxxy_0, \
                                         tg_xxxxyyzz_xxxz_0, tg_xxxxyyzz_xxyy_0, tg_xxxxyyzz_xxyz_0, tg_xxxxyyzz_xxzz_0, \
                                         tg_xxxxyyzz_xyyy_0, tg_xxxxyyzz_xyyz_0, tg_xxxxyyzz_xyzz_0, tg_xxxxyyzz_xzzz_0, \
                                         tg_xxxxyyzz_yyyy_0, tg_xxxxyyzz_yyyz_0, tg_xxxxyyzz_yyzz_0, tg_xxxxyyzz_yzzz_0, \
                                         tg_xxxxyyzz_zzzz_0, tg_xxxxzzz_yyzz_0, tg_xxxxzzz_yyzz_1, tg_xxxxzzz_yzzz_0, \
                                         tg_xxxxzzz_yzzz_1, tg_xxxxzzz_zzzz_0, tg_xxxxzzz_zzzz_1, tg_xxxyyyy_xxx_1, \
                                         tg_xxxyyyy_xxxx_0, tg_xxxyyyy_xxxx_1, tg_xxxyyyy_xxxy_0, tg_xxxyyyy_xxxy_1, \
                                         tg_xxxyyyy_xxxz_0, tg_xxxyyyy_xxxz_1, tg_xxxyyyy_xxy_1, tg_xxxyyyy_xxyy_0, \
                                         tg_xxxyyyy_xxyy_1, tg_xxxyyyy_xxyz_0, tg_xxxyyyy_xxyz_1, tg_xxxyyyy_xxz_1, \
                                         tg_xxxyyyy_xxzz_0, tg_xxxyyyy_xxzz_1, tg_xxxyyyy_xyy_1, tg_xxxyyyy_xyyy_0, \
                                         tg_xxxyyyy_xyyy_1, tg_xxxyyyy_xyyz_0, tg_xxxyyyy_xyyz_1, tg_xxxyyyy_xyz_1, \
                                         tg_xxxyyyy_xyzz_0, tg_xxxyyyy_xyzz_1, tg_xxxyyyy_xzz_1, tg_xxxyyyy_xzzz_0, \
                                         tg_xxxyyyy_xzzz_1, tg_xxxyyyy_yyy_1, tg_xxxyyyy_yyyy_0, tg_xxxyyyy_yyyy_1, \
                                         tg_xxxyyyy_yyyz_0, tg_xxxyyyy_yyyz_1, tg_xxxyyyy_yyz_1, tg_xxxyyyy_yyzz_0, \
                                         tg_xxxyyyy_yyzz_1, tg_xxxyyyy_yzz_1, tg_xxxyyyy_yzzz_0, tg_xxxyyyy_yzzz_1, \
                                         tg_xxxyyyy_zzz_1, tg_xxxyyyy_zzzz_0, tg_xxxyyyy_zzzz_1, tg_xxxyyyz_xxx_1, \
                                         tg_xxxyyyz_xxxx_0, tg_xxxyyyz_xxxx_1, tg_xxxyyyz_xxxy_0, tg_xxxyyyz_xxxy_1, \
                                         tg_xxxyyyz_xxxz_0, tg_xxxyyyz_xxxz_1, tg_xxxyyyz_xxy_1, tg_xxxyyyz_xxyy_0, \
                                         tg_xxxyyyz_xxyy_1, tg_xxxyyyz_xxyz_0, tg_xxxyyyz_xxyz_1, tg_xxxyyyz_xxz_1, \
                                         tg_xxxyyyz_xxzz_0, tg_xxxyyyz_xxzz_1, tg_xxxyyyz_xyy_1, tg_xxxyyyz_xyyy_0, \
                                         tg_xxxyyyz_xyyy_1, tg_xxxyyyz_xyyz_0, tg_xxxyyyz_xyyz_1, tg_xxxyyyz_xyz_1, \
                                         tg_xxxyyyz_xyzz_0, tg_xxxyyyz_xyzz_1, tg_xxxyyyz_xzz_1, tg_xxxyyyz_xzzz_0, \
                                         tg_xxxyyyz_xzzz_1, tg_xxxyyyz_yyy_1, tg_xxxyyyz_yyyy_0, tg_xxxyyyz_yyyy_1, \
                                         tg_xxxyyyz_yyyz_0, tg_xxxyyyz_yyyz_1, tg_xxxyyyz_yyz_1, tg_xxxyyyz_yyzz_0, \
                                         tg_xxxyyyz_yyzz_1, tg_xxxyyyz_yzz_1, tg_xxxyyyz_yzzz_0, tg_xxxyyyz_yzzz_1, \
                                         tg_xxxyyyz_zzz_1, tg_xxxyyyz_zzzz_0, tg_xxxyyyz_zzzz_1, tg_xxxyyzz_xxx_1, \
                                         tg_xxxyyzz_xxxx_0, tg_xxxyyzz_xxxx_1, tg_xxxyyzz_xxxy_0, tg_xxxyyzz_xxxy_1, \
                                         tg_xxxyyzz_xxxz_0, tg_xxxyyzz_xxxz_1, tg_xxxyyzz_xxy_1, tg_xxxyyzz_xxyy_0, \
                                         tg_xxxyyzz_xxyy_1, tg_xxxyyzz_xxyz_0, tg_xxxyyzz_xxyz_1, tg_xxxyyzz_xxz_1, \
                                         tg_xxxyyzz_xxzz_0, tg_xxxyyzz_xxzz_1, tg_xxxyyzz_xyy_1, tg_xxxyyzz_xyyy_0, \
                                         tg_xxxyyzz_xyyy_1, tg_xxxyyzz_xyyz_0, tg_xxxyyzz_xyyz_1, tg_xxxyyzz_xyz_1, \
                                         tg_xxxyyzz_xyzz_0, tg_xxxyyzz_xyzz_1, tg_xxxyyzz_xzz_1, tg_xxxyyzz_xzzz_0, \
                                         tg_xxxyyzz_xzzz_1, tg_xxxyyzz_yyy_1, tg_xxxyyzz_yyyy_0, tg_xxxyyzz_yyyy_1, \
                                         tg_xxxyyzz_yyyz_0, tg_xxxyyzz_yyyz_1, tg_xxxyyzz_yyz_1, tg_xxxyyzz_yyzz_0, \
                                         tg_xxxyyzz_yyzz_1, tg_xxxyyzz_yzz_1, tg_xxxyyzz_yzzz_0, tg_xxxyyzz_yzzz_1, \
                                         tg_xxxyyzz_zzz_1, tg_xxxyyzz_zzzz_0, tg_xxxyyzz_zzzz_1, tg_xxxzzz_yyzz_0, \
                                         tg_xxxzzz_yyzz_1, tg_xxxzzz_yzzz_0, tg_xxxzzz_yzzz_1, tg_xxxzzz_zzzz_0, \
                                         tg_xxxzzz_zzzz_1, tg_xxyyyy_xxxx_0, tg_xxyyyy_xxxx_1, tg_xxyyyy_xxxy_0, \
                                         tg_xxyyyy_xxxy_1, tg_xxyyyy_xxxz_0, tg_xxyyyy_xxxz_1, tg_xxyyyy_xxyy_0, \
                                         tg_xxyyyy_xxyy_1, tg_xxyyyy_xxyz_0, tg_xxyyyy_xxyz_1, tg_xxyyyy_xxzz_0, \
                                         tg_xxyyyy_xxzz_1, tg_xxyyyy_xyyy_0, tg_xxyyyy_xyyy_1, tg_xxyyyy_xyyz_0, \
                                         tg_xxyyyy_xyyz_1, tg_xxyyyy_xyzz_0, tg_xxyyyy_xyzz_1, tg_xxyyyy_xzzz_0, \
                                         tg_xxyyyy_xzzz_1, tg_xxyyyy_yyyy_0, tg_xxyyyy_yyyy_1, tg_xxyyyy_yyyz_0, \
                                         tg_xxyyyy_yyyz_1, tg_xxyyyy_yyzz_0, tg_xxyyyy_yyzz_1, tg_xxyyyy_yzzz_0, \
                                         tg_xxyyyy_yzzz_1, tg_xxyyyy_zzzz_0, tg_xxyyyy_zzzz_1, tg_xxyyyz_xxxx_0, \
                                         tg_xxyyyz_xxxx_1, tg_xxyyyz_xxxy_0, tg_xxyyyz_xxxy_1, tg_xxyyyz_xxxz_0, \
                                         tg_xxyyyz_xxxz_1, tg_xxyyyz_xxyy_0, tg_xxyyyz_xxyy_1, tg_xxyyyz_xxyz_0, \
                                         tg_xxyyyz_xxyz_1, tg_xxyyyz_xxzz_0, tg_xxyyyz_xxzz_1, tg_xxyyyz_xyyy_0, \
                                         tg_xxyyyz_xyyy_1, tg_xxyyyz_xyyz_0, tg_xxyyyz_xyyz_1, tg_xxyyyz_xyzz_0, \
                                         tg_xxyyyz_xyzz_1, tg_xxyyyz_xzzz_0, tg_xxyyyz_xzzz_1, tg_xxyyyz_yyyy_0, \
                                         tg_xxyyyz_yyyy_1, tg_xxyyyz_yyyz_0, tg_xxyyyz_yyyz_1, tg_xxyyyz_yyzz_0, \
                                         tg_xxyyyz_yyzz_1, tg_xxyyyz_yzzz_0, tg_xxyyyz_yzzz_1, tg_xxyyyz_zzzz_0, \
                                         tg_xxyyyz_zzzz_1, tg_xxyyzz_xxxx_0, tg_xxyyzz_xxxx_1, tg_xxyyzz_xxxy_0, \
                                         tg_xxyyzz_xxxy_1, tg_xxyyzz_xxxz_0, tg_xxyyzz_xxxz_1, tg_xxyyzz_xxyy_0, \
                                         tg_xxyyzz_xxyy_1, tg_xxyyzz_xxyz_0, tg_xxyyzz_xxyz_1, tg_xxyyzz_xxzz_0, \
                                         tg_xxyyzz_xxzz_1, tg_xxyyzz_xyyy_0, tg_xxyyzz_xyyy_1, tg_xxyyzz_xyyz_0, \
                                         tg_xxyyzz_xyyz_1, tg_xxyyzz_xyzz_0, tg_xxyyzz_xyzz_1, tg_xxyyzz_xzzz_0, \
                                         tg_xxyyzz_xzzz_1, tg_xxyyzz_yyyy_0, tg_xxyyzz_yyyy_1, tg_xxyyzz_yyyz_0, \
                                         tg_xxyyzz_yyyz_1, tg_xxyyzz_yyzz_0, tg_xxyyzz_yyzz_1, tg_xxyyzz_yzzz_0, \
                                         tg_xxyyzz_yzzz_1, tg_xxyyzz_zzzz_0, tg_xxyyzz_zzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxxxxzzz_yyzz_0[j] = pb_x * tg_xxxxzzz_yyzz_0[j] + wp_x[j] * tg_xxxxzzz_yyzz_1[j] + 2.0 * fl1_fx * tg_xxxzzz_yyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxzzz_yyzz_1[j];

                    tg_xxxxxzzz_yzzz_0[j] = pb_x * tg_xxxxzzz_yzzz_0[j] + wp_x[j] * tg_xxxxzzz_yzzz_1[j] + 2.0 * fl1_fx * tg_xxxzzz_yzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxzzz_yzzz_1[j];

                    tg_xxxxxzzz_zzzz_0[j] = pb_x * tg_xxxxzzz_zzzz_0[j] + wp_x[j] * tg_xxxxzzz_zzzz_1[j] + 2.0 * fl1_fx * tg_xxxzzz_zzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxzzz_zzzz_1[j];

                    tg_xxxxyyyy_xxxx_0[j] = pb_x * tg_xxxyyyy_xxxx_0[j] + wp_x[j] * tg_xxxyyyy_xxxx_1[j] + 1.5 * fl1_fx * tg_xxyyyy_xxxx_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyy_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxxyyyy_xxx_1[j];

                    tg_xxxxyyyy_xxxy_0[j] = pb_x * tg_xxxyyyy_xxxy_0[j] + wp_x[j] * tg_xxxyyyy_xxxy_1[j] + 1.5 * fl1_fx * tg_xxyyyy_xxxy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyy_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxxyyyy_xxy_1[j];

                    tg_xxxxyyyy_xxxz_0[j] = pb_x * tg_xxxyyyy_xxxz_0[j] + wp_x[j] * tg_xxxyyyy_xxxz_1[j] + 1.5 * fl1_fx * tg_xxyyyy_xxxz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyy_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxxyyyy_xxz_1[j];

                    tg_xxxxyyyy_xxyy_0[j] = pb_x * tg_xxxyyyy_xxyy_0[j] + wp_x[j] * tg_xxxyyyy_xxyy_1[j] + 1.5 * fl1_fx * tg_xxyyyy_xxyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyy_xxyy_1[j] + fl1_fxn * tg_xxxyyyy_xyy_1[j];

                    tg_xxxxyyyy_xxyz_0[j] = pb_x * tg_xxxyyyy_xxyz_0[j] + wp_x[j] * tg_xxxyyyy_xxyz_1[j] + 1.5 * fl1_fx * tg_xxyyyy_xxyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyy_xxyz_1[j] + fl1_fxn * tg_xxxyyyy_xyz_1[j];

                    tg_xxxxyyyy_xxzz_0[j] = pb_x * tg_xxxyyyy_xxzz_0[j] + wp_x[j] * tg_xxxyyyy_xxzz_1[j] + 1.5 * fl1_fx * tg_xxyyyy_xxzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyy_xxzz_1[j] + fl1_fxn * tg_xxxyyyy_xzz_1[j];

                    tg_xxxxyyyy_xyyy_0[j] = pb_x * tg_xxxyyyy_xyyy_0[j] + wp_x[j] * tg_xxxyyyy_xyyy_1[j] + 1.5 * fl1_fx * tg_xxyyyy_xyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyy_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxxyyyy_yyy_1[j];

                    tg_xxxxyyyy_xyyz_0[j] = pb_x * tg_xxxyyyy_xyyz_0[j] + wp_x[j] * tg_xxxyyyy_xyyz_1[j] + 1.5 * fl1_fx * tg_xxyyyy_xyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyy_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxxyyyy_yyz_1[j];

                    tg_xxxxyyyy_xyzz_0[j] = pb_x * tg_xxxyyyy_xyzz_0[j] + wp_x[j] * tg_xxxyyyy_xyzz_1[j] + 1.5 * fl1_fx * tg_xxyyyy_xyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyy_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxxyyyy_yzz_1[j];

                    tg_xxxxyyyy_xzzz_0[j] = pb_x * tg_xxxyyyy_xzzz_0[j] + wp_x[j] * tg_xxxyyyy_xzzz_1[j] + 1.5 * fl1_fx * tg_xxyyyy_xzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyy_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxxyyyy_zzz_1[j];

                    tg_xxxxyyyy_yyyy_0[j] = pb_x * tg_xxxyyyy_yyyy_0[j] + wp_x[j] * tg_xxxyyyy_yyyy_1[j] + 1.5 * fl1_fx * tg_xxyyyy_yyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyy_yyyy_1[j];

                    tg_xxxxyyyy_yyyz_0[j] = pb_x * tg_xxxyyyy_yyyz_0[j] + wp_x[j] * tg_xxxyyyy_yyyz_1[j] + 1.5 * fl1_fx * tg_xxyyyy_yyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyy_yyyz_1[j];

                    tg_xxxxyyyy_yyzz_0[j] = pb_x * tg_xxxyyyy_yyzz_0[j] + wp_x[j] * tg_xxxyyyy_yyzz_1[j] + 1.5 * fl1_fx * tg_xxyyyy_yyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyy_yyzz_1[j];

                    tg_xxxxyyyy_yzzz_0[j] = pb_x * tg_xxxyyyy_yzzz_0[j] + wp_x[j] * tg_xxxyyyy_yzzz_1[j] + 1.5 * fl1_fx * tg_xxyyyy_yzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyy_yzzz_1[j];

                    tg_xxxxyyyy_zzzz_0[j] = pb_x * tg_xxxyyyy_zzzz_0[j] + wp_x[j] * tg_xxxyyyy_zzzz_1[j] + 1.5 * fl1_fx * tg_xxyyyy_zzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyy_zzzz_1[j];

                    tg_xxxxyyyz_xxxx_0[j] = pb_x * tg_xxxyyyz_xxxx_0[j] + wp_x[j] * tg_xxxyyyz_xxxx_1[j] + 1.5 * fl1_fx * tg_xxyyyz_xxxx_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxxyyyz_xxx_1[j];

                    tg_xxxxyyyz_xxxy_0[j] = pb_x * tg_xxxyyyz_xxxy_0[j] + wp_x[j] * tg_xxxyyyz_xxxy_1[j] + 1.5 * fl1_fx * tg_xxyyyz_xxxy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxxyyyz_xxy_1[j];

                    tg_xxxxyyyz_xxxz_0[j] = pb_x * tg_xxxyyyz_xxxz_0[j] + wp_x[j] * tg_xxxyyyz_xxxz_1[j] + 1.5 * fl1_fx * tg_xxyyyz_xxxz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxxyyyz_xxz_1[j];

                    tg_xxxxyyyz_xxyy_0[j] = pb_x * tg_xxxyyyz_xxyy_0[j] + wp_x[j] * tg_xxxyyyz_xxyy_1[j] + 1.5 * fl1_fx * tg_xxyyyz_xxyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyz_xxyy_1[j] + fl1_fxn * tg_xxxyyyz_xyy_1[j];

                    tg_xxxxyyyz_xxyz_0[j] = pb_x * tg_xxxyyyz_xxyz_0[j] + wp_x[j] * tg_xxxyyyz_xxyz_1[j] + 1.5 * fl1_fx * tg_xxyyyz_xxyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyz_xxyz_1[j] + fl1_fxn * tg_xxxyyyz_xyz_1[j];

                    tg_xxxxyyyz_xxzz_0[j] = pb_x * tg_xxxyyyz_xxzz_0[j] + wp_x[j] * tg_xxxyyyz_xxzz_1[j] + 1.5 * fl1_fx * tg_xxyyyz_xxzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyz_xxzz_1[j] + fl1_fxn * tg_xxxyyyz_xzz_1[j];

                    tg_xxxxyyyz_xyyy_0[j] = pb_x * tg_xxxyyyz_xyyy_0[j] + wp_x[j] * tg_xxxyyyz_xyyy_1[j] + 1.5 * fl1_fx * tg_xxyyyz_xyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxxyyyz_yyy_1[j];

                    tg_xxxxyyyz_xyyz_0[j] = pb_x * tg_xxxyyyz_xyyz_0[j] + wp_x[j] * tg_xxxyyyz_xyyz_1[j] + 1.5 * fl1_fx * tg_xxyyyz_xyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxxyyyz_yyz_1[j];

                    tg_xxxxyyyz_xyzz_0[j] = pb_x * tg_xxxyyyz_xyzz_0[j] + wp_x[j] * tg_xxxyyyz_xyzz_1[j] + 1.5 * fl1_fx * tg_xxyyyz_xyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxxyyyz_yzz_1[j];

                    tg_xxxxyyyz_xzzz_0[j] = pb_x * tg_xxxyyyz_xzzz_0[j] + wp_x[j] * tg_xxxyyyz_xzzz_1[j] + 1.5 * fl1_fx * tg_xxyyyz_xzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxxyyyz_zzz_1[j];

                    tg_xxxxyyyz_yyyy_0[j] = pb_x * tg_xxxyyyz_yyyy_0[j] + wp_x[j] * tg_xxxyyyz_yyyy_1[j] + 1.5 * fl1_fx * tg_xxyyyz_yyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyz_yyyy_1[j];

                    tg_xxxxyyyz_yyyz_0[j] = pb_x * tg_xxxyyyz_yyyz_0[j] + wp_x[j] * tg_xxxyyyz_yyyz_1[j] + 1.5 * fl1_fx * tg_xxyyyz_yyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyz_yyyz_1[j];

                    tg_xxxxyyyz_yyzz_0[j] = pb_x * tg_xxxyyyz_yyzz_0[j] + wp_x[j] * tg_xxxyyyz_yyzz_1[j] + 1.5 * fl1_fx * tg_xxyyyz_yyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyz_yyzz_1[j];

                    tg_xxxxyyyz_yzzz_0[j] = pb_x * tg_xxxyyyz_yzzz_0[j] + wp_x[j] * tg_xxxyyyz_yzzz_1[j] + 1.5 * fl1_fx * tg_xxyyyz_yzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyz_yzzz_1[j];

                    tg_xxxxyyyz_zzzz_0[j] = pb_x * tg_xxxyyyz_zzzz_0[j] + wp_x[j] * tg_xxxyyyz_zzzz_1[j] + 1.5 * fl1_fx * tg_xxyyyz_zzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyz_zzzz_1[j];

                    tg_xxxxyyzz_xxxx_0[j] = pb_x * tg_xxxyyzz_xxxx_0[j] + wp_x[j] * tg_xxxyyzz_xxxx_1[j] + 1.5 * fl1_fx * tg_xxyyzz_xxxx_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxxyyzz_xxx_1[j];

                    tg_xxxxyyzz_xxxy_0[j] = pb_x * tg_xxxyyzz_xxxy_0[j] + wp_x[j] * tg_xxxyyzz_xxxy_1[j] + 1.5 * fl1_fx * tg_xxyyzz_xxxy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxxyyzz_xxy_1[j];

                    tg_xxxxyyzz_xxxz_0[j] = pb_x * tg_xxxyyzz_xxxz_0[j] + wp_x[j] * tg_xxxyyzz_xxxz_1[j] + 1.5 * fl1_fx * tg_xxyyzz_xxxz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxxyyzz_xxz_1[j];

                    tg_xxxxyyzz_xxyy_0[j] = pb_x * tg_xxxyyzz_xxyy_0[j] + wp_x[j] * tg_xxxyyzz_xxyy_1[j] + 1.5 * fl1_fx * tg_xxyyzz_xxyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyzz_xxyy_1[j] + fl1_fxn * tg_xxxyyzz_xyy_1[j];

                    tg_xxxxyyzz_xxyz_0[j] = pb_x * tg_xxxyyzz_xxyz_0[j] + wp_x[j] * tg_xxxyyzz_xxyz_1[j] + 1.5 * fl1_fx * tg_xxyyzz_xxyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyzz_xxyz_1[j] + fl1_fxn * tg_xxxyyzz_xyz_1[j];

                    tg_xxxxyyzz_xxzz_0[j] = pb_x * tg_xxxyyzz_xxzz_0[j] + wp_x[j] * tg_xxxyyzz_xxzz_1[j] + 1.5 * fl1_fx * tg_xxyyzz_xxzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyzz_xxzz_1[j] + fl1_fxn * tg_xxxyyzz_xzz_1[j];

                    tg_xxxxyyzz_xyyy_0[j] = pb_x * tg_xxxyyzz_xyyy_0[j] + wp_x[j] * tg_xxxyyzz_xyyy_1[j] + 1.5 * fl1_fx * tg_xxyyzz_xyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxxyyzz_yyy_1[j];

                    tg_xxxxyyzz_xyyz_0[j] = pb_x * tg_xxxyyzz_xyyz_0[j] + wp_x[j] * tg_xxxyyzz_xyyz_1[j] + 1.5 * fl1_fx * tg_xxyyzz_xyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxxyyzz_yyz_1[j];

                    tg_xxxxyyzz_xyzz_0[j] = pb_x * tg_xxxyyzz_xyzz_0[j] + wp_x[j] * tg_xxxyyzz_xyzz_1[j] + 1.5 * fl1_fx * tg_xxyyzz_xyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxxyyzz_yzz_1[j];

                    tg_xxxxyyzz_xzzz_0[j] = pb_x * tg_xxxyyzz_xzzz_0[j] + wp_x[j] * tg_xxxyyzz_xzzz_1[j] + 1.5 * fl1_fx * tg_xxyyzz_xzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxxyyzz_zzz_1[j];

                    tg_xxxxyyzz_yyyy_0[j] = pb_x * tg_xxxyyzz_yyyy_0[j] + wp_x[j] * tg_xxxyyzz_yyyy_1[j] + 1.5 * fl1_fx * tg_xxyyzz_yyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyzz_yyyy_1[j];

                    tg_xxxxyyzz_yyyz_0[j] = pb_x * tg_xxxyyzz_yyyz_0[j] + wp_x[j] * tg_xxxyyzz_yyyz_1[j] + 1.5 * fl1_fx * tg_xxyyzz_yyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyzz_yyyz_1[j];

                    tg_xxxxyyzz_yyzz_0[j] = pb_x * tg_xxxyyzz_yyzz_0[j] + wp_x[j] * tg_xxxyyzz_yyzz_1[j] + 1.5 * fl1_fx * tg_xxyyzz_yyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyzz_yyzz_1[j];

                    tg_xxxxyyzz_yzzz_0[j] = pb_x * tg_xxxyyzz_yzzz_0[j] + wp_x[j] * tg_xxxyyzz_yzzz_1[j] + 1.5 * fl1_fx * tg_xxyyzz_yzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyzz_yzzz_1[j];

                    tg_xxxxyyzz_zzzz_0[j] = pb_x * tg_xxxyyzz_zzzz_0[j] + wp_x[j] * tg_xxxyyzz_zzzz_1[j] + 1.5 * fl1_fx * tg_xxyyzz_zzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyzz_zzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSLSG_195_243(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (195,243)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {8, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {8, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_7_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_6_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_6_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_7_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_xxxyzzz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 195); 

                auto tg_xxxyzzz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 196); 

                auto tg_xxxyzzz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 197); 

                auto tg_xxxyzzz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 198); 

                auto tg_xxxyzzz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 199); 

                auto tg_xxxyzzz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 200); 

                auto tg_xxxyzzz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 201); 

                auto tg_xxxyzzz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 202); 

                auto tg_xxxyzzz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 203); 

                auto tg_xxxyzzz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 204); 

                auto tg_xxxyzzz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 205); 

                auto tg_xxxyzzz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 206); 

                auto tg_xxxyzzz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 207); 

                auto tg_xxxyzzz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 208); 

                auto tg_xxxyzzz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 209); 

                auto tg_xxxzzzz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 210); 

                auto tg_xxxzzzz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 211); 

                auto tg_xxxzzzz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 212); 

                auto tg_xxxzzzz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 213); 

                auto tg_xxxzzzz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 214); 

                auto tg_xxxzzzz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 215); 

                auto tg_xxxzzzz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 216); 

                auto tg_xxxzzzz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 217); 

                auto tg_xxxzzzz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 218); 

                auto tg_xxxzzzz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 219); 

                auto tg_xxxzzzz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 220); 

                auto tg_xxxzzzz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 221); 

                auto tg_xxxzzzz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 222); 

                auto tg_xxxzzzz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 223); 

                auto tg_xxxzzzz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 224); 

                auto tg_xxyyyyy_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 225); 

                auto tg_xxyyyyy_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 226); 

                auto tg_xxyyyyy_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 227); 

                auto tg_xxyyyyy_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 228); 

                auto tg_xxyyyyy_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 229); 

                auto tg_xxyyyyy_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 230); 

                auto tg_xxyyyyy_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 231); 

                auto tg_xxyyyyy_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 232); 

                auto tg_xxyyyyy_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 233); 

                auto tg_xxyyyyy_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 234); 

                auto tg_xxyyyyy_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 235); 

                auto tg_xxyyyyy_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 236); 

                auto tg_xxyyyyy_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 237); 

                auto tg_xxyyyyy_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 238); 

                auto tg_xxyyyyy_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 239); 

                auto tg_xxyyyyz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 240); 

                auto tg_xxyyyyz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 241); 

                auto tg_xxyyyyz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 242); 

                auto tg_xxxyzzz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 195); 

                auto tg_xxxyzzz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 196); 

                auto tg_xxxyzzz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 197); 

                auto tg_xxxyzzz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 198); 

                auto tg_xxxyzzz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 199); 

                auto tg_xxxyzzz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 200); 

                auto tg_xxxyzzz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 201); 

                auto tg_xxxyzzz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 202); 

                auto tg_xxxyzzz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 203); 

                auto tg_xxxyzzz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 204); 

                auto tg_xxxyzzz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 205); 

                auto tg_xxxyzzz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 206); 

                auto tg_xxxyzzz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 207); 

                auto tg_xxxyzzz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 208); 

                auto tg_xxxyzzz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 209); 

                auto tg_xxxzzzz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 210); 

                auto tg_xxxzzzz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 211); 

                auto tg_xxxzzzz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 212); 

                auto tg_xxxzzzz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 213); 

                auto tg_xxxzzzz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 214); 

                auto tg_xxxzzzz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 215); 

                auto tg_xxxzzzz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 216); 

                auto tg_xxxzzzz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 217); 

                auto tg_xxxzzzz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 218); 

                auto tg_xxxzzzz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 219); 

                auto tg_xxxzzzz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 220); 

                auto tg_xxxzzzz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 221); 

                auto tg_xxxzzzz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 222); 

                auto tg_xxxzzzz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 223); 

                auto tg_xxxzzzz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 224); 

                auto tg_xxyyyyy_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 225); 

                auto tg_xxyyyyy_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 226); 

                auto tg_xxyyyyy_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 227); 

                auto tg_xxyyyyy_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 228); 

                auto tg_xxyyyyy_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 229); 

                auto tg_xxyyyyy_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 230); 

                auto tg_xxyyyyy_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 231); 

                auto tg_xxyyyyy_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 232); 

                auto tg_xxyyyyy_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 233); 

                auto tg_xxyyyyy_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 234); 

                auto tg_xxyyyyy_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 235); 

                auto tg_xxyyyyy_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 236); 

                auto tg_xxyyyyy_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 237); 

                auto tg_xxyyyyy_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 238); 

                auto tg_xxyyyyy_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 239); 

                auto tg_xxyyyyz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 240); 

                auto tg_xxyyyyz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 241); 

                auto tg_xxyyyyz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 242); 

                auto tg_xxyzzz_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 195); 

                auto tg_xxyzzz_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 196); 

                auto tg_xxyzzz_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 197); 

                auto tg_xxyzzz_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 198); 

                auto tg_xxyzzz_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 199); 

                auto tg_xxyzzz_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 200); 

                auto tg_xxyzzz_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 201); 

                auto tg_xxyzzz_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 202); 

                auto tg_xxyzzz_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 203); 

                auto tg_xxyzzz_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 204); 

                auto tg_xxyzzz_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 205); 

                auto tg_xxyzzz_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 206); 

                auto tg_xxyzzz_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 207); 

                auto tg_xxyzzz_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 208); 

                auto tg_xxyzzz_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 209); 

                auto tg_xxzzzz_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 210); 

                auto tg_xxzzzz_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 211); 

                auto tg_xxzzzz_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 212); 

                auto tg_xxzzzz_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 213); 

                auto tg_xxzzzz_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 214); 

                auto tg_xxzzzz_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 215); 

                auto tg_xxzzzz_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 216); 

                auto tg_xxzzzz_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 217); 

                auto tg_xxzzzz_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 218); 

                auto tg_xxzzzz_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 219); 

                auto tg_xxzzzz_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 220); 

                auto tg_xxzzzz_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 221); 

                auto tg_xxzzzz_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 222); 

                auto tg_xxzzzz_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 223); 

                auto tg_xxzzzz_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 224); 

                auto tg_xyyyyy_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 225); 

                auto tg_xyyyyy_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 226); 

                auto tg_xyyyyy_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 227); 

                auto tg_xyyyyy_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 228); 

                auto tg_xyyyyy_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 229); 

                auto tg_xyyyyy_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 230); 

                auto tg_xyyyyy_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 231); 

                auto tg_xyyyyy_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 232); 

                auto tg_xyyyyy_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 233); 

                auto tg_xyyyyy_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 234); 

                auto tg_xyyyyy_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 235); 

                auto tg_xyyyyy_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 236); 

                auto tg_xyyyyy_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 237); 

                auto tg_xyyyyy_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 238); 

                auto tg_xyyyyy_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 239); 

                auto tg_xyyyyz_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 240); 

                auto tg_xyyyyz_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 241); 

                auto tg_xyyyyz_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 242); 

                auto tg_xxyzzz_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 195); 

                auto tg_xxyzzz_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 196); 

                auto tg_xxyzzz_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 197); 

                auto tg_xxyzzz_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 198); 

                auto tg_xxyzzz_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 199); 

                auto tg_xxyzzz_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 200); 

                auto tg_xxyzzz_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 201); 

                auto tg_xxyzzz_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 202); 

                auto tg_xxyzzz_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 203); 

                auto tg_xxyzzz_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 204); 

                auto tg_xxyzzz_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 205); 

                auto tg_xxyzzz_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 206); 

                auto tg_xxyzzz_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 207); 

                auto tg_xxyzzz_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 208); 

                auto tg_xxyzzz_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 209); 

                auto tg_xxzzzz_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 210); 

                auto tg_xxzzzz_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 211); 

                auto tg_xxzzzz_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 212); 

                auto tg_xxzzzz_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 213); 

                auto tg_xxzzzz_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 214); 

                auto tg_xxzzzz_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 215); 

                auto tg_xxzzzz_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 216); 

                auto tg_xxzzzz_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 217); 

                auto tg_xxzzzz_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 218); 

                auto tg_xxzzzz_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 219); 

                auto tg_xxzzzz_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 220); 

                auto tg_xxzzzz_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 221); 

                auto tg_xxzzzz_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 222); 

                auto tg_xxzzzz_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 223); 

                auto tg_xxzzzz_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 224); 

                auto tg_xyyyyy_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 225); 

                auto tg_xyyyyy_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 226); 

                auto tg_xyyyyy_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 227); 

                auto tg_xyyyyy_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 228); 

                auto tg_xyyyyy_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 229); 

                auto tg_xyyyyy_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 230); 

                auto tg_xyyyyy_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 231); 

                auto tg_xyyyyy_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 232); 

                auto tg_xyyyyy_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 233); 

                auto tg_xyyyyy_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 234); 

                auto tg_xyyyyy_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 235); 

                auto tg_xyyyyy_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 236); 

                auto tg_xyyyyy_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 237); 

                auto tg_xyyyyy_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 238); 

                auto tg_xyyyyy_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 239); 

                auto tg_xyyyyz_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 240); 

                auto tg_xyyyyz_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 241); 

                auto tg_xyyyyz_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 242); 

                auto tg_xxxyzzz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 130); 

                auto tg_xxxyzzz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 131); 

                auto tg_xxxyzzz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 132); 

                auto tg_xxxyzzz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 133); 

                auto tg_xxxyzzz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 134); 

                auto tg_xxxyzzz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 135); 

                auto tg_xxxyzzz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 136); 

                auto tg_xxxyzzz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 137); 

                auto tg_xxxyzzz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 138); 

                auto tg_xxxyzzz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 139); 

                auto tg_xxxzzzz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 140); 

                auto tg_xxxzzzz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 141); 

                auto tg_xxxzzzz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 142); 

                auto tg_xxxzzzz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 143); 

                auto tg_xxxzzzz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 144); 

                auto tg_xxxzzzz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 145); 

                auto tg_xxxzzzz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 146); 

                auto tg_xxxzzzz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 147); 

                auto tg_xxxzzzz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 148); 

                auto tg_xxxzzzz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 149); 

                auto tg_xxyyyyy_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 150); 

                auto tg_xxyyyyy_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 151); 

                auto tg_xxyyyyy_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 152); 

                auto tg_xxyyyyy_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 153); 

                auto tg_xxyyyyy_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 154); 

                auto tg_xxyyyyy_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 155); 

                auto tg_xxyyyyy_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 156); 

                auto tg_xxyyyyy_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 157); 

                auto tg_xxyyyyy_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 158); 

                auto tg_xxyyyyy_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 159); 

                auto tg_xxyyyyz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 160); 

                auto tg_xxyyyyz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 161); 

                auto tg_xxyyyyz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 162); 

                // set up pointers to integrals

                auto tg_xxxxyzzz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 195); 

                auto tg_xxxxyzzz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 196); 

                auto tg_xxxxyzzz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 197); 

                auto tg_xxxxyzzz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 198); 

                auto tg_xxxxyzzz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 199); 

                auto tg_xxxxyzzz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 200); 

                auto tg_xxxxyzzz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 201); 

                auto tg_xxxxyzzz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 202); 

                auto tg_xxxxyzzz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 203); 

                auto tg_xxxxyzzz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 204); 

                auto tg_xxxxyzzz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 205); 

                auto tg_xxxxyzzz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 206); 

                auto tg_xxxxyzzz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 207); 

                auto tg_xxxxyzzz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 208); 

                auto tg_xxxxyzzz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 209); 

                auto tg_xxxxzzzz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 210); 

                auto tg_xxxxzzzz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 211); 

                auto tg_xxxxzzzz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 212); 

                auto tg_xxxxzzzz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 213); 

                auto tg_xxxxzzzz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 214); 

                auto tg_xxxxzzzz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 215); 

                auto tg_xxxxzzzz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 216); 

                auto tg_xxxxzzzz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 217); 

                auto tg_xxxxzzzz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 218); 

                auto tg_xxxxzzzz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 219); 

                auto tg_xxxxzzzz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 220); 

                auto tg_xxxxzzzz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 221); 

                auto tg_xxxxzzzz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 222); 

                auto tg_xxxxzzzz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 223); 

                auto tg_xxxxzzzz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 224); 

                auto tg_xxxyyyyy_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 225); 

                auto tg_xxxyyyyy_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 226); 

                auto tg_xxxyyyyy_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 227); 

                auto tg_xxxyyyyy_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 228); 

                auto tg_xxxyyyyy_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 229); 

                auto tg_xxxyyyyy_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 230); 

                auto tg_xxxyyyyy_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 231); 

                auto tg_xxxyyyyy_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 232); 

                auto tg_xxxyyyyy_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 233); 

                auto tg_xxxyyyyy_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 234); 

                auto tg_xxxyyyyy_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 235); 

                auto tg_xxxyyyyy_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 236); 

                auto tg_xxxyyyyy_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 237); 

                auto tg_xxxyyyyy_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 238); 

                auto tg_xxxyyyyy_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 239); 

                auto tg_xxxyyyyz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 240); 

                auto tg_xxxyyyyz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 241); 

                auto tg_xxxyyyyz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 242); 

                // Batch of Integrals (195,243)

                #pragma omp simd aligned(fxn, fza, tg_xxxxyzzz_xxxx_0, tg_xxxxyzzz_xxxy_0, tg_xxxxyzzz_xxxz_0, \
                                         tg_xxxxyzzz_xxyy_0, tg_xxxxyzzz_xxyz_0, tg_xxxxyzzz_xxzz_0, tg_xxxxyzzz_xyyy_0, \
                                         tg_xxxxyzzz_xyyz_0, tg_xxxxyzzz_xyzz_0, tg_xxxxyzzz_xzzz_0, tg_xxxxyzzz_yyyy_0, \
                                         tg_xxxxyzzz_yyyz_0, tg_xxxxyzzz_yyzz_0, tg_xxxxyzzz_yzzz_0, tg_xxxxyzzz_zzzz_0, \
                                         tg_xxxxzzzz_xxxx_0, tg_xxxxzzzz_xxxy_0, tg_xxxxzzzz_xxxz_0, tg_xxxxzzzz_xxyy_0, \
                                         tg_xxxxzzzz_xxyz_0, tg_xxxxzzzz_xxzz_0, tg_xxxxzzzz_xyyy_0, tg_xxxxzzzz_xyyz_0, \
                                         tg_xxxxzzzz_xyzz_0, tg_xxxxzzzz_xzzz_0, tg_xxxxzzzz_yyyy_0, tg_xxxxzzzz_yyyz_0, \
                                         tg_xxxxzzzz_yyzz_0, tg_xxxxzzzz_yzzz_0, tg_xxxxzzzz_zzzz_0, tg_xxxyyyyy_xxxx_0, \
                                         tg_xxxyyyyy_xxxy_0, tg_xxxyyyyy_xxxz_0, tg_xxxyyyyy_xxyy_0, tg_xxxyyyyy_xxyz_0, \
                                         tg_xxxyyyyy_xxzz_0, tg_xxxyyyyy_xyyy_0, tg_xxxyyyyy_xyyz_0, tg_xxxyyyyy_xyzz_0, \
                                         tg_xxxyyyyy_xzzz_0, tg_xxxyyyyy_yyyy_0, tg_xxxyyyyy_yyyz_0, tg_xxxyyyyy_yyzz_0, \
                                         tg_xxxyyyyy_yzzz_0, tg_xxxyyyyy_zzzz_0, tg_xxxyyyyz_xxxx_0, tg_xxxyyyyz_xxxy_0, \
                                         tg_xxxyyyyz_xxxz_0, tg_xxxyzzz_xxx_1, tg_xxxyzzz_xxxx_0, tg_xxxyzzz_xxxx_1, \
                                         tg_xxxyzzz_xxxy_0, tg_xxxyzzz_xxxy_1, tg_xxxyzzz_xxxz_0, tg_xxxyzzz_xxxz_1, \
                                         tg_xxxyzzz_xxy_1, tg_xxxyzzz_xxyy_0, tg_xxxyzzz_xxyy_1, tg_xxxyzzz_xxyz_0, \
                                         tg_xxxyzzz_xxyz_1, tg_xxxyzzz_xxz_1, tg_xxxyzzz_xxzz_0, tg_xxxyzzz_xxzz_1, \
                                         tg_xxxyzzz_xyy_1, tg_xxxyzzz_xyyy_0, tg_xxxyzzz_xyyy_1, tg_xxxyzzz_xyyz_0, \
                                         tg_xxxyzzz_xyyz_1, tg_xxxyzzz_xyz_1, tg_xxxyzzz_xyzz_0, tg_xxxyzzz_xyzz_1, \
                                         tg_xxxyzzz_xzz_1, tg_xxxyzzz_xzzz_0, tg_xxxyzzz_xzzz_1, tg_xxxyzzz_yyy_1, \
                                         tg_xxxyzzz_yyyy_0, tg_xxxyzzz_yyyy_1, tg_xxxyzzz_yyyz_0, tg_xxxyzzz_yyyz_1, \
                                         tg_xxxyzzz_yyz_1, tg_xxxyzzz_yyzz_0, tg_xxxyzzz_yyzz_1, tg_xxxyzzz_yzz_1, \
                                         tg_xxxyzzz_yzzz_0, tg_xxxyzzz_yzzz_1, tg_xxxyzzz_zzz_1, tg_xxxyzzz_zzzz_0, \
                                         tg_xxxyzzz_zzzz_1, tg_xxxzzzz_xxx_1, tg_xxxzzzz_xxxx_0, tg_xxxzzzz_xxxx_1, \
                                         tg_xxxzzzz_xxxy_0, tg_xxxzzzz_xxxy_1, tg_xxxzzzz_xxxz_0, tg_xxxzzzz_xxxz_1, \
                                         tg_xxxzzzz_xxy_1, tg_xxxzzzz_xxyy_0, tg_xxxzzzz_xxyy_1, tg_xxxzzzz_xxyz_0, \
                                         tg_xxxzzzz_xxyz_1, tg_xxxzzzz_xxz_1, tg_xxxzzzz_xxzz_0, tg_xxxzzzz_xxzz_1, \
                                         tg_xxxzzzz_xyy_1, tg_xxxzzzz_xyyy_0, tg_xxxzzzz_xyyy_1, tg_xxxzzzz_xyyz_0, \
                                         tg_xxxzzzz_xyyz_1, tg_xxxzzzz_xyz_1, tg_xxxzzzz_xyzz_0, tg_xxxzzzz_xyzz_1, \
                                         tg_xxxzzzz_xzz_1, tg_xxxzzzz_xzzz_0, tg_xxxzzzz_xzzz_1, tg_xxxzzzz_yyy_1, \
                                         tg_xxxzzzz_yyyy_0, tg_xxxzzzz_yyyy_1, tg_xxxzzzz_yyyz_0, tg_xxxzzzz_yyyz_1, \
                                         tg_xxxzzzz_yyz_1, tg_xxxzzzz_yyzz_0, tg_xxxzzzz_yyzz_1, tg_xxxzzzz_yzz_1, \
                                         tg_xxxzzzz_yzzz_0, tg_xxxzzzz_yzzz_1, tg_xxxzzzz_zzz_1, tg_xxxzzzz_zzzz_0, \
                                         tg_xxxzzzz_zzzz_1, tg_xxyyyyy_xxx_1, tg_xxyyyyy_xxxx_0, tg_xxyyyyy_xxxx_1, \
                                         tg_xxyyyyy_xxxy_0, tg_xxyyyyy_xxxy_1, tg_xxyyyyy_xxxz_0, tg_xxyyyyy_xxxz_1, \
                                         tg_xxyyyyy_xxy_1, tg_xxyyyyy_xxyy_0, tg_xxyyyyy_xxyy_1, tg_xxyyyyy_xxyz_0, \
                                         tg_xxyyyyy_xxyz_1, tg_xxyyyyy_xxz_1, tg_xxyyyyy_xxzz_0, tg_xxyyyyy_xxzz_1, \
                                         tg_xxyyyyy_xyy_1, tg_xxyyyyy_xyyy_0, tg_xxyyyyy_xyyy_1, tg_xxyyyyy_xyyz_0, \
                                         tg_xxyyyyy_xyyz_1, tg_xxyyyyy_xyz_1, tg_xxyyyyy_xyzz_0, tg_xxyyyyy_xyzz_1, \
                                         tg_xxyyyyy_xzz_1, tg_xxyyyyy_xzzz_0, tg_xxyyyyy_xzzz_1, tg_xxyyyyy_yyy_1, \
                                         tg_xxyyyyy_yyyy_0, tg_xxyyyyy_yyyy_1, tg_xxyyyyy_yyyz_0, tg_xxyyyyy_yyyz_1, \
                                         tg_xxyyyyy_yyz_1, tg_xxyyyyy_yyzz_0, tg_xxyyyyy_yyzz_1, tg_xxyyyyy_yzz_1, \
                                         tg_xxyyyyy_yzzz_0, tg_xxyyyyy_yzzz_1, tg_xxyyyyy_zzz_1, tg_xxyyyyy_zzzz_0, \
                                         tg_xxyyyyy_zzzz_1, tg_xxyyyyz_xxx_1, tg_xxyyyyz_xxxx_0, tg_xxyyyyz_xxxx_1, \
                                         tg_xxyyyyz_xxxy_0, tg_xxyyyyz_xxxy_1, tg_xxyyyyz_xxxz_0, tg_xxyyyyz_xxxz_1, \
                                         tg_xxyyyyz_xxy_1, tg_xxyyyyz_xxz_1, tg_xxyzzz_xxxx_0, tg_xxyzzz_xxxx_1, \
                                         tg_xxyzzz_xxxy_0, tg_xxyzzz_xxxy_1, tg_xxyzzz_xxxz_0, tg_xxyzzz_xxxz_1, \
                                         tg_xxyzzz_xxyy_0, tg_xxyzzz_xxyy_1, tg_xxyzzz_xxyz_0, tg_xxyzzz_xxyz_1, \
                                         tg_xxyzzz_xxzz_0, tg_xxyzzz_xxzz_1, tg_xxyzzz_xyyy_0, tg_xxyzzz_xyyy_1, \
                                         tg_xxyzzz_xyyz_0, tg_xxyzzz_xyyz_1, tg_xxyzzz_xyzz_0, tg_xxyzzz_xyzz_1, \
                                         tg_xxyzzz_xzzz_0, tg_xxyzzz_xzzz_1, tg_xxyzzz_yyyy_0, tg_xxyzzz_yyyy_1, \
                                         tg_xxyzzz_yyyz_0, tg_xxyzzz_yyyz_1, tg_xxyzzz_yyzz_0, tg_xxyzzz_yyzz_1, \
                                         tg_xxyzzz_yzzz_0, tg_xxyzzz_yzzz_1, tg_xxyzzz_zzzz_0, tg_xxyzzz_zzzz_1, \
                                         tg_xxzzzz_xxxx_0, tg_xxzzzz_xxxx_1, tg_xxzzzz_xxxy_0, tg_xxzzzz_xxxy_1, \
                                         tg_xxzzzz_xxxz_0, tg_xxzzzz_xxxz_1, tg_xxzzzz_xxyy_0, tg_xxzzzz_xxyy_1, \
                                         tg_xxzzzz_xxyz_0, tg_xxzzzz_xxyz_1, tg_xxzzzz_xxzz_0, tg_xxzzzz_xxzz_1, \
                                         tg_xxzzzz_xyyy_0, tg_xxzzzz_xyyy_1, tg_xxzzzz_xyyz_0, tg_xxzzzz_xyyz_1, \
                                         tg_xxzzzz_xyzz_0, tg_xxzzzz_xyzz_1, tg_xxzzzz_xzzz_0, tg_xxzzzz_xzzz_1, \
                                         tg_xxzzzz_yyyy_0, tg_xxzzzz_yyyy_1, tg_xxzzzz_yyyz_0, tg_xxzzzz_yyyz_1, \
                                         tg_xxzzzz_yyzz_0, tg_xxzzzz_yyzz_1, tg_xxzzzz_yzzz_0, tg_xxzzzz_yzzz_1, \
                                         tg_xxzzzz_zzzz_0, tg_xxzzzz_zzzz_1, tg_xyyyyy_xxxx_0, tg_xyyyyy_xxxx_1, \
                                         tg_xyyyyy_xxxy_0, tg_xyyyyy_xxxy_1, tg_xyyyyy_xxxz_0, tg_xyyyyy_xxxz_1, \
                                         tg_xyyyyy_xxyy_0, tg_xyyyyy_xxyy_1, tg_xyyyyy_xxyz_0, tg_xyyyyy_xxyz_1, \
                                         tg_xyyyyy_xxzz_0, tg_xyyyyy_xxzz_1, tg_xyyyyy_xyyy_0, tg_xyyyyy_xyyy_1, \
                                         tg_xyyyyy_xyyz_0, tg_xyyyyy_xyyz_1, tg_xyyyyy_xyzz_0, tg_xyyyyy_xyzz_1, \
                                         tg_xyyyyy_xzzz_0, tg_xyyyyy_xzzz_1, tg_xyyyyy_yyyy_0, tg_xyyyyy_yyyy_1, \
                                         tg_xyyyyy_yyyz_0, tg_xyyyyy_yyyz_1, tg_xyyyyy_yyzz_0, tg_xyyyyy_yyzz_1, \
                                         tg_xyyyyy_yzzz_0, tg_xyyyyy_yzzz_1, tg_xyyyyy_zzzz_0, tg_xyyyyy_zzzz_1, \
                                         tg_xyyyyz_xxxx_0, tg_xyyyyz_xxxx_1, tg_xyyyyz_xxxy_0, tg_xyyyyz_xxxy_1, \
                                         tg_xyyyyz_xxxz_0, tg_xyyyyz_xxxz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxxxyzzz_xxxx_0[j] = pb_x * tg_xxxyzzz_xxxx_0[j] + wp_x[j] * tg_xxxyzzz_xxxx_1[j] + 1.5 * fl1_fx * tg_xxyzzz_xxxx_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyzzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxxyzzz_xxx_1[j];

                    tg_xxxxyzzz_xxxy_0[j] = pb_x * tg_xxxyzzz_xxxy_0[j] + wp_x[j] * tg_xxxyzzz_xxxy_1[j] + 1.5 * fl1_fx * tg_xxyzzz_xxxy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyzzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxxyzzz_xxy_1[j];

                    tg_xxxxyzzz_xxxz_0[j] = pb_x * tg_xxxyzzz_xxxz_0[j] + wp_x[j] * tg_xxxyzzz_xxxz_1[j] + 1.5 * fl1_fx * tg_xxyzzz_xxxz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyzzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxxyzzz_xxz_1[j];

                    tg_xxxxyzzz_xxyy_0[j] = pb_x * tg_xxxyzzz_xxyy_0[j] + wp_x[j] * tg_xxxyzzz_xxyy_1[j] + 1.5 * fl1_fx * tg_xxyzzz_xxyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyzzz_xxyy_1[j] + fl1_fxn * tg_xxxyzzz_xyy_1[j];

                    tg_xxxxyzzz_xxyz_0[j] = pb_x * tg_xxxyzzz_xxyz_0[j] + wp_x[j] * tg_xxxyzzz_xxyz_1[j] + 1.5 * fl1_fx * tg_xxyzzz_xxyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyzzz_xxyz_1[j] + fl1_fxn * tg_xxxyzzz_xyz_1[j];

                    tg_xxxxyzzz_xxzz_0[j] = pb_x * tg_xxxyzzz_xxzz_0[j] + wp_x[j] * tg_xxxyzzz_xxzz_1[j] + 1.5 * fl1_fx * tg_xxyzzz_xxzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyzzz_xxzz_1[j] + fl1_fxn * tg_xxxyzzz_xzz_1[j];

                    tg_xxxxyzzz_xyyy_0[j] = pb_x * tg_xxxyzzz_xyyy_0[j] + wp_x[j] * tg_xxxyzzz_xyyy_1[j] + 1.5 * fl1_fx * tg_xxyzzz_xyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyzzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxxyzzz_yyy_1[j];

                    tg_xxxxyzzz_xyyz_0[j] = pb_x * tg_xxxyzzz_xyyz_0[j] + wp_x[j] * tg_xxxyzzz_xyyz_1[j] + 1.5 * fl1_fx * tg_xxyzzz_xyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyzzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxxyzzz_yyz_1[j];

                    tg_xxxxyzzz_xyzz_0[j] = pb_x * tg_xxxyzzz_xyzz_0[j] + wp_x[j] * tg_xxxyzzz_xyzz_1[j] + 1.5 * fl1_fx * tg_xxyzzz_xyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxxyzzz_yzz_1[j];

                    tg_xxxxyzzz_xzzz_0[j] = pb_x * tg_xxxyzzz_xzzz_0[j] + wp_x[j] * tg_xxxyzzz_xzzz_1[j] + 1.5 * fl1_fx * tg_xxyzzz_xzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyzzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxxyzzz_zzz_1[j];

                    tg_xxxxyzzz_yyyy_0[j] = pb_x * tg_xxxyzzz_yyyy_0[j] + wp_x[j] * tg_xxxyzzz_yyyy_1[j] + 1.5 * fl1_fx * tg_xxyzzz_yyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyzzz_yyyy_1[j];

                    tg_xxxxyzzz_yyyz_0[j] = pb_x * tg_xxxyzzz_yyyz_0[j] + wp_x[j] * tg_xxxyzzz_yyyz_1[j] + 1.5 * fl1_fx * tg_xxyzzz_yyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyzzz_yyyz_1[j];

                    tg_xxxxyzzz_yyzz_0[j] = pb_x * tg_xxxyzzz_yyzz_0[j] + wp_x[j] * tg_xxxyzzz_yyzz_1[j] + 1.5 * fl1_fx * tg_xxyzzz_yyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyzzz_yyzz_1[j];

                    tg_xxxxyzzz_yzzz_0[j] = pb_x * tg_xxxyzzz_yzzz_0[j] + wp_x[j] * tg_xxxyzzz_yzzz_1[j] + 1.5 * fl1_fx * tg_xxyzzz_yzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyzzz_yzzz_1[j];

                    tg_xxxxyzzz_zzzz_0[j] = pb_x * tg_xxxyzzz_zzzz_0[j] + wp_x[j] * tg_xxxyzzz_zzzz_1[j] + 1.5 * fl1_fx * tg_xxyzzz_zzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyzzz_zzzz_1[j];

                    tg_xxxxzzzz_xxxx_0[j] = pb_x * tg_xxxzzzz_xxxx_0[j] + wp_x[j] * tg_xxxzzzz_xxxx_1[j] + 1.5 * fl1_fx * tg_xxzzzz_xxxx_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzzzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxxzzzz_xxx_1[j];

                    tg_xxxxzzzz_xxxy_0[j] = pb_x * tg_xxxzzzz_xxxy_0[j] + wp_x[j] * tg_xxxzzzz_xxxy_1[j] + 1.5 * fl1_fx * tg_xxzzzz_xxxy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzzzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxxzzzz_xxy_1[j];

                    tg_xxxxzzzz_xxxz_0[j] = pb_x * tg_xxxzzzz_xxxz_0[j] + wp_x[j] * tg_xxxzzzz_xxxz_1[j] + 1.5 * fl1_fx * tg_xxzzzz_xxxz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzzzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxxzzzz_xxz_1[j];

                    tg_xxxxzzzz_xxyy_0[j] = pb_x * tg_xxxzzzz_xxyy_0[j] + wp_x[j] * tg_xxxzzzz_xxyy_1[j] + 1.5 * fl1_fx * tg_xxzzzz_xxyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzzzz_xxyy_1[j] + fl1_fxn * tg_xxxzzzz_xyy_1[j];

                    tg_xxxxzzzz_xxyz_0[j] = pb_x * tg_xxxzzzz_xxyz_0[j] + wp_x[j] * tg_xxxzzzz_xxyz_1[j] + 1.5 * fl1_fx * tg_xxzzzz_xxyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzzzz_xxyz_1[j] + fl1_fxn * tg_xxxzzzz_xyz_1[j];

                    tg_xxxxzzzz_xxzz_0[j] = pb_x * tg_xxxzzzz_xxzz_0[j] + wp_x[j] * tg_xxxzzzz_xxzz_1[j] + 1.5 * fl1_fx * tg_xxzzzz_xxzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzzzz_xxzz_1[j] + fl1_fxn * tg_xxxzzzz_xzz_1[j];

                    tg_xxxxzzzz_xyyy_0[j] = pb_x * tg_xxxzzzz_xyyy_0[j] + wp_x[j] * tg_xxxzzzz_xyyy_1[j] + 1.5 * fl1_fx * tg_xxzzzz_xyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzzzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxxzzzz_yyy_1[j];

                    tg_xxxxzzzz_xyyz_0[j] = pb_x * tg_xxxzzzz_xyyz_0[j] + wp_x[j] * tg_xxxzzzz_xyyz_1[j] + 1.5 * fl1_fx * tg_xxzzzz_xyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzzzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxxzzzz_yyz_1[j];

                    tg_xxxxzzzz_xyzz_0[j] = pb_x * tg_xxxzzzz_xyzz_0[j] + wp_x[j] * tg_xxxzzzz_xyzz_1[j] + 1.5 * fl1_fx * tg_xxzzzz_xyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxxzzzz_yzz_1[j];

                    tg_xxxxzzzz_xzzz_0[j] = pb_x * tg_xxxzzzz_xzzz_0[j] + wp_x[j] * tg_xxxzzzz_xzzz_1[j] + 1.5 * fl1_fx * tg_xxzzzz_xzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzzzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxxzzzz_zzz_1[j];

                    tg_xxxxzzzz_yyyy_0[j] = pb_x * tg_xxxzzzz_yyyy_0[j] + wp_x[j] * tg_xxxzzzz_yyyy_1[j] + 1.5 * fl1_fx * tg_xxzzzz_yyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzzzz_yyyy_1[j];

                    tg_xxxxzzzz_yyyz_0[j] = pb_x * tg_xxxzzzz_yyyz_0[j] + wp_x[j] * tg_xxxzzzz_yyyz_1[j] + 1.5 * fl1_fx * tg_xxzzzz_yyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzzzz_yyyz_1[j];

                    tg_xxxxzzzz_yyzz_0[j] = pb_x * tg_xxxzzzz_yyzz_0[j] + wp_x[j] * tg_xxxzzzz_yyzz_1[j] + 1.5 * fl1_fx * tg_xxzzzz_yyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzzzz_yyzz_1[j];

                    tg_xxxxzzzz_yzzz_0[j] = pb_x * tg_xxxzzzz_yzzz_0[j] + wp_x[j] * tg_xxxzzzz_yzzz_1[j] + 1.5 * fl1_fx * tg_xxzzzz_yzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzzzz_yzzz_1[j];

                    tg_xxxxzzzz_zzzz_0[j] = pb_x * tg_xxxzzzz_zzzz_0[j] + wp_x[j] * tg_xxxzzzz_zzzz_1[j] + 1.5 * fl1_fx * tg_xxzzzz_zzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzzzz_zzzz_1[j];

                    tg_xxxyyyyy_xxxx_0[j] = pb_x * tg_xxyyyyy_xxxx_0[j] + wp_x[j] * tg_xxyyyyy_xxxx_1[j] + fl1_fx * tg_xyyyyy_xxxx_0[j] - fl1_fx * fl1_fza * tg_xyyyyy_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxyyyyy_xxx_1[j];

                    tg_xxxyyyyy_xxxy_0[j] = pb_x * tg_xxyyyyy_xxxy_0[j] + wp_x[j] * tg_xxyyyyy_xxxy_1[j] + fl1_fx * tg_xyyyyy_xxxy_0[j] - fl1_fx * fl1_fza * tg_xyyyyy_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxyyyyy_xxy_1[j];

                    tg_xxxyyyyy_xxxz_0[j] = pb_x * tg_xxyyyyy_xxxz_0[j] + wp_x[j] * tg_xxyyyyy_xxxz_1[j] + fl1_fx * tg_xyyyyy_xxxz_0[j] - fl1_fx * fl1_fza * tg_xyyyyy_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxyyyyy_xxz_1[j];

                    tg_xxxyyyyy_xxyy_0[j] = pb_x * tg_xxyyyyy_xxyy_0[j] + wp_x[j] * tg_xxyyyyy_xxyy_1[j] + fl1_fx * tg_xyyyyy_xxyy_0[j] - fl1_fx * fl1_fza * tg_xyyyyy_xxyy_1[j] + fl1_fxn * tg_xxyyyyy_xyy_1[j];

                    tg_xxxyyyyy_xxyz_0[j] = pb_x * tg_xxyyyyy_xxyz_0[j] + wp_x[j] * tg_xxyyyyy_xxyz_1[j] + fl1_fx * tg_xyyyyy_xxyz_0[j] - fl1_fx * fl1_fza * tg_xyyyyy_xxyz_1[j] + fl1_fxn * tg_xxyyyyy_xyz_1[j];

                    tg_xxxyyyyy_xxzz_0[j] = pb_x * tg_xxyyyyy_xxzz_0[j] + wp_x[j] * tg_xxyyyyy_xxzz_1[j] + fl1_fx * tg_xyyyyy_xxzz_0[j] - fl1_fx * fl1_fza * tg_xyyyyy_xxzz_1[j] + fl1_fxn * tg_xxyyyyy_xzz_1[j];

                    tg_xxxyyyyy_xyyy_0[j] = pb_x * tg_xxyyyyy_xyyy_0[j] + wp_x[j] * tg_xxyyyyy_xyyy_1[j] + fl1_fx * tg_xyyyyy_xyyy_0[j] - fl1_fx * fl1_fza * tg_xyyyyy_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxyyyyy_yyy_1[j];

                    tg_xxxyyyyy_xyyz_0[j] = pb_x * tg_xxyyyyy_xyyz_0[j] + wp_x[j] * tg_xxyyyyy_xyyz_1[j] + fl1_fx * tg_xyyyyy_xyyz_0[j] - fl1_fx * fl1_fza * tg_xyyyyy_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxyyyyy_yyz_1[j];

                    tg_xxxyyyyy_xyzz_0[j] = pb_x * tg_xxyyyyy_xyzz_0[j] + wp_x[j] * tg_xxyyyyy_xyzz_1[j] + fl1_fx * tg_xyyyyy_xyzz_0[j] - fl1_fx * fl1_fza * tg_xyyyyy_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxyyyyy_yzz_1[j];

                    tg_xxxyyyyy_xzzz_0[j] = pb_x * tg_xxyyyyy_xzzz_0[j] + wp_x[j] * tg_xxyyyyy_xzzz_1[j] + fl1_fx * tg_xyyyyy_xzzz_0[j] - fl1_fx * fl1_fza * tg_xyyyyy_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxyyyyy_zzz_1[j];

                    tg_xxxyyyyy_yyyy_0[j] = pb_x * tg_xxyyyyy_yyyy_0[j] + wp_x[j] * tg_xxyyyyy_yyyy_1[j] + fl1_fx * tg_xyyyyy_yyyy_0[j] - fl1_fx * fl1_fza * tg_xyyyyy_yyyy_1[j];

                    tg_xxxyyyyy_yyyz_0[j] = pb_x * tg_xxyyyyy_yyyz_0[j] + wp_x[j] * tg_xxyyyyy_yyyz_1[j] + fl1_fx * tg_xyyyyy_yyyz_0[j] - fl1_fx * fl1_fza * tg_xyyyyy_yyyz_1[j];

                    tg_xxxyyyyy_yyzz_0[j] = pb_x * tg_xxyyyyy_yyzz_0[j] + wp_x[j] * tg_xxyyyyy_yyzz_1[j] + fl1_fx * tg_xyyyyy_yyzz_0[j] - fl1_fx * fl1_fza * tg_xyyyyy_yyzz_1[j];

                    tg_xxxyyyyy_yzzz_0[j] = pb_x * tg_xxyyyyy_yzzz_0[j] + wp_x[j] * tg_xxyyyyy_yzzz_1[j] + fl1_fx * tg_xyyyyy_yzzz_0[j] - fl1_fx * fl1_fza * tg_xyyyyy_yzzz_1[j];

                    tg_xxxyyyyy_zzzz_0[j] = pb_x * tg_xxyyyyy_zzzz_0[j] + wp_x[j] * tg_xxyyyyy_zzzz_1[j] + fl1_fx * tg_xyyyyy_zzzz_0[j] - fl1_fx * fl1_fza * tg_xyyyyy_zzzz_1[j];

                    tg_xxxyyyyz_xxxx_0[j] = pb_x * tg_xxyyyyz_xxxx_0[j] + wp_x[j] * tg_xxyyyyz_xxxx_1[j] + fl1_fx * tg_xyyyyz_xxxx_0[j] - fl1_fx * fl1_fza * tg_xyyyyz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxyyyyz_xxx_1[j];

                    tg_xxxyyyyz_xxxy_0[j] = pb_x * tg_xxyyyyz_xxxy_0[j] + wp_x[j] * tg_xxyyyyz_xxxy_1[j] + fl1_fx * tg_xyyyyz_xxxy_0[j] - fl1_fx * fl1_fza * tg_xyyyyz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxyyyyz_xxy_1[j];

                    tg_xxxyyyyz_xxxz_0[j] = pb_x * tg_xxyyyyz_xxxz_0[j] + wp_x[j] * tg_xxyyyyz_xxxz_1[j] + fl1_fx * tg_xyyyyz_xxxz_0[j] - fl1_fx * fl1_fza * tg_xyyyyz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxyyyyz_xxz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSLSG_243_291(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (243,291)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {8, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {8, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_7_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_6_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_6_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_7_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_xxyyyyz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 243); 

                auto tg_xxyyyyz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 244); 

                auto tg_xxyyyyz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 245); 

                auto tg_xxyyyyz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 246); 

                auto tg_xxyyyyz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 247); 

                auto tg_xxyyyyz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 248); 

                auto tg_xxyyyyz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 249); 

                auto tg_xxyyyyz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 250); 

                auto tg_xxyyyyz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 251); 

                auto tg_xxyyyyz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 252); 

                auto tg_xxyyyyz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 253); 

                auto tg_xxyyyyz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 254); 

                auto tg_xxyyyzz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 255); 

                auto tg_xxyyyzz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 256); 

                auto tg_xxyyyzz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 257); 

                auto tg_xxyyyzz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 258); 

                auto tg_xxyyyzz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 259); 

                auto tg_xxyyyzz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 260); 

                auto tg_xxyyyzz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 261); 

                auto tg_xxyyyzz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 262); 

                auto tg_xxyyyzz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 263); 

                auto tg_xxyyyzz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 264); 

                auto tg_xxyyyzz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 265); 

                auto tg_xxyyyzz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 266); 

                auto tg_xxyyyzz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 267); 

                auto tg_xxyyyzz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 268); 

                auto tg_xxyyyzz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 269); 

                auto tg_xxyyzzz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 270); 

                auto tg_xxyyzzz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 271); 

                auto tg_xxyyzzz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 272); 

                auto tg_xxyyzzz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 273); 

                auto tg_xxyyzzz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 274); 

                auto tg_xxyyzzz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 275); 

                auto tg_xxyyzzz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 276); 

                auto tg_xxyyzzz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 277); 

                auto tg_xxyyzzz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 278); 

                auto tg_xxyyzzz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 279); 

                auto tg_xxyyzzz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 280); 

                auto tg_xxyyzzz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 281); 

                auto tg_xxyyzzz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 282); 

                auto tg_xxyyzzz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 283); 

                auto tg_xxyyzzz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 284); 

                auto tg_xxyzzzz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 285); 

                auto tg_xxyzzzz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 286); 

                auto tg_xxyzzzz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 287); 

                auto tg_xxyzzzz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 288); 

                auto tg_xxyzzzz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 289); 

                auto tg_xxyzzzz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 290); 

                auto tg_xxyyyyz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 243); 

                auto tg_xxyyyyz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 244); 

                auto tg_xxyyyyz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 245); 

                auto tg_xxyyyyz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 246); 

                auto tg_xxyyyyz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 247); 

                auto tg_xxyyyyz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 248); 

                auto tg_xxyyyyz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 249); 

                auto tg_xxyyyyz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 250); 

                auto tg_xxyyyyz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 251); 

                auto tg_xxyyyyz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 252); 

                auto tg_xxyyyyz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 253); 

                auto tg_xxyyyyz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 254); 

                auto tg_xxyyyzz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 255); 

                auto tg_xxyyyzz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 256); 

                auto tg_xxyyyzz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 257); 

                auto tg_xxyyyzz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 258); 

                auto tg_xxyyyzz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 259); 

                auto tg_xxyyyzz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 260); 

                auto tg_xxyyyzz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 261); 

                auto tg_xxyyyzz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 262); 

                auto tg_xxyyyzz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 263); 

                auto tg_xxyyyzz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 264); 

                auto tg_xxyyyzz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 265); 

                auto tg_xxyyyzz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 266); 

                auto tg_xxyyyzz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 267); 

                auto tg_xxyyyzz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 268); 

                auto tg_xxyyyzz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 269); 

                auto tg_xxyyzzz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 270); 

                auto tg_xxyyzzz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 271); 

                auto tg_xxyyzzz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 272); 

                auto tg_xxyyzzz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 273); 

                auto tg_xxyyzzz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 274); 

                auto tg_xxyyzzz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 275); 

                auto tg_xxyyzzz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 276); 

                auto tg_xxyyzzz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 277); 

                auto tg_xxyyzzz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 278); 

                auto tg_xxyyzzz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 279); 

                auto tg_xxyyzzz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 280); 

                auto tg_xxyyzzz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 281); 

                auto tg_xxyyzzz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 282); 

                auto tg_xxyyzzz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 283); 

                auto tg_xxyyzzz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 284); 

                auto tg_xxyzzzz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 285); 

                auto tg_xxyzzzz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 286); 

                auto tg_xxyzzzz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 287); 

                auto tg_xxyzzzz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 288); 

                auto tg_xxyzzzz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 289); 

                auto tg_xxyzzzz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 290); 

                auto tg_xyyyyz_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 243); 

                auto tg_xyyyyz_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 244); 

                auto tg_xyyyyz_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 245); 

                auto tg_xyyyyz_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 246); 

                auto tg_xyyyyz_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 247); 

                auto tg_xyyyyz_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 248); 

                auto tg_xyyyyz_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 249); 

                auto tg_xyyyyz_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 250); 

                auto tg_xyyyyz_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 251); 

                auto tg_xyyyyz_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 252); 

                auto tg_xyyyyz_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 253); 

                auto tg_xyyyyz_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 254); 

                auto tg_xyyyzz_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 255); 

                auto tg_xyyyzz_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 256); 

                auto tg_xyyyzz_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 257); 

                auto tg_xyyyzz_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 258); 

                auto tg_xyyyzz_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 259); 

                auto tg_xyyyzz_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 260); 

                auto tg_xyyyzz_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 261); 

                auto tg_xyyyzz_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 262); 

                auto tg_xyyyzz_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 263); 

                auto tg_xyyyzz_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 264); 

                auto tg_xyyyzz_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 265); 

                auto tg_xyyyzz_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 266); 

                auto tg_xyyyzz_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 267); 

                auto tg_xyyyzz_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 268); 

                auto tg_xyyyzz_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 269); 

                auto tg_xyyzzz_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 270); 

                auto tg_xyyzzz_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 271); 

                auto tg_xyyzzz_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 272); 

                auto tg_xyyzzz_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 273); 

                auto tg_xyyzzz_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 274); 

                auto tg_xyyzzz_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 275); 

                auto tg_xyyzzz_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 276); 

                auto tg_xyyzzz_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 277); 

                auto tg_xyyzzz_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 278); 

                auto tg_xyyzzz_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 279); 

                auto tg_xyyzzz_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 280); 

                auto tg_xyyzzz_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 281); 

                auto tg_xyyzzz_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 282); 

                auto tg_xyyzzz_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 283); 

                auto tg_xyyzzz_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 284); 

                auto tg_xyzzzz_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 285); 

                auto tg_xyzzzz_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 286); 

                auto tg_xyzzzz_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 287); 

                auto tg_xyzzzz_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 288); 

                auto tg_xyzzzz_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 289); 

                auto tg_xyzzzz_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 290); 

                auto tg_xyyyyz_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 243); 

                auto tg_xyyyyz_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 244); 

                auto tg_xyyyyz_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 245); 

                auto tg_xyyyyz_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 246); 

                auto tg_xyyyyz_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 247); 

                auto tg_xyyyyz_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 248); 

                auto tg_xyyyyz_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 249); 

                auto tg_xyyyyz_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 250); 

                auto tg_xyyyyz_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 251); 

                auto tg_xyyyyz_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 252); 

                auto tg_xyyyyz_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 253); 

                auto tg_xyyyyz_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 254); 

                auto tg_xyyyzz_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 255); 

                auto tg_xyyyzz_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 256); 

                auto tg_xyyyzz_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 257); 

                auto tg_xyyyzz_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 258); 

                auto tg_xyyyzz_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 259); 

                auto tg_xyyyzz_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 260); 

                auto tg_xyyyzz_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 261); 

                auto tg_xyyyzz_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 262); 

                auto tg_xyyyzz_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 263); 

                auto tg_xyyyzz_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 264); 

                auto tg_xyyyzz_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 265); 

                auto tg_xyyyzz_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 266); 

                auto tg_xyyyzz_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 267); 

                auto tg_xyyyzz_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 268); 

                auto tg_xyyyzz_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 269); 

                auto tg_xyyzzz_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 270); 

                auto tg_xyyzzz_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 271); 

                auto tg_xyyzzz_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 272); 

                auto tg_xyyzzz_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 273); 

                auto tg_xyyzzz_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 274); 

                auto tg_xyyzzz_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 275); 

                auto tg_xyyzzz_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 276); 

                auto tg_xyyzzz_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 277); 

                auto tg_xyyzzz_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 278); 

                auto tg_xyyzzz_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 279); 

                auto tg_xyyzzz_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 280); 

                auto tg_xyyzzz_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 281); 

                auto tg_xyyzzz_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 282); 

                auto tg_xyyzzz_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 283); 

                auto tg_xyyzzz_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 284); 

                auto tg_xyzzzz_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 285); 

                auto tg_xyzzzz_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 286); 

                auto tg_xyzzzz_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 287); 

                auto tg_xyzzzz_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 288); 

                auto tg_xyzzzz_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 289); 

                auto tg_xyzzzz_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 290); 

                auto tg_xxyyyyz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 163); 

                auto tg_xxyyyyz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 164); 

                auto tg_xxyyyyz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 165); 

                auto tg_xxyyyyz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 166); 

                auto tg_xxyyyyz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 167); 

                auto tg_xxyyyyz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 168); 

                auto tg_xxyyyyz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 169); 

                auto tg_xxyyyzz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 170); 

                auto tg_xxyyyzz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 171); 

                auto tg_xxyyyzz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 172); 

                auto tg_xxyyyzz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 173); 

                auto tg_xxyyyzz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 174); 

                auto tg_xxyyyzz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 175); 

                auto tg_xxyyyzz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 176); 

                auto tg_xxyyyzz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 177); 

                auto tg_xxyyyzz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 178); 

                auto tg_xxyyyzz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 179); 

                auto tg_xxyyzzz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 180); 

                auto tg_xxyyzzz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 181); 

                auto tg_xxyyzzz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 182); 

                auto tg_xxyyzzz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 183); 

                auto tg_xxyyzzz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 184); 

                auto tg_xxyyzzz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 185); 

                auto tg_xxyyzzz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 186); 

                auto tg_xxyyzzz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 187); 

                auto tg_xxyyzzz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 188); 

                auto tg_xxyyzzz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 189); 

                auto tg_xxyzzzz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 190); 

                auto tg_xxyzzzz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 191); 

                auto tg_xxyzzzz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 192); 

                auto tg_xxyzzzz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 193); 

                auto tg_xxyzzzz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 194); 

                auto tg_xxyzzzz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 195); 

                // set up pointers to integrals

                auto tg_xxxyyyyz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 243); 

                auto tg_xxxyyyyz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 244); 

                auto tg_xxxyyyyz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 245); 

                auto tg_xxxyyyyz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 246); 

                auto tg_xxxyyyyz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 247); 

                auto tg_xxxyyyyz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 248); 

                auto tg_xxxyyyyz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 249); 

                auto tg_xxxyyyyz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 250); 

                auto tg_xxxyyyyz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 251); 

                auto tg_xxxyyyyz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 252); 

                auto tg_xxxyyyyz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 253); 

                auto tg_xxxyyyyz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 254); 

                auto tg_xxxyyyzz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 255); 

                auto tg_xxxyyyzz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 256); 

                auto tg_xxxyyyzz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 257); 

                auto tg_xxxyyyzz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 258); 

                auto tg_xxxyyyzz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 259); 

                auto tg_xxxyyyzz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 260); 

                auto tg_xxxyyyzz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 261); 

                auto tg_xxxyyyzz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 262); 

                auto tg_xxxyyyzz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 263); 

                auto tg_xxxyyyzz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 264); 

                auto tg_xxxyyyzz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 265); 

                auto tg_xxxyyyzz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 266); 

                auto tg_xxxyyyzz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 267); 

                auto tg_xxxyyyzz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 268); 

                auto tg_xxxyyyzz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 269); 

                auto tg_xxxyyzzz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 270); 

                auto tg_xxxyyzzz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 271); 

                auto tg_xxxyyzzz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 272); 

                auto tg_xxxyyzzz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 273); 

                auto tg_xxxyyzzz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 274); 

                auto tg_xxxyyzzz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 275); 

                auto tg_xxxyyzzz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 276); 

                auto tg_xxxyyzzz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 277); 

                auto tg_xxxyyzzz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 278); 

                auto tg_xxxyyzzz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 279); 

                auto tg_xxxyyzzz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 280); 

                auto tg_xxxyyzzz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 281); 

                auto tg_xxxyyzzz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 282); 

                auto tg_xxxyyzzz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 283); 

                auto tg_xxxyyzzz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 284); 

                auto tg_xxxyzzzz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 285); 

                auto tg_xxxyzzzz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 286); 

                auto tg_xxxyzzzz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 287); 

                auto tg_xxxyzzzz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 288); 

                auto tg_xxxyzzzz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 289); 

                auto tg_xxxyzzzz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 290); 

                // Batch of Integrals (243,291)

                #pragma omp simd aligned(fxn, fza, tg_xxxyyyyz_xxyy_0, tg_xxxyyyyz_xxyz_0, tg_xxxyyyyz_xxzz_0, \
                                         tg_xxxyyyyz_xyyy_0, tg_xxxyyyyz_xyyz_0, tg_xxxyyyyz_xyzz_0, tg_xxxyyyyz_xzzz_0, \
                                         tg_xxxyyyyz_yyyy_0, tg_xxxyyyyz_yyyz_0, tg_xxxyyyyz_yyzz_0, tg_xxxyyyyz_yzzz_0, \
                                         tg_xxxyyyyz_zzzz_0, tg_xxxyyyzz_xxxx_0, tg_xxxyyyzz_xxxy_0, tg_xxxyyyzz_xxxz_0, \
                                         tg_xxxyyyzz_xxyy_0, tg_xxxyyyzz_xxyz_0, tg_xxxyyyzz_xxzz_0, tg_xxxyyyzz_xyyy_0, \
                                         tg_xxxyyyzz_xyyz_0, tg_xxxyyyzz_xyzz_0, tg_xxxyyyzz_xzzz_0, tg_xxxyyyzz_yyyy_0, \
                                         tg_xxxyyyzz_yyyz_0, tg_xxxyyyzz_yyzz_0, tg_xxxyyyzz_yzzz_0, tg_xxxyyyzz_zzzz_0, \
                                         tg_xxxyyzzz_xxxx_0, tg_xxxyyzzz_xxxy_0, tg_xxxyyzzz_xxxz_0, tg_xxxyyzzz_xxyy_0, \
                                         tg_xxxyyzzz_xxyz_0, tg_xxxyyzzz_xxzz_0, tg_xxxyyzzz_xyyy_0, tg_xxxyyzzz_xyyz_0, \
                                         tg_xxxyyzzz_xyzz_0, tg_xxxyyzzz_xzzz_0, tg_xxxyyzzz_yyyy_0, tg_xxxyyzzz_yyyz_0, \
                                         tg_xxxyyzzz_yyzz_0, tg_xxxyyzzz_yzzz_0, tg_xxxyyzzz_zzzz_0, tg_xxxyzzzz_xxxx_0, \
                                         tg_xxxyzzzz_xxxy_0, tg_xxxyzzzz_xxxz_0, tg_xxxyzzzz_xxyy_0, tg_xxxyzzzz_xxyz_0, \
                                         tg_xxxyzzzz_xxzz_0, tg_xxyyyyz_xxyy_0, tg_xxyyyyz_xxyy_1, tg_xxyyyyz_xxyz_0, \
                                         tg_xxyyyyz_xxyz_1, tg_xxyyyyz_xxzz_0, tg_xxyyyyz_xxzz_1, tg_xxyyyyz_xyy_1, \
                                         tg_xxyyyyz_xyyy_0, tg_xxyyyyz_xyyy_1, tg_xxyyyyz_xyyz_0, tg_xxyyyyz_xyyz_1, \
                                         tg_xxyyyyz_xyz_1, tg_xxyyyyz_xyzz_0, tg_xxyyyyz_xyzz_1, tg_xxyyyyz_xzz_1, \
                                         tg_xxyyyyz_xzzz_0, tg_xxyyyyz_xzzz_1, tg_xxyyyyz_yyy_1, tg_xxyyyyz_yyyy_0, \
                                         tg_xxyyyyz_yyyy_1, tg_xxyyyyz_yyyz_0, tg_xxyyyyz_yyyz_1, tg_xxyyyyz_yyz_1, \
                                         tg_xxyyyyz_yyzz_0, tg_xxyyyyz_yyzz_1, tg_xxyyyyz_yzz_1, tg_xxyyyyz_yzzz_0, \
                                         tg_xxyyyyz_yzzz_1, tg_xxyyyyz_zzz_1, tg_xxyyyyz_zzzz_0, tg_xxyyyyz_zzzz_1, \
                                         tg_xxyyyzz_xxx_1, tg_xxyyyzz_xxxx_0, tg_xxyyyzz_xxxx_1, tg_xxyyyzz_xxxy_0, \
                                         tg_xxyyyzz_xxxy_1, tg_xxyyyzz_xxxz_0, tg_xxyyyzz_xxxz_1, tg_xxyyyzz_xxy_1, \
                                         tg_xxyyyzz_xxyy_0, tg_xxyyyzz_xxyy_1, tg_xxyyyzz_xxyz_0, tg_xxyyyzz_xxyz_1, \
                                         tg_xxyyyzz_xxz_1, tg_xxyyyzz_xxzz_0, tg_xxyyyzz_xxzz_1, tg_xxyyyzz_xyy_1, \
                                         tg_xxyyyzz_xyyy_0, tg_xxyyyzz_xyyy_1, tg_xxyyyzz_xyyz_0, tg_xxyyyzz_xyyz_1, \
                                         tg_xxyyyzz_xyz_1, tg_xxyyyzz_xyzz_0, tg_xxyyyzz_xyzz_1, tg_xxyyyzz_xzz_1, \
                                         tg_xxyyyzz_xzzz_0, tg_xxyyyzz_xzzz_1, tg_xxyyyzz_yyy_1, tg_xxyyyzz_yyyy_0, \
                                         tg_xxyyyzz_yyyy_1, tg_xxyyyzz_yyyz_0, tg_xxyyyzz_yyyz_1, tg_xxyyyzz_yyz_1, \
                                         tg_xxyyyzz_yyzz_0, tg_xxyyyzz_yyzz_1, tg_xxyyyzz_yzz_1, tg_xxyyyzz_yzzz_0, \
                                         tg_xxyyyzz_yzzz_1, tg_xxyyyzz_zzz_1, tg_xxyyyzz_zzzz_0, tg_xxyyyzz_zzzz_1, \
                                         tg_xxyyzzz_xxx_1, tg_xxyyzzz_xxxx_0, tg_xxyyzzz_xxxx_1, tg_xxyyzzz_xxxy_0, \
                                         tg_xxyyzzz_xxxy_1, tg_xxyyzzz_xxxz_0, tg_xxyyzzz_xxxz_1, tg_xxyyzzz_xxy_1, \
                                         tg_xxyyzzz_xxyy_0, tg_xxyyzzz_xxyy_1, tg_xxyyzzz_xxyz_0, tg_xxyyzzz_xxyz_1, \
                                         tg_xxyyzzz_xxz_1, tg_xxyyzzz_xxzz_0, tg_xxyyzzz_xxzz_1, tg_xxyyzzz_xyy_1, \
                                         tg_xxyyzzz_xyyy_0, tg_xxyyzzz_xyyy_1, tg_xxyyzzz_xyyz_0, tg_xxyyzzz_xyyz_1, \
                                         tg_xxyyzzz_xyz_1, tg_xxyyzzz_xyzz_0, tg_xxyyzzz_xyzz_1, tg_xxyyzzz_xzz_1, \
                                         tg_xxyyzzz_xzzz_0, tg_xxyyzzz_xzzz_1, tg_xxyyzzz_yyy_1, tg_xxyyzzz_yyyy_0, \
                                         tg_xxyyzzz_yyyy_1, tg_xxyyzzz_yyyz_0, tg_xxyyzzz_yyyz_1, tg_xxyyzzz_yyz_1, \
                                         tg_xxyyzzz_yyzz_0, tg_xxyyzzz_yyzz_1, tg_xxyyzzz_yzz_1, tg_xxyyzzz_yzzz_0, \
                                         tg_xxyyzzz_yzzz_1, tg_xxyyzzz_zzz_1, tg_xxyyzzz_zzzz_0, tg_xxyyzzz_zzzz_1, \
                                         tg_xxyzzzz_xxx_1, tg_xxyzzzz_xxxx_0, tg_xxyzzzz_xxxx_1, tg_xxyzzzz_xxxy_0, \
                                         tg_xxyzzzz_xxxy_1, tg_xxyzzzz_xxxz_0, tg_xxyzzzz_xxxz_1, tg_xxyzzzz_xxy_1, \
                                         tg_xxyzzzz_xxyy_0, tg_xxyzzzz_xxyy_1, tg_xxyzzzz_xxyz_0, tg_xxyzzzz_xxyz_1, \
                                         tg_xxyzzzz_xxz_1, tg_xxyzzzz_xxzz_0, tg_xxyzzzz_xxzz_1, tg_xxyzzzz_xyy_1, \
                                         tg_xxyzzzz_xyz_1, tg_xxyzzzz_xzz_1, tg_xyyyyz_xxyy_0, tg_xyyyyz_xxyy_1, \
                                         tg_xyyyyz_xxyz_0, tg_xyyyyz_xxyz_1, tg_xyyyyz_xxzz_0, tg_xyyyyz_xxzz_1, \
                                         tg_xyyyyz_xyyy_0, tg_xyyyyz_xyyy_1, tg_xyyyyz_xyyz_0, tg_xyyyyz_xyyz_1, \
                                         tg_xyyyyz_xyzz_0, tg_xyyyyz_xyzz_1, tg_xyyyyz_xzzz_0, tg_xyyyyz_xzzz_1, \
                                         tg_xyyyyz_yyyy_0, tg_xyyyyz_yyyy_1, tg_xyyyyz_yyyz_0, tg_xyyyyz_yyyz_1, \
                                         tg_xyyyyz_yyzz_0, tg_xyyyyz_yyzz_1, tg_xyyyyz_yzzz_0, tg_xyyyyz_yzzz_1, \
                                         tg_xyyyyz_zzzz_0, tg_xyyyyz_zzzz_1, tg_xyyyzz_xxxx_0, tg_xyyyzz_xxxx_1, \
                                         tg_xyyyzz_xxxy_0, tg_xyyyzz_xxxy_1, tg_xyyyzz_xxxz_0, tg_xyyyzz_xxxz_1, \
                                         tg_xyyyzz_xxyy_0, tg_xyyyzz_xxyy_1, tg_xyyyzz_xxyz_0, tg_xyyyzz_xxyz_1, \
                                         tg_xyyyzz_xxzz_0, tg_xyyyzz_xxzz_1, tg_xyyyzz_xyyy_0, tg_xyyyzz_xyyy_1, \
                                         tg_xyyyzz_xyyz_0, tg_xyyyzz_xyyz_1, tg_xyyyzz_xyzz_0, tg_xyyyzz_xyzz_1, \
                                         tg_xyyyzz_xzzz_0, tg_xyyyzz_xzzz_1, tg_xyyyzz_yyyy_0, tg_xyyyzz_yyyy_1, \
                                         tg_xyyyzz_yyyz_0, tg_xyyyzz_yyyz_1, tg_xyyyzz_yyzz_0, tg_xyyyzz_yyzz_1, \
                                         tg_xyyyzz_yzzz_0, tg_xyyyzz_yzzz_1, tg_xyyyzz_zzzz_0, tg_xyyyzz_zzzz_1, \
                                         tg_xyyzzz_xxxx_0, tg_xyyzzz_xxxx_1, tg_xyyzzz_xxxy_0, tg_xyyzzz_xxxy_1, \
                                         tg_xyyzzz_xxxz_0, tg_xyyzzz_xxxz_1, tg_xyyzzz_xxyy_0, tg_xyyzzz_xxyy_1, \
                                         tg_xyyzzz_xxyz_0, tg_xyyzzz_xxyz_1, tg_xyyzzz_xxzz_0, tg_xyyzzz_xxzz_1, \
                                         tg_xyyzzz_xyyy_0, tg_xyyzzz_xyyy_1, tg_xyyzzz_xyyz_0, tg_xyyzzz_xyyz_1, \
                                         tg_xyyzzz_xyzz_0, tg_xyyzzz_xyzz_1, tg_xyyzzz_xzzz_0, tg_xyyzzz_xzzz_1, \
                                         tg_xyyzzz_yyyy_0, tg_xyyzzz_yyyy_1, tg_xyyzzz_yyyz_0, tg_xyyzzz_yyyz_1, \
                                         tg_xyyzzz_yyzz_0, tg_xyyzzz_yyzz_1, tg_xyyzzz_yzzz_0, tg_xyyzzz_yzzz_1, \
                                         tg_xyyzzz_zzzz_0, tg_xyyzzz_zzzz_1, tg_xyzzzz_xxxx_0, tg_xyzzzz_xxxx_1, \
                                         tg_xyzzzz_xxxy_0, tg_xyzzzz_xxxy_1, tg_xyzzzz_xxxz_0, tg_xyzzzz_xxxz_1, \
                                         tg_xyzzzz_xxyy_0, tg_xyzzzz_xxyy_1, tg_xyzzzz_xxyz_0, tg_xyzzzz_xxyz_1, \
                                         tg_xyzzzz_xxzz_0, tg_xyzzzz_xxzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxxyyyyz_xxyy_0[j] = pb_x * tg_xxyyyyz_xxyy_0[j] + wp_x[j] * tg_xxyyyyz_xxyy_1[j] + fl1_fx * tg_xyyyyz_xxyy_0[j] - fl1_fx * fl1_fza * tg_xyyyyz_xxyy_1[j] + fl1_fxn * tg_xxyyyyz_xyy_1[j];

                    tg_xxxyyyyz_xxyz_0[j] = pb_x * tg_xxyyyyz_xxyz_0[j] + wp_x[j] * tg_xxyyyyz_xxyz_1[j] + fl1_fx * tg_xyyyyz_xxyz_0[j] - fl1_fx * fl1_fza * tg_xyyyyz_xxyz_1[j] + fl1_fxn * tg_xxyyyyz_xyz_1[j];

                    tg_xxxyyyyz_xxzz_0[j] = pb_x * tg_xxyyyyz_xxzz_0[j] + wp_x[j] * tg_xxyyyyz_xxzz_1[j] + fl1_fx * tg_xyyyyz_xxzz_0[j] - fl1_fx * fl1_fza * tg_xyyyyz_xxzz_1[j] + fl1_fxn * tg_xxyyyyz_xzz_1[j];

                    tg_xxxyyyyz_xyyy_0[j] = pb_x * tg_xxyyyyz_xyyy_0[j] + wp_x[j] * tg_xxyyyyz_xyyy_1[j] + fl1_fx * tg_xyyyyz_xyyy_0[j] - fl1_fx * fl1_fza * tg_xyyyyz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxyyyyz_yyy_1[j];

                    tg_xxxyyyyz_xyyz_0[j] = pb_x * tg_xxyyyyz_xyyz_0[j] + wp_x[j] * tg_xxyyyyz_xyyz_1[j] + fl1_fx * tg_xyyyyz_xyyz_0[j] - fl1_fx * fl1_fza * tg_xyyyyz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxyyyyz_yyz_1[j];

                    tg_xxxyyyyz_xyzz_0[j] = pb_x * tg_xxyyyyz_xyzz_0[j] + wp_x[j] * tg_xxyyyyz_xyzz_1[j] + fl1_fx * tg_xyyyyz_xyzz_0[j] - fl1_fx * fl1_fza * tg_xyyyyz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxyyyyz_yzz_1[j];

                    tg_xxxyyyyz_xzzz_0[j] = pb_x * tg_xxyyyyz_xzzz_0[j] + wp_x[j] * tg_xxyyyyz_xzzz_1[j] + fl1_fx * tg_xyyyyz_xzzz_0[j] - fl1_fx * fl1_fza * tg_xyyyyz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxyyyyz_zzz_1[j];

                    tg_xxxyyyyz_yyyy_0[j] = pb_x * tg_xxyyyyz_yyyy_0[j] + wp_x[j] * tg_xxyyyyz_yyyy_1[j] + fl1_fx * tg_xyyyyz_yyyy_0[j] - fl1_fx * fl1_fza * tg_xyyyyz_yyyy_1[j];

                    tg_xxxyyyyz_yyyz_0[j] = pb_x * tg_xxyyyyz_yyyz_0[j] + wp_x[j] * tg_xxyyyyz_yyyz_1[j] + fl1_fx * tg_xyyyyz_yyyz_0[j] - fl1_fx * fl1_fza * tg_xyyyyz_yyyz_1[j];

                    tg_xxxyyyyz_yyzz_0[j] = pb_x * tg_xxyyyyz_yyzz_0[j] + wp_x[j] * tg_xxyyyyz_yyzz_1[j] + fl1_fx * tg_xyyyyz_yyzz_0[j] - fl1_fx * fl1_fza * tg_xyyyyz_yyzz_1[j];

                    tg_xxxyyyyz_yzzz_0[j] = pb_x * tg_xxyyyyz_yzzz_0[j] + wp_x[j] * tg_xxyyyyz_yzzz_1[j] + fl1_fx * tg_xyyyyz_yzzz_0[j] - fl1_fx * fl1_fza * tg_xyyyyz_yzzz_1[j];

                    tg_xxxyyyyz_zzzz_0[j] = pb_x * tg_xxyyyyz_zzzz_0[j] + wp_x[j] * tg_xxyyyyz_zzzz_1[j] + fl1_fx * tg_xyyyyz_zzzz_0[j] - fl1_fx * fl1_fza * tg_xyyyyz_zzzz_1[j];

                    tg_xxxyyyzz_xxxx_0[j] = pb_x * tg_xxyyyzz_xxxx_0[j] + wp_x[j] * tg_xxyyyzz_xxxx_1[j] + fl1_fx * tg_xyyyzz_xxxx_0[j] - fl1_fx * fl1_fza * tg_xyyyzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxyyyzz_xxx_1[j];

                    tg_xxxyyyzz_xxxy_0[j] = pb_x * tg_xxyyyzz_xxxy_0[j] + wp_x[j] * tg_xxyyyzz_xxxy_1[j] + fl1_fx * tg_xyyyzz_xxxy_0[j] - fl1_fx * fl1_fza * tg_xyyyzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxyyyzz_xxy_1[j];

                    tg_xxxyyyzz_xxxz_0[j] = pb_x * tg_xxyyyzz_xxxz_0[j] + wp_x[j] * tg_xxyyyzz_xxxz_1[j] + fl1_fx * tg_xyyyzz_xxxz_0[j] - fl1_fx * fl1_fza * tg_xyyyzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxyyyzz_xxz_1[j];

                    tg_xxxyyyzz_xxyy_0[j] = pb_x * tg_xxyyyzz_xxyy_0[j] + wp_x[j] * tg_xxyyyzz_xxyy_1[j] + fl1_fx * tg_xyyyzz_xxyy_0[j] - fl1_fx * fl1_fza * tg_xyyyzz_xxyy_1[j] + fl1_fxn * tg_xxyyyzz_xyy_1[j];

                    tg_xxxyyyzz_xxyz_0[j] = pb_x * tg_xxyyyzz_xxyz_0[j] + wp_x[j] * tg_xxyyyzz_xxyz_1[j] + fl1_fx * tg_xyyyzz_xxyz_0[j] - fl1_fx * fl1_fza * tg_xyyyzz_xxyz_1[j] + fl1_fxn * tg_xxyyyzz_xyz_1[j];

                    tg_xxxyyyzz_xxzz_0[j] = pb_x * tg_xxyyyzz_xxzz_0[j] + wp_x[j] * tg_xxyyyzz_xxzz_1[j] + fl1_fx * tg_xyyyzz_xxzz_0[j] - fl1_fx * fl1_fza * tg_xyyyzz_xxzz_1[j] + fl1_fxn * tg_xxyyyzz_xzz_1[j];

                    tg_xxxyyyzz_xyyy_0[j] = pb_x * tg_xxyyyzz_xyyy_0[j] + wp_x[j] * tg_xxyyyzz_xyyy_1[j] + fl1_fx * tg_xyyyzz_xyyy_0[j] - fl1_fx * fl1_fza * tg_xyyyzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxyyyzz_yyy_1[j];

                    tg_xxxyyyzz_xyyz_0[j] = pb_x * tg_xxyyyzz_xyyz_0[j] + wp_x[j] * tg_xxyyyzz_xyyz_1[j] + fl1_fx * tg_xyyyzz_xyyz_0[j] - fl1_fx * fl1_fza * tg_xyyyzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxyyyzz_yyz_1[j];

                    tg_xxxyyyzz_xyzz_0[j] = pb_x * tg_xxyyyzz_xyzz_0[j] + wp_x[j] * tg_xxyyyzz_xyzz_1[j] + fl1_fx * tg_xyyyzz_xyzz_0[j] - fl1_fx * fl1_fza * tg_xyyyzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxyyyzz_yzz_1[j];

                    tg_xxxyyyzz_xzzz_0[j] = pb_x * tg_xxyyyzz_xzzz_0[j] + wp_x[j] * tg_xxyyyzz_xzzz_1[j] + fl1_fx * tg_xyyyzz_xzzz_0[j] - fl1_fx * fl1_fza * tg_xyyyzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxyyyzz_zzz_1[j];

                    tg_xxxyyyzz_yyyy_0[j] = pb_x * tg_xxyyyzz_yyyy_0[j] + wp_x[j] * tg_xxyyyzz_yyyy_1[j] + fl1_fx * tg_xyyyzz_yyyy_0[j] - fl1_fx * fl1_fza * tg_xyyyzz_yyyy_1[j];

                    tg_xxxyyyzz_yyyz_0[j] = pb_x * tg_xxyyyzz_yyyz_0[j] + wp_x[j] * tg_xxyyyzz_yyyz_1[j] + fl1_fx * tg_xyyyzz_yyyz_0[j] - fl1_fx * fl1_fza * tg_xyyyzz_yyyz_1[j];

                    tg_xxxyyyzz_yyzz_0[j] = pb_x * tg_xxyyyzz_yyzz_0[j] + wp_x[j] * tg_xxyyyzz_yyzz_1[j] + fl1_fx * tg_xyyyzz_yyzz_0[j] - fl1_fx * fl1_fza * tg_xyyyzz_yyzz_1[j];

                    tg_xxxyyyzz_yzzz_0[j] = pb_x * tg_xxyyyzz_yzzz_0[j] + wp_x[j] * tg_xxyyyzz_yzzz_1[j] + fl1_fx * tg_xyyyzz_yzzz_0[j] - fl1_fx * fl1_fza * tg_xyyyzz_yzzz_1[j];

                    tg_xxxyyyzz_zzzz_0[j] = pb_x * tg_xxyyyzz_zzzz_0[j] + wp_x[j] * tg_xxyyyzz_zzzz_1[j] + fl1_fx * tg_xyyyzz_zzzz_0[j] - fl1_fx * fl1_fza * tg_xyyyzz_zzzz_1[j];

                    tg_xxxyyzzz_xxxx_0[j] = pb_x * tg_xxyyzzz_xxxx_0[j] + wp_x[j] * tg_xxyyzzz_xxxx_1[j] + fl1_fx * tg_xyyzzz_xxxx_0[j] - fl1_fx * fl1_fza * tg_xyyzzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxyyzzz_xxx_1[j];

                    tg_xxxyyzzz_xxxy_0[j] = pb_x * tg_xxyyzzz_xxxy_0[j] + wp_x[j] * tg_xxyyzzz_xxxy_1[j] + fl1_fx * tg_xyyzzz_xxxy_0[j] - fl1_fx * fl1_fza * tg_xyyzzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxyyzzz_xxy_1[j];

                    tg_xxxyyzzz_xxxz_0[j] = pb_x * tg_xxyyzzz_xxxz_0[j] + wp_x[j] * tg_xxyyzzz_xxxz_1[j] + fl1_fx * tg_xyyzzz_xxxz_0[j] - fl1_fx * fl1_fza * tg_xyyzzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxyyzzz_xxz_1[j];

                    tg_xxxyyzzz_xxyy_0[j] = pb_x * tg_xxyyzzz_xxyy_0[j] + wp_x[j] * tg_xxyyzzz_xxyy_1[j] + fl1_fx * tg_xyyzzz_xxyy_0[j] - fl1_fx * fl1_fza * tg_xyyzzz_xxyy_1[j] + fl1_fxn * tg_xxyyzzz_xyy_1[j];

                    tg_xxxyyzzz_xxyz_0[j] = pb_x * tg_xxyyzzz_xxyz_0[j] + wp_x[j] * tg_xxyyzzz_xxyz_1[j] + fl1_fx * tg_xyyzzz_xxyz_0[j] - fl1_fx * fl1_fza * tg_xyyzzz_xxyz_1[j] + fl1_fxn * tg_xxyyzzz_xyz_1[j];

                    tg_xxxyyzzz_xxzz_0[j] = pb_x * tg_xxyyzzz_xxzz_0[j] + wp_x[j] * tg_xxyyzzz_xxzz_1[j] + fl1_fx * tg_xyyzzz_xxzz_0[j] - fl1_fx * fl1_fza * tg_xyyzzz_xxzz_1[j] + fl1_fxn * tg_xxyyzzz_xzz_1[j];

                    tg_xxxyyzzz_xyyy_0[j] = pb_x * tg_xxyyzzz_xyyy_0[j] + wp_x[j] * tg_xxyyzzz_xyyy_1[j] + fl1_fx * tg_xyyzzz_xyyy_0[j] - fl1_fx * fl1_fza * tg_xyyzzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxyyzzz_yyy_1[j];

                    tg_xxxyyzzz_xyyz_0[j] = pb_x * tg_xxyyzzz_xyyz_0[j] + wp_x[j] * tg_xxyyzzz_xyyz_1[j] + fl1_fx * tg_xyyzzz_xyyz_0[j] - fl1_fx * fl1_fza * tg_xyyzzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxyyzzz_yyz_1[j];

                    tg_xxxyyzzz_xyzz_0[j] = pb_x * tg_xxyyzzz_xyzz_0[j] + wp_x[j] * tg_xxyyzzz_xyzz_1[j] + fl1_fx * tg_xyyzzz_xyzz_0[j] - fl1_fx * fl1_fza * tg_xyyzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxyyzzz_yzz_1[j];

                    tg_xxxyyzzz_xzzz_0[j] = pb_x * tg_xxyyzzz_xzzz_0[j] + wp_x[j] * tg_xxyyzzz_xzzz_1[j] + fl1_fx * tg_xyyzzz_xzzz_0[j] - fl1_fx * fl1_fza * tg_xyyzzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxyyzzz_zzz_1[j];

                    tg_xxxyyzzz_yyyy_0[j] = pb_x * tg_xxyyzzz_yyyy_0[j] + wp_x[j] * tg_xxyyzzz_yyyy_1[j] + fl1_fx * tg_xyyzzz_yyyy_0[j] - fl1_fx * fl1_fza * tg_xyyzzz_yyyy_1[j];

                    tg_xxxyyzzz_yyyz_0[j] = pb_x * tg_xxyyzzz_yyyz_0[j] + wp_x[j] * tg_xxyyzzz_yyyz_1[j] + fl1_fx * tg_xyyzzz_yyyz_0[j] - fl1_fx * fl1_fza * tg_xyyzzz_yyyz_1[j];

                    tg_xxxyyzzz_yyzz_0[j] = pb_x * tg_xxyyzzz_yyzz_0[j] + wp_x[j] * tg_xxyyzzz_yyzz_1[j] + fl1_fx * tg_xyyzzz_yyzz_0[j] - fl1_fx * fl1_fza * tg_xyyzzz_yyzz_1[j];

                    tg_xxxyyzzz_yzzz_0[j] = pb_x * tg_xxyyzzz_yzzz_0[j] + wp_x[j] * tg_xxyyzzz_yzzz_1[j] + fl1_fx * tg_xyyzzz_yzzz_0[j] - fl1_fx * fl1_fza * tg_xyyzzz_yzzz_1[j];

                    tg_xxxyyzzz_zzzz_0[j] = pb_x * tg_xxyyzzz_zzzz_0[j] + wp_x[j] * tg_xxyyzzz_zzzz_1[j] + fl1_fx * tg_xyyzzz_zzzz_0[j] - fl1_fx * fl1_fza * tg_xyyzzz_zzzz_1[j];

                    tg_xxxyzzzz_xxxx_0[j] = pb_x * tg_xxyzzzz_xxxx_0[j] + wp_x[j] * tg_xxyzzzz_xxxx_1[j] + fl1_fx * tg_xyzzzz_xxxx_0[j] - fl1_fx * fl1_fza * tg_xyzzzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxyzzzz_xxx_1[j];

                    tg_xxxyzzzz_xxxy_0[j] = pb_x * tg_xxyzzzz_xxxy_0[j] + wp_x[j] * tg_xxyzzzz_xxxy_1[j] + fl1_fx * tg_xyzzzz_xxxy_0[j] - fl1_fx * fl1_fza * tg_xyzzzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxyzzzz_xxy_1[j];

                    tg_xxxyzzzz_xxxz_0[j] = pb_x * tg_xxyzzzz_xxxz_0[j] + wp_x[j] * tg_xxyzzzz_xxxz_1[j] + fl1_fx * tg_xyzzzz_xxxz_0[j] - fl1_fx * fl1_fza * tg_xyzzzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxyzzzz_xxz_1[j];

                    tg_xxxyzzzz_xxyy_0[j] = pb_x * tg_xxyzzzz_xxyy_0[j] + wp_x[j] * tg_xxyzzzz_xxyy_1[j] + fl1_fx * tg_xyzzzz_xxyy_0[j] - fl1_fx * fl1_fza * tg_xyzzzz_xxyy_1[j] + fl1_fxn * tg_xxyzzzz_xyy_1[j];

                    tg_xxxyzzzz_xxyz_0[j] = pb_x * tg_xxyzzzz_xxyz_0[j] + wp_x[j] * tg_xxyzzzz_xxyz_1[j] + fl1_fx * tg_xyzzzz_xxyz_0[j] - fl1_fx * fl1_fza * tg_xyzzzz_xxyz_1[j] + fl1_fxn * tg_xxyzzzz_xyz_1[j];

                    tg_xxxyzzzz_xxzz_0[j] = pb_x * tg_xxyzzzz_xxzz_0[j] + wp_x[j] * tg_xxyzzzz_xxzz_1[j] + fl1_fx * tg_xyzzzz_xxzz_0[j] - fl1_fx * fl1_fza * tg_xyzzzz_xxzz_1[j] + fl1_fxn * tg_xxyzzzz_xzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSLSG_291_339(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (291,339)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {8, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {8, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_7_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_6_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_6_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_7_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_xxyzzzz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 291); 

                auto tg_xxyzzzz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 292); 

                auto tg_xxyzzzz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 293); 

                auto tg_xxyzzzz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 294); 

                auto tg_xxyzzzz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 295); 

                auto tg_xxyzzzz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 296); 

                auto tg_xxyzzzz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 297); 

                auto tg_xxyzzzz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 298); 

                auto tg_xxyzzzz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 299); 

                auto tg_xxzzzzz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 300); 

                auto tg_xxzzzzz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 301); 

                auto tg_xxzzzzz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 302); 

                auto tg_xxzzzzz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 303); 

                auto tg_xxzzzzz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 304); 

                auto tg_xxzzzzz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 305); 

                auto tg_xxzzzzz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 306); 

                auto tg_xxzzzzz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 307); 

                auto tg_xxzzzzz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 308); 

                auto tg_xxzzzzz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 309); 

                auto tg_xxzzzzz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 310); 

                auto tg_xxzzzzz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 311); 

                auto tg_xxzzzzz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 312); 

                auto tg_xxzzzzz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 313); 

                auto tg_xxzzzzz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 314); 

                auto tg_xyyyyyy_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 315); 

                auto tg_xyyyyyy_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 316); 

                auto tg_xyyyyyy_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 317); 

                auto tg_xyyyyyy_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 318); 

                auto tg_xyyyyyy_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 319); 

                auto tg_xyyyyyy_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 320); 

                auto tg_xyyyyyy_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 321); 

                auto tg_xyyyyyy_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 322); 

                auto tg_xyyyyyy_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 323); 

                auto tg_xyyyyyy_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 324); 

                auto tg_xyyyyyy_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 325); 

                auto tg_xyyyyyy_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 326); 

                auto tg_xyyyyyy_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 327); 

                auto tg_xyyyyyy_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 328); 

                auto tg_xyyyyyy_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 329); 

                auto tg_xyyyyyz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 330); 

                auto tg_xyyyyyz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 331); 

                auto tg_xyyyyyz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 332); 

                auto tg_xyyyyyz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 333); 

                auto tg_xyyyyyz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 334); 

                auto tg_xyyyyyz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 335); 

                auto tg_xyyyyyz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 336); 

                auto tg_xyyyyyz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 337); 

                auto tg_xyyyyyz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 338); 

                auto tg_xxyzzzz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 291); 

                auto tg_xxyzzzz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 292); 

                auto tg_xxyzzzz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 293); 

                auto tg_xxyzzzz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 294); 

                auto tg_xxyzzzz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 295); 

                auto tg_xxyzzzz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 296); 

                auto tg_xxyzzzz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 297); 

                auto tg_xxyzzzz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 298); 

                auto tg_xxyzzzz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 299); 

                auto tg_xxzzzzz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 300); 

                auto tg_xxzzzzz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 301); 

                auto tg_xxzzzzz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 302); 

                auto tg_xxzzzzz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 303); 

                auto tg_xxzzzzz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 304); 

                auto tg_xxzzzzz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 305); 

                auto tg_xxzzzzz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 306); 

                auto tg_xxzzzzz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 307); 

                auto tg_xxzzzzz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 308); 

                auto tg_xxzzzzz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 309); 

                auto tg_xxzzzzz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 310); 

                auto tg_xxzzzzz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 311); 

                auto tg_xxzzzzz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 312); 

                auto tg_xxzzzzz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 313); 

                auto tg_xxzzzzz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 314); 

                auto tg_xyyyyyy_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 315); 

                auto tg_xyyyyyy_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 316); 

                auto tg_xyyyyyy_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 317); 

                auto tg_xyyyyyy_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 318); 

                auto tg_xyyyyyy_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 319); 

                auto tg_xyyyyyy_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 320); 

                auto tg_xyyyyyy_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 321); 

                auto tg_xyyyyyy_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 322); 

                auto tg_xyyyyyy_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 323); 

                auto tg_xyyyyyy_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 324); 

                auto tg_xyyyyyy_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 325); 

                auto tg_xyyyyyy_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 326); 

                auto tg_xyyyyyy_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 327); 

                auto tg_xyyyyyy_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 328); 

                auto tg_xyyyyyy_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 329); 

                auto tg_xyyyyyz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 330); 

                auto tg_xyyyyyz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 331); 

                auto tg_xyyyyyz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 332); 

                auto tg_xyyyyyz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 333); 

                auto tg_xyyyyyz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 334); 

                auto tg_xyyyyyz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 335); 

                auto tg_xyyyyyz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 336); 

                auto tg_xyyyyyz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 337); 

                auto tg_xyyyyyz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 338); 

                auto tg_xyzzzz_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 291); 

                auto tg_xyzzzz_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 292); 

                auto tg_xyzzzz_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 293); 

                auto tg_xyzzzz_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 294); 

                auto tg_xyzzzz_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 295); 

                auto tg_xyzzzz_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 296); 

                auto tg_xyzzzz_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 297); 

                auto tg_xyzzzz_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 298); 

                auto tg_xyzzzz_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 299); 

                auto tg_xzzzzz_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 300); 

                auto tg_xzzzzz_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 301); 

                auto tg_xzzzzz_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 302); 

                auto tg_xzzzzz_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 303); 

                auto tg_xzzzzz_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 304); 

                auto tg_xzzzzz_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 305); 

                auto tg_xzzzzz_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 306); 

                auto tg_xzzzzz_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 307); 

                auto tg_xzzzzz_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 308); 

                auto tg_xzzzzz_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 309); 

                auto tg_xzzzzz_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 310); 

                auto tg_xzzzzz_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 311); 

                auto tg_xzzzzz_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 312); 

                auto tg_xzzzzz_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 313); 

                auto tg_xzzzzz_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 314); 

                auto tg_yyyyyy_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 315); 

                auto tg_yyyyyy_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 316); 

                auto tg_yyyyyy_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 317); 

                auto tg_yyyyyy_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 318); 

                auto tg_yyyyyy_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 319); 

                auto tg_yyyyyy_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 320); 

                auto tg_yyyyyy_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 321); 

                auto tg_yyyyyy_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 322); 

                auto tg_yyyyyy_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 323); 

                auto tg_yyyyyy_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 324); 

                auto tg_yyyyyy_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 325); 

                auto tg_yyyyyy_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 326); 

                auto tg_yyyyyy_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 327); 

                auto tg_yyyyyy_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 328); 

                auto tg_yyyyyy_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 329); 

                auto tg_yyyyyz_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 330); 

                auto tg_yyyyyz_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 331); 

                auto tg_yyyyyz_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 332); 

                auto tg_yyyyyz_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 333); 

                auto tg_yyyyyz_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 334); 

                auto tg_yyyyyz_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 335); 

                auto tg_yyyyyz_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 336); 

                auto tg_yyyyyz_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 337); 

                auto tg_yyyyyz_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 338); 

                auto tg_xyzzzz_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 291); 

                auto tg_xyzzzz_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 292); 

                auto tg_xyzzzz_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 293); 

                auto tg_xyzzzz_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 294); 

                auto tg_xyzzzz_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 295); 

                auto tg_xyzzzz_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 296); 

                auto tg_xyzzzz_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 297); 

                auto tg_xyzzzz_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 298); 

                auto tg_xyzzzz_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 299); 

                auto tg_xzzzzz_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 300); 

                auto tg_xzzzzz_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 301); 

                auto tg_xzzzzz_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 302); 

                auto tg_xzzzzz_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 303); 

                auto tg_xzzzzz_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 304); 

                auto tg_xzzzzz_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 305); 

                auto tg_xzzzzz_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 306); 

                auto tg_xzzzzz_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 307); 

                auto tg_xzzzzz_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 308); 

                auto tg_xzzzzz_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 309); 

                auto tg_xzzzzz_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 310); 

                auto tg_xzzzzz_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 311); 

                auto tg_xzzzzz_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 312); 

                auto tg_xzzzzz_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 313); 

                auto tg_xzzzzz_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 314); 

                auto tg_yyyyyy_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 315); 

                auto tg_yyyyyy_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 316); 

                auto tg_yyyyyy_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 317); 

                auto tg_yyyyyy_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 318); 

                auto tg_yyyyyy_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 319); 

                auto tg_yyyyyy_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 320); 

                auto tg_yyyyyy_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 321); 

                auto tg_yyyyyy_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 322); 

                auto tg_yyyyyy_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 323); 

                auto tg_yyyyyy_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 324); 

                auto tg_yyyyyy_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 325); 

                auto tg_yyyyyy_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 326); 

                auto tg_yyyyyy_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 327); 

                auto tg_yyyyyy_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 328); 

                auto tg_yyyyyy_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 329); 

                auto tg_yyyyyz_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 330); 

                auto tg_yyyyyz_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 331); 

                auto tg_yyyyyz_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 332); 

                auto tg_yyyyyz_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 333); 

                auto tg_yyyyyz_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 334); 

                auto tg_yyyyyz_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 335); 

                auto tg_yyyyyz_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 336); 

                auto tg_yyyyyz_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 337); 

                auto tg_yyyyyz_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 338); 

                auto tg_xxyzzzz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 196); 

                auto tg_xxyzzzz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 197); 

                auto tg_xxyzzzz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 198); 

                auto tg_xxyzzzz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 199); 

                auto tg_xxzzzzz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 200); 

                auto tg_xxzzzzz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 201); 

                auto tg_xxzzzzz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 202); 

                auto tg_xxzzzzz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 203); 

                auto tg_xxzzzzz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 204); 

                auto tg_xxzzzzz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 205); 

                auto tg_xxzzzzz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 206); 

                auto tg_xxzzzzz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 207); 

                auto tg_xxzzzzz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 208); 

                auto tg_xxzzzzz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 209); 

                auto tg_xyyyyyy_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 210); 

                auto tg_xyyyyyy_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 211); 

                auto tg_xyyyyyy_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 212); 

                auto tg_xyyyyyy_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 213); 

                auto tg_xyyyyyy_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 214); 

                auto tg_xyyyyyy_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 215); 

                auto tg_xyyyyyy_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 216); 

                auto tg_xyyyyyy_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 217); 

                auto tg_xyyyyyy_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 218); 

                auto tg_xyyyyyy_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 219); 

                auto tg_xyyyyyz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 220); 

                auto tg_xyyyyyz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 221); 

                auto tg_xyyyyyz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 222); 

                auto tg_xyyyyyz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 223); 

                auto tg_xyyyyyz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 224); 

                auto tg_xyyyyyz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 225); 

                auto tg_xyyyyyz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 226); 

                auto tg_xyyyyyz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 227); 

                auto tg_xyyyyyz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 228); 

                // set up pointers to integrals

                auto tg_xxxyzzzz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 291); 

                auto tg_xxxyzzzz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 292); 

                auto tg_xxxyzzzz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 293); 

                auto tg_xxxyzzzz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 294); 

                auto tg_xxxyzzzz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 295); 

                auto tg_xxxyzzzz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 296); 

                auto tg_xxxyzzzz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 297); 

                auto tg_xxxyzzzz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 298); 

                auto tg_xxxyzzzz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 299); 

                auto tg_xxxzzzzz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 300); 

                auto tg_xxxzzzzz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 301); 

                auto tg_xxxzzzzz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 302); 

                auto tg_xxxzzzzz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 303); 

                auto tg_xxxzzzzz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 304); 

                auto tg_xxxzzzzz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 305); 

                auto tg_xxxzzzzz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 306); 

                auto tg_xxxzzzzz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 307); 

                auto tg_xxxzzzzz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 308); 

                auto tg_xxxzzzzz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 309); 

                auto tg_xxxzzzzz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 310); 

                auto tg_xxxzzzzz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 311); 

                auto tg_xxxzzzzz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 312); 

                auto tg_xxxzzzzz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 313); 

                auto tg_xxxzzzzz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 314); 

                auto tg_xxyyyyyy_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 315); 

                auto tg_xxyyyyyy_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 316); 

                auto tg_xxyyyyyy_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 317); 

                auto tg_xxyyyyyy_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 318); 

                auto tg_xxyyyyyy_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 319); 

                auto tg_xxyyyyyy_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 320); 

                auto tg_xxyyyyyy_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 321); 

                auto tg_xxyyyyyy_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 322); 

                auto tg_xxyyyyyy_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 323); 

                auto tg_xxyyyyyy_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 324); 

                auto tg_xxyyyyyy_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 325); 

                auto tg_xxyyyyyy_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 326); 

                auto tg_xxyyyyyy_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 327); 

                auto tg_xxyyyyyy_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 328); 

                auto tg_xxyyyyyy_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 329); 

                auto tg_xxyyyyyz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 330); 

                auto tg_xxyyyyyz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 331); 

                auto tg_xxyyyyyz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 332); 

                auto tg_xxyyyyyz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 333); 

                auto tg_xxyyyyyz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 334); 

                auto tg_xxyyyyyz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 335); 

                auto tg_xxyyyyyz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 336); 

                auto tg_xxyyyyyz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 337); 

                auto tg_xxyyyyyz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 338); 

                // Batch of Integrals (291,339)

                #pragma omp simd aligned(fxn, fza, tg_xxxyzzzz_xyyy_0, tg_xxxyzzzz_xyyz_0, tg_xxxyzzzz_xyzz_0, \
                                         tg_xxxyzzzz_xzzz_0, tg_xxxyzzzz_yyyy_0, tg_xxxyzzzz_yyyz_0, tg_xxxyzzzz_yyzz_0, \
                                         tg_xxxyzzzz_yzzz_0, tg_xxxyzzzz_zzzz_0, tg_xxxzzzzz_xxxx_0, tg_xxxzzzzz_xxxy_0, \
                                         tg_xxxzzzzz_xxxz_0, tg_xxxzzzzz_xxyy_0, tg_xxxzzzzz_xxyz_0, tg_xxxzzzzz_xxzz_0, \
                                         tg_xxxzzzzz_xyyy_0, tg_xxxzzzzz_xyyz_0, tg_xxxzzzzz_xyzz_0, tg_xxxzzzzz_xzzz_0, \
                                         tg_xxxzzzzz_yyyy_0, tg_xxxzzzzz_yyyz_0, tg_xxxzzzzz_yyzz_0, tg_xxxzzzzz_yzzz_0, \
                                         tg_xxxzzzzz_zzzz_0, tg_xxyyyyyy_xxxx_0, tg_xxyyyyyy_xxxy_0, tg_xxyyyyyy_xxxz_0, \
                                         tg_xxyyyyyy_xxyy_0, tg_xxyyyyyy_xxyz_0, tg_xxyyyyyy_xxzz_0, tg_xxyyyyyy_xyyy_0, \
                                         tg_xxyyyyyy_xyyz_0, tg_xxyyyyyy_xyzz_0, tg_xxyyyyyy_xzzz_0, tg_xxyyyyyy_yyyy_0, \
                                         tg_xxyyyyyy_yyyz_0, tg_xxyyyyyy_yyzz_0, tg_xxyyyyyy_yzzz_0, tg_xxyyyyyy_zzzz_0, \
                                         tg_xxyyyyyz_xxxx_0, tg_xxyyyyyz_xxxy_0, tg_xxyyyyyz_xxxz_0, tg_xxyyyyyz_xxyy_0, \
                                         tg_xxyyyyyz_xxyz_0, tg_xxyyyyyz_xxzz_0, tg_xxyyyyyz_xyyy_0, tg_xxyyyyyz_xyyz_0, \
                                         tg_xxyyyyyz_xyzz_0, tg_xxyzzzz_xyyy_0, tg_xxyzzzz_xyyy_1, tg_xxyzzzz_xyyz_0, \
                                         tg_xxyzzzz_xyyz_1, tg_xxyzzzz_xyzz_0, tg_xxyzzzz_xyzz_1, tg_xxyzzzz_xzzz_0, \
                                         tg_xxyzzzz_xzzz_1, tg_xxyzzzz_yyy_1, tg_xxyzzzz_yyyy_0, tg_xxyzzzz_yyyy_1, \
                                         tg_xxyzzzz_yyyz_0, tg_xxyzzzz_yyyz_1, tg_xxyzzzz_yyz_1, tg_xxyzzzz_yyzz_0, \
                                         tg_xxyzzzz_yyzz_1, tg_xxyzzzz_yzz_1, tg_xxyzzzz_yzzz_0, tg_xxyzzzz_yzzz_1, \
                                         tg_xxyzzzz_zzz_1, tg_xxyzzzz_zzzz_0, tg_xxyzzzz_zzzz_1, tg_xxzzzzz_xxx_1, \
                                         tg_xxzzzzz_xxxx_0, tg_xxzzzzz_xxxx_1, tg_xxzzzzz_xxxy_0, tg_xxzzzzz_xxxy_1, \
                                         tg_xxzzzzz_xxxz_0, tg_xxzzzzz_xxxz_1, tg_xxzzzzz_xxy_1, tg_xxzzzzz_xxyy_0, \
                                         tg_xxzzzzz_xxyy_1, tg_xxzzzzz_xxyz_0, tg_xxzzzzz_xxyz_1, tg_xxzzzzz_xxz_1, \
                                         tg_xxzzzzz_xxzz_0, tg_xxzzzzz_xxzz_1, tg_xxzzzzz_xyy_1, tg_xxzzzzz_xyyy_0, \
                                         tg_xxzzzzz_xyyy_1, tg_xxzzzzz_xyyz_0, tg_xxzzzzz_xyyz_1, tg_xxzzzzz_xyz_1, \
                                         tg_xxzzzzz_xyzz_0, tg_xxzzzzz_xyzz_1, tg_xxzzzzz_xzz_1, tg_xxzzzzz_xzzz_0, \
                                         tg_xxzzzzz_xzzz_1, tg_xxzzzzz_yyy_1, tg_xxzzzzz_yyyy_0, tg_xxzzzzz_yyyy_1, \
                                         tg_xxzzzzz_yyyz_0, tg_xxzzzzz_yyyz_1, tg_xxzzzzz_yyz_1, tg_xxzzzzz_yyzz_0, \
                                         tg_xxzzzzz_yyzz_1, tg_xxzzzzz_yzz_1, tg_xxzzzzz_yzzz_0, tg_xxzzzzz_yzzz_1, \
                                         tg_xxzzzzz_zzz_1, tg_xxzzzzz_zzzz_0, tg_xxzzzzz_zzzz_1, tg_xyyyyyy_xxx_1, \
                                         tg_xyyyyyy_xxxx_0, tg_xyyyyyy_xxxx_1, tg_xyyyyyy_xxxy_0, tg_xyyyyyy_xxxy_1, \
                                         tg_xyyyyyy_xxxz_0, tg_xyyyyyy_xxxz_1, tg_xyyyyyy_xxy_1, tg_xyyyyyy_xxyy_0, \
                                         tg_xyyyyyy_xxyy_1, tg_xyyyyyy_xxyz_0, tg_xyyyyyy_xxyz_1, tg_xyyyyyy_xxz_1, \
                                         tg_xyyyyyy_xxzz_0, tg_xyyyyyy_xxzz_1, tg_xyyyyyy_xyy_1, tg_xyyyyyy_xyyy_0, \
                                         tg_xyyyyyy_xyyy_1, tg_xyyyyyy_xyyz_0, tg_xyyyyyy_xyyz_1, tg_xyyyyyy_xyz_1, \
                                         tg_xyyyyyy_xyzz_0, tg_xyyyyyy_xyzz_1, tg_xyyyyyy_xzz_1, tg_xyyyyyy_xzzz_0, \
                                         tg_xyyyyyy_xzzz_1, tg_xyyyyyy_yyy_1, tg_xyyyyyy_yyyy_0, tg_xyyyyyy_yyyy_1, \
                                         tg_xyyyyyy_yyyz_0, tg_xyyyyyy_yyyz_1, tg_xyyyyyy_yyz_1, tg_xyyyyyy_yyzz_0, \
                                         tg_xyyyyyy_yyzz_1, tg_xyyyyyy_yzz_1, tg_xyyyyyy_yzzz_0, tg_xyyyyyy_yzzz_1, \
                                         tg_xyyyyyy_zzz_1, tg_xyyyyyy_zzzz_0, tg_xyyyyyy_zzzz_1, tg_xyyyyyz_xxx_1, \
                                         tg_xyyyyyz_xxxx_0, tg_xyyyyyz_xxxx_1, tg_xyyyyyz_xxxy_0, tg_xyyyyyz_xxxy_1, \
                                         tg_xyyyyyz_xxxz_0, tg_xyyyyyz_xxxz_1, tg_xyyyyyz_xxy_1, tg_xyyyyyz_xxyy_0, \
                                         tg_xyyyyyz_xxyy_1, tg_xyyyyyz_xxyz_0, tg_xyyyyyz_xxyz_1, tg_xyyyyyz_xxz_1, \
                                         tg_xyyyyyz_xxzz_0, tg_xyyyyyz_xxzz_1, tg_xyyyyyz_xyy_1, tg_xyyyyyz_xyyy_0, \
                                         tg_xyyyyyz_xyyy_1, tg_xyyyyyz_xyyz_0, tg_xyyyyyz_xyyz_1, tg_xyyyyyz_xyz_1, \
                                         tg_xyyyyyz_xyzz_0, tg_xyyyyyz_xyzz_1, tg_xyyyyyz_xzz_1, tg_xyyyyyz_yyy_1, \
                                         tg_xyyyyyz_yyz_1, tg_xyyyyyz_yzz_1, tg_xyzzzz_xyyy_0, tg_xyzzzz_xyyy_1, \
                                         tg_xyzzzz_xyyz_0, tg_xyzzzz_xyyz_1, tg_xyzzzz_xyzz_0, tg_xyzzzz_xyzz_1, \
                                         tg_xyzzzz_xzzz_0, tg_xyzzzz_xzzz_1, tg_xyzzzz_yyyy_0, tg_xyzzzz_yyyy_1, \
                                         tg_xyzzzz_yyyz_0, tg_xyzzzz_yyyz_1, tg_xyzzzz_yyzz_0, tg_xyzzzz_yyzz_1, \
                                         tg_xyzzzz_yzzz_0, tg_xyzzzz_yzzz_1, tg_xyzzzz_zzzz_0, tg_xyzzzz_zzzz_1, \
                                         tg_xzzzzz_xxxx_0, tg_xzzzzz_xxxx_1, tg_xzzzzz_xxxy_0, tg_xzzzzz_xxxy_1, \
                                         tg_xzzzzz_xxxz_0, tg_xzzzzz_xxxz_1, tg_xzzzzz_xxyy_0, tg_xzzzzz_xxyy_1, \
                                         tg_xzzzzz_xxyz_0, tg_xzzzzz_xxyz_1, tg_xzzzzz_xxzz_0, tg_xzzzzz_xxzz_1, \
                                         tg_xzzzzz_xyyy_0, tg_xzzzzz_xyyy_1, tg_xzzzzz_xyyz_0, tg_xzzzzz_xyyz_1, \
                                         tg_xzzzzz_xyzz_0, tg_xzzzzz_xyzz_1, tg_xzzzzz_xzzz_0, tg_xzzzzz_xzzz_1, \
                                         tg_xzzzzz_yyyy_0, tg_xzzzzz_yyyy_1, tg_xzzzzz_yyyz_0, tg_xzzzzz_yyyz_1, \
                                         tg_xzzzzz_yyzz_0, tg_xzzzzz_yyzz_1, tg_xzzzzz_yzzz_0, tg_xzzzzz_yzzz_1, \
                                         tg_xzzzzz_zzzz_0, tg_xzzzzz_zzzz_1, tg_yyyyyy_xxxx_0, tg_yyyyyy_xxxx_1, \
                                         tg_yyyyyy_xxxy_0, tg_yyyyyy_xxxy_1, tg_yyyyyy_xxxz_0, tg_yyyyyy_xxxz_1, \
                                         tg_yyyyyy_xxyy_0, tg_yyyyyy_xxyy_1, tg_yyyyyy_xxyz_0, tg_yyyyyy_xxyz_1, \
                                         tg_yyyyyy_xxzz_0, tg_yyyyyy_xxzz_1, tg_yyyyyy_xyyy_0, tg_yyyyyy_xyyy_1, \
                                         tg_yyyyyy_xyyz_0, tg_yyyyyy_xyyz_1, tg_yyyyyy_xyzz_0, tg_yyyyyy_xyzz_1, \
                                         tg_yyyyyy_xzzz_0, tg_yyyyyy_xzzz_1, tg_yyyyyy_yyyy_0, tg_yyyyyy_yyyy_1, \
                                         tg_yyyyyy_yyyz_0, tg_yyyyyy_yyyz_1, tg_yyyyyy_yyzz_0, tg_yyyyyy_yyzz_1, \
                                         tg_yyyyyy_yzzz_0, tg_yyyyyy_yzzz_1, tg_yyyyyy_zzzz_0, tg_yyyyyy_zzzz_1, \
                                         tg_yyyyyz_xxxx_0, tg_yyyyyz_xxxx_1, tg_yyyyyz_xxxy_0, tg_yyyyyz_xxxy_1, \
                                         tg_yyyyyz_xxxz_0, tg_yyyyyz_xxxz_1, tg_yyyyyz_xxyy_0, tg_yyyyyz_xxyy_1, \
                                         tg_yyyyyz_xxyz_0, tg_yyyyyz_xxyz_1, tg_yyyyyz_xxzz_0, tg_yyyyyz_xxzz_1, \
                                         tg_yyyyyz_xyyy_0, tg_yyyyyz_xyyy_1, tg_yyyyyz_xyyz_0, tg_yyyyyz_xyyz_1, \
                                         tg_yyyyyz_xyzz_0, tg_yyyyyz_xyzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxxyzzzz_xyyy_0[j] = pb_x * tg_xxyzzzz_xyyy_0[j] + wp_x[j] * tg_xxyzzzz_xyyy_1[j] + fl1_fx * tg_xyzzzz_xyyy_0[j] - fl1_fx * fl1_fza * tg_xyzzzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxyzzzz_yyy_1[j];

                    tg_xxxyzzzz_xyyz_0[j] = pb_x * tg_xxyzzzz_xyyz_0[j] + wp_x[j] * tg_xxyzzzz_xyyz_1[j] + fl1_fx * tg_xyzzzz_xyyz_0[j] - fl1_fx * fl1_fza * tg_xyzzzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxyzzzz_yyz_1[j];

                    tg_xxxyzzzz_xyzz_0[j] = pb_x * tg_xxyzzzz_xyzz_0[j] + wp_x[j] * tg_xxyzzzz_xyzz_1[j] + fl1_fx * tg_xyzzzz_xyzz_0[j] - fl1_fx * fl1_fza * tg_xyzzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxyzzzz_yzz_1[j];

                    tg_xxxyzzzz_xzzz_0[j] = pb_x * tg_xxyzzzz_xzzz_0[j] + wp_x[j] * tg_xxyzzzz_xzzz_1[j] + fl1_fx * tg_xyzzzz_xzzz_0[j] - fl1_fx * fl1_fza * tg_xyzzzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxyzzzz_zzz_1[j];

                    tg_xxxyzzzz_yyyy_0[j] = pb_x * tg_xxyzzzz_yyyy_0[j] + wp_x[j] * tg_xxyzzzz_yyyy_1[j] + fl1_fx * tg_xyzzzz_yyyy_0[j] - fl1_fx * fl1_fza * tg_xyzzzz_yyyy_1[j];

                    tg_xxxyzzzz_yyyz_0[j] = pb_x * tg_xxyzzzz_yyyz_0[j] + wp_x[j] * tg_xxyzzzz_yyyz_1[j] + fl1_fx * tg_xyzzzz_yyyz_0[j] - fl1_fx * fl1_fza * tg_xyzzzz_yyyz_1[j];

                    tg_xxxyzzzz_yyzz_0[j] = pb_x * tg_xxyzzzz_yyzz_0[j] + wp_x[j] * tg_xxyzzzz_yyzz_1[j] + fl1_fx * tg_xyzzzz_yyzz_0[j] - fl1_fx * fl1_fza * tg_xyzzzz_yyzz_1[j];

                    tg_xxxyzzzz_yzzz_0[j] = pb_x * tg_xxyzzzz_yzzz_0[j] + wp_x[j] * tg_xxyzzzz_yzzz_1[j] + fl1_fx * tg_xyzzzz_yzzz_0[j] - fl1_fx * fl1_fza * tg_xyzzzz_yzzz_1[j];

                    tg_xxxyzzzz_zzzz_0[j] = pb_x * tg_xxyzzzz_zzzz_0[j] + wp_x[j] * tg_xxyzzzz_zzzz_1[j] + fl1_fx * tg_xyzzzz_zzzz_0[j] - fl1_fx * fl1_fza * tg_xyzzzz_zzzz_1[j];

                    tg_xxxzzzzz_xxxx_0[j] = pb_x * tg_xxzzzzz_xxxx_0[j] + wp_x[j] * tg_xxzzzzz_xxxx_1[j] + fl1_fx * tg_xzzzzz_xxxx_0[j] - fl1_fx * fl1_fza * tg_xzzzzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxzzzzz_xxx_1[j];

                    tg_xxxzzzzz_xxxy_0[j] = pb_x * tg_xxzzzzz_xxxy_0[j] + wp_x[j] * tg_xxzzzzz_xxxy_1[j] + fl1_fx * tg_xzzzzz_xxxy_0[j] - fl1_fx * fl1_fza * tg_xzzzzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxzzzzz_xxy_1[j];

                    tg_xxxzzzzz_xxxz_0[j] = pb_x * tg_xxzzzzz_xxxz_0[j] + wp_x[j] * tg_xxzzzzz_xxxz_1[j] + fl1_fx * tg_xzzzzz_xxxz_0[j] - fl1_fx * fl1_fza * tg_xzzzzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxzzzzz_xxz_1[j];

                    tg_xxxzzzzz_xxyy_0[j] = pb_x * tg_xxzzzzz_xxyy_0[j] + wp_x[j] * tg_xxzzzzz_xxyy_1[j] + fl1_fx * tg_xzzzzz_xxyy_0[j] - fl1_fx * fl1_fza * tg_xzzzzz_xxyy_1[j] + fl1_fxn * tg_xxzzzzz_xyy_1[j];

                    tg_xxxzzzzz_xxyz_0[j] = pb_x * tg_xxzzzzz_xxyz_0[j] + wp_x[j] * tg_xxzzzzz_xxyz_1[j] + fl1_fx * tg_xzzzzz_xxyz_0[j] - fl1_fx * fl1_fza * tg_xzzzzz_xxyz_1[j] + fl1_fxn * tg_xxzzzzz_xyz_1[j];

                    tg_xxxzzzzz_xxzz_0[j] = pb_x * tg_xxzzzzz_xxzz_0[j] + wp_x[j] * tg_xxzzzzz_xxzz_1[j] + fl1_fx * tg_xzzzzz_xxzz_0[j] - fl1_fx * fl1_fza * tg_xzzzzz_xxzz_1[j] + fl1_fxn * tg_xxzzzzz_xzz_1[j];

                    tg_xxxzzzzz_xyyy_0[j] = pb_x * tg_xxzzzzz_xyyy_0[j] + wp_x[j] * tg_xxzzzzz_xyyy_1[j] + fl1_fx * tg_xzzzzz_xyyy_0[j] - fl1_fx * fl1_fza * tg_xzzzzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxzzzzz_yyy_1[j];

                    tg_xxxzzzzz_xyyz_0[j] = pb_x * tg_xxzzzzz_xyyz_0[j] + wp_x[j] * tg_xxzzzzz_xyyz_1[j] + fl1_fx * tg_xzzzzz_xyyz_0[j] - fl1_fx * fl1_fza * tg_xzzzzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxzzzzz_yyz_1[j];

                    tg_xxxzzzzz_xyzz_0[j] = pb_x * tg_xxzzzzz_xyzz_0[j] + wp_x[j] * tg_xxzzzzz_xyzz_1[j] + fl1_fx * tg_xzzzzz_xyzz_0[j] - fl1_fx * fl1_fza * tg_xzzzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxzzzzz_yzz_1[j];

                    tg_xxxzzzzz_xzzz_0[j] = pb_x * tg_xxzzzzz_xzzz_0[j] + wp_x[j] * tg_xxzzzzz_xzzz_1[j] + fl1_fx * tg_xzzzzz_xzzz_0[j] - fl1_fx * fl1_fza * tg_xzzzzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxzzzzz_zzz_1[j];

                    tg_xxxzzzzz_yyyy_0[j] = pb_x * tg_xxzzzzz_yyyy_0[j] + wp_x[j] * tg_xxzzzzz_yyyy_1[j] + fl1_fx * tg_xzzzzz_yyyy_0[j] - fl1_fx * fl1_fza * tg_xzzzzz_yyyy_1[j];

                    tg_xxxzzzzz_yyyz_0[j] = pb_x * tg_xxzzzzz_yyyz_0[j] + wp_x[j] * tg_xxzzzzz_yyyz_1[j] + fl1_fx * tg_xzzzzz_yyyz_0[j] - fl1_fx * fl1_fza * tg_xzzzzz_yyyz_1[j];

                    tg_xxxzzzzz_yyzz_0[j] = pb_x * tg_xxzzzzz_yyzz_0[j] + wp_x[j] * tg_xxzzzzz_yyzz_1[j] + fl1_fx * tg_xzzzzz_yyzz_0[j] - fl1_fx * fl1_fza * tg_xzzzzz_yyzz_1[j];

                    tg_xxxzzzzz_yzzz_0[j] = pb_x * tg_xxzzzzz_yzzz_0[j] + wp_x[j] * tg_xxzzzzz_yzzz_1[j] + fl1_fx * tg_xzzzzz_yzzz_0[j] - fl1_fx * fl1_fza * tg_xzzzzz_yzzz_1[j];

                    tg_xxxzzzzz_zzzz_0[j] = pb_x * tg_xxzzzzz_zzzz_0[j] + wp_x[j] * tg_xxzzzzz_zzzz_1[j] + fl1_fx * tg_xzzzzz_zzzz_0[j] - fl1_fx * fl1_fza * tg_xzzzzz_zzzz_1[j];

                    tg_xxyyyyyy_xxxx_0[j] = pb_x * tg_xyyyyyy_xxxx_0[j] + wp_x[j] * tg_xyyyyyy_xxxx_1[j] + 0.5 * fl1_fx * tg_yyyyyy_xxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyy_xxxx_1[j] + 2.0 * fl1_fxn * tg_xyyyyyy_xxx_1[j];

                    tg_xxyyyyyy_xxxy_0[j] = pb_x * tg_xyyyyyy_xxxy_0[j] + wp_x[j] * tg_xyyyyyy_xxxy_1[j] + 0.5 * fl1_fx * tg_yyyyyy_xxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyy_xxxy_1[j] + 1.5 * fl1_fxn * tg_xyyyyyy_xxy_1[j];

                    tg_xxyyyyyy_xxxz_0[j] = pb_x * tg_xyyyyyy_xxxz_0[j] + wp_x[j] * tg_xyyyyyy_xxxz_1[j] + 0.5 * fl1_fx * tg_yyyyyy_xxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyy_xxxz_1[j] + 1.5 * fl1_fxn * tg_xyyyyyy_xxz_1[j];

                    tg_xxyyyyyy_xxyy_0[j] = pb_x * tg_xyyyyyy_xxyy_0[j] + wp_x[j] * tg_xyyyyyy_xxyy_1[j] + 0.5 * fl1_fx * tg_yyyyyy_xxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyy_xxyy_1[j] + fl1_fxn * tg_xyyyyyy_xyy_1[j];

                    tg_xxyyyyyy_xxyz_0[j] = pb_x * tg_xyyyyyy_xxyz_0[j] + wp_x[j] * tg_xyyyyyy_xxyz_1[j] + 0.5 * fl1_fx * tg_yyyyyy_xxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyy_xxyz_1[j] + fl1_fxn * tg_xyyyyyy_xyz_1[j];

                    tg_xxyyyyyy_xxzz_0[j] = pb_x * tg_xyyyyyy_xxzz_0[j] + wp_x[j] * tg_xyyyyyy_xxzz_1[j] + 0.5 * fl1_fx * tg_yyyyyy_xxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyy_xxzz_1[j] + fl1_fxn * tg_xyyyyyy_xzz_1[j];

                    tg_xxyyyyyy_xyyy_0[j] = pb_x * tg_xyyyyyy_xyyy_0[j] + wp_x[j] * tg_xyyyyyy_xyyy_1[j] + 0.5 * fl1_fx * tg_yyyyyy_xyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyy_xyyy_1[j] + 0.5 * fl1_fxn * tg_xyyyyyy_yyy_1[j];

                    tg_xxyyyyyy_xyyz_0[j] = pb_x * tg_xyyyyyy_xyyz_0[j] + wp_x[j] * tg_xyyyyyy_xyyz_1[j] + 0.5 * fl1_fx * tg_yyyyyy_xyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyy_xyyz_1[j] + 0.5 * fl1_fxn * tg_xyyyyyy_yyz_1[j];

                    tg_xxyyyyyy_xyzz_0[j] = pb_x * tg_xyyyyyy_xyzz_0[j] + wp_x[j] * tg_xyyyyyy_xyzz_1[j] + 0.5 * fl1_fx * tg_yyyyyy_xyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyy_xyzz_1[j] + 0.5 * fl1_fxn * tg_xyyyyyy_yzz_1[j];

                    tg_xxyyyyyy_xzzz_0[j] = pb_x * tg_xyyyyyy_xzzz_0[j] + wp_x[j] * tg_xyyyyyy_xzzz_1[j] + 0.5 * fl1_fx * tg_yyyyyy_xzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyy_xzzz_1[j] + 0.5 * fl1_fxn * tg_xyyyyyy_zzz_1[j];

                    tg_xxyyyyyy_yyyy_0[j] = pb_x * tg_xyyyyyy_yyyy_0[j] + wp_x[j] * tg_xyyyyyy_yyyy_1[j] + 0.5 * fl1_fx * tg_yyyyyy_yyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyy_yyyy_1[j];

                    tg_xxyyyyyy_yyyz_0[j] = pb_x * tg_xyyyyyy_yyyz_0[j] + wp_x[j] * tg_xyyyyyy_yyyz_1[j] + 0.5 * fl1_fx * tg_yyyyyy_yyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyy_yyyz_1[j];

                    tg_xxyyyyyy_yyzz_0[j] = pb_x * tg_xyyyyyy_yyzz_0[j] + wp_x[j] * tg_xyyyyyy_yyzz_1[j] + 0.5 * fl1_fx * tg_yyyyyy_yyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyy_yyzz_1[j];

                    tg_xxyyyyyy_yzzz_0[j] = pb_x * tg_xyyyyyy_yzzz_0[j] + wp_x[j] * tg_xyyyyyy_yzzz_1[j] + 0.5 * fl1_fx * tg_yyyyyy_yzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyy_yzzz_1[j];

                    tg_xxyyyyyy_zzzz_0[j] = pb_x * tg_xyyyyyy_zzzz_0[j] + wp_x[j] * tg_xyyyyyy_zzzz_1[j] + 0.5 * fl1_fx * tg_yyyyyy_zzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyy_zzzz_1[j];

                    tg_xxyyyyyz_xxxx_0[j] = pb_x * tg_xyyyyyz_xxxx_0[j] + wp_x[j] * tg_xyyyyyz_xxxx_1[j] + 0.5 * fl1_fx * tg_yyyyyz_xxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xyyyyyz_xxx_1[j];

                    tg_xxyyyyyz_xxxy_0[j] = pb_x * tg_xyyyyyz_xxxy_0[j] + wp_x[j] * tg_xyyyyyz_xxxy_1[j] + 0.5 * fl1_fx * tg_yyyyyz_xxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xyyyyyz_xxy_1[j];

                    tg_xxyyyyyz_xxxz_0[j] = pb_x * tg_xyyyyyz_xxxz_0[j] + wp_x[j] * tg_xyyyyyz_xxxz_1[j] + 0.5 * fl1_fx * tg_yyyyyz_xxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xyyyyyz_xxz_1[j];

                    tg_xxyyyyyz_xxyy_0[j] = pb_x * tg_xyyyyyz_xxyy_0[j] + wp_x[j] * tg_xyyyyyz_xxyy_1[j] + 0.5 * fl1_fx * tg_yyyyyz_xxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyz_xxyy_1[j] + fl1_fxn * tg_xyyyyyz_xyy_1[j];

                    tg_xxyyyyyz_xxyz_0[j] = pb_x * tg_xyyyyyz_xxyz_0[j] + wp_x[j] * tg_xyyyyyz_xxyz_1[j] + 0.5 * fl1_fx * tg_yyyyyz_xxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyz_xxyz_1[j] + fl1_fxn * tg_xyyyyyz_xyz_1[j];

                    tg_xxyyyyyz_xxzz_0[j] = pb_x * tg_xyyyyyz_xxzz_0[j] + wp_x[j] * tg_xyyyyyz_xxzz_1[j] + 0.5 * fl1_fx * tg_yyyyyz_xxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyz_xxzz_1[j] + fl1_fxn * tg_xyyyyyz_xzz_1[j];

                    tg_xxyyyyyz_xyyy_0[j] = pb_x * tg_xyyyyyz_xyyy_0[j] + wp_x[j] * tg_xyyyyyz_xyyy_1[j] + 0.5 * fl1_fx * tg_yyyyyz_xyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xyyyyyz_yyy_1[j];

                    tg_xxyyyyyz_xyyz_0[j] = pb_x * tg_xyyyyyz_xyyz_0[j] + wp_x[j] * tg_xyyyyyz_xyyz_1[j] + 0.5 * fl1_fx * tg_yyyyyz_xyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xyyyyyz_yyz_1[j];

                    tg_xxyyyyyz_xyzz_0[j] = pb_x * tg_xyyyyyz_xyzz_0[j] + wp_x[j] * tg_xyyyyyz_xyzz_1[j] + 0.5 * fl1_fx * tg_yyyyyz_xyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xyyyyyz_yzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSLSG_339_387(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (339,387)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {8, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {8, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_7_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_6_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_6_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_7_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_xyyyyyz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 339); 

                auto tg_xyyyyyz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 340); 

                auto tg_xyyyyyz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 341); 

                auto tg_xyyyyyz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 342); 

                auto tg_xyyyyyz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 343); 

                auto tg_xyyyyyz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 344); 

                auto tg_xyyyyzz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 345); 

                auto tg_xyyyyzz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 346); 

                auto tg_xyyyyzz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 347); 

                auto tg_xyyyyzz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 348); 

                auto tg_xyyyyzz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 349); 

                auto tg_xyyyyzz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 350); 

                auto tg_xyyyyzz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 351); 

                auto tg_xyyyyzz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 352); 

                auto tg_xyyyyzz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 353); 

                auto tg_xyyyyzz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 354); 

                auto tg_xyyyyzz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 355); 

                auto tg_xyyyyzz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 356); 

                auto tg_xyyyyzz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 357); 

                auto tg_xyyyyzz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 358); 

                auto tg_xyyyyzz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 359); 

                auto tg_xyyyzzz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 360); 

                auto tg_xyyyzzz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 361); 

                auto tg_xyyyzzz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 362); 

                auto tg_xyyyzzz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 363); 

                auto tg_xyyyzzz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 364); 

                auto tg_xyyyzzz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 365); 

                auto tg_xyyyzzz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 366); 

                auto tg_xyyyzzz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 367); 

                auto tg_xyyyzzz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 368); 

                auto tg_xyyyzzz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 369); 

                auto tg_xyyyzzz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 370); 

                auto tg_xyyyzzz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 371); 

                auto tg_xyyyzzz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 372); 

                auto tg_xyyyzzz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 373); 

                auto tg_xyyyzzz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 374); 

                auto tg_xyyzzzz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 375); 

                auto tg_xyyzzzz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 376); 

                auto tg_xyyzzzz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 377); 

                auto tg_xyyzzzz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 378); 

                auto tg_xyyzzzz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 379); 

                auto tg_xyyzzzz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 380); 

                auto tg_xyyzzzz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 381); 

                auto tg_xyyzzzz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 382); 

                auto tg_xyyzzzz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 383); 

                auto tg_xyyzzzz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 384); 

                auto tg_xyyzzzz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 385); 

                auto tg_xyyzzzz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 386); 

                auto tg_xyyyyyz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 339); 

                auto tg_xyyyyyz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 340); 

                auto tg_xyyyyyz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 341); 

                auto tg_xyyyyyz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 342); 

                auto tg_xyyyyyz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 343); 

                auto tg_xyyyyyz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 344); 

                auto tg_xyyyyzz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 345); 

                auto tg_xyyyyzz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 346); 

                auto tg_xyyyyzz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 347); 

                auto tg_xyyyyzz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 348); 

                auto tg_xyyyyzz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 349); 

                auto tg_xyyyyzz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 350); 

                auto tg_xyyyyzz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 351); 

                auto tg_xyyyyzz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 352); 

                auto tg_xyyyyzz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 353); 

                auto tg_xyyyyzz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 354); 

                auto tg_xyyyyzz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 355); 

                auto tg_xyyyyzz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 356); 

                auto tg_xyyyyzz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 357); 

                auto tg_xyyyyzz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 358); 

                auto tg_xyyyyzz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 359); 

                auto tg_xyyyzzz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 360); 

                auto tg_xyyyzzz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 361); 

                auto tg_xyyyzzz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 362); 

                auto tg_xyyyzzz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 363); 

                auto tg_xyyyzzz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 364); 

                auto tg_xyyyzzz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 365); 

                auto tg_xyyyzzz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 366); 

                auto tg_xyyyzzz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 367); 

                auto tg_xyyyzzz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 368); 

                auto tg_xyyyzzz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 369); 

                auto tg_xyyyzzz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 370); 

                auto tg_xyyyzzz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 371); 

                auto tg_xyyyzzz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 372); 

                auto tg_xyyyzzz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 373); 

                auto tg_xyyyzzz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 374); 

                auto tg_xyyzzzz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 375); 

                auto tg_xyyzzzz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 376); 

                auto tg_xyyzzzz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 377); 

                auto tg_xyyzzzz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 378); 

                auto tg_xyyzzzz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 379); 

                auto tg_xyyzzzz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 380); 

                auto tg_xyyzzzz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 381); 

                auto tg_xyyzzzz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 382); 

                auto tg_xyyzzzz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 383); 

                auto tg_xyyzzzz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 384); 

                auto tg_xyyzzzz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 385); 

                auto tg_xyyzzzz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 386); 

                auto tg_yyyyyz_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 339); 

                auto tg_yyyyyz_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 340); 

                auto tg_yyyyyz_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 341); 

                auto tg_yyyyyz_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 342); 

                auto tg_yyyyyz_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 343); 

                auto tg_yyyyyz_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 344); 

                auto tg_yyyyzz_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 345); 

                auto tg_yyyyzz_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 346); 

                auto tg_yyyyzz_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 347); 

                auto tg_yyyyzz_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 348); 

                auto tg_yyyyzz_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 349); 

                auto tg_yyyyzz_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 350); 

                auto tg_yyyyzz_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 351); 

                auto tg_yyyyzz_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 352); 

                auto tg_yyyyzz_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 353); 

                auto tg_yyyyzz_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 354); 

                auto tg_yyyyzz_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 355); 

                auto tg_yyyyzz_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 356); 

                auto tg_yyyyzz_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 357); 

                auto tg_yyyyzz_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 358); 

                auto tg_yyyyzz_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 359); 

                auto tg_yyyzzz_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 360); 

                auto tg_yyyzzz_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 361); 

                auto tg_yyyzzz_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 362); 

                auto tg_yyyzzz_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 363); 

                auto tg_yyyzzz_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 364); 

                auto tg_yyyzzz_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 365); 

                auto tg_yyyzzz_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 366); 

                auto tg_yyyzzz_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 367); 

                auto tg_yyyzzz_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 368); 

                auto tg_yyyzzz_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 369); 

                auto tg_yyyzzz_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 370); 

                auto tg_yyyzzz_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 371); 

                auto tg_yyyzzz_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 372); 

                auto tg_yyyzzz_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 373); 

                auto tg_yyyzzz_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 374); 

                auto tg_yyzzzz_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 375); 

                auto tg_yyzzzz_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 376); 

                auto tg_yyzzzz_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 377); 

                auto tg_yyzzzz_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 378); 

                auto tg_yyzzzz_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 379); 

                auto tg_yyzzzz_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 380); 

                auto tg_yyzzzz_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 381); 

                auto tg_yyzzzz_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 382); 

                auto tg_yyzzzz_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 383); 

                auto tg_yyzzzz_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 384); 

                auto tg_yyzzzz_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 385); 

                auto tg_yyzzzz_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 386); 

                auto tg_yyyyyz_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 339); 

                auto tg_yyyyyz_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 340); 

                auto tg_yyyyyz_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 341); 

                auto tg_yyyyyz_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 342); 

                auto tg_yyyyyz_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 343); 

                auto tg_yyyyyz_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 344); 

                auto tg_yyyyzz_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 345); 

                auto tg_yyyyzz_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 346); 

                auto tg_yyyyzz_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 347); 

                auto tg_yyyyzz_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 348); 

                auto tg_yyyyzz_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 349); 

                auto tg_yyyyzz_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 350); 

                auto tg_yyyyzz_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 351); 

                auto tg_yyyyzz_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 352); 

                auto tg_yyyyzz_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 353); 

                auto tg_yyyyzz_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 354); 

                auto tg_yyyyzz_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 355); 

                auto tg_yyyyzz_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 356); 

                auto tg_yyyyzz_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 357); 

                auto tg_yyyyzz_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 358); 

                auto tg_yyyyzz_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 359); 

                auto tg_yyyzzz_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 360); 

                auto tg_yyyzzz_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 361); 

                auto tg_yyyzzz_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 362); 

                auto tg_yyyzzz_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 363); 

                auto tg_yyyzzz_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 364); 

                auto tg_yyyzzz_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 365); 

                auto tg_yyyzzz_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 366); 

                auto tg_yyyzzz_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 367); 

                auto tg_yyyzzz_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 368); 

                auto tg_yyyzzz_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 369); 

                auto tg_yyyzzz_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 370); 

                auto tg_yyyzzz_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 371); 

                auto tg_yyyzzz_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 372); 

                auto tg_yyyzzz_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 373); 

                auto tg_yyyzzz_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 374); 

                auto tg_yyzzzz_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 375); 

                auto tg_yyzzzz_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 376); 

                auto tg_yyzzzz_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 377); 

                auto tg_yyzzzz_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 378); 

                auto tg_yyzzzz_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 379); 

                auto tg_yyzzzz_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 380); 

                auto tg_yyzzzz_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 381); 

                auto tg_yyzzzz_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 382); 

                auto tg_yyzzzz_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 383); 

                auto tg_yyzzzz_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 384); 

                auto tg_yyzzzz_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 385); 

                auto tg_yyzzzz_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 386); 

                auto tg_xyyyyyz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 229); 

                auto tg_xyyyyzz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 230); 

                auto tg_xyyyyzz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 231); 

                auto tg_xyyyyzz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 232); 

                auto tg_xyyyyzz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 233); 

                auto tg_xyyyyzz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 234); 

                auto tg_xyyyyzz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 235); 

                auto tg_xyyyyzz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 236); 

                auto tg_xyyyyzz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 237); 

                auto tg_xyyyyzz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 238); 

                auto tg_xyyyyzz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 239); 

                auto tg_xyyyzzz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 240); 

                auto tg_xyyyzzz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 241); 

                auto tg_xyyyzzz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 242); 

                auto tg_xyyyzzz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 243); 

                auto tg_xyyyzzz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 244); 

                auto tg_xyyyzzz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 245); 

                auto tg_xyyyzzz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 246); 

                auto tg_xyyyzzz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 247); 

                auto tg_xyyyzzz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 248); 

                auto tg_xyyyzzz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 249); 

                auto tg_xyyzzzz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 250); 

                auto tg_xyyzzzz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 251); 

                auto tg_xyyzzzz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 252); 

                auto tg_xyyzzzz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 253); 

                auto tg_xyyzzzz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 254); 

                auto tg_xyyzzzz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 255); 

                auto tg_xyyzzzz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 256); 

                auto tg_xyyzzzz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 257); 

                auto tg_xyyzzzz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 258); 

                auto tg_xyyzzzz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 259); 

                // set up pointers to integrals

                auto tg_xxyyyyyz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 339); 

                auto tg_xxyyyyyz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 340); 

                auto tg_xxyyyyyz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 341); 

                auto tg_xxyyyyyz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 342); 

                auto tg_xxyyyyyz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 343); 

                auto tg_xxyyyyyz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 344); 

                auto tg_xxyyyyzz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 345); 

                auto tg_xxyyyyzz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 346); 

                auto tg_xxyyyyzz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 347); 

                auto tg_xxyyyyzz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 348); 

                auto tg_xxyyyyzz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 349); 

                auto tg_xxyyyyzz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 350); 

                auto tg_xxyyyyzz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 351); 

                auto tg_xxyyyyzz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 352); 

                auto tg_xxyyyyzz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 353); 

                auto tg_xxyyyyzz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 354); 

                auto tg_xxyyyyzz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 355); 

                auto tg_xxyyyyzz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 356); 

                auto tg_xxyyyyzz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 357); 

                auto tg_xxyyyyzz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 358); 

                auto tg_xxyyyyzz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 359); 

                auto tg_xxyyyzzz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 360); 

                auto tg_xxyyyzzz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 361); 

                auto tg_xxyyyzzz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 362); 

                auto tg_xxyyyzzz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 363); 

                auto tg_xxyyyzzz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 364); 

                auto tg_xxyyyzzz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 365); 

                auto tg_xxyyyzzz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 366); 

                auto tg_xxyyyzzz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 367); 

                auto tg_xxyyyzzz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 368); 

                auto tg_xxyyyzzz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 369); 

                auto tg_xxyyyzzz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 370); 

                auto tg_xxyyyzzz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 371); 

                auto tg_xxyyyzzz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 372); 

                auto tg_xxyyyzzz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 373); 

                auto tg_xxyyyzzz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 374); 

                auto tg_xxyyzzzz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 375); 

                auto tg_xxyyzzzz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 376); 

                auto tg_xxyyzzzz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 377); 

                auto tg_xxyyzzzz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 378); 

                auto tg_xxyyzzzz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 379); 

                auto tg_xxyyzzzz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 380); 

                auto tg_xxyyzzzz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 381); 

                auto tg_xxyyzzzz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 382); 

                auto tg_xxyyzzzz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 383); 

                auto tg_xxyyzzzz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 384); 

                auto tg_xxyyzzzz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 385); 

                auto tg_xxyyzzzz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 386); 

                // Batch of Integrals (339,387)

                #pragma omp simd aligned(fxn, fza, tg_xxyyyyyz_xzzz_0, tg_xxyyyyyz_yyyy_0, tg_xxyyyyyz_yyyz_0, \
                                         tg_xxyyyyyz_yyzz_0, tg_xxyyyyyz_yzzz_0, tg_xxyyyyyz_zzzz_0, tg_xxyyyyzz_xxxx_0, \
                                         tg_xxyyyyzz_xxxy_0, tg_xxyyyyzz_xxxz_0, tg_xxyyyyzz_xxyy_0, tg_xxyyyyzz_xxyz_0, \
                                         tg_xxyyyyzz_xxzz_0, tg_xxyyyyzz_xyyy_0, tg_xxyyyyzz_xyyz_0, tg_xxyyyyzz_xyzz_0, \
                                         tg_xxyyyyzz_xzzz_0, tg_xxyyyyzz_yyyy_0, tg_xxyyyyzz_yyyz_0, tg_xxyyyyzz_yyzz_0, \
                                         tg_xxyyyyzz_yzzz_0, tg_xxyyyyzz_zzzz_0, tg_xxyyyzzz_xxxx_0, tg_xxyyyzzz_xxxy_0, \
                                         tg_xxyyyzzz_xxxz_0, tg_xxyyyzzz_xxyy_0, tg_xxyyyzzz_xxyz_0, tg_xxyyyzzz_xxzz_0, \
                                         tg_xxyyyzzz_xyyy_0, tg_xxyyyzzz_xyyz_0, tg_xxyyyzzz_xyzz_0, tg_xxyyyzzz_xzzz_0, \
                                         tg_xxyyyzzz_yyyy_0, tg_xxyyyzzz_yyyz_0, tg_xxyyyzzz_yyzz_0, tg_xxyyyzzz_yzzz_0, \
                                         tg_xxyyyzzz_zzzz_0, tg_xxyyzzzz_xxxx_0, tg_xxyyzzzz_xxxy_0, tg_xxyyzzzz_xxxz_0, \
                                         tg_xxyyzzzz_xxyy_0, tg_xxyyzzzz_xxyz_0, tg_xxyyzzzz_xxzz_0, tg_xxyyzzzz_xyyy_0, \
                                         tg_xxyyzzzz_xyyz_0, tg_xxyyzzzz_xyzz_0, tg_xxyyzzzz_xzzz_0, tg_xxyyzzzz_yyyy_0, \
                                         tg_xxyyzzzz_yyyz_0, tg_xyyyyyz_xzzz_0, tg_xyyyyyz_xzzz_1, tg_xyyyyyz_yyyy_0, \
                                         tg_xyyyyyz_yyyy_1, tg_xyyyyyz_yyyz_0, tg_xyyyyyz_yyyz_1, tg_xyyyyyz_yyzz_0, \
                                         tg_xyyyyyz_yyzz_1, tg_xyyyyyz_yzzz_0, tg_xyyyyyz_yzzz_1, tg_xyyyyyz_zzz_1, \
                                         tg_xyyyyyz_zzzz_0, tg_xyyyyyz_zzzz_1, tg_xyyyyzz_xxx_1, tg_xyyyyzz_xxxx_0, \
                                         tg_xyyyyzz_xxxx_1, tg_xyyyyzz_xxxy_0, tg_xyyyyzz_xxxy_1, tg_xyyyyzz_xxxz_0, \
                                         tg_xyyyyzz_xxxz_1, tg_xyyyyzz_xxy_1, tg_xyyyyzz_xxyy_0, tg_xyyyyzz_xxyy_1, \
                                         tg_xyyyyzz_xxyz_0, tg_xyyyyzz_xxyz_1, tg_xyyyyzz_xxz_1, tg_xyyyyzz_xxzz_0, \
                                         tg_xyyyyzz_xxzz_1, tg_xyyyyzz_xyy_1, tg_xyyyyzz_xyyy_0, tg_xyyyyzz_xyyy_1, \
                                         tg_xyyyyzz_xyyz_0, tg_xyyyyzz_xyyz_1, tg_xyyyyzz_xyz_1, tg_xyyyyzz_xyzz_0, \
                                         tg_xyyyyzz_xyzz_1, tg_xyyyyzz_xzz_1, tg_xyyyyzz_xzzz_0, tg_xyyyyzz_xzzz_1, \
                                         tg_xyyyyzz_yyy_1, tg_xyyyyzz_yyyy_0, tg_xyyyyzz_yyyy_1, tg_xyyyyzz_yyyz_0, \
                                         tg_xyyyyzz_yyyz_1, tg_xyyyyzz_yyz_1, tg_xyyyyzz_yyzz_0, tg_xyyyyzz_yyzz_1, \
                                         tg_xyyyyzz_yzz_1, tg_xyyyyzz_yzzz_0, tg_xyyyyzz_yzzz_1, tg_xyyyyzz_zzz_1, \
                                         tg_xyyyyzz_zzzz_0, tg_xyyyyzz_zzzz_1, tg_xyyyzzz_xxx_1, tg_xyyyzzz_xxxx_0, \
                                         tg_xyyyzzz_xxxx_1, tg_xyyyzzz_xxxy_0, tg_xyyyzzz_xxxy_1, tg_xyyyzzz_xxxz_0, \
                                         tg_xyyyzzz_xxxz_1, tg_xyyyzzz_xxy_1, tg_xyyyzzz_xxyy_0, tg_xyyyzzz_xxyy_1, \
                                         tg_xyyyzzz_xxyz_0, tg_xyyyzzz_xxyz_1, tg_xyyyzzz_xxz_1, tg_xyyyzzz_xxzz_0, \
                                         tg_xyyyzzz_xxzz_1, tg_xyyyzzz_xyy_1, tg_xyyyzzz_xyyy_0, tg_xyyyzzz_xyyy_1, \
                                         tg_xyyyzzz_xyyz_0, tg_xyyyzzz_xyyz_1, tg_xyyyzzz_xyz_1, tg_xyyyzzz_xyzz_0, \
                                         tg_xyyyzzz_xyzz_1, tg_xyyyzzz_xzz_1, tg_xyyyzzz_xzzz_0, tg_xyyyzzz_xzzz_1, \
                                         tg_xyyyzzz_yyy_1, tg_xyyyzzz_yyyy_0, tg_xyyyzzz_yyyy_1, tg_xyyyzzz_yyyz_0, \
                                         tg_xyyyzzz_yyyz_1, tg_xyyyzzz_yyz_1, tg_xyyyzzz_yyzz_0, tg_xyyyzzz_yyzz_1, \
                                         tg_xyyyzzz_yzz_1, tg_xyyyzzz_yzzz_0, tg_xyyyzzz_yzzz_1, tg_xyyyzzz_zzz_1, \
                                         tg_xyyyzzz_zzzz_0, tg_xyyyzzz_zzzz_1, tg_xyyzzzz_xxx_1, tg_xyyzzzz_xxxx_0, \
                                         tg_xyyzzzz_xxxx_1, tg_xyyzzzz_xxxy_0, tg_xyyzzzz_xxxy_1, tg_xyyzzzz_xxxz_0, \
                                         tg_xyyzzzz_xxxz_1, tg_xyyzzzz_xxy_1, tg_xyyzzzz_xxyy_0, tg_xyyzzzz_xxyy_1, \
                                         tg_xyyzzzz_xxyz_0, tg_xyyzzzz_xxyz_1, tg_xyyzzzz_xxz_1, tg_xyyzzzz_xxzz_0, \
                                         tg_xyyzzzz_xxzz_1, tg_xyyzzzz_xyy_1, tg_xyyzzzz_xyyy_0, tg_xyyzzzz_xyyy_1, \
                                         tg_xyyzzzz_xyyz_0, tg_xyyzzzz_xyyz_1, tg_xyyzzzz_xyz_1, tg_xyyzzzz_xyzz_0, \
                                         tg_xyyzzzz_xyzz_1, tg_xyyzzzz_xzz_1, tg_xyyzzzz_xzzz_0, tg_xyyzzzz_xzzz_1, \
                                         tg_xyyzzzz_yyy_1, tg_xyyzzzz_yyyy_0, tg_xyyzzzz_yyyy_1, tg_xyyzzzz_yyyz_0, \
                                         tg_xyyzzzz_yyyz_1, tg_xyyzzzz_yyz_1, tg_xyyzzzz_yzz_1, tg_xyyzzzz_zzz_1, \
                                         tg_yyyyyz_xzzz_0, tg_yyyyyz_xzzz_1, tg_yyyyyz_yyyy_0, tg_yyyyyz_yyyy_1, \
                                         tg_yyyyyz_yyyz_0, tg_yyyyyz_yyyz_1, tg_yyyyyz_yyzz_0, tg_yyyyyz_yyzz_1, \
                                         tg_yyyyyz_yzzz_0, tg_yyyyyz_yzzz_1, tg_yyyyyz_zzzz_0, tg_yyyyyz_zzzz_1, \
                                         tg_yyyyzz_xxxx_0, tg_yyyyzz_xxxx_1, tg_yyyyzz_xxxy_0, tg_yyyyzz_xxxy_1, \
                                         tg_yyyyzz_xxxz_0, tg_yyyyzz_xxxz_1, tg_yyyyzz_xxyy_0, tg_yyyyzz_xxyy_1, \
                                         tg_yyyyzz_xxyz_0, tg_yyyyzz_xxyz_1, tg_yyyyzz_xxzz_0, tg_yyyyzz_xxzz_1, \
                                         tg_yyyyzz_xyyy_0, tg_yyyyzz_xyyy_1, tg_yyyyzz_xyyz_0, tg_yyyyzz_xyyz_1, \
                                         tg_yyyyzz_xyzz_0, tg_yyyyzz_xyzz_1, tg_yyyyzz_xzzz_0, tg_yyyyzz_xzzz_1, \
                                         tg_yyyyzz_yyyy_0, tg_yyyyzz_yyyy_1, tg_yyyyzz_yyyz_0, tg_yyyyzz_yyyz_1, \
                                         tg_yyyyzz_yyzz_0, tg_yyyyzz_yyzz_1, tg_yyyyzz_yzzz_0, tg_yyyyzz_yzzz_1, \
                                         tg_yyyyzz_zzzz_0, tg_yyyyzz_zzzz_1, tg_yyyzzz_xxxx_0, tg_yyyzzz_xxxx_1, \
                                         tg_yyyzzz_xxxy_0, tg_yyyzzz_xxxy_1, tg_yyyzzz_xxxz_0, tg_yyyzzz_xxxz_1, \
                                         tg_yyyzzz_xxyy_0, tg_yyyzzz_xxyy_1, tg_yyyzzz_xxyz_0, tg_yyyzzz_xxyz_1, \
                                         tg_yyyzzz_xxzz_0, tg_yyyzzz_xxzz_1, tg_yyyzzz_xyyy_0, tg_yyyzzz_xyyy_1, \
                                         tg_yyyzzz_xyyz_0, tg_yyyzzz_xyyz_1, tg_yyyzzz_xyzz_0, tg_yyyzzz_xyzz_1, \
                                         tg_yyyzzz_xzzz_0, tg_yyyzzz_xzzz_1, tg_yyyzzz_yyyy_0, tg_yyyzzz_yyyy_1, \
                                         tg_yyyzzz_yyyz_0, tg_yyyzzz_yyyz_1, tg_yyyzzz_yyzz_0, tg_yyyzzz_yyzz_1, \
                                         tg_yyyzzz_yzzz_0, tg_yyyzzz_yzzz_1, tg_yyyzzz_zzzz_0, tg_yyyzzz_zzzz_1, \
                                         tg_yyzzzz_xxxx_0, tg_yyzzzz_xxxx_1, tg_yyzzzz_xxxy_0, tg_yyzzzz_xxxy_1, \
                                         tg_yyzzzz_xxxz_0, tg_yyzzzz_xxxz_1, tg_yyzzzz_xxyy_0, tg_yyzzzz_xxyy_1, \
                                         tg_yyzzzz_xxyz_0, tg_yyzzzz_xxyz_1, tg_yyzzzz_xxzz_0, tg_yyzzzz_xxzz_1, \
                                         tg_yyzzzz_xyyy_0, tg_yyzzzz_xyyy_1, tg_yyzzzz_xyyz_0, tg_yyzzzz_xyyz_1, \
                                         tg_yyzzzz_xyzz_0, tg_yyzzzz_xyzz_1, tg_yyzzzz_xzzz_0, tg_yyzzzz_xzzz_1, \
                                         tg_yyzzzz_yyyy_0, tg_yyzzzz_yyyy_1, tg_yyzzzz_yyyz_0, tg_yyzzzz_yyyz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxyyyyyz_xzzz_0[j] = pb_x * tg_xyyyyyz_xzzz_0[j] + wp_x[j] * tg_xyyyyyz_xzzz_1[j] + 0.5 * fl1_fx * tg_yyyyyz_xzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xyyyyyz_zzz_1[j];

                    tg_xxyyyyyz_yyyy_0[j] = pb_x * tg_xyyyyyz_yyyy_0[j] + wp_x[j] * tg_xyyyyyz_yyyy_1[j] + 0.5 * fl1_fx * tg_yyyyyz_yyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyz_yyyy_1[j];

                    tg_xxyyyyyz_yyyz_0[j] = pb_x * tg_xyyyyyz_yyyz_0[j] + wp_x[j] * tg_xyyyyyz_yyyz_1[j] + 0.5 * fl1_fx * tg_yyyyyz_yyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyz_yyyz_1[j];

                    tg_xxyyyyyz_yyzz_0[j] = pb_x * tg_xyyyyyz_yyzz_0[j] + wp_x[j] * tg_xyyyyyz_yyzz_1[j] + 0.5 * fl1_fx * tg_yyyyyz_yyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyz_yyzz_1[j];

                    tg_xxyyyyyz_yzzz_0[j] = pb_x * tg_xyyyyyz_yzzz_0[j] + wp_x[j] * tg_xyyyyyz_yzzz_1[j] + 0.5 * fl1_fx * tg_yyyyyz_yzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyz_yzzz_1[j];

                    tg_xxyyyyyz_zzzz_0[j] = pb_x * tg_xyyyyyz_zzzz_0[j] + wp_x[j] * tg_xyyyyyz_zzzz_1[j] + 0.5 * fl1_fx * tg_yyyyyz_zzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyz_zzzz_1[j];

                    tg_xxyyyyzz_xxxx_0[j] = pb_x * tg_xyyyyzz_xxxx_0[j] + wp_x[j] * tg_xyyyyzz_xxxx_1[j] + 0.5 * fl1_fx * tg_yyyyzz_xxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xyyyyzz_xxx_1[j];

                    tg_xxyyyyzz_xxxy_0[j] = pb_x * tg_xyyyyzz_xxxy_0[j] + wp_x[j] * tg_xyyyyzz_xxxy_1[j] + 0.5 * fl1_fx * tg_yyyyzz_xxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xyyyyzz_xxy_1[j];

                    tg_xxyyyyzz_xxxz_0[j] = pb_x * tg_xyyyyzz_xxxz_0[j] + wp_x[j] * tg_xyyyyzz_xxxz_1[j] + 0.5 * fl1_fx * tg_yyyyzz_xxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xyyyyzz_xxz_1[j];

                    tg_xxyyyyzz_xxyy_0[j] = pb_x * tg_xyyyyzz_xxyy_0[j] + wp_x[j] * tg_xyyyyzz_xxyy_1[j] + 0.5 * fl1_fx * tg_yyyyzz_xxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyzz_xxyy_1[j] + fl1_fxn * tg_xyyyyzz_xyy_1[j];

                    tg_xxyyyyzz_xxyz_0[j] = pb_x * tg_xyyyyzz_xxyz_0[j] + wp_x[j] * tg_xyyyyzz_xxyz_1[j] + 0.5 * fl1_fx * tg_yyyyzz_xxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyzz_xxyz_1[j] + fl1_fxn * tg_xyyyyzz_xyz_1[j];

                    tg_xxyyyyzz_xxzz_0[j] = pb_x * tg_xyyyyzz_xxzz_0[j] + wp_x[j] * tg_xyyyyzz_xxzz_1[j] + 0.5 * fl1_fx * tg_yyyyzz_xxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyzz_xxzz_1[j] + fl1_fxn * tg_xyyyyzz_xzz_1[j];

                    tg_xxyyyyzz_xyyy_0[j] = pb_x * tg_xyyyyzz_xyyy_0[j] + wp_x[j] * tg_xyyyyzz_xyyy_1[j] + 0.5 * fl1_fx * tg_yyyyzz_xyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xyyyyzz_yyy_1[j];

                    tg_xxyyyyzz_xyyz_0[j] = pb_x * tg_xyyyyzz_xyyz_0[j] + wp_x[j] * tg_xyyyyzz_xyyz_1[j] + 0.5 * fl1_fx * tg_yyyyzz_xyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xyyyyzz_yyz_1[j];

                    tg_xxyyyyzz_xyzz_0[j] = pb_x * tg_xyyyyzz_xyzz_0[j] + wp_x[j] * tg_xyyyyzz_xyzz_1[j] + 0.5 * fl1_fx * tg_yyyyzz_xyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xyyyyzz_yzz_1[j];

                    tg_xxyyyyzz_xzzz_0[j] = pb_x * tg_xyyyyzz_xzzz_0[j] + wp_x[j] * tg_xyyyyzz_xzzz_1[j] + 0.5 * fl1_fx * tg_yyyyzz_xzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xyyyyzz_zzz_1[j];

                    tg_xxyyyyzz_yyyy_0[j] = pb_x * tg_xyyyyzz_yyyy_0[j] + wp_x[j] * tg_xyyyyzz_yyyy_1[j] + 0.5 * fl1_fx * tg_yyyyzz_yyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyzz_yyyy_1[j];

                    tg_xxyyyyzz_yyyz_0[j] = pb_x * tg_xyyyyzz_yyyz_0[j] + wp_x[j] * tg_xyyyyzz_yyyz_1[j] + 0.5 * fl1_fx * tg_yyyyzz_yyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyzz_yyyz_1[j];

                    tg_xxyyyyzz_yyzz_0[j] = pb_x * tg_xyyyyzz_yyzz_0[j] + wp_x[j] * tg_xyyyyzz_yyzz_1[j] + 0.5 * fl1_fx * tg_yyyyzz_yyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyzz_yyzz_1[j];

                    tg_xxyyyyzz_yzzz_0[j] = pb_x * tg_xyyyyzz_yzzz_0[j] + wp_x[j] * tg_xyyyyzz_yzzz_1[j] + 0.5 * fl1_fx * tg_yyyyzz_yzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyzz_yzzz_1[j];

                    tg_xxyyyyzz_zzzz_0[j] = pb_x * tg_xyyyyzz_zzzz_0[j] + wp_x[j] * tg_xyyyyzz_zzzz_1[j] + 0.5 * fl1_fx * tg_yyyyzz_zzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyzz_zzzz_1[j];

                    tg_xxyyyzzz_xxxx_0[j] = pb_x * tg_xyyyzzz_xxxx_0[j] + wp_x[j] * tg_xyyyzzz_xxxx_1[j] + 0.5 * fl1_fx * tg_yyyzzz_xxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyzzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xyyyzzz_xxx_1[j];

                    tg_xxyyyzzz_xxxy_0[j] = pb_x * tg_xyyyzzz_xxxy_0[j] + wp_x[j] * tg_xyyyzzz_xxxy_1[j] + 0.5 * fl1_fx * tg_yyyzzz_xxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyzzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xyyyzzz_xxy_1[j];

                    tg_xxyyyzzz_xxxz_0[j] = pb_x * tg_xyyyzzz_xxxz_0[j] + wp_x[j] * tg_xyyyzzz_xxxz_1[j] + 0.5 * fl1_fx * tg_yyyzzz_xxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyzzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xyyyzzz_xxz_1[j];

                    tg_xxyyyzzz_xxyy_0[j] = pb_x * tg_xyyyzzz_xxyy_0[j] + wp_x[j] * tg_xyyyzzz_xxyy_1[j] + 0.5 * fl1_fx * tg_yyyzzz_xxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyzzz_xxyy_1[j] + fl1_fxn * tg_xyyyzzz_xyy_1[j];

                    tg_xxyyyzzz_xxyz_0[j] = pb_x * tg_xyyyzzz_xxyz_0[j] + wp_x[j] * tg_xyyyzzz_xxyz_1[j] + 0.5 * fl1_fx * tg_yyyzzz_xxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyzzz_xxyz_1[j] + fl1_fxn * tg_xyyyzzz_xyz_1[j];

                    tg_xxyyyzzz_xxzz_0[j] = pb_x * tg_xyyyzzz_xxzz_0[j] + wp_x[j] * tg_xyyyzzz_xxzz_1[j] + 0.5 * fl1_fx * tg_yyyzzz_xxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyzzz_xxzz_1[j] + fl1_fxn * tg_xyyyzzz_xzz_1[j];

                    tg_xxyyyzzz_xyyy_0[j] = pb_x * tg_xyyyzzz_xyyy_0[j] + wp_x[j] * tg_xyyyzzz_xyyy_1[j] + 0.5 * fl1_fx * tg_yyyzzz_xyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyzzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xyyyzzz_yyy_1[j];

                    tg_xxyyyzzz_xyyz_0[j] = pb_x * tg_xyyyzzz_xyyz_0[j] + wp_x[j] * tg_xyyyzzz_xyyz_1[j] + 0.5 * fl1_fx * tg_yyyzzz_xyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyzzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xyyyzzz_yyz_1[j];

                    tg_xxyyyzzz_xyzz_0[j] = pb_x * tg_xyyyzzz_xyzz_0[j] + wp_x[j] * tg_xyyyzzz_xyzz_1[j] + 0.5 * fl1_fx * tg_yyyzzz_xyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xyyyzzz_yzz_1[j];

                    tg_xxyyyzzz_xzzz_0[j] = pb_x * tg_xyyyzzz_xzzz_0[j] + wp_x[j] * tg_xyyyzzz_xzzz_1[j] + 0.5 * fl1_fx * tg_yyyzzz_xzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyzzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xyyyzzz_zzz_1[j];

                    tg_xxyyyzzz_yyyy_0[j] = pb_x * tg_xyyyzzz_yyyy_0[j] + wp_x[j] * tg_xyyyzzz_yyyy_1[j] + 0.5 * fl1_fx * tg_yyyzzz_yyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyzzz_yyyy_1[j];

                    tg_xxyyyzzz_yyyz_0[j] = pb_x * tg_xyyyzzz_yyyz_0[j] + wp_x[j] * tg_xyyyzzz_yyyz_1[j] + 0.5 * fl1_fx * tg_yyyzzz_yyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyzzz_yyyz_1[j];

                    tg_xxyyyzzz_yyzz_0[j] = pb_x * tg_xyyyzzz_yyzz_0[j] + wp_x[j] * tg_xyyyzzz_yyzz_1[j] + 0.5 * fl1_fx * tg_yyyzzz_yyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyzzz_yyzz_1[j];

                    tg_xxyyyzzz_yzzz_0[j] = pb_x * tg_xyyyzzz_yzzz_0[j] + wp_x[j] * tg_xyyyzzz_yzzz_1[j] + 0.5 * fl1_fx * tg_yyyzzz_yzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyzzz_yzzz_1[j];

                    tg_xxyyyzzz_zzzz_0[j] = pb_x * tg_xyyyzzz_zzzz_0[j] + wp_x[j] * tg_xyyyzzz_zzzz_1[j] + 0.5 * fl1_fx * tg_yyyzzz_zzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyzzz_zzzz_1[j];

                    tg_xxyyzzzz_xxxx_0[j] = pb_x * tg_xyyzzzz_xxxx_0[j] + wp_x[j] * tg_xyyzzzz_xxxx_1[j] + 0.5 * fl1_fx * tg_yyzzzz_xxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzzzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xyyzzzz_xxx_1[j];

                    tg_xxyyzzzz_xxxy_0[j] = pb_x * tg_xyyzzzz_xxxy_0[j] + wp_x[j] * tg_xyyzzzz_xxxy_1[j] + 0.5 * fl1_fx * tg_yyzzzz_xxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzzzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xyyzzzz_xxy_1[j];

                    tg_xxyyzzzz_xxxz_0[j] = pb_x * tg_xyyzzzz_xxxz_0[j] + wp_x[j] * tg_xyyzzzz_xxxz_1[j] + 0.5 * fl1_fx * tg_yyzzzz_xxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzzzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xyyzzzz_xxz_1[j];

                    tg_xxyyzzzz_xxyy_0[j] = pb_x * tg_xyyzzzz_xxyy_0[j] + wp_x[j] * tg_xyyzzzz_xxyy_1[j] + 0.5 * fl1_fx * tg_yyzzzz_xxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzzzz_xxyy_1[j] + fl1_fxn * tg_xyyzzzz_xyy_1[j];

                    tg_xxyyzzzz_xxyz_0[j] = pb_x * tg_xyyzzzz_xxyz_0[j] + wp_x[j] * tg_xyyzzzz_xxyz_1[j] + 0.5 * fl1_fx * tg_yyzzzz_xxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzzzz_xxyz_1[j] + fl1_fxn * tg_xyyzzzz_xyz_1[j];

                    tg_xxyyzzzz_xxzz_0[j] = pb_x * tg_xyyzzzz_xxzz_0[j] + wp_x[j] * tg_xyyzzzz_xxzz_1[j] + 0.5 * fl1_fx * tg_yyzzzz_xxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzzzz_xxzz_1[j] + fl1_fxn * tg_xyyzzzz_xzz_1[j];

                    tg_xxyyzzzz_xyyy_0[j] = pb_x * tg_xyyzzzz_xyyy_0[j] + wp_x[j] * tg_xyyzzzz_xyyy_1[j] + 0.5 * fl1_fx * tg_yyzzzz_xyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzzzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xyyzzzz_yyy_1[j];

                    tg_xxyyzzzz_xyyz_0[j] = pb_x * tg_xyyzzzz_xyyz_0[j] + wp_x[j] * tg_xyyzzzz_xyyz_1[j] + 0.5 * fl1_fx * tg_yyzzzz_xyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzzzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xyyzzzz_yyz_1[j];

                    tg_xxyyzzzz_xyzz_0[j] = pb_x * tg_xyyzzzz_xyzz_0[j] + wp_x[j] * tg_xyyzzzz_xyzz_1[j] + 0.5 * fl1_fx * tg_yyzzzz_xyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xyyzzzz_yzz_1[j];

                    tg_xxyyzzzz_xzzz_0[j] = pb_x * tg_xyyzzzz_xzzz_0[j] + wp_x[j] * tg_xyyzzzz_xzzz_1[j] + 0.5 * fl1_fx * tg_yyzzzz_xzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzzzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xyyzzzz_zzz_1[j];

                    tg_xxyyzzzz_yyyy_0[j] = pb_x * tg_xyyzzzz_yyyy_0[j] + wp_x[j] * tg_xyyzzzz_yyyy_1[j] + 0.5 * fl1_fx * tg_yyzzzz_yyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzzzz_yyyy_1[j];

                    tg_xxyyzzzz_yyyz_0[j] = pb_x * tg_xyyzzzz_yyyz_0[j] + wp_x[j] * tg_xyyzzzz_yyyz_1[j] + 0.5 * fl1_fx * tg_yyzzzz_yyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzzzz_yyyz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSLSG_387_435(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (387,435)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {8, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {8, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_7_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_6_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_6_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_7_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_xyyzzzz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 387); 

                auto tg_xyyzzzz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 388); 

                auto tg_xyyzzzz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 389); 

                auto tg_xyzzzzz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 390); 

                auto tg_xyzzzzz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 391); 

                auto tg_xyzzzzz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 392); 

                auto tg_xyzzzzz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 393); 

                auto tg_xyzzzzz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 394); 

                auto tg_xyzzzzz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 395); 

                auto tg_xyzzzzz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 396); 

                auto tg_xyzzzzz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 397); 

                auto tg_xyzzzzz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 398); 

                auto tg_xyzzzzz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 399); 

                auto tg_xyzzzzz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 400); 

                auto tg_xyzzzzz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 401); 

                auto tg_xyzzzzz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 402); 

                auto tg_xyzzzzz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 403); 

                auto tg_xyzzzzz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 404); 

                auto tg_xzzzzzz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 405); 

                auto tg_xzzzzzz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 406); 

                auto tg_xzzzzzz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 407); 

                auto tg_xzzzzzz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 408); 

                auto tg_xzzzzzz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 409); 

                auto tg_xzzzzzz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 410); 

                auto tg_xzzzzzz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 411); 

                auto tg_xzzzzzz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 412); 

                auto tg_xzzzzzz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 413); 

                auto tg_xzzzzzz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 414); 

                auto tg_xzzzzzz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 415); 

                auto tg_xzzzzzz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 416); 

                auto tg_xzzzzzz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 417); 

                auto tg_xzzzzzz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 418); 

                auto tg_xzzzzzz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 419); 

                auto tg_yyyyyyy_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 420); 

                auto tg_yyyyyyy_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 421); 

                auto tg_yyyyyyy_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 422); 

                auto tg_yyyyyyy_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 423); 

                auto tg_yyyyyyy_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 424); 

                auto tg_yyyyyyy_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 425); 

                auto tg_yyyyyyy_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 426); 

                auto tg_yyyyyyy_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 427); 

                auto tg_yyyyyyy_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 428); 

                auto tg_yyyyyyy_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 429); 

                auto tg_yyyyyyy_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 430); 

                auto tg_yyyyyyy_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 431); 

                auto tg_yyyyyyy_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 432); 

                auto tg_yyyyyyy_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 433); 

                auto tg_yyyyyyy_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 434); 

                auto tg_xyyzzzz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 387); 

                auto tg_xyyzzzz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 388); 

                auto tg_xyyzzzz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 389); 

                auto tg_xyzzzzz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 390); 

                auto tg_xyzzzzz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 391); 

                auto tg_xyzzzzz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 392); 

                auto tg_xyzzzzz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 393); 

                auto tg_xyzzzzz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 394); 

                auto tg_xyzzzzz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 395); 

                auto tg_xyzzzzz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 396); 

                auto tg_xyzzzzz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 397); 

                auto tg_xyzzzzz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 398); 

                auto tg_xyzzzzz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 399); 

                auto tg_xyzzzzz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 400); 

                auto tg_xyzzzzz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 401); 

                auto tg_xyzzzzz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 402); 

                auto tg_xyzzzzz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 403); 

                auto tg_xyzzzzz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 404); 

                auto tg_xzzzzzz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 405); 

                auto tg_xzzzzzz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 406); 

                auto tg_xzzzzzz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 407); 

                auto tg_xzzzzzz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 408); 

                auto tg_xzzzzzz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 409); 

                auto tg_xzzzzzz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 410); 

                auto tg_xzzzzzz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 411); 

                auto tg_xzzzzzz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 412); 

                auto tg_xzzzzzz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 413); 

                auto tg_xzzzzzz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 414); 

                auto tg_xzzzzzz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 415); 

                auto tg_xzzzzzz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 416); 

                auto tg_xzzzzzz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 417); 

                auto tg_xzzzzzz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 418); 

                auto tg_xzzzzzz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 419); 

                auto tg_yyyyyyy_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 420); 

                auto tg_yyyyyyy_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 421); 

                auto tg_yyyyyyy_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 422); 

                auto tg_yyyyyyy_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 423); 

                auto tg_yyyyyyy_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 424); 

                auto tg_yyyyyyy_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 425); 

                auto tg_yyyyyyy_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 426); 

                auto tg_yyyyyyy_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 427); 

                auto tg_yyyyyyy_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 428); 

                auto tg_yyyyyyy_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 429); 

                auto tg_yyyyyyy_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 430); 

                auto tg_yyyyyyy_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 431); 

                auto tg_yyyyyyy_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 432); 

                auto tg_yyyyyyy_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 433); 

                auto tg_yyyyyyy_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 434); 

                auto tg_yyzzzz_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 387); 

                auto tg_yyzzzz_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 388); 

                auto tg_yyzzzz_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 389); 

                auto tg_yzzzzz_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 390); 

                auto tg_yzzzzz_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 391); 

                auto tg_yzzzzz_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 392); 

                auto tg_yzzzzz_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 393); 

                auto tg_yzzzzz_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 394); 

                auto tg_yzzzzz_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 395); 

                auto tg_yzzzzz_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 396); 

                auto tg_yzzzzz_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 397); 

                auto tg_yzzzzz_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 398); 

                auto tg_yzzzzz_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 399); 

                auto tg_yzzzzz_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 400); 

                auto tg_yzzzzz_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 401); 

                auto tg_yzzzzz_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 402); 

                auto tg_yzzzzz_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 403); 

                auto tg_yzzzzz_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 404); 

                auto tg_zzzzzz_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 405); 

                auto tg_zzzzzz_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 406); 

                auto tg_zzzzzz_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 407); 

                auto tg_zzzzzz_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 408); 

                auto tg_zzzzzz_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 409); 

                auto tg_zzzzzz_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 410); 

                auto tg_zzzzzz_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 411); 

                auto tg_zzzzzz_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 412); 

                auto tg_zzzzzz_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 413); 

                auto tg_zzzzzz_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 414); 

                auto tg_zzzzzz_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 415); 

                auto tg_zzzzzz_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 416); 

                auto tg_zzzzzz_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 417); 

                auto tg_zzzzzz_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 418); 

                auto tg_zzzzzz_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 419); 

                auto tg_yyzzzz_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 387); 

                auto tg_yyzzzz_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 388); 

                auto tg_yyzzzz_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 389); 

                auto tg_yzzzzz_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 390); 

                auto tg_yzzzzz_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 391); 

                auto tg_yzzzzz_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 392); 

                auto tg_yzzzzz_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 393); 

                auto tg_yzzzzz_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 394); 

                auto tg_yzzzzz_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 395); 

                auto tg_yzzzzz_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 396); 

                auto tg_yzzzzz_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 397); 

                auto tg_yzzzzz_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 398); 

                auto tg_yzzzzz_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 399); 

                auto tg_yzzzzz_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 400); 

                auto tg_yzzzzz_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 401); 

                auto tg_yzzzzz_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 402); 

                auto tg_yzzzzz_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 403); 

                auto tg_yzzzzz_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 404); 

                auto tg_zzzzzz_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 405); 

                auto tg_zzzzzz_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 406); 

                auto tg_zzzzzz_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 407); 

                auto tg_zzzzzz_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 408); 

                auto tg_zzzzzz_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 409); 

                auto tg_zzzzzz_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 410); 

                auto tg_zzzzzz_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 411); 

                auto tg_zzzzzz_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 412); 

                auto tg_zzzzzz_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 413); 

                auto tg_zzzzzz_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 414); 

                auto tg_zzzzzz_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 415); 

                auto tg_zzzzzz_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 416); 

                auto tg_zzzzzz_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 417); 

                auto tg_zzzzzz_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 418); 

                auto tg_zzzzzz_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 419); 

                auto tg_xyzzzzz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 260); 

                auto tg_xyzzzzz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 261); 

                auto tg_xyzzzzz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 262); 

                auto tg_xyzzzzz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 263); 

                auto tg_xyzzzzz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 264); 

                auto tg_xyzzzzz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 265); 

                auto tg_xyzzzzz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 266); 

                auto tg_xyzzzzz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 267); 

                auto tg_xyzzzzz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 268); 

                auto tg_xyzzzzz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 269); 

                auto tg_xzzzzzz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 270); 

                auto tg_xzzzzzz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 271); 

                auto tg_xzzzzzz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 272); 

                auto tg_xzzzzzz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 273); 

                auto tg_xzzzzzz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 274); 

                auto tg_xzzzzzz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 275); 

                auto tg_xzzzzzz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 276); 

                auto tg_xzzzzzz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 277); 

                auto tg_xzzzzzz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 278); 

                auto tg_xzzzzzz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 279); 

                auto tg_yyyyyyy_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 280); 

                auto tg_yyyyyyy_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 281); 

                auto tg_yyyyyyy_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 282); 

                auto tg_yyyyyyy_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 283); 

                auto tg_yyyyyyy_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 284); 

                auto tg_yyyyyyy_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 285); 

                auto tg_yyyyyyy_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 286); 

                auto tg_yyyyyyy_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 287); 

                auto tg_yyyyyyy_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 288); 

                auto tg_yyyyyyy_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 289); 

                // set up pointers to integrals

                auto tg_xxyyzzzz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 387); 

                auto tg_xxyyzzzz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 388); 

                auto tg_xxyyzzzz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 389); 

                auto tg_xxyzzzzz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 390); 

                auto tg_xxyzzzzz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 391); 

                auto tg_xxyzzzzz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 392); 

                auto tg_xxyzzzzz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 393); 

                auto tg_xxyzzzzz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 394); 

                auto tg_xxyzzzzz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 395); 

                auto tg_xxyzzzzz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 396); 

                auto tg_xxyzzzzz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 397); 

                auto tg_xxyzzzzz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 398); 

                auto tg_xxyzzzzz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 399); 

                auto tg_xxyzzzzz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 400); 

                auto tg_xxyzzzzz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 401); 

                auto tg_xxyzzzzz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 402); 

                auto tg_xxyzzzzz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 403); 

                auto tg_xxyzzzzz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 404); 

                auto tg_xxzzzzzz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 405); 

                auto tg_xxzzzzzz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 406); 

                auto tg_xxzzzzzz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 407); 

                auto tg_xxzzzzzz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 408); 

                auto tg_xxzzzzzz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 409); 

                auto tg_xxzzzzzz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 410); 

                auto tg_xxzzzzzz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 411); 

                auto tg_xxzzzzzz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 412); 

                auto tg_xxzzzzzz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 413); 

                auto tg_xxzzzzzz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 414); 

                auto tg_xxzzzzzz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 415); 

                auto tg_xxzzzzzz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 416); 

                auto tg_xxzzzzzz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 417); 

                auto tg_xxzzzzzz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 418); 

                auto tg_xxzzzzzz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 419); 

                auto tg_xyyyyyyy_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 420); 

                auto tg_xyyyyyyy_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 421); 

                auto tg_xyyyyyyy_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 422); 

                auto tg_xyyyyyyy_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 423); 

                auto tg_xyyyyyyy_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 424); 

                auto tg_xyyyyyyy_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 425); 

                auto tg_xyyyyyyy_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 426); 

                auto tg_xyyyyyyy_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 427); 

                auto tg_xyyyyyyy_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 428); 

                auto tg_xyyyyyyy_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 429); 

                auto tg_xyyyyyyy_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 430); 

                auto tg_xyyyyyyy_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 431); 

                auto tg_xyyyyyyy_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 432); 

                auto tg_xyyyyyyy_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 433); 

                auto tg_xyyyyyyy_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 434); 

                // Batch of Integrals (387,435)

                #pragma omp simd aligned(fxn, fza, tg_xxyyzzzz_yyzz_0, tg_xxyyzzzz_yzzz_0, tg_xxyyzzzz_zzzz_0, \
                                         tg_xxyzzzzz_xxxx_0, tg_xxyzzzzz_xxxy_0, tg_xxyzzzzz_xxxz_0, tg_xxyzzzzz_xxyy_0, \
                                         tg_xxyzzzzz_xxyz_0, tg_xxyzzzzz_xxzz_0, tg_xxyzzzzz_xyyy_0, tg_xxyzzzzz_xyyz_0, \
                                         tg_xxyzzzzz_xyzz_0, tg_xxyzzzzz_xzzz_0, tg_xxyzzzzz_yyyy_0, tg_xxyzzzzz_yyyz_0, \
                                         tg_xxyzzzzz_yyzz_0, tg_xxyzzzzz_yzzz_0, tg_xxyzzzzz_zzzz_0, tg_xxzzzzzz_xxxx_0, \
                                         tg_xxzzzzzz_xxxy_0, tg_xxzzzzzz_xxxz_0, tg_xxzzzzzz_xxyy_0, tg_xxzzzzzz_xxyz_0, \
                                         tg_xxzzzzzz_xxzz_0, tg_xxzzzzzz_xyyy_0, tg_xxzzzzzz_xyyz_0, tg_xxzzzzzz_xyzz_0, \
                                         tg_xxzzzzzz_xzzz_0, tg_xxzzzzzz_yyyy_0, tg_xxzzzzzz_yyyz_0, tg_xxzzzzzz_yyzz_0, \
                                         tg_xxzzzzzz_yzzz_0, tg_xxzzzzzz_zzzz_0, tg_xyyyyyyy_xxxx_0, tg_xyyyyyyy_xxxy_0, \
                                         tg_xyyyyyyy_xxxz_0, tg_xyyyyyyy_xxyy_0, tg_xyyyyyyy_xxyz_0, tg_xyyyyyyy_xxzz_0, \
                                         tg_xyyyyyyy_xyyy_0, tg_xyyyyyyy_xyyz_0, tg_xyyyyyyy_xyzz_0, tg_xyyyyyyy_xzzz_0, \
                                         tg_xyyyyyyy_yyyy_0, tg_xyyyyyyy_yyyz_0, tg_xyyyyyyy_yyzz_0, tg_xyyyyyyy_yzzz_0, \
                                         tg_xyyyyyyy_zzzz_0, tg_xyyzzzz_yyzz_0, tg_xyyzzzz_yyzz_1, tg_xyyzzzz_yzzz_0, \
                                         tg_xyyzzzz_yzzz_1, tg_xyyzzzz_zzzz_0, tg_xyyzzzz_zzzz_1, tg_xyzzzzz_xxx_1, \
                                         tg_xyzzzzz_xxxx_0, tg_xyzzzzz_xxxx_1, tg_xyzzzzz_xxxy_0, tg_xyzzzzz_xxxy_1, \
                                         tg_xyzzzzz_xxxz_0, tg_xyzzzzz_xxxz_1, tg_xyzzzzz_xxy_1, tg_xyzzzzz_xxyy_0, \
                                         tg_xyzzzzz_xxyy_1, tg_xyzzzzz_xxyz_0, tg_xyzzzzz_xxyz_1, tg_xyzzzzz_xxz_1, \
                                         tg_xyzzzzz_xxzz_0, tg_xyzzzzz_xxzz_1, tg_xyzzzzz_xyy_1, tg_xyzzzzz_xyyy_0, \
                                         tg_xyzzzzz_xyyy_1, tg_xyzzzzz_xyyz_0, tg_xyzzzzz_xyyz_1, tg_xyzzzzz_xyz_1, \
                                         tg_xyzzzzz_xyzz_0, tg_xyzzzzz_xyzz_1, tg_xyzzzzz_xzz_1, tg_xyzzzzz_xzzz_0, \
                                         tg_xyzzzzz_xzzz_1, tg_xyzzzzz_yyy_1, tg_xyzzzzz_yyyy_0, tg_xyzzzzz_yyyy_1, \
                                         tg_xyzzzzz_yyyz_0, tg_xyzzzzz_yyyz_1, tg_xyzzzzz_yyz_1, tg_xyzzzzz_yyzz_0, \
                                         tg_xyzzzzz_yyzz_1, tg_xyzzzzz_yzz_1, tg_xyzzzzz_yzzz_0, tg_xyzzzzz_yzzz_1, \
                                         tg_xyzzzzz_zzz_1, tg_xyzzzzz_zzzz_0, tg_xyzzzzz_zzzz_1, tg_xzzzzzz_xxx_1, \
                                         tg_xzzzzzz_xxxx_0, tg_xzzzzzz_xxxx_1, tg_xzzzzzz_xxxy_0, tg_xzzzzzz_xxxy_1, \
                                         tg_xzzzzzz_xxxz_0, tg_xzzzzzz_xxxz_1, tg_xzzzzzz_xxy_1, tg_xzzzzzz_xxyy_0, \
                                         tg_xzzzzzz_xxyy_1, tg_xzzzzzz_xxyz_0, tg_xzzzzzz_xxyz_1, tg_xzzzzzz_xxz_1, \
                                         tg_xzzzzzz_xxzz_0, tg_xzzzzzz_xxzz_1, tg_xzzzzzz_xyy_1, tg_xzzzzzz_xyyy_0, \
                                         tg_xzzzzzz_xyyy_1, tg_xzzzzzz_xyyz_0, tg_xzzzzzz_xyyz_1, tg_xzzzzzz_xyz_1, \
                                         tg_xzzzzzz_xyzz_0, tg_xzzzzzz_xyzz_1, tg_xzzzzzz_xzz_1, tg_xzzzzzz_xzzz_0, \
                                         tg_xzzzzzz_xzzz_1, tg_xzzzzzz_yyy_1, tg_xzzzzzz_yyyy_0, tg_xzzzzzz_yyyy_1, \
                                         tg_xzzzzzz_yyyz_0, tg_xzzzzzz_yyyz_1, tg_xzzzzzz_yyz_1, tg_xzzzzzz_yyzz_0, \
                                         tg_xzzzzzz_yyzz_1, tg_xzzzzzz_yzz_1, tg_xzzzzzz_yzzz_0, tg_xzzzzzz_yzzz_1, \
                                         tg_xzzzzzz_zzz_1, tg_xzzzzzz_zzzz_0, tg_xzzzzzz_zzzz_1, tg_yyyyyyy_xxx_1, \
                                         tg_yyyyyyy_xxxx_0, tg_yyyyyyy_xxxx_1, tg_yyyyyyy_xxxy_0, tg_yyyyyyy_xxxy_1, \
                                         tg_yyyyyyy_xxxz_0, tg_yyyyyyy_xxxz_1, tg_yyyyyyy_xxy_1, tg_yyyyyyy_xxyy_0, \
                                         tg_yyyyyyy_xxyy_1, tg_yyyyyyy_xxyz_0, tg_yyyyyyy_xxyz_1, tg_yyyyyyy_xxz_1, \
                                         tg_yyyyyyy_xxzz_0, tg_yyyyyyy_xxzz_1, tg_yyyyyyy_xyy_1, tg_yyyyyyy_xyyy_0, \
                                         tg_yyyyyyy_xyyy_1, tg_yyyyyyy_xyyz_0, tg_yyyyyyy_xyyz_1, tg_yyyyyyy_xyz_1, \
                                         tg_yyyyyyy_xyzz_0, tg_yyyyyyy_xyzz_1, tg_yyyyyyy_xzz_1, tg_yyyyyyy_xzzz_0, \
                                         tg_yyyyyyy_xzzz_1, tg_yyyyyyy_yyy_1, tg_yyyyyyy_yyyy_0, tg_yyyyyyy_yyyy_1, \
                                         tg_yyyyyyy_yyyz_0, tg_yyyyyyy_yyyz_1, tg_yyyyyyy_yyz_1, tg_yyyyyyy_yyzz_0, \
                                         tg_yyyyyyy_yyzz_1, tg_yyyyyyy_yzz_1, tg_yyyyyyy_yzzz_0, tg_yyyyyyy_yzzz_1, \
                                         tg_yyyyyyy_zzz_1, tg_yyyyyyy_zzzz_0, tg_yyyyyyy_zzzz_1, tg_yyzzzz_yyzz_0, \
                                         tg_yyzzzz_yyzz_1, tg_yyzzzz_yzzz_0, tg_yyzzzz_yzzz_1, tg_yyzzzz_zzzz_0, \
                                         tg_yyzzzz_zzzz_1, tg_yzzzzz_xxxx_0, tg_yzzzzz_xxxx_1, tg_yzzzzz_xxxy_0, \
                                         tg_yzzzzz_xxxy_1, tg_yzzzzz_xxxz_0, tg_yzzzzz_xxxz_1, tg_yzzzzz_xxyy_0, \
                                         tg_yzzzzz_xxyy_1, tg_yzzzzz_xxyz_0, tg_yzzzzz_xxyz_1, tg_yzzzzz_xxzz_0, \
                                         tg_yzzzzz_xxzz_1, tg_yzzzzz_xyyy_0, tg_yzzzzz_xyyy_1, tg_yzzzzz_xyyz_0, \
                                         tg_yzzzzz_xyyz_1, tg_yzzzzz_xyzz_0, tg_yzzzzz_xyzz_1, tg_yzzzzz_xzzz_0, \
                                         tg_yzzzzz_xzzz_1, tg_yzzzzz_yyyy_0, tg_yzzzzz_yyyy_1, tg_yzzzzz_yyyz_0, \
                                         tg_yzzzzz_yyyz_1, tg_yzzzzz_yyzz_0, tg_yzzzzz_yyzz_1, tg_yzzzzz_yzzz_0, \
                                         tg_yzzzzz_yzzz_1, tg_yzzzzz_zzzz_0, tg_yzzzzz_zzzz_1, tg_zzzzzz_xxxx_0, \
                                         tg_zzzzzz_xxxx_1, tg_zzzzzz_xxxy_0, tg_zzzzzz_xxxy_1, tg_zzzzzz_xxxz_0, \
                                         tg_zzzzzz_xxxz_1, tg_zzzzzz_xxyy_0, tg_zzzzzz_xxyy_1, tg_zzzzzz_xxyz_0, \
                                         tg_zzzzzz_xxyz_1, tg_zzzzzz_xxzz_0, tg_zzzzzz_xxzz_1, tg_zzzzzz_xyyy_0, \
                                         tg_zzzzzz_xyyy_1, tg_zzzzzz_xyyz_0, tg_zzzzzz_xyyz_1, tg_zzzzzz_xyzz_0, \
                                         tg_zzzzzz_xyzz_1, tg_zzzzzz_xzzz_0, tg_zzzzzz_xzzz_1, tg_zzzzzz_yyyy_0, \
                                         tg_zzzzzz_yyyy_1, tg_zzzzzz_yyyz_0, tg_zzzzzz_yyyz_1, tg_zzzzzz_yyzz_0, \
                                         tg_zzzzzz_yyzz_1, tg_zzzzzz_yzzz_0, tg_zzzzzz_yzzz_1, tg_zzzzzz_zzzz_0, \
                                         tg_zzzzzz_zzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxyyzzzz_yyzz_0[j] = pb_x * tg_xyyzzzz_yyzz_0[j] + wp_x[j] * tg_xyyzzzz_yyzz_1[j] + 0.5 * fl1_fx * tg_yyzzzz_yyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzzzz_yyzz_1[j];

                    tg_xxyyzzzz_yzzz_0[j] = pb_x * tg_xyyzzzz_yzzz_0[j] + wp_x[j] * tg_xyyzzzz_yzzz_1[j] + 0.5 * fl1_fx * tg_yyzzzz_yzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzzzz_yzzz_1[j];

                    tg_xxyyzzzz_zzzz_0[j] = pb_x * tg_xyyzzzz_zzzz_0[j] + wp_x[j] * tg_xyyzzzz_zzzz_1[j] + 0.5 * fl1_fx * tg_yyzzzz_zzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzzzz_zzzz_1[j];

                    tg_xxyzzzzz_xxxx_0[j] = pb_x * tg_xyzzzzz_xxxx_0[j] + wp_x[j] * tg_xyzzzzz_xxxx_1[j] + 0.5 * fl1_fx * tg_yzzzzz_xxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzzzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xyzzzzz_xxx_1[j];

                    tg_xxyzzzzz_xxxy_0[j] = pb_x * tg_xyzzzzz_xxxy_0[j] + wp_x[j] * tg_xyzzzzz_xxxy_1[j] + 0.5 * fl1_fx * tg_yzzzzz_xxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzzzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xyzzzzz_xxy_1[j];

                    tg_xxyzzzzz_xxxz_0[j] = pb_x * tg_xyzzzzz_xxxz_0[j] + wp_x[j] * tg_xyzzzzz_xxxz_1[j] + 0.5 * fl1_fx * tg_yzzzzz_xxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzzzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xyzzzzz_xxz_1[j];

                    tg_xxyzzzzz_xxyy_0[j] = pb_x * tg_xyzzzzz_xxyy_0[j] + wp_x[j] * tg_xyzzzzz_xxyy_1[j] + 0.5 * fl1_fx * tg_yzzzzz_xxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzzzz_xxyy_1[j] + fl1_fxn * tg_xyzzzzz_xyy_1[j];

                    tg_xxyzzzzz_xxyz_0[j] = pb_x * tg_xyzzzzz_xxyz_0[j] + wp_x[j] * tg_xyzzzzz_xxyz_1[j] + 0.5 * fl1_fx * tg_yzzzzz_xxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzzzz_xxyz_1[j] + fl1_fxn * tg_xyzzzzz_xyz_1[j];

                    tg_xxyzzzzz_xxzz_0[j] = pb_x * tg_xyzzzzz_xxzz_0[j] + wp_x[j] * tg_xyzzzzz_xxzz_1[j] + 0.5 * fl1_fx * tg_yzzzzz_xxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzzzz_xxzz_1[j] + fl1_fxn * tg_xyzzzzz_xzz_1[j];

                    tg_xxyzzzzz_xyyy_0[j] = pb_x * tg_xyzzzzz_xyyy_0[j] + wp_x[j] * tg_xyzzzzz_xyyy_1[j] + 0.5 * fl1_fx * tg_yzzzzz_xyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzzzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xyzzzzz_yyy_1[j];

                    tg_xxyzzzzz_xyyz_0[j] = pb_x * tg_xyzzzzz_xyyz_0[j] + wp_x[j] * tg_xyzzzzz_xyyz_1[j] + 0.5 * fl1_fx * tg_yzzzzz_xyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzzzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xyzzzzz_yyz_1[j];

                    tg_xxyzzzzz_xyzz_0[j] = pb_x * tg_xyzzzzz_xyzz_0[j] + wp_x[j] * tg_xyzzzzz_xyzz_1[j] + 0.5 * fl1_fx * tg_yzzzzz_xyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xyzzzzz_yzz_1[j];

                    tg_xxyzzzzz_xzzz_0[j] = pb_x * tg_xyzzzzz_xzzz_0[j] + wp_x[j] * tg_xyzzzzz_xzzz_1[j] + 0.5 * fl1_fx * tg_yzzzzz_xzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzzzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xyzzzzz_zzz_1[j];

                    tg_xxyzzzzz_yyyy_0[j] = pb_x * tg_xyzzzzz_yyyy_0[j] + wp_x[j] * tg_xyzzzzz_yyyy_1[j] + 0.5 * fl1_fx * tg_yzzzzz_yyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzzzz_yyyy_1[j];

                    tg_xxyzzzzz_yyyz_0[j] = pb_x * tg_xyzzzzz_yyyz_0[j] + wp_x[j] * tg_xyzzzzz_yyyz_1[j] + 0.5 * fl1_fx * tg_yzzzzz_yyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzzzz_yyyz_1[j];

                    tg_xxyzzzzz_yyzz_0[j] = pb_x * tg_xyzzzzz_yyzz_0[j] + wp_x[j] * tg_xyzzzzz_yyzz_1[j] + 0.5 * fl1_fx * tg_yzzzzz_yyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzzzz_yyzz_1[j];

                    tg_xxyzzzzz_yzzz_0[j] = pb_x * tg_xyzzzzz_yzzz_0[j] + wp_x[j] * tg_xyzzzzz_yzzz_1[j] + 0.5 * fl1_fx * tg_yzzzzz_yzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzzzz_yzzz_1[j];

                    tg_xxyzzzzz_zzzz_0[j] = pb_x * tg_xyzzzzz_zzzz_0[j] + wp_x[j] * tg_xyzzzzz_zzzz_1[j] + 0.5 * fl1_fx * tg_yzzzzz_zzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzzzz_zzzz_1[j];

                    tg_xxzzzzzz_xxxx_0[j] = pb_x * tg_xzzzzzz_xxxx_0[j] + wp_x[j] * tg_xzzzzzz_xxxx_1[j] + 0.5 * fl1_fx * tg_zzzzzz_xxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xzzzzzz_xxx_1[j];

                    tg_xxzzzzzz_xxxy_0[j] = pb_x * tg_xzzzzzz_xxxy_0[j] + wp_x[j] * tg_xzzzzzz_xxxy_1[j] + 0.5 * fl1_fx * tg_zzzzzz_xxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xzzzzzz_xxy_1[j];

                    tg_xxzzzzzz_xxxz_0[j] = pb_x * tg_xzzzzzz_xxxz_0[j] + wp_x[j] * tg_xzzzzzz_xxxz_1[j] + 0.5 * fl1_fx * tg_zzzzzz_xxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xzzzzzz_xxz_1[j];

                    tg_xxzzzzzz_xxyy_0[j] = pb_x * tg_xzzzzzz_xxyy_0[j] + wp_x[j] * tg_xzzzzzz_xxyy_1[j] + 0.5 * fl1_fx * tg_zzzzzz_xxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_xxyy_1[j] + fl1_fxn * tg_xzzzzzz_xyy_1[j];

                    tg_xxzzzzzz_xxyz_0[j] = pb_x * tg_xzzzzzz_xxyz_0[j] + wp_x[j] * tg_xzzzzzz_xxyz_1[j] + 0.5 * fl1_fx * tg_zzzzzz_xxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_xxyz_1[j] + fl1_fxn * tg_xzzzzzz_xyz_1[j];

                    tg_xxzzzzzz_xxzz_0[j] = pb_x * tg_xzzzzzz_xxzz_0[j] + wp_x[j] * tg_xzzzzzz_xxzz_1[j] + 0.5 * fl1_fx * tg_zzzzzz_xxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_xxzz_1[j] + fl1_fxn * tg_xzzzzzz_xzz_1[j];

                    tg_xxzzzzzz_xyyy_0[j] = pb_x * tg_xzzzzzz_xyyy_0[j] + wp_x[j] * tg_xzzzzzz_xyyy_1[j] + 0.5 * fl1_fx * tg_zzzzzz_xyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xzzzzzz_yyy_1[j];

                    tg_xxzzzzzz_xyyz_0[j] = pb_x * tg_xzzzzzz_xyyz_0[j] + wp_x[j] * tg_xzzzzzz_xyyz_1[j] + 0.5 * fl1_fx * tg_zzzzzz_xyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xzzzzzz_yyz_1[j];

                    tg_xxzzzzzz_xyzz_0[j] = pb_x * tg_xzzzzzz_xyzz_0[j] + wp_x[j] * tg_xzzzzzz_xyzz_1[j] + 0.5 * fl1_fx * tg_zzzzzz_xyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xzzzzzz_yzz_1[j];

                    tg_xxzzzzzz_xzzz_0[j] = pb_x * tg_xzzzzzz_xzzz_0[j] + wp_x[j] * tg_xzzzzzz_xzzz_1[j] + 0.5 * fl1_fx * tg_zzzzzz_xzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xzzzzzz_zzz_1[j];

                    tg_xxzzzzzz_yyyy_0[j] = pb_x * tg_xzzzzzz_yyyy_0[j] + wp_x[j] * tg_xzzzzzz_yyyy_1[j] + 0.5 * fl1_fx * tg_zzzzzz_yyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_yyyy_1[j];

                    tg_xxzzzzzz_yyyz_0[j] = pb_x * tg_xzzzzzz_yyyz_0[j] + wp_x[j] * tg_xzzzzzz_yyyz_1[j] + 0.5 * fl1_fx * tg_zzzzzz_yyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_yyyz_1[j];

                    tg_xxzzzzzz_yyzz_0[j] = pb_x * tg_xzzzzzz_yyzz_0[j] + wp_x[j] * tg_xzzzzzz_yyzz_1[j] + 0.5 * fl1_fx * tg_zzzzzz_yyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_yyzz_1[j];

                    tg_xxzzzzzz_yzzz_0[j] = pb_x * tg_xzzzzzz_yzzz_0[j] + wp_x[j] * tg_xzzzzzz_yzzz_1[j] + 0.5 * fl1_fx * tg_zzzzzz_yzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_yzzz_1[j];

                    tg_xxzzzzzz_zzzz_0[j] = pb_x * tg_xzzzzzz_zzzz_0[j] + wp_x[j] * tg_xzzzzzz_zzzz_1[j] + 0.5 * fl1_fx * tg_zzzzzz_zzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_zzzz_1[j];

                    tg_xyyyyyyy_xxxx_0[j] = pb_x * tg_yyyyyyy_xxxx_0[j] + wp_x[j] * tg_yyyyyyy_xxxx_1[j] + 2.0 * fl1_fxn * tg_yyyyyyy_xxx_1[j];

                    tg_xyyyyyyy_xxxy_0[j] = pb_x * tg_yyyyyyy_xxxy_0[j] + wp_x[j] * tg_yyyyyyy_xxxy_1[j] + 1.5 * fl1_fxn * tg_yyyyyyy_xxy_1[j];

                    tg_xyyyyyyy_xxxz_0[j] = pb_x * tg_yyyyyyy_xxxz_0[j] + wp_x[j] * tg_yyyyyyy_xxxz_1[j] + 1.5 * fl1_fxn * tg_yyyyyyy_xxz_1[j];

                    tg_xyyyyyyy_xxyy_0[j] = pb_x * tg_yyyyyyy_xxyy_0[j] + wp_x[j] * tg_yyyyyyy_xxyy_1[j] + fl1_fxn * tg_yyyyyyy_xyy_1[j];

                    tg_xyyyyyyy_xxyz_0[j] = pb_x * tg_yyyyyyy_xxyz_0[j] + wp_x[j] * tg_yyyyyyy_xxyz_1[j] + fl1_fxn * tg_yyyyyyy_xyz_1[j];

                    tg_xyyyyyyy_xxzz_0[j] = pb_x * tg_yyyyyyy_xxzz_0[j] + wp_x[j] * tg_yyyyyyy_xxzz_1[j] + fl1_fxn * tg_yyyyyyy_xzz_1[j];

                    tg_xyyyyyyy_xyyy_0[j] = pb_x * tg_yyyyyyy_xyyy_0[j] + wp_x[j] * tg_yyyyyyy_xyyy_1[j] + 0.5 * fl1_fxn * tg_yyyyyyy_yyy_1[j];

                    tg_xyyyyyyy_xyyz_0[j] = pb_x * tg_yyyyyyy_xyyz_0[j] + wp_x[j] * tg_yyyyyyy_xyyz_1[j] + 0.5 * fl1_fxn * tg_yyyyyyy_yyz_1[j];

                    tg_xyyyyyyy_xyzz_0[j] = pb_x * tg_yyyyyyy_xyzz_0[j] + wp_x[j] * tg_yyyyyyy_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyyy_yzz_1[j];

                    tg_xyyyyyyy_xzzz_0[j] = pb_x * tg_yyyyyyy_xzzz_0[j] + wp_x[j] * tg_yyyyyyy_xzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyyy_zzz_1[j];

                    tg_xyyyyyyy_yyyy_0[j] = pb_x * tg_yyyyyyy_yyyy_0[j] + wp_x[j] * tg_yyyyyyy_yyyy_1[j];

                    tg_xyyyyyyy_yyyz_0[j] = pb_x * tg_yyyyyyy_yyyz_0[j] + wp_x[j] * tg_yyyyyyy_yyyz_1[j];

                    tg_xyyyyyyy_yyzz_0[j] = pb_x * tg_yyyyyyy_yyzz_0[j] + wp_x[j] * tg_yyyyyyy_yyzz_1[j];

                    tg_xyyyyyyy_yzzz_0[j] = pb_x * tg_yyyyyyy_yzzz_0[j] + wp_x[j] * tg_yyyyyyy_yzzz_1[j];

                    tg_xyyyyyyy_zzzz_0[j] = pb_x * tg_yyyyyyy_zzzz_0[j] + wp_x[j] * tg_yyyyyyy_zzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSLSG_435_483(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (435,483)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {8, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {8, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_7_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_7_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
            {
                // set up pointers to Obara-Saika factors

                auto fxn = osFactors.data(4 * idx);

                // set up distances R(PB) = P - B

                auto pb_x = r_pb_x[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                // set up pointers to auxilary integrals

                auto tg_yyyyyyz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 435); 

                auto tg_yyyyyyz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 436); 

                auto tg_yyyyyyz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 437); 

                auto tg_yyyyyyz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 438); 

                auto tg_yyyyyyz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 439); 

                auto tg_yyyyyyz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 440); 

                auto tg_yyyyyyz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 441); 

                auto tg_yyyyyyz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 442); 

                auto tg_yyyyyyz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 443); 

                auto tg_yyyyyyz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 444); 

                auto tg_yyyyyyz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 445); 

                auto tg_yyyyyyz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 446); 

                auto tg_yyyyyyz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 447); 

                auto tg_yyyyyyz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 448); 

                auto tg_yyyyyyz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 449); 

                auto tg_yyyyyzz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 450); 

                auto tg_yyyyyzz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 451); 

                auto tg_yyyyyzz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 452); 

                auto tg_yyyyyzz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 453); 

                auto tg_yyyyyzz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 454); 

                auto tg_yyyyyzz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 455); 

                auto tg_yyyyyzz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 456); 

                auto tg_yyyyyzz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 457); 

                auto tg_yyyyyzz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 458); 

                auto tg_yyyyyzz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 459); 

                auto tg_yyyyyzz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 460); 

                auto tg_yyyyyzz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 461); 

                auto tg_yyyyyzz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 462); 

                auto tg_yyyyyzz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 463); 

                auto tg_yyyyyzz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 464); 

                auto tg_yyyyzzz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 465); 

                auto tg_yyyyzzz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 466); 

                auto tg_yyyyzzz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 467); 

                auto tg_yyyyzzz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 468); 

                auto tg_yyyyzzz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 469); 

                auto tg_yyyyzzz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 470); 

                auto tg_yyyyzzz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 471); 

                auto tg_yyyyzzz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 472); 

                auto tg_yyyyzzz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 473); 

                auto tg_yyyyzzz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 474); 

                auto tg_yyyyzzz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 475); 

                auto tg_yyyyzzz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 476); 

                auto tg_yyyyzzz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 477); 

                auto tg_yyyyzzz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 478); 

                auto tg_yyyyzzz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 479); 

                auto tg_yyyzzzz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 480); 

                auto tg_yyyzzzz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 481); 

                auto tg_yyyzzzz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 482); 

                auto tg_yyyyyyz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 435); 

                auto tg_yyyyyyz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 436); 

                auto tg_yyyyyyz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 437); 

                auto tg_yyyyyyz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 438); 

                auto tg_yyyyyyz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 439); 

                auto tg_yyyyyyz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 440); 

                auto tg_yyyyyyz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 441); 

                auto tg_yyyyyyz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 442); 

                auto tg_yyyyyyz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 443); 

                auto tg_yyyyyyz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 444); 

                auto tg_yyyyyyz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 445); 

                auto tg_yyyyyyz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 446); 

                auto tg_yyyyyyz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 447); 

                auto tg_yyyyyyz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 448); 

                auto tg_yyyyyyz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 449); 

                auto tg_yyyyyzz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 450); 

                auto tg_yyyyyzz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 451); 

                auto tg_yyyyyzz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 452); 

                auto tg_yyyyyzz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 453); 

                auto tg_yyyyyzz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 454); 

                auto tg_yyyyyzz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 455); 

                auto tg_yyyyyzz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 456); 

                auto tg_yyyyyzz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 457); 

                auto tg_yyyyyzz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 458); 

                auto tg_yyyyyzz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 459); 

                auto tg_yyyyyzz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 460); 

                auto tg_yyyyyzz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 461); 

                auto tg_yyyyyzz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 462); 

                auto tg_yyyyyzz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 463); 

                auto tg_yyyyyzz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 464); 

                auto tg_yyyyzzz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 465); 

                auto tg_yyyyzzz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 466); 

                auto tg_yyyyzzz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 467); 

                auto tg_yyyyzzz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 468); 

                auto tg_yyyyzzz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 469); 

                auto tg_yyyyzzz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 470); 

                auto tg_yyyyzzz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 471); 

                auto tg_yyyyzzz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 472); 

                auto tg_yyyyzzz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 473); 

                auto tg_yyyyzzz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 474); 

                auto tg_yyyyzzz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 475); 

                auto tg_yyyyzzz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 476); 

                auto tg_yyyyzzz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 477); 

                auto tg_yyyyzzz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 478); 

                auto tg_yyyyzzz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 479); 

                auto tg_yyyzzzz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 480); 

                auto tg_yyyzzzz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 481); 

                auto tg_yyyzzzz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 482); 

                auto tg_yyyyyyz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 290); 

                auto tg_yyyyyyz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 291); 

                auto tg_yyyyyyz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 292); 

                auto tg_yyyyyyz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 293); 

                auto tg_yyyyyyz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 294); 

                auto tg_yyyyyyz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 295); 

                auto tg_yyyyyyz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 296); 

                auto tg_yyyyyyz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 297); 

                auto tg_yyyyyyz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 298); 

                auto tg_yyyyyyz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 299); 

                auto tg_yyyyyzz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 300); 

                auto tg_yyyyyzz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 301); 

                auto tg_yyyyyzz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 302); 

                auto tg_yyyyyzz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 303); 

                auto tg_yyyyyzz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 304); 

                auto tg_yyyyyzz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 305); 

                auto tg_yyyyyzz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 306); 

                auto tg_yyyyyzz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 307); 

                auto tg_yyyyyzz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 308); 

                auto tg_yyyyyzz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 309); 

                auto tg_yyyyzzz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 310); 

                auto tg_yyyyzzz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 311); 

                auto tg_yyyyzzz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 312); 

                auto tg_yyyyzzz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 313); 

                auto tg_yyyyzzz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 314); 

                auto tg_yyyyzzz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 315); 

                auto tg_yyyyzzz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 316); 

                auto tg_yyyyzzz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 317); 

                auto tg_yyyyzzz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 318); 

                auto tg_yyyyzzz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 319); 

                auto tg_yyyzzzz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 320); 

                auto tg_yyyzzzz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 321); 

                auto tg_yyyzzzz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 322); 

                // set up pointers to integrals

                auto tg_xyyyyyyz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 435); 

                auto tg_xyyyyyyz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 436); 

                auto tg_xyyyyyyz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 437); 

                auto tg_xyyyyyyz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 438); 

                auto tg_xyyyyyyz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 439); 

                auto tg_xyyyyyyz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 440); 

                auto tg_xyyyyyyz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 441); 

                auto tg_xyyyyyyz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 442); 

                auto tg_xyyyyyyz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 443); 

                auto tg_xyyyyyyz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 444); 

                auto tg_xyyyyyyz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 445); 

                auto tg_xyyyyyyz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 446); 

                auto tg_xyyyyyyz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 447); 

                auto tg_xyyyyyyz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 448); 

                auto tg_xyyyyyyz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 449); 

                auto tg_xyyyyyzz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 450); 

                auto tg_xyyyyyzz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 451); 

                auto tg_xyyyyyzz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 452); 

                auto tg_xyyyyyzz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 453); 

                auto tg_xyyyyyzz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 454); 

                auto tg_xyyyyyzz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 455); 

                auto tg_xyyyyyzz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 456); 

                auto tg_xyyyyyzz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 457); 

                auto tg_xyyyyyzz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 458); 

                auto tg_xyyyyyzz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 459); 

                auto tg_xyyyyyzz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 460); 

                auto tg_xyyyyyzz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 461); 

                auto tg_xyyyyyzz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 462); 

                auto tg_xyyyyyzz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 463); 

                auto tg_xyyyyyzz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 464); 

                auto tg_xyyyyzzz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 465); 

                auto tg_xyyyyzzz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 466); 

                auto tg_xyyyyzzz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 467); 

                auto tg_xyyyyzzz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 468); 

                auto tg_xyyyyzzz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 469); 

                auto tg_xyyyyzzz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 470); 

                auto tg_xyyyyzzz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 471); 

                auto tg_xyyyyzzz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 472); 

                auto tg_xyyyyzzz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 473); 

                auto tg_xyyyyzzz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 474); 

                auto tg_xyyyyzzz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 475); 

                auto tg_xyyyyzzz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 476); 

                auto tg_xyyyyzzz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 477); 

                auto tg_xyyyyzzz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 478); 

                auto tg_xyyyyzzz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 479); 

                auto tg_xyyyzzzz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 480); 

                auto tg_xyyyzzzz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 481); 

                auto tg_xyyyzzzz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 482); 

                // Batch of Integrals (435,483)

                #pragma omp simd aligned(fxn, tg_xyyyyyyz_xxxx_0, tg_xyyyyyyz_xxxy_0, tg_xyyyyyyz_xxxz_0, \
                                         tg_xyyyyyyz_xxyy_0, tg_xyyyyyyz_xxyz_0, tg_xyyyyyyz_xxzz_0, tg_xyyyyyyz_xyyy_0, \
                                         tg_xyyyyyyz_xyyz_0, tg_xyyyyyyz_xyzz_0, tg_xyyyyyyz_xzzz_0, tg_xyyyyyyz_yyyy_0, \
                                         tg_xyyyyyyz_yyyz_0, tg_xyyyyyyz_yyzz_0, tg_xyyyyyyz_yzzz_0, tg_xyyyyyyz_zzzz_0, \
                                         tg_xyyyyyzz_xxxx_0, tg_xyyyyyzz_xxxy_0, tg_xyyyyyzz_xxxz_0, tg_xyyyyyzz_xxyy_0, \
                                         tg_xyyyyyzz_xxyz_0, tg_xyyyyyzz_xxzz_0, tg_xyyyyyzz_xyyy_0, tg_xyyyyyzz_xyyz_0, \
                                         tg_xyyyyyzz_xyzz_0, tg_xyyyyyzz_xzzz_0, tg_xyyyyyzz_yyyy_0, tg_xyyyyyzz_yyyz_0, \
                                         tg_xyyyyyzz_yyzz_0, tg_xyyyyyzz_yzzz_0, tg_xyyyyyzz_zzzz_0, tg_xyyyyzzz_xxxx_0, \
                                         tg_xyyyyzzz_xxxy_0, tg_xyyyyzzz_xxxz_0, tg_xyyyyzzz_xxyy_0, tg_xyyyyzzz_xxyz_0, \
                                         tg_xyyyyzzz_xxzz_0, tg_xyyyyzzz_xyyy_0, tg_xyyyyzzz_xyyz_0, tg_xyyyyzzz_xyzz_0, \
                                         tg_xyyyyzzz_xzzz_0, tg_xyyyyzzz_yyyy_0, tg_xyyyyzzz_yyyz_0, tg_xyyyyzzz_yyzz_0, \
                                         tg_xyyyyzzz_yzzz_0, tg_xyyyyzzz_zzzz_0, tg_xyyyzzzz_xxxx_0, tg_xyyyzzzz_xxxy_0, \
                                         tg_xyyyzzzz_xxxz_0, tg_yyyyyyz_xxx_1, tg_yyyyyyz_xxxx_0, tg_yyyyyyz_xxxx_1, \
                                         tg_yyyyyyz_xxxy_0, tg_yyyyyyz_xxxy_1, tg_yyyyyyz_xxxz_0, tg_yyyyyyz_xxxz_1, \
                                         tg_yyyyyyz_xxy_1, tg_yyyyyyz_xxyy_0, tg_yyyyyyz_xxyy_1, tg_yyyyyyz_xxyz_0, \
                                         tg_yyyyyyz_xxyz_1, tg_yyyyyyz_xxz_1, tg_yyyyyyz_xxzz_0, tg_yyyyyyz_xxzz_1, \
                                         tg_yyyyyyz_xyy_1, tg_yyyyyyz_xyyy_0, tg_yyyyyyz_xyyy_1, tg_yyyyyyz_xyyz_0, \
                                         tg_yyyyyyz_xyyz_1, tg_yyyyyyz_xyz_1, tg_yyyyyyz_xyzz_0, tg_yyyyyyz_xyzz_1, \
                                         tg_yyyyyyz_xzz_1, tg_yyyyyyz_xzzz_0, tg_yyyyyyz_xzzz_1, tg_yyyyyyz_yyy_1, \
                                         tg_yyyyyyz_yyyy_0, tg_yyyyyyz_yyyy_1, tg_yyyyyyz_yyyz_0, tg_yyyyyyz_yyyz_1, \
                                         tg_yyyyyyz_yyz_1, tg_yyyyyyz_yyzz_0, tg_yyyyyyz_yyzz_1, tg_yyyyyyz_yzz_1, \
                                         tg_yyyyyyz_yzzz_0, tg_yyyyyyz_yzzz_1, tg_yyyyyyz_zzz_1, tg_yyyyyyz_zzzz_0, \
                                         tg_yyyyyyz_zzzz_1, tg_yyyyyzz_xxx_1, tg_yyyyyzz_xxxx_0, tg_yyyyyzz_xxxx_1, \
                                         tg_yyyyyzz_xxxy_0, tg_yyyyyzz_xxxy_1, tg_yyyyyzz_xxxz_0, tg_yyyyyzz_xxxz_1, \
                                         tg_yyyyyzz_xxy_1, tg_yyyyyzz_xxyy_0, tg_yyyyyzz_xxyy_1, tg_yyyyyzz_xxyz_0, \
                                         tg_yyyyyzz_xxyz_1, tg_yyyyyzz_xxz_1, tg_yyyyyzz_xxzz_0, tg_yyyyyzz_xxzz_1, \
                                         tg_yyyyyzz_xyy_1, tg_yyyyyzz_xyyy_0, tg_yyyyyzz_xyyy_1, tg_yyyyyzz_xyyz_0, \
                                         tg_yyyyyzz_xyyz_1, tg_yyyyyzz_xyz_1, tg_yyyyyzz_xyzz_0, tg_yyyyyzz_xyzz_1, \
                                         tg_yyyyyzz_xzz_1, tg_yyyyyzz_xzzz_0, tg_yyyyyzz_xzzz_1, tg_yyyyyzz_yyy_1, \
                                         tg_yyyyyzz_yyyy_0, tg_yyyyyzz_yyyy_1, tg_yyyyyzz_yyyz_0, tg_yyyyyzz_yyyz_1, \
                                         tg_yyyyyzz_yyz_1, tg_yyyyyzz_yyzz_0, tg_yyyyyzz_yyzz_1, tg_yyyyyzz_yzz_1, \
                                         tg_yyyyyzz_yzzz_0, tg_yyyyyzz_yzzz_1, tg_yyyyyzz_zzz_1, tg_yyyyyzz_zzzz_0, \
                                         tg_yyyyyzz_zzzz_1, tg_yyyyzzz_xxx_1, tg_yyyyzzz_xxxx_0, tg_yyyyzzz_xxxx_1, \
                                         tg_yyyyzzz_xxxy_0, tg_yyyyzzz_xxxy_1, tg_yyyyzzz_xxxz_0, tg_yyyyzzz_xxxz_1, \
                                         tg_yyyyzzz_xxy_1, tg_yyyyzzz_xxyy_0, tg_yyyyzzz_xxyy_1, tg_yyyyzzz_xxyz_0, \
                                         tg_yyyyzzz_xxyz_1, tg_yyyyzzz_xxz_1, tg_yyyyzzz_xxzz_0, tg_yyyyzzz_xxzz_1, \
                                         tg_yyyyzzz_xyy_1, tg_yyyyzzz_xyyy_0, tg_yyyyzzz_xyyy_1, tg_yyyyzzz_xyyz_0, \
                                         tg_yyyyzzz_xyyz_1, tg_yyyyzzz_xyz_1, tg_yyyyzzz_xyzz_0, tg_yyyyzzz_xyzz_1, \
                                         tg_yyyyzzz_xzz_1, tg_yyyyzzz_xzzz_0, tg_yyyyzzz_xzzz_1, tg_yyyyzzz_yyy_1, \
                                         tg_yyyyzzz_yyyy_0, tg_yyyyzzz_yyyy_1, tg_yyyyzzz_yyyz_0, tg_yyyyzzz_yyyz_1, \
                                         tg_yyyyzzz_yyz_1, tg_yyyyzzz_yyzz_0, tg_yyyyzzz_yyzz_1, tg_yyyyzzz_yzz_1, \
                                         tg_yyyyzzz_yzzz_0, tg_yyyyzzz_yzzz_1, tg_yyyyzzz_zzz_1, tg_yyyyzzz_zzzz_0, \
                                         tg_yyyyzzz_zzzz_1, tg_yyyzzzz_xxx_1, tg_yyyzzzz_xxxx_0, tg_yyyzzzz_xxxx_1, \
                                         tg_yyyzzzz_xxxy_0, tg_yyyzzzz_xxxy_1, tg_yyyzzzz_xxxz_0, tg_yyyzzzz_xxxz_1, \
                                         tg_yyyzzzz_xxy_1, tg_yyyzzzz_xxz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fxn = fxn[j];

                    tg_xyyyyyyz_xxxx_0[j] = pb_x * tg_yyyyyyz_xxxx_0[j] + wp_x[j] * tg_yyyyyyz_xxxx_1[j] + 2.0 * fl1_fxn * tg_yyyyyyz_xxx_1[j];

                    tg_xyyyyyyz_xxxy_0[j] = pb_x * tg_yyyyyyz_xxxy_0[j] + wp_x[j] * tg_yyyyyyz_xxxy_1[j] + 1.5 * fl1_fxn * tg_yyyyyyz_xxy_1[j];

                    tg_xyyyyyyz_xxxz_0[j] = pb_x * tg_yyyyyyz_xxxz_0[j] + wp_x[j] * tg_yyyyyyz_xxxz_1[j] + 1.5 * fl1_fxn * tg_yyyyyyz_xxz_1[j];

                    tg_xyyyyyyz_xxyy_0[j] = pb_x * tg_yyyyyyz_xxyy_0[j] + wp_x[j] * tg_yyyyyyz_xxyy_1[j] + fl1_fxn * tg_yyyyyyz_xyy_1[j];

                    tg_xyyyyyyz_xxyz_0[j] = pb_x * tg_yyyyyyz_xxyz_0[j] + wp_x[j] * tg_yyyyyyz_xxyz_1[j] + fl1_fxn * tg_yyyyyyz_xyz_1[j];

                    tg_xyyyyyyz_xxzz_0[j] = pb_x * tg_yyyyyyz_xxzz_0[j] + wp_x[j] * tg_yyyyyyz_xxzz_1[j] + fl1_fxn * tg_yyyyyyz_xzz_1[j];

                    tg_xyyyyyyz_xyyy_0[j] = pb_x * tg_yyyyyyz_xyyy_0[j] + wp_x[j] * tg_yyyyyyz_xyyy_1[j] + 0.5 * fl1_fxn * tg_yyyyyyz_yyy_1[j];

                    tg_xyyyyyyz_xyyz_0[j] = pb_x * tg_yyyyyyz_xyyz_0[j] + wp_x[j] * tg_yyyyyyz_xyyz_1[j] + 0.5 * fl1_fxn * tg_yyyyyyz_yyz_1[j];

                    tg_xyyyyyyz_xyzz_0[j] = pb_x * tg_yyyyyyz_xyzz_0[j] + wp_x[j] * tg_yyyyyyz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyyz_yzz_1[j];

                    tg_xyyyyyyz_xzzz_0[j] = pb_x * tg_yyyyyyz_xzzz_0[j] + wp_x[j] * tg_yyyyyyz_xzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyyz_zzz_1[j];

                    tg_xyyyyyyz_yyyy_0[j] = pb_x * tg_yyyyyyz_yyyy_0[j] + wp_x[j] * tg_yyyyyyz_yyyy_1[j];

                    tg_xyyyyyyz_yyyz_0[j] = pb_x * tg_yyyyyyz_yyyz_0[j] + wp_x[j] * tg_yyyyyyz_yyyz_1[j];

                    tg_xyyyyyyz_yyzz_0[j] = pb_x * tg_yyyyyyz_yyzz_0[j] + wp_x[j] * tg_yyyyyyz_yyzz_1[j];

                    tg_xyyyyyyz_yzzz_0[j] = pb_x * tg_yyyyyyz_yzzz_0[j] + wp_x[j] * tg_yyyyyyz_yzzz_1[j];

                    tg_xyyyyyyz_zzzz_0[j] = pb_x * tg_yyyyyyz_zzzz_0[j] + wp_x[j] * tg_yyyyyyz_zzzz_1[j];

                    tg_xyyyyyzz_xxxx_0[j] = pb_x * tg_yyyyyzz_xxxx_0[j] + wp_x[j] * tg_yyyyyzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_yyyyyzz_xxx_1[j];

                    tg_xyyyyyzz_xxxy_0[j] = pb_x * tg_yyyyyzz_xxxy_0[j] + wp_x[j] * tg_yyyyyzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_yyyyyzz_xxy_1[j];

                    tg_xyyyyyzz_xxxz_0[j] = pb_x * tg_yyyyyzz_xxxz_0[j] + wp_x[j] * tg_yyyyyzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_yyyyyzz_xxz_1[j];

                    tg_xyyyyyzz_xxyy_0[j] = pb_x * tg_yyyyyzz_xxyy_0[j] + wp_x[j] * tg_yyyyyzz_xxyy_1[j] + fl1_fxn * tg_yyyyyzz_xyy_1[j];

                    tg_xyyyyyzz_xxyz_0[j] = pb_x * tg_yyyyyzz_xxyz_0[j] + wp_x[j] * tg_yyyyyzz_xxyz_1[j] + fl1_fxn * tg_yyyyyzz_xyz_1[j];

                    tg_xyyyyyzz_xxzz_0[j] = pb_x * tg_yyyyyzz_xxzz_0[j] + wp_x[j] * tg_yyyyyzz_xxzz_1[j] + fl1_fxn * tg_yyyyyzz_xzz_1[j];

                    tg_xyyyyyzz_xyyy_0[j] = pb_x * tg_yyyyyzz_xyyy_0[j] + wp_x[j] * tg_yyyyyzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_yyyyyzz_yyy_1[j];

                    tg_xyyyyyzz_xyyz_0[j] = pb_x * tg_yyyyyzz_xyyz_0[j] + wp_x[j] * tg_yyyyyzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_yyyyyzz_yyz_1[j];

                    tg_xyyyyyzz_xyzz_0[j] = pb_x * tg_yyyyyzz_xyzz_0[j] + wp_x[j] * tg_yyyyyzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyzz_yzz_1[j];

                    tg_xyyyyyzz_xzzz_0[j] = pb_x * tg_yyyyyzz_xzzz_0[j] + wp_x[j] * tg_yyyyyzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyzz_zzz_1[j];

                    tg_xyyyyyzz_yyyy_0[j] = pb_x * tg_yyyyyzz_yyyy_0[j] + wp_x[j] * tg_yyyyyzz_yyyy_1[j];

                    tg_xyyyyyzz_yyyz_0[j] = pb_x * tg_yyyyyzz_yyyz_0[j] + wp_x[j] * tg_yyyyyzz_yyyz_1[j];

                    tg_xyyyyyzz_yyzz_0[j] = pb_x * tg_yyyyyzz_yyzz_0[j] + wp_x[j] * tg_yyyyyzz_yyzz_1[j];

                    tg_xyyyyyzz_yzzz_0[j] = pb_x * tg_yyyyyzz_yzzz_0[j] + wp_x[j] * tg_yyyyyzz_yzzz_1[j];

                    tg_xyyyyyzz_zzzz_0[j] = pb_x * tg_yyyyyzz_zzzz_0[j] + wp_x[j] * tg_yyyyyzz_zzzz_1[j];

                    tg_xyyyyzzz_xxxx_0[j] = pb_x * tg_yyyyzzz_xxxx_0[j] + wp_x[j] * tg_yyyyzzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_yyyyzzz_xxx_1[j];

                    tg_xyyyyzzz_xxxy_0[j] = pb_x * tg_yyyyzzz_xxxy_0[j] + wp_x[j] * tg_yyyyzzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_yyyyzzz_xxy_1[j];

                    tg_xyyyyzzz_xxxz_0[j] = pb_x * tg_yyyyzzz_xxxz_0[j] + wp_x[j] * tg_yyyyzzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_yyyyzzz_xxz_1[j];

                    tg_xyyyyzzz_xxyy_0[j] = pb_x * tg_yyyyzzz_xxyy_0[j] + wp_x[j] * tg_yyyyzzz_xxyy_1[j] + fl1_fxn * tg_yyyyzzz_xyy_1[j];

                    tg_xyyyyzzz_xxyz_0[j] = pb_x * tg_yyyyzzz_xxyz_0[j] + wp_x[j] * tg_yyyyzzz_xxyz_1[j] + fl1_fxn * tg_yyyyzzz_xyz_1[j];

                    tg_xyyyyzzz_xxzz_0[j] = pb_x * tg_yyyyzzz_xxzz_0[j] + wp_x[j] * tg_yyyyzzz_xxzz_1[j] + fl1_fxn * tg_yyyyzzz_xzz_1[j];

                    tg_xyyyyzzz_xyyy_0[j] = pb_x * tg_yyyyzzz_xyyy_0[j] + wp_x[j] * tg_yyyyzzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_yyyyzzz_yyy_1[j];

                    tg_xyyyyzzz_xyyz_0[j] = pb_x * tg_yyyyzzz_xyyz_0[j] + wp_x[j] * tg_yyyyzzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_yyyyzzz_yyz_1[j];

                    tg_xyyyyzzz_xyzz_0[j] = pb_x * tg_yyyyzzz_xyzz_0[j] + wp_x[j] * tg_yyyyzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyyyzzz_yzz_1[j];

                    tg_xyyyyzzz_xzzz_0[j] = pb_x * tg_yyyyzzz_xzzz_0[j] + wp_x[j] * tg_yyyyzzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyzzz_zzz_1[j];

                    tg_xyyyyzzz_yyyy_0[j] = pb_x * tg_yyyyzzz_yyyy_0[j] + wp_x[j] * tg_yyyyzzz_yyyy_1[j];

                    tg_xyyyyzzz_yyyz_0[j] = pb_x * tg_yyyyzzz_yyyz_0[j] + wp_x[j] * tg_yyyyzzz_yyyz_1[j];

                    tg_xyyyyzzz_yyzz_0[j] = pb_x * tg_yyyyzzz_yyzz_0[j] + wp_x[j] * tg_yyyyzzz_yyzz_1[j];

                    tg_xyyyyzzz_yzzz_0[j] = pb_x * tg_yyyyzzz_yzzz_0[j] + wp_x[j] * tg_yyyyzzz_yzzz_1[j];

                    tg_xyyyyzzz_zzzz_0[j] = pb_x * tg_yyyyzzz_zzzz_0[j] + wp_x[j] * tg_yyyyzzz_zzzz_1[j];

                    tg_xyyyzzzz_xxxx_0[j] = pb_x * tg_yyyzzzz_xxxx_0[j] + wp_x[j] * tg_yyyzzzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_yyyzzzz_xxx_1[j];

                    tg_xyyyzzzz_xxxy_0[j] = pb_x * tg_yyyzzzz_xxxy_0[j] + wp_x[j] * tg_yyyzzzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_yyyzzzz_xxy_1[j];

                    tg_xyyyzzzz_xxxz_0[j] = pb_x * tg_yyyzzzz_xxxz_0[j] + wp_x[j] * tg_yyyzzzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_yyyzzzz_xxz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSLSG_483_531(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (483,531)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {8, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {8, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_7_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_7_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
            {
                // set up pointers to Obara-Saika factors

                auto fxn = osFactors.data(4 * idx);

                // set up distances R(PB) = P - B

                auto pb_x = r_pb_x[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                // set up pointers to auxilary integrals

                auto tg_yyyzzzz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 483); 

                auto tg_yyyzzzz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 484); 

                auto tg_yyyzzzz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 485); 

                auto tg_yyyzzzz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 486); 

                auto tg_yyyzzzz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 487); 

                auto tg_yyyzzzz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 488); 

                auto tg_yyyzzzz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 489); 

                auto tg_yyyzzzz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 490); 

                auto tg_yyyzzzz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 491); 

                auto tg_yyyzzzz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 492); 

                auto tg_yyyzzzz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 493); 

                auto tg_yyyzzzz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 494); 

                auto tg_yyzzzzz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 495); 

                auto tg_yyzzzzz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 496); 

                auto tg_yyzzzzz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 497); 

                auto tg_yyzzzzz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 498); 

                auto tg_yyzzzzz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 499); 

                auto tg_yyzzzzz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 500); 

                auto tg_yyzzzzz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 501); 

                auto tg_yyzzzzz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 502); 

                auto tg_yyzzzzz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 503); 

                auto tg_yyzzzzz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 504); 

                auto tg_yyzzzzz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 505); 

                auto tg_yyzzzzz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 506); 

                auto tg_yyzzzzz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 507); 

                auto tg_yyzzzzz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 508); 

                auto tg_yyzzzzz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 509); 

                auto tg_yzzzzzz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 510); 

                auto tg_yzzzzzz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 511); 

                auto tg_yzzzzzz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 512); 

                auto tg_yzzzzzz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 513); 

                auto tg_yzzzzzz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 514); 

                auto tg_yzzzzzz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 515); 

                auto tg_yzzzzzz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 516); 

                auto tg_yzzzzzz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 517); 

                auto tg_yzzzzzz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 518); 

                auto tg_yzzzzzz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 519); 

                auto tg_yzzzzzz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 520); 

                auto tg_yzzzzzz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 521); 

                auto tg_yzzzzzz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 522); 

                auto tg_yzzzzzz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 523); 

                auto tg_yzzzzzz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 524); 

                auto tg_zzzzzzz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 525); 

                auto tg_zzzzzzz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 526); 

                auto tg_zzzzzzz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 527); 

                auto tg_zzzzzzz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 528); 

                auto tg_zzzzzzz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 529); 

                auto tg_zzzzzzz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 530); 

                auto tg_yyyzzzz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 483); 

                auto tg_yyyzzzz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 484); 

                auto tg_yyyzzzz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 485); 

                auto tg_yyyzzzz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 486); 

                auto tg_yyyzzzz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 487); 

                auto tg_yyyzzzz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 488); 

                auto tg_yyyzzzz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 489); 

                auto tg_yyyzzzz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 490); 

                auto tg_yyyzzzz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 491); 

                auto tg_yyyzzzz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 492); 

                auto tg_yyyzzzz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 493); 

                auto tg_yyyzzzz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 494); 

                auto tg_yyzzzzz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 495); 

                auto tg_yyzzzzz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 496); 

                auto tg_yyzzzzz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 497); 

                auto tg_yyzzzzz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 498); 

                auto tg_yyzzzzz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 499); 

                auto tg_yyzzzzz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 500); 

                auto tg_yyzzzzz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 501); 

                auto tg_yyzzzzz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 502); 

                auto tg_yyzzzzz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 503); 

                auto tg_yyzzzzz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 504); 

                auto tg_yyzzzzz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 505); 

                auto tg_yyzzzzz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 506); 

                auto tg_yyzzzzz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 507); 

                auto tg_yyzzzzz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 508); 

                auto tg_yyzzzzz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 509); 

                auto tg_yzzzzzz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 510); 

                auto tg_yzzzzzz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 511); 

                auto tg_yzzzzzz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 512); 

                auto tg_yzzzzzz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 513); 

                auto tg_yzzzzzz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 514); 

                auto tg_yzzzzzz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 515); 

                auto tg_yzzzzzz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 516); 

                auto tg_yzzzzzz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 517); 

                auto tg_yzzzzzz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 518); 

                auto tg_yzzzzzz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 519); 

                auto tg_yzzzzzz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 520); 

                auto tg_yzzzzzz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 521); 

                auto tg_yzzzzzz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 522); 

                auto tg_yzzzzzz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 523); 

                auto tg_yzzzzzz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 524); 

                auto tg_zzzzzzz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 525); 

                auto tg_zzzzzzz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 526); 

                auto tg_zzzzzzz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 527); 

                auto tg_zzzzzzz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 528); 

                auto tg_zzzzzzz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 529); 

                auto tg_zzzzzzz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 530); 

                auto tg_yyyzzzz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 323); 

                auto tg_yyyzzzz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 324); 

                auto tg_yyyzzzz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 325); 

                auto tg_yyyzzzz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 326); 

                auto tg_yyyzzzz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 327); 

                auto tg_yyyzzzz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 328); 

                auto tg_yyyzzzz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 329); 

                auto tg_yyzzzzz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 330); 

                auto tg_yyzzzzz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 331); 

                auto tg_yyzzzzz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 332); 

                auto tg_yyzzzzz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 333); 

                auto tg_yyzzzzz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 334); 

                auto tg_yyzzzzz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 335); 

                auto tg_yyzzzzz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 336); 

                auto tg_yyzzzzz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 337); 

                auto tg_yyzzzzz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 338); 

                auto tg_yyzzzzz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 339); 

                auto tg_yzzzzzz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 340); 

                auto tg_yzzzzzz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 341); 

                auto tg_yzzzzzz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 342); 

                auto tg_yzzzzzz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 343); 

                auto tg_yzzzzzz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 344); 

                auto tg_yzzzzzz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 345); 

                auto tg_yzzzzzz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 346); 

                auto tg_yzzzzzz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 347); 

                auto tg_yzzzzzz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 348); 

                auto tg_yzzzzzz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 349); 

                auto tg_zzzzzzz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 350); 

                auto tg_zzzzzzz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 351); 

                auto tg_zzzzzzz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 352); 

                auto tg_zzzzzzz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 353); 

                auto tg_zzzzzzz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 354); 

                auto tg_zzzzzzz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 355); 

                // set up pointers to integrals

                auto tg_xyyyzzzz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 483); 

                auto tg_xyyyzzzz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 484); 

                auto tg_xyyyzzzz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 485); 

                auto tg_xyyyzzzz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 486); 

                auto tg_xyyyzzzz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 487); 

                auto tg_xyyyzzzz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 488); 

                auto tg_xyyyzzzz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 489); 

                auto tg_xyyyzzzz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 490); 

                auto tg_xyyyzzzz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 491); 

                auto tg_xyyyzzzz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 492); 

                auto tg_xyyyzzzz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 493); 

                auto tg_xyyyzzzz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 494); 

                auto tg_xyyzzzzz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 495); 

                auto tg_xyyzzzzz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 496); 

                auto tg_xyyzzzzz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 497); 

                auto tg_xyyzzzzz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 498); 

                auto tg_xyyzzzzz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 499); 

                auto tg_xyyzzzzz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 500); 

                auto tg_xyyzzzzz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 501); 

                auto tg_xyyzzzzz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 502); 

                auto tg_xyyzzzzz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 503); 

                auto tg_xyyzzzzz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 504); 

                auto tg_xyyzzzzz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 505); 

                auto tg_xyyzzzzz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 506); 

                auto tg_xyyzzzzz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 507); 

                auto tg_xyyzzzzz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 508); 

                auto tg_xyyzzzzz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 509); 

                auto tg_xyzzzzzz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 510); 

                auto tg_xyzzzzzz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 511); 

                auto tg_xyzzzzzz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 512); 

                auto tg_xyzzzzzz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 513); 

                auto tg_xyzzzzzz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 514); 

                auto tg_xyzzzzzz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 515); 

                auto tg_xyzzzzzz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 516); 

                auto tg_xyzzzzzz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 517); 

                auto tg_xyzzzzzz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 518); 

                auto tg_xyzzzzzz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 519); 

                auto tg_xyzzzzzz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 520); 

                auto tg_xyzzzzzz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 521); 

                auto tg_xyzzzzzz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 522); 

                auto tg_xyzzzzzz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 523); 

                auto tg_xyzzzzzz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 524); 

                auto tg_xzzzzzzz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 525); 

                auto tg_xzzzzzzz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 526); 

                auto tg_xzzzzzzz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 527); 

                auto tg_xzzzzzzz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 528); 

                auto tg_xzzzzzzz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 529); 

                auto tg_xzzzzzzz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 530); 

                // Batch of Integrals (483,531)

                #pragma omp simd aligned(fxn, tg_xyyyzzzz_xxyy_0, tg_xyyyzzzz_xxyz_0, tg_xyyyzzzz_xxzz_0, \
                                         tg_xyyyzzzz_xyyy_0, tg_xyyyzzzz_xyyz_0, tg_xyyyzzzz_xyzz_0, tg_xyyyzzzz_xzzz_0, \
                                         tg_xyyyzzzz_yyyy_0, tg_xyyyzzzz_yyyz_0, tg_xyyyzzzz_yyzz_0, tg_xyyyzzzz_yzzz_0, \
                                         tg_xyyyzzzz_zzzz_0, tg_xyyzzzzz_xxxx_0, tg_xyyzzzzz_xxxy_0, tg_xyyzzzzz_xxxz_0, \
                                         tg_xyyzzzzz_xxyy_0, tg_xyyzzzzz_xxyz_0, tg_xyyzzzzz_xxzz_0, tg_xyyzzzzz_xyyy_0, \
                                         tg_xyyzzzzz_xyyz_0, tg_xyyzzzzz_xyzz_0, tg_xyyzzzzz_xzzz_0, tg_xyyzzzzz_yyyy_0, \
                                         tg_xyyzzzzz_yyyz_0, tg_xyyzzzzz_yyzz_0, tg_xyyzzzzz_yzzz_0, tg_xyyzzzzz_zzzz_0, \
                                         tg_xyzzzzzz_xxxx_0, tg_xyzzzzzz_xxxy_0, tg_xyzzzzzz_xxxz_0, tg_xyzzzzzz_xxyy_0, \
                                         tg_xyzzzzzz_xxyz_0, tg_xyzzzzzz_xxzz_0, tg_xyzzzzzz_xyyy_0, tg_xyzzzzzz_xyyz_0, \
                                         tg_xyzzzzzz_xyzz_0, tg_xyzzzzzz_xzzz_0, tg_xyzzzzzz_yyyy_0, tg_xyzzzzzz_yyyz_0, \
                                         tg_xyzzzzzz_yyzz_0, tg_xyzzzzzz_yzzz_0, tg_xyzzzzzz_zzzz_0, tg_xzzzzzzz_xxxx_0, \
                                         tg_xzzzzzzz_xxxy_0, tg_xzzzzzzz_xxxz_0, tg_xzzzzzzz_xxyy_0, tg_xzzzzzzz_xxyz_0, \
                                         tg_xzzzzzzz_xxzz_0, tg_yyyzzzz_xxyy_0, tg_yyyzzzz_xxyy_1, tg_yyyzzzz_xxyz_0, \
                                         tg_yyyzzzz_xxyz_1, tg_yyyzzzz_xxzz_0, tg_yyyzzzz_xxzz_1, tg_yyyzzzz_xyy_1, \
                                         tg_yyyzzzz_xyyy_0, tg_yyyzzzz_xyyy_1, tg_yyyzzzz_xyyz_0, tg_yyyzzzz_xyyz_1, \
                                         tg_yyyzzzz_xyz_1, tg_yyyzzzz_xyzz_0, tg_yyyzzzz_xyzz_1, tg_yyyzzzz_xzz_1, \
                                         tg_yyyzzzz_xzzz_0, tg_yyyzzzz_xzzz_1, tg_yyyzzzz_yyy_1, tg_yyyzzzz_yyyy_0, \
                                         tg_yyyzzzz_yyyy_1, tg_yyyzzzz_yyyz_0, tg_yyyzzzz_yyyz_1, tg_yyyzzzz_yyz_1, \
                                         tg_yyyzzzz_yyzz_0, tg_yyyzzzz_yyzz_1, tg_yyyzzzz_yzz_1, tg_yyyzzzz_yzzz_0, \
                                         tg_yyyzzzz_yzzz_1, tg_yyyzzzz_zzz_1, tg_yyyzzzz_zzzz_0, tg_yyyzzzz_zzzz_1, \
                                         tg_yyzzzzz_xxx_1, tg_yyzzzzz_xxxx_0, tg_yyzzzzz_xxxx_1, tg_yyzzzzz_xxxy_0, \
                                         tg_yyzzzzz_xxxy_1, tg_yyzzzzz_xxxz_0, tg_yyzzzzz_xxxz_1, tg_yyzzzzz_xxy_1, \
                                         tg_yyzzzzz_xxyy_0, tg_yyzzzzz_xxyy_1, tg_yyzzzzz_xxyz_0, tg_yyzzzzz_xxyz_1, \
                                         tg_yyzzzzz_xxz_1, tg_yyzzzzz_xxzz_0, tg_yyzzzzz_xxzz_1, tg_yyzzzzz_xyy_1, \
                                         tg_yyzzzzz_xyyy_0, tg_yyzzzzz_xyyy_1, tg_yyzzzzz_xyyz_0, tg_yyzzzzz_xyyz_1, \
                                         tg_yyzzzzz_xyz_1, tg_yyzzzzz_xyzz_0, tg_yyzzzzz_xyzz_1, tg_yyzzzzz_xzz_1, \
                                         tg_yyzzzzz_xzzz_0, tg_yyzzzzz_xzzz_1, tg_yyzzzzz_yyy_1, tg_yyzzzzz_yyyy_0, \
                                         tg_yyzzzzz_yyyy_1, tg_yyzzzzz_yyyz_0, tg_yyzzzzz_yyyz_1, tg_yyzzzzz_yyz_1, \
                                         tg_yyzzzzz_yyzz_0, tg_yyzzzzz_yyzz_1, tg_yyzzzzz_yzz_1, tg_yyzzzzz_yzzz_0, \
                                         tg_yyzzzzz_yzzz_1, tg_yyzzzzz_zzz_1, tg_yyzzzzz_zzzz_0, tg_yyzzzzz_zzzz_1, \
                                         tg_yzzzzzz_xxx_1, tg_yzzzzzz_xxxx_0, tg_yzzzzzz_xxxx_1, tg_yzzzzzz_xxxy_0, \
                                         tg_yzzzzzz_xxxy_1, tg_yzzzzzz_xxxz_0, tg_yzzzzzz_xxxz_1, tg_yzzzzzz_xxy_1, \
                                         tg_yzzzzzz_xxyy_0, tg_yzzzzzz_xxyy_1, tg_yzzzzzz_xxyz_0, tg_yzzzzzz_xxyz_1, \
                                         tg_yzzzzzz_xxz_1, tg_yzzzzzz_xxzz_0, tg_yzzzzzz_xxzz_1, tg_yzzzzzz_xyy_1, \
                                         tg_yzzzzzz_xyyy_0, tg_yzzzzzz_xyyy_1, tg_yzzzzzz_xyyz_0, tg_yzzzzzz_xyyz_1, \
                                         tg_yzzzzzz_xyz_1, tg_yzzzzzz_xyzz_0, tg_yzzzzzz_xyzz_1, tg_yzzzzzz_xzz_1, \
                                         tg_yzzzzzz_xzzz_0, tg_yzzzzzz_xzzz_1, tg_yzzzzzz_yyy_1, tg_yzzzzzz_yyyy_0, \
                                         tg_yzzzzzz_yyyy_1, tg_yzzzzzz_yyyz_0, tg_yzzzzzz_yyyz_1, tg_yzzzzzz_yyz_1, \
                                         tg_yzzzzzz_yyzz_0, tg_yzzzzzz_yyzz_1, tg_yzzzzzz_yzz_1, tg_yzzzzzz_yzzz_0, \
                                         tg_yzzzzzz_yzzz_1, tg_yzzzzzz_zzz_1, tg_yzzzzzz_zzzz_0, tg_yzzzzzz_zzzz_1, \
                                         tg_zzzzzzz_xxx_1, tg_zzzzzzz_xxxx_0, tg_zzzzzzz_xxxx_1, tg_zzzzzzz_xxxy_0, \
                                         tg_zzzzzzz_xxxy_1, tg_zzzzzzz_xxxz_0, tg_zzzzzzz_xxxz_1, tg_zzzzzzz_xxy_1, \
                                         tg_zzzzzzz_xxyy_0, tg_zzzzzzz_xxyy_1, tg_zzzzzzz_xxyz_0, tg_zzzzzzz_xxyz_1, \
                                         tg_zzzzzzz_xxz_1, tg_zzzzzzz_xxzz_0, tg_zzzzzzz_xxzz_1, tg_zzzzzzz_xyy_1, \
                                         tg_zzzzzzz_xyz_1, tg_zzzzzzz_xzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fxn = fxn[j];

                    tg_xyyyzzzz_xxyy_0[j] = pb_x * tg_yyyzzzz_xxyy_0[j] + wp_x[j] * tg_yyyzzzz_xxyy_1[j] + fl1_fxn * tg_yyyzzzz_xyy_1[j];

                    tg_xyyyzzzz_xxyz_0[j] = pb_x * tg_yyyzzzz_xxyz_0[j] + wp_x[j] * tg_yyyzzzz_xxyz_1[j] + fl1_fxn * tg_yyyzzzz_xyz_1[j];

                    tg_xyyyzzzz_xxzz_0[j] = pb_x * tg_yyyzzzz_xxzz_0[j] + wp_x[j] * tg_yyyzzzz_xxzz_1[j] + fl1_fxn * tg_yyyzzzz_xzz_1[j];

                    tg_xyyyzzzz_xyyy_0[j] = pb_x * tg_yyyzzzz_xyyy_0[j] + wp_x[j] * tg_yyyzzzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_yyyzzzz_yyy_1[j];

                    tg_xyyyzzzz_xyyz_0[j] = pb_x * tg_yyyzzzz_xyyz_0[j] + wp_x[j] * tg_yyyzzzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_yyyzzzz_yyz_1[j];

                    tg_xyyyzzzz_xyzz_0[j] = pb_x * tg_yyyzzzz_xyzz_0[j] + wp_x[j] * tg_yyyzzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyyzzzz_yzz_1[j];

                    tg_xyyyzzzz_xzzz_0[j] = pb_x * tg_yyyzzzz_xzzz_0[j] + wp_x[j] * tg_yyyzzzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_yyyzzzz_zzz_1[j];

                    tg_xyyyzzzz_yyyy_0[j] = pb_x * tg_yyyzzzz_yyyy_0[j] + wp_x[j] * tg_yyyzzzz_yyyy_1[j];

                    tg_xyyyzzzz_yyyz_0[j] = pb_x * tg_yyyzzzz_yyyz_0[j] + wp_x[j] * tg_yyyzzzz_yyyz_1[j];

                    tg_xyyyzzzz_yyzz_0[j] = pb_x * tg_yyyzzzz_yyzz_0[j] + wp_x[j] * tg_yyyzzzz_yyzz_1[j];

                    tg_xyyyzzzz_yzzz_0[j] = pb_x * tg_yyyzzzz_yzzz_0[j] + wp_x[j] * tg_yyyzzzz_yzzz_1[j];

                    tg_xyyyzzzz_zzzz_0[j] = pb_x * tg_yyyzzzz_zzzz_0[j] + wp_x[j] * tg_yyyzzzz_zzzz_1[j];

                    tg_xyyzzzzz_xxxx_0[j] = pb_x * tg_yyzzzzz_xxxx_0[j] + wp_x[j] * tg_yyzzzzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_yyzzzzz_xxx_1[j];

                    tg_xyyzzzzz_xxxy_0[j] = pb_x * tg_yyzzzzz_xxxy_0[j] + wp_x[j] * tg_yyzzzzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_yyzzzzz_xxy_1[j];

                    tg_xyyzzzzz_xxxz_0[j] = pb_x * tg_yyzzzzz_xxxz_0[j] + wp_x[j] * tg_yyzzzzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_yyzzzzz_xxz_1[j];

                    tg_xyyzzzzz_xxyy_0[j] = pb_x * tg_yyzzzzz_xxyy_0[j] + wp_x[j] * tg_yyzzzzz_xxyy_1[j] + fl1_fxn * tg_yyzzzzz_xyy_1[j];

                    tg_xyyzzzzz_xxyz_0[j] = pb_x * tg_yyzzzzz_xxyz_0[j] + wp_x[j] * tg_yyzzzzz_xxyz_1[j] + fl1_fxn * tg_yyzzzzz_xyz_1[j];

                    tg_xyyzzzzz_xxzz_0[j] = pb_x * tg_yyzzzzz_xxzz_0[j] + wp_x[j] * tg_yyzzzzz_xxzz_1[j] + fl1_fxn * tg_yyzzzzz_xzz_1[j];

                    tg_xyyzzzzz_xyyy_0[j] = pb_x * tg_yyzzzzz_xyyy_0[j] + wp_x[j] * tg_yyzzzzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_yyzzzzz_yyy_1[j];

                    tg_xyyzzzzz_xyyz_0[j] = pb_x * tg_yyzzzzz_xyyz_0[j] + wp_x[j] * tg_yyzzzzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_yyzzzzz_yyz_1[j];

                    tg_xyyzzzzz_xyzz_0[j] = pb_x * tg_yyzzzzz_xyzz_0[j] + wp_x[j] * tg_yyzzzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyzzzzz_yzz_1[j];

                    tg_xyyzzzzz_xzzz_0[j] = pb_x * tg_yyzzzzz_xzzz_0[j] + wp_x[j] * tg_yyzzzzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_yyzzzzz_zzz_1[j];

                    tg_xyyzzzzz_yyyy_0[j] = pb_x * tg_yyzzzzz_yyyy_0[j] + wp_x[j] * tg_yyzzzzz_yyyy_1[j];

                    tg_xyyzzzzz_yyyz_0[j] = pb_x * tg_yyzzzzz_yyyz_0[j] + wp_x[j] * tg_yyzzzzz_yyyz_1[j];

                    tg_xyyzzzzz_yyzz_0[j] = pb_x * tg_yyzzzzz_yyzz_0[j] + wp_x[j] * tg_yyzzzzz_yyzz_1[j];

                    tg_xyyzzzzz_yzzz_0[j] = pb_x * tg_yyzzzzz_yzzz_0[j] + wp_x[j] * tg_yyzzzzz_yzzz_1[j];

                    tg_xyyzzzzz_zzzz_0[j] = pb_x * tg_yyzzzzz_zzzz_0[j] + wp_x[j] * tg_yyzzzzz_zzzz_1[j];

                    tg_xyzzzzzz_xxxx_0[j] = pb_x * tg_yzzzzzz_xxxx_0[j] + wp_x[j] * tg_yzzzzzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_yzzzzzz_xxx_1[j];

                    tg_xyzzzzzz_xxxy_0[j] = pb_x * tg_yzzzzzz_xxxy_0[j] + wp_x[j] * tg_yzzzzzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_yzzzzzz_xxy_1[j];

                    tg_xyzzzzzz_xxxz_0[j] = pb_x * tg_yzzzzzz_xxxz_0[j] + wp_x[j] * tg_yzzzzzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_yzzzzzz_xxz_1[j];

                    tg_xyzzzzzz_xxyy_0[j] = pb_x * tg_yzzzzzz_xxyy_0[j] + wp_x[j] * tg_yzzzzzz_xxyy_1[j] + fl1_fxn * tg_yzzzzzz_xyy_1[j];

                    tg_xyzzzzzz_xxyz_0[j] = pb_x * tg_yzzzzzz_xxyz_0[j] + wp_x[j] * tg_yzzzzzz_xxyz_1[j] + fl1_fxn * tg_yzzzzzz_xyz_1[j];

                    tg_xyzzzzzz_xxzz_0[j] = pb_x * tg_yzzzzzz_xxzz_0[j] + wp_x[j] * tg_yzzzzzz_xxzz_1[j] + fl1_fxn * tg_yzzzzzz_xzz_1[j];

                    tg_xyzzzzzz_xyyy_0[j] = pb_x * tg_yzzzzzz_xyyy_0[j] + wp_x[j] * tg_yzzzzzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_yzzzzzz_yyy_1[j];

                    tg_xyzzzzzz_xyyz_0[j] = pb_x * tg_yzzzzzz_xyyz_0[j] + wp_x[j] * tg_yzzzzzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_yzzzzzz_yyz_1[j];

                    tg_xyzzzzzz_xyzz_0[j] = pb_x * tg_yzzzzzz_xyzz_0[j] + wp_x[j] * tg_yzzzzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yzzzzzz_yzz_1[j];

                    tg_xyzzzzzz_xzzz_0[j] = pb_x * tg_yzzzzzz_xzzz_0[j] + wp_x[j] * tg_yzzzzzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_yzzzzzz_zzz_1[j];

                    tg_xyzzzzzz_yyyy_0[j] = pb_x * tg_yzzzzzz_yyyy_0[j] + wp_x[j] * tg_yzzzzzz_yyyy_1[j];

                    tg_xyzzzzzz_yyyz_0[j] = pb_x * tg_yzzzzzz_yyyz_0[j] + wp_x[j] * tg_yzzzzzz_yyyz_1[j];

                    tg_xyzzzzzz_yyzz_0[j] = pb_x * tg_yzzzzzz_yyzz_0[j] + wp_x[j] * tg_yzzzzzz_yyzz_1[j];

                    tg_xyzzzzzz_yzzz_0[j] = pb_x * tg_yzzzzzz_yzzz_0[j] + wp_x[j] * tg_yzzzzzz_yzzz_1[j];

                    tg_xyzzzzzz_zzzz_0[j] = pb_x * tg_yzzzzzz_zzzz_0[j] + wp_x[j] * tg_yzzzzzz_zzzz_1[j];

                    tg_xzzzzzzz_xxxx_0[j] = pb_x * tg_zzzzzzz_xxxx_0[j] + wp_x[j] * tg_zzzzzzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_zzzzzzz_xxx_1[j];

                    tg_xzzzzzzz_xxxy_0[j] = pb_x * tg_zzzzzzz_xxxy_0[j] + wp_x[j] * tg_zzzzzzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_zzzzzzz_xxy_1[j];

                    tg_xzzzzzzz_xxxz_0[j] = pb_x * tg_zzzzzzz_xxxz_0[j] + wp_x[j] * tg_zzzzzzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_zzzzzzz_xxz_1[j];

                    tg_xzzzzzzz_xxyy_0[j] = pb_x * tg_zzzzzzz_xxyy_0[j] + wp_x[j] * tg_zzzzzzz_xxyy_1[j] + fl1_fxn * tg_zzzzzzz_xyy_1[j];

                    tg_xzzzzzzz_xxyz_0[j] = pb_x * tg_zzzzzzz_xxyz_0[j] + wp_x[j] * tg_zzzzzzz_xxyz_1[j] + fl1_fxn * tg_zzzzzzz_xyz_1[j];

                    tg_xzzzzzzz_xxzz_0[j] = pb_x * tg_zzzzzzz_xxzz_0[j] + wp_x[j] * tg_zzzzzzz_xxzz_1[j] + fl1_fxn * tg_zzzzzzz_xzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSLSG_531_579(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (531,579)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {8, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {8, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_7_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_7_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_6_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_6_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
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

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                auto wp_y = wpDistances.data(3 * idx + 1);

                // set up pointers to auxilary integrals

                auto tg_yyyyyyy_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 420); 

                auto tg_yyyyyyy_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 421); 

                auto tg_yyyyyyy_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 422); 

                auto tg_yyyyyyy_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 423); 

                auto tg_yyyyyyy_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 424); 

                auto tg_yyyyyyy_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 425); 

                auto tg_yyyyyyy_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 426); 

                auto tg_yyyyyyy_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 427); 

                auto tg_yyyyyyy_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 428); 

                auto tg_yyyyyyy_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 429); 

                auto tg_yyyyyyy_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 430); 

                auto tg_yyyyyyy_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 431); 

                auto tg_yyyyyyy_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 432); 

                auto tg_yyyyyyy_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 433); 

                auto tg_yyyyyyy_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 434); 

                auto tg_yyyyyyz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 435); 

                auto tg_yyyyyyz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 436); 

                auto tg_yyyyyyz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 437); 

                auto tg_yyyyyyz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 438); 

                auto tg_yyyyyyz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 439); 

                auto tg_yyyyyyz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 440); 

                auto tg_yyyyyyz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 441); 

                auto tg_yyyyyyz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 442); 

                auto tg_yyyyyyz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 443); 

                auto tg_yyyyyyz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 444); 

                auto tg_yyyyyyz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 445); 

                auto tg_yyyyyyz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 446); 

                auto tg_yyyyyyz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 447); 

                auto tg_yyyyyyz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 448); 

                auto tg_yyyyyyz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 449); 

                auto tg_yyyyyzz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 450); 

                auto tg_yyyyyzz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 451); 

                auto tg_yyyyyzz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 452); 

                auto tg_yyyyyzz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 453); 

                auto tg_yyyyyzz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 454); 

                auto tg_yyyyyzz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 455); 

                auto tg_yyyyyzz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 456); 

                auto tg_yyyyyzz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 457); 

                auto tg_yyyyyzz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 458); 

                auto tg_zzzzzzz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 531); 

                auto tg_zzzzzzz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 532); 

                auto tg_zzzzzzz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 533); 

                auto tg_zzzzzzz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 534); 

                auto tg_zzzzzzz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 535); 

                auto tg_zzzzzzz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 536); 

                auto tg_zzzzzzz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 537); 

                auto tg_zzzzzzz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 538); 

                auto tg_zzzzzzz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 539); 

                auto tg_yyyyyyy_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 420); 

                auto tg_yyyyyyy_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 421); 

                auto tg_yyyyyyy_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 422); 

                auto tg_yyyyyyy_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 423); 

                auto tg_yyyyyyy_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 424); 

                auto tg_yyyyyyy_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 425); 

                auto tg_yyyyyyy_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 426); 

                auto tg_yyyyyyy_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 427); 

                auto tg_yyyyyyy_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 428); 

                auto tg_yyyyyyy_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 429); 

                auto tg_yyyyyyy_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 430); 

                auto tg_yyyyyyy_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 431); 

                auto tg_yyyyyyy_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 432); 

                auto tg_yyyyyyy_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 433); 

                auto tg_yyyyyyy_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 434); 

                auto tg_yyyyyyz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 435); 

                auto tg_yyyyyyz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 436); 

                auto tg_yyyyyyz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 437); 

                auto tg_yyyyyyz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 438); 

                auto tg_yyyyyyz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 439); 

                auto tg_yyyyyyz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 440); 

                auto tg_yyyyyyz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 441); 

                auto tg_yyyyyyz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 442); 

                auto tg_yyyyyyz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 443); 

                auto tg_yyyyyyz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 444); 

                auto tg_yyyyyyz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 445); 

                auto tg_yyyyyyz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 446); 

                auto tg_yyyyyyz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 447); 

                auto tg_yyyyyyz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 448); 

                auto tg_yyyyyyz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 449); 

                auto tg_yyyyyzz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 450); 

                auto tg_yyyyyzz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 451); 

                auto tg_yyyyyzz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 452); 

                auto tg_yyyyyzz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 453); 

                auto tg_yyyyyzz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 454); 

                auto tg_yyyyyzz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 455); 

                auto tg_yyyyyzz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 456); 

                auto tg_yyyyyzz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 457); 

                auto tg_yyyyyzz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 458); 

                auto tg_zzzzzzz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 531); 

                auto tg_zzzzzzz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 532); 

                auto tg_zzzzzzz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 533); 

                auto tg_zzzzzzz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 534); 

                auto tg_zzzzzzz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 535); 

                auto tg_zzzzzzz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 536); 

                auto tg_zzzzzzz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 537); 

                auto tg_zzzzzzz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 538); 

                auto tg_zzzzzzz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 539); 

                auto tg_yyyyyy_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 315); 

                auto tg_yyyyyy_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 316); 

                auto tg_yyyyyy_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 317); 

                auto tg_yyyyyy_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 318); 

                auto tg_yyyyyy_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 319); 

                auto tg_yyyyyy_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 320); 

                auto tg_yyyyyy_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 321); 

                auto tg_yyyyyy_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 322); 

                auto tg_yyyyyy_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 323); 

                auto tg_yyyyyy_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 324); 

                auto tg_yyyyyy_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 325); 

                auto tg_yyyyyy_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 326); 

                auto tg_yyyyyy_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 327); 

                auto tg_yyyyyy_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 328); 

                auto tg_yyyyyy_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 329); 

                auto tg_yyyyyz_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 330); 

                auto tg_yyyyyz_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 331); 

                auto tg_yyyyyz_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 332); 

                auto tg_yyyyyz_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 333); 

                auto tg_yyyyyz_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 334); 

                auto tg_yyyyyz_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 335); 

                auto tg_yyyyyz_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 336); 

                auto tg_yyyyyz_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 337); 

                auto tg_yyyyyz_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 338); 

                auto tg_yyyyyz_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 339); 

                auto tg_yyyyyz_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 340); 

                auto tg_yyyyyz_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 341); 

                auto tg_yyyyyz_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 342); 

                auto tg_yyyyyz_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 343); 

                auto tg_yyyyyz_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 344); 

                auto tg_yyyyzz_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 345); 

                auto tg_yyyyzz_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 346); 

                auto tg_yyyyzz_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 347); 

                auto tg_yyyyzz_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 348); 

                auto tg_yyyyzz_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 349); 

                auto tg_yyyyzz_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 350); 

                auto tg_yyyyzz_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 351); 

                auto tg_yyyyzz_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 352); 

                auto tg_yyyyzz_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 353); 

                auto tg_yyyyyy_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 315); 

                auto tg_yyyyyy_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 316); 

                auto tg_yyyyyy_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 317); 

                auto tg_yyyyyy_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 318); 

                auto tg_yyyyyy_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 319); 

                auto tg_yyyyyy_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 320); 

                auto tg_yyyyyy_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 321); 

                auto tg_yyyyyy_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 322); 

                auto tg_yyyyyy_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 323); 

                auto tg_yyyyyy_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 324); 

                auto tg_yyyyyy_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 325); 

                auto tg_yyyyyy_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 326); 

                auto tg_yyyyyy_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 327); 

                auto tg_yyyyyy_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 328); 

                auto tg_yyyyyy_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 329); 

                auto tg_yyyyyz_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 330); 

                auto tg_yyyyyz_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 331); 

                auto tg_yyyyyz_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 332); 

                auto tg_yyyyyz_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 333); 

                auto tg_yyyyyz_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 334); 

                auto tg_yyyyyz_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 335); 

                auto tg_yyyyyz_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 336); 

                auto tg_yyyyyz_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 337); 

                auto tg_yyyyyz_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 338); 

                auto tg_yyyyyz_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 339); 

                auto tg_yyyyyz_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 340); 

                auto tg_yyyyyz_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 341); 

                auto tg_yyyyyz_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 342); 

                auto tg_yyyyyz_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 343); 

                auto tg_yyyyyz_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 344); 

                auto tg_yyyyzz_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 345); 

                auto tg_yyyyzz_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 346); 

                auto tg_yyyyzz_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 347); 

                auto tg_yyyyzz_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 348); 

                auto tg_yyyyzz_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 349); 

                auto tg_yyyyzz_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 350); 

                auto tg_yyyyzz_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 351); 

                auto tg_yyyyzz_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 352); 

                auto tg_yyyyzz_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 353); 

                auto tg_yyyyyyy_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 280); 

                auto tg_yyyyyyy_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 281); 

                auto tg_yyyyyyy_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 282); 

                auto tg_yyyyyyy_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 283); 

                auto tg_yyyyyyy_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 284); 

                auto tg_yyyyyyy_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 285); 

                auto tg_yyyyyyy_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 286); 

                auto tg_yyyyyyy_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 287); 

                auto tg_yyyyyyy_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 288); 

                auto tg_yyyyyyy_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 289); 

                auto tg_yyyyyyz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 290); 

                auto tg_yyyyyyz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 291); 

                auto tg_yyyyyyz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 292); 

                auto tg_yyyyyyz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 293); 

                auto tg_yyyyyyz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 294); 

                auto tg_yyyyyyz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 295); 

                auto tg_yyyyyyz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 296); 

                auto tg_yyyyyyz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 297); 

                auto tg_yyyyyyz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 298); 

                auto tg_yyyyyyz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 299); 

                auto tg_yyyyyzz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 300); 

                auto tg_yyyyyzz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 301); 

                auto tg_yyyyyzz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 302); 

                auto tg_yyyyyzz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 303); 

                auto tg_yyyyyzz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 304); 

                auto tg_yyyyyzz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 305); 

                auto tg_zzzzzzz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 356); 

                auto tg_zzzzzzz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 357); 

                auto tg_zzzzzzz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 358); 

                auto tg_zzzzzzz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 359); 

                // set up pointers to integrals

                auto tg_xzzzzzzz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 531); 

                auto tg_xzzzzzzz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 532); 

                auto tg_xzzzzzzz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 533); 

                auto tg_xzzzzzzz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 534); 

                auto tg_xzzzzzzz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 535); 

                auto tg_xzzzzzzz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 536); 

                auto tg_xzzzzzzz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 537); 

                auto tg_xzzzzzzz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 538); 

                auto tg_xzzzzzzz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 539); 

                auto tg_yyyyyyyy_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 540); 

                auto tg_yyyyyyyy_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 541); 

                auto tg_yyyyyyyy_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 542); 

                auto tg_yyyyyyyy_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 543); 

                auto tg_yyyyyyyy_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 544); 

                auto tg_yyyyyyyy_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 545); 

                auto tg_yyyyyyyy_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 546); 

                auto tg_yyyyyyyy_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 547); 

                auto tg_yyyyyyyy_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 548); 

                auto tg_yyyyyyyy_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 549); 

                auto tg_yyyyyyyy_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 550); 

                auto tg_yyyyyyyy_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 551); 

                auto tg_yyyyyyyy_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 552); 

                auto tg_yyyyyyyy_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 553); 

                auto tg_yyyyyyyy_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 554); 

                auto tg_yyyyyyyz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 555); 

                auto tg_yyyyyyyz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 556); 

                auto tg_yyyyyyyz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 557); 

                auto tg_yyyyyyyz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 558); 

                auto tg_yyyyyyyz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 559); 

                auto tg_yyyyyyyz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 560); 

                auto tg_yyyyyyyz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 561); 

                auto tg_yyyyyyyz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 562); 

                auto tg_yyyyyyyz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 563); 

                auto tg_yyyyyyyz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 564); 

                auto tg_yyyyyyyz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 565); 

                auto tg_yyyyyyyz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 566); 

                auto tg_yyyyyyyz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 567); 

                auto tg_yyyyyyyz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 568); 

                auto tg_yyyyyyyz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 569); 

                auto tg_yyyyyyzz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 570); 

                auto tg_yyyyyyzz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 571); 

                auto tg_yyyyyyzz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 572); 

                auto tg_yyyyyyzz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 573); 

                auto tg_yyyyyyzz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 574); 

                auto tg_yyyyyyzz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 575); 

                auto tg_yyyyyyzz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 576); 

                auto tg_yyyyyyzz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 577); 

                auto tg_yyyyyyzz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 578); 

                // Batch of Integrals (531,579)

                #pragma omp simd aligned(fxn, fza, tg_xzzzzzzz_xyyy_0, tg_xzzzzzzz_xyyz_0, tg_xzzzzzzz_xyzz_0, \
                                         tg_xzzzzzzz_xzzz_0, tg_xzzzzzzz_yyyy_0, tg_xzzzzzzz_yyyz_0, tg_xzzzzzzz_yyzz_0, \
                                         tg_xzzzzzzz_yzzz_0, tg_xzzzzzzz_zzzz_0, tg_yyyyyy_xxxx_0, tg_yyyyyy_xxxx_1, \
                                         tg_yyyyyy_xxxy_0, tg_yyyyyy_xxxy_1, tg_yyyyyy_xxxz_0, tg_yyyyyy_xxxz_1, \
                                         tg_yyyyyy_xxyy_0, tg_yyyyyy_xxyy_1, tg_yyyyyy_xxyz_0, tg_yyyyyy_xxyz_1, \
                                         tg_yyyyyy_xxzz_0, tg_yyyyyy_xxzz_1, tg_yyyyyy_xyyy_0, tg_yyyyyy_xyyy_1, \
                                         tg_yyyyyy_xyyz_0, tg_yyyyyy_xyyz_1, tg_yyyyyy_xyzz_0, tg_yyyyyy_xyzz_1, \
                                         tg_yyyyyy_xzzz_0, tg_yyyyyy_xzzz_1, tg_yyyyyy_yyyy_0, tg_yyyyyy_yyyy_1, \
                                         tg_yyyyyy_yyyz_0, tg_yyyyyy_yyyz_1, tg_yyyyyy_yyzz_0, tg_yyyyyy_yyzz_1, \
                                         tg_yyyyyy_yzzz_0, tg_yyyyyy_yzzz_1, tg_yyyyyy_zzzz_0, tg_yyyyyy_zzzz_1, \
                                         tg_yyyyyyy_xxx_1, tg_yyyyyyy_xxxx_0, tg_yyyyyyy_xxxx_1, tg_yyyyyyy_xxxy_0, \
                                         tg_yyyyyyy_xxxy_1, tg_yyyyyyy_xxxz_0, tg_yyyyyyy_xxxz_1, tg_yyyyyyy_xxy_1, \
                                         tg_yyyyyyy_xxyy_0, tg_yyyyyyy_xxyy_1, tg_yyyyyyy_xxyz_0, tg_yyyyyyy_xxyz_1, \
                                         tg_yyyyyyy_xxz_1, tg_yyyyyyy_xxzz_0, tg_yyyyyyy_xxzz_1, tg_yyyyyyy_xyy_1, \
                                         tg_yyyyyyy_xyyy_0, tg_yyyyyyy_xyyy_1, tg_yyyyyyy_xyyz_0, tg_yyyyyyy_xyyz_1, \
                                         tg_yyyyyyy_xyz_1, tg_yyyyyyy_xyzz_0, tg_yyyyyyy_xyzz_1, tg_yyyyyyy_xzz_1, \
                                         tg_yyyyyyy_xzzz_0, tg_yyyyyyy_xzzz_1, tg_yyyyyyy_yyy_1, tg_yyyyyyy_yyyy_0, \
                                         tg_yyyyyyy_yyyy_1, tg_yyyyyyy_yyyz_0, tg_yyyyyyy_yyyz_1, tg_yyyyyyy_yyz_1, \
                                         tg_yyyyyyy_yyzz_0, tg_yyyyyyy_yyzz_1, tg_yyyyyyy_yzz_1, tg_yyyyyyy_yzzz_0, \
                                         tg_yyyyyyy_yzzz_1, tg_yyyyyyy_zzz_1, tg_yyyyyyy_zzzz_0, tg_yyyyyyy_zzzz_1, \
                                         tg_yyyyyyyy_xxxx_0, tg_yyyyyyyy_xxxy_0, tg_yyyyyyyy_xxxz_0, tg_yyyyyyyy_xxyy_0, \
                                         tg_yyyyyyyy_xxyz_0, tg_yyyyyyyy_xxzz_0, tg_yyyyyyyy_xyyy_0, tg_yyyyyyyy_xyyz_0, \
                                         tg_yyyyyyyy_xyzz_0, tg_yyyyyyyy_xzzz_0, tg_yyyyyyyy_yyyy_0, tg_yyyyyyyy_yyyz_0, \
                                         tg_yyyyyyyy_yyzz_0, tg_yyyyyyyy_yzzz_0, tg_yyyyyyyy_zzzz_0, tg_yyyyyyyz_xxxx_0, \
                                         tg_yyyyyyyz_xxxy_0, tg_yyyyyyyz_xxxz_0, tg_yyyyyyyz_xxyy_0, tg_yyyyyyyz_xxyz_0, \
                                         tg_yyyyyyyz_xxzz_0, tg_yyyyyyyz_xyyy_0, tg_yyyyyyyz_xyyz_0, tg_yyyyyyyz_xyzz_0, \
                                         tg_yyyyyyyz_xzzz_0, tg_yyyyyyyz_yyyy_0, tg_yyyyyyyz_yyyz_0, tg_yyyyyyyz_yyzz_0, \
                                         tg_yyyyyyyz_yzzz_0, tg_yyyyyyyz_zzzz_0, tg_yyyyyyz_xxx_1, tg_yyyyyyz_xxxx_0, \
                                         tg_yyyyyyz_xxxx_1, tg_yyyyyyz_xxxy_0, tg_yyyyyyz_xxxy_1, tg_yyyyyyz_xxxz_0, \
                                         tg_yyyyyyz_xxxz_1, tg_yyyyyyz_xxy_1, tg_yyyyyyz_xxyy_0, tg_yyyyyyz_xxyy_1, \
                                         tg_yyyyyyz_xxyz_0, tg_yyyyyyz_xxyz_1, tg_yyyyyyz_xxz_1, tg_yyyyyyz_xxzz_0, \
                                         tg_yyyyyyz_xxzz_1, tg_yyyyyyz_xyy_1, tg_yyyyyyz_xyyy_0, tg_yyyyyyz_xyyy_1, \
                                         tg_yyyyyyz_xyyz_0, tg_yyyyyyz_xyyz_1, tg_yyyyyyz_xyz_1, tg_yyyyyyz_xyzz_0, \
                                         tg_yyyyyyz_xyzz_1, tg_yyyyyyz_xzz_1, tg_yyyyyyz_xzzz_0, tg_yyyyyyz_xzzz_1, \
                                         tg_yyyyyyz_yyy_1, tg_yyyyyyz_yyyy_0, tg_yyyyyyz_yyyy_1, tg_yyyyyyz_yyyz_0, \
                                         tg_yyyyyyz_yyyz_1, tg_yyyyyyz_yyz_1, tg_yyyyyyz_yyzz_0, tg_yyyyyyz_yyzz_1, \
                                         tg_yyyyyyz_yzz_1, tg_yyyyyyz_yzzz_0, tg_yyyyyyz_yzzz_1, tg_yyyyyyz_zzz_1, \
                                         tg_yyyyyyz_zzzz_0, tg_yyyyyyz_zzzz_1, tg_yyyyyyzz_xxxx_0, tg_yyyyyyzz_xxxy_0, \
                                         tg_yyyyyyzz_xxxz_0, tg_yyyyyyzz_xxyy_0, tg_yyyyyyzz_xxyz_0, tg_yyyyyyzz_xxzz_0, \
                                         tg_yyyyyyzz_xyyy_0, tg_yyyyyyzz_xyyz_0, tg_yyyyyyzz_xyzz_0, tg_yyyyyz_xxxx_0, \
                                         tg_yyyyyz_xxxx_1, tg_yyyyyz_xxxy_0, tg_yyyyyz_xxxy_1, tg_yyyyyz_xxxz_0, \
                                         tg_yyyyyz_xxxz_1, tg_yyyyyz_xxyy_0, tg_yyyyyz_xxyy_1, tg_yyyyyz_xxyz_0, \
                                         tg_yyyyyz_xxyz_1, tg_yyyyyz_xxzz_0, tg_yyyyyz_xxzz_1, tg_yyyyyz_xyyy_0, \
                                         tg_yyyyyz_xyyy_1, tg_yyyyyz_xyyz_0, tg_yyyyyz_xyyz_1, tg_yyyyyz_xyzz_0, \
                                         tg_yyyyyz_xyzz_1, tg_yyyyyz_xzzz_0, tg_yyyyyz_xzzz_1, tg_yyyyyz_yyyy_0, \
                                         tg_yyyyyz_yyyy_1, tg_yyyyyz_yyyz_0, tg_yyyyyz_yyyz_1, tg_yyyyyz_yyzz_0, \
                                         tg_yyyyyz_yyzz_1, tg_yyyyyz_yzzz_0, tg_yyyyyz_yzzz_1, tg_yyyyyz_zzzz_0, \
                                         tg_yyyyyz_zzzz_1, tg_yyyyyzz_xxx_1, tg_yyyyyzz_xxxx_0, tg_yyyyyzz_xxxx_1, \
                                         tg_yyyyyzz_xxxy_0, tg_yyyyyzz_xxxy_1, tg_yyyyyzz_xxxz_0, tg_yyyyyzz_xxxz_1, \
                                         tg_yyyyyzz_xxy_1, tg_yyyyyzz_xxyy_0, tg_yyyyyzz_xxyy_1, tg_yyyyyzz_xxyz_0, \
                                         tg_yyyyyzz_xxyz_1, tg_yyyyyzz_xxz_1, tg_yyyyyzz_xxzz_0, tg_yyyyyzz_xxzz_1, \
                                         tg_yyyyyzz_xyy_1, tg_yyyyyzz_xyyy_0, tg_yyyyyzz_xyyy_1, tg_yyyyyzz_xyyz_0, \
                                         tg_yyyyyzz_xyyz_1, tg_yyyyyzz_xyz_1, tg_yyyyyzz_xyzz_0, tg_yyyyyzz_xyzz_1, \
                                         tg_yyyyyzz_xzz_1, tg_yyyyzz_xxxx_0, tg_yyyyzz_xxxx_1, tg_yyyyzz_xxxy_0, \
                                         tg_yyyyzz_xxxy_1, tg_yyyyzz_xxxz_0, tg_yyyyzz_xxxz_1, tg_yyyyzz_xxyy_0, \
                                         tg_yyyyzz_xxyy_1, tg_yyyyzz_xxyz_0, tg_yyyyzz_xxyz_1, tg_yyyyzz_xxzz_0, \
                                         tg_yyyyzz_xxzz_1, tg_yyyyzz_xyyy_0, tg_yyyyzz_xyyy_1, tg_yyyyzz_xyyz_0, \
                                         tg_yyyyzz_xyyz_1, tg_yyyyzz_xyzz_0, tg_yyyyzz_xyzz_1, tg_zzzzzzz_xyyy_0, \
                                         tg_zzzzzzz_xyyy_1, tg_zzzzzzz_xyyz_0, tg_zzzzzzz_xyyz_1, tg_zzzzzzz_xyzz_0, \
                                         tg_zzzzzzz_xyzz_1, tg_zzzzzzz_xzzz_0, tg_zzzzzzz_xzzz_1, tg_zzzzzzz_yyy_1, \
                                         tg_zzzzzzz_yyyy_0, tg_zzzzzzz_yyyy_1, tg_zzzzzzz_yyyz_0, tg_zzzzzzz_yyyz_1, \
                                         tg_zzzzzzz_yyz_1, tg_zzzzzzz_yyzz_0, tg_zzzzzzz_yyzz_1, tg_zzzzzzz_yzz_1, \
                                         tg_zzzzzzz_yzzz_0, tg_zzzzzzz_yzzz_1, tg_zzzzzzz_zzz_1, tg_zzzzzzz_zzzz_0, \
                                         tg_zzzzzzz_zzzz_1, wp_x, wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xzzzzzzz_xyyy_0[j] = pb_x * tg_zzzzzzz_xyyy_0[j] + wp_x[j] * tg_zzzzzzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_yyy_1[j];

                    tg_xzzzzzzz_xyyz_0[j] = pb_x * tg_zzzzzzz_xyyz_0[j] + wp_x[j] * tg_zzzzzzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_yyz_1[j];

                    tg_xzzzzzzz_xyzz_0[j] = pb_x * tg_zzzzzzz_xyzz_0[j] + wp_x[j] * tg_zzzzzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_yzz_1[j];

                    tg_xzzzzzzz_xzzz_0[j] = pb_x * tg_zzzzzzz_xzzz_0[j] + wp_x[j] * tg_zzzzzzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_zzz_1[j];

                    tg_xzzzzzzz_yyyy_0[j] = pb_x * tg_zzzzzzz_yyyy_0[j] + wp_x[j] * tg_zzzzzzz_yyyy_1[j];

                    tg_xzzzzzzz_yyyz_0[j] = pb_x * tg_zzzzzzz_yyyz_0[j] + wp_x[j] * tg_zzzzzzz_yyyz_1[j];

                    tg_xzzzzzzz_yyzz_0[j] = pb_x * tg_zzzzzzz_yyzz_0[j] + wp_x[j] * tg_zzzzzzz_yyzz_1[j];

                    tg_xzzzzzzz_yzzz_0[j] = pb_x * tg_zzzzzzz_yzzz_0[j] + wp_x[j] * tg_zzzzzzz_yzzz_1[j];

                    tg_xzzzzzzz_zzzz_0[j] = pb_x * tg_zzzzzzz_zzzz_0[j] + wp_x[j] * tg_zzzzzzz_zzzz_1[j];

                    tg_yyyyyyyy_xxxx_0[j] = pb_y * tg_yyyyyyy_xxxx_0[j] + wp_y[j] * tg_yyyyyyy_xxxx_1[j] + 3.5 * fl1_fx * tg_yyyyyy_xxxx_0[j] - 3.5 * fl1_fx * fl1_fza * tg_yyyyyy_xxxx_1[j];

                    tg_yyyyyyyy_xxxy_0[j] = pb_y * tg_yyyyyyy_xxxy_0[j] + wp_y[j] * tg_yyyyyyy_xxxy_1[j] + 3.5 * fl1_fx * tg_yyyyyy_xxxy_0[j] - 3.5 * fl1_fx * fl1_fza * tg_yyyyyy_xxxy_1[j] + 0.5 * fl1_fxn * tg_yyyyyyy_xxx_1[j];

                    tg_yyyyyyyy_xxxz_0[j] = pb_y * tg_yyyyyyy_xxxz_0[j] + wp_y[j] * tg_yyyyyyy_xxxz_1[j] + 3.5 * fl1_fx * tg_yyyyyy_xxxz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_yyyyyy_xxxz_1[j];

                    tg_yyyyyyyy_xxyy_0[j] = pb_y * tg_yyyyyyy_xxyy_0[j] + wp_y[j] * tg_yyyyyyy_xxyy_1[j] + 3.5 * fl1_fx * tg_yyyyyy_xxyy_0[j] - 3.5 * fl1_fx * fl1_fza * tg_yyyyyy_xxyy_1[j] + fl1_fxn * tg_yyyyyyy_xxy_1[j];

                    tg_yyyyyyyy_xxyz_0[j] = pb_y * tg_yyyyyyy_xxyz_0[j] + wp_y[j] * tg_yyyyyyy_xxyz_1[j] + 3.5 * fl1_fx * tg_yyyyyy_xxyz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_yyyyyy_xxyz_1[j] + 0.5 * fl1_fxn * tg_yyyyyyy_xxz_1[j];

                    tg_yyyyyyyy_xxzz_0[j] = pb_y * tg_yyyyyyy_xxzz_0[j] + wp_y[j] * tg_yyyyyyy_xxzz_1[j] + 3.5 * fl1_fx * tg_yyyyyy_xxzz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_yyyyyy_xxzz_1[j];

                    tg_yyyyyyyy_xyyy_0[j] = pb_y * tg_yyyyyyy_xyyy_0[j] + wp_y[j] * tg_yyyyyyy_xyyy_1[j] + 3.5 * fl1_fx * tg_yyyyyy_xyyy_0[j] - 3.5 * fl1_fx * fl1_fza * tg_yyyyyy_xyyy_1[j] + 1.5 * fl1_fxn * tg_yyyyyyy_xyy_1[j];

                    tg_yyyyyyyy_xyyz_0[j] = pb_y * tg_yyyyyyy_xyyz_0[j] + wp_y[j] * tg_yyyyyyy_xyyz_1[j] + 3.5 * fl1_fx * tg_yyyyyy_xyyz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_yyyyyy_xyyz_1[j] + fl1_fxn * tg_yyyyyyy_xyz_1[j];

                    tg_yyyyyyyy_xyzz_0[j] = pb_y * tg_yyyyyyy_xyzz_0[j] + wp_y[j] * tg_yyyyyyy_xyzz_1[j] + 3.5 * fl1_fx * tg_yyyyyy_xyzz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_yyyyyy_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyyy_xzz_1[j];

                    tg_yyyyyyyy_xzzz_0[j] = pb_y * tg_yyyyyyy_xzzz_0[j] + wp_y[j] * tg_yyyyyyy_xzzz_1[j] + 3.5 * fl1_fx * tg_yyyyyy_xzzz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_yyyyyy_xzzz_1[j];

                    tg_yyyyyyyy_yyyy_0[j] = pb_y * tg_yyyyyyy_yyyy_0[j] + wp_y[j] * tg_yyyyyyy_yyyy_1[j] + 3.5 * fl1_fx * tg_yyyyyy_yyyy_0[j] - 3.5 * fl1_fx * fl1_fza * tg_yyyyyy_yyyy_1[j] + 2.0 * fl1_fxn * tg_yyyyyyy_yyy_1[j];

                    tg_yyyyyyyy_yyyz_0[j] = pb_y * tg_yyyyyyy_yyyz_0[j] + wp_y[j] * tg_yyyyyyy_yyyz_1[j] + 3.5 * fl1_fx * tg_yyyyyy_yyyz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_yyyyyy_yyyz_1[j] + 1.5 * fl1_fxn * tg_yyyyyyy_yyz_1[j];

                    tg_yyyyyyyy_yyzz_0[j] = pb_y * tg_yyyyyyy_yyzz_0[j] + wp_y[j] * tg_yyyyyyy_yyzz_1[j] + 3.5 * fl1_fx * tg_yyyyyy_yyzz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_yyyyyy_yyzz_1[j] + fl1_fxn * tg_yyyyyyy_yzz_1[j];

                    tg_yyyyyyyy_yzzz_0[j] = pb_y * tg_yyyyyyy_yzzz_0[j] + wp_y[j] * tg_yyyyyyy_yzzz_1[j] + 3.5 * fl1_fx * tg_yyyyyy_yzzz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_yyyyyy_yzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyyy_zzz_1[j];

                    tg_yyyyyyyy_zzzz_0[j] = pb_y * tg_yyyyyyy_zzzz_0[j] + wp_y[j] * tg_yyyyyyy_zzzz_1[j] + 3.5 * fl1_fx * tg_yyyyyy_zzzz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_yyyyyy_zzzz_1[j];

                    tg_yyyyyyyz_xxxx_0[j] = pb_y * tg_yyyyyyz_xxxx_0[j] + wp_y[j] * tg_yyyyyyz_xxxx_1[j] + 3.0 * fl1_fx * tg_yyyyyz_xxxx_0[j] - 3.0 * fl1_fx * fl1_fza * tg_yyyyyz_xxxx_1[j];

                    tg_yyyyyyyz_xxxy_0[j] = pb_y * tg_yyyyyyz_xxxy_0[j] + wp_y[j] * tg_yyyyyyz_xxxy_1[j] + 3.0 * fl1_fx * tg_yyyyyz_xxxy_0[j] - 3.0 * fl1_fx * fl1_fza * tg_yyyyyz_xxxy_1[j] + 0.5 * fl1_fxn * tg_yyyyyyz_xxx_1[j];

                    tg_yyyyyyyz_xxxz_0[j] = pb_y * tg_yyyyyyz_xxxz_0[j] + wp_y[j] * tg_yyyyyyz_xxxz_1[j] + 3.0 * fl1_fx * tg_yyyyyz_xxxz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_yyyyyz_xxxz_1[j];

                    tg_yyyyyyyz_xxyy_0[j] = pb_y * tg_yyyyyyz_xxyy_0[j] + wp_y[j] * tg_yyyyyyz_xxyy_1[j] + 3.0 * fl1_fx * tg_yyyyyz_xxyy_0[j] - 3.0 * fl1_fx * fl1_fza * tg_yyyyyz_xxyy_1[j] + fl1_fxn * tg_yyyyyyz_xxy_1[j];

                    tg_yyyyyyyz_xxyz_0[j] = pb_y * tg_yyyyyyz_xxyz_0[j] + wp_y[j] * tg_yyyyyyz_xxyz_1[j] + 3.0 * fl1_fx * tg_yyyyyz_xxyz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_yyyyyz_xxyz_1[j] + 0.5 * fl1_fxn * tg_yyyyyyz_xxz_1[j];

                    tg_yyyyyyyz_xxzz_0[j] = pb_y * tg_yyyyyyz_xxzz_0[j] + wp_y[j] * tg_yyyyyyz_xxzz_1[j] + 3.0 * fl1_fx * tg_yyyyyz_xxzz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_yyyyyz_xxzz_1[j];

                    tg_yyyyyyyz_xyyy_0[j] = pb_y * tg_yyyyyyz_xyyy_0[j] + wp_y[j] * tg_yyyyyyz_xyyy_1[j] + 3.0 * fl1_fx * tg_yyyyyz_xyyy_0[j] - 3.0 * fl1_fx * fl1_fza * tg_yyyyyz_xyyy_1[j] + 1.5 * fl1_fxn * tg_yyyyyyz_xyy_1[j];

                    tg_yyyyyyyz_xyyz_0[j] = pb_y * tg_yyyyyyz_xyyz_0[j] + wp_y[j] * tg_yyyyyyz_xyyz_1[j] + 3.0 * fl1_fx * tg_yyyyyz_xyyz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_yyyyyz_xyyz_1[j] + fl1_fxn * tg_yyyyyyz_xyz_1[j];

                    tg_yyyyyyyz_xyzz_0[j] = pb_y * tg_yyyyyyz_xyzz_0[j] + wp_y[j] * tg_yyyyyyz_xyzz_1[j] + 3.0 * fl1_fx * tg_yyyyyz_xyzz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_yyyyyz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyyz_xzz_1[j];

                    tg_yyyyyyyz_xzzz_0[j] = pb_y * tg_yyyyyyz_xzzz_0[j] + wp_y[j] * tg_yyyyyyz_xzzz_1[j] + 3.0 * fl1_fx * tg_yyyyyz_xzzz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_yyyyyz_xzzz_1[j];

                    tg_yyyyyyyz_yyyy_0[j] = pb_y * tg_yyyyyyz_yyyy_0[j] + wp_y[j] * tg_yyyyyyz_yyyy_1[j] + 3.0 * fl1_fx * tg_yyyyyz_yyyy_0[j] - 3.0 * fl1_fx * fl1_fza * tg_yyyyyz_yyyy_1[j] + 2.0 * fl1_fxn * tg_yyyyyyz_yyy_1[j];

                    tg_yyyyyyyz_yyyz_0[j] = pb_y * tg_yyyyyyz_yyyz_0[j] + wp_y[j] * tg_yyyyyyz_yyyz_1[j] + 3.0 * fl1_fx * tg_yyyyyz_yyyz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_yyyyyz_yyyz_1[j] + 1.5 * fl1_fxn * tg_yyyyyyz_yyz_1[j];

                    tg_yyyyyyyz_yyzz_0[j] = pb_y * tg_yyyyyyz_yyzz_0[j] + wp_y[j] * tg_yyyyyyz_yyzz_1[j] + 3.0 * fl1_fx * tg_yyyyyz_yyzz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_yyyyyz_yyzz_1[j] + fl1_fxn * tg_yyyyyyz_yzz_1[j];

                    tg_yyyyyyyz_yzzz_0[j] = pb_y * tg_yyyyyyz_yzzz_0[j] + wp_y[j] * tg_yyyyyyz_yzzz_1[j] + 3.0 * fl1_fx * tg_yyyyyz_yzzz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_yyyyyz_yzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyyz_zzz_1[j];

                    tg_yyyyyyyz_zzzz_0[j] = pb_y * tg_yyyyyyz_zzzz_0[j] + wp_y[j] * tg_yyyyyyz_zzzz_1[j] + 3.0 * fl1_fx * tg_yyyyyz_zzzz_0[j] - 3.0 * fl1_fx * fl1_fza * tg_yyyyyz_zzzz_1[j];

                    tg_yyyyyyzz_xxxx_0[j] = pb_y * tg_yyyyyzz_xxxx_0[j] + wp_y[j] * tg_yyyyyzz_xxxx_1[j] + 2.5 * fl1_fx * tg_yyyyzz_xxxx_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyyzz_xxxx_1[j];

                    tg_yyyyyyzz_xxxy_0[j] = pb_y * tg_yyyyyzz_xxxy_0[j] + wp_y[j] * tg_yyyyyzz_xxxy_1[j] + 2.5 * fl1_fx * tg_yyyyzz_xxxy_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyyzz_xxxy_1[j] + 0.5 * fl1_fxn * tg_yyyyyzz_xxx_1[j];

                    tg_yyyyyyzz_xxxz_0[j] = pb_y * tg_yyyyyzz_xxxz_0[j] + wp_y[j] * tg_yyyyyzz_xxxz_1[j] + 2.5 * fl1_fx * tg_yyyyzz_xxxz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyyzz_xxxz_1[j];

                    tg_yyyyyyzz_xxyy_0[j] = pb_y * tg_yyyyyzz_xxyy_0[j] + wp_y[j] * tg_yyyyyzz_xxyy_1[j] + 2.5 * fl1_fx * tg_yyyyzz_xxyy_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyyzz_xxyy_1[j] + fl1_fxn * tg_yyyyyzz_xxy_1[j];

                    tg_yyyyyyzz_xxyz_0[j] = pb_y * tg_yyyyyzz_xxyz_0[j] + wp_y[j] * tg_yyyyyzz_xxyz_1[j] + 2.5 * fl1_fx * tg_yyyyzz_xxyz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyyzz_xxyz_1[j] + 0.5 * fl1_fxn * tg_yyyyyzz_xxz_1[j];

                    tg_yyyyyyzz_xxzz_0[j] = pb_y * tg_yyyyyzz_xxzz_0[j] + wp_y[j] * tg_yyyyyzz_xxzz_1[j] + 2.5 * fl1_fx * tg_yyyyzz_xxzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyyzz_xxzz_1[j];

                    tg_yyyyyyzz_xyyy_0[j] = pb_y * tg_yyyyyzz_xyyy_0[j] + wp_y[j] * tg_yyyyyzz_xyyy_1[j] + 2.5 * fl1_fx * tg_yyyyzz_xyyy_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyyzz_xyyy_1[j] + 1.5 * fl1_fxn * tg_yyyyyzz_xyy_1[j];

                    tg_yyyyyyzz_xyyz_0[j] = pb_y * tg_yyyyyzz_xyyz_0[j] + wp_y[j] * tg_yyyyyzz_xyyz_1[j] + 2.5 * fl1_fx * tg_yyyyzz_xyyz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyyzz_xyyz_1[j] + fl1_fxn * tg_yyyyyzz_xyz_1[j];

                    tg_yyyyyyzz_xyzz_0[j] = pb_y * tg_yyyyyzz_xyzz_0[j] + wp_y[j] * tg_yyyyyzz_xyzz_1[j] + 2.5 * fl1_fx * tg_yyyyzz_xyzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyyzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyzz_xzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSLSG_579_627(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (579,627)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {8, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {8, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_7_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_6_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_6_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_7_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto pb_y = r_pb_y[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_y = wpDistances.data(3 * idx + 1);

                // set up pointers to auxilary integrals

                auto tg_yyyyyzz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 459); 

                auto tg_yyyyyzz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 460); 

                auto tg_yyyyyzz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 461); 

                auto tg_yyyyyzz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 462); 

                auto tg_yyyyyzz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 463); 

                auto tg_yyyyyzz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 464); 

                auto tg_yyyyzzz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 465); 

                auto tg_yyyyzzz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 466); 

                auto tg_yyyyzzz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 467); 

                auto tg_yyyyzzz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 468); 

                auto tg_yyyyzzz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 469); 

                auto tg_yyyyzzz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 470); 

                auto tg_yyyyzzz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 471); 

                auto tg_yyyyzzz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 472); 

                auto tg_yyyyzzz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 473); 

                auto tg_yyyyzzz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 474); 

                auto tg_yyyyzzz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 475); 

                auto tg_yyyyzzz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 476); 

                auto tg_yyyyzzz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 477); 

                auto tg_yyyyzzz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 478); 

                auto tg_yyyyzzz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 479); 

                auto tg_yyyzzzz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 480); 

                auto tg_yyyzzzz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 481); 

                auto tg_yyyzzzz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 482); 

                auto tg_yyyzzzz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 483); 

                auto tg_yyyzzzz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 484); 

                auto tg_yyyzzzz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 485); 

                auto tg_yyyzzzz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 486); 

                auto tg_yyyzzzz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 487); 

                auto tg_yyyzzzz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 488); 

                auto tg_yyyzzzz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 489); 

                auto tg_yyyzzzz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 490); 

                auto tg_yyyzzzz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 491); 

                auto tg_yyyzzzz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 492); 

                auto tg_yyyzzzz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 493); 

                auto tg_yyyzzzz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 494); 

                auto tg_yyzzzzz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 495); 

                auto tg_yyzzzzz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 496); 

                auto tg_yyzzzzz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 497); 

                auto tg_yyzzzzz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 498); 

                auto tg_yyzzzzz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 499); 

                auto tg_yyzzzzz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 500); 

                auto tg_yyzzzzz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 501); 

                auto tg_yyzzzzz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 502); 

                auto tg_yyzzzzz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 503); 

                auto tg_yyzzzzz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 504); 

                auto tg_yyzzzzz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 505); 

                auto tg_yyzzzzz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 506); 

                auto tg_yyyyyzz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 459); 

                auto tg_yyyyyzz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 460); 

                auto tg_yyyyyzz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 461); 

                auto tg_yyyyyzz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 462); 

                auto tg_yyyyyzz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 463); 

                auto tg_yyyyyzz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 464); 

                auto tg_yyyyzzz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 465); 

                auto tg_yyyyzzz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 466); 

                auto tg_yyyyzzz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 467); 

                auto tg_yyyyzzz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 468); 

                auto tg_yyyyzzz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 469); 

                auto tg_yyyyzzz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 470); 

                auto tg_yyyyzzz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 471); 

                auto tg_yyyyzzz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 472); 

                auto tg_yyyyzzz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 473); 

                auto tg_yyyyzzz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 474); 

                auto tg_yyyyzzz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 475); 

                auto tg_yyyyzzz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 476); 

                auto tg_yyyyzzz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 477); 

                auto tg_yyyyzzz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 478); 

                auto tg_yyyyzzz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 479); 

                auto tg_yyyzzzz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 480); 

                auto tg_yyyzzzz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 481); 

                auto tg_yyyzzzz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 482); 

                auto tg_yyyzzzz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 483); 

                auto tg_yyyzzzz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 484); 

                auto tg_yyyzzzz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 485); 

                auto tg_yyyzzzz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 486); 

                auto tg_yyyzzzz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 487); 

                auto tg_yyyzzzz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 488); 

                auto tg_yyyzzzz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 489); 

                auto tg_yyyzzzz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 490); 

                auto tg_yyyzzzz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 491); 

                auto tg_yyyzzzz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 492); 

                auto tg_yyyzzzz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 493); 

                auto tg_yyyzzzz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 494); 

                auto tg_yyzzzzz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 495); 

                auto tg_yyzzzzz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 496); 

                auto tg_yyzzzzz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 497); 

                auto tg_yyzzzzz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 498); 

                auto tg_yyzzzzz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 499); 

                auto tg_yyzzzzz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 500); 

                auto tg_yyzzzzz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 501); 

                auto tg_yyzzzzz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 502); 

                auto tg_yyzzzzz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 503); 

                auto tg_yyzzzzz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 504); 

                auto tg_yyzzzzz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 505); 

                auto tg_yyzzzzz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 506); 

                auto tg_yyyyzz_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 354); 

                auto tg_yyyyzz_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 355); 

                auto tg_yyyyzz_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 356); 

                auto tg_yyyyzz_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 357); 

                auto tg_yyyyzz_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 358); 

                auto tg_yyyyzz_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 359); 

                auto tg_yyyzzz_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 360); 

                auto tg_yyyzzz_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 361); 

                auto tg_yyyzzz_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 362); 

                auto tg_yyyzzz_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 363); 

                auto tg_yyyzzz_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 364); 

                auto tg_yyyzzz_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 365); 

                auto tg_yyyzzz_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 366); 

                auto tg_yyyzzz_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 367); 

                auto tg_yyyzzz_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 368); 

                auto tg_yyyzzz_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 369); 

                auto tg_yyyzzz_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 370); 

                auto tg_yyyzzz_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 371); 

                auto tg_yyyzzz_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 372); 

                auto tg_yyyzzz_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 373); 

                auto tg_yyyzzz_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 374); 

                auto tg_yyzzzz_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 375); 

                auto tg_yyzzzz_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 376); 

                auto tg_yyzzzz_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 377); 

                auto tg_yyzzzz_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 378); 

                auto tg_yyzzzz_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 379); 

                auto tg_yyzzzz_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 380); 

                auto tg_yyzzzz_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 381); 

                auto tg_yyzzzz_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 382); 

                auto tg_yyzzzz_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 383); 

                auto tg_yyzzzz_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 384); 

                auto tg_yyzzzz_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 385); 

                auto tg_yyzzzz_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 386); 

                auto tg_yyzzzz_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 387); 

                auto tg_yyzzzz_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 388); 

                auto tg_yyzzzz_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 389); 

                auto tg_yzzzzz_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 390); 

                auto tg_yzzzzz_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 391); 

                auto tg_yzzzzz_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 392); 

                auto tg_yzzzzz_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 393); 

                auto tg_yzzzzz_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 394); 

                auto tg_yzzzzz_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 395); 

                auto tg_yzzzzz_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 396); 

                auto tg_yzzzzz_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 397); 

                auto tg_yzzzzz_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 398); 

                auto tg_yzzzzz_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 399); 

                auto tg_yzzzzz_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 400); 

                auto tg_yzzzzz_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 401); 

                auto tg_yyyyzz_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 354); 

                auto tg_yyyyzz_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 355); 

                auto tg_yyyyzz_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 356); 

                auto tg_yyyyzz_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 357); 

                auto tg_yyyyzz_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 358); 

                auto tg_yyyyzz_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 359); 

                auto tg_yyyzzz_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 360); 

                auto tg_yyyzzz_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 361); 

                auto tg_yyyzzz_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 362); 

                auto tg_yyyzzz_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 363); 

                auto tg_yyyzzz_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 364); 

                auto tg_yyyzzz_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 365); 

                auto tg_yyyzzz_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 366); 

                auto tg_yyyzzz_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 367); 

                auto tg_yyyzzz_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 368); 

                auto tg_yyyzzz_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 369); 

                auto tg_yyyzzz_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 370); 

                auto tg_yyyzzz_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 371); 

                auto tg_yyyzzz_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 372); 

                auto tg_yyyzzz_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 373); 

                auto tg_yyyzzz_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 374); 

                auto tg_yyzzzz_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 375); 

                auto tg_yyzzzz_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 376); 

                auto tg_yyzzzz_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 377); 

                auto tg_yyzzzz_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 378); 

                auto tg_yyzzzz_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 379); 

                auto tg_yyzzzz_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 380); 

                auto tg_yyzzzz_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 381); 

                auto tg_yyzzzz_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 382); 

                auto tg_yyzzzz_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 383); 

                auto tg_yyzzzz_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 384); 

                auto tg_yyzzzz_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 385); 

                auto tg_yyzzzz_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 386); 

                auto tg_yyzzzz_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 387); 

                auto tg_yyzzzz_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 388); 

                auto tg_yyzzzz_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 389); 

                auto tg_yzzzzz_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 390); 

                auto tg_yzzzzz_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 391); 

                auto tg_yzzzzz_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 392); 

                auto tg_yzzzzz_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 393); 

                auto tg_yzzzzz_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 394); 

                auto tg_yzzzzz_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 395); 

                auto tg_yzzzzz_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 396); 

                auto tg_yzzzzz_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 397); 

                auto tg_yzzzzz_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 398); 

                auto tg_yzzzzz_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 399); 

                auto tg_yzzzzz_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 400); 

                auto tg_yzzzzz_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 401); 

                auto tg_yyyyyzz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 306); 

                auto tg_yyyyyzz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 307); 

                auto tg_yyyyyzz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 308); 

                auto tg_yyyyyzz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 309); 

                auto tg_yyyyzzz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 310); 

                auto tg_yyyyzzz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 311); 

                auto tg_yyyyzzz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 312); 

                auto tg_yyyyzzz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 313); 

                auto tg_yyyyzzz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 314); 

                auto tg_yyyyzzz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 315); 

                auto tg_yyyyzzz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 316); 

                auto tg_yyyyzzz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 317); 

                auto tg_yyyyzzz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 318); 

                auto tg_yyyyzzz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 319); 

                auto tg_yyyzzzz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 320); 

                auto tg_yyyzzzz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 321); 

                auto tg_yyyzzzz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 322); 

                auto tg_yyyzzzz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 323); 

                auto tg_yyyzzzz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 324); 

                auto tg_yyyzzzz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 325); 

                auto tg_yyyzzzz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 326); 

                auto tg_yyyzzzz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 327); 

                auto tg_yyyzzzz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 328); 

                auto tg_yyyzzzz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 329); 

                auto tg_yyzzzzz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 330); 

                auto tg_yyzzzzz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 331); 

                auto tg_yyzzzzz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 332); 

                auto tg_yyzzzzz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 333); 

                auto tg_yyzzzzz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 334); 

                auto tg_yyzzzzz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 335); 

                auto tg_yyzzzzz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 336); 

                auto tg_yyzzzzz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 337); 

                // set up pointers to integrals

                auto tg_yyyyyyzz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 579); 

                auto tg_yyyyyyzz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 580); 

                auto tg_yyyyyyzz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 581); 

                auto tg_yyyyyyzz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 582); 

                auto tg_yyyyyyzz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 583); 

                auto tg_yyyyyyzz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 584); 

                auto tg_yyyyyzzz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 585); 

                auto tg_yyyyyzzz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 586); 

                auto tg_yyyyyzzz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 587); 

                auto tg_yyyyyzzz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 588); 

                auto tg_yyyyyzzz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 589); 

                auto tg_yyyyyzzz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 590); 

                auto tg_yyyyyzzz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 591); 

                auto tg_yyyyyzzz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 592); 

                auto tg_yyyyyzzz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 593); 

                auto tg_yyyyyzzz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 594); 

                auto tg_yyyyyzzz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 595); 

                auto tg_yyyyyzzz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 596); 

                auto tg_yyyyyzzz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 597); 

                auto tg_yyyyyzzz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 598); 

                auto tg_yyyyyzzz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 599); 

                auto tg_yyyyzzzz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 600); 

                auto tg_yyyyzzzz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 601); 

                auto tg_yyyyzzzz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 602); 

                auto tg_yyyyzzzz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 603); 

                auto tg_yyyyzzzz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 604); 

                auto tg_yyyyzzzz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 605); 

                auto tg_yyyyzzzz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 606); 

                auto tg_yyyyzzzz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 607); 

                auto tg_yyyyzzzz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 608); 

                auto tg_yyyyzzzz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 609); 

                auto tg_yyyyzzzz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 610); 

                auto tg_yyyyzzzz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 611); 

                auto tg_yyyyzzzz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 612); 

                auto tg_yyyyzzzz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 613); 

                auto tg_yyyyzzzz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 614); 

                auto tg_yyyzzzzz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 615); 

                auto tg_yyyzzzzz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 616); 

                auto tg_yyyzzzzz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 617); 

                auto tg_yyyzzzzz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 618); 

                auto tg_yyyzzzzz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 619); 

                auto tg_yyyzzzzz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 620); 

                auto tg_yyyzzzzz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 621); 

                auto tg_yyyzzzzz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 622); 

                auto tg_yyyzzzzz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 623); 

                auto tg_yyyzzzzz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 624); 

                auto tg_yyyzzzzz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 625); 

                auto tg_yyyzzzzz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 626); 

                // Batch of Integrals (579,627)

                #pragma omp simd aligned(fxn, fza, tg_yyyyyyzz_xzzz_0, tg_yyyyyyzz_yyyy_0, tg_yyyyyyzz_yyyz_0, \
                                         tg_yyyyyyzz_yyzz_0, tg_yyyyyyzz_yzzz_0, tg_yyyyyyzz_zzzz_0, tg_yyyyyzz_xzzz_0, \
                                         tg_yyyyyzz_xzzz_1, tg_yyyyyzz_yyy_1, tg_yyyyyzz_yyyy_0, tg_yyyyyzz_yyyy_1, \
                                         tg_yyyyyzz_yyyz_0, tg_yyyyyzz_yyyz_1, tg_yyyyyzz_yyz_1, tg_yyyyyzz_yyzz_0, \
                                         tg_yyyyyzz_yyzz_1, tg_yyyyyzz_yzz_1, tg_yyyyyzz_yzzz_0, tg_yyyyyzz_yzzz_1, \
                                         tg_yyyyyzz_zzz_1, tg_yyyyyzz_zzzz_0, tg_yyyyyzz_zzzz_1, tg_yyyyyzzz_xxxx_0, \
                                         tg_yyyyyzzz_xxxy_0, tg_yyyyyzzz_xxxz_0, tg_yyyyyzzz_xxyy_0, tg_yyyyyzzz_xxyz_0, \
                                         tg_yyyyyzzz_xxzz_0, tg_yyyyyzzz_xyyy_0, tg_yyyyyzzz_xyyz_0, tg_yyyyyzzz_xyzz_0, \
                                         tg_yyyyyzzz_xzzz_0, tg_yyyyyzzz_yyyy_0, tg_yyyyyzzz_yyyz_0, tg_yyyyyzzz_yyzz_0, \
                                         tg_yyyyyzzz_yzzz_0, tg_yyyyyzzz_zzzz_0, tg_yyyyzz_xzzz_0, tg_yyyyzz_xzzz_1, \
                                         tg_yyyyzz_yyyy_0, tg_yyyyzz_yyyy_1, tg_yyyyzz_yyyz_0, tg_yyyyzz_yyyz_1, \
                                         tg_yyyyzz_yyzz_0, tg_yyyyzz_yyzz_1, tg_yyyyzz_yzzz_0, tg_yyyyzz_yzzz_1, \
                                         tg_yyyyzz_zzzz_0, tg_yyyyzz_zzzz_1, tg_yyyyzzz_xxx_1, tg_yyyyzzz_xxxx_0, \
                                         tg_yyyyzzz_xxxx_1, tg_yyyyzzz_xxxy_0, tg_yyyyzzz_xxxy_1, tg_yyyyzzz_xxxz_0, \
                                         tg_yyyyzzz_xxxz_1, tg_yyyyzzz_xxy_1, tg_yyyyzzz_xxyy_0, tg_yyyyzzz_xxyy_1, \
                                         tg_yyyyzzz_xxyz_0, tg_yyyyzzz_xxyz_1, tg_yyyyzzz_xxz_1, tg_yyyyzzz_xxzz_0, \
                                         tg_yyyyzzz_xxzz_1, tg_yyyyzzz_xyy_1, tg_yyyyzzz_xyyy_0, tg_yyyyzzz_xyyy_1, \
                                         tg_yyyyzzz_xyyz_0, tg_yyyyzzz_xyyz_1, tg_yyyyzzz_xyz_1, tg_yyyyzzz_xyzz_0, \
                                         tg_yyyyzzz_xyzz_1, tg_yyyyzzz_xzz_1, tg_yyyyzzz_xzzz_0, tg_yyyyzzz_xzzz_1, \
                                         tg_yyyyzzz_yyy_1, tg_yyyyzzz_yyyy_0, tg_yyyyzzz_yyyy_1, tg_yyyyzzz_yyyz_0, \
                                         tg_yyyyzzz_yyyz_1, tg_yyyyzzz_yyz_1, tg_yyyyzzz_yyzz_0, tg_yyyyzzz_yyzz_1, \
                                         tg_yyyyzzz_yzz_1, tg_yyyyzzz_yzzz_0, tg_yyyyzzz_yzzz_1, tg_yyyyzzz_zzz_1, \
                                         tg_yyyyzzz_zzzz_0, tg_yyyyzzz_zzzz_1, tg_yyyyzzzz_xxxx_0, tg_yyyyzzzz_xxxy_0, \
                                         tg_yyyyzzzz_xxxz_0, tg_yyyyzzzz_xxyy_0, tg_yyyyzzzz_xxyz_0, tg_yyyyzzzz_xxzz_0, \
                                         tg_yyyyzzzz_xyyy_0, tg_yyyyzzzz_xyyz_0, tg_yyyyzzzz_xyzz_0, tg_yyyyzzzz_xzzz_0, \
                                         tg_yyyyzzzz_yyyy_0, tg_yyyyzzzz_yyyz_0, tg_yyyyzzzz_yyzz_0, tg_yyyyzzzz_yzzz_0, \
                                         tg_yyyyzzzz_zzzz_0, tg_yyyzzz_xxxx_0, tg_yyyzzz_xxxx_1, tg_yyyzzz_xxxy_0, \
                                         tg_yyyzzz_xxxy_1, tg_yyyzzz_xxxz_0, tg_yyyzzz_xxxz_1, tg_yyyzzz_xxyy_0, \
                                         tg_yyyzzz_xxyy_1, tg_yyyzzz_xxyz_0, tg_yyyzzz_xxyz_1, tg_yyyzzz_xxzz_0, \
                                         tg_yyyzzz_xxzz_1, tg_yyyzzz_xyyy_0, tg_yyyzzz_xyyy_1, tg_yyyzzz_xyyz_0, \
                                         tg_yyyzzz_xyyz_1, tg_yyyzzz_xyzz_0, tg_yyyzzz_xyzz_1, tg_yyyzzz_xzzz_0, \
                                         tg_yyyzzz_xzzz_1, tg_yyyzzz_yyyy_0, tg_yyyzzz_yyyy_1, tg_yyyzzz_yyyz_0, \
                                         tg_yyyzzz_yyyz_1, tg_yyyzzz_yyzz_0, tg_yyyzzz_yyzz_1, tg_yyyzzz_yzzz_0, \
                                         tg_yyyzzz_yzzz_1, tg_yyyzzz_zzzz_0, tg_yyyzzz_zzzz_1, tg_yyyzzzz_xxx_1, \
                                         tg_yyyzzzz_xxxx_0, tg_yyyzzzz_xxxx_1, tg_yyyzzzz_xxxy_0, tg_yyyzzzz_xxxy_1, \
                                         tg_yyyzzzz_xxxz_0, tg_yyyzzzz_xxxz_1, tg_yyyzzzz_xxy_1, tg_yyyzzzz_xxyy_0, \
                                         tg_yyyzzzz_xxyy_1, tg_yyyzzzz_xxyz_0, tg_yyyzzzz_xxyz_1, tg_yyyzzzz_xxz_1, \
                                         tg_yyyzzzz_xxzz_0, tg_yyyzzzz_xxzz_1, tg_yyyzzzz_xyy_1, tg_yyyzzzz_xyyy_0, \
                                         tg_yyyzzzz_xyyy_1, tg_yyyzzzz_xyyz_0, tg_yyyzzzz_xyyz_1, tg_yyyzzzz_xyz_1, \
                                         tg_yyyzzzz_xyzz_0, tg_yyyzzzz_xyzz_1, tg_yyyzzzz_xzz_1, tg_yyyzzzz_xzzz_0, \
                                         tg_yyyzzzz_xzzz_1, tg_yyyzzzz_yyy_1, tg_yyyzzzz_yyyy_0, tg_yyyzzzz_yyyy_1, \
                                         tg_yyyzzzz_yyyz_0, tg_yyyzzzz_yyyz_1, tg_yyyzzzz_yyz_1, tg_yyyzzzz_yyzz_0, \
                                         tg_yyyzzzz_yyzz_1, tg_yyyzzzz_yzz_1, tg_yyyzzzz_yzzz_0, tg_yyyzzzz_yzzz_1, \
                                         tg_yyyzzzz_zzz_1, tg_yyyzzzz_zzzz_0, tg_yyyzzzz_zzzz_1, tg_yyyzzzzz_xxxx_0, \
                                         tg_yyyzzzzz_xxxy_0, tg_yyyzzzzz_xxxz_0, tg_yyyzzzzz_xxyy_0, tg_yyyzzzzz_xxyz_0, \
                                         tg_yyyzzzzz_xxzz_0, tg_yyyzzzzz_xyyy_0, tg_yyyzzzzz_xyyz_0, tg_yyyzzzzz_xyzz_0, \
                                         tg_yyyzzzzz_xzzz_0, tg_yyyzzzzz_yyyy_0, tg_yyyzzzzz_yyyz_0, tg_yyzzzz_xxxx_0, \
                                         tg_yyzzzz_xxxx_1, tg_yyzzzz_xxxy_0, tg_yyzzzz_xxxy_1, tg_yyzzzz_xxxz_0, \
                                         tg_yyzzzz_xxxz_1, tg_yyzzzz_xxyy_0, tg_yyzzzz_xxyy_1, tg_yyzzzz_xxyz_0, \
                                         tg_yyzzzz_xxyz_1, tg_yyzzzz_xxzz_0, tg_yyzzzz_xxzz_1, tg_yyzzzz_xyyy_0, \
                                         tg_yyzzzz_xyyy_1, tg_yyzzzz_xyyz_0, tg_yyzzzz_xyyz_1, tg_yyzzzz_xyzz_0, \
                                         tg_yyzzzz_xyzz_1, tg_yyzzzz_xzzz_0, tg_yyzzzz_xzzz_1, tg_yyzzzz_yyyy_0, \
                                         tg_yyzzzz_yyyy_1, tg_yyzzzz_yyyz_0, tg_yyzzzz_yyyz_1, tg_yyzzzz_yyzz_0, \
                                         tg_yyzzzz_yyzz_1, tg_yyzzzz_yzzz_0, tg_yyzzzz_yzzz_1, tg_yyzzzz_zzzz_0, \
                                         tg_yyzzzz_zzzz_1, tg_yyzzzzz_xxx_1, tg_yyzzzzz_xxxx_0, tg_yyzzzzz_xxxx_1, \
                                         tg_yyzzzzz_xxxy_0, tg_yyzzzzz_xxxy_1, tg_yyzzzzz_xxxz_0, tg_yyzzzzz_xxxz_1, \
                                         tg_yyzzzzz_xxy_1, tg_yyzzzzz_xxyy_0, tg_yyzzzzz_xxyy_1, tg_yyzzzzz_xxyz_0, \
                                         tg_yyzzzzz_xxyz_1, tg_yyzzzzz_xxz_1, tg_yyzzzzz_xxzz_0, tg_yyzzzzz_xxzz_1, \
                                         tg_yyzzzzz_xyy_1, tg_yyzzzzz_xyyy_0, tg_yyzzzzz_xyyy_1, tg_yyzzzzz_xyyz_0, \
                                         tg_yyzzzzz_xyyz_1, tg_yyzzzzz_xyz_1, tg_yyzzzzz_xyzz_0, tg_yyzzzzz_xyzz_1, \
                                         tg_yyzzzzz_xzz_1, tg_yyzzzzz_xzzz_0, tg_yyzzzzz_xzzz_1, tg_yyzzzzz_yyy_1, \
                                         tg_yyzzzzz_yyyy_0, tg_yyzzzzz_yyyy_1, tg_yyzzzzz_yyyz_0, tg_yyzzzzz_yyyz_1, \
                                         tg_yyzzzzz_yyz_1, tg_yzzzzz_xxxx_0, tg_yzzzzz_xxxx_1, tg_yzzzzz_xxxy_0, \
                                         tg_yzzzzz_xxxy_1, tg_yzzzzz_xxxz_0, tg_yzzzzz_xxxz_1, tg_yzzzzz_xxyy_0, \
                                         tg_yzzzzz_xxyy_1, tg_yzzzzz_xxyz_0, tg_yzzzzz_xxyz_1, tg_yzzzzz_xxzz_0, \
                                         tg_yzzzzz_xxzz_1, tg_yzzzzz_xyyy_0, tg_yzzzzz_xyyy_1, tg_yzzzzz_xyyz_0, \
                                         tg_yzzzzz_xyyz_1, tg_yzzzzz_xyzz_0, tg_yzzzzz_xyzz_1, tg_yzzzzz_xzzz_0, \
                                         tg_yzzzzz_xzzz_1, tg_yzzzzz_yyyy_0, tg_yzzzzz_yyyy_1, tg_yzzzzz_yyyz_0, \
                                         tg_yzzzzz_yyyz_1, wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_yyyyyyzz_xzzz_0[j] = pb_y * tg_yyyyyzz_xzzz_0[j] + wp_y[j] * tg_yyyyyzz_xzzz_1[j] + 2.5 * fl1_fx * tg_yyyyzz_xzzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyyzz_xzzz_1[j];

                    tg_yyyyyyzz_yyyy_0[j] = pb_y * tg_yyyyyzz_yyyy_0[j] + wp_y[j] * tg_yyyyyzz_yyyy_1[j] + 2.5 * fl1_fx * tg_yyyyzz_yyyy_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyyzz_yyyy_1[j] + 2.0 * fl1_fxn * tg_yyyyyzz_yyy_1[j];

                    tg_yyyyyyzz_yyyz_0[j] = pb_y * tg_yyyyyzz_yyyz_0[j] + wp_y[j] * tg_yyyyyzz_yyyz_1[j] + 2.5 * fl1_fx * tg_yyyyzz_yyyz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyyzz_yyyz_1[j] + 1.5 * fl1_fxn * tg_yyyyyzz_yyz_1[j];

                    tg_yyyyyyzz_yyzz_0[j] = pb_y * tg_yyyyyzz_yyzz_0[j] + wp_y[j] * tg_yyyyyzz_yyzz_1[j] + 2.5 * fl1_fx * tg_yyyyzz_yyzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyyzz_yyzz_1[j] + fl1_fxn * tg_yyyyyzz_yzz_1[j];

                    tg_yyyyyyzz_yzzz_0[j] = pb_y * tg_yyyyyzz_yzzz_0[j] + wp_y[j] * tg_yyyyyzz_yzzz_1[j] + 2.5 * fl1_fx * tg_yyyyzz_yzzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyyzz_yzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyzz_zzz_1[j];

                    tg_yyyyyyzz_zzzz_0[j] = pb_y * tg_yyyyyzz_zzzz_0[j] + wp_y[j] * tg_yyyyyzz_zzzz_1[j] + 2.5 * fl1_fx * tg_yyyyzz_zzzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyyzz_zzzz_1[j];

                    tg_yyyyyzzz_xxxx_0[j] = pb_y * tg_yyyyzzz_xxxx_0[j] + wp_y[j] * tg_yyyyzzz_xxxx_1[j] + 2.0 * fl1_fx * tg_yyyzzz_xxxx_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyzzz_xxxx_1[j];

                    tg_yyyyyzzz_xxxy_0[j] = pb_y * tg_yyyyzzz_xxxy_0[j] + wp_y[j] * tg_yyyyzzz_xxxy_1[j] + 2.0 * fl1_fx * tg_yyyzzz_xxxy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyzzz_xxxy_1[j] + 0.5 * fl1_fxn * tg_yyyyzzz_xxx_1[j];

                    tg_yyyyyzzz_xxxz_0[j] = pb_y * tg_yyyyzzz_xxxz_0[j] + wp_y[j] * tg_yyyyzzz_xxxz_1[j] + 2.0 * fl1_fx * tg_yyyzzz_xxxz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyzzz_xxxz_1[j];

                    tg_yyyyyzzz_xxyy_0[j] = pb_y * tg_yyyyzzz_xxyy_0[j] + wp_y[j] * tg_yyyyzzz_xxyy_1[j] + 2.0 * fl1_fx * tg_yyyzzz_xxyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyzzz_xxyy_1[j] + fl1_fxn * tg_yyyyzzz_xxy_1[j];

                    tg_yyyyyzzz_xxyz_0[j] = pb_y * tg_yyyyzzz_xxyz_0[j] + wp_y[j] * tg_yyyyzzz_xxyz_1[j] + 2.0 * fl1_fx * tg_yyyzzz_xxyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyzzz_xxyz_1[j] + 0.5 * fl1_fxn * tg_yyyyzzz_xxz_1[j];

                    tg_yyyyyzzz_xxzz_0[j] = pb_y * tg_yyyyzzz_xxzz_0[j] + wp_y[j] * tg_yyyyzzz_xxzz_1[j] + 2.0 * fl1_fx * tg_yyyzzz_xxzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyzzz_xxzz_1[j];

                    tg_yyyyyzzz_xyyy_0[j] = pb_y * tg_yyyyzzz_xyyy_0[j] + wp_y[j] * tg_yyyyzzz_xyyy_1[j] + 2.0 * fl1_fx * tg_yyyzzz_xyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyzzz_xyyy_1[j] + 1.5 * fl1_fxn * tg_yyyyzzz_xyy_1[j];

                    tg_yyyyyzzz_xyyz_0[j] = pb_y * tg_yyyyzzz_xyyz_0[j] + wp_y[j] * tg_yyyyzzz_xyyz_1[j] + 2.0 * fl1_fx * tg_yyyzzz_xyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyzzz_xyyz_1[j] + fl1_fxn * tg_yyyyzzz_xyz_1[j];

                    tg_yyyyyzzz_xyzz_0[j] = pb_y * tg_yyyyzzz_xyzz_0[j] + wp_y[j] * tg_yyyyzzz_xyzz_1[j] + 2.0 * fl1_fx * tg_yyyzzz_xyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyyyzzz_xzz_1[j];

                    tg_yyyyyzzz_xzzz_0[j] = pb_y * tg_yyyyzzz_xzzz_0[j] + wp_y[j] * tg_yyyyzzz_xzzz_1[j] + 2.0 * fl1_fx * tg_yyyzzz_xzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyzzz_xzzz_1[j];

                    tg_yyyyyzzz_yyyy_0[j] = pb_y * tg_yyyyzzz_yyyy_0[j] + wp_y[j] * tg_yyyyzzz_yyyy_1[j] + 2.0 * fl1_fx * tg_yyyzzz_yyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyzzz_yyyy_1[j] + 2.0 * fl1_fxn * tg_yyyyzzz_yyy_1[j];

                    tg_yyyyyzzz_yyyz_0[j] = pb_y * tg_yyyyzzz_yyyz_0[j] + wp_y[j] * tg_yyyyzzz_yyyz_1[j] + 2.0 * fl1_fx * tg_yyyzzz_yyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyzzz_yyyz_1[j] + 1.5 * fl1_fxn * tg_yyyyzzz_yyz_1[j];

                    tg_yyyyyzzz_yyzz_0[j] = pb_y * tg_yyyyzzz_yyzz_0[j] + wp_y[j] * tg_yyyyzzz_yyzz_1[j] + 2.0 * fl1_fx * tg_yyyzzz_yyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyzzz_yyzz_1[j] + fl1_fxn * tg_yyyyzzz_yzz_1[j];

                    tg_yyyyyzzz_yzzz_0[j] = pb_y * tg_yyyyzzz_yzzz_0[j] + wp_y[j] * tg_yyyyzzz_yzzz_1[j] + 2.0 * fl1_fx * tg_yyyzzz_yzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyzzz_yzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyzzz_zzz_1[j];

                    tg_yyyyyzzz_zzzz_0[j] = pb_y * tg_yyyyzzz_zzzz_0[j] + wp_y[j] * tg_yyyyzzz_zzzz_1[j] + 2.0 * fl1_fx * tg_yyyzzz_zzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyzzz_zzzz_1[j];

                    tg_yyyyzzzz_xxxx_0[j] = pb_y * tg_yyyzzzz_xxxx_0[j] + wp_y[j] * tg_yyyzzzz_xxxx_1[j] + 1.5 * fl1_fx * tg_yyzzzz_xxxx_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzzzz_xxxx_1[j];

                    tg_yyyyzzzz_xxxy_0[j] = pb_y * tg_yyyzzzz_xxxy_0[j] + wp_y[j] * tg_yyyzzzz_xxxy_1[j] + 1.5 * fl1_fx * tg_yyzzzz_xxxy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzzzz_xxxy_1[j] + 0.5 * fl1_fxn * tg_yyyzzzz_xxx_1[j];

                    tg_yyyyzzzz_xxxz_0[j] = pb_y * tg_yyyzzzz_xxxz_0[j] + wp_y[j] * tg_yyyzzzz_xxxz_1[j] + 1.5 * fl1_fx * tg_yyzzzz_xxxz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzzzz_xxxz_1[j];

                    tg_yyyyzzzz_xxyy_0[j] = pb_y * tg_yyyzzzz_xxyy_0[j] + wp_y[j] * tg_yyyzzzz_xxyy_1[j] + 1.5 * fl1_fx * tg_yyzzzz_xxyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzzzz_xxyy_1[j] + fl1_fxn * tg_yyyzzzz_xxy_1[j];

                    tg_yyyyzzzz_xxyz_0[j] = pb_y * tg_yyyzzzz_xxyz_0[j] + wp_y[j] * tg_yyyzzzz_xxyz_1[j] + 1.5 * fl1_fx * tg_yyzzzz_xxyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzzzz_xxyz_1[j] + 0.5 * fl1_fxn * tg_yyyzzzz_xxz_1[j];

                    tg_yyyyzzzz_xxzz_0[j] = pb_y * tg_yyyzzzz_xxzz_0[j] + wp_y[j] * tg_yyyzzzz_xxzz_1[j] + 1.5 * fl1_fx * tg_yyzzzz_xxzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzzzz_xxzz_1[j];

                    tg_yyyyzzzz_xyyy_0[j] = pb_y * tg_yyyzzzz_xyyy_0[j] + wp_y[j] * tg_yyyzzzz_xyyy_1[j] + 1.5 * fl1_fx * tg_yyzzzz_xyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzzzz_xyyy_1[j] + 1.5 * fl1_fxn * tg_yyyzzzz_xyy_1[j];

                    tg_yyyyzzzz_xyyz_0[j] = pb_y * tg_yyyzzzz_xyyz_0[j] + wp_y[j] * tg_yyyzzzz_xyyz_1[j] + 1.5 * fl1_fx * tg_yyzzzz_xyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzzzz_xyyz_1[j] + fl1_fxn * tg_yyyzzzz_xyz_1[j];

                    tg_yyyyzzzz_xyzz_0[j] = pb_y * tg_yyyzzzz_xyzz_0[j] + wp_y[j] * tg_yyyzzzz_xyzz_1[j] + 1.5 * fl1_fx * tg_yyzzzz_xyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyyzzzz_xzz_1[j];

                    tg_yyyyzzzz_xzzz_0[j] = pb_y * tg_yyyzzzz_xzzz_0[j] + wp_y[j] * tg_yyyzzzz_xzzz_1[j] + 1.5 * fl1_fx * tg_yyzzzz_xzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzzzz_xzzz_1[j];

                    tg_yyyyzzzz_yyyy_0[j] = pb_y * tg_yyyzzzz_yyyy_0[j] + wp_y[j] * tg_yyyzzzz_yyyy_1[j] + 1.5 * fl1_fx * tg_yyzzzz_yyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzzzz_yyyy_1[j] + 2.0 * fl1_fxn * tg_yyyzzzz_yyy_1[j];

                    tg_yyyyzzzz_yyyz_0[j] = pb_y * tg_yyyzzzz_yyyz_0[j] + wp_y[j] * tg_yyyzzzz_yyyz_1[j] + 1.5 * fl1_fx * tg_yyzzzz_yyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzzzz_yyyz_1[j] + 1.5 * fl1_fxn * tg_yyyzzzz_yyz_1[j];

                    tg_yyyyzzzz_yyzz_0[j] = pb_y * tg_yyyzzzz_yyzz_0[j] + wp_y[j] * tg_yyyzzzz_yyzz_1[j] + 1.5 * fl1_fx * tg_yyzzzz_yyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzzzz_yyzz_1[j] + fl1_fxn * tg_yyyzzzz_yzz_1[j];

                    tg_yyyyzzzz_yzzz_0[j] = pb_y * tg_yyyzzzz_yzzz_0[j] + wp_y[j] * tg_yyyzzzz_yzzz_1[j] + 1.5 * fl1_fx * tg_yyzzzz_yzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzzzz_yzzz_1[j] + 0.5 * fl1_fxn * tg_yyyzzzz_zzz_1[j];

                    tg_yyyyzzzz_zzzz_0[j] = pb_y * tg_yyyzzzz_zzzz_0[j] + wp_y[j] * tg_yyyzzzz_zzzz_1[j] + 1.5 * fl1_fx * tg_yyzzzz_zzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzzzz_zzzz_1[j];

                    tg_yyyzzzzz_xxxx_0[j] = pb_y * tg_yyzzzzz_xxxx_0[j] + wp_y[j] * tg_yyzzzzz_xxxx_1[j] + fl1_fx * tg_yzzzzz_xxxx_0[j] - fl1_fx * fl1_fza * tg_yzzzzz_xxxx_1[j];

                    tg_yyyzzzzz_xxxy_0[j] = pb_y * tg_yyzzzzz_xxxy_0[j] + wp_y[j] * tg_yyzzzzz_xxxy_1[j] + fl1_fx * tg_yzzzzz_xxxy_0[j] - fl1_fx * fl1_fza * tg_yzzzzz_xxxy_1[j] + 0.5 * fl1_fxn * tg_yyzzzzz_xxx_1[j];

                    tg_yyyzzzzz_xxxz_0[j] = pb_y * tg_yyzzzzz_xxxz_0[j] + wp_y[j] * tg_yyzzzzz_xxxz_1[j] + fl1_fx * tg_yzzzzz_xxxz_0[j] - fl1_fx * fl1_fza * tg_yzzzzz_xxxz_1[j];

                    tg_yyyzzzzz_xxyy_0[j] = pb_y * tg_yyzzzzz_xxyy_0[j] + wp_y[j] * tg_yyzzzzz_xxyy_1[j] + fl1_fx * tg_yzzzzz_xxyy_0[j] - fl1_fx * fl1_fza * tg_yzzzzz_xxyy_1[j] + fl1_fxn * tg_yyzzzzz_xxy_1[j];

                    tg_yyyzzzzz_xxyz_0[j] = pb_y * tg_yyzzzzz_xxyz_0[j] + wp_y[j] * tg_yyzzzzz_xxyz_1[j] + fl1_fx * tg_yzzzzz_xxyz_0[j] - fl1_fx * fl1_fza * tg_yzzzzz_xxyz_1[j] + 0.5 * fl1_fxn * tg_yyzzzzz_xxz_1[j];

                    tg_yyyzzzzz_xxzz_0[j] = pb_y * tg_yyzzzzz_xxzz_0[j] + wp_y[j] * tg_yyzzzzz_xxzz_1[j] + fl1_fx * tg_yzzzzz_xxzz_0[j] - fl1_fx * fl1_fza * tg_yzzzzz_xxzz_1[j];

                    tg_yyyzzzzz_xyyy_0[j] = pb_y * tg_yyzzzzz_xyyy_0[j] + wp_y[j] * tg_yyzzzzz_xyyy_1[j] + fl1_fx * tg_yzzzzz_xyyy_0[j] - fl1_fx * fl1_fza * tg_yzzzzz_xyyy_1[j] + 1.5 * fl1_fxn * tg_yyzzzzz_xyy_1[j];

                    tg_yyyzzzzz_xyyz_0[j] = pb_y * tg_yyzzzzz_xyyz_0[j] + wp_y[j] * tg_yyzzzzz_xyyz_1[j] + fl1_fx * tg_yzzzzz_xyyz_0[j] - fl1_fx * fl1_fza * tg_yzzzzz_xyyz_1[j] + fl1_fxn * tg_yyzzzzz_xyz_1[j];

                    tg_yyyzzzzz_xyzz_0[j] = pb_y * tg_yyzzzzz_xyzz_0[j] + wp_y[j] * tg_yyzzzzz_xyzz_1[j] + fl1_fx * tg_yzzzzz_xyzz_0[j] - fl1_fx * fl1_fza * tg_yzzzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyzzzzz_xzz_1[j];

                    tg_yyyzzzzz_xzzz_0[j] = pb_y * tg_yyzzzzz_xzzz_0[j] + wp_y[j] * tg_yyzzzzz_xzzz_1[j] + fl1_fx * tg_yzzzzz_xzzz_0[j] - fl1_fx * fl1_fza * tg_yzzzzz_xzzz_1[j];

                    tg_yyyzzzzz_yyyy_0[j] = pb_y * tg_yyzzzzz_yyyy_0[j] + wp_y[j] * tg_yyzzzzz_yyyy_1[j] + fl1_fx * tg_yzzzzz_yyyy_0[j] - fl1_fx * fl1_fza * tg_yzzzzz_yyyy_1[j] + 2.0 * fl1_fxn * tg_yyzzzzz_yyy_1[j];

                    tg_yyyzzzzz_yyyz_0[j] = pb_y * tg_yyzzzzz_yyyz_0[j] + wp_y[j] * tg_yyzzzzz_yyyz_1[j] + fl1_fx * tg_yzzzzz_yyyz_0[j] - fl1_fx * fl1_fza * tg_yzzzzz_yyyz_1[j] + 1.5 * fl1_fxn * tg_yyzzzzz_yyz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSLSG_627_675(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (627,675)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        auto r_pb_z = braGtoPairsBlock.getDistancesPBZ();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {8, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {8, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_7_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_6_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_6_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_7_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto pb_y = r_pb_y[i];

                auto pb_z = r_pb_z[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_y = wpDistances.data(3 * idx + 1);

                auto wp_z = wpDistances.data(3 * idx + 2);

                // set up pointers to auxilary integrals

                auto tg_yyzzzzz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 507); 

                auto tg_yyzzzzz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 508); 

                auto tg_yyzzzzz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 509); 

                auto tg_yzzzzzz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 510); 

                auto tg_yzzzzzz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 511); 

                auto tg_yzzzzzz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 512); 

                auto tg_yzzzzzz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 513); 

                auto tg_yzzzzzz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 514); 

                auto tg_yzzzzzz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 515); 

                auto tg_yzzzzzz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 516); 

                auto tg_yzzzzzz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 517); 

                auto tg_yzzzzzz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 518); 

                auto tg_yzzzzzz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 519); 

                auto tg_yzzzzzz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 520); 

                auto tg_yzzzzzz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 521); 

                auto tg_yzzzzzz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 522); 

                auto tg_yzzzzzz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 523); 

                auto tg_yzzzzzz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 524); 

                auto tg_zzzzzzz_xxxx_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 525); 

                auto tg_zzzzzzz_xxxy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 526); 

                auto tg_zzzzzzz_xxxz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 527); 

                auto tg_zzzzzzz_xxyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 528); 

                auto tg_zzzzzzz_xxyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 529); 

                auto tg_zzzzzzz_xxzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 530); 

                auto tg_zzzzzzz_xyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 531); 

                auto tg_zzzzzzz_xyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 532); 

                auto tg_zzzzzzz_xyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 533); 

                auto tg_zzzzzzz_xzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 534); 

                auto tg_zzzzzzz_yyyy_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 535); 

                auto tg_zzzzzzz_yyyz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 536); 

                auto tg_zzzzzzz_yyzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 537); 

                auto tg_zzzzzzz_yzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 538); 

                auto tg_zzzzzzz_zzzz_0 = primBuffer.data(pidx_g_7_4_m0 + 540 * idx + 539); 

                auto tg_yyzzzzz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 507); 

                auto tg_yyzzzzz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 508); 

                auto tg_yyzzzzz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 509); 

                auto tg_yzzzzzz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 510); 

                auto tg_yzzzzzz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 511); 

                auto tg_yzzzzzz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 512); 

                auto tg_yzzzzzz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 513); 

                auto tg_yzzzzzz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 514); 

                auto tg_yzzzzzz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 515); 

                auto tg_yzzzzzz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 516); 

                auto tg_yzzzzzz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 517); 

                auto tg_yzzzzzz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 518); 

                auto tg_yzzzzzz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 519); 

                auto tg_yzzzzzz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 520); 

                auto tg_yzzzzzz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 521); 

                auto tg_yzzzzzz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 522); 

                auto tg_yzzzzzz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 523); 

                auto tg_yzzzzzz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 524); 

                auto tg_zzzzzzz_xxxx_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 525); 

                auto tg_zzzzzzz_xxxy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 526); 

                auto tg_zzzzzzz_xxxz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 527); 

                auto tg_zzzzzzz_xxyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 528); 

                auto tg_zzzzzzz_xxyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 529); 

                auto tg_zzzzzzz_xxzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 530); 

                auto tg_zzzzzzz_xyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 531); 

                auto tg_zzzzzzz_xyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 532); 

                auto tg_zzzzzzz_xyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 533); 

                auto tg_zzzzzzz_xzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 534); 

                auto tg_zzzzzzz_yyyy_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 535); 

                auto tg_zzzzzzz_yyyz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 536); 

                auto tg_zzzzzzz_yyzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 537); 

                auto tg_zzzzzzz_yzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 538); 

                auto tg_zzzzzzz_zzzz_1 = primBuffer.data(pidx_g_7_4_m1 + 540 * idx + 539); 

                auto tg_yzzzzz_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 402); 

                auto tg_yzzzzz_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 403); 

                auto tg_yzzzzz_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 404); 

                auto tg_zzzzzz_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 405); 

                auto tg_zzzzzz_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 406); 

                auto tg_zzzzzz_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 407); 

                auto tg_zzzzzz_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 408); 

                auto tg_zzzzzz_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 409); 

                auto tg_zzzzzz_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 410); 

                auto tg_zzzzzz_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 411); 

                auto tg_zzzzzz_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 412); 

                auto tg_zzzzzz_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 413); 

                auto tg_zzzzzz_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 414); 

                auto tg_zzzzzz_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 415); 

                auto tg_zzzzzz_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 416); 

                auto tg_zzzzzz_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 417); 

                auto tg_zzzzzz_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 418); 

                auto tg_zzzzzz_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 419); 

                auto tg_yzzzzz_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 402); 

                auto tg_yzzzzz_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 403); 

                auto tg_yzzzzz_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 404); 

                auto tg_zzzzzz_xxxx_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 405); 

                auto tg_zzzzzz_xxxy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 406); 

                auto tg_zzzzzz_xxxz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 407); 

                auto tg_zzzzzz_xxyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 408); 

                auto tg_zzzzzz_xxyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 409); 

                auto tg_zzzzzz_xxzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 410); 

                auto tg_zzzzzz_xyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 411); 

                auto tg_zzzzzz_xyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 412); 

                auto tg_zzzzzz_xyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 413); 

                auto tg_zzzzzz_xzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 414); 

                auto tg_zzzzzz_yyyy_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 415); 

                auto tg_zzzzzz_yyyz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 416); 

                auto tg_zzzzzz_yyzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 417); 

                auto tg_zzzzzz_yzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 418); 

                auto tg_zzzzzz_zzzz_1 = primBuffer.data(pidx_g_6_4_m1 + 420 * idx + 419); 

                auto tg_yyzzzzz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 338); 

                auto tg_yyzzzzz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 339); 

                auto tg_yzzzzzz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 340); 

                auto tg_yzzzzzz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 341); 

                auto tg_yzzzzzz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 342); 

                auto tg_yzzzzzz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 343); 

                auto tg_yzzzzzz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 344); 

                auto tg_yzzzzzz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 345); 

                auto tg_yzzzzzz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 346); 

                auto tg_yzzzzzz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 347); 

                auto tg_yzzzzzz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 348); 

                auto tg_yzzzzzz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 349); 

                auto tg_zzzzzzz_xxx_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 350); 

                auto tg_zzzzzzz_xxy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 351); 

                auto tg_zzzzzzz_xxz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 352); 

                auto tg_zzzzzzz_xyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 353); 

                auto tg_zzzzzzz_xyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 354); 

                auto tg_zzzzzzz_xzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 355); 

                auto tg_zzzzzzz_yyy_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 356); 

                auto tg_zzzzzzz_yyz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 357); 

                auto tg_zzzzzzz_yzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 358); 

                auto tg_zzzzzzz_zzz_1 = primBuffer.data(pidx_g_7_3_m1 + 360 * idx + 359); 

                // set up pointers to integrals

                auto tg_yyyzzzzz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 627); 

                auto tg_yyyzzzzz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 628); 

                auto tg_yyyzzzzz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 629); 

                auto tg_yyzzzzzz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 630); 

                auto tg_yyzzzzzz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 631); 

                auto tg_yyzzzzzz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 632); 

                auto tg_yyzzzzzz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 633); 

                auto tg_yyzzzzzz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 634); 

                auto tg_yyzzzzzz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 635); 

                auto tg_yyzzzzzz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 636); 

                auto tg_yyzzzzzz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 637); 

                auto tg_yyzzzzzz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 638); 

                auto tg_yyzzzzzz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 639); 

                auto tg_yyzzzzzz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 640); 

                auto tg_yyzzzzzz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 641); 

                auto tg_yyzzzzzz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 642); 

                auto tg_yyzzzzzz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 643); 

                auto tg_yyzzzzzz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 644); 

                auto tg_yzzzzzzz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 645); 

                auto tg_yzzzzzzz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 646); 

                auto tg_yzzzzzzz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 647); 

                auto tg_yzzzzzzz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 648); 

                auto tg_yzzzzzzz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 649); 

                auto tg_yzzzzzzz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 650); 

                auto tg_yzzzzzzz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 651); 

                auto tg_yzzzzzzz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 652); 

                auto tg_yzzzzzzz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 653); 

                auto tg_yzzzzzzz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 654); 

                auto tg_yzzzzzzz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 655); 

                auto tg_yzzzzzzz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 656); 

                auto tg_yzzzzzzz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 657); 

                auto tg_yzzzzzzz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 658); 

                auto tg_yzzzzzzz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 659); 

                auto tg_zzzzzzzz_xxxx_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 660); 

                auto tg_zzzzzzzz_xxxy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 661); 

                auto tg_zzzzzzzz_xxxz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 662); 

                auto tg_zzzzzzzz_xxyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 663); 

                auto tg_zzzzzzzz_xxyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 664); 

                auto tg_zzzzzzzz_xxzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 665); 

                auto tg_zzzzzzzz_xyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 666); 

                auto tg_zzzzzzzz_xyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 667); 

                auto tg_zzzzzzzz_xyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 668); 

                auto tg_zzzzzzzz_xzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 669); 

                auto tg_zzzzzzzz_yyyy_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 670); 

                auto tg_zzzzzzzz_yyyz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 671); 

                auto tg_zzzzzzzz_yyzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 672); 

                auto tg_zzzzzzzz_yzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 673); 

                auto tg_zzzzzzzz_zzzz_0 = primBuffer.data(pidx_g_8_4_m0 + 675 * idx + 674); 

                // Batch of Integrals (627,675)

                #pragma omp simd aligned(fxn, fza, tg_yyyzzzzz_yyzz_0, tg_yyyzzzzz_yzzz_0, tg_yyyzzzzz_zzzz_0, \
                                         tg_yyzzzzz_yyzz_0, tg_yyzzzzz_yyzz_1, tg_yyzzzzz_yzz_1, tg_yyzzzzz_yzzz_0, \
                                         tg_yyzzzzz_yzzz_1, tg_yyzzzzz_zzz_1, tg_yyzzzzz_zzzz_0, tg_yyzzzzz_zzzz_1, \
                                         tg_yyzzzzzz_xxxx_0, tg_yyzzzzzz_xxxy_0, tg_yyzzzzzz_xxxz_0, tg_yyzzzzzz_xxyy_0, \
                                         tg_yyzzzzzz_xxyz_0, tg_yyzzzzzz_xxzz_0, tg_yyzzzzzz_xyyy_0, tg_yyzzzzzz_xyyz_0, \
                                         tg_yyzzzzzz_xyzz_0, tg_yyzzzzzz_xzzz_0, tg_yyzzzzzz_yyyy_0, tg_yyzzzzzz_yyyz_0, \
                                         tg_yyzzzzzz_yyzz_0, tg_yyzzzzzz_yzzz_0, tg_yyzzzzzz_zzzz_0, tg_yzzzzz_yyzz_0, \
                                         tg_yzzzzz_yyzz_1, tg_yzzzzz_yzzz_0, tg_yzzzzz_yzzz_1, tg_yzzzzz_zzzz_0, \
                                         tg_yzzzzz_zzzz_1, tg_yzzzzzz_xxx_1, tg_yzzzzzz_xxxx_0, tg_yzzzzzz_xxxx_1, \
                                         tg_yzzzzzz_xxxy_0, tg_yzzzzzz_xxxy_1, tg_yzzzzzz_xxxz_0, tg_yzzzzzz_xxxz_1, \
                                         tg_yzzzzzz_xxy_1, tg_yzzzzzz_xxyy_0, tg_yzzzzzz_xxyy_1, tg_yzzzzzz_xxyz_0, \
                                         tg_yzzzzzz_xxyz_1, tg_yzzzzzz_xxz_1, tg_yzzzzzz_xxzz_0, tg_yzzzzzz_xxzz_1, \
                                         tg_yzzzzzz_xyy_1, tg_yzzzzzz_xyyy_0, tg_yzzzzzz_xyyy_1, tg_yzzzzzz_xyyz_0, \
                                         tg_yzzzzzz_xyyz_1, tg_yzzzzzz_xyz_1, tg_yzzzzzz_xyzz_0, tg_yzzzzzz_xyzz_1, \
                                         tg_yzzzzzz_xzz_1, tg_yzzzzzz_xzzz_0, tg_yzzzzzz_xzzz_1, tg_yzzzzzz_yyy_1, \
                                         tg_yzzzzzz_yyyy_0, tg_yzzzzzz_yyyy_1, tg_yzzzzzz_yyyz_0, tg_yzzzzzz_yyyz_1, \
                                         tg_yzzzzzz_yyz_1, tg_yzzzzzz_yyzz_0, tg_yzzzzzz_yyzz_1, tg_yzzzzzz_yzz_1, \
                                         tg_yzzzzzz_yzzz_0, tg_yzzzzzz_yzzz_1, tg_yzzzzzz_zzz_1, tg_yzzzzzz_zzzz_0, \
                                         tg_yzzzzzz_zzzz_1, tg_yzzzzzzz_xxxx_0, tg_yzzzzzzz_xxxy_0, tg_yzzzzzzz_xxxz_0, \
                                         tg_yzzzzzzz_xxyy_0, tg_yzzzzzzz_xxyz_0, tg_yzzzzzzz_xxzz_0, tg_yzzzzzzz_xyyy_0, \
                                         tg_yzzzzzzz_xyyz_0, tg_yzzzzzzz_xyzz_0, tg_yzzzzzzz_xzzz_0, tg_yzzzzzzz_yyyy_0, \
                                         tg_yzzzzzzz_yyyz_0, tg_yzzzzzzz_yyzz_0, tg_yzzzzzzz_yzzz_0, tg_yzzzzzzz_zzzz_0, \
                                         tg_zzzzzz_xxxx_0, tg_zzzzzz_xxxx_1, tg_zzzzzz_xxxy_0, tg_zzzzzz_xxxy_1, \
                                         tg_zzzzzz_xxxz_0, tg_zzzzzz_xxxz_1, tg_zzzzzz_xxyy_0, tg_zzzzzz_xxyy_1, \
                                         tg_zzzzzz_xxyz_0, tg_zzzzzz_xxyz_1, tg_zzzzzz_xxzz_0, tg_zzzzzz_xxzz_1, \
                                         tg_zzzzzz_xyyy_0, tg_zzzzzz_xyyy_1, tg_zzzzzz_xyyz_0, tg_zzzzzz_xyyz_1, \
                                         tg_zzzzzz_xyzz_0, tg_zzzzzz_xyzz_1, tg_zzzzzz_xzzz_0, tg_zzzzzz_xzzz_1, \
                                         tg_zzzzzz_yyyy_0, tg_zzzzzz_yyyy_1, tg_zzzzzz_yyyz_0, tg_zzzzzz_yyyz_1, \
                                         tg_zzzzzz_yyzz_0, tg_zzzzzz_yyzz_1, tg_zzzzzz_yzzz_0, tg_zzzzzz_yzzz_1, \
                                         tg_zzzzzz_zzzz_0, tg_zzzzzz_zzzz_1, tg_zzzzzzz_xxx_1, tg_zzzzzzz_xxxx_0, \
                                         tg_zzzzzzz_xxxx_1, tg_zzzzzzz_xxxy_0, tg_zzzzzzz_xxxy_1, tg_zzzzzzz_xxxz_0, \
                                         tg_zzzzzzz_xxxz_1, tg_zzzzzzz_xxy_1, tg_zzzzzzz_xxyy_0, tg_zzzzzzz_xxyy_1, \
                                         tg_zzzzzzz_xxyz_0, tg_zzzzzzz_xxyz_1, tg_zzzzzzz_xxz_1, tg_zzzzzzz_xxzz_0, \
                                         tg_zzzzzzz_xxzz_1, tg_zzzzzzz_xyy_1, tg_zzzzzzz_xyyy_0, tg_zzzzzzz_xyyy_1, \
                                         tg_zzzzzzz_xyyz_0, tg_zzzzzzz_xyyz_1, tg_zzzzzzz_xyz_1, tg_zzzzzzz_xyzz_0, \
                                         tg_zzzzzzz_xyzz_1, tg_zzzzzzz_xzz_1, tg_zzzzzzz_xzzz_0, tg_zzzzzzz_xzzz_1, \
                                         tg_zzzzzzz_yyy_1, tg_zzzzzzz_yyyy_0, tg_zzzzzzz_yyyy_1, tg_zzzzzzz_yyyz_0, \
                                         tg_zzzzzzz_yyyz_1, tg_zzzzzzz_yyz_1, tg_zzzzzzz_yyzz_0, tg_zzzzzzz_yyzz_1, \
                                         tg_zzzzzzz_yzz_1, tg_zzzzzzz_yzzz_0, tg_zzzzzzz_yzzz_1, tg_zzzzzzz_zzz_1, \
                                         tg_zzzzzzz_zzzz_0, tg_zzzzzzz_zzzz_1, tg_zzzzzzzz_xxxx_0, tg_zzzzzzzz_xxxy_0, \
                                         tg_zzzzzzzz_xxxz_0, tg_zzzzzzzz_xxyy_0, tg_zzzzzzzz_xxyz_0, tg_zzzzzzzz_xxzz_0, \
                                         tg_zzzzzzzz_xyyy_0, tg_zzzzzzzz_xyyz_0, tg_zzzzzzzz_xyzz_0, tg_zzzzzzzz_xzzz_0, \
                                         tg_zzzzzzzz_yyyy_0, tg_zzzzzzzz_yyyz_0, tg_zzzzzzzz_yyzz_0, tg_zzzzzzzz_yzzz_0, \
                                         tg_zzzzzzzz_zzzz_0, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_yyyzzzzz_yyzz_0[j] = pb_y * tg_yyzzzzz_yyzz_0[j] + wp_y[j] * tg_yyzzzzz_yyzz_1[j] + fl1_fx * tg_yzzzzz_yyzz_0[j] - fl1_fx * fl1_fza * tg_yzzzzz_yyzz_1[j] + fl1_fxn * tg_yyzzzzz_yzz_1[j];

                    tg_yyyzzzzz_yzzz_0[j] = pb_y * tg_yyzzzzz_yzzz_0[j] + wp_y[j] * tg_yyzzzzz_yzzz_1[j] + fl1_fx * tg_yzzzzz_yzzz_0[j] - fl1_fx * fl1_fza * tg_yzzzzz_yzzz_1[j] + 0.5 * fl1_fxn * tg_yyzzzzz_zzz_1[j];

                    tg_yyyzzzzz_zzzz_0[j] = pb_y * tg_yyzzzzz_zzzz_0[j] + wp_y[j] * tg_yyzzzzz_zzzz_1[j] + fl1_fx * tg_yzzzzz_zzzz_0[j] - fl1_fx * fl1_fza * tg_yzzzzz_zzzz_1[j];

                    tg_yyzzzzzz_xxxx_0[j] = pb_y * tg_yzzzzzz_xxxx_0[j] + wp_y[j] * tg_yzzzzzz_xxxx_1[j] + 0.5 * fl1_fx * tg_zzzzzz_xxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_xxxx_1[j];

                    tg_yyzzzzzz_xxxy_0[j] = pb_y * tg_yzzzzzz_xxxy_0[j] + wp_y[j] * tg_yzzzzzz_xxxy_1[j] + 0.5 * fl1_fx * tg_zzzzzz_xxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_xxxy_1[j] + 0.5 * fl1_fxn * tg_yzzzzzz_xxx_1[j];

                    tg_yyzzzzzz_xxxz_0[j] = pb_y * tg_yzzzzzz_xxxz_0[j] + wp_y[j] * tg_yzzzzzz_xxxz_1[j] + 0.5 * fl1_fx * tg_zzzzzz_xxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_xxxz_1[j];

                    tg_yyzzzzzz_xxyy_0[j] = pb_y * tg_yzzzzzz_xxyy_0[j] + wp_y[j] * tg_yzzzzzz_xxyy_1[j] + 0.5 * fl1_fx * tg_zzzzzz_xxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_xxyy_1[j] + fl1_fxn * tg_yzzzzzz_xxy_1[j];

                    tg_yyzzzzzz_xxyz_0[j] = pb_y * tg_yzzzzzz_xxyz_0[j] + wp_y[j] * tg_yzzzzzz_xxyz_1[j] + 0.5 * fl1_fx * tg_zzzzzz_xxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_xxyz_1[j] + 0.5 * fl1_fxn * tg_yzzzzzz_xxz_1[j];

                    tg_yyzzzzzz_xxzz_0[j] = pb_y * tg_yzzzzzz_xxzz_0[j] + wp_y[j] * tg_yzzzzzz_xxzz_1[j] + 0.5 * fl1_fx * tg_zzzzzz_xxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_xxzz_1[j];

                    tg_yyzzzzzz_xyyy_0[j] = pb_y * tg_yzzzzzz_xyyy_0[j] + wp_y[j] * tg_yzzzzzz_xyyy_1[j] + 0.5 * fl1_fx * tg_zzzzzz_xyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_xyyy_1[j] + 1.5 * fl1_fxn * tg_yzzzzzz_xyy_1[j];

                    tg_yyzzzzzz_xyyz_0[j] = pb_y * tg_yzzzzzz_xyyz_0[j] + wp_y[j] * tg_yzzzzzz_xyyz_1[j] + 0.5 * fl1_fx * tg_zzzzzz_xyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_xyyz_1[j] + fl1_fxn * tg_yzzzzzz_xyz_1[j];

                    tg_yyzzzzzz_xyzz_0[j] = pb_y * tg_yzzzzzz_xyzz_0[j] + wp_y[j] * tg_yzzzzzz_xyzz_1[j] + 0.5 * fl1_fx * tg_zzzzzz_xyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yzzzzzz_xzz_1[j];

                    tg_yyzzzzzz_xzzz_0[j] = pb_y * tg_yzzzzzz_xzzz_0[j] + wp_y[j] * tg_yzzzzzz_xzzz_1[j] + 0.5 * fl1_fx * tg_zzzzzz_xzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_xzzz_1[j];

                    tg_yyzzzzzz_yyyy_0[j] = pb_y * tg_yzzzzzz_yyyy_0[j] + wp_y[j] * tg_yzzzzzz_yyyy_1[j] + 0.5 * fl1_fx * tg_zzzzzz_yyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_yyyy_1[j] + 2.0 * fl1_fxn * tg_yzzzzzz_yyy_1[j];

                    tg_yyzzzzzz_yyyz_0[j] = pb_y * tg_yzzzzzz_yyyz_0[j] + wp_y[j] * tg_yzzzzzz_yyyz_1[j] + 0.5 * fl1_fx * tg_zzzzzz_yyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_yyyz_1[j] + 1.5 * fl1_fxn * tg_yzzzzzz_yyz_1[j];

                    tg_yyzzzzzz_yyzz_0[j] = pb_y * tg_yzzzzzz_yyzz_0[j] + wp_y[j] * tg_yzzzzzz_yyzz_1[j] + 0.5 * fl1_fx * tg_zzzzzz_yyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_yyzz_1[j] + fl1_fxn * tg_yzzzzzz_yzz_1[j];

                    tg_yyzzzzzz_yzzz_0[j] = pb_y * tg_yzzzzzz_yzzz_0[j] + wp_y[j] * tg_yzzzzzz_yzzz_1[j] + 0.5 * fl1_fx * tg_zzzzzz_yzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_yzzz_1[j] + 0.5 * fl1_fxn * tg_yzzzzzz_zzz_1[j];

                    tg_yyzzzzzz_zzzz_0[j] = pb_y * tg_yzzzzzz_zzzz_0[j] + wp_y[j] * tg_yzzzzzz_zzzz_1[j] + 0.5 * fl1_fx * tg_zzzzzz_zzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_zzzz_1[j];

                    tg_yzzzzzzz_xxxx_0[j] = pb_y * tg_zzzzzzz_xxxx_0[j] + wp_y[j] * tg_zzzzzzz_xxxx_1[j];

                    tg_yzzzzzzz_xxxy_0[j] = pb_y * tg_zzzzzzz_xxxy_0[j] + wp_y[j] * tg_zzzzzzz_xxxy_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_xxx_1[j];

                    tg_yzzzzzzz_xxxz_0[j] = pb_y * tg_zzzzzzz_xxxz_0[j] + wp_y[j] * tg_zzzzzzz_xxxz_1[j];

                    tg_yzzzzzzz_xxyy_0[j] = pb_y * tg_zzzzzzz_xxyy_0[j] + wp_y[j] * tg_zzzzzzz_xxyy_1[j] + fl1_fxn * tg_zzzzzzz_xxy_1[j];

                    tg_yzzzzzzz_xxyz_0[j] = pb_y * tg_zzzzzzz_xxyz_0[j] + wp_y[j] * tg_zzzzzzz_xxyz_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_xxz_1[j];

                    tg_yzzzzzzz_xxzz_0[j] = pb_y * tg_zzzzzzz_xxzz_0[j] + wp_y[j] * tg_zzzzzzz_xxzz_1[j];

                    tg_yzzzzzzz_xyyy_0[j] = pb_y * tg_zzzzzzz_xyyy_0[j] + wp_y[j] * tg_zzzzzzz_xyyy_1[j] + 1.5 * fl1_fxn * tg_zzzzzzz_xyy_1[j];

                    tg_yzzzzzzz_xyyz_0[j] = pb_y * tg_zzzzzzz_xyyz_0[j] + wp_y[j] * tg_zzzzzzz_xyyz_1[j] + fl1_fxn * tg_zzzzzzz_xyz_1[j];

                    tg_yzzzzzzz_xyzz_0[j] = pb_y * tg_zzzzzzz_xyzz_0[j] + wp_y[j] * tg_zzzzzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_xzz_1[j];

                    tg_yzzzzzzz_xzzz_0[j] = pb_y * tg_zzzzzzz_xzzz_0[j] + wp_y[j] * tg_zzzzzzz_xzzz_1[j];

                    tg_yzzzzzzz_yyyy_0[j] = pb_y * tg_zzzzzzz_yyyy_0[j] + wp_y[j] * tg_zzzzzzz_yyyy_1[j] + 2.0 * fl1_fxn * tg_zzzzzzz_yyy_1[j];

                    tg_yzzzzzzz_yyyz_0[j] = pb_y * tg_zzzzzzz_yyyz_0[j] + wp_y[j] * tg_zzzzzzz_yyyz_1[j] + 1.5 * fl1_fxn * tg_zzzzzzz_yyz_1[j];

                    tg_yzzzzzzz_yyzz_0[j] = pb_y * tg_zzzzzzz_yyzz_0[j] + wp_y[j] * tg_zzzzzzz_yyzz_1[j] + fl1_fxn * tg_zzzzzzz_yzz_1[j];

                    tg_yzzzzzzz_yzzz_0[j] = pb_y * tg_zzzzzzz_yzzz_0[j] + wp_y[j] * tg_zzzzzzz_yzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_zzz_1[j];

                    tg_yzzzzzzz_zzzz_0[j] = pb_y * tg_zzzzzzz_zzzz_0[j] + wp_y[j] * tg_zzzzzzz_zzzz_1[j];

                    tg_zzzzzzzz_xxxx_0[j] = pb_z * tg_zzzzzzz_xxxx_0[j] + wp_z[j] * tg_zzzzzzz_xxxx_1[j] + 3.5 * fl1_fx * tg_zzzzzz_xxxx_0[j] - 3.5 * fl1_fx * fl1_fza * tg_zzzzzz_xxxx_1[j];

                    tg_zzzzzzzz_xxxy_0[j] = pb_z * tg_zzzzzzz_xxxy_0[j] + wp_z[j] * tg_zzzzzzz_xxxy_1[j] + 3.5 * fl1_fx * tg_zzzzzz_xxxy_0[j] - 3.5 * fl1_fx * fl1_fza * tg_zzzzzz_xxxy_1[j];

                    tg_zzzzzzzz_xxxz_0[j] = pb_z * tg_zzzzzzz_xxxz_0[j] + wp_z[j] * tg_zzzzzzz_xxxz_1[j] + 3.5 * fl1_fx * tg_zzzzzz_xxxz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_zzzzzz_xxxz_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_xxx_1[j];

                    tg_zzzzzzzz_xxyy_0[j] = pb_z * tg_zzzzzzz_xxyy_0[j] + wp_z[j] * tg_zzzzzzz_xxyy_1[j] + 3.5 * fl1_fx * tg_zzzzzz_xxyy_0[j] - 3.5 * fl1_fx * fl1_fza * tg_zzzzzz_xxyy_1[j];

                    tg_zzzzzzzz_xxyz_0[j] = pb_z * tg_zzzzzzz_xxyz_0[j] + wp_z[j] * tg_zzzzzzz_xxyz_1[j] + 3.5 * fl1_fx * tg_zzzzzz_xxyz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_zzzzzz_xxyz_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_xxy_1[j];

                    tg_zzzzzzzz_xxzz_0[j] = pb_z * tg_zzzzzzz_xxzz_0[j] + wp_z[j] * tg_zzzzzzz_xxzz_1[j] + 3.5 * fl1_fx * tg_zzzzzz_xxzz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_zzzzzz_xxzz_1[j] + fl1_fxn * tg_zzzzzzz_xxz_1[j];

                    tg_zzzzzzzz_xyyy_0[j] = pb_z * tg_zzzzzzz_xyyy_0[j] + wp_z[j] * tg_zzzzzzz_xyyy_1[j] + 3.5 * fl1_fx * tg_zzzzzz_xyyy_0[j] - 3.5 * fl1_fx * fl1_fza * tg_zzzzzz_xyyy_1[j];

                    tg_zzzzzzzz_xyyz_0[j] = pb_z * tg_zzzzzzz_xyyz_0[j] + wp_z[j] * tg_zzzzzzz_xyyz_1[j] + 3.5 * fl1_fx * tg_zzzzzz_xyyz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_zzzzzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_xyy_1[j];

                    tg_zzzzzzzz_xyzz_0[j] = pb_z * tg_zzzzzzz_xyzz_0[j] + wp_z[j] * tg_zzzzzzz_xyzz_1[j] + 3.5 * fl1_fx * tg_zzzzzz_xyzz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_zzzzzz_xyzz_1[j] + fl1_fxn * tg_zzzzzzz_xyz_1[j];

                    tg_zzzzzzzz_xzzz_0[j] = pb_z * tg_zzzzzzz_xzzz_0[j] + wp_z[j] * tg_zzzzzzz_xzzz_1[j] + 3.5 * fl1_fx * tg_zzzzzz_xzzz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_zzzzzz_xzzz_1[j] + 1.5 * fl1_fxn * tg_zzzzzzz_xzz_1[j];

                    tg_zzzzzzzz_yyyy_0[j] = pb_z * tg_zzzzzzz_yyyy_0[j] + wp_z[j] * tg_zzzzzzz_yyyy_1[j] + 3.5 * fl1_fx * tg_zzzzzz_yyyy_0[j] - 3.5 * fl1_fx * fl1_fza * tg_zzzzzz_yyyy_1[j];

                    tg_zzzzzzzz_yyyz_0[j] = pb_z * tg_zzzzzzz_yyyz_0[j] + wp_z[j] * tg_zzzzzzz_yyyz_1[j] + 3.5 * fl1_fx * tg_zzzzzz_yyyz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_zzzzzz_yyyz_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_yyy_1[j];

                    tg_zzzzzzzz_yyzz_0[j] = pb_z * tg_zzzzzzz_yyzz_0[j] + wp_z[j] * tg_zzzzzzz_yyzz_1[j] + 3.5 * fl1_fx * tg_zzzzzz_yyzz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_zzzzzz_yyzz_1[j] + fl1_fxn * tg_zzzzzzz_yyz_1[j];

                    tg_zzzzzzzz_yzzz_0[j] = pb_z * tg_zzzzzzz_yzzz_0[j] + wp_z[j] * tg_zzzzzzz_yzzz_1[j] + 3.5 * fl1_fx * tg_zzzzzz_yzzz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_zzzzzz_yzzz_1[j] + 1.5 * fl1_fxn * tg_zzzzzzz_yzz_1[j];

                    tg_zzzzzzzz_zzzz_0[j] = pb_z * tg_zzzzzzz_zzzz_0[j] + wp_z[j] * tg_zzzzzzz_zzzz_1[j] + 3.5 * fl1_fx * tg_zzzzzz_zzzz_0[j] - 3.5 * fl1_fx * fl1_fza * tg_zzzzzz_zzzz_1[j] + 2.0 * fl1_fxn * tg_zzzzzzz_zzz_1[j];
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

