//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectronRepulsionRecFuncForIG.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSISG(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSISG_0_47(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSISG_47_94(primBuffer,
                                                       recursionMap,
                                                       osFactors,
                                                       wpDistances, 
                                                       braGtoPairsBlock,
                                                       ketGtoPairsBlock,
                                                       nKetPrimPairs,
                                                       iContrPair); 

        erirecfunc::compElectronRepulsionForSISG_94_141(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 

        erirecfunc::compElectronRepulsionForSISG_141_188(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSISG_188_235(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSISG_235_282(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSISG_282_328(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSISG_328_374(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSISG_374_420(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 
    }

    void
    compElectronRepulsionForSISG_0_47(      CMemBlock2D<double>& primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& wpDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetPrimPairs,
                                      const int32_t              iContrPair)
    {
        // Batch of Integrals (0,47)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {6, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_6_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_6_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_5_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_5_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_4_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_5_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_xxxxx_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx); 

                auto tg_xxxxx_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 1); 

                auto tg_xxxxx_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 2); 

                auto tg_xxxxx_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 3); 

                auto tg_xxxxx_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 4); 

                auto tg_xxxxx_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 5); 

                auto tg_xxxxx_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 6); 

                auto tg_xxxxx_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 7); 

                auto tg_xxxxx_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 8); 

                auto tg_xxxxx_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 9); 

                auto tg_xxxxx_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 10); 

                auto tg_xxxxx_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 11); 

                auto tg_xxxxx_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 12); 

                auto tg_xxxxx_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 13); 

                auto tg_xxxxx_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 14); 

                auto tg_xxxxy_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 15); 

                auto tg_xxxxy_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 16); 

                auto tg_xxxxy_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 17); 

                auto tg_xxxxy_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 18); 

                auto tg_xxxxy_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 19); 

                auto tg_xxxxy_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 20); 

                auto tg_xxxxy_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 21); 

                auto tg_xxxxy_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 22); 

                auto tg_xxxxy_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 23); 

                auto tg_xxxxy_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 24); 

                auto tg_xxxxy_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 25); 

                auto tg_xxxxy_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 26); 

                auto tg_xxxxy_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 27); 

                auto tg_xxxxy_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 28); 

                auto tg_xxxxy_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 29); 

                auto tg_xxxxz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 30); 

                auto tg_xxxxz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 31); 

                auto tg_xxxxz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 32); 

                auto tg_xxxxz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 33); 

                auto tg_xxxxz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 34); 

                auto tg_xxxxz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 35); 

                auto tg_xxxxz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 36); 

                auto tg_xxxxz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 37); 

                auto tg_xxxxz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 38); 

                auto tg_xxxxz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 39); 

                auto tg_xxxxz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 40); 

                auto tg_xxxxz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 41); 

                auto tg_xxxxz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 42); 

                auto tg_xxxxz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 43); 

                auto tg_xxxxz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 44); 

                auto tg_xxxyy_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 45); 

                auto tg_xxxyy_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 46); 

                auto tg_xxxxx_xxxx_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx); 

                auto tg_xxxxx_xxxy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 1); 

                auto tg_xxxxx_xxxz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 2); 

                auto tg_xxxxx_xxyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 3); 

                auto tg_xxxxx_xxyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 4); 

                auto tg_xxxxx_xxzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 5); 

                auto tg_xxxxx_xyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 6); 

                auto tg_xxxxx_xyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 7); 

                auto tg_xxxxx_xyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 8); 

                auto tg_xxxxx_xzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 9); 

                auto tg_xxxxx_yyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 10); 

                auto tg_xxxxx_yyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 11); 

                auto tg_xxxxx_yyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 12); 

                auto tg_xxxxx_yzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 13); 

                auto tg_xxxxx_zzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 14); 

                auto tg_xxxxy_xxxx_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 15); 

                auto tg_xxxxy_xxxy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 16); 

                auto tg_xxxxy_xxxz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 17); 

                auto tg_xxxxy_xxyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 18); 

                auto tg_xxxxy_xxyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 19); 

                auto tg_xxxxy_xxzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 20); 

                auto tg_xxxxy_xyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 21); 

                auto tg_xxxxy_xyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 22); 

                auto tg_xxxxy_xyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 23); 

                auto tg_xxxxy_xzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 24); 

                auto tg_xxxxy_yyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 25); 

                auto tg_xxxxy_yyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 26); 

                auto tg_xxxxy_yyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 27); 

                auto tg_xxxxy_yzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 28); 

                auto tg_xxxxy_zzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 29); 

                auto tg_xxxxz_xxxx_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 30); 

                auto tg_xxxxz_xxxy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 31); 

                auto tg_xxxxz_xxxz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 32); 

                auto tg_xxxxz_xxyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 33); 

                auto tg_xxxxz_xxyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 34); 

                auto tg_xxxxz_xxzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 35); 

                auto tg_xxxxz_xyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 36); 

                auto tg_xxxxz_xyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 37); 

                auto tg_xxxxz_xyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 38); 

                auto tg_xxxxz_xzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 39); 

                auto tg_xxxxz_yyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 40); 

                auto tg_xxxxz_yyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 41); 

                auto tg_xxxxz_yyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 42); 

                auto tg_xxxxz_yzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 43); 

                auto tg_xxxxz_zzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 44); 

                auto tg_xxxyy_xxxx_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 45); 

                auto tg_xxxyy_xxxy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 46); 

                auto tg_xxxx_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx); 

                auto tg_xxxx_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 1); 

                auto tg_xxxx_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 2); 

                auto tg_xxxx_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 3); 

                auto tg_xxxx_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 4); 

                auto tg_xxxx_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 5); 

                auto tg_xxxx_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 6); 

                auto tg_xxxx_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 7); 

                auto tg_xxxx_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 8); 

                auto tg_xxxx_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 9); 

                auto tg_xxxx_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 10); 

                auto tg_xxxx_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 11); 

                auto tg_xxxx_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 12); 

                auto tg_xxxx_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 13); 

                auto tg_xxxx_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 14); 

                auto tg_xxxy_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 15); 

                auto tg_xxxy_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 16); 

                auto tg_xxxy_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 17); 

                auto tg_xxxy_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 18); 

                auto tg_xxxy_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 19); 

                auto tg_xxxy_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 20); 

                auto tg_xxxy_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 21); 

                auto tg_xxxy_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 22); 

                auto tg_xxxy_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 23); 

                auto tg_xxxy_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 24); 

                auto tg_xxxy_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 25); 

                auto tg_xxxy_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 26); 

                auto tg_xxxy_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 27); 

                auto tg_xxxy_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 28); 

                auto tg_xxxy_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 29); 

                auto tg_xxxz_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 30); 

                auto tg_xxxz_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 31); 

                auto tg_xxxz_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 32); 

                auto tg_xxxz_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 33); 

                auto tg_xxxz_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 34); 

                auto tg_xxxz_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 35); 

                auto tg_xxxz_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 36); 

                auto tg_xxxz_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 37); 

                auto tg_xxxz_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 38); 

                auto tg_xxxz_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 39); 

                auto tg_xxxz_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 40); 

                auto tg_xxxz_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 41); 

                auto tg_xxxz_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 42); 

                auto tg_xxxz_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 43); 

                auto tg_xxxz_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 44); 

                auto tg_xxyy_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 45); 

                auto tg_xxyy_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 46); 

                auto tg_xxxx_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx); 

                auto tg_xxxx_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 1); 

                auto tg_xxxx_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 2); 

                auto tg_xxxx_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 3); 

                auto tg_xxxx_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 4); 

                auto tg_xxxx_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 5); 

                auto tg_xxxx_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 6); 

                auto tg_xxxx_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 7); 

                auto tg_xxxx_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 8); 

                auto tg_xxxx_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 9); 

                auto tg_xxxx_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 10); 

                auto tg_xxxx_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 11); 

                auto tg_xxxx_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 12); 

                auto tg_xxxx_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 13); 

                auto tg_xxxx_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 14); 

                auto tg_xxxy_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 15); 

                auto tg_xxxy_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 16); 

                auto tg_xxxy_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 17); 

                auto tg_xxxy_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 18); 

                auto tg_xxxy_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 19); 

                auto tg_xxxy_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 20); 

                auto tg_xxxy_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 21); 

                auto tg_xxxy_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 22); 

                auto tg_xxxy_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 23); 

                auto tg_xxxy_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 24); 

                auto tg_xxxy_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 25); 

                auto tg_xxxy_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 26); 

                auto tg_xxxy_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 27); 

                auto tg_xxxy_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 28); 

                auto tg_xxxy_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 29); 

                auto tg_xxxz_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 30); 

                auto tg_xxxz_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 31); 

                auto tg_xxxz_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 32); 

                auto tg_xxxz_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 33); 

                auto tg_xxxz_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 34); 

                auto tg_xxxz_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 35); 

                auto tg_xxxz_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 36); 

                auto tg_xxxz_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 37); 

                auto tg_xxxz_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 38); 

                auto tg_xxxz_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 39); 

                auto tg_xxxz_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 40); 

                auto tg_xxxz_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 41); 

                auto tg_xxxz_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 42); 

                auto tg_xxxz_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 43); 

                auto tg_xxxz_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 44); 

                auto tg_xxyy_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 45); 

                auto tg_xxyy_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 46); 

                auto tg_xxxxx_xxx_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx); 

                auto tg_xxxxx_xxy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 1); 

                auto tg_xxxxx_xxz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 2); 

                auto tg_xxxxx_xyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 3); 

                auto tg_xxxxx_xyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 4); 

                auto tg_xxxxx_xzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 5); 

                auto tg_xxxxx_yyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 6); 

                auto tg_xxxxx_yyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 7); 

                auto tg_xxxxx_yzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 8); 

                auto tg_xxxxx_zzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 9); 

                auto tg_xxxxy_xxx_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 10); 

                auto tg_xxxxy_xxy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 11); 

                auto tg_xxxxy_xxz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 12); 

                auto tg_xxxxy_xyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 13); 

                auto tg_xxxxy_xyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 14); 

                auto tg_xxxxy_xzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 15); 

                auto tg_xxxxy_yyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 16); 

                auto tg_xxxxy_yyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 17); 

                auto tg_xxxxy_yzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 18); 

                auto tg_xxxxy_zzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 19); 

                auto tg_xxxxz_xxx_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 20); 

                auto tg_xxxxz_xxy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 21); 

                auto tg_xxxxz_xxz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 22); 

                auto tg_xxxxz_xyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 23); 

                auto tg_xxxxz_xyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 24); 

                auto tg_xxxxz_xzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 25); 

                auto tg_xxxxz_yyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 26); 

                auto tg_xxxxz_yyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 27); 

                auto tg_xxxxz_yzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 28); 

                auto tg_xxxxz_zzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 29); 

                auto tg_xxxyy_xxx_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 30); 

                auto tg_xxxyy_xxy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 31); 

                // set up pointers to integrals

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

                // Batch of Integrals (0,47)

                #pragma omp simd aligned(fxn, fza, tg_xxxx_xxxx_0, tg_xxxx_xxxx_1, tg_xxxx_xxxy_0, \
                                         tg_xxxx_xxxy_1, tg_xxxx_xxxz_0, tg_xxxx_xxxz_1, tg_xxxx_xxyy_0, tg_xxxx_xxyy_1, \
                                         tg_xxxx_xxyz_0, tg_xxxx_xxyz_1, tg_xxxx_xxzz_0, tg_xxxx_xxzz_1, tg_xxxx_xyyy_0, \
                                         tg_xxxx_xyyy_1, tg_xxxx_xyyz_0, tg_xxxx_xyyz_1, tg_xxxx_xyzz_0, tg_xxxx_xyzz_1, \
                                         tg_xxxx_xzzz_0, tg_xxxx_xzzz_1, tg_xxxx_yyyy_0, tg_xxxx_yyyy_1, tg_xxxx_yyyz_0, \
                                         tg_xxxx_yyyz_1, tg_xxxx_yyzz_0, tg_xxxx_yyzz_1, tg_xxxx_yzzz_0, tg_xxxx_yzzz_1, \
                                         tg_xxxx_zzzz_0, tg_xxxx_zzzz_1, tg_xxxxx_xxx_1, tg_xxxxx_xxxx_0, tg_xxxxx_xxxx_1, \
                                         tg_xxxxx_xxxy_0, tg_xxxxx_xxxy_1, tg_xxxxx_xxxz_0, tg_xxxxx_xxxz_1, tg_xxxxx_xxy_1, \
                                         tg_xxxxx_xxyy_0, tg_xxxxx_xxyy_1, tg_xxxxx_xxyz_0, tg_xxxxx_xxyz_1, tg_xxxxx_xxz_1, \
                                         tg_xxxxx_xxzz_0, tg_xxxxx_xxzz_1, tg_xxxxx_xyy_1, tg_xxxxx_xyyy_0, tg_xxxxx_xyyy_1, \
                                         tg_xxxxx_xyyz_0, tg_xxxxx_xyyz_1, tg_xxxxx_xyz_1, tg_xxxxx_xyzz_0, tg_xxxxx_xyzz_1, \
                                         tg_xxxxx_xzz_1, tg_xxxxx_xzzz_0, tg_xxxxx_xzzz_1, tg_xxxxx_yyy_1, tg_xxxxx_yyyy_0, \
                                         tg_xxxxx_yyyy_1, tg_xxxxx_yyyz_0, tg_xxxxx_yyyz_1, tg_xxxxx_yyz_1, tg_xxxxx_yyzz_0, \
                                         tg_xxxxx_yyzz_1, tg_xxxxx_yzz_1, tg_xxxxx_yzzz_0, tg_xxxxx_yzzz_1, tg_xxxxx_zzz_1, \
                                         tg_xxxxx_zzzz_0, tg_xxxxx_zzzz_1, tg_xxxxxx_xxxx_0, tg_xxxxxx_xxxy_0, \
                                         tg_xxxxxx_xxxz_0, tg_xxxxxx_xxyy_0, tg_xxxxxx_xxyz_0, tg_xxxxxx_xxzz_0, \
                                         tg_xxxxxx_xyyy_0, tg_xxxxxx_xyyz_0, tg_xxxxxx_xyzz_0, tg_xxxxxx_xzzz_0, \
                                         tg_xxxxxx_yyyy_0, tg_xxxxxx_yyyz_0, tg_xxxxxx_yyzz_0, tg_xxxxxx_yzzz_0, \
                                         tg_xxxxxx_zzzz_0, tg_xxxxxy_xxxx_0, tg_xxxxxy_xxxy_0, tg_xxxxxy_xxxz_0, \
                                         tg_xxxxxy_xxyy_0, tg_xxxxxy_xxyz_0, tg_xxxxxy_xxzz_0, tg_xxxxxy_xyyy_0, \
                                         tg_xxxxxy_xyyz_0, tg_xxxxxy_xyzz_0, tg_xxxxxy_xzzz_0, tg_xxxxxy_yyyy_0, \
                                         tg_xxxxxy_yyyz_0, tg_xxxxxy_yyzz_0, tg_xxxxxy_yzzz_0, tg_xxxxxy_zzzz_0, \
                                         tg_xxxxxz_xxxx_0, tg_xxxxxz_xxxy_0, tg_xxxxxz_xxxz_0, tg_xxxxxz_xxyy_0, \
                                         tg_xxxxxz_xxyz_0, tg_xxxxxz_xxzz_0, tg_xxxxxz_xyyy_0, tg_xxxxxz_xyyz_0, \
                                         tg_xxxxxz_xyzz_0, tg_xxxxxz_xzzz_0, tg_xxxxxz_yyyy_0, tg_xxxxxz_yyyz_0, \
                                         tg_xxxxxz_yyzz_0, tg_xxxxxz_yzzz_0, tg_xxxxxz_zzzz_0, tg_xxxxy_xxx_1, \
                                         tg_xxxxy_xxxx_0, tg_xxxxy_xxxx_1, tg_xxxxy_xxxy_0, tg_xxxxy_xxxy_1, tg_xxxxy_xxxz_0, \
                                         tg_xxxxy_xxxz_1, tg_xxxxy_xxy_1, tg_xxxxy_xxyy_0, tg_xxxxy_xxyy_1, tg_xxxxy_xxyz_0, \
                                         tg_xxxxy_xxyz_1, tg_xxxxy_xxz_1, tg_xxxxy_xxzz_0, tg_xxxxy_xxzz_1, tg_xxxxy_xyy_1, \
                                         tg_xxxxy_xyyy_0, tg_xxxxy_xyyy_1, tg_xxxxy_xyyz_0, tg_xxxxy_xyyz_1, tg_xxxxy_xyz_1, \
                                         tg_xxxxy_xyzz_0, tg_xxxxy_xyzz_1, tg_xxxxy_xzz_1, tg_xxxxy_xzzz_0, tg_xxxxy_xzzz_1, \
                                         tg_xxxxy_yyy_1, tg_xxxxy_yyyy_0, tg_xxxxy_yyyy_1, tg_xxxxy_yyyz_0, tg_xxxxy_yyyz_1, \
                                         tg_xxxxy_yyz_1, tg_xxxxy_yyzz_0, tg_xxxxy_yyzz_1, tg_xxxxy_yzz_1, tg_xxxxy_yzzz_0, \
                                         tg_xxxxy_yzzz_1, tg_xxxxy_zzz_1, tg_xxxxy_zzzz_0, tg_xxxxy_zzzz_1, tg_xxxxyy_xxxx_0, \
                                         tg_xxxxyy_xxxy_0, tg_xxxxz_xxx_1, tg_xxxxz_xxxx_0, tg_xxxxz_xxxx_1, tg_xxxxz_xxxy_0, \
                                         tg_xxxxz_xxxy_1, tg_xxxxz_xxxz_0, tg_xxxxz_xxxz_1, tg_xxxxz_xxy_1, tg_xxxxz_xxyy_0, \
                                         tg_xxxxz_xxyy_1, tg_xxxxz_xxyz_0, tg_xxxxz_xxyz_1, tg_xxxxz_xxz_1, tg_xxxxz_xxzz_0, \
                                         tg_xxxxz_xxzz_1, tg_xxxxz_xyy_1, tg_xxxxz_xyyy_0, tg_xxxxz_xyyy_1, tg_xxxxz_xyyz_0, \
                                         tg_xxxxz_xyyz_1, tg_xxxxz_xyz_1, tg_xxxxz_xyzz_0, tg_xxxxz_xyzz_1, tg_xxxxz_xzz_1, \
                                         tg_xxxxz_xzzz_0, tg_xxxxz_xzzz_1, tg_xxxxz_yyy_1, tg_xxxxz_yyyy_0, tg_xxxxz_yyyy_1, \
                                         tg_xxxxz_yyyz_0, tg_xxxxz_yyyz_1, tg_xxxxz_yyz_1, tg_xxxxz_yyzz_0, tg_xxxxz_yyzz_1, \
                                         tg_xxxxz_yzz_1, tg_xxxxz_yzzz_0, tg_xxxxz_yzzz_1, tg_xxxxz_zzz_1, tg_xxxxz_zzzz_0, \
                                         tg_xxxxz_zzzz_1, tg_xxxy_xxxx_0, tg_xxxy_xxxx_1, tg_xxxy_xxxy_0, tg_xxxy_xxxy_1, \
                                         tg_xxxy_xxxz_0, tg_xxxy_xxxz_1, tg_xxxy_xxyy_0, tg_xxxy_xxyy_1, tg_xxxy_xxyz_0, \
                                         tg_xxxy_xxyz_1, tg_xxxy_xxzz_0, tg_xxxy_xxzz_1, tg_xxxy_xyyy_0, tg_xxxy_xyyy_1, \
                                         tg_xxxy_xyyz_0, tg_xxxy_xyyz_1, tg_xxxy_xyzz_0, tg_xxxy_xyzz_1, tg_xxxy_xzzz_0, \
                                         tg_xxxy_xzzz_1, tg_xxxy_yyyy_0, tg_xxxy_yyyy_1, tg_xxxy_yyyz_0, tg_xxxy_yyyz_1, \
                                         tg_xxxy_yyzz_0, tg_xxxy_yyzz_1, tg_xxxy_yzzz_0, tg_xxxy_yzzz_1, tg_xxxy_zzzz_0, \
                                         tg_xxxy_zzzz_1, tg_xxxyy_xxx_1, tg_xxxyy_xxxx_0, tg_xxxyy_xxxx_1, tg_xxxyy_xxxy_0, \
                                         tg_xxxyy_xxxy_1, tg_xxxyy_xxy_1, tg_xxxz_xxxx_0, tg_xxxz_xxxx_1, tg_xxxz_xxxy_0, \
                                         tg_xxxz_xxxy_1, tg_xxxz_xxxz_0, tg_xxxz_xxxz_1, tg_xxxz_xxyy_0, tg_xxxz_xxyy_1, \
                                         tg_xxxz_xxyz_0, tg_xxxz_xxyz_1, tg_xxxz_xxzz_0, tg_xxxz_xxzz_1, tg_xxxz_xyyy_0, \
                                         tg_xxxz_xyyy_1, tg_xxxz_xyyz_0, tg_xxxz_xyyz_1, tg_xxxz_xyzz_0, tg_xxxz_xyzz_1, \
                                         tg_xxxz_xzzz_0, tg_xxxz_xzzz_1, tg_xxxz_yyyy_0, tg_xxxz_yyyy_1, tg_xxxz_yyyz_0, \
                                         tg_xxxz_yyyz_1, tg_xxxz_yyzz_0, tg_xxxz_yyzz_1, tg_xxxz_yzzz_0, tg_xxxz_yzzz_1, \
                                         tg_xxxz_zzzz_0, tg_xxxz_zzzz_1, tg_xxyy_xxxx_0, tg_xxyy_xxxx_1, tg_xxyy_xxxy_0, \
                                         tg_xxyy_xxxy_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxxxxx_xxxx_0[j] = pb_x * tg_xxxxx_xxxx_0[j] + wp_x[j] * tg_xxxxx_xxxx_1[j] + 2.5 * fl1_fx * tg_xxxx_xxxx_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxx_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxxxx_xxx_1[j];

                    tg_xxxxxx_xxxy_0[j] = pb_x * tg_xxxxx_xxxy_0[j] + wp_x[j] * tg_xxxxx_xxxy_1[j] + 2.5 * fl1_fx * tg_xxxx_xxxy_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxx_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxxxx_xxy_1[j];

                    tg_xxxxxx_xxxz_0[j] = pb_x * tg_xxxxx_xxxz_0[j] + wp_x[j] * tg_xxxxx_xxxz_1[j] + 2.5 * fl1_fx * tg_xxxx_xxxz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxx_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxxxx_xxz_1[j];

                    tg_xxxxxx_xxyy_0[j] = pb_x * tg_xxxxx_xxyy_0[j] + wp_x[j] * tg_xxxxx_xxyy_1[j] + 2.5 * fl1_fx * tg_xxxx_xxyy_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxx_xxyy_1[j] + fl1_fxn * tg_xxxxx_xyy_1[j];

                    tg_xxxxxx_xxyz_0[j] = pb_x * tg_xxxxx_xxyz_0[j] + wp_x[j] * tg_xxxxx_xxyz_1[j] + 2.5 * fl1_fx * tg_xxxx_xxyz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxx_xxyz_1[j] + fl1_fxn * tg_xxxxx_xyz_1[j];

                    tg_xxxxxx_xxzz_0[j] = pb_x * tg_xxxxx_xxzz_0[j] + wp_x[j] * tg_xxxxx_xxzz_1[j] + 2.5 * fl1_fx * tg_xxxx_xxzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxx_xxzz_1[j] + fl1_fxn * tg_xxxxx_xzz_1[j];

                    tg_xxxxxx_xyyy_0[j] = pb_x * tg_xxxxx_xyyy_0[j] + wp_x[j] * tg_xxxxx_xyyy_1[j] + 2.5 * fl1_fx * tg_xxxx_xyyy_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxx_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxxxx_yyy_1[j];

                    tg_xxxxxx_xyyz_0[j] = pb_x * tg_xxxxx_xyyz_0[j] + wp_x[j] * tg_xxxxx_xyyz_1[j] + 2.5 * fl1_fx * tg_xxxx_xyyz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxx_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxxxx_yyz_1[j];

                    tg_xxxxxx_xyzz_0[j] = pb_x * tg_xxxxx_xyzz_0[j] + wp_x[j] * tg_xxxxx_xyzz_1[j] + 2.5 * fl1_fx * tg_xxxx_xyzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxx_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxxxx_yzz_1[j];

                    tg_xxxxxx_xzzz_0[j] = pb_x * tg_xxxxx_xzzz_0[j] + wp_x[j] * tg_xxxxx_xzzz_1[j] + 2.5 * fl1_fx * tg_xxxx_xzzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxx_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxxxx_zzz_1[j];

                    tg_xxxxxx_yyyy_0[j] = pb_x * tg_xxxxx_yyyy_0[j] + wp_x[j] * tg_xxxxx_yyyy_1[j] + 2.5 * fl1_fx * tg_xxxx_yyyy_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxx_yyyy_1[j];

                    tg_xxxxxx_yyyz_0[j] = pb_x * tg_xxxxx_yyyz_0[j] + wp_x[j] * tg_xxxxx_yyyz_1[j] + 2.5 * fl1_fx * tg_xxxx_yyyz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxx_yyyz_1[j];

                    tg_xxxxxx_yyzz_0[j] = pb_x * tg_xxxxx_yyzz_0[j] + wp_x[j] * tg_xxxxx_yyzz_1[j] + 2.5 * fl1_fx * tg_xxxx_yyzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxx_yyzz_1[j];

                    tg_xxxxxx_yzzz_0[j] = pb_x * tg_xxxxx_yzzz_0[j] + wp_x[j] * tg_xxxxx_yzzz_1[j] + 2.5 * fl1_fx * tg_xxxx_yzzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxx_yzzz_1[j];

                    tg_xxxxxx_zzzz_0[j] = pb_x * tg_xxxxx_zzzz_0[j] + wp_x[j] * tg_xxxxx_zzzz_1[j] + 2.5 * fl1_fx * tg_xxxx_zzzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxx_zzzz_1[j];

                    tg_xxxxxy_xxxx_0[j] = pb_x * tg_xxxxy_xxxx_0[j] + wp_x[j] * tg_xxxxy_xxxx_1[j] + 2.0 * fl1_fx * tg_xxxy_xxxx_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxy_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxxxy_xxx_1[j];

                    tg_xxxxxy_xxxy_0[j] = pb_x * tg_xxxxy_xxxy_0[j] + wp_x[j] * tg_xxxxy_xxxy_1[j] + 2.0 * fl1_fx * tg_xxxy_xxxy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxy_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxxxy_xxy_1[j];

                    tg_xxxxxy_xxxz_0[j] = pb_x * tg_xxxxy_xxxz_0[j] + wp_x[j] * tg_xxxxy_xxxz_1[j] + 2.0 * fl1_fx * tg_xxxy_xxxz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxy_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxxxy_xxz_1[j];

                    tg_xxxxxy_xxyy_0[j] = pb_x * tg_xxxxy_xxyy_0[j] + wp_x[j] * tg_xxxxy_xxyy_1[j] + 2.0 * fl1_fx * tg_xxxy_xxyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxy_xxyy_1[j] + fl1_fxn * tg_xxxxy_xyy_1[j];

                    tg_xxxxxy_xxyz_0[j] = pb_x * tg_xxxxy_xxyz_0[j] + wp_x[j] * tg_xxxxy_xxyz_1[j] + 2.0 * fl1_fx * tg_xxxy_xxyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxy_xxyz_1[j] + fl1_fxn * tg_xxxxy_xyz_1[j];

                    tg_xxxxxy_xxzz_0[j] = pb_x * tg_xxxxy_xxzz_0[j] + wp_x[j] * tg_xxxxy_xxzz_1[j] + 2.0 * fl1_fx * tg_xxxy_xxzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxy_xxzz_1[j] + fl1_fxn * tg_xxxxy_xzz_1[j];

                    tg_xxxxxy_xyyy_0[j] = pb_x * tg_xxxxy_xyyy_0[j] + wp_x[j] * tg_xxxxy_xyyy_1[j] + 2.0 * fl1_fx * tg_xxxy_xyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxy_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxxxy_yyy_1[j];

                    tg_xxxxxy_xyyz_0[j] = pb_x * tg_xxxxy_xyyz_0[j] + wp_x[j] * tg_xxxxy_xyyz_1[j] + 2.0 * fl1_fx * tg_xxxy_xyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxy_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxxxy_yyz_1[j];

                    tg_xxxxxy_xyzz_0[j] = pb_x * tg_xxxxy_xyzz_0[j] + wp_x[j] * tg_xxxxy_xyzz_1[j] + 2.0 * fl1_fx * tg_xxxy_xyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxy_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxxxy_yzz_1[j];

                    tg_xxxxxy_xzzz_0[j] = pb_x * tg_xxxxy_xzzz_0[j] + wp_x[j] * tg_xxxxy_xzzz_1[j] + 2.0 * fl1_fx * tg_xxxy_xzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxy_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxxxy_zzz_1[j];

                    tg_xxxxxy_yyyy_0[j] = pb_x * tg_xxxxy_yyyy_0[j] + wp_x[j] * tg_xxxxy_yyyy_1[j] + 2.0 * fl1_fx * tg_xxxy_yyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxy_yyyy_1[j];

                    tg_xxxxxy_yyyz_0[j] = pb_x * tg_xxxxy_yyyz_0[j] + wp_x[j] * tg_xxxxy_yyyz_1[j] + 2.0 * fl1_fx * tg_xxxy_yyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxy_yyyz_1[j];

                    tg_xxxxxy_yyzz_0[j] = pb_x * tg_xxxxy_yyzz_0[j] + wp_x[j] * tg_xxxxy_yyzz_1[j] + 2.0 * fl1_fx * tg_xxxy_yyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxy_yyzz_1[j];

                    tg_xxxxxy_yzzz_0[j] = pb_x * tg_xxxxy_yzzz_0[j] + wp_x[j] * tg_xxxxy_yzzz_1[j] + 2.0 * fl1_fx * tg_xxxy_yzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxy_yzzz_1[j];

                    tg_xxxxxy_zzzz_0[j] = pb_x * tg_xxxxy_zzzz_0[j] + wp_x[j] * tg_xxxxy_zzzz_1[j] + 2.0 * fl1_fx * tg_xxxy_zzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxy_zzzz_1[j];

                    tg_xxxxxz_xxxx_0[j] = pb_x * tg_xxxxz_xxxx_0[j] + wp_x[j] * tg_xxxxz_xxxx_1[j] + 2.0 * fl1_fx * tg_xxxz_xxxx_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxxxz_xxx_1[j];

                    tg_xxxxxz_xxxy_0[j] = pb_x * tg_xxxxz_xxxy_0[j] + wp_x[j] * tg_xxxxz_xxxy_1[j] + 2.0 * fl1_fx * tg_xxxz_xxxy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxxxz_xxy_1[j];

                    tg_xxxxxz_xxxz_0[j] = pb_x * tg_xxxxz_xxxz_0[j] + wp_x[j] * tg_xxxxz_xxxz_1[j] + 2.0 * fl1_fx * tg_xxxz_xxxz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxxxz_xxz_1[j];

                    tg_xxxxxz_xxyy_0[j] = pb_x * tg_xxxxz_xxyy_0[j] + wp_x[j] * tg_xxxxz_xxyy_1[j] + 2.0 * fl1_fx * tg_xxxz_xxyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxz_xxyy_1[j] + fl1_fxn * tg_xxxxz_xyy_1[j];

                    tg_xxxxxz_xxyz_0[j] = pb_x * tg_xxxxz_xxyz_0[j] + wp_x[j] * tg_xxxxz_xxyz_1[j] + 2.0 * fl1_fx * tg_xxxz_xxyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxz_xxyz_1[j] + fl1_fxn * tg_xxxxz_xyz_1[j];

                    tg_xxxxxz_xxzz_0[j] = pb_x * tg_xxxxz_xxzz_0[j] + wp_x[j] * tg_xxxxz_xxzz_1[j] + 2.0 * fl1_fx * tg_xxxz_xxzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxz_xxzz_1[j] + fl1_fxn * tg_xxxxz_xzz_1[j];

                    tg_xxxxxz_xyyy_0[j] = pb_x * tg_xxxxz_xyyy_0[j] + wp_x[j] * tg_xxxxz_xyyy_1[j] + 2.0 * fl1_fx * tg_xxxz_xyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxxxz_yyy_1[j];

                    tg_xxxxxz_xyyz_0[j] = pb_x * tg_xxxxz_xyyz_0[j] + wp_x[j] * tg_xxxxz_xyyz_1[j] + 2.0 * fl1_fx * tg_xxxz_xyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxxxz_yyz_1[j];

                    tg_xxxxxz_xyzz_0[j] = pb_x * tg_xxxxz_xyzz_0[j] + wp_x[j] * tg_xxxxz_xyzz_1[j] + 2.0 * fl1_fx * tg_xxxz_xyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxxxz_yzz_1[j];

                    tg_xxxxxz_xzzz_0[j] = pb_x * tg_xxxxz_xzzz_0[j] + wp_x[j] * tg_xxxxz_xzzz_1[j] + 2.0 * fl1_fx * tg_xxxz_xzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxxxz_zzz_1[j];

                    tg_xxxxxz_yyyy_0[j] = pb_x * tg_xxxxz_yyyy_0[j] + wp_x[j] * tg_xxxxz_yyyy_1[j] + 2.0 * fl1_fx * tg_xxxz_yyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxz_yyyy_1[j];

                    tg_xxxxxz_yyyz_0[j] = pb_x * tg_xxxxz_yyyz_0[j] + wp_x[j] * tg_xxxxz_yyyz_1[j] + 2.0 * fl1_fx * tg_xxxz_yyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxz_yyyz_1[j];

                    tg_xxxxxz_yyzz_0[j] = pb_x * tg_xxxxz_yyzz_0[j] + wp_x[j] * tg_xxxxz_yyzz_1[j] + 2.0 * fl1_fx * tg_xxxz_yyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxz_yyzz_1[j];

                    tg_xxxxxz_yzzz_0[j] = pb_x * tg_xxxxz_yzzz_0[j] + wp_x[j] * tg_xxxxz_yzzz_1[j] + 2.0 * fl1_fx * tg_xxxz_yzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxz_yzzz_1[j];

                    tg_xxxxxz_zzzz_0[j] = pb_x * tg_xxxxz_zzzz_0[j] + wp_x[j] * tg_xxxxz_zzzz_1[j] + 2.0 * fl1_fx * tg_xxxz_zzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxz_zzzz_1[j];

                    tg_xxxxyy_xxxx_0[j] = pb_x * tg_xxxyy_xxxx_0[j] + wp_x[j] * tg_xxxyy_xxxx_1[j] + 1.5 * fl1_fx * tg_xxyy_xxxx_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyy_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxxyy_xxx_1[j];

                    tg_xxxxyy_xxxy_0[j] = pb_x * tg_xxxyy_xxxy_0[j] + wp_x[j] * tg_xxxyy_xxxy_1[j] + 1.5 * fl1_fx * tg_xxyy_xxxy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyy_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxxyy_xxy_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSISG_47_94(      CMemBlock2D<double>& primBuffer,
                                       const CRecursionMap&       recursionMap,
                                       const CMemBlock2D<double>& osFactors,
                                       const CMemBlock2D<double>& wpDistances,
                                       const CGtoPairsBlock&      braGtoPairsBlock,
                                       const CGtoPairsBlock&      ketGtoPairsBlock,
                                       const int32_t              nKetPrimPairs,
                                       const int32_t              iContrPair)
    {
        // Batch of Integrals (47,94)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {6, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_6_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_6_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_5_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_5_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_4_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_5_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_xxxyy_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 47); 

                auto tg_xxxyy_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 48); 

                auto tg_xxxyy_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 49); 

                auto tg_xxxyy_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 50); 

                auto tg_xxxyy_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 51); 

                auto tg_xxxyy_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 52); 

                auto tg_xxxyy_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 53); 

                auto tg_xxxyy_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 54); 

                auto tg_xxxyy_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 55); 

                auto tg_xxxyy_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 56); 

                auto tg_xxxyy_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 57); 

                auto tg_xxxyy_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 58); 

                auto tg_xxxyy_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 59); 

                auto tg_xxxyz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 60); 

                auto tg_xxxyz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 61); 

                auto tg_xxxyz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 62); 

                auto tg_xxxyz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 63); 

                auto tg_xxxyz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 64); 

                auto tg_xxxyz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 65); 

                auto tg_xxxyz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 66); 

                auto tg_xxxyz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 67); 

                auto tg_xxxyz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 68); 

                auto tg_xxxyz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 69); 

                auto tg_xxxyz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 70); 

                auto tg_xxxyz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 71); 

                auto tg_xxxyz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 72); 

                auto tg_xxxyz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 73); 

                auto tg_xxxyz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 74); 

                auto tg_xxxzz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 75); 

                auto tg_xxxzz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 76); 

                auto tg_xxxzz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 77); 

                auto tg_xxxzz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 78); 

                auto tg_xxxzz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 79); 

                auto tg_xxxzz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 80); 

                auto tg_xxxzz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 81); 

                auto tg_xxxzz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 82); 

                auto tg_xxxzz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 83); 

                auto tg_xxxzz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 84); 

                auto tg_xxxzz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 85); 

                auto tg_xxxzz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 86); 

                auto tg_xxxzz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 87); 

                auto tg_xxxzz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 88); 

                auto tg_xxxzz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 89); 

                auto tg_xxyyy_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 90); 

                auto tg_xxyyy_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 91); 

                auto tg_xxyyy_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 92); 

                auto tg_xxyyy_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 93); 

                auto tg_xxxyy_xxxz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 47); 

                auto tg_xxxyy_xxyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 48); 

                auto tg_xxxyy_xxyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 49); 

                auto tg_xxxyy_xxzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 50); 

                auto tg_xxxyy_xyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 51); 

                auto tg_xxxyy_xyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 52); 

                auto tg_xxxyy_xyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 53); 

                auto tg_xxxyy_xzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 54); 

                auto tg_xxxyy_yyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 55); 

                auto tg_xxxyy_yyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 56); 

                auto tg_xxxyy_yyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 57); 

                auto tg_xxxyy_yzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 58); 

                auto tg_xxxyy_zzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 59); 

                auto tg_xxxyz_xxxx_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 60); 

                auto tg_xxxyz_xxxy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 61); 

                auto tg_xxxyz_xxxz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 62); 

                auto tg_xxxyz_xxyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 63); 

                auto tg_xxxyz_xxyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 64); 

                auto tg_xxxyz_xxzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 65); 

                auto tg_xxxyz_xyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 66); 

                auto tg_xxxyz_xyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 67); 

                auto tg_xxxyz_xyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 68); 

                auto tg_xxxyz_xzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 69); 

                auto tg_xxxyz_yyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 70); 

                auto tg_xxxyz_yyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 71); 

                auto tg_xxxyz_yyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 72); 

                auto tg_xxxyz_yzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 73); 

                auto tg_xxxyz_zzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 74); 

                auto tg_xxxzz_xxxx_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 75); 

                auto tg_xxxzz_xxxy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 76); 

                auto tg_xxxzz_xxxz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 77); 

                auto tg_xxxzz_xxyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 78); 

                auto tg_xxxzz_xxyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 79); 

                auto tg_xxxzz_xxzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 80); 

                auto tg_xxxzz_xyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 81); 

                auto tg_xxxzz_xyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 82); 

                auto tg_xxxzz_xyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 83); 

                auto tg_xxxzz_xzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 84); 

                auto tg_xxxzz_yyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 85); 

                auto tg_xxxzz_yyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 86); 

                auto tg_xxxzz_yyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 87); 

                auto tg_xxxzz_yzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 88); 

                auto tg_xxxzz_zzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 89); 

                auto tg_xxyyy_xxxx_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 90); 

                auto tg_xxyyy_xxxy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 91); 

                auto tg_xxyyy_xxxz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 92); 

                auto tg_xxyyy_xxyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 93); 

                auto tg_xxyy_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 47); 

                auto tg_xxyy_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 48); 

                auto tg_xxyy_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 49); 

                auto tg_xxyy_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 50); 

                auto tg_xxyy_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 51); 

                auto tg_xxyy_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 52); 

                auto tg_xxyy_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 53); 

                auto tg_xxyy_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 54); 

                auto tg_xxyy_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 55); 

                auto tg_xxyy_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 56); 

                auto tg_xxyy_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 57); 

                auto tg_xxyy_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 58); 

                auto tg_xxyy_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 59); 

                auto tg_xxyz_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 60); 

                auto tg_xxyz_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 61); 

                auto tg_xxyz_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 62); 

                auto tg_xxyz_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 63); 

                auto tg_xxyz_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 64); 

                auto tg_xxyz_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 65); 

                auto tg_xxyz_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 66); 

                auto tg_xxyz_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 67); 

                auto tg_xxyz_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 68); 

                auto tg_xxyz_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 69); 

                auto tg_xxyz_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 70); 

                auto tg_xxyz_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 71); 

                auto tg_xxyz_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 72); 

                auto tg_xxyz_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 73); 

                auto tg_xxyz_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 74); 

                auto tg_xxzz_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 75); 

                auto tg_xxzz_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 76); 

                auto tg_xxzz_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 77); 

                auto tg_xxzz_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 78); 

                auto tg_xxzz_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 79); 

                auto tg_xxzz_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 80); 

                auto tg_xxzz_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 81); 

                auto tg_xxzz_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 82); 

                auto tg_xxzz_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 83); 

                auto tg_xxzz_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 84); 

                auto tg_xxzz_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 85); 

                auto tg_xxzz_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 86); 

                auto tg_xxzz_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 87); 

                auto tg_xxzz_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 88); 

                auto tg_xxzz_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 89); 

                auto tg_xyyy_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 90); 

                auto tg_xyyy_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 91); 

                auto tg_xyyy_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 92); 

                auto tg_xyyy_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 93); 

                auto tg_xxyy_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 47); 

                auto tg_xxyy_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 48); 

                auto tg_xxyy_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 49); 

                auto tg_xxyy_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 50); 

                auto tg_xxyy_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 51); 

                auto tg_xxyy_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 52); 

                auto tg_xxyy_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 53); 

                auto tg_xxyy_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 54); 

                auto tg_xxyy_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 55); 

                auto tg_xxyy_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 56); 

                auto tg_xxyy_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 57); 

                auto tg_xxyy_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 58); 

                auto tg_xxyy_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 59); 

                auto tg_xxyz_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 60); 

                auto tg_xxyz_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 61); 

                auto tg_xxyz_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 62); 

                auto tg_xxyz_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 63); 

                auto tg_xxyz_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 64); 

                auto tg_xxyz_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 65); 

                auto tg_xxyz_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 66); 

                auto tg_xxyz_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 67); 

                auto tg_xxyz_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 68); 

                auto tg_xxyz_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 69); 

                auto tg_xxyz_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 70); 

                auto tg_xxyz_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 71); 

                auto tg_xxyz_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 72); 

                auto tg_xxyz_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 73); 

                auto tg_xxyz_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 74); 

                auto tg_xxzz_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 75); 

                auto tg_xxzz_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 76); 

                auto tg_xxzz_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 77); 

                auto tg_xxzz_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 78); 

                auto tg_xxzz_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 79); 

                auto tg_xxzz_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 80); 

                auto tg_xxzz_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 81); 

                auto tg_xxzz_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 82); 

                auto tg_xxzz_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 83); 

                auto tg_xxzz_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 84); 

                auto tg_xxzz_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 85); 

                auto tg_xxzz_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 86); 

                auto tg_xxzz_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 87); 

                auto tg_xxzz_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 88); 

                auto tg_xxzz_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 89); 

                auto tg_xyyy_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 90); 

                auto tg_xyyy_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 91); 

                auto tg_xyyy_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 92); 

                auto tg_xyyy_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 93); 

                auto tg_xxxyy_xxz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 32); 

                auto tg_xxxyy_xyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 33); 

                auto tg_xxxyy_xyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 34); 

                auto tg_xxxyy_xzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 35); 

                auto tg_xxxyy_yyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 36); 

                auto tg_xxxyy_yyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 37); 

                auto tg_xxxyy_yzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 38); 

                auto tg_xxxyy_zzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 39); 

                auto tg_xxxyz_xxx_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 40); 

                auto tg_xxxyz_xxy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 41); 

                auto tg_xxxyz_xxz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 42); 

                auto tg_xxxyz_xyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 43); 

                auto tg_xxxyz_xyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 44); 

                auto tg_xxxyz_xzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 45); 

                auto tg_xxxyz_yyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 46); 

                auto tg_xxxyz_yyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 47); 

                auto tg_xxxyz_yzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 48); 

                auto tg_xxxyz_zzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 49); 

                auto tg_xxxzz_xxx_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 50); 

                auto tg_xxxzz_xxy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 51); 

                auto tg_xxxzz_xxz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 52); 

                auto tg_xxxzz_xyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 53); 

                auto tg_xxxzz_xyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 54); 

                auto tg_xxxzz_xzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 55); 

                auto tg_xxxzz_yyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 56); 

                auto tg_xxxzz_yyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 57); 

                auto tg_xxxzz_yzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 58); 

                auto tg_xxxzz_zzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 59); 

                auto tg_xxyyy_xxx_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 60); 

                auto tg_xxyyy_xxy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 61); 

                auto tg_xxyyy_xxz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 62); 

                auto tg_xxyyy_xyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 63); 

                // set up pointers to integrals

                auto tg_xxxxyy_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 47); 

                auto tg_xxxxyy_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 48); 

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

                // Batch of Integrals (47,94)

                #pragma omp simd aligned(fxn, fza, tg_xxxxyy_xxxz_0, tg_xxxxyy_xxyy_0, tg_xxxxyy_xxyz_0, \
                                         tg_xxxxyy_xxzz_0, tg_xxxxyy_xyyy_0, tg_xxxxyy_xyyz_0, tg_xxxxyy_xyzz_0, \
                                         tg_xxxxyy_xzzz_0, tg_xxxxyy_yyyy_0, tg_xxxxyy_yyyz_0, tg_xxxxyy_yyzz_0, \
                                         tg_xxxxyy_yzzz_0, tg_xxxxyy_zzzz_0, tg_xxxxyz_xxxx_0, tg_xxxxyz_xxxy_0, \
                                         tg_xxxxyz_xxxz_0, tg_xxxxyz_xxyy_0, tg_xxxxyz_xxyz_0, tg_xxxxyz_xxzz_0, \
                                         tg_xxxxyz_xyyy_0, tg_xxxxyz_xyyz_0, tg_xxxxyz_xyzz_0, tg_xxxxyz_xzzz_0, \
                                         tg_xxxxyz_yyyy_0, tg_xxxxyz_yyyz_0, tg_xxxxyz_yyzz_0, tg_xxxxyz_yzzz_0, \
                                         tg_xxxxyz_zzzz_0, tg_xxxxzz_xxxx_0, tg_xxxxzz_xxxy_0, tg_xxxxzz_xxxz_0, \
                                         tg_xxxxzz_xxyy_0, tg_xxxxzz_xxyz_0, tg_xxxxzz_xxzz_0, tg_xxxxzz_xyyy_0, \
                                         tg_xxxxzz_xyyz_0, tg_xxxxzz_xyzz_0, tg_xxxxzz_xzzz_0, tg_xxxxzz_yyyy_0, \
                                         tg_xxxxzz_yyyz_0, tg_xxxxzz_yyzz_0, tg_xxxxzz_yzzz_0, tg_xxxxzz_zzzz_0, \
                                         tg_xxxyy_xxxz_0, tg_xxxyy_xxxz_1, tg_xxxyy_xxyy_0, tg_xxxyy_xxyy_1, tg_xxxyy_xxyz_0, \
                                         tg_xxxyy_xxyz_1, tg_xxxyy_xxz_1, tg_xxxyy_xxzz_0, tg_xxxyy_xxzz_1, tg_xxxyy_xyy_1, \
                                         tg_xxxyy_xyyy_0, tg_xxxyy_xyyy_1, tg_xxxyy_xyyz_0, tg_xxxyy_xyyz_1, tg_xxxyy_xyz_1, \
                                         tg_xxxyy_xyzz_0, tg_xxxyy_xyzz_1, tg_xxxyy_xzz_1, tg_xxxyy_xzzz_0, tg_xxxyy_xzzz_1, \
                                         tg_xxxyy_yyy_1, tg_xxxyy_yyyy_0, tg_xxxyy_yyyy_1, tg_xxxyy_yyyz_0, tg_xxxyy_yyyz_1, \
                                         tg_xxxyy_yyz_1, tg_xxxyy_yyzz_0, tg_xxxyy_yyzz_1, tg_xxxyy_yzz_1, tg_xxxyy_yzzz_0, \
                                         tg_xxxyy_yzzz_1, tg_xxxyy_zzz_1, tg_xxxyy_zzzz_0, tg_xxxyy_zzzz_1, tg_xxxyyy_xxxx_0, \
                                         tg_xxxyyy_xxxy_0, tg_xxxyyy_xxxz_0, tg_xxxyyy_xxyy_0, tg_xxxyz_xxx_1, \
                                         tg_xxxyz_xxxx_0, tg_xxxyz_xxxx_1, tg_xxxyz_xxxy_0, tg_xxxyz_xxxy_1, tg_xxxyz_xxxz_0, \
                                         tg_xxxyz_xxxz_1, tg_xxxyz_xxy_1, tg_xxxyz_xxyy_0, tg_xxxyz_xxyy_1, tg_xxxyz_xxyz_0, \
                                         tg_xxxyz_xxyz_1, tg_xxxyz_xxz_1, tg_xxxyz_xxzz_0, tg_xxxyz_xxzz_1, tg_xxxyz_xyy_1, \
                                         tg_xxxyz_xyyy_0, tg_xxxyz_xyyy_1, tg_xxxyz_xyyz_0, tg_xxxyz_xyyz_1, tg_xxxyz_xyz_1, \
                                         tg_xxxyz_xyzz_0, tg_xxxyz_xyzz_1, tg_xxxyz_xzz_1, tg_xxxyz_xzzz_0, tg_xxxyz_xzzz_1, \
                                         tg_xxxyz_yyy_1, tg_xxxyz_yyyy_0, tg_xxxyz_yyyy_1, tg_xxxyz_yyyz_0, tg_xxxyz_yyyz_1, \
                                         tg_xxxyz_yyz_1, tg_xxxyz_yyzz_0, tg_xxxyz_yyzz_1, tg_xxxyz_yzz_1, tg_xxxyz_yzzz_0, \
                                         tg_xxxyz_yzzz_1, tg_xxxyz_zzz_1, tg_xxxyz_zzzz_0, tg_xxxyz_zzzz_1, tg_xxxzz_xxx_1, \
                                         tg_xxxzz_xxxx_0, tg_xxxzz_xxxx_1, tg_xxxzz_xxxy_0, tg_xxxzz_xxxy_1, tg_xxxzz_xxxz_0, \
                                         tg_xxxzz_xxxz_1, tg_xxxzz_xxy_1, tg_xxxzz_xxyy_0, tg_xxxzz_xxyy_1, tg_xxxzz_xxyz_0, \
                                         tg_xxxzz_xxyz_1, tg_xxxzz_xxz_1, tg_xxxzz_xxzz_0, tg_xxxzz_xxzz_1, tg_xxxzz_xyy_1, \
                                         tg_xxxzz_xyyy_0, tg_xxxzz_xyyy_1, tg_xxxzz_xyyz_0, tg_xxxzz_xyyz_1, tg_xxxzz_xyz_1, \
                                         tg_xxxzz_xyzz_0, tg_xxxzz_xyzz_1, tg_xxxzz_xzz_1, tg_xxxzz_xzzz_0, tg_xxxzz_xzzz_1, \
                                         tg_xxxzz_yyy_1, tg_xxxzz_yyyy_0, tg_xxxzz_yyyy_1, tg_xxxzz_yyyz_0, tg_xxxzz_yyyz_1, \
                                         tg_xxxzz_yyz_1, tg_xxxzz_yyzz_0, tg_xxxzz_yyzz_1, tg_xxxzz_yzz_1, tg_xxxzz_yzzz_0, \
                                         tg_xxxzz_yzzz_1, tg_xxxzz_zzz_1, tg_xxxzz_zzzz_0, tg_xxxzz_zzzz_1, tg_xxyy_xxxz_0, \
                                         tg_xxyy_xxxz_1, tg_xxyy_xxyy_0, tg_xxyy_xxyy_1, tg_xxyy_xxyz_0, tg_xxyy_xxyz_1, \
                                         tg_xxyy_xxzz_0, tg_xxyy_xxzz_1, tg_xxyy_xyyy_0, tg_xxyy_xyyy_1, tg_xxyy_xyyz_0, \
                                         tg_xxyy_xyyz_1, tg_xxyy_xyzz_0, tg_xxyy_xyzz_1, tg_xxyy_xzzz_0, tg_xxyy_xzzz_1, \
                                         tg_xxyy_yyyy_0, tg_xxyy_yyyy_1, tg_xxyy_yyyz_0, tg_xxyy_yyyz_1, tg_xxyy_yyzz_0, \
                                         tg_xxyy_yyzz_1, tg_xxyy_yzzz_0, tg_xxyy_yzzz_1, tg_xxyy_zzzz_0, tg_xxyy_zzzz_1, \
                                         tg_xxyyy_xxx_1, tg_xxyyy_xxxx_0, tg_xxyyy_xxxx_1, tg_xxyyy_xxxy_0, tg_xxyyy_xxxy_1, \
                                         tg_xxyyy_xxxz_0, tg_xxyyy_xxxz_1, tg_xxyyy_xxy_1, tg_xxyyy_xxyy_0, tg_xxyyy_xxyy_1, \
                                         tg_xxyyy_xxz_1, tg_xxyyy_xyy_1, tg_xxyz_xxxx_0, tg_xxyz_xxxx_1, tg_xxyz_xxxy_0, \
                                         tg_xxyz_xxxy_1, tg_xxyz_xxxz_0, tg_xxyz_xxxz_1, tg_xxyz_xxyy_0, tg_xxyz_xxyy_1, \
                                         tg_xxyz_xxyz_0, tg_xxyz_xxyz_1, tg_xxyz_xxzz_0, tg_xxyz_xxzz_1, tg_xxyz_xyyy_0, \
                                         tg_xxyz_xyyy_1, tg_xxyz_xyyz_0, tg_xxyz_xyyz_1, tg_xxyz_xyzz_0, tg_xxyz_xyzz_1, \
                                         tg_xxyz_xzzz_0, tg_xxyz_xzzz_1, tg_xxyz_yyyy_0, tg_xxyz_yyyy_1, tg_xxyz_yyyz_0, \
                                         tg_xxyz_yyyz_1, tg_xxyz_yyzz_0, tg_xxyz_yyzz_1, tg_xxyz_yzzz_0, tg_xxyz_yzzz_1, \
                                         tg_xxyz_zzzz_0, tg_xxyz_zzzz_1, tg_xxzz_xxxx_0, tg_xxzz_xxxx_1, tg_xxzz_xxxy_0, \
                                         tg_xxzz_xxxy_1, tg_xxzz_xxxz_0, tg_xxzz_xxxz_1, tg_xxzz_xxyy_0, tg_xxzz_xxyy_1, \
                                         tg_xxzz_xxyz_0, tg_xxzz_xxyz_1, tg_xxzz_xxzz_0, tg_xxzz_xxzz_1, tg_xxzz_xyyy_0, \
                                         tg_xxzz_xyyy_1, tg_xxzz_xyyz_0, tg_xxzz_xyyz_1, tg_xxzz_xyzz_0, tg_xxzz_xyzz_1, \
                                         tg_xxzz_xzzz_0, tg_xxzz_xzzz_1, tg_xxzz_yyyy_0, tg_xxzz_yyyy_1, tg_xxzz_yyyz_0, \
                                         tg_xxzz_yyyz_1, tg_xxzz_yyzz_0, tg_xxzz_yyzz_1, tg_xxzz_yzzz_0, tg_xxzz_yzzz_1, \
                                         tg_xxzz_zzzz_0, tg_xxzz_zzzz_1, tg_xyyy_xxxx_0, tg_xyyy_xxxx_1, tg_xyyy_xxxy_0, \
                                         tg_xyyy_xxxy_1, tg_xyyy_xxxz_0, tg_xyyy_xxxz_1, tg_xyyy_xxyy_0, tg_xyyy_xxyy_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxxxyy_xxxz_0[j] = pb_x * tg_xxxyy_xxxz_0[j] + wp_x[j] * tg_xxxyy_xxxz_1[j] + 1.5 * fl1_fx * tg_xxyy_xxxz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyy_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxxyy_xxz_1[j];

                    tg_xxxxyy_xxyy_0[j] = pb_x * tg_xxxyy_xxyy_0[j] + wp_x[j] * tg_xxxyy_xxyy_1[j] + 1.5 * fl1_fx * tg_xxyy_xxyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyy_xxyy_1[j] + fl1_fxn * tg_xxxyy_xyy_1[j];

                    tg_xxxxyy_xxyz_0[j] = pb_x * tg_xxxyy_xxyz_0[j] + wp_x[j] * tg_xxxyy_xxyz_1[j] + 1.5 * fl1_fx * tg_xxyy_xxyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyy_xxyz_1[j] + fl1_fxn * tg_xxxyy_xyz_1[j];

                    tg_xxxxyy_xxzz_0[j] = pb_x * tg_xxxyy_xxzz_0[j] + wp_x[j] * tg_xxxyy_xxzz_1[j] + 1.5 * fl1_fx * tg_xxyy_xxzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyy_xxzz_1[j] + fl1_fxn * tg_xxxyy_xzz_1[j];

                    tg_xxxxyy_xyyy_0[j] = pb_x * tg_xxxyy_xyyy_0[j] + wp_x[j] * tg_xxxyy_xyyy_1[j] + 1.5 * fl1_fx * tg_xxyy_xyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyy_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxxyy_yyy_1[j];

                    tg_xxxxyy_xyyz_0[j] = pb_x * tg_xxxyy_xyyz_0[j] + wp_x[j] * tg_xxxyy_xyyz_1[j] + 1.5 * fl1_fx * tg_xxyy_xyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyy_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxxyy_yyz_1[j];

                    tg_xxxxyy_xyzz_0[j] = pb_x * tg_xxxyy_xyzz_0[j] + wp_x[j] * tg_xxxyy_xyzz_1[j] + 1.5 * fl1_fx * tg_xxyy_xyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyy_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxxyy_yzz_1[j];

                    tg_xxxxyy_xzzz_0[j] = pb_x * tg_xxxyy_xzzz_0[j] + wp_x[j] * tg_xxxyy_xzzz_1[j] + 1.5 * fl1_fx * tg_xxyy_xzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyy_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxxyy_zzz_1[j];

                    tg_xxxxyy_yyyy_0[j] = pb_x * tg_xxxyy_yyyy_0[j] + wp_x[j] * tg_xxxyy_yyyy_1[j] + 1.5 * fl1_fx * tg_xxyy_yyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyy_yyyy_1[j];

                    tg_xxxxyy_yyyz_0[j] = pb_x * tg_xxxyy_yyyz_0[j] + wp_x[j] * tg_xxxyy_yyyz_1[j] + 1.5 * fl1_fx * tg_xxyy_yyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyy_yyyz_1[j];

                    tg_xxxxyy_yyzz_0[j] = pb_x * tg_xxxyy_yyzz_0[j] + wp_x[j] * tg_xxxyy_yyzz_1[j] + 1.5 * fl1_fx * tg_xxyy_yyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyy_yyzz_1[j];

                    tg_xxxxyy_yzzz_0[j] = pb_x * tg_xxxyy_yzzz_0[j] + wp_x[j] * tg_xxxyy_yzzz_1[j] + 1.5 * fl1_fx * tg_xxyy_yzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyy_yzzz_1[j];

                    tg_xxxxyy_zzzz_0[j] = pb_x * tg_xxxyy_zzzz_0[j] + wp_x[j] * tg_xxxyy_zzzz_1[j] + 1.5 * fl1_fx * tg_xxyy_zzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyy_zzzz_1[j];

                    tg_xxxxyz_xxxx_0[j] = pb_x * tg_xxxyz_xxxx_0[j] + wp_x[j] * tg_xxxyz_xxxx_1[j] + 1.5 * fl1_fx * tg_xxyz_xxxx_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxxyz_xxx_1[j];

                    tg_xxxxyz_xxxy_0[j] = pb_x * tg_xxxyz_xxxy_0[j] + wp_x[j] * tg_xxxyz_xxxy_1[j] + 1.5 * fl1_fx * tg_xxyz_xxxy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxxyz_xxy_1[j];

                    tg_xxxxyz_xxxz_0[j] = pb_x * tg_xxxyz_xxxz_0[j] + wp_x[j] * tg_xxxyz_xxxz_1[j] + 1.5 * fl1_fx * tg_xxyz_xxxz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxxyz_xxz_1[j];

                    tg_xxxxyz_xxyy_0[j] = pb_x * tg_xxxyz_xxyy_0[j] + wp_x[j] * tg_xxxyz_xxyy_1[j] + 1.5 * fl1_fx * tg_xxyz_xxyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyz_xxyy_1[j] + fl1_fxn * tg_xxxyz_xyy_1[j];

                    tg_xxxxyz_xxyz_0[j] = pb_x * tg_xxxyz_xxyz_0[j] + wp_x[j] * tg_xxxyz_xxyz_1[j] + 1.5 * fl1_fx * tg_xxyz_xxyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyz_xxyz_1[j] + fl1_fxn * tg_xxxyz_xyz_1[j];

                    tg_xxxxyz_xxzz_0[j] = pb_x * tg_xxxyz_xxzz_0[j] + wp_x[j] * tg_xxxyz_xxzz_1[j] + 1.5 * fl1_fx * tg_xxyz_xxzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyz_xxzz_1[j] + fl1_fxn * tg_xxxyz_xzz_1[j];

                    tg_xxxxyz_xyyy_0[j] = pb_x * tg_xxxyz_xyyy_0[j] + wp_x[j] * tg_xxxyz_xyyy_1[j] + 1.5 * fl1_fx * tg_xxyz_xyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxxyz_yyy_1[j];

                    tg_xxxxyz_xyyz_0[j] = pb_x * tg_xxxyz_xyyz_0[j] + wp_x[j] * tg_xxxyz_xyyz_1[j] + 1.5 * fl1_fx * tg_xxyz_xyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxxyz_yyz_1[j];

                    tg_xxxxyz_xyzz_0[j] = pb_x * tg_xxxyz_xyzz_0[j] + wp_x[j] * tg_xxxyz_xyzz_1[j] + 1.5 * fl1_fx * tg_xxyz_xyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxxyz_yzz_1[j];

                    tg_xxxxyz_xzzz_0[j] = pb_x * tg_xxxyz_xzzz_0[j] + wp_x[j] * tg_xxxyz_xzzz_1[j] + 1.5 * fl1_fx * tg_xxyz_xzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxxyz_zzz_1[j];

                    tg_xxxxyz_yyyy_0[j] = pb_x * tg_xxxyz_yyyy_0[j] + wp_x[j] * tg_xxxyz_yyyy_1[j] + 1.5 * fl1_fx * tg_xxyz_yyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyz_yyyy_1[j];

                    tg_xxxxyz_yyyz_0[j] = pb_x * tg_xxxyz_yyyz_0[j] + wp_x[j] * tg_xxxyz_yyyz_1[j] + 1.5 * fl1_fx * tg_xxyz_yyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyz_yyyz_1[j];

                    tg_xxxxyz_yyzz_0[j] = pb_x * tg_xxxyz_yyzz_0[j] + wp_x[j] * tg_xxxyz_yyzz_1[j] + 1.5 * fl1_fx * tg_xxyz_yyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyz_yyzz_1[j];

                    tg_xxxxyz_yzzz_0[j] = pb_x * tg_xxxyz_yzzz_0[j] + wp_x[j] * tg_xxxyz_yzzz_1[j] + 1.5 * fl1_fx * tg_xxyz_yzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyz_yzzz_1[j];

                    tg_xxxxyz_zzzz_0[j] = pb_x * tg_xxxyz_zzzz_0[j] + wp_x[j] * tg_xxxyz_zzzz_1[j] + 1.5 * fl1_fx * tg_xxyz_zzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyz_zzzz_1[j];

                    tg_xxxxzz_xxxx_0[j] = pb_x * tg_xxxzz_xxxx_0[j] + wp_x[j] * tg_xxxzz_xxxx_1[j] + 1.5 * fl1_fx * tg_xxzz_xxxx_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxxzz_xxx_1[j];

                    tg_xxxxzz_xxxy_0[j] = pb_x * tg_xxxzz_xxxy_0[j] + wp_x[j] * tg_xxxzz_xxxy_1[j] + 1.5 * fl1_fx * tg_xxzz_xxxy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxxzz_xxy_1[j];

                    tg_xxxxzz_xxxz_0[j] = pb_x * tg_xxxzz_xxxz_0[j] + wp_x[j] * tg_xxxzz_xxxz_1[j] + 1.5 * fl1_fx * tg_xxzz_xxxz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxxzz_xxz_1[j];

                    tg_xxxxzz_xxyy_0[j] = pb_x * tg_xxxzz_xxyy_0[j] + wp_x[j] * tg_xxxzz_xxyy_1[j] + 1.5 * fl1_fx * tg_xxzz_xxyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzz_xxyy_1[j] + fl1_fxn * tg_xxxzz_xyy_1[j];

                    tg_xxxxzz_xxyz_0[j] = pb_x * tg_xxxzz_xxyz_0[j] + wp_x[j] * tg_xxxzz_xxyz_1[j] + 1.5 * fl1_fx * tg_xxzz_xxyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzz_xxyz_1[j] + fl1_fxn * tg_xxxzz_xyz_1[j];

                    tg_xxxxzz_xxzz_0[j] = pb_x * tg_xxxzz_xxzz_0[j] + wp_x[j] * tg_xxxzz_xxzz_1[j] + 1.5 * fl1_fx * tg_xxzz_xxzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzz_xxzz_1[j] + fl1_fxn * tg_xxxzz_xzz_1[j];

                    tg_xxxxzz_xyyy_0[j] = pb_x * tg_xxxzz_xyyy_0[j] + wp_x[j] * tg_xxxzz_xyyy_1[j] + 1.5 * fl1_fx * tg_xxzz_xyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxxzz_yyy_1[j];

                    tg_xxxxzz_xyyz_0[j] = pb_x * tg_xxxzz_xyyz_0[j] + wp_x[j] * tg_xxxzz_xyyz_1[j] + 1.5 * fl1_fx * tg_xxzz_xyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxxzz_yyz_1[j];

                    tg_xxxxzz_xyzz_0[j] = pb_x * tg_xxxzz_xyzz_0[j] + wp_x[j] * tg_xxxzz_xyzz_1[j] + 1.5 * fl1_fx * tg_xxzz_xyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxxzz_yzz_1[j];

                    tg_xxxxzz_xzzz_0[j] = pb_x * tg_xxxzz_xzzz_0[j] + wp_x[j] * tg_xxxzz_xzzz_1[j] + 1.5 * fl1_fx * tg_xxzz_xzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxxzz_zzz_1[j];

                    tg_xxxxzz_yyyy_0[j] = pb_x * tg_xxxzz_yyyy_0[j] + wp_x[j] * tg_xxxzz_yyyy_1[j] + 1.5 * fl1_fx * tg_xxzz_yyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzz_yyyy_1[j];

                    tg_xxxxzz_yyyz_0[j] = pb_x * tg_xxxzz_yyyz_0[j] + wp_x[j] * tg_xxxzz_yyyz_1[j] + 1.5 * fl1_fx * tg_xxzz_yyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzz_yyyz_1[j];

                    tg_xxxxzz_yyzz_0[j] = pb_x * tg_xxxzz_yyzz_0[j] + wp_x[j] * tg_xxxzz_yyzz_1[j] + 1.5 * fl1_fx * tg_xxzz_yyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzz_yyzz_1[j];

                    tg_xxxxzz_yzzz_0[j] = pb_x * tg_xxxzz_yzzz_0[j] + wp_x[j] * tg_xxxzz_yzzz_1[j] + 1.5 * fl1_fx * tg_xxzz_yzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzz_yzzz_1[j];

                    tg_xxxxzz_zzzz_0[j] = pb_x * tg_xxxzz_zzzz_0[j] + wp_x[j] * tg_xxxzz_zzzz_1[j] + 1.5 * fl1_fx * tg_xxzz_zzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzz_zzzz_1[j];

                    tg_xxxyyy_xxxx_0[j] = pb_x * tg_xxyyy_xxxx_0[j] + wp_x[j] * tg_xxyyy_xxxx_1[j] + fl1_fx * tg_xyyy_xxxx_0[j] - fl1_fx * fl1_fza * tg_xyyy_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxyyy_xxx_1[j];

                    tg_xxxyyy_xxxy_0[j] = pb_x * tg_xxyyy_xxxy_0[j] + wp_x[j] * tg_xxyyy_xxxy_1[j] + fl1_fx * tg_xyyy_xxxy_0[j] - fl1_fx * fl1_fza * tg_xyyy_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxyyy_xxy_1[j];

                    tg_xxxyyy_xxxz_0[j] = pb_x * tg_xxyyy_xxxz_0[j] + wp_x[j] * tg_xxyyy_xxxz_1[j] + fl1_fx * tg_xyyy_xxxz_0[j] - fl1_fx * fl1_fza * tg_xyyy_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxyyy_xxz_1[j];

                    tg_xxxyyy_xxyy_0[j] = pb_x * tg_xxyyy_xxyy_0[j] + wp_x[j] * tg_xxyyy_xxyy_1[j] + fl1_fx * tg_xyyy_xxyy_0[j] - fl1_fx * fl1_fza * tg_xyyy_xxyy_1[j] + fl1_fxn * tg_xxyyy_xyy_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSISG_94_141(      CMemBlock2D<double>& primBuffer,
                                        const CRecursionMap&       recursionMap,
                                        const CMemBlock2D<double>& osFactors,
                                        const CMemBlock2D<double>& wpDistances,
                                        const CGtoPairsBlock&      braGtoPairsBlock,
                                        const CGtoPairsBlock&      ketGtoPairsBlock,
                                        const int32_t              nKetPrimPairs,
                                        const int32_t              iContrPair)
    {
        // Batch of Integrals (94,141)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {6, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_6_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_6_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_5_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_5_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_4_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_5_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_xxyyy_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 94); 

                auto tg_xxyyy_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 95); 

                auto tg_xxyyy_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 96); 

                auto tg_xxyyy_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 97); 

                auto tg_xxyyy_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 98); 

                auto tg_xxyyy_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 99); 

                auto tg_xxyyy_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 100); 

                auto tg_xxyyy_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 101); 

                auto tg_xxyyy_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 102); 

                auto tg_xxyyy_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 103); 

                auto tg_xxyyy_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 104); 

                auto tg_xxyyz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 105); 

                auto tg_xxyyz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 106); 

                auto tg_xxyyz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 107); 

                auto tg_xxyyz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 108); 

                auto tg_xxyyz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 109); 

                auto tg_xxyyz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 110); 

                auto tg_xxyyz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 111); 

                auto tg_xxyyz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 112); 

                auto tg_xxyyz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 113); 

                auto tg_xxyyz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 114); 

                auto tg_xxyyz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 115); 

                auto tg_xxyyz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 116); 

                auto tg_xxyyz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 117); 

                auto tg_xxyyz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 118); 

                auto tg_xxyyz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 119); 

                auto tg_xxyzz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 120); 

                auto tg_xxyzz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 121); 

                auto tg_xxyzz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 122); 

                auto tg_xxyzz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 123); 

                auto tg_xxyzz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 124); 

                auto tg_xxyzz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 125); 

                auto tg_xxyzz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 126); 

                auto tg_xxyzz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 127); 

                auto tg_xxyzz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 128); 

                auto tg_xxyzz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 129); 

                auto tg_xxyzz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 130); 

                auto tg_xxyzz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 131); 

                auto tg_xxyzz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 132); 

                auto tg_xxyzz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 133); 

                auto tg_xxyzz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 134); 

                auto tg_xxzzz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 135); 

                auto tg_xxzzz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 136); 

                auto tg_xxzzz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 137); 

                auto tg_xxzzz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 138); 

                auto tg_xxzzz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 139); 

                auto tg_xxzzz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 140); 

                auto tg_xxyyy_xxyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 94); 

                auto tg_xxyyy_xxzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 95); 

                auto tg_xxyyy_xyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 96); 

                auto tg_xxyyy_xyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 97); 

                auto tg_xxyyy_xyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 98); 

                auto tg_xxyyy_xzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 99); 

                auto tg_xxyyy_yyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 100); 

                auto tg_xxyyy_yyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 101); 

                auto tg_xxyyy_yyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 102); 

                auto tg_xxyyy_yzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 103); 

                auto tg_xxyyy_zzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 104); 

                auto tg_xxyyz_xxxx_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 105); 

                auto tg_xxyyz_xxxy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 106); 

                auto tg_xxyyz_xxxz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 107); 

                auto tg_xxyyz_xxyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 108); 

                auto tg_xxyyz_xxyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 109); 

                auto tg_xxyyz_xxzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 110); 

                auto tg_xxyyz_xyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 111); 

                auto tg_xxyyz_xyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 112); 

                auto tg_xxyyz_xyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 113); 

                auto tg_xxyyz_xzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 114); 

                auto tg_xxyyz_yyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 115); 

                auto tg_xxyyz_yyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 116); 

                auto tg_xxyyz_yyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 117); 

                auto tg_xxyyz_yzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 118); 

                auto tg_xxyyz_zzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 119); 

                auto tg_xxyzz_xxxx_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 120); 

                auto tg_xxyzz_xxxy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 121); 

                auto tg_xxyzz_xxxz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 122); 

                auto tg_xxyzz_xxyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 123); 

                auto tg_xxyzz_xxyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 124); 

                auto tg_xxyzz_xxzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 125); 

                auto tg_xxyzz_xyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 126); 

                auto tg_xxyzz_xyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 127); 

                auto tg_xxyzz_xyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 128); 

                auto tg_xxyzz_xzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 129); 

                auto tg_xxyzz_yyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 130); 

                auto tg_xxyzz_yyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 131); 

                auto tg_xxyzz_yyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 132); 

                auto tg_xxyzz_yzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 133); 

                auto tg_xxyzz_zzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 134); 

                auto tg_xxzzz_xxxx_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 135); 

                auto tg_xxzzz_xxxy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 136); 

                auto tg_xxzzz_xxxz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 137); 

                auto tg_xxzzz_xxyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 138); 

                auto tg_xxzzz_xxyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 139); 

                auto tg_xxzzz_xxzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 140); 

                auto tg_xyyy_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 94); 

                auto tg_xyyy_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 95); 

                auto tg_xyyy_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 96); 

                auto tg_xyyy_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 97); 

                auto tg_xyyy_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 98); 

                auto tg_xyyy_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 99); 

                auto tg_xyyy_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 100); 

                auto tg_xyyy_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 101); 

                auto tg_xyyy_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 102); 

                auto tg_xyyy_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 103); 

                auto tg_xyyy_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 104); 

                auto tg_xyyz_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 105); 

                auto tg_xyyz_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 106); 

                auto tg_xyyz_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 107); 

                auto tg_xyyz_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 108); 

                auto tg_xyyz_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 109); 

                auto tg_xyyz_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 110); 

                auto tg_xyyz_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 111); 

                auto tg_xyyz_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 112); 

                auto tg_xyyz_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 113); 

                auto tg_xyyz_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 114); 

                auto tg_xyyz_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 115); 

                auto tg_xyyz_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 116); 

                auto tg_xyyz_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 117); 

                auto tg_xyyz_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 118); 

                auto tg_xyyz_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 119); 

                auto tg_xyzz_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 120); 

                auto tg_xyzz_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 121); 

                auto tg_xyzz_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 122); 

                auto tg_xyzz_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 123); 

                auto tg_xyzz_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 124); 

                auto tg_xyzz_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 125); 

                auto tg_xyzz_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 126); 

                auto tg_xyzz_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 127); 

                auto tg_xyzz_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 128); 

                auto tg_xyzz_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 129); 

                auto tg_xyzz_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 130); 

                auto tg_xyzz_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 131); 

                auto tg_xyzz_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 132); 

                auto tg_xyzz_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 133); 

                auto tg_xyzz_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 134); 

                auto tg_xzzz_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 135); 

                auto tg_xzzz_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 136); 

                auto tg_xzzz_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 137); 

                auto tg_xzzz_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 138); 

                auto tg_xzzz_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 139); 

                auto tg_xzzz_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 140); 

                auto tg_xyyy_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 94); 

                auto tg_xyyy_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 95); 

                auto tg_xyyy_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 96); 

                auto tg_xyyy_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 97); 

                auto tg_xyyy_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 98); 

                auto tg_xyyy_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 99); 

                auto tg_xyyy_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 100); 

                auto tg_xyyy_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 101); 

                auto tg_xyyy_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 102); 

                auto tg_xyyy_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 103); 

                auto tg_xyyy_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 104); 

                auto tg_xyyz_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 105); 

                auto tg_xyyz_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 106); 

                auto tg_xyyz_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 107); 

                auto tg_xyyz_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 108); 

                auto tg_xyyz_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 109); 

                auto tg_xyyz_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 110); 

                auto tg_xyyz_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 111); 

                auto tg_xyyz_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 112); 

                auto tg_xyyz_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 113); 

                auto tg_xyyz_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 114); 

                auto tg_xyyz_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 115); 

                auto tg_xyyz_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 116); 

                auto tg_xyyz_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 117); 

                auto tg_xyyz_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 118); 

                auto tg_xyyz_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 119); 

                auto tg_xyzz_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 120); 

                auto tg_xyzz_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 121); 

                auto tg_xyzz_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 122); 

                auto tg_xyzz_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 123); 

                auto tg_xyzz_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 124); 

                auto tg_xyzz_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 125); 

                auto tg_xyzz_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 126); 

                auto tg_xyzz_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 127); 

                auto tg_xyzz_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 128); 

                auto tg_xyzz_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 129); 

                auto tg_xyzz_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 130); 

                auto tg_xyzz_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 131); 

                auto tg_xyzz_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 132); 

                auto tg_xyzz_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 133); 

                auto tg_xyzz_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 134); 

                auto tg_xzzz_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 135); 

                auto tg_xzzz_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 136); 

                auto tg_xzzz_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 137); 

                auto tg_xzzz_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 138); 

                auto tg_xzzz_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 139); 

                auto tg_xzzz_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 140); 

                auto tg_xxyyy_xyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 64); 

                auto tg_xxyyy_xzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 65); 

                auto tg_xxyyy_yyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 66); 

                auto tg_xxyyy_yyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 67); 

                auto tg_xxyyy_yzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 68); 

                auto tg_xxyyy_zzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 69); 

                auto tg_xxyyz_xxx_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 70); 

                auto tg_xxyyz_xxy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 71); 

                auto tg_xxyyz_xxz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 72); 

                auto tg_xxyyz_xyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 73); 

                auto tg_xxyyz_xyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 74); 

                auto tg_xxyyz_xzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 75); 

                auto tg_xxyyz_yyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 76); 

                auto tg_xxyyz_yyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 77); 

                auto tg_xxyyz_yzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 78); 

                auto tg_xxyyz_zzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 79); 

                auto tg_xxyzz_xxx_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 80); 

                auto tg_xxyzz_xxy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 81); 

                auto tg_xxyzz_xxz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 82); 

                auto tg_xxyzz_xyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 83); 

                auto tg_xxyzz_xyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 84); 

                auto tg_xxyzz_xzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 85); 

                auto tg_xxyzz_yyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 86); 

                auto tg_xxyzz_yyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 87); 

                auto tg_xxyzz_yzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 88); 

                auto tg_xxyzz_zzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 89); 

                auto tg_xxzzz_xxx_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 90); 

                auto tg_xxzzz_xxy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 91); 

                auto tg_xxzzz_xxz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 92); 

                auto tg_xxzzz_xyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 93); 

                auto tg_xxzzz_xyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 94); 

                auto tg_xxzzz_xzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 95); 

                // set up pointers to integrals

                auto tg_xxxyyy_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 94); 

                auto tg_xxxyyy_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 95); 

                auto tg_xxxyyy_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 96); 

                auto tg_xxxyyy_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 97); 

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

                // Batch of Integrals (94,141)

                #pragma omp simd aligned(fxn, fza, tg_xxxyyy_xxyz_0, tg_xxxyyy_xxzz_0, tg_xxxyyy_xyyy_0, \
                                         tg_xxxyyy_xyyz_0, tg_xxxyyy_xyzz_0, tg_xxxyyy_xzzz_0, tg_xxxyyy_yyyy_0, \
                                         tg_xxxyyy_yyyz_0, tg_xxxyyy_yyzz_0, tg_xxxyyy_yzzz_0, tg_xxxyyy_zzzz_0, \
                                         tg_xxxyyz_xxxx_0, tg_xxxyyz_xxxy_0, tg_xxxyyz_xxxz_0, tg_xxxyyz_xxyy_0, \
                                         tg_xxxyyz_xxyz_0, tg_xxxyyz_xxzz_0, tg_xxxyyz_xyyy_0, tg_xxxyyz_xyyz_0, \
                                         tg_xxxyyz_xyzz_0, tg_xxxyyz_xzzz_0, tg_xxxyyz_yyyy_0, tg_xxxyyz_yyyz_0, \
                                         tg_xxxyyz_yyzz_0, tg_xxxyyz_yzzz_0, tg_xxxyyz_zzzz_0, tg_xxxyzz_xxxx_0, \
                                         tg_xxxyzz_xxxy_0, tg_xxxyzz_xxxz_0, tg_xxxyzz_xxyy_0, tg_xxxyzz_xxyz_0, \
                                         tg_xxxyzz_xxzz_0, tg_xxxyzz_xyyy_0, tg_xxxyzz_xyyz_0, tg_xxxyzz_xyzz_0, \
                                         tg_xxxyzz_xzzz_0, tg_xxxyzz_yyyy_0, tg_xxxyzz_yyyz_0, tg_xxxyzz_yyzz_0, \
                                         tg_xxxyzz_yzzz_0, tg_xxxyzz_zzzz_0, tg_xxxzzz_xxxx_0, tg_xxxzzz_xxxy_0, \
                                         tg_xxxzzz_xxxz_0, tg_xxxzzz_xxyy_0, tg_xxxzzz_xxyz_0, tg_xxxzzz_xxzz_0, \
                                         tg_xxyyy_xxyz_0, tg_xxyyy_xxyz_1, tg_xxyyy_xxzz_0, tg_xxyyy_xxzz_1, tg_xxyyy_xyyy_0, \
                                         tg_xxyyy_xyyy_1, tg_xxyyy_xyyz_0, tg_xxyyy_xyyz_1, tg_xxyyy_xyz_1, tg_xxyyy_xyzz_0, \
                                         tg_xxyyy_xyzz_1, tg_xxyyy_xzz_1, tg_xxyyy_xzzz_0, tg_xxyyy_xzzz_1, tg_xxyyy_yyy_1, \
                                         tg_xxyyy_yyyy_0, tg_xxyyy_yyyy_1, tg_xxyyy_yyyz_0, tg_xxyyy_yyyz_1, tg_xxyyy_yyz_1, \
                                         tg_xxyyy_yyzz_0, tg_xxyyy_yyzz_1, tg_xxyyy_yzz_1, tg_xxyyy_yzzz_0, tg_xxyyy_yzzz_1, \
                                         tg_xxyyy_zzz_1, tg_xxyyy_zzzz_0, tg_xxyyy_zzzz_1, tg_xxyyz_xxx_1, tg_xxyyz_xxxx_0, \
                                         tg_xxyyz_xxxx_1, tg_xxyyz_xxxy_0, tg_xxyyz_xxxy_1, tg_xxyyz_xxxz_0, tg_xxyyz_xxxz_1, \
                                         tg_xxyyz_xxy_1, tg_xxyyz_xxyy_0, tg_xxyyz_xxyy_1, tg_xxyyz_xxyz_0, tg_xxyyz_xxyz_1, \
                                         tg_xxyyz_xxz_1, tg_xxyyz_xxzz_0, tg_xxyyz_xxzz_1, tg_xxyyz_xyy_1, tg_xxyyz_xyyy_0, \
                                         tg_xxyyz_xyyy_1, tg_xxyyz_xyyz_0, tg_xxyyz_xyyz_1, tg_xxyyz_xyz_1, tg_xxyyz_xyzz_0, \
                                         tg_xxyyz_xyzz_1, tg_xxyyz_xzz_1, tg_xxyyz_xzzz_0, tg_xxyyz_xzzz_1, tg_xxyyz_yyy_1, \
                                         tg_xxyyz_yyyy_0, tg_xxyyz_yyyy_1, tg_xxyyz_yyyz_0, tg_xxyyz_yyyz_1, tg_xxyyz_yyz_1, \
                                         tg_xxyyz_yyzz_0, tg_xxyyz_yyzz_1, tg_xxyyz_yzz_1, tg_xxyyz_yzzz_0, tg_xxyyz_yzzz_1, \
                                         tg_xxyyz_zzz_1, tg_xxyyz_zzzz_0, tg_xxyyz_zzzz_1, tg_xxyzz_xxx_1, tg_xxyzz_xxxx_0, \
                                         tg_xxyzz_xxxx_1, tg_xxyzz_xxxy_0, tg_xxyzz_xxxy_1, tg_xxyzz_xxxz_0, tg_xxyzz_xxxz_1, \
                                         tg_xxyzz_xxy_1, tg_xxyzz_xxyy_0, tg_xxyzz_xxyy_1, tg_xxyzz_xxyz_0, tg_xxyzz_xxyz_1, \
                                         tg_xxyzz_xxz_1, tg_xxyzz_xxzz_0, tg_xxyzz_xxzz_1, tg_xxyzz_xyy_1, tg_xxyzz_xyyy_0, \
                                         tg_xxyzz_xyyy_1, tg_xxyzz_xyyz_0, tg_xxyzz_xyyz_1, tg_xxyzz_xyz_1, tg_xxyzz_xyzz_0, \
                                         tg_xxyzz_xyzz_1, tg_xxyzz_xzz_1, tg_xxyzz_xzzz_0, tg_xxyzz_xzzz_1, tg_xxyzz_yyy_1, \
                                         tg_xxyzz_yyyy_0, tg_xxyzz_yyyy_1, tg_xxyzz_yyyz_0, tg_xxyzz_yyyz_1, tg_xxyzz_yyz_1, \
                                         tg_xxyzz_yyzz_0, tg_xxyzz_yyzz_1, tg_xxyzz_yzz_1, tg_xxyzz_yzzz_0, tg_xxyzz_yzzz_1, \
                                         tg_xxyzz_zzz_1, tg_xxyzz_zzzz_0, tg_xxyzz_zzzz_1, tg_xxzzz_xxx_1, tg_xxzzz_xxxx_0, \
                                         tg_xxzzz_xxxx_1, tg_xxzzz_xxxy_0, tg_xxzzz_xxxy_1, tg_xxzzz_xxxz_0, tg_xxzzz_xxxz_1, \
                                         tg_xxzzz_xxy_1, tg_xxzzz_xxyy_0, tg_xxzzz_xxyy_1, tg_xxzzz_xxyz_0, tg_xxzzz_xxyz_1, \
                                         tg_xxzzz_xxz_1, tg_xxzzz_xxzz_0, tg_xxzzz_xxzz_1, tg_xxzzz_xyy_1, tg_xxzzz_xyz_1, \
                                         tg_xxzzz_xzz_1, tg_xyyy_xxyz_0, tg_xyyy_xxyz_1, tg_xyyy_xxzz_0, tg_xyyy_xxzz_1, \
                                         tg_xyyy_xyyy_0, tg_xyyy_xyyy_1, tg_xyyy_xyyz_0, tg_xyyy_xyyz_1, tg_xyyy_xyzz_0, \
                                         tg_xyyy_xyzz_1, tg_xyyy_xzzz_0, tg_xyyy_xzzz_1, tg_xyyy_yyyy_0, tg_xyyy_yyyy_1, \
                                         tg_xyyy_yyyz_0, tg_xyyy_yyyz_1, tg_xyyy_yyzz_0, tg_xyyy_yyzz_1, tg_xyyy_yzzz_0, \
                                         tg_xyyy_yzzz_1, tg_xyyy_zzzz_0, tg_xyyy_zzzz_1, tg_xyyz_xxxx_0, tg_xyyz_xxxx_1, \
                                         tg_xyyz_xxxy_0, tg_xyyz_xxxy_1, tg_xyyz_xxxz_0, tg_xyyz_xxxz_1, tg_xyyz_xxyy_0, \
                                         tg_xyyz_xxyy_1, tg_xyyz_xxyz_0, tg_xyyz_xxyz_1, tg_xyyz_xxzz_0, tg_xyyz_xxzz_1, \
                                         tg_xyyz_xyyy_0, tg_xyyz_xyyy_1, tg_xyyz_xyyz_0, tg_xyyz_xyyz_1, tg_xyyz_xyzz_0, \
                                         tg_xyyz_xyzz_1, tg_xyyz_xzzz_0, tg_xyyz_xzzz_1, tg_xyyz_yyyy_0, tg_xyyz_yyyy_1, \
                                         tg_xyyz_yyyz_0, tg_xyyz_yyyz_1, tg_xyyz_yyzz_0, tg_xyyz_yyzz_1, tg_xyyz_yzzz_0, \
                                         tg_xyyz_yzzz_1, tg_xyyz_zzzz_0, tg_xyyz_zzzz_1, tg_xyzz_xxxx_0, tg_xyzz_xxxx_1, \
                                         tg_xyzz_xxxy_0, tg_xyzz_xxxy_1, tg_xyzz_xxxz_0, tg_xyzz_xxxz_1, tg_xyzz_xxyy_0, \
                                         tg_xyzz_xxyy_1, tg_xyzz_xxyz_0, tg_xyzz_xxyz_1, tg_xyzz_xxzz_0, tg_xyzz_xxzz_1, \
                                         tg_xyzz_xyyy_0, tg_xyzz_xyyy_1, tg_xyzz_xyyz_0, tg_xyzz_xyyz_1, tg_xyzz_xyzz_0, \
                                         tg_xyzz_xyzz_1, tg_xyzz_xzzz_0, tg_xyzz_xzzz_1, tg_xyzz_yyyy_0, tg_xyzz_yyyy_1, \
                                         tg_xyzz_yyyz_0, tg_xyzz_yyyz_1, tg_xyzz_yyzz_0, tg_xyzz_yyzz_1, tg_xyzz_yzzz_0, \
                                         tg_xyzz_yzzz_1, tg_xyzz_zzzz_0, tg_xyzz_zzzz_1, tg_xzzz_xxxx_0, tg_xzzz_xxxx_1, \
                                         tg_xzzz_xxxy_0, tg_xzzz_xxxy_1, tg_xzzz_xxxz_0, tg_xzzz_xxxz_1, tg_xzzz_xxyy_0, \
                                         tg_xzzz_xxyy_1, tg_xzzz_xxyz_0, tg_xzzz_xxyz_1, tg_xzzz_xxzz_0, tg_xzzz_xxzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxxyyy_xxyz_0[j] = pb_x * tg_xxyyy_xxyz_0[j] + wp_x[j] * tg_xxyyy_xxyz_1[j] + fl1_fx * tg_xyyy_xxyz_0[j] - fl1_fx * fl1_fza * tg_xyyy_xxyz_1[j] + fl1_fxn * tg_xxyyy_xyz_1[j];

                    tg_xxxyyy_xxzz_0[j] = pb_x * tg_xxyyy_xxzz_0[j] + wp_x[j] * tg_xxyyy_xxzz_1[j] + fl1_fx * tg_xyyy_xxzz_0[j] - fl1_fx * fl1_fza * tg_xyyy_xxzz_1[j] + fl1_fxn * tg_xxyyy_xzz_1[j];

                    tg_xxxyyy_xyyy_0[j] = pb_x * tg_xxyyy_xyyy_0[j] + wp_x[j] * tg_xxyyy_xyyy_1[j] + fl1_fx * tg_xyyy_xyyy_0[j] - fl1_fx * fl1_fza * tg_xyyy_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxyyy_yyy_1[j];

                    tg_xxxyyy_xyyz_0[j] = pb_x * tg_xxyyy_xyyz_0[j] + wp_x[j] * tg_xxyyy_xyyz_1[j] + fl1_fx * tg_xyyy_xyyz_0[j] - fl1_fx * fl1_fza * tg_xyyy_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxyyy_yyz_1[j];

                    tg_xxxyyy_xyzz_0[j] = pb_x * tg_xxyyy_xyzz_0[j] + wp_x[j] * tg_xxyyy_xyzz_1[j] + fl1_fx * tg_xyyy_xyzz_0[j] - fl1_fx * fl1_fza * tg_xyyy_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxyyy_yzz_1[j];

                    tg_xxxyyy_xzzz_0[j] = pb_x * tg_xxyyy_xzzz_0[j] + wp_x[j] * tg_xxyyy_xzzz_1[j] + fl1_fx * tg_xyyy_xzzz_0[j] - fl1_fx * fl1_fza * tg_xyyy_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxyyy_zzz_1[j];

                    tg_xxxyyy_yyyy_0[j] = pb_x * tg_xxyyy_yyyy_0[j] + wp_x[j] * tg_xxyyy_yyyy_1[j] + fl1_fx * tg_xyyy_yyyy_0[j] - fl1_fx * fl1_fza * tg_xyyy_yyyy_1[j];

                    tg_xxxyyy_yyyz_0[j] = pb_x * tg_xxyyy_yyyz_0[j] + wp_x[j] * tg_xxyyy_yyyz_1[j] + fl1_fx * tg_xyyy_yyyz_0[j] - fl1_fx * fl1_fza * tg_xyyy_yyyz_1[j];

                    tg_xxxyyy_yyzz_0[j] = pb_x * tg_xxyyy_yyzz_0[j] + wp_x[j] * tg_xxyyy_yyzz_1[j] + fl1_fx * tg_xyyy_yyzz_0[j] - fl1_fx * fl1_fza * tg_xyyy_yyzz_1[j];

                    tg_xxxyyy_yzzz_0[j] = pb_x * tg_xxyyy_yzzz_0[j] + wp_x[j] * tg_xxyyy_yzzz_1[j] + fl1_fx * tg_xyyy_yzzz_0[j] - fl1_fx * fl1_fza * tg_xyyy_yzzz_1[j];

                    tg_xxxyyy_zzzz_0[j] = pb_x * tg_xxyyy_zzzz_0[j] + wp_x[j] * tg_xxyyy_zzzz_1[j] + fl1_fx * tg_xyyy_zzzz_0[j] - fl1_fx * fl1_fza * tg_xyyy_zzzz_1[j];

                    tg_xxxyyz_xxxx_0[j] = pb_x * tg_xxyyz_xxxx_0[j] + wp_x[j] * tg_xxyyz_xxxx_1[j] + fl1_fx * tg_xyyz_xxxx_0[j] - fl1_fx * fl1_fza * tg_xyyz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxyyz_xxx_1[j];

                    tg_xxxyyz_xxxy_0[j] = pb_x * tg_xxyyz_xxxy_0[j] + wp_x[j] * tg_xxyyz_xxxy_1[j] + fl1_fx * tg_xyyz_xxxy_0[j] - fl1_fx * fl1_fza * tg_xyyz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxyyz_xxy_1[j];

                    tg_xxxyyz_xxxz_0[j] = pb_x * tg_xxyyz_xxxz_0[j] + wp_x[j] * tg_xxyyz_xxxz_1[j] + fl1_fx * tg_xyyz_xxxz_0[j] - fl1_fx * fl1_fza * tg_xyyz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxyyz_xxz_1[j];

                    tg_xxxyyz_xxyy_0[j] = pb_x * tg_xxyyz_xxyy_0[j] + wp_x[j] * tg_xxyyz_xxyy_1[j] + fl1_fx * tg_xyyz_xxyy_0[j] - fl1_fx * fl1_fza * tg_xyyz_xxyy_1[j] + fl1_fxn * tg_xxyyz_xyy_1[j];

                    tg_xxxyyz_xxyz_0[j] = pb_x * tg_xxyyz_xxyz_0[j] + wp_x[j] * tg_xxyyz_xxyz_1[j] + fl1_fx * tg_xyyz_xxyz_0[j] - fl1_fx * fl1_fza * tg_xyyz_xxyz_1[j] + fl1_fxn * tg_xxyyz_xyz_1[j];

                    tg_xxxyyz_xxzz_0[j] = pb_x * tg_xxyyz_xxzz_0[j] + wp_x[j] * tg_xxyyz_xxzz_1[j] + fl1_fx * tg_xyyz_xxzz_0[j] - fl1_fx * fl1_fza * tg_xyyz_xxzz_1[j] + fl1_fxn * tg_xxyyz_xzz_1[j];

                    tg_xxxyyz_xyyy_0[j] = pb_x * tg_xxyyz_xyyy_0[j] + wp_x[j] * tg_xxyyz_xyyy_1[j] + fl1_fx * tg_xyyz_xyyy_0[j] - fl1_fx * fl1_fza * tg_xyyz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxyyz_yyy_1[j];

                    tg_xxxyyz_xyyz_0[j] = pb_x * tg_xxyyz_xyyz_0[j] + wp_x[j] * tg_xxyyz_xyyz_1[j] + fl1_fx * tg_xyyz_xyyz_0[j] - fl1_fx * fl1_fza * tg_xyyz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxyyz_yyz_1[j];

                    tg_xxxyyz_xyzz_0[j] = pb_x * tg_xxyyz_xyzz_0[j] + wp_x[j] * tg_xxyyz_xyzz_1[j] + fl1_fx * tg_xyyz_xyzz_0[j] - fl1_fx * fl1_fza * tg_xyyz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxyyz_yzz_1[j];

                    tg_xxxyyz_xzzz_0[j] = pb_x * tg_xxyyz_xzzz_0[j] + wp_x[j] * tg_xxyyz_xzzz_1[j] + fl1_fx * tg_xyyz_xzzz_0[j] - fl1_fx * fl1_fza * tg_xyyz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxyyz_zzz_1[j];

                    tg_xxxyyz_yyyy_0[j] = pb_x * tg_xxyyz_yyyy_0[j] + wp_x[j] * tg_xxyyz_yyyy_1[j] + fl1_fx * tg_xyyz_yyyy_0[j] - fl1_fx * fl1_fza * tg_xyyz_yyyy_1[j];

                    tg_xxxyyz_yyyz_0[j] = pb_x * tg_xxyyz_yyyz_0[j] + wp_x[j] * tg_xxyyz_yyyz_1[j] + fl1_fx * tg_xyyz_yyyz_0[j] - fl1_fx * fl1_fza * tg_xyyz_yyyz_1[j];

                    tg_xxxyyz_yyzz_0[j] = pb_x * tg_xxyyz_yyzz_0[j] + wp_x[j] * tg_xxyyz_yyzz_1[j] + fl1_fx * tg_xyyz_yyzz_0[j] - fl1_fx * fl1_fza * tg_xyyz_yyzz_1[j];

                    tg_xxxyyz_yzzz_0[j] = pb_x * tg_xxyyz_yzzz_0[j] + wp_x[j] * tg_xxyyz_yzzz_1[j] + fl1_fx * tg_xyyz_yzzz_0[j] - fl1_fx * fl1_fza * tg_xyyz_yzzz_1[j];

                    tg_xxxyyz_zzzz_0[j] = pb_x * tg_xxyyz_zzzz_0[j] + wp_x[j] * tg_xxyyz_zzzz_1[j] + fl1_fx * tg_xyyz_zzzz_0[j] - fl1_fx * fl1_fza * tg_xyyz_zzzz_1[j];

                    tg_xxxyzz_xxxx_0[j] = pb_x * tg_xxyzz_xxxx_0[j] + wp_x[j] * tg_xxyzz_xxxx_1[j] + fl1_fx * tg_xyzz_xxxx_0[j] - fl1_fx * fl1_fza * tg_xyzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxyzz_xxx_1[j];

                    tg_xxxyzz_xxxy_0[j] = pb_x * tg_xxyzz_xxxy_0[j] + wp_x[j] * tg_xxyzz_xxxy_1[j] + fl1_fx * tg_xyzz_xxxy_0[j] - fl1_fx * fl1_fza * tg_xyzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxyzz_xxy_1[j];

                    tg_xxxyzz_xxxz_0[j] = pb_x * tg_xxyzz_xxxz_0[j] + wp_x[j] * tg_xxyzz_xxxz_1[j] + fl1_fx * tg_xyzz_xxxz_0[j] - fl1_fx * fl1_fza * tg_xyzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxyzz_xxz_1[j];

                    tg_xxxyzz_xxyy_0[j] = pb_x * tg_xxyzz_xxyy_0[j] + wp_x[j] * tg_xxyzz_xxyy_1[j] + fl1_fx * tg_xyzz_xxyy_0[j] - fl1_fx * fl1_fza * tg_xyzz_xxyy_1[j] + fl1_fxn * tg_xxyzz_xyy_1[j];

                    tg_xxxyzz_xxyz_0[j] = pb_x * tg_xxyzz_xxyz_0[j] + wp_x[j] * tg_xxyzz_xxyz_1[j] + fl1_fx * tg_xyzz_xxyz_0[j] - fl1_fx * fl1_fza * tg_xyzz_xxyz_1[j] + fl1_fxn * tg_xxyzz_xyz_1[j];

                    tg_xxxyzz_xxzz_0[j] = pb_x * tg_xxyzz_xxzz_0[j] + wp_x[j] * tg_xxyzz_xxzz_1[j] + fl1_fx * tg_xyzz_xxzz_0[j] - fl1_fx * fl1_fza * tg_xyzz_xxzz_1[j] + fl1_fxn * tg_xxyzz_xzz_1[j];

                    tg_xxxyzz_xyyy_0[j] = pb_x * tg_xxyzz_xyyy_0[j] + wp_x[j] * tg_xxyzz_xyyy_1[j] + fl1_fx * tg_xyzz_xyyy_0[j] - fl1_fx * fl1_fza * tg_xyzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxyzz_yyy_1[j];

                    tg_xxxyzz_xyyz_0[j] = pb_x * tg_xxyzz_xyyz_0[j] + wp_x[j] * tg_xxyzz_xyyz_1[j] + fl1_fx * tg_xyzz_xyyz_0[j] - fl1_fx * fl1_fza * tg_xyzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxyzz_yyz_1[j];

                    tg_xxxyzz_xyzz_0[j] = pb_x * tg_xxyzz_xyzz_0[j] + wp_x[j] * tg_xxyzz_xyzz_1[j] + fl1_fx * tg_xyzz_xyzz_0[j] - fl1_fx * fl1_fza * tg_xyzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxyzz_yzz_1[j];

                    tg_xxxyzz_xzzz_0[j] = pb_x * tg_xxyzz_xzzz_0[j] + wp_x[j] * tg_xxyzz_xzzz_1[j] + fl1_fx * tg_xyzz_xzzz_0[j] - fl1_fx * fl1_fza * tg_xyzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxyzz_zzz_1[j];

                    tg_xxxyzz_yyyy_0[j] = pb_x * tg_xxyzz_yyyy_0[j] + wp_x[j] * tg_xxyzz_yyyy_1[j] + fl1_fx * tg_xyzz_yyyy_0[j] - fl1_fx * fl1_fza * tg_xyzz_yyyy_1[j];

                    tg_xxxyzz_yyyz_0[j] = pb_x * tg_xxyzz_yyyz_0[j] + wp_x[j] * tg_xxyzz_yyyz_1[j] + fl1_fx * tg_xyzz_yyyz_0[j] - fl1_fx * fl1_fza * tg_xyzz_yyyz_1[j];

                    tg_xxxyzz_yyzz_0[j] = pb_x * tg_xxyzz_yyzz_0[j] + wp_x[j] * tg_xxyzz_yyzz_1[j] + fl1_fx * tg_xyzz_yyzz_0[j] - fl1_fx * fl1_fza * tg_xyzz_yyzz_1[j];

                    tg_xxxyzz_yzzz_0[j] = pb_x * tg_xxyzz_yzzz_0[j] + wp_x[j] * tg_xxyzz_yzzz_1[j] + fl1_fx * tg_xyzz_yzzz_0[j] - fl1_fx * fl1_fza * tg_xyzz_yzzz_1[j];

                    tg_xxxyzz_zzzz_0[j] = pb_x * tg_xxyzz_zzzz_0[j] + wp_x[j] * tg_xxyzz_zzzz_1[j] + fl1_fx * tg_xyzz_zzzz_0[j] - fl1_fx * fl1_fza * tg_xyzz_zzzz_1[j];

                    tg_xxxzzz_xxxx_0[j] = pb_x * tg_xxzzz_xxxx_0[j] + wp_x[j] * tg_xxzzz_xxxx_1[j] + fl1_fx * tg_xzzz_xxxx_0[j] - fl1_fx * fl1_fza * tg_xzzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxzzz_xxx_1[j];

                    tg_xxxzzz_xxxy_0[j] = pb_x * tg_xxzzz_xxxy_0[j] + wp_x[j] * tg_xxzzz_xxxy_1[j] + fl1_fx * tg_xzzz_xxxy_0[j] - fl1_fx * fl1_fza * tg_xzzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxzzz_xxy_1[j];

                    tg_xxxzzz_xxxz_0[j] = pb_x * tg_xxzzz_xxxz_0[j] + wp_x[j] * tg_xxzzz_xxxz_1[j] + fl1_fx * tg_xzzz_xxxz_0[j] - fl1_fx * fl1_fza * tg_xzzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxzzz_xxz_1[j];

                    tg_xxxzzz_xxyy_0[j] = pb_x * tg_xxzzz_xxyy_0[j] + wp_x[j] * tg_xxzzz_xxyy_1[j] + fl1_fx * tg_xzzz_xxyy_0[j] - fl1_fx * fl1_fza * tg_xzzz_xxyy_1[j] + fl1_fxn * tg_xxzzz_xyy_1[j];

                    tg_xxxzzz_xxyz_0[j] = pb_x * tg_xxzzz_xxyz_0[j] + wp_x[j] * tg_xxzzz_xxyz_1[j] + fl1_fx * tg_xzzz_xxyz_0[j] - fl1_fx * fl1_fza * tg_xzzz_xxyz_1[j] + fl1_fxn * tg_xxzzz_xyz_1[j];

                    tg_xxxzzz_xxzz_0[j] = pb_x * tg_xxzzz_xxzz_0[j] + wp_x[j] * tg_xxzzz_xxzz_1[j] + fl1_fx * tg_xzzz_xxzz_0[j] - fl1_fx * fl1_fza * tg_xzzz_xxzz_1[j] + fl1_fxn * tg_xxzzz_xzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSISG_141_188(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (141,188)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {6, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_6_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_6_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_5_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_5_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_4_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_5_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_xxzzz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 141); 

                auto tg_xxzzz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 142); 

                auto tg_xxzzz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 143); 

                auto tg_xxzzz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 144); 

                auto tg_xxzzz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 145); 

                auto tg_xxzzz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 146); 

                auto tg_xxzzz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 147); 

                auto tg_xxzzz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 148); 

                auto tg_xxzzz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 149); 

                auto tg_xyyyy_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 150); 

                auto tg_xyyyy_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 151); 

                auto tg_xyyyy_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 152); 

                auto tg_xyyyy_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 153); 

                auto tg_xyyyy_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 154); 

                auto tg_xyyyy_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 155); 

                auto tg_xyyyy_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 156); 

                auto tg_xyyyy_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 157); 

                auto tg_xyyyy_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 158); 

                auto tg_xyyyy_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 159); 

                auto tg_xyyyy_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 160); 

                auto tg_xyyyy_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 161); 

                auto tg_xyyyy_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 162); 

                auto tg_xyyyy_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 163); 

                auto tg_xyyyy_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 164); 

                auto tg_xyyyz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 165); 

                auto tg_xyyyz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 166); 

                auto tg_xyyyz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 167); 

                auto tg_xyyyz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 168); 

                auto tg_xyyyz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 169); 

                auto tg_xyyyz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 170); 

                auto tg_xyyyz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 171); 

                auto tg_xyyyz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 172); 

                auto tg_xyyyz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 173); 

                auto tg_xyyyz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 174); 

                auto tg_xyyyz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 175); 

                auto tg_xyyyz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 176); 

                auto tg_xyyyz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 177); 

                auto tg_xyyyz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 178); 

                auto tg_xyyyz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 179); 

                auto tg_xyyzz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 180); 

                auto tg_xyyzz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 181); 

                auto tg_xyyzz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 182); 

                auto tg_xyyzz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 183); 

                auto tg_xyyzz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 184); 

                auto tg_xyyzz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 185); 

                auto tg_xyyzz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 186); 

                auto tg_xyyzz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 187); 

                auto tg_xxzzz_xyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 141); 

                auto tg_xxzzz_xyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 142); 

                auto tg_xxzzz_xyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 143); 

                auto tg_xxzzz_xzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 144); 

                auto tg_xxzzz_yyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 145); 

                auto tg_xxzzz_yyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 146); 

                auto tg_xxzzz_yyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 147); 

                auto tg_xxzzz_yzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 148); 

                auto tg_xxzzz_zzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 149); 

                auto tg_xyyyy_xxxx_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 150); 

                auto tg_xyyyy_xxxy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 151); 

                auto tg_xyyyy_xxxz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 152); 

                auto tg_xyyyy_xxyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 153); 

                auto tg_xyyyy_xxyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 154); 

                auto tg_xyyyy_xxzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 155); 

                auto tg_xyyyy_xyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 156); 

                auto tg_xyyyy_xyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 157); 

                auto tg_xyyyy_xyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 158); 

                auto tg_xyyyy_xzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 159); 

                auto tg_xyyyy_yyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 160); 

                auto tg_xyyyy_yyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 161); 

                auto tg_xyyyy_yyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 162); 

                auto tg_xyyyy_yzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 163); 

                auto tg_xyyyy_zzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 164); 

                auto tg_xyyyz_xxxx_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 165); 

                auto tg_xyyyz_xxxy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 166); 

                auto tg_xyyyz_xxxz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 167); 

                auto tg_xyyyz_xxyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 168); 

                auto tg_xyyyz_xxyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 169); 

                auto tg_xyyyz_xxzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 170); 

                auto tg_xyyyz_xyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 171); 

                auto tg_xyyyz_xyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 172); 

                auto tg_xyyyz_xyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 173); 

                auto tg_xyyyz_xzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 174); 

                auto tg_xyyyz_yyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 175); 

                auto tg_xyyyz_yyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 176); 

                auto tg_xyyyz_yyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 177); 

                auto tg_xyyyz_yzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 178); 

                auto tg_xyyyz_zzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 179); 

                auto tg_xyyzz_xxxx_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 180); 

                auto tg_xyyzz_xxxy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 181); 

                auto tg_xyyzz_xxxz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 182); 

                auto tg_xyyzz_xxyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 183); 

                auto tg_xyyzz_xxyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 184); 

                auto tg_xyyzz_xxzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 185); 

                auto tg_xyyzz_xyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 186); 

                auto tg_xyyzz_xyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 187); 

                auto tg_xzzz_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 141); 

                auto tg_xzzz_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 142); 

                auto tg_xzzz_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 143); 

                auto tg_xzzz_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 144); 

                auto tg_xzzz_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 145); 

                auto tg_xzzz_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 146); 

                auto tg_xzzz_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 147); 

                auto tg_xzzz_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 148); 

                auto tg_xzzz_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 149); 

                auto tg_yyyy_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 150); 

                auto tg_yyyy_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 151); 

                auto tg_yyyy_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 152); 

                auto tg_yyyy_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 153); 

                auto tg_yyyy_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 154); 

                auto tg_yyyy_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 155); 

                auto tg_yyyy_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 156); 

                auto tg_yyyy_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 157); 

                auto tg_yyyy_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 158); 

                auto tg_yyyy_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 159); 

                auto tg_yyyy_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 160); 

                auto tg_yyyy_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 161); 

                auto tg_yyyy_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 162); 

                auto tg_yyyy_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 163); 

                auto tg_yyyy_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 164); 

                auto tg_yyyz_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 165); 

                auto tg_yyyz_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 166); 

                auto tg_yyyz_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 167); 

                auto tg_yyyz_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 168); 

                auto tg_yyyz_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 169); 

                auto tg_yyyz_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 170); 

                auto tg_yyyz_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 171); 

                auto tg_yyyz_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 172); 

                auto tg_yyyz_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 173); 

                auto tg_yyyz_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 174); 

                auto tg_yyyz_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 175); 

                auto tg_yyyz_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 176); 

                auto tg_yyyz_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 177); 

                auto tg_yyyz_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 178); 

                auto tg_yyyz_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 179); 

                auto tg_yyzz_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 180); 

                auto tg_yyzz_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 181); 

                auto tg_yyzz_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 182); 

                auto tg_yyzz_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 183); 

                auto tg_yyzz_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 184); 

                auto tg_yyzz_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 185); 

                auto tg_yyzz_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 186); 

                auto tg_yyzz_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 187); 

                auto tg_xzzz_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 141); 

                auto tg_xzzz_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 142); 

                auto tg_xzzz_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 143); 

                auto tg_xzzz_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 144); 

                auto tg_xzzz_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 145); 

                auto tg_xzzz_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 146); 

                auto tg_xzzz_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 147); 

                auto tg_xzzz_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 148); 

                auto tg_xzzz_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 149); 

                auto tg_yyyy_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 150); 

                auto tg_yyyy_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 151); 

                auto tg_yyyy_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 152); 

                auto tg_yyyy_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 153); 

                auto tg_yyyy_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 154); 

                auto tg_yyyy_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 155); 

                auto tg_yyyy_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 156); 

                auto tg_yyyy_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 157); 

                auto tg_yyyy_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 158); 

                auto tg_yyyy_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 159); 

                auto tg_yyyy_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 160); 

                auto tg_yyyy_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 161); 

                auto tg_yyyy_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 162); 

                auto tg_yyyy_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 163); 

                auto tg_yyyy_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 164); 

                auto tg_yyyz_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 165); 

                auto tg_yyyz_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 166); 

                auto tg_yyyz_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 167); 

                auto tg_yyyz_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 168); 

                auto tg_yyyz_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 169); 

                auto tg_yyyz_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 170); 

                auto tg_yyyz_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 171); 

                auto tg_yyyz_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 172); 

                auto tg_yyyz_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 173); 

                auto tg_yyyz_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 174); 

                auto tg_yyyz_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 175); 

                auto tg_yyyz_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 176); 

                auto tg_yyyz_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 177); 

                auto tg_yyyz_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 178); 

                auto tg_yyyz_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 179); 

                auto tg_yyzz_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 180); 

                auto tg_yyzz_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 181); 

                auto tg_yyzz_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 182); 

                auto tg_yyzz_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 183); 

                auto tg_yyzz_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 184); 

                auto tg_yyzz_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 185); 

                auto tg_yyzz_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 186); 

                auto tg_yyzz_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 187); 

                auto tg_xxzzz_yyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 96); 

                auto tg_xxzzz_yyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 97); 

                auto tg_xxzzz_yzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 98); 

                auto tg_xxzzz_zzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 99); 

                auto tg_xyyyy_xxx_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 100); 

                auto tg_xyyyy_xxy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 101); 

                auto tg_xyyyy_xxz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 102); 

                auto tg_xyyyy_xyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 103); 

                auto tg_xyyyy_xyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 104); 

                auto tg_xyyyy_xzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 105); 

                auto tg_xyyyy_yyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 106); 

                auto tg_xyyyy_yyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 107); 

                auto tg_xyyyy_yzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 108); 

                auto tg_xyyyy_zzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 109); 

                auto tg_xyyyz_xxx_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 110); 

                auto tg_xyyyz_xxy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 111); 

                auto tg_xyyyz_xxz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 112); 

                auto tg_xyyyz_xyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 113); 

                auto tg_xyyyz_xyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 114); 

                auto tg_xyyyz_xzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 115); 

                auto tg_xyyyz_yyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 116); 

                auto tg_xyyyz_yyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 117); 

                auto tg_xyyyz_yzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 118); 

                auto tg_xyyyz_zzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 119); 

                auto tg_xyyzz_xxx_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 120); 

                auto tg_xyyzz_xxy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 121); 

                auto tg_xyyzz_xxz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 122); 

                auto tg_xyyzz_xyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 123); 

                auto tg_xyyzz_xyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 124); 

                auto tg_xyyzz_xzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 125); 

                auto tg_xyyzz_yyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 126); 

                auto tg_xyyzz_yyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 127); 

                // set up pointers to integrals

                auto tg_xxxzzz_xyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 141); 

                auto tg_xxxzzz_xyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 142); 

                auto tg_xxxzzz_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 143); 

                auto tg_xxxzzz_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 144); 

                auto tg_xxxzzz_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 145); 

                auto tg_xxxzzz_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 146); 

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

                // Batch of Integrals (141,188)

                #pragma omp simd aligned(fxn, fza, tg_xxxzzz_xyyy_0, tg_xxxzzz_xyyz_0, tg_xxxzzz_xyzz_0, \
                                         tg_xxxzzz_xzzz_0, tg_xxxzzz_yyyy_0, tg_xxxzzz_yyyz_0, tg_xxxzzz_yyzz_0, \
                                         tg_xxxzzz_yzzz_0, tg_xxxzzz_zzzz_0, tg_xxyyyy_xxxx_0, tg_xxyyyy_xxxy_0, \
                                         tg_xxyyyy_xxxz_0, tg_xxyyyy_xxyy_0, tg_xxyyyy_xxyz_0, tg_xxyyyy_xxzz_0, \
                                         tg_xxyyyy_xyyy_0, tg_xxyyyy_xyyz_0, tg_xxyyyy_xyzz_0, tg_xxyyyy_xzzz_0, \
                                         tg_xxyyyy_yyyy_0, tg_xxyyyy_yyyz_0, tg_xxyyyy_yyzz_0, tg_xxyyyy_yzzz_0, \
                                         tg_xxyyyy_zzzz_0, tg_xxyyyz_xxxx_0, tg_xxyyyz_xxxy_0, tg_xxyyyz_xxxz_0, \
                                         tg_xxyyyz_xxyy_0, tg_xxyyyz_xxyz_0, tg_xxyyyz_xxzz_0, tg_xxyyyz_xyyy_0, \
                                         tg_xxyyyz_xyyz_0, tg_xxyyyz_xyzz_0, tg_xxyyyz_xzzz_0, tg_xxyyyz_yyyy_0, \
                                         tg_xxyyyz_yyyz_0, tg_xxyyyz_yyzz_0, tg_xxyyyz_yzzz_0, tg_xxyyyz_zzzz_0, \
                                         tg_xxyyzz_xxxx_0, tg_xxyyzz_xxxy_0, tg_xxyyzz_xxxz_0, tg_xxyyzz_xxyy_0, \
                                         tg_xxyyzz_xxyz_0, tg_xxyyzz_xxzz_0, tg_xxyyzz_xyyy_0, tg_xxyyzz_xyyz_0, \
                                         tg_xxzzz_xyyy_0, tg_xxzzz_xyyy_1, tg_xxzzz_xyyz_0, tg_xxzzz_xyyz_1, tg_xxzzz_xyzz_0, \
                                         tg_xxzzz_xyzz_1, tg_xxzzz_xzzz_0, tg_xxzzz_xzzz_1, tg_xxzzz_yyy_1, tg_xxzzz_yyyy_0, \
                                         tg_xxzzz_yyyy_1, tg_xxzzz_yyyz_0, tg_xxzzz_yyyz_1, tg_xxzzz_yyz_1, tg_xxzzz_yyzz_0, \
                                         tg_xxzzz_yyzz_1, tg_xxzzz_yzz_1, tg_xxzzz_yzzz_0, tg_xxzzz_yzzz_1, tg_xxzzz_zzz_1, \
                                         tg_xxzzz_zzzz_0, tg_xxzzz_zzzz_1, tg_xyyyy_xxx_1, tg_xyyyy_xxxx_0, tg_xyyyy_xxxx_1, \
                                         tg_xyyyy_xxxy_0, tg_xyyyy_xxxy_1, tg_xyyyy_xxxz_0, tg_xyyyy_xxxz_1, tg_xyyyy_xxy_1, \
                                         tg_xyyyy_xxyy_0, tg_xyyyy_xxyy_1, tg_xyyyy_xxyz_0, tg_xyyyy_xxyz_1, tg_xyyyy_xxz_1, \
                                         tg_xyyyy_xxzz_0, tg_xyyyy_xxzz_1, tg_xyyyy_xyy_1, tg_xyyyy_xyyy_0, tg_xyyyy_xyyy_1, \
                                         tg_xyyyy_xyyz_0, tg_xyyyy_xyyz_1, tg_xyyyy_xyz_1, tg_xyyyy_xyzz_0, tg_xyyyy_xyzz_1, \
                                         tg_xyyyy_xzz_1, tg_xyyyy_xzzz_0, tg_xyyyy_xzzz_1, tg_xyyyy_yyy_1, tg_xyyyy_yyyy_0, \
                                         tg_xyyyy_yyyy_1, tg_xyyyy_yyyz_0, tg_xyyyy_yyyz_1, tg_xyyyy_yyz_1, tg_xyyyy_yyzz_0, \
                                         tg_xyyyy_yyzz_1, tg_xyyyy_yzz_1, tg_xyyyy_yzzz_0, tg_xyyyy_yzzz_1, tg_xyyyy_zzz_1, \
                                         tg_xyyyy_zzzz_0, tg_xyyyy_zzzz_1, tg_xyyyz_xxx_1, tg_xyyyz_xxxx_0, tg_xyyyz_xxxx_1, \
                                         tg_xyyyz_xxxy_0, tg_xyyyz_xxxy_1, tg_xyyyz_xxxz_0, tg_xyyyz_xxxz_1, tg_xyyyz_xxy_1, \
                                         tg_xyyyz_xxyy_0, tg_xyyyz_xxyy_1, tg_xyyyz_xxyz_0, tg_xyyyz_xxyz_1, tg_xyyyz_xxz_1, \
                                         tg_xyyyz_xxzz_0, tg_xyyyz_xxzz_1, tg_xyyyz_xyy_1, tg_xyyyz_xyyy_0, tg_xyyyz_xyyy_1, \
                                         tg_xyyyz_xyyz_0, tg_xyyyz_xyyz_1, tg_xyyyz_xyz_1, tg_xyyyz_xyzz_0, tg_xyyyz_xyzz_1, \
                                         tg_xyyyz_xzz_1, tg_xyyyz_xzzz_0, tg_xyyyz_xzzz_1, tg_xyyyz_yyy_1, tg_xyyyz_yyyy_0, \
                                         tg_xyyyz_yyyy_1, tg_xyyyz_yyyz_0, tg_xyyyz_yyyz_1, tg_xyyyz_yyz_1, tg_xyyyz_yyzz_0, \
                                         tg_xyyyz_yyzz_1, tg_xyyyz_yzz_1, tg_xyyyz_yzzz_0, tg_xyyyz_yzzz_1, tg_xyyyz_zzz_1, \
                                         tg_xyyyz_zzzz_0, tg_xyyyz_zzzz_1, tg_xyyzz_xxx_1, tg_xyyzz_xxxx_0, tg_xyyzz_xxxx_1, \
                                         tg_xyyzz_xxxy_0, tg_xyyzz_xxxy_1, tg_xyyzz_xxxz_0, tg_xyyzz_xxxz_1, tg_xyyzz_xxy_1, \
                                         tg_xyyzz_xxyy_0, tg_xyyzz_xxyy_1, tg_xyyzz_xxyz_0, tg_xyyzz_xxyz_1, tg_xyyzz_xxz_1, \
                                         tg_xyyzz_xxzz_0, tg_xyyzz_xxzz_1, tg_xyyzz_xyy_1, tg_xyyzz_xyyy_0, tg_xyyzz_xyyy_1, \
                                         tg_xyyzz_xyyz_0, tg_xyyzz_xyyz_1, tg_xyyzz_xyz_1, tg_xyyzz_xzz_1, tg_xyyzz_yyy_1, \
                                         tg_xyyzz_yyz_1, tg_xzzz_xyyy_0, tg_xzzz_xyyy_1, tg_xzzz_xyyz_0, tg_xzzz_xyyz_1, \
                                         tg_xzzz_xyzz_0, tg_xzzz_xyzz_1, tg_xzzz_xzzz_0, tg_xzzz_xzzz_1, tg_xzzz_yyyy_0, \
                                         tg_xzzz_yyyy_1, tg_xzzz_yyyz_0, tg_xzzz_yyyz_1, tg_xzzz_yyzz_0, tg_xzzz_yyzz_1, \
                                         tg_xzzz_yzzz_0, tg_xzzz_yzzz_1, tg_xzzz_zzzz_0, tg_xzzz_zzzz_1, tg_yyyy_xxxx_0, \
                                         tg_yyyy_xxxx_1, tg_yyyy_xxxy_0, tg_yyyy_xxxy_1, tg_yyyy_xxxz_0, tg_yyyy_xxxz_1, \
                                         tg_yyyy_xxyy_0, tg_yyyy_xxyy_1, tg_yyyy_xxyz_0, tg_yyyy_xxyz_1, tg_yyyy_xxzz_0, \
                                         tg_yyyy_xxzz_1, tg_yyyy_xyyy_0, tg_yyyy_xyyy_1, tg_yyyy_xyyz_0, tg_yyyy_xyyz_1, \
                                         tg_yyyy_xyzz_0, tg_yyyy_xyzz_1, tg_yyyy_xzzz_0, tg_yyyy_xzzz_1, tg_yyyy_yyyy_0, \
                                         tg_yyyy_yyyy_1, tg_yyyy_yyyz_0, tg_yyyy_yyyz_1, tg_yyyy_yyzz_0, tg_yyyy_yyzz_1, \
                                         tg_yyyy_yzzz_0, tg_yyyy_yzzz_1, tg_yyyy_zzzz_0, tg_yyyy_zzzz_1, tg_yyyz_xxxx_0, \
                                         tg_yyyz_xxxx_1, tg_yyyz_xxxy_0, tg_yyyz_xxxy_1, tg_yyyz_xxxz_0, tg_yyyz_xxxz_1, \
                                         tg_yyyz_xxyy_0, tg_yyyz_xxyy_1, tg_yyyz_xxyz_0, tg_yyyz_xxyz_1, tg_yyyz_xxzz_0, \
                                         tg_yyyz_xxzz_1, tg_yyyz_xyyy_0, tg_yyyz_xyyy_1, tg_yyyz_xyyz_0, tg_yyyz_xyyz_1, \
                                         tg_yyyz_xyzz_0, tg_yyyz_xyzz_1, tg_yyyz_xzzz_0, tg_yyyz_xzzz_1, tg_yyyz_yyyy_0, \
                                         tg_yyyz_yyyy_1, tg_yyyz_yyyz_0, tg_yyyz_yyyz_1, tg_yyyz_yyzz_0, tg_yyyz_yyzz_1, \
                                         tg_yyyz_yzzz_0, tg_yyyz_yzzz_1, tg_yyyz_zzzz_0, tg_yyyz_zzzz_1, tg_yyzz_xxxx_0, \
                                         tg_yyzz_xxxx_1, tg_yyzz_xxxy_0, tg_yyzz_xxxy_1, tg_yyzz_xxxz_0, tg_yyzz_xxxz_1, \
                                         tg_yyzz_xxyy_0, tg_yyzz_xxyy_1, tg_yyzz_xxyz_0, tg_yyzz_xxyz_1, tg_yyzz_xxzz_0, \
                                         tg_yyzz_xxzz_1, tg_yyzz_xyyy_0, tg_yyzz_xyyy_1, tg_yyzz_xyyz_0, tg_yyzz_xyyz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxxzzz_xyyy_0[j] = pb_x * tg_xxzzz_xyyy_0[j] + wp_x[j] * tg_xxzzz_xyyy_1[j] + fl1_fx * tg_xzzz_xyyy_0[j] - fl1_fx * fl1_fza * tg_xzzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxzzz_yyy_1[j];

                    tg_xxxzzz_xyyz_0[j] = pb_x * tg_xxzzz_xyyz_0[j] + wp_x[j] * tg_xxzzz_xyyz_1[j] + fl1_fx * tg_xzzz_xyyz_0[j] - fl1_fx * fl1_fza * tg_xzzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxzzz_yyz_1[j];

                    tg_xxxzzz_xyzz_0[j] = pb_x * tg_xxzzz_xyzz_0[j] + wp_x[j] * tg_xxzzz_xyzz_1[j] + fl1_fx * tg_xzzz_xyzz_0[j] - fl1_fx * fl1_fza * tg_xzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxzzz_yzz_1[j];

                    tg_xxxzzz_xzzz_0[j] = pb_x * tg_xxzzz_xzzz_0[j] + wp_x[j] * tg_xxzzz_xzzz_1[j] + fl1_fx * tg_xzzz_xzzz_0[j] - fl1_fx * fl1_fza * tg_xzzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxzzz_zzz_1[j];

                    tg_xxxzzz_yyyy_0[j] = pb_x * tg_xxzzz_yyyy_0[j] + wp_x[j] * tg_xxzzz_yyyy_1[j] + fl1_fx * tg_xzzz_yyyy_0[j] - fl1_fx * fl1_fza * tg_xzzz_yyyy_1[j];

                    tg_xxxzzz_yyyz_0[j] = pb_x * tg_xxzzz_yyyz_0[j] + wp_x[j] * tg_xxzzz_yyyz_1[j] + fl1_fx * tg_xzzz_yyyz_0[j] - fl1_fx * fl1_fza * tg_xzzz_yyyz_1[j];

                    tg_xxxzzz_yyzz_0[j] = pb_x * tg_xxzzz_yyzz_0[j] + wp_x[j] * tg_xxzzz_yyzz_1[j] + fl1_fx * tg_xzzz_yyzz_0[j] - fl1_fx * fl1_fza * tg_xzzz_yyzz_1[j];

                    tg_xxxzzz_yzzz_0[j] = pb_x * tg_xxzzz_yzzz_0[j] + wp_x[j] * tg_xxzzz_yzzz_1[j] + fl1_fx * tg_xzzz_yzzz_0[j] - fl1_fx * fl1_fza * tg_xzzz_yzzz_1[j];

                    tg_xxxzzz_zzzz_0[j] = pb_x * tg_xxzzz_zzzz_0[j] + wp_x[j] * tg_xxzzz_zzzz_1[j] + fl1_fx * tg_xzzz_zzzz_0[j] - fl1_fx * fl1_fza * tg_xzzz_zzzz_1[j];

                    tg_xxyyyy_xxxx_0[j] = pb_x * tg_xyyyy_xxxx_0[j] + wp_x[j] * tg_xyyyy_xxxx_1[j] + 0.5 * fl1_fx * tg_yyyy_xxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyy_xxxx_1[j] + 2.0 * fl1_fxn * tg_xyyyy_xxx_1[j];

                    tg_xxyyyy_xxxy_0[j] = pb_x * tg_xyyyy_xxxy_0[j] + wp_x[j] * tg_xyyyy_xxxy_1[j] + 0.5 * fl1_fx * tg_yyyy_xxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyy_xxxy_1[j] + 1.5 * fl1_fxn * tg_xyyyy_xxy_1[j];

                    tg_xxyyyy_xxxz_0[j] = pb_x * tg_xyyyy_xxxz_0[j] + wp_x[j] * tg_xyyyy_xxxz_1[j] + 0.5 * fl1_fx * tg_yyyy_xxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyy_xxxz_1[j] + 1.5 * fl1_fxn * tg_xyyyy_xxz_1[j];

                    tg_xxyyyy_xxyy_0[j] = pb_x * tg_xyyyy_xxyy_0[j] + wp_x[j] * tg_xyyyy_xxyy_1[j] + 0.5 * fl1_fx * tg_yyyy_xxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyy_xxyy_1[j] + fl1_fxn * tg_xyyyy_xyy_1[j];

                    tg_xxyyyy_xxyz_0[j] = pb_x * tg_xyyyy_xxyz_0[j] + wp_x[j] * tg_xyyyy_xxyz_1[j] + 0.5 * fl1_fx * tg_yyyy_xxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyy_xxyz_1[j] + fl1_fxn * tg_xyyyy_xyz_1[j];

                    tg_xxyyyy_xxzz_0[j] = pb_x * tg_xyyyy_xxzz_0[j] + wp_x[j] * tg_xyyyy_xxzz_1[j] + 0.5 * fl1_fx * tg_yyyy_xxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyy_xxzz_1[j] + fl1_fxn * tg_xyyyy_xzz_1[j];

                    tg_xxyyyy_xyyy_0[j] = pb_x * tg_xyyyy_xyyy_0[j] + wp_x[j] * tg_xyyyy_xyyy_1[j] + 0.5 * fl1_fx * tg_yyyy_xyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyy_xyyy_1[j] + 0.5 * fl1_fxn * tg_xyyyy_yyy_1[j];

                    tg_xxyyyy_xyyz_0[j] = pb_x * tg_xyyyy_xyyz_0[j] + wp_x[j] * tg_xyyyy_xyyz_1[j] + 0.5 * fl1_fx * tg_yyyy_xyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyy_xyyz_1[j] + 0.5 * fl1_fxn * tg_xyyyy_yyz_1[j];

                    tg_xxyyyy_xyzz_0[j] = pb_x * tg_xyyyy_xyzz_0[j] + wp_x[j] * tg_xyyyy_xyzz_1[j] + 0.5 * fl1_fx * tg_yyyy_xyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyy_xyzz_1[j] + 0.5 * fl1_fxn * tg_xyyyy_yzz_1[j];

                    tg_xxyyyy_xzzz_0[j] = pb_x * tg_xyyyy_xzzz_0[j] + wp_x[j] * tg_xyyyy_xzzz_1[j] + 0.5 * fl1_fx * tg_yyyy_xzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyy_xzzz_1[j] + 0.5 * fl1_fxn * tg_xyyyy_zzz_1[j];

                    tg_xxyyyy_yyyy_0[j] = pb_x * tg_xyyyy_yyyy_0[j] + wp_x[j] * tg_xyyyy_yyyy_1[j] + 0.5 * fl1_fx * tg_yyyy_yyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyy_yyyy_1[j];

                    tg_xxyyyy_yyyz_0[j] = pb_x * tg_xyyyy_yyyz_0[j] + wp_x[j] * tg_xyyyy_yyyz_1[j] + 0.5 * fl1_fx * tg_yyyy_yyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyy_yyyz_1[j];

                    tg_xxyyyy_yyzz_0[j] = pb_x * tg_xyyyy_yyzz_0[j] + wp_x[j] * tg_xyyyy_yyzz_1[j] + 0.5 * fl1_fx * tg_yyyy_yyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyy_yyzz_1[j];

                    tg_xxyyyy_yzzz_0[j] = pb_x * tg_xyyyy_yzzz_0[j] + wp_x[j] * tg_xyyyy_yzzz_1[j] + 0.5 * fl1_fx * tg_yyyy_yzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyy_yzzz_1[j];

                    tg_xxyyyy_zzzz_0[j] = pb_x * tg_xyyyy_zzzz_0[j] + wp_x[j] * tg_xyyyy_zzzz_1[j] + 0.5 * fl1_fx * tg_yyyy_zzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyy_zzzz_1[j];

                    tg_xxyyyz_xxxx_0[j] = pb_x * tg_xyyyz_xxxx_0[j] + wp_x[j] * tg_xyyyz_xxxx_1[j] + 0.5 * fl1_fx * tg_yyyz_xxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xyyyz_xxx_1[j];

                    tg_xxyyyz_xxxy_0[j] = pb_x * tg_xyyyz_xxxy_0[j] + wp_x[j] * tg_xyyyz_xxxy_1[j] + 0.5 * fl1_fx * tg_yyyz_xxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xyyyz_xxy_1[j];

                    tg_xxyyyz_xxxz_0[j] = pb_x * tg_xyyyz_xxxz_0[j] + wp_x[j] * tg_xyyyz_xxxz_1[j] + 0.5 * fl1_fx * tg_yyyz_xxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xyyyz_xxz_1[j];

                    tg_xxyyyz_xxyy_0[j] = pb_x * tg_xyyyz_xxyy_0[j] + wp_x[j] * tg_xyyyz_xxyy_1[j] + 0.5 * fl1_fx * tg_yyyz_xxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyz_xxyy_1[j] + fl1_fxn * tg_xyyyz_xyy_1[j];

                    tg_xxyyyz_xxyz_0[j] = pb_x * tg_xyyyz_xxyz_0[j] + wp_x[j] * tg_xyyyz_xxyz_1[j] + 0.5 * fl1_fx * tg_yyyz_xxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyz_xxyz_1[j] + fl1_fxn * tg_xyyyz_xyz_1[j];

                    tg_xxyyyz_xxzz_0[j] = pb_x * tg_xyyyz_xxzz_0[j] + wp_x[j] * tg_xyyyz_xxzz_1[j] + 0.5 * fl1_fx * tg_yyyz_xxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyz_xxzz_1[j] + fl1_fxn * tg_xyyyz_xzz_1[j];

                    tg_xxyyyz_xyyy_0[j] = pb_x * tg_xyyyz_xyyy_0[j] + wp_x[j] * tg_xyyyz_xyyy_1[j] + 0.5 * fl1_fx * tg_yyyz_xyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xyyyz_yyy_1[j];

                    tg_xxyyyz_xyyz_0[j] = pb_x * tg_xyyyz_xyyz_0[j] + wp_x[j] * tg_xyyyz_xyyz_1[j] + 0.5 * fl1_fx * tg_yyyz_xyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xyyyz_yyz_1[j];

                    tg_xxyyyz_xyzz_0[j] = pb_x * tg_xyyyz_xyzz_0[j] + wp_x[j] * tg_xyyyz_xyzz_1[j] + 0.5 * fl1_fx * tg_yyyz_xyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xyyyz_yzz_1[j];

                    tg_xxyyyz_xzzz_0[j] = pb_x * tg_xyyyz_xzzz_0[j] + wp_x[j] * tg_xyyyz_xzzz_1[j] + 0.5 * fl1_fx * tg_yyyz_xzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xyyyz_zzz_1[j];

                    tg_xxyyyz_yyyy_0[j] = pb_x * tg_xyyyz_yyyy_0[j] + wp_x[j] * tg_xyyyz_yyyy_1[j] + 0.5 * fl1_fx * tg_yyyz_yyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyz_yyyy_1[j];

                    tg_xxyyyz_yyyz_0[j] = pb_x * tg_xyyyz_yyyz_0[j] + wp_x[j] * tg_xyyyz_yyyz_1[j] + 0.5 * fl1_fx * tg_yyyz_yyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyz_yyyz_1[j];

                    tg_xxyyyz_yyzz_0[j] = pb_x * tg_xyyyz_yyzz_0[j] + wp_x[j] * tg_xyyyz_yyzz_1[j] + 0.5 * fl1_fx * tg_yyyz_yyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyz_yyzz_1[j];

                    tg_xxyyyz_yzzz_0[j] = pb_x * tg_xyyyz_yzzz_0[j] + wp_x[j] * tg_xyyyz_yzzz_1[j] + 0.5 * fl1_fx * tg_yyyz_yzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyz_yzzz_1[j];

                    tg_xxyyyz_zzzz_0[j] = pb_x * tg_xyyyz_zzzz_0[j] + wp_x[j] * tg_xyyyz_zzzz_1[j] + 0.5 * fl1_fx * tg_yyyz_zzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyz_zzzz_1[j];

                    tg_xxyyzz_xxxx_0[j] = pb_x * tg_xyyzz_xxxx_0[j] + wp_x[j] * tg_xyyzz_xxxx_1[j] + 0.5 * fl1_fx * tg_yyzz_xxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xyyzz_xxx_1[j];

                    tg_xxyyzz_xxxy_0[j] = pb_x * tg_xyyzz_xxxy_0[j] + wp_x[j] * tg_xyyzz_xxxy_1[j] + 0.5 * fl1_fx * tg_yyzz_xxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xyyzz_xxy_1[j];

                    tg_xxyyzz_xxxz_0[j] = pb_x * tg_xyyzz_xxxz_0[j] + wp_x[j] * tg_xyyzz_xxxz_1[j] + 0.5 * fl1_fx * tg_yyzz_xxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xyyzz_xxz_1[j];

                    tg_xxyyzz_xxyy_0[j] = pb_x * tg_xyyzz_xxyy_0[j] + wp_x[j] * tg_xyyzz_xxyy_1[j] + 0.5 * fl1_fx * tg_yyzz_xxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzz_xxyy_1[j] + fl1_fxn * tg_xyyzz_xyy_1[j];

                    tg_xxyyzz_xxyz_0[j] = pb_x * tg_xyyzz_xxyz_0[j] + wp_x[j] * tg_xyyzz_xxyz_1[j] + 0.5 * fl1_fx * tg_yyzz_xxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzz_xxyz_1[j] + fl1_fxn * tg_xyyzz_xyz_1[j];

                    tg_xxyyzz_xxzz_0[j] = pb_x * tg_xyyzz_xxzz_0[j] + wp_x[j] * tg_xyyzz_xxzz_1[j] + 0.5 * fl1_fx * tg_yyzz_xxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzz_xxzz_1[j] + fl1_fxn * tg_xyyzz_xzz_1[j];

                    tg_xxyyzz_xyyy_0[j] = pb_x * tg_xyyzz_xyyy_0[j] + wp_x[j] * tg_xyyzz_xyyy_1[j] + 0.5 * fl1_fx * tg_yyzz_xyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xyyzz_yyy_1[j];

                    tg_xxyyzz_xyyz_0[j] = pb_x * tg_xyyzz_xyyz_0[j] + wp_x[j] * tg_xyyzz_xyyz_1[j] + 0.5 * fl1_fx * tg_yyzz_xyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xyyzz_yyz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSISG_188_235(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (188,235)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {6, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_6_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_6_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_5_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_5_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_4_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_5_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_xyyzz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 188); 

                auto tg_xyyzz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 189); 

                auto tg_xyyzz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 190); 

                auto tg_xyyzz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 191); 

                auto tg_xyyzz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 192); 

                auto tg_xyyzz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 193); 

                auto tg_xyyzz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 194); 

                auto tg_xyzzz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 195); 

                auto tg_xyzzz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 196); 

                auto tg_xyzzz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 197); 

                auto tg_xyzzz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 198); 

                auto tg_xyzzz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 199); 

                auto tg_xyzzz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 200); 

                auto tg_xyzzz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 201); 

                auto tg_xyzzz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 202); 

                auto tg_xyzzz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 203); 

                auto tg_xyzzz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 204); 

                auto tg_xyzzz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 205); 

                auto tg_xyzzz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 206); 

                auto tg_xyzzz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 207); 

                auto tg_xyzzz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 208); 

                auto tg_xyzzz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 209); 

                auto tg_xzzzz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 210); 

                auto tg_xzzzz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 211); 

                auto tg_xzzzz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 212); 

                auto tg_xzzzz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 213); 

                auto tg_xzzzz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 214); 

                auto tg_xzzzz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 215); 

                auto tg_xzzzz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 216); 

                auto tg_xzzzz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 217); 

                auto tg_xzzzz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 218); 

                auto tg_xzzzz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 219); 

                auto tg_xzzzz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 220); 

                auto tg_xzzzz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 221); 

                auto tg_xzzzz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 222); 

                auto tg_xzzzz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 223); 

                auto tg_xzzzz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 224); 

                auto tg_yyyyy_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 225); 

                auto tg_yyyyy_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 226); 

                auto tg_yyyyy_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 227); 

                auto tg_yyyyy_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 228); 

                auto tg_yyyyy_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 229); 

                auto tg_yyyyy_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 230); 

                auto tg_yyyyy_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 231); 

                auto tg_yyyyy_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 232); 

                auto tg_yyyyy_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 233); 

                auto tg_yyyyy_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 234); 

                auto tg_xyyzz_xyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 188); 

                auto tg_xyyzz_xzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 189); 

                auto tg_xyyzz_yyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 190); 

                auto tg_xyyzz_yyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 191); 

                auto tg_xyyzz_yyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 192); 

                auto tg_xyyzz_yzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 193); 

                auto tg_xyyzz_zzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 194); 

                auto tg_xyzzz_xxxx_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 195); 

                auto tg_xyzzz_xxxy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 196); 

                auto tg_xyzzz_xxxz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 197); 

                auto tg_xyzzz_xxyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 198); 

                auto tg_xyzzz_xxyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 199); 

                auto tg_xyzzz_xxzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 200); 

                auto tg_xyzzz_xyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 201); 

                auto tg_xyzzz_xyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 202); 

                auto tg_xyzzz_xyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 203); 

                auto tg_xyzzz_xzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 204); 

                auto tg_xyzzz_yyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 205); 

                auto tg_xyzzz_yyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 206); 

                auto tg_xyzzz_yyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 207); 

                auto tg_xyzzz_yzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 208); 

                auto tg_xyzzz_zzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 209); 

                auto tg_xzzzz_xxxx_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 210); 

                auto tg_xzzzz_xxxy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 211); 

                auto tg_xzzzz_xxxz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 212); 

                auto tg_xzzzz_xxyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 213); 

                auto tg_xzzzz_xxyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 214); 

                auto tg_xzzzz_xxzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 215); 

                auto tg_xzzzz_xyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 216); 

                auto tg_xzzzz_xyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 217); 

                auto tg_xzzzz_xyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 218); 

                auto tg_xzzzz_xzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 219); 

                auto tg_xzzzz_yyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 220); 

                auto tg_xzzzz_yyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 221); 

                auto tg_xzzzz_yyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 222); 

                auto tg_xzzzz_yzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 223); 

                auto tg_xzzzz_zzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 224); 

                auto tg_yyyyy_xxxx_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 225); 

                auto tg_yyyyy_xxxy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 226); 

                auto tg_yyyyy_xxxz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 227); 

                auto tg_yyyyy_xxyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 228); 

                auto tg_yyyyy_xxyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 229); 

                auto tg_yyyyy_xxzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 230); 

                auto tg_yyyyy_xyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 231); 

                auto tg_yyyyy_xyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 232); 

                auto tg_yyyyy_xyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 233); 

                auto tg_yyyyy_xzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 234); 

                auto tg_yyzz_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 188); 

                auto tg_yyzz_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 189); 

                auto tg_yyzz_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 190); 

                auto tg_yyzz_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 191); 

                auto tg_yyzz_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 192); 

                auto tg_yyzz_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 193); 

                auto tg_yyzz_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 194); 

                auto tg_yzzz_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 195); 

                auto tg_yzzz_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 196); 

                auto tg_yzzz_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 197); 

                auto tg_yzzz_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 198); 

                auto tg_yzzz_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 199); 

                auto tg_yzzz_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 200); 

                auto tg_yzzz_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 201); 

                auto tg_yzzz_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 202); 

                auto tg_yzzz_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 203); 

                auto tg_yzzz_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 204); 

                auto tg_yzzz_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 205); 

                auto tg_yzzz_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 206); 

                auto tg_yzzz_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 207); 

                auto tg_yzzz_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 208); 

                auto tg_yzzz_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 209); 

                auto tg_zzzz_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 210); 

                auto tg_zzzz_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 211); 

                auto tg_zzzz_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 212); 

                auto tg_zzzz_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 213); 

                auto tg_zzzz_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 214); 

                auto tg_zzzz_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 215); 

                auto tg_zzzz_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 216); 

                auto tg_zzzz_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 217); 

                auto tg_zzzz_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 218); 

                auto tg_zzzz_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 219); 

                auto tg_zzzz_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 220); 

                auto tg_zzzz_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 221); 

                auto tg_zzzz_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 222); 

                auto tg_zzzz_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 223); 

                auto tg_zzzz_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 224); 

                auto tg_yyzz_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 188); 

                auto tg_yyzz_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 189); 

                auto tg_yyzz_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 190); 

                auto tg_yyzz_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 191); 

                auto tg_yyzz_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 192); 

                auto tg_yyzz_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 193); 

                auto tg_yyzz_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 194); 

                auto tg_yzzz_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 195); 

                auto tg_yzzz_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 196); 

                auto tg_yzzz_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 197); 

                auto tg_yzzz_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 198); 

                auto tg_yzzz_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 199); 

                auto tg_yzzz_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 200); 

                auto tg_yzzz_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 201); 

                auto tg_yzzz_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 202); 

                auto tg_yzzz_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 203); 

                auto tg_yzzz_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 204); 

                auto tg_yzzz_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 205); 

                auto tg_yzzz_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 206); 

                auto tg_yzzz_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 207); 

                auto tg_yzzz_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 208); 

                auto tg_yzzz_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 209); 

                auto tg_zzzz_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 210); 

                auto tg_zzzz_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 211); 

                auto tg_zzzz_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 212); 

                auto tg_zzzz_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 213); 

                auto tg_zzzz_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 214); 

                auto tg_zzzz_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 215); 

                auto tg_zzzz_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 216); 

                auto tg_zzzz_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 217); 

                auto tg_zzzz_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 218); 

                auto tg_zzzz_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 219); 

                auto tg_zzzz_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 220); 

                auto tg_zzzz_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 221); 

                auto tg_zzzz_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 222); 

                auto tg_zzzz_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 223); 

                auto tg_zzzz_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 224); 

                auto tg_xyyzz_yzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 128); 

                auto tg_xyyzz_zzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 129); 

                auto tg_xyzzz_xxx_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 130); 

                auto tg_xyzzz_xxy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 131); 

                auto tg_xyzzz_xxz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 132); 

                auto tg_xyzzz_xyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 133); 

                auto tg_xyzzz_xyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 134); 

                auto tg_xyzzz_xzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 135); 

                auto tg_xyzzz_yyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 136); 

                auto tg_xyzzz_yyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 137); 

                auto tg_xyzzz_yzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 138); 

                auto tg_xyzzz_zzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 139); 

                auto tg_xzzzz_xxx_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 140); 

                auto tg_xzzzz_xxy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 141); 

                auto tg_xzzzz_xxz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 142); 

                auto tg_xzzzz_xyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 143); 

                auto tg_xzzzz_xyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 144); 

                auto tg_xzzzz_xzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 145); 

                auto tg_xzzzz_yyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 146); 

                auto tg_xzzzz_yyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 147); 

                auto tg_xzzzz_yzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 148); 

                auto tg_xzzzz_zzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 149); 

                auto tg_yyyyy_xxx_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 150); 

                auto tg_yyyyy_xxy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 151); 

                auto tg_yyyyy_xxz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 152); 

                auto tg_yyyyy_xyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 153); 

                auto tg_yyyyy_xyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 154); 

                auto tg_yyyyy_xzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 155); 

                auto tg_yyyyy_yyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 156); 

                auto tg_yyyyy_yyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 157); 

                auto tg_yyyyy_yzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 158); 

                auto tg_yyyyy_zzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 159); 

                // set up pointers to integrals

                auto tg_xxyyzz_xyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 188); 

                auto tg_xxyyzz_xzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 189); 

                auto tg_xxyyzz_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 190); 

                auto tg_xxyyzz_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 191); 

                auto tg_xxyyzz_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 192); 

                auto tg_xxyyzz_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 193); 

                auto tg_xxyyzz_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 194); 

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

                // Batch of Integrals (188,235)

                #pragma omp simd aligned(fxn, fza, tg_xxyyzz_xyzz_0, tg_xxyyzz_xzzz_0, tg_xxyyzz_yyyy_0, \
                                         tg_xxyyzz_yyyz_0, tg_xxyyzz_yyzz_0, tg_xxyyzz_yzzz_0, tg_xxyyzz_zzzz_0, \
                                         tg_xxyzzz_xxxx_0, tg_xxyzzz_xxxy_0, tg_xxyzzz_xxxz_0, tg_xxyzzz_xxyy_0, \
                                         tg_xxyzzz_xxyz_0, tg_xxyzzz_xxzz_0, tg_xxyzzz_xyyy_0, tg_xxyzzz_xyyz_0, \
                                         tg_xxyzzz_xyzz_0, tg_xxyzzz_xzzz_0, tg_xxyzzz_yyyy_0, tg_xxyzzz_yyyz_0, \
                                         tg_xxyzzz_yyzz_0, tg_xxyzzz_yzzz_0, tg_xxyzzz_zzzz_0, tg_xxzzzz_xxxx_0, \
                                         tg_xxzzzz_xxxy_0, tg_xxzzzz_xxxz_0, tg_xxzzzz_xxyy_0, tg_xxzzzz_xxyz_0, \
                                         tg_xxzzzz_xxzz_0, tg_xxzzzz_xyyy_0, tg_xxzzzz_xyyz_0, tg_xxzzzz_xyzz_0, \
                                         tg_xxzzzz_xzzz_0, tg_xxzzzz_yyyy_0, tg_xxzzzz_yyyz_0, tg_xxzzzz_yyzz_0, \
                                         tg_xxzzzz_yzzz_0, tg_xxzzzz_zzzz_0, tg_xyyyyy_xxxx_0, tg_xyyyyy_xxxy_0, \
                                         tg_xyyyyy_xxxz_0, tg_xyyyyy_xxyy_0, tg_xyyyyy_xxyz_0, tg_xyyyyy_xxzz_0, \
                                         tg_xyyyyy_xyyy_0, tg_xyyyyy_xyyz_0, tg_xyyyyy_xyzz_0, tg_xyyyyy_xzzz_0, \
                                         tg_xyyzz_xyzz_0, tg_xyyzz_xyzz_1, tg_xyyzz_xzzz_0, tg_xyyzz_xzzz_1, tg_xyyzz_yyyy_0, \
                                         tg_xyyzz_yyyy_1, tg_xyyzz_yyyz_0, tg_xyyzz_yyyz_1, tg_xyyzz_yyzz_0, tg_xyyzz_yyzz_1, \
                                         tg_xyyzz_yzz_1, tg_xyyzz_yzzz_0, tg_xyyzz_yzzz_1, tg_xyyzz_zzz_1, tg_xyyzz_zzzz_0, \
                                         tg_xyyzz_zzzz_1, tg_xyzzz_xxx_1, tg_xyzzz_xxxx_0, tg_xyzzz_xxxx_1, tg_xyzzz_xxxy_0, \
                                         tg_xyzzz_xxxy_1, tg_xyzzz_xxxz_0, tg_xyzzz_xxxz_1, tg_xyzzz_xxy_1, tg_xyzzz_xxyy_0, \
                                         tg_xyzzz_xxyy_1, tg_xyzzz_xxyz_0, tg_xyzzz_xxyz_1, tg_xyzzz_xxz_1, tg_xyzzz_xxzz_0, \
                                         tg_xyzzz_xxzz_1, tg_xyzzz_xyy_1, tg_xyzzz_xyyy_0, tg_xyzzz_xyyy_1, tg_xyzzz_xyyz_0, \
                                         tg_xyzzz_xyyz_1, tg_xyzzz_xyz_1, tg_xyzzz_xyzz_0, tg_xyzzz_xyzz_1, tg_xyzzz_xzz_1, \
                                         tg_xyzzz_xzzz_0, tg_xyzzz_xzzz_1, tg_xyzzz_yyy_1, tg_xyzzz_yyyy_0, tg_xyzzz_yyyy_1, \
                                         tg_xyzzz_yyyz_0, tg_xyzzz_yyyz_1, tg_xyzzz_yyz_1, tg_xyzzz_yyzz_0, tg_xyzzz_yyzz_1, \
                                         tg_xyzzz_yzz_1, tg_xyzzz_yzzz_0, tg_xyzzz_yzzz_1, tg_xyzzz_zzz_1, tg_xyzzz_zzzz_0, \
                                         tg_xyzzz_zzzz_1, tg_xzzzz_xxx_1, tg_xzzzz_xxxx_0, tg_xzzzz_xxxx_1, tg_xzzzz_xxxy_0, \
                                         tg_xzzzz_xxxy_1, tg_xzzzz_xxxz_0, tg_xzzzz_xxxz_1, tg_xzzzz_xxy_1, tg_xzzzz_xxyy_0, \
                                         tg_xzzzz_xxyy_1, tg_xzzzz_xxyz_0, tg_xzzzz_xxyz_1, tg_xzzzz_xxz_1, tg_xzzzz_xxzz_0, \
                                         tg_xzzzz_xxzz_1, tg_xzzzz_xyy_1, tg_xzzzz_xyyy_0, tg_xzzzz_xyyy_1, tg_xzzzz_xyyz_0, \
                                         tg_xzzzz_xyyz_1, tg_xzzzz_xyz_1, tg_xzzzz_xyzz_0, tg_xzzzz_xyzz_1, tg_xzzzz_xzz_1, \
                                         tg_xzzzz_xzzz_0, tg_xzzzz_xzzz_1, tg_xzzzz_yyy_1, tg_xzzzz_yyyy_0, tg_xzzzz_yyyy_1, \
                                         tg_xzzzz_yyyz_0, tg_xzzzz_yyyz_1, tg_xzzzz_yyz_1, tg_xzzzz_yyzz_0, tg_xzzzz_yyzz_1, \
                                         tg_xzzzz_yzz_1, tg_xzzzz_yzzz_0, tg_xzzzz_yzzz_1, tg_xzzzz_zzz_1, tg_xzzzz_zzzz_0, \
                                         tg_xzzzz_zzzz_1, tg_yyyyy_xxx_1, tg_yyyyy_xxxx_0, tg_yyyyy_xxxx_1, tg_yyyyy_xxxy_0, \
                                         tg_yyyyy_xxxy_1, tg_yyyyy_xxxz_0, tg_yyyyy_xxxz_1, tg_yyyyy_xxy_1, tg_yyyyy_xxyy_0, \
                                         tg_yyyyy_xxyy_1, tg_yyyyy_xxyz_0, tg_yyyyy_xxyz_1, tg_yyyyy_xxz_1, tg_yyyyy_xxzz_0, \
                                         tg_yyyyy_xxzz_1, tg_yyyyy_xyy_1, tg_yyyyy_xyyy_0, tg_yyyyy_xyyy_1, tg_yyyyy_xyyz_0, \
                                         tg_yyyyy_xyyz_1, tg_yyyyy_xyz_1, tg_yyyyy_xyzz_0, tg_yyyyy_xyzz_1, tg_yyyyy_xzz_1, \
                                         tg_yyyyy_xzzz_0, tg_yyyyy_xzzz_1, tg_yyyyy_yyy_1, tg_yyyyy_yyz_1, tg_yyyyy_yzz_1, \
                                         tg_yyyyy_zzz_1, tg_yyzz_xyzz_0, tg_yyzz_xyzz_1, tg_yyzz_xzzz_0, tg_yyzz_xzzz_1, \
                                         tg_yyzz_yyyy_0, tg_yyzz_yyyy_1, tg_yyzz_yyyz_0, tg_yyzz_yyyz_1, tg_yyzz_yyzz_0, \
                                         tg_yyzz_yyzz_1, tg_yyzz_yzzz_0, tg_yyzz_yzzz_1, tg_yyzz_zzzz_0, tg_yyzz_zzzz_1, \
                                         tg_yzzz_xxxx_0, tg_yzzz_xxxx_1, tg_yzzz_xxxy_0, tg_yzzz_xxxy_1, tg_yzzz_xxxz_0, \
                                         tg_yzzz_xxxz_1, tg_yzzz_xxyy_0, tg_yzzz_xxyy_1, tg_yzzz_xxyz_0, tg_yzzz_xxyz_1, \
                                         tg_yzzz_xxzz_0, tg_yzzz_xxzz_1, tg_yzzz_xyyy_0, tg_yzzz_xyyy_1, tg_yzzz_xyyz_0, \
                                         tg_yzzz_xyyz_1, tg_yzzz_xyzz_0, tg_yzzz_xyzz_1, tg_yzzz_xzzz_0, tg_yzzz_xzzz_1, \
                                         tg_yzzz_yyyy_0, tg_yzzz_yyyy_1, tg_yzzz_yyyz_0, tg_yzzz_yyyz_1, tg_yzzz_yyzz_0, \
                                         tg_yzzz_yyzz_1, tg_yzzz_yzzz_0, tg_yzzz_yzzz_1, tg_yzzz_zzzz_0, tg_yzzz_zzzz_1, \
                                         tg_zzzz_xxxx_0, tg_zzzz_xxxx_1, tg_zzzz_xxxy_0, tg_zzzz_xxxy_1, tg_zzzz_xxxz_0, \
                                         tg_zzzz_xxxz_1, tg_zzzz_xxyy_0, tg_zzzz_xxyy_1, tg_zzzz_xxyz_0, tg_zzzz_xxyz_1, \
                                         tg_zzzz_xxzz_0, tg_zzzz_xxzz_1, tg_zzzz_xyyy_0, tg_zzzz_xyyy_1, tg_zzzz_xyyz_0, \
                                         tg_zzzz_xyyz_1, tg_zzzz_xyzz_0, tg_zzzz_xyzz_1, tg_zzzz_xzzz_0, tg_zzzz_xzzz_1, \
                                         tg_zzzz_yyyy_0, tg_zzzz_yyyy_1, tg_zzzz_yyyz_0, tg_zzzz_yyyz_1, tg_zzzz_yyzz_0, \
                                         tg_zzzz_yyzz_1, tg_zzzz_yzzz_0, tg_zzzz_yzzz_1, tg_zzzz_zzzz_0, tg_zzzz_zzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxyyzz_xyzz_0[j] = pb_x * tg_xyyzz_xyzz_0[j] + wp_x[j] * tg_xyyzz_xyzz_1[j] + 0.5 * fl1_fx * tg_yyzz_xyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xyyzz_yzz_1[j];

                    tg_xxyyzz_xzzz_0[j] = pb_x * tg_xyyzz_xzzz_0[j] + wp_x[j] * tg_xyyzz_xzzz_1[j] + 0.5 * fl1_fx * tg_yyzz_xzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xyyzz_zzz_1[j];

                    tg_xxyyzz_yyyy_0[j] = pb_x * tg_xyyzz_yyyy_0[j] + wp_x[j] * tg_xyyzz_yyyy_1[j] + 0.5 * fl1_fx * tg_yyzz_yyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzz_yyyy_1[j];

                    tg_xxyyzz_yyyz_0[j] = pb_x * tg_xyyzz_yyyz_0[j] + wp_x[j] * tg_xyyzz_yyyz_1[j] + 0.5 * fl1_fx * tg_yyzz_yyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzz_yyyz_1[j];

                    tg_xxyyzz_yyzz_0[j] = pb_x * tg_xyyzz_yyzz_0[j] + wp_x[j] * tg_xyyzz_yyzz_1[j] + 0.5 * fl1_fx * tg_yyzz_yyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzz_yyzz_1[j];

                    tg_xxyyzz_yzzz_0[j] = pb_x * tg_xyyzz_yzzz_0[j] + wp_x[j] * tg_xyyzz_yzzz_1[j] + 0.5 * fl1_fx * tg_yyzz_yzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzz_yzzz_1[j];

                    tg_xxyyzz_zzzz_0[j] = pb_x * tg_xyyzz_zzzz_0[j] + wp_x[j] * tg_xyyzz_zzzz_1[j] + 0.5 * fl1_fx * tg_yyzz_zzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzz_zzzz_1[j];

                    tg_xxyzzz_xxxx_0[j] = pb_x * tg_xyzzz_xxxx_0[j] + wp_x[j] * tg_xyzzz_xxxx_1[j] + 0.5 * fl1_fx * tg_yzzz_xxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xyzzz_xxx_1[j];

                    tg_xxyzzz_xxxy_0[j] = pb_x * tg_xyzzz_xxxy_0[j] + wp_x[j] * tg_xyzzz_xxxy_1[j] + 0.5 * fl1_fx * tg_yzzz_xxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xyzzz_xxy_1[j];

                    tg_xxyzzz_xxxz_0[j] = pb_x * tg_xyzzz_xxxz_0[j] + wp_x[j] * tg_xyzzz_xxxz_1[j] + 0.5 * fl1_fx * tg_yzzz_xxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xyzzz_xxz_1[j];

                    tg_xxyzzz_xxyy_0[j] = pb_x * tg_xyzzz_xxyy_0[j] + wp_x[j] * tg_xyzzz_xxyy_1[j] + 0.5 * fl1_fx * tg_yzzz_xxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzz_xxyy_1[j] + fl1_fxn * tg_xyzzz_xyy_1[j];

                    tg_xxyzzz_xxyz_0[j] = pb_x * tg_xyzzz_xxyz_0[j] + wp_x[j] * tg_xyzzz_xxyz_1[j] + 0.5 * fl1_fx * tg_yzzz_xxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzz_xxyz_1[j] + fl1_fxn * tg_xyzzz_xyz_1[j];

                    tg_xxyzzz_xxzz_0[j] = pb_x * tg_xyzzz_xxzz_0[j] + wp_x[j] * tg_xyzzz_xxzz_1[j] + 0.5 * fl1_fx * tg_yzzz_xxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzz_xxzz_1[j] + fl1_fxn * tg_xyzzz_xzz_1[j];

                    tg_xxyzzz_xyyy_0[j] = pb_x * tg_xyzzz_xyyy_0[j] + wp_x[j] * tg_xyzzz_xyyy_1[j] + 0.5 * fl1_fx * tg_yzzz_xyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xyzzz_yyy_1[j];

                    tg_xxyzzz_xyyz_0[j] = pb_x * tg_xyzzz_xyyz_0[j] + wp_x[j] * tg_xyzzz_xyyz_1[j] + 0.5 * fl1_fx * tg_yzzz_xyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xyzzz_yyz_1[j];

                    tg_xxyzzz_xyzz_0[j] = pb_x * tg_xyzzz_xyzz_0[j] + wp_x[j] * tg_xyzzz_xyzz_1[j] + 0.5 * fl1_fx * tg_yzzz_xyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xyzzz_yzz_1[j];

                    tg_xxyzzz_xzzz_0[j] = pb_x * tg_xyzzz_xzzz_0[j] + wp_x[j] * tg_xyzzz_xzzz_1[j] + 0.5 * fl1_fx * tg_yzzz_xzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xyzzz_zzz_1[j];

                    tg_xxyzzz_yyyy_0[j] = pb_x * tg_xyzzz_yyyy_0[j] + wp_x[j] * tg_xyzzz_yyyy_1[j] + 0.5 * fl1_fx * tg_yzzz_yyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzz_yyyy_1[j];

                    tg_xxyzzz_yyyz_0[j] = pb_x * tg_xyzzz_yyyz_0[j] + wp_x[j] * tg_xyzzz_yyyz_1[j] + 0.5 * fl1_fx * tg_yzzz_yyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzz_yyyz_1[j];

                    tg_xxyzzz_yyzz_0[j] = pb_x * tg_xyzzz_yyzz_0[j] + wp_x[j] * tg_xyzzz_yyzz_1[j] + 0.5 * fl1_fx * tg_yzzz_yyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzz_yyzz_1[j];

                    tg_xxyzzz_yzzz_0[j] = pb_x * tg_xyzzz_yzzz_0[j] + wp_x[j] * tg_xyzzz_yzzz_1[j] + 0.5 * fl1_fx * tg_yzzz_yzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzz_yzzz_1[j];

                    tg_xxyzzz_zzzz_0[j] = pb_x * tg_xyzzz_zzzz_0[j] + wp_x[j] * tg_xyzzz_zzzz_1[j] + 0.5 * fl1_fx * tg_yzzz_zzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzz_zzzz_1[j];

                    tg_xxzzzz_xxxx_0[j] = pb_x * tg_xzzzz_xxxx_0[j] + wp_x[j] * tg_xzzzz_xxxx_1[j] + 0.5 * fl1_fx * tg_zzzz_xxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xzzzz_xxx_1[j];

                    tg_xxzzzz_xxxy_0[j] = pb_x * tg_xzzzz_xxxy_0[j] + wp_x[j] * tg_xzzzz_xxxy_1[j] + 0.5 * fl1_fx * tg_zzzz_xxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xzzzz_xxy_1[j];

                    tg_xxzzzz_xxxz_0[j] = pb_x * tg_xzzzz_xxxz_0[j] + wp_x[j] * tg_xzzzz_xxxz_1[j] + 0.5 * fl1_fx * tg_zzzz_xxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xzzzz_xxz_1[j];

                    tg_xxzzzz_xxyy_0[j] = pb_x * tg_xzzzz_xxyy_0[j] + wp_x[j] * tg_xzzzz_xxyy_1[j] + 0.5 * fl1_fx * tg_zzzz_xxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_xxyy_1[j] + fl1_fxn * tg_xzzzz_xyy_1[j];

                    tg_xxzzzz_xxyz_0[j] = pb_x * tg_xzzzz_xxyz_0[j] + wp_x[j] * tg_xzzzz_xxyz_1[j] + 0.5 * fl1_fx * tg_zzzz_xxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_xxyz_1[j] + fl1_fxn * tg_xzzzz_xyz_1[j];

                    tg_xxzzzz_xxzz_0[j] = pb_x * tg_xzzzz_xxzz_0[j] + wp_x[j] * tg_xzzzz_xxzz_1[j] + 0.5 * fl1_fx * tg_zzzz_xxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_xxzz_1[j] + fl1_fxn * tg_xzzzz_xzz_1[j];

                    tg_xxzzzz_xyyy_0[j] = pb_x * tg_xzzzz_xyyy_0[j] + wp_x[j] * tg_xzzzz_xyyy_1[j] + 0.5 * fl1_fx * tg_zzzz_xyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xzzzz_yyy_1[j];

                    tg_xxzzzz_xyyz_0[j] = pb_x * tg_xzzzz_xyyz_0[j] + wp_x[j] * tg_xzzzz_xyyz_1[j] + 0.5 * fl1_fx * tg_zzzz_xyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xzzzz_yyz_1[j];

                    tg_xxzzzz_xyzz_0[j] = pb_x * tg_xzzzz_xyzz_0[j] + wp_x[j] * tg_xzzzz_xyzz_1[j] + 0.5 * fl1_fx * tg_zzzz_xyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xzzzz_yzz_1[j];

                    tg_xxzzzz_xzzz_0[j] = pb_x * tg_xzzzz_xzzz_0[j] + wp_x[j] * tg_xzzzz_xzzz_1[j] + 0.5 * fl1_fx * tg_zzzz_xzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xzzzz_zzz_1[j];

                    tg_xxzzzz_yyyy_0[j] = pb_x * tg_xzzzz_yyyy_0[j] + wp_x[j] * tg_xzzzz_yyyy_1[j] + 0.5 * fl1_fx * tg_zzzz_yyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_yyyy_1[j];

                    tg_xxzzzz_yyyz_0[j] = pb_x * tg_xzzzz_yyyz_0[j] + wp_x[j] * tg_xzzzz_yyyz_1[j] + 0.5 * fl1_fx * tg_zzzz_yyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_yyyz_1[j];

                    tg_xxzzzz_yyzz_0[j] = pb_x * tg_xzzzz_yyzz_0[j] + wp_x[j] * tg_xzzzz_yyzz_1[j] + 0.5 * fl1_fx * tg_zzzz_yyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_yyzz_1[j];

                    tg_xxzzzz_yzzz_0[j] = pb_x * tg_xzzzz_yzzz_0[j] + wp_x[j] * tg_xzzzz_yzzz_1[j] + 0.5 * fl1_fx * tg_zzzz_yzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_yzzz_1[j];

                    tg_xxzzzz_zzzz_0[j] = pb_x * tg_xzzzz_zzzz_0[j] + wp_x[j] * tg_xzzzz_zzzz_1[j] + 0.5 * fl1_fx * tg_zzzz_zzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_zzzz_1[j];

                    tg_xyyyyy_xxxx_0[j] = pb_x * tg_yyyyy_xxxx_0[j] + wp_x[j] * tg_yyyyy_xxxx_1[j] + 2.0 * fl1_fxn * tg_yyyyy_xxx_1[j];

                    tg_xyyyyy_xxxy_0[j] = pb_x * tg_yyyyy_xxxy_0[j] + wp_x[j] * tg_yyyyy_xxxy_1[j] + 1.5 * fl1_fxn * tg_yyyyy_xxy_1[j];

                    tg_xyyyyy_xxxz_0[j] = pb_x * tg_yyyyy_xxxz_0[j] + wp_x[j] * tg_yyyyy_xxxz_1[j] + 1.5 * fl1_fxn * tg_yyyyy_xxz_1[j];

                    tg_xyyyyy_xxyy_0[j] = pb_x * tg_yyyyy_xxyy_0[j] + wp_x[j] * tg_yyyyy_xxyy_1[j] + fl1_fxn * tg_yyyyy_xyy_1[j];

                    tg_xyyyyy_xxyz_0[j] = pb_x * tg_yyyyy_xxyz_0[j] + wp_x[j] * tg_yyyyy_xxyz_1[j] + fl1_fxn * tg_yyyyy_xyz_1[j];

                    tg_xyyyyy_xxzz_0[j] = pb_x * tg_yyyyy_xxzz_0[j] + wp_x[j] * tg_yyyyy_xxzz_1[j] + fl1_fxn * tg_yyyyy_xzz_1[j];

                    tg_xyyyyy_xyyy_0[j] = pb_x * tg_yyyyy_xyyy_0[j] + wp_x[j] * tg_yyyyy_xyyy_1[j] + 0.5 * fl1_fxn * tg_yyyyy_yyy_1[j];

                    tg_xyyyyy_xyyz_0[j] = pb_x * tg_yyyyy_xyyz_0[j] + wp_x[j] * tg_yyyyy_xyyz_1[j] + 0.5 * fl1_fxn * tg_yyyyy_yyz_1[j];

                    tg_xyyyyy_xyzz_0[j] = pb_x * tg_yyyyy_xyzz_0[j] + wp_x[j] * tg_yyyyy_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyyyy_yzz_1[j];

                    tg_xyyyyy_xzzz_0[j] = pb_x * tg_yyyyy_xzzz_0[j] + wp_x[j] * tg_yyyyy_xzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyy_zzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSISG_235_282(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (235,282)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {6, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_6_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_6_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_5_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_5_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_5_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_yyyyy_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 235); 

                auto tg_yyyyy_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 236); 

                auto tg_yyyyy_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 237); 

                auto tg_yyyyy_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 238); 

                auto tg_yyyyy_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 239); 

                auto tg_yyyyz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 240); 

                auto tg_yyyyz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 241); 

                auto tg_yyyyz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 242); 

                auto tg_yyyyz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 243); 

                auto tg_yyyyz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 244); 

                auto tg_yyyyz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 245); 

                auto tg_yyyyz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 246); 

                auto tg_yyyyz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 247); 

                auto tg_yyyyz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 248); 

                auto tg_yyyyz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 249); 

                auto tg_yyyyz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 250); 

                auto tg_yyyyz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 251); 

                auto tg_yyyyz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 252); 

                auto tg_yyyyz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 253); 

                auto tg_yyyyz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 254); 

                auto tg_yyyzz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 255); 

                auto tg_yyyzz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 256); 

                auto tg_yyyzz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 257); 

                auto tg_yyyzz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 258); 

                auto tg_yyyzz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 259); 

                auto tg_yyyzz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 260); 

                auto tg_yyyzz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 261); 

                auto tg_yyyzz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 262); 

                auto tg_yyyzz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 263); 

                auto tg_yyyzz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 264); 

                auto tg_yyyzz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 265); 

                auto tg_yyyzz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 266); 

                auto tg_yyyzz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 267); 

                auto tg_yyyzz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 268); 

                auto tg_yyyzz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 269); 

                auto tg_yyzzz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 270); 

                auto tg_yyzzz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 271); 

                auto tg_yyzzz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 272); 

                auto tg_yyzzz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 273); 

                auto tg_yyzzz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 274); 

                auto tg_yyzzz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 275); 

                auto tg_yyzzz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 276); 

                auto tg_yyzzz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 277); 

                auto tg_yyzzz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 278); 

                auto tg_yyzzz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 279); 

                auto tg_yyzzz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 280); 

                auto tg_yyzzz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 281); 

                auto tg_yyyyy_yyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 235); 

                auto tg_yyyyy_yyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 236); 

                auto tg_yyyyy_yyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 237); 

                auto tg_yyyyy_yzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 238); 

                auto tg_yyyyy_zzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 239); 

                auto tg_yyyyz_xxxx_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 240); 

                auto tg_yyyyz_xxxy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 241); 

                auto tg_yyyyz_xxxz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 242); 

                auto tg_yyyyz_xxyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 243); 

                auto tg_yyyyz_xxyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 244); 

                auto tg_yyyyz_xxzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 245); 

                auto tg_yyyyz_xyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 246); 

                auto tg_yyyyz_xyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 247); 

                auto tg_yyyyz_xyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 248); 

                auto tg_yyyyz_xzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 249); 

                auto tg_yyyyz_yyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 250); 

                auto tg_yyyyz_yyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 251); 

                auto tg_yyyyz_yyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 252); 

                auto tg_yyyyz_yzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 253); 

                auto tg_yyyyz_zzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 254); 

                auto tg_yyyzz_xxxx_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 255); 

                auto tg_yyyzz_xxxy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 256); 

                auto tg_yyyzz_xxxz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 257); 

                auto tg_yyyzz_xxyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 258); 

                auto tg_yyyzz_xxyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 259); 

                auto tg_yyyzz_xxzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 260); 

                auto tg_yyyzz_xyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 261); 

                auto tg_yyyzz_xyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 262); 

                auto tg_yyyzz_xyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 263); 

                auto tg_yyyzz_xzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 264); 

                auto tg_yyyzz_yyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 265); 

                auto tg_yyyzz_yyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 266); 

                auto tg_yyyzz_yyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 267); 

                auto tg_yyyzz_yzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 268); 

                auto tg_yyyzz_zzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 269); 

                auto tg_yyzzz_xxxx_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 270); 

                auto tg_yyzzz_xxxy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 271); 

                auto tg_yyzzz_xxxz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 272); 

                auto tg_yyzzz_xxyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 273); 

                auto tg_yyzzz_xxyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 274); 

                auto tg_yyzzz_xxzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 275); 

                auto tg_yyzzz_xyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 276); 

                auto tg_yyzzz_xyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 277); 

                auto tg_yyzzz_xyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 278); 

                auto tg_yyzzz_xzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 279); 

                auto tg_yyzzz_yyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 280); 

                auto tg_yyzzz_yyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 281); 

                auto tg_yyyyz_xxx_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 160); 

                auto tg_yyyyz_xxy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 161); 

                auto tg_yyyyz_xxz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 162); 

                auto tg_yyyyz_xyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 163); 

                auto tg_yyyyz_xyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 164); 

                auto tg_yyyyz_xzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 165); 

                auto tg_yyyyz_yyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 166); 

                auto tg_yyyyz_yyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 167); 

                auto tg_yyyyz_yzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 168); 

                auto tg_yyyyz_zzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 169); 

                auto tg_yyyzz_xxx_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 170); 

                auto tg_yyyzz_xxy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 171); 

                auto tg_yyyzz_xxz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 172); 

                auto tg_yyyzz_xyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 173); 

                auto tg_yyyzz_xyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 174); 

                auto tg_yyyzz_xzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 175); 

                auto tg_yyyzz_yyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 176); 

                auto tg_yyyzz_yyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 177); 

                auto tg_yyyzz_yzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 178); 

                auto tg_yyyzz_zzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 179); 

                auto tg_yyzzz_xxx_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 180); 

                auto tg_yyzzz_xxy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 181); 

                auto tg_yyzzz_xxz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 182); 

                auto tg_yyzzz_xyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 183); 

                auto tg_yyzzz_xyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 184); 

                auto tg_yyzzz_xzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 185); 

                auto tg_yyzzz_yyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 186); 

                auto tg_yyzzz_yyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 187); 

                auto tg_yyzzz_yzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 188); 

                auto tg_yyzzz_zzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 189); 

                // set up pointers to integrals

                auto tg_xyyyyy_yyyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 235); 

                auto tg_xyyyyy_yyyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 236); 

                auto tg_xyyyyy_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 237); 

                auto tg_xyyyyy_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 238); 

                auto tg_xyyyyy_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 239); 

                auto tg_xyyyyz_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 240); 

                auto tg_xyyyyz_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 241); 

                auto tg_xyyyyz_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 242); 

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

                // Batch of Integrals (235,282)

                #pragma omp simd aligned(fxn, tg_xyyyyy_yyyy_0, tg_xyyyyy_yyyz_0, tg_xyyyyy_yyzz_0, \
                                         tg_xyyyyy_yzzz_0, tg_xyyyyy_zzzz_0, tg_xyyyyz_xxxx_0, tg_xyyyyz_xxxy_0, \
                                         tg_xyyyyz_xxxz_0, tg_xyyyyz_xxyy_0, tg_xyyyyz_xxyz_0, tg_xyyyyz_xxzz_0, \
                                         tg_xyyyyz_xyyy_0, tg_xyyyyz_xyyz_0, tg_xyyyyz_xyzz_0, tg_xyyyyz_xzzz_0, \
                                         tg_xyyyyz_yyyy_0, tg_xyyyyz_yyyz_0, tg_xyyyyz_yyzz_0, tg_xyyyyz_yzzz_0, \
                                         tg_xyyyyz_zzzz_0, tg_xyyyzz_xxxx_0, tg_xyyyzz_xxxy_0, tg_xyyyzz_xxxz_0, \
                                         tg_xyyyzz_xxyy_0, tg_xyyyzz_xxyz_0, tg_xyyyzz_xxzz_0, tg_xyyyzz_xyyy_0, \
                                         tg_xyyyzz_xyyz_0, tg_xyyyzz_xyzz_0, tg_xyyyzz_xzzz_0, tg_xyyyzz_yyyy_0, \
                                         tg_xyyyzz_yyyz_0, tg_xyyyzz_yyzz_0, tg_xyyyzz_yzzz_0, tg_xyyyzz_zzzz_0, \
                                         tg_xyyzzz_xxxx_0, tg_xyyzzz_xxxy_0, tg_xyyzzz_xxxz_0, tg_xyyzzz_xxyy_0, \
                                         tg_xyyzzz_xxyz_0, tg_xyyzzz_xxzz_0, tg_xyyzzz_xyyy_0, tg_xyyzzz_xyyz_0, \
                                         tg_xyyzzz_xyzz_0, tg_xyyzzz_xzzz_0, tg_xyyzzz_yyyy_0, tg_xyyzzz_yyyz_0, \
                                         tg_yyyyy_yyyy_0, tg_yyyyy_yyyy_1, tg_yyyyy_yyyz_0, tg_yyyyy_yyyz_1, tg_yyyyy_yyzz_0, \
                                         tg_yyyyy_yyzz_1, tg_yyyyy_yzzz_0, tg_yyyyy_yzzz_1, tg_yyyyy_zzzz_0, tg_yyyyy_zzzz_1, \
                                         tg_yyyyz_xxx_1, tg_yyyyz_xxxx_0, tg_yyyyz_xxxx_1, tg_yyyyz_xxxy_0, tg_yyyyz_xxxy_1, \
                                         tg_yyyyz_xxxz_0, tg_yyyyz_xxxz_1, tg_yyyyz_xxy_1, tg_yyyyz_xxyy_0, tg_yyyyz_xxyy_1, \
                                         tg_yyyyz_xxyz_0, tg_yyyyz_xxyz_1, tg_yyyyz_xxz_1, tg_yyyyz_xxzz_0, tg_yyyyz_xxzz_1, \
                                         tg_yyyyz_xyy_1, tg_yyyyz_xyyy_0, tg_yyyyz_xyyy_1, tg_yyyyz_xyyz_0, tg_yyyyz_xyyz_1, \
                                         tg_yyyyz_xyz_1, tg_yyyyz_xyzz_0, tg_yyyyz_xyzz_1, tg_yyyyz_xzz_1, tg_yyyyz_xzzz_0, \
                                         tg_yyyyz_xzzz_1, tg_yyyyz_yyy_1, tg_yyyyz_yyyy_0, tg_yyyyz_yyyy_1, tg_yyyyz_yyyz_0, \
                                         tg_yyyyz_yyyz_1, tg_yyyyz_yyz_1, tg_yyyyz_yyzz_0, tg_yyyyz_yyzz_1, tg_yyyyz_yzz_1, \
                                         tg_yyyyz_yzzz_0, tg_yyyyz_yzzz_1, tg_yyyyz_zzz_1, tg_yyyyz_zzzz_0, tg_yyyyz_zzzz_1, \
                                         tg_yyyzz_xxx_1, tg_yyyzz_xxxx_0, tg_yyyzz_xxxx_1, tg_yyyzz_xxxy_0, tg_yyyzz_xxxy_1, \
                                         tg_yyyzz_xxxz_0, tg_yyyzz_xxxz_1, tg_yyyzz_xxy_1, tg_yyyzz_xxyy_0, tg_yyyzz_xxyy_1, \
                                         tg_yyyzz_xxyz_0, tg_yyyzz_xxyz_1, tg_yyyzz_xxz_1, tg_yyyzz_xxzz_0, tg_yyyzz_xxzz_1, \
                                         tg_yyyzz_xyy_1, tg_yyyzz_xyyy_0, tg_yyyzz_xyyy_1, tg_yyyzz_xyyz_0, tg_yyyzz_xyyz_1, \
                                         tg_yyyzz_xyz_1, tg_yyyzz_xyzz_0, tg_yyyzz_xyzz_1, tg_yyyzz_xzz_1, tg_yyyzz_xzzz_0, \
                                         tg_yyyzz_xzzz_1, tg_yyyzz_yyy_1, tg_yyyzz_yyyy_0, tg_yyyzz_yyyy_1, tg_yyyzz_yyyz_0, \
                                         tg_yyyzz_yyyz_1, tg_yyyzz_yyz_1, tg_yyyzz_yyzz_0, tg_yyyzz_yyzz_1, tg_yyyzz_yzz_1, \
                                         tg_yyyzz_yzzz_0, tg_yyyzz_yzzz_1, tg_yyyzz_zzz_1, tg_yyyzz_zzzz_0, tg_yyyzz_zzzz_1, \
                                         tg_yyzzz_xxx_1, tg_yyzzz_xxxx_0, tg_yyzzz_xxxx_1, tg_yyzzz_xxxy_0, tg_yyzzz_xxxy_1, \
                                         tg_yyzzz_xxxz_0, tg_yyzzz_xxxz_1, tg_yyzzz_xxy_1, tg_yyzzz_xxyy_0, tg_yyzzz_xxyy_1, \
                                         tg_yyzzz_xxyz_0, tg_yyzzz_xxyz_1, tg_yyzzz_xxz_1, tg_yyzzz_xxzz_0, tg_yyzzz_xxzz_1, \
                                         tg_yyzzz_xyy_1, tg_yyzzz_xyyy_0, tg_yyzzz_xyyy_1, tg_yyzzz_xyyz_0, tg_yyzzz_xyyz_1, \
                                         tg_yyzzz_xyz_1, tg_yyzzz_xyzz_0, tg_yyzzz_xyzz_1, tg_yyzzz_xzz_1, tg_yyzzz_xzzz_0, \
                                         tg_yyzzz_xzzz_1, tg_yyzzz_yyy_1, tg_yyzzz_yyyy_0, tg_yyzzz_yyyy_1, tg_yyzzz_yyyz_0, \
                                         tg_yyzzz_yyyz_1, tg_yyzzz_yyz_1, tg_yyzzz_yzz_1, tg_yyzzz_zzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fxn = fxn[j];

                    tg_xyyyyy_yyyy_0[j] = pb_x * tg_yyyyy_yyyy_0[j] + wp_x[j] * tg_yyyyy_yyyy_1[j];

                    tg_xyyyyy_yyyz_0[j] = pb_x * tg_yyyyy_yyyz_0[j] + wp_x[j] * tg_yyyyy_yyyz_1[j];

                    tg_xyyyyy_yyzz_0[j] = pb_x * tg_yyyyy_yyzz_0[j] + wp_x[j] * tg_yyyyy_yyzz_1[j];

                    tg_xyyyyy_yzzz_0[j] = pb_x * tg_yyyyy_yzzz_0[j] + wp_x[j] * tg_yyyyy_yzzz_1[j];

                    tg_xyyyyy_zzzz_0[j] = pb_x * tg_yyyyy_zzzz_0[j] + wp_x[j] * tg_yyyyy_zzzz_1[j];

                    tg_xyyyyz_xxxx_0[j] = pb_x * tg_yyyyz_xxxx_0[j] + wp_x[j] * tg_yyyyz_xxxx_1[j] + 2.0 * fl1_fxn * tg_yyyyz_xxx_1[j];

                    tg_xyyyyz_xxxy_0[j] = pb_x * tg_yyyyz_xxxy_0[j] + wp_x[j] * tg_yyyyz_xxxy_1[j] + 1.5 * fl1_fxn * tg_yyyyz_xxy_1[j];

                    tg_xyyyyz_xxxz_0[j] = pb_x * tg_yyyyz_xxxz_0[j] + wp_x[j] * tg_yyyyz_xxxz_1[j] + 1.5 * fl1_fxn * tg_yyyyz_xxz_1[j];

                    tg_xyyyyz_xxyy_0[j] = pb_x * tg_yyyyz_xxyy_0[j] + wp_x[j] * tg_yyyyz_xxyy_1[j] + fl1_fxn * tg_yyyyz_xyy_1[j];

                    tg_xyyyyz_xxyz_0[j] = pb_x * tg_yyyyz_xxyz_0[j] + wp_x[j] * tg_yyyyz_xxyz_1[j] + fl1_fxn * tg_yyyyz_xyz_1[j];

                    tg_xyyyyz_xxzz_0[j] = pb_x * tg_yyyyz_xxzz_0[j] + wp_x[j] * tg_yyyyz_xxzz_1[j] + fl1_fxn * tg_yyyyz_xzz_1[j];

                    tg_xyyyyz_xyyy_0[j] = pb_x * tg_yyyyz_xyyy_0[j] + wp_x[j] * tg_yyyyz_xyyy_1[j] + 0.5 * fl1_fxn * tg_yyyyz_yyy_1[j];

                    tg_xyyyyz_xyyz_0[j] = pb_x * tg_yyyyz_xyyz_0[j] + wp_x[j] * tg_yyyyz_xyyz_1[j] + 0.5 * fl1_fxn * tg_yyyyz_yyz_1[j];

                    tg_xyyyyz_xyzz_0[j] = pb_x * tg_yyyyz_xyzz_0[j] + wp_x[j] * tg_yyyyz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyyyz_yzz_1[j];

                    tg_xyyyyz_xzzz_0[j] = pb_x * tg_yyyyz_xzzz_0[j] + wp_x[j] * tg_yyyyz_xzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyz_zzz_1[j];

                    tg_xyyyyz_yyyy_0[j] = pb_x * tg_yyyyz_yyyy_0[j] + wp_x[j] * tg_yyyyz_yyyy_1[j];

                    tg_xyyyyz_yyyz_0[j] = pb_x * tg_yyyyz_yyyz_0[j] + wp_x[j] * tg_yyyyz_yyyz_1[j];

                    tg_xyyyyz_yyzz_0[j] = pb_x * tg_yyyyz_yyzz_0[j] + wp_x[j] * tg_yyyyz_yyzz_1[j];

                    tg_xyyyyz_yzzz_0[j] = pb_x * tg_yyyyz_yzzz_0[j] + wp_x[j] * tg_yyyyz_yzzz_1[j];

                    tg_xyyyyz_zzzz_0[j] = pb_x * tg_yyyyz_zzzz_0[j] + wp_x[j] * tg_yyyyz_zzzz_1[j];

                    tg_xyyyzz_xxxx_0[j] = pb_x * tg_yyyzz_xxxx_0[j] + wp_x[j] * tg_yyyzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_yyyzz_xxx_1[j];

                    tg_xyyyzz_xxxy_0[j] = pb_x * tg_yyyzz_xxxy_0[j] + wp_x[j] * tg_yyyzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_yyyzz_xxy_1[j];

                    tg_xyyyzz_xxxz_0[j] = pb_x * tg_yyyzz_xxxz_0[j] + wp_x[j] * tg_yyyzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_yyyzz_xxz_1[j];

                    tg_xyyyzz_xxyy_0[j] = pb_x * tg_yyyzz_xxyy_0[j] + wp_x[j] * tg_yyyzz_xxyy_1[j] + fl1_fxn * tg_yyyzz_xyy_1[j];

                    tg_xyyyzz_xxyz_0[j] = pb_x * tg_yyyzz_xxyz_0[j] + wp_x[j] * tg_yyyzz_xxyz_1[j] + fl1_fxn * tg_yyyzz_xyz_1[j];

                    tg_xyyyzz_xxzz_0[j] = pb_x * tg_yyyzz_xxzz_0[j] + wp_x[j] * tg_yyyzz_xxzz_1[j] + fl1_fxn * tg_yyyzz_xzz_1[j];

                    tg_xyyyzz_xyyy_0[j] = pb_x * tg_yyyzz_xyyy_0[j] + wp_x[j] * tg_yyyzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_yyyzz_yyy_1[j];

                    tg_xyyyzz_xyyz_0[j] = pb_x * tg_yyyzz_xyyz_0[j] + wp_x[j] * tg_yyyzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_yyyzz_yyz_1[j];

                    tg_xyyyzz_xyzz_0[j] = pb_x * tg_yyyzz_xyzz_0[j] + wp_x[j] * tg_yyyzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyyzz_yzz_1[j];

                    tg_xyyyzz_xzzz_0[j] = pb_x * tg_yyyzz_xzzz_0[j] + wp_x[j] * tg_yyyzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_yyyzz_zzz_1[j];

                    tg_xyyyzz_yyyy_0[j] = pb_x * tg_yyyzz_yyyy_0[j] + wp_x[j] * tg_yyyzz_yyyy_1[j];

                    tg_xyyyzz_yyyz_0[j] = pb_x * tg_yyyzz_yyyz_0[j] + wp_x[j] * tg_yyyzz_yyyz_1[j];

                    tg_xyyyzz_yyzz_0[j] = pb_x * tg_yyyzz_yyzz_0[j] + wp_x[j] * tg_yyyzz_yyzz_1[j];

                    tg_xyyyzz_yzzz_0[j] = pb_x * tg_yyyzz_yzzz_0[j] + wp_x[j] * tg_yyyzz_yzzz_1[j];

                    tg_xyyyzz_zzzz_0[j] = pb_x * tg_yyyzz_zzzz_0[j] + wp_x[j] * tg_yyyzz_zzzz_1[j];

                    tg_xyyzzz_xxxx_0[j] = pb_x * tg_yyzzz_xxxx_0[j] + wp_x[j] * tg_yyzzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_yyzzz_xxx_1[j];

                    tg_xyyzzz_xxxy_0[j] = pb_x * tg_yyzzz_xxxy_0[j] + wp_x[j] * tg_yyzzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_yyzzz_xxy_1[j];

                    tg_xyyzzz_xxxz_0[j] = pb_x * tg_yyzzz_xxxz_0[j] + wp_x[j] * tg_yyzzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_yyzzz_xxz_1[j];

                    tg_xyyzzz_xxyy_0[j] = pb_x * tg_yyzzz_xxyy_0[j] + wp_x[j] * tg_yyzzz_xxyy_1[j] + fl1_fxn * tg_yyzzz_xyy_1[j];

                    tg_xyyzzz_xxyz_0[j] = pb_x * tg_yyzzz_xxyz_0[j] + wp_x[j] * tg_yyzzz_xxyz_1[j] + fl1_fxn * tg_yyzzz_xyz_1[j];

                    tg_xyyzzz_xxzz_0[j] = pb_x * tg_yyzzz_xxzz_0[j] + wp_x[j] * tg_yyzzz_xxzz_1[j] + fl1_fxn * tg_yyzzz_xzz_1[j];

                    tg_xyyzzz_xyyy_0[j] = pb_x * tg_yyzzz_xyyy_0[j] + wp_x[j] * tg_yyzzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_yyzzz_yyy_1[j];

                    tg_xyyzzz_xyyz_0[j] = pb_x * tg_yyzzz_xyyz_0[j] + wp_x[j] * tg_yyzzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_yyzzz_yyz_1[j];

                    tg_xyyzzz_xyzz_0[j] = pb_x * tg_yyzzz_xyzz_0[j] + wp_x[j] * tg_yyzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyzzz_yzz_1[j];

                    tg_xyyzzz_xzzz_0[j] = pb_x * tg_yyzzz_xzzz_0[j] + wp_x[j] * tg_yyzzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_yyzzz_zzz_1[j];

                    tg_xyyzzz_yyyy_0[j] = pb_x * tg_yyzzz_yyyy_0[j] + wp_x[j] * tg_yyzzz_yyyy_1[j];

                    tg_xyyzzz_yyyz_0[j] = pb_x * tg_yyzzz_yyyz_0[j] + wp_x[j] * tg_yyzzz_yyyz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSISG_282_328(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (282,328)

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
                                             {6, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_6_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_6_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_5_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_5_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_5_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_4_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tg_yyyyy_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 225); 

                auto tg_yyyyy_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 226); 

                auto tg_yyyyy_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 227); 

                auto tg_yyyyy_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 228); 

                auto tg_yyyyy_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 229); 

                auto tg_yyyyy_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 230); 

                auto tg_yyyyy_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 231); 

                auto tg_yyyyy_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 232); 

                auto tg_yyyyy_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 233); 

                auto tg_yyyyy_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 234); 

                auto tg_yyyyy_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 235); 

                auto tg_yyyyy_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 236); 

                auto tg_yyyyy_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 237); 

                auto tg_yyzzz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 282); 

                auto tg_yyzzz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 283); 

                auto tg_yyzzz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 284); 

                auto tg_yzzzz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 285); 

                auto tg_yzzzz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 286); 

                auto tg_yzzzz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 287); 

                auto tg_yzzzz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 288); 

                auto tg_yzzzz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 289); 

                auto tg_yzzzz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 290); 

                auto tg_yzzzz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 291); 

                auto tg_yzzzz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 292); 

                auto tg_yzzzz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 293); 

                auto tg_yzzzz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 294); 

                auto tg_yzzzz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 295); 

                auto tg_yzzzz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 296); 

                auto tg_yzzzz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 297); 

                auto tg_yzzzz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 298); 

                auto tg_yzzzz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 299); 

                auto tg_zzzzz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 300); 

                auto tg_zzzzz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 301); 

                auto tg_zzzzz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 302); 

                auto tg_zzzzz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 303); 

                auto tg_zzzzz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 304); 

                auto tg_zzzzz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 305); 

                auto tg_zzzzz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 306); 

                auto tg_zzzzz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 307); 

                auto tg_zzzzz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 308); 

                auto tg_zzzzz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 309); 

                auto tg_zzzzz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 310); 

                auto tg_zzzzz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 311); 

                auto tg_zzzzz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 312); 

                auto tg_zzzzz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 313); 

                auto tg_zzzzz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 314); 

                auto tg_yyyyy_xxxx_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 225); 

                auto tg_yyyyy_xxxy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 226); 

                auto tg_yyyyy_xxxz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 227); 

                auto tg_yyyyy_xxyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 228); 

                auto tg_yyyyy_xxyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 229); 

                auto tg_yyyyy_xxzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 230); 

                auto tg_yyyyy_xyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 231); 

                auto tg_yyyyy_xyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 232); 

                auto tg_yyyyy_xyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 233); 

                auto tg_yyyyy_xzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 234); 

                auto tg_yyyyy_yyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 235); 

                auto tg_yyyyy_yyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 236); 

                auto tg_yyyyy_yyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 237); 

                auto tg_yyzzz_yyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 282); 

                auto tg_yyzzz_yzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 283); 

                auto tg_yyzzz_zzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 284); 

                auto tg_yzzzz_xxxx_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 285); 

                auto tg_yzzzz_xxxy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 286); 

                auto tg_yzzzz_xxxz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 287); 

                auto tg_yzzzz_xxyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 288); 

                auto tg_yzzzz_xxyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 289); 

                auto tg_yzzzz_xxzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 290); 

                auto tg_yzzzz_xyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 291); 

                auto tg_yzzzz_xyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 292); 

                auto tg_yzzzz_xyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 293); 

                auto tg_yzzzz_xzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 294); 

                auto tg_yzzzz_yyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 295); 

                auto tg_yzzzz_yyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 296); 

                auto tg_yzzzz_yyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 297); 

                auto tg_yzzzz_yzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 298); 

                auto tg_yzzzz_zzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 299); 

                auto tg_zzzzz_xxxx_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 300); 

                auto tg_zzzzz_xxxy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 301); 

                auto tg_zzzzz_xxxz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 302); 

                auto tg_zzzzz_xxyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 303); 

                auto tg_zzzzz_xxyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 304); 

                auto tg_zzzzz_xxzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 305); 

                auto tg_zzzzz_xyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 306); 

                auto tg_zzzzz_xyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 307); 

                auto tg_zzzzz_xyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 308); 

                auto tg_zzzzz_xzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 309); 

                auto tg_zzzzz_yyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 310); 

                auto tg_zzzzz_yyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 311); 

                auto tg_zzzzz_yyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 312); 

                auto tg_zzzzz_yzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 313); 

                auto tg_zzzzz_zzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 314); 

                auto tg_yyyy_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 150); 

                auto tg_yyyy_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 151); 

                auto tg_yyyy_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 152); 

                auto tg_yyyy_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 153); 

                auto tg_yyyy_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 154); 

                auto tg_yyyy_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 155); 

                auto tg_yyyy_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 156); 

                auto tg_yyyy_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 157); 

                auto tg_yyyy_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 158); 

                auto tg_yyyy_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 159); 

                auto tg_yyyy_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 160); 

                auto tg_yyyy_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 161); 

                auto tg_yyyy_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 162); 

                auto tg_yyyy_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 150); 

                auto tg_yyyy_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 151); 

                auto tg_yyyy_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 152); 

                auto tg_yyyy_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 153); 

                auto tg_yyyy_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 154); 

                auto tg_yyyy_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 155); 

                auto tg_yyyy_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 156); 

                auto tg_yyyy_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 157); 

                auto tg_yyyy_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 158); 

                auto tg_yyyy_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 159); 

                auto tg_yyyy_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 160); 

                auto tg_yyyy_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 161); 

                auto tg_yyyy_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 162); 

                auto tg_yyyyy_xxx_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 150); 

                auto tg_yyyyy_xxy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 151); 

                auto tg_yyyyy_xxz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 152); 

                auto tg_yyyyy_xyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 153); 

                auto tg_yyyyy_xyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 154); 

                auto tg_yyyyy_xzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 155); 

                auto tg_yyyyy_yyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 156); 

                auto tg_yyyyy_yyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 157); 

                auto tg_yyyyy_yzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 158); 

                auto tg_yzzzz_xxx_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 190); 

                auto tg_yzzzz_xxy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 191); 

                auto tg_yzzzz_xxz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 192); 

                auto tg_yzzzz_xyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 193); 

                auto tg_yzzzz_xyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 194); 

                auto tg_yzzzz_xzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 195); 

                auto tg_yzzzz_yyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 196); 

                auto tg_yzzzz_yyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 197); 

                auto tg_yzzzz_yzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 198); 

                auto tg_yzzzz_zzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 199); 

                auto tg_zzzzz_xxx_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 200); 

                auto tg_zzzzz_xxy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 201); 

                auto tg_zzzzz_xxz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 202); 

                auto tg_zzzzz_xyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 203); 

                auto tg_zzzzz_xyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 204); 

                auto tg_zzzzz_xzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 205); 

                auto tg_zzzzz_yyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 206); 

                auto tg_zzzzz_yyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 207); 

                auto tg_zzzzz_yzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 208); 

                auto tg_zzzzz_zzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 209); 

                // set up pointers to integrals

                auto tg_xyyzzz_yyzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 282); 

                auto tg_xyyzzz_yzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 283); 

                auto tg_xyyzzz_zzzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 284); 

                auto tg_xyzzzz_xxxx_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 285); 

                auto tg_xyzzzz_xxxy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 286); 

                auto tg_xyzzzz_xxxz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 287); 

                auto tg_xyzzzz_xxyy_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 288); 

                auto tg_xyzzzz_xxyz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 289); 

                auto tg_xyzzzz_xxzz_0 = primBuffer.data(pidx_g_6_4_m0 + 420 * idx + 290); 

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

                // Batch of Integrals (282,328)

                #pragma omp simd aligned(fxn, fza, tg_xyyzzz_yyzz_0, tg_xyyzzz_yzzz_0, tg_xyyzzz_zzzz_0, \
                                         tg_xyzzzz_xxxx_0, tg_xyzzzz_xxxy_0, tg_xyzzzz_xxxz_0, tg_xyzzzz_xxyy_0, \
                                         tg_xyzzzz_xxyz_0, tg_xyzzzz_xxzz_0, tg_xyzzzz_xyyy_0, tg_xyzzzz_xyyz_0, \
                                         tg_xyzzzz_xyzz_0, tg_xyzzzz_xzzz_0, tg_xyzzzz_yyyy_0, tg_xyzzzz_yyyz_0, \
                                         tg_xyzzzz_yyzz_0, tg_xyzzzz_yzzz_0, tg_xyzzzz_zzzz_0, tg_xzzzzz_xxxx_0, \
                                         tg_xzzzzz_xxxy_0, tg_xzzzzz_xxxz_0, tg_xzzzzz_xxyy_0, tg_xzzzzz_xxyz_0, \
                                         tg_xzzzzz_xxzz_0, tg_xzzzzz_xyyy_0, tg_xzzzzz_xyyz_0, tg_xzzzzz_xyzz_0, \
                                         tg_xzzzzz_xzzz_0, tg_xzzzzz_yyyy_0, tg_xzzzzz_yyyz_0, tg_xzzzzz_yyzz_0, \
                                         tg_xzzzzz_yzzz_0, tg_xzzzzz_zzzz_0, tg_yyyy_xxxx_0, tg_yyyy_xxxx_1, tg_yyyy_xxxy_0, \
                                         tg_yyyy_xxxy_1, tg_yyyy_xxxz_0, tg_yyyy_xxxz_1, tg_yyyy_xxyy_0, tg_yyyy_xxyy_1, \
                                         tg_yyyy_xxyz_0, tg_yyyy_xxyz_1, tg_yyyy_xxzz_0, tg_yyyy_xxzz_1, tg_yyyy_xyyy_0, \
                                         tg_yyyy_xyyy_1, tg_yyyy_xyyz_0, tg_yyyy_xyyz_1, tg_yyyy_xyzz_0, tg_yyyy_xyzz_1, \
                                         tg_yyyy_xzzz_0, tg_yyyy_xzzz_1, tg_yyyy_yyyy_0, tg_yyyy_yyyy_1, tg_yyyy_yyyz_0, \
                                         tg_yyyy_yyyz_1, tg_yyyy_yyzz_0, tg_yyyy_yyzz_1, tg_yyyyy_xxx_1, tg_yyyyy_xxxx_0, \
                                         tg_yyyyy_xxxx_1, tg_yyyyy_xxxy_0, tg_yyyyy_xxxy_1, tg_yyyyy_xxxz_0, tg_yyyyy_xxxz_1, \
                                         tg_yyyyy_xxy_1, tg_yyyyy_xxyy_0, tg_yyyyy_xxyy_1, tg_yyyyy_xxyz_0, tg_yyyyy_xxyz_1, \
                                         tg_yyyyy_xxz_1, tg_yyyyy_xxzz_0, tg_yyyyy_xxzz_1, tg_yyyyy_xyy_1, tg_yyyyy_xyyy_0, \
                                         tg_yyyyy_xyyy_1, tg_yyyyy_xyyz_0, tg_yyyyy_xyyz_1, tg_yyyyy_xyz_1, tg_yyyyy_xyzz_0, \
                                         tg_yyyyy_xyzz_1, tg_yyyyy_xzz_1, tg_yyyyy_xzzz_0, tg_yyyyy_xzzz_1, tg_yyyyy_yyy_1, \
                                         tg_yyyyy_yyyy_0, tg_yyyyy_yyyy_1, tg_yyyyy_yyyz_0, tg_yyyyy_yyyz_1, tg_yyyyy_yyz_1, \
                                         tg_yyyyy_yyzz_0, tg_yyyyy_yyzz_1, tg_yyyyy_yzz_1, tg_yyyyyy_xxxx_0, \
                                         tg_yyyyyy_xxxy_0, tg_yyyyyy_xxxz_0, tg_yyyyyy_xxyy_0, tg_yyyyyy_xxyz_0, \
                                         tg_yyyyyy_xxzz_0, tg_yyyyyy_xyyy_0, tg_yyyyyy_xyyz_0, tg_yyyyyy_xyzz_0, \
                                         tg_yyyyyy_xzzz_0, tg_yyyyyy_yyyy_0, tg_yyyyyy_yyyz_0, tg_yyyyyy_yyzz_0, \
                                         tg_yyzzz_yyzz_0, tg_yyzzz_yyzz_1, tg_yyzzz_yzzz_0, tg_yyzzz_yzzz_1, tg_yyzzz_zzzz_0, \
                                         tg_yyzzz_zzzz_1, tg_yzzzz_xxx_1, tg_yzzzz_xxxx_0, tg_yzzzz_xxxx_1, tg_yzzzz_xxxy_0, \
                                         tg_yzzzz_xxxy_1, tg_yzzzz_xxxz_0, tg_yzzzz_xxxz_1, tg_yzzzz_xxy_1, tg_yzzzz_xxyy_0, \
                                         tg_yzzzz_xxyy_1, tg_yzzzz_xxyz_0, tg_yzzzz_xxyz_1, tg_yzzzz_xxz_1, tg_yzzzz_xxzz_0, \
                                         tg_yzzzz_xxzz_1, tg_yzzzz_xyy_1, tg_yzzzz_xyyy_0, tg_yzzzz_xyyy_1, tg_yzzzz_xyyz_0, \
                                         tg_yzzzz_xyyz_1, tg_yzzzz_xyz_1, tg_yzzzz_xyzz_0, tg_yzzzz_xyzz_1, tg_yzzzz_xzz_1, \
                                         tg_yzzzz_xzzz_0, tg_yzzzz_xzzz_1, tg_yzzzz_yyy_1, tg_yzzzz_yyyy_0, tg_yzzzz_yyyy_1, \
                                         tg_yzzzz_yyyz_0, tg_yzzzz_yyyz_1, tg_yzzzz_yyz_1, tg_yzzzz_yyzz_0, tg_yzzzz_yyzz_1, \
                                         tg_yzzzz_yzz_1, tg_yzzzz_yzzz_0, tg_yzzzz_yzzz_1, tg_yzzzz_zzz_1, tg_yzzzz_zzzz_0, \
                                         tg_yzzzz_zzzz_1, tg_zzzzz_xxx_1, tg_zzzzz_xxxx_0, tg_zzzzz_xxxx_1, tg_zzzzz_xxxy_0, \
                                         tg_zzzzz_xxxy_1, tg_zzzzz_xxxz_0, tg_zzzzz_xxxz_1, tg_zzzzz_xxy_1, tg_zzzzz_xxyy_0, \
                                         tg_zzzzz_xxyy_1, tg_zzzzz_xxyz_0, tg_zzzzz_xxyz_1, tg_zzzzz_xxz_1, tg_zzzzz_xxzz_0, \
                                         tg_zzzzz_xxzz_1, tg_zzzzz_xyy_1, tg_zzzzz_xyyy_0, tg_zzzzz_xyyy_1, tg_zzzzz_xyyz_0, \
                                         tg_zzzzz_xyyz_1, tg_zzzzz_xyz_1, tg_zzzzz_xyzz_0, tg_zzzzz_xyzz_1, tg_zzzzz_xzz_1, \
                                         tg_zzzzz_xzzz_0, tg_zzzzz_xzzz_1, tg_zzzzz_yyy_1, tg_zzzzz_yyyy_0, tg_zzzzz_yyyy_1, \
                                         tg_zzzzz_yyyz_0, tg_zzzzz_yyyz_1, tg_zzzzz_yyz_1, tg_zzzzz_yyzz_0, tg_zzzzz_yyzz_1, \
                                         tg_zzzzz_yzz_1, tg_zzzzz_yzzz_0, tg_zzzzz_yzzz_1, tg_zzzzz_zzz_1, tg_zzzzz_zzzz_0, \
                                         tg_zzzzz_zzzz_1, wp_x, wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xyyzzz_yyzz_0[j] = pb_x * tg_yyzzz_yyzz_0[j] + wp_x[j] * tg_yyzzz_yyzz_1[j];

                    tg_xyyzzz_yzzz_0[j] = pb_x * tg_yyzzz_yzzz_0[j] + wp_x[j] * tg_yyzzz_yzzz_1[j];

                    tg_xyyzzz_zzzz_0[j] = pb_x * tg_yyzzz_zzzz_0[j] + wp_x[j] * tg_yyzzz_zzzz_1[j];

                    tg_xyzzzz_xxxx_0[j] = pb_x * tg_yzzzz_xxxx_0[j] + wp_x[j] * tg_yzzzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_yzzzz_xxx_1[j];

                    tg_xyzzzz_xxxy_0[j] = pb_x * tg_yzzzz_xxxy_0[j] + wp_x[j] * tg_yzzzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_yzzzz_xxy_1[j];

                    tg_xyzzzz_xxxz_0[j] = pb_x * tg_yzzzz_xxxz_0[j] + wp_x[j] * tg_yzzzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_yzzzz_xxz_1[j];

                    tg_xyzzzz_xxyy_0[j] = pb_x * tg_yzzzz_xxyy_0[j] + wp_x[j] * tg_yzzzz_xxyy_1[j] + fl1_fxn * tg_yzzzz_xyy_1[j];

                    tg_xyzzzz_xxyz_0[j] = pb_x * tg_yzzzz_xxyz_0[j] + wp_x[j] * tg_yzzzz_xxyz_1[j] + fl1_fxn * tg_yzzzz_xyz_1[j];

                    tg_xyzzzz_xxzz_0[j] = pb_x * tg_yzzzz_xxzz_0[j] + wp_x[j] * tg_yzzzz_xxzz_1[j] + fl1_fxn * tg_yzzzz_xzz_1[j];

                    tg_xyzzzz_xyyy_0[j] = pb_x * tg_yzzzz_xyyy_0[j] + wp_x[j] * tg_yzzzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_yzzzz_yyy_1[j];

                    tg_xyzzzz_xyyz_0[j] = pb_x * tg_yzzzz_xyyz_0[j] + wp_x[j] * tg_yzzzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_yzzzz_yyz_1[j];

                    tg_xyzzzz_xyzz_0[j] = pb_x * tg_yzzzz_xyzz_0[j] + wp_x[j] * tg_yzzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yzzzz_yzz_1[j];

                    tg_xyzzzz_xzzz_0[j] = pb_x * tg_yzzzz_xzzz_0[j] + wp_x[j] * tg_yzzzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_yzzzz_zzz_1[j];

                    tg_xyzzzz_yyyy_0[j] = pb_x * tg_yzzzz_yyyy_0[j] + wp_x[j] * tg_yzzzz_yyyy_1[j];

                    tg_xyzzzz_yyyz_0[j] = pb_x * tg_yzzzz_yyyz_0[j] + wp_x[j] * tg_yzzzz_yyyz_1[j];

                    tg_xyzzzz_yyzz_0[j] = pb_x * tg_yzzzz_yyzz_0[j] + wp_x[j] * tg_yzzzz_yyzz_1[j];

                    tg_xyzzzz_yzzz_0[j] = pb_x * tg_yzzzz_yzzz_0[j] + wp_x[j] * tg_yzzzz_yzzz_1[j];

                    tg_xyzzzz_zzzz_0[j] = pb_x * tg_yzzzz_zzzz_0[j] + wp_x[j] * tg_yzzzz_zzzz_1[j];

                    tg_xzzzzz_xxxx_0[j] = pb_x * tg_zzzzz_xxxx_0[j] + wp_x[j] * tg_zzzzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_zzzzz_xxx_1[j];

                    tg_xzzzzz_xxxy_0[j] = pb_x * tg_zzzzz_xxxy_0[j] + wp_x[j] * tg_zzzzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_zzzzz_xxy_1[j];

                    tg_xzzzzz_xxxz_0[j] = pb_x * tg_zzzzz_xxxz_0[j] + wp_x[j] * tg_zzzzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_zzzzz_xxz_1[j];

                    tg_xzzzzz_xxyy_0[j] = pb_x * tg_zzzzz_xxyy_0[j] + wp_x[j] * tg_zzzzz_xxyy_1[j] + fl1_fxn * tg_zzzzz_xyy_1[j];

                    tg_xzzzzz_xxyz_0[j] = pb_x * tg_zzzzz_xxyz_0[j] + wp_x[j] * tg_zzzzz_xxyz_1[j] + fl1_fxn * tg_zzzzz_xyz_1[j];

                    tg_xzzzzz_xxzz_0[j] = pb_x * tg_zzzzz_xxzz_0[j] + wp_x[j] * tg_zzzzz_xxzz_1[j] + fl1_fxn * tg_zzzzz_xzz_1[j];

                    tg_xzzzzz_xyyy_0[j] = pb_x * tg_zzzzz_xyyy_0[j] + wp_x[j] * tg_zzzzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_zzzzz_yyy_1[j];

                    tg_xzzzzz_xyyz_0[j] = pb_x * tg_zzzzz_xyyz_0[j] + wp_x[j] * tg_zzzzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_zzzzz_yyz_1[j];

                    tg_xzzzzz_xyzz_0[j] = pb_x * tg_zzzzz_xyzz_0[j] + wp_x[j] * tg_zzzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_zzzzz_yzz_1[j];

                    tg_xzzzzz_xzzz_0[j] = pb_x * tg_zzzzz_xzzz_0[j] + wp_x[j] * tg_zzzzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzz_zzz_1[j];

                    tg_xzzzzz_yyyy_0[j] = pb_x * tg_zzzzz_yyyy_0[j] + wp_x[j] * tg_zzzzz_yyyy_1[j];

                    tg_xzzzzz_yyyz_0[j] = pb_x * tg_zzzzz_yyyz_0[j] + wp_x[j] * tg_zzzzz_yyyz_1[j];

                    tg_xzzzzz_yyzz_0[j] = pb_x * tg_zzzzz_yyzz_0[j] + wp_x[j] * tg_zzzzz_yyzz_1[j];

                    tg_xzzzzz_yzzz_0[j] = pb_x * tg_zzzzz_yzzz_0[j] + wp_x[j] * tg_zzzzz_yzzz_1[j];

                    tg_xzzzzz_zzzz_0[j] = pb_x * tg_zzzzz_zzzz_0[j] + wp_x[j] * tg_zzzzz_zzzz_1[j];

                    tg_yyyyyy_xxxx_0[j] = pb_y * tg_yyyyy_xxxx_0[j] + wp_y[j] * tg_yyyyy_xxxx_1[j] + 2.5 * fl1_fx * tg_yyyy_xxxx_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyy_xxxx_1[j];

                    tg_yyyyyy_xxxy_0[j] = pb_y * tg_yyyyy_xxxy_0[j] + wp_y[j] * tg_yyyyy_xxxy_1[j] + 2.5 * fl1_fx * tg_yyyy_xxxy_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyy_xxxy_1[j] + 0.5 * fl1_fxn * tg_yyyyy_xxx_1[j];

                    tg_yyyyyy_xxxz_0[j] = pb_y * tg_yyyyy_xxxz_0[j] + wp_y[j] * tg_yyyyy_xxxz_1[j] + 2.5 * fl1_fx * tg_yyyy_xxxz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyy_xxxz_1[j];

                    tg_yyyyyy_xxyy_0[j] = pb_y * tg_yyyyy_xxyy_0[j] + wp_y[j] * tg_yyyyy_xxyy_1[j] + 2.5 * fl1_fx * tg_yyyy_xxyy_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyy_xxyy_1[j] + fl1_fxn * tg_yyyyy_xxy_1[j];

                    tg_yyyyyy_xxyz_0[j] = pb_y * tg_yyyyy_xxyz_0[j] + wp_y[j] * tg_yyyyy_xxyz_1[j] + 2.5 * fl1_fx * tg_yyyy_xxyz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyy_xxyz_1[j] + 0.5 * fl1_fxn * tg_yyyyy_xxz_1[j];

                    tg_yyyyyy_xxzz_0[j] = pb_y * tg_yyyyy_xxzz_0[j] + wp_y[j] * tg_yyyyy_xxzz_1[j] + 2.5 * fl1_fx * tg_yyyy_xxzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyy_xxzz_1[j];

                    tg_yyyyyy_xyyy_0[j] = pb_y * tg_yyyyy_xyyy_0[j] + wp_y[j] * tg_yyyyy_xyyy_1[j] + 2.5 * fl1_fx * tg_yyyy_xyyy_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyy_xyyy_1[j] + 1.5 * fl1_fxn * tg_yyyyy_xyy_1[j];

                    tg_yyyyyy_xyyz_0[j] = pb_y * tg_yyyyy_xyyz_0[j] + wp_y[j] * tg_yyyyy_xyyz_1[j] + 2.5 * fl1_fx * tg_yyyy_xyyz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyy_xyyz_1[j] + fl1_fxn * tg_yyyyy_xyz_1[j];

                    tg_yyyyyy_xyzz_0[j] = pb_y * tg_yyyyy_xyzz_0[j] + wp_y[j] * tg_yyyyy_xyzz_1[j] + 2.5 * fl1_fx * tg_yyyy_xyzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyy_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyyyy_xzz_1[j];

                    tg_yyyyyy_xzzz_0[j] = pb_y * tg_yyyyy_xzzz_0[j] + wp_y[j] * tg_yyyyy_xzzz_1[j] + 2.5 * fl1_fx * tg_yyyy_xzzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyy_xzzz_1[j];

                    tg_yyyyyy_yyyy_0[j] = pb_y * tg_yyyyy_yyyy_0[j] + wp_y[j] * tg_yyyyy_yyyy_1[j] + 2.5 * fl1_fx * tg_yyyy_yyyy_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyy_yyyy_1[j] + 2.0 * fl1_fxn * tg_yyyyy_yyy_1[j];

                    tg_yyyyyy_yyyz_0[j] = pb_y * tg_yyyyy_yyyz_0[j] + wp_y[j] * tg_yyyyy_yyyz_1[j] + 2.5 * fl1_fx * tg_yyyy_yyyz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyy_yyyz_1[j] + 1.5 * fl1_fxn * tg_yyyyy_yyz_1[j];

                    tg_yyyyyy_yyzz_0[j] = pb_y * tg_yyyyy_yyzz_0[j] + wp_y[j] * tg_yyyyy_yyzz_1[j] + 2.5 * fl1_fx * tg_yyyy_yyzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyy_yyzz_1[j] + fl1_fxn * tg_yyyyy_yzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSISG_328_374(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (328,374)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {6, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_6_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_6_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_5_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_5_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_4_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_5_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_yyyyy_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 238); 

                auto tg_yyyyy_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 239); 

                auto tg_yyyyz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 240); 

                auto tg_yyyyz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 241); 

                auto tg_yyyyz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 242); 

                auto tg_yyyyz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 243); 

                auto tg_yyyyz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 244); 

                auto tg_yyyyz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 245); 

                auto tg_yyyyz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 246); 

                auto tg_yyyyz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 247); 

                auto tg_yyyyz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 248); 

                auto tg_yyyyz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 249); 

                auto tg_yyyyz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 250); 

                auto tg_yyyyz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 251); 

                auto tg_yyyyz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 252); 

                auto tg_yyyyz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 253); 

                auto tg_yyyyz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 254); 

                auto tg_yyyzz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 255); 

                auto tg_yyyzz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 256); 

                auto tg_yyyzz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 257); 

                auto tg_yyyzz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 258); 

                auto tg_yyyzz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 259); 

                auto tg_yyyzz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 260); 

                auto tg_yyyzz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 261); 

                auto tg_yyyzz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 262); 

                auto tg_yyyzz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 263); 

                auto tg_yyyzz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 264); 

                auto tg_yyyzz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 265); 

                auto tg_yyyzz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 266); 

                auto tg_yyyzz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 267); 

                auto tg_yyyzz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 268); 

                auto tg_yyyzz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 269); 

                auto tg_yyzzz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 270); 

                auto tg_yyzzz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 271); 

                auto tg_yyzzz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 272); 

                auto tg_yyzzz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 273); 

                auto tg_yyzzz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 274); 

                auto tg_yyzzz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 275); 

                auto tg_yyzzz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 276); 

                auto tg_yyzzz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 277); 

                auto tg_yyzzz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 278); 

                auto tg_yyzzz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 279); 

                auto tg_yyzzz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 280); 

                auto tg_yyzzz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 281); 

                auto tg_yyzzz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 282); 

                auto tg_yyzzz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 283); 

                auto tg_yyyyy_yzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 238); 

                auto tg_yyyyy_zzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 239); 

                auto tg_yyyyz_xxxx_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 240); 

                auto tg_yyyyz_xxxy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 241); 

                auto tg_yyyyz_xxxz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 242); 

                auto tg_yyyyz_xxyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 243); 

                auto tg_yyyyz_xxyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 244); 

                auto tg_yyyyz_xxzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 245); 

                auto tg_yyyyz_xyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 246); 

                auto tg_yyyyz_xyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 247); 

                auto tg_yyyyz_xyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 248); 

                auto tg_yyyyz_xzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 249); 

                auto tg_yyyyz_yyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 250); 

                auto tg_yyyyz_yyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 251); 

                auto tg_yyyyz_yyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 252); 

                auto tg_yyyyz_yzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 253); 

                auto tg_yyyyz_zzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 254); 

                auto tg_yyyzz_xxxx_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 255); 

                auto tg_yyyzz_xxxy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 256); 

                auto tg_yyyzz_xxxz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 257); 

                auto tg_yyyzz_xxyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 258); 

                auto tg_yyyzz_xxyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 259); 

                auto tg_yyyzz_xxzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 260); 

                auto tg_yyyzz_xyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 261); 

                auto tg_yyyzz_xyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 262); 

                auto tg_yyyzz_xyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 263); 

                auto tg_yyyzz_xzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 264); 

                auto tg_yyyzz_yyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 265); 

                auto tg_yyyzz_yyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 266); 

                auto tg_yyyzz_yyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 267); 

                auto tg_yyyzz_yzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 268); 

                auto tg_yyyzz_zzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 269); 

                auto tg_yyzzz_xxxx_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 270); 

                auto tg_yyzzz_xxxy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 271); 

                auto tg_yyzzz_xxxz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 272); 

                auto tg_yyzzz_xxyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 273); 

                auto tg_yyzzz_xxyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 274); 

                auto tg_yyzzz_xxzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 275); 

                auto tg_yyzzz_xyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 276); 

                auto tg_yyzzz_xyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 277); 

                auto tg_yyzzz_xyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 278); 

                auto tg_yyzzz_xzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 279); 

                auto tg_yyzzz_yyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 280); 

                auto tg_yyzzz_yyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 281); 

                auto tg_yyzzz_yyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 282); 

                auto tg_yyzzz_yzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 283); 

                auto tg_yyyy_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 163); 

                auto tg_yyyy_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 164); 

                auto tg_yyyz_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 165); 

                auto tg_yyyz_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 166); 

                auto tg_yyyz_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 167); 

                auto tg_yyyz_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 168); 

                auto tg_yyyz_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 169); 

                auto tg_yyyz_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 170); 

                auto tg_yyyz_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 171); 

                auto tg_yyyz_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 172); 

                auto tg_yyyz_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 173); 

                auto tg_yyyz_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 174); 

                auto tg_yyyz_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 175); 

                auto tg_yyyz_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 176); 

                auto tg_yyyz_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 177); 

                auto tg_yyyz_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 178); 

                auto tg_yyyz_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 179); 

                auto tg_yyzz_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 180); 

                auto tg_yyzz_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 181); 

                auto tg_yyzz_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 182); 

                auto tg_yyzz_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 183); 

                auto tg_yyzz_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 184); 

                auto tg_yyzz_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 185); 

                auto tg_yyzz_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 186); 

                auto tg_yyzz_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 187); 

                auto tg_yyzz_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 188); 

                auto tg_yyzz_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 189); 

                auto tg_yyzz_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 190); 

                auto tg_yyzz_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 191); 

                auto tg_yyzz_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 192); 

                auto tg_yyzz_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 193); 

                auto tg_yyzz_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 194); 

                auto tg_yzzz_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 195); 

                auto tg_yzzz_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 196); 

                auto tg_yzzz_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 197); 

                auto tg_yzzz_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 198); 

                auto tg_yzzz_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 199); 

                auto tg_yzzz_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 200); 

                auto tg_yzzz_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 201); 

                auto tg_yzzz_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 202); 

                auto tg_yzzz_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 203); 

                auto tg_yzzz_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 204); 

                auto tg_yzzz_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 205); 

                auto tg_yzzz_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 206); 

                auto tg_yzzz_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 207); 

                auto tg_yzzz_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 208); 

                auto tg_yyyy_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 163); 

                auto tg_yyyy_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 164); 

                auto tg_yyyz_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 165); 

                auto tg_yyyz_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 166); 

                auto tg_yyyz_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 167); 

                auto tg_yyyz_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 168); 

                auto tg_yyyz_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 169); 

                auto tg_yyyz_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 170); 

                auto tg_yyyz_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 171); 

                auto tg_yyyz_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 172); 

                auto tg_yyyz_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 173); 

                auto tg_yyyz_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 174); 

                auto tg_yyyz_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 175); 

                auto tg_yyyz_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 176); 

                auto tg_yyyz_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 177); 

                auto tg_yyyz_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 178); 

                auto tg_yyyz_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 179); 

                auto tg_yyzz_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 180); 

                auto tg_yyzz_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 181); 

                auto tg_yyzz_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 182); 

                auto tg_yyzz_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 183); 

                auto tg_yyzz_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 184); 

                auto tg_yyzz_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 185); 

                auto tg_yyzz_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 186); 

                auto tg_yyzz_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 187); 

                auto tg_yyzz_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 188); 

                auto tg_yyzz_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 189); 

                auto tg_yyzz_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 190); 

                auto tg_yyzz_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 191); 

                auto tg_yyzz_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 192); 

                auto tg_yyzz_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 193); 

                auto tg_yyzz_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 194); 

                auto tg_yzzz_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 195); 

                auto tg_yzzz_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 196); 

                auto tg_yzzz_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 197); 

                auto tg_yzzz_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 198); 

                auto tg_yzzz_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 199); 

                auto tg_yzzz_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 200); 

                auto tg_yzzz_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 201); 

                auto tg_yzzz_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 202); 

                auto tg_yzzz_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 203); 

                auto tg_yzzz_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 204); 

                auto tg_yzzz_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 205); 

                auto tg_yzzz_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 206); 

                auto tg_yzzz_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 207); 

                auto tg_yzzz_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 208); 

                auto tg_yyyyy_zzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 159); 

                auto tg_yyyyz_xxx_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 160); 

                auto tg_yyyyz_xxy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 161); 

                auto tg_yyyyz_xxz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 162); 

                auto tg_yyyyz_xyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 163); 

                auto tg_yyyyz_xyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 164); 

                auto tg_yyyyz_xzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 165); 

                auto tg_yyyyz_yyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 166); 

                auto tg_yyyyz_yyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 167); 

                auto tg_yyyyz_yzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 168); 

                auto tg_yyyyz_zzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 169); 

                auto tg_yyyzz_xxx_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 170); 

                auto tg_yyyzz_xxy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 171); 

                auto tg_yyyzz_xxz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 172); 

                auto tg_yyyzz_xyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 173); 

                auto tg_yyyzz_xyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 174); 

                auto tg_yyyzz_xzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 175); 

                auto tg_yyyzz_yyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 176); 

                auto tg_yyyzz_yyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 177); 

                auto tg_yyyzz_yzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 178); 

                auto tg_yyyzz_zzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 179); 

                auto tg_yyzzz_xxx_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 180); 

                auto tg_yyzzz_xxy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 181); 

                auto tg_yyzzz_xxz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 182); 

                auto tg_yyzzz_xyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 183); 

                auto tg_yyzzz_xyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 184); 

                auto tg_yyzzz_xzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 185); 

                auto tg_yyzzz_yyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 186); 

                auto tg_yyzzz_yyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 187); 

                auto tg_yyzzz_yzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 188); 

                auto tg_yyzzz_zzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 189); 

                // set up pointers to integrals

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

                // Batch of Integrals (328,374)

                #pragma omp simd aligned(fxn, fza, tg_yyyy_yzzz_0, tg_yyyy_yzzz_1, tg_yyyy_zzzz_0, \
                                         tg_yyyy_zzzz_1, tg_yyyyy_yzzz_0, tg_yyyyy_yzzz_1, tg_yyyyy_zzz_1, tg_yyyyy_zzzz_0, \
                                         tg_yyyyy_zzzz_1, tg_yyyyyy_yzzz_0, tg_yyyyyy_zzzz_0, tg_yyyyyz_xxxx_0, \
                                         tg_yyyyyz_xxxy_0, tg_yyyyyz_xxxz_0, tg_yyyyyz_xxyy_0, tg_yyyyyz_xxyz_0, \
                                         tg_yyyyyz_xxzz_0, tg_yyyyyz_xyyy_0, tg_yyyyyz_xyyz_0, tg_yyyyyz_xyzz_0, \
                                         tg_yyyyyz_xzzz_0, tg_yyyyyz_yyyy_0, tg_yyyyyz_yyyz_0, tg_yyyyyz_yyzz_0, \
                                         tg_yyyyyz_yzzz_0, tg_yyyyyz_zzzz_0, tg_yyyyz_xxx_1, tg_yyyyz_xxxx_0, tg_yyyyz_xxxx_1, \
                                         tg_yyyyz_xxxy_0, tg_yyyyz_xxxy_1, tg_yyyyz_xxxz_0, tg_yyyyz_xxxz_1, tg_yyyyz_xxy_1, \
                                         tg_yyyyz_xxyy_0, tg_yyyyz_xxyy_1, tg_yyyyz_xxyz_0, tg_yyyyz_xxyz_1, tg_yyyyz_xxz_1, \
                                         tg_yyyyz_xxzz_0, tg_yyyyz_xxzz_1, tg_yyyyz_xyy_1, tg_yyyyz_xyyy_0, tg_yyyyz_xyyy_1, \
                                         tg_yyyyz_xyyz_0, tg_yyyyz_xyyz_1, tg_yyyyz_xyz_1, tg_yyyyz_xyzz_0, tg_yyyyz_xyzz_1, \
                                         tg_yyyyz_xzz_1, tg_yyyyz_xzzz_0, tg_yyyyz_xzzz_1, tg_yyyyz_yyy_1, tg_yyyyz_yyyy_0, \
                                         tg_yyyyz_yyyy_1, tg_yyyyz_yyyz_0, tg_yyyyz_yyyz_1, tg_yyyyz_yyz_1, tg_yyyyz_yyzz_0, \
                                         tg_yyyyz_yyzz_1, tg_yyyyz_yzz_1, tg_yyyyz_yzzz_0, tg_yyyyz_yzzz_1, tg_yyyyz_zzz_1, \
                                         tg_yyyyz_zzzz_0, tg_yyyyz_zzzz_1, tg_yyyyzz_xxxx_0, tg_yyyyzz_xxxy_0, \
                                         tg_yyyyzz_xxxz_0, tg_yyyyzz_xxyy_0, tg_yyyyzz_xxyz_0, tg_yyyyzz_xxzz_0, \
                                         tg_yyyyzz_xyyy_0, tg_yyyyzz_xyyz_0, tg_yyyyzz_xyzz_0, tg_yyyyzz_xzzz_0, \
                                         tg_yyyyzz_yyyy_0, tg_yyyyzz_yyyz_0, tg_yyyyzz_yyzz_0, tg_yyyyzz_yzzz_0, \
                                         tg_yyyyzz_zzzz_0, tg_yyyz_xxxx_0, tg_yyyz_xxxx_1, tg_yyyz_xxxy_0, tg_yyyz_xxxy_1, \
                                         tg_yyyz_xxxz_0, tg_yyyz_xxxz_1, tg_yyyz_xxyy_0, tg_yyyz_xxyy_1, tg_yyyz_xxyz_0, \
                                         tg_yyyz_xxyz_1, tg_yyyz_xxzz_0, tg_yyyz_xxzz_1, tg_yyyz_xyyy_0, tg_yyyz_xyyy_1, \
                                         tg_yyyz_xyyz_0, tg_yyyz_xyyz_1, tg_yyyz_xyzz_0, tg_yyyz_xyzz_1, tg_yyyz_xzzz_0, \
                                         tg_yyyz_xzzz_1, tg_yyyz_yyyy_0, tg_yyyz_yyyy_1, tg_yyyz_yyyz_0, tg_yyyz_yyyz_1, \
                                         tg_yyyz_yyzz_0, tg_yyyz_yyzz_1, tg_yyyz_yzzz_0, tg_yyyz_yzzz_1, tg_yyyz_zzzz_0, \
                                         tg_yyyz_zzzz_1, tg_yyyzz_xxx_1, tg_yyyzz_xxxx_0, tg_yyyzz_xxxx_1, tg_yyyzz_xxxy_0, \
                                         tg_yyyzz_xxxy_1, tg_yyyzz_xxxz_0, tg_yyyzz_xxxz_1, tg_yyyzz_xxy_1, tg_yyyzz_xxyy_0, \
                                         tg_yyyzz_xxyy_1, tg_yyyzz_xxyz_0, tg_yyyzz_xxyz_1, tg_yyyzz_xxz_1, tg_yyyzz_xxzz_0, \
                                         tg_yyyzz_xxzz_1, tg_yyyzz_xyy_1, tg_yyyzz_xyyy_0, tg_yyyzz_xyyy_1, tg_yyyzz_xyyz_0, \
                                         tg_yyyzz_xyyz_1, tg_yyyzz_xyz_1, tg_yyyzz_xyzz_0, tg_yyyzz_xyzz_1, tg_yyyzz_xzz_1, \
                                         tg_yyyzz_xzzz_0, tg_yyyzz_xzzz_1, tg_yyyzz_yyy_1, tg_yyyzz_yyyy_0, tg_yyyzz_yyyy_1, \
                                         tg_yyyzz_yyyz_0, tg_yyyzz_yyyz_1, tg_yyyzz_yyz_1, tg_yyyzz_yyzz_0, tg_yyyzz_yyzz_1, \
                                         tg_yyyzz_yzz_1, tg_yyyzz_yzzz_0, tg_yyyzz_yzzz_1, tg_yyyzz_zzz_1, tg_yyyzz_zzzz_0, \
                                         tg_yyyzz_zzzz_1, tg_yyyzzz_xxxx_0, tg_yyyzzz_xxxy_0, tg_yyyzzz_xxxz_0, \
                                         tg_yyyzzz_xxyy_0, tg_yyyzzz_xxyz_0, tg_yyyzzz_xxzz_0, tg_yyyzzz_xyyy_0, \
                                         tg_yyyzzz_xyyz_0, tg_yyyzzz_xyzz_0, tg_yyyzzz_xzzz_0, tg_yyyzzz_yyyy_0, \
                                         tg_yyyzzz_yyyz_0, tg_yyyzzz_yyzz_0, tg_yyyzzz_yzzz_0, tg_yyzz_xxxx_0, tg_yyzz_xxxx_1, \
                                         tg_yyzz_xxxy_0, tg_yyzz_xxxy_1, tg_yyzz_xxxz_0, tg_yyzz_xxxz_1, tg_yyzz_xxyy_0, \
                                         tg_yyzz_xxyy_1, tg_yyzz_xxyz_0, tg_yyzz_xxyz_1, tg_yyzz_xxzz_0, tg_yyzz_xxzz_1, \
                                         tg_yyzz_xyyy_0, tg_yyzz_xyyy_1, tg_yyzz_xyyz_0, tg_yyzz_xyyz_1, tg_yyzz_xyzz_0, \
                                         tg_yyzz_xyzz_1, tg_yyzz_xzzz_0, tg_yyzz_xzzz_1, tg_yyzz_yyyy_0, tg_yyzz_yyyy_1, \
                                         tg_yyzz_yyyz_0, tg_yyzz_yyyz_1, tg_yyzz_yyzz_0, tg_yyzz_yyzz_1, tg_yyzz_yzzz_0, \
                                         tg_yyzz_yzzz_1, tg_yyzz_zzzz_0, tg_yyzz_zzzz_1, tg_yyzzz_xxx_1, tg_yyzzz_xxxx_0, \
                                         tg_yyzzz_xxxx_1, tg_yyzzz_xxxy_0, tg_yyzzz_xxxy_1, tg_yyzzz_xxxz_0, tg_yyzzz_xxxz_1, \
                                         tg_yyzzz_xxy_1, tg_yyzzz_xxyy_0, tg_yyzzz_xxyy_1, tg_yyzzz_xxyz_0, tg_yyzzz_xxyz_1, \
                                         tg_yyzzz_xxz_1, tg_yyzzz_xxzz_0, tg_yyzzz_xxzz_1, tg_yyzzz_xyy_1, tg_yyzzz_xyyy_0, \
                                         tg_yyzzz_xyyy_1, tg_yyzzz_xyyz_0, tg_yyzzz_xyyz_1, tg_yyzzz_xyz_1, tg_yyzzz_xyzz_0, \
                                         tg_yyzzz_xyzz_1, tg_yyzzz_xzz_1, tg_yyzzz_xzzz_0, tg_yyzzz_xzzz_1, tg_yyzzz_yyy_1, \
                                         tg_yyzzz_yyyy_0, tg_yyzzz_yyyy_1, tg_yyzzz_yyyz_0, tg_yyzzz_yyyz_1, tg_yyzzz_yyz_1, \
                                         tg_yyzzz_yyzz_0, tg_yyzzz_yyzz_1, tg_yyzzz_yzz_1, tg_yyzzz_yzzz_0, tg_yyzzz_yzzz_1, \
                                         tg_yyzzz_zzz_1, tg_yzzz_xxxx_0, tg_yzzz_xxxx_1, tg_yzzz_xxxy_0, tg_yzzz_xxxy_1, \
                                         tg_yzzz_xxxz_0, tg_yzzz_xxxz_1, tg_yzzz_xxyy_0, tg_yzzz_xxyy_1, tg_yzzz_xxyz_0, \
                                         tg_yzzz_xxyz_1, tg_yzzz_xxzz_0, tg_yzzz_xxzz_1, tg_yzzz_xyyy_0, tg_yzzz_xyyy_1, \
                                         tg_yzzz_xyyz_0, tg_yzzz_xyyz_1, tg_yzzz_xyzz_0, tg_yzzz_xyzz_1, tg_yzzz_xzzz_0, \
                                         tg_yzzz_xzzz_1, tg_yzzz_yyyy_0, tg_yzzz_yyyy_1, tg_yzzz_yyyz_0, tg_yzzz_yyyz_1, \
                                         tg_yzzz_yyzz_0, tg_yzzz_yyzz_1, tg_yzzz_yzzz_0, tg_yzzz_yzzz_1, wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_yyyyyy_yzzz_0[j] = pb_y * tg_yyyyy_yzzz_0[j] + wp_y[j] * tg_yyyyy_yzzz_1[j] + 2.5 * fl1_fx * tg_yyyy_yzzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyy_yzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyy_zzz_1[j];

                    tg_yyyyyy_zzzz_0[j] = pb_y * tg_yyyyy_zzzz_0[j] + wp_y[j] * tg_yyyyy_zzzz_1[j] + 2.5 * fl1_fx * tg_yyyy_zzzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyy_zzzz_1[j];

                    tg_yyyyyz_xxxx_0[j] = pb_y * tg_yyyyz_xxxx_0[j] + wp_y[j] * tg_yyyyz_xxxx_1[j] + 2.0 * fl1_fx * tg_yyyz_xxxx_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyz_xxxx_1[j];

                    tg_yyyyyz_xxxy_0[j] = pb_y * tg_yyyyz_xxxy_0[j] + wp_y[j] * tg_yyyyz_xxxy_1[j] + 2.0 * fl1_fx * tg_yyyz_xxxy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyz_xxxy_1[j] + 0.5 * fl1_fxn * tg_yyyyz_xxx_1[j];

                    tg_yyyyyz_xxxz_0[j] = pb_y * tg_yyyyz_xxxz_0[j] + wp_y[j] * tg_yyyyz_xxxz_1[j] + 2.0 * fl1_fx * tg_yyyz_xxxz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyz_xxxz_1[j];

                    tg_yyyyyz_xxyy_0[j] = pb_y * tg_yyyyz_xxyy_0[j] + wp_y[j] * tg_yyyyz_xxyy_1[j] + 2.0 * fl1_fx * tg_yyyz_xxyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyz_xxyy_1[j] + fl1_fxn * tg_yyyyz_xxy_1[j];

                    tg_yyyyyz_xxyz_0[j] = pb_y * tg_yyyyz_xxyz_0[j] + wp_y[j] * tg_yyyyz_xxyz_1[j] + 2.0 * fl1_fx * tg_yyyz_xxyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyz_xxyz_1[j] + 0.5 * fl1_fxn * tg_yyyyz_xxz_1[j];

                    tg_yyyyyz_xxzz_0[j] = pb_y * tg_yyyyz_xxzz_0[j] + wp_y[j] * tg_yyyyz_xxzz_1[j] + 2.0 * fl1_fx * tg_yyyz_xxzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyz_xxzz_1[j];

                    tg_yyyyyz_xyyy_0[j] = pb_y * tg_yyyyz_xyyy_0[j] + wp_y[j] * tg_yyyyz_xyyy_1[j] + 2.0 * fl1_fx * tg_yyyz_xyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyz_xyyy_1[j] + 1.5 * fl1_fxn * tg_yyyyz_xyy_1[j];

                    tg_yyyyyz_xyyz_0[j] = pb_y * tg_yyyyz_xyyz_0[j] + wp_y[j] * tg_yyyyz_xyyz_1[j] + 2.0 * fl1_fx * tg_yyyz_xyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyz_xyyz_1[j] + fl1_fxn * tg_yyyyz_xyz_1[j];

                    tg_yyyyyz_xyzz_0[j] = pb_y * tg_yyyyz_xyzz_0[j] + wp_y[j] * tg_yyyyz_xyzz_1[j] + 2.0 * fl1_fx * tg_yyyz_xyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyyyz_xzz_1[j];

                    tg_yyyyyz_xzzz_0[j] = pb_y * tg_yyyyz_xzzz_0[j] + wp_y[j] * tg_yyyyz_xzzz_1[j] + 2.0 * fl1_fx * tg_yyyz_xzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyz_xzzz_1[j];

                    tg_yyyyyz_yyyy_0[j] = pb_y * tg_yyyyz_yyyy_0[j] + wp_y[j] * tg_yyyyz_yyyy_1[j] + 2.0 * fl1_fx * tg_yyyz_yyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyz_yyyy_1[j] + 2.0 * fl1_fxn * tg_yyyyz_yyy_1[j];

                    tg_yyyyyz_yyyz_0[j] = pb_y * tg_yyyyz_yyyz_0[j] + wp_y[j] * tg_yyyyz_yyyz_1[j] + 2.0 * fl1_fx * tg_yyyz_yyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyz_yyyz_1[j] + 1.5 * fl1_fxn * tg_yyyyz_yyz_1[j];

                    tg_yyyyyz_yyzz_0[j] = pb_y * tg_yyyyz_yyzz_0[j] + wp_y[j] * tg_yyyyz_yyzz_1[j] + 2.0 * fl1_fx * tg_yyyz_yyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyz_yyzz_1[j] + fl1_fxn * tg_yyyyz_yzz_1[j];

                    tg_yyyyyz_yzzz_0[j] = pb_y * tg_yyyyz_yzzz_0[j] + wp_y[j] * tg_yyyyz_yzzz_1[j] + 2.0 * fl1_fx * tg_yyyz_yzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyz_yzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyz_zzz_1[j];

                    tg_yyyyyz_zzzz_0[j] = pb_y * tg_yyyyz_zzzz_0[j] + wp_y[j] * tg_yyyyz_zzzz_1[j] + 2.0 * fl1_fx * tg_yyyz_zzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyz_zzzz_1[j];

                    tg_yyyyzz_xxxx_0[j] = pb_y * tg_yyyzz_xxxx_0[j] + wp_y[j] * tg_yyyzz_xxxx_1[j] + 1.5 * fl1_fx * tg_yyzz_xxxx_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzz_xxxx_1[j];

                    tg_yyyyzz_xxxy_0[j] = pb_y * tg_yyyzz_xxxy_0[j] + wp_y[j] * tg_yyyzz_xxxy_1[j] + 1.5 * fl1_fx * tg_yyzz_xxxy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzz_xxxy_1[j] + 0.5 * fl1_fxn * tg_yyyzz_xxx_1[j];

                    tg_yyyyzz_xxxz_0[j] = pb_y * tg_yyyzz_xxxz_0[j] + wp_y[j] * tg_yyyzz_xxxz_1[j] + 1.5 * fl1_fx * tg_yyzz_xxxz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzz_xxxz_1[j];

                    tg_yyyyzz_xxyy_0[j] = pb_y * tg_yyyzz_xxyy_0[j] + wp_y[j] * tg_yyyzz_xxyy_1[j] + 1.5 * fl1_fx * tg_yyzz_xxyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzz_xxyy_1[j] + fl1_fxn * tg_yyyzz_xxy_1[j];

                    tg_yyyyzz_xxyz_0[j] = pb_y * tg_yyyzz_xxyz_0[j] + wp_y[j] * tg_yyyzz_xxyz_1[j] + 1.5 * fl1_fx * tg_yyzz_xxyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzz_xxyz_1[j] + 0.5 * fl1_fxn * tg_yyyzz_xxz_1[j];

                    tg_yyyyzz_xxzz_0[j] = pb_y * tg_yyyzz_xxzz_0[j] + wp_y[j] * tg_yyyzz_xxzz_1[j] + 1.5 * fl1_fx * tg_yyzz_xxzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzz_xxzz_1[j];

                    tg_yyyyzz_xyyy_0[j] = pb_y * tg_yyyzz_xyyy_0[j] + wp_y[j] * tg_yyyzz_xyyy_1[j] + 1.5 * fl1_fx * tg_yyzz_xyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzz_xyyy_1[j] + 1.5 * fl1_fxn * tg_yyyzz_xyy_1[j];

                    tg_yyyyzz_xyyz_0[j] = pb_y * tg_yyyzz_xyyz_0[j] + wp_y[j] * tg_yyyzz_xyyz_1[j] + 1.5 * fl1_fx * tg_yyzz_xyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzz_xyyz_1[j] + fl1_fxn * tg_yyyzz_xyz_1[j];

                    tg_yyyyzz_xyzz_0[j] = pb_y * tg_yyyzz_xyzz_0[j] + wp_y[j] * tg_yyyzz_xyzz_1[j] + 1.5 * fl1_fx * tg_yyzz_xyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyyzz_xzz_1[j];

                    tg_yyyyzz_xzzz_0[j] = pb_y * tg_yyyzz_xzzz_0[j] + wp_y[j] * tg_yyyzz_xzzz_1[j] + 1.5 * fl1_fx * tg_yyzz_xzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzz_xzzz_1[j];

                    tg_yyyyzz_yyyy_0[j] = pb_y * tg_yyyzz_yyyy_0[j] + wp_y[j] * tg_yyyzz_yyyy_1[j] + 1.5 * fl1_fx * tg_yyzz_yyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzz_yyyy_1[j] + 2.0 * fl1_fxn * tg_yyyzz_yyy_1[j];

                    tg_yyyyzz_yyyz_0[j] = pb_y * tg_yyyzz_yyyz_0[j] + wp_y[j] * tg_yyyzz_yyyz_1[j] + 1.5 * fl1_fx * tg_yyzz_yyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzz_yyyz_1[j] + 1.5 * fl1_fxn * tg_yyyzz_yyz_1[j];

                    tg_yyyyzz_yyzz_0[j] = pb_y * tg_yyyzz_yyzz_0[j] + wp_y[j] * tg_yyyzz_yyzz_1[j] + 1.5 * fl1_fx * tg_yyzz_yyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzz_yyzz_1[j] + fl1_fxn * tg_yyyzz_yzz_1[j];

                    tg_yyyyzz_yzzz_0[j] = pb_y * tg_yyyzz_yzzz_0[j] + wp_y[j] * tg_yyyzz_yzzz_1[j] + 1.5 * fl1_fx * tg_yyzz_yzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzz_yzzz_1[j] + 0.5 * fl1_fxn * tg_yyyzz_zzz_1[j];

                    tg_yyyyzz_zzzz_0[j] = pb_y * tg_yyyzz_zzzz_0[j] + wp_y[j] * tg_yyyzz_zzzz_1[j] + 1.5 * fl1_fx * tg_yyzz_zzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzz_zzzz_1[j];

                    tg_yyyzzz_xxxx_0[j] = pb_y * tg_yyzzz_xxxx_0[j] + wp_y[j] * tg_yyzzz_xxxx_1[j] + fl1_fx * tg_yzzz_xxxx_0[j] - fl1_fx * fl1_fza * tg_yzzz_xxxx_1[j];

                    tg_yyyzzz_xxxy_0[j] = pb_y * tg_yyzzz_xxxy_0[j] + wp_y[j] * tg_yyzzz_xxxy_1[j] + fl1_fx * tg_yzzz_xxxy_0[j] - fl1_fx * fl1_fza * tg_yzzz_xxxy_1[j] + 0.5 * fl1_fxn * tg_yyzzz_xxx_1[j];

                    tg_yyyzzz_xxxz_0[j] = pb_y * tg_yyzzz_xxxz_0[j] + wp_y[j] * tg_yyzzz_xxxz_1[j] + fl1_fx * tg_yzzz_xxxz_0[j] - fl1_fx * fl1_fza * tg_yzzz_xxxz_1[j];

                    tg_yyyzzz_xxyy_0[j] = pb_y * tg_yyzzz_xxyy_0[j] + wp_y[j] * tg_yyzzz_xxyy_1[j] + fl1_fx * tg_yzzz_xxyy_0[j] - fl1_fx * fl1_fza * tg_yzzz_xxyy_1[j] + fl1_fxn * tg_yyzzz_xxy_1[j];

                    tg_yyyzzz_xxyz_0[j] = pb_y * tg_yyzzz_xxyz_0[j] + wp_y[j] * tg_yyzzz_xxyz_1[j] + fl1_fx * tg_yzzz_xxyz_0[j] - fl1_fx * fl1_fza * tg_yzzz_xxyz_1[j] + 0.5 * fl1_fxn * tg_yyzzz_xxz_1[j];

                    tg_yyyzzz_xxzz_0[j] = pb_y * tg_yyzzz_xxzz_0[j] + wp_y[j] * tg_yyzzz_xxzz_1[j] + fl1_fx * tg_yzzz_xxzz_0[j] - fl1_fx * fl1_fza * tg_yzzz_xxzz_1[j];

                    tg_yyyzzz_xyyy_0[j] = pb_y * tg_yyzzz_xyyy_0[j] + wp_y[j] * tg_yyzzz_xyyy_1[j] + fl1_fx * tg_yzzz_xyyy_0[j] - fl1_fx * fl1_fza * tg_yzzz_xyyy_1[j] + 1.5 * fl1_fxn * tg_yyzzz_xyy_1[j];

                    tg_yyyzzz_xyyz_0[j] = pb_y * tg_yyzzz_xyyz_0[j] + wp_y[j] * tg_yyzzz_xyyz_1[j] + fl1_fx * tg_yzzz_xyyz_0[j] - fl1_fx * fl1_fza * tg_yzzz_xyyz_1[j] + fl1_fxn * tg_yyzzz_xyz_1[j];

                    tg_yyyzzz_xyzz_0[j] = pb_y * tg_yyzzz_xyzz_0[j] + wp_y[j] * tg_yyzzz_xyzz_1[j] + fl1_fx * tg_yzzz_xyzz_0[j] - fl1_fx * fl1_fza * tg_yzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyzzz_xzz_1[j];

                    tg_yyyzzz_xzzz_0[j] = pb_y * tg_yyzzz_xzzz_0[j] + wp_y[j] * tg_yyzzz_xzzz_1[j] + fl1_fx * tg_yzzz_xzzz_0[j] - fl1_fx * fl1_fza * tg_yzzz_xzzz_1[j];

                    tg_yyyzzz_yyyy_0[j] = pb_y * tg_yyzzz_yyyy_0[j] + wp_y[j] * tg_yyzzz_yyyy_1[j] + fl1_fx * tg_yzzz_yyyy_0[j] - fl1_fx * fl1_fza * tg_yzzz_yyyy_1[j] + 2.0 * fl1_fxn * tg_yyzzz_yyy_1[j];

                    tg_yyyzzz_yyyz_0[j] = pb_y * tg_yyzzz_yyyz_0[j] + wp_y[j] * tg_yyzzz_yyyz_1[j] + fl1_fx * tg_yzzz_yyyz_0[j] - fl1_fx * fl1_fza * tg_yzzz_yyyz_1[j] + 1.5 * fl1_fxn * tg_yyzzz_yyz_1[j];

                    tg_yyyzzz_yyzz_0[j] = pb_y * tg_yyzzz_yyzz_0[j] + wp_y[j] * tg_yyzzz_yyzz_1[j] + fl1_fx * tg_yzzz_yyzz_0[j] - fl1_fx * fl1_fza * tg_yzzz_yyzz_1[j] + fl1_fxn * tg_yyzzz_yzz_1[j];

                    tg_yyyzzz_yzzz_0[j] = pb_y * tg_yyzzz_yzzz_0[j] + wp_y[j] * tg_yyzzz_yzzz_1[j] + fl1_fx * tg_yzzz_yzzz_0[j] - fl1_fx * fl1_fza * tg_yzzz_yzzz_1[j] + 0.5 * fl1_fxn * tg_yyzzz_zzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSISG_374_420(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (374,420)

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
                                             {6, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_6_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_6_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_5_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_5_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_4_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_5_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_yyzzz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 284); 

                auto tg_yzzzz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 285); 

                auto tg_yzzzz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 286); 

                auto tg_yzzzz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 287); 

                auto tg_yzzzz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 288); 

                auto tg_yzzzz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 289); 

                auto tg_yzzzz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 290); 

                auto tg_yzzzz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 291); 

                auto tg_yzzzz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 292); 

                auto tg_yzzzz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 293); 

                auto tg_yzzzz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 294); 

                auto tg_yzzzz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 295); 

                auto tg_yzzzz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 296); 

                auto tg_yzzzz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 297); 

                auto tg_yzzzz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 298); 

                auto tg_yzzzz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 299); 

                auto tg_zzzzz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 300); 

                auto tg_zzzzz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 301); 

                auto tg_zzzzz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 302); 

                auto tg_zzzzz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 303); 

                auto tg_zzzzz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 304); 

                auto tg_zzzzz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 305); 

                auto tg_zzzzz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 306); 

                auto tg_zzzzz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 307); 

                auto tg_zzzzz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 308); 

                auto tg_zzzzz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 309); 

                auto tg_zzzzz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 310); 

                auto tg_zzzzz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 311); 

                auto tg_zzzzz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 312); 

                auto tg_zzzzz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 313); 

                auto tg_zzzzz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 314); 

                auto tg_yyzzz_zzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 284); 

                auto tg_yzzzz_xxxx_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 285); 

                auto tg_yzzzz_xxxy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 286); 

                auto tg_yzzzz_xxxz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 287); 

                auto tg_yzzzz_xxyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 288); 

                auto tg_yzzzz_xxyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 289); 

                auto tg_yzzzz_xxzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 290); 

                auto tg_yzzzz_xyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 291); 

                auto tg_yzzzz_xyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 292); 

                auto tg_yzzzz_xyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 293); 

                auto tg_yzzzz_xzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 294); 

                auto tg_yzzzz_yyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 295); 

                auto tg_yzzzz_yyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 296); 

                auto tg_yzzzz_yyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 297); 

                auto tg_yzzzz_yzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 298); 

                auto tg_yzzzz_zzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 299); 

                auto tg_zzzzz_xxxx_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 300); 

                auto tg_zzzzz_xxxy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 301); 

                auto tg_zzzzz_xxxz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 302); 

                auto tg_zzzzz_xxyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 303); 

                auto tg_zzzzz_xxyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 304); 

                auto tg_zzzzz_xxzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 305); 

                auto tg_zzzzz_xyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 306); 

                auto tg_zzzzz_xyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 307); 

                auto tg_zzzzz_xyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 308); 

                auto tg_zzzzz_xzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 309); 

                auto tg_zzzzz_yyyy_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 310); 

                auto tg_zzzzz_yyyz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 311); 

                auto tg_zzzzz_yyzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 312); 

                auto tg_zzzzz_yzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 313); 

                auto tg_zzzzz_zzzz_1 = primBuffer.data(pidx_g_5_4_m1 + 315 * idx + 314); 

                auto tg_yzzz_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 209); 

                auto tg_zzzz_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 210); 

                auto tg_zzzz_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 211); 

                auto tg_zzzz_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 212); 

                auto tg_zzzz_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 213); 

                auto tg_zzzz_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 214); 

                auto tg_zzzz_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 215); 

                auto tg_zzzz_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 216); 

                auto tg_zzzz_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 217); 

                auto tg_zzzz_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 218); 

                auto tg_zzzz_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 219); 

                auto tg_zzzz_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 220); 

                auto tg_zzzz_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 221); 

                auto tg_zzzz_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 222); 

                auto tg_zzzz_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 223); 

                auto tg_zzzz_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 224); 

                auto tg_yzzz_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 209); 

                auto tg_zzzz_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 210); 

                auto tg_zzzz_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 211); 

                auto tg_zzzz_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 212); 

                auto tg_zzzz_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 213); 

                auto tg_zzzz_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 214); 

                auto tg_zzzz_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 215); 

                auto tg_zzzz_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 216); 

                auto tg_zzzz_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 217); 

                auto tg_zzzz_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 218); 

                auto tg_zzzz_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 219); 

                auto tg_zzzz_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 220); 

                auto tg_zzzz_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 221); 

                auto tg_zzzz_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 222); 

                auto tg_zzzz_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 223); 

                auto tg_zzzz_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 224); 

                auto tg_yzzzz_xxx_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 190); 

                auto tg_yzzzz_xxy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 191); 

                auto tg_yzzzz_xxz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 192); 

                auto tg_yzzzz_xyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 193); 

                auto tg_yzzzz_xyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 194); 

                auto tg_yzzzz_xzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 195); 

                auto tg_yzzzz_yyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 196); 

                auto tg_yzzzz_yyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 197); 

                auto tg_yzzzz_yzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 198); 

                auto tg_yzzzz_zzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 199); 

                auto tg_zzzzz_xxx_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 200); 

                auto tg_zzzzz_xxy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 201); 

                auto tg_zzzzz_xxz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 202); 

                auto tg_zzzzz_xyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 203); 

                auto tg_zzzzz_xyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 204); 

                auto tg_zzzzz_xzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 205); 

                auto tg_zzzzz_yyy_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 206); 

                auto tg_zzzzz_yyz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 207); 

                auto tg_zzzzz_yzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 208); 

                auto tg_zzzzz_zzz_1 = primBuffer.data(pidx_g_5_3_m1 + 210 * idx + 209); 

                // set up pointers to integrals

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

                // Batch of Integrals (374,420)

                #pragma omp simd aligned(fxn, fza, tg_yyyzzz_zzzz_0, tg_yyzzz_zzzz_0, tg_yyzzz_zzzz_1, \
                                         tg_yyzzzz_xxxx_0, tg_yyzzzz_xxxy_0, tg_yyzzzz_xxxz_0, tg_yyzzzz_xxyy_0, \
                                         tg_yyzzzz_xxyz_0, tg_yyzzzz_xxzz_0, tg_yyzzzz_xyyy_0, tg_yyzzzz_xyyz_0, \
                                         tg_yyzzzz_xyzz_0, tg_yyzzzz_xzzz_0, tg_yyzzzz_yyyy_0, tg_yyzzzz_yyyz_0, \
                                         tg_yyzzzz_yyzz_0, tg_yyzzzz_yzzz_0, tg_yyzzzz_zzzz_0, tg_yzzz_zzzz_0, tg_yzzz_zzzz_1, \
                                         tg_yzzzz_xxx_1, tg_yzzzz_xxxx_0, tg_yzzzz_xxxx_1, tg_yzzzz_xxxy_0, tg_yzzzz_xxxy_1, \
                                         tg_yzzzz_xxxz_0, tg_yzzzz_xxxz_1, tg_yzzzz_xxy_1, tg_yzzzz_xxyy_0, tg_yzzzz_xxyy_1, \
                                         tg_yzzzz_xxyz_0, tg_yzzzz_xxyz_1, tg_yzzzz_xxz_1, tg_yzzzz_xxzz_0, tg_yzzzz_xxzz_1, \
                                         tg_yzzzz_xyy_1, tg_yzzzz_xyyy_0, tg_yzzzz_xyyy_1, tg_yzzzz_xyyz_0, tg_yzzzz_xyyz_1, \
                                         tg_yzzzz_xyz_1, tg_yzzzz_xyzz_0, tg_yzzzz_xyzz_1, tg_yzzzz_xzz_1, tg_yzzzz_xzzz_0, \
                                         tg_yzzzz_xzzz_1, tg_yzzzz_yyy_1, tg_yzzzz_yyyy_0, tg_yzzzz_yyyy_1, tg_yzzzz_yyyz_0, \
                                         tg_yzzzz_yyyz_1, tg_yzzzz_yyz_1, tg_yzzzz_yyzz_0, tg_yzzzz_yyzz_1, tg_yzzzz_yzz_1, \
                                         tg_yzzzz_yzzz_0, tg_yzzzz_yzzz_1, tg_yzzzz_zzz_1, tg_yzzzz_zzzz_0, tg_yzzzz_zzzz_1, \
                                         tg_yzzzzz_xxxx_0, tg_yzzzzz_xxxy_0, tg_yzzzzz_xxxz_0, tg_yzzzzz_xxyy_0, \
                                         tg_yzzzzz_xxyz_0, tg_yzzzzz_xxzz_0, tg_yzzzzz_xyyy_0, tg_yzzzzz_xyyz_0, \
                                         tg_yzzzzz_xyzz_0, tg_yzzzzz_xzzz_0, tg_yzzzzz_yyyy_0, tg_yzzzzz_yyyz_0, \
                                         tg_yzzzzz_yyzz_0, tg_yzzzzz_yzzz_0, tg_yzzzzz_zzzz_0, tg_zzzz_xxxx_0, tg_zzzz_xxxx_1, \
                                         tg_zzzz_xxxy_0, tg_zzzz_xxxy_1, tg_zzzz_xxxz_0, tg_zzzz_xxxz_1, tg_zzzz_xxyy_0, \
                                         tg_zzzz_xxyy_1, tg_zzzz_xxyz_0, tg_zzzz_xxyz_1, tg_zzzz_xxzz_0, tg_zzzz_xxzz_1, \
                                         tg_zzzz_xyyy_0, tg_zzzz_xyyy_1, tg_zzzz_xyyz_0, tg_zzzz_xyyz_1, tg_zzzz_xyzz_0, \
                                         tg_zzzz_xyzz_1, tg_zzzz_xzzz_0, tg_zzzz_xzzz_1, tg_zzzz_yyyy_0, tg_zzzz_yyyy_1, \
                                         tg_zzzz_yyyz_0, tg_zzzz_yyyz_1, tg_zzzz_yyzz_0, tg_zzzz_yyzz_1, tg_zzzz_yzzz_0, \
                                         tg_zzzz_yzzz_1, tg_zzzz_zzzz_0, tg_zzzz_zzzz_1, tg_zzzzz_xxx_1, tg_zzzzz_xxxx_0, \
                                         tg_zzzzz_xxxx_1, tg_zzzzz_xxxy_0, tg_zzzzz_xxxy_1, tg_zzzzz_xxxz_0, tg_zzzzz_xxxz_1, \
                                         tg_zzzzz_xxy_1, tg_zzzzz_xxyy_0, tg_zzzzz_xxyy_1, tg_zzzzz_xxyz_0, tg_zzzzz_xxyz_1, \
                                         tg_zzzzz_xxz_1, tg_zzzzz_xxzz_0, tg_zzzzz_xxzz_1, tg_zzzzz_xyy_1, tg_zzzzz_xyyy_0, \
                                         tg_zzzzz_xyyy_1, tg_zzzzz_xyyz_0, tg_zzzzz_xyyz_1, tg_zzzzz_xyz_1, tg_zzzzz_xyzz_0, \
                                         tg_zzzzz_xyzz_1, tg_zzzzz_xzz_1, tg_zzzzz_xzzz_0, tg_zzzzz_xzzz_1, tg_zzzzz_yyy_1, \
                                         tg_zzzzz_yyyy_0, tg_zzzzz_yyyy_1, tg_zzzzz_yyyz_0, tg_zzzzz_yyyz_1, tg_zzzzz_yyz_1, \
                                         tg_zzzzz_yyzz_0, tg_zzzzz_yyzz_1, tg_zzzzz_yzz_1, tg_zzzzz_yzzz_0, tg_zzzzz_yzzz_1, \
                                         tg_zzzzz_zzz_1, tg_zzzzz_zzzz_0, tg_zzzzz_zzzz_1, tg_zzzzzz_xxxx_0, \
                                         tg_zzzzzz_xxxy_0, tg_zzzzzz_xxxz_0, tg_zzzzzz_xxyy_0, tg_zzzzzz_xxyz_0, \
                                         tg_zzzzzz_xxzz_0, tg_zzzzzz_xyyy_0, tg_zzzzzz_xyyz_0, tg_zzzzzz_xyzz_0, \
                                         tg_zzzzzz_xzzz_0, tg_zzzzzz_yyyy_0, tg_zzzzzz_yyyz_0, tg_zzzzzz_yyzz_0, \
                                         tg_zzzzzz_yzzz_0, tg_zzzzzz_zzzz_0, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_yyyzzz_zzzz_0[j] = pb_y * tg_yyzzz_zzzz_0[j] + wp_y[j] * tg_yyzzz_zzzz_1[j] + fl1_fx * tg_yzzz_zzzz_0[j] - fl1_fx * fl1_fza * tg_yzzz_zzzz_1[j];

                    tg_yyzzzz_xxxx_0[j] = pb_y * tg_yzzzz_xxxx_0[j] + wp_y[j] * tg_yzzzz_xxxx_1[j] + 0.5 * fl1_fx * tg_zzzz_xxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_xxxx_1[j];

                    tg_yyzzzz_xxxy_0[j] = pb_y * tg_yzzzz_xxxy_0[j] + wp_y[j] * tg_yzzzz_xxxy_1[j] + 0.5 * fl1_fx * tg_zzzz_xxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_xxxy_1[j] + 0.5 * fl1_fxn * tg_yzzzz_xxx_1[j];

                    tg_yyzzzz_xxxz_0[j] = pb_y * tg_yzzzz_xxxz_0[j] + wp_y[j] * tg_yzzzz_xxxz_1[j] + 0.5 * fl1_fx * tg_zzzz_xxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_xxxz_1[j];

                    tg_yyzzzz_xxyy_0[j] = pb_y * tg_yzzzz_xxyy_0[j] + wp_y[j] * tg_yzzzz_xxyy_1[j] + 0.5 * fl1_fx * tg_zzzz_xxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_xxyy_1[j] + fl1_fxn * tg_yzzzz_xxy_1[j];

                    tg_yyzzzz_xxyz_0[j] = pb_y * tg_yzzzz_xxyz_0[j] + wp_y[j] * tg_yzzzz_xxyz_1[j] + 0.5 * fl1_fx * tg_zzzz_xxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_xxyz_1[j] + 0.5 * fl1_fxn * tg_yzzzz_xxz_1[j];

                    tg_yyzzzz_xxzz_0[j] = pb_y * tg_yzzzz_xxzz_0[j] + wp_y[j] * tg_yzzzz_xxzz_1[j] + 0.5 * fl1_fx * tg_zzzz_xxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_xxzz_1[j];

                    tg_yyzzzz_xyyy_0[j] = pb_y * tg_yzzzz_xyyy_0[j] + wp_y[j] * tg_yzzzz_xyyy_1[j] + 0.5 * fl1_fx * tg_zzzz_xyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_xyyy_1[j] + 1.5 * fl1_fxn * tg_yzzzz_xyy_1[j];

                    tg_yyzzzz_xyyz_0[j] = pb_y * tg_yzzzz_xyyz_0[j] + wp_y[j] * tg_yzzzz_xyyz_1[j] + 0.5 * fl1_fx * tg_zzzz_xyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_xyyz_1[j] + fl1_fxn * tg_yzzzz_xyz_1[j];

                    tg_yyzzzz_xyzz_0[j] = pb_y * tg_yzzzz_xyzz_0[j] + wp_y[j] * tg_yzzzz_xyzz_1[j] + 0.5 * fl1_fx * tg_zzzz_xyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yzzzz_xzz_1[j];

                    tg_yyzzzz_xzzz_0[j] = pb_y * tg_yzzzz_xzzz_0[j] + wp_y[j] * tg_yzzzz_xzzz_1[j] + 0.5 * fl1_fx * tg_zzzz_xzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_xzzz_1[j];

                    tg_yyzzzz_yyyy_0[j] = pb_y * tg_yzzzz_yyyy_0[j] + wp_y[j] * tg_yzzzz_yyyy_1[j] + 0.5 * fl1_fx * tg_zzzz_yyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_yyyy_1[j] + 2.0 * fl1_fxn * tg_yzzzz_yyy_1[j];

                    tg_yyzzzz_yyyz_0[j] = pb_y * tg_yzzzz_yyyz_0[j] + wp_y[j] * tg_yzzzz_yyyz_1[j] + 0.5 * fl1_fx * tg_zzzz_yyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_yyyz_1[j] + 1.5 * fl1_fxn * tg_yzzzz_yyz_1[j];

                    tg_yyzzzz_yyzz_0[j] = pb_y * tg_yzzzz_yyzz_0[j] + wp_y[j] * tg_yzzzz_yyzz_1[j] + 0.5 * fl1_fx * tg_zzzz_yyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_yyzz_1[j] + fl1_fxn * tg_yzzzz_yzz_1[j];

                    tg_yyzzzz_yzzz_0[j] = pb_y * tg_yzzzz_yzzz_0[j] + wp_y[j] * tg_yzzzz_yzzz_1[j] + 0.5 * fl1_fx * tg_zzzz_yzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_yzzz_1[j] + 0.5 * fl1_fxn * tg_yzzzz_zzz_1[j];

                    tg_yyzzzz_zzzz_0[j] = pb_y * tg_yzzzz_zzzz_0[j] + wp_y[j] * tg_yzzzz_zzzz_1[j] + 0.5 * fl1_fx * tg_zzzz_zzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_zzzz_1[j];

                    tg_yzzzzz_xxxx_0[j] = pb_y * tg_zzzzz_xxxx_0[j] + wp_y[j] * tg_zzzzz_xxxx_1[j];

                    tg_yzzzzz_xxxy_0[j] = pb_y * tg_zzzzz_xxxy_0[j] + wp_y[j] * tg_zzzzz_xxxy_1[j] + 0.5 * fl1_fxn * tg_zzzzz_xxx_1[j];

                    tg_yzzzzz_xxxz_0[j] = pb_y * tg_zzzzz_xxxz_0[j] + wp_y[j] * tg_zzzzz_xxxz_1[j];

                    tg_yzzzzz_xxyy_0[j] = pb_y * tg_zzzzz_xxyy_0[j] + wp_y[j] * tg_zzzzz_xxyy_1[j] + fl1_fxn * tg_zzzzz_xxy_1[j];

                    tg_yzzzzz_xxyz_0[j] = pb_y * tg_zzzzz_xxyz_0[j] + wp_y[j] * tg_zzzzz_xxyz_1[j] + 0.5 * fl1_fxn * tg_zzzzz_xxz_1[j];

                    tg_yzzzzz_xxzz_0[j] = pb_y * tg_zzzzz_xxzz_0[j] + wp_y[j] * tg_zzzzz_xxzz_1[j];

                    tg_yzzzzz_xyyy_0[j] = pb_y * tg_zzzzz_xyyy_0[j] + wp_y[j] * tg_zzzzz_xyyy_1[j] + 1.5 * fl1_fxn * tg_zzzzz_xyy_1[j];

                    tg_yzzzzz_xyyz_0[j] = pb_y * tg_zzzzz_xyyz_0[j] + wp_y[j] * tg_zzzzz_xyyz_1[j] + fl1_fxn * tg_zzzzz_xyz_1[j];

                    tg_yzzzzz_xyzz_0[j] = pb_y * tg_zzzzz_xyzz_0[j] + wp_y[j] * tg_zzzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_zzzzz_xzz_1[j];

                    tg_yzzzzz_xzzz_0[j] = pb_y * tg_zzzzz_xzzz_0[j] + wp_y[j] * tg_zzzzz_xzzz_1[j];

                    tg_yzzzzz_yyyy_0[j] = pb_y * tg_zzzzz_yyyy_0[j] + wp_y[j] * tg_zzzzz_yyyy_1[j] + 2.0 * fl1_fxn * tg_zzzzz_yyy_1[j];

                    tg_yzzzzz_yyyz_0[j] = pb_y * tg_zzzzz_yyyz_0[j] + wp_y[j] * tg_zzzzz_yyyz_1[j] + 1.5 * fl1_fxn * tg_zzzzz_yyz_1[j];

                    tg_yzzzzz_yyzz_0[j] = pb_y * tg_zzzzz_yyzz_0[j] + wp_y[j] * tg_zzzzz_yyzz_1[j] + fl1_fxn * tg_zzzzz_yzz_1[j];

                    tg_yzzzzz_yzzz_0[j] = pb_y * tg_zzzzz_yzzz_0[j] + wp_y[j] * tg_zzzzz_yzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzz_zzz_1[j];

                    tg_yzzzzz_zzzz_0[j] = pb_y * tg_zzzzz_zzzz_0[j] + wp_y[j] * tg_zzzzz_zzzz_1[j];

                    tg_zzzzzz_xxxx_0[j] = pb_z * tg_zzzzz_xxxx_0[j] + wp_z[j] * tg_zzzzz_xxxx_1[j] + 2.5 * fl1_fx * tg_zzzz_xxxx_0[j] - 2.5 * fl1_fx * fl1_fza * tg_zzzz_xxxx_1[j];

                    tg_zzzzzz_xxxy_0[j] = pb_z * tg_zzzzz_xxxy_0[j] + wp_z[j] * tg_zzzzz_xxxy_1[j] + 2.5 * fl1_fx * tg_zzzz_xxxy_0[j] - 2.5 * fl1_fx * fl1_fza * tg_zzzz_xxxy_1[j];

                    tg_zzzzzz_xxxz_0[j] = pb_z * tg_zzzzz_xxxz_0[j] + wp_z[j] * tg_zzzzz_xxxz_1[j] + 2.5 * fl1_fx * tg_zzzz_xxxz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_zzzz_xxxz_1[j] + 0.5 * fl1_fxn * tg_zzzzz_xxx_1[j];

                    tg_zzzzzz_xxyy_0[j] = pb_z * tg_zzzzz_xxyy_0[j] + wp_z[j] * tg_zzzzz_xxyy_1[j] + 2.5 * fl1_fx * tg_zzzz_xxyy_0[j] - 2.5 * fl1_fx * fl1_fza * tg_zzzz_xxyy_1[j];

                    tg_zzzzzz_xxyz_0[j] = pb_z * tg_zzzzz_xxyz_0[j] + wp_z[j] * tg_zzzzz_xxyz_1[j] + 2.5 * fl1_fx * tg_zzzz_xxyz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_zzzz_xxyz_1[j] + 0.5 * fl1_fxn * tg_zzzzz_xxy_1[j];

                    tg_zzzzzz_xxzz_0[j] = pb_z * tg_zzzzz_xxzz_0[j] + wp_z[j] * tg_zzzzz_xxzz_1[j] + 2.5 * fl1_fx * tg_zzzz_xxzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_zzzz_xxzz_1[j] + fl1_fxn * tg_zzzzz_xxz_1[j];

                    tg_zzzzzz_xyyy_0[j] = pb_z * tg_zzzzz_xyyy_0[j] + wp_z[j] * tg_zzzzz_xyyy_1[j] + 2.5 * fl1_fx * tg_zzzz_xyyy_0[j] - 2.5 * fl1_fx * fl1_fza * tg_zzzz_xyyy_1[j];

                    tg_zzzzzz_xyyz_0[j] = pb_z * tg_zzzzz_xyyz_0[j] + wp_z[j] * tg_zzzzz_xyyz_1[j] + 2.5 * fl1_fx * tg_zzzz_xyyz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_zzzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_zzzzz_xyy_1[j];

                    tg_zzzzzz_xyzz_0[j] = pb_z * tg_zzzzz_xyzz_0[j] + wp_z[j] * tg_zzzzz_xyzz_1[j] + 2.5 * fl1_fx * tg_zzzz_xyzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_zzzz_xyzz_1[j] + fl1_fxn * tg_zzzzz_xyz_1[j];

                    tg_zzzzzz_xzzz_0[j] = pb_z * tg_zzzzz_xzzz_0[j] + wp_z[j] * tg_zzzzz_xzzz_1[j] + 2.5 * fl1_fx * tg_zzzz_xzzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_zzzz_xzzz_1[j] + 1.5 * fl1_fxn * tg_zzzzz_xzz_1[j];

                    tg_zzzzzz_yyyy_0[j] = pb_z * tg_zzzzz_yyyy_0[j] + wp_z[j] * tg_zzzzz_yyyy_1[j] + 2.5 * fl1_fx * tg_zzzz_yyyy_0[j] - 2.5 * fl1_fx * fl1_fza * tg_zzzz_yyyy_1[j];

                    tg_zzzzzz_yyyz_0[j] = pb_z * tg_zzzzz_yyyz_0[j] + wp_z[j] * tg_zzzzz_yyyz_1[j] + 2.5 * fl1_fx * tg_zzzz_yyyz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_zzzz_yyyz_1[j] + 0.5 * fl1_fxn * tg_zzzzz_yyy_1[j];

                    tg_zzzzzz_yyzz_0[j] = pb_z * tg_zzzzz_yyzz_0[j] + wp_z[j] * tg_zzzzz_yyzz_1[j] + 2.5 * fl1_fx * tg_zzzz_yyzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_zzzz_yyzz_1[j] + fl1_fxn * tg_zzzzz_yyz_1[j];

                    tg_zzzzzz_yzzz_0[j] = pb_z * tg_zzzzz_yzzz_0[j] + wp_z[j] * tg_zzzzz_yzzz_1[j] + 2.5 * fl1_fx * tg_zzzz_yzzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_zzzz_yzzz_1[j] + 1.5 * fl1_fxn * tg_zzzzz_yzz_1[j];

                    tg_zzzzzz_zzzz_0[j] = pb_z * tg_zzzzz_zzzz_0[j] + wp_z[j] * tg_zzzzz_zzzz_1[j] + 2.5 * fl1_fx * tg_zzzz_zzzz_0[j] - 2.5 * fl1_fx * fl1_fza * tg_zzzz_zzzz_1[j] + 2.0 * fl1_fxn * tg_zzzzz_zzz_1[j];
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

