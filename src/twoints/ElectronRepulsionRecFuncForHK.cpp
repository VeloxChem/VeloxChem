//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectronRepulsionRecFuncForHK.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSHSK(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSHSK_0_95(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSHSK_95_190(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 

        erirecfunc::compElectronRepulsionForSHSK_190_285(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSHSK_285_380(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSHSK_380_474(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSHSK_474_568(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSHSK_568_662(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSHSK_662_756(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 
    }

    void
    compElectronRepulsionForSHSK_0_95(      CMemBlock2D<double>& primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& wpDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetPrimPairs,
                                      const int32_t              iContrPair)
    {
        // Batch of Integrals (0,95)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {5, -1, -1, -1},
                                             {7, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_7_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_4_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_3_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_4_6_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {6, -1, -1, -1}, 
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

                auto tg_xxxx_xxxxxxx_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx); 

                auto tg_xxxx_xxxxxxy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 1); 

                auto tg_xxxx_xxxxxxz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 2); 

                auto tg_xxxx_xxxxxyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 3); 

                auto tg_xxxx_xxxxxyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 4); 

                auto tg_xxxx_xxxxxzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 5); 

                auto tg_xxxx_xxxxyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 6); 

                auto tg_xxxx_xxxxyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 7); 

                auto tg_xxxx_xxxxyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 8); 

                auto tg_xxxx_xxxxzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 9); 

                auto tg_xxxx_xxxyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 10); 

                auto tg_xxxx_xxxyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 11); 

                auto tg_xxxx_xxxyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 12); 

                auto tg_xxxx_xxxyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 13); 

                auto tg_xxxx_xxxzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 14); 

                auto tg_xxxx_xxyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 15); 

                auto tg_xxxx_xxyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 16); 

                auto tg_xxxx_xxyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 17); 

                auto tg_xxxx_xxyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 18); 

                auto tg_xxxx_xxyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 19); 

                auto tg_xxxx_xxzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 20); 

                auto tg_xxxx_xyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 21); 

                auto tg_xxxx_xyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 22); 

                auto tg_xxxx_xyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 23); 

                auto tg_xxxx_xyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 24); 

                auto tg_xxxx_xyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 25); 

                auto tg_xxxx_xyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 26); 

                auto tg_xxxx_xzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 27); 

                auto tg_xxxx_yyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 28); 

                auto tg_xxxx_yyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 29); 

                auto tg_xxxx_yyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 30); 

                auto tg_xxxx_yyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 31); 

                auto tg_xxxx_yyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 32); 

                auto tg_xxxx_yyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 33); 

                auto tg_xxxx_yzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 34); 

                auto tg_xxxx_zzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 35); 

                auto tg_xxxy_xxxxxxx_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 36); 

                auto tg_xxxy_xxxxxxy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 37); 

                auto tg_xxxy_xxxxxxz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 38); 

                auto tg_xxxy_xxxxxyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 39); 

                auto tg_xxxy_xxxxxyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 40); 

                auto tg_xxxy_xxxxxzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 41); 

                auto tg_xxxy_xxxxyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 42); 

                auto tg_xxxy_xxxxyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 43); 

                auto tg_xxxy_xxxxyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 44); 

                auto tg_xxxy_xxxxzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 45); 

                auto tg_xxxy_xxxyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 46); 

                auto tg_xxxy_xxxyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 47); 

                auto tg_xxxy_xxxyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 48); 

                auto tg_xxxy_xxxyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 49); 

                auto tg_xxxy_xxxzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 50); 

                auto tg_xxxy_xxyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 51); 

                auto tg_xxxy_xxyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 52); 

                auto tg_xxxy_xxyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 53); 

                auto tg_xxxy_xxyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 54); 

                auto tg_xxxy_xxyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 55); 

                auto tg_xxxy_xxzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 56); 

                auto tg_xxxy_xyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 57); 

                auto tg_xxxy_xyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 58); 

                auto tg_xxxy_xyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 59); 

                auto tg_xxxy_xyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 60); 

                auto tg_xxxy_xyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 61); 

                auto tg_xxxy_xyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 62); 

                auto tg_xxxy_xzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 63); 

                auto tg_xxxy_yyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 64); 

                auto tg_xxxy_yyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 65); 

                auto tg_xxxy_yyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 66); 

                auto tg_xxxy_yyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 67); 

                auto tg_xxxy_yyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 68); 

                auto tg_xxxy_yyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 69); 

                auto tg_xxxy_yzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 70); 

                auto tg_xxxy_zzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 71); 

                auto tg_xxxz_xxxxxxx_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 72); 

                auto tg_xxxz_xxxxxxy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 73); 

                auto tg_xxxz_xxxxxxz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 74); 

                auto tg_xxxz_xxxxxyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 75); 

                auto tg_xxxz_xxxxxyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 76); 

                auto tg_xxxz_xxxxxzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 77); 

                auto tg_xxxz_xxxxyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 78); 

                auto tg_xxxz_xxxxyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 79); 

                auto tg_xxxz_xxxxyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 80); 

                auto tg_xxxz_xxxxzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 81); 

                auto tg_xxxz_xxxyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 82); 

                auto tg_xxxz_xxxyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 83); 

                auto tg_xxxz_xxxyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 84); 

                auto tg_xxxz_xxxyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 85); 

                auto tg_xxxz_xxxzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 86); 

                auto tg_xxxz_xxyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 87); 

                auto tg_xxxz_xxyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 88); 

                auto tg_xxxz_xxyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 89); 

                auto tg_xxxz_xxyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 90); 

                auto tg_xxxz_xxyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 91); 

                auto tg_xxxz_xxzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 92); 

                auto tg_xxxz_xyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 93); 

                auto tg_xxxz_xyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 94); 

                auto tg_xxxx_xxxxxxx_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx); 

                auto tg_xxxx_xxxxxxy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 1); 

                auto tg_xxxx_xxxxxxz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 2); 

                auto tg_xxxx_xxxxxyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 3); 

                auto tg_xxxx_xxxxxyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 4); 

                auto tg_xxxx_xxxxxzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 5); 

                auto tg_xxxx_xxxxyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 6); 

                auto tg_xxxx_xxxxyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 7); 

                auto tg_xxxx_xxxxyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 8); 

                auto tg_xxxx_xxxxzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 9); 

                auto tg_xxxx_xxxyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 10); 

                auto tg_xxxx_xxxyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 11); 

                auto tg_xxxx_xxxyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 12); 

                auto tg_xxxx_xxxyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 13); 

                auto tg_xxxx_xxxzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 14); 

                auto tg_xxxx_xxyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 15); 

                auto tg_xxxx_xxyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 16); 

                auto tg_xxxx_xxyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 17); 

                auto tg_xxxx_xxyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 18); 

                auto tg_xxxx_xxyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 19); 

                auto tg_xxxx_xxzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 20); 

                auto tg_xxxx_xyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 21); 

                auto tg_xxxx_xyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 22); 

                auto tg_xxxx_xyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 23); 

                auto tg_xxxx_xyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 24); 

                auto tg_xxxx_xyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 25); 

                auto tg_xxxx_xyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 26); 

                auto tg_xxxx_xzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 27); 

                auto tg_xxxx_yyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 28); 

                auto tg_xxxx_yyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 29); 

                auto tg_xxxx_yyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 30); 

                auto tg_xxxx_yyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 31); 

                auto tg_xxxx_yyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 32); 

                auto tg_xxxx_yyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 33); 

                auto tg_xxxx_yzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 34); 

                auto tg_xxxx_zzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 35); 

                auto tg_xxxy_xxxxxxx_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 36); 

                auto tg_xxxy_xxxxxxy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 37); 

                auto tg_xxxy_xxxxxxz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 38); 

                auto tg_xxxy_xxxxxyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 39); 

                auto tg_xxxy_xxxxxyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 40); 

                auto tg_xxxy_xxxxxzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 41); 

                auto tg_xxxy_xxxxyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 42); 

                auto tg_xxxy_xxxxyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 43); 

                auto tg_xxxy_xxxxyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 44); 

                auto tg_xxxy_xxxxzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 45); 

                auto tg_xxxy_xxxyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 46); 

                auto tg_xxxy_xxxyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 47); 

                auto tg_xxxy_xxxyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 48); 

                auto tg_xxxy_xxxyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 49); 

                auto tg_xxxy_xxxzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 50); 

                auto tg_xxxy_xxyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 51); 

                auto tg_xxxy_xxyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 52); 

                auto tg_xxxy_xxyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 53); 

                auto tg_xxxy_xxyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 54); 

                auto tg_xxxy_xxyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 55); 

                auto tg_xxxy_xxzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 56); 

                auto tg_xxxy_xyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 57); 

                auto tg_xxxy_xyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 58); 

                auto tg_xxxy_xyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 59); 

                auto tg_xxxy_xyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 60); 

                auto tg_xxxy_xyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 61); 

                auto tg_xxxy_xyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 62); 

                auto tg_xxxy_xzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 63); 

                auto tg_xxxy_yyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 64); 

                auto tg_xxxy_yyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 65); 

                auto tg_xxxy_yyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 66); 

                auto tg_xxxy_yyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 67); 

                auto tg_xxxy_yyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 68); 

                auto tg_xxxy_yyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 69); 

                auto tg_xxxy_yzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 70); 

                auto tg_xxxy_zzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 71); 

                auto tg_xxxz_xxxxxxx_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 72); 

                auto tg_xxxz_xxxxxxy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 73); 

                auto tg_xxxz_xxxxxxz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 74); 

                auto tg_xxxz_xxxxxyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 75); 

                auto tg_xxxz_xxxxxyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 76); 

                auto tg_xxxz_xxxxxzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 77); 

                auto tg_xxxz_xxxxyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 78); 

                auto tg_xxxz_xxxxyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 79); 

                auto tg_xxxz_xxxxyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 80); 

                auto tg_xxxz_xxxxzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 81); 

                auto tg_xxxz_xxxyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 82); 

                auto tg_xxxz_xxxyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 83); 

                auto tg_xxxz_xxxyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 84); 

                auto tg_xxxz_xxxyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 85); 

                auto tg_xxxz_xxxzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 86); 

                auto tg_xxxz_xxyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 87); 

                auto tg_xxxz_xxyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 88); 

                auto tg_xxxz_xxyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 89); 

                auto tg_xxxz_xxyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 90); 

                auto tg_xxxz_xxyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 91); 

                auto tg_xxxz_xxzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 92); 

                auto tg_xxxz_xyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 93); 

                auto tg_xxxz_xyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 94); 

                auto tg_xxx_xxxxxxx_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx); 

                auto tg_xxx_xxxxxxy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 1); 

                auto tg_xxx_xxxxxxz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 2); 

                auto tg_xxx_xxxxxyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 3); 

                auto tg_xxx_xxxxxyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 4); 

                auto tg_xxx_xxxxxzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 5); 

                auto tg_xxx_xxxxyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 6); 

                auto tg_xxx_xxxxyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 7); 

                auto tg_xxx_xxxxyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 8); 

                auto tg_xxx_xxxxzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 9); 

                auto tg_xxx_xxxyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 10); 

                auto tg_xxx_xxxyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 11); 

                auto tg_xxx_xxxyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 12); 

                auto tg_xxx_xxxyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 13); 

                auto tg_xxx_xxxzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 14); 

                auto tg_xxx_xxyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 15); 

                auto tg_xxx_xxyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 16); 

                auto tg_xxx_xxyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 17); 

                auto tg_xxx_xxyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 18); 

                auto tg_xxx_xxyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 19); 

                auto tg_xxx_xxzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 20); 

                auto tg_xxx_xyyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 21); 

                auto tg_xxx_xyyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 22); 

                auto tg_xxx_xyyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 23); 

                auto tg_xxx_xyyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 24); 

                auto tg_xxx_xyyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 25); 

                auto tg_xxx_xyzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 26); 

                auto tg_xxx_xzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 27); 

                auto tg_xxx_yyyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 28); 

                auto tg_xxx_yyyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 29); 

                auto tg_xxx_yyyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 30); 

                auto tg_xxx_yyyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 31); 

                auto tg_xxx_yyyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 32); 

                auto tg_xxx_yyzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 33); 

                auto tg_xxx_yzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 34); 

                auto tg_xxx_zzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 35); 

                auto tg_xxy_xxxxxxx_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 36); 

                auto tg_xxy_xxxxxxy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 37); 

                auto tg_xxy_xxxxxxz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 38); 

                auto tg_xxy_xxxxxyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 39); 

                auto tg_xxy_xxxxxyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 40); 

                auto tg_xxy_xxxxxzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 41); 

                auto tg_xxy_xxxxyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 42); 

                auto tg_xxy_xxxxyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 43); 

                auto tg_xxy_xxxxyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 44); 

                auto tg_xxy_xxxxzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 45); 

                auto tg_xxy_xxxyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 46); 

                auto tg_xxy_xxxyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 47); 

                auto tg_xxy_xxxyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 48); 

                auto tg_xxy_xxxyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 49); 

                auto tg_xxy_xxxzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 50); 

                auto tg_xxy_xxyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 51); 

                auto tg_xxy_xxyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 52); 

                auto tg_xxy_xxyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 53); 

                auto tg_xxy_xxyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 54); 

                auto tg_xxy_xxyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 55); 

                auto tg_xxy_xxzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 56); 

                auto tg_xxy_xyyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 57); 

                auto tg_xxy_xyyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 58); 

                auto tg_xxy_xyyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 59); 

                auto tg_xxy_xyyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 60); 

                auto tg_xxy_xyyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 61); 

                auto tg_xxy_xyzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 62); 

                auto tg_xxy_xzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 63); 

                auto tg_xxy_yyyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 64); 

                auto tg_xxy_yyyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 65); 

                auto tg_xxy_yyyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 66); 

                auto tg_xxy_yyyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 67); 

                auto tg_xxy_yyyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 68); 

                auto tg_xxy_yyzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 69); 

                auto tg_xxy_yzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 70); 

                auto tg_xxy_zzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 71); 

                auto tg_xxz_xxxxxxx_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 72); 

                auto tg_xxz_xxxxxxy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 73); 

                auto tg_xxz_xxxxxxz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 74); 

                auto tg_xxz_xxxxxyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 75); 

                auto tg_xxz_xxxxxyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 76); 

                auto tg_xxz_xxxxxzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 77); 

                auto tg_xxz_xxxxyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 78); 

                auto tg_xxz_xxxxyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 79); 

                auto tg_xxz_xxxxyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 80); 

                auto tg_xxz_xxxxzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 81); 

                auto tg_xxz_xxxyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 82); 

                auto tg_xxz_xxxyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 83); 

                auto tg_xxz_xxxyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 84); 

                auto tg_xxz_xxxyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 85); 

                auto tg_xxz_xxxzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 86); 

                auto tg_xxz_xxyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 87); 

                auto tg_xxz_xxyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 88); 

                auto tg_xxz_xxyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 89); 

                auto tg_xxz_xxyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 90); 

                auto tg_xxz_xxyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 91); 

                auto tg_xxz_xxzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 92); 

                auto tg_xxz_xyyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 93); 

                auto tg_xxz_xyyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 94); 

                auto tg_xxx_xxxxxxx_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx); 

                auto tg_xxx_xxxxxxy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 1); 

                auto tg_xxx_xxxxxxz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 2); 

                auto tg_xxx_xxxxxyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 3); 

                auto tg_xxx_xxxxxyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 4); 

                auto tg_xxx_xxxxxzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 5); 

                auto tg_xxx_xxxxyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 6); 

                auto tg_xxx_xxxxyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 7); 

                auto tg_xxx_xxxxyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 8); 

                auto tg_xxx_xxxxzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 9); 

                auto tg_xxx_xxxyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 10); 

                auto tg_xxx_xxxyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 11); 

                auto tg_xxx_xxxyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 12); 

                auto tg_xxx_xxxyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 13); 

                auto tg_xxx_xxxzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 14); 

                auto tg_xxx_xxyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 15); 

                auto tg_xxx_xxyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 16); 

                auto tg_xxx_xxyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 17); 

                auto tg_xxx_xxyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 18); 

                auto tg_xxx_xxyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 19); 

                auto tg_xxx_xxzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 20); 

                auto tg_xxx_xyyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 21); 

                auto tg_xxx_xyyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 22); 

                auto tg_xxx_xyyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 23); 

                auto tg_xxx_xyyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 24); 

                auto tg_xxx_xyyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 25); 

                auto tg_xxx_xyzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 26); 

                auto tg_xxx_xzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 27); 

                auto tg_xxx_yyyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 28); 

                auto tg_xxx_yyyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 29); 

                auto tg_xxx_yyyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 30); 

                auto tg_xxx_yyyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 31); 

                auto tg_xxx_yyyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 32); 

                auto tg_xxx_yyzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 33); 

                auto tg_xxx_yzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 34); 

                auto tg_xxx_zzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 35); 

                auto tg_xxy_xxxxxxx_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 36); 

                auto tg_xxy_xxxxxxy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 37); 

                auto tg_xxy_xxxxxxz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 38); 

                auto tg_xxy_xxxxxyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 39); 

                auto tg_xxy_xxxxxyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 40); 

                auto tg_xxy_xxxxxzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 41); 

                auto tg_xxy_xxxxyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 42); 

                auto tg_xxy_xxxxyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 43); 

                auto tg_xxy_xxxxyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 44); 

                auto tg_xxy_xxxxzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 45); 

                auto tg_xxy_xxxyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 46); 

                auto tg_xxy_xxxyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 47); 

                auto tg_xxy_xxxyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 48); 

                auto tg_xxy_xxxyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 49); 

                auto tg_xxy_xxxzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 50); 

                auto tg_xxy_xxyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 51); 

                auto tg_xxy_xxyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 52); 

                auto tg_xxy_xxyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 53); 

                auto tg_xxy_xxyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 54); 

                auto tg_xxy_xxyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 55); 

                auto tg_xxy_xxzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 56); 

                auto tg_xxy_xyyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 57); 

                auto tg_xxy_xyyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 58); 

                auto tg_xxy_xyyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 59); 

                auto tg_xxy_xyyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 60); 

                auto tg_xxy_xyyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 61); 

                auto tg_xxy_xyzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 62); 

                auto tg_xxy_xzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 63); 

                auto tg_xxy_yyyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 64); 

                auto tg_xxy_yyyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 65); 

                auto tg_xxy_yyyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 66); 

                auto tg_xxy_yyyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 67); 

                auto tg_xxy_yyyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 68); 

                auto tg_xxy_yyzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 69); 

                auto tg_xxy_yzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 70); 

                auto tg_xxy_zzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 71); 

                auto tg_xxz_xxxxxxx_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 72); 

                auto tg_xxz_xxxxxxy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 73); 

                auto tg_xxz_xxxxxxz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 74); 

                auto tg_xxz_xxxxxyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 75); 

                auto tg_xxz_xxxxxyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 76); 

                auto tg_xxz_xxxxxzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 77); 

                auto tg_xxz_xxxxyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 78); 

                auto tg_xxz_xxxxyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 79); 

                auto tg_xxz_xxxxyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 80); 

                auto tg_xxz_xxxxzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 81); 

                auto tg_xxz_xxxyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 82); 

                auto tg_xxz_xxxyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 83); 

                auto tg_xxz_xxxyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 84); 

                auto tg_xxz_xxxyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 85); 

                auto tg_xxz_xxxzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 86); 

                auto tg_xxz_xxyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 87); 

                auto tg_xxz_xxyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 88); 

                auto tg_xxz_xxyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 89); 

                auto tg_xxz_xxyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 90); 

                auto tg_xxz_xxyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 91); 

                auto tg_xxz_xxzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 92); 

                auto tg_xxz_xyyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 93); 

                auto tg_xxz_xyyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 94); 

                auto tg_xxxx_xxxxxx_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx); 

                auto tg_xxxx_xxxxxy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 1); 

                auto tg_xxxx_xxxxxz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 2); 

                auto tg_xxxx_xxxxyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 3); 

                auto tg_xxxx_xxxxyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 4); 

                auto tg_xxxx_xxxxzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 5); 

                auto tg_xxxx_xxxyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 6); 

                auto tg_xxxx_xxxyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 7); 

                auto tg_xxxx_xxxyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 8); 

                auto tg_xxxx_xxxzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 9); 

                auto tg_xxxx_xxyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 10); 

                auto tg_xxxx_xxyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 11); 

                auto tg_xxxx_xxyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 12); 

                auto tg_xxxx_xxyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 13); 

                auto tg_xxxx_xxzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 14); 

                auto tg_xxxx_xyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 15); 

                auto tg_xxxx_xyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 16); 

                auto tg_xxxx_xyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 17); 

                auto tg_xxxx_xyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 18); 

                auto tg_xxxx_xyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 19); 

                auto tg_xxxx_xzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 20); 

                auto tg_xxxx_yyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 21); 

                auto tg_xxxx_yyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 22); 

                auto tg_xxxx_yyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 23); 

                auto tg_xxxx_yyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 24); 

                auto tg_xxxx_yyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 25); 

                auto tg_xxxx_yzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 26); 

                auto tg_xxxx_zzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 27); 

                auto tg_xxxy_xxxxxx_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 28); 

                auto tg_xxxy_xxxxxy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 29); 

                auto tg_xxxy_xxxxxz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 30); 

                auto tg_xxxy_xxxxyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 31); 

                auto tg_xxxy_xxxxyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 32); 

                auto tg_xxxy_xxxxzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 33); 

                auto tg_xxxy_xxxyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 34); 

                auto tg_xxxy_xxxyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 35); 

                auto tg_xxxy_xxxyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 36); 

                auto tg_xxxy_xxxzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 37); 

                auto tg_xxxy_xxyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 38); 

                auto tg_xxxy_xxyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 39); 

                auto tg_xxxy_xxyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 40); 

                auto tg_xxxy_xxyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 41); 

                auto tg_xxxy_xxzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 42); 

                auto tg_xxxy_xyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 43); 

                auto tg_xxxy_xyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 44); 

                auto tg_xxxy_xyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 45); 

                auto tg_xxxy_xyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 46); 

                auto tg_xxxy_xyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 47); 

                auto tg_xxxy_xzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 48); 

                auto tg_xxxy_yyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 49); 

                auto tg_xxxy_yyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 50); 

                auto tg_xxxy_yyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 51); 

                auto tg_xxxy_yyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 52); 

                auto tg_xxxy_yyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 53); 

                auto tg_xxxy_yzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 54); 

                auto tg_xxxy_zzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 55); 

                auto tg_xxxz_xxxxxx_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 56); 

                auto tg_xxxz_xxxxxy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 57); 

                auto tg_xxxz_xxxxxz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 58); 

                auto tg_xxxz_xxxxyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 59); 

                auto tg_xxxz_xxxxyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 60); 

                auto tg_xxxz_xxxxzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 61); 

                auto tg_xxxz_xxxyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 62); 

                auto tg_xxxz_xxxyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 63); 

                auto tg_xxxz_xxxyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 64); 

                auto tg_xxxz_xxxzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 65); 

                auto tg_xxxz_xxyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 66); 

                auto tg_xxxz_xxyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 67); 

                auto tg_xxxz_xxyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 68); 

                auto tg_xxxz_xxyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 69); 

                auto tg_xxxz_xxzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 70); 

                auto tg_xxxz_xyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 71); 

                auto tg_xxxz_xyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 72); 

                auto tg_xxxz_xyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 73); 

                auto tg_xxxz_xyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 74); 

                auto tg_xxxz_xyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 75); 

                auto tg_xxxz_xzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 76); 

                auto tg_xxxz_yyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 77); 

                auto tg_xxxz_yyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 78); 

                // set up pointers to integrals

                auto tg_xxxxx_xxxxxxx_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx); 

                auto tg_xxxxx_xxxxxxy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 1); 

                auto tg_xxxxx_xxxxxxz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 2); 

                auto tg_xxxxx_xxxxxyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 3); 

                auto tg_xxxxx_xxxxxyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 4); 

                auto tg_xxxxx_xxxxxzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 5); 

                auto tg_xxxxx_xxxxyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 6); 

                auto tg_xxxxx_xxxxyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 7); 

                auto tg_xxxxx_xxxxyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 8); 

                auto tg_xxxxx_xxxxzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 9); 

                auto tg_xxxxx_xxxyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 10); 

                auto tg_xxxxx_xxxyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 11); 

                auto tg_xxxxx_xxxyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 12); 

                auto tg_xxxxx_xxxyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 13); 

                auto tg_xxxxx_xxxzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 14); 

                auto tg_xxxxx_xxyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 15); 

                auto tg_xxxxx_xxyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 16); 

                auto tg_xxxxx_xxyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 17); 

                auto tg_xxxxx_xxyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 18); 

                auto tg_xxxxx_xxyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 19); 

                auto tg_xxxxx_xxzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 20); 

                auto tg_xxxxx_xyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 21); 

                auto tg_xxxxx_xyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 22); 

                auto tg_xxxxx_xyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 23); 

                auto tg_xxxxx_xyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 24); 

                auto tg_xxxxx_xyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 25); 

                auto tg_xxxxx_xyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 26); 

                auto tg_xxxxx_xzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 27); 

                auto tg_xxxxx_yyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 28); 

                auto tg_xxxxx_yyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 29); 

                auto tg_xxxxx_yyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 30); 

                auto tg_xxxxx_yyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 31); 

                auto tg_xxxxx_yyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 32); 

                auto tg_xxxxx_yyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 33); 

                auto tg_xxxxx_yzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 34); 

                auto tg_xxxxx_zzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 35); 

                auto tg_xxxxy_xxxxxxx_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 36); 

                auto tg_xxxxy_xxxxxxy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 37); 

                auto tg_xxxxy_xxxxxxz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 38); 

                auto tg_xxxxy_xxxxxyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 39); 

                auto tg_xxxxy_xxxxxyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 40); 

                auto tg_xxxxy_xxxxxzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 41); 

                auto tg_xxxxy_xxxxyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 42); 

                auto tg_xxxxy_xxxxyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 43); 

                auto tg_xxxxy_xxxxyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 44); 

                auto tg_xxxxy_xxxxzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 45); 

                auto tg_xxxxy_xxxyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 46); 

                auto tg_xxxxy_xxxyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 47); 

                auto tg_xxxxy_xxxyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 48); 

                auto tg_xxxxy_xxxyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 49); 

                auto tg_xxxxy_xxxzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 50); 

                auto tg_xxxxy_xxyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 51); 

                auto tg_xxxxy_xxyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 52); 

                auto tg_xxxxy_xxyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 53); 

                auto tg_xxxxy_xxyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 54); 

                auto tg_xxxxy_xxyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 55); 

                auto tg_xxxxy_xxzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 56); 

                auto tg_xxxxy_xyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 57); 

                auto tg_xxxxy_xyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 58); 

                auto tg_xxxxy_xyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 59); 

                auto tg_xxxxy_xyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 60); 

                auto tg_xxxxy_xyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 61); 

                auto tg_xxxxy_xyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 62); 

                auto tg_xxxxy_xzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 63); 

                auto tg_xxxxy_yyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 64); 

                auto tg_xxxxy_yyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 65); 

                auto tg_xxxxy_yyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 66); 

                auto tg_xxxxy_yyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 67); 

                auto tg_xxxxy_yyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 68); 

                auto tg_xxxxy_yyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 69); 

                auto tg_xxxxy_yzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 70); 

                auto tg_xxxxy_zzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 71); 

                auto tg_xxxxz_xxxxxxx_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 72); 

                auto tg_xxxxz_xxxxxxy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 73); 

                auto tg_xxxxz_xxxxxxz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 74); 

                auto tg_xxxxz_xxxxxyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 75); 

                auto tg_xxxxz_xxxxxyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 76); 

                auto tg_xxxxz_xxxxxzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 77); 

                auto tg_xxxxz_xxxxyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 78); 

                auto tg_xxxxz_xxxxyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 79); 

                auto tg_xxxxz_xxxxyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 80); 

                auto tg_xxxxz_xxxxzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 81); 

                auto tg_xxxxz_xxxyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 82); 

                auto tg_xxxxz_xxxyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 83); 

                auto tg_xxxxz_xxxyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 84); 

                auto tg_xxxxz_xxxyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 85); 

                auto tg_xxxxz_xxxzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 86); 

                auto tg_xxxxz_xxyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 87); 

                auto tg_xxxxz_xxyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 88); 

                auto tg_xxxxz_xxyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 89); 

                auto tg_xxxxz_xxyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 90); 

                auto tg_xxxxz_xxyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 91); 

                auto tg_xxxxz_xxzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 92); 

                auto tg_xxxxz_xyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 93); 

                auto tg_xxxxz_xyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 94); 

                // Batch of Integrals (0,95)

                #pragma omp simd aligned(fxn, fza, tg_xxx_xxxxxxx_0, tg_xxx_xxxxxxx_1, tg_xxx_xxxxxxy_0, \
                                         tg_xxx_xxxxxxy_1, tg_xxx_xxxxxxz_0, tg_xxx_xxxxxxz_1, tg_xxx_xxxxxyy_0, \
                                         tg_xxx_xxxxxyy_1, tg_xxx_xxxxxyz_0, tg_xxx_xxxxxyz_1, tg_xxx_xxxxxzz_0, \
                                         tg_xxx_xxxxxzz_1, tg_xxx_xxxxyyy_0, tg_xxx_xxxxyyy_1, tg_xxx_xxxxyyz_0, \
                                         tg_xxx_xxxxyyz_1, tg_xxx_xxxxyzz_0, tg_xxx_xxxxyzz_1, tg_xxx_xxxxzzz_0, \
                                         tg_xxx_xxxxzzz_1, tg_xxx_xxxyyyy_0, tg_xxx_xxxyyyy_1, tg_xxx_xxxyyyz_0, \
                                         tg_xxx_xxxyyyz_1, tg_xxx_xxxyyzz_0, tg_xxx_xxxyyzz_1, tg_xxx_xxxyzzz_0, \
                                         tg_xxx_xxxyzzz_1, tg_xxx_xxxzzzz_0, tg_xxx_xxxzzzz_1, tg_xxx_xxyyyyy_0, \
                                         tg_xxx_xxyyyyy_1, tg_xxx_xxyyyyz_0, tg_xxx_xxyyyyz_1, tg_xxx_xxyyyzz_0, \
                                         tg_xxx_xxyyyzz_1, tg_xxx_xxyyzzz_0, tg_xxx_xxyyzzz_1, tg_xxx_xxyzzzz_0, \
                                         tg_xxx_xxyzzzz_1, tg_xxx_xxzzzzz_0, tg_xxx_xxzzzzz_1, tg_xxx_xyyyyyy_0, \
                                         tg_xxx_xyyyyyy_1, tg_xxx_xyyyyyz_0, tg_xxx_xyyyyyz_1, tg_xxx_xyyyyzz_0, \
                                         tg_xxx_xyyyyzz_1, tg_xxx_xyyyzzz_0, tg_xxx_xyyyzzz_1, tg_xxx_xyyzzzz_0, \
                                         tg_xxx_xyyzzzz_1, tg_xxx_xyzzzzz_0, tg_xxx_xyzzzzz_1, tg_xxx_xzzzzzz_0, \
                                         tg_xxx_xzzzzzz_1, tg_xxx_yyyyyyy_0, tg_xxx_yyyyyyy_1, tg_xxx_yyyyyyz_0, \
                                         tg_xxx_yyyyyyz_1, tg_xxx_yyyyyzz_0, tg_xxx_yyyyyzz_1, tg_xxx_yyyyzzz_0, \
                                         tg_xxx_yyyyzzz_1, tg_xxx_yyyzzzz_0, tg_xxx_yyyzzzz_1, tg_xxx_yyzzzzz_0, \
                                         tg_xxx_yyzzzzz_1, tg_xxx_yzzzzzz_0, tg_xxx_yzzzzzz_1, tg_xxx_zzzzzzz_0, \
                                         tg_xxx_zzzzzzz_1, tg_xxxx_xxxxxx_1, tg_xxxx_xxxxxxx_0, tg_xxxx_xxxxxxx_1, \
                                         tg_xxxx_xxxxxxy_0, tg_xxxx_xxxxxxy_1, tg_xxxx_xxxxxxz_0, tg_xxxx_xxxxxxz_1, \
                                         tg_xxxx_xxxxxy_1, tg_xxxx_xxxxxyy_0, tg_xxxx_xxxxxyy_1, tg_xxxx_xxxxxyz_0, \
                                         tg_xxxx_xxxxxyz_1, tg_xxxx_xxxxxz_1, tg_xxxx_xxxxxzz_0, tg_xxxx_xxxxxzz_1, \
                                         tg_xxxx_xxxxyy_1, tg_xxxx_xxxxyyy_0, tg_xxxx_xxxxyyy_1, tg_xxxx_xxxxyyz_0, \
                                         tg_xxxx_xxxxyyz_1, tg_xxxx_xxxxyz_1, tg_xxxx_xxxxyzz_0, tg_xxxx_xxxxyzz_1, \
                                         tg_xxxx_xxxxzz_1, tg_xxxx_xxxxzzz_0, tg_xxxx_xxxxzzz_1, tg_xxxx_xxxyyy_1, \
                                         tg_xxxx_xxxyyyy_0, tg_xxxx_xxxyyyy_1, tg_xxxx_xxxyyyz_0, tg_xxxx_xxxyyyz_1, \
                                         tg_xxxx_xxxyyz_1, tg_xxxx_xxxyyzz_0, tg_xxxx_xxxyyzz_1, tg_xxxx_xxxyzz_1, \
                                         tg_xxxx_xxxyzzz_0, tg_xxxx_xxxyzzz_1, tg_xxxx_xxxzzz_1, tg_xxxx_xxxzzzz_0, \
                                         tg_xxxx_xxxzzzz_1, tg_xxxx_xxyyyy_1, tg_xxxx_xxyyyyy_0, tg_xxxx_xxyyyyy_1, \
                                         tg_xxxx_xxyyyyz_0, tg_xxxx_xxyyyyz_1, tg_xxxx_xxyyyz_1, tg_xxxx_xxyyyzz_0, \
                                         tg_xxxx_xxyyyzz_1, tg_xxxx_xxyyzz_1, tg_xxxx_xxyyzzz_0, tg_xxxx_xxyyzzz_1, \
                                         tg_xxxx_xxyzzz_1, tg_xxxx_xxyzzzz_0, tg_xxxx_xxyzzzz_1, tg_xxxx_xxzzzz_1, \
                                         tg_xxxx_xxzzzzz_0, tg_xxxx_xxzzzzz_1, tg_xxxx_xyyyyy_1, tg_xxxx_xyyyyyy_0, \
                                         tg_xxxx_xyyyyyy_1, tg_xxxx_xyyyyyz_0, tg_xxxx_xyyyyyz_1, tg_xxxx_xyyyyz_1, \
                                         tg_xxxx_xyyyyzz_0, tg_xxxx_xyyyyzz_1, tg_xxxx_xyyyzz_1, tg_xxxx_xyyyzzz_0, \
                                         tg_xxxx_xyyyzzz_1, tg_xxxx_xyyzzz_1, tg_xxxx_xyyzzzz_0, tg_xxxx_xyyzzzz_1, \
                                         tg_xxxx_xyzzzz_1, tg_xxxx_xyzzzzz_0, tg_xxxx_xyzzzzz_1, tg_xxxx_xzzzzz_1, \
                                         tg_xxxx_xzzzzzz_0, tg_xxxx_xzzzzzz_1, tg_xxxx_yyyyyy_1, tg_xxxx_yyyyyyy_0, \
                                         tg_xxxx_yyyyyyy_1, tg_xxxx_yyyyyyz_0, tg_xxxx_yyyyyyz_1, tg_xxxx_yyyyyz_1, \
                                         tg_xxxx_yyyyyzz_0, tg_xxxx_yyyyyzz_1, tg_xxxx_yyyyzz_1, tg_xxxx_yyyyzzz_0, \
                                         tg_xxxx_yyyyzzz_1, tg_xxxx_yyyzzz_1, tg_xxxx_yyyzzzz_0, tg_xxxx_yyyzzzz_1, \
                                         tg_xxxx_yyzzzz_1, tg_xxxx_yyzzzzz_0, tg_xxxx_yyzzzzz_1, tg_xxxx_yzzzzz_1, \
                                         tg_xxxx_yzzzzzz_0, tg_xxxx_yzzzzzz_1, tg_xxxx_zzzzzz_1, tg_xxxx_zzzzzzz_0, \
                                         tg_xxxx_zzzzzzz_1, tg_xxxxx_xxxxxxx_0, tg_xxxxx_xxxxxxy_0, tg_xxxxx_xxxxxxz_0, \
                                         tg_xxxxx_xxxxxyy_0, tg_xxxxx_xxxxxyz_0, tg_xxxxx_xxxxxzz_0, tg_xxxxx_xxxxyyy_0, \
                                         tg_xxxxx_xxxxyyz_0, tg_xxxxx_xxxxyzz_0, tg_xxxxx_xxxxzzz_0, tg_xxxxx_xxxyyyy_0, \
                                         tg_xxxxx_xxxyyyz_0, tg_xxxxx_xxxyyzz_0, tg_xxxxx_xxxyzzz_0, tg_xxxxx_xxxzzzz_0, \
                                         tg_xxxxx_xxyyyyy_0, tg_xxxxx_xxyyyyz_0, tg_xxxxx_xxyyyzz_0, tg_xxxxx_xxyyzzz_0, \
                                         tg_xxxxx_xxyzzzz_0, tg_xxxxx_xxzzzzz_0, tg_xxxxx_xyyyyyy_0, tg_xxxxx_xyyyyyz_0, \
                                         tg_xxxxx_xyyyyzz_0, tg_xxxxx_xyyyzzz_0, tg_xxxxx_xyyzzzz_0, tg_xxxxx_xyzzzzz_0, \
                                         tg_xxxxx_xzzzzzz_0, tg_xxxxx_yyyyyyy_0, tg_xxxxx_yyyyyyz_0, tg_xxxxx_yyyyyzz_0, \
                                         tg_xxxxx_yyyyzzz_0, tg_xxxxx_yyyzzzz_0, tg_xxxxx_yyzzzzz_0, tg_xxxxx_yzzzzzz_0, \
                                         tg_xxxxx_zzzzzzz_0, tg_xxxxy_xxxxxxx_0, tg_xxxxy_xxxxxxy_0, tg_xxxxy_xxxxxxz_0, \
                                         tg_xxxxy_xxxxxyy_0, tg_xxxxy_xxxxxyz_0, tg_xxxxy_xxxxxzz_0, tg_xxxxy_xxxxyyy_0, \
                                         tg_xxxxy_xxxxyyz_0, tg_xxxxy_xxxxyzz_0, tg_xxxxy_xxxxzzz_0, tg_xxxxy_xxxyyyy_0, \
                                         tg_xxxxy_xxxyyyz_0, tg_xxxxy_xxxyyzz_0, tg_xxxxy_xxxyzzz_0, tg_xxxxy_xxxzzzz_0, \
                                         tg_xxxxy_xxyyyyy_0, tg_xxxxy_xxyyyyz_0, tg_xxxxy_xxyyyzz_0, tg_xxxxy_xxyyzzz_0, \
                                         tg_xxxxy_xxyzzzz_0, tg_xxxxy_xxzzzzz_0, tg_xxxxy_xyyyyyy_0, tg_xxxxy_xyyyyyz_0, \
                                         tg_xxxxy_xyyyyzz_0, tg_xxxxy_xyyyzzz_0, tg_xxxxy_xyyzzzz_0, tg_xxxxy_xyzzzzz_0, \
                                         tg_xxxxy_xzzzzzz_0, tg_xxxxy_yyyyyyy_0, tg_xxxxy_yyyyyyz_0, tg_xxxxy_yyyyyzz_0, \
                                         tg_xxxxy_yyyyzzz_0, tg_xxxxy_yyyzzzz_0, tg_xxxxy_yyzzzzz_0, tg_xxxxy_yzzzzzz_0, \
                                         tg_xxxxy_zzzzzzz_0, tg_xxxxz_xxxxxxx_0, tg_xxxxz_xxxxxxy_0, tg_xxxxz_xxxxxxz_0, \
                                         tg_xxxxz_xxxxxyy_0, tg_xxxxz_xxxxxyz_0, tg_xxxxz_xxxxxzz_0, tg_xxxxz_xxxxyyy_0, \
                                         tg_xxxxz_xxxxyyz_0, tg_xxxxz_xxxxyzz_0, tg_xxxxz_xxxxzzz_0, tg_xxxxz_xxxyyyy_0, \
                                         tg_xxxxz_xxxyyyz_0, tg_xxxxz_xxxyyzz_0, tg_xxxxz_xxxyzzz_0, tg_xxxxz_xxxzzzz_0, \
                                         tg_xxxxz_xxyyyyy_0, tg_xxxxz_xxyyyyz_0, tg_xxxxz_xxyyyzz_0, tg_xxxxz_xxyyzzz_0, \
                                         tg_xxxxz_xxyzzzz_0, tg_xxxxz_xxzzzzz_0, tg_xxxxz_xyyyyyy_0, tg_xxxxz_xyyyyyz_0, \
                                         tg_xxxy_xxxxxx_1, tg_xxxy_xxxxxxx_0, tg_xxxy_xxxxxxx_1, tg_xxxy_xxxxxxy_0, \
                                         tg_xxxy_xxxxxxy_1, tg_xxxy_xxxxxxz_0, tg_xxxy_xxxxxxz_1, tg_xxxy_xxxxxy_1, \
                                         tg_xxxy_xxxxxyy_0, tg_xxxy_xxxxxyy_1, tg_xxxy_xxxxxyz_0, tg_xxxy_xxxxxyz_1, \
                                         tg_xxxy_xxxxxz_1, tg_xxxy_xxxxxzz_0, tg_xxxy_xxxxxzz_1, tg_xxxy_xxxxyy_1, \
                                         tg_xxxy_xxxxyyy_0, tg_xxxy_xxxxyyy_1, tg_xxxy_xxxxyyz_0, tg_xxxy_xxxxyyz_1, \
                                         tg_xxxy_xxxxyz_1, tg_xxxy_xxxxyzz_0, tg_xxxy_xxxxyzz_1, tg_xxxy_xxxxzz_1, \
                                         tg_xxxy_xxxxzzz_0, tg_xxxy_xxxxzzz_1, tg_xxxy_xxxyyy_1, tg_xxxy_xxxyyyy_0, \
                                         tg_xxxy_xxxyyyy_1, tg_xxxy_xxxyyyz_0, tg_xxxy_xxxyyyz_1, tg_xxxy_xxxyyz_1, \
                                         tg_xxxy_xxxyyzz_0, tg_xxxy_xxxyyzz_1, tg_xxxy_xxxyzz_1, tg_xxxy_xxxyzzz_0, \
                                         tg_xxxy_xxxyzzz_1, tg_xxxy_xxxzzz_1, tg_xxxy_xxxzzzz_0, tg_xxxy_xxxzzzz_1, \
                                         tg_xxxy_xxyyyy_1, tg_xxxy_xxyyyyy_0, tg_xxxy_xxyyyyy_1, tg_xxxy_xxyyyyz_0, \
                                         tg_xxxy_xxyyyyz_1, tg_xxxy_xxyyyz_1, tg_xxxy_xxyyyzz_0, tg_xxxy_xxyyyzz_1, \
                                         tg_xxxy_xxyyzz_1, tg_xxxy_xxyyzzz_0, tg_xxxy_xxyyzzz_1, tg_xxxy_xxyzzz_1, \
                                         tg_xxxy_xxyzzzz_0, tg_xxxy_xxyzzzz_1, tg_xxxy_xxzzzz_1, tg_xxxy_xxzzzzz_0, \
                                         tg_xxxy_xxzzzzz_1, tg_xxxy_xyyyyy_1, tg_xxxy_xyyyyyy_0, tg_xxxy_xyyyyyy_1, \
                                         tg_xxxy_xyyyyyz_0, tg_xxxy_xyyyyyz_1, tg_xxxy_xyyyyz_1, tg_xxxy_xyyyyzz_0, \
                                         tg_xxxy_xyyyyzz_1, tg_xxxy_xyyyzz_1, tg_xxxy_xyyyzzz_0, tg_xxxy_xyyyzzz_1, \
                                         tg_xxxy_xyyzzz_1, tg_xxxy_xyyzzzz_0, tg_xxxy_xyyzzzz_1, tg_xxxy_xyzzzz_1, \
                                         tg_xxxy_xyzzzzz_0, tg_xxxy_xyzzzzz_1, tg_xxxy_xzzzzz_1, tg_xxxy_xzzzzzz_0, \
                                         tg_xxxy_xzzzzzz_1, tg_xxxy_yyyyyy_1, tg_xxxy_yyyyyyy_0, tg_xxxy_yyyyyyy_1, \
                                         tg_xxxy_yyyyyyz_0, tg_xxxy_yyyyyyz_1, tg_xxxy_yyyyyz_1, tg_xxxy_yyyyyzz_0, \
                                         tg_xxxy_yyyyyzz_1, tg_xxxy_yyyyzz_1, tg_xxxy_yyyyzzz_0, tg_xxxy_yyyyzzz_1, \
                                         tg_xxxy_yyyzzz_1, tg_xxxy_yyyzzzz_0, tg_xxxy_yyyzzzz_1, tg_xxxy_yyzzzz_1, \
                                         tg_xxxy_yyzzzzz_0, tg_xxxy_yyzzzzz_1, tg_xxxy_yzzzzz_1, tg_xxxy_yzzzzzz_0, \
                                         tg_xxxy_yzzzzzz_1, tg_xxxy_zzzzzz_1, tg_xxxy_zzzzzzz_0, tg_xxxy_zzzzzzz_1, \
                                         tg_xxxz_xxxxxx_1, tg_xxxz_xxxxxxx_0, tg_xxxz_xxxxxxx_1, tg_xxxz_xxxxxxy_0, \
                                         tg_xxxz_xxxxxxy_1, tg_xxxz_xxxxxxz_0, tg_xxxz_xxxxxxz_1, tg_xxxz_xxxxxy_1, \
                                         tg_xxxz_xxxxxyy_0, tg_xxxz_xxxxxyy_1, tg_xxxz_xxxxxyz_0, tg_xxxz_xxxxxyz_1, \
                                         tg_xxxz_xxxxxz_1, tg_xxxz_xxxxxzz_0, tg_xxxz_xxxxxzz_1, tg_xxxz_xxxxyy_1, \
                                         tg_xxxz_xxxxyyy_0, tg_xxxz_xxxxyyy_1, tg_xxxz_xxxxyyz_0, tg_xxxz_xxxxyyz_1, \
                                         tg_xxxz_xxxxyz_1, tg_xxxz_xxxxyzz_0, tg_xxxz_xxxxyzz_1, tg_xxxz_xxxxzz_1, \
                                         tg_xxxz_xxxxzzz_0, tg_xxxz_xxxxzzz_1, tg_xxxz_xxxyyy_1, tg_xxxz_xxxyyyy_0, \
                                         tg_xxxz_xxxyyyy_1, tg_xxxz_xxxyyyz_0, tg_xxxz_xxxyyyz_1, tg_xxxz_xxxyyz_1, \
                                         tg_xxxz_xxxyyzz_0, tg_xxxz_xxxyyzz_1, tg_xxxz_xxxyzz_1, tg_xxxz_xxxyzzz_0, \
                                         tg_xxxz_xxxyzzz_1, tg_xxxz_xxxzzz_1, tg_xxxz_xxxzzzz_0, tg_xxxz_xxxzzzz_1, \
                                         tg_xxxz_xxyyyy_1, tg_xxxz_xxyyyyy_0, tg_xxxz_xxyyyyy_1, tg_xxxz_xxyyyyz_0, \
                                         tg_xxxz_xxyyyyz_1, tg_xxxz_xxyyyz_1, tg_xxxz_xxyyyzz_0, tg_xxxz_xxyyyzz_1, \
                                         tg_xxxz_xxyyzz_1, tg_xxxz_xxyyzzz_0, tg_xxxz_xxyyzzz_1, tg_xxxz_xxyzzz_1, \
                                         tg_xxxz_xxyzzzz_0, tg_xxxz_xxyzzzz_1, tg_xxxz_xxzzzz_1, tg_xxxz_xxzzzzz_0, \
                                         tg_xxxz_xxzzzzz_1, tg_xxxz_xyyyyy_1, tg_xxxz_xyyyyyy_0, tg_xxxz_xyyyyyy_1, \
                                         tg_xxxz_xyyyyyz_0, tg_xxxz_xyyyyyz_1, tg_xxxz_xyyyyz_1, tg_xxxz_xyyyzz_1, \
                                         tg_xxxz_xyyzzz_1, tg_xxxz_xyzzzz_1, tg_xxxz_xzzzzz_1, tg_xxxz_yyyyyy_1, \
                                         tg_xxxz_yyyyyz_1, tg_xxy_xxxxxxx_0, tg_xxy_xxxxxxx_1, tg_xxy_xxxxxxy_0, \
                                         tg_xxy_xxxxxxy_1, tg_xxy_xxxxxxz_0, tg_xxy_xxxxxxz_1, tg_xxy_xxxxxyy_0, \
                                         tg_xxy_xxxxxyy_1, tg_xxy_xxxxxyz_0, tg_xxy_xxxxxyz_1, tg_xxy_xxxxxzz_0, \
                                         tg_xxy_xxxxxzz_1, tg_xxy_xxxxyyy_0, tg_xxy_xxxxyyy_1, tg_xxy_xxxxyyz_0, \
                                         tg_xxy_xxxxyyz_1, tg_xxy_xxxxyzz_0, tg_xxy_xxxxyzz_1, tg_xxy_xxxxzzz_0, \
                                         tg_xxy_xxxxzzz_1, tg_xxy_xxxyyyy_0, tg_xxy_xxxyyyy_1, tg_xxy_xxxyyyz_0, \
                                         tg_xxy_xxxyyyz_1, tg_xxy_xxxyyzz_0, tg_xxy_xxxyyzz_1, tg_xxy_xxxyzzz_0, \
                                         tg_xxy_xxxyzzz_1, tg_xxy_xxxzzzz_0, tg_xxy_xxxzzzz_1, tg_xxy_xxyyyyy_0, \
                                         tg_xxy_xxyyyyy_1, tg_xxy_xxyyyyz_0, tg_xxy_xxyyyyz_1, tg_xxy_xxyyyzz_0, \
                                         tg_xxy_xxyyyzz_1, tg_xxy_xxyyzzz_0, tg_xxy_xxyyzzz_1, tg_xxy_xxyzzzz_0, \
                                         tg_xxy_xxyzzzz_1, tg_xxy_xxzzzzz_0, tg_xxy_xxzzzzz_1, tg_xxy_xyyyyyy_0, \
                                         tg_xxy_xyyyyyy_1, tg_xxy_xyyyyyz_0, tg_xxy_xyyyyyz_1, tg_xxy_xyyyyzz_0, \
                                         tg_xxy_xyyyyzz_1, tg_xxy_xyyyzzz_0, tg_xxy_xyyyzzz_1, tg_xxy_xyyzzzz_0, \
                                         tg_xxy_xyyzzzz_1, tg_xxy_xyzzzzz_0, tg_xxy_xyzzzzz_1, tg_xxy_xzzzzzz_0, \
                                         tg_xxy_xzzzzzz_1, tg_xxy_yyyyyyy_0, tg_xxy_yyyyyyy_1, tg_xxy_yyyyyyz_0, \
                                         tg_xxy_yyyyyyz_1, tg_xxy_yyyyyzz_0, tg_xxy_yyyyyzz_1, tg_xxy_yyyyzzz_0, \
                                         tg_xxy_yyyyzzz_1, tg_xxy_yyyzzzz_0, tg_xxy_yyyzzzz_1, tg_xxy_yyzzzzz_0, \
                                         tg_xxy_yyzzzzz_1, tg_xxy_yzzzzzz_0, tg_xxy_yzzzzzz_1, tg_xxy_zzzzzzz_0, \
                                         tg_xxy_zzzzzzz_1, tg_xxz_xxxxxxx_0, tg_xxz_xxxxxxx_1, tg_xxz_xxxxxxy_0, \
                                         tg_xxz_xxxxxxy_1, tg_xxz_xxxxxxz_0, tg_xxz_xxxxxxz_1, tg_xxz_xxxxxyy_0, \
                                         tg_xxz_xxxxxyy_1, tg_xxz_xxxxxyz_0, tg_xxz_xxxxxyz_1, tg_xxz_xxxxxzz_0, \
                                         tg_xxz_xxxxxzz_1, tg_xxz_xxxxyyy_0, tg_xxz_xxxxyyy_1, tg_xxz_xxxxyyz_0, \
                                         tg_xxz_xxxxyyz_1, tg_xxz_xxxxyzz_0, tg_xxz_xxxxyzz_1, tg_xxz_xxxxzzz_0, \
                                         tg_xxz_xxxxzzz_1, tg_xxz_xxxyyyy_0, tg_xxz_xxxyyyy_1, tg_xxz_xxxyyyz_0, \
                                         tg_xxz_xxxyyyz_1, tg_xxz_xxxyyzz_0, tg_xxz_xxxyyzz_1, tg_xxz_xxxyzzz_0, \
                                         tg_xxz_xxxyzzz_1, tg_xxz_xxxzzzz_0, tg_xxz_xxxzzzz_1, tg_xxz_xxyyyyy_0, \
                                         tg_xxz_xxyyyyy_1, tg_xxz_xxyyyyz_0, tg_xxz_xxyyyyz_1, tg_xxz_xxyyyzz_0, \
                                         tg_xxz_xxyyyzz_1, tg_xxz_xxyyzzz_0, tg_xxz_xxyyzzz_1, tg_xxz_xxyzzzz_0, \
                                         tg_xxz_xxyzzzz_1, tg_xxz_xxzzzzz_0, tg_xxz_xxzzzzz_1, tg_xxz_xyyyyyy_0, \
                                         tg_xxz_xyyyyyy_1, tg_xxz_xyyyyyz_0, tg_xxz_xyyyyyz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxxxx_xxxxxxx_0[j] = pb_x * tg_xxxx_xxxxxxx_0[j] + wp_x[j] * tg_xxxx_xxxxxxx_1[j] + 2.0 * fl1_fx * tg_xxx_xxxxxxx_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxxxxxx_1[j] + 3.5 * fl1_fxn * tg_xxxx_xxxxxx_1[j];

                    tg_xxxxx_xxxxxxy_0[j] = pb_x * tg_xxxx_xxxxxxy_0[j] + wp_x[j] * tg_xxxx_xxxxxxy_1[j] + 2.0 * fl1_fx * tg_xxx_xxxxxxy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxxxxxy_1[j] + 3.0 * fl1_fxn * tg_xxxx_xxxxxy_1[j];

                    tg_xxxxx_xxxxxxz_0[j] = pb_x * tg_xxxx_xxxxxxz_0[j] + wp_x[j] * tg_xxxx_xxxxxxz_1[j] + 2.0 * fl1_fx * tg_xxx_xxxxxxz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxxxxxz_1[j] + 3.0 * fl1_fxn * tg_xxxx_xxxxxz_1[j];

                    tg_xxxxx_xxxxxyy_0[j] = pb_x * tg_xxxx_xxxxxyy_0[j] + wp_x[j] * tg_xxxx_xxxxxyy_1[j] + 2.0 * fl1_fx * tg_xxx_xxxxxyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxxxxyy_1[j] + 2.5 * fl1_fxn * tg_xxxx_xxxxyy_1[j];

                    tg_xxxxx_xxxxxyz_0[j] = pb_x * tg_xxxx_xxxxxyz_0[j] + wp_x[j] * tg_xxxx_xxxxxyz_1[j] + 2.0 * fl1_fx * tg_xxx_xxxxxyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxxxxyz_1[j] + 2.5 * fl1_fxn * tg_xxxx_xxxxyz_1[j];

                    tg_xxxxx_xxxxxzz_0[j] = pb_x * tg_xxxx_xxxxxzz_0[j] + wp_x[j] * tg_xxxx_xxxxxzz_1[j] + 2.0 * fl1_fx * tg_xxx_xxxxxzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxxxxzz_1[j] + 2.5 * fl1_fxn * tg_xxxx_xxxxzz_1[j];

                    tg_xxxxx_xxxxyyy_0[j] = pb_x * tg_xxxx_xxxxyyy_0[j] + wp_x[j] * tg_xxxx_xxxxyyy_1[j] + 2.0 * fl1_fx * tg_xxx_xxxxyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxxxyyy_1[j] + 2.0 * fl1_fxn * tg_xxxx_xxxyyy_1[j];

                    tg_xxxxx_xxxxyyz_0[j] = pb_x * tg_xxxx_xxxxyyz_0[j] + wp_x[j] * tg_xxxx_xxxxyyz_1[j] + 2.0 * fl1_fx * tg_xxx_xxxxyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxxxyyz_1[j] + 2.0 * fl1_fxn * tg_xxxx_xxxyyz_1[j];

                    tg_xxxxx_xxxxyzz_0[j] = pb_x * tg_xxxx_xxxxyzz_0[j] + wp_x[j] * tg_xxxx_xxxxyzz_1[j] + 2.0 * fl1_fx * tg_xxx_xxxxyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxxxyzz_1[j] + 2.0 * fl1_fxn * tg_xxxx_xxxyzz_1[j];

                    tg_xxxxx_xxxxzzz_0[j] = pb_x * tg_xxxx_xxxxzzz_0[j] + wp_x[j] * tg_xxxx_xxxxzzz_1[j] + 2.0 * fl1_fx * tg_xxx_xxxxzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxxxzzz_1[j] + 2.0 * fl1_fxn * tg_xxxx_xxxzzz_1[j];

                    tg_xxxxx_xxxyyyy_0[j] = pb_x * tg_xxxx_xxxyyyy_0[j] + wp_x[j] * tg_xxxx_xxxyyyy_1[j] + 2.0 * fl1_fx * tg_xxx_xxxyyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxxyyyy_1[j] + 1.5 * fl1_fxn * tg_xxxx_xxyyyy_1[j];

                    tg_xxxxx_xxxyyyz_0[j] = pb_x * tg_xxxx_xxxyyyz_0[j] + wp_x[j] * tg_xxxx_xxxyyyz_1[j] + 2.0 * fl1_fx * tg_xxx_xxxyyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_xxxx_xxyyyz_1[j];

                    tg_xxxxx_xxxyyzz_0[j] = pb_x * tg_xxxx_xxxyyzz_0[j] + wp_x[j] * tg_xxxx_xxxyyzz_1[j] + 2.0 * fl1_fx * tg_xxx_xxxyyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxxyyzz_1[j] + 1.5 * fl1_fxn * tg_xxxx_xxyyzz_1[j];

                    tg_xxxxx_xxxyzzz_0[j] = pb_x * tg_xxxx_xxxyzzz_0[j] + wp_x[j] * tg_xxxx_xxxyzzz_1[j] + 2.0 * fl1_fx * tg_xxx_xxxyzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxxyzzz_1[j] + 1.5 * fl1_fxn * tg_xxxx_xxyzzz_1[j];

                    tg_xxxxx_xxxzzzz_0[j] = pb_x * tg_xxxx_xxxzzzz_0[j] + wp_x[j] * tg_xxxx_xxxzzzz_1[j] + 2.0 * fl1_fx * tg_xxx_xxxzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxxzzzz_1[j] + 1.5 * fl1_fxn * tg_xxxx_xxzzzz_1[j];

                    tg_xxxxx_xxyyyyy_0[j] = pb_x * tg_xxxx_xxyyyyy_0[j] + wp_x[j] * tg_xxxx_xxyyyyy_1[j] + 2.0 * fl1_fx * tg_xxx_xxyyyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxyyyyy_1[j] + fl1_fxn * tg_xxxx_xyyyyy_1[j];

                    tg_xxxxx_xxyyyyz_0[j] = pb_x * tg_xxxx_xxyyyyz_0[j] + wp_x[j] * tg_xxxx_xxyyyyz_1[j] + 2.0 * fl1_fx * tg_xxx_xxyyyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxyyyyz_1[j] + fl1_fxn * tg_xxxx_xyyyyz_1[j];

                    tg_xxxxx_xxyyyzz_0[j] = pb_x * tg_xxxx_xxyyyzz_0[j] + wp_x[j] * tg_xxxx_xxyyyzz_1[j] + 2.0 * fl1_fx * tg_xxx_xxyyyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxyyyzz_1[j] + fl1_fxn * tg_xxxx_xyyyzz_1[j];

                    tg_xxxxx_xxyyzzz_0[j] = pb_x * tg_xxxx_xxyyzzz_0[j] + wp_x[j] * tg_xxxx_xxyyzzz_1[j] + 2.0 * fl1_fx * tg_xxx_xxyyzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxyyzzz_1[j] + fl1_fxn * tg_xxxx_xyyzzz_1[j];

                    tg_xxxxx_xxyzzzz_0[j] = pb_x * tg_xxxx_xxyzzzz_0[j] + wp_x[j] * tg_xxxx_xxyzzzz_1[j] + 2.0 * fl1_fx * tg_xxx_xxyzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxyzzzz_1[j] + fl1_fxn * tg_xxxx_xyzzzz_1[j];

                    tg_xxxxx_xxzzzzz_0[j] = pb_x * tg_xxxx_xxzzzzz_0[j] + wp_x[j] * tg_xxxx_xxzzzzz_1[j] + 2.0 * fl1_fx * tg_xxx_xxzzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxzzzzz_1[j] + fl1_fxn * tg_xxxx_xzzzzz_1[j];

                    tg_xxxxx_xyyyyyy_0[j] = pb_x * tg_xxxx_xyyyyyy_0[j] + wp_x[j] * tg_xxxx_xyyyyyy_1[j] + 2.0 * fl1_fx * tg_xxx_xyyyyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xyyyyyy_1[j] + 0.5 * fl1_fxn * tg_xxxx_yyyyyy_1[j];

                    tg_xxxxx_xyyyyyz_0[j] = pb_x * tg_xxxx_xyyyyyz_0[j] + wp_x[j] * tg_xxxx_xyyyyyz_1[j] + 2.0 * fl1_fx * tg_xxx_xyyyyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xyyyyyz_1[j] + 0.5 * fl1_fxn * tg_xxxx_yyyyyz_1[j];

                    tg_xxxxx_xyyyyzz_0[j] = pb_x * tg_xxxx_xyyyyzz_0[j] + wp_x[j] * tg_xxxx_xyyyyzz_1[j] + 2.0 * fl1_fx * tg_xxx_xyyyyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xyyyyzz_1[j] + 0.5 * fl1_fxn * tg_xxxx_yyyyzz_1[j];

                    tg_xxxxx_xyyyzzz_0[j] = pb_x * tg_xxxx_xyyyzzz_0[j] + wp_x[j] * tg_xxxx_xyyyzzz_1[j] + 2.0 * fl1_fx * tg_xxx_xyyyzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xyyyzzz_1[j] + 0.5 * fl1_fxn * tg_xxxx_yyyzzz_1[j];

                    tg_xxxxx_xyyzzzz_0[j] = pb_x * tg_xxxx_xyyzzzz_0[j] + wp_x[j] * tg_xxxx_xyyzzzz_1[j] + 2.0 * fl1_fx * tg_xxx_xyyzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xyyzzzz_1[j] + 0.5 * fl1_fxn * tg_xxxx_yyzzzz_1[j];

                    tg_xxxxx_xyzzzzz_0[j] = pb_x * tg_xxxx_xyzzzzz_0[j] + wp_x[j] * tg_xxxx_xyzzzzz_1[j] + 2.0 * fl1_fx * tg_xxx_xyzzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_xxxx_yzzzzz_1[j];

                    tg_xxxxx_xzzzzzz_0[j] = pb_x * tg_xxxx_xzzzzzz_0[j] + wp_x[j] * tg_xxxx_xzzzzzz_1[j] + 2.0 * fl1_fx * tg_xxx_xzzzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xxxx_zzzzzz_1[j];

                    tg_xxxxx_yyyyyyy_0[j] = pb_x * tg_xxxx_yyyyyyy_0[j] + wp_x[j] * tg_xxxx_yyyyyyy_1[j] + 2.0 * fl1_fx * tg_xxx_yyyyyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_yyyyyyy_1[j];

                    tg_xxxxx_yyyyyyz_0[j] = pb_x * tg_xxxx_yyyyyyz_0[j] + wp_x[j] * tg_xxxx_yyyyyyz_1[j] + 2.0 * fl1_fx * tg_xxx_yyyyyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_yyyyyyz_1[j];

                    tg_xxxxx_yyyyyzz_0[j] = pb_x * tg_xxxx_yyyyyzz_0[j] + wp_x[j] * tg_xxxx_yyyyyzz_1[j] + 2.0 * fl1_fx * tg_xxx_yyyyyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_yyyyyzz_1[j];

                    tg_xxxxx_yyyyzzz_0[j] = pb_x * tg_xxxx_yyyyzzz_0[j] + wp_x[j] * tg_xxxx_yyyyzzz_1[j] + 2.0 * fl1_fx * tg_xxx_yyyyzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_yyyyzzz_1[j];

                    tg_xxxxx_yyyzzzz_0[j] = pb_x * tg_xxxx_yyyzzzz_0[j] + wp_x[j] * tg_xxxx_yyyzzzz_1[j] + 2.0 * fl1_fx * tg_xxx_yyyzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_yyyzzzz_1[j];

                    tg_xxxxx_yyzzzzz_0[j] = pb_x * tg_xxxx_yyzzzzz_0[j] + wp_x[j] * tg_xxxx_yyzzzzz_1[j] + 2.0 * fl1_fx * tg_xxx_yyzzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_yyzzzzz_1[j];

                    tg_xxxxx_yzzzzzz_0[j] = pb_x * tg_xxxx_yzzzzzz_0[j] + wp_x[j] * tg_xxxx_yzzzzzz_1[j] + 2.0 * fl1_fx * tg_xxx_yzzzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_yzzzzzz_1[j];

                    tg_xxxxx_zzzzzzz_0[j] = pb_x * tg_xxxx_zzzzzzz_0[j] + wp_x[j] * tg_xxxx_zzzzzzz_1[j] + 2.0 * fl1_fx * tg_xxx_zzzzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_zzzzzzz_1[j];

                    tg_xxxxy_xxxxxxx_0[j] = pb_x * tg_xxxy_xxxxxxx_0[j] + wp_x[j] * tg_xxxy_xxxxxxx_1[j] + 1.5 * fl1_fx * tg_xxy_xxxxxxx_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxxxxxx_1[j] + 3.5 * fl1_fxn * tg_xxxy_xxxxxx_1[j];

                    tg_xxxxy_xxxxxxy_0[j] = pb_x * tg_xxxy_xxxxxxy_0[j] + wp_x[j] * tg_xxxy_xxxxxxy_1[j] + 1.5 * fl1_fx * tg_xxy_xxxxxxy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxxxxxy_1[j] + 3.0 * fl1_fxn * tg_xxxy_xxxxxy_1[j];

                    tg_xxxxy_xxxxxxz_0[j] = pb_x * tg_xxxy_xxxxxxz_0[j] + wp_x[j] * tg_xxxy_xxxxxxz_1[j] + 1.5 * fl1_fx * tg_xxy_xxxxxxz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxxxxxz_1[j] + 3.0 * fl1_fxn * tg_xxxy_xxxxxz_1[j];

                    tg_xxxxy_xxxxxyy_0[j] = pb_x * tg_xxxy_xxxxxyy_0[j] + wp_x[j] * tg_xxxy_xxxxxyy_1[j] + 1.5 * fl1_fx * tg_xxy_xxxxxyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxxxxyy_1[j] + 2.5 * fl1_fxn * tg_xxxy_xxxxyy_1[j];

                    tg_xxxxy_xxxxxyz_0[j] = pb_x * tg_xxxy_xxxxxyz_0[j] + wp_x[j] * tg_xxxy_xxxxxyz_1[j] + 1.5 * fl1_fx * tg_xxy_xxxxxyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxxxxyz_1[j] + 2.5 * fl1_fxn * tg_xxxy_xxxxyz_1[j];

                    tg_xxxxy_xxxxxzz_0[j] = pb_x * tg_xxxy_xxxxxzz_0[j] + wp_x[j] * tg_xxxy_xxxxxzz_1[j] + 1.5 * fl1_fx * tg_xxy_xxxxxzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxxxxzz_1[j] + 2.5 * fl1_fxn * tg_xxxy_xxxxzz_1[j];

                    tg_xxxxy_xxxxyyy_0[j] = pb_x * tg_xxxy_xxxxyyy_0[j] + wp_x[j] * tg_xxxy_xxxxyyy_1[j] + 1.5 * fl1_fx * tg_xxy_xxxxyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxxxyyy_1[j] + 2.0 * fl1_fxn * tg_xxxy_xxxyyy_1[j];

                    tg_xxxxy_xxxxyyz_0[j] = pb_x * tg_xxxy_xxxxyyz_0[j] + wp_x[j] * tg_xxxy_xxxxyyz_1[j] + 1.5 * fl1_fx * tg_xxy_xxxxyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxxxyyz_1[j] + 2.0 * fl1_fxn * tg_xxxy_xxxyyz_1[j];

                    tg_xxxxy_xxxxyzz_0[j] = pb_x * tg_xxxy_xxxxyzz_0[j] + wp_x[j] * tg_xxxy_xxxxyzz_1[j] + 1.5 * fl1_fx * tg_xxy_xxxxyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxxxyzz_1[j] + 2.0 * fl1_fxn * tg_xxxy_xxxyzz_1[j];

                    tg_xxxxy_xxxxzzz_0[j] = pb_x * tg_xxxy_xxxxzzz_0[j] + wp_x[j] * tg_xxxy_xxxxzzz_1[j] + 1.5 * fl1_fx * tg_xxy_xxxxzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxxxzzz_1[j] + 2.0 * fl1_fxn * tg_xxxy_xxxzzz_1[j];

                    tg_xxxxy_xxxyyyy_0[j] = pb_x * tg_xxxy_xxxyyyy_0[j] + wp_x[j] * tg_xxxy_xxxyyyy_1[j] + 1.5 * fl1_fx * tg_xxy_xxxyyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxxyyyy_1[j] + 1.5 * fl1_fxn * tg_xxxy_xxyyyy_1[j];

                    tg_xxxxy_xxxyyyz_0[j] = pb_x * tg_xxxy_xxxyyyz_0[j] + wp_x[j] * tg_xxxy_xxxyyyz_1[j] + 1.5 * fl1_fx * tg_xxy_xxxyyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_xxxy_xxyyyz_1[j];

                    tg_xxxxy_xxxyyzz_0[j] = pb_x * tg_xxxy_xxxyyzz_0[j] + wp_x[j] * tg_xxxy_xxxyyzz_1[j] + 1.5 * fl1_fx * tg_xxy_xxxyyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxxyyzz_1[j] + 1.5 * fl1_fxn * tg_xxxy_xxyyzz_1[j];

                    tg_xxxxy_xxxyzzz_0[j] = pb_x * tg_xxxy_xxxyzzz_0[j] + wp_x[j] * tg_xxxy_xxxyzzz_1[j] + 1.5 * fl1_fx * tg_xxy_xxxyzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxxyzzz_1[j] + 1.5 * fl1_fxn * tg_xxxy_xxyzzz_1[j];

                    tg_xxxxy_xxxzzzz_0[j] = pb_x * tg_xxxy_xxxzzzz_0[j] + wp_x[j] * tg_xxxy_xxxzzzz_1[j] + 1.5 * fl1_fx * tg_xxy_xxxzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxxzzzz_1[j] + 1.5 * fl1_fxn * tg_xxxy_xxzzzz_1[j];

                    tg_xxxxy_xxyyyyy_0[j] = pb_x * tg_xxxy_xxyyyyy_0[j] + wp_x[j] * tg_xxxy_xxyyyyy_1[j] + 1.5 * fl1_fx * tg_xxy_xxyyyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxyyyyy_1[j] + fl1_fxn * tg_xxxy_xyyyyy_1[j];

                    tg_xxxxy_xxyyyyz_0[j] = pb_x * tg_xxxy_xxyyyyz_0[j] + wp_x[j] * tg_xxxy_xxyyyyz_1[j] + 1.5 * fl1_fx * tg_xxy_xxyyyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxyyyyz_1[j] + fl1_fxn * tg_xxxy_xyyyyz_1[j];

                    tg_xxxxy_xxyyyzz_0[j] = pb_x * tg_xxxy_xxyyyzz_0[j] + wp_x[j] * tg_xxxy_xxyyyzz_1[j] + 1.5 * fl1_fx * tg_xxy_xxyyyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxyyyzz_1[j] + fl1_fxn * tg_xxxy_xyyyzz_1[j];

                    tg_xxxxy_xxyyzzz_0[j] = pb_x * tg_xxxy_xxyyzzz_0[j] + wp_x[j] * tg_xxxy_xxyyzzz_1[j] + 1.5 * fl1_fx * tg_xxy_xxyyzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxyyzzz_1[j] + fl1_fxn * tg_xxxy_xyyzzz_1[j];

                    tg_xxxxy_xxyzzzz_0[j] = pb_x * tg_xxxy_xxyzzzz_0[j] + wp_x[j] * tg_xxxy_xxyzzzz_1[j] + 1.5 * fl1_fx * tg_xxy_xxyzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxyzzzz_1[j] + fl1_fxn * tg_xxxy_xyzzzz_1[j];

                    tg_xxxxy_xxzzzzz_0[j] = pb_x * tg_xxxy_xxzzzzz_0[j] + wp_x[j] * tg_xxxy_xxzzzzz_1[j] + 1.5 * fl1_fx * tg_xxy_xxzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxzzzzz_1[j] + fl1_fxn * tg_xxxy_xzzzzz_1[j];

                    tg_xxxxy_xyyyyyy_0[j] = pb_x * tg_xxxy_xyyyyyy_0[j] + wp_x[j] * tg_xxxy_xyyyyyy_1[j] + 1.5 * fl1_fx * tg_xxy_xyyyyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xyyyyyy_1[j] + 0.5 * fl1_fxn * tg_xxxy_yyyyyy_1[j];

                    tg_xxxxy_xyyyyyz_0[j] = pb_x * tg_xxxy_xyyyyyz_0[j] + wp_x[j] * tg_xxxy_xyyyyyz_1[j] + 1.5 * fl1_fx * tg_xxy_xyyyyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xyyyyyz_1[j] + 0.5 * fl1_fxn * tg_xxxy_yyyyyz_1[j];

                    tg_xxxxy_xyyyyzz_0[j] = pb_x * tg_xxxy_xyyyyzz_0[j] + wp_x[j] * tg_xxxy_xyyyyzz_1[j] + 1.5 * fl1_fx * tg_xxy_xyyyyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xyyyyzz_1[j] + 0.5 * fl1_fxn * tg_xxxy_yyyyzz_1[j];

                    tg_xxxxy_xyyyzzz_0[j] = pb_x * tg_xxxy_xyyyzzz_0[j] + wp_x[j] * tg_xxxy_xyyyzzz_1[j] + 1.5 * fl1_fx * tg_xxy_xyyyzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xyyyzzz_1[j] + 0.5 * fl1_fxn * tg_xxxy_yyyzzz_1[j];

                    tg_xxxxy_xyyzzzz_0[j] = pb_x * tg_xxxy_xyyzzzz_0[j] + wp_x[j] * tg_xxxy_xyyzzzz_1[j] + 1.5 * fl1_fx * tg_xxy_xyyzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xyyzzzz_1[j] + 0.5 * fl1_fxn * tg_xxxy_yyzzzz_1[j];

                    tg_xxxxy_xyzzzzz_0[j] = pb_x * tg_xxxy_xyzzzzz_0[j] + wp_x[j] * tg_xxxy_xyzzzzz_1[j] + 1.5 * fl1_fx * tg_xxy_xyzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_xxxy_yzzzzz_1[j];

                    tg_xxxxy_xzzzzzz_0[j] = pb_x * tg_xxxy_xzzzzzz_0[j] + wp_x[j] * tg_xxxy_xzzzzzz_1[j] + 1.5 * fl1_fx * tg_xxy_xzzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xxxy_zzzzzz_1[j];

                    tg_xxxxy_yyyyyyy_0[j] = pb_x * tg_xxxy_yyyyyyy_0[j] + wp_x[j] * tg_xxxy_yyyyyyy_1[j] + 1.5 * fl1_fx * tg_xxy_yyyyyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_yyyyyyy_1[j];

                    tg_xxxxy_yyyyyyz_0[j] = pb_x * tg_xxxy_yyyyyyz_0[j] + wp_x[j] * tg_xxxy_yyyyyyz_1[j] + 1.5 * fl1_fx * tg_xxy_yyyyyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_yyyyyyz_1[j];

                    tg_xxxxy_yyyyyzz_0[j] = pb_x * tg_xxxy_yyyyyzz_0[j] + wp_x[j] * tg_xxxy_yyyyyzz_1[j] + 1.5 * fl1_fx * tg_xxy_yyyyyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_yyyyyzz_1[j];

                    tg_xxxxy_yyyyzzz_0[j] = pb_x * tg_xxxy_yyyyzzz_0[j] + wp_x[j] * tg_xxxy_yyyyzzz_1[j] + 1.5 * fl1_fx * tg_xxy_yyyyzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_yyyyzzz_1[j];

                    tg_xxxxy_yyyzzzz_0[j] = pb_x * tg_xxxy_yyyzzzz_0[j] + wp_x[j] * tg_xxxy_yyyzzzz_1[j] + 1.5 * fl1_fx * tg_xxy_yyyzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_yyyzzzz_1[j];

                    tg_xxxxy_yyzzzzz_0[j] = pb_x * tg_xxxy_yyzzzzz_0[j] + wp_x[j] * tg_xxxy_yyzzzzz_1[j] + 1.5 * fl1_fx * tg_xxy_yyzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_yyzzzzz_1[j];

                    tg_xxxxy_yzzzzzz_0[j] = pb_x * tg_xxxy_yzzzzzz_0[j] + wp_x[j] * tg_xxxy_yzzzzzz_1[j] + 1.5 * fl1_fx * tg_xxy_yzzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_yzzzzzz_1[j];

                    tg_xxxxy_zzzzzzz_0[j] = pb_x * tg_xxxy_zzzzzzz_0[j] + wp_x[j] * tg_xxxy_zzzzzzz_1[j] + 1.5 * fl1_fx * tg_xxy_zzzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_zzzzzzz_1[j];

                    tg_xxxxz_xxxxxxx_0[j] = pb_x * tg_xxxz_xxxxxxx_0[j] + wp_x[j] * tg_xxxz_xxxxxxx_1[j] + 1.5 * fl1_fx * tg_xxz_xxxxxxx_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxxxxxx_1[j] + 3.5 * fl1_fxn * tg_xxxz_xxxxxx_1[j];

                    tg_xxxxz_xxxxxxy_0[j] = pb_x * tg_xxxz_xxxxxxy_0[j] + wp_x[j] * tg_xxxz_xxxxxxy_1[j] + 1.5 * fl1_fx * tg_xxz_xxxxxxy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxxxxxy_1[j] + 3.0 * fl1_fxn * tg_xxxz_xxxxxy_1[j];

                    tg_xxxxz_xxxxxxz_0[j] = pb_x * tg_xxxz_xxxxxxz_0[j] + wp_x[j] * tg_xxxz_xxxxxxz_1[j] + 1.5 * fl1_fx * tg_xxz_xxxxxxz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxxxxxz_1[j] + 3.0 * fl1_fxn * tg_xxxz_xxxxxz_1[j];

                    tg_xxxxz_xxxxxyy_0[j] = pb_x * tg_xxxz_xxxxxyy_0[j] + wp_x[j] * tg_xxxz_xxxxxyy_1[j] + 1.5 * fl1_fx * tg_xxz_xxxxxyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxxxxyy_1[j] + 2.5 * fl1_fxn * tg_xxxz_xxxxyy_1[j];

                    tg_xxxxz_xxxxxyz_0[j] = pb_x * tg_xxxz_xxxxxyz_0[j] + wp_x[j] * tg_xxxz_xxxxxyz_1[j] + 1.5 * fl1_fx * tg_xxz_xxxxxyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxxxxyz_1[j] + 2.5 * fl1_fxn * tg_xxxz_xxxxyz_1[j];

                    tg_xxxxz_xxxxxzz_0[j] = pb_x * tg_xxxz_xxxxxzz_0[j] + wp_x[j] * tg_xxxz_xxxxxzz_1[j] + 1.5 * fl1_fx * tg_xxz_xxxxxzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxxxxzz_1[j] + 2.5 * fl1_fxn * tg_xxxz_xxxxzz_1[j];

                    tg_xxxxz_xxxxyyy_0[j] = pb_x * tg_xxxz_xxxxyyy_0[j] + wp_x[j] * tg_xxxz_xxxxyyy_1[j] + 1.5 * fl1_fx * tg_xxz_xxxxyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxxxyyy_1[j] + 2.0 * fl1_fxn * tg_xxxz_xxxyyy_1[j];

                    tg_xxxxz_xxxxyyz_0[j] = pb_x * tg_xxxz_xxxxyyz_0[j] + wp_x[j] * tg_xxxz_xxxxyyz_1[j] + 1.5 * fl1_fx * tg_xxz_xxxxyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxxxyyz_1[j] + 2.0 * fl1_fxn * tg_xxxz_xxxyyz_1[j];

                    tg_xxxxz_xxxxyzz_0[j] = pb_x * tg_xxxz_xxxxyzz_0[j] + wp_x[j] * tg_xxxz_xxxxyzz_1[j] + 1.5 * fl1_fx * tg_xxz_xxxxyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxxxyzz_1[j] + 2.0 * fl1_fxn * tg_xxxz_xxxyzz_1[j];

                    tg_xxxxz_xxxxzzz_0[j] = pb_x * tg_xxxz_xxxxzzz_0[j] + wp_x[j] * tg_xxxz_xxxxzzz_1[j] + 1.5 * fl1_fx * tg_xxz_xxxxzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxxxzzz_1[j] + 2.0 * fl1_fxn * tg_xxxz_xxxzzz_1[j];

                    tg_xxxxz_xxxyyyy_0[j] = pb_x * tg_xxxz_xxxyyyy_0[j] + wp_x[j] * tg_xxxz_xxxyyyy_1[j] + 1.5 * fl1_fx * tg_xxz_xxxyyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxxyyyy_1[j] + 1.5 * fl1_fxn * tg_xxxz_xxyyyy_1[j];

                    tg_xxxxz_xxxyyyz_0[j] = pb_x * tg_xxxz_xxxyyyz_0[j] + wp_x[j] * tg_xxxz_xxxyyyz_1[j] + 1.5 * fl1_fx * tg_xxz_xxxyyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_xxxz_xxyyyz_1[j];

                    tg_xxxxz_xxxyyzz_0[j] = pb_x * tg_xxxz_xxxyyzz_0[j] + wp_x[j] * tg_xxxz_xxxyyzz_1[j] + 1.5 * fl1_fx * tg_xxz_xxxyyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxxyyzz_1[j] + 1.5 * fl1_fxn * tg_xxxz_xxyyzz_1[j];

                    tg_xxxxz_xxxyzzz_0[j] = pb_x * tg_xxxz_xxxyzzz_0[j] + wp_x[j] * tg_xxxz_xxxyzzz_1[j] + 1.5 * fl1_fx * tg_xxz_xxxyzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxxyzzz_1[j] + 1.5 * fl1_fxn * tg_xxxz_xxyzzz_1[j];

                    tg_xxxxz_xxxzzzz_0[j] = pb_x * tg_xxxz_xxxzzzz_0[j] + wp_x[j] * tg_xxxz_xxxzzzz_1[j] + 1.5 * fl1_fx * tg_xxz_xxxzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxxzzzz_1[j] + 1.5 * fl1_fxn * tg_xxxz_xxzzzz_1[j];

                    tg_xxxxz_xxyyyyy_0[j] = pb_x * tg_xxxz_xxyyyyy_0[j] + wp_x[j] * tg_xxxz_xxyyyyy_1[j] + 1.5 * fl1_fx * tg_xxz_xxyyyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxyyyyy_1[j] + fl1_fxn * tg_xxxz_xyyyyy_1[j];

                    tg_xxxxz_xxyyyyz_0[j] = pb_x * tg_xxxz_xxyyyyz_0[j] + wp_x[j] * tg_xxxz_xxyyyyz_1[j] + 1.5 * fl1_fx * tg_xxz_xxyyyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxyyyyz_1[j] + fl1_fxn * tg_xxxz_xyyyyz_1[j];

                    tg_xxxxz_xxyyyzz_0[j] = pb_x * tg_xxxz_xxyyyzz_0[j] + wp_x[j] * tg_xxxz_xxyyyzz_1[j] + 1.5 * fl1_fx * tg_xxz_xxyyyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxyyyzz_1[j] + fl1_fxn * tg_xxxz_xyyyzz_1[j];

                    tg_xxxxz_xxyyzzz_0[j] = pb_x * tg_xxxz_xxyyzzz_0[j] + wp_x[j] * tg_xxxz_xxyyzzz_1[j] + 1.5 * fl1_fx * tg_xxz_xxyyzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxyyzzz_1[j] + fl1_fxn * tg_xxxz_xyyzzz_1[j];

                    tg_xxxxz_xxyzzzz_0[j] = pb_x * tg_xxxz_xxyzzzz_0[j] + wp_x[j] * tg_xxxz_xxyzzzz_1[j] + 1.5 * fl1_fx * tg_xxz_xxyzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxyzzzz_1[j] + fl1_fxn * tg_xxxz_xyzzzz_1[j];

                    tg_xxxxz_xxzzzzz_0[j] = pb_x * tg_xxxz_xxzzzzz_0[j] + wp_x[j] * tg_xxxz_xxzzzzz_1[j] + 1.5 * fl1_fx * tg_xxz_xxzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxzzzzz_1[j] + fl1_fxn * tg_xxxz_xzzzzz_1[j];

                    tg_xxxxz_xyyyyyy_0[j] = pb_x * tg_xxxz_xyyyyyy_0[j] + wp_x[j] * tg_xxxz_xyyyyyy_1[j] + 1.5 * fl1_fx * tg_xxz_xyyyyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xyyyyyy_1[j] + 0.5 * fl1_fxn * tg_xxxz_yyyyyy_1[j];

                    tg_xxxxz_xyyyyyz_0[j] = pb_x * tg_xxxz_xyyyyyz_0[j] + wp_x[j] * tg_xxxz_xyyyyyz_1[j] + 1.5 * fl1_fx * tg_xxz_xyyyyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xyyyyyz_1[j] + 0.5 * fl1_fxn * tg_xxxz_yyyyyz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSK_95_190(      CMemBlock2D<double>& primBuffer,
                                        const CRecursionMap&       recursionMap,
                                        const CMemBlock2D<double>& osFactors,
                                        const CMemBlock2D<double>& wpDistances,
                                        const CGtoPairsBlock&      braGtoPairsBlock,
                                        const CGtoPairsBlock&      ketGtoPairsBlock,
                                        const int32_t              nKetPrimPairs,
                                        const int32_t              iContrPair)
    {
        // Batch of Integrals (95,190)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {5, -1, -1, -1},
                                             {7, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_7_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_4_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_3_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_4_6_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {6, -1, -1, -1}, 
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

                auto tg_xxxz_xyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 95); 

                auto tg_xxxz_xyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 96); 

                auto tg_xxxz_xyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 97); 

                auto tg_xxxz_xyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 98); 

                auto tg_xxxz_xzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 99); 

                auto tg_xxxz_yyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 100); 

                auto tg_xxxz_yyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 101); 

                auto tg_xxxz_yyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 102); 

                auto tg_xxxz_yyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 103); 

                auto tg_xxxz_yyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 104); 

                auto tg_xxxz_yyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 105); 

                auto tg_xxxz_yzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 106); 

                auto tg_xxxz_zzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 107); 

                auto tg_xxyy_xxxxxxx_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 108); 

                auto tg_xxyy_xxxxxxy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 109); 

                auto tg_xxyy_xxxxxxz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 110); 

                auto tg_xxyy_xxxxxyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 111); 

                auto tg_xxyy_xxxxxyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 112); 

                auto tg_xxyy_xxxxxzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 113); 

                auto tg_xxyy_xxxxyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 114); 

                auto tg_xxyy_xxxxyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 115); 

                auto tg_xxyy_xxxxyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 116); 

                auto tg_xxyy_xxxxzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 117); 

                auto tg_xxyy_xxxyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 118); 

                auto tg_xxyy_xxxyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 119); 

                auto tg_xxyy_xxxyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 120); 

                auto tg_xxyy_xxxyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 121); 

                auto tg_xxyy_xxxzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 122); 

                auto tg_xxyy_xxyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 123); 

                auto tg_xxyy_xxyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 124); 

                auto tg_xxyy_xxyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 125); 

                auto tg_xxyy_xxyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 126); 

                auto tg_xxyy_xxyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 127); 

                auto tg_xxyy_xxzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 128); 

                auto tg_xxyy_xyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 129); 

                auto tg_xxyy_xyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 130); 

                auto tg_xxyy_xyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 131); 

                auto tg_xxyy_xyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 132); 

                auto tg_xxyy_xyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 133); 

                auto tg_xxyy_xyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 134); 

                auto tg_xxyy_xzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 135); 

                auto tg_xxyy_yyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 136); 

                auto tg_xxyy_yyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 137); 

                auto tg_xxyy_yyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 138); 

                auto tg_xxyy_yyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 139); 

                auto tg_xxyy_yyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 140); 

                auto tg_xxyy_yyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 141); 

                auto tg_xxyy_yzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 142); 

                auto tg_xxyy_zzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 143); 

                auto tg_xxyz_xxxxxxx_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 144); 

                auto tg_xxyz_xxxxxxy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 145); 

                auto tg_xxyz_xxxxxxz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 146); 

                auto tg_xxyz_xxxxxyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 147); 

                auto tg_xxyz_xxxxxyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 148); 

                auto tg_xxyz_xxxxxzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 149); 

                auto tg_xxyz_xxxxyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 150); 

                auto tg_xxyz_xxxxyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 151); 

                auto tg_xxyz_xxxxyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 152); 

                auto tg_xxyz_xxxxzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 153); 

                auto tg_xxyz_xxxyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 154); 

                auto tg_xxyz_xxxyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 155); 

                auto tg_xxyz_xxxyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 156); 

                auto tg_xxyz_xxxyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 157); 

                auto tg_xxyz_xxxzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 158); 

                auto tg_xxyz_xxyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 159); 

                auto tg_xxyz_xxyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 160); 

                auto tg_xxyz_xxyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 161); 

                auto tg_xxyz_xxyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 162); 

                auto tg_xxyz_xxyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 163); 

                auto tg_xxyz_xxzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 164); 

                auto tg_xxyz_xyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 165); 

                auto tg_xxyz_xyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 166); 

                auto tg_xxyz_xyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 167); 

                auto tg_xxyz_xyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 168); 

                auto tg_xxyz_xyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 169); 

                auto tg_xxyz_xyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 170); 

                auto tg_xxyz_xzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 171); 

                auto tg_xxyz_yyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 172); 

                auto tg_xxyz_yyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 173); 

                auto tg_xxyz_yyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 174); 

                auto tg_xxyz_yyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 175); 

                auto tg_xxyz_yyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 176); 

                auto tg_xxyz_yyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 177); 

                auto tg_xxyz_yzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 178); 

                auto tg_xxyz_zzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 179); 

                auto tg_xxzz_xxxxxxx_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 180); 

                auto tg_xxzz_xxxxxxy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 181); 

                auto tg_xxzz_xxxxxxz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 182); 

                auto tg_xxzz_xxxxxyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 183); 

                auto tg_xxzz_xxxxxyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 184); 

                auto tg_xxzz_xxxxxzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 185); 

                auto tg_xxzz_xxxxyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 186); 

                auto tg_xxzz_xxxxyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 187); 

                auto tg_xxzz_xxxxyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 188); 

                auto tg_xxzz_xxxxzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 189); 

                auto tg_xxxz_xyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 95); 

                auto tg_xxxz_xyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 96); 

                auto tg_xxxz_xyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 97); 

                auto tg_xxxz_xyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 98); 

                auto tg_xxxz_xzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 99); 

                auto tg_xxxz_yyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 100); 

                auto tg_xxxz_yyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 101); 

                auto tg_xxxz_yyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 102); 

                auto tg_xxxz_yyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 103); 

                auto tg_xxxz_yyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 104); 

                auto tg_xxxz_yyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 105); 

                auto tg_xxxz_yzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 106); 

                auto tg_xxxz_zzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 107); 

                auto tg_xxyy_xxxxxxx_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 108); 

                auto tg_xxyy_xxxxxxy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 109); 

                auto tg_xxyy_xxxxxxz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 110); 

                auto tg_xxyy_xxxxxyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 111); 

                auto tg_xxyy_xxxxxyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 112); 

                auto tg_xxyy_xxxxxzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 113); 

                auto tg_xxyy_xxxxyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 114); 

                auto tg_xxyy_xxxxyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 115); 

                auto tg_xxyy_xxxxyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 116); 

                auto tg_xxyy_xxxxzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 117); 

                auto tg_xxyy_xxxyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 118); 

                auto tg_xxyy_xxxyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 119); 

                auto tg_xxyy_xxxyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 120); 

                auto tg_xxyy_xxxyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 121); 

                auto tg_xxyy_xxxzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 122); 

                auto tg_xxyy_xxyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 123); 

                auto tg_xxyy_xxyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 124); 

                auto tg_xxyy_xxyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 125); 

                auto tg_xxyy_xxyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 126); 

                auto tg_xxyy_xxyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 127); 

                auto tg_xxyy_xxzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 128); 

                auto tg_xxyy_xyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 129); 

                auto tg_xxyy_xyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 130); 

                auto tg_xxyy_xyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 131); 

                auto tg_xxyy_xyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 132); 

                auto tg_xxyy_xyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 133); 

                auto tg_xxyy_xyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 134); 

                auto tg_xxyy_xzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 135); 

                auto tg_xxyy_yyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 136); 

                auto tg_xxyy_yyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 137); 

                auto tg_xxyy_yyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 138); 

                auto tg_xxyy_yyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 139); 

                auto tg_xxyy_yyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 140); 

                auto tg_xxyy_yyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 141); 

                auto tg_xxyy_yzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 142); 

                auto tg_xxyy_zzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 143); 

                auto tg_xxyz_xxxxxxx_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 144); 

                auto tg_xxyz_xxxxxxy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 145); 

                auto tg_xxyz_xxxxxxz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 146); 

                auto tg_xxyz_xxxxxyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 147); 

                auto tg_xxyz_xxxxxyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 148); 

                auto tg_xxyz_xxxxxzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 149); 

                auto tg_xxyz_xxxxyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 150); 

                auto tg_xxyz_xxxxyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 151); 

                auto tg_xxyz_xxxxyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 152); 

                auto tg_xxyz_xxxxzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 153); 

                auto tg_xxyz_xxxyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 154); 

                auto tg_xxyz_xxxyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 155); 

                auto tg_xxyz_xxxyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 156); 

                auto tg_xxyz_xxxyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 157); 

                auto tg_xxyz_xxxzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 158); 

                auto tg_xxyz_xxyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 159); 

                auto tg_xxyz_xxyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 160); 

                auto tg_xxyz_xxyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 161); 

                auto tg_xxyz_xxyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 162); 

                auto tg_xxyz_xxyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 163); 

                auto tg_xxyz_xxzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 164); 

                auto tg_xxyz_xyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 165); 

                auto tg_xxyz_xyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 166); 

                auto tg_xxyz_xyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 167); 

                auto tg_xxyz_xyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 168); 

                auto tg_xxyz_xyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 169); 

                auto tg_xxyz_xyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 170); 

                auto tg_xxyz_xzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 171); 

                auto tg_xxyz_yyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 172); 

                auto tg_xxyz_yyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 173); 

                auto tg_xxyz_yyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 174); 

                auto tg_xxyz_yyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 175); 

                auto tg_xxyz_yyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 176); 

                auto tg_xxyz_yyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 177); 

                auto tg_xxyz_yzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 178); 

                auto tg_xxyz_zzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 179); 

                auto tg_xxzz_xxxxxxx_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 180); 

                auto tg_xxzz_xxxxxxy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 181); 

                auto tg_xxzz_xxxxxxz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 182); 

                auto tg_xxzz_xxxxxyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 183); 

                auto tg_xxzz_xxxxxyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 184); 

                auto tg_xxzz_xxxxxzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 185); 

                auto tg_xxzz_xxxxyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 186); 

                auto tg_xxzz_xxxxyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 187); 

                auto tg_xxzz_xxxxyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 188); 

                auto tg_xxzz_xxxxzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 189); 

                auto tg_xxz_xyyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 95); 

                auto tg_xxz_xyyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 96); 

                auto tg_xxz_xyyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 97); 

                auto tg_xxz_xyzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 98); 

                auto tg_xxz_xzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 99); 

                auto tg_xxz_yyyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 100); 

                auto tg_xxz_yyyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 101); 

                auto tg_xxz_yyyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 102); 

                auto tg_xxz_yyyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 103); 

                auto tg_xxz_yyyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 104); 

                auto tg_xxz_yyzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 105); 

                auto tg_xxz_yzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 106); 

                auto tg_xxz_zzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 107); 

                auto tg_xyy_xxxxxxx_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 108); 

                auto tg_xyy_xxxxxxy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 109); 

                auto tg_xyy_xxxxxxz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 110); 

                auto tg_xyy_xxxxxyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 111); 

                auto tg_xyy_xxxxxyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 112); 

                auto tg_xyy_xxxxxzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 113); 

                auto tg_xyy_xxxxyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 114); 

                auto tg_xyy_xxxxyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 115); 

                auto tg_xyy_xxxxyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 116); 

                auto tg_xyy_xxxxzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 117); 

                auto tg_xyy_xxxyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 118); 

                auto tg_xyy_xxxyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 119); 

                auto tg_xyy_xxxyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 120); 

                auto tg_xyy_xxxyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 121); 

                auto tg_xyy_xxxzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 122); 

                auto tg_xyy_xxyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 123); 

                auto tg_xyy_xxyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 124); 

                auto tg_xyy_xxyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 125); 

                auto tg_xyy_xxyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 126); 

                auto tg_xyy_xxyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 127); 

                auto tg_xyy_xxzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 128); 

                auto tg_xyy_xyyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 129); 

                auto tg_xyy_xyyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 130); 

                auto tg_xyy_xyyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 131); 

                auto tg_xyy_xyyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 132); 

                auto tg_xyy_xyyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 133); 

                auto tg_xyy_xyzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 134); 

                auto tg_xyy_xzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 135); 

                auto tg_xyy_yyyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 136); 

                auto tg_xyy_yyyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 137); 

                auto tg_xyy_yyyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 138); 

                auto tg_xyy_yyyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 139); 

                auto tg_xyy_yyyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 140); 

                auto tg_xyy_yyzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 141); 

                auto tg_xyy_yzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 142); 

                auto tg_xyy_zzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 143); 

                auto tg_xyz_xxxxxxx_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 144); 

                auto tg_xyz_xxxxxxy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 145); 

                auto tg_xyz_xxxxxxz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 146); 

                auto tg_xyz_xxxxxyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 147); 

                auto tg_xyz_xxxxxyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 148); 

                auto tg_xyz_xxxxxzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 149); 

                auto tg_xyz_xxxxyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 150); 

                auto tg_xyz_xxxxyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 151); 

                auto tg_xyz_xxxxyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 152); 

                auto tg_xyz_xxxxzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 153); 

                auto tg_xyz_xxxyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 154); 

                auto tg_xyz_xxxyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 155); 

                auto tg_xyz_xxxyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 156); 

                auto tg_xyz_xxxyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 157); 

                auto tg_xyz_xxxzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 158); 

                auto tg_xyz_xxyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 159); 

                auto tg_xyz_xxyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 160); 

                auto tg_xyz_xxyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 161); 

                auto tg_xyz_xxyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 162); 

                auto tg_xyz_xxyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 163); 

                auto tg_xyz_xxzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 164); 

                auto tg_xyz_xyyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 165); 

                auto tg_xyz_xyyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 166); 

                auto tg_xyz_xyyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 167); 

                auto tg_xyz_xyyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 168); 

                auto tg_xyz_xyyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 169); 

                auto tg_xyz_xyzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 170); 

                auto tg_xyz_xzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 171); 

                auto tg_xyz_yyyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 172); 

                auto tg_xyz_yyyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 173); 

                auto tg_xyz_yyyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 174); 

                auto tg_xyz_yyyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 175); 

                auto tg_xyz_yyyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 176); 

                auto tg_xyz_yyzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 177); 

                auto tg_xyz_yzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 178); 

                auto tg_xyz_zzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 179); 

                auto tg_xzz_xxxxxxx_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 180); 

                auto tg_xzz_xxxxxxy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 181); 

                auto tg_xzz_xxxxxxz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 182); 

                auto tg_xzz_xxxxxyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 183); 

                auto tg_xzz_xxxxxyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 184); 

                auto tg_xzz_xxxxxzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 185); 

                auto tg_xzz_xxxxyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 186); 

                auto tg_xzz_xxxxyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 187); 

                auto tg_xzz_xxxxyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 188); 

                auto tg_xzz_xxxxzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 189); 

                auto tg_xxz_xyyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 95); 

                auto tg_xxz_xyyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 96); 

                auto tg_xxz_xyyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 97); 

                auto tg_xxz_xyzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 98); 

                auto tg_xxz_xzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 99); 

                auto tg_xxz_yyyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 100); 

                auto tg_xxz_yyyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 101); 

                auto tg_xxz_yyyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 102); 

                auto tg_xxz_yyyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 103); 

                auto tg_xxz_yyyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 104); 

                auto tg_xxz_yyzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 105); 

                auto tg_xxz_yzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 106); 

                auto tg_xxz_zzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 107); 

                auto tg_xyy_xxxxxxx_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 108); 

                auto tg_xyy_xxxxxxy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 109); 

                auto tg_xyy_xxxxxxz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 110); 

                auto tg_xyy_xxxxxyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 111); 

                auto tg_xyy_xxxxxyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 112); 

                auto tg_xyy_xxxxxzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 113); 

                auto tg_xyy_xxxxyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 114); 

                auto tg_xyy_xxxxyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 115); 

                auto tg_xyy_xxxxyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 116); 

                auto tg_xyy_xxxxzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 117); 

                auto tg_xyy_xxxyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 118); 

                auto tg_xyy_xxxyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 119); 

                auto tg_xyy_xxxyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 120); 

                auto tg_xyy_xxxyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 121); 

                auto tg_xyy_xxxzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 122); 

                auto tg_xyy_xxyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 123); 

                auto tg_xyy_xxyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 124); 

                auto tg_xyy_xxyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 125); 

                auto tg_xyy_xxyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 126); 

                auto tg_xyy_xxyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 127); 

                auto tg_xyy_xxzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 128); 

                auto tg_xyy_xyyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 129); 

                auto tg_xyy_xyyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 130); 

                auto tg_xyy_xyyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 131); 

                auto tg_xyy_xyyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 132); 

                auto tg_xyy_xyyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 133); 

                auto tg_xyy_xyzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 134); 

                auto tg_xyy_xzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 135); 

                auto tg_xyy_yyyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 136); 

                auto tg_xyy_yyyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 137); 

                auto tg_xyy_yyyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 138); 

                auto tg_xyy_yyyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 139); 

                auto tg_xyy_yyyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 140); 

                auto tg_xyy_yyzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 141); 

                auto tg_xyy_yzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 142); 

                auto tg_xyy_zzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 143); 

                auto tg_xyz_xxxxxxx_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 144); 

                auto tg_xyz_xxxxxxy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 145); 

                auto tg_xyz_xxxxxxz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 146); 

                auto tg_xyz_xxxxxyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 147); 

                auto tg_xyz_xxxxxyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 148); 

                auto tg_xyz_xxxxxzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 149); 

                auto tg_xyz_xxxxyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 150); 

                auto tg_xyz_xxxxyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 151); 

                auto tg_xyz_xxxxyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 152); 

                auto tg_xyz_xxxxzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 153); 

                auto tg_xyz_xxxyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 154); 

                auto tg_xyz_xxxyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 155); 

                auto tg_xyz_xxxyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 156); 

                auto tg_xyz_xxxyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 157); 

                auto tg_xyz_xxxzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 158); 

                auto tg_xyz_xxyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 159); 

                auto tg_xyz_xxyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 160); 

                auto tg_xyz_xxyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 161); 

                auto tg_xyz_xxyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 162); 

                auto tg_xyz_xxyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 163); 

                auto tg_xyz_xxzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 164); 

                auto tg_xyz_xyyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 165); 

                auto tg_xyz_xyyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 166); 

                auto tg_xyz_xyyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 167); 

                auto tg_xyz_xyyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 168); 

                auto tg_xyz_xyyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 169); 

                auto tg_xyz_xyzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 170); 

                auto tg_xyz_xzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 171); 

                auto tg_xyz_yyyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 172); 

                auto tg_xyz_yyyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 173); 

                auto tg_xyz_yyyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 174); 

                auto tg_xyz_yyyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 175); 

                auto tg_xyz_yyyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 176); 

                auto tg_xyz_yyzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 177); 

                auto tg_xyz_yzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 178); 

                auto tg_xyz_zzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 179); 

                auto tg_xzz_xxxxxxx_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 180); 

                auto tg_xzz_xxxxxxy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 181); 

                auto tg_xzz_xxxxxxz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 182); 

                auto tg_xzz_xxxxxyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 183); 

                auto tg_xzz_xxxxxyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 184); 

                auto tg_xzz_xxxxxzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 185); 

                auto tg_xzz_xxxxyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 186); 

                auto tg_xzz_xxxxyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 187); 

                auto tg_xzz_xxxxyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 188); 

                auto tg_xzz_xxxxzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 189); 

                auto tg_xxxz_yyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 79); 

                auto tg_xxxz_yyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 80); 

                auto tg_xxxz_yyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 81); 

                auto tg_xxxz_yzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 82); 

                auto tg_xxxz_zzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 83); 

                auto tg_xxyy_xxxxxx_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 84); 

                auto tg_xxyy_xxxxxy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 85); 

                auto tg_xxyy_xxxxxz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 86); 

                auto tg_xxyy_xxxxyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 87); 

                auto tg_xxyy_xxxxyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 88); 

                auto tg_xxyy_xxxxzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 89); 

                auto tg_xxyy_xxxyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 90); 

                auto tg_xxyy_xxxyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 91); 

                auto tg_xxyy_xxxyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 92); 

                auto tg_xxyy_xxxzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 93); 

                auto tg_xxyy_xxyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 94); 

                auto tg_xxyy_xxyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 95); 

                auto tg_xxyy_xxyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 96); 

                auto tg_xxyy_xxyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 97); 

                auto tg_xxyy_xxzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 98); 

                auto tg_xxyy_xyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 99); 

                auto tg_xxyy_xyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 100); 

                auto tg_xxyy_xyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 101); 

                auto tg_xxyy_xyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 102); 

                auto tg_xxyy_xyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 103); 

                auto tg_xxyy_xzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 104); 

                auto tg_xxyy_yyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 105); 

                auto tg_xxyy_yyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 106); 

                auto tg_xxyy_yyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 107); 

                auto tg_xxyy_yyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 108); 

                auto tg_xxyy_yyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 109); 

                auto tg_xxyy_yzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 110); 

                auto tg_xxyy_zzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 111); 

                auto tg_xxyz_xxxxxx_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 112); 

                auto tg_xxyz_xxxxxy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 113); 

                auto tg_xxyz_xxxxxz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 114); 

                auto tg_xxyz_xxxxyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 115); 

                auto tg_xxyz_xxxxyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 116); 

                auto tg_xxyz_xxxxzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 117); 

                auto tg_xxyz_xxxyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 118); 

                auto tg_xxyz_xxxyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 119); 

                auto tg_xxyz_xxxyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 120); 

                auto tg_xxyz_xxxzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 121); 

                auto tg_xxyz_xxyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 122); 

                auto tg_xxyz_xxyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 123); 

                auto tg_xxyz_xxyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 124); 

                auto tg_xxyz_xxyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 125); 

                auto tg_xxyz_xxzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 126); 

                auto tg_xxyz_xyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 127); 

                auto tg_xxyz_xyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 128); 

                auto tg_xxyz_xyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 129); 

                auto tg_xxyz_xyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 130); 

                auto tg_xxyz_xyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 131); 

                auto tg_xxyz_xzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 132); 

                auto tg_xxyz_yyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 133); 

                auto tg_xxyz_yyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 134); 

                auto tg_xxyz_yyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 135); 

                auto tg_xxyz_yyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 136); 

                auto tg_xxyz_yyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 137); 

                auto tg_xxyz_yzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 138); 

                auto tg_xxyz_zzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 139); 

                auto tg_xxzz_xxxxxx_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 140); 

                auto tg_xxzz_xxxxxy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 141); 

                auto tg_xxzz_xxxxxz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 142); 

                auto tg_xxzz_xxxxyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 143); 

                auto tg_xxzz_xxxxyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 144); 

                auto tg_xxzz_xxxxzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 145); 

                auto tg_xxzz_xxxyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 146); 

                auto tg_xxzz_xxxyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 147); 

                auto tg_xxzz_xxxyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 148); 

                auto tg_xxzz_xxxzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 149); 

                // set up pointers to integrals

                auto tg_xxxxz_xyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 95); 

                auto tg_xxxxz_xyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 96); 

                auto tg_xxxxz_xyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 97); 

                auto tg_xxxxz_xyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 98); 

                auto tg_xxxxz_xzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 99); 

                auto tg_xxxxz_yyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 100); 

                auto tg_xxxxz_yyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 101); 

                auto tg_xxxxz_yyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 102); 

                auto tg_xxxxz_yyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 103); 

                auto tg_xxxxz_yyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 104); 

                auto tg_xxxxz_yyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 105); 

                auto tg_xxxxz_yzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 106); 

                auto tg_xxxxz_zzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 107); 

                auto tg_xxxyy_xxxxxxx_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 108); 

                auto tg_xxxyy_xxxxxxy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 109); 

                auto tg_xxxyy_xxxxxxz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 110); 

                auto tg_xxxyy_xxxxxyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 111); 

                auto tg_xxxyy_xxxxxyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 112); 

                auto tg_xxxyy_xxxxxzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 113); 

                auto tg_xxxyy_xxxxyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 114); 

                auto tg_xxxyy_xxxxyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 115); 

                auto tg_xxxyy_xxxxyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 116); 

                auto tg_xxxyy_xxxxzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 117); 

                auto tg_xxxyy_xxxyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 118); 

                auto tg_xxxyy_xxxyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 119); 

                auto tg_xxxyy_xxxyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 120); 

                auto tg_xxxyy_xxxyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 121); 

                auto tg_xxxyy_xxxzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 122); 

                auto tg_xxxyy_xxyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 123); 

                auto tg_xxxyy_xxyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 124); 

                auto tg_xxxyy_xxyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 125); 

                auto tg_xxxyy_xxyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 126); 

                auto tg_xxxyy_xxyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 127); 

                auto tg_xxxyy_xxzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 128); 

                auto tg_xxxyy_xyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 129); 

                auto tg_xxxyy_xyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 130); 

                auto tg_xxxyy_xyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 131); 

                auto tg_xxxyy_xyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 132); 

                auto tg_xxxyy_xyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 133); 

                auto tg_xxxyy_xyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 134); 

                auto tg_xxxyy_xzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 135); 

                auto tg_xxxyy_yyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 136); 

                auto tg_xxxyy_yyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 137); 

                auto tg_xxxyy_yyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 138); 

                auto tg_xxxyy_yyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 139); 

                auto tg_xxxyy_yyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 140); 

                auto tg_xxxyy_yyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 141); 

                auto tg_xxxyy_yzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 142); 

                auto tg_xxxyy_zzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 143); 

                auto tg_xxxyz_xxxxxxx_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 144); 

                auto tg_xxxyz_xxxxxxy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 145); 

                auto tg_xxxyz_xxxxxxz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 146); 

                auto tg_xxxyz_xxxxxyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 147); 

                auto tg_xxxyz_xxxxxyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 148); 

                auto tg_xxxyz_xxxxxzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 149); 

                auto tg_xxxyz_xxxxyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 150); 

                auto tg_xxxyz_xxxxyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 151); 

                auto tg_xxxyz_xxxxyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 152); 

                auto tg_xxxyz_xxxxzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 153); 

                auto tg_xxxyz_xxxyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 154); 

                auto tg_xxxyz_xxxyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 155); 

                auto tg_xxxyz_xxxyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 156); 

                auto tg_xxxyz_xxxyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 157); 

                auto tg_xxxyz_xxxzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 158); 

                auto tg_xxxyz_xxyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 159); 

                auto tg_xxxyz_xxyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 160); 

                auto tg_xxxyz_xxyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 161); 

                auto tg_xxxyz_xxyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 162); 

                auto tg_xxxyz_xxyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 163); 

                auto tg_xxxyz_xxzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 164); 

                auto tg_xxxyz_xyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 165); 

                auto tg_xxxyz_xyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 166); 

                auto tg_xxxyz_xyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 167); 

                auto tg_xxxyz_xyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 168); 

                auto tg_xxxyz_xyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 169); 

                auto tg_xxxyz_xyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 170); 

                auto tg_xxxyz_xzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 171); 

                auto tg_xxxyz_yyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 172); 

                auto tg_xxxyz_yyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 173); 

                auto tg_xxxyz_yyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 174); 

                auto tg_xxxyz_yyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 175); 

                auto tg_xxxyz_yyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 176); 

                auto tg_xxxyz_yyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 177); 

                auto tg_xxxyz_yzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 178); 

                auto tg_xxxyz_zzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 179); 

                auto tg_xxxzz_xxxxxxx_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 180); 

                auto tg_xxxzz_xxxxxxy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 181); 

                auto tg_xxxzz_xxxxxxz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 182); 

                auto tg_xxxzz_xxxxxyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 183); 

                auto tg_xxxzz_xxxxxyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 184); 

                auto tg_xxxzz_xxxxxzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 185); 

                auto tg_xxxzz_xxxxyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 186); 

                auto tg_xxxzz_xxxxyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 187); 

                auto tg_xxxzz_xxxxyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 188); 

                auto tg_xxxzz_xxxxzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 189); 

                // Batch of Integrals (95,190)

                #pragma omp simd aligned(fxn, fza, tg_xxxxz_xyyyyzz_0, tg_xxxxz_xyyyzzz_0, tg_xxxxz_xyyzzzz_0, \
                                         tg_xxxxz_xyzzzzz_0, tg_xxxxz_xzzzzzz_0, tg_xxxxz_yyyyyyy_0, tg_xxxxz_yyyyyyz_0, \
                                         tg_xxxxz_yyyyyzz_0, tg_xxxxz_yyyyzzz_0, tg_xxxxz_yyyzzzz_0, tg_xxxxz_yyzzzzz_0, \
                                         tg_xxxxz_yzzzzzz_0, tg_xxxxz_zzzzzzz_0, tg_xxxyy_xxxxxxx_0, tg_xxxyy_xxxxxxy_0, \
                                         tg_xxxyy_xxxxxxz_0, tg_xxxyy_xxxxxyy_0, tg_xxxyy_xxxxxyz_0, tg_xxxyy_xxxxxzz_0, \
                                         tg_xxxyy_xxxxyyy_0, tg_xxxyy_xxxxyyz_0, tg_xxxyy_xxxxyzz_0, tg_xxxyy_xxxxzzz_0, \
                                         tg_xxxyy_xxxyyyy_0, tg_xxxyy_xxxyyyz_0, tg_xxxyy_xxxyyzz_0, tg_xxxyy_xxxyzzz_0, \
                                         tg_xxxyy_xxxzzzz_0, tg_xxxyy_xxyyyyy_0, tg_xxxyy_xxyyyyz_0, tg_xxxyy_xxyyyzz_0, \
                                         tg_xxxyy_xxyyzzz_0, tg_xxxyy_xxyzzzz_0, tg_xxxyy_xxzzzzz_0, tg_xxxyy_xyyyyyy_0, \
                                         tg_xxxyy_xyyyyyz_0, tg_xxxyy_xyyyyzz_0, tg_xxxyy_xyyyzzz_0, tg_xxxyy_xyyzzzz_0, \
                                         tg_xxxyy_xyzzzzz_0, tg_xxxyy_xzzzzzz_0, tg_xxxyy_yyyyyyy_0, tg_xxxyy_yyyyyyz_0, \
                                         tg_xxxyy_yyyyyzz_0, tg_xxxyy_yyyyzzz_0, tg_xxxyy_yyyzzzz_0, tg_xxxyy_yyzzzzz_0, \
                                         tg_xxxyy_yzzzzzz_0, tg_xxxyy_zzzzzzz_0, tg_xxxyz_xxxxxxx_0, tg_xxxyz_xxxxxxy_0, \
                                         tg_xxxyz_xxxxxxz_0, tg_xxxyz_xxxxxyy_0, tg_xxxyz_xxxxxyz_0, tg_xxxyz_xxxxxzz_0, \
                                         tg_xxxyz_xxxxyyy_0, tg_xxxyz_xxxxyyz_0, tg_xxxyz_xxxxyzz_0, tg_xxxyz_xxxxzzz_0, \
                                         tg_xxxyz_xxxyyyy_0, tg_xxxyz_xxxyyyz_0, tg_xxxyz_xxxyyzz_0, tg_xxxyz_xxxyzzz_0, \
                                         tg_xxxyz_xxxzzzz_0, tg_xxxyz_xxyyyyy_0, tg_xxxyz_xxyyyyz_0, tg_xxxyz_xxyyyzz_0, \
                                         tg_xxxyz_xxyyzzz_0, tg_xxxyz_xxyzzzz_0, tg_xxxyz_xxzzzzz_0, tg_xxxyz_xyyyyyy_0, \
                                         tg_xxxyz_xyyyyyz_0, tg_xxxyz_xyyyyzz_0, tg_xxxyz_xyyyzzz_0, tg_xxxyz_xyyzzzz_0, \
                                         tg_xxxyz_xyzzzzz_0, tg_xxxyz_xzzzzzz_0, tg_xxxyz_yyyyyyy_0, tg_xxxyz_yyyyyyz_0, \
                                         tg_xxxyz_yyyyyzz_0, tg_xxxyz_yyyyzzz_0, tg_xxxyz_yyyzzzz_0, tg_xxxyz_yyzzzzz_0, \
                                         tg_xxxyz_yzzzzzz_0, tg_xxxyz_zzzzzzz_0, tg_xxxz_xyyyyzz_0, tg_xxxz_xyyyyzz_1, \
                                         tg_xxxz_xyyyzzz_0, tg_xxxz_xyyyzzz_1, tg_xxxz_xyyzzzz_0, tg_xxxz_xyyzzzz_1, \
                                         tg_xxxz_xyzzzzz_0, tg_xxxz_xyzzzzz_1, tg_xxxz_xzzzzzz_0, tg_xxxz_xzzzzzz_1, \
                                         tg_xxxz_yyyyyyy_0, tg_xxxz_yyyyyyy_1, tg_xxxz_yyyyyyz_0, tg_xxxz_yyyyyyz_1, \
                                         tg_xxxz_yyyyyzz_0, tg_xxxz_yyyyyzz_1, tg_xxxz_yyyyzz_1, tg_xxxz_yyyyzzz_0, \
                                         tg_xxxz_yyyyzzz_1, tg_xxxz_yyyzzz_1, tg_xxxz_yyyzzzz_0, tg_xxxz_yyyzzzz_1, \
                                         tg_xxxz_yyzzzz_1, tg_xxxz_yyzzzzz_0, tg_xxxz_yyzzzzz_1, tg_xxxz_yzzzzz_1, \
                                         tg_xxxz_yzzzzzz_0, tg_xxxz_yzzzzzz_1, tg_xxxz_zzzzzz_1, tg_xxxz_zzzzzzz_0, \
                                         tg_xxxz_zzzzzzz_1, tg_xxxzz_xxxxxxx_0, tg_xxxzz_xxxxxxy_0, tg_xxxzz_xxxxxxz_0, \
                                         tg_xxxzz_xxxxxyy_0, tg_xxxzz_xxxxxyz_0, tg_xxxzz_xxxxxzz_0, tg_xxxzz_xxxxyyy_0, \
                                         tg_xxxzz_xxxxyyz_0, tg_xxxzz_xxxxyzz_0, tg_xxxzz_xxxxzzz_0, tg_xxyy_xxxxxx_1, \
                                         tg_xxyy_xxxxxxx_0, tg_xxyy_xxxxxxx_1, tg_xxyy_xxxxxxy_0, tg_xxyy_xxxxxxy_1, \
                                         tg_xxyy_xxxxxxz_0, tg_xxyy_xxxxxxz_1, tg_xxyy_xxxxxy_1, tg_xxyy_xxxxxyy_0, \
                                         tg_xxyy_xxxxxyy_1, tg_xxyy_xxxxxyz_0, tg_xxyy_xxxxxyz_1, tg_xxyy_xxxxxz_1, \
                                         tg_xxyy_xxxxxzz_0, tg_xxyy_xxxxxzz_1, tg_xxyy_xxxxyy_1, tg_xxyy_xxxxyyy_0, \
                                         tg_xxyy_xxxxyyy_1, tg_xxyy_xxxxyyz_0, tg_xxyy_xxxxyyz_1, tg_xxyy_xxxxyz_1, \
                                         tg_xxyy_xxxxyzz_0, tg_xxyy_xxxxyzz_1, tg_xxyy_xxxxzz_1, tg_xxyy_xxxxzzz_0, \
                                         tg_xxyy_xxxxzzz_1, tg_xxyy_xxxyyy_1, tg_xxyy_xxxyyyy_0, tg_xxyy_xxxyyyy_1, \
                                         tg_xxyy_xxxyyyz_0, tg_xxyy_xxxyyyz_1, tg_xxyy_xxxyyz_1, tg_xxyy_xxxyyzz_0, \
                                         tg_xxyy_xxxyyzz_1, tg_xxyy_xxxyzz_1, tg_xxyy_xxxyzzz_0, tg_xxyy_xxxyzzz_1, \
                                         tg_xxyy_xxxzzz_1, tg_xxyy_xxxzzzz_0, tg_xxyy_xxxzzzz_1, tg_xxyy_xxyyyy_1, \
                                         tg_xxyy_xxyyyyy_0, tg_xxyy_xxyyyyy_1, tg_xxyy_xxyyyyz_0, tg_xxyy_xxyyyyz_1, \
                                         tg_xxyy_xxyyyz_1, tg_xxyy_xxyyyzz_0, tg_xxyy_xxyyyzz_1, tg_xxyy_xxyyzz_1, \
                                         tg_xxyy_xxyyzzz_0, tg_xxyy_xxyyzzz_1, tg_xxyy_xxyzzz_1, tg_xxyy_xxyzzzz_0, \
                                         tg_xxyy_xxyzzzz_1, tg_xxyy_xxzzzz_1, tg_xxyy_xxzzzzz_0, tg_xxyy_xxzzzzz_1, \
                                         tg_xxyy_xyyyyy_1, tg_xxyy_xyyyyyy_0, tg_xxyy_xyyyyyy_1, tg_xxyy_xyyyyyz_0, \
                                         tg_xxyy_xyyyyyz_1, tg_xxyy_xyyyyz_1, tg_xxyy_xyyyyzz_0, tg_xxyy_xyyyyzz_1, \
                                         tg_xxyy_xyyyzz_1, tg_xxyy_xyyyzzz_0, tg_xxyy_xyyyzzz_1, tg_xxyy_xyyzzz_1, \
                                         tg_xxyy_xyyzzzz_0, tg_xxyy_xyyzzzz_1, tg_xxyy_xyzzzz_1, tg_xxyy_xyzzzzz_0, \
                                         tg_xxyy_xyzzzzz_1, tg_xxyy_xzzzzz_1, tg_xxyy_xzzzzzz_0, tg_xxyy_xzzzzzz_1, \
                                         tg_xxyy_yyyyyy_1, tg_xxyy_yyyyyyy_0, tg_xxyy_yyyyyyy_1, tg_xxyy_yyyyyyz_0, \
                                         tg_xxyy_yyyyyyz_1, tg_xxyy_yyyyyz_1, tg_xxyy_yyyyyzz_0, tg_xxyy_yyyyyzz_1, \
                                         tg_xxyy_yyyyzz_1, tg_xxyy_yyyyzzz_0, tg_xxyy_yyyyzzz_1, tg_xxyy_yyyzzz_1, \
                                         tg_xxyy_yyyzzzz_0, tg_xxyy_yyyzzzz_1, tg_xxyy_yyzzzz_1, tg_xxyy_yyzzzzz_0, \
                                         tg_xxyy_yyzzzzz_1, tg_xxyy_yzzzzz_1, tg_xxyy_yzzzzzz_0, tg_xxyy_yzzzzzz_1, \
                                         tg_xxyy_zzzzzz_1, tg_xxyy_zzzzzzz_0, tg_xxyy_zzzzzzz_1, tg_xxyz_xxxxxx_1, \
                                         tg_xxyz_xxxxxxx_0, tg_xxyz_xxxxxxx_1, tg_xxyz_xxxxxxy_0, tg_xxyz_xxxxxxy_1, \
                                         tg_xxyz_xxxxxxz_0, tg_xxyz_xxxxxxz_1, tg_xxyz_xxxxxy_1, tg_xxyz_xxxxxyy_0, \
                                         tg_xxyz_xxxxxyy_1, tg_xxyz_xxxxxyz_0, tg_xxyz_xxxxxyz_1, tg_xxyz_xxxxxz_1, \
                                         tg_xxyz_xxxxxzz_0, tg_xxyz_xxxxxzz_1, tg_xxyz_xxxxyy_1, tg_xxyz_xxxxyyy_0, \
                                         tg_xxyz_xxxxyyy_1, tg_xxyz_xxxxyyz_0, tg_xxyz_xxxxyyz_1, tg_xxyz_xxxxyz_1, \
                                         tg_xxyz_xxxxyzz_0, tg_xxyz_xxxxyzz_1, tg_xxyz_xxxxzz_1, tg_xxyz_xxxxzzz_0, \
                                         tg_xxyz_xxxxzzz_1, tg_xxyz_xxxyyy_1, tg_xxyz_xxxyyyy_0, tg_xxyz_xxxyyyy_1, \
                                         tg_xxyz_xxxyyyz_0, tg_xxyz_xxxyyyz_1, tg_xxyz_xxxyyz_1, tg_xxyz_xxxyyzz_0, \
                                         tg_xxyz_xxxyyzz_1, tg_xxyz_xxxyzz_1, tg_xxyz_xxxyzzz_0, tg_xxyz_xxxyzzz_1, \
                                         tg_xxyz_xxxzzz_1, tg_xxyz_xxxzzzz_0, tg_xxyz_xxxzzzz_1, tg_xxyz_xxyyyy_1, \
                                         tg_xxyz_xxyyyyy_0, tg_xxyz_xxyyyyy_1, tg_xxyz_xxyyyyz_0, tg_xxyz_xxyyyyz_1, \
                                         tg_xxyz_xxyyyz_1, tg_xxyz_xxyyyzz_0, tg_xxyz_xxyyyzz_1, tg_xxyz_xxyyzz_1, \
                                         tg_xxyz_xxyyzzz_0, tg_xxyz_xxyyzzz_1, tg_xxyz_xxyzzz_1, tg_xxyz_xxyzzzz_0, \
                                         tg_xxyz_xxyzzzz_1, tg_xxyz_xxzzzz_1, tg_xxyz_xxzzzzz_0, tg_xxyz_xxzzzzz_1, \
                                         tg_xxyz_xyyyyy_1, tg_xxyz_xyyyyyy_0, tg_xxyz_xyyyyyy_1, tg_xxyz_xyyyyyz_0, \
                                         tg_xxyz_xyyyyyz_1, tg_xxyz_xyyyyz_1, tg_xxyz_xyyyyzz_0, tg_xxyz_xyyyyzz_1, \
                                         tg_xxyz_xyyyzz_1, tg_xxyz_xyyyzzz_0, tg_xxyz_xyyyzzz_1, tg_xxyz_xyyzzz_1, \
                                         tg_xxyz_xyyzzzz_0, tg_xxyz_xyyzzzz_1, tg_xxyz_xyzzzz_1, tg_xxyz_xyzzzzz_0, \
                                         tg_xxyz_xyzzzzz_1, tg_xxyz_xzzzzz_1, tg_xxyz_xzzzzzz_0, tg_xxyz_xzzzzzz_1, \
                                         tg_xxyz_yyyyyy_1, tg_xxyz_yyyyyyy_0, tg_xxyz_yyyyyyy_1, tg_xxyz_yyyyyyz_0, \
                                         tg_xxyz_yyyyyyz_1, tg_xxyz_yyyyyz_1, tg_xxyz_yyyyyzz_0, tg_xxyz_yyyyyzz_1, \
                                         tg_xxyz_yyyyzz_1, tg_xxyz_yyyyzzz_0, tg_xxyz_yyyyzzz_1, tg_xxyz_yyyzzz_1, \
                                         tg_xxyz_yyyzzzz_0, tg_xxyz_yyyzzzz_1, tg_xxyz_yyzzzz_1, tg_xxyz_yyzzzzz_0, \
                                         tg_xxyz_yyzzzzz_1, tg_xxyz_yzzzzz_1, tg_xxyz_yzzzzzz_0, tg_xxyz_yzzzzzz_1, \
                                         tg_xxyz_zzzzzz_1, tg_xxyz_zzzzzzz_0, tg_xxyz_zzzzzzz_1, tg_xxz_xyyyyzz_0, \
                                         tg_xxz_xyyyyzz_1, tg_xxz_xyyyzzz_0, tg_xxz_xyyyzzz_1, tg_xxz_xyyzzzz_0, \
                                         tg_xxz_xyyzzzz_1, tg_xxz_xyzzzzz_0, tg_xxz_xyzzzzz_1, tg_xxz_xzzzzzz_0, \
                                         tg_xxz_xzzzzzz_1, tg_xxz_yyyyyyy_0, tg_xxz_yyyyyyy_1, tg_xxz_yyyyyyz_0, \
                                         tg_xxz_yyyyyyz_1, tg_xxz_yyyyyzz_0, tg_xxz_yyyyyzz_1, tg_xxz_yyyyzzz_0, \
                                         tg_xxz_yyyyzzz_1, tg_xxz_yyyzzzz_0, tg_xxz_yyyzzzz_1, tg_xxz_yyzzzzz_0, \
                                         tg_xxz_yyzzzzz_1, tg_xxz_yzzzzzz_0, tg_xxz_yzzzzzz_1, tg_xxz_zzzzzzz_0, \
                                         tg_xxz_zzzzzzz_1, tg_xxzz_xxxxxx_1, tg_xxzz_xxxxxxx_0, tg_xxzz_xxxxxxx_1, \
                                         tg_xxzz_xxxxxxy_0, tg_xxzz_xxxxxxy_1, tg_xxzz_xxxxxxz_0, tg_xxzz_xxxxxxz_1, \
                                         tg_xxzz_xxxxxy_1, tg_xxzz_xxxxxyy_0, tg_xxzz_xxxxxyy_1, tg_xxzz_xxxxxyz_0, \
                                         tg_xxzz_xxxxxyz_1, tg_xxzz_xxxxxz_1, tg_xxzz_xxxxxzz_0, tg_xxzz_xxxxxzz_1, \
                                         tg_xxzz_xxxxyy_1, tg_xxzz_xxxxyyy_0, tg_xxzz_xxxxyyy_1, tg_xxzz_xxxxyyz_0, \
                                         tg_xxzz_xxxxyyz_1, tg_xxzz_xxxxyz_1, tg_xxzz_xxxxyzz_0, tg_xxzz_xxxxyzz_1, \
                                         tg_xxzz_xxxxzz_1, tg_xxzz_xxxxzzz_0, tg_xxzz_xxxxzzz_1, tg_xxzz_xxxyyy_1, \
                                         tg_xxzz_xxxyyz_1, tg_xxzz_xxxyzz_1, tg_xxzz_xxxzzz_1, tg_xyy_xxxxxxx_0, \
                                         tg_xyy_xxxxxxx_1, tg_xyy_xxxxxxy_0, tg_xyy_xxxxxxy_1, tg_xyy_xxxxxxz_0, \
                                         tg_xyy_xxxxxxz_1, tg_xyy_xxxxxyy_0, tg_xyy_xxxxxyy_1, tg_xyy_xxxxxyz_0, \
                                         tg_xyy_xxxxxyz_1, tg_xyy_xxxxxzz_0, tg_xyy_xxxxxzz_1, tg_xyy_xxxxyyy_0, \
                                         tg_xyy_xxxxyyy_1, tg_xyy_xxxxyyz_0, tg_xyy_xxxxyyz_1, tg_xyy_xxxxyzz_0, \
                                         tg_xyy_xxxxyzz_1, tg_xyy_xxxxzzz_0, tg_xyy_xxxxzzz_1, tg_xyy_xxxyyyy_0, \
                                         tg_xyy_xxxyyyy_1, tg_xyy_xxxyyyz_0, tg_xyy_xxxyyyz_1, tg_xyy_xxxyyzz_0, \
                                         tg_xyy_xxxyyzz_1, tg_xyy_xxxyzzz_0, tg_xyy_xxxyzzz_1, tg_xyy_xxxzzzz_0, \
                                         tg_xyy_xxxzzzz_1, tg_xyy_xxyyyyy_0, tg_xyy_xxyyyyy_1, tg_xyy_xxyyyyz_0, \
                                         tg_xyy_xxyyyyz_1, tg_xyy_xxyyyzz_0, tg_xyy_xxyyyzz_1, tg_xyy_xxyyzzz_0, \
                                         tg_xyy_xxyyzzz_1, tg_xyy_xxyzzzz_0, tg_xyy_xxyzzzz_1, tg_xyy_xxzzzzz_0, \
                                         tg_xyy_xxzzzzz_1, tg_xyy_xyyyyyy_0, tg_xyy_xyyyyyy_1, tg_xyy_xyyyyyz_0, \
                                         tg_xyy_xyyyyyz_1, tg_xyy_xyyyyzz_0, tg_xyy_xyyyyzz_1, tg_xyy_xyyyzzz_0, \
                                         tg_xyy_xyyyzzz_1, tg_xyy_xyyzzzz_0, tg_xyy_xyyzzzz_1, tg_xyy_xyzzzzz_0, \
                                         tg_xyy_xyzzzzz_1, tg_xyy_xzzzzzz_0, tg_xyy_xzzzzzz_1, tg_xyy_yyyyyyy_0, \
                                         tg_xyy_yyyyyyy_1, tg_xyy_yyyyyyz_0, tg_xyy_yyyyyyz_1, tg_xyy_yyyyyzz_0, \
                                         tg_xyy_yyyyyzz_1, tg_xyy_yyyyzzz_0, tg_xyy_yyyyzzz_1, tg_xyy_yyyzzzz_0, \
                                         tg_xyy_yyyzzzz_1, tg_xyy_yyzzzzz_0, tg_xyy_yyzzzzz_1, tg_xyy_yzzzzzz_0, \
                                         tg_xyy_yzzzzzz_1, tg_xyy_zzzzzzz_0, tg_xyy_zzzzzzz_1, tg_xyz_xxxxxxx_0, \
                                         tg_xyz_xxxxxxx_1, tg_xyz_xxxxxxy_0, tg_xyz_xxxxxxy_1, tg_xyz_xxxxxxz_0, \
                                         tg_xyz_xxxxxxz_1, tg_xyz_xxxxxyy_0, tg_xyz_xxxxxyy_1, tg_xyz_xxxxxyz_0, \
                                         tg_xyz_xxxxxyz_1, tg_xyz_xxxxxzz_0, tg_xyz_xxxxxzz_1, tg_xyz_xxxxyyy_0, \
                                         tg_xyz_xxxxyyy_1, tg_xyz_xxxxyyz_0, tg_xyz_xxxxyyz_1, tg_xyz_xxxxyzz_0, \
                                         tg_xyz_xxxxyzz_1, tg_xyz_xxxxzzz_0, tg_xyz_xxxxzzz_1, tg_xyz_xxxyyyy_0, \
                                         tg_xyz_xxxyyyy_1, tg_xyz_xxxyyyz_0, tg_xyz_xxxyyyz_1, tg_xyz_xxxyyzz_0, \
                                         tg_xyz_xxxyyzz_1, tg_xyz_xxxyzzz_0, tg_xyz_xxxyzzz_1, tg_xyz_xxxzzzz_0, \
                                         tg_xyz_xxxzzzz_1, tg_xyz_xxyyyyy_0, tg_xyz_xxyyyyy_1, tg_xyz_xxyyyyz_0, \
                                         tg_xyz_xxyyyyz_1, tg_xyz_xxyyyzz_0, tg_xyz_xxyyyzz_1, tg_xyz_xxyyzzz_0, \
                                         tg_xyz_xxyyzzz_1, tg_xyz_xxyzzzz_0, tg_xyz_xxyzzzz_1, tg_xyz_xxzzzzz_0, \
                                         tg_xyz_xxzzzzz_1, tg_xyz_xyyyyyy_0, tg_xyz_xyyyyyy_1, tg_xyz_xyyyyyz_0, \
                                         tg_xyz_xyyyyyz_1, tg_xyz_xyyyyzz_0, tg_xyz_xyyyyzz_1, tg_xyz_xyyyzzz_0, \
                                         tg_xyz_xyyyzzz_1, tg_xyz_xyyzzzz_0, tg_xyz_xyyzzzz_1, tg_xyz_xyzzzzz_0, \
                                         tg_xyz_xyzzzzz_1, tg_xyz_xzzzzzz_0, tg_xyz_xzzzzzz_1, tg_xyz_yyyyyyy_0, \
                                         tg_xyz_yyyyyyy_1, tg_xyz_yyyyyyz_0, tg_xyz_yyyyyyz_1, tg_xyz_yyyyyzz_0, \
                                         tg_xyz_yyyyyzz_1, tg_xyz_yyyyzzz_0, tg_xyz_yyyyzzz_1, tg_xyz_yyyzzzz_0, \
                                         tg_xyz_yyyzzzz_1, tg_xyz_yyzzzzz_0, tg_xyz_yyzzzzz_1, tg_xyz_yzzzzzz_0, \
                                         tg_xyz_yzzzzzz_1, tg_xyz_zzzzzzz_0, tg_xyz_zzzzzzz_1, tg_xzz_xxxxxxx_0, \
                                         tg_xzz_xxxxxxx_1, tg_xzz_xxxxxxy_0, tg_xzz_xxxxxxy_1, tg_xzz_xxxxxxz_0, \
                                         tg_xzz_xxxxxxz_1, tg_xzz_xxxxxyy_0, tg_xzz_xxxxxyy_1, tg_xzz_xxxxxyz_0, \
                                         tg_xzz_xxxxxyz_1, tg_xzz_xxxxxzz_0, tg_xzz_xxxxxzz_1, tg_xzz_xxxxyyy_0, \
                                         tg_xzz_xxxxyyy_1, tg_xzz_xxxxyyz_0, tg_xzz_xxxxyyz_1, tg_xzz_xxxxyzz_0, \
                                         tg_xzz_xxxxyzz_1, tg_xzz_xxxxzzz_0, tg_xzz_xxxxzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxxxz_xyyyyzz_0[j] = pb_x * tg_xxxz_xyyyyzz_0[j] + wp_x[j] * tg_xxxz_xyyyyzz_1[j] + 1.5 * fl1_fx * tg_xxz_xyyyyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xyyyyzz_1[j] + 0.5 * fl1_fxn * tg_xxxz_yyyyzz_1[j];

                    tg_xxxxz_xyyyzzz_0[j] = pb_x * tg_xxxz_xyyyzzz_0[j] + wp_x[j] * tg_xxxz_xyyyzzz_1[j] + 1.5 * fl1_fx * tg_xxz_xyyyzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xyyyzzz_1[j] + 0.5 * fl1_fxn * tg_xxxz_yyyzzz_1[j];

                    tg_xxxxz_xyyzzzz_0[j] = pb_x * tg_xxxz_xyyzzzz_0[j] + wp_x[j] * tg_xxxz_xyyzzzz_1[j] + 1.5 * fl1_fx * tg_xxz_xyyzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xyyzzzz_1[j] + 0.5 * fl1_fxn * tg_xxxz_yyzzzz_1[j];

                    tg_xxxxz_xyzzzzz_0[j] = pb_x * tg_xxxz_xyzzzzz_0[j] + wp_x[j] * tg_xxxz_xyzzzzz_1[j] + 1.5 * fl1_fx * tg_xxz_xyzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_xxxz_yzzzzz_1[j];

                    tg_xxxxz_xzzzzzz_0[j] = pb_x * tg_xxxz_xzzzzzz_0[j] + wp_x[j] * tg_xxxz_xzzzzzz_1[j] + 1.5 * fl1_fx * tg_xxz_xzzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xxxz_zzzzzz_1[j];

                    tg_xxxxz_yyyyyyy_0[j] = pb_x * tg_xxxz_yyyyyyy_0[j] + wp_x[j] * tg_xxxz_yyyyyyy_1[j] + 1.5 * fl1_fx * tg_xxz_yyyyyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_yyyyyyy_1[j];

                    tg_xxxxz_yyyyyyz_0[j] = pb_x * tg_xxxz_yyyyyyz_0[j] + wp_x[j] * tg_xxxz_yyyyyyz_1[j] + 1.5 * fl1_fx * tg_xxz_yyyyyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_yyyyyyz_1[j];

                    tg_xxxxz_yyyyyzz_0[j] = pb_x * tg_xxxz_yyyyyzz_0[j] + wp_x[j] * tg_xxxz_yyyyyzz_1[j] + 1.5 * fl1_fx * tg_xxz_yyyyyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_yyyyyzz_1[j];

                    tg_xxxxz_yyyyzzz_0[j] = pb_x * tg_xxxz_yyyyzzz_0[j] + wp_x[j] * tg_xxxz_yyyyzzz_1[j] + 1.5 * fl1_fx * tg_xxz_yyyyzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_yyyyzzz_1[j];

                    tg_xxxxz_yyyzzzz_0[j] = pb_x * tg_xxxz_yyyzzzz_0[j] + wp_x[j] * tg_xxxz_yyyzzzz_1[j] + 1.5 * fl1_fx * tg_xxz_yyyzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_yyyzzzz_1[j];

                    tg_xxxxz_yyzzzzz_0[j] = pb_x * tg_xxxz_yyzzzzz_0[j] + wp_x[j] * tg_xxxz_yyzzzzz_1[j] + 1.5 * fl1_fx * tg_xxz_yyzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_yyzzzzz_1[j];

                    tg_xxxxz_yzzzzzz_0[j] = pb_x * tg_xxxz_yzzzzzz_0[j] + wp_x[j] * tg_xxxz_yzzzzzz_1[j] + 1.5 * fl1_fx * tg_xxz_yzzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_yzzzzzz_1[j];

                    tg_xxxxz_zzzzzzz_0[j] = pb_x * tg_xxxz_zzzzzzz_0[j] + wp_x[j] * tg_xxxz_zzzzzzz_1[j] + 1.5 * fl1_fx * tg_xxz_zzzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_zzzzzzz_1[j];

                    tg_xxxyy_xxxxxxx_0[j] = pb_x * tg_xxyy_xxxxxxx_0[j] + wp_x[j] * tg_xxyy_xxxxxxx_1[j] + fl1_fx * tg_xyy_xxxxxxx_0[j] - fl1_fx * fl1_fza * tg_xyy_xxxxxxx_1[j] + 3.5 * fl1_fxn * tg_xxyy_xxxxxx_1[j];

                    tg_xxxyy_xxxxxxy_0[j] = pb_x * tg_xxyy_xxxxxxy_0[j] + wp_x[j] * tg_xxyy_xxxxxxy_1[j] + fl1_fx * tg_xyy_xxxxxxy_0[j] - fl1_fx * fl1_fza * tg_xyy_xxxxxxy_1[j] + 3.0 * fl1_fxn * tg_xxyy_xxxxxy_1[j];

                    tg_xxxyy_xxxxxxz_0[j] = pb_x * tg_xxyy_xxxxxxz_0[j] + wp_x[j] * tg_xxyy_xxxxxxz_1[j] + fl1_fx * tg_xyy_xxxxxxz_0[j] - fl1_fx * fl1_fza * tg_xyy_xxxxxxz_1[j] + 3.0 * fl1_fxn * tg_xxyy_xxxxxz_1[j];

                    tg_xxxyy_xxxxxyy_0[j] = pb_x * tg_xxyy_xxxxxyy_0[j] + wp_x[j] * tg_xxyy_xxxxxyy_1[j] + fl1_fx * tg_xyy_xxxxxyy_0[j] - fl1_fx * fl1_fza * tg_xyy_xxxxxyy_1[j] + 2.5 * fl1_fxn * tg_xxyy_xxxxyy_1[j];

                    tg_xxxyy_xxxxxyz_0[j] = pb_x * tg_xxyy_xxxxxyz_0[j] + wp_x[j] * tg_xxyy_xxxxxyz_1[j] + fl1_fx * tg_xyy_xxxxxyz_0[j] - fl1_fx * fl1_fza * tg_xyy_xxxxxyz_1[j] + 2.5 * fl1_fxn * tg_xxyy_xxxxyz_1[j];

                    tg_xxxyy_xxxxxzz_0[j] = pb_x * tg_xxyy_xxxxxzz_0[j] + wp_x[j] * tg_xxyy_xxxxxzz_1[j] + fl1_fx * tg_xyy_xxxxxzz_0[j] - fl1_fx * fl1_fza * tg_xyy_xxxxxzz_1[j] + 2.5 * fl1_fxn * tg_xxyy_xxxxzz_1[j];

                    tg_xxxyy_xxxxyyy_0[j] = pb_x * tg_xxyy_xxxxyyy_0[j] + wp_x[j] * tg_xxyy_xxxxyyy_1[j] + fl1_fx * tg_xyy_xxxxyyy_0[j] - fl1_fx * fl1_fza * tg_xyy_xxxxyyy_1[j] + 2.0 * fl1_fxn * tg_xxyy_xxxyyy_1[j];

                    tg_xxxyy_xxxxyyz_0[j] = pb_x * tg_xxyy_xxxxyyz_0[j] + wp_x[j] * tg_xxyy_xxxxyyz_1[j] + fl1_fx * tg_xyy_xxxxyyz_0[j] - fl1_fx * fl1_fza * tg_xyy_xxxxyyz_1[j] + 2.0 * fl1_fxn * tg_xxyy_xxxyyz_1[j];

                    tg_xxxyy_xxxxyzz_0[j] = pb_x * tg_xxyy_xxxxyzz_0[j] + wp_x[j] * tg_xxyy_xxxxyzz_1[j] + fl1_fx * tg_xyy_xxxxyzz_0[j] - fl1_fx * fl1_fza * tg_xyy_xxxxyzz_1[j] + 2.0 * fl1_fxn * tg_xxyy_xxxyzz_1[j];

                    tg_xxxyy_xxxxzzz_0[j] = pb_x * tg_xxyy_xxxxzzz_0[j] + wp_x[j] * tg_xxyy_xxxxzzz_1[j] + fl1_fx * tg_xyy_xxxxzzz_0[j] - fl1_fx * fl1_fza * tg_xyy_xxxxzzz_1[j] + 2.0 * fl1_fxn * tg_xxyy_xxxzzz_1[j];

                    tg_xxxyy_xxxyyyy_0[j] = pb_x * tg_xxyy_xxxyyyy_0[j] + wp_x[j] * tg_xxyy_xxxyyyy_1[j] + fl1_fx * tg_xyy_xxxyyyy_0[j] - fl1_fx * fl1_fza * tg_xyy_xxxyyyy_1[j] + 1.5 * fl1_fxn * tg_xxyy_xxyyyy_1[j];

                    tg_xxxyy_xxxyyyz_0[j] = pb_x * tg_xxyy_xxxyyyz_0[j] + wp_x[j] * tg_xxyy_xxxyyyz_1[j] + fl1_fx * tg_xyy_xxxyyyz_0[j] - fl1_fx * fl1_fza * tg_xyy_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_xxyy_xxyyyz_1[j];

                    tg_xxxyy_xxxyyzz_0[j] = pb_x * tg_xxyy_xxxyyzz_0[j] + wp_x[j] * tg_xxyy_xxxyyzz_1[j] + fl1_fx * tg_xyy_xxxyyzz_0[j] - fl1_fx * fl1_fza * tg_xyy_xxxyyzz_1[j] + 1.5 * fl1_fxn * tg_xxyy_xxyyzz_1[j];

                    tg_xxxyy_xxxyzzz_0[j] = pb_x * tg_xxyy_xxxyzzz_0[j] + wp_x[j] * tg_xxyy_xxxyzzz_1[j] + fl1_fx * tg_xyy_xxxyzzz_0[j] - fl1_fx * fl1_fza * tg_xyy_xxxyzzz_1[j] + 1.5 * fl1_fxn * tg_xxyy_xxyzzz_1[j];

                    tg_xxxyy_xxxzzzz_0[j] = pb_x * tg_xxyy_xxxzzzz_0[j] + wp_x[j] * tg_xxyy_xxxzzzz_1[j] + fl1_fx * tg_xyy_xxxzzzz_0[j] - fl1_fx * fl1_fza * tg_xyy_xxxzzzz_1[j] + 1.5 * fl1_fxn * tg_xxyy_xxzzzz_1[j];

                    tg_xxxyy_xxyyyyy_0[j] = pb_x * tg_xxyy_xxyyyyy_0[j] + wp_x[j] * tg_xxyy_xxyyyyy_1[j] + fl1_fx * tg_xyy_xxyyyyy_0[j] - fl1_fx * fl1_fza * tg_xyy_xxyyyyy_1[j] + fl1_fxn * tg_xxyy_xyyyyy_1[j];

                    tg_xxxyy_xxyyyyz_0[j] = pb_x * tg_xxyy_xxyyyyz_0[j] + wp_x[j] * tg_xxyy_xxyyyyz_1[j] + fl1_fx * tg_xyy_xxyyyyz_0[j] - fl1_fx * fl1_fza * tg_xyy_xxyyyyz_1[j] + fl1_fxn * tg_xxyy_xyyyyz_1[j];

                    tg_xxxyy_xxyyyzz_0[j] = pb_x * tg_xxyy_xxyyyzz_0[j] + wp_x[j] * tg_xxyy_xxyyyzz_1[j] + fl1_fx * tg_xyy_xxyyyzz_0[j] - fl1_fx * fl1_fza * tg_xyy_xxyyyzz_1[j] + fl1_fxn * tg_xxyy_xyyyzz_1[j];

                    tg_xxxyy_xxyyzzz_0[j] = pb_x * tg_xxyy_xxyyzzz_0[j] + wp_x[j] * tg_xxyy_xxyyzzz_1[j] + fl1_fx * tg_xyy_xxyyzzz_0[j] - fl1_fx * fl1_fza * tg_xyy_xxyyzzz_1[j] + fl1_fxn * tg_xxyy_xyyzzz_1[j];

                    tg_xxxyy_xxyzzzz_0[j] = pb_x * tg_xxyy_xxyzzzz_0[j] + wp_x[j] * tg_xxyy_xxyzzzz_1[j] + fl1_fx * tg_xyy_xxyzzzz_0[j] - fl1_fx * fl1_fza * tg_xyy_xxyzzzz_1[j] + fl1_fxn * tg_xxyy_xyzzzz_1[j];

                    tg_xxxyy_xxzzzzz_0[j] = pb_x * tg_xxyy_xxzzzzz_0[j] + wp_x[j] * tg_xxyy_xxzzzzz_1[j] + fl1_fx * tg_xyy_xxzzzzz_0[j] - fl1_fx * fl1_fza * tg_xyy_xxzzzzz_1[j] + fl1_fxn * tg_xxyy_xzzzzz_1[j];

                    tg_xxxyy_xyyyyyy_0[j] = pb_x * tg_xxyy_xyyyyyy_0[j] + wp_x[j] * tg_xxyy_xyyyyyy_1[j] + fl1_fx * tg_xyy_xyyyyyy_0[j] - fl1_fx * fl1_fza * tg_xyy_xyyyyyy_1[j] + 0.5 * fl1_fxn * tg_xxyy_yyyyyy_1[j];

                    tg_xxxyy_xyyyyyz_0[j] = pb_x * tg_xxyy_xyyyyyz_0[j] + wp_x[j] * tg_xxyy_xyyyyyz_1[j] + fl1_fx * tg_xyy_xyyyyyz_0[j] - fl1_fx * fl1_fza * tg_xyy_xyyyyyz_1[j] + 0.5 * fl1_fxn * tg_xxyy_yyyyyz_1[j];

                    tg_xxxyy_xyyyyzz_0[j] = pb_x * tg_xxyy_xyyyyzz_0[j] + wp_x[j] * tg_xxyy_xyyyyzz_1[j] + fl1_fx * tg_xyy_xyyyyzz_0[j] - fl1_fx * fl1_fza * tg_xyy_xyyyyzz_1[j] + 0.5 * fl1_fxn * tg_xxyy_yyyyzz_1[j];

                    tg_xxxyy_xyyyzzz_0[j] = pb_x * tg_xxyy_xyyyzzz_0[j] + wp_x[j] * tg_xxyy_xyyyzzz_1[j] + fl1_fx * tg_xyy_xyyyzzz_0[j] - fl1_fx * fl1_fza * tg_xyy_xyyyzzz_1[j] + 0.5 * fl1_fxn * tg_xxyy_yyyzzz_1[j];

                    tg_xxxyy_xyyzzzz_0[j] = pb_x * tg_xxyy_xyyzzzz_0[j] + wp_x[j] * tg_xxyy_xyyzzzz_1[j] + fl1_fx * tg_xyy_xyyzzzz_0[j] - fl1_fx * fl1_fza * tg_xyy_xyyzzzz_1[j] + 0.5 * fl1_fxn * tg_xxyy_yyzzzz_1[j];

                    tg_xxxyy_xyzzzzz_0[j] = pb_x * tg_xxyy_xyzzzzz_0[j] + wp_x[j] * tg_xxyy_xyzzzzz_1[j] + fl1_fx * tg_xyy_xyzzzzz_0[j] - fl1_fx * fl1_fza * tg_xyy_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_xxyy_yzzzzz_1[j];

                    tg_xxxyy_xzzzzzz_0[j] = pb_x * tg_xxyy_xzzzzzz_0[j] + wp_x[j] * tg_xxyy_xzzzzzz_1[j] + fl1_fx * tg_xyy_xzzzzzz_0[j] - fl1_fx * fl1_fza * tg_xyy_xzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xxyy_zzzzzz_1[j];

                    tg_xxxyy_yyyyyyy_0[j] = pb_x * tg_xxyy_yyyyyyy_0[j] + wp_x[j] * tg_xxyy_yyyyyyy_1[j] + fl1_fx * tg_xyy_yyyyyyy_0[j] - fl1_fx * fl1_fza * tg_xyy_yyyyyyy_1[j];

                    tg_xxxyy_yyyyyyz_0[j] = pb_x * tg_xxyy_yyyyyyz_0[j] + wp_x[j] * tg_xxyy_yyyyyyz_1[j] + fl1_fx * tg_xyy_yyyyyyz_0[j] - fl1_fx * fl1_fza * tg_xyy_yyyyyyz_1[j];

                    tg_xxxyy_yyyyyzz_0[j] = pb_x * tg_xxyy_yyyyyzz_0[j] + wp_x[j] * tg_xxyy_yyyyyzz_1[j] + fl1_fx * tg_xyy_yyyyyzz_0[j] - fl1_fx * fl1_fza * tg_xyy_yyyyyzz_1[j];

                    tg_xxxyy_yyyyzzz_0[j] = pb_x * tg_xxyy_yyyyzzz_0[j] + wp_x[j] * tg_xxyy_yyyyzzz_1[j] + fl1_fx * tg_xyy_yyyyzzz_0[j] - fl1_fx * fl1_fza * tg_xyy_yyyyzzz_1[j];

                    tg_xxxyy_yyyzzzz_0[j] = pb_x * tg_xxyy_yyyzzzz_0[j] + wp_x[j] * tg_xxyy_yyyzzzz_1[j] + fl1_fx * tg_xyy_yyyzzzz_0[j] - fl1_fx * fl1_fza * tg_xyy_yyyzzzz_1[j];

                    tg_xxxyy_yyzzzzz_0[j] = pb_x * tg_xxyy_yyzzzzz_0[j] + wp_x[j] * tg_xxyy_yyzzzzz_1[j] + fl1_fx * tg_xyy_yyzzzzz_0[j] - fl1_fx * fl1_fza * tg_xyy_yyzzzzz_1[j];

                    tg_xxxyy_yzzzzzz_0[j] = pb_x * tg_xxyy_yzzzzzz_0[j] + wp_x[j] * tg_xxyy_yzzzzzz_1[j] + fl1_fx * tg_xyy_yzzzzzz_0[j] - fl1_fx * fl1_fza * tg_xyy_yzzzzzz_1[j];

                    tg_xxxyy_zzzzzzz_0[j] = pb_x * tg_xxyy_zzzzzzz_0[j] + wp_x[j] * tg_xxyy_zzzzzzz_1[j] + fl1_fx * tg_xyy_zzzzzzz_0[j] - fl1_fx * fl1_fza * tg_xyy_zzzzzzz_1[j];

                    tg_xxxyz_xxxxxxx_0[j] = pb_x * tg_xxyz_xxxxxxx_0[j] + wp_x[j] * tg_xxyz_xxxxxxx_1[j] + fl1_fx * tg_xyz_xxxxxxx_0[j] - fl1_fx * fl1_fza * tg_xyz_xxxxxxx_1[j] + 3.5 * fl1_fxn * tg_xxyz_xxxxxx_1[j];

                    tg_xxxyz_xxxxxxy_0[j] = pb_x * tg_xxyz_xxxxxxy_0[j] + wp_x[j] * tg_xxyz_xxxxxxy_1[j] + fl1_fx * tg_xyz_xxxxxxy_0[j] - fl1_fx * fl1_fza * tg_xyz_xxxxxxy_1[j] + 3.0 * fl1_fxn * tg_xxyz_xxxxxy_1[j];

                    tg_xxxyz_xxxxxxz_0[j] = pb_x * tg_xxyz_xxxxxxz_0[j] + wp_x[j] * tg_xxyz_xxxxxxz_1[j] + fl1_fx * tg_xyz_xxxxxxz_0[j] - fl1_fx * fl1_fza * tg_xyz_xxxxxxz_1[j] + 3.0 * fl1_fxn * tg_xxyz_xxxxxz_1[j];

                    tg_xxxyz_xxxxxyy_0[j] = pb_x * tg_xxyz_xxxxxyy_0[j] + wp_x[j] * tg_xxyz_xxxxxyy_1[j] + fl1_fx * tg_xyz_xxxxxyy_0[j] - fl1_fx * fl1_fza * tg_xyz_xxxxxyy_1[j] + 2.5 * fl1_fxn * tg_xxyz_xxxxyy_1[j];

                    tg_xxxyz_xxxxxyz_0[j] = pb_x * tg_xxyz_xxxxxyz_0[j] + wp_x[j] * tg_xxyz_xxxxxyz_1[j] + fl1_fx * tg_xyz_xxxxxyz_0[j] - fl1_fx * fl1_fza * tg_xyz_xxxxxyz_1[j] + 2.5 * fl1_fxn * tg_xxyz_xxxxyz_1[j];

                    tg_xxxyz_xxxxxzz_0[j] = pb_x * tg_xxyz_xxxxxzz_0[j] + wp_x[j] * tg_xxyz_xxxxxzz_1[j] + fl1_fx * tg_xyz_xxxxxzz_0[j] - fl1_fx * fl1_fza * tg_xyz_xxxxxzz_1[j] + 2.5 * fl1_fxn * tg_xxyz_xxxxzz_1[j];

                    tg_xxxyz_xxxxyyy_0[j] = pb_x * tg_xxyz_xxxxyyy_0[j] + wp_x[j] * tg_xxyz_xxxxyyy_1[j] + fl1_fx * tg_xyz_xxxxyyy_0[j] - fl1_fx * fl1_fza * tg_xyz_xxxxyyy_1[j] + 2.0 * fl1_fxn * tg_xxyz_xxxyyy_1[j];

                    tg_xxxyz_xxxxyyz_0[j] = pb_x * tg_xxyz_xxxxyyz_0[j] + wp_x[j] * tg_xxyz_xxxxyyz_1[j] + fl1_fx * tg_xyz_xxxxyyz_0[j] - fl1_fx * fl1_fza * tg_xyz_xxxxyyz_1[j] + 2.0 * fl1_fxn * tg_xxyz_xxxyyz_1[j];

                    tg_xxxyz_xxxxyzz_0[j] = pb_x * tg_xxyz_xxxxyzz_0[j] + wp_x[j] * tg_xxyz_xxxxyzz_1[j] + fl1_fx * tg_xyz_xxxxyzz_0[j] - fl1_fx * fl1_fza * tg_xyz_xxxxyzz_1[j] + 2.0 * fl1_fxn * tg_xxyz_xxxyzz_1[j];

                    tg_xxxyz_xxxxzzz_0[j] = pb_x * tg_xxyz_xxxxzzz_0[j] + wp_x[j] * tg_xxyz_xxxxzzz_1[j] + fl1_fx * tg_xyz_xxxxzzz_0[j] - fl1_fx * fl1_fza * tg_xyz_xxxxzzz_1[j] + 2.0 * fl1_fxn * tg_xxyz_xxxzzz_1[j];

                    tg_xxxyz_xxxyyyy_0[j] = pb_x * tg_xxyz_xxxyyyy_0[j] + wp_x[j] * tg_xxyz_xxxyyyy_1[j] + fl1_fx * tg_xyz_xxxyyyy_0[j] - fl1_fx * fl1_fza * tg_xyz_xxxyyyy_1[j] + 1.5 * fl1_fxn * tg_xxyz_xxyyyy_1[j];

                    tg_xxxyz_xxxyyyz_0[j] = pb_x * tg_xxyz_xxxyyyz_0[j] + wp_x[j] * tg_xxyz_xxxyyyz_1[j] + fl1_fx * tg_xyz_xxxyyyz_0[j] - fl1_fx * fl1_fza * tg_xyz_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_xxyz_xxyyyz_1[j];

                    tg_xxxyz_xxxyyzz_0[j] = pb_x * tg_xxyz_xxxyyzz_0[j] + wp_x[j] * tg_xxyz_xxxyyzz_1[j] + fl1_fx * tg_xyz_xxxyyzz_0[j] - fl1_fx * fl1_fza * tg_xyz_xxxyyzz_1[j] + 1.5 * fl1_fxn * tg_xxyz_xxyyzz_1[j];

                    tg_xxxyz_xxxyzzz_0[j] = pb_x * tg_xxyz_xxxyzzz_0[j] + wp_x[j] * tg_xxyz_xxxyzzz_1[j] + fl1_fx * tg_xyz_xxxyzzz_0[j] - fl1_fx * fl1_fza * tg_xyz_xxxyzzz_1[j] + 1.5 * fl1_fxn * tg_xxyz_xxyzzz_1[j];

                    tg_xxxyz_xxxzzzz_0[j] = pb_x * tg_xxyz_xxxzzzz_0[j] + wp_x[j] * tg_xxyz_xxxzzzz_1[j] + fl1_fx * tg_xyz_xxxzzzz_0[j] - fl1_fx * fl1_fza * tg_xyz_xxxzzzz_1[j] + 1.5 * fl1_fxn * tg_xxyz_xxzzzz_1[j];

                    tg_xxxyz_xxyyyyy_0[j] = pb_x * tg_xxyz_xxyyyyy_0[j] + wp_x[j] * tg_xxyz_xxyyyyy_1[j] + fl1_fx * tg_xyz_xxyyyyy_0[j] - fl1_fx * fl1_fza * tg_xyz_xxyyyyy_1[j] + fl1_fxn * tg_xxyz_xyyyyy_1[j];

                    tg_xxxyz_xxyyyyz_0[j] = pb_x * tg_xxyz_xxyyyyz_0[j] + wp_x[j] * tg_xxyz_xxyyyyz_1[j] + fl1_fx * tg_xyz_xxyyyyz_0[j] - fl1_fx * fl1_fza * tg_xyz_xxyyyyz_1[j] + fl1_fxn * tg_xxyz_xyyyyz_1[j];

                    tg_xxxyz_xxyyyzz_0[j] = pb_x * tg_xxyz_xxyyyzz_0[j] + wp_x[j] * tg_xxyz_xxyyyzz_1[j] + fl1_fx * tg_xyz_xxyyyzz_0[j] - fl1_fx * fl1_fza * tg_xyz_xxyyyzz_1[j] + fl1_fxn * tg_xxyz_xyyyzz_1[j];

                    tg_xxxyz_xxyyzzz_0[j] = pb_x * tg_xxyz_xxyyzzz_0[j] + wp_x[j] * tg_xxyz_xxyyzzz_1[j] + fl1_fx * tg_xyz_xxyyzzz_0[j] - fl1_fx * fl1_fza * tg_xyz_xxyyzzz_1[j] + fl1_fxn * tg_xxyz_xyyzzz_1[j];

                    tg_xxxyz_xxyzzzz_0[j] = pb_x * tg_xxyz_xxyzzzz_0[j] + wp_x[j] * tg_xxyz_xxyzzzz_1[j] + fl1_fx * tg_xyz_xxyzzzz_0[j] - fl1_fx * fl1_fza * tg_xyz_xxyzzzz_1[j] + fl1_fxn * tg_xxyz_xyzzzz_1[j];

                    tg_xxxyz_xxzzzzz_0[j] = pb_x * tg_xxyz_xxzzzzz_0[j] + wp_x[j] * tg_xxyz_xxzzzzz_1[j] + fl1_fx * tg_xyz_xxzzzzz_0[j] - fl1_fx * fl1_fza * tg_xyz_xxzzzzz_1[j] + fl1_fxn * tg_xxyz_xzzzzz_1[j];

                    tg_xxxyz_xyyyyyy_0[j] = pb_x * tg_xxyz_xyyyyyy_0[j] + wp_x[j] * tg_xxyz_xyyyyyy_1[j] + fl1_fx * tg_xyz_xyyyyyy_0[j] - fl1_fx * fl1_fza * tg_xyz_xyyyyyy_1[j] + 0.5 * fl1_fxn * tg_xxyz_yyyyyy_1[j];

                    tg_xxxyz_xyyyyyz_0[j] = pb_x * tg_xxyz_xyyyyyz_0[j] + wp_x[j] * tg_xxyz_xyyyyyz_1[j] + fl1_fx * tg_xyz_xyyyyyz_0[j] - fl1_fx * fl1_fza * tg_xyz_xyyyyyz_1[j] + 0.5 * fl1_fxn * tg_xxyz_yyyyyz_1[j];

                    tg_xxxyz_xyyyyzz_0[j] = pb_x * tg_xxyz_xyyyyzz_0[j] + wp_x[j] * tg_xxyz_xyyyyzz_1[j] + fl1_fx * tg_xyz_xyyyyzz_0[j] - fl1_fx * fl1_fza * tg_xyz_xyyyyzz_1[j] + 0.5 * fl1_fxn * tg_xxyz_yyyyzz_1[j];

                    tg_xxxyz_xyyyzzz_0[j] = pb_x * tg_xxyz_xyyyzzz_0[j] + wp_x[j] * tg_xxyz_xyyyzzz_1[j] + fl1_fx * tg_xyz_xyyyzzz_0[j] - fl1_fx * fl1_fza * tg_xyz_xyyyzzz_1[j] + 0.5 * fl1_fxn * tg_xxyz_yyyzzz_1[j];

                    tg_xxxyz_xyyzzzz_0[j] = pb_x * tg_xxyz_xyyzzzz_0[j] + wp_x[j] * tg_xxyz_xyyzzzz_1[j] + fl1_fx * tg_xyz_xyyzzzz_0[j] - fl1_fx * fl1_fza * tg_xyz_xyyzzzz_1[j] + 0.5 * fl1_fxn * tg_xxyz_yyzzzz_1[j];

                    tg_xxxyz_xyzzzzz_0[j] = pb_x * tg_xxyz_xyzzzzz_0[j] + wp_x[j] * tg_xxyz_xyzzzzz_1[j] + fl1_fx * tg_xyz_xyzzzzz_0[j] - fl1_fx * fl1_fza * tg_xyz_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_xxyz_yzzzzz_1[j];

                    tg_xxxyz_xzzzzzz_0[j] = pb_x * tg_xxyz_xzzzzzz_0[j] + wp_x[j] * tg_xxyz_xzzzzzz_1[j] + fl1_fx * tg_xyz_xzzzzzz_0[j] - fl1_fx * fl1_fza * tg_xyz_xzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xxyz_zzzzzz_1[j];

                    tg_xxxyz_yyyyyyy_0[j] = pb_x * tg_xxyz_yyyyyyy_0[j] + wp_x[j] * tg_xxyz_yyyyyyy_1[j] + fl1_fx * tg_xyz_yyyyyyy_0[j] - fl1_fx * fl1_fza * tg_xyz_yyyyyyy_1[j];

                    tg_xxxyz_yyyyyyz_0[j] = pb_x * tg_xxyz_yyyyyyz_0[j] + wp_x[j] * tg_xxyz_yyyyyyz_1[j] + fl1_fx * tg_xyz_yyyyyyz_0[j] - fl1_fx * fl1_fza * tg_xyz_yyyyyyz_1[j];

                    tg_xxxyz_yyyyyzz_0[j] = pb_x * tg_xxyz_yyyyyzz_0[j] + wp_x[j] * tg_xxyz_yyyyyzz_1[j] + fl1_fx * tg_xyz_yyyyyzz_0[j] - fl1_fx * fl1_fza * tg_xyz_yyyyyzz_1[j];

                    tg_xxxyz_yyyyzzz_0[j] = pb_x * tg_xxyz_yyyyzzz_0[j] + wp_x[j] * tg_xxyz_yyyyzzz_1[j] + fl1_fx * tg_xyz_yyyyzzz_0[j] - fl1_fx * fl1_fza * tg_xyz_yyyyzzz_1[j];

                    tg_xxxyz_yyyzzzz_0[j] = pb_x * tg_xxyz_yyyzzzz_0[j] + wp_x[j] * tg_xxyz_yyyzzzz_1[j] + fl1_fx * tg_xyz_yyyzzzz_0[j] - fl1_fx * fl1_fza * tg_xyz_yyyzzzz_1[j];

                    tg_xxxyz_yyzzzzz_0[j] = pb_x * tg_xxyz_yyzzzzz_0[j] + wp_x[j] * tg_xxyz_yyzzzzz_1[j] + fl1_fx * tg_xyz_yyzzzzz_0[j] - fl1_fx * fl1_fza * tg_xyz_yyzzzzz_1[j];

                    tg_xxxyz_yzzzzzz_0[j] = pb_x * tg_xxyz_yzzzzzz_0[j] + wp_x[j] * tg_xxyz_yzzzzzz_1[j] + fl1_fx * tg_xyz_yzzzzzz_0[j] - fl1_fx * fl1_fza * tg_xyz_yzzzzzz_1[j];

                    tg_xxxyz_zzzzzzz_0[j] = pb_x * tg_xxyz_zzzzzzz_0[j] + wp_x[j] * tg_xxyz_zzzzzzz_1[j] + fl1_fx * tg_xyz_zzzzzzz_0[j] - fl1_fx * fl1_fza * tg_xyz_zzzzzzz_1[j];

                    tg_xxxzz_xxxxxxx_0[j] = pb_x * tg_xxzz_xxxxxxx_0[j] + wp_x[j] * tg_xxzz_xxxxxxx_1[j] + fl1_fx * tg_xzz_xxxxxxx_0[j] - fl1_fx * fl1_fza * tg_xzz_xxxxxxx_1[j] + 3.5 * fl1_fxn * tg_xxzz_xxxxxx_1[j];

                    tg_xxxzz_xxxxxxy_0[j] = pb_x * tg_xxzz_xxxxxxy_0[j] + wp_x[j] * tg_xxzz_xxxxxxy_1[j] + fl1_fx * tg_xzz_xxxxxxy_0[j] - fl1_fx * fl1_fza * tg_xzz_xxxxxxy_1[j] + 3.0 * fl1_fxn * tg_xxzz_xxxxxy_1[j];

                    tg_xxxzz_xxxxxxz_0[j] = pb_x * tg_xxzz_xxxxxxz_0[j] + wp_x[j] * tg_xxzz_xxxxxxz_1[j] + fl1_fx * tg_xzz_xxxxxxz_0[j] - fl1_fx * fl1_fza * tg_xzz_xxxxxxz_1[j] + 3.0 * fl1_fxn * tg_xxzz_xxxxxz_1[j];

                    tg_xxxzz_xxxxxyy_0[j] = pb_x * tg_xxzz_xxxxxyy_0[j] + wp_x[j] * tg_xxzz_xxxxxyy_1[j] + fl1_fx * tg_xzz_xxxxxyy_0[j] - fl1_fx * fl1_fza * tg_xzz_xxxxxyy_1[j] + 2.5 * fl1_fxn * tg_xxzz_xxxxyy_1[j];

                    tg_xxxzz_xxxxxyz_0[j] = pb_x * tg_xxzz_xxxxxyz_0[j] + wp_x[j] * tg_xxzz_xxxxxyz_1[j] + fl1_fx * tg_xzz_xxxxxyz_0[j] - fl1_fx * fl1_fza * tg_xzz_xxxxxyz_1[j] + 2.5 * fl1_fxn * tg_xxzz_xxxxyz_1[j];

                    tg_xxxzz_xxxxxzz_0[j] = pb_x * tg_xxzz_xxxxxzz_0[j] + wp_x[j] * tg_xxzz_xxxxxzz_1[j] + fl1_fx * tg_xzz_xxxxxzz_0[j] - fl1_fx * fl1_fza * tg_xzz_xxxxxzz_1[j] + 2.5 * fl1_fxn * tg_xxzz_xxxxzz_1[j];

                    tg_xxxzz_xxxxyyy_0[j] = pb_x * tg_xxzz_xxxxyyy_0[j] + wp_x[j] * tg_xxzz_xxxxyyy_1[j] + fl1_fx * tg_xzz_xxxxyyy_0[j] - fl1_fx * fl1_fza * tg_xzz_xxxxyyy_1[j] + 2.0 * fl1_fxn * tg_xxzz_xxxyyy_1[j];

                    tg_xxxzz_xxxxyyz_0[j] = pb_x * tg_xxzz_xxxxyyz_0[j] + wp_x[j] * tg_xxzz_xxxxyyz_1[j] + fl1_fx * tg_xzz_xxxxyyz_0[j] - fl1_fx * fl1_fza * tg_xzz_xxxxyyz_1[j] + 2.0 * fl1_fxn * tg_xxzz_xxxyyz_1[j];

                    tg_xxxzz_xxxxyzz_0[j] = pb_x * tg_xxzz_xxxxyzz_0[j] + wp_x[j] * tg_xxzz_xxxxyzz_1[j] + fl1_fx * tg_xzz_xxxxyzz_0[j] - fl1_fx * fl1_fza * tg_xzz_xxxxyzz_1[j] + 2.0 * fl1_fxn * tg_xxzz_xxxyzz_1[j];

                    tg_xxxzz_xxxxzzz_0[j] = pb_x * tg_xxzz_xxxxzzz_0[j] + wp_x[j] * tg_xxzz_xxxxzzz_1[j] + fl1_fx * tg_xzz_xxxxzzz_0[j] - fl1_fx * fl1_fza * tg_xzz_xxxxzzz_1[j] + 2.0 * fl1_fxn * tg_xxzz_xxxzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSK_190_285(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (190,285)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {5, -1, -1, -1},
                                             {7, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_7_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_4_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_3_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_4_6_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {6, -1, -1, -1}, 
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

                auto tg_xxzz_xxxyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 190); 

                auto tg_xxzz_xxxyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 191); 

                auto tg_xxzz_xxxyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 192); 

                auto tg_xxzz_xxxyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 193); 

                auto tg_xxzz_xxxzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 194); 

                auto tg_xxzz_xxyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 195); 

                auto tg_xxzz_xxyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 196); 

                auto tg_xxzz_xxyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 197); 

                auto tg_xxzz_xxyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 198); 

                auto tg_xxzz_xxyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 199); 

                auto tg_xxzz_xxzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 200); 

                auto tg_xxzz_xyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 201); 

                auto tg_xxzz_xyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 202); 

                auto tg_xxzz_xyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 203); 

                auto tg_xxzz_xyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 204); 

                auto tg_xxzz_xyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 205); 

                auto tg_xxzz_xyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 206); 

                auto tg_xxzz_xzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 207); 

                auto tg_xxzz_yyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 208); 

                auto tg_xxzz_yyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 209); 

                auto tg_xxzz_yyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 210); 

                auto tg_xxzz_yyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 211); 

                auto tg_xxzz_yyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 212); 

                auto tg_xxzz_yyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 213); 

                auto tg_xxzz_yzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 214); 

                auto tg_xxzz_zzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 215); 

                auto tg_xyyy_xxxxxxx_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 216); 

                auto tg_xyyy_xxxxxxy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 217); 

                auto tg_xyyy_xxxxxxz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 218); 

                auto tg_xyyy_xxxxxyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 219); 

                auto tg_xyyy_xxxxxyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 220); 

                auto tg_xyyy_xxxxxzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 221); 

                auto tg_xyyy_xxxxyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 222); 

                auto tg_xyyy_xxxxyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 223); 

                auto tg_xyyy_xxxxyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 224); 

                auto tg_xyyy_xxxxzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 225); 

                auto tg_xyyy_xxxyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 226); 

                auto tg_xyyy_xxxyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 227); 

                auto tg_xyyy_xxxyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 228); 

                auto tg_xyyy_xxxyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 229); 

                auto tg_xyyy_xxxzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 230); 

                auto tg_xyyy_xxyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 231); 

                auto tg_xyyy_xxyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 232); 

                auto tg_xyyy_xxyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 233); 

                auto tg_xyyy_xxyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 234); 

                auto tg_xyyy_xxyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 235); 

                auto tg_xyyy_xxzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 236); 

                auto tg_xyyy_xyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 237); 

                auto tg_xyyy_xyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 238); 

                auto tg_xyyy_xyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 239); 

                auto tg_xyyy_xyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 240); 

                auto tg_xyyy_xyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 241); 

                auto tg_xyyy_xyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 242); 

                auto tg_xyyy_xzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 243); 

                auto tg_xyyy_yyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 244); 

                auto tg_xyyy_yyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 245); 

                auto tg_xyyy_yyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 246); 

                auto tg_xyyy_yyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 247); 

                auto tg_xyyy_yyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 248); 

                auto tg_xyyy_yyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 249); 

                auto tg_xyyy_yzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 250); 

                auto tg_xyyy_zzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 251); 

                auto tg_xyyz_xxxxxxx_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 252); 

                auto tg_xyyz_xxxxxxy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 253); 

                auto tg_xyyz_xxxxxxz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 254); 

                auto tg_xyyz_xxxxxyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 255); 

                auto tg_xyyz_xxxxxyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 256); 

                auto tg_xyyz_xxxxxzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 257); 

                auto tg_xyyz_xxxxyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 258); 

                auto tg_xyyz_xxxxyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 259); 

                auto tg_xyyz_xxxxyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 260); 

                auto tg_xyyz_xxxxzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 261); 

                auto tg_xyyz_xxxyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 262); 

                auto tg_xyyz_xxxyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 263); 

                auto tg_xyyz_xxxyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 264); 

                auto tg_xyyz_xxxyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 265); 

                auto tg_xyyz_xxxzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 266); 

                auto tg_xyyz_xxyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 267); 

                auto tg_xyyz_xxyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 268); 

                auto tg_xyyz_xxyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 269); 

                auto tg_xyyz_xxyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 270); 

                auto tg_xyyz_xxyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 271); 

                auto tg_xyyz_xxzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 272); 

                auto tg_xyyz_xyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 273); 

                auto tg_xyyz_xyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 274); 

                auto tg_xyyz_xyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 275); 

                auto tg_xyyz_xyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 276); 

                auto tg_xyyz_xyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 277); 

                auto tg_xyyz_xyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 278); 

                auto tg_xyyz_xzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 279); 

                auto tg_xyyz_yyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 280); 

                auto tg_xyyz_yyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 281); 

                auto tg_xyyz_yyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 282); 

                auto tg_xyyz_yyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 283); 

                auto tg_xyyz_yyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 284); 

                auto tg_xxzz_xxxyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 190); 

                auto tg_xxzz_xxxyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 191); 

                auto tg_xxzz_xxxyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 192); 

                auto tg_xxzz_xxxyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 193); 

                auto tg_xxzz_xxxzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 194); 

                auto tg_xxzz_xxyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 195); 

                auto tg_xxzz_xxyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 196); 

                auto tg_xxzz_xxyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 197); 

                auto tg_xxzz_xxyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 198); 

                auto tg_xxzz_xxyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 199); 

                auto tg_xxzz_xxzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 200); 

                auto tg_xxzz_xyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 201); 

                auto tg_xxzz_xyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 202); 

                auto tg_xxzz_xyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 203); 

                auto tg_xxzz_xyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 204); 

                auto tg_xxzz_xyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 205); 

                auto tg_xxzz_xyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 206); 

                auto tg_xxzz_xzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 207); 

                auto tg_xxzz_yyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 208); 

                auto tg_xxzz_yyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 209); 

                auto tg_xxzz_yyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 210); 

                auto tg_xxzz_yyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 211); 

                auto tg_xxzz_yyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 212); 

                auto tg_xxzz_yyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 213); 

                auto tg_xxzz_yzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 214); 

                auto tg_xxzz_zzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 215); 

                auto tg_xyyy_xxxxxxx_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 216); 

                auto tg_xyyy_xxxxxxy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 217); 

                auto tg_xyyy_xxxxxxz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 218); 

                auto tg_xyyy_xxxxxyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 219); 

                auto tg_xyyy_xxxxxyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 220); 

                auto tg_xyyy_xxxxxzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 221); 

                auto tg_xyyy_xxxxyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 222); 

                auto tg_xyyy_xxxxyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 223); 

                auto tg_xyyy_xxxxyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 224); 

                auto tg_xyyy_xxxxzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 225); 

                auto tg_xyyy_xxxyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 226); 

                auto tg_xyyy_xxxyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 227); 

                auto tg_xyyy_xxxyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 228); 

                auto tg_xyyy_xxxyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 229); 

                auto tg_xyyy_xxxzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 230); 

                auto tg_xyyy_xxyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 231); 

                auto tg_xyyy_xxyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 232); 

                auto tg_xyyy_xxyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 233); 

                auto tg_xyyy_xxyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 234); 

                auto tg_xyyy_xxyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 235); 

                auto tg_xyyy_xxzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 236); 

                auto tg_xyyy_xyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 237); 

                auto tg_xyyy_xyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 238); 

                auto tg_xyyy_xyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 239); 

                auto tg_xyyy_xyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 240); 

                auto tg_xyyy_xyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 241); 

                auto tg_xyyy_xyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 242); 

                auto tg_xyyy_xzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 243); 

                auto tg_xyyy_yyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 244); 

                auto tg_xyyy_yyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 245); 

                auto tg_xyyy_yyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 246); 

                auto tg_xyyy_yyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 247); 

                auto tg_xyyy_yyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 248); 

                auto tg_xyyy_yyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 249); 

                auto tg_xyyy_yzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 250); 

                auto tg_xyyy_zzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 251); 

                auto tg_xyyz_xxxxxxx_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 252); 

                auto tg_xyyz_xxxxxxy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 253); 

                auto tg_xyyz_xxxxxxz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 254); 

                auto tg_xyyz_xxxxxyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 255); 

                auto tg_xyyz_xxxxxyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 256); 

                auto tg_xyyz_xxxxxzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 257); 

                auto tg_xyyz_xxxxyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 258); 

                auto tg_xyyz_xxxxyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 259); 

                auto tg_xyyz_xxxxyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 260); 

                auto tg_xyyz_xxxxzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 261); 

                auto tg_xyyz_xxxyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 262); 

                auto tg_xyyz_xxxyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 263); 

                auto tg_xyyz_xxxyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 264); 

                auto tg_xyyz_xxxyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 265); 

                auto tg_xyyz_xxxzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 266); 

                auto tg_xyyz_xxyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 267); 

                auto tg_xyyz_xxyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 268); 

                auto tg_xyyz_xxyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 269); 

                auto tg_xyyz_xxyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 270); 

                auto tg_xyyz_xxyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 271); 

                auto tg_xyyz_xxzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 272); 

                auto tg_xyyz_xyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 273); 

                auto tg_xyyz_xyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 274); 

                auto tg_xyyz_xyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 275); 

                auto tg_xyyz_xyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 276); 

                auto tg_xyyz_xyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 277); 

                auto tg_xyyz_xyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 278); 

                auto tg_xyyz_xzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 279); 

                auto tg_xyyz_yyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 280); 

                auto tg_xyyz_yyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 281); 

                auto tg_xyyz_yyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 282); 

                auto tg_xyyz_yyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 283); 

                auto tg_xyyz_yyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 284); 

                auto tg_xzz_xxxyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 190); 

                auto tg_xzz_xxxyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 191); 

                auto tg_xzz_xxxyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 192); 

                auto tg_xzz_xxxyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 193); 

                auto tg_xzz_xxxzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 194); 

                auto tg_xzz_xxyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 195); 

                auto tg_xzz_xxyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 196); 

                auto tg_xzz_xxyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 197); 

                auto tg_xzz_xxyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 198); 

                auto tg_xzz_xxyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 199); 

                auto tg_xzz_xxzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 200); 

                auto tg_xzz_xyyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 201); 

                auto tg_xzz_xyyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 202); 

                auto tg_xzz_xyyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 203); 

                auto tg_xzz_xyyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 204); 

                auto tg_xzz_xyyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 205); 

                auto tg_xzz_xyzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 206); 

                auto tg_xzz_xzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 207); 

                auto tg_xzz_yyyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 208); 

                auto tg_xzz_yyyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 209); 

                auto tg_xzz_yyyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 210); 

                auto tg_xzz_yyyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 211); 

                auto tg_xzz_yyyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 212); 

                auto tg_xzz_yyzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 213); 

                auto tg_xzz_yzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 214); 

                auto tg_xzz_zzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 215); 

                auto tg_yyy_xxxxxxx_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 216); 

                auto tg_yyy_xxxxxxy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 217); 

                auto tg_yyy_xxxxxxz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 218); 

                auto tg_yyy_xxxxxyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 219); 

                auto tg_yyy_xxxxxyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 220); 

                auto tg_yyy_xxxxxzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 221); 

                auto tg_yyy_xxxxyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 222); 

                auto tg_yyy_xxxxyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 223); 

                auto tg_yyy_xxxxyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 224); 

                auto tg_yyy_xxxxzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 225); 

                auto tg_yyy_xxxyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 226); 

                auto tg_yyy_xxxyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 227); 

                auto tg_yyy_xxxyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 228); 

                auto tg_yyy_xxxyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 229); 

                auto tg_yyy_xxxzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 230); 

                auto tg_yyy_xxyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 231); 

                auto tg_yyy_xxyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 232); 

                auto tg_yyy_xxyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 233); 

                auto tg_yyy_xxyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 234); 

                auto tg_yyy_xxyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 235); 

                auto tg_yyy_xxzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 236); 

                auto tg_yyy_xyyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 237); 

                auto tg_yyy_xyyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 238); 

                auto tg_yyy_xyyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 239); 

                auto tg_yyy_xyyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 240); 

                auto tg_yyy_xyyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 241); 

                auto tg_yyy_xyzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 242); 

                auto tg_yyy_xzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 243); 

                auto tg_yyy_yyyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 244); 

                auto tg_yyy_yyyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 245); 

                auto tg_yyy_yyyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 246); 

                auto tg_yyy_yyyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 247); 

                auto tg_yyy_yyyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 248); 

                auto tg_yyy_yyzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 249); 

                auto tg_yyy_yzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 250); 

                auto tg_yyy_zzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 251); 

                auto tg_yyz_xxxxxxx_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 252); 

                auto tg_yyz_xxxxxxy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 253); 

                auto tg_yyz_xxxxxxz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 254); 

                auto tg_yyz_xxxxxyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 255); 

                auto tg_yyz_xxxxxyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 256); 

                auto tg_yyz_xxxxxzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 257); 

                auto tg_yyz_xxxxyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 258); 

                auto tg_yyz_xxxxyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 259); 

                auto tg_yyz_xxxxyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 260); 

                auto tg_yyz_xxxxzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 261); 

                auto tg_yyz_xxxyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 262); 

                auto tg_yyz_xxxyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 263); 

                auto tg_yyz_xxxyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 264); 

                auto tg_yyz_xxxyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 265); 

                auto tg_yyz_xxxzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 266); 

                auto tg_yyz_xxyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 267); 

                auto tg_yyz_xxyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 268); 

                auto tg_yyz_xxyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 269); 

                auto tg_yyz_xxyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 270); 

                auto tg_yyz_xxyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 271); 

                auto tg_yyz_xxzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 272); 

                auto tg_yyz_xyyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 273); 

                auto tg_yyz_xyyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 274); 

                auto tg_yyz_xyyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 275); 

                auto tg_yyz_xyyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 276); 

                auto tg_yyz_xyyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 277); 

                auto tg_yyz_xyzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 278); 

                auto tg_yyz_xzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 279); 

                auto tg_yyz_yyyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 280); 

                auto tg_yyz_yyyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 281); 

                auto tg_yyz_yyyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 282); 

                auto tg_yyz_yyyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 283); 

                auto tg_yyz_yyyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 284); 

                auto tg_xzz_xxxyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 190); 

                auto tg_xzz_xxxyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 191); 

                auto tg_xzz_xxxyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 192); 

                auto tg_xzz_xxxyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 193); 

                auto tg_xzz_xxxzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 194); 

                auto tg_xzz_xxyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 195); 

                auto tg_xzz_xxyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 196); 

                auto tg_xzz_xxyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 197); 

                auto tg_xzz_xxyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 198); 

                auto tg_xzz_xxyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 199); 

                auto tg_xzz_xxzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 200); 

                auto tg_xzz_xyyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 201); 

                auto tg_xzz_xyyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 202); 

                auto tg_xzz_xyyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 203); 

                auto tg_xzz_xyyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 204); 

                auto tg_xzz_xyyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 205); 

                auto tg_xzz_xyzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 206); 

                auto tg_xzz_xzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 207); 

                auto tg_xzz_yyyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 208); 

                auto tg_xzz_yyyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 209); 

                auto tg_xzz_yyyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 210); 

                auto tg_xzz_yyyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 211); 

                auto tg_xzz_yyyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 212); 

                auto tg_xzz_yyzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 213); 

                auto tg_xzz_yzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 214); 

                auto tg_xzz_zzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 215); 

                auto tg_yyy_xxxxxxx_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 216); 

                auto tg_yyy_xxxxxxy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 217); 

                auto tg_yyy_xxxxxxz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 218); 

                auto tg_yyy_xxxxxyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 219); 

                auto tg_yyy_xxxxxyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 220); 

                auto tg_yyy_xxxxxzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 221); 

                auto tg_yyy_xxxxyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 222); 

                auto tg_yyy_xxxxyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 223); 

                auto tg_yyy_xxxxyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 224); 

                auto tg_yyy_xxxxzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 225); 

                auto tg_yyy_xxxyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 226); 

                auto tg_yyy_xxxyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 227); 

                auto tg_yyy_xxxyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 228); 

                auto tg_yyy_xxxyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 229); 

                auto tg_yyy_xxxzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 230); 

                auto tg_yyy_xxyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 231); 

                auto tg_yyy_xxyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 232); 

                auto tg_yyy_xxyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 233); 

                auto tg_yyy_xxyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 234); 

                auto tg_yyy_xxyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 235); 

                auto tg_yyy_xxzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 236); 

                auto tg_yyy_xyyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 237); 

                auto tg_yyy_xyyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 238); 

                auto tg_yyy_xyyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 239); 

                auto tg_yyy_xyyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 240); 

                auto tg_yyy_xyyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 241); 

                auto tg_yyy_xyzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 242); 

                auto tg_yyy_xzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 243); 

                auto tg_yyy_yyyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 244); 

                auto tg_yyy_yyyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 245); 

                auto tg_yyy_yyyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 246); 

                auto tg_yyy_yyyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 247); 

                auto tg_yyy_yyyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 248); 

                auto tg_yyy_yyzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 249); 

                auto tg_yyy_yzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 250); 

                auto tg_yyy_zzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 251); 

                auto tg_yyz_xxxxxxx_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 252); 

                auto tg_yyz_xxxxxxy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 253); 

                auto tg_yyz_xxxxxxz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 254); 

                auto tg_yyz_xxxxxyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 255); 

                auto tg_yyz_xxxxxyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 256); 

                auto tg_yyz_xxxxxzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 257); 

                auto tg_yyz_xxxxyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 258); 

                auto tg_yyz_xxxxyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 259); 

                auto tg_yyz_xxxxyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 260); 

                auto tg_yyz_xxxxzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 261); 

                auto tg_yyz_xxxyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 262); 

                auto tg_yyz_xxxyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 263); 

                auto tg_yyz_xxxyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 264); 

                auto tg_yyz_xxxyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 265); 

                auto tg_yyz_xxxzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 266); 

                auto tg_yyz_xxyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 267); 

                auto tg_yyz_xxyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 268); 

                auto tg_yyz_xxyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 269); 

                auto tg_yyz_xxyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 270); 

                auto tg_yyz_xxyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 271); 

                auto tg_yyz_xxzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 272); 

                auto tg_yyz_xyyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 273); 

                auto tg_yyz_xyyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 274); 

                auto tg_yyz_xyyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 275); 

                auto tg_yyz_xyyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 276); 

                auto tg_yyz_xyyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 277); 

                auto tg_yyz_xyzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 278); 

                auto tg_yyz_xzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 279); 

                auto tg_yyz_yyyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 280); 

                auto tg_yyz_yyyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 281); 

                auto tg_yyz_yyyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 282); 

                auto tg_yyz_yyyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 283); 

                auto tg_yyz_yyyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 284); 

                auto tg_xxzz_xxyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 150); 

                auto tg_xxzz_xxyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 151); 

                auto tg_xxzz_xxyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 152); 

                auto tg_xxzz_xxyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 153); 

                auto tg_xxzz_xxzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 154); 

                auto tg_xxzz_xyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 155); 

                auto tg_xxzz_xyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 156); 

                auto tg_xxzz_xyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 157); 

                auto tg_xxzz_xyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 158); 

                auto tg_xxzz_xyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 159); 

                auto tg_xxzz_xzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 160); 

                auto tg_xxzz_yyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 161); 

                auto tg_xxzz_yyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 162); 

                auto tg_xxzz_yyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 163); 

                auto tg_xxzz_yyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 164); 

                auto tg_xxzz_yyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 165); 

                auto tg_xxzz_yzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 166); 

                auto tg_xxzz_zzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 167); 

                auto tg_xyyy_xxxxxx_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 168); 

                auto tg_xyyy_xxxxxy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 169); 

                auto tg_xyyy_xxxxxz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 170); 

                auto tg_xyyy_xxxxyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 171); 

                auto tg_xyyy_xxxxyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 172); 

                auto tg_xyyy_xxxxzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 173); 

                auto tg_xyyy_xxxyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 174); 

                auto tg_xyyy_xxxyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 175); 

                auto tg_xyyy_xxxyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 176); 

                auto tg_xyyy_xxxzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 177); 

                auto tg_xyyy_xxyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 178); 

                auto tg_xyyy_xxyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 179); 

                auto tg_xyyy_xxyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 180); 

                auto tg_xyyy_xxyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 181); 

                auto tg_xyyy_xxzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 182); 

                auto tg_xyyy_xyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 183); 

                auto tg_xyyy_xyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 184); 

                auto tg_xyyy_xyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 185); 

                auto tg_xyyy_xyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 186); 

                auto tg_xyyy_xyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 187); 

                auto tg_xyyy_xzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 188); 

                auto tg_xyyy_yyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 189); 

                auto tg_xyyy_yyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 190); 

                auto tg_xyyy_yyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 191); 

                auto tg_xyyy_yyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 192); 

                auto tg_xyyy_yyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 193); 

                auto tg_xyyy_yzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 194); 

                auto tg_xyyy_zzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 195); 

                auto tg_xyyz_xxxxxx_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 196); 

                auto tg_xyyz_xxxxxy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 197); 

                auto tg_xyyz_xxxxxz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 198); 

                auto tg_xyyz_xxxxyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 199); 

                auto tg_xyyz_xxxxyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 200); 

                auto tg_xyyz_xxxxzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 201); 

                auto tg_xyyz_xxxyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 202); 

                auto tg_xyyz_xxxyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 203); 

                auto tg_xyyz_xxxyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 204); 

                auto tg_xyyz_xxxzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 205); 

                auto tg_xyyz_xxyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 206); 

                auto tg_xyyz_xxyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 207); 

                auto tg_xyyz_xxyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 208); 

                auto tg_xyyz_xxyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 209); 

                auto tg_xyyz_xxzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 210); 

                auto tg_xyyz_xyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 211); 

                auto tg_xyyz_xyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 212); 

                auto tg_xyyz_xyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 213); 

                auto tg_xyyz_xyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 214); 

                auto tg_xyyz_xyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 215); 

                auto tg_xyyz_xzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 216); 

                auto tg_xyyz_yyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 217); 

                auto tg_xyyz_yyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 218); 

                auto tg_xyyz_yyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 219); 

                auto tg_xyyz_yyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 220); 

                auto tg_xyyz_yyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 221); 

                auto tg_xyyz_yzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 222); 

                auto tg_xyyz_zzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 223); 

                // set up pointers to integrals

                auto tg_xxxzz_xxxyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 190); 

                auto tg_xxxzz_xxxyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 191); 

                auto tg_xxxzz_xxxyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 192); 

                auto tg_xxxzz_xxxyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 193); 

                auto tg_xxxzz_xxxzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 194); 

                auto tg_xxxzz_xxyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 195); 

                auto tg_xxxzz_xxyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 196); 

                auto tg_xxxzz_xxyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 197); 

                auto tg_xxxzz_xxyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 198); 

                auto tg_xxxzz_xxyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 199); 

                auto tg_xxxzz_xxzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 200); 

                auto tg_xxxzz_xyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 201); 

                auto tg_xxxzz_xyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 202); 

                auto tg_xxxzz_xyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 203); 

                auto tg_xxxzz_xyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 204); 

                auto tg_xxxzz_xyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 205); 

                auto tg_xxxzz_xyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 206); 

                auto tg_xxxzz_xzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 207); 

                auto tg_xxxzz_yyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 208); 

                auto tg_xxxzz_yyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 209); 

                auto tg_xxxzz_yyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 210); 

                auto tg_xxxzz_yyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 211); 

                auto tg_xxxzz_yyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 212); 

                auto tg_xxxzz_yyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 213); 

                auto tg_xxxzz_yzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 214); 

                auto tg_xxxzz_zzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 215); 

                auto tg_xxyyy_xxxxxxx_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 216); 

                auto tg_xxyyy_xxxxxxy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 217); 

                auto tg_xxyyy_xxxxxxz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 218); 

                auto tg_xxyyy_xxxxxyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 219); 

                auto tg_xxyyy_xxxxxyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 220); 

                auto tg_xxyyy_xxxxxzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 221); 

                auto tg_xxyyy_xxxxyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 222); 

                auto tg_xxyyy_xxxxyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 223); 

                auto tg_xxyyy_xxxxyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 224); 

                auto tg_xxyyy_xxxxzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 225); 

                auto tg_xxyyy_xxxyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 226); 

                auto tg_xxyyy_xxxyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 227); 

                auto tg_xxyyy_xxxyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 228); 

                auto tg_xxyyy_xxxyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 229); 

                auto tg_xxyyy_xxxzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 230); 

                auto tg_xxyyy_xxyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 231); 

                auto tg_xxyyy_xxyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 232); 

                auto tg_xxyyy_xxyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 233); 

                auto tg_xxyyy_xxyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 234); 

                auto tg_xxyyy_xxyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 235); 

                auto tg_xxyyy_xxzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 236); 

                auto tg_xxyyy_xyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 237); 

                auto tg_xxyyy_xyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 238); 

                auto tg_xxyyy_xyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 239); 

                auto tg_xxyyy_xyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 240); 

                auto tg_xxyyy_xyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 241); 

                auto tg_xxyyy_xyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 242); 

                auto tg_xxyyy_xzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 243); 

                auto tg_xxyyy_yyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 244); 

                auto tg_xxyyy_yyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 245); 

                auto tg_xxyyy_yyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 246); 

                auto tg_xxyyy_yyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 247); 

                auto tg_xxyyy_yyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 248); 

                auto tg_xxyyy_yyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 249); 

                auto tg_xxyyy_yzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 250); 

                auto tg_xxyyy_zzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 251); 

                auto tg_xxyyz_xxxxxxx_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 252); 

                auto tg_xxyyz_xxxxxxy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 253); 

                auto tg_xxyyz_xxxxxxz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 254); 

                auto tg_xxyyz_xxxxxyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 255); 

                auto tg_xxyyz_xxxxxyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 256); 

                auto tg_xxyyz_xxxxxzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 257); 

                auto tg_xxyyz_xxxxyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 258); 

                auto tg_xxyyz_xxxxyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 259); 

                auto tg_xxyyz_xxxxyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 260); 

                auto tg_xxyyz_xxxxzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 261); 

                auto tg_xxyyz_xxxyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 262); 

                auto tg_xxyyz_xxxyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 263); 

                auto tg_xxyyz_xxxyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 264); 

                auto tg_xxyyz_xxxyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 265); 

                auto tg_xxyyz_xxxzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 266); 

                auto tg_xxyyz_xxyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 267); 

                auto tg_xxyyz_xxyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 268); 

                auto tg_xxyyz_xxyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 269); 

                auto tg_xxyyz_xxyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 270); 

                auto tg_xxyyz_xxyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 271); 

                auto tg_xxyyz_xxzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 272); 

                auto tg_xxyyz_xyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 273); 

                auto tg_xxyyz_xyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 274); 

                auto tg_xxyyz_xyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 275); 

                auto tg_xxyyz_xyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 276); 

                auto tg_xxyyz_xyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 277); 

                auto tg_xxyyz_xyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 278); 

                auto tg_xxyyz_xzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 279); 

                auto tg_xxyyz_yyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 280); 

                auto tg_xxyyz_yyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 281); 

                auto tg_xxyyz_yyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 282); 

                auto tg_xxyyz_yyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 283); 

                auto tg_xxyyz_yyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 284); 

                // Batch of Integrals (190,285)

                #pragma omp simd aligned(fxn, fza, tg_xxxzz_xxxyyyy_0, tg_xxxzz_xxxyyyz_0, tg_xxxzz_xxxyyzz_0, \
                                         tg_xxxzz_xxxyzzz_0, tg_xxxzz_xxxzzzz_0, tg_xxxzz_xxyyyyy_0, tg_xxxzz_xxyyyyz_0, \
                                         tg_xxxzz_xxyyyzz_0, tg_xxxzz_xxyyzzz_0, tg_xxxzz_xxyzzzz_0, tg_xxxzz_xxzzzzz_0, \
                                         tg_xxxzz_xyyyyyy_0, tg_xxxzz_xyyyyyz_0, tg_xxxzz_xyyyyzz_0, tg_xxxzz_xyyyzzz_0, \
                                         tg_xxxzz_xyyzzzz_0, tg_xxxzz_xyzzzzz_0, tg_xxxzz_xzzzzzz_0, tg_xxxzz_yyyyyyy_0, \
                                         tg_xxxzz_yyyyyyz_0, tg_xxxzz_yyyyyzz_0, tg_xxxzz_yyyyzzz_0, tg_xxxzz_yyyzzzz_0, \
                                         tg_xxxzz_yyzzzzz_0, tg_xxxzz_yzzzzzz_0, tg_xxxzz_zzzzzzz_0, tg_xxyyy_xxxxxxx_0, \
                                         tg_xxyyy_xxxxxxy_0, tg_xxyyy_xxxxxxz_0, tg_xxyyy_xxxxxyy_0, tg_xxyyy_xxxxxyz_0, \
                                         tg_xxyyy_xxxxxzz_0, tg_xxyyy_xxxxyyy_0, tg_xxyyy_xxxxyyz_0, tg_xxyyy_xxxxyzz_0, \
                                         tg_xxyyy_xxxxzzz_0, tg_xxyyy_xxxyyyy_0, tg_xxyyy_xxxyyyz_0, tg_xxyyy_xxxyyzz_0, \
                                         tg_xxyyy_xxxyzzz_0, tg_xxyyy_xxxzzzz_0, tg_xxyyy_xxyyyyy_0, tg_xxyyy_xxyyyyz_0, \
                                         tg_xxyyy_xxyyyzz_0, tg_xxyyy_xxyyzzz_0, tg_xxyyy_xxyzzzz_0, tg_xxyyy_xxzzzzz_0, \
                                         tg_xxyyy_xyyyyyy_0, tg_xxyyy_xyyyyyz_0, tg_xxyyy_xyyyyzz_0, tg_xxyyy_xyyyzzz_0, \
                                         tg_xxyyy_xyyzzzz_0, tg_xxyyy_xyzzzzz_0, tg_xxyyy_xzzzzzz_0, tg_xxyyy_yyyyyyy_0, \
                                         tg_xxyyy_yyyyyyz_0, tg_xxyyy_yyyyyzz_0, tg_xxyyy_yyyyzzz_0, tg_xxyyy_yyyzzzz_0, \
                                         tg_xxyyy_yyzzzzz_0, tg_xxyyy_yzzzzzz_0, tg_xxyyy_zzzzzzz_0, tg_xxyyz_xxxxxxx_0, \
                                         tg_xxyyz_xxxxxxy_0, tg_xxyyz_xxxxxxz_0, tg_xxyyz_xxxxxyy_0, tg_xxyyz_xxxxxyz_0, \
                                         tg_xxyyz_xxxxxzz_0, tg_xxyyz_xxxxyyy_0, tg_xxyyz_xxxxyyz_0, tg_xxyyz_xxxxyzz_0, \
                                         tg_xxyyz_xxxxzzz_0, tg_xxyyz_xxxyyyy_0, tg_xxyyz_xxxyyyz_0, tg_xxyyz_xxxyyzz_0, \
                                         tg_xxyyz_xxxyzzz_0, tg_xxyyz_xxxzzzz_0, tg_xxyyz_xxyyyyy_0, tg_xxyyz_xxyyyyz_0, \
                                         tg_xxyyz_xxyyyzz_0, tg_xxyyz_xxyyzzz_0, tg_xxyyz_xxyzzzz_0, tg_xxyyz_xxzzzzz_0, \
                                         tg_xxyyz_xyyyyyy_0, tg_xxyyz_xyyyyyz_0, tg_xxyyz_xyyyyzz_0, tg_xxyyz_xyyyzzz_0, \
                                         tg_xxyyz_xyyzzzz_0, tg_xxyyz_xyzzzzz_0, tg_xxyyz_xzzzzzz_0, tg_xxyyz_yyyyyyy_0, \
                                         tg_xxyyz_yyyyyyz_0, tg_xxyyz_yyyyyzz_0, tg_xxyyz_yyyyzzz_0, tg_xxyyz_yyyzzzz_0, \
                                         tg_xxzz_xxxyyyy_0, tg_xxzz_xxxyyyy_1, tg_xxzz_xxxyyyz_0, tg_xxzz_xxxyyyz_1, \
                                         tg_xxzz_xxxyyzz_0, tg_xxzz_xxxyyzz_1, tg_xxzz_xxxyzzz_0, tg_xxzz_xxxyzzz_1, \
                                         tg_xxzz_xxxzzzz_0, tg_xxzz_xxxzzzz_1, tg_xxzz_xxyyyy_1, tg_xxzz_xxyyyyy_0, \
                                         tg_xxzz_xxyyyyy_1, tg_xxzz_xxyyyyz_0, tg_xxzz_xxyyyyz_1, tg_xxzz_xxyyyz_1, \
                                         tg_xxzz_xxyyyzz_0, tg_xxzz_xxyyyzz_1, tg_xxzz_xxyyzz_1, tg_xxzz_xxyyzzz_0, \
                                         tg_xxzz_xxyyzzz_1, tg_xxzz_xxyzzz_1, tg_xxzz_xxyzzzz_0, tg_xxzz_xxyzzzz_1, \
                                         tg_xxzz_xxzzzz_1, tg_xxzz_xxzzzzz_0, tg_xxzz_xxzzzzz_1, tg_xxzz_xyyyyy_1, \
                                         tg_xxzz_xyyyyyy_0, tg_xxzz_xyyyyyy_1, tg_xxzz_xyyyyyz_0, tg_xxzz_xyyyyyz_1, \
                                         tg_xxzz_xyyyyz_1, tg_xxzz_xyyyyzz_0, tg_xxzz_xyyyyzz_1, tg_xxzz_xyyyzz_1, \
                                         tg_xxzz_xyyyzzz_0, tg_xxzz_xyyyzzz_1, tg_xxzz_xyyzzz_1, tg_xxzz_xyyzzzz_0, \
                                         tg_xxzz_xyyzzzz_1, tg_xxzz_xyzzzz_1, tg_xxzz_xyzzzzz_0, tg_xxzz_xyzzzzz_1, \
                                         tg_xxzz_xzzzzz_1, tg_xxzz_xzzzzzz_0, tg_xxzz_xzzzzzz_1, tg_xxzz_yyyyyy_1, \
                                         tg_xxzz_yyyyyyy_0, tg_xxzz_yyyyyyy_1, tg_xxzz_yyyyyyz_0, tg_xxzz_yyyyyyz_1, \
                                         tg_xxzz_yyyyyz_1, tg_xxzz_yyyyyzz_0, tg_xxzz_yyyyyzz_1, tg_xxzz_yyyyzz_1, \
                                         tg_xxzz_yyyyzzz_0, tg_xxzz_yyyyzzz_1, tg_xxzz_yyyzzz_1, tg_xxzz_yyyzzzz_0, \
                                         tg_xxzz_yyyzzzz_1, tg_xxzz_yyzzzz_1, tg_xxzz_yyzzzzz_0, tg_xxzz_yyzzzzz_1, \
                                         tg_xxzz_yzzzzz_1, tg_xxzz_yzzzzzz_0, tg_xxzz_yzzzzzz_1, tg_xxzz_zzzzzz_1, \
                                         tg_xxzz_zzzzzzz_0, tg_xxzz_zzzzzzz_1, tg_xyyy_xxxxxx_1, tg_xyyy_xxxxxxx_0, \
                                         tg_xyyy_xxxxxxx_1, tg_xyyy_xxxxxxy_0, tg_xyyy_xxxxxxy_1, tg_xyyy_xxxxxxz_0, \
                                         tg_xyyy_xxxxxxz_1, tg_xyyy_xxxxxy_1, tg_xyyy_xxxxxyy_0, tg_xyyy_xxxxxyy_1, \
                                         tg_xyyy_xxxxxyz_0, tg_xyyy_xxxxxyz_1, tg_xyyy_xxxxxz_1, tg_xyyy_xxxxxzz_0, \
                                         tg_xyyy_xxxxxzz_1, tg_xyyy_xxxxyy_1, tg_xyyy_xxxxyyy_0, tg_xyyy_xxxxyyy_1, \
                                         tg_xyyy_xxxxyyz_0, tg_xyyy_xxxxyyz_1, tg_xyyy_xxxxyz_1, tg_xyyy_xxxxyzz_0, \
                                         tg_xyyy_xxxxyzz_1, tg_xyyy_xxxxzz_1, tg_xyyy_xxxxzzz_0, tg_xyyy_xxxxzzz_1, \
                                         tg_xyyy_xxxyyy_1, tg_xyyy_xxxyyyy_0, tg_xyyy_xxxyyyy_1, tg_xyyy_xxxyyyz_0, \
                                         tg_xyyy_xxxyyyz_1, tg_xyyy_xxxyyz_1, tg_xyyy_xxxyyzz_0, tg_xyyy_xxxyyzz_1, \
                                         tg_xyyy_xxxyzz_1, tg_xyyy_xxxyzzz_0, tg_xyyy_xxxyzzz_1, tg_xyyy_xxxzzz_1, \
                                         tg_xyyy_xxxzzzz_0, tg_xyyy_xxxzzzz_1, tg_xyyy_xxyyyy_1, tg_xyyy_xxyyyyy_0, \
                                         tg_xyyy_xxyyyyy_1, tg_xyyy_xxyyyyz_0, tg_xyyy_xxyyyyz_1, tg_xyyy_xxyyyz_1, \
                                         tg_xyyy_xxyyyzz_0, tg_xyyy_xxyyyzz_1, tg_xyyy_xxyyzz_1, tg_xyyy_xxyyzzz_0, \
                                         tg_xyyy_xxyyzzz_1, tg_xyyy_xxyzzz_1, tg_xyyy_xxyzzzz_0, tg_xyyy_xxyzzzz_1, \
                                         tg_xyyy_xxzzzz_1, tg_xyyy_xxzzzzz_0, tg_xyyy_xxzzzzz_1, tg_xyyy_xyyyyy_1, \
                                         tg_xyyy_xyyyyyy_0, tg_xyyy_xyyyyyy_1, tg_xyyy_xyyyyyz_0, tg_xyyy_xyyyyyz_1, \
                                         tg_xyyy_xyyyyz_1, tg_xyyy_xyyyyzz_0, tg_xyyy_xyyyyzz_1, tg_xyyy_xyyyzz_1, \
                                         tg_xyyy_xyyyzzz_0, tg_xyyy_xyyyzzz_1, tg_xyyy_xyyzzz_1, tg_xyyy_xyyzzzz_0, \
                                         tg_xyyy_xyyzzzz_1, tg_xyyy_xyzzzz_1, tg_xyyy_xyzzzzz_0, tg_xyyy_xyzzzzz_1, \
                                         tg_xyyy_xzzzzz_1, tg_xyyy_xzzzzzz_0, tg_xyyy_xzzzzzz_1, tg_xyyy_yyyyyy_1, \
                                         tg_xyyy_yyyyyyy_0, tg_xyyy_yyyyyyy_1, tg_xyyy_yyyyyyz_0, tg_xyyy_yyyyyyz_1, \
                                         tg_xyyy_yyyyyz_1, tg_xyyy_yyyyyzz_0, tg_xyyy_yyyyyzz_1, tg_xyyy_yyyyzz_1, \
                                         tg_xyyy_yyyyzzz_0, tg_xyyy_yyyyzzz_1, tg_xyyy_yyyzzz_1, tg_xyyy_yyyzzzz_0, \
                                         tg_xyyy_yyyzzzz_1, tg_xyyy_yyzzzz_1, tg_xyyy_yyzzzzz_0, tg_xyyy_yyzzzzz_1, \
                                         tg_xyyy_yzzzzz_1, tg_xyyy_yzzzzzz_0, tg_xyyy_yzzzzzz_1, tg_xyyy_zzzzzz_1, \
                                         tg_xyyy_zzzzzzz_0, tg_xyyy_zzzzzzz_1, tg_xyyz_xxxxxx_1, tg_xyyz_xxxxxxx_0, \
                                         tg_xyyz_xxxxxxx_1, tg_xyyz_xxxxxxy_0, tg_xyyz_xxxxxxy_1, tg_xyyz_xxxxxxz_0, \
                                         tg_xyyz_xxxxxxz_1, tg_xyyz_xxxxxy_1, tg_xyyz_xxxxxyy_0, tg_xyyz_xxxxxyy_1, \
                                         tg_xyyz_xxxxxyz_0, tg_xyyz_xxxxxyz_1, tg_xyyz_xxxxxz_1, tg_xyyz_xxxxxzz_0, \
                                         tg_xyyz_xxxxxzz_1, tg_xyyz_xxxxyy_1, tg_xyyz_xxxxyyy_0, tg_xyyz_xxxxyyy_1, \
                                         tg_xyyz_xxxxyyz_0, tg_xyyz_xxxxyyz_1, tg_xyyz_xxxxyz_1, tg_xyyz_xxxxyzz_0, \
                                         tg_xyyz_xxxxyzz_1, tg_xyyz_xxxxzz_1, tg_xyyz_xxxxzzz_0, tg_xyyz_xxxxzzz_1, \
                                         tg_xyyz_xxxyyy_1, tg_xyyz_xxxyyyy_0, tg_xyyz_xxxyyyy_1, tg_xyyz_xxxyyyz_0, \
                                         tg_xyyz_xxxyyyz_1, tg_xyyz_xxxyyz_1, tg_xyyz_xxxyyzz_0, tg_xyyz_xxxyyzz_1, \
                                         tg_xyyz_xxxyzz_1, tg_xyyz_xxxyzzz_0, tg_xyyz_xxxyzzz_1, tg_xyyz_xxxzzz_1, \
                                         tg_xyyz_xxxzzzz_0, tg_xyyz_xxxzzzz_1, tg_xyyz_xxyyyy_1, tg_xyyz_xxyyyyy_0, \
                                         tg_xyyz_xxyyyyy_1, tg_xyyz_xxyyyyz_0, tg_xyyz_xxyyyyz_1, tg_xyyz_xxyyyz_1, \
                                         tg_xyyz_xxyyyzz_0, tg_xyyz_xxyyyzz_1, tg_xyyz_xxyyzz_1, tg_xyyz_xxyyzzz_0, \
                                         tg_xyyz_xxyyzzz_1, tg_xyyz_xxyzzz_1, tg_xyyz_xxyzzzz_0, tg_xyyz_xxyzzzz_1, \
                                         tg_xyyz_xxzzzz_1, tg_xyyz_xxzzzzz_0, tg_xyyz_xxzzzzz_1, tg_xyyz_xyyyyy_1, \
                                         tg_xyyz_xyyyyyy_0, tg_xyyz_xyyyyyy_1, tg_xyyz_xyyyyyz_0, tg_xyyz_xyyyyyz_1, \
                                         tg_xyyz_xyyyyz_1, tg_xyyz_xyyyyzz_0, tg_xyyz_xyyyyzz_1, tg_xyyz_xyyyzz_1, \
                                         tg_xyyz_xyyyzzz_0, tg_xyyz_xyyyzzz_1, tg_xyyz_xyyzzz_1, tg_xyyz_xyyzzzz_0, \
                                         tg_xyyz_xyyzzzz_1, tg_xyyz_xyzzzz_1, tg_xyyz_xyzzzzz_0, tg_xyyz_xyzzzzz_1, \
                                         tg_xyyz_xzzzzz_1, tg_xyyz_xzzzzzz_0, tg_xyyz_xzzzzzz_1, tg_xyyz_yyyyyy_1, \
                                         tg_xyyz_yyyyyyy_0, tg_xyyz_yyyyyyy_1, tg_xyyz_yyyyyyz_0, tg_xyyz_yyyyyyz_1, \
                                         tg_xyyz_yyyyyz_1, tg_xyyz_yyyyyzz_0, tg_xyyz_yyyyyzz_1, tg_xyyz_yyyyzz_1, \
                                         tg_xyyz_yyyyzzz_0, tg_xyyz_yyyyzzz_1, tg_xyyz_yyyzzz_1, tg_xyyz_yyyzzzz_0, \
                                         tg_xyyz_yyyzzzz_1, tg_xyyz_yyzzzz_1, tg_xyyz_yzzzzz_1, tg_xyyz_zzzzzz_1, \
                                         tg_xzz_xxxyyyy_0, tg_xzz_xxxyyyy_1, tg_xzz_xxxyyyz_0, tg_xzz_xxxyyyz_1, \
                                         tg_xzz_xxxyyzz_0, tg_xzz_xxxyyzz_1, tg_xzz_xxxyzzz_0, tg_xzz_xxxyzzz_1, \
                                         tg_xzz_xxxzzzz_0, tg_xzz_xxxzzzz_1, tg_xzz_xxyyyyy_0, tg_xzz_xxyyyyy_1, \
                                         tg_xzz_xxyyyyz_0, tg_xzz_xxyyyyz_1, tg_xzz_xxyyyzz_0, tg_xzz_xxyyyzz_1, \
                                         tg_xzz_xxyyzzz_0, tg_xzz_xxyyzzz_1, tg_xzz_xxyzzzz_0, tg_xzz_xxyzzzz_1, \
                                         tg_xzz_xxzzzzz_0, tg_xzz_xxzzzzz_1, tg_xzz_xyyyyyy_0, tg_xzz_xyyyyyy_1, \
                                         tg_xzz_xyyyyyz_0, tg_xzz_xyyyyyz_1, tg_xzz_xyyyyzz_0, tg_xzz_xyyyyzz_1, \
                                         tg_xzz_xyyyzzz_0, tg_xzz_xyyyzzz_1, tg_xzz_xyyzzzz_0, tg_xzz_xyyzzzz_1, \
                                         tg_xzz_xyzzzzz_0, tg_xzz_xyzzzzz_1, tg_xzz_xzzzzzz_0, tg_xzz_xzzzzzz_1, \
                                         tg_xzz_yyyyyyy_0, tg_xzz_yyyyyyy_1, tg_xzz_yyyyyyz_0, tg_xzz_yyyyyyz_1, \
                                         tg_xzz_yyyyyzz_0, tg_xzz_yyyyyzz_1, tg_xzz_yyyyzzz_0, tg_xzz_yyyyzzz_1, \
                                         tg_xzz_yyyzzzz_0, tg_xzz_yyyzzzz_1, tg_xzz_yyzzzzz_0, tg_xzz_yyzzzzz_1, \
                                         tg_xzz_yzzzzzz_0, tg_xzz_yzzzzzz_1, tg_xzz_zzzzzzz_0, tg_xzz_zzzzzzz_1, \
                                         tg_yyy_xxxxxxx_0, tg_yyy_xxxxxxx_1, tg_yyy_xxxxxxy_0, tg_yyy_xxxxxxy_1, \
                                         tg_yyy_xxxxxxz_0, tg_yyy_xxxxxxz_1, tg_yyy_xxxxxyy_0, tg_yyy_xxxxxyy_1, \
                                         tg_yyy_xxxxxyz_0, tg_yyy_xxxxxyz_1, tg_yyy_xxxxxzz_0, tg_yyy_xxxxxzz_1, \
                                         tg_yyy_xxxxyyy_0, tg_yyy_xxxxyyy_1, tg_yyy_xxxxyyz_0, tg_yyy_xxxxyyz_1, \
                                         tg_yyy_xxxxyzz_0, tg_yyy_xxxxyzz_1, tg_yyy_xxxxzzz_0, tg_yyy_xxxxzzz_1, \
                                         tg_yyy_xxxyyyy_0, tg_yyy_xxxyyyy_1, tg_yyy_xxxyyyz_0, tg_yyy_xxxyyyz_1, \
                                         tg_yyy_xxxyyzz_0, tg_yyy_xxxyyzz_1, tg_yyy_xxxyzzz_0, tg_yyy_xxxyzzz_1, \
                                         tg_yyy_xxxzzzz_0, tg_yyy_xxxzzzz_1, tg_yyy_xxyyyyy_0, tg_yyy_xxyyyyy_1, \
                                         tg_yyy_xxyyyyz_0, tg_yyy_xxyyyyz_1, tg_yyy_xxyyyzz_0, tg_yyy_xxyyyzz_1, \
                                         tg_yyy_xxyyzzz_0, tg_yyy_xxyyzzz_1, tg_yyy_xxyzzzz_0, tg_yyy_xxyzzzz_1, \
                                         tg_yyy_xxzzzzz_0, tg_yyy_xxzzzzz_1, tg_yyy_xyyyyyy_0, tg_yyy_xyyyyyy_1, \
                                         tg_yyy_xyyyyyz_0, tg_yyy_xyyyyyz_1, tg_yyy_xyyyyzz_0, tg_yyy_xyyyyzz_1, \
                                         tg_yyy_xyyyzzz_0, tg_yyy_xyyyzzz_1, tg_yyy_xyyzzzz_0, tg_yyy_xyyzzzz_1, \
                                         tg_yyy_xyzzzzz_0, tg_yyy_xyzzzzz_1, tg_yyy_xzzzzzz_0, tg_yyy_xzzzzzz_1, \
                                         tg_yyy_yyyyyyy_0, tg_yyy_yyyyyyy_1, tg_yyy_yyyyyyz_0, tg_yyy_yyyyyyz_1, \
                                         tg_yyy_yyyyyzz_0, tg_yyy_yyyyyzz_1, tg_yyy_yyyyzzz_0, tg_yyy_yyyyzzz_1, \
                                         tg_yyy_yyyzzzz_0, tg_yyy_yyyzzzz_1, tg_yyy_yyzzzzz_0, tg_yyy_yyzzzzz_1, \
                                         tg_yyy_yzzzzzz_0, tg_yyy_yzzzzzz_1, tg_yyy_zzzzzzz_0, tg_yyy_zzzzzzz_1, \
                                         tg_yyz_xxxxxxx_0, tg_yyz_xxxxxxx_1, tg_yyz_xxxxxxy_0, tg_yyz_xxxxxxy_1, \
                                         tg_yyz_xxxxxxz_0, tg_yyz_xxxxxxz_1, tg_yyz_xxxxxyy_0, tg_yyz_xxxxxyy_1, \
                                         tg_yyz_xxxxxyz_0, tg_yyz_xxxxxyz_1, tg_yyz_xxxxxzz_0, tg_yyz_xxxxxzz_1, \
                                         tg_yyz_xxxxyyy_0, tg_yyz_xxxxyyy_1, tg_yyz_xxxxyyz_0, tg_yyz_xxxxyyz_1, \
                                         tg_yyz_xxxxyzz_0, tg_yyz_xxxxyzz_1, tg_yyz_xxxxzzz_0, tg_yyz_xxxxzzz_1, \
                                         tg_yyz_xxxyyyy_0, tg_yyz_xxxyyyy_1, tg_yyz_xxxyyyz_0, tg_yyz_xxxyyyz_1, \
                                         tg_yyz_xxxyyzz_0, tg_yyz_xxxyyzz_1, tg_yyz_xxxyzzz_0, tg_yyz_xxxyzzz_1, \
                                         tg_yyz_xxxzzzz_0, tg_yyz_xxxzzzz_1, tg_yyz_xxyyyyy_0, tg_yyz_xxyyyyy_1, \
                                         tg_yyz_xxyyyyz_0, tg_yyz_xxyyyyz_1, tg_yyz_xxyyyzz_0, tg_yyz_xxyyyzz_1, \
                                         tg_yyz_xxyyzzz_0, tg_yyz_xxyyzzz_1, tg_yyz_xxyzzzz_0, tg_yyz_xxyzzzz_1, \
                                         tg_yyz_xxzzzzz_0, tg_yyz_xxzzzzz_1, tg_yyz_xyyyyyy_0, tg_yyz_xyyyyyy_1, \
                                         tg_yyz_xyyyyyz_0, tg_yyz_xyyyyyz_1, tg_yyz_xyyyyzz_0, tg_yyz_xyyyyzz_1, \
                                         tg_yyz_xyyyzzz_0, tg_yyz_xyyyzzz_1, tg_yyz_xyyzzzz_0, tg_yyz_xyyzzzz_1, \
                                         tg_yyz_xyzzzzz_0, tg_yyz_xyzzzzz_1, tg_yyz_xzzzzzz_0, tg_yyz_xzzzzzz_1, \
                                         tg_yyz_yyyyyyy_0, tg_yyz_yyyyyyy_1, tg_yyz_yyyyyyz_0, tg_yyz_yyyyyyz_1, \
                                         tg_yyz_yyyyyzz_0, tg_yyz_yyyyyzz_1, tg_yyz_yyyyzzz_0, tg_yyz_yyyyzzz_1, \
                                         tg_yyz_yyyzzzz_0, tg_yyz_yyyzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxxzz_xxxyyyy_0[j] = pb_x * tg_xxzz_xxxyyyy_0[j] + wp_x[j] * tg_xxzz_xxxyyyy_1[j] + fl1_fx * tg_xzz_xxxyyyy_0[j] - fl1_fx * fl1_fza * tg_xzz_xxxyyyy_1[j] + 1.5 * fl1_fxn * tg_xxzz_xxyyyy_1[j];

                    tg_xxxzz_xxxyyyz_0[j] = pb_x * tg_xxzz_xxxyyyz_0[j] + wp_x[j] * tg_xxzz_xxxyyyz_1[j] + fl1_fx * tg_xzz_xxxyyyz_0[j] - fl1_fx * fl1_fza * tg_xzz_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_xxzz_xxyyyz_1[j];

                    tg_xxxzz_xxxyyzz_0[j] = pb_x * tg_xxzz_xxxyyzz_0[j] + wp_x[j] * tg_xxzz_xxxyyzz_1[j] + fl1_fx * tg_xzz_xxxyyzz_0[j] - fl1_fx * fl1_fza * tg_xzz_xxxyyzz_1[j] + 1.5 * fl1_fxn * tg_xxzz_xxyyzz_1[j];

                    tg_xxxzz_xxxyzzz_0[j] = pb_x * tg_xxzz_xxxyzzz_0[j] + wp_x[j] * tg_xxzz_xxxyzzz_1[j] + fl1_fx * tg_xzz_xxxyzzz_0[j] - fl1_fx * fl1_fza * tg_xzz_xxxyzzz_1[j] + 1.5 * fl1_fxn * tg_xxzz_xxyzzz_1[j];

                    tg_xxxzz_xxxzzzz_0[j] = pb_x * tg_xxzz_xxxzzzz_0[j] + wp_x[j] * tg_xxzz_xxxzzzz_1[j] + fl1_fx * tg_xzz_xxxzzzz_0[j] - fl1_fx * fl1_fza * tg_xzz_xxxzzzz_1[j] + 1.5 * fl1_fxn * tg_xxzz_xxzzzz_1[j];

                    tg_xxxzz_xxyyyyy_0[j] = pb_x * tg_xxzz_xxyyyyy_0[j] + wp_x[j] * tg_xxzz_xxyyyyy_1[j] + fl1_fx * tg_xzz_xxyyyyy_0[j] - fl1_fx * fl1_fza * tg_xzz_xxyyyyy_1[j] + fl1_fxn * tg_xxzz_xyyyyy_1[j];

                    tg_xxxzz_xxyyyyz_0[j] = pb_x * tg_xxzz_xxyyyyz_0[j] + wp_x[j] * tg_xxzz_xxyyyyz_1[j] + fl1_fx * tg_xzz_xxyyyyz_0[j] - fl1_fx * fl1_fza * tg_xzz_xxyyyyz_1[j] + fl1_fxn * tg_xxzz_xyyyyz_1[j];

                    tg_xxxzz_xxyyyzz_0[j] = pb_x * tg_xxzz_xxyyyzz_0[j] + wp_x[j] * tg_xxzz_xxyyyzz_1[j] + fl1_fx * tg_xzz_xxyyyzz_0[j] - fl1_fx * fl1_fza * tg_xzz_xxyyyzz_1[j] + fl1_fxn * tg_xxzz_xyyyzz_1[j];

                    tg_xxxzz_xxyyzzz_0[j] = pb_x * tg_xxzz_xxyyzzz_0[j] + wp_x[j] * tg_xxzz_xxyyzzz_1[j] + fl1_fx * tg_xzz_xxyyzzz_0[j] - fl1_fx * fl1_fza * tg_xzz_xxyyzzz_1[j] + fl1_fxn * tg_xxzz_xyyzzz_1[j];

                    tg_xxxzz_xxyzzzz_0[j] = pb_x * tg_xxzz_xxyzzzz_0[j] + wp_x[j] * tg_xxzz_xxyzzzz_1[j] + fl1_fx * tg_xzz_xxyzzzz_0[j] - fl1_fx * fl1_fza * tg_xzz_xxyzzzz_1[j] + fl1_fxn * tg_xxzz_xyzzzz_1[j];

                    tg_xxxzz_xxzzzzz_0[j] = pb_x * tg_xxzz_xxzzzzz_0[j] + wp_x[j] * tg_xxzz_xxzzzzz_1[j] + fl1_fx * tg_xzz_xxzzzzz_0[j] - fl1_fx * fl1_fza * tg_xzz_xxzzzzz_1[j] + fl1_fxn * tg_xxzz_xzzzzz_1[j];

                    tg_xxxzz_xyyyyyy_0[j] = pb_x * tg_xxzz_xyyyyyy_0[j] + wp_x[j] * tg_xxzz_xyyyyyy_1[j] + fl1_fx * tg_xzz_xyyyyyy_0[j] - fl1_fx * fl1_fza * tg_xzz_xyyyyyy_1[j] + 0.5 * fl1_fxn * tg_xxzz_yyyyyy_1[j];

                    tg_xxxzz_xyyyyyz_0[j] = pb_x * tg_xxzz_xyyyyyz_0[j] + wp_x[j] * tg_xxzz_xyyyyyz_1[j] + fl1_fx * tg_xzz_xyyyyyz_0[j] - fl1_fx * fl1_fza * tg_xzz_xyyyyyz_1[j] + 0.5 * fl1_fxn * tg_xxzz_yyyyyz_1[j];

                    tg_xxxzz_xyyyyzz_0[j] = pb_x * tg_xxzz_xyyyyzz_0[j] + wp_x[j] * tg_xxzz_xyyyyzz_1[j] + fl1_fx * tg_xzz_xyyyyzz_0[j] - fl1_fx * fl1_fza * tg_xzz_xyyyyzz_1[j] + 0.5 * fl1_fxn * tg_xxzz_yyyyzz_1[j];

                    tg_xxxzz_xyyyzzz_0[j] = pb_x * tg_xxzz_xyyyzzz_0[j] + wp_x[j] * tg_xxzz_xyyyzzz_1[j] + fl1_fx * tg_xzz_xyyyzzz_0[j] - fl1_fx * fl1_fza * tg_xzz_xyyyzzz_1[j] + 0.5 * fl1_fxn * tg_xxzz_yyyzzz_1[j];

                    tg_xxxzz_xyyzzzz_0[j] = pb_x * tg_xxzz_xyyzzzz_0[j] + wp_x[j] * tg_xxzz_xyyzzzz_1[j] + fl1_fx * tg_xzz_xyyzzzz_0[j] - fl1_fx * fl1_fza * tg_xzz_xyyzzzz_1[j] + 0.5 * fl1_fxn * tg_xxzz_yyzzzz_1[j];

                    tg_xxxzz_xyzzzzz_0[j] = pb_x * tg_xxzz_xyzzzzz_0[j] + wp_x[j] * tg_xxzz_xyzzzzz_1[j] + fl1_fx * tg_xzz_xyzzzzz_0[j] - fl1_fx * fl1_fza * tg_xzz_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_xxzz_yzzzzz_1[j];

                    tg_xxxzz_xzzzzzz_0[j] = pb_x * tg_xxzz_xzzzzzz_0[j] + wp_x[j] * tg_xxzz_xzzzzzz_1[j] + fl1_fx * tg_xzz_xzzzzzz_0[j] - fl1_fx * fl1_fza * tg_xzz_xzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xxzz_zzzzzz_1[j];

                    tg_xxxzz_yyyyyyy_0[j] = pb_x * tg_xxzz_yyyyyyy_0[j] + wp_x[j] * tg_xxzz_yyyyyyy_1[j] + fl1_fx * tg_xzz_yyyyyyy_0[j] - fl1_fx * fl1_fza * tg_xzz_yyyyyyy_1[j];

                    tg_xxxzz_yyyyyyz_0[j] = pb_x * tg_xxzz_yyyyyyz_0[j] + wp_x[j] * tg_xxzz_yyyyyyz_1[j] + fl1_fx * tg_xzz_yyyyyyz_0[j] - fl1_fx * fl1_fza * tg_xzz_yyyyyyz_1[j];

                    tg_xxxzz_yyyyyzz_0[j] = pb_x * tg_xxzz_yyyyyzz_0[j] + wp_x[j] * tg_xxzz_yyyyyzz_1[j] + fl1_fx * tg_xzz_yyyyyzz_0[j] - fl1_fx * fl1_fza * tg_xzz_yyyyyzz_1[j];

                    tg_xxxzz_yyyyzzz_0[j] = pb_x * tg_xxzz_yyyyzzz_0[j] + wp_x[j] * tg_xxzz_yyyyzzz_1[j] + fl1_fx * tg_xzz_yyyyzzz_0[j] - fl1_fx * fl1_fza * tg_xzz_yyyyzzz_1[j];

                    tg_xxxzz_yyyzzzz_0[j] = pb_x * tg_xxzz_yyyzzzz_0[j] + wp_x[j] * tg_xxzz_yyyzzzz_1[j] + fl1_fx * tg_xzz_yyyzzzz_0[j] - fl1_fx * fl1_fza * tg_xzz_yyyzzzz_1[j];

                    tg_xxxzz_yyzzzzz_0[j] = pb_x * tg_xxzz_yyzzzzz_0[j] + wp_x[j] * tg_xxzz_yyzzzzz_1[j] + fl1_fx * tg_xzz_yyzzzzz_0[j] - fl1_fx * fl1_fza * tg_xzz_yyzzzzz_1[j];

                    tg_xxxzz_yzzzzzz_0[j] = pb_x * tg_xxzz_yzzzzzz_0[j] + wp_x[j] * tg_xxzz_yzzzzzz_1[j] + fl1_fx * tg_xzz_yzzzzzz_0[j] - fl1_fx * fl1_fza * tg_xzz_yzzzzzz_1[j];

                    tg_xxxzz_zzzzzzz_0[j] = pb_x * tg_xxzz_zzzzzzz_0[j] + wp_x[j] * tg_xxzz_zzzzzzz_1[j] + fl1_fx * tg_xzz_zzzzzzz_0[j] - fl1_fx * fl1_fza * tg_xzz_zzzzzzz_1[j];

                    tg_xxyyy_xxxxxxx_0[j] = pb_x * tg_xyyy_xxxxxxx_0[j] + wp_x[j] * tg_xyyy_xxxxxxx_1[j] + 0.5 * fl1_fx * tg_yyy_xxxxxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxxxxxx_1[j] + 3.5 * fl1_fxn * tg_xyyy_xxxxxx_1[j];

                    tg_xxyyy_xxxxxxy_0[j] = pb_x * tg_xyyy_xxxxxxy_0[j] + wp_x[j] * tg_xyyy_xxxxxxy_1[j] + 0.5 * fl1_fx * tg_yyy_xxxxxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxxxxxy_1[j] + 3.0 * fl1_fxn * tg_xyyy_xxxxxy_1[j];

                    tg_xxyyy_xxxxxxz_0[j] = pb_x * tg_xyyy_xxxxxxz_0[j] + wp_x[j] * tg_xyyy_xxxxxxz_1[j] + 0.5 * fl1_fx * tg_yyy_xxxxxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxxxxxz_1[j] + 3.0 * fl1_fxn * tg_xyyy_xxxxxz_1[j];

                    tg_xxyyy_xxxxxyy_0[j] = pb_x * tg_xyyy_xxxxxyy_0[j] + wp_x[j] * tg_xyyy_xxxxxyy_1[j] + 0.5 * fl1_fx * tg_yyy_xxxxxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxxxxyy_1[j] + 2.5 * fl1_fxn * tg_xyyy_xxxxyy_1[j];

                    tg_xxyyy_xxxxxyz_0[j] = pb_x * tg_xyyy_xxxxxyz_0[j] + wp_x[j] * tg_xyyy_xxxxxyz_1[j] + 0.5 * fl1_fx * tg_yyy_xxxxxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxxxxyz_1[j] + 2.5 * fl1_fxn * tg_xyyy_xxxxyz_1[j];

                    tg_xxyyy_xxxxxzz_0[j] = pb_x * tg_xyyy_xxxxxzz_0[j] + wp_x[j] * tg_xyyy_xxxxxzz_1[j] + 0.5 * fl1_fx * tg_yyy_xxxxxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxxxxzz_1[j] + 2.5 * fl1_fxn * tg_xyyy_xxxxzz_1[j];

                    tg_xxyyy_xxxxyyy_0[j] = pb_x * tg_xyyy_xxxxyyy_0[j] + wp_x[j] * tg_xyyy_xxxxyyy_1[j] + 0.5 * fl1_fx * tg_yyy_xxxxyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxxxyyy_1[j] + 2.0 * fl1_fxn * tg_xyyy_xxxyyy_1[j];

                    tg_xxyyy_xxxxyyz_0[j] = pb_x * tg_xyyy_xxxxyyz_0[j] + wp_x[j] * tg_xyyy_xxxxyyz_1[j] + 0.5 * fl1_fx * tg_yyy_xxxxyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxxxyyz_1[j] + 2.0 * fl1_fxn * tg_xyyy_xxxyyz_1[j];

                    tg_xxyyy_xxxxyzz_0[j] = pb_x * tg_xyyy_xxxxyzz_0[j] + wp_x[j] * tg_xyyy_xxxxyzz_1[j] + 0.5 * fl1_fx * tg_yyy_xxxxyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxxxyzz_1[j] + 2.0 * fl1_fxn * tg_xyyy_xxxyzz_1[j];

                    tg_xxyyy_xxxxzzz_0[j] = pb_x * tg_xyyy_xxxxzzz_0[j] + wp_x[j] * tg_xyyy_xxxxzzz_1[j] + 0.5 * fl1_fx * tg_yyy_xxxxzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxxxzzz_1[j] + 2.0 * fl1_fxn * tg_xyyy_xxxzzz_1[j];

                    tg_xxyyy_xxxyyyy_0[j] = pb_x * tg_xyyy_xxxyyyy_0[j] + wp_x[j] * tg_xyyy_xxxyyyy_1[j] + 0.5 * fl1_fx * tg_yyy_xxxyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxxyyyy_1[j] + 1.5 * fl1_fxn * tg_xyyy_xxyyyy_1[j];

                    tg_xxyyy_xxxyyyz_0[j] = pb_x * tg_xyyy_xxxyyyz_0[j] + wp_x[j] * tg_xyyy_xxxyyyz_1[j] + 0.5 * fl1_fx * tg_yyy_xxxyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_xyyy_xxyyyz_1[j];

                    tg_xxyyy_xxxyyzz_0[j] = pb_x * tg_xyyy_xxxyyzz_0[j] + wp_x[j] * tg_xyyy_xxxyyzz_1[j] + 0.5 * fl1_fx * tg_yyy_xxxyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxxyyzz_1[j] + 1.5 * fl1_fxn * tg_xyyy_xxyyzz_1[j];

                    tg_xxyyy_xxxyzzz_0[j] = pb_x * tg_xyyy_xxxyzzz_0[j] + wp_x[j] * tg_xyyy_xxxyzzz_1[j] + 0.5 * fl1_fx * tg_yyy_xxxyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxxyzzz_1[j] + 1.5 * fl1_fxn * tg_xyyy_xxyzzz_1[j];

                    tg_xxyyy_xxxzzzz_0[j] = pb_x * tg_xyyy_xxxzzzz_0[j] + wp_x[j] * tg_xyyy_xxxzzzz_1[j] + 0.5 * fl1_fx * tg_yyy_xxxzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxxzzzz_1[j] + 1.5 * fl1_fxn * tg_xyyy_xxzzzz_1[j];

                    tg_xxyyy_xxyyyyy_0[j] = pb_x * tg_xyyy_xxyyyyy_0[j] + wp_x[j] * tg_xyyy_xxyyyyy_1[j] + 0.5 * fl1_fx * tg_yyy_xxyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxyyyyy_1[j] + fl1_fxn * tg_xyyy_xyyyyy_1[j];

                    tg_xxyyy_xxyyyyz_0[j] = pb_x * tg_xyyy_xxyyyyz_0[j] + wp_x[j] * tg_xyyy_xxyyyyz_1[j] + 0.5 * fl1_fx * tg_yyy_xxyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxyyyyz_1[j] + fl1_fxn * tg_xyyy_xyyyyz_1[j];

                    tg_xxyyy_xxyyyzz_0[j] = pb_x * tg_xyyy_xxyyyzz_0[j] + wp_x[j] * tg_xyyy_xxyyyzz_1[j] + 0.5 * fl1_fx * tg_yyy_xxyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxyyyzz_1[j] + fl1_fxn * tg_xyyy_xyyyzz_1[j];

                    tg_xxyyy_xxyyzzz_0[j] = pb_x * tg_xyyy_xxyyzzz_0[j] + wp_x[j] * tg_xyyy_xxyyzzz_1[j] + 0.5 * fl1_fx * tg_yyy_xxyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxyyzzz_1[j] + fl1_fxn * tg_xyyy_xyyzzz_1[j];

                    tg_xxyyy_xxyzzzz_0[j] = pb_x * tg_xyyy_xxyzzzz_0[j] + wp_x[j] * tg_xyyy_xxyzzzz_1[j] + 0.5 * fl1_fx * tg_yyy_xxyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxyzzzz_1[j] + fl1_fxn * tg_xyyy_xyzzzz_1[j];

                    tg_xxyyy_xxzzzzz_0[j] = pb_x * tg_xyyy_xxzzzzz_0[j] + wp_x[j] * tg_xyyy_xxzzzzz_1[j] + 0.5 * fl1_fx * tg_yyy_xxzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxzzzzz_1[j] + fl1_fxn * tg_xyyy_xzzzzz_1[j];

                    tg_xxyyy_xyyyyyy_0[j] = pb_x * tg_xyyy_xyyyyyy_0[j] + wp_x[j] * tg_xyyy_xyyyyyy_1[j] + 0.5 * fl1_fx * tg_yyy_xyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xyyyyyy_1[j] + 0.5 * fl1_fxn * tg_xyyy_yyyyyy_1[j];

                    tg_xxyyy_xyyyyyz_0[j] = pb_x * tg_xyyy_xyyyyyz_0[j] + wp_x[j] * tg_xyyy_xyyyyyz_1[j] + 0.5 * fl1_fx * tg_yyy_xyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xyyyyyz_1[j] + 0.5 * fl1_fxn * tg_xyyy_yyyyyz_1[j];

                    tg_xxyyy_xyyyyzz_0[j] = pb_x * tg_xyyy_xyyyyzz_0[j] + wp_x[j] * tg_xyyy_xyyyyzz_1[j] + 0.5 * fl1_fx * tg_yyy_xyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xyyyyzz_1[j] + 0.5 * fl1_fxn * tg_xyyy_yyyyzz_1[j];

                    tg_xxyyy_xyyyzzz_0[j] = pb_x * tg_xyyy_xyyyzzz_0[j] + wp_x[j] * tg_xyyy_xyyyzzz_1[j] + 0.5 * fl1_fx * tg_yyy_xyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xyyyzzz_1[j] + 0.5 * fl1_fxn * tg_xyyy_yyyzzz_1[j];

                    tg_xxyyy_xyyzzzz_0[j] = pb_x * tg_xyyy_xyyzzzz_0[j] + wp_x[j] * tg_xyyy_xyyzzzz_1[j] + 0.5 * fl1_fx * tg_yyy_xyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xyyzzzz_1[j] + 0.5 * fl1_fxn * tg_xyyy_yyzzzz_1[j];

                    tg_xxyyy_xyzzzzz_0[j] = pb_x * tg_xyyy_xyzzzzz_0[j] + wp_x[j] * tg_xyyy_xyzzzzz_1[j] + 0.5 * fl1_fx * tg_yyy_xyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_xyyy_yzzzzz_1[j];

                    tg_xxyyy_xzzzzzz_0[j] = pb_x * tg_xyyy_xzzzzzz_0[j] + wp_x[j] * tg_xyyy_xzzzzzz_1[j] + 0.5 * fl1_fx * tg_yyy_xzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xyyy_zzzzzz_1[j];

                    tg_xxyyy_yyyyyyy_0[j] = pb_x * tg_xyyy_yyyyyyy_0[j] + wp_x[j] * tg_xyyy_yyyyyyy_1[j] + 0.5 * fl1_fx * tg_yyy_yyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_yyyyyyy_1[j];

                    tg_xxyyy_yyyyyyz_0[j] = pb_x * tg_xyyy_yyyyyyz_0[j] + wp_x[j] * tg_xyyy_yyyyyyz_1[j] + 0.5 * fl1_fx * tg_yyy_yyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_yyyyyyz_1[j];

                    tg_xxyyy_yyyyyzz_0[j] = pb_x * tg_xyyy_yyyyyzz_0[j] + wp_x[j] * tg_xyyy_yyyyyzz_1[j] + 0.5 * fl1_fx * tg_yyy_yyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_yyyyyzz_1[j];

                    tg_xxyyy_yyyyzzz_0[j] = pb_x * tg_xyyy_yyyyzzz_0[j] + wp_x[j] * tg_xyyy_yyyyzzz_1[j] + 0.5 * fl1_fx * tg_yyy_yyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_yyyyzzz_1[j];

                    tg_xxyyy_yyyzzzz_0[j] = pb_x * tg_xyyy_yyyzzzz_0[j] + wp_x[j] * tg_xyyy_yyyzzzz_1[j] + 0.5 * fl1_fx * tg_yyy_yyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_yyyzzzz_1[j];

                    tg_xxyyy_yyzzzzz_0[j] = pb_x * tg_xyyy_yyzzzzz_0[j] + wp_x[j] * tg_xyyy_yyzzzzz_1[j] + 0.5 * fl1_fx * tg_yyy_yyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_yyzzzzz_1[j];

                    tg_xxyyy_yzzzzzz_0[j] = pb_x * tg_xyyy_yzzzzzz_0[j] + wp_x[j] * tg_xyyy_yzzzzzz_1[j] + 0.5 * fl1_fx * tg_yyy_yzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_yzzzzzz_1[j];

                    tg_xxyyy_zzzzzzz_0[j] = pb_x * tg_xyyy_zzzzzzz_0[j] + wp_x[j] * tg_xyyy_zzzzzzz_1[j] + 0.5 * fl1_fx * tg_yyy_zzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_zzzzzzz_1[j];

                    tg_xxyyz_xxxxxxx_0[j] = pb_x * tg_xyyz_xxxxxxx_0[j] + wp_x[j] * tg_xyyz_xxxxxxx_1[j] + 0.5 * fl1_fx * tg_yyz_xxxxxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxxxxxx_1[j] + 3.5 * fl1_fxn * tg_xyyz_xxxxxx_1[j];

                    tg_xxyyz_xxxxxxy_0[j] = pb_x * tg_xyyz_xxxxxxy_0[j] + wp_x[j] * tg_xyyz_xxxxxxy_1[j] + 0.5 * fl1_fx * tg_yyz_xxxxxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxxxxxy_1[j] + 3.0 * fl1_fxn * tg_xyyz_xxxxxy_1[j];

                    tg_xxyyz_xxxxxxz_0[j] = pb_x * tg_xyyz_xxxxxxz_0[j] + wp_x[j] * tg_xyyz_xxxxxxz_1[j] + 0.5 * fl1_fx * tg_yyz_xxxxxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxxxxxz_1[j] + 3.0 * fl1_fxn * tg_xyyz_xxxxxz_1[j];

                    tg_xxyyz_xxxxxyy_0[j] = pb_x * tg_xyyz_xxxxxyy_0[j] + wp_x[j] * tg_xyyz_xxxxxyy_1[j] + 0.5 * fl1_fx * tg_yyz_xxxxxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxxxxyy_1[j] + 2.5 * fl1_fxn * tg_xyyz_xxxxyy_1[j];

                    tg_xxyyz_xxxxxyz_0[j] = pb_x * tg_xyyz_xxxxxyz_0[j] + wp_x[j] * tg_xyyz_xxxxxyz_1[j] + 0.5 * fl1_fx * tg_yyz_xxxxxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxxxxyz_1[j] + 2.5 * fl1_fxn * tg_xyyz_xxxxyz_1[j];

                    tg_xxyyz_xxxxxzz_0[j] = pb_x * tg_xyyz_xxxxxzz_0[j] + wp_x[j] * tg_xyyz_xxxxxzz_1[j] + 0.5 * fl1_fx * tg_yyz_xxxxxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxxxxzz_1[j] + 2.5 * fl1_fxn * tg_xyyz_xxxxzz_1[j];

                    tg_xxyyz_xxxxyyy_0[j] = pb_x * tg_xyyz_xxxxyyy_0[j] + wp_x[j] * tg_xyyz_xxxxyyy_1[j] + 0.5 * fl1_fx * tg_yyz_xxxxyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxxxyyy_1[j] + 2.0 * fl1_fxn * tg_xyyz_xxxyyy_1[j];

                    tg_xxyyz_xxxxyyz_0[j] = pb_x * tg_xyyz_xxxxyyz_0[j] + wp_x[j] * tg_xyyz_xxxxyyz_1[j] + 0.5 * fl1_fx * tg_yyz_xxxxyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxxxyyz_1[j] + 2.0 * fl1_fxn * tg_xyyz_xxxyyz_1[j];

                    tg_xxyyz_xxxxyzz_0[j] = pb_x * tg_xyyz_xxxxyzz_0[j] + wp_x[j] * tg_xyyz_xxxxyzz_1[j] + 0.5 * fl1_fx * tg_yyz_xxxxyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxxxyzz_1[j] + 2.0 * fl1_fxn * tg_xyyz_xxxyzz_1[j];

                    tg_xxyyz_xxxxzzz_0[j] = pb_x * tg_xyyz_xxxxzzz_0[j] + wp_x[j] * tg_xyyz_xxxxzzz_1[j] + 0.5 * fl1_fx * tg_yyz_xxxxzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxxxzzz_1[j] + 2.0 * fl1_fxn * tg_xyyz_xxxzzz_1[j];

                    tg_xxyyz_xxxyyyy_0[j] = pb_x * tg_xyyz_xxxyyyy_0[j] + wp_x[j] * tg_xyyz_xxxyyyy_1[j] + 0.5 * fl1_fx * tg_yyz_xxxyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxxyyyy_1[j] + 1.5 * fl1_fxn * tg_xyyz_xxyyyy_1[j];

                    tg_xxyyz_xxxyyyz_0[j] = pb_x * tg_xyyz_xxxyyyz_0[j] + wp_x[j] * tg_xyyz_xxxyyyz_1[j] + 0.5 * fl1_fx * tg_yyz_xxxyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_xyyz_xxyyyz_1[j];

                    tg_xxyyz_xxxyyzz_0[j] = pb_x * tg_xyyz_xxxyyzz_0[j] + wp_x[j] * tg_xyyz_xxxyyzz_1[j] + 0.5 * fl1_fx * tg_yyz_xxxyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxxyyzz_1[j] + 1.5 * fl1_fxn * tg_xyyz_xxyyzz_1[j];

                    tg_xxyyz_xxxyzzz_0[j] = pb_x * tg_xyyz_xxxyzzz_0[j] + wp_x[j] * tg_xyyz_xxxyzzz_1[j] + 0.5 * fl1_fx * tg_yyz_xxxyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxxyzzz_1[j] + 1.5 * fl1_fxn * tg_xyyz_xxyzzz_1[j];

                    tg_xxyyz_xxxzzzz_0[j] = pb_x * tg_xyyz_xxxzzzz_0[j] + wp_x[j] * tg_xyyz_xxxzzzz_1[j] + 0.5 * fl1_fx * tg_yyz_xxxzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxxzzzz_1[j] + 1.5 * fl1_fxn * tg_xyyz_xxzzzz_1[j];

                    tg_xxyyz_xxyyyyy_0[j] = pb_x * tg_xyyz_xxyyyyy_0[j] + wp_x[j] * tg_xyyz_xxyyyyy_1[j] + 0.5 * fl1_fx * tg_yyz_xxyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxyyyyy_1[j] + fl1_fxn * tg_xyyz_xyyyyy_1[j];

                    tg_xxyyz_xxyyyyz_0[j] = pb_x * tg_xyyz_xxyyyyz_0[j] + wp_x[j] * tg_xyyz_xxyyyyz_1[j] + 0.5 * fl1_fx * tg_yyz_xxyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxyyyyz_1[j] + fl1_fxn * tg_xyyz_xyyyyz_1[j];

                    tg_xxyyz_xxyyyzz_0[j] = pb_x * tg_xyyz_xxyyyzz_0[j] + wp_x[j] * tg_xyyz_xxyyyzz_1[j] + 0.5 * fl1_fx * tg_yyz_xxyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxyyyzz_1[j] + fl1_fxn * tg_xyyz_xyyyzz_1[j];

                    tg_xxyyz_xxyyzzz_0[j] = pb_x * tg_xyyz_xxyyzzz_0[j] + wp_x[j] * tg_xyyz_xxyyzzz_1[j] + 0.5 * fl1_fx * tg_yyz_xxyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxyyzzz_1[j] + fl1_fxn * tg_xyyz_xyyzzz_1[j];

                    tg_xxyyz_xxyzzzz_0[j] = pb_x * tg_xyyz_xxyzzzz_0[j] + wp_x[j] * tg_xyyz_xxyzzzz_1[j] + 0.5 * fl1_fx * tg_yyz_xxyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxyzzzz_1[j] + fl1_fxn * tg_xyyz_xyzzzz_1[j];

                    tg_xxyyz_xxzzzzz_0[j] = pb_x * tg_xyyz_xxzzzzz_0[j] + wp_x[j] * tg_xyyz_xxzzzzz_1[j] + 0.5 * fl1_fx * tg_yyz_xxzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxzzzzz_1[j] + fl1_fxn * tg_xyyz_xzzzzz_1[j];

                    tg_xxyyz_xyyyyyy_0[j] = pb_x * tg_xyyz_xyyyyyy_0[j] + wp_x[j] * tg_xyyz_xyyyyyy_1[j] + 0.5 * fl1_fx * tg_yyz_xyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xyyyyyy_1[j] + 0.5 * fl1_fxn * tg_xyyz_yyyyyy_1[j];

                    tg_xxyyz_xyyyyyz_0[j] = pb_x * tg_xyyz_xyyyyyz_0[j] + wp_x[j] * tg_xyyz_xyyyyyz_1[j] + 0.5 * fl1_fx * tg_yyz_xyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xyyyyyz_1[j] + 0.5 * fl1_fxn * tg_xyyz_yyyyyz_1[j];

                    tg_xxyyz_xyyyyzz_0[j] = pb_x * tg_xyyz_xyyyyzz_0[j] + wp_x[j] * tg_xyyz_xyyyyzz_1[j] + 0.5 * fl1_fx * tg_yyz_xyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xyyyyzz_1[j] + 0.5 * fl1_fxn * tg_xyyz_yyyyzz_1[j];

                    tg_xxyyz_xyyyzzz_0[j] = pb_x * tg_xyyz_xyyyzzz_0[j] + wp_x[j] * tg_xyyz_xyyyzzz_1[j] + 0.5 * fl1_fx * tg_yyz_xyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xyyyzzz_1[j] + 0.5 * fl1_fxn * tg_xyyz_yyyzzz_1[j];

                    tg_xxyyz_xyyzzzz_0[j] = pb_x * tg_xyyz_xyyzzzz_0[j] + wp_x[j] * tg_xyyz_xyyzzzz_1[j] + 0.5 * fl1_fx * tg_yyz_xyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xyyzzzz_1[j] + 0.5 * fl1_fxn * tg_xyyz_yyzzzz_1[j];

                    tg_xxyyz_xyzzzzz_0[j] = pb_x * tg_xyyz_xyzzzzz_0[j] + wp_x[j] * tg_xyyz_xyzzzzz_1[j] + 0.5 * fl1_fx * tg_yyz_xyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_xyyz_yzzzzz_1[j];

                    tg_xxyyz_xzzzzzz_0[j] = pb_x * tg_xyyz_xzzzzzz_0[j] + wp_x[j] * tg_xyyz_xzzzzzz_1[j] + 0.5 * fl1_fx * tg_yyz_xzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xyyz_zzzzzz_1[j];

                    tg_xxyyz_yyyyyyy_0[j] = pb_x * tg_xyyz_yyyyyyy_0[j] + wp_x[j] * tg_xyyz_yyyyyyy_1[j] + 0.5 * fl1_fx * tg_yyz_yyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_yyyyyyy_1[j];

                    tg_xxyyz_yyyyyyz_0[j] = pb_x * tg_xyyz_yyyyyyz_0[j] + wp_x[j] * tg_xyyz_yyyyyyz_1[j] + 0.5 * fl1_fx * tg_yyz_yyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_yyyyyyz_1[j];

                    tg_xxyyz_yyyyyzz_0[j] = pb_x * tg_xyyz_yyyyyzz_0[j] + wp_x[j] * tg_xyyz_yyyyyzz_1[j] + 0.5 * fl1_fx * tg_yyz_yyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_yyyyyzz_1[j];

                    tg_xxyyz_yyyyzzz_0[j] = pb_x * tg_xyyz_yyyyzzz_0[j] + wp_x[j] * tg_xyyz_yyyyzzz_1[j] + 0.5 * fl1_fx * tg_yyz_yyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_yyyyzzz_1[j];

                    tg_xxyyz_yyyzzzz_0[j] = pb_x * tg_xyyz_yyyzzzz_0[j] + wp_x[j] * tg_xyyz_yyyzzzz_1[j] + 0.5 * fl1_fx * tg_yyz_yyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_yyyzzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSK_285_380(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (285,380)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {5, -1, -1, -1},
                                             {7, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_7_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_4_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_3_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_4_6_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {6, -1, -1, -1}, 
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

                auto tg_xyyz_yyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 285); 

                auto tg_xyyz_yzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 286); 

                auto tg_xyyz_zzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 287); 

                auto tg_xyzz_xxxxxxx_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 288); 

                auto tg_xyzz_xxxxxxy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 289); 

                auto tg_xyzz_xxxxxxz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 290); 

                auto tg_xyzz_xxxxxyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 291); 

                auto tg_xyzz_xxxxxyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 292); 

                auto tg_xyzz_xxxxxzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 293); 

                auto tg_xyzz_xxxxyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 294); 

                auto tg_xyzz_xxxxyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 295); 

                auto tg_xyzz_xxxxyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 296); 

                auto tg_xyzz_xxxxzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 297); 

                auto tg_xyzz_xxxyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 298); 

                auto tg_xyzz_xxxyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 299); 

                auto tg_xyzz_xxxyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 300); 

                auto tg_xyzz_xxxyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 301); 

                auto tg_xyzz_xxxzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 302); 

                auto tg_xyzz_xxyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 303); 

                auto tg_xyzz_xxyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 304); 

                auto tg_xyzz_xxyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 305); 

                auto tg_xyzz_xxyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 306); 

                auto tg_xyzz_xxyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 307); 

                auto tg_xyzz_xxzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 308); 

                auto tg_xyzz_xyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 309); 

                auto tg_xyzz_xyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 310); 

                auto tg_xyzz_xyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 311); 

                auto tg_xyzz_xyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 312); 

                auto tg_xyzz_xyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 313); 

                auto tg_xyzz_xyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 314); 

                auto tg_xyzz_xzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 315); 

                auto tg_xyzz_yyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 316); 

                auto tg_xyzz_yyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 317); 

                auto tg_xyzz_yyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 318); 

                auto tg_xyzz_yyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 319); 

                auto tg_xyzz_yyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 320); 

                auto tg_xyzz_yyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 321); 

                auto tg_xyzz_yzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 322); 

                auto tg_xyzz_zzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 323); 

                auto tg_xzzz_xxxxxxx_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 324); 

                auto tg_xzzz_xxxxxxy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 325); 

                auto tg_xzzz_xxxxxxz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 326); 

                auto tg_xzzz_xxxxxyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 327); 

                auto tg_xzzz_xxxxxyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 328); 

                auto tg_xzzz_xxxxxzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 329); 

                auto tg_xzzz_xxxxyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 330); 

                auto tg_xzzz_xxxxyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 331); 

                auto tg_xzzz_xxxxyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 332); 

                auto tg_xzzz_xxxxzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 333); 

                auto tg_xzzz_xxxyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 334); 

                auto tg_xzzz_xxxyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 335); 

                auto tg_xzzz_xxxyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 336); 

                auto tg_xzzz_xxxyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 337); 

                auto tg_xzzz_xxxzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 338); 

                auto tg_xzzz_xxyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 339); 

                auto tg_xzzz_xxyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 340); 

                auto tg_xzzz_xxyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 341); 

                auto tg_xzzz_xxyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 342); 

                auto tg_xzzz_xxyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 343); 

                auto tg_xzzz_xxzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 344); 

                auto tg_xzzz_xyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 345); 

                auto tg_xzzz_xyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 346); 

                auto tg_xzzz_xyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 347); 

                auto tg_xzzz_xyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 348); 

                auto tg_xzzz_xyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 349); 

                auto tg_xzzz_xyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 350); 

                auto tg_xzzz_xzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 351); 

                auto tg_xzzz_yyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 352); 

                auto tg_xzzz_yyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 353); 

                auto tg_xzzz_yyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 354); 

                auto tg_xzzz_yyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 355); 

                auto tg_xzzz_yyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 356); 

                auto tg_xzzz_yyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 357); 

                auto tg_xzzz_yzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 358); 

                auto tg_xzzz_zzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 359); 

                auto tg_yyyy_xxxxxxx_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 360); 

                auto tg_yyyy_xxxxxxy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 361); 

                auto tg_yyyy_xxxxxxz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 362); 

                auto tg_yyyy_xxxxxyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 363); 

                auto tg_yyyy_xxxxxyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 364); 

                auto tg_yyyy_xxxxxzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 365); 

                auto tg_yyyy_xxxxyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 366); 

                auto tg_yyyy_xxxxyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 367); 

                auto tg_yyyy_xxxxyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 368); 

                auto tg_yyyy_xxxxzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 369); 

                auto tg_yyyy_xxxyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 370); 

                auto tg_yyyy_xxxyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 371); 

                auto tg_yyyy_xxxyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 372); 

                auto tg_yyyy_xxxyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 373); 

                auto tg_yyyy_xxxzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 374); 

                auto tg_yyyy_xxyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 375); 

                auto tg_yyyy_xxyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 376); 

                auto tg_yyyy_xxyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 377); 

                auto tg_yyyy_xxyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 378); 

                auto tg_yyyy_xxyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 379); 

                auto tg_xyyz_yyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 285); 

                auto tg_xyyz_yzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 286); 

                auto tg_xyyz_zzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 287); 

                auto tg_xyzz_xxxxxxx_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 288); 

                auto tg_xyzz_xxxxxxy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 289); 

                auto tg_xyzz_xxxxxxz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 290); 

                auto tg_xyzz_xxxxxyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 291); 

                auto tg_xyzz_xxxxxyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 292); 

                auto tg_xyzz_xxxxxzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 293); 

                auto tg_xyzz_xxxxyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 294); 

                auto tg_xyzz_xxxxyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 295); 

                auto tg_xyzz_xxxxyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 296); 

                auto tg_xyzz_xxxxzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 297); 

                auto tg_xyzz_xxxyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 298); 

                auto tg_xyzz_xxxyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 299); 

                auto tg_xyzz_xxxyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 300); 

                auto tg_xyzz_xxxyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 301); 

                auto tg_xyzz_xxxzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 302); 

                auto tg_xyzz_xxyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 303); 

                auto tg_xyzz_xxyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 304); 

                auto tg_xyzz_xxyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 305); 

                auto tg_xyzz_xxyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 306); 

                auto tg_xyzz_xxyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 307); 

                auto tg_xyzz_xxzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 308); 

                auto tg_xyzz_xyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 309); 

                auto tg_xyzz_xyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 310); 

                auto tg_xyzz_xyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 311); 

                auto tg_xyzz_xyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 312); 

                auto tg_xyzz_xyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 313); 

                auto tg_xyzz_xyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 314); 

                auto tg_xyzz_xzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 315); 

                auto tg_xyzz_yyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 316); 

                auto tg_xyzz_yyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 317); 

                auto tg_xyzz_yyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 318); 

                auto tg_xyzz_yyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 319); 

                auto tg_xyzz_yyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 320); 

                auto tg_xyzz_yyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 321); 

                auto tg_xyzz_yzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 322); 

                auto tg_xyzz_zzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 323); 

                auto tg_xzzz_xxxxxxx_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 324); 

                auto tg_xzzz_xxxxxxy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 325); 

                auto tg_xzzz_xxxxxxz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 326); 

                auto tg_xzzz_xxxxxyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 327); 

                auto tg_xzzz_xxxxxyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 328); 

                auto tg_xzzz_xxxxxzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 329); 

                auto tg_xzzz_xxxxyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 330); 

                auto tg_xzzz_xxxxyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 331); 

                auto tg_xzzz_xxxxyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 332); 

                auto tg_xzzz_xxxxzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 333); 

                auto tg_xzzz_xxxyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 334); 

                auto tg_xzzz_xxxyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 335); 

                auto tg_xzzz_xxxyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 336); 

                auto tg_xzzz_xxxyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 337); 

                auto tg_xzzz_xxxzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 338); 

                auto tg_xzzz_xxyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 339); 

                auto tg_xzzz_xxyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 340); 

                auto tg_xzzz_xxyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 341); 

                auto tg_xzzz_xxyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 342); 

                auto tg_xzzz_xxyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 343); 

                auto tg_xzzz_xxzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 344); 

                auto tg_xzzz_xyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 345); 

                auto tg_xzzz_xyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 346); 

                auto tg_xzzz_xyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 347); 

                auto tg_xzzz_xyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 348); 

                auto tg_xzzz_xyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 349); 

                auto tg_xzzz_xyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 350); 

                auto tg_xzzz_xzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 351); 

                auto tg_xzzz_yyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 352); 

                auto tg_xzzz_yyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 353); 

                auto tg_xzzz_yyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 354); 

                auto tg_xzzz_yyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 355); 

                auto tg_xzzz_yyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 356); 

                auto tg_xzzz_yyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 357); 

                auto tg_xzzz_yzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 358); 

                auto tg_xzzz_zzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 359); 

                auto tg_yyyy_xxxxxxx_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 360); 

                auto tg_yyyy_xxxxxxy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 361); 

                auto tg_yyyy_xxxxxxz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 362); 

                auto tg_yyyy_xxxxxyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 363); 

                auto tg_yyyy_xxxxxyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 364); 

                auto tg_yyyy_xxxxxzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 365); 

                auto tg_yyyy_xxxxyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 366); 

                auto tg_yyyy_xxxxyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 367); 

                auto tg_yyyy_xxxxyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 368); 

                auto tg_yyyy_xxxxzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 369); 

                auto tg_yyyy_xxxyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 370); 

                auto tg_yyyy_xxxyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 371); 

                auto tg_yyyy_xxxyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 372); 

                auto tg_yyyy_xxxyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 373); 

                auto tg_yyyy_xxxzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 374); 

                auto tg_yyyy_xxyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 375); 

                auto tg_yyyy_xxyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 376); 

                auto tg_yyyy_xxyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 377); 

                auto tg_yyyy_xxyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 378); 

                auto tg_yyyy_xxyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 379); 

                auto tg_yyz_yyzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 285); 

                auto tg_yyz_yzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 286); 

                auto tg_yyz_zzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 287); 

                auto tg_yzz_xxxxxxx_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 288); 

                auto tg_yzz_xxxxxxy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 289); 

                auto tg_yzz_xxxxxxz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 290); 

                auto tg_yzz_xxxxxyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 291); 

                auto tg_yzz_xxxxxyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 292); 

                auto tg_yzz_xxxxxzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 293); 

                auto tg_yzz_xxxxyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 294); 

                auto tg_yzz_xxxxyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 295); 

                auto tg_yzz_xxxxyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 296); 

                auto tg_yzz_xxxxzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 297); 

                auto tg_yzz_xxxyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 298); 

                auto tg_yzz_xxxyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 299); 

                auto tg_yzz_xxxyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 300); 

                auto tg_yzz_xxxyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 301); 

                auto tg_yzz_xxxzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 302); 

                auto tg_yzz_xxyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 303); 

                auto tg_yzz_xxyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 304); 

                auto tg_yzz_xxyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 305); 

                auto tg_yzz_xxyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 306); 

                auto tg_yzz_xxyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 307); 

                auto tg_yzz_xxzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 308); 

                auto tg_yzz_xyyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 309); 

                auto tg_yzz_xyyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 310); 

                auto tg_yzz_xyyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 311); 

                auto tg_yzz_xyyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 312); 

                auto tg_yzz_xyyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 313); 

                auto tg_yzz_xyzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 314); 

                auto tg_yzz_xzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 315); 

                auto tg_yzz_yyyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 316); 

                auto tg_yzz_yyyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 317); 

                auto tg_yzz_yyyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 318); 

                auto tg_yzz_yyyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 319); 

                auto tg_yzz_yyyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 320); 

                auto tg_yzz_yyzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 321); 

                auto tg_yzz_yzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 322); 

                auto tg_yzz_zzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 323); 

                auto tg_zzz_xxxxxxx_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 324); 

                auto tg_zzz_xxxxxxy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 325); 

                auto tg_zzz_xxxxxxz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 326); 

                auto tg_zzz_xxxxxyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 327); 

                auto tg_zzz_xxxxxyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 328); 

                auto tg_zzz_xxxxxzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 329); 

                auto tg_zzz_xxxxyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 330); 

                auto tg_zzz_xxxxyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 331); 

                auto tg_zzz_xxxxyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 332); 

                auto tg_zzz_xxxxzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 333); 

                auto tg_zzz_xxxyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 334); 

                auto tg_zzz_xxxyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 335); 

                auto tg_zzz_xxxyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 336); 

                auto tg_zzz_xxxyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 337); 

                auto tg_zzz_xxxzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 338); 

                auto tg_zzz_xxyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 339); 

                auto tg_zzz_xxyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 340); 

                auto tg_zzz_xxyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 341); 

                auto tg_zzz_xxyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 342); 

                auto tg_zzz_xxyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 343); 

                auto tg_zzz_xxzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 344); 

                auto tg_zzz_xyyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 345); 

                auto tg_zzz_xyyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 346); 

                auto tg_zzz_xyyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 347); 

                auto tg_zzz_xyyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 348); 

                auto tg_zzz_xyyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 349); 

                auto tg_zzz_xyzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 350); 

                auto tg_zzz_xzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 351); 

                auto tg_zzz_yyyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 352); 

                auto tg_zzz_yyyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 353); 

                auto tg_zzz_yyyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 354); 

                auto tg_zzz_yyyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 355); 

                auto tg_zzz_yyyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 356); 

                auto tg_zzz_yyzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 357); 

                auto tg_zzz_yzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 358); 

                auto tg_zzz_zzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 359); 

                auto tg_yyz_yyzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 285); 

                auto tg_yyz_yzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 286); 

                auto tg_yyz_zzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 287); 

                auto tg_yzz_xxxxxxx_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 288); 

                auto tg_yzz_xxxxxxy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 289); 

                auto tg_yzz_xxxxxxz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 290); 

                auto tg_yzz_xxxxxyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 291); 

                auto tg_yzz_xxxxxyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 292); 

                auto tg_yzz_xxxxxzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 293); 

                auto tg_yzz_xxxxyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 294); 

                auto tg_yzz_xxxxyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 295); 

                auto tg_yzz_xxxxyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 296); 

                auto tg_yzz_xxxxzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 297); 

                auto tg_yzz_xxxyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 298); 

                auto tg_yzz_xxxyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 299); 

                auto tg_yzz_xxxyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 300); 

                auto tg_yzz_xxxyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 301); 

                auto tg_yzz_xxxzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 302); 

                auto tg_yzz_xxyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 303); 

                auto tg_yzz_xxyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 304); 

                auto tg_yzz_xxyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 305); 

                auto tg_yzz_xxyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 306); 

                auto tg_yzz_xxyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 307); 

                auto tg_yzz_xxzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 308); 

                auto tg_yzz_xyyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 309); 

                auto tg_yzz_xyyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 310); 

                auto tg_yzz_xyyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 311); 

                auto tg_yzz_xyyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 312); 

                auto tg_yzz_xyyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 313); 

                auto tg_yzz_xyzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 314); 

                auto tg_yzz_xzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 315); 

                auto tg_yzz_yyyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 316); 

                auto tg_yzz_yyyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 317); 

                auto tg_yzz_yyyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 318); 

                auto tg_yzz_yyyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 319); 

                auto tg_yzz_yyyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 320); 

                auto tg_yzz_yyzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 321); 

                auto tg_yzz_yzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 322); 

                auto tg_yzz_zzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 323); 

                auto tg_zzz_xxxxxxx_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 324); 

                auto tg_zzz_xxxxxxy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 325); 

                auto tg_zzz_xxxxxxz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 326); 

                auto tg_zzz_xxxxxyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 327); 

                auto tg_zzz_xxxxxyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 328); 

                auto tg_zzz_xxxxxzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 329); 

                auto tg_zzz_xxxxyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 330); 

                auto tg_zzz_xxxxyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 331); 

                auto tg_zzz_xxxxyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 332); 

                auto tg_zzz_xxxxzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 333); 

                auto tg_zzz_xxxyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 334); 

                auto tg_zzz_xxxyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 335); 

                auto tg_zzz_xxxyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 336); 

                auto tg_zzz_xxxyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 337); 

                auto tg_zzz_xxxzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 338); 

                auto tg_zzz_xxyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 339); 

                auto tg_zzz_xxyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 340); 

                auto tg_zzz_xxyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 341); 

                auto tg_zzz_xxyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 342); 

                auto tg_zzz_xxyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 343); 

                auto tg_zzz_xxzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 344); 

                auto tg_zzz_xyyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 345); 

                auto tg_zzz_xyyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 346); 

                auto tg_zzz_xyyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 347); 

                auto tg_zzz_xyyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 348); 

                auto tg_zzz_xyyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 349); 

                auto tg_zzz_xyzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 350); 

                auto tg_zzz_xzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 351); 

                auto tg_zzz_yyyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 352); 

                auto tg_zzz_yyyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 353); 

                auto tg_zzz_yyyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 354); 

                auto tg_zzz_yyyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 355); 

                auto tg_zzz_yyyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 356); 

                auto tg_zzz_yyzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 357); 

                auto tg_zzz_yzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 358); 

                auto tg_zzz_zzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 359); 

                auto tg_xyzz_xxxxxx_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 224); 

                auto tg_xyzz_xxxxxy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 225); 

                auto tg_xyzz_xxxxxz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 226); 

                auto tg_xyzz_xxxxyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 227); 

                auto tg_xyzz_xxxxyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 228); 

                auto tg_xyzz_xxxxzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 229); 

                auto tg_xyzz_xxxyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 230); 

                auto tg_xyzz_xxxyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 231); 

                auto tg_xyzz_xxxyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 232); 

                auto tg_xyzz_xxxzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 233); 

                auto tg_xyzz_xxyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 234); 

                auto tg_xyzz_xxyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 235); 

                auto tg_xyzz_xxyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 236); 

                auto tg_xyzz_xxyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 237); 

                auto tg_xyzz_xxzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 238); 

                auto tg_xyzz_xyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 239); 

                auto tg_xyzz_xyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 240); 

                auto tg_xyzz_xyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 241); 

                auto tg_xyzz_xyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 242); 

                auto tg_xyzz_xyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 243); 

                auto tg_xyzz_xzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 244); 

                auto tg_xyzz_yyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 245); 

                auto tg_xyzz_yyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 246); 

                auto tg_xyzz_yyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 247); 

                auto tg_xyzz_yyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 248); 

                auto tg_xyzz_yyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 249); 

                auto tg_xyzz_yzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 250); 

                auto tg_xyzz_zzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 251); 

                auto tg_xzzz_xxxxxx_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 252); 

                auto tg_xzzz_xxxxxy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 253); 

                auto tg_xzzz_xxxxxz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 254); 

                auto tg_xzzz_xxxxyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 255); 

                auto tg_xzzz_xxxxyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 256); 

                auto tg_xzzz_xxxxzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 257); 

                auto tg_xzzz_xxxyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 258); 

                auto tg_xzzz_xxxyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 259); 

                auto tg_xzzz_xxxyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 260); 

                auto tg_xzzz_xxxzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 261); 

                auto tg_xzzz_xxyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 262); 

                auto tg_xzzz_xxyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 263); 

                auto tg_xzzz_xxyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 264); 

                auto tg_xzzz_xxyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 265); 

                auto tg_xzzz_xxzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 266); 

                auto tg_xzzz_xyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 267); 

                auto tg_xzzz_xyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 268); 

                auto tg_xzzz_xyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 269); 

                auto tg_xzzz_xyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 270); 

                auto tg_xzzz_xyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 271); 

                auto tg_xzzz_xzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 272); 

                auto tg_xzzz_yyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 273); 

                auto tg_xzzz_yyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 274); 

                auto tg_xzzz_yyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 275); 

                auto tg_xzzz_yyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 276); 

                auto tg_xzzz_yyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 277); 

                auto tg_xzzz_yzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 278); 

                auto tg_xzzz_zzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 279); 

                auto tg_yyyy_xxxxxx_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 280); 

                auto tg_yyyy_xxxxxy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 281); 

                auto tg_yyyy_xxxxxz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 282); 

                auto tg_yyyy_xxxxyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 283); 

                auto tg_yyyy_xxxxyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 284); 

                auto tg_yyyy_xxxxzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 285); 

                auto tg_yyyy_xxxyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 286); 

                auto tg_yyyy_xxxyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 287); 

                auto tg_yyyy_xxxyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 288); 

                auto tg_yyyy_xxxzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 289); 

                auto tg_yyyy_xxyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 290); 

                auto tg_yyyy_xxyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 291); 

                auto tg_yyyy_xxyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 292); 

                auto tg_yyyy_xxyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 293); 

                auto tg_yyyy_xxzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 294); 

                auto tg_yyyy_xyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 295); 

                auto tg_yyyy_xyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 296); 

                auto tg_yyyy_xyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 297); 

                auto tg_yyyy_xyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 298); 

                auto tg_yyyy_xyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 299); 

                // set up pointers to integrals

                auto tg_xxyyz_yyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 285); 

                auto tg_xxyyz_yzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 286); 

                auto tg_xxyyz_zzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 287); 

                auto tg_xxyzz_xxxxxxx_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 288); 

                auto tg_xxyzz_xxxxxxy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 289); 

                auto tg_xxyzz_xxxxxxz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 290); 

                auto tg_xxyzz_xxxxxyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 291); 

                auto tg_xxyzz_xxxxxyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 292); 

                auto tg_xxyzz_xxxxxzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 293); 

                auto tg_xxyzz_xxxxyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 294); 

                auto tg_xxyzz_xxxxyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 295); 

                auto tg_xxyzz_xxxxyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 296); 

                auto tg_xxyzz_xxxxzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 297); 

                auto tg_xxyzz_xxxyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 298); 

                auto tg_xxyzz_xxxyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 299); 

                auto tg_xxyzz_xxxyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 300); 

                auto tg_xxyzz_xxxyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 301); 

                auto tg_xxyzz_xxxzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 302); 

                auto tg_xxyzz_xxyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 303); 

                auto tg_xxyzz_xxyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 304); 

                auto tg_xxyzz_xxyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 305); 

                auto tg_xxyzz_xxyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 306); 

                auto tg_xxyzz_xxyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 307); 

                auto tg_xxyzz_xxzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 308); 

                auto tg_xxyzz_xyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 309); 

                auto tg_xxyzz_xyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 310); 

                auto tg_xxyzz_xyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 311); 

                auto tg_xxyzz_xyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 312); 

                auto tg_xxyzz_xyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 313); 

                auto tg_xxyzz_xyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 314); 

                auto tg_xxyzz_xzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 315); 

                auto tg_xxyzz_yyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 316); 

                auto tg_xxyzz_yyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 317); 

                auto tg_xxyzz_yyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 318); 

                auto tg_xxyzz_yyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 319); 

                auto tg_xxyzz_yyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 320); 

                auto tg_xxyzz_yyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 321); 

                auto tg_xxyzz_yzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 322); 

                auto tg_xxyzz_zzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 323); 

                auto tg_xxzzz_xxxxxxx_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 324); 

                auto tg_xxzzz_xxxxxxy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 325); 

                auto tg_xxzzz_xxxxxxz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 326); 

                auto tg_xxzzz_xxxxxyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 327); 

                auto tg_xxzzz_xxxxxyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 328); 

                auto tg_xxzzz_xxxxxzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 329); 

                auto tg_xxzzz_xxxxyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 330); 

                auto tg_xxzzz_xxxxyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 331); 

                auto tg_xxzzz_xxxxyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 332); 

                auto tg_xxzzz_xxxxzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 333); 

                auto tg_xxzzz_xxxyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 334); 

                auto tg_xxzzz_xxxyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 335); 

                auto tg_xxzzz_xxxyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 336); 

                auto tg_xxzzz_xxxyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 337); 

                auto tg_xxzzz_xxxzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 338); 

                auto tg_xxzzz_xxyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 339); 

                auto tg_xxzzz_xxyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 340); 

                auto tg_xxzzz_xxyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 341); 

                auto tg_xxzzz_xxyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 342); 

                auto tg_xxzzz_xxyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 343); 

                auto tg_xxzzz_xxzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 344); 

                auto tg_xxzzz_xyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 345); 

                auto tg_xxzzz_xyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 346); 

                auto tg_xxzzz_xyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 347); 

                auto tg_xxzzz_xyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 348); 

                auto tg_xxzzz_xyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 349); 

                auto tg_xxzzz_xyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 350); 

                auto tg_xxzzz_xzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 351); 

                auto tg_xxzzz_yyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 352); 

                auto tg_xxzzz_yyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 353); 

                auto tg_xxzzz_yyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 354); 

                auto tg_xxzzz_yyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 355); 

                auto tg_xxzzz_yyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 356); 

                auto tg_xxzzz_yyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 357); 

                auto tg_xxzzz_yzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 358); 

                auto tg_xxzzz_zzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 359); 

                auto tg_xyyyy_xxxxxxx_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 360); 

                auto tg_xyyyy_xxxxxxy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 361); 

                auto tg_xyyyy_xxxxxxz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 362); 

                auto tg_xyyyy_xxxxxyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 363); 

                auto tg_xyyyy_xxxxxyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 364); 

                auto tg_xyyyy_xxxxxzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 365); 

                auto tg_xyyyy_xxxxyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 366); 

                auto tg_xyyyy_xxxxyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 367); 

                auto tg_xyyyy_xxxxyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 368); 

                auto tg_xyyyy_xxxxzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 369); 

                auto tg_xyyyy_xxxyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 370); 

                auto tg_xyyyy_xxxyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 371); 

                auto tg_xyyyy_xxxyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 372); 

                auto tg_xyyyy_xxxyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 373); 

                auto tg_xyyyy_xxxzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 374); 

                auto tg_xyyyy_xxyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 375); 

                auto tg_xyyyy_xxyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 376); 

                auto tg_xyyyy_xxyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 377); 

                auto tg_xyyyy_xxyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 378); 

                auto tg_xyyyy_xxyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 379); 

                // Batch of Integrals (285,380)

                #pragma omp simd aligned(fxn, fza, tg_xxyyz_yyzzzzz_0, tg_xxyyz_yzzzzzz_0, tg_xxyyz_zzzzzzz_0, \
                                         tg_xxyzz_xxxxxxx_0, tg_xxyzz_xxxxxxy_0, tg_xxyzz_xxxxxxz_0, tg_xxyzz_xxxxxyy_0, \
                                         tg_xxyzz_xxxxxyz_0, tg_xxyzz_xxxxxzz_0, tg_xxyzz_xxxxyyy_0, tg_xxyzz_xxxxyyz_0, \
                                         tg_xxyzz_xxxxyzz_0, tg_xxyzz_xxxxzzz_0, tg_xxyzz_xxxyyyy_0, tg_xxyzz_xxxyyyz_0, \
                                         tg_xxyzz_xxxyyzz_0, tg_xxyzz_xxxyzzz_0, tg_xxyzz_xxxzzzz_0, tg_xxyzz_xxyyyyy_0, \
                                         tg_xxyzz_xxyyyyz_0, tg_xxyzz_xxyyyzz_0, tg_xxyzz_xxyyzzz_0, tg_xxyzz_xxyzzzz_0, \
                                         tg_xxyzz_xxzzzzz_0, tg_xxyzz_xyyyyyy_0, tg_xxyzz_xyyyyyz_0, tg_xxyzz_xyyyyzz_0, \
                                         tg_xxyzz_xyyyzzz_0, tg_xxyzz_xyyzzzz_0, tg_xxyzz_xyzzzzz_0, tg_xxyzz_xzzzzzz_0, \
                                         tg_xxyzz_yyyyyyy_0, tg_xxyzz_yyyyyyz_0, tg_xxyzz_yyyyyzz_0, tg_xxyzz_yyyyzzz_0, \
                                         tg_xxyzz_yyyzzzz_0, tg_xxyzz_yyzzzzz_0, tg_xxyzz_yzzzzzz_0, tg_xxyzz_zzzzzzz_0, \
                                         tg_xxzzz_xxxxxxx_0, tg_xxzzz_xxxxxxy_0, tg_xxzzz_xxxxxxz_0, tg_xxzzz_xxxxxyy_0, \
                                         tg_xxzzz_xxxxxyz_0, tg_xxzzz_xxxxxzz_0, tg_xxzzz_xxxxyyy_0, tg_xxzzz_xxxxyyz_0, \
                                         tg_xxzzz_xxxxyzz_0, tg_xxzzz_xxxxzzz_0, tg_xxzzz_xxxyyyy_0, tg_xxzzz_xxxyyyz_0, \
                                         tg_xxzzz_xxxyyzz_0, tg_xxzzz_xxxyzzz_0, tg_xxzzz_xxxzzzz_0, tg_xxzzz_xxyyyyy_0, \
                                         tg_xxzzz_xxyyyyz_0, tg_xxzzz_xxyyyzz_0, tg_xxzzz_xxyyzzz_0, tg_xxzzz_xxyzzzz_0, \
                                         tg_xxzzz_xxzzzzz_0, tg_xxzzz_xyyyyyy_0, tg_xxzzz_xyyyyyz_0, tg_xxzzz_xyyyyzz_0, \
                                         tg_xxzzz_xyyyzzz_0, tg_xxzzz_xyyzzzz_0, tg_xxzzz_xyzzzzz_0, tg_xxzzz_xzzzzzz_0, \
                                         tg_xxzzz_yyyyyyy_0, tg_xxzzz_yyyyyyz_0, tg_xxzzz_yyyyyzz_0, tg_xxzzz_yyyyzzz_0, \
                                         tg_xxzzz_yyyzzzz_0, tg_xxzzz_yyzzzzz_0, tg_xxzzz_yzzzzzz_0, tg_xxzzz_zzzzzzz_0, \
                                         tg_xyyyy_xxxxxxx_0, tg_xyyyy_xxxxxxy_0, tg_xyyyy_xxxxxxz_0, tg_xyyyy_xxxxxyy_0, \
                                         tg_xyyyy_xxxxxyz_0, tg_xyyyy_xxxxxzz_0, tg_xyyyy_xxxxyyy_0, tg_xyyyy_xxxxyyz_0, \
                                         tg_xyyyy_xxxxyzz_0, tg_xyyyy_xxxxzzz_0, tg_xyyyy_xxxyyyy_0, tg_xyyyy_xxxyyyz_0, \
                                         tg_xyyyy_xxxyyzz_0, tg_xyyyy_xxxyzzz_0, tg_xyyyy_xxxzzzz_0, tg_xyyyy_xxyyyyy_0, \
                                         tg_xyyyy_xxyyyyz_0, tg_xyyyy_xxyyyzz_0, tg_xyyyy_xxyyzzz_0, tg_xyyyy_xxyzzzz_0, \
                                         tg_xyyz_yyzzzzz_0, tg_xyyz_yyzzzzz_1, tg_xyyz_yzzzzzz_0, tg_xyyz_yzzzzzz_1, \
                                         tg_xyyz_zzzzzzz_0, tg_xyyz_zzzzzzz_1, tg_xyzz_xxxxxx_1, tg_xyzz_xxxxxxx_0, \
                                         tg_xyzz_xxxxxxx_1, tg_xyzz_xxxxxxy_0, tg_xyzz_xxxxxxy_1, tg_xyzz_xxxxxxz_0, \
                                         tg_xyzz_xxxxxxz_1, tg_xyzz_xxxxxy_1, tg_xyzz_xxxxxyy_0, tg_xyzz_xxxxxyy_1, \
                                         tg_xyzz_xxxxxyz_0, tg_xyzz_xxxxxyz_1, tg_xyzz_xxxxxz_1, tg_xyzz_xxxxxzz_0, \
                                         tg_xyzz_xxxxxzz_1, tg_xyzz_xxxxyy_1, tg_xyzz_xxxxyyy_0, tg_xyzz_xxxxyyy_1, \
                                         tg_xyzz_xxxxyyz_0, tg_xyzz_xxxxyyz_1, tg_xyzz_xxxxyz_1, tg_xyzz_xxxxyzz_0, \
                                         tg_xyzz_xxxxyzz_1, tg_xyzz_xxxxzz_1, tg_xyzz_xxxxzzz_0, tg_xyzz_xxxxzzz_1, \
                                         tg_xyzz_xxxyyy_1, tg_xyzz_xxxyyyy_0, tg_xyzz_xxxyyyy_1, tg_xyzz_xxxyyyz_0, \
                                         tg_xyzz_xxxyyyz_1, tg_xyzz_xxxyyz_1, tg_xyzz_xxxyyzz_0, tg_xyzz_xxxyyzz_1, \
                                         tg_xyzz_xxxyzz_1, tg_xyzz_xxxyzzz_0, tg_xyzz_xxxyzzz_1, tg_xyzz_xxxzzz_1, \
                                         tg_xyzz_xxxzzzz_0, tg_xyzz_xxxzzzz_1, tg_xyzz_xxyyyy_1, tg_xyzz_xxyyyyy_0, \
                                         tg_xyzz_xxyyyyy_1, tg_xyzz_xxyyyyz_0, tg_xyzz_xxyyyyz_1, tg_xyzz_xxyyyz_1, \
                                         tg_xyzz_xxyyyzz_0, tg_xyzz_xxyyyzz_1, tg_xyzz_xxyyzz_1, tg_xyzz_xxyyzzz_0, \
                                         tg_xyzz_xxyyzzz_1, tg_xyzz_xxyzzz_1, tg_xyzz_xxyzzzz_0, tg_xyzz_xxyzzzz_1, \
                                         tg_xyzz_xxzzzz_1, tg_xyzz_xxzzzzz_0, tg_xyzz_xxzzzzz_1, tg_xyzz_xyyyyy_1, \
                                         tg_xyzz_xyyyyyy_0, tg_xyzz_xyyyyyy_1, tg_xyzz_xyyyyyz_0, tg_xyzz_xyyyyyz_1, \
                                         tg_xyzz_xyyyyz_1, tg_xyzz_xyyyyzz_0, tg_xyzz_xyyyyzz_1, tg_xyzz_xyyyzz_1, \
                                         tg_xyzz_xyyyzzz_0, tg_xyzz_xyyyzzz_1, tg_xyzz_xyyzzz_1, tg_xyzz_xyyzzzz_0, \
                                         tg_xyzz_xyyzzzz_1, tg_xyzz_xyzzzz_1, tg_xyzz_xyzzzzz_0, tg_xyzz_xyzzzzz_1, \
                                         tg_xyzz_xzzzzz_1, tg_xyzz_xzzzzzz_0, tg_xyzz_xzzzzzz_1, tg_xyzz_yyyyyy_1, \
                                         tg_xyzz_yyyyyyy_0, tg_xyzz_yyyyyyy_1, tg_xyzz_yyyyyyz_0, tg_xyzz_yyyyyyz_1, \
                                         tg_xyzz_yyyyyz_1, tg_xyzz_yyyyyzz_0, tg_xyzz_yyyyyzz_1, tg_xyzz_yyyyzz_1, \
                                         tg_xyzz_yyyyzzz_0, tg_xyzz_yyyyzzz_1, tg_xyzz_yyyzzz_1, tg_xyzz_yyyzzzz_0, \
                                         tg_xyzz_yyyzzzz_1, tg_xyzz_yyzzzz_1, tg_xyzz_yyzzzzz_0, tg_xyzz_yyzzzzz_1, \
                                         tg_xyzz_yzzzzz_1, tg_xyzz_yzzzzzz_0, tg_xyzz_yzzzzzz_1, tg_xyzz_zzzzzz_1, \
                                         tg_xyzz_zzzzzzz_0, tg_xyzz_zzzzzzz_1, tg_xzzz_xxxxxx_1, tg_xzzz_xxxxxxx_0, \
                                         tg_xzzz_xxxxxxx_1, tg_xzzz_xxxxxxy_0, tg_xzzz_xxxxxxy_1, tg_xzzz_xxxxxxz_0, \
                                         tg_xzzz_xxxxxxz_1, tg_xzzz_xxxxxy_1, tg_xzzz_xxxxxyy_0, tg_xzzz_xxxxxyy_1, \
                                         tg_xzzz_xxxxxyz_0, tg_xzzz_xxxxxyz_1, tg_xzzz_xxxxxz_1, tg_xzzz_xxxxxzz_0, \
                                         tg_xzzz_xxxxxzz_1, tg_xzzz_xxxxyy_1, tg_xzzz_xxxxyyy_0, tg_xzzz_xxxxyyy_1, \
                                         tg_xzzz_xxxxyyz_0, tg_xzzz_xxxxyyz_1, tg_xzzz_xxxxyz_1, tg_xzzz_xxxxyzz_0, \
                                         tg_xzzz_xxxxyzz_1, tg_xzzz_xxxxzz_1, tg_xzzz_xxxxzzz_0, tg_xzzz_xxxxzzz_1, \
                                         tg_xzzz_xxxyyy_1, tg_xzzz_xxxyyyy_0, tg_xzzz_xxxyyyy_1, tg_xzzz_xxxyyyz_0, \
                                         tg_xzzz_xxxyyyz_1, tg_xzzz_xxxyyz_1, tg_xzzz_xxxyyzz_0, tg_xzzz_xxxyyzz_1, \
                                         tg_xzzz_xxxyzz_1, tg_xzzz_xxxyzzz_0, tg_xzzz_xxxyzzz_1, tg_xzzz_xxxzzz_1, \
                                         tg_xzzz_xxxzzzz_0, tg_xzzz_xxxzzzz_1, tg_xzzz_xxyyyy_1, tg_xzzz_xxyyyyy_0, \
                                         tg_xzzz_xxyyyyy_1, tg_xzzz_xxyyyyz_0, tg_xzzz_xxyyyyz_1, tg_xzzz_xxyyyz_1, \
                                         tg_xzzz_xxyyyzz_0, tg_xzzz_xxyyyzz_1, tg_xzzz_xxyyzz_1, tg_xzzz_xxyyzzz_0, \
                                         tg_xzzz_xxyyzzz_1, tg_xzzz_xxyzzz_1, tg_xzzz_xxyzzzz_0, tg_xzzz_xxyzzzz_1, \
                                         tg_xzzz_xxzzzz_1, tg_xzzz_xxzzzzz_0, tg_xzzz_xxzzzzz_1, tg_xzzz_xyyyyy_1, \
                                         tg_xzzz_xyyyyyy_0, tg_xzzz_xyyyyyy_1, tg_xzzz_xyyyyyz_0, tg_xzzz_xyyyyyz_1, \
                                         tg_xzzz_xyyyyz_1, tg_xzzz_xyyyyzz_0, tg_xzzz_xyyyyzz_1, tg_xzzz_xyyyzz_1, \
                                         tg_xzzz_xyyyzzz_0, tg_xzzz_xyyyzzz_1, tg_xzzz_xyyzzz_1, tg_xzzz_xyyzzzz_0, \
                                         tg_xzzz_xyyzzzz_1, tg_xzzz_xyzzzz_1, tg_xzzz_xyzzzzz_0, tg_xzzz_xyzzzzz_1, \
                                         tg_xzzz_xzzzzz_1, tg_xzzz_xzzzzzz_0, tg_xzzz_xzzzzzz_1, tg_xzzz_yyyyyy_1, \
                                         tg_xzzz_yyyyyyy_0, tg_xzzz_yyyyyyy_1, tg_xzzz_yyyyyyz_0, tg_xzzz_yyyyyyz_1, \
                                         tg_xzzz_yyyyyz_1, tg_xzzz_yyyyyzz_0, tg_xzzz_yyyyyzz_1, tg_xzzz_yyyyzz_1, \
                                         tg_xzzz_yyyyzzz_0, tg_xzzz_yyyyzzz_1, tg_xzzz_yyyzzz_1, tg_xzzz_yyyzzzz_0, \
                                         tg_xzzz_yyyzzzz_1, tg_xzzz_yyzzzz_1, tg_xzzz_yyzzzzz_0, tg_xzzz_yyzzzzz_1, \
                                         tg_xzzz_yzzzzz_1, tg_xzzz_yzzzzzz_0, tg_xzzz_yzzzzzz_1, tg_xzzz_zzzzzz_1, \
                                         tg_xzzz_zzzzzzz_0, tg_xzzz_zzzzzzz_1, tg_yyyy_xxxxxx_1, tg_yyyy_xxxxxxx_0, \
                                         tg_yyyy_xxxxxxx_1, tg_yyyy_xxxxxxy_0, tg_yyyy_xxxxxxy_1, tg_yyyy_xxxxxxz_0, \
                                         tg_yyyy_xxxxxxz_1, tg_yyyy_xxxxxy_1, tg_yyyy_xxxxxyy_0, tg_yyyy_xxxxxyy_1, \
                                         tg_yyyy_xxxxxyz_0, tg_yyyy_xxxxxyz_1, tg_yyyy_xxxxxz_1, tg_yyyy_xxxxxzz_0, \
                                         tg_yyyy_xxxxxzz_1, tg_yyyy_xxxxyy_1, tg_yyyy_xxxxyyy_0, tg_yyyy_xxxxyyy_1, \
                                         tg_yyyy_xxxxyyz_0, tg_yyyy_xxxxyyz_1, tg_yyyy_xxxxyz_1, tg_yyyy_xxxxyzz_0, \
                                         tg_yyyy_xxxxyzz_1, tg_yyyy_xxxxzz_1, tg_yyyy_xxxxzzz_0, tg_yyyy_xxxxzzz_1, \
                                         tg_yyyy_xxxyyy_1, tg_yyyy_xxxyyyy_0, tg_yyyy_xxxyyyy_1, tg_yyyy_xxxyyyz_0, \
                                         tg_yyyy_xxxyyyz_1, tg_yyyy_xxxyyz_1, tg_yyyy_xxxyyzz_0, tg_yyyy_xxxyyzz_1, \
                                         tg_yyyy_xxxyzz_1, tg_yyyy_xxxyzzz_0, tg_yyyy_xxxyzzz_1, tg_yyyy_xxxzzz_1, \
                                         tg_yyyy_xxxzzzz_0, tg_yyyy_xxxzzzz_1, tg_yyyy_xxyyyy_1, tg_yyyy_xxyyyyy_0, \
                                         tg_yyyy_xxyyyyy_1, tg_yyyy_xxyyyyz_0, tg_yyyy_xxyyyyz_1, tg_yyyy_xxyyyz_1, \
                                         tg_yyyy_xxyyyzz_0, tg_yyyy_xxyyyzz_1, tg_yyyy_xxyyzz_1, tg_yyyy_xxyyzzz_0, \
                                         tg_yyyy_xxyyzzz_1, tg_yyyy_xxyzzz_1, tg_yyyy_xxyzzzz_0, tg_yyyy_xxyzzzz_1, \
                                         tg_yyyy_xxzzzz_1, tg_yyyy_xyyyyy_1, tg_yyyy_xyyyyz_1, tg_yyyy_xyyyzz_1, \
                                         tg_yyyy_xyyzzz_1, tg_yyyy_xyzzzz_1, tg_yyz_yyzzzzz_0, tg_yyz_yyzzzzz_1, \
                                         tg_yyz_yzzzzzz_0, tg_yyz_yzzzzzz_1, tg_yyz_zzzzzzz_0, tg_yyz_zzzzzzz_1, \
                                         tg_yzz_xxxxxxx_0, tg_yzz_xxxxxxx_1, tg_yzz_xxxxxxy_0, tg_yzz_xxxxxxy_1, \
                                         tg_yzz_xxxxxxz_0, tg_yzz_xxxxxxz_1, tg_yzz_xxxxxyy_0, tg_yzz_xxxxxyy_1, \
                                         tg_yzz_xxxxxyz_0, tg_yzz_xxxxxyz_1, tg_yzz_xxxxxzz_0, tg_yzz_xxxxxzz_1, \
                                         tg_yzz_xxxxyyy_0, tg_yzz_xxxxyyy_1, tg_yzz_xxxxyyz_0, tg_yzz_xxxxyyz_1, \
                                         tg_yzz_xxxxyzz_0, tg_yzz_xxxxyzz_1, tg_yzz_xxxxzzz_0, tg_yzz_xxxxzzz_1, \
                                         tg_yzz_xxxyyyy_0, tg_yzz_xxxyyyy_1, tg_yzz_xxxyyyz_0, tg_yzz_xxxyyyz_1, \
                                         tg_yzz_xxxyyzz_0, tg_yzz_xxxyyzz_1, tg_yzz_xxxyzzz_0, tg_yzz_xxxyzzz_1, \
                                         tg_yzz_xxxzzzz_0, tg_yzz_xxxzzzz_1, tg_yzz_xxyyyyy_0, tg_yzz_xxyyyyy_1, \
                                         tg_yzz_xxyyyyz_0, tg_yzz_xxyyyyz_1, tg_yzz_xxyyyzz_0, tg_yzz_xxyyyzz_1, \
                                         tg_yzz_xxyyzzz_0, tg_yzz_xxyyzzz_1, tg_yzz_xxyzzzz_0, tg_yzz_xxyzzzz_1, \
                                         tg_yzz_xxzzzzz_0, tg_yzz_xxzzzzz_1, tg_yzz_xyyyyyy_0, tg_yzz_xyyyyyy_1, \
                                         tg_yzz_xyyyyyz_0, tg_yzz_xyyyyyz_1, tg_yzz_xyyyyzz_0, tg_yzz_xyyyyzz_1, \
                                         tg_yzz_xyyyzzz_0, tg_yzz_xyyyzzz_1, tg_yzz_xyyzzzz_0, tg_yzz_xyyzzzz_1, \
                                         tg_yzz_xyzzzzz_0, tg_yzz_xyzzzzz_1, tg_yzz_xzzzzzz_0, tg_yzz_xzzzzzz_1, \
                                         tg_yzz_yyyyyyy_0, tg_yzz_yyyyyyy_1, tg_yzz_yyyyyyz_0, tg_yzz_yyyyyyz_1, \
                                         tg_yzz_yyyyyzz_0, tg_yzz_yyyyyzz_1, tg_yzz_yyyyzzz_0, tg_yzz_yyyyzzz_1, \
                                         tg_yzz_yyyzzzz_0, tg_yzz_yyyzzzz_1, tg_yzz_yyzzzzz_0, tg_yzz_yyzzzzz_1, \
                                         tg_yzz_yzzzzzz_0, tg_yzz_yzzzzzz_1, tg_yzz_zzzzzzz_0, tg_yzz_zzzzzzz_1, \
                                         tg_zzz_xxxxxxx_0, tg_zzz_xxxxxxx_1, tg_zzz_xxxxxxy_0, tg_zzz_xxxxxxy_1, \
                                         tg_zzz_xxxxxxz_0, tg_zzz_xxxxxxz_1, tg_zzz_xxxxxyy_0, tg_zzz_xxxxxyy_1, \
                                         tg_zzz_xxxxxyz_0, tg_zzz_xxxxxyz_1, tg_zzz_xxxxxzz_0, tg_zzz_xxxxxzz_1, \
                                         tg_zzz_xxxxyyy_0, tg_zzz_xxxxyyy_1, tg_zzz_xxxxyyz_0, tg_zzz_xxxxyyz_1, \
                                         tg_zzz_xxxxyzz_0, tg_zzz_xxxxyzz_1, tg_zzz_xxxxzzz_0, tg_zzz_xxxxzzz_1, \
                                         tg_zzz_xxxyyyy_0, tg_zzz_xxxyyyy_1, tg_zzz_xxxyyyz_0, tg_zzz_xxxyyyz_1, \
                                         tg_zzz_xxxyyzz_0, tg_zzz_xxxyyzz_1, tg_zzz_xxxyzzz_0, tg_zzz_xxxyzzz_1, \
                                         tg_zzz_xxxzzzz_0, tg_zzz_xxxzzzz_1, tg_zzz_xxyyyyy_0, tg_zzz_xxyyyyy_1, \
                                         tg_zzz_xxyyyyz_0, tg_zzz_xxyyyyz_1, tg_zzz_xxyyyzz_0, tg_zzz_xxyyyzz_1, \
                                         tg_zzz_xxyyzzz_0, tg_zzz_xxyyzzz_1, tg_zzz_xxyzzzz_0, tg_zzz_xxyzzzz_1, \
                                         tg_zzz_xxzzzzz_0, tg_zzz_xxzzzzz_1, tg_zzz_xyyyyyy_0, tg_zzz_xyyyyyy_1, \
                                         tg_zzz_xyyyyyz_0, tg_zzz_xyyyyyz_1, tg_zzz_xyyyyzz_0, tg_zzz_xyyyyzz_1, \
                                         tg_zzz_xyyyzzz_0, tg_zzz_xyyyzzz_1, tg_zzz_xyyzzzz_0, tg_zzz_xyyzzzz_1, \
                                         tg_zzz_xyzzzzz_0, tg_zzz_xyzzzzz_1, tg_zzz_xzzzzzz_0, tg_zzz_xzzzzzz_1, \
                                         tg_zzz_yyyyyyy_0, tg_zzz_yyyyyyy_1, tg_zzz_yyyyyyz_0, tg_zzz_yyyyyyz_1, \
                                         tg_zzz_yyyyyzz_0, tg_zzz_yyyyyzz_1, tg_zzz_yyyyzzz_0, tg_zzz_yyyyzzz_1, \
                                         tg_zzz_yyyzzzz_0, tg_zzz_yyyzzzz_1, tg_zzz_yyzzzzz_0, tg_zzz_yyzzzzz_1, \
                                         tg_zzz_yzzzzzz_0, tg_zzz_yzzzzzz_1, tg_zzz_zzzzzzz_0, tg_zzz_zzzzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxyyz_yyzzzzz_0[j] = pb_x * tg_xyyz_yyzzzzz_0[j] + wp_x[j] * tg_xyyz_yyzzzzz_1[j] + 0.5 * fl1_fx * tg_yyz_yyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_yyzzzzz_1[j];

                    tg_xxyyz_yzzzzzz_0[j] = pb_x * tg_xyyz_yzzzzzz_0[j] + wp_x[j] * tg_xyyz_yzzzzzz_1[j] + 0.5 * fl1_fx * tg_yyz_yzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_yzzzzzz_1[j];

                    tg_xxyyz_zzzzzzz_0[j] = pb_x * tg_xyyz_zzzzzzz_0[j] + wp_x[j] * tg_xyyz_zzzzzzz_1[j] + 0.5 * fl1_fx * tg_yyz_zzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_zzzzzzz_1[j];

                    tg_xxyzz_xxxxxxx_0[j] = pb_x * tg_xyzz_xxxxxxx_0[j] + wp_x[j] * tg_xyzz_xxxxxxx_1[j] + 0.5 * fl1_fx * tg_yzz_xxxxxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxxxxxx_1[j] + 3.5 * fl1_fxn * tg_xyzz_xxxxxx_1[j];

                    tg_xxyzz_xxxxxxy_0[j] = pb_x * tg_xyzz_xxxxxxy_0[j] + wp_x[j] * tg_xyzz_xxxxxxy_1[j] + 0.5 * fl1_fx * tg_yzz_xxxxxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxxxxxy_1[j] + 3.0 * fl1_fxn * tg_xyzz_xxxxxy_1[j];

                    tg_xxyzz_xxxxxxz_0[j] = pb_x * tg_xyzz_xxxxxxz_0[j] + wp_x[j] * tg_xyzz_xxxxxxz_1[j] + 0.5 * fl1_fx * tg_yzz_xxxxxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxxxxxz_1[j] + 3.0 * fl1_fxn * tg_xyzz_xxxxxz_1[j];

                    tg_xxyzz_xxxxxyy_0[j] = pb_x * tg_xyzz_xxxxxyy_0[j] + wp_x[j] * tg_xyzz_xxxxxyy_1[j] + 0.5 * fl1_fx * tg_yzz_xxxxxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxxxxyy_1[j] + 2.5 * fl1_fxn * tg_xyzz_xxxxyy_1[j];

                    tg_xxyzz_xxxxxyz_0[j] = pb_x * tg_xyzz_xxxxxyz_0[j] + wp_x[j] * tg_xyzz_xxxxxyz_1[j] + 0.5 * fl1_fx * tg_yzz_xxxxxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxxxxyz_1[j] + 2.5 * fl1_fxn * tg_xyzz_xxxxyz_1[j];

                    tg_xxyzz_xxxxxzz_0[j] = pb_x * tg_xyzz_xxxxxzz_0[j] + wp_x[j] * tg_xyzz_xxxxxzz_1[j] + 0.5 * fl1_fx * tg_yzz_xxxxxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxxxxzz_1[j] + 2.5 * fl1_fxn * tg_xyzz_xxxxzz_1[j];

                    tg_xxyzz_xxxxyyy_0[j] = pb_x * tg_xyzz_xxxxyyy_0[j] + wp_x[j] * tg_xyzz_xxxxyyy_1[j] + 0.5 * fl1_fx * tg_yzz_xxxxyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxxxyyy_1[j] + 2.0 * fl1_fxn * tg_xyzz_xxxyyy_1[j];

                    tg_xxyzz_xxxxyyz_0[j] = pb_x * tg_xyzz_xxxxyyz_0[j] + wp_x[j] * tg_xyzz_xxxxyyz_1[j] + 0.5 * fl1_fx * tg_yzz_xxxxyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxxxyyz_1[j] + 2.0 * fl1_fxn * tg_xyzz_xxxyyz_1[j];

                    tg_xxyzz_xxxxyzz_0[j] = pb_x * tg_xyzz_xxxxyzz_0[j] + wp_x[j] * tg_xyzz_xxxxyzz_1[j] + 0.5 * fl1_fx * tg_yzz_xxxxyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxxxyzz_1[j] + 2.0 * fl1_fxn * tg_xyzz_xxxyzz_1[j];

                    tg_xxyzz_xxxxzzz_0[j] = pb_x * tg_xyzz_xxxxzzz_0[j] + wp_x[j] * tg_xyzz_xxxxzzz_1[j] + 0.5 * fl1_fx * tg_yzz_xxxxzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxxxzzz_1[j] + 2.0 * fl1_fxn * tg_xyzz_xxxzzz_1[j];

                    tg_xxyzz_xxxyyyy_0[j] = pb_x * tg_xyzz_xxxyyyy_0[j] + wp_x[j] * tg_xyzz_xxxyyyy_1[j] + 0.5 * fl1_fx * tg_yzz_xxxyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxxyyyy_1[j] + 1.5 * fl1_fxn * tg_xyzz_xxyyyy_1[j];

                    tg_xxyzz_xxxyyyz_0[j] = pb_x * tg_xyzz_xxxyyyz_0[j] + wp_x[j] * tg_xyzz_xxxyyyz_1[j] + 0.5 * fl1_fx * tg_yzz_xxxyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_xyzz_xxyyyz_1[j];

                    tg_xxyzz_xxxyyzz_0[j] = pb_x * tg_xyzz_xxxyyzz_0[j] + wp_x[j] * tg_xyzz_xxxyyzz_1[j] + 0.5 * fl1_fx * tg_yzz_xxxyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxxyyzz_1[j] + 1.5 * fl1_fxn * tg_xyzz_xxyyzz_1[j];

                    tg_xxyzz_xxxyzzz_0[j] = pb_x * tg_xyzz_xxxyzzz_0[j] + wp_x[j] * tg_xyzz_xxxyzzz_1[j] + 0.5 * fl1_fx * tg_yzz_xxxyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxxyzzz_1[j] + 1.5 * fl1_fxn * tg_xyzz_xxyzzz_1[j];

                    tg_xxyzz_xxxzzzz_0[j] = pb_x * tg_xyzz_xxxzzzz_0[j] + wp_x[j] * tg_xyzz_xxxzzzz_1[j] + 0.5 * fl1_fx * tg_yzz_xxxzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxxzzzz_1[j] + 1.5 * fl1_fxn * tg_xyzz_xxzzzz_1[j];

                    tg_xxyzz_xxyyyyy_0[j] = pb_x * tg_xyzz_xxyyyyy_0[j] + wp_x[j] * tg_xyzz_xxyyyyy_1[j] + 0.5 * fl1_fx * tg_yzz_xxyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxyyyyy_1[j] + fl1_fxn * tg_xyzz_xyyyyy_1[j];

                    tg_xxyzz_xxyyyyz_0[j] = pb_x * tg_xyzz_xxyyyyz_0[j] + wp_x[j] * tg_xyzz_xxyyyyz_1[j] + 0.5 * fl1_fx * tg_yzz_xxyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxyyyyz_1[j] + fl1_fxn * tg_xyzz_xyyyyz_1[j];

                    tg_xxyzz_xxyyyzz_0[j] = pb_x * tg_xyzz_xxyyyzz_0[j] + wp_x[j] * tg_xyzz_xxyyyzz_1[j] + 0.5 * fl1_fx * tg_yzz_xxyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxyyyzz_1[j] + fl1_fxn * tg_xyzz_xyyyzz_1[j];

                    tg_xxyzz_xxyyzzz_0[j] = pb_x * tg_xyzz_xxyyzzz_0[j] + wp_x[j] * tg_xyzz_xxyyzzz_1[j] + 0.5 * fl1_fx * tg_yzz_xxyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxyyzzz_1[j] + fl1_fxn * tg_xyzz_xyyzzz_1[j];

                    tg_xxyzz_xxyzzzz_0[j] = pb_x * tg_xyzz_xxyzzzz_0[j] + wp_x[j] * tg_xyzz_xxyzzzz_1[j] + 0.5 * fl1_fx * tg_yzz_xxyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxyzzzz_1[j] + fl1_fxn * tg_xyzz_xyzzzz_1[j];

                    tg_xxyzz_xxzzzzz_0[j] = pb_x * tg_xyzz_xxzzzzz_0[j] + wp_x[j] * tg_xyzz_xxzzzzz_1[j] + 0.5 * fl1_fx * tg_yzz_xxzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxzzzzz_1[j] + fl1_fxn * tg_xyzz_xzzzzz_1[j];

                    tg_xxyzz_xyyyyyy_0[j] = pb_x * tg_xyzz_xyyyyyy_0[j] + wp_x[j] * tg_xyzz_xyyyyyy_1[j] + 0.5 * fl1_fx * tg_yzz_xyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xyyyyyy_1[j] + 0.5 * fl1_fxn * tg_xyzz_yyyyyy_1[j];

                    tg_xxyzz_xyyyyyz_0[j] = pb_x * tg_xyzz_xyyyyyz_0[j] + wp_x[j] * tg_xyzz_xyyyyyz_1[j] + 0.5 * fl1_fx * tg_yzz_xyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xyyyyyz_1[j] + 0.5 * fl1_fxn * tg_xyzz_yyyyyz_1[j];

                    tg_xxyzz_xyyyyzz_0[j] = pb_x * tg_xyzz_xyyyyzz_0[j] + wp_x[j] * tg_xyzz_xyyyyzz_1[j] + 0.5 * fl1_fx * tg_yzz_xyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xyyyyzz_1[j] + 0.5 * fl1_fxn * tg_xyzz_yyyyzz_1[j];

                    tg_xxyzz_xyyyzzz_0[j] = pb_x * tg_xyzz_xyyyzzz_0[j] + wp_x[j] * tg_xyzz_xyyyzzz_1[j] + 0.5 * fl1_fx * tg_yzz_xyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xyyyzzz_1[j] + 0.5 * fl1_fxn * tg_xyzz_yyyzzz_1[j];

                    tg_xxyzz_xyyzzzz_0[j] = pb_x * tg_xyzz_xyyzzzz_0[j] + wp_x[j] * tg_xyzz_xyyzzzz_1[j] + 0.5 * fl1_fx * tg_yzz_xyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xyyzzzz_1[j] + 0.5 * fl1_fxn * tg_xyzz_yyzzzz_1[j];

                    tg_xxyzz_xyzzzzz_0[j] = pb_x * tg_xyzz_xyzzzzz_0[j] + wp_x[j] * tg_xyzz_xyzzzzz_1[j] + 0.5 * fl1_fx * tg_yzz_xyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_xyzz_yzzzzz_1[j];

                    tg_xxyzz_xzzzzzz_0[j] = pb_x * tg_xyzz_xzzzzzz_0[j] + wp_x[j] * tg_xyzz_xzzzzzz_1[j] + 0.5 * fl1_fx * tg_yzz_xzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xyzz_zzzzzz_1[j];

                    tg_xxyzz_yyyyyyy_0[j] = pb_x * tg_xyzz_yyyyyyy_0[j] + wp_x[j] * tg_xyzz_yyyyyyy_1[j] + 0.5 * fl1_fx * tg_yzz_yyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_yyyyyyy_1[j];

                    tg_xxyzz_yyyyyyz_0[j] = pb_x * tg_xyzz_yyyyyyz_0[j] + wp_x[j] * tg_xyzz_yyyyyyz_1[j] + 0.5 * fl1_fx * tg_yzz_yyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_yyyyyyz_1[j];

                    tg_xxyzz_yyyyyzz_0[j] = pb_x * tg_xyzz_yyyyyzz_0[j] + wp_x[j] * tg_xyzz_yyyyyzz_1[j] + 0.5 * fl1_fx * tg_yzz_yyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_yyyyyzz_1[j];

                    tg_xxyzz_yyyyzzz_0[j] = pb_x * tg_xyzz_yyyyzzz_0[j] + wp_x[j] * tg_xyzz_yyyyzzz_1[j] + 0.5 * fl1_fx * tg_yzz_yyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_yyyyzzz_1[j];

                    tg_xxyzz_yyyzzzz_0[j] = pb_x * tg_xyzz_yyyzzzz_0[j] + wp_x[j] * tg_xyzz_yyyzzzz_1[j] + 0.5 * fl1_fx * tg_yzz_yyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_yyyzzzz_1[j];

                    tg_xxyzz_yyzzzzz_0[j] = pb_x * tg_xyzz_yyzzzzz_0[j] + wp_x[j] * tg_xyzz_yyzzzzz_1[j] + 0.5 * fl1_fx * tg_yzz_yyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_yyzzzzz_1[j];

                    tg_xxyzz_yzzzzzz_0[j] = pb_x * tg_xyzz_yzzzzzz_0[j] + wp_x[j] * tg_xyzz_yzzzzzz_1[j] + 0.5 * fl1_fx * tg_yzz_yzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_yzzzzzz_1[j];

                    tg_xxyzz_zzzzzzz_0[j] = pb_x * tg_xyzz_zzzzzzz_0[j] + wp_x[j] * tg_xyzz_zzzzzzz_1[j] + 0.5 * fl1_fx * tg_yzz_zzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_zzzzzzz_1[j];

                    tg_xxzzz_xxxxxxx_0[j] = pb_x * tg_xzzz_xxxxxxx_0[j] + wp_x[j] * tg_xzzz_xxxxxxx_1[j] + 0.5 * fl1_fx * tg_zzz_xxxxxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxxxxx_1[j] + 3.5 * fl1_fxn * tg_xzzz_xxxxxx_1[j];

                    tg_xxzzz_xxxxxxy_0[j] = pb_x * tg_xzzz_xxxxxxy_0[j] + wp_x[j] * tg_xzzz_xxxxxxy_1[j] + 0.5 * fl1_fx * tg_zzz_xxxxxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxxxxy_1[j] + 3.0 * fl1_fxn * tg_xzzz_xxxxxy_1[j];

                    tg_xxzzz_xxxxxxz_0[j] = pb_x * tg_xzzz_xxxxxxz_0[j] + wp_x[j] * tg_xzzz_xxxxxxz_1[j] + 0.5 * fl1_fx * tg_zzz_xxxxxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxxxxz_1[j] + 3.0 * fl1_fxn * tg_xzzz_xxxxxz_1[j];

                    tg_xxzzz_xxxxxyy_0[j] = pb_x * tg_xzzz_xxxxxyy_0[j] + wp_x[j] * tg_xzzz_xxxxxyy_1[j] + 0.5 * fl1_fx * tg_zzz_xxxxxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxxxyy_1[j] + 2.5 * fl1_fxn * tg_xzzz_xxxxyy_1[j];

                    tg_xxzzz_xxxxxyz_0[j] = pb_x * tg_xzzz_xxxxxyz_0[j] + wp_x[j] * tg_xzzz_xxxxxyz_1[j] + 0.5 * fl1_fx * tg_zzz_xxxxxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxxxyz_1[j] + 2.5 * fl1_fxn * tg_xzzz_xxxxyz_1[j];

                    tg_xxzzz_xxxxxzz_0[j] = pb_x * tg_xzzz_xxxxxzz_0[j] + wp_x[j] * tg_xzzz_xxxxxzz_1[j] + 0.5 * fl1_fx * tg_zzz_xxxxxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxxxzz_1[j] + 2.5 * fl1_fxn * tg_xzzz_xxxxzz_1[j];

                    tg_xxzzz_xxxxyyy_0[j] = pb_x * tg_xzzz_xxxxyyy_0[j] + wp_x[j] * tg_xzzz_xxxxyyy_1[j] + 0.5 * fl1_fx * tg_zzz_xxxxyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxxyyy_1[j] + 2.0 * fl1_fxn * tg_xzzz_xxxyyy_1[j];

                    tg_xxzzz_xxxxyyz_0[j] = pb_x * tg_xzzz_xxxxyyz_0[j] + wp_x[j] * tg_xzzz_xxxxyyz_1[j] + 0.5 * fl1_fx * tg_zzz_xxxxyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxxyyz_1[j] + 2.0 * fl1_fxn * tg_xzzz_xxxyyz_1[j];

                    tg_xxzzz_xxxxyzz_0[j] = pb_x * tg_xzzz_xxxxyzz_0[j] + wp_x[j] * tg_xzzz_xxxxyzz_1[j] + 0.5 * fl1_fx * tg_zzz_xxxxyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxxyzz_1[j] + 2.0 * fl1_fxn * tg_xzzz_xxxyzz_1[j];

                    tg_xxzzz_xxxxzzz_0[j] = pb_x * tg_xzzz_xxxxzzz_0[j] + wp_x[j] * tg_xzzz_xxxxzzz_1[j] + 0.5 * fl1_fx * tg_zzz_xxxxzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxxzzz_1[j] + 2.0 * fl1_fxn * tg_xzzz_xxxzzz_1[j];

                    tg_xxzzz_xxxyyyy_0[j] = pb_x * tg_xzzz_xxxyyyy_0[j] + wp_x[j] * tg_xzzz_xxxyyyy_1[j] + 0.5 * fl1_fx * tg_zzz_xxxyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxyyyy_1[j] + 1.5 * fl1_fxn * tg_xzzz_xxyyyy_1[j];

                    tg_xxzzz_xxxyyyz_0[j] = pb_x * tg_xzzz_xxxyyyz_0[j] + wp_x[j] * tg_xzzz_xxxyyyz_1[j] + 0.5 * fl1_fx * tg_zzz_xxxyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_xzzz_xxyyyz_1[j];

                    tg_xxzzz_xxxyyzz_0[j] = pb_x * tg_xzzz_xxxyyzz_0[j] + wp_x[j] * tg_xzzz_xxxyyzz_1[j] + 0.5 * fl1_fx * tg_zzz_xxxyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxyyzz_1[j] + 1.5 * fl1_fxn * tg_xzzz_xxyyzz_1[j];

                    tg_xxzzz_xxxyzzz_0[j] = pb_x * tg_xzzz_xxxyzzz_0[j] + wp_x[j] * tg_xzzz_xxxyzzz_1[j] + 0.5 * fl1_fx * tg_zzz_xxxyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxyzzz_1[j] + 1.5 * fl1_fxn * tg_xzzz_xxyzzz_1[j];

                    tg_xxzzz_xxxzzzz_0[j] = pb_x * tg_xzzz_xxxzzzz_0[j] + wp_x[j] * tg_xzzz_xxxzzzz_1[j] + 0.5 * fl1_fx * tg_zzz_xxxzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxzzzz_1[j] + 1.5 * fl1_fxn * tg_xzzz_xxzzzz_1[j];

                    tg_xxzzz_xxyyyyy_0[j] = pb_x * tg_xzzz_xxyyyyy_0[j] + wp_x[j] * tg_xzzz_xxyyyyy_1[j] + 0.5 * fl1_fx * tg_zzz_xxyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxyyyyy_1[j] + fl1_fxn * tg_xzzz_xyyyyy_1[j];

                    tg_xxzzz_xxyyyyz_0[j] = pb_x * tg_xzzz_xxyyyyz_0[j] + wp_x[j] * tg_xzzz_xxyyyyz_1[j] + 0.5 * fl1_fx * tg_zzz_xxyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxyyyyz_1[j] + fl1_fxn * tg_xzzz_xyyyyz_1[j];

                    tg_xxzzz_xxyyyzz_0[j] = pb_x * tg_xzzz_xxyyyzz_0[j] + wp_x[j] * tg_xzzz_xxyyyzz_1[j] + 0.5 * fl1_fx * tg_zzz_xxyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxyyyzz_1[j] + fl1_fxn * tg_xzzz_xyyyzz_1[j];

                    tg_xxzzz_xxyyzzz_0[j] = pb_x * tg_xzzz_xxyyzzz_0[j] + wp_x[j] * tg_xzzz_xxyyzzz_1[j] + 0.5 * fl1_fx * tg_zzz_xxyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxyyzzz_1[j] + fl1_fxn * tg_xzzz_xyyzzz_1[j];

                    tg_xxzzz_xxyzzzz_0[j] = pb_x * tg_xzzz_xxyzzzz_0[j] + wp_x[j] * tg_xzzz_xxyzzzz_1[j] + 0.5 * fl1_fx * tg_zzz_xxyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxyzzzz_1[j] + fl1_fxn * tg_xzzz_xyzzzz_1[j];

                    tg_xxzzz_xxzzzzz_0[j] = pb_x * tg_xzzz_xxzzzzz_0[j] + wp_x[j] * tg_xzzz_xxzzzzz_1[j] + 0.5 * fl1_fx * tg_zzz_xxzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxzzzzz_1[j] + fl1_fxn * tg_xzzz_xzzzzz_1[j];

                    tg_xxzzz_xyyyyyy_0[j] = pb_x * tg_xzzz_xyyyyyy_0[j] + wp_x[j] * tg_xzzz_xyyyyyy_1[j] + 0.5 * fl1_fx * tg_zzz_xyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xyyyyyy_1[j] + 0.5 * fl1_fxn * tg_xzzz_yyyyyy_1[j];

                    tg_xxzzz_xyyyyyz_0[j] = pb_x * tg_xzzz_xyyyyyz_0[j] + wp_x[j] * tg_xzzz_xyyyyyz_1[j] + 0.5 * fl1_fx * tg_zzz_xyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xyyyyyz_1[j] + 0.5 * fl1_fxn * tg_xzzz_yyyyyz_1[j];

                    tg_xxzzz_xyyyyzz_0[j] = pb_x * tg_xzzz_xyyyyzz_0[j] + wp_x[j] * tg_xzzz_xyyyyzz_1[j] + 0.5 * fl1_fx * tg_zzz_xyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xyyyyzz_1[j] + 0.5 * fl1_fxn * tg_xzzz_yyyyzz_1[j];

                    tg_xxzzz_xyyyzzz_0[j] = pb_x * tg_xzzz_xyyyzzz_0[j] + wp_x[j] * tg_xzzz_xyyyzzz_1[j] + 0.5 * fl1_fx * tg_zzz_xyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xyyyzzz_1[j] + 0.5 * fl1_fxn * tg_xzzz_yyyzzz_1[j];

                    tg_xxzzz_xyyzzzz_0[j] = pb_x * tg_xzzz_xyyzzzz_0[j] + wp_x[j] * tg_xzzz_xyyzzzz_1[j] + 0.5 * fl1_fx * tg_zzz_xyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xyyzzzz_1[j] + 0.5 * fl1_fxn * tg_xzzz_yyzzzz_1[j];

                    tg_xxzzz_xyzzzzz_0[j] = pb_x * tg_xzzz_xyzzzzz_0[j] + wp_x[j] * tg_xzzz_xyzzzzz_1[j] + 0.5 * fl1_fx * tg_zzz_xyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_xzzz_yzzzzz_1[j];

                    tg_xxzzz_xzzzzzz_0[j] = pb_x * tg_xzzz_xzzzzzz_0[j] + wp_x[j] * tg_xzzz_xzzzzzz_1[j] + 0.5 * fl1_fx * tg_zzz_xzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xzzz_zzzzzz_1[j];

                    tg_xxzzz_yyyyyyy_0[j] = pb_x * tg_xzzz_yyyyyyy_0[j] + wp_x[j] * tg_xzzz_yyyyyyy_1[j] + 0.5 * fl1_fx * tg_zzz_yyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_yyyyyyy_1[j];

                    tg_xxzzz_yyyyyyz_0[j] = pb_x * tg_xzzz_yyyyyyz_0[j] + wp_x[j] * tg_xzzz_yyyyyyz_1[j] + 0.5 * fl1_fx * tg_zzz_yyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_yyyyyyz_1[j];

                    tg_xxzzz_yyyyyzz_0[j] = pb_x * tg_xzzz_yyyyyzz_0[j] + wp_x[j] * tg_xzzz_yyyyyzz_1[j] + 0.5 * fl1_fx * tg_zzz_yyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_yyyyyzz_1[j];

                    tg_xxzzz_yyyyzzz_0[j] = pb_x * tg_xzzz_yyyyzzz_0[j] + wp_x[j] * tg_xzzz_yyyyzzz_1[j] + 0.5 * fl1_fx * tg_zzz_yyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_yyyyzzz_1[j];

                    tg_xxzzz_yyyzzzz_0[j] = pb_x * tg_xzzz_yyyzzzz_0[j] + wp_x[j] * tg_xzzz_yyyzzzz_1[j] + 0.5 * fl1_fx * tg_zzz_yyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_yyyzzzz_1[j];

                    tg_xxzzz_yyzzzzz_0[j] = pb_x * tg_xzzz_yyzzzzz_0[j] + wp_x[j] * tg_xzzz_yyzzzzz_1[j] + 0.5 * fl1_fx * tg_zzz_yyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_yyzzzzz_1[j];

                    tg_xxzzz_yzzzzzz_0[j] = pb_x * tg_xzzz_yzzzzzz_0[j] + wp_x[j] * tg_xzzz_yzzzzzz_1[j] + 0.5 * fl1_fx * tg_zzz_yzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_yzzzzzz_1[j];

                    tg_xxzzz_zzzzzzz_0[j] = pb_x * tg_xzzz_zzzzzzz_0[j] + wp_x[j] * tg_xzzz_zzzzzzz_1[j] + 0.5 * fl1_fx * tg_zzz_zzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_zzzzzzz_1[j];

                    tg_xyyyy_xxxxxxx_0[j] = pb_x * tg_yyyy_xxxxxxx_0[j] + wp_x[j] * tg_yyyy_xxxxxxx_1[j] + 3.5 * fl1_fxn * tg_yyyy_xxxxxx_1[j];

                    tg_xyyyy_xxxxxxy_0[j] = pb_x * tg_yyyy_xxxxxxy_0[j] + wp_x[j] * tg_yyyy_xxxxxxy_1[j] + 3.0 * fl1_fxn * tg_yyyy_xxxxxy_1[j];

                    tg_xyyyy_xxxxxxz_0[j] = pb_x * tg_yyyy_xxxxxxz_0[j] + wp_x[j] * tg_yyyy_xxxxxxz_1[j] + 3.0 * fl1_fxn * tg_yyyy_xxxxxz_1[j];

                    tg_xyyyy_xxxxxyy_0[j] = pb_x * tg_yyyy_xxxxxyy_0[j] + wp_x[j] * tg_yyyy_xxxxxyy_1[j] + 2.5 * fl1_fxn * tg_yyyy_xxxxyy_1[j];

                    tg_xyyyy_xxxxxyz_0[j] = pb_x * tg_yyyy_xxxxxyz_0[j] + wp_x[j] * tg_yyyy_xxxxxyz_1[j] + 2.5 * fl1_fxn * tg_yyyy_xxxxyz_1[j];

                    tg_xyyyy_xxxxxzz_0[j] = pb_x * tg_yyyy_xxxxxzz_0[j] + wp_x[j] * tg_yyyy_xxxxxzz_1[j] + 2.5 * fl1_fxn * tg_yyyy_xxxxzz_1[j];

                    tg_xyyyy_xxxxyyy_0[j] = pb_x * tg_yyyy_xxxxyyy_0[j] + wp_x[j] * tg_yyyy_xxxxyyy_1[j] + 2.0 * fl1_fxn * tg_yyyy_xxxyyy_1[j];

                    tg_xyyyy_xxxxyyz_0[j] = pb_x * tg_yyyy_xxxxyyz_0[j] + wp_x[j] * tg_yyyy_xxxxyyz_1[j] + 2.0 * fl1_fxn * tg_yyyy_xxxyyz_1[j];

                    tg_xyyyy_xxxxyzz_0[j] = pb_x * tg_yyyy_xxxxyzz_0[j] + wp_x[j] * tg_yyyy_xxxxyzz_1[j] + 2.0 * fl1_fxn * tg_yyyy_xxxyzz_1[j];

                    tg_xyyyy_xxxxzzz_0[j] = pb_x * tg_yyyy_xxxxzzz_0[j] + wp_x[j] * tg_yyyy_xxxxzzz_1[j] + 2.0 * fl1_fxn * tg_yyyy_xxxzzz_1[j];

                    tg_xyyyy_xxxyyyy_0[j] = pb_x * tg_yyyy_xxxyyyy_0[j] + wp_x[j] * tg_yyyy_xxxyyyy_1[j] + 1.5 * fl1_fxn * tg_yyyy_xxyyyy_1[j];

                    tg_xyyyy_xxxyyyz_0[j] = pb_x * tg_yyyy_xxxyyyz_0[j] + wp_x[j] * tg_yyyy_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_yyyy_xxyyyz_1[j];

                    tg_xyyyy_xxxyyzz_0[j] = pb_x * tg_yyyy_xxxyyzz_0[j] + wp_x[j] * tg_yyyy_xxxyyzz_1[j] + 1.5 * fl1_fxn * tg_yyyy_xxyyzz_1[j];

                    tg_xyyyy_xxxyzzz_0[j] = pb_x * tg_yyyy_xxxyzzz_0[j] + wp_x[j] * tg_yyyy_xxxyzzz_1[j] + 1.5 * fl1_fxn * tg_yyyy_xxyzzz_1[j];

                    tg_xyyyy_xxxzzzz_0[j] = pb_x * tg_yyyy_xxxzzzz_0[j] + wp_x[j] * tg_yyyy_xxxzzzz_1[j] + 1.5 * fl1_fxn * tg_yyyy_xxzzzz_1[j];

                    tg_xyyyy_xxyyyyy_0[j] = pb_x * tg_yyyy_xxyyyyy_0[j] + wp_x[j] * tg_yyyy_xxyyyyy_1[j] + fl1_fxn * tg_yyyy_xyyyyy_1[j];

                    tg_xyyyy_xxyyyyz_0[j] = pb_x * tg_yyyy_xxyyyyz_0[j] + wp_x[j] * tg_yyyy_xxyyyyz_1[j] + fl1_fxn * tg_yyyy_xyyyyz_1[j];

                    tg_xyyyy_xxyyyzz_0[j] = pb_x * tg_yyyy_xxyyyzz_0[j] + wp_x[j] * tg_yyyy_xxyyyzz_1[j] + fl1_fxn * tg_yyyy_xyyyzz_1[j];

                    tg_xyyyy_xxyyzzz_0[j] = pb_x * tg_yyyy_xxyyzzz_0[j] + wp_x[j] * tg_yyyy_xxyyzzz_1[j] + fl1_fxn * tg_yyyy_xyyzzz_1[j];

                    tg_xyyyy_xxyzzzz_0[j] = pb_x * tg_yyyy_xxyzzzz_0[j] + wp_x[j] * tg_yyyy_xxyzzzz_1[j] + fl1_fxn * tg_yyyy_xyzzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSK_380_474(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (380,474)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {5, -1, -1, -1},
                                             {7, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_7_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_4_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_4_6_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {6, -1, -1, -1}, 
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

                auto tg_yyyy_xxzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 380); 

                auto tg_yyyy_xyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 381); 

                auto tg_yyyy_xyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 382); 

                auto tg_yyyy_xyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 383); 

                auto tg_yyyy_xyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 384); 

                auto tg_yyyy_xyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 385); 

                auto tg_yyyy_xyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 386); 

                auto tg_yyyy_xzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 387); 

                auto tg_yyyy_yyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 388); 

                auto tg_yyyy_yyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 389); 

                auto tg_yyyy_yyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 390); 

                auto tg_yyyy_yyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 391); 

                auto tg_yyyy_yyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 392); 

                auto tg_yyyy_yyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 393); 

                auto tg_yyyy_yzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 394); 

                auto tg_yyyy_zzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 395); 

                auto tg_yyyz_xxxxxxx_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 396); 

                auto tg_yyyz_xxxxxxy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 397); 

                auto tg_yyyz_xxxxxxz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 398); 

                auto tg_yyyz_xxxxxyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 399); 

                auto tg_yyyz_xxxxxyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 400); 

                auto tg_yyyz_xxxxxzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 401); 

                auto tg_yyyz_xxxxyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 402); 

                auto tg_yyyz_xxxxyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 403); 

                auto tg_yyyz_xxxxyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 404); 

                auto tg_yyyz_xxxxzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 405); 

                auto tg_yyyz_xxxyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 406); 

                auto tg_yyyz_xxxyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 407); 

                auto tg_yyyz_xxxyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 408); 

                auto tg_yyyz_xxxyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 409); 

                auto tg_yyyz_xxxzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 410); 

                auto tg_yyyz_xxyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 411); 

                auto tg_yyyz_xxyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 412); 

                auto tg_yyyz_xxyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 413); 

                auto tg_yyyz_xxyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 414); 

                auto tg_yyyz_xxyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 415); 

                auto tg_yyyz_xxzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 416); 

                auto tg_yyyz_xyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 417); 

                auto tg_yyyz_xyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 418); 

                auto tg_yyyz_xyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 419); 

                auto tg_yyyz_xyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 420); 

                auto tg_yyyz_xyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 421); 

                auto tg_yyyz_xyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 422); 

                auto tg_yyyz_xzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 423); 

                auto tg_yyyz_yyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 424); 

                auto tg_yyyz_yyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 425); 

                auto tg_yyyz_yyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 426); 

                auto tg_yyyz_yyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 427); 

                auto tg_yyyz_yyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 428); 

                auto tg_yyyz_yyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 429); 

                auto tg_yyyz_yzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 430); 

                auto tg_yyyz_zzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 431); 

                auto tg_yyzz_xxxxxxx_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 432); 

                auto tg_yyzz_xxxxxxy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 433); 

                auto tg_yyzz_xxxxxxz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 434); 

                auto tg_yyzz_xxxxxyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 435); 

                auto tg_yyzz_xxxxxyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 436); 

                auto tg_yyzz_xxxxxzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 437); 

                auto tg_yyzz_xxxxyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 438); 

                auto tg_yyzz_xxxxyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 439); 

                auto tg_yyzz_xxxxyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 440); 

                auto tg_yyzz_xxxxzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 441); 

                auto tg_yyzz_xxxyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 442); 

                auto tg_yyzz_xxxyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 443); 

                auto tg_yyzz_xxxyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 444); 

                auto tg_yyzz_xxxyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 445); 

                auto tg_yyzz_xxxzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 446); 

                auto tg_yyzz_xxyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 447); 

                auto tg_yyzz_xxyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 448); 

                auto tg_yyzz_xxyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 449); 

                auto tg_yyzz_xxyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 450); 

                auto tg_yyzz_xxyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 451); 

                auto tg_yyzz_xxzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 452); 

                auto tg_yyzz_xyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 453); 

                auto tg_yyzz_xyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 454); 

                auto tg_yyzz_xyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 455); 

                auto tg_yyzz_xyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 456); 

                auto tg_yyzz_xyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 457); 

                auto tg_yyzz_xyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 458); 

                auto tg_yyzz_xzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 459); 

                auto tg_yyzz_yyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 460); 

                auto tg_yyzz_yyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 461); 

                auto tg_yyzz_yyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 462); 

                auto tg_yyzz_yyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 463); 

                auto tg_yyzz_yyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 464); 

                auto tg_yyzz_yyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 465); 

                auto tg_yyzz_yzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 466); 

                auto tg_yyzz_zzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 467); 

                auto tg_yzzz_xxxxxxx_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 468); 

                auto tg_yzzz_xxxxxxy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 469); 

                auto tg_yzzz_xxxxxxz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 470); 

                auto tg_yzzz_xxxxxyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 471); 

                auto tg_yzzz_xxxxxyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 472); 

                auto tg_yzzz_xxxxxzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 473); 

                auto tg_yyyy_xxzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 380); 

                auto tg_yyyy_xyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 381); 

                auto tg_yyyy_xyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 382); 

                auto tg_yyyy_xyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 383); 

                auto tg_yyyy_xyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 384); 

                auto tg_yyyy_xyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 385); 

                auto tg_yyyy_xyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 386); 

                auto tg_yyyy_xzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 387); 

                auto tg_yyyy_yyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 388); 

                auto tg_yyyy_yyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 389); 

                auto tg_yyyy_yyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 390); 

                auto tg_yyyy_yyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 391); 

                auto tg_yyyy_yyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 392); 

                auto tg_yyyy_yyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 393); 

                auto tg_yyyy_yzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 394); 

                auto tg_yyyy_zzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 395); 

                auto tg_yyyz_xxxxxxx_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 396); 

                auto tg_yyyz_xxxxxxy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 397); 

                auto tg_yyyz_xxxxxxz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 398); 

                auto tg_yyyz_xxxxxyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 399); 

                auto tg_yyyz_xxxxxyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 400); 

                auto tg_yyyz_xxxxxzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 401); 

                auto tg_yyyz_xxxxyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 402); 

                auto tg_yyyz_xxxxyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 403); 

                auto tg_yyyz_xxxxyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 404); 

                auto tg_yyyz_xxxxzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 405); 

                auto tg_yyyz_xxxyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 406); 

                auto tg_yyyz_xxxyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 407); 

                auto tg_yyyz_xxxyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 408); 

                auto tg_yyyz_xxxyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 409); 

                auto tg_yyyz_xxxzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 410); 

                auto tg_yyyz_xxyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 411); 

                auto tg_yyyz_xxyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 412); 

                auto tg_yyyz_xxyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 413); 

                auto tg_yyyz_xxyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 414); 

                auto tg_yyyz_xxyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 415); 

                auto tg_yyyz_xxzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 416); 

                auto tg_yyyz_xyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 417); 

                auto tg_yyyz_xyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 418); 

                auto tg_yyyz_xyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 419); 

                auto tg_yyyz_xyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 420); 

                auto tg_yyyz_xyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 421); 

                auto tg_yyyz_xyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 422); 

                auto tg_yyyz_xzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 423); 

                auto tg_yyyz_yyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 424); 

                auto tg_yyyz_yyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 425); 

                auto tg_yyyz_yyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 426); 

                auto tg_yyyz_yyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 427); 

                auto tg_yyyz_yyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 428); 

                auto tg_yyyz_yyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 429); 

                auto tg_yyyz_yzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 430); 

                auto tg_yyyz_zzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 431); 

                auto tg_yyzz_xxxxxxx_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 432); 

                auto tg_yyzz_xxxxxxy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 433); 

                auto tg_yyzz_xxxxxxz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 434); 

                auto tg_yyzz_xxxxxyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 435); 

                auto tg_yyzz_xxxxxyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 436); 

                auto tg_yyzz_xxxxxzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 437); 

                auto tg_yyzz_xxxxyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 438); 

                auto tg_yyzz_xxxxyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 439); 

                auto tg_yyzz_xxxxyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 440); 

                auto tg_yyzz_xxxxzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 441); 

                auto tg_yyzz_xxxyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 442); 

                auto tg_yyzz_xxxyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 443); 

                auto tg_yyzz_xxxyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 444); 

                auto tg_yyzz_xxxyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 445); 

                auto tg_yyzz_xxxzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 446); 

                auto tg_yyzz_xxyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 447); 

                auto tg_yyzz_xxyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 448); 

                auto tg_yyzz_xxyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 449); 

                auto tg_yyzz_xxyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 450); 

                auto tg_yyzz_xxyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 451); 

                auto tg_yyzz_xxzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 452); 

                auto tg_yyzz_xyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 453); 

                auto tg_yyzz_xyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 454); 

                auto tg_yyzz_xyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 455); 

                auto tg_yyzz_xyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 456); 

                auto tg_yyzz_xyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 457); 

                auto tg_yyzz_xyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 458); 

                auto tg_yyzz_xzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 459); 

                auto tg_yyzz_yyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 460); 

                auto tg_yyzz_yyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 461); 

                auto tg_yyzz_yyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 462); 

                auto tg_yyzz_yyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 463); 

                auto tg_yyzz_yyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 464); 

                auto tg_yyzz_yyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 465); 

                auto tg_yyzz_yzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 466); 

                auto tg_yyzz_zzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 467); 

                auto tg_yzzz_xxxxxxx_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 468); 

                auto tg_yzzz_xxxxxxy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 469); 

                auto tg_yzzz_xxxxxxz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 470); 

                auto tg_yzzz_xxxxxyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 471); 

                auto tg_yzzz_xxxxxyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 472); 

                auto tg_yzzz_xxxxxzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 473); 

                auto tg_yyyy_xzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 300); 

                auto tg_yyyy_yyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 301); 

                auto tg_yyyy_yyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 302); 

                auto tg_yyyy_yyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 303); 

                auto tg_yyyy_yyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 304); 

                auto tg_yyyy_yyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 305); 

                auto tg_yyyy_yzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 306); 

                auto tg_yyyy_zzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 307); 

                auto tg_yyyz_xxxxxx_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 308); 

                auto tg_yyyz_xxxxxy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 309); 

                auto tg_yyyz_xxxxxz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 310); 

                auto tg_yyyz_xxxxyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 311); 

                auto tg_yyyz_xxxxyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 312); 

                auto tg_yyyz_xxxxzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 313); 

                auto tg_yyyz_xxxyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 314); 

                auto tg_yyyz_xxxyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 315); 

                auto tg_yyyz_xxxyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 316); 

                auto tg_yyyz_xxxzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 317); 

                auto tg_yyyz_xxyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 318); 

                auto tg_yyyz_xxyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 319); 

                auto tg_yyyz_xxyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 320); 

                auto tg_yyyz_xxyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 321); 

                auto tg_yyyz_xxzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 322); 

                auto tg_yyyz_xyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 323); 

                auto tg_yyyz_xyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 324); 

                auto tg_yyyz_xyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 325); 

                auto tg_yyyz_xyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 326); 

                auto tg_yyyz_xyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 327); 

                auto tg_yyyz_xzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 328); 

                auto tg_yyyz_yyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 329); 

                auto tg_yyyz_yyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 330); 

                auto tg_yyyz_yyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 331); 

                auto tg_yyyz_yyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 332); 

                auto tg_yyyz_yyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 333); 

                auto tg_yyyz_yzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 334); 

                auto tg_yyyz_zzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 335); 

                auto tg_yyzz_xxxxxx_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 336); 

                auto tg_yyzz_xxxxxy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 337); 

                auto tg_yyzz_xxxxxz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 338); 

                auto tg_yyzz_xxxxyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 339); 

                auto tg_yyzz_xxxxyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 340); 

                auto tg_yyzz_xxxxzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 341); 

                auto tg_yyzz_xxxyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 342); 

                auto tg_yyzz_xxxyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 343); 

                auto tg_yyzz_xxxyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 344); 

                auto tg_yyzz_xxxzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 345); 

                auto tg_yyzz_xxyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 346); 

                auto tg_yyzz_xxyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 347); 

                auto tg_yyzz_xxyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 348); 

                auto tg_yyzz_xxyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 349); 

                auto tg_yyzz_xxzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 350); 

                auto tg_yyzz_xyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 351); 

                auto tg_yyzz_xyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 352); 

                auto tg_yyzz_xyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 353); 

                auto tg_yyzz_xyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 354); 

                auto tg_yyzz_xyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 355); 

                auto tg_yyzz_xzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 356); 

                auto tg_yyzz_yyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 357); 

                auto tg_yyzz_yyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 358); 

                auto tg_yyzz_yyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 359); 

                auto tg_yyzz_yyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 360); 

                auto tg_yyzz_yyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 361); 

                auto tg_yyzz_yzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 362); 

                auto tg_yyzz_zzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 363); 

                auto tg_yzzz_xxxxxx_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 364); 

                auto tg_yzzz_xxxxxy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 365); 

                auto tg_yzzz_xxxxxz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 366); 

                auto tg_yzzz_xxxxyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 367); 

                auto tg_yzzz_xxxxyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 368); 

                auto tg_yzzz_xxxxzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 369); 

                // set up pointers to integrals

                auto tg_xyyyy_xxzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 380); 

                auto tg_xyyyy_xyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 381); 

                auto tg_xyyyy_xyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 382); 

                auto tg_xyyyy_xyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 383); 

                auto tg_xyyyy_xyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 384); 

                auto tg_xyyyy_xyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 385); 

                auto tg_xyyyy_xyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 386); 

                auto tg_xyyyy_xzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 387); 

                auto tg_xyyyy_yyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 388); 

                auto tg_xyyyy_yyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 389); 

                auto tg_xyyyy_yyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 390); 

                auto tg_xyyyy_yyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 391); 

                auto tg_xyyyy_yyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 392); 

                auto tg_xyyyy_yyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 393); 

                auto tg_xyyyy_yzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 394); 

                auto tg_xyyyy_zzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 395); 

                auto tg_xyyyz_xxxxxxx_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 396); 

                auto tg_xyyyz_xxxxxxy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 397); 

                auto tg_xyyyz_xxxxxxz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 398); 

                auto tg_xyyyz_xxxxxyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 399); 

                auto tg_xyyyz_xxxxxyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 400); 

                auto tg_xyyyz_xxxxxzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 401); 

                auto tg_xyyyz_xxxxyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 402); 

                auto tg_xyyyz_xxxxyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 403); 

                auto tg_xyyyz_xxxxyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 404); 

                auto tg_xyyyz_xxxxzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 405); 

                auto tg_xyyyz_xxxyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 406); 

                auto tg_xyyyz_xxxyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 407); 

                auto tg_xyyyz_xxxyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 408); 

                auto tg_xyyyz_xxxyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 409); 

                auto tg_xyyyz_xxxzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 410); 

                auto tg_xyyyz_xxyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 411); 

                auto tg_xyyyz_xxyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 412); 

                auto tg_xyyyz_xxyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 413); 

                auto tg_xyyyz_xxyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 414); 

                auto tg_xyyyz_xxyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 415); 

                auto tg_xyyyz_xxzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 416); 

                auto tg_xyyyz_xyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 417); 

                auto tg_xyyyz_xyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 418); 

                auto tg_xyyyz_xyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 419); 

                auto tg_xyyyz_xyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 420); 

                auto tg_xyyyz_xyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 421); 

                auto tg_xyyyz_xyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 422); 

                auto tg_xyyyz_xzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 423); 

                auto tg_xyyyz_yyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 424); 

                auto tg_xyyyz_yyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 425); 

                auto tg_xyyyz_yyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 426); 

                auto tg_xyyyz_yyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 427); 

                auto tg_xyyyz_yyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 428); 

                auto tg_xyyyz_yyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 429); 

                auto tg_xyyyz_yzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 430); 

                auto tg_xyyyz_zzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 431); 

                auto tg_xyyzz_xxxxxxx_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 432); 

                auto tg_xyyzz_xxxxxxy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 433); 

                auto tg_xyyzz_xxxxxxz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 434); 

                auto tg_xyyzz_xxxxxyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 435); 

                auto tg_xyyzz_xxxxxyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 436); 

                auto tg_xyyzz_xxxxxzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 437); 

                auto tg_xyyzz_xxxxyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 438); 

                auto tg_xyyzz_xxxxyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 439); 

                auto tg_xyyzz_xxxxyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 440); 

                auto tg_xyyzz_xxxxzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 441); 

                auto tg_xyyzz_xxxyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 442); 

                auto tg_xyyzz_xxxyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 443); 

                auto tg_xyyzz_xxxyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 444); 

                auto tg_xyyzz_xxxyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 445); 

                auto tg_xyyzz_xxxzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 446); 

                auto tg_xyyzz_xxyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 447); 

                auto tg_xyyzz_xxyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 448); 

                auto tg_xyyzz_xxyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 449); 

                auto tg_xyyzz_xxyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 450); 

                auto tg_xyyzz_xxyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 451); 

                auto tg_xyyzz_xxzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 452); 

                auto tg_xyyzz_xyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 453); 

                auto tg_xyyzz_xyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 454); 

                auto tg_xyyzz_xyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 455); 

                auto tg_xyyzz_xyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 456); 

                auto tg_xyyzz_xyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 457); 

                auto tg_xyyzz_xyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 458); 

                auto tg_xyyzz_xzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 459); 

                auto tg_xyyzz_yyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 460); 

                auto tg_xyyzz_yyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 461); 

                auto tg_xyyzz_yyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 462); 

                auto tg_xyyzz_yyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 463); 

                auto tg_xyyzz_yyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 464); 

                auto tg_xyyzz_yyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 465); 

                auto tg_xyyzz_yzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 466); 

                auto tg_xyyzz_zzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 467); 

                auto tg_xyzzz_xxxxxxx_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 468); 

                auto tg_xyzzz_xxxxxxy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 469); 

                auto tg_xyzzz_xxxxxxz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 470); 

                auto tg_xyzzz_xxxxxyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 471); 

                auto tg_xyzzz_xxxxxyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 472); 

                auto tg_xyzzz_xxxxxzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 473); 

                // Batch of Integrals (380,474)

                #pragma omp simd aligned(fxn, tg_xyyyy_xxzzzzz_0, tg_xyyyy_xyyyyyy_0, tg_xyyyy_xyyyyyz_0, \
                                         tg_xyyyy_xyyyyzz_0, tg_xyyyy_xyyyzzz_0, tg_xyyyy_xyyzzzz_0, tg_xyyyy_xyzzzzz_0, \
                                         tg_xyyyy_xzzzzzz_0, tg_xyyyy_yyyyyyy_0, tg_xyyyy_yyyyyyz_0, tg_xyyyy_yyyyyzz_0, \
                                         tg_xyyyy_yyyyzzz_0, tg_xyyyy_yyyzzzz_0, tg_xyyyy_yyzzzzz_0, tg_xyyyy_yzzzzzz_0, \
                                         tg_xyyyy_zzzzzzz_0, tg_xyyyz_xxxxxxx_0, tg_xyyyz_xxxxxxy_0, tg_xyyyz_xxxxxxz_0, \
                                         tg_xyyyz_xxxxxyy_0, tg_xyyyz_xxxxxyz_0, tg_xyyyz_xxxxxzz_0, tg_xyyyz_xxxxyyy_0, \
                                         tg_xyyyz_xxxxyyz_0, tg_xyyyz_xxxxyzz_0, tg_xyyyz_xxxxzzz_0, tg_xyyyz_xxxyyyy_0, \
                                         tg_xyyyz_xxxyyyz_0, tg_xyyyz_xxxyyzz_0, tg_xyyyz_xxxyzzz_0, tg_xyyyz_xxxzzzz_0, \
                                         tg_xyyyz_xxyyyyy_0, tg_xyyyz_xxyyyyz_0, tg_xyyyz_xxyyyzz_0, tg_xyyyz_xxyyzzz_0, \
                                         tg_xyyyz_xxyzzzz_0, tg_xyyyz_xxzzzzz_0, tg_xyyyz_xyyyyyy_0, tg_xyyyz_xyyyyyz_0, \
                                         tg_xyyyz_xyyyyzz_0, tg_xyyyz_xyyyzzz_0, tg_xyyyz_xyyzzzz_0, tg_xyyyz_xyzzzzz_0, \
                                         tg_xyyyz_xzzzzzz_0, tg_xyyyz_yyyyyyy_0, tg_xyyyz_yyyyyyz_0, tg_xyyyz_yyyyyzz_0, \
                                         tg_xyyyz_yyyyzzz_0, tg_xyyyz_yyyzzzz_0, tg_xyyyz_yyzzzzz_0, tg_xyyyz_yzzzzzz_0, \
                                         tg_xyyyz_zzzzzzz_0, tg_xyyzz_xxxxxxx_0, tg_xyyzz_xxxxxxy_0, tg_xyyzz_xxxxxxz_0, \
                                         tg_xyyzz_xxxxxyy_0, tg_xyyzz_xxxxxyz_0, tg_xyyzz_xxxxxzz_0, tg_xyyzz_xxxxyyy_0, \
                                         tg_xyyzz_xxxxyyz_0, tg_xyyzz_xxxxyzz_0, tg_xyyzz_xxxxzzz_0, tg_xyyzz_xxxyyyy_0, \
                                         tg_xyyzz_xxxyyyz_0, tg_xyyzz_xxxyyzz_0, tg_xyyzz_xxxyzzz_0, tg_xyyzz_xxxzzzz_0, \
                                         tg_xyyzz_xxyyyyy_0, tg_xyyzz_xxyyyyz_0, tg_xyyzz_xxyyyzz_0, tg_xyyzz_xxyyzzz_0, \
                                         tg_xyyzz_xxyzzzz_0, tg_xyyzz_xxzzzzz_0, tg_xyyzz_xyyyyyy_0, tg_xyyzz_xyyyyyz_0, \
                                         tg_xyyzz_xyyyyzz_0, tg_xyyzz_xyyyzzz_0, tg_xyyzz_xyyzzzz_0, tg_xyyzz_xyzzzzz_0, \
                                         tg_xyyzz_xzzzzzz_0, tg_xyyzz_yyyyyyy_0, tg_xyyzz_yyyyyyz_0, tg_xyyzz_yyyyyzz_0, \
                                         tg_xyyzz_yyyyzzz_0, tg_xyyzz_yyyzzzz_0, tg_xyyzz_yyzzzzz_0, tg_xyyzz_yzzzzzz_0, \
                                         tg_xyyzz_zzzzzzz_0, tg_xyzzz_xxxxxxx_0, tg_xyzzz_xxxxxxy_0, tg_xyzzz_xxxxxxz_0, \
                                         tg_xyzzz_xxxxxyy_0, tg_xyzzz_xxxxxyz_0, tg_xyzzz_xxxxxzz_0, tg_yyyy_xxzzzzz_0, \
                                         tg_yyyy_xxzzzzz_1, tg_yyyy_xyyyyyy_0, tg_yyyy_xyyyyyy_1, tg_yyyy_xyyyyyz_0, \
                                         tg_yyyy_xyyyyyz_1, tg_yyyy_xyyyyzz_0, tg_yyyy_xyyyyzz_1, tg_yyyy_xyyyzzz_0, \
                                         tg_yyyy_xyyyzzz_1, tg_yyyy_xyyzzzz_0, tg_yyyy_xyyzzzz_1, tg_yyyy_xyzzzzz_0, \
                                         tg_yyyy_xyzzzzz_1, tg_yyyy_xzzzzz_1, tg_yyyy_xzzzzzz_0, tg_yyyy_xzzzzzz_1, \
                                         tg_yyyy_yyyyyy_1, tg_yyyy_yyyyyyy_0, tg_yyyy_yyyyyyy_1, tg_yyyy_yyyyyyz_0, \
                                         tg_yyyy_yyyyyyz_1, tg_yyyy_yyyyyz_1, tg_yyyy_yyyyyzz_0, tg_yyyy_yyyyyzz_1, \
                                         tg_yyyy_yyyyzz_1, tg_yyyy_yyyyzzz_0, tg_yyyy_yyyyzzz_1, tg_yyyy_yyyzzz_1, \
                                         tg_yyyy_yyyzzzz_0, tg_yyyy_yyyzzzz_1, tg_yyyy_yyzzzz_1, tg_yyyy_yyzzzzz_0, \
                                         tg_yyyy_yyzzzzz_1, tg_yyyy_yzzzzz_1, tg_yyyy_yzzzzzz_0, tg_yyyy_yzzzzzz_1, \
                                         tg_yyyy_zzzzzz_1, tg_yyyy_zzzzzzz_0, tg_yyyy_zzzzzzz_1, tg_yyyz_xxxxxx_1, \
                                         tg_yyyz_xxxxxxx_0, tg_yyyz_xxxxxxx_1, tg_yyyz_xxxxxxy_0, tg_yyyz_xxxxxxy_1, \
                                         tg_yyyz_xxxxxxz_0, tg_yyyz_xxxxxxz_1, tg_yyyz_xxxxxy_1, tg_yyyz_xxxxxyy_0, \
                                         tg_yyyz_xxxxxyy_1, tg_yyyz_xxxxxyz_0, tg_yyyz_xxxxxyz_1, tg_yyyz_xxxxxz_1, \
                                         tg_yyyz_xxxxxzz_0, tg_yyyz_xxxxxzz_1, tg_yyyz_xxxxyy_1, tg_yyyz_xxxxyyy_0, \
                                         tg_yyyz_xxxxyyy_1, tg_yyyz_xxxxyyz_0, tg_yyyz_xxxxyyz_1, tg_yyyz_xxxxyz_1, \
                                         tg_yyyz_xxxxyzz_0, tg_yyyz_xxxxyzz_1, tg_yyyz_xxxxzz_1, tg_yyyz_xxxxzzz_0, \
                                         tg_yyyz_xxxxzzz_1, tg_yyyz_xxxyyy_1, tg_yyyz_xxxyyyy_0, tg_yyyz_xxxyyyy_1, \
                                         tg_yyyz_xxxyyyz_0, tg_yyyz_xxxyyyz_1, tg_yyyz_xxxyyz_1, tg_yyyz_xxxyyzz_0, \
                                         tg_yyyz_xxxyyzz_1, tg_yyyz_xxxyzz_1, tg_yyyz_xxxyzzz_0, tg_yyyz_xxxyzzz_1, \
                                         tg_yyyz_xxxzzz_1, tg_yyyz_xxxzzzz_0, tg_yyyz_xxxzzzz_1, tg_yyyz_xxyyyy_1, \
                                         tg_yyyz_xxyyyyy_0, tg_yyyz_xxyyyyy_1, tg_yyyz_xxyyyyz_0, tg_yyyz_xxyyyyz_1, \
                                         tg_yyyz_xxyyyz_1, tg_yyyz_xxyyyzz_0, tg_yyyz_xxyyyzz_1, tg_yyyz_xxyyzz_1, \
                                         tg_yyyz_xxyyzzz_0, tg_yyyz_xxyyzzz_1, tg_yyyz_xxyzzz_1, tg_yyyz_xxyzzzz_0, \
                                         tg_yyyz_xxyzzzz_1, tg_yyyz_xxzzzz_1, tg_yyyz_xxzzzzz_0, tg_yyyz_xxzzzzz_1, \
                                         tg_yyyz_xyyyyy_1, tg_yyyz_xyyyyyy_0, tg_yyyz_xyyyyyy_1, tg_yyyz_xyyyyyz_0, \
                                         tg_yyyz_xyyyyyz_1, tg_yyyz_xyyyyz_1, tg_yyyz_xyyyyzz_0, tg_yyyz_xyyyyzz_1, \
                                         tg_yyyz_xyyyzz_1, tg_yyyz_xyyyzzz_0, tg_yyyz_xyyyzzz_1, tg_yyyz_xyyzzz_1, \
                                         tg_yyyz_xyyzzzz_0, tg_yyyz_xyyzzzz_1, tg_yyyz_xyzzzz_1, tg_yyyz_xyzzzzz_0, \
                                         tg_yyyz_xyzzzzz_1, tg_yyyz_xzzzzz_1, tg_yyyz_xzzzzzz_0, tg_yyyz_xzzzzzz_1, \
                                         tg_yyyz_yyyyyy_1, tg_yyyz_yyyyyyy_0, tg_yyyz_yyyyyyy_1, tg_yyyz_yyyyyyz_0, \
                                         tg_yyyz_yyyyyyz_1, tg_yyyz_yyyyyz_1, tg_yyyz_yyyyyzz_0, tg_yyyz_yyyyyzz_1, \
                                         tg_yyyz_yyyyzz_1, tg_yyyz_yyyyzzz_0, tg_yyyz_yyyyzzz_1, tg_yyyz_yyyzzz_1, \
                                         tg_yyyz_yyyzzzz_0, tg_yyyz_yyyzzzz_1, tg_yyyz_yyzzzz_1, tg_yyyz_yyzzzzz_0, \
                                         tg_yyyz_yyzzzzz_1, tg_yyyz_yzzzzz_1, tg_yyyz_yzzzzzz_0, tg_yyyz_yzzzzzz_1, \
                                         tg_yyyz_zzzzzz_1, tg_yyyz_zzzzzzz_0, tg_yyyz_zzzzzzz_1, tg_yyzz_xxxxxx_1, \
                                         tg_yyzz_xxxxxxx_0, tg_yyzz_xxxxxxx_1, tg_yyzz_xxxxxxy_0, tg_yyzz_xxxxxxy_1, \
                                         tg_yyzz_xxxxxxz_0, tg_yyzz_xxxxxxz_1, tg_yyzz_xxxxxy_1, tg_yyzz_xxxxxyy_0, \
                                         tg_yyzz_xxxxxyy_1, tg_yyzz_xxxxxyz_0, tg_yyzz_xxxxxyz_1, tg_yyzz_xxxxxz_1, \
                                         tg_yyzz_xxxxxzz_0, tg_yyzz_xxxxxzz_1, tg_yyzz_xxxxyy_1, tg_yyzz_xxxxyyy_0, \
                                         tg_yyzz_xxxxyyy_1, tg_yyzz_xxxxyyz_0, tg_yyzz_xxxxyyz_1, tg_yyzz_xxxxyz_1, \
                                         tg_yyzz_xxxxyzz_0, tg_yyzz_xxxxyzz_1, tg_yyzz_xxxxzz_1, tg_yyzz_xxxxzzz_0, \
                                         tg_yyzz_xxxxzzz_1, tg_yyzz_xxxyyy_1, tg_yyzz_xxxyyyy_0, tg_yyzz_xxxyyyy_1, \
                                         tg_yyzz_xxxyyyz_0, tg_yyzz_xxxyyyz_1, tg_yyzz_xxxyyz_1, tg_yyzz_xxxyyzz_0, \
                                         tg_yyzz_xxxyyzz_1, tg_yyzz_xxxyzz_1, tg_yyzz_xxxyzzz_0, tg_yyzz_xxxyzzz_1, \
                                         tg_yyzz_xxxzzz_1, tg_yyzz_xxxzzzz_0, tg_yyzz_xxxzzzz_1, tg_yyzz_xxyyyy_1, \
                                         tg_yyzz_xxyyyyy_0, tg_yyzz_xxyyyyy_1, tg_yyzz_xxyyyyz_0, tg_yyzz_xxyyyyz_1, \
                                         tg_yyzz_xxyyyz_1, tg_yyzz_xxyyyzz_0, tg_yyzz_xxyyyzz_1, tg_yyzz_xxyyzz_1, \
                                         tg_yyzz_xxyyzzz_0, tg_yyzz_xxyyzzz_1, tg_yyzz_xxyzzz_1, tg_yyzz_xxyzzzz_0, \
                                         tg_yyzz_xxyzzzz_1, tg_yyzz_xxzzzz_1, tg_yyzz_xxzzzzz_0, tg_yyzz_xxzzzzz_1, \
                                         tg_yyzz_xyyyyy_1, tg_yyzz_xyyyyyy_0, tg_yyzz_xyyyyyy_1, tg_yyzz_xyyyyyz_0, \
                                         tg_yyzz_xyyyyyz_1, tg_yyzz_xyyyyz_1, tg_yyzz_xyyyyzz_0, tg_yyzz_xyyyyzz_1, \
                                         tg_yyzz_xyyyzz_1, tg_yyzz_xyyyzzz_0, tg_yyzz_xyyyzzz_1, tg_yyzz_xyyzzz_1, \
                                         tg_yyzz_xyyzzzz_0, tg_yyzz_xyyzzzz_1, tg_yyzz_xyzzzz_1, tg_yyzz_xyzzzzz_0, \
                                         tg_yyzz_xyzzzzz_1, tg_yyzz_xzzzzz_1, tg_yyzz_xzzzzzz_0, tg_yyzz_xzzzzzz_1, \
                                         tg_yyzz_yyyyyy_1, tg_yyzz_yyyyyyy_0, tg_yyzz_yyyyyyy_1, tg_yyzz_yyyyyyz_0, \
                                         tg_yyzz_yyyyyyz_1, tg_yyzz_yyyyyz_1, tg_yyzz_yyyyyzz_0, tg_yyzz_yyyyyzz_1, \
                                         tg_yyzz_yyyyzz_1, tg_yyzz_yyyyzzz_0, tg_yyzz_yyyyzzz_1, tg_yyzz_yyyzzz_1, \
                                         tg_yyzz_yyyzzzz_0, tg_yyzz_yyyzzzz_1, tg_yyzz_yyzzzz_1, tg_yyzz_yyzzzzz_0, \
                                         tg_yyzz_yyzzzzz_1, tg_yyzz_yzzzzz_1, tg_yyzz_yzzzzzz_0, tg_yyzz_yzzzzzz_1, \
                                         tg_yyzz_zzzzzz_1, tg_yyzz_zzzzzzz_0, tg_yyzz_zzzzzzz_1, tg_yzzz_xxxxxx_1, \
                                         tg_yzzz_xxxxxxx_0, tg_yzzz_xxxxxxx_1, tg_yzzz_xxxxxxy_0, tg_yzzz_xxxxxxy_1, \
                                         tg_yzzz_xxxxxxz_0, tg_yzzz_xxxxxxz_1, tg_yzzz_xxxxxy_1, tg_yzzz_xxxxxyy_0, \
                                         tg_yzzz_xxxxxyy_1, tg_yzzz_xxxxxyz_0, tg_yzzz_xxxxxyz_1, tg_yzzz_xxxxxz_1, \
                                         tg_yzzz_xxxxxzz_0, tg_yzzz_xxxxxzz_1, tg_yzzz_xxxxyy_1, tg_yzzz_xxxxyz_1, \
                                         tg_yzzz_xxxxzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fxn = fxn[j];

                    tg_xyyyy_xxzzzzz_0[j] = pb_x * tg_yyyy_xxzzzzz_0[j] + wp_x[j] * tg_yyyy_xxzzzzz_1[j] + fl1_fxn * tg_yyyy_xzzzzz_1[j];

                    tg_xyyyy_xyyyyyy_0[j] = pb_x * tg_yyyy_xyyyyyy_0[j] + wp_x[j] * tg_yyyy_xyyyyyy_1[j] + 0.5 * fl1_fxn * tg_yyyy_yyyyyy_1[j];

                    tg_xyyyy_xyyyyyz_0[j] = pb_x * tg_yyyy_xyyyyyz_0[j] + wp_x[j] * tg_yyyy_xyyyyyz_1[j] + 0.5 * fl1_fxn * tg_yyyy_yyyyyz_1[j];

                    tg_xyyyy_xyyyyzz_0[j] = pb_x * tg_yyyy_xyyyyzz_0[j] + wp_x[j] * tg_yyyy_xyyyyzz_1[j] + 0.5 * fl1_fxn * tg_yyyy_yyyyzz_1[j];

                    tg_xyyyy_xyyyzzz_0[j] = pb_x * tg_yyyy_xyyyzzz_0[j] + wp_x[j] * tg_yyyy_xyyyzzz_1[j] + 0.5 * fl1_fxn * tg_yyyy_yyyzzz_1[j];

                    tg_xyyyy_xyyzzzz_0[j] = pb_x * tg_yyyy_xyyzzzz_0[j] + wp_x[j] * tg_yyyy_xyyzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyy_yyzzzz_1[j];

                    tg_xyyyy_xyzzzzz_0[j] = pb_x * tg_yyyy_xyzzzzz_0[j] + wp_x[j] * tg_yyyy_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyy_yzzzzz_1[j];

                    tg_xyyyy_xzzzzzz_0[j] = pb_x * tg_yyyy_xzzzzzz_0[j] + wp_x[j] * tg_yyyy_xzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyy_zzzzzz_1[j];

                    tg_xyyyy_yyyyyyy_0[j] = pb_x * tg_yyyy_yyyyyyy_0[j] + wp_x[j] * tg_yyyy_yyyyyyy_1[j];

                    tg_xyyyy_yyyyyyz_0[j] = pb_x * tg_yyyy_yyyyyyz_0[j] + wp_x[j] * tg_yyyy_yyyyyyz_1[j];

                    tg_xyyyy_yyyyyzz_0[j] = pb_x * tg_yyyy_yyyyyzz_0[j] + wp_x[j] * tg_yyyy_yyyyyzz_1[j];

                    tg_xyyyy_yyyyzzz_0[j] = pb_x * tg_yyyy_yyyyzzz_0[j] + wp_x[j] * tg_yyyy_yyyyzzz_1[j];

                    tg_xyyyy_yyyzzzz_0[j] = pb_x * tg_yyyy_yyyzzzz_0[j] + wp_x[j] * tg_yyyy_yyyzzzz_1[j];

                    tg_xyyyy_yyzzzzz_0[j] = pb_x * tg_yyyy_yyzzzzz_0[j] + wp_x[j] * tg_yyyy_yyzzzzz_1[j];

                    tg_xyyyy_yzzzzzz_0[j] = pb_x * tg_yyyy_yzzzzzz_0[j] + wp_x[j] * tg_yyyy_yzzzzzz_1[j];

                    tg_xyyyy_zzzzzzz_0[j] = pb_x * tg_yyyy_zzzzzzz_0[j] + wp_x[j] * tg_yyyy_zzzzzzz_1[j];

                    tg_xyyyz_xxxxxxx_0[j] = pb_x * tg_yyyz_xxxxxxx_0[j] + wp_x[j] * tg_yyyz_xxxxxxx_1[j] + 3.5 * fl1_fxn * tg_yyyz_xxxxxx_1[j];

                    tg_xyyyz_xxxxxxy_0[j] = pb_x * tg_yyyz_xxxxxxy_0[j] + wp_x[j] * tg_yyyz_xxxxxxy_1[j] + 3.0 * fl1_fxn * tg_yyyz_xxxxxy_1[j];

                    tg_xyyyz_xxxxxxz_0[j] = pb_x * tg_yyyz_xxxxxxz_0[j] + wp_x[j] * tg_yyyz_xxxxxxz_1[j] + 3.0 * fl1_fxn * tg_yyyz_xxxxxz_1[j];

                    tg_xyyyz_xxxxxyy_0[j] = pb_x * tg_yyyz_xxxxxyy_0[j] + wp_x[j] * tg_yyyz_xxxxxyy_1[j] + 2.5 * fl1_fxn * tg_yyyz_xxxxyy_1[j];

                    tg_xyyyz_xxxxxyz_0[j] = pb_x * tg_yyyz_xxxxxyz_0[j] + wp_x[j] * tg_yyyz_xxxxxyz_1[j] + 2.5 * fl1_fxn * tg_yyyz_xxxxyz_1[j];

                    tg_xyyyz_xxxxxzz_0[j] = pb_x * tg_yyyz_xxxxxzz_0[j] + wp_x[j] * tg_yyyz_xxxxxzz_1[j] + 2.5 * fl1_fxn * tg_yyyz_xxxxzz_1[j];

                    tg_xyyyz_xxxxyyy_0[j] = pb_x * tg_yyyz_xxxxyyy_0[j] + wp_x[j] * tg_yyyz_xxxxyyy_1[j] + 2.0 * fl1_fxn * tg_yyyz_xxxyyy_1[j];

                    tg_xyyyz_xxxxyyz_0[j] = pb_x * tg_yyyz_xxxxyyz_0[j] + wp_x[j] * tg_yyyz_xxxxyyz_1[j] + 2.0 * fl1_fxn * tg_yyyz_xxxyyz_1[j];

                    tg_xyyyz_xxxxyzz_0[j] = pb_x * tg_yyyz_xxxxyzz_0[j] + wp_x[j] * tg_yyyz_xxxxyzz_1[j] + 2.0 * fl1_fxn * tg_yyyz_xxxyzz_1[j];

                    tg_xyyyz_xxxxzzz_0[j] = pb_x * tg_yyyz_xxxxzzz_0[j] + wp_x[j] * tg_yyyz_xxxxzzz_1[j] + 2.0 * fl1_fxn * tg_yyyz_xxxzzz_1[j];

                    tg_xyyyz_xxxyyyy_0[j] = pb_x * tg_yyyz_xxxyyyy_0[j] + wp_x[j] * tg_yyyz_xxxyyyy_1[j] + 1.5 * fl1_fxn * tg_yyyz_xxyyyy_1[j];

                    tg_xyyyz_xxxyyyz_0[j] = pb_x * tg_yyyz_xxxyyyz_0[j] + wp_x[j] * tg_yyyz_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_yyyz_xxyyyz_1[j];

                    tg_xyyyz_xxxyyzz_0[j] = pb_x * tg_yyyz_xxxyyzz_0[j] + wp_x[j] * tg_yyyz_xxxyyzz_1[j] + 1.5 * fl1_fxn * tg_yyyz_xxyyzz_1[j];

                    tg_xyyyz_xxxyzzz_0[j] = pb_x * tg_yyyz_xxxyzzz_0[j] + wp_x[j] * tg_yyyz_xxxyzzz_1[j] + 1.5 * fl1_fxn * tg_yyyz_xxyzzz_1[j];

                    tg_xyyyz_xxxzzzz_0[j] = pb_x * tg_yyyz_xxxzzzz_0[j] + wp_x[j] * tg_yyyz_xxxzzzz_1[j] + 1.5 * fl1_fxn * tg_yyyz_xxzzzz_1[j];

                    tg_xyyyz_xxyyyyy_0[j] = pb_x * tg_yyyz_xxyyyyy_0[j] + wp_x[j] * tg_yyyz_xxyyyyy_1[j] + fl1_fxn * tg_yyyz_xyyyyy_1[j];

                    tg_xyyyz_xxyyyyz_0[j] = pb_x * tg_yyyz_xxyyyyz_0[j] + wp_x[j] * tg_yyyz_xxyyyyz_1[j] + fl1_fxn * tg_yyyz_xyyyyz_1[j];

                    tg_xyyyz_xxyyyzz_0[j] = pb_x * tg_yyyz_xxyyyzz_0[j] + wp_x[j] * tg_yyyz_xxyyyzz_1[j] + fl1_fxn * tg_yyyz_xyyyzz_1[j];

                    tg_xyyyz_xxyyzzz_0[j] = pb_x * tg_yyyz_xxyyzzz_0[j] + wp_x[j] * tg_yyyz_xxyyzzz_1[j] + fl1_fxn * tg_yyyz_xyyzzz_1[j];

                    tg_xyyyz_xxyzzzz_0[j] = pb_x * tg_yyyz_xxyzzzz_0[j] + wp_x[j] * tg_yyyz_xxyzzzz_1[j] + fl1_fxn * tg_yyyz_xyzzzz_1[j];

                    tg_xyyyz_xxzzzzz_0[j] = pb_x * tg_yyyz_xxzzzzz_0[j] + wp_x[j] * tg_yyyz_xxzzzzz_1[j] + fl1_fxn * tg_yyyz_xzzzzz_1[j];

                    tg_xyyyz_xyyyyyy_0[j] = pb_x * tg_yyyz_xyyyyyy_0[j] + wp_x[j] * tg_yyyz_xyyyyyy_1[j] + 0.5 * fl1_fxn * tg_yyyz_yyyyyy_1[j];

                    tg_xyyyz_xyyyyyz_0[j] = pb_x * tg_yyyz_xyyyyyz_0[j] + wp_x[j] * tg_yyyz_xyyyyyz_1[j] + 0.5 * fl1_fxn * tg_yyyz_yyyyyz_1[j];

                    tg_xyyyz_xyyyyzz_0[j] = pb_x * tg_yyyz_xyyyyzz_0[j] + wp_x[j] * tg_yyyz_xyyyyzz_1[j] + 0.5 * fl1_fxn * tg_yyyz_yyyyzz_1[j];

                    tg_xyyyz_xyyyzzz_0[j] = pb_x * tg_yyyz_xyyyzzz_0[j] + wp_x[j] * tg_yyyz_xyyyzzz_1[j] + 0.5 * fl1_fxn * tg_yyyz_yyyzzz_1[j];

                    tg_xyyyz_xyyzzzz_0[j] = pb_x * tg_yyyz_xyyzzzz_0[j] + wp_x[j] * tg_yyyz_xyyzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyz_yyzzzz_1[j];

                    tg_xyyyz_xyzzzzz_0[j] = pb_x * tg_yyyz_xyzzzzz_0[j] + wp_x[j] * tg_yyyz_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyz_yzzzzz_1[j];

                    tg_xyyyz_xzzzzzz_0[j] = pb_x * tg_yyyz_xzzzzzz_0[j] + wp_x[j] * tg_yyyz_xzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyz_zzzzzz_1[j];

                    tg_xyyyz_yyyyyyy_0[j] = pb_x * tg_yyyz_yyyyyyy_0[j] + wp_x[j] * tg_yyyz_yyyyyyy_1[j];

                    tg_xyyyz_yyyyyyz_0[j] = pb_x * tg_yyyz_yyyyyyz_0[j] + wp_x[j] * tg_yyyz_yyyyyyz_1[j];

                    tg_xyyyz_yyyyyzz_0[j] = pb_x * tg_yyyz_yyyyyzz_0[j] + wp_x[j] * tg_yyyz_yyyyyzz_1[j];

                    tg_xyyyz_yyyyzzz_0[j] = pb_x * tg_yyyz_yyyyzzz_0[j] + wp_x[j] * tg_yyyz_yyyyzzz_1[j];

                    tg_xyyyz_yyyzzzz_0[j] = pb_x * tg_yyyz_yyyzzzz_0[j] + wp_x[j] * tg_yyyz_yyyzzzz_1[j];

                    tg_xyyyz_yyzzzzz_0[j] = pb_x * tg_yyyz_yyzzzzz_0[j] + wp_x[j] * tg_yyyz_yyzzzzz_1[j];

                    tg_xyyyz_yzzzzzz_0[j] = pb_x * tg_yyyz_yzzzzzz_0[j] + wp_x[j] * tg_yyyz_yzzzzzz_1[j];

                    tg_xyyyz_zzzzzzz_0[j] = pb_x * tg_yyyz_zzzzzzz_0[j] + wp_x[j] * tg_yyyz_zzzzzzz_1[j];

                    tg_xyyzz_xxxxxxx_0[j] = pb_x * tg_yyzz_xxxxxxx_0[j] + wp_x[j] * tg_yyzz_xxxxxxx_1[j] + 3.5 * fl1_fxn * tg_yyzz_xxxxxx_1[j];

                    tg_xyyzz_xxxxxxy_0[j] = pb_x * tg_yyzz_xxxxxxy_0[j] + wp_x[j] * tg_yyzz_xxxxxxy_1[j] + 3.0 * fl1_fxn * tg_yyzz_xxxxxy_1[j];

                    tg_xyyzz_xxxxxxz_0[j] = pb_x * tg_yyzz_xxxxxxz_0[j] + wp_x[j] * tg_yyzz_xxxxxxz_1[j] + 3.0 * fl1_fxn * tg_yyzz_xxxxxz_1[j];

                    tg_xyyzz_xxxxxyy_0[j] = pb_x * tg_yyzz_xxxxxyy_0[j] + wp_x[j] * tg_yyzz_xxxxxyy_1[j] + 2.5 * fl1_fxn * tg_yyzz_xxxxyy_1[j];

                    tg_xyyzz_xxxxxyz_0[j] = pb_x * tg_yyzz_xxxxxyz_0[j] + wp_x[j] * tg_yyzz_xxxxxyz_1[j] + 2.5 * fl1_fxn * tg_yyzz_xxxxyz_1[j];

                    tg_xyyzz_xxxxxzz_0[j] = pb_x * tg_yyzz_xxxxxzz_0[j] + wp_x[j] * tg_yyzz_xxxxxzz_1[j] + 2.5 * fl1_fxn * tg_yyzz_xxxxzz_1[j];

                    tg_xyyzz_xxxxyyy_0[j] = pb_x * tg_yyzz_xxxxyyy_0[j] + wp_x[j] * tg_yyzz_xxxxyyy_1[j] + 2.0 * fl1_fxn * tg_yyzz_xxxyyy_1[j];

                    tg_xyyzz_xxxxyyz_0[j] = pb_x * tg_yyzz_xxxxyyz_0[j] + wp_x[j] * tg_yyzz_xxxxyyz_1[j] + 2.0 * fl1_fxn * tg_yyzz_xxxyyz_1[j];

                    tg_xyyzz_xxxxyzz_0[j] = pb_x * tg_yyzz_xxxxyzz_0[j] + wp_x[j] * tg_yyzz_xxxxyzz_1[j] + 2.0 * fl1_fxn * tg_yyzz_xxxyzz_1[j];

                    tg_xyyzz_xxxxzzz_0[j] = pb_x * tg_yyzz_xxxxzzz_0[j] + wp_x[j] * tg_yyzz_xxxxzzz_1[j] + 2.0 * fl1_fxn * tg_yyzz_xxxzzz_1[j];

                    tg_xyyzz_xxxyyyy_0[j] = pb_x * tg_yyzz_xxxyyyy_0[j] + wp_x[j] * tg_yyzz_xxxyyyy_1[j] + 1.5 * fl1_fxn * tg_yyzz_xxyyyy_1[j];

                    tg_xyyzz_xxxyyyz_0[j] = pb_x * tg_yyzz_xxxyyyz_0[j] + wp_x[j] * tg_yyzz_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_yyzz_xxyyyz_1[j];

                    tg_xyyzz_xxxyyzz_0[j] = pb_x * tg_yyzz_xxxyyzz_0[j] + wp_x[j] * tg_yyzz_xxxyyzz_1[j] + 1.5 * fl1_fxn * tg_yyzz_xxyyzz_1[j];

                    tg_xyyzz_xxxyzzz_0[j] = pb_x * tg_yyzz_xxxyzzz_0[j] + wp_x[j] * tg_yyzz_xxxyzzz_1[j] + 1.5 * fl1_fxn * tg_yyzz_xxyzzz_1[j];

                    tg_xyyzz_xxxzzzz_0[j] = pb_x * tg_yyzz_xxxzzzz_0[j] + wp_x[j] * tg_yyzz_xxxzzzz_1[j] + 1.5 * fl1_fxn * tg_yyzz_xxzzzz_1[j];

                    tg_xyyzz_xxyyyyy_0[j] = pb_x * tg_yyzz_xxyyyyy_0[j] + wp_x[j] * tg_yyzz_xxyyyyy_1[j] + fl1_fxn * tg_yyzz_xyyyyy_1[j];

                    tg_xyyzz_xxyyyyz_0[j] = pb_x * tg_yyzz_xxyyyyz_0[j] + wp_x[j] * tg_yyzz_xxyyyyz_1[j] + fl1_fxn * tg_yyzz_xyyyyz_1[j];

                    tg_xyyzz_xxyyyzz_0[j] = pb_x * tg_yyzz_xxyyyzz_0[j] + wp_x[j] * tg_yyzz_xxyyyzz_1[j] + fl1_fxn * tg_yyzz_xyyyzz_1[j];

                    tg_xyyzz_xxyyzzz_0[j] = pb_x * tg_yyzz_xxyyzzz_0[j] + wp_x[j] * tg_yyzz_xxyyzzz_1[j] + fl1_fxn * tg_yyzz_xyyzzz_1[j];

                    tg_xyyzz_xxyzzzz_0[j] = pb_x * tg_yyzz_xxyzzzz_0[j] + wp_x[j] * tg_yyzz_xxyzzzz_1[j] + fl1_fxn * tg_yyzz_xyzzzz_1[j];

                    tg_xyyzz_xxzzzzz_0[j] = pb_x * tg_yyzz_xxzzzzz_0[j] + wp_x[j] * tg_yyzz_xxzzzzz_1[j] + fl1_fxn * tg_yyzz_xzzzzz_1[j];

                    tg_xyyzz_xyyyyyy_0[j] = pb_x * tg_yyzz_xyyyyyy_0[j] + wp_x[j] * tg_yyzz_xyyyyyy_1[j] + 0.5 * fl1_fxn * tg_yyzz_yyyyyy_1[j];

                    tg_xyyzz_xyyyyyz_0[j] = pb_x * tg_yyzz_xyyyyyz_0[j] + wp_x[j] * tg_yyzz_xyyyyyz_1[j] + 0.5 * fl1_fxn * tg_yyzz_yyyyyz_1[j];

                    tg_xyyzz_xyyyyzz_0[j] = pb_x * tg_yyzz_xyyyyzz_0[j] + wp_x[j] * tg_yyzz_xyyyyzz_1[j] + 0.5 * fl1_fxn * tg_yyzz_yyyyzz_1[j];

                    tg_xyyzz_xyyyzzz_0[j] = pb_x * tg_yyzz_xyyyzzz_0[j] + wp_x[j] * tg_yyzz_xyyyzzz_1[j] + 0.5 * fl1_fxn * tg_yyzz_yyyzzz_1[j];

                    tg_xyyzz_xyyzzzz_0[j] = pb_x * tg_yyzz_xyyzzzz_0[j] + wp_x[j] * tg_yyzz_xyyzzzz_1[j] + 0.5 * fl1_fxn * tg_yyzz_yyzzzz_1[j];

                    tg_xyyzz_xyzzzzz_0[j] = pb_x * tg_yyzz_xyzzzzz_0[j] + wp_x[j] * tg_yyzz_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyzz_yzzzzz_1[j];

                    tg_xyyzz_xzzzzzz_0[j] = pb_x * tg_yyzz_xzzzzzz_0[j] + wp_x[j] * tg_yyzz_xzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyzz_zzzzzz_1[j];

                    tg_xyyzz_yyyyyyy_0[j] = pb_x * tg_yyzz_yyyyyyy_0[j] + wp_x[j] * tg_yyzz_yyyyyyy_1[j];

                    tg_xyyzz_yyyyyyz_0[j] = pb_x * tg_yyzz_yyyyyyz_0[j] + wp_x[j] * tg_yyzz_yyyyyyz_1[j];

                    tg_xyyzz_yyyyyzz_0[j] = pb_x * tg_yyzz_yyyyyzz_0[j] + wp_x[j] * tg_yyzz_yyyyyzz_1[j];

                    tg_xyyzz_yyyyzzz_0[j] = pb_x * tg_yyzz_yyyyzzz_0[j] + wp_x[j] * tg_yyzz_yyyyzzz_1[j];

                    tg_xyyzz_yyyzzzz_0[j] = pb_x * tg_yyzz_yyyzzzz_0[j] + wp_x[j] * tg_yyzz_yyyzzzz_1[j];

                    tg_xyyzz_yyzzzzz_0[j] = pb_x * tg_yyzz_yyzzzzz_0[j] + wp_x[j] * tg_yyzz_yyzzzzz_1[j];

                    tg_xyyzz_yzzzzzz_0[j] = pb_x * tg_yyzz_yzzzzzz_0[j] + wp_x[j] * tg_yyzz_yzzzzzz_1[j];

                    tg_xyyzz_zzzzzzz_0[j] = pb_x * tg_yyzz_zzzzzzz_0[j] + wp_x[j] * tg_yyzz_zzzzzzz_1[j];

                    tg_xyzzz_xxxxxxx_0[j] = pb_x * tg_yzzz_xxxxxxx_0[j] + wp_x[j] * tg_yzzz_xxxxxxx_1[j] + 3.5 * fl1_fxn * tg_yzzz_xxxxxx_1[j];

                    tg_xyzzz_xxxxxxy_0[j] = pb_x * tg_yzzz_xxxxxxy_0[j] + wp_x[j] * tg_yzzz_xxxxxxy_1[j] + 3.0 * fl1_fxn * tg_yzzz_xxxxxy_1[j];

                    tg_xyzzz_xxxxxxz_0[j] = pb_x * tg_yzzz_xxxxxxz_0[j] + wp_x[j] * tg_yzzz_xxxxxxz_1[j] + 3.0 * fl1_fxn * tg_yzzz_xxxxxz_1[j];

                    tg_xyzzz_xxxxxyy_0[j] = pb_x * tg_yzzz_xxxxxyy_0[j] + wp_x[j] * tg_yzzz_xxxxxyy_1[j] + 2.5 * fl1_fxn * tg_yzzz_xxxxyy_1[j];

                    tg_xyzzz_xxxxxyz_0[j] = pb_x * tg_yzzz_xxxxxyz_0[j] + wp_x[j] * tg_yzzz_xxxxxyz_1[j] + 2.5 * fl1_fxn * tg_yzzz_xxxxyz_1[j];

                    tg_xyzzz_xxxxxzz_0[j] = pb_x * tg_yzzz_xxxxxzz_0[j] + wp_x[j] * tg_yzzz_xxxxxzz_1[j] + 2.5 * fl1_fxn * tg_yzzz_xxxxzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSK_474_568(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (474,568)

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
                                             {5, -1, -1, -1},
                                             {7, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_7_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_4_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_4_6_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_3_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {7, -1, -1, -1}, 
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

                auto tg_yyyy_xxxxxxx_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 360); 

                auto tg_yyyy_xxxxxxy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 361); 

                auto tg_yyyy_xxxxxxz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 362); 

                auto tg_yyyy_xxxxxyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 363); 

                auto tg_yyyy_xxxxxyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 364); 

                auto tg_yyyy_xxxxxzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 365); 

                auto tg_yyyy_xxxxyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 366); 

                auto tg_yyyy_xxxxyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 367); 

                auto tg_yyyy_xxxxyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 368); 

                auto tg_yyyy_xxxxzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 369); 

                auto tg_yyyy_xxxyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 370); 

                auto tg_yyyy_xxxyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 371); 

                auto tg_yyyy_xxxyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 372); 

                auto tg_yyyy_xxxyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 373); 

                auto tg_yyyy_xxxzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 374); 

                auto tg_yyyy_xxyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 375); 

                auto tg_yyyy_xxyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 376); 

                auto tg_yyyy_xxyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 377); 

                auto tg_yyyy_xxyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 378); 

                auto tg_yyyy_xxyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 379); 

                auto tg_yyyy_xxzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 380); 

                auto tg_yyyy_xyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 381); 

                auto tg_yyyy_xyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 382); 

                auto tg_yyyy_xyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 383); 

                auto tg_yyyy_xyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 384); 

                auto tg_yyyy_xyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 385); 

                auto tg_yyyy_xyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 386); 

                auto tg_yyyy_xzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 387); 

                auto tg_yzzz_xxxxyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 474); 

                auto tg_yzzz_xxxxyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 475); 

                auto tg_yzzz_xxxxyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 476); 

                auto tg_yzzz_xxxxzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 477); 

                auto tg_yzzz_xxxyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 478); 

                auto tg_yzzz_xxxyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 479); 

                auto tg_yzzz_xxxyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 480); 

                auto tg_yzzz_xxxyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 481); 

                auto tg_yzzz_xxxzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 482); 

                auto tg_yzzz_xxyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 483); 

                auto tg_yzzz_xxyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 484); 

                auto tg_yzzz_xxyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 485); 

                auto tg_yzzz_xxyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 486); 

                auto tg_yzzz_xxyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 487); 

                auto tg_yzzz_xxzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 488); 

                auto tg_yzzz_xyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 489); 

                auto tg_yzzz_xyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 490); 

                auto tg_yzzz_xyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 491); 

                auto tg_yzzz_xyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 492); 

                auto tg_yzzz_xyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 493); 

                auto tg_yzzz_xyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 494); 

                auto tg_yzzz_xzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 495); 

                auto tg_yzzz_yyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 496); 

                auto tg_yzzz_yyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 497); 

                auto tg_yzzz_yyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 498); 

                auto tg_yzzz_yyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 499); 

                auto tg_yzzz_yyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 500); 

                auto tg_yzzz_yyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 501); 

                auto tg_yzzz_yzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 502); 

                auto tg_yzzz_zzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 503); 

                auto tg_zzzz_xxxxxxx_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 504); 

                auto tg_zzzz_xxxxxxy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 505); 

                auto tg_zzzz_xxxxxxz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 506); 

                auto tg_zzzz_xxxxxyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 507); 

                auto tg_zzzz_xxxxxyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 508); 

                auto tg_zzzz_xxxxxzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 509); 

                auto tg_zzzz_xxxxyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 510); 

                auto tg_zzzz_xxxxyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 511); 

                auto tg_zzzz_xxxxyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 512); 

                auto tg_zzzz_xxxxzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 513); 

                auto tg_zzzz_xxxyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 514); 

                auto tg_zzzz_xxxyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 515); 

                auto tg_zzzz_xxxyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 516); 

                auto tg_zzzz_xxxyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 517); 

                auto tg_zzzz_xxxzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 518); 

                auto tg_zzzz_xxyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 519); 

                auto tg_zzzz_xxyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 520); 

                auto tg_zzzz_xxyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 521); 

                auto tg_zzzz_xxyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 522); 

                auto tg_zzzz_xxyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 523); 

                auto tg_zzzz_xxzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 524); 

                auto tg_zzzz_xyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 525); 

                auto tg_zzzz_xyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 526); 

                auto tg_zzzz_xyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 527); 

                auto tg_zzzz_xyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 528); 

                auto tg_zzzz_xyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 529); 

                auto tg_zzzz_xyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 530); 

                auto tg_zzzz_xzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 531); 

                auto tg_zzzz_yyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 532); 

                auto tg_zzzz_yyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 533); 

                auto tg_zzzz_yyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 534); 

                auto tg_zzzz_yyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 535); 

                auto tg_zzzz_yyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 536); 

                auto tg_zzzz_yyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 537); 

                auto tg_zzzz_yzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 538); 

                auto tg_zzzz_zzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 539); 

                auto tg_yyyy_xxxxxxx_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 360); 

                auto tg_yyyy_xxxxxxy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 361); 

                auto tg_yyyy_xxxxxxz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 362); 

                auto tg_yyyy_xxxxxyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 363); 

                auto tg_yyyy_xxxxxyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 364); 

                auto tg_yyyy_xxxxxzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 365); 

                auto tg_yyyy_xxxxyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 366); 

                auto tg_yyyy_xxxxyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 367); 

                auto tg_yyyy_xxxxyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 368); 

                auto tg_yyyy_xxxxzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 369); 

                auto tg_yyyy_xxxyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 370); 

                auto tg_yyyy_xxxyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 371); 

                auto tg_yyyy_xxxyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 372); 

                auto tg_yyyy_xxxyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 373); 

                auto tg_yyyy_xxxzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 374); 

                auto tg_yyyy_xxyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 375); 

                auto tg_yyyy_xxyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 376); 

                auto tg_yyyy_xxyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 377); 

                auto tg_yyyy_xxyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 378); 

                auto tg_yyyy_xxyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 379); 

                auto tg_yyyy_xxzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 380); 

                auto tg_yyyy_xyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 381); 

                auto tg_yyyy_xyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 382); 

                auto tg_yyyy_xyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 383); 

                auto tg_yyyy_xyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 384); 

                auto tg_yyyy_xyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 385); 

                auto tg_yyyy_xyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 386); 

                auto tg_yyyy_xzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 387); 

                auto tg_yzzz_xxxxyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 474); 

                auto tg_yzzz_xxxxyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 475); 

                auto tg_yzzz_xxxxyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 476); 

                auto tg_yzzz_xxxxzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 477); 

                auto tg_yzzz_xxxyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 478); 

                auto tg_yzzz_xxxyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 479); 

                auto tg_yzzz_xxxyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 480); 

                auto tg_yzzz_xxxyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 481); 

                auto tg_yzzz_xxxzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 482); 

                auto tg_yzzz_xxyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 483); 

                auto tg_yzzz_xxyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 484); 

                auto tg_yzzz_xxyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 485); 

                auto tg_yzzz_xxyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 486); 

                auto tg_yzzz_xxyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 487); 

                auto tg_yzzz_xxzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 488); 

                auto tg_yzzz_xyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 489); 

                auto tg_yzzz_xyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 490); 

                auto tg_yzzz_xyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 491); 

                auto tg_yzzz_xyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 492); 

                auto tg_yzzz_xyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 493); 

                auto tg_yzzz_xyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 494); 

                auto tg_yzzz_xzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 495); 

                auto tg_yzzz_yyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 496); 

                auto tg_yzzz_yyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 497); 

                auto tg_yzzz_yyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 498); 

                auto tg_yzzz_yyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 499); 

                auto tg_yzzz_yyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 500); 

                auto tg_yzzz_yyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 501); 

                auto tg_yzzz_yzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 502); 

                auto tg_yzzz_zzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 503); 

                auto tg_zzzz_xxxxxxx_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 504); 

                auto tg_zzzz_xxxxxxy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 505); 

                auto tg_zzzz_xxxxxxz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 506); 

                auto tg_zzzz_xxxxxyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 507); 

                auto tg_zzzz_xxxxxyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 508); 

                auto tg_zzzz_xxxxxzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 509); 

                auto tg_zzzz_xxxxyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 510); 

                auto tg_zzzz_xxxxyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 511); 

                auto tg_zzzz_xxxxyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 512); 

                auto tg_zzzz_xxxxzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 513); 

                auto tg_zzzz_xxxyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 514); 

                auto tg_zzzz_xxxyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 515); 

                auto tg_zzzz_xxxyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 516); 

                auto tg_zzzz_xxxyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 517); 

                auto tg_zzzz_xxxzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 518); 

                auto tg_zzzz_xxyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 519); 

                auto tg_zzzz_xxyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 520); 

                auto tg_zzzz_xxyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 521); 

                auto tg_zzzz_xxyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 522); 

                auto tg_zzzz_xxyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 523); 

                auto tg_zzzz_xxzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 524); 

                auto tg_zzzz_xyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 525); 

                auto tg_zzzz_xyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 526); 

                auto tg_zzzz_xyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 527); 

                auto tg_zzzz_xyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 528); 

                auto tg_zzzz_xyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 529); 

                auto tg_zzzz_xyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 530); 

                auto tg_zzzz_xzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 531); 

                auto tg_zzzz_yyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 532); 

                auto tg_zzzz_yyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 533); 

                auto tg_zzzz_yyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 534); 

                auto tg_zzzz_yyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 535); 

                auto tg_zzzz_yyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 536); 

                auto tg_zzzz_yyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 537); 

                auto tg_zzzz_yzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 538); 

                auto tg_zzzz_zzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 539); 

                auto tg_yyy_xxxxxxx_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 216); 

                auto tg_yyy_xxxxxxy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 217); 

                auto tg_yyy_xxxxxxz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 218); 

                auto tg_yyy_xxxxxyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 219); 

                auto tg_yyy_xxxxxyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 220); 

                auto tg_yyy_xxxxxzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 221); 

                auto tg_yyy_xxxxyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 222); 

                auto tg_yyy_xxxxyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 223); 

                auto tg_yyy_xxxxyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 224); 

                auto tg_yyy_xxxxzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 225); 

                auto tg_yyy_xxxyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 226); 

                auto tg_yyy_xxxyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 227); 

                auto tg_yyy_xxxyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 228); 

                auto tg_yyy_xxxyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 229); 

                auto tg_yyy_xxxzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 230); 

                auto tg_yyy_xxyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 231); 

                auto tg_yyy_xxyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 232); 

                auto tg_yyy_xxyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 233); 

                auto tg_yyy_xxyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 234); 

                auto tg_yyy_xxyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 235); 

                auto tg_yyy_xxzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 236); 

                auto tg_yyy_xyyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 237); 

                auto tg_yyy_xyyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 238); 

                auto tg_yyy_xyyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 239); 

                auto tg_yyy_xyyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 240); 

                auto tg_yyy_xyyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 241); 

                auto tg_yyy_xyzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 242); 

                auto tg_yyy_xzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 243); 

                auto tg_yyy_xxxxxxx_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 216); 

                auto tg_yyy_xxxxxxy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 217); 

                auto tg_yyy_xxxxxxz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 218); 

                auto tg_yyy_xxxxxyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 219); 

                auto tg_yyy_xxxxxyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 220); 

                auto tg_yyy_xxxxxzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 221); 

                auto tg_yyy_xxxxyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 222); 

                auto tg_yyy_xxxxyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 223); 

                auto tg_yyy_xxxxyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 224); 

                auto tg_yyy_xxxxzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 225); 

                auto tg_yyy_xxxyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 226); 

                auto tg_yyy_xxxyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 227); 

                auto tg_yyy_xxxyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 228); 

                auto tg_yyy_xxxyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 229); 

                auto tg_yyy_xxxzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 230); 

                auto tg_yyy_xxyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 231); 

                auto tg_yyy_xxyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 232); 

                auto tg_yyy_xxyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 233); 

                auto tg_yyy_xxyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 234); 

                auto tg_yyy_xxyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 235); 

                auto tg_yyy_xxzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 236); 

                auto tg_yyy_xyyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 237); 

                auto tg_yyy_xyyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 238); 

                auto tg_yyy_xyyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 239); 

                auto tg_yyy_xyyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 240); 

                auto tg_yyy_xyyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 241); 

                auto tg_yyy_xyzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 242); 

                auto tg_yyy_xzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 243); 

                auto tg_yyyy_xxxxxx_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 280); 

                auto tg_yyyy_xxxxxy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 281); 

                auto tg_yyyy_xxxxxz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 282); 

                auto tg_yyyy_xxxxyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 283); 

                auto tg_yyyy_xxxxyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 284); 

                auto tg_yyyy_xxxxzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 285); 

                auto tg_yyyy_xxxyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 286); 

                auto tg_yyyy_xxxyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 287); 

                auto tg_yyyy_xxxyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 288); 

                auto tg_yyyy_xxxzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 289); 

                auto tg_yyyy_xxyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 290); 

                auto tg_yyyy_xxyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 291); 

                auto tg_yyyy_xxyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 292); 

                auto tg_yyyy_xxyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 293); 

                auto tg_yyyy_xxzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 294); 

                auto tg_yyyy_xyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 295); 

                auto tg_yyyy_xyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 296); 

                auto tg_yyyy_xyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 297); 

                auto tg_yyyy_xyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 298); 

                auto tg_yyyy_xyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 299); 

                auto tg_yyyy_xzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 300); 

                auto tg_yzzz_xxxyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 370); 

                auto tg_yzzz_xxxyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 371); 

                auto tg_yzzz_xxxyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 372); 

                auto tg_yzzz_xxxzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 373); 

                auto tg_yzzz_xxyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 374); 

                auto tg_yzzz_xxyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 375); 

                auto tg_yzzz_xxyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 376); 

                auto tg_yzzz_xxyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 377); 

                auto tg_yzzz_xxzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 378); 

                auto tg_yzzz_xyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 379); 

                auto tg_yzzz_xyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 380); 

                auto tg_yzzz_xyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 381); 

                auto tg_yzzz_xyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 382); 

                auto tg_yzzz_xyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 383); 

                auto tg_yzzz_xzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 384); 

                auto tg_yzzz_yyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 385); 

                auto tg_yzzz_yyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 386); 

                auto tg_yzzz_yyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 387); 

                auto tg_yzzz_yyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 388); 

                auto tg_yzzz_yyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 389); 

                auto tg_yzzz_yzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 390); 

                auto tg_yzzz_zzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 391); 

                auto tg_zzzz_xxxxxx_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 392); 

                auto tg_zzzz_xxxxxy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 393); 

                auto tg_zzzz_xxxxxz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 394); 

                auto tg_zzzz_xxxxyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 395); 

                auto tg_zzzz_xxxxyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 396); 

                auto tg_zzzz_xxxxzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 397); 

                auto tg_zzzz_xxxyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 398); 

                auto tg_zzzz_xxxyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 399); 

                auto tg_zzzz_xxxyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 400); 

                auto tg_zzzz_xxxzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 401); 

                auto tg_zzzz_xxyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 402); 

                auto tg_zzzz_xxyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 403); 

                auto tg_zzzz_xxyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 404); 

                auto tg_zzzz_xxyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 405); 

                auto tg_zzzz_xxzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 406); 

                auto tg_zzzz_xyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 407); 

                auto tg_zzzz_xyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 408); 

                auto tg_zzzz_xyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 409); 

                auto tg_zzzz_xyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 410); 

                auto tg_zzzz_xyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 411); 

                auto tg_zzzz_xzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 412); 

                auto tg_zzzz_yyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 413); 

                auto tg_zzzz_yyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 414); 

                auto tg_zzzz_yyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 415); 

                auto tg_zzzz_yyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 416); 

                auto tg_zzzz_yyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 417); 

                auto tg_zzzz_yzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 418); 

                auto tg_zzzz_zzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 419); 

                // set up pointers to integrals

                auto tg_xyzzz_xxxxyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 474); 

                auto tg_xyzzz_xxxxyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 475); 

                auto tg_xyzzz_xxxxyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 476); 

                auto tg_xyzzz_xxxxzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 477); 

                auto tg_xyzzz_xxxyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 478); 

                auto tg_xyzzz_xxxyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 479); 

                auto tg_xyzzz_xxxyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 480); 

                auto tg_xyzzz_xxxyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 481); 

                auto tg_xyzzz_xxxzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 482); 

                auto tg_xyzzz_xxyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 483); 

                auto tg_xyzzz_xxyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 484); 

                auto tg_xyzzz_xxyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 485); 

                auto tg_xyzzz_xxyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 486); 

                auto tg_xyzzz_xxyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 487); 

                auto tg_xyzzz_xxzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 488); 

                auto tg_xyzzz_xyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 489); 

                auto tg_xyzzz_xyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 490); 

                auto tg_xyzzz_xyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 491); 

                auto tg_xyzzz_xyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 492); 

                auto tg_xyzzz_xyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 493); 

                auto tg_xyzzz_xyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 494); 

                auto tg_xyzzz_xzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 495); 

                auto tg_xyzzz_yyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 496); 

                auto tg_xyzzz_yyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 497); 

                auto tg_xyzzz_yyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 498); 

                auto tg_xyzzz_yyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 499); 

                auto tg_xyzzz_yyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 500); 

                auto tg_xyzzz_yyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 501); 

                auto tg_xyzzz_yzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 502); 

                auto tg_xyzzz_zzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 503); 

                auto tg_xzzzz_xxxxxxx_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 504); 

                auto tg_xzzzz_xxxxxxy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 505); 

                auto tg_xzzzz_xxxxxxz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 506); 

                auto tg_xzzzz_xxxxxyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 507); 

                auto tg_xzzzz_xxxxxyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 508); 

                auto tg_xzzzz_xxxxxzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 509); 

                auto tg_xzzzz_xxxxyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 510); 

                auto tg_xzzzz_xxxxyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 511); 

                auto tg_xzzzz_xxxxyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 512); 

                auto tg_xzzzz_xxxxzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 513); 

                auto tg_xzzzz_xxxyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 514); 

                auto tg_xzzzz_xxxyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 515); 

                auto tg_xzzzz_xxxyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 516); 

                auto tg_xzzzz_xxxyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 517); 

                auto tg_xzzzz_xxxzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 518); 

                auto tg_xzzzz_xxyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 519); 

                auto tg_xzzzz_xxyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 520); 

                auto tg_xzzzz_xxyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 521); 

                auto tg_xzzzz_xxyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 522); 

                auto tg_xzzzz_xxyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 523); 

                auto tg_xzzzz_xxzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 524); 

                auto tg_xzzzz_xyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 525); 

                auto tg_xzzzz_xyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 526); 

                auto tg_xzzzz_xyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 527); 

                auto tg_xzzzz_xyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 528); 

                auto tg_xzzzz_xyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 529); 

                auto tg_xzzzz_xyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 530); 

                auto tg_xzzzz_xzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 531); 

                auto tg_xzzzz_yyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 532); 

                auto tg_xzzzz_yyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 533); 

                auto tg_xzzzz_yyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 534); 

                auto tg_xzzzz_yyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 535); 

                auto tg_xzzzz_yyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 536); 

                auto tg_xzzzz_yyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 537); 

                auto tg_xzzzz_yzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 538); 

                auto tg_xzzzz_zzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 539); 

                auto tg_yyyyy_xxxxxxx_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 540); 

                auto tg_yyyyy_xxxxxxy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 541); 

                auto tg_yyyyy_xxxxxxz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 542); 

                auto tg_yyyyy_xxxxxyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 543); 

                auto tg_yyyyy_xxxxxyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 544); 

                auto tg_yyyyy_xxxxxzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 545); 

                auto tg_yyyyy_xxxxyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 546); 

                auto tg_yyyyy_xxxxyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 547); 

                auto tg_yyyyy_xxxxyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 548); 

                auto tg_yyyyy_xxxxzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 549); 

                auto tg_yyyyy_xxxyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 550); 

                auto tg_yyyyy_xxxyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 551); 

                auto tg_yyyyy_xxxyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 552); 

                auto tg_yyyyy_xxxyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 553); 

                auto tg_yyyyy_xxxzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 554); 

                auto tg_yyyyy_xxyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 555); 

                auto tg_yyyyy_xxyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 556); 

                auto tg_yyyyy_xxyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 557); 

                auto tg_yyyyy_xxyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 558); 

                auto tg_yyyyy_xxyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 559); 

                auto tg_yyyyy_xxzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 560); 

                auto tg_yyyyy_xyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 561); 

                auto tg_yyyyy_xyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 562); 

                auto tg_yyyyy_xyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 563); 

                auto tg_yyyyy_xyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 564); 

                auto tg_yyyyy_xyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 565); 

                auto tg_yyyyy_xyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 566); 

                auto tg_yyyyy_xzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 567); 

                // Batch of Integrals (474,568)

                #pragma omp simd aligned(fxn, fza, tg_xyzzz_xxxxyyy_0, tg_xyzzz_xxxxyyz_0, tg_xyzzz_xxxxyzz_0, \
                                         tg_xyzzz_xxxxzzz_0, tg_xyzzz_xxxyyyy_0, tg_xyzzz_xxxyyyz_0, tg_xyzzz_xxxyyzz_0, \
                                         tg_xyzzz_xxxyzzz_0, tg_xyzzz_xxxzzzz_0, tg_xyzzz_xxyyyyy_0, tg_xyzzz_xxyyyyz_0, \
                                         tg_xyzzz_xxyyyzz_0, tg_xyzzz_xxyyzzz_0, tg_xyzzz_xxyzzzz_0, tg_xyzzz_xxzzzzz_0, \
                                         tg_xyzzz_xyyyyyy_0, tg_xyzzz_xyyyyyz_0, tg_xyzzz_xyyyyzz_0, tg_xyzzz_xyyyzzz_0, \
                                         tg_xyzzz_xyyzzzz_0, tg_xyzzz_xyzzzzz_0, tg_xyzzz_xzzzzzz_0, tg_xyzzz_yyyyyyy_0, \
                                         tg_xyzzz_yyyyyyz_0, tg_xyzzz_yyyyyzz_0, tg_xyzzz_yyyyzzz_0, tg_xyzzz_yyyzzzz_0, \
                                         tg_xyzzz_yyzzzzz_0, tg_xyzzz_yzzzzzz_0, tg_xyzzz_zzzzzzz_0, tg_xzzzz_xxxxxxx_0, \
                                         tg_xzzzz_xxxxxxy_0, tg_xzzzz_xxxxxxz_0, tg_xzzzz_xxxxxyy_0, tg_xzzzz_xxxxxyz_0, \
                                         tg_xzzzz_xxxxxzz_0, tg_xzzzz_xxxxyyy_0, tg_xzzzz_xxxxyyz_0, tg_xzzzz_xxxxyzz_0, \
                                         tg_xzzzz_xxxxzzz_0, tg_xzzzz_xxxyyyy_0, tg_xzzzz_xxxyyyz_0, tg_xzzzz_xxxyyzz_0, \
                                         tg_xzzzz_xxxyzzz_0, tg_xzzzz_xxxzzzz_0, tg_xzzzz_xxyyyyy_0, tg_xzzzz_xxyyyyz_0, \
                                         tg_xzzzz_xxyyyzz_0, tg_xzzzz_xxyyzzz_0, tg_xzzzz_xxyzzzz_0, tg_xzzzz_xxzzzzz_0, \
                                         tg_xzzzz_xyyyyyy_0, tg_xzzzz_xyyyyyz_0, tg_xzzzz_xyyyyzz_0, tg_xzzzz_xyyyzzz_0, \
                                         tg_xzzzz_xyyzzzz_0, tg_xzzzz_xyzzzzz_0, tg_xzzzz_xzzzzzz_0, tg_xzzzz_yyyyyyy_0, \
                                         tg_xzzzz_yyyyyyz_0, tg_xzzzz_yyyyyzz_0, tg_xzzzz_yyyyzzz_0, tg_xzzzz_yyyzzzz_0, \
                                         tg_xzzzz_yyzzzzz_0, tg_xzzzz_yzzzzzz_0, tg_xzzzz_zzzzzzz_0, tg_yyy_xxxxxxx_0, \
                                         tg_yyy_xxxxxxx_1, tg_yyy_xxxxxxy_0, tg_yyy_xxxxxxy_1, tg_yyy_xxxxxxz_0, \
                                         tg_yyy_xxxxxxz_1, tg_yyy_xxxxxyy_0, tg_yyy_xxxxxyy_1, tg_yyy_xxxxxyz_0, \
                                         tg_yyy_xxxxxyz_1, tg_yyy_xxxxxzz_0, tg_yyy_xxxxxzz_1, tg_yyy_xxxxyyy_0, \
                                         tg_yyy_xxxxyyy_1, tg_yyy_xxxxyyz_0, tg_yyy_xxxxyyz_1, tg_yyy_xxxxyzz_0, \
                                         tg_yyy_xxxxyzz_1, tg_yyy_xxxxzzz_0, tg_yyy_xxxxzzz_1, tg_yyy_xxxyyyy_0, \
                                         tg_yyy_xxxyyyy_1, tg_yyy_xxxyyyz_0, tg_yyy_xxxyyyz_1, tg_yyy_xxxyyzz_0, \
                                         tg_yyy_xxxyyzz_1, tg_yyy_xxxyzzz_0, tg_yyy_xxxyzzz_1, tg_yyy_xxxzzzz_0, \
                                         tg_yyy_xxxzzzz_1, tg_yyy_xxyyyyy_0, tg_yyy_xxyyyyy_1, tg_yyy_xxyyyyz_0, \
                                         tg_yyy_xxyyyyz_1, tg_yyy_xxyyyzz_0, tg_yyy_xxyyyzz_1, tg_yyy_xxyyzzz_0, \
                                         tg_yyy_xxyyzzz_1, tg_yyy_xxyzzzz_0, tg_yyy_xxyzzzz_1, tg_yyy_xxzzzzz_0, \
                                         tg_yyy_xxzzzzz_1, tg_yyy_xyyyyyy_0, tg_yyy_xyyyyyy_1, tg_yyy_xyyyyyz_0, \
                                         tg_yyy_xyyyyyz_1, tg_yyy_xyyyyzz_0, tg_yyy_xyyyyzz_1, tg_yyy_xyyyzzz_0, \
                                         tg_yyy_xyyyzzz_1, tg_yyy_xyyzzzz_0, tg_yyy_xyyzzzz_1, tg_yyy_xyzzzzz_0, \
                                         tg_yyy_xyzzzzz_1, tg_yyy_xzzzzzz_0, tg_yyy_xzzzzzz_1, tg_yyyy_xxxxxx_1, \
                                         tg_yyyy_xxxxxxx_0, tg_yyyy_xxxxxxx_1, tg_yyyy_xxxxxxy_0, tg_yyyy_xxxxxxy_1, \
                                         tg_yyyy_xxxxxxz_0, tg_yyyy_xxxxxxz_1, tg_yyyy_xxxxxy_1, tg_yyyy_xxxxxyy_0, \
                                         tg_yyyy_xxxxxyy_1, tg_yyyy_xxxxxyz_0, tg_yyyy_xxxxxyz_1, tg_yyyy_xxxxxz_1, \
                                         tg_yyyy_xxxxxzz_0, tg_yyyy_xxxxxzz_1, tg_yyyy_xxxxyy_1, tg_yyyy_xxxxyyy_0, \
                                         tg_yyyy_xxxxyyy_1, tg_yyyy_xxxxyyz_0, tg_yyyy_xxxxyyz_1, tg_yyyy_xxxxyz_1, \
                                         tg_yyyy_xxxxyzz_0, tg_yyyy_xxxxyzz_1, tg_yyyy_xxxxzz_1, tg_yyyy_xxxxzzz_0, \
                                         tg_yyyy_xxxxzzz_1, tg_yyyy_xxxyyy_1, tg_yyyy_xxxyyyy_0, tg_yyyy_xxxyyyy_1, \
                                         tg_yyyy_xxxyyyz_0, tg_yyyy_xxxyyyz_1, tg_yyyy_xxxyyz_1, tg_yyyy_xxxyyzz_0, \
                                         tg_yyyy_xxxyyzz_1, tg_yyyy_xxxyzz_1, tg_yyyy_xxxyzzz_0, tg_yyyy_xxxyzzz_1, \
                                         tg_yyyy_xxxzzz_1, tg_yyyy_xxxzzzz_0, tg_yyyy_xxxzzzz_1, tg_yyyy_xxyyyy_1, \
                                         tg_yyyy_xxyyyyy_0, tg_yyyy_xxyyyyy_1, tg_yyyy_xxyyyyz_0, tg_yyyy_xxyyyyz_1, \
                                         tg_yyyy_xxyyyz_1, tg_yyyy_xxyyyzz_0, tg_yyyy_xxyyyzz_1, tg_yyyy_xxyyzz_1, \
                                         tg_yyyy_xxyyzzz_0, tg_yyyy_xxyyzzz_1, tg_yyyy_xxyzzz_1, tg_yyyy_xxyzzzz_0, \
                                         tg_yyyy_xxyzzzz_1, tg_yyyy_xxzzzz_1, tg_yyyy_xxzzzzz_0, tg_yyyy_xxzzzzz_1, \
                                         tg_yyyy_xyyyyy_1, tg_yyyy_xyyyyyy_0, tg_yyyy_xyyyyyy_1, tg_yyyy_xyyyyyz_0, \
                                         tg_yyyy_xyyyyyz_1, tg_yyyy_xyyyyz_1, tg_yyyy_xyyyyzz_0, tg_yyyy_xyyyyzz_1, \
                                         tg_yyyy_xyyyzz_1, tg_yyyy_xyyyzzz_0, tg_yyyy_xyyyzzz_1, tg_yyyy_xyyzzz_1, \
                                         tg_yyyy_xyyzzzz_0, tg_yyyy_xyyzzzz_1, tg_yyyy_xyzzzz_1, tg_yyyy_xyzzzzz_0, \
                                         tg_yyyy_xyzzzzz_1, tg_yyyy_xzzzzz_1, tg_yyyy_xzzzzzz_0, tg_yyyy_xzzzzzz_1, \
                                         tg_yyyyy_xxxxxxx_0, tg_yyyyy_xxxxxxy_0, tg_yyyyy_xxxxxxz_0, tg_yyyyy_xxxxxyy_0, \
                                         tg_yyyyy_xxxxxyz_0, tg_yyyyy_xxxxxzz_0, tg_yyyyy_xxxxyyy_0, tg_yyyyy_xxxxyyz_0, \
                                         tg_yyyyy_xxxxyzz_0, tg_yyyyy_xxxxzzz_0, tg_yyyyy_xxxyyyy_0, tg_yyyyy_xxxyyyz_0, \
                                         tg_yyyyy_xxxyyzz_0, tg_yyyyy_xxxyzzz_0, tg_yyyyy_xxxzzzz_0, tg_yyyyy_xxyyyyy_0, \
                                         tg_yyyyy_xxyyyyz_0, tg_yyyyy_xxyyyzz_0, tg_yyyyy_xxyyzzz_0, tg_yyyyy_xxyzzzz_0, \
                                         tg_yyyyy_xxzzzzz_0, tg_yyyyy_xyyyyyy_0, tg_yyyyy_xyyyyyz_0, tg_yyyyy_xyyyyzz_0, \
                                         tg_yyyyy_xyyyzzz_0, tg_yyyyy_xyyzzzz_0, tg_yyyyy_xyzzzzz_0, tg_yyyyy_xzzzzzz_0, \
                                         tg_yzzz_xxxxyyy_0, tg_yzzz_xxxxyyy_1, tg_yzzz_xxxxyyz_0, tg_yzzz_xxxxyyz_1, \
                                         tg_yzzz_xxxxyzz_0, tg_yzzz_xxxxyzz_1, tg_yzzz_xxxxzzz_0, tg_yzzz_xxxxzzz_1, \
                                         tg_yzzz_xxxyyy_1, tg_yzzz_xxxyyyy_0, tg_yzzz_xxxyyyy_1, tg_yzzz_xxxyyyz_0, \
                                         tg_yzzz_xxxyyyz_1, tg_yzzz_xxxyyz_1, tg_yzzz_xxxyyzz_0, tg_yzzz_xxxyyzz_1, \
                                         tg_yzzz_xxxyzz_1, tg_yzzz_xxxyzzz_0, tg_yzzz_xxxyzzz_1, tg_yzzz_xxxzzz_1, \
                                         tg_yzzz_xxxzzzz_0, tg_yzzz_xxxzzzz_1, tg_yzzz_xxyyyy_1, tg_yzzz_xxyyyyy_0, \
                                         tg_yzzz_xxyyyyy_1, tg_yzzz_xxyyyyz_0, tg_yzzz_xxyyyyz_1, tg_yzzz_xxyyyz_1, \
                                         tg_yzzz_xxyyyzz_0, tg_yzzz_xxyyyzz_1, tg_yzzz_xxyyzz_1, tg_yzzz_xxyyzzz_0, \
                                         tg_yzzz_xxyyzzz_1, tg_yzzz_xxyzzz_1, tg_yzzz_xxyzzzz_0, tg_yzzz_xxyzzzz_1, \
                                         tg_yzzz_xxzzzz_1, tg_yzzz_xxzzzzz_0, tg_yzzz_xxzzzzz_1, tg_yzzz_xyyyyy_1, \
                                         tg_yzzz_xyyyyyy_0, tg_yzzz_xyyyyyy_1, tg_yzzz_xyyyyyz_0, tg_yzzz_xyyyyyz_1, \
                                         tg_yzzz_xyyyyz_1, tg_yzzz_xyyyyzz_0, tg_yzzz_xyyyyzz_1, tg_yzzz_xyyyzz_1, \
                                         tg_yzzz_xyyyzzz_0, tg_yzzz_xyyyzzz_1, tg_yzzz_xyyzzz_1, tg_yzzz_xyyzzzz_0, \
                                         tg_yzzz_xyyzzzz_1, tg_yzzz_xyzzzz_1, tg_yzzz_xyzzzzz_0, tg_yzzz_xyzzzzz_1, \
                                         tg_yzzz_xzzzzz_1, tg_yzzz_xzzzzzz_0, tg_yzzz_xzzzzzz_1, tg_yzzz_yyyyyy_1, \
                                         tg_yzzz_yyyyyyy_0, tg_yzzz_yyyyyyy_1, tg_yzzz_yyyyyyz_0, tg_yzzz_yyyyyyz_1, \
                                         tg_yzzz_yyyyyz_1, tg_yzzz_yyyyyzz_0, tg_yzzz_yyyyyzz_1, tg_yzzz_yyyyzz_1, \
                                         tg_yzzz_yyyyzzz_0, tg_yzzz_yyyyzzz_1, tg_yzzz_yyyzzz_1, tg_yzzz_yyyzzzz_0, \
                                         tg_yzzz_yyyzzzz_1, tg_yzzz_yyzzzz_1, tg_yzzz_yyzzzzz_0, tg_yzzz_yyzzzzz_1, \
                                         tg_yzzz_yzzzzz_1, tg_yzzz_yzzzzzz_0, tg_yzzz_yzzzzzz_1, tg_yzzz_zzzzzz_1, \
                                         tg_yzzz_zzzzzzz_0, tg_yzzz_zzzzzzz_1, tg_zzzz_xxxxxx_1, tg_zzzz_xxxxxxx_0, \
                                         tg_zzzz_xxxxxxx_1, tg_zzzz_xxxxxxy_0, tg_zzzz_xxxxxxy_1, tg_zzzz_xxxxxxz_0, \
                                         tg_zzzz_xxxxxxz_1, tg_zzzz_xxxxxy_1, tg_zzzz_xxxxxyy_0, tg_zzzz_xxxxxyy_1, \
                                         tg_zzzz_xxxxxyz_0, tg_zzzz_xxxxxyz_1, tg_zzzz_xxxxxz_1, tg_zzzz_xxxxxzz_0, \
                                         tg_zzzz_xxxxxzz_1, tg_zzzz_xxxxyy_1, tg_zzzz_xxxxyyy_0, tg_zzzz_xxxxyyy_1, \
                                         tg_zzzz_xxxxyyz_0, tg_zzzz_xxxxyyz_1, tg_zzzz_xxxxyz_1, tg_zzzz_xxxxyzz_0, \
                                         tg_zzzz_xxxxyzz_1, tg_zzzz_xxxxzz_1, tg_zzzz_xxxxzzz_0, tg_zzzz_xxxxzzz_1, \
                                         tg_zzzz_xxxyyy_1, tg_zzzz_xxxyyyy_0, tg_zzzz_xxxyyyy_1, tg_zzzz_xxxyyyz_0, \
                                         tg_zzzz_xxxyyyz_1, tg_zzzz_xxxyyz_1, tg_zzzz_xxxyyzz_0, tg_zzzz_xxxyyzz_1, \
                                         tg_zzzz_xxxyzz_1, tg_zzzz_xxxyzzz_0, tg_zzzz_xxxyzzz_1, tg_zzzz_xxxzzz_1, \
                                         tg_zzzz_xxxzzzz_0, tg_zzzz_xxxzzzz_1, tg_zzzz_xxyyyy_1, tg_zzzz_xxyyyyy_0, \
                                         tg_zzzz_xxyyyyy_1, tg_zzzz_xxyyyyz_0, tg_zzzz_xxyyyyz_1, tg_zzzz_xxyyyz_1, \
                                         tg_zzzz_xxyyyzz_0, tg_zzzz_xxyyyzz_1, tg_zzzz_xxyyzz_1, tg_zzzz_xxyyzzz_0, \
                                         tg_zzzz_xxyyzzz_1, tg_zzzz_xxyzzz_1, tg_zzzz_xxyzzzz_0, tg_zzzz_xxyzzzz_1, \
                                         tg_zzzz_xxzzzz_1, tg_zzzz_xxzzzzz_0, tg_zzzz_xxzzzzz_1, tg_zzzz_xyyyyy_1, \
                                         tg_zzzz_xyyyyyy_0, tg_zzzz_xyyyyyy_1, tg_zzzz_xyyyyyz_0, tg_zzzz_xyyyyyz_1, \
                                         tg_zzzz_xyyyyz_1, tg_zzzz_xyyyyzz_0, tg_zzzz_xyyyyzz_1, tg_zzzz_xyyyzz_1, \
                                         tg_zzzz_xyyyzzz_0, tg_zzzz_xyyyzzz_1, tg_zzzz_xyyzzz_1, tg_zzzz_xyyzzzz_0, \
                                         tg_zzzz_xyyzzzz_1, tg_zzzz_xyzzzz_1, tg_zzzz_xyzzzzz_0, tg_zzzz_xyzzzzz_1, \
                                         tg_zzzz_xzzzzz_1, tg_zzzz_xzzzzzz_0, tg_zzzz_xzzzzzz_1, tg_zzzz_yyyyyy_1, \
                                         tg_zzzz_yyyyyyy_0, tg_zzzz_yyyyyyy_1, tg_zzzz_yyyyyyz_0, tg_zzzz_yyyyyyz_1, \
                                         tg_zzzz_yyyyyz_1, tg_zzzz_yyyyyzz_0, tg_zzzz_yyyyyzz_1, tg_zzzz_yyyyzz_1, \
                                         tg_zzzz_yyyyzzz_0, tg_zzzz_yyyyzzz_1, tg_zzzz_yyyzzz_1, tg_zzzz_yyyzzzz_0, \
                                         tg_zzzz_yyyzzzz_1, tg_zzzz_yyzzzz_1, tg_zzzz_yyzzzzz_0, tg_zzzz_yyzzzzz_1, \
                                         tg_zzzz_yzzzzz_1, tg_zzzz_yzzzzzz_0, tg_zzzz_yzzzzzz_1, tg_zzzz_zzzzzz_1, \
                                         tg_zzzz_zzzzzzz_0, tg_zzzz_zzzzzzz_1, wp_x, wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xyzzz_xxxxyyy_0[j] = pb_x * tg_yzzz_xxxxyyy_0[j] + wp_x[j] * tg_yzzz_xxxxyyy_1[j] + 2.0 * fl1_fxn * tg_yzzz_xxxyyy_1[j];

                    tg_xyzzz_xxxxyyz_0[j] = pb_x * tg_yzzz_xxxxyyz_0[j] + wp_x[j] * tg_yzzz_xxxxyyz_1[j] + 2.0 * fl1_fxn * tg_yzzz_xxxyyz_1[j];

                    tg_xyzzz_xxxxyzz_0[j] = pb_x * tg_yzzz_xxxxyzz_0[j] + wp_x[j] * tg_yzzz_xxxxyzz_1[j] + 2.0 * fl1_fxn * tg_yzzz_xxxyzz_1[j];

                    tg_xyzzz_xxxxzzz_0[j] = pb_x * tg_yzzz_xxxxzzz_0[j] + wp_x[j] * tg_yzzz_xxxxzzz_1[j] + 2.0 * fl1_fxn * tg_yzzz_xxxzzz_1[j];

                    tg_xyzzz_xxxyyyy_0[j] = pb_x * tg_yzzz_xxxyyyy_0[j] + wp_x[j] * tg_yzzz_xxxyyyy_1[j] + 1.5 * fl1_fxn * tg_yzzz_xxyyyy_1[j];

                    tg_xyzzz_xxxyyyz_0[j] = pb_x * tg_yzzz_xxxyyyz_0[j] + wp_x[j] * tg_yzzz_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_yzzz_xxyyyz_1[j];

                    tg_xyzzz_xxxyyzz_0[j] = pb_x * tg_yzzz_xxxyyzz_0[j] + wp_x[j] * tg_yzzz_xxxyyzz_1[j] + 1.5 * fl1_fxn * tg_yzzz_xxyyzz_1[j];

                    tg_xyzzz_xxxyzzz_0[j] = pb_x * tg_yzzz_xxxyzzz_0[j] + wp_x[j] * tg_yzzz_xxxyzzz_1[j] + 1.5 * fl1_fxn * tg_yzzz_xxyzzz_1[j];

                    tg_xyzzz_xxxzzzz_0[j] = pb_x * tg_yzzz_xxxzzzz_0[j] + wp_x[j] * tg_yzzz_xxxzzzz_1[j] + 1.5 * fl1_fxn * tg_yzzz_xxzzzz_1[j];

                    tg_xyzzz_xxyyyyy_0[j] = pb_x * tg_yzzz_xxyyyyy_0[j] + wp_x[j] * tg_yzzz_xxyyyyy_1[j] + fl1_fxn * tg_yzzz_xyyyyy_1[j];

                    tg_xyzzz_xxyyyyz_0[j] = pb_x * tg_yzzz_xxyyyyz_0[j] + wp_x[j] * tg_yzzz_xxyyyyz_1[j] + fl1_fxn * tg_yzzz_xyyyyz_1[j];

                    tg_xyzzz_xxyyyzz_0[j] = pb_x * tg_yzzz_xxyyyzz_0[j] + wp_x[j] * tg_yzzz_xxyyyzz_1[j] + fl1_fxn * tg_yzzz_xyyyzz_1[j];

                    tg_xyzzz_xxyyzzz_0[j] = pb_x * tg_yzzz_xxyyzzz_0[j] + wp_x[j] * tg_yzzz_xxyyzzz_1[j] + fl1_fxn * tg_yzzz_xyyzzz_1[j];

                    tg_xyzzz_xxyzzzz_0[j] = pb_x * tg_yzzz_xxyzzzz_0[j] + wp_x[j] * tg_yzzz_xxyzzzz_1[j] + fl1_fxn * tg_yzzz_xyzzzz_1[j];

                    tg_xyzzz_xxzzzzz_0[j] = pb_x * tg_yzzz_xxzzzzz_0[j] + wp_x[j] * tg_yzzz_xxzzzzz_1[j] + fl1_fxn * tg_yzzz_xzzzzz_1[j];

                    tg_xyzzz_xyyyyyy_0[j] = pb_x * tg_yzzz_xyyyyyy_0[j] + wp_x[j] * tg_yzzz_xyyyyyy_1[j] + 0.5 * fl1_fxn * tg_yzzz_yyyyyy_1[j];

                    tg_xyzzz_xyyyyyz_0[j] = pb_x * tg_yzzz_xyyyyyz_0[j] + wp_x[j] * tg_yzzz_xyyyyyz_1[j] + 0.5 * fl1_fxn * tg_yzzz_yyyyyz_1[j];

                    tg_xyzzz_xyyyyzz_0[j] = pb_x * tg_yzzz_xyyyyzz_0[j] + wp_x[j] * tg_yzzz_xyyyyzz_1[j] + 0.5 * fl1_fxn * tg_yzzz_yyyyzz_1[j];

                    tg_xyzzz_xyyyzzz_0[j] = pb_x * tg_yzzz_xyyyzzz_0[j] + wp_x[j] * tg_yzzz_xyyyzzz_1[j] + 0.5 * fl1_fxn * tg_yzzz_yyyzzz_1[j];

                    tg_xyzzz_xyyzzzz_0[j] = pb_x * tg_yzzz_xyyzzzz_0[j] + wp_x[j] * tg_yzzz_xyyzzzz_1[j] + 0.5 * fl1_fxn * tg_yzzz_yyzzzz_1[j];

                    tg_xyzzz_xyzzzzz_0[j] = pb_x * tg_yzzz_xyzzzzz_0[j] + wp_x[j] * tg_yzzz_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_yzzz_yzzzzz_1[j];

                    tg_xyzzz_xzzzzzz_0[j] = pb_x * tg_yzzz_xzzzzzz_0[j] + wp_x[j] * tg_yzzz_xzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yzzz_zzzzzz_1[j];

                    tg_xyzzz_yyyyyyy_0[j] = pb_x * tg_yzzz_yyyyyyy_0[j] + wp_x[j] * tg_yzzz_yyyyyyy_1[j];

                    tg_xyzzz_yyyyyyz_0[j] = pb_x * tg_yzzz_yyyyyyz_0[j] + wp_x[j] * tg_yzzz_yyyyyyz_1[j];

                    tg_xyzzz_yyyyyzz_0[j] = pb_x * tg_yzzz_yyyyyzz_0[j] + wp_x[j] * tg_yzzz_yyyyyzz_1[j];

                    tg_xyzzz_yyyyzzz_0[j] = pb_x * tg_yzzz_yyyyzzz_0[j] + wp_x[j] * tg_yzzz_yyyyzzz_1[j];

                    tg_xyzzz_yyyzzzz_0[j] = pb_x * tg_yzzz_yyyzzzz_0[j] + wp_x[j] * tg_yzzz_yyyzzzz_1[j];

                    tg_xyzzz_yyzzzzz_0[j] = pb_x * tg_yzzz_yyzzzzz_0[j] + wp_x[j] * tg_yzzz_yyzzzzz_1[j];

                    tg_xyzzz_yzzzzzz_0[j] = pb_x * tg_yzzz_yzzzzzz_0[j] + wp_x[j] * tg_yzzz_yzzzzzz_1[j];

                    tg_xyzzz_zzzzzzz_0[j] = pb_x * tg_yzzz_zzzzzzz_0[j] + wp_x[j] * tg_yzzz_zzzzzzz_1[j];

                    tg_xzzzz_xxxxxxx_0[j] = pb_x * tg_zzzz_xxxxxxx_0[j] + wp_x[j] * tg_zzzz_xxxxxxx_1[j] + 3.5 * fl1_fxn * tg_zzzz_xxxxxx_1[j];

                    tg_xzzzz_xxxxxxy_0[j] = pb_x * tg_zzzz_xxxxxxy_0[j] + wp_x[j] * tg_zzzz_xxxxxxy_1[j] + 3.0 * fl1_fxn * tg_zzzz_xxxxxy_1[j];

                    tg_xzzzz_xxxxxxz_0[j] = pb_x * tg_zzzz_xxxxxxz_0[j] + wp_x[j] * tg_zzzz_xxxxxxz_1[j] + 3.0 * fl1_fxn * tg_zzzz_xxxxxz_1[j];

                    tg_xzzzz_xxxxxyy_0[j] = pb_x * tg_zzzz_xxxxxyy_0[j] + wp_x[j] * tg_zzzz_xxxxxyy_1[j] + 2.5 * fl1_fxn * tg_zzzz_xxxxyy_1[j];

                    tg_xzzzz_xxxxxyz_0[j] = pb_x * tg_zzzz_xxxxxyz_0[j] + wp_x[j] * tg_zzzz_xxxxxyz_1[j] + 2.5 * fl1_fxn * tg_zzzz_xxxxyz_1[j];

                    tg_xzzzz_xxxxxzz_0[j] = pb_x * tg_zzzz_xxxxxzz_0[j] + wp_x[j] * tg_zzzz_xxxxxzz_1[j] + 2.5 * fl1_fxn * tg_zzzz_xxxxzz_1[j];

                    tg_xzzzz_xxxxyyy_0[j] = pb_x * tg_zzzz_xxxxyyy_0[j] + wp_x[j] * tg_zzzz_xxxxyyy_1[j] + 2.0 * fl1_fxn * tg_zzzz_xxxyyy_1[j];

                    tg_xzzzz_xxxxyyz_0[j] = pb_x * tg_zzzz_xxxxyyz_0[j] + wp_x[j] * tg_zzzz_xxxxyyz_1[j] + 2.0 * fl1_fxn * tg_zzzz_xxxyyz_1[j];

                    tg_xzzzz_xxxxyzz_0[j] = pb_x * tg_zzzz_xxxxyzz_0[j] + wp_x[j] * tg_zzzz_xxxxyzz_1[j] + 2.0 * fl1_fxn * tg_zzzz_xxxyzz_1[j];

                    tg_xzzzz_xxxxzzz_0[j] = pb_x * tg_zzzz_xxxxzzz_0[j] + wp_x[j] * tg_zzzz_xxxxzzz_1[j] + 2.0 * fl1_fxn * tg_zzzz_xxxzzz_1[j];

                    tg_xzzzz_xxxyyyy_0[j] = pb_x * tg_zzzz_xxxyyyy_0[j] + wp_x[j] * tg_zzzz_xxxyyyy_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxyyyy_1[j];

                    tg_xzzzz_xxxyyyz_0[j] = pb_x * tg_zzzz_xxxyyyz_0[j] + wp_x[j] * tg_zzzz_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxyyyz_1[j];

                    tg_xzzzz_xxxyyzz_0[j] = pb_x * tg_zzzz_xxxyyzz_0[j] + wp_x[j] * tg_zzzz_xxxyyzz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxyyzz_1[j];

                    tg_xzzzz_xxxyzzz_0[j] = pb_x * tg_zzzz_xxxyzzz_0[j] + wp_x[j] * tg_zzzz_xxxyzzz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxyzzz_1[j];

                    tg_xzzzz_xxxzzzz_0[j] = pb_x * tg_zzzz_xxxzzzz_0[j] + wp_x[j] * tg_zzzz_xxxzzzz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxzzzz_1[j];

                    tg_xzzzz_xxyyyyy_0[j] = pb_x * tg_zzzz_xxyyyyy_0[j] + wp_x[j] * tg_zzzz_xxyyyyy_1[j] + fl1_fxn * tg_zzzz_xyyyyy_1[j];

                    tg_xzzzz_xxyyyyz_0[j] = pb_x * tg_zzzz_xxyyyyz_0[j] + wp_x[j] * tg_zzzz_xxyyyyz_1[j] + fl1_fxn * tg_zzzz_xyyyyz_1[j];

                    tg_xzzzz_xxyyyzz_0[j] = pb_x * tg_zzzz_xxyyyzz_0[j] + wp_x[j] * tg_zzzz_xxyyyzz_1[j] + fl1_fxn * tg_zzzz_xyyyzz_1[j];

                    tg_xzzzz_xxyyzzz_0[j] = pb_x * tg_zzzz_xxyyzzz_0[j] + wp_x[j] * tg_zzzz_xxyyzzz_1[j] + fl1_fxn * tg_zzzz_xyyzzz_1[j];

                    tg_xzzzz_xxyzzzz_0[j] = pb_x * tg_zzzz_xxyzzzz_0[j] + wp_x[j] * tg_zzzz_xxyzzzz_1[j] + fl1_fxn * tg_zzzz_xyzzzz_1[j];

                    tg_xzzzz_xxzzzzz_0[j] = pb_x * tg_zzzz_xxzzzzz_0[j] + wp_x[j] * tg_zzzz_xxzzzzz_1[j] + fl1_fxn * tg_zzzz_xzzzzz_1[j];

                    tg_xzzzz_xyyyyyy_0[j] = pb_x * tg_zzzz_xyyyyyy_0[j] + wp_x[j] * tg_zzzz_xyyyyyy_1[j] + 0.5 * fl1_fxn * tg_zzzz_yyyyyy_1[j];

                    tg_xzzzz_xyyyyyz_0[j] = pb_x * tg_zzzz_xyyyyyz_0[j] + wp_x[j] * tg_zzzz_xyyyyyz_1[j] + 0.5 * fl1_fxn * tg_zzzz_yyyyyz_1[j];

                    tg_xzzzz_xyyyyzz_0[j] = pb_x * tg_zzzz_xyyyyzz_0[j] + wp_x[j] * tg_zzzz_xyyyyzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_yyyyzz_1[j];

                    tg_xzzzz_xyyyzzz_0[j] = pb_x * tg_zzzz_xyyyzzz_0[j] + wp_x[j] * tg_zzzz_xyyyzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_yyyzzz_1[j];

                    tg_xzzzz_xyyzzzz_0[j] = pb_x * tg_zzzz_xyyzzzz_0[j] + wp_x[j] * tg_zzzz_xyyzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_yyzzzz_1[j];

                    tg_xzzzz_xyzzzzz_0[j] = pb_x * tg_zzzz_xyzzzzz_0[j] + wp_x[j] * tg_zzzz_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_yzzzzz_1[j];

                    tg_xzzzz_xzzzzzz_0[j] = pb_x * tg_zzzz_xzzzzzz_0[j] + wp_x[j] * tg_zzzz_xzzzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_zzzzzz_1[j];

                    tg_xzzzz_yyyyyyy_0[j] = pb_x * tg_zzzz_yyyyyyy_0[j] + wp_x[j] * tg_zzzz_yyyyyyy_1[j];

                    tg_xzzzz_yyyyyyz_0[j] = pb_x * tg_zzzz_yyyyyyz_0[j] + wp_x[j] * tg_zzzz_yyyyyyz_1[j];

                    tg_xzzzz_yyyyyzz_0[j] = pb_x * tg_zzzz_yyyyyzz_0[j] + wp_x[j] * tg_zzzz_yyyyyzz_1[j];

                    tg_xzzzz_yyyyzzz_0[j] = pb_x * tg_zzzz_yyyyzzz_0[j] + wp_x[j] * tg_zzzz_yyyyzzz_1[j];

                    tg_xzzzz_yyyzzzz_0[j] = pb_x * tg_zzzz_yyyzzzz_0[j] + wp_x[j] * tg_zzzz_yyyzzzz_1[j];

                    tg_xzzzz_yyzzzzz_0[j] = pb_x * tg_zzzz_yyzzzzz_0[j] + wp_x[j] * tg_zzzz_yyzzzzz_1[j];

                    tg_xzzzz_yzzzzzz_0[j] = pb_x * tg_zzzz_yzzzzzz_0[j] + wp_x[j] * tg_zzzz_yzzzzzz_1[j];

                    tg_xzzzz_zzzzzzz_0[j] = pb_x * tg_zzzz_zzzzzzz_0[j] + wp_x[j] * tg_zzzz_zzzzzzz_1[j];

                    tg_yyyyy_xxxxxxx_0[j] = pb_y * tg_yyyy_xxxxxxx_0[j] + wp_y[j] * tg_yyyy_xxxxxxx_1[j] + 2.0 * fl1_fx * tg_yyy_xxxxxxx_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxxxxxx_1[j];

                    tg_yyyyy_xxxxxxy_0[j] = pb_y * tg_yyyy_xxxxxxy_0[j] + wp_y[j] * tg_yyyy_xxxxxxy_1[j] + 2.0 * fl1_fx * tg_yyy_xxxxxxy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxxxxxy_1[j] + 0.5 * fl1_fxn * tg_yyyy_xxxxxx_1[j];

                    tg_yyyyy_xxxxxxz_0[j] = pb_y * tg_yyyy_xxxxxxz_0[j] + wp_y[j] * tg_yyyy_xxxxxxz_1[j] + 2.0 * fl1_fx * tg_yyy_xxxxxxz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxxxxxz_1[j];

                    tg_yyyyy_xxxxxyy_0[j] = pb_y * tg_yyyy_xxxxxyy_0[j] + wp_y[j] * tg_yyyy_xxxxxyy_1[j] + 2.0 * fl1_fx * tg_yyy_xxxxxyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxxxxyy_1[j] + fl1_fxn * tg_yyyy_xxxxxy_1[j];

                    tg_yyyyy_xxxxxyz_0[j] = pb_y * tg_yyyy_xxxxxyz_0[j] + wp_y[j] * tg_yyyy_xxxxxyz_1[j] + 2.0 * fl1_fx * tg_yyy_xxxxxyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxxxxyz_1[j] + 0.5 * fl1_fxn * tg_yyyy_xxxxxz_1[j];

                    tg_yyyyy_xxxxxzz_0[j] = pb_y * tg_yyyy_xxxxxzz_0[j] + wp_y[j] * tg_yyyy_xxxxxzz_1[j] + 2.0 * fl1_fx * tg_yyy_xxxxxzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxxxxzz_1[j];

                    tg_yyyyy_xxxxyyy_0[j] = pb_y * tg_yyyy_xxxxyyy_0[j] + wp_y[j] * tg_yyyy_xxxxyyy_1[j] + 2.0 * fl1_fx * tg_yyy_xxxxyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxxxyyy_1[j] + 1.5 * fl1_fxn * tg_yyyy_xxxxyy_1[j];

                    tg_yyyyy_xxxxyyz_0[j] = pb_y * tg_yyyy_xxxxyyz_0[j] + wp_y[j] * tg_yyyy_xxxxyyz_1[j] + 2.0 * fl1_fx * tg_yyy_xxxxyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxxxyyz_1[j] + fl1_fxn * tg_yyyy_xxxxyz_1[j];

                    tg_yyyyy_xxxxyzz_0[j] = pb_y * tg_yyyy_xxxxyzz_0[j] + wp_y[j] * tg_yyyy_xxxxyzz_1[j] + 2.0 * fl1_fx * tg_yyy_xxxxyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxxxyzz_1[j] + 0.5 * fl1_fxn * tg_yyyy_xxxxzz_1[j];

                    tg_yyyyy_xxxxzzz_0[j] = pb_y * tg_yyyy_xxxxzzz_0[j] + wp_y[j] * tg_yyyy_xxxxzzz_1[j] + 2.0 * fl1_fx * tg_yyy_xxxxzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxxxzzz_1[j];

                    tg_yyyyy_xxxyyyy_0[j] = pb_y * tg_yyyy_xxxyyyy_0[j] + wp_y[j] * tg_yyyy_xxxyyyy_1[j] + 2.0 * fl1_fx * tg_yyy_xxxyyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxxyyyy_1[j] + 2.0 * fl1_fxn * tg_yyyy_xxxyyy_1[j];

                    tg_yyyyy_xxxyyyz_0[j] = pb_y * tg_yyyy_xxxyyyz_0[j] + wp_y[j] * tg_yyyy_xxxyyyz_1[j] + 2.0 * fl1_fx * tg_yyy_xxxyyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_yyyy_xxxyyz_1[j];

                    tg_yyyyy_xxxyyzz_0[j] = pb_y * tg_yyyy_xxxyyzz_0[j] + wp_y[j] * tg_yyyy_xxxyyzz_1[j] + 2.0 * fl1_fx * tg_yyy_xxxyyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxxyyzz_1[j] + fl1_fxn * tg_yyyy_xxxyzz_1[j];

                    tg_yyyyy_xxxyzzz_0[j] = pb_y * tg_yyyy_xxxyzzz_0[j] + wp_y[j] * tg_yyyy_xxxyzzz_1[j] + 2.0 * fl1_fx * tg_yyy_xxxyzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxxyzzz_1[j] + 0.5 * fl1_fxn * tg_yyyy_xxxzzz_1[j];

                    tg_yyyyy_xxxzzzz_0[j] = pb_y * tg_yyyy_xxxzzzz_0[j] + wp_y[j] * tg_yyyy_xxxzzzz_1[j] + 2.0 * fl1_fx * tg_yyy_xxxzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxxzzzz_1[j];

                    tg_yyyyy_xxyyyyy_0[j] = pb_y * tg_yyyy_xxyyyyy_0[j] + wp_y[j] * tg_yyyy_xxyyyyy_1[j] + 2.0 * fl1_fx * tg_yyy_xxyyyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxyyyyy_1[j] + 2.5 * fl1_fxn * tg_yyyy_xxyyyy_1[j];

                    tg_yyyyy_xxyyyyz_0[j] = pb_y * tg_yyyy_xxyyyyz_0[j] + wp_y[j] * tg_yyyy_xxyyyyz_1[j] + 2.0 * fl1_fx * tg_yyy_xxyyyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxyyyyz_1[j] + 2.0 * fl1_fxn * tg_yyyy_xxyyyz_1[j];

                    tg_yyyyy_xxyyyzz_0[j] = pb_y * tg_yyyy_xxyyyzz_0[j] + wp_y[j] * tg_yyyy_xxyyyzz_1[j] + 2.0 * fl1_fx * tg_yyy_xxyyyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxyyyzz_1[j] + 1.5 * fl1_fxn * tg_yyyy_xxyyzz_1[j];

                    tg_yyyyy_xxyyzzz_0[j] = pb_y * tg_yyyy_xxyyzzz_0[j] + wp_y[j] * tg_yyyy_xxyyzzz_1[j] + 2.0 * fl1_fx * tg_yyy_xxyyzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxyyzzz_1[j] + fl1_fxn * tg_yyyy_xxyzzz_1[j];

                    tg_yyyyy_xxyzzzz_0[j] = pb_y * tg_yyyy_xxyzzzz_0[j] + wp_y[j] * tg_yyyy_xxyzzzz_1[j] + 2.0 * fl1_fx * tg_yyy_xxyzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxyzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyy_xxzzzz_1[j];

                    tg_yyyyy_xxzzzzz_0[j] = pb_y * tg_yyyy_xxzzzzz_0[j] + wp_y[j] * tg_yyyy_xxzzzzz_1[j] + 2.0 * fl1_fx * tg_yyy_xxzzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxzzzzz_1[j];

                    tg_yyyyy_xyyyyyy_0[j] = pb_y * tg_yyyy_xyyyyyy_0[j] + wp_y[j] * tg_yyyy_xyyyyyy_1[j] + 2.0 * fl1_fx * tg_yyy_xyyyyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xyyyyyy_1[j] + 3.0 * fl1_fxn * tg_yyyy_xyyyyy_1[j];

                    tg_yyyyy_xyyyyyz_0[j] = pb_y * tg_yyyy_xyyyyyz_0[j] + wp_y[j] * tg_yyyy_xyyyyyz_1[j] + 2.0 * fl1_fx * tg_yyy_xyyyyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xyyyyyz_1[j] + 2.5 * fl1_fxn * tg_yyyy_xyyyyz_1[j];

                    tg_yyyyy_xyyyyzz_0[j] = pb_y * tg_yyyy_xyyyyzz_0[j] + wp_y[j] * tg_yyyy_xyyyyzz_1[j] + 2.0 * fl1_fx * tg_yyy_xyyyyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xyyyyzz_1[j] + 2.0 * fl1_fxn * tg_yyyy_xyyyzz_1[j];

                    tg_yyyyy_xyyyzzz_0[j] = pb_y * tg_yyyy_xyyyzzz_0[j] + wp_y[j] * tg_yyyy_xyyyzzz_1[j] + 2.0 * fl1_fx * tg_yyy_xyyyzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xyyyzzz_1[j] + 1.5 * fl1_fxn * tg_yyyy_xyyzzz_1[j];

                    tg_yyyyy_xyyzzzz_0[j] = pb_y * tg_yyyy_xyyzzzz_0[j] + wp_y[j] * tg_yyyy_xyyzzzz_1[j] + 2.0 * fl1_fx * tg_yyy_xyyzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xyyzzzz_1[j] + fl1_fxn * tg_yyyy_xyzzzz_1[j];

                    tg_yyyyy_xyzzzzz_0[j] = pb_y * tg_yyyy_xyzzzzz_0[j] + wp_y[j] * tg_yyyy_xyzzzzz_1[j] + 2.0 * fl1_fx * tg_yyy_xyzzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyy_xzzzzz_1[j];

                    tg_yyyyy_xzzzzzz_0[j] = pb_y * tg_yyyy_xzzzzzz_0[j] + wp_y[j] * tg_yyyy_xzzzzzz_1[j] + 2.0 * fl1_fx * tg_yyy_xzzzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xzzzzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSK_568_662(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (568,662)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {5, -1, -1, -1},
                                             {7, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_7_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_4_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_3_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_4_6_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {6, -1, -1, -1}, 
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

                auto tg_yyyy_yyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 388); 

                auto tg_yyyy_yyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 389); 

                auto tg_yyyy_yyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 390); 

                auto tg_yyyy_yyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 391); 

                auto tg_yyyy_yyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 392); 

                auto tg_yyyy_yyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 393); 

                auto tg_yyyy_yzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 394); 

                auto tg_yyyy_zzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 395); 

                auto tg_yyyz_xxxxxxx_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 396); 

                auto tg_yyyz_xxxxxxy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 397); 

                auto tg_yyyz_xxxxxxz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 398); 

                auto tg_yyyz_xxxxxyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 399); 

                auto tg_yyyz_xxxxxyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 400); 

                auto tg_yyyz_xxxxxzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 401); 

                auto tg_yyyz_xxxxyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 402); 

                auto tg_yyyz_xxxxyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 403); 

                auto tg_yyyz_xxxxyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 404); 

                auto tg_yyyz_xxxxzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 405); 

                auto tg_yyyz_xxxyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 406); 

                auto tg_yyyz_xxxyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 407); 

                auto tg_yyyz_xxxyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 408); 

                auto tg_yyyz_xxxyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 409); 

                auto tg_yyyz_xxxzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 410); 

                auto tg_yyyz_xxyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 411); 

                auto tg_yyyz_xxyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 412); 

                auto tg_yyyz_xxyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 413); 

                auto tg_yyyz_xxyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 414); 

                auto tg_yyyz_xxyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 415); 

                auto tg_yyyz_xxzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 416); 

                auto tg_yyyz_xyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 417); 

                auto tg_yyyz_xyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 418); 

                auto tg_yyyz_xyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 419); 

                auto tg_yyyz_xyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 420); 

                auto tg_yyyz_xyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 421); 

                auto tg_yyyz_xyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 422); 

                auto tg_yyyz_xzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 423); 

                auto tg_yyyz_yyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 424); 

                auto tg_yyyz_yyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 425); 

                auto tg_yyyz_yyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 426); 

                auto tg_yyyz_yyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 427); 

                auto tg_yyyz_yyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 428); 

                auto tg_yyyz_yyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 429); 

                auto tg_yyyz_yzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 430); 

                auto tg_yyyz_zzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 431); 

                auto tg_yyzz_xxxxxxx_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 432); 

                auto tg_yyzz_xxxxxxy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 433); 

                auto tg_yyzz_xxxxxxz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 434); 

                auto tg_yyzz_xxxxxyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 435); 

                auto tg_yyzz_xxxxxyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 436); 

                auto tg_yyzz_xxxxxzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 437); 

                auto tg_yyzz_xxxxyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 438); 

                auto tg_yyzz_xxxxyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 439); 

                auto tg_yyzz_xxxxyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 440); 

                auto tg_yyzz_xxxxzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 441); 

                auto tg_yyzz_xxxyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 442); 

                auto tg_yyzz_xxxyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 443); 

                auto tg_yyzz_xxxyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 444); 

                auto tg_yyzz_xxxyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 445); 

                auto tg_yyzz_xxxzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 446); 

                auto tg_yyzz_xxyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 447); 

                auto tg_yyzz_xxyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 448); 

                auto tg_yyzz_xxyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 449); 

                auto tg_yyzz_xxyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 450); 

                auto tg_yyzz_xxyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 451); 

                auto tg_yyzz_xxzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 452); 

                auto tg_yyzz_xyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 453); 

                auto tg_yyzz_xyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 454); 

                auto tg_yyzz_xyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 455); 

                auto tg_yyzz_xyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 456); 

                auto tg_yyzz_xyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 457); 

                auto tg_yyzz_xyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 458); 

                auto tg_yyzz_xzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 459); 

                auto tg_yyzz_yyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 460); 

                auto tg_yyzz_yyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 461); 

                auto tg_yyzz_yyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 462); 

                auto tg_yyzz_yyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 463); 

                auto tg_yyzz_yyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 464); 

                auto tg_yyzz_yyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 465); 

                auto tg_yyzz_yzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 466); 

                auto tg_yyzz_zzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 467); 

                auto tg_yzzz_xxxxxxx_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 468); 

                auto tg_yzzz_xxxxxxy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 469); 

                auto tg_yzzz_xxxxxxz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 470); 

                auto tg_yzzz_xxxxxyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 471); 

                auto tg_yzzz_xxxxxyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 472); 

                auto tg_yzzz_xxxxxzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 473); 

                auto tg_yzzz_xxxxyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 474); 

                auto tg_yzzz_xxxxyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 475); 

                auto tg_yzzz_xxxxyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 476); 

                auto tg_yzzz_xxxxzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 477); 

                auto tg_yzzz_xxxyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 478); 

                auto tg_yzzz_xxxyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 479); 

                auto tg_yzzz_xxxyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 480); 

                auto tg_yzzz_xxxyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 481); 

                auto tg_yyyy_yyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 388); 

                auto tg_yyyy_yyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 389); 

                auto tg_yyyy_yyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 390); 

                auto tg_yyyy_yyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 391); 

                auto tg_yyyy_yyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 392); 

                auto tg_yyyy_yyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 393); 

                auto tg_yyyy_yzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 394); 

                auto tg_yyyy_zzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 395); 

                auto tg_yyyz_xxxxxxx_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 396); 

                auto tg_yyyz_xxxxxxy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 397); 

                auto tg_yyyz_xxxxxxz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 398); 

                auto tg_yyyz_xxxxxyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 399); 

                auto tg_yyyz_xxxxxyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 400); 

                auto tg_yyyz_xxxxxzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 401); 

                auto tg_yyyz_xxxxyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 402); 

                auto tg_yyyz_xxxxyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 403); 

                auto tg_yyyz_xxxxyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 404); 

                auto tg_yyyz_xxxxzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 405); 

                auto tg_yyyz_xxxyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 406); 

                auto tg_yyyz_xxxyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 407); 

                auto tg_yyyz_xxxyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 408); 

                auto tg_yyyz_xxxyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 409); 

                auto tg_yyyz_xxxzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 410); 

                auto tg_yyyz_xxyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 411); 

                auto tg_yyyz_xxyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 412); 

                auto tg_yyyz_xxyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 413); 

                auto tg_yyyz_xxyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 414); 

                auto tg_yyyz_xxyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 415); 

                auto tg_yyyz_xxzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 416); 

                auto tg_yyyz_xyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 417); 

                auto tg_yyyz_xyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 418); 

                auto tg_yyyz_xyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 419); 

                auto tg_yyyz_xyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 420); 

                auto tg_yyyz_xyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 421); 

                auto tg_yyyz_xyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 422); 

                auto tg_yyyz_xzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 423); 

                auto tg_yyyz_yyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 424); 

                auto tg_yyyz_yyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 425); 

                auto tg_yyyz_yyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 426); 

                auto tg_yyyz_yyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 427); 

                auto tg_yyyz_yyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 428); 

                auto tg_yyyz_yyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 429); 

                auto tg_yyyz_yzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 430); 

                auto tg_yyyz_zzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 431); 

                auto tg_yyzz_xxxxxxx_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 432); 

                auto tg_yyzz_xxxxxxy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 433); 

                auto tg_yyzz_xxxxxxz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 434); 

                auto tg_yyzz_xxxxxyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 435); 

                auto tg_yyzz_xxxxxyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 436); 

                auto tg_yyzz_xxxxxzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 437); 

                auto tg_yyzz_xxxxyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 438); 

                auto tg_yyzz_xxxxyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 439); 

                auto tg_yyzz_xxxxyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 440); 

                auto tg_yyzz_xxxxzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 441); 

                auto tg_yyzz_xxxyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 442); 

                auto tg_yyzz_xxxyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 443); 

                auto tg_yyzz_xxxyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 444); 

                auto tg_yyzz_xxxyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 445); 

                auto tg_yyzz_xxxzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 446); 

                auto tg_yyzz_xxyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 447); 

                auto tg_yyzz_xxyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 448); 

                auto tg_yyzz_xxyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 449); 

                auto tg_yyzz_xxyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 450); 

                auto tg_yyzz_xxyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 451); 

                auto tg_yyzz_xxzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 452); 

                auto tg_yyzz_xyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 453); 

                auto tg_yyzz_xyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 454); 

                auto tg_yyzz_xyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 455); 

                auto tg_yyzz_xyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 456); 

                auto tg_yyzz_xyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 457); 

                auto tg_yyzz_xyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 458); 

                auto tg_yyzz_xzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 459); 

                auto tg_yyzz_yyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 460); 

                auto tg_yyzz_yyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 461); 

                auto tg_yyzz_yyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 462); 

                auto tg_yyzz_yyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 463); 

                auto tg_yyzz_yyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 464); 

                auto tg_yyzz_yyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 465); 

                auto tg_yyzz_yzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 466); 

                auto tg_yyzz_zzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 467); 

                auto tg_yzzz_xxxxxxx_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 468); 

                auto tg_yzzz_xxxxxxy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 469); 

                auto tg_yzzz_xxxxxxz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 470); 

                auto tg_yzzz_xxxxxyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 471); 

                auto tg_yzzz_xxxxxyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 472); 

                auto tg_yzzz_xxxxxzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 473); 

                auto tg_yzzz_xxxxyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 474); 

                auto tg_yzzz_xxxxyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 475); 

                auto tg_yzzz_xxxxyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 476); 

                auto tg_yzzz_xxxxzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 477); 

                auto tg_yzzz_xxxyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 478); 

                auto tg_yzzz_xxxyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 479); 

                auto tg_yzzz_xxxyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 480); 

                auto tg_yzzz_xxxyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 481); 

                auto tg_yyy_yyyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 244); 

                auto tg_yyy_yyyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 245); 

                auto tg_yyy_yyyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 246); 

                auto tg_yyy_yyyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 247); 

                auto tg_yyy_yyyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 248); 

                auto tg_yyy_yyzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 249); 

                auto tg_yyy_yzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 250); 

                auto tg_yyy_zzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 251); 

                auto tg_yyz_xxxxxxx_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 252); 

                auto tg_yyz_xxxxxxy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 253); 

                auto tg_yyz_xxxxxxz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 254); 

                auto tg_yyz_xxxxxyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 255); 

                auto tg_yyz_xxxxxyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 256); 

                auto tg_yyz_xxxxxzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 257); 

                auto tg_yyz_xxxxyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 258); 

                auto tg_yyz_xxxxyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 259); 

                auto tg_yyz_xxxxyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 260); 

                auto tg_yyz_xxxxzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 261); 

                auto tg_yyz_xxxyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 262); 

                auto tg_yyz_xxxyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 263); 

                auto tg_yyz_xxxyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 264); 

                auto tg_yyz_xxxyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 265); 

                auto tg_yyz_xxxzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 266); 

                auto tg_yyz_xxyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 267); 

                auto tg_yyz_xxyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 268); 

                auto tg_yyz_xxyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 269); 

                auto tg_yyz_xxyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 270); 

                auto tg_yyz_xxyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 271); 

                auto tg_yyz_xxzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 272); 

                auto tg_yyz_xyyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 273); 

                auto tg_yyz_xyyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 274); 

                auto tg_yyz_xyyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 275); 

                auto tg_yyz_xyyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 276); 

                auto tg_yyz_xyyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 277); 

                auto tg_yyz_xyzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 278); 

                auto tg_yyz_xzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 279); 

                auto tg_yyz_yyyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 280); 

                auto tg_yyz_yyyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 281); 

                auto tg_yyz_yyyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 282); 

                auto tg_yyz_yyyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 283); 

                auto tg_yyz_yyyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 284); 

                auto tg_yyz_yyzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 285); 

                auto tg_yyz_yzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 286); 

                auto tg_yyz_zzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 287); 

                auto tg_yzz_xxxxxxx_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 288); 

                auto tg_yzz_xxxxxxy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 289); 

                auto tg_yzz_xxxxxxz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 290); 

                auto tg_yzz_xxxxxyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 291); 

                auto tg_yzz_xxxxxyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 292); 

                auto tg_yzz_xxxxxzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 293); 

                auto tg_yzz_xxxxyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 294); 

                auto tg_yzz_xxxxyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 295); 

                auto tg_yzz_xxxxyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 296); 

                auto tg_yzz_xxxxzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 297); 

                auto tg_yzz_xxxyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 298); 

                auto tg_yzz_xxxyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 299); 

                auto tg_yzz_xxxyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 300); 

                auto tg_yzz_xxxyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 301); 

                auto tg_yzz_xxxzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 302); 

                auto tg_yzz_xxyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 303); 

                auto tg_yzz_xxyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 304); 

                auto tg_yzz_xxyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 305); 

                auto tg_yzz_xxyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 306); 

                auto tg_yzz_xxyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 307); 

                auto tg_yzz_xxzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 308); 

                auto tg_yzz_xyyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 309); 

                auto tg_yzz_xyyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 310); 

                auto tg_yzz_xyyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 311); 

                auto tg_yzz_xyyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 312); 

                auto tg_yzz_xyyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 313); 

                auto tg_yzz_xyzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 314); 

                auto tg_yzz_xzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 315); 

                auto tg_yzz_yyyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 316); 

                auto tg_yzz_yyyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 317); 

                auto tg_yzz_yyyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 318); 

                auto tg_yzz_yyyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 319); 

                auto tg_yzz_yyyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 320); 

                auto tg_yzz_yyzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 321); 

                auto tg_yzz_yzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 322); 

                auto tg_yzz_zzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 323); 

                auto tg_zzz_xxxxxxx_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 324); 

                auto tg_zzz_xxxxxxy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 325); 

                auto tg_zzz_xxxxxxz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 326); 

                auto tg_zzz_xxxxxyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 327); 

                auto tg_zzz_xxxxxyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 328); 

                auto tg_zzz_xxxxxzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 329); 

                auto tg_zzz_xxxxyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 330); 

                auto tg_zzz_xxxxyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 331); 

                auto tg_zzz_xxxxyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 332); 

                auto tg_zzz_xxxxzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 333); 

                auto tg_zzz_xxxyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 334); 

                auto tg_zzz_xxxyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 335); 

                auto tg_zzz_xxxyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 336); 

                auto tg_zzz_xxxyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 337); 

                auto tg_yyy_yyyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 244); 

                auto tg_yyy_yyyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 245); 

                auto tg_yyy_yyyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 246); 

                auto tg_yyy_yyyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 247); 

                auto tg_yyy_yyyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 248); 

                auto tg_yyy_yyzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 249); 

                auto tg_yyy_yzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 250); 

                auto tg_yyy_zzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 251); 

                auto tg_yyz_xxxxxxx_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 252); 

                auto tg_yyz_xxxxxxy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 253); 

                auto tg_yyz_xxxxxxz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 254); 

                auto tg_yyz_xxxxxyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 255); 

                auto tg_yyz_xxxxxyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 256); 

                auto tg_yyz_xxxxxzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 257); 

                auto tg_yyz_xxxxyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 258); 

                auto tg_yyz_xxxxyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 259); 

                auto tg_yyz_xxxxyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 260); 

                auto tg_yyz_xxxxzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 261); 

                auto tg_yyz_xxxyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 262); 

                auto tg_yyz_xxxyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 263); 

                auto tg_yyz_xxxyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 264); 

                auto tg_yyz_xxxyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 265); 

                auto tg_yyz_xxxzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 266); 

                auto tg_yyz_xxyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 267); 

                auto tg_yyz_xxyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 268); 

                auto tg_yyz_xxyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 269); 

                auto tg_yyz_xxyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 270); 

                auto tg_yyz_xxyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 271); 

                auto tg_yyz_xxzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 272); 

                auto tg_yyz_xyyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 273); 

                auto tg_yyz_xyyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 274); 

                auto tg_yyz_xyyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 275); 

                auto tg_yyz_xyyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 276); 

                auto tg_yyz_xyyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 277); 

                auto tg_yyz_xyzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 278); 

                auto tg_yyz_xzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 279); 

                auto tg_yyz_yyyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 280); 

                auto tg_yyz_yyyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 281); 

                auto tg_yyz_yyyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 282); 

                auto tg_yyz_yyyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 283); 

                auto tg_yyz_yyyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 284); 

                auto tg_yyz_yyzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 285); 

                auto tg_yyz_yzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 286); 

                auto tg_yyz_zzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 287); 

                auto tg_yzz_xxxxxxx_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 288); 

                auto tg_yzz_xxxxxxy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 289); 

                auto tg_yzz_xxxxxxz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 290); 

                auto tg_yzz_xxxxxyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 291); 

                auto tg_yzz_xxxxxyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 292); 

                auto tg_yzz_xxxxxzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 293); 

                auto tg_yzz_xxxxyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 294); 

                auto tg_yzz_xxxxyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 295); 

                auto tg_yzz_xxxxyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 296); 

                auto tg_yzz_xxxxzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 297); 

                auto tg_yzz_xxxyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 298); 

                auto tg_yzz_xxxyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 299); 

                auto tg_yzz_xxxyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 300); 

                auto tg_yzz_xxxyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 301); 

                auto tg_yzz_xxxzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 302); 

                auto tg_yzz_xxyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 303); 

                auto tg_yzz_xxyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 304); 

                auto tg_yzz_xxyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 305); 

                auto tg_yzz_xxyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 306); 

                auto tg_yzz_xxyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 307); 

                auto tg_yzz_xxzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 308); 

                auto tg_yzz_xyyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 309); 

                auto tg_yzz_xyyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 310); 

                auto tg_yzz_xyyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 311); 

                auto tg_yzz_xyyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 312); 

                auto tg_yzz_xyyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 313); 

                auto tg_yzz_xyzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 314); 

                auto tg_yzz_xzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 315); 

                auto tg_yzz_yyyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 316); 

                auto tg_yzz_yyyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 317); 

                auto tg_yzz_yyyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 318); 

                auto tg_yzz_yyyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 319); 

                auto tg_yzz_yyyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 320); 

                auto tg_yzz_yyzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 321); 

                auto tg_yzz_yzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 322); 

                auto tg_yzz_zzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 323); 

                auto tg_zzz_xxxxxxx_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 324); 

                auto tg_zzz_xxxxxxy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 325); 

                auto tg_zzz_xxxxxxz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 326); 

                auto tg_zzz_xxxxxyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 327); 

                auto tg_zzz_xxxxxyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 328); 

                auto tg_zzz_xxxxxzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 329); 

                auto tg_zzz_xxxxyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 330); 

                auto tg_zzz_xxxxyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 331); 

                auto tg_zzz_xxxxyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 332); 

                auto tg_zzz_xxxxzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 333); 

                auto tg_zzz_xxxyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 334); 

                auto tg_zzz_xxxyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 335); 

                auto tg_zzz_xxxyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 336); 

                auto tg_zzz_xxxyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 337); 

                auto tg_yyyy_yyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 301); 

                auto tg_yyyy_yyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 302); 

                auto tg_yyyy_yyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 303); 

                auto tg_yyyy_yyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 304); 

                auto tg_yyyy_yyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 305); 

                auto tg_yyyy_yzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 306); 

                auto tg_yyyy_zzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 307); 

                auto tg_yyyz_xxxxxx_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 308); 

                auto tg_yyyz_xxxxxy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 309); 

                auto tg_yyyz_xxxxxz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 310); 

                auto tg_yyyz_xxxxyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 311); 

                auto tg_yyyz_xxxxyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 312); 

                auto tg_yyyz_xxxxzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 313); 

                auto tg_yyyz_xxxyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 314); 

                auto tg_yyyz_xxxyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 315); 

                auto tg_yyyz_xxxyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 316); 

                auto tg_yyyz_xxxzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 317); 

                auto tg_yyyz_xxyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 318); 

                auto tg_yyyz_xxyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 319); 

                auto tg_yyyz_xxyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 320); 

                auto tg_yyyz_xxyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 321); 

                auto tg_yyyz_xxzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 322); 

                auto tg_yyyz_xyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 323); 

                auto tg_yyyz_xyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 324); 

                auto tg_yyyz_xyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 325); 

                auto tg_yyyz_xyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 326); 

                auto tg_yyyz_xyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 327); 

                auto tg_yyyz_xzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 328); 

                auto tg_yyyz_yyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 329); 

                auto tg_yyyz_yyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 330); 

                auto tg_yyyz_yyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 331); 

                auto tg_yyyz_yyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 332); 

                auto tg_yyyz_yyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 333); 

                auto tg_yyyz_yzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 334); 

                auto tg_yyyz_zzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 335); 

                auto tg_yyzz_xxxxxx_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 336); 

                auto tg_yyzz_xxxxxy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 337); 

                auto tg_yyzz_xxxxxz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 338); 

                auto tg_yyzz_xxxxyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 339); 

                auto tg_yyzz_xxxxyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 340); 

                auto tg_yyzz_xxxxzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 341); 

                auto tg_yyzz_xxxyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 342); 

                auto tg_yyzz_xxxyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 343); 

                auto tg_yyzz_xxxyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 344); 

                auto tg_yyzz_xxxzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 345); 

                auto tg_yyzz_xxyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 346); 

                auto tg_yyzz_xxyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 347); 

                auto tg_yyzz_xxyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 348); 

                auto tg_yyzz_xxyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 349); 

                auto tg_yyzz_xxzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 350); 

                auto tg_yyzz_xyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 351); 

                auto tg_yyzz_xyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 352); 

                auto tg_yyzz_xyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 353); 

                auto tg_yyzz_xyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 354); 

                auto tg_yyzz_xyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 355); 

                auto tg_yyzz_xzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 356); 

                auto tg_yyzz_yyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 357); 

                auto tg_yyzz_yyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 358); 

                auto tg_yyzz_yyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 359); 

                auto tg_yyzz_yyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 360); 

                auto tg_yyzz_yyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 361); 

                auto tg_yyzz_yzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 362); 

                auto tg_yyzz_zzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 363); 

                auto tg_yzzz_xxxxxx_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 364); 

                auto tg_yzzz_xxxxxy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 365); 

                auto tg_yzzz_xxxxxz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 366); 

                auto tg_yzzz_xxxxyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 367); 

                auto tg_yzzz_xxxxyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 368); 

                auto tg_yzzz_xxxxzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 369); 

                auto tg_yzzz_xxxyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 370); 

                auto tg_yzzz_xxxyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 371); 

                auto tg_yzzz_xxxyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 372); 

                auto tg_yzzz_xxxzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 373); 

                // set up pointers to integrals

                auto tg_yyyyy_yyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 568); 

                auto tg_yyyyy_yyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 569); 

                auto tg_yyyyy_yyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 570); 

                auto tg_yyyyy_yyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 571); 

                auto tg_yyyyy_yyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 572); 

                auto tg_yyyyy_yyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 573); 

                auto tg_yyyyy_yzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 574); 

                auto tg_yyyyy_zzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 575); 

                auto tg_yyyyz_xxxxxxx_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 576); 

                auto tg_yyyyz_xxxxxxy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 577); 

                auto tg_yyyyz_xxxxxxz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 578); 

                auto tg_yyyyz_xxxxxyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 579); 

                auto tg_yyyyz_xxxxxyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 580); 

                auto tg_yyyyz_xxxxxzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 581); 

                auto tg_yyyyz_xxxxyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 582); 

                auto tg_yyyyz_xxxxyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 583); 

                auto tg_yyyyz_xxxxyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 584); 

                auto tg_yyyyz_xxxxzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 585); 

                auto tg_yyyyz_xxxyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 586); 

                auto tg_yyyyz_xxxyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 587); 

                auto tg_yyyyz_xxxyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 588); 

                auto tg_yyyyz_xxxyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 589); 

                auto tg_yyyyz_xxxzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 590); 

                auto tg_yyyyz_xxyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 591); 

                auto tg_yyyyz_xxyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 592); 

                auto tg_yyyyz_xxyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 593); 

                auto tg_yyyyz_xxyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 594); 

                auto tg_yyyyz_xxyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 595); 

                auto tg_yyyyz_xxzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 596); 

                auto tg_yyyyz_xyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 597); 

                auto tg_yyyyz_xyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 598); 

                auto tg_yyyyz_xyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 599); 

                auto tg_yyyyz_xyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 600); 

                auto tg_yyyyz_xyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 601); 

                auto tg_yyyyz_xyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 602); 

                auto tg_yyyyz_xzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 603); 

                auto tg_yyyyz_yyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 604); 

                auto tg_yyyyz_yyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 605); 

                auto tg_yyyyz_yyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 606); 

                auto tg_yyyyz_yyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 607); 

                auto tg_yyyyz_yyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 608); 

                auto tg_yyyyz_yyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 609); 

                auto tg_yyyyz_yzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 610); 

                auto tg_yyyyz_zzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 611); 

                auto tg_yyyzz_xxxxxxx_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 612); 

                auto tg_yyyzz_xxxxxxy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 613); 

                auto tg_yyyzz_xxxxxxz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 614); 

                auto tg_yyyzz_xxxxxyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 615); 

                auto tg_yyyzz_xxxxxyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 616); 

                auto tg_yyyzz_xxxxxzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 617); 

                auto tg_yyyzz_xxxxyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 618); 

                auto tg_yyyzz_xxxxyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 619); 

                auto tg_yyyzz_xxxxyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 620); 

                auto tg_yyyzz_xxxxzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 621); 

                auto tg_yyyzz_xxxyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 622); 

                auto tg_yyyzz_xxxyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 623); 

                auto tg_yyyzz_xxxyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 624); 

                auto tg_yyyzz_xxxyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 625); 

                auto tg_yyyzz_xxxzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 626); 

                auto tg_yyyzz_xxyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 627); 

                auto tg_yyyzz_xxyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 628); 

                auto tg_yyyzz_xxyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 629); 

                auto tg_yyyzz_xxyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 630); 

                auto tg_yyyzz_xxyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 631); 

                auto tg_yyyzz_xxzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 632); 

                auto tg_yyyzz_xyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 633); 

                auto tg_yyyzz_xyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 634); 

                auto tg_yyyzz_xyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 635); 

                auto tg_yyyzz_xyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 636); 

                auto tg_yyyzz_xyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 637); 

                auto tg_yyyzz_xyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 638); 

                auto tg_yyyzz_xzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 639); 

                auto tg_yyyzz_yyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 640); 

                auto tg_yyyzz_yyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 641); 

                auto tg_yyyzz_yyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 642); 

                auto tg_yyyzz_yyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 643); 

                auto tg_yyyzz_yyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 644); 

                auto tg_yyyzz_yyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 645); 

                auto tg_yyyzz_yzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 646); 

                auto tg_yyyzz_zzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 647); 

                auto tg_yyzzz_xxxxxxx_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 648); 

                auto tg_yyzzz_xxxxxxy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 649); 

                auto tg_yyzzz_xxxxxxz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 650); 

                auto tg_yyzzz_xxxxxyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 651); 

                auto tg_yyzzz_xxxxxyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 652); 

                auto tg_yyzzz_xxxxxzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 653); 

                auto tg_yyzzz_xxxxyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 654); 

                auto tg_yyzzz_xxxxyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 655); 

                auto tg_yyzzz_xxxxyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 656); 

                auto tg_yyzzz_xxxxzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 657); 

                auto tg_yyzzz_xxxyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 658); 

                auto tg_yyzzz_xxxyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 659); 

                auto tg_yyzzz_xxxyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 660); 

                auto tg_yyzzz_xxxyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 661); 

                // Batch of Integrals (568,662)

                #pragma omp simd aligned(fxn, fza, tg_yyy_yyyyyyy_0, tg_yyy_yyyyyyy_1, tg_yyy_yyyyyyz_0, \
                                         tg_yyy_yyyyyyz_1, tg_yyy_yyyyyzz_0, tg_yyy_yyyyyzz_1, tg_yyy_yyyyzzz_0, \
                                         tg_yyy_yyyyzzz_1, tg_yyy_yyyzzzz_0, tg_yyy_yyyzzzz_1, tg_yyy_yyzzzzz_0, \
                                         tg_yyy_yyzzzzz_1, tg_yyy_yzzzzzz_0, tg_yyy_yzzzzzz_1, tg_yyy_zzzzzzz_0, \
                                         tg_yyy_zzzzzzz_1, tg_yyyy_yyyyyy_1, tg_yyyy_yyyyyyy_0, tg_yyyy_yyyyyyy_1, \
                                         tg_yyyy_yyyyyyz_0, tg_yyyy_yyyyyyz_1, tg_yyyy_yyyyyz_1, tg_yyyy_yyyyyzz_0, \
                                         tg_yyyy_yyyyyzz_1, tg_yyyy_yyyyzz_1, tg_yyyy_yyyyzzz_0, tg_yyyy_yyyyzzz_1, \
                                         tg_yyyy_yyyzzz_1, tg_yyyy_yyyzzzz_0, tg_yyyy_yyyzzzz_1, tg_yyyy_yyzzzz_1, \
                                         tg_yyyy_yyzzzzz_0, tg_yyyy_yyzzzzz_1, tg_yyyy_yzzzzz_1, tg_yyyy_yzzzzzz_0, \
                                         tg_yyyy_yzzzzzz_1, tg_yyyy_zzzzzz_1, tg_yyyy_zzzzzzz_0, tg_yyyy_zzzzzzz_1, \
                                         tg_yyyyy_yyyyyyy_0, tg_yyyyy_yyyyyyz_0, tg_yyyyy_yyyyyzz_0, tg_yyyyy_yyyyzzz_0, \
                                         tg_yyyyy_yyyzzzz_0, tg_yyyyy_yyzzzzz_0, tg_yyyyy_yzzzzzz_0, tg_yyyyy_zzzzzzz_0, \
                                         tg_yyyyz_xxxxxxx_0, tg_yyyyz_xxxxxxy_0, tg_yyyyz_xxxxxxz_0, tg_yyyyz_xxxxxyy_0, \
                                         tg_yyyyz_xxxxxyz_0, tg_yyyyz_xxxxxzz_0, tg_yyyyz_xxxxyyy_0, tg_yyyyz_xxxxyyz_0, \
                                         tg_yyyyz_xxxxyzz_0, tg_yyyyz_xxxxzzz_0, tg_yyyyz_xxxyyyy_0, tg_yyyyz_xxxyyyz_0, \
                                         tg_yyyyz_xxxyyzz_0, tg_yyyyz_xxxyzzz_0, tg_yyyyz_xxxzzzz_0, tg_yyyyz_xxyyyyy_0, \
                                         tg_yyyyz_xxyyyyz_0, tg_yyyyz_xxyyyzz_0, tg_yyyyz_xxyyzzz_0, tg_yyyyz_xxyzzzz_0, \
                                         tg_yyyyz_xxzzzzz_0, tg_yyyyz_xyyyyyy_0, tg_yyyyz_xyyyyyz_0, tg_yyyyz_xyyyyzz_0, \
                                         tg_yyyyz_xyyyzzz_0, tg_yyyyz_xyyzzzz_0, tg_yyyyz_xyzzzzz_0, tg_yyyyz_xzzzzzz_0, \
                                         tg_yyyyz_yyyyyyy_0, tg_yyyyz_yyyyyyz_0, tg_yyyyz_yyyyyzz_0, tg_yyyyz_yyyyzzz_0, \
                                         tg_yyyyz_yyyzzzz_0, tg_yyyyz_yyzzzzz_0, tg_yyyyz_yzzzzzz_0, tg_yyyyz_zzzzzzz_0, \
                                         tg_yyyz_xxxxxx_1, tg_yyyz_xxxxxxx_0, tg_yyyz_xxxxxxx_1, tg_yyyz_xxxxxxy_0, \
                                         tg_yyyz_xxxxxxy_1, tg_yyyz_xxxxxxz_0, tg_yyyz_xxxxxxz_1, tg_yyyz_xxxxxy_1, \
                                         tg_yyyz_xxxxxyy_0, tg_yyyz_xxxxxyy_1, tg_yyyz_xxxxxyz_0, tg_yyyz_xxxxxyz_1, \
                                         tg_yyyz_xxxxxz_1, tg_yyyz_xxxxxzz_0, tg_yyyz_xxxxxzz_1, tg_yyyz_xxxxyy_1, \
                                         tg_yyyz_xxxxyyy_0, tg_yyyz_xxxxyyy_1, tg_yyyz_xxxxyyz_0, tg_yyyz_xxxxyyz_1, \
                                         tg_yyyz_xxxxyz_1, tg_yyyz_xxxxyzz_0, tg_yyyz_xxxxyzz_1, tg_yyyz_xxxxzz_1, \
                                         tg_yyyz_xxxxzzz_0, tg_yyyz_xxxxzzz_1, tg_yyyz_xxxyyy_1, tg_yyyz_xxxyyyy_0, \
                                         tg_yyyz_xxxyyyy_1, tg_yyyz_xxxyyyz_0, tg_yyyz_xxxyyyz_1, tg_yyyz_xxxyyz_1, \
                                         tg_yyyz_xxxyyzz_0, tg_yyyz_xxxyyzz_1, tg_yyyz_xxxyzz_1, tg_yyyz_xxxyzzz_0, \
                                         tg_yyyz_xxxyzzz_1, tg_yyyz_xxxzzz_1, tg_yyyz_xxxzzzz_0, tg_yyyz_xxxzzzz_1, \
                                         tg_yyyz_xxyyyy_1, tg_yyyz_xxyyyyy_0, tg_yyyz_xxyyyyy_1, tg_yyyz_xxyyyyz_0, \
                                         tg_yyyz_xxyyyyz_1, tg_yyyz_xxyyyz_1, tg_yyyz_xxyyyzz_0, tg_yyyz_xxyyyzz_1, \
                                         tg_yyyz_xxyyzz_1, tg_yyyz_xxyyzzz_0, tg_yyyz_xxyyzzz_1, tg_yyyz_xxyzzz_1, \
                                         tg_yyyz_xxyzzzz_0, tg_yyyz_xxyzzzz_1, tg_yyyz_xxzzzz_1, tg_yyyz_xxzzzzz_0, \
                                         tg_yyyz_xxzzzzz_1, tg_yyyz_xyyyyy_1, tg_yyyz_xyyyyyy_0, tg_yyyz_xyyyyyy_1, \
                                         tg_yyyz_xyyyyyz_0, tg_yyyz_xyyyyyz_1, tg_yyyz_xyyyyz_1, tg_yyyz_xyyyyzz_0, \
                                         tg_yyyz_xyyyyzz_1, tg_yyyz_xyyyzz_1, tg_yyyz_xyyyzzz_0, tg_yyyz_xyyyzzz_1, \
                                         tg_yyyz_xyyzzz_1, tg_yyyz_xyyzzzz_0, tg_yyyz_xyyzzzz_1, tg_yyyz_xyzzzz_1, \
                                         tg_yyyz_xyzzzzz_0, tg_yyyz_xyzzzzz_1, tg_yyyz_xzzzzz_1, tg_yyyz_xzzzzzz_0, \
                                         tg_yyyz_xzzzzzz_1, tg_yyyz_yyyyyy_1, tg_yyyz_yyyyyyy_0, tg_yyyz_yyyyyyy_1, \
                                         tg_yyyz_yyyyyyz_0, tg_yyyz_yyyyyyz_1, tg_yyyz_yyyyyz_1, tg_yyyz_yyyyyzz_0, \
                                         tg_yyyz_yyyyyzz_1, tg_yyyz_yyyyzz_1, tg_yyyz_yyyyzzz_0, tg_yyyz_yyyyzzz_1, \
                                         tg_yyyz_yyyzzz_1, tg_yyyz_yyyzzzz_0, tg_yyyz_yyyzzzz_1, tg_yyyz_yyzzzz_1, \
                                         tg_yyyz_yyzzzzz_0, tg_yyyz_yyzzzzz_1, tg_yyyz_yzzzzz_1, tg_yyyz_yzzzzzz_0, \
                                         tg_yyyz_yzzzzzz_1, tg_yyyz_zzzzzz_1, tg_yyyz_zzzzzzz_0, tg_yyyz_zzzzzzz_1, \
                                         tg_yyyzz_xxxxxxx_0, tg_yyyzz_xxxxxxy_0, tg_yyyzz_xxxxxxz_0, tg_yyyzz_xxxxxyy_0, \
                                         tg_yyyzz_xxxxxyz_0, tg_yyyzz_xxxxxzz_0, tg_yyyzz_xxxxyyy_0, tg_yyyzz_xxxxyyz_0, \
                                         tg_yyyzz_xxxxyzz_0, tg_yyyzz_xxxxzzz_0, tg_yyyzz_xxxyyyy_0, tg_yyyzz_xxxyyyz_0, \
                                         tg_yyyzz_xxxyyzz_0, tg_yyyzz_xxxyzzz_0, tg_yyyzz_xxxzzzz_0, tg_yyyzz_xxyyyyy_0, \
                                         tg_yyyzz_xxyyyyz_0, tg_yyyzz_xxyyyzz_0, tg_yyyzz_xxyyzzz_0, tg_yyyzz_xxyzzzz_0, \
                                         tg_yyyzz_xxzzzzz_0, tg_yyyzz_xyyyyyy_0, tg_yyyzz_xyyyyyz_0, tg_yyyzz_xyyyyzz_0, \
                                         tg_yyyzz_xyyyzzz_0, tg_yyyzz_xyyzzzz_0, tg_yyyzz_xyzzzzz_0, tg_yyyzz_xzzzzzz_0, \
                                         tg_yyyzz_yyyyyyy_0, tg_yyyzz_yyyyyyz_0, tg_yyyzz_yyyyyzz_0, tg_yyyzz_yyyyzzz_0, \
                                         tg_yyyzz_yyyzzzz_0, tg_yyyzz_yyzzzzz_0, tg_yyyzz_yzzzzzz_0, tg_yyyzz_zzzzzzz_0, \
                                         tg_yyz_xxxxxxx_0, tg_yyz_xxxxxxx_1, tg_yyz_xxxxxxy_0, tg_yyz_xxxxxxy_1, \
                                         tg_yyz_xxxxxxz_0, tg_yyz_xxxxxxz_1, tg_yyz_xxxxxyy_0, tg_yyz_xxxxxyy_1, \
                                         tg_yyz_xxxxxyz_0, tg_yyz_xxxxxyz_1, tg_yyz_xxxxxzz_0, tg_yyz_xxxxxzz_1, \
                                         tg_yyz_xxxxyyy_0, tg_yyz_xxxxyyy_1, tg_yyz_xxxxyyz_0, tg_yyz_xxxxyyz_1, \
                                         tg_yyz_xxxxyzz_0, tg_yyz_xxxxyzz_1, tg_yyz_xxxxzzz_0, tg_yyz_xxxxzzz_1, \
                                         tg_yyz_xxxyyyy_0, tg_yyz_xxxyyyy_1, tg_yyz_xxxyyyz_0, tg_yyz_xxxyyyz_1, \
                                         tg_yyz_xxxyyzz_0, tg_yyz_xxxyyzz_1, tg_yyz_xxxyzzz_0, tg_yyz_xxxyzzz_1, \
                                         tg_yyz_xxxzzzz_0, tg_yyz_xxxzzzz_1, tg_yyz_xxyyyyy_0, tg_yyz_xxyyyyy_1, \
                                         tg_yyz_xxyyyyz_0, tg_yyz_xxyyyyz_1, tg_yyz_xxyyyzz_0, tg_yyz_xxyyyzz_1, \
                                         tg_yyz_xxyyzzz_0, tg_yyz_xxyyzzz_1, tg_yyz_xxyzzzz_0, tg_yyz_xxyzzzz_1, \
                                         tg_yyz_xxzzzzz_0, tg_yyz_xxzzzzz_1, tg_yyz_xyyyyyy_0, tg_yyz_xyyyyyy_1, \
                                         tg_yyz_xyyyyyz_0, tg_yyz_xyyyyyz_1, tg_yyz_xyyyyzz_0, tg_yyz_xyyyyzz_1, \
                                         tg_yyz_xyyyzzz_0, tg_yyz_xyyyzzz_1, tg_yyz_xyyzzzz_0, tg_yyz_xyyzzzz_1, \
                                         tg_yyz_xyzzzzz_0, tg_yyz_xyzzzzz_1, tg_yyz_xzzzzzz_0, tg_yyz_xzzzzzz_1, \
                                         tg_yyz_yyyyyyy_0, tg_yyz_yyyyyyy_1, tg_yyz_yyyyyyz_0, tg_yyz_yyyyyyz_1, \
                                         tg_yyz_yyyyyzz_0, tg_yyz_yyyyyzz_1, tg_yyz_yyyyzzz_0, tg_yyz_yyyyzzz_1, \
                                         tg_yyz_yyyzzzz_0, tg_yyz_yyyzzzz_1, tg_yyz_yyzzzzz_0, tg_yyz_yyzzzzz_1, \
                                         tg_yyz_yzzzzzz_0, tg_yyz_yzzzzzz_1, tg_yyz_zzzzzzz_0, tg_yyz_zzzzzzz_1, \
                                         tg_yyzz_xxxxxx_1, tg_yyzz_xxxxxxx_0, tg_yyzz_xxxxxxx_1, tg_yyzz_xxxxxxy_0, \
                                         tg_yyzz_xxxxxxy_1, tg_yyzz_xxxxxxz_0, tg_yyzz_xxxxxxz_1, tg_yyzz_xxxxxy_1, \
                                         tg_yyzz_xxxxxyy_0, tg_yyzz_xxxxxyy_1, tg_yyzz_xxxxxyz_0, tg_yyzz_xxxxxyz_1, \
                                         tg_yyzz_xxxxxz_1, tg_yyzz_xxxxxzz_0, tg_yyzz_xxxxxzz_1, tg_yyzz_xxxxyy_1, \
                                         tg_yyzz_xxxxyyy_0, tg_yyzz_xxxxyyy_1, tg_yyzz_xxxxyyz_0, tg_yyzz_xxxxyyz_1, \
                                         tg_yyzz_xxxxyz_1, tg_yyzz_xxxxyzz_0, tg_yyzz_xxxxyzz_1, tg_yyzz_xxxxzz_1, \
                                         tg_yyzz_xxxxzzz_0, tg_yyzz_xxxxzzz_1, tg_yyzz_xxxyyy_1, tg_yyzz_xxxyyyy_0, \
                                         tg_yyzz_xxxyyyy_1, tg_yyzz_xxxyyyz_0, tg_yyzz_xxxyyyz_1, tg_yyzz_xxxyyz_1, \
                                         tg_yyzz_xxxyyzz_0, tg_yyzz_xxxyyzz_1, tg_yyzz_xxxyzz_1, tg_yyzz_xxxyzzz_0, \
                                         tg_yyzz_xxxyzzz_1, tg_yyzz_xxxzzz_1, tg_yyzz_xxxzzzz_0, tg_yyzz_xxxzzzz_1, \
                                         tg_yyzz_xxyyyy_1, tg_yyzz_xxyyyyy_0, tg_yyzz_xxyyyyy_1, tg_yyzz_xxyyyyz_0, \
                                         tg_yyzz_xxyyyyz_1, tg_yyzz_xxyyyz_1, tg_yyzz_xxyyyzz_0, tg_yyzz_xxyyyzz_1, \
                                         tg_yyzz_xxyyzz_1, tg_yyzz_xxyyzzz_0, tg_yyzz_xxyyzzz_1, tg_yyzz_xxyzzz_1, \
                                         tg_yyzz_xxyzzzz_0, tg_yyzz_xxyzzzz_1, tg_yyzz_xxzzzz_1, tg_yyzz_xxzzzzz_0, \
                                         tg_yyzz_xxzzzzz_1, tg_yyzz_xyyyyy_1, tg_yyzz_xyyyyyy_0, tg_yyzz_xyyyyyy_1, \
                                         tg_yyzz_xyyyyyz_0, tg_yyzz_xyyyyyz_1, tg_yyzz_xyyyyz_1, tg_yyzz_xyyyyzz_0, \
                                         tg_yyzz_xyyyyzz_1, tg_yyzz_xyyyzz_1, tg_yyzz_xyyyzzz_0, tg_yyzz_xyyyzzz_1, \
                                         tg_yyzz_xyyzzz_1, tg_yyzz_xyyzzzz_0, tg_yyzz_xyyzzzz_1, tg_yyzz_xyzzzz_1, \
                                         tg_yyzz_xyzzzzz_0, tg_yyzz_xyzzzzz_1, tg_yyzz_xzzzzz_1, tg_yyzz_xzzzzzz_0, \
                                         tg_yyzz_xzzzzzz_1, tg_yyzz_yyyyyy_1, tg_yyzz_yyyyyyy_0, tg_yyzz_yyyyyyy_1, \
                                         tg_yyzz_yyyyyyz_0, tg_yyzz_yyyyyyz_1, tg_yyzz_yyyyyz_1, tg_yyzz_yyyyyzz_0, \
                                         tg_yyzz_yyyyyzz_1, tg_yyzz_yyyyzz_1, tg_yyzz_yyyyzzz_0, tg_yyzz_yyyyzzz_1, \
                                         tg_yyzz_yyyzzz_1, tg_yyzz_yyyzzzz_0, tg_yyzz_yyyzzzz_1, tg_yyzz_yyzzzz_1, \
                                         tg_yyzz_yyzzzzz_0, tg_yyzz_yyzzzzz_1, tg_yyzz_yzzzzz_1, tg_yyzz_yzzzzzz_0, \
                                         tg_yyzz_yzzzzzz_1, tg_yyzz_zzzzzz_1, tg_yyzz_zzzzzzz_0, tg_yyzz_zzzzzzz_1, \
                                         tg_yyzzz_xxxxxxx_0, tg_yyzzz_xxxxxxy_0, tg_yyzzz_xxxxxxz_0, tg_yyzzz_xxxxxyy_0, \
                                         tg_yyzzz_xxxxxyz_0, tg_yyzzz_xxxxxzz_0, tg_yyzzz_xxxxyyy_0, tg_yyzzz_xxxxyyz_0, \
                                         tg_yyzzz_xxxxyzz_0, tg_yyzzz_xxxxzzz_0, tg_yyzzz_xxxyyyy_0, tg_yyzzz_xxxyyyz_0, \
                                         tg_yyzzz_xxxyyzz_0, tg_yyzzz_xxxyzzz_0, tg_yzz_xxxxxxx_0, tg_yzz_xxxxxxx_1, \
                                         tg_yzz_xxxxxxy_0, tg_yzz_xxxxxxy_1, tg_yzz_xxxxxxz_0, tg_yzz_xxxxxxz_1, \
                                         tg_yzz_xxxxxyy_0, tg_yzz_xxxxxyy_1, tg_yzz_xxxxxyz_0, tg_yzz_xxxxxyz_1, \
                                         tg_yzz_xxxxxzz_0, tg_yzz_xxxxxzz_1, tg_yzz_xxxxyyy_0, tg_yzz_xxxxyyy_1, \
                                         tg_yzz_xxxxyyz_0, tg_yzz_xxxxyyz_1, tg_yzz_xxxxyzz_0, tg_yzz_xxxxyzz_1, \
                                         tg_yzz_xxxxzzz_0, tg_yzz_xxxxzzz_1, tg_yzz_xxxyyyy_0, tg_yzz_xxxyyyy_1, \
                                         tg_yzz_xxxyyyz_0, tg_yzz_xxxyyyz_1, tg_yzz_xxxyyzz_0, tg_yzz_xxxyyzz_1, \
                                         tg_yzz_xxxyzzz_0, tg_yzz_xxxyzzz_1, tg_yzz_xxxzzzz_0, tg_yzz_xxxzzzz_1, \
                                         tg_yzz_xxyyyyy_0, tg_yzz_xxyyyyy_1, tg_yzz_xxyyyyz_0, tg_yzz_xxyyyyz_1, \
                                         tg_yzz_xxyyyzz_0, tg_yzz_xxyyyzz_1, tg_yzz_xxyyzzz_0, tg_yzz_xxyyzzz_1, \
                                         tg_yzz_xxyzzzz_0, tg_yzz_xxyzzzz_1, tg_yzz_xxzzzzz_0, tg_yzz_xxzzzzz_1, \
                                         tg_yzz_xyyyyyy_0, tg_yzz_xyyyyyy_1, tg_yzz_xyyyyyz_0, tg_yzz_xyyyyyz_1, \
                                         tg_yzz_xyyyyzz_0, tg_yzz_xyyyyzz_1, tg_yzz_xyyyzzz_0, tg_yzz_xyyyzzz_1, \
                                         tg_yzz_xyyzzzz_0, tg_yzz_xyyzzzz_1, tg_yzz_xyzzzzz_0, tg_yzz_xyzzzzz_1, \
                                         tg_yzz_xzzzzzz_0, tg_yzz_xzzzzzz_1, tg_yzz_yyyyyyy_0, tg_yzz_yyyyyyy_1, \
                                         tg_yzz_yyyyyyz_0, tg_yzz_yyyyyyz_1, tg_yzz_yyyyyzz_0, tg_yzz_yyyyyzz_1, \
                                         tg_yzz_yyyyzzz_0, tg_yzz_yyyyzzz_1, tg_yzz_yyyzzzz_0, tg_yzz_yyyzzzz_1, \
                                         tg_yzz_yyzzzzz_0, tg_yzz_yyzzzzz_1, tg_yzz_yzzzzzz_0, tg_yzz_yzzzzzz_1, \
                                         tg_yzz_zzzzzzz_0, tg_yzz_zzzzzzz_1, tg_yzzz_xxxxxx_1, tg_yzzz_xxxxxxx_0, \
                                         tg_yzzz_xxxxxxx_1, tg_yzzz_xxxxxxy_0, tg_yzzz_xxxxxxy_1, tg_yzzz_xxxxxxz_0, \
                                         tg_yzzz_xxxxxxz_1, tg_yzzz_xxxxxy_1, tg_yzzz_xxxxxyy_0, tg_yzzz_xxxxxyy_1, \
                                         tg_yzzz_xxxxxyz_0, tg_yzzz_xxxxxyz_1, tg_yzzz_xxxxxz_1, tg_yzzz_xxxxxzz_0, \
                                         tg_yzzz_xxxxxzz_1, tg_yzzz_xxxxyy_1, tg_yzzz_xxxxyyy_0, tg_yzzz_xxxxyyy_1, \
                                         tg_yzzz_xxxxyyz_0, tg_yzzz_xxxxyyz_1, tg_yzzz_xxxxyz_1, tg_yzzz_xxxxyzz_0, \
                                         tg_yzzz_xxxxyzz_1, tg_yzzz_xxxxzz_1, tg_yzzz_xxxxzzz_0, tg_yzzz_xxxxzzz_1, \
                                         tg_yzzz_xxxyyy_1, tg_yzzz_xxxyyyy_0, tg_yzzz_xxxyyyy_1, tg_yzzz_xxxyyyz_0, \
                                         tg_yzzz_xxxyyyz_1, tg_yzzz_xxxyyz_1, tg_yzzz_xxxyyzz_0, tg_yzzz_xxxyyzz_1, \
                                         tg_yzzz_xxxyzz_1, tg_yzzz_xxxyzzz_0, tg_yzzz_xxxyzzz_1, tg_yzzz_xxxzzz_1, \
                                         tg_zzz_xxxxxxx_0, tg_zzz_xxxxxxx_1, tg_zzz_xxxxxxy_0, tg_zzz_xxxxxxy_1, \
                                         tg_zzz_xxxxxxz_0, tg_zzz_xxxxxxz_1, tg_zzz_xxxxxyy_0, tg_zzz_xxxxxyy_1, \
                                         tg_zzz_xxxxxyz_0, tg_zzz_xxxxxyz_1, tg_zzz_xxxxxzz_0, tg_zzz_xxxxxzz_1, \
                                         tg_zzz_xxxxyyy_0, tg_zzz_xxxxyyy_1, tg_zzz_xxxxyyz_0, tg_zzz_xxxxyyz_1, \
                                         tg_zzz_xxxxyzz_0, tg_zzz_xxxxyzz_1, tg_zzz_xxxxzzz_0, tg_zzz_xxxxzzz_1, \
                                         tg_zzz_xxxyyyy_0, tg_zzz_xxxyyyy_1, tg_zzz_xxxyyyz_0, tg_zzz_xxxyyyz_1, \
                                         tg_zzz_xxxyyzz_0, tg_zzz_xxxyyzz_1, tg_zzz_xxxyzzz_0, tg_zzz_xxxyzzz_1, wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_yyyyy_yyyyyyy_0[j] = pb_y * tg_yyyy_yyyyyyy_0[j] + wp_y[j] * tg_yyyy_yyyyyyy_1[j] + 2.0 * fl1_fx * tg_yyy_yyyyyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_yyyyyyy_1[j] + 3.5 * fl1_fxn * tg_yyyy_yyyyyy_1[j];

                    tg_yyyyy_yyyyyyz_0[j] = pb_y * tg_yyyy_yyyyyyz_0[j] + wp_y[j] * tg_yyyy_yyyyyyz_1[j] + 2.0 * fl1_fx * tg_yyy_yyyyyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_yyyyyyz_1[j] + 3.0 * fl1_fxn * tg_yyyy_yyyyyz_1[j];

                    tg_yyyyy_yyyyyzz_0[j] = pb_y * tg_yyyy_yyyyyzz_0[j] + wp_y[j] * tg_yyyy_yyyyyzz_1[j] + 2.0 * fl1_fx * tg_yyy_yyyyyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_yyyyyzz_1[j] + 2.5 * fl1_fxn * tg_yyyy_yyyyzz_1[j];

                    tg_yyyyy_yyyyzzz_0[j] = pb_y * tg_yyyy_yyyyzzz_0[j] + wp_y[j] * tg_yyyy_yyyyzzz_1[j] + 2.0 * fl1_fx * tg_yyy_yyyyzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_yyyyzzz_1[j] + 2.0 * fl1_fxn * tg_yyyy_yyyzzz_1[j];

                    tg_yyyyy_yyyzzzz_0[j] = pb_y * tg_yyyy_yyyzzzz_0[j] + wp_y[j] * tg_yyyy_yyyzzzz_1[j] + 2.0 * fl1_fx * tg_yyy_yyyzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_yyyzzzz_1[j] + 1.5 * fl1_fxn * tg_yyyy_yyzzzz_1[j];

                    tg_yyyyy_yyzzzzz_0[j] = pb_y * tg_yyyy_yyzzzzz_0[j] + wp_y[j] * tg_yyyy_yyzzzzz_1[j] + 2.0 * fl1_fx * tg_yyy_yyzzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_yyzzzzz_1[j] + fl1_fxn * tg_yyyy_yzzzzz_1[j];

                    tg_yyyyy_yzzzzzz_0[j] = pb_y * tg_yyyy_yzzzzzz_0[j] + wp_y[j] * tg_yyyy_yzzzzzz_1[j] + 2.0 * fl1_fx * tg_yyy_yzzzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_yzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyy_zzzzzz_1[j];

                    tg_yyyyy_zzzzzzz_0[j] = pb_y * tg_yyyy_zzzzzzz_0[j] + wp_y[j] * tg_yyyy_zzzzzzz_1[j] + 2.0 * fl1_fx * tg_yyy_zzzzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_zzzzzzz_1[j];

                    tg_yyyyz_xxxxxxx_0[j] = pb_y * tg_yyyz_xxxxxxx_0[j] + wp_y[j] * tg_yyyz_xxxxxxx_1[j] + 1.5 * fl1_fx * tg_yyz_xxxxxxx_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxxxxxx_1[j];

                    tg_yyyyz_xxxxxxy_0[j] = pb_y * tg_yyyz_xxxxxxy_0[j] + wp_y[j] * tg_yyyz_xxxxxxy_1[j] + 1.5 * fl1_fx * tg_yyz_xxxxxxy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxxxxxy_1[j] + 0.5 * fl1_fxn * tg_yyyz_xxxxxx_1[j];

                    tg_yyyyz_xxxxxxz_0[j] = pb_y * tg_yyyz_xxxxxxz_0[j] + wp_y[j] * tg_yyyz_xxxxxxz_1[j] + 1.5 * fl1_fx * tg_yyz_xxxxxxz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxxxxxz_1[j];

                    tg_yyyyz_xxxxxyy_0[j] = pb_y * tg_yyyz_xxxxxyy_0[j] + wp_y[j] * tg_yyyz_xxxxxyy_1[j] + 1.5 * fl1_fx * tg_yyz_xxxxxyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxxxxyy_1[j] + fl1_fxn * tg_yyyz_xxxxxy_1[j];

                    tg_yyyyz_xxxxxyz_0[j] = pb_y * tg_yyyz_xxxxxyz_0[j] + wp_y[j] * tg_yyyz_xxxxxyz_1[j] + 1.5 * fl1_fx * tg_yyz_xxxxxyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxxxxyz_1[j] + 0.5 * fl1_fxn * tg_yyyz_xxxxxz_1[j];

                    tg_yyyyz_xxxxxzz_0[j] = pb_y * tg_yyyz_xxxxxzz_0[j] + wp_y[j] * tg_yyyz_xxxxxzz_1[j] + 1.5 * fl1_fx * tg_yyz_xxxxxzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxxxxzz_1[j];

                    tg_yyyyz_xxxxyyy_0[j] = pb_y * tg_yyyz_xxxxyyy_0[j] + wp_y[j] * tg_yyyz_xxxxyyy_1[j] + 1.5 * fl1_fx * tg_yyz_xxxxyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxxxyyy_1[j] + 1.5 * fl1_fxn * tg_yyyz_xxxxyy_1[j];

                    tg_yyyyz_xxxxyyz_0[j] = pb_y * tg_yyyz_xxxxyyz_0[j] + wp_y[j] * tg_yyyz_xxxxyyz_1[j] + 1.5 * fl1_fx * tg_yyz_xxxxyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxxxyyz_1[j] + fl1_fxn * tg_yyyz_xxxxyz_1[j];

                    tg_yyyyz_xxxxyzz_0[j] = pb_y * tg_yyyz_xxxxyzz_0[j] + wp_y[j] * tg_yyyz_xxxxyzz_1[j] + 1.5 * fl1_fx * tg_yyz_xxxxyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxxxyzz_1[j] + 0.5 * fl1_fxn * tg_yyyz_xxxxzz_1[j];

                    tg_yyyyz_xxxxzzz_0[j] = pb_y * tg_yyyz_xxxxzzz_0[j] + wp_y[j] * tg_yyyz_xxxxzzz_1[j] + 1.5 * fl1_fx * tg_yyz_xxxxzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxxxzzz_1[j];

                    tg_yyyyz_xxxyyyy_0[j] = pb_y * tg_yyyz_xxxyyyy_0[j] + wp_y[j] * tg_yyyz_xxxyyyy_1[j] + 1.5 * fl1_fx * tg_yyz_xxxyyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxxyyyy_1[j] + 2.0 * fl1_fxn * tg_yyyz_xxxyyy_1[j];

                    tg_yyyyz_xxxyyyz_0[j] = pb_y * tg_yyyz_xxxyyyz_0[j] + wp_y[j] * tg_yyyz_xxxyyyz_1[j] + 1.5 * fl1_fx * tg_yyz_xxxyyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_yyyz_xxxyyz_1[j];

                    tg_yyyyz_xxxyyzz_0[j] = pb_y * tg_yyyz_xxxyyzz_0[j] + wp_y[j] * tg_yyyz_xxxyyzz_1[j] + 1.5 * fl1_fx * tg_yyz_xxxyyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxxyyzz_1[j] + fl1_fxn * tg_yyyz_xxxyzz_1[j];

                    tg_yyyyz_xxxyzzz_0[j] = pb_y * tg_yyyz_xxxyzzz_0[j] + wp_y[j] * tg_yyyz_xxxyzzz_1[j] + 1.5 * fl1_fx * tg_yyz_xxxyzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxxyzzz_1[j] + 0.5 * fl1_fxn * tg_yyyz_xxxzzz_1[j];

                    tg_yyyyz_xxxzzzz_0[j] = pb_y * tg_yyyz_xxxzzzz_0[j] + wp_y[j] * tg_yyyz_xxxzzzz_1[j] + 1.5 * fl1_fx * tg_yyz_xxxzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxxzzzz_1[j];

                    tg_yyyyz_xxyyyyy_0[j] = pb_y * tg_yyyz_xxyyyyy_0[j] + wp_y[j] * tg_yyyz_xxyyyyy_1[j] + 1.5 * fl1_fx * tg_yyz_xxyyyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxyyyyy_1[j] + 2.5 * fl1_fxn * tg_yyyz_xxyyyy_1[j];

                    tg_yyyyz_xxyyyyz_0[j] = pb_y * tg_yyyz_xxyyyyz_0[j] + wp_y[j] * tg_yyyz_xxyyyyz_1[j] + 1.5 * fl1_fx * tg_yyz_xxyyyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxyyyyz_1[j] + 2.0 * fl1_fxn * tg_yyyz_xxyyyz_1[j];

                    tg_yyyyz_xxyyyzz_0[j] = pb_y * tg_yyyz_xxyyyzz_0[j] + wp_y[j] * tg_yyyz_xxyyyzz_1[j] + 1.5 * fl1_fx * tg_yyz_xxyyyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxyyyzz_1[j] + 1.5 * fl1_fxn * tg_yyyz_xxyyzz_1[j];

                    tg_yyyyz_xxyyzzz_0[j] = pb_y * tg_yyyz_xxyyzzz_0[j] + wp_y[j] * tg_yyyz_xxyyzzz_1[j] + 1.5 * fl1_fx * tg_yyz_xxyyzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxyyzzz_1[j] + fl1_fxn * tg_yyyz_xxyzzz_1[j];

                    tg_yyyyz_xxyzzzz_0[j] = pb_y * tg_yyyz_xxyzzzz_0[j] + wp_y[j] * tg_yyyz_xxyzzzz_1[j] + 1.5 * fl1_fx * tg_yyz_xxyzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxyzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyz_xxzzzz_1[j];

                    tg_yyyyz_xxzzzzz_0[j] = pb_y * tg_yyyz_xxzzzzz_0[j] + wp_y[j] * tg_yyyz_xxzzzzz_1[j] + 1.5 * fl1_fx * tg_yyz_xxzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxzzzzz_1[j];

                    tg_yyyyz_xyyyyyy_0[j] = pb_y * tg_yyyz_xyyyyyy_0[j] + wp_y[j] * tg_yyyz_xyyyyyy_1[j] + 1.5 * fl1_fx * tg_yyz_xyyyyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xyyyyyy_1[j] + 3.0 * fl1_fxn * tg_yyyz_xyyyyy_1[j];

                    tg_yyyyz_xyyyyyz_0[j] = pb_y * tg_yyyz_xyyyyyz_0[j] + wp_y[j] * tg_yyyz_xyyyyyz_1[j] + 1.5 * fl1_fx * tg_yyz_xyyyyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xyyyyyz_1[j] + 2.5 * fl1_fxn * tg_yyyz_xyyyyz_1[j];

                    tg_yyyyz_xyyyyzz_0[j] = pb_y * tg_yyyz_xyyyyzz_0[j] + wp_y[j] * tg_yyyz_xyyyyzz_1[j] + 1.5 * fl1_fx * tg_yyz_xyyyyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xyyyyzz_1[j] + 2.0 * fl1_fxn * tg_yyyz_xyyyzz_1[j];

                    tg_yyyyz_xyyyzzz_0[j] = pb_y * tg_yyyz_xyyyzzz_0[j] + wp_y[j] * tg_yyyz_xyyyzzz_1[j] + 1.5 * fl1_fx * tg_yyz_xyyyzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xyyyzzz_1[j] + 1.5 * fl1_fxn * tg_yyyz_xyyzzz_1[j];

                    tg_yyyyz_xyyzzzz_0[j] = pb_y * tg_yyyz_xyyzzzz_0[j] + wp_y[j] * tg_yyyz_xyyzzzz_1[j] + 1.5 * fl1_fx * tg_yyz_xyyzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xyyzzzz_1[j] + fl1_fxn * tg_yyyz_xyzzzz_1[j];

                    tg_yyyyz_xyzzzzz_0[j] = pb_y * tg_yyyz_xyzzzzz_0[j] + wp_y[j] * tg_yyyz_xyzzzzz_1[j] + 1.5 * fl1_fx * tg_yyz_xyzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyz_xzzzzz_1[j];

                    tg_yyyyz_xzzzzzz_0[j] = pb_y * tg_yyyz_xzzzzzz_0[j] + wp_y[j] * tg_yyyz_xzzzzzz_1[j] + 1.5 * fl1_fx * tg_yyz_xzzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xzzzzzz_1[j];

                    tg_yyyyz_yyyyyyy_0[j] = pb_y * tg_yyyz_yyyyyyy_0[j] + wp_y[j] * tg_yyyz_yyyyyyy_1[j] + 1.5 * fl1_fx * tg_yyz_yyyyyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_yyyyyyy_1[j] + 3.5 * fl1_fxn * tg_yyyz_yyyyyy_1[j];

                    tg_yyyyz_yyyyyyz_0[j] = pb_y * tg_yyyz_yyyyyyz_0[j] + wp_y[j] * tg_yyyz_yyyyyyz_1[j] + 1.5 * fl1_fx * tg_yyz_yyyyyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_yyyyyyz_1[j] + 3.0 * fl1_fxn * tg_yyyz_yyyyyz_1[j];

                    tg_yyyyz_yyyyyzz_0[j] = pb_y * tg_yyyz_yyyyyzz_0[j] + wp_y[j] * tg_yyyz_yyyyyzz_1[j] + 1.5 * fl1_fx * tg_yyz_yyyyyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_yyyyyzz_1[j] + 2.5 * fl1_fxn * tg_yyyz_yyyyzz_1[j];

                    tg_yyyyz_yyyyzzz_0[j] = pb_y * tg_yyyz_yyyyzzz_0[j] + wp_y[j] * tg_yyyz_yyyyzzz_1[j] + 1.5 * fl1_fx * tg_yyz_yyyyzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_yyyyzzz_1[j] + 2.0 * fl1_fxn * tg_yyyz_yyyzzz_1[j];

                    tg_yyyyz_yyyzzzz_0[j] = pb_y * tg_yyyz_yyyzzzz_0[j] + wp_y[j] * tg_yyyz_yyyzzzz_1[j] + 1.5 * fl1_fx * tg_yyz_yyyzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_yyyzzzz_1[j] + 1.5 * fl1_fxn * tg_yyyz_yyzzzz_1[j];

                    tg_yyyyz_yyzzzzz_0[j] = pb_y * tg_yyyz_yyzzzzz_0[j] + wp_y[j] * tg_yyyz_yyzzzzz_1[j] + 1.5 * fl1_fx * tg_yyz_yyzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_yyzzzzz_1[j] + fl1_fxn * tg_yyyz_yzzzzz_1[j];

                    tg_yyyyz_yzzzzzz_0[j] = pb_y * tg_yyyz_yzzzzzz_0[j] + wp_y[j] * tg_yyyz_yzzzzzz_1[j] + 1.5 * fl1_fx * tg_yyz_yzzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_yzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyz_zzzzzz_1[j];

                    tg_yyyyz_zzzzzzz_0[j] = pb_y * tg_yyyz_zzzzzzz_0[j] + wp_y[j] * tg_yyyz_zzzzzzz_1[j] + 1.5 * fl1_fx * tg_yyz_zzzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_zzzzzzz_1[j];

                    tg_yyyzz_xxxxxxx_0[j] = pb_y * tg_yyzz_xxxxxxx_0[j] + wp_y[j] * tg_yyzz_xxxxxxx_1[j] + fl1_fx * tg_yzz_xxxxxxx_0[j] - fl1_fx * fl1_fza * tg_yzz_xxxxxxx_1[j];

                    tg_yyyzz_xxxxxxy_0[j] = pb_y * tg_yyzz_xxxxxxy_0[j] + wp_y[j] * tg_yyzz_xxxxxxy_1[j] + fl1_fx * tg_yzz_xxxxxxy_0[j] - fl1_fx * fl1_fza * tg_yzz_xxxxxxy_1[j] + 0.5 * fl1_fxn * tg_yyzz_xxxxxx_1[j];

                    tg_yyyzz_xxxxxxz_0[j] = pb_y * tg_yyzz_xxxxxxz_0[j] + wp_y[j] * tg_yyzz_xxxxxxz_1[j] + fl1_fx * tg_yzz_xxxxxxz_0[j] - fl1_fx * fl1_fza * tg_yzz_xxxxxxz_1[j];

                    tg_yyyzz_xxxxxyy_0[j] = pb_y * tg_yyzz_xxxxxyy_0[j] + wp_y[j] * tg_yyzz_xxxxxyy_1[j] + fl1_fx * tg_yzz_xxxxxyy_0[j] - fl1_fx * fl1_fza * tg_yzz_xxxxxyy_1[j] + fl1_fxn * tg_yyzz_xxxxxy_1[j];

                    tg_yyyzz_xxxxxyz_0[j] = pb_y * tg_yyzz_xxxxxyz_0[j] + wp_y[j] * tg_yyzz_xxxxxyz_1[j] + fl1_fx * tg_yzz_xxxxxyz_0[j] - fl1_fx * fl1_fza * tg_yzz_xxxxxyz_1[j] + 0.5 * fl1_fxn * tg_yyzz_xxxxxz_1[j];

                    tg_yyyzz_xxxxxzz_0[j] = pb_y * tg_yyzz_xxxxxzz_0[j] + wp_y[j] * tg_yyzz_xxxxxzz_1[j] + fl1_fx * tg_yzz_xxxxxzz_0[j] - fl1_fx * fl1_fza * tg_yzz_xxxxxzz_1[j];

                    tg_yyyzz_xxxxyyy_0[j] = pb_y * tg_yyzz_xxxxyyy_0[j] + wp_y[j] * tg_yyzz_xxxxyyy_1[j] + fl1_fx * tg_yzz_xxxxyyy_0[j] - fl1_fx * fl1_fza * tg_yzz_xxxxyyy_1[j] + 1.5 * fl1_fxn * tg_yyzz_xxxxyy_1[j];

                    tg_yyyzz_xxxxyyz_0[j] = pb_y * tg_yyzz_xxxxyyz_0[j] + wp_y[j] * tg_yyzz_xxxxyyz_1[j] + fl1_fx * tg_yzz_xxxxyyz_0[j] - fl1_fx * fl1_fza * tg_yzz_xxxxyyz_1[j] + fl1_fxn * tg_yyzz_xxxxyz_1[j];

                    tg_yyyzz_xxxxyzz_0[j] = pb_y * tg_yyzz_xxxxyzz_0[j] + wp_y[j] * tg_yyzz_xxxxyzz_1[j] + fl1_fx * tg_yzz_xxxxyzz_0[j] - fl1_fx * fl1_fza * tg_yzz_xxxxyzz_1[j] + 0.5 * fl1_fxn * tg_yyzz_xxxxzz_1[j];

                    tg_yyyzz_xxxxzzz_0[j] = pb_y * tg_yyzz_xxxxzzz_0[j] + wp_y[j] * tg_yyzz_xxxxzzz_1[j] + fl1_fx * tg_yzz_xxxxzzz_0[j] - fl1_fx * fl1_fza * tg_yzz_xxxxzzz_1[j];

                    tg_yyyzz_xxxyyyy_0[j] = pb_y * tg_yyzz_xxxyyyy_0[j] + wp_y[j] * tg_yyzz_xxxyyyy_1[j] + fl1_fx * tg_yzz_xxxyyyy_0[j] - fl1_fx * fl1_fza * tg_yzz_xxxyyyy_1[j] + 2.0 * fl1_fxn * tg_yyzz_xxxyyy_1[j];

                    tg_yyyzz_xxxyyyz_0[j] = pb_y * tg_yyzz_xxxyyyz_0[j] + wp_y[j] * tg_yyzz_xxxyyyz_1[j] + fl1_fx * tg_yzz_xxxyyyz_0[j] - fl1_fx * fl1_fza * tg_yzz_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_yyzz_xxxyyz_1[j];

                    tg_yyyzz_xxxyyzz_0[j] = pb_y * tg_yyzz_xxxyyzz_0[j] + wp_y[j] * tg_yyzz_xxxyyzz_1[j] + fl1_fx * tg_yzz_xxxyyzz_0[j] - fl1_fx * fl1_fza * tg_yzz_xxxyyzz_1[j] + fl1_fxn * tg_yyzz_xxxyzz_1[j];

                    tg_yyyzz_xxxyzzz_0[j] = pb_y * tg_yyzz_xxxyzzz_0[j] + wp_y[j] * tg_yyzz_xxxyzzz_1[j] + fl1_fx * tg_yzz_xxxyzzz_0[j] - fl1_fx * fl1_fza * tg_yzz_xxxyzzz_1[j] + 0.5 * fl1_fxn * tg_yyzz_xxxzzz_1[j];

                    tg_yyyzz_xxxzzzz_0[j] = pb_y * tg_yyzz_xxxzzzz_0[j] + wp_y[j] * tg_yyzz_xxxzzzz_1[j] + fl1_fx * tg_yzz_xxxzzzz_0[j] - fl1_fx * fl1_fza * tg_yzz_xxxzzzz_1[j];

                    tg_yyyzz_xxyyyyy_0[j] = pb_y * tg_yyzz_xxyyyyy_0[j] + wp_y[j] * tg_yyzz_xxyyyyy_1[j] + fl1_fx * tg_yzz_xxyyyyy_0[j] - fl1_fx * fl1_fza * tg_yzz_xxyyyyy_1[j] + 2.5 * fl1_fxn * tg_yyzz_xxyyyy_1[j];

                    tg_yyyzz_xxyyyyz_0[j] = pb_y * tg_yyzz_xxyyyyz_0[j] + wp_y[j] * tg_yyzz_xxyyyyz_1[j] + fl1_fx * tg_yzz_xxyyyyz_0[j] - fl1_fx * fl1_fza * tg_yzz_xxyyyyz_1[j] + 2.0 * fl1_fxn * tg_yyzz_xxyyyz_1[j];

                    tg_yyyzz_xxyyyzz_0[j] = pb_y * tg_yyzz_xxyyyzz_0[j] + wp_y[j] * tg_yyzz_xxyyyzz_1[j] + fl1_fx * tg_yzz_xxyyyzz_0[j] - fl1_fx * fl1_fza * tg_yzz_xxyyyzz_1[j] + 1.5 * fl1_fxn * tg_yyzz_xxyyzz_1[j];

                    tg_yyyzz_xxyyzzz_0[j] = pb_y * tg_yyzz_xxyyzzz_0[j] + wp_y[j] * tg_yyzz_xxyyzzz_1[j] + fl1_fx * tg_yzz_xxyyzzz_0[j] - fl1_fx * fl1_fza * tg_yzz_xxyyzzz_1[j] + fl1_fxn * tg_yyzz_xxyzzz_1[j];

                    tg_yyyzz_xxyzzzz_0[j] = pb_y * tg_yyzz_xxyzzzz_0[j] + wp_y[j] * tg_yyzz_xxyzzzz_1[j] + fl1_fx * tg_yzz_xxyzzzz_0[j] - fl1_fx * fl1_fza * tg_yzz_xxyzzzz_1[j] + 0.5 * fl1_fxn * tg_yyzz_xxzzzz_1[j];

                    tg_yyyzz_xxzzzzz_0[j] = pb_y * tg_yyzz_xxzzzzz_0[j] + wp_y[j] * tg_yyzz_xxzzzzz_1[j] + fl1_fx * tg_yzz_xxzzzzz_0[j] - fl1_fx * fl1_fza * tg_yzz_xxzzzzz_1[j];

                    tg_yyyzz_xyyyyyy_0[j] = pb_y * tg_yyzz_xyyyyyy_0[j] + wp_y[j] * tg_yyzz_xyyyyyy_1[j] + fl1_fx * tg_yzz_xyyyyyy_0[j] - fl1_fx * fl1_fza * tg_yzz_xyyyyyy_1[j] + 3.0 * fl1_fxn * tg_yyzz_xyyyyy_1[j];

                    tg_yyyzz_xyyyyyz_0[j] = pb_y * tg_yyzz_xyyyyyz_0[j] + wp_y[j] * tg_yyzz_xyyyyyz_1[j] + fl1_fx * tg_yzz_xyyyyyz_0[j] - fl1_fx * fl1_fza * tg_yzz_xyyyyyz_1[j] + 2.5 * fl1_fxn * tg_yyzz_xyyyyz_1[j];

                    tg_yyyzz_xyyyyzz_0[j] = pb_y * tg_yyzz_xyyyyzz_0[j] + wp_y[j] * tg_yyzz_xyyyyzz_1[j] + fl1_fx * tg_yzz_xyyyyzz_0[j] - fl1_fx * fl1_fza * tg_yzz_xyyyyzz_1[j] + 2.0 * fl1_fxn * tg_yyzz_xyyyzz_1[j];

                    tg_yyyzz_xyyyzzz_0[j] = pb_y * tg_yyzz_xyyyzzz_0[j] + wp_y[j] * tg_yyzz_xyyyzzz_1[j] + fl1_fx * tg_yzz_xyyyzzz_0[j] - fl1_fx * fl1_fza * tg_yzz_xyyyzzz_1[j] + 1.5 * fl1_fxn * tg_yyzz_xyyzzz_1[j];

                    tg_yyyzz_xyyzzzz_0[j] = pb_y * tg_yyzz_xyyzzzz_0[j] + wp_y[j] * tg_yyzz_xyyzzzz_1[j] + fl1_fx * tg_yzz_xyyzzzz_0[j] - fl1_fx * fl1_fza * tg_yzz_xyyzzzz_1[j] + fl1_fxn * tg_yyzz_xyzzzz_1[j];

                    tg_yyyzz_xyzzzzz_0[j] = pb_y * tg_yyzz_xyzzzzz_0[j] + wp_y[j] * tg_yyzz_xyzzzzz_1[j] + fl1_fx * tg_yzz_xyzzzzz_0[j] - fl1_fx * fl1_fza * tg_yzz_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyzz_xzzzzz_1[j];

                    tg_yyyzz_xzzzzzz_0[j] = pb_y * tg_yyzz_xzzzzzz_0[j] + wp_y[j] * tg_yyzz_xzzzzzz_1[j] + fl1_fx * tg_yzz_xzzzzzz_0[j] - fl1_fx * fl1_fza * tg_yzz_xzzzzzz_1[j];

                    tg_yyyzz_yyyyyyy_0[j] = pb_y * tg_yyzz_yyyyyyy_0[j] + wp_y[j] * tg_yyzz_yyyyyyy_1[j] + fl1_fx * tg_yzz_yyyyyyy_0[j] - fl1_fx * fl1_fza * tg_yzz_yyyyyyy_1[j] + 3.5 * fl1_fxn * tg_yyzz_yyyyyy_1[j];

                    tg_yyyzz_yyyyyyz_0[j] = pb_y * tg_yyzz_yyyyyyz_0[j] + wp_y[j] * tg_yyzz_yyyyyyz_1[j] + fl1_fx * tg_yzz_yyyyyyz_0[j] - fl1_fx * fl1_fza * tg_yzz_yyyyyyz_1[j] + 3.0 * fl1_fxn * tg_yyzz_yyyyyz_1[j];

                    tg_yyyzz_yyyyyzz_0[j] = pb_y * tg_yyzz_yyyyyzz_0[j] + wp_y[j] * tg_yyzz_yyyyyzz_1[j] + fl1_fx * tg_yzz_yyyyyzz_0[j] - fl1_fx * fl1_fza * tg_yzz_yyyyyzz_1[j] + 2.5 * fl1_fxn * tg_yyzz_yyyyzz_1[j];

                    tg_yyyzz_yyyyzzz_0[j] = pb_y * tg_yyzz_yyyyzzz_0[j] + wp_y[j] * tg_yyzz_yyyyzzz_1[j] + fl1_fx * tg_yzz_yyyyzzz_0[j] - fl1_fx * fl1_fza * tg_yzz_yyyyzzz_1[j] + 2.0 * fl1_fxn * tg_yyzz_yyyzzz_1[j];

                    tg_yyyzz_yyyzzzz_0[j] = pb_y * tg_yyzz_yyyzzzz_0[j] + wp_y[j] * tg_yyzz_yyyzzzz_1[j] + fl1_fx * tg_yzz_yyyzzzz_0[j] - fl1_fx * fl1_fza * tg_yzz_yyyzzzz_1[j] + 1.5 * fl1_fxn * tg_yyzz_yyzzzz_1[j];

                    tg_yyyzz_yyzzzzz_0[j] = pb_y * tg_yyzz_yyzzzzz_0[j] + wp_y[j] * tg_yyzz_yyzzzzz_1[j] + fl1_fx * tg_yzz_yyzzzzz_0[j] - fl1_fx * fl1_fza * tg_yzz_yyzzzzz_1[j] + fl1_fxn * tg_yyzz_yzzzzz_1[j];

                    tg_yyyzz_yzzzzzz_0[j] = pb_y * tg_yyzz_yzzzzzz_0[j] + wp_y[j] * tg_yyzz_yzzzzzz_1[j] + fl1_fx * tg_yzz_yzzzzzz_0[j] - fl1_fx * fl1_fza * tg_yzz_yzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyzz_zzzzzz_1[j];

                    tg_yyyzz_zzzzzzz_0[j] = pb_y * tg_yyzz_zzzzzzz_0[j] + wp_y[j] * tg_yyzz_zzzzzzz_1[j] + fl1_fx * tg_yzz_zzzzzzz_0[j] - fl1_fx * fl1_fza * tg_yzz_zzzzzzz_1[j];

                    tg_yyzzz_xxxxxxx_0[j] = pb_y * tg_yzzz_xxxxxxx_0[j] + wp_y[j] * tg_yzzz_xxxxxxx_1[j] + 0.5 * fl1_fx * tg_zzz_xxxxxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxxxxx_1[j];

                    tg_yyzzz_xxxxxxy_0[j] = pb_y * tg_yzzz_xxxxxxy_0[j] + wp_y[j] * tg_yzzz_xxxxxxy_1[j] + 0.5 * fl1_fx * tg_zzz_xxxxxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxxxxy_1[j] + 0.5 * fl1_fxn * tg_yzzz_xxxxxx_1[j];

                    tg_yyzzz_xxxxxxz_0[j] = pb_y * tg_yzzz_xxxxxxz_0[j] + wp_y[j] * tg_yzzz_xxxxxxz_1[j] + 0.5 * fl1_fx * tg_zzz_xxxxxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxxxxz_1[j];

                    tg_yyzzz_xxxxxyy_0[j] = pb_y * tg_yzzz_xxxxxyy_0[j] + wp_y[j] * tg_yzzz_xxxxxyy_1[j] + 0.5 * fl1_fx * tg_zzz_xxxxxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxxxyy_1[j] + fl1_fxn * tg_yzzz_xxxxxy_1[j];

                    tg_yyzzz_xxxxxyz_0[j] = pb_y * tg_yzzz_xxxxxyz_0[j] + wp_y[j] * tg_yzzz_xxxxxyz_1[j] + 0.5 * fl1_fx * tg_zzz_xxxxxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxxxyz_1[j] + 0.5 * fl1_fxn * tg_yzzz_xxxxxz_1[j];

                    tg_yyzzz_xxxxxzz_0[j] = pb_y * tg_yzzz_xxxxxzz_0[j] + wp_y[j] * tg_yzzz_xxxxxzz_1[j] + 0.5 * fl1_fx * tg_zzz_xxxxxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxxxzz_1[j];

                    tg_yyzzz_xxxxyyy_0[j] = pb_y * tg_yzzz_xxxxyyy_0[j] + wp_y[j] * tg_yzzz_xxxxyyy_1[j] + 0.5 * fl1_fx * tg_zzz_xxxxyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxxyyy_1[j] + 1.5 * fl1_fxn * tg_yzzz_xxxxyy_1[j];

                    tg_yyzzz_xxxxyyz_0[j] = pb_y * tg_yzzz_xxxxyyz_0[j] + wp_y[j] * tg_yzzz_xxxxyyz_1[j] + 0.5 * fl1_fx * tg_zzz_xxxxyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxxyyz_1[j] + fl1_fxn * tg_yzzz_xxxxyz_1[j];

                    tg_yyzzz_xxxxyzz_0[j] = pb_y * tg_yzzz_xxxxyzz_0[j] + wp_y[j] * tg_yzzz_xxxxyzz_1[j] + 0.5 * fl1_fx * tg_zzz_xxxxyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxxyzz_1[j] + 0.5 * fl1_fxn * tg_yzzz_xxxxzz_1[j];

                    tg_yyzzz_xxxxzzz_0[j] = pb_y * tg_yzzz_xxxxzzz_0[j] + wp_y[j] * tg_yzzz_xxxxzzz_1[j] + 0.5 * fl1_fx * tg_zzz_xxxxzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxxzzz_1[j];

                    tg_yyzzz_xxxyyyy_0[j] = pb_y * tg_yzzz_xxxyyyy_0[j] + wp_y[j] * tg_yzzz_xxxyyyy_1[j] + 0.5 * fl1_fx * tg_zzz_xxxyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxyyyy_1[j] + 2.0 * fl1_fxn * tg_yzzz_xxxyyy_1[j];

                    tg_yyzzz_xxxyyyz_0[j] = pb_y * tg_yzzz_xxxyyyz_0[j] + wp_y[j] * tg_yzzz_xxxyyyz_1[j] + 0.5 * fl1_fx * tg_zzz_xxxyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_yzzz_xxxyyz_1[j];

                    tg_yyzzz_xxxyyzz_0[j] = pb_y * tg_yzzz_xxxyyzz_0[j] + wp_y[j] * tg_yzzz_xxxyyzz_1[j] + 0.5 * fl1_fx * tg_zzz_xxxyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxyyzz_1[j] + fl1_fxn * tg_yzzz_xxxyzz_1[j];

                    tg_yyzzz_xxxyzzz_0[j] = pb_y * tg_yzzz_xxxyzzz_0[j] + wp_y[j] * tg_yzzz_xxxyzzz_1[j] + 0.5 * fl1_fx * tg_zzz_xxxyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxyzzz_1[j] + 0.5 * fl1_fxn * tg_yzzz_xxxzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSK_662_756(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (662,756)

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
                                             {5, -1, -1, -1},
                                             {7, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_7_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_4_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_3_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_4_6_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {6, -1, -1, -1}, 
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

                auto tg_yzzz_xxxzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 482); 

                auto tg_yzzz_xxyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 483); 

                auto tg_yzzz_xxyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 484); 

                auto tg_yzzz_xxyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 485); 

                auto tg_yzzz_xxyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 486); 

                auto tg_yzzz_xxyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 487); 

                auto tg_yzzz_xxzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 488); 

                auto tg_yzzz_xyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 489); 

                auto tg_yzzz_xyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 490); 

                auto tg_yzzz_xyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 491); 

                auto tg_yzzz_xyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 492); 

                auto tg_yzzz_xyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 493); 

                auto tg_yzzz_xyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 494); 

                auto tg_yzzz_xzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 495); 

                auto tg_yzzz_yyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 496); 

                auto tg_yzzz_yyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 497); 

                auto tg_yzzz_yyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 498); 

                auto tg_yzzz_yyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 499); 

                auto tg_yzzz_yyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 500); 

                auto tg_yzzz_yyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 501); 

                auto tg_yzzz_yzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 502); 

                auto tg_yzzz_zzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 503); 

                auto tg_zzzz_xxxxxxx_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 504); 

                auto tg_zzzz_xxxxxxy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 505); 

                auto tg_zzzz_xxxxxxz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 506); 

                auto tg_zzzz_xxxxxyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 507); 

                auto tg_zzzz_xxxxxyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 508); 

                auto tg_zzzz_xxxxxzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 509); 

                auto tg_zzzz_xxxxyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 510); 

                auto tg_zzzz_xxxxyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 511); 

                auto tg_zzzz_xxxxyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 512); 

                auto tg_zzzz_xxxxzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 513); 

                auto tg_zzzz_xxxyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 514); 

                auto tg_zzzz_xxxyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 515); 

                auto tg_zzzz_xxxyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 516); 

                auto tg_zzzz_xxxyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 517); 

                auto tg_zzzz_xxxzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 518); 

                auto tg_zzzz_xxyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 519); 

                auto tg_zzzz_xxyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 520); 

                auto tg_zzzz_xxyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 521); 

                auto tg_zzzz_xxyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 522); 

                auto tg_zzzz_xxyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 523); 

                auto tg_zzzz_xxzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 524); 

                auto tg_zzzz_xyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 525); 

                auto tg_zzzz_xyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 526); 

                auto tg_zzzz_xyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 527); 

                auto tg_zzzz_xyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 528); 

                auto tg_zzzz_xyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 529); 

                auto tg_zzzz_xyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 530); 

                auto tg_zzzz_xzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 531); 

                auto tg_zzzz_yyyyyyy_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 532); 

                auto tg_zzzz_yyyyyyz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 533); 

                auto tg_zzzz_yyyyyzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 534); 

                auto tg_zzzz_yyyyzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 535); 

                auto tg_zzzz_yyyzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 536); 

                auto tg_zzzz_yyzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 537); 

                auto tg_zzzz_yzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 538); 

                auto tg_zzzz_zzzzzzz_0 = primBuffer.data(pidx_g_4_7_m0 + 540 * idx + 539); 

                auto tg_yzzz_xxxzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 482); 

                auto tg_yzzz_xxyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 483); 

                auto tg_yzzz_xxyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 484); 

                auto tg_yzzz_xxyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 485); 

                auto tg_yzzz_xxyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 486); 

                auto tg_yzzz_xxyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 487); 

                auto tg_yzzz_xxzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 488); 

                auto tg_yzzz_xyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 489); 

                auto tg_yzzz_xyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 490); 

                auto tg_yzzz_xyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 491); 

                auto tg_yzzz_xyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 492); 

                auto tg_yzzz_xyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 493); 

                auto tg_yzzz_xyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 494); 

                auto tg_yzzz_xzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 495); 

                auto tg_yzzz_yyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 496); 

                auto tg_yzzz_yyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 497); 

                auto tg_yzzz_yyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 498); 

                auto tg_yzzz_yyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 499); 

                auto tg_yzzz_yyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 500); 

                auto tg_yzzz_yyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 501); 

                auto tg_yzzz_yzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 502); 

                auto tg_yzzz_zzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 503); 

                auto tg_zzzz_xxxxxxx_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 504); 

                auto tg_zzzz_xxxxxxy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 505); 

                auto tg_zzzz_xxxxxxz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 506); 

                auto tg_zzzz_xxxxxyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 507); 

                auto tg_zzzz_xxxxxyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 508); 

                auto tg_zzzz_xxxxxzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 509); 

                auto tg_zzzz_xxxxyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 510); 

                auto tg_zzzz_xxxxyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 511); 

                auto tg_zzzz_xxxxyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 512); 

                auto tg_zzzz_xxxxzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 513); 

                auto tg_zzzz_xxxyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 514); 

                auto tg_zzzz_xxxyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 515); 

                auto tg_zzzz_xxxyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 516); 

                auto tg_zzzz_xxxyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 517); 

                auto tg_zzzz_xxxzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 518); 

                auto tg_zzzz_xxyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 519); 

                auto tg_zzzz_xxyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 520); 

                auto tg_zzzz_xxyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 521); 

                auto tg_zzzz_xxyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 522); 

                auto tg_zzzz_xxyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 523); 

                auto tg_zzzz_xxzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 524); 

                auto tg_zzzz_xyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 525); 

                auto tg_zzzz_xyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 526); 

                auto tg_zzzz_xyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 527); 

                auto tg_zzzz_xyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 528); 

                auto tg_zzzz_xyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 529); 

                auto tg_zzzz_xyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 530); 

                auto tg_zzzz_xzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 531); 

                auto tg_zzzz_yyyyyyy_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 532); 

                auto tg_zzzz_yyyyyyz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 533); 

                auto tg_zzzz_yyyyyzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 534); 

                auto tg_zzzz_yyyyzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 535); 

                auto tg_zzzz_yyyzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 536); 

                auto tg_zzzz_yyzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 537); 

                auto tg_zzzz_yzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 538); 

                auto tg_zzzz_zzzzzzz_1 = primBuffer.data(pidx_g_4_7_m1 + 540 * idx + 539); 

                auto tg_zzz_xxxxxxx_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 324); 

                auto tg_zzz_xxxxxxy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 325); 

                auto tg_zzz_xxxxxxz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 326); 

                auto tg_zzz_xxxxxyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 327); 

                auto tg_zzz_xxxxxyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 328); 

                auto tg_zzz_xxxxxzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 329); 

                auto tg_zzz_xxxxyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 330); 

                auto tg_zzz_xxxxyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 331); 

                auto tg_zzz_xxxxyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 332); 

                auto tg_zzz_xxxxzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 333); 

                auto tg_zzz_xxxyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 334); 

                auto tg_zzz_xxxyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 335); 

                auto tg_zzz_xxxyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 336); 

                auto tg_zzz_xxxyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 337); 

                auto tg_zzz_xxxzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 338); 

                auto tg_zzz_xxyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 339); 

                auto tg_zzz_xxyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 340); 

                auto tg_zzz_xxyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 341); 

                auto tg_zzz_xxyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 342); 

                auto tg_zzz_xxyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 343); 

                auto tg_zzz_xxzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 344); 

                auto tg_zzz_xyyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 345); 

                auto tg_zzz_xyyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 346); 

                auto tg_zzz_xyyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 347); 

                auto tg_zzz_xyyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 348); 

                auto tg_zzz_xyyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 349); 

                auto tg_zzz_xyzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 350); 

                auto tg_zzz_xzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 351); 

                auto tg_zzz_yyyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 352); 

                auto tg_zzz_yyyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 353); 

                auto tg_zzz_yyyyyzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 354); 

                auto tg_zzz_yyyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 355); 

                auto tg_zzz_yyyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 356); 

                auto tg_zzz_yyzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 357); 

                auto tg_zzz_yzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 358); 

                auto tg_zzz_zzzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 359); 

                auto tg_zzz_xxxxxxx_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 324); 

                auto tg_zzz_xxxxxxy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 325); 

                auto tg_zzz_xxxxxxz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 326); 

                auto tg_zzz_xxxxxyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 327); 

                auto tg_zzz_xxxxxyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 328); 

                auto tg_zzz_xxxxxzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 329); 

                auto tg_zzz_xxxxyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 330); 

                auto tg_zzz_xxxxyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 331); 

                auto tg_zzz_xxxxyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 332); 

                auto tg_zzz_xxxxzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 333); 

                auto tg_zzz_xxxyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 334); 

                auto tg_zzz_xxxyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 335); 

                auto tg_zzz_xxxyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 336); 

                auto tg_zzz_xxxyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 337); 

                auto tg_zzz_xxxzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 338); 

                auto tg_zzz_xxyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 339); 

                auto tg_zzz_xxyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 340); 

                auto tg_zzz_xxyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 341); 

                auto tg_zzz_xxyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 342); 

                auto tg_zzz_xxyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 343); 

                auto tg_zzz_xxzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 344); 

                auto tg_zzz_xyyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 345); 

                auto tg_zzz_xyyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 346); 

                auto tg_zzz_xyyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 347); 

                auto tg_zzz_xyyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 348); 

                auto tg_zzz_xyyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 349); 

                auto tg_zzz_xyzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 350); 

                auto tg_zzz_xzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 351); 

                auto tg_zzz_yyyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 352); 

                auto tg_zzz_yyyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 353); 

                auto tg_zzz_yyyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 354); 

                auto tg_zzz_yyyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 355); 

                auto tg_zzz_yyyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 356); 

                auto tg_zzz_yyzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 357); 

                auto tg_zzz_yzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 358); 

                auto tg_zzz_zzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 359); 

                auto tg_yzzz_xxyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 374); 

                auto tg_yzzz_xxyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 375); 

                auto tg_yzzz_xxyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 376); 

                auto tg_yzzz_xxyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 377); 

                auto tg_yzzz_xxzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 378); 

                auto tg_yzzz_xyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 379); 

                auto tg_yzzz_xyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 380); 

                auto tg_yzzz_xyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 381); 

                auto tg_yzzz_xyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 382); 

                auto tg_yzzz_xyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 383); 

                auto tg_yzzz_xzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 384); 

                auto tg_yzzz_yyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 385); 

                auto tg_yzzz_yyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 386); 

                auto tg_yzzz_yyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 387); 

                auto tg_yzzz_yyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 388); 

                auto tg_yzzz_yyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 389); 

                auto tg_yzzz_yzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 390); 

                auto tg_yzzz_zzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 391); 

                auto tg_zzzz_xxxxxx_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 392); 

                auto tg_zzzz_xxxxxy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 393); 

                auto tg_zzzz_xxxxxz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 394); 

                auto tg_zzzz_xxxxyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 395); 

                auto tg_zzzz_xxxxyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 396); 

                auto tg_zzzz_xxxxzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 397); 

                auto tg_zzzz_xxxyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 398); 

                auto tg_zzzz_xxxyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 399); 

                auto tg_zzzz_xxxyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 400); 

                auto tg_zzzz_xxxzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 401); 

                auto tg_zzzz_xxyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 402); 

                auto tg_zzzz_xxyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 403); 

                auto tg_zzzz_xxyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 404); 

                auto tg_zzzz_xxyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 405); 

                auto tg_zzzz_xxzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 406); 

                auto tg_zzzz_xyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 407); 

                auto tg_zzzz_xyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 408); 

                auto tg_zzzz_xyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 409); 

                auto tg_zzzz_xyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 410); 

                auto tg_zzzz_xyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 411); 

                auto tg_zzzz_xzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 412); 

                auto tg_zzzz_yyyyyy_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 413); 

                auto tg_zzzz_yyyyyz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 414); 

                auto tg_zzzz_yyyyzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 415); 

                auto tg_zzzz_yyyzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 416); 

                auto tg_zzzz_yyzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 417); 

                auto tg_zzzz_yzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 418); 

                auto tg_zzzz_zzzzzz_1 = primBuffer.data(pidx_g_4_6_m1 + 420 * idx + 419); 

                // set up pointers to integrals

                auto tg_yyzzz_xxxzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 662); 

                auto tg_yyzzz_xxyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 663); 

                auto tg_yyzzz_xxyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 664); 

                auto tg_yyzzz_xxyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 665); 

                auto tg_yyzzz_xxyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 666); 

                auto tg_yyzzz_xxyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 667); 

                auto tg_yyzzz_xxzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 668); 

                auto tg_yyzzz_xyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 669); 

                auto tg_yyzzz_xyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 670); 

                auto tg_yyzzz_xyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 671); 

                auto tg_yyzzz_xyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 672); 

                auto tg_yyzzz_xyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 673); 

                auto tg_yyzzz_xyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 674); 

                auto tg_yyzzz_xzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 675); 

                auto tg_yyzzz_yyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 676); 

                auto tg_yyzzz_yyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 677); 

                auto tg_yyzzz_yyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 678); 

                auto tg_yyzzz_yyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 679); 

                auto tg_yyzzz_yyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 680); 

                auto tg_yyzzz_yyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 681); 

                auto tg_yyzzz_yzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 682); 

                auto tg_yyzzz_zzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 683); 

                auto tg_yzzzz_xxxxxxx_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 684); 

                auto tg_yzzzz_xxxxxxy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 685); 

                auto tg_yzzzz_xxxxxxz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 686); 

                auto tg_yzzzz_xxxxxyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 687); 

                auto tg_yzzzz_xxxxxyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 688); 

                auto tg_yzzzz_xxxxxzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 689); 

                auto tg_yzzzz_xxxxyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 690); 

                auto tg_yzzzz_xxxxyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 691); 

                auto tg_yzzzz_xxxxyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 692); 

                auto tg_yzzzz_xxxxzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 693); 

                auto tg_yzzzz_xxxyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 694); 

                auto tg_yzzzz_xxxyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 695); 

                auto tg_yzzzz_xxxyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 696); 

                auto tg_yzzzz_xxxyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 697); 

                auto tg_yzzzz_xxxzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 698); 

                auto tg_yzzzz_xxyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 699); 

                auto tg_yzzzz_xxyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 700); 

                auto tg_yzzzz_xxyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 701); 

                auto tg_yzzzz_xxyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 702); 

                auto tg_yzzzz_xxyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 703); 

                auto tg_yzzzz_xxzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 704); 

                auto tg_yzzzz_xyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 705); 

                auto tg_yzzzz_xyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 706); 

                auto tg_yzzzz_xyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 707); 

                auto tg_yzzzz_xyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 708); 

                auto tg_yzzzz_xyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 709); 

                auto tg_yzzzz_xyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 710); 

                auto tg_yzzzz_xzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 711); 

                auto tg_yzzzz_yyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 712); 

                auto tg_yzzzz_yyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 713); 

                auto tg_yzzzz_yyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 714); 

                auto tg_yzzzz_yyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 715); 

                auto tg_yzzzz_yyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 716); 

                auto tg_yzzzz_yyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 717); 

                auto tg_yzzzz_yzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 718); 

                auto tg_yzzzz_zzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 719); 

                auto tg_zzzzz_xxxxxxx_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 720); 

                auto tg_zzzzz_xxxxxxy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 721); 

                auto tg_zzzzz_xxxxxxz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 722); 

                auto tg_zzzzz_xxxxxyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 723); 

                auto tg_zzzzz_xxxxxyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 724); 

                auto tg_zzzzz_xxxxxzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 725); 

                auto tg_zzzzz_xxxxyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 726); 

                auto tg_zzzzz_xxxxyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 727); 

                auto tg_zzzzz_xxxxyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 728); 

                auto tg_zzzzz_xxxxzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 729); 

                auto tg_zzzzz_xxxyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 730); 

                auto tg_zzzzz_xxxyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 731); 

                auto tg_zzzzz_xxxyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 732); 

                auto tg_zzzzz_xxxyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 733); 

                auto tg_zzzzz_xxxzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 734); 

                auto tg_zzzzz_xxyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 735); 

                auto tg_zzzzz_xxyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 736); 

                auto tg_zzzzz_xxyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 737); 

                auto tg_zzzzz_xxyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 738); 

                auto tg_zzzzz_xxyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 739); 

                auto tg_zzzzz_xxzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 740); 

                auto tg_zzzzz_xyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 741); 

                auto tg_zzzzz_xyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 742); 

                auto tg_zzzzz_xyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 743); 

                auto tg_zzzzz_xyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 744); 

                auto tg_zzzzz_xyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 745); 

                auto tg_zzzzz_xyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 746); 

                auto tg_zzzzz_xzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 747); 

                auto tg_zzzzz_yyyyyyy_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 748); 

                auto tg_zzzzz_yyyyyyz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 749); 

                auto tg_zzzzz_yyyyyzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 750); 

                auto tg_zzzzz_yyyyzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 751); 

                auto tg_zzzzz_yyyzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 752); 

                auto tg_zzzzz_yyzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 753); 

                auto tg_zzzzz_yzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 754); 

                auto tg_zzzzz_zzzzzzz_0 = primBuffer.data(pidx_g_5_7_m0 + 756 * idx + 755); 

                // Batch of Integrals (662,756)

                #pragma omp simd aligned(fxn, fza, tg_yyzzz_xxxzzzz_0, tg_yyzzz_xxyyyyy_0, tg_yyzzz_xxyyyyz_0, \
                                         tg_yyzzz_xxyyyzz_0, tg_yyzzz_xxyyzzz_0, tg_yyzzz_xxyzzzz_0, tg_yyzzz_xxzzzzz_0, \
                                         tg_yyzzz_xyyyyyy_0, tg_yyzzz_xyyyyyz_0, tg_yyzzz_xyyyyzz_0, tg_yyzzz_xyyyzzz_0, \
                                         tg_yyzzz_xyyzzzz_0, tg_yyzzz_xyzzzzz_0, tg_yyzzz_xzzzzzz_0, tg_yyzzz_yyyyyyy_0, \
                                         tg_yyzzz_yyyyyyz_0, tg_yyzzz_yyyyyzz_0, tg_yyzzz_yyyyzzz_0, tg_yyzzz_yyyzzzz_0, \
                                         tg_yyzzz_yyzzzzz_0, tg_yyzzz_yzzzzzz_0, tg_yyzzz_zzzzzzz_0, tg_yzzz_xxxzzzz_0, \
                                         tg_yzzz_xxxzzzz_1, tg_yzzz_xxyyyy_1, tg_yzzz_xxyyyyy_0, tg_yzzz_xxyyyyy_1, \
                                         tg_yzzz_xxyyyyz_0, tg_yzzz_xxyyyyz_1, tg_yzzz_xxyyyz_1, tg_yzzz_xxyyyzz_0, \
                                         tg_yzzz_xxyyyzz_1, tg_yzzz_xxyyzz_1, tg_yzzz_xxyyzzz_0, tg_yzzz_xxyyzzz_1, \
                                         tg_yzzz_xxyzzz_1, tg_yzzz_xxyzzzz_0, tg_yzzz_xxyzzzz_1, tg_yzzz_xxzzzz_1, \
                                         tg_yzzz_xxzzzzz_0, tg_yzzz_xxzzzzz_1, tg_yzzz_xyyyyy_1, tg_yzzz_xyyyyyy_0, \
                                         tg_yzzz_xyyyyyy_1, tg_yzzz_xyyyyyz_0, tg_yzzz_xyyyyyz_1, tg_yzzz_xyyyyz_1, \
                                         tg_yzzz_xyyyyzz_0, tg_yzzz_xyyyyzz_1, tg_yzzz_xyyyzz_1, tg_yzzz_xyyyzzz_0, \
                                         tg_yzzz_xyyyzzz_1, tg_yzzz_xyyzzz_1, tg_yzzz_xyyzzzz_0, tg_yzzz_xyyzzzz_1, \
                                         tg_yzzz_xyzzzz_1, tg_yzzz_xyzzzzz_0, tg_yzzz_xyzzzzz_1, tg_yzzz_xzzzzz_1, \
                                         tg_yzzz_xzzzzzz_0, tg_yzzz_xzzzzzz_1, tg_yzzz_yyyyyy_1, tg_yzzz_yyyyyyy_0, \
                                         tg_yzzz_yyyyyyy_1, tg_yzzz_yyyyyyz_0, tg_yzzz_yyyyyyz_1, tg_yzzz_yyyyyz_1, \
                                         tg_yzzz_yyyyyzz_0, tg_yzzz_yyyyyzz_1, tg_yzzz_yyyyzz_1, tg_yzzz_yyyyzzz_0, \
                                         tg_yzzz_yyyyzzz_1, tg_yzzz_yyyzzz_1, tg_yzzz_yyyzzzz_0, tg_yzzz_yyyzzzz_1, \
                                         tg_yzzz_yyzzzz_1, tg_yzzz_yyzzzzz_0, tg_yzzz_yyzzzzz_1, tg_yzzz_yzzzzz_1, \
                                         tg_yzzz_yzzzzzz_0, tg_yzzz_yzzzzzz_1, tg_yzzz_zzzzzz_1, tg_yzzz_zzzzzzz_0, \
                                         tg_yzzz_zzzzzzz_1, tg_yzzzz_xxxxxxx_0, tg_yzzzz_xxxxxxy_0, tg_yzzzz_xxxxxxz_0, \
                                         tg_yzzzz_xxxxxyy_0, tg_yzzzz_xxxxxyz_0, tg_yzzzz_xxxxxzz_0, tg_yzzzz_xxxxyyy_0, \
                                         tg_yzzzz_xxxxyyz_0, tg_yzzzz_xxxxyzz_0, tg_yzzzz_xxxxzzz_0, tg_yzzzz_xxxyyyy_0, \
                                         tg_yzzzz_xxxyyyz_0, tg_yzzzz_xxxyyzz_0, tg_yzzzz_xxxyzzz_0, tg_yzzzz_xxxzzzz_0, \
                                         tg_yzzzz_xxyyyyy_0, tg_yzzzz_xxyyyyz_0, tg_yzzzz_xxyyyzz_0, tg_yzzzz_xxyyzzz_0, \
                                         tg_yzzzz_xxyzzzz_0, tg_yzzzz_xxzzzzz_0, tg_yzzzz_xyyyyyy_0, tg_yzzzz_xyyyyyz_0, \
                                         tg_yzzzz_xyyyyzz_0, tg_yzzzz_xyyyzzz_0, tg_yzzzz_xyyzzzz_0, tg_yzzzz_xyzzzzz_0, \
                                         tg_yzzzz_xzzzzzz_0, tg_yzzzz_yyyyyyy_0, tg_yzzzz_yyyyyyz_0, tg_yzzzz_yyyyyzz_0, \
                                         tg_yzzzz_yyyyzzz_0, tg_yzzzz_yyyzzzz_0, tg_yzzzz_yyzzzzz_0, tg_yzzzz_yzzzzzz_0, \
                                         tg_yzzzz_zzzzzzz_0, tg_zzz_xxxxxxx_0, tg_zzz_xxxxxxx_1, tg_zzz_xxxxxxy_0, \
                                         tg_zzz_xxxxxxy_1, tg_zzz_xxxxxxz_0, tg_zzz_xxxxxxz_1, tg_zzz_xxxxxyy_0, \
                                         tg_zzz_xxxxxyy_1, tg_zzz_xxxxxyz_0, tg_zzz_xxxxxyz_1, tg_zzz_xxxxxzz_0, \
                                         tg_zzz_xxxxxzz_1, tg_zzz_xxxxyyy_0, tg_zzz_xxxxyyy_1, tg_zzz_xxxxyyz_0, \
                                         tg_zzz_xxxxyyz_1, tg_zzz_xxxxyzz_0, tg_zzz_xxxxyzz_1, tg_zzz_xxxxzzz_0, \
                                         tg_zzz_xxxxzzz_1, tg_zzz_xxxyyyy_0, tg_zzz_xxxyyyy_1, tg_zzz_xxxyyyz_0, \
                                         tg_zzz_xxxyyyz_1, tg_zzz_xxxyyzz_0, tg_zzz_xxxyyzz_1, tg_zzz_xxxyzzz_0, \
                                         tg_zzz_xxxyzzz_1, tg_zzz_xxxzzzz_0, tg_zzz_xxxzzzz_1, tg_zzz_xxyyyyy_0, \
                                         tg_zzz_xxyyyyy_1, tg_zzz_xxyyyyz_0, tg_zzz_xxyyyyz_1, tg_zzz_xxyyyzz_0, \
                                         tg_zzz_xxyyyzz_1, tg_zzz_xxyyzzz_0, tg_zzz_xxyyzzz_1, tg_zzz_xxyzzzz_0, \
                                         tg_zzz_xxyzzzz_1, tg_zzz_xxzzzzz_0, tg_zzz_xxzzzzz_1, tg_zzz_xyyyyyy_0, \
                                         tg_zzz_xyyyyyy_1, tg_zzz_xyyyyyz_0, tg_zzz_xyyyyyz_1, tg_zzz_xyyyyzz_0, \
                                         tg_zzz_xyyyyzz_1, tg_zzz_xyyyzzz_0, tg_zzz_xyyyzzz_1, tg_zzz_xyyzzzz_0, \
                                         tg_zzz_xyyzzzz_1, tg_zzz_xyzzzzz_0, tg_zzz_xyzzzzz_1, tg_zzz_xzzzzzz_0, \
                                         tg_zzz_xzzzzzz_1, tg_zzz_yyyyyyy_0, tg_zzz_yyyyyyy_1, tg_zzz_yyyyyyz_0, \
                                         tg_zzz_yyyyyyz_1, tg_zzz_yyyyyzz_0, tg_zzz_yyyyyzz_1, tg_zzz_yyyyzzz_0, \
                                         tg_zzz_yyyyzzz_1, tg_zzz_yyyzzzz_0, tg_zzz_yyyzzzz_1, tg_zzz_yyzzzzz_0, \
                                         tg_zzz_yyzzzzz_1, tg_zzz_yzzzzzz_0, tg_zzz_yzzzzzz_1, tg_zzz_zzzzzzz_0, \
                                         tg_zzz_zzzzzzz_1, tg_zzzz_xxxxxx_1, tg_zzzz_xxxxxxx_0, tg_zzzz_xxxxxxx_1, \
                                         tg_zzzz_xxxxxxy_0, tg_zzzz_xxxxxxy_1, tg_zzzz_xxxxxxz_0, tg_zzzz_xxxxxxz_1, \
                                         tg_zzzz_xxxxxy_1, tg_zzzz_xxxxxyy_0, tg_zzzz_xxxxxyy_1, tg_zzzz_xxxxxyz_0, \
                                         tg_zzzz_xxxxxyz_1, tg_zzzz_xxxxxz_1, tg_zzzz_xxxxxzz_0, tg_zzzz_xxxxxzz_1, \
                                         tg_zzzz_xxxxyy_1, tg_zzzz_xxxxyyy_0, tg_zzzz_xxxxyyy_1, tg_zzzz_xxxxyyz_0, \
                                         tg_zzzz_xxxxyyz_1, tg_zzzz_xxxxyz_1, tg_zzzz_xxxxyzz_0, tg_zzzz_xxxxyzz_1, \
                                         tg_zzzz_xxxxzz_1, tg_zzzz_xxxxzzz_0, tg_zzzz_xxxxzzz_1, tg_zzzz_xxxyyy_1, \
                                         tg_zzzz_xxxyyyy_0, tg_zzzz_xxxyyyy_1, tg_zzzz_xxxyyyz_0, tg_zzzz_xxxyyyz_1, \
                                         tg_zzzz_xxxyyz_1, tg_zzzz_xxxyyzz_0, tg_zzzz_xxxyyzz_1, tg_zzzz_xxxyzz_1, \
                                         tg_zzzz_xxxyzzz_0, tg_zzzz_xxxyzzz_1, tg_zzzz_xxxzzz_1, tg_zzzz_xxxzzzz_0, \
                                         tg_zzzz_xxxzzzz_1, tg_zzzz_xxyyyy_1, tg_zzzz_xxyyyyy_0, tg_zzzz_xxyyyyy_1, \
                                         tg_zzzz_xxyyyyz_0, tg_zzzz_xxyyyyz_1, tg_zzzz_xxyyyz_1, tg_zzzz_xxyyyzz_0, \
                                         tg_zzzz_xxyyyzz_1, tg_zzzz_xxyyzz_1, tg_zzzz_xxyyzzz_0, tg_zzzz_xxyyzzz_1, \
                                         tg_zzzz_xxyzzz_1, tg_zzzz_xxyzzzz_0, tg_zzzz_xxyzzzz_1, tg_zzzz_xxzzzz_1, \
                                         tg_zzzz_xxzzzzz_0, tg_zzzz_xxzzzzz_1, tg_zzzz_xyyyyy_1, tg_zzzz_xyyyyyy_0, \
                                         tg_zzzz_xyyyyyy_1, tg_zzzz_xyyyyyz_0, tg_zzzz_xyyyyyz_1, tg_zzzz_xyyyyz_1, \
                                         tg_zzzz_xyyyyzz_0, tg_zzzz_xyyyyzz_1, tg_zzzz_xyyyzz_1, tg_zzzz_xyyyzzz_0, \
                                         tg_zzzz_xyyyzzz_1, tg_zzzz_xyyzzz_1, tg_zzzz_xyyzzzz_0, tg_zzzz_xyyzzzz_1, \
                                         tg_zzzz_xyzzzz_1, tg_zzzz_xyzzzzz_0, tg_zzzz_xyzzzzz_1, tg_zzzz_xzzzzz_1, \
                                         tg_zzzz_xzzzzzz_0, tg_zzzz_xzzzzzz_1, tg_zzzz_yyyyyy_1, tg_zzzz_yyyyyyy_0, \
                                         tg_zzzz_yyyyyyy_1, tg_zzzz_yyyyyyz_0, tg_zzzz_yyyyyyz_1, tg_zzzz_yyyyyz_1, \
                                         tg_zzzz_yyyyyzz_0, tg_zzzz_yyyyyzz_1, tg_zzzz_yyyyzz_1, tg_zzzz_yyyyzzz_0, \
                                         tg_zzzz_yyyyzzz_1, tg_zzzz_yyyzzz_1, tg_zzzz_yyyzzzz_0, tg_zzzz_yyyzzzz_1, \
                                         tg_zzzz_yyzzzz_1, tg_zzzz_yyzzzzz_0, tg_zzzz_yyzzzzz_1, tg_zzzz_yzzzzz_1, \
                                         tg_zzzz_yzzzzzz_0, tg_zzzz_yzzzzzz_1, tg_zzzz_zzzzzz_1, tg_zzzz_zzzzzzz_0, \
                                         tg_zzzz_zzzzzzz_1, tg_zzzzz_xxxxxxx_0, tg_zzzzz_xxxxxxy_0, tg_zzzzz_xxxxxxz_0, \
                                         tg_zzzzz_xxxxxyy_0, tg_zzzzz_xxxxxyz_0, tg_zzzzz_xxxxxzz_0, tg_zzzzz_xxxxyyy_0, \
                                         tg_zzzzz_xxxxyyz_0, tg_zzzzz_xxxxyzz_0, tg_zzzzz_xxxxzzz_0, tg_zzzzz_xxxyyyy_0, \
                                         tg_zzzzz_xxxyyyz_0, tg_zzzzz_xxxyyzz_0, tg_zzzzz_xxxyzzz_0, tg_zzzzz_xxxzzzz_0, \
                                         tg_zzzzz_xxyyyyy_0, tg_zzzzz_xxyyyyz_0, tg_zzzzz_xxyyyzz_0, tg_zzzzz_xxyyzzz_0, \
                                         tg_zzzzz_xxyzzzz_0, tg_zzzzz_xxzzzzz_0, tg_zzzzz_xyyyyyy_0, tg_zzzzz_xyyyyyz_0, \
                                         tg_zzzzz_xyyyyzz_0, tg_zzzzz_xyyyzzz_0, tg_zzzzz_xyyzzzz_0, tg_zzzzz_xyzzzzz_0, \
                                         tg_zzzzz_xzzzzzz_0, tg_zzzzz_yyyyyyy_0, tg_zzzzz_yyyyyyz_0, tg_zzzzz_yyyyyzz_0, \
                                         tg_zzzzz_yyyyzzz_0, tg_zzzzz_yyyzzzz_0, tg_zzzzz_yyzzzzz_0, tg_zzzzz_yzzzzzz_0, \
                                         tg_zzzzz_zzzzzzz_0, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_yyzzz_xxxzzzz_0[j] = pb_y * tg_yzzz_xxxzzzz_0[j] + wp_y[j] * tg_yzzz_xxxzzzz_1[j] + 0.5 * fl1_fx * tg_zzz_xxxzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxzzzz_1[j];

                    tg_yyzzz_xxyyyyy_0[j] = pb_y * tg_yzzz_xxyyyyy_0[j] + wp_y[j] * tg_yzzz_xxyyyyy_1[j] + 0.5 * fl1_fx * tg_zzz_xxyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxyyyyy_1[j] + 2.5 * fl1_fxn * tg_yzzz_xxyyyy_1[j];

                    tg_yyzzz_xxyyyyz_0[j] = pb_y * tg_yzzz_xxyyyyz_0[j] + wp_y[j] * tg_yzzz_xxyyyyz_1[j] + 0.5 * fl1_fx * tg_zzz_xxyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxyyyyz_1[j] + 2.0 * fl1_fxn * tg_yzzz_xxyyyz_1[j];

                    tg_yyzzz_xxyyyzz_0[j] = pb_y * tg_yzzz_xxyyyzz_0[j] + wp_y[j] * tg_yzzz_xxyyyzz_1[j] + 0.5 * fl1_fx * tg_zzz_xxyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxyyyzz_1[j] + 1.5 * fl1_fxn * tg_yzzz_xxyyzz_1[j];

                    tg_yyzzz_xxyyzzz_0[j] = pb_y * tg_yzzz_xxyyzzz_0[j] + wp_y[j] * tg_yzzz_xxyyzzz_1[j] + 0.5 * fl1_fx * tg_zzz_xxyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxyyzzz_1[j] + fl1_fxn * tg_yzzz_xxyzzz_1[j];

                    tg_yyzzz_xxyzzzz_0[j] = pb_y * tg_yzzz_xxyzzzz_0[j] + wp_y[j] * tg_yzzz_xxyzzzz_1[j] + 0.5 * fl1_fx * tg_zzz_xxyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxyzzzz_1[j] + 0.5 * fl1_fxn * tg_yzzz_xxzzzz_1[j];

                    tg_yyzzz_xxzzzzz_0[j] = pb_y * tg_yzzz_xxzzzzz_0[j] + wp_y[j] * tg_yzzz_xxzzzzz_1[j] + 0.5 * fl1_fx * tg_zzz_xxzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxzzzzz_1[j];

                    tg_yyzzz_xyyyyyy_0[j] = pb_y * tg_yzzz_xyyyyyy_0[j] + wp_y[j] * tg_yzzz_xyyyyyy_1[j] + 0.5 * fl1_fx * tg_zzz_xyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xyyyyyy_1[j] + 3.0 * fl1_fxn * tg_yzzz_xyyyyy_1[j];

                    tg_yyzzz_xyyyyyz_0[j] = pb_y * tg_yzzz_xyyyyyz_0[j] + wp_y[j] * tg_yzzz_xyyyyyz_1[j] + 0.5 * fl1_fx * tg_zzz_xyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xyyyyyz_1[j] + 2.5 * fl1_fxn * tg_yzzz_xyyyyz_1[j];

                    tg_yyzzz_xyyyyzz_0[j] = pb_y * tg_yzzz_xyyyyzz_0[j] + wp_y[j] * tg_yzzz_xyyyyzz_1[j] + 0.5 * fl1_fx * tg_zzz_xyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xyyyyzz_1[j] + 2.0 * fl1_fxn * tg_yzzz_xyyyzz_1[j];

                    tg_yyzzz_xyyyzzz_0[j] = pb_y * tg_yzzz_xyyyzzz_0[j] + wp_y[j] * tg_yzzz_xyyyzzz_1[j] + 0.5 * fl1_fx * tg_zzz_xyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xyyyzzz_1[j] + 1.5 * fl1_fxn * tg_yzzz_xyyzzz_1[j];

                    tg_yyzzz_xyyzzzz_0[j] = pb_y * tg_yzzz_xyyzzzz_0[j] + wp_y[j] * tg_yzzz_xyyzzzz_1[j] + 0.5 * fl1_fx * tg_zzz_xyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xyyzzzz_1[j] + fl1_fxn * tg_yzzz_xyzzzz_1[j];

                    tg_yyzzz_xyzzzzz_0[j] = pb_y * tg_yzzz_xyzzzzz_0[j] + wp_y[j] * tg_yzzz_xyzzzzz_1[j] + 0.5 * fl1_fx * tg_zzz_xyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_yzzz_xzzzzz_1[j];

                    tg_yyzzz_xzzzzzz_0[j] = pb_y * tg_yzzz_xzzzzzz_0[j] + wp_y[j] * tg_yzzz_xzzzzzz_1[j] + 0.5 * fl1_fx * tg_zzz_xzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xzzzzzz_1[j];

                    tg_yyzzz_yyyyyyy_0[j] = pb_y * tg_yzzz_yyyyyyy_0[j] + wp_y[j] * tg_yzzz_yyyyyyy_1[j] + 0.5 * fl1_fx * tg_zzz_yyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_yyyyyyy_1[j] + 3.5 * fl1_fxn * tg_yzzz_yyyyyy_1[j];

                    tg_yyzzz_yyyyyyz_0[j] = pb_y * tg_yzzz_yyyyyyz_0[j] + wp_y[j] * tg_yzzz_yyyyyyz_1[j] + 0.5 * fl1_fx * tg_zzz_yyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_yyyyyyz_1[j] + 3.0 * fl1_fxn * tg_yzzz_yyyyyz_1[j];

                    tg_yyzzz_yyyyyzz_0[j] = pb_y * tg_yzzz_yyyyyzz_0[j] + wp_y[j] * tg_yzzz_yyyyyzz_1[j] + 0.5 * fl1_fx * tg_zzz_yyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_yyyyyzz_1[j] + 2.5 * fl1_fxn * tg_yzzz_yyyyzz_1[j];

                    tg_yyzzz_yyyyzzz_0[j] = pb_y * tg_yzzz_yyyyzzz_0[j] + wp_y[j] * tg_yzzz_yyyyzzz_1[j] + 0.5 * fl1_fx * tg_zzz_yyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_yyyyzzz_1[j] + 2.0 * fl1_fxn * tg_yzzz_yyyzzz_1[j];

                    tg_yyzzz_yyyzzzz_0[j] = pb_y * tg_yzzz_yyyzzzz_0[j] + wp_y[j] * tg_yzzz_yyyzzzz_1[j] + 0.5 * fl1_fx * tg_zzz_yyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_yyyzzzz_1[j] + 1.5 * fl1_fxn * tg_yzzz_yyzzzz_1[j];

                    tg_yyzzz_yyzzzzz_0[j] = pb_y * tg_yzzz_yyzzzzz_0[j] + wp_y[j] * tg_yzzz_yyzzzzz_1[j] + 0.5 * fl1_fx * tg_zzz_yyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_yyzzzzz_1[j] + fl1_fxn * tg_yzzz_yzzzzz_1[j];

                    tg_yyzzz_yzzzzzz_0[j] = pb_y * tg_yzzz_yzzzzzz_0[j] + wp_y[j] * tg_yzzz_yzzzzzz_1[j] + 0.5 * fl1_fx * tg_zzz_yzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_yzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yzzz_zzzzzz_1[j];

                    tg_yyzzz_zzzzzzz_0[j] = pb_y * tg_yzzz_zzzzzzz_0[j] + wp_y[j] * tg_yzzz_zzzzzzz_1[j] + 0.5 * fl1_fx * tg_zzz_zzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_zzzzzzz_1[j];

                    tg_yzzzz_xxxxxxx_0[j] = pb_y * tg_zzzz_xxxxxxx_0[j] + wp_y[j] * tg_zzzz_xxxxxxx_1[j];

                    tg_yzzzz_xxxxxxy_0[j] = pb_y * tg_zzzz_xxxxxxy_0[j] + wp_y[j] * tg_zzzz_xxxxxxy_1[j] + 0.5 * fl1_fxn * tg_zzzz_xxxxxx_1[j];

                    tg_yzzzz_xxxxxxz_0[j] = pb_y * tg_zzzz_xxxxxxz_0[j] + wp_y[j] * tg_zzzz_xxxxxxz_1[j];

                    tg_yzzzz_xxxxxyy_0[j] = pb_y * tg_zzzz_xxxxxyy_0[j] + wp_y[j] * tg_zzzz_xxxxxyy_1[j] + fl1_fxn * tg_zzzz_xxxxxy_1[j];

                    tg_yzzzz_xxxxxyz_0[j] = pb_y * tg_zzzz_xxxxxyz_0[j] + wp_y[j] * tg_zzzz_xxxxxyz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xxxxxz_1[j];

                    tg_yzzzz_xxxxxzz_0[j] = pb_y * tg_zzzz_xxxxxzz_0[j] + wp_y[j] * tg_zzzz_xxxxxzz_1[j];

                    tg_yzzzz_xxxxyyy_0[j] = pb_y * tg_zzzz_xxxxyyy_0[j] + wp_y[j] * tg_zzzz_xxxxyyy_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxxxyy_1[j];

                    tg_yzzzz_xxxxyyz_0[j] = pb_y * tg_zzzz_xxxxyyz_0[j] + wp_y[j] * tg_zzzz_xxxxyyz_1[j] + fl1_fxn * tg_zzzz_xxxxyz_1[j];

                    tg_yzzzz_xxxxyzz_0[j] = pb_y * tg_zzzz_xxxxyzz_0[j] + wp_y[j] * tg_zzzz_xxxxyzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xxxxzz_1[j];

                    tg_yzzzz_xxxxzzz_0[j] = pb_y * tg_zzzz_xxxxzzz_0[j] + wp_y[j] * tg_zzzz_xxxxzzz_1[j];

                    tg_yzzzz_xxxyyyy_0[j] = pb_y * tg_zzzz_xxxyyyy_0[j] + wp_y[j] * tg_zzzz_xxxyyyy_1[j] + 2.0 * fl1_fxn * tg_zzzz_xxxyyy_1[j];

                    tg_yzzzz_xxxyyyz_0[j] = pb_y * tg_zzzz_xxxyyyz_0[j] + wp_y[j] * tg_zzzz_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxxyyz_1[j];

                    tg_yzzzz_xxxyyzz_0[j] = pb_y * tg_zzzz_xxxyyzz_0[j] + wp_y[j] * tg_zzzz_xxxyyzz_1[j] + fl1_fxn * tg_zzzz_xxxyzz_1[j];

                    tg_yzzzz_xxxyzzz_0[j] = pb_y * tg_zzzz_xxxyzzz_0[j] + wp_y[j] * tg_zzzz_xxxyzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xxxzzz_1[j];

                    tg_yzzzz_xxxzzzz_0[j] = pb_y * tg_zzzz_xxxzzzz_0[j] + wp_y[j] * tg_zzzz_xxxzzzz_1[j];

                    tg_yzzzz_xxyyyyy_0[j] = pb_y * tg_zzzz_xxyyyyy_0[j] + wp_y[j] * tg_zzzz_xxyyyyy_1[j] + 2.5 * fl1_fxn * tg_zzzz_xxyyyy_1[j];

                    tg_yzzzz_xxyyyyz_0[j] = pb_y * tg_zzzz_xxyyyyz_0[j] + wp_y[j] * tg_zzzz_xxyyyyz_1[j] + 2.0 * fl1_fxn * tg_zzzz_xxyyyz_1[j];

                    tg_yzzzz_xxyyyzz_0[j] = pb_y * tg_zzzz_xxyyyzz_0[j] + wp_y[j] * tg_zzzz_xxyyyzz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxyyzz_1[j];

                    tg_yzzzz_xxyyzzz_0[j] = pb_y * tg_zzzz_xxyyzzz_0[j] + wp_y[j] * tg_zzzz_xxyyzzz_1[j] + fl1_fxn * tg_zzzz_xxyzzz_1[j];

                    tg_yzzzz_xxyzzzz_0[j] = pb_y * tg_zzzz_xxyzzzz_0[j] + wp_y[j] * tg_zzzz_xxyzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xxzzzz_1[j];

                    tg_yzzzz_xxzzzzz_0[j] = pb_y * tg_zzzz_xxzzzzz_0[j] + wp_y[j] * tg_zzzz_xxzzzzz_1[j];

                    tg_yzzzz_xyyyyyy_0[j] = pb_y * tg_zzzz_xyyyyyy_0[j] + wp_y[j] * tg_zzzz_xyyyyyy_1[j] + 3.0 * fl1_fxn * tg_zzzz_xyyyyy_1[j];

                    tg_yzzzz_xyyyyyz_0[j] = pb_y * tg_zzzz_xyyyyyz_0[j] + wp_y[j] * tg_zzzz_xyyyyyz_1[j] + 2.5 * fl1_fxn * tg_zzzz_xyyyyz_1[j];

                    tg_yzzzz_xyyyyzz_0[j] = pb_y * tg_zzzz_xyyyyzz_0[j] + wp_y[j] * tg_zzzz_xyyyyzz_1[j] + 2.0 * fl1_fxn * tg_zzzz_xyyyzz_1[j];

                    tg_yzzzz_xyyyzzz_0[j] = pb_y * tg_zzzz_xyyyzzz_0[j] + wp_y[j] * tg_zzzz_xyyyzzz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xyyzzz_1[j];

                    tg_yzzzz_xyyzzzz_0[j] = pb_y * tg_zzzz_xyyzzzz_0[j] + wp_y[j] * tg_zzzz_xyyzzzz_1[j] + fl1_fxn * tg_zzzz_xyzzzz_1[j];

                    tg_yzzzz_xyzzzzz_0[j] = pb_y * tg_zzzz_xyzzzzz_0[j] + wp_y[j] * tg_zzzz_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xzzzzz_1[j];

                    tg_yzzzz_xzzzzzz_0[j] = pb_y * tg_zzzz_xzzzzzz_0[j] + wp_y[j] * tg_zzzz_xzzzzzz_1[j];

                    tg_yzzzz_yyyyyyy_0[j] = pb_y * tg_zzzz_yyyyyyy_0[j] + wp_y[j] * tg_zzzz_yyyyyyy_1[j] + 3.5 * fl1_fxn * tg_zzzz_yyyyyy_1[j];

                    tg_yzzzz_yyyyyyz_0[j] = pb_y * tg_zzzz_yyyyyyz_0[j] + wp_y[j] * tg_zzzz_yyyyyyz_1[j] + 3.0 * fl1_fxn * tg_zzzz_yyyyyz_1[j];

                    tg_yzzzz_yyyyyzz_0[j] = pb_y * tg_zzzz_yyyyyzz_0[j] + wp_y[j] * tg_zzzz_yyyyyzz_1[j] + 2.5 * fl1_fxn * tg_zzzz_yyyyzz_1[j];

                    tg_yzzzz_yyyyzzz_0[j] = pb_y * tg_zzzz_yyyyzzz_0[j] + wp_y[j] * tg_zzzz_yyyyzzz_1[j] + 2.0 * fl1_fxn * tg_zzzz_yyyzzz_1[j];

                    tg_yzzzz_yyyzzzz_0[j] = pb_y * tg_zzzz_yyyzzzz_0[j] + wp_y[j] * tg_zzzz_yyyzzzz_1[j] + 1.5 * fl1_fxn * tg_zzzz_yyzzzz_1[j];

                    tg_yzzzz_yyzzzzz_0[j] = pb_y * tg_zzzz_yyzzzzz_0[j] + wp_y[j] * tg_zzzz_yyzzzzz_1[j] + fl1_fxn * tg_zzzz_yzzzzz_1[j];

                    tg_yzzzz_yzzzzzz_0[j] = pb_y * tg_zzzz_yzzzzzz_0[j] + wp_y[j] * tg_zzzz_yzzzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_zzzzzz_1[j];

                    tg_yzzzz_zzzzzzz_0[j] = pb_y * tg_zzzz_zzzzzzz_0[j] + wp_y[j] * tg_zzzz_zzzzzzz_1[j];

                    tg_zzzzz_xxxxxxx_0[j] = pb_z * tg_zzzz_xxxxxxx_0[j] + wp_z[j] * tg_zzzz_xxxxxxx_1[j] + 2.0 * fl1_fx * tg_zzz_xxxxxxx_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxxxxxx_1[j];

                    tg_zzzzz_xxxxxxy_0[j] = pb_z * tg_zzzz_xxxxxxy_0[j] + wp_z[j] * tg_zzzz_xxxxxxy_1[j] + 2.0 * fl1_fx * tg_zzz_xxxxxxy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxxxxxy_1[j];

                    tg_zzzzz_xxxxxxz_0[j] = pb_z * tg_zzzz_xxxxxxz_0[j] + wp_z[j] * tg_zzzz_xxxxxxz_1[j] + 2.0 * fl1_fx * tg_zzz_xxxxxxz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxxxxxz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xxxxxx_1[j];

                    tg_zzzzz_xxxxxyy_0[j] = pb_z * tg_zzzz_xxxxxyy_0[j] + wp_z[j] * tg_zzzz_xxxxxyy_1[j] + 2.0 * fl1_fx * tg_zzz_xxxxxyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxxxxyy_1[j];

                    tg_zzzzz_xxxxxyz_0[j] = pb_z * tg_zzzz_xxxxxyz_0[j] + wp_z[j] * tg_zzzz_xxxxxyz_1[j] + 2.0 * fl1_fx * tg_zzz_xxxxxyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxxxxyz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xxxxxy_1[j];

                    tg_zzzzz_xxxxxzz_0[j] = pb_z * tg_zzzz_xxxxxzz_0[j] + wp_z[j] * tg_zzzz_xxxxxzz_1[j] + 2.0 * fl1_fx * tg_zzz_xxxxxzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxxxxzz_1[j] + fl1_fxn * tg_zzzz_xxxxxz_1[j];

                    tg_zzzzz_xxxxyyy_0[j] = pb_z * tg_zzzz_xxxxyyy_0[j] + wp_z[j] * tg_zzzz_xxxxyyy_1[j] + 2.0 * fl1_fx * tg_zzz_xxxxyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxxxyyy_1[j];

                    tg_zzzzz_xxxxyyz_0[j] = pb_z * tg_zzzz_xxxxyyz_0[j] + wp_z[j] * tg_zzzz_xxxxyyz_1[j] + 2.0 * fl1_fx * tg_zzz_xxxxyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxxxyyz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xxxxyy_1[j];

                    tg_zzzzz_xxxxyzz_0[j] = pb_z * tg_zzzz_xxxxyzz_0[j] + wp_z[j] * tg_zzzz_xxxxyzz_1[j] + 2.0 * fl1_fx * tg_zzz_xxxxyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxxxyzz_1[j] + fl1_fxn * tg_zzzz_xxxxyz_1[j];

                    tg_zzzzz_xxxxzzz_0[j] = pb_z * tg_zzzz_xxxxzzz_0[j] + wp_z[j] * tg_zzzz_xxxxzzz_1[j] + 2.0 * fl1_fx * tg_zzz_xxxxzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxxxzzz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxxxzz_1[j];

                    tg_zzzzz_xxxyyyy_0[j] = pb_z * tg_zzzz_xxxyyyy_0[j] + wp_z[j] * tg_zzzz_xxxyyyy_1[j] + 2.0 * fl1_fx * tg_zzz_xxxyyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxxyyyy_1[j];

                    tg_zzzzz_xxxyyyz_0[j] = pb_z * tg_zzzz_xxxyyyz_0[j] + wp_z[j] * tg_zzzz_xxxyyyz_1[j] + 2.0 * fl1_fx * tg_zzz_xxxyyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxxyyyz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xxxyyy_1[j];

                    tg_zzzzz_xxxyyzz_0[j] = pb_z * tg_zzzz_xxxyyzz_0[j] + wp_z[j] * tg_zzzz_xxxyyzz_1[j] + 2.0 * fl1_fx * tg_zzz_xxxyyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxxyyzz_1[j] + fl1_fxn * tg_zzzz_xxxyyz_1[j];

                    tg_zzzzz_xxxyzzz_0[j] = pb_z * tg_zzzz_xxxyzzz_0[j] + wp_z[j] * tg_zzzz_xxxyzzz_1[j] + 2.0 * fl1_fx * tg_zzz_xxxyzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxxyzzz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxxyzz_1[j];

                    tg_zzzzz_xxxzzzz_0[j] = pb_z * tg_zzzz_xxxzzzz_0[j] + wp_z[j] * tg_zzzz_xxxzzzz_1[j] + 2.0 * fl1_fx * tg_zzz_xxxzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxxzzzz_1[j] + 2.0 * fl1_fxn * tg_zzzz_xxxzzz_1[j];

                    tg_zzzzz_xxyyyyy_0[j] = pb_z * tg_zzzz_xxyyyyy_0[j] + wp_z[j] * tg_zzzz_xxyyyyy_1[j] + 2.0 * fl1_fx * tg_zzz_xxyyyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxyyyyy_1[j];

                    tg_zzzzz_xxyyyyz_0[j] = pb_z * tg_zzzz_xxyyyyz_0[j] + wp_z[j] * tg_zzzz_xxyyyyz_1[j] + 2.0 * fl1_fx * tg_zzz_xxyyyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxyyyyz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xxyyyy_1[j];

                    tg_zzzzz_xxyyyzz_0[j] = pb_z * tg_zzzz_xxyyyzz_0[j] + wp_z[j] * tg_zzzz_xxyyyzz_1[j] + 2.0 * fl1_fx * tg_zzz_xxyyyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxyyyzz_1[j] + fl1_fxn * tg_zzzz_xxyyyz_1[j];

                    tg_zzzzz_xxyyzzz_0[j] = pb_z * tg_zzzz_xxyyzzz_0[j] + wp_z[j] * tg_zzzz_xxyyzzz_1[j] + 2.0 * fl1_fx * tg_zzz_xxyyzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxyyzzz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxyyzz_1[j];

                    tg_zzzzz_xxyzzzz_0[j] = pb_z * tg_zzzz_xxyzzzz_0[j] + wp_z[j] * tg_zzzz_xxyzzzz_1[j] + 2.0 * fl1_fx * tg_zzz_xxyzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxyzzzz_1[j] + 2.0 * fl1_fxn * tg_zzzz_xxyzzz_1[j];

                    tg_zzzzz_xxzzzzz_0[j] = pb_z * tg_zzzz_xxzzzzz_0[j] + wp_z[j] * tg_zzzz_xxzzzzz_1[j] + 2.0 * fl1_fx * tg_zzz_xxzzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxzzzzz_1[j] + 2.5 * fl1_fxn * tg_zzzz_xxzzzz_1[j];

                    tg_zzzzz_xyyyyyy_0[j] = pb_z * tg_zzzz_xyyyyyy_0[j] + wp_z[j] * tg_zzzz_xyyyyyy_1[j] + 2.0 * fl1_fx * tg_zzz_xyyyyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xyyyyyy_1[j];

                    tg_zzzzz_xyyyyyz_0[j] = pb_z * tg_zzzz_xyyyyyz_0[j] + wp_z[j] * tg_zzzz_xyyyyyz_1[j] + 2.0 * fl1_fx * tg_zzz_xyyyyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xyyyyyz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xyyyyy_1[j];

                    tg_zzzzz_xyyyyzz_0[j] = pb_z * tg_zzzz_xyyyyzz_0[j] + wp_z[j] * tg_zzzz_xyyyyzz_1[j] + 2.0 * fl1_fx * tg_zzz_xyyyyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xyyyyzz_1[j] + fl1_fxn * tg_zzzz_xyyyyz_1[j];

                    tg_zzzzz_xyyyzzz_0[j] = pb_z * tg_zzzz_xyyyzzz_0[j] + wp_z[j] * tg_zzzz_xyyyzzz_1[j] + 2.0 * fl1_fx * tg_zzz_xyyyzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xyyyzzz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xyyyzz_1[j];

                    tg_zzzzz_xyyzzzz_0[j] = pb_z * tg_zzzz_xyyzzzz_0[j] + wp_z[j] * tg_zzzz_xyyzzzz_1[j] + 2.0 * fl1_fx * tg_zzz_xyyzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xyyzzzz_1[j] + 2.0 * fl1_fxn * tg_zzzz_xyyzzz_1[j];

                    tg_zzzzz_xyzzzzz_0[j] = pb_z * tg_zzzz_xyzzzzz_0[j] + wp_z[j] * tg_zzzz_xyzzzzz_1[j] + 2.0 * fl1_fx * tg_zzz_xyzzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xyzzzzz_1[j] + 2.5 * fl1_fxn * tg_zzzz_xyzzzz_1[j];

                    tg_zzzzz_xzzzzzz_0[j] = pb_z * tg_zzzz_xzzzzzz_0[j] + wp_z[j] * tg_zzzz_xzzzzzz_1[j] + 2.0 * fl1_fx * tg_zzz_xzzzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xzzzzzz_1[j] + 3.0 * fl1_fxn * tg_zzzz_xzzzzz_1[j];

                    tg_zzzzz_yyyyyyy_0[j] = pb_z * tg_zzzz_yyyyyyy_0[j] + wp_z[j] * tg_zzzz_yyyyyyy_1[j] + 2.0 * fl1_fx * tg_zzz_yyyyyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_yyyyyyy_1[j];

                    tg_zzzzz_yyyyyyz_0[j] = pb_z * tg_zzzz_yyyyyyz_0[j] + wp_z[j] * tg_zzzz_yyyyyyz_1[j] + 2.0 * fl1_fx * tg_zzz_yyyyyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_yyyyyyz_1[j] + 0.5 * fl1_fxn * tg_zzzz_yyyyyy_1[j];

                    tg_zzzzz_yyyyyzz_0[j] = pb_z * tg_zzzz_yyyyyzz_0[j] + wp_z[j] * tg_zzzz_yyyyyzz_1[j] + 2.0 * fl1_fx * tg_zzz_yyyyyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_yyyyyzz_1[j] + fl1_fxn * tg_zzzz_yyyyyz_1[j];

                    tg_zzzzz_yyyyzzz_0[j] = pb_z * tg_zzzz_yyyyzzz_0[j] + wp_z[j] * tg_zzzz_yyyyzzz_1[j] + 2.0 * fl1_fx * tg_zzz_yyyyzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_yyyyzzz_1[j] + 1.5 * fl1_fxn * tg_zzzz_yyyyzz_1[j];

                    tg_zzzzz_yyyzzzz_0[j] = pb_z * tg_zzzz_yyyzzzz_0[j] + wp_z[j] * tg_zzzz_yyyzzzz_1[j] + 2.0 * fl1_fx * tg_zzz_yyyzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_yyyzzzz_1[j] + 2.0 * fl1_fxn * tg_zzzz_yyyzzz_1[j];

                    tg_zzzzz_yyzzzzz_0[j] = pb_z * tg_zzzz_yyzzzzz_0[j] + wp_z[j] * tg_zzzz_yyzzzzz_1[j] + 2.0 * fl1_fx * tg_zzz_yyzzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_yyzzzzz_1[j] + 2.5 * fl1_fxn * tg_zzzz_yyzzzz_1[j];

                    tg_zzzzz_yzzzzzz_0[j] = pb_z * tg_zzzz_yzzzzzz_0[j] + wp_z[j] * tg_zzzz_yzzzzzz_1[j] + 2.0 * fl1_fx * tg_zzz_yzzzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_yzzzzzz_1[j] + 3.0 * fl1_fxn * tg_zzzz_yzzzzz_1[j];

                    tg_zzzzz_zzzzzzz_0[j] = pb_z * tg_zzzz_zzzzzzz_0[j] + wp_z[j] * tg_zzzz_zzzzzzz_1[j] + 2.0 * fl1_fx * tg_zzz_zzzzzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_zzzzzzz_1[j] + 3.5 * fl1_fxn * tg_zzzz_zzzzzz_1[j];
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

