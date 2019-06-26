//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectronRepulsionRecFuncForGK.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSGSK(      CMemBlock2D<double>* primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSGSK_0_90(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSGSK_90_180(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 

        erirecfunc::compElectronRepulsionForSGSK_180_270(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSGSK_270_360(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSGSK_360_450(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSGSK_450_540(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 
    }

    void
    compElectronRepulsionForSGSK_0_90(      CMemBlock2D<double>* primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& wpDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetPrimPairs,
                                      const int32_t              iContrPair)
    {
        // Batch of Integrals (0,90)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {4, -1, -1, -1},
                                             {7, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_7_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {7, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_7_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_7_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_7_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_2_7_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {7, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_2_7_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {7, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {6, -1, -1, -1}, 
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

                auto tg_xxx_xxxxxxx_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx); 

                auto tg_xxx_xxxxxxy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 1); 

                auto tg_xxx_xxxxxxz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 2); 

                auto tg_xxx_xxxxxyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 3); 

                auto tg_xxx_xxxxxyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 4); 

                auto tg_xxx_xxxxxzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 5); 

                auto tg_xxx_xxxxyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 6); 

                auto tg_xxx_xxxxyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 7); 

                auto tg_xxx_xxxxyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 8); 

                auto tg_xxx_xxxxzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 9); 

                auto tg_xxx_xxxyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 10); 

                auto tg_xxx_xxxyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 11); 

                auto tg_xxx_xxxyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 12); 

                auto tg_xxx_xxxyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 13); 

                auto tg_xxx_xxxzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 14); 

                auto tg_xxx_xxyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 15); 

                auto tg_xxx_xxyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 16); 

                auto tg_xxx_xxyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 17); 

                auto tg_xxx_xxyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 18); 

                auto tg_xxx_xxyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 19); 

                auto tg_xxx_xxzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 20); 

                auto tg_xxx_xyyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 21); 

                auto tg_xxx_xyyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 22); 

                auto tg_xxx_xyyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 23); 

                auto tg_xxx_xyyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 24); 

                auto tg_xxx_xyyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 25); 

                auto tg_xxx_xyzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 26); 

                auto tg_xxx_xzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 27); 

                auto tg_xxx_yyyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 28); 

                auto tg_xxx_yyyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 29); 

                auto tg_xxx_yyyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 30); 

                auto tg_xxx_yyyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 31); 

                auto tg_xxx_yyyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 32); 

                auto tg_xxx_yyzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 33); 

                auto tg_xxx_yzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 34); 

                auto tg_xxx_zzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 35); 

                auto tg_xxy_xxxxxxx_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 36); 

                auto tg_xxy_xxxxxxy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 37); 

                auto tg_xxy_xxxxxxz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 38); 

                auto tg_xxy_xxxxxyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 39); 

                auto tg_xxy_xxxxxyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 40); 

                auto tg_xxy_xxxxxzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 41); 

                auto tg_xxy_xxxxyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 42); 

                auto tg_xxy_xxxxyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 43); 

                auto tg_xxy_xxxxyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 44); 

                auto tg_xxy_xxxxzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 45); 

                auto tg_xxy_xxxyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 46); 

                auto tg_xxy_xxxyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 47); 

                auto tg_xxy_xxxyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 48); 

                auto tg_xxy_xxxyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 49); 

                auto tg_xxy_xxxzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 50); 

                auto tg_xxy_xxyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 51); 

                auto tg_xxy_xxyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 52); 

                auto tg_xxy_xxyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 53); 

                auto tg_xxy_xxyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 54); 

                auto tg_xxy_xxyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 55); 

                auto tg_xxy_xxzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 56); 

                auto tg_xxy_xyyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 57); 

                auto tg_xxy_xyyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 58); 

                auto tg_xxy_xyyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 59); 

                auto tg_xxy_xyyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 60); 

                auto tg_xxy_xyyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 61); 

                auto tg_xxy_xyzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 62); 

                auto tg_xxy_xzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 63); 

                auto tg_xxy_yyyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 64); 

                auto tg_xxy_yyyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 65); 

                auto tg_xxy_yyyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 66); 

                auto tg_xxy_yyyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 67); 

                auto tg_xxy_yyyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 68); 

                auto tg_xxy_yyzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 69); 

                auto tg_xxy_yzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 70); 

                auto tg_xxy_zzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 71); 

                auto tg_xxz_xxxxxxx_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 72); 

                auto tg_xxz_xxxxxxy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 73); 

                auto tg_xxz_xxxxxxz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 74); 

                auto tg_xxz_xxxxxyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 75); 

                auto tg_xxz_xxxxxyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 76); 

                auto tg_xxz_xxxxxzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 77); 

                auto tg_xxz_xxxxyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 78); 

                auto tg_xxz_xxxxyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 79); 

                auto tg_xxz_xxxxyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 80); 

                auto tg_xxz_xxxxzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 81); 

                auto tg_xxz_xxxyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 82); 

                auto tg_xxz_xxxyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 83); 

                auto tg_xxz_xxxyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 84); 

                auto tg_xxz_xxxyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 85); 

                auto tg_xxz_xxxzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 86); 

                auto tg_xxz_xxyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 87); 

                auto tg_xxz_xxyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 88); 

                auto tg_xxz_xxyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 89); 

                auto tg_xxx_xxxxxxx_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx); 

                auto tg_xxx_xxxxxxy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 1); 

                auto tg_xxx_xxxxxxz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 2); 

                auto tg_xxx_xxxxxyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 3); 

                auto tg_xxx_xxxxxyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 4); 

                auto tg_xxx_xxxxxzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 5); 

                auto tg_xxx_xxxxyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 6); 

                auto tg_xxx_xxxxyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 7); 

                auto tg_xxx_xxxxyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 8); 

                auto tg_xxx_xxxxzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 9); 

                auto tg_xxx_xxxyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 10); 

                auto tg_xxx_xxxyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 11); 

                auto tg_xxx_xxxyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 12); 

                auto tg_xxx_xxxyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 13); 

                auto tg_xxx_xxxzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 14); 

                auto tg_xxx_xxyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 15); 

                auto tg_xxx_xxyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 16); 

                auto tg_xxx_xxyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 17); 

                auto tg_xxx_xxyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 18); 

                auto tg_xxx_xxyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 19); 

                auto tg_xxx_xxzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 20); 

                auto tg_xxx_xyyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 21); 

                auto tg_xxx_xyyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 22); 

                auto tg_xxx_xyyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 23); 

                auto tg_xxx_xyyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 24); 

                auto tg_xxx_xyyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 25); 

                auto tg_xxx_xyzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 26); 

                auto tg_xxx_xzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 27); 

                auto tg_xxx_yyyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 28); 

                auto tg_xxx_yyyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 29); 

                auto tg_xxx_yyyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 30); 

                auto tg_xxx_yyyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 31); 

                auto tg_xxx_yyyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 32); 

                auto tg_xxx_yyzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 33); 

                auto tg_xxx_yzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 34); 

                auto tg_xxx_zzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 35); 

                auto tg_xxy_xxxxxxx_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 36); 

                auto tg_xxy_xxxxxxy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 37); 

                auto tg_xxy_xxxxxxz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 38); 

                auto tg_xxy_xxxxxyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 39); 

                auto tg_xxy_xxxxxyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 40); 

                auto tg_xxy_xxxxxzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 41); 

                auto tg_xxy_xxxxyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 42); 

                auto tg_xxy_xxxxyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 43); 

                auto tg_xxy_xxxxyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 44); 

                auto tg_xxy_xxxxzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 45); 

                auto tg_xxy_xxxyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 46); 

                auto tg_xxy_xxxyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 47); 

                auto tg_xxy_xxxyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 48); 

                auto tg_xxy_xxxyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 49); 

                auto tg_xxy_xxxzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 50); 

                auto tg_xxy_xxyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 51); 

                auto tg_xxy_xxyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 52); 

                auto tg_xxy_xxyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 53); 

                auto tg_xxy_xxyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 54); 

                auto tg_xxy_xxyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 55); 

                auto tg_xxy_xxzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 56); 

                auto tg_xxy_xyyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 57); 

                auto tg_xxy_xyyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 58); 

                auto tg_xxy_xyyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 59); 

                auto tg_xxy_xyyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 60); 

                auto tg_xxy_xyyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 61); 

                auto tg_xxy_xyzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 62); 

                auto tg_xxy_xzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 63); 

                auto tg_xxy_yyyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 64); 

                auto tg_xxy_yyyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 65); 

                auto tg_xxy_yyyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 66); 

                auto tg_xxy_yyyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 67); 

                auto tg_xxy_yyyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 68); 

                auto tg_xxy_yyzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 69); 

                auto tg_xxy_yzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 70); 

                auto tg_xxy_zzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 71); 

                auto tg_xxz_xxxxxxx_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 72); 

                auto tg_xxz_xxxxxxy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 73); 

                auto tg_xxz_xxxxxxz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 74); 

                auto tg_xxz_xxxxxyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 75); 

                auto tg_xxz_xxxxxyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 76); 

                auto tg_xxz_xxxxxzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 77); 

                auto tg_xxz_xxxxyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 78); 

                auto tg_xxz_xxxxyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 79); 

                auto tg_xxz_xxxxyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 80); 

                auto tg_xxz_xxxxzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 81); 

                auto tg_xxz_xxxyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 82); 

                auto tg_xxz_xxxyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 83); 

                auto tg_xxz_xxxyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 84); 

                auto tg_xxz_xxxyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 85); 

                auto tg_xxz_xxxzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 86); 

                auto tg_xxz_xxyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 87); 

                auto tg_xxz_xxyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 88); 

                auto tg_xxz_xxyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 89); 

                auto tg_xx_xxxxxxx_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx); 

                auto tg_xx_xxxxxxy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 1); 

                auto tg_xx_xxxxxxz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 2); 

                auto tg_xx_xxxxxyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 3); 

                auto tg_xx_xxxxxyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 4); 

                auto tg_xx_xxxxxzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 5); 

                auto tg_xx_xxxxyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 6); 

                auto tg_xx_xxxxyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 7); 

                auto tg_xx_xxxxyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 8); 

                auto tg_xx_xxxxzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 9); 

                auto tg_xx_xxxyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 10); 

                auto tg_xx_xxxyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 11); 

                auto tg_xx_xxxyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 12); 

                auto tg_xx_xxxyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 13); 

                auto tg_xx_xxxzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 14); 

                auto tg_xx_xxyyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 15); 

                auto tg_xx_xxyyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 16); 

                auto tg_xx_xxyyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 17); 

                auto tg_xx_xxyyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 18); 

                auto tg_xx_xxyzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 19); 

                auto tg_xx_xxzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 20); 

                auto tg_xx_xyyyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 21); 

                auto tg_xx_xyyyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 22); 

                auto tg_xx_xyyyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 23); 

                auto tg_xx_xyyyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 24); 

                auto tg_xx_xyyzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 25); 

                auto tg_xx_xyzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 26); 

                auto tg_xx_xzzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 27); 

                auto tg_xx_yyyyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 28); 

                auto tg_xx_yyyyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 29); 

                auto tg_xx_yyyyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 30); 

                auto tg_xx_yyyyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 31); 

                auto tg_xx_yyyzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 32); 

                auto tg_xx_yyzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 33); 

                auto tg_xx_yzzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 34); 

                auto tg_xx_zzzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 35); 

                auto tg_xy_xxxxxxx_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 36); 

                auto tg_xy_xxxxxxy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 37); 

                auto tg_xy_xxxxxxz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 38); 

                auto tg_xy_xxxxxyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 39); 

                auto tg_xy_xxxxxyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 40); 

                auto tg_xy_xxxxxzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 41); 

                auto tg_xy_xxxxyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 42); 

                auto tg_xy_xxxxyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 43); 

                auto tg_xy_xxxxyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 44); 

                auto tg_xy_xxxxzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 45); 

                auto tg_xy_xxxyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 46); 

                auto tg_xy_xxxyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 47); 

                auto tg_xy_xxxyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 48); 

                auto tg_xy_xxxyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 49); 

                auto tg_xy_xxxzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 50); 

                auto tg_xy_xxyyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 51); 

                auto tg_xy_xxyyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 52); 

                auto tg_xy_xxyyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 53); 

                auto tg_xy_xxyyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 54); 

                auto tg_xy_xxyzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 55); 

                auto tg_xy_xxzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 56); 

                auto tg_xy_xyyyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 57); 

                auto tg_xy_xyyyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 58); 

                auto tg_xy_xyyyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 59); 

                auto tg_xy_xyyyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 60); 

                auto tg_xy_xyyzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 61); 

                auto tg_xy_xyzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 62); 

                auto tg_xy_xzzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 63); 

                auto tg_xy_yyyyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 64); 

                auto tg_xy_yyyyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 65); 

                auto tg_xy_yyyyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 66); 

                auto tg_xy_yyyyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 67); 

                auto tg_xy_yyyzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 68); 

                auto tg_xy_yyzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 69); 

                auto tg_xy_yzzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 70); 

                auto tg_xy_zzzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 71); 

                auto tg_xz_xxxxxxx_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 72); 

                auto tg_xz_xxxxxxy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 73); 

                auto tg_xz_xxxxxxz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 74); 

                auto tg_xz_xxxxxyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 75); 

                auto tg_xz_xxxxxyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 76); 

                auto tg_xz_xxxxxzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 77); 

                auto tg_xz_xxxxyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 78); 

                auto tg_xz_xxxxyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 79); 

                auto tg_xz_xxxxyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 80); 

                auto tg_xz_xxxxzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 81); 

                auto tg_xz_xxxyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 82); 

                auto tg_xz_xxxyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 83); 

                auto tg_xz_xxxyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 84); 

                auto tg_xz_xxxyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 85); 

                auto tg_xz_xxxzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 86); 

                auto tg_xz_xxyyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 87); 

                auto tg_xz_xxyyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 88); 

                auto tg_xz_xxyyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 89); 

                auto tg_xx_xxxxxxx_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx); 

                auto tg_xx_xxxxxxy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 1); 

                auto tg_xx_xxxxxxz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 2); 

                auto tg_xx_xxxxxyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 3); 

                auto tg_xx_xxxxxyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 4); 

                auto tg_xx_xxxxxzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 5); 

                auto tg_xx_xxxxyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 6); 

                auto tg_xx_xxxxyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 7); 

                auto tg_xx_xxxxyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 8); 

                auto tg_xx_xxxxzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 9); 

                auto tg_xx_xxxyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 10); 

                auto tg_xx_xxxyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 11); 

                auto tg_xx_xxxyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 12); 

                auto tg_xx_xxxyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 13); 

                auto tg_xx_xxxzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 14); 

                auto tg_xx_xxyyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 15); 

                auto tg_xx_xxyyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 16); 

                auto tg_xx_xxyyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 17); 

                auto tg_xx_xxyyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 18); 

                auto tg_xx_xxyzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 19); 

                auto tg_xx_xxzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 20); 

                auto tg_xx_xyyyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 21); 

                auto tg_xx_xyyyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 22); 

                auto tg_xx_xyyyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 23); 

                auto tg_xx_xyyyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 24); 

                auto tg_xx_xyyzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 25); 

                auto tg_xx_xyzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 26); 

                auto tg_xx_xzzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 27); 

                auto tg_xx_yyyyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 28); 

                auto tg_xx_yyyyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 29); 

                auto tg_xx_yyyyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 30); 

                auto tg_xx_yyyyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 31); 

                auto tg_xx_yyyzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 32); 

                auto tg_xx_yyzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 33); 

                auto tg_xx_yzzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 34); 

                auto tg_xx_zzzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 35); 

                auto tg_xy_xxxxxxx_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 36); 

                auto tg_xy_xxxxxxy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 37); 

                auto tg_xy_xxxxxxz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 38); 

                auto tg_xy_xxxxxyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 39); 

                auto tg_xy_xxxxxyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 40); 

                auto tg_xy_xxxxxzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 41); 

                auto tg_xy_xxxxyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 42); 

                auto tg_xy_xxxxyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 43); 

                auto tg_xy_xxxxyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 44); 

                auto tg_xy_xxxxzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 45); 

                auto tg_xy_xxxyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 46); 

                auto tg_xy_xxxyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 47); 

                auto tg_xy_xxxyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 48); 

                auto tg_xy_xxxyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 49); 

                auto tg_xy_xxxzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 50); 

                auto tg_xy_xxyyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 51); 

                auto tg_xy_xxyyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 52); 

                auto tg_xy_xxyyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 53); 

                auto tg_xy_xxyyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 54); 

                auto tg_xy_xxyzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 55); 

                auto tg_xy_xxzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 56); 

                auto tg_xy_xyyyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 57); 

                auto tg_xy_xyyyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 58); 

                auto tg_xy_xyyyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 59); 

                auto tg_xy_xyyyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 60); 

                auto tg_xy_xyyzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 61); 

                auto tg_xy_xyzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 62); 

                auto tg_xy_xzzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 63); 

                auto tg_xy_yyyyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 64); 

                auto tg_xy_yyyyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 65); 

                auto tg_xy_yyyyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 66); 

                auto tg_xy_yyyyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 67); 

                auto tg_xy_yyyzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 68); 

                auto tg_xy_yyzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 69); 

                auto tg_xy_yzzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 70); 

                auto tg_xy_zzzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 71); 

                auto tg_xz_xxxxxxx_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 72); 

                auto tg_xz_xxxxxxy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 73); 

                auto tg_xz_xxxxxxz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 74); 

                auto tg_xz_xxxxxyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 75); 

                auto tg_xz_xxxxxyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 76); 

                auto tg_xz_xxxxxzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 77); 

                auto tg_xz_xxxxyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 78); 

                auto tg_xz_xxxxyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 79); 

                auto tg_xz_xxxxyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 80); 

                auto tg_xz_xxxxzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 81); 

                auto tg_xz_xxxyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 82); 

                auto tg_xz_xxxyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 83); 

                auto tg_xz_xxxyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 84); 

                auto tg_xz_xxxyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 85); 

                auto tg_xz_xxxzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 86); 

                auto tg_xz_xxyyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 87); 

                auto tg_xz_xxyyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 88); 

                auto tg_xz_xxyyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 89); 

                auto tg_xxx_xxxxxx_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx); 

                auto tg_xxx_xxxxxy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 1); 

                auto tg_xxx_xxxxxz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 2); 

                auto tg_xxx_xxxxyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 3); 

                auto tg_xxx_xxxxyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 4); 

                auto tg_xxx_xxxxzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 5); 

                auto tg_xxx_xxxyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 6); 

                auto tg_xxx_xxxyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 7); 

                auto tg_xxx_xxxyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 8); 

                auto tg_xxx_xxxzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 9); 

                auto tg_xxx_xxyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 10); 

                auto tg_xxx_xxyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 11); 

                auto tg_xxx_xxyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 12); 

                auto tg_xxx_xxyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 13); 

                auto tg_xxx_xxzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 14); 

                auto tg_xxx_xyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 15); 

                auto tg_xxx_xyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 16); 

                auto tg_xxx_xyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 17); 

                auto tg_xxx_xyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 18); 

                auto tg_xxx_xyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 19); 

                auto tg_xxx_xzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 20); 

                auto tg_xxx_yyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 21); 

                auto tg_xxx_yyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 22); 

                auto tg_xxx_yyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 23); 

                auto tg_xxx_yyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 24); 

                auto tg_xxx_yyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 25); 

                auto tg_xxx_yzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 26); 

                auto tg_xxx_zzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 27); 

                auto tg_xxy_xxxxxx_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 28); 

                auto tg_xxy_xxxxxy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 29); 

                auto tg_xxy_xxxxxz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 30); 

                auto tg_xxy_xxxxyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 31); 

                auto tg_xxy_xxxxyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 32); 

                auto tg_xxy_xxxxzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 33); 

                auto tg_xxy_xxxyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 34); 

                auto tg_xxy_xxxyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 35); 

                auto tg_xxy_xxxyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 36); 

                auto tg_xxy_xxxzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 37); 

                auto tg_xxy_xxyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 38); 

                auto tg_xxy_xxyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 39); 

                auto tg_xxy_xxyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 40); 

                auto tg_xxy_xxyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 41); 

                auto tg_xxy_xxzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 42); 

                auto tg_xxy_xyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 43); 

                auto tg_xxy_xyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 44); 

                auto tg_xxy_xyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 45); 

                auto tg_xxy_xyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 46); 

                auto tg_xxy_xyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 47); 

                auto tg_xxy_xzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 48); 

                auto tg_xxy_yyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 49); 

                auto tg_xxy_yyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 50); 

                auto tg_xxy_yyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 51); 

                auto tg_xxy_yyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 52); 

                auto tg_xxy_yyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 53); 

                auto tg_xxy_yzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 54); 

                auto tg_xxy_zzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 55); 

                auto tg_xxz_xxxxxx_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 56); 

                auto tg_xxz_xxxxxy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 57); 

                auto tg_xxz_xxxxxz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 58); 

                auto tg_xxz_xxxxyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 59); 

                auto tg_xxz_xxxxyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 60); 

                auto tg_xxz_xxxxzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 61); 

                auto tg_xxz_xxxyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 62); 

                auto tg_xxz_xxxyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 63); 

                auto tg_xxz_xxxyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 64); 

                auto tg_xxz_xxxzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 65); 

                auto tg_xxz_xxyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 66); 

                auto tg_xxz_xxyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 67); 

                auto tg_xxz_xxyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 68); 

                auto tg_xxz_xxyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 69); 

                auto tg_xxz_xxzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 70); 

                auto tg_xxz_xyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 71); 

                auto tg_xxz_xyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 72); 

                auto tg_xxz_xyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 73); 

                // set up pointers to integrals

                auto tg_xxxx_xxxxxxx_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx); 

                auto tg_xxxx_xxxxxxy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 1); 

                auto tg_xxxx_xxxxxxz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 2); 

                auto tg_xxxx_xxxxxyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 3); 

                auto tg_xxxx_xxxxxyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 4); 

                auto tg_xxxx_xxxxxzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 5); 

                auto tg_xxxx_xxxxyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 6); 

                auto tg_xxxx_xxxxyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 7); 

                auto tg_xxxx_xxxxyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 8); 

                auto tg_xxxx_xxxxzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 9); 

                auto tg_xxxx_xxxyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 10); 

                auto tg_xxxx_xxxyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 11); 

                auto tg_xxxx_xxxyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 12); 

                auto tg_xxxx_xxxyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 13); 

                auto tg_xxxx_xxxzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 14); 

                auto tg_xxxx_xxyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 15); 

                auto tg_xxxx_xxyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 16); 

                auto tg_xxxx_xxyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 17); 

                auto tg_xxxx_xxyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 18); 

                auto tg_xxxx_xxyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 19); 

                auto tg_xxxx_xxzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 20); 

                auto tg_xxxx_xyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 21); 

                auto tg_xxxx_xyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 22); 

                auto tg_xxxx_xyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 23); 

                auto tg_xxxx_xyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 24); 

                auto tg_xxxx_xyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 25); 

                auto tg_xxxx_xyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 26); 

                auto tg_xxxx_xzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 27); 

                auto tg_xxxx_yyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 28); 

                auto tg_xxxx_yyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 29); 

                auto tg_xxxx_yyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 30); 

                auto tg_xxxx_yyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 31); 

                auto tg_xxxx_yyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 32); 

                auto tg_xxxx_yyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 33); 

                auto tg_xxxx_yzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 34); 

                auto tg_xxxx_zzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 35); 

                auto tg_xxxy_xxxxxxx_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 36); 

                auto tg_xxxy_xxxxxxy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 37); 

                auto tg_xxxy_xxxxxxz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 38); 

                auto tg_xxxy_xxxxxyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 39); 

                auto tg_xxxy_xxxxxyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 40); 

                auto tg_xxxy_xxxxxzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 41); 

                auto tg_xxxy_xxxxyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 42); 

                auto tg_xxxy_xxxxyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 43); 

                auto tg_xxxy_xxxxyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 44); 

                auto tg_xxxy_xxxxzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 45); 

                auto tg_xxxy_xxxyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 46); 

                auto tg_xxxy_xxxyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 47); 

                auto tg_xxxy_xxxyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 48); 

                auto tg_xxxy_xxxyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 49); 

                auto tg_xxxy_xxxzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 50); 

                auto tg_xxxy_xxyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 51); 

                auto tg_xxxy_xxyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 52); 

                auto tg_xxxy_xxyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 53); 

                auto tg_xxxy_xxyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 54); 

                auto tg_xxxy_xxyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 55); 

                auto tg_xxxy_xxzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 56); 

                auto tg_xxxy_xyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 57); 

                auto tg_xxxy_xyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 58); 

                auto tg_xxxy_xyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 59); 

                auto tg_xxxy_xyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 60); 

                auto tg_xxxy_xyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 61); 

                auto tg_xxxy_xyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 62); 

                auto tg_xxxy_xzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 63); 

                auto tg_xxxy_yyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 64); 

                auto tg_xxxy_yyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 65); 

                auto tg_xxxy_yyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 66); 

                auto tg_xxxy_yyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 67); 

                auto tg_xxxy_yyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 68); 

                auto tg_xxxy_yyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 69); 

                auto tg_xxxy_yzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 70); 

                auto tg_xxxy_zzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 71); 

                auto tg_xxxz_xxxxxxx_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 72); 

                auto tg_xxxz_xxxxxxy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 73); 

                auto tg_xxxz_xxxxxxz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 74); 

                auto tg_xxxz_xxxxxyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 75); 

                auto tg_xxxz_xxxxxyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 76); 

                auto tg_xxxz_xxxxxzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 77); 

                auto tg_xxxz_xxxxyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 78); 

                auto tg_xxxz_xxxxyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 79); 

                auto tg_xxxz_xxxxyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 80); 

                auto tg_xxxz_xxxxzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 81); 

                auto tg_xxxz_xxxyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 82); 

                auto tg_xxxz_xxxyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 83); 

                auto tg_xxxz_xxxyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 84); 

                auto tg_xxxz_xxxyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 85); 

                auto tg_xxxz_xxxzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 86); 

                auto tg_xxxz_xxyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 87); 

                auto tg_xxxz_xxyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 88); 

                auto tg_xxxz_xxyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 89); 

                // Batch of Integrals (0,90)

                #pragma omp simd aligned(fxn, fza, tg_xx_xxxxxxx_0, tg_xx_xxxxxxx_1, tg_xx_xxxxxxy_0, \
                                         tg_xx_xxxxxxy_1, tg_xx_xxxxxxz_0, tg_xx_xxxxxxz_1, tg_xx_xxxxxyy_0, tg_xx_xxxxxyy_1, \
                                         tg_xx_xxxxxyz_0, tg_xx_xxxxxyz_1, tg_xx_xxxxxzz_0, tg_xx_xxxxxzz_1, tg_xx_xxxxyyy_0, \
                                         tg_xx_xxxxyyy_1, tg_xx_xxxxyyz_0, tg_xx_xxxxyyz_1, tg_xx_xxxxyzz_0, tg_xx_xxxxyzz_1, \
                                         tg_xx_xxxxzzz_0, tg_xx_xxxxzzz_1, tg_xx_xxxyyyy_0, tg_xx_xxxyyyy_1, tg_xx_xxxyyyz_0, \
                                         tg_xx_xxxyyyz_1, tg_xx_xxxyyzz_0, tg_xx_xxxyyzz_1, tg_xx_xxxyzzz_0, tg_xx_xxxyzzz_1, \
                                         tg_xx_xxxzzzz_0, tg_xx_xxxzzzz_1, tg_xx_xxyyyyy_0, tg_xx_xxyyyyy_1, tg_xx_xxyyyyz_0, \
                                         tg_xx_xxyyyyz_1, tg_xx_xxyyyzz_0, tg_xx_xxyyyzz_1, tg_xx_xxyyzzz_0, tg_xx_xxyyzzz_1, \
                                         tg_xx_xxyzzzz_0, tg_xx_xxyzzzz_1, tg_xx_xxzzzzz_0, tg_xx_xxzzzzz_1, tg_xx_xyyyyyy_0, \
                                         tg_xx_xyyyyyy_1, tg_xx_xyyyyyz_0, tg_xx_xyyyyyz_1, tg_xx_xyyyyzz_0, tg_xx_xyyyyzz_1, \
                                         tg_xx_xyyyzzz_0, tg_xx_xyyyzzz_1, tg_xx_xyyzzzz_0, tg_xx_xyyzzzz_1, tg_xx_xyzzzzz_0, \
                                         tg_xx_xyzzzzz_1, tg_xx_xzzzzzz_0, tg_xx_xzzzzzz_1, tg_xx_yyyyyyy_0, tg_xx_yyyyyyy_1, \
                                         tg_xx_yyyyyyz_0, tg_xx_yyyyyyz_1, tg_xx_yyyyyzz_0, tg_xx_yyyyyzz_1, tg_xx_yyyyzzz_0, \
                                         tg_xx_yyyyzzz_1, tg_xx_yyyzzzz_0, tg_xx_yyyzzzz_1, tg_xx_yyzzzzz_0, tg_xx_yyzzzzz_1, \
                                         tg_xx_yzzzzzz_0, tg_xx_yzzzzzz_1, tg_xx_zzzzzzz_0, tg_xx_zzzzzzz_1, tg_xxx_xxxxxx_1, \
                                         tg_xxx_xxxxxxx_0, tg_xxx_xxxxxxx_1, tg_xxx_xxxxxxy_0, tg_xxx_xxxxxxy_1, \
                                         tg_xxx_xxxxxxz_0, tg_xxx_xxxxxxz_1, tg_xxx_xxxxxy_1, tg_xxx_xxxxxyy_0, \
                                         tg_xxx_xxxxxyy_1, tg_xxx_xxxxxyz_0, tg_xxx_xxxxxyz_1, tg_xxx_xxxxxz_1, \
                                         tg_xxx_xxxxxzz_0, tg_xxx_xxxxxzz_1, tg_xxx_xxxxyy_1, tg_xxx_xxxxyyy_0, \
                                         tg_xxx_xxxxyyy_1, tg_xxx_xxxxyyz_0, tg_xxx_xxxxyyz_1, tg_xxx_xxxxyz_1, \
                                         tg_xxx_xxxxyzz_0, tg_xxx_xxxxyzz_1, tg_xxx_xxxxzz_1, tg_xxx_xxxxzzz_0, \
                                         tg_xxx_xxxxzzz_1, tg_xxx_xxxyyy_1, tg_xxx_xxxyyyy_0, tg_xxx_xxxyyyy_1, \
                                         tg_xxx_xxxyyyz_0, tg_xxx_xxxyyyz_1, tg_xxx_xxxyyz_1, tg_xxx_xxxyyzz_0, \
                                         tg_xxx_xxxyyzz_1, tg_xxx_xxxyzz_1, tg_xxx_xxxyzzz_0, tg_xxx_xxxyzzz_1, \
                                         tg_xxx_xxxzzz_1, tg_xxx_xxxzzzz_0, tg_xxx_xxxzzzz_1, tg_xxx_xxyyyy_1, \
                                         tg_xxx_xxyyyyy_0, tg_xxx_xxyyyyy_1, tg_xxx_xxyyyyz_0, tg_xxx_xxyyyyz_1, \
                                         tg_xxx_xxyyyz_1, tg_xxx_xxyyyzz_0, tg_xxx_xxyyyzz_1, tg_xxx_xxyyzz_1, \
                                         tg_xxx_xxyyzzz_0, tg_xxx_xxyyzzz_1, tg_xxx_xxyzzz_1, tg_xxx_xxyzzzz_0, \
                                         tg_xxx_xxyzzzz_1, tg_xxx_xxzzzz_1, tg_xxx_xxzzzzz_0, tg_xxx_xxzzzzz_1, \
                                         tg_xxx_xyyyyy_1, tg_xxx_xyyyyyy_0, tg_xxx_xyyyyyy_1, tg_xxx_xyyyyyz_0, \
                                         tg_xxx_xyyyyyz_1, tg_xxx_xyyyyz_1, tg_xxx_xyyyyzz_0, tg_xxx_xyyyyzz_1, \
                                         tg_xxx_xyyyzz_1, tg_xxx_xyyyzzz_0, tg_xxx_xyyyzzz_1, tg_xxx_xyyzzz_1, \
                                         tg_xxx_xyyzzzz_0, tg_xxx_xyyzzzz_1, tg_xxx_xyzzzz_1, tg_xxx_xyzzzzz_0, \
                                         tg_xxx_xyzzzzz_1, tg_xxx_xzzzzz_1, tg_xxx_xzzzzzz_0, tg_xxx_xzzzzzz_1, \
                                         tg_xxx_yyyyyy_1, tg_xxx_yyyyyyy_0, tg_xxx_yyyyyyy_1, tg_xxx_yyyyyyz_0, \
                                         tg_xxx_yyyyyyz_1, tg_xxx_yyyyyz_1, tg_xxx_yyyyyzz_0, tg_xxx_yyyyyzz_1, \
                                         tg_xxx_yyyyzz_1, tg_xxx_yyyyzzz_0, tg_xxx_yyyyzzz_1, tg_xxx_yyyzzz_1, \
                                         tg_xxx_yyyzzzz_0, tg_xxx_yyyzzzz_1, tg_xxx_yyzzzz_1, tg_xxx_yyzzzzz_0, \
                                         tg_xxx_yyzzzzz_1, tg_xxx_yzzzzz_1, tg_xxx_yzzzzzz_0, tg_xxx_yzzzzzz_1, \
                                         tg_xxx_zzzzzz_1, tg_xxx_zzzzzzz_0, tg_xxx_zzzzzzz_1, tg_xxxx_xxxxxxx_0, \
                                         tg_xxxx_xxxxxxy_0, tg_xxxx_xxxxxxz_0, tg_xxxx_xxxxxyy_0, tg_xxxx_xxxxxyz_0, \
                                         tg_xxxx_xxxxxzz_0, tg_xxxx_xxxxyyy_0, tg_xxxx_xxxxyyz_0, tg_xxxx_xxxxyzz_0, \
                                         tg_xxxx_xxxxzzz_0, tg_xxxx_xxxyyyy_0, tg_xxxx_xxxyyyz_0, tg_xxxx_xxxyyzz_0, \
                                         tg_xxxx_xxxyzzz_0, tg_xxxx_xxxzzzz_0, tg_xxxx_xxyyyyy_0, tg_xxxx_xxyyyyz_0, \
                                         tg_xxxx_xxyyyzz_0, tg_xxxx_xxyyzzz_0, tg_xxxx_xxyzzzz_0, tg_xxxx_xxzzzzz_0, \
                                         tg_xxxx_xyyyyyy_0, tg_xxxx_xyyyyyz_0, tg_xxxx_xyyyyzz_0, tg_xxxx_xyyyzzz_0, \
                                         tg_xxxx_xyyzzzz_0, tg_xxxx_xyzzzzz_0, tg_xxxx_xzzzzzz_0, tg_xxxx_yyyyyyy_0, \
                                         tg_xxxx_yyyyyyz_0, tg_xxxx_yyyyyzz_0, tg_xxxx_yyyyzzz_0, tg_xxxx_yyyzzzz_0, \
                                         tg_xxxx_yyzzzzz_0, tg_xxxx_yzzzzzz_0, tg_xxxx_zzzzzzz_0, tg_xxxy_xxxxxxx_0, \
                                         tg_xxxy_xxxxxxy_0, tg_xxxy_xxxxxxz_0, tg_xxxy_xxxxxyy_0, tg_xxxy_xxxxxyz_0, \
                                         tg_xxxy_xxxxxzz_0, tg_xxxy_xxxxyyy_0, tg_xxxy_xxxxyyz_0, tg_xxxy_xxxxyzz_0, \
                                         tg_xxxy_xxxxzzz_0, tg_xxxy_xxxyyyy_0, tg_xxxy_xxxyyyz_0, tg_xxxy_xxxyyzz_0, \
                                         tg_xxxy_xxxyzzz_0, tg_xxxy_xxxzzzz_0, tg_xxxy_xxyyyyy_0, tg_xxxy_xxyyyyz_0, \
                                         tg_xxxy_xxyyyzz_0, tg_xxxy_xxyyzzz_0, tg_xxxy_xxyzzzz_0, tg_xxxy_xxzzzzz_0, \
                                         tg_xxxy_xyyyyyy_0, tg_xxxy_xyyyyyz_0, tg_xxxy_xyyyyzz_0, tg_xxxy_xyyyzzz_0, \
                                         tg_xxxy_xyyzzzz_0, tg_xxxy_xyzzzzz_0, tg_xxxy_xzzzzzz_0, tg_xxxy_yyyyyyy_0, \
                                         tg_xxxy_yyyyyyz_0, tg_xxxy_yyyyyzz_0, tg_xxxy_yyyyzzz_0, tg_xxxy_yyyzzzz_0, \
                                         tg_xxxy_yyzzzzz_0, tg_xxxy_yzzzzzz_0, tg_xxxy_zzzzzzz_0, tg_xxxz_xxxxxxx_0, \
                                         tg_xxxz_xxxxxxy_0, tg_xxxz_xxxxxxz_0, tg_xxxz_xxxxxyy_0, tg_xxxz_xxxxxyz_0, \
                                         tg_xxxz_xxxxxzz_0, tg_xxxz_xxxxyyy_0, tg_xxxz_xxxxyyz_0, tg_xxxz_xxxxyzz_0, \
                                         tg_xxxz_xxxxzzz_0, tg_xxxz_xxxyyyy_0, tg_xxxz_xxxyyyz_0, tg_xxxz_xxxyyzz_0, \
                                         tg_xxxz_xxxyzzz_0, tg_xxxz_xxxzzzz_0, tg_xxxz_xxyyyyy_0, tg_xxxz_xxyyyyz_0, \
                                         tg_xxxz_xxyyyzz_0, tg_xxy_xxxxxx_1, tg_xxy_xxxxxxx_0, tg_xxy_xxxxxxx_1, \
                                         tg_xxy_xxxxxxy_0, tg_xxy_xxxxxxy_1, tg_xxy_xxxxxxz_0, tg_xxy_xxxxxxz_1, \
                                         tg_xxy_xxxxxy_1, tg_xxy_xxxxxyy_0, tg_xxy_xxxxxyy_1, tg_xxy_xxxxxyz_0, \
                                         tg_xxy_xxxxxyz_1, tg_xxy_xxxxxz_1, tg_xxy_xxxxxzz_0, tg_xxy_xxxxxzz_1, \
                                         tg_xxy_xxxxyy_1, tg_xxy_xxxxyyy_0, tg_xxy_xxxxyyy_1, tg_xxy_xxxxyyz_0, \
                                         tg_xxy_xxxxyyz_1, tg_xxy_xxxxyz_1, tg_xxy_xxxxyzz_0, tg_xxy_xxxxyzz_1, \
                                         tg_xxy_xxxxzz_1, tg_xxy_xxxxzzz_0, tg_xxy_xxxxzzz_1, tg_xxy_xxxyyy_1, \
                                         tg_xxy_xxxyyyy_0, tg_xxy_xxxyyyy_1, tg_xxy_xxxyyyz_0, tg_xxy_xxxyyyz_1, \
                                         tg_xxy_xxxyyz_1, tg_xxy_xxxyyzz_0, tg_xxy_xxxyyzz_1, tg_xxy_xxxyzz_1, \
                                         tg_xxy_xxxyzzz_0, tg_xxy_xxxyzzz_1, tg_xxy_xxxzzz_1, tg_xxy_xxxzzzz_0, \
                                         tg_xxy_xxxzzzz_1, tg_xxy_xxyyyy_1, tg_xxy_xxyyyyy_0, tg_xxy_xxyyyyy_1, \
                                         tg_xxy_xxyyyyz_0, tg_xxy_xxyyyyz_1, tg_xxy_xxyyyz_1, tg_xxy_xxyyyzz_0, \
                                         tg_xxy_xxyyyzz_1, tg_xxy_xxyyzz_1, tg_xxy_xxyyzzz_0, tg_xxy_xxyyzzz_1, \
                                         tg_xxy_xxyzzz_1, tg_xxy_xxyzzzz_0, tg_xxy_xxyzzzz_1, tg_xxy_xxzzzz_1, \
                                         tg_xxy_xxzzzzz_0, tg_xxy_xxzzzzz_1, tg_xxy_xyyyyy_1, tg_xxy_xyyyyyy_0, \
                                         tg_xxy_xyyyyyy_1, tg_xxy_xyyyyyz_0, tg_xxy_xyyyyyz_1, tg_xxy_xyyyyz_1, \
                                         tg_xxy_xyyyyzz_0, tg_xxy_xyyyyzz_1, tg_xxy_xyyyzz_1, tg_xxy_xyyyzzz_0, \
                                         tg_xxy_xyyyzzz_1, tg_xxy_xyyzzz_1, tg_xxy_xyyzzzz_0, tg_xxy_xyyzzzz_1, \
                                         tg_xxy_xyzzzz_1, tg_xxy_xyzzzzz_0, tg_xxy_xyzzzzz_1, tg_xxy_xzzzzz_1, \
                                         tg_xxy_xzzzzzz_0, tg_xxy_xzzzzzz_1, tg_xxy_yyyyyy_1, tg_xxy_yyyyyyy_0, \
                                         tg_xxy_yyyyyyy_1, tg_xxy_yyyyyyz_0, tg_xxy_yyyyyyz_1, tg_xxy_yyyyyz_1, \
                                         tg_xxy_yyyyyzz_0, tg_xxy_yyyyyzz_1, tg_xxy_yyyyzz_1, tg_xxy_yyyyzzz_0, \
                                         tg_xxy_yyyyzzz_1, tg_xxy_yyyzzz_1, tg_xxy_yyyzzzz_0, tg_xxy_yyyzzzz_1, \
                                         tg_xxy_yyzzzz_1, tg_xxy_yyzzzzz_0, tg_xxy_yyzzzzz_1, tg_xxy_yzzzzz_1, \
                                         tg_xxy_yzzzzzz_0, tg_xxy_yzzzzzz_1, tg_xxy_zzzzzz_1, tg_xxy_zzzzzzz_0, \
                                         tg_xxy_zzzzzzz_1, tg_xxz_xxxxxx_1, tg_xxz_xxxxxxx_0, tg_xxz_xxxxxxx_1, \
                                         tg_xxz_xxxxxxy_0, tg_xxz_xxxxxxy_1, tg_xxz_xxxxxxz_0, tg_xxz_xxxxxxz_1, \
                                         tg_xxz_xxxxxy_1, tg_xxz_xxxxxyy_0, tg_xxz_xxxxxyy_1, tg_xxz_xxxxxyz_0, \
                                         tg_xxz_xxxxxyz_1, tg_xxz_xxxxxz_1, tg_xxz_xxxxxzz_0, tg_xxz_xxxxxzz_1, \
                                         tg_xxz_xxxxyy_1, tg_xxz_xxxxyyy_0, tg_xxz_xxxxyyy_1, tg_xxz_xxxxyyz_0, \
                                         tg_xxz_xxxxyyz_1, tg_xxz_xxxxyz_1, tg_xxz_xxxxyzz_0, tg_xxz_xxxxyzz_1, \
                                         tg_xxz_xxxxzz_1, tg_xxz_xxxxzzz_0, tg_xxz_xxxxzzz_1, tg_xxz_xxxyyy_1, \
                                         tg_xxz_xxxyyyy_0, tg_xxz_xxxyyyy_1, tg_xxz_xxxyyyz_0, tg_xxz_xxxyyyz_1, \
                                         tg_xxz_xxxyyz_1, tg_xxz_xxxyyzz_0, tg_xxz_xxxyyzz_1, tg_xxz_xxxyzz_1, \
                                         tg_xxz_xxxyzzz_0, tg_xxz_xxxyzzz_1, tg_xxz_xxxzzz_1, tg_xxz_xxxzzzz_0, \
                                         tg_xxz_xxxzzzz_1, tg_xxz_xxyyyy_1, tg_xxz_xxyyyyy_0, tg_xxz_xxyyyyy_1, \
                                         tg_xxz_xxyyyyz_0, tg_xxz_xxyyyyz_1, tg_xxz_xxyyyz_1, tg_xxz_xxyyyzz_0, \
                                         tg_xxz_xxyyyzz_1, tg_xxz_xxyyzz_1, tg_xxz_xxyzzz_1, tg_xxz_xxzzzz_1, tg_xxz_xyyyyy_1, \
                                         tg_xxz_xyyyyz_1, tg_xxz_xyyyzz_1, tg_xy_xxxxxxx_0, tg_xy_xxxxxxx_1, tg_xy_xxxxxxy_0, \
                                         tg_xy_xxxxxxy_1, tg_xy_xxxxxxz_0, tg_xy_xxxxxxz_1, tg_xy_xxxxxyy_0, tg_xy_xxxxxyy_1, \
                                         tg_xy_xxxxxyz_0, tg_xy_xxxxxyz_1, tg_xy_xxxxxzz_0, tg_xy_xxxxxzz_1, tg_xy_xxxxyyy_0, \
                                         tg_xy_xxxxyyy_1, tg_xy_xxxxyyz_0, tg_xy_xxxxyyz_1, tg_xy_xxxxyzz_0, tg_xy_xxxxyzz_1, \
                                         tg_xy_xxxxzzz_0, tg_xy_xxxxzzz_1, tg_xy_xxxyyyy_0, tg_xy_xxxyyyy_1, tg_xy_xxxyyyz_0, \
                                         tg_xy_xxxyyyz_1, tg_xy_xxxyyzz_0, tg_xy_xxxyyzz_1, tg_xy_xxxyzzz_0, tg_xy_xxxyzzz_1, \
                                         tg_xy_xxxzzzz_0, tg_xy_xxxzzzz_1, tg_xy_xxyyyyy_0, tg_xy_xxyyyyy_1, tg_xy_xxyyyyz_0, \
                                         tg_xy_xxyyyyz_1, tg_xy_xxyyyzz_0, tg_xy_xxyyyzz_1, tg_xy_xxyyzzz_0, tg_xy_xxyyzzz_1, \
                                         tg_xy_xxyzzzz_0, tg_xy_xxyzzzz_1, tg_xy_xxzzzzz_0, tg_xy_xxzzzzz_1, tg_xy_xyyyyyy_0, \
                                         tg_xy_xyyyyyy_1, tg_xy_xyyyyyz_0, tg_xy_xyyyyyz_1, tg_xy_xyyyyzz_0, tg_xy_xyyyyzz_1, \
                                         tg_xy_xyyyzzz_0, tg_xy_xyyyzzz_1, tg_xy_xyyzzzz_0, tg_xy_xyyzzzz_1, tg_xy_xyzzzzz_0, \
                                         tg_xy_xyzzzzz_1, tg_xy_xzzzzzz_0, tg_xy_xzzzzzz_1, tg_xy_yyyyyyy_0, tg_xy_yyyyyyy_1, \
                                         tg_xy_yyyyyyz_0, tg_xy_yyyyyyz_1, tg_xy_yyyyyzz_0, tg_xy_yyyyyzz_1, tg_xy_yyyyzzz_0, \
                                         tg_xy_yyyyzzz_1, tg_xy_yyyzzzz_0, tg_xy_yyyzzzz_1, tg_xy_yyzzzzz_0, tg_xy_yyzzzzz_1, \
                                         tg_xy_yzzzzzz_0, tg_xy_yzzzzzz_1, tg_xy_zzzzzzz_0, tg_xy_zzzzzzz_1, tg_xz_xxxxxxx_0, \
                                         tg_xz_xxxxxxx_1, tg_xz_xxxxxxy_0, tg_xz_xxxxxxy_1, tg_xz_xxxxxxz_0, tg_xz_xxxxxxz_1, \
                                         tg_xz_xxxxxyy_0, tg_xz_xxxxxyy_1, tg_xz_xxxxxyz_0, tg_xz_xxxxxyz_1, tg_xz_xxxxxzz_0, \
                                         tg_xz_xxxxxzz_1, tg_xz_xxxxyyy_0, tg_xz_xxxxyyy_1, tg_xz_xxxxyyz_0, tg_xz_xxxxyyz_1, \
                                         tg_xz_xxxxyzz_0, tg_xz_xxxxyzz_1, tg_xz_xxxxzzz_0, tg_xz_xxxxzzz_1, tg_xz_xxxyyyy_0, \
                                         tg_xz_xxxyyyy_1, tg_xz_xxxyyyz_0, tg_xz_xxxyyyz_1, tg_xz_xxxyyzz_0, tg_xz_xxxyyzz_1, \
                                         tg_xz_xxxyzzz_0, tg_xz_xxxyzzz_1, tg_xz_xxxzzzz_0, tg_xz_xxxzzzz_1, tg_xz_xxyyyyy_0, \
                                         tg_xz_xxyyyyy_1, tg_xz_xxyyyyz_0, tg_xz_xxyyyyz_1, tg_xz_xxyyyzz_0, tg_xz_xxyyyzz_1, \
                                         wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxx_xxxxxxx_0[j] = pb_x * tg_xxx_xxxxxxx_0[j] + fr * tg_xxx_xxxxxxx_1[j] + 1.5 * fl1_fx * (tg_xx_xxxxxxx_0[j] - tg_xx_xxxxxxx_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_xxx_xxxxxx_1[j];

                    tg_xxxx_xxxxxxy_0[j] = pb_x * tg_xxx_xxxxxxy_0[j] + fr * tg_xxx_xxxxxxy_1[j] + 1.5 * fl1_fx * (tg_xx_xxxxxxy_0[j] - tg_xx_xxxxxxy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxx_xxxxxy_1[j];

                    tg_xxxx_xxxxxxz_0[j] = pb_x * tg_xxx_xxxxxxz_0[j] + fr * tg_xxx_xxxxxxz_1[j] + 1.5 * fl1_fx * (tg_xx_xxxxxxz_0[j] - tg_xx_xxxxxxz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxx_xxxxxz_1[j];

                    tg_xxxx_xxxxxyy_0[j] = pb_x * tg_xxx_xxxxxyy_0[j] + fr * tg_xxx_xxxxxyy_1[j] + 1.5 * fl1_fx * (tg_xx_xxxxxyy_0[j] - tg_xx_xxxxxyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxx_xxxxyy_1[j];

                    tg_xxxx_xxxxxyz_0[j] = pb_x * tg_xxx_xxxxxyz_0[j] + fr * tg_xxx_xxxxxyz_1[j] + 1.5 * fl1_fx * (tg_xx_xxxxxyz_0[j] - tg_xx_xxxxxyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxx_xxxxyz_1[j];

                    tg_xxxx_xxxxxzz_0[j] = pb_x * tg_xxx_xxxxxzz_0[j] + fr * tg_xxx_xxxxxzz_1[j] + 1.5 * fl1_fx * (tg_xx_xxxxxzz_0[j] - tg_xx_xxxxxzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxx_xxxxzz_1[j];

                    tg_xxxx_xxxxyyy_0[j] = pb_x * tg_xxx_xxxxyyy_0[j] + fr * tg_xxx_xxxxyyy_1[j] + 1.5 * fl1_fx * (tg_xx_xxxxyyy_0[j] - tg_xx_xxxxyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxx_xxxyyy_1[j];

                    tg_xxxx_xxxxyyz_0[j] = pb_x * tg_xxx_xxxxyyz_0[j] + fr * tg_xxx_xxxxyyz_1[j] + 1.5 * fl1_fx * (tg_xx_xxxxyyz_0[j] - tg_xx_xxxxyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxx_xxxyyz_1[j];

                    tg_xxxx_xxxxyzz_0[j] = pb_x * tg_xxx_xxxxyzz_0[j] + fr * tg_xxx_xxxxyzz_1[j] + 1.5 * fl1_fx * (tg_xx_xxxxyzz_0[j] - tg_xx_xxxxyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxx_xxxyzz_1[j];

                    tg_xxxx_xxxxzzz_0[j] = pb_x * tg_xxx_xxxxzzz_0[j] + fr * tg_xxx_xxxxzzz_1[j] + 1.5 * fl1_fx * (tg_xx_xxxxzzz_0[j] - tg_xx_xxxxzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxx_xxxzzz_1[j];

                    tg_xxxx_xxxyyyy_0[j] = pb_x * tg_xxx_xxxyyyy_0[j] + fr * tg_xxx_xxxyyyy_1[j] + 1.5 * fl1_fx * (tg_xx_xxxyyyy_0[j] - tg_xx_xxxyyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxx_xxyyyy_1[j];

                    tg_xxxx_xxxyyyz_0[j] = pb_x * tg_xxx_xxxyyyz_0[j] + fr * tg_xxx_xxxyyyz_1[j] + 1.5 * fl1_fx * (tg_xx_xxxyyyz_0[j] - tg_xx_xxxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxx_xxyyyz_1[j];

                    tg_xxxx_xxxyyzz_0[j] = pb_x * tg_xxx_xxxyyzz_0[j] + fr * tg_xxx_xxxyyzz_1[j] + 1.5 * fl1_fx * (tg_xx_xxxyyzz_0[j] - tg_xx_xxxyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxx_xxyyzz_1[j];

                    tg_xxxx_xxxyzzz_0[j] = pb_x * tg_xxx_xxxyzzz_0[j] + fr * tg_xxx_xxxyzzz_1[j] + 1.5 * fl1_fx * (tg_xx_xxxyzzz_0[j] - tg_xx_xxxyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxx_xxyzzz_1[j];

                    tg_xxxx_xxxzzzz_0[j] = pb_x * tg_xxx_xxxzzzz_0[j] + fr * tg_xxx_xxxzzzz_1[j] + 1.5 * fl1_fx * (tg_xx_xxxzzzz_0[j] - tg_xx_xxxzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxx_xxzzzz_1[j];

                    tg_xxxx_xxyyyyy_0[j] = pb_x * tg_xxx_xxyyyyy_0[j] + fr * tg_xxx_xxyyyyy_1[j] + 1.5 * fl1_fx * (tg_xx_xxyyyyy_0[j] - tg_xx_xxyyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxx_xyyyyy_1[j];

                    tg_xxxx_xxyyyyz_0[j] = pb_x * tg_xxx_xxyyyyz_0[j] + fr * tg_xxx_xxyyyyz_1[j] + 1.5 * fl1_fx * (tg_xx_xxyyyyz_0[j] - tg_xx_xxyyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxx_xyyyyz_1[j];

                    tg_xxxx_xxyyyzz_0[j] = pb_x * tg_xxx_xxyyyzz_0[j] + fr * tg_xxx_xxyyyzz_1[j] + 1.5 * fl1_fx * (tg_xx_xxyyyzz_0[j] - tg_xx_xxyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxx_xyyyzz_1[j];

                    tg_xxxx_xxyyzzz_0[j] = pb_x * tg_xxx_xxyyzzz_0[j] + fr * tg_xxx_xxyyzzz_1[j] + 1.5 * fl1_fx * (tg_xx_xxyyzzz_0[j] - tg_xx_xxyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxx_xyyzzz_1[j];

                    tg_xxxx_xxyzzzz_0[j] = pb_x * tg_xxx_xxyzzzz_0[j] + fr * tg_xxx_xxyzzzz_1[j] + 1.5 * fl1_fx * (tg_xx_xxyzzzz_0[j] - tg_xx_xxyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxx_xyzzzz_1[j];

                    tg_xxxx_xxzzzzz_0[j] = pb_x * tg_xxx_xxzzzzz_0[j] + fr * tg_xxx_xxzzzzz_1[j] + 1.5 * fl1_fx * (tg_xx_xxzzzzz_0[j] - tg_xx_xxzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxx_xzzzzz_1[j];

                    tg_xxxx_xyyyyyy_0[j] = pb_x * tg_xxx_xyyyyyy_0[j] + fr * tg_xxx_xyyyyyy_1[j] + 1.5 * fl1_fx * (tg_xx_xyyyyyy_0[j] - tg_xx_xyyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxx_yyyyyy_1[j];

                    tg_xxxx_xyyyyyz_0[j] = pb_x * tg_xxx_xyyyyyz_0[j] + fr * tg_xxx_xyyyyyz_1[j] + 1.5 * fl1_fx * (tg_xx_xyyyyyz_0[j] - tg_xx_xyyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxx_yyyyyz_1[j];

                    tg_xxxx_xyyyyzz_0[j] = pb_x * tg_xxx_xyyyyzz_0[j] + fr * tg_xxx_xyyyyzz_1[j] + 1.5 * fl1_fx * (tg_xx_xyyyyzz_0[j] - tg_xx_xyyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxx_yyyyzz_1[j];

                    tg_xxxx_xyyyzzz_0[j] = pb_x * tg_xxx_xyyyzzz_0[j] + fr * tg_xxx_xyyyzzz_1[j] + 1.5 * fl1_fx * (tg_xx_xyyyzzz_0[j] - tg_xx_xyyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxx_yyyzzz_1[j];

                    tg_xxxx_xyyzzzz_0[j] = pb_x * tg_xxx_xyyzzzz_0[j] + fr * tg_xxx_xyyzzzz_1[j] + 1.5 * fl1_fx * (tg_xx_xyyzzzz_0[j] - tg_xx_xyyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxx_yyzzzz_1[j];

                    tg_xxxx_xyzzzzz_0[j] = pb_x * tg_xxx_xyzzzzz_0[j] + fr * tg_xxx_xyzzzzz_1[j] + 1.5 * fl1_fx * (tg_xx_xyzzzzz_0[j] - tg_xx_xyzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxx_yzzzzz_1[j];

                    tg_xxxx_xzzzzzz_0[j] = pb_x * tg_xxx_xzzzzzz_0[j] + fr * tg_xxx_xzzzzzz_1[j] + 1.5 * fl1_fx * (tg_xx_xzzzzzz_0[j] - tg_xx_xzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxx_zzzzzz_1[j];

                    tg_xxxx_yyyyyyy_0[j] = pb_x * tg_xxx_yyyyyyy_0[j] + fr * tg_xxx_yyyyyyy_1[j] + 1.5 * fl1_fx * (tg_xx_yyyyyyy_0[j] - tg_xx_yyyyyyy_1[j] * fl1_fza);

                    tg_xxxx_yyyyyyz_0[j] = pb_x * tg_xxx_yyyyyyz_0[j] + fr * tg_xxx_yyyyyyz_1[j] + 1.5 * fl1_fx * (tg_xx_yyyyyyz_0[j] - tg_xx_yyyyyyz_1[j] * fl1_fza);

                    tg_xxxx_yyyyyzz_0[j] = pb_x * tg_xxx_yyyyyzz_0[j] + fr * tg_xxx_yyyyyzz_1[j] + 1.5 * fl1_fx * (tg_xx_yyyyyzz_0[j] - tg_xx_yyyyyzz_1[j] * fl1_fza);

                    tg_xxxx_yyyyzzz_0[j] = pb_x * tg_xxx_yyyyzzz_0[j] + fr * tg_xxx_yyyyzzz_1[j] + 1.5 * fl1_fx * (tg_xx_yyyyzzz_0[j] - tg_xx_yyyyzzz_1[j] * fl1_fza);

                    tg_xxxx_yyyzzzz_0[j] = pb_x * tg_xxx_yyyzzzz_0[j] + fr * tg_xxx_yyyzzzz_1[j] + 1.5 * fl1_fx * (tg_xx_yyyzzzz_0[j] - tg_xx_yyyzzzz_1[j] * fl1_fza);

                    tg_xxxx_yyzzzzz_0[j] = pb_x * tg_xxx_yyzzzzz_0[j] + fr * tg_xxx_yyzzzzz_1[j] + 1.5 * fl1_fx * (tg_xx_yyzzzzz_0[j] - tg_xx_yyzzzzz_1[j] * fl1_fza);

                    tg_xxxx_yzzzzzz_0[j] = pb_x * tg_xxx_yzzzzzz_0[j] + fr * tg_xxx_yzzzzzz_1[j] + 1.5 * fl1_fx * (tg_xx_yzzzzzz_0[j] - tg_xx_yzzzzzz_1[j] * fl1_fza);

                    tg_xxxx_zzzzzzz_0[j] = pb_x * tg_xxx_zzzzzzz_0[j] + fr * tg_xxx_zzzzzzz_1[j] + 1.5 * fl1_fx * (tg_xx_zzzzzzz_0[j] - tg_xx_zzzzzzz_1[j] * fl1_fza);

                    tg_xxxy_xxxxxxx_0[j] = pb_x * tg_xxy_xxxxxxx_0[j] + fr * tg_xxy_xxxxxxx_1[j] + fl1_fx * (tg_xy_xxxxxxx_0[j] - tg_xy_xxxxxxx_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_xxy_xxxxxx_1[j];

                    tg_xxxy_xxxxxxy_0[j] = pb_x * tg_xxy_xxxxxxy_0[j] + fr * tg_xxy_xxxxxxy_1[j] + fl1_fx * (tg_xy_xxxxxxy_0[j] - tg_xy_xxxxxxy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxy_xxxxxy_1[j];

                    tg_xxxy_xxxxxxz_0[j] = pb_x * tg_xxy_xxxxxxz_0[j] + fr * tg_xxy_xxxxxxz_1[j] + fl1_fx * (tg_xy_xxxxxxz_0[j] - tg_xy_xxxxxxz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxy_xxxxxz_1[j];

                    tg_xxxy_xxxxxyy_0[j] = pb_x * tg_xxy_xxxxxyy_0[j] + fr * tg_xxy_xxxxxyy_1[j] + fl1_fx * (tg_xy_xxxxxyy_0[j] - tg_xy_xxxxxyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxy_xxxxyy_1[j];

                    tg_xxxy_xxxxxyz_0[j] = pb_x * tg_xxy_xxxxxyz_0[j] + fr * tg_xxy_xxxxxyz_1[j] + fl1_fx * (tg_xy_xxxxxyz_0[j] - tg_xy_xxxxxyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxy_xxxxyz_1[j];

                    tg_xxxy_xxxxxzz_0[j] = pb_x * tg_xxy_xxxxxzz_0[j] + fr * tg_xxy_xxxxxzz_1[j] + fl1_fx * (tg_xy_xxxxxzz_0[j] - tg_xy_xxxxxzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxy_xxxxzz_1[j];

                    tg_xxxy_xxxxyyy_0[j] = pb_x * tg_xxy_xxxxyyy_0[j] + fr * tg_xxy_xxxxyyy_1[j] + fl1_fx * (tg_xy_xxxxyyy_0[j] - tg_xy_xxxxyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxy_xxxyyy_1[j];

                    tg_xxxy_xxxxyyz_0[j] = pb_x * tg_xxy_xxxxyyz_0[j] + fr * tg_xxy_xxxxyyz_1[j] + fl1_fx * (tg_xy_xxxxyyz_0[j] - tg_xy_xxxxyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxy_xxxyyz_1[j];

                    tg_xxxy_xxxxyzz_0[j] = pb_x * tg_xxy_xxxxyzz_0[j] + fr * tg_xxy_xxxxyzz_1[j] + fl1_fx * (tg_xy_xxxxyzz_0[j] - tg_xy_xxxxyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxy_xxxyzz_1[j];

                    tg_xxxy_xxxxzzz_0[j] = pb_x * tg_xxy_xxxxzzz_0[j] + fr * tg_xxy_xxxxzzz_1[j] + fl1_fx * (tg_xy_xxxxzzz_0[j] - tg_xy_xxxxzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxy_xxxzzz_1[j];

                    tg_xxxy_xxxyyyy_0[j] = pb_x * tg_xxy_xxxyyyy_0[j] + fr * tg_xxy_xxxyyyy_1[j] + fl1_fx * (tg_xy_xxxyyyy_0[j] - tg_xy_xxxyyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxy_xxyyyy_1[j];

                    tg_xxxy_xxxyyyz_0[j] = pb_x * tg_xxy_xxxyyyz_0[j] + fr * tg_xxy_xxxyyyz_1[j] + fl1_fx * (tg_xy_xxxyyyz_0[j] - tg_xy_xxxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxy_xxyyyz_1[j];

                    tg_xxxy_xxxyyzz_0[j] = pb_x * tg_xxy_xxxyyzz_0[j] + fr * tg_xxy_xxxyyzz_1[j] + fl1_fx * (tg_xy_xxxyyzz_0[j] - tg_xy_xxxyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxy_xxyyzz_1[j];

                    tg_xxxy_xxxyzzz_0[j] = pb_x * tg_xxy_xxxyzzz_0[j] + fr * tg_xxy_xxxyzzz_1[j] + fl1_fx * (tg_xy_xxxyzzz_0[j] - tg_xy_xxxyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxy_xxyzzz_1[j];

                    tg_xxxy_xxxzzzz_0[j] = pb_x * tg_xxy_xxxzzzz_0[j] + fr * tg_xxy_xxxzzzz_1[j] + fl1_fx * (tg_xy_xxxzzzz_0[j] - tg_xy_xxxzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxy_xxzzzz_1[j];

                    tg_xxxy_xxyyyyy_0[j] = pb_x * tg_xxy_xxyyyyy_0[j] + fr * tg_xxy_xxyyyyy_1[j] + fl1_fx * (tg_xy_xxyyyyy_0[j] - tg_xy_xxyyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxy_xyyyyy_1[j];

                    tg_xxxy_xxyyyyz_0[j] = pb_x * tg_xxy_xxyyyyz_0[j] + fr * tg_xxy_xxyyyyz_1[j] + fl1_fx * (tg_xy_xxyyyyz_0[j] - tg_xy_xxyyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxy_xyyyyz_1[j];

                    tg_xxxy_xxyyyzz_0[j] = pb_x * tg_xxy_xxyyyzz_0[j] + fr * tg_xxy_xxyyyzz_1[j] + fl1_fx * (tg_xy_xxyyyzz_0[j] - tg_xy_xxyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxy_xyyyzz_1[j];

                    tg_xxxy_xxyyzzz_0[j] = pb_x * tg_xxy_xxyyzzz_0[j] + fr * tg_xxy_xxyyzzz_1[j] + fl1_fx * (tg_xy_xxyyzzz_0[j] - tg_xy_xxyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxy_xyyzzz_1[j];

                    tg_xxxy_xxyzzzz_0[j] = pb_x * tg_xxy_xxyzzzz_0[j] + fr * tg_xxy_xxyzzzz_1[j] + fl1_fx * (tg_xy_xxyzzzz_0[j] - tg_xy_xxyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxy_xyzzzz_1[j];

                    tg_xxxy_xxzzzzz_0[j] = pb_x * tg_xxy_xxzzzzz_0[j] + fr * tg_xxy_xxzzzzz_1[j] + fl1_fx * (tg_xy_xxzzzzz_0[j] - tg_xy_xxzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxy_xzzzzz_1[j];

                    tg_xxxy_xyyyyyy_0[j] = pb_x * tg_xxy_xyyyyyy_0[j] + fr * tg_xxy_xyyyyyy_1[j] + fl1_fx * (tg_xy_xyyyyyy_0[j] - tg_xy_xyyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxy_yyyyyy_1[j];

                    tg_xxxy_xyyyyyz_0[j] = pb_x * tg_xxy_xyyyyyz_0[j] + fr * tg_xxy_xyyyyyz_1[j] + fl1_fx * (tg_xy_xyyyyyz_0[j] - tg_xy_xyyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxy_yyyyyz_1[j];

                    tg_xxxy_xyyyyzz_0[j] = pb_x * tg_xxy_xyyyyzz_0[j] + fr * tg_xxy_xyyyyzz_1[j] + fl1_fx * (tg_xy_xyyyyzz_0[j] - tg_xy_xyyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxy_yyyyzz_1[j];

                    tg_xxxy_xyyyzzz_0[j] = pb_x * tg_xxy_xyyyzzz_0[j] + fr * tg_xxy_xyyyzzz_1[j] + fl1_fx * (tg_xy_xyyyzzz_0[j] - tg_xy_xyyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxy_yyyzzz_1[j];

                    tg_xxxy_xyyzzzz_0[j] = pb_x * tg_xxy_xyyzzzz_0[j] + fr * tg_xxy_xyyzzzz_1[j] + fl1_fx * (tg_xy_xyyzzzz_0[j] - tg_xy_xyyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxy_yyzzzz_1[j];

                    tg_xxxy_xyzzzzz_0[j] = pb_x * tg_xxy_xyzzzzz_0[j] + fr * tg_xxy_xyzzzzz_1[j] + fl1_fx * (tg_xy_xyzzzzz_0[j] - tg_xy_xyzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxy_yzzzzz_1[j];

                    tg_xxxy_xzzzzzz_0[j] = pb_x * tg_xxy_xzzzzzz_0[j] + fr * tg_xxy_xzzzzzz_1[j] + fl1_fx * (tg_xy_xzzzzzz_0[j] - tg_xy_xzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxy_zzzzzz_1[j];

                    tg_xxxy_yyyyyyy_0[j] = pb_x * tg_xxy_yyyyyyy_0[j] + fr * tg_xxy_yyyyyyy_1[j] + fl1_fx * (tg_xy_yyyyyyy_0[j] - tg_xy_yyyyyyy_1[j] * fl1_fza);

                    tg_xxxy_yyyyyyz_0[j] = pb_x * tg_xxy_yyyyyyz_0[j] + fr * tg_xxy_yyyyyyz_1[j] + fl1_fx * (tg_xy_yyyyyyz_0[j] - tg_xy_yyyyyyz_1[j] * fl1_fza);

                    tg_xxxy_yyyyyzz_0[j] = pb_x * tg_xxy_yyyyyzz_0[j] + fr * tg_xxy_yyyyyzz_1[j] + fl1_fx * (tg_xy_yyyyyzz_0[j] - tg_xy_yyyyyzz_1[j] * fl1_fza);

                    tg_xxxy_yyyyzzz_0[j] = pb_x * tg_xxy_yyyyzzz_0[j] + fr * tg_xxy_yyyyzzz_1[j] + fl1_fx * (tg_xy_yyyyzzz_0[j] - tg_xy_yyyyzzz_1[j] * fl1_fza);

                    tg_xxxy_yyyzzzz_0[j] = pb_x * tg_xxy_yyyzzzz_0[j] + fr * tg_xxy_yyyzzzz_1[j] + fl1_fx * (tg_xy_yyyzzzz_0[j] - tg_xy_yyyzzzz_1[j] * fl1_fza);

                    tg_xxxy_yyzzzzz_0[j] = pb_x * tg_xxy_yyzzzzz_0[j] + fr * tg_xxy_yyzzzzz_1[j] + fl1_fx * (tg_xy_yyzzzzz_0[j] - tg_xy_yyzzzzz_1[j] * fl1_fza);

                    tg_xxxy_yzzzzzz_0[j] = pb_x * tg_xxy_yzzzzzz_0[j] + fr * tg_xxy_yzzzzzz_1[j] + fl1_fx * (tg_xy_yzzzzzz_0[j] - tg_xy_yzzzzzz_1[j] * fl1_fza);

                    tg_xxxy_zzzzzzz_0[j] = pb_x * tg_xxy_zzzzzzz_0[j] + fr * tg_xxy_zzzzzzz_1[j] + fl1_fx * (tg_xy_zzzzzzz_0[j] - tg_xy_zzzzzzz_1[j] * fl1_fza);

                    tg_xxxz_xxxxxxx_0[j] = pb_x * tg_xxz_xxxxxxx_0[j] + fr * tg_xxz_xxxxxxx_1[j] + fl1_fx * (tg_xz_xxxxxxx_0[j] - tg_xz_xxxxxxx_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_xxz_xxxxxx_1[j];

                    tg_xxxz_xxxxxxy_0[j] = pb_x * tg_xxz_xxxxxxy_0[j] + fr * tg_xxz_xxxxxxy_1[j] + fl1_fx * (tg_xz_xxxxxxy_0[j] - tg_xz_xxxxxxy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxz_xxxxxy_1[j];

                    tg_xxxz_xxxxxxz_0[j] = pb_x * tg_xxz_xxxxxxz_0[j] + fr * tg_xxz_xxxxxxz_1[j] + fl1_fx * (tg_xz_xxxxxxz_0[j] - tg_xz_xxxxxxz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxz_xxxxxz_1[j];

                    tg_xxxz_xxxxxyy_0[j] = pb_x * tg_xxz_xxxxxyy_0[j] + fr * tg_xxz_xxxxxyy_1[j] + fl1_fx * (tg_xz_xxxxxyy_0[j] - tg_xz_xxxxxyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxz_xxxxyy_1[j];

                    tg_xxxz_xxxxxyz_0[j] = pb_x * tg_xxz_xxxxxyz_0[j] + fr * tg_xxz_xxxxxyz_1[j] + fl1_fx * (tg_xz_xxxxxyz_0[j] - tg_xz_xxxxxyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxz_xxxxyz_1[j];

                    tg_xxxz_xxxxxzz_0[j] = pb_x * tg_xxz_xxxxxzz_0[j] + fr * tg_xxz_xxxxxzz_1[j] + fl1_fx * (tg_xz_xxxxxzz_0[j] - tg_xz_xxxxxzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxz_xxxxzz_1[j];

                    tg_xxxz_xxxxyyy_0[j] = pb_x * tg_xxz_xxxxyyy_0[j] + fr * tg_xxz_xxxxyyy_1[j] + fl1_fx * (tg_xz_xxxxyyy_0[j] - tg_xz_xxxxyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxz_xxxyyy_1[j];

                    tg_xxxz_xxxxyyz_0[j] = pb_x * tg_xxz_xxxxyyz_0[j] + fr * tg_xxz_xxxxyyz_1[j] + fl1_fx * (tg_xz_xxxxyyz_0[j] - tg_xz_xxxxyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxz_xxxyyz_1[j];

                    tg_xxxz_xxxxyzz_0[j] = pb_x * tg_xxz_xxxxyzz_0[j] + fr * tg_xxz_xxxxyzz_1[j] + fl1_fx * (tg_xz_xxxxyzz_0[j] - tg_xz_xxxxyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxz_xxxyzz_1[j];

                    tg_xxxz_xxxxzzz_0[j] = pb_x * tg_xxz_xxxxzzz_0[j] + fr * tg_xxz_xxxxzzz_1[j] + fl1_fx * (tg_xz_xxxxzzz_0[j] - tg_xz_xxxxzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxz_xxxzzz_1[j];

                    tg_xxxz_xxxyyyy_0[j] = pb_x * tg_xxz_xxxyyyy_0[j] + fr * tg_xxz_xxxyyyy_1[j] + fl1_fx * (tg_xz_xxxyyyy_0[j] - tg_xz_xxxyyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxz_xxyyyy_1[j];

                    tg_xxxz_xxxyyyz_0[j] = pb_x * tg_xxz_xxxyyyz_0[j] + fr * tg_xxz_xxxyyyz_1[j] + fl1_fx * (tg_xz_xxxyyyz_0[j] - tg_xz_xxxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxz_xxyyyz_1[j];

                    tg_xxxz_xxxyyzz_0[j] = pb_x * tg_xxz_xxxyyzz_0[j] + fr * tg_xxz_xxxyyzz_1[j] + fl1_fx * (tg_xz_xxxyyzz_0[j] - tg_xz_xxxyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxz_xxyyzz_1[j];

                    tg_xxxz_xxxyzzz_0[j] = pb_x * tg_xxz_xxxyzzz_0[j] + fr * tg_xxz_xxxyzzz_1[j] + fl1_fx * (tg_xz_xxxyzzz_0[j] - tg_xz_xxxyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxz_xxyzzz_1[j];

                    tg_xxxz_xxxzzzz_0[j] = pb_x * tg_xxz_xxxzzzz_0[j] + fr * tg_xxz_xxxzzzz_1[j] + fl1_fx * (tg_xz_xxxzzzz_0[j] - tg_xz_xxxzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxz_xxzzzz_1[j];

                    tg_xxxz_xxyyyyy_0[j] = pb_x * tg_xxz_xxyyyyy_0[j] + fr * tg_xxz_xxyyyyy_1[j] + fl1_fx * (tg_xz_xxyyyyy_0[j] - tg_xz_xxyyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxz_xyyyyy_1[j];

                    tg_xxxz_xxyyyyz_0[j] = pb_x * tg_xxz_xxyyyyz_0[j] + fr * tg_xxz_xxyyyyz_1[j] + fl1_fx * (tg_xz_xxyyyyz_0[j] - tg_xz_xxyyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxz_xyyyyz_1[j];

                    tg_xxxz_xxyyyzz_0[j] = pb_x * tg_xxz_xxyyyzz_0[j] + fr * tg_xxz_xxyyyzz_1[j] + fl1_fx * (tg_xz_xxyyyzz_0[j] - tg_xz_xxyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxz_xyyyzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSK_90_180(      CMemBlock2D<double>* primBuffer,
                                        const CRecursionMap&       recursionMap,
                                        const CMemBlock2D<double>& osFactors,
                                        const CMemBlock2D<double>& wpDistances,
                                        const CGtoPairsBlock&      braGtoPairsBlock,
                                        const CGtoPairsBlock&      ketGtoPairsBlock,
                                        const int32_t              nKetPrimPairs,
                                        const int32_t              iContrPair)
    {
        // Batch of Integrals (90,180)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {4, -1, -1, -1},
                                             {7, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_7_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {7, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_7_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_7_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_7_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_2_7_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {7, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_2_7_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {7, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {6, -1, -1, -1}, 
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

                auto tg_xxz_xxyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 90); 

                auto tg_xxz_xxyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 91); 

                auto tg_xxz_xxzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 92); 

                auto tg_xxz_xyyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 93); 

                auto tg_xxz_xyyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 94); 

                auto tg_xxz_xyyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 95); 

                auto tg_xxz_xyyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 96); 

                auto tg_xxz_xyyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 97); 

                auto tg_xxz_xyzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 98); 

                auto tg_xxz_xzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 99); 

                auto tg_xxz_yyyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 100); 

                auto tg_xxz_yyyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 101); 

                auto tg_xxz_yyyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 102); 

                auto tg_xxz_yyyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 103); 

                auto tg_xxz_yyyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 104); 

                auto tg_xxz_yyzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 105); 

                auto tg_xxz_yzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 106); 

                auto tg_xxz_zzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 107); 

                auto tg_xyy_xxxxxxx_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 108); 

                auto tg_xyy_xxxxxxy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 109); 

                auto tg_xyy_xxxxxxz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 110); 

                auto tg_xyy_xxxxxyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 111); 

                auto tg_xyy_xxxxxyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 112); 

                auto tg_xyy_xxxxxzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 113); 

                auto tg_xyy_xxxxyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 114); 

                auto tg_xyy_xxxxyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 115); 

                auto tg_xyy_xxxxyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 116); 

                auto tg_xyy_xxxxzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 117); 

                auto tg_xyy_xxxyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 118); 

                auto tg_xyy_xxxyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 119); 

                auto tg_xyy_xxxyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 120); 

                auto tg_xyy_xxxyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 121); 

                auto tg_xyy_xxxzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 122); 

                auto tg_xyy_xxyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 123); 

                auto tg_xyy_xxyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 124); 

                auto tg_xyy_xxyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 125); 

                auto tg_xyy_xxyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 126); 

                auto tg_xyy_xxyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 127); 

                auto tg_xyy_xxzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 128); 

                auto tg_xyy_xyyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 129); 

                auto tg_xyy_xyyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 130); 

                auto tg_xyy_xyyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 131); 

                auto tg_xyy_xyyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 132); 

                auto tg_xyy_xyyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 133); 

                auto tg_xyy_xyzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 134); 

                auto tg_xyy_xzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 135); 

                auto tg_xyy_yyyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 136); 

                auto tg_xyy_yyyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 137); 

                auto tg_xyy_yyyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 138); 

                auto tg_xyy_yyyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 139); 

                auto tg_xyy_yyyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 140); 

                auto tg_xyy_yyzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 141); 

                auto tg_xyy_yzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 142); 

                auto tg_xyy_zzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 143); 

                auto tg_xyz_xxxxxxx_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 144); 

                auto tg_xyz_xxxxxxy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 145); 

                auto tg_xyz_xxxxxxz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 146); 

                auto tg_xyz_xxxxxyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 147); 

                auto tg_xyz_xxxxxyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 148); 

                auto tg_xyz_xxxxxzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 149); 

                auto tg_xyz_xxxxyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 150); 

                auto tg_xyz_xxxxyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 151); 

                auto tg_xyz_xxxxyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 152); 

                auto tg_xyz_xxxxzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 153); 

                auto tg_xyz_xxxyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 154); 

                auto tg_xyz_xxxyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 155); 

                auto tg_xyz_xxxyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 156); 

                auto tg_xyz_xxxyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 157); 

                auto tg_xyz_xxxzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 158); 

                auto tg_xyz_xxyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 159); 

                auto tg_xyz_xxyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 160); 

                auto tg_xyz_xxyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 161); 

                auto tg_xyz_xxyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 162); 

                auto tg_xyz_xxyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 163); 

                auto tg_xyz_xxzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 164); 

                auto tg_xyz_xyyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 165); 

                auto tg_xyz_xyyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 166); 

                auto tg_xyz_xyyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 167); 

                auto tg_xyz_xyyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 168); 

                auto tg_xyz_xyyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 169); 

                auto tg_xyz_xyzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 170); 

                auto tg_xyz_xzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 171); 

                auto tg_xyz_yyyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 172); 

                auto tg_xyz_yyyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 173); 

                auto tg_xyz_yyyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 174); 

                auto tg_xyz_yyyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 175); 

                auto tg_xyz_yyyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 176); 

                auto tg_xyz_yyzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 177); 

                auto tg_xyz_yzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 178); 

                auto tg_xyz_zzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 179); 

                auto tg_xxz_xxyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 90); 

                auto tg_xxz_xxyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 91); 

                auto tg_xxz_xxzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 92); 

                auto tg_xxz_xyyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 93); 

                auto tg_xxz_xyyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 94); 

                auto tg_xxz_xyyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 95); 

                auto tg_xxz_xyyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 96); 

                auto tg_xxz_xyyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 97); 

                auto tg_xxz_xyzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 98); 

                auto tg_xxz_xzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 99); 

                auto tg_xxz_yyyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 100); 

                auto tg_xxz_yyyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 101); 

                auto tg_xxz_yyyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 102); 

                auto tg_xxz_yyyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 103); 

                auto tg_xxz_yyyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 104); 

                auto tg_xxz_yyzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 105); 

                auto tg_xxz_yzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 106); 

                auto tg_xxz_zzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 107); 

                auto tg_xyy_xxxxxxx_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 108); 

                auto tg_xyy_xxxxxxy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 109); 

                auto tg_xyy_xxxxxxz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 110); 

                auto tg_xyy_xxxxxyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 111); 

                auto tg_xyy_xxxxxyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 112); 

                auto tg_xyy_xxxxxzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 113); 

                auto tg_xyy_xxxxyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 114); 

                auto tg_xyy_xxxxyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 115); 

                auto tg_xyy_xxxxyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 116); 

                auto tg_xyy_xxxxzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 117); 

                auto tg_xyy_xxxyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 118); 

                auto tg_xyy_xxxyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 119); 

                auto tg_xyy_xxxyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 120); 

                auto tg_xyy_xxxyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 121); 

                auto tg_xyy_xxxzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 122); 

                auto tg_xyy_xxyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 123); 

                auto tg_xyy_xxyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 124); 

                auto tg_xyy_xxyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 125); 

                auto tg_xyy_xxyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 126); 

                auto tg_xyy_xxyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 127); 

                auto tg_xyy_xxzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 128); 

                auto tg_xyy_xyyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 129); 

                auto tg_xyy_xyyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 130); 

                auto tg_xyy_xyyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 131); 

                auto tg_xyy_xyyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 132); 

                auto tg_xyy_xyyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 133); 

                auto tg_xyy_xyzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 134); 

                auto tg_xyy_xzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 135); 

                auto tg_xyy_yyyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 136); 

                auto tg_xyy_yyyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 137); 

                auto tg_xyy_yyyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 138); 

                auto tg_xyy_yyyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 139); 

                auto tg_xyy_yyyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 140); 

                auto tg_xyy_yyzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 141); 

                auto tg_xyy_yzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 142); 

                auto tg_xyy_zzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 143); 

                auto tg_xyz_xxxxxxx_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 144); 

                auto tg_xyz_xxxxxxy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 145); 

                auto tg_xyz_xxxxxxz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 146); 

                auto tg_xyz_xxxxxyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 147); 

                auto tg_xyz_xxxxxyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 148); 

                auto tg_xyz_xxxxxzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 149); 

                auto tg_xyz_xxxxyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 150); 

                auto tg_xyz_xxxxyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 151); 

                auto tg_xyz_xxxxyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 152); 

                auto tg_xyz_xxxxzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 153); 

                auto tg_xyz_xxxyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 154); 

                auto tg_xyz_xxxyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 155); 

                auto tg_xyz_xxxyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 156); 

                auto tg_xyz_xxxyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 157); 

                auto tg_xyz_xxxzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 158); 

                auto tg_xyz_xxyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 159); 

                auto tg_xyz_xxyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 160); 

                auto tg_xyz_xxyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 161); 

                auto tg_xyz_xxyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 162); 

                auto tg_xyz_xxyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 163); 

                auto tg_xyz_xxzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 164); 

                auto tg_xyz_xyyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 165); 

                auto tg_xyz_xyyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 166); 

                auto tg_xyz_xyyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 167); 

                auto tg_xyz_xyyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 168); 

                auto tg_xyz_xyyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 169); 

                auto tg_xyz_xyzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 170); 

                auto tg_xyz_xzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 171); 

                auto tg_xyz_yyyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 172); 

                auto tg_xyz_yyyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 173); 

                auto tg_xyz_yyyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 174); 

                auto tg_xyz_yyyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 175); 

                auto tg_xyz_yyyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 176); 

                auto tg_xyz_yyzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 177); 

                auto tg_xyz_yzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 178); 

                auto tg_xyz_zzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 179); 

                auto tg_xz_xxyyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 90); 

                auto tg_xz_xxyzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 91); 

                auto tg_xz_xxzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 92); 

                auto tg_xz_xyyyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 93); 

                auto tg_xz_xyyyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 94); 

                auto tg_xz_xyyyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 95); 

                auto tg_xz_xyyyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 96); 

                auto tg_xz_xyyzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 97); 

                auto tg_xz_xyzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 98); 

                auto tg_xz_xzzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 99); 

                auto tg_xz_yyyyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 100); 

                auto tg_xz_yyyyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 101); 

                auto tg_xz_yyyyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 102); 

                auto tg_xz_yyyyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 103); 

                auto tg_xz_yyyzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 104); 

                auto tg_xz_yyzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 105); 

                auto tg_xz_yzzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 106); 

                auto tg_xz_zzzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 107); 

                auto tg_yy_xxxxxxx_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 108); 

                auto tg_yy_xxxxxxy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 109); 

                auto tg_yy_xxxxxxz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 110); 

                auto tg_yy_xxxxxyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 111); 

                auto tg_yy_xxxxxyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 112); 

                auto tg_yy_xxxxxzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 113); 

                auto tg_yy_xxxxyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 114); 

                auto tg_yy_xxxxyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 115); 

                auto tg_yy_xxxxyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 116); 

                auto tg_yy_xxxxzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 117); 

                auto tg_yy_xxxyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 118); 

                auto tg_yy_xxxyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 119); 

                auto tg_yy_xxxyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 120); 

                auto tg_yy_xxxyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 121); 

                auto tg_yy_xxxzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 122); 

                auto tg_yy_xxyyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 123); 

                auto tg_yy_xxyyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 124); 

                auto tg_yy_xxyyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 125); 

                auto tg_yy_xxyyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 126); 

                auto tg_yy_xxyzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 127); 

                auto tg_yy_xxzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 128); 

                auto tg_yy_xyyyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 129); 

                auto tg_yy_xyyyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 130); 

                auto tg_yy_xyyyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 131); 

                auto tg_yy_xyyyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 132); 

                auto tg_yy_xyyzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 133); 

                auto tg_yy_xyzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 134); 

                auto tg_yy_xzzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 135); 

                auto tg_yy_yyyyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 136); 

                auto tg_yy_yyyyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 137); 

                auto tg_yy_yyyyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 138); 

                auto tg_yy_yyyyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 139); 

                auto tg_yy_yyyzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 140); 

                auto tg_yy_yyzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 141); 

                auto tg_yy_yzzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 142); 

                auto tg_yy_zzzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 143); 

                auto tg_yz_xxxxxxx_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 144); 

                auto tg_yz_xxxxxxy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 145); 

                auto tg_yz_xxxxxxz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 146); 

                auto tg_yz_xxxxxyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 147); 

                auto tg_yz_xxxxxyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 148); 

                auto tg_yz_xxxxxzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 149); 

                auto tg_yz_xxxxyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 150); 

                auto tg_yz_xxxxyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 151); 

                auto tg_yz_xxxxyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 152); 

                auto tg_yz_xxxxzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 153); 

                auto tg_yz_xxxyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 154); 

                auto tg_yz_xxxyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 155); 

                auto tg_yz_xxxyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 156); 

                auto tg_yz_xxxyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 157); 

                auto tg_yz_xxxzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 158); 

                auto tg_yz_xxyyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 159); 

                auto tg_yz_xxyyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 160); 

                auto tg_yz_xxyyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 161); 

                auto tg_yz_xxyyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 162); 

                auto tg_yz_xxyzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 163); 

                auto tg_yz_xxzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 164); 

                auto tg_yz_xyyyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 165); 

                auto tg_yz_xyyyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 166); 

                auto tg_yz_xyyyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 167); 

                auto tg_yz_xyyyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 168); 

                auto tg_yz_xyyzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 169); 

                auto tg_yz_xyzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 170); 

                auto tg_yz_xzzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 171); 

                auto tg_yz_yyyyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 172); 

                auto tg_yz_yyyyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 173); 

                auto tg_yz_yyyyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 174); 

                auto tg_yz_yyyyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 175); 

                auto tg_yz_yyyzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 176); 

                auto tg_yz_yyzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 177); 

                auto tg_yz_yzzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 178); 

                auto tg_yz_zzzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 179); 

                auto tg_xz_xxyyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 90); 

                auto tg_xz_xxyzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 91); 

                auto tg_xz_xxzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 92); 

                auto tg_xz_xyyyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 93); 

                auto tg_xz_xyyyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 94); 

                auto tg_xz_xyyyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 95); 

                auto tg_xz_xyyyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 96); 

                auto tg_xz_xyyzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 97); 

                auto tg_xz_xyzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 98); 

                auto tg_xz_xzzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 99); 

                auto tg_xz_yyyyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 100); 

                auto tg_xz_yyyyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 101); 

                auto tg_xz_yyyyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 102); 

                auto tg_xz_yyyyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 103); 

                auto tg_xz_yyyzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 104); 

                auto tg_xz_yyzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 105); 

                auto tg_xz_yzzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 106); 

                auto tg_xz_zzzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 107); 

                auto tg_yy_xxxxxxx_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 108); 

                auto tg_yy_xxxxxxy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 109); 

                auto tg_yy_xxxxxxz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 110); 

                auto tg_yy_xxxxxyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 111); 

                auto tg_yy_xxxxxyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 112); 

                auto tg_yy_xxxxxzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 113); 

                auto tg_yy_xxxxyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 114); 

                auto tg_yy_xxxxyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 115); 

                auto tg_yy_xxxxyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 116); 

                auto tg_yy_xxxxzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 117); 

                auto tg_yy_xxxyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 118); 

                auto tg_yy_xxxyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 119); 

                auto tg_yy_xxxyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 120); 

                auto tg_yy_xxxyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 121); 

                auto tg_yy_xxxzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 122); 

                auto tg_yy_xxyyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 123); 

                auto tg_yy_xxyyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 124); 

                auto tg_yy_xxyyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 125); 

                auto tg_yy_xxyyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 126); 

                auto tg_yy_xxyzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 127); 

                auto tg_yy_xxzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 128); 

                auto tg_yy_xyyyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 129); 

                auto tg_yy_xyyyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 130); 

                auto tg_yy_xyyyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 131); 

                auto tg_yy_xyyyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 132); 

                auto tg_yy_xyyzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 133); 

                auto tg_yy_xyzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 134); 

                auto tg_yy_xzzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 135); 

                auto tg_yy_yyyyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 136); 

                auto tg_yy_yyyyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 137); 

                auto tg_yy_yyyyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 138); 

                auto tg_yy_yyyyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 139); 

                auto tg_yy_yyyzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 140); 

                auto tg_yy_yyzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 141); 

                auto tg_yy_yzzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 142); 

                auto tg_yy_zzzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 143); 

                auto tg_yz_xxxxxxx_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 144); 

                auto tg_yz_xxxxxxy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 145); 

                auto tg_yz_xxxxxxz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 146); 

                auto tg_yz_xxxxxyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 147); 

                auto tg_yz_xxxxxyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 148); 

                auto tg_yz_xxxxxzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 149); 

                auto tg_yz_xxxxyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 150); 

                auto tg_yz_xxxxyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 151); 

                auto tg_yz_xxxxyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 152); 

                auto tg_yz_xxxxzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 153); 

                auto tg_yz_xxxyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 154); 

                auto tg_yz_xxxyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 155); 

                auto tg_yz_xxxyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 156); 

                auto tg_yz_xxxyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 157); 

                auto tg_yz_xxxzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 158); 

                auto tg_yz_xxyyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 159); 

                auto tg_yz_xxyyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 160); 

                auto tg_yz_xxyyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 161); 

                auto tg_yz_xxyyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 162); 

                auto tg_yz_xxyzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 163); 

                auto tg_yz_xxzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 164); 

                auto tg_yz_xyyyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 165); 

                auto tg_yz_xyyyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 166); 

                auto tg_yz_xyyyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 167); 

                auto tg_yz_xyyyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 168); 

                auto tg_yz_xyyzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 169); 

                auto tg_yz_xyzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 170); 

                auto tg_yz_xzzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 171); 

                auto tg_yz_yyyyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 172); 

                auto tg_yz_yyyyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 173); 

                auto tg_yz_yyyyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 174); 

                auto tg_yz_yyyyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 175); 

                auto tg_yz_yyyzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 176); 

                auto tg_yz_yyzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 177); 

                auto tg_yz_yzzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 178); 

                auto tg_yz_zzzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 179); 

                auto tg_xxz_xyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 74); 

                auto tg_xxz_xyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 75); 

                auto tg_xxz_xzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 76); 

                auto tg_xxz_yyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 77); 

                auto tg_xxz_yyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 78); 

                auto tg_xxz_yyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 79); 

                auto tg_xxz_yyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 80); 

                auto tg_xxz_yyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 81); 

                auto tg_xxz_yzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 82); 

                auto tg_xxz_zzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 83); 

                auto tg_xyy_xxxxxx_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 84); 

                auto tg_xyy_xxxxxy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 85); 

                auto tg_xyy_xxxxxz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 86); 

                auto tg_xyy_xxxxyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 87); 

                auto tg_xyy_xxxxyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 88); 

                auto tg_xyy_xxxxzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 89); 

                auto tg_xyy_xxxyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 90); 

                auto tg_xyy_xxxyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 91); 

                auto tg_xyy_xxxyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 92); 

                auto tg_xyy_xxxzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 93); 

                auto tg_xyy_xxyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 94); 

                auto tg_xyy_xxyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 95); 

                auto tg_xyy_xxyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 96); 

                auto tg_xyy_xxyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 97); 

                auto tg_xyy_xxzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 98); 

                auto tg_xyy_xyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 99); 

                auto tg_xyy_xyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 100); 

                auto tg_xyy_xyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 101); 

                auto tg_xyy_xyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 102); 

                auto tg_xyy_xyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 103); 

                auto tg_xyy_xzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 104); 

                auto tg_xyy_yyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 105); 

                auto tg_xyy_yyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 106); 

                auto tg_xyy_yyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 107); 

                auto tg_xyy_yyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 108); 

                auto tg_xyy_yyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 109); 

                auto tg_xyy_yzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 110); 

                auto tg_xyy_zzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 111); 

                auto tg_xyz_xxxxxx_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 112); 

                auto tg_xyz_xxxxxy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 113); 

                auto tg_xyz_xxxxxz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 114); 

                auto tg_xyz_xxxxyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 115); 

                auto tg_xyz_xxxxyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 116); 

                auto tg_xyz_xxxxzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 117); 

                auto tg_xyz_xxxyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 118); 

                auto tg_xyz_xxxyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 119); 

                auto tg_xyz_xxxyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 120); 

                auto tg_xyz_xxxzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 121); 

                auto tg_xyz_xxyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 122); 

                auto tg_xyz_xxyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 123); 

                auto tg_xyz_xxyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 124); 

                auto tg_xyz_xxyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 125); 

                auto tg_xyz_xxzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 126); 

                auto tg_xyz_xyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 127); 

                auto tg_xyz_xyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 128); 

                auto tg_xyz_xyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 129); 

                auto tg_xyz_xyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 130); 

                auto tg_xyz_xyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 131); 

                auto tg_xyz_xzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 132); 

                auto tg_xyz_yyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 133); 

                auto tg_xyz_yyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 134); 

                auto tg_xyz_yyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 135); 

                auto tg_xyz_yyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 136); 

                auto tg_xyz_yyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 137); 

                auto tg_xyz_yzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 138); 

                auto tg_xyz_zzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 139); 

                // set up pointers to integrals

                auto tg_xxxz_xxyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 90); 

                auto tg_xxxz_xxyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 91); 

                auto tg_xxxz_xxzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 92); 

                auto tg_xxxz_xyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 93); 

                auto tg_xxxz_xyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 94); 

                auto tg_xxxz_xyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 95); 

                auto tg_xxxz_xyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 96); 

                auto tg_xxxz_xyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 97); 

                auto tg_xxxz_xyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 98); 

                auto tg_xxxz_xzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 99); 

                auto tg_xxxz_yyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 100); 

                auto tg_xxxz_yyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 101); 

                auto tg_xxxz_yyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 102); 

                auto tg_xxxz_yyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 103); 

                auto tg_xxxz_yyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 104); 

                auto tg_xxxz_yyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 105); 

                auto tg_xxxz_yzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 106); 

                auto tg_xxxz_zzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 107); 

                auto tg_xxyy_xxxxxxx_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 108); 

                auto tg_xxyy_xxxxxxy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 109); 

                auto tg_xxyy_xxxxxxz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 110); 

                auto tg_xxyy_xxxxxyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 111); 

                auto tg_xxyy_xxxxxyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 112); 

                auto tg_xxyy_xxxxxzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 113); 

                auto tg_xxyy_xxxxyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 114); 

                auto tg_xxyy_xxxxyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 115); 

                auto tg_xxyy_xxxxyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 116); 

                auto tg_xxyy_xxxxzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 117); 

                auto tg_xxyy_xxxyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 118); 

                auto tg_xxyy_xxxyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 119); 

                auto tg_xxyy_xxxyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 120); 

                auto tg_xxyy_xxxyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 121); 

                auto tg_xxyy_xxxzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 122); 

                auto tg_xxyy_xxyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 123); 

                auto tg_xxyy_xxyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 124); 

                auto tg_xxyy_xxyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 125); 

                auto tg_xxyy_xxyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 126); 

                auto tg_xxyy_xxyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 127); 

                auto tg_xxyy_xxzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 128); 

                auto tg_xxyy_xyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 129); 

                auto tg_xxyy_xyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 130); 

                auto tg_xxyy_xyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 131); 

                auto tg_xxyy_xyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 132); 

                auto tg_xxyy_xyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 133); 

                auto tg_xxyy_xyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 134); 

                auto tg_xxyy_xzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 135); 

                auto tg_xxyy_yyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 136); 

                auto tg_xxyy_yyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 137); 

                auto tg_xxyy_yyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 138); 

                auto tg_xxyy_yyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 139); 

                auto tg_xxyy_yyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 140); 

                auto tg_xxyy_yyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 141); 

                auto tg_xxyy_yzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 142); 

                auto tg_xxyy_zzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 143); 

                auto tg_xxyz_xxxxxxx_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 144); 

                auto tg_xxyz_xxxxxxy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 145); 

                auto tg_xxyz_xxxxxxz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 146); 

                auto tg_xxyz_xxxxxyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 147); 

                auto tg_xxyz_xxxxxyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 148); 

                auto tg_xxyz_xxxxxzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 149); 

                auto tg_xxyz_xxxxyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 150); 

                auto tg_xxyz_xxxxyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 151); 

                auto tg_xxyz_xxxxyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 152); 

                auto tg_xxyz_xxxxzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 153); 

                auto tg_xxyz_xxxyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 154); 

                auto tg_xxyz_xxxyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 155); 

                auto tg_xxyz_xxxyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 156); 

                auto tg_xxyz_xxxyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 157); 

                auto tg_xxyz_xxxzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 158); 

                auto tg_xxyz_xxyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 159); 

                auto tg_xxyz_xxyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 160); 

                auto tg_xxyz_xxyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 161); 

                auto tg_xxyz_xxyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 162); 

                auto tg_xxyz_xxyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 163); 

                auto tg_xxyz_xxzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 164); 

                auto tg_xxyz_xyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 165); 

                auto tg_xxyz_xyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 166); 

                auto tg_xxyz_xyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 167); 

                auto tg_xxyz_xyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 168); 

                auto tg_xxyz_xyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 169); 

                auto tg_xxyz_xyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 170); 

                auto tg_xxyz_xzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 171); 

                auto tg_xxyz_yyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 172); 

                auto tg_xxyz_yyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 173); 

                auto tg_xxyz_yyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 174); 

                auto tg_xxyz_yyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 175); 

                auto tg_xxyz_yyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 176); 

                auto tg_xxyz_yyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 177); 

                auto tg_xxyz_yzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 178); 

                auto tg_xxyz_zzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 179); 

                // Batch of Integrals (90,180)

                #pragma omp simd aligned(fxn, fza, tg_xxxz_xxyyzzz_0, tg_xxxz_xxyzzzz_0, tg_xxxz_xxzzzzz_0, \
                                         tg_xxxz_xyyyyyy_0, tg_xxxz_xyyyyyz_0, tg_xxxz_xyyyyzz_0, tg_xxxz_xyyyzzz_0, \
                                         tg_xxxz_xyyzzzz_0, tg_xxxz_xyzzzzz_0, tg_xxxz_xzzzzzz_0, tg_xxxz_yyyyyyy_0, \
                                         tg_xxxz_yyyyyyz_0, tg_xxxz_yyyyyzz_0, tg_xxxz_yyyyzzz_0, tg_xxxz_yyyzzzz_0, \
                                         tg_xxxz_yyzzzzz_0, tg_xxxz_yzzzzzz_0, tg_xxxz_zzzzzzz_0, tg_xxyy_xxxxxxx_0, \
                                         tg_xxyy_xxxxxxy_0, tg_xxyy_xxxxxxz_0, tg_xxyy_xxxxxyy_0, tg_xxyy_xxxxxyz_0, \
                                         tg_xxyy_xxxxxzz_0, tg_xxyy_xxxxyyy_0, tg_xxyy_xxxxyyz_0, tg_xxyy_xxxxyzz_0, \
                                         tg_xxyy_xxxxzzz_0, tg_xxyy_xxxyyyy_0, tg_xxyy_xxxyyyz_0, tg_xxyy_xxxyyzz_0, \
                                         tg_xxyy_xxxyzzz_0, tg_xxyy_xxxzzzz_0, tg_xxyy_xxyyyyy_0, tg_xxyy_xxyyyyz_0, \
                                         tg_xxyy_xxyyyzz_0, tg_xxyy_xxyyzzz_0, tg_xxyy_xxyzzzz_0, tg_xxyy_xxzzzzz_0, \
                                         tg_xxyy_xyyyyyy_0, tg_xxyy_xyyyyyz_0, tg_xxyy_xyyyyzz_0, tg_xxyy_xyyyzzz_0, \
                                         tg_xxyy_xyyzzzz_0, tg_xxyy_xyzzzzz_0, tg_xxyy_xzzzzzz_0, tg_xxyy_yyyyyyy_0, \
                                         tg_xxyy_yyyyyyz_0, tg_xxyy_yyyyyzz_0, tg_xxyy_yyyyzzz_0, tg_xxyy_yyyzzzz_0, \
                                         tg_xxyy_yyzzzzz_0, tg_xxyy_yzzzzzz_0, tg_xxyy_zzzzzzz_0, tg_xxyz_xxxxxxx_0, \
                                         tg_xxyz_xxxxxxy_0, tg_xxyz_xxxxxxz_0, tg_xxyz_xxxxxyy_0, tg_xxyz_xxxxxyz_0, \
                                         tg_xxyz_xxxxxzz_0, tg_xxyz_xxxxyyy_0, tg_xxyz_xxxxyyz_0, tg_xxyz_xxxxyzz_0, \
                                         tg_xxyz_xxxxzzz_0, tg_xxyz_xxxyyyy_0, tg_xxyz_xxxyyyz_0, tg_xxyz_xxxyyzz_0, \
                                         tg_xxyz_xxxyzzz_0, tg_xxyz_xxxzzzz_0, tg_xxyz_xxyyyyy_0, tg_xxyz_xxyyyyz_0, \
                                         tg_xxyz_xxyyyzz_0, tg_xxyz_xxyyzzz_0, tg_xxyz_xxyzzzz_0, tg_xxyz_xxzzzzz_0, \
                                         tg_xxyz_xyyyyyy_0, tg_xxyz_xyyyyyz_0, tg_xxyz_xyyyyzz_0, tg_xxyz_xyyyzzz_0, \
                                         tg_xxyz_xyyzzzz_0, tg_xxyz_xyzzzzz_0, tg_xxyz_xzzzzzz_0, tg_xxyz_yyyyyyy_0, \
                                         tg_xxyz_yyyyyyz_0, tg_xxyz_yyyyyzz_0, tg_xxyz_yyyyzzz_0, tg_xxyz_yyyzzzz_0, \
                                         tg_xxyz_yyzzzzz_0, tg_xxyz_yzzzzzz_0, tg_xxyz_zzzzzzz_0, tg_xxz_xxyyzzz_0, \
                                         tg_xxz_xxyyzzz_1, tg_xxz_xxyzzzz_0, tg_xxz_xxyzzzz_1, tg_xxz_xxzzzzz_0, \
                                         tg_xxz_xxzzzzz_1, tg_xxz_xyyyyyy_0, tg_xxz_xyyyyyy_1, tg_xxz_xyyyyyz_0, \
                                         tg_xxz_xyyyyyz_1, tg_xxz_xyyyyzz_0, tg_xxz_xyyyyzz_1, tg_xxz_xyyyzzz_0, \
                                         tg_xxz_xyyyzzz_1, tg_xxz_xyyzzz_1, tg_xxz_xyyzzzz_0, tg_xxz_xyyzzzz_1, \
                                         tg_xxz_xyzzzz_1, tg_xxz_xyzzzzz_0, tg_xxz_xyzzzzz_1, tg_xxz_xzzzzz_1, \
                                         tg_xxz_xzzzzzz_0, tg_xxz_xzzzzzz_1, tg_xxz_yyyyyy_1, tg_xxz_yyyyyyy_0, \
                                         tg_xxz_yyyyyyy_1, tg_xxz_yyyyyyz_0, tg_xxz_yyyyyyz_1, tg_xxz_yyyyyz_1, \
                                         tg_xxz_yyyyyzz_0, tg_xxz_yyyyyzz_1, tg_xxz_yyyyzz_1, tg_xxz_yyyyzzz_0, \
                                         tg_xxz_yyyyzzz_1, tg_xxz_yyyzzz_1, tg_xxz_yyyzzzz_0, tg_xxz_yyyzzzz_1, \
                                         tg_xxz_yyzzzz_1, tg_xxz_yyzzzzz_0, tg_xxz_yyzzzzz_1, tg_xxz_yzzzzz_1, \
                                         tg_xxz_yzzzzzz_0, tg_xxz_yzzzzzz_1, tg_xxz_zzzzzz_1, tg_xxz_zzzzzzz_0, \
                                         tg_xxz_zzzzzzz_1, tg_xyy_xxxxxx_1, tg_xyy_xxxxxxx_0, tg_xyy_xxxxxxx_1, \
                                         tg_xyy_xxxxxxy_0, tg_xyy_xxxxxxy_1, tg_xyy_xxxxxxz_0, tg_xyy_xxxxxxz_1, \
                                         tg_xyy_xxxxxy_1, tg_xyy_xxxxxyy_0, tg_xyy_xxxxxyy_1, tg_xyy_xxxxxyz_0, \
                                         tg_xyy_xxxxxyz_1, tg_xyy_xxxxxz_1, tg_xyy_xxxxxzz_0, tg_xyy_xxxxxzz_1, \
                                         tg_xyy_xxxxyy_1, tg_xyy_xxxxyyy_0, tg_xyy_xxxxyyy_1, tg_xyy_xxxxyyz_0, \
                                         tg_xyy_xxxxyyz_1, tg_xyy_xxxxyz_1, tg_xyy_xxxxyzz_0, tg_xyy_xxxxyzz_1, \
                                         tg_xyy_xxxxzz_1, tg_xyy_xxxxzzz_0, tg_xyy_xxxxzzz_1, tg_xyy_xxxyyy_1, \
                                         tg_xyy_xxxyyyy_0, tg_xyy_xxxyyyy_1, tg_xyy_xxxyyyz_0, tg_xyy_xxxyyyz_1, \
                                         tg_xyy_xxxyyz_1, tg_xyy_xxxyyzz_0, tg_xyy_xxxyyzz_1, tg_xyy_xxxyzz_1, \
                                         tg_xyy_xxxyzzz_0, tg_xyy_xxxyzzz_1, tg_xyy_xxxzzz_1, tg_xyy_xxxzzzz_0, \
                                         tg_xyy_xxxzzzz_1, tg_xyy_xxyyyy_1, tg_xyy_xxyyyyy_0, tg_xyy_xxyyyyy_1, \
                                         tg_xyy_xxyyyyz_0, tg_xyy_xxyyyyz_1, tg_xyy_xxyyyz_1, tg_xyy_xxyyyzz_0, \
                                         tg_xyy_xxyyyzz_1, tg_xyy_xxyyzz_1, tg_xyy_xxyyzzz_0, tg_xyy_xxyyzzz_1, \
                                         tg_xyy_xxyzzz_1, tg_xyy_xxyzzzz_0, tg_xyy_xxyzzzz_1, tg_xyy_xxzzzz_1, \
                                         tg_xyy_xxzzzzz_0, tg_xyy_xxzzzzz_1, tg_xyy_xyyyyy_1, tg_xyy_xyyyyyy_0, \
                                         tg_xyy_xyyyyyy_1, tg_xyy_xyyyyyz_0, tg_xyy_xyyyyyz_1, tg_xyy_xyyyyz_1, \
                                         tg_xyy_xyyyyzz_0, tg_xyy_xyyyyzz_1, tg_xyy_xyyyzz_1, tg_xyy_xyyyzzz_0, \
                                         tg_xyy_xyyyzzz_1, tg_xyy_xyyzzz_1, tg_xyy_xyyzzzz_0, tg_xyy_xyyzzzz_1, \
                                         tg_xyy_xyzzzz_1, tg_xyy_xyzzzzz_0, tg_xyy_xyzzzzz_1, tg_xyy_xzzzzz_1, \
                                         tg_xyy_xzzzzzz_0, tg_xyy_xzzzzzz_1, tg_xyy_yyyyyy_1, tg_xyy_yyyyyyy_0, \
                                         tg_xyy_yyyyyyy_1, tg_xyy_yyyyyyz_0, tg_xyy_yyyyyyz_1, tg_xyy_yyyyyz_1, \
                                         tg_xyy_yyyyyzz_0, tg_xyy_yyyyyzz_1, tg_xyy_yyyyzz_1, tg_xyy_yyyyzzz_0, \
                                         tg_xyy_yyyyzzz_1, tg_xyy_yyyzzz_1, tg_xyy_yyyzzzz_0, tg_xyy_yyyzzzz_1, \
                                         tg_xyy_yyzzzz_1, tg_xyy_yyzzzzz_0, tg_xyy_yyzzzzz_1, tg_xyy_yzzzzz_1, \
                                         tg_xyy_yzzzzzz_0, tg_xyy_yzzzzzz_1, tg_xyy_zzzzzz_1, tg_xyy_zzzzzzz_0, \
                                         tg_xyy_zzzzzzz_1, tg_xyz_xxxxxx_1, tg_xyz_xxxxxxx_0, tg_xyz_xxxxxxx_1, \
                                         tg_xyz_xxxxxxy_0, tg_xyz_xxxxxxy_1, tg_xyz_xxxxxxz_0, tg_xyz_xxxxxxz_1, \
                                         tg_xyz_xxxxxy_1, tg_xyz_xxxxxyy_0, tg_xyz_xxxxxyy_1, tg_xyz_xxxxxyz_0, \
                                         tg_xyz_xxxxxyz_1, tg_xyz_xxxxxz_1, tg_xyz_xxxxxzz_0, tg_xyz_xxxxxzz_1, \
                                         tg_xyz_xxxxyy_1, tg_xyz_xxxxyyy_0, tg_xyz_xxxxyyy_1, tg_xyz_xxxxyyz_0, \
                                         tg_xyz_xxxxyyz_1, tg_xyz_xxxxyz_1, tg_xyz_xxxxyzz_0, tg_xyz_xxxxyzz_1, \
                                         tg_xyz_xxxxzz_1, tg_xyz_xxxxzzz_0, tg_xyz_xxxxzzz_1, tg_xyz_xxxyyy_1, \
                                         tg_xyz_xxxyyyy_0, tg_xyz_xxxyyyy_1, tg_xyz_xxxyyyz_0, tg_xyz_xxxyyyz_1, \
                                         tg_xyz_xxxyyz_1, tg_xyz_xxxyyzz_0, tg_xyz_xxxyyzz_1, tg_xyz_xxxyzz_1, \
                                         tg_xyz_xxxyzzz_0, tg_xyz_xxxyzzz_1, tg_xyz_xxxzzz_1, tg_xyz_xxxzzzz_0, \
                                         tg_xyz_xxxzzzz_1, tg_xyz_xxyyyy_1, tg_xyz_xxyyyyy_0, tg_xyz_xxyyyyy_1, \
                                         tg_xyz_xxyyyyz_0, tg_xyz_xxyyyyz_1, tg_xyz_xxyyyz_1, tg_xyz_xxyyyzz_0, \
                                         tg_xyz_xxyyyzz_1, tg_xyz_xxyyzz_1, tg_xyz_xxyyzzz_0, tg_xyz_xxyyzzz_1, \
                                         tg_xyz_xxyzzz_1, tg_xyz_xxyzzzz_0, tg_xyz_xxyzzzz_1, tg_xyz_xxzzzz_1, \
                                         tg_xyz_xxzzzzz_0, tg_xyz_xxzzzzz_1, tg_xyz_xyyyyy_1, tg_xyz_xyyyyyy_0, \
                                         tg_xyz_xyyyyyy_1, tg_xyz_xyyyyyz_0, tg_xyz_xyyyyyz_1, tg_xyz_xyyyyz_1, \
                                         tg_xyz_xyyyyzz_0, tg_xyz_xyyyyzz_1, tg_xyz_xyyyzz_1, tg_xyz_xyyyzzz_0, \
                                         tg_xyz_xyyyzzz_1, tg_xyz_xyyzzz_1, tg_xyz_xyyzzzz_0, tg_xyz_xyyzzzz_1, \
                                         tg_xyz_xyzzzz_1, tg_xyz_xyzzzzz_0, tg_xyz_xyzzzzz_1, tg_xyz_xzzzzz_1, \
                                         tg_xyz_xzzzzzz_0, tg_xyz_xzzzzzz_1, tg_xyz_yyyyyy_1, tg_xyz_yyyyyyy_0, \
                                         tg_xyz_yyyyyyy_1, tg_xyz_yyyyyyz_0, tg_xyz_yyyyyyz_1, tg_xyz_yyyyyz_1, \
                                         tg_xyz_yyyyyzz_0, tg_xyz_yyyyyzz_1, tg_xyz_yyyyzz_1, tg_xyz_yyyyzzz_0, \
                                         tg_xyz_yyyyzzz_1, tg_xyz_yyyzzz_1, tg_xyz_yyyzzzz_0, tg_xyz_yyyzzzz_1, \
                                         tg_xyz_yyzzzz_1, tg_xyz_yyzzzzz_0, tg_xyz_yyzzzzz_1, tg_xyz_yzzzzz_1, \
                                         tg_xyz_yzzzzzz_0, tg_xyz_yzzzzzz_1, tg_xyz_zzzzzz_1, tg_xyz_zzzzzzz_0, \
                                         tg_xyz_zzzzzzz_1, tg_xz_xxyyzzz_0, tg_xz_xxyyzzz_1, tg_xz_xxyzzzz_0, tg_xz_xxyzzzz_1, \
                                         tg_xz_xxzzzzz_0, tg_xz_xxzzzzz_1, tg_xz_xyyyyyy_0, tg_xz_xyyyyyy_1, tg_xz_xyyyyyz_0, \
                                         tg_xz_xyyyyyz_1, tg_xz_xyyyyzz_0, tg_xz_xyyyyzz_1, tg_xz_xyyyzzz_0, tg_xz_xyyyzzz_1, \
                                         tg_xz_xyyzzzz_0, tg_xz_xyyzzzz_1, tg_xz_xyzzzzz_0, tg_xz_xyzzzzz_1, tg_xz_xzzzzzz_0, \
                                         tg_xz_xzzzzzz_1, tg_xz_yyyyyyy_0, tg_xz_yyyyyyy_1, tg_xz_yyyyyyz_0, tg_xz_yyyyyyz_1, \
                                         tg_xz_yyyyyzz_0, tg_xz_yyyyyzz_1, tg_xz_yyyyzzz_0, tg_xz_yyyyzzz_1, tg_xz_yyyzzzz_0, \
                                         tg_xz_yyyzzzz_1, tg_xz_yyzzzzz_0, tg_xz_yyzzzzz_1, tg_xz_yzzzzzz_0, tg_xz_yzzzzzz_1, \
                                         tg_xz_zzzzzzz_0, tg_xz_zzzzzzz_1, tg_yy_xxxxxxx_0, tg_yy_xxxxxxx_1, tg_yy_xxxxxxy_0, \
                                         tg_yy_xxxxxxy_1, tg_yy_xxxxxxz_0, tg_yy_xxxxxxz_1, tg_yy_xxxxxyy_0, tg_yy_xxxxxyy_1, \
                                         tg_yy_xxxxxyz_0, tg_yy_xxxxxyz_1, tg_yy_xxxxxzz_0, tg_yy_xxxxxzz_1, tg_yy_xxxxyyy_0, \
                                         tg_yy_xxxxyyy_1, tg_yy_xxxxyyz_0, tg_yy_xxxxyyz_1, tg_yy_xxxxyzz_0, tg_yy_xxxxyzz_1, \
                                         tg_yy_xxxxzzz_0, tg_yy_xxxxzzz_1, tg_yy_xxxyyyy_0, tg_yy_xxxyyyy_1, tg_yy_xxxyyyz_0, \
                                         tg_yy_xxxyyyz_1, tg_yy_xxxyyzz_0, tg_yy_xxxyyzz_1, tg_yy_xxxyzzz_0, tg_yy_xxxyzzz_1, \
                                         tg_yy_xxxzzzz_0, tg_yy_xxxzzzz_1, tg_yy_xxyyyyy_0, tg_yy_xxyyyyy_1, tg_yy_xxyyyyz_0, \
                                         tg_yy_xxyyyyz_1, tg_yy_xxyyyzz_0, tg_yy_xxyyyzz_1, tg_yy_xxyyzzz_0, tg_yy_xxyyzzz_1, \
                                         tg_yy_xxyzzzz_0, tg_yy_xxyzzzz_1, tg_yy_xxzzzzz_0, tg_yy_xxzzzzz_1, tg_yy_xyyyyyy_0, \
                                         tg_yy_xyyyyyy_1, tg_yy_xyyyyyz_0, tg_yy_xyyyyyz_1, tg_yy_xyyyyzz_0, tg_yy_xyyyyzz_1, \
                                         tg_yy_xyyyzzz_0, tg_yy_xyyyzzz_1, tg_yy_xyyzzzz_0, tg_yy_xyyzzzz_1, tg_yy_xyzzzzz_0, \
                                         tg_yy_xyzzzzz_1, tg_yy_xzzzzzz_0, tg_yy_xzzzzzz_1, tg_yy_yyyyyyy_0, tg_yy_yyyyyyy_1, \
                                         tg_yy_yyyyyyz_0, tg_yy_yyyyyyz_1, tg_yy_yyyyyzz_0, tg_yy_yyyyyzz_1, tg_yy_yyyyzzz_0, \
                                         tg_yy_yyyyzzz_1, tg_yy_yyyzzzz_0, tg_yy_yyyzzzz_1, tg_yy_yyzzzzz_0, tg_yy_yyzzzzz_1, \
                                         tg_yy_yzzzzzz_0, tg_yy_yzzzzzz_1, tg_yy_zzzzzzz_0, tg_yy_zzzzzzz_1, tg_yz_xxxxxxx_0, \
                                         tg_yz_xxxxxxx_1, tg_yz_xxxxxxy_0, tg_yz_xxxxxxy_1, tg_yz_xxxxxxz_0, tg_yz_xxxxxxz_1, \
                                         tg_yz_xxxxxyy_0, tg_yz_xxxxxyy_1, tg_yz_xxxxxyz_0, tg_yz_xxxxxyz_1, tg_yz_xxxxxzz_0, \
                                         tg_yz_xxxxxzz_1, tg_yz_xxxxyyy_0, tg_yz_xxxxyyy_1, tg_yz_xxxxyyz_0, tg_yz_xxxxyyz_1, \
                                         tg_yz_xxxxyzz_0, tg_yz_xxxxyzz_1, tg_yz_xxxxzzz_0, tg_yz_xxxxzzz_1, tg_yz_xxxyyyy_0, \
                                         tg_yz_xxxyyyy_1, tg_yz_xxxyyyz_0, tg_yz_xxxyyyz_1, tg_yz_xxxyyzz_0, tg_yz_xxxyyzz_1, \
                                         tg_yz_xxxyzzz_0, tg_yz_xxxyzzz_1, tg_yz_xxxzzzz_0, tg_yz_xxxzzzz_1, tg_yz_xxyyyyy_0, \
                                         tg_yz_xxyyyyy_1, tg_yz_xxyyyyz_0, tg_yz_xxyyyyz_1, tg_yz_xxyyyzz_0, tg_yz_xxyyyzz_1, \
                                         tg_yz_xxyyzzz_0, tg_yz_xxyyzzz_1, tg_yz_xxyzzzz_0, tg_yz_xxyzzzz_1, tg_yz_xxzzzzz_0, \
                                         tg_yz_xxzzzzz_1, tg_yz_xyyyyyy_0, tg_yz_xyyyyyy_1, tg_yz_xyyyyyz_0, tg_yz_xyyyyyz_1, \
                                         tg_yz_xyyyyzz_0, tg_yz_xyyyyzz_1, tg_yz_xyyyzzz_0, tg_yz_xyyyzzz_1, tg_yz_xyyzzzz_0, \
                                         tg_yz_xyyzzzz_1, tg_yz_xyzzzzz_0, tg_yz_xyzzzzz_1, tg_yz_xzzzzzz_0, tg_yz_xzzzzzz_1, \
                                         tg_yz_yyyyyyy_0, tg_yz_yyyyyyy_1, tg_yz_yyyyyyz_0, tg_yz_yyyyyyz_1, tg_yz_yyyyyzz_0, \
                                         tg_yz_yyyyyzz_1, tg_yz_yyyyzzz_0, tg_yz_yyyyzzz_1, tg_yz_yyyzzzz_0, tg_yz_yyyzzzz_1, \
                                         tg_yz_yyzzzzz_0, tg_yz_yyzzzzz_1, tg_yz_yzzzzzz_0, tg_yz_yzzzzzz_1, tg_yz_zzzzzzz_0, \
                                         tg_yz_zzzzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxz_xxyyzzz_0[j] = pb_x * tg_xxz_xxyyzzz_0[j] + fr * tg_xxz_xxyyzzz_1[j] + fl1_fx * (tg_xz_xxyyzzz_0[j] - tg_xz_xxyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxz_xyyzzz_1[j];

                    tg_xxxz_xxyzzzz_0[j] = pb_x * tg_xxz_xxyzzzz_0[j] + fr * tg_xxz_xxyzzzz_1[j] + fl1_fx * (tg_xz_xxyzzzz_0[j] - tg_xz_xxyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxz_xyzzzz_1[j];

                    tg_xxxz_xxzzzzz_0[j] = pb_x * tg_xxz_xxzzzzz_0[j] + fr * tg_xxz_xxzzzzz_1[j] + fl1_fx * (tg_xz_xxzzzzz_0[j] - tg_xz_xxzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxz_xzzzzz_1[j];

                    tg_xxxz_xyyyyyy_0[j] = pb_x * tg_xxz_xyyyyyy_0[j] + fr * tg_xxz_xyyyyyy_1[j] + fl1_fx * (tg_xz_xyyyyyy_0[j] - tg_xz_xyyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxz_yyyyyy_1[j];

                    tg_xxxz_xyyyyyz_0[j] = pb_x * tg_xxz_xyyyyyz_0[j] + fr * tg_xxz_xyyyyyz_1[j] + fl1_fx * (tg_xz_xyyyyyz_0[j] - tg_xz_xyyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxz_yyyyyz_1[j];

                    tg_xxxz_xyyyyzz_0[j] = pb_x * tg_xxz_xyyyyzz_0[j] + fr * tg_xxz_xyyyyzz_1[j] + fl1_fx * (tg_xz_xyyyyzz_0[j] - tg_xz_xyyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxz_yyyyzz_1[j];

                    tg_xxxz_xyyyzzz_0[j] = pb_x * tg_xxz_xyyyzzz_0[j] + fr * tg_xxz_xyyyzzz_1[j] + fl1_fx * (tg_xz_xyyyzzz_0[j] - tg_xz_xyyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxz_yyyzzz_1[j];

                    tg_xxxz_xyyzzzz_0[j] = pb_x * tg_xxz_xyyzzzz_0[j] + fr * tg_xxz_xyyzzzz_1[j] + fl1_fx * (tg_xz_xyyzzzz_0[j] - tg_xz_xyyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxz_yyzzzz_1[j];

                    tg_xxxz_xyzzzzz_0[j] = pb_x * tg_xxz_xyzzzzz_0[j] + fr * tg_xxz_xyzzzzz_1[j] + fl1_fx * (tg_xz_xyzzzzz_0[j] - tg_xz_xyzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxz_yzzzzz_1[j];

                    tg_xxxz_xzzzzzz_0[j] = pb_x * tg_xxz_xzzzzzz_0[j] + fr * tg_xxz_xzzzzzz_1[j] + fl1_fx * (tg_xz_xzzzzzz_0[j] - tg_xz_xzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxz_zzzzzz_1[j];

                    tg_xxxz_yyyyyyy_0[j] = pb_x * tg_xxz_yyyyyyy_0[j] + fr * tg_xxz_yyyyyyy_1[j] + fl1_fx * (tg_xz_yyyyyyy_0[j] - tg_xz_yyyyyyy_1[j] * fl1_fza);

                    tg_xxxz_yyyyyyz_0[j] = pb_x * tg_xxz_yyyyyyz_0[j] + fr * tg_xxz_yyyyyyz_1[j] + fl1_fx * (tg_xz_yyyyyyz_0[j] - tg_xz_yyyyyyz_1[j] * fl1_fza);

                    tg_xxxz_yyyyyzz_0[j] = pb_x * tg_xxz_yyyyyzz_0[j] + fr * tg_xxz_yyyyyzz_1[j] + fl1_fx * (tg_xz_yyyyyzz_0[j] - tg_xz_yyyyyzz_1[j] * fl1_fza);

                    tg_xxxz_yyyyzzz_0[j] = pb_x * tg_xxz_yyyyzzz_0[j] + fr * tg_xxz_yyyyzzz_1[j] + fl1_fx * (tg_xz_yyyyzzz_0[j] - tg_xz_yyyyzzz_1[j] * fl1_fza);

                    tg_xxxz_yyyzzzz_0[j] = pb_x * tg_xxz_yyyzzzz_0[j] + fr * tg_xxz_yyyzzzz_1[j] + fl1_fx * (tg_xz_yyyzzzz_0[j] - tg_xz_yyyzzzz_1[j] * fl1_fza);

                    tg_xxxz_yyzzzzz_0[j] = pb_x * tg_xxz_yyzzzzz_0[j] + fr * tg_xxz_yyzzzzz_1[j] + fl1_fx * (tg_xz_yyzzzzz_0[j] - tg_xz_yyzzzzz_1[j] * fl1_fza);

                    tg_xxxz_yzzzzzz_0[j] = pb_x * tg_xxz_yzzzzzz_0[j] + fr * tg_xxz_yzzzzzz_1[j] + fl1_fx * (tg_xz_yzzzzzz_0[j] - tg_xz_yzzzzzz_1[j] * fl1_fza);

                    tg_xxxz_zzzzzzz_0[j] = pb_x * tg_xxz_zzzzzzz_0[j] + fr * tg_xxz_zzzzzzz_1[j] + fl1_fx * (tg_xz_zzzzzzz_0[j] - tg_xz_zzzzzzz_1[j] * fl1_fza);

                    tg_xxyy_xxxxxxx_0[j] = pb_x * tg_xyy_xxxxxxx_0[j] + fr * tg_xyy_xxxxxxx_1[j] + 0.5 * fl1_fx * (tg_yy_xxxxxxx_0[j] - tg_yy_xxxxxxx_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_xyy_xxxxxx_1[j];

                    tg_xxyy_xxxxxxy_0[j] = pb_x * tg_xyy_xxxxxxy_0[j] + fr * tg_xyy_xxxxxxy_1[j] + 0.5 * fl1_fx * (tg_yy_xxxxxxy_0[j] - tg_yy_xxxxxxy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xyy_xxxxxy_1[j];

                    tg_xxyy_xxxxxxz_0[j] = pb_x * tg_xyy_xxxxxxz_0[j] + fr * tg_xyy_xxxxxxz_1[j] + 0.5 * fl1_fx * (tg_yy_xxxxxxz_0[j] - tg_yy_xxxxxxz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xyy_xxxxxz_1[j];

                    tg_xxyy_xxxxxyy_0[j] = pb_x * tg_xyy_xxxxxyy_0[j] + fr * tg_xyy_xxxxxyy_1[j] + 0.5 * fl1_fx * (tg_yy_xxxxxyy_0[j] - tg_yy_xxxxxyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyy_xxxxyy_1[j];

                    tg_xxyy_xxxxxyz_0[j] = pb_x * tg_xyy_xxxxxyz_0[j] + fr * tg_xyy_xxxxxyz_1[j] + 0.5 * fl1_fx * (tg_yy_xxxxxyz_0[j] - tg_yy_xxxxxyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyy_xxxxyz_1[j];

                    tg_xxyy_xxxxxzz_0[j] = pb_x * tg_xyy_xxxxxzz_0[j] + fr * tg_xyy_xxxxxzz_1[j] + 0.5 * fl1_fx * (tg_yy_xxxxxzz_0[j] - tg_yy_xxxxxzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyy_xxxxzz_1[j];

                    tg_xxyy_xxxxyyy_0[j] = pb_x * tg_xyy_xxxxyyy_0[j] + fr * tg_xyy_xxxxyyy_1[j] + 0.5 * fl1_fx * (tg_yy_xxxxyyy_0[j] - tg_yy_xxxxyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyy_xxxyyy_1[j];

                    tg_xxyy_xxxxyyz_0[j] = pb_x * tg_xyy_xxxxyyz_0[j] + fr * tg_xyy_xxxxyyz_1[j] + 0.5 * fl1_fx * (tg_yy_xxxxyyz_0[j] - tg_yy_xxxxyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyy_xxxyyz_1[j];

                    tg_xxyy_xxxxyzz_0[j] = pb_x * tg_xyy_xxxxyzz_0[j] + fr * tg_xyy_xxxxyzz_1[j] + 0.5 * fl1_fx * (tg_yy_xxxxyzz_0[j] - tg_yy_xxxxyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyy_xxxyzz_1[j];

                    tg_xxyy_xxxxzzz_0[j] = pb_x * tg_xyy_xxxxzzz_0[j] + fr * tg_xyy_xxxxzzz_1[j] + 0.5 * fl1_fx * (tg_yy_xxxxzzz_0[j] - tg_yy_xxxxzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyy_xxxzzz_1[j];

                    tg_xxyy_xxxyyyy_0[j] = pb_x * tg_xyy_xxxyyyy_0[j] + fr * tg_xyy_xxxyyyy_1[j] + 0.5 * fl1_fx * (tg_yy_xxxyyyy_0[j] - tg_yy_xxxyyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyy_xxyyyy_1[j];

                    tg_xxyy_xxxyyyz_0[j] = pb_x * tg_xyy_xxxyyyz_0[j] + fr * tg_xyy_xxxyyyz_1[j] + 0.5 * fl1_fx * (tg_yy_xxxyyyz_0[j] - tg_yy_xxxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyy_xxyyyz_1[j];

                    tg_xxyy_xxxyyzz_0[j] = pb_x * tg_xyy_xxxyyzz_0[j] + fr * tg_xyy_xxxyyzz_1[j] + 0.5 * fl1_fx * (tg_yy_xxxyyzz_0[j] - tg_yy_xxxyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyy_xxyyzz_1[j];

                    tg_xxyy_xxxyzzz_0[j] = pb_x * tg_xyy_xxxyzzz_0[j] + fr * tg_xyy_xxxyzzz_1[j] + 0.5 * fl1_fx * (tg_yy_xxxyzzz_0[j] - tg_yy_xxxyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyy_xxyzzz_1[j];

                    tg_xxyy_xxxzzzz_0[j] = pb_x * tg_xyy_xxxzzzz_0[j] + fr * tg_xyy_xxxzzzz_1[j] + 0.5 * fl1_fx * (tg_yy_xxxzzzz_0[j] - tg_yy_xxxzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyy_xxzzzz_1[j];

                    tg_xxyy_xxyyyyy_0[j] = pb_x * tg_xyy_xxyyyyy_0[j] + fr * tg_xyy_xxyyyyy_1[j] + 0.5 * fl1_fx * (tg_yy_xxyyyyy_0[j] - tg_yy_xxyyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyy_xyyyyy_1[j];

                    tg_xxyy_xxyyyyz_0[j] = pb_x * tg_xyy_xxyyyyz_0[j] + fr * tg_xyy_xxyyyyz_1[j] + 0.5 * fl1_fx * (tg_yy_xxyyyyz_0[j] - tg_yy_xxyyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyy_xyyyyz_1[j];

                    tg_xxyy_xxyyyzz_0[j] = pb_x * tg_xyy_xxyyyzz_0[j] + fr * tg_xyy_xxyyyzz_1[j] + 0.5 * fl1_fx * (tg_yy_xxyyyzz_0[j] - tg_yy_xxyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyy_xyyyzz_1[j];

                    tg_xxyy_xxyyzzz_0[j] = pb_x * tg_xyy_xxyyzzz_0[j] + fr * tg_xyy_xxyyzzz_1[j] + 0.5 * fl1_fx * (tg_yy_xxyyzzz_0[j] - tg_yy_xxyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyy_xyyzzz_1[j];

                    tg_xxyy_xxyzzzz_0[j] = pb_x * tg_xyy_xxyzzzz_0[j] + fr * tg_xyy_xxyzzzz_1[j] + 0.5 * fl1_fx * (tg_yy_xxyzzzz_0[j] - tg_yy_xxyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyy_xyzzzz_1[j];

                    tg_xxyy_xxzzzzz_0[j] = pb_x * tg_xyy_xxzzzzz_0[j] + fr * tg_xyy_xxzzzzz_1[j] + 0.5 * fl1_fx * (tg_yy_xxzzzzz_0[j] - tg_yy_xxzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyy_xzzzzz_1[j];

                    tg_xxyy_xyyyyyy_0[j] = pb_x * tg_xyy_xyyyyyy_0[j] + fr * tg_xyy_xyyyyyy_1[j] + 0.5 * fl1_fx * (tg_yy_xyyyyyy_0[j] - tg_yy_xyyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyy_yyyyyy_1[j];

                    tg_xxyy_xyyyyyz_0[j] = pb_x * tg_xyy_xyyyyyz_0[j] + fr * tg_xyy_xyyyyyz_1[j] + 0.5 * fl1_fx * (tg_yy_xyyyyyz_0[j] - tg_yy_xyyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyy_yyyyyz_1[j];

                    tg_xxyy_xyyyyzz_0[j] = pb_x * tg_xyy_xyyyyzz_0[j] + fr * tg_xyy_xyyyyzz_1[j] + 0.5 * fl1_fx * (tg_yy_xyyyyzz_0[j] - tg_yy_xyyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyy_yyyyzz_1[j];

                    tg_xxyy_xyyyzzz_0[j] = pb_x * tg_xyy_xyyyzzz_0[j] + fr * tg_xyy_xyyyzzz_1[j] + 0.5 * fl1_fx * (tg_yy_xyyyzzz_0[j] - tg_yy_xyyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyy_yyyzzz_1[j];

                    tg_xxyy_xyyzzzz_0[j] = pb_x * tg_xyy_xyyzzzz_0[j] + fr * tg_xyy_xyyzzzz_1[j] + 0.5 * fl1_fx * (tg_yy_xyyzzzz_0[j] - tg_yy_xyyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyy_yyzzzz_1[j];

                    tg_xxyy_xyzzzzz_0[j] = pb_x * tg_xyy_xyzzzzz_0[j] + fr * tg_xyy_xyzzzzz_1[j] + 0.5 * fl1_fx * (tg_yy_xyzzzzz_0[j] - tg_yy_xyzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyy_yzzzzz_1[j];

                    tg_xxyy_xzzzzzz_0[j] = pb_x * tg_xyy_xzzzzzz_0[j] + fr * tg_xyy_xzzzzzz_1[j] + 0.5 * fl1_fx * (tg_yy_xzzzzzz_0[j] - tg_yy_xzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyy_zzzzzz_1[j];

                    tg_xxyy_yyyyyyy_0[j] = pb_x * tg_xyy_yyyyyyy_0[j] + fr * tg_xyy_yyyyyyy_1[j] + 0.5 * fl1_fx * (tg_yy_yyyyyyy_0[j] - tg_yy_yyyyyyy_1[j] * fl1_fza);

                    tg_xxyy_yyyyyyz_0[j] = pb_x * tg_xyy_yyyyyyz_0[j] + fr * tg_xyy_yyyyyyz_1[j] + 0.5 * fl1_fx * (tg_yy_yyyyyyz_0[j] - tg_yy_yyyyyyz_1[j] * fl1_fza);

                    tg_xxyy_yyyyyzz_0[j] = pb_x * tg_xyy_yyyyyzz_0[j] + fr * tg_xyy_yyyyyzz_1[j] + 0.5 * fl1_fx * (tg_yy_yyyyyzz_0[j] - tg_yy_yyyyyzz_1[j] * fl1_fza);

                    tg_xxyy_yyyyzzz_0[j] = pb_x * tg_xyy_yyyyzzz_0[j] + fr * tg_xyy_yyyyzzz_1[j] + 0.5 * fl1_fx * (tg_yy_yyyyzzz_0[j] - tg_yy_yyyyzzz_1[j] * fl1_fza);

                    tg_xxyy_yyyzzzz_0[j] = pb_x * tg_xyy_yyyzzzz_0[j] + fr * tg_xyy_yyyzzzz_1[j] + 0.5 * fl1_fx * (tg_yy_yyyzzzz_0[j] - tg_yy_yyyzzzz_1[j] * fl1_fza);

                    tg_xxyy_yyzzzzz_0[j] = pb_x * tg_xyy_yyzzzzz_0[j] + fr * tg_xyy_yyzzzzz_1[j] + 0.5 * fl1_fx * (tg_yy_yyzzzzz_0[j] - tg_yy_yyzzzzz_1[j] * fl1_fza);

                    tg_xxyy_yzzzzzz_0[j] = pb_x * tg_xyy_yzzzzzz_0[j] + fr * tg_xyy_yzzzzzz_1[j] + 0.5 * fl1_fx * (tg_yy_yzzzzzz_0[j] - tg_yy_yzzzzzz_1[j] * fl1_fza);

                    tg_xxyy_zzzzzzz_0[j] = pb_x * tg_xyy_zzzzzzz_0[j] + fr * tg_xyy_zzzzzzz_1[j] + 0.5 * fl1_fx * (tg_yy_zzzzzzz_0[j] - tg_yy_zzzzzzz_1[j] * fl1_fza);

                    tg_xxyz_xxxxxxx_0[j] = pb_x * tg_xyz_xxxxxxx_0[j] + fr * tg_xyz_xxxxxxx_1[j] + 0.5 * fl1_fx * (tg_yz_xxxxxxx_0[j] - tg_yz_xxxxxxx_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_xyz_xxxxxx_1[j];

                    tg_xxyz_xxxxxxy_0[j] = pb_x * tg_xyz_xxxxxxy_0[j] + fr * tg_xyz_xxxxxxy_1[j] + 0.5 * fl1_fx * (tg_yz_xxxxxxy_0[j] - tg_yz_xxxxxxy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xyz_xxxxxy_1[j];

                    tg_xxyz_xxxxxxz_0[j] = pb_x * tg_xyz_xxxxxxz_0[j] + fr * tg_xyz_xxxxxxz_1[j] + 0.5 * fl1_fx * (tg_yz_xxxxxxz_0[j] - tg_yz_xxxxxxz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xyz_xxxxxz_1[j];

                    tg_xxyz_xxxxxyy_0[j] = pb_x * tg_xyz_xxxxxyy_0[j] + fr * tg_xyz_xxxxxyy_1[j] + 0.5 * fl1_fx * (tg_yz_xxxxxyy_0[j] - tg_yz_xxxxxyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyz_xxxxyy_1[j];

                    tg_xxyz_xxxxxyz_0[j] = pb_x * tg_xyz_xxxxxyz_0[j] + fr * tg_xyz_xxxxxyz_1[j] + 0.5 * fl1_fx * (tg_yz_xxxxxyz_0[j] - tg_yz_xxxxxyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyz_xxxxyz_1[j];

                    tg_xxyz_xxxxxzz_0[j] = pb_x * tg_xyz_xxxxxzz_0[j] + fr * tg_xyz_xxxxxzz_1[j] + 0.5 * fl1_fx * (tg_yz_xxxxxzz_0[j] - tg_yz_xxxxxzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyz_xxxxzz_1[j];

                    tg_xxyz_xxxxyyy_0[j] = pb_x * tg_xyz_xxxxyyy_0[j] + fr * tg_xyz_xxxxyyy_1[j] + 0.5 * fl1_fx * (tg_yz_xxxxyyy_0[j] - tg_yz_xxxxyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyz_xxxyyy_1[j];

                    tg_xxyz_xxxxyyz_0[j] = pb_x * tg_xyz_xxxxyyz_0[j] + fr * tg_xyz_xxxxyyz_1[j] + 0.5 * fl1_fx * (tg_yz_xxxxyyz_0[j] - tg_yz_xxxxyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyz_xxxyyz_1[j];

                    tg_xxyz_xxxxyzz_0[j] = pb_x * tg_xyz_xxxxyzz_0[j] + fr * tg_xyz_xxxxyzz_1[j] + 0.5 * fl1_fx * (tg_yz_xxxxyzz_0[j] - tg_yz_xxxxyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyz_xxxyzz_1[j];

                    tg_xxyz_xxxxzzz_0[j] = pb_x * tg_xyz_xxxxzzz_0[j] + fr * tg_xyz_xxxxzzz_1[j] + 0.5 * fl1_fx * (tg_yz_xxxxzzz_0[j] - tg_yz_xxxxzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyz_xxxzzz_1[j];

                    tg_xxyz_xxxyyyy_0[j] = pb_x * tg_xyz_xxxyyyy_0[j] + fr * tg_xyz_xxxyyyy_1[j] + 0.5 * fl1_fx * (tg_yz_xxxyyyy_0[j] - tg_yz_xxxyyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyz_xxyyyy_1[j];

                    tg_xxyz_xxxyyyz_0[j] = pb_x * tg_xyz_xxxyyyz_0[j] + fr * tg_xyz_xxxyyyz_1[j] + 0.5 * fl1_fx * (tg_yz_xxxyyyz_0[j] - tg_yz_xxxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyz_xxyyyz_1[j];

                    tg_xxyz_xxxyyzz_0[j] = pb_x * tg_xyz_xxxyyzz_0[j] + fr * tg_xyz_xxxyyzz_1[j] + 0.5 * fl1_fx * (tg_yz_xxxyyzz_0[j] - tg_yz_xxxyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyz_xxyyzz_1[j];

                    tg_xxyz_xxxyzzz_0[j] = pb_x * tg_xyz_xxxyzzz_0[j] + fr * tg_xyz_xxxyzzz_1[j] + 0.5 * fl1_fx * (tg_yz_xxxyzzz_0[j] - tg_yz_xxxyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyz_xxyzzz_1[j];

                    tg_xxyz_xxxzzzz_0[j] = pb_x * tg_xyz_xxxzzzz_0[j] + fr * tg_xyz_xxxzzzz_1[j] + 0.5 * fl1_fx * (tg_yz_xxxzzzz_0[j] - tg_yz_xxxzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyz_xxzzzz_1[j];

                    tg_xxyz_xxyyyyy_0[j] = pb_x * tg_xyz_xxyyyyy_0[j] + fr * tg_xyz_xxyyyyy_1[j] + 0.5 * fl1_fx * (tg_yz_xxyyyyy_0[j] - tg_yz_xxyyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyz_xyyyyy_1[j];

                    tg_xxyz_xxyyyyz_0[j] = pb_x * tg_xyz_xxyyyyz_0[j] + fr * tg_xyz_xxyyyyz_1[j] + 0.5 * fl1_fx * (tg_yz_xxyyyyz_0[j] - tg_yz_xxyyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyz_xyyyyz_1[j];

                    tg_xxyz_xxyyyzz_0[j] = pb_x * tg_xyz_xxyyyzz_0[j] + fr * tg_xyz_xxyyyzz_1[j] + 0.5 * fl1_fx * (tg_yz_xxyyyzz_0[j] - tg_yz_xxyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyz_xyyyzz_1[j];

                    tg_xxyz_xxyyzzz_0[j] = pb_x * tg_xyz_xxyyzzz_0[j] + fr * tg_xyz_xxyyzzz_1[j] + 0.5 * fl1_fx * (tg_yz_xxyyzzz_0[j] - tg_yz_xxyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyz_xyyzzz_1[j];

                    tg_xxyz_xxyzzzz_0[j] = pb_x * tg_xyz_xxyzzzz_0[j] + fr * tg_xyz_xxyzzzz_1[j] + 0.5 * fl1_fx * (tg_yz_xxyzzzz_0[j] - tg_yz_xxyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyz_xyzzzz_1[j];

                    tg_xxyz_xxzzzzz_0[j] = pb_x * tg_xyz_xxzzzzz_0[j] + fr * tg_xyz_xxzzzzz_1[j] + 0.5 * fl1_fx * (tg_yz_xxzzzzz_0[j] - tg_yz_xxzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyz_xzzzzz_1[j];

                    tg_xxyz_xyyyyyy_0[j] = pb_x * tg_xyz_xyyyyyy_0[j] + fr * tg_xyz_xyyyyyy_1[j] + 0.5 * fl1_fx * (tg_yz_xyyyyyy_0[j] - tg_yz_xyyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyz_yyyyyy_1[j];

                    tg_xxyz_xyyyyyz_0[j] = pb_x * tg_xyz_xyyyyyz_0[j] + fr * tg_xyz_xyyyyyz_1[j] + 0.5 * fl1_fx * (tg_yz_xyyyyyz_0[j] - tg_yz_xyyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyz_yyyyyz_1[j];

                    tg_xxyz_xyyyyzz_0[j] = pb_x * tg_xyz_xyyyyzz_0[j] + fr * tg_xyz_xyyyyzz_1[j] + 0.5 * fl1_fx * (tg_yz_xyyyyzz_0[j] - tg_yz_xyyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyz_yyyyzz_1[j];

                    tg_xxyz_xyyyzzz_0[j] = pb_x * tg_xyz_xyyyzzz_0[j] + fr * tg_xyz_xyyyzzz_1[j] + 0.5 * fl1_fx * (tg_yz_xyyyzzz_0[j] - tg_yz_xyyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyz_yyyzzz_1[j];

                    tg_xxyz_xyyzzzz_0[j] = pb_x * tg_xyz_xyyzzzz_0[j] + fr * tg_xyz_xyyzzzz_1[j] + 0.5 * fl1_fx * (tg_yz_xyyzzzz_0[j] - tg_yz_xyyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyz_yyzzzz_1[j];

                    tg_xxyz_xyzzzzz_0[j] = pb_x * tg_xyz_xyzzzzz_0[j] + fr * tg_xyz_xyzzzzz_1[j] + 0.5 * fl1_fx * (tg_yz_xyzzzzz_0[j] - tg_yz_xyzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyz_yzzzzz_1[j];

                    tg_xxyz_xzzzzzz_0[j] = pb_x * tg_xyz_xzzzzzz_0[j] + fr * tg_xyz_xzzzzzz_1[j] + 0.5 * fl1_fx * (tg_yz_xzzzzzz_0[j] - tg_yz_xzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyz_zzzzzz_1[j];

                    tg_xxyz_yyyyyyy_0[j] = pb_x * tg_xyz_yyyyyyy_0[j] + fr * tg_xyz_yyyyyyy_1[j] + 0.5 * fl1_fx * (tg_yz_yyyyyyy_0[j] - tg_yz_yyyyyyy_1[j] * fl1_fza);

                    tg_xxyz_yyyyyyz_0[j] = pb_x * tg_xyz_yyyyyyz_0[j] + fr * tg_xyz_yyyyyyz_1[j] + 0.5 * fl1_fx * (tg_yz_yyyyyyz_0[j] - tg_yz_yyyyyyz_1[j] * fl1_fza);

                    tg_xxyz_yyyyyzz_0[j] = pb_x * tg_xyz_yyyyyzz_0[j] + fr * tg_xyz_yyyyyzz_1[j] + 0.5 * fl1_fx * (tg_yz_yyyyyzz_0[j] - tg_yz_yyyyyzz_1[j] * fl1_fza);

                    tg_xxyz_yyyyzzz_0[j] = pb_x * tg_xyz_yyyyzzz_0[j] + fr * tg_xyz_yyyyzzz_1[j] + 0.5 * fl1_fx * (tg_yz_yyyyzzz_0[j] - tg_yz_yyyyzzz_1[j] * fl1_fza);

                    tg_xxyz_yyyzzzz_0[j] = pb_x * tg_xyz_yyyzzzz_0[j] + fr * tg_xyz_yyyzzzz_1[j] + 0.5 * fl1_fx * (tg_yz_yyyzzzz_0[j] - tg_yz_yyyzzzz_1[j] * fl1_fza);

                    tg_xxyz_yyzzzzz_0[j] = pb_x * tg_xyz_yyzzzzz_0[j] + fr * tg_xyz_yyzzzzz_1[j] + 0.5 * fl1_fx * (tg_yz_yyzzzzz_0[j] - tg_yz_yyzzzzz_1[j] * fl1_fza);

                    tg_xxyz_yzzzzzz_0[j] = pb_x * tg_xyz_yzzzzzz_0[j] + fr * tg_xyz_yzzzzzz_1[j] + 0.5 * fl1_fx * (tg_yz_yzzzzzz_0[j] - tg_yz_yzzzzzz_1[j] * fl1_fza);

                    tg_xxyz_zzzzzzz_0[j] = pb_x * tg_xyz_zzzzzzz_0[j] + fr * tg_xyz_zzzzzzz_1[j] + 0.5 * fl1_fx * (tg_yz_zzzzzzz_0[j] - tg_yz_zzzzzzz_1[j] * fl1_fza);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSK_180_270(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (180,270)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {4, -1, -1, -1},
                                             {7, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_7_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {7, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_7_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_7_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_7_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_2_7_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {7, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_2_7_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {7, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {6, -1, -1, -1}, 
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

                auto tg_xzz_xxxxxxx_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 180); 

                auto tg_xzz_xxxxxxy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 181); 

                auto tg_xzz_xxxxxxz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 182); 

                auto tg_xzz_xxxxxyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 183); 

                auto tg_xzz_xxxxxyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 184); 

                auto tg_xzz_xxxxxzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 185); 

                auto tg_xzz_xxxxyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 186); 

                auto tg_xzz_xxxxyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 187); 

                auto tg_xzz_xxxxyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 188); 

                auto tg_xzz_xxxxzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 189); 

                auto tg_xzz_xxxyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 190); 

                auto tg_xzz_xxxyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 191); 

                auto tg_xzz_xxxyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 192); 

                auto tg_xzz_xxxyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 193); 

                auto tg_xzz_xxxzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 194); 

                auto tg_xzz_xxyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 195); 

                auto tg_xzz_xxyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 196); 

                auto tg_xzz_xxyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 197); 

                auto tg_xzz_xxyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 198); 

                auto tg_xzz_xxyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 199); 

                auto tg_xzz_xxzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 200); 

                auto tg_xzz_xyyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 201); 

                auto tg_xzz_xyyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 202); 

                auto tg_xzz_xyyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 203); 

                auto tg_xzz_xyyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 204); 

                auto tg_xzz_xyyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 205); 

                auto tg_xzz_xyzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 206); 

                auto tg_xzz_xzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 207); 

                auto tg_xzz_yyyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 208); 

                auto tg_xzz_yyyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 209); 

                auto tg_xzz_yyyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 210); 

                auto tg_xzz_yyyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 211); 

                auto tg_xzz_yyyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 212); 

                auto tg_xzz_yyzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 213); 

                auto tg_xzz_yzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 214); 

                auto tg_xzz_zzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 215); 

                auto tg_yyy_xxxxxxx_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 216); 

                auto tg_yyy_xxxxxxy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 217); 

                auto tg_yyy_xxxxxxz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 218); 

                auto tg_yyy_xxxxxyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 219); 

                auto tg_yyy_xxxxxyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 220); 

                auto tg_yyy_xxxxxzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 221); 

                auto tg_yyy_xxxxyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 222); 

                auto tg_yyy_xxxxyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 223); 

                auto tg_yyy_xxxxyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 224); 

                auto tg_yyy_xxxxzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 225); 

                auto tg_yyy_xxxyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 226); 

                auto tg_yyy_xxxyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 227); 

                auto tg_yyy_xxxyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 228); 

                auto tg_yyy_xxxyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 229); 

                auto tg_yyy_xxxzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 230); 

                auto tg_yyy_xxyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 231); 

                auto tg_yyy_xxyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 232); 

                auto tg_yyy_xxyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 233); 

                auto tg_yyy_xxyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 234); 

                auto tg_yyy_xxyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 235); 

                auto tg_yyy_xxzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 236); 

                auto tg_yyy_xyyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 237); 

                auto tg_yyy_xyyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 238); 

                auto tg_yyy_xyyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 239); 

                auto tg_yyy_xyyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 240); 

                auto tg_yyy_xyyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 241); 

                auto tg_yyy_xyzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 242); 

                auto tg_yyy_xzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 243); 

                auto tg_yyy_yyyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 244); 

                auto tg_yyy_yyyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 245); 

                auto tg_yyy_yyyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 246); 

                auto tg_yyy_yyyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 247); 

                auto tg_yyy_yyyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 248); 

                auto tg_yyy_yyzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 249); 

                auto tg_yyy_yzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 250); 

                auto tg_yyy_zzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 251); 

                auto tg_yyz_xxxxxxx_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 252); 

                auto tg_yyz_xxxxxxy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 253); 

                auto tg_yyz_xxxxxxz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 254); 

                auto tg_yyz_xxxxxyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 255); 

                auto tg_yyz_xxxxxyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 256); 

                auto tg_yyz_xxxxxzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 257); 

                auto tg_yyz_xxxxyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 258); 

                auto tg_yyz_xxxxyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 259); 

                auto tg_yyz_xxxxyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 260); 

                auto tg_yyz_xxxxzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 261); 

                auto tg_yyz_xxxyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 262); 

                auto tg_yyz_xxxyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 263); 

                auto tg_yyz_xxxyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 264); 

                auto tg_yyz_xxxyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 265); 

                auto tg_yyz_xxxzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 266); 

                auto tg_yyz_xxyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 267); 

                auto tg_yyz_xxyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 268); 

                auto tg_yyz_xxyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 269); 

                auto tg_xzz_xxxxxxx_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 180); 

                auto tg_xzz_xxxxxxy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 181); 

                auto tg_xzz_xxxxxxz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 182); 

                auto tg_xzz_xxxxxyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 183); 

                auto tg_xzz_xxxxxyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 184); 

                auto tg_xzz_xxxxxzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 185); 

                auto tg_xzz_xxxxyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 186); 

                auto tg_xzz_xxxxyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 187); 

                auto tg_xzz_xxxxyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 188); 

                auto tg_xzz_xxxxzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 189); 

                auto tg_xzz_xxxyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 190); 

                auto tg_xzz_xxxyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 191); 

                auto tg_xzz_xxxyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 192); 

                auto tg_xzz_xxxyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 193); 

                auto tg_xzz_xxxzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 194); 

                auto tg_xzz_xxyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 195); 

                auto tg_xzz_xxyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 196); 

                auto tg_xzz_xxyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 197); 

                auto tg_xzz_xxyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 198); 

                auto tg_xzz_xxyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 199); 

                auto tg_xzz_xxzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 200); 

                auto tg_xzz_xyyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 201); 

                auto tg_xzz_xyyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 202); 

                auto tg_xzz_xyyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 203); 

                auto tg_xzz_xyyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 204); 

                auto tg_xzz_xyyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 205); 

                auto tg_xzz_xyzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 206); 

                auto tg_xzz_xzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 207); 

                auto tg_xzz_yyyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 208); 

                auto tg_xzz_yyyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 209); 

                auto tg_xzz_yyyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 210); 

                auto tg_xzz_yyyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 211); 

                auto tg_xzz_yyyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 212); 

                auto tg_xzz_yyzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 213); 

                auto tg_xzz_yzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 214); 

                auto tg_xzz_zzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 215); 

                auto tg_yyy_xxxxxxx_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 216); 

                auto tg_yyy_xxxxxxy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 217); 

                auto tg_yyy_xxxxxxz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 218); 

                auto tg_yyy_xxxxxyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 219); 

                auto tg_yyy_xxxxxyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 220); 

                auto tg_yyy_xxxxxzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 221); 

                auto tg_yyy_xxxxyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 222); 

                auto tg_yyy_xxxxyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 223); 

                auto tg_yyy_xxxxyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 224); 

                auto tg_yyy_xxxxzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 225); 

                auto tg_yyy_xxxyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 226); 

                auto tg_yyy_xxxyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 227); 

                auto tg_yyy_xxxyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 228); 

                auto tg_yyy_xxxyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 229); 

                auto tg_yyy_xxxzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 230); 

                auto tg_yyy_xxyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 231); 

                auto tg_yyy_xxyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 232); 

                auto tg_yyy_xxyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 233); 

                auto tg_yyy_xxyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 234); 

                auto tg_yyy_xxyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 235); 

                auto tg_yyy_xxzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 236); 

                auto tg_yyy_xyyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 237); 

                auto tg_yyy_xyyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 238); 

                auto tg_yyy_xyyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 239); 

                auto tg_yyy_xyyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 240); 

                auto tg_yyy_xyyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 241); 

                auto tg_yyy_xyzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 242); 

                auto tg_yyy_xzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 243); 

                auto tg_yyy_yyyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 244); 

                auto tg_yyy_yyyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 245); 

                auto tg_yyy_yyyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 246); 

                auto tg_yyy_yyyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 247); 

                auto tg_yyy_yyyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 248); 

                auto tg_yyy_yyzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 249); 

                auto tg_yyy_yzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 250); 

                auto tg_yyy_zzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 251); 

                auto tg_yyz_xxxxxxx_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 252); 

                auto tg_yyz_xxxxxxy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 253); 

                auto tg_yyz_xxxxxxz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 254); 

                auto tg_yyz_xxxxxyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 255); 

                auto tg_yyz_xxxxxyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 256); 

                auto tg_yyz_xxxxxzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 257); 

                auto tg_yyz_xxxxyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 258); 

                auto tg_yyz_xxxxyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 259); 

                auto tg_yyz_xxxxyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 260); 

                auto tg_yyz_xxxxzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 261); 

                auto tg_yyz_xxxyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 262); 

                auto tg_yyz_xxxyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 263); 

                auto tg_yyz_xxxyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 264); 

                auto tg_yyz_xxxyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 265); 

                auto tg_yyz_xxxzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 266); 

                auto tg_yyz_xxyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 267); 

                auto tg_yyz_xxyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 268); 

                auto tg_yyz_xxyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 269); 

                auto tg_zz_xxxxxxx_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 180); 

                auto tg_zz_xxxxxxy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 181); 

                auto tg_zz_xxxxxxz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 182); 

                auto tg_zz_xxxxxyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 183); 

                auto tg_zz_xxxxxyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 184); 

                auto tg_zz_xxxxxzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 185); 

                auto tg_zz_xxxxyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 186); 

                auto tg_zz_xxxxyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 187); 

                auto tg_zz_xxxxyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 188); 

                auto tg_zz_xxxxzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 189); 

                auto tg_zz_xxxyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 190); 

                auto tg_zz_xxxyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 191); 

                auto tg_zz_xxxyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 192); 

                auto tg_zz_xxxyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 193); 

                auto tg_zz_xxxzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 194); 

                auto tg_zz_xxyyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 195); 

                auto tg_zz_xxyyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 196); 

                auto tg_zz_xxyyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 197); 

                auto tg_zz_xxyyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 198); 

                auto tg_zz_xxyzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 199); 

                auto tg_zz_xxzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 200); 

                auto tg_zz_xyyyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 201); 

                auto tg_zz_xyyyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 202); 

                auto tg_zz_xyyyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 203); 

                auto tg_zz_xyyyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 204); 

                auto tg_zz_xyyzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 205); 

                auto tg_zz_xyzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 206); 

                auto tg_zz_xzzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 207); 

                auto tg_zz_yyyyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 208); 

                auto tg_zz_yyyyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 209); 

                auto tg_zz_yyyyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 210); 

                auto tg_zz_yyyyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 211); 

                auto tg_zz_yyyzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 212); 

                auto tg_zz_yyzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 213); 

                auto tg_zz_yzzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 214); 

                auto tg_zz_zzzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 215); 

                auto tg_zz_xxxxxxx_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 180); 

                auto tg_zz_xxxxxxy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 181); 

                auto tg_zz_xxxxxxz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 182); 

                auto tg_zz_xxxxxyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 183); 

                auto tg_zz_xxxxxyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 184); 

                auto tg_zz_xxxxxzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 185); 

                auto tg_zz_xxxxyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 186); 

                auto tg_zz_xxxxyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 187); 

                auto tg_zz_xxxxyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 188); 

                auto tg_zz_xxxxzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 189); 

                auto tg_zz_xxxyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 190); 

                auto tg_zz_xxxyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 191); 

                auto tg_zz_xxxyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 192); 

                auto tg_zz_xxxyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 193); 

                auto tg_zz_xxxzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 194); 

                auto tg_zz_xxyyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 195); 

                auto tg_zz_xxyyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 196); 

                auto tg_zz_xxyyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 197); 

                auto tg_zz_xxyyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 198); 

                auto tg_zz_xxyzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 199); 

                auto tg_zz_xxzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 200); 

                auto tg_zz_xyyyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 201); 

                auto tg_zz_xyyyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 202); 

                auto tg_zz_xyyyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 203); 

                auto tg_zz_xyyyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 204); 

                auto tg_zz_xyyzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 205); 

                auto tg_zz_xyzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 206); 

                auto tg_zz_xzzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 207); 

                auto tg_zz_yyyyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 208); 

                auto tg_zz_yyyyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 209); 

                auto tg_zz_yyyyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 210); 

                auto tg_zz_yyyyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 211); 

                auto tg_zz_yyyzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 212); 

                auto tg_zz_yyzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 213); 

                auto tg_zz_yzzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 214); 

                auto tg_zz_zzzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 215); 

                auto tg_xzz_xxxxxx_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 140); 

                auto tg_xzz_xxxxxy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 141); 

                auto tg_xzz_xxxxxz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 142); 

                auto tg_xzz_xxxxyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 143); 

                auto tg_xzz_xxxxyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 144); 

                auto tg_xzz_xxxxzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 145); 

                auto tg_xzz_xxxyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 146); 

                auto tg_xzz_xxxyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 147); 

                auto tg_xzz_xxxyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 148); 

                auto tg_xzz_xxxzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 149); 

                auto tg_xzz_xxyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 150); 

                auto tg_xzz_xxyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 151); 

                auto tg_xzz_xxyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 152); 

                auto tg_xzz_xxyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 153); 

                auto tg_xzz_xxzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 154); 

                auto tg_xzz_xyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 155); 

                auto tg_xzz_xyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 156); 

                auto tg_xzz_xyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 157); 

                auto tg_xzz_xyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 158); 

                auto tg_xzz_xyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 159); 

                auto tg_xzz_xzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 160); 

                auto tg_xzz_yyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 161); 

                auto tg_xzz_yyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 162); 

                auto tg_xzz_yyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 163); 

                auto tg_xzz_yyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 164); 

                auto tg_xzz_yyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 165); 

                auto tg_xzz_yzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 166); 

                auto tg_xzz_zzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 167); 

                auto tg_yyy_xxxxxx_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 168); 

                auto tg_yyy_xxxxxy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 169); 

                auto tg_yyy_xxxxxz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 170); 

                auto tg_yyy_xxxxyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 171); 

                auto tg_yyy_xxxxyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 172); 

                auto tg_yyy_xxxxzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 173); 

                auto tg_yyy_xxxyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 174); 

                auto tg_yyy_xxxyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 175); 

                auto tg_yyy_xxxyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 176); 

                auto tg_yyy_xxxzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 177); 

                auto tg_yyy_xxyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 178); 

                auto tg_yyy_xxyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 179); 

                auto tg_yyy_xxyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 180); 

                auto tg_yyy_xxyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 181); 

                auto tg_yyy_xxzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 182); 

                auto tg_yyy_xyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 183); 

                auto tg_yyy_xyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 184); 

                auto tg_yyy_xyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 185); 

                auto tg_yyy_xyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 186); 

                auto tg_yyy_xyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 187); 

                auto tg_yyy_xzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 188); 

                auto tg_yyy_yyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 189); 

                auto tg_yyy_yyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 190); 

                auto tg_yyy_yyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 191); 

                auto tg_yyy_yyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 192); 

                auto tg_yyy_yyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 193); 

                auto tg_yyy_yzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 194); 

                auto tg_yyy_zzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 195); 

                auto tg_yyz_xxxxxx_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 196); 

                auto tg_yyz_xxxxxy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 197); 

                auto tg_yyz_xxxxxz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 198); 

                auto tg_yyz_xxxxyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 199); 

                auto tg_yyz_xxxxyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 200); 

                auto tg_yyz_xxxxzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 201); 

                auto tg_yyz_xxxyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 202); 

                auto tg_yyz_xxxyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 203); 

                auto tg_yyz_xxxyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 204); 

                auto tg_yyz_xxxzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 205); 

                auto tg_yyz_xxyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 206); 

                auto tg_yyz_xxyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 207); 

                auto tg_yyz_xxyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 208); 

                auto tg_yyz_xxyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 209); 

                auto tg_yyz_xxzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 210); 

                auto tg_yyz_xyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 211); 

                auto tg_yyz_xyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 212); 

                auto tg_yyz_xyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 213); 

                // set up pointers to integrals

                auto tg_xxzz_xxxxxxx_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 180); 

                auto tg_xxzz_xxxxxxy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 181); 

                auto tg_xxzz_xxxxxxz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 182); 

                auto tg_xxzz_xxxxxyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 183); 

                auto tg_xxzz_xxxxxyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 184); 

                auto tg_xxzz_xxxxxzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 185); 

                auto tg_xxzz_xxxxyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 186); 

                auto tg_xxzz_xxxxyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 187); 

                auto tg_xxzz_xxxxyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 188); 

                auto tg_xxzz_xxxxzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 189); 

                auto tg_xxzz_xxxyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 190); 

                auto tg_xxzz_xxxyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 191); 

                auto tg_xxzz_xxxyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 192); 

                auto tg_xxzz_xxxyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 193); 

                auto tg_xxzz_xxxzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 194); 

                auto tg_xxzz_xxyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 195); 

                auto tg_xxzz_xxyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 196); 

                auto tg_xxzz_xxyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 197); 

                auto tg_xxzz_xxyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 198); 

                auto tg_xxzz_xxyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 199); 

                auto tg_xxzz_xxzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 200); 

                auto tg_xxzz_xyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 201); 

                auto tg_xxzz_xyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 202); 

                auto tg_xxzz_xyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 203); 

                auto tg_xxzz_xyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 204); 

                auto tg_xxzz_xyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 205); 

                auto tg_xxzz_xyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 206); 

                auto tg_xxzz_xzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 207); 

                auto tg_xxzz_yyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 208); 

                auto tg_xxzz_yyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 209); 

                auto tg_xxzz_yyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 210); 

                auto tg_xxzz_yyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 211); 

                auto tg_xxzz_yyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 212); 

                auto tg_xxzz_yyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 213); 

                auto tg_xxzz_yzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 214); 

                auto tg_xxzz_zzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 215); 

                auto tg_xyyy_xxxxxxx_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 216); 

                auto tg_xyyy_xxxxxxy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 217); 

                auto tg_xyyy_xxxxxxz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 218); 

                auto tg_xyyy_xxxxxyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 219); 

                auto tg_xyyy_xxxxxyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 220); 

                auto tg_xyyy_xxxxxzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 221); 

                auto tg_xyyy_xxxxyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 222); 

                auto tg_xyyy_xxxxyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 223); 

                auto tg_xyyy_xxxxyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 224); 

                auto tg_xyyy_xxxxzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 225); 

                auto tg_xyyy_xxxyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 226); 

                auto tg_xyyy_xxxyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 227); 

                auto tg_xyyy_xxxyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 228); 

                auto tg_xyyy_xxxyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 229); 

                auto tg_xyyy_xxxzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 230); 

                auto tg_xyyy_xxyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 231); 

                auto tg_xyyy_xxyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 232); 

                auto tg_xyyy_xxyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 233); 

                auto tg_xyyy_xxyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 234); 

                auto tg_xyyy_xxyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 235); 

                auto tg_xyyy_xxzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 236); 

                auto tg_xyyy_xyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 237); 

                auto tg_xyyy_xyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 238); 

                auto tg_xyyy_xyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 239); 

                auto tg_xyyy_xyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 240); 

                auto tg_xyyy_xyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 241); 

                auto tg_xyyy_xyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 242); 

                auto tg_xyyy_xzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 243); 

                auto tg_xyyy_yyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 244); 

                auto tg_xyyy_yyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 245); 

                auto tg_xyyy_yyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 246); 

                auto tg_xyyy_yyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 247); 

                auto tg_xyyy_yyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 248); 

                auto tg_xyyy_yyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 249); 

                auto tg_xyyy_yzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 250); 

                auto tg_xyyy_zzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 251); 

                auto tg_xyyz_xxxxxxx_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 252); 

                auto tg_xyyz_xxxxxxy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 253); 

                auto tg_xyyz_xxxxxxz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 254); 

                auto tg_xyyz_xxxxxyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 255); 

                auto tg_xyyz_xxxxxyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 256); 

                auto tg_xyyz_xxxxxzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 257); 

                auto tg_xyyz_xxxxyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 258); 

                auto tg_xyyz_xxxxyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 259); 

                auto tg_xyyz_xxxxyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 260); 

                auto tg_xyyz_xxxxzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 261); 

                auto tg_xyyz_xxxyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 262); 

                auto tg_xyyz_xxxyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 263); 

                auto tg_xyyz_xxxyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 264); 

                auto tg_xyyz_xxxyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 265); 

                auto tg_xyyz_xxxzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 266); 

                auto tg_xyyz_xxyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 267); 

                auto tg_xyyz_xxyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 268); 

                auto tg_xyyz_xxyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 269); 

                // Batch of Integrals (180,270)

                #pragma omp simd aligned(fxn, fza, tg_xxzz_xxxxxxx_0, tg_xxzz_xxxxxxy_0, tg_xxzz_xxxxxxz_0, \
                                         tg_xxzz_xxxxxyy_0, tg_xxzz_xxxxxyz_0, tg_xxzz_xxxxxzz_0, tg_xxzz_xxxxyyy_0, \
                                         tg_xxzz_xxxxyyz_0, tg_xxzz_xxxxyzz_0, tg_xxzz_xxxxzzz_0, tg_xxzz_xxxyyyy_0, \
                                         tg_xxzz_xxxyyyz_0, tg_xxzz_xxxyyzz_0, tg_xxzz_xxxyzzz_0, tg_xxzz_xxxzzzz_0, \
                                         tg_xxzz_xxyyyyy_0, tg_xxzz_xxyyyyz_0, tg_xxzz_xxyyyzz_0, tg_xxzz_xxyyzzz_0, \
                                         tg_xxzz_xxyzzzz_0, tg_xxzz_xxzzzzz_0, tg_xxzz_xyyyyyy_0, tg_xxzz_xyyyyyz_0, \
                                         tg_xxzz_xyyyyzz_0, tg_xxzz_xyyyzzz_0, tg_xxzz_xyyzzzz_0, tg_xxzz_xyzzzzz_0, \
                                         tg_xxzz_xzzzzzz_0, tg_xxzz_yyyyyyy_0, tg_xxzz_yyyyyyz_0, tg_xxzz_yyyyyzz_0, \
                                         tg_xxzz_yyyyzzz_0, tg_xxzz_yyyzzzz_0, tg_xxzz_yyzzzzz_0, tg_xxzz_yzzzzzz_0, \
                                         tg_xxzz_zzzzzzz_0, tg_xyyy_xxxxxxx_0, tg_xyyy_xxxxxxy_0, tg_xyyy_xxxxxxz_0, \
                                         tg_xyyy_xxxxxyy_0, tg_xyyy_xxxxxyz_0, tg_xyyy_xxxxxzz_0, tg_xyyy_xxxxyyy_0, \
                                         tg_xyyy_xxxxyyz_0, tg_xyyy_xxxxyzz_0, tg_xyyy_xxxxzzz_0, tg_xyyy_xxxyyyy_0, \
                                         tg_xyyy_xxxyyyz_0, tg_xyyy_xxxyyzz_0, tg_xyyy_xxxyzzz_0, tg_xyyy_xxxzzzz_0, \
                                         tg_xyyy_xxyyyyy_0, tg_xyyy_xxyyyyz_0, tg_xyyy_xxyyyzz_0, tg_xyyy_xxyyzzz_0, \
                                         tg_xyyy_xxyzzzz_0, tg_xyyy_xxzzzzz_0, tg_xyyy_xyyyyyy_0, tg_xyyy_xyyyyyz_0, \
                                         tg_xyyy_xyyyyzz_0, tg_xyyy_xyyyzzz_0, tg_xyyy_xyyzzzz_0, tg_xyyy_xyzzzzz_0, \
                                         tg_xyyy_xzzzzzz_0, tg_xyyy_yyyyyyy_0, tg_xyyy_yyyyyyz_0, tg_xyyy_yyyyyzz_0, \
                                         tg_xyyy_yyyyzzz_0, tg_xyyy_yyyzzzz_0, tg_xyyy_yyzzzzz_0, tg_xyyy_yzzzzzz_0, \
                                         tg_xyyy_zzzzzzz_0, tg_xyyz_xxxxxxx_0, tg_xyyz_xxxxxxy_0, tg_xyyz_xxxxxxz_0, \
                                         tg_xyyz_xxxxxyy_0, tg_xyyz_xxxxxyz_0, tg_xyyz_xxxxxzz_0, tg_xyyz_xxxxyyy_0, \
                                         tg_xyyz_xxxxyyz_0, tg_xyyz_xxxxyzz_0, tg_xyyz_xxxxzzz_0, tg_xyyz_xxxyyyy_0, \
                                         tg_xyyz_xxxyyyz_0, tg_xyyz_xxxyyzz_0, tg_xyyz_xxxyzzz_0, tg_xyyz_xxxzzzz_0, \
                                         tg_xyyz_xxyyyyy_0, tg_xyyz_xxyyyyz_0, tg_xyyz_xxyyyzz_0, tg_xzz_xxxxxx_1, \
                                         tg_xzz_xxxxxxx_0, tg_xzz_xxxxxxx_1, tg_xzz_xxxxxxy_0, tg_xzz_xxxxxxy_1, \
                                         tg_xzz_xxxxxxz_0, tg_xzz_xxxxxxz_1, tg_xzz_xxxxxy_1, tg_xzz_xxxxxyy_0, \
                                         tg_xzz_xxxxxyy_1, tg_xzz_xxxxxyz_0, tg_xzz_xxxxxyz_1, tg_xzz_xxxxxz_1, \
                                         tg_xzz_xxxxxzz_0, tg_xzz_xxxxxzz_1, tg_xzz_xxxxyy_1, tg_xzz_xxxxyyy_0, \
                                         tg_xzz_xxxxyyy_1, tg_xzz_xxxxyyz_0, tg_xzz_xxxxyyz_1, tg_xzz_xxxxyz_1, \
                                         tg_xzz_xxxxyzz_0, tg_xzz_xxxxyzz_1, tg_xzz_xxxxzz_1, tg_xzz_xxxxzzz_0, \
                                         tg_xzz_xxxxzzz_1, tg_xzz_xxxyyy_1, tg_xzz_xxxyyyy_0, tg_xzz_xxxyyyy_1, \
                                         tg_xzz_xxxyyyz_0, tg_xzz_xxxyyyz_1, tg_xzz_xxxyyz_1, tg_xzz_xxxyyzz_0, \
                                         tg_xzz_xxxyyzz_1, tg_xzz_xxxyzz_1, tg_xzz_xxxyzzz_0, tg_xzz_xxxyzzz_1, \
                                         tg_xzz_xxxzzz_1, tg_xzz_xxxzzzz_0, tg_xzz_xxxzzzz_1, tg_xzz_xxyyyy_1, \
                                         tg_xzz_xxyyyyy_0, tg_xzz_xxyyyyy_1, tg_xzz_xxyyyyz_0, tg_xzz_xxyyyyz_1, \
                                         tg_xzz_xxyyyz_1, tg_xzz_xxyyyzz_0, tg_xzz_xxyyyzz_1, tg_xzz_xxyyzz_1, \
                                         tg_xzz_xxyyzzz_0, tg_xzz_xxyyzzz_1, tg_xzz_xxyzzz_1, tg_xzz_xxyzzzz_0, \
                                         tg_xzz_xxyzzzz_1, tg_xzz_xxzzzz_1, tg_xzz_xxzzzzz_0, tg_xzz_xxzzzzz_1, \
                                         tg_xzz_xyyyyy_1, tg_xzz_xyyyyyy_0, tg_xzz_xyyyyyy_1, tg_xzz_xyyyyyz_0, \
                                         tg_xzz_xyyyyyz_1, tg_xzz_xyyyyz_1, tg_xzz_xyyyyzz_0, tg_xzz_xyyyyzz_1, \
                                         tg_xzz_xyyyzz_1, tg_xzz_xyyyzzz_0, tg_xzz_xyyyzzz_1, tg_xzz_xyyzzz_1, \
                                         tg_xzz_xyyzzzz_0, tg_xzz_xyyzzzz_1, tg_xzz_xyzzzz_1, tg_xzz_xyzzzzz_0, \
                                         tg_xzz_xyzzzzz_1, tg_xzz_xzzzzz_1, tg_xzz_xzzzzzz_0, tg_xzz_xzzzzzz_1, \
                                         tg_xzz_yyyyyy_1, tg_xzz_yyyyyyy_0, tg_xzz_yyyyyyy_1, tg_xzz_yyyyyyz_0, \
                                         tg_xzz_yyyyyyz_1, tg_xzz_yyyyyz_1, tg_xzz_yyyyyzz_0, tg_xzz_yyyyyzz_1, \
                                         tg_xzz_yyyyzz_1, tg_xzz_yyyyzzz_0, tg_xzz_yyyyzzz_1, tg_xzz_yyyzzz_1, \
                                         tg_xzz_yyyzzzz_0, tg_xzz_yyyzzzz_1, tg_xzz_yyzzzz_1, tg_xzz_yyzzzzz_0, \
                                         tg_xzz_yyzzzzz_1, tg_xzz_yzzzzz_1, tg_xzz_yzzzzzz_0, tg_xzz_yzzzzzz_1, \
                                         tg_xzz_zzzzzz_1, tg_xzz_zzzzzzz_0, tg_xzz_zzzzzzz_1, tg_yyy_xxxxxx_1, \
                                         tg_yyy_xxxxxxx_0, tg_yyy_xxxxxxx_1, tg_yyy_xxxxxxy_0, tg_yyy_xxxxxxy_1, \
                                         tg_yyy_xxxxxxz_0, tg_yyy_xxxxxxz_1, tg_yyy_xxxxxy_1, tg_yyy_xxxxxyy_0, \
                                         tg_yyy_xxxxxyy_1, tg_yyy_xxxxxyz_0, tg_yyy_xxxxxyz_1, tg_yyy_xxxxxz_1, \
                                         tg_yyy_xxxxxzz_0, tg_yyy_xxxxxzz_1, tg_yyy_xxxxyy_1, tg_yyy_xxxxyyy_0, \
                                         tg_yyy_xxxxyyy_1, tg_yyy_xxxxyyz_0, tg_yyy_xxxxyyz_1, tg_yyy_xxxxyz_1, \
                                         tg_yyy_xxxxyzz_0, tg_yyy_xxxxyzz_1, tg_yyy_xxxxzz_1, tg_yyy_xxxxzzz_0, \
                                         tg_yyy_xxxxzzz_1, tg_yyy_xxxyyy_1, tg_yyy_xxxyyyy_0, tg_yyy_xxxyyyy_1, \
                                         tg_yyy_xxxyyyz_0, tg_yyy_xxxyyyz_1, tg_yyy_xxxyyz_1, tg_yyy_xxxyyzz_0, \
                                         tg_yyy_xxxyyzz_1, tg_yyy_xxxyzz_1, tg_yyy_xxxyzzz_0, tg_yyy_xxxyzzz_1, \
                                         tg_yyy_xxxzzz_1, tg_yyy_xxxzzzz_0, tg_yyy_xxxzzzz_1, tg_yyy_xxyyyy_1, \
                                         tg_yyy_xxyyyyy_0, tg_yyy_xxyyyyy_1, tg_yyy_xxyyyyz_0, tg_yyy_xxyyyyz_1, \
                                         tg_yyy_xxyyyz_1, tg_yyy_xxyyyzz_0, tg_yyy_xxyyyzz_1, tg_yyy_xxyyzz_1, \
                                         tg_yyy_xxyyzzz_0, tg_yyy_xxyyzzz_1, tg_yyy_xxyzzz_1, tg_yyy_xxyzzzz_0, \
                                         tg_yyy_xxyzzzz_1, tg_yyy_xxzzzz_1, tg_yyy_xxzzzzz_0, tg_yyy_xxzzzzz_1, \
                                         tg_yyy_xyyyyy_1, tg_yyy_xyyyyyy_0, tg_yyy_xyyyyyy_1, tg_yyy_xyyyyyz_0, \
                                         tg_yyy_xyyyyyz_1, tg_yyy_xyyyyz_1, tg_yyy_xyyyyzz_0, tg_yyy_xyyyyzz_1, \
                                         tg_yyy_xyyyzz_1, tg_yyy_xyyyzzz_0, tg_yyy_xyyyzzz_1, tg_yyy_xyyzzz_1, \
                                         tg_yyy_xyyzzzz_0, tg_yyy_xyyzzzz_1, tg_yyy_xyzzzz_1, tg_yyy_xyzzzzz_0, \
                                         tg_yyy_xyzzzzz_1, tg_yyy_xzzzzz_1, tg_yyy_xzzzzzz_0, tg_yyy_xzzzzzz_1, \
                                         tg_yyy_yyyyyy_1, tg_yyy_yyyyyyy_0, tg_yyy_yyyyyyy_1, tg_yyy_yyyyyyz_0, \
                                         tg_yyy_yyyyyyz_1, tg_yyy_yyyyyz_1, tg_yyy_yyyyyzz_0, tg_yyy_yyyyyzz_1, \
                                         tg_yyy_yyyyzz_1, tg_yyy_yyyyzzz_0, tg_yyy_yyyyzzz_1, tg_yyy_yyyzzz_1, \
                                         tg_yyy_yyyzzzz_0, tg_yyy_yyyzzzz_1, tg_yyy_yyzzzz_1, tg_yyy_yyzzzzz_0, \
                                         tg_yyy_yyzzzzz_1, tg_yyy_yzzzzz_1, tg_yyy_yzzzzzz_0, tg_yyy_yzzzzzz_1, \
                                         tg_yyy_zzzzzz_1, tg_yyy_zzzzzzz_0, tg_yyy_zzzzzzz_1, tg_yyz_xxxxxx_1, \
                                         tg_yyz_xxxxxxx_0, tg_yyz_xxxxxxx_1, tg_yyz_xxxxxxy_0, tg_yyz_xxxxxxy_1, \
                                         tg_yyz_xxxxxxz_0, tg_yyz_xxxxxxz_1, tg_yyz_xxxxxy_1, tg_yyz_xxxxxyy_0, \
                                         tg_yyz_xxxxxyy_1, tg_yyz_xxxxxyz_0, tg_yyz_xxxxxyz_1, tg_yyz_xxxxxz_1, \
                                         tg_yyz_xxxxxzz_0, tg_yyz_xxxxxzz_1, tg_yyz_xxxxyy_1, tg_yyz_xxxxyyy_0, \
                                         tg_yyz_xxxxyyy_1, tg_yyz_xxxxyyz_0, tg_yyz_xxxxyyz_1, tg_yyz_xxxxyz_1, \
                                         tg_yyz_xxxxyzz_0, tg_yyz_xxxxyzz_1, tg_yyz_xxxxzz_1, tg_yyz_xxxxzzz_0, \
                                         tg_yyz_xxxxzzz_1, tg_yyz_xxxyyy_1, tg_yyz_xxxyyyy_0, tg_yyz_xxxyyyy_1, \
                                         tg_yyz_xxxyyyz_0, tg_yyz_xxxyyyz_1, tg_yyz_xxxyyz_1, tg_yyz_xxxyyzz_0, \
                                         tg_yyz_xxxyyzz_1, tg_yyz_xxxyzz_1, tg_yyz_xxxyzzz_0, tg_yyz_xxxyzzz_1, \
                                         tg_yyz_xxxzzz_1, tg_yyz_xxxzzzz_0, tg_yyz_xxxzzzz_1, tg_yyz_xxyyyy_1, \
                                         tg_yyz_xxyyyyy_0, tg_yyz_xxyyyyy_1, tg_yyz_xxyyyyz_0, tg_yyz_xxyyyyz_1, \
                                         tg_yyz_xxyyyz_1, tg_yyz_xxyyyzz_0, tg_yyz_xxyyyzz_1, tg_yyz_xxyyzz_1, \
                                         tg_yyz_xxyzzz_1, tg_yyz_xxzzzz_1, tg_yyz_xyyyyy_1, tg_yyz_xyyyyz_1, tg_yyz_xyyyzz_1, \
                                         tg_zz_xxxxxxx_0, tg_zz_xxxxxxx_1, tg_zz_xxxxxxy_0, tg_zz_xxxxxxy_1, tg_zz_xxxxxxz_0, \
                                         tg_zz_xxxxxxz_1, tg_zz_xxxxxyy_0, tg_zz_xxxxxyy_1, tg_zz_xxxxxyz_0, tg_zz_xxxxxyz_1, \
                                         tg_zz_xxxxxzz_0, tg_zz_xxxxxzz_1, tg_zz_xxxxyyy_0, tg_zz_xxxxyyy_1, tg_zz_xxxxyyz_0, \
                                         tg_zz_xxxxyyz_1, tg_zz_xxxxyzz_0, tg_zz_xxxxyzz_1, tg_zz_xxxxzzz_0, tg_zz_xxxxzzz_1, \
                                         tg_zz_xxxyyyy_0, tg_zz_xxxyyyy_1, tg_zz_xxxyyyz_0, tg_zz_xxxyyyz_1, tg_zz_xxxyyzz_0, \
                                         tg_zz_xxxyyzz_1, tg_zz_xxxyzzz_0, tg_zz_xxxyzzz_1, tg_zz_xxxzzzz_0, tg_zz_xxxzzzz_1, \
                                         tg_zz_xxyyyyy_0, tg_zz_xxyyyyy_1, tg_zz_xxyyyyz_0, tg_zz_xxyyyyz_1, tg_zz_xxyyyzz_0, \
                                         tg_zz_xxyyyzz_1, tg_zz_xxyyzzz_0, tg_zz_xxyyzzz_1, tg_zz_xxyzzzz_0, tg_zz_xxyzzzz_1, \
                                         tg_zz_xxzzzzz_0, tg_zz_xxzzzzz_1, tg_zz_xyyyyyy_0, tg_zz_xyyyyyy_1, tg_zz_xyyyyyz_0, \
                                         tg_zz_xyyyyyz_1, tg_zz_xyyyyzz_0, tg_zz_xyyyyzz_1, tg_zz_xyyyzzz_0, tg_zz_xyyyzzz_1, \
                                         tg_zz_xyyzzzz_0, tg_zz_xyyzzzz_1, tg_zz_xyzzzzz_0, tg_zz_xyzzzzz_1, tg_zz_xzzzzzz_0, \
                                         tg_zz_xzzzzzz_1, tg_zz_yyyyyyy_0, tg_zz_yyyyyyy_1, tg_zz_yyyyyyz_0, tg_zz_yyyyyyz_1, \
                                         tg_zz_yyyyyzz_0, tg_zz_yyyyyzz_1, tg_zz_yyyyzzz_0, tg_zz_yyyyzzz_1, tg_zz_yyyzzzz_0, \
                                         tg_zz_yyyzzzz_1, tg_zz_yyzzzzz_0, tg_zz_yyzzzzz_1, tg_zz_yzzzzzz_0, tg_zz_yzzzzzz_1, \
                                         tg_zz_zzzzzzz_0, tg_zz_zzzzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxzz_xxxxxxx_0[j] = pb_x * tg_xzz_xxxxxxx_0[j] + fr * tg_xzz_xxxxxxx_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxxxx_0[j] - tg_zz_xxxxxxx_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_xzz_xxxxxx_1[j];

                    tg_xxzz_xxxxxxy_0[j] = pb_x * tg_xzz_xxxxxxy_0[j] + fr * tg_xzz_xxxxxxy_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxxxy_0[j] - tg_zz_xxxxxxy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xzz_xxxxxy_1[j];

                    tg_xxzz_xxxxxxz_0[j] = pb_x * tg_xzz_xxxxxxz_0[j] + fr * tg_xzz_xxxxxxz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxxxz_0[j] - tg_zz_xxxxxxz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xzz_xxxxxz_1[j];

                    tg_xxzz_xxxxxyy_0[j] = pb_x * tg_xzz_xxxxxyy_0[j] + fr * tg_xzz_xxxxxyy_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxxyy_0[j] - tg_zz_xxxxxyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xzz_xxxxyy_1[j];

                    tg_xxzz_xxxxxyz_0[j] = pb_x * tg_xzz_xxxxxyz_0[j] + fr * tg_xzz_xxxxxyz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxxyz_0[j] - tg_zz_xxxxxyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xzz_xxxxyz_1[j];

                    tg_xxzz_xxxxxzz_0[j] = pb_x * tg_xzz_xxxxxzz_0[j] + fr * tg_xzz_xxxxxzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxxzz_0[j] - tg_zz_xxxxxzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xzz_xxxxzz_1[j];

                    tg_xxzz_xxxxyyy_0[j] = pb_x * tg_xzz_xxxxyyy_0[j] + fr * tg_xzz_xxxxyyy_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxyyy_0[j] - tg_zz_xxxxyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzz_xxxyyy_1[j];

                    tg_xxzz_xxxxyyz_0[j] = pb_x * tg_xzz_xxxxyyz_0[j] + fr * tg_xzz_xxxxyyz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxyyz_0[j] - tg_zz_xxxxyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzz_xxxyyz_1[j];

                    tg_xxzz_xxxxyzz_0[j] = pb_x * tg_xzz_xxxxyzz_0[j] + fr * tg_xzz_xxxxyzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxyzz_0[j] - tg_zz_xxxxyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzz_xxxyzz_1[j];

                    tg_xxzz_xxxxzzz_0[j] = pb_x * tg_xzz_xxxxzzz_0[j] + fr * tg_xzz_xxxxzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxzzz_0[j] - tg_zz_xxxxzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzz_xxxzzz_1[j];

                    tg_xxzz_xxxyyyy_0[j] = pb_x * tg_xzz_xxxyyyy_0[j] + fr * tg_xzz_xxxyyyy_1[j] + 0.5 * fl1_fx * (tg_zz_xxxyyyy_0[j] - tg_zz_xxxyyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzz_xxyyyy_1[j];

                    tg_xxzz_xxxyyyz_0[j] = pb_x * tg_xzz_xxxyyyz_0[j] + fr * tg_xzz_xxxyyyz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxyyyz_0[j] - tg_zz_xxxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzz_xxyyyz_1[j];

                    tg_xxzz_xxxyyzz_0[j] = pb_x * tg_xzz_xxxyyzz_0[j] + fr * tg_xzz_xxxyyzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxyyzz_0[j] - tg_zz_xxxyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzz_xxyyzz_1[j];

                    tg_xxzz_xxxyzzz_0[j] = pb_x * tg_xzz_xxxyzzz_0[j] + fr * tg_xzz_xxxyzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxyzzz_0[j] - tg_zz_xxxyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzz_xxyzzz_1[j];

                    tg_xxzz_xxxzzzz_0[j] = pb_x * tg_xzz_xxxzzzz_0[j] + fr * tg_xzz_xxxzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxzzzz_0[j] - tg_zz_xxxzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzz_xxzzzz_1[j];

                    tg_xxzz_xxyyyyy_0[j] = pb_x * tg_xzz_xxyyyyy_0[j] + fr * tg_xzz_xxyyyyy_1[j] + 0.5 * fl1_fx * (tg_zz_xxyyyyy_0[j] - tg_zz_xxyyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xzz_xyyyyy_1[j];

                    tg_xxzz_xxyyyyz_0[j] = pb_x * tg_xzz_xxyyyyz_0[j] + fr * tg_xzz_xxyyyyz_1[j] + 0.5 * fl1_fx * (tg_zz_xxyyyyz_0[j] - tg_zz_xxyyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xzz_xyyyyz_1[j];

                    tg_xxzz_xxyyyzz_0[j] = pb_x * tg_xzz_xxyyyzz_0[j] + fr * tg_xzz_xxyyyzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxyyyzz_0[j] - tg_zz_xxyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xzz_xyyyzz_1[j];

                    tg_xxzz_xxyyzzz_0[j] = pb_x * tg_xzz_xxyyzzz_0[j] + fr * tg_xzz_xxyyzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxyyzzz_0[j] - tg_zz_xxyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xzz_xyyzzz_1[j];

                    tg_xxzz_xxyzzzz_0[j] = pb_x * tg_xzz_xxyzzzz_0[j] + fr * tg_xzz_xxyzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxyzzzz_0[j] - tg_zz_xxyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xzz_xyzzzz_1[j];

                    tg_xxzz_xxzzzzz_0[j] = pb_x * tg_xzz_xxzzzzz_0[j] + fr * tg_xzz_xxzzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxzzzzz_0[j] - tg_zz_xxzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xzz_xzzzzz_1[j];

                    tg_xxzz_xyyyyyy_0[j] = pb_x * tg_xzz_xyyyyyy_0[j] + fr * tg_xzz_xyyyyyy_1[j] + 0.5 * fl1_fx * (tg_zz_xyyyyyy_0[j] - tg_zz_xyyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzz_yyyyyy_1[j];

                    tg_xxzz_xyyyyyz_0[j] = pb_x * tg_xzz_xyyyyyz_0[j] + fr * tg_xzz_xyyyyyz_1[j] + 0.5 * fl1_fx * (tg_zz_xyyyyyz_0[j] - tg_zz_xyyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzz_yyyyyz_1[j];

                    tg_xxzz_xyyyyzz_0[j] = pb_x * tg_xzz_xyyyyzz_0[j] + fr * tg_xzz_xyyyyzz_1[j] + 0.5 * fl1_fx * (tg_zz_xyyyyzz_0[j] - tg_zz_xyyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzz_yyyyzz_1[j];

                    tg_xxzz_xyyyzzz_0[j] = pb_x * tg_xzz_xyyyzzz_0[j] + fr * tg_xzz_xyyyzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xyyyzzz_0[j] - tg_zz_xyyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzz_yyyzzz_1[j];

                    tg_xxzz_xyyzzzz_0[j] = pb_x * tg_xzz_xyyzzzz_0[j] + fr * tg_xzz_xyyzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xyyzzzz_0[j] - tg_zz_xyyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzz_yyzzzz_1[j];

                    tg_xxzz_xyzzzzz_0[j] = pb_x * tg_xzz_xyzzzzz_0[j] + fr * tg_xzz_xyzzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xyzzzzz_0[j] - tg_zz_xyzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzz_yzzzzz_1[j];

                    tg_xxzz_xzzzzzz_0[j] = pb_x * tg_xzz_xzzzzzz_0[j] + fr * tg_xzz_xzzzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xzzzzzz_0[j] - tg_zz_xzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzz_zzzzzz_1[j];

                    tg_xxzz_yyyyyyy_0[j] = pb_x * tg_xzz_yyyyyyy_0[j] + fr * tg_xzz_yyyyyyy_1[j] + 0.5 * fl1_fx * (tg_zz_yyyyyyy_0[j] - tg_zz_yyyyyyy_1[j] * fl1_fza);

                    tg_xxzz_yyyyyyz_0[j] = pb_x * tg_xzz_yyyyyyz_0[j] + fr * tg_xzz_yyyyyyz_1[j] + 0.5 * fl1_fx * (tg_zz_yyyyyyz_0[j] - tg_zz_yyyyyyz_1[j] * fl1_fza);

                    tg_xxzz_yyyyyzz_0[j] = pb_x * tg_xzz_yyyyyzz_0[j] + fr * tg_xzz_yyyyyzz_1[j] + 0.5 * fl1_fx * (tg_zz_yyyyyzz_0[j] - tg_zz_yyyyyzz_1[j] * fl1_fza);

                    tg_xxzz_yyyyzzz_0[j] = pb_x * tg_xzz_yyyyzzz_0[j] + fr * tg_xzz_yyyyzzz_1[j] + 0.5 * fl1_fx * (tg_zz_yyyyzzz_0[j] - tg_zz_yyyyzzz_1[j] * fl1_fza);

                    tg_xxzz_yyyzzzz_0[j] = pb_x * tg_xzz_yyyzzzz_0[j] + fr * tg_xzz_yyyzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_yyyzzzz_0[j] - tg_zz_yyyzzzz_1[j] * fl1_fza);

                    tg_xxzz_yyzzzzz_0[j] = pb_x * tg_xzz_yyzzzzz_0[j] + fr * tg_xzz_yyzzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_yyzzzzz_0[j] - tg_zz_yyzzzzz_1[j] * fl1_fza);

                    tg_xxzz_yzzzzzz_0[j] = pb_x * tg_xzz_yzzzzzz_0[j] + fr * tg_xzz_yzzzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_yzzzzzz_0[j] - tg_zz_yzzzzzz_1[j] * fl1_fza);

                    tg_xxzz_zzzzzzz_0[j] = pb_x * tg_xzz_zzzzzzz_0[j] + fr * tg_xzz_zzzzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_zzzzzzz_0[j] - tg_zz_zzzzzzz_1[j] * fl1_fza);

                    tg_xyyy_xxxxxxx_0[j] = pb_x * tg_yyy_xxxxxxx_0[j] + fr * tg_yyy_xxxxxxx_1[j] + 3.5 * fl1_fxn * tg_yyy_xxxxxx_1[j];

                    tg_xyyy_xxxxxxy_0[j] = pb_x * tg_yyy_xxxxxxy_0[j] + fr * tg_yyy_xxxxxxy_1[j] + 3.0 * fl1_fxn * tg_yyy_xxxxxy_1[j];

                    tg_xyyy_xxxxxxz_0[j] = pb_x * tg_yyy_xxxxxxz_0[j] + fr * tg_yyy_xxxxxxz_1[j] + 3.0 * fl1_fxn * tg_yyy_xxxxxz_1[j];

                    tg_xyyy_xxxxxyy_0[j] = pb_x * tg_yyy_xxxxxyy_0[j] + fr * tg_yyy_xxxxxyy_1[j] + 2.5 * fl1_fxn * tg_yyy_xxxxyy_1[j];

                    tg_xyyy_xxxxxyz_0[j] = pb_x * tg_yyy_xxxxxyz_0[j] + fr * tg_yyy_xxxxxyz_1[j] + 2.5 * fl1_fxn * tg_yyy_xxxxyz_1[j];

                    tg_xyyy_xxxxxzz_0[j] = pb_x * tg_yyy_xxxxxzz_0[j] + fr * tg_yyy_xxxxxzz_1[j] + 2.5 * fl1_fxn * tg_yyy_xxxxzz_1[j];

                    tg_xyyy_xxxxyyy_0[j] = pb_x * tg_yyy_xxxxyyy_0[j] + fr * tg_yyy_xxxxyyy_1[j] + 2.0 * fl1_fxn * tg_yyy_xxxyyy_1[j];

                    tg_xyyy_xxxxyyz_0[j] = pb_x * tg_yyy_xxxxyyz_0[j] + fr * tg_yyy_xxxxyyz_1[j] + 2.0 * fl1_fxn * tg_yyy_xxxyyz_1[j];

                    tg_xyyy_xxxxyzz_0[j] = pb_x * tg_yyy_xxxxyzz_0[j] + fr * tg_yyy_xxxxyzz_1[j] + 2.0 * fl1_fxn * tg_yyy_xxxyzz_1[j];

                    tg_xyyy_xxxxzzz_0[j] = pb_x * tg_yyy_xxxxzzz_0[j] + fr * tg_yyy_xxxxzzz_1[j] + 2.0 * fl1_fxn * tg_yyy_xxxzzz_1[j];

                    tg_xyyy_xxxyyyy_0[j] = pb_x * tg_yyy_xxxyyyy_0[j] + fr * tg_yyy_xxxyyyy_1[j] + 1.5 * fl1_fxn * tg_yyy_xxyyyy_1[j];

                    tg_xyyy_xxxyyyz_0[j] = pb_x * tg_yyy_xxxyyyz_0[j] + fr * tg_yyy_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_yyy_xxyyyz_1[j];

                    tg_xyyy_xxxyyzz_0[j] = pb_x * tg_yyy_xxxyyzz_0[j] + fr * tg_yyy_xxxyyzz_1[j] + 1.5 * fl1_fxn * tg_yyy_xxyyzz_1[j];

                    tg_xyyy_xxxyzzz_0[j] = pb_x * tg_yyy_xxxyzzz_0[j] + fr * tg_yyy_xxxyzzz_1[j] + 1.5 * fl1_fxn * tg_yyy_xxyzzz_1[j];

                    tg_xyyy_xxxzzzz_0[j] = pb_x * tg_yyy_xxxzzzz_0[j] + fr * tg_yyy_xxxzzzz_1[j] + 1.5 * fl1_fxn * tg_yyy_xxzzzz_1[j];

                    tg_xyyy_xxyyyyy_0[j] = pb_x * tg_yyy_xxyyyyy_0[j] + fr * tg_yyy_xxyyyyy_1[j] + fl1_fxn * tg_yyy_xyyyyy_1[j];

                    tg_xyyy_xxyyyyz_0[j] = pb_x * tg_yyy_xxyyyyz_0[j] + fr * tg_yyy_xxyyyyz_1[j] + fl1_fxn * tg_yyy_xyyyyz_1[j];

                    tg_xyyy_xxyyyzz_0[j] = pb_x * tg_yyy_xxyyyzz_0[j] + fr * tg_yyy_xxyyyzz_1[j] + fl1_fxn * tg_yyy_xyyyzz_1[j];

                    tg_xyyy_xxyyzzz_0[j] = pb_x * tg_yyy_xxyyzzz_0[j] + fr * tg_yyy_xxyyzzz_1[j] + fl1_fxn * tg_yyy_xyyzzz_1[j];

                    tg_xyyy_xxyzzzz_0[j] = pb_x * tg_yyy_xxyzzzz_0[j] + fr * tg_yyy_xxyzzzz_1[j] + fl1_fxn * tg_yyy_xyzzzz_1[j];

                    tg_xyyy_xxzzzzz_0[j] = pb_x * tg_yyy_xxzzzzz_0[j] + fr * tg_yyy_xxzzzzz_1[j] + fl1_fxn * tg_yyy_xzzzzz_1[j];

                    tg_xyyy_xyyyyyy_0[j] = pb_x * tg_yyy_xyyyyyy_0[j] + fr * tg_yyy_xyyyyyy_1[j] + 0.5 * fl1_fxn * tg_yyy_yyyyyy_1[j];

                    tg_xyyy_xyyyyyz_0[j] = pb_x * tg_yyy_xyyyyyz_0[j] + fr * tg_yyy_xyyyyyz_1[j] + 0.5 * fl1_fxn * tg_yyy_yyyyyz_1[j];

                    tg_xyyy_xyyyyzz_0[j] = pb_x * tg_yyy_xyyyyzz_0[j] + fr * tg_yyy_xyyyyzz_1[j] + 0.5 * fl1_fxn * tg_yyy_yyyyzz_1[j];

                    tg_xyyy_xyyyzzz_0[j] = pb_x * tg_yyy_xyyyzzz_0[j] + fr * tg_yyy_xyyyzzz_1[j] + 0.5 * fl1_fxn * tg_yyy_yyyzzz_1[j];

                    tg_xyyy_xyyzzzz_0[j] = pb_x * tg_yyy_xyyzzzz_0[j] + fr * tg_yyy_xyyzzzz_1[j] + 0.5 * fl1_fxn * tg_yyy_yyzzzz_1[j];

                    tg_xyyy_xyzzzzz_0[j] = pb_x * tg_yyy_xyzzzzz_0[j] + fr * tg_yyy_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyy_yzzzzz_1[j];

                    tg_xyyy_xzzzzzz_0[j] = pb_x * tg_yyy_xzzzzzz_0[j] + fr * tg_yyy_xzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyy_zzzzzz_1[j];

                    tg_xyyy_yyyyyyy_0[j] = pb_x * tg_yyy_yyyyyyy_0[j] + fr * tg_yyy_yyyyyyy_1[j];

                    tg_xyyy_yyyyyyz_0[j] = pb_x * tg_yyy_yyyyyyz_0[j] + fr * tg_yyy_yyyyyyz_1[j];

                    tg_xyyy_yyyyyzz_0[j] = pb_x * tg_yyy_yyyyyzz_0[j] + fr * tg_yyy_yyyyyzz_1[j];

                    tg_xyyy_yyyyzzz_0[j] = pb_x * tg_yyy_yyyyzzz_0[j] + fr * tg_yyy_yyyyzzz_1[j];

                    tg_xyyy_yyyzzzz_0[j] = pb_x * tg_yyy_yyyzzzz_0[j] + fr * tg_yyy_yyyzzzz_1[j];

                    tg_xyyy_yyzzzzz_0[j] = pb_x * tg_yyy_yyzzzzz_0[j] + fr * tg_yyy_yyzzzzz_1[j];

                    tg_xyyy_yzzzzzz_0[j] = pb_x * tg_yyy_yzzzzzz_0[j] + fr * tg_yyy_yzzzzzz_1[j];

                    tg_xyyy_zzzzzzz_0[j] = pb_x * tg_yyy_zzzzzzz_0[j] + fr * tg_yyy_zzzzzzz_1[j];

                    tg_xyyz_xxxxxxx_0[j] = pb_x * tg_yyz_xxxxxxx_0[j] + fr * tg_yyz_xxxxxxx_1[j] + 3.5 * fl1_fxn * tg_yyz_xxxxxx_1[j];

                    tg_xyyz_xxxxxxy_0[j] = pb_x * tg_yyz_xxxxxxy_0[j] + fr * tg_yyz_xxxxxxy_1[j] + 3.0 * fl1_fxn * tg_yyz_xxxxxy_1[j];

                    tg_xyyz_xxxxxxz_0[j] = pb_x * tg_yyz_xxxxxxz_0[j] + fr * tg_yyz_xxxxxxz_1[j] + 3.0 * fl1_fxn * tg_yyz_xxxxxz_1[j];

                    tg_xyyz_xxxxxyy_0[j] = pb_x * tg_yyz_xxxxxyy_0[j] + fr * tg_yyz_xxxxxyy_1[j] + 2.5 * fl1_fxn * tg_yyz_xxxxyy_1[j];

                    tg_xyyz_xxxxxyz_0[j] = pb_x * tg_yyz_xxxxxyz_0[j] + fr * tg_yyz_xxxxxyz_1[j] + 2.5 * fl1_fxn * tg_yyz_xxxxyz_1[j];

                    tg_xyyz_xxxxxzz_0[j] = pb_x * tg_yyz_xxxxxzz_0[j] + fr * tg_yyz_xxxxxzz_1[j] + 2.5 * fl1_fxn * tg_yyz_xxxxzz_1[j];

                    tg_xyyz_xxxxyyy_0[j] = pb_x * tg_yyz_xxxxyyy_0[j] + fr * tg_yyz_xxxxyyy_1[j] + 2.0 * fl1_fxn * tg_yyz_xxxyyy_1[j];

                    tg_xyyz_xxxxyyz_0[j] = pb_x * tg_yyz_xxxxyyz_0[j] + fr * tg_yyz_xxxxyyz_1[j] + 2.0 * fl1_fxn * tg_yyz_xxxyyz_1[j];

                    tg_xyyz_xxxxyzz_0[j] = pb_x * tg_yyz_xxxxyzz_0[j] + fr * tg_yyz_xxxxyzz_1[j] + 2.0 * fl1_fxn * tg_yyz_xxxyzz_1[j];

                    tg_xyyz_xxxxzzz_0[j] = pb_x * tg_yyz_xxxxzzz_0[j] + fr * tg_yyz_xxxxzzz_1[j] + 2.0 * fl1_fxn * tg_yyz_xxxzzz_1[j];

                    tg_xyyz_xxxyyyy_0[j] = pb_x * tg_yyz_xxxyyyy_0[j] + fr * tg_yyz_xxxyyyy_1[j] + 1.5 * fl1_fxn * tg_yyz_xxyyyy_1[j];

                    tg_xyyz_xxxyyyz_0[j] = pb_x * tg_yyz_xxxyyyz_0[j] + fr * tg_yyz_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_yyz_xxyyyz_1[j];

                    tg_xyyz_xxxyyzz_0[j] = pb_x * tg_yyz_xxxyyzz_0[j] + fr * tg_yyz_xxxyyzz_1[j] + 1.5 * fl1_fxn * tg_yyz_xxyyzz_1[j];

                    tg_xyyz_xxxyzzz_0[j] = pb_x * tg_yyz_xxxyzzz_0[j] + fr * tg_yyz_xxxyzzz_1[j] + 1.5 * fl1_fxn * tg_yyz_xxyzzz_1[j];

                    tg_xyyz_xxxzzzz_0[j] = pb_x * tg_yyz_xxxzzzz_0[j] + fr * tg_yyz_xxxzzzz_1[j] + 1.5 * fl1_fxn * tg_yyz_xxzzzz_1[j];

                    tg_xyyz_xxyyyyy_0[j] = pb_x * tg_yyz_xxyyyyy_0[j] + fr * tg_yyz_xxyyyyy_1[j] + fl1_fxn * tg_yyz_xyyyyy_1[j];

                    tg_xyyz_xxyyyyz_0[j] = pb_x * tg_yyz_xxyyyyz_0[j] + fr * tg_yyz_xxyyyyz_1[j] + fl1_fxn * tg_yyz_xyyyyz_1[j];

                    tg_xyyz_xxyyyzz_0[j] = pb_x * tg_yyz_xxyyyzz_0[j] + fr * tg_yyz_xxyyyzz_1[j] + fl1_fxn * tg_yyz_xyyyzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSK_270_360(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (270,360)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {4, -1, -1, -1},
                                             {7, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_7_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {7, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_7_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_7_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_7_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {6, -1, -1, -1}, 
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

                auto tg_yyz_xxyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 270); 

                auto tg_yyz_xxyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 271); 

                auto tg_yyz_xxzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 272); 

                auto tg_yyz_xyyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 273); 

                auto tg_yyz_xyyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 274); 

                auto tg_yyz_xyyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 275); 

                auto tg_yyz_xyyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 276); 

                auto tg_yyz_xyyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 277); 

                auto tg_yyz_xyzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 278); 

                auto tg_yyz_xzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 279); 

                auto tg_yyz_yyyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 280); 

                auto tg_yyz_yyyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 281); 

                auto tg_yyz_yyyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 282); 

                auto tg_yyz_yyyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 283); 

                auto tg_yyz_yyyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 284); 

                auto tg_yyz_yyzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 285); 

                auto tg_yyz_yzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 286); 

                auto tg_yyz_zzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 287); 

                auto tg_yzz_xxxxxxx_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 288); 

                auto tg_yzz_xxxxxxy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 289); 

                auto tg_yzz_xxxxxxz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 290); 

                auto tg_yzz_xxxxxyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 291); 

                auto tg_yzz_xxxxxyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 292); 

                auto tg_yzz_xxxxxzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 293); 

                auto tg_yzz_xxxxyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 294); 

                auto tg_yzz_xxxxyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 295); 

                auto tg_yzz_xxxxyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 296); 

                auto tg_yzz_xxxxzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 297); 

                auto tg_yzz_xxxyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 298); 

                auto tg_yzz_xxxyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 299); 

                auto tg_yzz_xxxyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 300); 

                auto tg_yzz_xxxyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 301); 

                auto tg_yzz_xxxzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 302); 

                auto tg_yzz_xxyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 303); 

                auto tg_yzz_xxyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 304); 

                auto tg_yzz_xxyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 305); 

                auto tg_yzz_xxyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 306); 

                auto tg_yzz_xxyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 307); 

                auto tg_yzz_xxzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 308); 

                auto tg_yzz_xyyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 309); 

                auto tg_yzz_xyyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 310); 

                auto tg_yzz_xyyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 311); 

                auto tg_yzz_xyyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 312); 

                auto tg_yzz_xyyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 313); 

                auto tg_yzz_xyzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 314); 

                auto tg_yzz_xzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 315); 

                auto tg_yzz_yyyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 316); 

                auto tg_yzz_yyyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 317); 

                auto tg_yzz_yyyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 318); 

                auto tg_yzz_yyyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 319); 

                auto tg_yzz_yyyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 320); 

                auto tg_yzz_yyzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 321); 

                auto tg_yzz_yzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 322); 

                auto tg_yzz_zzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 323); 

                auto tg_zzz_xxxxxxx_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 324); 

                auto tg_zzz_xxxxxxy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 325); 

                auto tg_zzz_xxxxxxz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 326); 

                auto tg_zzz_xxxxxyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 327); 

                auto tg_zzz_xxxxxyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 328); 

                auto tg_zzz_xxxxxzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 329); 

                auto tg_zzz_xxxxyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 330); 

                auto tg_zzz_xxxxyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 331); 

                auto tg_zzz_xxxxyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 332); 

                auto tg_zzz_xxxxzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 333); 

                auto tg_zzz_xxxyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 334); 

                auto tg_zzz_xxxyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 335); 

                auto tg_zzz_xxxyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 336); 

                auto tg_zzz_xxxyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 337); 

                auto tg_zzz_xxxzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 338); 

                auto tg_zzz_xxyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 339); 

                auto tg_zzz_xxyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 340); 

                auto tg_zzz_xxyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 341); 

                auto tg_zzz_xxyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 342); 

                auto tg_zzz_xxyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 343); 

                auto tg_zzz_xxzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 344); 

                auto tg_zzz_xyyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 345); 

                auto tg_zzz_xyyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 346); 

                auto tg_zzz_xyyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 347); 

                auto tg_zzz_xyyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 348); 

                auto tg_zzz_xyyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 349); 

                auto tg_zzz_xyzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 350); 

                auto tg_zzz_xzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 351); 

                auto tg_zzz_yyyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 352); 

                auto tg_zzz_yyyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 353); 

                auto tg_zzz_yyyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 354); 

                auto tg_zzz_yyyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 355); 

                auto tg_zzz_yyyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 356); 

                auto tg_zzz_yyzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 357); 

                auto tg_zzz_yzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 358); 

                auto tg_zzz_zzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 359); 

                auto tg_yyz_xxyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 270); 

                auto tg_yyz_xxyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 271); 

                auto tg_yyz_xxzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 272); 

                auto tg_yyz_xyyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 273); 

                auto tg_yyz_xyyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 274); 

                auto tg_yyz_xyyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 275); 

                auto tg_yyz_xyyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 276); 

                auto tg_yyz_xyyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 277); 

                auto tg_yyz_xyzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 278); 

                auto tg_yyz_xzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 279); 

                auto tg_yyz_yyyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 280); 

                auto tg_yyz_yyyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 281); 

                auto tg_yyz_yyyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 282); 

                auto tg_yyz_yyyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 283); 

                auto tg_yyz_yyyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 284); 

                auto tg_yyz_yyzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 285); 

                auto tg_yyz_yzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 286); 

                auto tg_yyz_zzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 287); 

                auto tg_yzz_xxxxxxx_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 288); 

                auto tg_yzz_xxxxxxy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 289); 

                auto tg_yzz_xxxxxxz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 290); 

                auto tg_yzz_xxxxxyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 291); 

                auto tg_yzz_xxxxxyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 292); 

                auto tg_yzz_xxxxxzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 293); 

                auto tg_yzz_xxxxyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 294); 

                auto tg_yzz_xxxxyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 295); 

                auto tg_yzz_xxxxyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 296); 

                auto tg_yzz_xxxxzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 297); 

                auto tg_yzz_xxxyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 298); 

                auto tg_yzz_xxxyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 299); 

                auto tg_yzz_xxxyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 300); 

                auto tg_yzz_xxxyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 301); 

                auto tg_yzz_xxxzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 302); 

                auto tg_yzz_xxyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 303); 

                auto tg_yzz_xxyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 304); 

                auto tg_yzz_xxyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 305); 

                auto tg_yzz_xxyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 306); 

                auto tg_yzz_xxyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 307); 

                auto tg_yzz_xxzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 308); 

                auto tg_yzz_xyyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 309); 

                auto tg_yzz_xyyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 310); 

                auto tg_yzz_xyyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 311); 

                auto tg_yzz_xyyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 312); 

                auto tg_yzz_xyyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 313); 

                auto tg_yzz_xyzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 314); 

                auto tg_yzz_xzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 315); 

                auto tg_yzz_yyyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 316); 

                auto tg_yzz_yyyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 317); 

                auto tg_yzz_yyyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 318); 

                auto tg_yzz_yyyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 319); 

                auto tg_yzz_yyyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 320); 

                auto tg_yzz_yyzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 321); 

                auto tg_yzz_yzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 322); 

                auto tg_yzz_zzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 323); 

                auto tg_zzz_xxxxxxx_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 324); 

                auto tg_zzz_xxxxxxy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 325); 

                auto tg_zzz_xxxxxxz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 326); 

                auto tg_zzz_xxxxxyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 327); 

                auto tg_zzz_xxxxxyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 328); 

                auto tg_zzz_xxxxxzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 329); 

                auto tg_zzz_xxxxyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 330); 

                auto tg_zzz_xxxxyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 331); 

                auto tg_zzz_xxxxyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 332); 

                auto tg_zzz_xxxxzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 333); 

                auto tg_zzz_xxxyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 334); 

                auto tg_zzz_xxxyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 335); 

                auto tg_zzz_xxxyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 336); 

                auto tg_zzz_xxxyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 337); 

                auto tg_zzz_xxxzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 338); 

                auto tg_zzz_xxyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 339); 

                auto tg_zzz_xxyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 340); 

                auto tg_zzz_xxyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 341); 

                auto tg_zzz_xxyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 342); 

                auto tg_zzz_xxyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 343); 

                auto tg_zzz_xxzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 344); 

                auto tg_zzz_xyyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 345); 

                auto tg_zzz_xyyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 346); 

                auto tg_zzz_xyyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 347); 

                auto tg_zzz_xyyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 348); 

                auto tg_zzz_xyyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 349); 

                auto tg_zzz_xyzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 350); 

                auto tg_zzz_xzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 351); 

                auto tg_zzz_yyyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 352); 

                auto tg_zzz_yyyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 353); 

                auto tg_zzz_yyyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 354); 

                auto tg_zzz_yyyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 355); 

                auto tg_zzz_yyyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 356); 

                auto tg_zzz_yyzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 357); 

                auto tg_zzz_yzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 358); 

                auto tg_zzz_zzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 359); 

                auto tg_yyz_xyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 214); 

                auto tg_yyz_xyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 215); 

                auto tg_yyz_xzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 216); 

                auto tg_yyz_yyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 217); 

                auto tg_yyz_yyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 218); 

                auto tg_yyz_yyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 219); 

                auto tg_yyz_yyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 220); 

                auto tg_yyz_yyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 221); 

                auto tg_yyz_yzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 222); 

                auto tg_yyz_zzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 223); 

                auto tg_yzz_xxxxxx_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 224); 

                auto tg_yzz_xxxxxy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 225); 

                auto tg_yzz_xxxxxz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 226); 

                auto tg_yzz_xxxxyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 227); 

                auto tg_yzz_xxxxyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 228); 

                auto tg_yzz_xxxxzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 229); 

                auto tg_yzz_xxxyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 230); 

                auto tg_yzz_xxxyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 231); 

                auto tg_yzz_xxxyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 232); 

                auto tg_yzz_xxxzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 233); 

                auto tg_yzz_xxyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 234); 

                auto tg_yzz_xxyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 235); 

                auto tg_yzz_xxyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 236); 

                auto tg_yzz_xxyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 237); 

                auto tg_yzz_xxzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 238); 

                auto tg_yzz_xyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 239); 

                auto tg_yzz_xyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 240); 

                auto tg_yzz_xyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 241); 

                auto tg_yzz_xyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 242); 

                auto tg_yzz_xyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 243); 

                auto tg_yzz_xzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 244); 

                auto tg_yzz_yyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 245); 

                auto tg_yzz_yyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 246); 

                auto tg_yzz_yyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 247); 

                auto tg_yzz_yyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 248); 

                auto tg_yzz_yyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 249); 

                auto tg_yzz_yzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 250); 

                auto tg_yzz_zzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 251); 

                auto tg_zzz_xxxxxx_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 252); 

                auto tg_zzz_xxxxxy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 253); 

                auto tg_zzz_xxxxxz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 254); 

                auto tg_zzz_xxxxyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 255); 

                auto tg_zzz_xxxxyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 256); 

                auto tg_zzz_xxxxzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 257); 

                auto tg_zzz_xxxyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 258); 

                auto tg_zzz_xxxyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 259); 

                auto tg_zzz_xxxyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 260); 

                auto tg_zzz_xxxzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 261); 

                auto tg_zzz_xxyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 262); 

                auto tg_zzz_xxyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 263); 

                auto tg_zzz_xxyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 264); 

                auto tg_zzz_xxyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 265); 

                auto tg_zzz_xxzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 266); 

                auto tg_zzz_xyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 267); 

                auto tg_zzz_xyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 268); 

                auto tg_zzz_xyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 269); 

                auto tg_zzz_xyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 270); 

                auto tg_zzz_xyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 271); 

                auto tg_zzz_xzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 272); 

                auto tg_zzz_yyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 273); 

                auto tg_zzz_yyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 274); 

                auto tg_zzz_yyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 275); 

                auto tg_zzz_yyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 276); 

                auto tg_zzz_yyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 277); 

                auto tg_zzz_yzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 278); 

                auto tg_zzz_zzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 279); 

                // set up pointers to integrals

                auto tg_xyyz_xxyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 270); 

                auto tg_xyyz_xxyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 271); 

                auto tg_xyyz_xxzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 272); 

                auto tg_xyyz_xyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 273); 

                auto tg_xyyz_xyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 274); 

                auto tg_xyyz_xyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 275); 

                auto tg_xyyz_xyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 276); 

                auto tg_xyyz_xyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 277); 

                auto tg_xyyz_xyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 278); 

                auto tg_xyyz_xzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 279); 

                auto tg_xyyz_yyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 280); 

                auto tg_xyyz_yyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 281); 

                auto tg_xyyz_yyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 282); 

                auto tg_xyyz_yyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 283); 

                auto tg_xyyz_yyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 284); 

                auto tg_xyyz_yyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 285); 

                auto tg_xyyz_yzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 286); 

                auto tg_xyyz_zzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 287); 

                auto tg_xyzz_xxxxxxx_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 288); 

                auto tg_xyzz_xxxxxxy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 289); 

                auto tg_xyzz_xxxxxxz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 290); 

                auto tg_xyzz_xxxxxyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 291); 

                auto tg_xyzz_xxxxxyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 292); 

                auto tg_xyzz_xxxxxzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 293); 

                auto tg_xyzz_xxxxyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 294); 

                auto tg_xyzz_xxxxyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 295); 

                auto tg_xyzz_xxxxyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 296); 

                auto tg_xyzz_xxxxzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 297); 

                auto tg_xyzz_xxxyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 298); 

                auto tg_xyzz_xxxyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 299); 

                auto tg_xyzz_xxxyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 300); 

                auto tg_xyzz_xxxyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 301); 

                auto tg_xyzz_xxxzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 302); 

                auto tg_xyzz_xxyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 303); 

                auto tg_xyzz_xxyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 304); 

                auto tg_xyzz_xxyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 305); 

                auto tg_xyzz_xxyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 306); 

                auto tg_xyzz_xxyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 307); 

                auto tg_xyzz_xxzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 308); 

                auto tg_xyzz_xyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 309); 

                auto tg_xyzz_xyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 310); 

                auto tg_xyzz_xyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 311); 

                auto tg_xyzz_xyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 312); 

                auto tg_xyzz_xyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 313); 

                auto tg_xyzz_xyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 314); 

                auto tg_xyzz_xzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 315); 

                auto tg_xyzz_yyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 316); 

                auto tg_xyzz_yyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 317); 

                auto tg_xyzz_yyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 318); 

                auto tg_xyzz_yyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 319); 

                auto tg_xyzz_yyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 320); 

                auto tg_xyzz_yyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 321); 

                auto tg_xyzz_yzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 322); 

                auto tg_xyzz_zzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 323); 

                auto tg_xzzz_xxxxxxx_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 324); 

                auto tg_xzzz_xxxxxxy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 325); 

                auto tg_xzzz_xxxxxxz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 326); 

                auto tg_xzzz_xxxxxyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 327); 

                auto tg_xzzz_xxxxxyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 328); 

                auto tg_xzzz_xxxxxzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 329); 

                auto tg_xzzz_xxxxyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 330); 

                auto tg_xzzz_xxxxyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 331); 

                auto tg_xzzz_xxxxyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 332); 

                auto tg_xzzz_xxxxzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 333); 

                auto tg_xzzz_xxxyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 334); 

                auto tg_xzzz_xxxyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 335); 

                auto tg_xzzz_xxxyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 336); 

                auto tg_xzzz_xxxyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 337); 

                auto tg_xzzz_xxxzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 338); 

                auto tg_xzzz_xxyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 339); 

                auto tg_xzzz_xxyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 340); 

                auto tg_xzzz_xxyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 341); 

                auto tg_xzzz_xxyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 342); 

                auto tg_xzzz_xxyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 343); 

                auto tg_xzzz_xxzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 344); 

                auto tg_xzzz_xyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 345); 

                auto tg_xzzz_xyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 346); 

                auto tg_xzzz_xyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 347); 

                auto tg_xzzz_xyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 348); 

                auto tg_xzzz_xyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 349); 

                auto tg_xzzz_xyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 350); 

                auto tg_xzzz_xzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 351); 

                auto tg_xzzz_yyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 352); 

                auto tg_xzzz_yyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 353); 

                auto tg_xzzz_yyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 354); 

                auto tg_xzzz_yyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 355); 

                auto tg_xzzz_yyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 356); 

                auto tg_xzzz_yyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 357); 

                auto tg_xzzz_yzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 358); 

                auto tg_xzzz_zzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 359); 

                // Batch of Integrals (270,360)

                #pragma omp simd aligned(fxn, tg_xyyz_xxyyzzz_0, tg_xyyz_xxyzzzz_0, tg_xyyz_xxzzzzz_0, \
                                         tg_xyyz_xyyyyyy_0, tg_xyyz_xyyyyyz_0, tg_xyyz_xyyyyzz_0, tg_xyyz_xyyyzzz_0, \
                                         tg_xyyz_xyyzzzz_0, tg_xyyz_xyzzzzz_0, tg_xyyz_xzzzzzz_0, tg_xyyz_yyyyyyy_0, \
                                         tg_xyyz_yyyyyyz_0, tg_xyyz_yyyyyzz_0, tg_xyyz_yyyyzzz_0, tg_xyyz_yyyzzzz_0, \
                                         tg_xyyz_yyzzzzz_0, tg_xyyz_yzzzzzz_0, tg_xyyz_zzzzzzz_0, tg_xyzz_xxxxxxx_0, \
                                         tg_xyzz_xxxxxxy_0, tg_xyzz_xxxxxxz_0, tg_xyzz_xxxxxyy_0, tg_xyzz_xxxxxyz_0, \
                                         tg_xyzz_xxxxxzz_0, tg_xyzz_xxxxyyy_0, tg_xyzz_xxxxyyz_0, tg_xyzz_xxxxyzz_0, \
                                         tg_xyzz_xxxxzzz_0, tg_xyzz_xxxyyyy_0, tg_xyzz_xxxyyyz_0, tg_xyzz_xxxyyzz_0, \
                                         tg_xyzz_xxxyzzz_0, tg_xyzz_xxxzzzz_0, tg_xyzz_xxyyyyy_0, tg_xyzz_xxyyyyz_0, \
                                         tg_xyzz_xxyyyzz_0, tg_xyzz_xxyyzzz_0, tg_xyzz_xxyzzzz_0, tg_xyzz_xxzzzzz_0, \
                                         tg_xyzz_xyyyyyy_0, tg_xyzz_xyyyyyz_0, tg_xyzz_xyyyyzz_0, tg_xyzz_xyyyzzz_0, \
                                         tg_xyzz_xyyzzzz_0, tg_xyzz_xyzzzzz_0, tg_xyzz_xzzzzzz_0, tg_xyzz_yyyyyyy_0, \
                                         tg_xyzz_yyyyyyz_0, tg_xyzz_yyyyyzz_0, tg_xyzz_yyyyzzz_0, tg_xyzz_yyyzzzz_0, \
                                         tg_xyzz_yyzzzzz_0, tg_xyzz_yzzzzzz_0, tg_xyzz_zzzzzzz_0, tg_xzzz_xxxxxxx_0, \
                                         tg_xzzz_xxxxxxy_0, tg_xzzz_xxxxxxz_0, tg_xzzz_xxxxxyy_0, tg_xzzz_xxxxxyz_0, \
                                         tg_xzzz_xxxxxzz_0, tg_xzzz_xxxxyyy_0, tg_xzzz_xxxxyyz_0, tg_xzzz_xxxxyzz_0, \
                                         tg_xzzz_xxxxzzz_0, tg_xzzz_xxxyyyy_0, tg_xzzz_xxxyyyz_0, tg_xzzz_xxxyyzz_0, \
                                         tg_xzzz_xxxyzzz_0, tg_xzzz_xxxzzzz_0, tg_xzzz_xxyyyyy_0, tg_xzzz_xxyyyyz_0, \
                                         tg_xzzz_xxyyyzz_0, tg_xzzz_xxyyzzz_0, tg_xzzz_xxyzzzz_0, tg_xzzz_xxzzzzz_0, \
                                         tg_xzzz_xyyyyyy_0, tg_xzzz_xyyyyyz_0, tg_xzzz_xyyyyzz_0, tg_xzzz_xyyyzzz_0, \
                                         tg_xzzz_xyyzzzz_0, tg_xzzz_xyzzzzz_0, tg_xzzz_xzzzzzz_0, tg_xzzz_yyyyyyy_0, \
                                         tg_xzzz_yyyyyyz_0, tg_xzzz_yyyyyzz_0, tg_xzzz_yyyyzzz_0, tg_xzzz_yyyzzzz_0, \
                                         tg_xzzz_yyzzzzz_0, tg_xzzz_yzzzzzz_0, tg_xzzz_zzzzzzz_0, tg_yyz_xxyyzzz_0, \
                                         tg_yyz_xxyyzzz_1, tg_yyz_xxyzzzz_0, tg_yyz_xxyzzzz_1, tg_yyz_xxzzzzz_0, \
                                         tg_yyz_xxzzzzz_1, tg_yyz_xyyyyyy_0, tg_yyz_xyyyyyy_1, tg_yyz_xyyyyyz_0, \
                                         tg_yyz_xyyyyyz_1, tg_yyz_xyyyyzz_0, tg_yyz_xyyyyzz_1, tg_yyz_xyyyzzz_0, \
                                         tg_yyz_xyyyzzz_1, tg_yyz_xyyzzz_1, tg_yyz_xyyzzzz_0, tg_yyz_xyyzzzz_1, \
                                         tg_yyz_xyzzzz_1, tg_yyz_xyzzzzz_0, tg_yyz_xyzzzzz_1, tg_yyz_xzzzzz_1, \
                                         tg_yyz_xzzzzzz_0, tg_yyz_xzzzzzz_1, tg_yyz_yyyyyy_1, tg_yyz_yyyyyyy_0, \
                                         tg_yyz_yyyyyyy_1, tg_yyz_yyyyyyz_0, tg_yyz_yyyyyyz_1, tg_yyz_yyyyyz_1, \
                                         tg_yyz_yyyyyzz_0, tg_yyz_yyyyyzz_1, tg_yyz_yyyyzz_1, tg_yyz_yyyyzzz_0, \
                                         tg_yyz_yyyyzzz_1, tg_yyz_yyyzzz_1, tg_yyz_yyyzzzz_0, tg_yyz_yyyzzzz_1, \
                                         tg_yyz_yyzzzz_1, tg_yyz_yyzzzzz_0, tg_yyz_yyzzzzz_1, tg_yyz_yzzzzz_1, \
                                         tg_yyz_yzzzzzz_0, tg_yyz_yzzzzzz_1, tg_yyz_zzzzzz_1, tg_yyz_zzzzzzz_0, \
                                         tg_yyz_zzzzzzz_1, tg_yzz_xxxxxx_1, tg_yzz_xxxxxxx_0, tg_yzz_xxxxxxx_1, \
                                         tg_yzz_xxxxxxy_0, tg_yzz_xxxxxxy_1, tg_yzz_xxxxxxz_0, tg_yzz_xxxxxxz_1, \
                                         tg_yzz_xxxxxy_1, tg_yzz_xxxxxyy_0, tg_yzz_xxxxxyy_1, tg_yzz_xxxxxyz_0, \
                                         tg_yzz_xxxxxyz_1, tg_yzz_xxxxxz_1, tg_yzz_xxxxxzz_0, tg_yzz_xxxxxzz_1, \
                                         tg_yzz_xxxxyy_1, tg_yzz_xxxxyyy_0, tg_yzz_xxxxyyy_1, tg_yzz_xxxxyyz_0, \
                                         tg_yzz_xxxxyyz_1, tg_yzz_xxxxyz_1, tg_yzz_xxxxyzz_0, tg_yzz_xxxxyzz_1, \
                                         tg_yzz_xxxxzz_1, tg_yzz_xxxxzzz_0, tg_yzz_xxxxzzz_1, tg_yzz_xxxyyy_1, \
                                         tg_yzz_xxxyyyy_0, tg_yzz_xxxyyyy_1, tg_yzz_xxxyyyz_0, tg_yzz_xxxyyyz_1, \
                                         tg_yzz_xxxyyz_1, tg_yzz_xxxyyzz_0, tg_yzz_xxxyyzz_1, tg_yzz_xxxyzz_1, \
                                         tg_yzz_xxxyzzz_0, tg_yzz_xxxyzzz_1, tg_yzz_xxxzzz_1, tg_yzz_xxxzzzz_0, \
                                         tg_yzz_xxxzzzz_1, tg_yzz_xxyyyy_1, tg_yzz_xxyyyyy_0, tg_yzz_xxyyyyy_1, \
                                         tg_yzz_xxyyyyz_0, tg_yzz_xxyyyyz_1, tg_yzz_xxyyyz_1, tg_yzz_xxyyyzz_0, \
                                         tg_yzz_xxyyyzz_1, tg_yzz_xxyyzz_1, tg_yzz_xxyyzzz_0, tg_yzz_xxyyzzz_1, \
                                         tg_yzz_xxyzzz_1, tg_yzz_xxyzzzz_0, tg_yzz_xxyzzzz_1, tg_yzz_xxzzzz_1, \
                                         tg_yzz_xxzzzzz_0, tg_yzz_xxzzzzz_1, tg_yzz_xyyyyy_1, tg_yzz_xyyyyyy_0, \
                                         tg_yzz_xyyyyyy_1, tg_yzz_xyyyyyz_0, tg_yzz_xyyyyyz_1, tg_yzz_xyyyyz_1, \
                                         tg_yzz_xyyyyzz_0, tg_yzz_xyyyyzz_1, tg_yzz_xyyyzz_1, tg_yzz_xyyyzzz_0, \
                                         tg_yzz_xyyyzzz_1, tg_yzz_xyyzzz_1, tg_yzz_xyyzzzz_0, tg_yzz_xyyzzzz_1, \
                                         tg_yzz_xyzzzz_1, tg_yzz_xyzzzzz_0, tg_yzz_xyzzzzz_1, tg_yzz_xzzzzz_1, \
                                         tg_yzz_xzzzzzz_0, tg_yzz_xzzzzzz_1, tg_yzz_yyyyyy_1, tg_yzz_yyyyyyy_0, \
                                         tg_yzz_yyyyyyy_1, tg_yzz_yyyyyyz_0, tg_yzz_yyyyyyz_1, tg_yzz_yyyyyz_1, \
                                         tg_yzz_yyyyyzz_0, tg_yzz_yyyyyzz_1, tg_yzz_yyyyzz_1, tg_yzz_yyyyzzz_0, \
                                         tg_yzz_yyyyzzz_1, tg_yzz_yyyzzz_1, tg_yzz_yyyzzzz_0, tg_yzz_yyyzzzz_1, \
                                         tg_yzz_yyzzzz_1, tg_yzz_yyzzzzz_0, tg_yzz_yyzzzzz_1, tg_yzz_yzzzzz_1, \
                                         tg_yzz_yzzzzzz_0, tg_yzz_yzzzzzz_1, tg_yzz_zzzzzz_1, tg_yzz_zzzzzzz_0, \
                                         tg_yzz_zzzzzzz_1, tg_zzz_xxxxxx_1, tg_zzz_xxxxxxx_0, tg_zzz_xxxxxxx_1, \
                                         tg_zzz_xxxxxxy_0, tg_zzz_xxxxxxy_1, tg_zzz_xxxxxxz_0, tg_zzz_xxxxxxz_1, \
                                         tg_zzz_xxxxxy_1, tg_zzz_xxxxxyy_0, tg_zzz_xxxxxyy_1, tg_zzz_xxxxxyz_0, \
                                         tg_zzz_xxxxxyz_1, tg_zzz_xxxxxz_1, tg_zzz_xxxxxzz_0, tg_zzz_xxxxxzz_1, \
                                         tg_zzz_xxxxyy_1, tg_zzz_xxxxyyy_0, tg_zzz_xxxxyyy_1, tg_zzz_xxxxyyz_0, \
                                         tg_zzz_xxxxyyz_1, tg_zzz_xxxxyz_1, tg_zzz_xxxxyzz_0, tg_zzz_xxxxyzz_1, \
                                         tg_zzz_xxxxzz_1, tg_zzz_xxxxzzz_0, tg_zzz_xxxxzzz_1, tg_zzz_xxxyyy_1, \
                                         tg_zzz_xxxyyyy_0, tg_zzz_xxxyyyy_1, tg_zzz_xxxyyyz_0, tg_zzz_xxxyyyz_1, \
                                         tg_zzz_xxxyyz_1, tg_zzz_xxxyyzz_0, tg_zzz_xxxyyzz_1, tg_zzz_xxxyzz_1, \
                                         tg_zzz_xxxyzzz_0, tg_zzz_xxxyzzz_1, tg_zzz_xxxzzz_1, tg_zzz_xxxzzzz_0, \
                                         tg_zzz_xxxzzzz_1, tg_zzz_xxyyyy_1, tg_zzz_xxyyyyy_0, tg_zzz_xxyyyyy_1, \
                                         tg_zzz_xxyyyyz_0, tg_zzz_xxyyyyz_1, tg_zzz_xxyyyz_1, tg_zzz_xxyyyzz_0, \
                                         tg_zzz_xxyyyzz_1, tg_zzz_xxyyzz_1, tg_zzz_xxyyzzz_0, tg_zzz_xxyyzzz_1, \
                                         tg_zzz_xxyzzz_1, tg_zzz_xxyzzzz_0, tg_zzz_xxyzzzz_1, tg_zzz_xxzzzz_1, \
                                         tg_zzz_xxzzzzz_0, tg_zzz_xxzzzzz_1, tg_zzz_xyyyyy_1, tg_zzz_xyyyyyy_0, \
                                         tg_zzz_xyyyyyy_1, tg_zzz_xyyyyyz_0, tg_zzz_xyyyyyz_1, tg_zzz_xyyyyz_1, \
                                         tg_zzz_xyyyyzz_0, tg_zzz_xyyyyzz_1, tg_zzz_xyyyzz_1, tg_zzz_xyyyzzz_0, \
                                         tg_zzz_xyyyzzz_1, tg_zzz_xyyzzz_1, tg_zzz_xyyzzzz_0, tg_zzz_xyyzzzz_1, \
                                         tg_zzz_xyzzzz_1, tg_zzz_xyzzzzz_0, tg_zzz_xyzzzzz_1, tg_zzz_xzzzzz_1, \
                                         tg_zzz_xzzzzzz_0, tg_zzz_xzzzzzz_1, tg_zzz_yyyyyy_1, tg_zzz_yyyyyyy_0, \
                                         tg_zzz_yyyyyyy_1, tg_zzz_yyyyyyz_0, tg_zzz_yyyyyyz_1, tg_zzz_yyyyyz_1, \
                                         tg_zzz_yyyyyzz_0, tg_zzz_yyyyyzz_1, tg_zzz_yyyyzz_1, tg_zzz_yyyyzzz_0, \
                                         tg_zzz_yyyyzzz_1, tg_zzz_yyyzzz_1, tg_zzz_yyyzzzz_0, tg_zzz_yyyzzzz_1, \
                                         tg_zzz_yyzzzz_1, tg_zzz_yyzzzzz_0, tg_zzz_yyzzzzz_1, tg_zzz_yzzzzz_1, \
                                         tg_zzz_yzzzzzz_0, tg_zzz_yzzzzzz_1, tg_zzz_zzzzzz_1, tg_zzz_zzzzzzz_0, \
                                         tg_zzz_zzzzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fxn = fxn[j];

                    double fr = wp_x[j]; 

                    tg_xyyz_xxyyzzz_0[j] = pb_x * tg_yyz_xxyyzzz_0[j] + fr * tg_yyz_xxyyzzz_1[j] + fl1_fxn * tg_yyz_xyyzzz_1[j];

                    tg_xyyz_xxyzzzz_0[j] = pb_x * tg_yyz_xxyzzzz_0[j] + fr * tg_yyz_xxyzzzz_1[j] + fl1_fxn * tg_yyz_xyzzzz_1[j];

                    tg_xyyz_xxzzzzz_0[j] = pb_x * tg_yyz_xxzzzzz_0[j] + fr * tg_yyz_xxzzzzz_1[j] + fl1_fxn * tg_yyz_xzzzzz_1[j];

                    tg_xyyz_xyyyyyy_0[j] = pb_x * tg_yyz_xyyyyyy_0[j] + fr * tg_yyz_xyyyyyy_1[j] + 0.5 * fl1_fxn * tg_yyz_yyyyyy_1[j];

                    tg_xyyz_xyyyyyz_0[j] = pb_x * tg_yyz_xyyyyyz_0[j] + fr * tg_yyz_xyyyyyz_1[j] + 0.5 * fl1_fxn * tg_yyz_yyyyyz_1[j];

                    tg_xyyz_xyyyyzz_0[j] = pb_x * tg_yyz_xyyyyzz_0[j] + fr * tg_yyz_xyyyyzz_1[j] + 0.5 * fl1_fxn * tg_yyz_yyyyzz_1[j];

                    tg_xyyz_xyyyzzz_0[j] = pb_x * tg_yyz_xyyyzzz_0[j] + fr * tg_yyz_xyyyzzz_1[j] + 0.5 * fl1_fxn * tg_yyz_yyyzzz_1[j];

                    tg_xyyz_xyyzzzz_0[j] = pb_x * tg_yyz_xyyzzzz_0[j] + fr * tg_yyz_xyyzzzz_1[j] + 0.5 * fl1_fxn * tg_yyz_yyzzzz_1[j];

                    tg_xyyz_xyzzzzz_0[j] = pb_x * tg_yyz_xyzzzzz_0[j] + fr * tg_yyz_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyz_yzzzzz_1[j];

                    tg_xyyz_xzzzzzz_0[j] = pb_x * tg_yyz_xzzzzzz_0[j] + fr * tg_yyz_xzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyz_zzzzzz_1[j];

                    tg_xyyz_yyyyyyy_0[j] = pb_x * tg_yyz_yyyyyyy_0[j] + fr * tg_yyz_yyyyyyy_1[j];

                    tg_xyyz_yyyyyyz_0[j] = pb_x * tg_yyz_yyyyyyz_0[j] + fr * tg_yyz_yyyyyyz_1[j];

                    tg_xyyz_yyyyyzz_0[j] = pb_x * tg_yyz_yyyyyzz_0[j] + fr * tg_yyz_yyyyyzz_1[j];

                    tg_xyyz_yyyyzzz_0[j] = pb_x * tg_yyz_yyyyzzz_0[j] + fr * tg_yyz_yyyyzzz_1[j];

                    tg_xyyz_yyyzzzz_0[j] = pb_x * tg_yyz_yyyzzzz_0[j] + fr * tg_yyz_yyyzzzz_1[j];

                    tg_xyyz_yyzzzzz_0[j] = pb_x * tg_yyz_yyzzzzz_0[j] + fr * tg_yyz_yyzzzzz_1[j];

                    tg_xyyz_yzzzzzz_0[j] = pb_x * tg_yyz_yzzzzzz_0[j] + fr * tg_yyz_yzzzzzz_1[j];

                    tg_xyyz_zzzzzzz_0[j] = pb_x * tg_yyz_zzzzzzz_0[j] + fr * tg_yyz_zzzzzzz_1[j];

                    tg_xyzz_xxxxxxx_0[j] = pb_x * tg_yzz_xxxxxxx_0[j] + fr * tg_yzz_xxxxxxx_1[j] + 3.5 * fl1_fxn * tg_yzz_xxxxxx_1[j];

                    tg_xyzz_xxxxxxy_0[j] = pb_x * tg_yzz_xxxxxxy_0[j] + fr * tg_yzz_xxxxxxy_1[j] + 3.0 * fl1_fxn * tg_yzz_xxxxxy_1[j];

                    tg_xyzz_xxxxxxz_0[j] = pb_x * tg_yzz_xxxxxxz_0[j] + fr * tg_yzz_xxxxxxz_1[j] + 3.0 * fl1_fxn * tg_yzz_xxxxxz_1[j];

                    tg_xyzz_xxxxxyy_0[j] = pb_x * tg_yzz_xxxxxyy_0[j] + fr * tg_yzz_xxxxxyy_1[j] + 2.5 * fl1_fxn * tg_yzz_xxxxyy_1[j];

                    tg_xyzz_xxxxxyz_0[j] = pb_x * tg_yzz_xxxxxyz_0[j] + fr * tg_yzz_xxxxxyz_1[j] + 2.5 * fl1_fxn * tg_yzz_xxxxyz_1[j];

                    tg_xyzz_xxxxxzz_0[j] = pb_x * tg_yzz_xxxxxzz_0[j] + fr * tg_yzz_xxxxxzz_1[j] + 2.5 * fl1_fxn * tg_yzz_xxxxzz_1[j];

                    tg_xyzz_xxxxyyy_0[j] = pb_x * tg_yzz_xxxxyyy_0[j] + fr * tg_yzz_xxxxyyy_1[j] + 2.0 * fl1_fxn * tg_yzz_xxxyyy_1[j];

                    tg_xyzz_xxxxyyz_0[j] = pb_x * tg_yzz_xxxxyyz_0[j] + fr * tg_yzz_xxxxyyz_1[j] + 2.0 * fl1_fxn * tg_yzz_xxxyyz_1[j];

                    tg_xyzz_xxxxyzz_0[j] = pb_x * tg_yzz_xxxxyzz_0[j] + fr * tg_yzz_xxxxyzz_1[j] + 2.0 * fl1_fxn * tg_yzz_xxxyzz_1[j];

                    tg_xyzz_xxxxzzz_0[j] = pb_x * tg_yzz_xxxxzzz_0[j] + fr * tg_yzz_xxxxzzz_1[j] + 2.0 * fl1_fxn * tg_yzz_xxxzzz_1[j];

                    tg_xyzz_xxxyyyy_0[j] = pb_x * tg_yzz_xxxyyyy_0[j] + fr * tg_yzz_xxxyyyy_1[j] + 1.5 * fl1_fxn * tg_yzz_xxyyyy_1[j];

                    tg_xyzz_xxxyyyz_0[j] = pb_x * tg_yzz_xxxyyyz_0[j] + fr * tg_yzz_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_yzz_xxyyyz_1[j];

                    tg_xyzz_xxxyyzz_0[j] = pb_x * tg_yzz_xxxyyzz_0[j] + fr * tg_yzz_xxxyyzz_1[j] + 1.5 * fl1_fxn * tg_yzz_xxyyzz_1[j];

                    tg_xyzz_xxxyzzz_0[j] = pb_x * tg_yzz_xxxyzzz_0[j] + fr * tg_yzz_xxxyzzz_1[j] + 1.5 * fl1_fxn * tg_yzz_xxyzzz_1[j];

                    tg_xyzz_xxxzzzz_0[j] = pb_x * tg_yzz_xxxzzzz_0[j] + fr * tg_yzz_xxxzzzz_1[j] + 1.5 * fl1_fxn * tg_yzz_xxzzzz_1[j];

                    tg_xyzz_xxyyyyy_0[j] = pb_x * tg_yzz_xxyyyyy_0[j] + fr * tg_yzz_xxyyyyy_1[j] + fl1_fxn * tg_yzz_xyyyyy_1[j];

                    tg_xyzz_xxyyyyz_0[j] = pb_x * tg_yzz_xxyyyyz_0[j] + fr * tg_yzz_xxyyyyz_1[j] + fl1_fxn * tg_yzz_xyyyyz_1[j];

                    tg_xyzz_xxyyyzz_0[j] = pb_x * tg_yzz_xxyyyzz_0[j] + fr * tg_yzz_xxyyyzz_1[j] + fl1_fxn * tg_yzz_xyyyzz_1[j];

                    tg_xyzz_xxyyzzz_0[j] = pb_x * tg_yzz_xxyyzzz_0[j] + fr * tg_yzz_xxyyzzz_1[j] + fl1_fxn * tg_yzz_xyyzzz_1[j];

                    tg_xyzz_xxyzzzz_0[j] = pb_x * tg_yzz_xxyzzzz_0[j] + fr * tg_yzz_xxyzzzz_1[j] + fl1_fxn * tg_yzz_xyzzzz_1[j];

                    tg_xyzz_xxzzzzz_0[j] = pb_x * tg_yzz_xxzzzzz_0[j] + fr * tg_yzz_xxzzzzz_1[j] + fl1_fxn * tg_yzz_xzzzzz_1[j];

                    tg_xyzz_xyyyyyy_0[j] = pb_x * tg_yzz_xyyyyyy_0[j] + fr * tg_yzz_xyyyyyy_1[j] + 0.5 * fl1_fxn * tg_yzz_yyyyyy_1[j];

                    tg_xyzz_xyyyyyz_0[j] = pb_x * tg_yzz_xyyyyyz_0[j] + fr * tg_yzz_xyyyyyz_1[j] + 0.5 * fl1_fxn * tg_yzz_yyyyyz_1[j];

                    tg_xyzz_xyyyyzz_0[j] = pb_x * tg_yzz_xyyyyzz_0[j] + fr * tg_yzz_xyyyyzz_1[j] + 0.5 * fl1_fxn * tg_yzz_yyyyzz_1[j];

                    tg_xyzz_xyyyzzz_0[j] = pb_x * tg_yzz_xyyyzzz_0[j] + fr * tg_yzz_xyyyzzz_1[j] + 0.5 * fl1_fxn * tg_yzz_yyyzzz_1[j];

                    tg_xyzz_xyyzzzz_0[j] = pb_x * tg_yzz_xyyzzzz_0[j] + fr * tg_yzz_xyyzzzz_1[j] + 0.5 * fl1_fxn * tg_yzz_yyzzzz_1[j];

                    tg_xyzz_xyzzzzz_0[j] = pb_x * tg_yzz_xyzzzzz_0[j] + fr * tg_yzz_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_yzz_yzzzzz_1[j];

                    tg_xyzz_xzzzzzz_0[j] = pb_x * tg_yzz_xzzzzzz_0[j] + fr * tg_yzz_xzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yzz_zzzzzz_1[j];

                    tg_xyzz_yyyyyyy_0[j] = pb_x * tg_yzz_yyyyyyy_0[j] + fr * tg_yzz_yyyyyyy_1[j];

                    tg_xyzz_yyyyyyz_0[j] = pb_x * tg_yzz_yyyyyyz_0[j] + fr * tg_yzz_yyyyyyz_1[j];

                    tg_xyzz_yyyyyzz_0[j] = pb_x * tg_yzz_yyyyyzz_0[j] + fr * tg_yzz_yyyyyzz_1[j];

                    tg_xyzz_yyyyzzz_0[j] = pb_x * tg_yzz_yyyyzzz_0[j] + fr * tg_yzz_yyyyzzz_1[j];

                    tg_xyzz_yyyzzzz_0[j] = pb_x * tg_yzz_yyyzzzz_0[j] + fr * tg_yzz_yyyzzzz_1[j];

                    tg_xyzz_yyzzzzz_0[j] = pb_x * tg_yzz_yyzzzzz_0[j] + fr * tg_yzz_yyzzzzz_1[j];

                    tg_xyzz_yzzzzzz_0[j] = pb_x * tg_yzz_yzzzzzz_0[j] + fr * tg_yzz_yzzzzzz_1[j];

                    tg_xyzz_zzzzzzz_0[j] = pb_x * tg_yzz_zzzzzzz_0[j] + fr * tg_yzz_zzzzzzz_1[j];

                    tg_xzzz_xxxxxxx_0[j] = pb_x * tg_zzz_xxxxxxx_0[j] + fr * tg_zzz_xxxxxxx_1[j] + 3.5 * fl1_fxn * tg_zzz_xxxxxx_1[j];

                    tg_xzzz_xxxxxxy_0[j] = pb_x * tg_zzz_xxxxxxy_0[j] + fr * tg_zzz_xxxxxxy_1[j] + 3.0 * fl1_fxn * tg_zzz_xxxxxy_1[j];

                    tg_xzzz_xxxxxxz_0[j] = pb_x * tg_zzz_xxxxxxz_0[j] + fr * tg_zzz_xxxxxxz_1[j] + 3.0 * fl1_fxn * tg_zzz_xxxxxz_1[j];

                    tg_xzzz_xxxxxyy_0[j] = pb_x * tg_zzz_xxxxxyy_0[j] + fr * tg_zzz_xxxxxyy_1[j] + 2.5 * fl1_fxn * tg_zzz_xxxxyy_1[j];

                    tg_xzzz_xxxxxyz_0[j] = pb_x * tg_zzz_xxxxxyz_0[j] + fr * tg_zzz_xxxxxyz_1[j] + 2.5 * fl1_fxn * tg_zzz_xxxxyz_1[j];

                    tg_xzzz_xxxxxzz_0[j] = pb_x * tg_zzz_xxxxxzz_0[j] + fr * tg_zzz_xxxxxzz_1[j] + 2.5 * fl1_fxn * tg_zzz_xxxxzz_1[j];

                    tg_xzzz_xxxxyyy_0[j] = pb_x * tg_zzz_xxxxyyy_0[j] + fr * tg_zzz_xxxxyyy_1[j] + 2.0 * fl1_fxn * tg_zzz_xxxyyy_1[j];

                    tg_xzzz_xxxxyyz_0[j] = pb_x * tg_zzz_xxxxyyz_0[j] + fr * tg_zzz_xxxxyyz_1[j] + 2.0 * fl1_fxn * tg_zzz_xxxyyz_1[j];

                    tg_xzzz_xxxxyzz_0[j] = pb_x * tg_zzz_xxxxyzz_0[j] + fr * tg_zzz_xxxxyzz_1[j] + 2.0 * fl1_fxn * tg_zzz_xxxyzz_1[j];

                    tg_xzzz_xxxxzzz_0[j] = pb_x * tg_zzz_xxxxzzz_0[j] + fr * tg_zzz_xxxxzzz_1[j] + 2.0 * fl1_fxn * tg_zzz_xxxzzz_1[j];

                    tg_xzzz_xxxyyyy_0[j] = pb_x * tg_zzz_xxxyyyy_0[j] + fr * tg_zzz_xxxyyyy_1[j] + 1.5 * fl1_fxn * tg_zzz_xxyyyy_1[j];

                    tg_xzzz_xxxyyyz_0[j] = pb_x * tg_zzz_xxxyyyz_0[j] + fr * tg_zzz_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_zzz_xxyyyz_1[j];

                    tg_xzzz_xxxyyzz_0[j] = pb_x * tg_zzz_xxxyyzz_0[j] + fr * tg_zzz_xxxyyzz_1[j] + 1.5 * fl1_fxn * tg_zzz_xxyyzz_1[j];

                    tg_xzzz_xxxyzzz_0[j] = pb_x * tg_zzz_xxxyzzz_0[j] + fr * tg_zzz_xxxyzzz_1[j] + 1.5 * fl1_fxn * tg_zzz_xxyzzz_1[j];

                    tg_xzzz_xxxzzzz_0[j] = pb_x * tg_zzz_xxxzzzz_0[j] + fr * tg_zzz_xxxzzzz_1[j] + 1.5 * fl1_fxn * tg_zzz_xxzzzz_1[j];

                    tg_xzzz_xxyyyyy_0[j] = pb_x * tg_zzz_xxyyyyy_0[j] + fr * tg_zzz_xxyyyyy_1[j] + fl1_fxn * tg_zzz_xyyyyy_1[j];

                    tg_xzzz_xxyyyyz_0[j] = pb_x * tg_zzz_xxyyyyz_0[j] + fr * tg_zzz_xxyyyyz_1[j] + fl1_fxn * tg_zzz_xyyyyz_1[j];

                    tg_xzzz_xxyyyzz_0[j] = pb_x * tg_zzz_xxyyyzz_0[j] + fr * tg_zzz_xxyyyzz_1[j] + fl1_fxn * tg_zzz_xyyyzz_1[j];

                    tg_xzzz_xxyyzzz_0[j] = pb_x * tg_zzz_xxyyzzz_0[j] + fr * tg_zzz_xxyyzzz_1[j] + fl1_fxn * tg_zzz_xyyzzz_1[j];

                    tg_xzzz_xxyzzzz_0[j] = pb_x * tg_zzz_xxyzzzz_0[j] + fr * tg_zzz_xxyzzzz_1[j] + fl1_fxn * tg_zzz_xyzzzz_1[j];

                    tg_xzzz_xxzzzzz_0[j] = pb_x * tg_zzz_xxzzzzz_0[j] + fr * tg_zzz_xxzzzzz_1[j] + fl1_fxn * tg_zzz_xzzzzz_1[j];

                    tg_xzzz_xyyyyyy_0[j] = pb_x * tg_zzz_xyyyyyy_0[j] + fr * tg_zzz_xyyyyyy_1[j] + 0.5 * fl1_fxn * tg_zzz_yyyyyy_1[j];

                    tg_xzzz_xyyyyyz_0[j] = pb_x * tg_zzz_xyyyyyz_0[j] + fr * tg_zzz_xyyyyyz_1[j] + 0.5 * fl1_fxn * tg_zzz_yyyyyz_1[j];

                    tg_xzzz_xyyyyzz_0[j] = pb_x * tg_zzz_xyyyyzz_0[j] + fr * tg_zzz_xyyyyzz_1[j] + 0.5 * fl1_fxn * tg_zzz_yyyyzz_1[j];

                    tg_xzzz_xyyyzzz_0[j] = pb_x * tg_zzz_xyyyzzz_0[j] + fr * tg_zzz_xyyyzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_yyyzzz_1[j];

                    tg_xzzz_xyyzzzz_0[j] = pb_x * tg_zzz_xyyzzzz_0[j] + fr * tg_zzz_xyyzzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_yyzzzz_1[j];

                    tg_xzzz_xyzzzzz_0[j] = pb_x * tg_zzz_xyzzzzz_0[j] + fr * tg_zzz_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_yzzzzz_1[j];

                    tg_xzzz_xzzzzzz_0[j] = pb_x * tg_zzz_xzzzzzz_0[j] + fr * tg_zzz_xzzzzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_zzzzzz_1[j];

                    tg_xzzz_yyyyyyy_0[j] = pb_x * tg_zzz_yyyyyyy_0[j] + fr * tg_zzz_yyyyyyy_1[j];

                    tg_xzzz_yyyyyyz_0[j] = pb_x * tg_zzz_yyyyyyz_0[j] + fr * tg_zzz_yyyyyyz_1[j];

                    tg_xzzz_yyyyyzz_0[j] = pb_x * tg_zzz_yyyyyzz_0[j] + fr * tg_zzz_yyyyyzz_1[j];

                    tg_xzzz_yyyyzzz_0[j] = pb_x * tg_zzz_yyyyzzz_0[j] + fr * tg_zzz_yyyyzzz_1[j];

                    tg_xzzz_yyyzzzz_0[j] = pb_x * tg_zzz_yyyzzzz_0[j] + fr * tg_zzz_yyyzzzz_1[j];

                    tg_xzzz_yyzzzzz_0[j] = pb_x * tg_zzz_yyzzzzz_0[j] + fr * tg_zzz_yyzzzzz_1[j];

                    tg_xzzz_yzzzzzz_0[j] = pb_x * tg_zzz_yzzzzzz_0[j] + fr * tg_zzz_yzzzzzz_1[j];

                    tg_xzzz_zzzzzzz_0[j] = pb_x * tg_zzz_zzzzzzz_0[j] + fr * tg_zzz_zzzzzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSK_360_450(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (360,450)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {4, -1, -1, -1},
                                             {7, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_7_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {7, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_7_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_7_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_7_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_2_7_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {7, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_2_7_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {7, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {6, -1, -1, -1}, 
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

                auto tg_yyy_xxxxxxx_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 216); 

                auto tg_yyy_xxxxxxy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 217); 

                auto tg_yyy_xxxxxxz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 218); 

                auto tg_yyy_xxxxxyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 219); 

                auto tg_yyy_xxxxxyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 220); 

                auto tg_yyy_xxxxxzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 221); 

                auto tg_yyy_xxxxyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 222); 

                auto tg_yyy_xxxxyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 223); 

                auto tg_yyy_xxxxyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 224); 

                auto tg_yyy_xxxxzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 225); 

                auto tg_yyy_xxxyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 226); 

                auto tg_yyy_xxxyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 227); 

                auto tg_yyy_xxxyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 228); 

                auto tg_yyy_xxxyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 229); 

                auto tg_yyy_xxxzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 230); 

                auto tg_yyy_xxyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 231); 

                auto tg_yyy_xxyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 232); 

                auto tg_yyy_xxyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 233); 

                auto tg_yyy_xxyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 234); 

                auto tg_yyy_xxyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 235); 

                auto tg_yyy_xxzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 236); 

                auto tg_yyy_xyyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 237); 

                auto tg_yyy_xyyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 238); 

                auto tg_yyy_xyyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 239); 

                auto tg_yyy_xyyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 240); 

                auto tg_yyy_xyyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 241); 

                auto tg_yyy_xyzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 242); 

                auto tg_yyy_xzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 243); 

                auto tg_yyy_yyyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 244); 

                auto tg_yyy_yyyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 245); 

                auto tg_yyy_yyyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 246); 

                auto tg_yyy_yyyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 247); 

                auto tg_yyy_yyyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 248); 

                auto tg_yyy_yyzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 249); 

                auto tg_yyy_yzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 250); 

                auto tg_yyy_zzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 251); 

                auto tg_yyz_xxxxxxx_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 252); 

                auto tg_yyz_xxxxxxy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 253); 

                auto tg_yyz_xxxxxxz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 254); 

                auto tg_yyz_xxxxxyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 255); 

                auto tg_yyz_xxxxxyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 256); 

                auto tg_yyz_xxxxxzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 257); 

                auto tg_yyz_xxxxyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 258); 

                auto tg_yyz_xxxxyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 259); 

                auto tg_yyz_xxxxyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 260); 

                auto tg_yyz_xxxxzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 261); 

                auto tg_yyz_xxxyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 262); 

                auto tg_yyz_xxxyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 263); 

                auto tg_yyz_xxxyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 264); 

                auto tg_yyz_xxxyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 265); 

                auto tg_yyz_xxxzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 266); 

                auto tg_yyz_xxyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 267); 

                auto tg_yyz_xxyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 268); 

                auto tg_yyz_xxyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 269); 

                auto tg_yyz_xxyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 270); 

                auto tg_yyz_xxyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 271); 

                auto tg_yyz_xxzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 272); 

                auto tg_yyz_xyyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 273); 

                auto tg_yyz_xyyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 274); 

                auto tg_yyz_xyyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 275); 

                auto tg_yyz_xyyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 276); 

                auto tg_yyz_xyyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 277); 

                auto tg_yyz_xyzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 278); 

                auto tg_yyz_xzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 279); 

                auto tg_yyz_yyyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 280); 

                auto tg_yyz_yyyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 281); 

                auto tg_yyz_yyyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 282); 

                auto tg_yyz_yyyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 283); 

                auto tg_yyz_yyyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 284); 

                auto tg_yyz_yyzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 285); 

                auto tg_yyz_yzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 286); 

                auto tg_yyz_zzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 287); 

                auto tg_yzz_xxxxxxx_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 288); 

                auto tg_yzz_xxxxxxy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 289); 

                auto tg_yzz_xxxxxxz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 290); 

                auto tg_yzz_xxxxxyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 291); 

                auto tg_yzz_xxxxxyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 292); 

                auto tg_yzz_xxxxxzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 293); 

                auto tg_yzz_xxxxyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 294); 

                auto tg_yzz_xxxxyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 295); 

                auto tg_yzz_xxxxyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 296); 

                auto tg_yzz_xxxxzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 297); 

                auto tg_yzz_xxxyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 298); 

                auto tg_yzz_xxxyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 299); 

                auto tg_yzz_xxxyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 300); 

                auto tg_yzz_xxxyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 301); 

                auto tg_yzz_xxxzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 302); 

                auto tg_yzz_xxyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 303); 

                auto tg_yzz_xxyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 304); 

                auto tg_yzz_xxyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 305); 

                auto tg_yyy_xxxxxxx_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 216); 

                auto tg_yyy_xxxxxxy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 217); 

                auto tg_yyy_xxxxxxz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 218); 

                auto tg_yyy_xxxxxyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 219); 

                auto tg_yyy_xxxxxyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 220); 

                auto tg_yyy_xxxxxzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 221); 

                auto tg_yyy_xxxxyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 222); 

                auto tg_yyy_xxxxyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 223); 

                auto tg_yyy_xxxxyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 224); 

                auto tg_yyy_xxxxzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 225); 

                auto tg_yyy_xxxyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 226); 

                auto tg_yyy_xxxyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 227); 

                auto tg_yyy_xxxyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 228); 

                auto tg_yyy_xxxyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 229); 

                auto tg_yyy_xxxzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 230); 

                auto tg_yyy_xxyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 231); 

                auto tg_yyy_xxyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 232); 

                auto tg_yyy_xxyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 233); 

                auto tg_yyy_xxyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 234); 

                auto tg_yyy_xxyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 235); 

                auto tg_yyy_xxzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 236); 

                auto tg_yyy_xyyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 237); 

                auto tg_yyy_xyyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 238); 

                auto tg_yyy_xyyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 239); 

                auto tg_yyy_xyyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 240); 

                auto tg_yyy_xyyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 241); 

                auto tg_yyy_xyzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 242); 

                auto tg_yyy_xzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 243); 

                auto tg_yyy_yyyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 244); 

                auto tg_yyy_yyyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 245); 

                auto tg_yyy_yyyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 246); 

                auto tg_yyy_yyyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 247); 

                auto tg_yyy_yyyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 248); 

                auto tg_yyy_yyzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 249); 

                auto tg_yyy_yzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 250); 

                auto tg_yyy_zzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 251); 

                auto tg_yyz_xxxxxxx_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 252); 

                auto tg_yyz_xxxxxxy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 253); 

                auto tg_yyz_xxxxxxz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 254); 

                auto tg_yyz_xxxxxyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 255); 

                auto tg_yyz_xxxxxyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 256); 

                auto tg_yyz_xxxxxzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 257); 

                auto tg_yyz_xxxxyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 258); 

                auto tg_yyz_xxxxyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 259); 

                auto tg_yyz_xxxxyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 260); 

                auto tg_yyz_xxxxzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 261); 

                auto tg_yyz_xxxyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 262); 

                auto tg_yyz_xxxyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 263); 

                auto tg_yyz_xxxyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 264); 

                auto tg_yyz_xxxyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 265); 

                auto tg_yyz_xxxzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 266); 

                auto tg_yyz_xxyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 267); 

                auto tg_yyz_xxyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 268); 

                auto tg_yyz_xxyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 269); 

                auto tg_yyz_xxyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 270); 

                auto tg_yyz_xxyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 271); 

                auto tg_yyz_xxzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 272); 

                auto tg_yyz_xyyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 273); 

                auto tg_yyz_xyyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 274); 

                auto tg_yyz_xyyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 275); 

                auto tg_yyz_xyyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 276); 

                auto tg_yyz_xyyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 277); 

                auto tg_yyz_xyzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 278); 

                auto tg_yyz_xzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 279); 

                auto tg_yyz_yyyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 280); 

                auto tg_yyz_yyyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 281); 

                auto tg_yyz_yyyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 282); 

                auto tg_yyz_yyyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 283); 

                auto tg_yyz_yyyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 284); 

                auto tg_yyz_yyzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 285); 

                auto tg_yyz_yzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 286); 

                auto tg_yyz_zzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 287); 

                auto tg_yzz_xxxxxxx_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 288); 

                auto tg_yzz_xxxxxxy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 289); 

                auto tg_yzz_xxxxxxz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 290); 

                auto tg_yzz_xxxxxyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 291); 

                auto tg_yzz_xxxxxyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 292); 

                auto tg_yzz_xxxxxzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 293); 

                auto tg_yzz_xxxxyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 294); 

                auto tg_yzz_xxxxyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 295); 

                auto tg_yzz_xxxxyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 296); 

                auto tg_yzz_xxxxzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 297); 

                auto tg_yzz_xxxyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 298); 

                auto tg_yzz_xxxyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 299); 

                auto tg_yzz_xxxyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 300); 

                auto tg_yzz_xxxyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 301); 

                auto tg_yzz_xxxzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 302); 

                auto tg_yzz_xxyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 303); 

                auto tg_yzz_xxyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 304); 

                auto tg_yzz_xxyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 305); 

                auto tg_yy_xxxxxxx_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 108); 

                auto tg_yy_xxxxxxy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 109); 

                auto tg_yy_xxxxxxz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 110); 

                auto tg_yy_xxxxxyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 111); 

                auto tg_yy_xxxxxyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 112); 

                auto tg_yy_xxxxxzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 113); 

                auto tg_yy_xxxxyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 114); 

                auto tg_yy_xxxxyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 115); 

                auto tg_yy_xxxxyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 116); 

                auto tg_yy_xxxxzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 117); 

                auto tg_yy_xxxyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 118); 

                auto tg_yy_xxxyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 119); 

                auto tg_yy_xxxyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 120); 

                auto tg_yy_xxxyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 121); 

                auto tg_yy_xxxzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 122); 

                auto tg_yy_xxyyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 123); 

                auto tg_yy_xxyyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 124); 

                auto tg_yy_xxyyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 125); 

                auto tg_yy_xxyyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 126); 

                auto tg_yy_xxyzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 127); 

                auto tg_yy_xxzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 128); 

                auto tg_yy_xyyyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 129); 

                auto tg_yy_xyyyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 130); 

                auto tg_yy_xyyyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 131); 

                auto tg_yy_xyyyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 132); 

                auto tg_yy_xyyzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 133); 

                auto tg_yy_xyzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 134); 

                auto tg_yy_xzzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 135); 

                auto tg_yy_yyyyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 136); 

                auto tg_yy_yyyyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 137); 

                auto tg_yy_yyyyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 138); 

                auto tg_yy_yyyyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 139); 

                auto tg_yy_yyyzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 140); 

                auto tg_yy_yyzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 141); 

                auto tg_yy_yzzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 142); 

                auto tg_yy_zzzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 143); 

                auto tg_yz_xxxxxxx_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 144); 

                auto tg_yz_xxxxxxy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 145); 

                auto tg_yz_xxxxxxz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 146); 

                auto tg_yz_xxxxxyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 147); 

                auto tg_yz_xxxxxyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 148); 

                auto tg_yz_xxxxxzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 149); 

                auto tg_yz_xxxxyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 150); 

                auto tg_yz_xxxxyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 151); 

                auto tg_yz_xxxxyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 152); 

                auto tg_yz_xxxxzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 153); 

                auto tg_yz_xxxyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 154); 

                auto tg_yz_xxxyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 155); 

                auto tg_yz_xxxyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 156); 

                auto tg_yz_xxxyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 157); 

                auto tg_yz_xxxzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 158); 

                auto tg_yz_xxyyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 159); 

                auto tg_yz_xxyyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 160); 

                auto tg_yz_xxyyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 161); 

                auto tg_yz_xxyyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 162); 

                auto tg_yz_xxyzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 163); 

                auto tg_yz_xxzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 164); 

                auto tg_yz_xyyyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 165); 

                auto tg_yz_xyyyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 166); 

                auto tg_yz_xyyyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 167); 

                auto tg_yz_xyyyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 168); 

                auto tg_yz_xyyzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 169); 

                auto tg_yz_xyzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 170); 

                auto tg_yz_xzzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 171); 

                auto tg_yz_yyyyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 172); 

                auto tg_yz_yyyyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 173); 

                auto tg_yz_yyyyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 174); 

                auto tg_yz_yyyyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 175); 

                auto tg_yz_yyyzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 176); 

                auto tg_yz_yyzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 177); 

                auto tg_yz_yzzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 178); 

                auto tg_yz_zzzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 179); 

                auto tg_zz_xxxxxxx_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 180); 

                auto tg_zz_xxxxxxy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 181); 

                auto tg_zz_xxxxxxz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 182); 

                auto tg_zz_xxxxxyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 183); 

                auto tg_zz_xxxxxyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 184); 

                auto tg_zz_xxxxxzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 185); 

                auto tg_zz_xxxxyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 186); 

                auto tg_zz_xxxxyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 187); 

                auto tg_zz_xxxxyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 188); 

                auto tg_zz_xxxxzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 189); 

                auto tg_zz_xxxyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 190); 

                auto tg_zz_xxxyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 191); 

                auto tg_zz_xxxyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 192); 

                auto tg_zz_xxxyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 193); 

                auto tg_zz_xxxzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 194); 

                auto tg_zz_xxyyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 195); 

                auto tg_zz_xxyyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 196); 

                auto tg_zz_xxyyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 197); 

                auto tg_yy_xxxxxxx_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 108); 

                auto tg_yy_xxxxxxy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 109); 

                auto tg_yy_xxxxxxz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 110); 

                auto tg_yy_xxxxxyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 111); 

                auto tg_yy_xxxxxyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 112); 

                auto tg_yy_xxxxxzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 113); 

                auto tg_yy_xxxxyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 114); 

                auto tg_yy_xxxxyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 115); 

                auto tg_yy_xxxxyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 116); 

                auto tg_yy_xxxxzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 117); 

                auto tg_yy_xxxyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 118); 

                auto tg_yy_xxxyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 119); 

                auto tg_yy_xxxyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 120); 

                auto tg_yy_xxxyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 121); 

                auto tg_yy_xxxzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 122); 

                auto tg_yy_xxyyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 123); 

                auto tg_yy_xxyyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 124); 

                auto tg_yy_xxyyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 125); 

                auto tg_yy_xxyyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 126); 

                auto tg_yy_xxyzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 127); 

                auto tg_yy_xxzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 128); 

                auto tg_yy_xyyyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 129); 

                auto tg_yy_xyyyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 130); 

                auto tg_yy_xyyyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 131); 

                auto tg_yy_xyyyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 132); 

                auto tg_yy_xyyzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 133); 

                auto tg_yy_xyzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 134); 

                auto tg_yy_xzzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 135); 

                auto tg_yy_yyyyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 136); 

                auto tg_yy_yyyyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 137); 

                auto tg_yy_yyyyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 138); 

                auto tg_yy_yyyyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 139); 

                auto tg_yy_yyyzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 140); 

                auto tg_yy_yyzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 141); 

                auto tg_yy_yzzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 142); 

                auto tg_yy_zzzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 143); 

                auto tg_yz_xxxxxxx_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 144); 

                auto tg_yz_xxxxxxy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 145); 

                auto tg_yz_xxxxxxz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 146); 

                auto tg_yz_xxxxxyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 147); 

                auto tg_yz_xxxxxyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 148); 

                auto tg_yz_xxxxxzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 149); 

                auto tg_yz_xxxxyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 150); 

                auto tg_yz_xxxxyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 151); 

                auto tg_yz_xxxxyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 152); 

                auto tg_yz_xxxxzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 153); 

                auto tg_yz_xxxyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 154); 

                auto tg_yz_xxxyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 155); 

                auto tg_yz_xxxyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 156); 

                auto tg_yz_xxxyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 157); 

                auto tg_yz_xxxzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 158); 

                auto tg_yz_xxyyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 159); 

                auto tg_yz_xxyyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 160); 

                auto tg_yz_xxyyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 161); 

                auto tg_yz_xxyyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 162); 

                auto tg_yz_xxyzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 163); 

                auto tg_yz_xxzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 164); 

                auto tg_yz_xyyyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 165); 

                auto tg_yz_xyyyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 166); 

                auto tg_yz_xyyyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 167); 

                auto tg_yz_xyyyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 168); 

                auto tg_yz_xyyzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 169); 

                auto tg_yz_xyzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 170); 

                auto tg_yz_xzzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 171); 

                auto tg_yz_yyyyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 172); 

                auto tg_yz_yyyyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 173); 

                auto tg_yz_yyyyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 174); 

                auto tg_yz_yyyyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 175); 

                auto tg_yz_yyyzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 176); 

                auto tg_yz_yyzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 177); 

                auto tg_yz_yzzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 178); 

                auto tg_yz_zzzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 179); 

                auto tg_zz_xxxxxxx_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 180); 

                auto tg_zz_xxxxxxy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 181); 

                auto tg_zz_xxxxxxz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 182); 

                auto tg_zz_xxxxxyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 183); 

                auto tg_zz_xxxxxyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 184); 

                auto tg_zz_xxxxxzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 185); 

                auto tg_zz_xxxxyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 186); 

                auto tg_zz_xxxxyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 187); 

                auto tg_zz_xxxxyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 188); 

                auto tg_zz_xxxxzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 189); 

                auto tg_zz_xxxyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 190); 

                auto tg_zz_xxxyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 191); 

                auto tg_zz_xxxyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 192); 

                auto tg_zz_xxxyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 193); 

                auto tg_zz_xxxzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 194); 

                auto tg_zz_xxyyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 195); 

                auto tg_zz_xxyyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 196); 

                auto tg_zz_xxyyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 197); 

                auto tg_yyy_xxxxxx_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 168); 

                auto tg_yyy_xxxxxy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 169); 

                auto tg_yyy_xxxxxz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 170); 

                auto tg_yyy_xxxxyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 171); 

                auto tg_yyy_xxxxyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 172); 

                auto tg_yyy_xxxxzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 173); 

                auto tg_yyy_xxxyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 174); 

                auto tg_yyy_xxxyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 175); 

                auto tg_yyy_xxxyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 176); 

                auto tg_yyy_xxxzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 177); 

                auto tg_yyy_xxyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 178); 

                auto tg_yyy_xxyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 179); 

                auto tg_yyy_xxyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 180); 

                auto tg_yyy_xxyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 181); 

                auto tg_yyy_xxzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 182); 

                auto tg_yyy_xyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 183); 

                auto tg_yyy_xyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 184); 

                auto tg_yyy_xyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 185); 

                auto tg_yyy_xyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 186); 

                auto tg_yyy_xyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 187); 

                auto tg_yyy_xzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 188); 

                auto tg_yyy_yyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 189); 

                auto tg_yyy_yyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 190); 

                auto tg_yyy_yyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 191); 

                auto tg_yyy_yyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 192); 

                auto tg_yyy_yyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 193); 

                auto tg_yyy_yzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 194); 

                auto tg_yyy_zzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 195); 

                auto tg_yyz_xxxxxx_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 196); 

                auto tg_yyz_xxxxxy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 197); 

                auto tg_yyz_xxxxxz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 198); 

                auto tg_yyz_xxxxyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 199); 

                auto tg_yyz_xxxxyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 200); 

                auto tg_yyz_xxxxzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 201); 

                auto tg_yyz_xxxyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 202); 

                auto tg_yyz_xxxyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 203); 

                auto tg_yyz_xxxyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 204); 

                auto tg_yyz_xxxzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 205); 

                auto tg_yyz_xxyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 206); 

                auto tg_yyz_xxyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 207); 

                auto tg_yyz_xxyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 208); 

                auto tg_yyz_xxyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 209); 

                auto tg_yyz_xxzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 210); 

                auto tg_yyz_xyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 211); 

                auto tg_yyz_xyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 212); 

                auto tg_yyz_xyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 213); 

                auto tg_yyz_xyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 214); 

                auto tg_yyz_xyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 215); 

                auto tg_yyz_xzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 216); 

                auto tg_yyz_yyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 217); 

                auto tg_yyz_yyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 218); 

                auto tg_yyz_yyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 219); 

                auto tg_yyz_yyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 220); 

                auto tg_yyz_yyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 221); 

                auto tg_yyz_yzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 222); 

                auto tg_yyz_zzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 223); 

                auto tg_yzz_xxxxxx_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 224); 

                auto tg_yzz_xxxxxy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 225); 

                auto tg_yzz_xxxxxz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 226); 

                auto tg_yzz_xxxxyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 227); 

                auto tg_yzz_xxxxyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 228); 

                auto tg_yzz_xxxxzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 229); 

                auto tg_yzz_xxxyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 230); 

                auto tg_yzz_xxxyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 231); 

                auto tg_yzz_xxxyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 232); 

                auto tg_yzz_xxxzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 233); 

                auto tg_yzz_xxyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 234); 

                auto tg_yzz_xxyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 235); 

                auto tg_yzz_xxyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 236); 

                // set up pointers to integrals

                auto tg_yyyy_xxxxxxx_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 360); 

                auto tg_yyyy_xxxxxxy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 361); 

                auto tg_yyyy_xxxxxxz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 362); 

                auto tg_yyyy_xxxxxyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 363); 

                auto tg_yyyy_xxxxxyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 364); 

                auto tg_yyyy_xxxxxzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 365); 

                auto tg_yyyy_xxxxyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 366); 

                auto tg_yyyy_xxxxyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 367); 

                auto tg_yyyy_xxxxyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 368); 

                auto tg_yyyy_xxxxzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 369); 

                auto tg_yyyy_xxxyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 370); 

                auto tg_yyyy_xxxyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 371); 

                auto tg_yyyy_xxxyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 372); 

                auto tg_yyyy_xxxyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 373); 

                auto tg_yyyy_xxxzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 374); 

                auto tg_yyyy_xxyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 375); 

                auto tg_yyyy_xxyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 376); 

                auto tg_yyyy_xxyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 377); 

                auto tg_yyyy_xxyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 378); 

                auto tg_yyyy_xxyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 379); 

                auto tg_yyyy_xxzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 380); 

                auto tg_yyyy_xyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 381); 

                auto tg_yyyy_xyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 382); 

                auto tg_yyyy_xyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 383); 

                auto tg_yyyy_xyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 384); 

                auto tg_yyyy_xyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 385); 

                auto tg_yyyy_xyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 386); 

                auto tg_yyyy_xzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 387); 

                auto tg_yyyy_yyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 388); 

                auto tg_yyyy_yyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 389); 

                auto tg_yyyy_yyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 390); 

                auto tg_yyyy_yyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 391); 

                auto tg_yyyy_yyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 392); 

                auto tg_yyyy_yyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 393); 

                auto tg_yyyy_yzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 394); 

                auto tg_yyyy_zzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 395); 

                auto tg_yyyz_xxxxxxx_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 396); 

                auto tg_yyyz_xxxxxxy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 397); 

                auto tg_yyyz_xxxxxxz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 398); 

                auto tg_yyyz_xxxxxyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 399); 

                auto tg_yyyz_xxxxxyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 400); 

                auto tg_yyyz_xxxxxzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 401); 

                auto tg_yyyz_xxxxyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 402); 

                auto tg_yyyz_xxxxyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 403); 

                auto tg_yyyz_xxxxyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 404); 

                auto tg_yyyz_xxxxzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 405); 

                auto tg_yyyz_xxxyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 406); 

                auto tg_yyyz_xxxyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 407); 

                auto tg_yyyz_xxxyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 408); 

                auto tg_yyyz_xxxyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 409); 

                auto tg_yyyz_xxxzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 410); 

                auto tg_yyyz_xxyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 411); 

                auto tg_yyyz_xxyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 412); 

                auto tg_yyyz_xxyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 413); 

                auto tg_yyyz_xxyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 414); 

                auto tg_yyyz_xxyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 415); 

                auto tg_yyyz_xxzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 416); 

                auto tg_yyyz_xyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 417); 

                auto tg_yyyz_xyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 418); 

                auto tg_yyyz_xyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 419); 

                auto tg_yyyz_xyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 420); 

                auto tg_yyyz_xyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 421); 

                auto tg_yyyz_xyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 422); 

                auto tg_yyyz_xzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 423); 

                auto tg_yyyz_yyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 424); 

                auto tg_yyyz_yyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 425); 

                auto tg_yyyz_yyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 426); 

                auto tg_yyyz_yyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 427); 

                auto tg_yyyz_yyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 428); 

                auto tg_yyyz_yyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 429); 

                auto tg_yyyz_yzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 430); 

                auto tg_yyyz_zzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 431); 

                auto tg_yyzz_xxxxxxx_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 432); 

                auto tg_yyzz_xxxxxxy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 433); 

                auto tg_yyzz_xxxxxxz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 434); 

                auto tg_yyzz_xxxxxyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 435); 

                auto tg_yyzz_xxxxxyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 436); 

                auto tg_yyzz_xxxxxzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 437); 

                auto tg_yyzz_xxxxyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 438); 

                auto tg_yyzz_xxxxyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 439); 

                auto tg_yyzz_xxxxyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 440); 

                auto tg_yyzz_xxxxzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 441); 

                auto tg_yyzz_xxxyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 442); 

                auto tg_yyzz_xxxyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 443); 

                auto tg_yyzz_xxxyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 444); 

                auto tg_yyzz_xxxyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 445); 

                auto tg_yyzz_xxxzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 446); 

                auto tg_yyzz_xxyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 447); 

                auto tg_yyzz_xxyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 448); 

                auto tg_yyzz_xxyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 449); 

                // Batch of Integrals (360,450)

                #pragma omp simd aligned(fxn, fza, tg_yy_xxxxxxx_0, tg_yy_xxxxxxx_1, tg_yy_xxxxxxy_0, \
                                         tg_yy_xxxxxxy_1, tg_yy_xxxxxxz_0, tg_yy_xxxxxxz_1, tg_yy_xxxxxyy_0, tg_yy_xxxxxyy_1, \
                                         tg_yy_xxxxxyz_0, tg_yy_xxxxxyz_1, tg_yy_xxxxxzz_0, tg_yy_xxxxxzz_1, tg_yy_xxxxyyy_0, \
                                         tg_yy_xxxxyyy_1, tg_yy_xxxxyyz_0, tg_yy_xxxxyyz_1, tg_yy_xxxxyzz_0, tg_yy_xxxxyzz_1, \
                                         tg_yy_xxxxzzz_0, tg_yy_xxxxzzz_1, tg_yy_xxxyyyy_0, tg_yy_xxxyyyy_1, tg_yy_xxxyyyz_0, \
                                         tg_yy_xxxyyyz_1, tg_yy_xxxyyzz_0, tg_yy_xxxyyzz_1, tg_yy_xxxyzzz_0, tg_yy_xxxyzzz_1, \
                                         tg_yy_xxxzzzz_0, tg_yy_xxxzzzz_1, tg_yy_xxyyyyy_0, tg_yy_xxyyyyy_1, tg_yy_xxyyyyz_0, \
                                         tg_yy_xxyyyyz_1, tg_yy_xxyyyzz_0, tg_yy_xxyyyzz_1, tg_yy_xxyyzzz_0, tg_yy_xxyyzzz_1, \
                                         tg_yy_xxyzzzz_0, tg_yy_xxyzzzz_1, tg_yy_xxzzzzz_0, tg_yy_xxzzzzz_1, tg_yy_xyyyyyy_0, \
                                         tg_yy_xyyyyyy_1, tg_yy_xyyyyyz_0, tg_yy_xyyyyyz_1, tg_yy_xyyyyzz_0, tg_yy_xyyyyzz_1, \
                                         tg_yy_xyyyzzz_0, tg_yy_xyyyzzz_1, tg_yy_xyyzzzz_0, tg_yy_xyyzzzz_1, tg_yy_xyzzzzz_0, \
                                         tg_yy_xyzzzzz_1, tg_yy_xzzzzzz_0, tg_yy_xzzzzzz_1, tg_yy_yyyyyyy_0, tg_yy_yyyyyyy_1, \
                                         tg_yy_yyyyyyz_0, tg_yy_yyyyyyz_1, tg_yy_yyyyyzz_0, tg_yy_yyyyyzz_1, tg_yy_yyyyzzz_0, \
                                         tg_yy_yyyyzzz_1, tg_yy_yyyzzzz_0, tg_yy_yyyzzzz_1, tg_yy_yyzzzzz_0, tg_yy_yyzzzzz_1, \
                                         tg_yy_yzzzzzz_0, tg_yy_yzzzzzz_1, tg_yy_zzzzzzz_0, tg_yy_zzzzzzz_1, tg_yyy_xxxxxx_1, \
                                         tg_yyy_xxxxxxx_0, tg_yyy_xxxxxxx_1, tg_yyy_xxxxxxy_0, tg_yyy_xxxxxxy_1, \
                                         tg_yyy_xxxxxxz_0, tg_yyy_xxxxxxz_1, tg_yyy_xxxxxy_1, tg_yyy_xxxxxyy_0, \
                                         tg_yyy_xxxxxyy_1, tg_yyy_xxxxxyz_0, tg_yyy_xxxxxyz_1, tg_yyy_xxxxxz_1, \
                                         tg_yyy_xxxxxzz_0, tg_yyy_xxxxxzz_1, tg_yyy_xxxxyy_1, tg_yyy_xxxxyyy_0, \
                                         tg_yyy_xxxxyyy_1, tg_yyy_xxxxyyz_0, tg_yyy_xxxxyyz_1, tg_yyy_xxxxyz_1, \
                                         tg_yyy_xxxxyzz_0, tg_yyy_xxxxyzz_1, tg_yyy_xxxxzz_1, tg_yyy_xxxxzzz_0, \
                                         tg_yyy_xxxxzzz_1, tg_yyy_xxxyyy_1, tg_yyy_xxxyyyy_0, tg_yyy_xxxyyyy_1, \
                                         tg_yyy_xxxyyyz_0, tg_yyy_xxxyyyz_1, tg_yyy_xxxyyz_1, tg_yyy_xxxyyzz_0, \
                                         tg_yyy_xxxyyzz_1, tg_yyy_xxxyzz_1, tg_yyy_xxxyzzz_0, tg_yyy_xxxyzzz_1, \
                                         tg_yyy_xxxzzz_1, tg_yyy_xxxzzzz_0, tg_yyy_xxxzzzz_1, tg_yyy_xxyyyy_1, \
                                         tg_yyy_xxyyyyy_0, tg_yyy_xxyyyyy_1, tg_yyy_xxyyyyz_0, tg_yyy_xxyyyyz_1, \
                                         tg_yyy_xxyyyz_1, tg_yyy_xxyyyzz_0, tg_yyy_xxyyyzz_1, tg_yyy_xxyyzz_1, \
                                         tg_yyy_xxyyzzz_0, tg_yyy_xxyyzzz_1, tg_yyy_xxyzzz_1, tg_yyy_xxyzzzz_0, \
                                         tg_yyy_xxyzzzz_1, tg_yyy_xxzzzz_1, tg_yyy_xxzzzzz_0, tg_yyy_xxzzzzz_1, \
                                         tg_yyy_xyyyyy_1, tg_yyy_xyyyyyy_0, tg_yyy_xyyyyyy_1, tg_yyy_xyyyyyz_0, \
                                         tg_yyy_xyyyyyz_1, tg_yyy_xyyyyz_1, tg_yyy_xyyyyzz_0, tg_yyy_xyyyyzz_1, \
                                         tg_yyy_xyyyzz_1, tg_yyy_xyyyzzz_0, tg_yyy_xyyyzzz_1, tg_yyy_xyyzzz_1, \
                                         tg_yyy_xyyzzzz_0, tg_yyy_xyyzzzz_1, tg_yyy_xyzzzz_1, tg_yyy_xyzzzzz_0, \
                                         tg_yyy_xyzzzzz_1, tg_yyy_xzzzzz_1, tg_yyy_xzzzzzz_0, tg_yyy_xzzzzzz_1, \
                                         tg_yyy_yyyyyy_1, tg_yyy_yyyyyyy_0, tg_yyy_yyyyyyy_1, tg_yyy_yyyyyyz_0, \
                                         tg_yyy_yyyyyyz_1, tg_yyy_yyyyyz_1, tg_yyy_yyyyyzz_0, tg_yyy_yyyyyzz_1, \
                                         tg_yyy_yyyyzz_1, tg_yyy_yyyyzzz_0, tg_yyy_yyyyzzz_1, tg_yyy_yyyzzz_1, \
                                         tg_yyy_yyyzzzz_0, tg_yyy_yyyzzzz_1, tg_yyy_yyzzzz_1, tg_yyy_yyzzzzz_0, \
                                         tg_yyy_yyzzzzz_1, tg_yyy_yzzzzz_1, tg_yyy_yzzzzzz_0, tg_yyy_yzzzzzz_1, \
                                         tg_yyy_zzzzzz_1, tg_yyy_zzzzzzz_0, tg_yyy_zzzzzzz_1, tg_yyyy_xxxxxxx_0, \
                                         tg_yyyy_xxxxxxy_0, tg_yyyy_xxxxxxz_0, tg_yyyy_xxxxxyy_0, tg_yyyy_xxxxxyz_0, \
                                         tg_yyyy_xxxxxzz_0, tg_yyyy_xxxxyyy_0, tg_yyyy_xxxxyyz_0, tg_yyyy_xxxxyzz_0, \
                                         tg_yyyy_xxxxzzz_0, tg_yyyy_xxxyyyy_0, tg_yyyy_xxxyyyz_0, tg_yyyy_xxxyyzz_0, \
                                         tg_yyyy_xxxyzzz_0, tg_yyyy_xxxzzzz_0, tg_yyyy_xxyyyyy_0, tg_yyyy_xxyyyyz_0, \
                                         tg_yyyy_xxyyyzz_0, tg_yyyy_xxyyzzz_0, tg_yyyy_xxyzzzz_0, tg_yyyy_xxzzzzz_0, \
                                         tg_yyyy_xyyyyyy_0, tg_yyyy_xyyyyyz_0, tg_yyyy_xyyyyzz_0, tg_yyyy_xyyyzzz_0, \
                                         tg_yyyy_xyyzzzz_0, tg_yyyy_xyzzzzz_0, tg_yyyy_xzzzzzz_0, tg_yyyy_yyyyyyy_0, \
                                         tg_yyyy_yyyyyyz_0, tg_yyyy_yyyyyzz_0, tg_yyyy_yyyyzzz_0, tg_yyyy_yyyzzzz_0, \
                                         tg_yyyy_yyzzzzz_0, tg_yyyy_yzzzzzz_0, tg_yyyy_zzzzzzz_0, tg_yyyz_xxxxxxx_0, \
                                         tg_yyyz_xxxxxxy_0, tg_yyyz_xxxxxxz_0, tg_yyyz_xxxxxyy_0, tg_yyyz_xxxxxyz_0, \
                                         tg_yyyz_xxxxxzz_0, tg_yyyz_xxxxyyy_0, tg_yyyz_xxxxyyz_0, tg_yyyz_xxxxyzz_0, \
                                         tg_yyyz_xxxxzzz_0, tg_yyyz_xxxyyyy_0, tg_yyyz_xxxyyyz_0, tg_yyyz_xxxyyzz_0, \
                                         tg_yyyz_xxxyzzz_0, tg_yyyz_xxxzzzz_0, tg_yyyz_xxyyyyy_0, tg_yyyz_xxyyyyz_0, \
                                         tg_yyyz_xxyyyzz_0, tg_yyyz_xxyyzzz_0, tg_yyyz_xxyzzzz_0, tg_yyyz_xxzzzzz_0, \
                                         tg_yyyz_xyyyyyy_0, tg_yyyz_xyyyyyz_0, tg_yyyz_xyyyyzz_0, tg_yyyz_xyyyzzz_0, \
                                         tg_yyyz_xyyzzzz_0, tg_yyyz_xyzzzzz_0, tg_yyyz_xzzzzzz_0, tg_yyyz_yyyyyyy_0, \
                                         tg_yyyz_yyyyyyz_0, tg_yyyz_yyyyyzz_0, tg_yyyz_yyyyzzz_0, tg_yyyz_yyyzzzz_0, \
                                         tg_yyyz_yyzzzzz_0, tg_yyyz_yzzzzzz_0, tg_yyyz_zzzzzzz_0, tg_yyz_xxxxxx_1, \
                                         tg_yyz_xxxxxxx_0, tg_yyz_xxxxxxx_1, tg_yyz_xxxxxxy_0, tg_yyz_xxxxxxy_1, \
                                         tg_yyz_xxxxxxz_0, tg_yyz_xxxxxxz_1, tg_yyz_xxxxxy_1, tg_yyz_xxxxxyy_0, \
                                         tg_yyz_xxxxxyy_1, tg_yyz_xxxxxyz_0, tg_yyz_xxxxxyz_1, tg_yyz_xxxxxz_1, \
                                         tg_yyz_xxxxxzz_0, tg_yyz_xxxxxzz_1, tg_yyz_xxxxyy_1, tg_yyz_xxxxyyy_0, \
                                         tg_yyz_xxxxyyy_1, tg_yyz_xxxxyyz_0, tg_yyz_xxxxyyz_1, tg_yyz_xxxxyz_1, \
                                         tg_yyz_xxxxyzz_0, tg_yyz_xxxxyzz_1, tg_yyz_xxxxzz_1, tg_yyz_xxxxzzz_0, \
                                         tg_yyz_xxxxzzz_1, tg_yyz_xxxyyy_1, tg_yyz_xxxyyyy_0, tg_yyz_xxxyyyy_1, \
                                         tg_yyz_xxxyyyz_0, tg_yyz_xxxyyyz_1, tg_yyz_xxxyyz_1, tg_yyz_xxxyyzz_0, \
                                         tg_yyz_xxxyyzz_1, tg_yyz_xxxyzz_1, tg_yyz_xxxyzzz_0, tg_yyz_xxxyzzz_1, \
                                         tg_yyz_xxxzzz_1, tg_yyz_xxxzzzz_0, tg_yyz_xxxzzzz_1, tg_yyz_xxyyyy_1, \
                                         tg_yyz_xxyyyyy_0, tg_yyz_xxyyyyy_1, tg_yyz_xxyyyyz_0, tg_yyz_xxyyyyz_1, \
                                         tg_yyz_xxyyyz_1, tg_yyz_xxyyyzz_0, tg_yyz_xxyyyzz_1, tg_yyz_xxyyzz_1, \
                                         tg_yyz_xxyyzzz_0, tg_yyz_xxyyzzz_1, tg_yyz_xxyzzz_1, tg_yyz_xxyzzzz_0, \
                                         tg_yyz_xxyzzzz_1, tg_yyz_xxzzzz_1, tg_yyz_xxzzzzz_0, tg_yyz_xxzzzzz_1, \
                                         tg_yyz_xyyyyy_1, tg_yyz_xyyyyyy_0, tg_yyz_xyyyyyy_1, tg_yyz_xyyyyyz_0, \
                                         tg_yyz_xyyyyyz_1, tg_yyz_xyyyyz_1, tg_yyz_xyyyyzz_0, tg_yyz_xyyyyzz_1, \
                                         tg_yyz_xyyyzz_1, tg_yyz_xyyyzzz_0, tg_yyz_xyyyzzz_1, tg_yyz_xyyzzz_1, \
                                         tg_yyz_xyyzzzz_0, tg_yyz_xyyzzzz_1, tg_yyz_xyzzzz_1, tg_yyz_xyzzzzz_0, \
                                         tg_yyz_xyzzzzz_1, tg_yyz_xzzzzz_1, tg_yyz_xzzzzzz_0, tg_yyz_xzzzzzz_1, \
                                         tg_yyz_yyyyyy_1, tg_yyz_yyyyyyy_0, tg_yyz_yyyyyyy_1, tg_yyz_yyyyyyz_0, \
                                         tg_yyz_yyyyyyz_1, tg_yyz_yyyyyz_1, tg_yyz_yyyyyzz_0, tg_yyz_yyyyyzz_1, \
                                         tg_yyz_yyyyzz_1, tg_yyz_yyyyzzz_0, tg_yyz_yyyyzzz_1, tg_yyz_yyyzzz_1, \
                                         tg_yyz_yyyzzzz_0, tg_yyz_yyyzzzz_1, tg_yyz_yyzzzz_1, tg_yyz_yyzzzzz_0, \
                                         tg_yyz_yyzzzzz_1, tg_yyz_yzzzzz_1, tg_yyz_yzzzzzz_0, tg_yyz_yzzzzzz_1, \
                                         tg_yyz_zzzzzz_1, tg_yyz_zzzzzzz_0, tg_yyz_zzzzzzz_1, tg_yyzz_xxxxxxx_0, \
                                         tg_yyzz_xxxxxxy_0, tg_yyzz_xxxxxxz_0, tg_yyzz_xxxxxyy_0, tg_yyzz_xxxxxyz_0, \
                                         tg_yyzz_xxxxxzz_0, tg_yyzz_xxxxyyy_0, tg_yyzz_xxxxyyz_0, tg_yyzz_xxxxyzz_0, \
                                         tg_yyzz_xxxxzzz_0, tg_yyzz_xxxyyyy_0, tg_yyzz_xxxyyyz_0, tg_yyzz_xxxyyzz_0, \
                                         tg_yyzz_xxxyzzz_0, tg_yyzz_xxxzzzz_0, tg_yyzz_xxyyyyy_0, tg_yyzz_xxyyyyz_0, \
                                         tg_yyzz_xxyyyzz_0, tg_yz_xxxxxxx_0, tg_yz_xxxxxxx_1, tg_yz_xxxxxxy_0, tg_yz_xxxxxxy_1, \
                                         tg_yz_xxxxxxz_0, tg_yz_xxxxxxz_1, tg_yz_xxxxxyy_0, tg_yz_xxxxxyy_1, tg_yz_xxxxxyz_0, \
                                         tg_yz_xxxxxyz_1, tg_yz_xxxxxzz_0, tg_yz_xxxxxzz_1, tg_yz_xxxxyyy_0, tg_yz_xxxxyyy_1, \
                                         tg_yz_xxxxyyz_0, tg_yz_xxxxyyz_1, tg_yz_xxxxyzz_0, tg_yz_xxxxyzz_1, tg_yz_xxxxzzz_0, \
                                         tg_yz_xxxxzzz_1, tg_yz_xxxyyyy_0, tg_yz_xxxyyyy_1, tg_yz_xxxyyyz_0, tg_yz_xxxyyyz_1, \
                                         tg_yz_xxxyyzz_0, tg_yz_xxxyyzz_1, tg_yz_xxxyzzz_0, tg_yz_xxxyzzz_1, tg_yz_xxxzzzz_0, \
                                         tg_yz_xxxzzzz_1, tg_yz_xxyyyyy_0, tg_yz_xxyyyyy_1, tg_yz_xxyyyyz_0, tg_yz_xxyyyyz_1, \
                                         tg_yz_xxyyyzz_0, tg_yz_xxyyyzz_1, tg_yz_xxyyzzz_0, tg_yz_xxyyzzz_1, tg_yz_xxyzzzz_0, \
                                         tg_yz_xxyzzzz_1, tg_yz_xxzzzzz_0, tg_yz_xxzzzzz_1, tg_yz_xyyyyyy_0, tg_yz_xyyyyyy_1, \
                                         tg_yz_xyyyyyz_0, tg_yz_xyyyyyz_1, tg_yz_xyyyyzz_0, tg_yz_xyyyyzz_1, tg_yz_xyyyzzz_0, \
                                         tg_yz_xyyyzzz_1, tg_yz_xyyzzzz_0, tg_yz_xyyzzzz_1, tg_yz_xyzzzzz_0, tg_yz_xyzzzzz_1, \
                                         tg_yz_xzzzzzz_0, tg_yz_xzzzzzz_1, tg_yz_yyyyyyy_0, tg_yz_yyyyyyy_1, tg_yz_yyyyyyz_0, \
                                         tg_yz_yyyyyyz_1, tg_yz_yyyyyzz_0, tg_yz_yyyyyzz_1, tg_yz_yyyyzzz_0, tg_yz_yyyyzzz_1, \
                                         tg_yz_yyyzzzz_0, tg_yz_yyyzzzz_1, tg_yz_yyzzzzz_0, tg_yz_yyzzzzz_1, tg_yz_yzzzzzz_0, \
                                         tg_yz_yzzzzzz_1, tg_yz_zzzzzzz_0, tg_yz_zzzzzzz_1, tg_yzz_xxxxxx_1, \
                                         tg_yzz_xxxxxxx_0, tg_yzz_xxxxxxx_1, tg_yzz_xxxxxxy_0, tg_yzz_xxxxxxy_1, \
                                         tg_yzz_xxxxxxz_0, tg_yzz_xxxxxxz_1, tg_yzz_xxxxxy_1, tg_yzz_xxxxxyy_0, \
                                         tg_yzz_xxxxxyy_1, tg_yzz_xxxxxyz_0, tg_yzz_xxxxxyz_1, tg_yzz_xxxxxz_1, \
                                         tg_yzz_xxxxxzz_0, tg_yzz_xxxxxzz_1, tg_yzz_xxxxyy_1, tg_yzz_xxxxyyy_0, \
                                         tg_yzz_xxxxyyy_1, tg_yzz_xxxxyyz_0, tg_yzz_xxxxyyz_1, tg_yzz_xxxxyz_1, \
                                         tg_yzz_xxxxyzz_0, tg_yzz_xxxxyzz_1, tg_yzz_xxxxzz_1, tg_yzz_xxxxzzz_0, \
                                         tg_yzz_xxxxzzz_1, tg_yzz_xxxyyy_1, tg_yzz_xxxyyyy_0, tg_yzz_xxxyyyy_1, \
                                         tg_yzz_xxxyyyz_0, tg_yzz_xxxyyyz_1, tg_yzz_xxxyyz_1, tg_yzz_xxxyyzz_0, \
                                         tg_yzz_xxxyyzz_1, tg_yzz_xxxyzz_1, tg_yzz_xxxyzzz_0, tg_yzz_xxxyzzz_1, \
                                         tg_yzz_xxxzzz_1, tg_yzz_xxxzzzz_0, tg_yzz_xxxzzzz_1, tg_yzz_xxyyyy_1, \
                                         tg_yzz_xxyyyyy_0, tg_yzz_xxyyyyy_1, tg_yzz_xxyyyyz_0, tg_yzz_xxyyyyz_1, \
                                         tg_yzz_xxyyyz_1, tg_yzz_xxyyyzz_0, tg_yzz_xxyyyzz_1, tg_yzz_xxyyzz_1, \
                                         tg_zz_xxxxxxx_0, tg_zz_xxxxxxx_1, tg_zz_xxxxxxy_0, tg_zz_xxxxxxy_1, tg_zz_xxxxxxz_0, \
                                         tg_zz_xxxxxxz_1, tg_zz_xxxxxyy_0, tg_zz_xxxxxyy_1, tg_zz_xxxxxyz_0, tg_zz_xxxxxyz_1, \
                                         tg_zz_xxxxxzz_0, tg_zz_xxxxxzz_1, tg_zz_xxxxyyy_0, tg_zz_xxxxyyy_1, tg_zz_xxxxyyz_0, \
                                         tg_zz_xxxxyyz_1, tg_zz_xxxxyzz_0, tg_zz_xxxxyzz_1, tg_zz_xxxxzzz_0, tg_zz_xxxxzzz_1, \
                                         tg_zz_xxxyyyy_0, tg_zz_xxxyyyy_1, tg_zz_xxxyyyz_0, tg_zz_xxxyyyz_1, tg_zz_xxxyyzz_0, \
                                         tg_zz_xxxyyzz_1, tg_zz_xxxyzzz_0, tg_zz_xxxyzzz_1, tg_zz_xxxzzzz_0, tg_zz_xxxzzzz_1, \
                                         tg_zz_xxyyyyy_0, tg_zz_xxyyyyy_1, tg_zz_xxyyyyz_0, tg_zz_xxyyyyz_1, tg_zz_xxyyyzz_0, \
                                         tg_zz_xxyyyzz_1, wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_y[j]; 

                    tg_yyyy_xxxxxxx_0[j] = pb_y * tg_yyy_xxxxxxx_0[j] + fr * tg_yyy_xxxxxxx_1[j] + 1.5 * fl1_fx * (tg_yy_xxxxxxx_0[j] - tg_yy_xxxxxxx_1[j] * fl1_fza);

                    tg_yyyy_xxxxxxy_0[j] = pb_y * tg_yyy_xxxxxxy_0[j] + fr * tg_yyy_xxxxxxy_1[j] + 1.5 * fl1_fx * (tg_yy_xxxxxxy_0[j] - tg_yy_xxxxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyy_xxxxxx_1[j];

                    tg_yyyy_xxxxxxz_0[j] = pb_y * tg_yyy_xxxxxxz_0[j] + fr * tg_yyy_xxxxxxz_1[j] + 1.5 * fl1_fx * (tg_yy_xxxxxxz_0[j] - tg_yy_xxxxxxz_1[j] * fl1_fza);

                    tg_yyyy_xxxxxyy_0[j] = pb_y * tg_yyy_xxxxxyy_0[j] + fr * tg_yyy_xxxxxyy_1[j] + 1.5 * fl1_fx * (tg_yy_xxxxxyy_0[j] - tg_yy_xxxxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyy_xxxxxy_1[j];

                    tg_yyyy_xxxxxyz_0[j] = pb_y * tg_yyy_xxxxxyz_0[j] + fr * tg_yyy_xxxxxyz_1[j] + 1.5 * fl1_fx * (tg_yy_xxxxxyz_0[j] - tg_yy_xxxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyy_xxxxxz_1[j];

                    tg_yyyy_xxxxxzz_0[j] = pb_y * tg_yyy_xxxxxzz_0[j] + fr * tg_yyy_xxxxxzz_1[j] + 1.5 * fl1_fx * (tg_yy_xxxxxzz_0[j] - tg_yy_xxxxxzz_1[j] * fl1_fza);

                    tg_yyyy_xxxxyyy_0[j] = pb_y * tg_yyy_xxxxyyy_0[j] + fr * tg_yyy_xxxxyyy_1[j] + 1.5 * fl1_fx * (tg_yy_xxxxyyy_0[j] - tg_yy_xxxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyy_xxxxyy_1[j];

                    tg_yyyy_xxxxyyz_0[j] = pb_y * tg_yyy_xxxxyyz_0[j] + fr * tg_yyy_xxxxyyz_1[j] + 1.5 * fl1_fx * (tg_yy_xxxxyyz_0[j] - tg_yy_xxxxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyy_xxxxyz_1[j];

                    tg_yyyy_xxxxyzz_0[j] = pb_y * tg_yyy_xxxxyzz_0[j] + fr * tg_yyy_xxxxyzz_1[j] + 1.5 * fl1_fx * (tg_yy_xxxxyzz_0[j] - tg_yy_xxxxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyy_xxxxzz_1[j];

                    tg_yyyy_xxxxzzz_0[j] = pb_y * tg_yyy_xxxxzzz_0[j] + fr * tg_yyy_xxxxzzz_1[j] + 1.5 * fl1_fx * (tg_yy_xxxxzzz_0[j] - tg_yy_xxxxzzz_1[j] * fl1_fza);

                    tg_yyyy_xxxyyyy_0[j] = pb_y * tg_yyy_xxxyyyy_0[j] + fr * tg_yyy_xxxyyyy_1[j] + 1.5 * fl1_fx * (tg_yy_xxxyyyy_0[j] - tg_yy_xxxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyy_xxxyyy_1[j];

                    tg_yyyy_xxxyyyz_0[j] = pb_y * tg_yyy_xxxyyyz_0[j] + fr * tg_yyy_xxxyyyz_1[j] + 1.5 * fl1_fx * (tg_yy_xxxyyyz_0[j] - tg_yy_xxxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyy_xxxyyz_1[j];

                    tg_yyyy_xxxyyzz_0[j] = pb_y * tg_yyy_xxxyyzz_0[j] + fr * tg_yyy_xxxyyzz_1[j] + 1.5 * fl1_fx * (tg_yy_xxxyyzz_0[j] - tg_yy_xxxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyy_xxxyzz_1[j];

                    tg_yyyy_xxxyzzz_0[j] = pb_y * tg_yyy_xxxyzzz_0[j] + fr * tg_yyy_xxxyzzz_1[j] + 1.5 * fl1_fx * (tg_yy_xxxyzzz_0[j] - tg_yy_xxxyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyy_xxxzzz_1[j];

                    tg_yyyy_xxxzzzz_0[j] = pb_y * tg_yyy_xxxzzzz_0[j] + fr * tg_yyy_xxxzzzz_1[j] + 1.5 * fl1_fx * (tg_yy_xxxzzzz_0[j] - tg_yy_xxxzzzz_1[j] * fl1_fza);

                    tg_yyyy_xxyyyyy_0[j] = pb_y * tg_yyy_xxyyyyy_0[j] + fr * tg_yyy_xxyyyyy_1[j] + 1.5 * fl1_fx * (tg_yy_xxyyyyy_0[j] - tg_yy_xxyyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyy_xxyyyy_1[j];

                    tg_yyyy_xxyyyyz_0[j] = pb_y * tg_yyy_xxyyyyz_0[j] + fr * tg_yyy_xxyyyyz_1[j] + 1.5 * fl1_fx * (tg_yy_xxyyyyz_0[j] - tg_yy_xxyyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyy_xxyyyz_1[j];

                    tg_yyyy_xxyyyzz_0[j] = pb_y * tg_yyy_xxyyyzz_0[j] + fr * tg_yyy_xxyyyzz_1[j] + 1.5 * fl1_fx * (tg_yy_xxyyyzz_0[j] - tg_yy_xxyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyy_xxyyzz_1[j];

                    tg_yyyy_xxyyzzz_0[j] = pb_y * tg_yyy_xxyyzzz_0[j] + fr * tg_yyy_xxyyzzz_1[j] + 1.5 * fl1_fx * (tg_yy_xxyyzzz_0[j] - tg_yy_xxyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyy_xxyzzz_1[j];

                    tg_yyyy_xxyzzzz_0[j] = pb_y * tg_yyy_xxyzzzz_0[j] + fr * tg_yyy_xxyzzzz_1[j] + 1.5 * fl1_fx * (tg_yy_xxyzzzz_0[j] - tg_yy_xxyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyy_xxzzzz_1[j];

                    tg_yyyy_xxzzzzz_0[j] = pb_y * tg_yyy_xxzzzzz_0[j] + fr * tg_yyy_xxzzzzz_1[j] + 1.5 * fl1_fx * (tg_yy_xxzzzzz_0[j] - tg_yy_xxzzzzz_1[j] * fl1_fza);

                    tg_yyyy_xyyyyyy_0[j] = pb_y * tg_yyy_xyyyyyy_0[j] + fr * tg_yyy_xyyyyyy_1[j] + 1.5 * fl1_fx * (tg_yy_xyyyyyy_0[j] - tg_yy_xyyyyyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yyy_xyyyyy_1[j];

                    tg_yyyy_xyyyyyz_0[j] = pb_y * tg_yyy_xyyyyyz_0[j] + fr * tg_yyy_xyyyyyz_1[j] + 1.5 * fl1_fx * (tg_yy_xyyyyyz_0[j] - tg_yy_xyyyyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyy_xyyyyz_1[j];

                    tg_yyyy_xyyyyzz_0[j] = pb_y * tg_yyy_xyyyyzz_0[j] + fr * tg_yyy_xyyyyzz_1[j] + 1.5 * fl1_fx * (tg_yy_xyyyyzz_0[j] - tg_yy_xyyyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyy_xyyyzz_1[j];

                    tg_yyyy_xyyyzzz_0[j] = pb_y * tg_yyy_xyyyzzz_0[j] + fr * tg_yyy_xyyyzzz_1[j] + 1.5 * fl1_fx * (tg_yy_xyyyzzz_0[j] - tg_yy_xyyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyy_xyyzzz_1[j];

                    tg_yyyy_xyyzzzz_0[j] = pb_y * tg_yyy_xyyzzzz_0[j] + fr * tg_yyy_xyyzzzz_1[j] + 1.5 * fl1_fx * (tg_yy_xyyzzzz_0[j] - tg_yy_xyyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyy_xyzzzz_1[j];

                    tg_yyyy_xyzzzzz_0[j] = pb_y * tg_yyy_xyzzzzz_0[j] + fr * tg_yyy_xyzzzzz_1[j] + 1.5 * fl1_fx * (tg_yy_xyzzzzz_0[j] - tg_yy_xyzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyy_xzzzzz_1[j];

                    tg_yyyy_xzzzzzz_0[j] = pb_y * tg_yyy_xzzzzzz_0[j] + fr * tg_yyy_xzzzzzz_1[j] + 1.5 * fl1_fx * (tg_yy_xzzzzzz_0[j] - tg_yy_xzzzzzz_1[j] * fl1_fza);

                    tg_yyyy_yyyyyyy_0[j] = pb_y * tg_yyy_yyyyyyy_0[j] + fr * tg_yyy_yyyyyyy_1[j] + 1.5 * fl1_fx * (tg_yy_yyyyyyy_0[j] - tg_yy_yyyyyyy_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_yyy_yyyyyy_1[j];

                    tg_yyyy_yyyyyyz_0[j] = pb_y * tg_yyy_yyyyyyz_0[j] + fr * tg_yyy_yyyyyyz_1[j] + 1.5 * fl1_fx * (tg_yy_yyyyyyz_0[j] - tg_yy_yyyyyyz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yyy_yyyyyz_1[j];

                    tg_yyyy_yyyyyzz_0[j] = pb_y * tg_yyy_yyyyyzz_0[j] + fr * tg_yyy_yyyyyzz_1[j] + 1.5 * fl1_fx * (tg_yy_yyyyyzz_0[j] - tg_yy_yyyyyzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyy_yyyyzz_1[j];

                    tg_yyyy_yyyyzzz_0[j] = pb_y * tg_yyy_yyyyzzz_0[j] + fr * tg_yyy_yyyyzzz_1[j] + 1.5 * fl1_fx * (tg_yy_yyyyzzz_0[j] - tg_yy_yyyyzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyy_yyyzzz_1[j];

                    tg_yyyy_yyyzzzz_0[j] = pb_y * tg_yyy_yyyzzzz_0[j] + fr * tg_yyy_yyyzzzz_1[j] + 1.5 * fl1_fx * (tg_yy_yyyzzzz_0[j] - tg_yy_yyyzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyy_yyzzzz_1[j];

                    tg_yyyy_yyzzzzz_0[j] = pb_y * tg_yyy_yyzzzzz_0[j] + fr * tg_yyy_yyzzzzz_1[j] + 1.5 * fl1_fx * (tg_yy_yyzzzzz_0[j] - tg_yy_yyzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyy_yzzzzz_1[j];

                    tg_yyyy_yzzzzzz_0[j] = pb_y * tg_yyy_yzzzzzz_0[j] + fr * tg_yyy_yzzzzzz_1[j] + 1.5 * fl1_fx * (tg_yy_yzzzzzz_0[j] - tg_yy_yzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyy_zzzzzz_1[j];

                    tg_yyyy_zzzzzzz_0[j] = pb_y * tg_yyy_zzzzzzz_0[j] + fr * tg_yyy_zzzzzzz_1[j] + 1.5 * fl1_fx * (tg_yy_zzzzzzz_0[j] - tg_yy_zzzzzzz_1[j] * fl1_fza);

                    tg_yyyz_xxxxxxx_0[j] = pb_y * tg_yyz_xxxxxxx_0[j] + fr * tg_yyz_xxxxxxx_1[j] + fl1_fx * (tg_yz_xxxxxxx_0[j] - tg_yz_xxxxxxx_1[j] * fl1_fza);

                    tg_yyyz_xxxxxxy_0[j] = pb_y * tg_yyz_xxxxxxy_0[j] + fr * tg_yyz_xxxxxxy_1[j] + fl1_fx * (tg_yz_xxxxxxy_0[j] - tg_yz_xxxxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyz_xxxxxx_1[j];

                    tg_yyyz_xxxxxxz_0[j] = pb_y * tg_yyz_xxxxxxz_0[j] + fr * tg_yyz_xxxxxxz_1[j] + fl1_fx * (tg_yz_xxxxxxz_0[j] - tg_yz_xxxxxxz_1[j] * fl1_fza);

                    tg_yyyz_xxxxxyy_0[j] = pb_y * tg_yyz_xxxxxyy_0[j] + fr * tg_yyz_xxxxxyy_1[j] + fl1_fx * (tg_yz_xxxxxyy_0[j] - tg_yz_xxxxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyz_xxxxxy_1[j];

                    tg_yyyz_xxxxxyz_0[j] = pb_y * tg_yyz_xxxxxyz_0[j] + fr * tg_yyz_xxxxxyz_1[j] + fl1_fx * (tg_yz_xxxxxyz_0[j] - tg_yz_xxxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyz_xxxxxz_1[j];

                    tg_yyyz_xxxxxzz_0[j] = pb_y * tg_yyz_xxxxxzz_0[j] + fr * tg_yyz_xxxxxzz_1[j] + fl1_fx * (tg_yz_xxxxxzz_0[j] - tg_yz_xxxxxzz_1[j] * fl1_fza);

                    tg_yyyz_xxxxyyy_0[j] = pb_y * tg_yyz_xxxxyyy_0[j] + fr * tg_yyz_xxxxyyy_1[j] + fl1_fx * (tg_yz_xxxxyyy_0[j] - tg_yz_xxxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyz_xxxxyy_1[j];

                    tg_yyyz_xxxxyyz_0[j] = pb_y * tg_yyz_xxxxyyz_0[j] + fr * tg_yyz_xxxxyyz_1[j] + fl1_fx * (tg_yz_xxxxyyz_0[j] - tg_yz_xxxxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyz_xxxxyz_1[j];

                    tg_yyyz_xxxxyzz_0[j] = pb_y * tg_yyz_xxxxyzz_0[j] + fr * tg_yyz_xxxxyzz_1[j] + fl1_fx * (tg_yz_xxxxyzz_0[j] - tg_yz_xxxxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyz_xxxxzz_1[j];

                    tg_yyyz_xxxxzzz_0[j] = pb_y * tg_yyz_xxxxzzz_0[j] + fr * tg_yyz_xxxxzzz_1[j] + fl1_fx * (tg_yz_xxxxzzz_0[j] - tg_yz_xxxxzzz_1[j] * fl1_fza);

                    tg_yyyz_xxxyyyy_0[j] = pb_y * tg_yyz_xxxyyyy_0[j] + fr * tg_yyz_xxxyyyy_1[j] + fl1_fx * (tg_yz_xxxyyyy_0[j] - tg_yz_xxxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyz_xxxyyy_1[j];

                    tg_yyyz_xxxyyyz_0[j] = pb_y * tg_yyz_xxxyyyz_0[j] + fr * tg_yyz_xxxyyyz_1[j] + fl1_fx * (tg_yz_xxxyyyz_0[j] - tg_yz_xxxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyz_xxxyyz_1[j];

                    tg_yyyz_xxxyyzz_0[j] = pb_y * tg_yyz_xxxyyzz_0[j] + fr * tg_yyz_xxxyyzz_1[j] + fl1_fx * (tg_yz_xxxyyzz_0[j] - tg_yz_xxxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyz_xxxyzz_1[j];

                    tg_yyyz_xxxyzzz_0[j] = pb_y * tg_yyz_xxxyzzz_0[j] + fr * tg_yyz_xxxyzzz_1[j] + fl1_fx * (tg_yz_xxxyzzz_0[j] - tg_yz_xxxyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyz_xxxzzz_1[j];

                    tg_yyyz_xxxzzzz_0[j] = pb_y * tg_yyz_xxxzzzz_0[j] + fr * tg_yyz_xxxzzzz_1[j] + fl1_fx * (tg_yz_xxxzzzz_0[j] - tg_yz_xxxzzzz_1[j] * fl1_fza);

                    tg_yyyz_xxyyyyy_0[j] = pb_y * tg_yyz_xxyyyyy_0[j] + fr * tg_yyz_xxyyyyy_1[j] + fl1_fx * (tg_yz_xxyyyyy_0[j] - tg_yz_xxyyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyz_xxyyyy_1[j];

                    tg_yyyz_xxyyyyz_0[j] = pb_y * tg_yyz_xxyyyyz_0[j] + fr * tg_yyz_xxyyyyz_1[j] + fl1_fx * (tg_yz_xxyyyyz_0[j] - tg_yz_xxyyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyz_xxyyyz_1[j];

                    tg_yyyz_xxyyyzz_0[j] = pb_y * tg_yyz_xxyyyzz_0[j] + fr * tg_yyz_xxyyyzz_1[j] + fl1_fx * (tg_yz_xxyyyzz_0[j] - tg_yz_xxyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyz_xxyyzz_1[j];

                    tg_yyyz_xxyyzzz_0[j] = pb_y * tg_yyz_xxyyzzz_0[j] + fr * tg_yyz_xxyyzzz_1[j] + fl1_fx * (tg_yz_xxyyzzz_0[j] - tg_yz_xxyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyz_xxyzzz_1[j];

                    tg_yyyz_xxyzzzz_0[j] = pb_y * tg_yyz_xxyzzzz_0[j] + fr * tg_yyz_xxyzzzz_1[j] + fl1_fx * (tg_yz_xxyzzzz_0[j] - tg_yz_xxyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyz_xxzzzz_1[j];

                    tg_yyyz_xxzzzzz_0[j] = pb_y * tg_yyz_xxzzzzz_0[j] + fr * tg_yyz_xxzzzzz_1[j] + fl1_fx * (tg_yz_xxzzzzz_0[j] - tg_yz_xxzzzzz_1[j] * fl1_fza);

                    tg_yyyz_xyyyyyy_0[j] = pb_y * tg_yyz_xyyyyyy_0[j] + fr * tg_yyz_xyyyyyy_1[j] + fl1_fx * (tg_yz_xyyyyyy_0[j] - tg_yz_xyyyyyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yyz_xyyyyy_1[j];

                    tg_yyyz_xyyyyyz_0[j] = pb_y * tg_yyz_xyyyyyz_0[j] + fr * tg_yyz_xyyyyyz_1[j] + fl1_fx * (tg_yz_xyyyyyz_0[j] - tg_yz_xyyyyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyz_xyyyyz_1[j];

                    tg_yyyz_xyyyyzz_0[j] = pb_y * tg_yyz_xyyyyzz_0[j] + fr * tg_yyz_xyyyyzz_1[j] + fl1_fx * (tg_yz_xyyyyzz_0[j] - tg_yz_xyyyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyz_xyyyzz_1[j];

                    tg_yyyz_xyyyzzz_0[j] = pb_y * tg_yyz_xyyyzzz_0[j] + fr * tg_yyz_xyyyzzz_1[j] + fl1_fx * (tg_yz_xyyyzzz_0[j] - tg_yz_xyyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyz_xyyzzz_1[j];

                    tg_yyyz_xyyzzzz_0[j] = pb_y * tg_yyz_xyyzzzz_0[j] + fr * tg_yyz_xyyzzzz_1[j] + fl1_fx * (tg_yz_xyyzzzz_0[j] - tg_yz_xyyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyz_xyzzzz_1[j];

                    tg_yyyz_xyzzzzz_0[j] = pb_y * tg_yyz_xyzzzzz_0[j] + fr * tg_yyz_xyzzzzz_1[j] + fl1_fx * (tg_yz_xyzzzzz_0[j] - tg_yz_xyzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyz_xzzzzz_1[j];

                    tg_yyyz_xzzzzzz_0[j] = pb_y * tg_yyz_xzzzzzz_0[j] + fr * tg_yyz_xzzzzzz_1[j] + fl1_fx * (tg_yz_xzzzzzz_0[j] - tg_yz_xzzzzzz_1[j] * fl1_fza);

                    tg_yyyz_yyyyyyy_0[j] = pb_y * tg_yyz_yyyyyyy_0[j] + fr * tg_yyz_yyyyyyy_1[j] + fl1_fx * (tg_yz_yyyyyyy_0[j] - tg_yz_yyyyyyy_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_yyz_yyyyyy_1[j];

                    tg_yyyz_yyyyyyz_0[j] = pb_y * tg_yyz_yyyyyyz_0[j] + fr * tg_yyz_yyyyyyz_1[j] + fl1_fx * (tg_yz_yyyyyyz_0[j] - tg_yz_yyyyyyz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yyz_yyyyyz_1[j];

                    tg_yyyz_yyyyyzz_0[j] = pb_y * tg_yyz_yyyyyzz_0[j] + fr * tg_yyz_yyyyyzz_1[j] + fl1_fx * (tg_yz_yyyyyzz_0[j] - tg_yz_yyyyyzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyz_yyyyzz_1[j];

                    tg_yyyz_yyyyzzz_0[j] = pb_y * tg_yyz_yyyyzzz_0[j] + fr * tg_yyz_yyyyzzz_1[j] + fl1_fx * (tg_yz_yyyyzzz_0[j] - tg_yz_yyyyzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyz_yyyzzz_1[j];

                    tg_yyyz_yyyzzzz_0[j] = pb_y * tg_yyz_yyyzzzz_0[j] + fr * tg_yyz_yyyzzzz_1[j] + fl1_fx * (tg_yz_yyyzzzz_0[j] - tg_yz_yyyzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyz_yyzzzz_1[j];

                    tg_yyyz_yyzzzzz_0[j] = pb_y * tg_yyz_yyzzzzz_0[j] + fr * tg_yyz_yyzzzzz_1[j] + fl1_fx * (tg_yz_yyzzzzz_0[j] - tg_yz_yyzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyz_yzzzzz_1[j];

                    tg_yyyz_yzzzzzz_0[j] = pb_y * tg_yyz_yzzzzzz_0[j] + fr * tg_yyz_yzzzzzz_1[j] + fl1_fx * (tg_yz_yzzzzzz_0[j] - tg_yz_yzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyz_zzzzzz_1[j];

                    tg_yyyz_zzzzzzz_0[j] = pb_y * tg_yyz_zzzzzzz_0[j] + fr * tg_yyz_zzzzzzz_1[j] + fl1_fx * (tg_yz_zzzzzzz_0[j] - tg_yz_zzzzzzz_1[j] * fl1_fza);

                    tg_yyzz_xxxxxxx_0[j] = pb_y * tg_yzz_xxxxxxx_0[j] + fr * tg_yzz_xxxxxxx_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxxxx_0[j] - tg_zz_xxxxxxx_1[j] * fl1_fza);

                    tg_yyzz_xxxxxxy_0[j] = pb_y * tg_yzz_xxxxxxy_0[j] + fr * tg_yzz_xxxxxxy_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxxxy_0[j] - tg_zz_xxxxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzz_xxxxxx_1[j];

                    tg_yyzz_xxxxxxz_0[j] = pb_y * tg_yzz_xxxxxxz_0[j] + fr * tg_yzz_xxxxxxz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxxxz_0[j] - tg_zz_xxxxxxz_1[j] * fl1_fza);

                    tg_yyzz_xxxxxyy_0[j] = pb_y * tg_yzz_xxxxxyy_0[j] + fr * tg_yzz_xxxxxyy_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxxyy_0[j] - tg_zz_xxxxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yzz_xxxxxy_1[j];

                    tg_yyzz_xxxxxyz_0[j] = pb_y * tg_yzz_xxxxxyz_0[j] + fr * tg_yzz_xxxxxyz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxxyz_0[j] - tg_zz_xxxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzz_xxxxxz_1[j];

                    tg_yyzz_xxxxxzz_0[j] = pb_y * tg_yzz_xxxxxzz_0[j] + fr * tg_yzz_xxxxxzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxxzz_0[j] - tg_zz_xxxxxzz_1[j] * fl1_fza);

                    tg_yyzz_xxxxyyy_0[j] = pb_y * tg_yzz_xxxxyyy_0[j] + fr * tg_yzz_xxxxyyy_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxyyy_0[j] - tg_zz_xxxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzz_xxxxyy_1[j];

                    tg_yyzz_xxxxyyz_0[j] = pb_y * tg_yzz_xxxxyyz_0[j] + fr * tg_yzz_xxxxyyz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxyyz_0[j] - tg_zz_xxxxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yzz_xxxxyz_1[j];

                    tg_yyzz_xxxxyzz_0[j] = pb_y * tg_yzz_xxxxyzz_0[j] + fr * tg_yzz_xxxxyzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxyzz_0[j] - tg_zz_xxxxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzz_xxxxzz_1[j];

                    tg_yyzz_xxxxzzz_0[j] = pb_y * tg_yzz_xxxxzzz_0[j] + fr * tg_yzz_xxxxzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxzzz_0[j] - tg_zz_xxxxzzz_1[j] * fl1_fza);

                    tg_yyzz_xxxyyyy_0[j] = pb_y * tg_yzz_xxxyyyy_0[j] + fr * tg_yzz_xxxyyyy_1[j] + 0.5 * fl1_fx * (tg_zz_xxxyyyy_0[j] - tg_zz_xxxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzz_xxxyyy_1[j];

                    tg_yyzz_xxxyyyz_0[j] = pb_y * tg_yzz_xxxyyyz_0[j] + fr * tg_yzz_xxxyyyz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxyyyz_0[j] - tg_zz_xxxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzz_xxxyyz_1[j];

                    tg_yyzz_xxxyyzz_0[j] = pb_y * tg_yzz_xxxyyzz_0[j] + fr * tg_yzz_xxxyyzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxyyzz_0[j] - tg_zz_xxxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yzz_xxxyzz_1[j];

                    tg_yyzz_xxxyzzz_0[j] = pb_y * tg_yzz_xxxyzzz_0[j] + fr * tg_yzz_xxxyzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxyzzz_0[j] - tg_zz_xxxyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzz_xxxzzz_1[j];

                    tg_yyzz_xxxzzzz_0[j] = pb_y * tg_yzz_xxxzzzz_0[j] + fr * tg_yzz_xxxzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxzzzz_0[j] - tg_zz_xxxzzzz_1[j] * fl1_fza);

                    tg_yyzz_xxyyyyy_0[j] = pb_y * tg_yzz_xxyyyyy_0[j] + fr * tg_yzz_xxyyyyy_1[j] + 0.5 * fl1_fx * (tg_zz_xxyyyyy_0[j] - tg_zz_xxyyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yzz_xxyyyy_1[j];

                    tg_yyzz_xxyyyyz_0[j] = pb_y * tg_yzz_xxyyyyz_0[j] + fr * tg_yzz_xxyyyyz_1[j] + 0.5 * fl1_fx * (tg_zz_xxyyyyz_0[j] - tg_zz_xxyyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzz_xxyyyz_1[j];

                    tg_yyzz_xxyyyzz_0[j] = pb_y * tg_yzz_xxyyyzz_0[j] + fr * tg_yzz_xxyyyzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxyyyzz_0[j] - tg_zz_xxyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzz_xxyyzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSK_450_540(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (450,540)

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
                                             {4, -1, -1, -1},
                                             {7, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_7_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {7, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_7_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_7_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_7_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_2_7_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {7, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_2_7_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {7, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {6, -1, -1, -1}, 
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

                auto tg_yzz_xxyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 306); 

                auto tg_yzz_xxyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 307); 

                auto tg_yzz_xxzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 308); 

                auto tg_yzz_xyyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 309); 

                auto tg_yzz_xyyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 310); 

                auto tg_yzz_xyyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 311); 

                auto tg_yzz_xyyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 312); 

                auto tg_yzz_xyyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 313); 

                auto tg_yzz_xyzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 314); 

                auto tg_yzz_xzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 315); 

                auto tg_yzz_yyyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 316); 

                auto tg_yzz_yyyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 317); 

                auto tg_yzz_yyyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 318); 

                auto tg_yzz_yyyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 319); 

                auto tg_yzz_yyyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 320); 

                auto tg_yzz_yyzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 321); 

                auto tg_yzz_yzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 322); 

                auto tg_yzz_zzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 323); 

                auto tg_zzz_xxxxxxx_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 324); 

                auto tg_zzz_xxxxxxy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 325); 

                auto tg_zzz_xxxxxxz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 326); 

                auto tg_zzz_xxxxxyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 327); 

                auto tg_zzz_xxxxxyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 328); 

                auto tg_zzz_xxxxxzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 329); 

                auto tg_zzz_xxxxyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 330); 

                auto tg_zzz_xxxxyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 331); 

                auto tg_zzz_xxxxyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 332); 

                auto tg_zzz_xxxxzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 333); 

                auto tg_zzz_xxxyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 334); 

                auto tg_zzz_xxxyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 335); 

                auto tg_zzz_xxxyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 336); 

                auto tg_zzz_xxxyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 337); 

                auto tg_zzz_xxxzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 338); 

                auto tg_zzz_xxyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 339); 

                auto tg_zzz_xxyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 340); 

                auto tg_zzz_xxyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 341); 

                auto tg_zzz_xxyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 342); 

                auto tg_zzz_xxyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 343); 

                auto tg_zzz_xxzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 344); 

                auto tg_zzz_xyyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 345); 

                auto tg_zzz_xyyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 346); 

                auto tg_zzz_xyyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 347); 

                auto tg_zzz_xyyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 348); 

                auto tg_zzz_xyyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 349); 

                auto tg_zzz_xyzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 350); 

                auto tg_zzz_xzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 351); 

                auto tg_zzz_yyyyyyy_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 352); 

                auto tg_zzz_yyyyyyz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 353); 

                auto tg_zzz_yyyyyzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 354); 

                auto tg_zzz_yyyyzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 355); 

                auto tg_zzz_yyyzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 356); 

                auto tg_zzz_yyzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 357); 

                auto tg_zzz_yzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 358); 

                auto tg_zzz_zzzzzzz_0 = primBuffer[pidx_g_3_7_m0].data(360 * idx + 359); 

                auto tg_yzz_xxyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 306); 

                auto tg_yzz_xxyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 307); 

                auto tg_yzz_xxzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 308); 

                auto tg_yzz_xyyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 309); 

                auto tg_yzz_xyyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 310); 

                auto tg_yzz_xyyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 311); 

                auto tg_yzz_xyyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 312); 

                auto tg_yzz_xyyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 313); 

                auto tg_yzz_xyzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 314); 

                auto tg_yzz_xzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 315); 

                auto tg_yzz_yyyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 316); 

                auto tg_yzz_yyyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 317); 

                auto tg_yzz_yyyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 318); 

                auto tg_yzz_yyyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 319); 

                auto tg_yzz_yyyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 320); 

                auto tg_yzz_yyzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 321); 

                auto tg_yzz_yzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 322); 

                auto tg_yzz_zzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 323); 

                auto tg_zzz_xxxxxxx_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 324); 

                auto tg_zzz_xxxxxxy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 325); 

                auto tg_zzz_xxxxxxz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 326); 

                auto tg_zzz_xxxxxyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 327); 

                auto tg_zzz_xxxxxyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 328); 

                auto tg_zzz_xxxxxzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 329); 

                auto tg_zzz_xxxxyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 330); 

                auto tg_zzz_xxxxyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 331); 

                auto tg_zzz_xxxxyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 332); 

                auto tg_zzz_xxxxzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 333); 

                auto tg_zzz_xxxyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 334); 

                auto tg_zzz_xxxyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 335); 

                auto tg_zzz_xxxyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 336); 

                auto tg_zzz_xxxyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 337); 

                auto tg_zzz_xxxzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 338); 

                auto tg_zzz_xxyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 339); 

                auto tg_zzz_xxyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 340); 

                auto tg_zzz_xxyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 341); 

                auto tg_zzz_xxyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 342); 

                auto tg_zzz_xxyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 343); 

                auto tg_zzz_xxzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 344); 

                auto tg_zzz_xyyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 345); 

                auto tg_zzz_xyyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 346); 

                auto tg_zzz_xyyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 347); 

                auto tg_zzz_xyyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 348); 

                auto tg_zzz_xyyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 349); 

                auto tg_zzz_xyzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 350); 

                auto tg_zzz_xzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 351); 

                auto tg_zzz_yyyyyyy_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 352); 

                auto tg_zzz_yyyyyyz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 353); 

                auto tg_zzz_yyyyyzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 354); 

                auto tg_zzz_yyyyzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 355); 

                auto tg_zzz_yyyzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 356); 

                auto tg_zzz_yyzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 357); 

                auto tg_zzz_yzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 358); 

                auto tg_zzz_zzzzzzz_1 = primBuffer[pidx_g_3_7_m1].data(360 * idx + 359); 

                auto tg_zz_xxxxxxx_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 180); 

                auto tg_zz_xxxxxxy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 181); 

                auto tg_zz_xxxxxxz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 182); 

                auto tg_zz_xxxxxyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 183); 

                auto tg_zz_xxxxxyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 184); 

                auto tg_zz_xxxxxzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 185); 

                auto tg_zz_xxxxyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 186); 

                auto tg_zz_xxxxyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 187); 

                auto tg_zz_xxxxyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 188); 

                auto tg_zz_xxxxzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 189); 

                auto tg_zz_xxxyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 190); 

                auto tg_zz_xxxyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 191); 

                auto tg_zz_xxxyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 192); 

                auto tg_zz_xxxyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 193); 

                auto tg_zz_xxxzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 194); 

                auto tg_zz_xxyyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 195); 

                auto tg_zz_xxyyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 196); 

                auto tg_zz_xxyyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 197); 

                auto tg_zz_xxyyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 198); 

                auto tg_zz_xxyzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 199); 

                auto tg_zz_xxzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 200); 

                auto tg_zz_xyyyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 201); 

                auto tg_zz_xyyyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 202); 

                auto tg_zz_xyyyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 203); 

                auto tg_zz_xyyyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 204); 

                auto tg_zz_xyyzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 205); 

                auto tg_zz_xyzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 206); 

                auto tg_zz_xzzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 207); 

                auto tg_zz_yyyyyyy_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 208); 

                auto tg_zz_yyyyyyz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 209); 

                auto tg_zz_yyyyyzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 210); 

                auto tg_zz_yyyyzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 211); 

                auto tg_zz_yyyzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 212); 

                auto tg_zz_yyzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 213); 

                auto tg_zz_yzzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 214); 

                auto tg_zz_zzzzzzz_0 = primBuffer[pidx_g_2_7_m0].data(216 * idx + 215); 

                auto tg_zz_xxxxxxx_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 180); 

                auto tg_zz_xxxxxxy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 181); 

                auto tg_zz_xxxxxxz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 182); 

                auto tg_zz_xxxxxyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 183); 

                auto tg_zz_xxxxxyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 184); 

                auto tg_zz_xxxxxzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 185); 

                auto tg_zz_xxxxyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 186); 

                auto tg_zz_xxxxyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 187); 

                auto tg_zz_xxxxyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 188); 

                auto tg_zz_xxxxzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 189); 

                auto tg_zz_xxxyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 190); 

                auto tg_zz_xxxyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 191); 

                auto tg_zz_xxxyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 192); 

                auto tg_zz_xxxyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 193); 

                auto tg_zz_xxxzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 194); 

                auto tg_zz_xxyyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 195); 

                auto tg_zz_xxyyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 196); 

                auto tg_zz_xxyyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 197); 

                auto tg_zz_xxyyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 198); 

                auto tg_zz_xxyzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 199); 

                auto tg_zz_xxzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 200); 

                auto tg_zz_xyyyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 201); 

                auto tg_zz_xyyyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 202); 

                auto tg_zz_xyyyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 203); 

                auto tg_zz_xyyyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 204); 

                auto tg_zz_xyyzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 205); 

                auto tg_zz_xyzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 206); 

                auto tg_zz_xzzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 207); 

                auto tg_zz_yyyyyyy_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 208); 

                auto tg_zz_yyyyyyz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 209); 

                auto tg_zz_yyyyyzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 210); 

                auto tg_zz_yyyyzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 211); 

                auto tg_zz_yyyzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 212); 

                auto tg_zz_yyzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 213); 

                auto tg_zz_yzzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 214); 

                auto tg_zz_zzzzzzz_1 = primBuffer[pidx_g_2_7_m1].data(216 * idx + 215); 

                auto tg_yzz_xxyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 237); 

                auto tg_yzz_xxzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 238); 

                auto tg_yzz_xyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 239); 

                auto tg_yzz_xyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 240); 

                auto tg_yzz_xyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 241); 

                auto tg_yzz_xyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 242); 

                auto tg_yzz_xyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 243); 

                auto tg_yzz_xzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 244); 

                auto tg_yzz_yyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 245); 

                auto tg_yzz_yyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 246); 

                auto tg_yzz_yyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 247); 

                auto tg_yzz_yyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 248); 

                auto tg_yzz_yyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 249); 

                auto tg_yzz_yzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 250); 

                auto tg_yzz_zzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 251); 

                auto tg_zzz_xxxxxx_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 252); 

                auto tg_zzz_xxxxxy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 253); 

                auto tg_zzz_xxxxxz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 254); 

                auto tg_zzz_xxxxyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 255); 

                auto tg_zzz_xxxxyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 256); 

                auto tg_zzz_xxxxzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 257); 

                auto tg_zzz_xxxyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 258); 

                auto tg_zzz_xxxyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 259); 

                auto tg_zzz_xxxyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 260); 

                auto tg_zzz_xxxzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 261); 

                auto tg_zzz_xxyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 262); 

                auto tg_zzz_xxyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 263); 

                auto tg_zzz_xxyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 264); 

                auto tg_zzz_xxyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 265); 

                auto tg_zzz_xxzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 266); 

                auto tg_zzz_xyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 267); 

                auto tg_zzz_xyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 268); 

                auto tg_zzz_xyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 269); 

                auto tg_zzz_xyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 270); 

                auto tg_zzz_xyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 271); 

                auto tg_zzz_xzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 272); 

                auto tg_zzz_yyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 273); 

                auto tg_zzz_yyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 274); 

                auto tg_zzz_yyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 275); 

                auto tg_zzz_yyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 276); 

                auto tg_zzz_yyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 277); 

                auto tg_zzz_yzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 278); 

                auto tg_zzz_zzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 279); 

                // set up pointers to integrals

                auto tg_yyzz_xxyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 450); 

                auto tg_yyzz_xxyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 451); 

                auto tg_yyzz_xxzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 452); 

                auto tg_yyzz_xyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 453); 

                auto tg_yyzz_xyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 454); 

                auto tg_yyzz_xyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 455); 

                auto tg_yyzz_xyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 456); 

                auto tg_yyzz_xyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 457); 

                auto tg_yyzz_xyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 458); 

                auto tg_yyzz_xzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 459); 

                auto tg_yyzz_yyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 460); 

                auto tg_yyzz_yyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 461); 

                auto tg_yyzz_yyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 462); 

                auto tg_yyzz_yyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 463); 

                auto tg_yyzz_yyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 464); 

                auto tg_yyzz_yyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 465); 

                auto tg_yyzz_yzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 466); 

                auto tg_yyzz_zzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 467); 

                auto tg_yzzz_xxxxxxx_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 468); 

                auto tg_yzzz_xxxxxxy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 469); 

                auto tg_yzzz_xxxxxxz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 470); 

                auto tg_yzzz_xxxxxyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 471); 

                auto tg_yzzz_xxxxxyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 472); 

                auto tg_yzzz_xxxxxzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 473); 

                auto tg_yzzz_xxxxyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 474); 

                auto tg_yzzz_xxxxyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 475); 

                auto tg_yzzz_xxxxyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 476); 

                auto tg_yzzz_xxxxzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 477); 

                auto tg_yzzz_xxxyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 478); 

                auto tg_yzzz_xxxyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 479); 

                auto tg_yzzz_xxxyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 480); 

                auto tg_yzzz_xxxyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 481); 

                auto tg_yzzz_xxxzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 482); 

                auto tg_yzzz_xxyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 483); 

                auto tg_yzzz_xxyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 484); 

                auto tg_yzzz_xxyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 485); 

                auto tg_yzzz_xxyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 486); 

                auto tg_yzzz_xxyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 487); 

                auto tg_yzzz_xxzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 488); 

                auto tg_yzzz_xyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 489); 

                auto tg_yzzz_xyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 490); 

                auto tg_yzzz_xyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 491); 

                auto tg_yzzz_xyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 492); 

                auto tg_yzzz_xyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 493); 

                auto tg_yzzz_xyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 494); 

                auto tg_yzzz_xzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 495); 

                auto tg_yzzz_yyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 496); 

                auto tg_yzzz_yyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 497); 

                auto tg_yzzz_yyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 498); 

                auto tg_yzzz_yyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 499); 

                auto tg_yzzz_yyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 500); 

                auto tg_yzzz_yyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 501); 

                auto tg_yzzz_yzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 502); 

                auto tg_yzzz_zzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 503); 

                auto tg_zzzz_xxxxxxx_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 504); 

                auto tg_zzzz_xxxxxxy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 505); 

                auto tg_zzzz_xxxxxxz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 506); 

                auto tg_zzzz_xxxxxyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 507); 

                auto tg_zzzz_xxxxxyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 508); 

                auto tg_zzzz_xxxxxzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 509); 

                auto tg_zzzz_xxxxyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 510); 

                auto tg_zzzz_xxxxyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 511); 

                auto tg_zzzz_xxxxyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 512); 

                auto tg_zzzz_xxxxzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 513); 

                auto tg_zzzz_xxxyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 514); 

                auto tg_zzzz_xxxyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 515); 

                auto tg_zzzz_xxxyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 516); 

                auto tg_zzzz_xxxyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 517); 

                auto tg_zzzz_xxxzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 518); 

                auto tg_zzzz_xxyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 519); 

                auto tg_zzzz_xxyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 520); 

                auto tg_zzzz_xxyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 521); 

                auto tg_zzzz_xxyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 522); 

                auto tg_zzzz_xxyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 523); 

                auto tg_zzzz_xxzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 524); 

                auto tg_zzzz_xyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 525); 

                auto tg_zzzz_xyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 526); 

                auto tg_zzzz_xyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 527); 

                auto tg_zzzz_xyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 528); 

                auto tg_zzzz_xyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 529); 

                auto tg_zzzz_xyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 530); 

                auto tg_zzzz_xzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 531); 

                auto tg_zzzz_yyyyyyy_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 532); 

                auto tg_zzzz_yyyyyyz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 533); 

                auto tg_zzzz_yyyyyzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 534); 

                auto tg_zzzz_yyyyzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 535); 

                auto tg_zzzz_yyyzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 536); 

                auto tg_zzzz_yyzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 537); 

                auto tg_zzzz_yzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 538); 

                auto tg_zzzz_zzzzzzz_0 = primBuffer[pidx_g_4_7_m0].data(540 * idx + 539); 

                // Batch of Integrals (450,540)

                #pragma omp simd aligned(fxn, fza, tg_yyzz_xxyyzzz_0, tg_yyzz_xxyzzzz_0, tg_yyzz_xxzzzzz_0, \
                                         tg_yyzz_xyyyyyy_0, tg_yyzz_xyyyyyz_0, tg_yyzz_xyyyyzz_0, tg_yyzz_xyyyzzz_0, \
                                         tg_yyzz_xyyzzzz_0, tg_yyzz_xyzzzzz_0, tg_yyzz_xzzzzzz_0, tg_yyzz_yyyyyyy_0, \
                                         tg_yyzz_yyyyyyz_0, tg_yyzz_yyyyyzz_0, tg_yyzz_yyyyzzz_0, tg_yyzz_yyyzzzz_0, \
                                         tg_yyzz_yyzzzzz_0, tg_yyzz_yzzzzzz_0, tg_yyzz_zzzzzzz_0, tg_yzz_xxyyzzz_0, \
                                         tg_yzz_xxyyzzz_1, tg_yzz_xxyzzz_1, tg_yzz_xxyzzzz_0, tg_yzz_xxyzzzz_1, \
                                         tg_yzz_xxzzzz_1, tg_yzz_xxzzzzz_0, tg_yzz_xxzzzzz_1, tg_yzz_xyyyyy_1, \
                                         tg_yzz_xyyyyyy_0, tg_yzz_xyyyyyy_1, tg_yzz_xyyyyyz_0, tg_yzz_xyyyyyz_1, \
                                         tg_yzz_xyyyyz_1, tg_yzz_xyyyyzz_0, tg_yzz_xyyyyzz_1, tg_yzz_xyyyzz_1, \
                                         tg_yzz_xyyyzzz_0, tg_yzz_xyyyzzz_1, tg_yzz_xyyzzz_1, tg_yzz_xyyzzzz_0, \
                                         tg_yzz_xyyzzzz_1, tg_yzz_xyzzzz_1, tg_yzz_xyzzzzz_0, tg_yzz_xyzzzzz_1, \
                                         tg_yzz_xzzzzz_1, tg_yzz_xzzzzzz_0, tg_yzz_xzzzzzz_1, tg_yzz_yyyyyy_1, \
                                         tg_yzz_yyyyyyy_0, tg_yzz_yyyyyyy_1, tg_yzz_yyyyyyz_0, tg_yzz_yyyyyyz_1, \
                                         tg_yzz_yyyyyz_1, tg_yzz_yyyyyzz_0, tg_yzz_yyyyyzz_1, tg_yzz_yyyyzz_1, \
                                         tg_yzz_yyyyzzz_0, tg_yzz_yyyyzzz_1, tg_yzz_yyyzzz_1, tg_yzz_yyyzzzz_0, \
                                         tg_yzz_yyyzzzz_1, tg_yzz_yyzzzz_1, tg_yzz_yyzzzzz_0, tg_yzz_yyzzzzz_1, \
                                         tg_yzz_yzzzzz_1, tg_yzz_yzzzzzz_0, tg_yzz_yzzzzzz_1, tg_yzz_zzzzzz_1, \
                                         tg_yzz_zzzzzzz_0, tg_yzz_zzzzzzz_1, tg_yzzz_xxxxxxx_0, tg_yzzz_xxxxxxy_0, \
                                         tg_yzzz_xxxxxxz_0, tg_yzzz_xxxxxyy_0, tg_yzzz_xxxxxyz_0, tg_yzzz_xxxxxzz_0, \
                                         tg_yzzz_xxxxyyy_0, tg_yzzz_xxxxyyz_0, tg_yzzz_xxxxyzz_0, tg_yzzz_xxxxzzz_0, \
                                         tg_yzzz_xxxyyyy_0, tg_yzzz_xxxyyyz_0, tg_yzzz_xxxyyzz_0, tg_yzzz_xxxyzzz_0, \
                                         tg_yzzz_xxxzzzz_0, tg_yzzz_xxyyyyy_0, tg_yzzz_xxyyyyz_0, tg_yzzz_xxyyyzz_0, \
                                         tg_yzzz_xxyyzzz_0, tg_yzzz_xxyzzzz_0, tg_yzzz_xxzzzzz_0, tg_yzzz_xyyyyyy_0, \
                                         tg_yzzz_xyyyyyz_0, tg_yzzz_xyyyyzz_0, tg_yzzz_xyyyzzz_0, tg_yzzz_xyyzzzz_0, \
                                         tg_yzzz_xyzzzzz_0, tg_yzzz_xzzzzzz_0, tg_yzzz_yyyyyyy_0, tg_yzzz_yyyyyyz_0, \
                                         tg_yzzz_yyyyyzz_0, tg_yzzz_yyyyzzz_0, tg_yzzz_yyyzzzz_0, tg_yzzz_yyzzzzz_0, \
                                         tg_yzzz_yzzzzzz_0, tg_yzzz_zzzzzzz_0, tg_zz_xxxxxxx_0, tg_zz_xxxxxxx_1, \
                                         tg_zz_xxxxxxy_0, tg_zz_xxxxxxy_1, tg_zz_xxxxxxz_0, tg_zz_xxxxxxz_1, tg_zz_xxxxxyy_0, \
                                         tg_zz_xxxxxyy_1, tg_zz_xxxxxyz_0, tg_zz_xxxxxyz_1, tg_zz_xxxxxzz_0, tg_zz_xxxxxzz_1, \
                                         tg_zz_xxxxyyy_0, tg_zz_xxxxyyy_1, tg_zz_xxxxyyz_0, tg_zz_xxxxyyz_1, tg_zz_xxxxyzz_0, \
                                         tg_zz_xxxxyzz_1, tg_zz_xxxxzzz_0, tg_zz_xxxxzzz_1, tg_zz_xxxyyyy_0, tg_zz_xxxyyyy_1, \
                                         tg_zz_xxxyyyz_0, tg_zz_xxxyyyz_1, tg_zz_xxxyyzz_0, tg_zz_xxxyyzz_1, tg_zz_xxxyzzz_0, \
                                         tg_zz_xxxyzzz_1, tg_zz_xxxzzzz_0, tg_zz_xxxzzzz_1, tg_zz_xxyyyyy_0, tg_zz_xxyyyyy_1, \
                                         tg_zz_xxyyyyz_0, tg_zz_xxyyyyz_1, tg_zz_xxyyyzz_0, tg_zz_xxyyyzz_1, tg_zz_xxyyzzz_0, \
                                         tg_zz_xxyyzzz_1, tg_zz_xxyzzzz_0, tg_zz_xxyzzzz_1, tg_zz_xxzzzzz_0, tg_zz_xxzzzzz_1, \
                                         tg_zz_xyyyyyy_0, tg_zz_xyyyyyy_1, tg_zz_xyyyyyz_0, tg_zz_xyyyyyz_1, tg_zz_xyyyyzz_0, \
                                         tg_zz_xyyyyzz_1, tg_zz_xyyyzzz_0, tg_zz_xyyyzzz_1, tg_zz_xyyzzzz_0, tg_zz_xyyzzzz_1, \
                                         tg_zz_xyzzzzz_0, tg_zz_xyzzzzz_1, tg_zz_xzzzzzz_0, tg_zz_xzzzzzz_1, tg_zz_yyyyyyy_0, \
                                         tg_zz_yyyyyyy_1, tg_zz_yyyyyyz_0, tg_zz_yyyyyyz_1, tg_zz_yyyyyzz_0, tg_zz_yyyyyzz_1, \
                                         tg_zz_yyyyzzz_0, tg_zz_yyyyzzz_1, tg_zz_yyyzzzz_0, tg_zz_yyyzzzz_1, tg_zz_yyzzzzz_0, \
                                         tg_zz_yyzzzzz_1, tg_zz_yzzzzzz_0, tg_zz_yzzzzzz_1, tg_zz_zzzzzzz_0, tg_zz_zzzzzzz_1, \
                                         tg_zzz_xxxxxx_1, tg_zzz_xxxxxxx_0, tg_zzz_xxxxxxx_1, tg_zzz_xxxxxxy_0, \
                                         tg_zzz_xxxxxxy_1, tg_zzz_xxxxxxz_0, tg_zzz_xxxxxxz_1, tg_zzz_xxxxxy_1, \
                                         tg_zzz_xxxxxyy_0, tg_zzz_xxxxxyy_1, tg_zzz_xxxxxyz_0, tg_zzz_xxxxxyz_1, \
                                         tg_zzz_xxxxxz_1, tg_zzz_xxxxxzz_0, tg_zzz_xxxxxzz_1, tg_zzz_xxxxyy_1, \
                                         tg_zzz_xxxxyyy_0, tg_zzz_xxxxyyy_1, tg_zzz_xxxxyyz_0, tg_zzz_xxxxyyz_1, \
                                         tg_zzz_xxxxyz_1, tg_zzz_xxxxyzz_0, tg_zzz_xxxxyzz_1, tg_zzz_xxxxzz_1, \
                                         tg_zzz_xxxxzzz_0, tg_zzz_xxxxzzz_1, tg_zzz_xxxyyy_1, tg_zzz_xxxyyyy_0, \
                                         tg_zzz_xxxyyyy_1, tg_zzz_xxxyyyz_0, tg_zzz_xxxyyyz_1, tg_zzz_xxxyyz_1, \
                                         tg_zzz_xxxyyzz_0, tg_zzz_xxxyyzz_1, tg_zzz_xxxyzz_1, tg_zzz_xxxyzzz_0, \
                                         tg_zzz_xxxyzzz_1, tg_zzz_xxxzzz_1, tg_zzz_xxxzzzz_0, tg_zzz_xxxzzzz_1, \
                                         tg_zzz_xxyyyy_1, tg_zzz_xxyyyyy_0, tg_zzz_xxyyyyy_1, tg_zzz_xxyyyyz_0, \
                                         tg_zzz_xxyyyyz_1, tg_zzz_xxyyyz_1, tg_zzz_xxyyyzz_0, tg_zzz_xxyyyzz_1, \
                                         tg_zzz_xxyyzz_1, tg_zzz_xxyyzzz_0, tg_zzz_xxyyzzz_1, tg_zzz_xxyzzz_1, \
                                         tg_zzz_xxyzzzz_0, tg_zzz_xxyzzzz_1, tg_zzz_xxzzzz_1, tg_zzz_xxzzzzz_0, \
                                         tg_zzz_xxzzzzz_1, tg_zzz_xyyyyy_1, tg_zzz_xyyyyyy_0, tg_zzz_xyyyyyy_1, \
                                         tg_zzz_xyyyyyz_0, tg_zzz_xyyyyyz_1, tg_zzz_xyyyyz_1, tg_zzz_xyyyyzz_0, \
                                         tg_zzz_xyyyyzz_1, tg_zzz_xyyyzz_1, tg_zzz_xyyyzzz_0, tg_zzz_xyyyzzz_1, \
                                         tg_zzz_xyyzzz_1, tg_zzz_xyyzzzz_0, tg_zzz_xyyzzzz_1, tg_zzz_xyzzzz_1, \
                                         tg_zzz_xyzzzzz_0, tg_zzz_xyzzzzz_1, tg_zzz_xzzzzz_1, tg_zzz_xzzzzzz_0, \
                                         tg_zzz_xzzzzzz_1, tg_zzz_yyyyyy_1, tg_zzz_yyyyyyy_0, tg_zzz_yyyyyyy_1, \
                                         tg_zzz_yyyyyyz_0, tg_zzz_yyyyyyz_1, tg_zzz_yyyyyz_1, tg_zzz_yyyyyzz_0, \
                                         tg_zzz_yyyyyzz_1, tg_zzz_yyyyzz_1, tg_zzz_yyyyzzz_0, tg_zzz_yyyyzzz_1, \
                                         tg_zzz_yyyzzz_1, tg_zzz_yyyzzzz_0, tg_zzz_yyyzzzz_1, tg_zzz_yyzzzz_1, \
                                         tg_zzz_yyzzzzz_0, tg_zzz_yyzzzzz_1, tg_zzz_yzzzzz_1, tg_zzz_yzzzzzz_0, \
                                         tg_zzz_yzzzzzz_1, tg_zzz_zzzzzz_1, tg_zzz_zzzzzzz_0, tg_zzz_zzzzzzz_1, \
                                         tg_zzzz_xxxxxxx_0, tg_zzzz_xxxxxxy_0, tg_zzzz_xxxxxxz_0, tg_zzzz_xxxxxyy_0, \
                                         tg_zzzz_xxxxxyz_0, tg_zzzz_xxxxxzz_0, tg_zzzz_xxxxyyy_0, tg_zzzz_xxxxyyz_0, \
                                         tg_zzzz_xxxxyzz_0, tg_zzzz_xxxxzzz_0, tg_zzzz_xxxyyyy_0, tg_zzzz_xxxyyyz_0, \
                                         tg_zzzz_xxxyyzz_0, tg_zzzz_xxxyzzz_0, tg_zzzz_xxxzzzz_0, tg_zzzz_xxyyyyy_0, \
                                         tg_zzzz_xxyyyyz_0, tg_zzzz_xxyyyzz_0, tg_zzzz_xxyyzzz_0, tg_zzzz_xxyzzzz_0, \
                                         tg_zzzz_xxzzzzz_0, tg_zzzz_xyyyyyy_0, tg_zzzz_xyyyyyz_0, tg_zzzz_xyyyyzz_0, \
                                         tg_zzzz_xyyyzzz_0, tg_zzzz_xyyzzzz_0, tg_zzzz_xyzzzzz_0, tg_zzzz_xzzzzzz_0, \
                                         tg_zzzz_yyyyyyy_0, tg_zzzz_yyyyyyz_0, tg_zzzz_yyyyyzz_0, tg_zzzz_yyyyzzz_0, \
                                         tg_zzzz_yyyzzzz_0, tg_zzzz_yyzzzzz_0, tg_zzzz_yzzzzzz_0, tg_zzzz_zzzzzzz_0, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_y[j]; 

                    tg_yyzz_xxyyzzz_0[j] = pb_y * tg_yzz_xxyyzzz_0[j] + fr * tg_yzz_xxyyzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxyyzzz_0[j] - tg_zz_xxyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yzz_xxyzzz_1[j];

                    tg_yyzz_xxyzzzz_0[j] = pb_y * tg_yzz_xxyzzzz_0[j] + fr * tg_yzz_xxyzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxyzzzz_0[j] - tg_zz_xxyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzz_xxzzzz_1[j];

                    tg_yyzz_xxzzzzz_0[j] = pb_y * tg_yzz_xxzzzzz_0[j] + fr * tg_yzz_xxzzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxzzzzz_0[j] - tg_zz_xxzzzzz_1[j] * fl1_fza);

                    tg_yyzz_xyyyyyy_0[j] = pb_y * tg_yzz_xyyyyyy_0[j] + fr * tg_yzz_xyyyyyy_1[j] + 0.5 * fl1_fx * (tg_zz_xyyyyyy_0[j] - tg_zz_xyyyyyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yzz_xyyyyy_1[j];

                    tg_yyzz_xyyyyyz_0[j] = pb_y * tg_yzz_xyyyyyz_0[j] + fr * tg_yzz_xyyyyyz_1[j] + 0.5 * fl1_fx * (tg_zz_xyyyyyz_0[j] - tg_zz_xyyyyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yzz_xyyyyz_1[j];

                    tg_yyzz_xyyyyzz_0[j] = pb_y * tg_yzz_xyyyyzz_0[j] + fr * tg_yzz_xyyyyzz_1[j] + 0.5 * fl1_fx * (tg_zz_xyyyyzz_0[j] - tg_zz_xyyyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzz_xyyyzz_1[j];

                    tg_yyzz_xyyyzzz_0[j] = pb_y * tg_yzz_xyyyzzz_0[j] + fr * tg_yzz_xyyyzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xyyyzzz_0[j] - tg_zz_xyyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzz_xyyzzz_1[j];

                    tg_yyzz_xyyzzzz_0[j] = pb_y * tg_yzz_xyyzzzz_0[j] + fr * tg_yzz_xyyzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xyyzzzz_0[j] - tg_zz_xyyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yzz_xyzzzz_1[j];

                    tg_yyzz_xyzzzzz_0[j] = pb_y * tg_yzz_xyzzzzz_0[j] + fr * tg_yzz_xyzzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xyzzzzz_0[j] - tg_zz_xyzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzz_xzzzzz_1[j];

                    tg_yyzz_xzzzzzz_0[j] = pb_y * tg_yzz_xzzzzzz_0[j] + fr * tg_yzz_xzzzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xzzzzzz_0[j] - tg_zz_xzzzzzz_1[j] * fl1_fza);

                    tg_yyzz_yyyyyyy_0[j] = pb_y * tg_yzz_yyyyyyy_0[j] + fr * tg_yzz_yyyyyyy_1[j] + 0.5 * fl1_fx * (tg_zz_yyyyyyy_0[j] - tg_zz_yyyyyyy_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_yzz_yyyyyy_1[j];

                    tg_yyzz_yyyyyyz_0[j] = pb_y * tg_yzz_yyyyyyz_0[j] + fr * tg_yzz_yyyyyyz_1[j] + 0.5 * fl1_fx * (tg_zz_yyyyyyz_0[j] - tg_zz_yyyyyyz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yzz_yyyyyz_1[j];

                    tg_yyzz_yyyyyzz_0[j] = pb_y * tg_yzz_yyyyyzz_0[j] + fr * tg_yzz_yyyyyzz_1[j] + 0.5 * fl1_fx * (tg_zz_yyyyyzz_0[j] - tg_zz_yyyyyzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yzz_yyyyzz_1[j];

                    tg_yyzz_yyyyzzz_0[j] = pb_y * tg_yzz_yyyyzzz_0[j] + fr * tg_yzz_yyyyzzz_1[j] + 0.5 * fl1_fx * (tg_zz_yyyyzzz_0[j] - tg_zz_yyyyzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzz_yyyzzz_1[j];

                    tg_yyzz_yyyzzzz_0[j] = pb_y * tg_yzz_yyyzzzz_0[j] + fr * tg_yzz_yyyzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_yyyzzzz_0[j] - tg_zz_yyyzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzz_yyzzzz_1[j];

                    tg_yyzz_yyzzzzz_0[j] = pb_y * tg_yzz_yyzzzzz_0[j] + fr * tg_yzz_yyzzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_yyzzzzz_0[j] - tg_zz_yyzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yzz_yzzzzz_1[j];

                    tg_yyzz_yzzzzzz_0[j] = pb_y * tg_yzz_yzzzzzz_0[j] + fr * tg_yzz_yzzzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_yzzzzzz_0[j] - tg_zz_yzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzz_zzzzzz_1[j];

                    tg_yyzz_zzzzzzz_0[j] = pb_y * tg_yzz_zzzzzzz_0[j] + fr * tg_yzz_zzzzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_zzzzzzz_0[j] - tg_zz_zzzzzzz_1[j] * fl1_fza);

                    tg_yzzz_xxxxxxx_0[j] = pb_y * tg_zzz_xxxxxxx_0[j] + fr * tg_zzz_xxxxxxx_1[j];

                    tg_yzzz_xxxxxxy_0[j] = pb_y * tg_zzz_xxxxxxy_0[j] + fr * tg_zzz_xxxxxxy_1[j] + 0.5 * fl1_fxn * tg_zzz_xxxxxx_1[j];

                    tg_yzzz_xxxxxxz_0[j] = pb_y * tg_zzz_xxxxxxz_0[j] + fr * tg_zzz_xxxxxxz_1[j];

                    tg_yzzz_xxxxxyy_0[j] = pb_y * tg_zzz_xxxxxyy_0[j] + fr * tg_zzz_xxxxxyy_1[j] + fl1_fxn * tg_zzz_xxxxxy_1[j];

                    tg_yzzz_xxxxxyz_0[j] = pb_y * tg_zzz_xxxxxyz_0[j] + fr * tg_zzz_xxxxxyz_1[j] + 0.5 * fl1_fxn * tg_zzz_xxxxxz_1[j];

                    tg_yzzz_xxxxxzz_0[j] = pb_y * tg_zzz_xxxxxzz_0[j] + fr * tg_zzz_xxxxxzz_1[j];

                    tg_yzzz_xxxxyyy_0[j] = pb_y * tg_zzz_xxxxyyy_0[j] + fr * tg_zzz_xxxxyyy_1[j] + 1.5 * fl1_fxn * tg_zzz_xxxxyy_1[j];

                    tg_yzzz_xxxxyyz_0[j] = pb_y * tg_zzz_xxxxyyz_0[j] + fr * tg_zzz_xxxxyyz_1[j] + fl1_fxn * tg_zzz_xxxxyz_1[j];

                    tg_yzzz_xxxxyzz_0[j] = pb_y * tg_zzz_xxxxyzz_0[j] + fr * tg_zzz_xxxxyzz_1[j] + 0.5 * fl1_fxn * tg_zzz_xxxxzz_1[j];

                    tg_yzzz_xxxxzzz_0[j] = pb_y * tg_zzz_xxxxzzz_0[j] + fr * tg_zzz_xxxxzzz_1[j];

                    tg_yzzz_xxxyyyy_0[j] = pb_y * tg_zzz_xxxyyyy_0[j] + fr * tg_zzz_xxxyyyy_1[j] + 2.0 * fl1_fxn * tg_zzz_xxxyyy_1[j];

                    tg_yzzz_xxxyyyz_0[j] = pb_y * tg_zzz_xxxyyyz_0[j] + fr * tg_zzz_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_zzz_xxxyyz_1[j];

                    tg_yzzz_xxxyyzz_0[j] = pb_y * tg_zzz_xxxyyzz_0[j] + fr * tg_zzz_xxxyyzz_1[j] + fl1_fxn * tg_zzz_xxxyzz_1[j];

                    tg_yzzz_xxxyzzz_0[j] = pb_y * tg_zzz_xxxyzzz_0[j] + fr * tg_zzz_xxxyzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_xxxzzz_1[j];

                    tg_yzzz_xxxzzzz_0[j] = pb_y * tg_zzz_xxxzzzz_0[j] + fr * tg_zzz_xxxzzzz_1[j];

                    tg_yzzz_xxyyyyy_0[j] = pb_y * tg_zzz_xxyyyyy_0[j] + fr * tg_zzz_xxyyyyy_1[j] + 2.5 * fl1_fxn * tg_zzz_xxyyyy_1[j];

                    tg_yzzz_xxyyyyz_0[j] = pb_y * tg_zzz_xxyyyyz_0[j] + fr * tg_zzz_xxyyyyz_1[j] + 2.0 * fl1_fxn * tg_zzz_xxyyyz_1[j];

                    tg_yzzz_xxyyyzz_0[j] = pb_y * tg_zzz_xxyyyzz_0[j] + fr * tg_zzz_xxyyyzz_1[j] + 1.5 * fl1_fxn * tg_zzz_xxyyzz_1[j];

                    tg_yzzz_xxyyzzz_0[j] = pb_y * tg_zzz_xxyyzzz_0[j] + fr * tg_zzz_xxyyzzz_1[j] + fl1_fxn * tg_zzz_xxyzzz_1[j];

                    tg_yzzz_xxyzzzz_0[j] = pb_y * tg_zzz_xxyzzzz_0[j] + fr * tg_zzz_xxyzzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_xxzzzz_1[j];

                    tg_yzzz_xxzzzzz_0[j] = pb_y * tg_zzz_xxzzzzz_0[j] + fr * tg_zzz_xxzzzzz_1[j];

                    tg_yzzz_xyyyyyy_0[j] = pb_y * tg_zzz_xyyyyyy_0[j] + fr * tg_zzz_xyyyyyy_1[j] + 3.0 * fl1_fxn * tg_zzz_xyyyyy_1[j];

                    tg_yzzz_xyyyyyz_0[j] = pb_y * tg_zzz_xyyyyyz_0[j] + fr * tg_zzz_xyyyyyz_1[j] + 2.5 * fl1_fxn * tg_zzz_xyyyyz_1[j];

                    tg_yzzz_xyyyyzz_0[j] = pb_y * tg_zzz_xyyyyzz_0[j] + fr * tg_zzz_xyyyyzz_1[j] + 2.0 * fl1_fxn * tg_zzz_xyyyzz_1[j];

                    tg_yzzz_xyyyzzz_0[j] = pb_y * tg_zzz_xyyyzzz_0[j] + fr * tg_zzz_xyyyzzz_1[j] + 1.5 * fl1_fxn * tg_zzz_xyyzzz_1[j];

                    tg_yzzz_xyyzzzz_0[j] = pb_y * tg_zzz_xyyzzzz_0[j] + fr * tg_zzz_xyyzzzz_1[j] + fl1_fxn * tg_zzz_xyzzzz_1[j];

                    tg_yzzz_xyzzzzz_0[j] = pb_y * tg_zzz_xyzzzzz_0[j] + fr * tg_zzz_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_xzzzzz_1[j];

                    tg_yzzz_xzzzzzz_0[j] = pb_y * tg_zzz_xzzzzzz_0[j] + fr * tg_zzz_xzzzzzz_1[j];

                    tg_yzzz_yyyyyyy_0[j] = pb_y * tg_zzz_yyyyyyy_0[j] + fr * tg_zzz_yyyyyyy_1[j] + 3.5 * fl1_fxn * tg_zzz_yyyyyy_1[j];

                    tg_yzzz_yyyyyyz_0[j] = pb_y * tg_zzz_yyyyyyz_0[j] + fr * tg_zzz_yyyyyyz_1[j] + 3.0 * fl1_fxn * tg_zzz_yyyyyz_1[j];

                    tg_yzzz_yyyyyzz_0[j] = pb_y * tg_zzz_yyyyyzz_0[j] + fr * tg_zzz_yyyyyzz_1[j] + 2.5 * fl1_fxn * tg_zzz_yyyyzz_1[j];

                    tg_yzzz_yyyyzzz_0[j] = pb_y * tg_zzz_yyyyzzz_0[j] + fr * tg_zzz_yyyyzzz_1[j] + 2.0 * fl1_fxn * tg_zzz_yyyzzz_1[j];

                    tg_yzzz_yyyzzzz_0[j] = pb_y * tg_zzz_yyyzzzz_0[j] + fr * tg_zzz_yyyzzzz_1[j] + 1.5 * fl1_fxn * tg_zzz_yyzzzz_1[j];

                    tg_yzzz_yyzzzzz_0[j] = pb_y * tg_zzz_yyzzzzz_0[j] + fr * tg_zzz_yyzzzzz_1[j] + fl1_fxn * tg_zzz_yzzzzz_1[j];

                    tg_yzzz_yzzzzzz_0[j] = pb_y * tg_zzz_yzzzzzz_0[j] + fr * tg_zzz_yzzzzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_zzzzzz_1[j];

                    tg_yzzz_zzzzzzz_0[j] = pb_y * tg_zzz_zzzzzzz_0[j] + fr * tg_zzz_zzzzzzz_1[j];

                    fr = wp_z[j]; 

                    tg_zzzz_xxxxxxx_0[j] = pb_z * tg_zzz_xxxxxxx_0[j] + fr * tg_zzz_xxxxxxx_1[j] + 1.5 * fl1_fx * (tg_zz_xxxxxxx_0[j] - tg_zz_xxxxxxx_1[j] * fl1_fza);

                    tg_zzzz_xxxxxxy_0[j] = pb_z * tg_zzz_xxxxxxy_0[j] + fr * tg_zzz_xxxxxxy_1[j] + 1.5 * fl1_fx * (tg_zz_xxxxxxy_0[j] - tg_zz_xxxxxxy_1[j] * fl1_fza);

                    tg_zzzz_xxxxxxz_0[j] = pb_z * tg_zzz_xxxxxxz_0[j] + fr * tg_zzz_xxxxxxz_1[j] + 1.5 * fl1_fx * (tg_zz_xxxxxxz_0[j] - tg_zz_xxxxxxz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzz_xxxxxx_1[j];

                    tg_zzzz_xxxxxyy_0[j] = pb_z * tg_zzz_xxxxxyy_0[j] + fr * tg_zzz_xxxxxyy_1[j] + 1.5 * fl1_fx * (tg_zz_xxxxxyy_0[j] - tg_zz_xxxxxyy_1[j] * fl1_fza);

                    tg_zzzz_xxxxxyz_0[j] = pb_z * tg_zzz_xxxxxyz_0[j] + fr * tg_zzz_xxxxxyz_1[j] + 1.5 * fl1_fx * (tg_zz_xxxxxyz_0[j] - tg_zz_xxxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzz_xxxxxy_1[j];

                    tg_zzzz_xxxxxzz_0[j] = pb_z * tg_zzz_xxxxxzz_0[j] + fr * tg_zzz_xxxxxzz_1[j] + 1.5 * fl1_fx * (tg_zz_xxxxxzz_0[j] - tg_zz_xxxxxzz_1[j] * fl1_fza) + fl1_fxn * tg_zzz_xxxxxz_1[j];

                    tg_zzzz_xxxxyyy_0[j] = pb_z * tg_zzz_xxxxyyy_0[j] + fr * tg_zzz_xxxxyyy_1[j] + 1.5 * fl1_fx * (tg_zz_xxxxyyy_0[j] - tg_zz_xxxxyyy_1[j] * fl1_fza);

                    tg_zzzz_xxxxyyz_0[j] = pb_z * tg_zzz_xxxxyyz_0[j] + fr * tg_zzz_xxxxyyz_1[j] + 1.5 * fl1_fx * (tg_zz_xxxxyyz_0[j] - tg_zz_xxxxyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzz_xxxxyy_1[j];

                    tg_zzzz_xxxxyzz_0[j] = pb_z * tg_zzz_xxxxyzz_0[j] + fr * tg_zzz_xxxxyzz_1[j] + 1.5 * fl1_fx * (tg_zz_xxxxyzz_0[j] - tg_zz_xxxxyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzz_xxxxyz_1[j];

                    tg_zzzz_xxxxzzz_0[j] = pb_z * tg_zzz_xxxxzzz_0[j] + fr * tg_zzz_xxxxzzz_1[j] + 1.5 * fl1_fx * (tg_zz_xxxxzzz_0[j] - tg_zz_xxxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzz_xxxxzz_1[j];

                    tg_zzzz_xxxyyyy_0[j] = pb_z * tg_zzz_xxxyyyy_0[j] + fr * tg_zzz_xxxyyyy_1[j] + 1.5 * fl1_fx * (tg_zz_xxxyyyy_0[j] - tg_zz_xxxyyyy_1[j] * fl1_fza);

                    tg_zzzz_xxxyyyz_0[j] = pb_z * tg_zzz_xxxyyyz_0[j] + fr * tg_zzz_xxxyyyz_1[j] + 1.5 * fl1_fx * (tg_zz_xxxyyyz_0[j] - tg_zz_xxxyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzz_xxxyyy_1[j];

                    tg_zzzz_xxxyyzz_0[j] = pb_z * tg_zzz_xxxyyzz_0[j] + fr * tg_zzz_xxxyyzz_1[j] + 1.5 * fl1_fx * (tg_zz_xxxyyzz_0[j] - tg_zz_xxxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzz_xxxyyz_1[j];

                    tg_zzzz_xxxyzzz_0[j] = pb_z * tg_zzz_xxxyzzz_0[j] + fr * tg_zzz_xxxyzzz_1[j] + 1.5 * fl1_fx * (tg_zz_xxxyzzz_0[j] - tg_zz_xxxyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzz_xxxyzz_1[j];

                    tg_zzzz_xxxzzzz_0[j] = pb_z * tg_zzz_xxxzzzz_0[j] + fr * tg_zzz_xxxzzzz_1[j] + 1.5 * fl1_fx * (tg_zz_xxxzzzz_0[j] - tg_zz_xxxzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzz_xxxzzz_1[j];

                    tg_zzzz_xxyyyyy_0[j] = pb_z * tg_zzz_xxyyyyy_0[j] + fr * tg_zzz_xxyyyyy_1[j] + 1.5 * fl1_fx * (tg_zz_xxyyyyy_0[j] - tg_zz_xxyyyyy_1[j] * fl1_fza);

                    tg_zzzz_xxyyyyz_0[j] = pb_z * tg_zzz_xxyyyyz_0[j] + fr * tg_zzz_xxyyyyz_1[j] + 1.5 * fl1_fx * (tg_zz_xxyyyyz_0[j] - tg_zz_xxyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzz_xxyyyy_1[j];

                    tg_zzzz_xxyyyzz_0[j] = pb_z * tg_zzz_xxyyyzz_0[j] + fr * tg_zzz_xxyyyzz_1[j] + 1.5 * fl1_fx * (tg_zz_xxyyyzz_0[j] - tg_zz_xxyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzz_xxyyyz_1[j];

                    tg_zzzz_xxyyzzz_0[j] = pb_z * tg_zzz_xxyyzzz_0[j] + fr * tg_zzz_xxyyzzz_1[j] + 1.5 * fl1_fx * (tg_zz_xxyyzzz_0[j] - tg_zz_xxyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzz_xxyyzz_1[j];

                    tg_zzzz_xxyzzzz_0[j] = pb_z * tg_zzz_xxyzzzz_0[j] + fr * tg_zzz_xxyzzzz_1[j] + 1.5 * fl1_fx * (tg_zz_xxyzzzz_0[j] - tg_zz_xxyzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzz_xxyzzz_1[j];

                    tg_zzzz_xxzzzzz_0[j] = pb_z * tg_zzz_xxzzzzz_0[j] + fr * tg_zzz_xxzzzzz_1[j] + 1.5 * fl1_fx * (tg_zz_xxzzzzz_0[j] - tg_zz_xxzzzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_zzz_xxzzzz_1[j];

                    tg_zzzz_xyyyyyy_0[j] = pb_z * tg_zzz_xyyyyyy_0[j] + fr * tg_zzz_xyyyyyy_1[j] + 1.5 * fl1_fx * (tg_zz_xyyyyyy_0[j] - tg_zz_xyyyyyy_1[j] * fl1_fza);

                    tg_zzzz_xyyyyyz_0[j] = pb_z * tg_zzz_xyyyyyz_0[j] + fr * tg_zzz_xyyyyyz_1[j] + 1.5 * fl1_fx * (tg_zz_xyyyyyz_0[j] - tg_zz_xyyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzz_xyyyyy_1[j];

                    tg_zzzz_xyyyyzz_0[j] = pb_z * tg_zzz_xyyyyzz_0[j] + fr * tg_zzz_xyyyyzz_1[j] + 1.5 * fl1_fx * (tg_zz_xyyyyzz_0[j] - tg_zz_xyyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzz_xyyyyz_1[j];

                    tg_zzzz_xyyyzzz_0[j] = pb_z * tg_zzz_xyyyzzz_0[j] + fr * tg_zzz_xyyyzzz_1[j] + 1.5 * fl1_fx * (tg_zz_xyyyzzz_0[j] - tg_zz_xyyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzz_xyyyzz_1[j];

                    tg_zzzz_xyyzzzz_0[j] = pb_z * tg_zzz_xyyzzzz_0[j] + fr * tg_zzz_xyyzzzz_1[j] + 1.5 * fl1_fx * (tg_zz_xyyzzzz_0[j] - tg_zz_xyyzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzz_xyyzzz_1[j];

                    tg_zzzz_xyzzzzz_0[j] = pb_z * tg_zzz_xyzzzzz_0[j] + fr * tg_zzz_xyzzzzz_1[j] + 1.5 * fl1_fx * (tg_zz_xyzzzzz_0[j] - tg_zz_xyzzzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_zzz_xyzzzz_1[j];

                    tg_zzzz_xzzzzzz_0[j] = pb_z * tg_zzz_xzzzzzz_0[j] + fr * tg_zzz_xzzzzzz_1[j] + 1.5 * fl1_fx * (tg_zz_xzzzzzz_0[j] - tg_zz_xzzzzzz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_zzz_xzzzzz_1[j];

                    tg_zzzz_yyyyyyy_0[j] = pb_z * tg_zzz_yyyyyyy_0[j] + fr * tg_zzz_yyyyyyy_1[j] + 1.5 * fl1_fx * (tg_zz_yyyyyyy_0[j] - tg_zz_yyyyyyy_1[j] * fl1_fza);

                    tg_zzzz_yyyyyyz_0[j] = pb_z * tg_zzz_yyyyyyz_0[j] + fr * tg_zzz_yyyyyyz_1[j] + 1.5 * fl1_fx * (tg_zz_yyyyyyz_0[j] - tg_zz_yyyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzz_yyyyyy_1[j];

                    tg_zzzz_yyyyyzz_0[j] = pb_z * tg_zzz_yyyyyzz_0[j] + fr * tg_zzz_yyyyyzz_1[j] + 1.5 * fl1_fx * (tg_zz_yyyyyzz_0[j] - tg_zz_yyyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzz_yyyyyz_1[j];

                    tg_zzzz_yyyyzzz_0[j] = pb_z * tg_zzz_yyyyzzz_0[j] + fr * tg_zzz_yyyyzzz_1[j] + 1.5 * fl1_fx * (tg_zz_yyyyzzz_0[j] - tg_zz_yyyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzz_yyyyzz_1[j];

                    tg_zzzz_yyyzzzz_0[j] = pb_z * tg_zzz_yyyzzzz_0[j] + fr * tg_zzz_yyyzzzz_1[j] + 1.5 * fl1_fx * (tg_zz_yyyzzzz_0[j] - tg_zz_yyyzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzz_yyyzzz_1[j];

                    tg_zzzz_yyzzzzz_0[j] = pb_z * tg_zzz_yyzzzzz_0[j] + fr * tg_zzz_yyzzzzz_1[j] + 1.5 * fl1_fx * (tg_zz_yyzzzzz_0[j] - tg_zz_yyzzzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_zzz_yyzzzz_1[j];

                    tg_zzzz_yzzzzzz_0[j] = pb_z * tg_zzz_yzzzzzz_0[j] + fr * tg_zzz_yzzzzzz_1[j] + 1.5 * fl1_fx * (tg_zz_yzzzzzz_0[j] - tg_zz_yzzzzzz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_zzz_yzzzzz_1[j];

                    tg_zzzz_zzzzzzz_0[j] = pb_z * tg_zzz_zzzzzzz_0[j] + fr * tg_zzz_zzzzzzz_1[j] + 1.5 * fl1_fx * (tg_zz_zzzzzzz_0[j] - tg_zz_zzzzzzz_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_zzz_zzzzzz_1[j];
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

