//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectronRepulsionRecFuncForFK.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSFSK(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSFSK_0_90(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSFSK_90_180(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 

        erirecfunc::compElectronRepulsionForSFSK_180_270(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSFSK_270_360(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 
    }

    void
    compElectronRepulsionForSFSK_0_90(      CMemBlock2D<double>& primBuffer,
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
                                             {3, -1, -1, -1},
                                             {7, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_3_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_3_7_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_2_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_2_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_1_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_1_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_2_6_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {6, -1, -1, -1}, 
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

                auto tg_xx_xxxxxxx_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx); 

                auto tg_xx_xxxxxxy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 1); 

                auto tg_xx_xxxxxxz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 2); 

                auto tg_xx_xxxxxyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 3); 

                auto tg_xx_xxxxxyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 4); 

                auto tg_xx_xxxxxzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 5); 

                auto tg_xx_xxxxyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 6); 

                auto tg_xx_xxxxyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 7); 

                auto tg_xx_xxxxyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 8); 

                auto tg_xx_xxxxzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 9); 

                auto tg_xx_xxxyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 10); 

                auto tg_xx_xxxyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 11); 

                auto tg_xx_xxxyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 12); 

                auto tg_xx_xxxyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 13); 

                auto tg_xx_xxxzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 14); 

                auto tg_xx_xxyyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 15); 

                auto tg_xx_xxyyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 16); 

                auto tg_xx_xxyyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 17); 

                auto tg_xx_xxyyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 18); 

                auto tg_xx_xxyzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 19); 

                auto tg_xx_xxzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 20); 

                auto tg_xx_xyyyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 21); 

                auto tg_xx_xyyyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 22); 

                auto tg_xx_xyyyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 23); 

                auto tg_xx_xyyyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 24); 

                auto tg_xx_xyyzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 25); 

                auto tg_xx_xyzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 26); 

                auto tg_xx_xzzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 27); 

                auto tg_xx_yyyyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 28); 

                auto tg_xx_yyyyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 29); 

                auto tg_xx_yyyyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 30); 

                auto tg_xx_yyyyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 31); 

                auto tg_xx_yyyzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 32); 

                auto tg_xx_yyzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 33); 

                auto tg_xx_yzzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 34); 

                auto tg_xx_zzzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 35); 

                auto tg_xy_xxxxxxx_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 36); 

                auto tg_xy_xxxxxxy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 37); 

                auto tg_xy_xxxxxxz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 38); 

                auto tg_xy_xxxxxyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 39); 

                auto tg_xy_xxxxxyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 40); 

                auto tg_xy_xxxxxzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 41); 

                auto tg_xy_xxxxyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 42); 

                auto tg_xy_xxxxyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 43); 

                auto tg_xy_xxxxyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 44); 

                auto tg_xy_xxxxzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 45); 

                auto tg_xy_xxxyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 46); 

                auto tg_xy_xxxyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 47); 

                auto tg_xy_xxxyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 48); 

                auto tg_xy_xxxyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 49); 

                auto tg_xy_xxxzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 50); 

                auto tg_xy_xxyyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 51); 

                auto tg_xy_xxyyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 52); 

                auto tg_xy_xxyyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 53); 

                auto tg_xy_xxyyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 54); 

                auto tg_xy_xxyzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 55); 

                auto tg_xy_xxzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 56); 

                auto tg_xy_xyyyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 57); 

                auto tg_xy_xyyyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 58); 

                auto tg_xy_xyyyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 59); 

                auto tg_xy_xyyyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 60); 

                auto tg_xy_xyyzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 61); 

                auto tg_xy_xyzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 62); 

                auto tg_xy_xzzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 63); 

                auto tg_xy_yyyyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 64); 

                auto tg_xy_yyyyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 65); 

                auto tg_xy_yyyyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 66); 

                auto tg_xy_yyyyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 67); 

                auto tg_xy_yyyzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 68); 

                auto tg_xy_yyzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 69); 

                auto tg_xy_yzzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 70); 

                auto tg_xy_zzzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 71); 

                auto tg_xz_xxxxxxx_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 72); 

                auto tg_xz_xxxxxxy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 73); 

                auto tg_xz_xxxxxxz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 74); 

                auto tg_xz_xxxxxyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 75); 

                auto tg_xz_xxxxxyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 76); 

                auto tg_xz_xxxxxzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 77); 

                auto tg_xz_xxxxyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 78); 

                auto tg_xz_xxxxyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 79); 

                auto tg_xz_xxxxyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 80); 

                auto tg_xz_xxxxzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 81); 

                auto tg_xz_xxxyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 82); 

                auto tg_xz_xxxyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 83); 

                auto tg_xz_xxxyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 84); 

                auto tg_xz_xxxyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 85); 

                auto tg_xz_xxxzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 86); 

                auto tg_xz_xxyyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 87); 

                auto tg_xz_xxyyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 88); 

                auto tg_xz_xxyyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 89); 

                auto tg_xx_xxxxxxx_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx); 

                auto tg_xx_xxxxxxy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 1); 

                auto tg_xx_xxxxxxz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 2); 

                auto tg_xx_xxxxxyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 3); 

                auto tg_xx_xxxxxyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 4); 

                auto tg_xx_xxxxxzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 5); 

                auto tg_xx_xxxxyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 6); 

                auto tg_xx_xxxxyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 7); 

                auto tg_xx_xxxxyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 8); 

                auto tg_xx_xxxxzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 9); 

                auto tg_xx_xxxyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 10); 

                auto tg_xx_xxxyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 11); 

                auto tg_xx_xxxyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 12); 

                auto tg_xx_xxxyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 13); 

                auto tg_xx_xxxzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 14); 

                auto tg_xx_xxyyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 15); 

                auto tg_xx_xxyyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 16); 

                auto tg_xx_xxyyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 17); 

                auto tg_xx_xxyyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 18); 

                auto tg_xx_xxyzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 19); 

                auto tg_xx_xxzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 20); 

                auto tg_xx_xyyyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 21); 

                auto tg_xx_xyyyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 22); 

                auto tg_xx_xyyyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 23); 

                auto tg_xx_xyyyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 24); 

                auto tg_xx_xyyzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 25); 

                auto tg_xx_xyzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 26); 

                auto tg_xx_xzzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 27); 

                auto tg_xx_yyyyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 28); 

                auto tg_xx_yyyyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 29); 

                auto tg_xx_yyyyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 30); 

                auto tg_xx_yyyyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 31); 

                auto tg_xx_yyyzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 32); 

                auto tg_xx_yyzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 33); 

                auto tg_xx_yzzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 34); 

                auto tg_xx_zzzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 35); 

                auto tg_xy_xxxxxxx_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 36); 

                auto tg_xy_xxxxxxy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 37); 

                auto tg_xy_xxxxxxz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 38); 

                auto tg_xy_xxxxxyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 39); 

                auto tg_xy_xxxxxyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 40); 

                auto tg_xy_xxxxxzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 41); 

                auto tg_xy_xxxxyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 42); 

                auto tg_xy_xxxxyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 43); 

                auto tg_xy_xxxxyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 44); 

                auto tg_xy_xxxxzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 45); 

                auto tg_xy_xxxyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 46); 

                auto tg_xy_xxxyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 47); 

                auto tg_xy_xxxyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 48); 

                auto tg_xy_xxxyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 49); 

                auto tg_xy_xxxzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 50); 

                auto tg_xy_xxyyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 51); 

                auto tg_xy_xxyyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 52); 

                auto tg_xy_xxyyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 53); 

                auto tg_xy_xxyyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 54); 

                auto tg_xy_xxyzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 55); 

                auto tg_xy_xxzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 56); 

                auto tg_xy_xyyyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 57); 

                auto tg_xy_xyyyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 58); 

                auto tg_xy_xyyyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 59); 

                auto tg_xy_xyyyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 60); 

                auto tg_xy_xyyzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 61); 

                auto tg_xy_xyzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 62); 

                auto tg_xy_xzzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 63); 

                auto tg_xy_yyyyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 64); 

                auto tg_xy_yyyyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 65); 

                auto tg_xy_yyyyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 66); 

                auto tg_xy_yyyyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 67); 

                auto tg_xy_yyyzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 68); 

                auto tg_xy_yyzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 69); 

                auto tg_xy_yzzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 70); 

                auto tg_xy_zzzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 71); 

                auto tg_xz_xxxxxxx_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 72); 

                auto tg_xz_xxxxxxy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 73); 

                auto tg_xz_xxxxxxz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 74); 

                auto tg_xz_xxxxxyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 75); 

                auto tg_xz_xxxxxyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 76); 

                auto tg_xz_xxxxxzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 77); 

                auto tg_xz_xxxxyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 78); 

                auto tg_xz_xxxxyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 79); 

                auto tg_xz_xxxxyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 80); 

                auto tg_xz_xxxxzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 81); 

                auto tg_xz_xxxyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 82); 

                auto tg_xz_xxxyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 83); 

                auto tg_xz_xxxyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 84); 

                auto tg_xz_xxxyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 85); 

                auto tg_xz_xxxzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 86); 

                auto tg_xz_xxyyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 87); 

                auto tg_xz_xxyyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 88); 

                auto tg_xz_xxyyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 89); 

                auto tg_x_xxxxxxx_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx); 

                auto tg_x_xxxxxxy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 1); 

                auto tg_x_xxxxxxz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 2); 

                auto tg_x_xxxxxyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 3); 

                auto tg_x_xxxxxyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 4); 

                auto tg_x_xxxxxzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 5); 

                auto tg_x_xxxxyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 6); 

                auto tg_x_xxxxyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 7); 

                auto tg_x_xxxxyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 8); 

                auto tg_x_xxxxzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 9); 

                auto tg_x_xxxyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 10); 

                auto tg_x_xxxyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 11); 

                auto tg_x_xxxyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 12); 

                auto tg_x_xxxyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 13); 

                auto tg_x_xxxzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 14); 

                auto tg_x_xxyyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 15); 

                auto tg_x_xxyyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 16); 

                auto tg_x_xxyyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 17); 

                auto tg_x_xxyyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 18); 

                auto tg_x_xxyzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 19); 

                auto tg_x_xxzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 20); 

                auto tg_x_xyyyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 21); 

                auto tg_x_xyyyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 22); 

                auto tg_x_xyyyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 23); 

                auto tg_x_xyyyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 24); 

                auto tg_x_xyyzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 25); 

                auto tg_x_xyzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 26); 

                auto tg_x_xzzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 27); 

                auto tg_x_yyyyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 28); 

                auto tg_x_yyyyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 29); 

                auto tg_x_yyyyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 30); 

                auto tg_x_yyyyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 31); 

                auto tg_x_yyyzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 32); 

                auto tg_x_yyzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 33); 

                auto tg_x_yzzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 34); 

                auto tg_x_zzzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 35); 

                auto tg_y_xxxxxxx_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 36); 

                auto tg_y_xxxxxxy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 37); 

                auto tg_y_xxxxxxz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 38); 

                auto tg_y_xxxxxyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 39); 

                auto tg_y_xxxxxyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 40); 

                auto tg_y_xxxxxzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 41); 

                auto tg_y_xxxxyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 42); 

                auto tg_y_xxxxyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 43); 

                auto tg_y_xxxxyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 44); 

                auto tg_y_xxxxzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 45); 

                auto tg_y_xxxyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 46); 

                auto tg_y_xxxyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 47); 

                auto tg_y_xxxyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 48); 

                auto tg_y_xxxyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 49); 

                auto tg_y_xxxzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 50); 

                auto tg_y_xxyyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 51); 

                auto tg_y_xxyyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 52); 

                auto tg_y_xxyyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 53); 

                auto tg_y_xxyyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 54); 

                auto tg_y_xxyzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 55); 

                auto tg_y_xxzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 56); 

                auto tg_y_xyyyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 57); 

                auto tg_y_xyyyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 58); 

                auto tg_y_xyyyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 59); 

                auto tg_y_xyyyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 60); 

                auto tg_y_xyyzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 61); 

                auto tg_y_xyzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 62); 

                auto tg_y_xzzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 63); 

                auto tg_y_yyyyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 64); 

                auto tg_y_yyyyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 65); 

                auto tg_y_yyyyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 66); 

                auto tg_y_yyyyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 67); 

                auto tg_y_yyyzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 68); 

                auto tg_y_yyzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 69); 

                auto tg_y_yzzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 70); 

                auto tg_y_zzzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 71); 

                auto tg_z_xxxxxxx_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 72); 

                auto tg_z_xxxxxxy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 73); 

                auto tg_z_xxxxxxz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 74); 

                auto tg_z_xxxxxyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 75); 

                auto tg_z_xxxxxyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 76); 

                auto tg_z_xxxxxzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 77); 

                auto tg_z_xxxxyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 78); 

                auto tg_z_xxxxyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 79); 

                auto tg_z_xxxxyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 80); 

                auto tg_z_xxxxzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 81); 

                auto tg_z_xxxyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 82); 

                auto tg_z_xxxyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 83); 

                auto tg_z_xxxyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 84); 

                auto tg_z_xxxyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 85); 

                auto tg_z_xxxzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 86); 

                auto tg_z_xxyyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 87); 

                auto tg_z_xxyyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 88); 

                auto tg_z_xxyyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 89); 

                auto tg_x_xxxxxxx_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx); 

                auto tg_x_xxxxxxy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 1); 

                auto tg_x_xxxxxxz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 2); 

                auto tg_x_xxxxxyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 3); 

                auto tg_x_xxxxxyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 4); 

                auto tg_x_xxxxxzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 5); 

                auto tg_x_xxxxyyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 6); 

                auto tg_x_xxxxyyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 7); 

                auto tg_x_xxxxyzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 8); 

                auto tg_x_xxxxzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 9); 

                auto tg_x_xxxyyyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 10); 

                auto tg_x_xxxyyyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 11); 

                auto tg_x_xxxyyzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 12); 

                auto tg_x_xxxyzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 13); 

                auto tg_x_xxxzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 14); 

                auto tg_x_xxyyyyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 15); 

                auto tg_x_xxyyyyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 16); 

                auto tg_x_xxyyyzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 17); 

                auto tg_x_xxyyzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 18); 

                auto tg_x_xxyzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 19); 

                auto tg_x_xxzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 20); 

                auto tg_x_xyyyyyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 21); 

                auto tg_x_xyyyyyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 22); 

                auto tg_x_xyyyyzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 23); 

                auto tg_x_xyyyzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 24); 

                auto tg_x_xyyzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 25); 

                auto tg_x_xyzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 26); 

                auto tg_x_xzzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 27); 

                auto tg_x_yyyyyyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 28); 

                auto tg_x_yyyyyyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 29); 

                auto tg_x_yyyyyzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 30); 

                auto tg_x_yyyyzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 31); 

                auto tg_x_yyyzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 32); 

                auto tg_x_yyzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 33); 

                auto tg_x_yzzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 34); 

                auto tg_x_zzzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 35); 

                auto tg_y_xxxxxxx_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 36); 

                auto tg_y_xxxxxxy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 37); 

                auto tg_y_xxxxxxz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 38); 

                auto tg_y_xxxxxyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 39); 

                auto tg_y_xxxxxyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 40); 

                auto tg_y_xxxxxzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 41); 

                auto tg_y_xxxxyyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 42); 

                auto tg_y_xxxxyyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 43); 

                auto tg_y_xxxxyzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 44); 

                auto tg_y_xxxxzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 45); 

                auto tg_y_xxxyyyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 46); 

                auto tg_y_xxxyyyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 47); 

                auto tg_y_xxxyyzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 48); 

                auto tg_y_xxxyzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 49); 

                auto tg_y_xxxzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 50); 

                auto tg_y_xxyyyyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 51); 

                auto tg_y_xxyyyyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 52); 

                auto tg_y_xxyyyzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 53); 

                auto tg_y_xxyyzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 54); 

                auto tg_y_xxyzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 55); 

                auto tg_y_xxzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 56); 

                auto tg_y_xyyyyyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 57); 

                auto tg_y_xyyyyyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 58); 

                auto tg_y_xyyyyzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 59); 

                auto tg_y_xyyyzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 60); 

                auto tg_y_xyyzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 61); 

                auto tg_y_xyzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 62); 

                auto tg_y_xzzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 63); 

                auto tg_y_yyyyyyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 64); 

                auto tg_y_yyyyyyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 65); 

                auto tg_y_yyyyyzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 66); 

                auto tg_y_yyyyzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 67); 

                auto tg_y_yyyzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 68); 

                auto tg_y_yyzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 69); 

                auto tg_y_yzzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 70); 

                auto tg_y_zzzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 71); 

                auto tg_z_xxxxxxx_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 72); 

                auto tg_z_xxxxxxy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 73); 

                auto tg_z_xxxxxxz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 74); 

                auto tg_z_xxxxxyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 75); 

                auto tg_z_xxxxxyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 76); 

                auto tg_z_xxxxxzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 77); 

                auto tg_z_xxxxyyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 78); 

                auto tg_z_xxxxyyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 79); 

                auto tg_z_xxxxyzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 80); 

                auto tg_z_xxxxzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 81); 

                auto tg_z_xxxyyyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 82); 

                auto tg_z_xxxyyyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 83); 

                auto tg_z_xxxyyzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 84); 

                auto tg_z_xxxyzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 85); 

                auto tg_z_xxxzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 86); 

                auto tg_z_xxyyyyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 87); 

                auto tg_z_xxyyyyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 88); 

                auto tg_z_xxyyyzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 89); 

                auto tg_xx_xxxxxx_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx); 

                auto tg_xx_xxxxxy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 1); 

                auto tg_xx_xxxxxz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 2); 

                auto tg_xx_xxxxyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 3); 

                auto tg_xx_xxxxyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 4); 

                auto tg_xx_xxxxzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 5); 

                auto tg_xx_xxxyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 6); 

                auto tg_xx_xxxyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 7); 

                auto tg_xx_xxxyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 8); 

                auto tg_xx_xxxzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 9); 

                auto tg_xx_xxyyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 10); 

                auto tg_xx_xxyyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 11); 

                auto tg_xx_xxyyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 12); 

                auto tg_xx_xxyzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 13); 

                auto tg_xx_xxzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 14); 

                auto tg_xx_xyyyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 15); 

                auto tg_xx_xyyyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 16); 

                auto tg_xx_xyyyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 17); 

                auto tg_xx_xyyzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 18); 

                auto tg_xx_xyzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 19); 

                auto tg_xx_xzzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 20); 

                auto tg_xx_yyyyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 21); 

                auto tg_xx_yyyyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 22); 

                auto tg_xx_yyyyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 23); 

                auto tg_xx_yyyzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 24); 

                auto tg_xx_yyzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 25); 

                auto tg_xx_yzzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 26); 

                auto tg_xx_zzzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 27); 

                auto tg_xy_xxxxxx_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 28); 

                auto tg_xy_xxxxxy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 29); 

                auto tg_xy_xxxxxz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 30); 

                auto tg_xy_xxxxyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 31); 

                auto tg_xy_xxxxyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 32); 

                auto tg_xy_xxxxzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 33); 

                auto tg_xy_xxxyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 34); 

                auto tg_xy_xxxyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 35); 

                auto tg_xy_xxxyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 36); 

                auto tg_xy_xxxzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 37); 

                auto tg_xy_xxyyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 38); 

                auto tg_xy_xxyyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 39); 

                auto tg_xy_xxyyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 40); 

                auto tg_xy_xxyzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 41); 

                auto tg_xy_xxzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 42); 

                auto tg_xy_xyyyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 43); 

                auto tg_xy_xyyyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 44); 

                auto tg_xy_xyyyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 45); 

                auto tg_xy_xyyzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 46); 

                auto tg_xy_xyzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 47); 

                auto tg_xy_xzzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 48); 

                auto tg_xy_yyyyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 49); 

                auto tg_xy_yyyyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 50); 

                auto tg_xy_yyyyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 51); 

                auto tg_xy_yyyzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 52); 

                auto tg_xy_yyzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 53); 

                auto tg_xy_yzzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 54); 

                auto tg_xy_zzzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 55); 

                auto tg_xz_xxxxxx_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 56); 

                auto tg_xz_xxxxxy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 57); 

                auto tg_xz_xxxxxz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 58); 

                auto tg_xz_xxxxyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 59); 

                auto tg_xz_xxxxyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 60); 

                auto tg_xz_xxxxzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 61); 

                auto tg_xz_xxxyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 62); 

                auto tg_xz_xxxyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 63); 

                auto tg_xz_xxxyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 64); 

                auto tg_xz_xxxzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 65); 

                auto tg_xz_xxyyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 66); 

                auto tg_xz_xxyyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 67); 

                auto tg_xz_xxyyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 68); 

                auto tg_xz_xxyzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 69); 

                auto tg_xz_xxzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 70); 

                auto tg_xz_xyyyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 71); 

                auto tg_xz_xyyyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 72); 

                auto tg_xz_xyyyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 73); 

                // set up pointers to integrals

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

                // Batch of Integrals (0,90)

                #pragma omp simd aligned(fxn, fza, tg_x_xxxxxxx_0, tg_x_xxxxxxx_1, tg_x_xxxxxxy_0, \
                                         tg_x_xxxxxxy_1, tg_x_xxxxxxz_0, tg_x_xxxxxxz_1, tg_x_xxxxxyy_0, tg_x_xxxxxyy_1, \
                                         tg_x_xxxxxyz_0, tg_x_xxxxxyz_1, tg_x_xxxxxzz_0, tg_x_xxxxxzz_1, tg_x_xxxxyyy_0, \
                                         tg_x_xxxxyyy_1, tg_x_xxxxyyz_0, tg_x_xxxxyyz_1, tg_x_xxxxyzz_0, tg_x_xxxxyzz_1, \
                                         tg_x_xxxxzzz_0, tg_x_xxxxzzz_1, tg_x_xxxyyyy_0, tg_x_xxxyyyy_1, tg_x_xxxyyyz_0, \
                                         tg_x_xxxyyyz_1, tg_x_xxxyyzz_0, tg_x_xxxyyzz_1, tg_x_xxxyzzz_0, tg_x_xxxyzzz_1, \
                                         tg_x_xxxzzzz_0, tg_x_xxxzzzz_1, tg_x_xxyyyyy_0, tg_x_xxyyyyy_1, tg_x_xxyyyyz_0, \
                                         tg_x_xxyyyyz_1, tg_x_xxyyyzz_0, tg_x_xxyyyzz_1, tg_x_xxyyzzz_0, tg_x_xxyyzzz_1, \
                                         tg_x_xxyzzzz_0, tg_x_xxyzzzz_1, tg_x_xxzzzzz_0, tg_x_xxzzzzz_1, tg_x_xyyyyyy_0, \
                                         tg_x_xyyyyyy_1, tg_x_xyyyyyz_0, tg_x_xyyyyyz_1, tg_x_xyyyyzz_0, tg_x_xyyyyzz_1, \
                                         tg_x_xyyyzzz_0, tg_x_xyyyzzz_1, tg_x_xyyzzzz_0, tg_x_xyyzzzz_1, tg_x_xyzzzzz_0, \
                                         tg_x_xyzzzzz_1, tg_x_xzzzzzz_0, tg_x_xzzzzzz_1, tg_x_yyyyyyy_0, tg_x_yyyyyyy_1, \
                                         tg_x_yyyyyyz_0, tg_x_yyyyyyz_1, tg_x_yyyyyzz_0, tg_x_yyyyyzz_1, tg_x_yyyyzzz_0, \
                                         tg_x_yyyyzzz_1, tg_x_yyyzzzz_0, tg_x_yyyzzzz_1, tg_x_yyzzzzz_0, tg_x_yyzzzzz_1, \
                                         tg_x_yzzzzzz_0, tg_x_yzzzzzz_1, tg_x_zzzzzzz_0, tg_x_zzzzzzz_1, tg_xx_xxxxxx_1, \
                                         tg_xx_xxxxxxx_0, tg_xx_xxxxxxx_1, tg_xx_xxxxxxy_0, tg_xx_xxxxxxy_1, tg_xx_xxxxxxz_0, \
                                         tg_xx_xxxxxxz_1, tg_xx_xxxxxy_1, tg_xx_xxxxxyy_0, tg_xx_xxxxxyy_1, tg_xx_xxxxxyz_0, \
                                         tg_xx_xxxxxyz_1, tg_xx_xxxxxz_1, tg_xx_xxxxxzz_0, tg_xx_xxxxxzz_1, tg_xx_xxxxyy_1, \
                                         tg_xx_xxxxyyy_0, tg_xx_xxxxyyy_1, tg_xx_xxxxyyz_0, tg_xx_xxxxyyz_1, tg_xx_xxxxyz_1, \
                                         tg_xx_xxxxyzz_0, tg_xx_xxxxyzz_1, tg_xx_xxxxzz_1, tg_xx_xxxxzzz_0, tg_xx_xxxxzzz_1, \
                                         tg_xx_xxxyyy_1, tg_xx_xxxyyyy_0, tg_xx_xxxyyyy_1, tg_xx_xxxyyyz_0, tg_xx_xxxyyyz_1, \
                                         tg_xx_xxxyyz_1, tg_xx_xxxyyzz_0, tg_xx_xxxyyzz_1, tg_xx_xxxyzz_1, tg_xx_xxxyzzz_0, \
                                         tg_xx_xxxyzzz_1, tg_xx_xxxzzz_1, tg_xx_xxxzzzz_0, tg_xx_xxxzzzz_1, tg_xx_xxyyyy_1, \
                                         tg_xx_xxyyyyy_0, tg_xx_xxyyyyy_1, tg_xx_xxyyyyz_0, tg_xx_xxyyyyz_1, tg_xx_xxyyyz_1, \
                                         tg_xx_xxyyyzz_0, tg_xx_xxyyyzz_1, tg_xx_xxyyzz_1, tg_xx_xxyyzzz_0, tg_xx_xxyyzzz_1, \
                                         tg_xx_xxyzzz_1, tg_xx_xxyzzzz_0, tg_xx_xxyzzzz_1, tg_xx_xxzzzz_1, tg_xx_xxzzzzz_0, \
                                         tg_xx_xxzzzzz_1, tg_xx_xyyyyy_1, tg_xx_xyyyyyy_0, tg_xx_xyyyyyy_1, tg_xx_xyyyyyz_0, \
                                         tg_xx_xyyyyyz_1, tg_xx_xyyyyz_1, tg_xx_xyyyyzz_0, tg_xx_xyyyyzz_1, tg_xx_xyyyzz_1, \
                                         tg_xx_xyyyzzz_0, tg_xx_xyyyzzz_1, tg_xx_xyyzzz_1, tg_xx_xyyzzzz_0, tg_xx_xyyzzzz_1, \
                                         tg_xx_xyzzzz_1, tg_xx_xyzzzzz_0, tg_xx_xyzzzzz_1, tg_xx_xzzzzz_1, tg_xx_xzzzzzz_0, \
                                         tg_xx_xzzzzzz_1, tg_xx_yyyyyy_1, tg_xx_yyyyyyy_0, tg_xx_yyyyyyy_1, tg_xx_yyyyyyz_0, \
                                         tg_xx_yyyyyyz_1, tg_xx_yyyyyz_1, tg_xx_yyyyyzz_0, tg_xx_yyyyyzz_1, tg_xx_yyyyzz_1, \
                                         tg_xx_yyyyzzz_0, tg_xx_yyyyzzz_1, tg_xx_yyyzzz_1, tg_xx_yyyzzzz_0, tg_xx_yyyzzzz_1, \
                                         tg_xx_yyzzzz_1, tg_xx_yyzzzzz_0, tg_xx_yyzzzzz_1, tg_xx_yzzzzz_1, tg_xx_yzzzzzz_0, \
                                         tg_xx_yzzzzzz_1, tg_xx_zzzzzz_1, tg_xx_zzzzzzz_0, tg_xx_zzzzzzz_1, tg_xxx_xxxxxxx_0, \
                                         tg_xxx_xxxxxxy_0, tg_xxx_xxxxxxz_0, tg_xxx_xxxxxyy_0, tg_xxx_xxxxxyz_0, \
                                         tg_xxx_xxxxxzz_0, tg_xxx_xxxxyyy_0, tg_xxx_xxxxyyz_0, tg_xxx_xxxxyzz_0, \
                                         tg_xxx_xxxxzzz_0, tg_xxx_xxxyyyy_0, tg_xxx_xxxyyyz_0, tg_xxx_xxxyyzz_0, \
                                         tg_xxx_xxxyzzz_0, tg_xxx_xxxzzzz_0, tg_xxx_xxyyyyy_0, tg_xxx_xxyyyyz_0, \
                                         tg_xxx_xxyyyzz_0, tg_xxx_xxyyzzz_0, tg_xxx_xxyzzzz_0, tg_xxx_xxzzzzz_0, \
                                         tg_xxx_xyyyyyy_0, tg_xxx_xyyyyyz_0, tg_xxx_xyyyyzz_0, tg_xxx_xyyyzzz_0, \
                                         tg_xxx_xyyzzzz_0, tg_xxx_xyzzzzz_0, tg_xxx_xzzzzzz_0, tg_xxx_yyyyyyy_0, \
                                         tg_xxx_yyyyyyz_0, tg_xxx_yyyyyzz_0, tg_xxx_yyyyzzz_0, tg_xxx_yyyzzzz_0, \
                                         tg_xxx_yyzzzzz_0, tg_xxx_yzzzzzz_0, tg_xxx_zzzzzzz_0, tg_xxy_xxxxxxx_0, \
                                         tg_xxy_xxxxxxy_0, tg_xxy_xxxxxxz_0, tg_xxy_xxxxxyy_0, tg_xxy_xxxxxyz_0, \
                                         tg_xxy_xxxxxzz_0, tg_xxy_xxxxyyy_0, tg_xxy_xxxxyyz_0, tg_xxy_xxxxyzz_0, \
                                         tg_xxy_xxxxzzz_0, tg_xxy_xxxyyyy_0, tg_xxy_xxxyyyz_0, tg_xxy_xxxyyzz_0, \
                                         tg_xxy_xxxyzzz_0, tg_xxy_xxxzzzz_0, tg_xxy_xxyyyyy_0, tg_xxy_xxyyyyz_0, \
                                         tg_xxy_xxyyyzz_0, tg_xxy_xxyyzzz_0, tg_xxy_xxyzzzz_0, tg_xxy_xxzzzzz_0, \
                                         tg_xxy_xyyyyyy_0, tg_xxy_xyyyyyz_0, tg_xxy_xyyyyzz_0, tg_xxy_xyyyzzz_0, \
                                         tg_xxy_xyyzzzz_0, tg_xxy_xyzzzzz_0, tg_xxy_xzzzzzz_0, tg_xxy_yyyyyyy_0, \
                                         tg_xxy_yyyyyyz_0, tg_xxy_yyyyyzz_0, tg_xxy_yyyyzzz_0, tg_xxy_yyyzzzz_0, \
                                         tg_xxy_yyzzzzz_0, tg_xxy_yzzzzzz_0, tg_xxy_zzzzzzz_0, tg_xxz_xxxxxxx_0, \
                                         tg_xxz_xxxxxxy_0, tg_xxz_xxxxxxz_0, tg_xxz_xxxxxyy_0, tg_xxz_xxxxxyz_0, \
                                         tg_xxz_xxxxxzz_0, tg_xxz_xxxxyyy_0, tg_xxz_xxxxyyz_0, tg_xxz_xxxxyzz_0, \
                                         tg_xxz_xxxxzzz_0, tg_xxz_xxxyyyy_0, tg_xxz_xxxyyyz_0, tg_xxz_xxxyyzz_0, \
                                         tg_xxz_xxxyzzz_0, tg_xxz_xxxzzzz_0, tg_xxz_xxyyyyy_0, tg_xxz_xxyyyyz_0, \
                                         tg_xxz_xxyyyzz_0, tg_xy_xxxxxx_1, tg_xy_xxxxxxx_0, tg_xy_xxxxxxx_1, tg_xy_xxxxxxy_0, \
                                         tg_xy_xxxxxxy_1, tg_xy_xxxxxxz_0, tg_xy_xxxxxxz_1, tg_xy_xxxxxy_1, tg_xy_xxxxxyy_0, \
                                         tg_xy_xxxxxyy_1, tg_xy_xxxxxyz_0, tg_xy_xxxxxyz_1, tg_xy_xxxxxz_1, tg_xy_xxxxxzz_0, \
                                         tg_xy_xxxxxzz_1, tg_xy_xxxxyy_1, tg_xy_xxxxyyy_0, tg_xy_xxxxyyy_1, tg_xy_xxxxyyz_0, \
                                         tg_xy_xxxxyyz_1, tg_xy_xxxxyz_1, tg_xy_xxxxyzz_0, tg_xy_xxxxyzz_1, tg_xy_xxxxzz_1, \
                                         tg_xy_xxxxzzz_0, tg_xy_xxxxzzz_1, tg_xy_xxxyyy_1, tg_xy_xxxyyyy_0, tg_xy_xxxyyyy_1, \
                                         tg_xy_xxxyyyz_0, tg_xy_xxxyyyz_1, tg_xy_xxxyyz_1, tg_xy_xxxyyzz_0, tg_xy_xxxyyzz_1, \
                                         tg_xy_xxxyzz_1, tg_xy_xxxyzzz_0, tg_xy_xxxyzzz_1, tg_xy_xxxzzz_1, tg_xy_xxxzzzz_0, \
                                         tg_xy_xxxzzzz_1, tg_xy_xxyyyy_1, tg_xy_xxyyyyy_0, tg_xy_xxyyyyy_1, tg_xy_xxyyyyz_0, \
                                         tg_xy_xxyyyyz_1, tg_xy_xxyyyz_1, tg_xy_xxyyyzz_0, tg_xy_xxyyyzz_1, tg_xy_xxyyzz_1, \
                                         tg_xy_xxyyzzz_0, tg_xy_xxyyzzz_1, tg_xy_xxyzzz_1, tg_xy_xxyzzzz_0, tg_xy_xxyzzzz_1, \
                                         tg_xy_xxzzzz_1, tg_xy_xxzzzzz_0, tg_xy_xxzzzzz_1, tg_xy_xyyyyy_1, tg_xy_xyyyyyy_0, \
                                         tg_xy_xyyyyyy_1, tg_xy_xyyyyyz_0, tg_xy_xyyyyyz_1, tg_xy_xyyyyz_1, tg_xy_xyyyyzz_0, \
                                         tg_xy_xyyyyzz_1, tg_xy_xyyyzz_1, tg_xy_xyyyzzz_0, tg_xy_xyyyzzz_1, tg_xy_xyyzzz_1, \
                                         tg_xy_xyyzzzz_0, tg_xy_xyyzzzz_1, tg_xy_xyzzzz_1, tg_xy_xyzzzzz_0, tg_xy_xyzzzzz_1, \
                                         tg_xy_xzzzzz_1, tg_xy_xzzzzzz_0, tg_xy_xzzzzzz_1, tg_xy_yyyyyy_1, tg_xy_yyyyyyy_0, \
                                         tg_xy_yyyyyyy_1, tg_xy_yyyyyyz_0, tg_xy_yyyyyyz_1, tg_xy_yyyyyz_1, tg_xy_yyyyyzz_0, \
                                         tg_xy_yyyyyzz_1, tg_xy_yyyyzz_1, tg_xy_yyyyzzz_0, tg_xy_yyyyzzz_1, tg_xy_yyyzzz_1, \
                                         tg_xy_yyyzzzz_0, tg_xy_yyyzzzz_1, tg_xy_yyzzzz_1, tg_xy_yyzzzzz_0, tg_xy_yyzzzzz_1, \
                                         tg_xy_yzzzzz_1, tg_xy_yzzzzzz_0, tg_xy_yzzzzzz_1, tg_xy_zzzzzz_1, tg_xy_zzzzzzz_0, \
                                         tg_xy_zzzzzzz_1, tg_xz_xxxxxx_1, tg_xz_xxxxxxx_0, tg_xz_xxxxxxx_1, tg_xz_xxxxxxy_0, \
                                         tg_xz_xxxxxxy_1, tg_xz_xxxxxxz_0, tg_xz_xxxxxxz_1, tg_xz_xxxxxy_1, tg_xz_xxxxxyy_0, \
                                         tg_xz_xxxxxyy_1, tg_xz_xxxxxyz_0, tg_xz_xxxxxyz_1, tg_xz_xxxxxz_1, tg_xz_xxxxxzz_0, \
                                         tg_xz_xxxxxzz_1, tg_xz_xxxxyy_1, tg_xz_xxxxyyy_0, tg_xz_xxxxyyy_1, tg_xz_xxxxyyz_0, \
                                         tg_xz_xxxxyyz_1, tg_xz_xxxxyz_1, tg_xz_xxxxyzz_0, tg_xz_xxxxyzz_1, tg_xz_xxxxzz_1, \
                                         tg_xz_xxxxzzz_0, tg_xz_xxxxzzz_1, tg_xz_xxxyyy_1, tg_xz_xxxyyyy_0, tg_xz_xxxyyyy_1, \
                                         tg_xz_xxxyyyz_0, tg_xz_xxxyyyz_1, tg_xz_xxxyyz_1, tg_xz_xxxyyzz_0, tg_xz_xxxyyzz_1, \
                                         tg_xz_xxxyzz_1, tg_xz_xxxyzzz_0, tg_xz_xxxyzzz_1, tg_xz_xxxzzz_1, tg_xz_xxxzzzz_0, \
                                         tg_xz_xxxzzzz_1, tg_xz_xxyyyy_1, tg_xz_xxyyyyy_0, tg_xz_xxyyyyy_1, tg_xz_xxyyyyz_0, \
                                         tg_xz_xxyyyyz_1, tg_xz_xxyyyz_1, tg_xz_xxyyyzz_0, tg_xz_xxyyyzz_1, tg_xz_xxyyzz_1, \
                                         tg_xz_xxyzzz_1, tg_xz_xxzzzz_1, tg_xz_xyyyyy_1, tg_xz_xyyyyz_1, tg_xz_xyyyzz_1, \
                                         tg_y_xxxxxxx_0, tg_y_xxxxxxx_1, tg_y_xxxxxxy_0, tg_y_xxxxxxy_1, tg_y_xxxxxxz_0, \
                                         tg_y_xxxxxxz_1, tg_y_xxxxxyy_0, tg_y_xxxxxyy_1, tg_y_xxxxxyz_0, tg_y_xxxxxyz_1, \
                                         tg_y_xxxxxzz_0, tg_y_xxxxxzz_1, tg_y_xxxxyyy_0, tg_y_xxxxyyy_1, tg_y_xxxxyyz_0, \
                                         tg_y_xxxxyyz_1, tg_y_xxxxyzz_0, tg_y_xxxxyzz_1, tg_y_xxxxzzz_0, tg_y_xxxxzzz_1, \
                                         tg_y_xxxyyyy_0, tg_y_xxxyyyy_1, tg_y_xxxyyyz_0, tg_y_xxxyyyz_1, tg_y_xxxyyzz_0, \
                                         tg_y_xxxyyzz_1, tg_y_xxxyzzz_0, tg_y_xxxyzzz_1, tg_y_xxxzzzz_0, tg_y_xxxzzzz_1, \
                                         tg_y_xxyyyyy_0, tg_y_xxyyyyy_1, tg_y_xxyyyyz_0, tg_y_xxyyyyz_1, tg_y_xxyyyzz_0, \
                                         tg_y_xxyyyzz_1, tg_y_xxyyzzz_0, tg_y_xxyyzzz_1, tg_y_xxyzzzz_0, tg_y_xxyzzzz_1, \
                                         tg_y_xxzzzzz_0, tg_y_xxzzzzz_1, tg_y_xyyyyyy_0, tg_y_xyyyyyy_1, tg_y_xyyyyyz_0, \
                                         tg_y_xyyyyyz_1, tg_y_xyyyyzz_0, tg_y_xyyyyzz_1, tg_y_xyyyzzz_0, tg_y_xyyyzzz_1, \
                                         tg_y_xyyzzzz_0, tg_y_xyyzzzz_1, tg_y_xyzzzzz_0, tg_y_xyzzzzz_1, tg_y_xzzzzzz_0, \
                                         tg_y_xzzzzzz_1, tg_y_yyyyyyy_0, tg_y_yyyyyyy_1, tg_y_yyyyyyz_0, tg_y_yyyyyyz_1, \
                                         tg_y_yyyyyzz_0, tg_y_yyyyyzz_1, tg_y_yyyyzzz_0, tg_y_yyyyzzz_1, tg_y_yyyzzzz_0, \
                                         tg_y_yyyzzzz_1, tg_y_yyzzzzz_0, tg_y_yyzzzzz_1, tg_y_yzzzzzz_0, tg_y_yzzzzzz_1, \
                                         tg_y_zzzzzzz_0, tg_y_zzzzzzz_1, tg_z_xxxxxxx_0, tg_z_xxxxxxx_1, tg_z_xxxxxxy_0, \
                                         tg_z_xxxxxxy_1, tg_z_xxxxxxz_0, tg_z_xxxxxxz_1, tg_z_xxxxxyy_0, tg_z_xxxxxyy_1, \
                                         tg_z_xxxxxyz_0, tg_z_xxxxxyz_1, tg_z_xxxxxzz_0, tg_z_xxxxxzz_1, tg_z_xxxxyyy_0, \
                                         tg_z_xxxxyyy_1, tg_z_xxxxyyz_0, tg_z_xxxxyyz_1, tg_z_xxxxyzz_0, tg_z_xxxxyzz_1, \
                                         tg_z_xxxxzzz_0, tg_z_xxxxzzz_1, tg_z_xxxyyyy_0, tg_z_xxxyyyy_1, tg_z_xxxyyyz_0, \
                                         tg_z_xxxyyyz_1, tg_z_xxxyyzz_0, tg_z_xxxyyzz_1, tg_z_xxxyzzz_0, tg_z_xxxyzzz_1, \
                                         tg_z_xxxzzzz_0, tg_z_xxxzzzz_1, tg_z_xxyyyyy_0, tg_z_xxyyyyy_1, tg_z_xxyyyyz_0, \
                                         tg_z_xxyyyyz_1, tg_z_xxyyyzz_0, tg_z_xxyyyzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxx_xxxxxxx_0[j] = pb_x * tg_xx_xxxxxxx_0[j] + wp_x[j] * tg_xx_xxxxxxx_1[j] + fl1_fx * tg_x_xxxxxxx_0[j] - fl1_fx * fl1_fza * tg_x_xxxxxxx_1[j] + 3.5 * fl1_fxn * tg_xx_xxxxxx_1[j];

                    tg_xxx_xxxxxxy_0[j] = pb_x * tg_xx_xxxxxxy_0[j] + wp_x[j] * tg_xx_xxxxxxy_1[j] + fl1_fx * tg_x_xxxxxxy_0[j] - fl1_fx * fl1_fza * tg_x_xxxxxxy_1[j] + 3.0 * fl1_fxn * tg_xx_xxxxxy_1[j];

                    tg_xxx_xxxxxxz_0[j] = pb_x * tg_xx_xxxxxxz_0[j] + wp_x[j] * tg_xx_xxxxxxz_1[j] + fl1_fx * tg_x_xxxxxxz_0[j] - fl1_fx * fl1_fza * tg_x_xxxxxxz_1[j] + 3.0 * fl1_fxn * tg_xx_xxxxxz_1[j];

                    tg_xxx_xxxxxyy_0[j] = pb_x * tg_xx_xxxxxyy_0[j] + wp_x[j] * tg_xx_xxxxxyy_1[j] + fl1_fx * tg_x_xxxxxyy_0[j] - fl1_fx * fl1_fza * tg_x_xxxxxyy_1[j] + 2.5 * fl1_fxn * tg_xx_xxxxyy_1[j];

                    tg_xxx_xxxxxyz_0[j] = pb_x * tg_xx_xxxxxyz_0[j] + wp_x[j] * tg_xx_xxxxxyz_1[j] + fl1_fx * tg_x_xxxxxyz_0[j] - fl1_fx * fl1_fza * tg_x_xxxxxyz_1[j] + 2.5 * fl1_fxn * tg_xx_xxxxyz_1[j];

                    tg_xxx_xxxxxzz_0[j] = pb_x * tg_xx_xxxxxzz_0[j] + wp_x[j] * tg_xx_xxxxxzz_1[j] + fl1_fx * tg_x_xxxxxzz_0[j] - fl1_fx * fl1_fza * tg_x_xxxxxzz_1[j] + 2.5 * fl1_fxn * tg_xx_xxxxzz_1[j];

                    tg_xxx_xxxxyyy_0[j] = pb_x * tg_xx_xxxxyyy_0[j] + wp_x[j] * tg_xx_xxxxyyy_1[j] + fl1_fx * tg_x_xxxxyyy_0[j] - fl1_fx * fl1_fza * tg_x_xxxxyyy_1[j] + 2.0 * fl1_fxn * tg_xx_xxxyyy_1[j];

                    tg_xxx_xxxxyyz_0[j] = pb_x * tg_xx_xxxxyyz_0[j] + wp_x[j] * tg_xx_xxxxyyz_1[j] + fl1_fx * tg_x_xxxxyyz_0[j] - fl1_fx * fl1_fza * tg_x_xxxxyyz_1[j] + 2.0 * fl1_fxn * tg_xx_xxxyyz_1[j];

                    tg_xxx_xxxxyzz_0[j] = pb_x * tg_xx_xxxxyzz_0[j] + wp_x[j] * tg_xx_xxxxyzz_1[j] + fl1_fx * tg_x_xxxxyzz_0[j] - fl1_fx * fl1_fza * tg_x_xxxxyzz_1[j] + 2.0 * fl1_fxn * tg_xx_xxxyzz_1[j];

                    tg_xxx_xxxxzzz_0[j] = pb_x * tg_xx_xxxxzzz_0[j] + wp_x[j] * tg_xx_xxxxzzz_1[j] + fl1_fx * tg_x_xxxxzzz_0[j] - fl1_fx * fl1_fza * tg_x_xxxxzzz_1[j] + 2.0 * fl1_fxn * tg_xx_xxxzzz_1[j];

                    tg_xxx_xxxyyyy_0[j] = pb_x * tg_xx_xxxyyyy_0[j] + wp_x[j] * tg_xx_xxxyyyy_1[j] + fl1_fx * tg_x_xxxyyyy_0[j] - fl1_fx * fl1_fza * tg_x_xxxyyyy_1[j] + 1.5 * fl1_fxn * tg_xx_xxyyyy_1[j];

                    tg_xxx_xxxyyyz_0[j] = pb_x * tg_xx_xxxyyyz_0[j] + wp_x[j] * tg_xx_xxxyyyz_1[j] + fl1_fx * tg_x_xxxyyyz_0[j] - fl1_fx * fl1_fza * tg_x_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_xx_xxyyyz_1[j];

                    tg_xxx_xxxyyzz_0[j] = pb_x * tg_xx_xxxyyzz_0[j] + wp_x[j] * tg_xx_xxxyyzz_1[j] + fl1_fx * tg_x_xxxyyzz_0[j] - fl1_fx * fl1_fza * tg_x_xxxyyzz_1[j] + 1.5 * fl1_fxn * tg_xx_xxyyzz_1[j];

                    tg_xxx_xxxyzzz_0[j] = pb_x * tg_xx_xxxyzzz_0[j] + wp_x[j] * tg_xx_xxxyzzz_1[j] + fl1_fx * tg_x_xxxyzzz_0[j] - fl1_fx * fl1_fza * tg_x_xxxyzzz_1[j] + 1.5 * fl1_fxn * tg_xx_xxyzzz_1[j];

                    tg_xxx_xxxzzzz_0[j] = pb_x * tg_xx_xxxzzzz_0[j] + wp_x[j] * tg_xx_xxxzzzz_1[j] + fl1_fx * tg_x_xxxzzzz_0[j] - fl1_fx * fl1_fza * tg_x_xxxzzzz_1[j] + 1.5 * fl1_fxn * tg_xx_xxzzzz_1[j];

                    tg_xxx_xxyyyyy_0[j] = pb_x * tg_xx_xxyyyyy_0[j] + wp_x[j] * tg_xx_xxyyyyy_1[j] + fl1_fx * tg_x_xxyyyyy_0[j] - fl1_fx * fl1_fza * tg_x_xxyyyyy_1[j] + fl1_fxn * tg_xx_xyyyyy_1[j];

                    tg_xxx_xxyyyyz_0[j] = pb_x * tg_xx_xxyyyyz_0[j] + wp_x[j] * tg_xx_xxyyyyz_1[j] + fl1_fx * tg_x_xxyyyyz_0[j] - fl1_fx * fl1_fza * tg_x_xxyyyyz_1[j] + fl1_fxn * tg_xx_xyyyyz_1[j];

                    tg_xxx_xxyyyzz_0[j] = pb_x * tg_xx_xxyyyzz_0[j] + wp_x[j] * tg_xx_xxyyyzz_1[j] + fl1_fx * tg_x_xxyyyzz_0[j] - fl1_fx * fl1_fza * tg_x_xxyyyzz_1[j] + fl1_fxn * tg_xx_xyyyzz_1[j];

                    tg_xxx_xxyyzzz_0[j] = pb_x * tg_xx_xxyyzzz_0[j] + wp_x[j] * tg_xx_xxyyzzz_1[j] + fl1_fx * tg_x_xxyyzzz_0[j] - fl1_fx * fl1_fza * tg_x_xxyyzzz_1[j] + fl1_fxn * tg_xx_xyyzzz_1[j];

                    tg_xxx_xxyzzzz_0[j] = pb_x * tg_xx_xxyzzzz_0[j] + wp_x[j] * tg_xx_xxyzzzz_1[j] + fl1_fx * tg_x_xxyzzzz_0[j] - fl1_fx * fl1_fza * tg_x_xxyzzzz_1[j] + fl1_fxn * tg_xx_xyzzzz_1[j];

                    tg_xxx_xxzzzzz_0[j] = pb_x * tg_xx_xxzzzzz_0[j] + wp_x[j] * tg_xx_xxzzzzz_1[j] + fl1_fx * tg_x_xxzzzzz_0[j] - fl1_fx * fl1_fza * tg_x_xxzzzzz_1[j] + fl1_fxn * tg_xx_xzzzzz_1[j];

                    tg_xxx_xyyyyyy_0[j] = pb_x * tg_xx_xyyyyyy_0[j] + wp_x[j] * tg_xx_xyyyyyy_1[j] + fl1_fx * tg_x_xyyyyyy_0[j] - fl1_fx * fl1_fza * tg_x_xyyyyyy_1[j] + 0.5 * fl1_fxn * tg_xx_yyyyyy_1[j];

                    tg_xxx_xyyyyyz_0[j] = pb_x * tg_xx_xyyyyyz_0[j] + wp_x[j] * tg_xx_xyyyyyz_1[j] + fl1_fx * tg_x_xyyyyyz_0[j] - fl1_fx * fl1_fza * tg_x_xyyyyyz_1[j] + 0.5 * fl1_fxn * tg_xx_yyyyyz_1[j];

                    tg_xxx_xyyyyzz_0[j] = pb_x * tg_xx_xyyyyzz_0[j] + wp_x[j] * tg_xx_xyyyyzz_1[j] + fl1_fx * tg_x_xyyyyzz_0[j] - fl1_fx * fl1_fza * tg_x_xyyyyzz_1[j] + 0.5 * fl1_fxn * tg_xx_yyyyzz_1[j];

                    tg_xxx_xyyyzzz_0[j] = pb_x * tg_xx_xyyyzzz_0[j] + wp_x[j] * tg_xx_xyyyzzz_1[j] + fl1_fx * tg_x_xyyyzzz_0[j] - fl1_fx * fl1_fza * tg_x_xyyyzzz_1[j] + 0.5 * fl1_fxn * tg_xx_yyyzzz_1[j];

                    tg_xxx_xyyzzzz_0[j] = pb_x * tg_xx_xyyzzzz_0[j] + wp_x[j] * tg_xx_xyyzzzz_1[j] + fl1_fx * tg_x_xyyzzzz_0[j] - fl1_fx * fl1_fza * tg_x_xyyzzzz_1[j] + 0.5 * fl1_fxn * tg_xx_yyzzzz_1[j];

                    tg_xxx_xyzzzzz_0[j] = pb_x * tg_xx_xyzzzzz_0[j] + wp_x[j] * tg_xx_xyzzzzz_1[j] + fl1_fx * tg_x_xyzzzzz_0[j] - fl1_fx * fl1_fza * tg_x_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_xx_yzzzzz_1[j];

                    tg_xxx_xzzzzzz_0[j] = pb_x * tg_xx_xzzzzzz_0[j] + wp_x[j] * tg_xx_xzzzzzz_1[j] + fl1_fx * tg_x_xzzzzzz_0[j] - fl1_fx * fl1_fza * tg_x_xzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xx_zzzzzz_1[j];

                    tg_xxx_yyyyyyy_0[j] = pb_x * tg_xx_yyyyyyy_0[j] + wp_x[j] * tg_xx_yyyyyyy_1[j] + fl1_fx * tg_x_yyyyyyy_0[j] - fl1_fx * fl1_fza * tg_x_yyyyyyy_1[j];

                    tg_xxx_yyyyyyz_0[j] = pb_x * tg_xx_yyyyyyz_0[j] + wp_x[j] * tg_xx_yyyyyyz_1[j] + fl1_fx * tg_x_yyyyyyz_0[j] - fl1_fx * fl1_fza * tg_x_yyyyyyz_1[j];

                    tg_xxx_yyyyyzz_0[j] = pb_x * tg_xx_yyyyyzz_0[j] + wp_x[j] * tg_xx_yyyyyzz_1[j] + fl1_fx * tg_x_yyyyyzz_0[j] - fl1_fx * fl1_fza * tg_x_yyyyyzz_1[j];

                    tg_xxx_yyyyzzz_0[j] = pb_x * tg_xx_yyyyzzz_0[j] + wp_x[j] * tg_xx_yyyyzzz_1[j] + fl1_fx * tg_x_yyyyzzz_0[j] - fl1_fx * fl1_fza * tg_x_yyyyzzz_1[j];

                    tg_xxx_yyyzzzz_0[j] = pb_x * tg_xx_yyyzzzz_0[j] + wp_x[j] * tg_xx_yyyzzzz_1[j] + fl1_fx * tg_x_yyyzzzz_0[j] - fl1_fx * fl1_fza * tg_x_yyyzzzz_1[j];

                    tg_xxx_yyzzzzz_0[j] = pb_x * tg_xx_yyzzzzz_0[j] + wp_x[j] * tg_xx_yyzzzzz_1[j] + fl1_fx * tg_x_yyzzzzz_0[j] - fl1_fx * fl1_fza * tg_x_yyzzzzz_1[j];

                    tg_xxx_yzzzzzz_0[j] = pb_x * tg_xx_yzzzzzz_0[j] + wp_x[j] * tg_xx_yzzzzzz_1[j] + fl1_fx * tg_x_yzzzzzz_0[j] - fl1_fx * fl1_fza * tg_x_yzzzzzz_1[j];

                    tg_xxx_zzzzzzz_0[j] = pb_x * tg_xx_zzzzzzz_0[j] + wp_x[j] * tg_xx_zzzzzzz_1[j] + fl1_fx * tg_x_zzzzzzz_0[j] - fl1_fx * fl1_fza * tg_x_zzzzzzz_1[j];

                    tg_xxy_xxxxxxx_0[j] = pb_x * tg_xy_xxxxxxx_0[j] + wp_x[j] * tg_xy_xxxxxxx_1[j] + 0.5 * fl1_fx * tg_y_xxxxxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxxxxx_1[j] + 3.5 * fl1_fxn * tg_xy_xxxxxx_1[j];

                    tg_xxy_xxxxxxy_0[j] = pb_x * tg_xy_xxxxxxy_0[j] + wp_x[j] * tg_xy_xxxxxxy_1[j] + 0.5 * fl1_fx * tg_y_xxxxxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxxxxy_1[j] + 3.0 * fl1_fxn * tg_xy_xxxxxy_1[j];

                    tg_xxy_xxxxxxz_0[j] = pb_x * tg_xy_xxxxxxz_0[j] + wp_x[j] * tg_xy_xxxxxxz_1[j] + 0.5 * fl1_fx * tg_y_xxxxxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxxxxz_1[j] + 3.0 * fl1_fxn * tg_xy_xxxxxz_1[j];

                    tg_xxy_xxxxxyy_0[j] = pb_x * tg_xy_xxxxxyy_0[j] + wp_x[j] * tg_xy_xxxxxyy_1[j] + 0.5 * fl1_fx * tg_y_xxxxxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxxxyy_1[j] + 2.5 * fl1_fxn * tg_xy_xxxxyy_1[j];

                    tg_xxy_xxxxxyz_0[j] = pb_x * tg_xy_xxxxxyz_0[j] + wp_x[j] * tg_xy_xxxxxyz_1[j] + 0.5 * fl1_fx * tg_y_xxxxxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxxxyz_1[j] + 2.5 * fl1_fxn * tg_xy_xxxxyz_1[j];

                    tg_xxy_xxxxxzz_0[j] = pb_x * tg_xy_xxxxxzz_0[j] + wp_x[j] * tg_xy_xxxxxzz_1[j] + 0.5 * fl1_fx * tg_y_xxxxxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxxxzz_1[j] + 2.5 * fl1_fxn * tg_xy_xxxxzz_1[j];

                    tg_xxy_xxxxyyy_0[j] = pb_x * tg_xy_xxxxyyy_0[j] + wp_x[j] * tg_xy_xxxxyyy_1[j] + 0.5 * fl1_fx * tg_y_xxxxyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxxyyy_1[j] + 2.0 * fl1_fxn * tg_xy_xxxyyy_1[j];

                    tg_xxy_xxxxyyz_0[j] = pb_x * tg_xy_xxxxyyz_0[j] + wp_x[j] * tg_xy_xxxxyyz_1[j] + 0.5 * fl1_fx * tg_y_xxxxyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxxyyz_1[j] + 2.0 * fl1_fxn * tg_xy_xxxyyz_1[j];

                    tg_xxy_xxxxyzz_0[j] = pb_x * tg_xy_xxxxyzz_0[j] + wp_x[j] * tg_xy_xxxxyzz_1[j] + 0.5 * fl1_fx * tg_y_xxxxyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxxyzz_1[j] + 2.0 * fl1_fxn * tg_xy_xxxyzz_1[j];

                    tg_xxy_xxxxzzz_0[j] = pb_x * tg_xy_xxxxzzz_0[j] + wp_x[j] * tg_xy_xxxxzzz_1[j] + 0.5 * fl1_fx * tg_y_xxxxzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxxzzz_1[j] + 2.0 * fl1_fxn * tg_xy_xxxzzz_1[j];

                    tg_xxy_xxxyyyy_0[j] = pb_x * tg_xy_xxxyyyy_0[j] + wp_x[j] * tg_xy_xxxyyyy_1[j] + 0.5 * fl1_fx * tg_y_xxxyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxyyyy_1[j] + 1.5 * fl1_fxn * tg_xy_xxyyyy_1[j];

                    tg_xxy_xxxyyyz_0[j] = pb_x * tg_xy_xxxyyyz_0[j] + wp_x[j] * tg_xy_xxxyyyz_1[j] + 0.5 * fl1_fx * tg_y_xxxyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_xy_xxyyyz_1[j];

                    tg_xxy_xxxyyzz_0[j] = pb_x * tg_xy_xxxyyzz_0[j] + wp_x[j] * tg_xy_xxxyyzz_1[j] + 0.5 * fl1_fx * tg_y_xxxyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxyyzz_1[j] + 1.5 * fl1_fxn * tg_xy_xxyyzz_1[j];

                    tg_xxy_xxxyzzz_0[j] = pb_x * tg_xy_xxxyzzz_0[j] + wp_x[j] * tg_xy_xxxyzzz_1[j] + 0.5 * fl1_fx * tg_y_xxxyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxyzzz_1[j] + 1.5 * fl1_fxn * tg_xy_xxyzzz_1[j];

                    tg_xxy_xxxzzzz_0[j] = pb_x * tg_xy_xxxzzzz_0[j] + wp_x[j] * tg_xy_xxxzzzz_1[j] + 0.5 * fl1_fx * tg_y_xxxzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxzzzz_1[j] + 1.5 * fl1_fxn * tg_xy_xxzzzz_1[j];

                    tg_xxy_xxyyyyy_0[j] = pb_x * tg_xy_xxyyyyy_0[j] + wp_x[j] * tg_xy_xxyyyyy_1[j] + 0.5 * fl1_fx * tg_y_xxyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxyyyyy_1[j] + fl1_fxn * tg_xy_xyyyyy_1[j];

                    tg_xxy_xxyyyyz_0[j] = pb_x * tg_xy_xxyyyyz_0[j] + wp_x[j] * tg_xy_xxyyyyz_1[j] + 0.5 * fl1_fx * tg_y_xxyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxyyyyz_1[j] + fl1_fxn * tg_xy_xyyyyz_1[j];

                    tg_xxy_xxyyyzz_0[j] = pb_x * tg_xy_xxyyyzz_0[j] + wp_x[j] * tg_xy_xxyyyzz_1[j] + 0.5 * fl1_fx * tg_y_xxyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxyyyzz_1[j] + fl1_fxn * tg_xy_xyyyzz_1[j];

                    tg_xxy_xxyyzzz_0[j] = pb_x * tg_xy_xxyyzzz_0[j] + wp_x[j] * tg_xy_xxyyzzz_1[j] + 0.5 * fl1_fx * tg_y_xxyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxyyzzz_1[j] + fl1_fxn * tg_xy_xyyzzz_1[j];

                    tg_xxy_xxyzzzz_0[j] = pb_x * tg_xy_xxyzzzz_0[j] + wp_x[j] * tg_xy_xxyzzzz_1[j] + 0.5 * fl1_fx * tg_y_xxyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxyzzzz_1[j] + fl1_fxn * tg_xy_xyzzzz_1[j];

                    tg_xxy_xxzzzzz_0[j] = pb_x * tg_xy_xxzzzzz_0[j] + wp_x[j] * tg_xy_xxzzzzz_1[j] + 0.5 * fl1_fx * tg_y_xxzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxzzzzz_1[j] + fl1_fxn * tg_xy_xzzzzz_1[j];

                    tg_xxy_xyyyyyy_0[j] = pb_x * tg_xy_xyyyyyy_0[j] + wp_x[j] * tg_xy_xyyyyyy_1[j] + 0.5 * fl1_fx * tg_y_xyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xyyyyyy_1[j] + 0.5 * fl1_fxn * tg_xy_yyyyyy_1[j];

                    tg_xxy_xyyyyyz_0[j] = pb_x * tg_xy_xyyyyyz_0[j] + wp_x[j] * tg_xy_xyyyyyz_1[j] + 0.5 * fl1_fx * tg_y_xyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xyyyyyz_1[j] + 0.5 * fl1_fxn * tg_xy_yyyyyz_1[j];

                    tg_xxy_xyyyyzz_0[j] = pb_x * tg_xy_xyyyyzz_0[j] + wp_x[j] * tg_xy_xyyyyzz_1[j] + 0.5 * fl1_fx * tg_y_xyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xyyyyzz_1[j] + 0.5 * fl1_fxn * tg_xy_yyyyzz_1[j];

                    tg_xxy_xyyyzzz_0[j] = pb_x * tg_xy_xyyyzzz_0[j] + wp_x[j] * tg_xy_xyyyzzz_1[j] + 0.5 * fl1_fx * tg_y_xyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xyyyzzz_1[j] + 0.5 * fl1_fxn * tg_xy_yyyzzz_1[j];

                    tg_xxy_xyyzzzz_0[j] = pb_x * tg_xy_xyyzzzz_0[j] + wp_x[j] * tg_xy_xyyzzzz_1[j] + 0.5 * fl1_fx * tg_y_xyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xyyzzzz_1[j] + 0.5 * fl1_fxn * tg_xy_yyzzzz_1[j];

                    tg_xxy_xyzzzzz_0[j] = pb_x * tg_xy_xyzzzzz_0[j] + wp_x[j] * tg_xy_xyzzzzz_1[j] + 0.5 * fl1_fx * tg_y_xyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_xy_yzzzzz_1[j];

                    tg_xxy_xzzzzzz_0[j] = pb_x * tg_xy_xzzzzzz_0[j] + wp_x[j] * tg_xy_xzzzzzz_1[j] + 0.5 * fl1_fx * tg_y_xzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xy_zzzzzz_1[j];

                    tg_xxy_yyyyyyy_0[j] = pb_x * tg_xy_yyyyyyy_0[j] + wp_x[j] * tg_xy_yyyyyyy_1[j] + 0.5 * fl1_fx * tg_y_yyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_yyyyyyy_1[j];

                    tg_xxy_yyyyyyz_0[j] = pb_x * tg_xy_yyyyyyz_0[j] + wp_x[j] * tg_xy_yyyyyyz_1[j] + 0.5 * fl1_fx * tg_y_yyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_yyyyyyz_1[j];

                    tg_xxy_yyyyyzz_0[j] = pb_x * tg_xy_yyyyyzz_0[j] + wp_x[j] * tg_xy_yyyyyzz_1[j] + 0.5 * fl1_fx * tg_y_yyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_yyyyyzz_1[j];

                    tg_xxy_yyyyzzz_0[j] = pb_x * tg_xy_yyyyzzz_0[j] + wp_x[j] * tg_xy_yyyyzzz_1[j] + 0.5 * fl1_fx * tg_y_yyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_yyyyzzz_1[j];

                    tg_xxy_yyyzzzz_0[j] = pb_x * tg_xy_yyyzzzz_0[j] + wp_x[j] * tg_xy_yyyzzzz_1[j] + 0.5 * fl1_fx * tg_y_yyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_yyyzzzz_1[j];

                    tg_xxy_yyzzzzz_0[j] = pb_x * tg_xy_yyzzzzz_0[j] + wp_x[j] * tg_xy_yyzzzzz_1[j] + 0.5 * fl1_fx * tg_y_yyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_yyzzzzz_1[j];

                    tg_xxy_yzzzzzz_0[j] = pb_x * tg_xy_yzzzzzz_0[j] + wp_x[j] * tg_xy_yzzzzzz_1[j] + 0.5 * fl1_fx * tg_y_yzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_yzzzzzz_1[j];

                    tg_xxy_zzzzzzz_0[j] = pb_x * tg_xy_zzzzzzz_0[j] + wp_x[j] * tg_xy_zzzzzzz_1[j] + 0.5 * fl1_fx * tg_y_zzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_zzzzzzz_1[j];

                    tg_xxz_xxxxxxx_0[j] = pb_x * tg_xz_xxxxxxx_0[j] + wp_x[j] * tg_xz_xxxxxxx_1[j] + 0.5 * fl1_fx * tg_z_xxxxxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxxx_1[j] + 3.5 * fl1_fxn * tg_xz_xxxxxx_1[j];

                    tg_xxz_xxxxxxy_0[j] = pb_x * tg_xz_xxxxxxy_0[j] + wp_x[j] * tg_xz_xxxxxxy_1[j] + 0.5 * fl1_fx * tg_z_xxxxxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxxy_1[j] + 3.0 * fl1_fxn * tg_xz_xxxxxy_1[j];

                    tg_xxz_xxxxxxz_0[j] = pb_x * tg_xz_xxxxxxz_0[j] + wp_x[j] * tg_xz_xxxxxxz_1[j] + 0.5 * fl1_fx * tg_z_xxxxxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxxz_1[j] + 3.0 * fl1_fxn * tg_xz_xxxxxz_1[j];

                    tg_xxz_xxxxxyy_0[j] = pb_x * tg_xz_xxxxxyy_0[j] + wp_x[j] * tg_xz_xxxxxyy_1[j] + 0.5 * fl1_fx * tg_z_xxxxxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxyy_1[j] + 2.5 * fl1_fxn * tg_xz_xxxxyy_1[j];

                    tg_xxz_xxxxxyz_0[j] = pb_x * tg_xz_xxxxxyz_0[j] + wp_x[j] * tg_xz_xxxxxyz_1[j] + 0.5 * fl1_fx * tg_z_xxxxxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxyz_1[j] + 2.5 * fl1_fxn * tg_xz_xxxxyz_1[j];

                    tg_xxz_xxxxxzz_0[j] = pb_x * tg_xz_xxxxxzz_0[j] + wp_x[j] * tg_xz_xxxxxzz_1[j] + 0.5 * fl1_fx * tg_z_xxxxxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxzz_1[j] + 2.5 * fl1_fxn * tg_xz_xxxxzz_1[j];

                    tg_xxz_xxxxyyy_0[j] = pb_x * tg_xz_xxxxyyy_0[j] + wp_x[j] * tg_xz_xxxxyyy_1[j] + 0.5 * fl1_fx * tg_z_xxxxyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxyyy_1[j] + 2.0 * fl1_fxn * tg_xz_xxxyyy_1[j];

                    tg_xxz_xxxxyyz_0[j] = pb_x * tg_xz_xxxxyyz_0[j] + wp_x[j] * tg_xz_xxxxyyz_1[j] + 0.5 * fl1_fx * tg_z_xxxxyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxyyz_1[j] + 2.0 * fl1_fxn * tg_xz_xxxyyz_1[j];

                    tg_xxz_xxxxyzz_0[j] = pb_x * tg_xz_xxxxyzz_0[j] + wp_x[j] * tg_xz_xxxxyzz_1[j] + 0.5 * fl1_fx * tg_z_xxxxyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxyzz_1[j] + 2.0 * fl1_fxn * tg_xz_xxxyzz_1[j];

                    tg_xxz_xxxxzzz_0[j] = pb_x * tg_xz_xxxxzzz_0[j] + wp_x[j] * tg_xz_xxxxzzz_1[j] + 0.5 * fl1_fx * tg_z_xxxxzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxzzz_1[j] + 2.0 * fl1_fxn * tg_xz_xxxzzz_1[j];

                    tg_xxz_xxxyyyy_0[j] = pb_x * tg_xz_xxxyyyy_0[j] + wp_x[j] * tg_xz_xxxyyyy_1[j] + 0.5 * fl1_fx * tg_z_xxxyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxyyyy_1[j] + 1.5 * fl1_fxn * tg_xz_xxyyyy_1[j];

                    tg_xxz_xxxyyyz_0[j] = pb_x * tg_xz_xxxyyyz_0[j] + wp_x[j] * tg_xz_xxxyyyz_1[j] + 0.5 * fl1_fx * tg_z_xxxyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_xz_xxyyyz_1[j];

                    tg_xxz_xxxyyzz_0[j] = pb_x * tg_xz_xxxyyzz_0[j] + wp_x[j] * tg_xz_xxxyyzz_1[j] + 0.5 * fl1_fx * tg_z_xxxyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxyyzz_1[j] + 1.5 * fl1_fxn * tg_xz_xxyyzz_1[j];

                    tg_xxz_xxxyzzz_0[j] = pb_x * tg_xz_xxxyzzz_0[j] + wp_x[j] * tg_xz_xxxyzzz_1[j] + 0.5 * fl1_fx * tg_z_xxxyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxyzzz_1[j] + 1.5 * fl1_fxn * tg_xz_xxyzzz_1[j];

                    tg_xxz_xxxzzzz_0[j] = pb_x * tg_xz_xxxzzzz_0[j] + wp_x[j] * tg_xz_xxxzzzz_1[j] + 0.5 * fl1_fx * tg_z_xxxzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxzzzz_1[j] + 1.5 * fl1_fxn * tg_xz_xxzzzz_1[j];

                    tg_xxz_xxyyyyy_0[j] = pb_x * tg_xz_xxyyyyy_0[j] + wp_x[j] * tg_xz_xxyyyyy_1[j] + 0.5 * fl1_fx * tg_z_xxyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxyyyyy_1[j] + fl1_fxn * tg_xz_xyyyyy_1[j];

                    tg_xxz_xxyyyyz_0[j] = pb_x * tg_xz_xxyyyyz_0[j] + wp_x[j] * tg_xz_xxyyyyz_1[j] + 0.5 * fl1_fx * tg_z_xxyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxyyyyz_1[j] + fl1_fxn * tg_xz_xyyyyz_1[j];

                    tg_xxz_xxyyyzz_0[j] = pb_x * tg_xz_xxyyyzz_0[j] + wp_x[j] * tg_xz_xxyyyzz_1[j] + 0.5 * fl1_fx * tg_z_xxyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxyyyzz_1[j] + fl1_fxn * tg_xz_xyyyzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSFSK_90_180(      CMemBlock2D<double>& primBuffer,
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
                                             {3, -1, -1, -1},
                                             {7, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_3_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_3_7_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_2_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_2_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_1_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_1_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_2_6_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {6, -1, -1, -1}, 
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

                auto tg_xz_xxyyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 90); 

                auto tg_xz_xxyzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 91); 

                auto tg_xz_xxzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 92); 

                auto tg_xz_xyyyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 93); 

                auto tg_xz_xyyyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 94); 

                auto tg_xz_xyyyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 95); 

                auto tg_xz_xyyyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 96); 

                auto tg_xz_xyyzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 97); 

                auto tg_xz_xyzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 98); 

                auto tg_xz_xzzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 99); 

                auto tg_xz_yyyyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 100); 

                auto tg_xz_yyyyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 101); 

                auto tg_xz_yyyyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 102); 

                auto tg_xz_yyyyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 103); 

                auto tg_xz_yyyzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 104); 

                auto tg_xz_yyzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 105); 

                auto tg_xz_yzzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 106); 

                auto tg_xz_zzzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 107); 

                auto tg_yy_xxxxxxx_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 108); 

                auto tg_yy_xxxxxxy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 109); 

                auto tg_yy_xxxxxxz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 110); 

                auto tg_yy_xxxxxyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 111); 

                auto tg_yy_xxxxxyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 112); 

                auto tg_yy_xxxxxzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 113); 

                auto tg_yy_xxxxyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 114); 

                auto tg_yy_xxxxyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 115); 

                auto tg_yy_xxxxyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 116); 

                auto tg_yy_xxxxzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 117); 

                auto tg_yy_xxxyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 118); 

                auto tg_yy_xxxyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 119); 

                auto tg_yy_xxxyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 120); 

                auto tg_yy_xxxyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 121); 

                auto tg_yy_xxxzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 122); 

                auto tg_yy_xxyyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 123); 

                auto tg_yy_xxyyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 124); 

                auto tg_yy_xxyyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 125); 

                auto tg_yy_xxyyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 126); 

                auto tg_yy_xxyzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 127); 

                auto tg_yy_xxzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 128); 

                auto tg_yy_xyyyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 129); 

                auto tg_yy_xyyyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 130); 

                auto tg_yy_xyyyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 131); 

                auto tg_yy_xyyyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 132); 

                auto tg_yy_xyyzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 133); 

                auto tg_yy_xyzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 134); 

                auto tg_yy_xzzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 135); 

                auto tg_yy_yyyyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 136); 

                auto tg_yy_yyyyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 137); 

                auto tg_yy_yyyyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 138); 

                auto tg_yy_yyyyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 139); 

                auto tg_yy_yyyzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 140); 

                auto tg_yy_yyzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 141); 

                auto tg_yy_yzzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 142); 

                auto tg_yy_zzzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 143); 

                auto tg_yz_xxxxxxx_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 144); 

                auto tg_yz_xxxxxxy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 145); 

                auto tg_yz_xxxxxxz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 146); 

                auto tg_yz_xxxxxyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 147); 

                auto tg_yz_xxxxxyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 148); 

                auto tg_yz_xxxxxzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 149); 

                auto tg_yz_xxxxyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 150); 

                auto tg_yz_xxxxyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 151); 

                auto tg_yz_xxxxyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 152); 

                auto tg_yz_xxxxzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 153); 

                auto tg_yz_xxxyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 154); 

                auto tg_yz_xxxyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 155); 

                auto tg_yz_xxxyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 156); 

                auto tg_yz_xxxyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 157); 

                auto tg_yz_xxxzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 158); 

                auto tg_yz_xxyyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 159); 

                auto tg_yz_xxyyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 160); 

                auto tg_yz_xxyyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 161); 

                auto tg_yz_xxyyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 162); 

                auto tg_yz_xxyzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 163); 

                auto tg_yz_xxzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 164); 

                auto tg_yz_xyyyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 165); 

                auto tg_yz_xyyyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 166); 

                auto tg_yz_xyyyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 167); 

                auto tg_yz_xyyyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 168); 

                auto tg_yz_xyyzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 169); 

                auto tg_yz_xyzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 170); 

                auto tg_yz_xzzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 171); 

                auto tg_yz_yyyyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 172); 

                auto tg_yz_yyyyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 173); 

                auto tg_yz_yyyyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 174); 

                auto tg_yz_yyyyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 175); 

                auto tg_yz_yyyzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 176); 

                auto tg_yz_yyzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 177); 

                auto tg_yz_yzzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 178); 

                auto tg_yz_zzzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 179); 

                auto tg_xz_xxyyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 90); 

                auto tg_xz_xxyzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 91); 

                auto tg_xz_xxzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 92); 

                auto tg_xz_xyyyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 93); 

                auto tg_xz_xyyyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 94); 

                auto tg_xz_xyyyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 95); 

                auto tg_xz_xyyyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 96); 

                auto tg_xz_xyyzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 97); 

                auto tg_xz_xyzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 98); 

                auto tg_xz_xzzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 99); 

                auto tg_xz_yyyyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 100); 

                auto tg_xz_yyyyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 101); 

                auto tg_xz_yyyyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 102); 

                auto tg_xz_yyyyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 103); 

                auto tg_xz_yyyzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 104); 

                auto tg_xz_yyzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 105); 

                auto tg_xz_yzzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 106); 

                auto tg_xz_zzzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 107); 

                auto tg_yy_xxxxxxx_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 108); 

                auto tg_yy_xxxxxxy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 109); 

                auto tg_yy_xxxxxxz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 110); 

                auto tg_yy_xxxxxyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 111); 

                auto tg_yy_xxxxxyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 112); 

                auto tg_yy_xxxxxzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 113); 

                auto tg_yy_xxxxyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 114); 

                auto tg_yy_xxxxyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 115); 

                auto tg_yy_xxxxyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 116); 

                auto tg_yy_xxxxzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 117); 

                auto tg_yy_xxxyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 118); 

                auto tg_yy_xxxyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 119); 

                auto tg_yy_xxxyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 120); 

                auto tg_yy_xxxyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 121); 

                auto tg_yy_xxxzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 122); 

                auto tg_yy_xxyyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 123); 

                auto tg_yy_xxyyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 124); 

                auto tg_yy_xxyyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 125); 

                auto tg_yy_xxyyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 126); 

                auto tg_yy_xxyzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 127); 

                auto tg_yy_xxzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 128); 

                auto tg_yy_xyyyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 129); 

                auto tg_yy_xyyyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 130); 

                auto tg_yy_xyyyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 131); 

                auto tg_yy_xyyyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 132); 

                auto tg_yy_xyyzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 133); 

                auto tg_yy_xyzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 134); 

                auto tg_yy_xzzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 135); 

                auto tg_yy_yyyyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 136); 

                auto tg_yy_yyyyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 137); 

                auto tg_yy_yyyyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 138); 

                auto tg_yy_yyyyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 139); 

                auto tg_yy_yyyzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 140); 

                auto tg_yy_yyzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 141); 

                auto tg_yy_yzzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 142); 

                auto tg_yy_zzzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 143); 

                auto tg_yz_xxxxxxx_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 144); 

                auto tg_yz_xxxxxxy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 145); 

                auto tg_yz_xxxxxxz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 146); 

                auto tg_yz_xxxxxyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 147); 

                auto tg_yz_xxxxxyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 148); 

                auto tg_yz_xxxxxzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 149); 

                auto tg_yz_xxxxyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 150); 

                auto tg_yz_xxxxyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 151); 

                auto tg_yz_xxxxyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 152); 

                auto tg_yz_xxxxzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 153); 

                auto tg_yz_xxxyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 154); 

                auto tg_yz_xxxyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 155); 

                auto tg_yz_xxxyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 156); 

                auto tg_yz_xxxyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 157); 

                auto tg_yz_xxxzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 158); 

                auto tg_yz_xxyyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 159); 

                auto tg_yz_xxyyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 160); 

                auto tg_yz_xxyyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 161); 

                auto tg_yz_xxyyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 162); 

                auto tg_yz_xxyzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 163); 

                auto tg_yz_xxzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 164); 

                auto tg_yz_xyyyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 165); 

                auto tg_yz_xyyyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 166); 

                auto tg_yz_xyyyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 167); 

                auto tg_yz_xyyyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 168); 

                auto tg_yz_xyyzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 169); 

                auto tg_yz_xyzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 170); 

                auto tg_yz_xzzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 171); 

                auto tg_yz_yyyyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 172); 

                auto tg_yz_yyyyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 173); 

                auto tg_yz_yyyyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 174); 

                auto tg_yz_yyyyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 175); 

                auto tg_yz_yyyzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 176); 

                auto tg_yz_yyzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 177); 

                auto tg_yz_yzzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 178); 

                auto tg_yz_zzzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 179); 

                auto tg_z_xxyyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 90); 

                auto tg_z_xxyzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 91); 

                auto tg_z_xxzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 92); 

                auto tg_z_xyyyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 93); 

                auto tg_z_xyyyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 94); 

                auto tg_z_xyyyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 95); 

                auto tg_z_xyyyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 96); 

                auto tg_z_xyyzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 97); 

                auto tg_z_xyzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 98); 

                auto tg_z_xzzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 99); 

                auto tg_z_yyyyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 100); 

                auto tg_z_yyyyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 101); 

                auto tg_z_yyyyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 102); 

                auto tg_z_yyyyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 103); 

                auto tg_z_yyyzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 104); 

                auto tg_z_yyzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 105); 

                auto tg_z_yzzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 106); 

                auto tg_z_zzzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 107); 

                auto tg_z_xxyyzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 90); 

                auto tg_z_xxyzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 91); 

                auto tg_z_xxzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 92); 

                auto tg_z_xyyyyyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 93); 

                auto tg_z_xyyyyyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 94); 

                auto tg_z_xyyyyzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 95); 

                auto tg_z_xyyyzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 96); 

                auto tg_z_xyyzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 97); 

                auto tg_z_xyzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 98); 

                auto tg_z_xzzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 99); 

                auto tg_z_yyyyyyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 100); 

                auto tg_z_yyyyyyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 101); 

                auto tg_z_yyyyyzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 102); 

                auto tg_z_yyyyzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 103); 

                auto tg_z_yyyzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 104); 

                auto tg_z_yyzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 105); 

                auto tg_z_yzzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 106); 

                auto tg_z_zzzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 107); 

                auto tg_xz_xyyzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 74); 

                auto tg_xz_xyzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 75); 

                auto tg_xz_xzzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 76); 

                auto tg_xz_yyyyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 77); 

                auto tg_xz_yyyyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 78); 

                auto tg_xz_yyyyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 79); 

                auto tg_xz_yyyzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 80); 

                auto tg_xz_yyzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 81); 

                auto tg_xz_yzzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 82); 

                auto tg_xz_zzzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 83); 

                auto tg_yy_xxxxxx_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 84); 

                auto tg_yy_xxxxxy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 85); 

                auto tg_yy_xxxxxz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 86); 

                auto tg_yy_xxxxyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 87); 

                auto tg_yy_xxxxyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 88); 

                auto tg_yy_xxxxzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 89); 

                auto tg_yy_xxxyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 90); 

                auto tg_yy_xxxyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 91); 

                auto tg_yy_xxxyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 92); 

                auto tg_yy_xxxzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 93); 

                auto tg_yy_xxyyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 94); 

                auto tg_yy_xxyyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 95); 

                auto tg_yy_xxyyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 96); 

                auto tg_yy_xxyzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 97); 

                auto tg_yy_xxzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 98); 

                auto tg_yy_xyyyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 99); 

                auto tg_yy_xyyyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 100); 

                auto tg_yy_xyyyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 101); 

                auto tg_yy_xyyzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 102); 

                auto tg_yy_xyzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 103); 

                auto tg_yy_xzzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 104); 

                auto tg_yy_yyyyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 105); 

                auto tg_yy_yyyyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 106); 

                auto tg_yy_yyyyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 107); 

                auto tg_yy_yyyzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 108); 

                auto tg_yy_yyzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 109); 

                auto tg_yy_yzzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 110); 

                auto tg_yy_zzzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 111); 

                auto tg_yz_xxxxxx_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 112); 

                auto tg_yz_xxxxxy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 113); 

                auto tg_yz_xxxxxz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 114); 

                auto tg_yz_xxxxyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 115); 

                auto tg_yz_xxxxyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 116); 

                auto tg_yz_xxxxzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 117); 

                auto tg_yz_xxxyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 118); 

                auto tg_yz_xxxyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 119); 

                auto tg_yz_xxxyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 120); 

                auto tg_yz_xxxzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 121); 

                auto tg_yz_xxyyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 122); 

                auto tg_yz_xxyyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 123); 

                auto tg_yz_xxyyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 124); 

                auto tg_yz_xxyzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 125); 

                auto tg_yz_xxzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 126); 

                auto tg_yz_xyyyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 127); 

                auto tg_yz_xyyyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 128); 

                auto tg_yz_xyyyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 129); 

                auto tg_yz_xyyzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 130); 

                auto tg_yz_xyzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 131); 

                auto tg_yz_xzzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 132); 

                auto tg_yz_yyyyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 133); 

                auto tg_yz_yyyyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 134); 

                auto tg_yz_yyyyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 135); 

                auto tg_yz_yyyzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 136); 

                auto tg_yz_yyzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 137); 

                auto tg_yz_yzzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 138); 

                auto tg_yz_zzzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 139); 

                // set up pointers to integrals

                auto tg_xxz_xxyyzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 90); 

                auto tg_xxz_xxyzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 91); 

                auto tg_xxz_xxzzzzz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 92); 

                auto tg_xxz_xyyyyyy_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 93); 

                auto tg_xxz_xyyyyyz_0 = primBuffer.data(pidx_g_3_7_m0 + 360 * idx + 94); 

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

                // Batch of Integrals (90,180)

                #pragma omp simd aligned(fxn, fza, tg_xxz_xxyyzzz_0, tg_xxz_xxyzzzz_0, tg_xxz_xxzzzzz_0, \
                                         tg_xxz_xyyyyyy_0, tg_xxz_xyyyyyz_0, tg_xxz_xyyyyzz_0, tg_xxz_xyyyzzz_0, \
                                         tg_xxz_xyyzzzz_0, tg_xxz_xyzzzzz_0, tg_xxz_xzzzzzz_0, tg_xxz_yyyyyyy_0, \
                                         tg_xxz_yyyyyyz_0, tg_xxz_yyyyyzz_0, tg_xxz_yyyyzzz_0, tg_xxz_yyyzzzz_0, \
                                         tg_xxz_yyzzzzz_0, tg_xxz_yzzzzzz_0, tg_xxz_zzzzzzz_0, tg_xyy_xxxxxxx_0, \
                                         tg_xyy_xxxxxxy_0, tg_xyy_xxxxxxz_0, tg_xyy_xxxxxyy_0, tg_xyy_xxxxxyz_0, \
                                         tg_xyy_xxxxxzz_0, tg_xyy_xxxxyyy_0, tg_xyy_xxxxyyz_0, tg_xyy_xxxxyzz_0, \
                                         tg_xyy_xxxxzzz_0, tg_xyy_xxxyyyy_0, tg_xyy_xxxyyyz_0, tg_xyy_xxxyyzz_0, \
                                         tg_xyy_xxxyzzz_0, tg_xyy_xxxzzzz_0, tg_xyy_xxyyyyy_0, tg_xyy_xxyyyyz_0, \
                                         tg_xyy_xxyyyzz_0, tg_xyy_xxyyzzz_0, tg_xyy_xxyzzzz_0, tg_xyy_xxzzzzz_0, \
                                         tg_xyy_xyyyyyy_0, tg_xyy_xyyyyyz_0, tg_xyy_xyyyyzz_0, tg_xyy_xyyyzzz_0, \
                                         tg_xyy_xyyzzzz_0, tg_xyy_xyzzzzz_0, tg_xyy_xzzzzzz_0, tg_xyy_yyyyyyy_0, \
                                         tg_xyy_yyyyyyz_0, tg_xyy_yyyyyzz_0, tg_xyy_yyyyzzz_0, tg_xyy_yyyzzzz_0, \
                                         tg_xyy_yyzzzzz_0, tg_xyy_yzzzzzz_0, tg_xyy_zzzzzzz_0, tg_xyz_xxxxxxx_0, \
                                         tg_xyz_xxxxxxy_0, tg_xyz_xxxxxxz_0, tg_xyz_xxxxxyy_0, tg_xyz_xxxxxyz_0, \
                                         tg_xyz_xxxxxzz_0, tg_xyz_xxxxyyy_0, tg_xyz_xxxxyyz_0, tg_xyz_xxxxyzz_0, \
                                         tg_xyz_xxxxzzz_0, tg_xyz_xxxyyyy_0, tg_xyz_xxxyyyz_0, tg_xyz_xxxyyzz_0, \
                                         tg_xyz_xxxyzzz_0, tg_xyz_xxxzzzz_0, tg_xyz_xxyyyyy_0, tg_xyz_xxyyyyz_0, \
                                         tg_xyz_xxyyyzz_0, tg_xyz_xxyyzzz_0, tg_xyz_xxyzzzz_0, tg_xyz_xxzzzzz_0, \
                                         tg_xyz_xyyyyyy_0, tg_xyz_xyyyyyz_0, tg_xyz_xyyyyzz_0, tg_xyz_xyyyzzz_0, \
                                         tg_xyz_xyyzzzz_0, tg_xyz_xyzzzzz_0, tg_xyz_xzzzzzz_0, tg_xyz_yyyyyyy_0, \
                                         tg_xyz_yyyyyyz_0, tg_xyz_yyyyyzz_0, tg_xyz_yyyyzzz_0, tg_xyz_yyyzzzz_0, \
                                         tg_xyz_yyzzzzz_0, tg_xyz_yzzzzzz_0, tg_xyz_zzzzzzz_0, tg_xz_xxyyzzz_0, \
                                         tg_xz_xxyyzzz_1, tg_xz_xxyzzzz_0, tg_xz_xxyzzzz_1, tg_xz_xxzzzzz_0, tg_xz_xxzzzzz_1, \
                                         tg_xz_xyyyyyy_0, tg_xz_xyyyyyy_1, tg_xz_xyyyyyz_0, tg_xz_xyyyyyz_1, tg_xz_xyyyyzz_0, \
                                         tg_xz_xyyyyzz_1, tg_xz_xyyyzzz_0, tg_xz_xyyyzzz_1, tg_xz_xyyzzz_1, tg_xz_xyyzzzz_0, \
                                         tg_xz_xyyzzzz_1, tg_xz_xyzzzz_1, tg_xz_xyzzzzz_0, tg_xz_xyzzzzz_1, tg_xz_xzzzzz_1, \
                                         tg_xz_xzzzzzz_0, tg_xz_xzzzzzz_1, tg_xz_yyyyyy_1, tg_xz_yyyyyyy_0, tg_xz_yyyyyyy_1, \
                                         tg_xz_yyyyyyz_0, tg_xz_yyyyyyz_1, tg_xz_yyyyyz_1, tg_xz_yyyyyzz_0, tg_xz_yyyyyzz_1, \
                                         tg_xz_yyyyzz_1, tg_xz_yyyyzzz_0, tg_xz_yyyyzzz_1, tg_xz_yyyzzz_1, tg_xz_yyyzzzz_0, \
                                         tg_xz_yyyzzzz_1, tg_xz_yyzzzz_1, tg_xz_yyzzzzz_0, tg_xz_yyzzzzz_1, tg_xz_yzzzzz_1, \
                                         tg_xz_yzzzzzz_0, tg_xz_yzzzzzz_1, tg_xz_zzzzzz_1, tg_xz_zzzzzzz_0, tg_xz_zzzzzzz_1, \
                                         tg_yy_xxxxxx_1, tg_yy_xxxxxxx_0, tg_yy_xxxxxxx_1, tg_yy_xxxxxxy_0, tg_yy_xxxxxxy_1, \
                                         tg_yy_xxxxxxz_0, tg_yy_xxxxxxz_1, tg_yy_xxxxxy_1, tg_yy_xxxxxyy_0, tg_yy_xxxxxyy_1, \
                                         tg_yy_xxxxxyz_0, tg_yy_xxxxxyz_1, tg_yy_xxxxxz_1, tg_yy_xxxxxzz_0, tg_yy_xxxxxzz_1, \
                                         tg_yy_xxxxyy_1, tg_yy_xxxxyyy_0, tg_yy_xxxxyyy_1, tg_yy_xxxxyyz_0, tg_yy_xxxxyyz_1, \
                                         tg_yy_xxxxyz_1, tg_yy_xxxxyzz_0, tg_yy_xxxxyzz_1, tg_yy_xxxxzz_1, tg_yy_xxxxzzz_0, \
                                         tg_yy_xxxxzzz_1, tg_yy_xxxyyy_1, tg_yy_xxxyyyy_0, tg_yy_xxxyyyy_1, tg_yy_xxxyyyz_0, \
                                         tg_yy_xxxyyyz_1, tg_yy_xxxyyz_1, tg_yy_xxxyyzz_0, tg_yy_xxxyyzz_1, tg_yy_xxxyzz_1, \
                                         tg_yy_xxxyzzz_0, tg_yy_xxxyzzz_1, tg_yy_xxxzzz_1, tg_yy_xxxzzzz_0, tg_yy_xxxzzzz_1, \
                                         tg_yy_xxyyyy_1, tg_yy_xxyyyyy_0, tg_yy_xxyyyyy_1, tg_yy_xxyyyyz_0, tg_yy_xxyyyyz_1, \
                                         tg_yy_xxyyyz_1, tg_yy_xxyyyzz_0, tg_yy_xxyyyzz_1, tg_yy_xxyyzz_1, tg_yy_xxyyzzz_0, \
                                         tg_yy_xxyyzzz_1, tg_yy_xxyzzz_1, tg_yy_xxyzzzz_0, tg_yy_xxyzzzz_1, tg_yy_xxzzzz_1, \
                                         tg_yy_xxzzzzz_0, tg_yy_xxzzzzz_1, tg_yy_xyyyyy_1, tg_yy_xyyyyyy_0, tg_yy_xyyyyyy_1, \
                                         tg_yy_xyyyyyz_0, tg_yy_xyyyyyz_1, tg_yy_xyyyyz_1, tg_yy_xyyyyzz_0, tg_yy_xyyyyzz_1, \
                                         tg_yy_xyyyzz_1, tg_yy_xyyyzzz_0, tg_yy_xyyyzzz_1, tg_yy_xyyzzz_1, tg_yy_xyyzzzz_0, \
                                         tg_yy_xyyzzzz_1, tg_yy_xyzzzz_1, tg_yy_xyzzzzz_0, tg_yy_xyzzzzz_1, tg_yy_xzzzzz_1, \
                                         tg_yy_xzzzzzz_0, tg_yy_xzzzzzz_1, tg_yy_yyyyyy_1, tg_yy_yyyyyyy_0, tg_yy_yyyyyyy_1, \
                                         tg_yy_yyyyyyz_0, tg_yy_yyyyyyz_1, tg_yy_yyyyyz_1, tg_yy_yyyyyzz_0, tg_yy_yyyyyzz_1, \
                                         tg_yy_yyyyzz_1, tg_yy_yyyyzzz_0, tg_yy_yyyyzzz_1, tg_yy_yyyzzz_1, tg_yy_yyyzzzz_0, \
                                         tg_yy_yyyzzzz_1, tg_yy_yyzzzz_1, tg_yy_yyzzzzz_0, tg_yy_yyzzzzz_1, tg_yy_yzzzzz_1, \
                                         tg_yy_yzzzzzz_0, tg_yy_yzzzzzz_1, tg_yy_zzzzzz_1, tg_yy_zzzzzzz_0, tg_yy_zzzzzzz_1, \
                                         tg_yz_xxxxxx_1, tg_yz_xxxxxxx_0, tg_yz_xxxxxxx_1, tg_yz_xxxxxxy_0, tg_yz_xxxxxxy_1, \
                                         tg_yz_xxxxxxz_0, tg_yz_xxxxxxz_1, tg_yz_xxxxxy_1, tg_yz_xxxxxyy_0, tg_yz_xxxxxyy_1, \
                                         tg_yz_xxxxxyz_0, tg_yz_xxxxxyz_1, tg_yz_xxxxxz_1, tg_yz_xxxxxzz_0, tg_yz_xxxxxzz_1, \
                                         tg_yz_xxxxyy_1, tg_yz_xxxxyyy_0, tg_yz_xxxxyyy_1, tg_yz_xxxxyyz_0, tg_yz_xxxxyyz_1, \
                                         tg_yz_xxxxyz_1, tg_yz_xxxxyzz_0, tg_yz_xxxxyzz_1, tg_yz_xxxxzz_1, tg_yz_xxxxzzz_0, \
                                         tg_yz_xxxxzzz_1, tg_yz_xxxyyy_1, tg_yz_xxxyyyy_0, tg_yz_xxxyyyy_1, tg_yz_xxxyyyz_0, \
                                         tg_yz_xxxyyyz_1, tg_yz_xxxyyz_1, tg_yz_xxxyyzz_0, tg_yz_xxxyyzz_1, tg_yz_xxxyzz_1, \
                                         tg_yz_xxxyzzz_0, tg_yz_xxxyzzz_1, tg_yz_xxxzzz_1, tg_yz_xxxzzzz_0, tg_yz_xxxzzzz_1, \
                                         tg_yz_xxyyyy_1, tg_yz_xxyyyyy_0, tg_yz_xxyyyyy_1, tg_yz_xxyyyyz_0, tg_yz_xxyyyyz_1, \
                                         tg_yz_xxyyyz_1, tg_yz_xxyyyzz_0, tg_yz_xxyyyzz_1, tg_yz_xxyyzz_1, tg_yz_xxyyzzz_0, \
                                         tg_yz_xxyyzzz_1, tg_yz_xxyzzz_1, tg_yz_xxyzzzz_0, tg_yz_xxyzzzz_1, tg_yz_xxzzzz_1, \
                                         tg_yz_xxzzzzz_0, tg_yz_xxzzzzz_1, tg_yz_xyyyyy_1, tg_yz_xyyyyyy_0, tg_yz_xyyyyyy_1, \
                                         tg_yz_xyyyyyz_0, tg_yz_xyyyyyz_1, tg_yz_xyyyyz_1, tg_yz_xyyyyzz_0, tg_yz_xyyyyzz_1, \
                                         tg_yz_xyyyzz_1, tg_yz_xyyyzzz_0, tg_yz_xyyyzzz_1, tg_yz_xyyzzz_1, tg_yz_xyyzzzz_0, \
                                         tg_yz_xyyzzzz_1, tg_yz_xyzzzz_1, tg_yz_xyzzzzz_0, tg_yz_xyzzzzz_1, tg_yz_xzzzzz_1, \
                                         tg_yz_xzzzzzz_0, tg_yz_xzzzzzz_1, tg_yz_yyyyyy_1, tg_yz_yyyyyyy_0, tg_yz_yyyyyyy_1, \
                                         tg_yz_yyyyyyz_0, tg_yz_yyyyyyz_1, tg_yz_yyyyyz_1, tg_yz_yyyyyzz_0, tg_yz_yyyyyzz_1, \
                                         tg_yz_yyyyzz_1, tg_yz_yyyyzzz_0, tg_yz_yyyyzzz_1, tg_yz_yyyzzz_1, tg_yz_yyyzzzz_0, \
                                         tg_yz_yyyzzzz_1, tg_yz_yyzzzz_1, tg_yz_yyzzzzz_0, tg_yz_yyzzzzz_1, tg_yz_yzzzzz_1, \
                                         tg_yz_yzzzzzz_0, tg_yz_yzzzzzz_1, tg_yz_zzzzzz_1, tg_yz_zzzzzzz_0, tg_yz_zzzzzzz_1, \
                                         tg_z_xxyyzzz_0, tg_z_xxyyzzz_1, tg_z_xxyzzzz_0, tg_z_xxyzzzz_1, tg_z_xxzzzzz_0, \
                                         tg_z_xxzzzzz_1, tg_z_xyyyyyy_0, tg_z_xyyyyyy_1, tg_z_xyyyyyz_0, tg_z_xyyyyyz_1, \
                                         tg_z_xyyyyzz_0, tg_z_xyyyyzz_1, tg_z_xyyyzzz_0, tg_z_xyyyzzz_1, tg_z_xyyzzzz_0, \
                                         tg_z_xyyzzzz_1, tg_z_xyzzzzz_0, tg_z_xyzzzzz_1, tg_z_xzzzzzz_0, tg_z_xzzzzzz_1, \
                                         tg_z_yyyyyyy_0, tg_z_yyyyyyy_1, tg_z_yyyyyyz_0, tg_z_yyyyyyz_1, tg_z_yyyyyzz_0, \
                                         tg_z_yyyyyzz_1, tg_z_yyyyzzz_0, tg_z_yyyyzzz_1, tg_z_yyyzzzz_0, tg_z_yyyzzzz_1, \
                                         tg_z_yyzzzzz_0, tg_z_yyzzzzz_1, tg_z_yzzzzzz_0, tg_z_yzzzzzz_1, tg_z_zzzzzzz_0, \
                                         tg_z_zzzzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxz_xxyyzzz_0[j] = pb_x * tg_xz_xxyyzzz_0[j] + wp_x[j] * tg_xz_xxyyzzz_1[j] + 0.5 * fl1_fx * tg_z_xxyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxyyzzz_1[j] + fl1_fxn * tg_xz_xyyzzz_1[j];

                    tg_xxz_xxyzzzz_0[j] = pb_x * tg_xz_xxyzzzz_0[j] + wp_x[j] * tg_xz_xxyzzzz_1[j] + 0.5 * fl1_fx * tg_z_xxyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxyzzzz_1[j] + fl1_fxn * tg_xz_xyzzzz_1[j];

                    tg_xxz_xxzzzzz_0[j] = pb_x * tg_xz_xxzzzzz_0[j] + wp_x[j] * tg_xz_xxzzzzz_1[j] + 0.5 * fl1_fx * tg_z_xxzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxzzzzz_1[j] + fl1_fxn * tg_xz_xzzzzz_1[j];

                    tg_xxz_xyyyyyy_0[j] = pb_x * tg_xz_xyyyyyy_0[j] + wp_x[j] * tg_xz_xyyyyyy_1[j] + 0.5 * fl1_fx * tg_z_xyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyyyyyy_1[j] + 0.5 * fl1_fxn * tg_xz_yyyyyy_1[j];

                    tg_xxz_xyyyyyz_0[j] = pb_x * tg_xz_xyyyyyz_0[j] + wp_x[j] * tg_xz_xyyyyyz_1[j] + 0.5 * fl1_fx * tg_z_xyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyyyyyz_1[j] + 0.5 * fl1_fxn * tg_xz_yyyyyz_1[j];

                    tg_xxz_xyyyyzz_0[j] = pb_x * tg_xz_xyyyyzz_0[j] + wp_x[j] * tg_xz_xyyyyzz_1[j] + 0.5 * fl1_fx * tg_z_xyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyyyyzz_1[j] + 0.5 * fl1_fxn * tg_xz_yyyyzz_1[j];

                    tg_xxz_xyyyzzz_0[j] = pb_x * tg_xz_xyyyzzz_0[j] + wp_x[j] * tg_xz_xyyyzzz_1[j] + 0.5 * fl1_fx * tg_z_xyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyyyzzz_1[j] + 0.5 * fl1_fxn * tg_xz_yyyzzz_1[j];

                    tg_xxz_xyyzzzz_0[j] = pb_x * tg_xz_xyyzzzz_0[j] + wp_x[j] * tg_xz_xyyzzzz_1[j] + 0.5 * fl1_fx * tg_z_xyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyyzzzz_1[j] + 0.5 * fl1_fxn * tg_xz_yyzzzz_1[j];

                    tg_xxz_xyzzzzz_0[j] = pb_x * tg_xz_xyzzzzz_0[j] + wp_x[j] * tg_xz_xyzzzzz_1[j] + 0.5 * fl1_fx * tg_z_xyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_xz_yzzzzz_1[j];

                    tg_xxz_xzzzzzz_0[j] = pb_x * tg_xz_xzzzzzz_0[j] + wp_x[j] * tg_xz_xzzzzzz_1[j] + 0.5 * fl1_fx * tg_z_xzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xz_zzzzzz_1[j];

                    tg_xxz_yyyyyyy_0[j] = pb_x * tg_xz_yyyyyyy_0[j] + wp_x[j] * tg_xz_yyyyyyy_1[j] + 0.5 * fl1_fx * tg_z_yyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyyyyyy_1[j];

                    tg_xxz_yyyyyyz_0[j] = pb_x * tg_xz_yyyyyyz_0[j] + wp_x[j] * tg_xz_yyyyyyz_1[j] + 0.5 * fl1_fx * tg_z_yyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyyyyyz_1[j];

                    tg_xxz_yyyyyzz_0[j] = pb_x * tg_xz_yyyyyzz_0[j] + wp_x[j] * tg_xz_yyyyyzz_1[j] + 0.5 * fl1_fx * tg_z_yyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyyyyzz_1[j];

                    tg_xxz_yyyyzzz_0[j] = pb_x * tg_xz_yyyyzzz_0[j] + wp_x[j] * tg_xz_yyyyzzz_1[j] + 0.5 * fl1_fx * tg_z_yyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyyyzzz_1[j];

                    tg_xxz_yyyzzzz_0[j] = pb_x * tg_xz_yyyzzzz_0[j] + wp_x[j] * tg_xz_yyyzzzz_1[j] + 0.5 * fl1_fx * tg_z_yyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyyzzzz_1[j];

                    tg_xxz_yyzzzzz_0[j] = pb_x * tg_xz_yyzzzzz_0[j] + wp_x[j] * tg_xz_yyzzzzz_1[j] + 0.5 * fl1_fx * tg_z_yyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyzzzzz_1[j];

                    tg_xxz_yzzzzzz_0[j] = pb_x * tg_xz_yzzzzzz_0[j] + wp_x[j] * tg_xz_yzzzzzz_1[j] + 0.5 * fl1_fx * tg_z_yzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yzzzzzz_1[j];

                    tg_xxz_zzzzzzz_0[j] = pb_x * tg_xz_zzzzzzz_0[j] + wp_x[j] * tg_xz_zzzzzzz_1[j] + 0.5 * fl1_fx * tg_z_zzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_zzzzzzz_1[j];

                    tg_xyy_xxxxxxx_0[j] = pb_x * tg_yy_xxxxxxx_0[j] + wp_x[j] * tg_yy_xxxxxxx_1[j] + 3.5 * fl1_fxn * tg_yy_xxxxxx_1[j];

                    tg_xyy_xxxxxxy_0[j] = pb_x * tg_yy_xxxxxxy_0[j] + wp_x[j] * tg_yy_xxxxxxy_1[j] + 3.0 * fl1_fxn * tg_yy_xxxxxy_1[j];

                    tg_xyy_xxxxxxz_0[j] = pb_x * tg_yy_xxxxxxz_0[j] + wp_x[j] * tg_yy_xxxxxxz_1[j] + 3.0 * fl1_fxn * tg_yy_xxxxxz_1[j];

                    tg_xyy_xxxxxyy_0[j] = pb_x * tg_yy_xxxxxyy_0[j] + wp_x[j] * tg_yy_xxxxxyy_1[j] + 2.5 * fl1_fxn * tg_yy_xxxxyy_1[j];

                    tg_xyy_xxxxxyz_0[j] = pb_x * tg_yy_xxxxxyz_0[j] + wp_x[j] * tg_yy_xxxxxyz_1[j] + 2.5 * fl1_fxn * tg_yy_xxxxyz_1[j];

                    tg_xyy_xxxxxzz_0[j] = pb_x * tg_yy_xxxxxzz_0[j] + wp_x[j] * tg_yy_xxxxxzz_1[j] + 2.5 * fl1_fxn * tg_yy_xxxxzz_1[j];

                    tg_xyy_xxxxyyy_0[j] = pb_x * tg_yy_xxxxyyy_0[j] + wp_x[j] * tg_yy_xxxxyyy_1[j] + 2.0 * fl1_fxn * tg_yy_xxxyyy_1[j];

                    tg_xyy_xxxxyyz_0[j] = pb_x * tg_yy_xxxxyyz_0[j] + wp_x[j] * tg_yy_xxxxyyz_1[j] + 2.0 * fl1_fxn * tg_yy_xxxyyz_1[j];

                    tg_xyy_xxxxyzz_0[j] = pb_x * tg_yy_xxxxyzz_0[j] + wp_x[j] * tg_yy_xxxxyzz_1[j] + 2.0 * fl1_fxn * tg_yy_xxxyzz_1[j];

                    tg_xyy_xxxxzzz_0[j] = pb_x * tg_yy_xxxxzzz_0[j] + wp_x[j] * tg_yy_xxxxzzz_1[j] + 2.0 * fl1_fxn * tg_yy_xxxzzz_1[j];

                    tg_xyy_xxxyyyy_0[j] = pb_x * tg_yy_xxxyyyy_0[j] + wp_x[j] * tg_yy_xxxyyyy_1[j] + 1.5 * fl1_fxn * tg_yy_xxyyyy_1[j];

                    tg_xyy_xxxyyyz_0[j] = pb_x * tg_yy_xxxyyyz_0[j] + wp_x[j] * tg_yy_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_yy_xxyyyz_1[j];

                    tg_xyy_xxxyyzz_0[j] = pb_x * tg_yy_xxxyyzz_0[j] + wp_x[j] * tg_yy_xxxyyzz_1[j] + 1.5 * fl1_fxn * tg_yy_xxyyzz_1[j];

                    tg_xyy_xxxyzzz_0[j] = pb_x * tg_yy_xxxyzzz_0[j] + wp_x[j] * tg_yy_xxxyzzz_1[j] + 1.5 * fl1_fxn * tg_yy_xxyzzz_1[j];

                    tg_xyy_xxxzzzz_0[j] = pb_x * tg_yy_xxxzzzz_0[j] + wp_x[j] * tg_yy_xxxzzzz_1[j] + 1.5 * fl1_fxn * tg_yy_xxzzzz_1[j];

                    tg_xyy_xxyyyyy_0[j] = pb_x * tg_yy_xxyyyyy_0[j] + wp_x[j] * tg_yy_xxyyyyy_1[j] + fl1_fxn * tg_yy_xyyyyy_1[j];

                    tg_xyy_xxyyyyz_0[j] = pb_x * tg_yy_xxyyyyz_0[j] + wp_x[j] * tg_yy_xxyyyyz_1[j] + fl1_fxn * tg_yy_xyyyyz_1[j];

                    tg_xyy_xxyyyzz_0[j] = pb_x * tg_yy_xxyyyzz_0[j] + wp_x[j] * tg_yy_xxyyyzz_1[j] + fl1_fxn * tg_yy_xyyyzz_1[j];

                    tg_xyy_xxyyzzz_0[j] = pb_x * tg_yy_xxyyzzz_0[j] + wp_x[j] * tg_yy_xxyyzzz_1[j] + fl1_fxn * tg_yy_xyyzzz_1[j];

                    tg_xyy_xxyzzzz_0[j] = pb_x * tg_yy_xxyzzzz_0[j] + wp_x[j] * tg_yy_xxyzzzz_1[j] + fl1_fxn * tg_yy_xyzzzz_1[j];

                    tg_xyy_xxzzzzz_0[j] = pb_x * tg_yy_xxzzzzz_0[j] + wp_x[j] * tg_yy_xxzzzzz_1[j] + fl1_fxn * tg_yy_xzzzzz_1[j];

                    tg_xyy_xyyyyyy_0[j] = pb_x * tg_yy_xyyyyyy_0[j] + wp_x[j] * tg_yy_xyyyyyy_1[j] + 0.5 * fl1_fxn * tg_yy_yyyyyy_1[j];

                    tg_xyy_xyyyyyz_0[j] = pb_x * tg_yy_xyyyyyz_0[j] + wp_x[j] * tg_yy_xyyyyyz_1[j] + 0.5 * fl1_fxn * tg_yy_yyyyyz_1[j];

                    tg_xyy_xyyyyzz_0[j] = pb_x * tg_yy_xyyyyzz_0[j] + wp_x[j] * tg_yy_xyyyyzz_1[j] + 0.5 * fl1_fxn * tg_yy_yyyyzz_1[j];

                    tg_xyy_xyyyzzz_0[j] = pb_x * tg_yy_xyyyzzz_0[j] + wp_x[j] * tg_yy_xyyyzzz_1[j] + 0.5 * fl1_fxn * tg_yy_yyyzzz_1[j];

                    tg_xyy_xyyzzzz_0[j] = pb_x * tg_yy_xyyzzzz_0[j] + wp_x[j] * tg_yy_xyyzzzz_1[j] + 0.5 * fl1_fxn * tg_yy_yyzzzz_1[j];

                    tg_xyy_xyzzzzz_0[j] = pb_x * tg_yy_xyzzzzz_0[j] + wp_x[j] * tg_yy_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_yy_yzzzzz_1[j];

                    tg_xyy_xzzzzzz_0[j] = pb_x * tg_yy_xzzzzzz_0[j] + wp_x[j] * tg_yy_xzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yy_zzzzzz_1[j];

                    tg_xyy_yyyyyyy_0[j] = pb_x * tg_yy_yyyyyyy_0[j] + wp_x[j] * tg_yy_yyyyyyy_1[j];

                    tg_xyy_yyyyyyz_0[j] = pb_x * tg_yy_yyyyyyz_0[j] + wp_x[j] * tg_yy_yyyyyyz_1[j];

                    tg_xyy_yyyyyzz_0[j] = pb_x * tg_yy_yyyyyzz_0[j] + wp_x[j] * tg_yy_yyyyyzz_1[j];

                    tg_xyy_yyyyzzz_0[j] = pb_x * tg_yy_yyyyzzz_0[j] + wp_x[j] * tg_yy_yyyyzzz_1[j];

                    tg_xyy_yyyzzzz_0[j] = pb_x * tg_yy_yyyzzzz_0[j] + wp_x[j] * tg_yy_yyyzzzz_1[j];

                    tg_xyy_yyzzzzz_0[j] = pb_x * tg_yy_yyzzzzz_0[j] + wp_x[j] * tg_yy_yyzzzzz_1[j];

                    tg_xyy_yzzzzzz_0[j] = pb_x * tg_yy_yzzzzzz_0[j] + wp_x[j] * tg_yy_yzzzzzz_1[j];

                    tg_xyy_zzzzzzz_0[j] = pb_x * tg_yy_zzzzzzz_0[j] + wp_x[j] * tg_yy_zzzzzzz_1[j];

                    tg_xyz_xxxxxxx_0[j] = pb_x * tg_yz_xxxxxxx_0[j] + wp_x[j] * tg_yz_xxxxxxx_1[j] + 3.5 * fl1_fxn * tg_yz_xxxxxx_1[j];

                    tg_xyz_xxxxxxy_0[j] = pb_x * tg_yz_xxxxxxy_0[j] + wp_x[j] * tg_yz_xxxxxxy_1[j] + 3.0 * fl1_fxn * tg_yz_xxxxxy_1[j];

                    tg_xyz_xxxxxxz_0[j] = pb_x * tg_yz_xxxxxxz_0[j] + wp_x[j] * tg_yz_xxxxxxz_1[j] + 3.0 * fl1_fxn * tg_yz_xxxxxz_1[j];

                    tg_xyz_xxxxxyy_0[j] = pb_x * tg_yz_xxxxxyy_0[j] + wp_x[j] * tg_yz_xxxxxyy_1[j] + 2.5 * fl1_fxn * tg_yz_xxxxyy_1[j];

                    tg_xyz_xxxxxyz_0[j] = pb_x * tg_yz_xxxxxyz_0[j] + wp_x[j] * tg_yz_xxxxxyz_1[j] + 2.5 * fl1_fxn * tg_yz_xxxxyz_1[j];

                    tg_xyz_xxxxxzz_0[j] = pb_x * tg_yz_xxxxxzz_0[j] + wp_x[j] * tg_yz_xxxxxzz_1[j] + 2.5 * fl1_fxn * tg_yz_xxxxzz_1[j];

                    tg_xyz_xxxxyyy_0[j] = pb_x * tg_yz_xxxxyyy_0[j] + wp_x[j] * tg_yz_xxxxyyy_1[j] + 2.0 * fl1_fxn * tg_yz_xxxyyy_1[j];

                    tg_xyz_xxxxyyz_0[j] = pb_x * tg_yz_xxxxyyz_0[j] + wp_x[j] * tg_yz_xxxxyyz_1[j] + 2.0 * fl1_fxn * tg_yz_xxxyyz_1[j];

                    tg_xyz_xxxxyzz_0[j] = pb_x * tg_yz_xxxxyzz_0[j] + wp_x[j] * tg_yz_xxxxyzz_1[j] + 2.0 * fl1_fxn * tg_yz_xxxyzz_1[j];

                    tg_xyz_xxxxzzz_0[j] = pb_x * tg_yz_xxxxzzz_0[j] + wp_x[j] * tg_yz_xxxxzzz_1[j] + 2.0 * fl1_fxn * tg_yz_xxxzzz_1[j];

                    tg_xyz_xxxyyyy_0[j] = pb_x * tg_yz_xxxyyyy_0[j] + wp_x[j] * tg_yz_xxxyyyy_1[j] + 1.5 * fl1_fxn * tg_yz_xxyyyy_1[j];

                    tg_xyz_xxxyyyz_0[j] = pb_x * tg_yz_xxxyyyz_0[j] + wp_x[j] * tg_yz_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_yz_xxyyyz_1[j];

                    tg_xyz_xxxyyzz_0[j] = pb_x * tg_yz_xxxyyzz_0[j] + wp_x[j] * tg_yz_xxxyyzz_1[j] + 1.5 * fl1_fxn * tg_yz_xxyyzz_1[j];

                    tg_xyz_xxxyzzz_0[j] = pb_x * tg_yz_xxxyzzz_0[j] + wp_x[j] * tg_yz_xxxyzzz_1[j] + 1.5 * fl1_fxn * tg_yz_xxyzzz_1[j];

                    tg_xyz_xxxzzzz_0[j] = pb_x * tg_yz_xxxzzzz_0[j] + wp_x[j] * tg_yz_xxxzzzz_1[j] + 1.5 * fl1_fxn * tg_yz_xxzzzz_1[j];

                    tg_xyz_xxyyyyy_0[j] = pb_x * tg_yz_xxyyyyy_0[j] + wp_x[j] * tg_yz_xxyyyyy_1[j] + fl1_fxn * tg_yz_xyyyyy_1[j];

                    tg_xyz_xxyyyyz_0[j] = pb_x * tg_yz_xxyyyyz_0[j] + wp_x[j] * tg_yz_xxyyyyz_1[j] + fl1_fxn * tg_yz_xyyyyz_1[j];

                    tg_xyz_xxyyyzz_0[j] = pb_x * tg_yz_xxyyyzz_0[j] + wp_x[j] * tg_yz_xxyyyzz_1[j] + fl1_fxn * tg_yz_xyyyzz_1[j];

                    tg_xyz_xxyyzzz_0[j] = pb_x * tg_yz_xxyyzzz_0[j] + wp_x[j] * tg_yz_xxyyzzz_1[j] + fl1_fxn * tg_yz_xyyzzz_1[j];

                    tg_xyz_xxyzzzz_0[j] = pb_x * tg_yz_xxyzzzz_0[j] + wp_x[j] * tg_yz_xxyzzzz_1[j] + fl1_fxn * tg_yz_xyzzzz_1[j];

                    tg_xyz_xxzzzzz_0[j] = pb_x * tg_yz_xxzzzzz_0[j] + wp_x[j] * tg_yz_xxzzzzz_1[j] + fl1_fxn * tg_yz_xzzzzz_1[j];

                    tg_xyz_xyyyyyy_0[j] = pb_x * tg_yz_xyyyyyy_0[j] + wp_x[j] * tg_yz_xyyyyyy_1[j] + 0.5 * fl1_fxn * tg_yz_yyyyyy_1[j];

                    tg_xyz_xyyyyyz_0[j] = pb_x * tg_yz_xyyyyyz_0[j] + wp_x[j] * tg_yz_xyyyyyz_1[j] + 0.5 * fl1_fxn * tg_yz_yyyyyz_1[j];

                    tg_xyz_xyyyyzz_0[j] = pb_x * tg_yz_xyyyyzz_0[j] + wp_x[j] * tg_yz_xyyyyzz_1[j] + 0.5 * fl1_fxn * tg_yz_yyyyzz_1[j];

                    tg_xyz_xyyyzzz_0[j] = pb_x * tg_yz_xyyyzzz_0[j] + wp_x[j] * tg_yz_xyyyzzz_1[j] + 0.5 * fl1_fxn * tg_yz_yyyzzz_1[j];

                    tg_xyz_xyyzzzz_0[j] = pb_x * tg_yz_xyyzzzz_0[j] + wp_x[j] * tg_yz_xyyzzzz_1[j] + 0.5 * fl1_fxn * tg_yz_yyzzzz_1[j];

                    tg_xyz_xyzzzzz_0[j] = pb_x * tg_yz_xyzzzzz_0[j] + wp_x[j] * tg_yz_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_yz_yzzzzz_1[j];

                    tg_xyz_xzzzzzz_0[j] = pb_x * tg_yz_xzzzzzz_0[j] + wp_x[j] * tg_yz_xzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yz_zzzzzz_1[j];

                    tg_xyz_yyyyyyy_0[j] = pb_x * tg_yz_yyyyyyy_0[j] + wp_x[j] * tg_yz_yyyyyyy_1[j];

                    tg_xyz_yyyyyyz_0[j] = pb_x * tg_yz_yyyyyyz_0[j] + wp_x[j] * tg_yz_yyyyyyz_1[j];

                    tg_xyz_yyyyyzz_0[j] = pb_x * tg_yz_yyyyyzz_0[j] + wp_x[j] * tg_yz_yyyyyzz_1[j];

                    tg_xyz_yyyyzzz_0[j] = pb_x * tg_yz_yyyyzzz_0[j] + wp_x[j] * tg_yz_yyyyzzz_1[j];

                    tg_xyz_yyyzzzz_0[j] = pb_x * tg_yz_yyyzzzz_0[j] + wp_x[j] * tg_yz_yyyzzzz_1[j];

                    tg_xyz_yyzzzzz_0[j] = pb_x * tg_yz_yyzzzzz_0[j] + wp_x[j] * tg_yz_yyzzzzz_1[j];

                    tg_xyz_yzzzzzz_0[j] = pb_x * tg_yz_yzzzzzz_0[j] + wp_x[j] * tg_yz_yzzzzzz_1[j];

                    tg_xyz_zzzzzzz_0[j] = pb_x * tg_yz_zzzzzzz_0[j] + wp_x[j] * tg_yz_zzzzzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSFSK_180_270(      CMemBlock2D<double>& primBuffer,
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

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {3, -1, -1, -1},
                                             {7, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_3_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_3_7_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_2_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_2_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_2_6_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {6, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_1_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_1_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {7, -1, -1, -1}, 
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

                auto tg_yy_xxxxxxx_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 108); 

                auto tg_yy_xxxxxxy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 109); 

                auto tg_yy_xxxxxxz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 110); 

                auto tg_yy_xxxxxyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 111); 

                auto tg_yy_xxxxxyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 112); 

                auto tg_yy_xxxxxzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 113); 

                auto tg_yy_xxxxyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 114); 

                auto tg_yy_xxxxyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 115); 

                auto tg_yy_xxxxyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 116); 

                auto tg_yy_xxxxzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 117); 

                auto tg_yy_xxxyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 118); 

                auto tg_yy_xxxyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 119); 

                auto tg_yy_xxxyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 120); 

                auto tg_yy_xxxyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 121); 

                auto tg_yy_xxxzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 122); 

                auto tg_yy_xxyyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 123); 

                auto tg_yy_xxyyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 124); 

                auto tg_yy_xxyyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 125); 

                auto tg_yy_xxyyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 126); 

                auto tg_yy_xxyzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 127); 

                auto tg_yy_xxzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 128); 

                auto tg_yy_xyyyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 129); 

                auto tg_yy_xyyyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 130); 

                auto tg_yy_xyyyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 131); 

                auto tg_yy_xyyyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 132); 

                auto tg_yy_xyyzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 133); 

                auto tg_yy_xyzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 134); 

                auto tg_yy_xzzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 135); 

                auto tg_yy_yyyyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 136); 

                auto tg_yy_yyyyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 137); 

                auto tg_yy_yyyyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 138); 

                auto tg_yy_yyyyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 139); 

                auto tg_yy_yyyzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 140); 

                auto tg_yy_yyzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 141); 

                auto tg_yy_yzzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 142); 

                auto tg_yy_zzzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 143); 

                auto tg_yz_xxxxxxx_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 144); 

                auto tg_yz_xxxxxxy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 145); 

                auto tg_yz_xxxxxxz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 146); 

                auto tg_yz_xxxxxyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 147); 

                auto tg_yz_xxxxxyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 148); 

                auto tg_yz_xxxxxzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 149); 

                auto tg_yz_xxxxyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 150); 

                auto tg_yz_xxxxyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 151); 

                auto tg_yz_xxxxyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 152); 

                auto tg_yz_xxxxzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 153); 

                auto tg_yz_xxxyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 154); 

                auto tg_yz_xxxyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 155); 

                auto tg_yz_xxxyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 156); 

                auto tg_yz_xxxyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 157); 

                auto tg_yz_xxxzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 158); 

                auto tg_yz_xxyyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 159); 

                auto tg_yz_xxyyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 160); 

                auto tg_yz_xxyyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 161); 

                auto tg_zz_xxxxxxx_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 180); 

                auto tg_zz_xxxxxxy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 181); 

                auto tg_zz_xxxxxxz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 182); 

                auto tg_zz_xxxxxyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 183); 

                auto tg_zz_xxxxxyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 184); 

                auto tg_zz_xxxxxzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 185); 

                auto tg_zz_xxxxyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 186); 

                auto tg_zz_xxxxyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 187); 

                auto tg_zz_xxxxyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 188); 

                auto tg_zz_xxxxzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 189); 

                auto tg_zz_xxxyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 190); 

                auto tg_zz_xxxyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 191); 

                auto tg_zz_xxxyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 192); 

                auto tg_zz_xxxyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 193); 

                auto tg_zz_xxxzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 194); 

                auto tg_zz_xxyyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 195); 

                auto tg_zz_xxyyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 196); 

                auto tg_zz_xxyyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 197); 

                auto tg_zz_xxyyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 198); 

                auto tg_zz_xxyzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 199); 

                auto tg_zz_xxzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 200); 

                auto tg_zz_xyyyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 201); 

                auto tg_zz_xyyyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 202); 

                auto tg_zz_xyyyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 203); 

                auto tg_zz_xyyyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 204); 

                auto tg_zz_xyyzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 205); 

                auto tg_zz_xyzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 206); 

                auto tg_zz_xzzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 207); 

                auto tg_zz_yyyyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 208); 

                auto tg_zz_yyyyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 209); 

                auto tg_zz_yyyyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 210); 

                auto tg_zz_yyyyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 211); 

                auto tg_zz_yyyzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 212); 

                auto tg_zz_yyzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 213); 

                auto tg_zz_yzzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 214); 

                auto tg_zz_zzzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 215); 

                auto tg_yy_xxxxxxx_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 108); 

                auto tg_yy_xxxxxxy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 109); 

                auto tg_yy_xxxxxxz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 110); 

                auto tg_yy_xxxxxyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 111); 

                auto tg_yy_xxxxxyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 112); 

                auto tg_yy_xxxxxzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 113); 

                auto tg_yy_xxxxyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 114); 

                auto tg_yy_xxxxyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 115); 

                auto tg_yy_xxxxyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 116); 

                auto tg_yy_xxxxzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 117); 

                auto tg_yy_xxxyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 118); 

                auto tg_yy_xxxyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 119); 

                auto tg_yy_xxxyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 120); 

                auto tg_yy_xxxyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 121); 

                auto tg_yy_xxxzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 122); 

                auto tg_yy_xxyyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 123); 

                auto tg_yy_xxyyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 124); 

                auto tg_yy_xxyyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 125); 

                auto tg_yy_xxyyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 126); 

                auto tg_yy_xxyzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 127); 

                auto tg_yy_xxzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 128); 

                auto tg_yy_xyyyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 129); 

                auto tg_yy_xyyyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 130); 

                auto tg_yy_xyyyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 131); 

                auto tg_yy_xyyyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 132); 

                auto tg_yy_xyyzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 133); 

                auto tg_yy_xyzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 134); 

                auto tg_yy_xzzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 135); 

                auto tg_yy_yyyyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 136); 

                auto tg_yy_yyyyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 137); 

                auto tg_yy_yyyyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 138); 

                auto tg_yy_yyyyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 139); 

                auto tg_yy_yyyzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 140); 

                auto tg_yy_yyzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 141); 

                auto tg_yy_yzzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 142); 

                auto tg_yy_zzzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 143); 

                auto tg_yz_xxxxxxx_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 144); 

                auto tg_yz_xxxxxxy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 145); 

                auto tg_yz_xxxxxxz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 146); 

                auto tg_yz_xxxxxyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 147); 

                auto tg_yz_xxxxxyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 148); 

                auto tg_yz_xxxxxzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 149); 

                auto tg_yz_xxxxyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 150); 

                auto tg_yz_xxxxyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 151); 

                auto tg_yz_xxxxyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 152); 

                auto tg_yz_xxxxzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 153); 

                auto tg_yz_xxxyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 154); 

                auto tg_yz_xxxyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 155); 

                auto tg_yz_xxxyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 156); 

                auto tg_yz_xxxyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 157); 

                auto tg_yz_xxxzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 158); 

                auto tg_yz_xxyyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 159); 

                auto tg_yz_xxyyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 160); 

                auto tg_yz_xxyyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 161); 

                auto tg_zz_xxxxxxx_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 180); 

                auto tg_zz_xxxxxxy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 181); 

                auto tg_zz_xxxxxxz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 182); 

                auto tg_zz_xxxxxyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 183); 

                auto tg_zz_xxxxxyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 184); 

                auto tg_zz_xxxxxzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 185); 

                auto tg_zz_xxxxyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 186); 

                auto tg_zz_xxxxyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 187); 

                auto tg_zz_xxxxyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 188); 

                auto tg_zz_xxxxzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 189); 

                auto tg_zz_xxxyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 190); 

                auto tg_zz_xxxyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 191); 

                auto tg_zz_xxxyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 192); 

                auto tg_zz_xxxyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 193); 

                auto tg_zz_xxxzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 194); 

                auto tg_zz_xxyyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 195); 

                auto tg_zz_xxyyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 196); 

                auto tg_zz_xxyyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 197); 

                auto tg_zz_xxyyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 198); 

                auto tg_zz_xxyzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 199); 

                auto tg_zz_xxzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 200); 

                auto tg_zz_xyyyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 201); 

                auto tg_zz_xyyyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 202); 

                auto tg_zz_xyyyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 203); 

                auto tg_zz_xyyyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 204); 

                auto tg_zz_xyyzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 205); 

                auto tg_zz_xyzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 206); 

                auto tg_zz_xzzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 207); 

                auto tg_zz_yyyyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 208); 

                auto tg_zz_yyyyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 209); 

                auto tg_zz_yyyyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 210); 

                auto tg_zz_yyyyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 211); 

                auto tg_zz_yyyzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 212); 

                auto tg_zz_yyzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 213); 

                auto tg_zz_yzzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 214); 

                auto tg_zz_zzzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 215); 

                auto tg_y_xxxxxxx_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 36); 

                auto tg_y_xxxxxxy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 37); 

                auto tg_y_xxxxxxz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 38); 

                auto tg_y_xxxxxyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 39); 

                auto tg_y_xxxxxyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 40); 

                auto tg_y_xxxxxzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 41); 

                auto tg_y_xxxxyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 42); 

                auto tg_y_xxxxyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 43); 

                auto tg_y_xxxxyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 44); 

                auto tg_y_xxxxzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 45); 

                auto tg_y_xxxyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 46); 

                auto tg_y_xxxyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 47); 

                auto tg_y_xxxyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 48); 

                auto tg_y_xxxyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 49); 

                auto tg_y_xxxzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 50); 

                auto tg_y_xxyyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 51); 

                auto tg_y_xxyyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 52); 

                auto tg_y_xxyyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 53); 

                auto tg_y_xxyyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 54); 

                auto tg_y_xxyzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 55); 

                auto tg_y_xxzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 56); 

                auto tg_y_xyyyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 57); 

                auto tg_y_xyyyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 58); 

                auto tg_y_xyyyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 59); 

                auto tg_y_xyyyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 60); 

                auto tg_y_xyyzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 61); 

                auto tg_y_xyzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 62); 

                auto tg_y_xzzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 63); 

                auto tg_y_yyyyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 64); 

                auto tg_y_yyyyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 65); 

                auto tg_y_yyyyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 66); 

                auto tg_y_yyyyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 67); 

                auto tg_y_yyyzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 68); 

                auto tg_y_yyzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 69); 

                auto tg_y_yzzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 70); 

                auto tg_y_zzzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 71); 

                auto tg_z_xxxxxxx_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 72); 

                auto tg_z_xxxxxxy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 73); 

                auto tg_z_xxxxxxz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 74); 

                auto tg_z_xxxxxyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 75); 

                auto tg_z_xxxxxyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 76); 

                auto tg_z_xxxxxzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 77); 

                auto tg_z_xxxxyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 78); 

                auto tg_z_xxxxyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 79); 

                auto tg_z_xxxxyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 80); 

                auto tg_z_xxxxzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 81); 

                auto tg_z_xxxyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 82); 

                auto tg_z_xxxyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 83); 

                auto tg_z_xxxyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 84); 

                auto tg_z_xxxyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 85); 

                auto tg_z_xxxzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 86); 

                auto tg_z_xxyyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 87); 

                auto tg_z_xxyyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 88); 

                auto tg_z_xxyyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 89); 

                auto tg_y_xxxxxxx_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 36); 

                auto tg_y_xxxxxxy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 37); 

                auto tg_y_xxxxxxz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 38); 

                auto tg_y_xxxxxyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 39); 

                auto tg_y_xxxxxyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 40); 

                auto tg_y_xxxxxzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 41); 

                auto tg_y_xxxxyyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 42); 

                auto tg_y_xxxxyyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 43); 

                auto tg_y_xxxxyzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 44); 

                auto tg_y_xxxxzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 45); 

                auto tg_y_xxxyyyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 46); 

                auto tg_y_xxxyyyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 47); 

                auto tg_y_xxxyyzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 48); 

                auto tg_y_xxxyzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 49); 

                auto tg_y_xxxzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 50); 

                auto tg_y_xxyyyyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 51); 

                auto tg_y_xxyyyyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 52); 

                auto tg_y_xxyyyzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 53); 

                auto tg_y_xxyyzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 54); 

                auto tg_y_xxyzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 55); 

                auto tg_y_xxzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 56); 

                auto tg_y_xyyyyyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 57); 

                auto tg_y_xyyyyyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 58); 

                auto tg_y_xyyyyzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 59); 

                auto tg_y_xyyyzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 60); 

                auto tg_y_xyyzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 61); 

                auto tg_y_xyzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 62); 

                auto tg_y_xzzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 63); 

                auto tg_y_yyyyyyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 64); 

                auto tg_y_yyyyyyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 65); 

                auto tg_y_yyyyyzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 66); 

                auto tg_y_yyyyzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 67); 

                auto tg_y_yyyzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 68); 

                auto tg_y_yyzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 69); 

                auto tg_y_yzzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 70); 

                auto tg_y_zzzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 71); 

                auto tg_z_xxxxxxx_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 72); 

                auto tg_z_xxxxxxy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 73); 

                auto tg_z_xxxxxxz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 74); 

                auto tg_z_xxxxxyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 75); 

                auto tg_z_xxxxxyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 76); 

                auto tg_z_xxxxxzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 77); 

                auto tg_z_xxxxyyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 78); 

                auto tg_z_xxxxyyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 79); 

                auto tg_z_xxxxyzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 80); 

                auto tg_z_xxxxzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 81); 

                auto tg_z_xxxyyyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 82); 

                auto tg_z_xxxyyyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 83); 

                auto tg_z_xxxyyzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 84); 

                auto tg_z_xxxyzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 85); 

                auto tg_z_xxxzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 86); 

                auto tg_z_xxyyyyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 87); 

                auto tg_z_xxyyyyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 88); 

                auto tg_z_xxyyyzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 89); 

                auto tg_yy_xxxxxx_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 84); 

                auto tg_yy_xxxxxy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 85); 

                auto tg_yy_xxxxxz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 86); 

                auto tg_yy_xxxxyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 87); 

                auto tg_yy_xxxxyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 88); 

                auto tg_yy_xxxxzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 89); 

                auto tg_yy_xxxyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 90); 

                auto tg_yy_xxxyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 91); 

                auto tg_yy_xxxyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 92); 

                auto tg_yy_xxxzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 93); 

                auto tg_yy_xxyyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 94); 

                auto tg_yy_xxyyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 95); 

                auto tg_yy_xxyyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 96); 

                auto tg_yy_xxyzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 97); 

                auto tg_yy_xxzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 98); 

                auto tg_yy_xyyyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 99); 

                auto tg_yy_xyyyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 100); 

                auto tg_yy_xyyyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 101); 

                auto tg_yy_xyyzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 102); 

                auto tg_yy_xyzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 103); 

                auto tg_yy_xzzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 104); 

                auto tg_yy_yyyyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 105); 

                auto tg_yy_yyyyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 106); 

                auto tg_yy_yyyyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 107); 

                auto tg_yy_yyyzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 108); 

                auto tg_yy_yyzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 109); 

                auto tg_yy_yzzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 110); 

                auto tg_yy_zzzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 111); 

                auto tg_yz_xxxxxx_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 112); 

                auto tg_yz_xxxxxy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 113); 

                auto tg_yz_xxxxxz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 114); 

                auto tg_yz_xxxxyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 115); 

                auto tg_yz_xxxxyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 116); 

                auto tg_yz_xxxxzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 117); 

                auto tg_yz_xxxyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 118); 

                auto tg_yz_xxxyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 119); 

                auto tg_yz_xxxyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 120); 

                auto tg_yz_xxxzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 121); 

                auto tg_yz_xxyyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 122); 

                auto tg_yz_xxyyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 123); 

                auto tg_yz_xxyyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 124); 

                auto tg_zz_xxxxxx_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 140); 

                auto tg_zz_xxxxxy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 141); 

                auto tg_zz_xxxxxz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 142); 

                auto tg_zz_xxxxyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 143); 

                auto tg_zz_xxxxyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 144); 

                auto tg_zz_xxxxzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 145); 

                auto tg_zz_xxxyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 146); 

                auto tg_zz_xxxyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 147); 

                auto tg_zz_xxxyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 148); 

                auto tg_zz_xxxzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 149); 

                auto tg_zz_xxyyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 150); 

                auto tg_zz_xxyyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 151); 

                auto tg_zz_xxyyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 152); 

                auto tg_zz_xxyzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 153); 

                auto tg_zz_xxzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 154); 

                auto tg_zz_xyyyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 155); 

                auto tg_zz_xyyyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 156); 

                auto tg_zz_xyyyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 157); 

                auto tg_zz_xyyzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 158); 

                auto tg_zz_xyzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 159); 

                auto tg_zz_xzzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 160); 

                auto tg_zz_yyyyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 161); 

                auto tg_zz_yyyyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 162); 

                auto tg_zz_yyyyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 163); 

                auto tg_zz_yyyzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 164); 

                auto tg_zz_yyzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 165); 

                auto tg_zz_yzzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 166); 

                auto tg_zz_zzzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 167); 

                // set up pointers to integrals

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

                // Batch of Integrals (180,270)

                #pragma omp simd aligned(fxn, fza, tg_xzz_xxxxxxx_0, tg_xzz_xxxxxxy_0, tg_xzz_xxxxxxz_0, \
                                         tg_xzz_xxxxxyy_0, tg_xzz_xxxxxyz_0, tg_xzz_xxxxxzz_0, tg_xzz_xxxxyyy_0, \
                                         tg_xzz_xxxxyyz_0, tg_xzz_xxxxyzz_0, tg_xzz_xxxxzzz_0, tg_xzz_xxxyyyy_0, \
                                         tg_xzz_xxxyyyz_0, tg_xzz_xxxyyzz_0, tg_xzz_xxxyzzz_0, tg_xzz_xxxzzzz_0, \
                                         tg_xzz_xxyyyyy_0, tg_xzz_xxyyyyz_0, tg_xzz_xxyyyzz_0, tg_xzz_xxyyzzz_0, \
                                         tg_xzz_xxyzzzz_0, tg_xzz_xxzzzzz_0, tg_xzz_xyyyyyy_0, tg_xzz_xyyyyyz_0, \
                                         tg_xzz_xyyyyzz_0, tg_xzz_xyyyzzz_0, tg_xzz_xyyzzzz_0, tg_xzz_xyzzzzz_0, \
                                         tg_xzz_xzzzzzz_0, tg_xzz_yyyyyyy_0, tg_xzz_yyyyyyz_0, tg_xzz_yyyyyzz_0, \
                                         tg_xzz_yyyyzzz_0, tg_xzz_yyyzzzz_0, tg_xzz_yyzzzzz_0, tg_xzz_yzzzzzz_0, \
                                         tg_xzz_zzzzzzz_0, tg_y_xxxxxxx_0, tg_y_xxxxxxx_1, tg_y_xxxxxxy_0, tg_y_xxxxxxy_1, \
                                         tg_y_xxxxxxz_0, tg_y_xxxxxxz_1, tg_y_xxxxxyy_0, tg_y_xxxxxyy_1, tg_y_xxxxxyz_0, \
                                         tg_y_xxxxxyz_1, tg_y_xxxxxzz_0, tg_y_xxxxxzz_1, tg_y_xxxxyyy_0, tg_y_xxxxyyy_1, \
                                         tg_y_xxxxyyz_0, tg_y_xxxxyyz_1, tg_y_xxxxyzz_0, tg_y_xxxxyzz_1, tg_y_xxxxzzz_0, \
                                         tg_y_xxxxzzz_1, tg_y_xxxyyyy_0, tg_y_xxxyyyy_1, tg_y_xxxyyyz_0, tg_y_xxxyyyz_1, \
                                         tg_y_xxxyyzz_0, tg_y_xxxyyzz_1, tg_y_xxxyzzz_0, tg_y_xxxyzzz_1, tg_y_xxxzzzz_0, \
                                         tg_y_xxxzzzz_1, tg_y_xxyyyyy_0, tg_y_xxyyyyy_1, tg_y_xxyyyyz_0, tg_y_xxyyyyz_1, \
                                         tg_y_xxyyyzz_0, tg_y_xxyyyzz_1, tg_y_xxyyzzz_0, tg_y_xxyyzzz_1, tg_y_xxyzzzz_0, \
                                         tg_y_xxyzzzz_1, tg_y_xxzzzzz_0, tg_y_xxzzzzz_1, tg_y_xyyyyyy_0, tg_y_xyyyyyy_1, \
                                         tg_y_xyyyyyz_0, tg_y_xyyyyyz_1, tg_y_xyyyyzz_0, tg_y_xyyyyzz_1, tg_y_xyyyzzz_0, \
                                         tg_y_xyyyzzz_1, tg_y_xyyzzzz_0, tg_y_xyyzzzz_1, tg_y_xyzzzzz_0, tg_y_xyzzzzz_1, \
                                         tg_y_xzzzzzz_0, tg_y_xzzzzzz_1, tg_y_yyyyyyy_0, tg_y_yyyyyyy_1, tg_y_yyyyyyz_0, \
                                         tg_y_yyyyyyz_1, tg_y_yyyyyzz_0, tg_y_yyyyyzz_1, tg_y_yyyyzzz_0, tg_y_yyyyzzz_1, \
                                         tg_y_yyyzzzz_0, tg_y_yyyzzzz_1, tg_y_yyzzzzz_0, tg_y_yyzzzzz_1, tg_y_yzzzzzz_0, \
                                         tg_y_yzzzzzz_1, tg_y_zzzzzzz_0, tg_y_zzzzzzz_1, tg_yy_xxxxxx_1, tg_yy_xxxxxxx_0, \
                                         tg_yy_xxxxxxx_1, tg_yy_xxxxxxy_0, tg_yy_xxxxxxy_1, tg_yy_xxxxxxz_0, tg_yy_xxxxxxz_1, \
                                         tg_yy_xxxxxy_1, tg_yy_xxxxxyy_0, tg_yy_xxxxxyy_1, tg_yy_xxxxxyz_0, tg_yy_xxxxxyz_1, \
                                         tg_yy_xxxxxz_1, tg_yy_xxxxxzz_0, tg_yy_xxxxxzz_1, tg_yy_xxxxyy_1, tg_yy_xxxxyyy_0, \
                                         tg_yy_xxxxyyy_1, tg_yy_xxxxyyz_0, tg_yy_xxxxyyz_1, tg_yy_xxxxyz_1, tg_yy_xxxxyzz_0, \
                                         tg_yy_xxxxyzz_1, tg_yy_xxxxzz_1, tg_yy_xxxxzzz_0, tg_yy_xxxxzzz_1, tg_yy_xxxyyy_1, \
                                         tg_yy_xxxyyyy_0, tg_yy_xxxyyyy_1, tg_yy_xxxyyyz_0, tg_yy_xxxyyyz_1, tg_yy_xxxyyz_1, \
                                         tg_yy_xxxyyzz_0, tg_yy_xxxyyzz_1, tg_yy_xxxyzz_1, tg_yy_xxxyzzz_0, tg_yy_xxxyzzz_1, \
                                         tg_yy_xxxzzz_1, tg_yy_xxxzzzz_0, tg_yy_xxxzzzz_1, tg_yy_xxyyyy_1, tg_yy_xxyyyyy_0, \
                                         tg_yy_xxyyyyy_1, tg_yy_xxyyyyz_0, tg_yy_xxyyyyz_1, tg_yy_xxyyyz_1, tg_yy_xxyyyzz_0, \
                                         tg_yy_xxyyyzz_1, tg_yy_xxyyzz_1, tg_yy_xxyyzzz_0, tg_yy_xxyyzzz_1, tg_yy_xxyzzz_1, \
                                         tg_yy_xxyzzzz_0, tg_yy_xxyzzzz_1, tg_yy_xxzzzz_1, tg_yy_xxzzzzz_0, tg_yy_xxzzzzz_1, \
                                         tg_yy_xyyyyy_1, tg_yy_xyyyyyy_0, tg_yy_xyyyyyy_1, tg_yy_xyyyyyz_0, tg_yy_xyyyyyz_1, \
                                         tg_yy_xyyyyz_1, tg_yy_xyyyyzz_0, tg_yy_xyyyyzz_1, tg_yy_xyyyzz_1, tg_yy_xyyyzzz_0, \
                                         tg_yy_xyyyzzz_1, tg_yy_xyyzzz_1, tg_yy_xyyzzzz_0, tg_yy_xyyzzzz_1, tg_yy_xyzzzz_1, \
                                         tg_yy_xyzzzzz_0, tg_yy_xyzzzzz_1, tg_yy_xzzzzz_1, tg_yy_xzzzzzz_0, tg_yy_xzzzzzz_1, \
                                         tg_yy_yyyyyy_1, tg_yy_yyyyyyy_0, tg_yy_yyyyyyy_1, tg_yy_yyyyyyz_0, tg_yy_yyyyyyz_1, \
                                         tg_yy_yyyyyz_1, tg_yy_yyyyyzz_0, tg_yy_yyyyyzz_1, tg_yy_yyyyzz_1, tg_yy_yyyyzzz_0, \
                                         tg_yy_yyyyzzz_1, tg_yy_yyyzzz_1, tg_yy_yyyzzzz_0, tg_yy_yyyzzzz_1, tg_yy_yyzzzz_1, \
                                         tg_yy_yyzzzzz_0, tg_yy_yyzzzzz_1, tg_yy_yzzzzz_1, tg_yy_yzzzzzz_0, tg_yy_yzzzzzz_1, \
                                         tg_yy_zzzzzz_1, tg_yy_zzzzzzz_0, tg_yy_zzzzzzz_1, tg_yyy_xxxxxxx_0, \
                                         tg_yyy_xxxxxxy_0, tg_yyy_xxxxxxz_0, tg_yyy_xxxxxyy_0, tg_yyy_xxxxxyz_0, \
                                         tg_yyy_xxxxxzz_0, tg_yyy_xxxxyyy_0, tg_yyy_xxxxyyz_0, tg_yyy_xxxxyzz_0, \
                                         tg_yyy_xxxxzzz_0, tg_yyy_xxxyyyy_0, tg_yyy_xxxyyyz_0, tg_yyy_xxxyyzz_0, \
                                         tg_yyy_xxxyzzz_0, tg_yyy_xxxzzzz_0, tg_yyy_xxyyyyy_0, tg_yyy_xxyyyyz_0, \
                                         tg_yyy_xxyyyzz_0, tg_yyy_xxyyzzz_0, tg_yyy_xxyzzzz_0, tg_yyy_xxzzzzz_0, \
                                         tg_yyy_xyyyyyy_0, tg_yyy_xyyyyyz_0, tg_yyy_xyyyyzz_0, tg_yyy_xyyyzzz_0, \
                                         tg_yyy_xyyzzzz_0, tg_yyy_xyzzzzz_0, tg_yyy_xzzzzzz_0, tg_yyy_yyyyyyy_0, \
                                         tg_yyy_yyyyyyz_0, tg_yyy_yyyyyzz_0, tg_yyy_yyyyzzz_0, tg_yyy_yyyzzzz_0, \
                                         tg_yyy_yyzzzzz_0, tg_yyy_yzzzzzz_0, tg_yyy_zzzzzzz_0, tg_yyz_xxxxxxx_0, \
                                         tg_yyz_xxxxxxy_0, tg_yyz_xxxxxxz_0, tg_yyz_xxxxxyy_0, tg_yyz_xxxxxyz_0, \
                                         tg_yyz_xxxxxzz_0, tg_yyz_xxxxyyy_0, tg_yyz_xxxxyyz_0, tg_yyz_xxxxyzz_0, \
                                         tg_yyz_xxxxzzz_0, tg_yyz_xxxyyyy_0, tg_yyz_xxxyyyz_0, tg_yyz_xxxyyzz_0, \
                                         tg_yyz_xxxyzzz_0, tg_yyz_xxxzzzz_0, tg_yyz_xxyyyyy_0, tg_yyz_xxyyyyz_0, \
                                         tg_yyz_xxyyyzz_0, tg_yz_xxxxxx_1, tg_yz_xxxxxxx_0, tg_yz_xxxxxxx_1, tg_yz_xxxxxxy_0, \
                                         tg_yz_xxxxxxy_1, tg_yz_xxxxxxz_0, tg_yz_xxxxxxz_1, tg_yz_xxxxxy_1, tg_yz_xxxxxyy_0, \
                                         tg_yz_xxxxxyy_1, tg_yz_xxxxxyz_0, tg_yz_xxxxxyz_1, tg_yz_xxxxxz_1, tg_yz_xxxxxzz_0, \
                                         tg_yz_xxxxxzz_1, tg_yz_xxxxyy_1, tg_yz_xxxxyyy_0, tg_yz_xxxxyyy_1, tg_yz_xxxxyyz_0, \
                                         tg_yz_xxxxyyz_1, tg_yz_xxxxyz_1, tg_yz_xxxxyzz_0, tg_yz_xxxxyzz_1, tg_yz_xxxxzz_1, \
                                         tg_yz_xxxxzzz_0, tg_yz_xxxxzzz_1, tg_yz_xxxyyy_1, tg_yz_xxxyyyy_0, tg_yz_xxxyyyy_1, \
                                         tg_yz_xxxyyyz_0, tg_yz_xxxyyyz_1, tg_yz_xxxyyz_1, tg_yz_xxxyyzz_0, tg_yz_xxxyyzz_1, \
                                         tg_yz_xxxyzz_1, tg_yz_xxxyzzz_0, tg_yz_xxxyzzz_1, tg_yz_xxxzzz_1, tg_yz_xxxzzzz_0, \
                                         tg_yz_xxxzzzz_1, tg_yz_xxyyyy_1, tg_yz_xxyyyyy_0, tg_yz_xxyyyyy_1, tg_yz_xxyyyyz_0, \
                                         tg_yz_xxyyyyz_1, tg_yz_xxyyyz_1, tg_yz_xxyyyzz_0, tg_yz_xxyyyzz_1, tg_yz_xxyyzz_1, \
                                         tg_z_xxxxxxx_0, tg_z_xxxxxxx_1, tg_z_xxxxxxy_0, tg_z_xxxxxxy_1, tg_z_xxxxxxz_0, \
                                         tg_z_xxxxxxz_1, tg_z_xxxxxyy_0, tg_z_xxxxxyy_1, tg_z_xxxxxyz_0, tg_z_xxxxxyz_1, \
                                         tg_z_xxxxxzz_0, tg_z_xxxxxzz_1, tg_z_xxxxyyy_0, tg_z_xxxxyyy_1, tg_z_xxxxyyz_0, \
                                         tg_z_xxxxyyz_1, tg_z_xxxxyzz_0, tg_z_xxxxyzz_1, tg_z_xxxxzzz_0, tg_z_xxxxzzz_1, \
                                         tg_z_xxxyyyy_0, tg_z_xxxyyyy_1, tg_z_xxxyyyz_0, tg_z_xxxyyyz_1, tg_z_xxxyyzz_0, \
                                         tg_z_xxxyyzz_1, tg_z_xxxyzzz_0, tg_z_xxxyzzz_1, tg_z_xxxzzzz_0, tg_z_xxxzzzz_1, \
                                         tg_z_xxyyyyy_0, tg_z_xxyyyyy_1, tg_z_xxyyyyz_0, tg_z_xxyyyyz_1, tg_z_xxyyyzz_0, \
                                         tg_z_xxyyyzz_1, tg_zz_xxxxxx_1, tg_zz_xxxxxxx_0, tg_zz_xxxxxxx_1, tg_zz_xxxxxxy_0, \
                                         tg_zz_xxxxxxy_1, tg_zz_xxxxxxz_0, tg_zz_xxxxxxz_1, tg_zz_xxxxxy_1, tg_zz_xxxxxyy_0, \
                                         tg_zz_xxxxxyy_1, tg_zz_xxxxxyz_0, tg_zz_xxxxxyz_1, tg_zz_xxxxxz_1, tg_zz_xxxxxzz_0, \
                                         tg_zz_xxxxxzz_1, tg_zz_xxxxyy_1, tg_zz_xxxxyyy_0, tg_zz_xxxxyyy_1, tg_zz_xxxxyyz_0, \
                                         tg_zz_xxxxyyz_1, tg_zz_xxxxyz_1, tg_zz_xxxxyzz_0, tg_zz_xxxxyzz_1, tg_zz_xxxxzz_1, \
                                         tg_zz_xxxxzzz_0, tg_zz_xxxxzzz_1, tg_zz_xxxyyy_1, tg_zz_xxxyyyy_0, tg_zz_xxxyyyy_1, \
                                         tg_zz_xxxyyyz_0, tg_zz_xxxyyyz_1, tg_zz_xxxyyz_1, tg_zz_xxxyyzz_0, tg_zz_xxxyyzz_1, \
                                         tg_zz_xxxyzz_1, tg_zz_xxxyzzz_0, tg_zz_xxxyzzz_1, tg_zz_xxxzzz_1, tg_zz_xxxzzzz_0, \
                                         tg_zz_xxxzzzz_1, tg_zz_xxyyyy_1, tg_zz_xxyyyyy_0, tg_zz_xxyyyyy_1, tg_zz_xxyyyyz_0, \
                                         tg_zz_xxyyyyz_1, tg_zz_xxyyyz_1, tg_zz_xxyyyzz_0, tg_zz_xxyyyzz_1, tg_zz_xxyyzz_1, \
                                         tg_zz_xxyyzzz_0, tg_zz_xxyyzzz_1, tg_zz_xxyzzz_1, tg_zz_xxyzzzz_0, tg_zz_xxyzzzz_1, \
                                         tg_zz_xxzzzz_1, tg_zz_xxzzzzz_0, tg_zz_xxzzzzz_1, tg_zz_xyyyyy_1, tg_zz_xyyyyyy_0, \
                                         tg_zz_xyyyyyy_1, tg_zz_xyyyyyz_0, tg_zz_xyyyyyz_1, tg_zz_xyyyyz_1, tg_zz_xyyyyzz_0, \
                                         tg_zz_xyyyyzz_1, tg_zz_xyyyzz_1, tg_zz_xyyyzzz_0, tg_zz_xyyyzzz_1, tg_zz_xyyzzz_1, \
                                         tg_zz_xyyzzzz_0, tg_zz_xyyzzzz_1, tg_zz_xyzzzz_1, tg_zz_xyzzzzz_0, tg_zz_xyzzzzz_1, \
                                         tg_zz_xzzzzz_1, tg_zz_xzzzzzz_0, tg_zz_xzzzzzz_1, tg_zz_yyyyyy_1, tg_zz_yyyyyyy_0, \
                                         tg_zz_yyyyyyy_1, tg_zz_yyyyyyz_0, tg_zz_yyyyyyz_1, tg_zz_yyyyyz_1, tg_zz_yyyyyzz_0, \
                                         tg_zz_yyyyyzz_1, tg_zz_yyyyzz_1, tg_zz_yyyyzzz_0, tg_zz_yyyyzzz_1, tg_zz_yyyzzz_1, \
                                         tg_zz_yyyzzzz_0, tg_zz_yyyzzzz_1, tg_zz_yyzzzz_1, tg_zz_yyzzzzz_0, tg_zz_yyzzzzz_1, \
                                         tg_zz_yzzzzz_1, tg_zz_yzzzzzz_0, tg_zz_yzzzzzz_1, tg_zz_zzzzzz_1, tg_zz_zzzzzzz_0, \
                                         tg_zz_zzzzzzz_1, wp_x, wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xzz_xxxxxxx_0[j] = pb_x * tg_zz_xxxxxxx_0[j] + wp_x[j] * tg_zz_xxxxxxx_1[j] + 3.5 * fl1_fxn * tg_zz_xxxxxx_1[j];

                    tg_xzz_xxxxxxy_0[j] = pb_x * tg_zz_xxxxxxy_0[j] + wp_x[j] * tg_zz_xxxxxxy_1[j] + 3.0 * fl1_fxn * tg_zz_xxxxxy_1[j];

                    tg_xzz_xxxxxxz_0[j] = pb_x * tg_zz_xxxxxxz_0[j] + wp_x[j] * tg_zz_xxxxxxz_1[j] + 3.0 * fl1_fxn * tg_zz_xxxxxz_1[j];

                    tg_xzz_xxxxxyy_0[j] = pb_x * tg_zz_xxxxxyy_0[j] + wp_x[j] * tg_zz_xxxxxyy_1[j] + 2.5 * fl1_fxn * tg_zz_xxxxyy_1[j];

                    tg_xzz_xxxxxyz_0[j] = pb_x * tg_zz_xxxxxyz_0[j] + wp_x[j] * tg_zz_xxxxxyz_1[j] + 2.5 * fl1_fxn * tg_zz_xxxxyz_1[j];

                    tg_xzz_xxxxxzz_0[j] = pb_x * tg_zz_xxxxxzz_0[j] + wp_x[j] * tg_zz_xxxxxzz_1[j] + 2.5 * fl1_fxn * tg_zz_xxxxzz_1[j];

                    tg_xzz_xxxxyyy_0[j] = pb_x * tg_zz_xxxxyyy_0[j] + wp_x[j] * tg_zz_xxxxyyy_1[j] + 2.0 * fl1_fxn * tg_zz_xxxyyy_1[j];

                    tg_xzz_xxxxyyz_0[j] = pb_x * tg_zz_xxxxyyz_0[j] + wp_x[j] * tg_zz_xxxxyyz_1[j] + 2.0 * fl1_fxn * tg_zz_xxxyyz_1[j];

                    tg_xzz_xxxxyzz_0[j] = pb_x * tg_zz_xxxxyzz_0[j] + wp_x[j] * tg_zz_xxxxyzz_1[j] + 2.0 * fl1_fxn * tg_zz_xxxyzz_1[j];

                    tg_xzz_xxxxzzz_0[j] = pb_x * tg_zz_xxxxzzz_0[j] + wp_x[j] * tg_zz_xxxxzzz_1[j] + 2.0 * fl1_fxn * tg_zz_xxxzzz_1[j];

                    tg_xzz_xxxyyyy_0[j] = pb_x * tg_zz_xxxyyyy_0[j] + wp_x[j] * tg_zz_xxxyyyy_1[j] + 1.5 * fl1_fxn * tg_zz_xxyyyy_1[j];

                    tg_xzz_xxxyyyz_0[j] = pb_x * tg_zz_xxxyyyz_0[j] + wp_x[j] * tg_zz_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_zz_xxyyyz_1[j];

                    tg_xzz_xxxyyzz_0[j] = pb_x * tg_zz_xxxyyzz_0[j] + wp_x[j] * tg_zz_xxxyyzz_1[j] + 1.5 * fl1_fxn * tg_zz_xxyyzz_1[j];

                    tg_xzz_xxxyzzz_0[j] = pb_x * tg_zz_xxxyzzz_0[j] + wp_x[j] * tg_zz_xxxyzzz_1[j] + 1.5 * fl1_fxn * tg_zz_xxyzzz_1[j];

                    tg_xzz_xxxzzzz_0[j] = pb_x * tg_zz_xxxzzzz_0[j] + wp_x[j] * tg_zz_xxxzzzz_1[j] + 1.5 * fl1_fxn * tg_zz_xxzzzz_1[j];

                    tg_xzz_xxyyyyy_0[j] = pb_x * tg_zz_xxyyyyy_0[j] + wp_x[j] * tg_zz_xxyyyyy_1[j] + fl1_fxn * tg_zz_xyyyyy_1[j];

                    tg_xzz_xxyyyyz_0[j] = pb_x * tg_zz_xxyyyyz_0[j] + wp_x[j] * tg_zz_xxyyyyz_1[j] + fl1_fxn * tg_zz_xyyyyz_1[j];

                    tg_xzz_xxyyyzz_0[j] = pb_x * tg_zz_xxyyyzz_0[j] + wp_x[j] * tg_zz_xxyyyzz_1[j] + fl1_fxn * tg_zz_xyyyzz_1[j];

                    tg_xzz_xxyyzzz_0[j] = pb_x * tg_zz_xxyyzzz_0[j] + wp_x[j] * tg_zz_xxyyzzz_1[j] + fl1_fxn * tg_zz_xyyzzz_1[j];

                    tg_xzz_xxyzzzz_0[j] = pb_x * tg_zz_xxyzzzz_0[j] + wp_x[j] * tg_zz_xxyzzzz_1[j] + fl1_fxn * tg_zz_xyzzzz_1[j];

                    tg_xzz_xxzzzzz_0[j] = pb_x * tg_zz_xxzzzzz_0[j] + wp_x[j] * tg_zz_xxzzzzz_1[j] + fl1_fxn * tg_zz_xzzzzz_1[j];

                    tg_xzz_xyyyyyy_0[j] = pb_x * tg_zz_xyyyyyy_0[j] + wp_x[j] * tg_zz_xyyyyyy_1[j] + 0.5 * fl1_fxn * tg_zz_yyyyyy_1[j];

                    tg_xzz_xyyyyyz_0[j] = pb_x * tg_zz_xyyyyyz_0[j] + wp_x[j] * tg_zz_xyyyyyz_1[j] + 0.5 * fl1_fxn * tg_zz_yyyyyz_1[j];

                    tg_xzz_xyyyyzz_0[j] = pb_x * tg_zz_xyyyyzz_0[j] + wp_x[j] * tg_zz_xyyyyzz_1[j] + 0.5 * fl1_fxn * tg_zz_yyyyzz_1[j];

                    tg_xzz_xyyyzzz_0[j] = pb_x * tg_zz_xyyyzzz_0[j] + wp_x[j] * tg_zz_xyyyzzz_1[j] + 0.5 * fl1_fxn * tg_zz_yyyzzz_1[j];

                    tg_xzz_xyyzzzz_0[j] = pb_x * tg_zz_xyyzzzz_0[j] + wp_x[j] * tg_zz_xyyzzzz_1[j] + 0.5 * fl1_fxn * tg_zz_yyzzzz_1[j];

                    tg_xzz_xyzzzzz_0[j] = pb_x * tg_zz_xyzzzzz_0[j] + wp_x[j] * tg_zz_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_zz_yzzzzz_1[j];

                    tg_xzz_xzzzzzz_0[j] = pb_x * tg_zz_xzzzzzz_0[j] + wp_x[j] * tg_zz_xzzzzzz_1[j] + 0.5 * fl1_fxn * tg_zz_zzzzzz_1[j];

                    tg_xzz_yyyyyyy_0[j] = pb_x * tg_zz_yyyyyyy_0[j] + wp_x[j] * tg_zz_yyyyyyy_1[j];

                    tg_xzz_yyyyyyz_0[j] = pb_x * tg_zz_yyyyyyz_0[j] + wp_x[j] * tg_zz_yyyyyyz_1[j];

                    tg_xzz_yyyyyzz_0[j] = pb_x * tg_zz_yyyyyzz_0[j] + wp_x[j] * tg_zz_yyyyyzz_1[j];

                    tg_xzz_yyyyzzz_0[j] = pb_x * tg_zz_yyyyzzz_0[j] + wp_x[j] * tg_zz_yyyyzzz_1[j];

                    tg_xzz_yyyzzzz_0[j] = pb_x * tg_zz_yyyzzzz_0[j] + wp_x[j] * tg_zz_yyyzzzz_1[j];

                    tg_xzz_yyzzzzz_0[j] = pb_x * tg_zz_yyzzzzz_0[j] + wp_x[j] * tg_zz_yyzzzzz_1[j];

                    tg_xzz_yzzzzzz_0[j] = pb_x * tg_zz_yzzzzzz_0[j] + wp_x[j] * tg_zz_yzzzzzz_1[j];

                    tg_xzz_zzzzzzz_0[j] = pb_x * tg_zz_zzzzzzz_0[j] + wp_x[j] * tg_zz_zzzzzzz_1[j];

                    tg_yyy_xxxxxxx_0[j] = pb_y * tg_yy_xxxxxxx_0[j] + wp_y[j] * tg_yy_xxxxxxx_1[j] + fl1_fx * tg_y_xxxxxxx_0[j] - fl1_fx * fl1_fza * tg_y_xxxxxxx_1[j];

                    tg_yyy_xxxxxxy_0[j] = pb_y * tg_yy_xxxxxxy_0[j] + wp_y[j] * tg_yy_xxxxxxy_1[j] + fl1_fx * tg_y_xxxxxxy_0[j] - fl1_fx * fl1_fza * tg_y_xxxxxxy_1[j] + 0.5 * fl1_fxn * tg_yy_xxxxxx_1[j];

                    tg_yyy_xxxxxxz_0[j] = pb_y * tg_yy_xxxxxxz_0[j] + wp_y[j] * tg_yy_xxxxxxz_1[j] + fl1_fx * tg_y_xxxxxxz_0[j] - fl1_fx * fl1_fza * tg_y_xxxxxxz_1[j];

                    tg_yyy_xxxxxyy_0[j] = pb_y * tg_yy_xxxxxyy_0[j] + wp_y[j] * tg_yy_xxxxxyy_1[j] + fl1_fx * tg_y_xxxxxyy_0[j] - fl1_fx * fl1_fza * tg_y_xxxxxyy_1[j] + fl1_fxn * tg_yy_xxxxxy_1[j];

                    tg_yyy_xxxxxyz_0[j] = pb_y * tg_yy_xxxxxyz_0[j] + wp_y[j] * tg_yy_xxxxxyz_1[j] + fl1_fx * tg_y_xxxxxyz_0[j] - fl1_fx * fl1_fza * tg_y_xxxxxyz_1[j] + 0.5 * fl1_fxn * tg_yy_xxxxxz_1[j];

                    tg_yyy_xxxxxzz_0[j] = pb_y * tg_yy_xxxxxzz_0[j] + wp_y[j] * tg_yy_xxxxxzz_1[j] + fl1_fx * tg_y_xxxxxzz_0[j] - fl1_fx * fl1_fza * tg_y_xxxxxzz_1[j];

                    tg_yyy_xxxxyyy_0[j] = pb_y * tg_yy_xxxxyyy_0[j] + wp_y[j] * tg_yy_xxxxyyy_1[j] + fl1_fx * tg_y_xxxxyyy_0[j] - fl1_fx * fl1_fza * tg_y_xxxxyyy_1[j] + 1.5 * fl1_fxn * tg_yy_xxxxyy_1[j];

                    tg_yyy_xxxxyyz_0[j] = pb_y * tg_yy_xxxxyyz_0[j] + wp_y[j] * tg_yy_xxxxyyz_1[j] + fl1_fx * tg_y_xxxxyyz_0[j] - fl1_fx * fl1_fza * tg_y_xxxxyyz_1[j] + fl1_fxn * tg_yy_xxxxyz_1[j];

                    tg_yyy_xxxxyzz_0[j] = pb_y * tg_yy_xxxxyzz_0[j] + wp_y[j] * tg_yy_xxxxyzz_1[j] + fl1_fx * tg_y_xxxxyzz_0[j] - fl1_fx * fl1_fza * tg_y_xxxxyzz_1[j] + 0.5 * fl1_fxn * tg_yy_xxxxzz_1[j];

                    tg_yyy_xxxxzzz_0[j] = pb_y * tg_yy_xxxxzzz_0[j] + wp_y[j] * tg_yy_xxxxzzz_1[j] + fl1_fx * tg_y_xxxxzzz_0[j] - fl1_fx * fl1_fza * tg_y_xxxxzzz_1[j];

                    tg_yyy_xxxyyyy_0[j] = pb_y * tg_yy_xxxyyyy_0[j] + wp_y[j] * tg_yy_xxxyyyy_1[j] + fl1_fx * tg_y_xxxyyyy_0[j] - fl1_fx * fl1_fza * tg_y_xxxyyyy_1[j] + 2.0 * fl1_fxn * tg_yy_xxxyyy_1[j];

                    tg_yyy_xxxyyyz_0[j] = pb_y * tg_yy_xxxyyyz_0[j] + wp_y[j] * tg_yy_xxxyyyz_1[j] + fl1_fx * tg_y_xxxyyyz_0[j] - fl1_fx * fl1_fza * tg_y_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_yy_xxxyyz_1[j];

                    tg_yyy_xxxyyzz_0[j] = pb_y * tg_yy_xxxyyzz_0[j] + wp_y[j] * tg_yy_xxxyyzz_1[j] + fl1_fx * tg_y_xxxyyzz_0[j] - fl1_fx * fl1_fza * tg_y_xxxyyzz_1[j] + fl1_fxn * tg_yy_xxxyzz_1[j];

                    tg_yyy_xxxyzzz_0[j] = pb_y * tg_yy_xxxyzzz_0[j] + wp_y[j] * tg_yy_xxxyzzz_1[j] + fl1_fx * tg_y_xxxyzzz_0[j] - fl1_fx * fl1_fza * tg_y_xxxyzzz_1[j] + 0.5 * fl1_fxn * tg_yy_xxxzzz_1[j];

                    tg_yyy_xxxzzzz_0[j] = pb_y * tg_yy_xxxzzzz_0[j] + wp_y[j] * tg_yy_xxxzzzz_1[j] + fl1_fx * tg_y_xxxzzzz_0[j] - fl1_fx * fl1_fza * tg_y_xxxzzzz_1[j];

                    tg_yyy_xxyyyyy_0[j] = pb_y * tg_yy_xxyyyyy_0[j] + wp_y[j] * tg_yy_xxyyyyy_1[j] + fl1_fx * tg_y_xxyyyyy_0[j] - fl1_fx * fl1_fza * tg_y_xxyyyyy_1[j] + 2.5 * fl1_fxn * tg_yy_xxyyyy_1[j];

                    tg_yyy_xxyyyyz_0[j] = pb_y * tg_yy_xxyyyyz_0[j] + wp_y[j] * tg_yy_xxyyyyz_1[j] + fl1_fx * tg_y_xxyyyyz_0[j] - fl1_fx * fl1_fza * tg_y_xxyyyyz_1[j] + 2.0 * fl1_fxn * tg_yy_xxyyyz_1[j];

                    tg_yyy_xxyyyzz_0[j] = pb_y * tg_yy_xxyyyzz_0[j] + wp_y[j] * tg_yy_xxyyyzz_1[j] + fl1_fx * tg_y_xxyyyzz_0[j] - fl1_fx * fl1_fza * tg_y_xxyyyzz_1[j] + 1.5 * fl1_fxn * tg_yy_xxyyzz_1[j];

                    tg_yyy_xxyyzzz_0[j] = pb_y * tg_yy_xxyyzzz_0[j] + wp_y[j] * tg_yy_xxyyzzz_1[j] + fl1_fx * tg_y_xxyyzzz_0[j] - fl1_fx * fl1_fza * tg_y_xxyyzzz_1[j] + fl1_fxn * tg_yy_xxyzzz_1[j];

                    tg_yyy_xxyzzzz_0[j] = pb_y * tg_yy_xxyzzzz_0[j] + wp_y[j] * tg_yy_xxyzzzz_1[j] + fl1_fx * tg_y_xxyzzzz_0[j] - fl1_fx * fl1_fza * tg_y_xxyzzzz_1[j] + 0.5 * fl1_fxn * tg_yy_xxzzzz_1[j];

                    tg_yyy_xxzzzzz_0[j] = pb_y * tg_yy_xxzzzzz_0[j] + wp_y[j] * tg_yy_xxzzzzz_1[j] + fl1_fx * tg_y_xxzzzzz_0[j] - fl1_fx * fl1_fza * tg_y_xxzzzzz_1[j];

                    tg_yyy_xyyyyyy_0[j] = pb_y * tg_yy_xyyyyyy_0[j] + wp_y[j] * tg_yy_xyyyyyy_1[j] + fl1_fx * tg_y_xyyyyyy_0[j] - fl1_fx * fl1_fza * tg_y_xyyyyyy_1[j] + 3.0 * fl1_fxn * tg_yy_xyyyyy_1[j];

                    tg_yyy_xyyyyyz_0[j] = pb_y * tg_yy_xyyyyyz_0[j] + wp_y[j] * tg_yy_xyyyyyz_1[j] + fl1_fx * tg_y_xyyyyyz_0[j] - fl1_fx * fl1_fza * tg_y_xyyyyyz_1[j] + 2.5 * fl1_fxn * tg_yy_xyyyyz_1[j];

                    tg_yyy_xyyyyzz_0[j] = pb_y * tg_yy_xyyyyzz_0[j] + wp_y[j] * tg_yy_xyyyyzz_1[j] + fl1_fx * tg_y_xyyyyzz_0[j] - fl1_fx * fl1_fza * tg_y_xyyyyzz_1[j] + 2.0 * fl1_fxn * tg_yy_xyyyzz_1[j];

                    tg_yyy_xyyyzzz_0[j] = pb_y * tg_yy_xyyyzzz_0[j] + wp_y[j] * tg_yy_xyyyzzz_1[j] + fl1_fx * tg_y_xyyyzzz_0[j] - fl1_fx * fl1_fza * tg_y_xyyyzzz_1[j] + 1.5 * fl1_fxn * tg_yy_xyyzzz_1[j];

                    tg_yyy_xyyzzzz_0[j] = pb_y * tg_yy_xyyzzzz_0[j] + wp_y[j] * tg_yy_xyyzzzz_1[j] + fl1_fx * tg_y_xyyzzzz_0[j] - fl1_fx * fl1_fza * tg_y_xyyzzzz_1[j] + fl1_fxn * tg_yy_xyzzzz_1[j];

                    tg_yyy_xyzzzzz_0[j] = pb_y * tg_yy_xyzzzzz_0[j] + wp_y[j] * tg_yy_xyzzzzz_1[j] + fl1_fx * tg_y_xyzzzzz_0[j] - fl1_fx * fl1_fza * tg_y_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_yy_xzzzzz_1[j];

                    tg_yyy_xzzzzzz_0[j] = pb_y * tg_yy_xzzzzzz_0[j] + wp_y[j] * tg_yy_xzzzzzz_1[j] + fl1_fx * tg_y_xzzzzzz_0[j] - fl1_fx * fl1_fza * tg_y_xzzzzzz_1[j];

                    tg_yyy_yyyyyyy_0[j] = pb_y * tg_yy_yyyyyyy_0[j] + wp_y[j] * tg_yy_yyyyyyy_1[j] + fl1_fx * tg_y_yyyyyyy_0[j] - fl1_fx * fl1_fza * tg_y_yyyyyyy_1[j] + 3.5 * fl1_fxn * tg_yy_yyyyyy_1[j];

                    tg_yyy_yyyyyyz_0[j] = pb_y * tg_yy_yyyyyyz_0[j] + wp_y[j] * tg_yy_yyyyyyz_1[j] + fl1_fx * tg_y_yyyyyyz_0[j] - fl1_fx * fl1_fza * tg_y_yyyyyyz_1[j] + 3.0 * fl1_fxn * tg_yy_yyyyyz_1[j];

                    tg_yyy_yyyyyzz_0[j] = pb_y * tg_yy_yyyyyzz_0[j] + wp_y[j] * tg_yy_yyyyyzz_1[j] + fl1_fx * tg_y_yyyyyzz_0[j] - fl1_fx * fl1_fza * tg_y_yyyyyzz_1[j] + 2.5 * fl1_fxn * tg_yy_yyyyzz_1[j];

                    tg_yyy_yyyyzzz_0[j] = pb_y * tg_yy_yyyyzzz_0[j] + wp_y[j] * tg_yy_yyyyzzz_1[j] + fl1_fx * tg_y_yyyyzzz_0[j] - fl1_fx * fl1_fza * tg_y_yyyyzzz_1[j] + 2.0 * fl1_fxn * tg_yy_yyyzzz_1[j];

                    tg_yyy_yyyzzzz_0[j] = pb_y * tg_yy_yyyzzzz_0[j] + wp_y[j] * tg_yy_yyyzzzz_1[j] + fl1_fx * tg_y_yyyzzzz_0[j] - fl1_fx * fl1_fza * tg_y_yyyzzzz_1[j] + 1.5 * fl1_fxn * tg_yy_yyzzzz_1[j];

                    tg_yyy_yyzzzzz_0[j] = pb_y * tg_yy_yyzzzzz_0[j] + wp_y[j] * tg_yy_yyzzzzz_1[j] + fl1_fx * tg_y_yyzzzzz_0[j] - fl1_fx * fl1_fza * tg_y_yyzzzzz_1[j] + fl1_fxn * tg_yy_yzzzzz_1[j];

                    tg_yyy_yzzzzzz_0[j] = pb_y * tg_yy_yzzzzzz_0[j] + wp_y[j] * tg_yy_yzzzzzz_1[j] + fl1_fx * tg_y_yzzzzzz_0[j] - fl1_fx * fl1_fza * tg_y_yzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yy_zzzzzz_1[j];

                    tg_yyy_zzzzzzz_0[j] = pb_y * tg_yy_zzzzzzz_0[j] + wp_y[j] * tg_yy_zzzzzzz_1[j] + fl1_fx * tg_y_zzzzzzz_0[j] - fl1_fx * fl1_fza * tg_y_zzzzzzz_1[j];

                    tg_yyz_xxxxxxx_0[j] = pb_y * tg_yz_xxxxxxx_0[j] + wp_y[j] * tg_yz_xxxxxxx_1[j] + 0.5 * fl1_fx * tg_z_xxxxxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxxx_1[j];

                    tg_yyz_xxxxxxy_0[j] = pb_y * tg_yz_xxxxxxy_0[j] + wp_y[j] * tg_yz_xxxxxxy_1[j] + 0.5 * fl1_fx * tg_z_xxxxxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxxy_1[j] + 0.5 * fl1_fxn * tg_yz_xxxxxx_1[j];

                    tg_yyz_xxxxxxz_0[j] = pb_y * tg_yz_xxxxxxz_0[j] + wp_y[j] * tg_yz_xxxxxxz_1[j] + 0.5 * fl1_fx * tg_z_xxxxxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxxz_1[j];

                    tg_yyz_xxxxxyy_0[j] = pb_y * tg_yz_xxxxxyy_0[j] + wp_y[j] * tg_yz_xxxxxyy_1[j] + 0.5 * fl1_fx * tg_z_xxxxxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxyy_1[j] + fl1_fxn * tg_yz_xxxxxy_1[j];

                    tg_yyz_xxxxxyz_0[j] = pb_y * tg_yz_xxxxxyz_0[j] + wp_y[j] * tg_yz_xxxxxyz_1[j] + 0.5 * fl1_fx * tg_z_xxxxxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxyz_1[j] + 0.5 * fl1_fxn * tg_yz_xxxxxz_1[j];

                    tg_yyz_xxxxxzz_0[j] = pb_y * tg_yz_xxxxxzz_0[j] + wp_y[j] * tg_yz_xxxxxzz_1[j] + 0.5 * fl1_fx * tg_z_xxxxxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxzz_1[j];

                    tg_yyz_xxxxyyy_0[j] = pb_y * tg_yz_xxxxyyy_0[j] + wp_y[j] * tg_yz_xxxxyyy_1[j] + 0.5 * fl1_fx * tg_z_xxxxyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxyyy_1[j] + 1.5 * fl1_fxn * tg_yz_xxxxyy_1[j];

                    tg_yyz_xxxxyyz_0[j] = pb_y * tg_yz_xxxxyyz_0[j] + wp_y[j] * tg_yz_xxxxyyz_1[j] + 0.5 * fl1_fx * tg_z_xxxxyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxyyz_1[j] + fl1_fxn * tg_yz_xxxxyz_1[j];

                    tg_yyz_xxxxyzz_0[j] = pb_y * tg_yz_xxxxyzz_0[j] + wp_y[j] * tg_yz_xxxxyzz_1[j] + 0.5 * fl1_fx * tg_z_xxxxyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxyzz_1[j] + 0.5 * fl1_fxn * tg_yz_xxxxzz_1[j];

                    tg_yyz_xxxxzzz_0[j] = pb_y * tg_yz_xxxxzzz_0[j] + wp_y[j] * tg_yz_xxxxzzz_1[j] + 0.5 * fl1_fx * tg_z_xxxxzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxzzz_1[j];

                    tg_yyz_xxxyyyy_0[j] = pb_y * tg_yz_xxxyyyy_0[j] + wp_y[j] * tg_yz_xxxyyyy_1[j] + 0.5 * fl1_fx * tg_z_xxxyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxyyyy_1[j] + 2.0 * fl1_fxn * tg_yz_xxxyyy_1[j];

                    tg_yyz_xxxyyyz_0[j] = pb_y * tg_yz_xxxyyyz_0[j] + wp_y[j] * tg_yz_xxxyyyz_1[j] + 0.5 * fl1_fx * tg_z_xxxyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_yz_xxxyyz_1[j];

                    tg_yyz_xxxyyzz_0[j] = pb_y * tg_yz_xxxyyzz_0[j] + wp_y[j] * tg_yz_xxxyyzz_1[j] + 0.5 * fl1_fx * tg_z_xxxyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxyyzz_1[j] + fl1_fxn * tg_yz_xxxyzz_1[j];

                    tg_yyz_xxxyzzz_0[j] = pb_y * tg_yz_xxxyzzz_0[j] + wp_y[j] * tg_yz_xxxyzzz_1[j] + 0.5 * fl1_fx * tg_z_xxxyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxyzzz_1[j] + 0.5 * fl1_fxn * tg_yz_xxxzzz_1[j];

                    tg_yyz_xxxzzzz_0[j] = pb_y * tg_yz_xxxzzzz_0[j] + wp_y[j] * tg_yz_xxxzzzz_1[j] + 0.5 * fl1_fx * tg_z_xxxzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxzzzz_1[j];

                    tg_yyz_xxyyyyy_0[j] = pb_y * tg_yz_xxyyyyy_0[j] + wp_y[j] * tg_yz_xxyyyyy_1[j] + 0.5 * fl1_fx * tg_z_xxyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxyyyyy_1[j] + 2.5 * fl1_fxn * tg_yz_xxyyyy_1[j];

                    tg_yyz_xxyyyyz_0[j] = pb_y * tg_yz_xxyyyyz_0[j] + wp_y[j] * tg_yz_xxyyyyz_1[j] + 0.5 * fl1_fx * tg_z_xxyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxyyyyz_1[j] + 2.0 * fl1_fxn * tg_yz_xxyyyz_1[j];

                    tg_yyz_xxyyyzz_0[j] = pb_y * tg_yz_xxyyyzz_0[j] + wp_y[j] * tg_yz_xxyyyzz_1[j] + 0.5 * fl1_fx * tg_z_xxyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxyyyzz_1[j] + 1.5 * fl1_fxn * tg_yz_xxyyzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSFSK_270_360(      CMemBlock2D<double>& primBuffer,
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

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        auto r_pb_z = braGtoPairsBlock.getDistancesPBZ();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {3, -1, -1, -1},
                                             {7, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_3_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_3_7_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_2_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_2_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_1_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_1_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_2_6_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {6, -1, -1, -1}, 
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

                auto tg_yz_xxyyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 162); 

                auto tg_yz_xxyzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 163); 

                auto tg_yz_xxzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 164); 

                auto tg_yz_xyyyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 165); 

                auto tg_yz_xyyyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 166); 

                auto tg_yz_xyyyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 167); 

                auto tg_yz_xyyyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 168); 

                auto tg_yz_xyyzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 169); 

                auto tg_yz_xyzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 170); 

                auto tg_yz_xzzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 171); 

                auto tg_yz_yyyyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 172); 

                auto tg_yz_yyyyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 173); 

                auto tg_yz_yyyyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 174); 

                auto tg_yz_yyyyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 175); 

                auto tg_yz_yyyzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 176); 

                auto tg_yz_yyzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 177); 

                auto tg_yz_yzzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 178); 

                auto tg_yz_zzzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 179); 

                auto tg_zz_xxxxxxx_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 180); 

                auto tg_zz_xxxxxxy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 181); 

                auto tg_zz_xxxxxxz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 182); 

                auto tg_zz_xxxxxyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 183); 

                auto tg_zz_xxxxxyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 184); 

                auto tg_zz_xxxxxzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 185); 

                auto tg_zz_xxxxyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 186); 

                auto tg_zz_xxxxyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 187); 

                auto tg_zz_xxxxyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 188); 

                auto tg_zz_xxxxzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 189); 

                auto tg_zz_xxxyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 190); 

                auto tg_zz_xxxyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 191); 

                auto tg_zz_xxxyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 192); 

                auto tg_zz_xxxyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 193); 

                auto tg_zz_xxxzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 194); 

                auto tg_zz_xxyyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 195); 

                auto tg_zz_xxyyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 196); 

                auto tg_zz_xxyyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 197); 

                auto tg_zz_xxyyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 198); 

                auto tg_zz_xxyzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 199); 

                auto tg_zz_xxzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 200); 

                auto tg_zz_xyyyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 201); 

                auto tg_zz_xyyyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 202); 

                auto tg_zz_xyyyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 203); 

                auto tg_zz_xyyyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 204); 

                auto tg_zz_xyyzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 205); 

                auto tg_zz_xyzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 206); 

                auto tg_zz_xzzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 207); 

                auto tg_zz_yyyyyyy_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 208); 

                auto tg_zz_yyyyyyz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 209); 

                auto tg_zz_yyyyyzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 210); 

                auto tg_zz_yyyyzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 211); 

                auto tg_zz_yyyzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 212); 

                auto tg_zz_yyzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 213); 

                auto tg_zz_yzzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 214); 

                auto tg_zz_zzzzzzz_0 = primBuffer.data(pidx_g_2_7_m0 + 216 * idx + 215); 

                auto tg_yz_xxyyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 162); 

                auto tg_yz_xxyzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 163); 

                auto tg_yz_xxzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 164); 

                auto tg_yz_xyyyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 165); 

                auto tg_yz_xyyyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 166); 

                auto tg_yz_xyyyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 167); 

                auto tg_yz_xyyyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 168); 

                auto tg_yz_xyyzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 169); 

                auto tg_yz_xyzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 170); 

                auto tg_yz_xzzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 171); 

                auto tg_yz_yyyyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 172); 

                auto tg_yz_yyyyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 173); 

                auto tg_yz_yyyyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 174); 

                auto tg_yz_yyyyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 175); 

                auto tg_yz_yyyzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 176); 

                auto tg_yz_yyzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 177); 

                auto tg_yz_yzzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 178); 

                auto tg_yz_zzzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 179); 

                auto tg_zz_xxxxxxx_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 180); 

                auto tg_zz_xxxxxxy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 181); 

                auto tg_zz_xxxxxxz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 182); 

                auto tg_zz_xxxxxyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 183); 

                auto tg_zz_xxxxxyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 184); 

                auto tg_zz_xxxxxzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 185); 

                auto tg_zz_xxxxyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 186); 

                auto tg_zz_xxxxyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 187); 

                auto tg_zz_xxxxyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 188); 

                auto tg_zz_xxxxzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 189); 

                auto tg_zz_xxxyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 190); 

                auto tg_zz_xxxyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 191); 

                auto tg_zz_xxxyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 192); 

                auto tg_zz_xxxyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 193); 

                auto tg_zz_xxxzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 194); 

                auto tg_zz_xxyyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 195); 

                auto tg_zz_xxyyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 196); 

                auto tg_zz_xxyyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 197); 

                auto tg_zz_xxyyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 198); 

                auto tg_zz_xxyzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 199); 

                auto tg_zz_xxzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 200); 

                auto tg_zz_xyyyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 201); 

                auto tg_zz_xyyyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 202); 

                auto tg_zz_xyyyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 203); 

                auto tg_zz_xyyyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 204); 

                auto tg_zz_xyyzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 205); 

                auto tg_zz_xyzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 206); 

                auto tg_zz_xzzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 207); 

                auto tg_zz_yyyyyyy_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 208); 

                auto tg_zz_yyyyyyz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 209); 

                auto tg_zz_yyyyyzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 210); 

                auto tg_zz_yyyyzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 211); 

                auto tg_zz_yyyzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 212); 

                auto tg_zz_yyzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 213); 

                auto tg_zz_yzzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 214); 

                auto tg_zz_zzzzzzz_1 = primBuffer.data(pidx_g_2_7_m1 + 216 * idx + 215); 

                auto tg_z_xxxxxxx_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 72); 

                auto tg_z_xxxxxxy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 73); 

                auto tg_z_xxxxxxz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 74); 

                auto tg_z_xxxxxyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 75); 

                auto tg_z_xxxxxyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 76); 

                auto tg_z_xxxxxzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 77); 

                auto tg_z_xxxxyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 78); 

                auto tg_z_xxxxyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 79); 

                auto tg_z_xxxxyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 80); 

                auto tg_z_xxxxzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 81); 

                auto tg_z_xxxyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 82); 

                auto tg_z_xxxyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 83); 

                auto tg_z_xxxyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 84); 

                auto tg_z_xxxyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 85); 

                auto tg_z_xxxzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 86); 

                auto tg_z_xxyyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 87); 

                auto tg_z_xxyyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 88); 

                auto tg_z_xxyyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 89); 

                auto tg_z_xxyyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 90); 

                auto tg_z_xxyzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 91); 

                auto tg_z_xxzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 92); 

                auto tg_z_xyyyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 93); 

                auto tg_z_xyyyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 94); 

                auto tg_z_xyyyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 95); 

                auto tg_z_xyyyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 96); 

                auto tg_z_xyyzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 97); 

                auto tg_z_xyzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 98); 

                auto tg_z_xzzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 99); 

                auto tg_z_yyyyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 100); 

                auto tg_z_yyyyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 101); 

                auto tg_z_yyyyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 102); 

                auto tg_z_yyyyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 103); 

                auto tg_z_yyyzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 104); 

                auto tg_z_yyzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 105); 

                auto tg_z_yzzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 106); 

                auto tg_z_zzzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 107); 

                auto tg_z_xxxxxxx_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 72); 

                auto tg_z_xxxxxxy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 73); 

                auto tg_z_xxxxxxz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 74); 

                auto tg_z_xxxxxyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 75); 

                auto tg_z_xxxxxyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 76); 

                auto tg_z_xxxxxzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 77); 

                auto tg_z_xxxxyyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 78); 

                auto tg_z_xxxxyyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 79); 

                auto tg_z_xxxxyzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 80); 

                auto tg_z_xxxxzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 81); 

                auto tg_z_xxxyyyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 82); 

                auto tg_z_xxxyyyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 83); 

                auto tg_z_xxxyyzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 84); 

                auto tg_z_xxxyzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 85); 

                auto tg_z_xxxzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 86); 

                auto tg_z_xxyyyyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 87); 

                auto tg_z_xxyyyyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 88); 

                auto tg_z_xxyyyzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 89); 

                auto tg_z_xxyyzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 90); 

                auto tg_z_xxyzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 91); 

                auto tg_z_xxzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 92); 

                auto tg_z_xyyyyyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 93); 

                auto tg_z_xyyyyyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 94); 

                auto tg_z_xyyyyzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 95); 

                auto tg_z_xyyyzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 96); 

                auto tg_z_xyyzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 97); 

                auto tg_z_xyzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 98); 

                auto tg_z_xzzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 99); 

                auto tg_z_yyyyyyy_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 100); 

                auto tg_z_yyyyyyz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 101); 

                auto tg_z_yyyyyzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 102); 

                auto tg_z_yyyyzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 103); 

                auto tg_z_yyyzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 104); 

                auto tg_z_yyzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 105); 

                auto tg_z_yzzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 106); 

                auto tg_z_zzzzzzz_1 = primBuffer.data(pidx_g_1_7_m1 + 108 * idx + 107); 

                auto tg_yz_xxyzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 125); 

                auto tg_yz_xxzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 126); 

                auto tg_yz_xyyyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 127); 

                auto tg_yz_xyyyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 128); 

                auto tg_yz_xyyyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 129); 

                auto tg_yz_xyyzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 130); 

                auto tg_yz_xyzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 131); 

                auto tg_yz_xzzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 132); 

                auto tg_yz_yyyyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 133); 

                auto tg_yz_yyyyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 134); 

                auto tg_yz_yyyyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 135); 

                auto tg_yz_yyyzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 136); 

                auto tg_yz_yyzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 137); 

                auto tg_yz_yzzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 138); 

                auto tg_yz_zzzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 139); 

                auto tg_zz_xxxxxx_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 140); 

                auto tg_zz_xxxxxy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 141); 

                auto tg_zz_xxxxxz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 142); 

                auto tg_zz_xxxxyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 143); 

                auto tg_zz_xxxxyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 144); 

                auto tg_zz_xxxxzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 145); 

                auto tg_zz_xxxyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 146); 

                auto tg_zz_xxxyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 147); 

                auto tg_zz_xxxyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 148); 

                auto tg_zz_xxxzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 149); 

                auto tg_zz_xxyyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 150); 

                auto tg_zz_xxyyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 151); 

                auto tg_zz_xxyyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 152); 

                auto tg_zz_xxyzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 153); 

                auto tg_zz_xxzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 154); 

                auto tg_zz_xyyyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 155); 

                auto tg_zz_xyyyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 156); 

                auto tg_zz_xyyyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 157); 

                auto tg_zz_xyyzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 158); 

                auto tg_zz_xyzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 159); 

                auto tg_zz_xzzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 160); 

                auto tg_zz_yyyyyy_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 161); 

                auto tg_zz_yyyyyz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 162); 

                auto tg_zz_yyyyzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 163); 

                auto tg_zz_yyyzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 164); 

                auto tg_zz_yyzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 165); 

                auto tg_zz_yzzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 166); 

                auto tg_zz_zzzzzz_1 = primBuffer.data(pidx_g_2_6_m1 + 168 * idx + 167); 

                // set up pointers to integrals

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

                // Batch of Integrals (270,360)

                #pragma omp simd aligned(fxn, fza, tg_yyz_xxyyzzz_0, tg_yyz_xxyzzzz_0, tg_yyz_xxzzzzz_0, \
                                         tg_yyz_xyyyyyy_0, tg_yyz_xyyyyyz_0, tg_yyz_xyyyyzz_0, tg_yyz_xyyyzzz_0, \
                                         tg_yyz_xyyzzzz_0, tg_yyz_xyzzzzz_0, tg_yyz_xzzzzzz_0, tg_yyz_yyyyyyy_0, \
                                         tg_yyz_yyyyyyz_0, tg_yyz_yyyyyzz_0, tg_yyz_yyyyzzz_0, tg_yyz_yyyzzzz_0, \
                                         tg_yyz_yyzzzzz_0, tg_yyz_yzzzzzz_0, tg_yyz_zzzzzzz_0, tg_yz_xxyyzzz_0, \
                                         tg_yz_xxyyzzz_1, tg_yz_xxyzzz_1, tg_yz_xxyzzzz_0, tg_yz_xxyzzzz_1, tg_yz_xxzzzz_1, \
                                         tg_yz_xxzzzzz_0, tg_yz_xxzzzzz_1, tg_yz_xyyyyy_1, tg_yz_xyyyyyy_0, tg_yz_xyyyyyy_1, \
                                         tg_yz_xyyyyyz_0, tg_yz_xyyyyyz_1, tg_yz_xyyyyz_1, tg_yz_xyyyyzz_0, tg_yz_xyyyyzz_1, \
                                         tg_yz_xyyyzz_1, tg_yz_xyyyzzz_0, tg_yz_xyyyzzz_1, tg_yz_xyyzzz_1, tg_yz_xyyzzzz_0, \
                                         tg_yz_xyyzzzz_1, tg_yz_xyzzzz_1, tg_yz_xyzzzzz_0, tg_yz_xyzzzzz_1, tg_yz_xzzzzz_1, \
                                         tg_yz_xzzzzzz_0, tg_yz_xzzzzzz_1, tg_yz_yyyyyy_1, tg_yz_yyyyyyy_0, tg_yz_yyyyyyy_1, \
                                         tg_yz_yyyyyyz_0, tg_yz_yyyyyyz_1, tg_yz_yyyyyz_1, tg_yz_yyyyyzz_0, tg_yz_yyyyyzz_1, \
                                         tg_yz_yyyyzz_1, tg_yz_yyyyzzz_0, tg_yz_yyyyzzz_1, tg_yz_yyyzzz_1, tg_yz_yyyzzzz_0, \
                                         tg_yz_yyyzzzz_1, tg_yz_yyzzzz_1, tg_yz_yyzzzzz_0, tg_yz_yyzzzzz_1, tg_yz_yzzzzz_1, \
                                         tg_yz_yzzzzzz_0, tg_yz_yzzzzzz_1, tg_yz_zzzzzz_1, tg_yz_zzzzzzz_0, tg_yz_zzzzzzz_1, \
                                         tg_yzz_xxxxxxx_0, tg_yzz_xxxxxxy_0, tg_yzz_xxxxxxz_0, tg_yzz_xxxxxyy_0, \
                                         tg_yzz_xxxxxyz_0, tg_yzz_xxxxxzz_0, tg_yzz_xxxxyyy_0, tg_yzz_xxxxyyz_0, \
                                         tg_yzz_xxxxyzz_0, tg_yzz_xxxxzzz_0, tg_yzz_xxxyyyy_0, tg_yzz_xxxyyyz_0, \
                                         tg_yzz_xxxyyzz_0, tg_yzz_xxxyzzz_0, tg_yzz_xxxzzzz_0, tg_yzz_xxyyyyy_0, \
                                         tg_yzz_xxyyyyz_0, tg_yzz_xxyyyzz_0, tg_yzz_xxyyzzz_0, tg_yzz_xxyzzzz_0, \
                                         tg_yzz_xxzzzzz_0, tg_yzz_xyyyyyy_0, tg_yzz_xyyyyyz_0, tg_yzz_xyyyyzz_0, \
                                         tg_yzz_xyyyzzz_0, tg_yzz_xyyzzzz_0, tg_yzz_xyzzzzz_0, tg_yzz_xzzzzzz_0, \
                                         tg_yzz_yyyyyyy_0, tg_yzz_yyyyyyz_0, tg_yzz_yyyyyzz_0, tg_yzz_yyyyzzz_0, \
                                         tg_yzz_yyyzzzz_0, tg_yzz_yyzzzzz_0, tg_yzz_yzzzzzz_0, tg_yzz_zzzzzzz_0, \
                                         tg_z_xxxxxxx_0, tg_z_xxxxxxx_1, tg_z_xxxxxxy_0, tg_z_xxxxxxy_1, tg_z_xxxxxxz_0, \
                                         tg_z_xxxxxxz_1, tg_z_xxxxxyy_0, tg_z_xxxxxyy_1, tg_z_xxxxxyz_0, tg_z_xxxxxyz_1, \
                                         tg_z_xxxxxzz_0, tg_z_xxxxxzz_1, tg_z_xxxxyyy_0, tg_z_xxxxyyy_1, tg_z_xxxxyyz_0, \
                                         tg_z_xxxxyyz_1, tg_z_xxxxyzz_0, tg_z_xxxxyzz_1, tg_z_xxxxzzz_0, tg_z_xxxxzzz_1, \
                                         tg_z_xxxyyyy_0, tg_z_xxxyyyy_1, tg_z_xxxyyyz_0, tg_z_xxxyyyz_1, tg_z_xxxyyzz_0, \
                                         tg_z_xxxyyzz_1, tg_z_xxxyzzz_0, tg_z_xxxyzzz_1, tg_z_xxxzzzz_0, tg_z_xxxzzzz_1, \
                                         tg_z_xxyyyyy_0, tg_z_xxyyyyy_1, tg_z_xxyyyyz_0, tg_z_xxyyyyz_1, tg_z_xxyyyzz_0, \
                                         tg_z_xxyyyzz_1, tg_z_xxyyzzz_0, tg_z_xxyyzzz_1, tg_z_xxyzzzz_0, tg_z_xxyzzzz_1, \
                                         tg_z_xxzzzzz_0, tg_z_xxzzzzz_1, tg_z_xyyyyyy_0, tg_z_xyyyyyy_1, tg_z_xyyyyyz_0, \
                                         tg_z_xyyyyyz_1, tg_z_xyyyyzz_0, tg_z_xyyyyzz_1, tg_z_xyyyzzz_0, tg_z_xyyyzzz_1, \
                                         tg_z_xyyzzzz_0, tg_z_xyyzzzz_1, tg_z_xyzzzzz_0, tg_z_xyzzzzz_1, tg_z_xzzzzzz_0, \
                                         tg_z_xzzzzzz_1, tg_z_yyyyyyy_0, tg_z_yyyyyyy_1, tg_z_yyyyyyz_0, tg_z_yyyyyyz_1, \
                                         tg_z_yyyyyzz_0, tg_z_yyyyyzz_1, tg_z_yyyyzzz_0, tg_z_yyyyzzz_1, tg_z_yyyzzzz_0, \
                                         tg_z_yyyzzzz_1, tg_z_yyzzzzz_0, tg_z_yyzzzzz_1, tg_z_yzzzzzz_0, tg_z_yzzzzzz_1, \
                                         tg_z_zzzzzzz_0, tg_z_zzzzzzz_1, tg_zz_xxxxxx_1, tg_zz_xxxxxxx_0, tg_zz_xxxxxxx_1, \
                                         tg_zz_xxxxxxy_0, tg_zz_xxxxxxy_1, tg_zz_xxxxxxz_0, tg_zz_xxxxxxz_1, tg_zz_xxxxxy_1, \
                                         tg_zz_xxxxxyy_0, tg_zz_xxxxxyy_1, tg_zz_xxxxxyz_0, tg_zz_xxxxxyz_1, tg_zz_xxxxxz_1, \
                                         tg_zz_xxxxxzz_0, tg_zz_xxxxxzz_1, tg_zz_xxxxyy_1, tg_zz_xxxxyyy_0, tg_zz_xxxxyyy_1, \
                                         tg_zz_xxxxyyz_0, tg_zz_xxxxyyz_1, tg_zz_xxxxyz_1, tg_zz_xxxxyzz_0, tg_zz_xxxxyzz_1, \
                                         tg_zz_xxxxzz_1, tg_zz_xxxxzzz_0, tg_zz_xxxxzzz_1, tg_zz_xxxyyy_1, tg_zz_xxxyyyy_0, \
                                         tg_zz_xxxyyyy_1, tg_zz_xxxyyyz_0, tg_zz_xxxyyyz_1, tg_zz_xxxyyz_1, tg_zz_xxxyyzz_0, \
                                         tg_zz_xxxyyzz_1, tg_zz_xxxyzz_1, tg_zz_xxxyzzz_0, tg_zz_xxxyzzz_1, tg_zz_xxxzzz_1, \
                                         tg_zz_xxxzzzz_0, tg_zz_xxxzzzz_1, tg_zz_xxyyyy_1, tg_zz_xxyyyyy_0, tg_zz_xxyyyyy_1, \
                                         tg_zz_xxyyyyz_0, tg_zz_xxyyyyz_1, tg_zz_xxyyyz_1, tg_zz_xxyyyzz_0, tg_zz_xxyyyzz_1, \
                                         tg_zz_xxyyzz_1, tg_zz_xxyyzzz_0, tg_zz_xxyyzzz_1, tg_zz_xxyzzz_1, tg_zz_xxyzzzz_0, \
                                         tg_zz_xxyzzzz_1, tg_zz_xxzzzz_1, tg_zz_xxzzzzz_0, tg_zz_xxzzzzz_1, tg_zz_xyyyyy_1, \
                                         tg_zz_xyyyyyy_0, tg_zz_xyyyyyy_1, tg_zz_xyyyyyz_0, tg_zz_xyyyyyz_1, tg_zz_xyyyyz_1, \
                                         tg_zz_xyyyyzz_0, tg_zz_xyyyyzz_1, tg_zz_xyyyzz_1, tg_zz_xyyyzzz_0, tg_zz_xyyyzzz_1, \
                                         tg_zz_xyyzzz_1, tg_zz_xyyzzzz_0, tg_zz_xyyzzzz_1, tg_zz_xyzzzz_1, tg_zz_xyzzzzz_0, \
                                         tg_zz_xyzzzzz_1, tg_zz_xzzzzz_1, tg_zz_xzzzzzz_0, tg_zz_xzzzzzz_1, tg_zz_yyyyyy_1, \
                                         tg_zz_yyyyyyy_0, tg_zz_yyyyyyy_1, tg_zz_yyyyyyz_0, tg_zz_yyyyyyz_1, tg_zz_yyyyyz_1, \
                                         tg_zz_yyyyyzz_0, tg_zz_yyyyyzz_1, tg_zz_yyyyzz_1, tg_zz_yyyyzzz_0, tg_zz_yyyyzzz_1, \
                                         tg_zz_yyyzzz_1, tg_zz_yyyzzzz_0, tg_zz_yyyzzzz_1, tg_zz_yyzzzz_1, tg_zz_yyzzzzz_0, \
                                         tg_zz_yyzzzzz_1, tg_zz_yzzzzz_1, tg_zz_yzzzzzz_0, tg_zz_yzzzzzz_1, tg_zz_zzzzzz_1, \
                                         tg_zz_zzzzzzz_0, tg_zz_zzzzzzz_1, tg_zzz_xxxxxxx_0, tg_zzz_xxxxxxy_0, \
                                         tg_zzz_xxxxxxz_0, tg_zzz_xxxxxyy_0, tg_zzz_xxxxxyz_0, tg_zzz_xxxxxzz_0, \
                                         tg_zzz_xxxxyyy_0, tg_zzz_xxxxyyz_0, tg_zzz_xxxxyzz_0, tg_zzz_xxxxzzz_0, \
                                         tg_zzz_xxxyyyy_0, tg_zzz_xxxyyyz_0, tg_zzz_xxxyyzz_0, tg_zzz_xxxyzzz_0, \
                                         tg_zzz_xxxzzzz_0, tg_zzz_xxyyyyy_0, tg_zzz_xxyyyyz_0, tg_zzz_xxyyyzz_0, \
                                         tg_zzz_xxyyzzz_0, tg_zzz_xxyzzzz_0, tg_zzz_xxzzzzz_0, tg_zzz_xyyyyyy_0, \
                                         tg_zzz_xyyyyyz_0, tg_zzz_xyyyyzz_0, tg_zzz_xyyyzzz_0, tg_zzz_xyyzzzz_0, \
                                         tg_zzz_xyzzzzz_0, tg_zzz_xzzzzzz_0, tg_zzz_yyyyyyy_0, tg_zzz_yyyyyyz_0, \
                                         tg_zzz_yyyyyzz_0, tg_zzz_yyyyzzz_0, tg_zzz_yyyzzzz_0, tg_zzz_yyzzzzz_0, \
                                         tg_zzz_yzzzzzz_0, tg_zzz_zzzzzzz_0, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_yyz_xxyyzzz_0[j] = pb_y * tg_yz_xxyyzzz_0[j] + wp_y[j] * tg_yz_xxyyzzz_1[j] + 0.5 * fl1_fx * tg_z_xxyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxyyzzz_1[j] + fl1_fxn * tg_yz_xxyzzz_1[j];

                    tg_yyz_xxyzzzz_0[j] = pb_y * tg_yz_xxyzzzz_0[j] + wp_y[j] * tg_yz_xxyzzzz_1[j] + 0.5 * fl1_fx * tg_z_xxyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxyzzzz_1[j] + 0.5 * fl1_fxn * tg_yz_xxzzzz_1[j];

                    tg_yyz_xxzzzzz_0[j] = pb_y * tg_yz_xxzzzzz_0[j] + wp_y[j] * tg_yz_xxzzzzz_1[j] + 0.5 * fl1_fx * tg_z_xxzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxzzzzz_1[j];

                    tg_yyz_xyyyyyy_0[j] = pb_y * tg_yz_xyyyyyy_0[j] + wp_y[j] * tg_yz_xyyyyyy_1[j] + 0.5 * fl1_fx * tg_z_xyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyyyyyy_1[j] + 3.0 * fl1_fxn * tg_yz_xyyyyy_1[j];

                    tg_yyz_xyyyyyz_0[j] = pb_y * tg_yz_xyyyyyz_0[j] + wp_y[j] * tg_yz_xyyyyyz_1[j] + 0.5 * fl1_fx * tg_z_xyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyyyyyz_1[j] + 2.5 * fl1_fxn * tg_yz_xyyyyz_1[j];

                    tg_yyz_xyyyyzz_0[j] = pb_y * tg_yz_xyyyyzz_0[j] + wp_y[j] * tg_yz_xyyyyzz_1[j] + 0.5 * fl1_fx * tg_z_xyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyyyyzz_1[j] + 2.0 * fl1_fxn * tg_yz_xyyyzz_1[j];

                    tg_yyz_xyyyzzz_0[j] = pb_y * tg_yz_xyyyzzz_0[j] + wp_y[j] * tg_yz_xyyyzzz_1[j] + 0.5 * fl1_fx * tg_z_xyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyyyzzz_1[j] + 1.5 * fl1_fxn * tg_yz_xyyzzz_1[j];

                    tg_yyz_xyyzzzz_0[j] = pb_y * tg_yz_xyyzzzz_0[j] + wp_y[j] * tg_yz_xyyzzzz_1[j] + 0.5 * fl1_fx * tg_z_xyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyyzzzz_1[j] + fl1_fxn * tg_yz_xyzzzz_1[j];

                    tg_yyz_xyzzzzz_0[j] = pb_y * tg_yz_xyzzzzz_0[j] + wp_y[j] * tg_yz_xyzzzzz_1[j] + 0.5 * fl1_fx * tg_z_xyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_yz_xzzzzz_1[j];

                    tg_yyz_xzzzzzz_0[j] = pb_y * tg_yz_xzzzzzz_0[j] + wp_y[j] * tg_yz_xzzzzzz_1[j] + 0.5 * fl1_fx * tg_z_xzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xzzzzzz_1[j];

                    tg_yyz_yyyyyyy_0[j] = pb_y * tg_yz_yyyyyyy_0[j] + wp_y[j] * tg_yz_yyyyyyy_1[j] + 0.5 * fl1_fx * tg_z_yyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyyyyyy_1[j] + 3.5 * fl1_fxn * tg_yz_yyyyyy_1[j];

                    tg_yyz_yyyyyyz_0[j] = pb_y * tg_yz_yyyyyyz_0[j] + wp_y[j] * tg_yz_yyyyyyz_1[j] + 0.5 * fl1_fx * tg_z_yyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyyyyyz_1[j] + 3.0 * fl1_fxn * tg_yz_yyyyyz_1[j];

                    tg_yyz_yyyyyzz_0[j] = pb_y * tg_yz_yyyyyzz_0[j] + wp_y[j] * tg_yz_yyyyyzz_1[j] + 0.5 * fl1_fx * tg_z_yyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyyyyzz_1[j] + 2.5 * fl1_fxn * tg_yz_yyyyzz_1[j];

                    tg_yyz_yyyyzzz_0[j] = pb_y * tg_yz_yyyyzzz_0[j] + wp_y[j] * tg_yz_yyyyzzz_1[j] + 0.5 * fl1_fx * tg_z_yyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyyyzzz_1[j] + 2.0 * fl1_fxn * tg_yz_yyyzzz_1[j];

                    tg_yyz_yyyzzzz_0[j] = pb_y * tg_yz_yyyzzzz_0[j] + wp_y[j] * tg_yz_yyyzzzz_1[j] + 0.5 * fl1_fx * tg_z_yyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyyzzzz_1[j] + 1.5 * fl1_fxn * tg_yz_yyzzzz_1[j];

                    tg_yyz_yyzzzzz_0[j] = pb_y * tg_yz_yyzzzzz_0[j] + wp_y[j] * tg_yz_yyzzzzz_1[j] + 0.5 * fl1_fx * tg_z_yyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyzzzzz_1[j] + fl1_fxn * tg_yz_yzzzzz_1[j];

                    tg_yyz_yzzzzzz_0[j] = pb_y * tg_yz_yzzzzzz_0[j] + wp_y[j] * tg_yz_yzzzzzz_1[j] + 0.5 * fl1_fx * tg_z_yzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yz_zzzzzz_1[j];

                    tg_yyz_zzzzzzz_0[j] = pb_y * tg_yz_zzzzzzz_0[j] + wp_y[j] * tg_yz_zzzzzzz_1[j] + 0.5 * fl1_fx * tg_z_zzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_zzzzzzz_1[j];

                    tg_yzz_xxxxxxx_0[j] = pb_y * tg_zz_xxxxxxx_0[j] + wp_y[j] * tg_zz_xxxxxxx_1[j];

                    tg_yzz_xxxxxxy_0[j] = pb_y * tg_zz_xxxxxxy_0[j] + wp_y[j] * tg_zz_xxxxxxy_1[j] + 0.5 * fl1_fxn * tg_zz_xxxxxx_1[j];

                    tg_yzz_xxxxxxz_0[j] = pb_y * tg_zz_xxxxxxz_0[j] + wp_y[j] * tg_zz_xxxxxxz_1[j];

                    tg_yzz_xxxxxyy_0[j] = pb_y * tg_zz_xxxxxyy_0[j] + wp_y[j] * tg_zz_xxxxxyy_1[j] + fl1_fxn * tg_zz_xxxxxy_1[j];

                    tg_yzz_xxxxxyz_0[j] = pb_y * tg_zz_xxxxxyz_0[j] + wp_y[j] * tg_zz_xxxxxyz_1[j] + 0.5 * fl1_fxn * tg_zz_xxxxxz_1[j];

                    tg_yzz_xxxxxzz_0[j] = pb_y * tg_zz_xxxxxzz_0[j] + wp_y[j] * tg_zz_xxxxxzz_1[j];

                    tg_yzz_xxxxyyy_0[j] = pb_y * tg_zz_xxxxyyy_0[j] + wp_y[j] * tg_zz_xxxxyyy_1[j] + 1.5 * fl1_fxn * tg_zz_xxxxyy_1[j];

                    tg_yzz_xxxxyyz_0[j] = pb_y * tg_zz_xxxxyyz_0[j] + wp_y[j] * tg_zz_xxxxyyz_1[j] + fl1_fxn * tg_zz_xxxxyz_1[j];

                    tg_yzz_xxxxyzz_0[j] = pb_y * tg_zz_xxxxyzz_0[j] + wp_y[j] * tg_zz_xxxxyzz_1[j] + 0.5 * fl1_fxn * tg_zz_xxxxzz_1[j];

                    tg_yzz_xxxxzzz_0[j] = pb_y * tg_zz_xxxxzzz_0[j] + wp_y[j] * tg_zz_xxxxzzz_1[j];

                    tg_yzz_xxxyyyy_0[j] = pb_y * tg_zz_xxxyyyy_0[j] + wp_y[j] * tg_zz_xxxyyyy_1[j] + 2.0 * fl1_fxn * tg_zz_xxxyyy_1[j];

                    tg_yzz_xxxyyyz_0[j] = pb_y * tg_zz_xxxyyyz_0[j] + wp_y[j] * tg_zz_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_zz_xxxyyz_1[j];

                    tg_yzz_xxxyyzz_0[j] = pb_y * tg_zz_xxxyyzz_0[j] + wp_y[j] * tg_zz_xxxyyzz_1[j] + fl1_fxn * tg_zz_xxxyzz_1[j];

                    tg_yzz_xxxyzzz_0[j] = pb_y * tg_zz_xxxyzzz_0[j] + wp_y[j] * tg_zz_xxxyzzz_1[j] + 0.5 * fl1_fxn * tg_zz_xxxzzz_1[j];

                    tg_yzz_xxxzzzz_0[j] = pb_y * tg_zz_xxxzzzz_0[j] + wp_y[j] * tg_zz_xxxzzzz_1[j];

                    tg_yzz_xxyyyyy_0[j] = pb_y * tg_zz_xxyyyyy_0[j] + wp_y[j] * tg_zz_xxyyyyy_1[j] + 2.5 * fl1_fxn * tg_zz_xxyyyy_1[j];

                    tg_yzz_xxyyyyz_0[j] = pb_y * tg_zz_xxyyyyz_0[j] + wp_y[j] * tg_zz_xxyyyyz_1[j] + 2.0 * fl1_fxn * tg_zz_xxyyyz_1[j];

                    tg_yzz_xxyyyzz_0[j] = pb_y * tg_zz_xxyyyzz_0[j] + wp_y[j] * tg_zz_xxyyyzz_1[j] + 1.5 * fl1_fxn * tg_zz_xxyyzz_1[j];

                    tg_yzz_xxyyzzz_0[j] = pb_y * tg_zz_xxyyzzz_0[j] + wp_y[j] * tg_zz_xxyyzzz_1[j] + fl1_fxn * tg_zz_xxyzzz_1[j];

                    tg_yzz_xxyzzzz_0[j] = pb_y * tg_zz_xxyzzzz_0[j] + wp_y[j] * tg_zz_xxyzzzz_1[j] + 0.5 * fl1_fxn * tg_zz_xxzzzz_1[j];

                    tg_yzz_xxzzzzz_0[j] = pb_y * tg_zz_xxzzzzz_0[j] + wp_y[j] * tg_zz_xxzzzzz_1[j];

                    tg_yzz_xyyyyyy_0[j] = pb_y * tg_zz_xyyyyyy_0[j] + wp_y[j] * tg_zz_xyyyyyy_1[j] + 3.0 * fl1_fxn * tg_zz_xyyyyy_1[j];

                    tg_yzz_xyyyyyz_0[j] = pb_y * tg_zz_xyyyyyz_0[j] + wp_y[j] * tg_zz_xyyyyyz_1[j] + 2.5 * fl1_fxn * tg_zz_xyyyyz_1[j];

                    tg_yzz_xyyyyzz_0[j] = pb_y * tg_zz_xyyyyzz_0[j] + wp_y[j] * tg_zz_xyyyyzz_1[j] + 2.0 * fl1_fxn * tg_zz_xyyyzz_1[j];

                    tg_yzz_xyyyzzz_0[j] = pb_y * tg_zz_xyyyzzz_0[j] + wp_y[j] * tg_zz_xyyyzzz_1[j] + 1.5 * fl1_fxn * tg_zz_xyyzzz_1[j];

                    tg_yzz_xyyzzzz_0[j] = pb_y * tg_zz_xyyzzzz_0[j] + wp_y[j] * tg_zz_xyyzzzz_1[j] + fl1_fxn * tg_zz_xyzzzz_1[j];

                    tg_yzz_xyzzzzz_0[j] = pb_y * tg_zz_xyzzzzz_0[j] + wp_y[j] * tg_zz_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_zz_xzzzzz_1[j];

                    tg_yzz_xzzzzzz_0[j] = pb_y * tg_zz_xzzzzzz_0[j] + wp_y[j] * tg_zz_xzzzzzz_1[j];

                    tg_yzz_yyyyyyy_0[j] = pb_y * tg_zz_yyyyyyy_0[j] + wp_y[j] * tg_zz_yyyyyyy_1[j] + 3.5 * fl1_fxn * tg_zz_yyyyyy_1[j];

                    tg_yzz_yyyyyyz_0[j] = pb_y * tg_zz_yyyyyyz_0[j] + wp_y[j] * tg_zz_yyyyyyz_1[j] + 3.0 * fl1_fxn * tg_zz_yyyyyz_1[j];

                    tg_yzz_yyyyyzz_0[j] = pb_y * tg_zz_yyyyyzz_0[j] + wp_y[j] * tg_zz_yyyyyzz_1[j] + 2.5 * fl1_fxn * tg_zz_yyyyzz_1[j];

                    tg_yzz_yyyyzzz_0[j] = pb_y * tg_zz_yyyyzzz_0[j] + wp_y[j] * tg_zz_yyyyzzz_1[j] + 2.0 * fl1_fxn * tg_zz_yyyzzz_1[j];

                    tg_yzz_yyyzzzz_0[j] = pb_y * tg_zz_yyyzzzz_0[j] + wp_y[j] * tg_zz_yyyzzzz_1[j] + 1.5 * fl1_fxn * tg_zz_yyzzzz_1[j];

                    tg_yzz_yyzzzzz_0[j] = pb_y * tg_zz_yyzzzzz_0[j] + wp_y[j] * tg_zz_yyzzzzz_1[j] + fl1_fxn * tg_zz_yzzzzz_1[j];

                    tg_yzz_yzzzzzz_0[j] = pb_y * tg_zz_yzzzzzz_0[j] + wp_y[j] * tg_zz_yzzzzzz_1[j] + 0.5 * fl1_fxn * tg_zz_zzzzzz_1[j];

                    tg_yzz_zzzzzzz_0[j] = pb_y * tg_zz_zzzzzzz_0[j] + wp_y[j] * tg_zz_zzzzzzz_1[j];

                    tg_zzz_xxxxxxx_0[j] = pb_z * tg_zz_xxxxxxx_0[j] + wp_z[j] * tg_zz_xxxxxxx_1[j] + fl1_fx * tg_z_xxxxxxx_0[j] - fl1_fx * fl1_fza * tg_z_xxxxxxx_1[j];

                    tg_zzz_xxxxxxy_0[j] = pb_z * tg_zz_xxxxxxy_0[j] + wp_z[j] * tg_zz_xxxxxxy_1[j] + fl1_fx * tg_z_xxxxxxy_0[j] - fl1_fx * fl1_fza * tg_z_xxxxxxy_1[j];

                    tg_zzz_xxxxxxz_0[j] = pb_z * tg_zz_xxxxxxz_0[j] + wp_z[j] * tg_zz_xxxxxxz_1[j] + fl1_fx * tg_z_xxxxxxz_0[j] - fl1_fx * fl1_fza * tg_z_xxxxxxz_1[j] + 0.5 * fl1_fxn * tg_zz_xxxxxx_1[j];

                    tg_zzz_xxxxxyy_0[j] = pb_z * tg_zz_xxxxxyy_0[j] + wp_z[j] * tg_zz_xxxxxyy_1[j] + fl1_fx * tg_z_xxxxxyy_0[j] - fl1_fx * fl1_fza * tg_z_xxxxxyy_1[j];

                    tg_zzz_xxxxxyz_0[j] = pb_z * tg_zz_xxxxxyz_0[j] + wp_z[j] * tg_zz_xxxxxyz_1[j] + fl1_fx * tg_z_xxxxxyz_0[j] - fl1_fx * fl1_fza * tg_z_xxxxxyz_1[j] + 0.5 * fl1_fxn * tg_zz_xxxxxy_1[j];

                    tg_zzz_xxxxxzz_0[j] = pb_z * tg_zz_xxxxxzz_0[j] + wp_z[j] * tg_zz_xxxxxzz_1[j] + fl1_fx * tg_z_xxxxxzz_0[j] - fl1_fx * fl1_fza * tg_z_xxxxxzz_1[j] + fl1_fxn * tg_zz_xxxxxz_1[j];

                    tg_zzz_xxxxyyy_0[j] = pb_z * tg_zz_xxxxyyy_0[j] + wp_z[j] * tg_zz_xxxxyyy_1[j] + fl1_fx * tg_z_xxxxyyy_0[j] - fl1_fx * fl1_fza * tg_z_xxxxyyy_1[j];

                    tg_zzz_xxxxyyz_0[j] = pb_z * tg_zz_xxxxyyz_0[j] + wp_z[j] * tg_zz_xxxxyyz_1[j] + fl1_fx * tg_z_xxxxyyz_0[j] - fl1_fx * fl1_fza * tg_z_xxxxyyz_1[j] + 0.5 * fl1_fxn * tg_zz_xxxxyy_1[j];

                    tg_zzz_xxxxyzz_0[j] = pb_z * tg_zz_xxxxyzz_0[j] + wp_z[j] * tg_zz_xxxxyzz_1[j] + fl1_fx * tg_z_xxxxyzz_0[j] - fl1_fx * fl1_fza * tg_z_xxxxyzz_1[j] + fl1_fxn * tg_zz_xxxxyz_1[j];

                    tg_zzz_xxxxzzz_0[j] = pb_z * tg_zz_xxxxzzz_0[j] + wp_z[j] * tg_zz_xxxxzzz_1[j] + fl1_fx * tg_z_xxxxzzz_0[j] - fl1_fx * fl1_fza * tg_z_xxxxzzz_1[j] + 1.5 * fl1_fxn * tg_zz_xxxxzz_1[j];

                    tg_zzz_xxxyyyy_0[j] = pb_z * tg_zz_xxxyyyy_0[j] + wp_z[j] * tg_zz_xxxyyyy_1[j] + fl1_fx * tg_z_xxxyyyy_0[j] - fl1_fx * fl1_fza * tg_z_xxxyyyy_1[j];

                    tg_zzz_xxxyyyz_0[j] = pb_z * tg_zz_xxxyyyz_0[j] + wp_z[j] * tg_zz_xxxyyyz_1[j] + fl1_fx * tg_z_xxxyyyz_0[j] - fl1_fx * fl1_fza * tg_z_xxxyyyz_1[j] + 0.5 * fl1_fxn * tg_zz_xxxyyy_1[j];

                    tg_zzz_xxxyyzz_0[j] = pb_z * tg_zz_xxxyyzz_0[j] + wp_z[j] * tg_zz_xxxyyzz_1[j] + fl1_fx * tg_z_xxxyyzz_0[j] - fl1_fx * fl1_fza * tg_z_xxxyyzz_1[j] + fl1_fxn * tg_zz_xxxyyz_1[j];

                    tg_zzz_xxxyzzz_0[j] = pb_z * tg_zz_xxxyzzz_0[j] + wp_z[j] * tg_zz_xxxyzzz_1[j] + fl1_fx * tg_z_xxxyzzz_0[j] - fl1_fx * fl1_fza * tg_z_xxxyzzz_1[j] + 1.5 * fl1_fxn * tg_zz_xxxyzz_1[j];

                    tg_zzz_xxxzzzz_0[j] = pb_z * tg_zz_xxxzzzz_0[j] + wp_z[j] * tg_zz_xxxzzzz_1[j] + fl1_fx * tg_z_xxxzzzz_0[j] - fl1_fx * fl1_fza * tg_z_xxxzzzz_1[j] + 2.0 * fl1_fxn * tg_zz_xxxzzz_1[j];

                    tg_zzz_xxyyyyy_0[j] = pb_z * tg_zz_xxyyyyy_0[j] + wp_z[j] * tg_zz_xxyyyyy_1[j] + fl1_fx * tg_z_xxyyyyy_0[j] - fl1_fx * fl1_fza * tg_z_xxyyyyy_1[j];

                    tg_zzz_xxyyyyz_0[j] = pb_z * tg_zz_xxyyyyz_0[j] + wp_z[j] * tg_zz_xxyyyyz_1[j] + fl1_fx * tg_z_xxyyyyz_0[j] - fl1_fx * fl1_fza * tg_z_xxyyyyz_1[j] + 0.5 * fl1_fxn * tg_zz_xxyyyy_1[j];

                    tg_zzz_xxyyyzz_0[j] = pb_z * tg_zz_xxyyyzz_0[j] + wp_z[j] * tg_zz_xxyyyzz_1[j] + fl1_fx * tg_z_xxyyyzz_0[j] - fl1_fx * fl1_fza * tg_z_xxyyyzz_1[j] + fl1_fxn * tg_zz_xxyyyz_1[j];

                    tg_zzz_xxyyzzz_0[j] = pb_z * tg_zz_xxyyzzz_0[j] + wp_z[j] * tg_zz_xxyyzzz_1[j] + fl1_fx * tg_z_xxyyzzz_0[j] - fl1_fx * fl1_fza * tg_z_xxyyzzz_1[j] + 1.5 * fl1_fxn * tg_zz_xxyyzz_1[j];

                    tg_zzz_xxyzzzz_0[j] = pb_z * tg_zz_xxyzzzz_0[j] + wp_z[j] * tg_zz_xxyzzzz_1[j] + fl1_fx * tg_z_xxyzzzz_0[j] - fl1_fx * fl1_fza * tg_z_xxyzzzz_1[j] + 2.0 * fl1_fxn * tg_zz_xxyzzz_1[j];

                    tg_zzz_xxzzzzz_0[j] = pb_z * tg_zz_xxzzzzz_0[j] + wp_z[j] * tg_zz_xxzzzzz_1[j] + fl1_fx * tg_z_xxzzzzz_0[j] - fl1_fx * fl1_fza * tg_z_xxzzzzz_1[j] + 2.5 * fl1_fxn * tg_zz_xxzzzz_1[j];

                    tg_zzz_xyyyyyy_0[j] = pb_z * tg_zz_xyyyyyy_0[j] + wp_z[j] * tg_zz_xyyyyyy_1[j] + fl1_fx * tg_z_xyyyyyy_0[j] - fl1_fx * fl1_fza * tg_z_xyyyyyy_1[j];

                    tg_zzz_xyyyyyz_0[j] = pb_z * tg_zz_xyyyyyz_0[j] + wp_z[j] * tg_zz_xyyyyyz_1[j] + fl1_fx * tg_z_xyyyyyz_0[j] - fl1_fx * fl1_fza * tg_z_xyyyyyz_1[j] + 0.5 * fl1_fxn * tg_zz_xyyyyy_1[j];

                    tg_zzz_xyyyyzz_0[j] = pb_z * tg_zz_xyyyyzz_0[j] + wp_z[j] * tg_zz_xyyyyzz_1[j] + fl1_fx * tg_z_xyyyyzz_0[j] - fl1_fx * fl1_fza * tg_z_xyyyyzz_1[j] + fl1_fxn * tg_zz_xyyyyz_1[j];

                    tg_zzz_xyyyzzz_0[j] = pb_z * tg_zz_xyyyzzz_0[j] + wp_z[j] * tg_zz_xyyyzzz_1[j] + fl1_fx * tg_z_xyyyzzz_0[j] - fl1_fx * fl1_fza * tg_z_xyyyzzz_1[j] + 1.5 * fl1_fxn * tg_zz_xyyyzz_1[j];

                    tg_zzz_xyyzzzz_0[j] = pb_z * tg_zz_xyyzzzz_0[j] + wp_z[j] * tg_zz_xyyzzzz_1[j] + fl1_fx * tg_z_xyyzzzz_0[j] - fl1_fx * fl1_fza * tg_z_xyyzzzz_1[j] + 2.0 * fl1_fxn * tg_zz_xyyzzz_1[j];

                    tg_zzz_xyzzzzz_0[j] = pb_z * tg_zz_xyzzzzz_0[j] + wp_z[j] * tg_zz_xyzzzzz_1[j] + fl1_fx * tg_z_xyzzzzz_0[j] - fl1_fx * fl1_fza * tg_z_xyzzzzz_1[j] + 2.5 * fl1_fxn * tg_zz_xyzzzz_1[j];

                    tg_zzz_xzzzzzz_0[j] = pb_z * tg_zz_xzzzzzz_0[j] + wp_z[j] * tg_zz_xzzzzzz_1[j] + fl1_fx * tg_z_xzzzzzz_0[j] - fl1_fx * fl1_fza * tg_z_xzzzzzz_1[j] + 3.0 * fl1_fxn * tg_zz_xzzzzz_1[j];

                    tg_zzz_yyyyyyy_0[j] = pb_z * tg_zz_yyyyyyy_0[j] + wp_z[j] * tg_zz_yyyyyyy_1[j] + fl1_fx * tg_z_yyyyyyy_0[j] - fl1_fx * fl1_fza * tg_z_yyyyyyy_1[j];

                    tg_zzz_yyyyyyz_0[j] = pb_z * tg_zz_yyyyyyz_0[j] + wp_z[j] * tg_zz_yyyyyyz_1[j] + fl1_fx * tg_z_yyyyyyz_0[j] - fl1_fx * fl1_fza * tg_z_yyyyyyz_1[j] + 0.5 * fl1_fxn * tg_zz_yyyyyy_1[j];

                    tg_zzz_yyyyyzz_0[j] = pb_z * tg_zz_yyyyyzz_0[j] + wp_z[j] * tg_zz_yyyyyzz_1[j] + fl1_fx * tg_z_yyyyyzz_0[j] - fl1_fx * fl1_fza * tg_z_yyyyyzz_1[j] + fl1_fxn * tg_zz_yyyyyz_1[j];

                    tg_zzz_yyyyzzz_0[j] = pb_z * tg_zz_yyyyzzz_0[j] + wp_z[j] * tg_zz_yyyyzzz_1[j] + fl1_fx * tg_z_yyyyzzz_0[j] - fl1_fx * fl1_fza * tg_z_yyyyzzz_1[j] + 1.5 * fl1_fxn * tg_zz_yyyyzz_1[j];

                    tg_zzz_yyyzzzz_0[j] = pb_z * tg_zz_yyyzzzz_0[j] + wp_z[j] * tg_zz_yyyzzzz_1[j] + fl1_fx * tg_z_yyyzzzz_0[j] - fl1_fx * fl1_fza * tg_z_yyyzzzz_1[j] + 2.0 * fl1_fxn * tg_zz_yyyzzz_1[j];

                    tg_zzz_yyzzzzz_0[j] = pb_z * tg_zz_yyzzzzz_0[j] + wp_z[j] * tg_zz_yyzzzzz_1[j] + fl1_fx * tg_z_yyzzzzz_0[j] - fl1_fx * fl1_fza * tg_z_yyzzzzz_1[j] + 2.5 * fl1_fxn * tg_zz_yyzzzz_1[j];

                    tg_zzz_yzzzzzz_0[j] = pb_z * tg_zz_yzzzzzz_0[j] + wp_z[j] * tg_zz_yzzzzzz_1[j] + fl1_fx * tg_z_yzzzzzz_0[j] - fl1_fx * fl1_fza * tg_z_yzzzzzz_1[j] + 3.0 * fl1_fxn * tg_zz_yzzzzz_1[j];

                    tg_zzz_zzzzzzz_0[j] = pb_z * tg_zz_zzzzzzz_0[j] + wp_z[j] * tg_zz_zzzzzzz_1[j] + fl1_fx * tg_z_zzzzzzz_0[j] - fl1_fx * fl1_fza * tg_z_zzzzzzz_1[j] + 3.5 * fl1_fxn * tg_zz_zzzzzz_1[j];
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

