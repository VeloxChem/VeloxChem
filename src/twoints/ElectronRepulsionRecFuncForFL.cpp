//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectronRepulsionRecFuncForFL.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSFSL(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSFSL_0_90(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSFSL_90_180(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 

        erirecfunc::compElectronRepulsionForSFSL_180_270(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSFSL_270_360(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSFSL_360_450(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 
    }

    void
    compElectronRepulsionForSFSL_0_90(      CMemBlock2D<double>& primBuffer,
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
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_3_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_3_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_2_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_2_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_1_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_1_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_2_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {7, -1, -1, -1}, 
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

                auto tg_xx_xxxxxxxx_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx); 

                auto tg_xx_xxxxxxxy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 1); 

                auto tg_xx_xxxxxxxz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 2); 

                auto tg_xx_xxxxxxyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 3); 

                auto tg_xx_xxxxxxyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 4); 

                auto tg_xx_xxxxxxzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 5); 

                auto tg_xx_xxxxxyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 6); 

                auto tg_xx_xxxxxyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 7); 

                auto tg_xx_xxxxxyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 8); 

                auto tg_xx_xxxxxzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 9); 

                auto tg_xx_xxxxyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 10); 

                auto tg_xx_xxxxyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 11); 

                auto tg_xx_xxxxyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 12); 

                auto tg_xx_xxxxyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 13); 

                auto tg_xx_xxxxzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 14); 

                auto tg_xx_xxxyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 15); 

                auto tg_xx_xxxyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 16); 

                auto tg_xx_xxxyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 17); 

                auto tg_xx_xxxyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 18); 

                auto tg_xx_xxxyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 19); 

                auto tg_xx_xxxzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 20); 

                auto tg_xx_xxyyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 21); 

                auto tg_xx_xxyyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 22); 

                auto tg_xx_xxyyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 23); 

                auto tg_xx_xxyyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 24); 

                auto tg_xx_xxyyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 25); 

                auto tg_xx_xxyzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 26); 

                auto tg_xx_xxzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 27); 

                auto tg_xx_xyyyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 28); 

                auto tg_xx_xyyyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 29); 

                auto tg_xx_xyyyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 30); 

                auto tg_xx_xyyyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 31); 

                auto tg_xx_xyyyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 32); 

                auto tg_xx_xyyzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 33); 

                auto tg_xx_xyzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 34); 

                auto tg_xx_xzzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 35); 

                auto tg_xx_yyyyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 36); 

                auto tg_xx_yyyyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 37); 

                auto tg_xx_yyyyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 38); 

                auto tg_xx_yyyyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 39); 

                auto tg_xx_yyyyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 40); 

                auto tg_xx_yyyzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 41); 

                auto tg_xx_yyzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 42); 

                auto tg_xx_yzzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 43); 

                auto tg_xx_zzzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 44); 

                auto tg_xy_xxxxxxxx_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 45); 

                auto tg_xy_xxxxxxxy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 46); 

                auto tg_xy_xxxxxxxz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 47); 

                auto tg_xy_xxxxxxyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 48); 

                auto tg_xy_xxxxxxyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 49); 

                auto tg_xy_xxxxxxzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 50); 

                auto tg_xy_xxxxxyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 51); 

                auto tg_xy_xxxxxyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 52); 

                auto tg_xy_xxxxxyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 53); 

                auto tg_xy_xxxxxzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 54); 

                auto tg_xy_xxxxyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 55); 

                auto tg_xy_xxxxyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 56); 

                auto tg_xy_xxxxyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 57); 

                auto tg_xy_xxxxyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 58); 

                auto tg_xy_xxxxzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 59); 

                auto tg_xy_xxxyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 60); 

                auto tg_xy_xxxyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 61); 

                auto tg_xy_xxxyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 62); 

                auto tg_xy_xxxyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 63); 

                auto tg_xy_xxxyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 64); 

                auto tg_xy_xxxzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 65); 

                auto tg_xy_xxyyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 66); 

                auto tg_xy_xxyyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 67); 

                auto tg_xy_xxyyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 68); 

                auto tg_xy_xxyyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 69); 

                auto tg_xy_xxyyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 70); 

                auto tg_xy_xxyzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 71); 

                auto tg_xy_xxzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 72); 

                auto tg_xy_xyyyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 73); 

                auto tg_xy_xyyyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 74); 

                auto tg_xy_xyyyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 75); 

                auto tg_xy_xyyyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 76); 

                auto tg_xy_xyyyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 77); 

                auto tg_xy_xyyzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 78); 

                auto tg_xy_xyzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 79); 

                auto tg_xy_xzzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 80); 

                auto tg_xy_yyyyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 81); 

                auto tg_xy_yyyyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 82); 

                auto tg_xy_yyyyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 83); 

                auto tg_xy_yyyyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 84); 

                auto tg_xy_yyyyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 85); 

                auto tg_xy_yyyzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 86); 

                auto tg_xy_yyzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 87); 

                auto tg_xy_yzzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 88); 

                auto tg_xy_zzzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 89); 

                auto tg_xx_xxxxxxxx_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx); 

                auto tg_xx_xxxxxxxy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 1); 

                auto tg_xx_xxxxxxxz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 2); 

                auto tg_xx_xxxxxxyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 3); 

                auto tg_xx_xxxxxxyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 4); 

                auto tg_xx_xxxxxxzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 5); 

                auto tg_xx_xxxxxyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 6); 

                auto tg_xx_xxxxxyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 7); 

                auto tg_xx_xxxxxyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 8); 

                auto tg_xx_xxxxxzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 9); 

                auto tg_xx_xxxxyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 10); 

                auto tg_xx_xxxxyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 11); 

                auto tg_xx_xxxxyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 12); 

                auto tg_xx_xxxxyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 13); 

                auto tg_xx_xxxxzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 14); 

                auto tg_xx_xxxyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 15); 

                auto tg_xx_xxxyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 16); 

                auto tg_xx_xxxyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 17); 

                auto tg_xx_xxxyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 18); 

                auto tg_xx_xxxyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 19); 

                auto tg_xx_xxxzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 20); 

                auto tg_xx_xxyyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 21); 

                auto tg_xx_xxyyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 22); 

                auto tg_xx_xxyyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 23); 

                auto tg_xx_xxyyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 24); 

                auto tg_xx_xxyyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 25); 

                auto tg_xx_xxyzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 26); 

                auto tg_xx_xxzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 27); 

                auto tg_xx_xyyyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 28); 

                auto tg_xx_xyyyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 29); 

                auto tg_xx_xyyyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 30); 

                auto tg_xx_xyyyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 31); 

                auto tg_xx_xyyyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 32); 

                auto tg_xx_xyyzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 33); 

                auto tg_xx_xyzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 34); 

                auto tg_xx_xzzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 35); 

                auto tg_xx_yyyyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 36); 

                auto tg_xx_yyyyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 37); 

                auto tg_xx_yyyyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 38); 

                auto tg_xx_yyyyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 39); 

                auto tg_xx_yyyyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 40); 

                auto tg_xx_yyyzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 41); 

                auto tg_xx_yyzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 42); 

                auto tg_xx_yzzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 43); 

                auto tg_xx_zzzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 44); 

                auto tg_xy_xxxxxxxx_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 45); 

                auto tg_xy_xxxxxxxy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 46); 

                auto tg_xy_xxxxxxxz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 47); 

                auto tg_xy_xxxxxxyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 48); 

                auto tg_xy_xxxxxxyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 49); 

                auto tg_xy_xxxxxxzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 50); 

                auto tg_xy_xxxxxyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 51); 

                auto tg_xy_xxxxxyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 52); 

                auto tg_xy_xxxxxyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 53); 

                auto tg_xy_xxxxxzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 54); 

                auto tg_xy_xxxxyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 55); 

                auto tg_xy_xxxxyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 56); 

                auto tg_xy_xxxxyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 57); 

                auto tg_xy_xxxxyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 58); 

                auto tg_xy_xxxxzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 59); 

                auto tg_xy_xxxyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 60); 

                auto tg_xy_xxxyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 61); 

                auto tg_xy_xxxyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 62); 

                auto tg_xy_xxxyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 63); 

                auto tg_xy_xxxyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 64); 

                auto tg_xy_xxxzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 65); 

                auto tg_xy_xxyyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 66); 

                auto tg_xy_xxyyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 67); 

                auto tg_xy_xxyyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 68); 

                auto tg_xy_xxyyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 69); 

                auto tg_xy_xxyyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 70); 

                auto tg_xy_xxyzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 71); 

                auto tg_xy_xxzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 72); 

                auto tg_xy_xyyyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 73); 

                auto tg_xy_xyyyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 74); 

                auto tg_xy_xyyyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 75); 

                auto tg_xy_xyyyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 76); 

                auto tg_xy_xyyyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 77); 

                auto tg_xy_xyyzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 78); 

                auto tg_xy_xyzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 79); 

                auto tg_xy_xzzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 80); 

                auto tg_xy_yyyyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 81); 

                auto tg_xy_yyyyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 82); 

                auto tg_xy_yyyyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 83); 

                auto tg_xy_yyyyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 84); 

                auto tg_xy_yyyyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 85); 

                auto tg_xy_yyyzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 86); 

                auto tg_xy_yyzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 87); 

                auto tg_xy_yzzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 88); 

                auto tg_xy_zzzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 89); 

                auto tg_x_xxxxxxxx_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx); 

                auto tg_x_xxxxxxxy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 1); 

                auto tg_x_xxxxxxxz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 2); 

                auto tg_x_xxxxxxyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 3); 

                auto tg_x_xxxxxxyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 4); 

                auto tg_x_xxxxxxzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 5); 

                auto tg_x_xxxxxyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 6); 

                auto tg_x_xxxxxyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 7); 

                auto tg_x_xxxxxyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 8); 

                auto tg_x_xxxxxzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 9); 

                auto tg_x_xxxxyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 10); 

                auto tg_x_xxxxyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 11); 

                auto tg_x_xxxxyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 12); 

                auto tg_x_xxxxyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 13); 

                auto tg_x_xxxxzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 14); 

                auto tg_x_xxxyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 15); 

                auto tg_x_xxxyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 16); 

                auto tg_x_xxxyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 17); 

                auto tg_x_xxxyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 18); 

                auto tg_x_xxxyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 19); 

                auto tg_x_xxxzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 20); 

                auto tg_x_xxyyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 21); 

                auto tg_x_xxyyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 22); 

                auto tg_x_xxyyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 23); 

                auto tg_x_xxyyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 24); 

                auto tg_x_xxyyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 25); 

                auto tg_x_xxyzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 26); 

                auto tg_x_xxzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 27); 

                auto tg_x_xyyyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 28); 

                auto tg_x_xyyyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 29); 

                auto tg_x_xyyyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 30); 

                auto tg_x_xyyyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 31); 

                auto tg_x_xyyyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 32); 

                auto tg_x_xyyzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 33); 

                auto tg_x_xyzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 34); 

                auto tg_x_xzzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 35); 

                auto tg_x_yyyyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 36); 

                auto tg_x_yyyyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 37); 

                auto tg_x_yyyyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 38); 

                auto tg_x_yyyyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 39); 

                auto tg_x_yyyyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 40); 

                auto tg_x_yyyzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 41); 

                auto tg_x_yyzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 42); 

                auto tg_x_yzzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 43); 

                auto tg_x_zzzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 44); 

                auto tg_y_xxxxxxxx_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 45); 

                auto tg_y_xxxxxxxy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 46); 

                auto tg_y_xxxxxxxz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 47); 

                auto tg_y_xxxxxxyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 48); 

                auto tg_y_xxxxxxyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 49); 

                auto tg_y_xxxxxxzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 50); 

                auto tg_y_xxxxxyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 51); 

                auto tg_y_xxxxxyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 52); 

                auto tg_y_xxxxxyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 53); 

                auto tg_y_xxxxxzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 54); 

                auto tg_y_xxxxyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 55); 

                auto tg_y_xxxxyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 56); 

                auto tg_y_xxxxyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 57); 

                auto tg_y_xxxxyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 58); 

                auto tg_y_xxxxzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 59); 

                auto tg_y_xxxyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 60); 

                auto tg_y_xxxyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 61); 

                auto tg_y_xxxyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 62); 

                auto tg_y_xxxyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 63); 

                auto tg_y_xxxyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 64); 

                auto tg_y_xxxzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 65); 

                auto tg_y_xxyyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 66); 

                auto tg_y_xxyyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 67); 

                auto tg_y_xxyyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 68); 

                auto tg_y_xxyyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 69); 

                auto tg_y_xxyyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 70); 

                auto tg_y_xxyzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 71); 

                auto tg_y_xxzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 72); 

                auto tg_y_xyyyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 73); 

                auto tg_y_xyyyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 74); 

                auto tg_y_xyyyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 75); 

                auto tg_y_xyyyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 76); 

                auto tg_y_xyyyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 77); 

                auto tg_y_xyyzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 78); 

                auto tg_y_xyzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 79); 

                auto tg_y_xzzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 80); 

                auto tg_y_yyyyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 81); 

                auto tg_y_yyyyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 82); 

                auto tg_y_yyyyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 83); 

                auto tg_y_yyyyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 84); 

                auto tg_y_yyyyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 85); 

                auto tg_y_yyyzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 86); 

                auto tg_y_yyzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 87); 

                auto tg_y_yzzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 88); 

                auto tg_y_zzzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 89); 

                auto tg_x_xxxxxxxx_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx); 

                auto tg_x_xxxxxxxy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 1); 

                auto tg_x_xxxxxxxz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 2); 

                auto tg_x_xxxxxxyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 3); 

                auto tg_x_xxxxxxyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 4); 

                auto tg_x_xxxxxxzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 5); 

                auto tg_x_xxxxxyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 6); 

                auto tg_x_xxxxxyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 7); 

                auto tg_x_xxxxxyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 8); 

                auto tg_x_xxxxxzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 9); 

                auto tg_x_xxxxyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 10); 

                auto tg_x_xxxxyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 11); 

                auto tg_x_xxxxyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 12); 

                auto tg_x_xxxxyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 13); 

                auto tg_x_xxxxzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 14); 

                auto tg_x_xxxyyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 15); 

                auto tg_x_xxxyyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 16); 

                auto tg_x_xxxyyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 17); 

                auto tg_x_xxxyyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 18); 

                auto tg_x_xxxyzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 19); 

                auto tg_x_xxxzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 20); 

                auto tg_x_xxyyyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 21); 

                auto tg_x_xxyyyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 22); 

                auto tg_x_xxyyyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 23); 

                auto tg_x_xxyyyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 24); 

                auto tg_x_xxyyzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 25); 

                auto tg_x_xxyzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 26); 

                auto tg_x_xxzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 27); 

                auto tg_x_xyyyyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 28); 

                auto tg_x_xyyyyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 29); 

                auto tg_x_xyyyyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 30); 

                auto tg_x_xyyyyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 31); 

                auto tg_x_xyyyzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 32); 

                auto tg_x_xyyzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 33); 

                auto tg_x_xyzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 34); 

                auto tg_x_xzzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 35); 

                auto tg_x_yyyyyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 36); 

                auto tg_x_yyyyyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 37); 

                auto tg_x_yyyyyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 38); 

                auto tg_x_yyyyyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 39); 

                auto tg_x_yyyyzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 40); 

                auto tg_x_yyyzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 41); 

                auto tg_x_yyzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 42); 

                auto tg_x_yzzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 43); 

                auto tg_x_zzzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 44); 

                auto tg_y_xxxxxxxx_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 45); 

                auto tg_y_xxxxxxxy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 46); 

                auto tg_y_xxxxxxxz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 47); 

                auto tg_y_xxxxxxyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 48); 

                auto tg_y_xxxxxxyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 49); 

                auto tg_y_xxxxxxzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 50); 

                auto tg_y_xxxxxyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 51); 

                auto tg_y_xxxxxyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 52); 

                auto tg_y_xxxxxyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 53); 

                auto tg_y_xxxxxzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 54); 

                auto tg_y_xxxxyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 55); 

                auto tg_y_xxxxyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 56); 

                auto tg_y_xxxxyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 57); 

                auto tg_y_xxxxyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 58); 

                auto tg_y_xxxxzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 59); 

                auto tg_y_xxxyyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 60); 

                auto tg_y_xxxyyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 61); 

                auto tg_y_xxxyyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 62); 

                auto tg_y_xxxyyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 63); 

                auto tg_y_xxxyzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 64); 

                auto tg_y_xxxzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 65); 

                auto tg_y_xxyyyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 66); 

                auto tg_y_xxyyyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 67); 

                auto tg_y_xxyyyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 68); 

                auto tg_y_xxyyyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 69); 

                auto tg_y_xxyyzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 70); 

                auto tg_y_xxyzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 71); 

                auto tg_y_xxzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 72); 

                auto tg_y_xyyyyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 73); 

                auto tg_y_xyyyyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 74); 

                auto tg_y_xyyyyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 75); 

                auto tg_y_xyyyyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 76); 

                auto tg_y_xyyyzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 77); 

                auto tg_y_xyyzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 78); 

                auto tg_y_xyzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 79); 

                auto tg_y_xzzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 80); 

                auto tg_y_yyyyyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 81); 

                auto tg_y_yyyyyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 82); 

                auto tg_y_yyyyyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 83); 

                auto tg_y_yyyyyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 84); 

                auto tg_y_yyyyzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 85); 

                auto tg_y_yyyzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 86); 

                auto tg_y_yyzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 87); 

                auto tg_y_yzzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 88); 

                auto tg_y_zzzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 89); 

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

                // set up pointers to integrals

                auto tg_xxx_xxxxxxxx_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx); 

                auto tg_xxx_xxxxxxxy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 1); 

                auto tg_xxx_xxxxxxxz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 2); 

                auto tg_xxx_xxxxxxyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 3); 

                auto tg_xxx_xxxxxxyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 4); 

                auto tg_xxx_xxxxxxzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 5); 

                auto tg_xxx_xxxxxyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 6); 

                auto tg_xxx_xxxxxyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 7); 

                auto tg_xxx_xxxxxyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 8); 

                auto tg_xxx_xxxxxzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 9); 

                auto tg_xxx_xxxxyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 10); 

                auto tg_xxx_xxxxyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 11); 

                auto tg_xxx_xxxxyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 12); 

                auto tg_xxx_xxxxyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 13); 

                auto tg_xxx_xxxxzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 14); 

                auto tg_xxx_xxxyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 15); 

                auto tg_xxx_xxxyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 16); 

                auto tg_xxx_xxxyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 17); 

                auto tg_xxx_xxxyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 18); 

                auto tg_xxx_xxxyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 19); 

                auto tg_xxx_xxxzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 20); 

                auto tg_xxx_xxyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 21); 

                auto tg_xxx_xxyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 22); 

                auto tg_xxx_xxyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 23); 

                auto tg_xxx_xxyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 24); 

                auto tg_xxx_xxyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 25); 

                auto tg_xxx_xxyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 26); 

                auto tg_xxx_xxzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 27); 

                auto tg_xxx_xyyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 28); 

                auto tg_xxx_xyyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 29); 

                auto tg_xxx_xyyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 30); 

                auto tg_xxx_xyyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 31); 

                auto tg_xxx_xyyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 32); 

                auto tg_xxx_xyyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 33); 

                auto tg_xxx_xyzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 34); 

                auto tg_xxx_xzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 35); 

                auto tg_xxx_yyyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 36); 

                auto tg_xxx_yyyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 37); 

                auto tg_xxx_yyyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 38); 

                auto tg_xxx_yyyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 39); 

                auto tg_xxx_yyyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 40); 

                auto tg_xxx_yyyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 41); 

                auto tg_xxx_yyzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 42); 

                auto tg_xxx_yzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 43); 

                auto tg_xxx_zzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 44); 

                auto tg_xxy_xxxxxxxx_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 45); 

                auto tg_xxy_xxxxxxxy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 46); 

                auto tg_xxy_xxxxxxxz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 47); 

                auto tg_xxy_xxxxxxyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 48); 

                auto tg_xxy_xxxxxxyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 49); 

                auto tg_xxy_xxxxxxzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 50); 

                auto tg_xxy_xxxxxyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 51); 

                auto tg_xxy_xxxxxyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 52); 

                auto tg_xxy_xxxxxyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 53); 

                auto tg_xxy_xxxxxzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 54); 

                auto tg_xxy_xxxxyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 55); 

                auto tg_xxy_xxxxyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 56); 

                auto tg_xxy_xxxxyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 57); 

                auto tg_xxy_xxxxyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 58); 

                auto tg_xxy_xxxxzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 59); 

                auto tg_xxy_xxxyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 60); 

                auto tg_xxy_xxxyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 61); 

                auto tg_xxy_xxxyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 62); 

                auto tg_xxy_xxxyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 63); 

                auto tg_xxy_xxxyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 64); 

                auto tg_xxy_xxxzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 65); 

                auto tg_xxy_xxyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 66); 

                auto tg_xxy_xxyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 67); 

                auto tg_xxy_xxyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 68); 

                auto tg_xxy_xxyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 69); 

                auto tg_xxy_xxyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 70); 

                auto tg_xxy_xxyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 71); 

                auto tg_xxy_xxzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 72); 

                auto tg_xxy_xyyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 73); 

                auto tg_xxy_xyyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 74); 

                auto tg_xxy_xyyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 75); 

                auto tg_xxy_xyyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 76); 

                auto tg_xxy_xyyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 77); 

                auto tg_xxy_xyyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 78); 

                auto tg_xxy_xyzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 79); 

                auto tg_xxy_xzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 80); 

                auto tg_xxy_yyyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 81); 

                auto tg_xxy_yyyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 82); 

                auto tg_xxy_yyyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 83); 

                auto tg_xxy_yyyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 84); 

                auto tg_xxy_yyyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 85); 

                auto tg_xxy_yyyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 86); 

                auto tg_xxy_yyzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 87); 

                auto tg_xxy_yzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 88); 

                auto tg_xxy_zzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 89); 

                // Batch of Integrals (0,90)

                #pragma omp simd aligned(fxn, fza, tg_x_xxxxxxxx_0, tg_x_xxxxxxxx_1, tg_x_xxxxxxxy_0, \
                                         tg_x_xxxxxxxy_1, tg_x_xxxxxxxz_0, tg_x_xxxxxxxz_1, tg_x_xxxxxxyy_0, tg_x_xxxxxxyy_1, \
                                         tg_x_xxxxxxyz_0, tg_x_xxxxxxyz_1, tg_x_xxxxxxzz_0, tg_x_xxxxxxzz_1, tg_x_xxxxxyyy_0, \
                                         tg_x_xxxxxyyy_1, tg_x_xxxxxyyz_0, tg_x_xxxxxyyz_1, tg_x_xxxxxyzz_0, tg_x_xxxxxyzz_1, \
                                         tg_x_xxxxxzzz_0, tg_x_xxxxxzzz_1, tg_x_xxxxyyyy_0, tg_x_xxxxyyyy_1, tg_x_xxxxyyyz_0, \
                                         tg_x_xxxxyyyz_1, tg_x_xxxxyyzz_0, tg_x_xxxxyyzz_1, tg_x_xxxxyzzz_0, tg_x_xxxxyzzz_1, \
                                         tg_x_xxxxzzzz_0, tg_x_xxxxzzzz_1, tg_x_xxxyyyyy_0, tg_x_xxxyyyyy_1, tg_x_xxxyyyyz_0, \
                                         tg_x_xxxyyyyz_1, tg_x_xxxyyyzz_0, tg_x_xxxyyyzz_1, tg_x_xxxyyzzz_0, tg_x_xxxyyzzz_1, \
                                         tg_x_xxxyzzzz_0, tg_x_xxxyzzzz_1, tg_x_xxxzzzzz_0, tg_x_xxxzzzzz_1, tg_x_xxyyyyyy_0, \
                                         tg_x_xxyyyyyy_1, tg_x_xxyyyyyz_0, tg_x_xxyyyyyz_1, tg_x_xxyyyyzz_0, tg_x_xxyyyyzz_1, \
                                         tg_x_xxyyyzzz_0, tg_x_xxyyyzzz_1, tg_x_xxyyzzzz_0, tg_x_xxyyzzzz_1, tg_x_xxyzzzzz_0, \
                                         tg_x_xxyzzzzz_1, tg_x_xxzzzzzz_0, tg_x_xxzzzzzz_1, tg_x_xyyyyyyy_0, tg_x_xyyyyyyy_1, \
                                         tg_x_xyyyyyyz_0, tg_x_xyyyyyyz_1, tg_x_xyyyyyzz_0, tg_x_xyyyyyzz_1, tg_x_xyyyyzzz_0, \
                                         tg_x_xyyyyzzz_1, tg_x_xyyyzzzz_0, tg_x_xyyyzzzz_1, tg_x_xyyzzzzz_0, tg_x_xyyzzzzz_1, \
                                         tg_x_xyzzzzzz_0, tg_x_xyzzzzzz_1, tg_x_xzzzzzzz_0, tg_x_xzzzzzzz_1, tg_x_yyyyyyyy_0, \
                                         tg_x_yyyyyyyy_1, tg_x_yyyyyyyz_0, tg_x_yyyyyyyz_1, tg_x_yyyyyyzz_0, tg_x_yyyyyyzz_1, \
                                         tg_x_yyyyyzzz_0, tg_x_yyyyyzzz_1, tg_x_yyyyzzzz_0, tg_x_yyyyzzzz_1, tg_x_yyyzzzzz_0, \
                                         tg_x_yyyzzzzz_1, tg_x_yyzzzzzz_0, tg_x_yyzzzzzz_1, tg_x_yzzzzzzz_0, tg_x_yzzzzzzz_1, \
                                         tg_x_zzzzzzzz_0, tg_x_zzzzzzzz_1, tg_xx_xxxxxxx_1, tg_xx_xxxxxxxx_0, \
                                         tg_xx_xxxxxxxx_1, tg_xx_xxxxxxxy_0, tg_xx_xxxxxxxy_1, tg_xx_xxxxxxxz_0, \
                                         tg_xx_xxxxxxxz_1, tg_xx_xxxxxxy_1, tg_xx_xxxxxxyy_0, tg_xx_xxxxxxyy_1, \
                                         tg_xx_xxxxxxyz_0, tg_xx_xxxxxxyz_1, tg_xx_xxxxxxz_1, tg_xx_xxxxxxzz_0, \
                                         tg_xx_xxxxxxzz_1, tg_xx_xxxxxyy_1, tg_xx_xxxxxyyy_0, tg_xx_xxxxxyyy_1, \
                                         tg_xx_xxxxxyyz_0, tg_xx_xxxxxyyz_1, tg_xx_xxxxxyz_1, tg_xx_xxxxxyzz_0, \
                                         tg_xx_xxxxxyzz_1, tg_xx_xxxxxzz_1, tg_xx_xxxxxzzz_0, tg_xx_xxxxxzzz_1, \
                                         tg_xx_xxxxyyy_1, tg_xx_xxxxyyyy_0, tg_xx_xxxxyyyy_1, tg_xx_xxxxyyyz_0, \
                                         tg_xx_xxxxyyyz_1, tg_xx_xxxxyyz_1, tg_xx_xxxxyyzz_0, tg_xx_xxxxyyzz_1, \
                                         tg_xx_xxxxyzz_1, tg_xx_xxxxyzzz_0, tg_xx_xxxxyzzz_1, tg_xx_xxxxzzz_1, \
                                         tg_xx_xxxxzzzz_0, tg_xx_xxxxzzzz_1, tg_xx_xxxyyyy_1, tg_xx_xxxyyyyy_0, \
                                         tg_xx_xxxyyyyy_1, tg_xx_xxxyyyyz_0, tg_xx_xxxyyyyz_1, tg_xx_xxxyyyz_1, \
                                         tg_xx_xxxyyyzz_0, tg_xx_xxxyyyzz_1, tg_xx_xxxyyzz_1, tg_xx_xxxyyzzz_0, \
                                         tg_xx_xxxyyzzz_1, tg_xx_xxxyzzz_1, tg_xx_xxxyzzzz_0, tg_xx_xxxyzzzz_1, \
                                         tg_xx_xxxzzzz_1, tg_xx_xxxzzzzz_0, tg_xx_xxxzzzzz_1, tg_xx_xxyyyyy_1, \
                                         tg_xx_xxyyyyyy_0, tg_xx_xxyyyyyy_1, tg_xx_xxyyyyyz_0, tg_xx_xxyyyyyz_1, \
                                         tg_xx_xxyyyyz_1, tg_xx_xxyyyyzz_0, tg_xx_xxyyyyzz_1, tg_xx_xxyyyzz_1, \
                                         tg_xx_xxyyyzzz_0, tg_xx_xxyyyzzz_1, tg_xx_xxyyzzz_1, tg_xx_xxyyzzzz_0, \
                                         tg_xx_xxyyzzzz_1, tg_xx_xxyzzzz_1, tg_xx_xxyzzzzz_0, tg_xx_xxyzzzzz_1, \
                                         tg_xx_xxzzzzz_1, tg_xx_xxzzzzzz_0, tg_xx_xxzzzzzz_1, tg_xx_xyyyyyy_1, \
                                         tg_xx_xyyyyyyy_0, tg_xx_xyyyyyyy_1, tg_xx_xyyyyyyz_0, tg_xx_xyyyyyyz_1, \
                                         tg_xx_xyyyyyz_1, tg_xx_xyyyyyzz_0, tg_xx_xyyyyyzz_1, tg_xx_xyyyyzz_1, \
                                         tg_xx_xyyyyzzz_0, tg_xx_xyyyyzzz_1, tg_xx_xyyyzzz_1, tg_xx_xyyyzzzz_0, \
                                         tg_xx_xyyyzzzz_1, tg_xx_xyyzzzz_1, tg_xx_xyyzzzzz_0, tg_xx_xyyzzzzz_1, \
                                         tg_xx_xyzzzzz_1, tg_xx_xyzzzzzz_0, tg_xx_xyzzzzzz_1, tg_xx_xzzzzzz_1, \
                                         tg_xx_xzzzzzzz_0, tg_xx_xzzzzzzz_1, tg_xx_yyyyyyy_1, tg_xx_yyyyyyyy_0, \
                                         tg_xx_yyyyyyyy_1, tg_xx_yyyyyyyz_0, tg_xx_yyyyyyyz_1, tg_xx_yyyyyyz_1, \
                                         tg_xx_yyyyyyzz_0, tg_xx_yyyyyyzz_1, tg_xx_yyyyyzz_1, tg_xx_yyyyyzzz_0, \
                                         tg_xx_yyyyyzzz_1, tg_xx_yyyyzzz_1, tg_xx_yyyyzzzz_0, tg_xx_yyyyzzzz_1, \
                                         tg_xx_yyyzzzz_1, tg_xx_yyyzzzzz_0, tg_xx_yyyzzzzz_1, tg_xx_yyzzzzz_1, \
                                         tg_xx_yyzzzzzz_0, tg_xx_yyzzzzzz_1, tg_xx_yzzzzzz_1, tg_xx_yzzzzzzz_0, \
                                         tg_xx_yzzzzzzz_1, tg_xx_zzzzzzz_1, tg_xx_zzzzzzzz_0, tg_xx_zzzzzzzz_1, \
                                         tg_xxx_xxxxxxxx_0, tg_xxx_xxxxxxxy_0, tg_xxx_xxxxxxxz_0, tg_xxx_xxxxxxyy_0, \
                                         tg_xxx_xxxxxxyz_0, tg_xxx_xxxxxxzz_0, tg_xxx_xxxxxyyy_0, tg_xxx_xxxxxyyz_0, \
                                         tg_xxx_xxxxxyzz_0, tg_xxx_xxxxxzzz_0, tg_xxx_xxxxyyyy_0, tg_xxx_xxxxyyyz_0, \
                                         tg_xxx_xxxxyyzz_0, tg_xxx_xxxxyzzz_0, tg_xxx_xxxxzzzz_0, tg_xxx_xxxyyyyy_0, \
                                         tg_xxx_xxxyyyyz_0, tg_xxx_xxxyyyzz_0, tg_xxx_xxxyyzzz_0, tg_xxx_xxxyzzzz_0, \
                                         tg_xxx_xxxzzzzz_0, tg_xxx_xxyyyyyy_0, tg_xxx_xxyyyyyz_0, tg_xxx_xxyyyyzz_0, \
                                         tg_xxx_xxyyyzzz_0, tg_xxx_xxyyzzzz_0, tg_xxx_xxyzzzzz_0, tg_xxx_xxzzzzzz_0, \
                                         tg_xxx_xyyyyyyy_0, tg_xxx_xyyyyyyz_0, tg_xxx_xyyyyyzz_0, tg_xxx_xyyyyzzz_0, \
                                         tg_xxx_xyyyzzzz_0, tg_xxx_xyyzzzzz_0, tg_xxx_xyzzzzzz_0, tg_xxx_xzzzzzzz_0, \
                                         tg_xxx_yyyyyyyy_0, tg_xxx_yyyyyyyz_0, tg_xxx_yyyyyyzz_0, tg_xxx_yyyyyzzz_0, \
                                         tg_xxx_yyyyzzzz_0, tg_xxx_yyyzzzzz_0, tg_xxx_yyzzzzzz_0, tg_xxx_yzzzzzzz_0, \
                                         tg_xxx_zzzzzzzz_0, tg_xxy_xxxxxxxx_0, tg_xxy_xxxxxxxy_0, tg_xxy_xxxxxxxz_0, \
                                         tg_xxy_xxxxxxyy_0, tg_xxy_xxxxxxyz_0, tg_xxy_xxxxxxzz_0, tg_xxy_xxxxxyyy_0, \
                                         tg_xxy_xxxxxyyz_0, tg_xxy_xxxxxyzz_0, tg_xxy_xxxxxzzz_0, tg_xxy_xxxxyyyy_0, \
                                         tg_xxy_xxxxyyyz_0, tg_xxy_xxxxyyzz_0, tg_xxy_xxxxyzzz_0, tg_xxy_xxxxzzzz_0, \
                                         tg_xxy_xxxyyyyy_0, tg_xxy_xxxyyyyz_0, tg_xxy_xxxyyyzz_0, tg_xxy_xxxyyzzz_0, \
                                         tg_xxy_xxxyzzzz_0, tg_xxy_xxxzzzzz_0, tg_xxy_xxyyyyyy_0, tg_xxy_xxyyyyyz_0, \
                                         tg_xxy_xxyyyyzz_0, tg_xxy_xxyyyzzz_0, tg_xxy_xxyyzzzz_0, tg_xxy_xxyzzzzz_0, \
                                         tg_xxy_xxzzzzzz_0, tg_xxy_xyyyyyyy_0, tg_xxy_xyyyyyyz_0, tg_xxy_xyyyyyzz_0, \
                                         tg_xxy_xyyyyzzz_0, tg_xxy_xyyyzzzz_0, tg_xxy_xyyzzzzz_0, tg_xxy_xyzzzzzz_0, \
                                         tg_xxy_xzzzzzzz_0, tg_xxy_yyyyyyyy_0, tg_xxy_yyyyyyyz_0, tg_xxy_yyyyyyzz_0, \
                                         tg_xxy_yyyyyzzz_0, tg_xxy_yyyyzzzz_0, tg_xxy_yyyzzzzz_0, tg_xxy_yyzzzzzz_0, \
                                         tg_xxy_yzzzzzzz_0, tg_xxy_zzzzzzzz_0, tg_xy_xxxxxxx_1, tg_xy_xxxxxxxx_0, \
                                         tg_xy_xxxxxxxx_1, tg_xy_xxxxxxxy_0, tg_xy_xxxxxxxy_1, tg_xy_xxxxxxxz_0, \
                                         tg_xy_xxxxxxxz_1, tg_xy_xxxxxxy_1, tg_xy_xxxxxxyy_0, tg_xy_xxxxxxyy_1, \
                                         tg_xy_xxxxxxyz_0, tg_xy_xxxxxxyz_1, tg_xy_xxxxxxz_1, tg_xy_xxxxxxzz_0, \
                                         tg_xy_xxxxxxzz_1, tg_xy_xxxxxyy_1, tg_xy_xxxxxyyy_0, tg_xy_xxxxxyyy_1, \
                                         tg_xy_xxxxxyyz_0, tg_xy_xxxxxyyz_1, tg_xy_xxxxxyz_1, tg_xy_xxxxxyzz_0, \
                                         tg_xy_xxxxxyzz_1, tg_xy_xxxxxzz_1, tg_xy_xxxxxzzz_0, tg_xy_xxxxxzzz_1, \
                                         tg_xy_xxxxyyy_1, tg_xy_xxxxyyyy_0, tg_xy_xxxxyyyy_1, tg_xy_xxxxyyyz_0, \
                                         tg_xy_xxxxyyyz_1, tg_xy_xxxxyyz_1, tg_xy_xxxxyyzz_0, tg_xy_xxxxyyzz_1, \
                                         tg_xy_xxxxyzz_1, tg_xy_xxxxyzzz_0, tg_xy_xxxxyzzz_1, tg_xy_xxxxzzz_1, \
                                         tg_xy_xxxxzzzz_0, tg_xy_xxxxzzzz_1, tg_xy_xxxyyyy_1, tg_xy_xxxyyyyy_0, \
                                         tg_xy_xxxyyyyy_1, tg_xy_xxxyyyyz_0, tg_xy_xxxyyyyz_1, tg_xy_xxxyyyz_1, \
                                         tg_xy_xxxyyyzz_0, tg_xy_xxxyyyzz_1, tg_xy_xxxyyzz_1, tg_xy_xxxyyzzz_0, \
                                         tg_xy_xxxyyzzz_1, tg_xy_xxxyzzz_1, tg_xy_xxxyzzzz_0, tg_xy_xxxyzzzz_1, \
                                         tg_xy_xxxzzzz_1, tg_xy_xxxzzzzz_0, tg_xy_xxxzzzzz_1, tg_xy_xxyyyyy_1, \
                                         tg_xy_xxyyyyyy_0, tg_xy_xxyyyyyy_1, tg_xy_xxyyyyyz_0, tg_xy_xxyyyyyz_1, \
                                         tg_xy_xxyyyyz_1, tg_xy_xxyyyyzz_0, tg_xy_xxyyyyzz_1, tg_xy_xxyyyzz_1, \
                                         tg_xy_xxyyyzzz_0, tg_xy_xxyyyzzz_1, tg_xy_xxyyzzz_1, tg_xy_xxyyzzzz_0, \
                                         tg_xy_xxyyzzzz_1, tg_xy_xxyzzzz_1, tg_xy_xxyzzzzz_0, tg_xy_xxyzzzzz_1, \
                                         tg_xy_xxzzzzz_1, tg_xy_xxzzzzzz_0, tg_xy_xxzzzzzz_1, tg_xy_xyyyyyy_1, \
                                         tg_xy_xyyyyyyy_0, tg_xy_xyyyyyyy_1, tg_xy_xyyyyyyz_0, tg_xy_xyyyyyyz_1, \
                                         tg_xy_xyyyyyz_1, tg_xy_xyyyyyzz_0, tg_xy_xyyyyyzz_1, tg_xy_xyyyyzz_1, \
                                         tg_xy_xyyyyzzz_0, tg_xy_xyyyyzzz_1, tg_xy_xyyyzzz_1, tg_xy_xyyyzzzz_0, \
                                         tg_xy_xyyyzzzz_1, tg_xy_xyyzzzz_1, tg_xy_xyyzzzzz_0, tg_xy_xyyzzzzz_1, \
                                         tg_xy_xyzzzzz_1, tg_xy_xyzzzzzz_0, tg_xy_xyzzzzzz_1, tg_xy_xzzzzzz_1, \
                                         tg_xy_xzzzzzzz_0, tg_xy_xzzzzzzz_1, tg_xy_yyyyyyy_1, tg_xy_yyyyyyyy_0, \
                                         tg_xy_yyyyyyyy_1, tg_xy_yyyyyyyz_0, tg_xy_yyyyyyyz_1, tg_xy_yyyyyyz_1, \
                                         tg_xy_yyyyyyzz_0, tg_xy_yyyyyyzz_1, tg_xy_yyyyyzz_1, tg_xy_yyyyyzzz_0, \
                                         tg_xy_yyyyyzzz_1, tg_xy_yyyyzzz_1, tg_xy_yyyyzzzz_0, tg_xy_yyyyzzzz_1, \
                                         tg_xy_yyyzzzz_1, tg_xy_yyyzzzzz_0, tg_xy_yyyzzzzz_1, tg_xy_yyzzzzz_1, \
                                         tg_xy_yyzzzzzz_0, tg_xy_yyzzzzzz_1, tg_xy_yzzzzzz_1, tg_xy_yzzzzzzz_0, \
                                         tg_xy_yzzzzzzz_1, tg_xy_zzzzzzz_1, tg_xy_zzzzzzzz_0, tg_xy_zzzzzzzz_1, \
                                         tg_y_xxxxxxxx_0, tg_y_xxxxxxxx_1, tg_y_xxxxxxxy_0, tg_y_xxxxxxxy_1, tg_y_xxxxxxxz_0, \
                                         tg_y_xxxxxxxz_1, tg_y_xxxxxxyy_0, tg_y_xxxxxxyy_1, tg_y_xxxxxxyz_0, tg_y_xxxxxxyz_1, \
                                         tg_y_xxxxxxzz_0, tg_y_xxxxxxzz_1, tg_y_xxxxxyyy_0, tg_y_xxxxxyyy_1, tg_y_xxxxxyyz_0, \
                                         tg_y_xxxxxyyz_1, tg_y_xxxxxyzz_0, tg_y_xxxxxyzz_1, tg_y_xxxxxzzz_0, tg_y_xxxxxzzz_1, \
                                         tg_y_xxxxyyyy_0, tg_y_xxxxyyyy_1, tg_y_xxxxyyyz_0, tg_y_xxxxyyyz_1, tg_y_xxxxyyzz_0, \
                                         tg_y_xxxxyyzz_1, tg_y_xxxxyzzz_0, tg_y_xxxxyzzz_1, tg_y_xxxxzzzz_0, tg_y_xxxxzzzz_1, \
                                         tg_y_xxxyyyyy_0, tg_y_xxxyyyyy_1, tg_y_xxxyyyyz_0, tg_y_xxxyyyyz_1, tg_y_xxxyyyzz_0, \
                                         tg_y_xxxyyyzz_1, tg_y_xxxyyzzz_0, tg_y_xxxyyzzz_1, tg_y_xxxyzzzz_0, tg_y_xxxyzzzz_1, \
                                         tg_y_xxxzzzzz_0, tg_y_xxxzzzzz_1, tg_y_xxyyyyyy_0, tg_y_xxyyyyyy_1, tg_y_xxyyyyyz_0, \
                                         tg_y_xxyyyyyz_1, tg_y_xxyyyyzz_0, tg_y_xxyyyyzz_1, tg_y_xxyyyzzz_0, tg_y_xxyyyzzz_1, \
                                         tg_y_xxyyzzzz_0, tg_y_xxyyzzzz_1, tg_y_xxyzzzzz_0, tg_y_xxyzzzzz_1, tg_y_xxzzzzzz_0, \
                                         tg_y_xxzzzzzz_1, tg_y_xyyyyyyy_0, tg_y_xyyyyyyy_1, tg_y_xyyyyyyz_0, tg_y_xyyyyyyz_1, \
                                         tg_y_xyyyyyzz_0, tg_y_xyyyyyzz_1, tg_y_xyyyyzzz_0, tg_y_xyyyyzzz_1, tg_y_xyyyzzzz_0, \
                                         tg_y_xyyyzzzz_1, tg_y_xyyzzzzz_0, tg_y_xyyzzzzz_1, tg_y_xyzzzzzz_0, tg_y_xyzzzzzz_1, \
                                         tg_y_xzzzzzzz_0, tg_y_xzzzzzzz_1, tg_y_yyyyyyyy_0, tg_y_yyyyyyyy_1, tg_y_yyyyyyyz_0, \
                                         tg_y_yyyyyyyz_1, tg_y_yyyyyyzz_0, tg_y_yyyyyyzz_1, tg_y_yyyyyzzz_0, tg_y_yyyyyzzz_1, \
                                         tg_y_yyyyzzzz_0, tg_y_yyyyzzzz_1, tg_y_yyyzzzzz_0, tg_y_yyyzzzzz_1, tg_y_yyzzzzzz_0, \
                                         tg_y_yyzzzzzz_1, tg_y_yzzzzzzz_0, tg_y_yzzzzzzz_1, tg_y_zzzzzzzz_0, tg_y_zzzzzzzz_1, \
                                         wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxx_xxxxxxxx_0[j] = pb_x * tg_xx_xxxxxxxx_0[j] + wp_x[j] * tg_xx_xxxxxxxx_1[j] + fl1_fx * tg_x_xxxxxxxx_0[j] - fl1_fx * fl1_fza * tg_x_xxxxxxxx_1[j] + 4.0 * fl1_fxn * tg_xx_xxxxxxx_1[j];

                    tg_xxx_xxxxxxxy_0[j] = pb_x * tg_xx_xxxxxxxy_0[j] + wp_x[j] * tg_xx_xxxxxxxy_1[j] + fl1_fx * tg_x_xxxxxxxy_0[j] - fl1_fx * fl1_fza * tg_x_xxxxxxxy_1[j] + 3.5 * fl1_fxn * tg_xx_xxxxxxy_1[j];

                    tg_xxx_xxxxxxxz_0[j] = pb_x * tg_xx_xxxxxxxz_0[j] + wp_x[j] * tg_xx_xxxxxxxz_1[j] + fl1_fx * tg_x_xxxxxxxz_0[j] - fl1_fx * fl1_fza * tg_x_xxxxxxxz_1[j] + 3.5 * fl1_fxn * tg_xx_xxxxxxz_1[j];

                    tg_xxx_xxxxxxyy_0[j] = pb_x * tg_xx_xxxxxxyy_0[j] + wp_x[j] * tg_xx_xxxxxxyy_1[j] + fl1_fx * tg_x_xxxxxxyy_0[j] - fl1_fx * fl1_fza * tg_x_xxxxxxyy_1[j] + 3.0 * fl1_fxn * tg_xx_xxxxxyy_1[j];

                    tg_xxx_xxxxxxyz_0[j] = pb_x * tg_xx_xxxxxxyz_0[j] + wp_x[j] * tg_xx_xxxxxxyz_1[j] + fl1_fx * tg_x_xxxxxxyz_0[j] - fl1_fx * fl1_fza * tg_x_xxxxxxyz_1[j] + 3.0 * fl1_fxn * tg_xx_xxxxxyz_1[j];

                    tg_xxx_xxxxxxzz_0[j] = pb_x * tg_xx_xxxxxxzz_0[j] + wp_x[j] * tg_xx_xxxxxxzz_1[j] + fl1_fx * tg_x_xxxxxxzz_0[j] - fl1_fx * fl1_fza * tg_x_xxxxxxzz_1[j] + 3.0 * fl1_fxn * tg_xx_xxxxxzz_1[j];

                    tg_xxx_xxxxxyyy_0[j] = pb_x * tg_xx_xxxxxyyy_0[j] + wp_x[j] * tg_xx_xxxxxyyy_1[j] + fl1_fx * tg_x_xxxxxyyy_0[j] - fl1_fx * fl1_fza * tg_x_xxxxxyyy_1[j] + 2.5 * fl1_fxn * tg_xx_xxxxyyy_1[j];

                    tg_xxx_xxxxxyyz_0[j] = pb_x * tg_xx_xxxxxyyz_0[j] + wp_x[j] * tg_xx_xxxxxyyz_1[j] + fl1_fx * tg_x_xxxxxyyz_0[j] - fl1_fx * fl1_fza * tg_x_xxxxxyyz_1[j] + 2.5 * fl1_fxn * tg_xx_xxxxyyz_1[j];

                    tg_xxx_xxxxxyzz_0[j] = pb_x * tg_xx_xxxxxyzz_0[j] + wp_x[j] * tg_xx_xxxxxyzz_1[j] + fl1_fx * tg_x_xxxxxyzz_0[j] - fl1_fx * fl1_fza * tg_x_xxxxxyzz_1[j] + 2.5 * fl1_fxn * tg_xx_xxxxyzz_1[j];

                    tg_xxx_xxxxxzzz_0[j] = pb_x * tg_xx_xxxxxzzz_0[j] + wp_x[j] * tg_xx_xxxxxzzz_1[j] + fl1_fx * tg_x_xxxxxzzz_0[j] - fl1_fx * fl1_fza * tg_x_xxxxxzzz_1[j] + 2.5 * fl1_fxn * tg_xx_xxxxzzz_1[j];

                    tg_xxx_xxxxyyyy_0[j] = pb_x * tg_xx_xxxxyyyy_0[j] + wp_x[j] * tg_xx_xxxxyyyy_1[j] + fl1_fx * tg_x_xxxxyyyy_0[j] - fl1_fx * fl1_fza * tg_x_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_xx_xxxyyyy_1[j];

                    tg_xxx_xxxxyyyz_0[j] = pb_x * tg_xx_xxxxyyyz_0[j] + wp_x[j] * tg_xx_xxxxyyyz_1[j] + fl1_fx * tg_x_xxxxyyyz_0[j] - fl1_fx * fl1_fza * tg_x_xxxxyyyz_1[j] + 2.0 * fl1_fxn * tg_xx_xxxyyyz_1[j];

                    tg_xxx_xxxxyyzz_0[j] = pb_x * tg_xx_xxxxyyzz_0[j] + wp_x[j] * tg_xx_xxxxyyzz_1[j] + fl1_fx * tg_x_xxxxyyzz_0[j] - fl1_fx * fl1_fza * tg_x_xxxxyyzz_1[j] + 2.0 * fl1_fxn * tg_xx_xxxyyzz_1[j];

                    tg_xxx_xxxxyzzz_0[j] = pb_x * tg_xx_xxxxyzzz_0[j] + wp_x[j] * tg_xx_xxxxyzzz_1[j] + fl1_fx * tg_x_xxxxyzzz_0[j] - fl1_fx * fl1_fza * tg_x_xxxxyzzz_1[j] + 2.0 * fl1_fxn * tg_xx_xxxyzzz_1[j];

                    tg_xxx_xxxxzzzz_0[j] = pb_x * tg_xx_xxxxzzzz_0[j] + wp_x[j] * tg_xx_xxxxzzzz_1[j] + fl1_fx * tg_x_xxxxzzzz_0[j] - fl1_fx * fl1_fza * tg_x_xxxxzzzz_1[j] + 2.0 * fl1_fxn * tg_xx_xxxzzzz_1[j];

                    tg_xxx_xxxyyyyy_0[j] = pb_x * tg_xx_xxxyyyyy_0[j] + wp_x[j] * tg_xx_xxxyyyyy_1[j] + fl1_fx * tg_x_xxxyyyyy_0[j] - fl1_fx * fl1_fza * tg_x_xxxyyyyy_1[j] + 1.5 * fl1_fxn * tg_xx_xxyyyyy_1[j];

                    tg_xxx_xxxyyyyz_0[j] = pb_x * tg_xx_xxxyyyyz_0[j] + wp_x[j] * tg_xx_xxxyyyyz_1[j] + fl1_fx * tg_x_xxxyyyyz_0[j] - fl1_fx * fl1_fza * tg_x_xxxyyyyz_1[j] + 1.5 * fl1_fxn * tg_xx_xxyyyyz_1[j];

                    tg_xxx_xxxyyyzz_0[j] = pb_x * tg_xx_xxxyyyzz_0[j] + wp_x[j] * tg_xx_xxxyyyzz_1[j] + fl1_fx * tg_x_xxxyyyzz_0[j] - fl1_fx * fl1_fza * tg_x_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_xx_xxyyyzz_1[j];

                    tg_xxx_xxxyyzzz_0[j] = pb_x * tg_xx_xxxyyzzz_0[j] + wp_x[j] * tg_xx_xxxyyzzz_1[j] + fl1_fx * tg_x_xxxyyzzz_0[j] - fl1_fx * fl1_fza * tg_x_xxxyyzzz_1[j] + 1.5 * fl1_fxn * tg_xx_xxyyzzz_1[j];

                    tg_xxx_xxxyzzzz_0[j] = pb_x * tg_xx_xxxyzzzz_0[j] + wp_x[j] * tg_xx_xxxyzzzz_1[j] + fl1_fx * tg_x_xxxyzzzz_0[j] - fl1_fx * fl1_fza * tg_x_xxxyzzzz_1[j] + 1.5 * fl1_fxn * tg_xx_xxyzzzz_1[j];

                    tg_xxx_xxxzzzzz_0[j] = pb_x * tg_xx_xxxzzzzz_0[j] + wp_x[j] * tg_xx_xxxzzzzz_1[j] + fl1_fx * tg_x_xxxzzzzz_0[j] - fl1_fx * fl1_fza * tg_x_xxxzzzzz_1[j] + 1.5 * fl1_fxn * tg_xx_xxzzzzz_1[j];

                    tg_xxx_xxyyyyyy_0[j] = pb_x * tg_xx_xxyyyyyy_0[j] + wp_x[j] * tg_xx_xxyyyyyy_1[j] + fl1_fx * tg_x_xxyyyyyy_0[j] - fl1_fx * fl1_fza * tg_x_xxyyyyyy_1[j] + fl1_fxn * tg_xx_xyyyyyy_1[j];

                    tg_xxx_xxyyyyyz_0[j] = pb_x * tg_xx_xxyyyyyz_0[j] + wp_x[j] * tg_xx_xxyyyyyz_1[j] + fl1_fx * tg_x_xxyyyyyz_0[j] - fl1_fx * fl1_fza * tg_x_xxyyyyyz_1[j] + fl1_fxn * tg_xx_xyyyyyz_1[j];

                    tg_xxx_xxyyyyzz_0[j] = pb_x * tg_xx_xxyyyyzz_0[j] + wp_x[j] * tg_xx_xxyyyyzz_1[j] + fl1_fx * tg_x_xxyyyyzz_0[j] - fl1_fx * fl1_fza * tg_x_xxyyyyzz_1[j] + fl1_fxn * tg_xx_xyyyyzz_1[j];

                    tg_xxx_xxyyyzzz_0[j] = pb_x * tg_xx_xxyyyzzz_0[j] + wp_x[j] * tg_xx_xxyyyzzz_1[j] + fl1_fx * tg_x_xxyyyzzz_0[j] - fl1_fx * fl1_fza * tg_x_xxyyyzzz_1[j] + fl1_fxn * tg_xx_xyyyzzz_1[j];

                    tg_xxx_xxyyzzzz_0[j] = pb_x * tg_xx_xxyyzzzz_0[j] + wp_x[j] * tg_xx_xxyyzzzz_1[j] + fl1_fx * tg_x_xxyyzzzz_0[j] - fl1_fx * fl1_fza * tg_x_xxyyzzzz_1[j] + fl1_fxn * tg_xx_xyyzzzz_1[j];

                    tg_xxx_xxyzzzzz_0[j] = pb_x * tg_xx_xxyzzzzz_0[j] + wp_x[j] * tg_xx_xxyzzzzz_1[j] + fl1_fx * tg_x_xxyzzzzz_0[j] - fl1_fx * fl1_fza * tg_x_xxyzzzzz_1[j] + fl1_fxn * tg_xx_xyzzzzz_1[j];

                    tg_xxx_xxzzzzzz_0[j] = pb_x * tg_xx_xxzzzzzz_0[j] + wp_x[j] * tg_xx_xxzzzzzz_1[j] + fl1_fx * tg_x_xxzzzzzz_0[j] - fl1_fx * fl1_fza * tg_x_xxzzzzzz_1[j] + fl1_fxn * tg_xx_xzzzzzz_1[j];

                    tg_xxx_xyyyyyyy_0[j] = pb_x * tg_xx_xyyyyyyy_0[j] + wp_x[j] * tg_xx_xyyyyyyy_1[j] + fl1_fx * tg_x_xyyyyyyy_0[j] - fl1_fx * fl1_fza * tg_x_xyyyyyyy_1[j] + 0.5 * fl1_fxn * tg_xx_yyyyyyy_1[j];

                    tg_xxx_xyyyyyyz_0[j] = pb_x * tg_xx_xyyyyyyz_0[j] + wp_x[j] * tg_xx_xyyyyyyz_1[j] + fl1_fx * tg_x_xyyyyyyz_0[j] - fl1_fx * fl1_fza * tg_x_xyyyyyyz_1[j] + 0.5 * fl1_fxn * tg_xx_yyyyyyz_1[j];

                    tg_xxx_xyyyyyzz_0[j] = pb_x * tg_xx_xyyyyyzz_0[j] + wp_x[j] * tg_xx_xyyyyyzz_1[j] + fl1_fx * tg_x_xyyyyyzz_0[j] - fl1_fx * fl1_fza * tg_x_xyyyyyzz_1[j] + 0.5 * fl1_fxn * tg_xx_yyyyyzz_1[j];

                    tg_xxx_xyyyyzzz_0[j] = pb_x * tg_xx_xyyyyzzz_0[j] + wp_x[j] * tg_xx_xyyyyzzz_1[j] + fl1_fx * tg_x_xyyyyzzz_0[j] - fl1_fx * fl1_fza * tg_x_xyyyyzzz_1[j] + 0.5 * fl1_fxn * tg_xx_yyyyzzz_1[j];

                    tg_xxx_xyyyzzzz_0[j] = pb_x * tg_xx_xyyyzzzz_0[j] + wp_x[j] * tg_xx_xyyyzzzz_1[j] + fl1_fx * tg_x_xyyyzzzz_0[j] - fl1_fx * fl1_fza * tg_x_xyyyzzzz_1[j] + 0.5 * fl1_fxn * tg_xx_yyyzzzz_1[j];

                    tg_xxx_xyyzzzzz_0[j] = pb_x * tg_xx_xyyzzzzz_0[j] + wp_x[j] * tg_xx_xyyzzzzz_1[j] + fl1_fx * tg_x_xyyzzzzz_0[j] - fl1_fx * fl1_fza * tg_x_xyyzzzzz_1[j] + 0.5 * fl1_fxn * tg_xx_yyzzzzz_1[j];

                    tg_xxx_xyzzzzzz_0[j] = pb_x * tg_xx_xyzzzzzz_0[j] + wp_x[j] * tg_xx_xyzzzzzz_1[j] + fl1_fx * tg_x_xyzzzzzz_0[j] - fl1_fx * fl1_fza * tg_x_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xx_yzzzzzz_1[j];

                    tg_xxx_xzzzzzzz_0[j] = pb_x * tg_xx_xzzzzzzz_0[j] + wp_x[j] * tg_xx_xzzzzzzz_1[j] + fl1_fx * tg_x_xzzzzzzz_0[j] - fl1_fx * fl1_fza * tg_x_xzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xx_zzzzzzz_1[j];

                    tg_xxx_yyyyyyyy_0[j] = pb_x * tg_xx_yyyyyyyy_0[j] + wp_x[j] * tg_xx_yyyyyyyy_1[j] + fl1_fx * tg_x_yyyyyyyy_0[j] - fl1_fx * fl1_fza * tg_x_yyyyyyyy_1[j];

                    tg_xxx_yyyyyyyz_0[j] = pb_x * tg_xx_yyyyyyyz_0[j] + wp_x[j] * tg_xx_yyyyyyyz_1[j] + fl1_fx * tg_x_yyyyyyyz_0[j] - fl1_fx * fl1_fza * tg_x_yyyyyyyz_1[j];

                    tg_xxx_yyyyyyzz_0[j] = pb_x * tg_xx_yyyyyyzz_0[j] + wp_x[j] * tg_xx_yyyyyyzz_1[j] + fl1_fx * tg_x_yyyyyyzz_0[j] - fl1_fx * fl1_fza * tg_x_yyyyyyzz_1[j];

                    tg_xxx_yyyyyzzz_0[j] = pb_x * tg_xx_yyyyyzzz_0[j] + wp_x[j] * tg_xx_yyyyyzzz_1[j] + fl1_fx * tg_x_yyyyyzzz_0[j] - fl1_fx * fl1_fza * tg_x_yyyyyzzz_1[j];

                    tg_xxx_yyyyzzzz_0[j] = pb_x * tg_xx_yyyyzzzz_0[j] + wp_x[j] * tg_xx_yyyyzzzz_1[j] + fl1_fx * tg_x_yyyyzzzz_0[j] - fl1_fx * fl1_fza * tg_x_yyyyzzzz_1[j];

                    tg_xxx_yyyzzzzz_0[j] = pb_x * tg_xx_yyyzzzzz_0[j] + wp_x[j] * tg_xx_yyyzzzzz_1[j] + fl1_fx * tg_x_yyyzzzzz_0[j] - fl1_fx * fl1_fza * tg_x_yyyzzzzz_1[j];

                    tg_xxx_yyzzzzzz_0[j] = pb_x * tg_xx_yyzzzzzz_0[j] + wp_x[j] * tg_xx_yyzzzzzz_1[j] + fl1_fx * tg_x_yyzzzzzz_0[j] - fl1_fx * fl1_fza * tg_x_yyzzzzzz_1[j];

                    tg_xxx_yzzzzzzz_0[j] = pb_x * tg_xx_yzzzzzzz_0[j] + wp_x[j] * tg_xx_yzzzzzzz_1[j] + fl1_fx * tg_x_yzzzzzzz_0[j] - fl1_fx * fl1_fza * tg_x_yzzzzzzz_1[j];

                    tg_xxx_zzzzzzzz_0[j] = pb_x * tg_xx_zzzzzzzz_0[j] + wp_x[j] * tg_xx_zzzzzzzz_1[j] + fl1_fx * tg_x_zzzzzzzz_0[j] - fl1_fx * fl1_fza * tg_x_zzzzzzzz_1[j];

                    tg_xxy_xxxxxxxx_0[j] = pb_x * tg_xy_xxxxxxxx_0[j] + wp_x[j] * tg_xy_xxxxxxxx_1[j] + 0.5 * fl1_fx * tg_y_xxxxxxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxxxxxx_1[j] + 4.0 * fl1_fxn * tg_xy_xxxxxxx_1[j];

                    tg_xxy_xxxxxxxy_0[j] = pb_x * tg_xy_xxxxxxxy_0[j] + wp_x[j] * tg_xy_xxxxxxxy_1[j] + 0.5 * fl1_fx * tg_y_xxxxxxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxxxxxy_1[j] + 3.5 * fl1_fxn * tg_xy_xxxxxxy_1[j];

                    tg_xxy_xxxxxxxz_0[j] = pb_x * tg_xy_xxxxxxxz_0[j] + wp_x[j] * tg_xy_xxxxxxxz_1[j] + 0.5 * fl1_fx * tg_y_xxxxxxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxxxxxz_1[j] + 3.5 * fl1_fxn * tg_xy_xxxxxxz_1[j];

                    tg_xxy_xxxxxxyy_0[j] = pb_x * tg_xy_xxxxxxyy_0[j] + wp_x[j] * tg_xy_xxxxxxyy_1[j] + 0.5 * fl1_fx * tg_y_xxxxxxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxxxxyy_1[j] + 3.0 * fl1_fxn * tg_xy_xxxxxyy_1[j];

                    tg_xxy_xxxxxxyz_0[j] = pb_x * tg_xy_xxxxxxyz_0[j] + wp_x[j] * tg_xy_xxxxxxyz_1[j] + 0.5 * fl1_fx * tg_y_xxxxxxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxxxxyz_1[j] + 3.0 * fl1_fxn * tg_xy_xxxxxyz_1[j];

                    tg_xxy_xxxxxxzz_0[j] = pb_x * tg_xy_xxxxxxzz_0[j] + wp_x[j] * tg_xy_xxxxxxzz_1[j] + 0.5 * fl1_fx * tg_y_xxxxxxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxxxxzz_1[j] + 3.0 * fl1_fxn * tg_xy_xxxxxzz_1[j];

                    tg_xxy_xxxxxyyy_0[j] = pb_x * tg_xy_xxxxxyyy_0[j] + wp_x[j] * tg_xy_xxxxxyyy_1[j] + 0.5 * fl1_fx * tg_y_xxxxxyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxxxyyy_1[j] + 2.5 * fl1_fxn * tg_xy_xxxxyyy_1[j];

                    tg_xxy_xxxxxyyz_0[j] = pb_x * tg_xy_xxxxxyyz_0[j] + wp_x[j] * tg_xy_xxxxxyyz_1[j] + 0.5 * fl1_fx * tg_y_xxxxxyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxxxyyz_1[j] + 2.5 * fl1_fxn * tg_xy_xxxxyyz_1[j];

                    tg_xxy_xxxxxyzz_0[j] = pb_x * tg_xy_xxxxxyzz_0[j] + wp_x[j] * tg_xy_xxxxxyzz_1[j] + 0.5 * fl1_fx * tg_y_xxxxxyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxxxyzz_1[j] + 2.5 * fl1_fxn * tg_xy_xxxxyzz_1[j];

                    tg_xxy_xxxxxzzz_0[j] = pb_x * tg_xy_xxxxxzzz_0[j] + wp_x[j] * tg_xy_xxxxxzzz_1[j] + 0.5 * fl1_fx * tg_y_xxxxxzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxxxzzz_1[j] + 2.5 * fl1_fxn * tg_xy_xxxxzzz_1[j];

                    tg_xxy_xxxxyyyy_0[j] = pb_x * tg_xy_xxxxyyyy_0[j] + wp_x[j] * tg_xy_xxxxyyyy_1[j] + 0.5 * fl1_fx * tg_y_xxxxyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_xy_xxxyyyy_1[j];

                    tg_xxy_xxxxyyyz_0[j] = pb_x * tg_xy_xxxxyyyz_0[j] + wp_x[j] * tg_xy_xxxxyyyz_1[j] + 0.5 * fl1_fx * tg_y_xxxxyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxxyyyz_1[j] + 2.0 * fl1_fxn * tg_xy_xxxyyyz_1[j];

                    tg_xxy_xxxxyyzz_0[j] = pb_x * tg_xy_xxxxyyzz_0[j] + wp_x[j] * tg_xy_xxxxyyzz_1[j] + 0.5 * fl1_fx * tg_y_xxxxyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxxyyzz_1[j] + 2.0 * fl1_fxn * tg_xy_xxxyyzz_1[j];

                    tg_xxy_xxxxyzzz_0[j] = pb_x * tg_xy_xxxxyzzz_0[j] + wp_x[j] * tg_xy_xxxxyzzz_1[j] + 0.5 * fl1_fx * tg_y_xxxxyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxxyzzz_1[j] + 2.0 * fl1_fxn * tg_xy_xxxyzzz_1[j];

                    tg_xxy_xxxxzzzz_0[j] = pb_x * tg_xy_xxxxzzzz_0[j] + wp_x[j] * tg_xy_xxxxzzzz_1[j] + 0.5 * fl1_fx * tg_y_xxxxzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxxzzzz_1[j] + 2.0 * fl1_fxn * tg_xy_xxxzzzz_1[j];

                    tg_xxy_xxxyyyyy_0[j] = pb_x * tg_xy_xxxyyyyy_0[j] + wp_x[j] * tg_xy_xxxyyyyy_1[j] + 0.5 * fl1_fx * tg_y_xxxyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxyyyyy_1[j] + 1.5 * fl1_fxn * tg_xy_xxyyyyy_1[j];

                    tg_xxy_xxxyyyyz_0[j] = pb_x * tg_xy_xxxyyyyz_0[j] + wp_x[j] * tg_xy_xxxyyyyz_1[j] + 0.5 * fl1_fx * tg_y_xxxyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxyyyyz_1[j] + 1.5 * fl1_fxn * tg_xy_xxyyyyz_1[j];

                    tg_xxy_xxxyyyzz_0[j] = pb_x * tg_xy_xxxyyyzz_0[j] + wp_x[j] * tg_xy_xxxyyyzz_1[j] + 0.5 * fl1_fx * tg_y_xxxyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_xy_xxyyyzz_1[j];

                    tg_xxy_xxxyyzzz_0[j] = pb_x * tg_xy_xxxyyzzz_0[j] + wp_x[j] * tg_xy_xxxyyzzz_1[j] + 0.5 * fl1_fx * tg_y_xxxyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxyyzzz_1[j] + 1.5 * fl1_fxn * tg_xy_xxyyzzz_1[j];

                    tg_xxy_xxxyzzzz_0[j] = pb_x * tg_xy_xxxyzzzz_0[j] + wp_x[j] * tg_xy_xxxyzzzz_1[j] + 0.5 * fl1_fx * tg_y_xxxyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxyzzzz_1[j] + 1.5 * fl1_fxn * tg_xy_xxyzzzz_1[j];

                    tg_xxy_xxxzzzzz_0[j] = pb_x * tg_xy_xxxzzzzz_0[j] + wp_x[j] * tg_xy_xxxzzzzz_1[j] + 0.5 * fl1_fx * tg_y_xxxzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxzzzzz_1[j] + 1.5 * fl1_fxn * tg_xy_xxzzzzz_1[j];

                    tg_xxy_xxyyyyyy_0[j] = pb_x * tg_xy_xxyyyyyy_0[j] + wp_x[j] * tg_xy_xxyyyyyy_1[j] + 0.5 * fl1_fx * tg_y_xxyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxyyyyyy_1[j] + fl1_fxn * tg_xy_xyyyyyy_1[j];

                    tg_xxy_xxyyyyyz_0[j] = pb_x * tg_xy_xxyyyyyz_0[j] + wp_x[j] * tg_xy_xxyyyyyz_1[j] + 0.5 * fl1_fx * tg_y_xxyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxyyyyyz_1[j] + fl1_fxn * tg_xy_xyyyyyz_1[j];

                    tg_xxy_xxyyyyzz_0[j] = pb_x * tg_xy_xxyyyyzz_0[j] + wp_x[j] * tg_xy_xxyyyyzz_1[j] + 0.5 * fl1_fx * tg_y_xxyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxyyyyzz_1[j] + fl1_fxn * tg_xy_xyyyyzz_1[j];

                    tg_xxy_xxyyyzzz_0[j] = pb_x * tg_xy_xxyyyzzz_0[j] + wp_x[j] * tg_xy_xxyyyzzz_1[j] + 0.5 * fl1_fx * tg_y_xxyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxyyyzzz_1[j] + fl1_fxn * tg_xy_xyyyzzz_1[j];

                    tg_xxy_xxyyzzzz_0[j] = pb_x * tg_xy_xxyyzzzz_0[j] + wp_x[j] * tg_xy_xxyyzzzz_1[j] + 0.5 * fl1_fx * tg_y_xxyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxyyzzzz_1[j] + fl1_fxn * tg_xy_xyyzzzz_1[j];

                    tg_xxy_xxyzzzzz_0[j] = pb_x * tg_xy_xxyzzzzz_0[j] + wp_x[j] * tg_xy_xxyzzzzz_1[j] + 0.5 * fl1_fx * tg_y_xxyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxyzzzzz_1[j] + fl1_fxn * tg_xy_xyzzzzz_1[j];

                    tg_xxy_xxzzzzzz_0[j] = pb_x * tg_xy_xxzzzzzz_0[j] + wp_x[j] * tg_xy_xxzzzzzz_1[j] + 0.5 * fl1_fx * tg_y_xxzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxzzzzzz_1[j] + fl1_fxn * tg_xy_xzzzzzz_1[j];

                    tg_xxy_xyyyyyyy_0[j] = pb_x * tg_xy_xyyyyyyy_0[j] + wp_x[j] * tg_xy_xyyyyyyy_1[j] + 0.5 * fl1_fx * tg_y_xyyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xyyyyyyy_1[j] + 0.5 * fl1_fxn * tg_xy_yyyyyyy_1[j];

                    tg_xxy_xyyyyyyz_0[j] = pb_x * tg_xy_xyyyyyyz_0[j] + wp_x[j] * tg_xy_xyyyyyyz_1[j] + 0.5 * fl1_fx * tg_y_xyyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xyyyyyyz_1[j] + 0.5 * fl1_fxn * tg_xy_yyyyyyz_1[j];

                    tg_xxy_xyyyyyzz_0[j] = pb_x * tg_xy_xyyyyyzz_0[j] + wp_x[j] * tg_xy_xyyyyyzz_1[j] + 0.5 * fl1_fx * tg_y_xyyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xyyyyyzz_1[j] + 0.5 * fl1_fxn * tg_xy_yyyyyzz_1[j];

                    tg_xxy_xyyyyzzz_0[j] = pb_x * tg_xy_xyyyyzzz_0[j] + wp_x[j] * tg_xy_xyyyyzzz_1[j] + 0.5 * fl1_fx * tg_y_xyyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xyyyyzzz_1[j] + 0.5 * fl1_fxn * tg_xy_yyyyzzz_1[j];

                    tg_xxy_xyyyzzzz_0[j] = pb_x * tg_xy_xyyyzzzz_0[j] + wp_x[j] * tg_xy_xyyyzzzz_1[j] + 0.5 * fl1_fx * tg_y_xyyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xyyyzzzz_1[j] + 0.5 * fl1_fxn * tg_xy_yyyzzzz_1[j];

                    tg_xxy_xyyzzzzz_0[j] = pb_x * tg_xy_xyyzzzzz_0[j] + wp_x[j] * tg_xy_xyyzzzzz_1[j] + 0.5 * fl1_fx * tg_y_xyyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xyyzzzzz_1[j] + 0.5 * fl1_fxn * tg_xy_yyzzzzz_1[j];

                    tg_xxy_xyzzzzzz_0[j] = pb_x * tg_xy_xyzzzzzz_0[j] + wp_x[j] * tg_xy_xyzzzzzz_1[j] + 0.5 * fl1_fx * tg_y_xyzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xy_yzzzzzz_1[j];

                    tg_xxy_xzzzzzzz_0[j] = pb_x * tg_xy_xzzzzzzz_0[j] + wp_x[j] * tg_xy_xzzzzzzz_1[j] + 0.5 * fl1_fx * tg_y_xzzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xy_zzzzzzz_1[j];

                    tg_xxy_yyyyyyyy_0[j] = pb_x * tg_xy_yyyyyyyy_0[j] + wp_x[j] * tg_xy_yyyyyyyy_1[j] + 0.5 * fl1_fx * tg_y_yyyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_yyyyyyyy_1[j];

                    tg_xxy_yyyyyyyz_0[j] = pb_x * tg_xy_yyyyyyyz_0[j] + wp_x[j] * tg_xy_yyyyyyyz_1[j] + 0.5 * fl1_fx * tg_y_yyyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_yyyyyyyz_1[j];

                    tg_xxy_yyyyyyzz_0[j] = pb_x * tg_xy_yyyyyyzz_0[j] + wp_x[j] * tg_xy_yyyyyyzz_1[j] + 0.5 * fl1_fx * tg_y_yyyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_yyyyyyzz_1[j];

                    tg_xxy_yyyyyzzz_0[j] = pb_x * tg_xy_yyyyyzzz_0[j] + wp_x[j] * tg_xy_yyyyyzzz_1[j] + 0.5 * fl1_fx * tg_y_yyyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_yyyyyzzz_1[j];

                    tg_xxy_yyyyzzzz_0[j] = pb_x * tg_xy_yyyyzzzz_0[j] + wp_x[j] * tg_xy_yyyyzzzz_1[j] + 0.5 * fl1_fx * tg_y_yyyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_yyyyzzzz_1[j];

                    tg_xxy_yyyzzzzz_0[j] = pb_x * tg_xy_yyyzzzzz_0[j] + wp_x[j] * tg_xy_yyyzzzzz_1[j] + 0.5 * fl1_fx * tg_y_yyyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_yyyzzzzz_1[j];

                    tg_xxy_yyzzzzzz_0[j] = pb_x * tg_xy_yyzzzzzz_0[j] + wp_x[j] * tg_xy_yyzzzzzz_1[j] + 0.5 * fl1_fx * tg_y_yyzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_yyzzzzzz_1[j];

                    tg_xxy_yzzzzzzz_0[j] = pb_x * tg_xy_yzzzzzzz_0[j] + wp_x[j] * tg_xy_yzzzzzzz_1[j] + 0.5 * fl1_fx * tg_y_yzzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_yzzzzzzz_1[j];

                    tg_xxy_zzzzzzzz_0[j] = pb_x * tg_xy_zzzzzzzz_0[j] + wp_x[j] * tg_xy_zzzzzzzz_1[j] + 0.5 * fl1_fx * tg_y_zzzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_zzzzzzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSFSL_90_180(      CMemBlock2D<double>& primBuffer,
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
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_3_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_3_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_2_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_2_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_1_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_1_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_2_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {7, -1, -1, -1}, 
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

                auto tg_xz_xxxxxxxx_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 90); 

                auto tg_xz_xxxxxxxy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 91); 

                auto tg_xz_xxxxxxxz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 92); 

                auto tg_xz_xxxxxxyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 93); 

                auto tg_xz_xxxxxxyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 94); 

                auto tg_xz_xxxxxxzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 95); 

                auto tg_xz_xxxxxyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 96); 

                auto tg_xz_xxxxxyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 97); 

                auto tg_xz_xxxxxyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 98); 

                auto tg_xz_xxxxxzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 99); 

                auto tg_xz_xxxxyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 100); 

                auto tg_xz_xxxxyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 101); 

                auto tg_xz_xxxxyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 102); 

                auto tg_xz_xxxxyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 103); 

                auto tg_xz_xxxxzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 104); 

                auto tg_xz_xxxyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 105); 

                auto tg_xz_xxxyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 106); 

                auto tg_xz_xxxyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 107); 

                auto tg_xz_xxxyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 108); 

                auto tg_xz_xxxyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 109); 

                auto tg_xz_xxxzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 110); 

                auto tg_xz_xxyyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 111); 

                auto tg_xz_xxyyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 112); 

                auto tg_xz_xxyyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 113); 

                auto tg_xz_xxyyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 114); 

                auto tg_xz_xxyyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 115); 

                auto tg_xz_xxyzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 116); 

                auto tg_xz_xxzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 117); 

                auto tg_xz_xyyyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 118); 

                auto tg_xz_xyyyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 119); 

                auto tg_xz_xyyyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 120); 

                auto tg_xz_xyyyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 121); 

                auto tg_xz_xyyyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 122); 

                auto tg_xz_xyyzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 123); 

                auto tg_xz_xyzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 124); 

                auto tg_xz_xzzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 125); 

                auto tg_xz_yyyyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 126); 

                auto tg_xz_yyyyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 127); 

                auto tg_xz_yyyyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 128); 

                auto tg_xz_yyyyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 129); 

                auto tg_xz_yyyyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 130); 

                auto tg_xz_yyyzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 131); 

                auto tg_xz_yyzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 132); 

                auto tg_xz_yzzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 133); 

                auto tg_xz_zzzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 134); 

                auto tg_yy_xxxxxxxx_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 135); 

                auto tg_yy_xxxxxxxy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 136); 

                auto tg_yy_xxxxxxxz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 137); 

                auto tg_yy_xxxxxxyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 138); 

                auto tg_yy_xxxxxxyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 139); 

                auto tg_yy_xxxxxxzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 140); 

                auto tg_yy_xxxxxyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 141); 

                auto tg_yy_xxxxxyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 142); 

                auto tg_yy_xxxxxyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 143); 

                auto tg_yy_xxxxxzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 144); 

                auto tg_yy_xxxxyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 145); 

                auto tg_yy_xxxxyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 146); 

                auto tg_yy_xxxxyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 147); 

                auto tg_yy_xxxxyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 148); 

                auto tg_yy_xxxxzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 149); 

                auto tg_yy_xxxyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 150); 

                auto tg_yy_xxxyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 151); 

                auto tg_yy_xxxyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 152); 

                auto tg_yy_xxxyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 153); 

                auto tg_yy_xxxyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 154); 

                auto tg_yy_xxxzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 155); 

                auto tg_yy_xxyyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 156); 

                auto tg_yy_xxyyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 157); 

                auto tg_yy_xxyyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 158); 

                auto tg_yy_xxyyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 159); 

                auto tg_yy_xxyyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 160); 

                auto tg_yy_xxyzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 161); 

                auto tg_yy_xxzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 162); 

                auto tg_yy_xyyyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 163); 

                auto tg_yy_xyyyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 164); 

                auto tg_yy_xyyyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 165); 

                auto tg_yy_xyyyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 166); 

                auto tg_yy_xyyyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 167); 

                auto tg_yy_xyyzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 168); 

                auto tg_yy_xyzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 169); 

                auto tg_yy_xzzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 170); 

                auto tg_yy_yyyyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 171); 

                auto tg_yy_yyyyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 172); 

                auto tg_yy_yyyyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 173); 

                auto tg_yy_yyyyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 174); 

                auto tg_yy_yyyyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 175); 

                auto tg_yy_yyyzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 176); 

                auto tg_yy_yyzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 177); 

                auto tg_yy_yzzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 178); 

                auto tg_yy_zzzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 179); 

                auto tg_xz_xxxxxxxx_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 90); 

                auto tg_xz_xxxxxxxy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 91); 

                auto tg_xz_xxxxxxxz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 92); 

                auto tg_xz_xxxxxxyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 93); 

                auto tg_xz_xxxxxxyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 94); 

                auto tg_xz_xxxxxxzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 95); 

                auto tg_xz_xxxxxyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 96); 

                auto tg_xz_xxxxxyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 97); 

                auto tg_xz_xxxxxyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 98); 

                auto tg_xz_xxxxxzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 99); 

                auto tg_xz_xxxxyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 100); 

                auto tg_xz_xxxxyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 101); 

                auto tg_xz_xxxxyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 102); 

                auto tg_xz_xxxxyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 103); 

                auto tg_xz_xxxxzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 104); 

                auto tg_xz_xxxyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 105); 

                auto tg_xz_xxxyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 106); 

                auto tg_xz_xxxyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 107); 

                auto tg_xz_xxxyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 108); 

                auto tg_xz_xxxyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 109); 

                auto tg_xz_xxxzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 110); 

                auto tg_xz_xxyyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 111); 

                auto tg_xz_xxyyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 112); 

                auto tg_xz_xxyyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 113); 

                auto tg_xz_xxyyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 114); 

                auto tg_xz_xxyyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 115); 

                auto tg_xz_xxyzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 116); 

                auto tg_xz_xxzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 117); 

                auto tg_xz_xyyyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 118); 

                auto tg_xz_xyyyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 119); 

                auto tg_xz_xyyyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 120); 

                auto tg_xz_xyyyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 121); 

                auto tg_xz_xyyyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 122); 

                auto tg_xz_xyyzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 123); 

                auto tg_xz_xyzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 124); 

                auto tg_xz_xzzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 125); 

                auto tg_xz_yyyyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 126); 

                auto tg_xz_yyyyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 127); 

                auto tg_xz_yyyyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 128); 

                auto tg_xz_yyyyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 129); 

                auto tg_xz_yyyyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 130); 

                auto tg_xz_yyyzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 131); 

                auto tg_xz_yyzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 132); 

                auto tg_xz_yzzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 133); 

                auto tg_xz_zzzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 134); 

                auto tg_yy_xxxxxxxx_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 135); 

                auto tg_yy_xxxxxxxy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 136); 

                auto tg_yy_xxxxxxxz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 137); 

                auto tg_yy_xxxxxxyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 138); 

                auto tg_yy_xxxxxxyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 139); 

                auto tg_yy_xxxxxxzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 140); 

                auto tg_yy_xxxxxyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 141); 

                auto tg_yy_xxxxxyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 142); 

                auto tg_yy_xxxxxyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 143); 

                auto tg_yy_xxxxxzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 144); 

                auto tg_yy_xxxxyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 145); 

                auto tg_yy_xxxxyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 146); 

                auto tg_yy_xxxxyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 147); 

                auto tg_yy_xxxxyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 148); 

                auto tg_yy_xxxxzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 149); 

                auto tg_yy_xxxyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 150); 

                auto tg_yy_xxxyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 151); 

                auto tg_yy_xxxyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 152); 

                auto tg_yy_xxxyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 153); 

                auto tg_yy_xxxyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 154); 

                auto tg_yy_xxxzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 155); 

                auto tg_yy_xxyyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 156); 

                auto tg_yy_xxyyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 157); 

                auto tg_yy_xxyyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 158); 

                auto tg_yy_xxyyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 159); 

                auto tg_yy_xxyyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 160); 

                auto tg_yy_xxyzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 161); 

                auto tg_yy_xxzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 162); 

                auto tg_yy_xyyyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 163); 

                auto tg_yy_xyyyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 164); 

                auto tg_yy_xyyyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 165); 

                auto tg_yy_xyyyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 166); 

                auto tg_yy_xyyyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 167); 

                auto tg_yy_xyyzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 168); 

                auto tg_yy_xyzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 169); 

                auto tg_yy_xzzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 170); 

                auto tg_yy_yyyyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 171); 

                auto tg_yy_yyyyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 172); 

                auto tg_yy_yyyyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 173); 

                auto tg_yy_yyyyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 174); 

                auto tg_yy_yyyyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 175); 

                auto tg_yy_yyyzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 176); 

                auto tg_yy_yyzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 177); 

                auto tg_yy_yzzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 178); 

                auto tg_yy_zzzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 179); 

                auto tg_z_xxxxxxxx_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 90); 

                auto tg_z_xxxxxxxy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 91); 

                auto tg_z_xxxxxxxz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 92); 

                auto tg_z_xxxxxxyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 93); 

                auto tg_z_xxxxxxyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 94); 

                auto tg_z_xxxxxxzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 95); 

                auto tg_z_xxxxxyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 96); 

                auto tg_z_xxxxxyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 97); 

                auto tg_z_xxxxxyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 98); 

                auto tg_z_xxxxxzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 99); 

                auto tg_z_xxxxyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 100); 

                auto tg_z_xxxxyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 101); 

                auto tg_z_xxxxyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 102); 

                auto tg_z_xxxxyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 103); 

                auto tg_z_xxxxzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 104); 

                auto tg_z_xxxyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 105); 

                auto tg_z_xxxyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 106); 

                auto tg_z_xxxyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 107); 

                auto tg_z_xxxyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 108); 

                auto tg_z_xxxyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 109); 

                auto tg_z_xxxzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 110); 

                auto tg_z_xxyyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 111); 

                auto tg_z_xxyyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 112); 

                auto tg_z_xxyyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 113); 

                auto tg_z_xxyyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 114); 

                auto tg_z_xxyyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 115); 

                auto tg_z_xxyzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 116); 

                auto tg_z_xxzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 117); 

                auto tg_z_xyyyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 118); 

                auto tg_z_xyyyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 119); 

                auto tg_z_xyyyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 120); 

                auto tg_z_xyyyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 121); 

                auto tg_z_xyyyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 122); 

                auto tg_z_xyyzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 123); 

                auto tg_z_xyzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 124); 

                auto tg_z_xzzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 125); 

                auto tg_z_yyyyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 126); 

                auto tg_z_yyyyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 127); 

                auto tg_z_yyyyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 128); 

                auto tg_z_yyyyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 129); 

                auto tg_z_yyyyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 130); 

                auto tg_z_yyyzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 131); 

                auto tg_z_yyzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 132); 

                auto tg_z_yzzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 133); 

                auto tg_z_zzzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 134); 

                auto tg_z_xxxxxxxx_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 90); 

                auto tg_z_xxxxxxxy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 91); 

                auto tg_z_xxxxxxxz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 92); 

                auto tg_z_xxxxxxyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 93); 

                auto tg_z_xxxxxxyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 94); 

                auto tg_z_xxxxxxzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 95); 

                auto tg_z_xxxxxyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 96); 

                auto tg_z_xxxxxyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 97); 

                auto tg_z_xxxxxyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 98); 

                auto tg_z_xxxxxzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 99); 

                auto tg_z_xxxxyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 100); 

                auto tg_z_xxxxyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 101); 

                auto tg_z_xxxxyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 102); 

                auto tg_z_xxxxyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 103); 

                auto tg_z_xxxxzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 104); 

                auto tg_z_xxxyyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 105); 

                auto tg_z_xxxyyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 106); 

                auto tg_z_xxxyyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 107); 

                auto tg_z_xxxyyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 108); 

                auto tg_z_xxxyzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 109); 

                auto tg_z_xxxzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 110); 

                auto tg_z_xxyyyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 111); 

                auto tg_z_xxyyyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 112); 

                auto tg_z_xxyyyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 113); 

                auto tg_z_xxyyyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 114); 

                auto tg_z_xxyyzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 115); 

                auto tg_z_xxyzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 116); 

                auto tg_z_xxzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 117); 

                auto tg_z_xyyyyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 118); 

                auto tg_z_xyyyyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 119); 

                auto tg_z_xyyyyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 120); 

                auto tg_z_xyyyyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 121); 

                auto tg_z_xyyyzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 122); 

                auto tg_z_xyyzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 123); 

                auto tg_z_xyzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 124); 

                auto tg_z_xzzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 125); 

                auto tg_z_yyyyyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 126); 

                auto tg_z_yyyyyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 127); 

                auto tg_z_yyyyyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 128); 

                auto tg_z_yyyyyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 129); 

                auto tg_z_yyyyzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 130); 

                auto tg_z_yyyzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 131); 

                auto tg_z_yyzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 132); 

                auto tg_z_yzzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 133); 

                auto tg_z_zzzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 134); 

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

                // set up pointers to integrals

                auto tg_xxz_xxxxxxxx_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 90); 

                auto tg_xxz_xxxxxxxy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 91); 

                auto tg_xxz_xxxxxxxz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 92); 

                auto tg_xxz_xxxxxxyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 93); 

                auto tg_xxz_xxxxxxyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 94); 

                auto tg_xxz_xxxxxxzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 95); 

                auto tg_xxz_xxxxxyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 96); 

                auto tg_xxz_xxxxxyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 97); 

                auto tg_xxz_xxxxxyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 98); 

                auto tg_xxz_xxxxxzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 99); 

                auto tg_xxz_xxxxyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 100); 

                auto tg_xxz_xxxxyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 101); 

                auto tg_xxz_xxxxyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 102); 

                auto tg_xxz_xxxxyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 103); 

                auto tg_xxz_xxxxzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 104); 

                auto tg_xxz_xxxyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 105); 

                auto tg_xxz_xxxyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 106); 

                auto tg_xxz_xxxyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 107); 

                auto tg_xxz_xxxyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 108); 

                auto tg_xxz_xxxyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 109); 

                auto tg_xxz_xxxzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 110); 

                auto tg_xxz_xxyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 111); 

                auto tg_xxz_xxyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 112); 

                auto tg_xxz_xxyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 113); 

                auto tg_xxz_xxyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 114); 

                auto tg_xxz_xxyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 115); 

                auto tg_xxz_xxyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 116); 

                auto tg_xxz_xxzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 117); 

                auto tg_xxz_xyyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 118); 

                auto tg_xxz_xyyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 119); 

                auto tg_xxz_xyyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 120); 

                auto tg_xxz_xyyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 121); 

                auto tg_xxz_xyyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 122); 

                auto tg_xxz_xyyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 123); 

                auto tg_xxz_xyzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 124); 

                auto tg_xxz_xzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 125); 

                auto tg_xxz_yyyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 126); 

                auto tg_xxz_yyyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 127); 

                auto tg_xxz_yyyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 128); 

                auto tg_xxz_yyyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 129); 

                auto tg_xxz_yyyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 130); 

                auto tg_xxz_yyyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 131); 

                auto tg_xxz_yyzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 132); 

                auto tg_xxz_yzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 133); 

                auto tg_xxz_zzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 134); 

                auto tg_xyy_xxxxxxxx_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 135); 

                auto tg_xyy_xxxxxxxy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 136); 

                auto tg_xyy_xxxxxxxz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 137); 

                auto tg_xyy_xxxxxxyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 138); 

                auto tg_xyy_xxxxxxyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 139); 

                auto tg_xyy_xxxxxxzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 140); 

                auto tg_xyy_xxxxxyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 141); 

                auto tg_xyy_xxxxxyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 142); 

                auto tg_xyy_xxxxxyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 143); 

                auto tg_xyy_xxxxxzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 144); 

                auto tg_xyy_xxxxyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 145); 

                auto tg_xyy_xxxxyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 146); 

                auto tg_xyy_xxxxyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 147); 

                auto tg_xyy_xxxxyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 148); 

                auto tg_xyy_xxxxzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 149); 

                auto tg_xyy_xxxyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 150); 

                auto tg_xyy_xxxyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 151); 

                auto tg_xyy_xxxyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 152); 

                auto tg_xyy_xxxyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 153); 

                auto tg_xyy_xxxyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 154); 

                auto tg_xyy_xxxzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 155); 

                auto tg_xyy_xxyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 156); 

                auto tg_xyy_xxyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 157); 

                auto tg_xyy_xxyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 158); 

                auto tg_xyy_xxyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 159); 

                auto tg_xyy_xxyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 160); 

                auto tg_xyy_xxyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 161); 

                auto tg_xyy_xxzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 162); 

                auto tg_xyy_xyyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 163); 

                auto tg_xyy_xyyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 164); 

                auto tg_xyy_xyyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 165); 

                auto tg_xyy_xyyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 166); 

                auto tg_xyy_xyyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 167); 

                auto tg_xyy_xyyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 168); 

                auto tg_xyy_xyzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 169); 

                auto tg_xyy_xzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 170); 

                auto tg_xyy_yyyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 171); 

                auto tg_xyy_yyyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 172); 

                auto tg_xyy_yyyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 173); 

                auto tg_xyy_yyyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 174); 

                auto tg_xyy_yyyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 175); 

                auto tg_xyy_yyyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 176); 

                auto tg_xyy_yyzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 177); 

                auto tg_xyy_yzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 178); 

                auto tg_xyy_zzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 179); 

                // Batch of Integrals (90,180)

                #pragma omp simd aligned(fxn, fza, tg_xxz_xxxxxxxx_0, tg_xxz_xxxxxxxy_0, tg_xxz_xxxxxxxz_0, \
                                         tg_xxz_xxxxxxyy_0, tg_xxz_xxxxxxyz_0, tg_xxz_xxxxxxzz_0, tg_xxz_xxxxxyyy_0, \
                                         tg_xxz_xxxxxyyz_0, tg_xxz_xxxxxyzz_0, tg_xxz_xxxxxzzz_0, tg_xxz_xxxxyyyy_0, \
                                         tg_xxz_xxxxyyyz_0, tg_xxz_xxxxyyzz_0, tg_xxz_xxxxyzzz_0, tg_xxz_xxxxzzzz_0, \
                                         tg_xxz_xxxyyyyy_0, tg_xxz_xxxyyyyz_0, tg_xxz_xxxyyyzz_0, tg_xxz_xxxyyzzz_0, \
                                         tg_xxz_xxxyzzzz_0, tg_xxz_xxxzzzzz_0, tg_xxz_xxyyyyyy_0, tg_xxz_xxyyyyyz_0, \
                                         tg_xxz_xxyyyyzz_0, tg_xxz_xxyyyzzz_0, tg_xxz_xxyyzzzz_0, tg_xxz_xxyzzzzz_0, \
                                         tg_xxz_xxzzzzzz_0, tg_xxz_xyyyyyyy_0, tg_xxz_xyyyyyyz_0, tg_xxz_xyyyyyzz_0, \
                                         tg_xxz_xyyyyzzz_0, tg_xxz_xyyyzzzz_0, tg_xxz_xyyzzzzz_0, tg_xxz_xyzzzzzz_0, \
                                         tg_xxz_xzzzzzzz_0, tg_xxz_yyyyyyyy_0, tg_xxz_yyyyyyyz_0, tg_xxz_yyyyyyzz_0, \
                                         tg_xxz_yyyyyzzz_0, tg_xxz_yyyyzzzz_0, tg_xxz_yyyzzzzz_0, tg_xxz_yyzzzzzz_0, \
                                         tg_xxz_yzzzzzzz_0, tg_xxz_zzzzzzzz_0, tg_xyy_xxxxxxxx_0, tg_xyy_xxxxxxxy_0, \
                                         tg_xyy_xxxxxxxz_0, tg_xyy_xxxxxxyy_0, tg_xyy_xxxxxxyz_0, tg_xyy_xxxxxxzz_0, \
                                         tg_xyy_xxxxxyyy_0, tg_xyy_xxxxxyyz_0, tg_xyy_xxxxxyzz_0, tg_xyy_xxxxxzzz_0, \
                                         tg_xyy_xxxxyyyy_0, tg_xyy_xxxxyyyz_0, tg_xyy_xxxxyyzz_0, tg_xyy_xxxxyzzz_0, \
                                         tg_xyy_xxxxzzzz_0, tg_xyy_xxxyyyyy_0, tg_xyy_xxxyyyyz_0, tg_xyy_xxxyyyzz_0, \
                                         tg_xyy_xxxyyzzz_0, tg_xyy_xxxyzzzz_0, tg_xyy_xxxzzzzz_0, tg_xyy_xxyyyyyy_0, \
                                         tg_xyy_xxyyyyyz_0, tg_xyy_xxyyyyzz_0, tg_xyy_xxyyyzzz_0, tg_xyy_xxyyzzzz_0, \
                                         tg_xyy_xxyzzzzz_0, tg_xyy_xxzzzzzz_0, tg_xyy_xyyyyyyy_0, tg_xyy_xyyyyyyz_0, \
                                         tg_xyy_xyyyyyzz_0, tg_xyy_xyyyyzzz_0, tg_xyy_xyyyzzzz_0, tg_xyy_xyyzzzzz_0, \
                                         tg_xyy_xyzzzzzz_0, tg_xyy_xzzzzzzz_0, tg_xyy_yyyyyyyy_0, tg_xyy_yyyyyyyz_0, \
                                         tg_xyy_yyyyyyzz_0, tg_xyy_yyyyyzzz_0, tg_xyy_yyyyzzzz_0, tg_xyy_yyyzzzzz_0, \
                                         tg_xyy_yyzzzzzz_0, tg_xyy_yzzzzzzz_0, tg_xyy_zzzzzzzz_0, tg_xz_xxxxxxx_1, \
                                         tg_xz_xxxxxxxx_0, tg_xz_xxxxxxxx_1, tg_xz_xxxxxxxy_0, tg_xz_xxxxxxxy_1, \
                                         tg_xz_xxxxxxxz_0, tg_xz_xxxxxxxz_1, tg_xz_xxxxxxy_1, tg_xz_xxxxxxyy_0, \
                                         tg_xz_xxxxxxyy_1, tg_xz_xxxxxxyz_0, tg_xz_xxxxxxyz_1, tg_xz_xxxxxxz_1, \
                                         tg_xz_xxxxxxzz_0, tg_xz_xxxxxxzz_1, tg_xz_xxxxxyy_1, tg_xz_xxxxxyyy_0, \
                                         tg_xz_xxxxxyyy_1, tg_xz_xxxxxyyz_0, tg_xz_xxxxxyyz_1, tg_xz_xxxxxyz_1, \
                                         tg_xz_xxxxxyzz_0, tg_xz_xxxxxyzz_1, tg_xz_xxxxxzz_1, tg_xz_xxxxxzzz_0, \
                                         tg_xz_xxxxxzzz_1, tg_xz_xxxxyyy_1, tg_xz_xxxxyyyy_0, tg_xz_xxxxyyyy_1, \
                                         tg_xz_xxxxyyyz_0, tg_xz_xxxxyyyz_1, tg_xz_xxxxyyz_1, tg_xz_xxxxyyzz_0, \
                                         tg_xz_xxxxyyzz_1, tg_xz_xxxxyzz_1, tg_xz_xxxxyzzz_0, tg_xz_xxxxyzzz_1, \
                                         tg_xz_xxxxzzz_1, tg_xz_xxxxzzzz_0, tg_xz_xxxxzzzz_1, tg_xz_xxxyyyy_1, \
                                         tg_xz_xxxyyyyy_0, tg_xz_xxxyyyyy_1, tg_xz_xxxyyyyz_0, tg_xz_xxxyyyyz_1, \
                                         tg_xz_xxxyyyz_1, tg_xz_xxxyyyzz_0, tg_xz_xxxyyyzz_1, tg_xz_xxxyyzz_1, \
                                         tg_xz_xxxyyzzz_0, tg_xz_xxxyyzzz_1, tg_xz_xxxyzzz_1, tg_xz_xxxyzzzz_0, \
                                         tg_xz_xxxyzzzz_1, tg_xz_xxxzzzz_1, tg_xz_xxxzzzzz_0, tg_xz_xxxzzzzz_1, \
                                         tg_xz_xxyyyyy_1, tg_xz_xxyyyyyy_0, tg_xz_xxyyyyyy_1, tg_xz_xxyyyyyz_0, \
                                         tg_xz_xxyyyyyz_1, tg_xz_xxyyyyz_1, tg_xz_xxyyyyzz_0, tg_xz_xxyyyyzz_1, \
                                         tg_xz_xxyyyzz_1, tg_xz_xxyyyzzz_0, tg_xz_xxyyyzzz_1, tg_xz_xxyyzzz_1, \
                                         tg_xz_xxyyzzzz_0, tg_xz_xxyyzzzz_1, tg_xz_xxyzzzz_1, tg_xz_xxyzzzzz_0, \
                                         tg_xz_xxyzzzzz_1, tg_xz_xxzzzzz_1, tg_xz_xxzzzzzz_0, tg_xz_xxzzzzzz_1, \
                                         tg_xz_xyyyyyy_1, tg_xz_xyyyyyyy_0, tg_xz_xyyyyyyy_1, tg_xz_xyyyyyyz_0, \
                                         tg_xz_xyyyyyyz_1, tg_xz_xyyyyyz_1, tg_xz_xyyyyyzz_0, tg_xz_xyyyyyzz_1, \
                                         tg_xz_xyyyyzz_1, tg_xz_xyyyyzzz_0, tg_xz_xyyyyzzz_1, tg_xz_xyyyzzz_1, \
                                         tg_xz_xyyyzzzz_0, tg_xz_xyyyzzzz_1, tg_xz_xyyzzzz_1, tg_xz_xyyzzzzz_0, \
                                         tg_xz_xyyzzzzz_1, tg_xz_xyzzzzz_1, tg_xz_xyzzzzzz_0, tg_xz_xyzzzzzz_1, \
                                         tg_xz_xzzzzzz_1, tg_xz_xzzzzzzz_0, tg_xz_xzzzzzzz_1, tg_xz_yyyyyyy_1, \
                                         tg_xz_yyyyyyyy_0, tg_xz_yyyyyyyy_1, tg_xz_yyyyyyyz_0, tg_xz_yyyyyyyz_1, \
                                         tg_xz_yyyyyyz_1, tg_xz_yyyyyyzz_0, tg_xz_yyyyyyzz_1, tg_xz_yyyyyzz_1, \
                                         tg_xz_yyyyyzzz_0, tg_xz_yyyyyzzz_1, tg_xz_yyyyzzz_1, tg_xz_yyyyzzzz_0, \
                                         tg_xz_yyyyzzzz_1, tg_xz_yyyzzzz_1, tg_xz_yyyzzzzz_0, tg_xz_yyyzzzzz_1, \
                                         tg_xz_yyzzzzz_1, tg_xz_yyzzzzzz_0, tg_xz_yyzzzzzz_1, tg_xz_yzzzzzz_1, \
                                         tg_xz_yzzzzzzz_0, tg_xz_yzzzzzzz_1, tg_xz_zzzzzzz_1, tg_xz_zzzzzzzz_0, \
                                         tg_xz_zzzzzzzz_1, tg_yy_xxxxxxx_1, tg_yy_xxxxxxxx_0, tg_yy_xxxxxxxx_1, \
                                         tg_yy_xxxxxxxy_0, tg_yy_xxxxxxxy_1, tg_yy_xxxxxxxz_0, tg_yy_xxxxxxxz_1, \
                                         tg_yy_xxxxxxy_1, tg_yy_xxxxxxyy_0, tg_yy_xxxxxxyy_1, tg_yy_xxxxxxyz_0, \
                                         tg_yy_xxxxxxyz_1, tg_yy_xxxxxxz_1, tg_yy_xxxxxxzz_0, tg_yy_xxxxxxzz_1, \
                                         tg_yy_xxxxxyy_1, tg_yy_xxxxxyyy_0, tg_yy_xxxxxyyy_1, tg_yy_xxxxxyyz_0, \
                                         tg_yy_xxxxxyyz_1, tg_yy_xxxxxyz_1, tg_yy_xxxxxyzz_0, tg_yy_xxxxxyzz_1, \
                                         tg_yy_xxxxxzz_1, tg_yy_xxxxxzzz_0, tg_yy_xxxxxzzz_1, tg_yy_xxxxyyy_1, \
                                         tg_yy_xxxxyyyy_0, tg_yy_xxxxyyyy_1, tg_yy_xxxxyyyz_0, tg_yy_xxxxyyyz_1, \
                                         tg_yy_xxxxyyz_1, tg_yy_xxxxyyzz_0, tg_yy_xxxxyyzz_1, tg_yy_xxxxyzz_1, \
                                         tg_yy_xxxxyzzz_0, tg_yy_xxxxyzzz_1, tg_yy_xxxxzzz_1, tg_yy_xxxxzzzz_0, \
                                         tg_yy_xxxxzzzz_1, tg_yy_xxxyyyy_1, tg_yy_xxxyyyyy_0, tg_yy_xxxyyyyy_1, \
                                         tg_yy_xxxyyyyz_0, tg_yy_xxxyyyyz_1, tg_yy_xxxyyyz_1, tg_yy_xxxyyyzz_0, \
                                         tg_yy_xxxyyyzz_1, tg_yy_xxxyyzz_1, tg_yy_xxxyyzzz_0, tg_yy_xxxyyzzz_1, \
                                         tg_yy_xxxyzzz_1, tg_yy_xxxyzzzz_0, tg_yy_xxxyzzzz_1, tg_yy_xxxzzzz_1, \
                                         tg_yy_xxxzzzzz_0, tg_yy_xxxzzzzz_1, tg_yy_xxyyyyy_1, tg_yy_xxyyyyyy_0, \
                                         tg_yy_xxyyyyyy_1, tg_yy_xxyyyyyz_0, tg_yy_xxyyyyyz_1, tg_yy_xxyyyyz_1, \
                                         tg_yy_xxyyyyzz_0, tg_yy_xxyyyyzz_1, tg_yy_xxyyyzz_1, tg_yy_xxyyyzzz_0, \
                                         tg_yy_xxyyyzzz_1, tg_yy_xxyyzzz_1, tg_yy_xxyyzzzz_0, tg_yy_xxyyzzzz_1, \
                                         tg_yy_xxyzzzz_1, tg_yy_xxyzzzzz_0, tg_yy_xxyzzzzz_1, tg_yy_xxzzzzz_1, \
                                         tg_yy_xxzzzzzz_0, tg_yy_xxzzzzzz_1, tg_yy_xyyyyyy_1, tg_yy_xyyyyyyy_0, \
                                         tg_yy_xyyyyyyy_1, tg_yy_xyyyyyyz_0, tg_yy_xyyyyyyz_1, tg_yy_xyyyyyz_1, \
                                         tg_yy_xyyyyyzz_0, tg_yy_xyyyyyzz_1, tg_yy_xyyyyzz_1, tg_yy_xyyyyzzz_0, \
                                         tg_yy_xyyyyzzz_1, tg_yy_xyyyzzz_1, tg_yy_xyyyzzzz_0, tg_yy_xyyyzzzz_1, \
                                         tg_yy_xyyzzzz_1, tg_yy_xyyzzzzz_0, tg_yy_xyyzzzzz_1, tg_yy_xyzzzzz_1, \
                                         tg_yy_xyzzzzzz_0, tg_yy_xyzzzzzz_1, tg_yy_xzzzzzz_1, tg_yy_xzzzzzzz_0, \
                                         tg_yy_xzzzzzzz_1, tg_yy_yyyyyyy_1, tg_yy_yyyyyyyy_0, tg_yy_yyyyyyyy_1, \
                                         tg_yy_yyyyyyyz_0, tg_yy_yyyyyyyz_1, tg_yy_yyyyyyz_1, tg_yy_yyyyyyzz_0, \
                                         tg_yy_yyyyyyzz_1, tg_yy_yyyyyzz_1, tg_yy_yyyyyzzz_0, tg_yy_yyyyyzzz_1, \
                                         tg_yy_yyyyzzz_1, tg_yy_yyyyzzzz_0, tg_yy_yyyyzzzz_1, tg_yy_yyyzzzz_1, \
                                         tg_yy_yyyzzzzz_0, tg_yy_yyyzzzzz_1, tg_yy_yyzzzzz_1, tg_yy_yyzzzzzz_0, \
                                         tg_yy_yyzzzzzz_1, tg_yy_yzzzzzz_1, tg_yy_yzzzzzzz_0, tg_yy_yzzzzzzz_1, \
                                         tg_yy_zzzzzzz_1, tg_yy_zzzzzzzz_0, tg_yy_zzzzzzzz_1, tg_z_xxxxxxxx_0, \
                                         tg_z_xxxxxxxx_1, tg_z_xxxxxxxy_0, tg_z_xxxxxxxy_1, tg_z_xxxxxxxz_0, tg_z_xxxxxxxz_1, \
                                         tg_z_xxxxxxyy_0, tg_z_xxxxxxyy_1, tg_z_xxxxxxyz_0, tg_z_xxxxxxyz_1, tg_z_xxxxxxzz_0, \
                                         tg_z_xxxxxxzz_1, tg_z_xxxxxyyy_0, tg_z_xxxxxyyy_1, tg_z_xxxxxyyz_0, tg_z_xxxxxyyz_1, \
                                         tg_z_xxxxxyzz_0, tg_z_xxxxxyzz_1, tg_z_xxxxxzzz_0, tg_z_xxxxxzzz_1, tg_z_xxxxyyyy_0, \
                                         tg_z_xxxxyyyy_1, tg_z_xxxxyyyz_0, tg_z_xxxxyyyz_1, tg_z_xxxxyyzz_0, tg_z_xxxxyyzz_1, \
                                         tg_z_xxxxyzzz_0, tg_z_xxxxyzzz_1, tg_z_xxxxzzzz_0, tg_z_xxxxzzzz_1, tg_z_xxxyyyyy_0, \
                                         tg_z_xxxyyyyy_1, tg_z_xxxyyyyz_0, tg_z_xxxyyyyz_1, tg_z_xxxyyyzz_0, tg_z_xxxyyyzz_1, \
                                         tg_z_xxxyyzzz_0, tg_z_xxxyyzzz_1, tg_z_xxxyzzzz_0, tg_z_xxxyzzzz_1, tg_z_xxxzzzzz_0, \
                                         tg_z_xxxzzzzz_1, tg_z_xxyyyyyy_0, tg_z_xxyyyyyy_1, tg_z_xxyyyyyz_0, tg_z_xxyyyyyz_1, \
                                         tg_z_xxyyyyzz_0, tg_z_xxyyyyzz_1, tg_z_xxyyyzzz_0, tg_z_xxyyyzzz_1, tg_z_xxyyzzzz_0, \
                                         tg_z_xxyyzzzz_1, tg_z_xxyzzzzz_0, tg_z_xxyzzzzz_1, tg_z_xxzzzzzz_0, tg_z_xxzzzzzz_1, \
                                         tg_z_xyyyyyyy_0, tg_z_xyyyyyyy_1, tg_z_xyyyyyyz_0, tg_z_xyyyyyyz_1, tg_z_xyyyyyzz_0, \
                                         tg_z_xyyyyyzz_1, tg_z_xyyyyzzz_0, tg_z_xyyyyzzz_1, tg_z_xyyyzzzz_0, tg_z_xyyyzzzz_1, \
                                         tg_z_xyyzzzzz_0, tg_z_xyyzzzzz_1, tg_z_xyzzzzzz_0, tg_z_xyzzzzzz_1, tg_z_xzzzzzzz_0, \
                                         tg_z_xzzzzzzz_1, tg_z_yyyyyyyy_0, tg_z_yyyyyyyy_1, tg_z_yyyyyyyz_0, tg_z_yyyyyyyz_1, \
                                         tg_z_yyyyyyzz_0, tg_z_yyyyyyzz_1, tg_z_yyyyyzzz_0, tg_z_yyyyyzzz_1, tg_z_yyyyzzzz_0, \
                                         tg_z_yyyyzzzz_1, tg_z_yyyzzzzz_0, tg_z_yyyzzzzz_1, tg_z_yyzzzzzz_0, tg_z_yyzzzzzz_1, \
                                         tg_z_yzzzzzzz_0, tg_z_yzzzzzzz_1, tg_z_zzzzzzzz_0, tg_z_zzzzzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxz_xxxxxxxx_0[j] = pb_x * tg_xz_xxxxxxxx_0[j] + wp_x[j] * tg_xz_xxxxxxxx_1[j] + 0.5 * fl1_fx * tg_z_xxxxxxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxxxx_1[j] + 4.0 * fl1_fxn * tg_xz_xxxxxxx_1[j];

                    tg_xxz_xxxxxxxy_0[j] = pb_x * tg_xz_xxxxxxxy_0[j] + wp_x[j] * tg_xz_xxxxxxxy_1[j] + 0.5 * fl1_fx * tg_z_xxxxxxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxxxy_1[j] + 3.5 * fl1_fxn * tg_xz_xxxxxxy_1[j];

                    tg_xxz_xxxxxxxz_0[j] = pb_x * tg_xz_xxxxxxxz_0[j] + wp_x[j] * tg_xz_xxxxxxxz_1[j] + 0.5 * fl1_fx * tg_z_xxxxxxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxxxz_1[j] + 3.5 * fl1_fxn * tg_xz_xxxxxxz_1[j];

                    tg_xxz_xxxxxxyy_0[j] = pb_x * tg_xz_xxxxxxyy_0[j] + wp_x[j] * tg_xz_xxxxxxyy_1[j] + 0.5 * fl1_fx * tg_z_xxxxxxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxxyy_1[j] + 3.0 * fl1_fxn * tg_xz_xxxxxyy_1[j];

                    tg_xxz_xxxxxxyz_0[j] = pb_x * tg_xz_xxxxxxyz_0[j] + wp_x[j] * tg_xz_xxxxxxyz_1[j] + 0.5 * fl1_fx * tg_z_xxxxxxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxxyz_1[j] + 3.0 * fl1_fxn * tg_xz_xxxxxyz_1[j];

                    tg_xxz_xxxxxxzz_0[j] = pb_x * tg_xz_xxxxxxzz_0[j] + wp_x[j] * tg_xz_xxxxxxzz_1[j] + 0.5 * fl1_fx * tg_z_xxxxxxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxxzz_1[j] + 3.0 * fl1_fxn * tg_xz_xxxxxzz_1[j];

                    tg_xxz_xxxxxyyy_0[j] = pb_x * tg_xz_xxxxxyyy_0[j] + wp_x[j] * tg_xz_xxxxxyyy_1[j] + 0.5 * fl1_fx * tg_z_xxxxxyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxyyy_1[j] + 2.5 * fl1_fxn * tg_xz_xxxxyyy_1[j];

                    tg_xxz_xxxxxyyz_0[j] = pb_x * tg_xz_xxxxxyyz_0[j] + wp_x[j] * tg_xz_xxxxxyyz_1[j] + 0.5 * fl1_fx * tg_z_xxxxxyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxyyz_1[j] + 2.5 * fl1_fxn * tg_xz_xxxxyyz_1[j];

                    tg_xxz_xxxxxyzz_0[j] = pb_x * tg_xz_xxxxxyzz_0[j] + wp_x[j] * tg_xz_xxxxxyzz_1[j] + 0.5 * fl1_fx * tg_z_xxxxxyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxyzz_1[j] + 2.5 * fl1_fxn * tg_xz_xxxxyzz_1[j];

                    tg_xxz_xxxxxzzz_0[j] = pb_x * tg_xz_xxxxxzzz_0[j] + wp_x[j] * tg_xz_xxxxxzzz_1[j] + 0.5 * fl1_fx * tg_z_xxxxxzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxzzz_1[j] + 2.5 * fl1_fxn * tg_xz_xxxxzzz_1[j];

                    tg_xxz_xxxxyyyy_0[j] = pb_x * tg_xz_xxxxyyyy_0[j] + wp_x[j] * tg_xz_xxxxyyyy_1[j] + 0.5 * fl1_fx * tg_z_xxxxyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_xz_xxxyyyy_1[j];

                    tg_xxz_xxxxyyyz_0[j] = pb_x * tg_xz_xxxxyyyz_0[j] + wp_x[j] * tg_xz_xxxxyyyz_1[j] + 0.5 * fl1_fx * tg_z_xxxxyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxyyyz_1[j] + 2.0 * fl1_fxn * tg_xz_xxxyyyz_1[j];

                    tg_xxz_xxxxyyzz_0[j] = pb_x * tg_xz_xxxxyyzz_0[j] + wp_x[j] * tg_xz_xxxxyyzz_1[j] + 0.5 * fl1_fx * tg_z_xxxxyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxyyzz_1[j] + 2.0 * fl1_fxn * tg_xz_xxxyyzz_1[j];

                    tg_xxz_xxxxyzzz_0[j] = pb_x * tg_xz_xxxxyzzz_0[j] + wp_x[j] * tg_xz_xxxxyzzz_1[j] + 0.5 * fl1_fx * tg_z_xxxxyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxyzzz_1[j] + 2.0 * fl1_fxn * tg_xz_xxxyzzz_1[j];

                    tg_xxz_xxxxzzzz_0[j] = pb_x * tg_xz_xxxxzzzz_0[j] + wp_x[j] * tg_xz_xxxxzzzz_1[j] + 0.5 * fl1_fx * tg_z_xxxxzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxzzzz_1[j] + 2.0 * fl1_fxn * tg_xz_xxxzzzz_1[j];

                    tg_xxz_xxxyyyyy_0[j] = pb_x * tg_xz_xxxyyyyy_0[j] + wp_x[j] * tg_xz_xxxyyyyy_1[j] + 0.5 * fl1_fx * tg_z_xxxyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxyyyyy_1[j] + 1.5 * fl1_fxn * tg_xz_xxyyyyy_1[j];

                    tg_xxz_xxxyyyyz_0[j] = pb_x * tg_xz_xxxyyyyz_0[j] + wp_x[j] * tg_xz_xxxyyyyz_1[j] + 0.5 * fl1_fx * tg_z_xxxyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxyyyyz_1[j] + 1.5 * fl1_fxn * tg_xz_xxyyyyz_1[j];

                    tg_xxz_xxxyyyzz_0[j] = pb_x * tg_xz_xxxyyyzz_0[j] + wp_x[j] * tg_xz_xxxyyyzz_1[j] + 0.5 * fl1_fx * tg_z_xxxyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_xz_xxyyyzz_1[j];

                    tg_xxz_xxxyyzzz_0[j] = pb_x * tg_xz_xxxyyzzz_0[j] + wp_x[j] * tg_xz_xxxyyzzz_1[j] + 0.5 * fl1_fx * tg_z_xxxyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxyyzzz_1[j] + 1.5 * fl1_fxn * tg_xz_xxyyzzz_1[j];

                    tg_xxz_xxxyzzzz_0[j] = pb_x * tg_xz_xxxyzzzz_0[j] + wp_x[j] * tg_xz_xxxyzzzz_1[j] + 0.5 * fl1_fx * tg_z_xxxyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxyzzzz_1[j] + 1.5 * fl1_fxn * tg_xz_xxyzzzz_1[j];

                    tg_xxz_xxxzzzzz_0[j] = pb_x * tg_xz_xxxzzzzz_0[j] + wp_x[j] * tg_xz_xxxzzzzz_1[j] + 0.5 * fl1_fx * tg_z_xxxzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxzzzzz_1[j] + 1.5 * fl1_fxn * tg_xz_xxzzzzz_1[j];

                    tg_xxz_xxyyyyyy_0[j] = pb_x * tg_xz_xxyyyyyy_0[j] + wp_x[j] * tg_xz_xxyyyyyy_1[j] + 0.5 * fl1_fx * tg_z_xxyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxyyyyyy_1[j] + fl1_fxn * tg_xz_xyyyyyy_1[j];

                    tg_xxz_xxyyyyyz_0[j] = pb_x * tg_xz_xxyyyyyz_0[j] + wp_x[j] * tg_xz_xxyyyyyz_1[j] + 0.5 * fl1_fx * tg_z_xxyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxyyyyyz_1[j] + fl1_fxn * tg_xz_xyyyyyz_1[j];

                    tg_xxz_xxyyyyzz_0[j] = pb_x * tg_xz_xxyyyyzz_0[j] + wp_x[j] * tg_xz_xxyyyyzz_1[j] + 0.5 * fl1_fx * tg_z_xxyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxyyyyzz_1[j] + fl1_fxn * tg_xz_xyyyyzz_1[j];

                    tg_xxz_xxyyyzzz_0[j] = pb_x * tg_xz_xxyyyzzz_0[j] + wp_x[j] * tg_xz_xxyyyzzz_1[j] + 0.5 * fl1_fx * tg_z_xxyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxyyyzzz_1[j] + fl1_fxn * tg_xz_xyyyzzz_1[j];

                    tg_xxz_xxyyzzzz_0[j] = pb_x * tg_xz_xxyyzzzz_0[j] + wp_x[j] * tg_xz_xxyyzzzz_1[j] + 0.5 * fl1_fx * tg_z_xxyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxyyzzzz_1[j] + fl1_fxn * tg_xz_xyyzzzz_1[j];

                    tg_xxz_xxyzzzzz_0[j] = pb_x * tg_xz_xxyzzzzz_0[j] + wp_x[j] * tg_xz_xxyzzzzz_1[j] + 0.5 * fl1_fx * tg_z_xxyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxyzzzzz_1[j] + fl1_fxn * tg_xz_xyzzzzz_1[j];

                    tg_xxz_xxzzzzzz_0[j] = pb_x * tg_xz_xxzzzzzz_0[j] + wp_x[j] * tg_xz_xxzzzzzz_1[j] + 0.5 * fl1_fx * tg_z_xxzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxzzzzzz_1[j] + fl1_fxn * tg_xz_xzzzzzz_1[j];

                    tg_xxz_xyyyyyyy_0[j] = pb_x * tg_xz_xyyyyyyy_0[j] + wp_x[j] * tg_xz_xyyyyyyy_1[j] + 0.5 * fl1_fx * tg_z_xyyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyyyyyyy_1[j] + 0.5 * fl1_fxn * tg_xz_yyyyyyy_1[j];

                    tg_xxz_xyyyyyyz_0[j] = pb_x * tg_xz_xyyyyyyz_0[j] + wp_x[j] * tg_xz_xyyyyyyz_1[j] + 0.5 * fl1_fx * tg_z_xyyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyyyyyyz_1[j] + 0.5 * fl1_fxn * tg_xz_yyyyyyz_1[j];

                    tg_xxz_xyyyyyzz_0[j] = pb_x * tg_xz_xyyyyyzz_0[j] + wp_x[j] * tg_xz_xyyyyyzz_1[j] + 0.5 * fl1_fx * tg_z_xyyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyyyyyzz_1[j] + 0.5 * fl1_fxn * tg_xz_yyyyyzz_1[j];

                    tg_xxz_xyyyyzzz_0[j] = pb_x * tg_xz_xyyyyzzz_0[j] + wp_x[j] * tg_xz_xyyyyzzz_1[j] + 0.5 * fl1_fx * tg_z_xyyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyyyyzzz_1[j] + 0.5 * fl1_fxn * tg_xz_yyyyzzz_1[j];

                    tg_xxz_xyyyzzzz_0[j] = pb_x * tg_xz_xyyyzzzz_0[j] + wp_x[j] * tg_xz_xyyyzzzz_1[j] + 0.5 * fl1_fx * tg_z_xyyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyyyzzzz_1[j] + 0.5 * fl1_fxn * tg_xz_yyyzzzz_1[j];

                    tg_xxz_xyyzzzzz_0[j] = pb_x * tg_xz_xyyzzzzz_0[j] + wp_x[j] * tg_xz_xyyzzzzz_1[j] + 0.5 * fl1_fx * tg_z_xyyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyyzzzzz_1[j] + 0.5 * fl1_fxn * tg_xz_yyzzzzz_1[j];

                    tg_xxz_xyzzzzzz_0[j] = pb_x * tg_xz_xyzzzzzz_0[j] + wp_x[j] * tg_xz_xyzzzzzz_1[j] + 0.5 * fl1_fx * tg_z_xyzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xz_yzzzzzz_1[j];

                    tg_xxz_xzzzzzzz_0[j] = pb_x * tg_xz_xzzzzzzz_0[j] + wp_x[j] * tg_xz_xzzzzzzz_1[j] + 0.5 * fl1_fx * tg_z_xzzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xz_zzzzzzz_1[j];

                    tg_xxz_yyyyyyyy_0[j] = pb_x * tg_xz_yyyyyyyy_0[j] + wp_x[j] * tg_xz_yyyyyyyy_1[j] + 0.5 * fl1_fx * tg_z_yyyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyyyyyyy_1[j];

                    tg_xxz_yyyyyyyz_0[j] = pb_x * tg_xz_yyyyyyyz_0[j] + wp_x[j] * tg_xz_yyyyyyyz_1[j] + 0.5 * fl1_fx * tg_z_yyyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyyyyyyz_1[j];

                    tg_xxz_yyyyyyzz_0[j] = pb_x * tg_xz_yyyyyyzz_0[j] + wp_x[j] * tg_xz_yyyyyyzz_1[j] + 0.5 * fl1_fx * tg_z_yyyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyyyyyzz_1[j];

                    tg_xxz_yyyyyzzz_0[j] = pb_x * tg_xz_yyyyyzzz_0[j] + wp_x[j] * tg_xz_yyyyyzzz_1[j] + 0.5 * fl1_fx * tg_z_yyyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyyyyzzz_1[j];

                    tg_xxz_yyyyzzzz_0[j] = pb_x * tg_xz_yyyyzzzz_0[j] + wp_x[j] * tg_xz_yyyyzzzz_1[j] + 0.5 * fl1_fx * tg_z_yyyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyyyzzzz_1[j];

                    tg_xxz_yyyzzzzz_0[j] = pb_x * tg_xz_yyyzzzzz_0[j] + wp_x[j] * tg_xz_yyyzzzzz_1[j] + 0.5 * fl1_fx * tg_z_yyyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyyzzzzz_1[j];

                    tg_xxz_yyzzzzzz_0[j] = pb_x * tg_xz_yyzzzzzz_0[j] + wp_x[j] * tg_xz_yyzzzzzz_1[j] + 0.5 * fl1_fx * tg_z_yyzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyzzzzzz_1[j];

                    tg_xxz_yzzzzzzz_0[j] = pb_x * tg_xz_yzzzzzzz_0[j] + wp_x[j] * tg_xz_yzzzzzzz_1[j] + 0.5 * fl1_fx * tg_z_yzzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yzzzzzzz_1[j];

                    tg_xxz_zzzzzzzz_0[j] = pb_x * tg_xz_zzzzzzzz_0[j] + wp_x[j] * tg_xz_zzzzzzzz_1[j] + 0.5 * fl1_fx * tg_z_zzzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_zzzzzzzz_1[j];

                    tg_xyy_xxxxxxxx_0[j] = pb_x * tg_yy_xxxxxxxx_0[j] + wp_x[j] * tg_yy_xxxxxxxx_1[j] + 4.0 * fl1_fxn * tg_yy_xxxxxxx_1[j];

                    tg_xyy_xxxxxxxy_0[j] = pb_x * tg_yy_xxxxxxxy_0[j] + wp_x[j] * tg_yy_xxxxxxxy_1[j] + 3.5 * fl1_fxn * tg_yy_xxxxxxy_1[j];

                    tg_xyy_xxxxxxxz_0[j] = pb_x * tg_yy_xxxxxxxz_0[j] + wp_x[j] * tg_yy_xxxxxxxz_1[j] + 3.5 * fl1_fxn * tg_yy_xxxxxxz_1[j];

                    tg_xyy_xxxxxxyy_0[j] = pb_x * tg_yy_xxxxxxyy_0[j] + wp_x[j] * tg_yy_xxxxxxyy_1[j] + 3.0 * fl1_fxn * tg_yy_xxxxxyy_1[j];

                    tg_xyy_xxxxxxyz_0[j] = pb_x * tg_yy_xxxxxxyz_0[j] + wp_x[j] * tg_yy_xxxxxxyz_1[j] + 3.0 * fl1_fxn * tg_yy_xxxxxyz_1[j];

                    tg_xyy_xxxxxxzz_0[j] = pb_x * tg_yy_xxxxxxzz_0[j] + wp_x[j] * tg_yy_xxxxxxzz_1[j] + 3.0 * fl1_fxn * tg_yy_xxxxxzz_1[j];

                    tg_xyy_xxxxxyyy_0[j] = pb_x * tg_yy_xxxxxyyy_0[j] + wp_x[j] * tg_yy_xxxxxyyy_1[j] + 2.5 * fl1_fxn * tg_yy_xxxxyyy_1[j];

                    tg_xyy_xxxxxyyz_0[j] = pb_x * tg_yy_xxxxxyyz_0[j] + wp_x[j] * tg_yy_xxxxxyyz_1[j] + 2.5 * fl1_fxn * tg_yy_xxxxyyz_1[j];

                    tg_xyy_xxxxxyzz_0[j] = pb_x * tg_yy_xxxxxyzz_0[j] + wp_x[j] * tg_yy_xxxxxyzz_1[j] + 2.5 * fl1_fxn * tg_yy_xxxxyzz_1[j];

                    tg_xyy_xxxxxzzz_0[j] = pb_x * tg_yy_xxxxxzzz_0[j] + wp_x[j] * tg_yy_xxxxxzzz_1[j] + 2.5 * fl1_fxn * tg_yy_xxxxzzz_1[j];

                    tg_xyy_xxxxyyyy_0[j] = pb_x * tg_yy_xxxxyyyy_0[j] + wp_x[j] * tg_yy_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_yy_xxxyyyy_1[j];

                    tg_xyy_xxxxyyyz_0[j] = pb_x * tg_yy_xxxxyyyz_0[j] + wp_x[j] * tg_yy_xxxxyyyz_1[j] + 2.0 * fl1_fxn * tg_yy_xxxyyyz_1[j];

                    tg_xyy_xxxxyyzz_0[j] = pb_x * tg_yy_xxxxyyzz_0[j] + wp_x[j] * tg_yy_xxxxyyzz_1[j] + 2.0 * fl1_fxn * tg_yy_xxxyyzz_1[j];

                    tg_xyy_xxxxyzzz_0[j] = pb_x * tg_yy_xxxxyzzz_0[j] + wp_x[j] * tg_yy_xxxxyzzz_1[j] + 2.0 * fl1_fxn * tg_yy_xxxyzzz_1[j];

                    tg_xyy_xxxxzzzz_0[j] = pb_x * tg_yy_xxxxzzzz_0[j] + wp_x[j] * tg_yy_xxxxzzzz_1[j] + 2.0 * fl1_fxn * tg_yy_xxxzzzz_1[j];

                    tg_xyy_xxxyyyyy_0[j] = pb_x * tg_yy_xxxyyyyy_0[j] + wp_x[j] * tg_yy_xxxyyyyy_1[j] + 1.5 * fl1_fxn * tg_yy_xxyyyyy_1[j];

                    tg_xyy_xxxyyyyz_0[j] = pb_x * tg_yy_xxxyyyyz_0[j] + wp_x[j] * tg_yy_xxxyyyyz_1[j] + 1.5 * fl1_fxn * tg_yy_xxyyyyz_1[j];

                    tg_xyy_xxxyyyzz_0[j] = pb_x * tg_yy_xxxyyyzz_0[j] + wp_x[j] * tg_yy_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_yy_xxyyyzz_1[j];

                    tg_xyy_xxxyyzzz_0[j] = pb_x * tg_yy_xxxyyzzz_0[j] + wp_x[j] * tg_yy_xxxyyzzz_1[j] + 1.5 * fl1_fxn * tg_yy_xxyyzzz_1[j];

                    tg_xyy_xxxyzzzz_0[j] = pb_x * tg_yy_xxxyzzzz_0[j] + wp_x[j] * tg_yy_xxxyzzzz_1[j] + 1.5 * fl1_fxn * tg_yy_xxyzzzz_1[j];

                    tg_xyy_xxxzzzzz_0[j] = pb_x * tg_yy_xxxzzzzz_0[j] + wp_x[j] * tg_yy_xxxzzzzz_1[j] + 1.5 * fl1_fxn * tg_yy_xxzzzzz_1[j];

                    tg_xyy_xxyyyyyy_0[j] = pb_x * tg_yy_xxyyyyyy_0[j] + wp_x[j] * tg_yy_xxyyyyyy_1[j] + fl1_fxn * tg_yy_xyyyyyy_1[j];

                    tg_xyy_xxyyyyyz_0[j] = pb_x * tg_yy_xxyyyyyz_0[j] + wp_x[j] * tg_yy_xxyyyyyz_1[j] + fl1_fxn * tg_yy_xyyyyyz_1[j];

                    tg_xyy_xxyyyyzz_0[j] = pb_x * tg_yy_xxyyyyzz_0[j] + wp_x[j] * tg_yy_xxyyyyzz_1[j] + fl1_fxn * tg_yy_xyyyyzz_1[j];

                    tg_xyy_xxyyyzzz_0[j] = pb_x * tg_yy_xxyyyzzz_0[j] + wp_x[j] * tg_yy_xxyyyzzz_1[j] + fl1_fxn * tg_yy_xyyyzzz_1[j];

                    tg_xyy_xxyyzzzz_0[j] = pb_x * tg_yy_xxyyzzzz_0[j] + wp_x[j] * tg_yy_xxyyzzzz_1[j] + fl1_fxn * tg_yy_xyyzzzz_1[j];

                    tg_xyy_xxyzzzzz_0[j] = pb_x * tg_yy_xxyzzzzz_0[j] + wp_x[j] * tg_yy_xxyzzzzz_1[j] + fl1_fxn * tg_yy_xyzzzzz_1[j];

                    tg_xyy_xxzzzzzz_0[j] = pb_x * tg_yy_xxzzzzzz_0[j] + wp_x[j] * tg_yy_xxzzzzzz_1[j] + fl1_fxn * tg_yy_xzzzzzz_1[j];

                    tg_xyy_xyyyyyyy_0[j] = pb_x * tg_yy_xyyyyyyy_0[j] + wp_x[j] * tg_yy_xyyyyyyy_1[j] + 0.5 * fl1_fxn * tg_yy_yyyyyyy_1[j];

                    tg_xyy_xyyyyyyz_0[j] = pb_x * tg_yy_xyyyyyyz_0[j] + wp_x[j] * tg_yy_xyyyyyyz_1[j] + 0.5 * fl1_fxn * tg_yy_yyyyyyz_1[j];

                    tg_xyy_xyyyyyzz_0[j] = pb_x * tg_yy_xyyyyyzz_0[j] + wp_x[j] * tg_yy_xyyyyyzz_1[j] + 0.5 * fl1_fxn * tg_yy_yyyyyzz_1[j];

                    tg_xyy_xyyyyzzz_0[j] = pb_x * tg_yy_xyyyyzzz_0[j] + wp_x[j] * tg_yy_xyyyyzzz_1[j] + 0.5 * fl1_fxn * tg_yy_yyyyzzz_1[j];

                    tg_xyy_xyyyzzzz_0[j] = pb_x * tg_yy_xyyyzzzz_0[j] + wp_x[j] * tg_yy_xyyyzzzz_1[j] + 0.5 * fl1_fxn * tg_yy_yyyzzzz_1[j];

                    tg_xyy_xyyzzzzz_0[j] = pb_x * tg_yy_xyyzzzzz_0[j] + wp_x[j] * tg_yy_xyyzzzzz_1[j] + 0.5 * fl1_fxn * tg_yy_yyzzzzz_1[j];

                    tg_xyy_xyzzzzzz_0[j] = pb_x * tg_yy_xyzzzzzz_0[j] + wp_x[j] * tg_yy_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yy_yzzzzzz_1[j];

                    tg_xyy_xzzzzzzz_0[j] = pb_x * tg_yy_xzzzzzzz_0[j] + wp_x[j] * tg_yy_xzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yy_zzzzzzz_1[j];

                    tg_xyy_yyyyyyyy_0[j] = pb_x * tg_yy_yyyyyyyy_0[j] + wp_x[j] * tg_yy_yyyyyyyy_1[j];

                    tg_xyy_yyyyyyyz_0[j] = pb_x * tg_yy_yyyyyyyz_0[j] + wp_x[j] * tg_yy_yyyyyyyz_1[j];

                    tg_xyy_yyyyyyzz_0[j] = pb_x * tg_yy_yyyyyyzz_0[j] + wp_x[j] * tg_yy_yyyyyyzz_1[j];

                    tg_xyy_yyyyyzzz_0[j] = pb_x * tg_yy_yyyyyzzz_0[j] + wp_x[j] * tg_yy_yyyyyzzz_1[j];

                    tg_xyy_yyyyzzzz_0[j] = pb_x * tg_yy_yyyyzzzz_0[j] + wp_x[j] * tg_yy_yyyyzzzz_1[j];

                    tg_xyy_yyyzzzzz_0[j] = pb_x * tg_yy_yyyzzzzz_0[j] + wp_x[j] * tg_yy_yyyzzzzz_1[j];

                    tg_xyy_yyzzzzzz_0[j] = pb_x * tg_yy_yyzzzzzz_0[j] + wp_x[j] * tg_yy_yyzzzzzz_1[j];

                    tg_xyy_yzzzzzzz_0[j] = pb_x * tg_yy_yzzzzzzz_0[j] + wp_x[j] * tg_yy_yzzzzzzz_1[j];

                    tg_xyy_zzzzzzzz_0[j] = pb_x * tg_yy_zzzzzzzz_0[j] + wp_x[j] * tg_yy_zzzzzzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSFSL_180_270(      CMemBlock2D<double>& primBuffer,
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

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {3, -1, -1, -1},
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_3_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_3_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_2_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_2_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_2_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {7, -1, -1, -1}, 
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

                auto tg_yz_xxxxxxxx_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 180); 

                auto tg_yz_xxxxxxxy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 181); 

                auto tg_yz_xxxxxxxz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 182); 

                auto tg_yz_xxxxxxyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 183); 

                auto tg_yz_xxxxxxyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 184); 

                auto tg_yz_xxxxxxzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 185); 

                auto tg_yz_xxxxxyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 186); 

                auto tg_yz_xxxxxyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 187); 

                auto tg_yz_xxxxxyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 188); 

                auto tg_yz_xxxxxzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 189); 

                auto tg_yz_xxxxyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 190); 

                auto tg_yz_xxxxyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 191); 

                auto tg_yz_xxxxyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 192); 

                auto tg_yz_xxxxyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 193); 

                auto tg_yz_xxxxzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 194); 

                auto tg_yz_xxxyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 195); 

                auto tg_yz_xxxyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 196); 

                auto tg_yz_xxxyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 197); 

                auto tg_yz_xxxyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 198); 

                auto tg_yz_xxxyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 199); 

                auto tg_yz_xxxzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 200); 

                auto tg_yz_xxyyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 201); 

                auto tg_yz_xxyyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 202); 

                auto tg_yz_xxyyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 203); 

                auto tg_yz_xxyyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 204); 

                auto tg_yz_xxyyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 205); 

                auto tg_yz_xxyzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 206); 

                auto tg_yz_xxzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 207); 

                auto tg_yz_xyyyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 208); 

                auto tg_yz_xyyyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 209); 

                auto tg_yz_xyyyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 210); 

                auto tg_yz_xyyyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 211); 

                auto tg_yz_xyyyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 212); 

                auto tg_yz_xyyzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 213); 

                auto tg_yz_xyzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 214); 

                auto tg_yz_xzzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 215); 

                auto tg_yz_yyyyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 216); 

                auto tg_yz_yyyyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 217); 

                auto tg_yz_yyyyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 218); 

                auto tg_yz_yyyyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 219); 

                auto tg_yz_yyyyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 220); 

                auto tg_yz_yyyzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 221); 

                auto tg_yz_yyzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 222); 

                auto tg_yz_yzzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 223); 

                auto tg_yz_zzzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 224); 

                auto tg_zz_xxxxxxxx_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 225); 

                auto tg_zz_xxxxxxxy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 226); 

                auto tg_zz_xxxxxxxz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 227); 

                auto tg_zz_xxxxxxyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 228); 

                auto tg_zz_xxxxxxyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 229); 

                auto tg_zz_xxxxxxzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 230); 

                auto tg_zz_xxxxxyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 231); 

                auto tg_zz_xxxxxyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 232); 

                auto tg_zz_xxxxxyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 233); 

                auto tg_zz_xxxxxzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 234); 

                auto tg_zz_xxxxyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 235); 

                auto tg_zz_xxxxyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 236); 

                auto tg_zz_xxxxyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 237); 

                auto tg_zz_xxxxyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 238); 

                auto tg_zz_xxxxzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 239); 

                auto tg_zz_xxxyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 240); 

                auto tg_zz_xxxyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 241); 

                auto tg_zz_xxxyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 242); 

                auto tg_zz_xxxyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 243); 

                auto tg_zz_xxxyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 244); 

                auto tg_zz_xxxzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 245); 

                auto tg_zz_xxyyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 246); 

                auto tg_zz_xxyyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 247); 

                auto tg_zz_xxyyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 248); 

                auto tg_zz_xxyyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 249); 

                auto tg_zz_xxyyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 250); 

                auto tg_zz_xxyzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 251); 

                auto tg_zz_xxzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 252); 

                auto tg_zz_xyyyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 253); 

                auto tg_zz_xyyyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 254); 

                auto tg_zz_xyyyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 255); 

                auto tg_zz_xyyyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 256); 

                auto tg_zz_xyyyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 257); 

                auto tg_zz_xyyzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 258); 

                auto tg_zz_xyzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 259); 

                auto tg_zz_xzzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 260); 

                auto tg_zz_yyyyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 261); 

                auto tg_zz_yyyyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 262); 

                auto tg_zz_yyyyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 263); 

                auto tg_zz_yyyyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 264); 

                auto tg_zz_yyyyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 265); 

                auto tg_zz_yyyzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 266); 

                auto tg_zz_yyzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 267); 

                auto tg_zz_yzzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 268); 

                auto tg_zz_zzzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 269); 

                auto tg_yz_xxxxxxxx_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 180); 

                auto tg_yz_xxxxxxxy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 181); 

                auto tg_yz_xxxxxxxz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 182); 

                auto tg_yz_xxxxxxyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 183); 

                auto tg_yz_xxxxxxyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 184); 

                auto tg_yz_xxxxxxzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 185); 

                auto tg_yz_xxxxxyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 186); 

                auto tg_yz_xxxxxyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 187); 

                auto tg_yz_xxxxxyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 188); 

                auto tg_yz_xxxxxzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 189); 

                auto tg_yz_xxxxyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 190); 

                auto tg_yz_xxxxyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 191); 

                auto tg_yz_xxxxyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 192); 

                auto tg_yz_xxxxyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 193); 

                auto tg_yz_xxxxzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 194); 

                auto tg_yz_xxxyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 195); 

                auto tg_yz_xxxyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 196); 

                auto tg_yz_xxxyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 197); 

                auto tg_yz_xxxyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 198); 

                auto tg_yz_xxxyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 199); 

                auto tg_yz_xxxzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 200); 

                auto tg_yz_xxyyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 201); 

                auto tg_yz_xxyyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 202); 

                auto tg_yz_xxyyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 203); 

                auto tg_yz_xxyyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 204); 

                auto tg_yz_xxyyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 205); 

                auto tg_yz_xxyzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 206); 

                auto tg_yz_xxzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 207); 

                auto tg_yz_xyyyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 208); 

                auto tg_yz_xyyyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 209); 

                auto tg_yz_xyyyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 210); 

                auto tg_yz_xyyyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 211); 

                auto tg_yz_xyyyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 212); 

                auto tg_yz_xyyzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 213); 

                auto tg_yz_xyzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 214); 

                auto tg_yz_xzzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 215); 

                auto tg_yz_yyyyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 216); 

                auto tg_yz_yyyyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 217); 

                auto tg_yz_yyyyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 218); 

                auto tg_yz_yyyyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 219); 

                auto tg_yz_yyyyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 220); 

                auto tg_yz_yyyzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 221); 

                auto tg_yz_yyzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 222); 

                auto tg_yz_yzzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 223); 

                auto tg_yz_zzzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 224); 

                auto tg_zz_xxxxxxxx_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 225); 

                auto tg_zz_xxxxxxxy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 226); 

                auto tg_zz_xxxxxxxz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 227); 

                auto tg_zz_xxxxxxyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 228); 

                auto tg_zz_xxxxxxyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 229); 

                auto tg_zz_xxxxxxzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 230); 

                auto tg_zz_xxxxxyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 231); 

                auto tg_zz_xxxxxyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 232); 

                auto tg_zz_xxxxxyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 233); 

                auto tg_zz_xxxxxzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 234); 

                auto tg_zz_xxxxyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 235); 

                auto tg_zz_xxxxyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 236); 

                auto tg_zz_xxxxyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 237); 

                auto tg_zz_xxxxyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 238); 

                auto tg_zz_xxxxzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 239); 

                auto tg_zz_xxxyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 240); 

                auto tg_zz_xxxyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 241); 

                auto tg_zz_xxxyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 242); 

                auto tg_zz_xxxyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 243); 

                auto tg_zz_xxxyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 244); 

                auto tg_zz_xxxzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 245); 

                auto tg_zz_xxyyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 246); 

                auto tg_zz_xxyyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 247); 

                auto tg_zz_xxyyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 248); 

                auto tg_zz_xxyyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 249); 

                auto tg_zz_xxyyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 250); 

                auto tg_zz_xxyzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 251); 

                auto tg_zz_xxzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 252); 

                auto tg_zz_xyyyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 253); 

                auto tg_zz_xyyyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 254); 

                auto tg_zz_xyyyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 255); 

                auto tg_zz_xyyyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 256); 

                auto tg_zz_xyyyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 257); 

                auto tg_zz_xyyzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 258); 

                auto tg_zz_xyzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 259); 

                auto tg_zz_xzzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 260); 

                auto tg_zz_yyyyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 261); 

                auto tg_zz_yyyyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 262); 

                auto tg_zz_yyyyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 263); 

                auto tg_zz_yyyyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 264); 

                auto tg_zz_yyyyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 265); 

                auto tg_zz_yyyzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 266); 

                auto tg_zz_yyzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 267); 

                auto tg_zz_yzzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 268); 

                auto tg_zz_zzzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 269); 

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

                // set up pointers to integrals

                auto tg_xyz_xxxxxxxx_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 180); 

                auto tg_xyz_xxxxxxxy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 181); 

                auto tg_xyz_xxxxxxxz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 182); 

                auto tg_xyz_xxxxxxyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 183); 

                auto tg_xyz_xxxxxxyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 184); 

                auto tg_xyz_xxxxxxzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 185); 

                auto tg_xyz_xxxxxyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 186); 

                auto tg_xyz_xxxxxyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 187); 

                auto tg_xyz_xxxxxyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 188); 

                auto tg_xyz_xxxxxzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 189); 

                auto tg_xyz_xxxxyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 190); 

                auto tg_xyz_xxxxyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 191); 

                auto tg_xyz_xxxxyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 192); 

                auto tg_xyz_xxxxyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 193); 

                auto tg_xyz_xxxxzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 194); 

                auto tg_xyz_xxxyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 195); 

                auto tg_xyz_xxxyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 196); 

                auto tg_xyz_xxxyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 197); 

                auto tg_xyz_xxxyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 198); 

                auto tg_xyz_xxxyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 199); 

                auto tg_xyz_xxxzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 200); 

                auto tg_xyz_xxyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 201); 

                auto tg_xyz_xxyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 202); 

                auto tg_xyz_xxyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 203); 

                auto tg_xyz_xxyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 204); 

                auto tg_xyz_xxyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 205); 

                auto tg_xyz_xxyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 206); 

                auto tg_xyz_xxzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 207); 

                auto tg_xyz_xyyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 208); 

                auto tg_xyz_xyyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 209); 

                auto tg_xyz_xyyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 210); 

                auto tg_xyz_xyyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 211); 

                auto tg_xyz_xyyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 212); 

                auto tg_xyz_xyyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 213); 

                auto tg_xyz_xyzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 214); 

                auto tg_xyz_xzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 215); 

                auto tg_xyz_yyyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 216); 

                auto tg_xyz_yyyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 217); 

                auto tg_xyz_yyyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 218); 

                auto tg_xyz_yyyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 219); 

                auto tg_xyz_yyyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 220); 

                auto tg_xyz_yyyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 221); 

                auto tg_xyz_yyzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 222); 

                auto tg_xyz_yzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 223); 

                auto tg_xyz_zzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 224); 

                auto tg_xzz_xxxxxxxx_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 225); 

                auto tg_xzz_xxxxxxxy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 226); 

                auto tg_xzz_xxxxxxxz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 227); 

                auto tg_xzz_xxxxxxyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 228); 

                auto tg_xzz_xxxxxxyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 229); 

                auto tg_xzz_xxxxxxzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 230); 

                auto tg_xzz_xxxxxyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 231); 

                auto tg_xzz_xxxxxyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 232); 

                auto tg_xzz_xxxxxyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 233); 

                auto tg_xzz_xxxxxzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 234); 

                auto tg_xzz_xxxxyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 235); 

                auto tg_xzz_xxxxyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 236); 

                auto tg_xzz_xxxxyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 237); 

                auto tg_xzz_xxxxyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 238); 

                auto tg_xzz_xxxxzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 239); 

                auto tg_xzz_xxxyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 240); 

                auto tg_xzz_xxxyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 241); 

                auto tg_xzz_xxxyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 242); 

                auto tg_xzz_xxxyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 243); 

                auto tg_xzz_xxxyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 244); 

                auto tg_xzz_xxxzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 245); 

                auto tg_xzz_xxyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 246); 

                auto tg_xzz_xxyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 247); 

                auto tg_xzz_xxyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 248); 

                auto tg_xzz_xxyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 249); 

                auto tg_xzz_xxyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 250); 

                auto tg_xzz_xxyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 251); 

                auto tg_xzz_xxzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 252); 

                auto tg_xzz_xyyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 253); 

                auto tg_xzz_xyyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 254); 

                auto tg_xzz_xyyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 255); 

                auto tg_xzz_xyyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 256); 

                auto tg_xzz_xyyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 257); 

                auto tg_xzz_xyyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 258); 

                auto tg_xzz_xyzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 259); 

                auto tg_xzz_xzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 260); 

                auto tg_xzz_yyyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 261); 

                auto tg_xzz_yyyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 262); 

                auto tg_xzz_yyyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 263); 

                auto tg_xzz_yyyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 264); 

                auto tg_xzz_yyyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 265); 

                auto tg_xzz_yyyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 266); 

                auto tg_xzz_yyzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 267); 

                auto tg_xzz_yzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 268); 

                auto tg_xzz_zzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 269); 

                // Batch of Integrals (180,270)

                #pragma omp simd aligned(fxn, tg_xyz_xxxxxxxx_0, tg_xyz_xxxxxxxy_0, tg_xyz_xxxxxxxz_0, \
                                         tg_xyz_xxxxxxyy_0, tg_xyz_xxxxxxyz_0, tg_xyz_xxxxxxzz_0, tg_xyz_xxxxxyyy_0, \
                                         tg_xyz_xxxxxyyz_0, tg_xyz_xxxxxyzz_0, tg_xyz_xxxxxzzz_0, tg_xyz_xxxxyyyy_0, \
                                         tg_xyz_xxxxyyyz_0, tg_xyz_xxxxyyzz_0, tg_xyz_xxxxyzzz_0, tg_xyz_xxxxzzzz_0, \
                                         tg_xyz_xxxyyyyy_0, tg_xyz_xxxyyyyz_0, tg_xyz_xxxyyyzz_0, tg_xyz_xxxyyzzz_0, \
                                         tg_xyz_xxxyzzzz_0, tg_xyz_xxxzzzzz_0, tg_xyz_xxyyyyyy_0, tg_xyz_xxyyyyyz_0, \
                                         tg_xyz_xxyyyyzz_0, tg_xyz_xxyyyzzz_0, tg_xyz_xxyyzzzz_0, tg_xyz_xxyzzzzz_0, \
                                         tg_xyz_xxzzzzzz_0, tg_xyz_xyyyyyyy_0, tg_xyz_xyyyyyyz_0, tg_xyz_xyyyyyzz_0, \
                                         tg_xyz_xyyyyzzz_0, tg_xyz_xyyyzzzz_0, tg_xyz_xyyzzzzz_0, tg_xyz_xyzzzzzz_0, \
                                         tg_xyz_xzzzzzzz_0, tg_xyz_yyyyyyyy_0, tg_xyz_yyyyyyyz_0, tg_xyz_yyyyyyzz_0, \
                                         tg_xyz_yyyyyzzz_0, tg_xyz_yyyyzzzz_0, tg_xyz_yyyzzzzz_0, tg_xyz_yyzzzzzz_0, \
                                         tg_xyz_yzzzzzzz_0, tg_xyz_zzzzzzzz_0, tg_xzz_xxxxxxxx_0, tg_xzz_xxxxxxxy_0, \
                                         tg_xzz_xxxxxxxz_0, tg_xzz_xxxxxxyy_0, tg_xzz_xxxxxxyz_0, tg_xzz_xxxxxxzz_0, \
                                         tg_xzz_xxxxxyyy_0, tg_xzz_xxxxxyyz_0, tg_xzz_xxxxxyzz_0, tg_xzz_xxxxxzzz_0, \
                                         tg_xzz_xxxxyyyy_0, tg_xzz_xxxxyyyz_0, tg_xzz_xxxxyyzz_0, tg_xzz_xxxxyzzz_0, \
                                         tg_xzz_xxxxzzzz_0, tg_xzz_xxxyyyyy_0, tg_xzz_xxxyyyyz_0, tg_xzz_xxxyyyzz_0, \
                                         tg_xzz_xxxyyzzz_0, tg_xzz_xxxyzzzz_0, tg_xzz_xxxzzzzz_0, tg_xzz_xxyyyyyy_0, \
                                         tg_xzz_xxyyyyyz_0, tg_xzz_xxyyyyzz_0, tg_xzz_xxyyyzzz_0, tg_xzz_xxyyzzzz_0, \
                                         tg_xzz_xxyzzzzz_0, tg_xzz_xxzzzzzz_0, tg_xzz_xyyyyyyy_0, tg_xzz_xyyyyyyz_0, \
                                         tg_xzz_xyyyyyzz_0, tg_xzz_xyyyyzzz_0, tg_xzz_xyyyzzzz_0, tg_xzz_xyyzzzzz_0, \
                                         tg_xzz_xyzzzzzz_0, tg_xzz_xzzzzzzz_0, tg_xzz_yyyyyyyy_0, tg_xzz_yyyyyyyz_0, \
                                         tg_xzz_yyyyyyzz_0, tg_xzz_yyyyyzzz_0, tg_xzz_yyyyzzzz_0, tg_xzz_yyyzzzzz_0, \
                                         tg_xzz_yyzzzzzz_0, tg_xzz_yzzzzzzz_0, tg_xzz_zzzzzzzz_0, tg_yz_xxxxxxx_1, \
                                         tg_yz_xxxxxxxx_0, tg_yz_xxxxxxxx_1, tg_yz_xxxxxxxy_0, tg_yz_xxxxxxxy_1, \
                                         tg_yz_xxxxxxxz_0, tg_yz_xxxxxxxz_1, tg_yz_xxxxxxy_1, tg_yz_xxxxxxyy_0, \
                                         tg_yz_xxxxxxyy_1, tg_yz_xxxxxxyz_0, tg_yz_xxxxxxyz_1, tg_yz_xxxxxxz_1, \
                                         tg_yz_xxxxxxzz_0, tg_yz_xxxxxxzz_1, tg_yz_xxxxxyy_1, tg_yz_xxxxxyyy_0, \
                                         tg_yz_xxxxxyyy_1, tg_yz_xxxxxyyz_0, tg_yz_xxxxxyyz_1, tg_yz_xxxxxyz_1, \
                                         tg_yz_xxxxxyzz_0, tg_yz_xxxxxyzz_1, tg_yz_xxxxxzz_1, tg_yz_xxxxxzzz_0, \
                                         tg_yz_xxxxxzzz_1, tg_yz_xxxxyyy_1, tg_yz_xxxxyyyy_0, tg_yz_xxxxyyyy_1, \
                                         tg_yz_xxxxyyyz_0, tg_yz_xxxxyyyz_1, tg_yz_xxxxyyz_1, tg_yz_xxxxyyzz_0, \
                                         tg_yz_xxxxyyzz_1, tg_yz_xxxxyzz_1, tg_yz_xxxxyzzz_0, tg_yz_xxxxyzzz_1, \
                                         tg_yz_xxxxzzz_1, tg_yz_xxxxzzzz_0, tg_yz_xxxxzzzz_1, tg_yz_xxxyyyy_1, \
                                         tg_yz_xxxyyyyy_0, tg_yz_xxxyyyyy_1, tg_yz_xxxyyyyz_0, tg_yz_xxxyyyyz_1, \
                                         tg_yz_xxxyyyz_1, tg_yz_xxxyyyzz_0, tg_yz_xxxyyyzz_1, tg_yz_xxxyyzz_1, \
                                         tg_yz_xxxyyzzz_0, tg_yz_xxxyyzzz_1, tg_yz_xxxyzzz_1, tg_yz_xxxyzzzz_0, \
                                         tg_yz_xxxyzzzz_1, tg_yz_xxxzzzz_1, tg_yz_xxxzzzzz_0, tg_yz_xxxzzzzz_1, \
                                         tg_yz_xxyyyyy_1, tg_yz_xxyyyyyy_0, tg_yz_xxyyyyyy_1, tg_yz_xxyyyyyz_0, \
                                         tg_yz_xxyyyyyz_1, tg_yz_xxyyyyz_1, tg_yz_xxyyyyzz_0, tg_yz_xxyyyyzz_1, \
                                         tg_yz_xxyyyzz_1, tg_yz_xxyyyzzz_0, tg_yz_xxyyyzzz_1, tg_yz_xxyyzzz_1, \
                                         tg_yz_xxyyzzzz_0, tg_yz_xxyyzzzz_1, tg_yz_xxyzzzz_1, tg_yz_xxyzzzzz_0, \
                                         tg_yz_xxyzzzzz_1, tg_yz_xxzzzzz_1, tg_yz_xxzzzzzz_0, tg_yz_xxzzzzzz_1, \
                                         tg_yz_xyyyyyy_1, tg_yz_xyyyyyyy_0, tg_yz_xyyyyyyy_1, tg_yz_xyyyyyyz_0, \
                                         tg_yz_xyyyyyyz_1, tg_yz_xyyyyyz_1, tg_yz_xyyyyyzz_0, tg_yz_xyyyyyzz_1, \
                                         tg_yz_xyyyyzz_1, tg_yz_xyyyyzzz_0, tg_yz_xyyyyzzz_1, tg_yz_xyyyzzz_1, \
                                         tg_yz_xyyyzzzz_0, tg_yz_xyyyzzzz_1, tg_yz_xyyzzzz_1, tg_yz_xyyzzzzz_0, \
                                         tg_yz_xyyzzzzz_1, tg_yz_xyzzzzz_1, tg_yz_xyzzzzzz_0, tg_yz_xyzzzzzz_1, \
                                         tg_yz_xzzzzzz_1, tg_yz_xzzzzzzz_0, tg_yz_xzzzzzzz_1, tg_yz_yyyyyyy_1, \
                                         tg_yz_yyyyyyyy_0, tg_yz_yyyyyyyy_1, tg_yz_yyyyyyyz_0, tg_yz_yyyyyyyz_1, \
                                         tg_yz_yyyyyyz_1, tg_yz_yyyyyyzz_0, tg_yz_yyyyyyzz_1, tg_yz_yyyyyzz_1, \
                                         tg_yz_yyyyyzzz_0, tg_yz_yyyyyzzz_1, tg_yz_yyyyzzz_1, tg_yz_yyyyzzzz_0, \
                                         tg_yz_yyyyzzzz_1, tg_yz_yyyzzzz_1, tg_yz_yyyzzzzz_0, tg_yz_yyyzzzzz_1, \
                                         tg_yz_yyzzzzz_1, tg_yz_yyzzzzzz_0, tg_yz_yyzzzzzz_1, tg_yz_yzzzzzz_1, \
                                         tg_yz_yzzzzzzz_0, tg_yz_yzzzzzzz_1, tg_yz_zzzzzzz_1, tg_yz_zzzzzzzz_0, \
                                         tg_yz_zzzzzzzz_1, tg_zz_xxxxxxx_1, tg_zz_xxxxxxxx_0, tg_zz_xxxxxxxx_1, \
                                         tg_zz_xxxxxxxy_0, tg_zz_xxxxxxxy_1, tg_zz_xxxxxxxz_0, tg_zz_xxxxxxxz_1, \
                                         tg_zz_xxxxxxy_1, tg_zz_xxxxxxyy_0, tg_zz_xxxxxxyy_1, tg_zz_xxxxxxyz_0, \
                                         tg_zz_xxxxxxyz_1, tg_zz_xxxxxxz_1, tg_zz_xxxxxxzz_0, tg_zz_xxxxxxzz_1, \
                                         tg_zz_xxxxxyy_1, tg_zz_xxxxxyyy_0, tg_zz_xxxxxyyy_1, tg_zz_xxxxxyyz_0, \
                                         tg_zz_xxxxxyyz_1, tg_zz_xxxxxyz_1, tg_zz_xxxxxyzz_0, tg_zz_xxxxxyzz_1, \
                                         tg_zz_xxxxxzz_1, tg_zz_xxxxxzzz_0, tg_zz_xxxxxzzz_1, tg_zz_xxxxyyy_1, \
                                         tg_zz_xxxxyyyy_0, tg_zz_xxxxyyyy_1, tg_zz_xxxxyyyz_0, tg_zz_xxxxyyyz_1, \
                                         tg_zz_xxxxyyz_1, tg_zz_xxxxyyzz_0, tg_zz_xxxxyyzz_1, tg_zz_xxxxyzz_1, \
                                         tg_zz_xxxxyzzz_0, tg_zz_xxxxyzzz_1, tg_zz_xxxxzzz_1, tg_zz_xxxxzzzz_0, \
                                         tg_zz_xxxxzzzz_1, tg_zz_xxxyyyy_1, tg_zz_xxxyyyyy_0, tg_zz_xxxyyyyy_1, \
                                         tg_zz_xxxyyyyz_0, tg_zz_xxxyyyyz_1, tg_zz_xxxyyyz_1, tg_zz_xxxyyyzz_0, \
                                         tg_zz_xxxyyyzz_1, tg_zz_xxxyyzz_1, tg_zz_xxxyyzzz_0, tg_zz_xxxyyzzz_1, \
                                         tg_zz_xxxyzzz_1, tg_zz_xxxyzzzz_0, tg_zz_xxxyzzzz_1, tg_zz_xxxzzzz_1, \
                                         tg_zz_xxxzzzzz_0, tg_zz_xxxzzzzz_1, tg_zz_xxyyyyy_1, tg_zz_xxyyyyyy_0, \
                                         tg_zz_xxyyyyyy_1, tg_zz_xxyyyyyz_0, tg_zz_xxyyyyyz_1, tg_zz_xxyyyyz_1, \
                                         tg_zz_xxyyyyzz_0, tg_zz_xxyyyyzz_1, tg_zz_xxyyyzz_1, tg_zz_xxyyyzzz_0, \
                                         tg_zz_xxyyyzzz_1, tg_zz_xxyyzzz_1, tg_zz_xxyyzzzz_0, tg_zz_xxyyzzzz_1, \
                                         tg_zz_xxyzzzz_1, tg_zz_xxyzzzzz_0, tg_zz_xxyzzzzz_1, tg_zz_xxzzzzz_1, \
                                         tg_zz_xxzzzzzz_0, tg_zz_xxzzzzzz_1, tg_zz_xyyyyyy_1, tg_zz_xyyyyyyy_0, \
                                         tg_zz_xyyyyyyy_1, tg_zz_xyyyyyyz_0, tg_zz_xyyyyyyz_1, tg_zz_xyyyyyz_1, \
                                         tg_zz_xyyyyyzz_0, tg_zz_xyyyyyzz_1, tg_zz_xyyyyzz_1, tg_zz_xyyyyzzz_0, \
                                         tg_zz_xyyyyzzz_1, tg_zz_xyyyzzz_1, tg_zz_xyyyzzzz_0, tg_zz_xyyyzzzz_1, \
                                         tg_zz_xyyzzzz_1, tg_zz_xyyzzzzz_0, tg_zz_xyyzzzzz_1, tg_zz_xyzzzzz_1, \
                                         tg_zz_xyzzzzzz_0, tg_zz_xyzzzzzz_1, tg_zz_xzzzzzz_1, tg_zz_xzzzzzzz_0, \
                                         tg_zz_xzzzzzzz_1, tg_zz_yyyyyyy_1, tg_zz_yyyyyyyy_0, tg_zz_yyyyyyyy_1, \
                                         tg_zz_yyyyyyyz_0, tg_zz_yyyyyyyz_1, tg_zz_yyyyyyz_1, tg_zz_yyyyyyzz_0, \
                                         tg_zz_yyyyyyzz_1, tg_zz_yyyyyzz_1, tg_zz_yyyyyzzz_0, tg_zz_yyyyyzzz_1, \
                                         tg_zz_yyyyzzz_1, tg_zz_yyyyzzzz_0, tg_zz_yyyyzzzz_1, tg_zz_yyyzzzz_1, \
                                         tg_zz_yyyzzzzz_0, tg_zz_yyyzzzzz_1, tg_zz_yyzzzzz_1, tg_zz_yyzzzzzz_0, \
                                         tg_zz_yyzzzzzz_1, tg_zz_yzzzzzz_1, tg_zz_yzzzzzzz_0, tg_zz_yzzzzzzz_1, \
                                         tg_zz_zzzzzzz_1, tg_zz_zzzzzzzz_0, tg_zz_zzzzzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fxn = fxn[j];

                    tg_xyz_xxxxxxxx_0[j] = pb_x * tg_yz_xxxxxxxx_0[j] + wp_x[j] * tg_yz_xxxxxxxx_1[j] + 4.0 * fl1_fxn * tg_yz_xxxxxxx_1[j];

                    tg_xyz_xxxxxxxy_0[j] = pb_x * tg_yz_xxxxxxxy_0[j] + wp_x[j] * tg_yz_xxxxxxxy_1[j] + 3.5 * fl1_fxn * tg_yz_xxxxxxy_1[j];

                    tg_xyz_xxxxxxxz_0[j] = pb_x * tg_yz_xxxxxxxz_0[j] + wp_x[j] * tg_yz_xxxxxxxz_1[j] + 3.5 * fl1_fxn * tg_yz_xxxxxxz_1[j];

                    tg_xyz_xxxxxxyy_0[j] = pb_x * tg_yz_xxxxxxyy_0[j] + wp_x[j] * tg_yz_xxxxxxyy_1[j] + 3.0 * fl1_fxn * tg_yz_xxxxxyy_1[j];

                    tg_xyz_xxxxxxyz_0[j] = pb_x * tg_yz_xxxxxxyz_0[j] + wp_x[j] * tg_yz_xxxxxxyz_1[j] + 3.0 * fl1_fxn * tg_yz_xxxxxyz_1[j];

                    tg_xyz_xxxxxxzz_0[j] = pb_x * tg_yz_xxxxxxzz_0[j] + wp_x[j] * tg_yz_xxxxxxzz_1[j] + 3.0 * fl1_fxn * tg_yz_xxxxxzz_1[j];

                    tg_xyz_xxxxxyyy_0[j] = pb_x * tg_yz_xxxxxyyy_0[j] + wp_x[j] * tg_yz_xxxxxyyy_1[j] + 2.5 * fl1_fxn * tg_yz_xxxxyyy_1[j];

                    tg_xyz_xxxxxyyz_0[j] = pb_x * tg_yz_xxxxxyyz_0[j] + wp_x[j] * tg_yz_xxxxxyyz_1[j] + 2.5 * fl1_fxn * tg_yz_xxxxyyz_1[j];

                    tg_xyz_xxxxxyzz_0[j] = pb_x * tg_yz_xxxxxyzz_0[j] + wp_x[j] * tg_yz_xxxxxyzz_1[j] + 2.5 * fl1_fxn * tg_yz_xxxxyzz_1[j];

                    tg_xyz_xxxxxzzz_0[j] = pb_x * tg_yz_xxxxxzzz_0[j] + wp_x[j] * tg_yz_xxxxxzzz_1[j] + 2.5 * fl1_fxn * tg_yz_xxxxzzz_1[j];

                    tg_xyz_xxxxyyyy_0[j] = pb_x * tg_yz_xxxxyyyy_0[j] + wp_x[j] * tg_yz_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_yz_xxxyyyy_1[j];

                    tg_xyz_xxxxyyyz_0[j] = pb_x * tg_yz_xxxxyyyz_0[j] + wp_x[j] * tg_yz_xxxxyyyz_1[j] + 2.0 * fl1_fxn * tg_yz_xxxyyyz_1[j];

                    tg_xyz_xxxxyyzz_0[j] = pb_x * tg_yz_xxxxyyzz_0[j] + wp_x[j] * tg_yz_xxxxyyzz_1[j] + 2.0 * fl1_fxn * tg_yz_xxxyyzz_1[j];

                    tg_xyz_xxxxyzzz_0[j] = pb_x * tg_yz_xxxxyzzz_0[j] + wp_x[j] * tg_yz_xxxxyzzz_1[j] + 2.0 * fl1_fxn * tg_yz_xxxyzzz_1[j];

                    tg_xyz_xxxxzzzz_0[j] = pb_x * tg_yz_xxxxzzzz_0[j] + wp_x[j] * tg_yz_xxxxzzzz_1[j] + 2.0 * fl1_fxn * tg_yz_xxxzzzz_1[j];

                    tg_xyz_xxxyyyyy_0[j] = pb_x * tg_yz_xxxyyyyy_0[j] + wp_x[j] * tg_yz_xxxyyyyy_1[j] + 1.5 * fl1_fxn * tg_yz_xxyyyyy_1[j];

                    tg_xyz_xxxyyyyz_0[j] = pb_x * tg_yz_xxxyyyyz_0[j] + wp_x[j] * tg_yz_xxxyyyyz_1[j] + 1.5 * fl1_fxn * tg_yz_xxyyyyz_1[j];

                    tg_xyz_xxxyyyzz_0[j] = pb_x * tg_yz_xxxyyyzz_0[j] + wp_x[j] * tg_yz_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_yz_xxyyyzz_1[j];

                    tg_xyz_xxxyyzzz_0[j] = pb_x * tg_yz_xxxyyzzz_0[j] + wp_x[j] * tg_yz_xxxyyzzz_1[j] + 1.5 * fl1_fxn * tg_yz_xxyyzzz_1[j];

                    tg_xyz_xxxyzzzz_0[j] = pb_x * tg_yz_xxxyzzzz_0[j] + wp_x[j] * tg_yz_xxxyzzzz_1[j] + 1.5 * fl1_fxn * tg_yz_xxyzzzz_1[j];

                    tg_xyz_xxxzzzzz_0[j] = pb_x * tg_yz_xxxzzzzz_0[j] + wp_x[j] * tg_yz_xxxzzzzz_1[j] + 1.5 * fl1_fxn * tg_yz_xxzzzzz_1[j];

                    tg_xyz_xxyyyyyy_0[j] = pb_x * tg_yz_xxyyyyyy_0[j] + wp_x[j] * tg_yz_xxyyyyyy_1[j] + fl1_fxn * tg_yz_xyyyyyy_1[j];

                    tg_xyz_xxyyyyyz_0[j] = pb_x * tg_yz_xxyyyyyz_0[j] + wp_x[j] * tg_yz_xxyyyyyz_1[j] + fl1_fxn * tg_yz_xyyyyyz_1[j];

                    tg_xyz_xxyyyyzz_0[j] = pb_x * tg_yz_xxyyyyzz_0[j] + wp_x[j] * tg_yz_xxyyyyzz_1[j] + fl1_fxn * tg_yz_xyyyyzz_1[j];

                    tg_xyz_xxyyyzzz_0[j] = pb_x * tg_yz_xxyyyzzz_0[j] + wp_x[j] * tg_yz_xxyyyzzz_1[j] + fl1_fxn * tg_yz_xyyyzzz_1[j];

                    tg_xyz_xxyyzzzz_0[j] = pb_x * tg_yz_xxyyzzzz_0[j] + wp_x[j] * tg_yz_xxyyzzzz_1[j] + fl1_fxn * tg_yz_xyyzzzz_1[j];

                    tg_xyz_xxyzzzzz_0[j] = pb_x * tg_yz_xxyzzzzz_0[j] + wp_x[j] * tg_yz_xxyzzzzz_1[j] + fl1_fxn * tg_yz_xyzzzzz_1[j];

                    tg_xyz_xxzzzzzz_0[j] = pb_x * tg_yz_xxzzzzzz_0[j] + wp_x[j] * tg_yz_xxzzzzzz_1[j] + fl1_fxn * tg_yz_xzzzzzz_1[j];

                    tg_xyz_xyyyyyyy_0[j] = pb_x * tg_yz_xyyyyyyy_0[j] + wp_x[j] * tg_yz_xyyyyyyy_1[j] + 0.5 * fl1_fxn * tg_yz_yyyyyyy_1[j];

                    tg_xyz_xyyyyyyz_0[j] = pb_x * tg_yz_xyyyyyyz_0[j] + wp_x[j] * tg_yz_xyyyyyyz_1[j] + 0.5 * fl1_fxn * tg_yz_yyyyyyz_1[j];

                    tg_xyz_xyyyyyzz_0[j] = pb_x * tg_yz_xyyyyyzz_0[j] + wp_x[j] * tg_yz_xyyyyyzz_1[j] + 0.5 * fl1_fxn * tg_yz_yyyyyzz_1[j];

                    tg_xyz_xyyyyzzz_0[j] = pb_x * tg_yz_xyyyyzzz_0[j] + wp_x[j] * tg_yz_xyyyyzzz_1[j] + 0.5 * fl1_fxn * tg_yz_yyyyzzz_1[j];

                    tg_xyz_xyyyzzzz_0[j] = pb_x * tg_yz_xyyyzzzz_0[j] + wp_x[j] * tg_yz_xyyyzzzz_1[j] + 0.5 * fl1_fxn * tg_yz_yyyzzzz_1[j];

                    tg_xyz_xyyzzzzz_0[j] = pb_x * tg_yz_xyyzzzzz_0[j] + wp_x[j] * tg_yz_xyyzzzzz_1[j] + 0.5 * fl1_fxn * tg_yz_yyzzzzz_1[j];

                    tg_xyz_xyzzzzzz_0[j] = pb_x * tg_yz_xyzzzzzz_0[j] + wp_x[j] * tg_yz_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yz_yzzzzzz_1[j];

                    tg_xyz_xzzzzzzz_0[j] = pb_x * tg_yz_xzzzzzzz_0[j] + wp_x[j] * tg_yz_xzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yz_zzzzzzz_1[j];

                    tg_xyz_yyyyyyyy_0[j] = pb_x * tg_yz_yyyyyyyy_0[j] + wp_x[j] * tg_yz_yyyyyyyy_1[j];

                    tg_xyz_yyyyyyyz_0[j] = pb_x * tg_yz_yyyyyyyz_0[j] + wp_x[j] * tg_yz_yyyyyyyz_1[j];

                    tg_xyz_yyyyyyzz_0[j] = pb_x * tg_yz_yyyyyyzz_0[j] + wp_x[j] * tg_yz_yyyyyyzz_1[j];

                    tg_xyz_yyyyyzzz_0[j] = pb_x * tg_yz_yyyyyzzz_0[j] + wp_x[j] * tg_yz_yyyyyzzz_1[j];

                    tg_xyz_yyyyzzzz_0[j] = pb_x * tg_yz_yyyyzzzz_0[j] + wp_x[j] * tg_yz_yyyyzzzz_1[j];

                    tg_xyz_yyyzzzzz_0[j] = pb_x * tg_yz_yyyzzzzz_0[j] + wp_x[j] * tg_yz_yyyzzzzz_1[j];

                    tg_xyz_yyzzzzzz_0[j] = pb_x * tg_yz_yyzzzzzz_0[j] + wp_x[j] * tg_yz_yyzzzzzz_1[j];

                    tg_xyz_yzzzzzzz_0[j] = pb_x * tg_yz_yzzzzzzz_0[j] + wp_x[j] * tg_yz_yzzzzzzz_1[j];

                    tg_xyz_zzzzzzzz_0[j] = pb_x * tg_yz_zzzzzzzz_0[j] + wp_x[j] * tg_yz_zzzzzzzz_1[j];

                    tg_xzz_xxxxxxxx_0[j] = pb_x * tg_zz_xxxxxxxx_0[j] + wp_x[j] * tg_zz_xxxxxxxx_1[j] + 4.0 * fl1_fxn * tg_zz_xxxxxxx_1[j];

                    tg_xzz_xxxxxxxy_0[j] = pb_x * tg_zz_xxxxxxxy_0[j] + wp_x[j] * tg_zz_xxxxxxxy_1[j] + 3.5 * fl1_fxn * tg_zz_xxxxxxy_1[j];

                    tg_xzz_xxxxxxxz_0[j] = pb_x * tg_zz_xxxxxxxz_0[j] + wp_x[j] * tg_zz_xxxxxxxz_1[j] + 3.5 * fl1_fxn * tg_zz_xxxxxxz_1[j];

                    tg_xzz_xxxxxxyy_0[j] = pb_x * tg_zz_xxxxxxyy_0[j] + wp_x[j] * tg_zz_xxxxxxyy_1[j] + 3.0 * fl1_fxn * tg_zz_xxxxxyy_1[j];

                    tg_xzz_xxxxxxyz_0[j] = pb_x * tg_zz_xxxxxxyz_0[j] + wp_x[j] * tg_zz_xxxxxxyz_1[j] + 3.0 * fl1_fxn * tg_zz_xxxxxyz_1[j];

                    tg_xzz_xxxxxxzz_0[j] = pb_x * tg_zz_xxxxxxzz_0[j] + wp_x[j] * tg_zz_xxxxxxzz_1[j] + 3.0 * fl1_fxn * tg_zz_xxxxxzz_1[j];

                    tg_xzz_xxxxxyyy_0[j] = pb_x * tg_zz_xxxxxyyy_0[j] + wp_x[j] * tg_zz_xxxxxyyy_1[j] + 2.5 * fl1_fxn * tg_zz_xxxxyyy_1[j];

                    tg_xzz_xxxxxyyz_0[j] = pb_x * tg_zz_xxxxxyyz_0[j] + wp_x[j] * tg_zz_xxxxxyyz_1[j] + 2.5 * fl1_fxn * tg_zz_xxxxyyz_1[j];

                    tg_xzz_xxxxxyzz_0[j] = pb_x * tg_zz_xxxxxyzz_0[j] + wp_x[j] * tg_zz_xxxxxyzz_1[j] + 2.5 * fl1_fxn * tg_zz_xxxxyzz_1[j];

                    tg_xzz_xxxxxzzz_0[j] = pb_x * tg_zz_xxxxxzzz_0[j] + wp_x[j] * tg_zz_xxxxxzzz_1[j] + 2.5 * fl1_fxn * tg_zz_xxxxzzz_1[j];

                    tg_xzz_xxxxyyyy_0[j] = pb_x * tg_zz_xxxxyyyy_0[j] + wp_x[j] * tg_zz_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_zz_xxxyyyy_1[j];

                    tg_xzz_xxxxyyyz_0[j] = pb_x * tg_zz_xxxxyyyz_0[j] + wp_x[j] * tg_zz_xxxxyyyz_1[j] + 2.0 * fl1_fxn * tg_zz_xxxyyyz_1[j];

                    tg_xzz_xxxxyyzz_0[j] = pb_x * tg_zz_xxxxyyzz_0[j] + wp_x[j] * tg_zz_xxxxyyzz_1[j] + 2.0 * fl1_fxn * tg_zz_xxxyyzz_1[j];

                    tg_xzz_xxxxyzzz_0[j] = pb_x * tg_zz_xxxxyzzz_0[j] + wp_x[j] * tg_zz_xxxxyzzz_1[j] + 2.0 * fl1_fxn * tg_zz_xxxyzzz_1[j];

                    tg_xzz_xxxxzzzz_0[j] = pb_x * tg_zz_xxxxzzzz_0[j] + wp_x[j] * tg_zz_xxxxzzzz_1[j] + 2.0 * fl1_fxn * tg_zz_xxxzzzz_1[j];

                    tg_xzz_xxxyyyyy_0[j] = pb_x * tg_zz_xxxyyyyy_0[j] + wp_x[j] * tg_zz_xxxyyyyy_1[j] + 1.5 * fl1_fxn * tg_zz_xxyyyyy_1[j];

                    tg_xzz_xxxyyyyz_0[j] = pb_x * tg_zz_xxxyyyyz_0[j] + wp_x[j] * tg_zz_xxxyyyyz_1[j] + 1.5 * fl1_fxn * tg_zz_xxyyyyz_1[j];

                    tg_xzz_xxxyyyzz_0[j] = pb_x * tg_zz_xxxyyyzz_0[j] + wp_x[j] * tg_zz_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_zz_xxyyyzz_1[j];

                    tg_xzz_xxxyyzzz_0[j] = pb_x * tg_zz_xxxyyzzz_0[j] + wp_x[j] * tg_zz_xxxyyzzz_1[j] + 1.5 * fl1_fxn * tg_zz_xxyyzzz_1[j];

                    tg_xzz_xxxyzzzz_0[j] = pb_x * tg_zz_xxxyzzzz_0[j] + wp_x[j] * tg_zz_xxxyzzzz_1[j] + 1.5 * fl1_fxn * tg_zz_xxyzzzz_1[j];

                    tg_xzz_xxxzzzzz_0[j] = pb_x * tg_zz_xxxzzzzz_0[j] + wp_x[j] * tg_zz_xxxzzzzz_1[j] + 1.5 * fl1_fxn * tg_zz_xxzzzzz_1[j];

                    tg_xzz_xxyyyyyy_0[j] = pb_x * tg_zz_xxyyyyyy_0[j] + wp_x[j] * tg_zz_xxyyyyyy_1[j] + fl1_fxn * tg_zz_xyyyyyy_1[j];

                    tg_xzz_xxyyyyyz_0[j] = pb_x * tg_zz_xxyyyyyz_0[j] + wp_x[j] * tg_zz_xxyyyyyz_1[j] + fl1_fxn * tg_zz_xyyyyyz_1[j];

                    tg_xzz_xxyyyyzz_0[j] = pb_x * tg_zz_xxyyyyzz_0[j] + wp_x[j] * tg_zz_xxyyyyzz_1[j] + fl1_fxn * tg_zz_xyyyyzz_1[j];

                    tg_xzz_xxyyyzzz_0[j] = pb_x * tg_zz_xxyyyzzz_0[j] + wp_x[j] * tg_zz_xxyyyzzz_1[j] + fl1_fxn * tg_zz_xyyyzzz_1[j];

                    tg_xzz_xxyyzzzz_0[j] = pb_x * tg_zz_xxyyzzzz_0[j] + wp_x[j] * tg_zz_xxyyzzzz_1[j] + fl1_fxn * tg_zz_xyyzzzz_1[j];

                    tg_xzz_xxyzzzzz_0[j] = pb_x * tg_zz_xxyzzzzz_0[j] + wp_x[j] * tg_zz_xxyzzzzz_1[j] + fl1_fxn * tg_zz_xyzzzzz_1[j];

                    tg_xzz_xxzzzzzz_0[j] = pb_x * tg_zz_xxzzzzzz_0[j] + wp_x[j] * tg_zz_xxzzzzzz_1[j] + fl1_fxn * tg_zz_xzzzzzz_1[j];

                    tg_xzz_xyyyyyyy_0[j] = pb_x * tg_zz_xyyyyyyy_0[j] + wp_x[j] * tg_zz_xyyyyyyy_1[j] + 0.5 * fl1_fxn * tg_zz_yyyyyyy_1[j];

                    tg_xzz_xyyyyyyz_0[j] = pb_x * tg_zz_xyyyyyyz_0[j] + wp_x[j] * tg_zz_xyyyyyyz_1[j] + 0.5 * fl1_fxn * tg_zz_yyyyyyz_1[j];

                    tg_xzz_xyyyyyzz_0[j] = pb_x * tg_zz_xyyyyyzz_0[j] + wp_x[j] * tg_zz_xyyyyyzz_1[j] + 0.5 * fl1_fxn * tg_zz_yyyyyzz_1[j];

                    tg_xzz_xyyyyzzz_0[j] = pb_x * tg_zz_xyyyyzzz_0[j] + wp_x[j] * tg_zz_xyyyyzzz_1[j] + 0.5 * fl1_fxn * tg_zz_yyyyzzz_1[j];

                    tg_xzz_xyyyzzzz_0[j] = pb_x * tg_zz_xyyyzzzz_0[j] + wp_x[j] * tg_zz_xyyyzzzz_1[j] + 0.5 * fl1_fxn * tg_zz_yyyzzzz_1[j];

                    tg_xzz_xyyzzzzz_0[j] = pb_x * tg_zz_xyyzzzzz_0[j] + wp_x[j] * tg_zz_xyyzzzzz_1[j] + 0.5 * fl1_fxn * tg_zz_yyzzzzz_1[j];

                    tg_xzz_xyzzzzzz_0[j] = pb_x * tg_zz_xyzzzzzz_0[j] + wp_x[j] * tg_zz_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_zz_yzzzzzz_1[j];

                    tg_xzz_xzzzzzzz_0[j] = pb_x * tg_zz_xzzzzzzz_0[j] + wp_x[j] * tg_zz_xzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_zz_zzzzzzz_1[j];

                    tg_xzz_yyyyyyyy_0[j] = pb_x * tg_zz_yyyyyyyy_0[j] + wp_x[j] * tg_zz_yyyyyyyy_1[j];

                    tg_xzz_yyyyyyyz_0[j] = pb_x * tg_zz_yyyyyyyz_0[j] + wp_x[j] * tg_zz_yyyyyyyz_1[j];

                    tg_xzz_yyyyyyzz_0[j] = pb_x * tg_zz_yyyyyyzz_0[j] + wp_x[j] * tg_zz_yyyyyyzz_1[j];

                    tg_xzz_yyyyyzzz_0[j] = pb_x * tg_zz_yyyyyzzz_0[j] + wp_x[j] * tg_zz_yyyyyzzz_1[j];

                    tg_xzz_yyyyzzzz_0[j] = pb_x * tg_zz_yyyyzzzz_0[j] + wp_x[j] * tg_zz_yyyyzzzz_1[j];

                    tg_xzz_yyyzzzzz_0[j] = pb_x * tg_zz_yyyzzzzz_0[j] + wp_x[j] * tg_zz_yyyzzzzz_1[j];

                    tg_xzz_yyzzzzzz_0[j] = pb_x * tg_zz_yyzzzzzz_0[j] + wp_x[j] * tg_zz_yyzzzzzz_1[j];

                    tg_xzz_yzzzzzzz_0[j] = pb_x * tg_zz_yzzzzzzz_0[j] + wp_x[j] * tg_zz_yzzzzzzz_1[j];

                    tg_xzz_zzzzzzzz_0[j] = pb_x * tg_zz_zzzzzzzz_0[j] + wp_x[j] * tg_zz_zzzzzzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSFSL_270_360(      CMemBlock2D<double>& primBuffer,
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

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {3, -1, -1, -1},
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_3_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_3_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_2_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_2_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_1_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_1_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_2_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {7, -1, -1, -1}, 
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

                auto tg_yy_xxxxxxxx_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 135); 

                auto tg_yy_xxxxxxxy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 136); 

                auto tg_yy_xxxxxxxz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 137); 

                auto tg_yy_xxxxxxyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 138); 

                auto tg_yy_xxxxxxyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 139); 

                auto tg_yy_xxxxxxzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 140); 

                auto tg_yy_xxxxxyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 141); 

                auto tg_yy_xxxxxyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 142); 

                auto tg_yy_xxxxxyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 143); 

                auto tg_yy_xxxxxzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 144); 

                auto tg_yy_xxxxyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 145); 

                auto tg_yy_xxxxyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 146); 

                auto tg_yy_xxxxyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 147); 

                auto tg_yy_xxxxyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 148); 

                auto tg_yy_xxxxzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 149); 

                auto tg_yy_xxxyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 150); 

                auto tg_yy_xxxyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 151); 

                auto tg_yy_xxxyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 152); 

                auto tg_yy_xxxyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 153); 

                auto tg_yy_xxxyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 154); 

                auto tg_yy_xxxzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 155); 

                auto tg_yy_xxyyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 156); 

                auto tg_yy_xxyyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 157); 

                auto tg_yy_xxyyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 158); 

                auto tg_yy_xxyyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 159); 

                auto tg_yy_xxyyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 160); 

                auto tg_yy_xxyzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 161); 

                auto tg_yy_xxzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 162); 

                auto tg_yy_xyyyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 163); 

                auto tg_yy_xyyyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 164); 

                auto tg_yy_xyyyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 165); 

                auto tg_yy_xyyyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 166); 

                auto tg_yy_xyyyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 167); 

                auto tg_yy_xyyzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 168); 

                auto tg_yy_xyzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 169); 

                auto tg_yy_xzzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 170); 

                auto tg_yy_yyyyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 171); 

                auto tg_yy_yyyyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 172); 

                auto tg_yy_yyyyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 173); 

                auto tg_yy_yyyyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 174); 

                auto tg_yy_yyyyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 175); 

                auto tg_yy_yyyzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 176); 

                auto tg_yy_yyzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 177); 

                auto tg_yy_yzzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 178); 

                auto tg_yy_zzzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 179); 

                auto tg_yz_xxxxxxxx_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 180); 

                auto tg_yz_xxxxxxxy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 181); 

                auto tg_yz_xxxxxxxz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 182); 

                auto tg_yz_xxxxxxyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 183); 

                auto tg_yz_xxxxxxyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 184); 

                auto tg_yz_xxxxxxzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 185); 

                auto tg_yz_xxxxxyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 186); 

                auto tg_yz_xxxxxyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 187); 

                auto tg_yz_xxxxxyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 188); 

                auto tg_yz_xxxxxzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 189); 

                auto tg_yz_xxxxyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 190); 

                auto tg_yz_xxxxyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 191); 

                auto tg_yz_xxxxyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 192); 

                auto tg_yz_xxxxyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 193); 

                auto tg_yz_xxxxzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 194); 

                auto tg_yz_xxxyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 195); 

                auto tg_yz_xxxyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 196); 

                auto tg_yz_xxxyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 197); 

                auto tg_yz_xxxyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 198); 

                auto tg_yz_xxxyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 199); 

                auto tg_yz_xxxzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 200); 

                auto tg_yz_xxyyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 201); 

                auto tg_yz_xxyyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 202); 

                auto tg_yz_xxyyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 203); 

                auto tg_yz_xxyyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 204); 

                auto tg_yz_xxyyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 205); 

                auto tg_yz_xxyzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 206); 

                auto tg_yz_xxzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 207); 

                auto tg_yz_xyyyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 208); 

                auto tg_yz_xyyyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 209); 

                auto tg_yz_xyyyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 210); 

                auto tg_yz_xyyyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 211); 

                auto tg_yz_xyyyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 212); 

                auto tg_yz_xyyzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 213); 

                auto tg_yz_xyzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 214); 

                auto tg_yz_xzzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 215); 

                auto tg_yz_yyyyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 216); 

                auto tg_yz_yyyyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 217); 

                auto tg_yz_yyyyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 218); 

                auto tg_yz_yyyyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 219); 

                auto tg_yz_yyyyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 220); 

                auto tg_yz_yyyzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 221); 

                auto tg_yz_yyzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 222); 

                auto tg_yz_yzzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 223); 

                auto tg_yz_zzzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 224); 

                auto tg_yy_xxxxxxxx_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 135); 

                auto tg_yy_xxxxxxxy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 136); 

                auto tg_yy_xxxxxxxz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 137); 

                auto tg_yy_xxxxxxyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 138); 

                auto tg_yy_xxxxxxyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 139); 

                auto tg_yy_xxxxxxzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 140); 

                auto tg_yy_xxxxxyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 141); 

                auto tg_yy_xxxxxyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 142); 

                auto tg_yy_xxxxxyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 143); 

                auto tg_yy_xxxxxzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 144); 

                auto tg_yy_xxxxyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 145); 

                auto tg_yy_xxxxyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 146); 

                auto tg_yy_xxxxyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 147); 

                auto tg_yy_xxxxyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 148); 

                auto tg_yy_xxxxzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 149); 

                auto tg_yy_xxxyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 150); 

                auto tg_yy_xxxyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 151); 

                auto tg_yy_xxxyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 152); 

                auto tg_yy_xxxyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 153); 

                auto tg_yy_xxxyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 154); 

                auto tg_yy_xxxzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 155); 

                auto tg_yy_xxyyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 156); 

                auto tg_yy_xxyyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 157); 

                auto tg_yy_xxyyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 158); 

                auto tg_yy_xxyyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 159); 

                auto tg_yy_xxyyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 160); 

                auto tg_yy_xxyzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 161); 

                auto tg_yy_xxzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 162); 

                auto tg_yy_xyyyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 163); 

                auto tg_yy_xyyyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 164); 

                auto tg_yy_xyyyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 165); 

                auto tg_yy_xyyyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 166); 

                auto tg_yy_xyyyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 167); 

                auto tg_yy_xyyzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 168); 

                auto tg_yy_xyzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 169); 

                auto tg_yy_xzzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 170); 

                auto tg_yy_yyyyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 171); 

                auto tg_yy_yyyyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 172); 

                auto tg_yy_yyyyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 173); 

                auto tg_yy_yyyyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 174); 

                auto tg_yy_yyyyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 175); 

                auto tg_yy_yyyzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 176); 

                auto tg_yy_yyzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 177); 

                auto tg_yy_yzzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 178); 

                auto tg_yy_zzzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 179); 

                auto tg_yz_xxxxxxxx_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 180); 

                auto tg_yz_xxxxxxxy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 181); 

                auto tg_yz_xxxxxxxz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 182); 

                auto tg_yz_xxxxxxyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 183); 

                auto tg_yz_xxxxxxyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 184); 

                auto tg_yz_xxxxxxzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 185); 

                auto tg_yz_xxxxxyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 186); 

                auto tg_yz_xxxxxyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 187); 

                auto tg_yz_xxxxxyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 188); 

                auto tg_yz_xxxxxzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 189); 

                auto tg_yz_xxxxyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 190); 

                auto tg_yz_xxxxyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 191); 

                auto tg_yz_xxxxyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 192); 

                auto tg_yz_xxxxyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 193); 

                auto tg_yz_xxxxzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 194); 

                auto tg_yz_xxxyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 195); 

                auto tg_yz_xxxyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 196); 

                auto tg_yz_xxxyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 197); 

                auto tg_yz_xxxyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 198); 

                auto tg_yz_xxxyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 199); 

                auto tg_yz_xxxzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 200); 

                auto tg_yz_xxyyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 201); 

                auto tg_yz_xxyyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 202); 

                auto tg_yz_xxyyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 203); 

                auto tg_yz_xxyyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 204); 

                auto tg_yz_xxyyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 205); 

                auto tg_yz_xxyzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 206); 

                auto tg_yz_xxzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 207); 

                auto tg_yz_xyyyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 208); 

                auto tg_yz_xyyyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 209); 

                auto tg_yz_xyyyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 210); 

                auto tg_yz_xyyyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 211); 

                auto tg_yz_xyyyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 212); 

                auto tg_yz_xyyzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 213); 

                auto tg_yz_xyzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 214); 

                auto tg_yz_xzzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 215); 

                auto tg_yz_yyyyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 216); 

                auto tg_yz_yyyyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 217); 

                auto tg_yz_yyyyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 218); 

                auto tg_yz_yyyyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 219); 

                auto tg_yz_yyyyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 220); 

                auto tg_yz_yyyzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 221); 

                auto tg_yz_yyzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 222); 

                auto tg_yz_yzzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 223); 

                auto tg_yz_zzzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 224); 

                auto tg_y_xxxxxxxx_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 45); 

                auto tg_y_xxxxxxxy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 46); 

                auto tg_y_xxxxxxxz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 47); 

                auto tg_y_xxxxxxyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 48); 

                auto tg_y_xxxxxxyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 49); 

                auto tg_y_xxxxxxzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 50); 

                auto tg_y_xxxxxyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 51); 

                auto tg_y_xxxxxyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 52); 

                auto tg_y_xxxxxyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 53); 

                auto tg_y_xxxxxzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 54); 

                auto tg_y_xxxxyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 55); 

                auto tg_y_xxxxyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 56); 

                auto tg_y_xxxxyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 57); 

                auto tg_y_xxxxyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 58); 

                auto tg_y_xxxxzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 59); 

                auto tg_y_xxxyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 60); 

                auto tg_y_xxxyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 61); 

                auto tg_y_xxxyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 62); 

                auto tg_y_xxxyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 63); 

                auto tg_y_xxxyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 64); 

                auto tg_y_xxxzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 65); 

                auto tg_y_xxyyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 66); 

                auto tg_y_xxyyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 67); 

                auto tg_y_xxyyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 68); 

                auto tg_y_xxyyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 69); 

                auto tg_y_xxyyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 70); 

                auto tg_y_xxyzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 71); 

                auto tg_y_xxzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 72); 

                auto tg_y_xyyyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 73); 

                auto tg_y_xyyyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 74); 

                auto tg_y_xyyyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 75); 

                auto tg_y_xyyyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 76); 

                auto tg_y_xyyyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 77); 

                auto tg_y_xyyzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 78); 

                auto tg_y_xyzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 79); 

                auto tg_y_xzzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 80); 

                auto tg_y_yyyyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 81); 

                auto tg_y_yyyyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 82); 

                auto tg_y_yyyyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 83); 

                auto tg_y_yyyyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 84); 

                auto tg_y_yyyyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 85); 

                auto tg_y_yyyzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 86); 

                auto tg_y_yyzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 87); 

                auto tg_y_yzzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 88); 

                auto tg_y_zzzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 89); 

                auto tg_z_xxxxxxxx_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 90); 

                auto tg_z_xxxxxxxy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 91); 

                auto tg_z_xxxxxxxz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 92); 

                auto tg_z_xxxxxxyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 93); 

                auto tg_z_xxxxxxyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 94); 

                auto tg_z_xxxxxxzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 95); 

                auto tg_z_xxxxxyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 96); 

                auto tg_z_xxxxxyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 97); 

                auto tg_z_xxxxxyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 98); 

                auto tg_z_xxxxxzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 99); 

                auto tg_z_xxxxyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 100); 

                auto tg_z_xxxxyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 101); 

                auto tg_z_xxxxyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 102); 

                auto tg_z_xxxxyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 103); 

                auto tg_z_xxxxzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 104); 

                auto tg_z_xxxyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 105); 

                auto tg_z_xxxyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 106); 

                auto tg_z_xxxyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 107); 

                auto tg_z_xxxyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 108); 

                auto tg_z_xxxyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 109); 

                auto tg_z_xxxzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 110); 

                auto tg_z_xxyyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 111); 

                auto tg_z_xxyyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 112); 

                auto tg_z_xxyyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 113); 

                auto tg_z_xxyyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 114); 

                auto tg_z_xxyyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 115); 

                auto tg_z_xxyzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 116); 

                auto tg_z_xxzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 117); 

                auto tg_z_xyyyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 118); 

                auto tg_z_xyyyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 119); 

                auto tg_z_xyyyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 120); 

                auto tg_z_xyyyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 121); 

                auto tg_z_xyyyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 122); 

                auto tg_z_xyyzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 123); 

                auto tg_z_xyzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 124); 

                auto tg_z_xzzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 125); 

                auto tg_z_yyyyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 126); 

                auto tg_z_yyyyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 127); 

                auto tg_z_yyyyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 128); 

                auto tg_z_yyyyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 129); 

                auto tg_z_yyyyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 130); 

                auto tg_z_yyyzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 131); 

                auto tg_z_yyzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 132); 

                auto tg_z_yzzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 133); 

                auto tg_z_zzzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 134); 

                auto tg_y_xxxxxxxx_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 45); 

                auto tg_y_xxxxxxxy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 46); 

                auto tg_y_xxxxxxxz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 47); 

                auto tg_y_xxxxxxyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 48); 

                auto tg_y_xxxxxxyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 49); 

                auto tg_y_xxxxxxzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 50); 

                auto tg_y_xxxxxyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 51); 

                auto tg_y_xxxxxyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 52); 

                auto tg_y_xxxxxyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 53); 

                auto tg_y_xxxxxzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 54); 

                auto tg_y_xxxxyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 55); 

                auto tg_y_xxxxyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 56); 

                auto tg_y_xxxxyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 57); 

                auto tg_y_xxxxyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 58); 

                auto tg_y_xxxxzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 59); 

                auto tg_y_xxxyyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 60); 

                auto tg_y_xxxyyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 61); 

                auto tg_y_xxxyyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 62); 

                auto tg_y_xxxyyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 63); 

                auto tg_y_xxxyzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 64); 

                auto tg_y_xxxzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 65); 

                auto tg_y_xxyyyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 66); 

                auto tg_y_xxyyyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 67); 

                auto tg_y_xxyyyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 68); 

                auto tg_y_xxyyyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 69); 

                auto tg_y_xxyyzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 70); 

                auto tg_y_xxyzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 71); 

                auto tg_y_xxzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 72); 

                auto tg_y_xyyyyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 73); 

                auto tg_y_xyyyyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 74); 

                auto tg_y_xyyyyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 75); 

                auto tg_y_xyyyyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 76); 

                auto tg_y_xyyyzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 77); 

                auto tg_y_xyyzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 78); 

                auto tg_y_xyzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 79); 

                auto tg_y_xzzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 80); 

                auto tg_y_yyyyyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 81); 

                auto tg_y_yyyyyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 82); 

                auto tg_y_yyyyyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 83); 

                auto tg_y_yyyyyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 84); 

                auto tg_y_yyyyzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 85); 

                auto tg_y_yyyzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 86); 

                auto tg_y_yyzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 87); 

                auto tg_y_yzzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 88); 

                auto tg_y_zzzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 89); 

                auto tg_z_xxxxxxxx_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 90); 

                auto tg_z_xxxxxxxy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 91); 

                auto tg_z_xxxxxxxz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 92); 

                auto tg_z_xxxxxxyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 93); 

                auto tg_z_xxxxxxyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 94); 

                auto tg_z_xxxxxxzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 95); 

                auto tg_z_xxxxxyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 96); 

                auto tg_z_xxxxxyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 97); 

                auto tg_z_xxxxxyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 98); 

                auto tg_z_xxxxxzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 99); 

                auto tg_z_xxxxyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 100); 

                auto tg_z_xxxxyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 101); 

                auto tg_z_xxxxyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 102); 

                auto tg_z_xxxxyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 103); 

                auto tg_z_xxxxzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 104); 

                auto tg_z_xxxyyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 105); 

                auto tg_z_xxxyyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 106); 

                auto tg_z_xxxyyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 107); 

                auto tg_z_xxxyyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 108); 

                auto tg_z_xxxyzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 109); 

                auto tg_z_xxxzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 110); 

                auto tg_z_xxyyyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 111); 

                auto tg_z_xxyyyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 112); 

                auto tg_z_xxyyyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 113); 

                auto tg_z_xxyyyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 114); 

                auto tg_z_xxyyzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 115); 

                auto tg_z_xxyzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 116); 

                auto tg_z_xxzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 117); 

                auto tg_z_xyyyyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 118); 

                auto tg_z_xyyyyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 119); 

                auto tg_z_xyyyyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 120); 

                auto tg_z_xyyyyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 121); 

                auto tg_z_xyyyzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 122); 

                auto tg_z_xyyzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 123); 

                auto tg_z_xyzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 124); 

                auto tg_z_xzzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 125); 

                auto tg_z_yyyyyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 126); 

                auto tg_z_yyyyyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 127); 

                auto tg_z_yyyyyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 128); 

                auto tg_z_yyyyyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 129); 

                auto tg_z_yyyyzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 130); 

                auto tg_z_yyyzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 131); 

                auto tg_z_yyzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 132); 

                auto tg_z_yzzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 133); 

                auto tg_z_zzzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 134); 

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

                // set up pointers to integrals

                auto tg_yyy_xxxxxxxx_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 270); 

                auto tg_yyy_xxxxxxxy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 271); 

                auto tg_yyy_xxxxxxxz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 272); 

                auto tg_yyy_xxxxxxyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 273); 

                auto tg_yyy_xxxxxxyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 274); 

                auto tg_yyy_xxxxxxzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 275); 

                auto tg_yyy_xxxxxyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 276); 

                auto tg_yyy_xxxxxyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 277); 

                auto tg_yyy_xxxxxyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 278); 

                auto tg_yyy_xxxxxzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 279); 

                auto tg_yyy_xxxxyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 280); 

                auto tg_yyy_xxxxyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 281); 

                auto tg_yyy_xxxxyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 282); 

                auto tg_yyy_xxxxyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 283); 

                auto tg_yyy_xxxxzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 284); 

                auto tg_yyy_xxxyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 285); 

                auto tg_yyy_xxxyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 286); 

                auto tg_yyy_xxxyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 287); 

                auto tg_yyy_xxxyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 288); 

                auto tg_yyy_xxxyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 289); 

                auto tg_yyy_xxxzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 290); 

                auto tg_yyy_xxyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 291); 

                auto tg_yyy_xxyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 292); 

                auto tg_yyy_xxyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 293); 

                auto tg_yyy_xxyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 294); 

                auto tg_yyy_xxyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 295); 

                auto tg_yyy_xxyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 296); 

                auto tg_yyy_xxzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 297); 

                auto tg_yyy_xyyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 298); 

                auto tg_yyy_xyyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 299); 

                auto tg_yyy_xyyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 300); 

                auto tg_yyy_xyyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 301); 

                auto tg_yyy_xyyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 302); 

                auto tg_yyy_xyyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 303); 

                auto tg_yyy_xyzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 304); 

                auto tg_yyy_xzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 305); 

                auto tg_yyy_yyyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 306); 

                auto tg_yyy_yyyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 307); 

                auto tg_yyy_yyyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 308); 

                auto tg_yyy_yyyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 309); 

                auto tg_yyy_yyyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 310); 

                auto tg_yyy_yyyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 311); 

                auto tg_yyy_yyzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 312); 

                auto tg_yyy_yzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 313); 

                auto tg_yyy_zzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 314); 

                auto tg_yyz_xxxxxxxx_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 315); 

                auto tg_yyz_xxxxxxxy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 316); 

                auto tg_yyz_xxxxxxxz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 317); 

                auto tg_yyz_xxxxxxyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 318); 

                auto tg_yyz_xxxxxxyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 319); 

                auto tg_yyz_xxxxxxzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 320); 

                auto tg_yyz_xxxxxyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 321); 

                auto tg_yyz_xxxxxyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 322); 

                auto tg_yyz_xxxxxyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 323); 

                auto tg_yyz_xxxxxzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 324); 

                auto tg_yyz_xxxxyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 325); 

                auto tg_yyz_xxxxyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 326); 

                auto tg_yyz_xxxxyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 327); 

                auto tg_yyz_xxxxyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 328); 

                auto tg_yyz_xxxxzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 329); 

                auto tg_yyz_xxxyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 330); 

                auto tg_yyz_xxxyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 331); 

                auto tg_yyz_xxxyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 332); 

                auto tg_yyz_xxxyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 333); 

                auto tg_yyz_xxxyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 334); 

                auto tg_yyz_xxxzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 335); 

                auto tg_yyz_xxyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 336); 

                auto tg_yyz_xxyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 337); 

                auto tg_yyz_xxyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 338); 

                auto tg_yyz_xxyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 339); 

                auto tg_yyz_xxyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 340); 

                auto tg_yyz_xxyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 341); 

                auto tg_yyz_xxzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 342); 

                auto tg_yyz_xyyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 343); 

                auto tg_yyz_xyyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 344); 

                auto tg_yyz_xyyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 345); 

                auto tg_yyz_xyyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 346); 

                auto tg_yyz_xyyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 347); 

                auto tg_yyz_xyyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 348); 

                auto tg_yyz_xyzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 349); 

                auto tg_yyz_xzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 350); 

                auto tg_yyz_yyyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 351); 

                auto tg_yyz_yyyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 352); 

                auto tg_yyz_yyyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 353); 

                auto tg_yyz_yyyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 354); 

                auto tg_yyz_yyyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 355); 

                auto tg_yyz_yyyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 356); 

                auto tg_yyz_yyzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 357); 

                auto tg_yyz_yzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 358); 

                auto tg_yyz_zzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 359); 

                // Batch of Integrals (270,360)

                #pragma omp simd aligned(fxn, fza, tg_y_xxxxxxxx_0, tg_y_xxxxxxxx_1, tg_y_xxxxxxxy_0, \
                                         tg_y_xxxxxxxy_1, tg_y_xxxxxxxz_0, tg_y_xxxxxxxz_1, tg_y_xxxxxxyy_0, tg_y_xxxxxxyy_1, \
                                         tg_y_xxxxxxyz_0, tg_y_xxxxxxyz_1, tg_y_xxxxxxzz_0, tg_y_xxxxxxzz_1, tg_y_xxxxxyyy_0, \
                                         tg_y_xxxxxyyy_1, tg_y_xxxxxyyz_0, tg_y_xxxxxyyz_1, tg_y_xxxxxyzz_0, tg_y_xxxxxyzz_1, \
                                         tg_y_xxxxxzzz_0, tg_y_xxxxxzzz_1, tg_y_xxxxyyyy_0, tg_y_xxxxyyyy_1, tg_y_xxxxyyyz_0, \
                                         tg_y_xxxxyyyz_1, tg_y_xxxxyyzz_0, tg_y_xxxxyyzz_1, tg_y_xxxxyzzz_0, tg_y_xxxxyzzz_1, \
                                         tg_y_xxxxzzzz_0, tg_y_xxxxzzzz_1, tg_y_xxxyyyyy_0, tg_y_xxxyyyyy_1, tg_y_xxxyyyyz_0, \
                                         tg_y_xxxyyyyz_1, tg_y_xxxyyyzz_0, tg_y_xxxyyyzz_1, tg_y_xxxyyzzz_0, tg_y_xxxyyzzz_1, \
                                         tg_y_xxxyzzzz_0, tg_y_xxxyzzzz_1, tg_y_xxxzzzzz_0, tg_y_xxxzzzzz_1, tg_y_xxyyyyyy_0, \
                                         tg_y_xxyyyyyy_1, tg_y_xxyyyyyz_0, tg_y_xxyyyyyz_1, tg_y_xxyyyyzz_0, tg_y_xxyyyyzz_1, \
                                         tg_y_xxyyyzzz_0, tg_y_xxyyyzzz_1, tg_y_xxyyzzzz_0, tg_y_xxyyzzzz_1, tg_y_xxyzzzzz_0, \
                                         tg_y_xxyzzzzz_1, tg_y_xxzzzzzz_0, tg_y_xxzzzzzz_1, tg_y_xyyyyyyy_0, tg_y_xyyyyyyy_1, \
                                         tg_y_xyyyyyyz_0, tg_y_xyyyyyyz_1, tg_y_xyyyyyzz_0, tg_y_xyyyyyzz_1, tg_y_xyyyyzzz_0, \
                                         tg_y_xyyyyzzz_1, tg_y_xyyyzzzz_0, tg_y_xyyyzzzz_1, tg_y_xyyzzzzz_0, tg_y_xyyzzzzz_1, \
                                         tg_y_xyzzzzzz_0, tg_y_xyzzzzzz_1, tg_y_xzzzzzzz_0, tg_y_xzzzzzzz_1, tg_y_yyyyyyyy_0, \
                                         tg_y_yyyyyyyy_1, tg_y_yyyyyyyz_0, tg_y_yyyyyyyz_1, tg_y_yyyyyyzz_0, tg_y_yyyyyyzz_1, \
                                         tg_y_yyyyyzzz_0, tg_y_yyyyyzzz_1, tg_y_yyyyzzzz_0, tg_y_yyyyzzzz_1, tg_y_yyyzzzzz_0, \
                                         tg_y_yyyzzzzz_1, tg_y_yyzzzzzz_0, tg_y_yyzzzzzz_1, tg_y_yzzzzzzz_0, tg_y_yzzzzzzz_1, \
                                         tg_y_zzzzzzzz_0, tg_y_zzzzzzzz_1, tg_yy_xxxxxxx_1, tg_yy_xxxxxxxx_0, \
                                         tg_yy_xxxxxxxx_1, tg_yy_xxxxxxxy_0, tg_yy_xxxxxxxy_1, tg_yy_xxxxxxxz_0, \
                                         tg_yy_xxxxxxxz_1, tg_yy_xxxxxxy_1, tg_yy_xxxxxxyy_0, tg_yy_xxxxxxyy_1, \
                                         tg_yy_xxxxxxyz_0, tg_yy_xxxxxxyz_1, tg_yy_xxxxxxz_1, tg_yy_xxxxxxzz_0, \
                                         tg_yy_xxxxxxzz_1, tg_yy_xxxxxyy_1, tg_yy_xxxxxyyy_0, tg_yy_xxxxxyyy_1, \
                                         tg_yy_xxxxxyyz_0, tg_yy_xxxxxyyz_1, tg_yy_xxxxxyz_1, tg_yy_xxxxxyzz_0, \
                                         tg_yy_xxxxxyzz_1, tg_yy_xxxxxzz_1, tg_yy_xxxxxzzz_0, tg_yy_xxxxxzzz_1, \
                                         tg_yy_xxxxyyy_1, tg_yy_xxxxyyyy_0, tg_yy_xxxxyyyy_1, tg_yy_xxxxyyyz_0, \
                                         tg_yy_xxxxyyyz_1, tg_yy_xxxxyyz_1, tg_yy_xxxxyyzz_0, tg_yy_xxxxyyzz_1, \
                                         tg_yy_xxxxyzz_1, tg_yy_xxxxyzzz_0, tg_yy_xxxxyzzz_1, tg_yy_xxxxzzz_1, \
                                         tg_yy_xxxxzzzz_0, tg_yy_xxxxzzzz_1, tg_yy_xxxyyyy_1, tg_yy_xxxyyyyy_0, \
                                         tg_yy_xxxyyyyy_1, tg_yy_xxxyyyyz_0, tg_yy_xxxyyyyz_1, tg_yy_xxxyyyz_1, \
                                         tg_yy_xxxyyyzz_0, tg_yy_xxxyyyzz_1, tg_yy_xxxyyzz_1, tg_yy_xxxyyzzz_0, \
                                         tg_yy_xxxyyzzz_1, tg_yy_xxxyzzz_1, tg_yy_xxxyzzzz_0, tg_yy_xxxyzzzz_1, \
                                         tg_yy_xxxzzzz_1, tg_yy_xxxzzzzz_0, tg_yy_xxxzzzzz_1, tg_yy_xxyyyyy_1, \
                                         tg_yy_xxyyyyyy_0, tg_yy_xxyyyyyy_1, tg_yy_xxyyyyyz_0, tg_yy_xxyyyyyz_1, \
                                         tg_yy_xxyyyyz_1, tg_yy_xxyyyyzz_0, tg_yy_xxyyyyzz_1, tg_yy_xxyyyzz_1, \
                                         tg_yy_xxyyyzzz_0, tg_yy_xxyyyzzz_1, tg_yy_xxyyzzz_1, tg_yy_xxyyzzzz_0, \
                                         tg_yy_xxyyzzzz_1, tg_yy_xxyzzzz_1, tg_yy_xxyzzzzz_0, tg_yy_xxyzzzzz_1, \
                                         tg_yy_xxzzzzz_1, tg_yy_xxzzzzzz_0, tg_yy_xxzzzzzz_1, tg_yy_xyyyyyy_1, \
                                         tg_yy_xyyyyyyy_0, tg_yy_xyyyyyyy_1, tg_yy_xyyyyyyz_0, tg_yy_xyyyyyyz_1, \
                                         tg_yy_xyyyyyz_1, tg_yy_xyyyyyzz_0, tg_yy_xyyyyyzz_1, tg_yy_xyyyyzz_1, \
                                         tg_yy_xyyyyzzz_0, tg_yy_xyyyyzzz_1, tg_yy_xyyyzzz_1, tg_yy_xyyyzzzz_0, \
                                         tg_yy_xyyyzzzz_1, tg_yy_xyyzzzz_1, tg_yy_xyyzzzzz_0, tg_yy_xyyzzzzz_1, \
                                         tg_yy_xyzzzzz_1, tg_yy_xyzzzzzz_0, tg_yy_xyzzzzzz_1, tg_yy_xzzzzzz_1, \
                                         tg_yy_xzzzzzzz_0, tg_yy_xzzzzzzz_1, tg_yy_yyyyyyy_1, tg_yy_yyyyyyyy_0, \
                                         tg_yy_yyyyyyyy_1, tg_yy_yyyyyyyz_0, tg_yy_yyyyyyyz_1, tg_yy_yyyyyyz_1, \
                                         tg_yy_yyyyyyzz_0, tg_yy_yyyyyyzz_1, tg_yy_yyyyyzz_1, tg_yy_yyyyyzzz_0, \
                                         tg_yy_yyyyyzzz_1, tg_yy_yyyyzzz_1, tg_yy_yyyyzzzz_0, tg_yy_yyyyzzzz_1, \
                                         tg_yy_yyyzzzz_1, tg_yy_yyyzzzzz_0, tg_yy_yyyzzzzz_1, tg_yy_yyzzzzz_1, \
                                         tg_yy_yyzzzzzz_0, tg_yy_yyzzzzzz_1, tg_yy_yzzzzzz_1, tg_yy_yzzzzzzz_0, \
                                         tg_yy_yzzzzzzz_1, tg_yy_zzzzzzz_1, tg_yy_zzzzzzzz_0, tg_yy_zzzzzzzz_1, \
                                         tg_yyy_xxxxxxxx_0, tg_yyy_xxxxxxxy_0, tg_yyy_xxxxxxxz_0, tg_yyy_xxxxxxyy_0, \
                                         tg_yyy_xxxxxxyz_0, tg_yyy_xxxxxxzz_0, tg_yyy_xxxxxyyy_0, tg_yyy_xxxxxyyz_0, \
                                         tg_yyy_xxxxxyzz_0, tg_yyy_xxxxxzzz_0, tg_yyy_xxxxyyyy_0, tg_yyy_xxxxyyyz_0, \
                                         tg_yyy_xxxxyyzz_0, tg_yyy_xxxxyzzz_0, tg_yyy_xxxxzzzz_0, tg_yyy_xxxyyyyy_0, \
                                         tg_yyy_xxxyyyyz_0, tg_yyy_xxxyyyzz_0, tg_yyy_xxxyyzzz_0, tg_yyy_xxxyzzzz_0, \
                                         tg_yyy_xxxzzzzz_0, tg_yyy_xxyyyyyy_0, tg_yyy_xxyyyyyz_0, tg_yyy_xxyyyyzz_0, \
                                         tg_yyy_xxyyyzzz_0, tg_yyy_xxyyzzzz_0, tg_yyy_xxyzzzzz_0, tg_yyy_xxzzzzzz_0, \
                                         tg_yyy_xyyyyyyy_0, tg_yyy_xyyyyyyz_0, tg_yyy_xyyyyyzz_0, tg_yyy_xyyyyzzz_0, \
                                         tg_yyy_xyyyzzzz_0, tg_yyy_xyyzzzzz_0, tg_yyy_xyzzzzzz_0, tg_yyy_xzzzzzzz_0, \
                                         tg_yyy_yyyyyyyy_0, tg_yyy_yyyyyyyz_0, tg_yyy_yyyyyyzz_0, tg_yyy_yyyyyzzz_0, \
                                         tg_yyy_yyyyzzzz_0, tg_yyy_yyyzzzzz_0, tg_yyy_yyzzzzzz_0, tg_yyy_yzzzzzzz_0, \
                                         tg_yyy_zzzzzzzz_0, tg_yyz_xxxxxxxx_0, tg_yyz_xxxxxxxy_0, tg_yyz_xxxxxxxz_0, \
                                         tg_yyz_xxxxxxyy_0, tg_yyz_xxxxxxyz_0, tg_yyz_xxxxxxzz_0, tg_yyz_xxxxxyyy_0, \
                                         tg_yyz_xxxxxyyz_0, tg_yyz_xxxxxyzz_0, tg_yyz_xxxxxzzz_0, tg_yyz_xxxxyyyy_0, \
                                         tg_yyz_xxxxyyyz_0, tg_yyz_xxxxyyzz_0, tg_yyz_xxxxyzzz_0, tg_yyz_xxxxzzzz_0, \
                                         tg_yyz_xxxyyyyy_0, tg_yyz_xxxyyyyz_0, tg_yyz_xxxyyyzz_0, tg_yyz_xxxyyzzz_0, \
                                         tg_yyz_xxxyzzzz_0, tg_yyz_xxxzzzzz_0, tg_yyz_xxyyyyyy_0, tg_yyz_xxyyyyyz_0, \
                                         tg_yyz_xxyyyyzz_0, tg_yyz_xxyyyzzz_0, tg_yyz_xxyyzzzz_0, tg_yyz_xxyzzzzz_0, \
                                         tg_yyz_xxzzzzzz_0, tg_yyz_xyyyyyyy_0, tg_yyz_xyyyyyyz_0, tg_yyz_xyyyyyzz_0, \
                                         tg_yyz_xyyyyzzz_0, tg_yyz_xyyyzzzz_0, tg_yyz_xyyzzzzz_0, tg_yyz_xyzzzzzz_0, \
                                         tg_yyz_xzzzzzzz_0, tg_yyz_yyyyyyyy_0, tg_yyz_yyyyyyyz_0, tg_yyz_yyyyyyzz_0, \
                                         tg_yyz_yyyyyzzz_0, tg_yyz_yyyyzzzz_0, tg_yyz_yyyzzzzz_0, tg_yyz_yyzzzzzz_0, \
                                         tg_yyz_yzzzzzzz_0, tg_yyz_zzzzzzzz_0, tg_yz_xxxxxxx_1, tg_yz_xxxxxxxx_0, \
                                         tg_yz_xxxxxxxx_1, tg_yz_xxxxxxxy_0, tg_yz_xxxxxxxy_1, tg_yz_xxxxxxxz_0, \
                                         tg_yz_xxxxxxxz_1, tg_yz_xxxxxxy_1, tg_yz_xxxxxxyy_0, tg_yz_xxxxxxyy_1, \
                                         tg_yz_xxxxxxyz_0, tg_yz_xxxxxxyz_1, tg_yz_xxxxxxz_1, tg_yz_xxxxxxzz_0, \
                                         tg_yz_xxxxxxzz_1, tg_yz_xxxxxyy_1, tg_yz_xxxxxyyy_0, tg_yz_xxxxxyyy_1, \
                                         tg_yz_xxxxxyyz_0, tg_yz_xxxxxyyz_1, tg_yz_xxxxxyz_1, tg_yz_xxxxxyzz_0, \
                                         tg_yz_xxxxxyzz_1, tg_yz_xxxxxzz_1, tg_yz_xxxxxzzz_0, tg_yz_xxxxxzzz_1, \
                                         tg_yz_xxxxyyy_1, tg_yz_xxxxyyyy_0, tg_yz_xxxxyyyy_1, tg_yz_xxxxyyyz_0, \
                                         tg_yz_xxxxyyyz_1, tg_yz_xxxxyyz_1, tg_yz_xxxxyyzz_0, tg_yz_xxxxyyzz_1, \
                                         tg_yz_xxxxyzz_1, tg_yz_xxxxyzzz_0, tg_yz_xxxxyzzz_1, tg_yz_xxxxzzz_1, \
                                         tg_yz_xxxxzzzz_0, tg_yz_xxxxzzzz_1, tg_yz_xxxyyyy_1, tg_yz_xxxyyyyy_0, \
                                         tg_yz_xxxyyyyy_1, tg_yz_xxxyyyyz_0, tg_yz_xxxyyyyz_1, tg_yz_xxxyyyz_1, \
                                         tg_yz_xxxyyyzz_0, tg_yz_xxxyyyzz_1, tg_yz_xxxyyzz_1, tg_yz_xxxyyzzz_0, \
                                         tg_yz_xxxyyzzz_1, tg_yz_xxxyzzz_1, tg_yz_xxxyzzzz_0, tg_yz_xxxyzzzz_1, \
                                         tg_yz_xxxzzzz_1, tg_yz_xxxzzzzz_0, tg_yz_xxxzzzzz_1, tg_yz_xxyyyyy_1, \
                                         tg_yz_xxyyyyyy_0, tg_yz_xxyyyyyy_1, tg_yz_xxyyyyyz_0, tg_yz_xxyyyyyz_1, \
                                         tg_yz_xxyyyyz_1, tg_yz_xxyyyyzz_0, tg_yz_xxyyyyzz_1, tg_yz_xxyyyzz_1, \
                                         tg_yz_xxyyyzzz_0, tg_yz_xxyyyzzz_1, tg_yz_xxyyzzz_1, tg_yz_xxyyzzzz_0, \
                                         tg_yz_xxyyzzzz_1, tg_yz_xxyzzzz_1, tg_yz_xxyzzzzz_0, tg_yz_xxyzzzzz_1, \
                                         tg_yz_xxzzzzz_1, tg_yz_xxzzzzzz_0, tg_yz_xxzzzzzz_1, tg_yz_xyyyyyy_1, \
                                         tg_yz_xyyyyyyy_0, tg_yz_xyyyyyyy_1, tg_yz_xyyyyyyz_0, tg_yz_xyyyyyyz_1, \
                                         tg_yz_xyyyyyz_1, tg_yz_xyyyyyzz_0, tg_yz_xyyyyyzz_1, tg_yz_xyyyyzz_1, \
                                         tg_yz_xyyyyzzz_0, tg_yz_xyyyyzzz_1, tg_yz_xyyyzzz_1, tg_yz_xyyyzzzz_0, \
                                         tg_yz_xyyyzzzz_1, tg_yz_xyyzzzz_1, tg_yz_xyyzzzzz_0, tg_yz_xyyzzzzz_1, \
                                         tg_yz_xyzzzzz_1, tg_yz_xyzzzzzz_0, tg_yz_xyzzzzzz_1, tg_yz_xzzzzzz_1, \
                                         tg_yz_xzzzzzzz_0, tg_yz_xzzzzzzz_1, tg_yz_yyyyyyy_1, tg_yz_yyyyyyyy_0, \
                                         tg_yz_yyyyyyyy_1, tg_yz_yyyyyyyz_0, tg_yz_yyyyyyyz_1, tg_yz_yyyyyyz_1, \
                                         tg_yz_yyyyyyzz_0, tg_yz_yyyyyyzz_1, tg_yz_yyyyyzz_1, tg_yz_yyyyyzzz_0, \
                                         tg_yz_yyyyyzzz_1, tg_yz_yyyyzzz_1, tg_yz_yyyyzzzz_0, tg_yz_yyyyzzzz_1, \
                                         tg_yz_yyyzzzz_1, tg_yz_yyyzzzzz_0, tg_yz_yyyzzzzz_1, tg_yz_yyzzzzz_1, \
                                         tg_yz_yyzzzzzz_0, tg_yz_yyzzzzzz_1, tg_yz_yzzzzzz_1, tg_yz_yzzzzzzz_0, \
                                         tg_yz_yzzzzzzz_1, tg_yz_zzzzzzz_1, tg_yz_zzzzzzzz_0, tg_yz_zzzzzzzz_1, \
                                         tg_z_xxxxxxxx_0, tg_z_xxxxxxxx_1, tg_z_xxxxxxxy_0, tg_z_xxxxxxxy_1, tg_z_xxxxxxxz_0, \
                                         tg_z_xxxxxxxz_1, tg_z_xxxxxxyy_0, tg_z_xxxxxxyy_1, tg_z_xxxxxxyz_0, tg_z_xxxxxxyz_1, \
                                         tg_z_xxxxxxzz_0, tg_z_xxxxxxzz_1, tg_z_xxxxxyyy_0, tg_z_xxxxxyyy_1, tg_z_xxxxxyyz_0, \
                                         tg_z_xxxxxyyz_1, tg_z_xxxxxyzz_0, tg_z_xxxxxyzz_1, tg_z_xxxxxzzz_0, tg_z_xxxxxzzz_1, \
                                         tg_z_xxxxyyyy_0, tg_z_xxxxyyyy_1, tg_z_xxxxyyyz_0, tg_z_xxxxyyyz_1, tg_z_xxxxyyzz_0, \
                                         tg_z_xxxxyyzz_1, tg_z_xxxxyzzz_0, tg_z_xxxxyzzz_1, tg_z_xxxxzzzz_0, tg_z_xxxxzzzz_1, \
                                         tg_z_xxxyyyyy_0, tg_z_xxxyyyyy_1, tg_z_xxxyyyyz_0, tg_z_xxxyyyyz_1, tg_z_xxxyyyzz_0, \
                                         tg_z_xxxyyyzz_1, tg_z_xxxyyzzz_0, tg_z_xxxyyzzz_1, tg_z_xxxyzzzz_0, tg_z_xxxyzzzz_1, \
                                         tg_z_xxxzzzzz_0, tg_z_xxxzzzzz_1, tg_z_xxyyyyyy_0, tg_z_xxyyyyyy_1, tg_z_xxyyyyyz_0, \
                                         tg_z_xxyyyyyz_1, tg_z_xxyyyyzz_0, tg_z_xxyyyyzz_1, tg_z_xxyyyzzz_0, tg_z_xxyyyzzz_1, \
                                         tg_z_xxyyzzzz_0, tg_z_xxyyzzzz_1, tg_z_xxyzzzzz_0, tg_z_xxyzzzzz_1, tg_z_xxzzzzzz_0, \
                                         tg_z_xxzzzzzz_1, tg_z_xyyyyyyy_0, tg_z_xyyyyyyy_1, tg_z_xyyyyyyz_0, tg_z_xyyyyyyz_1, \
                                         tg_z_xyyyyyzz_0, tg_z_xyyyyyzz_1, tg_z_xyyyyzzz_0, tg_z_xyyyyzzz_1, tg_z_xyyyzzzz_0, \
                                         tg_z_xyyyzzzz_1, tg_z_xyyzzzzz_0, tg_z_xyyzzzzz_1, tg_z_xyzzzzzz_0, tg_z_xyzzzzzz_1, \
                                         tg_z_xzzzzzzz_0, tg_z_xzzzzzzz_1, tg_z_yyyyyyyy_0, tg_z_yyyyyyyy_1, tg_z_yyyyyyyz_0, \
                                         tg_z_yyyyyyyz_1, tg_z_yyyyyyzz_0, tg_z_yyyyyyzz_1, tg_z_yyyyyzzz_0, tg_z_yyyyyzzz_1, \
                                         tg_z_yyyyzzzz_0, tg_z_yyyyzzzz_1, tg_z_yyyzzzzz_0, tg_z_yyyzzzzz_1, tg_z_yyzzzzzz_0, \
                                         tg_z_yyzzzzzz_1, tg_z_yzzzzzzz_0, tg_z_yzzzzzzz_1, tg_z_zzzzzzzz_0, tg_z_zzzzzzzz_1, \
                                         wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_yyy_xxxxxxxx_0[j] = pb_y * tg_yy_xxxxxxxx_0[j] + wp_y[j] * tg_yy_xxxxxxxx_1[j] + fl1_fx * tg_y_xxxxxxxx_0[j] - fl1_fx * fl1_fza * tg_y_xxxxxxxx_1[j];

                    tg_yyy_xxxxxxxy_0[j] = pb_y * tg_yy_xxxxxxxy_0[j] + wp_y[j] * tg_yy_xxxxxxxy_1[j] + fl1_fx * tg_y_xxxxxxxy_0[j] - fl1_fx * fl1_fza * tg_y_xxxxxxxy_1[j] + 0.5 * fl1_fxn * tg_yy_xxxxxxx_1[j];

                    tg_yyy_xxxxxxxz_0[j] = pb_y * tg_yy_xxxxxxxz_0[j] + wp_y[j] * tg_yy_xxxxxxxz_1[j] + fl1_fx * tg_y_xxxxxxxz_0[j] - fl1_fx * fl1_fza * tg_y_xxxxxxxz_1[j];

                    tg_yyy_xxxxxxyy_0[j] = pb_y * tg_yy_xxxxxxyy_0[j] + wp_y[j] * tg_yy_xxxxxxyy_1[j] + fl1_fx * tg_y_xxxxxxyy_0[j] - fl1_fx * fl1_fza * tg_y_xxxxxxyy_1[j] + fl1_fxn * tg_yy_xxxxxxy_1[j];

                    tg_yyy_xxxxxxyz_0[j] = pb_y * tg_yy_xxxxxxyz_0[j] + wp_y[j] * tg_yy_xxxxxxyz_1[j] + fl1_fx * tg_y_xxxxxxyz_0[j] - fl1_fx * fl1_fza * tg_y_xxxxxxyz_1[j] + 0.5 * fl1_fxn * tg_yy_xxxxxxz_1[j];

                    tg_yyy_xxxxxxzz_0[j] = pb_y * tg_yy_xxxxxxzz_0[j] + wp_y[j] * tg_yy_xxxxxxzz_1[j] + fl1_fx * tg_y_xxxxxxzz_0[j] - fl1_fx * fl1_fza * tg_y_xxxxxxzz_1[j];

                    tg_yyy_xxxxxyyy_0[j] = pb_y * tg_yy_xxxxxyyy_0[j] + wp_y[j] * tg_yy_xxxxxyyy_1[j] + fl1_fx * tg_y_xxxxxyyy_0[j] - fl1_fx * fl1_fza * tg_y_xxxxxyyy_1[j] + 1.5 * fl1_fxn * tg_yy_xxxxxyy_1[j];

                    tg_yyy_xxxxxyyz_0[j] = pb_y * tg_yy_xxxxxyyz_0[j] + wp_y[j] * tg_yy_xxxxxyyz_1[j] + fl1_fx * tg_y_xxxxxyyz_0[j] - fl1_fx * fl1_fza * tg_y_xxxxxyyz_1[j] + fl1_fxn * tg_yy_xxxxxyz_1[j];

                    tg_yyy_xxxxxyzz_0[j] = pb_y * tg_yy_xxxxxyzz_0[j] + wp_y[j] * tg_yy_xxxxxyzz_1[j] + fl1_fx * tg_y_xxxxxyzz_0[j] - fl1_fx * fl1_fza * tg_y_xxxxxyzz_1[j] + 0.5 * fl1_fxn * tg_yy_xxxxxzz_1[j];

                    tg_yyy_xxxxxzzz_0[j] = pb_y * tg_yy_xxxxxzzz_0[j] + wp_y[j] * tg_yy_xxxxxzzz_1[j] + fl1_fx * tg_y_xxxxxzzz_0[j] - fl1_fx * fl1_fza * tg_y_xxxxxzzz_1[j];

                    tg_yyy_xxxxyyyy_0[j] = pb_y * tg_yy_xxxxyyyy_0[j] + wp_y[j] * tg_yy_xxxxyyyy_1[j] + fl1_fx * tg_y_xxxxyyyy_0[j] - fl1_fx * fl1_fza * tg_y_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_yy_xxxxyyy_1[j];

                    tg_yyy_xxxxyyyz_0[j] = pb_y * tg_yy_xxxxyyyz_0[j] + wp_y[j] * tg_yy_xxxxyyyz_1[j] + fl1_fx * tg_y_xxxxyyyz_0[j] - fl1_fx * fl1_fza * tg_y_xxxxyyyz_1[j] + 1.5 * fl1_fxn * tg_yy_xxxxyyz_1[j];

                    tg_yyy_xxxxyyzz_0[j] = pb_y * tg_yy_xxxxyyzz_0[j] + wp_y[j] * tg_yy_xxxxyyzz_1[j] + fl1_fx * tg_y_xxxxyyzz_0[j] - fl1_fx * fl1_fza * tg_y_xxxxyyzz_1[j] + fl1_fxn * tg_yy_xxxxyzz_1[j];

                    tg_yyy_xxxxyzzz_0[j] = pb_y * tg_yy_xxxxyzzz_0[j] + wp_y[j] * tg_yy_xxxxyzzz_1[j] + fl1_fx * tg_y_xxxxyzzz_0[j] - fl1_fx * fl1_fza * tg_y_xxxxyzzz_1[j] + 0.5 * fl1_fxn * tg_yy_xxxxzzz_1[j];

                    tg_yyy_xxxxzzzz_0[j] = pb_y * tg_yy_xxxxzzzz_0[j] + wp_y[j] * tg_yy_xxxxzzzz_1[j] + fl1_fx * tg_y_xxxxzzzz_0[j] - fl1_fx * fl1_fza * tg_y_xxxxzzzz_1[j];

                    tg_yyy_xxxyyyyy_0[j] = pb_y * tg_yy_xxxyyyyy_0[j] + wp_y[j] * tg_yy_xxxyyyyy_1[j] + fl1_fx * tg_y_xxxyyyyy_0[j] - fl1_fx * fl1_fza * tg_y_xxxyyyyy_1[j] + 2.5 * fl1_fxn * tg_yy_xxxyyyy_1[j];

                    tg_yyy_xxxyyyyz_0[j] = pb_y * tg_yy_xxxyyyyz_0[j] + wp_y[j] * tg_yy_xxxyyyyz_1[j] + fl1_fx * tg_y_xxxyyyyz_0[j] - fl1_fx * fl1_fza * tg_y_xxxyyyyz_1[j] + 2.0 * fl1_fxn * tg_yy_xxxyyyz_1[j];

                    tg_yyy_xxxyyyzz_0[j] = pb_y * tg_yy_xxxyyyzz_0[j] + wp_y[j] * tg_yy_xxxyyyzz_1[j] + fl1_fx * tg_y_xxxyyyzz_0[j] - fl1_fx * fl1_fza * tg_y_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_yy_xxxyyzz_1[j];

                    tg_yyy_xxxyyzzz_0[j] = pb_y * tg_yy_xxxyyzzz_0[j] + wp_y[j] * tg_yy_xxxyyzzz_1[j] + fl1_fx * tg_y_xxxyyzzz_0[j] - fl1_fx * fl1_fza * tg_y_xxxyyzzz_1[j] + fl1_fxn * tg_yy_xxxyzzz_1[j];

                    tg_yyy_xxxyzzzz_0[j] = pb_y * tg_yy_xxxyzzzz_0[j] + wp_y[j] * tg_yy_xxxyzzzz_1[j] + fl1_fx * tg_y_xxxyzzzz_0[j] - fl1_fx * fl1_fza * tg_y_xxxyzzzz_1[j] + 0.5 * fl1_fxn * tg_yy_xxxzzzz_1[j];

                    tg_yyy_xxxzzzzz_0[j] = pb_y * tg_yy_xxxzzzzz_0[j] + wp_y[j] * tg_yy_xxxzzzzz_1[j] + fl1_fx * tg_y_xxxzzzzz_0[j] - fl1_fx * fl1_fza * tg_y_xxxzzzzz_1[j];

                    tg_yyy_xxyyyyyy_0[j] = pb_y * tg_yy_xxyyyyyy_0[j] + wp_y[j] * tg_yy_xxyyyyyy_1[j] + fl1_fx * tg_y_xxyyyyyy_0[j] - fl1_fx * fl1_fza * tg_y_xxyyyyyy_1[j] + 3.0 * fl1_fxn * tg_yy_xxyyyyy_1[j];

                    tg_yyy_xxyyyyyz_0[j] = pb_y * tg_yy_xxyyyyyz_0[j] + wp_y[j] * tg_yy_xxyyyyyz_1[j] + fl1_fx * tg_y_xxyyyyyz_0[j] - fl1_fx * fl1_fza * tg_y_xxyyyyyz_1[j] + 2.5 * fl1_fxn * tg_yy_xxyyyyz_1[j];

                    tg_yyy_xxyyyyzz_0[j] = pb_y * tg_yy_xxyyyyzz_0[j] + wp_y[j] * tg_yy_xxyyyyzz_1[j] + fl1_fx * tg_y_xxyyyyzz_0[j] - fl1_fx * fl1_fza * tg_y_xxyyyyzz_1[j] + 2.0 * fl1_fxn * tg_yy_xxyyyzz_1[j];

                    tg_yyy_xxyyyzzz_0[j] = pb_y * tg_yy_xxyyyzzz_0[j] + wp_y[j] * tg_yy_xxyyyzzz_1[j] + fl1_fx * tg_y_xxyyyzzz_0[j] - fl1_fx * fl1_fza * tg_y_xxyyyzzz_1[j] + 1.5 * fl1_fxn * tg_yy_xxyyzzz_1[j];

                    tg_yyy_xxyyzzzz_0[j] = pb_y * tg_yy_xxyyzzzz_0[j] + wp_y[j] * tg_yy_xxyyzzzz_1[j] + fl1_fx * tg_y_xxyyzzzz_0[j] - fl1_fx * fl1_fza * tg_y_xxyyzzzz_1[j] + fl1_fxn * tg_yy_xxyzzzz_1[j];

                    tg_yyy_xxyzzzzz_0[j] = pb_y * tg_yy_xxyzzzzz_0[j] + wp_y[j] * tg_yy_xxyzzzzz_1[j] + fl1_fx * tg_y_xxyzzzzz_0[j] - fl1_fx * fl1_fza * tg_y_xxyzzzzz_1[j] + 0.5 * fl1_fxn * tg_yy_xxzzzzz_1[j];

                    tg_yyy_xxzzzzzz_0[j] = pb_y * tg_yy_xxzzzzzz_0[j] + wp_y[j] * tg_yy_xxzzzzzz_1[j] + fl1_fx * tg_y_xxzzzzzz_0[j] - fl1_fx * fl1_fza * tg_y_xxzzzzzz_1[j];

                    tg_yyy_xyyyyyyy_0[j] = pb_y * tg_yy_xyyyyyyy_0[j] + wp_y[j] * tg_yy_xyyyyyyy_1[j] + fl1_fx * tg_y_xyyyyyyy_0[j] - fl1_fx * fl1_fza * tg_y_xyyyyyyy_1[j] + 3.5 * fl1_fxn * tg_yy_xyyyyyy_1[j];

                    tg_yyy_xyyyyyyz_0[j] = pb_y * tg_yy_xyyyyyyz_0[j] + wp_y[j] * tg_yy_xyyyyyyz_1[j] + fl1_fx * tg_y_xyyyyyyz_0[j] - fl1_fx * fl1_fza * tg_y_xyyyyyyz_1[j] + 3.0 * fl1_fxn * tg_yy_xyyyyyz_1[j];

                    tg_yyy_xyyyyyzz_0[j] = pb_y * tg_yy_xyyyyyzz_0[j] + wp_y[j] * tg_yy_xyyyyyzz_1[j] + fl1_fx * tg_y_xyyyyyzz_0[j] - fl1_fx * fl1_fza * tg_y_xyyyyyzz_1[j] + 2.5 * fl1_fxn * tg_yy_xyyyyzz_1[j];

                    tg_yyy_xyyyyzzz_0[j] = pb_y * tg_yy_xyyyyzzz_0[j] + wp_y[j] * tg_yy_xyyyyzzz_1[j] + fl1_fx * tg_y_xyyyyzzz_0[j] - fl1_fx * fl1_fza * tg_y_xyyyyzzz_1[j] + 2.0 * fl1_fxn * tg_yy_xyyyzzz_1[j];

                    tg_yyy_xyyyzzzz_0[j] = pb_y * tg_yy_xyyyzzzz_0[j] + wp_y[j] * tg_yy_xyyyzzzz_1[j] + fl1_fx * tg_y_xyyyzzzz_0[j] - fl1_fx * fl1_fza * tg_y_xyyyzzzz_1[j] + 1.5 * fl1_fxn * tg_yy_xyyzzzz_1[j];

                    tg_yyy_xyyzzzzz_0[j] = pb_y * tg_yy_xyyzzzzz_0[j] + wp_y[j] * tg_yy_xyyzzzzz_1[j] + fl1_fx * tg_y_xyyzzzzz_0[j] - fl1_fx * fl1_fza * tg_y_xyyzzzzz_1[j] + fl1_fxn * tg_yy_xyzzzzz_1[j];

                    tg_yyy_xyzzzzzz_0[j] = pb_y * tg_yy_xyzzzzzz_0[j] + wp_y[j] * tg_yy_xyzzzzzz_1[j] + fl1_fx * tg_y_xyzzzzzz_0[j] - fl1_fx * fl1_fza * tg_y_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yy_xzzzzzz_1[j];

                    tg_yyy_xzzzzzzz_0[j] = pb_y * tg_yy_xzzzzzzz_0[j] + wp_y[j] * tg_yy_xzzzzzzz_1[j] + fl1_fx * tg_y_xzzzzzzz_0[j] - fl1_fx * fl1_fza * tg_y_xzzzzzzz_1[j];

                    tg_yyy_yyyyyyyy_0[j] = pb_y * tg_yy_yyyyyyyy_0[j] + wp_y[j] * tg_yy_yyyyyyyy_1[j] + fl1_fx * tg_y_yyyyyyyy_0[j] - fl1_fx * fl1_fza * tg_y_yyyyyyyy_1[j] + 4.0 * fl1_fxn * tg_yy_yyyyyyy_1[j];

                    tg_yyy_yyyyyyyz_0[j] = pb_y * tg_yy_yyyyyyyz_0[j] + wp_y[j] * tg_yy_yyyyyyyz_1[j] + fl1_fx * tg_y_yyyyyyyz_0[j] - fl1_fx * fl1_fza * tg_y_yyyyyyyz_1[j] + 3.5 * fl1_fxn * tg_yy_yyyyyyz_1[j];

                    tg_yyy_yyyyyyzz_0[j] = pb_y * tg_yy_yyyyyyzz_0[j] + wp_y[j] * tg_yy_yyyyyyzz_1[j] + fl1_fx * tg_y_yyyyyyzz_0[j] - fl1_fx * fl1_fza * tg_y_yyyyyyzz_1[j] + 3.0 * fl1_fxn * tg_yy_yyyyyzz_1[j];

                    tg_yyy_yyyyyzzz_0[j] = pb_y * tg_yy_yyyyyzzz_0[j] + wp_y[j] * tg_yy_yyyyyzzz_1[j] + fl1_fx * tg_y_yyyyyzzz_0[j] - fl1_fx * fl1_fza * tg_y_yyyyyzzz_1[j] + 2.5 * fl1_fxn * tg_yy_yyyyzzz_1[j];

                    tg_yyy_yyyyzzzz_0[j] = pb_y * tg_yy_yyyyzzzz_0[j] + wp_y[j] * tg_yy_yyyyzzzz_1[j] + fl1_fx * tg_y_yyyyzzzz_0[j] - fl1_fx * fl1_fza * tg_y_yyyyzzzz_1[j] + 2.0 * fl1_fxn * tg_yy_yyyzzzz_1[j];

                    tg_yyy_yyyzzzzz_0[j] = pb_y * tg_yy_yyyzzzzz_0[j] + wp_y[j] * tg_yy_yyyzzzzz_1[j] + fl1_fx * tg_y_yyyzzzzz_0[j] - fl1_fx * fl1_fza * tg_y_yyyzzzzz_1[j] + 1.5 * fl1_fxn * tg_yy_yyzzzzz_1[j];

                    tg_yyy_yyzzzzzz_0[j] = pb_y * tg_yy_yyzzzzzz_0[j] + wp_y[j] * tg_yy_yyzzzzzz_1[j] + fl1_fx * tg_y_yyzzzzzz_0[j] - fl1_fx * fl1_fza * tg_y_yyzzzzzz_1[j] + fl1_fxn * tg_yy_yzzzzzz_1[j];

                    tg_yyy_yzzzzzzz_0[j] = pb_y * tg_yy_yzzzzzzz_0[j] + wp_y[j] * tg_yy_yzzzzzzz_1[j] + fl1_fx * tg_y_yzzzzzzz_0[j] - fl1_fx * fl1_fza * tg_y_yzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yy_zzzzzzz_1[j];

                    tg_yyy_zzzzzzzz_0[j] = pb_y * tg_yy_zzzzzzzz_0[j] + wp_y[j] * tg_yy_zzzzzzzz_1[j] + fl1_fx * tg_y_zzzzzzzz_0[j] - fl1_fx * fl1_fza * tg_y_zzzzzzzz_1[j];

                    tg_yyz_xxxxxxxx_0[j] = pb_y * tg_yz_xxxxxxxx_0[j] + wp_y[j] * tg_yz_xxxxxxxx_1[j] + 0.5 * fl1_fx * tg_z_xxxxxxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxxxx_1[j];

                    tg_yyz_xxxxxxxy_0[j] = pb_y * tg_yz_xxxxxxxy_0[j] + wp_y[j] * tg_yz_xxxxxxxy_1[j] + 0.5 * fl1_fx * tg_z_xxxxxxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxxxy_1[j] + 0.5 * fl1_fxn * tg_yz_xxxxxxx_1[j];

                    tg_yyz_xxxxxxxz_0[j] = pb_y * tg_yz_xxxxxxxz_0[j] + wp_y[j] * tg_yz_xxxxxxxz_1[j] + 0.5 * fl1_fx * tg_z_xxxxxxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxxxz_1[j];

                    tg_yyz_xxxxxxyy_0[j] = pb_y * tg_yz_xxxxxxyy_0[j] + wp_y[j] * tg_yz_xxxxxxyy_1[j] + 0.5 * fl1_fx * tg_z_xxxxxxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxxyy_1[j] + fl1_fxn * tg_yz_xxxxxxy_1[j];

                    tg_yyz_xxxxxxyz_0[j] = pb_y * tg_yz_xxxxxxyz_0[j] + wp_y[j] * tg_yz_xxxxxxyz_1[j] + 0.5 * fl1_fx * tg_z_xxxxxxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxxyz_1[j] + 0.5 * fl1_fxn * tg_yz_xxxxxxz_1[j];

                    tg_yyz_xxxxxxzz_0[j] = pb_y * tg_yz_xxxxxxzz_0[j] + wp_y[j] * tg_yz_xxxxxxzz_1[j] + 0.5 * fl1_fx * tg_z_xxxxxxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxxzz_1[j];

                    tg_yyz_xxxxxyyy_0[j] = pb_y * tg_yz_xxxxxyyy_0[j] + wp_y[j] * tg_yz_xxxxxyyy_1[j] + 0.5 * fl1_fx * tg_z_xxxxxyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxyyy_1[j] + 1.5 * fl1_fxn * tg_yz_xxxxxyy_1[j];

                    tg_yyz_xxxxxyyz_0[j] = pb_y * tg_yz_xxxxxyyz_0[j] + wp_y[j] * tg_yz_xxxxxyyz_1[j] + 0.5 * fl1_fx * tg_z_xxxxxyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxyyz_1[j] + fl1_fxn * tg_yz_xxxxxyz_1[j];

                    tg_yyz_xxxxxyzz_0[j] = pb_y * tg_yz_xxxxxyzz_0[j] + wp_y[j] * tg_yz_xxxxxyzz_1[j] + 0.5 * fl1_fx * tg_z_xxxxxyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxyzz_1[j] + 0.5 * fl1_fxn * tg_yz_xxxxxzz_1[j];

                    tg_yyz_xxxxxzzz_0[j] = pb_y * tg_yz_xxxxxzzz_0[j] + wp_y[j] * tg_yz_xxxxxzzz_1[j] + 0.5 * fl1_fx * tg_z_xxxxxzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxxzzz_1[j];

                    tg_yyz_xxxxyyyy_0[j] = pb_y * tg_yz_xxxxyyyy_0[j] + wp_y[j] * tg_yz_xxxxyyyy_1[j] + 0.5 * fl1_fx * tg_z_xxxxyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_yz_xxxxyyy_1[j];

                    tg_yyz_xxxxyyyz_0[j] = pb_y * tg_yz_xxxxyyyz_0[j] + wp_y[j] * tg_yz_xxxxyyyz_1[j] + 0.5 * fl1_fx * tg_z_xxxxyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxyyyz_1[j] + 1.5 * fl1_fxn * tg_yz_xxxxyyz_1[j];

                    tg_yyz_xxxxyyzz_0[j] = pb_y * tg_yz_xxxxyyzz_0[j] + wp_y[j] * tg_yz_xxxxyyzz_1[j] + 0.5 * fl1_fx * tg_z_xxxxyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxyyzz_1[j] + fl1_fxn * tg_yz_xxxxyzz_1[j];

                    tg_yyz_xxxxyzzz_0[j] = pb_y * tg_yz_xxxxyzzz_0[j] + wp_y[j] * tg_yz_xxxxyzzz_1[j] + 0.5 * fl1_fx * tg_z_xxxxyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxyzzz_1[j] + 0.5 * fl1_fxn * tg_yz_xxxxzzz_1[j];

                    tg_yyz_xxxxzzzz_0[j] = pb_y * tg_yz_xxxxzzzz_0[j] + wp_y[j] * tg_yz_xxxxzzzz_1[j] + 0.5 * fl1_fx * tg_z_xxxxzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxzzzz_1[j];

                    tg_yyz_xxxyyyyy_0[j] = pb_y * tg_yz_xxxyyyyy_0[j] + wp_y[j] * tg_yz_xxxyyyyy_1[j] + 0.5 * fl1_fx * tg_z_xxxyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxyyyyy_1[j] + 2.5 * fl1_fxn * tg_yz_xxxyyyy_1[j];

                    tg_yyz_xxxyyyyz_0[j] = pb_y * tg_yz_xxxyyyyz_0[j] + wp_y[j] * tg_yz_xxxyyyyz_1[j] + 0.5 * fl1_fx * tg_z_xxxyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxyyyyz_1[j] + 2.0 * fl1_fxn * tg_yz_xxxyyyz_1[j];

                    tg_yyz_xxxyyyzz_0[j] = pb_y * tg_yz_xxxyyyzz_0[j] + wp_y[j] * tg_yz_xxxyyyzz_1[j] + 0.5 * fl1_fx * tg_z_xxxyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_yz_xxxyyzz_1[j];

                    tg_yyz_xxxyyzzz_0[j] = pb_y * tg_yz_xxxyyzzz_0[j] + wp_y[j] * tg_yz_xxxyyzzz_1[j] + 0.5 * fl1_fx * tg_z_xxxyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxyyzzz_1[j] + fl1_fxn * tg_yz_xxxyzzz_1[j];

                    tg_yyz_xxxyzzzz_0[j] = pb_y * tg_yz_xxxyzzzz_0[j] + wp_y[j] * tg_yz_xxxyzzzz_1[j] + 0.5 * fl1_fx * tg_z_xxxyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxyzzzz_1[j] + 0.5 * fl1_fxn * tg_yz_xxxzzzz_1[j];

                    tg_yyz_xxxzzzzz_0[j] = pb_y * tg_yz_xxxzzzzz_0[j] + wp_y[j] * tg_yz_xxxzzzzz_1[j] + 0.5 * fl1_fx * tg_z_xxxzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxzzzzz_1[j];

                    tg_yyz_xxyyyyyy_0[j] = pb_y * tg_yz_xxyyyyyy_0[j] + wp_y[j] * tg_yz_xxyyyyyy_1[j] + 0.5 * fl1_fx * tg_z_xxyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxyyyyyy_1[j] + 3.0 * fl1_fxn * tg_yz_xxyyyyy_1[j];

                    tg_yyz_xxyyyyyz_0[j] = pb_y * tg_yz_xxyyyyyz_0[j] + wp_y[j] * tg_yz_xxyyyyyz_1[j] + 0.5 * fl1_fx * tg_z_xxyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxyyyyyz_1[j] + 2.5 * fl1_fxn * tg_yz_xxyyyyz_1[j];

                    tg_yyz_xxyyyyzz_0[j] = pb_y * tg_yz_xxyyyyzz_0[j] + wp_y[j] * tg_yz_xxyyyyzz_1[j] + 0.5 * fl1_fx * tg_z_xxyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxyyyyzz_1[j] + 2.0 * fl1_fxn * tg_yz_xxyyyzz_1[j];

                    tg_yyz_xxyyyzzz_0[j] = pb_y * tg_yz_xxyyyzzz_0[j] + wp_y[j] * tg_yz_xxyyyzzz_1[j] + 0.5 * fl1_fx * tg_z_xxyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxyyyzzz_1[j] + 1.5 * fl1_fxn * tg_yz_xxyyzzz_1[j];

                    tg_yyz_xxyyzzzz_0[j] = pb_y * tg_yz_xxyyzzzz_0[j] + wp_y[j] * tg_yz_xxyyzzzz_1[j] + 0.5 * fl1_fx * tg_z_xxyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxyyzzzz_1[j] + fl1_fxn * tg_yz_xxyzzzz_1[j];

                    tg_yyz_xxyzzzzz_0[j] = pb_y * tg_yz_xxyzzzzz_0[j] + wp_y[j] * tg_yz_xxyzzzzz_1[j] + 0.5 * fl1_fx * tg_z_xxyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxyzzzzz_1[j] + 0.5 * fl1_fxn * tg_yz_xxzzzzz_1[j];

                    tg_yyz_xxzzzzzz_0[j] = pb_y * tg_yz_xxzzzzzz_0[j] + wp_y[j] * tg_yz_xxzzzzzz_1[j] + 0.5 * fl1_fx * tg_z_xxzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxzzzzzz_1[j];

                    tg_yyz_xyyyyyyy_0[j] = pb_y * tg_yz_xyyyyyyy_0[j] + wp_y[j] * tg_yz_xyyyyyyy_1[j] + 0.5 * fl1_fx * tg_z_xyyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyyyyyyy_1[j] + 3.5 * fl1_fxn * tg_yz_xyyyyyy_1[j];

                    tg_yyz_xyyyyyyz_0[j] = pb_y * tg_yz_xyyyyyyz_0[j] + wp_y[j] * tg_yz_xyyyyyyz_1[j] + 0.5 * fl1_fx * tg_z_xyyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyyyyyyz_1[j] + 3.0 * fl1_fxn * tg_yz_xyyyyyz_1[j];

                    tg_yyz_xyyyyyzz_0[j] = pb_y * tg_yz_xyyyyyzz_0[j] + wp_y[j] * tg_yz_xyyyyyzz_1[j] + 0.5 * fl1_fx * tg_z_xyyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyyyyyzz_1[j] + 2.5 * fl1_fxn * tg_yz_xyyyyzz_1[j];

                    tg_yyz_xyyyyzzz_0[j] = pb_y * tg_yz_xyyyyzzz_0[j] + wp_y[j] * tg_yz_xyyyyzzz_1[j] + 0.5 * fl1_fx * tg_z_xyyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyyyyzzz_1[j] + 2.0 * fl1_fxn * tg_yz_xyyyzzz_1[j];

                    tg_yyz_xyyyzzzz_0[j] = pb_y * tg_yz_xyyyzzzz_0[j] + wp_y[j] * tg_yz_xyyyzzzz_1[j] + 0.5 * fl1_fx * tg_z_xyyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyyyzzzz_1[j] + 1.5 * fl1_fxn * tg_yz_xyyzzzz_1[j];

                    tg_yyz_xyyzzzzz_0[j] = pb_y * tg_yz_xyyzzzzz_0[j] + wp_y[j] * tg_yz_xyyzzzzz_1[j] + 0.5 * fl1_fx * tg_z_xyyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyyzzzzz_1[j] + fl1_fxn * tg_yz_xyzzzzz_1[j];

                    tg_yyz_xyzzzzzz_0[j] = pb_y * tg_yz_xyzzzzzz_0[j] + wp_y[j] * tg_yz_xyzzzzzz_1[j] + 0.5 * fl1_fx * tg_z_xyzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yz_xzzzzzz_1[j];

                    tg_yyz_xzzzzzzz_0[j] = pb_y * tg_yz_xzzzzzzz_0[j] + wp_y[j] * tg_yz_xzzzzzzz_1[j] + 0.5 * fl1_fx * tg_z_xzzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xzzzzzzz_1[j];

                    tg_yyz_yyyyyyyy_0[j] = pb_y * tg_yz_yyyyyyyy_0[j] + wp_y[j] * tg_yz_yyyyyyyy_1[j] + 0.5 * fl1_fx * tg_z_yyyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyyyyyyy_1[j] + 4.0 * fl1_fxn * tg_yz_yyyyyyy_1[j];

                    tg_yyz_yyyyyyyz_0[j] = pb_y * tg_yz_yyyyyyyz_0[j] + wp_y[j] * tg_yz_yyyyyyyz_1[j] + 0.5 * fl1_fx * tg_z_yyyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyyyyyyz_1[j] + 3.5 * fl1_fxn * tg_yz_yyyyyyz_1[j];

                    tg_yyz_yyyyyyzz_0[j] = pb_y * tg_yz_yyyyyyzz_0[j] + wp_y[j] * tg_yz_yyyyyyzz_1[j] + 0.5 * fl1_fx * tg_z_yyyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyyyyyzz_1[j] + 3.0 * fl1_fxn * tg_yz_yyyyyzz_1[j];

                    tg_yyz_yyyyyzzz_0[j] = pb_y * tg_yz_yyyyyzzz_0[j] + wp_y[j] * tg_yz_yyyyyzzz_1[j] + 0.5 * fl1_fx * tg_z_yyyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyyyyzzz_1[j] + 2.5 * fl1_fxn * tg_yz_yyyyzzz_1[j];

                    tg_yyz_yyyyzzzz_0[j] = pb_y * tg_yz_yyyyzzzz_0[j] + wp_y[j] * tg_yz_yyyyzzzz_1[j] + 0.5 * fl1_fx * tg_z_yyyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyyyzzzz_1[j] + 2.0 * fl1_fxn * tg_yz_yyyzzzz_1[j];

                    tg_yyz_yyyzzzzz_0[j] = pb_y * tg_yz_yyyzzzzz_0[j] + wp_y[j] * tg_yz_yyyzzzzz_1[j] + 0.5 * fl1_fx * tg_z_yyyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyyzzzzz_1[j] + 1.5 * fl1_fxn * tg_yz_yyzzzzz_1[j];

                    tg_yyz_yyzzzzzz_0[j] = pb_y * tg_yz_yyzzzzzz_0[j] + wp_y[j] * tg_yz_yyzzzzzz_1[j] + 0.5 * fl1_fx * tg_z_yyzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyzzzzzz_1[j] + fl1_fxn * tg_yz_yzzzzzz_1[j];

                    tg_yyz_yzzzzzzz_0[j] = pb_y * tg_yz_yzzzzzzz_0[j] + wp_y[j] * tg_yz_yzzzzzzz_1[j] + 0.5 * fl1_fx * tg_z_yzzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yz_zzzzzzz_1[j];

                    tg_yyz_zzzzzzzz_0[j] = pb_y * tg_yz_zzzzzzzz_0[j] + wp_y[j] * tg_yz_zzzzzzzz_1[j] + 0.5 * fl1_fx * tg_z_zzzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_zzzzzzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSFSL_360_450(      CMemBlock2D<double>& primBuffer,
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

        auto r_pb_z = braGtoPairsBlock.getDistancesPBZ();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {3, -1, -1, -1},
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_3_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_3_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_2_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_2_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_2_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_1_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_1_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {8, -1, -1, -1}, 
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

                auto tg_zz_xxxxxxxx_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 225); 

                auto tg_zz_xxxxxxxy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 226); 

                auto tg_zz_xxxxxxxz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 227); 

                auto tg_zz_xxxxxxyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 228); 

                auto tg_zz_xxxxxxyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 229); 

                auto tg_zz_xxxxxxzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 230); 

                auto tg_zz_xxxxxyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 231); 

                auto tg_zz_xxxxxyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 232); 

                auto tg_zz_xxxxxyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 233); 

                auto tg_zz_xxxxxzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 234); 

                auto tg_zz_xxxxyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 235); 

                auto tg_zz_xxxxyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 236); 

                auto tg_zz_xxxxyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 237); 

                auto tg_zz_xxxxyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 238); 

                auto tg_zz_xxxxzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 239); 

                auto tg_zz_xxxyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 240); 

                auto tg_zz_xxxyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 241); 

                auto tg_zz_xxxyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 242); 

                auto tg_zz_xxxyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 243); 

                auto tg_zz_xxxyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 244); 

                auto tg_zz_xxxzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 245); 

                auto tg_zz_xxyyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 246); 

                auto tg_zz_xxyyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 247); 

                auto tg_zz_xxyyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 248); 

                auto tg_zz_xxyyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 249); 

                auto tg_zz_xxyyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 250); 

                auto tg_zz_xxyzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 251); 

                auto tg_zz_xxzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 252); 

                auto tg_zz_xyyyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 253); 

                auto tg_zz_xyyyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 254); 

                auto tg_zz_xyyyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 255); 

                auto tg_zz_xyyyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 256); 

                auto tg_zz_xyyyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 257); 

                auto tg_zz_xyyzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 258); 

                auto tg_zz_xyzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 259); 

                auto tg_zz_xzzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 260); 

                auto tg_zz_yyyyyyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 261); 

                auto tg_zz_yyyyyyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 262); 

                auto tg_zz_yyyyyyzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 263); 

                auto tg_zz_yyyyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 264); 

                auto tg_zz_yyyyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 265); 

                auto tg_zz_yyyzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 266); 

                auto tg_zz_yyzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 267); 

                auto tg_zz_yzzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 268); 

                auto tg_zz_zzzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 269); 

                auto tg_zz_xxxxxxxx_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 225); 

                auto tg_zz_xxxxxxxy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 226); 

                auto tg_zz_xxxxxxxz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 227); 

                auto tg_zz_xxxxxxyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 228); 

                auto tg_zz_xxxxxxyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 229); 

                auto tg_zz_xxxxxxzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 230); 

                auto tg_zz_xxxxxyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 231); 

                auto tg_zz_xxxxxyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 232); 

                auto tg_zz_xxxxxyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 233); 

                auto tg_zz_xxxxxzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 234); 

                auto tg_zz_xxxxyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 235); 

                auto tg_zz_xxxxyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 236); 

                auto tg_zz_xxxxyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 237); 

                auto tg_zz_xxxxyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 238); 

                auto tg_zz_xxxxzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 239); 

                auto tg_zz_xxxyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 240); 

                auto tg_zz_xxxyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 241); 

                auto tg_zz_xxxyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 242); 

                auto tg_zz_xxxyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 243); 

                auto tg_zz_xxxyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 244); 

                auto tg_zz_xxxzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 245); 

                auto tg_zz_xxyyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 246); 

                auto tg_zz_xxyyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 247); 

                auto tg_zz_xxyyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 248); 

                auto tg_zz_xxyyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 249); 

                auto tg_zz_xxyyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 250); 

                auto tg_zz_xxyzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 251); 

                auto tg_zz_xxzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 252); 

                auto tg_zz_xyyyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 253); 

                auto tg_zz_xyyyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 254); 

                auto tg_zz_xyyyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 255); 

                auto tg_zz_xyyyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 256); 

                auto tg_zz_xyyyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 257); 

                auto tg_zz_xyyzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 258); 

                auto tg_zz_xyzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 259); 

                auto tg_zz_xzzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 260); 

                auto tg_zz_yyyyyyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 261); 

                auto tg_zz_yyyyyyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 262); 

                auto tg_zz_yyyyyyzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 263); 

                auto tg_zz_yyyyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 264); 

                auto tg_zz_yyyyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 265); 

                auto tg_zz_yyyzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 266); 

                auto tg_zz_yyzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 267); 

                auto tg_zz_yzzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 268); 

                auto tg_zz_zzzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 269); 

                auto tg_z_xxxxxxxx_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 90); 

                auto tg_z_xxxxxxxy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 91); 

                auto tg_z_xxxxxxxz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 92); 

                auto tg_z_xxxxxxyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 93); 

                auto tg_z_xxxxxxyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 94); 

                auto tg_z_xxxxxxzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 95); 

                auto tg_z_xxxxxyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 96); 

                auto tg_z_xxxxxyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 97); 

                auto tg_z_xxxxxyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 98); 

                auto tg_z_xxxxxzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 99); 

                auto tg_z_xxxxyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 100); 

                auto tg_z_xxxxyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 101); 

                auto tg_z_xxxxyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 102); 

                auto tg_z_xxxxyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 103); 

                auto tg_z_xxxxzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 104); 

                auto tg_z_xxxyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 105); 

                auto tg_z_xxxyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 106); 

                auto tg_z_xxxyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 107); 

                auto tg_z_xxxyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 108); 

                auto tg_z_xxxyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 109); 

                auto tg_z_xxxzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 110); 

                auto tg_z_xxyyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 111); 

                auto tg_z_xxyyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 112); 

                auto tg_z_xxyyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 113); 

                auto tg_z_xxyyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 114); 

                auto tg_z_xxyyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 115); 

                auto tg_z_xxyzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 116); 

                auto tg_z_xxzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 117); 

                auto tg_z_xyyyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 118); 

                auto tg_z_xyyyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 119); 

                auto tg_z_xyyyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 120); 

                auto tg_z_xyyyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 121); 

                auto tg_z_xyyyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 122); 

                auto tg_z_xyyzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 123); 

                auto tg_z_xyzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 124); 

                auto tg_z_xzzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 125); 

                auto tg_z_yyyyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 126); 

                auto tg_z_yyyyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 127); 

                auto tg_z_yyyyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 128); 

                auto tg_z_yyyyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 129); 

                auto tg_z_yyyyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 130); 

                auto tg_z_yyyzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 131); 

                auto tg_z_yyzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 132); 

                auto tg_z_yzzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 133); 

                auto tg_z_zzzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 134); 

                auto tg_z_xxxxxxxx_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 90); 

                auto tg_z_xxxxxxxy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 91); 

                auto tg_z_xxxxxxxz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 92); 

                auto tg_z_xxxxxxyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 93); 

                auto tg_z_xxxxxxyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 94); 

                auto tg_z_xxxxxxzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 95); 

                auto tg_z_xxxxxyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 96); 

                auto tg_z_xxxxxyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 97); 

                auto tg_z_xxxxxyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 98); 

                auto tg_z_xxxxxzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 99); 

                auto tg_z_xxxxyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 100); 

                auto tg_z_xxxxyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 101); 

                auto tg_z_xxxxyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 102); 

                auto tg_z_xxxxyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 103); 

                auto tg_z_xxxxzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 104); 

                auto tg_z_xxxyyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 105); 

                auto tg_z_xxxyyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 106); 

                auto tg_z_xxxyyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 107); 

                auto tg_z_xxxyyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 108); 

                auto tg_z_xxxyzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 109); 

                auto tg_z_xxxzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 110); 

                auto tg_z_xxyyyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 111); 

                auto tg_z_xxyyyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 112); 

                auto tg_z_xxyyyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 113); 

                auto tg_z_xxyyyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 114); 

                auto tg_z_xxyyzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 115); 

                auto tg_z_xxyzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 116); 

                auto tg_z_xxzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 117); 

                auto tg_z_xyyyyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 118); 

                auto tg_z_xyyyyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 119); 

                auto tg_z_xyyyyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 120); 

                auto tg_z_xyyyyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 121); 

                auto tg_z_xyyyzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 122); 

                auto tg_z_xyyzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 123); 

                auto tg_z_xyzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 124); 

                auto tg_z_xzzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 125); 

                auto tg_z_yyyyyyyy_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 126); 

                auto tg_z_yyyyyyyz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 127); 

                auto tg_z_yyyyyyzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 128); 

                auto tg_z_yyyyyzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 129); 

                auto tg_z_yyyyzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 130); 

                auto tg_z_yyyzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 131); 

                auto tg_z_yyzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 132); 

                auto tg_z_yzzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 133); 

                auto tg_z_zzzzzzzz_1 = primBuffer.data(pidx_g_1_8_m1 + 135 * idx + 134); 

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

                // set up pointers to integrals

                auto tg_yzz_xxxxxxxx_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 360); 

                auto tg_yzz_xxxxxxxy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 361); 

                auto tg_yzz_xxxxxxxz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 362); 

                auto tg_yzz_xxxxxxyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 363); 

                auto tg_yzz_xxxxxxyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 364); 

                auto tg_yzz_xxxxxxzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 365); 

                auto tg_yzz_xxxxxyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 366); 

                auto tg_yzz_xxxxxyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 367); 

                auto tg_yzz_xxxxxyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 368); 

                auto tg_yzz_xxxxxzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 369); 

                auto tg_yzz_xxxxyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 370); 

                auto tg_yzz_xxxxyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 371); 

                auto tg_yzz_xxxxyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 372); 

                auto tg_yzz_xxxxyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 373); 

                auto tg_yzz_xxxxzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 374); 

                auto tg_yzz_xxxyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 375); 

                auto tg_yzz_xxxyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 376); 

                auto tg_yzz_xxxyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 377); 

                auto tg_yzz_xxxyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 378); 

                auto tg_yzz_xxxyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 379); 

                auto tg_yzz_xxxzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 380); 

                auto tg_yzz_xxyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 381); 

                auto tg_yzz_xxyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 382); 

                auto tg_yzz_xxyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 383); 

                auto tg_yzz_xxyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 384); 

                auto tg_yzz_xxyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 385); 

                auto tg_yzz_xxyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 386); 

                auto tg_yzz_xxzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 387); 

                auto tg_yzz_xyyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 388); 

                auto tg_yzz_xyyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 389); 

                auto tg_yzz_xyyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 390); 

                auto tg_yzz_xyyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 391); 

                auto tg_yzz_xyyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 392); 

                auto tg_yzz_xyyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 393); 

                auto tg_yzz_xyzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 394); 

                auto tg_yzz_xzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 395); 

                auto tg_yzz_yyyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 396); 

                auto tg_yzz_yyyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 397); 

                auto tg_yzz_yyyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 398); 

                auto tg_yzz_yyyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 399); 

                auto tg_yzz_yyyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 400); 

                auto tg_yzz_yyyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 401); 

                auto tg_yzz_yyzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 402); 

                auto tg_yzz_yzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 403); 

                auto tg_yzz_zzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 404); 

                auto tg_zzz_xxxxxxxx_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 405); 

                auto tg_zzz_xxxxxxxy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 406); 

                auto tg_zzz_xxxxxxxz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 407); 

                auto tg_zzz_xxxxxxyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 408); 

                auto tg_zzz_xxxxxxyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 409); 

                auto tg_zzz_xxxxxxzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 410); 

                auto tg_zzz_xxxxxyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 411); 

                auto tg_zzz_xxxxxyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 412); 

                auto tg_zzz_xxxxxyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 413); 

                auto tg_zzz_xxxxxzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 414); 

                auto tg_zzz_xxxxyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 415); 

                auto tg_zzz_xxxxyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 416); 

                auto tg_zzz_xxxxyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 417); 

                auto tg_zzz_xxxxyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 418); 

                auto tg_zzz_xxxxzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 419); 

                auto tg_zzz_xxxyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 420); 

                auto tg_zzz_xxxyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 421); 

                auto tg_zzz_xxxyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 422); 

                auto tg_zzz_xxxyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 423); 

                auto tg_zzz_xxxyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 424); 

                auto tg_zzz_xxxzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 425); 

                auto tg_zzz_xxyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 426); 

                auto tg_zzz_xxyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 427); 

                auto tg_zzz_xxyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 428); 

                auto tg_zzz_xxyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 429); 

                auto tg_zzz_xxyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 430); 

                auto tg_zzz_xxyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 431); 

                auto tg_zzz_xxzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 432); 

                auto tg_zzz_xyyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 433); 

                auto tg_zzz_xyyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 434); 

                auto tg_zzz_xyyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 435); 

                auto tg_zzz_xyyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 436); 

                auto tg_zzz_xyyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 437); 

                auto tg_zzz_xyyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 438); 

                auto tg_zzz_xyzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 439); 

                auto tg_zzz_xzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 440); 

                auto tg_zzz_yyyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 441); 

                auto tg_zzz_yyyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 442); 

                auto tg_zzz_yyyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 443); 

                auto tg_zzz_yyyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 444); 

                auto tg_zzz_yyyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 445); 

                auto tg_zzz_yyyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 446); 

                auto tg_zzz_yyzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 447); 

                auto tg_zzz_yzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 448); 

                auto tg_zzz_zzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 449); 

                // Batch of Integrals (360,450)

                #pragma omp simd aligned(fxn, fza, tg_yzz_xxxxxxxx_0, tg_yzz_xxxxxxxy_0, tg_yzz_xxxxxxxz_0, \
                                         tg_yzz_xxxxxxyy_0, tg_yzz_xxxxxxyz_0, tg_yzz_xxxxxxzz_0, tg_yzz_xxxxxyyy_0, \
                                         tg_yzz_xxxxxyyz_0, tg_yzz_xxxxxyzz_0, tg_yzz_xxxxxzzz_0, tg_yzz_xxxxyyyy_0, \
                                         tg_yzz_xxxxyyyz_0, tg_yzz_xxxxyyzz_0, tg_yzz_xxxxyzzz_0, tg_yzz_xxxxzzzz_0, \
                                         tg_yzz_xxxyyyyy_0, tg_yzz_xxxyyyyz_0, tg_yzz_xxxyyyzz_0, tg_yzz_xxxyyzzz_0, \
                                         tg_yzz_xxxyzzzz_0, tg_yzz_xxxzzzzz_0, tg_yzz_xxyyyyyy_0, tg_yzz_xxyyyyyz_0, \
                                         tg_yzz_xxyyyyzz_0, tg_yzz_xxyyyzzz_0, tg_yzz_xxyyzzzz_0, tg_yzz_xxyzzzzz_0, \
                                         tg_yzz_xxzzzzzz_0, tg_yzz_xyyyyyyy_0, tg_yzz_xyyyyyyz_0, tg_yzz_xyyyyyzz_0, \
                                         tg_yzz_xyyyyzzz_0, tg_yzz_xyyyzzzz_0, tg_yzz_xyyzzzzz_0, tg_yzz_xyzzzzzz_0, \
                                         tg_yzz_xzzzzzzz_0, tg_yzz_yyyyyyyy_0, tg_yzz_yyyyyyyz_0, tg_yzz_yyyyyyzz_0, \
                                         tg_yzz_yyyyyzzz_0, tg_yzz_yyyyzzzz_0, tg_yzz_yyyzzzzz_0, tg_yzz_yyzzzzzz_0, \
                                         tg_yzz_yzzzzzzz_0, tg_yzz_zzzzzzzz_0, tg_z_xxxxxxxx_0, tg_z_xxxxxxxx_1, \
                                         tg_z_xxxxxxxy_0, tg_z_xxxxxxxy_1, tg_z_xxxxxxxz_0, tg_z_xxxxxxxz_1, tg_z_xxxxxxyy_0, \
                                         tg_z_xxxxxxyy_1, tg_z_xxxxxxyz_0, tg_z_xxxxxxyz_1, tg_z_xxxxxxzz_0, tg_z_xxxxxxzz_1, \
                                         tg_z_xxxxxyyy_0, tg_z_xxxxxyyy_1, tg_z_xxxxxyyz_0, tg_z_xxxxxyyz_1, tg_z_xxxxxyzz_0, \
                                         tg_z_xxxxxyzz_1, tg_z_xxxxxzzz_0, tg_z_xxxxxzzz_1, tg_z_xxxxyyyy_0, tg_z_xxxxyyyy_1, \
                                         tg_z_xxxxyyyz_0, tg_z_xxxxyyyz_1, tg_z_xxxxyyzz_0, tg_z_xxxxyyzz_1, tg_z_xxxxyzzz_0, \
                                         tg_z_xxxxyzzz_1, tg_z_xxxxzzzz_0, tg_z_xxxxzzzz_1, tg_z_xxxyyyyy_0, tg_z_xxxyyyyy_1, \
                                         tg_z_xxxyyyyz_0, tg_z_xxxyyyyz_1, tg_z_xxxyyyzz_0, tg_z_xxxyyyzz_1, tg_z_xxxyyzzz_0, \
                                         tg_z_xxxyyzzz_1, tg_z_xxxyzzzz_0, tg_z_xxxyzzzz_1, tg_z_xxxzzzzz_0, tg_z_xxxzzzzz_1, \
                                         tg_z_xxyyyyyy_0, tg_z_xxyyyyyy_1, tg_z_xxyyyyyz_0, tg_z_xxyyyyyz_1, tg_z_xxyyyyzz_0, \
                                         tg_z_xxyyyyzz_1, tg_z_xxyyyzzz_0, tg_z_xxyyyzzz_1, tg_z_xxyyzzzz_0, tg_z_xxyyzzzz_1, \
                                         tg_z_xxyzzzzz_0, tg_z_xxyzzzzz_1, tg_z_xxzzzzzz_0, tg_z_xxzzzzzz_1, tg_z_xyyyyyyy_0, \
                                         tg_z_xyyyyyyy_1, tg_z_xyyyyyyz_0, tg_z_xyyyyyyz_1, tg_z_xyyyyyzz_0, tg_z_xyyyyyzz_1, \
                                         tg_z_xyyyyzzz_0, tg_z_xyyyyzzz_1, tg_z_xyyyzzzz_0, tg_z_xyyyzzzz_1, tg_z_xyyzzzzz_0, \
                                         tg_z_xyyzzzzz_1, tg_z_xyzzzzzz_0, tg_z_xyzzzzzz_1, tg_z_xzzzzzzz_0, tg_z_xzzzzzzz_1, \
                                         tg_z_yyyyyyyy_0, tg_z_yyyyyyyy_1, tg_z_yyyyyyyz_0, tg_z_yyyyyyyz_1, tg_z_yyyyyyzz_0, \
                                         tg_z_yyyyyyzz_1, tg_z_yyyyyzzz_0, tg_z_yyyyyzzz_1, tg_z_yyyyzzzz_0, tg_z_yyyyzzzz_1, \
                                         tg_z_yyyzzzzz_0, tg_z_yyyzzzzz_1, tg_z_yyzzzzzz_0, tg_z_yyzzzzzz_1, tg_z_yzzzzzzz_0, \
                                         tg_z_yzzzzzzz_1, tg_z_zzzzzzzz_0, tg_z_zzzzzzzz_1, tg_zz_xxxxxxx_1, \
                                         tg_zz_xxxxxxxx_0, tg_zz_xxxxxxxx_1, tg_zz_xxxxxxxy_0, tg_zz_xxxxxxxy_1, \
                                         tg_zz_xxxxxxxz_0, tg_zz_xxxxxxxz_1, tg_zz_xxxxxxy_1, tg_zz_xxxxxxyy_0, \
                                         tg_zz_xxxxxxyy_1, tg_zz_xxxxxxyz_0, tg_zz_xxxxxxyz_1, tg_zz_xxxxxxz_1, \
                                         tg_zz_xxxxxxzz_0, tg_zz_xxxxxxzz_1, tg_zz_xxxxxyy_1, tg_zz_xxxxxyyy_0, \
                                         tg_zz_xxxxxyyy_1, tg_zz_xxxxxyyz_0, tg_zz_xxxxxyyz_1, tg_zz_xxxxxyz_1, \
                                         tg_zz_xxxxxyzz_0, tg_zz_xxxxxyzz_1, tg_zz_xxxxxzz_1, tg_zz_xxxxxzzz_0, \
                                         tg_zz_xxxxxzzz_1, tg_zz_xxxxyyy_1, tg_zz_xxxxyyyy_0, tg_zz_xxxxyyyy_1, \
                                         tg_zz_xxxxyyyz_0, tg_zz_xxxxyyyz_1, tg_zz_xxxxyyz_1, tg_zz_xxxxyyzz_0, \
                                         tg_zz_xxxxyyzz_1, tg_zz_xxxxyzz_1, tg_zz_xxxxyzzz_0, tg_zz_xxxxyzzz_1, \
                                         tg_zz_xxxxzzz_1, tg_zz_xxxxzzzz_0, tg_zz_xxxxzzzz_1, tg_zz_xxxyyyy_1, \
                                         tg_zz_xxxyyyyy_0, tg_zz_xxxyyyyy_1, tg_zz_xxxyyyyz_0, tg_zz_xxxyyyyz_1, \
                                         tg_zz_xxxyyyz_1, tg_zz_xxxyyyzz_0, tg_zz_xxxyyyzz_1, tg_zz_xxxyyzz_1, \
                                         tg_zz_xxxyyzzz_0, tg_zz_xxxyyzzz_1, tg_zz_xxxyzzz_1, tg_zz_xxxyzzzz_0, \
                                         tg_zz_xxxyzzzz_1, tg_zz_xxxzzzz_1, tg_zz_xxxzzzzz_0, tg_zz_xxxzzzzz_1, \
                                         tg_zz_xxyyyyy_1, tg_zz_xxyyyyyy_0, tg_zz_xxyyyyyy_1, tg_zz_xxyyyyyz_0, \
                                         tg_zz_xxyyyyyz_1, tg_zz_xxyyyyz_1, tg_zz_xxyyyyzz_0, tg_zz_xxyyyyzz_1, \
                                         tg_zz_xxyyyzz_1, tg_zz_xxyyyzzz_0, tg_zz_xxyyyzzz_1, tg_zz_xxyyzzz_1, \
                                         tg_zz_xxyyzzzz_0, tg_zz_xxyyzzzz_1, tg_zz_xxyzzzz_1, tg_zz_xxyzzzzz_0, \
                                         tg_zz_xxyzzzzz_1, tg_zz_xxzzzzz_1, tg_zz_xxzzzzzz_0, tg_zz_xxzzzzzz_1, \
                                         tg_zz_xyyyyyy_1, tg_zz_xyyyyyyy_0, tg_zz_xyyyyyyy_1, tg_zz_xyyyyyyz_0, \
                                         tg_zz_xyyyyyyz_1, tg_zz_xyyyyyz_1, tg_zz_xyyyyyzz_0, tg_zz_xyyyyyzz_1, \
                                         tg_zz_xyyyyzz_1, tg_zz_xyyyyzzz_0, tg_zz_xyyyyzzz_1, tg_zz_xyyyzzz_1, \
                                         tg_zz_xyyyzzzz_0, tg_zz_xyyyzzzz_1, tg_zz_xyyzzzz_1, tg_zz_xyyzzzzz_0, \
                                         tg_zz_xyyzzzzz_1, tg_zz_xyzzzzz_1, tg_zz_xyzzzzzz_0, tg_zz_xyzzzzzz_1, \
                                         tg_zz_xzzzzzz_1, tg_zz_xzzzzzzz_0, tg_zz_xzzzzzzz_1, tg_zz_yyyyyyy_1, \
                                         tg_zz_yyyyyyyy_0, tg_zz_yyyyyyyy_1, tg_zz_yyyyyyyz_0, tg_zz_yyyyyyyz_1, \
                                         tg_zz_yyyyyyz_1, tg_zz_yyyyyyzz_0, tg_zz_yyyyyyzz_1, tg_zz_yyyyyzz_1, \
                                         tg_zz_yyyyyzzz_0, tg_zz_yyyyyzzz_1, tg_zz_yyyyzzz_1, tg_zz_yyyyzzzz_0, \
                                         tg_zz_yyyyzzzz_1, tg_zz_yyyzzzz_1, tg_zz_yyyzzzzz_0, tg_zz_yyyzzzzz_1, \
                                         tg_zz_yyzzzzz_1, tg_zz_yyzzzzzz_0, tg_zz_yyzzzzzz_1, tg_zz_yzzzzzz_1, \
                                         tg_zz_yzzzzzzz_0, tg_zz_yzzzzzzz_1, tg_zz_zzzzzzz_1, tg_zz_zzzzzzzz_0, \
                                         tg_zz_zzzzzzzz_1, tg_zzz_xxxxxxxx_0, tg_zzz_xxxxxxxy_0, tg_zzz_xxxxxxxz_0, \
                                         tg_zzz_xxxxxxyy_0, tg_zzz_xxxxxxyz_0, tg_zzz_xxxxxxzz_0, tg_zzz_xxxxxyyy_0, \
                                         tg_zzz_xxxxxyyz_0, tg_zzz_xxxxxyzz_0, tg_zzz_xxxxxzzz_0, tg_zzz_xxxxyyyy_0, \
                                         tg_zzz_xxxxyyyz_0, tg_zzz_xxxxyyzz_0, tg_zzz_xxxxyzzz_0, tg_zzz_xxxxzzzz_0, \
                                         tg_zzz_xxxyyyyy_0, tg_zzz_xxxyyyyz_0, tg_zzz_xxxyyyzz_0, tg_zzz_xxxyyzzz_0, \
                                         tg_zzz_xxxyzzzz_0, tg_zzz_xxxzzzzz_0, tg_zzz_xxyyyyyy_0, tg_zzz_xxyyyyyz_0, \
                                         tg_zzz_xxyyyyzz_0, tg_zzz_xxyyyzzz_0, tg_zzz_xxyyzzzz_0, tg_zzz_xxyzzzzz_0, \
                                         tg_zzz_xxzzzzzz_0, tg_zzz_xyyyyyyy_0, tg_zzz_xyyyyyyz_0, tg_zzz_xyyyyyzz_0, \
                                         tg_zzz_xyyyyzzz_0, tg_zzz_xyyyzzzz_0, tg_zzz_xyyzzzzz_0, tg_zzz_xyzzzzzz_0, \
                                         tg_zzz_xzzzzzzz_0, tg_zzz_yyyyyyyy_0, tg_zzz_yyyyyyyz_0, tg_zzz_yyyyyyzz_0, \
                                         tg_zzz_yyyyyzzz_0, tg_zzz_yyyyzzzz_0, tg_zzz_yyyzzzzz_0, tg_zzz_yyzzzzzz_0, \
                                         tg_zzz_yzzzzzzz_0, tg_zzz_zzzzzzzz_0, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_yzz_xxxxxxxx_0[j] = pb_y * tg_zz_xxxxxxxx_0[j] + wp_y[j] * tg_zz_xxxxxxxx_1[j];

                    tg_yzz_xxxxxxxy_0[j] = pb_y * tg_zz_xxxxxxxy_0[j] + wp_y[j] * tg_zz_xxxxxxxy_1[j] + 0.5 * fl1_fxn * tg_zz_xxxxxxx_1[j];

                    tg_yzz_xxxxxxxz_0[j] = pb_y * tg_zz_xxxxxxxz_0[j] + wp_y[j] * tg_zz_xxxxxxxz_1[j];

                    tg_yzz_xxxxxxyy_0[j] = pb_y * tg_zz_xxxxxxyy_0[j] + wp_y[j] * tg_zz_xxxxxxyy_1[j] + fl1_fxn * tg_zz_xxxxxxy_1[j];

                    tg_yzz_xxxxxxyz_0[j] = pb_y * tg_zz_xxxxxxyz_0[j] + wp_y[j] * tg_zz_xxxxxxyz_1[j] + 0.5 * fl1_fxn * tg_zz_xxxxxxz_1[j];

                    tg_yzz_xxxxxxzz_0[j] = pb_y * tg_zz_xxxxxxzz_0[j] + wp_y[j] * tg_zz_xxxxxxzz_1[j];

                    tg_yzz_xxxxxyyy_0[j] = pb_y * tg_zz_xxxxxyyy_0[j] + wp_y[j] * tg_zz_xxxxxyyy_1[j] + 1.5 * fl1_fxn * tg_zz_xxxxxyy_1[j];

                    tg_yzz_xxxxxyyz_0[j] = pb_y * tg_zz_xxxxxyyz_0[j] + wp_y[j] * tg_zz_xxxxxyyz_1[j] + fl1_fxn * tg_zz_xxxxxyz_1[j];

                    tg_yzz_xxxxxyzz_0[j] = pb_y * tg_zz_xxxxxyzz_0[j] + wp_y[j] * tg_zz_xxxxxyzz_1[j] + 0.5 * fl1_fxn * tg_zz_xxxxxzz_1[j];

                    tg_yzz_xxxxxzzz_0[j] = pb_y * tg_zz_xxxxxzzz_0[j] + wp_y[j] * tg_zz_xxxxxzzz_1[j];

                    tg_yzz_xxxxyyyy_0[j] = pb_y * tg_zz_xxxxyyyy_0[j] + wp_y[j] * tg_zz_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_zz_xxxxyyy_1[j];

                    tg_yzz_xxxxyyyz_0[j] = pb_y * tg_zz_xxxxyyyz_0[j] + wp_y[j] * tg_zz_xxxxyyyz_1[j] + 1.5 * fl1_fxn * tg_zz_xxxxyyz_1[j];

                    tg_yzz_xxxxyyzz_0[j] = pb_y * tg_zz_xxxxyyzz_0[j] + wp_y[j] * tg_zz_xxxxyyzz_1[j] + fl1_fxn * tg_zz_xxxxyzz_1[j];

                    tg_yzz_xxxxyzzz_0[j] = pb_y * tg_zz_xxxxyzzz_0[j] + wp_y[j] * tg_zz_xxxxyzzz_1[j] + 0.5 * fl1_fxn * tg_zz_xxxxzzz_1[j];

                    tg_yzz_xxxxzzzz_0[j] = pb_y * tg_zz_xxxxzzzz_0[j] + wp_y[j] * tg_zz_xxxxzzzz_1[j];

                    tg_yzz_xxxyyyyy_0[j] = pb_y * tg_zz_xxxyyyyy_0[j] + wp_y[j] * tg_zz_xxxyyyyy_1[j] + 2.5 * fl1_fxn * tg_zz_xxxyyyy_1[j];

                    tg_yzz_xxxyyyyz_0[j] = pb_y * tg_zz_xxxyyyyz_0[j] + wp_y[j] * tg_zz_xxxyyyyz_1[j] + 2.0 * fl1_fxn * tg_zz_xxxyyyz_1[j];

                    tg_yzz_xxxyyyzz_0[j] = pb_y * tg_zz_xxxyyyzz_0[j] + wp_y[j] * tg_zz_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_zz_xxxyyzz_1[j];

                    tg_yzz_xxxyyzzz_0[j] = pb_y * tg_zz_xxxyyzzz_0[j] + wp_y[j] * tg_zz_xxxyyzzz_1[j] + fl1_fxn * tg_zz_xxxyzzz_1[j];

                    tg_yzz_xxxyzzzz_0[j] = pb_y * tg_zz_xxxyzzzz_0[j] + wp_y[j] * tg_zz_xxxyzzzz_1[j] + 0.5 * fl1_fxn * tg_zz_xxxzzzz_1[j];

                    tg_yzz_xxxzzzzz_0[j] = pb_y * tg_zz_xxxzzzzz_0[j] + wp_y[j] * tg_zz_xxxzzzzz_1[j];

                    tg_yzz_xxyyyyyy_0[j] = pb_y * tg_zz_xxyyyyyy_0[j] + wp_y[j] * tg_zz_xxyyyyyy_1[j] + 3.0 * fl1_fxn * tg_zz_xxyyyyy_1[j];

                    tg_yzz_xxyyyyyz_0[j] = pb_y * tg_zz_xxyyyyyz_0[j] + wp_y[j] * tg_zz_xxyyyyyz_1[j] + 2.5 * fl1_fxn * tg_zz_xxyyyyz_1[j];

                    tg_yzz_xxyyyyzz_0[j] = pb_y * tg_zz_xxyyyyzz_0[j] + wp_y[j] * tg_zz_xxyyyyzz_1[j] + 2.0 * fl1_fxn * tg_zz_xxyyyzz_1[j];

                    tg_yzz_xxyyyzzz_0[j] = pb_y * tg_zz_xxyyyzzz_0[j] + wp_y[j] * tg_zz_xxyyyzzz_1[j] + 1.5 * fl1_fxn * tg_zz_xxyyzzz_1[j];

                    tg_yzz_xxyyzzzz_0[j] = pb_y * tg_zz_xxyyzzzz_0[j] + wp_y[j] * tg_zz_xxyyzzzz_1[j] + fl1_fxn * tg_zz_xxyzzzz_1[j];

                    tg_yzz_xxyzzzzz_0[j] = pb_y * tg_zz_xxyzzzzz_0[j] + wp_y[j] * tg_zz_xxyzzzzz_1[j] + 0.5 * fl1_fxn * tg_zz_xxzzzzz_1[j];

                    tg_yzz_xxzzzzzz_0[j] = pb_y * tg_zz_xxzzzzzz_0[j] + wp_y[j] * tg_zz_xxzzzzzz_1[j];

                    tg_yzz_xyyyyyyy_0[j] = pb_y * tg_zz_xyyyyyyy_0[j] + wp_y[j] * tg_zz_xyyyyyyy_1[j] + 3.5 * fl1_fxn * tg_zz_xyyyyyy_1[j];

                    tg_yzz_xyyyyyyz_0[j] = pb_y * tg_zz_xyyyyyyz_0[j] + wp_y[j] * tg_zz_xyyyyyyz_1[j] + 3.0 * fl1_fxn * tg_zz_xyyyyyz_1[j];

                    tg_yzz_xyyyyyzz_0[j] = pb_y * tg_zz_xyyyyyzz_0[j] + wp_y[j] * tg_zz_xyyyyyzz_1[j] + 2.5 * fl1_fxn * tg_zz_xyyyyzz_1[j];

                    tg_yzz_xyyyyzzz_0[j] = pb_y * tg_zz_xyyyyzzz_0[j] + wp_y[j] * tg_zz_xyyyyzzz_1[j] + 2.0 * fl1_fxn * tg_zz_xyyyzzz_1[j];

                    tg_yzz_xyyyzzzz_0[j] = pb_y * tg_zz_xyyyzzzz_0[j] + wp_y[j] * tg_zz_xyyyzzzz_1[j] + 1.5 * fl1_fxn * tg_zz_xyyzzzz_1[j];

                    tg_yzz_xyyzzzzz_0[j] = pb_y * tg_zz_xyyzzzzz_0[j] + wp_y[j] * tg_zz_xyyzzzzz_1[j] + fl1_fxn * tg_zz_xyzzzzz_1[j];

                    tg_yzz_xyzzzzzz_0[j] = pb_y * tg_zz_xyzzzzzz_0[j] + wp_y[j] * tg_zz_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_zz_xzzzzzz_1[j];

                    tg_yzz_xzzzzzzz_0[j] = pb_y * tg_zz_xzzzzzzz_0[j] + wp_y[j] * tg_zz_xzzzzzzz_1[j];

                    tg_yzz_yyyyyyyy_0[j] = pb_y * tg_zz_yyyyyyyy_0[j] + wp_y[j] * tg_zz_yyyyyyyy_1[j] + 4.0 * fl1_fxn * tg_zz_yyyyyyy_1[j];

                    tg_yzz_yyyyyyyz_0[j] = pb_y * tg_zz_yyyyyyyz_0[j] + wp_y[j] * tg_zz_yyyyyyyz_1[j] + 3.5 * fl1_fxn * tg_zz_yyyyyyz_1[j];

                    tg_yzz_yyyyyyzz_0[j] = pb_y * tg_zz_yyyyyyzz_0[j] + wp_y[j] * tg_zz_yyyyyyzz_1[j] + 3.0 * fl1_fxn * tg_zz_yyyyyzz_1[j];

                    tg_yzz_yyyyyzzz_0[j] = pb_y * tg_zz_yyyyyzzz_0[j] + wp_y[j] * tg_zz_yyyyyzzz_1[j] + 2.5 * fl1_fxn * tg_zz_yyyyzzz_1[j];

                    tg_yzz_yyyyzzzz_0[j] = pb_y * tg_zz_yyyyzzzz_0[j] + wp_y[j] * tg_zz_yyyyzzzz_1[j] + 2.0 * fl1_fxn * tg_zz_yyyzzzz_1[j];

                    tg_yzz_yyyzzzzz_0[j] = pb_y * tg_zz_yyyzzzzz_0[j] + wp_y[j] * tg_zz_yyyzzzzz_1[j] + 1.5 * fl1_fxn * tg_zz_yyzzzzz_1[j];

                    tg_yzz_yyzzzzzz_0[j] = pb_y * tg_zz_yyzzzzzz_0[j] + wp_y[j] * tg_zz_yyzzzzzz_1[j] + fl1_fxn * tg_zz_yzzzzzz_1[j];

                    tg_yzz_yzzzzzzz_0[j] = pb_y * tg_zz_yzzzzzzz_0[j] + wp_y[j] * tg_zz_yzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_zz_zzzzzzz_1[j];

                    tg_yzz_zzzzzzzz_0[j] = pb_y * tg_zz_zzzzzzzz_0[j] + wp_y[j] * tg_zz_zzzzzzzz_1[j];

                    tg_zzz_xxxxxxxx_0[j] = pb_z * tg_zz_xxxxxxxx_0[j] + wp_z[j] * tg_zz_xxxxxxxx_1[j] + fl1_fx * tg_z_xxxxxxxx_0[j] - fl1_fx * fl1_fza * tg_z_xxxxxxxx_1[j];

                    tg_zzz_xxxxxxxy_0[j] = pb_z * tg_zz_xxxxxxxy_0[j] + wp_z[j] * tg_zz_xxxxxxxy_1[j] + fl1_fx * tg_z_xxxxxxxy_0[j] - fl1_fx * fl1_fza * tg_z_xxxxxxxy_1[j];

                    tg_zzz_xxxxxxxz_0[j] = pb_z * tg_zz_xxxxxxxz_0[j] + wp_z[j] * tg_zz_xxxxxxxz_1[j] + fl1_fx * tg_z_xxxxxxxz_0[j] - fl1_fx * fl1_fza * tg_z_xxxxxxxz_1[j] + 0.5 * fl1_fxn * tg_zz_xxxxxxx_1[j];

                    tg_zzz_xxxxxxyy_0[j] = pb_z * tg_zz_xxxxxxyy_0[j] + wp_z[j] * tg_zz_xxxxxxyy_1[j] + fl1_fx * tg_z_xxxxxxyy_0[j] - fl1_fx * fl1_fza * tg_z_xxxxxxyy_1[j];

                    tg_zzz_xxxxxxyz_0[j] = pb_z * tg_zz_xxxxxxyz_0[j] + wp_z[j] * tg_zz_xxxxxxyz_1[j] + fl1_fx * tg_z_xxxxxxyz_0[j] - fl1_fx * fl1_fza * tg_z_xxxxxxyz_1[j] + 0.5 * fl1_fxn * tg_zz_xxxxxxy_1[j];

                    tg_zzz_xxxxxxzz_0[j] = pb_z * tg_zz_xxxxxxzz_0[j] + wp_z[j] * tg_zz_xxxxxxzz_1[j] + fl1_fx * tg_z_xxxxxxzz_0[j] - fl1_fx * fl1_fza * tg_z_xxxxxxzz_1[j] + fl1_fxn * tg_zz_xxxxxxz_1[j];

                    tg_zzz_xxxxxyyy_0[j] = pb_z * tg_zz_xxxxxyyy_0[j] + wp_z[j] * tg_zz_xxxxxyyy_1[j] + fl1_fx * tg_z_xxxxxyyy_0[j] - fl1_fx * fl1_fza * tg_z_xxxxxyyy_1[j];

                    tg_zzz_xxxxxyyz_0[j] = pb_z * tg_zz_xxxxxyyz_0[j] + wp_z[j] * tg_zz_xxxxxyyz_1[j] + fl1_fx * tg_z_xxxxxyyz_0[j] - fl1_fx * fl1_fza * tg_z_xxxxxyyz_1[j] + 0.5 * fl1_fxn * tg_zz_xxxxxyy_1[j];

                    tg_zzz_xxxxxyzz_0[j] = pb_z * tg_zz_xxxxxyzz_0[j] + wp_z[j] * tg_zz_xxxxxyzz_1[j] + fl1_fx * tg_z_xxxxxyzz_0[j] - fl1_fx * fl1_fza * tg_z_xxxxxyzz_1[j] + fl1_fxn * tg_zz_xxxxxyz_1[j];

                    tg_zzz_xxxxxzzz_0[j] = pb_z * tg_zz_xxxxxzzz_0[j] + wp_z[j] * tg_zz_xxxxxzzz_1[j] + fl1_fx * tg_z_xxxxxzzz_0[j] - fl1_fx * fl1_fza * tg_z_xxxxxzzz_1[j] + 1.5 * fl1_fxn * tg_zz_xxxxxzz_1[j];

                    tg_zzz_xxxxyyyy_0[j] = pb_z * tg_zz_xxxxyyyy_0[j] + wp_z[j] * tg_zz_xxxxyyyy_1[j] + fl1_fx * tg_z_xxxxyyyy_0[j] - fl1_fx * fl1_fza * tg_z_xxxxyyyy_1[j];

                    tg_zzz_xxxxyyyz_0[j] = pb_z * tg_zz_xxxxyyyz_0[j] + wp_z[j] * tg_zz_xxxxyyyz_1[j] + fl1_fx * tg_z_xxxxyyyz_0[j] - fl1_fx * fl1_fza * tg_z_xxxxyyyz_1[j] + 0.5 * fl1_fxn * tg_zz_xxxxyyy_1[j];

                    tg_zzz_xxxxyyzz_0[j] = pb_z * tg_zz_xxxxyyzz_0[j] + wp_z[j] * tg_zz_xxxxyyzz_1[j] + fl1_fx * tg_z_xxxxyyzz_0[j] - fl1_fx * fl1_fza * tg_z_xxxxyyzz_1[j] + fl1_fxn * tg_zz_xxxxyyz_1[j];

                    tg_zzz_xxxxyzzz_0[j] = pb_z * tg_zz_xxxxyzzz_0[j] + wp_z[j] * tg_zz_xxxxyzzz_1[j] + fl1_fx * tg_z_xxxxyzzz_0[j] - fl1_fx * fl1_fza * tg_z_xxxxyzzz_1[j] + 1.5 * fl1_fxn * tg_zz_xxxxyzz_1[j];

                    tg_zzz_xxxxzzzz_0[j] = pb_z * tg_zz_xxxxzzzz_0[j] + wp_z[j] * tg_zz_xxxxzzzz_1[j] + fl1_fx * tg_z_xxxxzzzz_0[j] - fl1_fx * fl1_fza * tg_z_xxxxzzzz_1[j] + 2.0 * fl1_fxn * tg_zz_xxxxzzz_1[j];

                    tg_zzz_xxxyyyyy_0[j] = pb_z * tg_zz_xxxyyyyy_0[j] + wp_z[j] * tg_zz_xxxyyyyy_1[j] + fl1_fx * tg_z_xxxyyyyy_0[j] - fl1_fx * fl1_fza * tg_z_xxxyyyyy_1[j];

                    tg_zzz_xxxyyyyz_0[j] = pb_z * tg_zz_xxxyyyyz_0[j] + wp_z[j] * tg_zz_xxxyyyyz_1[j] + fl1_fx * tg_z_xxxyyyyz_0[j] - fl1_fx * fl1_fza * tg_z_xxxyyyyz_1[j] + 0.5 * fl1_fxn * tg_zz_xxxyyyy_1[j];

                    tg_zzz_xxxyyyzz_0[j] = pb_z * tg_zz_xxxyyyzz_0[j] + wp_z[j] * tg_zz_xxxyyyzz_1[j] + fl1_fx * tg_z_xxxyyyzz_0[j] - fl1_fx * fl1_fza * tg_z_xxxyyyzz_1[j] + fl1_fxn * tg_zz_xxxyyyz_1[j];

                    tg_zzz_xxxyyzzz_0[j] = pb_z * tg_zz_xxxyyzzz_0[j] + wp_z[j] * tg_zz_xxxyyzzz_1[j] + fl1_fx * tg_z_xxxyyzzz_0[j] - fl1_fx * fl1_fza * tg_z_xxxyyzzz_1[j] + 1.5 * fl1_fxn * tg_zz_xxxyyzz_1[j];

                    tg_zzz_xxxyzzzz_0[j] = pb_z * tg_zz_xxxyzzzz_0[j] + wp_z[j] * tg_zz_xxxyzzzz_1[j] + fl1_fx * tg_z_xxxyzzzz_0[j] - fl1_fx * fl1_fza * tg_z_xxxyzzzz_1[j] + 2.0 * fl1_fxn * tg_zz_xxxyzzz_1[j];

                    tg_zzz_xxxzzzzz_0[j] = pb_z * tg_zz_xxxzzzzz_0[j] + wp_z[j] * tg_zz_xxxzzzzz_1[j] + fl1_fx * tg_z_xxxzzzzz_0[j] - fl1_fx * fl1_fza * tg_z_xxxzzzzz_1[j] + 2.5 * fl1_fxn * tg_zz_xxxzzzz_1[j];

                    tg_zzz_xxyyyyyy_0[j] = pb_z * tg_zz_xxyyyyyy_0[j] + wp_z[j] * tg_zz_xxyyyyyy_1[j] + fl1_fx * tg_z_xxyyyyyy_0[j] - fl1_fx * fl1_fza * tg_z_xxyyyyyy_1[j];

                    tg_zzz_xxyyyyyz_0[j] = pb_z * tg_zz_xxyyyyyz_0[j] + wp_z[j] * tg_zz_xxyyyyyz_1[j] + fl1_fx * tg_z_xxyyyyyz_0[j] - fl1_fx * fl1_fza * tg_z_xxyyyyyz_1[j] + 0.5 * fl1_fxn * tg_zz_xxyyyyy_1[j];

                    tg_zzz_xxyyyyzz_0[j] = pb_z * tg_zz_xxyyyyzz_0[j] + wp_z[j] * tg_zz_xxyyyyzz_1[j] + fl1_fx * tg_z_xxyyyyzz_0[j] - fl1_fx * fl1_fza * tg_z_xxyyyyzz_1[j] + fl1_fxn * tg_zz_xxyyyyz_1[j];

                    tg_zzz_xxyyyzzz_0[j] = pb_z * tg_zz_xxyyyzzz_0[j] + wp_z[j] * tg_zz_xxyyyzzz_1[j] + fl1_fx * tg_z_xxyyyzzz_0[j] - fl1_fx * fl1_fza * tg_z_xxyyyzzz_1[j] + 1.5 * fl1_fxn * tg_zz_xxyyyzz_1[j];

                    tg_zzz_xxyyzzzz_0[j] = pb_z * tg_zz_xxyyzzzz_0[j] + wp_z[j] * tg_zz_xxyyzzzz_1[j] + fl1_fx * tg_z_xxyyzzzz_0[j] - fl1_fx * fl1_fza * tg_z_xxyyzzzz_1[j] + 2.0 * fl1_fxn * tg_zz_xxyyzzz_1[j];

                    tg_zzz_xxyzzzzz_0[j] = pb_z * tg_zz_xxyzzzzz_0[j] + wp_z[j] * tg_zz_xxyzzzzz_1[j] + fl1_fx * tg_z_xxyzzzzz_0[j] - fl1_fx * fl1_fza * tg_z_xxyzzzzz_1[j] + 2.5 * fl1_fxn * tg_zz_xxyzzzz_1[j];

                    tg_zzz_xxzzzzzz_0[j] = pb_z * tg_zz_xxzzzzzz_0[j] + wp_z[j] * tg_zz_xxzzzzzz_1[j] + fl1_fx * tg_z_xxzzzzzz_0[j] - fl1_fx * fl1_fza * tg_z_xxzzzzzz_1[j] + 3.0 * fl1_fxn * tg_zz_xxzzzzz_1[j];

                    tg_zzz_xyyyyyyy_0[j] = pb_z * tg_zz_xyyyyyyy_0[j] + wp_z[j] * tg_zz_xyyyyyyy_1[j] + fl1_fx * tg_z_xyyyyyyy_0[j] - fl1_fx * fl1_fza * tg_z_xyyyyyyy_1[j];

                    tg_zzz_xyyyyyyz_0[j] = pb_z * tg_zz_xyyyyyyz_0[j] + wp_z[j] * tg_zz_xyyyyyyz_1[j] + fl1_fx * tg_z_xyyyyyyz_0[j] - fl1_fx * fl1_fza * tg_z_xyyyyyyz_1[j] + 0.5 * fl1_fxn * tg_zz_xyyyyyy_1[j];

                    tg_zzz_xyyyyyzz_0[j] = pb_z * tg_zz_xyyyyyzz_0[j] + wp_z[j] * tg_zz_xyyyyyzz_1[j] + fl1_fx * tg_z_xyyyyyzz_0[j] - fl1_fx * fl1_fza * tg_z_xyyyyyzz_1[j] + fl1_fxn * tg_zz_xyyyyyz_1[j];

                    tg_zzz_xyyyyzzz_0[j] = pb_z * tg_zz_xyyyyzzz_0[j] + wp_z[j] * tg_zz_xyyyyzzz_1[j] + fl1_fx * tg_z_xyyyyzzz_0[j] - fl1_fx * fl1_fza * tg_z_xyyyyzzz_1[j] + 1.5 * fl1_fxn * tg_zz_xyyyyzz_1[j];

                    tg_zzz_xyyyzzzz_0[j] = pb_z * tg_zz_xyyyzzzz_0[j] + wp_z[j] * tg_zz_xyyyzzzz_1[j] + fl1_fx * tg_z_xyyyzzzz_0[j] - fl1_fx * fl1_fza * tg_z_xyyyzzzz_1[j] + 2.0 * fl1_fxn * tg_zz_xyyyzzz_1[j];

                    tg_zzz_xyyzzzzz_0[j] = pb_z * tg_zz_xyyzzzzz_0[j] + wp_z[j] * tg_zz_xyyzzzzz_1[j] + fl1_fx * tg_z_xyyzzzzz_0[j] - fl1_fx * fl1_fza * tg_z_xyyzzzzz_1[j] + 2.5 * fl1_fxn * tg_zz_xyyzzzz_1[j];

                    tg_zzz_xyzzzzzz_0[j] = pb_z * tg_zz_xyzzzzzz_0[j] + wp_z[j] * tg_zz_xyzzzzzz_1[j] + fl1_fx * tg_z_xyzzzzzz_0[j] - fl1_fx * fl1_fza * tg_z_xyzzzzzz_1[j] + 3.0 * fl1_fxn * tg_zz_xyzzzzz_1[j];

                    tg_zzz_xzzzzzzz_0[j] = pb_z * tg_zz_xzzzzzzz_0[j] + wp_z[j] * tg_zz_xzzzzzzz_1[j] + fl1_fx * tg_z_xzzzzzzz_0[j] - fl1_fx * fl1_fza * tg_z_xzzzzzzz_1[j] + 3.5 * fl1_fxn * tg_zz_xzzzzzz_1[j];

                    tg_zzz_yyyyyyyy_0[j] = pb_z * tg_zz_yyyyyyyy_0[j] + wp_z[j] * tg_zz_yyyyyyyy_1[j] + fl1_fx * tg_z_yyyyyyyy_0[j] - fl1_fx * fl1_fza * tg_z_yyyyyyyy_1[j];

                    tg_zzz_yyyyyyyz_0[j] = pb_z * tg_zz_yyyyyyyz_0[j] + wp_z[j] * tg_zz_yyyyyyyz_1[j] + fl1_fx * tg_z_yyyyyyyz_0[j] - fl1_fx * fl1_fza * tg_z_yyyyyyyz_1[j] + 0.5 * fl1_fxn * tg_zz_yyyyyyy_1[j];

                    tg_zzz_yyyyyyzz_0[j] = pb_z * tg_zz_yyyyyyzz_0[j] + wp_z[j] * tg_zz_yyyyyyzz_1[j] + fl1_fx * tg_z_yyyyyyzz_0[j] - fl1_fx * fl1_fza * tg_z_yyyyyyzz_1[j] + fl1_fxn * tg_zz_yyyyyyz_1[j];

                    tg_zzz_yyyyyzzz_0[j] = pb_z * tg_zz_yyyyyzzz_0[j] + wp_z[j] * tg_zz_yyyyyzzz_1[j] + fl1_fx * tg_z_yyyyyzzz_0[j] - fl1_fx * fl1_fza * tg_z_yyyyyzzz_1[j] + 1.5 * fl1_fxn * tg_zz_yyyyyzz_1[j];

                    tg_zzz_yyyyzzzz_0[j] = pb_z * tg_zz_yyyyzzzz_0[j] + wp_z[j] * tg_zz_yyyyzzzz_1[j] + fl1_fx * tg_z_yyyyzzzz_0[j] - fl1_fx * fl1_fza * tg_z_yyyyzzzz_1[j] + 2.0 * fl1_fxn * tg_zz_yyyyzzz_1[j];

                    tg_zzz_yyyzzzzz_0[j] = pb_z * tg_zz_yyyzzzzz_0[j] + wp_z[j] * tg_zz_yyyzzzzz_1[j] + fl1_fx * tg_z_yyyzzzzz_0[j] - fl1_fx * fl1_fza * tg_z_yyyzzzzz_1[j] + 2.5 * fl1_fxn * tg_zz_yyyzzzz_1[j];

                    tg_zzz_yyzzzzzz_0[j] = pb_z * tg_zz_yyzzzzzz_0[j] + wp_z[j] * tg_zz_yyzzzzzz_1[j] + fl1_fx * tg_z_yyzzzzzz_0[j] - fl1_fx * fl1_fza * tg_z_yyzzzzzz_1[j] + 3.0 * fl1_fxn * tg_zz_yyzzzzz_1[j];

                    tg_zzz_yzzzzzzz_0[j] = pb_z * tg_zz_yzzzzzzz_0[j] + wp_z[j] * tg_zz_yzzzzzzz_1[j] + fl1_fx * tg_z_yzzzzzzz_0[j] - fl1_fx * fl1_fza * tg_z_yzzzzzzz_1[j] + 3.5 * fl1_fxn * tg_zz_yzzzzzz_1[j];

                    tg_zzz_zzzzzzzz_0[j] = pb_z * tg_zz_zzzzzzzz_0[j] + wp_z[j] * tg_zz_zzzzzzzz_1[j] + fl1_fx * tg_z_zzzzzzzz_0[j] - fl1_fx * fl1_fza * tg_z_zzzzzzzz_1[j] + 4.0 * fl1_fxn * tg_zz_zzzzzzz_1[j];
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

