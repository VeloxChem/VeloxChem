//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectronRepulsionRecFuncForGL.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSGSL(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSGSL_0_49(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSGSL_49_98(primBuffer,
                                                       recursionMap,
                                                       osFactors,
                                                       wpDistances, 
                                                       braGtoPairsBlock,
                                                       ketGtoPairsBlock,
                                                       nKetPrimPairs,
                                                       iContrPair); 

        erirecfunc::compElectronRepulsionForSGSL_98_147(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 

        erirecfunc::compElectronRepulsionForSGSL_147_195(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSGSL_195_243(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSGSL_243_291(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSGSL_291_339(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSGSL_339_387(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSGSL_387_435(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSGSL_435_483(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSGSL_483_531(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSGSL_531_579(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSGSL_579_627(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSGSL_627_675(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 
    }

    void
    compElectronRepulsionForSGSL_0_49(      CMemBlock2D<double>& primBuffer,
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
                                             {4, -1, -1, -1},
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_2_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_2_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

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

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                // set up pointers to auxilary integrals

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

                auto tg_xxx_xxxxxxxx_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx); 

                auto tg_xxx_xxxxxxxy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 1); 

                auto tg_xxx_xxxxxxxz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 2); 

                auto tg_xxx_xxxxxxyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 3); 

                auto tg_xxx_xxxxxxyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 4); 

                auto tg_xxx_xxxxxxzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 5); 

                auto tg_xxx_xxxxxyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 6); 

                auto tg_xxx_xxxxxyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 7); 

                auto tg_xxx_xxxxxyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 8); 

                auto tg_xxx_xxxxxzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 9); 

                auto tg_xxx_xxxxyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 10); 

                auto tg_xxx_xxxxyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 11); 

                auto tg_xxx_xxxxyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 12); 

                auto tg_xxx_xxxxyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 13); 

                auto tg_xxx_xxxxzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 14); 

                auto tg_xxx_xxxyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 15); 

                auto tg_xxx_xxxyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 16); 

                auto tg_xxx_xxxyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 17); 

                auto tg_xxx_xxxyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 18); 

                auto tg_xxx_xxxyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 19); 

                auto tg_xxx_xxxzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 20); 

                auto tg_xxx_xxyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 21); 

                auto tg_xxx_xxyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 22); 

                auto tg_xxx_xxyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 23); 

                auto tg_xxx_xxyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 24); 

                auto tg_xxx_xxyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 25); 

                auto tg_xxx_xxyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 26); 

                auto tg_xxx_xxzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 27); 

                auto tg_xxx_xyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 28); 

                auto tg_xxx_xyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 29); 

                auto tg_xxx_xyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 30); 

                auto tg_xxx_xyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 31); 

                auto tg_xxx_xyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 32); 

                auto tg_xxx_xyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 33); 

                auto tg_xxx_xyzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 34); 

                auto tg_xxx_xzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 35); 

                auto tg_xxx_yyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 36); 

                auto tg_xxx_yyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 37); 

                auto tg_xxx_yyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 38); 

                auto tg_xxx_yyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 39); 

                auto tg_xxx_yyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 40); 

                auto tg_xxx_yyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 41); 

                auto tg_xxx_yyzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 42); 

                auto tg_xxx_yzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 43); 

                auto tg_xxx_zzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 44); 

                auto tg_xxy_xxxxxxxx_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 45); 

                auto tg_xxy_xxxxxxxy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 46); 

                auto tg_xxy_xxxxxxxz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 47); 

                auto tg_xxy_xxxxxxyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 48); 

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

                // set up pointers to integrals

                auto tg_xxxx_xxxxxxxx_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx); 

                auto tg_xxxx_xxxxxxxy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 1); 

                auto tg_xxxx_xxxxxxxz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 2); 

                auto tg_xxxx_xxxxxxyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 3); 

                auto tg_xxxx_xxxxxxyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 4); 

                auto tg_xxxx_xxxxxxzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 5); 

                auto tg_xxxx_xxxxxyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 6); 

                auto tg_xxxx_xxxxxyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 7); 

                auto tg_xxxx_xxxxxyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 8); 

                auto tg_xxxx_xxxxxzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 9); 

                auto tg_xxxx_xxxxyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 10); 

                auto tg_xxxx_xxxxyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 11); 

                auto tg_xxxx_xxxxyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 12); 

                auto tg_xxxx_xxxxyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 13); 

                auto tg_xxxx_xxxxzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 14); 

                auto tg_xxxx_xxxyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 15); 

                auto tg_xxxx_xxxyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 16); 

                auto tg_xxxx_xxxyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 17); 

                auto tg_xxxx_xxxyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 18); 

                auto tg_xxxx_xxxyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 19); 

                auto tg_xxxx_xxxzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 20); 

                auto tg_xxxx_xxyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 21); 

                auto tg_xxxx_xxyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 22); 

                auto tg_xxxx_xxyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 23); 

                auto tg_xxxx_xxyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 24); 

                auto tg_xxxx_xxyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 25); 

                auto tg_xxxx_xxyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 26); 

                auto tg_xxxx_xxzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 27); 

                auto tg_xxxx_xyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 28); 

                auto tg_xxxx_xyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 29); 

                auto tg_xxxx_xyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 30); 

                auto tg_xxxx_xyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 31); 

                auto tg_xxxx_xyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 32); 

                auto tg_xxxx_xyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 33); 

                auto tg_xxxx_xyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 34); 

                auto tg_xxxx_xzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 35); 

                auto tg_xxxx_yyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 36); 

                auto tg_xxxx_yyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 37); 

                auto tg_xxxx_yyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 38); 

                auto tg_xxxx_yyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 39); 

                auto tg_xxxx_yyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 40); 

                auto tg_xxxx_yyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 41); 

                auto tg_xxxx_yyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 42); 

                auto tg_xxxx_yzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 43); 

                auto tg_xxxx_zzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 44); 

                auto tg_xxxy_xxxxxxxx_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 45); 

                auto tg_xxxy_xxxxxxxy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 46); 

                auto tg_xxxy_xxxxxxxz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 47); 

                auto tg_xxxy_xxxxxxyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 48); 

                // Batch of Integrals (0,49)

                #pragma omp simd aligned(fxn, fza, tg_xx_xxxxxxxx_0, tg_xx_xxxxxxxx_1, tg_xx_xxxxxxxy_0, \
                                         tg_xx_xxxxxxxy_1, tg_xx_xxxxxxxz_0, tg_xx_xxxxxxxz_1, tg_xx_xxxxxxyy_0, \
                                         tg_xx_xxxxxxyy_1, tg_xx_xxxxxxyz_0, tg_xx_xxxxxxyz_1, tg_xx_xxxxxxzz_0, \
                                         tg_xx_xxxxxxzz_1, tg_xx_xxxxxyyy_0, tg_xx_xxxxxyyy_1, tg_xx_xxxxxyyz_0, \
                                         tg_xx_xxxxxyyz_1, tg_xx_xxxxxyzz_0, tg_xx_xxxxxyzz_1, tg_xx_xxxxxzzz_0, \
                                         tg_xx_xxxxxzzz_1, tg_xx_xxxxyyyy_0, tg_xx_xxxxyyyy_1, tg_xx_xxxxyyyz_0, \
                                         tg_xx_xxxxyyyz_1, tg_xx_xxxxyyzz_0, tg_xx_xxxxyyzz_1, tg_xx_xxxxyzzz_0, \
                                         tg_xx_xxxxyzzz_1, tg_xx_xxxxzzzz_0, tg_xx_xxxxzzzz_1, tg_xx_xxxyyyyy_0, \
                                         tg_xx_xxxyyyyy_1, tg_xx_xxxyyyyz_0, tg_xx_xxxyyyyz_1, tg_xx_xxxyyyzz_0, \
                                         tg_xx_xxxyyyzz_1, tg_xx_xxxyyzzz_0, tg_xx_xxxyyzzz_1, tg_xx_xxxyzzzz_0, \
                                         tg_xx_xxxyzzzz_1, tg_xx_xxxzzzzz_0, tg_xx_xxxzzzzz_1, tg_xx_xxyyyyyy_0, \
                                         tg_xx_xxyyyyyy_1, tg_xx_xxyyyyyz_0, tg_xx_xxyyyyyz_1, tg_xx_xxyyyyzz_0, \
                                         tg_xx_xxyyyyzz_1, tg_xx_xxyyyzzz_0, tg_xx_xxyyyzzz_1, tg_xx_xxyyzzzz_0, \
                                         tg_xx_xxyyzzzz_1, tg_xx_xxyzzzzz_0, tg_xx_xxyzzzzz_1, tg_xx_xxzzzzzz_0, \
                                         tg_xx_xxzzzzzz_1, tg_xx_xyyyyyyy_0, tg_xx_xyyyyyyy_1, tg_xx_xyyyyyyz_0, \
                                         tg_xx_xyyyyyyz_1, tg_xx_xyyyyyzz_0, tg_xx_xyyyyyzz_1, tg_xx_xyyyyzzz_0, \
                                         tg_xx_xyyyyzzz_1, tg_xx_xyyyzzzz_0, tg_xx_xyyyzzzz_1, tg_xx_xyyzzzzz_0, \
                                         tg_xx_xyyzzzzz_1, tg_xx_xyzzzzzz_0, tg_xx_xyzzzzzz_1, tg_xx_xzzzzzzz_0, \
                                         tg_xx_xzzzzzzz_1, tg_xx_yyyyyyyy_0, tg_xx_yyyyyyyy_1, tg_xx_yyyyyyyz_0, \
                                         tg_xx_yyyyyyyz_1, tg_xx_yyyyyyzz_0, tg_xx_yyyyyyzz_1, tg_xx_yyyyyzzz_0, \
                                         tg_xx_yyyyyzzz_1, tg_xx_yyyyzzzz_0, tg_xx_yyyyzzzz_1, tg_xx_yyyzzzzz_0, \
                                         tg_xx_yyyzzzzz_1, tg_xx_yyzzzzzz_0, tg_xx_yyzzzzzz_1, tg_xx_yzzzzzzz_0, \
                                         tg_xx_yzzzzzzz_1, tg_xx_zzzzzzzz_0, tg_xx_zzzzzzzz_1, tg_xxx_xxxxxxx_1, \
                                         tg_xxx_xxxxxxxx_0, tg_xxx_xxxxxxxx_1, tg_xxx_xxxxxxxy_0, tg_xxx_xxxxxxxy_1, \
                                         tg_xxx_xxxxxxxz_0, tg_xxx_xxxxxxxz_1, tg_xxx_xxxxxxy_1, tg_xxx_xxxxxxyy_0, \
                                         tg_xxx_xxxxxxyy_1, tg_xxx_xxxxxxyz_0, tg_xxx_xxxxxxyz_1, tg_xxx_xxxxxxz_1, \
                                         tg_xxx_xxxxxxzz_0, tg_xxx_xxxxxxzz_1, tg_xxx_xxxxxyy_1, tg_xxx_xxxxxyyy_0, \
                                         tg_xxx_xxxxxyyy_1, tg_xxx_xxxxxyyz_0, tg_xxx_xxxxxyyz_1, tg_xxx_xxxxxyz_1, \
                                         tg_xxx_xxxxxyzz_0, tg_xxx_xxxxxyzz_1, tg_xxx_xxxxxzz_1, tg_xxx_xxxxxzzz_0, \
                                         tg_xxx_xxxxxzzz_1, tg_xxx_xxxxyyy_1, tg_xxx_xxxxyyyy_0, tg_xxx_xxxxyyyy_1, \
                                         tg_xxx_xxxxyyyz_0, tg_xxx_xxxxyyyz_1, tg_xxx_xxxxyyz_1, tg_xxx_xxxxyyzz_0, \
                                         tg_xxx_xxxxyyzz_1, tg_xxx_xxxxyzz_1, tg_xxx_xxxxyzzz_0, tg_xxx_xxxxyzzz_1, \
                                         tg_xxx_xxxxzzz_1, tg_xxx_xxxxzzzz_0, tg_xxx_xxxxzzzz_1, tg_xxx_xxxyyyy_1, \
                                         tg_xxx_xxxyyyyy_0, tg_xxx_xxxyyyyy_1, tg_xxx_xxxyyyyz_0, tg_xxx_xxxyyyyz_1, \
                                         tg_xxx_xxxyyyz_1, tg_xxx_xxxyyyzz_0, tg_xxx_xxxyyyzz_1, tg_xxx_xxxyyzz_1, \
                                         tg_xxx_xxxyyzzz_0, tg_xxx_xxxyyzzz_1, tg_xxx_xxxyzzz_1, tg_xxx_xxxyzzzz_0, \
                                         tg_xxx_xxxyzzzz_1, tg_xxx_xxxzzzz_1, tg_xxx_xxxzzzzz_0, tg_xxx_xxxzzzzz_1, \
                                         tg_xxx_xxyyyyy_1, tg_xxx_xxyyyyyy_0, tg_xxx_xxyyyyyy_1, tg_xxx_xxyyyyyz_0, \
                                         tg_xxx_xxyyyyyz_1, tg_xxx_xxyyyyz_1, tg_xxx_xxyyyyzz_0, tg_xxx_xxyyyyzz_1, \
                                         tg_xxx_xxyyyzz_1, tg_xxx_xxyyyzzz_0, tg_xxx_xxyyyzzz_1, tg_xxx_xxyyzzz_1, \
                                         tg_xxx_xxyyzzzz_0, tg_xxx_xxyyzzzz_1, tg_xxx_xxyzzzz_1, tg_xxx_xxyzzzzz_0, \
                                         tg_xxx_xxyzzzzz_1, tg_xxx_xxzzzzz_1, tg_xxx_xxzzzzzz_0, tg_xxx_xxzzzzzz_1, \
                                         tg_xxx_xyyyyyy_1, tg_xxx_xyyyyyyy_0, tg_xxx_xyyyyyyy_1, tg_xxx_xyyyyyyz_0, \
                                         tg_xxx_xyyyyyyz_1, tg_xxx_xyyyyyz_1, tg_xxx_xyyyyyzz_0, tg_xxx_xyyyyyzz_1, \
                                         tg_xxx_xyyyyzz_1, tg_xxx_xyyyyzzz_0, tg_xxx_xyyyyzzz_1, tg_xxx_xyyyzzz_1, \
                                         tg_xxx_xyyyzzzz_0, tg_xxx_xyyyzzzz_1, tg_xxx_xyyzzzz_1, tg_xxx_xyyzzzzz_0, \
                                         tg_xxx_xyyzzzzz_1, tg_xxx_xyzzzzz_1, tg_xxx_xyzzzzzz_0, tg_xxx_xyzzzzzz_1, \
                                         tg_xxx_xzzzzzz_1, tg_xxx_xzzzzzzz_0, tg_xxx_xzzzzzzz_1, tg_xxx_yyyyyyy_1, \
                                         tg_xxx_yyyyyyyy_0, tg_xxx_yyyyyyyy_1, tg_xxx_yyyyyyyz_0, tg_xxx_yyyyyyyz_1, \
                                         tg_xxx_yyyyyyz_1, tg_xxx_yyyyyyzz_0, tg_xxx_yyyyyyzz_1, tg_xxx_yyyyyzz_1, \
                                         tg_xxx_yyyyyzzz_0, tg_xxx_yyyyyzzz_1, tg_xxx_yyyyzzz_1, tg_xxx_yyyyzzzz_0, \
                                         tg_xxx_yyyyzzzz_1, tg_xxx_yyyzzzz_1, tg_xxx_yyyzzzzz_0, tg_xxx_yyyzzzzz_1, \
                                         tg_xxx_yyzzzzz_1, tg_xxx_yyzzzzzz_0, tg_xxx_yyzzzzzz_1, tg_xxx_yzzzzzz_1, \
                                         tg_xxx_yzzzzzzz_0, tg_xxx_yzzzzzzz_1, tg_xxx_zzzzzzz_1, tg_xxx_zzzzzzzz_0, \
                                         tg_xxx_zzzzzzzz_1, tg_xxxx_xxxxxxxx_0, tg_xxxx_xxxxxxxy_0, tg_xxxx_xxxxxxxz_0, \
                                         tg_xxxx_xxxxxxyy_0, tg_xxxx_xxxxxxyz_0, tg_xxxx_xxxxxxzz_0, tg_xxxx_xxxxxyyy_0, \
                                         tg_xxxx_xxxxxyyz_0, tg_xxxx_xxxxxyzz_0, tg_xxxx_xxxxxzzz_0, tg_xxxx_xxxxyyyy_0, \
                                         tg_xxxx_xxxxyyyz_0, tg_xxxx_xxxxyyzz_0, tg_xxxx_xxxxyzzz_0, tg_xxxx_xxxxzzzz_0, \
                                         tg_xxxx_xxxyyyyy_0, tg_xxxx_xxxyyyyz_0, tg_xxxx_xxxyyyzz_0, tg_xxxx_xxxyyzzz_0, \
                                         tg_xxxx_xxxyzzzz_0, tg_xxxx_xxxzzzzz_0, tg_xxxx_xxyyyyyy_0, tg_xxxx_xxyyyyyz_0, \
                                         tg_xxxx_xxyyyyzz_0, tg_xxxx_xxyyyzzz_0, tg_xxxx_xxyyzzzz_0, tg_xxxx_xxyzzzzz_0, \
                                         tg_xxxx_xxzzzzzz_0, tg_xxxx_xyyyyyyy_0, tg_xxxx_xyyyyyyz_0, tg_xxxx_xyyyyyzz_0, \
                                         tg_xxxx_xyyyyzzz_0, tg_xxxx_xyyyzzzz_0, tg_xxxx_xyyzzzzz_0, tg_xxxx_xyzzzzzz_0, \
                                         tg_xxxx_xzzzzzzz_0, tg_xxxx_yyyyyyyy_0, tg_xxxx_yyyyyyyz_0, tg_xxxx_yyyyyyzz_0, \
                                         tg_xxxx_yyyyyzzz_0, tg_xxxx_yyyyzzzz_0, tg_xxxx_yyyzzzzz_0, tg_xxxx_yyzzzzzz_0, \
                                         tg_xxxx_yzzzzzzz_0, tg_xxxx_zzzzzzzz_0, tg_xxxy_xxxxxxxx_0, tg_xxxy_xxxxxxxy_0, \
                                         tg_xxxy_xxxxxxxz_0, tg_xxxy_xxxxxxyy_0, tg_xxy_xxxxxxx_1, tg_xxy_xxxxxxxx_0, \
                                         tg_xxy_xxxxxxxx_1, tg_xxy_xxxxxxxy_0, tg_xxy_xxxxxxxy_1, tg_xxy_xxxxxxxz_0, \
                                         tg_xxy_xxxxxxxz_1, tg_xxy_xxxxxxy_1, tg_xxy_xxxxxxyy_0, tg_xxy_xxxxxxyy_1, \
                                         tg_xxy_xxxxxxz_1, tg_xxy_xxxxxyy_1, tg_xy_xxxxxxxx_0, tg_xy_xxxxxxxx_1, \
                                         tg_xy_xxxxxxxy_0, tg_xy_xxxxxxxy_1, tg_xy_xxxxxxxz_0, tg_xy_xxxxxxxz_1, \
                                         tg_xy_xxxxxxyy_0, tg_xy_xxxxxxyy_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxxx_xxxxxxxx_0[j] = pb_x * tg_xxx_xxxxxxxx_0[j] + wp_x[j] * tg_xxx_xxxxxxxx_1[j] + 1.5 * fl1_fx * tg_xx_xxxxxxxx_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxxxxxxx_1[j] + 4.0 * fl1_fxn * tg_xxx_xxxxxxx_1[j];

                    tg_xxxx_xxxxxxxy_0[j] = pb_x * tg_xxx_xxxxxxxy_0[j] + wp_x[j] * tg_xxx_xxxxxxxy_1[j] + 1.5 * fl1_fx * tg_xx_xxxxxxxy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxxxxxxy_1[j] + 3.5 * fl1_fxn * tg_xxx_xxxxxxy_1[j];

                    tg_xxxx_xxxxxxxz_0[j] = pb_x * tg_xxx_xxxxxxxz_0[j] + wp_x[j] * tg_xxx_xxxxxxxz_1[j] + 1.5 * fl1_fx * tg_xx_xxxxxxxz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxxxxxxz_1[j] + 3.5 * fl1_fxn * tg_xxx_xxxxxxz_1[j];

                    tg_xxxx_xxxxxxyy_0[j] = pb_x * tg_xxx_xxxxxxyy_0[j] + wp_x[j] * tg_xxx_xxxxxxyy_1[j] + 1.5 * fl1_fx * tg_xx_xxxxxxyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxxxxxyy_1[j] + 3.0 * fl1_fxn * tg_xxx_xxxxxyy_1[j];

                    tg_xxxx_xxxxxxyz_0[j] = pb_x * tg_xxx_xxxxxxyz_0[j] + wp_x[j] * tg_xxx_xxxxxxyz_1[j] + 1.5 * fl1_fx * tg_xx_xxxxxxyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxxxxxyz_1[j] + 3.0 * fl1_fxn * tg_xxx_xxxxxyz_1[j];

                    tg_xxxx_xxxxxxzz_0[j] = pb_x * tg_xxx_xxxxxxzz_0[j] + wp_x[j] * tg_xxx_xxxxxxzz_1[j] + 1.5 * fl1_fx * tg_xx_xxxxxxzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxxxxxzz_1[j] + 3.0 * fl1_fxn * tg_xxx_xxxxxzz_1[j];

                    tg_xxxx_xxxxxyyy_0[j] = pb_x * tg_xxx_xxxxxyyy_0[j] + wp_x[j] * tg_xxx_xxxxxyyy_1[j] + 1.5 * fl1_fx * tg_xx_xxxxxyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxxxxyyy_1[j] + 2.5 * fl1_fxn * tg_xxx_xxxxyyy_1[j];

                    tg_xxxx_xxxxxyyz_0[j] = pb_x * tg_xxx_xxxxxyyz_0[j] + wp_x[j] * tg_xxx_xxxxxyyz_1[j] + 1.5 * fl1_fx * tg_xx_xxxxxyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxxxxyyz_1[j] + 2.5 * fl1_fxn * tg_xxx_xxxxyyz_1[j];

                    tg_xxxx_xxxxxyzz_0[j] = pb_x * tg_xxx_xxxxxyzz_0[j] + wp_x[j] * tg_xxx_xxxxxyzz_1[j] + 1.5 * fl1_fx * tg_xx_xxxxxyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxxxxyzz_1[j] + 2.5 * fl1_fxn * tg_xxx_xxxxyzz_1[j];

                    tg_xxxx_xxxxxzzz_0[j] = pb_x * tg_xxx_xxxxxzzz_0[j] + wp_x[j] * tg_xxx_xxxxxzzz_1[j] + 1.5 * fl1_fx * tg_xx_xxxxxzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxxxxzzz_1[j] + 2.5 * fl1_fxn * tg_xxx_xxxxzzz_1[j];

                    tg_xxxx_xxxxyyyy_0[j] = pb_x * tg_xxx_xxxxyyyy_0[j] + wp_x[j] * tg_xxx_xxxxyyyy_1[j] + 1.5 * fl1_fx * tg_xx_xxxxyyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_xxx_xxxyyyy_1[j];

                    tg_xxxx_xxxxyyyz_0[j] = pb_x * tg_xxx_xxxxyyyz_0[j] + wp_x[j] * tg_xxx_xxxxyyyz_1[j] + 1.5 * fl1_fx * tg_xx_xxxxyyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxxxyyyz_1[j] + 2.0 * fl1_fxn * tg_xxx_xxxyyyz_1[j];

                    tg_xxxx_xxxxyyzz_0[j] = pb_x * tg_xxx_xxxxyyzz_0[j] + wp_x[j] * tg_xxx_xxxxyyzz_1[j] + 1.5 * fl1_fx * tg_xx_xxxxyyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxxxyyzz_1[j] + 2.0 * fl1_fxn * tg_xxx_xxxyyzz_1[j];

                    tg_xxxx_xxxxyzzz_0[j] = pb_x * tg_xxx_xxxxyzzz_0[j] + wp_x[j] * tg_xxx_xxxxyzzz_1[j] + 1.5 * fl1_fx * tg_xx_xxxxyzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxxxyzzz_1[j] + 2.0 * fl1_fxn * tg_xxx_xxxyzzz_1[j];

                    tg_xxxx_xxxxzzzz_0[j] = pb_x * tg_xxx_xxxxzzzz_0[j] + wp_x[j] * tg_xxx_xxxxzzzz_1[j] + 1.5 * fl1_fx * tg_xx_xxxxzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxxxzzzz_1[j] + 2.0 * fl1_fxn * tg_xxx_xxxzzzz_1[j];

                    tg_xxxx_xxxyyyyy_0[j] = pb_x * tg_xxx_xxxyyyyy_0[j] + wp_x[j] * tg_xxx_xxxyyyyy_1[j] + 1.5 * fl1_fx * tg_xx_xxxyyyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxxyyyyy_1[j] + 1.5 * fl1_fxn * tg_xxx_xxyyyyy_1[j];

                    tg_xxxx_xxxyyyyz_0[j] = pb_x * tg_xxx_xxxyyyyz_0[j] + wp_x[j] * tg_xxx_xxxyyyyz_1[j] + 1.5 * fl1_fx * tg_xx_xxxyyyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxxyyyyz_1[j] + 1.5 * fl1_fxn * tg_xxx_xxyyyyz_1[j];

                    tg_xxxx_xxxyyyzz_0[j] = pb_x * tg_xxx_xxxyyyzz_0[j] + wp_x[j] * tg_xxx_xxxyyyzz_1[j] + 1.5 * fl1_fx * tg_xx_xxxyyyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_xxx_xxyyyzz_1[j];

                    tg_xxxx_xxxyyzzz_0[j] = pb_x * tg_xxx_xxxyyzzz_0[j] + wp_x[j] * tg_xxx_xxxyyzzz_1[j] + 1.5 * fl1_fx * tg_xx_xxxyyzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxxyyzzz_1[j] + 1.5 * fl1_fxn * tg_xxx_xxyyzzz_1[j];

                    tg_xxxx_xxxyzzzz_0[j] = pb_x * tg_xxx_xxxyzzzz_0[j] + wp_x[j] * tg_xxx_xxxyzzzz_1[j] + 1.5 * fl1_fx * tg_xx_xxxyzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxxyzzzz_1[j] + 1.5 * fl1_fxn * tg_xxx_xxyzzzz_1[j];

                    tg_xxxx_xxxzzzzz_0[j] = pb_x * tg_xxx_xxxzzzzz_0[j] + wp_x[j] * tg_xxx_xxxzzzzz_1[j] + 1.5 * fl1_fx * tg_xx_xxxzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxxzzzzz_1[j] + 1.5 * fl1_fxn * tg_xxx_xxzzzzz_1[j];

                    tg_xxxx_xxyyyyyy_0[j] = pb_x * tg_xxx_xxyyyyyy_0[j] + wp_x[j] * tg_xxx_xxyyyyyy_1[j] + 1.5 * fl1_fx * tg_xx_xxyyyyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxyyyyyy_1[j] + fl1_fxn * tg_xxx_xyyyyyy_1[j];

                    tg_xxxx_xxyyyyyz_0[j] = pb_x * tg_xxx_xxyyyyyz_0[j] + wp_x[j] * tg_xxx_xxyyyyyz_1[j] + 1.5 * fl1_fx * tg_xx_xxyyyyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxyyyyyz_1[j] + fl1_fxn * tg_xxx_xyyyyyz_1[j];

                    tg_xxxx_xxyyyyzz_0[j] = pb_x * tg_xxx_xxyyyyzz_0[j] + wp_x[j] * tg_xxx_xxyyyyzz_1[j] + 1.5 * fl1_fx * tg_xx_xxyyyyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxyyyyzz_1[j] + fl1_fxn * tg_xxx_xyyyyzz_1[j];

                    tg_xxxx_xxyyyzzz_0[j] = pb_x * tg_xxx_xxyyyzzz_0[j] + wp_x[j] * tg_xxx_xxyyyzzz_1[j] + 1.5 * fl1_fx * tg_xx_xxyyyzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxyyyzzz_1[j] + fl1_fxn * tg_xxx_xyyyzzz_1[j];

                    tg_xxxx_xxyyzzzz_0[j] = pb_x * tg_xxx_xxyyzzzz_0[j] + wp_x[j] * tg_xxx_xxyyzzzz_1[j] + 1.5 * fl1_fx * tg_xx_xxyyzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxyyzzzz_1[j] + fl1_fxn * tg_xxx_xyyzzzz_1[j];

                    tg_xxxx_xxyzzzzz_0[j] = pb_x * tg_xxx_xxyzzzzz_0[j] + wp_x[j] * tg_xxx_xxyzzzzz_1[j] + 1.5 * fl1_fx * tg_xx_xxyzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxyzzzzz_1[j] + fl1_fxn * tg_xxx_xyzzzzz_1[j];

                    tg_xxxx_xxzzzzzz_0[j] = pb_x * tg_xxx_xxzzzzzz_0[j] + wp_x[j] * tg_xxx_xxzzzzzz_1[j] + 1.5 * fl1_fx * tg_xx_xxzzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxzzzzzz_1[j] + fl1_fxn * tg_xxx_xzzzzzz_1[j];

                    tg_xxxx_xyyyyyyy_0[j] = pb_x * tg_xxx_xyyyyyyy_0[j] + wp_x[j] * tg_xxx_xyyyyyyy_1[j] + 1.5 * fl1_fx * tg_xx_xyyyyyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xyyyyyyy_1[j] + 0.5 * fl1_fxn * tg_xxx_yyyyyyy_1[j];

                    tg_xxxx_xyyyyyyz_0[j] = pb_x * tg_xxx_xyyyyyyz_0[j] + wp_x[j] * tg_xxx_xyyyyyyz_1[j] + 1.5 * fl1_fx * tg_xx_xyyyyyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xyyyyyyz_1[j] + 0.5 * fl1_fxn * tg_xxx_yyyyyyz_1[j];

                    tg_xxxx_xyyyyyzz_0[j] = pb_x * tg_xxx_xyyyyyzz_0[j] + wp_x[j] * tg_xxx_xyyyyyzz_1[j] + 1.5 * fl1_fx * tg_xx_xyyyyyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xyyyyyzz_1[j] + 0.5 * fl1_fxn * tg_xxx_yyyyyzz_1[j];

                    tg_xxxx_xyyyyzzz_0[j] = pb_x * tg_xxx_xyyyyzzz_0[j] + wp_x[j] * tg_xxx_xyyyyzzz_1[j] + 1.5 * fl1_fx * tg_xx_xyyyyzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xyyyyzzz_1[j] + 0.5 * fl1_fxn * tg_xxx_yyyyzzz_1[j];

                    tg_xxxx_xyyyzzzz_0[j] = pb_x * tg_xxx_xyyyzzzz_0[j] + wp_x[j] * tg_xxx_xyyyzzzz_1[j] + 1.5 * fl1_fx * tg_xx_xyyyzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xyyyzzzz_1[j] + 0.5 * fl1_fxn * tg_xxx_yyyzzzz_1[j];

                    tg_xxxx_xyyzzzzz_0[j] = pb_x * tg_xxx_xyyzzzzz_0[j] + wp_x[j] * tg_xxx_xyyzzzzz_1[j] + 1.5 * fl1_fx * tg_xx_xyyzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xyyzzzzz_1[j] + 0.5 * fl1_fxn * tg_xxx_yyzzzzz_1[j];

                    tg_xxxx_xyzzzzzz_0[j] = pb_x * tg_xxx_xyzzzzzz_0[j] + wp_x[j] * tg_xxx_xyzzzzzz_1[j] + 1.5 * fl1_fx * tg_xx_xyzzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xxx_yzzzzzz_1[j];

                    tg_xxxx_xzzzzzzz_0[j] = pb_x * tg_xxx_xzzzzzzz_0[j] + wp_x[j] * tg_xxx_xzzzzzzz_1[j] + 1.5 * fl1_fx * tg_xx_xzzzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xxx_zzzzzzz_1[j];

                    tg_xxxx_yyyyyyyy_0[j] = pb_x * tg_xxx_yyyyyyyy_0[j] + wp_x[j] * tg_xxx_yyyyyyyy_1[j] + 1.5 * fl1_fx * tg_xx_yyyyyyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_yyyyyyyy_1[j];

                    tg_xxxx_yyyyyyyz_0[j] = pb_x * tg_xxx_yyyyyyyz_0[j] + wp_x[j] * tg_xxx_yyyyyyyz_1[j] + 1.5 * fl1_fx * tg_xx_yyyyyyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_yyyyyyyz_1[j];

                    tg_xxxx_yyyyyyzz_0[j] = pb_x * tg_xxx_yyyyyyzz_0[j] + wp_x[j] * tg_xxx_yyyyyyzz_1[j] + 1.5 * fl1_fx * tg_xx_yyyyyyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_yyyyyyzz_1[j];

                    tg_xxxx_yyyyyzzz_0[j] = pb_x * tg_xxx_yyyyyzzz_0[j] + wp_x[j] * tg_xxx_yyyyyzzz_1[j] + 1.5 * fl1_fx * tg_xx_yyyyyzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_yyyyyzzz_1[j];

                    tg_xxxx_yyyyzzzz_0[j] = pb_x * tg_xxx_yyyyzzzz_0[j] + wp_x[j] * tg_xxx_yyyyzzzz_1[j] + 1.5 * fl1_fx * tg_xx_yyyyzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_yyyyzzzz_1[j];

                    tg_xxxx_yyyzzzzz_0[j] = pb_x * tg_xxx_yyyzzzzz_0[j] + wp_x[j] * tg_xxx_yyyzzzzz_1[j] + 1.5 * fl1_fx * tg_xx_yyyzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_yyyzzzzz_1[j];

                    tg_xxxx_yyzzzzzz_0[j] = pb_x * tg_xxx_yyzzzzzz_0[j] + wp_x[j] * tg_xxx_yyzzzzzz_1[j] + 1.5 * fl1_fx * tg_xx_yyzzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_yyzzzzzz_1[j];

                    tg_xxxx_yzzzzzzz_0[j] = pb_x * tg_xxx_yzzzzzzz_0[j] + wp_x[j] * tg_xxx_yzzzzzzz_1[j] + 1.5 * fl1_fx * tg_xx_yzzzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_yzzzzzzz_1[j];

                    tg_xxxx_zzzzzzzz_0[j] = pb_x * tg_xxx_zzzzzzzz_0[j] + wp_x[j] * tg_xxx_zzzzzzzz_1[j] + 1.5 * fl1_fx * tg_xx_zzzzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_zzzzzzzz_1[j];

                    tg_xxxy_xxxxxxxx_0[j] = pb_x * tg_xxy_xxxxxxxx_0[j] + wp_x[j] * tg_xxy_xxxxxxxx_1[j] + fl1_fx * tg_xy_xxxxxxxx_0[j] - fl1_fx * fl1_fza * tg_xy_xxxxxxxx_1[j] + 4.0 * fl1_fxn * tg_xxy_xxxxxxx_1[j];

                    tg_xxxy_xxxxxxxy_0[j] = pb_x * tg_xxy_xxxxxxxy_0[j] + wp_x[j] * tg_xxy_xxxxxxxy_1[j] + fl1_fx * tg_xy_xxxxxxxy_0[j] - fl1_fx * fl1_fza * tg_xy_xxxxxxxy_1[j] + 3.5 * fl1_fxn * tg_xxy_xxxxxxy_1[j];

                    tg_xxxy_xxxxxxxz_0[j] = pb_x * tg_xxy_xxxxxxxz_0[j] + wp_x[j] * tg_xxy_xxxxxxxz_1[j] + fl1_fx * tg_xy_xxxxxxxz_0[j] - fl1_fx * fl1_fza * tg_xy_xxxxxxxz_1[j] + 3.5 * fl1_fxn * tg_xxy_xxxxxxz_1[j];

                    tg_xxxy_xxxxxxyy_0[j] = pb_x * tg_xxy_xxxxxxyy_0[j] + wp_x[j] * tg_xxy_xxxxxxyy_1[j] + fl1_fx * tg_xy_xxxxxxyy_0[j] - fl1_fx * fl1_fza * tg_xy_xxxxxxyy_1[j] + 3.0 * fl1_fxn * tg_xxy_xxxxxyy_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSL_49_98(      CMemBlock2D<double>& primBuffer,
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
                                             {4, -1, -1, -1},
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_2_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_2_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

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

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                // set up pointers to auxilary integrals

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

                auto tg_xxz_xxxxxxxx_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 90); 

                auto tg_xxz_xxxxxxxy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 91); 

                auto tg_xxz_xxxxxxxz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 92); 

                auto tg_xxz_xxxxxxyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 93); 

                auto tg_xxz_xxxxxxyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 94); 

                auto tg_xxz_xxxxxxzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 95); 

                auto tg_xxz_xxxxxyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 96); 

                auto tg_xxz_xxxxxyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 97); 

                auto tg_xxy_xxxxxxyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 49); 

                auto tg_xxy_xxxxxxzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 50); 

                auto tg_xxy_xxxxxyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 51); 

                auto tg_xxy_xxxxxyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 52); 

                auto tg_xxy_xxxxxyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 53); 

                auto tg_xxy_xxxxxzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 54); 

                auto tg_xxy_xxxxyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 55); 

                auto tg_xxy_xxxxyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 56); 

                auto tg_xxy_xxxxyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 57); 

                auto tg_xxy_xxxxyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 58); 

                auto tg_xxy_xxxxzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 59); 

                auto tg_xxy_xxxyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 60); 

                auto tg_xxy_xxxyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 61); 

                auto tg_xxy_xxxyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 62); 

                auto tg_xxy_xxxyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 63); 

                auto tg_xxy_xxxyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 64); 

                auto tg_xxy_xxxzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 65); 

                auto tg_xxy_xxyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 66); 

                auto tg_xxy_xxyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 67); 

                auto tg_xxy_xxyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 68); 

                auto tg_xxy_xxyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 69); 

                auto tg_xxy_xxyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 70); 

                auto tg_xxy_xxyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 71); 

                auto tg_xxy_xxzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 72); 

                auto tg_xxy_xyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 73); 

                auto tg_xxy_xyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 74); 

                auto tg_xxy_xyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 75); 

                auto tg_xxy_xyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 76); 

                auto tg_xxy_xyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 77); 

                auto tg_xxy_xyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 78); 

                auto tg_xxy_xyzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 79); 

                auto tg_xxy_xzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 80); 

                auto tg_xxy_yyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 81); 

                auto tg_xxy_yyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 82); 

                auto tg_xxy_yyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 83); 

                auto tg_xxy_yyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 84); 

                auto tg_xxy_yyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 85); 

                auto tg_xxy_yyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 86); 

                auto tg_xxy_yyzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 87); 

                auto tg_xxy_yzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 88); 

                auto tg_xxy_zzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 89); 

                auto tg_xxz_xxxxxxxx_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 90); 

                auto tg_xxz_xxxxxxxy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 91); 

                auto tg_xxz_xxxxxxxz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 92); 

                auto tg_xxz_xxxxxxyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 93); 

                auto tg_xxz_xxxxxxyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 94); 

                auto tg_xxz_xxxxxxzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 95); 

                auto tg_xxz_xxxxxyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 96); 

                auto tg_xxz_xxxxxyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 97); 

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

                auto tg_xz_xxxxxxxx_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 90); 

                auto tg_xz_xxxxxxxy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 91); 

                auto tg_xz_xxxxxxxz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 92); 

                auto tg_xz_xxxxxxyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 93); 

                auto tg_xz_xxxxxxyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 94); 

                auto tg_xz_xxxxxxzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 95); 

                auto tg_xz_xxxxxyyy_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 96); 

                auto tg_xz_xxxxxyyz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 97); 

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

                auto tg_xz_xxxxxxxx_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 90); 

                auto tg_xz_xxxxxxxy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 91); 

                auto tg_xz_xxxxxxxz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 92); 

                auto tg_xz_xxxxxxyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 93); 

                auto tg_xz_xxxxxxyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 94); 

                auto tg_xz_xxxxxxzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 95); 

                auto tg_xz_xxxxxyyy_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 96); 

                auto tg_xz_xxxxxyyz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 97); 

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

                // set up pointers to integrals

                auto tg_xxxy_xxxxxxyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 49); 

                auto tg_xxxy_xxxxxxzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 50); 

                auto tg_xxxy_xxxxxyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 51); 

                auto tg_xxxy_xxxxxyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 52); 

                auto tg_xxxy_xxxxxyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 53); 

                auto tg_xxxy_xxxxxzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 54); 

                auto tg_xxxy_xxxxyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 55); 

                auto tg_xxxy_xxxxyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 56); 

                auto tg_xxxy_xxxxyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 57); 

                auto tg_xxxy_xxxxyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 58); 

                auto tg_xxxy_xxxxzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 59); 

                auto tg_xxxy_xxxyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 60); 

                auto tg_xxxy_xxxyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 61); 

                auto tg_xxxy_xxxyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 62); 

                auto tg_xxxy_xxxyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 63); 

                auto tg_xxxy_xxxyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 64); 

                auto tg_xxxy_xxxzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 65); 

                auto tg_xxxy_xxyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 66); 

                auto tg_xxxy_xxyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 67); 

                auto tg_xxxy_xxyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 68); 

                auto tg_xxxy_xxyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 69); 

                auto tg_xxxy_xxyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 70); 

                auto tg_xxxy_xxyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 71); 

                auto tg_xxxy_xxzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 72); 

                auto tg_xxxy_xyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 73); 

                auto tg_xxxy_xyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 74); 

                auto tg_xxxy_xyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 75); 

                auto tg_xxxy_xyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 76); 

                auto tg_xxxy_xyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 77); 

                auto tg_xxxy_xyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 78); 

                auto tg_xxxy_xyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 79); 

                auto tg_xxxy_xzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 80); 

                auto tg_xxxy_yyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 81); 

                auto tg_xxxy_yyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 82); 

                auto tg_xxxy_yyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 83); 

                auto tg_xxxy_yyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 84); 

                auto tg_xxxy_yyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 85); 

                auto tg_xxxy_yyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 86); 

                auto tg_xxxy_yyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 87); 

                auto tg_xxxy_yzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 88); 

                auto tg_xxxy_zzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 89); 

                auto tg_xxxz_xxxxxxxx_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 90); 

                auto tg_xxxz_xxxxxxxy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 91); 

                auto tg_xxxz_xxxxxxxz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 92); 

                auto tg_xxxz_xxxxxxyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 93); 

                auto tg_xxxz_xxxxxxyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 94); 

                auto tg_xxxz_xxxxxxzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 95); 

                auto tg_xxxz_xxxxxyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 96); 

                auto tg_xxxz_xxxxxyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 97); 

                // Batch of Integrals (49,98)

                #pragma omp simd aligned(fxn, fza, tg_xxxy_xxxxxxyz_0, tg_xxxy_xxxxxxzz_0, tg_xxxy_xxxxxyyy_0, \
                                         tg_xxxy_xxxxxyyz_0, tg_xxxy_xxxxxyzz_0, tg_xxxy_xxxxxzzz_0, tg_xxxy_xxxxyyyy_0, \
                                         tg_xxxy_xxxxyyyz_0, tg_xxxy_xxxxyyzz_0, tg_xxxy_xxxxyzzz_0, tg_xxxy_xxxxzzzz_0, \
                                         tg_xxxy_xxxyyyyy_0, tg_xxxy_xxxyyyyz_0, tg_xxxy_xxxyyyzz_0, tg_xxxy_xxxyyzzz_0, \
                                         tg_xxxy_xxxyzzzz_0, tg_xxxy_xxxzzzzz_0, tg_xxxy_xxyyyyyy_0, tg_xxxy_xxyyyyyz_0, \
                                         tg_xxxy_xxyyyyzz_0, tg_xxxy_xxyyyzzz_0, tg_xxxy_xxyyzzzz_0, tg_xxxy_xxyzzzzz_0, \
                                         tg_xxxy_xxzzzzzz_0, tg_xxxy_xyyyyyyy_0, tg_xxxy_xyyyyyyz_0, tg_xxxy_xyyyyyzz_0, \
                                         tg_xxxy_xyyyyzzz_0, tg_xxxy_xyyyzzzz_0, tg_xxxy_xyyzzzzz_0, tg_xxxy_xyzzzzzz_0, \
                                         tg_xxxy_xzzzzzzz_0, tg_xxxy_yyyyyyyy_0, tg_xxxy_yyyyyyyz_0, tg_xxxy_yyyyyyzz_0, \
                                         tg_xxxy_yyyyyzzz_0, tg_xxxy_yyyyzzzz_0, tg_xxxy_yyyzzzzz_0, tg_xxxy_yyzzzzzz_0, \
                                         tg_xxxy_yzzzzzzz_0, tg_xxxy_zzzzzzzz_0, tg_xxxz_xxxxxxxx_0, tg_xxxz_xxxxxxxy_0, \
                                         tg_xxxz_xxxxxxxz_0, tg_xxxz_xxxxxxyy_0, tg_xxxz_xxxxxxyz_0, tg_xxxz_xxxxxxzz_0, \
                                         tg_xxxz_xxxxxyyy_0, tg_xxxz_xxxxxyyz_0, tg_xxy_xxxxxxyz_0, tg_xxy_xxxxxxyz_1, \
                                         tg_xxy_xxxxxxzz_0, tg_xxy_xxxxxxzz_1, tg_xxy_xxxxxyyy_0, tg_xxy_xxxxxyyy_1, \
                                         tg_xxy_xxxxxyyz_0, tg_xxy_xxxxxyyz_1, tg_xxy_xxxxxyz_1, tg_xxy_xxxxxyzz_0, \
                                         tg_xxy_xxxxxyzz_1, tg_xxy_xxxxxzz_1, tg_xxy_xxxxxzzz_0, tg_xxy_xxxxxzzz_1, \
                                         tg_xxy_xxxxyyy_1, tg_xxy_xxxxyyyy_0, tg_xxy_xxxxyyyy_1, tg_xxy_xxxxyyyz_0, \
                                         tg_xxy_xxxxyyyz_1, tg_xxy_xxxxyyz_1, tg_xxy_xxxxyyzz_0, tg_xxy_xxxxyyzz_1, \
                                         tg_xxy_xxxxyzz_1, tg_xxy_xxxxyzzz_0, tg_xxy_xxxxyzzz_1, tg_xxy_xxxxzzz_1, \
                                         tg_xxy_xxxxzzzz_0, tg_xxy_xxxxzzzz_1, tg_xxy_xxxyyyy_1, tg_xxy_xxxyyyyy_0, \
                                         tg_xxy_xxxyyyyy_1, tg_xxy_xxxyyyyz_0, tg_xxy_xxxyyyyz_1, tg_xxy_xxxyyyz_1, \
                                         tg_xxy_xxxyyyzz_0, tg_xxy_xxxyyyzz_1, tg_xxy_xxxyyzz_1, tg_xxy_xxxyyzzz_0, \
                                         tg_xxy_xxxyyzzz_1, tg_xxy_xxxyzzz_1, tg_xxy_xxxyzzzz_0, tg_xxy_xxxyzzzz_1, \
                                         tg_xxy_xxxzzzz_1, tg_xxy_xxxzzzzz_0, tg_xxy_xxxzzzzz_1, tg_xxy_xxyyyyy_1, \
                                         tg_xxy_xxyyyyyy_0, tg_xxy_xxyyyyyy_1, tg_xxy_xxyyyyyz_0, tg_xxy_xxyyyyyz_1, \
                                         tg_xxy_xxyyyyz_1, tg_xxy_xxyyyyzz_0, tg_xxy_xxyyyyzz_1, tg_xxy_xxyyyzz_1, \
                                         tg_xxy_xxyyyzzz_0, tg_xxy_xxyyyzzz_1, tg_xxy_xxyyzzz_1, tg_xxy_xxyyzzzz_0, \
                                         tg_xxy_xxyyzzzz_1, tg_xxy_xxyzzzz_1, tg_xxy_xxyzzzzz_0, tg_xxy_xxyzzzzz_1, \
                                         tg_xxy_xxzzzzz_1, tg_xxy_xxzzzzzz_0, tg_xxy_xxzzzzzz_1, tg_xxy_xyyyyyy_1, \
                                         tg_xxy_xyyyyyyy_0, tg_xxy_xyyyyyyy_1, tg_xxy_xyyyyyyz_0, tg_xxy_xyyyyyyz_1, \
                                         tg_xxy_xyyyyyz_1, tg_xxy_xyyyyyzz_0, tg_xxy_xyyyyyzz_1, tg_xxy_xyyyyzz_1, \
                                         tg_xxy_xyyyyzzz_0, tg_xxy_xyyyyzzz_1, tg_xxy_xyyyzzz_1, tg_xxy_xyyyzzzz_0, \
                                         tg_xxy_xyyyzzzz_1, tg_xxy_xyyzzzz_1, tg_xxy_xyyzzzzz_0, tg_xxy_xyyzzzzz_1, \
                                         tg_xxy_xyzzzzz_1, tg_xxy_xyzzzzzz_0, tg_xxy_xyzzzzzz_1, tg_xxy_xzzzzzz_1, \
                                         tg_xxy_xzzzzzzz_0, tg_xxy_xzzzzzzz_1, tg_xxy_yyyyyyy_1, tg_xxy_yyyyyyyy_0, \
                                         tg_xxy_yyyyyyyy_1, tg_xxy_yyyyyyyz_0, tg_xxy_yyyyyyyz_1, tg_xxy_yyyyyyz_1, \
                                         tg_xxy_yyyyyyzz_0, tg_xxy_yyyyyyzz_1, tg_xxy_yyyyyzz_1, tg_xxy_yyyyyzzz_0, \
                                         tg_xxy_yyyyyzzz_1, tg_xxy_yyyyzzz_1, tg_xxy_yyyyzzzz_0, tg_xxy_yyyyzzzz_1, \
                                         tg_xxy_yyyzzzz_1, tg_xxy_yyyzzzzz_0, tg_xxy_yyyzzzzz_1, tg_xxy_yyzzzzz_1, \
                                         tg_xxy_yyzzzzzz_0, tg_xxy_yyzzzzzz_1, tg_xxy_yzzzzzz_1, tg_xxy_yzzzzzzz_0, \
                                         tg_xxy_yzzzzzzz_1, tg_xxy_zzzzzzz_1, tg_xxy_zzzzzzzz_0, tg_xxy_zzzzzzzz_1, \
                                         tg_xxz_xxxxxxx_1, tg_xxz_xxxxxxxx_0, tg_xxz_xxxxxxxx_1, tg_xxz_xxxxxxxy_0, \
                                         tg_xxz_xxxxxxxy_1, tg_xxz_xxxxxxxz_0, tg_xxz_xxxxxxxz_1, tg_xxz_xxxxxxy_1, \
                                         tg_xxz_xxxxxxyy_0, tg_xxz_xxxxxxyy_1, tg_xxz_xxxxxxyz_0, tg_xxz_xxxxxxyz_1, \
                                         tg_xxz_xxxxxxz_1, tg_xxz_xxxxxxzz_0, tg_xxz_xxxxxxzz_1, tg_xxz_xxxxxyy_1, \
                                         tg_xxz_xxxxxyyy_0, tg_xxz_xxxxxyyy_1, tg_xxz_xxxxxyyz_0, tg_xxz_xxxxxyyz_1, \
                                         tg_xxz_xxxxxyz_1, tg_xxz_xxxxxzz_1, tg_xxz_xxxxyyy_1, tg_xxz_xxxxyyz_1, \
                                         tg_xy_xxxxxxyz_0, tg_xy_xxxxxxyz_1, tg_xy_xxxxxxzz_0, tg_xy_xxxxxxzz_1, \
                                         tg_xy_xxxxxyyy_0, tg_xy_xxxxxyyy_1, tg_xy_xxxxxyyz_0, tg_xy_xxxxxyyz_1, \
                                         tg_xy_xxxxxyzz_0, tg_xy_xxxxxyzz_1, tg_xy_xxxxxzzz_0, tg_xy_xxxxxzzz_1, \
                                         tg_xy_xxxxyyyy_0, tg_xy_xxxxyyyy_1, tg_xy_xxxxyyyz_0, tg_xy_xxxxyyyz_1, \
                                         tg_xy_xxxxyyzz_0, tg_xy_xxxxyyzz_1, tg_xy_xxxxyzzz_0, tg_xy_xxxxyzzz_1, \
                                         tg_xy_xxxxzzzz_0, tg_xy_xxxxzzzz_1, tg_xy_xxxyyyyy_0, tg_xy_xxxyyyyy_1, \
                                         tg_xy_xxxyyyyz_0, tg_xy_xxxyyyyz_1, tg_xy_xxxyyyzz_0, tg_xy_xxxyyyzz_1, \
                                         tg_xy_xxxyyzzz_0, tg_xy_xxxyyzzz_1, tg_xy_xxxyzzzz_0, tg_xy_xxxyzzzz_1, \
                                         tg_xy_xxxzzzzz_0, tg_xy_xxxzzzzz_1, tg_xy_xxyyyyyy_0, tg_xy_xxyyyyyy_1, \
                                         tg_xy_xxyyyyyz_0, tg_xy_xxyyyyyz_1, tg_xy_xxyyyyzz_0, tg_xy_xxyyyyzz_1, \
                                         tg_xy_xxyyyzzz_0, tg_xy_xxyyyzzz_1, tg_xy_xxyyzzzz_0, tg_xy_xxyyzzzz_1, \
                                         tg_xy_xxyzzzzz_0, tg_xy_xxyzzzzz_1, tg_xy_xxzzzzzz_0, tg_xy_xxzzzzzz_1, \
                                         tg_xy_xyyyyyyy_0, tg_xy_xyyyyyyy_1, tg_xy_xyyyyyyz_0, tg_xy_xyyyyyyz_1, \
                                         tg_xy_xyyyyyzz_0, tg_xy_xyyyyyzz_1, tg_xy_xyyyyzzz_0, tg_xy_xyyyyzzz_1, \
                                         tg_xy_xyyyzzzz_0, tg_xy_xyyyzzzz_1, tg_xy_xyyzzzzz_0, tg_xy_xyyzzzzz_1, \
                                         tg_xy_xyzzzzzz_0, tg_xy_xyzzzzzz_1, tg_xy_xzzzzzzz_0, tg_xy_xzzzzzzz_1, \
                                         tg_xy_yyyyyyyy_0, tg_xy_yyyyyyyy_1, tg_xy_yyyyyyyz_0, tg_xy_yyyyyyyz_1, \
                                         tg_xy_yyyyyyzz_0, tg_xy_yyyyyyzz_1, tg_xy_yyyyyzzz_0, tg_xy_yyyyyzzz_1, \
                                         tg_xy_yyyyzzzz_0, tg_xy_yyyyzzzz_1, tg_xy_yyyzzzzz_0, tg_xy_yyyzzzzz_1, \
                                         tg_xy_yyzzzzzz_0, tg_xy_yyzzzzzz_1, tg_xy_yzzzzzzz_0, tg_xy_yzzzzzzz_1, \
                                         tg_xy_zzzzzzzz_0, tg_xy_zzzzzzzz_1, tg_xz_xxxxxxxx_0, tg_xz_xxxxxxxx_1, \
                                         tg_xz_xxxxxxxy_0, tg_xz_xxxxxxxy_1, tg_xz_xxxxxxxz_0, tg_xz_xxxxxxxz_1, \
                                         tg_xz_xxxxxxyy_0, tg_xz_xxxxxxyy_1, tg_xz_xxxxxxyz_0, tg_xz_xxxxxxyz_1, \
                                         tg_xz_xxxxxxzz_0, tg_xz_xxxxxxzz_1, tg_xz_xxxxxyyy_0, tg_xz_xxxxxyyy_1, \
                                         tg_xz_xxxxxyyz_0, tg_xz_xxxxxyyz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxxy_xxxxxxyz_0[j] = pb_x * tg_xxy_xxxxxxyz_0[j] + wp_x[j] * tg_xxy_xxxxxxyz_1[j] + fl1_fx * tg_xy_xxxxxxyz_0[j] - fl1_fx * fl1_fza * tg_xy_xxxxxxyz_1[j] + 3.0 * fl1_fxn * tg_xxy_xxxxxyz_1[j];

                    tg_xxxy_xxxxxxzz_0[j] = pb_x * tg_xxy_xxxxxxzz_0[j] + wp_x[j] * tg_xxy_xxxxxxzz_1[j] + fl1_fx * tg_xy_xxxxxxzz_0[j] - fl1_fx * fl1_fza * tg_xy_xxxxxxzz_1[j] + 3.0 * fl1_fxn * tg_xxy_xxxxxzz_1[j];

                    tg_xxxy_xxxxxyyy_0[j] = pb_x * tg_xxy_xxxxxyyy_0[j] + wp_x[j] * tg_xxy_xxxxxyyy_1[j] + fl1_fx * tg_xy_xxxxxyyy_0[j] - fl1_fx * fl1_fza * tg_xy_xxxxxyyy_1[j] + 2.5 * fl1_fxn * tg_xxy_xxxxyyy_1[j];

                    tg_xxxy_xxxxxyyz_0[j] = pb_x * tg_xxy_xxxxxyyz_0[j] + wp_x[j] * tg_xxy_xxxxxyyz_1[j] + fl1_fx * tg_xy_xxxxxyyz_0[j] - fl1_fx * fl1_fza * tg_xy_xxxxxyyz_1[j] + 2.5 * fl1_fxn * tg_xxy_xxxxyyz_1[j];

                    tg_xxxy_xxxxxyzz_0[j] = pb_x * tg_xxy_xxxxxyzz_0[j] + wp_x[j] * tg_xxy_xxxxxyzz_1[j] + fl1_fx * tg_xy_xxxxxyzz_0[j] - fl1_fx * fl1_fza * tg_xy_xxxxxyzz_1[j] + 2.5 * fl1_fxn * tg_xxy_xxxxyzz_1[j];

                    tg_xxxy_xxxxxzzz_0[j] = pb_x * tg_xxy_xxxxxzzz_0[j] + wp_x[j] * tg_xxy_xxxxxzzz_1[j] + fl1_fx * tg_xy_xxxxxzzz_0[j] - fl1_fx * fl1_fza * tg_xy_xxxxxzzz_1[j] + 2.5 * fl1_fxn * tg_xxy_xxxxzzz_1[j];

                    tg_xxxy_xxxxyyyy_0[j] = pb_x * tg_xxy_xxxxyyyy_0[j] + wp_x[j] * tg_xxy_xxxxyyyy_1[j] + fl1_fx * tg_xy_xxxxyyyy_0[j] - fl1_fx * fl1_fza * tg_xy_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_xxy_xxxyyyy_1[j];

                    tg_xxxy_xxxxyyyz_0[j] = pb_x * tg_xxy_xxxxyyyz_0[j] + wp_x[j] * tg_xxy_xxxxyyyz_1[j] + fl1_fx * tg_xy_xxxxyyyz_0[j] - fl1_fx * fl1_fza * tg_xy_xxxxyyyz_1[j] + 2.0 * fl1_fxn * tg_xxy_xxxyyyz_1[j];

                    tg_xxxy_xxxxyyzz_0[j] = pb_x * tg_xxy_xxxxyyzz_0[j] + wp_x[j] * tg_xxy_xxxxyyzz_1[j] + fl1_fx * tg_xy_xxxxyyzz_0[j] - fl1_fx * fl1_fza * tg_xy_xxxxyyzz_1[j] + 2.0 * fl1_fxn * tg_xxy_xxxyyzz_1[j];

                    tg_xxxy_xxxxyzzz_0[j] = pb_x * tg_xxy_xxxxyzzz_0[j] + wp_x[j] * tg_xxy_xxxxyzzz_1[j] + fl1_fx * tg_xy_xxxxyzzz_0[j] - fl1_fx * fl1_fza * tg_xy_xxxxyzzz_1[j] + 2.0 * fl1_fxn * tg_xxy_xxxyzzz_1[j];

                    tg_xxxy_xxxxzzzz_0[j] = pb_x * tg_xxy_xxxxzzzz_0[j] + wp_x[j] * tg_xxy_xxxxzzzz_1[j] + fl1_fx * tg_xy_xxxxzzzz_0[j] - fl1_fx * fl1_fza * tg_xy_xxxxzzzz_1[j] + 2.0 * fl1_fxn * tg_xxy_xxxzzzz_1[j];

                    tg_xxxy_xxxyyyyy_0[j] = pb_x * tg_xxy_xxxyyyyy_0[j] + wp_x[j] * tg_xxy_xxxyyyyy_1[j] + fl1_fx * tg_xy_xxxyyyyy_0[j] - fl1_fx * fl1_fza * tg_xy_xxxyyyyy_1[j] + 1.5 * fl1_fxn * tg_xxy_xxyyyyy_1[j];

                    tg_xxxy_xxxyyyyz_0[j] = pb_x * tg_xxy_xxxyyyyz_0[j] + wp_x[j] * tg_xxy_xxxyyyyz_1[j] + fl1_fx * tg_xy_xxxyyyyz_0[j] - fl1_fx * fl1_fza * tg_xy_xxxyyyyz_1[j] + 1.5 * fl1_fxn * tg_xxy_xxyyyyz_1[j];

                    tg_xxxy_xxxyyyzz_0[j] = pb_x * tg_xxy_xxxyyyzz_0[j] + wp_x[j] * tg_xxy_xxxyyyzz_1[j] + fl1_fx * tg_xy_xxxyyyzz_0[j] - fl1_fx * fl1_fza * tg_xy_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_xxy_xxyyyzz_1[j];

                    tg_xxxy_xxxyyzzz_0[j] = pb_x * tg_xxy_xxxyyzzz_0[j] + wp_x[j] * tg_xxy_xxxyyzzz_1[j] + fl1_fx * tg_xy_xxxyyzzz_0[j] - fl1_fx * fl1_fza * tg_xy_xxxyyzzz_1[j] + 1.5 * fl1_fxn * tg_xxy_xxyyzzz_1[j];

                    tg_xxxy_xxxyzzzz_0[j] = pb_x * tg_xxy_xxxyzzzz_0[j] + wp_x[j] * tg_xxy_xxxyzzzz_1[j] + fl1_fx * tg_xy_xxxyzzzz_0[j] - fl1_fx * fl1_fza * tg_xy_xxxyzzzz_1[j] + 1.5 * fl1_fxn * tg_xxy_xxyzzzz_1[j];

                    tg_xxxy_xxxzzzzz_0[j] = pb_x * tg_xxy_xxxzzzzz_0[j] + wp_x[j] * tg_xxy_xxxzzzzz_1[j] + fl1_fx * tg_xy_xxxzzzzz_0[j] - fl1_fx * fl1_fza * tg_xy_xxxzzzzz_1[j] + 1.5 * fl1_fxn * tg_xxy_xxzzzzz_1[j];

                    tg_xxxy_xxyyyyyy_0[j] = pb_x * tg_xxy_xxyyyyyy_0[j] + wp_x[j] * tg_xxy_xxyyyyyy_1[j] + fl1_fx * tg_xy_xxyyyyyy_0[j] - fl1_fx * fl1_fza * tg_xy_xxyyyyyy_1[j] + fl1_fxn * tg_xxy_xyyyyyy_1[j];

                    tg_xxxy_xxyyyyyz_0[j] = pb_x * tg_xxy_xxyyyyyz_0[j] + wp_x[j] * tg_xxy_xxyyyyyz_1[j] + fl1_fx * tg_xy_xxyyyyyz_0[j] - fl1_fx * fl1_fza * tg_xy_xxyyyyyz_1[j] + fl1_fxn * tg_xxy_xyyyyyz_1[j];

                    tg_xxxy_xxyyyyzz_0[j] = pb_x * tg_xxy_xxyyyyzz_0[j] + wp_x[j] * tg_xxy_xxyyyyzz_1[j] + fl1_fx * tg_xy_xxyyyyzz_0[j] - fl1_fx * fl1_fza * tg_xy_xxyyyyzz_1[j] + fl1_fxn * tg_xxy_xyyyyzz_1[j];

                    tg_xxxy_xxyyyzzz_0[j] = pb_x * tg_xxy_xxyyyzzz_0[j] + wp_x[j] * tg_xxy_xxyyyzzz_1[j] + fl1_fx * tg_xy_xxyyyzzz_0[j] - fl1_fx * fl1_fza * tg_xy_xxyyyzzz_1[j] + fl1_fxn * tg_xxy_xyyyzzz_1[j];

                    tg_xxxy_xxyyzzzz_0[j] = pb_x * tg_xxy_xxyyzzzz_0[j] + wp_x[j] * tg_xxy_xxyyzzzz_1[j] + fl1_fx * tg_xy_xxyyzzzz_0[j] - fl1_fx * fl1_fza * tg_xy_xxyyzzzz_1[j] + fl1_fxn * tg_xxy_xyyzzzz_1[j];

                    tg_xxxy_xxyzzzzz_0[j] = pb_x * tg_xxy_xxyzzzzz_0[j] + wp_x[j] * tg_xxy_xxyzzzzz_1[j] + fl1_fx * tg_xy_xxyzzzzz_0[j] - fl1_fx * fl1_fza * tg_xy_xxyzzzzz_1[j] + fl1_fxn * tg_xxy_xyzzzzz_1[j];

                    tg_xxxy_xxzzzzzz_0[j] = pb_x * tg_xxy_xxzzzzzz_0[j] + wp_x[j] * tg_xxy_xxzzzzzz_1[j] + fl1_fx * tg_xy_xxzzzzzz_0[j] - fl1_fx * fl1_fza * tg_xy_xxzzzzzz_1[j] + fl1_fxn * tg_xxy_xzzzzzz_1[j];

                    tg_xxxy_xyyyyyyy_0[j] = pb_x * tg_xxy_xyyyyyyy_0[j] + wp_x[j] * tg_xxy_xyyyyyyy_1[j] + fl1_fx * tg_xy_xyyyyyyy_0[j] - fl1_fx * fl1_fza * tg_xy_xyyyyyyy_1[j] + 0.5 * fl1_fxn * tg_xxy_yyyyyyy_1[j];

                    tg_xxxy_xyyyyyyz_0[j] = pb_x * tg_xxy_xyyyyyyz_0[j] + wp_x[j] * tg_xxy_xyyyyyyz_1[j] + fl1_fx * tg_xy_xyyyyyyz_0[j] - fl1_fx * fl1_fza * tg_xy_xyyyyyyz_1[j] + 0.5 * fl1_fxn * tg_xxy_yyyyyyz_1[j];

                    tg_xxxy_xyyyyyzz_0[j] = pb_x * tg_xxy_xyyyyyzz_0[j] + wp_x[j] * tg_xxy_xyyyyyzz_1[j] + fl1_fx * tg_xy_xyyyyyzz_0[j] - fl1_fx * fl1_fza * tg_xy_xyyyyyzz_1[j] + 0.5 * fl1_fxn * tg_xxy_yyyyyzz_1[j];

                    tg_xxxy_xyyyyzzz_0[j] = pb_x * tg_xxy_xyyyyzzz_0[j] + wp_x[j] * tg_xxy_xyyyyzzz_1[j] + fl1_fx * tg_xy_xyyyyzzz_0[j] - fl1_fx * fl1_fza * tg_xy_xyyyyzzz_1[j] + 0.5 * fl1_fxn * tg_xxy_yyyyzzz_1[j];

                    tg_xxxy_xyyyzzzz_0[j] = pb_x * tg_xxy_xyyyzzzz_0[j] + wp_x[j] * tg_xxy_xyyyzzzz_1[j] + fl1_fx * tg_xy_xyyyzzzz_0[j] - fl1_fx * fl1_fza * tg_xy_xyyyzzzz_1[j] + 0.5 * fl1_fxn * tg_xxy_yyyzzzz_1[j];

                    tg_xxxy_xyyzzzzz_0[j] = pb_x * tg_xxy_xyyzzzzz_0[j] + wp_x[j] * tg_xxy_xyyzzzzz_1[j] + fl1_fx * tg_xy_xyyzzzzz_0[j] - fl1_fx * fl1_fza * tg_xy_xyyzzzzz_1[j] + 0.5 * fl1_fxn * tg_xxy_yyzzzzz_1[j];

                    tg_xxxy_xyzzzzzz_0[j] = pb_x * tg_xxy_xyzzzzzz_0[j] + wp_x[j] * tg_xxy_xyzzzzzz_1[j] + fl1_fx * tg_xy_xyzzzzzz_0[j] - fl1_fx * fl1_fza * tg_xy_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xxy_yzzzzzz_1[j];

                    tg_xxxy_xzzzzzzz_0[j] = pb_x * tg_xxy_xzzzzzzz_0[j] + wp_x[j] * tg_xxy_xzzzzzzz_1[j] + fl1_fx * tg_xy_xzzzzzzz_0[j] - fl1_fx * fl1_fza * tg_xy_xzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xxy_zzzzzzz_1[j];

                    tg_xxxy_yyyyyyyy_0[j] = pb_x * tg_xxy_yyyyyyyy_0[j] + wp_x[j] * tg_xxy_yyyyyyyy_1[j] + fl1_fx * tg_xy_yyyyyyyy_0[j] - fl1_fx * fl1_fza * tg_xy_yyyyyyyy_1[j];

                    tg_xxxy_yyyyyyyz_0[j] = pb_x * tg_xxy_yyyyyyyz_0[j] + wp_x[j] * tg_xxy_yyyyyyyz_1[j] + fl1_fx * tg_xy_yyyyyyyz_0[j] - fl1_fx * fl1_fza * tg_xy_yyyyyyyz_1[j];

                    tg_xxxy_yyyyyyzz_0[j] = pb_x * tg_xxy_yyyyyyzz_0[j] + wp_x[j] * tg_xxy_yyyyyyzz_1[j] + fl1_fx * tg_xy_yyyyyyzz_0[j] - fl1_fx * fl1_fza * tg_xy_yyyyyyzz_1[j];

                    tg_xxxy_yyyyyzzz_0[j] = pb_x * tg_xxy_yyyyyzzz_0[j] + wp_x[j] * tg_xxy_yyyyyzzz_1[j] + fl1_fx * tg_xy_yyyyyzzz_0[j] - fl1_fx * fl1_fza * tg_xy_yyyyyzzz_1[j];

                    tg_xxxy_yyyyzzzz_0[j] = pb_x * tg_xxy_yyyyzzzz_0[j] + wp_x[j] * tg_xxy_yyyyzzzz_1[j] + fl1_fx * tg_xy_yyyyzzzz_0[j] - fl1_fx * fl1_fza * tg_xy_yyyyzzzz_1[j];

                    tg_xxxy_yyyzzzzz_0[j] = pb_x * tg_xxy_yyyzzzzz_0[j] + wp_x[j] * tg_xxy_yyyzzzzz_1[j] + fl1_fx * tg_xy_yyyzzzzz_0[j] - fl1_fx * fl1_fza * tg_xy_yyyzzzzz_1[j];

                    tg_xxxy_yyzzzzzz_0[j] = pb_x * tg_xxy_yyzzzzzz_0[j] + wp_x[j] * tg_xxy_yyzzzzzz_1[j] + fl1_fx * tg_xy_yyzzzzzz_0[j] - fl1_fx * fl1_fza * tg_xy_yyzzzzzz_1[j];

                    tg_xxxy_yzzzzzzz_0[j] = pb_x * tg_xxy_yzzzzzzz_0[j] + wp_x[j] * tg_xxy_yzzzzzzz_1[j] + fl1_fx * tg_xy_yzzzzzzz_0[j] - fl1_fx * fl1_fza * tg_xy_yzzzzzzz_1[j];

                    tg_xxxy_zzzzzzzz_0[j] = pb_x * tg_xxy_zzzzzzzz_0[j] + wp_x[j] * tg_xxy_zzzzzzzz_1[j] + fl1_fx * tg_xy_zzzzzzzz_0[j] - fl1_fx * fl1_fza * tg_xy_zzzzzzzz_1[j];

                    tg_xxxz_xxxxxxxx_0[j] = pb_x * tg_xxz_xxxxxxxx_0[j] + wp_x[j] * tg_xxz_xxxxxxxx_1[j] + fl1_fx * tg_xz_xxxxxxxx_0[j] - fl1_fx * fl1_fza * tg_xz_xxxxxxxx_1[j] + 4.0 * fl1_fxn * tg_xxz_xxxxxxx_1[j];

                    tg_xxxz_xxxxxxxy_0[j] = pb_x * tg_xxz_xxxxxxxy_0[j] + wp_x[j] * tg_xxz_xxxxxxxy_1[j] + fl1_fx * tg_xz_xxxxxxxy_0[j] - fl1_fx * fl1_fza * tg_xz_xxxxxxxy_1[j] + 3.5 * fl1_fxn * tg_xxz_xxxxxxy_1[j];

                    tg_xxxz_xxxxxxxz_0[j] = pb_x * tg_xxz_xxxxxxxz_0[j] + wp_x[j] * tg_xxz_xxxxxxxz_1[j] + fl1_fx * tg_xz_xxxxxxxz_0[j] - fl1_fx * fl1_fza * tg_xz_xxxxxxxz_1[j] + 3.5 * fl1_fxn * tg_xxz_xxxxxxz_1[j];

                    tg_xxxz_xxxxxxyy_0[j] = pb_x * tg_xxz_xxxxxxyy_0[j] + wp_x[j] * tg_xxz_xxxxxxyy_1[j] + fl1_fx * tg_xz_xxxxxxyy_0[j] - fl1_fx * fl1_fza * tg_xz_xxxxxxyy_1[j] + 3.0 * fl1_fxn * tg_xxz_xxxxxyy_1[j];

                    tg_xxxz_xxxxxxyz_0[j] = pb_x * tg_xxz_xxxxxxyz_0[j] + wp_x[j] * tg_xxz_xxxxxxyz_1[j] + fl1_fx * tg_xz_xxxxxxyz_0[j] - fl1_fx * fl1_fza * tg_xz_xxxxxxyz_1[j] + 3.0 * fl1_fxn * tg_xxz_xxxxxyz_1[j];

                    tg_xxxz_xxxxxxzz_0[j] = pb_x * tg_xxz_xxxxxxzz_0[j] + wp_x[j] * tg_xxz_xxxxxxzz_1[j] + fl1_fx * tg_xz_xxxxxxzz_0[j] - fl1_fx * fl1_fza * tg_xz_xxxxxxzz_1[j] + 3.0 * fl1_fxn * tg_xxz_xxxxxzz_1[j];

                    tg_xxxz_xxxxxyyy_0[j] = pb_x * tg_xxz_xxxxxyyy_0[j] + wp_x[j] * tg_xxz_xxxxxyyy_1[j] + fl1_fx * tg_xz_xxxxxyyy_0[j] - fl1_fx * fl1_fza * tg_xz_xxxxxyyy_1[j] + 2.5 * fl1_fxn * tg_xxz_xxxxyyy_1[j];

                    tg_xxxz_xxxxxyyz_0[j] = pb_x * tg_xxz_xxxxxyyz_0[j] + wp_x[j] * tg_xxz_xxxxxyyz_1[j] + fl1_fx * tg_xz_xxxxxyyz_0[j] - fl1_fx * fl1_fza * tg_xz_xxxxxyyz_1[j] + 2.5 * fl1_fxn * tg_xxz_xxxxyyz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSL_98_147(      CMemBlock2D<double>& primBuffer,
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
                                             {4, -1, -1, -1},
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_2_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_2_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

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

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                // set up pointers to auxilary integrals

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

                auto tg_xxz_xxxxxyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 98); 

                auto tg_xxz_xxxxxzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 99); 

                auto tg_xxz_xxxxyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 100); 

                auto tg_xxz_xxxxyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 101); 

                auto tg_xxz_xxxxyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 102); 

                auto tg_xxz_xxxxyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 103); 

                auto tg_xxz_xxxxzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 104); 

                auto tg_xxz_xxxyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 105); 

                auto tg_xxz_xxxyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 106); 

                auto tg_xxz_xxxyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 107); 

                auto tg_xxz_xxxyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 108); 

                auto tg_xxz_xxxyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 109); 

                auto tg_xxz_xxxzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 110); 

                auto tg_xxz_xxyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 111); 

                auto tg_xxz_xxyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 112); 

                auto tg_xxz_xxyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 113); 

                auto tg_xxz_xxyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 114); 

                auto tg_xxz_xxyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 115); 

                auto tg_xxz_xxyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 116); 

                auto tg_xxz_xxzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 117); 

                auto tg_xxz_xyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 118); 

                auto tg_xxz_xyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 119); 

                auto tg_xxz_xyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 120); 

                auto tg_xxz_xyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 121); 

                auto tg_xxz_xyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 122); 

                auto tg_xxz_xyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 123); 

                auto tg_xxz_xyzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 124); 

                auto tg_xxz_xzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 125); 

                auto tg_xxz_yyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 126); 

                auto tg_xxz_yyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 127); 

                auto tg_xxz_yyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 128); 

                auto tg_xxz_yyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 129); 

                auto tg_xxz_yyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 130); 

                auto tg_xxz_yyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 131); 

                auto tg_xxz_yyzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 132); 

                auto tg_xxz_yzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 133); 

                auto tg_xxz_zzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 134); 

                auto tg_xyy_xxxxxxxx_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 135); 

                auto tg_xyy_xxxxxxxy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 136); 

                auto tg_xyy_xxxxxxxz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 137); 

                auto tg_xyy_xxxxxxyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 138); 

                auto tg_xyy_xxxxxxyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 139); 

                auto tg_xyy_xxxxxxzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 140); 

                auto tg_xyy_xxxxxyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 141); 

                auto tg_xyy_xxxxxyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 142); 

                auto tg_xyy_xxxxxyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 143); 

                auto tg_xyy_xxxxxzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 144); 

                auto tg_xyy_xxxxyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 145); 

                auto tg_xyy_xxxxyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 146); 

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

                // set up pointers to integrals

                auto tg_xxxz_xxxxxyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 98); 

                auto tg_xxxz_xxxxxzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 99); 

                auto tg_xxxz_xxxxyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 100); 

                auto tg_xxxz_xxxxyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 101); 

                auto tg_xxxz_xxxxyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 102); 

                auto tg_xxxz_xxxxyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 103); 

                auto tg_xxxz_xxxxzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 104); 

                auto tg_xxxz_xxxyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 105); 

                auto tg_xxxz_xxxyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 106); 

                auto tg_xxxz_xxxyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 107); 

                auto tg_xxxz_xxxyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 108); 

                auto tg_xxxz_xxxyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 109); 

                auto tg_xxxz_xxxzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 110); 

                auto tg_xxxz_xxyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 111); 

                auto tg_xxxz_xxyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 112); 

                auto tg_xxxz_xxyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 113); 

                auto tg_xxxz_xxyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 114); 

                auto tg_xxxz_xxyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 115); 

                auto tg_xxxz_xxyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 116); 

                auto tg_xxxz_xxzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 117); 

                auto tg_xxxz_xyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 118); 

                auto tg_xxxz_xyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 119); 

                auto tg_xxxz_xyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 120); 

                auto tg_xxxz_xyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 121); 

                auto tg_xxxz_xyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 122); 

                auto tg_xxxz_xyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 123); 

                auto tg_xxxz_xyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 124); 

                auto tg_xxxz_xzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 125); 

                auto tg_xxxz_yyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 126); 

                auto tg_xxxz_yyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 127); 

                auto tg_xxxz_yyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 128); 

                auto tg_xxxz_yyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 129); 

                auto tg_xxxz_yyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 130); 

                auto tg_xxxz_yyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 131); 

                auto tg_xxxz_yyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 132); 

                auto tg_xxxz_yzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 133); 

                auto tg_xxxz_zzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 134); 

                auto tg_xxyy_xxxxxxxx_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 135); 

                auto tg_xxyy_xxxxxxxy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 136); 

                auto tg_xxyy_xxxxxxxz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 137); 

                auto tg_xxyy_xxxxxxyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 138); 

                auto tg_xxyy_xxxxxxyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 139); 

                auto tg_xxyy_xxxxxxzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 140); 

                auto tg_xxyy_xxxxxyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 141); 

                auto tg_xxyy_xxxxxyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 142); 

                auto tg_xxyy_xxxxxyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 143); 

                auto tg_xxyy_xxxxxzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 144); 

                auto tg_xxyy_xxxxyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 145); 

                auto tg_xxyy_xxxxyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 146); 

                // Batch of Integrals (98,147)

                #pragma omp simd aligned(fxn, fza, tg_xxxz_xxxxxyzz_0, tg_xxxz_xxxxxzzz_0, tg_xxxz_xxxxyyyy_0, \
                                         tg_xxxz_xxxxyyyz_0, tg_xxxz_xxxxyyzz_0, tg_xxxz_xxxxyzzz_0, tg_xxxz_xxxxzzzz_0, \
                                         tg_xxxz_xxxyyyyy_0, tg_xxxz_xxxyyyyz_0, tg_xxxz_xxxyyyzz_0, tg_xxxz_xxxyyzzz_0, \
                                         tg_xxxz_xxxyzzzz_0, tg_xxxz_xxxzzzzz_0, tg_xxxz_xxyyyyyy_0, tg_xxxz_xxyyyyyz_0, \
                                         tg_xxxz_xxyyyyzz_0, tg_xxxz_xxyyyzzz_0, tg_xxxz_xxyyzzzz_0, tg_xxxz_xxyzzzzz_0, \
                                         tg_xxxz_xxzzzzzz_0, tg_xxxz_xyyyyyyy_0, tg_xxxz_xyyyyyyz_0, tg_xxxz_xyyyyyzz_0, \
                                         tg_xxxz_xyyyyzzz_0, tg_xxxz_xyyyzzzz_0, tg_xxxz_xyyzzzzz_0, tg_xxxz_xyzzzzzz_0, \
                                         tg_xxxz_xzzzzzzz_0, tg_xxxz_yyyyyyyy_0, tg_xxxz_yyyyyyyz_0, tg_xxxz_yyyyyyzz_0, \
                                         tg_xxxz_yyyyyzzz_0, tg_xxxz_yyyyzzzz_0, tg_xxxz_yyyzzzzz_0, tg_xxxz_yyzzzzzz_0, \
                                         tg_xxxz_yzzzzzzz_0, tg_xxxz_zzzzzzzz_0, tg_xxyy_xxxxxxxx_0, tg_xxyy_xxxxxxxy_0, \
                                         tg_xxyy_xxxxxxxz_0, tg_xxyy_xxxxxxyy_0, tg_xxyy_xxxxxxyz_0, tg_xxyy_xxxxxxzz_0, \
                                         tg_xxyy_xxxxxyyy_0, tg_xxyy_xxxxxyyz_0, tg_xxyy_xxxxxyzz_0, tg_xxyy_xxxxxzzz_0, \
                                         tg_xxyy_xxxxyyyy_0, tg_xxyy_xxxxyyyz_0, tg_xxz_xxxxxyzz_0, tg_xxz_xxxxxyzz_1, \
                                         tg_xxz_xxxxxzzz_0, tg_xxz_xxxxxzzz_1, tg_xxz_xxxxyyyy_0, tg_xxz_xxxxyyyy_1, \
                                         tg_xxz_xxxxyyyz_0, tg_xxz_xxxxyyyz_1, tg_xxz_xxxxyyzz_0, tg_xxz_xxxxyyzz_1, \
                                         tg_xxz_xxxxyzz_1, tg_xxz_xxxxyzzz_0, tg_xxz_xxxxyzzz_1, tg_xxz_xxxxzzz_1, \
                                         tg_xxz_xxxxzzzz_0, tg_xxz_xxxxzzzz_1, tg_xxz_xxxyyyy_1, tg_xxz_xxxyyyyy_0, \
                                         tg_xxz_xxxyyyyy_1, tg_xxz_xxxyyyyz_0, tg_xxz_xxxyyyyz_1, tg_xxz_xxxyyyz_1, \
                                         tg_xxz_xxxyyyzz_0, tg_xxz_xxxyyyzz_1, tg_xxz_xxxyyzz_1, tg_xxz_xxxyyzzz_0, \
                                         tg_xxz_xxxyyzzz_1, tg_xxz_xxxyzzz_1, tg_xxz_xxxyzzzz_0, tg_xxz_xxxyzzzz_1, \
                                         tg_xxz_xxxzzzz_1, tg_xxz_xxxzzzzz_0, tg_xxz_xxxzzzzz_1, tg_xxz_xxyyyyy_1, \
                                         tg_xxz_xxyyyyyy_0, tg_xxz_xxyyyyyy_1, tg_xxz_xxyyyyyz_0, tg_xxz_xxyyyyyz_1, \
                                         tg_xxz_xxyyyyz_1, tg_xxz_xxyyyyzz_0, tg_xxz_xxyyyyzz_1, tg_xxz_xxyyyzz_1, \
                                         tg_xxz_xxyyyzzz_0, tg_xxz_xxyyyzzz_1, tg_xxz_xxyyzzz_1, tg_xxz_xxyyzzzz_0, \
                                         tg_xxz_xxyyzzzz_1, tg_xxz_xxyzzzz_1, tg_xxz_xxyzzzzz_0, tg_xxz_xxyzzzzz_1, \
                                         tg_xxz_xxzzzzz_1, tg_xxz_xxzzzzzz_0, tg_xxz_xxzzzzzz_1, tg_xxz_xyyyyyy_1, \
                                         tg_xxz_xyyyyyyy_0, tg_xxz_xyyyyyyy_1, tg_xxz_xyyyyyyz_0, tg_xxz_xyyyyyyz_1, \
                                         tg_xxz_xyyyyyz_1, tg_xxz_xyyyyyzz_0, tg_xxz_xyyyyyzz_1, tg_xxz_xyyyyzz_1, \
                                         tg_xxz_xyyyyzzz_0, tg_xxz_xyyyyzzz_1, tg_xxz_xyyyzzz_1, tg_xxz_xyyyzzzz_0, \
                                         tg_xxz_xyyyzzzz_1, tg_xxz_xyyzzzz_1, tg_xxz_xyyzzzzz_0, tg_xxz_xyyzzzzz_1, \
                                         tg_xxz_xyzzzzz_1, tg_xxz_xyzzzzzz_0, tg_xxz_xyzzzzzz_1, tg_xxz_xzzzzzz_1, \
                                         tg_xxz_xzzzzzzz_0, tg_xxz_xzzzzzzz_1, tg_xxz_yyyyyyy_1, tg_xxz_yyyyyyyy_0, \
                                         tg_xxz_yyyyyyyy_1, tg_xxz_yyyyyyyz_0, tg_xxz_yyyyyyyz_1, tg_xxz_yyyyyyz_1, \
                                         tg_xxz_yyyyyyzz_0, tg_xxz_yyyyyyzz_1, tg_xxz_yyyyyzz_1, tg_xxz_yyyyyzzz_0, \
                                         tg_xxz_yyyyyzzz_1, tg_xxz_yyyyzzz_1, tg_xxz_yyyyzzzz_0, tg_xxz_yyyyzzzz_1, \
                                         tg_xxz_yyyzzzz_1, tg_xxz_yyyzzzzz_0, tg_xxz_yyyzzzzz_1, tg_xxz_yyzzzzz_1, \
                                         tg_xxz_yyzzzzzz_0, tg_xxz_yyzzzzzz_1, tg_xxz_yzzzzzz_1, tg_xxz_yzzzzzzz_0, \
                                         tg_xxz_yzzzzzzz_1, tg_xxz_zzzzzzz_1, tg_xxz_zzzzzzzz_0, tg_xxz_zzzzzzzz_1, \
                                         tg_xyy_xxxxxxx_1, tg_xyy_xxxxxxxx_0, tg_xyy_xxxxxxxx_1, tg_xyy_xxxxxxxy_0, \
                                         tg_xyy_xxxxxxxy_1, tg_xyy_xxxxxxxz_0, tg_xyy_xxxxxxxz_1, tg_xyy_xxxxxxy_1, \
                                         tg_xyy_xxxxxxyy_0, tg_xyy_xxxxxxyy_1, tg_xyy_xxxxxxyz_0, tg_xyy_xxxxxxyz_1, \
                                         tg_xyy_xxxxxxz_1, tg_xyy_xxxxxxzz_0, tg_xyy_xxxxxxzz_1, tg_xyy_xxxxxyy_1, \
                                         tg_xyy_xxxxxyyy_0, tg_xyy_xxxxxyyy_1, tg_xyy_xxxxxyyz_0, tg_xyy_xxxxxyyz_1, \
                                         tg_xyy_xxxxxyz_1, tg_xyy_xxxxxyzz_0, tg_xyy_xxxxxyzz_1, tg_xyy_xxxxxzz_1, \
                                         tg_xyy_xxxxxzzz_0, tg_xyy_xxxxxzzz_1, tg_xyy_xxxxyyy_1, tg_xyy_xxxxyyyy_0, \
                                         tg_xyy_xxxxyyyy_1, tg_xyy_xxxxyyyz_0, tg_xyy_xxxxyyyz_1, tg_xyy_xxxxyyz_1, \
                                         tg_xyy_xxxxyzz_1, tg_xyy_xxxxzzz_1, tg_xyy_xxxyyyy_1, tg_xyy_xxxyyyz_1, \
                                         tg_xz_xxxxxyzz_0, tg_xz_xxxxxyzz_1, tg_xz_xxxxxzzz_0, tg_xz_xxxxxzzz_1, \
                                         tg_xz_xxxxyyyy_0, tg_xz_xxxxyyyy_1, tg_xz_xxxxyyyz_0, tg_xz_xxxxyyyz_1, \
                                         tg_xz_xxxxyyzz_0, tg_xz_xxxxyyzz_1, tg_xz_xxxxyzzz_0, tg_xz_xxxxyzzz_1, \
                                         tg_xz_xxxxzzzz_0, tg_xz_xxxxzzzz_1, tg_xz_xxxyyyyy_0, tg_xz_xxxyyyyy_1, \
                                         tg_xz_xxxyyyyz_0, tg_xz_xxxyyyyz_1, tg_xz_xxxyyyzz_0, tg_xz_xxxyyyzz_1, \
                                         tg_xz_xxxyyzzz_0, tg_xz_xxxyyzzz_1, tg_xz_xxxyzzzz_0, tg_xz_xxxyzzzz_1, \
                                         tg_xz_xxxzzzzz_0, tg_xz_xxxzzzzz_1, tg_xz_xxyyyyyy_0, tg_xz_xxyyyyyy_1, \
                                         tg_xz_xxyyyyyz_0, tg_xz_xxyyyyyz_1, tg_xz_xxyyyyzz_0, tg_xz_xxyyyyzz_1, \
                                         tg_xz_xxyyyzzz_0, tg_xz_xxyyyzzz_1, tg_xz_xxyyzzzz_0, tg_xz_xxyyzzzz_1, \
                                         tg_xz_xxyzzzzz_0, tg_xz_xxyzzzzz_1, tg_xz_xxzzzzzz_0, tg_xz_xxzzzzzz_1, \
                                         tg_xz_xyyyyyyy_0, tg_xz_xyyyyyyy_1, tg_xz_xyyyyyyz_0, tg_xz_xyyyyyyz_1, \
                                         tg_xz_xyyyyyzz_0, tg_xz_xyyyyyzz_1, tg_xz_xyyyyzzz_0, tg_xz_xyyyyzzz_1, \
                                         tg_xz_xyyyzzzz_0, tg_xz_xyyyzzzz_1, tg_xz_xyyzzzzz_0, tg_xz_xyyzzzzz_1, \
                                         tg_xz_xyzzzzzz_0, tg_xz_xyzzzzzz_1, tg_xz_xzzzzzzz_0, tg_xz_xzzzzzzz_1, \
                                         tg_xz_yyyyyyyy_0, tg_xz_yyyyyyyy_1, tg_xz_yyyyyyyz_0, tg_xz_yyyyyyyz_1, \
                                         tg_xz_yyyyyyzz_0, tg_xz_yyyyyyzz_1, tg_xz_yyyyyzzz_0, tg_xz_yyyyyzzz_1, \
                                         tg_xz_yyyyzzzz_0, tg_xz_yyyyzzzz_1, tg_xz_yyyzzzzz_0, tg_xz_yyyzzzzz_1, \
                                         tg_xz_yyzzzzzz_0, tg_xz_yyzzzzzz_1, tg_xz_yzzzzzzz_0, tg_xz_yzzzzzzz_1, \
                                         tg_xz_zzzzzzzz_0, tg_xz_zzzzzzzz_1, tg_yy_xxxxxxxx_0, tg_yy_xxxxxxxx_1, \
                                         tg_yy_xxxxxxxy_0, tg_yy_xxxxxxxy_1, tg_yy_xxxxxxxz_0, tg_yy_xxxxxxxz_1, \
                                         tg_yy_xxxxxxyy_0, tg_yy_xxxxxxyy_1, tg_yy_xxxxxxyz_0, tg_yy_xxxxxxyz_1, \
                                         tg_yy_xxxxxxzz_0, tg_yy_xxxxxxzz_1, tg_yy_xxxxxyyy_0, tg_yy_xxxxxyyy_1, \
                                         tg_yy_xxxxxyyz_0, tg_yy_xxxxxyyz_1, tg_yy_xxxxxyzz_0, tg_yy_xxxxxyzz_1, \
                                         tg_yy_xxxxxzzz_0, tg_yy_xxxxxzzz_1, tg_yy_xxxxyyyy_0, tg_yy_xxxxyyyy_1, \
                                         tg_yy_xxxxyyyz_0, tg_yy_xxxxyyyz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxxz_xxxxxyzz_0[j] = pb_x * tg_xxz_xxxxxyzz_0[j] + wp_x[j] * tg_xxz_xxxxxyzz_1[j] + fl1_fx * tg_xz_xxxxxyzz_0[j] - fl1_fx * fl1_fza * tg_xz_xxxxxyzz_1[j] + 2.5 * fl1_fxn * tg_xxz_xxxxyzz_1[j];

                    tg_xxxz_xxxxxzzz_0[j] = pb_x * tg_xxz_xxxxxzzz_0[j] + wp_x[j] * tg_xxz_xxxxxzzz_1[j] + fl1_fx * tg_xz_xxxxxzzz_0[j] - fl1_fx * fl1_fza * tg_xz_xxxxxzzz_1[j] + 2.5 * fl1_fxn * tg_xxz_xxxxzzz_1[j];

                    tg_xxxz_xxxxyyyy_0[j] = pb_x * tg_xxz_xxxxyyyy_0[j] + wp_x[j] * tg_xxz_xxxxyyyy_1[j] + fl1_fx * tg_xz_xxxxyyyy_0[j] - fl1_fx * fl1_fza * tg_xz_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_xxz_xxxyyyy_1[j];

                    tg_xxxz_xxxxyyyz_0[j] = pb_x * tg_xxz_xxxxyyyz_0[j] + wp_x[j] * tg_xxz_xxxxyyyz_1[j] + fl1_fx * tg_xz_xxxxyyyz_0[j] - fl1_fx * fl1_fza * tg_xz_xxxxyyyz_1[j] + 2.0 * fl1_fxn * tg_xxz_xxxyyyz_1[j];

                    tg_xxxz_xxxxyyzz_0[j] = pb_x * tg_xxz_xxxxyyzz_0[j] + wp_x[j] * tg_xxz_xxxxyyzz_1[j] + fl1_fx * tg_xz_xxxxyyzz_0[j] - fl1_fx * fl1_fza * tg_xz_xxxxyyzz_1[j] + 2.0 * fl1_fxn * tg_xxz_xxxyyzz_1[j];

                    tg_xxxz_xxxxyzzz_0[j] = pb_x * tg_xxz_xxxxyzzz_0[j] + wp_x[j] * tg_xxz_xxxxyzzz_1[j] + fl1_fx * tg_xz_xxxxyzzz_0[j] - fl1_fx * fl1_fza * tg_xz_xxxxyzzz_1[j] + 2.0 * fl1_fxn * tg_xxz_xxxyzzz_1[j];

                    tg_xxxz_xxxxzzzz_0[j] = pb_x * tg_xxz_xxxxzzzz_0[j] + wp_x[j] * tg_xxz_xxxxzzzz_1[j] + fl1_fx * tg_xz_xxxxzzzz_0[j] - fl1_fx * fl1_fza * tg_xz_xxxxzzzz_1[j] + 2.0 * fl1_fxn * tg_xxz_xxxzzzz_1[j];

                    tg_xxxz_xxxyyyyy_0[j] = pb_x * tg_xxz_xxxyyyyy_0[j] + wp_x[j] * tg_xxz_xxxyyyyy_1[j] + fl1_fx * tg_xz_xxxyyyyy_0[j] - fl1_fx * fl1_fza * tg_xz_xxxyyyyy_1[j] + 1.5 * fl1_fxn * tg_xxz_xxyyyyy_1[j];

                    tg_xxxz_xxxyyyyz_0[j] = pb_x * tg_xxz_xxxyyyyz_0[j] + wp_x[j] * tg_xxz_xxxyyyyz_1[j] + fl1_fx * tg_xz_xxxyyyyz_0[j] - fl1_fx * fl1_fza * tg_xz_xxxyyyyz_1[j] + 1.5 * fl1_fxn * tg_xxz_xxyyyyz_1[j];

                    tg_xxxz_xxxyyyzz_0[j] = pb_x * tg_xxz_xxxyyyzz_0[j] + wp_x[j] * tg_xxz_xxxyyyzz_1[j] + fl1_fx * tg_xz_xxxyyyzz_0[j] - fl1_fx * fl1_fza * tg_xz_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_xxz_xxyyyzz_1[j];

                    tg_xxxz_xxxyyzzz_0[j] = pb_x * tg_xxz_xxxyyzzz_0[j] + wp_x[j] * tg_xxz_xxxyyzzz_1[j] + fl1_fx * tg_xz_xxxyyzzz_0[j] - fl1_fx * fl1_fza * tg_xz_xxxyyzzz_1[j] + 1.5 * fl1_fxn * tg_xxz_xxyyzzz_1[j];

                    tg_xxxz_xxxyzzzz_0[j] = pb_x * tg_xxz_xxxyzzzz_0[j] + wp_x[j] * tg_xxz_xxxyzzzz_1[j] + fl1_fx * tg_xz_xxxyzzzz_0[j] - fl1_fx * fl1_fza * tg_xz_xxxyzzzz_1[j] + 1.5 * fl1_fxn * tg_xxz_xxyzzzz_1[j];

                    tg_xxxz_xxxzzzzz_0[j] = pb_x * tg_xxz_xxxzzzzz_0[j] + wp_x[j] * tg_xxz_xxxzzzzz_1[j] + fl1_fx * tg_xz_xxxzzzzz_0[j] - fl1_fx * fl1_fza * tg_xz_xxxzzzzz_1[j] + 1.5 * fl1_fxn * tg_xxz_xxzzzzz_1[j];

                    tg_xxxz_xxyyyyyy_0[j] = pb_x * tg_xxz_xxyyyyyy_0[j] + wp_x[j] * tg_xxz_xxyyyyyy_1[j] + fl1_fx * tg_xz_xxyyyyyy_0[j] - fl1_fx * fl1_fza * tg_xz_xxyyyyyy_1[j] + fl1_fxn * tg_xxz_xyyyyyy_1[j];

                    tg_xxxz_xxyyyyyz_0[j] = pb_x * tg_xxz_xxyyyyyz_0[j] + wp_x[j] * tg_xxz_xxyyyyyz_1[j] + fl1_fx * tg_xz_xxyyyyyz_0[j] - fl1_fx * fl1_fza * tg_xz_xxyyyyyz_1[j] + fl1_fxn * tg_xxz_xyyyyyz_1[j];

                    tg_xxxz_xxyyyyzz_0[j] = pb_x * tg_xxz_xxyyyyzz_0[j] + wp_x[j] * tg_xxz_xxyyyyzz_1[j] + fl1_fx * tg_xz_xxyyyyzz_0[j] - fl1_fx * fl1_fza * tg_xz_xxyyyyzz_1[j] + fl1_fxn * tg_xxz_xyyyyzz_1[j];

                    tg_xxxz_xxyyyzzz_0[j] = pb_x * tg_xxz_xxyyyzzz_0[j] + wp_x[j] * tg_xxz_xxyyyzzz_1[j] + fl1_fx * tg_xz_xxyyyzzz_0[j] - fl1_fx * fl1_fza * tg_xz_xxyyyzzz_1[j] + fl1_fxn * tg_xxz_xyyyzzz_1[j];

                    tg_xxxz_xxyyzzzz_0[j] = pb_x * tg_xxz_xxyyzzzz_0[j] + wp_x[j] * tg_xxz_xxyyzzzz_1[j] + fl1_fx * tg_xz_xxyyzzzz_0[j] - fl1_fx * fl1_fza * tg_xz_xxyyzzzz_1[j] + fl1_fxn * tg_xxz_xyyzzzz_1[j];

                    tg_xxxz_xxyzzzzz_0[j] = pb_x * tg_xxz_xxyzzzzz_0[j] + wp_x[j] * tg_xxz_xxyzzzzz_1[j] + fl1_fx * tg_xz_xxyzzzzz_0[j] - fl1_fx * fl1_fza * tg_xz_xxyzzzzz_1[j] + fl1_fxn * tg_xxz_xyzzzzz_1[j];

                    tg_xxxz_xxzzzzzz_0[j] = pb_x * tg_xxz_xxzzzzzz_0[j] + wp_x[j] * tg_xxz_xxzzzzzz_1[j] + fl1_fx * tg_xz_xxzzzzzz_0[j] - fl1_fx * fl1_fza * tg_xz_xxzzzzzz_1[j] + fl1_fxn * tg_xxz_xzzzzzz_1[j];

                    tg_xxxz_xyyyyyyy_0[j] = pb_x * tg_xxz_xyyyyyyy_0[j] + wp_x[j] * tg_xxz_xyyyyyyy_1[j] + fl1_fx * tg_xz_xyyyyyyy_0[j] - fl1_fx * fl1_fza * tg_xz_xyyyyyyy_1[j] + 0.5 * fl1_fxn * tg_xxz_yyyyyyy_1[j];

                    tg_xxxz_xyyyyyyz_0[j] = pb_x * tg_xxz_xyyyyyyz_0[j] + wp_x[j] * tg_xxz_xyyyyyyz_1[j] + fl1_fx * tg_xz_xyyyyyyz_0[j] - fl1_fx * fl1_fza * tg_xz_xyyyyyyz_1[j] + 0.5 * fl1_fxn * tg_xxz_yyyyyyz_1[j];

                    tg_xxxz_xyyyyyzz_0[j] = pb_x * tg_xxz_xyyyyyzz_0[j] + wp_x[j] * tg_xxz_xyyyyyzz_1[j] + fl1_fx * tg_xz_xyyyyyzz_0[j] - fl1_fx * fl1_fza * tg_xz_xyyyyyzz_1[j] + 0.5 * fl1_fxn * tg_xxz_yyyyyzz_1[j];

                    tg_xxxz_xyyyyzzz_0[j] = pb_x * tg_xxz_xyyyyzzz_0[j] + wp_x[j] * tg_xxz_xyyyyzzz_1[j] + fl1_fx * tg_xz_xyyyyzzz_0[j] - fl1_fx * fl1_fza * tg_xz_xyyyyzzz_1[j] + 0.5 * fl1_fxn * tg_xxz_yyyyzzz_1[j];

                    tg_xxxz_xyyyzzzz_0[j] = pb_x * tg_xxz_xyyyzzzz_0[j] + wp_x[j] * tg_xxz_xyyyzzzz_1[j] + fl1_fx * tg_xz_xyyyzzzz_0[j] - fl1_fx * fl1_fza * tg_xz_xyyyzzzz_1[j] + 0.5 * fl1_fxn * tg_xxz_yyyzzzz_1[j];

                    tg_xxxz_xyyzzzzz_0[j] = pb_x * tg_xxz_xyyzzzzz_0[j] + wp_x[j] * tg_xxz_xyyzzzzz_1[j] + fl1_fx * tg_xz_xyyzzzzz_0[j] - fl1_fx * fl1_fza * tg_xz_xyyzzzzz_1[j] + 0.5 * fl1_fxn * tg_xxz_yyzzzzz_1[j];

                    tg_xxxz_xyzzzzzz_0[j] = pb_x * tg_xxz_xyzzzzzz_0[j] + wp_x[j] * tg_xxz_xyzzzzzz_1[j] + fl1_fx * tg_xz_xyzzzzzz_0[j] - fl1_fx * fl1_fza * tg_xz_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xxz_yzzzzzz_1[j];

                    tg_xxxz_xzzzzzzz_0[j] = pb_x * tg_xxz_xzzzzzzz_0[j] + wp_x[j] * tg_xxz_xzzzzzzz_1[j] + fl1_fx * tg_xz_xzzzzzzz_0[j] - fl1_fx * fl1_fza * tg_xz_xzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xxz_zzzzzzz_1[j];

                    tg_xxxz_yyyyyyyy_0[j] = pb_x * tg_xxz_yyyyyyyy_0[j] + wp_x[j] * tg_xxz_yyyyyyyy_1[j] + fl1_fx * tg_xz_yyyyyyyy_0[j] - fl1_fx * fl1_fza * tg_xz_yyyyyyyy_1[j];

                    tg_xxxz_yyyyyyyz_0[j] = pb_x * tg_xxz_yyyyyyyz_0[j] + wp_x[j] * tg_xxz_yyyyyyyz_1[j] + fl1_fx * tg_xz_yyyyyyyz_0[j] - fl1_fx * fl1_fza * tg_xz_yyyyyyyz_1[j];

                    tg_xxxz_yyyyyyzz_0[j] = pb_x * tg_xxz_yyyyyyzz_0[j] + wp_x[j] * tg_xxz_yyyyyyzz_1[j] + fl1_fx * tg_xz_yyyyyyzz_0[j] - fl1_fx * fl1_fza * tg_xz_yyyyyyzz_1[j];

                    tg_xxxz_yyyyyzzz_0[j] = pb_x * tg_xxz_yyyyyzzz_0[j] + wp_x[j] * tg_xxz_yyyyyzzz_1[j] + fl1_fx * tg_xz_yyyyyzzz_0[j] - fl1_fx * fl1_fza * tg_xz_yyyyyzzz_1[j];

                    tg_xxxz_yyyyzzzz_0[j] = pb_x * tg_xxz_yyyyzzzz_0[j] + wp_x[j] * tg_xxz_yyyyzzzz_1[j] + fl1_fx * tg_xz_yyyyzzzz_0[j] - fl1_fx * fl1_fza * tg_xz_yyyyzzzz_1[j];

                    tg_xxxz_yyyzzzzz_0[j] = pb_x * tg_xxz_yyyzzzzz_0[j] + wp_x[j] * tg_xxz_yyyzzzzz_1[j] + fl1_fx * tg_xz_yyyzzzzz_0[j] - fl1_fx * fl1_fza * tg_xz_yyyzzzzz_1[j];

                    tg_xxxz_yyzzzzzz_0[j] = pb_x * tg_xxz_yyzzzzzz_0[j] + wp_x[j] * tg_xxz_yyzzzzzz_1[j] + fl1_fx * tg_xz_yyzzzzzz_0[j] - fl1_fx * fl1_fza * tg_xz_yyzzzzzz_1[j];

                    tg_xxxz_yzzzzzzz_0[j] = pb_x * tg_xxz_yzzzzzzz_0[j] + wp_x[j] * tg_xxz_yzzzzzzz_1[j] + fl1_fx * tg_xz_yzzzzzzz_0[j] - fl1_fx * fl1_fza * tg_xz_yzzzzzzz_1[j];

                    tg_xxxz_zzzzzzzz_0[j] = pb_x * tg_xxz_zzzzzzzz_0[j] + wp_x[j] * tg_xxz_zzzzzzzz_1[j] + fl1_fx * tg_xz_zzzzzzzz_0[j] - fl1_fx * fl1_fza * tg_xz_zzzzzzzz_1[j];

                    tg_xxyy_xxxxxxxx_0[j] = pb_x * tg_xyy_xxxxxxxx_0[j] + wp_x[j] * tg_xyy_xxxxxxxx_1[j] + 0.5 * fl1_fx * tg_yy_xxxxxxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxxxxxxx_1[j] + 4.0 * fl1_fxn * tg_xyy_xxxxxxx_1[j];

                    tg_xxyy_xxxxxxxy_0[j] = pb_x * tg_xyy_xxxxxxxy_0[j] + wp_x[j] * tg_xyy_xxxxxxxy_1[j] + 0.5 * fl1_fx * tg_yy_xxxxxxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxxxxxxy_1[j] + 3.5 * fl1_fxn * tg_xyy_xxxxxxy_1[j];

                    tg_xxyy_xxxxxxxz_0[j] = pb_x * tg_xyy_xxxxxxxz_0[j] + wp_x[j] * tg_xyy_xxxxxxxz_1[j] + 0.5 * fl1_fx * tg_yy_xxxxxxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxxxxxxz_1[j] + 3.5 * fl1_fxn * tg_xyy_xxxxxxz_1[j];

                    tg_xxyy_xxxxxxyy_0[j] = pb_x * tg_xyy_xxxxxxyy_0[j] + wp_x[j] * tg_xyy_xxxxxxyy_1[j] + 0.5 * fl1_fx * tg_yy_xxxxxxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxxxxxyy_1[j] + 3.0 * fl1_fxn * tg_xyy_xxxxxyy_1[j];

                    tg_xxyy_xxxxxxyz_0[j] = pb_x * tg_xyy_xxxxxxyz_0[j] + wp_x[j] * tg_xyy_xxxxxxyz_1[j] + 0.5 * fl1_fx * tg_yy_xxxxxxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxxxxxyz_1[j] + 3.0 * fl1_fxn * tg_xyy_xxxxxyz_1[j];

                    tg_xxyy_xxxxxxzz_0[j] = pb_x * tg_xyy_xxxxxxzz_0[j] + wp_x[j] * tg_xyy_xxxxxxzz_1[j] + 0.5 * fl1_fx * tg_yy_xxxxxxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxxxxxzz_1[j] + 3.0 * fl1_fxn * tg_xyy_xxxxxzz_1[j];

                    tg_xxyy_xxxxxyyy_0[j] = pb_x * tg_xyy_xxxxxyyy_0[j] + wp_x[j] * tg_xyy_xxxxxyyy_1[j] + 0.5 * fl1_fx * tg_yy_xxxxxyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxxxxyyy_1[j] + 2.5 * fl1_fxn * tg_xyy_xxxxyyy_1[j];

                    tg_xxyy_xxxxxyyz_0[j] = pb_x * tg_xyy_xxxxxyyz_0[j] + wp_x[j] * tg_xyy_xxxxxyyz_1[j] + 0.5 * fl1_fx * tg_yy_xxxxxyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxxxxyyz_1[j] + 2.5 * fl1_fxn * tg_xyy_xxxxyyz_1[j];

                    tg_xxyy_xxxxxyzz_0[j] = pb_x * tg_xyy_xxxxxyzz_0[j] + wp_x[j] * tg_xyy_xxxxxyzz_1[j] + 0.5 * fl1_fx * tg_yy_xxxxxyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxxxxyzz_1[j] + 2.5 * fl1_fxn * tg_xyy_xxxxyzz_1[j];

                    tg_xxyy_xxxxxzzz_0[j] = pb_x * tg_xyy_xxxxxzzz_0[j] + wp_x[j] * tg_xyy_xxxxxzzz_1[j] + 0.5 * fl1_fx * tg_yy_xxxxxzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxxxxzzz_1[j] + 2.5 * fl1_fxn * tg_xyy_xxxxzzz_1[j];

                    tg_xxyy_xxxxyyyy_0[j] = pb_x * tg_xyy_xxxxyyyy_0[j] + wp_x[j] * tg_xyy_xxxxyyyy_1[j] + 0.5 * fl1_fx * tg_yy_xxxxyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_xyy_xxxyyyy_1[j];

                    tg_xxyy_xxxxyyyz_0[j] = pb_x * tg_xyy_xxxxyyyz_0[j] + wp_x[j] * tg_xyy_xxxxyyyz_1[j] + 0.5 * fl1_fx * tg_yy_xxxxyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxxxyyyz_1[j] + 2.0 * fl1_fxn * tg_xyy_xxxyyyz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSL_147_195(      CMemBlock2D<double>& primBuffer,
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
                                             {4, -1, -1, -1},
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_2_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_2_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

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

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                // set up pointers to auxilary integrals

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

                auto tg_xyy_xxxxyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 147); 

                auto tg_xyy_xxxxyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 148); 

                auto tg_xyy_xxxxzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 149); 

                auto tg_xyy_xxxyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 150); 

                auto tg_xyy_xxxyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 151); 

                auto tg_xyy_xxxyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 152); 

                auto tg_xyy_xxxyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 153); 

                auto tg_xyy_xxxyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 154); 

                auto tg_xyy_xxxzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 155); 

                auto tg_xyy_xxyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 156); 

                auto tg_xyy_xxyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 157); 

                auto tg_xyy_xxyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 158); 

                auto tg_xyy_xxyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 159); 

                auto tg_xyy_xxyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 160); 

                auto tg_xyy_xxyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 161); 

                auto tg_xyy_xxzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 162); 

                auto tg_xyy_xyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 163); 

                auto tg_xyy_xyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 164); 

                auto tg_xyy_xyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 165); 

                auto tg_xyy_xyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 166); 

                auto tg_xyy_xyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 167); 

                auto tg_xyy_xyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 168); 

                auto tg_xyy_xyzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 169); 

                auto tg_xyy_xzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 170); 

                auto tg_xyy_yyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 171); 

                auto tg_xyy_yyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 172); 

                auto tg_xyy_yyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 173); 

                auto tg_xyy_yyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 174); 

                auto tg_xyy_yyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 175); 

                auto tg_xyy_yyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 176); 

                auto tg_xyy_yyzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 177); 

                auto tg_xyy_yzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 178); 

                auto tg_xyy_zzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 179); 

                auto tg_xyz_xxxxxxxx_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 180); 

                auto tg_xyz_xxxxxxxy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 181); 

                auto tg_xyz_xxxxxxxz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 182); 

                auto tg_xyz_xxxxxxyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 183); 

                auto tg_xyz_xxxxxxyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 184); 

                auto tg_xyz_xxxxxxzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 185); 

                auto tg_xyz_xxxxxyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 186); 

                auto tg_xyz_xxxxxyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 187); 

                auto tg_xyz_xxxxxyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 188); 

                auto tg_xyz_xxxxxzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 189); 

                auto tg_xyz_xxxxyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 190); 

                auto tg_xyz_xxxxyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 191); 

                auto tg_xyz_xxxxyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 192); 

                auto tg_xyz_xxxxyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 193); 

                auto tg_xyz_xxxxzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 194); 

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

                // set up pointers to integrals

                auto tg_xxyy_xxxxyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 147); 

                auto tg_xxyy_xxxxyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 148); 

                auto tg_xxyy_xxxxzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 149); 

                auto tg_xxyy_xxxyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 150); 

                auto tg_xxyy_xxxyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 151); 

                auto tg_xxyy_xxxyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 152); 

                auto tg_xxyy_xxxyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 153); 

                auto tg_xxyy_xxxyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 154); 

                auto tg_xxyy_xxxzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 155); 

                auto tg_xxyy_xxyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 156); 

                auto tg_xxyy_xxyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 157); 

                auto tg_xxyy_xxyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 158); 

                auto tg_xxyy_xxyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 159); 

                auto tg_xxyy_xxyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 160); 

                auto tg_xxyy_xxyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 161); 

                auto tg_xxyy_xxzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 162); 

                auto tg_xxyy_xyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 163); 

                auto tg_xxyy_xyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 164); 

                auto tg_xxyy_xyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 165); 

                auto tg_xxyy_xyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 166); 

                auto tg_xxyy_xyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 167); 

                auto tg_xxyy_xyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 168); 

                auto tg_xxyy_xyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 169); 

                auto tg_xxyy_xzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 170); 

                auto tg_xxyy_yyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 171); 

                auto tg_xxyy_yyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 172); 

                auto tg_xxyy_yyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 173); 

                auto tg_xxyy_yyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 174); 

                auto tg_xxyy_yyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 175); 

                auto tg_xxyy_yyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 176); 

                auto tg_xxyy_yyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 177); 

                auto tg_xxyy_yzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 178); 

                auto tg_xxyy_zzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 179); 

                auto tg_xxyz_xxxxxxxx_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 180); 

                auto tg_xxyz_xxxxxxxy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 181); 

                auto tg_xxyz_xxxxxxxz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 182); 

                auto tg_xxyz_xxxxxxyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 183); 

                auto tg_xxyz_xxxxxxyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 184); 

                auto tg_xxyz_xxxxxxzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 185); 

                auto tg_xxyz_xxxxxyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 186); 

                auto tg_xxyz_xxxxxyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 187); 

                auto tg_xxyz_xxxxxyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 188); 

                auto tg_xxyz_xxxxxzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 189); 

                auto tg_xxyz_xxxxyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 190); 

                auto tg_xxyz_xxxxyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 191); 

                auto tg_xxyz_xxxxyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 192); 

                auto tg_xxyz_xxxxyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 193); 

                auto tg_xxyz_xxxxzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 194); 

                // Batch of Integrals (147,195)

                #pragma omp simd aligned(fxn, fza, tg_xxyy_xxxxyyzz_0, tg_xxyy_xxxxyzzz_0, tg_xxyy_xxxxzzzz_0, \
                                         tg_xxyy_xxxyyyyy_0, tg_xxyy_xxxyyyyz_0, tg_xxyy_xxxyyyzz_0, tg_xxyy_xxxyyzzz_0, \
                                         tg_xxyy_xxxyzzzz_0, tg_xxyy_xxxzzzzz_0, tg_xxyy_xxyyyyyy_0, tg_xxyy_xxyyyyyz_0, \
                                         tg_xxyy_xxyyyyzz_0, tg_xxyy_xxyyyzzz_0, tg_xxyy_xxyyzzzz_0, tg_xxyy_xxyzzzzz_0, \
                                         tg_xxyy_xxzzzzzz_0, tg_xxyy_xyyyyyyy_0, tg_xxyy_xyyyyyyz_0, tg_xxyy_xyyyyyzz_0, \
                                         tg_xxyy_xyyyyzzz_0, tg_xxyy_xyyyzzzz_0, tg_xxyy_xyyzzzzz_0, tg_xxyy_xyzzzzzz_0, \
                                         tg_xxyy_xzzzzzzz_0, tg_xxyy_yyyyyyyy_0, tg_xxyy_yyyyyyyz_0, tg_xxyy_yyyyyyzz_0, \
                                         tg_xxyy_yyyyyzzz_0, tg_xxyy_yyyyzzzz_0, tg_xxyy_yyyzzzzz_0, tg_xxyy_yyzzzzzz_0, \
                                         tg_xxyy_yzzzzzzz_0, tg_xxyy_zzzzzzzz_0, tg_xxyz_xxxxxxxx_0, tg_xxyz_xxxxxxxy_0, \
                                         tg_xxyz_xxxxxxxz_0, tg_xxyz_xxxxxxyy_0, tg_xxyz_xxxxxxyz_0, tg_xxyz_xxxxxxzz_0, \
                                         tg_xxyz_xxxxxyyy_0, tg_xxyz_xxxxxyyz_0, tg_xxyz_xxxxxyzz_0, tg_xxyz_xxxxxzzz_0, \
                                         tg_xxyz_xxxxyyyy_0, tg_xxyz_xxxxyyyz_0, tg_xxyz_xxxxyyzz_0, tg_xxyz_xxxxyzzz_0, \
                                         tg_xxyz_xxxxzzzz_0, tg_xyy_xxxxyyzz_0, tg_xyy_xxxxyyzz_1, tg_xyy_xxxxyzzz_0, \
                                         tg_xyy_xxxxyzzz_1, tg_xyy_xxxxzzzz_0, tg_xyy_xxxxzzzz_1, tg_xyy_xxxyyyyy_0, \
                                         tg_xyy_xxxyyyyy_1, tg_xyy_xxxyyyyz_0, tg_xyy_xxxyyyyz_1, tg_xyy_xxxyyyzz_0, \
                                         tg_xyy_xxxyyyzz_1, tg_xyy_xxxyyzz_1, tg_xyy_xxxyyzzz_0, tg_xyy_xxxyyzzz_1, \
                                         tg_xyy_xxxyzzz_1, tg_xyy_xxxyzzzz_0, tg_xyy_xxxyzzzz_1, tg_xyy_xxxzzzz_1, \
                                         tg_xyy_xxxzzzzz_0, tg_xyy_xxxzzzzz_1, tg_xyy_xxyyyyy_1, tg_xyy_xxyyyyyy_0, \
                                         tg_xyy_xxyyyyyy_1, tg_xyy_xxyyyyyz_0, tg_xyy_xxyyyyyz_1, tg_xyy_xxyyyyz_1, \
                                         tg_xyy_xxyyyyzz_0, tg_xyy_xxyyyyzz_1, tg_xyy_xxyyyzz_1, tg_xyy_xxyyyzzz_0, \
                                         tg_xyy_xxyyyzzz_1, tg_xyy_xxyyzzz_1, tg_xyy_xxyyzzzz_0, tg_xyy_xxyyzzzz_1, \
                                         tg_xyy_xxyzzzz_1, tg_xyy_xxyzzzzz_0, tg_xyy_xxyzzzzz_1, tg_xyy_xxzzzzz_1, \
                                         tg_xyy_xxzzzzzz_0, tg_xyy_xxzzzzzz_1, tg_xyy_xyyyyyy_1, tg_xyy_xyyyyyyy_0, \
                                         tg_xyy_xyyyyyyy_1, tg_xyy_xyyyyyyz_0, tg_xyy_xyyyyyyz_1, tg_xyy_xyyyyyz_1, \
                                         tg_xyy_xyyyyyzz_0, tg_xyy_xyyyyyzz_1, tg_xyy_xyyyyzz_1, tg_xyy_xyyyyzzz_0, \
                                         tg_xyy_xyyyyzzz_1, tg_xyy_xyyyzzz_1, tg_xyy_xyyyzzzz_0, tg_xyy_xyyyzzzz_1, \
                                         tg_xyy_xyyzzzz_1, tg_xyy_xyyzzzzz_0, tg_xyy_xyyzzzzz_1, tg_xyy_xyzzzzz_1, \
                                         tg_xyy_xyzzzzzz_0, tg_xyy_xyzzzzzz_1, tg_xyy_xzzzzzz_1, tg_xyy_xzzzzzzz_0, \
                                         tg_xyy_xzzzzzzz_1, tg_xyy_yyyyyyy_1, tg_xyy_yyyyyyyy_0, tg_xyy_yyyyyyyy_1, \
                                         tg_xyy_yyyyyyyz_0, tg_xyy_yyyyyyyz_1, tg_xyy_yyyyyyz_1, tg_xyy_yyyyyyzz_0, \
                                         tg_xyy_yyyyyyzz_1, tg_xyy_yyyyyzz_1, tg_xyy_yyyyyzzz_0, tg_xyy_yyyyyzzz_1, \
                                         tg_xyy_yyyyzzz_1, tg_xyy_yyyyzzzz_0, tg_xyy_yyyyzzzz_1, tg_xyy_yyyzzzz_1, \
                                         tg_xyy_yyyzzzzz_0, tg_xyy_yyyzzzzz_1, tg_xyy_yyzzzzz_1, tg_xyy_yyzzzzzz_0, \
                                         tg_xyy_yyzzzzzz_1, tg_xyy_yzzzzzz_1, tg_xyy_yzzzzzzz_0, tg_xyy_yzzzzzzz_1, \
                                         tg_xyy_zzzzzzz_1, tg_xyy_zzzzzzzz_0, tg_xyy_zzzzzzzz_1, tg_xyz_xxxxxxx_1, \
                                         tg_xyz_xxxxxxxx_0, tg_xyz_xxxxxxxx_1, tg_xyz_xxxxxxxy_0, tg_xyz_xxxxxxxy_1, \
                                         tg_xyz_xxxxxxxz_0, tg_xyz_xxxxxxxz_1, tg_xyz_xxxxxxy_1, tg_xyz_xxxxxxyy_0, \
                                         tg_xyz_xxxxxxyy_1, tg_xyz_xxxxxxyz_0, tg_xyz_xxxxxxyz_1, tg_xyz_xxxxxxz_1, \
                                         tg_xyz_xxxxxxzz_0, tg_xyz_xxxxxxzz_1, tg_xyz_xxxxxyy_1, tg_xyz_xxxxxyyy_0, \
                                         tg_xyz_xxxxxyyy_1, tg_xyz_xxxxxyyz_0, tg_xyz_xxxxxyyz_1, tg_xyz_xxxxxyz_1, \
                                         tg_xyz_xxxxxyzz_0, tg_xyz_xxxxxyzz_1, tg_xyz_xxxxxzz_1, tg_xyz_xxxxxzzz_0, \
                                         tg_xyz_xxxxxzzz_1, tg_xyz_xxxxyyy_1, tg_xyz_xxxxyyyy_0, tg_xyz_xxxxyyyy_1, \
                                         tg_xyz_xxxxyyyz_0, tg_xyz_xxxxyyyz_1, tg_xyz_xxxxyyz_1, tg_xyz_xxxxyyzz_0, \
                                         tg_xyz_xxxxyyzz_1, tg_xyz_xxxxyzz_1, tg_xyz_xxxxyzzz_0, tg_xyz_xxxxyzzz_1, \
                                         tg_xyz_xxxxzzz_1, tg_xyz_xxxxzzzz_0, tg_xyz_xxxxzzzz_1, tg_xyz_xxxyyyy_1, \
                                         tg_xyz_xxxyyyz_1, tg_xyz_xxxyyzz_1, tg_xyz_xxxyzzz_1, tg_xyz_xxxzzzz_1, \
                                         tg_yy_xxxxyyzz_0, tg_yy_xxxxyyzz_1, tg_yy_xxxxyzzz_0, tg_yy_xxxxyzzz_1, \
                                         tg_yy_xxxxzzzz_0, tg_yy_xxxxzzzz_1, tg_yy_xxxyyyyy_0, tg_yy_xxxyyyyy_1, \
                                         tg_yy_xxxyyyyz_0, tg_yy_xxxyyyyz_1, tg_yy_xxxyyyzz_0, tg_yy_xxxyyyzz_1, \
                                         tg_yy_xxxyyzzz_0, tg_yy_xxxyyzzz_1, tg_yy_xxxyzzzz_0, tg_yy_xxxyzzzz_1, \
                                         tg_yy_xxxzzzzz_0, tg_yy_xxxzzzzz_1, tg_yy_xxyyyyyy_0, tg_yy_xxyyyyyy_1, \
                                         tg_yy_xxyyyyyz_0, tg_yy_xxyyyyyz_1, tg_yy_xxyyyyzz_0, tg_yy_xxyyyyzz_1, \
                                         tg_yy_xxyyyzzz_0, tg_yy_xxyyyzzz_1, tg_yy_xxyyzzzz_0, tg_yy_xxyyzzzz_1, \
                                         tg_yy_xxyzzzzz_0, tg_yy_xxyzzzzz_1, tg_yy_xxzzzzzz_0, tg_yy_xxzzzzzz_1, \
                                         tg_yy_xyyyyyyy_0, tg_yy_xyyyyyyy_1, tg_yy_xyyyyyyz_0, tg_yy_xyyyyyyz_1, \
                                         tg_yy_xyyyyyzz_0, tg_yy_xyyyyyzz_1, tg_yy_xyyyyzzz_0, tg_yy_xyyyyzzz_1, \
                                         tg_yy_xyyyzzzz_0, tg_yy_xyyyzzzz_1, tg_yy_xyyzzzzz_0, tg_yy_xyyzzzzz_1, \
                                         tg_yy_xyzzzzzz_0, tg_yy_xyzzzzzz_1, tg_yy_xzzzzzzz_0, tg_yy_xzzzzzzz_1, \
                                         tg_yy_yyyyyyyy_0, tg_yy_yyyyyyyy_1, tg_yy_yyyyyyyz_0, tg_yy_yyyyyyyz_1, \
                                         tg_yy_yyyyyyzz_0, tg_yy_yyyyyyzz_1, tg_yy_yyyyyzzz_0, tg_yy_yyyyyzzz_1, \
                                         tg_yy_yyyyzzzz_0, tg_yy_yyyyzzzz_1, tg_yy_yyyzzzzz_0, tg_yy_yyyzzzzz_1, \
                                         tg_yy_yyzzzzzz_0, tg_yy_yyzzzzzz_1, tg_yy_yzzzzzzz_0, tg_yy_yzzzzzzz_1, \
                                         tg_yy_zzzzzzzz_0, tg_yy_zzzzzzzz_1, tg_yz_xxxxxxxx_0, tg_yz_xxxxxxxx_1, \
                                         tg_yz_xxxxxxxy_0, tg_yz_xxxxxxxy_1, tg_yz_xxxxxxxz_0, tg_yz_xxxxxxxz_1, \
                                         tg_yz_xxxxxxyy_0, tg_yz_xxxxxxyy_1, tg_yz_xxxxxxyz_0, tg_yz_xxxxxxyz_1, \
                                         tg_yz_xxxxxxzz_0, tg_yz_xxxxxxzz_1, tg_yz_xxxxxyyy_0, tg_yz_xxxxxyyy_1, \
                                         tg_yz_xxxxxyyz_0, tg_yz_xxxxxyyz_1, tg_yz_xxxxxyzz_0, tg_yz_xxxxxyzz_1, \
                                         tg_yz_xxxxxzzz_0, tg_yz_xxxxxzzz_1, tg_yz_xxxxyyyy_0, tg_yz_xxxxyyyy_1, \
                                         tg_yz_xxxxyyyz_0, tg_yz_xxxxyyyz_1, tg_yz_xxxxyyzz_0, tg_yz_xxxxyyzz_1, \
                                         tg_yz_xxxxyzzz_0, tg_yz_xxxxyzzz_1, tg_yz_xxxxzzzz_0, tg_yz_xxxxzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxyy_xxxxyyzz_0[j] = pb_x * tg_xyy_xxxxyyzz_0[j] + wp_x[j] * tg_xyy_xxxxyyzz_1[j] + 0.5 * fl1_fx * tg_yy_xxxxyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxxxyyzz_1[j] + 2.0 * fl1_fxn * tg_xyy_xxxyyzz_1[j];

                    tg_xxyy_xxxxyzzz_0[j] = pb_x * tg_xyy_xxxxyzzz_0[j] + wp_x[j] * tg_xyy_xxxxyzzz_1[j] + 0.5 * fl1_fx * tg_yy_xxxxyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxxxyzzz_1[j] + 2.0 * fl1_fxn * tg_xyy_xxxyzzz_1[j];

                    tg_xxyy_xxxxzzzz_0[j] = pb_x * tg_xyy_xxxxzzzz_0[j] + wp_x[j] * tg_xyy_xxxxzzzz_1[j] + 0.5 * fl1_fx * tg_yy_xxxxzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxxxzzzz_1[j] + 2.0 * fl1_fxn * tg_xyy_xxxzzzz_1[j];

                    tg_xxyy_xxxyyyyy_0[j] = pb_x * tg_xyy_xxxyyyyy_0[j] + wp_x[j] * tg_xyy_xxxyyyyy_1[j] + 0.5 * fl1_fx * tg_yy_xxxyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxxyyyyy_1[j] + 1.5 * fl1_fxn * tg_xyy_xxyyyyy_1[j];

                    tg_xxyy_xxxyyyyz_0[j] = pb_x * tg_xyy_xxxyyyyz_0[j] + wp_x[j] * tg_xyy_xxxyyyyz_1[j] + 0.5 * fl1_fx * tg_yy_xxxyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxxyyyyz_1[j] + 1.5 * fl1_fxn * tg_xyy_xxyyyyz_1[j];

                    tg_xxyy_xxxyyyzz_0[j] = pb_x * tg_xyy_xxxyyyzz_0[j] + wp_x[j] * tg_xyy_xxxyyyzz_1[j] + 0.5 * fl1_fx * tg_yy_xxxyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_xyy_xxyyyzz_1[j];

                    tg_xxyy_xxxyyzzz_0[j] = pb_x * tg_xyy_xxxyyzzz_0[j] + wp_x[j] * tg_xyy_xxxyyzzz_1[j] + 0.5 * fl1_fx * tg_yy_xxxyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxxyyzzz_1[j] + 1.5 * fl1_fxn * tg_xyy_xxyyzzz_1[j];

                    tg_xxyy_xxxyzzzz_0[j] = pb_x * tg_xyy_xxxyzzzz_0[j] + wp_x[j] * tg_xyy_xxxyzzzz_1[j] + 0.5 * fl1_fx * tg_yy_xxxyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxxyzzzz_1[j] + 1.5 * fl1_fxn * tg_xyy_xxyzzzz_1[j];

                    tg_xxyy_xxxzzzzz_0[j] = pb_x * tg_xyy_xxxzzzzz_0[j] + wp_x[j] * tg_xyy_xxxzzzzz_1[j] + 0.5 * fl1_fx * tg_yy_xxxzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxxzzzzz_1[j] + 1.5 * fl1_fxn * tg_xyy_xxzzzzz_1[j];

                    tg_xxyy_xxyyyyyy_0[j] = pb_x * tg_xyy_xxyyyyyy_0[j] + wp_x[j] * tg_xyy_xxyyyyyy_1[j] + 0.5 * fl1_fx * tg_yy_xxyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxyyyyyy_1[j] + fl1_fxn * tg_xyy_xyyyyyy_1[j];

                    tg_xxyy_xxyyyyyz_0[j] = pb_x * tg_xyy_xxyyyyyz_0[j] + wp_x[j] * tg_xyy_xxyyyyyz_1[j] + 0.5 * fl1_fx * tg_yy_xxyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxyyyyyz_1[j] + fl1_fxn * tg_xyy_xyyyyyz_1[j];

                    tg_xxyy_xxyyyyzz_0[j] = pb_x * tg_xyy_xxyyyyzz_0[j] + wp_x[j] * tg_xyy_xxyyyyzz_1[j] + 0.5 * fl1_fx * tg_yy_xxyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxyyyyzz_1[j] + fl1_fxn * tg_xyy_xyyyyzz_1[j];

                    tg_xxyy_xxyyyzzz_0[j] = pb_x * tg_xyy_xxyyyzzz_0[j] + wp_x[j] * tg_xyy_xxyyyzzz_1[j] + 0.5 * fl1_fx * tg_yy_xxyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxyyyzzz_1[j] + fl1_fxn * tg_xyy_xyyyzzz_1[j];

                    tg_xxyy_xxyyzzzz_0[j] = pb_x * tg_xyy_xxyyzzzz_0[j] + wp_x[j] * tg_xyy_xxyyzzzz_1[j] + 0.5 * fl1_fx * tg_yy_xxyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxyyzzzz_1[j] + fl1_fxn * tg_xyy_xyyzzzz_1[j];

                    tg_xxyy_xxyzzzzz_0[j] = pb_x * tg_xyy_xxyzzzzz_0[j] + wp_x[j] * tg_xyy_xxyzzzzz_1[j] + 0.5 * fl1_fx * tg_yy_xxyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxyzzzzz_1[j] + fl1_fxn * tg_xyy_xyzzzzz_1[j];

                    tg_xxyy_xxzzzzzz_0[j] = pb_x * tg_xyy_xxzzzzzz_0[j] + wp_x[j] * tg_xyy_xxzzzzzz_1[j] + 0.5 * fl1_fx * tg_yy_xxzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxzzzzzz_1[j] + fl1_fxn * tg_xyy_xzzzzzz_1[j];

                    tg_xxyy_xyyyyyyy_0[j] = pb_x * tg_xyy_xyyyyyyy_0[j] + wp_x[j] * tg_xyy_xyyyyyyy_1[j] + 0.5 * fl1_fx * tg_yy_xyyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xyyyyyyy_1[j] + 0.5 * fl1_fxn * tg_xyy_yyyyyyy_1[j];

                    tg_xxyy_xyyyyyyz_0[j] = pb_x * tg_xyy_xyyyyyyz_0[j] + wp_x[j] * tg_xyy_xyyyyyyz_1[j] + 0.5 * fl1_fx * tg_yy_xyyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xyyyyyyz_1[j] + 0.5 * fl1_fxn * tg_xyy_yyyyyyz_1[j];

                    tg_xxyy_xyyyyyzz_0[j] = pb_x * tg_xyy_xyyyyyzz_0[j] + wp_x[j] * tg_xyy_xyyyyyzz_1[j] + 0.5 * fl1_fx * tg_yy_xyyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xyyyyyzz_1[j] + 0.5 * fl1_fxn * tg_xyy_yyyyyzz_1[j];

                    tg_xxyy_xyyyyzzz_0[j] = pb_x * tg_xyy_xyyyyzzz_0[j] + wp_x[j] * tg_xyy_xyyyyzzz_1[j] + 0.5 * fl1_fx * tg_yy_xyyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xyyyyzzz_1[j] + 0.5 * fl1_fxn * tg_xyy_yyyyzzz_1[j];

                    tg_xxyy_xyyyzzzz_0[j] = pb_x * tg_xyy_xyyyzzzz_0[j] + wp_x[j] * tg_xyy_xyyyzzzz_1[j] + 0.5 * fl1_fx * tg_yy_xyyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xyyyzzzz_1[j] + 0.5 * fl1_fxn * tg_xyy_yyyzzzz_1[j];

                    tg_xxyy_xyyzzzzz_0[j] = pb_x * tg_xyy_xyyzzzzz_0[j] + wp_x[j] * tg_xyy_xyyzzzzz_1[j] + 0.5 * fl1_fx * tg_yy_xyyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xyyzzzzz_1[j] + 0.5 * fl1_fxn * tg_xyy_yyzzzzz_1[j];

                    tg_xxyy_xyzzzzzz_0[j] = pb_x * tg_xyy_xyzzzzzz_0[j] + wp_x[j] * tg_xyy_xyzzzzzz_1[j] + 0.5 * fl1_fx * tg_yy_xyzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xyy_yzzzzzz_1[j];

                    tg_xxyy_xzzzzzzz_0[j] = pb_x * tg_xyy_xzzzzzzz_0[j] + wp_x[j] * tg_xyy_xzzzzzzz_1[j] + 0.5 * fl1_fx * tg_yy_xzzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xyy_zzzzzzz_1[j];

                    tg_xxyy_yyyyyyyy_0[j] = pb_x * tg_xyy_yyyyyyyy_0[j] + wp_x[j] * tg_xyy_yyyyyyyy_1[j] + 0.5 * fl1_fx * tg_yy_yyyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_yyyyyyyy_1[j];

                    tg_xxyy_yyyyyyyz_0[j] = pb_x * tg_xyy_yyyyyyyz_0[j] + wp_x[j] * tg_xyy_yyyyyyyz_1[j] + 0.5 * fl1_fx * tg_yy_yyyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_yyyyyyyz_1[j];

                    tg_xxyy_yyyyyyzz_0[j] = pb_x * tg_xyy_yyyyyyzz_0[j] + wp_x[j] * tg_xyy_yyyyyyzz_1[j] + 0.5 * fl1_fx * tg_yy_yyyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_yyyyyyzz_1[j];

                    tg_xxyy_yyyyyzzz_0[j] = pb_x * tg_xyy_yyyyyzzz_0[j] + wp_x[j] * tg_xyy_yyyyyzzz_1[j] + 0.5 * fl1_fx * tg_yy_yyyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_yyyyyzzz_1[j];

                    tg_xxyy_yyyyzzzz_0[j] = pb_x * tg_xyy_yyyyzzzz_0[j] + wp_x[j] * tg_xyy_yyyyzzzz_1[j] + 0.5 * fl1_fx * tg_yy_yyyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_yyyyzzzz_1[j];

                    tg_xxyy_yyyzzzzz_0[j] = pb_x * tg_xyy_yyyzzzzz_0[j] + wp_x[j] * tg_xyy_yyyzzzzz_1[j] + 0.5 * fl1_fx * tg_yy_yyyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_yyyzzzzz_1[j];

                    tg_xxyy_yyzzzzzz_0[j] = pb_x * tg_xyy_yyzzzzzz_0[j] + wp_x[j] * tg_xyy_yyzzzzzz_1[j] + 0.5 * fl1_fx * tg_yy_yyzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_yyzzzzzz_1[j];

                    tg_xxyy_yzzzzzzz_0[j] = pb_x * tg_xyy_yzzzzzzz_0[j] + wp_x[j] * tg_xyy_yzzzzzzz_1[j] + 0.5 * fl1_fx * tg_yy_yzzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_yzzzzzzz_1[j];

                    tg_xxyy_zzzzzzzz_0[j] = pb_x * tg_xyy_zzzzzzzz_0[j] + wp_x[j] * tg_xyy_zzzzzzzz_1[j] + 0.5 * fl1_fx * tg_yy_zzzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_zzzzzzzz_1[j];

                    tg_xxyz_xxxxxxxx_0[j] = pb_x * tg_xyz_xxxxxxxx_0[j] + wp_x[j] * tg_xyz_xxxxxxxx_1[j] + 0.5 * fl1_fx * tg_yz_xxxxxxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxxxxxxx_1[j] + 4.0 * fl1_fxn * tg_xyz_xxxxxxx_1[j];

                    tg_xxyz_xxxxxxxy_0[j] = pb_x * tg_xyz_xxxxxxxy_0[j] + wp_x[j] * tg_xyz_xxxxxxxy_1[j] + 0.5 * fl1_fx * tg_yz_xxxxxxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxxxxxxy_1[j] + 3.5 * fl1_fxn * tg_xyz_xxxxxxy_1[j];

                    tg_xxyz_xxxxxxxz_0[j] = pb_x * tg_xyz_xxxxxxxz_0[j] + wp_x[j] * tg_xyz_xxxxxxxz_1[j] + 0.5 * fl1_fx * tg_yz_xxxxxxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxxxxxxz_1[j] + 3.5 * fl1_fxn * tg_xyz_xxxxxxz_1[j];

                    tg_xxyz_xxxxxxyy_0[j] = pb_x * tg_xyz_xxxxxxyy_0[j] + wp_x[j] * tg_xyz_xxxxxxyy_1[j] + 0.5 * fl1_fx * tg_yz_xxxxxxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxxxxxyy_1[j] + 3.0 * fl1_fxn * tg_xyz_xxxxxyy_1[j];

                    tg_xxyz_xxxxxxyz_0[j] = pb_x * tg_xyz_xxxxxxyz_0[j] + wp_x[j] * tg_xyz_xxxxxxyz_1[j] + 0.5 * fl1_fx * tg_yz_xxxxxxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxxxxxyz_1[j] + 3.0 * fl1_fxn * tg_xyz_xxxxxyz_1[j];

                    tg_xxyz_xxxxxxzz_0[j] = pb_x * tg_xyz_xxxxxxzz_0[j] + wp_x[j] * tg_xyz_xxxxxxzz_1[j] + 0.5 * fl1_fx * tg_yz_xxxxxxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxxxxxzz_1[j] + 3.0 * fl1_fxn * tg_xyz_xxxxxzz_1[j];

                    tg_xxyz_xxxxxyyy_0[j] = pb_x * tg_xyz_xxxxxyyy_0[j] + wp_x[j] * tg_xyz_xxxxxyyy_1[j] + 0.5 * fl1_fx * tg_yz_xxxxxyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxxxxyyy_1[j] + 2.5 * fl1_fxn * tg_xyz_xxxxyyy_1[j];

                    tg_xxyz_xxxxxyyz_0[j] = pb_x * tg_xyz_xxxxxyyz_0[j] + wp_x[j] * tg_xyz_xxxxxyyz_1[j] + 0.5 * fl1_fx * tg_yz_xxxxxyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxxxxyyz_1[j] + 2.5 * fl1_fxn * tg_xyz_xxxxyyz_1[j];

                    tg_xxyz_xxxxxyzz_0[j] = pb_x * tg_xyz_xxxxxyzz_0[j] + wp_x[j] * tg_xyz_xxxxxyzz_1[j] + 0.5 * fl1_fx * tg_yz_xxxxxyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxxxxyzz_1[j] + 2.5 * fl1_fxn * tg_xyz_xxxxyzz_1[j];

                    tg_xxyz_xxxxxzzz_0[j] = pb_x * tg_xyz_xxxxxzzz_0[j] + wp_x[j] * tg_xyz_xxxxxzzz_1[j] + 0.5 * fl1_fx * tg_yz_xxxxxzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxxxxzzz_1[j] + 2.5 * fl1_fxn * tg_xyz_xxxxzzz_1[j];

                    tg_xxyz_xxxxyyyy_0[j] = pb_x * tg_xyz_xxxxyyyy_0[j] + wp_x[j] * tg_xyz_xxxxyyyy_1[j] + 0.5 * fl1_fx * tg_yz_xxxxyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_xyz_xxxyyyy_1[j];

                    tg_xxyz_xxxxyyyz_0[j] = pb_x * tg_xyz_xxxxyyyz_0[j] + wp_x[j] * tg_xyz_xxxxyyyz_1[j] + 0.5 * fl1_fx * tg_yz_xxxxyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxxxyyyz_1[j] + 2.0 * fl1_fxn * tg_xyz_xxxyyyz_1[j];

                    tg_xxyz_xxxxyyzz_0[j] = pb_x * tg_xyz_xxxxyyzz_0[j] + wp_x[j] * tg_xyz_xxxxyyzz_1[j] + 0.5 * fl1_fx * tg_yz_xxxxyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxxxyyzz_1[j] + 2.0 * fl1_fxn * tg_xyz_xxxyyzz_1[j];

                    tg_xxyz_xxxxyzzz_0[j] = pb_x * tg_xyz_xxxxyzzz_0[j] + wp_x[j] * tg_xyz_xxxxyzzz_1[j] + 0.5 * fl1_fx * tg_yz_xxxxyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxxxyzzz_1[j] + 2.0 * fl1_fxn * tg_xyz_xxxyzzz_1[j];

                    tg_xxyz_xxxxzzzz_0[j] = pb_x * tg_xyz_xxxxzzzz_0[j] + wp_x[j] * tg_xyz_xxxxzzzz_1[j] + 0.5 * fl1_fx * tg_yz_xxxxzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxxxzzzz_1[j] + 2.0 * fl1_fxn * tg_xyz_xxxzzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSL_195_243(      CMemBlock2D<double>& primBuffer,
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
                                             {4, -1, -1, -1},
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_2_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_2_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

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

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                // set up pointers to auxilary integrals

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

                auto tg_xyz_xxxyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 195); 

                auto tg_xyz_xxxyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 196); 

                auto tg_xyz_xxxyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 197); 

                auto tg_xyz_xxxyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 198); 

                auto tg_xyz_xxxyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 199); 

                auto tg_xyz_xxxzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 200); 

                auto tg_xyz_xxyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 201); 

                auto tg_xyz_xxyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 202); 

                auto tg_xyz_xxyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 203); 

                auto tg_xyz_xxyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 204); 

                auto tg_xyz_xxyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 205); 

                auto tg_xyz_xxyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 206); 

                auto tg_xyz_xxzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 207); 

                auto tg_xyz_xyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 208); 

                auto tg_xyz_xyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 209); 

                auto tg_xyz_xyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 210); 

                auto tg_xyz_xyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 211); 

                auto tg_xyz_xyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 212); 

                auto tg_xyz_xyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 213); 

                auto tg_xyz_xyzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 214); 

                auto tg_xyz_xzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 215); 

                auto tg_xyz_yyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 216); 

                auto tg_xyz_yyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 217); 

                auto tg_xyz_yyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 218); 

                auto tg_xyz_yyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 219); 

                auto tg_xyz_yyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 220); 

                auto tg_xyz_yyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 221); 

                auto tg_xyz_yyzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 222); 

                auto tg_xyz_yzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 223); 

                auto tg_xyz_zzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 224); 

                auto tg_xzz_xxxxxxxx_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 225); 

                auto tg_xzz_xxxxxxxy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 226); 

                auto tg_xzz_xxxxxxxz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 227); 

                auto tg_xzz_xxxxxxyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 228); 

                auto tg_xzz_xxxxxxyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 229); 

                auto tg_xzz_xxxxxxzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 230); 

                auto tg_xzz_xxxxxyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 231); 

                auto tg_xzz_xxxxxyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 232); 

                auto tg_xzz_xxxxxyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 233); 

                auto tg_xzz_xxxxxzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 234); 

                auto tg_xzz_xxxxyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 235); 

                auto tg_xzz_xxxxyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 236); 

                auto tg_xzz_xxxxyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 237); 

                auto tg_xzz_xxxxyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 238); 

                auto tg_xzz_xxxxzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 239); 

                auto tg_xzz_xxxyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 240); 

                auto tg_xzz_xxxyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 241); 

                auto tg_xzz_xxxyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 242); 

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

                auto tg_xzz_xxxyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 190); 

                auto tg_xzz_xxxyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 191); 

                auto tg_xzz_xxxyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 192); 

                auto tg_xzz_xxxyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 193); 

                auto tg_xzz_xxxzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 194); 

                auto tg_xzz_xxyyyyy_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 195); 

                auto tg_xzz_xxyyyyz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 196); 

                auto tg_xzz_xxyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 197); 

                // set up pointers to integrals

                auto tg_xxyz_xxxyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 195); 

                auto tg_xxyz_xxxyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 196); 

                auto tg_xxyz_xxxyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 197); 

                auto tg_xxyz_xxxyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 198); 

                auto tg_xxyz_xxxyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 199); 

                auto tg_xxyz_xxxzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 200); 

                auto tg_xxyz_xxyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 201); 

                auto tg_xxyz_xxyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 202); 

                auto tg_xxyz_xxyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 203); 

                auto tg_xxyz_xxyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 204); 

                auto tg_xxyz_xxyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 205); 

                auto tg_xxyz_xxyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 206); 

                auto tg_xxyz_xxzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 207); 

                auto tg_xxyz_xyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 208); 

                auto tg_xxyz_xyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 209); 

                auto tg_xxyz_xyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 210); 

                auto tg_xxyz_xyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 211); 

                auto tg_xxyz_xyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 212); 

                auto tg_xxyz_xyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 213); 

                auto tg_xxyz_xyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 214); 

                auto tg_xxyz_xzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 215); 

                auto tg_xxyz_yyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 216); 

                auto tg_xxyz_yyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 217); 

                auto tg_xxyz_yyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 218); 

                auto tg_xxyz_yyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 219); 

                auto tg_xxyz_yyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 220); 

                auto tg_xxyz_yyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 221); 

                auto tg_xxyz_yyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 222); 

                auto tg_xxyz_yzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 223); 

                auto tg_xxyz_zzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 224); 

                auto tg_xxzz_xxxxxxxx_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 225); 

                auto tg_xxzz_xxxxxxxy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 226); 

                auto tg_xxzz_xxxxxxxz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 227); 

                auto tg_xxzz_xxxxxxyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 228); 

                auto tg_xxzz_xxxxxxyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 229); 

                auto tg_xxzz_xxxxxxzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 230); 

                auto tg_xxzz_xxxxxyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 231); 

                auto tg_xxzz_xxxxxyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 232); 

                auto tg_xxzz_xxxxxyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 233); 

                auto tg_xxzz_xxxxxzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 234); 

                auto tg_xxzz_xxxxyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 235); 

                auto tg_xxzz_xxxxyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 236); 

                auto tg_xxzz_xxxxyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 237); 

                auto tg_xxzz_xxxxyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 238); 

                auto tg_xxzz_xxxxzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 239); 

                auto tg_xxzz_xxxyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 240); 

                auto tg_xxzz_xxxyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 241); 

                auto tg_xxzz_xxxyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 242); 

                // Batch of Integrals (195,243)

                #pragma omp simd aligned(fxn, fza, tg_xxyz_xxxyyyyy_0, tg_xxyz_xxxyyyyz_0, tg_xxyz_xxxyyyzz_0, \
                                         tg_xxyz_xxxyyzzz_0, tg_xxyz_xxxyzzzz_0, tg_xxyz_xxxzzzzz_0, tg_xxyz_xxyyyyyy_0, \
                                         tg_xxyz_xxyyyyyz_0, tg_xxyz_xxyyyyzz_0, tg_xxyz_xxyyyzzz_0, tg_xxyz_xxyyzzzz_0, \
                                         tg_xxyz_xxyzzzzz_0, tg_xxyz_xxzzzzzz_0, tg_xxyz_xyyyyyyy_0, tg_xxyz_xyyyyyyz_0, \
                                         tg_xxyz_xyyyyyzz_0, tg_xxyz_xyyyyzzz_0, tg_xxyz_xyyyzzzz_0, tg_xxyz_xyyzzzzz_0, \
                                         tg_xxyz_xyzzzzzz_0, tg_xxyz_xzzzzzzz_0, tg_xxyz_yyyyyyyy_0, tg_xxyz_yyyyyyyz_0, \
                                         tg_xxyz_yyyyyyzz_0, tg_xxyz_yyyyyzzz_0, tg_xxyz_yyyyzzzz_0, tg_xxyz_yyyzzzzz_0, \
                                         tg_xxyz_yyzzzzzz_0, tg_xxyz_yzzzzzzz_0, tg_xxyz_zzzzzzzz_0, tg_xxzz_xxxxxxxx_0, \
                                         tg_xxzz_xxxxxxxy_0, tg_xxzz_xxxxxxxz_0, tg_xxzz_xxxxxxyy_0, tg_xxzz_xxxxxxyz_0, \
                                         tg_xxzz_xxxxxxzz_0, tg_xxzz_xxxxxyyy_0, tg_xxzz_xxxxxyyz_0, tg_xxzz_xxxxxyzz_0, \
                                         tg_xxzz_xxxxxzzz_0, tg_xxzz_xxxxyyyy_0, tg_xxzz_xxxxyyyz_0, tg_xxzz_xxxxyyzz_0, \
                                         tg_xxzz_xxxxyzzz_0, tg_xxzz_xxxxzzzz_0, tg_xxzz_xxxyyyyy_0, tg_xxzz_xxxyyyyz_0, \
                                         tg_xxzz_xxxyyyzz_0, tg_xyz_xxxyyyyy_0, tg_xyz_xxxyyyyy_1, tg_xyz_xxxyyyyz_0, \
                                         tg_xyz_xxxyyyyz_1, tg_xyz_xxxyyyzz_0, tg_xyz_xxxyyyzz_1, tg_xyz_xxxyyzzz_0, \
                                         tg_xyz_xxxyyzzz_1, tg_xyz_xxxyzzzz_0, tg_xyz_xxxyzzzz_1, tg_xyz_xxxzzzzz_0, \
                                         tg_xyz_xxxzzzzz_1, tg_xyz_xxyyyyy_1, tg_xyz_xxyyyyyy_0, tg_xyz_xxyyyyyy_1, \
                                         tg_xyz_xxyyyyyz_0, tg_xyz_xxyyyyyz_1, tg_xyz_xxyyyyz_1, tg_xyz_xxyyyyzz_0, \
                                         tg_xyz_xxyyyyzz_1, tg_xyz_xxyyyzz_1, tg_xyz_xxyyyzzz_0, tg_xyz_xxyyyzzz_1, \
                                         tg_xyz_xxyyzzz_1, tg_xyz_xxyyzzzz_0, tg_xyz_xxyyzzzz_1, tg_xyz_xxyzzzz_1, \
                                         tg_xyz_xxyzzzzz_0, tg_xyz_xxyzzzzz_1, tg_xyz_xxzzzzz_1, tg_xyz_xxzzzzzz_0, \
                                         tg_xyz_xxzzzzzz_1, tg_xyz_xyyyyyy_1, tg_xyz_xyyyyyyy_0, tg_xyz_xyyyyyyy_1, \
                                         tg_xyz_xyyyyyyz_0, tg_xyz_xyyyyyyz_1, tg_xyz_xyyyyyz_1, tg_xyz_xyyyyyzz_0, \
                                         tg_xyz_xyyyyyzz_1, tg_xyz_xyyyyzz_1, tg_xyz_xyyyyzzz_0, tg_xyz_xyyyyzzz_1, \
                                         tg_xyz_xyyyzzz_1, tg_xyz_xyyyzzzz_0, tg_xyz_xyyyzzzz_1, tg_xyz_xyyzzzz_1, \
                                         tg_xyz_xyyzzzzz_0, tg_xyz_xyyzzzzz_1, tg_xyz_xyzzzzz_1, tg_xyz_xyzzzzzz_0, \
                                         tg_xyz_xyzzzzzz_1, tg_xyz_xzzzzzz_1, tg_xyz_xzzzzzzz_0, tg_xyz_xzzzzzzz_1, \
                                         tg_xyz_yyyyyyy_1, tg_xyz_yyyyyyyy_0, tg_xyz_yyyyyyyy_1, tg_xyz_yyyyyyyz_0, \
                                         tg_xyz_yyyyyyyz_1, tg_xyz_yyyyyyz_1, tg_xyz_yyyyyyzz_0, tg_xyz_yyyyyyzz_1, \
                                         tg_xyz_yyyyyzz_1, tg_xyz_yyyyyzzz_0, tg_xyz_yyyyyzzz_1, tg_xyz_yyyyzzz_1, \
                                         tg_xyz_yyyyzzzz_0, tg_xyz_yyyyzzzz_1, tg_xyz_yyyzzzz_1, tg_xyz_yyyzzzzz_0, \
                                         tg_xyz_yyyzzzzz_1, tg_xyz_yyzzzzz_1, tg_xyz_yyzzzzzz_0, tg_xyz_yyzzzzzz_1, \
                                         tg_xyz_yzzzzzz_1, tg_xyz_yzzzzzzz_0, tg_xyz_yzzzzzzz_1, tg_xyz_zzzzzzz_1, \
                                         tg_xyz_zzzzzzzz_0, tg_xyz_zzzzzzzz_1, tg_xzz_xxxxxxx_1, tg_xzz_xxxxxxxx_0, \
                                         tg_xzz_xxxxxxxx_1, tg_xzz_xxxxxxxy_0, tg_xzz_xxxxxxxy_1, tg_xzz_xxxxxxxz_0, \
                                         tg_xzz_xxxxxxxz_1, tg_xzz_xxxxxxy_1, tg_xzz_xxxxxxyy_0, tg_xzz_xxxxxxyy_1, \
                                         tg_xzz_xxxxxxyz_0, tg_xzz_xxxxxxyz_1, tg_xzz_xxxxxxz_1, tg_xzz_xxxxxxzz_0, \
                                         tg_xzz_xxxxxxzz_1, tg_xzz_xxxxxyy_1, tg_xzz_xxxxxyyy_0, tg_xzz_xxxxxyyy_1, \
                                         tg_xzz_xxxxxyyz_0, tg_xzz_xxxxxyyz_1, tg_xzz_xxxxxyz_1, tg_xzz_xxxxxyzz_0, \
                                         tg_xzz_xxxxxyzz_1, tg_xzz_xxxxxzz_1, tg_xzz_xxxxxzzz_0, tg_xzz_xxxxxzzz_1, \
                                         tg_xzz_xxxxyyy_1, tg_xzz_xxxxyyyy_0, tg_xzz_xxxxyyyy_1, tg_xzz_xxxxyyyz_0, \
                                         tg_xzz_xxxxyyyz_1, tg_xzz_xxxxyyz_1, tg_xzz_xxxxyyzz_0, tg_xzz_xxxxyyzz_1, \
                                         tg_xzz_xxxxyzz_1, tg_xzz_xxxxyzzz_0, tg_xzz_xxxxyzzz_1, tg_xzz_xxxxzzz_1, \
                                         tg_xzz_xxxxzzzz_0, tg_xzz_xxxxzzzz_1, tg_xzz_xxxyyyy_1, tg_xzz_xxxyyyyy_0, \
                                         tg_xzz_xxxyyyyy_1, tg_xzz_xxxyyyyz_0, tg_xzz_xxxyyyyz_1, tg_xzz_xxxyyyz_1, \
                                         tg_xzz_xxxyyyzz_0, tg_xzz_xxxyyyzz_1, tg_xzz_xxxyyzz_1, tg_xzz_xxxyzzz_1, \
                                         tg_xzz_xxxzzzz_1, tg_xzz_xxyyyyy_1, tg_xzz_xxyyyyz_1, tg_xzz_xxyyyzz_1, \
                                         tg_yz_xxxyyyyy_0, tg_yz_xxxyyyyy_1, tg_yz_xxxyyyyz_0, tg_yz_xxxyyyyz_1, \
                                         tg_yz_xxxyyyzz_0, tg_yz_xxxyyyzz_1, tg_yz_xxxyyzzz_0, tg_yz_xxxyyzzz_1, \
                                         tg_yz_xxxyzzzz_0, tg_yz_xxxyzzzz_1, tg_yz_xxxzzzzz_0, tg_yz_xxxzzzzz_1, \
                                         tg_yz_xxyyyyyy_0, tg_yz_xxyyyyyy_1, tg_yz_xxyyyyyz_0, tg_yz_xxyyyyyz_1, \
                                         tg_yz_xxyyyyzz_0, tg_yz_xxyyyyzz_1, tg_yz_xxyyyzzz_0, tg_yz_xxyyyzzz_1, \
                                         tg_yz_xxyyzzzz_0, tg_yz_xxyyzzzz_1, tg_yz_xxyzzzzz_0, tg_yz_xxyzzzzz_1, \
                                         tg_yz_xxzzzzzz_0, tg_yz_xxzzzzzz_1, tg_yz_xyyyyyyy_0, tg_yz_xyyyyyyy_1, \
                                         tg_yz_xyyyyyyz_0, tg_yz_xyyyyyyz_1, tg_yz_xyyyyyzz_0, tg_yz_xyyyyyzz_1, \
                                         tg_yz_xyyyyzzz_0, tg_yz_xyyyyzzz_1, tg_yz_xyyyzzzz_0, tg_yz_xyyyzzzz_1, \
                                         tg_yz_xyyzzzzz_0, tg_yz_xyyzzzzz_1, tg_yz_xyzzzzzz_0, tg_yz_xyzzzzzz_1, \
                                         tg_yz_xzzzzzzz_0, tg_yz_xzzzzzzz_1, tg_yz_yyyyyyyy_0, tg_yz_yyyyyyyy_1, \
                                         tg_yz_yyyyyyyz_0, tg_yz_yyyyyyyz_1, tg_yz_yyyyyyzz_0, tg_yz_yyyyyyzz_1, \
                                         tg_yz_yyyyyzzz_0, tg_yz_yyyyyzzz_1, tg_yz_yyyyzzzz_0, tg_yz_yyyyzzzz_1, \
                                         tg_yz_yyyzzzzz_0, tg_yz_yyyzzzzz_1, tg_yz_yyzzzzzz_0, tg_yz_yyzzzzzz_1, \
                                         tg_yz_yzzzzzzz_0, tg_yz_yzzzzzzz_1, tg_yz_zzzzzzzz_0, tg_yz_zzzzzzzz_1, \
                                         tg_zz_xxxxxxxx_0, tg_zz_xxxxxxxx_1, tg_zz_xxxxxxxy_0, tg_zz_xxxxxxxy_1, \
                                         tg_zz_xxxxxxxz_0, tg_zz_xxxxxxxz_1, tg_zz_xxxxxxyy_0, tg_zz_xxxxxxyy_1, \
                                         tg_zz_xxxxxxyz_0, tg_zz_xxxxxxyz_1, tg_zz_xxxxxxzz_0, tg_zz_xxxxxxzz_1, \
                                         tg_zz_xxxxxyyy_0, tg_zz_xxxxxyyy_1, tg_zz_xxxxxyyz_0, tg_zz_xxxxxyyz_1, \
                                         tg_zz_xxxxxyzz_0, tg_zz_xxxxxyzz_1, tg_zz_xxxxxzzz_0, tg_zz_xxxxxzzz_1, \
                                         tg_zz_xxxxyyyy_0, tg_zz_xxxxyyyy_1, tg_zz_xxxxyyyz_0, tg_zz_xxxxyyyz_1, \
                                         tg_zz_xxxxyyzz_0, tg_zz_xxxxyyzz_1, tg_zz_xxxxyzzz_0, tg_zz_xxxxyzzz_1, \
                                         tg_zz_xxxxzzzz_0, tg_zz_xxxxzzzz_1, tg_zz_xxxyyyyy_0, tg_zz_xxxyyyyy_1, \
                                         tg_zz_xxxyyyyz_0, tg_zz_xxxyyyyz_1, tg_zz_xxxyyyzz_0, tg_zz_xxxyyyzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxyz_xxxyyyyy_0[j] = pb_x * tg_xyz_xxxyyyyy_0[j] + wp_x[j] * tg_xyz_xxxyyyyy_1[j] + 0.5 * fl1_fx * tg_yz_xxxyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxxyyyyy_1[j] + 1.5 * fl1_fxn * tg_xyz_xxyyyyy_1[j];

                    tg_xxyz_xxxyyyyz_0[j] = pb_x * tg_xyz_xxxyyyyz_0[j] + wp_x[j] * tg_xyz_xxxyyyyz_1[j] + 0.5 * fl1_fx * tg_yz_xxxyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxxyyyyz_1[j] + 1.5 * fl1_fxn * tg_xyz_xxyyyyz_1[j];

                    tg_xxyz_xxxyyyzz_0[j] = pb_x * tg_xyz_xxxyyyzz_0[j] + wp_x[j] * tg_xyz_xxxyyyzz_1[j] + 0.5 * fl1_fx * tg_yz_xxxyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_xyz_xxyyyzz_1[j];

                    tg_xxyz_xxxyyzzz_0[j] = pb_x * tg_xyz_xxxyyzzz_0[j] + wp_x[j] * tg_xyz_xxxyyzzz_1[j] + 0.5 * fl1_fx * tg_yz_xxxyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxxyyzzz_1[j] + 1.5 * fl1_fxn * tg_xyz_xxyyzzz_1[j];

                    tg_xxyz_xxxyzzzz_0[j] = pb_x * tg_xyz_xxxyzzzz_0[j] + wp_x[j] * tg_xyz_xxxyzzzz_1[j] + 0.5 * fl1_fx * tg_yz_xxxyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxxyzzzz_1[j] + 1.5 * fl1_fxn * tg_xyz_xxyzzzz_1[j];

                    tg_xxyz_xxxzzzzz_0[j] = pb_x * tg_xyz_xxxzzzzz_0[j] + wp_x[j] * tg_xyz_xxxzzzzz_1[j] + 0.5 * fl1_fx * tg_yz_xxxzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxxzzzzz_1[j] + 1.5 * fl1_fxn * tg_xyz_xxzzzzz_1[j];

                    tg_xxyz_xxyyyyyy_0[j] = pb_x * tg_xyz_xxyyyyyy_0[j] + wp_x[j] * tg_xyz_xxyyyyyy_1[j] + 0.5 * fl1_fx * tg_yz_xxyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxyyyyyy_1[j] + fl1_fxn * tg_xyz_xyyyyyy_1[j];

                    tg_xxyz_xxyyyyyz_0[j] = pb_x * tg_xyz_xxyyyyyz_0[j] + wp_x[j] * tg_xyz_xxyyyyyz_1[j] + 0.5 * fl1_fx * tg_yz_xxyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxyyyyyz_1[j] + fl1_fxn * tg_xyz_xyyyyyz_1[j];

                    tg_xxyz_xxyyyyzz_0[j] = pb_x * tg_xyz_xxyyyyzz_0[j] + wp_x[j] * tg_xyz_xxyyyyzz_1[j] + 0.5 * fl1_fx * tg_yz_xxyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxyyyyzz_1[j] + fl1_fxn * tg_xyz_xyyyyzz_1[j];

                    tg_xxyz_xxyyyzzz_0[j] = pb_x * tg_xyz_xxyyyzzz_0[j] + wp_x[j] * tg_xyz_xxyyyzzz_1[j] + 0.5 * fl1_fx * tg_yz_xxyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxyyyzzz_1[j] + fl1_fxn * tg_xyz_xyyyzzz_1[j];

                    tg_xxyz_xxyyzzzz_0[j] = pb_x * tg_xyz_xxyyzzzz_0[j] + wp_x[j] * tg_xyz_xxyyzzzz_1[j] + 0.5 * fl1_fx * tg_yz_xxyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxyyzzzz_1[j] + fl1_fxn * tg_xyz_xyyzzzz_1[j];

                    tg_xxyz_xxyzzzzz_0[j] = pb_x * tg_xyz_xxyzzzzz_0[j] + wp_x[j] * tg_xyz_xxyzzzzz_1[j] + 0.5 * fl1_fx * tg_yz_xxyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxyzzzzz_1[j] + fl1_fxn * tg_xyz_xyzzzzz_1[j];

                    tg_xxyz_xxzzzzzz_0[j] = pb_x * tg_xyz_xxzzzzzz_0[j] + wp_x[j] * tg_xyz_xxzzzzzz_1[j] + 0.5 * fl1_fx * tg_yz_xxzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxzzzzzz_1[j] + fl1_fxn * tg_xyz_xzzzzzz_1[j];

                    tg_xxyz_xyyyyyyy_0[j] = pb_x * tg_xyz_xyyyyyyy_0[j] + wp_x[j] * tg_xyz_xyyyyyyy_1[j] + 0.5 * fl1_fx * tg_yz_xyyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xyyyyyyy_1[j] + 0.5 * fl1_fxn * tg_xyz_yyyyyyy_1[j];

                    tg_xxyz_xyyyyyyz_0[j] = pb_x * tg_xyz_xyyyyyyz_0[j] + wp_x[j] * tg_xyz_xyyyyyyz_1[j] + 0.5 * fl1_fx * tg_yz_xyyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xyyyyyyz_1[j] + 0.5 * fl1_fxn * tg_xyz_yyyyyyz_1[j];

                    tg_xxyz_xyyyyyzz_0[j] = pb_x * tg_xyz_xyyyyyzz_0[j] + wp_x[j] * tg_xyz_xyyyyyzz_1[j] + 0.5 * fl1_fx * tg_yz_xyyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xyyyyyzz_1[j] + 0.5 * fl1_fxn * tg_xyz_yyyyyzz_1[j];

                    tg_xxyz_xyyyyzzz_0[j] = pb_x * tg_xyz_xyyyyzzz_0[j] + wp_x[j] * tg_xyz_xyyyyzzz_1[j] + 0.5 * fl1_fx * tg_yz_xyyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xyyyyzzz_1[j] + 0.5 * fl1_fxn * tg_xyz_yyyyzzz_1[j];

                    tg_xxyz_xyyyzzzz_0[j] = pb_x * tg_xyz_xyyyzzzz_0[j] + wp_x[j] * tg_xyz_xyyyzzzz_1[j] + 0.5 * fl1_fx * tg_yz_xyyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xyyyzzzz_1[j] + 0.5 * fl1_fxn * tg_xyz_yyyzzzz_1[j];

                    tg_xxyz_xyyzzzzz_0[j] = pb_x * tg_xyz_xyyzzzzz_0[j] + wp_x[j] * tg_xyz_xyyzzzzz_1[j] + 0.5 * fl1_fx * tg_yz_xyyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xyyzzzzz_1[j] + 0.5 * fl1_fxn * tg_xyz_yyzzzzz_1[j];

                    tg_xxyz_xyzzzzzz_0[j] = pb_x * tg_xyz_xyzzzzzz_0[j] + wp_x[j] * tg_xyz_xyzzzzzz_1[j] + 0.5 * fl1_fx * tg_yz_xyzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xyz_yzzzzzz_1[j];

                    tg_xxyz_xzzzzzzz_0[j] = pb_x * tg_xyz_xzzzzzzz_0[j] + wp_x[j] * tg_xyz_xzzzzzzz_1[j] + 0.5 * fl1_fx * tg_yz_xzzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xyz_zzzzzzz_1[j];

                    tg_xxyz_yyyyyyyy_0[j] = pb_x * tg_xyz_yyyyyyyy_0[j] + wp_x[j] * tg_xyz_yyyyyyyy_1[j] + 0.5 * fl1_fx * tg_yz_yyyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_yyyyyyyy_1[j];

                    tg_xxyz_yyyyyyyz_0[j] = pb_x * tg_xyz_yyyyyyyz_0[j] + wp_x[j] * tg_xyz_yyyyyyyz_1[j] + 0.5 * fl1_fx * tg_yz_yyyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_yyyyyyyz_1[j];

                    tg_xxyz_yyyyyyzz_0[j] = pb_x * tg_xyz_yyyyyyzz_0[j] + wp_x[j] * tg_xyz_yyyyyyzz_1[j] + 0.5 * fl1_fx * tg_yz_yyyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_yyyyyyzz_1[j];

                    tg_xxyz_yyyyyzzz_0[j] = pb_x * tg_xyz_yyyyyzzz_0[j] + wp_x[j] * tg_xyz_yyyyyzzz_1[j] + 0.5 * fl1_fx * tg_yz_yyyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_yyyyyzzz_1[j];

                    tg_xxyz_yyyyzzzz_0[j] = pb_x * tg_xyz_yyyyzzzz_0[j] + wp_x[j] * tg_xyz_yyyyzzzz_1[j] + 0.5 * fl1_fx * tg_yz_yyyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_yyyyzzzz_1[j];

                    tg_xxyz_yyyzzzzz_0[j] = pb_x * tg_xyz_yyyzzzzz_0[j] + wp_x[j] * tg_xyz_yyyzzzzz_1[j] + 0.5 * fl1_fx * tg_yz_yyyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_yyyzzzzz_1[j];

                    tg_xxyz_yyzzzzzz_0[j] = pb_x * tg_xyz_yyzzzzzz_0[j] + wp_x[j] * tg_xyz_yyzzzzzz_1[j] + 0.5 * fl1_fx * tg_yz_yyzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_yyzzzzzz_1[j];

                    tg_xxyz_yzzzzzzz_0[j] = pb_x * tg_xyz_yzzzzzzz_0[j] + wp_x[j] * tg_xyz_yzzzzzzz_1[j] + 0.5 * fl1_fx * tg_yz_yzzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_yzzzzzzz_1[j];

                    tg_xxyz_zzzzzzzz_0[j] = pb_x * tg_xyz_zzzzzzzz_0[j] + wp_x[j] * tg_xyz_zzzzzzzz_1[j] + 0.5 * fl1_fx * tg_yz_zzzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_zzzzzzzz_1[j];

                    tg_xxzz_xxxxxxxx_0[j] = pb_x * tg_xzz_xxxxxxxx_0[j] + wp_x[j] * tg_xzz_xxxxxxxx_1[j] + 0.5 * fl1_fx * tg_zz_xxxxxxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxxxxx_1[j] + 4.0 * fl1_fxn * tg_xzz_xxxxxxx_1[j];

                    tg_xxzz_xxxxxxxy_0[j] = pb_x * tg_xzz_xxxxxxxy_0[j] + wp_x[j] * tg_xzz_xxxxxxxy_1[j] + 0.5 * fl1_fx * tg_zz_xxxxxxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxxxxy_1[j] + 3.5 * fl1_fxn * tg_xzz_xxxxxxy_1[j];

                    tg_xxzz_xxxxxxxz_0[j] = pb_x * tg_xzz_xxxxxxxz_0[j] + wp_x[j] * tg_xzz_xxxxxxxz_1[j] + 0.5 * fl1_fx * tg_zz_xxxxxxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxxxxz_1[j] + 3.5 * fl1_fxn * tg_xzz_xxxxxxz_1[j];

                    tg_xxzz_xxxxxxyy_0[j] = pb_x * tg_xzz_xxxxxxyy_0[j] + wp_x[j] * tg_xzz_xxxxxxyy_1[j] + 0.5 * fl1_fx * tg_zz_xxxxxxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxxxyy_1[j] + 3.0 * fl1_fxn * tg_xzz_xxxxxyy_1[j];

                    tg_xxzz_xxxxxxyz_0[j] = pb_x * tg_xzz_xxxxxxyz_0[j] + wp_x[j] * tg_xzz_xxxxxxyz_1[j] + 0.5 * fl1_fx * tg_zz_xxxxxxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxxxyz_1[j] + 3.0 * fl1_fxn * tg_xzz_xxxxxyz_1[j];

                    tg_xxzz_xxxxxxzz_0[j] = pb_x * tg_xzz_xxxxxxzz_0[j] + wp_x[j] * tg_xzz_xxxxxxzz_1[j] + 0.5 * fl1_fx * tg_zz_xxxxxxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxxxzz_1[j] + 3.0 * fl1_fxn * tg_xzz_xxxxxzz_1[j];

                    tg_xxzz_xxxxxyyy_0[j] = pb_x * tg_xzz_xxxxxyyy_0[j] + wp_x[j] * tg_xzz_xxxxxyyy_1[j] + 0.5 * fl1_fx * tg_zz_xxxxxyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxxyyy_1[j] + 2.5 * fl1_fxn * tg_xzz_xxxxyyy_1[j];

                    tg_xxzz_xxxxxyyz_0[j] = pb_x * tg_xzz_xxxxxyyz_0[j] + wp_x[j] * tg_xzz_xxxxxyyz_1[j] + 0.5 * fl1_fx * tg_zz_xxxxxyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxxyyz_1[j] + 2.5 * fl1_fxn * tg_xzz_xxxxyyz_1[j];

                    tg_xxzz_xxxxxyzz_0[j] = pb_x * tg_xzz_xxxxxyzz_0[j] + wp_x[j] * tg_xzz_xxxxxyzz_1[j] + 0.5 * fl1_fx * tg_zz_xxxxxyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxxyzz_1[j] + 2.5 * fl1_fxn * tg_xzz_xxxxyzz_1[j];

                    tg_xxzz_xxxxxzzz_0[j] = pb_x * tg_xzz_xxxxxzzz_0[j] + wp_x[j] * tg_xzz_xxxxxzzz_1[j] + 0.5 * fl1_fx * tg_zz_xxxxxzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxxzzz_1[j] + 2.5 * fl1_fxn * tg_xzz_xxxxzzz_1[j];

                    tg_xxzz_xxxxyyyy_0[j] = pb_x * tg_xzz_xxxxyyyy_0[j] + wp_x[j] * tg_xzz_xxxxyyyy_1[j] + 0.5 * fl1_fx * tg_zz_xxxxyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_xzz_xxxyyyy_1[j];

                    tg_xxzz_xxxxyyyz_0[j] = pb_x * tg_xzz_xxxxyyyz_0[j] + wp_x[j] * tg_xzz_xxxxyyyz_1[j] + 0.5 * fl1_fx * tg_zz_xxxxyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxyyyz_1[j] + 2.0 * fl1_fxn * tg_xzz_xxxyyyz_1[j];

                    tg_xxzz_xxxxyyzz_0[j] = pb_x * tg_xzz_xxxxyyzz_0[j] + wp_x[j] * tg_xzz_xxxxyyzz_1[j] + 0.5 * fl1_fx * tg_zz_xxxxyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxyyzz_1[j] + 2.0 * fl1_fxn * tg_xzz_xxxyyzz_1[j];

                    tg_xxzz_xxxxyzzz_0[j] = pb_x * tg_xzz_xxxxyzzz_0[j] + wp_x[j] * tg_xzz_xxxxyzzz_1[j] + 0.5 * fl1_fx * tg_zz_xxxxyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxyzzz_1[j] + 2.0 * fl1_fxn * tg_xzz_xxxyzzz_1[j];

                    tg_xxzz_xxxxzzzz_0[j] = pb_x * tg_xzz_xxxxzzzz_0[j] + wp_x[j] * tg_xzz_xxxxzzzz_1[j] + 0.5 * fl1_fx * tg_zz_xxxxzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxzzzz_1[j] + 2.0 * fl1_fxn * tg_xzz_xxxzzzz_1[j];

                    tg_xxzz_xxxyyyyy_0[j] = pb_x * tg_xzz_xxxyyyyy_0[j] + wp_x[j] * tg_xzz_xxxyyyyy_1[j] + 0.5 * fl1_fx * tg_zz_xxxyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxyyyyy_1[j] + 1.5 * fl1_fxn * tg_xzz_xxyyyyy_1[j];

                    tg_xxzz_xxxyyyyz_0[j] = pb_x * tg_xzz_xxxyyyyz_0[j] + wp_x[j] * tg_xzz_xxxyyyyz_1[j] + 0.5 * fl1_fx * tg_zz_xxxyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxyyyyz_1[j] + 1.5 * fl1_fxn * tg_xzz_xxyyyyz_1[j];

                    tg_xxzz_xxxyyyzz_0[j] = pb_x * tg_xzz_xxxyyyzz_0[j] + wp_x[j] * tg_xzz_xxxyyyzz_1[j] + 0.5 * fl1_fx * tg_zz_xxxyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_xzz_xxyyyzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSL_243_291(      CMemBlock2D<double>& primBuffer,
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
                                             {4, -1, -1, -1},
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_2_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_2_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

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

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                // set up pointers to auxilary integrals

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

                auto tg_xzz_xxxyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 243); 

                auto tg_xzz_xxxyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 244); 

                auto tg_xzz_xxxzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 245); 

                auto tg_xzz_xxyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 246); 

                auto tg_xzz_xxyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 247); 

                auto tg_xzz_xxyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 248); 

                auto tg_xzz_xxyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 249); 

                auto tg_xzz_xxyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 250); 

                auto tg_xzz_xxyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 251); 

                auto tg_xzz_xxzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 252); 

                auto tg_xzz_xyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 253); 

                auto tg_xzz_xyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 254); 

                auto tg_xzz_xyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 255); 

                auto tg_xzz_xyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 256); 

                auto tg_xzz_xyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 257); 

                auto tg_xzz_xyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 258); 

                auto tg_xzz_xyzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 259); 

                auto tg_xzz_xzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 260); 

                auto tg_xzz_yyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 261); 

                auto tg_xzz_yyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 262); 

                auto tg_xzz_yyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 263); 

                auto tg_xzz_yyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 264); 

                auto tg_xzz_yyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 265); 

                auto tg_xzz_yyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 266); 

                auto tg_xzz_yyzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 267); 

                auto tg_xzz_yzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 268); 

                auto tg_xzz_zzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 269); 

                auto tg_yyy_xxxxxxxx_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 270); 

                auto tg_yyy_xxxxxxxy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 271); 

                auto tg_yyy_xxxxxxxz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 272); 

                auto tg_yyy_xxxxxxyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 273); 

                auto tg_yyy_xxxxxxyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 274); 

                auto tg_yyy_xxxxxxzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 275); 

                auto tg_yyy_xxxxxyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 276); 

                auto tg_yyy_xxxxxyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 277); 

                auto tg_yyy_xxxxxyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 278); 

                auto tg_yyy_xxxxxzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 279); 

                auto tg_yyy_xxxxyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 280); 

                auto tg_yyy_xxxxyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 281); 

                auto tg_yyy_xxxxyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 282); 

                auto tg_yyy_xxxxyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 283); 

                auto tg_yyy_xxxxzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 284); 

                auto tg_yyy_xxxyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 285); 

                auto tg_yyy_xxxyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 286); 

                auto tg_yyy_xxxyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 287); 

                auto tg_yyy_xxxyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 288); 

                auto tg_yyy_xxxyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 289); 

                auto tg_yyy_xxxzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 290); 

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

                // set up pointers to integrals

                auto tg_xxzz_xxxyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 243); 

                auto tg_xxzz_xxxyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 244); 

                auto tg_xxzz_xxxzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 245); 

                auto tg_xxzz_xxyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 246); 

                auto tg_xxzz_xxyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 247); 

                auto tg_xxzz_xxyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 248); 

                auto tg_xxzz_xxyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 249); 

                auto tg_xxzz_xxyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 250); 

                auto tg_xxzz_xxyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 251); 

                auto tg_xxzz_xxzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 252); 

                auto tg_xxzz_xyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 253); 

                auto tg_xxzz_xyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 254); 

                auto tg_xxzz_xyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 255); 

                auto tg_xxzz_xyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 256); 

                auto tg_xxzz_xyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 257); 

                auto tg_xxzz_xyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 258); 

                auto tg_xxzz_xyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 259); 

                auto tg_xxzz_xzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 260); 

                auto tg_xxzz_yyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 261); 

                auto tg_xxzz_yyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 262); 

                auto tg_xxzz_yyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 263); 

                auto tg_xxzz_yyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 264); 

                auto tg_xxzz_yyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 265); 

                auto tg_xxzz_yyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 266); 

                auto tg_xxzz_yyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 267); 

                auto tg_xxzz_yzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 268); 

                auto tg_xxzz_zzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 269); 

                auto tg_xyyy_xxxxxxxx_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 270); 

                auto tg_xyyy_xxxxxxxy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 271); 

                auto tg_xyyy_xxxxxxxz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 272); 

                auto tg_xyyy_xxxxxxyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 273); 

                auto tg_xyyy_xxxxxxyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 274); 

                auto tg_xyyy_xxxxxxzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 275); 

                auto tg_xyyy_xxxxxyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 276); 

                auto tg_xyyy_xxxxxyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 277); 

                auto tg_xyyy_xxxxxyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 278); 

                auto tg_xyyy_xxxxxzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 279); 

                auto tg_xyyy_xxxxyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 280); 

                auto tg_xyyy_xxxxyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 281); 

                auto tg_xyyy_xxxxyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 282); 

                auto tg_xyyy_xxxxyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 283); 

                auto tg_xyyy_xxxxzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 284); 

                auto tg_xyyy_xxxyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 285); 

                auto tg_xyyy_xxxyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 286); 

                auto tg_xyyy_xxxyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 287); 

                auto tg_xyyy_xxxyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 288); 

                auto tg_xyyy_xxxyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 289); 

                auto tg_xyyy_xxxzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 290); 

                // Batch of Integrals (243,291)

                #pragma omp simd aligned(fxn, fza, tg_xxzz_xxxyyzzz_0, tg_xxzz_xxxyzzzz_0, tg_xxzz_xxxzzzzz_0, \
                                         tg_xxzz_xxyyyyyy_0, tg_xxzz_xxyyyyyz_0, tg_xxzz_xxyyyyzz_0, tg_xxzz_xxyyyzzz_0, \
                                         tg_xxzz_xxyyzzzz_0, tg_xxzz_xxyzzzzz_0, tg_xxzz_xxzzzzzz_0, tg_xxzz_xyyyyyyy_0, \
                                         tg_xxzz_xyyyyyyz_0, tg_xxzz_xyyyyyzz_0, tg_xxzz_xyyyyzzz_0, tg_xxzz_xyyyzzzz_0, \
                                         tg_xxzz_xyyzzzzz_0, tg_xxzz_xyzzzzzz_0, tg_xxzz_xzzzzzzz_0, tg_xxzz_yyyyyyyy_0, \
                                         tg_xxzz_yyyyyyyz_0, tg_xxzz_yyyyyyzz_0, tg_xxzz_yyyyyzzz_0, tg_xxzz_yyyyzzzz_0, \
                                         tg_xxzz_yyyzzzzz_0, tg_xxzz_yyzzzzzz_0, tg_xxzz_yzzzzzzz_0, tg_xxzz_zzzzzzzz_0, \
                                         tg_xyyy_xxxxxxxx_0, tg_xyyy_xxxxxxxy_0, tg_xyyy_xxxxxxxz_0, tg_xyyy_xxxxxxyy_0, \
                                         tg_xyyy_xxxxxxyz_0, tg_xyyy_xxxxxxzz_0, tg_xyyy_xxxxxyyy_0, tg_xyyy_xxxxxyyz_0, \
                                         tg_xyyy_xxxxxyzz_0, tg_xyyy_xxxxxzzz_0, tg_xyyy_xxxxyyyy_0, tg_xyyy_xxxxyyyz_0, \
                                         tg_xyyy_xxxxyyzz_0, tg_xyyy_xxxxyzzz_0, tg_xyyy_xxxxzzzz_0, tg_xyyy_xxxyyyyy_0, \
                                         tg_xyyy_xxxyyyyz_0, tg_xyyy_xxxyyyzz_0, tg_xyyy_xxxyyzzz_0, tg_xyyy_xxxyzzzz_0, \
                                         tg_xyyy_xxxzzzzz_0, tg_xzz_xxxyyzzz_0, tg_xzz_xxxyyzzz_1, tg_xzz_xxxyzzzz_0, \
                                         tg_xzz_xxxyzzzz_1, tg_xzz_xxxzzzzz_0, tg_xzz_xxxzzzzz_1, tg_xzz_xxyyyyyy_0, \
                                         tg_xzz_xxyyyyyy_1, tg_xzz_xxyyyyyz_0, tg_xzz_xxyyyyyz_1, tg_xzz_xxyyyyzz_0, \
                                         tg_xzz_xxyyyyzz_1, tg_xzz_xxyyyzzz_0, tg_xzz_xxyyyzzz_1, tg_xzz_xxyyzzz_1, \
                                         tg_xzz_xxyyzzzz_0, tg_xzz_xxyyzzzz_1, tg_xzz_xxyzzzz_1, tg_xzz_xxyzzzzz_0, \
                                         tg_xzz_xxyzzzzz_1, tg_xzz_xxzzzzz_1, tg_xzz_xxzzzzzz_0, tg_xzz_xxzzzzzz_1, \
                                         tg_xzz_xyyyyyy_1, tg_xzz_xyyyyyyy_0, tg_xzz_xyyyyyyy_1, tg_xzz_xyyyyyyz_0, \
                                         tg_xzz_xyyyyyyz_1, tg_xzz_xyyyyyz_1, tg_xzz_xyyyyyzz_0, tg_xzz_xyyyyyzz_1, \
                                         tg_xzz_xyyyyzz_1, tg_xzz_xyyyyzzz_0, tg_xzz_xyyyyzzz_1, tg_xzz_xyyyzzz_1, \
                                         tg_xzz_xyyyzzzz_0, tg_xzz_xyyyzzzz_1, tg_xzz_xyyzzzz_1, tg_xzz_xyyzzzzz_0, \
                                         tg_xzz_xyyzzzzz_1, tg_xzz_xyzzzzz_1, tg_xzz_xyzzzzzz_0, tg_xzz_xyzzzzzz_1, \
                                         tg_xzz_xzzzzzz_1, tg_xzz_xzzzzzzz_0, tg_xzz_xzzzzzzz_1, tg_xzz_yyyyyyy_1, \
                                         tg_xzz_yyyyyyyy_0, tg_xzz_yyyyyyyy_1, tg_xzz_yyyyyyyz_0, tg_xzz_yyyyyyyz_1, \
                                         tg_xzz_yyyyyyz_1, tg_xzz_yyyyyyzz_0, tg_xzz_yyyyyyzz_1, tg_xzz_yyyyyzz_1, \
                                         tg_xzz_yyyyyzzz_0, tg_xzz_yyyyyzzz_1, tg_xzz_yyyyzzz_1, tg_xzz_yyyyzzzz_0, \
                                         tg_xzz_yyyyzzzz_1, tg_xzz_yyyzzzz_1, tg_xzz_yyyzzzzz_0, tg_xzz_yyyzzzzz_1, \
                                         tg_xzz_yyzzzzz_1, tg_xzz_yyzzzzzz_0, tg_xzz_yyzzzzzz_1, tg_xzz_yzzzzzz_1, \
                                         tg_xzz_yzzzzzzz_0, tg_xzz_yzzzzzzz_1, tg_xzz_zzzzzzz_1, tg_xzz_zzzzzzzz_0, \
                                         tg_xzz_zzzzzzzz_1, tg_yyy_xxxxxxx_1, tg_yyy_xxxxxxxx_0, tg_yyy_xxxxxxxx_1, \
                                         tg_yyy_xxxxxxxy_0, tg_yyy_xxxxxxxy_1, tg_yyy_xxxxxxxz_0, tg_yyy_xxxxxxxz_1, \
                                         tg_yyy_xxxxxxy_1, tg_yyy_xxxxxxyy_0, tg_yyy_xxxxxxyy_1, tg_yyy_xxxxxxyz_0, \
                                         tg_yyy_xxxxxxyz_1, tg_yyy_xxxxxxz_1, tg_yyy_xxxxxxzz_0, tg_yyy_xxxxxxzz_1, \
                                         tg_yyy_xxxxxyy_1, tg_yyy_xxxxxyyy_0, tg_yyy_xxxxxyyy_1, tg_yyy_xxxxxyyz_0, \
                                         tg_yyy_xxxxxyyz_1, tg_yyy_xxxxxyz_1, tg_yyy_xxxxxyzz_0, tg_yyy_xxxxxyzz_1, \
                                         tg_yyy_xxxxxzz_1, tg_yyy_xxxxxzzz_0, tg_yyy_xxxxxzzz_1, tg_yyy_xxxxyyy_1, \
                                         tg_yyy_xxxxyyyy_0, tg_yyy_xxxxyyyy_1, tg_yyy_xxxxyyyz_0, tg_yyy_xxxxyyyz_1, \
                                         tg_yyy_xxxxyyz_1, tg_yyy_xxxxyyzz_0, tg_yyy_xxxxyyzz_1, tg_yyy_xxxxyzz_1, \
                                         tg_yyy_xxxxyzzz_0, tg_yyy_xxxxyzzz_1, tg_yyy_xxxxzzz_1, tg_yyy_xxxxzzzz_0, \
                                         tg_yyy_xxxxzzzz_1, tg_yyy_xxxyyyy_1, tg_yyy_xxxyyyyy_0, tg_yyy_xxxyyyyy_1, \
                                         tg_yyy_xxxyyyyz_0, tg_yyy_xxxyyyyz_1, tg_yyy_xxxyyyz_1, tg_yyy_xxxyyyzz_0, \
                                         tg_yyy_xxxyyyzz_1, tg_yyy_xxxyyzz_1, tg_yyy_xxxyyzzz_0, tg_yyy_xxxyyzzz_1, \
                                         tg_yyy_xxxyzzz_1, tg_yyy_xxxyzzzz_0, tg_yyy_xxxyzzzz_1, tg_yyy_xxxzzzz_1, \
                                         tg_yyy_xxxzzzzz_0, tg_yyy_xxxzzzzz_1, tg_yyy_xxyyyyy_1, tg_yyy_xxyyyyz_1, \
                                         tg_yyy_xxyyyzz_1, tg_yyy_xxyyzzz_1, tg_yyy_xxyzzzz_1, tg_yyy_xxzzzzz_1, \
                                         tg_zz_xxxyyzzz_0, tg_zz_xxxyyzzz_1, tg_zz_xxxyzzzz_0, tg_zz_xxxyzzzz_1, \
                                         tg_zz_xxxzzzzz_0, tg_zz_xxxzzzzz_1, tg_zz_xxyyyyyy_0, tg_zz_xxyyyyyy_1, \
                                         tg_zz_xxyyyyyz_0, tg_zz_xxyyyyyz_1, tg_zz_xxyyyyzz_0, tg_zz_xxyyyyzz_1, \
                                         tg_zz_xxyyyzzz_0, tg_zz_xxyyyzzz_1, tg_zz_xxyyzzzz_0, tg_zz_xxyyzzzz_1, \
                                         tg_zz_xxyzzzzz_0, tg_zz_xxyzzzzz_1, tg_zz_xxzzzzzz_0, tg_zz_xxzzzzzz_1, \
                                         tg_zz_xyyyyyyy_0, tg_zz_xyyyyyyy_1, tg_zz_xyyyyyyz_0, tg_zz_xyyyyyyz_1, \
                                         tg_zz_xyyyyyzz_0, tg_zz_xyyyyyzz_1, tg_zz_xyyyyzzz_0, tg_zz_xyyyyzzz_1, \
                                         tg_zz_xyyyzzzz_0, tg_zz_xyyyzzzz_1, tg_zz_xyyzzzzz_0, tg_zz_xyyzzzzz_1, \
                                         tg_zz_xyzzzzzz_0, tg_zz_xyzzzzzz_1, tg_zz_xzzzzzzz_0, tg_zz_xzzzzzzz_1, \
                                         tg_zz_yyyyyyyy_0, tg_zz_yyyyyyyy_1, tg_zz_yyyyyyyz_0, tg_zz_yyyyyyyz_1, \
                                         tg_zz_yyyyyyzz_0, tg_zz_yyyyyyzz_1, tg_zz_yyyyyzzz_0, tg_zz_yyyyyzzz_1, \
                                         tg_zz_yyyyzzzz_0, tg_zz_yyyyzzzz_1, tg_zz_yyyzzzzz_0, tg_zz_yyyzzzzz_1, \
                                         tg_zz_yyzzzzzz_0, tg_zz_yyzzzzzz_1, tg_zz_yzzzzzzz_0, tg_zz_yzzzzzzz_1, \
                                         tg_zz_zzzzzzzz_0, tg_zz_zzzzzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxzz_xxxyyzzz_0[j] = pb_x * tg_xzz_xxxyyzzz_0[j] + wp_x[j] * tg_xzz_xxxyyzzz_1[j] + 0.5 * fl1_fx * tg_zz_xxxyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxyyzzz_1[j] + 1.5 * fl1_fxn * tg_xzz_xxyyzzz_1[j];

                    tg_xxzz_xxxyzzzz_0[j] = pb_x * tg_xzz_xxxyzzzz_0[j] + wp_x[j] * tg_xzz_xxxyzzzz_1[j] + 0.5 * fl1_fx * tg_zz_xxxyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxyzzzz_1[j] + 1.5 * fl1_fxn * tg_xzz_xxyzzzz_1[j];

                    tg_xxzz_xxxzzzzz_0[j] = pb_x * tg_xzz_xxxzzzzz_0[j] + wp_x[j] * tg_xzz_xxxzzzzz_1[j] + 0.5 * fl1_fx * tg_zz_xxxzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxzzzzz_1[j] + 1.5 * fl1_fxn * tg_xzz_xxzzzzz_1[j];

                    tg_xxzz_xxyyyyyy_0[j] = pb_x * tg_xzz_xxyyyyyy_0[j] + wp_x[j] * tg_xzz_xxyyyyyy_1[j] + 0.5 * fl1_fx * tg_zz_xxyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxyyyyyy_1[j] + fl1_fxn * tg_xzz_xyyyyyy_1[j];

                    tg_xxzz_xxyyyyyz_0[j] = pb_x * tg_xzz_xxyyyyyz_0[j] + wp_x[j] * tg_xzz_xxyyyyyz_1[j] + 0.5 * fl1_fx * tg_zz_xxyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxyyyyyz_1[j] + fl1_fxn * tg_xzz_xyyyyyz_1[j];

                    tg_xxzz_xxyyyyzz_0[j] = pb_x * tg_xzz_xxyyyyzz_0[j] + wp_x[j] * tg_xzz_xxyyyyzz_1[j] + 0.5 * fl1_fx * tg_zz_xxyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxyyyyzz_1[j] + fl1_fxn * tg_xzz_xyyyyzz_1[j];

                    tg_xxzz_xxyyyzzz_0[j] = pb_x * tg_xzz_xxyyyzzz_0[j] + wp_x[j] * tg_xzz_xxyyyzzz_1[j] + 0.5 * fl1_fx * tg_zz_xxyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxyyyzzz_1[j] + fl1_fxn * tg_xzz_xyyyzzz_1[j];

                    tg_xxzz_xxyyzzzz_0[j] = pb_x * tg_xzz_xxyyzzzz_0[j] + wp_x[j] * tg_xzz_xxyyzzzz_1[j] + 0.5 * fl1_fx * tg_zz_xxyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxyyzzzz_1[j] + fl1_fxn * tg_xzz_xyyzzzz_1[j];

                    tg_xxzz_xxyzzzzz_0[j] = pb_x * tg_xzz_xxyzzzzz_0[j] + wp_x[j] * tg_xzz_xxyzzzzz_1[j] + 0.5 * fl1_fx * tg_zz_xxyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxyzzzzz_1[j] + fl1_fxn * tg_xzz_xyzzzzz_1[j];

                    tg_xxzz_xxzzzzzz_0[j] = pb_x * tg_xzz_xxzzzzzz_0[j] + wp_x[j] * tg_xzz_xxzzzzzz_1[j] + 0.5 * fl1_fx * tg_zz_xxzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxzzzzzz_1[j] + fl1_fxn * tg_xzz_xzzzzzz_1[j];

                    tg_xxzz_xyyyyyyy_0[j] = pb_x * tg_xzz_xyyyyyyy_0[j] + wp_x[j] * tg_xzz_xyyyyyyy_1[j] + 0.5 * fl1_fx * tg_zz_xyyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xyyyyyyy_1[j] + 0.5 * fl1_fxn * tg_xzz_yyyyyyy_1[j];

                    tg_xxzz_xyyyyyyz_0[j] = pb_x * tg_xzz_xyyyyyyz_0[j] + wp_x[j] * tg_xzz_xyyyyyyz_1[j] + 0.5 * fl1_fx * tg_zz_xyyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xyyyyyyz_1[j] + 0.5 * fl1_fxn * tg_xzz_yyyyyyz_1[j];

                    tg_xxzz_xyyyyyzz_0[j] = pb_x * tg_xzz_xyyyyyzz_0[j] + wp_x[j] * tg_xzz_xyyyyyzz_1[j] + 0.5 * fl1_fx * tg_zz_xyyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xyyyyyzz_1[j] + 0.5 * fl1_fxn * tg_xzz_yyyyyzz_1[j];

                    tg_xxzz_xyyyyzzz_0[j] = pb_x * tg_xzz_xyyyyzzz_0[j] + wp_x[j] * tg_xzz_xyyyyzzz_1[j] + 0.5 * fl1_fx * tg_zz_xyyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xyyyyzzz_1[j] + 0.5 * fl1_fxn * tg_xzz_yyyyzzz_1[j];

                    tg_xxzz_xyyyzzzz_0[j] = pb_x * tg_xzz_xyyyzzzz_0[j] + wp_x[j] * tg_xzz_xyyyzzzz_1[j] + 0.5 * fl1_fx * tg_zz_xyyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xyyyzzzz_1[j] + 0.5 * fl1_fxn * tg_xzz_yyyzzzz_1[j];

                    tg_xxzz_xyyzzzzz_0[j] = pb_x * tg_xzz_xyyzzzzz_0[j] + wp_x[j] * tg_xzz_xyyzzzzz_1[j] + 0.5 * fl1_fx * tg_zz_xyyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xyyzzzzz_1[j] + 0.5 * fl1_fxn * tg_xzz_yyzzzzz_1[j];

                    tg_xxzz_xyzzzzzz_0[j] = pb_x * tg_xzz_xyzzzzzz_0[j] + wp_x[j] * tg_xzz_xyzzzzzz_1[j] + 0.5 * fl1_fx * tg_zz_xyzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xzz_yzzzzzz_1[j];

                    tg_xxzz_xzzzzzzz_0[j] = pb_x * tg_xzz_xzzzzzzz_0[j] + wp_x[j] * tg_xzz_xzzzzzzz_1[j] + 0.5 * fl1_fx * tg_zz_xzzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_xzz_zzzzzzz_1[j];

                    tg_xxzz_yyyyyyyy_0[j] = pb_x * tg_xzz_yyyyyyyy_0[j] + wp_x[j] * tg_xzz_yyyyyyyy_1[j] + 0.5 * fl1_fx * tg_zz_yyyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_yyyyyyyy_1[j];

                    tg_xxzz_yyyyyyyz_0[j] = pb_x * tg_xzz_yyyyyyyz_0[j] + wp_x[j] * tg_xzz_yyyyyyyz_1[j] + 0.5 * fl1_fx * tg_zz_yyyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_yyyyyyyz_1[j];

                    tg_xxzz_yyyyyyzz_0[j] = pb_x * tg_xzz_yyyyyyzz_0[j] + wp_x[j] * tg_xzz_yyyyyyzz_1[j] + 0.5 * fl1_fx * tg_zz_yyyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_yyyyyyzz_1[j];

                    tg_xxzz_yyyyyzzz_0[j] = pb_x * tg_xzz_yyyyyzzz_0[j] + wp_x[j] * tg_xzz_yyyyyzzz_1[j] + 0.5 * fl1_fx * tg_zz_yyyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_yyyyyzzz_1[j];

                    tg_xxzz_yyyyzzzz_0[j] = pb_x * tg_xzz_yyyyzzzz_0[j] + wp_x[j] * tg_xzz_yyyyzzzz_1[j] + 0.5 * fl1_fx * tg_zz_yyyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_yyyyzzzz_1[j];

                    tg_xxzz_yyyzzzzz_0[j] = pb_x * tg_xzz_yyyzzzzz_0[j] + wp_x[j] * tg_xzz_yyyzzzzz_1[j] + 0.5 * fl1_fx * tg_zz_yyyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_yyyzzzzz_1[j];

                    tg_xxzz_yyzzzzzz_0[j] = pb_x * tg_xzz_yyzzzzzz_0[j] + wp_x[j] * tg_xzz_yyzzzzzz_1[j] + 0.5 * fl1_fx * tg_zz_yyzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_yyzzzzzz_1[j];

                    tg_xxzz_yzzzzzzz_0[j] = pb_x * tg_xzz_yzzzzzzz_0[j] + wp_x[j] * tg_xzz_yzzzzzzz_1[j] + 0.5 * fl1_fx * tg_zz_yzzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_yzzzzzzz_1[j];

                    tg_xxzz_zzzzzzzz_0[j] = pb_x * tg_xzz_zzzzzzzz_0[j] + wp_x[j] * tg_xzz_zzzzzzzz_1[j] + 0.5 * fl1_fx * tg_zz_zzzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_zzzzzzzz_1[j];

                    tg_xyyy_xxxxxxxx_0[j] = pb_x * tg_yyy_xxxxxxxx_0[j] + wp_x[j] * tg_yyy_xxxxxxxx_1[j] + 4.0 * fl1_fxn * tg_yyy_xxxxxxx_1[j];

                    tg_xyyy_xxxxxxxy_0[j] = pb_x * tg_yyy_xxxxxxxy_0[j] + wp_x[j] * tg_yyy_xxxxxxxy_1[j] + 3.5 * fl1_fxn * tg_yyy_xxxxxxy_1[j];

                    tg_xyyy_xxxxxxxz_0[j] = pb_x * tg_yyy_xxxxxxxz_0[j] + wp_x[j] * tg_yyy_xxxxxxxz_1[j] + 3.5 * fl1_fxn * tg_yyy_xxxxxxz_1[j];

                    tg_xyyy_xxxxxxyy_0[j] = pb_x * tg_yyy_xxxxxxyy_0[j] + wp_x[j] * tg_yyy_xxxxxxyy_1[j] + 3.0 * fl1_fxn * tg_yyy_xxxxxyy_1[j];

                    tg_xyyy_xxxxxxyz_0[j] = pb_x * tg_yyy_xxxxxxyz_0[j] + wp_x[j] * tg_yyy_xxxxxxyz_1[j] + 3.0 * fl1_fxn * tg_yyy_xxxxxyz_1[j];

                    tg_xyyy_xxxxxxzz_0[j] = pb_x * tg_yyy_xxxxxxzz_0[j] + wp_x[j] * tg_yyy_xxxxxxzz_1[j] + 3.0 * fl1_fxn * tg_yyy_xxxxxzz_1[j];

                    tg_xyyy_xxxxxyyy_0[j] = pb_x * tg_yyy_xxxxxyyy_0[j] + wp_x[j] * tg_yyy_xxxxxyyy_1[j] + 2.5 * fl1_fxn * tg_yyy_xxxxyyy_1[j];

                    tg_xyyy_xxxxxyyz_0[j] = pb_x * tg_yyy_xxxxxyyz_0[j] + wp_x[j] * tg_yyy_xxxxxyyz_1[j] + 2.5 * fl1_fxn * tg_yyy_xxxxyyz_1[j];

                    tg_xyyy_xxxxxyzz_0[j] = pb_x * tg_yyy_xxxxxyzz_0[j] + wp_x[j] * tg_yyy_xxxxxyzz_1[j] + 2.5 * fl1_fxn * tg_yyy_xxxxyzz_1[j];

                    tg_xyyy_xxxxxzzz_0[j] = pb_x * tg_yyy_xxxxxzzz_0[j] + wp_x[j] * tg_yyy_xxxxxzzz_1[j] + 2.5 * fl1_fxn * tg_yyy_xxxxzzz_1[j];

                    tg_xyyy_xxxxyyyy_0[j] = pb_x * tg_yyy_xxxxyyyy_0[j] + wp_x[j] * tg_yyy_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_yyy_xxxyyyy_1[j];

                    tg_xyyy_xxxxyyyz_0[j] = pb_x * tg_yyy_xxxxyyyz_0[j] + wp_x[j] * tg_yyy_xxxxyyyz_1[j] + 2.0 * fl1_fxn * tg_yyy_xxxyyyz_1[j];

                    tg_xyyy_xxxxyyzz_0[j] = pb_x * tg_yyy_xxxxyyzz_0[j] + wp_x[j] * tg_yyy_xxxxyyzz_1[j] + 2.0 * fl1_fxn * tg_yyy_xxxyyzz_1[j];

                    tg_xyyy_xxxxyzzz_0[j] = pb_x * tg_yyy_xxxxyzzz_0[j] + wp_x[j] * tg_yyy_xxxxyzzz_1[j] + 2.0 * fl1_fxn * tg_yyy_xxxyzzz_1[j];

                    tg_xyyy_xxxxzzzz_0[j] = pb_x * tg_yyy_xxxxzzzz_0[j] + wp_x[j] * tg_yyy_xxxxzzzz_1[j] + 2.0 * fl1_fxn * tg_yyy_xxxzzzz_1[j];

                    tg_xyyy_xxxyyyyy_0[j] = pb_x * tg_yyy_xxxyyyyy_0[j] + wp_x[j] * tg_yyy_xxxyyyyy_1[j] + 1.5 * fl1_fxn * tg_yyy_xxyyyyy_1[j];

                    tg_xyyy_xxxyyyyz_0[j] = pb_x * tg_yyy_xxxyyyyz_0[j] + wp_x[j] * tg_yyy_xxxyyyyz_1[j] + 1.5 * fl1_fxn * tg_yyy_xxyyyyz_1[j];

                    tg_xyyy_xxxyyyzz_0[j] = pb_x * tg_yyy_xxxyyyzz_0[j] + wp_x[j] * tg_yyy_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_yyy_xxyyyzz_1[j];

                    tg_xyyy_xxxyyzzz_0[j] = pb_x * tg_yyy_xxxyyzzz_0[j] + wp_x[j] * tg_yyy_xxxyyzzz_1[j] + 1.5 * fl1_fxn * tg_yyy_xxyyzzz_1[j];

                    tg_xyyy_xxxyzzzz_0[j] = pb_x * tg_yyy_xxxyzzzz_0[j] + wp_x[j] * tg_yyy_xxxyzzzz_1[j] + 1.5 * fl1_fxn * tg_yyy_xxyzzzz_1[j];

                    tg_xyyy_xxxzzzzz_0[j] = pb_x * tg_yyy_xxxzzzzz_0[j] + wp_x[j] * tg_yyy_xxxzzzzz_1[j] + 1.5 * fl1_fxn * tg_yyy_xxzzzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSL_291_339(      CMemBlock2D<double>& primBuffer,
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

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {4, -1, -1, -1},
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_3_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {7, -1, -1, -1}, 
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

                auto tg_yyy_xxyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 291); 

                auto tg_yyy_xxyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 292); 

                auto tg_yyy_xxyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 293); 

                auto tg_yyy_xxyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 294); 

                auto tg_yyy_xxyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 295); 

                auto tg_yyy_xxyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 296); 

                auto tg_yyy_xxzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 297); 

                auto tg_yyy_xyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 298); 

                auto tg_yyy_xyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 299); 

                auto tg_yyy_xyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 300); 

                auto tg_yyy_xyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 301); 

                auto tg_yyy_xyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 302); 

                auto tg_yyy_xyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 303); 

                auto tg_yyy_xyzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 304); 

                auto tg_yyy_xzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 305); 

                auto tg_yyy_yyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 306); 

                auto tg_yyy_yyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 307); 

                auto tg_yyy_yyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 308); 

                auto tg_yyy_yyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 309); 

                auto tg_yyy_yyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 310); 

                auto tg_yyy_yyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 311); 

                auto tg_yyy_yyzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 312); 

                auto tg_yyy_yzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 313); 

                auto tg_yyy_zzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 314); 

                auto tg_yyz_xxxxxxxx_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 315); 

                auto tg_yyz_xxxxxxxy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 316); 

                auto tg_yyz_xxxxxxxz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 317); 

                auto tg_yyz_xxxxxxyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 318); 

                auto tg_yyz_xxxxxxyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 319); 

                auto tg_yyz_xxxxxxzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 320); 

                auto tg_yyz_xxxxxyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 321); 

                auto tg_yyz_xxxxxyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 322); 

                auto tg_yyz_xxxxxyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 323); 

                auto tg_yyz_xxxxxzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 324); 

                auto tg_yyz_xxxxyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 325); 

                auto tg_yyz_xxxxyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 326); 

                auto tg_yyz_xxxxyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 327); 

                auto tg_yyz_xxxxyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 328); 

                auto tg_yyz_xxxxzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 329); 

                auto tg_yyz_xxxyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 330); 

                auto tg_yyz_xxxyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 331); 

                auto tg_yyz_xxxyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 332); 

                auto tg_yyz_xxxyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 333); 

                auto tg_yyz_xxxyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 334); 

                auto tg_yyz_xxxzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 335); 

                auto tg_yyz_xxyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 336); 

                auto tg_yyz_xxyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 337); 

                auto tg_yyz_xxyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 338); 

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

                // set up pointers to integrals

                auto tg_xyyy_xxyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 291); 

                auto tg_xyyy_xxyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 292); 

                auto tg_xyyy_xxyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 293); 

                auto tg_xyyy_xxyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 294); 

                auto tg_xyyy_xxyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 295); 

                auto tg_xyyy_xxyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 296); 

                auto tg_xyyy_xxzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 297); 

                auto tg_xyyy_xyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 298); 

                auto tg_xyyy_xyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 299); 

                auto tg_xyyy_xyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 300); 

                auto tg_xyyy_xyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 301); 

                auto tg_xyyy_xyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 302); 

                auto tg_xyyy_xyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 303); 

                auto tg_xyyy_xyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 304); 

                auto tg_xyyy_xzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 305); 

                auto tg_xyyy_yyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 306); 

                auto tg_xyyy_yyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 307); 

                auto tg_xyyy_yyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 308); 

                auto tg_xyyy_yyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 309); 

                auto tg_xyyy_yyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 310); 

                auto tg_xyyy_yyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 311); 

                auto tg_xyyy_yyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 312); 

                auto tg_xyyy_yzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 313); 

                auto tg_xyyy_zzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 314); 

                auto tg_xyyz_xxxxxxxx_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 315); 

                auto tg_xyyz_xxxxxxxy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 316); 

                auto tg_xyyz_xxxxxxxz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 317); 

                auto tg_xyyz_xxxxxxyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 318); 

                auto tg_xyyz_xxxxxxyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 319); 

                auto tg_xyyz_xxxxxxzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 320); 

                auto tg_xyyz_xxxxxyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 321); 

                auto tg_xyyz_xxxxxyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 322); 

                auto tg_xyyz_xxxxxyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 323); 

                auto tg_xyyz_xxxxxzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 324); 

                auto tg_xyyz_xxxxyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 325); 

                auto tg_xyyz_xxxxyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 326); 

                auto tg_xyyz_xxxxyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 327); 

                auto tg_xyyz_xxxxyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 328); 

                auto tg_xyyz_xxxxzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 329); 

                auto tg_xyyz_xxxyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 330); 

                auto tg_xyyz_xxxyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 331); 

                auto tg_xyyz_xxxyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 332); 

                auto tg_xyyz_xxxyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 333); 

                auto tg_xyyz_xxxyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 334); 

                auto tg_xyyz_xxxzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 335); 

                auto tg_xyyz_xxyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 336); 

                auto tg_xyyz_xxyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 337); 

                auto tg_xyyz_xxyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 338); 

                // Batch of Integrals (291,339)

                #pragma omp simd aligned(fxn, tg_xyyy_xxyyyyyy_0, tg_xyyy_xxyyyyyz_0, tg_xyyy_xxyyyyzz_0, \
                                         tg_xyyy_xxyyyzzz_0, tg_xyyy_xxyyzzzz_0, tg_xyyy_xxyzzzzz_0, tg_xyyy_xxzzzzzz_0, \
                                         tg_xyyy_xyyyyyyy_0, tg_xyyy_xyyyyyyz_0, tg_xyyy_xyyyyyzz_0, tg_xyyy_xyyyyzzz_0, \
                                         tg_xyyy_xyyyzzzz_0, tg_xyyy_xyyzzzzz_0, tg_xyyy_xyzzzzzz_0, tg_xyyy_xzzzzzzz_0, \
                                         tg_xyyy_yyyyyyyy_0, tg_xyyy_yyyyyyyz_0, tg_xyyy_yyyyyyzz_0, tg_xyyy_yyyyyzzz_0, \
                                         tg_xyyy_yyyyzzzz_0, tg_xyyy_yyyzzzzz_0, tg_xyyy_yyzzzzzz_0, tg_xyyy_yzzzzzzz_0, \
                                         tg_xyyy_zzzzzzzz_0, tg_xyyz_xxxxxxxx_0, tg_xyyz_xxxxxxxy_0, tg_xyyz_xxxxxxxz_0, \
                                         tg_xyyz_xxxxxxyy_0, tg_xyyz_xxxxxxyz_0, tg_xyyz_xxxxxxzz_0, tg_xyyz_xxxxxyyy_0, \
                                         tg_xyyz_xxxxxyyz_0, tg_xyyz_xxxxxyzz_0, tg_xyyz_xxxxxzzz_0, tg_xyyz_xxxxyyyy_0, \
                                         tg_xyyz_xxxxyyyz_0, tg_xyyz_xxxxyyzz_0, tg_xyyz_xxxxyzzz_0, tg_xyyz_xxxxzzzz_0, \
                                         tg_xyyz_xxxyyyyy_0, tg_xyyz_xxxyyyyz_0, tg_xyyz_xxxyyyzz_0, tg_xyyz_xxxyyzzz_0, \
                                         tg_xyyz_xxxyzzzz_0, tg_xyyz_xxxzzzzz_0, tg_xyyz_xxyyyyyy_0, tg_xyyz_xxyyyyyz_0, \
                                         tg_xyyz_xxyyyyzz_0, tg_yyy_xxyyyyyy_0, tg_yyy_xxyyyyyy_1, tg_yyy_xxyyyyyz_0, \
                                         tg_yyy_xxyyyyyz_1, tg_yyy_xxyyyyzz_0, tg_yyy_xxyyyyzz_1, tg_yyy_xxyyyzzz_0, \
                                         tg_yyy_xxyyyzzz_1, tg_yyy_xxyyzzzz_0, tg_yyy_xxyyzzzz_1, tg_yyy_xxyzzzzz_0, \
                                         tg_yyy_xxyzzzzz_1, tg_yyy_xxzzzzzz_0, tg_yyy_xxzzzzzz_1, tg_yyy_xyyyyyy_1, \
                                         tg_yyy_xyyyyyyy_0, tg_yyy_xyyyyyyy_1, tg_yyy_xyyyyyyz_0, tg_yyy_xyyyyyyz_1, \
                                         tg_yyy_xyyyyyz_1, tg_yyy_xyyyyyzz_0, tg_yyy_xyyyyyzz_1, tg_yyy_xyyyyzz_1, \
                                         tg_yyy_xyyyyzzz_0, tg_yyy_xyyyyzzz_1, tg_yyy_xyyyzzz_1, tg_yyy_xyyyzzzz_0, \
                                         tg_yyy_xyyyzzzz_1, tg_yyy_xyyzzzz_1, tg_yyy_xyyzzzzz_0, tg_yyy_xyyzzzzz_1, \
                                         tg_yyy_xyzzzzz_1, tg_yyy_xyzzzzzz_0, tg_yyy_xyzzzzzz_1, tg_yyy_xzzzzzz_1, \
                                         tg_yyy_xzzzzzzz_0, tg_yyy_xzzzzzzz_1, tg_yyy_yyyyyyy_1, tg_yyy_yyyyyyyy_0, \
                                         tg_yyy_yyyyyyyy_1, tg_yyy_yyyyyyyz_0, tg_yyy_yyyyyyyz_1, tg_yyy_yyyyyyz_1, \
                                         tg_yyy_yyyyyyzz_0, tg_yyy_yyyyyyzz_1, tg_yyy_yyyyyzz_1, tg_yyy_yyyyyzzz_0, \
                                         tg_yyy_yyyyyzzz_1, tg_yyy_yyyyzzz_1, tg_yyy_yyyyzzzz_0, tg_yyy_yyyyzzzz_1, \
                                         tg_yyy_yyyzzzz_1, tg_yyy_yyyzzzzz_0, tg_yyy_yyyzzzzz_1, tg_yyy_yyzzzzz_1, \
                                         tg_yyy_yyzzzzzz_0, tg_yyy_yyzzzzzz_1, tg_yyy_yzzzzzz_1, tg_yyy_yzzzzzzz_0, \
                                         tg_yyy_yzzzzzzz_1, tg_yyy_zzzzzzz_1, tg_yyy_zzzzzzzz_0, tg_yyy_zzzzzzzz_1, \
                                         tg_yyz_xxxxxxx_1, tg_yyz_xxxxxxxx_0, tg_yyz_xxxxxxxx_1, tg_yyz_xxxxxxxy_0, \
                                         tg_yyz_xxxxxxxy_1, tg_yyz_xxxxxxxz_0, tg_yyz_xxxxxxxz_1, tg_yyz_xxxxxxy_1, \
                                         tg_yyz_xxxxxxyy_0, tg_yyz_xxxxxxyy_1, tg_yyz_xxxxxxyz_0, tg_yyz_xxxxxxyz_1, \
                                         tg_yyz_xxxxxxz_1, tg_yyz_xxxxxxzz_0, tg_yyz_xxxxxxzz_1, tg_yyz_xxxxxyy_1, \
                                         tg_yyz_xxxxxyyy_0, tg_yyz_xxxxxyyy_1, tg_yyz_xxxxxyyz_0, tg_yyz_xxxxxyyz_1, \
                                         tg_yyz_xxxxxyz_1, tg_yyz_xxxxxyzz_0, tg_yyz_xxxxxyzz_1, tg_yyz_xxxxxzz_1, \
                                         tg_yyz_xxxxxzzz_0, tg_yyz_xxxxxzzz_1, tg_yyz_xxxxyyy_1, tg_yyz_xxxxyyyy_0, \
                                         tg_yyz_xxxxyyyy_1, tg_yyz_xxxxyyyz_0, tg_yyz_xxxxyyyz_1, tg_yyz_xxxxyyz_1, \
                                         tg_yyz_xxxxyyzz_0, tg_yyz_xxxxyyzz_1, tg_yyz_xxxxyzz_1, tg_yyz_xxxxyzzz_0, \
                                         tg_yyz_xxxxyzzz_1, tg_yyz_xxxxzzz_1, tg_yyz_xxxxzzzz_0, tg_yyz_xxxxzzzz_1, \
                                         tg_yyz_xxxyyyy_1, tg_yyz_xxxyyyyy_0, tg_yyz_xxxyyyyy_1, tg_yyz_xxxyyyyz_0, \
                                         tg_yyz_xxxyyyyz_1, tg_yyz_xxxyyyz_1, tg_yyz_xxxyyyzz_0, tg_yyz_xxxyyyzz_1, \
                                         tg_yyz_xxxyyzz_1, tg_yyz_xxxyyzzz_0, tg_yyz_xxxyyzzz_1, tg_yyz_xxxyzzz_1, \
                                         tg_yyz_xxxyzzzz_0, tg_yyz_xxxyzzzz_1, tg_yyz_xxxzzzz_1, tg_yyz_xxxzzzzz_0, \
                                         tg_yyz_xxxzzzzz_1, tg_yyz_xxyyyyy_1, tg_yyz_xxyyyyyy_0, tg_yyz_xxyyyyyy_1, \
                                         tg_yyz_xxyyyyyz_0, tg_yyz_xxyyyyyz_1, tg_yyz_xxyyyyz_1, tg_yyz_xxyyyyzz_0, \
                                         tg_yyz_xxyyyyzz_1, tg_yyz_xxyyyzz_1, tg_yyz_xxyyzzz_1, tg_yyz_xxyzzzz_1, \
                                         tg_yyz_xxzzzzz_1, tg_yyz_xyyyyyy_1, tg_yyz_xyyyyyz_1, tg_yyz_xyyyyzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fxn = fxn[j];

                    tg_xyyy_xxyyyyyy_0[j] = pb_x * tg_yyy_xxyyyyyy_0[j] + wp_x[j] * tg_yyy_xxyyyyyy_1[j] + fl1_fxn * tg_yyy_xyyyyyy_1[j];

                    tg_xyyy_xxyyyyyz_0[j] = pb_x * tg_yyy_xxyyyyyz_0[j] + wp_x[j] * tg_yyy_xxyyyyyz_1[j] + fl1_fxn * tg_yyy_xyyyyyz_1[j];

                    tg_xyyy_xxyyyyzz_0[j] = pb_x * tg_yyy_xxyyyyzz_0[j] + wp_x[j] * tg_yyy_xxyyyyzz_1[j] + fl1_fxn * tg_yyy_xyyyyzz_1[j];

                    tg_xyyy_xxyyyzzz_0[j] = pb_x * tg_yyy_xxyyyzzz_0[j] + wp_x[j] * tg_yyy_xxyyyzzz_1[j] + fl1_fxn * tg_yyy_xyyyzzz_1[j];

                    tg_xyyy_xxyyzzzz_0[j] = pb_x * tg_yyy_xxyyzzzz_0[j] + wp_x[j] * tg_yyy_xxyyzzzz_1[j] + fl1_fxn * tg_yyy_xyyzzzz_1[j];

                    tg_xyyy_xxyzzzzz_0[j] = pb_x * tg_yyy_xxyzzzzz_0[j] + wp_x[j] * tg_yyy_xxyzzzzz_1[j] + fl1_fxn * tg_yyy_xyzzzzz_1[j];

                    tg_xyyy_xxzzzzzz_0[j] = pb_x * tg_yyy_xxzzzzzz_0[j] + wp_x[j] * tg_yyy_xxzzzzzz_1[j] + fl1_fxn * tg_yyy_xzzzzzz_1[j];

                    tg_xyyy_xyyyyyyy_0[j] = pb_x * tg_yyy_xyyyyyyy_0[j] + wp_x[j] * tg_yyy_xyyyyyyy_1[j] + 0.5 * fl1_fxn * tg_yyy_yyyyyyy_1[j];

                    tg_xyyy_xyyyyyyz_0[j] = pb_x * tg_yyy_xyyyyyyz_0[j] + wp_x[j] * tg_yyy_xyyyyyyz_1[j] + 0.5 * fl1_fxn * tg_yyy_yyyyyyz_1[j];

                    tg_xyyy_xyyyyyzz_0[j] = pb_x * tg_yyy_xyyyyyzz_0[j] + wp_x[j] * tg_yyy_xyyyyyzz_1[j] + 0.5 * fl1_fxn * tg_yyy_yyyyyzz_1[j];

                    tg_xyyy_xyyyyzzz_0[j] = pb_x * tg_yyy_xyyyyzzz_0[j] + wp_x[j] * tg_yyy_xyyyyzzz_1[j] + 0.5 * fl1_fxn * tg_yyy_yyyyzzz_1[j];

                    tg_xyyy_xyyyzzzz_0[j] = pb_x * tg_yyy_xyyyzzzz_0[j] + wp_x[j] * tg_yyy_xyyyzzzz_1[j] + 0.5 * fl1_fxn * tg_yyy_yyyzzzz_1[j];

                    tg_xyyy_xyyzzzzz_0[j] = pb_x * tg_yyy_xyyzzzzz_0[j] + wp_x[j] * tg_yyy_xyyzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyy_yyzzzzz_1[j];

                    tg_xyyy_xyzzzzzz_0[j] = pb_x * tg_yyy_xyzzzzzz_0[j] + wp_x[j] * tg_yyy_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyy_yzzzzzz_1[j];

                    tg_xyyy_xzzzzzzz_0[j] = pb_x * tg_yyy_xzzzzzzz_0[j] + wp_x[j] * tg_yyy_xzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyy_zzzzzzz_1[j];

                    tg_xyyy_yyyyyyyy_0[j] = pb_x * tg_yyy_yyyyyyyy_0[j] + wp_x[j] * tg_yyy_yyyyyyyy_1[j];

                    tg_xyyy_yyyyyyyz_0[j] = pb_x * tg_yyy_yyyyyyyz_0[j] + wp_x[j] * tg_yyy_yyyyyyyz_1[j];

                    tg_xyyy_yyyyyyzz_0[j] = pb_x * tg_yyy_yyyyyyzz_0[j] + wp_x[j] * tg_yyy_yyyyyyzz_1[j];

                    tg_xyyy_yyyyyzzz_0[j] = pb_x * tg_yyy_yyyyyzzz_0[j] + wp_x[j] * tg_yyy_yyyyyzzz_1[j];

                    tg_xyyy_yyyyzzzz_0[j] = pb_x * tg_yyy_yyyyzzzz_0[j] + wp_x[j] * tg_yyy_yyyyzzzz_1[j];

                    tg_xyyy_yyyzzzzz_0[j] = pb_x * tg_yyy_yyyzzzzz_0[j] + wp_x[j] * tg_yyy_yyyzzzzz_1[j];

                    tg_xyyy_yyzzzzzz_0[j] = pb_x * tg_yyy_yyzzzzzz_0[j] + wp_x[j] * tg_yyy_yyzzzzzz_1[j];

                    tg_xyyy_yzzzzzzz_0[j] = pb_x * tg_yyy_yzzzzzzz_0[j] + wp_x[j] * tg_yyy_yzzzzzzz_1[j];

                    tg_xyyy_zzzzzzzz_0[j] = pb_x * tg_yyy_zzzzzzzz_0[j] + wp_x[j] * tg_yyy_zzzzzzzz_1[j];

                    tg_xyyz_xxxxxxxx_0[j] = pb_x * tg_yyz_xxxxxxxx_0[j] + wp_x[j] * tg_yyz_xxxxxxxx_1[j] + 4.0 * fl1_fxn * tg_yyz_xxxxxxx_1[j];

                    tg_xyyz_xxxxxxxy_0[j] = pb_x * tg_yyz_xxxxxxxy_0[j] + wp_x[j] * tg_yyz_xxxxxxxy_1[j] + 3.5 * fl1_fxn * tg_yyz_xxxxxxy_1[j];

                    tg_xyyz_xxxxxxxz_0[j] = pb_x * tg_yyz_xxxxxxxz_0[j] + wp_x[j] * tg_yyz_xxxxxxxz_1[j] + 3.5 * fl1_fxn * tg_yyz_xxxxxxz_1[j];

                    tg_xyyz_xxxxxxyy_0[j] = pb_x * tg_yyz_xxxxxxyy_0[j] + wp_x[j] * tg_yyz_xxxxxxyy_1[j] + 3.0 * fl1_fxn * tg_yyz_xxxxxyy_1[j];

                    tg_xyyz_xxxxxxyz_0[j] = pb_x * tg_yyz_xxxxxxyz_0[j] + wp_x[j] * tg_yyz_xxxxxxyz_1[j] + 3.0 * fl1_fxn * tg_yyz_xxxxxyz_1[j];

                    tg_xyyz_xxxxxxzz_0[j] = pb_x * tg_yyz_xxxxxxzz_0[j] + wp_x[j] * tg_yyz_xxxxxxzz_1[j] + 3.0 * fl1_fxn * tg_yyz_xxxxxzz_1[j];

                    tg_xyyz_xxxxxyyy_0[j] = pb_x * tg_yyz_xxxxxyyy_0[j] + wp_x[j] * tg_yyz_xxxxxyyy_1[j] + 2.5 * fl1_fxn * tg_yyz_xxxxyyy_1[j];

                    tg_xyyz_xxxxxyyz_0[j] = pb_x * tg_yyz_xxxxxyyz_0[j] + wp_x[j] * tg_yyz_xxxxxyyz_1[j] + 2.5 * fl1_fxn * tg_yyz_xxxxyyz_1[j];

                    tg_xyyz_xxxxxyzz_0[j] = pb_x * tg_yyz_xxxxxyzz_0[j] + wp_x[j] * tg_yyz_xxxxxyzz_1[j] + 2.5 * fl1_fxn * tg_yyz_xxxxyzz_1[j];

                    tg_xyyz_xxxxxzzz_0[j] = pb_x * tg_yyz_xxxxxzzz_0[j] + wp_x[j] * tg_yyz_xxxxxzzz_1[j] + 2.5 * fl1_fxn * tg_yyz_xxxxzzz_1[j];

                    tg_xyyz_xxxxyyyy_0[j] = pb_x * tg_yyz_xxxxyyyy_0[j] + wp_x[j] * tg_yyz_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_yyz_xxxyyyy_1[j];

                    tg_xyyz_xxxxyyyz_0[j] = pb_x * tg_yyz_xxxxyyyz_0[j] + wp_x[j] * tg_yyz_xxxxyyyz_1[j] + 2.0 * fl1_fxn * tg_yyz_xxxyyyz_1[j];

                    tg_xyyz_xxxxyyzz_0[j] = pb_x * tg_yyz_xxxxyyzz_0[j] + wp_x[j] * tg_yyz_xxxxyyzz_1[j] + 2.0 * fl1_fxn * tg_yyz_xxxyyzz_1[j];

                    tg_xyyz_xxxxyzzz_0[j] = pb_x * tg_yyz_xxxxyzzz_0[j] + wp_x[j] * tg_yyz_xxxxyzzz_1[j] + 2.0 * fl1_fxn * tg_yyz_xxxyzzz_1[j];

                    tg_xyyz_xxxxzzzz_0[j] = pb_x * tg_yyz_xxxxzzzz_0[j] + wp_x[j] * tg_yyz_xxxxzzzz_1[j] + 2.0 * fl1_fxn * tg_yyz_xxxzzzz_1[j];

                    tg_xyyz_xxxyyyyy_0[j] = pb_x * tg_yyz_xxxyyyyy_0[j] + wp_x[j] * tg_yyz_xxxyyyyy_1[j] + 1.5 * fl1_fxn * tg_yyz_xxyyyyy_1[j];

                    tg_xyyz_xxxyyyyz_0[j] = pb_x * tg_yyz_xxxyyyyz_0[j] + wp_x[j] * tg_yyz_xxxyyyyz_1[j] + 1.5 * fl1_fxn * tg_yyz_xxyyyyz_1[j];

                    tg_xyyz_xxxyyyzz_0[j] = pb_x * tg_yyz_xxxyyyzz_0[j] + wp_x[j] * tg_yyz_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_yyz_xxyyyzz_1[j];

                    tg_xyyz_xxxyyzzz_0[j] = pb_x * tg_yyz_xxxyyzzz_0[j] + wp_x[j] * tg_yyz_xxxyyzzz_1[j] + 1.5 * fl1_fxn * tg_yyz_xxyyzzz_1[j];

                    tg_xyyz_xxxyzzzz_0[j] = pb_x * tg_yyz_xxxyzzzz_0[j] + wp_x[j] * tg_yyz_xxxyzzzz_1[j] + 1.5 * fl1_fxn * tg_yyz_xxyzzzz_1[j];

                    tg_xyyz_xxxzzzzz_0[j] = pb_x * tg_yyz_xxxzzzzz_0[j] + wp_x[j] * tg_yyz_xxxzzzzz_1[j] + 1.5 * fl1_fxn * tg_yyz_xxzzzzz_1[j];

                    tg_xyyz_xxyyyyyy_0[j] = pb_x * tg_yyz_xxyyyyyy_0[j] + wp_x[j] * tg_yyz_xxyyyyyy_1[j] + fl1_fxn * tg_yyz_xyyyyyy_1[j];

                    tg_xyyz_xxyyyyyz_0[j] = pb_x * tg_yyz_xxyyyyyz_0[j] + wp_x[j] * tg_yyz_xxyyyyyz_1[j] + fl1_fxn * tg_yyz_xyyyyyz_1[j];

                    tg_xyyz_xxyyyyzz_0[j] = pb_x * tg_yyz_xxyyyyzz_0[j] + wp_x[j] * tg_yyz_xxyyyyzz_1[j] + fl1_fxn * tg_yyz_xyyyyzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSL_339_387(      CMemBlock2D<double>& primBuffer,
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

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {4, -1, -1, -1},
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_3_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {7, -1, -1, -1}, 
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

                auto tg_yyz_xxyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 339); 

                auto tg_yyz_xxyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 340); 

                auto tg_yyz_xxyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 341); 

                auto tg_yyz_xxzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 342); 

                auto tg_yyz_xyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 343); 

                auto tg_yyz_xyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 344); 

                auto tg_yyz_xyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 345); 

                auto tg_yyz_xyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 346); 

                auto tg_yyz_xyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 347); 

                auto tg_yyz_xyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 348); 

                auto tg_yyz_xyzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 349); 

                auto tg_yyz_xzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 350); 

                auto tg_yyz_yyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 351); 

                auto tg_yyz_yyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 352); 

                auto tg_yyz_yyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 353); 

                auto tg_yyz_yyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 354); 

                auto tg_yyz_yyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 355); 

                auto tg_yyz_yyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 356); 

                auto tg_yyz_yyzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 357); 

                auto tg_yyz_yzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 358); 

                auto tg_yyz_zzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 359); 

                auto tg_yzz_xxxxxxxx_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 360); 

                auto tg_yzz_xxxxxxxy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 361); 

                auto tg_yzz_xxxxxxxz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 362); 

                auto tg_yzz_xxxxxxyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 363); 

                auto tg_yzz_xxxxxxyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 364); 

                auto tg_yzz_xxxxxxzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 365); 

                auto tg_yzz_xxxxxyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 366); 

                auto tg_yzz_xxxxxyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 367); 

                auto tg_yzz_xxxxxyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 368); 

                auto tg_yzz_xxxxxzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 369); 

                auto tg_yzz_xxxxyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 370); 

                auto tg_yzz_xxxxyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 371); 

                auto tg_yzz_xxxxyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 372); 

                auto tg_yzz_xxxxyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 373); 

                auto tg_yzz_xxxxzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 374); 

                auto tg_yzz_xxxyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 375); 

                auto tg_yzz_xxxyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 376); 

                auto tg_yzz_xxxyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 377); 

                auto tg_yzz_xxxyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 378); 

                auto tg_yzz_xxxyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 379); 

                auto tg_yzz_xxxzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 380); 

                auto tg_yzz_xxyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 381); 

                auto tg_yzz_xxyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 382); 

                auto tg_yzz_xxyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 383); 

                auto tg_yzz_xxyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 384); 

                auto tg_yzz_xxyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 385); 

                auto tg_yzz_xxyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 386); 

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

                // set up pointers to integrals

                auto tg_xyyz_xxyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 339); 

                auto tg_xyyz_xxyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 340); 

                auto tg_xyyz_xxyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 341); 

                auto tg_xyyz_xxzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 342); 

                auto tg_xyyz_xyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 343); 

                auto tg_xyyz_xyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 344); 

                auto tg_xyyz_xyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 345); 

                auto tg_xyyz_xyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 346); 

                auto tg_xyyz_xyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 347); 

                auto tg_xyyz_xyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 348); 

                auto tg_xyyz_xyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 349); 

                auto tg_xyyz_xzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 350); 

                auto tg_xyyz_yyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 351); 

                auto tg_xyyz_yyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 352); 

                auto tg_xyyz_yyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 353); 

                auto tg_xyyz_yyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 354); 

                auto tg_xyyz_yyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 355); 

                auto tg_xyyz_yyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 356); 

                auto tg_xyyz_yyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 357); 

                auto tg_xyyz_yzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 358); 

                auto tg_xyyz_zzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 359); 

                auto tg_xyzz_xxxxxxxx_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 360); 

                auto tg_xyzz_xxxxxxxy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 361); 

                auto tg_xyzz_xxxxxxxz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 362); 

                auto tg_xyzz_xxxxxxyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 363); 

                auto tg_xyzz_xxxxxxyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 364); 

                auto tg_xyzz_xxxxxxzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 365); 

                auto tg_xyzz_xxxxxyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 366); 

                auto tg_xyzz_xxxxxyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 367); 

                auto tg_xyzz_xxxxxyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 368); 

                auto tg_xyzz_xxxxxzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 369); 

                auto tg_xyzz_xxxxyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 370); 

                auto tg_xyzz_xxxxyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 371); 

                auto tg_xyzz_xxxxyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 372); 

                auto tg_xyzz_xxxxyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 373); 

                auto tg_xyzz_xxxxzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 374); 

                auto tg_xyzz_xxxyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 375); 

                auto tg_xyzz_xxxyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 376); 

                auto tg_xyzz_xxxyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 377); 

                auto tg_xyzz_xxxyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 378); 

                auto tg_xyzz_xxxyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 379); 

                auto tg_xyzz_xxxzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 380); 

                auto tg_xyzz_xxyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 381); 

                auto tg_xyzz_xxyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 382); 

                auto tg_xyzz_xxyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 383); 

                auto tg_xyzz_xxyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 384); 

                auto tg_xyzz_xxyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 385); 

                auto tg_xyzz_xxyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 386); 

                // Batch of Integrals (339,387)

                #pragma omp simd aligned(fxn, tg_xyyz_xxyyyzzz_0, tg_xyyz_xxyyzzzz_0, tg_xyyz_xxyzzzzz_0, \
                                         tg_xyyz_xxzzzzzz_0, tg_xyyz_xyyyyyyy_0, tg_xyyz_xyyyyyyz_0, tg_xyyz_xyyyyyzz_0, \
                                         tg_xyyz_xyyyyzzz_0, tg_xyyz_xyyyzzzz_0, tg_xyyz_xyyzzzzz_0, tg_xyyz_xyzzzzzz_0, \
                                         tg_xyyz_xzzzzzzz_0, tg_xyyz_yyyyyyyy_0, tg_xyyz_yyyyyyyz_0, tg_xyyz_yyyyyyzz_0, \
                                         tg_xyyz_yyyyyzzz_0, tg_xyyz_yyyyzzzz_0, tg_xyyz_yyyzzzzz_0, tg_xyyz_yyzzzzzz_0, \
                                         tg_xyyz_yzzzzzzz_0, tg_xyyz_zzzzzzzz_0, tg_xyzz_xxxxxxxx_0, tg_xyzz_xxxxxxxy_0, \
                                         tg_xyzz_xxxxxxxz_0, tg_xyzz_xxxxxxyy_0, tg_xyzz_xxxxxxyz_0, tg_xyzz_xxxxxxzz_0, \
                                         tg_xyzz_xxxxxyyy_0, tg_xyzz_xxxxxyyz_0, tg_xyzz_xxxxxyzz_0, tg_xyzz_xxxxxzzz_0, \
                                         tg_xyzz_xxxxyyyy_0, tg_xyzz_xxxxyyyz_0, tg_xyzz_xxxxyyzz_0, tg_xyzz_xxxxyzzz_0, \
                                         tg_xyzz_xxxxzzzz_0, tg_xyzz_xxxyyyyy_0, tg_xyzz_xxxyyyyz_0, tg_xyzz_xxxyyyzz_0, \
                                         tg_xyzz_xxxyyzzz_0, tg_xyzz_xxxyzzzz_0, tg_xyzz_xxxzzzzz_0, tg_xyzz_xxyyyyyy_0, \
                                         tg_xyzz_xxyyyyyz_0, tg_xyzz_xxyyyyzz_0, tg_xyzz_xxyyyzzz_0, tg_xyzz_xxyyzzzz_0, \
                                         tg_xyzz_xxyzzzzz_0, tg_yyz_xxyyyzzz_0, tg_yyz_xxyyyzzz_1, tg_yyz_xxyyzzzz_0, \
                                         tg_yyz_xxyyzzzz_1, tg_yyz_xxyzzzzz_0, tg_yyz_xxyzzzzz_1, tg_yyz_xxzzzzzz_0, \
                                         tg_yyz_xxzzzzzz_1, tg_yyz_xyyyyyyy_0, tg_yyz_xyyyyyyy_1, tg_yyz_xyyyyyyz_0, \
                                         tg_yyz_xyyyyyyz_1, tg_yyz_xyyyyyzz_0, tg_yyz_xyyyyyzz_1, tg_yyz_xyyyyzzz_0, \
                                         tg_yyz_xyyyyzzz_1, tg_yyz_xyyyzzz_1, tg_yyz_xyyyzzzz_0, tg_yyz_xyyyzzzz_1, \
                                         tg_yyz_xyyzzzz_1, tg_yyz_xyyzzzzz_0, tg_yyz_xyyzzzzz_1, tg_yyz_xyzzzzz_1, \
                                         tg_yyz_xyzzzzzz_0, tg_yyz_xyzzzzzz_1, tg_yyz_xzzzzzz_1, tg_yyz_xzzzzzzz_0, \
                                         tg_yyz_xzzzzzzz_1, tg_yyz_yyyyyyy_1, tg_yyz_yyyyyyyy_0, tg_yyz_yyyyyyyy_1, \
                                         tg_yyz_yyyyyyyz_0, tg_yyz_yyyyyyyz_1, tg_yyz_yyyyyyz_1, tg_yyz_yyyyyyzz_0, \
                                         tg_yyz_yyyyyyzz_1, tg_yyz_yyyyyzz_1, tg_yyz_yyyyyzzz_0, tg_yyz_yyyyyzzz_1, \
                                         tg_yyz_yyyyzzz_1, tg_yyz_yyyyzzzz_0, tg_yyz_yyyyzzzz_1, tg_yyz_yyyzzzz_1, \
                                         tg_yyz_yyyzzzzz_0, tg_yyz_yyyzzzzz_1, tg_yyz_yyzzzzz_1, tg_yyz_yyzzzzzz_0, \
                                         tg_yyz_yyzzzzzz_1, tg_yyz_yzzzzzz_1, tg_yyz_yzzzzzzz_0, tg_yyz_yzzzzzzz_1, \
                                         tg_yyz_zzzzzzz_1, tg_yyz_zzzzzzzz_0, tg_yyz_zzzzzzzz_1, tg_yzz_xxxxxxx_1, \
                                         tg_yzz_xxxxxxxx_0, tg_yzz_xxxxxxxx_1, tg_yzz_xxxxxxxy_0, tg_yzz_xxxxxxxy_1, \
                                         tg_yzz_xxxxxxxz_0, tg_yzz_xxxxxxxz_1, tg_yzz_xxxxxxy_1, tg_yzz_xxxxxxyy_0, \
                                         tg_yzz_xxxxxxyy_1, tg_yzz_xxxxxxyz_0, tg_yzz_xxxxxxyz_1, tg_yzz_xxxxxxz_1, \
                                         tg_yzz_xxxxxxzz_0, tg_yzz_xxxxxxzz_1, tg_yzz_xxxxxyy_1, tg_yzz_xxxxxyyy_0, \
                                         tg_yzz_xxxxxyyy_1, tg_yzz_xxxxxyyz_0, tg_yzz_xxxxxyyz_1, tg_yzz_xxxxxyz_1, \
                                         tg_yzz_xxxxxyzz_0, tg_yzz_xxxxxyzz_1, tg_yzz_xxxxxzz_1, tg_yzz_xxxxxzzz_0, \
                                         tg_yzz_xxxxxzzz_1, tg_yzz_xxxxyyy_1, tg_yzz_xxxxyyyy_0, tg_yzz_xxxxyyyy_1, \
                                         tg_yzz_xxxxyyyz_0, tg_yzz_xxxxyyyz_1, tg_yzz_xxxxyyz_1, tg_yzz_xxxxyyzz_0, \
                                         tg_yzz_xxxxyyzz_1, tg_yzz_xxxxyzz_1, tg_yzz_xxxxyzzz_0, tg_yzz_xxxxyzzz_1, \
                                         tg_yzz_xxxxzzz_1, tg_yzz_xxxxzzzz_0, tg_yzz_xxxxzzzz_1, tg_yzz_xxxyyyy_1, \
                                         tg_yzz_xxxyyyyy_0, tg_yzz_xxxyyyyy_1, tg_yzz_xxxyyyyz_0, tg_yzz_xxxyyyyz_1, \
                                         tg_yzz_xxxyyyz_1, tg_yzz_xxxyyyzz_0, tg_yzz_xxxyyyzz_1, tg_yzz_xxxyyzz_1, \
                                         tg_yzz_xxxyyzzz_0, tg_yzz_xxxyyzzz_1, tg_yzz_xxxyzzz_1, tg_yzz_xxxyzzzz_0, \
                                         tg_yzz_xxxyzzzz_1, tg_yzz_xxxzzzz_1, tg_yzz_xxxzzzzz_0, tg_yzz_xxxzzzzz_1, \
                                         tg_yzz_xxyyyyy_1, tg_yzz_xxyyyyyy_0, tg_yzz_xxyyyyyy_1, tg_yzz_xxyyyyyz_0, \
                                         tg_yzz_xxyyyyyz_1, tg_yzz_xxyyyyz_1, tg_yzz_xxyyyyzz_0, tg_yzz_xxyyyyzz_1, \
                                         tg_yzz_xxyyyzz_1, tg_yzz_xxyyyzzz_0, tg_yzz_xxyyyzzz_1, tg_yzz_xxyyzzz_1, \
                                         tg_yzz_xxyyzzzz_0, tg_yzz_xxyyzzzz_1, tg_yzz_xxyzzzz_1, tg_yzz_xxyzzzzz_0, \
                                         tg_yzz_xxyzzzzz_1, tg_yzz_xxzzzzz_1, tg_yzz_xyyyyyy_1, tg_yzz_xyyyyyz_1, \
                                         tg_yzz_xyyyyzz_1, tg_yzz_xyyyzzz_1, tg_yzz_xyyzzzz_1, tg_yzz_xyzzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fxn = fxn[j];

                    tg_xyyz_xxyyyzzz_0[j] = pb_x * tg_yyz_xxyyyzzz_0[j] + wp_x[j] * tg_yyz_xxyyyzzz_1[j] + fl1_fxn * tg_yyz_xyyyzzz_1[j];

                    tg_xyyz_xxyyzzzz_0[j] = pb_x * tg_yyz_xxyyzzzz_0[j] + wp_x[j] * tg_yyz_xxyyzzzz_1[j] + fl1_fxn * tg_yyz_xyyzzzz_1[j];

                    tg_xyyz_xxyzzzzz_0[j] = pb_x * tg_yyz_xxyzzzzz_0[j] + wp_x[j] * tg_yyz_xxyzzzzz_1[j] + fl1_fxn * tg_yyz_xyzzzzz_1[j];

                    tg_xyyz_xxzzzzzz_0[j] = pb_x * tg_yyz_xxzzzzzz_0[j] + wp_x[j] * tg_yyz_xxzzzzzz_1[j] + fl1_fxn * tg_yyz_xzzzzzz_1[j];

                    tg_xyyz_xyyyyyyy_0[j] = pb_x * tg_yyz_xyyyyyyy_0[j] + wp_x[j] * tg_yyz_xyyyyyyy_1[j] + 0.5 * fl1_fxn * tg_yyz_yyyyyyy_1[j];

                    tg_xyyz_xyyyyyyz_0[j] = pb_x * tg_yyz_xyyyyyyz_0[j] + wp_x[j] * tg_yyz_xyyyyyyz_1[j] + 0.5 * fl1_fxn * tg_yyz_yyyyyyz_1[j];

                    tg_xyyz_xyyyyyzz_0[j] = pb_x * tg_yyz_xyyyyyzz_0[j] + wp_x[j] * tg_yyz_xyyyyyzz_1[j] + 0.5 * fl1_fxn * tg_yyz_yyyyyzz_1[j];

                    tg_xyyz_xyyyyzzz_0[j] = pb_x * tg_yyz_xyyyyzzz_0[j] + wp_x[j] * tg_yyz_xyyyyzzz_1[j] + 0.5 * fl1_fxn * tg_yyz_yyyyzzz_1[j];

                    tg_xyyz_xyyyzzzz_0[j] = pb_x * tg_yyz_xyyyzzzz_0[j] + wp_x[j] * tg_yyz_xyyyzzzz_1[j] + 0.5 * fl1_fxn * tg_yyz_yyyzzzz_1[j];

                    tg_xyyz_xyyzzzzz_0[j] = pb_x * tg_yyz_xyyzzzzz_0[j] + wp_x[j] * tg_yyz_xyyzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyz_yyzzzzz_1[j];

                    tg_xyyz_xyzzzzzz_0[j] = pb_x * tg_yyz_xyzzzzzz_0[j] + wp_x[j] * tg_yyz_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyz_yzzzzzz_1[j];

                    tg_xyyz_xzzzzzzz_0[j] = pb_x * tg_yyz_xzzzzzzz_0[j] + wp_x[j] * tg_yyz_xzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyz_zzzzzzz_1[j];

                    tg_xyyz_yyyyyyyy_0[j] = pb_x * tg_yyz_yyyyyyyy_0[j] + wp_x[j] * tg_yyz_yyyyyyyy_1[j];

                    tg_xyyz_yyyyyyyz_0[j] = pb_x * tg_yyz_yyyyyyyz_0[j] + wp_x[j] * tg_yyz_yyyyyyyz_1[j];

                    tg_xyyz_yyyyyyzz_0[j] = pb_x * tg_yyz_yyyyyyzz_0[j] + wp_x[j] * tg_yyz_yyyyyyzz_1[j];

                    tg_xyyz_yyyyyzzz_0[j] = pb_x * tg_yyz_yyyyyzzz_0[j] + wp_x[j] * tg_yyz_yyyyyzzz_1[j];

                    tg_xyyz_yyyyzzzz_0[j] = pb_x * tg_yyz_yyyyzzzz_0[j] + wp_x[j] * tg_yyz_yyyyzzzz_1[j];

                    tg_xyyz_yyyzzzzz_0[j] = pb_x * tg_yyz_yyyzzzzz_0[j] + wp_x[j] * tg_yyz_yyyzzzzz_1[j];

                    tg_xyyz_yyzzzzzz_0[j] = pb_x * tg_yyz_yyzzzzzz_0[j] + wp_x[j] * tg_yyz_yyzzzzzz_1[j];

                    tg_xyyz_yzzzzzzz_0[j] = pb_x * tg_yyz_yzzzzzzz_0[j] + wp_x[j] * tg_yyz_yzzzzzzz_1[j];

                    tg_xyyz_zzzzzzzz_0[j] = pb_x * tg_yyz_zzzzzzzz_0[j] + wp_x[j] * tg_yyz_zzzzzzzz_1[j];

                    tg_xyzz_xxxxxxxx_0[j] = pb_x * tg_yzz_xxxxxxxx_0[j] + wp_x[j] * tg_yzz_xxxxxxxx_1[j] + 4.0 * fl1_fxn * tg_yzz_xxxxxxx_1[j];

                    tg_xyzz_xxxxxxxy_0[j] = pb_x * tg_yzz_xxxxxxxy_0[j] + wp_x[j] * tg_yzz_xxxxxxxy_1[j] + 3.5 * fl1_fxn * tg_yzz_xxxxxxy_1[j];

                    tg_xyzz_xxxxxxxz_0[j] = pb_x * tg_yzz_xxxxxxxz_0[j] + wp_x[j] * tg_yzz_xxxxxxxz_1[j] + 3.5 * fl1_fxn * tg_yzz_xxxxxxz_1[j];

                    tg_xyzz_xxxxxxyy_0[j] = pb_x * tg_yzz_xxxxxxyy_0[j] + wp_x[j] * tg_yzz_xxxxxxyy_1[j] + 3.0 * fl1_fxn * tg_yzz_xxxxxyy_1[j];

                    tg_xyzz_xxxxxxyz_0[j] = pb_x * tg_yzz_xxxxxxyz_0[j] + wp_x[j] * tg_yzz_xxxxxxyz_1[j] + 3.0 * fl1_fxn * tg_yzz_xxxxxyz_1[j];

                    tg_xyzz_xxxxxxzz_0[j] = pb_x * tg_yzz_xxxxxxzz_0[j] + wp_x[j] * tg_yzz_xxxxxxzz_1[j] + 3.0 * fl1_fxn * tg_yzz_xxxxxzz_1[j];

                    tg_xyzz_xxxxxyyy_0[j] = pb_x * tg_yzz_xxxxxyyy_0[j] + wp_x[j] * tg_yzz_xxxxxyyy_1[j] + 2.5 * fl1_fxn * tg_yzz_xxxxyyy_1[j];

                    tg_xyzz_xxxxxyyz_0[j] = pb_x * tg_yzz_xxxxxyyz_0[j] + wp_x[j] * tg_yzz_xxxxxyyz_1[j] + 2.5 * fl1_fxn * tg_yzz_xxxxyyz_1[j];

                    tg_xyzz_xxxxxyzz_0[j] = pb_x * tg_yzz_xxxxxyzz_0[j] + wp_x[j] * tg_yzz_xxxxxyzz_1[j] + 2.5 * fl1_fxn * tg_yzz_xxxxyzz_1[j];

                    tg_xyzz_xxxxxzzz_0[j] = pb_x * tg_yzz_xxxxxzzz_0[j] + wp_x[j] * tg_yzz_xxxxxzzz_1[j] + 2.5 * fl1_fxn * tg_yzz_xxxxzzz_1[j];

                    tg_xyzz_xxxxyyyy_0[j] = pb_x * tg_yzz_xxxxyyyy_0[j] + wp_x[j] * tg_yzz_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_yzz_xxxyyyy_1[j];

                    tg_xyzz_xxxxyyyz_0[j] = pb_x * tg_yzz_xxxxyyyz_0[j] + wp_x[j] * tg_yzz_xxxxyyyz_1[j] + 2.0 * fl1_fxn * tg_yzz_xxxyyyz_1[j];

                    tg_xyzz_xxxxyyzz_0[j] = pb_x * tg_yzz_xxxxyyzz_0[j] + wp_x[j] * tg_yzz_xxxxyyzz_1[j] + 2.0 * fl1_fxn * tg_yzz_xxxyyzz_1[j];

                    tg_xyzz_xxxxyzzz_0[j] = pb_x * tg_yzz_xxxxyzzz_0[j] + wp_x[j] * tg_yzz_xxxxyzzz_1[j] + 2.0 * fl1_fxn * tg_yzz_xxxyzzz_1[j];

                    tg_xyzz_xxxxzzzz_0[j] = pb_x * tg_yzz_xxxxzzzz_0[j] + wp_x[j] * tg_yzz_xxxxzzzz_1[j] + 2.0 * fl1_fxn * tg_yzz_xxxzzzz_1[j];

                    tg_xyzz_xxxyyyyy_0[j] = pb_x * tg_yzz_xxxyyyyy_0[j] + wp_x[j] * tg_yzz_xxxyyyyy_1[j] + 1.5 * fl1_fxn * tg_yzz_xxyyyyy_1[j];

                    tg_xyzz_xxxyyyyz_0[j] = pb_x * tg_yzz_xxxyyyyz_0[j] + wp_x[j] * tg_yzz_xxxyyyyz_1[j] + 1.5 * fl1_fxn * tg_yzz_xxyyyyz_1[j];

                    tg_xyzz_xxxyyyzz_0[j] = pb_x * tg_yzz_xxxyyyzz_0[j] + wp_x[j] * tg_yzz_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_yzz_xxyyyzz_1[j];

                    tg_xyzz_xxxyyzzz_0[j] = pb_x * tg_yzz_xxxyyzzz_0[j] + wp_x[j] * tg_yzz_xxxyyzzz_1[j] + 1.5 * fl1_fxn * tg_yzz_xxyyzzz_1[j];

                    tg_xyzz_xxxyzzzz_0[j] = pb_x * tg_yzz_xxxyzzzz_0[j] + wp_x[j] * tg_yzz_xxxyzzzz_1[j] + 1.5 * fl1_fxn * tg_yzz_xxyzzzz_1[j];

                    tg_xyzz_xxxzzzzz_0[j] = pb_x * tg_yzz_xxxzzzzz_0[j] + wp_x[j] * tg_yzz_xxxzzzzz_1[j] + 1.5 * fl1_fxn * tg_yzz_xxzzzzz_1[j];

                    tg_xyzz_xxyyyyyy_0[j] = pb_x * tg_yzz_xxyyyyyy_0[j] + wp_x[j] * tg_yzz_xxyyyyyy_1[j] + fl1_fxn * tg_yzz_xyyyyyy_1[j];

                    tg_xyzz_xxyyyyyz_0[j] = pb_x * tg_yzz_xxyyyyyz_0[j] + wp_x[j] * tg_yzz_xxyyyyyz_1[j] + fl1_fxn * tg_yzz_xyyyyyz_1[j];

                    tg_xyzz_xxyyyyzz_0[j] = pb_x * tg_yzz_xxyyyyzz_0[j] + wp_x[j] * tg_yzz_xxyyyyzz_1[j] + fl1_fxn * tg_yzz_xyyyyzz_1[j];

                    tg_xyzz_xxyyyzzz_0[j] = pb_x * tg_yzz_xxyyyzzz_0[j] + wp_x[j] * tg_yzz_xxyyyzzz_1[j] + fl1_fxn * tg_yzz_xyyyzzz_1[j];

                    tg_xyzz_xxyyzzzz_0[j] = pb_x * tg_yzz_xxyyzzzz_0[j] + wp_x[j] * tg_yzz_xxyyzzzz_1[j] + fl1_fxn * tg_yzz_xyyzzzz_1[j];

                    tg_xyzz_xxyzzzzz_0[j] = pb_x * tg_yzz_xxyzzzzz_0[j] + wp_x[j] * tg_yzz_xxyzzzzz_1[j] + fl1_fxn * tg_yzz_xyzzzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSL_387_435(      CMemBlock2D<double>& primBuffer,
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

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {4, -1, -1, -1},
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_3_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {7, -1, -1, -1}, 
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

                auto tg_yzz_xxzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 387); 

                auto tg_yzz_xyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 388); 

                auto tg_yzz_xyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 389); 

                auto tg_yzz_xyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 390); 

                auto tg_yzz_xyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 391); 

                auto tg_yzz_xyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 392); 

                auto tg_yzz_xyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 393); 

                auto tg_yzz_xyzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 394); 

                auto tg_yzz_xzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 395); 

                auto tg_yzz_yyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 396); 

                auto tg_yzz_yyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 397); 

                auto tg_yzz_yyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 398); 

                auto tg_yzz_yyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 399); 

                auto tg_yzz_yyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 400); 

                auto tg_yzz_yyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 401); 

                auto tg_yzz_yyzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 402); 

                auto tg_yzz_yzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 403); 

                auto tg_yzz_zzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 404); 

                auto tg_zzz_xxxxxxxx_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 405); 

                auto tg_zzz_xxxxxxxy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 406); 

                auto tg_zzz_xxxxxxxz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 407); 

                auto tg_zzz_xxxxxxyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 408); 

                auto tg_zzz_xxxxxxyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 409); 

                auto tg_zzz_xxxxxxzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 410); 

                auto tg_zzz_xxxxxyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 411); 

                auto tg_zzz_xxxxxyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 412); 

                auto tg_zzz_xxxxxyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 413); 

                auto tg_zzz_xxxxxzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 414); 

                auto tg_zzz_xxxxyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 415); 

                auto tg_zzz_xxxxyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 416); 

                auto tg_zzz_xxxxyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 417); 

                auto tg_zzz_xxxxyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 418); 

                auto tg_zzz_xxxxzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 419); 

                auto tg_zzz_xxxyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 420); 

                auto tg_zzz_xxxyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 421); 

                auto tg_zzz_xxxyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 422); 

                auto tg_zzz_xxxyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 423); 

                auto tg_zzz_xxxyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 424); 

                auto tg_zzz_xxxzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 425); 

                auto tg_zzz_xxyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 426); 

                auto tg_zzz_xxyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 427); 

                auto tg_zzz_xxyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 428); 

                auto tg_zzz_xxyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 429); 

                auto tg_zzz_xxyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 430); 

                auto tg_zzz_xxyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 431); 

                auto tg_zzz_xxzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 432); 

                auto tg_zzz_xyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 433); 

                auto tg_zzz_xyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 434); 

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

                // set up pointers to integrals

                auto tg_xyzz_xxzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 387); 

                auto tg_xyzz_xyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 388); 

                auto tg_xyzz_xyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 389); 

                auto tg_xyzz_xyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 390); 

                auto tg_xyzz_xyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 391); 

                auto tg_xyzz_xyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 392); 

                auto tg_xyzz_xyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 393); 

                auto tg_xyzz_xyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 394); 

                auto tg_xyzz_xzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 395); 

                auto tg_xyzz_yyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 396); 

                auto tg_xyzz_yyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 397); 

                auto tg_xyzz_yyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 398); 

                auto tg_xyzz_yyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 399); 

                auto tg_xyzz_yyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 400); 

                auto tg_xyzz_yyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 401); 

                auto tg_xyzz_yyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 402); 

                auto tg_xyzz_yzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 403); 

                auto tg_xyzz_zzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 404); 

                auto tg_xzzz_xxxxxxxx_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 405); 

                auto tg_xzzz_xxxxxxxy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 406); 

                auto tg_xzzz_xxxxxxxz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 407); 

                auto tg_xzzz_xxxxxxyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 408); 

                auto tg_xzzz_xxxxxxyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 409); 

                auto tg_xzzz_xxxxxxzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 410); 

                auto tg_xzzz_xxxxxyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 411); 

                auto tg_xzzz_xxxxxyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 412); 

                auto tg_xzzz_xxxxxyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 413); 

                auto tg_xzzz_xxxxxzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 414); 

                auto tg_xzzz_xxxxyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 415); 

                auto tg_xzzz_xxxxyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 416); 

                auto tg_xzzz_xxxxyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 417); 

                auto tg_xzzz_xxxxyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 418); 

                auto tg_xzzz_xxxxzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 419); 

                auto tg_xzzz_xxxyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 420); 

                auto tg_xzzz_xxxyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 421); 

                auto tg_xzzz_xxxyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 422); 

                auto tg_xzzz_xxxyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 423); 

                auto tg_xzzz_xxxyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 424); 

                auto tg_xzzz_xxxzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 425); 

                auto tg_xzzz_xxyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 426); 

                auto tg_xzzz_xxyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 427); 

                auto tg_xzzz_xxyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 428); 

                auto tg_xzzz_xxyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 429); 

                auto tg_xzzz_xxyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 430); 

                auto tg_xzzz_xxyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 431); 

                auto tg_xzzz_xxzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 432); 

                auto tg_xzzz_xyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 433); 

                auto tg_xzzz_xyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 434); 

                // Batch of Integrals (387,435)

                #pragma omp simd aligned(fxn, tg_xyzz_xxzzzzzz_0, tg_xyzz_xyyyyyyy_0, tg_xyzz_xyyyyyyz_0, \
                                         tg_xyzz_xyyyyyzz_0, tg_xyzz_xyyyyzzz_0, tg_xyzz_xyyyzzzz_0, tg_xyzz_xyyzzzzz_0, \
                                         tg_xyzz_xyzzzzzz_0, tg_xyzz_xzzzzzzz_0, tg_xyzz_yyyyyyyy_0, tg_xyzz_yyyyyyyz_0, \
                                         tg_xyzz_yyyyyyzz_0, tg_xyzz_yyyyyzzz_0, tg_xyzz_yyyyzzzz_0, tg_xyzz_yyyzzzzz_0, \
                                         tg_xyzz_yyzzzzzz_0, tg_xyzz_yzzzzzzz_0, tg_xyzz_zzzzzzzz_0, tg_xzzz_xxxxxxxx_0, \
                                         tg_xzzz_xxxxxxxy_0, tg_xzzz_xxxxxxxz_0, tg_xzzz_xxxxxxyy_0, tg_xzzz_xxxxxxyz_0, \
                                         tg_xzzz_xxxxxxzz_0, tg_xzzz_xxxxxyyy_0, tg_xzzz_xxxxxyyz_0, tg_xzzz_xxxxxyzz_0, \
                                         tg_xzzz_xxxxxzzz_0, tg_xzzz_xxxxyyyy_0, tg_xzzz_xxxxyyyz_0, tg_xzzz_xxxxyyzz_0, \
                                         tg_xzzz_xxxxyzzz_0, tg_xzzz_xxxxzzzz_0, tg_xzzz_xxxyyyyy_0, tg_xzzz_xxxyyyyz_0, \
                                         tg_xzzz_xxxyyyzz_0, tg_xzzz_xxxyyzzz_0, tg_xzzz_xxxyzzzz_0, tg_xzzz_xxxzzzzz_0, \
                                         tg_xzzz_xxyyyyyy_0, tg_xzzz_xxyyyyyz_0, tg_xzzz_xxyyyyzz_0, tg_xzzz_xxyyyzzz_0, \
                                         tg_xzzz_xxyyzzzz_0, tg_xzzz_xxyzzzzz_0, tg_xzzz_xxzzzzzz_0, tg_xzzz_xyyyyyyy_0, \
                                         tg_xzzz_xyyyyyyz_0, tg_yzz_xxzzzzzz_0, tg_yzz_xxzzzzzz_1, tg_yzz_xyyyyyyy_0, \
                                         tg_yzz_xyyyyyyy_1, tg_yzz_xyyyyyyz_0, tg_yzz_xyyyyyyz_1, tg_yzz_xyyyyyzz_0, \
                                         tg_yzz_xyyyyyzz_1, tg_yzz_xyyyyzzz_0, tg_yzz_xyyyyzzz_1, tg_yzz_xyyyzzzz_0, \
                                         tg_yzz_xyyyzzzz_1, tg_yzz_xyyzzzzz_0, tg_yzz_xyyzzzzz_1, tg_yzz_xyzzzzzz_0, \
                                         tg_yzz_xyzzzzzz_1, tg_yzz_xzzzzzz_1, tg_yzz_xzzzzzzz_0, tg_yzz_xzzzzzzz_1, \
                                         tg_yzz_yyyyyyy_1, tg_yzz_yyyyyyyy_0, tg_yzz_yyyyyyyy_1, tg_yzz_yyyyyyyz_0, \
                                         tg_yzz_yyyyyyyz_1, tg_yzz_yyyyyyz_1, tg_yzz_yyyyyyzz_0, tg_yzz_yyyyyyzz_1, \
                                         tg_yzz_yyyyyzz_1, tg_yzz_yyyyyzzz_0, tg_yzz_yyyyyzzz_1, tg_yzz_yyyyzzz_1, \
                                         tg_yzz_yyyyzzzz_0, tg_yzz_yyyyzzzz_1, tg_yzz_yyyzzzz_1, tg_yzz_yyyzzzzz_0, \
                                         tg_yzz_yyyzzzzz_1, tg_yzz_yyzzzzz_1, tg_yzz_yyzzzzzz_0, tg_yzz_yyzzzzzz_1, \
                                         tg_yzz_yzzzzzz_1, tg_yzz_yzzzzzzz_0, tg_yzz_yzzzzzzz_1, tg_yzz_zzzzzzz_1, \
                                         tg_yzz_zzzzzzzz_0, tg_yzz_zzzzzzzz_1, tg_zzz_xxxxxxx_1, tg_zzz_xxxxxxxx_0, \
                                         tg_zzz_xxxxxxxx_1, tg_zzz_xxxxxxxy_0, tg_zzz_xxxxxxxy_1, tg_zzz_xxxxxxxz_0, \
                                         tg_zzz_xxxxxxxz_1, tg_zzz_xxxxxxy_1, tg_zzz_xxxxxxyy_0, tg_zzz_xxxxxxyy_1, \
                                         tg_zzz_xxxxxxyz_0, tg_zzz_xxxxxxyz_1, tg_zzz_xxxxxxz_1, tg_zzz_xxxxxxzz_0, \
                                         tg_zzz_xxxxxxzz_1, tg_zzz_xxxxxyy_1, tg_zzz_xxxxxyyy_0, tg_zzz_xxxxxyyy_1, \
                                         tg_zzz_xxxxxyyz_0, tg_zzz_xxxxxyyz_1, tg_zzz_xxxxxyz_1, tg_zzz_xxxxxyzz_0, \
                                         tg_zzz_xxxxxyzz_1, tg_zzz_xxxxxzz_1, tg_zzz_xxxxxzzz_0, tg_zzz_xxxxxzzz_1, \
                                         tg_zzz_xxxxyyy_1, tg_zzz_xxxxyyyy_0, tg_zzz_xxxxyyyy_1, tg_zzz_xxxxyyyz_0, \
                                         tg_zzz_xxxxyyyz_1, tg_zzz_xxxxyyz_1, tg_zzz_xxxxyyzz_0, tg_zzz_xxxxyyzz_1, \
                                         tg_zzz_xxxxyzz_1, tg_zzz_xxxxyzzz_0, tg_zzz_xxxxyzzz_1, tg_zzz_xxxxzzz_1, \
                                         tg_zzz_xxxxzzzz_0, tg_zzz_xxxxzzzz_1, tg_zzz_xxxyyyy_1, tg_zzz_xxxyyyyy_0, \
                                         tg_zzz_xxxyyyyy_1, tg_zzz_xxxyyyyz_0, tg_zzz_xxxyyyyz_1, tg_zzz_xxxyyyz_1, \
                                         tg_zzz_xxxyyyzz_0, tg_zzz_xxxyyyzz_1, tg_zzz_xxxyyzz_1, tg_zzz_xxxyyzzz_0, \
                                         tg_zzz_xxxyyzzz_1, tg_zzz_xxxyzzz_1, tg_zzz_xxxyzzzz_0, tg_zzz_xxxyzzzz_1, \
                                         tg_zzz_xxxzzzz_1, tg_zzz_xxxzzzzz_0, tg_zzz_xxxzzzzz_1, tg_zzz_xxyyyyy_1, \
                                         tg_zzz_xxyyyyyy_0, tg_zzz_xxyyyyyy_1, tg_zzz_xxyyyyyz_0, tg_zzz_xxyyyyyz_1, \
                                         tg_zzz_xxyyyyz_1, tg_zzz_xxyyyyzz_0, tg_zzz_xxyyyyzz_1, tg_zzz_xxyyyzz_1, \
                                         tg_zzz_xxyyyzzz_0, tg_zzz_xxyyyzzz_1, tg_zzz_xxyyzzz_1, tg_zzz_xxyyzzzz_0, \
                                         tg_zzz_xxyyzzzz_1, tg_zzz_xxyzzzz_1, tg_zzz_xxyzzzzz_0, tg_zzz_xxyzzzzz_1, \
                                         tg_zzz_xxzzzzz_1, tg_zzz_xxzzzzzz_0, tg_zzz_xxzzzzzz_1, tg_zzz_xyyyyyy_1, \
                                         tg_zzz_xyyyyyyy_0, tg_zzz_xyyyyyyy_1, tg_zzz_xyyyyyyz_0, tg_zzz_xyyyyyyz_1, \
                                         tg_zzz_xyyyyyz_1, tg_zzz_xyyyyzz_1, tg_zzz_xyyyzzz_1, tg_zzz_xyyzzzz_1, \
                                         tg_zzz_xyzzzzz_1, tg_zzz_xzzzzzz_1, tg_zzz_yyyyyyy_1, tg_zzz_yyyyyyz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fxn = fxn[j];

                    tg_xyzz_xxzzzzzz_0[j] = pb_x * tg_yzz_xxzzzzzz_0[j] + wp_x[j] * tg_yzz_xxzzzzzz_1[j] + fl1_fxn * tg_yzz_xzzzzzz_1[j];

                    tg_xyzz_xyyyyyyy_0[j] = pb_x * tg_yzz_xyyyyyyy_0[j] + wp_x[j] * tg_yzz_xyyyyyyy_1[j] + 0.5 * fl1_fxn * tg_yzz_yyyyyyy_1[j];

                    tg_xyzz_xyyyyyyz_0[j] = pb_x * tg_yzz_xyyyyyyz_0[j] + wp_x[j] * tg_yzz_xyyyyyyz_1[j] + 0.5 * fl1_fxn * tg_yzz_yyyyyyz_1[j];

                    tg_xyzz_xyyyyyzz_0[j] = pb_x * tg_yzz_xyyyyyzz_0[j] + wp_x[j] * tg_yzz_xyyyyyzz_1[j] + 0.5 * fl1_fxn * tg_yzz_yyyyyzz_1[j];

                    tg_xyzz_xyyyyzzz_0[j] = pb_x * tg_yzz_xyyyyzzz_0[j] + wp_x[j] * tg_yzz_xyyyyzzz_1[j] + 0.5 * fl1_fxn * tg_yzz_yyyyzzz_1[j];

                    tg_xyzz_xyyyzzzz_0[j] = pb_x * tg_yzz_xyyyzzzz_0[j] + wp_x[j] * tg_yzz_xyyyzzzz_1[j] + 0.5 * fl1_fxn * tg_yzz_yyyzzzz_1[j];

                    tg_xyzz_xyyzzzzz_0[j] = pb_x * tg_yzz_xyyzzzzz_0[j] + wp_x[j] * tg_yzz_xyyzzzzz_1[j] + 0.5 * fl1_fxn * tg_yzz_yyzzzzz_1[j];

                    tg_xyzz_xyzzzzzz_0[j] = pb_x * tg_yzz_xyzzzzzz_0[j] + wp_x[j] * tg_yzz_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yzz_yzzzzzz_1[j];

                    tg_xyzz_xzzzzzzz_0[j] = pb_x * tg_yzz_xzzzzzzz_0[j] + wp_x[j] * tg_yzz_xzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yzz_zzzzzzz_1[j];

                    tg_xyzz_yyyyyyyy_0[j] = pb_x * tg_yzz_yyyyyyyy_0[j] + wp_x[j] * tg_yzz_yyyyyyyy_1[j];

                    tg_xyzz_yyyyyyyz_0[j] = pb_x * tg_yzz_yyyyyyyz_0[j] + wp_x[j] * tg_yzz_yyyyyyyz_1[j];

                    tg_xyzz_yyyyyyzz_0[j] = pb_x * tg_yzz_yyyyyyzz_0[j] + wp_x[j] * tg_yzz_yyyyyyzz_1[j];

                    tg_xyzz_yyyyyzzz_0[j] = pb_x * tg_yzz_yyyyyzzz_0[j] + wp_x[j] * tg_yzz_yyyyyzzz_1[j];

                    tg_xyzz_yyyyzzzz_0[j] = pb_x * tg_yzz_yyyyzzzz_0[j] + wp_x[j] * tg_yzz_yyyyzzzz_1[j];

                    tg_xyzz_yyyzzzzz_0[j] = pb_x * tg_yzz_yyyzzzzz_0[j] + wp_x[j] * tg_yzz_yyyzzzzz_1[j];

                    tg_xyzz_yyzzzzzz_0[j] = pb_x * tg_yzz_yyzzzzzz_0[j] + wp_x[j] * tg_yzz_yyzzzzzz_1[j];

                    tg_xyzz_yzzzzzzz_0[j] = pb_x * tg_yzz_yzzzzzzz_0[j] + wp_x[j] * tg_yzz_yzzzzzzz_1[j];

                    tg_xyzz_zzzzzzzz_0[j] = pb_x * tg_yzz_zzzzzzzz_0[j] + wp_x[j] * tg_yzz_zzzzzzzz_1[j];

                    tg_xzzz_xxxxxxxx_0[j] = pb_x * tg_zzz_xxxxxxxx_0[j] + wp_x[j] * tg_zzz_xxxxxxxx_1[j] + 4.0 * fl1_fxn * tg_zzz_xxxxxxx_1[j];

                    tg_xzzz_xxxxxxxy_0[j] = pb_x * tg_zzz_xxxxxxxy_0[j] + wp_x[j] * tg_zzz_xxxxxxxy_1[j] + 3.5 * fl1_fxn * tg_zzz_xxxxxxy_1[j];

                    tg_xzzz_xxxxxxxz_0[j] = pb_x * tg_zzz_xxxxxxxz_0[j] + wp_x[j] * tg_zzz_xxxxxxxz_1[j] + 3.5 * fl1_fxn * tg_zzz_xxxxxxz_1[j];

                    tg_xzzz_xxxxxxyy_0[j] = pb_x * tg_zzz_xxxxxxyy_0[j] + wp_x[j] * tg_zzz_xxxxxxyy_1[j] + 3.0 * fl1_fxn * tg_zzz_xxxxxyy_1[j];

                    tg_xzzz_xxxxxxyz_0[j] = pb_x * tg_zzz_xxxxxxyz_0[j] + wp_x[j] * tg_zzz_xxxxxxyz_1[j] + 3.0 * fl1_fxn * tg_zzz_xxxxxyz_1[j];

                    tg_xzzz_xxxxxxzz_0[j] = pb_x * tg_zzz_xxxxxxzz_0[j] + wp_x[j] * tg_zzz_xxxxxxzz_1[j] + 3.0 * fl1_fxn * tg_zzz_xxxxxzz_1[j];

                    tg_xzzz_xxxxxyyy_0[j] = pb_x * tg_zzz_xxxxxyyy_0[j] + wp_x[j] * tg_zzz_xxxxxyyy_1[j] + 2.5 * fl1_fxn * tg_zzz_xxxxyyy_1[j];

                    tg_xzzz_xxxxxyyz_0[j] = pb_x * tg_zzz_xxxxxyyz_0[j] + wp_x[j] * tg_zzz_xxxxxyyz_1[j] + 2.5 * fl1_fxn * tg_zzz_xxxxyyz_1[j];

                    tg_xzzz_xxxxxyzz_0[j] = pb_x * tg_zzz_xxxxxyzz_0[j] + wp_x[j] * tg_zzz_xxxxxyzz_1[j] + 2.5 * fl1_fxn * tg_zzz_xxxxyzz_1[j];

                    tg_xzzz_xxxxxzzz_0[j] = pb_x * tg_zzz_xxxxxzzz_0[j] + wp_x[j] * tg_zzz_xxxxxzzz_1[j] + 2.5 * fl1_fxn * tg_zzz_xxxxzzz_1[j];

                    tg_xzzz_xxxxyyyy_0[j] = pb_x * tg_zzz_xxxxyyyy_0[j] + wp_x[j] * tg_zzz_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_zzz_xxxyyyy_1[j];

                    tg_xzzz_xxxxyyyz_0[j] = pb_x * tg_zzz_xxxxyyyz_0[j] + wp_x[j] * tg_zzz_xxxxyyyz_1[j] + 2.0 * fl1_fxn * tg_zzz_xxxyyyz_1[j];

                    tg_xzzz_xxxxyyzz_0[j] = pb_x * tg_zzz_xxxxyyzz_0[j] + wp_x[j] * tg_zzz_xxxxyyzz_1[j] + 2.0 * fl1_fxn * tg_zzz_xxxyyzz_1[j];

                    tg_xzzz_xxxxyzzz_0[j] = pb_x * tg_zzz_xxxxyzzz_0[j] + wp_x[j] * tg_zzz_xxxxyzzz_1[j] + 2.0 * fl1_fxn * tg_zzz_xxxyzzz_1[j];

                    tg_xzzz_xxxxzzzz_0[j] = pb_x * tg_zzz_xxxxzzzz_0[j] + wp_x[j] * tg_zzz_xxxxzzzz_1[j] + 2.0 * fl1_fxn * tg_zzz_xxxzzzz_1[j];

                    tg_xzzz_xxxyyyyy_0[j] = pb_x * tg_zzz_xxxyyyyy_0[j] + wp_x[j] * tg_zzz_xxxyyyyy_1[j] + 1.5 * fl1_fxn * tg_zzz_xxyyyyy_1[j];

                    tg_xzzz_xxxyyyyz_0[j] = pb_x * tg_zzz_xxxyyyyz_0[j] + wp_x[j] * tg_zzz_xxxyyyyz_1[j] + 1.5 * fl1_fxn * tg_zzz_xxyyyyz_1[j];

                    tg_xzzz_xxxyyyzz_0[j] = pb_x * tg_zzz_xxxyyyzz_0[j] + wp_x[j] * tg_zzz_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_zzz_xxyyyzz_1[j];

                    tg_xzzz_xxxyyzzz_0[j] = pb_x * tg_zzz_xxxyyzzz_0[j] + wp_x[j] * tg_zzz_xxxyyzzz_1[j] + 1.5 * fl1_fxn * tg_zzz_xxyyzzz_1[j];

                    tg_xzzz_xxxyzzzz_0[j] = pb_x * tg_zzz_xxxyzzzz_0[j] + wp_x[j] * tg_zzz_xxxyzzzz_1[j] + 1.5 * fl1_fxn * tg_zzz_xxyzzzz_1[j];

                    tg_xzzz_xxxzzzzz_0[j] = pb_x * tg_zzz_xxxzzzzz_0[j] + wp_x[j] * tg_zzz_xxxzzzzz_1[j] + 1.5 * fl1_fxn * tg_zzz_xxzzzzz_1[j];

                    tg_xzzz_xxyyyyyy_0[j] = pb_x * tg_zzz_xxyyyyyy_0[j] + wp_x[j] * tg_zzz_xxyyyyyy_1[j] + fl1_fxn * tg_zzz_xyyyyyy_1[j];

                    tg_xzzz_xxyyyyyz_0[j] = pb_x * tg_zzz_xxyyyyyz_0[j] + wp_x[j] * tg_zzz_xxyyyyyz_1[j] + fl1_fxn * tg_zzz_xyyyyyz_1[j];

                    tg_xzzz_xxyyyyzz_0[j] = pb_x * tg_zzz_xxyyyyzz_0[j] + wp_x[j] * tg_zzz_xxyyyyzz_1[j] + fl1_fxn * tg_zzz_xyyyyzz_1[j];

                    tg_xzzz_xxyyyzzz_0[j] = pb_x * tg_zzz_xxyyyzzz_0[j] + wp_x[j] * tg_zzz_xxyyyzzz_1[j] + fl1_fxn * tg_zzz_xyyyzzz_1[j];

                    tg_xzzz_xxyyzzzz_0[j] = pb_x * tg_zzz_xxyyzzzz_0[j] + wp_x[j] * tg_zzz_xxyyzzzz_1[j] + fl1_fxn * tg_zzz_xyyzzzz_1[j];

                    tg_xzzz_xxyzzzzz_0[j] = pb_x * tg_zzz_xxyzzzzz_0[j] + wp_x[j] * tg_zzz_xxyzzzzz_1[j] + fl1_fxn * tg_zzz_xyzzzzz_1[j];

                    tg_xzzz_xxzzzzzz_0[j] = pb_x * tg_zzz_xxzzzzzz_0[j] + wp_x[j] * tg_zzz_xxzzzzzz_1[j] + fl1_fxn * tg_zzz_xzzzzzz_1[j];

                    tg_xzzz_xyyyyyyy_0[j] = pb_x * tg_zzz_xyyyyyyy_0[j] + wp_x[j] * tg_zzz_xyyyyyyy_1[j] + 0.5 * fl1_fxn * tg_zzz_yyyyyyy_1[j];

                    tg_xzzz_xyyyyyyz_0[j] = pb_x * tg_zzz_xyyyyyyz_0[j] + wp_x[j] * tg_zzz_xyyyyyyz_1[j] + 0.5 * fl1_fxn * tg_zzz_yyyyyyz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSL_435_483(      CMemBlock2D<double>& primBuffer,
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

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {4, -1, -1, -1},
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_3_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_2_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_2_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
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

                auto tg_yyy_xxxxxxxx_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 270); 

                auto tg_yyy_xxxxxxxy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 271); 

                auto tg_yyy_xxxxxxxz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 272); 

                auto tg_yyy_xxxxxxyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 273); 

                auto tg_yyy_xxxxxxyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 274); 

                auto tg_yyy_xxxxxxzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 275); 

                auto tg_yyy_xxxxxyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 276); 

                auto tg_yyy_xxxxxyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 277); 

                auto tg_yyy_xxxxxyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 278); 

                auto tg_yyy_xxxxxzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 279); 

                auto tg_yyy_xxxxyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 280); 

                auto tg_yyy_xxxxyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 281); 

                auto tg_yyy_xxxxyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 282); 

                auto tg_yyy_xxxxyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 283); 

                auto tg_yyy_xxxxzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 284); 

                auto tg_yyy_xxxyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 285); 

                auto tg_yyy_xxxyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 286); 

                auto tg_yyy_xxxyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 287); 

                auto tg_yyy_xxxyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 288); 

                auto tg_yyy_xxxyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 289); 

                auto tg_yyy_xxxzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 290); 

                auto tg_yyy_xxyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 291); 

                auto tg_yyy_xxyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 292); 

                auto tg_yyy_xxyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 293); 

                auto tg_yyy_xxyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 294); 

                auto tg_yyy_xxyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 295); 

                auto tg_yyy_xxyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 296); 

                auto tg_yyy_xxzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 297); 

                auto tg_yyy_xyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 298); 

                auto tg_yyy_xyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 299); 

                auto tg_yyy_xyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 300); 

                auto tg_yyy_xyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 301); 

                auto tg_yyy_xyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 302); 

                auto tg_zzz_xyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 435); 

                auto tg_zzz_xyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 436); 

                auto tg_zzz_xyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 437); 

                auto tg_zzz_xyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 438); 

                auto tg_zzz_xyzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 439); 

                auto tg_zzz_xzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 440); 

                auto tg_zzz_yyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 441); 

                auto tg_zzz_yyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 442); 

                auto tg_zzz_yyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 443); 

                auto tg_zzz_yyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 444); 

                auto tg_zzz_yyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 445); 

                auto tg_zzz_yyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 446); 

                auto tg_zzz_yyzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 447); 

                auto tg_zzz_yzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 448); 

                auto tg_zzz_zzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 449); 

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

                auto tg_zzz_yyyyyzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 354); 

                auto tg_zzz_yyyyzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 355); 

                auto tg_zzz_yyyzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 356); 

                auto tg_zzz_yyzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 357); 

                auto tg_zzz_yzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 358); 

                auto tg_zzz_zzzzzzz_1 = primBuffer.data(pidx_g_3_7_m1 + 360 * idx + 359); 

                // set up pointers to integrals

                auto tg_xzzz_xyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 435); 

                auto tg_xzzz_xyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 436); 

                auto tg_xzzz_xyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 437); 

                auto tg_xzzz_xyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 438); 

                auto tg_xzzz_xyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 439); 

                auto tg_xzzz_xzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 440); 

                auto tg_xzzz_yyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 441); 

                auto tg_xzzz_yyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 442); 

                auto tg_xzzz_yyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 443); 

                auto tg_xzzz_yyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 444); 

                auto tg_xzzz_yyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 445); 

                auto tg_xzzz_yyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 446); 

                auto tg_xzzz_yyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 447); 

                auto tg_xzzz_yzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 448); 

                auto tg_xzzz_zzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 449); 

                auto tg_yyyy_xxxxxxxx_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 450); 

                auto tg_yyyy_xxxxxxxy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 451); 

                auto tg_yyyy_xxxxxxxz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 452); 

                auto tg_yyyy_xxxxxxyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 453); 

                auto tg_yyyy_xxxxxxyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 454); 

                auto tg_yyyy_xxxxxxzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 455); 

                auto tg_yyyy_xxxxxyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 456); 

                auto tg_yyyy_xxxxxyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 457); 

                auto tg_yyyy_xxxxxyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 458); 

                auto tg_yyyy_xxxxxzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 459); 

                auto tg_yyyy_xxxxyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 460); 

                auto tg_yyyy_xxxxyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 461); 

                auto tg_yyyy_xxxxyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 462); 

                auto tg_yyyy_xxxxyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 463); 

                auto tg_yyyy_xxxxzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 464); 

                auto tg_yyyy_xxxyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 465); 

                auto tg_yyyy_xxxyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 466); 

                auto tg_yyyy_xxxyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 467); 

                auto tg_yyyy_xxxyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 468); 

                auto tg_yyyy_xxxyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 469); 

                auto tg_yyyy_xxxzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 470); 

                auto tg_yyyy_xxyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 471); 

                auto tg_yyyy_xxyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 472); 

                auto tg_yyyy_xxyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 473); 

                auto tg_yyyy_xxyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 474); 

                auto tg_yyyy_xxyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 475); 

                auto tg_yyyy_xxyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 476); 

                auto tg_yyyy_xxzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 477); 

                auto tg_yyyy_xyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 478); 

                auto tg_yyyy_xyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 479); 

                auto tg_yyyy_xyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 480); 

                auto tg_yyyy_xyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 481); 

                auto tg_yyyy_xyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 482); 

                // Batch of Integrals (435,483)

                #pragma omp simd aligned(fxn, fza, tg_xzzz_xyyyyyzz_0, tg_xzzz_xyyyyzzz_0, tg_xzzz_xyyyzzzz_0, \
                                         tg_xzzz_xyyzzzzz_0, tg_xzzz_xyzzzzzz_0, tg_xzzz_xzzzzzzz_0, tg_xzzz_yyyyyyyy_0, \
                                         tg_xzzz_yyyyyyyz_0, tg_xzzz_yyyyyyzz_0, tg_xzzz_yyyyyzzz_0, tg_xzzz_yyyyzzzz_0, \
                                         tg_xzzz_yyyzzzzz_0, tg_xzzz_yyzzzzzz_0, tg_xzzz_yzzzzzzz_0, tg_xzzz_zzzzzzzz_0, \
                                         tg_yy_xxxxxxxx_0, tg_yy_xxxxxxxx_1, tg_yy_xxxxxxxy_0, tg_yy_xxxxxxxy_1, \
                                         tg_yy_xxxxxxxz_0, tg_yy_xxxxxxxz_1, tg_yy_xxxxxxyy_0, tg_yy_xxxxxxyy_1, \
                                         tg_yy_xxxxxxyz_0, tg_yy_xxxxxxyz_1, tg_yy_xxxxxxzz_0, tg_yy_xxxxxxzz_1, \
                                         tg_yy_xxxxxyyy_0, tg_yy_xxxxxyyy_1, tg_yy_xxxxxyyz_0, tg_yy_xxxxxyyz_1, \
                                         tg_yy_xxxxxyzz_0, tg_yy_xxxxxyzz_1, tg_yy_xxxxxzzz_0, tg_yy_xxxxxzzz_1, \
                                         tg_yy_xxxxyyyy_0, tg_yy_xxxxyyyy_1, tg_yy_xxxxyyyz_0, tg_yy_xxxxyyyz_1, \
                                         tg_yy_xxxxyyzz_0, tg_yy_xxxxyyzz_1, tg_yy_xxxxyzzz_0, tg_yy_xxxxyzzz_1, \
                                         tg_yy_xxxxzzzz_0, tg_yy_xxxxzzzz_1, tg_yy_xxxyyyyy_0, tg_yy_xxxyyyyy_1, \
                                         tg_yy_xxxyyyyz_0, tg_yy_xxxyyyyz_1, tg_yy_xxxyyyzz_0, tg_yy_xxxyyyzz_1, \
                                         tg_yy_xxxyyzzz_0, tg_yy_xxxyyzzz_1, tg_yy_xxxyzzzz_0, tg_yy_xxxyzzzz_1, \
                                         tg_yy_xxxzzzzz_0, tg_yy_xxxzzzzz_1, tg_yy_xxyyyyyy_0, tg_yy_xxyyyyyy_1, \
                                         tg_yy_xxyyyyyz_0, tg_yy_xxyyyyyz_1, tg_yy_xxyyyyzz_0, tg_yy_xxyyyyzz_1, \
                                         tg_yy_xxyyyzzz_0, tg_yy_xxyyyzzz_1, tg_yy_xxyyzzzz_0, tg_yy_xxyyzzzz_1, \
                                         tg_yy_xxyzzzzz_0, tg_yy_xxyzzzzz_1, tg_yy_xxzzzzzz_0, tg_yy_xxzzzzzz_1, \
                                         tg_yy_xyyyyyyy_0, tg_yy_xyyyyyyy_1, tg_yy_xyyyyyyz_0, tg_yy_xyyyyyyz_1, \
                                         tg_yy_xyyyyyzz_0, tg_yy_xyyyyyzz_1, tg_yy_xyyyyzzz_0, tg_yy_xyyyyzzz_1, \
                                         tg_yy_xyyyzzzz_0, tg_yy_xyyyzzzz_1, tg_yyy_xxxxxxx_1, tg_yyy_xxxxxxxx_0, \
                                         tg_yyy_xxxxxxxx_1, tg_yyy_xxxxxxxy_0, tg_yyy_xxxxxxxy_1, tg_yyy_xxxxxxxz_0, \
                                         tg_yyy_xxxxxxxz_1, tg_yyy_xxxxxxy_1, tg_yyy_xxxxxxyy_0, tg_yyy_xxxxxxyy_1, \
                                         tg_yyy_xxxxxxyz_0, tg_yyy_xxxxxxyz_1, tg_yyy_xxxxxxz_1, tg_yyy_xxxxxxzz_0, \
                                         tg_yyy_xxxxxxzz_1, tg_yyy_xxxxxyy_1, tg_yyy_xxxxxyyy_0, tg_yyy_xxxxxyyy_1, \
                                         tg_yyy_xxxxxyyz_0, tg_yyy_xxxxxyyz_1, tg_yyy_xxxxxyz_1, tg_yyy_xxxxxyzz_0, \
                                         tg_yyy_xxxxxyzz_1, tg_yyy_xxxxxzz_1, tg_yyy_xxxxxzzz_0, tg_yyy_xxxxxzzz_1, \
                                         tg_yyy_xxxxyyy_1, tg_yyy_xxxxyyyy_0, tg_yyy_xxxxyyyy_1, tg_yyy_xxxxyyyz_0, \
                                         tg_yyy_xxxxyyyz_1, tg_yyy_xxxxyyz_1, tg_yyy_xxxxyyzz_0, tg_yyy_xxxxyyzz_1, \
                                         tg_yyy_xxxxyzz_1, tg_yyy_xxxxyzzz_0, tg_yyy_xxxxyzzz_1, tg_yyy_xxxxzzz_1, \
                                         tg_yyy_xxxxzzzz_0, tg_yyy_xxxxzzzz_1, tg_yyy_xxxyyyy_1, tg_yyy_xxxyyyyy_0, \
                                         tg_yyy_xxxyyyyy_1, tg_yyy_xxxyyyyz_0, tg_yyy_xxxyyyyz_1, tg_yyy_xxxyyyz_1, \
                                         tg_yyy_xxxyyyzz_0, tg_yyy_xxxyyyzz_1, tg_yyy_xxxyyzz_1, tg_yyy_xxxyyzzz_0, \
                                         tg_yyy_xxxyyzzz_1, tg_yyy_xxxyzzz_1, tg_yyy_xxxyzzzz_0, tg_yyy_xxxyzzzz_1, \
                                         tg_yyy_xxxzzzz_1, tg_yyy_xxxzzzzz_0, tg_yyy_xxxzzzzz_1, tg_yyy_xxyyyyy_1, \
                                         tg_yyy_xxyyyyyy_0, tg_yyy_xxyyyyyy_1, tg_yyy_xxyyyyyz_0, tg_yyy_xxyyyyyz_1, \
                                         tg_yyy_xxyyyyz_1, tg_yyy_xxyyyyzz_0, tg_yyy_xxyyyyzz_1, tg_yyy_xxyyyzz_1, \
                                         tg_yyy_xxyyyzzz_0, tg_yyy_xxyyyzzz_1, tg_yyy_xxyyzzz_1, tg_yyy_xxyyzzzz_0, \
                                         tg_yyy_xxyyzzzz_1, tg_yyy_xxyzzzz_1, tg_yyy_xxyzzzzz_0, tg_yyy_xxyzzzzz_1, \
                                         tg_yyy_xxzzzzz_1, tg_yyy_xxzzzzzz_0, tg_yyy_xxzzzzzz_1, tg_yyy_xyyyyyy_1, \
                                         tg_yyy_xyyyyyyy_0, tg_yyy_xyyyyyyy_1, tg_yyy_xyyyyyyz_0, tg_yyy_xyyyyyyz_1, \
                                         tg_yyy_xyyyyyz_1, tg_yyy_xyyyyyzz_0, tg_yyy_xyyyyyzz_1, tg_yyy_xyyyyzz_1, \
                                         tg_yyy_xyyyyzzz_0, tg_yyy_xyyyyzzz_1, tg_yyy_xyyyzzz_1, tg_yyy_xyyyzzzz_0, \
                                         tg_yyy_xyyyzzzz_1, tg_yyy_xyyzzzz_1, tg_yyyy_xxxxxxxx_0, tg_yyyy_xxxxxxxy_0, \
                                         tg_yyyy_xxxxxxxz_0, tg_yyyy_xxxxxxyy_0, tg_yyyy_xxxxxxyz_0, tg_yyyy_xxxxxxzz_0, \
                                         tg_yyyy_xxxxxyyy_0, tg_yyyy_xxxxxyyz_0, tg_yyyy_xxxxxyzz_0, tg_yyyy_xxxxxzzz_0, \
                                         tg_yyyy_xxxxyyyy_0, tg_yyyy_xxxxyyyz_0, tg_yyyy_xxxxyyzz_0, tg_yyyy_xxxxyzzz_0, \
                                         tg_yyyy_xxxxzzzz_0, tg_yyyy_xxxyyyyy_0, tg_yyyy_xxxyyyyz_0, tg_yyyy_xxxyyyzz_0, \
                                         tg_yyyy_xxxyyzzz_0, tg_yyyy_xxxyzzzz_0, tg_yyyy_xxxzzzzz_0, tg_yyyy_xxyyyyyy_0, \
                                         tg_yyyy_xxyyyyyz_0, tg_yyyy_xxyyyyzz_0, tg_yyyy_xxyyyzzz_0, tg_yyyy_xxyyzzzz_0, \
                                         tg_yyyy_xxyzzzzz_0, tg_yyyy_xxzzzzzz_0, tg_yyyy_xyyyyyyy_0, tg_yyyy_xyyyyyyz_0, \
                                         tg_yyyy_xyyyyyzz_0, tg_yyyy_xyyyyzzz_0, tg_yyyy_xyyyzzzz_0, tg_zzz_xyyyyyzz_0, \
                                         tg_zzz_xyyyyyzz_1, tg_zzz_xyyyyzzz_0, tg_zzz_xyyyyzzz_1, tg_zzz_xyyyzzzz_0, \
                                         tg_zzz_xyyyzzzz_1, tg_zzz_xyyzzzzz_0, tg_zzz_xyyzzzzz_1, tg_zzz_xyzzzzzz_0, \
                                         tg_zzz_xyzzzzzz_1, tg_zzz_xzzzzzzz_0, tg_zzz_xzzzzzzz_1, tg_zzz_yyyyyyyy_0, \
                                         tg_zzz_yyyyyyyy_1, tg_zzz_yyyyyyyz_0, tg_zzz_yyyyyyyz_1, tg_zzz_yyyyyyzz_0, \
                                         tg_zzz_yyyyyyzz_1, tg_zzz_yyyyyzz_1, tg_zzz_yyyyyzzz_0, tg_zzz_yyyyyzzz_1, \
                                         tg_zzz_yyyyzzz_1, tg_zzz_yyyyzzzz_0, tg_zzz_yyyyzzzz_1, tg_zzz_yyyzzzz_1, \
                                         tg_zzz_yyyzzzzz_0, tg_zzz_yyyzzzzz_1, tg_zzz_yyzzzzz_1, tg_zzz_yyzzzzzz_0, \
                                         tg_zzz_yyzzzzzz_1, tg_zzz_yzzzzzz_1, tg_zzz_yzzzzzzz_0, tg_zzz_yzzzzzzz_1, \
                                         tg_zzz_zzzzzzz_1, tg_zzz_zzzzzzzz_0, tg_zzz_zzzzzzzz_1, wp_x, wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xzzz_xyyyyyzz_0[j] = pb_x * tg_zzz_xyyyyyzz_0[j] + wp_x[j] * tg_zzz_xyyyyyzz_1[j] + 0.5 * fl1_fxn * tg_zzz_yyyyyzz_1[j];

                    tg_xzzz_xyyyyzzz_0[j] = pb_x * tg_zzz_xyyyyzzz_0[j] + wp_x[j] * tg_zzz_xyyyyzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_yyyyzzz_1[j];

                    tg_xzzz_xyyyzzzz_0[j] = pb_x * tg_zzz_xyyyzzzz_0[j] + wp_x[j] * tg_zzz_xyyyzzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_yyyzzzz_1[j];

                    tg_xzzz_xyyzzzzz_0[j] = pb_x * tg_zzz_xyyzzzzz_0[j] + wp_x[j] * tg_zzz_xyyzzzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_yyzzzzz_1[j];

                    tg_xzzz_xyzzzzzz_0[j] = pb_x * tg_zzz_xyzzzzzz_0[j] + wp_x[j] * tg_zzz_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_yzzzzzz_1[j];

                    tg_xzzz_xzzzzzzz_0[j] = pb_x * tg_zzz_xzzzzzzz_0[j] + wp_x[j] * tg_zzz_xzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_zzzzzzz_1[j];

                    tg_xzzz_yyyyyyyy_0[j] = pb_x * tg_zzz_yyyyyyyy_0[j] + wp_x[j] * tg_zzz_yyyyyyyy_1[j];

                    tg_xzzz_yyyyyyyz_0[j] = pb_x * tg_zzz_yyyyyyyz_0[j] + wp_x[j] * tg_zzz_yyyyyyyz_1[j];

                    tg_xzzz_yyyyyyzz_0[j] = pb_x * tg_zzz_yyyyyyzz_0[j] + wp_x[j] * tg_zzz_yyyyyyzz_1[j];

                    tg_xzzz_yyyyyzzz_0[j] = pb_x * tg_zzz_yyyyyzzz_0[j] + wp_x[j] * tg_zzz_yyyyyzzz_1[j];

                    tg_xzzz_yyyyzzzz_0[j] = pb_x * tg_zzz_yyyyzzzz_0[j] + wp_x[j] * tg_zzz_yyyyzzzz_1[j];

                    tg_xzzz_yyyzzzzz_0[j] = pb_x * tg_zzz_yyyzzzzz_0[j] + wp_x[j] * tg_zzz_yyyzzzzz_1[j];

                    tg_xzzz_yyzzzzzz_0[j] = pb_x * tg_zzz_yyzzzzzz_0[j] + wp_x[j] * tg_zzz_yyzzzzzz_1[j];

                    tg_xzzz_yzzzzzzz_0[j] = pb_x * tg_zzz_yzzzzzzz_0[j] + wp_x[j] * tg_zzz_yzzzzzzz_1[j];

                    tg_xzzz_zzzzzzzz_0[j] = pb_x * tg_zzz_zzzzzzzz_0[j] + wp_x[j] * tg_zzz_zzzzzzzz_1[j];

                    tg_yyyy_xxxxxxxx_0[j] = pb_y * tg_yyy_xxxxxxxx_0[j] + wp_y[j] * tg_yyy_xxxxxxxx_1[j] + 1.5 * fl1_fx * tg_yy_xxxxxxxx_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxxxxxxx_1[j];

                    tg_yyyy_xxxxxxxy_0[j] = pb_y * tg_yyy_xxxxxxxy_0[j] + wp_y[j] * tg_yyy_xxxxxxxy_1[j] + 1.5 * fl1_fx * tg_yy_xxxxxxxy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxxxxxxy_1[j] + 0.5 * fl1_fxn * tg_yyy_xxxxxxx_1[j];

                    tg_yyyy_xxxxxxxz_0[j] = pb_y * tg_yyy_xxxxxxxz_0[j] + wp_y[j] * tg_yyy_xxxxxxxz_1[j] + 1.5 * fl1_fx * tg_yy_xxxxxxxz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxxxxxxz_1[j];

                    tg_yyyy_xxxxxxyy_0[j] = pb_y * tg_yyy_xxxxxxyy_0[j] + wp_y[j] * tg_yyy_xxxxxxyy_1[j] + 1.5 * fl1_fx * tg_yy_xxxxxxyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxxxxxyy_1[j] + fl1_fxn * tg_yyy_xxxxxxy_1[j];

                    tg_yyyy_xxxxxxyz_0[j] = pb_y * tg_yyy_xxxxxxyz_0[j] + wp_y[j] * tg_yyy_xxxxxxyz_1[j] + 1.5 * fl1_fx * tg_yy_xxxxxxyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxxxxxyz_1[j] + 0.5 * fl1_fxn * tg_yyy_xxxxxxz_1[j];

                    tg_yyyy_xxxxxxzz_0[j] = pb_y * tg_yyy_xxxxxxzz_0[j] + wp_y[j] * tg_yyy_xxxxxxzz_1[j] + 1.5 * fl1_fx * tg_yy_xxxxxxzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxxxxxzz_1[j];

                    tg_yyyy_xxxxxyyy_0[j] = pb_y * tg_yyy_xxxxxyyy_0[j] + wp_y[j] * tg_yyy_xxxxxyyy_1[j] + 1.5 * fl1_fx * tg_yy_xxxxxyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxxxxyyy_1[j] + 1.5 * fl1_fxn * tg_yyy_xxxxxyy_1[j];

                    tg_yyyy_xxxxxyyz_0[j] = pb_y * tg_yyy_xxxxxyyz_0[j] + wp_y[j] * tg_yyy_xxxxxyyz_1[j] + 1.5 * fl1_fx * tg_yy_xxxxxyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxxxxyyz_1[j] + fl1_fxn * tg_yyy_xxxxxyz_1[j];

                    tg_yyyy_xxxxxyzz_0[j] = pb_y * tg_yyy_xxxxxyzz_0[j] + wp_y[j] * tg_yyy_xxxxxyzz_1[j] + 1.5 * fl1_fx * tg_yy_xxxxxyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxxxxyzz_1[j] + 0.5 * fl1_fxn * tg_yyy_xxxxxzz_1[j];

                    tg_yyyy_xxxxxzzz_0[j] = pb_y * tg_yyy_xxxxxzzz_0[j] + wp_y[j] * tg_yyy_xxxxxzzz_1[j] + 1.5 * fl1_fx * tg_yy_xxxxxzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxxxxzzz_1[j];

                    tg_yyyy_xxxxyyyy_0[j] = pb_y * tg_yyy_xxxxyyyy_0[j] + wp_y[j] * tg_yyy_xxxxyyyy_1[j] + 1.5 * fl1_fx * tg_yy_xxxxyyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_yyy_xxxxyyy_1[j];

                    tg_yyyy_xxxxyyyz_0[j] = pb_y * tg_yyy_xxxxyyyz_0[j] + wp_y[j] * tg_yyy_xxxxyyyz_1[j] + 1.5 * fl1_fx * tg_yy_xxxxyyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxxxyyyz_1[j] + 1.5 * fl1_fxn * tg_yyy_xxxxyyz_1[j];

                    tg_yyyy_xxxxyyzz_0[j] = pb_y * tg_yyy_xxxxyyzz_0[j] + wp_y[j] * tg_yyy_xxxxyyzz_1[j] + 1.5 * fl1_fx * tg_yy_xxxxyyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxxxyyzz_1[j] + fl1_fxn * tg_yyy_xxxxyzz_1[j];

                    tg_yyyy_xxxxyzzz_0[j] = pb_y * tg_yyy_xxxxyzzz_0[j] + wp_y[j] * tg_yyy_xxxxyzzz_1[j] + 1.5 * fl1_fx * tg_yy_xxxxyzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxxxyzzz_1[j] + 0.5 * fl1_fxn * tg_yyy_xxxxzzz_1[j];

                    tg_yyyy_xxxxzzzz_0[j] = pb_y * tg_yyy_xxxxzzzz_0[j] + wp_y[j] * tg_yyy_xxxxzzzz_1[j] + 1.5 * fl1_fx * tg_yy_xxxxzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxxxzzzz_1[j];

                    tg_yyyy_xxxyyyyy_0[j] = pb_y * tg_yyy_xxxyyyyy_0[j] + wp_y[j] * tg_yyy_xxxyyyyy_1[j] + 1.5 * fl1_fx * tg_yy_xxxyyyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxxyyyyy_1[j] + 2.5 * fl1_fxn * tg_yyy_xxxyyyy_1[j];

                    tg_yyyy_xxxyyyyz_0[j] = pb_y * tg_yyy_xxxyyyyz_0[j] + wp_y[j] * tg_yyy_xxxyyyyz_1[j] + 1.5 * fl1_fx * tg_yy_xxxyyyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxxyyyyz_1[j] + 2.0 * fl1_fxn * tg_yyy_xxxyyyz_1[j];

                    tg_yyyy_xxxyyyzz_0[j] = pb_y * tg_yyy_xxxyyyzz_0[j] + wp_y[j] * tg_yyy_xxxyyyzz_1[j] + 1.5 * fl1_fx * tg_yy_xxxyyyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_yyy_xxxyyzz_1[j];

                    tg_yyyy_xxxyyzzz_0[j] = pb_y * tg_yyy_xxxyyzzz_0[j] + wp_y[j] * tg_yyy_xxxyyzzz_1[j] + 1.5 * fl1_fx * tg_yy_xxxyyzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxxyyzzz_1[j] + fl1_fxn * tg_yyy_xxxyzzz_1[j];

                    tg_yyyy_xxxyzzzz_0[j] = pb_y * tg_yyy_xxxyzzzz_0[j] + wp_y[j] * tg_yyy_xxxyzzzz_1[j] + 1.5 * fl1_fx * tg_yy_xxxyzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxxyzzzz_1[j] + 0.5 * fl1_fxn * tg_yyy_xxxzzzz_1[j];

                    tg_yyyy_xxxzzzzz_0[j] = pb_y * tg_yyy_xxxzzzzz_0[j] + wp_y[j] * tg_yyy_xxxzzzzz_1[j] + 1.5 * fl1_fx * tg_yy_xxxzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxxzzzzz_1[j];

                    tg_yyyy_xxyyyyyy_0[j] = pb_y * tg_yyy_xxyyyyyy_0[j] + wp_y[j] * tg_yyy_xxyyyyyy_1[j] + 1.5 * fl1_fx * tg_yy_xxyyyyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxyyyyyy_1[j] + 3.0 * fl1_fxn * tg_yyy_xxyyyyy_1[j];

                    tg_yyyy_xxyyyyyz_0[j] = pb_y * tg_yyy_xxyyyyyz_0[j] + wp_y[j] * tg_yyy_xxyyyyyz_1[j] + 1.5 * fl1_fx * tg_yy_xxyyyyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxyyyyyz_1[j] + 2.5 * fl1_fxn * tg_yyy_xxyyyyz_1[j];

                    tg_yyyy_xxyyyyzz_0[j] = pb_y * tg_yyy_xxyyyyzz_0[j] + wp_y[j] * tg_yyy_xxyyyyzz_1[j] + 1.5 * fl1_fx * tg_yy_xxyyyyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxyyyyzz_1[j] + 2.0 * fl1_fxn * tg_yyy_xxyyyzz_1[j];

                    tg_yyyy_xxyyyzzz_0[j] = pb_y * tg_yyy_xxyyyzzz_0[j] + wp_y[j] * tg_yyy_xxyyyzzz_1[j] + 1.5 * fl1_fx * tg_yy_xxyyyzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxyyyzzz_1[j] + 1.5 * fl1_fxn * tg_yyy_xxyyzzz_1[j];

                    tg_yyyy_xxyyzzzz_0[j] = pb_y * tg_yyy_xxyyzzzz_0[j] + wp_y[j] * tg_yyy_xxyyzzzz_1[j] + 1.5 * fl1_fx * tg_yy_xxyyzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxyyzzzz_1[j] + fl1_fxn * tg_yyy_xxyzzzz_1[j];

                    tg_yyyy_xxyzzzzz_0[j] = pb_y * tg_yyy_xxyzzzzz_0[j] + wp_y[j] * tg_yyy_xxyzzzzz_1[j] + 1.5 * fl1_fx * tg_yy_xxyzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxyzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyy_xxzzzzz_1[j];

                    tg_yyyy_xxzzzzzz_0[j] = pb_y * tg_yyy_xxzzzzzz_0[j] + wp_y[j] * tg_yyy_xxzzzzzz_1[j] + 1.5 * fl1_fx * tg_yy_xxzzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxzzzzzz_1[j];

                    tg_yyyy_xyyyyyyy_0[j] = pb_y * tg_yyy_xyyyyyyy_0[j] + wp_y[j] * tg_yyy_xyyyyyyy_1[j] + 1.5 * fl1_fx * tg_yy_xyyyyyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xyyyyyyy_1[j] + 3.5 * fl1_fxn * tg_yyy_xyyyyyy_1[j];

                    tg_yyyy_xyyyyyyz_0[j] = pb_y * tg_yyy_xyyyyyyz_0[j] + wp_y[j] * tg_yyy_xyyyyyyz_1[j] + 1.5 * fl1_fx * tg_yy_xyyyyyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xyyyyyyz_1[j] + 3.0 * fl1_fxn * tg_yyy_xyyyyyz_1[j];

                    tg_yyyy_xyyyyyzz_0[j] = pb_y * tg_yyy_xyyyyyzz_0[j] + wp_y[j] * tg_yyy_xyyyyyzz_1[j] + 1.5 * fl1_fx * tg_yy_xyyyyyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xyyyyyzz_1[j] + 2.5 * fl1_fxn * tg_yyy_xyyyyzz_1[j];

                    tg_yyyy_xyyyyzzz_0[j] = pb_y * tg_yyy_xyyyyzzz_0[j] + wp_y[j] * tg_yyy_xyyyyzzz_1[j] + 1.5 * fl1_fx * tg_yy_xyyyyzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xyyyyzzz_1[j] + 2.0 * fl1_fxn * tg_yyy_xyyyzzz_1[j];

                    tg_yyyy_xyyyzzzz_0[j] = pb_y * tg_yyy_xyyyzzzz_0[j] + wp_y[j] * tg_yyy_xyyyzzzz_1[j] + 1.5 * fl1_fx * tg_yy_xyyyzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xyyyzzzz_1[j] + 1.5 * fl1_fxn * tg_yyy_xyyzzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSL_483_531(      CMemBlock2D<double>& primBuffer,
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

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {4, -1, -1, -1},
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_2_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_2_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

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

                auto pb_y = r_pb_y[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_y = wpDistances.data(3 * idx + 1);

                // set up pointers to auxilary integrals

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

                auto tg_yyy_xyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 303); 

                auto tg_yyy_xyzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 304); 

                auto tg_yyy_xzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 305); 

                auto tg_yyy_yyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 306); 

                auto tg_yyy_yyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 307); 

                auto tg_yyy_yyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 308); 

                auto tg_yyy_yyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 309); 

                auto tg_yyy_yyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 310); 

                auto tg_yyy_yyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 311); 

                auto tg_yyy_yyzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 312); 

                auto tg_yyy_yzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 313); 

                auto tg_yyy_zzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 314); 

                auto tg_yyz_xxxxxxxx_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 315); 

                auto tg_yyz_xxxxxxxy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 316); 

                auto tg_yyz_xxxxxxxz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 317); 

                auto tg_yyz_xxxxxxyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 318); 

                auto tg_yyz_xxxxxxyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 319); 

                auto tg_yyz_xxxxxxzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 320); 

                auto tg_yyz_xxxxxyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 321); 

                auto tg_yyz_xxxxxyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 322); 

                auto tg_yyz_xxxxxyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 323); 

                auto tg_yyz_xxxxxzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 324); 

                auto tg_yyz_xxxxyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 325); 

                auto tg_yyz_xxxxyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 326); 

                auto tg_yyz_xxxxyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 327); 

                auto tg_yyz_xxxxyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 328); 

                auto tg_yyz_xxxxzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 329); 

                auto tg_yyz_xxxyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 330); 

                auto tg_yyz_xxxyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 331); 

                auto tg_yyz_xxxyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 332); 

                auto tg_yyz_xxxyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 333); 

                auto tg_yyz_xxxyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 334); 

                auto tg_yyz_xxxzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 335); 

                auto tg_yyz_xxyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 336); 

                auto tg_yyz_xxyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 337); 

                auto tg_yyz_xxyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 338); 

                auto tg_yyz_xxyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 339); 

                auto tg_yyz_xxyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 340); 

                auto tg_yyz_xxyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 341); 

                auto tg_yyz_xxzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 342); 

                auto tg_yyz_xyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 343); 

                auto tg_yyz_xyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 344); 

                auto tg_yyz_xyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 345); 

                auto tg_yyz_xyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 346); 

                auto tg_yyz_xyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 347); 

                auto tg_yyz_xyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 348); 

                auto tg_yyz_xyzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 349); 

                auto tg_yyz_xzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 350); 

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

                // set up pointers to integrals

                auto tg_yyyy_xyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 483); 

                auto tg_yyyy_xyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 484); 

                auto tg_yyyy_xzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 485); 

                auto tg_yyyy_yyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 486); 

                auto tg_yyyy_yyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 487); 

                auto tg_yyyy_yyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 488); 

                auto tg_yyyy_yyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 489); 

                auto tg_yyyy_yyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 490); 

                auto tg_yyyy_yyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 491); 

                auto tg_yyyy_yyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 492); 

                auto tg_yyyy_yzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 493); 

                auto tg_yyyy_zzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 494); 

                auto tg_yyyz_xxxxxxxx_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 495); 

                auto tg_yyyz_xxxxxxxy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 496); 

                auto tg_yyyz_xxxxxxxz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 497); 

                auto tg_yyyz_xxxxxxyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 498); 

                auto tg_yyyz_xxxxxxyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 499); 

                auto tg_yyyz_xxxxxxzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 500); 

                auto tg_yyyz_xxxxxyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 501); 

                auto tg_yyyz_xxxxxyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 502); 

                auto tg_yyyz_xxxxxyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 503); 

                auto tg_yyyz_xxxxxzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 504); 

                auto tg_yyyz_xxxxyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 505); 

                auto tg_yyyz_xxxxyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 506); 

                auto tg_yyyz_xxxxyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 507); 

                auto tg_yyyz_xxxxyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 508); 

                auto tg_yyyz_xxxxzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 509); 

                auto tg_yyyz_xxxyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 510); 

                auto tg_yyyz_xxxyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 511); 

                auto tg_yyyz_xxxyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 512); 

                auto tg_yyyz_xxxyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 513); 

                auto tg_yyyz_xxxyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 514); 

                auto tg_yyyz_xxxzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 515); 

                auto tg_yyyz_xxyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 516); 

                auto tg_yyyz_xxyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 517); 

                auto tg_yyyz_xxyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 518); 

                auto tg_yyyz_xxyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 519); 

                auto tg_yyyz_xxyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 520); 

                auto tg_yyyz_xxyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 521); 

                auto tg_yyyz_xxzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 522); 

                auto tg_yyyz_xyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 523); 

                auto tg_yyyz_xyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 524); 

                auto tg_yyyz_xyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 525); 

                auto tg_yyyz_xyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 526); 

                auto tg_yyyz_xyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 527); 

                auto tg_yyyz_xyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 528); 

                auto tg_yyyz_xyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 529); 

                auto tg_yyyz_xzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 530); 

                // Batch of Integrals (483,531)

                #pragma omp simd aligned(fxn, fza, tg_yy_xyyzzzzz_0, tg_yy_xyyzzzzz_1, tg_yy_xyzzzzzz_0, \
                                         tg_yy_xyzzzzzz_1, tg_yy_xzzzzzzz_0, tg_yy_xzzzzzzz_1, tg_yy_yyyyyyyy_0, \
                                         tg_yy_yyyyyyyy_1, tg_yy_yyyyyyyz_0, tg_yy_yyyyyyyz_1, tg_yy_yyyyyyzz_0, \
                                         tg_yy_yyyyyyzz_1, tg_yy_yyyyyzzz_0, tg_yy_yyyyyzzz_1, tg_yy_yyyyzzzz_0, \
                                         tg_yy_yyyyzzzz_1, tg_yy_yyyzzzzz_0, tg_yy_yyyzzzzz_1, tg_yy_yyzzzzzz_0, \
                                         tg_yy_yyzzzzzz_1, tg_yy_yzzzzzzz_0, tg_yy_yzzzzzzz_1, tg_yy_zzzzzzzz_0, \
                                         tg_yy_zzzzzzzz_1, tg_yyy_xyyzzzzz_0, tg_yyy_xyyzzzzz_1, tg_yyy_xyzzzzz_1, \
                                         tg_yyy_xyzzzzzz_0, tg_yyy_xyzzzzzz_1, tg_yyy_xzzzzzz_1, tg_yyy_xzzzzzzz_0, \
                                         tg_yyy_xzzzzzzz_1, tg_yyy_yyyyyyy_1, tg_yyy_yyyyyyyy_0, tg_yyy_yyyyyyyy_1, \
                                         tg_yyy_yyyyyyyz_0, tg_yyy_yyyyyyyz_1, tg_yyy_yyyyyyz_1, tg_yyy_yyyyyyzz_0, \
                                         tg_yyy_yyyyyyzz_1, tg_yyy_yyyyyzz_1, tg_yyy_yyyyyzzz_0, tg_yyy_yyyyyzzz_1, \
                                         tg_yyy_yyyyzzz_1, tg_yyy_yyyyzzzz_0, tg_yyy_yyyyzzzz_1, tg_yyy_yyyzzzz_1, \
                                         tg_yyy_yyyzzzzz_0, tg_yyy_yyyzzzzz_1, tg_yyy_yyzzzzz_1, tg_yyy_yyzzzzzz_0, \
                                         tg_yyy_yyzzzzzz_1, tg_yyy_yzzzzzz_1, tg_yyy_yzzzzzzz_0, tg_yyy_yzzzzzzz_1, \
                                         tg_yyy_zzzzzzz_1, tg_yyy_zzzzzzzz_0, tg_yyy_zzzzzzzz_1, tg_yyyy_xyyzzzzz_0, \
                                         tg_yyyy_xyzzzzzz_0, tg_yyyy_xzzzzzzz_0, tg_yyyy_yyyyyyyy_0, tg_yyyy_yyyyyyyz_0, \
                                         tg_yyyy_yyyyyyzz_0, tg_yyyy_yyyyyzzz_0, tg_yyyy_yyyyzzzz_0, tg_yyyy_yyyzzzzz_0, \
                                         tg_yyyy_yyzzzzzz_0, tg_yyyy_yzzzzzzz_0, tg_yyyy_zzzzzzzz_0, tg_yyyz_xxxxxxxx_0, \
                                         tg_yyyz_xxxxxxxy_0, tg_yyyz_xxxxxxxz_0, tg_yyyz_xxxxxxyy_0, tg_yyyz_xxxxxxyz_0, \
                                         tg_yyyz_xxxxxxzz_0, tg_yyyz_xxxxxyyy_0, tg_yyyz_xxxxxyyz_0, tg_yyyz_xxxxxyzz_0, \
                                         tg_yyyz_xxxxxzzz_0, tg_yyyz_xxxxyyyy_0, tg_yyyz_xxxxyyyz_0, tg_yyyz_xxxxyyzz_0, \
                                         tg_yyyz_xxxxyzzz_0, tg_yyyz_xxxxzzzz_0, tg_yyyz_xxxyyyyy_0, tg_yyyz_xxxyyyyz_0, \
                                         tg_yyyz_xxxyyyzz_0, tg_yyyz_xxxyyzzz_0, tg_yyyz_xxxyzzzz_0, tg_yyyz_xxxzzzzz_0, \
                                         tg_yyyz_xxyyyyyy_0, tg_yyyz_xxyyyyyz_0, tg_yyyz_xxyyyyzz_0, tg_yyyz_xxyyyzzz_0, \
                                         tg_yyyz_xxyyzzzz_0, tg_yyyz_xxyzzzzz_0, tg_yyyz_xxzzzzzz_0, tg_yyyz_xyyyyyyy_0, \
                                         tg_yyyz_xyyyyyyz_0, tg_yyyz_xyyyyyzz_0, tg_yyyz_xyyyyzzz_0, tg_yyyz_xyyyzzzz_0, \
                                         tg_yyyz_xyyzzzzz_0, tg_yyyz_xyzzzzzz_0, tg_yyyz_xzzzzzzz_0, tg_yyz_xxxxxxx_1, \
                                         tg_yyz_xxxxxxxx_0, tg_yyz_xxxxxxxx_1, tg_yyz_xxxxxxxy_0, tg_yyz_xxxxxxxy_1, \
                                         tg_yyz_xxxxxxxz_0, tg_yyz_xxxxxxxz_1, tg_yyz_xxxxxxy_1, tg_yyz_xxxxxxyy_0, \
                                         tg_yyz_xxxxxxyy_1, tg_yyz_xxxxxxyz_0, tg_yyz_xxxxxxyz_1, tg_yyz_xxxxxxz_1, \
                                         tg_yyz_xxxxxxzz_0, tg_yyz_xxxxxxzz_1, tg_yyz_xxxxxyy_1, tg_yyz_xxxxxyyy_0, \
                                         tg_yyz_xxxxxyyy_1, tg_yyz_xxxxxyyz_0, tg_yyz_xxxxxyyz_1, tg_yyz_xxxxxyz_1, \
                                         tg_yyz_xxxxxyzz_0, tg_yyz_xxxxxyzz_1, tg_yyz_xxxxxzz_1, tg_yyz_xxxxxzzz_0, \
                                         tg_yyz_xxxxxzzz_1, tg_yyz_xxxxyyy_1, tg_yyz_xxxxyyyy_0, tg_yyz_xxxxyyyy_1, \
                                         tg_yyz_xxxxyyyz_0, tg_yyz_xxxxyyyz_1, tg_yyz_xxxxyyz_1, tg_yyz_xxxxyyzz_0, \
                                         tg_yyz_xxxxyyzz_1, tg_yyz_xxxxyzz_1, tg_yyz_xxxxyzzz_0, tg_yyz_xxxxyzzz_1, \
                                         tg_yyz_xxxxzzz_1, tg_yyz_xxxxzzzz_0, tg_yyz_xxxxzzzz_1, tg_yyz_xxxyyyy_1, \
                                         tg_yyz_xxxyyyyy_0, tg_yyz_xxxyyyyy_1, tg_yyz_xxxyyyyz_0, tg_yyz_xxxyyyyz_1, \
                                         tg_yyz_xxxyyyz_1, tg_yyz_xxxyyyzz_0, tg_yyz_xxxyyyzz_1, tg_yyz_xxxyyzz_1, \
                                         tg_yyz_xxxyyzzz_0, tg_yyz_xxxyyzzz_1, tg_yyz_xxxyzzz_1, tg_yyz_xxxyzzzz_0, \
                                         tg_yyz_xxxyzzzz_1, tg_yyz_xxxzzzz_1, tg_yyz_xxxzzzzz_0, tg_yyz_xxxzzzzz_1, \
                                         tg_yyz_xxyyyyy_1, tg_yyz_xxyyyyyy_0, tg_yyz_xxyyyyyy_1, tg_yyz_xxyyyyyz_0, \
                                         tg_yyz_xxyyyyyz_1, tg_yyz_xxyyyyz_1, tg_yyz_xxyyyyzz_0, tg_yyz_xxyyyyzz_1, \
                                         tg_yyz_xxyyyzz_1, tg_yyz_xxyyyzzz_0, tg_yyz_xxyyyzzz_1, tg_yyz_xxyyzzz_1, \
                                         tg_yyz_xxyyzzzz_0, tg_yyz_xxyyzzzz_1, tg_yyz_xxyzzzz_1, tg_yyz_xxyzzzzz_0, \
                                         tg_yyz_xxyzzzzz_1, tg_yyz_xxzzzzz_1, tg_yyz_xxzzzzzz_0, tg_yyz_xxzzzzzz_1, \
                                         tg_yyz_xyyyyyy_1, tg_yyz_xyyyyyyy_0, tg_yyz_xyyyyyyy_1, tg_yyz_xyyyyyyz_0, \
                                         tg_yyz_xyyyyyyz_1, tg_yyz_xyyyyyz_1, tg_yyz_xyyyyyzz_0, tg_yyz_xyyyyyzz_1, \
                                         tg_yyz_xyyyyzz_1, tg_yyz_xyyyyzzz_0, tg_yyz_xyyyyzzz_1, tg_yyz_xyyyzzz_1, \
                                         tg_yyz_xyyyzzzz_0, tg_yyz_xyyyzzzz_1, tg_yyz_xyyzzzz_1, tg_yyz_xyyzzzzz_0, \
                                         tg_yyz_xyyzzzzz_1, tg_yyz_xyzzzzz_1, tg_yyz_xyzzzzzz_0, tg_yyz_xyzzzzzz_1, \
                                         tg_yyz_xzzzzzz_1, tg_yyz_xzzzzzzz_0, tg_yyz_xzzzzzzz_1, tg_yz_xxxxxxxx_0, \
                                         tg_yz_xxxxxxxx_1, tg_yz_xxxxxxxy_0, tg_yz_xxxxxxxy_1, tg_yz_xxxxxxxz_0, \
                                         tg_yz_xxxxxxxz_1, tg_yz_xxxxxxyy_0, tg_yz_xxxxxxyy_1, tg_yz_xxxxxxyz_0, \
                                         tg_yz_xxxxxxyz_1, tg_yz_xxxxxxzz_0, tg_yz_xxxxxxzz_1, tg_yz_xxxxxyyy_0, \
                                         tg_yz_xxxxxyyy_1, tg_yz_xxxxxyyz_0, tg_yz_xxxxxyyz_1, tg_yz_xxxxxyzz_0, \
                                         tg_yz_xxxxxyzz_1, tg_yz_xxxxxzzz_0, tg_yz_xxxxxzzz_1, tg_yz_xxxxyyyy_0, \
                                         tg_yz_xxxxyyyy_1, tg_yz_xxxxyyyz_0, tg_yz_xxxxyyyz_1, tg_yz_xxxxyyzz_0, \
                                         tg_yz_xxxxyyzz_1, tg_yz_xxxxyzzz_0, tg_yz_xxxxyzzz_1, tg_yz_xxxxzzzz_0, \
                                         tg_yz_xxxxzzzz_1, tg_yz_xxxyyyyy_0, tg_yz_xxxyyyyy_1, tg_yz_xxxyyyyz_0, \
                                         tg_yz_xxxyyyyz_1, tg_yz_xxxyyyzz_0, tg_yz_xxxyyyzz_1, tg_yz_xxxyyzzz_0, \
                                         tg_yz_xxxyyzzz_1, tg_yz_xxxyzzzz_0, tg_yz_xxxyzzzz_1, tg_yz_xxxzzzzz_0, \
                                         tg_yz_xxxzzzzz_1, tg_yz_xxyyyyyy_0, tg_yz_xxyyyyyy_1, tg_yz_xxyyyyyz_0, \
                                         tg_yz_xxyyyyyz_1, tg_yz_xxyyyyzz_0, tg_yz_xxyyyyzz_1, tg_yz_xxyyyzzz_0, \
                                         tg_yz_xxyyyzzz_1, tg_yz_xxyyzzzz_0, tg_yz_xxyyzzzz_1, tg_yz_xxyzzzzz_0, \
                                         tg_yz_xxyzzzzz_1, tg_yz_xxzzzzzz_0, tg_yz_xxzzzzzz_1, tg_yz_xyyyyyyy_0, \
                                         tg_yz_xyyyyyyy_1, tg_yz_xyyyyyyz_0, tg_yz_xyyyyyyz_1, tg_yz_xyyyyyzz_0, \
                                         tg_yz_xyyyyyzz_1, tg_yz_xyyyyzzz_0, tg_yz_xyyyyzzz_1, tg_yz_xyyyzzzz_0, \
                                         tg_yz_xyyyzzzz_1, tg_yz_xyyzzzzz_0, tg_yz_xyyzzzzz_1, tg_yz_xyzzzzzz_0, \
                                         tg_yz_xyzzzzzz_1, tg_yz_xzzzzzzz_0, tg_yz_xzzzzzzz_1, wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_yyyy_xyyzzzzz_0[j] = pb_y * tg_yyy_xyyzzzzz_0[j] + wp_y[j] * tg_yyy_xyyzzzzz_1[j] + 1.5 * fl1_fx * tg_yy_xyyzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xyyzzzzz_1[j] + fl1_fxn * tg_yyy_xyzzzzz_1[j];

                    tg_yyyy_xyzzzzzz_0[j] = pb_y * tg_yyy_xyzzzzzz_0[j] + wp_y[j] * tg_yyy_xyzzzzzz_1[j] + 1.5 * fl1_fx * tg_yy_xyzzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyy_xzzzzzz_1[j];

                    tg_yyyy_xzzzzzzz_0[j] = pb_y * tg_yyy_xzzzzzzz_0[j] + wp_y[j] * tg_yyy_xzzzzzzz_1[j] + 1.5 * fl1_fx * tg_yy_xzzzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xzzzzzzz_1[j];

                    tg_yyyy_yyyyyyyy_0[j] = pb_y * tg_yyy_yyyyyyyy_0[j] + wp_y[j] * tg_yyy_yyyyyyyy_1[j] + 1.5 * fl1_fx * tg_yy_yyyyyyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_yyyyyyyy_1[j] + 4.0 * fl1_fxn * tg_yyy_yyyyyyy_1[j];

                    tg_yyyy_yyyyyyyz_0[j] = pb_y * tg_yyy_yyyyyyyz_0[j] + wp_y[j] * tg_yyy_yyyyyyyz_1[j] + 1.5 * fl1_fx * tg_yy_yyyyyyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_yyyyyyyz_1[j] + 3.5 * fl1_fxn * tg_yyy_yyyyyyz_1[j];

                    tg_yyyy_yyyyyyzz_0[j] = pb_y * tg_yyy_yyyyyyzz_0[j] + wp_y[j] * tg_yyy_yyyyyyzz_1[j] + 1.5 * fl1_fx * tg_yy_yyyyyyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_yyyyyyzz_1[j] + 3.0 * fl1_fxn * tg_yyy_yyyyyzz_1[j];

                    tg_yyyy_yyyyyzzz_0[j] = pb_y * tg_yyy_yyyyyzzz_0[j] + wp_y[j] * tg_yyy_yyyyyzzz_1[j] + 1.5 * fl1_fx * tg_yy_yyyyyzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_yyyyyzzz_1[j] + 2.5 * fl1_fxn * tg_yyy_yyyyzzz_1[j];

                    tg_yyyy_yyyyzzzz_0[j] = pb_y * tg_yyy_yyyyzzzz_0[j] + wp_y[j] * tg_yyy_yyyyzzzz_1[j] + 1.5 * fl1_fx * tg_yy_yyyyzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_yyyyzzzz_1[j] + 2.0 * fl1_fxn * tg_yyy_yyyzzzz_1[j];

                    tg_yyyy_yyyzzzzz_0[j] = pb_y * tg_yyy_yyyzzzzz_0[j] + wp_y[j] * tg_yyy_yyyzzzzz_1[j] + 1.5 * fl1_fx * tg_yy_yyyzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_yyyzzzzz_1[j] + 1.5 * fl1_fxn * tg_yyy_yyzzzzz_1[j];

                    tg_yyyy_yyzzzzzz_0[j] = pb_y * tg_yyy_yyzzzzzz_0[j] + wp_y[j] * tg_yyy_yyzzzzzz_1[j] + 1.5 * fl1_fx * tg_yy_yyzzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_yyzzzzzz_1[j] + fl1_fxn * tg_yyy_yzzzzzz_1[j];

                    tg_yyyy_yzzzzzzz_0[j] = pb_y * tg_yyy_yzzzzzzz_0[j] + wp_y[j] * tg_yyy_yzzzzzzz_1[j] + 1.5 * fl1_fx * tg_yy_yzzzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_yzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyy_zzzzzzz_1[j];

                    tg_yyyy_zzzzzzzz_0[j] = pb_y * tg_yyy_zzzzzzzz_0[j] + wp_y[j] * tg_yyy_zzzzzzzz_1[j] + 1.5 * fl1_fx * tg_yy_zzzzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_zzzzzzzz_1[j];

                    tg_yyyz_xxxxxxxx_0[j] = pb_y * tg_yyz_xxxxxxxx_0[j] + wp_y[j] * tg_yyz_xxxxxxxx_1[j] + fl1_fx * tg_yz_xxxxxxxx_0[j] - fl1_fx * fl1_fza * tg_yz_xxxxxxxx_1[j];

                    tg_yyyz_xxxxxxxy_0[j] = pb_y * tg_yyz_xxxxxxxy_0[j] + wp_y[j] * tg_yyz_xxxxxxxy_1[j] + fl1_fx * tg_yz_xxxxxxxy_0[j] - fl1_fx * fl1_fza * tg_yz_xxxxxxxy_1[j] + 0.5 * fl1_fxn * tg_yyz_xxxxxxx_1[j];

                    tg_yyyz_xxxxxxxz_0[j] = pb_y * tg_yyz_xxxxxxxz_0[j] + wp_y[j] * tg_yyz_xxxxxxxz_1[j] + fl1_fx * tg_yz_xxxxxxxz_0[j] - fl1_fx * fl1_fza * tg_yz_xxxxxxxz_1[j];

                    tg_yyyz_xxxxxxyy_0[j] = pb_y * tg_yyz_xxxxxxyy_0[j] + wp_y[j] * tg_yyz_xxxxxxyy_1[j] + fl1_fx * tg_yz_xxxxxxyy_0[j] - fl1_fx * fl1_fza * tg_yz_xxxxxxyy_1[j] + fl1_fxn * tg_yyz_xxxxxxy_1[j];

                    tg_yyyz_xxxxxxyz_0[j] = pb_y * tg_yyz_xxxxxxyz_0[j] + wp_y[j] * tg_yyz_xxxxxxyz_1[j] + fl1_fx * tg_yz_xxxxxxyz_0[j] - fl1_fx * fl1_fza * tg_yz_xxxxxxyz_1[j] + 0.5 * fl1_fxn * tg_yyz_xxxxxxz_1[j];

                    tg_yyyz_xxxxxxzz_0[j] = pb_y * tg_yyz_xxxxxxzz_0[j] + wp_y[j] * tg_yyz_xxxxxxzz_1[j] + fl1_fx * tg_yz_xxxxxxzz_0[j] - fl1_fx * fl1_fza * tg_yz_xxxxxxzz_1[j];

                    tg_yyyz_xxxxxyyy_0[j] = pb_y * tg_yyz_xxxxxyyy_0[j] + wp_y[j] * tg_yyz_xxxxxyyy_1[j] + fl1_fx * tg_yz_xxxxxyyy_0[j] - fl1_fx * fl1_fza * tg_yz_xxxxxyyy_1[j] + 1.5 * fl1_fxn * tg_yyz_xxxxxyy_1[j];

                    tg_yyyz_xxxxxyyz_0[j] = pb_y * tg_yyz_xxxxxyyz_0[j] + wp_y[j] * tg_yyz_xxxxxyyz_1[j] + fl1_fx * tg_yz_xxxxxyyz_0[j] - fl1_fx * fl1_fza * tg_yz_xxxxxyyz_1[j] + fl1_fxn * tg_yyz_xxxxxyz_1[j];

                    tg_yyyz_xxxxxyzz_0[j] = pb_y * tg_yyz_xxxxxyzz_0[j] + wp_y[j] * tg_yyz_xxxxxyzz_1[j] + fl1_fx * tg_yz_xxxxxyzz_0[j] - fl1_fx * fl1_fza * tg_yz_xxxxxyzz_1[j] + 0.5 * fl1_fxn * tg_yyz_xxxxxzz_1[j];

                    tg_yyyz_xxxxxzzz_0[j] = pb_y * tg_yyz_xxxxxzzz_0[j] + wp_y[j] * tg_yyz_xxxxxzzz_1[j] + fl1_fx * tg_yz_xxxxxzzz_0[j] - fl1_fx * fl1_fza * tg_yz_xxxxxzzz_1[j];

                    tg_yyyz_xxxxyyyy_0[j] = pb_y * tg_yyz_xxxxyyyy_0[j] + wp_y[j] * tg_yyz_xxxxyyyy_1[j] + fl1_fx * tg_yz_xxxxyyyy_0[j] - fl1_fx * fl1_fza * tg_yz_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_yyz_xxxxyyy_1[j];

                    tg_yyyz_xxxxyyyz_0[j] = pb_y * tg_yyz_xxxxyyyz_0[j] + wp_y[j] * tg_yyz_xxxxyyyz_1[j] + fl1_fx * tg_yz_xxxxyyyz_0[j] - fl1_fx * fl1_fza * tg_yz_xxxxyyyz_1[j] + 1.5 * fl1_fxn * tg_yyz_xxxxyyz_1[j];

                    tg_yyyz_xxxxyyzz_0[j] = pb_y * tg_yyz_xxxxyyzz_0[j] + wp_y[j] * tg_yyz_xxxxyyzz_1[j] + fl1_fx * tg_yz_xxxxyyzz_0[j] - fl1_fx * fl1_fza * tg_yz_xxxxyyzz_1[j] + fl1_fxn * tg_yyz_xxxxyzz_1[j];

                    tg_yyyz_xxxxyzzz_0[j] = pb_y * tg_yyz_xxxxyzzz_0[j] + wp_y[j] * tg_yyz_xxxxyzzz_1[j] + fl1_fx * tg_yz_xxxxyzzz_0[j] - fl1_fx * fl1_fza * tg_yz_xxxxyzzz_1[j] + 0.5 * fl1_fxn * tg_yyz_xxxxzzz_1[j];

                    tg_yyyz_xxxxzzzz_0[j] = pb_y * tg_yyz_xxxxzzzz_0[j] + wp_y[j] * tg_yyz_xxxxzzzz_1[j] + fl1_fx * tg_yz_xxxxzzzz_0[j] - fl1_fx * fl1_fza * tg_yz_xxxxzzzz_1[j];

                    tg_yyyz_xxxyyyyy_0[j] = pb_y * tg_yyz_xxxyyyyy_0[j] + wp_y[j] * tg_yyz_xxxyyyyy_1[j] + fl1_fx * tg_yz_xxxyyyyy_0[j] - fl1_fx * fl1_fza * tg_yz_xxxyyyyy_1[j] + 2.5 * fl1_fxn * tg_yyz_xxxyyyy_1[j];

                    tg_yyyz_xxxyyyyz_0[j] = pb_y * tg_yyz_xxxyyyyz_0[j] + wp_y[j] * tg_yyz_xxxyyyyz_1[j] + fl1_fx * tg_yz_xxxyyyyz_0[j] - fl1_fx * fl1_fza * tg_yz_xxxyyyyz_1[j] + 2.0 * fl1_fxn * tg_yyz_xxxyyyz_1[j];

                    tg_yyyz_xxxyyyzz_0[j] = pb_y * tg_yyz_xxxyyyzz_0[j] + wp_y[j] * tg_yyz_xxxyyyzz_1[j] + fl1_fx * tg_yz_xxxyyyzz_0[j] - fl1_fx * fl1_fza * tg_yz_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_yyz_xxxyyzz_1[j];

                    tg_yyyz_xxxyyzzz_0[j] = pb_y * tg_yyz_xxxyyzzz_0[j] + wp_y[j] * tg_yyz_xxxyyzzz_1[j] + fl1_fx * tg_yz_xxxyyzzz_0[j] - fl1_fx * fl1_fza * tg_yz_xxxyyzzz_1[j] + fl1_fxn * tg_yyz_xxxyzzz_1[j];

                    tg_yyyz_xxxyzzzz_0[j] = pb_y * tg_yyz_xxxyzzzz_0[j] + wp_y[j] * tg_yyz_xxxyzzzz_1[j] + fl1_fx * tg_yz_xxxyzzzz_0[j] - fl1_fx * fl1_fza * tg_yz_xxxyzzzz_1[j] + 0.5 * fl1_fxn * tg_yyz_xxxzzzz_1[j];

                    tg_yyyz_xxxzzzzz_0[j] = pb_y * tg_yyz_xxxzzzzz_0[j] + wp_y[j] * tg_yyz_xxxzzzzz_1[j] + fl1_fx * tg_yz_xxxzzzzz_0[j] - fl1_fx * fl1_fza * tg_yz_xxxzzzzz_1[j];

                    tg_yyyz_xxyyyyyy_0[j] = pb_y * tg_yyz_xxyyyyyy_0[j] + wp_y[j] * tg_yyz_xxyyyyyy_1[j] + fl1_fx * tg_yz_xxyyyyyy_0[j] - fl1_fx * fl1_fza * tg_yz_xxyyyyyy_1[j] + 3.0 * fl1_fxn * tg_yyz_xxyyyyy_1[j];

                    tg_yyyz_xxyyyyyz_0[j] = pb_y * tg_yyz_xxyyyyyz_0[j] + wp_y[j] * tg_yyz_xxyyyyyz_1[j] + fl1_fx * tg_yz_xxyyyyyz_0[j] - fl1_fx * fl1_fza * tg_yz_xxyyyyyz_1[j] + 2.5 * fl1_fxn * tg_yyz_xxyyyyz_1[j];

                    tg_yyyz_xxyyyyzz_0[j] = pb_y * tg_yyz_xxyyyyzz_0[j] + wp_y[j] * tg_yyz_xxyyyyzz_1[j] + fl1_fx * tg_yz_xxyyyyzz_0[j] - fl1_fx * fl1_fza * tg_yz_xxyyyyzz_1[j] + 2.0 * fl1_fxn * tg_yyz_xxyyyzz_1[j];

                    tg_yyyz_xxyyyzzz_0[j] = pb_y * tg_yyz_xxyyyzzz_0[j] + wp_y[j] * tg_yyz_xxyyyzzz_1[j] + fl1_fx * tg_yz_xxyyyzzz_0[j] - fl1_fx * fl1_fza * tg_yz_xxyyyzzz_1[j] + 1.5 * fl1_fxn * tg_yyz_xxyyzzz_1[j];

                    tg_yyyz_xxyyzzzz_0[j] = pb_y * tg_yyz_xxyyzzzz_0[j] + wp_y[j] * tg_yyz_xxyyzzzz_1[j] + fl1_fx * tg_yz_xxyyzzzz_0[j] - fl1_fx * fl1_fza * tg_yz_xxyyzzzz_1[j] + fl1_fxn * tg_yyz_xxyzzzz_1[j];

                    tg_yyyz_xxyzzzzz_0[j] = pb_y * tg_yyz_xxyzzzzz_0[j] + wp_y[j] * tg_yyz_xxyzzzzz_1[j] + fl1_fx * tg_yz_xxyzzzzz_0[j] - fl1_fx * fl1_fza * tg_yz_xxyzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyz_xxzzzzz_1[j];

                    tg_yyyz_xxzzzzzz_0[j] = pb_y * tg_yyz_xxzzzzzz_0[j] + wp_y[j] * tg_yyz_xxzzzzzz_1[j] + fl1_fx * tg_yz_xxzzzzzz_0[j] - fl1_fx * fl1_fza * tg_yz_xxzzzzzz_1[j];

                    tg_yyyz_xyyyyyyy_0[j] = pb_y * tg_yyz_xyyyyyyy_0[j] + wp_y[j] * tg_yyz_xyyyyyyy_1[j] + fl1_fx * tg_yz_xyyyyyyy_0[j] - fl1_fx * fl1_fza * tg_yz_xyyyyyyy_1[j] + 3.5 * fl1_fxn * tg_yyz_xyyyyyy_1[j];

                    tg_yyyz_xyyyyyyz_0[j] = pb_y * tg_yyz_xyyyyyyz_0[j] + wp_y[j] * tg_yyz_xyyyyyyz_1[j] + fl1_fx * tg_yz_xyyyyyyz_0[j] - fl1_fx * fl1_fza * tg_yz_xyyyyyyz_1[j] + 3.0 * fl1_fxn * tg_yyz_xyyyyyz_1[j];

                    tg_yyyz_xyyyyyzz_0[j] = pb_y * tg_yyz_xyyyyyzz_0[j] + wp_y[j] * tg_yyz_xyyyyyzz_1[j] + fl1_fx * tg_yz_xyyyyyzz_0[j] - fl1_fx * fl1_fza * tg_yz_xyyyyyzz_1[j] + 2.5 * fl1_fxn * tg_yyz_xyyyyzz_1[j];

                    tg_yyyz_xyyyyzzz_0[j] = pb_y * tg_yyz_xyyyyzzz_0[j] + wp_y[j] * tg_yyz_xyyyyzzz_1[j] + fl1_fx * tg_yz_xyyyyzzz_0[j] - fl1_fx * fl1_fza * tg_yz_xyyyyzzz_1[j] + 2.0 * fl1_fxn * tg_yyz_xyyyzzz_1[j];

                    tg_yyyz_xyyyzzzz_0[j] = pb_y * tg_yyz_xyyyzzzz_0[j] + wp_y[j] * tg_yyz_xyyyzzzz_1[j] + fl1_fx * tg_yz_xyyyzzzz_0[j] - fl1_fx * fl1_fza * tg_yz_xyyyzzzz_1[j] + 1.5 * fl1_fxn * tg_yyz_xyyzzzz_1[j];

                    tg_yyyz_xyyzzzzz_0[j] = pb_y * tg_yyz_xyyzzzzz_0[j] + wp_y[j] * tg_yyz_xyyzzzzz_1[j] + fl1_fx * tg_yz_xyyzzzzz_0[j] - fl1_fx * fl1_fza * tg_yz_xyyzzzzz_1[j] + fl1_fxn * tg_yyz_xyzzzzz_1[j];

                    tg_yyyz_xyzzzzzz_0[j] = pb_y * tg_yyz_xyzzzzzz_0[j] + wp_y[j] * tg_yyz_xyzzzzzz_1[j] + fl1_fx * tg_yz_xyzzzzzz_0[j] - fl1_fx * fl1_fza * tg_yz_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyz_xzzzzzz_1[j];

                    tg_yyyz_xzzzzzzz_0[j] = pb_y * tg_yyz_xzzzzzzz_0[j] + wp_y[j] * tg_yyz_xzzzzzzz_1[j] + fl1_fx * tg_yz_xzzzzzzz_0[j] - fl1_fx * fl1_fza * tg_yz_xzzzzzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSL_531_579(      CMemBlock2D<double>& primBuffer,
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

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {4, -1, -1, -1},
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_2_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_2_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

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

                auto pb_y = r_pb_y[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_y = wpDistances.data(3 * idx + 1);

                // set up pointers to auxilary integrals

                auto tg_yyz_yyyyyyyy_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 351); 

                auto tg_yyz_yyyyyyyz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 352); 

                auto tg_yyz_yyyyyyzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 353); 

                auto tg_yyz_yyyyyzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 354); 

                auto tg_yyz_yyyyzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 355); 

                auto tg_yyz_yyyzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 356); 

                auto tg_yyz_yyzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 357); 

                auto tg_yyz_yzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 358); 

                auto tg_yyz_zzzzzzzz_0 = primBuffer.data(pidx_g_3_8_m0 + 450 * idx + 359); 

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

                auto tg_yyz_yyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 351); 

                auto tg_yyz_yyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 352); 

                auto tg_yyz_yyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 353); 

                auto tg_yyz_yyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 354); 

                auto tg_yyz_yyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 355); 

                auto tg_yyz_yyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 356); 

                auto tg_yyz_yyzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 357); 

                auto tg_yyz_yzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 358); 

                auto tg_yyz_zzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 359); 

                auto tg_yzz_xxxxxxxx_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 360); 

                auto tg_yzz_xxxxxxxy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 361); 

                auto tg_yzz_xxxxxxxz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 362); 

                auto tg_yzz_xxxxxxyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 363); 

                auto tg_yzz_xxxxxxyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 364); 

                auto tg_yzz_xxxxxxzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 365); 

                auto tg_yzz_xxxxxyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 366); 

                auto tg_yzz_xxxxxyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 367); 

                auto tg_yzz_xxxxxyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 368); 

                auto tg_yzz_xxxxxzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 369); 

                auto tg_yzz_xxxxyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 370); 

                auto tg_yzz_xxxxyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 371); 

                auto tg_yzz_xxxxyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 372); 

                auto tg_yzz_xxxxyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 373); 

                auto tg_yzz_xxxxzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 374); 

                auto tg_yzz_xxxyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 375); 

                auto tg_yzz_xxxyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 376); 

                auto tg_yzz_xxxyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 377); 

                auto tg_yzz_xxxyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 378); 

                auto tg_yzz_xxxyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 379); 

                auto tg_yzz_xxxzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 380); 

                auto tg_yzz_xxyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 381); 

                auto tg_yzz_xxyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 382); 

                auto tg_yzz_xxyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 383); 

                auto tg_yzz_xxyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 384); 

                auto tg_yzz_xxyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 385); 

                auto tg_yzz_xxyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 386); 

                auto tg_yzz_xxzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 387); 

                auto tg_yzz_xyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 388); 

                auto tg_yzz_xyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 389); 

                auto tg_yzz_xyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 390); 

                auto tg_yzz_xyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 391); 

                auto tg_yzz_xyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 392); 

                auto tg_yzz_xyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 393); 

                auto tg_yzz_xyzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 394); 

                auto tg_yzz_xzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 395); 

                auto tg_yzz_yyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 396); 

                auto tg_yzz_yyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 397); 

                auto tg_yzz_yyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 398); 

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

                // set up pointers to integrals

                auto tg_yyyz_yyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 531); 

                auto tg_yyyz_yyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 532); 

                auto tg_yyyz_yyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 533); 

                auto tg_yyyz_yyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 534); 

                auto tg_yyyz_yyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 535); 

                auto tg_yyyz_yyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 536); 

                auto tg_yyyz_yyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 537); 

                auto tg_yyyz_yzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 538); 

                auto tg_yyyz_zzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 539); 

                auto tg_yyzz_xxxxxxxx_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 540); 

                auto tg_yyzz_xxxxxxxy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 541); 

                auto tg_yyzz_xxxxxxxz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 542); 

                auto tg_yyzz_xxxxxxyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 543); 

                auto tg_yyzz_xxxxxxyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 544); 

                auto tg_yyzz_xxxxxxzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 545); 

                auto tg_yyzz_xxxxxyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 546); 

                auto tg_yyzz_xxxxxyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 547); 

                auto tg_yyzz_xxxxxyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 548); 

                auto tg_yyzz_xxxxxzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 549); 

                auto tg_yyzz_xxxxyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 550); 

                auto tg_yyzz_xxxxyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 551); 

                auto tg_yyzz_xxxxyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 552); 

                auto tg_yyzz_xxxxyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 553); 

                auto tg_yyzz_xxxxzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 554); 

                auto tg_yyzz_xxxyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 555); 

                auto tg_yyzz_xxxyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 556); 

                auto tg_yyzz_xxxyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 557); 

                auto tg_yyzz_xxxyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 558); 

                auto tg_yyzz_xxxyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 559); 

                auto tg_yyzz_xxxzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 560); 

                auto tg_yyzz_xxyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 561); 

                auto tg_yyzz_xxyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 562); 

                auto tg_yyzz_xxyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 563); 

                auto tg_yyzz_xxyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 564); 

                auto tg_yyzz_xxyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 565); 

                auto tg_yyzz_xxyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 566); 

                auto tg_yyzz_xxzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 567); 

                auto tg_yyzz_xyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 568); 

                auto tg_yyzz_xyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 569); 

                auto tg_yyzz_xyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 570); 

                auto tg_yyzz_xyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 571); 

                auto tg_yyzz_xyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 572); 

                auto tg_yyzz_xyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 573); 

                auto tg_yyzz_xyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 574); 

                auto tg_yyzz_xzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 575); 

                auto tg_yyzz_yyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 576); 

                auto tg_yyzz_yyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 577); 

                auto tg_yyzz_yyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 578); 

                // Batch of Integrals (531,579)

                #pragma omp simd aligned(fxn, fza, tg_yyyz_yyyyyyyy_0, tg_yyyz_yyyyyyyz_0, tg_yyyz_yyyyyyzz_0, \
                                         tg_yyyz_yyyyyzzz_0, tg_yyyz_yyyyzzzz_0, tg_yyyz_yyyzzzzz_0, tg_yyyz_yyzzzzzz_0, \
                                         tg_yyyz_yzzzzzzz_0, tg_yyyz_zzzzzzzz_0, tg_yyz_yyyyyyy_1, tg_yyz_yyyyyyyy_0, \
                                         tg_yyz_yyyyyyyy_1, tg_yyz_yyyyyyyz_0, tg_yyz_yyyyyyyz_1, tg_yyz_yyyyyyz_1, \
                                         tg_yyz_yyyyyyzz_0, tg_yyz_yyyyyyzz_1, tg_yyz_yyyyyzz_1, tg_yyz_yyyyyzzz_0, \
                                         tg_yyz_yyyyyzzz_1, tg_yyz_yyyyzzz_1, tg_yyz_yyyyzzzz_0, tg_yyz_yyyyzzzz_1, \
                                         tg_yyz_yyyzzzz_1, tg_yyz_yyyzzzzz_0, tg_yyz_yyyzzzzz_1, tg_yyz_yyzzzzz_1, \
                                         tg_yyz_yyzzzzzz_0, tg_yyz_yyzzzzzz_1, tg_yyz_yzzzzzz_1, tg_yyz_yzzzzzzz_0, \
                                         tg_yyz_yzzzzzzz_1, tg_yyz_zzzzzzz_1, tg_yyz_zzzzzzzz_0, tg_yyz_zzzzzzzz_1, \
                                         tg_yyzz_xxxxxxxx_0, tg_yyzz_xxxxxxxy_0, tg_yyzz_xxxxxxxz_0, tg_yyzz_xxxxxxyy_0, \
                                         tg_yyzz_xxxxxxyz_0, tg_yyzz_xxxxxxzz_0, tg_yyzz_xxxxxyyy_0, tg_yyzz_xxxxxyyz_0, \
                                         tg_yyzz_xxxxxyzz_0, tg_yyzz_xxxxxzzz_0, tg_yyzz_xxxxyyyy_0, tg_yyzz_xxxxyyyz_0, \
                                         tg_yyzz_xxxxyyzz_0, tg_yyzz_xxxxyzzz_0, tg_yyzz_xxxxzzzz_0, tg_yyzz_xxxyyyyy_0, \
                                         tg_yyzz_xxxyyyyz_0, tg_yyzz_xxxyyyzz_0, tg_yyzz_xxxyyzzz_0, tg_yyzz_xxxyzzzz_0, \
                                         tg_yyzz_xxxzzzzz_0, tg_yyzz_xxyyyyyy_0, tg_yyzz_xxyyyyyz_0, tg_yyzz_xxyyyyzz_0, \
                                         tg_yyzz_xxyyyzzz_0, tg_yyzz_xxyyzzzz_0, tg_yyzz_xxyzzzzz_0, tg_yyzz_xxzzzzzz_0, \
                                         tg_yyzz_xyyyyyyy_0, tg_yyzz_xyyyyyyz_0, tg_yyzz_xyyyyyzz_0, tg_yyzz_xyyyyzzz_0, \
                                         tg_yyzz_xyyyzzzz_0, tg_yyzz_xyyzzzzz_0, tg_yyzz_xyzzzzzz_0, tg_yyzz_xzzzzzzz_0, \
                                         tg_yyzz_yyyyyyyy_0, tg_yyzz_yyyyyyyz_0, tg_yyzz_yyyyyyzz_0, tg_yz_yyyyyyyy_0, \
                                         tg_yz_yyyyyyyy_1, tg_yz_yyyyyyyz_0, tg_yz_yyyyyyyz_1, tg_yz_yyyyyyzz_0, \
                                         tg_yz_yyyyyyzz_1, tg_yz_yyyyyzzz_0, tg_yz_yyyyyzzz_1, tg_yz_yyyyzzzz_0, \
                                         tg_yz_yyyyzzzz_1, tg_yz_yyyzzzzz_0, tg_yz_yyyzzzzz_1, tg_yz_yyzzzzzz_0, \
                                         tg_yz_yyzzzzzz_1, tg_yz_yzzzzzzz_0, tg_yz_yzzzzzzz_1, tg_yz_zzzzzzzz_0, \
                                         tg_yz_zzzzzzzz_1, tg_yzz_xxxxxxx_1, tg_yzz_xxxxxxxx_0, tg_yzz_xxxxxxxx_1, \
                                         tg_yzz_xxxxxxxy_0, tg_yzz_xxxxxxxy_1, tg_yzz_xxxxxxxz_0, tg_yzz_xxxxxxxz_1, \
                                         tg_yzz_xxxxxxy_1, tg_yzz_xxxxxxyy_0, tg_yzz_xxxxxxyy_1, tg_yzz_xxxxxxyz_0, \
                                         tg_yzz_xxxxxxyz_1, tg_yzz_xxxxxxz_1, tg_yzz_xxxxxxzz_0, tg_yzz_xxxxxxzz_1, \
                                         tg_yzz_xxxxxyy_1, tg_yzz_xxxxxyyy_0, tg_yzz_xxxxxyyy_1, tg_yzz_xxxxxyyz_0, \
                                         tg_yzz_xxxxxyyz_1, tg_yzz_xxxxxyz_1, tg_yzz_xxxxxyzz_0, tg_yzz_xxxxxyzz_1, \
                                         tg_yzz_xxxxxzz_1, tg_yzz_xxxxxzzz_0, tg_yzz_xxxxxzzz_1, tg_yzz_xxxxyyy_1, \
                                         tg_yzz_xxxxyyyy_0, tg_yzz_xxxxyyyy_1, tg_yzz_xxxxyyyz_0, tg_yzz_xxxxyyyz_1, \
                                         tg_yzz_xxxxyyz_1, tg_yzz_xxxxyyzz_0, tg_yzz_xxxxyyzz_1, tg_yzz_xxxxyzz_1, \
                                         tg_yzz_xxxxyzzz_0, tg_yzz_xxxxyzzz_1, tg_yzz_xxxxzzz_1, tg_yzz_xxxxzzzz_0, \
                                         tg_yzz_xxxxzzzz_1, tg_yzz_xxxyyyy_1, tg_yzz_xxxyyyyy_0, tg_yzz_xxxyyyyy_1, \
                                         tg_yzz_xxxyyyyz_0, tg_yzz_xxxyyyyz_1, tg_yzz_xxxyyyz_1, tg_yzz_xxxyyyzz_0, \
                                         tg_yzz_xxxyyyzz_1, tg_yzz_xxxyyzz_1, tg_yzz_xxxyyzzz_0, tg_yzz_xxxyyzzz_1, \
                                         tg_yzz_xxxyzzz_1, tg_yzz_xxxyzzzz_0, tg_yzz_xxxyzzzz_1, tg_yzz_xxxzzzz_1, \
                                         tg_yzz_xxxzzzzz_0, tg_yzz_xxxzzzzz_1, tg_yzz_xxyyyyy_1, tg_yzz_xxyyyyyy_0, \
                                         tg_yzz_xxyyyyyy_1, tg_yzz_xxyyyyyz_0, tg_yzz_xxyyyyyz_1, tg_yzz_xxyyyyz_1, \
                                         tg_yzz_xxyyyyzz_0, tg_yzz_xxyyyyzz_1, tg_yzz_xxyyyzz_1, tg_yzz_xxyyyzzz_0, \
                                         tg_yzz_xxyyyzzz_1, tg_yzz_xxyyzzz_1, tg_yzz_xxyyzzzz_0, tg_yzz_xxyyzzzz_1, \
                                         tg_yzz_xxyzzzz_1, tg_yzz_xxyzzzzz_0, tg_yzz_xxyzzzzz_1, tg_yzz_xxzzzzz_1, \
                                         tg_yzz_xxzzzzzz_0, tg_yzz_xxzzzzzz_1, tg_yzz_xyyyyyy_1, tg_yzz_xyyyyyyy_0, \
                                         tg_yzz_xyyyyyyy_1, tg_yzz_xyyyyyyz_0, tg_yzz_xyyyyyyz_1, tg_yzz_xyyyyyz_1, \
                                         tg_yzz_xyyyyyzz_0, tg_yzz_xyyyyyzz_1, tg_yzz_xyyyyzz_1, tg_yzz_xyyyyzzz_0, \
                                         tg_yzz_xyyyyzzz_1, tg_yzz_xyyyzzz_1, tg_yzz_xyyyzzzz_0, tg_yzz_xyyyzzzz_1, \
                                         tg_yzz_xyyzzzz_1, tg_yzz_xyyzzzzz_0, tg_yzz_xyyzzzzz_1, tg_yzz_xyzzzzz_1, \
                                         tg_yzz_xyzzzzzz_0, tg_yzz_xyzzzzzz_1, tg_yzz_xzzzzzz_1, tg_yzz_xzzzzzzz_0, \
                                         tg_yzz_xzzzzzzz_1, tg_yzz_yyyyyyy_1, tg_yzz_yyyyyyyy_0, tg_yzz_yyyyyyyy_1, \
                                         tg_yzz_yyyyyyyz_0, tg_yzz_yyyyyyyz_1, tg_yzz_yyyyyyz_1, tg_yzz_yyyyyyzz_0, \
                                         tg_yzz_yyyyyyzz_1, tg_yzz_yyyyyzz_1, tg_zz_xxxxxxxx_0, tg_zz_xxxxxxxx_1, \
                                         tg_zz_xxxxxxxy_0, tg_zz_xxxxxxxy_1, tg_zz_xxxxxxxz_0, tg_zz_xxxxxxxz_1, \
                                         tg_zz_xxxxxxyy_0, tg_zz_xxxxxxyy_1, tg_zz_xxxxxxyz_0, tg_zz_xxxxxxyz_1, \
                                         tg_zz_xxxxxxzz_0, tg_zz_xxxxxxzz_1, tg_zz_xxxxxyyy_0, tg_zz_xxxxxyyy_1, \
                                         tg_zz_xxxxxyyz_0, tg_zz_xxxxxyyz_1, tg_zz_xxxxxyzz_0, tg_zz_xxxxxyzz_1, \
                                         tg_zz_xxxxxzzz_0, tg_zz_xxxxxzzz_1, tg_zz_xxxxyyyy_0, tg_zz_xxxxyyyy_1, \
                                         tg_zz_xxxxyyyz_0, tg_zz_xxxxyyyz_1, tg_zz_xxxxyyzz_0, tg_zz_xxxxyyzz_1, \
                                         tg_zz_xxxxyzzz_0, tg_zz_xxxxyzzz_1, tg_zz_xxxxzzzz_0, tg_zz_xxxxzzzz_1, \
                                         tg_zz_xxxyyyyy_0, tg_zz_xxxyyyyy_1, tg_zz_xxxyyyyz_0, tg_zz_xxxyyyyz_1, \
                                         tg_zz_xxxyyyzz_0, tg_zz_xxxyyyzz_1, tg_zz_xxxyyzzz_0, tg_zz_xxxyyzzz_1, \
                                         tg_zz_xxxyzzzz_0, tg_zz_xxxyzzzz_1, tg_zz_xxxzzzzz_0, tg_zz_xxxzzzzz_1, \
                                         tg_zz_xxyyyyyy_0, tg_zz_xxyyyyyy_1, tg_zz_xxyyyyyz_0, tg_zz_xxyyyyyz_1, \
                                         tg_zz_xxyyyyzz_0, tg_zz_xxyyyyzz_1, tg_zz_xxyyyzzz_0, tg_zz_xxyyyzzz_1, \
                                         tg_zz_xxyyzzzz_0, tg_zz_xxyyzzzz_1, tg_zz_xxyzzzzz_0, tg_zz_xxyzzzzz_1, \
                                         tg_zz_xxzzzzzz_0, tg_zz_xxzzzzzz_1, tg_zz_xyyyyyyy_0, tg_zz_xyyyyyyy_1, \
                                         tg_zz_xyyyyyyz_0, tg_zz_xyyyyyyz_1, tg_zz_xyyyyyzz_0, tg_zz_xyyyyyzz_1, \
                                         tg_zz_xyyyyzzz_0, tg_zz_xyyyyzzz_1, tg_zz_xyyyzzzz_0, tg_zz_xyyyzzzz_1, \
                                         tg_zz_xyyzzzzz_0, tg_zz_xyyzzzzz_1, tg_zz_xyzzzzzz_0, tg_zz_xyzzzzzz_1, \
                                         tg_zz_xzzzzzzz_0, tg_zz_xzzzzzzz_1, tg_zz_yyyyyyyy_0, tg_zz_yyyyyyyy_1, \
                                         tg_zz_yyyyyyyz_0, tg_zz_yyyyyyyz_1, tg_zz_yyyyyyzz_0, tg_zz_yyyyyyzz_1, wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_yyyz_yyyyyyyy_0[j] = pb_y * tg_yyz_yyyyyyyy_0[j] + wp_y[j] * tg_yyz_yyyyyyyy_1[j] + fl1_fx * tg_yz_yyyyyyyy_0[j] - fl1_fx * fl1_fza * tg_yz_yyyyyyyy_1[j] + 4.0 * fl1_fxn * tg_yyz_yyyyyyy_1[j];

                    tg_yyyz_yyyyyyyz_0[j] = pb_y * tg_yyz_yyyyyyyz_0[j] + wp_y[j] * tg_yyz_yyyyyyyz_1[j] + fl1_fx * tg_yz_yyyyyyyz_0[j] - fl1_fx * fl1_fza * tg_yz_yyyyyyyz_1[j] + 3.5 * fl1_fxn * tg_yyz_yyyyyyz_1[j];

                    tg_yyyz_yyyyyyzz_0[j] = pb_y * tg_yyz_yyyyyyzz_0[j] + wp_y[j] * tg_yyz_yyyyyyzz_1[j] + fl1_fx * tg_yz_yyyyyyzz_0[j] - fl1_fx * fl1_fza * tg_yz_yyyyyyzz_1[j] + 3.0 * fl1_fxn * tg_yyz_yyyyyzz_1[j];

                    tg_yyyz_yyyyyzzz_0[j] = pb_y * tg_yyz_yyyyyzzz_0[j] + wp_y[j] * tg_yyz_yyyyyzzz_1[j] + fl1_fx * tg_yz_yyyyyzzz_0[j] - fl1_fx * fl1_fza * tg_yz_yyyyyzzz_1[j] + 2.5 * fl1_fxn * tg_yyz_yyyyzzz_1[j];

                    tg_yyyz_yyyyzzzz_0[j] = pb_y * tg_yyz_yyyyzzzz_0[j] + wp_y[j] * tg_yyz_yyyyzzzz_1[j] + fl1_fx * tg_yz_yyyyzzzz_0[j] - fl1_fx * fl1_fza * tg_yz_yyyyzzzz_1[j] + 2.0 * fl1_fxn * tg_yyz_yyyzzzz_1[j];

                    tg_yyyz_yyyzzzzz_0[j] = pb_y * tg_yyz_yyyzzzzz_0[j] + wp_y[j] * tg_yyz_yyyzzzzz_1[j] + fl1_fx * tg_yz_yyyzzzzz_0[j] - fl1_fx * fl1_fza * tg_yz_yyyzzzzz_1[j] + 1.5 * fl1_fxn * tg_yyz_yyzzzzz_1[j];

                    tg_yyyz_yyzzzzzz_0[j] = pb_y * tg_yyz_yyzzzzzz_0[j] + wp_y[j] * tg_yyz_yyzzzzzz_1[j] + fl1_fx * tg_yz_yyzzzzzz_0[j] - fl1_fx * fl1_fza * tg_yz_yyzzzzzz_1[j] + fl1_fxn * tg_yyz_yzzzzzz_1[j];

                    tg_yyyz_yzzzzzzz_0[j] = pb_y * tg_yyz_yzzzzzzz_0[j] + wp_y[j] * tg_yyz_yzzzzzzz_1[j] + fl1_fx * tg_yz_yzzzzzzz_0[j] - fl1_fx * fl1_fza * tg_yz_yzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyz_zzzzzzz_1[j];

                    tg_yyyz_zzzzzzzz_0[j] = pb_y * tg_yyz_zzzzzzzz_0[j] + wp_y[j] * tg_yyz_zzzzzzzz_1[j] + fl1_fx * tg_yz_zzzzzzzz_0[j] - fl1_fx * fl1_fza * tg_yz_zzzzzzzz_1[j];

                    tg_yyzz_xxxxxxxx_0[j] = pb_y * tg_yzz_xxxxxxxx_0[j] + wp_y[j] * tg_yzz_xxxxxxxx_1[j] + 0.5 * fl1_fx * tg_zz_xxxxxxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxxxxx_1[j];

                    tg_yyzz_xxxxxxxy_0[j] = pb_y * tg_yzz_xxxxxxxy_0[j] + wp_y[j] * tg_yzz_xxxxxxxy_1[j] + 0.5 * fl1_fx * tg_zz_xxxxxxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxxxxy_1[j] + 0.5 * fl1_fxn * tg_yzz_xxxxxxx_1[j];

                    tg_yyzz_xxxxxxxz_0[j] = pb_y * tg_yzz_xxxxxxxz_0[j] + wp_y[j] * tg_yzz_xxxxxxxz_1[j] + 0.5 * fl1_fx * tg_zz_xxxxxxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxxxxz_1[j];

                    tg_yyzz_xxxxxxyy_0[j] = pb_y * tg_yzz_xxxxxxyy_0[j] + wp_y[j] * tg_yzz_xxxxxxyy_1[j] + 0.5 * fl1_fx * tg_zz_xxxxxxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxxxyy_1[j] + fl1_fxn * tg_yzz_xxxxxxy_1[j];

                    tg_yyzz_xxxxxxyz_0[j] = pb_y * tg_yzz_xxxxxxyz_0[j] + wp_y[j] * tg_yzz_xxxxxxyz_1[j] + 0.5 * fl1_fx * tg_zz_xxxxxxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxxxyz_1[j] + 0.5 * fl1_fxn * tg_yzz_xxxxxxz_1[j];

                    tg_yyzz_xxxxxxzz_0[j] = pb_y * tg_yzz_xxxxxxzz_0[j] + wp_y[j] * tg_yzz_xxxxxxzz_1[j] + 0.5 * fl1_fx * tg_zz_xxxxxxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxxxzz_1[j];

                    tg_yyzz_xxxxxyyy_0[j] = pb_y * tg_yzz_xxxxxyyy_0[j] + wp_y[j] * tg_yzz_xxxxxyyy_1[j] + 0.5 * fl1_fx * tg_zz_xxxxxyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxxyyy_1[j] + 1.5 * fl1_fxn * tg_yzz_xxxxxyy_1[j];

                    tg_yyzz_xxxxxyyz_0[j] = pb_y * tg_yzz_xxxxxyyz_0[j] + wp_y[j] * tg_yzz_xxxxxyyz_1[j] + 0.5 * fl1_fx * tg_zz_xxxxxyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxxyyz_1[j] + fl1_fxn * tg_yzz_xxxxxyz_1[j];

                    tg_yyzz_xxxxxyzz_0[j] = pb_y * tg_yzz_xxxxxyzz_0[j] + wp_y[j] * tg_yzz_xxxxxyzz_1[j] + 0.5 * fl1_fx * tg_zz_xxxxxyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxxyzz_1[j] + 0.5 * fl1_fxn * tg_yzz_xxxxxzz_1[j];

                    tg_yyzz_xxxxxzzz_0[j] = pb_y * tg_yzz_xxxxxzzz_0[j] + wp_y[j] * tg_yzz_xxxxxzzz_1[j] + 0.5 * fl1_fx * tg_zz_xxxxxzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxxzzz_1[j];

                    tg_yyzz_xxxxyyyy_0[j] = pb_y * tg_yzz_xxxxyyyy_0[j] + wp_y[j] * tg_yzz_xxxxyyyy_1[j] + 0.5 * fl1_fx * tg_zz_xxxxyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_yzz_xxxxyyy_1[j];

                    tg_yyzz_xxxxyyyz_0[j] = pb_y * tg_yzz_xxxxyyyz_0[j] + wp_y[j] * tg_yzz_xxxxyyyz_1[j] + 0.5 * fl1_fx * tg_zz_xxxxyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxyyyz_1[j] + 1.5 * fl1_fxn * tg_yzz_xxxxyyz_1[j];

                    tg_yyzz_xxxxyyzz_0[j] = pb_y * tg_yzz_xxxxyyzz_0[j] + wp_y[j] * tg_yzz_xxxxyyzz_1[j] + 0.5 * fl1_fx * tg_zz_xxxxyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxyyzz_1[j] + fl1_fxn * tg_yzz_xxxxyzz_1[j];

                    tg_yyzz_xxxxyzzz_0[j] = pb_y * tg_yzz_xxxxyzzz_0[j] + wp_y[j] * tg_yzz_xxxxyzzz_1[j] + 0.5 * fl1_fx * tg_zz_xxxxyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxyzzz_1[j] + 0.5 * fl1_fxn * tg_yzz_xxxxzzz_1[j];

                    tg_yyzz_xxxxzzzz_0[j] = pb_y * tg_yzz_xxxxzzzz_0[j] + wp_y[j] * tg_yzz_xxxxzzzz_1[j] + 0.5 * fl1_fx * tg_zz_xxxxzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxxzzzz_1[j];

                    tg_yyzz_xxxyyyyy_0[j] = pb_y * tg_yzz_xxxyyyyy_0[j] + wp_y[j] * tg_yzz_xxxyyyyy_1[j] + 0.5 * fl1_fx * tg_zz_xxxyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxyyyyy_1[j] + 2.5 * fl1_fxn * tg_yzz_xxxyyyy_1[j];

                    tg_yyzz_xxxyyyyz_0[j] = pb_y * tg_yzz_xxxyyyyz_0[j] + wp_y[j] * tg_yzz_xxxyyyyz_1[j] + 0.5 * fl1_fx * tg_zz_xxxyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxyyyyz_1[j] + 2.0 * fl1_fxn * tg_yzz_xxxyyyz_1[j];

                    tg_yyzz_xxxyyyzz_0[j] = pb_y * tg_yzz_xxxyyyzz_0[j] + wp_y[j] * tg_yzz_xxxyyyzz_1[j] + 0.5 * fl1_fx * tg_zz_xxxyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_yzz_xxxyyzz_1[j];

                    tg_yyzz_xxxyyzzz_0[j] = pb_y * tg_yzz_xxxyyzzz_0[j] + wp_y[j] * tg_yzz_xxxyyzzz_1[j] + 0.5 * fl1_fx * tg_zz_xxxyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxyyzzz_1[j] + fl1_fxn * tg_yzz_xxxyzzz_1[j];

                    tg_yyzz_xxxyzzzz_0[j] = pb_y * tg_yzz_xxxyzzzz_0[j] + wp_y[j] * tg_yzz_xxxyzzzz_1[j] + 0.5 * fl1_fx * tg_zz_xxxyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxyzzzz_1[j] + 0.5 * fl1_fxn * tg_yzz_xxxzzzz_1[j];

                    tg_yyzz_xxxzzzzz_0[j] = pb_y * tg_yzz_xxxzzzzz_0[j] + wp_y[j] * tg_yzz_xxxzzzzz_1[j] + 0.5 * fl1_fx * tg_zz_xxxzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxzzzzz_1[j];

                    tg_yyzz_xxyyyyyy_0[j] = pb_y * tg_yzz_xxyyyyyy_0[j] + wp_y[j] * tg_yzz_xxyyyyyy_1[j] + 0.5 * fl1_fx * tg_zz_xxyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxyyyyyy_1[j] + 3.0 * fl1_fxn * tg_yzz_xxyyyyy_1[j];

                    tg_yyzz_xxyyyyyz_0[j] = pb_y * tg_yzz_xxyyyyyz_0[j] + wp_y[j] * tg_yzz_xxyyyyyz_1[j] + 0.5 * fl1_fx * tg_zz_xxyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxyyyyyz_1[j] + 2.5 * fl1_fxn * tg_yzz_xxyyyyz_1[j];

                    tg_yyzz_xxyyyyzz_0[j] = pb_y * tg_yzz_xxyyyyzz_0[j] + wp_y[j] * tg_yzz_xxyyyyzz_1[j] + 0.5 * fl1_fx * tg_zz_xxyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxyyyyzz_1[j] + 2.0 * fl1_fxn * tg_yzz_xxyyyzz_1[j];

                    tg_yyzz_xxyyyzzz_0[j] = pb_y * tg_yzz_xxyyyzzz_0[j] + wp_y[j] * tg_yzz_xxyyyzzz_1[j] + 0.5 * fl1_fx * tg_zz_xxyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxyyyzzz_1[j] + 1.5 * fl1_fxn * tg_yzz_xxyyzzz_1[j];

                    tg_yyzz_xxyyzzzz_0[j] = pb_y * tg_yzz_xxyyzzzz_0[j] + wp_y[j] * tg_yzz_xxyyzzzz_1[j] + 0.5 * fl1_fx * tg_zz_xxyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxyyzzzz_1[j] + fl1_fxn * tg_yzz_xxyzzzz_1[j];

                    tg_yyzz_xxyzzzzz_0[j] = pb_y * tg_yzz_xxyzzzzz_0[j] + wp_y[j] * tg_yzz_xxyzzzzz_1[j] + 0.5 * fl1_fx * tg_zz_xxyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxyzzzzz_1[j] + 0.5 * fl1_fxn * tg_yzz_xxzzzzz_1[j];

                    tg_yyzz_xxzzzzzz_0[j] = pb_y * tg_yzz_xxzzzzzz_0[j] + wp_y[j] * tg_yzz_xxzzzzzz_1[j] + 0.5 * fl1_fx * tg_zz_xxzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxzzzzzz_1[j];

                    tg_yyzz_xyyyyyyy_0[j] = pb_y * tg_yzz_xyyyyyyy_0[j] + wp_y[j] * tg_yzz_xyyyyyyy_1[j] + 0.5 * fl1_fx * tg_zz_xyyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xyyyyyyy_1[j] + 3.5 * fl1_fxn * tg_yzz_xyyyyyy_1[j];

                    tg_yyzz_xyyyyyyz_0[j] = pb_y * tg_yzz_xyyyyyyz_0[j] + wp_y[j] * tg_yzz_xyyyyyyz_1[j] + 0.5 * fl1_fx * tg_zz_xyyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xyyyyyyz_1[j] + 3.0 * fl1_fxn * tg_yzz_xyyyyyz_1[j];

                    tg_yyzz_xyyyyyzz_0[j] = pb_y * tg_yzz_xyyyyyzz_0[j] + wp_y[j] * tg_yzz_xyyyyyzz_1[j] + 0.5 * fl1_fx * tg_zz_xyyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xyyyyyzz_1[j] + 2.5 * fl1_fxn * tg_yzz_xyyyyzz_1[j];

                    tg_yyzz_xyyyyzzz_0[j] = pb_y * tg_yzz_xyyyyzzz_0[j] + wp_y[j] * tg_yzz_xyyyyzzz_1[j] + 0.5 * fl1_fx * tg_zz_xyyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xyyyyzzz_1[j] + 2.0 * fl1_fxn * tg_yzz_xyyyzzz_1[j];

                    tg_yyzz_xyyyzzzz_0[j] = pb_y * tg_yzz_xyyyzzzz_0[j] + wp_y[j] * tg_yzz_xyyyzzzz_1[j] + 0.5 * fl1_fx * tg_zz_xyyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xyyyzzzz_1[j] + 1.5 * fl1_fxn * tg_yzz_xyyzzzz_1[j];

                    tg_yyzz_xyyzzzzz_0[j] = pb_y * tg_yzz_xyyzzzzz_0[j] + wp_y[j] * tg_yzz_xyyzzzzz_1[j] + 0.5 * fl1_fx * tg_zz_xyyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xyyzzzzz_1[j] + fl1_fxn * tg_yzz_xyzzzzz_1[j];

                    tg_yyzz_xyzzzzzz_0[j] = pb_y * tg_yzz_xyzzzzzz_0[j] + wp_y[j] * tg_yzz_xyzzzzzz_1[j] + 0.5 * fl1_fx * tg_zz_xyzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yzz_xzzzzzz_1[j];

                    tg_yyzz_xzzzzzzz_0[j] = pb_y * tg_yzz_xzzzzzzz_0[j] + wp_y[j] * tg_yzz_xzzzzzzz_1[j] + 0.5 * fl1_fx * tg_zz_xzzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xzzzzzzz_1[j];

                    tg_yyzz_yyyyyyyy_0[j] = pb_y * tg_yzz_yyyyyyyy_0[j] + wp_y[j] * tg_yzz_yyyyyyyy_1[j] + 0.5 * fl1_fx * tg_zz_yyyyyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_yyyyyyyy_1[j] + 4.0 * fl1_fxn * tg_yzz_yyyyyyy_1[j];

                    tg_yyzz_yyyyyyyz_0[j] = pb_y * tg_yzz_yyyyyyyz_0[j] + wp_y[j] * tg_yzz_yyyyyyyz_1[j] + 0.5 * fl1_fx * tg_zz_yyyyyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_yyyyyyyz_1[j] + 3.5 * fl1_fxn * tg_yzz_yyyyyyz_1[j];

                    tg_yyzz_yyyyyyzz_0[j] = pb_y * tg_yzz_yyyyyyzz_0[j] + wp_y[j] * tg_yzz_yyyyyyzz_1[j] + 0.5 * fl1_fx * tg_zz_yyyyyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_yyyyyyzz_1[j] + 3.0 * fl1_fxn * tg_yzz_yyyyyzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSL_579_627(      CMemBlock2D<double>& primBuffer,
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
                                             {4, -1, -1, -1},
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_2_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_2_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

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

                auto pb_y = r_pb_y[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_y = wpDistances.data(3 * idx + 1);

                // set up pointers to auxilary integrals

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

                auto tg_yzz_yyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 399); 

                auto tg_yzz_yyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 400); 

                auto tg_yzz_yyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 401); 

                auto tg_yzz_yyzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 402); 

                auto tg_yzz_yzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 403); 

                auto tg_yzz_zzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 404); 

                auto tg_zzz_xxxxxxxx_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 405); 

                auto tg_zzz_xxxxxxxy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 406); 

                auto tg_zzz_xxxxxxxz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 407); 

                auto tg_zzz_xxxxxxyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 408); 

                auto tg_zzz_xxxxxxyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 409); 

                auto tg_zzz_xxxxxxzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 410); 

                auto tg_zzz_xxxxxyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 411); 

                auto tg_zzz_xxxxxyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 412); 

                auto tg_zzz_xxxxxyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 413); 

                auto tg_zzz_xxxxxzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 414); 

                auto tg_zzz_xxxxyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 415); 

                auto tg_zzz_xxxxyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 416); 

                auto tg_zzz_xxxxyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 417); 

                auto tg_zzz_xxxxyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 418); 

                auto tg_zzz_xxxxzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 419); 

                auto tg_zzz_xxxyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 420); 

                auto tg_zzz_xxxyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 421); 

                auto tg_zzz_xxxyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 422); 

                auto tg_zzz_xxxyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 423); 

                auto tg_zzz_xxxyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 424); 

                auto tg_zzz_xxxzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 425); 

                auto tg_zzz_xxyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 426); 

                auto tg_zzz_xxyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 427); 

                auto tg_zzz_xxyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 428); 

                auto tg_zzz_xxyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 429); 

                auto tg_zzz_xxyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 430); 

                auto tg_zzz_xxyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 431); 

                auto tg_zzz_xxzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 432); 

                auto tg_zzz_xyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 433); 

                auto tg_zzz_xyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 434); 

                auto tg_zzz_xyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 435); 

                auto tg_zzz_xyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 436); 

                auto tg_zzz_xyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 437); 

                auto tg_zzz_xyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 438); 

                auto tg_zzz_xyzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 439); 

                auto tg_zzz_xzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 440); 

                auto tg_zzz_yyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 441); 

                auto tg_zzz_yyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 442); 

                auto tg_zzz_yyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 443); 

                auto tg_zzz_yyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 444); 

                auto tg_zzz_yyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 445); 

                auto tg_zzz_yyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 446); 

                auto tg_zz_yyyyyzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 264); 

                auto tg_zz_yyyyzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 265); 

                auto tg_zz_yyyzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 266); 

                auto tg_zz_yyzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 267); 

                auto tg_zz_yzzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 268); 

                auto tg_zz_zzzzzzzz_0 = primBuffer.data(pidx_g_2_8_m0 + 270 * idx + 269); 

                auto tg_zz_yyyyyzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 264); 

                auto tg_zz_yyyyzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 265); 

                auto tg_zz_yyyzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 266); 

                auto tg_zz_yyzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 267); 

                auto tg_zz_yzzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 268); 

                auto tg_zz_zzzzzzzz_1 = primBuffer.data(pidx_g_2_8_m1 + 270 * idx + 269); 

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

                // set up pointers to integrals

                auto tg_yyzz_yyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 579); 

                auto tg_yyzz_yyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 580); 

                auto tg_yyzz_yyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 581); 

                auto tg_yyzz_yyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 582); 

                auto tg_yyzz_yzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 583); 

                auto tg_yyzz_zzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 584); 

                auto tg_yzzz_xxxxxxxx_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 585); 

                auto tg_yzzz_xxxxxxxy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 586); 

                auto tg_yzzz_xxxxxxxz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 587); 

                auto tg_yzzz_xxxxxxyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 588); 

                auto tg_yzzz_xxxxxxyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 589); 

                auto tg_yzzz_xxxxxxzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 590); 

                auto tg_yzzz_xxxxxyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 591); 

                auto tg_yzzz_xxxxxyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 592); 

                auto tg_yzzz_xxxxxyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 593); 

                auto tg_yzzz_xxxxxzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 594); 

                auto tg_yzzz_xxxxyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 595); 

                auto tg_yzzz_xxxxyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 596); 

                auto tg_yzzz_xxxxyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 597); 

                auto tg_yzzz_xxxxyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 598); 

                auto tg_yzzz_xxxxzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 599); 

                auto tg_yzzz_xxxyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 600); 

                auto tg_yzzz_xxxyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 601); 

                auto tg_yzzz_xxxyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 602); 

                auto tg_yzzz_xxxyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 603); 

                auto tg_yzzz_xxxyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 604); 

                auto tg_yzzz_xxxzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 605); 

                auto tg_yzzz_xxyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 606); 

                auto tg_yzzz_xxyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 607); 

                auto tg_yzzz_xxyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 608); 

                auto tg_yzzz_xxyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 609); 

                auto tg_yzzz_xxyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 610); 

                auto tg_yzzz_xxyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 611); 

                auto tg_yzzz_xxzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 612); 

                auto tg_yzzz_xyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 613); 

                auto tg_yzzz_xyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 614); 

                auto tg_yzzz_xyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 615); 

                auto tg_yzzz_xyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 616); 

                auto tg_yzzz_xyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 617); 

                auto tg_yzzz_xyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 618); 

                auto tg_yzzz_xyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 619); 

                auto tg_yzzz_xzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 620); 

                auto tg_yzzz_yyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 621); 

                auto tg_yzzz_yyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 622); 

                auto tg_yzzz_yyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 623); 

                auto tg_yzzz_yyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 624); 

                auto tg_yzzz_yyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 625); 

                auto tg_yzzz_yyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 626); 

                // Batch of Integrals (579,627)

                #pragma omp simd aligned(fxn, fza, tg_yyzz_yyyyyzzz_0, tg_yyzz_yyyyzzzz_0, tg_yyzz_yyyzzzzz_0, \
                                         tg_yyzz_yyzzzzzz_0, tg_yyzz_yzzzzzzz_0, tg_yyzz_zzzzzzzz_0, tg_yzz_yyyyyzzz_0, \
                                         tg_yzz_yyyyyzzz_1, tg_yzz_yyyyzzz_1, tg_yzz_yyyyzzzz_0, tg_yzz_yyyyzzzz_1, \
                                         tg_yzz_yyyzzzz_1, tg_yzz_yyyzzzzz_0, tg_yzz_yyyzzzzz_1, tg_yzz_yyzzzzz_1, \
                                         tg_yzz_yyzzzzzz_0, tg_yzz_yyzzzzzz_1, tg_yzz_yzzzzzz_1, tg_yzz_yzzzzzzz_0, \
                                         tg_yzz_yzzzzzzz_1, tg_yzz_zzzzzzz_1, tg_yzz_zzzzzzzz_0, tg_yzz_zzzzzzzz_1, \
                                         tg_yzzz_xxxxxxxx_0, tg_yzzz_xxxxxxxy_0, tg_yzzz_xxxxxxxz_0, tg_yzzz_xxxxxxyy_0, \
                                         tg_yzzz_xxxxxxyz_0, tg_yzzz_xxxxxxzz_0, tg_yzzz_xxxxxyyy_0, tg_yzzz_xxxxxyyz_0, \
                                         tg_yzzz_xxxxxyzz_0, tg_yzzz_xxxxxzzz_0, tg_yzzz_xxxxyyyy_0, tg_yzzz_xxxxyyyz_0, \
                                         tg_yzzz_xxxxyyzz_0, tg_yzzz_xxxxyzzz_0, tg_yzzz_xxxxzzzz_0, tg_yzzz_xxxyyyyy_0, \
                                         tg_yzzz_xxxyyyyz_0, tg_yzzz_xxxyyyzz_0, tg_yzzz_xxxyyzzz_0, tg_yzzz_xxxyzzzz_0, \
                                         tg_yzzz_xxxzzzzz_0, tg_yzzz_xxyyyyyy_0, tg_yzzz_xxyyyyyz_0, tg_yzzz_xxyyyyzz_0, \
                                         tg_yzzz_xxyyyzzz_0, tg_yzzz_xxyyzzzz_0, tg_yzzz_xxyzzzzz_0, tg_yzzz_xxzzzzzz_0, \
                                         tg_yzzz_xyyyyyyy_0, tg_yzzz_xyyyyyyz_0, tg_yzzz_xyyyyyzz_0, tg_yzzz_xyyyyzzz_0, \
                                         tg_yzzz_xyyyzzzz_0, tg_yzzz_xyyzzzzz_0, tg_yzzz_xyzzzzzz_0, tg_yzzz_xzzzzzzz_0, \
                                         tg_yzzz_yyyyyyyy_0, tg_yzzz_yyyyyyyz_0, tg_yzzz_yyyyyyzz_0, tg_yzzz_yyyyyzzz_0, \
                                         tg_yzzz_yyyyzzzz_0, tg_yzzz_yyyzzzzz_0, tg_zz_yyyyyzzz_0, tg_zz_yyyyyzzz_1, \
                                         tg_zz_yyyyzzzz_0, tg_zz_yyyyzzzz_1, tg_zz_yyyzzzzz_0, tg_zz_yyyzzzzz_1, \
                                         tg_zz_yyzzzzzz_0, tg_zz_yyzzzzzz_1, tg_zz_yzzzzzzz_0, tg_zz_yzzzzzzz_1, \
                                         tg_zz_zzzzzzzz_0, tg_zz_zzzzzzzz_1, tg_zzz_xxxxxxx_1, tg_zzz_xxxxxxxx_0, \
                                         tg_zzz_xxxxxxxx_1, tg_zzz_xxxxxxxy_0, tg_zzz_xxxxxxxy_1, tg_zzz_xxxxxxxz_0, \
                                         tg_zzz_xxxxxxxz_1, tg_zzz_xxxxxxy_1, tg_zzz_xxxxxxyy_0, tg_zzz_xxxxxxyy_1, \
                                         tg_zzz_xxxxxxyz_0, tg_zzz_xxxxxxyz_1, tg_zzz_xxxxxxz_1, tg_zzz_xxxxxxzz_0, \
                                         tg_zzz_xxxxxxzz_1, tg_zzz_xxxxxyy_1, tg_zzz_xxxxxyyy_0, tg_zzz_xxxxxyyy_1, \
                                         tg_zzz_xxxxxyyz_0, tg_zzz_xxxxxyyz_1, tg_zzz_xxxxxyz_1, tg_zzz_xxxxxyzz_0, \
                                         tg_zzz_xxxxxyzz_1, tg_zzz_xxxxxzz_1, tg_zzz_xxxxxzzz_0, tg_zzz_xxxxxzzz_1, \
                                         tg_zzz_xxxxyyy_1, tg_zzz_xxxxyyyy_0, tg_zzz_xxxxyyyy_1, tg_zzz_xxxxyyyz_0, \
                                         tg_zzz_xxxxyyyz_1, tg_zzz_xxxxyyz_1, tg_zzz_xxxxyyzz_0, tg_zzz_xxxxyyzz_1, \
                                         tg_zzz_xxxxyzz_1, tg_zzz_xxxxyzzz_0, tg_zzz_xxxxyzzz_1, tg_zzz_xxxxzzz_1, \
                                         tg_zzz_xxxxzzzz_0, tg_zzz_xxxxzzzz_1, tg_zzz_xxxyyyy_1, tg_zzz_xxxyyyyy_0, \
                                         tg_zzz_xxxyyyyy_1, tg_zzz_xxxyyyyz_0, tg_zzz_xxxyyyyz_1, tg_zzz_xxxyyyz_1, \
                                         tg_zzz_xxxyyyzz_0, tg_zzz_xxxyyyzz_1, tg_zzz_xxxyyzz_1, tg_zzz_xxxyyzzz_0, \
                                         tg_zzz_xxxyyzzz_1, tg_zzz_xxxyzzz_1, tg_zzz_xxxyzzzz_0, tg_zzz_xxxyzzzz_1, \
                                         tg_zzz_xxxzzzz_1, tg_zzz_xxxzzzzz_0, tg_zzz_xxxzzzzz_1, tg_zzz_xxyyyyy_1, \
                                         tg_zzz_xxyyyyyy_0, tg_zzz_xxyyyyyy_1, tg_zzz_xxyyyyyz_0, tg_zzz_xxyyyyyz_1, \
                                         tg_zzz_xxyyyyz_1, tg_zzz_xxyyyyzz_0, tg_zzz_xxyyyyzz_1, tg_zzz_xxyyyzz_1, \
                                         tg_zzz_xxyyyzzz_0, tg_zzz_xxyyyzzz_1, tg_zzz_xxyyzzz_1, tg_zzz_xxyyzzzz_0, \
                                         tg_zzz_xxyyzzzz_1, tg_zzz_xxyzzzz_1, tg_zzz_xxyzzzzz_0, tg_zzz_xxyzzzzz_1, \
                                         tg_zzz_xxzzzzz_1, tg_zzz_xxzzzzzz_0, tg_zzz_xxzzzzzz_1, tg_zzz_xyyyyyy_1, \
                                         tg_zzz_xyyyyyyy_0, tg_zzz_xyyyyyyy_1, tg_zzz_xyyyyyyz_0, tg_zzz_xyyyyyyz_1, \
                                         tg_zzz_xyyyyyz_1, tg_zzz_xyyyyyzz_0, tg_zzz_xyyyyyzz_1, tg_zzz_xyyyyzz_1, \
                                         tg_zzz_xyyyyzzz_0, tg_zzz_xyyyyzzz_1, tg_zzz_xyyyzzz_1, tg_zzz_xyyyzzzz_0, \
                                         tg_zzz_xyyyzzzz_1, tg_zzz_xyyzzzz_1, tg_zzz_xyyzzzzz_0, tg_zzz_xyyzzzzz_1, \
                                         tg_zzz_xyzzzzz_1, tg_zzz_xyzzzzzz_0, tg_zzz_xyzzzzzz_1, tg_zzz_xzzzzzz_1, \
                                         tg_zzz_xzzzzzzz_0, tg_zzz_xzzzzzzz_1, tg_zzz_yyyyyyy_1, tg_zzz_yyyyyyyy_0, \
                                         tg_zzz_yyyyyyyy_1, tg_zzz_yyyyyyyz_0, tg_zzz_yyyyyyyz_1, tg_zzz_yyyyyyz_1, \
                                         tg_zzz_yyyyyyzz_0, tg_zzz_yyyyyyzz_1, tg_zzz_yyyyyzz_1, tg_zzz_yyyyyzzz_0, \
                                         tg_zzz_yyyyyzzz_1, tg_zzz_yyyyzzz_1, tg_zzz_yyyyzzzz_0, tg_zzz_yyyyzzzz_1, \
                                         tg_zzz_yyyzzzz_1, tg_zzz_yyyzzzzz_0, tg_zzz_yyyzzzzz_1, tg_zzz_yyzzzzz_1, wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_yyzz_yyyyyzzz_0[j] = pb_y * tg_yzz_yyyyyzzz_0[j] + wp_y[j] * tg_yzz_yyyyyzzz_1[j] + 0.5 * fl1_fx * tg_zz_yyyyyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_yyyyyzzz_1[j] + 2.5 * fl1_fxn * tg_yzz_yyyyzzz_1[j];

                    tg_yyzz_yyyyzzzz_0[j] = pb_y * tg_yzz_yyyyzzzz_0[j] + wp_y[j] * tg_yzz_yyyyzzzz_1[j] + 0.5 * fl1_fx * tg_zz_yyyyzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_yyyyzzzz_1[j] + 2.0 * fl1_fxn * tg_yzz_yyyzzzz_1[j];

                    tg_yyzz_yyyzzzzz_0[j] = pb_y * tg_yzz_yyyzzzzz_0[j] + wp_y[j] * tg_yzz_yyyzzzzz_1[j] + 0.5 * fl1_fx * tg_zz_yyyzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_yyyzzzzz_1[j] + 1.5 * fl1_fxn * tg_yzz_yyzzzzz_1[j];

                    tg_yyzz_yyzzzzzz_0[j] = pb_y * tg_yzz_yyzzzzzz_0[j] + wp_y[j] * tg_yzz_yyzzzzzz_1[j] + 0.5 * fl1_fx * tg_zz_yyzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_yyzzzzzz_1[j] + fl1_fxn * tg_yzz_yzzzzzz_1[j];

                    tg_yyzz_yzzzzzzz_0[j] = pb_y * tg_yzz_yzzzzzzz_0[j] + wp_y[j] * tg_yzz_yzzzzzzz_1[j] + 0.5 * fl1_fx * tg_zz_yzzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_yzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yzz_zzzzzzz_1[j];

                    tg_yyzz_zzzzzzzz_0[j] = pb_y * tg_yzz_zzzzzzzz_0[j] + wp_y[j] * tg_yzz_zzzzzzzz_1[j] + 0.5 * fl1_fx * tg_zz_zzzzzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_zzzzzzzz_1[j];

                    tg_yzzz_xxxxxxxx_0[j] = pb_y * tg_zzz_xxxxxxxx_0[j] + wp_y[j] * tg_zzz_xxxxxxxx_1[j];

                    tg_yzzz_xxxxxxxy_0[j] = pb_y * tg_zzz_xxxxxxxy_0[j] + wp_y[j] * tg_zzz_xxxxxxxy_1[j] + 0.5 * fl1_fxn * tg_zzz_xxxxxxx_1[j];

                    tg_yzzz_xxxxxxxz_0[j] = pb_y * tg_zzz_xxxxxxxz_0[j] + wp_y[j] * tg_zzz_xxxxxxxz_1[j];

                    tg_yzzz_xxxxxxyy_0[j] = pb_y * tg_zzz_xxxxxxyy_0[j] + wp_y[j] * tg_zzz_xxxxxxyy_1[j] + fl1_fxn * tg_zzz_xxxxxxy_1[j];

                    tg_yzzz_xxxxxxyz_0[j] = pb_y * tg_zzz_xxxxxxyz_0[j] + wp_y[j] * tg_zzz_xxxxxxyz_1[j] + 0.5 * fl1_fxn * tg_zzz_xxxxxxz_1[j];

                    tg_yzzz_xxxxxxzz_0[j] = pb_y * tg_zzz_xxxxxxzz_0[j] + wp_y[j] * tg_zzz_xxxxxxzz_1[j];

                    tg_yzzz_xxxxxyyy_0[j] = pb_y * tg_zzz_xxxxxyyy_0[j] + wp_y[j] * tg_zzz_xxxxxyyy_1[j] + 1.5 * fl1_fxn * tg_zzz_xxxxxyy_1[j];

                    tg_yzzz_xxxxxyyz_0[j] = pb_y * tg_zzz_xxxxxyyz_0[j] + wp_y[j] * tg_zzz_xxxxxyyz_1[j] + fl1_fxn * tg_zzz_xxxxxyz_1[j];

                    tg_yzzz_xxxxxyzz_0[j] = pb_y * tg_zzz_xxxxxyzz_0[j] + wp_y[j] * tg_zzz_xxxxxyzz_1[j] + 0.5 * fl1_fxn * tg_zzz_xxxxxzz_1[j];

                    tg_yzzz_xxxxxzzz_0[j] = pb_y * tg_zzz_xxxxxzzz_0[j] + wp_y[j] * tg_zzz_xxxxxzzz_1[j];

                    tg_yzzz_xxxxyyyy_0[j] = pb_y * tg_zzz_xxxxyyyy_0[j] + wp_y[j] * tg_zzz_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_zzz_xxxxyyy_1[j];

                    tg_yzzz_xxxxyyyz_0[j] = pb_y * tg_zzz_xxxxyyyz_0[j] + wp_y[j] * tg_zzz_xxxxyyyz_1[j] + 1.5 * fl1_fxn * tg_zzz_xxxxyyz_1[j];

                    tg_yzzz_xxxxyyzz_0[j] = pb_y * tg_zzz_xxxxyyzz_0[j] + wp_y[j] * tg_zzz_xxxxyyzz_1[j] + fl1_fxn * tg_zzz_xxxxyzz_1[j];

                    tg_yzzz_xxxxyzzz_0[j] = pb_y * tg_zzz_xxxxyzzz_0[j] + wp_y[j] * tg_zzz_xxxxyzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_xxxxzzz_1[j];

                    tg_yzzz_xxxxzzzz_0[j] = pb_y * tg_zzz_xxxxzzzz_0[j] + wp_y[j] * tg_zzz_xxxxzzzz_1[j];

                    tg_yzzz_xxxyyyyy_0[j] = pb_y * tg_zzz_xxxyyyyy_0[j] + wp_y[j] * tg_zzz_xxxyyyyy_1[j] + 2.5 * fl1_fxn * tg_zzz_xxxyyyy_1[j];

                    tg_yzzz_xxxyyyyz_0[j] = pb_y * tg_zzz_xxxyyyyz_0[j] + wp_y[j] * tg_zzz_xxxyyyyz_1[j] + 2.0 * fl1_fxn * tg_zzz_xxxyyyz_1[j];

                    tg_yzzz_xxxyyyzz_0[j] = pb_y * tg_zzz_xxxyyyzz_0[j] + wp_y[j] * tg_zzz_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_zzz_xxxyyzz_1[j];

                    tg_yzzz_xxxyyzzz_0[j] = pb_y * tg_zzz_xxxyyzzz_0[j] + wp_y[j] * tg_zzz_xxxyyzzz_1[j] + fl1_fxn * tg_zzz_xxxyzzz_1[j];

                    tg_yzzz_xxxyzzzz_0[j] = pb_y * tg_zzz_xxxyzzzz_0[j] + wp_y[j] * tg_zzz_xxxyzzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_xxxzzzz_1[j];

                    tg_yzzz_xxxzzzzz_0[j] = pb_y * tg_zzz_xxxzzzzz_0[j] + wp_y[j] * tg_zzz_xxxzzzzz_1[j];

                    tg_yzzz_xxyyyyyy_0[j] = pb_y * tg_zzz_xxyyyyyy_0[j] + wp_y[j] * tg_zzz_xxyyyyyy_1[j] + 3.0 * fl1_fxn * tg_zzz_xxyyyyy_1[j];

                    tg_yzzz_xxyyyyyz_0[j] = pb_y * tg_zzz_xxyyyyyz_0[j] + wp_y[j] * tg_zzz_xxyyyyyz_1[j] + 2.5 * fl1_fxn * tg_zzz_xxyyyyz_1[j];

                    tg_yzzz_xxyyyyzz_0[j] = pb_y * tg_zzz_xxyyyyzz_0[j] + wp_y[j] * tg_zzz_xxyyyyzz_1[j] + 2.0 * fl1_fxn * tg_zzz_xxyyyzz_1[j];

                    tg_yzzz_xxyyyzzz_0[j] = pb_y * tg_zzz_xxyyyzzz_0[j] + wp_y[j] * tg_zzz_xxyyyzzz_1[j] + 1.5 * fl1_fxn * tg_zzz_xxyyzzz_1[j];

                    tg_yzzz_xxyyzzzz_0[j] = pb_y * tg_zzz_xxyyzzzz_0[j] + wp_y[j] * tg_zzz_xxyyzzzz_1[j] + fl1_fxn * tg_zzz_xxyzzzz_1[j];

                    tg_yzzz_xxyzzzzz_0[j] = pb_y * tg_zzz_xxyzzzzz_0[j] + wp_y[j] * tg_zzz_xxyzzzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_xxzzzzz_1[j];

                    tg_yzzz_xxzzzzzz_0[j] = pb_y * tg_zzz_xxzzzzzz_0[j] + wp_y[j] * tg_zzz_xxzzzzzz_1[j];

                    tg_yzzz_xyyyyyyy_0[j] = pb_y * tg_zzz_xyyyyyyy_0[j] + wp_y[j] * tg_zzz_xyyyyyyy_1[j] + 3.5 * fl1_fxn * tg_zzz_xyyyyyy_1[j];

                    tg_yzzz_xyyyyyyz_0[j] = pb_y * tg_zzz_xyyyyyyz_0[j] + wp_y[j] * tg_zzz_xyyyyyyz_1[j] + 3.0 * fl1_fxn * tg_zzz_xyyyyyz_1[j];

                    tg_yzzz_xyyyyyzz_0[j] = pb_y * tg_zzz_xyyyyyzz_0[j] + wp_y[j] * tg_zzz_xyyyyyzz_1[j] + 2.5 * fl1_fxn * tg_zzz_xyyyyzz_1[j];

                    tg_yzzz_xyyyyzzz_0[j] = pb_y * tg_zzz_xyyyyzzz_0[j] + wp_y[j] * tg_zzz_xyyyyzzz_1[j] + 2.0 * fl1_fxn * tg_zzz_xyyyzzz_1[j];

                    tg_yzzz_xyyyzzzz_0[j] = pb_y * tg_zzz_xyyyzzzz_0[j] + wp_y[j] * tg_zzz_xyyyzzzz_1[j] + 1.5 * fl1_fxn * tg_zzz_xyyzzzz_1[j];

                    tg_yzzz_xyyzzzzz_0[j] = pb_y * tg_zzz_xyyzzzzz_0[j] + wp_y[j] * tg_zzz_xyyzzzzz_1[j] + fl1_fxn * tg_zzz_xyzzzzz_1[j];

                    tg_yzzz_xyzzzzzz_0[j] = pb_y * tg_zzz_xyzzzzzz_0[j] + wp_y[j] * tg_zzz_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_xzzzzzz_1[j];

                    tg_yzzz_xzzzzzzz_0[j] = pb_y * tg_zzz_xzzzzzzz_0[j] + wp_y[j] * tg_zzz_xzzzzzzz_1[j];

                    tg_yzzz_yyyyyyyy_0[j] = pb_y * tg_zzz_yyyyyyyy_0[j] + wp_y[j] * tg_zzz_yyyyyyyy_1[j] + 4.0 * fl1_fxn * tg_zzz_yyyyyyy_1[j];

                    tg_yzzz_yyyyyyyz_0[j] = pb_y * tg_zzz_yyyyyyyz_0[j] + wp_y[j] * tg_zzz_yyyyyyyz_1[j] + 3.5 * fl1_fxn * tg_zzz_yyyyyyz_1[j];

                    tg_yzzz_yyyyyyzz_0[j] = pb_y * tg_zzz_yyyyyyzz_0[j] + wp_y[j] * tg_zzz_yyyyyyzz_1[j] + 3.0 * fl1_fxn * tg_zzz_yyyyyzz_1[j];

                    tg_yzzz_yyyyyzzz_0[j] = pb_y * tg_zzz_yyyyyzzz_0[j] + wp_y[j] * tg_zzz_yyyyyzzz_1[j] + 2.5 * fl1_fxn * tg_zzz_yyyyzzz_1[j];

                    tg_yzzz_yyyyzzzz_0[j] = pb_y * tg_zzz_yyyyzzzz_0[j] + wp_y[j] * tg_zzz_yyyyzzzz_1[j] + 2.0 * fl1_fxn * tg_zzz_yyyzzzz_1[j];

                    tg_yzzz_yyyzzzzz_0[j] = pb_y * tg_zzz_yyyzzzzz_0[j] + wp_y[j] * tg_zzz_yyyzzzzz_1[j] + 1.5 * fl1_fxn * tg_zzz_yyzzzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSL_627_675(      CMemBlock2D<double>& primBuffer,
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
                                             {4, -1, -1, -1},
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_3_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_2_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_2_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {8, -1, -1, -1}, 
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

                auto tg_zzz_xxxxxxxx_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 405); 

                auto tg_zzz_xxxxxxxy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 406); 

                auto tg_zzz_xxxxxxxz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 407); 

                auto tg_zzz_xxxxxxyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 408); 

                auto tg_zzz_xxxxxxyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 409); 

                auto tg_zzz_xxxxxxzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 410); 

                auto tg_zzz_xxxxxyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 411); 

                auto tg_zzz_xxxxxyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 412); 

                auto tg_zzz_xxxxxyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 413); 

                auto tg_zzz_xxxxxzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 414); 

                auto tg_zzz_xxxxyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 415); 

                auto tg_zzz_xxxxyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 416); 

                auto tg_zzz_xxxxyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 417); 

                auto tg_zzz_xxxxyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 418); 

                auto tg_zzz_xxxxzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 419); 

                auto tg_zzz_xxxyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 420); 

                auto tg_zzz_xxxyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 421); 

                auto tg_zzz_xxxyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 422); 

                auto tg_zzz_xxxyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 423); 

                auto tg_zzz_xxxyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 424); 

                auto tg_zzz_xxxzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 425); 

                auto tg_zzz_xxyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 426); 

                auto tg_zzz_xxyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 427); 

                auto tg_zzz_xxyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 428); 

                auto tg_zzz_xxyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 429); 

                auto tg_zzz_xxyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 430); 

                auto tg_zzz_xxyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 431); 

                auto tg_zzz_xxzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 432); 

                auto tg_zzz_xyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 433); 

                auto tg_zzz_xyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 434); 

                auto tg_zzz_xyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 435); 

                auto tg_zzz_xyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 436); 

                auto tg_zzz_xyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 437); 

                auto tg_zzz_xyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 438); 

                auto tg_zzz_xyzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 439); 

                auto tg_zzz_xzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 440); 

                auto tg_zzz_yyyyyyyy_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 441); 

                auto tg_zzz_yyyyyyyz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 442); 

                auto tg_zzz_yyyyyyzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 443); 

                auto tg_zzz_yyyyyzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 444); 

                auto tg_zzz_yyyyzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 445); 

                auto tg_zzz_yyyzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 446); 

                auto tg_zzz_yyzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 447); 

                auto tg_zzz_yzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 448); 

                auto tg_zzz_zzzzzzzz_1 = primBuffer.data(pidx_g_3_8_m1 + 450 * idx + 449); 

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

                // set up pointers to integrals

                auto tg_yzzz_yyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 627); 

                auto tg_yzzz_yzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 628); 

                auto tg_yzzz_zzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 629); 

                auto tg_zzzz_xxxxxxxx_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 630); 

                auto tg_zzzz_xxxxxxxy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 631); 

                auto tg_zzzz_xxxxxxxz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 632); 

                auto tg_zzzz_xxxxxxyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 633); 

                auto tg_zzzz_xxxxxxyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 634); 

                auto tg_zzzz_xxxxxxzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 635); 

                auto tg_zzzz_xxxxxyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 636); 

                auto tg_zzzz_xxxxxyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 637); 

                auto tg_zzzz_xxxxxyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 638); 

                auto tg_zzzz_xxxxxzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 639); 

                auto tg_zzzz_xxxxyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 640); 

                auto tg_zzzz_xxxxyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 641); 

                auto tg_zzzz_xxxxyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 642); 

                auto tg_zzzz_xxxxyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 643); 

                auto tg_zzzz_xxxxzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 644); 

                auto tg_zzzz_xxxyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 645); 

                auto tg_zzzz_xxxyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 646); 

                auto tg_zzzz_xxxyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 647); 

                auto tg_zzzz_xxxyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 648); 

                auto tg_zzzz_xxxyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 649); 

                auto tg_zzzz_xxxzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 650); 

                auto tg_zzzz_xxyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 651); 

                auto tg_zzzz_xxyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 652); 

                auto tg_zzzz_xxyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 653); 

                auto tg_zzzz_xxyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 654); 

                auto tg_zzzz_xxyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 655); 

                auto tg_zzzz_xxyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 656); 

                auto tg_zzzz_xxzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 657); 

                auto tg_zzzz_xyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 658); 

                auto tg_zzzz_xyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 659); 

                auto tg_zzzz_xyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 660); 

                auto tg_zzzz_xyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 661); 

                auto tg_zzzz_xyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 662); 

                auto tg_zzzz_xyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 663); 

                auto tg_zzzz_xyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 664); 

                auto tg_zzzz_xzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 665); 

                auto tg_zzzz_yyyyyyyy_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 666); 

                auto tg_zzzz_yyyyyyyz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 667); 

                auto tg_zzzz_yyyyyyzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 668); 

                auto tg_zzzz_yyyyyzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 669); 

                auto tg_zzzz_yyyyzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 670); 

                auto tg_zzzz_yyyzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 671); 

                auto tg_zzzz_yyzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 672); 

                auto tg_zzzz_yzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 673); 

                auto tg_zzzz_zzzzzzzz_0 = primBuffer.data(pidx_g_4_8_m0 + 675 * idx + 674); 

                // Batch of Integrals (627,675)

                #pragma omp simd aligned(fxn, fza, tg_yzzz_yyzzzzzz_0, tg_yzzz_yzzzzzzz_0, tg_yzzz_zzzzzzzz_0, \
                                         tg_zz_xxxxxxxx_0, tg_zz_xxxxxxxx_1, tg_zz_xxxxxxxy_0, tg_zz_xxxxxxxy_1, \
                                         tg_zz_xxxxxxxz_0, tg_zz_xxxxxxxz_1, tg_zz_xxxxxxyy_0, tg_zz_xxxxxxyy_1, \
                                         tg_zz_xxxxxxyz_0, tg_zz_xxxxxxyz_1, tg_zz_xxxxxxzz_0, tg_zz_xxxxxxzz_1, \
                                         tg_zz_xxxxxyyy_0, tg_zz_xxxxxyyy_1, tg_zz_xxxxxyyz_0, tg_zz_xxxxxyyz_1, \
                                         tg_zz_xxxxxyzz_0, tg_zz_xxxxxyzz_1, tg_zz_xxxxxzzz_0, tg_zz_xxxxxzzz_1, \
                                         tg_zz_xxxxyyyy_0, tg_zz_xxxxyyyy_1, tg_zz_xxxxyyyz_0, tg_zz_xxxxyyyz_1, \
                                         tg_zz_xxxxyyzz_0, tg_zz_xxxxyyzz_1, tg_zz_xxxxyzzz_0, tg_zz_xxxxyzzz_1, \
                                         tg_zz_xxxxzzzz_0, tg_zz_xxxxzzzz_1, tg_zz_xxxyyyyy_0, tg_zz_xxxyyyyy_1, \
                                         tg_zz_xxxyyyyz_0, tg_zz_xxxyyyyz_1, tg_zz_xxxyyyzz_0, tg_zz_xxxyyyzz_1, \
                                         tg_zz_xxxyyzzz_0, tg_zz_xxxyyzzz_1, tg_zz_xxxyzzzz_0, tg_zz_xxxyzzzz_1, \
                                         tg_zz_xxxzzzzz_0, tg_zz_xxxzzzzz_1, tg_zz_xxyyyyyy_0, tg_zz_xxyyyyyy_1, \
                                         tg_zz_xxyyyyyz_0, tg_zz_xxyyyyyz_1, tg_zz_xxyyyyzz_0, tg_zz_xxyyyyzz_1, \
                                         tg_zz_xxyyyzzz_0, tg_zz_xxyyyzzz_1, tg_zz_xxyyzzzz_0, tg_zz_xxyyzzzz_1, \
                                         tg_zz_xxyzzzzz_0, tg_zz_xxyzzzzz_1, tg_zz_xxzzzzzz_0, tg_zz_xxzzzzzz_1, \
                                         tg_zz_xyyyyyyy_0, tg_zz_xyyyyyyy_1, tg_zz_xyyyyyyz_0, tg_zz_xyyyyyyz_1, \
                                         tg_zz_xyyyyyzz_0, tg_zz_xyyyyyzz_1, tg_zz_xyyyyzzz_0, tg_zz_xyyyyzzz_1, \
                                         tg_zz_xyyyzzzz_0, tg_zz_xyyyzzzz_1, tg_zz_xyyzzzzz_0, tg_zz_xyyzzzzz_1, \
                                         tg_zz_xyzzzzzz_0, tg_zz_xyzzzzzz_1, tg_zz_xzzzzzzz_0, tg_zz_xzzzzzzz_1, \
                                         tg_zz_yyyyyyyy_0, tg_zz_yyyyyyyy_1, tg_zz_yyyyyyyz_0, tg_zz_yyyyyyyz_1, \
                                         tg_zz_yyyyyyzz_0, tg_zz_yyyyyyzz_1, tg_zz_yyyyyzzz_0, tg_zz_yyyyyzzz_1, \
                                         tg_zz_yyyyzzzz_0, tg_zz_yyyyzzzz_1, tg_zz_yyyzzzzz_0, tg_zz_yyyzzzzz_1, \
                                         tg_zz_yyzzzzzz_0, tg_zz_yyzzzzzz_1, tg_zz_yzzzzzzz_0, tg_zz_yzzzzzzz_1, \
                                         tg_zz_zzzzzzzz_0, tg_zz_zzzzzzzz_1, tg_zzz_xxxxxxx_1, tg_zzz_xxxxxxxx_0, \
                                         tg_zzz_xxxxxxxx_1, tg_zzz_xxxxxxxy_0, tg_zzz_xxxxxxxy_1, tg_zzz_xxxxxxxz_0, \
                                         tg_zzz_xxxxxxxz_1, tg_zzz_xxxxxxy_1, tg_zzz_xxxxxxyy_0, tg_zzz_xxxxxxyy_1, \
                                         tg_zzz_xxxxxxyz_0, tg_zzz_xxxxxxyz_1, tg_zzz_xxxxxxz_1, tg_zzz_xxxxxxzz_0, \
                                         tg_zzz_xxxxxxzz_1, tg_zzz_xxxxxyy_1, tg_zzz_xxxxxyyy_0, tg_zzz_xxxxxyyy_1, \
                                         tg_zzz_xxxxxyyz_0, tg_zzz_xxxxxyyz_1, tg_zzz_xxxxxyz_1, tg_zzz_xxxxxyzz_0, \
                                         tg_zzz_xxxxxyzz_1, tg_zzz_xxxxxzz_1, tg_zzz_xxxxxzzz_0, tg_zzz_xxxxxzzz_1, \
                                         tg_zzz_xxxxyyy_1, tg_zzz_xxxxyyyy_0, tg_zzz_xxxxyyyy_1, tg_zzz_xxxxyyyz_0, \
                                         tg_zzz_xxxxyyyz_1, tg_zzz_xxxxyyz_1, tg_zzz_xxxxyyzz_0, tg_zzz_xxxxyyzz_1, \
                                         tg_zzz_xxxxyzz_1, tg_zzz_xxxxyzzz_0, tg_zzz_xxxxyzzz_1, tg_zzz_xxxxzzz_1, \
                                         tg_zzz_xxxxzzzz_0, tg_zzz_xxxxzzzz_1, tg_zzz_xxxyyyy_1, tg_zzz_xxxyyyyy_0, \
                                         tg_zzz_xxxyyyyy_1, tg_zzz_xxxyyyyz_0, tg_zzz_xxxyyyyz_1, tg_zzz_xxxyyyz_1, \
                                         tg_zzz_xxxyyyzz_0, tg_zzz_xxxyyyzz_1, tg_zzz_xxxyyzz_1, tg_zzz_xxxyyzzz_0, \
                                         tg_zzz_xxxyyzzz_1, tg_zzz_xxxyzzz_1, tg_zzz_xxxyzzzz_0, tg_zzz_xxxyzzzz_1, \
                                         tg_zzz_xxxzzzz_1, tg_zzz_xxxzzzzz_0, tg_zzz_xxxzzzzz_1, tg_zzz_xxyyyyy_1, \
                                         tg_zzz_xxyyyyyy_0, tg_zzz_xxyyyyyy_1, tg_zzz_xxyyyyyz_0, tg_zzz_xxyyyyyz_1, \
                                         tg_zzz_xxyyyyz_1, tg_zzz_xxyyyyzz_0, tg_zzz_xxyyyyzz_1, tg_zzz_xxyyyzz_1, \
                                         tg_zzz_xxyyyzzz_0, tg_zzz_xxyyyzzz_1, tg_zzz_xxyyzzz_1, tg_zzz_xxyyzzzz_0, \
                                         tg_zzz_xxyyzzzz_1, tg_zzz_xxyzzzz_1, tg_zzz_xxyzzzzz_0, tg_zzz_xxyzzzzz_1, \
                                         tg_zzz_xxzzzzz_1, tg_zzz_xxzzzzzz_0, tg_zzz_xxzzzzzz_1, tg_zzz_xyyyyyy_1, \
                                         tg_zzz_xyyyyyyy_0, tg_zzz_xyyyyyyy_1, tg_zzz_xyyyyyyz_0, tg_zzz_xyyyyyyz_1, \
                                         tg_zzz_xyyyyyz_1, tg_zzz_xyyyyyzz_0, tg_zzz_xyyyyyzz_1, tg_zzz_xyyyyzz_1, \
                                         tg_zzz_xyyyyzzz_0, tg_zzz_xyyyyzzz_1, tg_zzz_xyyyzzz_1, tg_zzz_xyyyzzzz_0, \
                                         tg_zzz_xyyyzzzz_1, tg_zzz_xyyzzzz_1, tg_zzz_xyyzzzzz_0, tg_zzz_xyyzzzzz_1, \
                                         tg_zzz_xyzzzzz_1, tg_zzz_xyzzzzzz_0, tg_zzz_xyzzzzzz_1, tg_zzz_xzzzzzz_1, \
                                         tg_zzz_xzzzzzzz_0, tg_zzz_xzzzzzzz_1, tg_zzz_yyyyyyy_1, tg_zzz_yyyyyyyy_0, \
                                         tg_zzz_yyyyyyyy_1, tg_zzz_yyyyyyyz_0, tg_zzz_yyyyyyyz_1, tg_zzz_yyyyyyz_1, \
                                         tg_zzz_yyyyyyzz_0, tg_zzz_yyyyyyzz_1, tg_zzz_yyyyyzz_1, tg_zzz_yyyyyzzz_0, \
                                         tg_zzz_yyyyyzzz_1, tg_zzz_yyyyzzz_1, tg_zzz_yyyyzzzz_0, tg_zzz_yyyyzzzz_1, \
                                         tg_zzz_yyyzzzz_1, tg_zzz_yyyzzzzz_0, tg_zzz_yyyzzzzz_1, tg_zzz_yyzzzzz_1, \
                                         tg_zzz_yyzzzzzz_0, tg_zzz_yyzzzzzz_1, tg_zzz_yzzzzzz_1, tg_zzz_yzzzzzzz_0, \
                                         tg_zzz_yzzzzzzz_1, tg_zzz_zzzzzzz_1, tg_zzz_zzzzzzzz_0, tg_zzz_zzzzzzzz_1, \
                                         tg_zzzz_xxxxxxxx_0, tg_zzzz_xxxxxxxy_0, tg_zzzz_xxxxxxxz_0, tg_zzzz_xxxxxxyy_0, \
                                         tg_zzzz_xxxxxxyz_0, tg_zzzz_xxxxxxzz_0, tg_zzzz_xxxxxyyy_0, tg_zzzz_xxxxxyyz_0, \
                                         tg_zzzz_xxxxxyzz_0, tg_zzzz_xxxxxzzz_0, tg_zzzz_xxxxyyyy_0, tg_zzzz_xxxxyyyz_0, \
                                         tg_zzzz_xxxxyyzz_0, tg_zzzz_xxxxyzzz_0, tg_zzzz_xxxxzzzz_0, tg_zzzz_xxxyyyyy_0, \
                                         tg_zzzz_xxxyyyyz_0, tg_zzzz_xxxyyyzz_0, tg_zzzz_xxxyyzzz_0, tg_zzzz_xxxyzzzz_0, \
                                         tg_zzzz_xxxzzzzz_0, tg_zzzz_xxyyyyyy_0, tg_zzzz_xxyyyyyz_0, tg_zzzz_xxyyyyzz_0, \
                                         tg_zzzz_xxyyyzzz_0, tg_zzzz_xxyyzzzz_0, tg_zzzz_xxyzzzzz_0, tg_zzzz_xxzzzzzz_0, \
                                         tg_zzzz_xyyyyyyy_0, tg_zzzz_xyyyyyyz_0, tg_zzzz_xyyyyyzz_0, tg_zzzz_xyyyyzzz_0, \
                                         tg_zzzz_xyyyzzzz_0, tg_zzzz_xyyzzzzz_0, tg_zzzz_xyzzzzzz_0, tg_zzzz_xzzzzzzz_0, \
                                         tg_zzzz_yyyyyyyy_0, tg_zzzz_yyyyyyyz_0, tg_zzzz_yyyyyyzz_0, tg_zzzz_yyyyyzzz_0, \
                                         tg_zzzz_yyyyzzzz_0, tg_zzzz_yyyzzzzz_0, tg_zzzz_yyzzzzzz_0, tg_zzzz_yzzzzzzz_0, \
                                         tg_zzzz_zzzzzzzz_0, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_yzzz_yyzzzzzz_0[j] = pb_y * tg_zzz_yyzzzzzz_0[j] + wp_y[j] * tg_zzz_yyzzzzzz_1[j] + fl1_fxn * tg_zzz_yzzzzzz_1[j];

                    tg_yzzz_yzzzzzzz_0[j] = pb_y * tg_zzz_yzzzzzzz_0[j] + wp_y[j] * tg_zzz_yzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_zzzzzzz_1[j];

                    tg_yzzz_zzzzzzzz_0[j] = pb_y * tg_zzz_zzzzzzzz_0[j] + wp_y[j] * tg_zzz_zzzzzzzz_1[j];

                    tg_zzzz_xxxxxxxx_0[j] = pb_z * tg_zzz_xxxxxxxx_0[j] + wp_z[j] * tg_zzz_xxxxxxxx_1[j] + 1.5 * fl1_fx * tg_zz_xxxxxxxx_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxxxxxxx_1[j];

                    tg_zzzz_xxxxxxxy_0[j] = pb_z * tg_zzz_xxxxxxxy_0[j] + wp_z[j] * tg_zzz_xxxxxxxy_1[j] + 1.5 * fl1_fx * tg_zz_xxxxxxxy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxxxxxxy_1[j];

                    tg_zzzz_xxxxxxxz_0[j] = pb_z * tg_zzz_xxxxxxxz_0[j] + wp_z[j] * tg_zzz_xxxxxxxz_1[j] + 1.5 * fl1_fx * tg_zz_xxxxxxxz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxxxxxxz_1[j] + 0.5 * fl1_fxn * tg_zzz_xxxxxxx_1[j];

                    tg_zzzz_xxxxxxyy_0[j] = pb_z * tg_zzz_xxxxxxyy_0[j] + wp_z[j] * tg_zzz_xxxxxxyy_1[j] + 1.5 * fl1_fx * tg_zz_xxxxxxyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxxxxxyy_1[j];

                    tg_zzzz_xxxxxxyz_0[j] = pb_z * tg_zzz_xxxxxxyz_0[j] + wp_z[j] * tg_zzz_xxxxxxyz_1[j] + 1.5 * fl1_fx * tg_zz_xxxxxxyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxxxxxyz_1[j] + 0.5 * fl1_fxn * tg_zzz_xxxxxxy_1[j];

                    tg_zzzz_xxxxxxzz_0[j] = pb_z * tg_zzz_xxxxxxzz_0[j] + wp_z[j] * tg_zzz_xxxxxxzz_1[j] + 1.5 * fl1_fx * tg_zz_xxxxxxzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxxxxxzz_1[j] + fl1_fxn * tg_zzz_xxxxxxz_1[j];

                    tg_zzzz_xxxxxyyy_0[j] = pb_z * tg_zzz_xxxxxyyy_0[j] + wp_z[j] * tg_zzz_xxxxxyyy_1[j] + 1.5 * fl1_fx * tg_zz_xxxxxyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxxxxyyy_1[j];

                    tg_zzzz_xxxxxyyz_0[j] = pb_z * tg_zzz_xxxxxyyz_0[j] + wp_z[j] * tg_zzz_xxxxxyyz_1[j] + 1.5 * fl1_fx * tg_zz_xxxxxyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxxxxyyz_1[j] + 0.5 * fl1_fxn * tg_zzz_xxxxxyy_1[j];

                    tg_zzzz_xxxxxyzz_0[j] = pb_z * tg_zzz_xxxxxyzz_0[j] + wp_z[j] * tg_zzz_xxxxxyzz_1[j] + 1.5 * fl1_fx * tg_zz_xxxxxyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxxxxyzz_1[j] + fl1_fxn * tg_zzz_xxxxxyz_1[j];

                    tg_zzzz_xxxxxzzz_0[j] = pb_z * tg_zzz_xxxxxzzz_0[j] + wp_z[j] * tg_zzz_xxxxxzzz_1[j] + 1.5 * fl1_fx * tg_zz_xxxxxzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxxxxzzz_1[j] + 1.5 * fl1_fxn * tg_zzz_xxxxxzz_1[j];

                    tg_zzzz_xxxxyyyy_0[j] = pb_z * tg_zzz_xxxxyyyy_0[j] + wp_z[j] * tg_zzz_xxxxyyyy_1[j] + 1.5 * fl1_fx * tg_zz_xxxxyyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxxxyyyy_1[j];

                    tg_zzzz_xxxxyyyz_0[j] = pb_z * tg_zzz_xxxxyyyz_0[j] + wp_z[j] * tg_zzz_xxxxyyyz_1[j] + 1.5 * fl1_fx * tg_zz_xxxxyyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxxxyyyz_1[j] + 0.5 * fl1_fxn * tg_zzz_xxxxyyy_1[j];

                    tg_zzzz_xxxxyyzz_0[j] = pb_z * tg_zzz_xxxxyyzz_0[j] + wp_z[j] * tg_zzz_xxxxyyzz_1[j] + 1.5 * fl1_fx * tg_zz_xxxxyyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxxxyyzz_1[j] + fl1_fxn * tg_zzz_xxxxyyz_1[j];

                    tg_zzzz_xxxxyzzz_0[j] = pb_z * tg_zzz_xxxxyzzz_0[j] + wp_z[j] * tg_zzz_xxxxyzzz_1[j] + 1.5 * fl1_fx * tg_zz_xxxxyzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxxxyzzz_1[j] + 1.5 * fl1_fxn * tg_zzz_xxxxyzz_1[j];

                    tg_zzzz_xxxxzzzz_0[j] = pb_z * tg_zzz_xxxxzzzz_0[j] + wp_z[j] * tg_zzz_xxxxzzzz_1[j] + 1.5 * fl1_fx * tg_zz_xxxxzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxxxzzzz_1[j] + 2.0 * fl1_fxn * tg_zzz_xxxxzzz_1[j];

                    tg_zzzz_xxxyyyyy_0[j] = pb_z * tg_zzz_xxxyyyyy_0[j] + wp_z[j] * tg_zzz_xxxyyyyy_1[j] + 1.5 * fl1_fx * tg_zz_xxxyyyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxxyyyyy_1[j];

                    tg_zzzz_xxxyyyyz_0[j] = pb_z * tg_zzz_xxxyyyyz_0[j] + wp_z[j] * tg_zzz_xxxyyyyz_1[j] + 1.5 * fl1_fx * tg_zz_xxxyyyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxxyyyyz_1[j] + 0.5 * fl1_fxn * tg_zzz_xxxyyyy_1[j];

                    tg_zzzz_xxxyyyzz_0[j] = pb_z * tg_zzz_xxxyyyzz_0[j] + wp_z[j] * tg_zzz_xxxyyyzz_1[j] + 1.5 * fl1_fx * tg_zz_xxxyyyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxxyyyzz_1[j] + fl1_fxn * tg_zzz_xxxyyyz_1[j];

                    tg_zzzz_xxxyyzzz_0[j] = pb_z * tg_zzz_xxxyyzzz_0[j] + wp_z[j] * tg_zzz_xxxyyzzz_1[j] + 1.5 * fl1_fx * tg_zz_xxxyyzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxxyyzzz_1[j] + 1.5 * fl1_fxn * tg_zzz_xxxyyzz_1[j];

                    tg_zzzz_xxxyzzzz_0[j] = pb_z * tg_zzz_xxxyzzzz_0[j] + wp_z[j] * tg_zzz_xxxyzzzz_1[j] + 1.5 * fl1_fx * tg_zz_xxxyzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxxyzzzz_1[j] + 2.0 * fl1_fxn * tg_zzz_xxxyzzz_1[j];

                    tg_zzzz_xxxzzzzz_0[j] = pb_z * tg_zzz_xxxzzzzz_0[j] + wp_z[j] * tg_zzz_xxxzzzzz_1[j] + 1.5 * fl1_fx * tg_zz_xxxzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxxzzzzz_1[j] + 2.5 * fl1_fxn * tg_zzz_xxxzzzz_1[j];

                    tg_zzzz_xxyyyyyy_0[j] = pb_z * tg_zzz_xxyyyyyy_0[j] + wp_z[j] * tg_zzz_xxyyyyyy_1[j] + 1.5 * fl1_fx * tg_zz_xxyyyyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxyyyyyy_1[j];

                    tg_zzzz_xxyyyyyz_0[j] = pb_z * tg_zzz_xxyyyyyz_0[j] + wp_z[j] * tg_zzz_xxyyyyyz_1[j] + 1.5 * fl1_fx * tg_zz_xxyyyyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxyyyyyz_1[j] + 0.5 * fl1_fxn * tg_zzz_xxyyyyy_1[j];

                    tg_zzzz_xxyyyyzz_0[j] = pb_z * tg_zzz_xxyyyyzz_0[j] + wp_z[j] * tg_zzz_xxyyyyzz_1[j] + 1.5 * fl1_fx * tg_zz_xxyyyyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxyyyyzz_1[j] + fl1_fxn * tg_zzz_xxyyyyz_1[j];

                    tg_zzzz_xxyyyzzz_0[j] = pb_z * tg_zzz_xxyyyzzz_0[j] + wp_z[j] * tg_zzz_xxyyyzzz_1[j] + 1.5 * fl1_fx * tg_zz_xxyyyzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxyyyzzz_1[j] + 1.5 * fl1_fxn * tg_zzz_xxyyyzz_1[j];

                    tg_zzzz_xxyyzzzz_0[j] = pb_z * tg_zzz_xxyyzzzz_0[j] + wp_z[j] * tg_zzz_xxyyzzzz_1[j] + 1.5 * fl1_fx * tg_zz_xxyyzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxyyzzzz_1[j] + 2.0 * fl1_fxn * tg_zzz_xxyyzzz_1[j];

                    tg_zzzz_xxyzzzzz_0[j] = pb_z * tg_zzz_xxyzzzzz_0[j] + wp_z[j] * tg_zzz_xxyzzzzz_1[j] + 1.5 * fl1_fx * tg_zz_xxyzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxyzzzzz_1[j] + 2.5 * fl1_fxn * tg_zzz_xxyzzzz_1[j];

                    tg_zzzz_xxzzzzzz_0[j] = pb_z * tg_zzz_xxzzzzzz_0[j] + wp_z[j] * tg_zzz_xxzzzzzz_1[j] + 1.5 * fl1_fx * tg_zz_xxzzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxzzzzzz_1[j] + 3.0 * fl1_fxn * tg_zzz_xxzzzzz_1[j];

                    tg_zzzz_xyyyyyyy_0[j] = pb_z * tg_zzz_xyyyyyyy_0[j] + wp_z[j] * tg_zzz_xyyyyyyy_1[j] + 1.5 * fl1_fx * tg_zz_xyyyyyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xyyyyyyy_1[j];

                    tg_zzzz_xyyyyyyz_0[j] = pb_z * tg_zzz_xyyyyyyz_0[j] + wp_z[j] * tg_zzz_xyyyyyyz_1[j] + 1.5 * fl1_fx * tg_zz_xyyyyyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xyyyyyyz_1[j] + 0.5 * fl1_fxn * tg_zzz_xyyyyyy_1[j];

                    tg_zzzz_xyyyyyzz_0[j] = pb_z * tg_zzz_xyyyyyzz_0[j] + wp_z[j] * tg_zzz_xyyyyyzz_1[j] + 1.5 * fl1_fx * tg_zz_xyyyyyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xyyyyyzz_1[j] + fl1_fxn * tg_zzz_xyyyyyz_1[j];

                    tg_zzzz_xyyyyzzz_0[j] = pb_z * tg_zzz_xyyyyzzz_0[j] + wp_z[j] * tg_zzz_xyyyyzzz_1[j] + 1.5 * fl1_fx * tg_zz_xyyyyzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xyyyyzzz_1[j] + 1.5 * fl1_fxn * tg_zzz_xyyyyzz_1[j];

                    tg_zzzz_xyyyzzzz_0[j] = pb_z * tg_zzz_xyyyzzzz_0[j] + wp_z[j] * tg_zzz_xyyyzzzz_1[j] + 1.5 * fl1_fx * tg_zz_xyyyzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xyyyzzzz_1[j] + 2.0 * fl1_fxn * tg_zzz_xyyyzzz_1[j];

                    tg_zzzz_xyyzzzzz_0[j] = pb_z * tg_zzz_xyyzzzzz_0[j] + wp_z[j] * tg_zzz_xyyzzzzz_1[j] + 1.5 * fl1_fx * tg_zz_xyyzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xyyzzzzz_1[j] + 2.5 * fl1_fxn * tg_zzz_xyyzzzz_1[j];

                    tg_zzzz_xyzzzzzz_0[j] = pb_z * tg_zzz_xyzzzzzz_0[j] + wp_z[j] * tg_zzz_xyzzzzzz_1[j] + 1.5 * fl1_fx * tg_zz_xyzzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xyzzzzzz_1[j] + 3.0 * fl1_fxn * tg_zzz_xyzzzzz_1[j];

                    tg_zzzz_xzzzzzzz_0[j] = pb_z * tg_zzz_xzzzzzzz_0[j] + wp_z[j] * tg_zzz_xzzzzzzz_1[j] + 1.5 * fl1_fx * tg_zz_xzzzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xzzzzzzz_1[j] + 3.5 * fl1_fxn * tg_zzz_xzzzzzz_1[j];

                    tg_zzzz_yyyyyyyy_0[j] = pb_z * tg_zzz_yyyyyyyy_0[j] + wp_z[j] * tg_zzz_yyyyyyyy_1[j] + 1.5 * fl1_fx * tg_zz_yyyyyyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_yyyyyyyy_1[j];

                    tg_zzzz_yyyyyyyz_0[j] = pb_z * tg_zzz_yyyyyyyz_0[j] + wp_z[j] * tg_zzz_yyyyyyyz_1[j] + 1.5 * fl1_fx * tg_zz_yyyyyyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_yyyyyyyz_1[j] + 0.5 * fl1_fxn * tg_zzz_yyyyyyy_1[j];

                    tg_zzzz_yyyyyyzz_0[j] = pb_z * tg_zzz_yyyyyyzz_0[j] + wp_z[j] * tg_zzz_yyyyyyzz_1[j] + 1.5 * fl1_fx * tg_zz_yyyyyyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_yyyyyyzz_1[j] + fl1_fxn * tg_zzz_yyyyyyz_1[j];

                    tg_zzzz_yyyyyzzz_0[j] = pb_z * tg_zzz_yyyyyzzz_0[j] + wp_z[j] * tg_zzz_yyyyyzzz_1[j] + 1.5 * fl1_fx * tg_zz_yyyyyzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_yyyyyzzz_1[j] + 1.5 * fl1_fxn * tg_zzz_yyyyyzz_1[j];

                    tg_zzzz_yyyyzzzz_0[j] = pb_z * tg_zzz_yyyyzzzz_0[j] + wp_z[j] * tg_zzz_yyyyzzzz_1[j] + 1.5 * fl1_fx * tg_zz_yyyyzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_yyyyzzzz_1[j] + 2.0 * fl1_fxn * tg_zzz_yyyyzzz_1[j];

                    tg_zzzz_yyyzzzzz_0[j] = pb_z * tg_zzz_yyyzzzzz_0[j] + wp_z[j] * tg_zzz_yyyzzzzz_1[j] + 1.5 * fl1_fx * tg_zz_yyyzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_yyyzzzzz_1[j] + 2.5 * fl1_fxn * tg_zzz_yyyzzzz_1[j];

                    tg_zzzz_yyzzzzzz_0[j] = pb_z * tg_zzz_yyzzzzzz_0[j] + wp_z[j] * tg_zzz_yyzzzzzz_1[j] + 1.5 * fl1_fx * tg_zz_yyzzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_yyzzzzzz_1[j] + 3.0 * fl1_fxn * tg_zzz_yyzzzzz_1[j];

                    tg_zzzz_yzzzzzzz_0[j] = pb_z * tg_zzz_yzzzzzzz_0[j] + wp_z[j] * tg_zzz_yzzzzzzz_1[j] + 1.5 * fl1_fx * tg_zz_yzzzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_yzzzzzzz_1[j] + 3.5 * fl1_fxn * tg_zzz_yzzzzzz_1[j];

                    tg_zzzz_zzzzzzzz_0[j] = pb_z * tg_zzz_zzzzzzzz_0[j] + wp_z[j] * tg_zzz_zzzzzzzz_1[j] + 1.5 * fl1_fx * tg_zz_zzzzzzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_zzzzzzzz_1[j] + 4.0 * fl1_fxn * tg_zzz_zzzzzzz_1[j];
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

