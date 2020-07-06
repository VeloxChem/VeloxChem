//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ElectronRepulsionRecFuncForFI.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSFSI(      CMemBlock2D<double>* primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSFSI_0_94(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSFSI_94_187(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 

        erirecfunc::compElectronRepulsionForSFSI_187_280(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 
    }

    void
    compElectronRepulsionForSFSI_0_94(      CMemBlock2D<double>* primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& wpDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetPrimPairs,
                                      const int32_t              iContrPair)
    {
        // Batch of Integrals (0,94)

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
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_3_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_3_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_2_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_2_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_1_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {1, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_1_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {1, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_2_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {5, -1, -1, -1}, 
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

                auto tg_xx_xxxxxx_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx); 

                auto tg_xx_xxxxxy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 1); 

                auto tg_xx_xxxxxz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 2); 

                auto tg_xx_xxxxyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 3); 

                auto tg_xx_xxxxyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 4); 

                auto tg_xx_xxxxzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 5); 

                auto tg_xx_xxxyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 6); 

                auto tg_xx_xxxyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 7); 

                auto tg_xx_xxxyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 8); 

                auto tg_xx_xxxzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 9); 

                auto tg_xx_xxyyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 10); 

                auto tg_xx_xxyyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 11); 

                auto tg_xx_xxyyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 12); 

                auto tg_xx_xxyzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 13); 

                auto tg_xx_xxzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 14); 

                auto tg_xx_xyyyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 15); 

                auto tg_xx_xyyyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 16); 

                auto tg_xx_xyyyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 17); 

                auto tg_xx_xyyzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 18); 

                auto tg_xx_xyzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 19); 

                auto tg_xx_xzzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 20); 

                auto tg_xx_yyyyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 21); 

                auto tg_xx_yyyyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 22); 

                auto tg_xx_yyyyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 23); 

                auto tg_xx_yyyzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 24); 

                auto tg_xx_yyzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 25); 

                auto tg_xx_yzzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 26); 

                auto tg_xx_zzzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 27); 

                auto tg_xy_xxxxxx_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 28); 

                auto tg_xy_xxxxxy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 29); 

                auto tg_xy_xxxxxz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 30); 

                auto tg_xy_xxxxyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 31); 

                auto tg_xy_xxxxyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 32); 

                auto tg_xy_xxxxzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 33); 

                auto tg_xy_xxxyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 34); 

                auto tg_xy_xxxyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 35); 

                auto tg_xy_xxxyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 36); 

                auto tg_xy_xxxzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 37); 

                auto tg_xy_xxyyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 38); 

                auto tg_xy_xxyyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 39); 

                auto tg_xy_xxyyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 40); 

                auto tg_xy_xxyzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 41); 

                auto tg_xy_xxzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 42); 

                auto tg_xy_xyyyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 43); 

                auto tg_xy_xyyyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 44); 

                auto tg_xy_xyyyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 45); 

                auto tg_xy_xyyzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 46); 

                auto tg_xy_xyzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 47); 

                auto tg_xy_xzzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 48); 

                auto tg_xy_yyyyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 49); 

                auto tg_xy_yyyyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 50); 

                auto tg_xy_yyyyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 51); 

                auto tg_xy_yyyzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 52); 

                auto tg_xy_yyzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 53); 

                auto tg_xy_yzzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 54); 

                auto tg_xy_zzzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 55); 

                auto tg_xz_xxxxxx_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 56); 

                auto tg_xz_xxxxxy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 57); 

                auto tg_xz_xxxxxz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 58); 

                auto tg_xz_xxxxyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 59); 

                auto tg_xz_xxxxyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 60); 

                auto tg_xz_xxxxzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 61); 

                auto tg_xz_xxxyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 62); 

                auto tg_xz_xxxyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 63); 

                auto tg_xz_xxxyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 64); 

                auto tg_xz_xxxzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 65); 

                auto tg_xz_xxyyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 66); 

                auto tg_xz_xxyyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 67); 

                auto tg_xz_xxyyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 68); 

                auto tg_xz_xxyzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 69); 

                auto tg_xz_xxzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 70); 

                auto tg_xz_xyyyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 71); 

                auto tg_xz_xyyyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 72); 

                auto tg_xz_xyyyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 73); 

                auto tg_xz_xyyzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 74); 

                auto tg_xz_xyzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 75); 

                auto tg_xz_xzzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 76); 

                auto tg_xz_yyyyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 77); 

                auto tg_xz_yyyyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 78); 

                auto tg_xz_yyyyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 79); 

                auto tg_xz_yyyzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 80); 

                auto tg_xz_yyzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 81); 

                auto tg_xz_yzzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 82); 

                auto tg_xz_zzzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 83); 

                auto tg_yy_xxxxxx_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 84); 

                auto tg_yy_xxxxxy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 85); 

                auto tg_yy_xxxxxz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 86); 

                auto tg_yy_xxxxyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 87); 

                auto tg_yy_xxxxyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 88); 

                auto tg_yy_xxxxzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 89); 

                auto tg_yy_xxxyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 90); 

                auto tg_yy_xxxyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 91); 

                auto tg_yy_xxxyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 92); 

                auto tg_yy_xxxzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 93); 

                auto tg_xx_xxxxxx_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx); 

                auto tg_xx_xxxxxy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 1); 

                auto tg_xx_xxxxxz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 2); 

                auto tg_xx_xxxxyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 3); 

                auto tg_xx_xxxxyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 4); 

                auto tg_xx_xxxxzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 5); 

                auto tg_xx_xxxyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 6); 

                auto tg_xx_xxxyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 7); 

                auto tg_xx_xxxyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 8); 

                auto tg_xx_xxxzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 9); 

                auto tg_xx_xxyyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 10); 

                auto tg_xx_xxyyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 11); 

                auto tg_xx_xxyyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 12); 

                auto tg_xx_xxyzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 13); 

                auto tg_xx_xxzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 14); 

                auto tg_xx_xyyyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 15); 

                auto tg_xx_xyyyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 16); 

                auto tg_xx_xyyyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 17); 

                auto tg_xx_xyyzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 18); 

                auto tg_xx_xyzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 19); 

                auto tg_xx_xzzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 20); 

                auto tg_xx_yyyyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 21); 

                auto tg_xx_yyyyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 22); 

                auto tg_xx_yyyyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 23); 

                auto tg_xx_yyyzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 24); 

                auto tg_xx_yyzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 25); 

                auto tg_xx_yzzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 26); 

                auto tg_xx_zzzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 27); 

                auto tg_xy_xxxxxx_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 28); 

                auto tg_xy_xxxxxy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 29); 

                auto tg_xy_xxxxxz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 30); 

                auto tg_xy_xxxxyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 31); 

                auto tg_xy_xxxxyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 32); 

                auto tg_xy_xxxxzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 33); 

                auto tg_xy_xxxyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 34); 

                auto tg_xy_xxxyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 35); 

                auto tg_xy_xxxyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 36); 

                auto tg_xy_xxxzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 37); 

                auto tg_xy_xxyyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 38); 

                auto tg_xy_xxyyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 39); 

                auto tg_xy_xxyyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 40); 

                auto tg_xy_xxyzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 41); 

                auto tg_xy_xxzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 42); 

                auto tg_xy_xyyyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 43); 

                auto tg_xy_xyyyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 44); 

                auto tg_xy_xyyyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 45); 

                auto tg_xy_xyyzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 46); 

                auto tg_xy_xyzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 47); 

                auto tg_xy_xzzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 48); 

                auto tg_xy_yyyyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 49); 

                auto tg_xy_yyyyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 50); 

                auto tg_xy_yyyyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 51); 

                auto tg_xy_yyyzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 52); 

                auto tg_xy_yyzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 53); 

                auto tg_xy_yzzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 54); 

                auto tg_xy_zzzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 55); 

                auto tg_xz_xxxxxx_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 56); 

                auto tg_xz_xxxxxy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 57); 

                auto tg_xz_xxxxxz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 58); 

                auto tg_xz_xxxxyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 59); 

                auto tg_xz_xxxxyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 60); 

                auto tg_xz_xxxxzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 61); 

                auto tg_xz_xxxyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 62); 

                auto tg_xz_xxxyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 63); 

                auto tg_xz_xxxyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 64); 

                auto tg_xz_xxxzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 65); 

                auto tg_xz_xxyyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 66); 

                auto tg_xz_xxyyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 67); 

                auto tg_xz_xxyyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 68); 

                auto tg_xz_xxyzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 69); 

                auto tg_xz_xxzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 70); 

                auto tg_xz_xyyyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 71); 

                auto tg_xz_xyyyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 72); 

                auto tg_xz_xyyyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 73); 

                auto tg_xz_xyyzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 74); 

                auto tg_xz_xyzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 75); 

                auto tg_xz_xzzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 76); 

                auto tg_xz_yyyyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 77); 

                auto tg_xz_yyyyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 78); 

                auto tg_xz_yyyyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 79); 

                auto tg_xz_yyyzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 80); 

                auto tg_xz_yyzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 81); 

                auto tg_xz_yzzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 82); 

                auto tg_xz_zzzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 83); 

                auto tg_yy_xxxxxx_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 84); 

                auto tg_yy_xxxxxy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 85); 

                auto tg_yy_xxxxxz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 86); 

                auto tg_yy_xxxxyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 87); 

                auto tg_yy_xxxxyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 88); 

                auto tg_yy_xxxxzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 89); 

                auto tg_yy_xxxyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 90); 

                auto tg_yy_xxxyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 91); 

                auto tg_yy_xxxyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 92); 

                auto tg_yy_xxxzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 93); 

                auto tg_x_xxxxxx_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx); 

                auto tg_x_xxxxxy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 1); 

                auto tg_x_xxxxxz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 2); 

                auto tg_x_xxxxyy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 3); 

                auto tg_x_xxxxyz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 4); 

                auto tg_x_xxxxzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 5); 

                auto tg_x_xxxyyy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 6); 

                auto tg_x_xxxyyz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 7); 

                auto tg_x_xxxyzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 8); 

                auto tg_x_xxxzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 9); 

                auto tg_x_xxyyyy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 10); 

                auto tg_x_xxyyyz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 11); 

                auto tg_x_xxyyzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 12); 

                auto tg_x_xxyzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 13); 

                auto tg_x_xxzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 14); 

                auto tg_x_xyyyyy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 15); 

                auto tg_x_xyyyyz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 16); 

                auto tg_x_xyyyzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 17); 

                auto tg_x_xyyzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 18); 

                auto tg_x_xyzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 19); 

                auto tg_x_xzzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 20); 

                auto tg_x_yyyyyy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 21); 

                auto tg_x_yyyyyz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 22); 

                auto tg_x_yyyyzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 23); 

                auto tg_x_yyyzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 24); 

                auto tg_x_yyzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 25); 

                auto tg_x_yzzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 26); 

                auto tg_x_zzzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 27); 

                auto tg_y_xxxxxx_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 28); 

                auto tg_y_xxxxxy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 29); 

                auto tg_y_xxxxxz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 30); 

                auto tg_y_xxxxyy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 31); 

                auto tg_y_xxxxyz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 32); 

                auto tg_y_xxxxzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 33); 

                auto tg_y_xxxyyy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 34); 

                auto tg_y_xxxyyz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 35); 

                auto tg_y_xxxyzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 36); 

                auto tg_y_xxxzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 37); 

                auto tg_y_xxyyyy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 38); 

                auto tg_y_xxyyyz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 39); 

                auto tg_y_xxyyzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 40); 

                auto tg_y_xxyzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 41); 

                auto tg_y_xxzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 42); 

                auto tg_y_xyyyyy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 43); 

                auto tg_y_xyyyyz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 44); 

                auto tg_y_xyyyzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 45); 

                auto tg_y_xyyzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 46); 

                auto tg_y_xyzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 47); 

                auto tg_y_xzzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 48); 

                auto tg_y_yyyyyy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 49); 

                auto tg_y_yyyyyz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 50); 

                auto tg_y_yyyyzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 51); 

                auto tg_y_yyyzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 52); 

                auto tg_y_yyzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 53); 

                auto tg_y_yzzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 54); 

                auto tg_y_zzzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 55); 

                auto tg_z_xxxxxx_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 56); 

                auto tg_z_xxxxxy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 57); 

                auto tg_z_xxxxxz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 58); 

                auto tg_z_xxxxyy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 59); 

                auto tg_z_xxxxyz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 60); 

                auto tg_z_xxxxzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 61); 

                auto tg_z_xxxyyy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 62); 

                auto tg_z_xxxyyz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 63); 

                auto tg_z_xxxyzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 64); 

                auto tg_z_xxxzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 65); 

                auto tg_z_xxyyyy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 66); 

                auto tg_z_xxyyyz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 67); 

                auto tg_z_xxyyzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 68); 

                auto tg_z_xxyzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 69); 

                auto tg_z_xxzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 70); 

                auto tg_z_xyyyyy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 71); 

                auto tg_z_xyyyyz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 72); 

                auto tg_z_xyyyzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 73); 

                auto tg_z_xyyzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 74); 

                auto tg_z_xyzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 75); 

                auto tg_z_xzzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 76); 

                auto tg_z_yyyyyy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 77); 

                auto tg_z_yyyyyz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 78); 

                auto tg_z_yyyyzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 79); 

                auto tg_z_yyyzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 80); 

                auto tg_z_yyzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 81); 

                auto tg_z_yzzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 82); 

                auto tg_z_zzzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 83); 

                auto tg_x_xxxxxx_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx); 

                auto tg_x_xxxxxy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 1); 

                auto tg_x_xxxxxz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 2); 

                auto tg_x_xxxxyy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 3); 

                auto tg_x_xxxxyz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 4); 

                auto tg_x_xxxxzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 5); 

                auto tg_x_xxxyyy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 6); 

                auto tg_x_xxxyyz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 7); 

                auto tg_x_xxxyzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 8); 

                auto tg_x_xxxzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 9); 

                auto tg_x_xxyyyy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 10); 

                auto tg_x_xxyyyz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 11); 

                auto tg_x_xxyyzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 12); 

                auto tg_x_xxyzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 13); 

                auto tg_x_xxzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 14); 

                auto tg_x_xyyyyy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 15); 

                auto tg_x_xyyyyz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 16); 

                auto tg_x_xyyyzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 17); 

                auto tg_x_xyyzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 18); 

                auto tg_x_xyzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 19); 

                auto tg_x_xzzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 20); 

                auto tg_x_yyyyyy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 21); 

                auto tg_x_yyyyyz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 22); 

                auto tg_x_yyyyzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 23); 

                auto tg_x_yyyzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 24); 

                auto tg_x_yyzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 25); 

                auto tg_x_yzzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 26); 

                auto tg_x_zzzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 27); 

                auto tg_y_xxxxxx_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 28); 

                auto tg_y_xxxxxy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 29); 

                auto tg_y_xxxxxz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 30); 

                auto tg_y_xxxxyy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 31); 

                auto tg_y_xxxxyz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 32); 

                auto tg_y_xxxxzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 33); 

                auto tg_y_xxxyyy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 34); 

                auto tg_y_xxxyyz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 35); 

                auto tg_y_xxxyzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 36); 

                auto tg_y_xxxzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 37); 

                auto tg_y_xxyyyy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 38); 

                auto tg_y_xxyyyz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 39); 

                auto tg_y_xxyyzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 40); 

                auto tg_y_xxyzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 41); 

                auto tg_y_xxzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 42); 

                auto tg_y_xyyyyy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 43); 

                auto tg_y_xyyyyz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 44); 

                auto tg_y_xyyyzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 45); 

                auto tg_y_xyyzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 46); 

                auto tg_y_xyzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 47); 

                auto tg_y_xzzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 48); 

                auto tg_y_yyyyyy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 49); 

                auto tg_y_yyyyyz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 50); 

                auto tg_y_yyyyzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 51); 

                auto tg_y_yyyzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 52); 

                auto tg_y_yyzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 53); 

                auto tg_y_yzzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 54); 

                auto tg_y_zzzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 55); 

                auto tg_z_xxxxxx_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 56); 

                auto tg_z_xxxxxy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 57); 

                auto tg_z_xxxxxz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 58); 

                auto tg_z_xxxxyy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 59); 

                auto tg_z_xxxxyz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 60); 

                auto tg_z_xxxxzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 61); 

                auto tg_z_xxxyyy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 62); 

                auto tg_z_xxxyyz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 63); 

                auto tg_z_xxxyzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 64); 

                auto tg_z_xxxzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 65); 

                auto tg_z_xxyyyy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 66); 

                auto tg_z_xxyyyz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 67); 

                auto tg_z_xxyyzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 68); 

                auto tg_z_xxyzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 69); 

                auto tg_z_xxzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 70); 

                auto tg_z_xyyyyy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 71); 

                auto tg_z_xyyyyz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 72); 

                auto tg_z_xyyyzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 73); 

                auto tg_z_xyyzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 74); 

                auto tg_z_xyzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 75); 

                auto tg_z_xzzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 76); 

                auto tg_z_yyyyyy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 77); 

                auto tg_z_yyyyyz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 78); 

                auto tg_z_yyyyzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 79); 

                auto tg_z_yyyzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 80); 

                auto tg_z_yyzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 81); 

                auto tg_z_yzzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 82); 

                auto tg_z_zzzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 83); 

                auto tg_xx_xxxxx_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx); 

                auto tg_xx_xxxxy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 1); 

                auto tg_xx_xxxxz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 2); 

                auto tg_xx_xxxyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 3); 

                auto tg_xx_xxxyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 4); 

                auto tg_xx_xxxzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 5); 

                auto tg_xx_xxyyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 6); 

                auto tg_xx_xxyyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 7); 

                auto tg_xx_xxyzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 8); 

                auto tg_xx_xxzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 9); 

                auto tg_xx_xyyyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 10); 

                auto tg_xx_xyyyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 11); 

                auto tg_xx_xyyzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 12); 

                auto tg_xx_xyzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 13); 

                auto tg_xx_xzzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 14); 

                auto tg_xx_yyyyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 15); 

                auto tg_xx_yyyyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 16); 

                auto tg_xx_yyyzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 17); 

                auto tg_xx_yyzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 18); 

                auto tg_xx_yzzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 19); 

                auto tg_xx_zzzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 20); 

                auto tg_xy_xxxxx_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 21); 

                auto tg_xy_xxxxy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 22); 

                auto tg_xy_xxxxz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 23); 

                auto tg_xy_xxxyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 24); 

                auto tg_xy_xxxyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 25); 

                auto tg_xy_xxxzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 26); 

                auto tg_xy_xxyyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 27); 

                auto tg_xy_xxyyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 28); 

                auto tg_xy_xxyzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 29); 

                auto tg_xy_xxzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 30); 

                auto tg_xy_xyyyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 31); 

                auto tg_xy_xyyyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 32); 

                auto tg_xy_xyyzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 33); 

                auto tg_xy_xyzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 34); 

                auto tg_xy_xzzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 35); 

                auto tg_xy_yyyyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 36); 

                auto tg_xy_yyyyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 37); 

                auto tg_xy_yyyzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 38); 

                auto tg_xy_yyzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 39); 

                auto tg_xy_yzzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 40); 

                auto tg_xy_zzzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 41); 

                auto tg_xz_xxxxx_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 42); 

                auto tg_xz_xxxxy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 43); 

                auto tg_xz_xxxxz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 44); 

                auto tg_xz_xxxyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 45); 

                auto tg_xz_xxxyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 46); 

                auto tg_xz_xxxzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 47); 

                auto tg_xz_xxyyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 48); 

                auto tg_xz_xxyyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 49); 

                auto tg_xz_xxyzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 50); 

                auto tg_xz_xxzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 51); 

                auto tg_xz_xyyyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 52); 

                auto tg_xz_xyyyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 53); 

                auto tg_xz_xyyzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 54); 

                auto tg_xz_xyzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 55); 

                auto tg_xz_xzzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 56); 

                auto tg_xz_yyyyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 57); 

                auto tg_xz_yyyyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 58); 

                auto tg_xz_yyyzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 59); 

                auto tg_xz_yyzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 60); 

                auto tg_xz_yzzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 61); 

                auto tg_xz_zzzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 62); 

                auto tg_yy_xxxxx_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 63); 

                auto tg_yy_xxxxy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 64); 

                auto tg_yy_xxxxz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 65); 

                auto tg_yy_xxxyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 66); 

                auto tg_yy_xxxyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 67); 

                auto tg_yy_xxxzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 68); 

                auto tg_yy_xxyyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 69); 

                auto tg_yy_xxyyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 70); 

                auto tg_yy_xxyzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 71); 

                auto tg_yy_xxzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 72); 

                // set up pointers to integrals

                auto tg_xxx_xxxxxx_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx); 

                auto tg_xxx_xxxxxy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 1); 

                auto tg_xxx_xxxxxz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 2); 

                auto tg_xxx_xxxxyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 3); 

                auto tg_xxx_xxxxyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 4); 

                auto tg_xxx_xxxxzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 5); 

                auto tg_xxx_xxxyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 6); 

                auto tg_xxx_xxxyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 7); 

                auto tg_xxx_xxxyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 8); 

                auto tg_xxx_xxxzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 9); 

                auto tg_xxx_xxyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 10); 

                auto tg_xxx_xxyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 11); 

                auto tg_xxx_xxyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 12); 

                auto tg_xxx_xxyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 13); 

                auto tg_xxx_xxzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 14); 

                auto tg_xxx_xyyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 15); 

                auto tg_xxx_xyyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 16); 

                auto tg_xxx_xyyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 17); 

                auto tg_xxx_xyyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 18); 

                auto tg_xxx_xyzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 19); 

                auto tg_xxx_xzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 20); 

                auto tg_xxx_yyyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 21); 

                auto tg_xxx_yyyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 22); 

                auto tg_xxx_yyyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 23); 

                auto tg_xxx_yyyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 24); 

                auto tg_xxx_yyzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 25); 

                auto tg_xxx_yzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 26); 

                auto tg_xxx_zzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 27); 

                auto tg_xxy_xxxxxx_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 28); 

                auto tg_xxy_xxxxxy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 29); 

                auto tg_xxy_xxxxxz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 30); 

                auto tg_xxy_xxxxyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 31); 

                auto tg_xxy_xxxxyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 32); 

                auto tg_xxy_xxxxzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 33); 

                auto tg_xxy_xxxyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 34); 

                auto tg_xxy_xxxyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 35); 

                auto tg_xxy_xxxyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 36); 

                auto tg_xxy_xxxzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 37); 

                auto tg_xxy_xxyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 38); 

                auto tg_xxy_xxyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 39); 

                auto tg_xxy_xxyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 40); 

                auto tg_xxy_xxyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 41); 

                auto tg_xxy_xxzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 42); 

                auto tg_xxy_xyyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 43); 

                auto tg_xxy_xyyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 44); 

                auto tg_xxy_xyyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 45); 

                auto tg_xxy_xyyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 46); 

                auto tg_xxy_xyzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 47); 

                auto tg_xxy_xzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 48); 

                auto tg_xxy_yyyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 49); 

                auto tg_xxy_yyyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 50); 

                auto tg_xxy_yyyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 51); 

                auto tg_xxy_yyyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 52); 

                auto tg_xxy_yyzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 53); 

                auto tg_xxy_yzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 54); 

                auto tg_xxy_zzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 55); 

                auto tg_xxz_xxxxxx_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 56); 

                auto tg_xxz_xxxxxy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 57); 

                auto tg_xxz_xxxxxz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 58); 

                auto tg_xxz_xxxxyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 59); 

                auto tg_xxz_xxxxyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 60); 

                auto tg_xxz_xxxxzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 61); 

                auto tg_xxz_xxxyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 62); 

                auto tg_xxz_xxxyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 63); 

                auto tg_xxz_xxxyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 64); 

                auto tg_xxz_xxxzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 65); 

                auto tg_xxz_xxyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 66); 

                auto tg_xxz_xxyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 67); 

                auto tg_xxz_xxyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 68); 

                auto tg_xxz_xxyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 69); 

                auto tg_xxz_xxzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 70); 

                auto tg_xxz_xyyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 71); 

                auto tg_xxz_xyyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 72); 

                auto tg_xxz_xyyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 73); 

                auto tg_xxz_xyyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 74); 

                auto tg_xxz_xyzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 75); 

                auto tg_xxz_xzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 76); 

                auto tg_xxz_yyyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 77); 

                auto tg_xxz_yyyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 78); 

                auto tg_xxz_yyyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 79); 

                auto tg_xxz_yyyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 80); 

                auto tg_xxz_yyzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 81); 

                auto tg_xxz_yzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 82); 

                auto tg_xxz_zzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 83); 

                auto tg_xyy_xxxxxx_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 84); 

                auto tg_xyy_xxxxxy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 85); 

                auto tg_xyy_xxxxxz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 86); 

                auto tg_xyy_xxxxyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 87); 

                auto tg_xyy_xxxxyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 88); 

                auto tg_xyy_xxxxzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 89); 

                auto tg_xyy_xxxyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 90); 

                auto tg_xyy_xxxyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 91); 

                auto tg_xyy_xxxyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 92); 

                auto tg_xyy_xxxzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 93); 

                // Batch of Integrals (0,94)

                #pragma omp simd aligned(fxn, fza, tg_x_xxxxxx_0, tg_x_xxxxxx_1, tg_x_xxxxxy_0, tg_x_xxxxxy_1, \
                                         tg_x_xxxxxz_0, tg_x_xxxxxz_1, tg_x_xxxxyy_0, tg_x_xxxxyy_1, tg_x_xxxxyz_0, \
                                         tg_x_xxxxyz_1, tg_x_xxxxzz_0, tg_x_xxxxzz_1, tg_x_xxxyyy_0, tg_x_xxxyyy_1, \
                                         tg_x_xxxyyz_0, tg_x_xxxyyz_1, tg_x_xxxyzz_0, tg_x_xxxyzz_1, tg_x_xxxzzz_0, \
                                         tg_x_xxxzzz_1, tg_x_xxyyyy_0, tg_x_xxyyyy_1, tg_x_xxyyyz_0, tg_x_xxyyyz_1, \
                                         tg_x_xxyyzz_0, tg_x_xxyyzz_1, tg_x_xxyzzz_0, tg_x_xxyzzz_1, tg_x_xxzzzz_0, \
                                         tg_x_xxzzzz_1, tg_x_xyyyyy_0, tg_x_xyyyyy_1, tg_x_xyyyyz_0, tg_x_xyyyyz_1, \
                                         tg_x_xyyyzz_0, tg_x_xyyyzz_1, tg_x_xyyzzz_0, tg_x_xyyzzz_1, tg_x_xyzzzz_0, \
                                         tg_x_xyzzzz_1, tg_x_xzzzzz_0, tg_x_xzzzzz_1, tg_x_yyyyyy_0, tg_x_yyyyyy_1, \
                                         tg_x_yyyyyz_0, tg_x_yyyyyz_1, tg_x_yyyyzz_0, tg_x_yyyyzz_1, tg_x_yyyzzz_0, \
                                         tg_x_yyyzzz_1, tg_x_yyzzzz_0, tg_x_yyzzzz_1, tg_x_yzzzzz_0, tg_x_yzzzzz_1, \
                                         tg_x_zzzzzz_0, tg_x_zzzzzz_1, tg_xx_xxxxx_1, tg_xx_xxxxxx_0, tg_xx_xxxxxx_1, \
                                         tg_xx_xxxxxy_0, tg_xx_xxxxxy_1, tg_xx_xxxxxz_0, tg_xx_xxxxxz_1, tg_xx_xxxxy_1, \
                                         tg_xx_xxxxyy_0, tg_xx_xxxxyy_1, tg_xx_xxxxyz_0, tg_xx_xxxxyz_1, tg_xx_xxxxz_1, \
                                         tg_xx_xxxxzz_0, tg_xx_xxxxzz_1, tg_xx_xxxyy_1, tg_xx_xxxyyy_0, tg_xx_xxxyyy_1, \
                                         tg_xx_xxxyyz_0, tg_xx_xxxyyz_1, tg_xx_xxxyz_1, tg_xx_xxxyzz_0, tg_xx_xxxyzz_1, \
                                         tg_xx_xxxzz_1, tg_xx_xxxzzz_0, tg_xx_xxxzzz_1, tg_xx_xxyyy_1, tg_xx_xxyyyy_0, \
                                         tg_xx_xxyyyy_1, tg_xx_xxyyyz_0, tg_xx_xxyyyz_1, tg_xx_xxyyz_1, tg_xx_xxyyzz_0, \
                                         tg_xx_xxyyzz_1, tg_xx_xxyzz_1, tg_xx_xxyzzz_0, tg_xx_xxyzzz_1, tg_xx_xxzzz_1, \
                                         tg_xx_xxzzzz_0, tg_xx_xxzzzz_1, tg_xx_xyyyy_1, tg_xx_xyyyyy_0, tg_xx_xyyyyy_1, \
                                         tg_xx_xyyyyz_0, tg_xx_xyyyyz_1, tg_xx_xyyyz_1, tg_xx_xyyyzz_0, tg_xx_xyyyzz_1, \
                                         tg_xx_xyyzz_1, tg_xx_xyyzzz_0, tg_xx_xyyzzz_1, tg_xx_xyzzz_1, tg_xx_xyzzzz_0, \
                                         tg_xx_xyzzzz_1, tg_xx_xzzzz_1, tg_xx_xzzzzz_0, tg_xx_xzzzzz_1, tg_xx_yyyyy_1, \
                                         tg_xx_yyyyyy_0, tg_xx_yyyyyy_1, tg_xx_yyyyyz_0, tg_xx_yyyyyz_1, tg_xx_yyyyz_1, \
                                         tg_xx_yyyyzz_0, tg_xx_yyyyzz_1, tg_xx_yyyzz_1, tg_xx_yyyzzz_0, tg_xx_yyyzzz_1, \
                                         tg_xx_yyzzz_1, tg_xx_yyzzzz_0, tg_xx_yyzzzz_1, tg_xx_yzzzz_1, tg_xx_yzzzzz_0, \
                                         tg_xx_yzzzzz_1, tg_xx_zzzzz_1, tg_xx_zzzzzz_0, tg_xx_zzzzzz_1, tg_xxx_xxxxxx_0, \
                                         tg_xxx_xxxxxy_0, tg_xxx_xxxxxz_0, tg_xxx_xxxxyy_0, tg_xxx_xxxxyz_0, tg_xxx_xxxxzz_0, \
                                         tg_xxx_xxxyyy_0, tg_xxx_xxxyyz_0, tg_xxx_xxxyzz_0, tg_xxx_xxxzzz_0, tg_xxx_xxyyyy_0, \
                                         tg_xxx_xxyyyz_0, tg_xxx_xxyyzz_0, tg_xxx_xxyzzz_0, tg_xxx_xxzzzz_0, tg_xxx_xyyyyy_0, \
                                         tg_xxx_xyyyyz_0, tg_xxx_xyyyzz_0, tg_xxx_xyyzzz_0, tg_xxx_xyzzzz_0, tg_xxx_xzzzzz_0, \
                                         tg_xxx_yyyyyy_0, tg_xxx_yyyyyz_0, tg_xxx_yyyyzz_0, tg_xxx_yyyzzz_0, tg_xxx_yyzzzz_0, \
                                         tg_xxx_yzzzzz_0, tg_xxx_zzzzzz_0, tg_xxy_xxxxxx_0, tg_xxy_xxxxxy_0, tg_xxy_xxxxxz_0, \
                                         tg_xxy_xxxxyy_0, tg_xxy_xxxxyz_0, tg_xxy_xxxxzz_0, tg_xxy_xxxyyy_0, tg_xxy_xxxyyz_0, \
                                         tg_xxy_xxxyzz_0, tg_xxy_xxxzzz_0, tg_xxy_xxyyyy_0, tg_xxy_xxyyyz_0, tg_xxy_xxyyzz_0, \
                                         tg_xxy_xxyzzz_0, tg_xxy_xxzzzz_0, tg_xxy_xyyyyy_0, tg_xxy_xyyyyz_0, tg_xxy_xyyyzz_0, \
                                         tg_xxy_xyyzzz_0, tg_xxy_xyzzzz_0, tg_xxy_xzzzzz_0, tg_xxy_yyyyyy_0, tg_xxy_yyyyyz_0, \
                                         tg_xxy_yyyyzz_0, tg_xxy_yyyzzz_0, tg_xxy_yyzzzz_0, tg_xxy_yzzzzz_0, tg_xxy_zzzzzz_0, \
                                         tg_xxz_xxxxxx_0, tg_xxz_xxxxxy_0, tg_xxz_xxxxxz_0, tg_xxz_xxxxyy_0, tg_xxz_xxxxyz_0, \
                                         tg_xxz_xxxxzz_0, tg_xxz_xxxyyy_0, tg_xxz_xxxyyz_0, tg_xxz_xxxyzz_0, tg_xxz_xxxzzz_0, \
                                         tg_xxz_xxyyyy_0, tg_xxz_xxyyyz_0, tg_xxz_xxyyzz_0, tg_xxz_xxyzzz_0, tg_xxz_xxzzzz_0, \
                                         tg_xxz_xyyyyy_0, tg_xxz_xyyyyz_0, tg_xxz_xyyyzz_0, tg_xxz_xyyzzz_0, tg_xxz_xyzzzz_0, \
                                         tg_xxz_xzzzzz_0, tg_xxz_yyyyyy_0, tg_xxz_yyyyyz_0, tg_xxz_yyyyzz_0, tg_xxz_yyyzzz_0, \
                                         tg_xxz_yyzzzz_0, tg_xxz_yzzzzz_0, tg_xxz_zzzzzz_0, tg_xy_xxxxx_1, tg_xy_xxxxxx_0, \
                                         tg_xy_xxxxxx_1, tg_xy_xxxxxy_0, tg_xy_xxxxxy_1, tg_xy_xxxxxz_0, tg_xy_xxxxxz_1, \
                                         tg_xy_xxxxy_1, tg_xy_xxxxyy_0, tg_xy_xxxxyy_1, tg_xy_xxxxyz_0, tg_xy_xxxxyz_1, \
                                         tg_xy_xxxxz_1, tg_xy_xxxxzz_0, tg_xy_xxxxzz_1, tg_xy_xxxyy_1, tg_xy_xxxyyy_0, \
                                         tg_xy_xxxyyy_1, tg_xy_xxxyyz_0, tg_xy_xxxyyz_1, tg_xy_xxxyz_1, tg_xy_xxxyzz_0, \
                                         tg_xy_xxxyzz_1, tg_xy_xxxzz_1, tg_xy_xxxzzz_0, tg_xy_xxxzzz_1, tg_xy_xxyyy_1, \
                                         tg_xy_xxyyyy_0, tg_xy_xxyyyy_1, tg_xy_xxyyyz_0, tg_xy_xxyyyz_1, tg_xy_xxyyz_1, \
                                         tg_xy_xxyyzz_0, tg_xy_xxyyzz_1, tg_xy_xxyzz_1, tg_xy_xxyzzz_0, tg_xy_xxyzzz_1, \
                                         tg_xy_xxzzz_1, tg_xy_xxzzzz_0, tg_xy_xxzzzz_1, tg_xy_xyyyy_1, tg_xy_xyyyyy_0, \
                                         tg_xy_xyyyyy_1, tg_xy_xyyyyz_0, tg_xy_xyyyyz_1, tg_xy_xyyyz_1, tg_xy_xyyyzz_0, \
                                         tg_xy_xyyyzz_1, tg_xy_xyyzz_1, tg_xy_xyyzzz_0, tg_xy_xyyzzz_1, tg_xy_xyzzz_1, \
                                         tg_xy_xyzzzz_0, tg_xy_xyzzzz_1, tg_xy_xzzzz_1, tg_xy_xzzzzz_0, tg_xy_xzzzzz_1, \
                                         tg_xy_yyyyy_1, tg_xy_yyyyyy_0, tg_xy_yyyyyy_1, tg_xy_yyyyyz_0, tg_xy_yyyyyz_1, \
                                         tg_xy_yyyyz_1, tg_xy_yyyyzz_0, tg_xy_yyyyzz_1, tg_xy_yyyzz_1, tg_xy_yyyzzz_0, \
                                         tg_xy_yyyzzz_1, tg_xy_yyzzz_1, tg_xy_yyzzzz_0, tg_xy_yyzzzz_1, tg_xy_yzzzz_1, \
                                         tg_xy_yzzzzz_0, tg_xy_yzzzzz_1, tg_xy_zzzzz_1, tg_xy_zzzzzz_0, tg_xy_zzzzzz_1, \
                                         tg_xyy_xxxxxx_0, tg_xyy_xxxxxy_0, tg_xyy_xxxxxz_0, tg_xyy_xxxxyy_0, tg_xyy_xxxxyz_0, \
                                         tg_xyy_xxxxzz_0, tg_xyy_xxxyyy_0, tg_xyy_xxxyyz_0, tg_xyy_xxxyzz_0, tg_xyy_xxxzzz_0, \
                                         tg_xz_xxxxx_1, tg_xz_xxxxxx_0, tg_xz_xxxxxx_1, tg_xz_xxxxxy_0, tg_xz_xxxxxy_1, \
                                         tg_xz_xxxxxz_0, tg_xz_xxxxxz_1, tg_xz_xxxxy_1, tg_xz_xxxxyy_0, tg_xz_xxxxyy_1, \
                                         tg_xz_xxxxyz_0, tg_xz_xxxxyz_1, tg_xz_xxxxz_1, tg_xz_xxxxzz_0, tg_xz_xxxxzz_1, \
                                         tg_xz_xxxyy_1, tg_xz_xxxyyy_0, tg_xz_xxxyyy_1, tg_xz_xxxyyz_0, tg_xz_xxxyyz_1, \
                                         tg_xz_xxxyz_1, tg_xz_xxxyzz_0, tg_xz_xxxyzz_1, tg_xz_xxxzz_1, tg_xz_xxxzzz_0, \
                                         tg_xz_xxxzzz_1, tg_xz_xxyyy_1, tg_xz_xxyyyy_0, tg_xz_xxyyyy_1, tg_xz_xxyyyz_0, \
                                         tg_xz_xxyyyz_1, tg_xz_xxyyz_1, tg_xz_xxyyzz_0, tg_xz_xxyyzz_1, tg_xz_xxyzz_1, \
                                         tg_xz_xxyzzz_0, tg_xz_xxyzzz_1, tg_xz_xxzzz_1, tg_xz_xxzzzz_0, tg_xz_xxzzzz_1, \
                                         tg_xz_xyyyy_1, tg_xz_xyyyyy_0, tg_xz_xyyyyy_1, tg_xz_xyyyyz_0, tg_xz_xyyyyz_1, \
                                         tg_xz_xyyyz_1, tg_xz_xyyyzz_0, tg_xz_xyyyzz_1, tg_xz_xyyzz_1, tg_xz_xyyzzz_0, \
                                         tg_xz_xyyzzz_1, tg_xz_xyzzz_1, tg_xz_xyzzzz_0, tg_xz_xyzzzz_1, tg_xz_xzzzz_1, \
                                         tg_xz_xzzzzz_0, tg_xz_xzzzzz_1, tg_xz_yyyyy_1, tg_xz_yyyyyy_0, tg_xz_yyyyyy_1, \
                                         tg_xz_yyyyyz_0, tg_xz_yyyyyz_1, tg_xz_yyyyz_1, tg_xz_yyyyzz_0, tg_xz_yyyyzz_1, \
                                         tg_xz_yyyzz_1, tg_xz_yyyzzz_0, tg_xz_yyyzzz_1, tg_xz_yyzzz_1, tg_xz_yyzzzz_0, \
                                         tg_xz_yyzzzz_1, tg_xz_yzzzz_1, tg_xz_yzzzzz_0, tg_xz_yzzzzz_1, tg_xz_zzzzz_1, \
                                         tg_xz_zzzzzz_0, tg_xz_zzzzzz_1, tg_y_xxxxxx_0, tg_y_xxxxxx_1, tg_y_xxxxxy_0, \
                                         tg_y_xxxxxy_1, tg_y_xxxxxz_0, tg_y_xxxxxz_1, tg_y_xxxxyy_0, tg_y_xxxxyy_1, \
                                         tg_y_xxxxyz_0, tg_y_xxxxyz_1, tg_y_xxxxzz_0, tg_y_xxxxzz_1, tg_y_xxxyyy_0, \
                                         tg_y_xxxyyy_1, tg_y_xxxyyz_0, tg_y_xxxyyz_1, tg_y_xxxyzz_0, tg_y_xxxyzz_1, \
                                         tg_y_xxxzzz_0, tg_y_xxxzzz_1, tg_y_xxyyyy_0, tg_y_xxyyyy_1, tg_y_xxyyyz_0, \
                                         tg_y_xxyyyz_1, tg_y_xxyyzz_0, tg_y_xxyyzz_1, tg_y_xxyzzz_0, tg_y_xxyzzz_1, \
                                         tg_y_xxzzzz_0, tg_y_xxzzzz_1, tg_y_xyyyyy_0, tg_y_xyyyyy_1, tg_y_xyyyyz_0, \
                                         tg_y_xyyyyz_1, tg_y_xyyyzz_0, tg_y_xyyyzz_1, tg_y_xyyzzz_0, tg_y_xyyzzz_1, \
                                         tg_y_xyzzzz_0, tg_y_xyzzzz_1, tg_y_xzzzzz_0, tg_y_xzzzzz_1, tg_y_yyyyyy_0, \
                                         tg_y_yyyyyy_1, tg_y_yyyyyz_0, tg_y_yyyyyz_1, tg_y_yyyyzz_0, tg_y_yyyyzz_1, \
                                         tg_y_yyyzzz_0, tg_y_yyyzzz_1, tg_y_yyzzzz_0, tg_y_yyzzzz_1, tg_y_yzzzzz_0, \
                                         tg_y_yzzzzz_1, tg_y_zzzzzz_0, tg_y_zzzzzz_1, tg_yy_xxxxx_1, tg_yy_xxxxxx_0, \
                                         tg_yy_xxxxxx_1, tg_yy_xxxxxy_0, tg_yy_xxxxxy_1, tg_yy_xxxxxz_0, tg_yy_xxxxxz_1, \
                                         tg_yy_xxxxy_1, tg_yy_xxxxyy_0, tg_yy_xxxxyy_1, tg_yy_xxxxyz_0, tg_yy_xxxxyz_1, \
                                         tg_yy_xxxxz_1, tg_yy_xxxxzz_0, tg_yy_xxxxzz_1, tg_yy_xxxyy_1, tg_yy_xxxyyy_0, \
                                         tg_yy_xxxyyy_1, tg_yy_xxxyyz_0, tg_yy_xxxyyz_1, tg_yy_xxxyz_1, tg_yy_xxxyzz_0, \
                                         tg_yy_xxxyzz_1, tg_yy_xxxzz_1, tg_yy_xxxzzz_0, tg_yy_xxxzzz_1, tg_yy_xxyyy_1, \
                                         tg_yy_xxyyz_1, tg_yy_xxyzz_1, tg_yy_xxzzz_1, tg_z_xxxxxx_0, tg_z_xxxxxx_1, \
                                         tg_z_xxxxxy_0, tg_z_xxxxxy_1, tg_z_xxxxxz_0, tg_z_xxxxxz_1, tg_z_xxxxyy_0, \
                                         tg_z_xxxxyy_1, tg_z_xxxxyz_0, tg_z_xxxxyz_1, tg_z_xxxxzz_0, tg_z_xxxxzz_1, \
                                         tg_z_xxxyyy_0, tg_z_xxxyyy_1, tg_z_xxxyyz_0, tg_z_xxxyyz_1, tg_z_xxxyzz_0, \
                                         tg_z_xxxyzz_1, tg_z_xxxzzz_0, tg_z_xxxzzz_1, tg_z_xxyyyy_0, tg_z_xxyyyy_1, \
                                         tg_z_xxyyyz_0, tg_z_xxyyyz_1, tg_z_xxyyzz_0, tg_z_xxyyzz_1, tg_z_xxyzzz_0, \
                                         tg_z_xxyzzz_1, tg_z_xxzzzz_0, tg_z_xxzzzz_1, tg_z_xyyyyy_0, tg_z_xyyyyy_1, \
                                         tg_z_xyyyyz_0, tg_z_xyyyyz_1, tg_z_xyyyzz_0, tg_z_xyyyzz_1, tg_z_xyyzzz_0, \
                                         tg_z_xyyzzz_1, tg_z_xyzzzz_0, tg_z_xyzzzz_1, tg_z_xzzzzz_0, tg_z_xzzzzz_1, \
                                         tg_z_yyyyyy_0, tg_z_yyyyyy_1, tg_z_yyyyyz_0, tg_z_yyyyyz_1, tg_z_yyyyzz_0, \
                                         tg_z_yyyyzz_1, tg_z_yyyzzz_0, tg_z_yyyzzz_1, tg_z_yyzzzz_0, tg_z_yyzzzz_1, \
                                         tg_z_yzzzzz_0, tg_z_yzzzzz_1, tg_z_zzzzzz_0, tg_z_zzzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxx_xxxxxx_0[j] = pb_x * tg_xx_xxxxxx_0[j] + fr * tg_xx_xxxxxx_1[j] + fl1_fx * (tg_x_xxxxxx_0[j] - tg_x_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xx_xxxxx_1[j];

                    tg_xxx_xxxxxy_0[j] = pb_x * tg_xx_xxxxxy_0[j] + fr * tg_xx_xxxxxy_1[j] + fl1_fx * (tg_x_xxxxxy_0[j] - tg_x_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xx_xxxxy_1[j];

                    tg_xxx_xxxxxz_0[j] = pb_x * tg_xx_xxxxxz_0[j] + fr * tg_xx_xxxxxz_1[j] + fl1_fx * (tg_x_xxxxxz_0[j] - tg_x_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xx_xxxxz_1[j];

                    tg_xxx_xxxxyy_0[j] = pb_x * tg_xx_xxxxyy_0[j] + fr * tg_xx_xxxxyy_1[j] + fl1_fx * (tg_x_xxxxyy_0[j] - tg_x_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xx_xxxyy_1[j];

                    tg_xxx_xxxxyz_0[j] = pb_x * tg_xx_xxxxyz_0[j] + fr * tg_xx_xxxxyz_1[j] + fl1_fx * (tg_x_xxxxyz_0[j] - tg_x_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xx_xxxyz_1[j];

                    tg_xxx_xxxxzz_0[j] = pb_x * tg_xx_xxxxzz_0[j] + fr * tg_xx_xxxxzz_1[j] + fl1_fx * (tg_x_xxxxzz_0[j] - tg_x_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xx_xxxzz_1[j];

                    tg_xxx_xxxyyy_0[j] = pb_x * tg_xx_xxxyyy_0[j] + fr * tg_xx_xxxyyy_1[j] + fl1_fx * (tg_x_xxxyyy_0[j] - tg_x_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xx_xxyyy_1[j];

                    tg_xxx_xxxyyz_0[j] = pb_x * tg_xx_xxxyyz_0[j] + fr * tg_xx_xxxyyz_1[j] + fl1_fx * (tg_x_xxxyyz_0[j] - tg_x_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xx_xxyyz_1[j];

                    tg_xxx_xxxyzz_0[j] = pb_x * tg_xx_xxxyzz_0[j] + fr * tg_xx_xxxyzz_1[j] + fl1_fx * (tg_x_xxxyzz_0[j] - tg_x_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xx_xxyzz_1[j];

                    tg_xxx_xxxzzz_0[j] = pb_x * tg_xx_xxxzzz_0[j] + fr * tg_xx_xxxzzz_1[j] + fl1_fx * (tg_x_xxxzzz_0[j] - tg_x_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xx_xxzzz_1[j];

                    tg_xxx_xxyyyy_0[j] = pb_x * tg_xx_xxyyyy_0[j] + fr * tg_xx_xxyyyy_1[j] + fl1_fx * (tg_x_xxyyyy_0[j] - tg_x_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xx_xyyyy_1[j];

                    tg_xxx_xxyyyz_0[j] = pb_x * tg_xx_xxyyyz_0[j] + fr * tg_xx_xxyyyz_1[j] + fl1_fx * (tg_x_xxyyyz_0[j] - tg_x_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xx_xyyyz_1[j];

                    tg_xxx_xxyyzz_0[j] = pb_x * tg_xx_xxyyzz_0[j] + fr * tg_xx_xxyyzz_1[j] + fl1_fx * (tg_x_xxyyzz_0[j] - tg_x_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xx_xyyzz_1[j];

                    tg_xxx_xxyzzz_0[j] = pb_x * tg_xx_xxyzzz_0[j] + fr * tg_xx_xxyzzz_1[j] + fl1_fx * (tg_x_xxyzzz_0[j] - tg_x_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xx_xyzzz_1[j];

                    tg_xxx_xxzzzz_0[j] = pb_x * tg_xx_xxzzzz_0[j] + fr * tg_xx_xxzzzz_1[j] + fl1_fx * (tg_x_xxzzzz_0[j] - tg_x_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xx_xzzzz_1[j];

                    tg_xxx_xyyyyy_0[j] = pb_x * tg_xx_xyyyyy_0[j] + fr * tg_xx_xyyyyy_1[j] + fl1_fx * (tg_x_xyyyyy_0[j] - tg_x_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xx_yyyyy_1[j];

                    tg_xxx_xyyyyz_0[j] = pb_x * tg_xx_xyyyyz_0[j] + fr * tg_xx_xyyyyz_1[j] + fl1_fx * (tg_x_xyyyyz_0[j] - tg_x_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xx_yyyyz_1[j];

                    tg_xxx_xyyyzz_0[j] = pb_x * tg_xx_xyyyzz_0[j] + fr * tg_xx_xyyyzz_1[j] + fl1_fx * (tg_x_xyyyzz_0[j] - tg_x_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xx_yyyzz_1[j];

                    tg_xxx_xyyzzz_0[j] = pb_x * tg_xx_xyyzzz_0[j] + fr * tg_xx_xyyzzz_1[j] + fl1_fx * (tg_x_xyyzzz_0[j] - tg_x_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xx_yyzzz_1[j];

                    tg_xxx_xyzzzz_0[j] = pb_x * tg_xx_xyzzzz_0[j] + fr * tg_xx_xyzzzz_1[j] + fl1_fx * (tg_x_xyzzzz_0[j] - tg_x_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xx_yzzzz_1[j];

                    tg_xxx_xzzzzz_0[j] = pb_x * tg_xx_xzzzzz_0[j] + fr * tg_xx_xzzzzz_1[j] + fl1_fx * (tg_x_xzzzzz_0[j] - tg_x_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xx_zzzzz_1[j];

                    tg_xxx_yyyyyy_0[j] = pb_x * tg_xx_yyyyyy_0[j] + fr * tg_xx_yyyyyy_1[j] + fl1_fx * (tg_x_yyyyyy_0[j] - tg_x_yyyyyy_1[j] * fl1_fza);

                    tg_xxx_yyyyyz_0[j] = pb_x * tg_xx_yyyyyz_0[j] + fr * tg_xx_yyyyyz_1[j] + fl1_fx * (tg_x_yyyyyz_0[j] - tg_x_yyyyyz_1[j] * fl1_fza);

                    tg_xxx_yyyyzz_0[j] = pb_x * tg_xx_yyyyzz_0[j] + fr * tg_xx_yyyyzz_1[j] + fl1_fx * (tg_x_yyyyzz_0[j] - tg_x_yyyyzz_1[j] * fl1_fza);

                    tg_xxx_yyyzzz_0[j] = pb_x * tg_xx_yyyzzz_0[j] + fr * tg_xx_yyyzzz_1[j] + fl1_fx * (tg_x_yyyzzz_0[j] - tg_x_yyyzzz_1[j] * fl1_fza);

                    tg_xxx_yyzzzz_0[j] = pb_x * tg_xx_yyzzzz_0[j] + fr * tg_xx_yyzzzz_1[j] + fl1_fx * (tg_x_yyzzzz_0[j] - tg_x_yyzzzz_1[j] * fl1_fza);

                    tg_xxx_yzzzzz_0[j] = pb_x * tg_xx_yzzzzz_0[j] + fr * tg_xx_yzzzzz_1[j] + fl1_fx * (tg_x_yzzzzz_0[j] - tg_x_yzzzzz_1[j] * fl1_fza);

                    tg_xxx_zzzzzz_0[j] = pb_x * tg_xx_zzzzzz_0[j] + fr * tg_xx_zzzzzz_1[j] + fl1_fx * (tg_x_zzzzzz_0[j] - tg_x_zzzzzz_1[j] * fl1_fza);

                    tg_xxy_xxxxxx_0[j] = pb_x * tg_xy_xxxxxx_0[j] + fr * tg_xy_xxxxxx_1[j] + 0.5 * fl1_fx * (tg_y_xxxxxx_0[j] - tg_y_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xy_xxxxx_1[j];

                    tg_xxy_xxxxxy_0[j] = pb_x * tg_xy_xxxxxy_0[j] + fr * tg_xy_xxxxxy_1[j] + 0.5 * fl1_fx * (tg_y_xxxxxy_0[j] - tg_y_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xy_xxxxy_1[j];

                    tg_xxy_xxxxxz_0[j] = pb_x * tg_xy_xxxxxz_0[j] + fr * tg_xy_xxxxxz_1[j] + 0.5 * fl1_fx * (tg_y_xxxxxz_0[j] - tg_y_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xy_xxxxz_1[j];

                    tg_xxy_xxxxyy_0[j] = pb_x * tg_xy_xxxxyy_0[j] + fr * tg_xy_xxxxyy_1[j] + 0.5 * fl1_fx * (tg_y_xxxxyy_0[j] - tg_y_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xy_xxxyy_1[j];

                    tg_xxy_xxxxyz_0[j] = pb_x * tg_xy_xxxxyz_0[j] + fr * tg_xy_xxxxyz_1[j] + 0.5 * fl1_fx * (tg_y_xxxxyz_0[j] - tg_y_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xy_xxxyz_1[j];

                    tg_xxy_xxxxzz_0[j] = pb_x * tg_xy_xxxxzz_0[j] + fr * tg_xy_xxxxzz_1[j] + 0.5 * fl1_fx * (tg_y_xxxxzz_0[j] - tg_y_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xy_xxxzz_1[j];

                    tg_xxy_xxxyyy_0[j] = pb_x * tg_xy_xxxyyy_0[j] + fr * tg_xy_xxxyyy_1[j] + 0.5 * fl1_fx * (tg_y_xxxyyy_0[j] - tg_y_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xy_xxyyy_1[j];

                    tg_xxy_xxxyyz_0[j] = pb_x * tg_xy_xxxyyz_0[j] + fr * tg_xy_xxxyyz_1[j] + 0.5 * fl1_fx * (tg_y_xxxyyz_0[j] - tg_y_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xy_xxyyz_1[j];

                    tg_xxy_xxxyzz_0[j] = pb_x * tg_xy_xxxyzz_0[j] + fr * tg_xy_xxxyzz_1[j] + 0.5 * fl1_fx * (tg_y_xxxyzz_0[j] - tg_y_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xy_xxyzz_1[j];

                    tg_xxy_xxxzzz_0[j] = pb_x * tg_xy_xxxzzz_0[j] + fr * tg_xy_xxxzzz_1[j] + 0.5 * fl1_fx * (tg_y_xxxzzz_0[j] - tg_y_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xy_xxzzz_1[j];

                    tg_xxy_xxyyyy_0[j] = pb_x * tg_xy_xxyyyy_0[j] + fr * tg_xy_xxyyyy_1[j] + 0.5 * fl1_fx * (tg_y_xxyyyy_0[j] - tg_y_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xy_xyyyy_1[j];

                    tg_xxy_xxyyyz_0[j] = pb_x * tg_xy_xxyyyz_0[j] + fr * tg_xy_xxyyyz_1[j] + 0.5 * fl1_fx * (tg_y_xxyyyz_0[j] - tg_y_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xy_xyyyz_1[j];

                    tg_xxy_xxyyzz_0[j] = pb_x * tg_xy_xxyyzz_0[j] + fr * tg_xy_xxyyzz_1[j] + 0.5 * fl1_fx * (tg_y_xxyyzz_0[j] - tg_y_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xy_xyyzz_1[j];

                    tg_xxy_xxyzzz_0[j] = pb_x * tg_xy_xxyzzz_0[j] + fr * tg_xy_xxyzzz_1[j] + 0.5 * fl1_fx * (tg_y_xxyzzz_0[j] - tg_y_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xy_xyzzz_1[j];

                    tg_xxy_xxzzzz_0[j] = pb_x * tg_xy_xxzzzz_0[j] + fr * tg_xy_xxzzzz_1[j] + 0.5 * fl1_fx * (tg_y_xxzzzz_0[j] - tg_y_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xy_xzzzz_1[j];

                    tg_xxy_xyyyyy_0[j] = pb_x * tg_xy_xyyyyy_0[j] + fr * tg_xy_xyyyyy_1[j] + 0.5 * fl1_fx * (tg_y_xyyyyy_0[j] - tg_y_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xy_yyyyy_1[j];

                    tg_xxy_xyyyyz_0[j] = pb_x * tg_xy_xyyyyz_0[j] + fr * tg_xy_xyyyyz_1[j] + 0.5 * fl1_fx * (tg_y_xyyyyz_0[j] - tg_y_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xy_yyyyz_1[j];

                    tg_xxy_xyyyzz_0[j] = pb_x * tg_xy_xyyyzz_0[j] + fr * tg_xy_xyyyzz_1[j] + 0.5 * fl1_fx * (tg_y_xyyyzz_0[j] - tg_y_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xy_yyyzz_1[j];

                    tg_xxy_xyyzzz_0[j] = pb_x * tg_xy_xyyzzz_0[j] + fr * tg_xy_xyyzzz_1[j] + 0.5 * fl1_fx * (tg_y_xyyzzz_0[j] - tg_y_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xy_yyzzz_1[j];

                    tg_xxy_xyzzzz_0[j] = pb_x * tg_xy_xyzzzz_0[j] + fr * tg_xy_xyzzzz_1[j] + 0.5 * fl1_fx * (tg_y_xyzzzz_0[j] - tg_y_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xy_yzzzz_1[j];

                    tg_xxy_xzzzzz_0[j] = pb_x * tg_xy_xzzzzz_0[j] + fr * tg_xy_xzzzzz_1[j] + 0.5 * fl1_fx * (tg_y_xzzzzz_0[j] - tg_y_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xy_zzzzz_1[j];

                    tg_xxy_yyyyyy_0[j] = pb_x * tg_xy_yyyyyy_0[j] + fr * tg_xy_yyyyyy_1[j] + 0.5 * fl1_fx * (tg_y_yyyyyy_0[j] - tg_y_yyyyyy_1[j] * fl1_fza);

                    tg_xxy_yyyyyz_0[j] = pb_x * tg_xy_yyyyyz_0[j] + fr * tg_xy_yyyyyz_1[j] + 0.5 * fl1_fx * (tg_y_yyyyyz_0[j] - tg_y_yyyyyz_1[j] * fl1_fza);

                    tg_xxy_yyyyzz_0[j] = pb_x * tg_xy_yyyyzz_0[j] + fr * tg_xy_yyyyzz_1[j] + 0.5 * fl1_fx * (tg_y_yyyyzz_0[j] - tg_y_yyyyzz_1[j] * fl1_fza);

                    tg_xxy_yyyzzz_0[j] = pb_x * tg_xy_yyyzzz_0[j] + fr * tg_xy_yyyzzz_1[j] + 0.5 * fl1_fx * (tg_y_yyyzzz_0[j] - tg_y_yyyzzz_1[j] * fl1_fza);

                    tg_xxy_yyzzzz_0[j] = pb_x * tg_xy_yyzzzz_0[j] + fr * tg_xy_yyzzzz_1[j] + 0.5 * fl1_fx * (tg_y_yyzzzz_0[j] - tg_y_yyzzzz_1[j] * fl1_fza);

                    tg_xxy_yzzzzz_0[j] = pb_x * tg_xy_yzzzzz_0[j] + fr * tg_xy_yzzzzz_1[j] + 0.5 * fl1_fx * (tg_y_yzzzzz_0[j] - tg_y_yzzzzz_1[j] * fl1_fza);

                    tg_xxy_zzzzzz_0[j] = pb_x * tg_xy_zzzzzz_0[j] + fr * tg_xy_zzzzzz_1[j] + 0.5 * fl1_fx * (tg_y_zzzzzz_0[j] - tg_y_zzzzzz_1[j] * fl1_fza);

                    tg_xxz_xxxxxx_0[j] = pb_x * tg_xz_xxxxxx_0[j] + fr * tg_xz_xxxxxx_1[j] + 0.5 * fl1_fx * (tg_z_xxxxxx_0[j] - tg_z_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xz_xxxxx_1[j];

                    tg_xxz_xxxxxy_0[j] = pb_x * tg_xz_xxxxxy_0[j] + fr * tg_xz_xxxxxy_1[j] + 0.5 * fl1_fx * (tg_z_xxxxxy_0[j] - tg_z_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xz_xxxxy_1[j];

                    tg_xxz_xxxxxz_0[j] = pb_x * tg_xz_xxxxxz_0[j] + fr * tg_xz_xxxxxz_1[j] + 0.5 * fl1_fx * (tg_z_xxxxxz_0[j] - tg_z_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xz_xxxxz_1[j];

                    tg_xxz_xxxxyy_0[j] = pb_x * tg_xz_xxxxyy_0[j] + fr * tg_xz_xxxxyy_1[j] + 0.5 * fl1_fx * (tg_z_xxxxyy_0[j] - tg_z_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xz_xxxyy_1[j];

                    tg_xxz_xxxxyz_0[j] = pb_x * tg_xz_xxxxyz_0[j] + fr * tg_xz_xxxxyz_1[j] + 0.5 * fl1_fx * (tg_z_xxxxyz_0[j] - tg_z_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xz_xxxyz_1[j];

                    tg_xxz_xxxxzz_0[j] = pb_x * tg_xz_xxxxzz_0[j] + fr * tg_xz_xxxxzz_1[j] + 0.5 * fl1_fx * (tg_z_xxxxzz_0[j] - tg_z_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xz_xxxzz_1[j];

                    tg_xxz_xxxyyy_0[j] = pb_x * tg_xz_xxxyyy_0[j] + fr * tg_xz_xxxyyy_1[j] + 0.5 * fl1_fx * (tg_z_xxxyyy_0[j] - tg_z_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xz_xxyyy_1[j];

                    tg_xxz_xxxyyz_0[j] = pb_x * tg_xz_xxxyyz_0[j] + fr * tg_xz_xxxyyz_1[j] + 0.5 * fl1_fx * (tg_z_xxxyyz_0[j] - tg_z_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xz_xxyyz_1[j];

                    tg_xxz_xxxyzz_0[j] = pb_x * tg_xz_xxxyzz_0[j] + fr * tg_xz_xxxyzz_1[j] + 0.5 * fl1_fx * (tg_z_xxxyzz_0[j] - tg_z_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xz_xxyzz_1[j];

                    tg_xxz_xxxzzz_0[j] = pb_x * tg_xz_xxxzzz_0[j] + fr * tg_xz_xxxzzz_1[j] + 0.5 * fl1_fx * (tg_z_xxxzzz_0[j] - tg_z_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xz_xxzzz_1[j];

                    tg_xxz_xxyyyy_0[j] = pb_x * tg_xz_xxyyyy_0[j] + fr * tg_xz_xxyyyy_1[j] + 0.5 * fl1_fx * (tg_z_xxyyyy_0[j] - tg_z_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xz_xyyyy_1[j];

                    tg_xxz_xxyyyz_0[j] = pb_x * tg_xz_xxyyyz_0[j] + fr * tg_xz_xxyyyz_1[j] + 0.5 * fl1_fx * (tg_z_xxyyyz_0[j] - tg_z_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xz_xyyyz_1[j];

                    tg_xxz_xxyyzz_0[j] = pb_x * tg_xz_xxyyzz_0[j] + fr * tg_xz_xxyyzz_1[j] + 0.5 * fl1_fx * (tg_z_xxyyzz_0[j] - tg_z_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xz_xyyzz_1[j];

                    tg_xxz_xxyzzz_0[j] = pb_x * tg_xz_xxyzzz_0[j] + fr * tg_xz_xxyzzz_1[j] + 0.5 * fl1_fx * (tg_z_xxyzzz_0[j] - tg_z_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xz_xyzzz_1[j];

                    tg_xxz_xxzzzz_0[j] = pb_x * tg_xz_xxzzzz_0[j] + fr * tg_xz_xxzzzz_1[j] + 0.5 * fl1_fx * (tg_z_xxzzzz_0[j] - tg_z_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xz_xzzzz_1[j];

                    tg_xxz_xyyyyy_0[j] = pb_x * tg_xz_xyyyyy_0[j] + fr * tg_xz_xyyyyy_1[j] + 0.5 * fl1_fx * (tg_z_xyyyyy_0[j] - tg_z_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xz_yyyyy_1[j];

                    tg_xxz_xyyyyz_0[j] = pb_x * tg_xz_xyyyyz_0[j] + fr * tg_xz_xyyyyz_1[j] + 0.5 * fl1_fx * (tg_z_xyyyyz_0[j] - tg_z_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xz_yyyyz_1[j];

                    tg_xxz_xyyyzz_0[j] = pb_x * tg_xz_xyyyzz_0[j] + fr * tg_xz_xyyyzz_1[j] + 0.5 * fl1_fx * (tg_z_xyyyzz_0[j] - tg_z_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xz_yyyzz_1[j];

                    tg_xxz_xyyzzz_0[j] = pb_x * tg_xz_xyyzzz_0[j] + fr * tg_xz_xyyzzz_1[j] + 0.5 * fl1_fx * (tg_z_xyyzzz_0[j] - tg_z_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xz_yyzzz_1[j];

                    tg_xxz_xyzzzz_0[j] = pb_x * tg_xz_xyzzzz_0[j] + fr * tg_xz_xyzzzz_1[j] + 0.5 * fl1_fx * (tg_z_xyzzzz_0[j] - tg_z_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xz_yzzzz_1[j];

                    tg_xxz_xzzzzz_0[j] = pb_x * tg_xz_xzzzzz_0[j] + fr * tg_xz_xzzzzz_1[j] + 0.5 * fl1_fx * (tg_z_xzzzzz_0[j] - tg_z_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xz_zzzzz_1[j];

                    tg_xxz_yyyyyy_0[j] = pb_x * tg_xz_yyyyyy_0[j] + fr * tg_xz_yyyyyy_1[j] + 0.5 * fl1_fx * (tg_z_yyyyyy_0[j] - tg_z_yyyyyy_1[j] * fl1_fza);

                    tg_xxz_yyyyyz_0[j] = pb_x * tg_xz_yyyyyz_0[j] + fr * tg_xz_yyyyyz_1[j] + 0.5 * fl1_fx * (tg_z_yyyyyz_0[j] - tg_z_yyyyyz_1[j] * fl1_fza);

                    tg_xxz_yyyyzz_0[j] = pb_x * tg_xz_yyyyzz_0[j] + fr * tg_xz_yyyyzz_1[j] + 0.5 * fl1_fx * (tg_z_yyyyzz_0[j] - tg_z_yyyyzz_1[j] * fl1_fza);

                    tg_xxz_yyyzzz_0[j] = pb_x * tg_xz_yyyzzz_0[j] + fr * tg_xz_yyyzzz_1[j] + 0.5 * fl1_fx * (tg_z_yyyzzz_0[j] - tg_z_yyyzzz_1[j] * fl1_fza);

                    tg_xxz_yyzzzz_0[j] = pb_x * tg_xz_yyzzzz_0[j] + fr * tg_xz_yyzzzz_1[j] + 0.5 * fl1_fx * (tg_z_yyzzzz_0[j] - tg_z_yyzzzz_1[j] * fl1_fza);

                    tg_xxz_yzzzzz_0[j] = pb_x * tg_xz_yzzzzz_0[j] + fr * tg_xz_yzzzzz_1[j] + 0.5 * fl1_fx * (tg_z_yzzzzz_0[j] - tg_z_yzzzzz_1[j] * fl1_fza);

                    tg_xxz_zzzzzz_0[j] = pb_x * tg_xz_zzzzzz_0[j] + fr * tg_xz_zzzzzz_1[j] + 0.5 * fl1_fx * (tg_z_zzzzzz_0[j] - tg_z_zzzzzz_1[j] * fl1_fza);

                    tg_xyy_xxxxxx_0[j] = pb_x * tg_yy_xxxxxx_0[j] + fr * tg_yy_xxxxxx_1[j] + 3.0 * fl1_fxn * tg_yy_xxxxx_1[j];

                    tg_xyy_xxxxxy_0[j] = pb_x * tg_yy_xxxxxy_0[j] + fr * tg_yy_xxxxxy_1[j] + 2.5 * fl1_fxn * tg_yy_xxxxy_1[j];

                    tg_xyy_xxxxxz_0[j] = pb_x * tg_yy_xxxxxz_0[j] + fr * tg_yy_xxxxxz_1[j] + 2.5 * fl1_fxn * tg_yy_xxxxz_1[j];

                    tg_xyy_xxxxyy_0[j] = pb_x * tg_yy_xxxxyy_0[j] + fr * tg_yy_xxxxyy_1[j] + 2.0 * fl1_fxn * tg_yy_xxxyy_1[j];

                    tg_xyy_xxxxyz_0[j] = pb_x * tg_yy_xxxxyz_0[j] + fr * tg_yy_xxxxyz_1[j] + 2.0 * fl1_fxn * tg_yy_xxxyz_1[j];

                    tg_xyy_xxxxzz_0[j] = pb_x * tg_yy_xxxxzz_0[j] + fr * tg_yy_xxxxzz_1[j] + 2.0 * fl1_fxn * tg_yy_xxxzz_1[j];

                    tg_xyy_xxxyyy_0[j] = pb_x * tg_yy_xxxyyy_0[j] + fr * tg_yy_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_yy_xxyyy_1[j];

                    tg_xyy_xxxyyz_0[j] = pb_x * tg_yy_xxxyyz_0[j] + fr * tg_yy_xxxyyz_1[j] + 1.5 * fl1_fxn * tg_yy_xxyyz_1[j];

                    tg_xyy_xxxyzz_0[j] = pb_x * tg_yy_xxxyzz_0[j] + fr * tg_yy_xxxyzz_1[j] + 1.5 * fl1_fxn * tg_yy_xxyzz_1[j];

                    tg_xyy_xxxzzz_0[j] = pb_x * tg_yy_xxxzzz_0[j] + fr * tg_yy_xxxzzz_1[j] + 1.5 * fl1_fxn * tg_yy_xxzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSFSI_94_187(      CMemBlock2D<double>* primBuffer,
                                        const CRecursionMap&       recursionMap,
                                        const CMemBlock2D<double>& osFactors,
                                        const CMemBlock2D<double>& wpDistances,
                                        const CGtoPairsBlock&      braGtoPairsBlock,
                                        const CGtoPairsBlock&      ketGtoPairsBlock,
                                        const int32_t              nKetPrimPairs,
                                        const int32_t              iContrPair)
    {
        // Batch of Integrals (94,187)

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
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_3_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_3_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_2_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_2_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_2_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_1_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {1, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_1_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {1, -1, -1, -1}, {6, -1, -1, -1}, 
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

                auto tg_yy_xxxxxx_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 84); 

                auto tg_yy_xxxxxy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 85); 

                auto tg_yy_xxxxxz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 86); 

                auto tg_yy_xxxxyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 87); 

                auto tg_yy_xxxxyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 88); 

                auto tg_yy_xxxxzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 89); 

                auto tg_yy_xxxyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 90); 

                auto tg_yy_xxxyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 91); 

                auto tg_yy_xxxyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 92); 

                auto tg_yy_xxxzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 93); 

                auto tg_yy_xxyyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 94); 

                auto tg_yy_xxyyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 95); 

                auto tg_yy_xxyyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 96); 

                auto tg_yy_xxyzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 97); 

                auto tg_yy_xxzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 98); 

                auto tg_yy_xyyyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 99); 

                auto tg_yy_xyyyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 100); 

                auto tg_yy_xyyyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 101); 

                auto tg_yy_xyyzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 102); 

                auto tg_yy_xyzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 103); 

                auto tg_yy_xzzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 104); 

                auto tg_yy_yyyyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 105); 

                auto tg_yy_yyyyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 106); 

                auto tg_yy_yyyyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 107); 

                auto tg_yy_yyyzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 108); 

                auto tg_yy_yyzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 109); 

                auto tg_yy_yzzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 110); 

                auto tg_yy_zzzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 111); 

                auto tg_yz_xxxxxx_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 112); 

                auto tg_yz_xxxxxy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 113); 

                auto tg_yz_xxxxxz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 114); 

                auto tg_yz_xxxxyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 115); 

                auto tg_yz_xxxxyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 116); 

                auto tg_yz_xxxxzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 117); 

                auto tg_yz_xxxyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 118); 

                auto tg_yz_xxxyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 119); 

                auto tg_yz_xxxyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 120); 

                auto tg_yz_xxxzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 121); 

                auto tg_yz_xxyyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 122); 

                auto tg_yz_xxyyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 123); 

                auto tg_yz_xxyyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 124); 

                auto tg_yz_xxyzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 125); 

                auto tg_yz_xxzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 126); 

                auto tg_yz_xyyyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 127); 

                auto tg_yz_xyyyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 128); 

                auto tg_yz_xyyyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 129); 

                auto tg_yz_xyyzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 130); 

                auto tg_yz_xyzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 131); 

                auto tg_yz_xzzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 132); 

                auto tg_yz_yyyyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 133); 

                auto tg_yz_yyyyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 134); 

                auto tg_yz_yyyyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 135); 

                auto tg_yz_yyyzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 136); 

                auto tg_yz_yyzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 137); 

                auto tg_yz_yzzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 138); 

                auto tg_yz_zzzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 139); 

                auto tg_zz_xxxxxx_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 140); 

                auto tg_zz_xxxxxy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 141); 

                auto tg_zz_xxxxxz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 142); 

                auto tg_zz_xxxxyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 143); 

                auto tg_zz_xxxxyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 144); 

                auto tg_zz_xxxxzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 145); 

                auto tg_zz_xxxyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 146); 

                auto tg_zz_xxxyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 147); 

                auto tg_zz_xxxyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 148); 

                auto tg_zz_xxxzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 149); 

                auto tg_zz_xxyyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 150); 

                auto tg_zz_xxyyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 151); 

                auto tg_zz_xxyyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 152); 

                auto tg_zz_xxyzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 153); 

                auto tg_zz_xxzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 154); 

                auto tg_zz_xyyyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 155); 

                auto tg_zz_xyyyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 156); 

                auto tg_zz_xyyyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 157); 

                auto tg_zz_xyyzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 158); 

                auto tg_zz_xyzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 159); 

                auto tg_zz_xzzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 160); 

                auto tg_zz_yyyyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 161); 

                auto tg_zz_yyyyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 162); 

                auto tg_zz_yyyyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 163); 

                auto tg_zz_yyyzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 164); 

                auto tg_zz_yyzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 165); 

                auto tg_zz_yzzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 166); 

                auto tg_zz_zzzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 167); 

                auto tg_yy_xxxxxx_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 84); 

                auto tg_yy_xxxxxy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 85); 

                auto tg_yy_xxxxxz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 86); 

                auto tg_yy_xxxxyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 87); 

                auto tg_yy_xxxxyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 88); 

                auto tg_yy_xxxxzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 89); 

                auto tg_yy_xxxyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 90); 

                auto tg_yy_xxxyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 91); 

                auto tg_yy_xxxyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 92); 

                auto tg_yy_xxxzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 93); 

                auto tg_yy_xxyyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 94); 

                auto tg_yy_xxyyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 95); 

                auto tg_yy_xxyyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 96); 

                auto tg_yy_xxyzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 97); 

                auto tg_yy_xxzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 98); 

                auto tg_yy_xyyyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 99); 

                auto tg_yy_xyyyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 100); 

                auto tg_yy_xyyyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 101); 

                auto tg_yy_xyyzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 102); 

                auto tg_yy_xyzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 103); 

                auto tg_yy_xzzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 104); 

                auto tg_yy_yyyyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 105); 

                auto tg_yy_yyyyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 106); 

                auto tg_yy_yyyyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 107); 

                auto tg_yy_yyyzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 108); 

                auto tg_yy_yyzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 109); 

                auto tg_yy_yzzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 110); 

                auto tg_yy_zzzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 111); 

                auto tg_yz_xxxxxx_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 112); 

                auto tg_yz_xxxxxy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 113); 

                auto tg_yz_xxxxxz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 114); 

                auto tg_yz_xxxxyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 115); 

                auto tg_yz_xxxxyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 116); 

                auto tg_yz_xxxxzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 117); 

                auto tg_yz_xxxyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 118); 

                auto tg_yz_xxxyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 119); 

                auto tg_yz_xxxyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 120); 

                auto tg_yz_xxxzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 121); 

                auto tg_yz_xxyyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 122); 

                auto tg_yz_xxyyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 123); 

                auto tg_yz_xxyyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 124); 

                auto tg_yz_xxyzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 125); 

                auto tg_yz_xxzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 126); 

                auto tg_yz_xyyyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 127); 

                auto tg_yz_xyyyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 128); 

                auto tg_yz_xyyyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 129); 

                auto tg_yz_xyyzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 130); 

                auto tg_yz_xyzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 131); 

                auto tg_yz_xzzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 132); 

                auto tg_yz_yyyyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 133); 

                auto tg_yz_yyyyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 134); 

                auto tg_yz_yyyyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 135); 

                auto tg_yz_yyyzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 136); 

                auto tg_yz_yyzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 137); 

                auto tg_yz_yzzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 138); 

                auto tg_yz_zzzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 139); 

                auto tg_zz_xxxxxx_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 140); 

                auto tg_zz_xxxxxy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 141); 

                auto tg_zz_xxxxxz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 142); 

                auto tg_zz_xxxxyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 143); 

                auto tg_zz_xxxxyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 144); 

                auto tg_zz_xxxxzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 145); 

                auto tg_zz_xxxyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 146); 

                auto tg_zz_xxxyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 147); 

                auto tg_zz_xxxyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 148); 

                auto tg_zz_xxxzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 149); 

                auto tg_zz_xxyyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 150); 

                auto tg_zz_xxyyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 151); 

                auto tg_zz_xxyyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 152); 

                auto tg_zz_xxyzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 153); 

                auto tg_zz_xxzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 154); 

                auto tg_zz_xyyyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 155); 

                auto tg_zz_xyyyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 156); 

                auto tg_zz_xyyyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 157); 

                auto tg_zz_xyyzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 158); 

                auto tg_zz_xyzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 159); 

                auto tg_zz_xzzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 160); 

                auto tg_zz_yyyyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 161); 

                auto tg_zz_yyyyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 162); 

                auto tg_zz_yyyyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 163); 

                auto tg_zz_yyyzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 164); 

                auto tg_zz_yyzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 165); 

                auto tg_zz_yzzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 166); 

                auto tg_zz_zzzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 167); 

                auto tg_y_xxxxxx_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 28); 

                auto tg_y_xxxxxy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 29); 

                auto tg_y_xxxxxz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 30); 

                auto tg_y_xxxxyy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 31); 

                auto tg_y_xxxxyz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 32); 

                auto tg_y_xxxxzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 33); 

                auto tg_y_xxxyyy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 34); 

                auto tg_y_xxxyyz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 35); 

                auto tg_y_xxxyzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 36); 

                auto tg_y_xxxzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 37); 

                auto tg_y_xxyyyy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 38); 

                auto tg_y_xxyyyz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 39); 

                auto tg_y_xxyyzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 40); 

                auto tg_y_xxyzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 41); 

                auto tg_y_xxzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 42); 

                auto tg_y_xyyyyy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 43); 

                auto tg_y_xyyyyz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 44); 

                auto tg_y_xyyyzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 45); 

                auto tg_y_xyyzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 46); 

                auto tg_y_xxxxxx_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 28); 

                auto tg_y_xxxxxy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 29); 

                auto tg_y_xxxxxz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 30); 

                auto tg_y_xxxxyy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 31); 

                auto tg_y_xxxxyz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 32); 

                auto tg_y_xxxxzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 33); 

                auto tg_y_xxxyyy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 34); 

                auto tg_y_xxxyyz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 35); 

                auto tg_y_xxxyzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 36); 

                auto tg_y_xxxzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 37); 

                auto tg_y_xxyyyy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 38); 

                auto tg_y_xxyyyz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 39); 

                auto tg_y_xxyyzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 40); 

                auto tg_y_xxyzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 41); 

                auto tg_y_xxzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 42); 

                auto tg_y_xyyyyy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 43); 

                auto tg_y_xyyyyz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 44); 

                auto tg_y_xyyyzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 45); 

                auto tg_y_xyyzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 46); 

                auto tg_yy_xxxxx_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 63); 

                auto tg_yy_xxxxy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 64); 

                auto tg_yy_xxxxz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 65); 

                auto tg_yy_xxxyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 66); 

                auto tg_yy_xxxyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 67); 

                auto tg_yy_xxxzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 68); 

                auto tg_yy_xxyyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 69); 

                auto tg_yy_xxyyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 70); 

                auto tg_yy_xxyzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 71); 

                auto tg_yy_xxzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 72); 

                auto tg_yy_xyyyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 73); 

                auto tg_yy_xyyyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 74); 

                auto tg_yy_xyyzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 75); 

                auto tg_yy_xyzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 76); 

                auto tg_yy_xzzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 77); 

                auto tg_yy_yyyyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 78); 

                auto tg_yy_yyyyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 79); 

                auto tg_yy_yyyzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 80); 

                auto tg_yy_yyzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 81); 

                auto tg_yy_yzzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 82); 

                auto tg_yy_zzzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 83); 

                auto tg_yz_xxxxx_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 84); 

                auto tg_yz_xxxxy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 85); 

                auto tg_yz_xxxxz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 86); 

                auto tg_yz_xxxyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 87); 

                auto tg_yz_xxxyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 88); 

                auto tg_yz_xxxzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 89); 

                auto tg_yz_xxyyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 90); 

                auto tg_yz_xxyyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 91); 

                auto tg_yz_xxyzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 92); 

                auto tg_yz_xxzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 93); 

                auto tg_yz_xyyyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 94); 

                auto tg_yz_xyyyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 95); 

                auto tg_yz_xyyzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 96); 

                auto tg_yz_xyzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 97); 

                auto tg_yz_xzzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 98); 

                auto tg_yz_yyyyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 99); 

                auto tg_yz_yyyyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 100); 

                auto tg_yz_yyyzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 101); 

                auto tg_yz_yyzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 102); 

                auto tg_yz_yzzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 103); 

                auto tg_yz_zzzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 104); 

                auto tg_zz_xxxxx_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 105); 

                auto tg_zz_xxxxy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 106); 

                auto tg_zz_xxxxz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 107); 

                auto tg_zz_xxxyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 108); 

                auto tg_zz_xxxyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 109); 

                auto tg_zz_xxxzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 110); 

                auto tg_zz_xxyyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 111); 

                auto tg_zz_xxyyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 112); 

                auto tg_zz_xxyzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 113); 

                auto tg_zz_xxzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 114); 

                auto tg_zz_xyyyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 115); 

                auto tg_zz_xyyyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 116); 

                auto tg_zz_xyyzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 117); 

                auto tg_zz_xyzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 118); 

                auto tg_zz_xzzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 119); 

                auto tg_zz_yyyyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 120); 

                auto tg_zz_yyyyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 121); 

                auto tg_zz_yyyzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 122); 

                auto tg_zz_yyzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 123); 

                auto tg_zz_yzzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 124); 

                auto tg_zz_zzzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 125); 

                // set up pointers to integrals

                auto tg_xyy_xxyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 94); 

                auto tg_xyy_xxyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 95); 

                auto tg_xyy_xxyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 96); 

                auto tg_xyy_xxyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 97); 

                auto tg_xyy_xxzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 98); 

                auto tg_xyy_xyyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 99); 

                auto tg_xyy_xyyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 100); 

                auto tg_xyy_xyyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 101); 

                auto tg_xyy_xyyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 102); 

                auto tg_xyy_xyzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 103); 

                auto tg_xyy_xzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 104); 

                auto tg_xyy_yyyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 105); 

                auto tg_xyy_yyyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 106); 

                auto tg_xyy_yyyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 107); 

                auto tg_xyy_yyyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 108); 

                auto tg_xyy_yyzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 109); 

                auto tg_xyy_yzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 110); 

                auto tg_xyy_zzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 111); 

                auto tg_xyz_xxxxxx_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 112); 

                auto tg_xyz_xxxxxy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 113); 

                auto tg_xyz_xxxxxz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 114); 

                auto tg_xyz_xxxxyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 115); 

                auto tg_xyz_xxxxyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 116); 

                auto tg_xyz_xxxxzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 117); 

                auto tg_xyz_xxxyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 118); 

                auto tg_xyz_xxxyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 119); 

                auto tg_xyz_xxxyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 120); 

                auto tg_xyz_xxxzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 121); 

                auto tg_xyz_xxyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 122); 

                auto tg_xyz_xxyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 123); 

                auto tg_xyz_xxyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 124); 

                auto tg_xyz_xxyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 125); 

                auto tg_xyz_xxzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 126); 

                auto tg_xyz_xyyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 127); 

                auto tg_xyz_xyyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 128); 

                auto tg_xyz_xyyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 129); 

                auto tg_xyz_xyyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 130); 

                auto tg_xyz_xyzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 131); 

                auto tg_xyz_xzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 132); 

                auto tg_xyz_yyyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 133); 

                auto tg_xyz_yyyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 134); 

                auto tg_xyz_yyyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 135); 

                auto tg_xyz_yyyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 136); 

                auto tg_xyz_yyzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 137); 

                auto tg_xyz_yzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 138); 

                auto tg_xyz_zzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 139); 

                auto tg_xzz_xxxxxx_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 140); 

                auto tg_xzz_xxxxxy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 141); 

                auto tg_xzz_xxxxxz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 142); 

                auto tg_xzz_xxxxyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 143); 

                auto tg_xzz_xxxxyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 144); 

                auto tg_xzz_xxxxzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 145); 

                auto tg_xzz_xxxyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 146); 

                auto tg_xzz_xxxyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 147); 

                auto tg_xzz_xxxyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 148); 

                auto tg_xzz_xxxzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 149); 

                auto tg_xzz_xxyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 150); 

                auto tg_xzz_xxyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 151); 

                auto tg_xzz_xxyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 152); 

                auto tg_xzz_xxyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 153); 

                auto tg_xzz_xxzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 154); 

                auto tg_xzz_xyyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 155); 

                auto tg_xzz_xyyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 156); 

                auto tg_xzz_xyyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 157); 

                auto tg_xzz_xyyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 158); 

                auto tg_xzz_xyzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 159); 

                auto tg_xzz_xzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 160); 

                auto tg_xzz_yyyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 161); 

                auto tg_xzz_yyyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 162); 

                auto tg_xzz_yyyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 163); 

                auto tg_xzz_yyyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 164); 

                auto tg_xzz_yyzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 165); 

                auto tg_xzz_yzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 166); 

                auto tg_xzz_zzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 167); 

                auto tg_yyy_xxxxxx_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 168); 

                auto tg_yyy_xxxxxy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 169); 

                auto tg_yyy_xxxxxz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 170); 

                auto tg_yyy_xxxxyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 171); 

                auto tg_yyy_xxxxyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 172); 

                auto tg_yyy_xxxxzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 173); 

                auto tg_yyy_xxxyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 174); 

                auto tg_yyy_xxxyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 175); 

                auto tg_yyy_xxxyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 176); 

                auto tg_yyy_xxxzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 177); 

                auto tg_yyy_xxyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 178); 

                auto tg_yyy_xxyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 179); 

                auto tg_yyy_xxyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 180); 

                auto tg_yyy_xxyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 181); 

                auto tg_yyy_xxzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 182); 

                auto tg_yyy_xyyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 183); 

                auto tg_yyy_xyyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 184); 

                auto tg_yyy_xyyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 185); 

                auto tg_yyy_xyyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 186); 

                // Batch of Integrals (94,187)

                #pragma omp simd aligned(fxn, fza, tg_xyy_xxyyyy_0, tg_xyy_xxyyyz_0, tg_xyy_xxyyzz_0, \
                                         tg_xyy_xxyzzz_0, tg_xyy_xxzzzz_0, tg_xyy_xyyyyy_0, tg_xyy_xyyyyz_0, tg_xyy_xyyyzz_0, \
                                         tg_xyy_xyyzzz_0, tg_xyy_xyzzzz_0, tg_xyy_xzzzzz_0, tg_xyy_yyyyyy_0, tg_xyy_yyyyyz_0, \
                                         tg_xyy_yyyyzz_0, tg_xyy_yyyzzz_0, tg_xyy_yyzzzz_0, tg_xyy_yzzzzz_0, tg_xyy_zzzzzz_0, \
                                         tg_xyz_xxxxxx_0, tg_xyz_xxxxxy_0, tg_xyz_xxxxxz_0, tg_xyz_xxxxyy_0, tg_xyz_xxxxyz_0, \
                                         tg_xyz_xxxxzz_0, tg_xyz_xxxyyy_0, tg_xyz_xxxyyz_0, tg_xyz_xxxyzz_0, tg_xyz_xxxzzz_0, \
                                         tg_xyz_xxyyyy_0, tg_xyz_xxyyyz_0, tg_xyz_xxyyzz_0, tg_xyz_xxyzzz_0, tg_xyz_xxzzzz_0, \
                                         tg_xyz_xyyyyy_0, tg_xyz_xyyyyz_0, tg_xyz_xyyyzz_0, tg_xyz_xyyzzz_0, tg_xyz_xyzzzz_0, \
                                         tg_xyz_xzzzzz_0, tg_xyz_yyyyyy_0, tg_xyz_yyyyyz_0, tg_xyz_yyyyzz_0, tg_xyz_yyyzzz_0, \
                                         tg_xyz_yyzzzz_0, tg_xyz_yzzzzz_0, tg_xyz_zzzzzz_0, tg_xzz_xxxxxx_0, tg_xzz_xxxxxy_0, \
                                         tg_xzz_xxxxxz_0, tg_xzz_xxxxyy_0, tg_xzz_xxxxyz_0, tg_xzz_xxxxzz_0, tg_xzz_xxxyyy_0, \
                                         tg_xzz_xxxyyz_0, tg_xzz_xxxyzz_0, tg_xzz_xxxzzz_0, tg_xzz_xxyyyy_0, tg_xzz_xxyyyz_0, \
                                         tg_xzz_xxyyzz_0, tg_xzz_xxyzzz_0, tg_xzz_xxzzzz_0, tg_xzz_xyyyyy_0, tg_xzz_xyyyyz_0, \
                                         tg_xzz_xyyyzz_0, tg_xzz_xyyzzz_0, tg_xzz_xyzzzz_0, tg_xzz_xzzzzz_0, tg_xzz_yyyyyy_0, \
                                         tg_xzz_yyyyyz_0, tg_xzz_yyyyzz_0, tg_xzz_yyyzzz_0, tg_xzz_yyzzzz_0, tg_xzz_yzzzzz_0, \
                                         tg_xzz_zzzzzz_0, tg_y_xxxxxx_0, tg_y_xxxxxx_1, tg_y_xxxxxy_0, tg_y_xxxxxy_1, \
                                         tg_y_xxxxxz_0, tg_y_xxxxxz_1, tg_y_xxxxyy_0, tg_y_xxxxyy_1, tg_y_xxxxyz_0, \
                                         tg_y_xxxxyz_1, tg_y_xxxxzz_0, tg_y_xxxxzz_1, tg_y_xxxyyy_0, tg_y_xxxyyy_1, \
                                         tg_y_xxxyyz_0, tg_y_xxxyyz_1, tg_y_xxxyzz_0, tg_y_xxxyzz_1, tg_y_xxxzzz_0, \
                                         tg_y_xxxzzz_1, tg_y_xxyyyy_0, tg_y_xxyyyy_1, tg_y_xxyyyz_0, tg_y_xxyyyz_1, \
                                         tg_y_xxyyzz_0, tg_y_xxyyzz_1, tg_y_xxyzzz_0, tg_y_xxyzzz_1, tg_y_xxzzzz_0, \
                                         tg_y_xxzzzz_1, tg_y_xyyyyy_0, tg_y_xyyyyy_1, tg_y_xyyyyz_0, tg_y_xyyyyz_1, \
                                         tg_y_xyyyzz_0, tg_y_xyyyzz_1, tg_y_xyyzzz_0, tg_y_xyyzzz_1, tg_yy_xxxxx_1, \
                                         tg_yy_xxxxxx_0, tg_yy_xxxxxx_1, tg_yy_xxxxxy_0, tg_yy_xxxxxy_1, tg_yy_xxxxxz_0, \
                                         tg_yy_xxxxxz_1, tg_yy_xxxxy_1, tg_yy_xxxxyy_0, tg_yy_xxxxyy_1, tg_yy_xxxxyz_0, \
                                         tg_yy_xxxxyz_1, tg_yy_xxxxz_1, tg_yy_xxxxzz_0, tg_yy_xxxxzz_1, tg_yy_xxxyy_1, \
                                         tg_yy_xxxyyy_0, tg_yy_xxxyyy_1, tg_yy_xxxyyz_0, tg_yy_xxxyyz_1, tg_yy_xxxyz_1, \
                                         tg_yy_xxxyzz_0, tg_yy_xxxyzz_1, tg_yy_xxxzz_1, tg_yy_xxxzzz_0, tg_yy_xxxzzz_1, \
                                         tg_yy_xxyyy_1, tg_yy_xxyyyy_0, tg_yy_xxyyyy_1, tg_yy_xxyyyz_0, tg_yy_xxyyyz_1, \
                                         tg_yy_xxyyz_1, tg_yy_xxyyzz_0, tg_yy_xxyyzz_1, tg_yy_xxyzz_1, tg_yy_xxyzzz_0, \
                                         tg_yy_xxyzzz_1, tg_yy_xxzzz_1, tg_yy_xxzzzz_0, tg_yy_xxzzzz_1, tg_yy_xyyyy_1, \
                                         tg_yy_xyyyyy_0, tg_yy_xyyyyy_1, tg_yy_xyyyyz_0, tg_yy_xyyyyz_1, tg_yy_xyyyz_1, \
                                         tg_yy_xyyyzz_0, tg_yy_xyyyzz_1, tg_yy_xyyzz_1, tg_yy_xyyzzz_0, tg_yy_xyyzzz_1, \
                                         tg_yy_xyzzz_1, tg_yy_xyzzzz_0, tg_yy_xyzzzz_1, tg_yy_xzzzz_1, tg_yy_xzzzzz_0, \
                                         tg_yy_xzzzzz_1, tg_yy_yyyyy_1, tg_yy_yyyyyy_0, tg_yy_yyyyyy_1, tg_yy_yyyyyz_0, \
                                         tg_yy_yyyyyz_1, tg_yy_yyyyz_1, tg_yy_yyyyzz_0, tg_yy_yyyyzz_1, tg_yy_yyyzz_1, \
                                         tg_yy_yyyzzz_0, tg_yy_yyyzzz_1, tg_yy_yyzzz_1, tg_yy_yyzzzz_0, tg_yy_yyzzzz_1, \
                                         tg_yy_yzzzz_1, tg_yy_yzzzzz_0, tg_yy_yzzzzz_1, tg_yy_zzzzz_1, tg_yy_zzzzzz_0, \
                                         tg_yy_zzzzzz_1, tg_yyy_xxxxxx_0, tg_yyy_xxxxxy_0, tg_yyy_xxxxxz_0, tg_yyy_xxxxyy_0, \
                                         tg_yyy_xxxxyz_0, tg_yyy_xxxxzz_0, tg_yyy_xxxyyy_0, tg_yyy_xxxyyz_0, tg_yyy_xxxyzz_0, \
                                         tg_yyy_xxxzzz_0, tg_yyy_xxyyyy_0, tg_yyy_xxyyyz_0, tg_yyy_xxyyzz_0, tg_yyy_xxyzzz_0, \
                                         tg_yyy_xxzzzz_0, tg_yyy_xyyyyy_0, tg_yyy_xyyyyz_0, tg_yyy_xyyyzz_0, tg_yyy_xyyzzz_0, \
                                         tg_yz_xxxxx_1, tg_yz_xxxxxx_0, tg_yz_xxxxxx_1, tg_yz_xxxxxy_0, tg_yz_xxxxxy_1, \
                                         tg_yz_xxxxxz_0, tg_yz_xxxxxz_1, tg_yz_xxxxy_1, tg_yz_xxxxyy_0, tg_yz_xxxxyy_1, \
                                         tg_yz_xxxxyz_0, tg_yz_xxxxyz_1, tg_yz_xxxxz_1, tg_yz_xxxxzz_0, tg_yz_xxxxzz_1, \
                                         tg_yz_xxxyy_1, tg_yz_xxxyyy_0, tg_yz_xxxyyy_1, tg_yz_xxxyyz_0, tg_yz_xxxyyz_1, \
                                         tg_yz_xxxyz_1, tg_yz_xxxyzz_0, tg_yz_xxxyzz_1, tg_yz_xxxzz_1, tg_yz_xxxzzz_0, \
                                         tg_yz_xxxzzz_1, tg_yz_xxyyy_1, tg_yz_xxyyyy_0, tg_yz_xxyyyy_1, tg_yz_xxyyyz_0, \
                                         tg_yz_xxyyyz_1, tg_yz_xxyyz_1, tg_yz_xxyyzz_0, tg_yz_xxyyzz_1, tg_yz_xxyzz_1, \
                                         tg_yz_xxyzzz_0, tg_yz_xxyzzz_1, tg_yz_xxzzz_1, tg_yz_xxzzzz_0, tg_yz_xxzzzz_1, \
                                         tg_yz_xyyyy_1, tg_yz_xyyyyy_0, tg_yz_xyyyyy_1, tg_yz_xyyyyz_0, tg_yz_xyyyyz_1, \
                                         tg_yz_xyyyz_1, tg_yz_xyyyzz_0, tg_yz_xyyyzz_1, tg_yz_xyyzz_1, tg_yz_xyyzzz_0, \
                                         tg_yz_xyyzzz_1, tg_yz_xyzzz_1, tg_yz_xyzzzz_0, tg_yz_xyzzzz_1, tg_yz_xzzzz_1, \
                                         tg_yz_xzzzzz_0, tg_yz_xzzzzz_1, tg_yz_yyyyy_1, tg_yz_yyyyyy_0, tg_yz_yyyyyy_1, \
                                         tg_yz_yyyyyz_0, tg_yz_yyyyyz_1, tg_yz_yyyyz_1, tg_yz_yyyyzz_0, tg_yz_yyyyzz_1, \
                                         tg_yz_yyyzz_1, tg_yz_yyyzzz_0, tg_yz_yyyzzz_1, tg_yz_yyzzz_1, tg_yz_yyzzzz_0, \
                                         tg_yz_yyzzzz_1, tg_yz_yzzzz_1, tg_yz_yzzzzz_0, tg_yz_yzzzzz_1, tg_yz_zzzzz_1, \
                                         tg_yz_zzzzzz_0, tg_yz_zzzzzz_1, tg_zz_xxxxx_1, tg_zz_xxxxxx_0, tg_zz_xxxxxx_1, \
                                         tg_zz_xxxxxy_0, tg_zz_xxxxxy_1, tg_zz_xxxxxz_0, tg_zz_xxxxxz_1, tg_zz_xxxxy_1, \
                                         tg_zz_xxxxyy_0, tg_zz_xxxxyy_1, tg_zz_xxxxyz_0, tg_zz_xxxxyz_1, tg_zz_xxxxz_1, \
                                         tg_zz_xxxxzz_0, tg_zz_xxxxzz_1, tg_zz_xxxyy_1, tg_zz_xxxyyy_0, tg_zz_xxxyyy_1, \
                                         tg_zz_xxxyyz_0, tg_zz_xxxyyz_1, tg_zz_xxxyz_1, tg_zz_xxxyzz_0, tg_zz_xxxyzz_1, \
                                         tg_zz_xxxzz_1, tg_zz_xxxzzz_0, tg_zz_xxxzzz_1, tg_zz_xxyyy_1, tg_zz_xxyyyy_0, \
                                         tg_zz_xxyyyy_1, tg_zz_xxyyyz_0, tg_zz_xxyyyz_1, tg_zz_xxyyz_1, tg_zz_xxyyzz_0, \
                                         tg_zz_xxyyzz_1, tg_zz_xxyzz_1, tg_zz_xxyzzz_0, tg_zz_xxyzzz_1, tg_zz_xxzzz_1, \
                                         tg_zz_xxzzzz_0, tg_zz_xxzzzz_1, tg_zz_xyyyy_1, tg_zz_xyyyyy_0, tg_zz_xyyyyy_1, \
                                         tg_zz_xyyyyz_0, tg_zz_xyyyyz_1, tg_zz_xyyyz_1, tg_zz_xyyyzz_0, tg_zz_xyyyzz_1, \
                                         tg_zz_xyyzz_1, tg_zz_xyyzzz_0, tg_zz_xyyzzz_1, tg_zz_xyzzz_1, tg_zz_xyzzzz_0, \
                                         tg_zz_xyzzzz_1, tg_zz_xzzzz_1, tg_zz_xzzzzz_0, tg_zz_xzzzzz_1, tg_zz_yyyyy_1, \
                                         tg_zz_yyyyyy_0, tg_zz_yyyyyy_1, tg_zz_yyyyyz_0, tg_zz_yyyyyz_1, tg_zz_yyyyz_1, \
                                         tg_zz_yyyyzz_0, tg_zz_yyyyzz_1, tg_zz_yyyzz_1, tg_zz_yyyzzz_0, tg_zz_yyyzzz_1, \
                                         tg_zz_yyzzz_1, tg_zz_yyzzzz_0, tg_zz_yyzzzz_1, tg_zz_yzzzz_1, tg_zz_yzzzzz_0, \
                                         tg_zz_yzzzzz_1, tg_zz_zzzzz_1, tg_zz_zzzzzz_0, tg_zz_zzzzzz_1, wp_x, wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xyy_xxyyyy_0[j] = pb_x * tg_yy_xxyyyy_0[j] + fr * tg_yy_xxyyyy_1[j] + fl1_fxn * tg_yy_xyyyy_1[j];

                    tg_xyy_xxyyyz_0[j] = pb_x * tg_yy_xxyyyz_0[j] + fr * tg_yy_xxyyyz_1[j] + fl1_fxn * tg_yy_xyyyz_1[j];

                    tg_xyy_xxyyzz_0[j] = pb_x * tg_yy_xxyyzz_0[j] + fr * tg_yy_xxyyzz_1[j] + fl1_fxn * tg_yy_xyyzz_1[j];

                    tg_xyy_xxyzzz_0[j] = pb_x * tg_yy_xxyzzz_0[j] + fr * tg_yy_xxyzzz_1[j] + fl1_fxn * tg_yy_xyzzz_1[j];

                    tg_xyy_xxzzzz_0[j] = pb_x * tg_yy_xxzzzz_0[j] + fr * tg_yy_xxzzzz_1[j] + fl1_fxn * tg_yy_xzzzz_1[j];

                    tg_xyy_xyyyyy_0[j] = pb_x * tg_yy_xyyyyy_0[j] + fr * tg_yy_xyyyyy_1[j] + 0.5 * fl1_fxn * tg_yy_yyyyy_1[j];

                    tg_xyy_xyyyyz_0[j] = pb_x * tg_yy_xyyyyz_0[j] + fr * tg_yy_xyyyyz_1[j] + 0.5 * fl1_fxn * tg_yy_yyyyz_1[j];

                    tg_xyy_xyyyzz_0[j] = pb_x * tg_yy_xyyyzz_0[j] + fr * tg_yy_xyyyzz_1[j] + 0.5 * fl1_fxn * tg_yy_yyyzz_1[j];

                    tg_xyy_xyyzzz_0[j] = pb_x * tg_yy_xyyzzz_0[j] + fr * tg_yy_xyyzzz_1[j] + 0.5 * fl1_fxn * tg_yy_yyzzz_1[j];

                    tg_xyy_xyzzzz_0[j] = pb_x * tg_yy_xyzzzz_0[j] + fr * tg_yy_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_yy_yzzzz_1[j];

                    tg_xyy_xzzzzz_0[j] = pb_x * tg_yy_xzzzzz_0[j] + fr * tg_yy_xzzzzz_1[j] + 0.5 * fl1_fxn * tg_yy_zzzzz_1[j];

                    tg_xyy_yyyyyy_0[j] = pb_x * tg_yy_yyyyyy_0[j] + fr * tg_yy_yyyyyy_1[j];

                    tg_xyy_yyyyyz_0[j] = pb_x * tg_yy_yyyyyz_0[j] + fr * tg_yy_yyyyyz_1[j];

                    tg_xyy_yyyyzz_0[j] = pb_x * tg_yy_yyyyzz_0[j] + fr * tg_yy_yyyyzz_1[j];

                    tg_xyy_yyyzzz_0[j] = pb_x * tg_yy_yyyzzz_0[j] + fr * tg_yy_yyyzzz_1[j];

                    tg_xyy_yyzzzz_0[j] = pb_x * tg_yy_yyzzzz_0[j] + fr * tg_yy_yyzzzz_1[j];

                    tg_xyy_yzzzzz_0[j] = pb_x * tg_yy_yzzzzz_0[j] + fr * tg_yy_yzzzzz_1[j];

                    tg_xyy_zzzzzz_0[j] = pb_x * tg_yy_zzzzzz_0[j] + fr * tg_yy_zzzzzz_1[j];

                    tg_xyz_xxxxxx_0[j] = pb_x * tg_yz_xxxxxx_0[j] + fr * tg_yz_xxxxxx_1[j] + 3.0 * fl1_fxn * tg_yz_xxxxx_1[j];

                    tg_xyz_xxxxxy_0[j] = pb_x * tg_yz_xxxxxy_0[j] + fr * tg_yz_xxxxxy_1[j] + 2.5 * fl1_fxn * tg_yz_xxxxy_1[j];

                    tg_xyz_xxxxxz_0[j] = pb_x * tg_yz_xxxxxz_0[j] + fr * tg_yz_xxxxxz_1[j] + 2.5 * fl1_fxn * tg_yz_xxxxz_1[j];

                    tg_xyz_xxxxyy_0[j] = pb_x * tg_yz_xxxxyy_0[j] + fr * tg_yz_xxxxyy_1[j] + 2.0 * fl1_fxn * tg_yz_xxxyy_1[j];

                    tg_xyz_xxxxyz_0[j] = pb_x * tg_yz_xxxxyz_0[j] + fr * tg_yz_xxxxyz_1[j] + 2.0 * fl1_fxn * tg_yz_xxxyz_1[j];

                    tg_xyz_xxxxzz_0[j] = pb_x * tg_yz_xxxxzz_0[j] + fr * tg_yz_xxxxzz_1[j] + 2.0 * fl1_fxn * tg_yz_xxxzz_1[j];

                    tg_xyz_xxxyyy_0[j] = pb_x * tg_yz_xxxyyy_0[j] + fr * tg_yz_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_yz_xxyyy_1[j];

                    tg_xyz_xxxyyz_0[j] = pb_x * tg_yz_xxxyyz_0[j] + fr * tg_yz_xxxyyz_1[j] + 1.5 * fl1_fxn * tg_yz_xxyyz_1[j];

                    tg_xyz_xxxyzz_0[j] = pb_x * tg_yz_xxxyzz_0[j] + fr * tg_yz_xxxyzz_1[j] + 1.5 * fl1_fxn * tg_yz_xxyzz_1[j];

                    tg_xyz_xxxzzz_0[j] = pb_x * tg_yz_xxxzzz_0[j] + fr * tg_yz_xxxzzz_1[j] + 1.5 * fl1_fxn * tg_yz_xxzzz_1[j];

                    tg_xyz_xxyyyy_0[j] = pb_x * tg_yz_xxyyyy_0[j] + fr * tg_yz_xxyyyy_1[j] + fl1_fxn * tg_yz_xyyyy_1[j];

                    tg_xyz_xxyyyz_0[j] = pb_x * tg_yz_xxyyyz_0[j] + fr * tg_yz_xxyyyz_1[j] + fl1_fxn * tg_yz_xyyyz_1[j];

                    tg_xyz_xxyyzz_0[j] = pb_x * tg_yz_xxyyzz_0[j] + fr * tg_yz_xxyyzz_1[j] + fl1_fxn * tg_yz_xyyzz_1[j];

                    tg_xyz_xxyzzz_0[j] = pb_x * tg_yz_xxyzzz_0[j] + fr * tg_yz_xxyzzz_1[j] + fl1_fxn * tg_yz_xyzzz_1[j];

                    tg_xyz_xxzzzz_0[j] = pb_x * tg_yz_xxzzzz_0[j] + fr * tg_yz_xxzzzz_1[j] + fl1_fxn * tg_yz_xzzzz_1[j];

                    tg_xyz_xyyyyy_0[j] = pb_x * tg_yz_xyyyyy_0[j] + fr * tg_yz_xyyyyy_1[j] + 0.5 * fl1_fxn * tg_yz_yyyyy_1[j];

                    tg_xyz_xyyyyz_0[j] = pb_x * tg_yz_xyyyyz_0[j] + fr * tg_yz_xyyyyz_1[j] + 0.5 * fl1_fxn * tg_yz_yyyyz_1[j];

                    tg_xyz_xyyyzz_0[j] = pb_x * tg_yz_xyyyzz_0[j] + fr * tg_yz_xyyyzz_1[j] + 0.5 * fl1_fxn * tg_yz_yyyzz_1[j];

                    tg_xyz_xyyzzz_0[j] = pb_x * tg_yz_xyyzzz_0[j] + fr * tg_yz_xyyzzz_1[j] + 0.5 * fl1_fxn * tg_yz_yyzzz_1[j];

                    tg_xyz_xyzzzz_0[j] = pb_x * tg_yz_xyzzzz_0[j] + fr * tg_yz_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_yz_yzzzz_1[j];

                    tg_xyz_xzzzzz_0[j] = pb_x * tg_yz_xzzzzz_0[j] + fr * tg_yz_xzzzzz_1[j] + 0.5 * fl1_fxn * tg_yz_zzzzz_1[j];

                    tg_xyz_yyyyyy_0[j] = pb_x * tg_yz_yyyyyy_0[j] + fr * tg_yz_yyyyyy_1[j];

                    tg_xyz_yyyyyz_0[j] = pb_x * tg_yz_yyyyyz_0[j] + fr * tg_yz_yyyyyz_1[j];

                    tg_xyz_yyyyzz_0[j] = pb_x * tg_yz_yyyyzz_0[j] + fr * tg_yz_yyyyzz_1[j];

                    tg_xyz_yyyzzz_0[j] = pb_x * tg_yz_yyyzzz_0[j] + fr * tg_yz_yyyzzz_1[j];

                    tg_xyz_yyzzzz_0[j] = pb_x * tg_yz_yyzzzz_0[j] + fr * tg_yz_yyzzzz_1[j];

                    tg_xyz_yzzzzz_0[j] = pb_x * tg_yz_yzzzzz_0[j] + fr * tg_yz_yzzzzz_1[j];

                    tg_xyz_zzzzzz_0[j] = pb_x * tg_yz_zzzzzz_0[j] + fr * tg_yz_zzzzzz_1[j];

                    tg_xzz_xxxxxx_0[j] = pb_x * tg_zz_xxxxxx_0[j] + fr * tg_zz_xxxxxx_1[j] + 3.0 * fl1_fxn * tg_zz_xxxxx_1[j];

                    tg_xzz_xxxxxy_0[j] = pb_x * tg_zz_xxxxxy_0[j] + fr * tg_zz_xxxxxy_1[j] + 2.5 * fl1_fxn * tg_zz_xxxxy_1[j];

                    tg_xzz_xxxxxz_0[j] = pb_x * tg_zz_xxxxxz_0[j] + fr * tg_zz_xxxxxz_1[j] + 2.5 * fl1_fxn * tg_zz_xxxxz_1[j];

                    tg_xzz_xxxxyy_0[j] = pb_x * tg_zz_xxxxyy_0[j] + fr * tg_zz_xxxxyy_1[j] + 2.0 * fl1_fxn * tg_zz_xxxyy_1[j];

                    tg_xzz_xxxxyz_0[j] = pb_x * tg_zz_xxxxyz_0[j] + fr * tg_zz_xxxxyz_1[j] + 2.0 * fl1_fxn * tg_zz_xxxyz_1[j];

                    tg_xzz_xxxxzz_0[j] = pb_x * tg_zz_xxxxzz_0[j] + fr * tg_zz_xxxxzz_1[j] + 2.0 * fl1_fxn * tg_zz_xxxzz_1[j];

                    tg_xzz_xxxyyy_0[j] = pb_x * tg_zz_xxxyyy_0[j] + fr * tg_zz_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_zz_xxyyy_1[j];

                    tg_xzz_xxxyyz_0[j] = pb_x * tg_zz_xxxyyz_0[j] + fr * tg_zz_xxxyyz_1[j] + 1.5 * fl1_fxn * tg_zz_xxyyz_1[j];

                    tg_xzz_xxxyzz_0[j] = pb_x * tg_zz_xxxyzz_0[j] + fr * tg_zz_xxxyzz_1[j] + 1.5 * fl1_fxn * tg_zz_xxyzz_1[j];

                    tg_xzz_xxxzzz_0[j] = pb_x * tg_zz_xxxzzz_0[j] + fr * tg_zz_xxxzzz_1[j] + 1.5 * fl1_fxn * tg_zz_xxzzz_1[j];

                    tg_xzz_xxyyyy_0[j] = pb_x * tg_zz_xxyyyy_0[j] + fr * tg_zz_xxyyyy_1[j] + fl1_fxn * tg_zz_xyyyy_1[j];

                    tg_xzz_xxyyyz_0[j] = pb_x * tg_zz_xxyyyz_0[j] + fr * tg_zz_xxyyyz_1[j] + fl1_fxn * tg_zz_xyyyz_1[j];

                    tg_xzz_xxyyzz_0[j] = pb_x * tg_zz_xxyyzz_0[j] + fr * tg_zz_xxyyzz_1[j] + fl1_fxn * tg_zz_xyyzz_1[j];

                    tg_xzz_xxyzzz_0[j] = pb_x * tg_zz_xxyzzz_0[j] + fr * tg_zz_xxyzzz_1[j] + fl1_fxn * tg_zz_xyzzz_1[j];

                    tg_xzz_xxzzzz_0[j] = pb_x * tg_zz_xxzzzz_0[j] + fr * tg_zz_xxzzzz_1[j] + fl1_fxn * tg_zz_xzzzz_1[j];

                    tg_xzz_xyyyyy_0[j] = pb_x * tg_zz_xyyyyy_0[j] + fr * tg_zz_xyyyyy_1[j] + 0.5 * fl1_fxn * tg_zz_yyyyy_1[j];

                    tg_xzz_xyyyyz_0[j] = pb_x * tg_zz_xyyyyz_0[j] + fr * tg_zz_xyyyyz_1[j] + 0.5 * fl1_fxn * tg_zz_yyyyz_1[j];

                    tg_xzz_xyyyzz_0[j] = pb_x * tg_zz_xyyyzz_0[j] + fr * tg_zz_xyyyzz_1[j] + 0.5 * fl1_fxn * tg_zz_yyyzz_1[j];

                    tg_xzz_xyyzzz_0[j] = pb_x * tg_zz_xyyzzz_0[j] + fr * tg_zz_xyyzzz_1[j] + 0.5 * fl1_fxn * tg_zz_yyzzz_1[j];

                    tg_xzz_xyzzzz_0[j] = pb_x * tg_zz_xyzzzz_0[j] + fr * tg_zz_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_zz_yzzzz_1[j];

                    tg_xzz_xzzzzz_0[j] = pb_x * tg_zz_xzzzzz_0[j] + fr * tg_zz_xzzzzz_1[j] + 0.5 * fl1_fxn * tg_zz_zzzzz_1[j];

                    tg_xzz_yyyyyy_0[j] = pb_x * tg_zz_yyyyyy_0[j] + fr * tg_zz_yyyyyy_1[j];

                    tg_xzz_yyyyyz_0[j] = pb_x * tg_zz_yyyyyz_0[j] + fr * tg_zz_yyyyyz_1[j];

                    tg_xzz_yyyyzz_0[j] = pb_x * tg_zz_yyyyzz_0[j] + fr * tg_zz_yyyyzz_1[j];

                    tg_xzz_yyyzzz_0[j] = pb_x * tg_zz_yyyzzz_0[j] + fr * tg_zz_yyyzzz_1[j];

                    tg_xzz_yyzzzz_0[j] = pb_x * tg_zz_yyzzzz_0[j] + fr * tg_zz_yyzzzz_1[j];

                    tg_xzz_yzzzzz_0[j] = pb_x * tg_zz_yzzzzz_0[j] + fr * tg_zz_yzzzzz_1[j];

                    tg_xzz_zzzzzz_0[j] = pb_x * tg_zz_zzzzzz_0[j] + fr * tg_zz_zzzzzz_1[j];

                    fr = wp_y[j]; 

                    tg_yyy_xxxxxx_0[j] = pb_y * tg_yy_xxxxxx_0[j] + fr * tg_yy_xxxxxx_1[j] + fl1_fx * (tg_y_xxxxxx_0[j] - tg_y_xxxxxx_1[j] * fl1_fza);

                    tg_yyy_xxxxxy_0[j] = pb_y * tg_yy_xxxxxy_0[j] + fr * tg_yy_xxxxxy_1[j] + fl1_fx * (tg_y_xxxxxy_0[j] - tg_y_xxxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yy_xxxxx_1[j];

                    tg_yyy_xxxxxz_0[j] = pb_y * tg_yy_xxxxxz_0[j] + fr * tg_yy_xxxxxz_1[j] + fl1_fx * (tg_y_xxxxxz_0[j] - tg_y_xxxxxz_1[j] * fl1_fza);

                    tg_yyy_xxxxyy_0[j] = pb_y * tg_yy_xxxxyy_0[j] + fr * tg_yy_xxxxyy_1[j] + fl1_fx * (tg_y_xxxxyy_0[j] - tg_y_xxxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yy_xxxxy_1[j];

                    tg_yyy_xxxxyz_0[j] = pb_y * tg_yy_xxxxyz_0[j] + fr * tg_yy_xxxxyz_1[j] + fl1_fx * (tg_y_xxxxyz_0[j] - tg_y_xxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yy_xxxxz_1[j];

                    tg_yyy_xxxxzz_0[j] = pb_y * tg_yy_xxxxzz_0[j] + fr * tg_yy_xxxxzz_1[j] + fl1_fx * (tg_y_xxxxzz_0[j] - tg_y_xxxxzz_1[j] * fl1_fza);

                    tg_yyy_xxxyyy_0[j] = pb_y * tg_yy_xxxyyy_0[j] + fr * tg_yy_xxxyyy_1[j] + fl1_fx * (tg_y_xxxyyy_0[j] - tg_y_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yy_xxxyy_1[j];

                    tg_yyy_xxxyyz_0[j] = pb_y * tg_yy_xxxyyz_0[j] + fr * tg_yy_xxxyyz_1[j] + fl1_fx * (tg_y_xxxyyz_0[j] - tg_y_xxxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yy_xxxyz_1[j];

                    tg_yyy_xxxyzz_0[j] = pb_y * tg_yy_xxxyzz_0[j] + fr * tg_yy_xxxyzz_1[j] + fl1_fx * (tg_y_xxxyzz_0[j] - tg_y_xxxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yy_xxxzz_1[j];

                    tg_yyy_xxxzzz_0[j] = pb_y * tg_yy_xxxzzz_0[j] + fr * tg_yy_xxxzzz_1[j] + fl1_fx * (tg_y_xxxzzz_0[j] - tg_y_xxxzzz_1[j] * fl1_fza);

                    tg_yyy_xxyyyy_0[j] = pb_y * tg_yy_xxyyyy_0[j] + fr * tg_yy_xxyyyy_1[j] + fl1_fx * (tg_y_xxyyyy_0[j] - tg_y_xxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yy_xxyyy_1[j];

                    tg_yyy_xxyyyz_0[j] = pb_y * tg_yy_xxyyyz_0[j] + fr * tg_yy_xxyyyz_1[j] + fl1_fx * (tg_y_xxyyyz_0[j] - tg_y_xxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yy_xxyyz_1[j];

                    tg_yyy_xxyyzz_0[j] = pb_y * tg_yy_xxyyzz_0[j] + fr * tg_yy_xxyyzz_1[j] + fl1_fx * (tg_y_xxyyzz_0[j] - tg_y_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yy_xxyzz_1[j];

                    tg_yyy_xxyzzz_0[j] = pb_y * tg_yy_xxyzzz_0[j] + fr * tg_yy_xxyzzz_1[j] + fl1_fx * (tg_y_xxyzzz_0[j] - tg_y_xxyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yy_xxzzz_1[j];

                    tg_yyy_xxzzzz_0[j] = pb_y * tg_yy_xxzzzz_0[j] + fr * tg_yy_xxzzzz_1[j] + fl1_fx * (tg_y_xxzzzz_0[j] - tg_y_xxzzzz_1[j] * fl1_fza);

                    tg_yyy_xyyyyy_0[j] = pb_y * tg_yy_xyyyyy_0[j] + fr * tg_yy_xyyyyy_1[j] + fl1_fx * (tg_y_xyyyyy_0[j] - tg_y_xyyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yy_xyyyy_1[j];

                    tg_yyy_xyyyyz_0[j] = pb_y * tg_yy_xyyyyz_0[j] + fr * tg_yy_xyyyyz_1[j] + fl1_fx * (tg_y_xyyyyz_0[j] - tg_y_xyyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yy_xyyyz_1[j];

                    tg_yyy_xyyyzz_0[j] = pb_y * tg_yy_xyyyzz_0[j] + fr * tg_yy_xyyyzz_1[j] + fl1_fx * (tg_y_xyyyzz_0[j] - tg_y_xyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yy_xyyzz_1[j];

                    tg_yyy_xyyzzz_0[j] = pb_y * tg_yy_xyyzzz_0[j] + fr * tg_yy_xyyzzz_1[j] + fl1_fx * (tg_y_xyyzzz_0[j] - tg_y_xyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yy_xyzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSFSI_187_280(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (187,280)

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
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_3_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_3_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_2_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_2_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_1_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {1, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_1_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {1, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_2_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {5, -1, -1, -1}, 
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

                auto tg_yy_xyzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 103); 

                auto tg_yy_xzzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 104); 

                auto tg_yy_yyyyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 105); 

                auto tg_yy_yyyyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 106); 

                auto tg_yy_yyyyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 107); 

                auto tg_yy_yyyzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 108); 

                auto tg_yy_yyzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 109); 

                auto tg_yy_yzzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 110); 

                auto tg_yy_zzzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 111); 

                auto tg_yz_xxxxxx_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 112); 

                auto tg_yz_xxxxxy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 113); 

                auto tg_yz_xxxxxz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 114); 

                auto tg_yz_xxxxyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 115); 

                auto tg_yz_xxxxyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 116); 

                auto tg_yz_xxxxzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 117); 

                auto tg_yz_xxxyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 118); 

                auto tg_yz_xxxyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 119); 

                auto tg_yz_xxxyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 120); 

                auto tg_yz_xxxzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 121); 

                auto tg_yz_xxyyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 122); 

                auto tg_yz_xxyyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 123); 

                auto tg_yz_xxyyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 124); 

                auto tg_yz_xxyzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 125); 

                auto tg_yz_xxzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 126); 

                auto tg_yz_xyyyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 127); 

                auto tg_yz_xyyyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 128); 

                auto tg_yz_xyyyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 129); 

                auto tg_yz_xyyzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 130); 

                auto tg_yz_xyzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 131); 

                auto tg_yz_xzzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 132); 

                auto tg_yz_yyyyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 133); 

                auto tg_yz_yyyyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 134); 

                auto tg_yz_yyyyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 135); 

                auto tg_yz_yyyzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 136); 

                auto tg_yz_yyzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 137); 

                auto tg_yz_yzzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 138); 

                auto tg_yz_zzzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 139); 

                auto tg_zz_xxxxxx_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 140); 

                auto tg_zz_xxxxxy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 141); 

                auto tg_zz_xxxxxz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 142); 

                auto tg_zz_xxxxyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 143); 

                auto tg_zz_xxxxyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 144); 

                auto tg_zz_xxxxzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 145); 

                auto tg_zz_xxxyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 146); 

                auto tg_zz_xxxyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 147); 

                auto tg_zz_xxxyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 148); 

                auto tg_zz_xxxzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 149); 

                auto tg_zz_xxyyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 150); 

                auto tg_zz_xxyyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 151); 

                auto tg_zz_xxyyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 152); 

                auto tg_zz_xxyzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 153); 

                auto tg_zz_xxzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 154); 

                auto tg_zz_xyyyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 155); 

                auto tg_zz_xyyyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 156); 

                auto tg_zz_xyyyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 157); 

                auto tg_zz_xyyzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 158); 

                auto tg_zz_xyzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 159); 

                auto tg_zz_xzzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 160); 

                auto tg_zz_yyyyyy_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 161); 

                auto tg_zz_yyyyyz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 162); 

                auto tg_zz_yyyyzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 163); 

                auto tg_zz_yyyzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 164); 

                auto tg_zz_yyzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 165); 

                auto tg_zz_yzzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 166); 

                auto tg_zz_zzzzzz_0 = primBuffer[pidx_g_2_6_m0].data(168 * idx + 167); 

                auto tg_yy_xyzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 103); 

                auto tg_yy_xzzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 104); 

                auto tg_yy_yyyyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 105); 

                auto tg_yy_yyyyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 106); 

                auto tg_yy_yyyyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 107); 

                auto tg_yy_yyyzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 108); 

                auto tg_yy_yyzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 109); 

                auto tg_yy_yzzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 110); 

                auto tg_yy_zzzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 111); 

                auto tg_yz_xxxxxx_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 112); 

                auto tg_yz_xxxxxy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 113); 

                auto tg_yz_xxxxxz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 114); 

                auto tg_yz_xxxxyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 115); 

                auto tg_yz_xxxxyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 116); 

                auto tg_yz_xxxxzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 117); 

                auto tg_yz_xxxyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 118); 

                auto tg_yz_xxxyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 119); 

                auto tg_yz_xxxyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 120); 

                auto tg_yz_xxxzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 121); 

                auto tg_yz_xxyyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 122); 

                auto tg_yz_xxyyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 123); 

                auto tg_yz_xxyyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 124); 

                auto tg_yz_xxyzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 125); 

                auto tg_yz_xxzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 126); 

                auto tg_yz_xyyyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 127); 

                auto tg_yz_xyyyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 128); 

                auto tg_yz_xyyyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 129); 

                auto tg_yz_xyyzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 130); 

                auto tg_yz_xyzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 131); 

                auto tg_yz_xzzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 132); 

                auto tg_yz_yyyyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 133); 

                auto tg_yz_yyyyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 134); 

                auto tg_yz_yyyyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 135); 

                auto tg_yz_yyyzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 136); 

                auto tg_yz_yyzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 137); 

                auto tg_yz_yzzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 138); 

                auto tg_yz_zzzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 139); 

                auto tg_zz_xxxxxx_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 140); 

                auto tg_zz_xxxxxy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 141); 

                auto tg_zz_xxxxxz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 142); 

                auto tg_zz_xxxxyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 143); 

                auto tg_zz_xxxxyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 144); 

                auto tg_zz_xxxxzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 145); 

                auto tg_zz_xxxyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 146); 

                auto tg_zz_xxxyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 147); 

                auto tg_zz_xxxyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 148); 

                auto tg_zz_xxxzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 149); 

                auto tg_zz_xxyyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 150); 

                auto tg_zz_xxyyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 151); 

                auto tg_zz_xxyyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 152); 

                auto tg_zz_xxyzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 153); 

                auto tg_zz_xxzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 154); 

                auto tg_zz_xyyyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 155); 

                auto tg_zz_xyyyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 156); 

                auto tg_zz_xyyyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 157); 

                auto tg_zz_xyyzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 158); 

                auto tg_zz_xyzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 159); 

                auto tg_zz_xzzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 160); 

                auto tg_zz_yyyyyy_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 161); 

                auto tg_zz_yyyyyz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 162); 

                auto tg_zz_yyyyzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 163); 

                auto tg_zz_yyyzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 164); 

                auto tg_zz_yyzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 165); 

                auto tg_zz_yzzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 166); 

                auto tg_zz_zzzzzz_1 = primBuffer[pidx_g_2_6_m1].data(168 * idx + 167); 

                auto tg_y_xyzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 47); 

                auto tg_y_xzzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 48); 

                auto tg_y_yyyyyy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 49); 

                auto tg_y_yyyyyz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 50); 

                auto tg_y_yyyyzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 51); 

                auto tg_y_yyyzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 52); 

                auto tg_y_yyzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 53); 

                auto tg_y_yzzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 54); 

                auto tg_y_zzzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 55); 

                auto tg_z_xxxxxx_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 56); 

                auto tg_z_xxxxxy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 57); 

                auto tg_z_xxxxxz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 58); 

                auto tg_z_xxxxyy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 59); 

                auto tg_z_xxxxyz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 60); 

                auto tg_z_xxxxzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 61); 

                auto tg_z_xxxyyy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 62); 

                auto tg_z_xxxyyz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 63); 

                auto tg_z_xxxyzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 64); 

                auto tg_z_xxxzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 65); 

                auto tg_z_xxyyyy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 66); 

                auto tg_z_xxyyyz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 67); 

                auto tg_z_xxyyzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 68); 

                auto tg_z_xxyzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 69); 

                auto tg_z_xxzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 70); 

                auto tg_z_xyyyyy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 71); 

                auto tg_z_xyyyyz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 72); 

                auto tg_z_xyyyzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 73); 

                auto tg_z_xyyzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 74); 

                auto tg_z_xyzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 75); 

                auto tg_z_xzzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 76); 

                auto tg_z_yyyyyy_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 77); 

                auto tg_z_yyyyyz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 78); 

                auto tg_z_yyyyzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 79); 

                auto tg_z_yyyzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 80); 

                auto tg_z_yyzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 81); 

                auto tg_z_yzzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 82); 

                auto tg_z_zzzzzz_0 = primBuffer[pidx_g_1_6_m0].data(84 * idx + 83); 

                auto tg_y_xyzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 47); 

                auto tg_y_xzzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 48); 

                auto tg_y_yyyyyy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 49); 

                auto tg_y_yyyyyz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 50); 

                auto tg_y_yyyyzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 51); 

                auto tg_y_yyyzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 52); 

                auto tg_y_yyzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 53); 

                auto tg_y_yzzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 54); 

                auto tg_y_zzzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 55); 

                auto tg_z_xxxxxx_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 56); 

                auto tg_z_xxxxxy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 57); 

                auto tg_z_xxxxxz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 58); 

                auto tg_z_xxxxyy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 59); 

                auto tg_z_xxxxyz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 60); 

                auto tg_z_xxxxzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 61); 

                auto tg_z_xxxyyy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 62); 

                auto tg_z_xxxyyz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 63); 

                auto tg_z_xxxyzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 64); 

                auto tg_z_xxxzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 65); 

                auto tg_z_xxyyyy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 66); 

                auto tg_z_xxyyyz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 67); 

                auto tg_z_xxyyzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 68); 

                auto tg_z_xxyzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 69); 

                auto tg_z_xxzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 70); 

                auto tg_z_xyyyyy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 71); 

                auto tg_z_xyyyyz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 72); 

                auto tg_z_xyyyzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 73); 

                auto tg_z_xyyzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 74); 

                auto tg_z_xyzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 75); 

                auto tg_z_xzzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 76); 

                auto tg_z_yyyyyy_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 77); 

                auto tg_z_yyyyyz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 78); 

                auto tg_z_yyyyzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 79); 

                auto tg_z_yyyzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 80); 

                auto tg_z_yyzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 81); 

                auto tg_z_yzzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 82); 

                auto tg_z_zzzzzz_1 = primBuffer[pidx_g_1_6_m1].data(84 * idx + 83); 

                auto tg_yy_xzzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 77); 

                auto tg_yy_yyyyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 78); 

                auto tg_yy_yyyyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 79); 

                auto tg_yy_yyyzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 80); 

                auto tg_yy_yyzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 81); 

                auto tg_yy_yzzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 82); 

                auto tg_yy_zzzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 83); 

                auto tg_yz_xxxxx_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 84); 

                auto tg_yz_xxxxy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 85); 

                auto tg_yz_xxxxz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 86); 

                auto tg_yz_xxxyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 87); 

                auto tg_yz_xxxyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 88); 

                auto tg_yz_xxxzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 89); 

                auto tg_yz_xxyyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 90); 

                auto tg_yz_xxyyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 91); 

                auto tg_yz_xxyzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 92); 

                auto tg_yz_xxzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 93); 

                auto tg_yz_xyyyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 94); 

                auto tg_yz_xyyyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 95); 

                auto tg_yz_xyyzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 96); 

                auto tg_yz_xyzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 97); 

                auto tg_yz_xzzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 98); 

                auto tg_yz_yyyyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 99); 

                auto tg_yz_yyyyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 100); 

                auto tg_yz_yyyzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 101); 

                auto tg_yz_yyzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 102); 

                auto tg_yz_yzzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 103); 

                auto tg_yz_zzzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 104); 

                auto tg_zz_xxxxx_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 105); 

                auto tg_zz_xxxxy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 106); 

                auto tg_zz_xxxxz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 107); 

                auto tg_zz_xxxyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 108); 

                auto tg_zz_xxxyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 109); 

                auto tg_zz_xxxzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 110); 

                auto tg_zz_xxyyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 111); 

                auto tg_zz_xxyyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 112); 

                auto tg_zz_xxyzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 113); 

                auto tg_zz_xxzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 114); 

                auto tg_zz_xyyyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 115); 

                auto tg_zz_xyyyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 116); 

                auto tg_zz_xyyzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 117); 

                auto tg_zz_xyzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 118); 

                auto tg_zz_xzzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 119); 

                auto tg_zz_yyyyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 120); 

                auto tg_zz_yyyyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 121); 

                auto tg_zz_yyyzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 122); 

                auto tg_zz_yyzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 123); 

                auto tg_zz_yzzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 124); 

                auto tg_zz_zzzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 125); 

                // set up pointers to integrals

                auto tg_yyy_xyzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 187); 

                auto tg_yyy_xzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 188); 

                auto tg_yyy_yyyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 189); 

                auto tg_yyy_yyyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 190); 

                auto tg_yyy_yyyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 191); 

                auto tg_yyy_yyyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 192); 

                auto tg_yyy_yyzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 193); 

                auto tg_yyy_yzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 194); 

                auto tg_yyy_zzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 195); 

                auto tg_yyz_xxxxxx_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 196); 

                auto tg_yyz_xxxxxy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 197); 

                auto tg_yyz_xxxxxz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 198); 

                auto tg_yyz_xxxxyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 199); 

                auto tg_yyz_xxxxyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 200); 

                auto tg_yyz_xxxxzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 201); 

                auto tg_yyz_xxxyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 202); 

                auto tg_yyz_xxxyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 203); 

                auto tg_yyz_xxxyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 204); 

                auto tg_yyz_xxxzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 205); 

                auto tg_yyz_xxyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 206); 

                auto tg_yyz_xxyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 207); 

                auto tg_yyz_xxyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 208); 

                auto tg_yyz_xxyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 209); 

                auto tg_yyz_xxzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 210); 

                auto tg_yyz_xyyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 211); 

                auto tg_yyz_xyyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 212); 

                auto tg_yyz_xyyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 213); 

                auto tg_yyz_xyyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 214); 

                auto tg_yyz_xyzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 215); 

                auto tg_yyz_xzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 216); 

                auto tg_yyz_yyyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 217); 

                auto tg_yyz_yyyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 218); 

                auto tg_yyz_yyyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 219); 

                auto tg_yyz_yyyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 220); 

                auto tg_yyz_yyzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 221); 

                auto tg_yyz_yzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 222); 

                auto tg_yyz_zzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 223); 

                auto tg_yzz_xxxxxx_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 224); 

                auto tg_yzz_xxxxxy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 225); 

                auto tg_yzz_xxxxxz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 226); 

                auto tg_yzz_xxxxyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 227); 

                auto tg_yzz_xxxxyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 228); 

                auto tg_yzz_xxxxzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 229); 

                auto tg_yzz_xxxyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 230); 

                auto tg_yzz_xxxyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 231); 

                auto tg_yzz_xxxyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 232); 

                auto tg_yzz_xxxzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 233); 

                auto tg_yzz_xxyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 234); 

                auto tg_yzz_xxyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 235); 

                auto tg_yzz_xxyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 236); 

                auto tg_yzz_xxyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 237); 

                auto tg_yzz_xxzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 238); 

                auto tg_yzz_xyyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 239); 

                auto tg_yzz_xyyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 240); 

                auto tg_yzz_xyyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 241); 

                auto tg_yzz_xyyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 242); 

                auto tg_yzz_xyzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 243); 

                auto tg_yzz_xzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 244); 

                auto tg_yzz_yyyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 245); 

                auto tg_yzz_yyyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 246); 

                auto tg_yzz_yyyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 247); 

                auto tg_yzz_yyyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 248); 

                auto tg_yzz_yyzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 249); 

                auto tg_yzz_yzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 250); 

                auto tg_yzz_zzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 251); 

                auto tg_zzz_xxxxxx_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 252); 

                auto tg_zzz_xxxxxy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 253); 

                auto tg_zzz_xxxxxz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 254); 

                auto tg_zzz_xxxxyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 255); 

                auto tg_zzz_xxxxyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 256); 

                auto tg_zzz_xxxxzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 257); 

                auto tg_zzz_xxxyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 258); 

                auto tg_zzz_xxxyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 259); 

                auto tg_zzz_xxxyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 260); 

                auto tg_zzz_xxxzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 261); 

                auto tg_zzz_xxyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 262); 

                auto tg_zzz_xxyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 263); 

                auto tg_zzz_xxyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 264); 

                auto tg_zzz_xxyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 265); 

                auto tg_zzz_xxzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 266); 

                auto tg_zzz_xyyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 267); 

                auto tg_zzz_xyyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 268); 

                auto tg_zzz_xyyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 269); 

                auto tg_zzz_xyyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 270); 

                auto tg_zzz_xyzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 271); 

                auto tg_zzz_xzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 272); 

                auto tg_zzz_yyyyyy_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 273); 

                auto tg_zzz_yyyyyz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 274); 

                auto tg_zzz_yyyyzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 275); 

                auto tg_zzz_yyyzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 276); 

                auto tg_zzz_yyzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 277); 

                auto tg_zzz_yzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 278); 

                auto tg_zzz_zzzzzz_0 = primBuffer[pidx_g_3_6_m0].data(280 * idx + 279); 

                // Batch of Integrals (187,280)

                #pragma omp simd aligned(fxn, fza, tg_y_xyzzzz_0, tg_y_xyzzzz_1, tg_y_xzzzzz_0, tg_y_xzzzzz_1, \
                                         tg_y_yyyyyy_0, tg_y_yyyyyy_1, tg_y_yyyyyz_0, tg_y_yyyyyz_1, tg_y_yyyyzz_0, \
                                         tg_y_yyyyzz_1, tg_y_yyyzzz_0, tg_y_yyyzzz_1, tg_y_yyzzzz_0, tg_y_yyzzzz_1, \
                                         tg_y_yzzzzz_0, tg_y_yzzzzz_1, tg_y_zzzzzz_0, tg_y_zzzzzz_1, tg_yy_xyzzzz_0, \
                                         tg_yy_xyzzzz_1, tg_yy_xzzzz_1, tg_yy_xzzzzz_0, tg_yy_xzzzzz_1, tg_yy_yyyyy_1, \
                                         tg_yy_yyyyyy_0, tg_yy_yyyyyy_1, tg_yy_yyyyyz_0, tg_yy_yyyyyz_1, tg_yy_yyyyz_1, \
                                         tg_yy_yyyyzz_0, tg_yy_yyyyzz_1, tg_yy_yyyzz_1, tg_yy_yyyzzz_0, tg_yy_yyyzzz_1, \
                                         tg_yy_yyzzz_1, tg_yy_yyzzzz_0, tg_yy_yyzzzz_1, tg_yy_yzzzz_1, tg_yy_yzzzzz_0, \
                                         tg_yy_yzzzzz_1, tg_yy_zzzzz_1, tg_yy_zzzzzz_0, tg_yy_zzzzzz_1, tg_yyy_xyzzzz_0, \
                                         tg_yyy_xzzzzz_0, tg_yyy_yyyyyy_0, tg_yyy_yyyyyz_0, tg_yyy_yyyyzz_0, tg_yyy_yyyzzz_0, \
                                         tg_yyy_yyzzzz_0, tg_yyy_yzzzzz_0, tg_yyy_zzzzzz_0, tg_yyz_xxxxxx_0, tg_yyz_xxxxxy_0, \
                                         tg_yyz_xxxxxz_0, tg_yyz_xxxxyy_0, tg_yyz_xxxxyz_0, tg_yyz_xxxxzz_0, tg_yyz_xxxyyy_0, \
                                         tg_yyz_xxxyyz_0, tg_yyz_xxxyzz_0, tg_yyz_xxxzzz_0, tg_yyz_xxyyyy_0, tg_yyz_xxyyyz_0, \
                                         tg_yyz_xxyyzz_0, tg_yyz_xxyzzz_0, tg_yyz_xxzzzz_0, tg_yyz_xyyyyy_0, tg_yyz_xyyyyz_0, \
                                         tg_yyz_xyyyzz_0, tg_yyz_xyyzzz_0, tg_yyz_xyzzzz_0, tg_yyz_xzzzzz_0, tg_yyz_yyyyyy_0, \
                                         tg_yyz_yyyyyz_0, tg_yyz_yyyyzz_0, tg_yyz_yyyzzz_0, tg_yyz_yyzzzz_0, tg_yyz_yzzzzz_0, \
                                         tg_yyz_zzzzzz_0, tg_yz_xxxxx_1, tg_yz_xxxxxx_0, tg_yz_xxxxxx_1, tg_yz_xxxxxy_0, \
                                         tg_yz_xxxxxy_1, tg_yz_xxxxxz_0, tg_yz_xxxxxz_1, tg_yz_xxxxy_1, tg_yz_xxxxyy_0, \
                                         tg_yz_xxxxyy_1, tg_yz_xxxxyz_0, tg_yz_xxxxyz_1, tg_yz_xxxxz_1, tg_yz_xxxxzz_0, \
                                         tg_yz_xxxxzz_1, tg_yz_xxxyy_1, tg_yz_xxxyyy_0, tg_yz_xxxyyy_1, tg_yz_xxxyyz_0, \
                                         tg_yz_xxxyyz_1, tg_yz_xxxyz_1, tg_yz_xxxyzz_0, tg_yz_xxxyzz_1, tg_yz_xxxzz_1, \
                                         tg_yz_xxxzzz_0, tg_yz_xxxzzz_1, tg_yz_xxyyy_1, tg_yz_xxyyyy_0, tg_yz_xxyyyy_1, \
                                         tg_yz_xxyyyz_0, tg_yz_xxyyyz_1, tg_yz_xxyyz_1, tg_yz_xxyyzz_0, tg_yz_xxyyzz_1, \
                                         tg_yz_xxyzz_1, tg_yz_xxyzzz_0, tg_yz_xxyzzz_1, tg_yz_xxzzz_1, tg_yz_xxzzzz_0, \
                                         tg_yz_xxzzzz_1, tg_yz_xyyyy_1, tg_yz_xyyyyy_0, tg_yz_xyyyyy_1, tg_yz_xyyyyz_0, \
                                         tg_yz_xyyyyz_1, tg_yz_xyyyz_1, tg_yz_xyyyzz_0, tg_yz_xyyyzz_1, tg_yz_xyyzz_1, \
                                         tg_yz_xyyzzz_0, tg_yz_xyyzzz_1, tg_yz_xyzzz_1, tg_yz_xyzzzz_0, tg_yz_xyzzzz_1, \
                                         tg_yz_xzzzz_1, tg_yz_xzzzzz_0, tg_yz_xzzzzz_1, tg_yz_yyyyy_1, tg_yz_yyyyyy_0, \
                                         tg_yz_yyyyyy_1, tg_yz_yyyyyz_0, tg_yz_yyyyyz_1, tg_yz_yyyyz_1, tg_yz_yyyyzz_0, \
                                         tg_yz_yyyyzz_1, tg_yz_yyyzz_1, tg_yz_yyyzzz_0, tg_yz_yyyzzz_1, tg_yz_yyzzz_1, \
                                         tg_yz_yyzzzz_0, tg_yz_yyzzzz_1, tg_yz_yzzzz_1, tg_yz_yzzzzz_0, tg_yz_yzzzzz_1, \
                                         tg_yz_zzzzz_1, tg_yz_zzzzzz_0, tg_yz_zzzzzz_1, tg_yzz_xxxxxx_0, tg_yzz_xxxxxy_0, \
                                         tg_yzz_xxxxxz_0, tg_yzz_xxxxyy_0, tg_yzz_xxxxyz_0, tg_yzz_xxxxzz_0, tg_yzz_xxxyyy_0, \
                                         tg_yzz_xxxyyz_0, tg_yzz_xxxyzz_0, tg_yzz_xxxzzz_0, tg_yzz_xxyyyy_0, tg_yzz_xxyyyz_0, \
                                         tg_yzz_xxyyzz_0, tg_yzz_xxyzzz_0, tg_yzz_xxzzzz_0, tg_yzz_xyyyyy_0, tg_yzz_xyyyyz_0, \
                                         tg_yzz_xyyyzz_0, tg_yzz_xyyzzz_0, tg_yzz_xyzzzz_0, tg_yzz_xzzzzz_0, tg_yzz_yyyyyy_0, \
                                         tg_yzz_yyyyyz_0, tg_yzz_yyyyzz_0, tg_yzz_yyyzzz_0, tg_yzz_yyzzzz_0, tg_yzz_yzzzzz_0, \
                                         tg_yzz_zzzzzz_0, tg_z_xxxxxx_0, tg_z_xxxxxx_1, tg_z_xxxxxy_0, tg_z_xxxxxy_1, \
                                         tg_z_xxxxxz_0, tg_z_xxxxxz_1, tg_z_xxxxyy_0, tg_z_xxxxyy_1, tg_z_xxxxyz_0, \
                                         tg_z_xxxxyz_1, tg_z_xxxxzz_0, tg_z_xxxxzz_1, tg_z_xxxyyy_0, tg_z_xxxyyy_1, \
                                         tg_z_xxxyyz_0, tg_z_xxxyyz_1, tg_z_xxxyzz_0, tg_z_xxxyzz_1, tg_z_xxxzzz_0, \
                                         tg_z_xxxzzz_1, tg_z_xxyyyy_0, tg_z_xxyyyy_1, tg_z_xxyyyz_0, tg_z_xxyyyz_1, \
                                         tg_z_xxyyzz_0, tg_z_xxyyzz_1, tg_z_xxyzzz_0, tg_z_xxyzzz_1, tg_z_xxzzzz_0, \
                                         tg_z_xxzzzz_1, tg_z_xyyyyy_0, tg_z_xyyyyy_1, tg_z_xyyyyz_0, tg_z_xyyyyz_1, \
                                         tg_z_xyyyzz_0, tg_z_xyyyzz_1, tg_z_xyyzzz_0, tg_z_xyyzzz_1, tg_z_xyzzzz_0, \
                                         tg_z_xyzzzz_1, tg_z_xzzzzz_0, tg_z_xzzzzz_1, tg_z_yyyyyy_0, tg_z_yyyyyy_1, \
                                         tg_z_yyyyyz_0, tg_z_yyyyyz_1, tg_z_yyyyzz_0, tg_z_yyyyzz_1, tg_z_yyyzzz_0, \
                                         tg_z_yyyzzz_1, tg_z_yyzzzz_0, tg_z_yyzzzz_1, tg_z_yzzzzz_0, tg_z_yzzzzz_1, \
                                         tg_z_zzzzzz_0, tg_z_zzzzzz_1, tg_zz_xxxxx_1, tg_zz_xxxxxx_0, tg_zz_xxxxxx_1, \
                                         tg_zz_xxxxxy_0, tg_zz_xxxxxy_1, tg_zz_xxxxxz_0, tg_zz_xxxxxz_1, tg_zz_xxxxy_1, \
                                         tg_zz_xxxxyy_0, tg_zz_xxxxyy_1, tg_zz_xxxxyz_0, tg_zz_xxxxyz_1, tg_zz_xxxxz_1, \
                                         tg_zz_xxxxzz_0, tg_zz_xxxxzz_1, tg_zz_xxxyy_1, tg_zz_xxxyyy_0, tg_zz_xxxyyy_1, \
                                         tg_zz_xxxyyz_0, tg_zz_xxxyyz_1, tg_zz_xxxyz_1, tg_zz_xxxyzz_0, tg_zz_xxxyzz_1, \
                                         tg_zz_xxxzz_1, tg_zz_xxxzzz_0, tg_zz_xxxzzz_1, tg_zz_xxyyy_1, tg_zz_xxyyyy_0, \
                                         tg_zz_xxyyyy_1, tg_zz_xxyyyz_0, tg_zz_xxyyyz_1, tg_zz_xxyyz_1, tg_zz_xxyyzz_0, \
                                         tg_zz_xxyyzz_1, tg_zz_xxyzz_1, tg_zz_xxyzzz_0, tg_zz_xxyzzz_1, tg_zz_xxzzz_1, \
                                         tg_zz_xxzzzz_0, tg_zz_xxzzzz_1, tg_zz_xyyyy_1, tg_zz_xyyyyy_0, tg_zz_xyyyyy_1, \
                                         tg_zz_xyyyyz_0, tg_zz_xyyyyz_1, tg_zz_xyyyz_1, tg_zz_xyyyzz_0, tg_zz_xyyyzz_1, \
                                         tg_zz_xyyzz_1, tg_zz_xyyzzz_0, tg_zz_xyyzzz_1, tg_zz_xyzzz_1, tg_zz_xyzzzz_0, \
                                         tg_zz_xyzzzz_1, tg_zz_xzzzz_1, tg_zz_xzzzzz_0, tg_zz_xzzzzz_1, tg_zz_yyyyy_1, \
                                         tg_zz_yyyyyy_0, tg_zz_yyyyyy_1, tg_zz_yyyyyz_0, tg_zz_yyyyyz_1, tg_zz_yyyyz_1, \
                                         tg_zz_yyyyzz_0, tg_zz_yyyyzz_1, tg_zz_yyyzz_1, tg_zz_yyyzzz_0, tg_zz_yyyzzz_1, \
                                         tg_zz_yyzzz_1, tg_zz_yyzzzz_0, tg_zz_yyzzzz_1, tg_zz_yzzzz_1, tg_zz_yzzzzz_0, \
                                         tg_zz_yzzzzz_1, tg_zz_zzzzz_1, tg_zz_zzzzzz_0, tg_zz_zzzzzz_1, tg_zzz_xxxxxx_0, \
                                         tg_zzz_xxxxxy_0, tg_zzz_xxxxxz_0, tg_zzz_xxxxyy_0, tg_zzz_xxxxyz_0, tg_zzz_xxxxzz_0, \
                                         tg_zzz_xxxyyy_0, tg_zzz_xxxyyz_0, tg_zzz_xxxyzz_0, tg_zzz_xxxzzz_0, tg_zzz_xxyyyy_0, \
                                         tg_zzz_xxyyyz_0, tg_zzz_xxyyzz_0, tg_zzz_xxyzzz_0, tg_zzz_xxzzzz_0, tg_zzz_xyyyyy_0, \
                                         tg_zzz_xyyyyz_0, tg_zzz_xyyyzz_0, tg_zzz_xyyzzz_0, tg_zzz_xyzzzz_0, tg_zzz_xzzzzz_0, \
                                         tg_zzz_yyyyyy_0, tg_zzz_yyyyyz_0, tg_zzz_yyyyzz_0, tg_zzz_yyyzzz_0, tg_zzz_yyzzzz_0, \
                                         tg_zzz_yzzzzz_0, tg_zzz_zzzzzz_0, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_y[j]; 

                    tg_yyy_xyzzzz_0[j] = pb_y * tg_yy_xyzzzz_0[j] + fr * tg_yy_xyzzzz_1[j] + fl1_fx * (tg_y_xyzzzz_0[j] - tg_y_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yy_xzzzz_1[j];

                    tg_yyy_xzzzzz_0[j] = pb_y * tg_yy_xzzzzz_0[j] + fr * tg_yy_xzzzzz_1[j] + fl1_fx * (tg_y_xzzzzz_0[j] - tg_y_xzzzzz_1[j] * fl1_fza);

                    tg_yyy_yyyyyy_0[j] = pb_y * tg_yy_yyyyyy_0[j] + fr * tg_yy_yyyyyy_1[j] + fl1_fx * (tg_y_yyyyyy_0[j] - tg_y_yyyyyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yy_yyyyy_1[j];

                    tg_yyy_yyyyyz_0[j] = pb_y * tg_yy_yyyyyz_0[j] + fr * tg_yy_yyyyyz_1[j] + fl1_fx * (tg_y_yyyyyz_0[j] - tg_y_yyyyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yy_yyyyz_1[j];

                    tg_yyy_yyyyzz_0[j] = pb_y * tg_yy_yyyyzz_0[j] + fr * tg_yy_yyyyzz_1[j] + fl1_fx * (tg_y_yyyyzz_0[j] - tg_y_yyyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yy_yyyzz_1[j];

                    tg_yyy_yyyzzz_0[j] = pb_y * tg_yy_yyyzzz_0[j] + fr * tg_yy_yyyzzz_1[j] + fl1_fx * (tg_y_yyyzzz_0[j] - tg_y_yyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yy_yyzzz_1[j];

                    tg_yyy_yyzzzz_0[j] = pb_y * tg_yy_yyzzzz_0[j] + fr * tg_yy_yyzzzz_1[j] + fl1_fx * (tg_y_yyzzzz_0[j] - tg_y_yyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yy_yzzzz_1[j];

                    tg_yyy_yzzzzz_0[j] = pb_y * tg_yy_yzzzzz_0[j] + fr * tg_yy_yzzzzz_1[j] + fl1_fx * (tg_y_yzzzzz_0[j] - tg_y_yzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yy_zzzzz_1[j];

                    tg_yyy_zzzzzz_0[j] = pb_y * tg_yy_zzzzzz_0[j] + fr * tg_yy_zzzzzz_1[j] + fl1_fx * (tg_y_zzzzzz_0[j] - tg_y_zzzzzz_1[j] * fl1_fza);

                    tg_yyz_xxxxxx_0[j] = pb_y * tg_yz_xxxxxx_0[j] + fr * tg_yz_xxxxxx_1[j] + 0.5 * fl1_fx * (tg_z_xxxxxx_0[j] - tg_z_xxxxxx_1[j] * fl1_fza);

                    tg_yyz_xxxxxy_0[j] = pb_y * tg_yz_xxxxxy_0[j] + fr * tg_yz_xxxxxy_1[j] + 0.5 * fl1_fx * (tg_z_xxxxxy_0[j] - tg_z_xxxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yz_xxxxx_1[j];

                    tg_yyz_xxxxxz_0[j] = pb_y * tg_yz_xxxxxz_0[j] + fr * tg_yz_xxxxxz_1[j] + 0.5 * fl1_fx * (tg_z_xxxxxz_0[j] - tg_z_xxxxxz_1[j] * fl1_fza);

                    tg_yyz_xxxxyy_0[j] = pb_y * tg_yz_xxxxyy_0[j] + fr * tg_yz_xxxxyy_1[j] + 0.5 * fl1_fx * (tg_z_xxxxyy_0[j] - tg_z_xxxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yz_xxxxy_1[j];

                    tg_yyz_xxxxyz_0[j] = pb_y * tg_yz_xxxxyz_0[j] + fr * tg_yz_xxxxyz_1[j] + 0.5 * fl1_fx * (tg_z_xxxxyz_0[j] - tg_z_xxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yz_xxxxz_1[j];

                    tg_yyz_xxxxzz_0[j] = pb_y * tg_yz_xxxxzz_0[j] + fr * tg_yz_xxxxzz_1[j] + 0.5 * fl1_fx * (tg_z_xxxxzz_0[j] - tg_z_xxxxzz_1[j] * fl1_fza);

                    tg_yyz_xxxyyy_0[j] = pb_y * tg_yz_xxxyyy_0[j] + fr * tg_yz_xxxyyy_1[j] + 0.5 * fl1_fx * (tg_z_xxxyyy_0[j] - tg_z_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yz_xxxyy_1[j];

                    tg_yyz_xxxyyz_0[j] = pb_y * tg_yz_xxxyyz_0[j] + fr * tg_yz_xxxyyz_1[j] + 0.5 * fl1_fx * (tg_z_xxxyyz_0[j] - tg_z_xxxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yz_xxxyz_1[j];

                    tg_yyz_xxxyzz_0[j] = pb_y * tg_yz_xxxyzz_0[j] + fr * tg_yz_xxxyzz_1[j] + 0.5 * fl1_fx * (tg_z_xxxyzz_0[j] - tg_z_xxxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yz_xxxzz_1[j];

                    tg_yyz_xxxzzz_0[j] = pb_y * tg_yz_xxxzzz_0[j] + fr * tg_yz_xxxzzz_1[j] + 0.5 * fl1_fx * (tg_z_xxxzzz_0[j] - tg_z_xxxzzz_1[j] * fl1_fza);

                    tg_yyz_xxyyyy_0[j] = pb_y * tg_yz_xxyyyy_0[j] + fr * tg_yz_xxyyyy_1[j] + 0.5 * fl1_fx * (tg_z_xxyyyy_0[j] - tg_z_xxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yz_xxyyy_1[j];

                    tg_yyz_xxyyyz_0[j] = pb_y * tg_yz_xxyyyz_0[j] + fr * tg_yz_xxyyyz_1[j] + 0.5 * fl1_fx * (tg_z_xxyyyz_0[j] - tg_z_xxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yz_xxyyz_1[j];

                    tg_yyz_xxyyzz_0[j] = pb_y * tg_yz_xxyyzz_0[j] + fr * tg_yz_xxyyzz_1[j] + 0.5 * fl1_fx * (tg_z_xxyyzz_0[j] - tg_z_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yz_xxyzz_1[j];

                    tg_yyz_xxyzzz_0[j] = pb_y * tg_yz_xxyzzz_0[j] + fr * tg_yz_xxyzzz_1[j] + 0.5 * fl1_fx * (tg_z_xxyzzz_0[j] - tg_z_xxyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yz_xxzzz_1[j];

                    tg_yyz_xxzzzz_0[j] = pb_y * tg_yz_xxzzzz_0[j] + fr * tg_yz_xxzzzz_1[j] + 0.5 * fl1_fx * (tg_z_xxzzzz_0[j] - tg_z_xxzzzz_1[j] * fl1_fza);

                    tg_yyz_xyyyyy_0[j] = pb_y * tg_yz_xyyyyy_0[j] + fr * tg_yz_xyyyyy_1[j] + 0.5 * fl1_fx * (tg_z_xyyyyy_0[j] - tg_z_xyyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yz_xyyyy_1[j];

                    tg_yyz_xyyyyz_0[j] = pb_y * tg_yz_xyyyyz_0[j] + fr * tg_yz_xyyyyz_1[j] + 0.5 * fl1_fx * (tg_z_xyyyyz_0[j] - tg_z_xyyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yz_xyyyz_1[j];

                    tg_yyz_xyyyzz_0[j] = pb_y * tg_yz_xyyyzz_0[j] + fr * tg_yz_xyyyzz_1[j] + 0.5 * fl1_fx * (tg_z_xyyyzz_0[j] - tg_z_xyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yz_xyyzz_1[j];

                    tg_yyz_xyyzzz_0[j] = pb_y * tg_yz_xyyzzz_0[j] + fr * tg_yz_xyyzzz_1[j] + 0.5 * fl1_fx * (tg_z_xyyzzz_0[j] - tg_z_xyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yz_xyzzz_1[j];

                    tg_yyz_xyzzzz_0[j] = pb_y * tg_yz_xyzzzz_0[j] + fr * tg_yz_xyzzzz_1[j] + 0.5 * fl1_fx * (tg_z_xyzzzz_0[j] - tg_z_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yz_xzzzz_1[j];

                    tg_yyz_xzzzzz_0[j] = pb_y * tg_yz_xzzzzz_0[j] + fr * tg_yz_xzzzzz_1[j] + 0.5 * fl1_fx * (tg_z_xzzzzz_0[j] - tg_z_xzzzzz_1[j] * fl1_fza);

                    tg_yyz_yyyyyy_0[j] = pb_y * tg_yz_yyyyyy_0[j] + fr * tg_yz_yyyyyy_1[j] + 0.5 * fl1_fx * (tg_z_yyyyyy_0[j] - tg_z_yyyyyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yz_yyyyy_1[j];

                    tg_yyz_yyyyyz_0[j] = pb_y * tg_yz_yyyyyz_0[j] + fr * tg_yz_yyyyyz_1[j] + 0.5 * fl1_fx * (tg_z_yyyyyz_0[j] - tg_z_yyyyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yz_yyyyz_1[j];

                    tg_yyz_yyyyzz_0[j] = pb_y * tg_yz_yyyyzz_0[j] + fr * tg_yz_yyyyzz_1[j] + 0.5 * fl1_fx * (tg_z_yyyyzz_0[j] - tg_z_yyyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yz_yyyzz_1[j];

                    tg_yyz_yyyzzz_0[j] = pb_y * tg_yz_yyyzzz_0[j] + fr * tg_yz_yyyzzz_1[j] + 0.5 * fl1_fx * (tg_z_yyyzzz_0[j] - tg_z_yyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yz_yyzzz_1[j];

                    tg_yyz_yyzzzz_0[j] = pb_y * tg_yz_yyzzzz_0[j] + fr * tg_yz_yyzzzz_1[j] + 0.5 * fl1_fx * (tg_z_yyzzzz_0[j] - tg_z_yyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yz_yzzzz_1[j];

                    tg_yyz_yzzzzz_0[j] = pb_y * tg_yz_yzzzzz_0[j] + fr * tg_yz_yzzzzz_1[j] + 0.5 * fl1_fx * (tg_z_yzzzzz_0[j] - tg_z_yzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yz_zzzzz_1[j];

                    tg_yyz_zzzzzz_0[j] = pb_y * tg_yz_zzzzzz_0[j] + fr * tg_yz_zzzzzz_1[j] + 0.5 * fl1_fx * (tg_z_zzzzzz_0[j] - tg_z_zzzzzz_1[j] * fl1_fza);

                    tg_yzz_xxxxxx_0[j] = pb_y * tg_zz_xxxxxx_0[j] + fr * tg_zz_xxxxxx_1[j];

                    tg_yzz_xxxxxy_0[j] = pb_y * tg_zz_xxxxxy_0[j] + fr * tg_zz_xxxxxy_1[j] + 0.5 * fl1_fxn * tg_zz_xxxxx_1[j];

                    tg_yzz_xxxxxz_0[j] = pb_y * tg_zz_xxxxxz_0[j] + fr * tg_zz_xxxxxz_1[j];

                    tg_yzz_xxxxyy_0[j] = pb_y * tg_zz_xxxxyy_0[j] + fr * tg_zz_xxxxyy_1[j] + fl1_fxn * tg_zz_xxxxy_1[j];

                    tg_yzz_xxxxyz_0[j] = pb_y * tg_zz_xxxxyz_0[j] + fr * tg_zz_xxxxyz_1[j] + 0.5 * fl1_fxn * tg_zz_xxxxz_1[j];

                    tg_yzz_xxxxzz_0[j] = pb_y * tg_zz_xxxxzz_0[j] + fr * tg_zz_xxxxzz_1[j];

                    tg_yzz_xxxyyy_0[j] = pb_y * tg_zz_xxxyyy_0[j] + fr * tg_zz_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_zz_xxxyy_1[j];

                    tg_yzz_xxxyyz_0[j] = pb_y * tg_zz_xxxyyz_0[j] + fr * tg_zz_xxxyyz_1[j] + fl1_fxn * tg_zz_xxxyz_1[j];

                    tg_yzz_xxxyzz_0[j] = pb_y * tg_zz_xxxyzz_0[j] + fr * tg_zz_xxxyzz_1[j] + 0.5 * fl1_fxn * tg_zz_xxxzz_1[j];

                    tg_yzz_xxxzzz_0[j] = pb_y * tg_zz_xxxzzz_0[j] + fr * tg_zz_xxxzzz_1[j];

                    tg_yzz_xxyyyy_0[j] = pb_y * tg_zz_xxyyyy_0[j] + fr * tg_zz_xxyyyy_1[j] + 2.0 * fl1_fxn * tg_zz_xxyyy_1[j];

                    tg_yzz_xxyyyz_0[j] = pb_y * tg_zz_xxyyyz_0[j] + fr * tg_zz_xxyyyz_1[j] + 1.5 * fl1_fxn * tg_zz_xxyyz_1[j];

                    tg_yzz_xxyyzz_0[j] = pb_y * tg_zz_xxyyzz_0[j] + fr * tg_zz_xxyyzz_1[j] + fl1_fxn * tg_zz_xxyzz_1[j];

                    tg_yzz_xxyzzz_0[j] = pb_y * tg_zz_xxyzzz_0[j] + fr * tg_zz_xxyzzz_1[j] + 0.5 * fl1_fxn * tg_zz_xxzzz_1[j];

                    tg_yzz_xxzzzz_0[j] = pb_y * tg_zz_xxzzzz_0[j] + fr * tg_zz_xxzzzz_1[j];

                    tg_yzz_xyyyyy_0[j] = pb_y * tg_zz_xyyyyy_0[j] + fr * tg_zz_xyyyyy_1[j] + 2.5 * fl1_fxn * tg_zz_xyyyy_1[j];

                    tg_yzz_xyyyyz_0[j] = pb_y * tg_zz_xyyyyz_0[j] + fr * tg_zz_xyyyyz_1[j] + 2.0 * fl1_fxn * tg_zz_xyyyz_1[j];

                    tg_yzz_xyyyzz_0[j] = pb_y * tg_zz_xyyyzz_0[j] + fr * tg_zz_xyyyzz_1[j] + 1.5 * fl1_fxn * tg_zz_xyyzz_1[j];

                    tg_yzz_xyyzzz_0[j] = pb_y * tg_zz_xyyzzz_0[j] + fr * tg_zz_xyyzzz_1[j] + fl1_fxn * tg_zz_xyzzz_1[j];

                    tg_yzz_xyzzzz_0[j] = pb_y * tg_zz_xyzzzz_0[j] + fr * tg_zz_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_zz_xzzzz_1[j];

                    tg_yzz_xzzzzz_0[j] = pb_y * tg_zz_xzzzzz_0[j] + fr * tg_zz_xzzzzz_1[j];

                    tg_yzz_yyyyyy_0[j] = pb_y * tg_zz_yyyyyy_0[j] + fr * tg_zz_yyyyyy_1[j] + 3.0 * fl1_fxn * tg_zz_yyyyy_1[j];

                    tg_yzz_yyyyyz_0[j] = pb_y * tg_zz_yyyyyz_0[j] + fr * tg_zz_yyyyyz_1[j] + 2.5 * fl1_fxn * tg_zz_yyyyz_1[j];

                    tg_yzz_yyyyzz_0[j] = pb_y * tg_zz_yyyyzz_0[j] + fr * tg_zz_yyyyzz_1[j] + 2.0 * fl1_fxn * tg_zz_yyyzz_1[j];

                    tg_yzz_yyyzzz_0[j] = pb_y * tg_zz_yyyzzz_0[j] + fr * tg_zz_yyyzzz_1[j] + 1.5 * fl1_fxn * tg_zz_yyzzz_1[j];

                    tg_yzz_yyzzzz_0[j] = pb_y * tg_zz_yyzzzz_0[j] + fr * tg_zz_yyzzzz_1[j] + fl1_fxn * tg_zz_yzzzz_1[j];

                    tg_yzz_yzzzzz_0[j] = pb_y * tg_zz_yzzzzz_0[j] + fr * tg_zz_yzzzzz_1[j] + 0.5 * fl1_fxn * tg_zz_zzzzz_1[j];

                    tg_yzz_zzzzzz_0[j] = pb_y * tg_zz_zzzzzz_0[j] + fr * tg_zz_zzzzzz_1[j];

                    fr = wp_z[j]; 

                    tg_zzz_xxxxxx_0[j] = pb_z * tg_zz_xxxxxx_0[j] + fr * tg_zz_xxxxxx_1[j] + fl1_fx * (tg_z_xxxxxx_0[j] - tg_z_xxxxxx_1[j] * fl1_fza);

                    tg_zzz_xxxxxy_0[j] = pb_z * tg_zz_xxxxxy_0[j] + fr * tg_zz_xxxxxy_1[j] + fl1_fx * (tg_z_xxxxxy_0[j] - tg_z_xxxxxy_1[j] * fl1_fza);

                    tg_zzz_xxxxxz_0[j] = pb_z * tg_zz_xxxxxz_0[j] + fr * tg_zz_xxxxxz_1[j] + fl1_fx * (tg_z_xxxxxz_0[j] - tg_z_xxxxxz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zz_xxxxx_1[j];

                    tg_zzz_xxxxyy_0[j] = pb_z * tg_zz_xxxxyy_0[j] + fr * tg_zz_xxxxyy_1[j] + fl1_fx * (tg_z_xxxxyy_0[j] - tg_z_xxxxyy_1[j] * fl1_fza);

                    tg_zzz_xxxxyz_0[j] = pb_z * tg_zz_xxxxyz_0[j] + fr * tg_zz_xxxxyz_1[j] + fl1_fx * (tg_z_xxxxyz_0[j] - tg_z_xxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zz_xxxxy_1[j];

                    tg_zzz_xxxxzz_0[j] = pb_z * tg_zz_xxxxzz_0[j] + fr * tg_zz_xxxxzz_1[j] + fl1_fx * (tg_z_xxxxzz_0[j] - tg_z_xxxxzz_1[j] * fl1_fza) + fl1_fxn * tg_zz_xxxxz_1[j];

                    tg_zzz_xxxyyy_0[j] = pb_z * tg_zz_xxxyyy_0[j] + fr * tg_zz_xxxyyy_1[j] + fl1_fx * (tg_z_xxxyyy_0[j] - tg_z_xxxyyy_1[j] * fl1_fza);

                    tg_zzz_xxxyyz_0[j] = pb_z * tg_zz_xxxyyz_0[j] + fr * tg_zz_xxxyyz_1[j] + fl1_fx * (tg_z_xxxyyz_0[j] - tg_z_xxxyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zz_xxxyy_1[j];

                    tg_zzz_xxxyzz_0[j] = pb_z * tg_zz_xxxyzz_0[j] + fr * tg_zz_xxxyzz_1[j] + fl1_fx * (tg_z_xxxyzz_0[j] - tg_z_xxxyzz_1[j] * fl1_fza) + fl1_fxn * tg_zz_xxxyz_1[j];

                    tg_zzz_xxxzzz_0[j] = pb_z * tg_zz_xxxzzz_0[j] + fr * tg_zz_xxxzzz_1[j] + fl1_fx * (tg_z_xxxzzz_0[j] - tg_z_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zz_xxxzz_1[j];

                    tg_zzz_xxyyyy_0[j] = pb_z * tg_zz_xxyyyy_0[j] + fr * tg_zz_xxyyyy_1[j] + fl1_fx * (tg_z_xxyyyy_0[j] - tg_z_xxyyyy_1[j] * fl1_fza);

                    tg_zzz_xxyyyz_0[j] = pb_z * tg_zz_xxyyyz_0[j] + fr * tg_zz_xxyyyz_1[j] + fl1_fx * (tg_z_xxyyyz_0[j] - tg_z_xxyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zz_xxyyy_1[j];

                    tg_zzz_xxyyzz_0[j] = pb_z * tg_zz_xxyyzz_0[j] + fr * tg_zz_xxyyzz_1[j] + fl1_fx * (tg_z_xxyyzz_0[j] - tg_z_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zz_xxyyz_1[j];

                    tg_zzz_xxyzzz_0[j] = pb_z * tg_zz_xxyzzz_0[j] + fr * tg_zz_xxyzzz_1[j] + fl1_fx * (tg_z_xxyzzz_0[j] - tg_z_xxyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zz_xxyzz_1[j];

                    tg_zzz_xxzzzz_0[j] = pb_z * tg_zz_xxzzzz_0[j] + fr * tg_zz_xxzzzz_1[j] + fl1_fx * (tg_z_xxzzzz_0[j] - tg_z_xxzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zz_xxzzz_1[j];

                    tg_zzz_xyyyyy_0[j] = pb_z * tg_zz_xyyyyy_0[j] + fr * tg_zz_xyyyyy_1[j] + fl1_fx * (tg_z_xyyyyy_0[j] - tg_z_xyyyyy_1[j] * fl1_fza);

                    tg_zzz_xyyyyz_0[j] = pb_z * tg_zz_xyyyyz_0[j] + fr * tg_zz_xyyyyz_1[j] + fl1_fx * (tg_z_xyyyyz_0[j] - tg_z_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zz_xyyyy_1[j];

                    tg_zzz_xyyyzz_0[j] = pb_z * tg_zz_xyyyzz_0[j] + fr * tg_zz_xyyyzz_1[j] + fl1_fx * (tg_z_xyyyzz_0[j] - tg_z_xyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zz_xyyyz_1[j];

                    tg_zzz_xyyzzz_0[j] = pb_z * tg_zz_xyyzzz_0[j] + fr * tg_zz_xyyzzz_1[j] + fl1_fx * (tg_z_xyyzzz_0[j] - tg_z_xyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zz_xyyzz_1[j];

                    tg_zzz_xyzzzz_0[j] = pb_z * tg_zz_xyzzzz_0[j] + fr * tg_zz_xyzzzz_1[j] + fl1_fx * (tg_z_xyzzzz_0[j] - tg_z_xyzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zz_xyzzz_1[j];

                    tg_zzz_xzzzzz_0[j] = pb_z * tg_zz_xzzzzz_0[j] + fr * tg_zz_xzzzzz_1[j] + fl1_fx * (tg_z_xzzzzz_0[j] - tg_z_xzzzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_zz_xzzzz_1[j];

                    tg_zzz_yyyyyy_0[j] = pb_z * tg_zz_yyyyyy_0[j] + fr * tg_zz_yyyyyy_1[j] + fl1_fx * (tg_z_yyyyyy_0[j] - tg_z_yyyyyy_1[j] * fl1_fza);

                    tg_zzz_yyyyyz_0[j] = pb_z * tg_zz_yyyyyz_0[j] + fr * tg_zz_yyyyyz_1[j] + fl1_fx * (tg_z_yyyyyz_0[j] - tg_z_yyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zz_yyyyy_1[j];

                    tg_zzz_yyyyzz_0[j] = pb_z * tg_zz_yyyyzz_0[j] + fr * tg_zz_yyyyzz_1[j] + fl1_fx * (tg_z_yyyyzz_0[j] - tg_z_yyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zz_yyyyz_1[j];

                    tg_zzz_yyyzzz_0[j] = pb_z * tg_zz_yyyzzz_0[j] + fr * tg_zz_yyyzzz_1[j] + fl1_fx * (tg_z_yyyzzz_0[j] - tg_z_yyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zz_yyyzz_1[j];

                    tg_zzz_yyzzzz_0[j] = pb_z * tg_zz_yyzzzz_0[j] + fr * tg_zz_yyzzzz_1[j] + fl1_fx * (tg_z_yyzzzz_0[j] - tg_z_yyzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zz_yyzzz_1[j];

                    tg_zzz_yzzzzz_0[j] = pb_z * tg_zz_yzzzzz_0[j] + fr * tg_zz_yzzzzz_1[j] + fl1_fx * (tg_z_yzzzzz_0[j] - tg_z_yzzzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_zz_yzzzz_1[j];

                    tg_zzz_zzzzzz_0[j] = pb_z * tg_zz_zzzzzz_0[j] + fr * tg_zz_zzzzzz_1[j] + fl1_fx * (tg_z_zzzzzz_0[j] - tg_z_zzzzzz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_zz_zzzzz_1[j];
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

