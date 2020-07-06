//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ElectronRepulsionRecFuncForII.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSISI(      CMemBlock2D<double>* primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSISI_0_98(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSISI_98_196(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 

        erirecfunc::compElectronRepulsionForSISI_196_294(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSISI_294_392(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSISI_392_490(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSISI_490_588(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSISI_588_686(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSISI_686_784(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 
    }

    void
    compElectronRepulsionForSISI_0_98(      CMemBlock2D<double>* primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& wpDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetPrimPairs,
                                      const int32_t              iContrPair)
    {
        // Batch of Integrals (0,98)

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
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_6_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_6_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_5_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
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

                auto tg_xxxxx_xxxxxx_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx); 

                auto tg_xxxxx_xxxxxy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 1); 

                auto tg_xxxxx_xxxxxz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 2); 

                auto tg_xxxxx_xxxxyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 3); 

                auto tg_xxxxx_xxxxyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 4); 

                auto tg_xxxxx_xxxxzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 5); 

                auto tg_xxxxx_xxxyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 6); 

                auto tg_xxxxx_xxxyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 7); 

                auto tg_xxxxx_xxxyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 8); 

                auto tg_xxxxx_xxxzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 9); 

                auto tg_xxxxx_xxyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 10); 

                auto tg_xxxxx_xxyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 11); 

                auto tg_xxxxx_xxyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 12); 

                auto tg_xxxxx_xxyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 13); 

                auto tg_xxxxx_xxzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 14); 

                auto tg_xxxxx_xyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 15); 

                auto tg_xxxxx_xyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 16); 

                auto tg_xxxxx_xyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 17); 

                auto tg_xxxxx_xyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 18); 

                auto tg_xxxxx_xyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 19); 

                auto tg_xxxxx_xzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 20); 

                auto tg_xxxxx_yyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 21); 

                auto tg_xxxxx_yyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 22); 

                auto tg_xxxxx_yyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 23); 

                auto tg_xxxxx_yyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 24); 

                auto tg_xxxxx_yyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 25); 

                auto tg_xxxxx_yzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 26); 

                auto tg_xxxxx_zzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 27); 

                auto tg_xxxxy_xxxxxx_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 28); 

                auto tg_xxxxy_xxxxxy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 29); 

                auto tg_xxxxy_xxxxxz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 30); 

                auto tg_xxxxy_xxxxyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 31); 

                auto tg_xxxxy_xxxxyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 32); 

                auto tg_xxxxy_xxxxzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 33); 

                auto tg_xxxxy_xxxyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 34); 

                auto tg_xxxxy_xxxyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 35); 

                auto tg_xxxxy_xxxyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 36); 

                auto tg_xxxxy_xxxzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 37); 

                auto tg_xxxxy_xxyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 38); 

                auto tg_xxxxy_xxyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 39); 

                auto tg_xxxxy_xxyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 40); 

                auto tg_xxxxy_xxyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 41); 

                auto tg_xxxxy_xxzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 42); 

                auto tg_xxxxy_xyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 43); 

                auto tg_xxxxy_xyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 44); 

                auto tg_xxxxy_xyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 45); 

                auto tg_xxxxy_xyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 46); 

                auto tg_xxxxy_xyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 47); 

                auto tg_xxxxy_xzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 48); 

                auto tg_xxxxy_yyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 49); 

                auto tg_xxxxy_yyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 50); 

                auto tg_xxxxy_yyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 51); 

                auto tg_xxxxy_yyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 52); 

                auto tg_xxxxy_yyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 53); 

                auto tg_xxxxy_yzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 54); 

                auto tg_xxxxy_zzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 55); 

                auto tg_xxxxz_xxxxxx_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 56); 

                auto tg_xxxxz_xxxxxy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 57); 

                auto tg_xxxxz_xxxxxz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 58); 

                auto tg_xxxxz_xxxxyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 59); 

                auto tg_xxxxz_xxxxyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 60); 

                auto tg_xxxxz_xxxxzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 61); 

                auto tg_xxxxz_xxxyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 62); 

                auto tg_xxxxz_xxxyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 63); 

                auto tg_xxxxz_xxxyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 64); 

                auto tg_xxxxz_xxxzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 65); 

                auto tg_xxxxz_xxyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 66); 

                auto tg_xxxxz_xxyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 67); 

                auto tg_xxxxz_xxyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 68); 

                auto tg_xxxxz_xxyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 69); 

                auto tg_xxxxz_xxzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 70); 

                auto tg_xxxxz_xyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 71); 

                auto tg_xxxxz_xyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 72); 

                auto tg_xxxxz_xyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 73); 

                auto tg_xxxxz_xyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 74); 

                auto tg_xxxxz_xyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 75); 

                auto tg_xxxxz_xzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 76); 

                auto tg_xxxxz_yyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 77); 

                auto tg_xxxxz_yyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 78); 

                auto tg_xxxxz_yyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 79); 

                auto tg_xxxxz_yyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 80); 

                auto tg_xxxxz_yyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 81); 

                auto tg_xxxxz_yzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 82); 

                auto tg_xxxxz_zzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 83); 

                auto tg_xxxyy_xxxxxx_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 84); 

                auto tg_xxxyy_xxxxxy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 85); 

                auto tg_xxxyy_xxxxxz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 86); 

                auto tg_xxxyy_xxxxyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 87); 

                auto tg_xxxyy_xxxxyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 88); 

                auto tg_xxxyy_xxxxzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 89); 

                auto tg_xxxyy_xxxyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 90); 

                auto tg_xxxyy_xxxyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 91); 

                auto tg_xxxyy_xxxyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 92); 

                auto tg_xxxyy_xxxzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 93); 

                auto tg_xxxyy_xxyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 94); 

                auto tg_xxxyy_xxyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 95); 

                auto tg_xxxyy_xxyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 96); 

                auto tg_xxxyy_xxyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 97); 

                auto tg_xxxxx_xxxxxx_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx); 

                auto tg_xxxxx_xxxxxy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 1); 

                auto tg_xxxxx_xxxxxz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 2); 

                auto tg_xxxxx_xxxxyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 3); 

                auto tg_xxxxx_xxxxyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 4); 

                auto tg_xxxxx_xxxxzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 5); 

                auto tg_xxxxx_xxxyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 6); 

                auto tg_xxxxx_xxxyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 7); 

                auto tg_xxxxx_xxxyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 8); 

                auto tg_xxxxx_xxxzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 9); 

                auto tg_xxxxx_xxyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 10); 

                auto tg_xxxxx_xxyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 11); 

                auto tg_xxxxx_xxyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 12); 

                auto tg_xxxxx_xxyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 13); 

                auto tg_xxxxx_xxzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 14); 

                auto tg_xxxxx_xyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 15); 

                auto tg_xxxxx_xyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 16); 

                auto tg_xxxxx_xyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 17); 

                auto tg_xxxxx_xyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 18); 

                auto tg_xxxxx_xyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 19); 

                auto tg_xxxxx_xzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 20); 

                auto tg_xxxxx_yyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 21); 

                auto tg_xxxxx_yyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 22); 

                auto tg_xxxxx_yyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 23); 

                auto tg_xxxxx_yyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 24); 

                auto tg_xxxxx_yyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 25); 

                auto tg_xxxxx_yzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 26); 

                auto tg_xxxxx_zzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 27); 

                auto tg_xxxxy_xxxxxx_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 28); 

                auto tg_xxxxy_xxxxxy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 29); 

                auto tg_xxxxy_xxxxxz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 30); 

                auto tg_xxxxy_xxxxyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 31); 

                auto tg_xxxxy_xxxxyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 32); 

                auto tg_xxxxy_xxxxzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 33); 

                auto tg_xxxxy_xxxyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 34); 

                auto tg_xxxxy_xxxyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 35); 

                auto tg_xxxxy_xxxyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 36); 

                auto tg_xxxxy_xxxzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 37); 

                auto tg_xxxxy_xxyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 38); 

                auto tg_xxxxy_xxyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 39); 

                auto tg_xxxxy_xxyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 40); 

                auto tg_xxxxy_xxyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 41); 

                auto tg_xxxxy_xxzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 42); 

                auto tg_xxxxy_xyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 43); 

                auto tg_xxxxy_xyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 44); 

                auto tg_xxxxy_xyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 45); 

                auto tg_xxxxy_xyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 46); 

                auto tg_xxxxy_xyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 47); 

                auto tg_xxxxy_xzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 48); 

                auto tg_xxxxy_yyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 49); 

                auto tg_xxxxy_yyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 50); 

                auto tg_xxxxy_yyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 51); 

                auto tg_xxxxy_yyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 52); 

                auto tg_xxxxy_yyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 53); 

                auto tg_xxxxy_yzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 54); 

                auto tg_xxxxy_zzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 55); 

                auto tg_xxxxz_xxxxxx_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 56); 

                auto tg_xxxxz_xxxxxy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 57); 

                auto tg_xxxxz_xxxxxz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 58); 

                auto tg_xxxxz_xxxxyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 59); 

                auto tg_xxxxz_xxxxyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 60); 

                auto tg_xxxxz_xxxxzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 61); 

                auto tg_xxxxz_xxxyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 62); 

                auto tg_xxxxz_xxxyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 63); 

                auto tg_xxxxz_xxxyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 64); 

                auto tg_xxxxz_xxxzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 65); 

                auto tg_xxxxz_xxyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 66); 

                auto tg_xxxxz_xxyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 67); 

                auto tg_xxxxz_xxyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 68); 

                auto tg_xxxxz_xxyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 69); 

                auto tg_xxxxz_xxzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 70); 

                auto tg_xxxxz_xyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 71); 

                auto tg_xxxxz_xyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 72); 

                auto tg_xxxxz_xyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 73); 

                auto tg_xxxxz_xyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 74); 

                auto tg_xxxxz_xyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 75); 

                auto tg_xxxxz_xzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 76); 

                auto tg_xxxxz_yyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 77); 

                auto tg_xxxxz_yyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 78); 

                auto tg_xxxxz_yyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 79); 

                auto tg_xxxxz_yyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 80); 

                auto tg_xxxxz_yyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 81); 

                auto tg_xxxxz_yzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 82); 

                auto tg_xxxxz_zzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 83); 

                auto tg_xxxyy_xxxxxx_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 84); 

                auto tg_xxxyy_xxxxxy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 85); 

                auto tg_xxxyy_xxxxxz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 86); 

                auto tg_xxxyy_xxxxyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 87); 

                auto tg_xxxyy_xxxxyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 88); 

                auto tg_xxxyy_xxxxzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 89); 

                auto tg_xxxyy_xxxyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 90); 

                auto tg_xxxyy_xxxyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 91); 

                auto tg_xxxyy_xxxyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 92); 

                auto tg_xxxyy_xxxzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 93); 

                auto tg_xxxyy_xxyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 94); 

                auto tg_xxxyy_xxyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 95); 

                auto tg_xxxyy_xxyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 96); 

                auto tg_xxxyy_xxyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 97); 

                auto tg_xxxx_xxxxxx_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx); 

                auto tg_xxxx_xxxxxy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 1); 

                auto tg_xxxx_xxxxxz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 2); 

                auto tg_xxxx_xxxxyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 3); 

                auto tg_xxxx_xxxxyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 4); 

                auto tg_xxxx_xxxxzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 5); 

                auto tg_xxxx_xxxyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 6); 

                auto tg_xxxx_xxxyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 7); 

                auto tg_xxxx_xxxyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 8); 

                auto tg_xxxx_xxxzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 9); 

                auto tg_xxxx_xxyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 10); 

                auto tg_xxxx_xxyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 11); 

                auto tg_xxxx_xxyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 12); 

                auto tg_xxxx_xxyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 13); 

                auto tg_xxxx_xxzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 14); 

                auto tg_xxxx_xyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 15); 

                auto tg_xxxx_xyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 16); 

                auto tg_xxxx_xyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 17); 

                auto tg_xxxx_xyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 18); 

                auto tg_xxxx_xyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 19); 

                auto tg_xxxx_xzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 20); 

                auto tg_xxxx_yyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 21); 

                auto tg_xxxx_yyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 22); 

                auto tg_xxxx_yyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 23); 

                auto tg_xxxx_yyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 24); 

                auto tg_xxxx_yyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 25); 

                auto tg_xxxx_yzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 26); 

                auto tg_xxxx_zzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 27); 

                auto tg_xxxy_xxxxxx_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 28); 

                auto tg_xxxy_xxxxxy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 29); 

                auto tg_xxxy_xxxxxz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 30); 

                auto tg_xxxy_xxxxyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 31); 

                auto tg_xxxy_xxxxyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 32); 

                auto tg_xxxy_xxxxzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 33); 

                auto tg_xxxy_xxxyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 34); 

                auto tg_xxxy_xxxyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 35); 

                auto tg_xxxy_xxxyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 36); 

                auto tg_xxxy_xxxzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 37); 

                auto tg_xxxy_xxyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 38); 

                auto tg_xxxy_xxyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 39); 

                auto tg_xxxy_xxyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 40); 

                auto tg_xxxy_xxyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 41); 

                auto tg_xxxy_xxzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 42); 

                auto tg_xxxy_xyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 43); 

                auto tg_xxxy_xyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 44); 

                auto tg_xxxy_xyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 45); 

                auto tg_xxxy_xyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 46); 

                auto tg_xxxy_xyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 47); 

                auto tg_xxxy_xzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 48); 

                auto tg_xxxy_yyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 49); 

                auto tg_xxxy_yyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 50); 

                auto tg_xxxy_yyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 51); 

                auto tg_xxxy_yyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 52); 

                auto tg_xxxy_yyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 53); 

                auto tg_xxxy_yzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 54); 

                auto tg_xxxy_zzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 55); 

                auto tg_xxxz_xxxxxx_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 56); 

                auto tg_xxxz_xxxxxy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 57); 

                auto tg_xxxz_xxxxxz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 58); 

                auto tg_xxxz_xxxxyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 59); 

                auto tg_xxxz_xxxxyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 60); 

                auto tg_xxxz_xxxxzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 61); 

                auto tg_xxxz_xxxyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 62); 

                auto tg_xxxz_xxxyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 63); 

                auto tg_xxxz_xxxyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 64); 

                auto tg_xxxz_xxxzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 65); 

                auto tg_xxxz_xxyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 66); 

                auto tg_xxxz_xxyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 67); 

                auto tg_xxxz_xxyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 68); 

                auto tg_xxxz_xxyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 69); 

                auto tg_xxxz_xxzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 70); 

                auto tg_xxxz_xyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 71); 

                auto tg_xxxz_xyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 72); 

                auto tg_xxxz_xyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 73); 

                auto tg_xxxz_xyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 74); 

                auto tg_xxxz_xyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 75); 

                auto tg_xxxz_xzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 76); 

                auto tg_xxxz_yyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 77); 

                auto tg_xxxz_yyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 78); 

                auto tg_xxxz_yyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 79); 

                auto tg_xxxz_yyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 80); 

                auto tg_xxxz_yyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 81); 

                auto tg_xxxz_yzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 82); 

                auto tg_xxxz_zzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 83); 

                auto tg_xxyy_xxxxxx_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 84); 

                auto tg_xxyy_xxxxxy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 85); 

                auto tg_xxyy_xxxxxz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 86); 

                auto tg_xxyy_xxxxyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 87); 

                auto tg_xxyy_xxxxyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 88); 

                auto tg_xxyy_xxxxzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 89); 

                auto tg_xxyy_xxxyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 90); 

                auto tg_xxyy_xxxyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 91); 

                auto tg_xxyy_xxxyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 92); 

                auto tg_xxyy_xxxzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 93); 

                auto tg_xxyy_xxyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 94); 

                auto tg_xxyy_xxyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 95); 

                auto tg_xxyy_xxyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 96); 

                auto tg_xxyy_xxyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 97); 

                auto tg_xxxx_xxxxxx_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx); 

                auto tg_xxxx_xxxxxy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 1); 

                auto tg_xxxx_xxxxxz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 2); 

                auto tg_xxxx_xxxxyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 3); 

                auto tg_xxxx_xxxxyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 4); 

                auto tg_xxxx_xxxxzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 5); 

                auto tg_xxxx_xxxyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 6); 

                auto tg_xxxx_xxxyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 7); 

                auto tg_xxxx_xxxyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 8); 

                auto tg_xxxx_xxxzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 9); 

                auto tg_xxxx_xxyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 10); 

                auto tg_xxxx_xxyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 11); 

                auto tg_xxxx_xxyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 12); 

                auto tg_xxxx_xxyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 13); 

                auto tg_xxxx_xxzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 14); 

                auto tg_xxxx_xyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 15); 

                auto tg_xxxx_xyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 16); 

                auto tg_xxxx_xyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 17); 

                auto tg_xxxx_xyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 18); 

                auto tg_xxxx_xyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 19); 

                auto tg_xxxx_xzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 20); 

                auto tg_xxxx_yyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 21); 

                auto tg_xxxx_yyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 22); 

                auto tg_xxxx_yyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 23); 

                auto tg_xxxx_yyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 24); 

                auto tg_xxxx_yyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 25); 

                auto tg_xxxx_yzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 26); 

                auto tg_xxxx_zzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 27); 

                auto tg_xxxy_xxxxxx_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 28); 

                auto tg_xxxy_xxxxxy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 29); 

                auto tg_xxxy_xxxxxz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 30); 

                auto tg_xxxy_xxxxyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 31); 

                auto tg_xxxy_xxxxyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 32); 

                auto tg_xxxy_xxxxzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 33); 

                auto tg_xxxy_xxxyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 34); 

                auto tg_xxxy_xxxyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 35); 

                auto tg_xxxy_xxxyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 36); 

                auto tg_xxxy_xxxzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 37); 

                auto tg_xxxy_xxyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 38); 

                auto tg_xxxy_xxyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 39); 

                auto tg_xxxy_xxyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 40); 

                auto tg_xxxy_xxyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 41); 

                auto tg_xxxy_xxzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 42); 

                auto tg_xxxy_xyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 43); 

                auto tg_xxxy_xyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 44); 

                auto tg_xxxy_xyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 45); 

                auto tg_xxxy_xyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 46); 

                auto tg_xxxy_xyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 47); 

                auto tg_xxxy_xzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 48); 

                auto tg_xxxy_yyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 49); 

                auto tg_xxxy_yyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 50); 

                auto tg_xxxy_yyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 51); 

                auto tg_xxxy_yyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 52); 

                auto tg_xxxy_yyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 53); 

                auto tg_xxxy_yzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 54); 

                auto tg_xxxy_zzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 55); 

                auto tg_xxxz_xxxxxx_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 56); 

                auto tg_xxxz_xxxxxy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 57); 

                auto tg_xxxz_xxxxxz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 58); 

                auto tg_xxxz_xxxxyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 59); 

                auto tg_xxxz_xxxxyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 60); 

                auto tg_xxxz_xxxxzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 61); 

                auto tg_xxxz_xxxyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 62); 

                auto tg_xxxz_xxxyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 63); 

                auto tg_xxxz_xxxyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 64); 

                auto tg_xxxz_xxxzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 65); 

                auto tg_xxxz_xxyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 66); 

                auto tg_xxxz_xxyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 67); 

                auto tg_xxxz_xxyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 68); 

                auto tg_xxxz_xxyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 69); 

                auto tg_xxxz_xxzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 70); 

                auto tg_xxxz_xyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 71); 

                auto tg_xxxz_xyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 72); 

                auto tg_xxxz_xyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 73); 

                auto tg_xxxz_xyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 74); 

                auto tg_xxxz_xyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 75); 

                auto tg_xxxz_xzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 76); 

                auto tg_xxxz_yyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 77); 

                auto tg_xxxz_yyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 78); 

                auto tg_xxxz_yyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 79); 

                auto tg_xxxz_yyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 80); 

                auto tg_xxxz_yyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 81); 

                auto tg_xxxz_yzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 82); 

                auto tg_xxxz_zzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 83); 

                auto tg_xxyy_xxxxxx_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 84); 

                auto tg_xxyy_xxxxxy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 85); 

                auto tg_xxyy_xxxxxz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 86); 

                auto tg_xxyy_xxxxyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 87); 

                auto tg_xxyy_xxxxyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 88); 

                auto tg_xxyy_xxxxzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 89); 

                auto tg_xxyy_xxxyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 90); 

                auto tg_xxyy_xxxyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 91); 

                auto tg_xxyy_xxxyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 92); 

                auto tg_xxyy_xxxzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 93); 

                auto tg_xxyy_xxyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 94); 

                auto tg_xxyy_xxyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 95); 

                auto tg_xxyy_xxyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 96); 

                auto tg_xxyy_xxyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 97); 

                auto tg_xxxxx_xxxxx_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx); 

                auto tg_xxxxx_xxxxy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 1); 

                auto tg_xxxxx_xxxxz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 2); 

                auto tg_xxxxx_xxxyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 3); 

                auto tg_xxxxx_xxxyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 4); 

                auto tg_xxxxx_xxxzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 5); 

                auto tg_xxxxx_xxyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 6); 

                auto tg_xxxxx_xxyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 7); 

                auto tg_xxxxx_xxyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 8); 

                auto tg_xxxxx_xxzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 9); 

                auto tg_xxxxx_xyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 10); 

                auto tg_xxxxx_xyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 11); 

                auto tg_xxxxx_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 12); 

                auto tg_xxxxx_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 13); 

                auto tg_xxxxx_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 14); 

                auto tg_xxxxx_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 15); 

                auto tg_xxxxx_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 16); 

                auto tg_xxxxx_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 17); 

                auto tg_xxxxx_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 18); 

                auto tg_xxxxx_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 19); 

                auto tg_xxxxx_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 20); 

                auto tg_xxxxy_xxxxx_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 21); 

                auto tg_xxxxy_xxxxy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 22); 

                auto tg_xxxxy_xxxxz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 23); 

                auto tg_xxxxy_xxxyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 24); 

                auto tg_xxxxy_xxxyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 25); 

                auto tg_xxxxy_xxxzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 26); 

                auto tg_xxxxy_xxyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 27); 

                auto tg_xxxxy_xxyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 28); 

                auto tg_xxxxy_xxyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 29); 

                auto tg_xxxxy_xxzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 30); 

                auto tg_xxxxy_xyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 31); 

                auto tg_xxxxy_xyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 32); 

                auto tg_xxxxy_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 33); 

                auto tg_xxxxy_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 34); 

                auto tg_xxxxy_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 35); 

                auto tg_xxxxy_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 36); 

                auto tg_xxxxy_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 37); 

                auto tg_xxxxy_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 38); 

                auto tg_xxxxy_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 39); 

                auto tg_xxxxy_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 40); 

                auto tg_xxxxy_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 41); 

                auto tg_xxxxz_xxxxx_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 42); 

                auto tg_xxxxz_xxxxy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 43); 

                auto tg_xxxxz_xxxxz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 44); 

                auto tg_xxxxz_xxxyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 45); 

                auto tg_xxxxz_xxxyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 46); 

                auto tg_xxxxz_xxxzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 47); 

                auto tg_xxxxz_xxyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 48); 

                auto tg_xxxxz_xxyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 49); 

                auto tg_xxxxz_xxyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 50); 

                auto tg_xxxxz_xxzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 51); 

                auto tg_xxxxz_xyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 52); 

                auto tg_xxxxz_xyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 53); 

                auto tg_xxxxz_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 54); 

                auto tg_xxxxz_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 55); 

                auto tg_xxxxz_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 56); 

                auto tg_xxxxz_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 57); 

                auto tg_xxxxz_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 58); 

                auto tg_xxxxz_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 59); 

                auto tg_xxxxz_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 60); 

                auto tg_xxxxz_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 61); 

                auto tg_xxxxz_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 62); 

                auto tg_xxxyy_xxxxx_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 63); 

                auto tg_xxxyy_xxxxy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 64); 

                auto tg_xxxyy_xxxxz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 65); 

                auto tg_xxxyy_xxxyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 66); 

                auto tg_xxxyy_xxxyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 67); 

                auto tg_xxxyy_xxxzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 68); 

                auto tg_xxxyy_xxyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 69); 

                auto tg_xxxyy_xxyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 70); 

                auto tg_xxxyy_xxyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 71); 

                auto tg_xxxyy_xxzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 72); 

                auto tg_xxxyy_xyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 73); 

                auto tg_xxxyy_xyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 74); 

                auto tg_xxxyy_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 75); 

                auto tg_xxxyy_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 76); 

                // set up pointers to integrals

                auto tg_xxxxxx_xxxxxx_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx); 

                auto tg_xxxxxx_xxxxxy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 1); 

                auto tg_xxxxxx_xxxxxz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 2); 

                auto tg_xxxxxx_xxxxyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 3); 

                auto tg_xxxxxx_xxxxyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 4); 

                auto tg_xxxxxx_xxxxzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 5); 

                auto tg_xxxxxx_xxxyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 6); 

                auto tg_xxxxxx_xxxyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 7); 

                auto tg_xxxxxx_xxxyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 8); 

                auto tg_xxxxxx_xxxzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 9); 

                auto tg_xxxxxx_xxyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 10); 

                auto tg_xxxxxx_xxyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 11); 

                auto tg_xxxxxx_xxyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 12); 

                auto tg_xxxxxx_xxyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 13); 

                auto tg_xxxxxx_xxzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 14); 

                auto tg_xxxxxx_xyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 15); 

                auto tg_xxxxxx_xyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 16); 

                auto tg_xxxxxx_xyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 17); 

                auto tg_xxxxxx_xyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 18); 

                auto tg_xxxxxx_xyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 19); 

                auto tg_xxxxxx_xzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 20); 

                auto tg_xxxxxx_yyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 21); 

                auto tg_xxxxxx_yyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 22); 

                auto tg_xxxxxx_yyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 23); 

                auto tg_xxxxxx_yyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 24); 

                auto tg_xxxxxx_yyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 25); 

                auto tg_xxxxxx_yzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 26); 

                auto tg_xxxxxx_zzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 27); 

                auto tg_xxxxxy_xxxxxx_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 28); 

                auto tg_xxxxxy_xxxxxy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 29); 

                auto tg_xxxxxy_xxxxxz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 30); 

                auto tg_xxxxxy_xxxxyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 31); 

                auto tg_xxxxxy_xxxxyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 32); 

                auto tg_xxxxxy_xxxxzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 33); 

                auto tg_xxxxxy_xxxyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 34); 

                auto tg_xxxxxy_xxxyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 35); 

                auto tg_xxxxxy_xxxyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 36); 

                auto tg_xxxxxy_xxxzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 37); 

                auto tg_xxxxxy_xxyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 38); 

                auto tg_xxxxxy_xxyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 39); 

                auto tg_xxxxxy_xxyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 40); 

                auto tg_xxxxxy_xxyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 41); 

                auto tg_xxxxxy_xxzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 42); 

                auto tg_xxxxxy_xyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 43); 

                auto tg_xxxxxy_xyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 44); 

                auto tg_xxxxxy_xyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 45); 

                auto tg_xxxxxy_xyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 46); 

                auto tg_xxxxxy_xyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 47); 

                auto tg_xxxxxy_xzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 48); 

                auto tg_xxxxxy_yyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 49); 

                auto tg_xxxxxy_yyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 50); 

                auto tg_xxxxxy_yyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 51); 

                auto tg_xxxxxy_yyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 52); 

                auto tg_xxxxxy_yyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 53); 

                auto tg_xxxxxy_yzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 54); 

                auto tg_xxxxxy_zzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 55); 

                auto tg_xxxxxz_xxxxxx_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 56); 

                auto tg_xxxxxz_xxxxxy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 57); 

                auto tg_xxxxxz_xxxxxz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 58); 

                auto tg_xxxxxz_xxxxyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 59); 

                auto tg_xxxxxz_xxxxyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 60); 

                auto tg_xxxxxz_xxxxzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 61); 

                auto tg_xxxxxz_xxxyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 62); 

                auto tg_xxxxxz_xxxyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 63); 

                auto tg_xxxxxz_xxxyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 64); 

                auto tg_xxxxxz_xxxzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 65); 

                auto tg_xxxxxz_xxyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 66); 

                auto tg_xxxxxz_xxyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 67); 

                auto tg_xxxxxz_xxyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 68); 

                auto tg_xxxxxz_xxyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 69); 

                auto tg_xxxxxz_xxzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 70); 

                auto tg_xxxxxz_xyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 71); 

                auto tg_xxxxxz_xyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 72); 

                auto tg_xxxxxz_xyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 73); 

                auto tg_xxxxxz_xyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 74); 

                auto tg_xxxxxz_xyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 75); 

                auto tg_xxxxxz_xzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 76); 

                auto tg_xxxxxz_yyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 77); 

                auto tg_xxxxxz_yyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 78); 

                auto tg_xxxxxz_yyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 79); 

                auto tg_xxxxxz_yyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 80); 

                auto tg_xxxxxz_yyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 81); 

                auto tg_xxxxxz_yzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 82); 

                auto tg_xxxxxz_zzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 83); 

                auto tg_xxxxyy_xxxxxx_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 84); 

                auto tg_xxxxyy_xxxxxy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 85); 

                auto tg_xxxxyy_xxxxxz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 86); 

                auto tg_xxxxyy_xxxxyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 87); 

                auto tg_xxxxyy_xxxxyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 88); 

                auto tg_xxxxyy_xxxxzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 89); 

                auto tg_xxxxyy_xxxyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 90); 

                auto tg_xxxxyy_xxxyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 91); 

                auto tg_xxxxyy_xxxyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 92); 

                auto tg_xxxxyy_xxxzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 93); 

                auto tg_xxxxyy_xxyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 94); 

                auto tg_xxxxyy_xxyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 95); 

                auto tg_xxxxyy_xxyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 96); 

                auto tg_xxxxyy_xxyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 97); 

                // Batch of Integrals (0,98)

                #pragma omp simd aligned(fxn, fza, tg_xxxx_xxxxxx_0, tg_xxxx_xxxxxx_1, tg_xxxx_xxxxxy_0, \
                                         tg_xxxx_xxxxxy_1, tg_xxxx_xxxxxz_0, tg_xxxx_xxxxxz_1, tg_xxxx_xxxxyy_0, \
                                         tg_xxxx_xxxxyy_1, tg_xxxx_xxxxyz_0, tg_xxxx_xxxxyz_1, tg_xxxx_xxxxzz_0, \
                                         tg_xxxx_xxxxzz_1, tg_xxxx_xxxyyy_0, tg_xxxx_xxxyyy_1, tg_xxxx_xxxyyz_0, \
                                         tg_xxxx_xxxyyz_1, tg_xxxx_xxxyzz_0, tg_xxxx_xxxyzz_1, tg_xxxx_xxxzzz_0, \
                                         tg_xxxx_xxxzzz_1, tg_xxxx_xxyyyy_0, tg_xxxx_xxyyyy_1, tg_xxxx_xxyyyz_0, \
                                         tg_xxxx_xxyyyz_1, tg_xxxx_xxyyzz_0, tg_xxxx_xxyyzz_1, tg_xxxx_xxyzzz_0, \
                                         tg_xxxx_xxyzzz_1, tg_xxxx_xxzzzz_0, tg_xxxx_xxzzzz_1, tg_xxxx_xyyyyy_0, \
                                         tg_xxxx_xyyyyy_1, tg_xxxx_xyyyyz_0, tg_xxxx_xyyyyz_1, tg_xxxx_xyyyzz_0, \
                                         tg_xxxx_xyyyzz_1, tg_xxxx_xyyzzz_0, tg_xxxx_xyyzzz_1, tg_xxxx_xyzzzz_0, \
                                         tg_xxxx_xyzzzz_1, tg_xxxx_xzzzzz_0, tg_xxxx_xzzzzz_1, tg_xxxx_yyyyyy_0, \
                                         tg_xxxx_yyyyyy_1, tg_xxxx_yyyyyz_0, tg_xxxx_yyyyyz_1, tg_xxxx_yyyyzz_0, \
                                         tg_xxxx_yyyyzz_1, tg_xxxx_yyyzzz_0, tg_xxxx_yyyzzz_1, tg_xxxx_yyzzzz_0, \
                                         tg_xxxx_yyzzzz_1, tg_xxxx_yzzzzz_0, tg_xxxx_yzzzzz_1, tg_xxxx_zzzzzz_0, \
                                         tg_xxxx_zzzzzz_1, tg_xxxxx_xxxxx_1, tg_xxxxx_xxxxxx_0, tg_xxxxx_xxxxxx_1, \
                                         tg_xxxxx_xxxxxy_0, tg_xxxxx_xxxxxy_1, tg_xxxxx_xxxxxz_0, tg_xxxxx_xxxxxz_1, \
                                         tg_xxxxx_xxxxy_1, tg_xxxxx_xxxxyy_0, tg_xxxxx_xxxxyy_1, tg_xxxxx_xxxxyz_0, \
                                         tg_xxxxx_xxxxyz_1, tg_xxxxx_xxxxz_1, tg_xxxxx_xxxxzz_0, tg_xxxxx_xxxxzz_1, \
                                         tg_xxxxx_xxxyy_1, tg_xxxxx_xxxyyy_0, tg_xxxxx_xxxyyy_1, tg_xxxxx_xxxyyz_0, \
                                         tg_xxxxx_xxxyyz_1, tg_xxxxx_xxxyz_1, tg_xxxxx_xxxyzz_0, tg_xxxxx_xxxyzz_1, \
                                         tg_xxxxx_xxxzz_1, tg_xxxxx_xxxzzz_0, tg_xxxxx_xxxzzz_1, tg_xxxxx_xxyyy_1, \
                                         tg_xxxxx_xxyyyy_0, tg_xxxxx_xxyyyy_1, tg_xxxxx_xxyyyz_0, tg_xxxxx_xxyyyz_1, \
                                         tg_xxxxx_xxyyz_1, tg_xxxxx_xxyyzz_0, tg_xxxxx_xxyyzz_1, tg_xxxxx_xxyzz_1, \
                                         tg_xxxxx_xxyzzz_0, tg_xxxxx_xxyzzz_1, tg_xxxxx_xxzzz_1, tg_xxxxx_xxzzzz_0, \
                                         tg_xxxxx_xxzzzz_1, tg_xxxxx_xyyyy_1, tg_xxxxx_xyyyyy_0, tg_xxxxx_xyyyyy_1, \
                                         tg_xxxxx_xyyyyz_0, tg_xxxxx_xyyyyz_1, tg_xxxxx_xyyyz_1, tg_xxxxx_xyyyzz_0, \
                                         tg_xxxxx_xyyyzz_1, tg_xxxxx_xyyzz_1, tg_xxxxx_xyyzzz_0, tg_xxxxx_xyyzzz_1, \
                                         tg_xxxxx_xyzzz_1, tg_xxxxx_xyzzzz_0, tg_xxxxx_xyzzzz_1, tg_xxxxx_xzzzz_1, \
                                         tg_xxxxx_xzzzzz_0, tg_xxxxx_xzzzzz_1, tg_xxxxx_yyyyy_1, tg_xxxxx_yyyyyy_0, \
                                         tg_xxxxx_yyyyyy_1, tg_xxxxx_yyyyyz_0, tg_xxxxx_yyyyyz_1, tg_xxxxx_yyyyz_1, \
                                         tg_xxxxx_yyyyzz_0, tg_xxxxx_yyyyzz_1, tg_xxxxx_yyyzz_1, tg_xxxxx_yyyzzz_0, \
                                         tg_xxxxx_yyyzzz_1, tg_xxxxx_yyzzz_1, tg_xxxxx_yyzzzz_0, tg_xxxxx_yyzzzz_1, \
                                         tg_xxxxx_yzzzz_1, tg_xxxxx_yzzzzz_0, tg_xxxxx_yzzzzz_1, tg_xxxxx_zzzzz_1, \
                                         tg_xxxxx_zzzzzz_0, tg_xxxxx_zzzzzz_1, tg_xxxxxx_xxxxxx_0, tg_xxxxxx_xxxxxy_0, \
                                         tg_xxxxxx_xxxxxz_0, tg_xxxxxx_xxxxyy_0, tg_xxxxxx_xxxxyz_0, tg_xxxxxx_xxxxzz_0, \
                                         tg_xxxxxx_xxxyyy_0, tg_xxxxxx_xxxyyz_0, tg_xxxxxx_xxxyzz_0, tg_xxxxxx_xxxzzz_0, \
                                         tg_xxxxxx_xxyyyy_0, tg_xxxxxx_xxyyyz_0, tg_xxxxxx_xxyyzz_0, tg_xxxxxx_xxyzzz_0, \
                                         tg_xxxxxx_xxzzzz_0, tg_xxxxxx_xyyyyy_0, tg_xxxxxx_xyyyyz_0, tg_xxxxxx_xyyyzz_0, \
                                         tg_xxxxxx_xyyzzz_0, tg_xxxxxx_xyzzzz_0, tg_xxxxxx_xzzzzz_0, tg_xxxxxx_yyyyyy_0, \
                                         tg_xxxxxx_yyyyyz_0, tg_xxxxxx_yyyyzz_0, tg_xxxxxx_yyyzzz_0, tg_xxxxxx_yyzzzz_0, \
                                         tg_xxxxxx_yzzzzz_0, tg_xxxxxx_zzzzzz_0, tg_xxxxxy_xxxxxx_0, tg_xxxxxy_xxxxxy_0, \
                                         tg_xxxxxy_xxxxxz_0, tg_xxxxxy_xxxxyy_0, tg_xxxxxy_xxxxyz_0, tg_xxxxxy_xxxxzz_0, \
                                         tg_xxxxxy_xxxyyy_0, tg_xxxxxy_xxxyyz_0, tg_xxxxxy_xxxyzz_0, tg_xxxxxy_xxxzzz_0, \
                                         tg_xxxxxy_xxyyyy_0, tg_xxxxxy_xxyyyz_0, tg_xxxxxy_xxyyzz_0, tg_xxxxxy_xxyzzz_0, \
                                         tg_xxxxxy_xxzzzz_0, tg_xxxxxy_xyyyyy_0, tg_xxxxxy_xyyyyz_0, tg_xxxxxy_xyyyzz_0, \
                                         tg_xxxxxy_xyyzzz_0, tg_xxxxxy_xyzzzz_0, tg_xxxxxy_xzzzzz_0, tg_xxxxxy_yyyyyy_0, \
                                         tg_xxxxxy_yyyyyz_0, tg_xxxxxy_yyyyzz_0, tg_xxxxxy_yyyzzz_0, tg_xxxxxy_yyzzzz_0, \
                                         tg_xxxxxy_yzzzzz_0, tg_xxxxxy_zzzzzz_0, tg_xxxxxz_xxxxxx_0, tg_xxxxxz_xxxxxy_0, \
                                         tg_xxxxxz_xxxxxz_0, tg_xxxxxz_xxxxyy_0, tg_xxxxxz_xxxxyz_0, tg_xxxxxz_xxxxzz_0, \
                                         tg_xxxxxz_xxxyyy_0, tg_xxxxxz_xxxyyz_0, tg_xxxxxz_xxxyzz_0, tg_xxxxxz_xxxzzz_0, \
                                         tg_xxxxxz_xxyyyy_0, tg_xxxxxz_xxyyyz_0, tg_xxxxxz_xxyyzz_0, tg_xxxxxz_xxyzzz_0, \
                                         tg_xxxxxz_xxzzzz_0, tg_xxxxxz_xyyyyy_0, tg_xxxxxz_xyyyyz_0, tg_xxxxxz_xyyyzz_0, \
                                         tg_xxxxxz_xyyzzz_0, tg_xxxxxz_xyzzzz_0, tg_xxxxxz_xzzzzz_0, tg_xxxxxz_yyyyyy_0, \
                                         tg_xxxxxz_yyyyyz_0, tg_xxxxxz_yyyyzz_0, tg_xxxxxz_yyyzzz_0, tg_xxxxxz_yyzzzz_0, \
                                         tg_xxxxxz_yzzzzz_0, tg_xxxxxz_zzzzzz_0, tg_xxxxy_xxxxx_1, tg_xxxxy_xxxxxx_0, \
                                         tg_xxxxy_xxxxxx_1, tg_xxxxy_xxxxxy_0, tg_xxxxy_xxxxxy_1, tg_xxxxy_xxxxxz_0, \
                                         tg_xxxxy_xxxxxz_1, tg_xxxxy_xxxxy_1, tg_xxxxy_xxxxyy_0, tg_xxxxy_xxxxyy_1, \
                                         tg_xxxxy_xxxxyz_0, tg_xxxxy_xxxxyz_1, tg_xxxxy_xxxxz_1, tg_xxxxy_xxxxzz_0, \
                                         tg_xxxxy_xxxxzz_1, tg_xxxxy_xxxyy_1, tg_xxxxy_xxxyyy_0, tg_xxxxy_xxxyyy_1, \
                                         tg_xxxxy_xxxyyz_0, tg_xxxxy_xxxyyz_1, tg_xxxxy_xxxyz_1, tg_xxxxy_xxxyzz_0, \
                                         tg_xxxxy_xxxyzz_1, tg_xxxxy_xxxzz_1, tg_xxxxy_xxxzzz_0, tg_xxxxy_xxxzzz_1, \
                                         tg_xxxxy_xxyyy_1, tg_xxxxy_xxyyyy_0, tg_xxxxy_xxyyyy_1, tg_xxxxy_xxyyyz_0, \
                                         tg_xxxxy_xxyyyz_1, tg_xxxxy_xxyyz_1, tg_xxxxy_xxyyzz_0, tg_xxxxy_xxyyzz_1, \
                                         tg_xxxxy_xxyzz_1, tg_xxxxy_xxyzzz_0, tg_xxxxy_xxyzzz_1, tg_xxxxy_xxzzz_1, \
                                         tg_xxxxy_xxzzzz_0, tg_xxxxy_xxzzzz_1, tg_xxxxy_xyyyy_1, tg_xxxxy_xyyyyy_0, \
                                         tg_xxxxy_xyyyyy_1, tg_xxxxy_xyyyyz_0, tg_xxxxy_xyyyyz_1, tg_xxxxy_xyyyz_1, \
                                         tg_xxxxy_xyyyzz_0, tg_xxxxy_xyyyzz_1, tg_xxxxy_xyyzz_1, tg_xxxxy_xyyzzz_0, \
                                         tg_xxxxy_xyyzzz_1, tg_xxxxy_xyzzz_1, tg_xxxxy_xyzzzz_0, tg_xxxxy_xyzzzz_1, \
                                         tg_xxxxy_xzzzz_1, tg_xxxxy_xzzzzz_0, tg_xxxxy_xzzzzz_1, tg_xxxxy_yyyyy_1, \
                                         tg_xxxxy_yyyyyy_0, tg_xxxxy_yyyyyy_1, tg_xxxxy_yyyyyz_0, tg_xxxxy_yyyyyz_1, \
                                         tg_xxxxy_yyyyz_1, tg_xxxxy_yyyyzz_0, tg_xxxxy_yyyyzz_1, tg_xxxxy_yyyzz_1, \
                                         tg_xxxxy_yyyzzz_0, tg_xxxxy_yyyzzz_1, tg_xxxxy_yyzzz_1, tg_xxxxy_yyzzzz_0, \
                                         tg_xxxxy_yyzzzz_1, tg_xxxxy_yzzzz_1, tg_xxxxy_yzzzzz_0, tg_xxxxy_yzzzzz_1, \
                                         tg_xxxxy_zzzzz_1, tg_xxxxy_zzzzzz_0, tg_xxxxy_zzzzzz_1, tg_xxxxyy_xxxxxx_0, \
                                         tg_xxxxyy_xxxxxy_0, tg_xxxxyy_xxxxxz_0, tg_xxxxyy_xxxxyy_0, tg_xxxxyy_xxxxyz_0, \
                                         tg_xxxxyy_xxxxzz_0, tg_xxxxyy_xxxyyy_0, tg_xxxxyy_xxxyyz_0, tg_xxxxyy_xxxyzz_0, \
                                         tg_xxxxyy_xxxzzz_0, tg_xxxxyy_xxyyyy_0, tg_xxxxyy_xxyyyz_0, tg_xxxxyy_xxyyzz_0, \
                                         tg_xxxxyy_xxyzzz_0, tg_xxxxz_xxxxx_1, tg_xxxxz_xxxxxx_0, tg_xxxxz_xxxxxx_1, \
                                         tg_xxxxz_xxxxxy_0, tg_xxxxz_xxxxxy_1, tg_xxxxz_xxxxxz_0, tg_xxxxz_xxxxxz_1, \
                                         tg_xxxxz_xxxxy_1, tg_xxxxz_xxxxyy_0, tg_xxxxz_xxxxyy_1, tg_xxxxz_xxxxyz_0, \
                                         tg_xxxxz_xxxxyz_1, tg_xxxxz_xxxxz_1, tg_xxxxz_xxxxzz_0, tg_xxxxz_xxxxzz_1, \
                                         tg_xxxxz_xxxyy_1, tg_xxxxz_xxxyyy_0, tg_xxxxz_xxxyyy_1, tg_xxxxz_xxxyyz_0, \
                                         tg_xxxxz_xxxyyz_1, tg_xxxxz_xxxyz_1, tg_xxxxz_xxxyzz_0, tg_xxxxz_xxxyzz_1, \
                                         tg_xxxxz_xxxzz_1, tg_xxxxz_xxxzzz_0, tg_xxxxz_xxxzzz_1, tg_xxxxz_xxyyy_1, \
                                         tg_xxxxz_xxyyyy_0, tg_xxxxz_xxyyyy_1, tg_xxxxz_xxyyyz_0, tg_xxxxz_xxyyyz_1, \
                                         tg_xxxxz_xxyyz_1, tg_xxxxz_xxyyzz_0, tg_xxxxz_xxyyzz_1, tg_xxxxz_xxyzz_1, \
                                         tg_xxxxz_xxyzzz_0, tg_xxxxz_xxyzzz_1, tg_xxxxz_xxzzz_1, tg_xxxxz_xxzzzz_0, \
                                         tg_xxxxz_xxzzzz_1, tg_xxxxz_xyyyy_1, tg_xxxxz_xyyyyy_0, tg_xxxxz_xyyyyy_1, \
                                         tg_xxxxz_xyyyyz_0, tg_xxxxz_xyyyyz_1, tg_xxxxz_xyyyz_1, tg_xxxxz_xyyyzz_0, \
                                         tg_xxxxz_xyyyzz_1, tg_xxxxz_xyyzz_1, tg_xxxxz_xyyzzz_0, tg_xxxxz_xyyzzz_1, \
                                         tg_xxxxz_xyzzz_1, tg_xxxxz_xyzzzz_0, tg_xxxxz_xyzzzz_1, tg_xxxxz_xzzzz_1, \
                                         tg_xxxxz_xzzzzz_0, tg_xxxxz_xzzzzz_1, tg_xxxxz_yyyyy_1, tg_xxxxz_yyyyyy_0, \
                                         tg_xxxxz_yyyyyy_1, tg_xxxxz_yyyyyz_0, tg_xxxxz_yyyyyz_1, tg_xxxxz_yyyyz_1, \
                                         tg_xxxxz_yyyyzz_0, tg_xxxxz_yyyyzz_1, tg_xxxxz_yyyzz_1, tg_xxxxz_yyyzzz_0, \
                                         tg_xxxxz_yyyzzz_1, tg_xxxxz_yyzzz_1, tg_xxxxz_yyzzzz_0, tg_xxxxz_yyzzzz_1, \
                                         tg_xxxxz_yzzzz_1, tg_xxxxz_yzzzzz_0, tg_xxxxz_yzzzzz_1, tg_xxxxz_zzzzz_1, \
                                         tg_xxxxz_zzzzzz_0, tg_xxxxz_zzzzzz_1, tg_xxxy_xxxxxx_0, tg_xxxy_xxxxxx_1, \
                                         tg_xxxy_xxxxxy_0, tg_xxxy_xxxxxy_1, tg_xxxy_xxxxxz_0, tg_xxxy_xxxxxz_1, \
                                         tg_xxxy_xxxxyy_0, tg_xxxy_xxxxyy_1, tg_xxxy_xxxxyz_0, tg_xxxy_xxxxyz_1, \
                                         tg_xxxy_xxxxzz_0, tg_xxxy_xxxxzz_1, tg_xxxy_xxxyyy_0, tg_xxxy_xxxyyy_1, \
                                         tg_xxxy_xxxyyz_0, tg_xxxy_xxxyyz_1, tg_xxxy_xxxyzz_0, tg_xxxy_xxxyzz_1, \
                                         tg_xxxy_xxxzzz_0, tg_xxxy_xxxzzz_1, tg_xxxy_xxyyyy_0, tg_xxxy_xxyyyy_1, \
                                         tg_xxxy_xxyyyz_0, tg_xxxy_xxyyyz_1, tg_xxxy_xxyyzz_0, tg_xxxy_xxyyzz_1, \
                                         tg_xxxy_xxyzzz_0, tg_xxxy_xxyzzz_1, tg_xxxy_xxzzzz_0, tg_xxxy_xxzzzz_1, \
                                         tg_xxxy_xyyyyy_0, tg_xxxy_xyyyyy_1, tg_xxxy_xyyyyz_0, tg_xxxy_xyyyyz_1, \
                                         tg_xxxy_xyyyzz_0, tg_xxxy_xyyyzz_1, tg_xxxy_xyyzzz_0, tg_xxxy_xyyzzz_1, \
                                         tg_xxxy_xyzzzz_0, tg_xxxy_xyzzzz_1, tg_xxxy_xzzzzz_0, tg_xxxy_xzzzzz_1, \
                                         tg_xxxy_yyyyyy_0, tg_xxxy_yyyyyy_1, tg_xxxy_yyyyyz_0, tg_xxxy_yyyyyz_1, \
                                         tg_xxxy_yyyyzz_0, tg_xxxy_yyyyzz_1, tg_xxxy_yyyzzz_0, tg_xxxy_yyyzzz_1, \
                                         tg_xxxy_yyzzzz_0, tg_xxxy_yyzzzz_1, tg_xxxy_yzzzzz_0, tg_xxxy_yzzzzz_1, \
                                         tg_xxxy_zzzzzz_0, tg_xxxy_zzzzzz_1, tg_xxxyy_xxxxx_1, tg_xxxyy_xxxxxx_0, \
                                         tg_xxxyy_xxxxxx_1, tg_xxxyy_xxxxxy_0, tg_xxxyy_xxxxxy_1, tg_xxxyy_xxxxxz_0, \
                                         tg_xxxyy_xxxxxz_1, tg_xxxyy_xxxxy_1, tg_xxxyy_xxxxyy_0, tg_xxxyy_xxxxyy_1, \
                                         tg_xxxyy_xxxxyz_0, tg_xxxyy_xxxxyz_1, tg_xxxyy_xxxxz_1, tg_xxxyy_xxxxzz_0, \
                                         tg_xxxyy_xxxxzz_1, tg_xxxyy_xxxyy_1, tg_xxxyy_xxxyyy_0, tg_xxxyy_xxxyyy_1, \
                                         tg_xxxyy_xxxyyz_0, tg_xxxyy_xxxyyz_1, tg_xxxyy_xxxyz_1, tg_xxxyy_xxxyzz_0, \
                                         tg_xxxyy_xxxyzz_1, tg_xxxyy_xxxzz_1, tg_xxxyy_xxxzzz_0, tg_xxxyy_xxxzzz_1, \
                                         tg_xxxyy_xxyyy_1, tg_xxxyy_xxyyyy_0, tg_xxxyy_xxyyyy_1, tg_xxxyy_xxyyyz_0, \
                                         tg_xxxyy_xxyyyz_1, tg_xxxyy_xxyyz_1, tg_xxxyy_xxyyzz_0, tg_xxxyy_xxyyzz_1, \
                                         tg_xxxyy_xxyzz_1, tg_xxxyy_xxyzzz_0, tg_xxxyy_xxyzzz_1, tg_xxxyy_xxzzz_1, \
                                         tg_xxxyy_xyyyy_1, tg_xxxyy_xyyyz_1, tg_xxxyy_xyyzz_1, tg_xxxyy_xyzzz_1, \
                                         tg_xxxz_xxxxxx_0, tg_xxxz_xxxxxx_1, tg_xxxz_xxxxxy_0, tg_xxxz_xxxxxy_1, \
                                         tg_xxxz_xxxxxz_0, tg_xxxz_xxxxxz_1, tg_xxxz_xxxxyy_0, tg_xxxz_xxxxyy_1, \
                                         tg_xxxz_xxxxyz_0, tg_xxxz_xxxxyz_1, tg_xxxz_xxxxzz_0, tg_xxxz_xxxxzz_1, \
                                         tg_xxxz_xxxyyy_0, tg_xxxz_xxxyyy_1, tg_xxxz_xxxyyz_0, tg_xxxz_xxxyyz_1, \
                                         tg_xxxz_xxxyzz_0, tg_xxxz_xxxyzz_1, tg_xxxz_xxxzzz_0, tg_xxxz_xxxzzz_1, \
                                         tg_xxxz_xxyyyy_0, tg_xxxz_xxyyyy_1, tg_xxxz_xxyyyz_0, tg_xxxz_xxyyyz_1, \
                                         tg_xxxz_xxyyzz_0, tg_xxxz_xxyyzz_1, tg_xxxz_xxyzzz_0, tg_xxxz_xxyzzz_1, \
                                         tg_xxxz_xxzzzz_0, tg_xxxz_xxzzzz_1, tg_xxxz_xyyyyy_0, tg_xxxz_xyyyyy_1, \
                                         tg_xxxz_xyyyyz_0, tg_xxxz_xyyyyz_1, tg_xxxz_xyyyzz_0, tg_xxxz_xyyyzz_1, \
                                         tg_xxxz_xyyzzz_0, tg_xxxz_xyyzzz_1, tg_xxxz_xyzzzz_0, tg_xxxz_xyzzzz_1, \
                                         tg_xxxz_xzzzzz_0, tg_xxxz_xzzzzz_1, tg_xxxz_yyyyyy_0, tg_xxxz_yyyyyy_1, \
                                         tg_xxxz_yyyyyz_0, tg_xxxz_yyyyyz_1, tg_xxxz_yyyyzz_0, tg_xxxz_yyyyzz_1, \
                                         tg_xxxz_yyyzzz_0, tg_xxxz_yyyzzz_1, tg_xxxz_yyzzzz_0, tg_xxxz_yyzzzz_1, \
                                         tg_xxxz_yzzzzz_0, tg_xxxz_yzzzzz_1, tg_xxxz_zzzzzz_0, tg_xxxz_zzzzzz_1, \
                                         tg_xxyy_xxxxxx_0, tg_xxyy_xxxxxx_1, tg_xxyy_xxxxxy_0, tg_xxyy_xxxxxy_1, \
                                         tg_xxyy_xxxxxz_0, tg_xxyy_xxxxxz_1, tg_xxyy_xxxxyy_0, tg_xxyy_xxxxyy_1, \
                                         tg_xxyy_xxxxyz_0, tg_xxyy_xxxxyz_1, tg_xxyy_xxxxzz_0, tg_xxyy_xxxxzz_1, \
                                         tg_xxyy_xxxyyy_0, tg_xxyy_xxxyyy_1, tg_xxyy_xxxyyz_0, tg_xxyy_xxxyyz_1, \
                                         tg_xxyy_xxxyzz_0, tg_xxyy_xxxyzz_1, tg_xxyy_xxxzzz_0, tg_xxyy_xxxzzz_1, \
                                         tg_xxyy_xxyyyy_0, tg_xxyy_xxyyyy_1, tg_xxyy_xxyyyz_0, tg_xxyy_xxyyyz_1, \
                                         tg_xxyy_xxyyzz_0, tg_xxyy_xxyyzz_1, tg_xxyy_xxyzzz_0, tg_xxyy_xxyzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxxxx_xxxxxx_0[j] = pb_x * tg_xxxxx_xxxxxx_0[j] + fr * tg_xxxxx_xxxxxx_1[j] + 2.5 * fl1_fx * (tg_xxxx_xxxxxx_0[j] - tg_xxxx_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxxxx_xxxxx_1[j];

                    tg_xxxxxx_xxxxxy_0[j] = pb_x * tg_xxxxx_xxxxxy_0[j] + fr * tg_xxxxx_xxxxxy_1[j] + 2.5 * fl1_fx * (tg_xxxx_xxxxxy_0[j] - tg_xxxx_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxx_xxxxy_1[j];

                    tg_xxxxxx_xxxxxz_0[j] = pb_x * tg_xxxxx_xxxxxz_0[j] + fr * tg_xxxxx_xxxxxz_1[j] + 2.5 * fl1_fx * (tg_xxxx_xxxxxz_0[j] - tg_xxxx_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxx_xxxxz_1[j];

                    tg_xxxxxx_xxxxyy_0[j] = pb_x * tg_xxxxx_xxxxyy_0[j] + fr * tg_xxxxx_xxxxyy_1[j] + 2.5 * fl1_fx * (tg_xxxx_xxxxyy_0[j] - tg_xxxx_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxx_xxxyy_1[j];

                    tg_xxxxxx_xxxxyz_0[j] = pb_x * tg_xxxxx_xxxxyz_0[j] + fr * tg_xxxxx_xxxxyz_1[j] + 2.5 * fl1_fx * (tg_xxxx_xxxxyz_0[j] - tg_xxxx_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxx_xxxyz_1[j];

                    tg_xxxxxx_xxxxzz_0[j] = pb_x * tg_xxxxx_xxxxzz_0[j] + fr * tg_xxxxx_xxxxzz_1[j] + 2.5 * fl1_fx * (tg_xxxx_xxxxzz_0[j] - tg_xxxx_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxx_xxxzz_1[j];

                    tg_xxxxxx_xxxyyy_0[j] = pb_x * tg_xxxxx_xxxyyy_0[j] + fr * tg_xxxxx_xxxyyy_1[j] + 2.5 * fl1_fx * (tg_xxxx_xxxyyy_0[j] - tg_xxxx_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxx_xxyyy_1[j];

                    tg_xxxxxx_xxxyyz_0[j] = pb_x * tg_xxxxx_xxxyyz_0[j] + fr * tg_xxxxx_xxxyyz_1[j] + 2.5 * fl1_fx * (tg_xxxx_xxxyyz_0[j] - tg_xxxx_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxx_xxyyz_1[j];

                    tg_xxxxxx_xxxyzz_0[j] = pb_x * tg_xxxxx_xxxyzz_0[j] + fr * tg_xxxxx_xxxyzz_1[j] + 2.5 * fl1_fx * (tg_xxxx_xxxyzz_0[j] - tg_xxxx_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxx_xxyzz_1[j];

                    tg_xxxxxx_xxxzzz_0[j] = pb_x * tg_xxxxx_xxxzzz_0[j] + fr * tg_xxxxx_xxxzzz_1[j] + 2.5 * fl1_fx * (tg_xxxx_xxxzzz_0[j] - tg_xxxx_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxx_xxzzz_1[j];

                    tg_xxxxxx_xxyyyy_0[j] = pb_x * tg_xxxxx_xxyyyy_0[j] + fr * tg_xxxxx_xxyyyy_1[j] + 2.5 * fl1_fx * (tg_xxxx_xxyyyy_0[j] - tg_xxxx_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxx_xyyyy_1[j];

                    tg_xxxxxx_xxyyyz_0[j] = pb_x * tg_xxxxx_xxyyyz_0[j] + fr * tg_xxxxx_xxyyyz_1[j] + 2.5 * fl1_fx * (tg_xxxx_xxyyyz_0[j] - tg_xxxx_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxx_xyyyz_1[j];

                    tg_xxxxxx_xxyyzz_0[j] = pb_x * tg_xxxxx_xxyyzz_0[j] + fr * tg_xxxxx_xxyyzz_1[j] + 2.5 * fl1_fx * (tg_xxxx_xxyyzz_0[j] - tg_xxxx_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxx_xyyzz_1[j];

                    tg_xxxxxx_xxyzzz_0[j] = pb_x * tg_xxxxx_xxyzzz_0[j] + fr * tg_xxxxx_xxyzzz_1[j] + 2.5 * fl1_fx * (tg_xxxx_xxyzzz_0[j] - tg_xxxx_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxx_xyzzz_1[j];

                    tg_xxxxxx_xxzzzz_0[j] = pb_x * tg_xxxxx_xxzzzz_0[j] + fr * tg_xxxxx_xxzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxx_xxzzzz_0[j] - tg_xxxx_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxx_xzzzz_1[j];

                    tg_xxxxxx_xyyyyy_0[j] = pb_x * tg_xxxxx_xyyyyy_0[j] + fr * tg_xxxxx_xyyyyy_1[j] + 2.5 * fl1_fx * (tg_xxxx_xyyyyy_0[j] - tg_xxxx_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxx_yyyyy_1[j];

                    tg_xxxxxx_xyyyyz_0[j] = pb_x * tg_xxxxx_xyyyyz_0[j] + fr * tg_xxxxx_xyyyyz_1[j] + 2.5 * fl1_fx * (tg_xxxx_xyyyyz_0[j] - tg_xxxx_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxx_yyyyz_1[j];

                    tg_xxxxxx_xyyyzz_0[j] = pb_x * tg_xxxxx_xyyyzz_0[j] + fr * tg_xxxxx_xyyyzz_1[j] + 2.5 * fl1_fx * (tg_xxxx_xyyyzz_0[j] - tg_xxxx_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxx_yyyzz_1[j];

                    tg_xxxxxx_xyyzzz_0[j] = pb_x * tg_xxxxx_xyyzzz_0[j] + fr * tg_xxxxx_xyyzzz_1[j] + 2.5 * fl1_fx * (tg_xxxx_xyyzzz_0[j] - tg_xxxx_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxx_yyzzz_1[j];

                    tg_xxxxxx_xyzzzz_0[j] = pb_x * tg_xxxxx_xyzzzz_0[j] + fr * tg_xxxxx_xyzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxx_xyzzzz_0[j] - tg_xxxx_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxx_yzzzz_1[j];

                    tg_xxxxxx_xzzzzz_0[j] = pb_x * tg_xxxxx_xzzzzz_0[j] + fr * tg_xxxxx_xzzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxx_xzzzzz_0[j] - tg_xxxx_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxx_zzzzz_1[j];

                    tg_xxxxxx_yyyyyy_0[j] = pb_x * tg_xxxxx_yyyyyy_0[j] + fr * tg_xxxxx_yyyyyy_1[j] + 2.5 * fl1_fx * (tg_xxxx_yyyyyy_0[j] - tg_xxxx_yyyyyy_1[j] * fl1_fza);

                    tg_xxxxxx_yyyyyz_0[j] = pb_x * tg_xxxxx_yyyyyz_0[j] + fr * tg_xxxxx_yyyyyz_1[j] + 2.5 * fl1_fx * (tg_xxxx_yyyyyz_0[j] - tg_xxxx_yyyyyz_1[j] * fl1_fza);

                    tg_xxxxxx_yyyyzz_0[j] = pb_x * tg_xxxxx_yyyyzz_0[j] + fr * tg_xxxxx_yyyyzz_1[j] + 2.5 * fl1_fx * (tg_xxxx_yyyyzz_0[j] - tg_xxxx_yyyyzz_1[j] * fl1_fza);

                    tg_xxxxxx_yyyzzz_0[j] = pb_x * tg_xxxxx_yyyzzz_0[j] + fr * tg_xxxxx_yyyzzz_1[j] + 2.5 * fl1_fx * (tg_xxxx_yyyzzz_0[j] - tg_xxxx_yyyzzz_1[j] * fl1_fza);

                    tg_xxxxxx_yyzzzz_0[j] = pb_x * tg_xxxxx_yyzzzz_0[j] + fr * tg_xxxxx_yyzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxx_yyzzzz_0[j] - tg_xxxx_yyzzzz_1[j] * fl1_fza);

                    tg_xxxxxx_yzzzzz_0[j] = pb_x * tg_xxxxx_yzzzzz_0[j] + fr * tg_xxxxx_yzzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxx_yzzzzz_0[j] - tg_xxxx_yzzzzz_1[j] * fl1_fza);

                    tg_xxxxxx_zzzzzz_0[j] = pb_x * tg_xxxxx_zzzzzz_0[j] + fr * tg_xxxxx_zzzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxx_zzzzzz_0[j] - tg_xxxx_zzzzzz_1[j] * fl1_fza);

                    tg_xxxxxy_xxxxxx_0[j] = pb_x * tg_xxxxy_xxxxxx_0[j] + fr * tg_xxxxy_xxxxxx_1[j] + 2.0 * fl1_fx * (tg_xxxy_xxxxxx_0[j] - tg_xxxy_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxxxy_xxxxx_1[j];

                    tg_xxxxxy_xxxxxy_0[j] = pb_x * tg_xxxxy_xxxxxy_0[j] + fr * tg_xxxxy_xxxxxy_1[j] + 2.0 * fl1_fx * (tg_xxxy_xxxxxy_0[j] - tg_xxxy_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxy_xxxxy_1[j];

                    tg_xxxxxy_xxxxxz_0[j] = pb_x * tg_xxxxy_xxxxxz_0[j] + fr * tg_xxxxy_xxxxxz_1[j] + 2.0 * fl1_fx * (tg_xxxy_xxxxxz_0[j] - tg_xxxy_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxy_xxxxz_1[j];

                    tg_xxxxxy_xxxxyy_0[j] = pb_x * tg_xxxxy_xxxxyy_0[j] + fr * tg_xxxxy_xxxxyy_1[j] + 2.0 * fl1_fx * (tg_xxxy_xxxxyy_0[j] - tg_xxxy_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxy_xxxyy_1[j];

                    tg_xxxxxy_xxxxyz_0[j] = pb_x * tg_xxxxy_xxxxyz_0[j] + fr * tg_xxxxy_xxxxyz_1[j] + 2.0 * fl1_fx * (tg_xxxy_xxxxyz_0[j] - tg_xxxy_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxy_xxxyz_1[j];

                    tg_xxxxxy_xxxxzz_0[j] = pb_x * tg_xxxxy_xxxxzz_0[j] + fr * tg_xxxxy_xxxxzz_1[j] + 2.0 * fl1_fx * (tg_xxxy_xxxxzz_0[j] - tg_xxxy_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxy_xxxzz_1[j];

                    tg_xxxxxy_xxxyyy_0[j] = pb_x * tg_xxxxy_xxxyyy_0[j] + fr * tg_xxxxy_xxxyyy_1[j] + 2.0 * fl1_fx * (tg_xxxy_xxxyyy_0[j] - tg_xxxy_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxy_xxyyy_1[j];

                    tg_xxxxxy_xxxyyz_0[j] = pb_x * tg_xxxxy_xxxyyz_0[j] + fr * tg_xxxxy_xxxyyz_1[j] + 2.0 * fl1_fx * (tg_xxxy_xxxyyz_0[j] - tg_xxxy_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxy_xxyyz_1[j];

                    tg_xxxxxy_xxxyzz_0[j] = pb_x * tg_xxxxy_xxxyzz_0[j] + fr * tg_xxxxy_xxxyzz_1[j] + 2.0 * fl1_fx * (tg_xxxy_xxxyzz_0[j] - tg_xxxy_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxy_xxyzz_1[j];

                    tg_xxxxxy_xxxzzz_0[j] = pb_x * tg_xxxxy_xxxzzz_0[j] + fr * tg_xxxxy_xxxzzz_1[j] + 2.0 * fl1_fx * (tg_xxxy_xxxzzz_0[j] - tg_xxxy_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxy_xxzzz_1[j];

                    tg_xxxxxy_xxyyyy_0[j] = pb_x * tg_xxxxy_xxyyyy_0[j] + fr * tg_xxxxy_xxyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxy_xxyyyy_0[j] - tg_xxxy_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxy_xyyyy_1[j];

                    tg_xxxxxy_xxyyyz_0[j] = pb_x * tg_xxxxy_xxyyyz_0[j] + fr * tg_xxxxy_xxyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxy_xxyyyz_0[j] - tg_xxxy_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxy_xyyyz_1[j];

                    tg_xxxxxy_xxyyzz_0[j] = pb_x * tg_xxxxy_xxyyzz_0[j] + fr * tg_xxxxy_xxyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxy_xxyyzz_0[j] - tg_xxxy_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxy_xyyzz_1[j];

                    tg_xxxxxy_xxyzzz_0[j] = pb_x * tg_xxxxy_xxyzzz_0[j] + fr * tg_xxxxy_xxyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxy_xxyzzz_0[j] - tg_xxxy_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxy_xyzzz_1[j];

                    tg_xxxxxy_xxzzzz_0[j] = pb_x * tg_xxxxy_xxzzzz_0[j] + fr * tg_xxxxy_xxzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxy_xxzzzz_0[j] - tg_xxxy_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxy_xzzzz_1[j];

                    tg_xxxxxy_xyyyyy_0[j] = pb_x * tg_xxxxy_xyyyyy_0[j] + fr * tg_xxxxy_xyyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxy_xyyyyy_0[j] - tg_xxxy_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxy_yyyyy_1[j];

                    tg_xxxxxy_xyyyyz_0[j] = pb_x * tg_xxxxy_xyyyyz_0[j] + fr * tg_xxxxy_xyyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxy_xyyyyz_0[j] - tg_xxxy_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxy_yyyyz_1[j];

                    tg_xxxxxy_xyyyzz_0[j] = pb_x * tg_xxxxy_xyyyzz_0[j] + fr * tg_xxxxy_xyyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxy_xyyyzz_0[j] - tg_xxxy_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxy_yyyzz_1[j];

                    tg_xxxxxy_xyyzzz_0[j] = pb_x * tg_xxxxy_xyyzzz_0[j] + fr * tg_xxxxy_xyyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxy_xyyzzz_0[j] - tg_xxxy_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxy_yyzzz_1[j];

                    tg_xxxxxy_xyzzzz_0[j] = pb_x * tg_xxxxy_xyzzzz_0[j] + fr * tg_xxxxy_xyzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxy_xyzzzz_0[j] - tg_xxxy_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxy_yzzzz_1[j];

                    tg_xxxxxy_xzzzzz_0[j] = pb_x * tg_xxxxy_xzzzzz_0[j] + fr * tg_xxxxy_xzzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxy_xzzzzz_0[j] - tg_xxxy_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxy_zzzzz_1[j];

                    tg_xxxxxy_yyyyyy_0[j] = pb_x * tg_xxxxy_yyyyyy_0[j] + fr * tg_xxxxy_yyyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxy_yyyyyy_0[j] - tg_xxxy_yyyyyy_1[j] * fl1_fza);

                    tg_xxxxxy_yyyyyz_0[j] = pb_x * tg_xxxxy_yyyyyz_0[j] + fr * tg_xxxxy_yyyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxy_yyyyyz_0[j] - tg_xxxy_yyyyyz_1[j] * fl1_fza);

                    tg_xxxxxy_yyyyzz_0[j] = pb_x * tg_xxxxy_yyyyzz_0[j] + fr * tg_xxxxy_yyyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxy_yyyyzz_0[j] - tg_xxxy_yyyyzz_1[j] * fl1_fza);

                    tg_xxxxxy_yyyzzz_0[j] = pb_x * tg_xxxxy_yyyzzz_0[j] + fr * tg_xxxxy_yyyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxy_yyyzzz_0[j] - tg_xxxy_yyyzzz_1[j] * fl1_fza);

                    tg_xxxxxy_yyzzzz_0[j] = pb_x * tg_xxxxy_yyzzzz_0[j] + fr * tg_xxxxy_yyzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxy_yyzzzz_0[j] - tg_xxxy_yyzzzz_1[j] * fl1_fza);

                    tg_xxxxxy_yzzzzz_0[j] = pb_x * tg_xxxxy_yzzzzz_0[j] + fr * tg_xxxxy_yzzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxy_yzzzzz_0[j] - tg_xxxy_yzzzzz_1[j] * fl1_fza);

                    tg_xxxxxy_zzzzzz_0[j] = pb_x * tg_xxxxy_zzzzzz_0[j] + fr * tg_xxxxy_zzzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxy_zzzzzz_0[j] - tg_xxxy_zzzzzz_1[j] * fl1_fza);

                    tg_xxxxxz_xxxxxx_0[j] = pb_x * tg_xxxxz_xxxxxx_0[j] + fr * tg_xxxxz_xxxxxx_1[j] + 2.0 * fl1_fx * (tg_xxxz_xxxxxx_0[j] - tg_xxxz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxxxz_xxxxx_1[j];

                    tg_xxxxxz_xxxxxy_0[j] = pb_x * tg_xxxxz_xxxxxy_0[j] + fr * tg_xxxxz_xxxxxy_1[j] + 2.0 * fl1_fx * (tg_xxxz_xxxxxy_0[j] - tg_xxxz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxz_xxxxy_1[j];

                    tg_xxxxxz_xxxxxz_0[j] = pb_x * tg_xxxxz_xxxxxz_0[j] + fr * tg_xxxxz_xxxxxz_1[j] + 2.0 * fl1_fx * (tg_xxxz_xxxxxz_0[j] - tg_xxxz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxz_xxxxz_1[j];

                    tg_xxxxxz_xxxxyy_0[j] = pb_x * tg_xxxxz_xxxxyy_0[j] + fr * tg_xxxxz_xxxxyy_1[j] + 2.0 * fl1_fx * (tg_xxxz_xxxxyy_0[j] - tg_xxxz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxz_xxxyy_1[j];

                    tg_xxxxxz_xxxxyz_0[j] = pb_x * tg_xxxxz_xxxxyz_0[j] + fr * tg_xxxxz_xxxxyz_1[j] + 2.0 * fl1_fx * (tg_xxxz_xxxxyz_0[j] - tg_xxxz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxz_xxxyz_1[j];

                    tg_xxxxxz_xxxxzz_0[j] = pb_x * tg_xxxxz_xxxxzz_0[j] + fr * tg_xxxxz_xxxxzz_1[j] + 2.0 * fl1_fx * (tg_xxxz_xxxxzz_0[j] - tg_xxxz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxz_xxxzz_1[j];

                    tg_xxxxxz_xxxyyy_0[j] = pb_x * tg_xxxxz_xxxyyy_0[j] + fr * tg_xxxxz_xxxyyy_1[j] + 2.0 * fl1_fx * (tg_xxxz_xxxyyy_0[j] - tg_xxxz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxz_xxyyy_1[j];

                    tg_xxxxxz_xxxyyz_0[j] = pb_x * tg_xxxxz_xxxyyz_0[j] + fr * tg_xxxxz_xxxyyz_1[j] + 2.0 * fl1_fx * (tg_xxxz_xxxyyz_0[j] - tg_xxxz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxz_xxyyz_1[j];

                    tg_xxxxxz_xxxyzz_0[j] = pb_x * tg_xxxxz_xxxyzz_0[j] + fr * tg_xxxxz_xxxyzz_1[j] + 2.0 * fl1_fx * (tg_xxxz_xxxyzz_0[j] - tg_xxxz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxz_xxyzz_1[j];

                    tg_xxxxxz_xxxzzz_0[j] = pb_x * tg_xxxxz_xxxzzz_0[j] + fr * tg_xxxxz_xxxzzz_1[j] + 2.0 * fl1_fx * (tg_xxxz_xxxzzz_0[j] - tg_xxxz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxz_xxzzz_1[j];

                    tg_xxxxxz_xxyyyy_0[j] = pb_x * tg_xxxxz_xxyyyy_0[j] + fr * tg_xxxxz_xxyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxz_xxyyyy_0[j] - tg_xxxz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxz_xyyyy_1[j];

                    tg_xxxxxz_xxyyyz_0[j] = pb_x * tg_xxxxz_xxyyyz_0[j] + fr * tg_xxxxz_xxyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxz_xxyyyz_0[j] - tg_xxxz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxz_xyyyz_1[j];

                    tg_xxxxxz_xxyyzz_0[j] = pb_x * tg_xxxxz_xxyyzz_0[j] + fr * tg_xxxxz_xxyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxz_xxyyzz_0[j] - tg_xxxz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxz_xyyzz_1[j];

                    tg_xxxxxz_xxyzzz_0[j] = pb_x * tg_xxxxz_xxyzzz_0[j] + fr * tg_xxxxz_xxyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxz_xxyzzz_0[j] - tg_xxxz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxz_xyzzz_1[j];

                    tg_xxxxxz_xxzzzz_0[j] = pb_x * tg_xxxxz_xxzzzz_0[j] + fr * tg_xxxxz_xxzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxz_xxzzzz_0[j] - tg_xxxz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxz_xzzzz_1[j];

                    tg_xxxxxz_xyyyyy_0[j] = pb_x * tg_xxxxz_xyyyyy_0[j] + fr * tg_xxxxz_xyyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxz_xyyyyy_0[j] - tg_xxxz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxz_yyyyy_1[j];

                    tg_xxxxxz_xyyyyz_0[j] = pb_x * tg_xxxxz_xyyyyz_0[j] + fr * tg_xxxxz_xyyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxz_xyyyyz_0[j] - tg_xxxz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxz_yyyyz_1[j];

                    tg_xxxxxz_xyyyzz_0[j] = pb_x * tg_xxxxz_xyyyzz_0[j] + fr * tg_xxxxz_xyyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxz_xyyyzz_0[j] - tg_xxxz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxz_yyyzz_1[j];

                    tg_xxxxxz_xyyzzz_0[j] = pb_x * tg_xxxxz_xyyzzz_0[j] + fr * tg_xxxxz_xyyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxz_xyyzzz_0[j] - tg_xxxz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxz_yyzzz_1[j];

                    tg_xxxxxz_xyzzzz_0[j] = pb_x * tg_xxxxz_xyzzzz_0[j] + fr * tg_xxxxz_xyzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxz_xyzzzz_0[j] - tg_xxxz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxz_yzzzz_1[j];

                    tg_xxxxxz_xzzzzz_0[j] = pb_x * tg_xxxxz_xzzzzz_0[j] + fr * tg_xxxxz_xzzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxz_xzzzzz_0[j] - tg_xxxz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxz_zzzzz_1[j];

                    tg_xxxxxz_yyyyyy_0[j] = pb_x * tg_xxxxz_yyyyyy_0[j] + fr * tg_xxxxz_yyyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxz_yyyyyy_0[j] - tg_xxxz_yyyyyy_1[j] * fl1_fza);

                    tg_xxxxxz_yyyyyz_0[j] = pb_x * tg_xxxxz_yyyyyz_0[j] + fr * tg_xxxxz_yyyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxz_yyyyyz_0[j] - tg_xxxz_yyyyyz_1[j] * fl1_fza);

                    tg_xxxxxz_yyyyzz_0[j] = pb_x * tg_xxxxz_yyyyzz_0[j] + fr * tg_xxxxz_yyyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxz_yyyyzz_0[j] - tg_xxxz_yyyyzz_1[j] * fl1_fza);

                    tg_xxxxxz_yyyzzz_0[j] = pb_x * tg_xxxxz_yyyzzz_0[j] + fr * tg_xxxxz_yyyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxz_yyyzzz_0[j] - tg_xxxz_yyyzzz_1[j] * fl1_fza);

                    tg_xxxxxz_yyzzzz_0[j] = pb_x * tg_xxxxz_yyzzzz_0[j] + fr * tg_xxxxz_yyzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxz_yyzzzz_0[j] - tg_xxxz_yyzzzz_1[j] * fl1_fza);

                    tg_xxxxxz_yzzzzz_0[j] = pb_x * tg_xxxxz_yzzzzz_0[j] + fr * tg_xxxxz_yzzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxz_yzzzzz_0[j] - tg_xxxz_yzzzzz_1[j] * fl1_fza);

                    tg_xxxxxz_zzzzzz_0[j] = pb_x * tg_xxxxz_zzzzzz_0[j] + fr * tg_xxxxz_zzzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxz_zzzzzz_0[j] - tg_xxxz_zzzzzz_1[j] * fl1_fza);

                    tg_xxxxyy_xxxxxx_0[j] = pb_x * tg_xxxyy_xxxxxx_0[j] + fr * tg_xxxyy_xxxxxx_1[j] + 1.5 * fl1_fx * (tg_xxyy_xxxxxx_0[j] - tg_xxyy_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxxyy_xxxxx_1[j];

                    tg_xxxxyy_xxxxxy_0[j] = pb_x * tg_xxxyy_xxxxxy_0[j] + fr * tg_xxxyy_xxxxxy_1[j] + 1.5 * fl1_fx * (tg_xxyy_xxxxxy_0[j] - tg_xxyy_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxyy_xxxxy_1[j];

                    tg_xxxxyy_xxxxxz_0[j] = pb_x * tg_xxxyy_xxxxxz_0[j] + fr * tg_xxxyy_xxxxxz_1[j] + 1.5 * fl1_fx * (tg_xxyy_xxxxxz_0[j] - tg_xxyy_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxyy_xxxxz_1[j];

                    tg_xxxxyy_xxxxyy_0[j] = pb_x * tg_xxxyy_xxxxyy_0[j] + fr * tg_xxxyy_xxxxyy_1[j] + 1.5 * fl1_fx * (tg_xxyy_xxxxyy_0[j] - tg_xxyy_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyy_xxxyy_1[j];

                    tg_xxxxyy_xxxxyz_0[j] = pb_x * tg_xxxyy_xxxxyz_0[j] + fr * tg_xxxyy_xxxxyz_1[j] + 1.5 * fl1_fx * (tg_xxyy_xxxxyz_0[j] - tg_xxyy_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyy_xxxyz_1[j];

                    tg_xxxxyy_xxxxzz_0[j] = pb_x * tg_xxxyy_xxxxzz_0[j] + fr * tg_xxxyy_xxxxzz_1[j] + 1.5 * fl1_fx * (tg_xxyy_xxxxzz_0[j] - tg_xxyy_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyy_xxxzz_1[j];

                    tg_xxxxyy_xxxyyy_0[j] = pb_x * tg_xxxyy_xxxyyy_0[j] + fr * tg_xxxyy_xxxyyy_1[j] + 1.5 * fl1_fx * (tg_xxyy_xxxyyy_0[j] - tg_xxyy_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyy_xxyyy_1[j];

                    tg_xxxxyy_xxxyyz_0[j] = pb_x * tg_xxxyy_xxxyyz_0[j] + fr * tg_xxxyy_xxxyyz_1[j] + 1.5 * fl1_fx * (tg_xxyy_xxxyyz_0[j] - tg_xxyy_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyy_xxyyz_1[j];

                    tg_xxxxyy_xxxyzz_0[j] = pb_x * tg_xxxyy_xxxyzz_0[j] + fr * tg_xxxyy_xxxyzz_1[j] + 1.5 * fl1_fx * (tg_xxyy_xxxyzz_0[j] - tg_xxyy_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyy_xxyzz_1[j];

                    tg_xxxxyy_xxxzzz_0[j] = pb_x * tg_xxxyy_xxxzzz_0[j] + fr * tg_xxxyy_xxxzzz_1[j] + 1.5 * fl1_fx * (tg_xxyy_xxxzzz_0[j] - tg_xxyy_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyy_xxzzz_1[j];

                    tg_xxxxyy_xxyyyy_0[j] = pb_x * tg_xxxyy_xxyyyy_0[j] + fr * tg_xxxyy_xxyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyy_xxyyyy_0[j] - tg_xxyy_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxyy_xyyyy_1[j];

                    tg_xxxxyy_xxyyyz_0[j] = pb_x * tg_xxxyy_xxyyyz_0[j] + fr * tg_xxxyy_xxyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyy_xxyyyz_0[j] - tg_xxyy_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyy_xyyyz_1[j];

                    tg_xxxxyy_xxyyzz_0[j] = pb_x * tg_xxxyy_xxyyzz_0[j] + fr * tg_xxxyy_xxyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyy_xxyyzz_0[j] - tg_xxyy_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyy_xyyzz_1[j];

                    tg_xxxxyy_xxyzzz_0[j] = pb_x * tg_xxxyy_xxyzzz_0[j] + fr * tg_xxxyy_xxyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyy_xxyzzz_0[j] - tg_xxyy_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyy_xyzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSISI_98_196(      CMemBlock2D<double>* primBuffer,
                                        const CRecursionMap&       recursionMap,
                                        const CMemBlock2D<double>& osFactors,
                                        const CMemBlock2D<double>& wpDistances,
                                        const CGtoPairsBlock&      braGtoPairsBlock,
                                        const CGtoPairsBlock&      ketGtoPairsBlock,
                                        const int32_t              nKetPrimPairs,
                                        const int32_t              iContrPair)
    {
        // Batch of Integrals (98,196)

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
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_6_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_6_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_5_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
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

                auto tg_xxxyy_xxzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 98); 

                auto tg_xxxyy_xyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 99); 

                auto tg_xxxyy_xyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 100); 

                auto tg_xxxyy_xyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 101); 

                auto tg_xxxyy_xyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 102); 

                auto tg_xxxyy_xyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 103); 

                auto tg_xxxyy_xzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 104); 

                auto tg_xxxyy_yyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 105); 

                auto tg_xxxyy_yyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 106); 

                auto tg_xxxyy_yyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 107); 

                auto tg_xxxyy_yyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 108); 

                auto tg_xxxyy_yyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 109); 

                auto tg_xxxyy_yzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 110); 

                auto tg_xxxyy_zzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 111); 

                auto tg_xxxyz_xxxxxx_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 112); 

                auto tg_xxxyz_xxxxxy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 113); 

                auto tg_xxxyz_xxxxxz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 114); 

                auto tg_xxxyz_xxxxyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 115); 

                auto tg_xxxyz_xxxxyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 116); 

                auto tg_xxxyz_xxxxzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 117); 

                auto tg_xxxyz_xxxyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 118); 

                auto tg_xxxyz_xxxyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 119); 

                auto tg_xxxyz_xxxyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 120); 

                auto tg_xxxyz_xxxzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 121); 

                auto tg_xxxyz_xxyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 122); 

                auto tg_xxxyz_xxyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 123); 

                auto tg_xxxyz_xxyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 124); 

                auto tg_xxxyz_xxyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 125); 

                auto tg_xxxyz_xxzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 126); 

                auto tg_xxxyz_xyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 127); 

                auto tg_xxxyz_xyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 128); 

                auto tg_xxxyz_xyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 129); 

                auto tg_xxxyz_xyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 130); 

                auto tg_xxxyz_xyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 131); 

                auto tg_xxxyz_xzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 132); 

                auto tg_xxxyz_yyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 133); 

                auto tg_xxxyz_yyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 134); 

                auto tg_xxxyz_yyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 135); 

                auto tg_xxxyz_yyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 136); 

                auto tg_xxxyz_yyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 137); 

                auto tg_xxxyz_yzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 138); 

                auto tg_xxxyz_zzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 139); 

                auto tg_xxxzz_xxxxxx_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 140); 

                auto tg_xxxzz_xxxxxy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 141); 

                auto tg_xxxzz_xxxxxz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 142); 

                auto tg_xxxzz_xxxxyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 143); 

                auto tg_xxxzz_xxxxyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 144); 

                auto tg_xxxzz_xxxxzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 145); 

                auto tg_xxxzz_xxxyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 146); 

                auto tg_xxxzz_xxxyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 147); 

                auto tg_xxxzz_xxxyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 148); 

                auto tg_xxxzz_xxxzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 149); 

                auto tg_xxxzz_xxyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 150); 

                auto tg_xxxzz_xxyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 151); 

                auto tg_xxxzz_xxyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 152); 

                auto tg_xxxzz_xxyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 153); 

                auto tg_xxxzz_xxzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 154); 

                auto tg_xxxzz_xyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 155); 

                auto tg_xxxzz_xyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 156); 

                auto tg_xxxzz_xyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 157); 

                auto tg_xxxzz_xyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 158); 

                auto tg_xxxzz_xyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 159); 

                auto tg_xxxzz_xzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 160); 

                auto tg_xxxzz_yyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 161); 

                auto tg_xxxzz_yyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 162); 

                auto tg_xxxzz_yyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 163); 

                auto tg_xxxzz_yyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 164); 

                auto tg_xxxzz_yyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 165); 

                auto tg_xxxzz_yzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 166); 

                auto tg_xxxzz_zzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 167); 

                auto tg_xxyyy_xxxxxx_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 168); 

                auto tg_xxyyy_xxxxxy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 169); 

                auto tg_xxyyy_xxxxxz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 170); 

                auto tg_xxyyy_xxxxyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 171); 

                auto tg_xxyyy_xxxxyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 172); 

                auto tg_xxyyy_xxxxzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 173); 

                auto tg_xxyyy_xxxyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 174); 

                auto tg_xxyyy_xxxyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 175); 

                auto tg_xxyyy_xxxyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 176); 

                auto tg_xxyyy_xxxzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 177); 

                auto tg_xxyyy_xxyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 178); 

                auto tg_xxyyy_xxyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 179); 

                auto tg_xxyyy_xxyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 180); 

                auto tg_xxyyy_xxyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 181); 

                auto tg_xxyyy_xxzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 182); 

                auto tg_xxyyy_xyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 183); 

                auto tg_xxyyy_xyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 184); 

                auto tg_xxyyy_xyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 185); 

                auto tg_xxyyy_xyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 186); 

                auto tg_xxyyy_xyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 187); 

                auto tg_xxyyy_xzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 188); 

                auto tg_xxyyy_yyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 189); 

                auto tg_xxyyy_yyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 190); 

                auto tg_xxyyy_yyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 191); 

                auto tg_xxyyy_yyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 192); 

                auto tg_xxyyy_yyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 193); 

                auto tg_xxyyy_yzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 194); 

                auto tg_xxyyy_zzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 195); 

                auto tg_xxxyy_xxzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 98); 

                auto tg_xxxyy_xyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 99); 

                auto tg_xxxyy_xyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 100); 

                auto tg_xxxyy_xyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 101); 

                auto tg_xxxyy_xyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 102); 

                auto tg_xxxyy_xyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 103); 

                auto tg_xxxyy_xzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 104); 

                auto tg_xxxyy_yyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 105); 

                auto tg_xxxyy_yyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 106); 

                auto tg_xxxyy_yyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 107); 

                auto tg_xxxyy_yyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 108); 

                auto tg_xxxyy_yyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 109); 

                auto tg_xxxyy_yzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 110); 

                auto tg_xxxyy_zzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 111); 

                auto tg_xxxyz_xxxxxx_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 112); 

                auto tg_xxxyz_xxxxxy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 113); 

                auto tg_xxxyz_xxxxxz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 114); 

                auto tg_xxxyz_xxxxyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 115); 

                auto tg_xxxyz_xxxxyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 116); 

                auto tg_xxxyz_xxxxzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 117); 

                auto tg_xxxyz_xxxyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 118); 

                auto tg_xxxyz_xxxyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 119); 

                auto tg_xxxyz_xxxyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 120); 

                auto tg_xxxyz_xxxzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 121); 

                auto tg_xxxyz_xxyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 122); 

                auto tg_xxxyz_xxyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 123); 

                auto tg_xxxyz_xxyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 124); 

                auto tg_xxxyz_xxyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 125); 

                auto tg_xxxyz_xxzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 126); 

                auto tg_xxxyz_xyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 127); 

                auto tg_xxxyz_xyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 128); 

                auto tg_xxxyz_xyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 129); 

                auto tg_xxxyz_xyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 130); 

                auto tg_xxxyz_xyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 131); 

                auto tg_xxxyz_xzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 132); 

                auto tg_xxxyz_yyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 133); 

                auto tg_xxxyz_yyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 134); 

                auto tg_xxxyz_yyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 135); 

                auto tg_xxxyz_yyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 136); 

                auto tg_xxxyz_yyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 137); 

                auto tg_xxxyz_yzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 138); 

                auto tg_xxxyz_zzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 139); 

                auto tg_xxxzz_xxxxxx_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 140); 

                auto tg_xxxzz_xxxxxy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 141); 

                auto tg_xxxzz_xxxxxz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 142); 

                auto tg_xxxzz_xxxxyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 143); 

                auto tg_xxxzz_xxxxyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 144); 

                auto tg_xxxzz_xxxxzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 145); 

                auto tg_xxxzz_xxxyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 146); 

                auto tg_xxxzz_xxxyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 147); 

                auto tg_xxxzz_xxxyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 148); 

                auto tg_xxxzz_xxxzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 149); 

                auto tg_xxxzz_xxyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 150); 

                auto tg_xxxzz_xxyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 151); 

                auto tg_xxxzz_xxyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 152); 

                auto tg_xxxzz_xxyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 153); 

                auto tg_xxxzz_xxzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 154); 

                auto tg_xxxzz_xyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 155); 

                auto tg_xxxzz_xyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 156); 

                auto tg_xxxzz_xyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 157); 

                auto tg_xxxzz_xyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 158); 

                auto tg_xxxzz_xyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 159); 

                auto tg_xxxzz_xzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 160); 

                auto tg_xxxzz_yyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 161); 

                auto tg_xxxzz_yyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 162); 

                auto tg_xxxzz_yyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 163); 

                auto tg_xxxzz_yyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 164); 

                auto tg_xxxzz_yyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 165); 

                auto tg_xxxzz_yzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 166); 

                auto tg_xxxzz_zzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 167); 

                auto tg_xxyyy_xxxxxx_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 168); 

                auto tg_xxyyy_xxxxxy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 169); 

                auto tg_xxyyy_xxxxxz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 170); 

                auto tg_xxyyy_xxxxyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 171); 

                auto tg_xxyyy_xxxxyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 172); 

                auto tg_xxyyy_xxxxzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 173); 

                auto tg_xxyyy_xxxyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 174); 

                auto tg_xxyyy_xxxyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 175); 

                auto tg_xxyyy_xxxyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 176); 

                auto tg_xxyyy_xxxzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 177); 

                auto tg_xxyyy_xxyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 178); 

                auto tg_xxyyy_xxyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 179); 

                auto tg_xxyyy_xxyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 180); 

                auto tg_xxyyy_xxyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 181); 

                auto tg_xxyyy_xxzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 182); 

                auto tg_xxyyy_xyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 183); 

                auto tg_xxyyy_xyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 184); 

                auto tg_xxyyy_xyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 185); 

                auto tg_xxyyy_xyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 186); 

                auto tg_xxyyy_xyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 187); 

                auto tg_xxyyy_xzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 188); 

                auto tg_xxyyy_yyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 189); 

                auto tg_xxyyy_yyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 190); 

                auto tg_xxyyy_yyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 191); 

                auto tg_xxyyy_yyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 192); 

                auto tg_xxyyy_yyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 193); 

                auto tg_xxyyy_yzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 194); 

                auto tg_xxyyy_zzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 195); 

                auto tg_xxyy_xxzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 98); 

                auto tg_xxyy_xyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 99); 

                auto tg_xxyy_xyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 100); 

                auto tg_xxyy_xyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 101); 

                auto tg_xxyy_xyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 102); 

                auto tg_xxyy_xyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 103); 

                auto tg_xxyy_xzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 104); 

                auto tg_xxyy_yyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 105); 

                auto tg_xxyy_yyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 106); 

                auto tg_xxyy_yyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 107); 

                auto tg_xxyy_yyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 108); 

                auto tg_xxyy_yyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 109); 

                auto tg_xxyy_yzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 110); 

                auto tg_xxyy_zzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 111); 

                auto tg_xxyz_xxxxxx_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 112); 

                auto tg_xxyz_xxxxxy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 113); 

                auto tg_xxyz_xxxxxz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 114); 

                auto tg_xxyz_xxxxyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 115); 

                auto tg_xxyz_xxxxyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 116); 

                auto tg_xxyz_xxxxzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 117); 

                auto tg_xxyz_xxxyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 118); 

                auto tg_xxyz_xxxyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 119); 

                auto tg_xxyz_xxxyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 120); 

                auto tg_xxyz_xxxzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 121); 

                auto tg_xxyz_xxyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 122); 

                auto tg_xxyz_xxyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 123); 

                auto tg_xxyz_xxyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 124); 

                auto tg_xxyz_xxyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 125); 

                auto tg_xxyz_xxzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 126); 

                auto tg_xxyz_xyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 127); 

                auto tg_xxyz_xyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 128); 

                auto tg_xxyz_xyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 129); 

                auto tg_xxyz_xyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 130); 

                auto tg_xxyz_xyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 131); 

                auto tg_xxyz_xzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 132); 

                auto tg_xxyz_yyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 133); 

                auto tg_xxyz_yyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 134); 

                auto tg_xxyz_yyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 135); 

                auto tg_xxyz_yyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 136); 

                auto tg_xxyz_yyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 137); 

                auto tg_xxyz_yzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 138); 

                auto tg_xxyz_zzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 139); 

                auto tg_xxzz_xxxxxx_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 140); 

                auto tg_xxzz_xxxxxy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 141); 

                auto tg_xxzz_xxxxxz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 142); 

                auto tg_xxzz_xxxxyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 143); 

                auto tg_xxzz_xxxxyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 144); 

                auto tg_xxzz_xxxxzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 145); 

                auto tg_xxzz_xxxyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 146); 

                auto tg_xxzz_xxxyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 147); 

                auto tg_xxzz_xxxyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 148); 

                auto tg_xxzz_xxxzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 149); 

                auto tg_xxzz_xxyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 150); 

                auto tg_xxzz_xxyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 151); 

                auto tg_xxzz_xxyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 152); 

                auto tg_xxzz_xxyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 153); 

                auto tg_xxzz_xxzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 154); 

                auto tg_xxzz_xyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 155); 

                auto tg_xxzz_xyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 156); 

                auto tg_xxzz_xyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 157); 

                auto tg_xxzz_xyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 158); 

                auto tg_xxzz_xyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 159); 

                auto tg_xxzz_xzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 160); 

                auto tg_xxzz_yyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 161); 

                auto tg_xxzz_yyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 162); 

                auto tg_xxzz_yyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 163); 

                auto tg_xxzz_yyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 164); 

                auto tg_xxzz_yyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 165); 

                auto tg_xxzz_yzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 166); 

                auto tg_xxzz_zzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 167); 

                auto tg_xyyy_xxxxxx_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 168); 

                auto tg_xyyy_xxxxxy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 169); 

                auto tg_xyyy_xxxxxz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 170); 

                auto tg_xyyy_xxxxyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 171); 

                auto tg_xyyy_xxxxyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 172); 

                auto tg_xyyy_xxxxzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 173); 

                auto tg_xyyy_xxxyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 174); 

                auto tg_xyyy_xxxyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 175); 

                auto tg_xyyy_xxxyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 176); 

                auto tg_xyyy_xxxzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 177); 

                auto tg_xyyy_xxyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 178); 

                auto tg_xyyy_xxyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 179); 

                auto tg_xyyy_xxyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 180); 

                auto tg_xyyy_xxyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 181); 

                auto tg_xyyy_xxzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 182); 

                auto tg_xyyy_xyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 183); 

                auto tg_xyyy_xyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 184); 

                auto tg_xyyy_xyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 185); 

                auto tg_xyyy_xyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 186); 

                auto tg_xyyy_xyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 187); 

                auto tg_xyyy_xzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 188); 

                auto tg_xyyy_yyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 189); 

                auto tg_xyyy_yyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 190); 

                auto tg_xyyy_yyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 191); 

                auto tg_xyyy_yyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 192); 

                auto tg_xyyy_yyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 193); 

                auto tg_xyyy_yzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 194); 

                auto tg_xyyy_zzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 195); 

                auto tg_xxyy_xxzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 98); 

                auto tg_xxyy_xyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 99); 

                auto tg_xxyy_xyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 100); 

                auto tg_xxyy_xyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 101); 

                auto tg_xxyy_xyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 102); 

                auto tg_xxyy_xyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 103); 

                auto tg_xxyy_xzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 104); 

                auto tg_xxyy_yyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 105); 

                auto tg_xxyy_yyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 106); 

                auto tg_xxyy_yyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 107); 

                auto tg_xxyy_yyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 108); 

                auto tg_xxyy_yyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 109); 

                auto tg_xxyy_yzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 110); 

                auto tg_xxyy_zzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 111); 

                auto tg_xxyz_xxxxxx_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 112); 

                auto tg_xxyz_xxxxxy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 113); 

                auto tg_xxyz_xxxxxz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 114); 

                auto tg_xxyz_xxxxyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 115); 

                auto tg_xxyz_xxxxyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 116); 

                auto tg_xxyz_xxxxzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 117); 

                auto tg_xxyz_xxxyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 118); 

                auto tg_xxyz_xxxyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 119); 

                auto tg_xxyz_xxxyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 120); 

                auto tg_xxyz_xxxzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 121); 

                auto tg_xxyz_xxyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 122); 

                auto tg_xxyz_xxyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 123); 

                auto tg_xxyz_xxyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 124); 

                auto tg_xxyz_xxyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 125); 

                auto tg_xxyz_xxzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 126); 

                auto tg_xxyz_xyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 127); 

                auto tg_xxyz_xyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 128); 

                auto tg_xxyz_xyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 129); 

                auto tg_xxyz_xyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 130); 

                auto tg_xxyz_xyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 131); 

                auto tg_xxyz_xzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 132); 

                auto tg_xxyz_yyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 133); 

                auto tg_xxyz_yyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 134); 

                auto tg_xxyz_yyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 135); 

                auto tg_xxyz_yyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 136); 

                auto tg_xxyz_yyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 137); 

                auto tg_xxyz_yzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 138); 

                auto tg_xxyz_zzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 139); 

                auto tg_xxzz_xxxxxx_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 140); 

                auto tg_xxzz_xxxxxy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 141); 

                auto tg_xxzz_xxxxxz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 142); 

                auto tg_xxzz_xxxxyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 143); 

                auto tg_xxzz_xxxxyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 144); 

                auto tg_xxzz_xxxxzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 145); 

                auto tg_xxzz_xxxyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 146); 

                auto tg_xxzz_xxxyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 147); 

                auto tg_xxzz_xxxyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 148); 

                auto tg_xxzz_xxxzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 149); 

                auto tg_xxzz_xxyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 150); 

                auto tg_xxzz_xxyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 151); 

                auto tg_xxzz_xxyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 152); 

                auto tg_xxzz_xxyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 153); 

                auto tg_xxzz_xxzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 154); 

                auto tg_xxzz_xyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 155); 

                auto tg_xxzz_xyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 156); 

                auto tg_xxzz_xyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 157); 

                auto tg_xxzz_xyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 158); 

                auto tg_xxzz_xyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 159); 

                auto tg_xxzz_xzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 160); 

                auto tg_xxzz_yyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 161); 

                auto tg_xxzz_yyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 162); 

                auto tg_xxzz_yyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 163); 

                auto tg_xxzz_yyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 164); 

                auto tg_xxzz_yyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 165); 

                auto tg_xxzz_yzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 166); 

                auto tg_xxzz_zzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 167); 

                auto tg_xyyy_xxxxxx_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 168); 

                auto tg_xyyy_xxxxxy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 169); 

                auto tg_xyyy_xxxxxz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 170); 

                auto tg_xyyy_xxxxyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 171); 

                auto tg_xyyy_xxxxyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 172); 

                auto tg_xyyy_xxxxzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 173); 

                auto tg_xyyy_xxxyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 174); 

                auto tg_xyyy_xxxyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 175); 

                auto tg_xyyy_xxxyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 176); 

                auto tg_xyyy_xxxzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 177); 

                auto tg_xyyy_xxyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 178); 

                auto tg_xyyy_xxyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 179); 

                auto tg_xyyy_xxyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 180); 

                auto tg_xyyy_xxyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 181); 

                auto tg_xyyy_xxzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 182); 

                auto tg_xyyy_xyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 183); 

                auto tg_xyyy_xyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 184); 

                auto tg_xyyy_xyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 185); 

                auto tg_xyyy_xyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 186); 

                auto tg_xyyy_xyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 187); 

                auto tg_xyyy_xzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 188); 

                auto tg_xyyy_yyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 189); 

                auto tg_xyyy_yyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 190); 

                auto tg_xyyy_yyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 191); 

                auto tg_xyyy_yyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 192); 

                auto tg_xyyy_yyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 193); 

                auto tg_xyyy_yzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 194); 

                auto tg_xyyy_zzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 195); 

                auto tg_xxxyy_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 77); 

                auto tg_xxxyy_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 78); 

                auto tg_xxxyy_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 79); 

                auto tg_xxxyy_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 80); 

                auto tg_xxxyy_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 81); 

                auto tg_xxxyy_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 82); 

                auto tg_xxxyy_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 83); 

                auto tg_xxxyz_xxxxx_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 84); 

                auto tg_xxxyz_xxxxy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 85); 

                auto tg_xxxyz_xxxxz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 86); 

                auto tg_xxxyz_xxxyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 87); 

                auto tg_xxxyz_xxxyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 88); 

                auto tg_xxxyz_xxxzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 89); 

                auto tg_xxxyz_xxyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 90); 

                auto tg_xxxyz_xxyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 91); 

                auto tg_xxxyz_xxyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 92); 

                auto tg_xxxyz_xxzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 93); 

                auto tg_xxxyz_xyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 94); 

                auto tg_xxxyz_xyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 95); 

                auto tg_xxxyz_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 96); 

                auto tg_xxxyz_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 97); 

                auto tg_xxxyz_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 98); 

                auto tg_xxxyz_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 99); 

                auto tg_xxxyz_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 100); 

                auto tg_xxxyz_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 101); 

                auto tg_xxxyz_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 102); 

                auto tg_xxxyz_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 103); 

                auto tg_xxxyz_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 104); 

                auto tg_xxxzz_xxxxx_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 105); 

                auto tg_xxxzz_xxxxy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 106); 

                auto tg_xxxzz_xxxxz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 107); 

                auto tg_xxxzz_xxxyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 108); 

                auto tg_xxxzz_xxxyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 109); 

                auto tg_xxxzz_xxxzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 110); 

                auto tg_xxxzz_xxyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 111); 

                auto tg_xxxzz_xxyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 112); 

                auto tg_xxxzz_xxyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 113); 

                auto tg_xxxzz_xxzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 114); 

                auto tg_xxxzz_xyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 115); 

                auto tg_xxxzz_xyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 116); 

                auto tg_xxxzz_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 117); 

                auto tg_xxxzz_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 118); 

                auto tg_xxxzz_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 119); 

                auto tg_xxxzz_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 120); 

                auto tg_xxxzz_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 121); 

                auto tg_xxxzz_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 122); 

                auto tg_xxxzz_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 123); 

                auto tg_xxxzz_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 124); 

                auto tg_xxxzz_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 125); 

                auto tg_xxyyy_xxxxx_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 126); 

                auto tg_xxyyy_xxxxy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 127); 

                auto tg_xxyyy_xxxxz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 128); 

                auto tg_xxyyy_xxxyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 129); 

                auto tg_xxyyy_xxxyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 130); 

                auto tg_xxyyy_xxxzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 131); 

                auto tg_xxyyy_xxyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 132); 

                auto tg_xxyyy_xxyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 133); 

                auto tg_xxyyy_xxyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 134); 

                auto tg_xxyyy_xxzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 135); 

                auto tg_xxyyy_xyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 136); 

                auto tg_xxyyy_xyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 137); 

                auto tg_xxyyy_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 138); 

                auto tg_xxyyy_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 139); 

                auto tg_xxyyy_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 140); 

                auto tg_xxyyy_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 141); 

                auto tg_xxyyy_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 142); 

                auto tg_xxyyy_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 143); 

                auto tg_xxyyy_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 144); 

                auto tg_xxyyy_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 145); 

                auto tg_xxyyy_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 146); 

                // set up pointers to integrals

                auto tg_xxxxyy_xxzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 98); 

                auto tg_xxxxyy_xyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 99); 

                auto tg_xxxxyy_xyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 100); 

                auto tg_xxxxyy_xyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 101); 

                auto tg_xxxxyy_xyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 102); 

                auto tg_xxxxyy_xyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 103); 

                auto tg_xxxxyy_xzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 104); 

                auto tg_xxxxyy_yyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 105); 

                auto tg_xxxxyy_yyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 106); 

                auto tg_xxxxyy_yyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 107); 

                auto tg_xxxxyy_yyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 108); 

                auto tg_xxxxyy_yyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 109); 

                auto tg_xxxxyy_yzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 110); 

                auto tg_xxxxyy_zzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 111); 

                auto tg_xxxxyz_xxxxxx_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 112); 

                auto tg_xxxxyz_xxxxxy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 113); 

                auto tg_xxxxyz_xxxxxz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 114); 

                auto tg_xxxxyz_xxxxyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 115); 

                auto tg_xxxxyz_xxxxyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 116); 

                auto tg_xxxxyz_xxxxzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 117); 

                auto tg_xxxxyz_xxxyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 118); 

                auto tg_xxxxyz_xxxyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 119); 

                auto tg_xxxxyz_xxxyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 120); 

                auto tg_xxxxyz_xxxzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 121); 

                auto tg_xxxxyz_xxyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 122); 

                auto tg_xxxxyz_xxyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 123); 

                auto tg_xxxxyz_xxyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 124); 

                auto tg_xxxxyz_xxyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 125); 

                auto tg_xxxxyz_xxzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 126); 

                auto tg_xxxxyz_xyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 127); 

                auto tg_xxxxyz_xyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 128); 

                auto tg_xxxxyz_xyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 129); 

                auto tg_xxxxyz_xyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 130); 

                auto tg_xxxxyz_xyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 131); 

                auto tg_xxxxyz_xzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 132); 

                auto tg_xxxxyz_yyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 133); 

                auto tg_xxxxyz_yyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 134); 

                auto tg_xxxxyz_yyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 135); 

                auto tg_xxxxyz_yyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 136); 

                auto tg_xxxxyz_yyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 137); 

                auto tg_xxxxyz_yzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 138); 

                auto tg_xxxxyz_zzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 139); 

                auto tg_xxxxzz_xxxxxx_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 140); 

                auto tg_xxxxzz_xxxxxy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 141); 

                auto tg_xxxxzz_xxxxxz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 142); 

                auto tg_xxxxzz_xxxxyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 143); 

                auto tg_xxxxzz_xxxxyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 144); 

                auto tg_xxxxzz_xxxxzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 145); 

                auto tg_xxxxzz_xxxyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 146); 

                auto tg_xxxxzz_xxxyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 147); 

                auto tg_xxxxzz_xxxyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 148); 

                auto tg_xxxxzz_xxxzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 149); 

                auto tg_xxxxzz_xxyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 150); 

                auto tg_xxxxzz_xxyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 151); 

                auto tg_xxxxzz_xxyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 152); 

                auto tg_xxxxzz_xxyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 153); 

                auto tg_xxxxzz_xxzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 154); 

                auto tg_xxxxzz_xyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 155); 

                auto tg_xxxxzz_xyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 156); 

                auto tg_xxxxzz_xyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 157); 

                auto tg_xxxxzz_xyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 158); 

                auto tg_xxxxzz_xyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 159); 

                auto tg_xxxxzz_xzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 160); 

                auto tg_xxxxzz_yyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 161); 

                auto tg_xxxxzz_yyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 162); 

                auto tg_xxxxzz_yyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 163); 

                auto tg_xxxxzz_yyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 164); 

                auto tg_xxxxzz_yyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 165); 

                auto tg_xxxxzz_yzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 166); 

                auto tg_xxxxzz_zzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 167); 

                auto tg_xxxyyy_xxxxxx_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 168); 

                auto tg_xxxyyy_xxxxxy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 169); 

                auto tg_xxxyyy_xxxxxz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 170); 

                auto tg_xxxyyy_xxxxyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 171); 

                auto tg_xxxyyy_xxxxyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 172); 

                auto tg_xxxyyy_xxxxzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 173); 

                auto tg_xxxyyy_xxxyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 174); 

                auto tg_xxxyyy_xxxyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 175); 

                auto tg_xxxyyy_xxxyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 176); 

                auto tg_xxxyyy_xxxzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 177); 

                auto tg_xxxyyy_xxyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 178); 

                auto tg_xxxyyy_xxyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 179); 

                auto tg_xxxyyy_xxyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 180); 

                auto tg_xxxyyy_xxyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 181); 

                auto tg_xxxyyy_xxzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 182); 

                auto tg_xxxyyy_xyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 183); 

                auto tg_xxxyyy_xyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 184); 

                auto tg_xxxyyy_xyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 185); 

                auto tg_xxxyyy_xyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 186); 

                auto tg_xxxyyy_xyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 187); 

                auto tg_xxxyyy_xzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 188); 

                auto tg_xxxyyy_yyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 189); 

                auto tg_xxxyyy_yyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 190); 

                auto tg_xxxyyy_yyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 191); 

                auto tg_xxxyyy_yyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 192); 

                auto tg_xxxyyy_yyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 193); 

                auto tg_xxxyyy_yzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 194); 

                auto tg_xxxyyy_zzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 195); 

                // Batch of Integrals (98,196)

                #pragma omp simd aligned(fxn, fza, tg_xxxxyy_xxzzzz_0, tg_xxxxyy_xyyyyy_0, tg_xxxxyy_xyyyyz_0, \
                                         tg_xxxxyy_xyyyzz_0, tg_xxxxyy_xyyzzz_0, tg_xxxxyy_xyzzzz_0, tg_xxxxyy_xzzzzz_0, \
                                         tg_xxxxyy_yyyyyy_0, tg_xxxxyy_yyyyyz_0, tg_xxxxyy_yyyyzz_0, tg_xxxxyy_yyyzzz_0, \
                                         tg_xxxxyy_yyzzzz_0, tg_xxxxyy_yzzzzz_0, tg_xxxxyy_zzzzzz_0, tg_xxxxyz_xxxxxx_0, \
                                         tg_xxxxyz_xxxxxy_0, tg_xxxxyz_xxxxxz_0, tg_xxxxyz_xxxxyy_0, tg_xxxxyz_xxxxyz_0, \
                                         tg_xxxxyz_xxxxzz_0, tg_xxxxyz_xxxyyy_0, tg_xxxxyz_xxxyyz_0, tg_xxxxyz_xxxyzz_0, \
                                         tg_xxxxyz_xxxzzz_0, tg_xxxxyz_xxyyyy_0, tg_xxxxyz_xxyyyz_0, tg_xxxxyz_xxyyzz_0, \
                                         tg_xxxxyz_xxyzzz_0, tg_xxxxyz_xxzzzz_0, tg_xxxxyz_xyyyyy_0, tg_xxxxyz_xyyyyz_0, \
                                         tg_xxxxyz_xyyyzz_0, tg_xxxxyz_xyyzzz_0, tg_xxxxyz_xyzzzz_0, tg_xxxxyz_xzzzzz_0, \
                                         tg_xxxxyz_yyyyyy_0, tg_xxxxyz_yyyyyz_0, tg_xxxxyz_yyyyzz_0, tg_xxxxyz_yyyzzz_0, \
                                         tg_xxxxyz_yyzzzz_0, tg_xxxxyz_yzzzzz_0, tg_xxxxyz_zzzzzz_0, tg_xxxxzz_xxxxxx_0, \
                                         tg_xxxxzz_xxxxxy_0, tg_xxxxzz_xxxxxz_0, tg_xxxxzz_xxxxyy_0, tg_xxxxzz_xxxxyz_0, \
                                         tg_xxxxzz_xxxxzz_0, tg_xxxxzz_xxxyyy_0, tg_xxxxzz_xxxyyz_0, tg_xxxxzz_xxxyzz_0, \
                                         tg_xxxxzz_xxxzzz_0, tg_xxxxzz_xxyyyy_0, tg_xxxxzz_xxyyyz_0, tg_xxxxzz_xxyyzz_0, \
                                         tg_xxxxzz_xxyzzz_0, tg_xxxxzz_xxzzzz_0, tg_xxxxzz_xyyyyy_0, tg_xxxxzz_xyyyyz_0, \
                                         tg_xxxxzz_xyyyzz_0, tg_xxxxzz_xyyzzz_0, tg_xxxxzz_xyzzzz_0, tg_xxxxzz_xzzzzz_0, \
                                         tg_xxxxzz_yyyyyy_0, tg_xxxxzz_yyyyyz_0, tg_xxxxzz_yyyyzz_0, tg_xxxxzz_yyyzzz_0, \
                                         tg_xxxxzz_yyzzzz_0, tg_xxxxzz_yzzzzz_0, tg_xxxxzz_zzzzzz_0, tg_xxxyy_xxzzzz_0, \
                                         tg_xxxyy_xxzzzz_1, tg_xxxyy_xyyyyy_0, tg_xxxyy_xyyyyy_1, tg_xxxyy_xyyyyz_0, \
                                         tg_xxxyy_xyyyyz_1, tg_xxxyy_xyyyzz_0, tg_xxxyy_xyyyzz_1, tg_xxxyy_xyyzzz_0, \
                                         tg_xxxyy_xyyzzz_1, tg_xxxyy_xyzzzz_0, tg_xxxyy_xyzzzz_1, tg_xxxyy_xzzzz_1, \
                                         tg_xxxyy_xzzzzz_0, tg_xxxyy_xzzzzz_1, tg_xxxyy_yyyyy_1, tg_xxxyy_yyyyyy_0, \
                                         tg_xxxyy_yyyyyy_1, tg_xxxyy_yyyyyz_0, tg_xxxyy_yyyyyz_1, tg_xxxyy_yyyyz_1, \
                                         tg_xxxyy_yyyyzz_0, tg_xxxyy_yyyyzz_1, tg_xxxyy_yyyzz_1, tg_xxxyy_yyyzzz_0, \
                                         tg_xxxyy_yyyzzz_1, tg_xxxyy_yyzzz_1, tg_xxxyy_yyzzzz_0, tg_xxxyy_yyzzzz_1, \
                                         tg_xxxyy_yzzzz_1, tg_xxxyy_yzzzzz_0, tg_xxxyy_yzzzzz_1, tg_xxxyy_zzzzz_1, \
                                         tg_xxxyy_zzzzzz_0, tg_xxxyy_zzzzzz_1, tg_xxxyyy_xxxxxx_0, tg_xxxyyy_xxxxxy_0, \
                                         tg_xxxyyy_xxxxxz_0, tg_xxxyyy_xxxxyy_0, tg_xxxyyy_xxxxyz_0, tg_xxxyyy_xxxxzz_0, \
                                         tg_xxxyyy_xxxyyy_0, tg_xxxyyy_xxxyyz_0, tg_xxxyyy_xxxyzz_0, tg_xxxyyy_xxxzzz_0, \
                                         tg_xxxyyy_xxyyyy_0, tg_xxxyyy_xxyyyz_0, tg_xxxyyy_xxyyzz_0, tg_xxxyyy_xxyzzz_0, \
                                         tg_xxxyyy_xxzzzz_0, tg_xxxyyy_xyyyyy_0, tg_xxxyyy_xyyyyz_0, tg_xxxyyy_xyyyzz_0, \
                                         tg_xxxyyy_xyyzzz_0, tg_xxxyyy_xyzzzz_0, tg_xxxyyy_xzzzzz_0, tg_xxxyyy_yyyyyy_0, \
                                         tg_xxxyyy_yyyyyz_0, tg_xxxyyy_yyyyzz_0, tg_xxxyyy_yyyzzz_0, tg_xxxyyy_yyzzzz_0, \
                                         tg_xxxyyy_yzzzzz_0, tg_xxxyyy_zzzzzz_0, tg_xxxyz_xxxxx_1, tg_xxxyz_xxxxxx_0, \
                                         tg_xxxyz_xxxxxx_1, tg_xxxyz_xxxxxy_0, tg_xxxyz_xxxxxy_1, tg_xxxyz_xxxxxz_0, \
                                         tg_xxxyz_xxxxxz_1, tg_xxxyz_xxxxy_1, tg_xxxyz_xxxxyy_0, tg_xxxyz_xxxxyy_1, \
                                         tg_xxxyz_xxxxyz_0, tg_xxxyz_xxxxyz_1, tg_xxxyz_xxxxz_1, tg_xxxyz_xxxxzz_0, \
                                         tg_xxxyz_xxxxzz_1, tg_xxxyz_xxxyy_1, tg_xxxyz_xxxyyy_0, tg_xxxyz_xxxyyy_1, \
                                         tg_xxxyz_xxxyyz_0, tg_xxxyz_xxxyyz_1, tg_xxxyz_xxxyz_1, tg_xxxyz_xxxyzz_0, \
                                         tg_xxxyz_xxxyzz_1, tg_xxxyz_xxxzz_1, tg_xxxyz_xxxzzz_0, tg_xxxyz_xxxzzz_1, \
                                         tg_xxxyz_xxyyy_1, tg_xxxyz_xxyyyy_0, tg_xxxyz_xxyyyy_1, tg_xxxyz_xxyyyz_0, \
                                         tg_xxxyz_xxyyyz_1, tg_xxxyz_xxyyz_1, tg_xxxyz_xxyyzz_0, tg_xxxyz_xxyyzz_1, \
                                         tg_xxxyz_xxyzz_1, tg_xxxyz_xxyzzz_0, tg_xxxyz_xxyzzz_1, tg_xxxyz_xxzzz_1, \
                                         tg_xxxyz_xxzzzz_0, tg_xxxyz_xxzzzz_1, tg_xxxyz_xyyyy_1, tg_xxxyz_xyyyyy_0, \
                                         tg_xxxyz_xyyyyy_1, tg_xxxyz_xyyyyz_0, tg_xxxyz_xyyyyz_1, tg_xxxyz_xyyyz_1, \
                                         tg_xxxyz_xyyyzz_0, tg_xxxyz_xyyyzz_1, tg_xxxyz_xyyzz_1, tg_xxxyz_xyyzzz_0, \
                                         tg_xxxyz_xyyzzz_1, tg_xxxyz_xyzzz_1, tg_xxxyz_xyzzzz_0, tg_xxxyz_xyzzzz_1, \
                                         tg_xxxyz_xzzzz_1, tg_xxxyz_xzzzzz_0, tg_xxxyz_xzzzzz_1, tg_xxxyz_yyyyy_1, \
                                         tg_xxxyz_yyyyyy_0, tg_xxxyz_yyyyyy_1, tg_xxxyz_yyyyyz_0, tg_xxxyz_yyyyyz_1, \
                                         tg_xxxyz_yyyyz_1, tg_xxxyz_yyyyzz_0, tg_xxxyz_yyyyzz_1, tg_xxxyz_yyyzz_1, \
                                         tg_xxxyz_yyyzzz_0, tg_xxxyz_yyyzzz_1, tg_xxxyz_yyzzz_1, tg_xxxyz_yyzzzz_0, \
                                         tg_xxxyz_yyzzzz_1, tg_xxxyz_yzzzz_1, tg_xxxyz_yzzzzz_0, tg_xxxyz_yzzzzz_1, \
                                         tg_xxxyz_zzzzz_1, tg_xxxyz_zzzzzz_0, tg_xxxyz_zzzzzz_1, tg_xxxzz_xxxxx_1, \
                                         tg_xxxzz_xxxxxx_0, tg_xxxzz_xxxxxx_1, tg_xxxzz_xxxxxy_0, tg_xxxzz_xxxxxy_1, \
                                         tg_xxxzz_xxxxxz_0, tg_xxxzz_xxxxxz_1, tg_xxxzz_xxxxy_1, tg_xxxzz_xxxxyy_0, \
                                         tg_xxxzz_xxxxyy_1, tg_xxxzz_xxxxyz_0, tg_xxxzz_xxxxyz_1, tg_xxxzz_xxxxz_1, \
                                         tg_xxxzz_xxxxzz_0, tg_xxxzz_xxxxzz_1, tg_xxxzz_xxxyy_1, tg_xxxzz_xxxyyy_0, \
                                         tg_xxxzz_xxxyyy_1, tg_xxxzz_xxxyyz_0, tg_xxxzz_xxxyyz_1, tg_xxxzz_xxxyz_1, \
                                         tg_xxxzz_xxxyzz_0, tg_xxxzz_xxxyzz_1, tg_xxxzz_xxxzz_1, tg_xxxzz_xxxzzz_0, \
                                         tg_xxxzz_xxxzzz_1, tg_xxxzz_xxyyy_1, tg_xxxzz_xxyyyy_0, tg_xxxzz_xxyyyy_1, \
                                         tg_xxxzz_xxyyyz_0, tg_xxxzz_xxyyyz_1, tg_xxxzz_xxyyz_1, tg_xxxzz_xxyyzz_0, \
                                         tg_xxxzz_xxyyzz_1, tg_xxxzz_xxyzz_1, tg_xxxzz_xxyzzz_0, tg_xxxzz_xxyzzz_1, \
                                         tg_xxxzz_xxzzz_1, tg_xxxzz_xxzzzz_0, tg_xxxzz_xxzzzz_1, tg_xxxzz_xyyyy_1, \
                                         tg_xxxzz_xyyyyy_0, tg_xxxzz_xyyyyy_1, tg_xxxzz_xyyyyz_0, tg_xxxzz_xyyyyz_1, \
                                         tg_xxxzz_xyyyz_1, tg_xxxzz_xyyyzz_0, tg_xxxzz_xyyyzz_1, tg_xxxzz_xyyzz_1, \
                                         tg_xxxzz_xyyzzz_0, tg_xxxzz_xyyzzz_1, tg_xxxzz_xyzzz_1, tg_xxxzz_xyzzzz_0, \
                                         tg_xxxzz_xyzzzz_1, tg_xxxzz_xzzzz_1, tg_xxxzz_xzzzzz_0, tg_xxxzz_xzzzzz_1, \
                                         tg_xxxzz_yyyyy_1, tg_xxxzz_yyyyyy_0, tg_xxxzz_yyyyyy_1, tg_xxxzz_yyyyyz_0, \
                                         tg_xxxzz_yyyyyz_1, tg_xxxzz_yyyyz_1, tg_xxxzz_yyyyzz_0, tg_xxxzz_yyyyzz_1, \
                                         tg_xxxzz_yyyzz_1, tg_xxxzz_yyyzzz_0, tg_xxxzz_yyyzzz_1, tg_xxxzz_yyzzz_1, \
                                         tg_xxxzz_yyzzzz_0, tg_xxxzz_yyzzzz_1, tg_xxxzz_yzzzz_1, tg_xxxzz_yzzzzz_0, \
                                         tg_xxxzz_yzzzzz_1, tg_xxxzz_zzzzz_1, tg_xxxzz_zzzzzz_0, tg_xxxzz_zzzzzz_1, \
                                         tg_xxyy_xxzzzz_0, tg_xxyy_xxzzzz_1, tg_xxyy_xyyyyy_0, tg_xxyy_xyyyyy_1, \
                                         tg_xxyy_xyyyyz_0, tg_xxyy_xyyyyz_1, tg_xxyy_xyyyzz_0, tg_xxyy_xyyyzz_1, \
                                         tg_xxyy_xyyzzz_0, tg_xxyy_xyyzzz_1, tg_xxyy_xyzzzz_0, tg_xxyy_xyzzzz_1, \
                                         tg_xxyy_xzzzzz_0, tg_xxyy_xzzzzz_1, tg_xxyy_yyyyyy_0, tg_xxyy_yyyyyy_1, \
                                         tg_xxyy_yyyyyz_0, tg_xxyy_yyyyyz_1, tg_xxyy_yyyyzz_0, tg_xxyy_yyyyzz_1, \
                                         tg_xxyy_yyyzzz_0, tg_xxyy_yyyzzz_1, tg_xxyy_yyzzzz_0, tg_xxyy_yyzzzz_1, \
                                         tg_xxyy_yzzzzz_0, tg_xxyy_yzzzzz_1, tg_xxyy_zzzzzz_0, tg_xxyy_zzzzzz_1, \
                                         tg_xxyyy_xxxxx_1, tg_xxyyy_xxxxxx_0, tg_xxyyy_xxxxxx_1, tg_xxyyy_xxxxxy_0, \
                                         tg_xxyyy_xxxxxy_1, tg_xxyyy_xxxxxz_0, tg_xxyyy_xxxxxz_1, tg_xxyyy_xxxxy_1, \
                                         tg_xxyyy_xxxxyy_0, tg_xxyyy_xxxxyy_1, tg_xxyyy_xxxxyz_0, tg_xxyyy_xxxxyz_1, \
                                         tg_xxyyy_xxxxz_1, tg_xxyyy_xxxxzz_0, tg_xxyyy_xxxxzz_1, tg_xxyyy_xxxyy_1, \
                                         tg_xxyyy_xxxyyy_0, tg_xxyyy_xxxyyy_1, tg_xxyyy_xxxyyz_0, tg_xxyyy_xxxyyz_1, \
                                         tg_xxyyy_xxxyz_1, tg_xxyyy_xxxyzz_0, tg_xxyyy_xxxyzz_1, tg_xxyyy_xxxzz_1, \
                                         tg_xxyyy_xxxzzz_0, tg_xxyyy_xxxzzz_1, tg_xxyyy_xxyyy_1, tg_xxyyy_xxyyyy_0, \
                                         tg_xxyyy_xxyyyy_1, tg_xxyyy_xxyyyz_0, tg_xxyyy_xxyyyz_1, tg_xxyyy_xxyyz_1, \
                                         tg_xxyyy_xxyyzz_0, tg_xxyyy_xxyyzz_1, tg_xxyyy_xxyzz_1, tg_xxyyy_xxyzzz_0, \
                                         tg_xxyyy_xxyzzz_1, tg_xxyyy_xxzzz_1, tg_xxyyy_xxzzzz_0, tg_xxyyy_xxzzzz_1, \
                                         tg_xxyyy_xyyyy_1, tg_xxyyy_xyyyyy_0, tg_xxyyy_xyyyyy_1, tg_xxyyy_xyyyyz_0, \
                                         tg_xxyyy_xyyyyz_1, tg_xxyyy_xyyyz_1, tg_xxyyy_xyyyzz_0, tg_xxyyy_xyyyzz_1, \
                                         tg_xxyyy_xyyzz_1, tg_xxyyy_xyyzzz_0, tg_xxyyy_xyyzzz_1, tg_xxyyy_xyzzz_1, \
                                         tg_xxyyy_xyzzzz_0, tg_xxyyy_xyzzzz_1, tg_xxyyy_xzzzz_1, tg_xxyyy_xzzzzz_0, \
                                         tg_xxyyy_xzzzzz_1, tg_xxyyy_yyyyy_1, tg_xxyyy_yyyyyy_0, tg_xxyyy_yyyyyy_1, \
                                         tg_xxyyy_yyyyyz_0, tg_xxyyy_yyyyyz_1, tg_xxyyy_yyyyz_1, tg_xxyyy_yyyyzz_0, \
                                         tg_xxyyy_yyyyzz_1, tg_xxyyy_yyyzz_1, tg_xxyyy_yyyzzz_0, tg_xxyyy_yyyzzz_1, \
                                         tg_xxyyy_yyzzz_1, tg_xxyyy_yyzzzz_0, tg_xxyyy_yyzzzz_1, tg_xxyyy_yzzzz_1, \
                                         tg_xxyyy_yzzzzz_0, tg_xxyyy_yzzzzz_1, tg_xxyyy_zzzzz_1, tg_xxyyy_zzzzzz_0, \
                                         tg_xxyyy_zzzzzz_1, tg_xxyz_xxxxxx_0, tg_xxyz_xxxxxx_1, tg_xxyz_xxxxxy_0, \
                                         tg_xxyz_xxxxxy_1, tg_xxyz_xxxxxz_0, tg_xxyz_xxxxxz_1, tg_xxyz_xxxxyy_0, \
                                         tg_xxyz_xxxxyy_1, tg_xxyz_xxxxyz_0, tg_xxyz_xxxxyz_1, tg_xxyz_xxxxzz_0, \
                                         tg_xxyz_xxxxzz_1, tg_xxyz_xxxyyy_0, tg_xxyz_xxxyyy_1, tg_xxyz_xxxyyz_0, \
                                         tg_xxyz_xxxyyz_1, tg_xxyz_xxxyzz_0, tg_xxyz_xxxyzz_1, tg_xxyz_xxxzzz_0, \
                                         tg_xxyz_xxxzzz_1, tg_xxyz_xxyyyy_0, tg_xxyz_xxyyyy_1, tg_xxyz_xxyyyz_0, \
                                         tg_xxyz_xxyyyz_1, tg_xxyz_xxyyzz_0, tg_xxyz_xxyyzz_1, tg_xxyz_xxyzzz_0, \
                                         tg_xxyz_xxyzzz_1, tg_xxyz_xxzzzz_0, tg_xxyz_xxzzzz_1, tg_xxyz_xyyyyy_0, \
                                         tg_xxyz_xyyyyy_1, tg_xxyz_xyyyyz_0, tg_xxyz_xyyyyz_1, tg_xxyz_xyyyzz_0, \
                                         tg_xxyz_xyyyzz_1, tg_xxyz_xyyzzz_0, tg_xxyz_xyyzzz_1, tg_xxyz_xyzzzz_0, \
                                         tg_xxyz_xyzzzz_1, tg_xxyz_xzzzzz_0, tg_xxyz_xzzzzz_1, tg_xxyz_yyyyyy_0, \
                                         tg_xxyz_yyyyyy_1, tg_xxyz_yyyyyz_0, tg_xxyz_yyyyyz_1, tg_xxyz_yyyyzz_0, \
                                         tg_xxyz_yyyyzz_1, tg_xxyz_yyyzzz_0, tg_xxyz_yyyzzz_1, tg_xxyz_yyzzzz_0, \
                                         tg_xxyz_yyzzzz_1, tg_xxyz_yzzzzz_0, tg_xxyz_yzzzzz_1, tg_xxyz_zzzzzz_0, \
                                         tg_xxyz_zzzzzz_1, tg_xxzz_xxxxxx_0, tg_xxzz_xxxxxx_1, tg_xxzz_xxxxxy_0, \
                                         tg_xxzz_xxxxxy_1, tg_xxzz_xxxxxz_0, tg_xxzz_xxxxxz_1, tg_xxzz_xxxxyy_0, \
                                         tg_xxzz_xxxxyy_1, tg_xxzz_xxxxyz_0, tg_xxzz_xxxxyz_1, tg_xxzz_xxxxzz_0, \
                                         tg_xxzz_xxxxzz_1, tg_xxzz_xxxyyy_0, tg_xxzz_xxxyyy_1, tg_xxzz_xxxyyz_0, \
                                         tg_xxzz_xxxyyz_1, tg_xxzz_xxxyzz_0, tg_xxzz_xxxyzz_1, tg_xxzz_xxxzzz_0, \
                                         tg_xxzz_xxxzzz_1, tg_xxzz_xxyyyy_0, tg_xxzz_xxyyyy_1, tg_xxzz_xxyyyz_0, \
                                         tg_xxzz_xxyyyz_1, tg_xxzz_xxyyzz_0, tg_xxzz_xxyyzz_1, tg_xxzz_xxyzzz_0, \
                                         tg_xxzz_xxyzzz_1, tg_xxzz_xxzzzz_0, tg_xxzz_xxzzzz_1, tg_xxzz_xyyyyy_0, \
                                         tg_xxzz_xyyyyy_1, tg_xxzz_xyyyyz_0, tg_xxzz_xyyyyz_1, tg_xxzz_xyyyzz_0, \
                                         tg_xxzz_xyyyzz_1, tg_xxzz_xyyzzz_0, tg_xxzz_xyyzzz_1, tg_xxzz_xyzzzz_0, \
                                         tg_xxzz_xyzzzz_1, tg_xxzz_xzzzzz_0, tg_xxzz_xzzzzz_1, tg_xxzz_yyyyyy_0, \
                                         tg_xxzz_yyyyyy_1, tg_xxzz_yyyyyz_0, tg_xxzz_yyyyyz_1, tg_xxzz_yyyyzz_0, \
                                         tg_xxzz_yyyyzz_1, tg_xxzz_yyyzzz_0, tg_xxzz_yyyzzz_1, tg_xxzz_yyzzzz_0, \
                                         tg_xxzz_yyzzzz_1, tg_xxzz_yzzzzz_0, tg_xxzz_yzzzzz_1, tg_xxzz_zzzzzz_0, \
                                         tg_xxzz_zzzzzz_1, tg_xyyy_xxxxxx_0, tg_xyyy_xxxxxx_1, tg_xyyy_xxxxxy_0, \
                                         tg_xyyy_xxxxxy_1, tg_xyyy_xxxxxz_0, tg_xyyy_xxxxxz_1, tg_xyyy_xxxxyy_0, \
                                         tg_xyyy_xxxxyy_1, tg_xyyy_xxxxyz_0, tg_xyyy_xxxxyz_1, tg_xyyy_xxxxzz_0, \
                                         tg_xyyy_xxxxzz_1, tg_xyyy_xxxyyy_0, tg_xyyy_xxxyyy_1, tg_xyyy_xxxyyz_0, \
                                         tg_xyyy_xxxyyz_1, tg_xyyy_xxxyzz_0, tg_xyyy_xxxyzz_1, tg_xyyy_xxxzzz_0, \
                                         tg_xyyy_xxxzzz_1, tg_xyyy_xxyyyy_0, tg_xyyy_xxyyyy_1, tg_xyyy_xxyyyz_0, \
                                         tg_xyyy_xxyyyz_1, tg_xyyy_xxyyzz_0, tg_xyyy_xxyyzz_1, tg_xyyy_xxyzzz_0, \
                                         tg_xyyy_xxyzzz_1, tg_xyyy_xxzzzz_0, tg_xyyy_xxzzzz_1, tg_xyyy_xyyyyy_0, \
                                         tg_xyyy_xyyyyy_1, tg_xyyy_xyyyyz_0, tg_xyyy_xyyyyz_1, tg_xyyy_xyyyzz_0, \
                                         tg_xyyy_xyyyzz_1, tg_xyyy_xyyzzz_0, tg_xyyy_xyyzzz_1, tg_xyyy_xyzzzz_0, \
                                         tg_xyyy_xyzzzz_1, tg_xyyy_xzzzzz_0, tg_xyyy_xzzzzz_1, tg_xyyy_yyyyyy_0, \
                                         tg_xyyy_yyyyyy_1, tg_xyyy_yyyyyz_0, tg_xyyy_yyyyyz_1, tg_xyyy_yyyyzz_0, \
                                         tg_xyyy_yyyyzz_1, tg_xyyy_yyyzzz_0, tg_xyyy_yyyzzz_1, tg_xyyy_yyzzzz_0, \
                                         tg_xyyy_yyzzzz_1, tg_xyyy_yzzzzz_0, tg_xyyy_yzzzzz_1, tg_xyyy_zzzzzz_0, \
                                         tg_xyyy_zzzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxxyy_xxzzzz_0[j] = pb_x * tg_xxxyy_xxzzzz_0[j] + fr * tg_xxxyy_xxzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyy_xxzzzz_0[j] - tg_xxyy_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyy_xzzzz_1[j];

                    tg_xxxxyy_xyyyyy_0[j] = pb_x * tg_xxxyy_xyyyyy_0[j] + fr * tg_xxxyy_xyyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyy_xyyyyy_0[j] - tg_xxyy_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyy_yyyyy_1[j];

                    tg_xxxxyy_xyyyyz_0[j] = pb_x * tg_xxxyy_xyyyyz_0[j] + fr * tg_xxxyy_xyyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyy_xyyyyz_0[j] - tg_xxyy_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyy_yyyyz_1[j];

                    tg_xxxxyy_xyyyzz_0[j] = pb_x * tg_xxxyy_xyyyzz_0[j] + fr * tg_xxxyy_xyyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyy_xyyyzz_0[j] - tg_xxyy_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyy_yyyzz_1[j];

                    tg_xxxxyy_xyyzzz_0[j] = pb_x * tg_xxxyy_xyyzzz_0[j] + fr * tg_xxxyy_xyyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyy_xyyzzz_0[j] - tg_xxyy_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyy_yyzzz_1[j];

                    tg_xxxxyy_xyzzzz_0[j] = pb_x * tg_xxxyy_xyzzzz_0[j] + fr * tg_xxxyy_xyzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyy_xyzzzz_0[j] - tg_xxyy_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyy_yzzzz_1[j];

                    tg_xxxxyy_xzzzzz_0[j] = pb_x * tg_xxxyy_xzzzzz_0[j] + fr * tg_xxxyy_xzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyy_xzzzzz_0[j] - tg_xxyy_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyy_zzzzz_1[j];

                    tg_xxxxyy_yyyyyy_0[j] = pb_x * tg_xxxyy_yyyyyy_0[j] + fr * tg_xxxyy_yyyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyy_yyyyyy_0[j] - tg_xxyy_yyyyyy_1[j] * fl1_fza);

                    tg_xxxxyy_yyyyyz_0[j] = pb_x * tg_xxxyy_yyyyyz_0[j] + fr * tg_xxxyy_yyyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyy_yyyyyz_0[j] - tg_xxyy_yyyyyz_1[j] * fl1_fza);

                    tg_xxxxyy_yyyyzz_0[j] = pb_x * tg_xxxyy_yyyyzz_0[j] + fr * tg_xxxyy_yyyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyy_yyyyzz_0[j] - tg_xxyy_yyyyzz_1[j] * fl1_fza);

                    tg_xxxxyy_yyyzzz_0[j] = pb_x * tg_xxxyy_yyyzzz_0[j] + fr * tg_xxxyy_yyyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyy_yyyzzz_0[j] - tg_xxyy_yyyzzz_1[j] * fl1_fza);

                    tg_xxxxyy_yyzzzz_0[j] = pb_x * tg_xxxyy_yyzzzz_0[j] + fr * tg_xxxyy_yyzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyy_yyzzzz_0[j] - tg_xxyy_yyzzzz_1[j] * fl1_fza);

                    tg_xxxxyy_yzzzzz_0[j] = pb_x * tg_xxxyy_yzzzzz_0[j] + fr * tg_xxxyy_yzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyy_yzzzzz_0[j] - tg_xxyy_yzzzzz_1[j] * fl1_fza);

                    tg_xxxxyy_zzzzzz_0[j] = pb_x * tg_xxxyy_zzzzzz_0[j] + fr * tg_xxxyy_zzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyy_zzzzzz_0[j] - tg_xxyy_zzzzzz_1[j] * fl1_fza);

                    tg_xxxxyz_xxxxxx_0[j] = pb_x * tg_xxxyz_xxxxxx_0[j] + fr * tg_xxxyz_xxxxxx_1[j] + 1.5 * fl1_fx * (tg_xxyz_xxxxxx_0[j] - tg_xxyz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxxyz_xxxxx_1[j];

                    tg_xxxxyz_xxxxxy_0[j] = pb_x * tg_xxxyz_xxxxxy_0[j] + fr * tg_xxxyz_xxxxxy_1[j] + 1.5 * fl1_fx * (tg_xxyz_xxxxxy_0[j] - tg_xxyz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxyz_xxxxy_1[j];

                    tg_xxxxyz_xxxxxz_0[j] = pb_x * tg_xxxyz_xxxxxz_0[j] + fr * tg_xxxyz_xxxxxz_1[j] + 1.5 * fl1_fx * (tg_xxyz_xxxxxz_0[j] - tg_xxyz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxyz_xxxxz_1[j];

                    tg_xxxxyz_xxxxyy_0[j] = pb_x * tg_xxxyz_xxxxyy_0[j] + fr * tg_xxxyz_xxxxyy_1[j] + 1.5 * fl1_fx * (tg_xxyz_xxxxyy_0[j] - tg_xxyz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyz_xxxyy_1[j];

                    tg_xxxxyz_xxxxyz_0[j] = pb_x * tg_xxxyz_xxxxyz_0[j] + fr * tg_xxxyz_xxxxyz_1[j] + 1.5 * fl1_fx * (tg_xxyz_xxxxyz_0[j] - tg_xxyz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyz_xxxyz_1[j];

                    tg_xxxxyz_xxxxzz_0[j] = pb_x * tg_xxxyz_xxxxzz_0[j] + fr * tg_xxxyz_xxxxzz_1[j] + 1.5 * fl1_fx * (tg_xxyz_xxxxzz_0[j] - tg_xxyz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyz_xxxzz_1[j];

                    tg_xxxxyz_xxxyyy_0[j] = pb_x * tg_xxxyz_xxxyyy_0[j] + fr * tg_xxxyz_xxxyyy_1[j] + 1.5 * fl1_fx * (tg_xxyz_xxxyyy_0[j] - tg_xxyz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyz_xxyyy_1[j];

                    tg_xxxxyz_xxxyyz_0[j] = pb_x * tg_xxxyz_xxxyyz_0[j] + fr * tg_xxxyz_xxxyyz_1[j] + 1.5 * fl1_fx * (tg_xxyz_xxxyyz_0[j] - tg_xxyz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyz_xxyyz_1[j];

                    tg_xxxxyz_xxxyzz_0[j] = pb_x * tg_xxxyz_xxxyzz_0[j] + fr * tg_xxxyz_xxxyzz_1[j] + 1.5 * fl1_fx * (tg_xxyz_xxxyzz_0[j] - tg_xxyz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyz_xxyzz_1[j];

                    tg_xxxxyz_xxxzzz_0[j] = pb_x * tg_xxxyz_xxxzzz_0[j] + fr * tg_xxxyz_xxxzzz_1[j] + 1.5 * fl1_fx * (tg_xxyz_xxxzzz_0[j] - tg_xxyz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyz_xxzzz_1[j];

                    tg_xxxxyz_xxyyyy_0[j] = pb_x * tg_xxxyz_xxyyyy_0[j] + fr * tg_xxxyz_xxyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyz_xxyyyy_0[j] - tg_xxyz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxyz_xyyyy_1[j];

                    tg_xxxxyz_xxyyyz_0[j] = pb_x * tg_xxxyz_xxyyyz_0[j] + fr * tg_xxxyz_xxyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyz_xxyyyz_0[j] - tg_xxyz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyz_xyyyz_1[j];

                    tg_xxxxyz_xxyyzz_0[j] = pb_x * tg_xxxyz_xxyyzz_0[j] + fr * tg_xxxyz_xxyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyz_xxyyzz_0[j] - tg_xxyz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyz_xyyzz_1[j];

                    tg_xxxxyz_xxyzzz_0[j] = pb_x * tg_xxxyz_xxyzzz_0[j] + fr * tg_xxxyz_xxyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyz_xxyzzz_0[j] - tg_xxyz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyz_xyzzz_1[j];

                    tg_xxxxyz_xxzzzz_0[j] = pb_x * tg_xxxyz_xxzzzz_0[j] + fr * tg_xxxyz_xxzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyz_xxzzzz_0[j] - tg_xxyz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyz_xzzzz_1[j];

                    tg_xxxxyz_xyyyyy_0[j] = pb_x * tg_xxxyz_xyyyyy_0[j] + fr * tg_xxxyz_xyyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyz_xyyyyy_0[j] - tg_xxyz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyz_yyyyy_1[j];

                    tg_xxxxyz_xyyyyz_0[j] = pb_x * tg_xxxyz_xyyyyz_0[j] + fr * tg_xxxyz_xyyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyz_xyyyyz_0[j] - tg_xxyz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyz_yyyyz_1[j];

                    tg_xxxxyz_xyyyzz_0[j] = pb_x * tg_xxxyz_xyyyzz_0[j] + fr * tg_xxxyz_xyyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyz_xyyyzz_0[j] - tg_xxyz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyz_yyyzz_1[j];

                    tg_xxxxyz_xyyzzz_0[j] = pb_x * tg_xxxyz_xyyzzz_0[j] + fr * tg_xxxyz_xyyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyz_xyyzzz_0[j] - tg_xxyz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyz_yyzzz_1[j];

                    tg_xxxxyz_xyzzzz_0[j] = pb_x * tg_xxxyz_xyzzzz_0[j] + fr * tg_xxxyz_xyzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyz_xyzzzz_0[j] - tg_xxyz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyz_yzzzz_1[j];

                    tg_xxxxyz_xzzzzz_0[j] = pb_x * tg_xxxyz_xzzzzz_0[j] + fr * tg_xxxyz_xzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyz_xzzzzz_0[j] - tg_xxyz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyz_zzzzz_1[j];

                    tg_xxxxyz_yyyyyy_0[j] = pb_x * tg_xxxyz_yyyyyy_0[j] + fr * tg_xxxyz_yyyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyz_yyyyyy_0[j] - tg_xxyz_yyyyyy_1[j] * fl1_fza);

                    tg_xxxxyz_yyyyyz_0[j] = pb_x * tg_xxxyz_yyyyyz_0[j] + fr * tg_xxxyz_yyyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyz_yyyyyz_0[j] - tg_xxyz_yyyyyz_1[j] * fl1_fza);

                    tg_xxxxyz_yyyyzz_0[j] = pb_x * tg_xxxyz_yyyyzz_0[j] + fr * tg_xxxyz_yyyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyz_yyyyzz_0[j] - tg_xxyz_yyyyzz_1[j] * fl1_fza);

                    tg_xxxxyz_yyyzzz_0[j] = pb_x * tg_xxxyz_yyyzzz_0[j] + fr * tg_xxxyz_yyyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyz_yyyzzz_0[j] - tg_xxyz_yyyzzz_1[j] * fl1_fza);

                    tg_xxxxyz_yyzzzz_0[j] = pb_x * tg_xxxyz_yyzzzz_0[j] + fr * tg_xxxyz_yyzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyz_yyzzzz_0[j] - tg_xxyz_yyzzzz_1[j] * fl1_fza);

                    tg_xxxxyz_yzzzzz_0[j] = pb_x * tg_xxxyz_yzzzzz_0[j] + fr * tg_xxxyz_yzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyz_yzzzzz_0[j] - tg_xxyz_yzzzzz_1[j] * fl1_fza);

                    tg_xxxxyz_zzzzzz_0[j] = pb_x * tg_xxxyz_zzzzzz_0[j] + fr * tg_xxxyz_zzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyz_zzzzzz_0[j] - tg_xxyz_zzzzzz_1[j] * fl1_fza);

                    tg_xxxxzz_xxxxxx_0[j] = pb_x * tg_xxxzz_xxxxxx_0[j] + fr * tg_xxxzz_xxxxxx_1[j] + 1.5 * fl1_fx * (tg_xxzz_xxxxxx_0[j] - tg_xxzz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxxzz_xxxxx_1[j];

                    tg_xxxxzz_xxxxxy_0[j] = pb_x * tg_xxxzz_xxxxxy_0[j] + fr * tg_xxxzz_xxxxxy_1[j] + 1.5 * fl1_fx * (tg_xxzz_xxxxxy_0[j] - tg_xxzz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxzz_xxxxy_1[j];

                    tg_xxxxzz_xxxxxz_0[j] = pb_x * tg_xxxzz_xxxxxz_0[j] + fr * tg_xxxzz_xxxxxz_1[j] + 1.5 * fl1_fx * (tg_xxzz_xxxxxz_0[j] - tg_xxzz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxzz_xxxxz_1[j];

                    tg_xxxxzz_xxxxyy_0[j] = pb_x * tg_xxxzz_xxxxyy_0[j] + fr * tg_xxxzz_xxxxyy_1[j] + 1.5 * fl1_fx * (tg_xxzz_xxxxyy_0[j] - tg_xxzz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxzz_xxxyy_1[j];

                    tg_xxxxzz_xxxxyz_0[j] = pb_x * tg_xxxzz_xxxxyz_0[j] + fr * tg_xxxzz_xxxxyz_1[j] + 1.5 * fl1_fx * (tg_xxzz_xxxxyz_0[j] - tg_xxzz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxzz_xxxyz_1[j];

                    tg_xxxxzz_xxxxzz_0[j] = pb_x * tg_xxxzz_xxxxzz_0[j] + fr * tg_xxxzz_xxxxzz_1[j] + 1.5 * fl1_fx * (tg_xxzz_xxxxzz_0[j] - tg_xxzz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxzz_xxxzz_1[j];

                    tg_xxxxzz_xxxyyy_0[j] = pb_x * tg_xxxzz_xxxyyy_0[j] + fr * tg_xxxzz_xxxyyy_1[j] + 1.5 * fl1_fx * (tg_xxzz_xxxyyy_0[j] - tg_xxzz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxzz_xxyyy_1[j];

                    tg_xxxxzz_xxxyyz_0[j] = pb_x * tg_xxxzz_xxxyyz_0[j] + fr * tg_xxxzz_xxxyyz_1[j] + 1.5 * fl1_fx * (tg_xxzz_xxxyyz_0[j] - tg_xxzz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxzz_xxyyz_1[j];

                    tg_xxxxzz_xxxyzz_0[j] = pb_x * tg_xxxzz_xxxyzz_0[j] + fr * tg_xxxzz_xxxyzz_1[j] + 1.5 * fl1_fx * (tg_xxzz_xxxyzz_0[j] - tg_xxzz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxzz_xxyzz_1[j];

                    tg_xxxxzz_xxxzzz_0[j] = pb_x * tg_xxxzz_xxxzzz_0[j] + fr * tg_xxxzz_xxxzzz_1[j] + 1.5 * fl1_fx * (tg_xxzz_xxxzzz_0[j] - tg_xxzz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxzz_xxzzz_1[j];

                    tg_xxxxzz_xxyyyy_0[j] = pb_x * tg_xxxzz_xxyyyy_0[j] + fr * tg_xxxzz_xxyyyy_1[j] + 1.5 * fl1_fx * (tg_xxzz_xxyyyy_0[j] - tg_xxzz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxzz_xyyyy_1[j];

                    tg_xxxxzz_xxyyyz_0[j] = pb_x * tg_xxxzz_xxyyyz_0[j] + fr * tg_xxxzz_xxyyyz_1[j] + 1.5 * fl1_fx * (tg_xxzz_xxyyyz_0[j] - tg_xxzz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxzz_xyyyz_1[j];

                    tg_xxxxzz_xxyyzz_0[j] = pb_x * tg_xxxzz_xxyyzz_0[j] + fr * tg_xxxzz_xxyyzz_1[j] + 1.5 * fl1_fx * (tg_xxzz_xxyyzz_0[j] - tg_xxzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxzz_xyyzz_1[j];

                    tg_xxxxzz_xxyzzz_0[j] = pb_x * tg_xxxzz_xxyzzz_0[j] + fr * tg_xxxzz_xxyzzz_1[j] + 1.5 * fl1_fx * (tg_xxzz_xxyzzz_0[j] - tg_xxzz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxzz_xyzzz_1[j];

                    tg_xxxxzz_xxzzzz_0[j] = pb_x * tg_xxxzz_xxzzzz_0[j] + fr * tg_xxxzz_xxzzzz_1[j] + 1.5 * fl1_fx * (tg_xxzz_xxzzzz_0[j] - tg_xxzz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxzz_xzzzz_1[j];

                    tg_xxxxzz_xyyyyy_0[j] = pb_x * tg_xxxzz_xyyyyy_0[j] + fr * tg_xxxzz_xyyyyy_1[j] + 1.5 * fl1_fx * (tg_xxzz_xyyyyy_0[j] - tg_xxzz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzz_yyyyy_1[j];

                    tg_xxxxzz_xyyyyz_0[j] = pb_x * tg_xxxzz_xyyyyz_0[j] + fr * tg_xxxzz_xyyyyz_1[j] + 1.5 * fl1_fx * (tg_xxzz_xyyyyz_0[j] - tg_xxzz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzz_yyyyz_1[j];

                    tg_xxxxzz_xyyyzz_0[j] = pb_x * tg_xxxzz_xyyyzz_0[j] + fr * tg_xxxzz_xyyyzz_1[j] + 1.5 * fl1_fx * (tg_xxzz_xyyyzz_0[j] - tg_xxzz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzz_yyyzz_1[j];

                    tg_xxxxzz_xyyzzz_0[j] = pb_x * tg_xxxzz_xyyzzz_0[j] + fr * tg_xxxzz_xyyzzz_1[j] + 1.5 * fl1_fx * (tg_xxzz_xyyzzz_0[j] - tg_xxzz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzz_yyzzz_1[j];

                    tg_xxxxzz_xyzzzz_0[j] = pb_x * tg_xxxzz_xyzzzz_0[j] + fr * tg_xxxzz_xyzzzz_1[j] + 1.5 * fl1_fx * (tg_xxzz_xyzzzz_0[j] - tg_xxzz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzz_yzzzz_1[j];

                    tg_xxxxzz_xzzzzz_0[j] = pb_x * tg_xxxzz_xzzzzz_0[j] + fr * tg_xxxzz_xzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxzz_xzzzzz_0[j] - tg_xxzz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzz_zzzzz_1[j];

                    tg_xxxxzz_yyyyyy_0[j] = pb_x * tg_xxxzz_yyyyyy_0[j] + fr * tg_xxxzz_yyyyyy_1[j] + 1.5 * fl1_fx * (tg_xxzz_yyyyyy_0[j] - tg_xxzz_yyyyyy_1[j] * fl1_fza);

                    tg_xxxxzz_yyyyyz_0[j] = pb_x * tg_xxxzz_yyyyyz_0[j] + fr * tg_xxxzz_yyyyyz_1[j] + 1.5 * fl1_fx * (tg_xxzz_yyyyyz_0[j] - tg_xxzz_yyyyyz_1[j] * fl1_fza);

                    tg_xxxxzz_yyyyzz_0[j] = pb_x * tg_xxxzz_yyyyzz_0[j] + fr * tg_xxxzz_yyyyzz_1[j] + 1.5 * fl1_fx * (tg_xxzz_yyyyzz_0[j] - tg_xxzz_yyyyzz_1[j] * fl1_fza);

                    tg_xxxxzz_yyyzzz_0[j] = pb_x * tg_xxxzz_yyyzzz_0[j] + fr * tg_xxxzz_yyyzzz_1[j] + 1.5 * fl1_fx * (tg_xxzz_yyyzzz_0[j] - tg_xxzz_yyyzzz_1[j] * fl1_fza);

                    tg_xxxxzz_yyzzzz_0[j] = pb_x * tg_xxxzz_yyzzzz_0[j] + fr * tg_xxxzz_yyzzzz_1[j] + 1.5 * fl1_fx * (tg_xxzz_yyzzzz_0[j] - tg_xxzz_yyzzzz_1[j] * fl1_fza);

                    tg_xxxxzz_yzzzzz_0[j] = pb_x * tg_xxxzz_yzzzzz_0[j] + fr * tg_xxxzz_yzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxzz_yzzzzz_0[j] - tg_xxzz_yzzzzz_1[j] * fl1_fza);

                    tg_xxxxzz_zzzzzz_0[j] = pb_x * tg_xxxzz_zzzzzz_0[j] + fr * tg_xxxzz_zzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxzz_zzzzzz_0[j] - tg_xxzz_zzzzzz_1[j] * fl1_fza);

                    tg_xxxyyy_xxxxxx_0[j] = pb_x * tg_xxyyy_xxxxxx_0[j] + fr * tg_xxyyy_xxxxxx_1[j] + fl1_fx * (tg_xyyy_xxxxxx_0[j] - tg_xyyy_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxyyy_xxxxx_1[j];

                    tg_xxxyyy_xxxxxy_0[j] = pb_x * tg_xxyyy_xxxxxy_0[j] + fr * tg_xxyyy_xxxxxy_1[j] + fl1_fx * (tg_xyyy_xxxxxy_0[j] - tg_xyyy_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyyy_xxxxy_1[j];

                    tg_xxxyyy_xxxxxz_0[j] = pb_x * tg_xxyyy_xxxxxz_0[j] + fr * tg_xxyyy_xxxxxz_1[j] + fl1_fx * (tg_xyyy_xxxxxz_0[j] - tg_xyyy_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyyy_xxxxz_1[j];

                    tg_xxxyyy_xxxxyy_0[j] = pb_x * tg_xxyyy_xxxxyy_0[j] + fr * tg_xxyyy_xxxxyy_1[j] + fl1_fx * (tg_xyyy_xxxxyy_0[j] - tg_xyyy_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyy_xxxyy_1[j];

                    tg_xxxyyy_xxxxyz_0[j] = pb_x * tg_xxyyy_xxxxyz_0[j] + fr * tg_xxyyy_xxxxyz_1[j] + fl1_fx * (tg_xyyy_xxxxyz_0[j] - tg_xyyy_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyy_xxxyz_1[j];

                    tg_xxxyyy_xxxxzz_0[j] = pb_x * tg_xxyyy_xxxxzz_0[j] + fr * tg_xxyyy_xxxxzz_1[j] + fl1_fx * (tg_xyyy_xxxxzz_0[j] - tg_xyyy_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyy_xxxzz_1[j];

                    tg_xxxyyy_xxxyyy_0[j] = pb_x * tg_xxyyy_xxxyyy_0[j] + fr * tg_xxyyy_xxxyyy_1[j] + fl1_fx * (tg_xyyy_xxxyyy_0[j] - tg_xyyy_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyy_xxyyy_1[j];

                    tg_xxxyyy_xxxyyz_0[j] = pb_x * tg_xxyyy_xxxyyz_0[j] + fr * tg_xxyyy_xxxyyz_1[j] + fl1_fx * (tg_xyyy_xxxyyz_0[j] - tg_xyyy_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyy_xxyyz_1[j];

                    tg_xxxyyy_xxxyzz_0[j] = pb_x * tg_xxyyy_xxxyzz_0[j] + fr * tg_xxyyy_xxxyzz_1[j] + fl1_fx * (tg_xyyy_xxxyzz_0[j] - tg_xyyy_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyy_xxyzz_1[j];

                    tg_xxxyyy_xxxzzz_0[j] = pb_x * tg_xxyyy_xxxzzz_0[j] + fr * tg_xxyyy_xxxzzz_1[j] + fl1_fx * (tg_xyyy_xxxzzz_0[j] - tg_xyyy_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyy_xxzzz_1[j];

                    tg_xxxyyy_xxyyyy_0[j] = pb_x * tg_xxyyy_xxyyyy_0[j] + fr * tg_xxyyy_xxyyyy_1[j] + fl1_fx * (tg_xyyy_xxyyyy_0[j] - tg_xyyy_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxyyy_xyyyy_1[j];

                    tg_xxxyyy_xxyyyz_0[j] = pb_x * tg_xxyyy_xxyyyz_0[j] + fr * tg_xxyyy_xxyyyz_1[j] + fl1_fx * (tg_xyyy_xxyyyz_0[j] - tg_xyyy_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyy_xyyyz_1[j];

                    tg_xxxyyy_xxyyzz_0[j] = pb_x * tg_xxyyy_xxyyzz_0[j] + fr * tg_xxyyy_xxyyzz_1[j] + fl1_fx * (tg_xyyy_xxyyzz_0[j] - tg_xyyy_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyy_xyyzz_1[j];

                    tg_xxxyyy_xxyzzz_0[j] = pb_x * tg_xxyyy_xxyzzz_0[j] + fr * tg_xxyyy_xxyzzz_1[j] + fl1_fx * (tg_xyyy_xxyzzz_0[j] - tg_xyyy_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyy_xyzzz_1[j];

                    tg_xxxyyy_xxzzzz_0[j] = pb_x * tg_xxyyy_xxzzzz_0[j] + fr * tg_xxyyy_xxzzzz_1[j] + fl1_fx * (tg_xyyy_xxzzzz_0[j] - tg_xyyy_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyy_xzzzz_1[j];

                    tg_xxxyyy_xyyyyy_0[j] = pb_x * tg_xxyyy_xyyyyy_0[j] + fr * tg_xxyyy_xyyyyy_1[j] + fl1_fx * (tg_xyyy_xyyyyy_0[j] - tg_xyyy_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyy_yyyyy_1[j];

                    tg_xxxyyy_xyyyyz_0[j] = pb_x * tg_xxyyy_xyyyyz_0[j] + fr * tg_xxyyy_xyyyyz_1[j] + fl1_fx * (tg_xyyy_xyyyyz_0[j] - tg_xyyy_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyy_yyyyz_1[j];

                    tg_xxxyyy_xyyyzz_0[j] = pb_x * tg_xxyyy_xyyyzz_0[j] + fr * tg_xxyyy_xyyyzz_1[j] + fl1_fx * (tg_xyyy_xyyyzz_0[j] - tg_xyyy_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyy_yyyzz_1[j];

                    tg_xxxyyy_xyyzzz_0[j] = pb_x * tg_xxyyy_xyyzzz_0[j] + fr * tg_xxyyy_xyyzzz_1[j] + fl1_fx * (tg_xyyy_xyyzzz_0[j] - tg_xyyy_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyy_yyzzz_1[j];

                    tg_xxxyyy_xyzzzz_0[j] = pb_x * tg_xxyyy_xyzzzz_0[j] + fr * tg_xxyyy_xyzzzz_1[j] + fl1_fx * (tg_xyyy_xyzzzz_0[j] - tg_xyyy_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyy_yzzzz_1[j];

                    tg_xxxyyy_xzzzzz_0[j] = pb_x * tg_xxyyy_xzzzzz_0[j] + fr * tg_xxyyy_xzzzzz_1[j] + fl1_fx * (tg_xyyy_xzzzzz_0[j] - tg_xyyy_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyy_zzzzz_1[j];

                    tg_xxxyyy_yyyyyy_0[j] = pb_x * tg_xxyyy_yyyyyy_0[j] + fr * tg_xxyyy_yyyyyy_1[j] + fl1_fx * (tg_xyyy_yyyyyy_0[j] - tg_xyyy_yyyyyy_1[j] * fl1_fza);

                    tg_xxxyyy_yyyyyz_0[j] = pb_x * tg_xxyyy_yyyyyz_0[j] + fr * tg_xxyyy_yyyyyz_1[j] + fl1_fx * (tg_xyyy_yyyyyz_0[j] - tg_xyyy_yyyyyz_1[j] * fl1_fza);

                    tg_xxxyyy_yyyyzz_0[j] = pb_x * tg_xxyyy_yyyyzz_0[j] + fr * tg_xxyyy_yyyyzz_1[j] + fl1_fx * (tg_xyyy_yyyyzz_0[j] - tg_xyyy_yyyyzz_1[j] * fl1_fza);

                    tg_xxxyyy_yyyzzz_0[j] = pb_x * tg_xxyyy_yyyzzz_0[j] + fr * tg_xxyyy_yyyzzz_1[j] + fl1_fx * (tg_xyyy_yyyzzz_0[j] - tg_xyyy_yyyzzz_1[j] * fl1_fza);

                    tg_xxxyyy_yyzzzz_0[j] = pb_x * tg_xxyyy_yyzzzz_0[j] + fr * tg_xxyyy_yyzzzz_1[j] + fl1_fx * (tg_xyyy_yyzzzz_0[j] - tg_xyyy_yyzzzz_1[j] * fl1_fza);

                    tg_xxxyyy_yzzzzz_0[j] = pb_x * tg_xxyyy_yzzzzz_0[j] + fr * tg_xxyyy_yzzzzz_1[j] + fl1_fx * (tg_xyyy_yzzzzz_0[j] - tg_xyyy_yzzzzz_1[j] * fl1_fza);

                    tg_xxxyyy_zzzzzz_0[j] = pb_x * tg_xxyyy_zzzzzz_0[j] + fr * tg_xxyyy_zzzzzz_1[j] + fl1_fx * (tg_xyyy_zzzzzz_0[j] - tg_xyyy_zzzzzz_1[j] * fl1_fza);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSISI_196_294(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (196,294)

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
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_6_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_6_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_5_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
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

                auto tg_xxyyz_xxxxxx_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 196); 

                auto tg_xxyyz_xxxxxy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 197); 

                auto tg_xxyyz_xxxxxz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 198); 

                auto tg_xxyyz_xxxxyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 199); 

                auto tg_xxyyz_xxxxyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 200); 

                auto tg_xxyyz_xxxxzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 201); 

                auto tg_xxyyz_xxxyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 202); 

                auto tg_xxyyz_xxxyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 203); 

                auto tg_xxyyz_xxxyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 204); 

                auto tg_xxyyz_xxxzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 205); 

                auto tg_xxyyz_xxyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 206); 

                auto tg_xxyyz_xxyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 207); 

                auto tg_xxyyz_xxyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 208); 

                auto tg_xxyyz_xxyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 209); 

                auto tg_xxyyz_xxzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 210); 

                auto tg_xxyyz_xyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 211); 

                auto tg_xxyyz_xyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 212); 

                auto tg_xxyyz_xyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 213); 

                auto tg_xxyyz_xyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 214); 

                auto tg_xxyyz_xyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 215); 

                auto tg_xxyyz_xzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 216); 

                auto tg_xxyyz_yyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 217); 

                auto tg_xxyyz_yyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 218); 

                auto tg_xxyyz_yyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 219); 

                auto tg_xxyyz_yyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 220); 

                auto tg_xxyyz_yyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 221); 

                auto tg_xxyyz_yzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 222); 

                auto tg_xxyyz_zzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 223); 

                auto tg_xxyzz_xxxxxx_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 224); 

                auto tg_xxyzz_xxxxxy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 225); 

                auto tg_xxyzz_xxxxxz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 226); 

                auto tg_xxyzz_xxxxyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 227); 

                auto tg_xxyzz_xxxxyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 228); 

                auto tg_xxyzz_xxxxzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 229); 

                auto tg_xxyzz_xxxyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 230); 

                auto tg_xxyzz_xxxyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 231); 

                auto tg_xxyzz_xxxyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 232); 

                auto tg_xxyzz_xxxzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 233); 

                auto tg_xxyzz_xxyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 234); 

                auto tg_xxyzz_xxyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 235); 

                auto tg_xxyzz_xxyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 236); 

                auto tg_xxyzz_xxyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 237); 

                auto tg_xxyzz_xxzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 238); 

                auto tg_xxyzz_xyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 239); 

                auto tg_xxyzz_xyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 240); 

                auto tg_xxyzz_xyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 241); 

                auto tg_xxyzz_xyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 242); 

                auto tg_xxyzz_xyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 243); 

                auto tg_xxyzz_xzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 244); 

                auto tg_xxyzz_yyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 245); 

                auto tg_xxyzz_yyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 246); 

                auto tg_xxyzz_yyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 247); 

                auto tg_xxyzz_yyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 248); 

                auto tg_xxyzz_yyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 249); 

                auto tg_xxyzz_yzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 250); 

                auto tg_xxyzz_zzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 251); 

                auto tg_xxzzz_xxxxxx_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 252); 

                auto tg_xxzzz_xxxxxy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 253); 

                auto tg_xxzzz_xxxxxz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 254); 

                auto tg_xxzzz_xxxxyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 255); 

                auto tg_xxzzz_xxxxyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 256); 

                auto tg_xxzzz_xxxxzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 257); 

                auto tg_xxzzz_xxxyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 258); 

                auto tg_xxzzz_xxxyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 259); 

                auto tg_xxzzz_xxxyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 260); 

                auto tg_xxzzz_xxxzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 261); 

                auto tg_xxzzz_xxyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 262); 

                auto tg_xxzzz_xxyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 263); 

                auto tg_xxzzz_xxyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 264); 

                auto tg_xxzzz_xxyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 265); 

                auto tg_xxzzz_xxzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 266); 

                auto tg_xxzzz_xyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 267); 

                auto tg_xxzzz_xyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 268); 

                auto tg_xxzzz_xyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 269); 

                auto tg_xxzzz_xyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 270); 

                auto tg_xxzzz_xyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 271); 

                auto tg_xxzzz_xzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 272); 

                auto tg_xxzzz_yyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 273); 

                auto tg_xxzzz_yyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 274); 

                auto tg_xxzzz_yyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 275); 

                auto tg_xxzzz_yyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 276); 

                auto tg_xxzzz_yyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 277); 

                auto tg_xxzzz_yzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 278); 

                auto tg_xxzzz_zzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 279); 

                auto tg_xyyyy_xxxxxx_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 280); 

                auto tg_xyyyy_xxxxxy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 281); 

                auto tg_xyyyy_xxxxxz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 282); 

                auto tg_xyyyy_xxxxyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 283); 

                auto tg_xyyyy_xxxxyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 284); 

                auto tg_xyyyy_xxxxzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 285); 

                auto tg_xyyyy_xxxyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 286); 

                auto tg_xyyyy_xxxyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 287); 

                auto tg_xyyyy_xxxyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 288); 

                auto tg_xyyyy_xxxzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 289); 

                auto tg_xyyyy_xxyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 290); 

                auto tg_xyyyy_xxyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 291); 

                auto tg_xyyyy_xxyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 292); 

                auto tg_xyyyy_xxyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 293); 

                auto tg_xxyyz_xxxxxx_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 196); 

                auto tg_xxyyz_xxxxxy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 197); 

                auto tg_xxyyz_xxxxxz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 198); 

                auto tg_xxyyz_xxxxyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 199); 

                auto tg_xxyyz_xxxxyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 200); 

                auto tg_xxyyz_xxxxzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 201); 

                auto tg_xxyyz_xxxyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 202); 

                auto tg_xxyyz_xxxyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 203); 

                auto tg_xxyyz_xxxyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 204); 

                auto tg_xxyyz_xxxzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 205); 

                auto tg_xxyyz_xxyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 206); 

                auto tg_xxyyz_xxyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 207); 

                auto tg_xxyyz_xxyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 208); 

                auto tg_xxyyz_xxyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 209); 

                auto tg_xxyyz_xxzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 210); 

                auto tg_xxyyz_xyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 211); 

                auto tg_xxyyz_xyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 212); 

                auto tg_xxyyz_xyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 213); 

                auto tg_xxyyz_xyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 214); 

                auto tg_xxyyz_xyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 215); 

                auto tg_xxyyz_xzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 216); 

                auto tg_xxyyz_yyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 217); 

                auto tg_xxyyz_yyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 218); 

                auto tg_xxyyz_yyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 219); 

                auto tg_xxyyz_yyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 220); 

                auto tg_xxyyz_yyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 221); 

                auto tg_xxyyz_yzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 222); 

                auto tg_xxyyz_zzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 223); 

                auto tg_xxyzz_xxxxxx_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 224); 

                auto tg_xxyzz_xxxxxy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 225); 

                auto tg_xxyzz_xxxxxz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 226); 

                auto tg_xxyzz_xxxxyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 227); 

                auto tg_xxyzz_xxxxyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 228); 

                auto tg_xxyzz_xxxxzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 229); 

                auto tg_xxyzz_xxxyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 230); 

                auto tg_xxyzz_xxxyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 231); 

                auto tg_xxyzz_xxxyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 232); 

                auto tg_xxyzz_xxxzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 233); 

                auto tg_xxyzz_xxyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 234); 

                auto tg_xxyzz_xxyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 235); 

                auto tg_xxyzz_xxyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 236); 

                auto tg_xxyzz_xxyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 237); 

                auto tg_xxyzz_xxzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 238); 

                auto tg_xxyzz_xyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 239); 

                auto tg_xxyzz_xyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 240); 

                auto tg_xxyzz_xyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 241); 

                auto tg_xxyzz_xyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 242); 

                auto tg_xxyzz_xyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 243); 

                auto tg_xxyzz_xzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 244); 

                auto tg_xxyzz_yyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 245); 

                auto tg_xxyzz_yyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 246); 

                auto tg_xxyzz_yyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 247); 

                auto tg_xxyzz_yyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 248); 

                auto tg_xxyzz_yyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 249); 

                auto tg_xxyzz_yzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 250); 

                auto tg_xxyzz_zzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 251); 

                auto tg_xxzzz_xxxxxx_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 252); 

                auto tg_xxzzz_xxxxxy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 253); 

                auto tg_xxzzz_xxxxxz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 254); 

                auto tg_xxzzz_xxxxyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 255); 

                auto tg_xxzzz_xxxxyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 256); 

                auto tg_xxzzz_xxxxzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 257); 

                auto tg_xxzzz_xxxyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 258); 

                auto tg_xxzzz_xxxyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 259); 

                auto tg_xxzzz_xxxyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 260); 

                auto tg_xxzzz_xxxzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 261); 

                auto tg_xxzzz_xxyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 262); 

                auto tg_xxzzz_xxyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 263); 

                auto tg_xxzzz_xxyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 264); 

                auto tg_xxzzz_xxyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 265); 

                auto tg_xxzzz_xxzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 266); 

                auto tg_xxzzz_xyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 267); 

                auto tg_xxzzz_xyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 268); 

                auto tg_xxzzz_xyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 269); 

                auto tg_xxzzz_xyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 270); 

                auto tg_xxzzz_xyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 271); 

                auto tg_xxzzz_xzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 272); 

                auto tg_xxzzz_yyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 273); 

                auto tg_xxzzz_yyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 274); 

                auto tg_xxzzz_yyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 275); 

                auto tg_xxzzz_yyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 276); 

                auto tg_xxzzz_yyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 277); 

                auto tg_xxzzz_yzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 278); 

                auto tg_xxzzz_zzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 279); 

                auto tg_xyyyy_xxxxxx_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 280); 

                auto tg_xyyyy_xxxxxy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 281); 

                auto tg_xyyyy_xxxxxz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 282); 

                auto tg_xyyyy_xxxxyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 283); 

                auto tg_xyyyy_xxxxyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 284); 

                auto tg_xyyyy_xxxxzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 285); 

                auto tg_xyyyy_xxxyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 286); 

                auto tg_xyyyy_xxxyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 287); 

                auto tg_xyyyy_xxxyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 288); 

                auto tg_xyyyy_xxxzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 289); 

                auto tg_xyyyy_xxyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 290); 

                auto tg_xyyyy_xxyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 291); 

                auto tg_xyyyy_xxyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 292); 

                auto tg_xyyyy_xxyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 293); 

                auto tg_xyyz_xxxxxx_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 196); 

                auto tg_xyyz_xxxxxy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 197); 

                auto tg_xyyz_xxxxxz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 198); 

                auto tg_xyyz_xxxxyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 199); 

                auto tg_xyyz_xxxxyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 200); 

                auto tg_xyyz_xxxxzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 201); 

                auto tg_xyyz_xxxyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 202); 

                auto tg_xyyz_xxxyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 203); 

                auto tg_xyyz_xxxyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 204); 

                auto tg_xyyz_xxxzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 205); 

                auto tg_xyyz_xxyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 206); 

                auto tg_xyyz_xxyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 207); 

                auto tg_xyyz_xxyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 208); 

                auto tg_xyyz_xxyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 209); 

                auto tg_xyyz_xxzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 210); 

                auto tg_xyyz_xyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 211); 

                auto tg_xyyz_xyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 212); 

                auto tg_xyyz_xyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 213); 

                auto tg_xyyz_xyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 214); 

                auto tg_xyyz_xyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 215); 

                auto tg_xyyz_xzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 216); 

                auto tg_xyyz_yyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 217); 

                auto tg_xyyz_yyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 218); 

                auto tg_xyyz_yyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 219); 

                auto tg_xyyz_yyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 220); 

                auto tg_xyyz_yyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 221); 

                auto tg_xyyz_yzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 222); 

                auto tg_xyyz_zzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 223); 

                auto tg_xyzz_xxxxxx_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 224); 

                auto tg_xyzz_xxxxxy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 225); 

                auto tg_xyzz_xxxxxz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 226); 

                auto tg_xyzz_xxxxyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 227); 

                auto tg_xyzz_xxxxyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 228); 

                auto tg_xyzz_xxxxzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 229); 

                auto tg_xyzz_xxxyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 230); 

                auto tg_xyzz_xxxyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 231); 

                auto tg_xyzz_xxxyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 232); 

                auto tg_xyzz_xxxzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 233); 

                auto tg_xyzz_xxyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 234); 

                auto tg_xyzz_xxyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 235); 

                auto tg_xyzz_xxyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 236); 

                auto tg_xyzz_xxyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 237); 

                auto tg_xyzz_xxzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 238); 

                auto tg_xyzz_xyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 239); 

                auto tg_xyzz_xyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 240); 

                auto tg_xyzz_xyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 241); 

                auto tg_xyzz_xyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 242); 

                auto tg_xyzz_xyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 243); 

                auto tg_xyzz_xzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 244); 

                auto tg_xyzz_yyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 245); 

                auto tg_xyzz_yyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 246); 

                auto tg_xyzz_yyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 247); 

                auto tg_xyzz_yyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 248); 

                auto tg_xyzz_yyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 249); 

                auto tg_xyzz_yzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 250); 

                auto tg_xyzz_zzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 251); 

                auto tg_xzzz_xxxxxx_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 252); 

                auto tg_xzzz_xxxxxy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 253); 

                auto tg_xzzz_xxxxxz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 254); 

                auto tg_xzzz_xxxxyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 255); 

                auto tg_xzzz_xxxxyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 256); 

                auto tg_xzzz_xxxxzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 257); 

                auto tg_xzzz_xxxyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 258); 

                auto tg_xzzz_xxxyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 259); 

                auto tg_xzzz_xxxyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 260); 

                auto tg_xzzz_xxxzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 261); 

                auto tg_xzzz_xxyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 262); 

                auto tg_xzzz_xxyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 263); 

                auto tg_xzzz_xxyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 264); 

                auto tg_xzzz_xxyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 265); 

                auto tg_xzzz_xxzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 266); 

                auto tg_xzzz_xyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 267); 

                auto tg_xzzz_xyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 268); 

                auto tg_xzzz_xyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 269); 

                auto tg_xzzz_xyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 270); 

                auto tg_xzzz_xyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 271); 

                auto tg_xzzz_xzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 272); 

                auto tg_xzzz_yyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 273); 

                auto tg_xzzz_yyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 274); 

                auto tg_xzzz_yyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 275); 

                auto tg_xzzz_yyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 276); 

                auto tg_xzzz_yyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 277); 

                auto tg_xzzz_yzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 278); 

                auto tg_xzzz_zzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 279); 

                auto tg_yyyy_xxxxxx_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 280); 

                auto tg_yyyy_xxxxxy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 281); 

                auto tg_yyyy_xxxxxz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 282); 

                auto tg_yyyy_xxxxyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 283); 

                auto tg_yyyy_xxxxyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 284); 

                auto tg_yyyy_xxxxzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 285); 

                auto tg_yyyy_xxxyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 286); 

                auto tg_yyyy_xxxyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 287); 

                auto tg_yyyy_xxxyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 288); 

                auto tg_yyyy_xxxzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 289); 

                auto tg_yyyy_xxyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 290); 

                auto tg_yyyy_xxyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 291); 

                auto tg_yyyy_xxyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 292); 

                auto tg_yyyy_xxyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 293); 

                auto tg_xyyz_xxxxxx_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 196); 

                auto tg_xyyz_xxxxxy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 197); 

                auto tg_xyyz_xxxxxz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 198); 

                auto tg_xyyz_xxxxyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 199); 

                auto tg_xyyz_xxxxyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 200); 

                auto tg_xyyz_xxxxzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 201); 

                auto tg_xyyz_xxxyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 202); 

                auto tg_xyyz_xxxyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 203); 

                auto tg_xyyz_xxxyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 204); 

                auto tg_xyyz_xxxzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 205); 

                auto tg_xyyz_xxyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 206); 

                auto tg_xyyz_xxyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 207); 

                auto tg_xyyz_xxyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 208); 

                auto tg_xyyz_xxyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 209); 

                auto tg_xyyz_xxzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 210); 

                auto tg_xyyz_xyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 211); 

                auto tg_xyyz_xyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 212); 

                auto tg_xyyz_xyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 213); 

                auto tg_xyyz_xyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 214); 

                auto tg_xyyz_xyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 215); 

                auto tg_xyyz_xzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 216); 

                auto tg_xyyz_yyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 217); 

                auto tg_xyyz_yyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 218); 

                auto tg_xyyz_yyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 219); 

                auto tg_xyyz_yyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 220); 

                auto tg_xyyz_yyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 221); 

                auto tg_xyyz_yzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 222); 

                auto tg_xyyz_zzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 223); 

                auto tg_xyzz_xxxxxx_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 224); 

                auto tg_xyzz_xxxxxy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 225); 

                auto tg_xyzz_xxxxxz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 226); 

                auto tg_xyzz_xxxxyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 227); 

                auto tg_xyzz_xxxxyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 228); 

                auto tg_xyzz_xxxxzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 229); 

                auto tg_xyzz_xxxyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 230); 

                auto tg_xyzz_xxxyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 231); 

                auto tg_xyzz_xxxyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 232); 

                auto tg_xyzz_xxxzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 233); 

                auto tg_xyzz_xxyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 234); 

                auto tg_xyzz_xxyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 235); 

                auto tg_xyzz_xxyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 236); 

                auto tg_xyzz_xxyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 237); 

                auto tg_xyzz_xxzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 238); 

                auto tg_xyzz_xyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 239); 

                auto tg_xyzz_xyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 240); 

                auto tg_xyzz_xyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 241); 

                auto tg_xyzz_xyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 242); 

                auto tg_xyzz_xyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 243); 

                auto tg_xyzz_xzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 244); 

                auto tg_xyzz_yyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 245); 

                auto tg_xyzz_yyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 246); 

                auto tg_xyzz_yyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 247); 

                auto tg_xyzz_yyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 248); 

                auto tg_xyzz_yyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 249); 

                auto tg_xyzz_yzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 250); 

                auto tg_xyzz_zzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 251); 

                auto tg_xzzz_xxxxxx_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 252); 

                auto tg_xzzz_xxxxxy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 253); 

                auto tg_xzzz_xxxxxz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 254); 

                auto tg_xzzz_xxxxyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 255); 

                auto tg_xzzz_xxxxyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 256); 

                auto tg_xzzz_xxxxzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 257); 

                auto tg_xzzz_xxxyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 258); 

                auto tg_xzzz_xxxyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 259); 

                auto tg_xzzz_xxxyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 260); 

                auto tg_xzzz_xxxzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 261); 

                auto tg_xzzz_xxyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 262); 

                auto tg_xzzz_xxyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 263); 

                auto tg_xzzz_xxyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 264); 

                auto tg_xzzz_xxyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 265); 

                auto tg_xzzz_xxzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 266); 

                auto tg_xzzz_xyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 267); 

                auto tg_xzzz_xyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 268); 

                auto tg_xzzz_xyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 269); 

                auto tg_xzzz_xyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 270); 

                auto tg_xzzz_xyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 271); 

                auto tg_xzzz_xzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 272); 

                auto tg_xzzz_yyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 273); 

                auto tg_xzzz_yyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 274); 

                auto tg_xzzz_yyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 275); 

                auto tg_xzzz_yyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 276); 

                auto tg_xzzz_yyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 277); 

                auto tg_xzzz_yzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 278); 

                auto tg_xzzz_zzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 279); 

                auto tg_yyyy_xxxxxx_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 280); 

                auto tg_yyyy_xxxxxy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 281); 

                auto tg_yyyy_xxxxxz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 282); 

                auto tg_yyyy_xxxxyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 283); 

                auto tg_yyyy_xxxxyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 284); 

                auto tg_yyyy_xxxxzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 285); 

                auto tg_yyyy_xxxyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 286); 

                auto tg_yyyy_xxxyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 287); 

                auto tg_yyyy_xxxyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 288); 

                auto tg_yyyy_xxxzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 289); 

                auto tg_yyyy_xxyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 290); 

                auto tg_yyyy_xxyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 291); 

                auto tg_yyyy_xxyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 292); 

                auto tg_yyyy_xxyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 293); 

                auto tg_xxyyz_xxxxx_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 147); 

                auto tg_xxyyz_xxxxy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 148); 

                auto tg_xxyyz_xxxxz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 149); 

                auto tg_xxyyz_xxxyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 150); 

                auto tg_xxyyz_xxxyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 151); 

                auto tg_xxyyz_xxxzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 152); 

                auto tg_xxyyz_xxyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 153); 

                auto tg_xxyyz_xxyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 154); 

                auto tg_xxyyz_xxyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 155); 

                auto tg_xxyyz_xxzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 156); 

                auto tg_xxyyz_xyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 157); 

                auto tg_xxyyz_xyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 158); 

                auto tg_xxyyz_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 159); 

                auto tg_xxyyz_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 160); 

                auto tg_xxyyz_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 161); 

                auto tg_xxyyz_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 162); 

                auto tg_xxyyz_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 163); 

                auto tg_xxyyz_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 164); 

                auto tg_xxyyz_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 165); 

                auto tg_xxyyz_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 166); 

                auto tg_xxyyz_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 167); 

                auto tg_xxyzz_xxxxx_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 168); 

                auto tg_xxyzz_xxxxy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 169); 

                auto tg_xxyzz_xxxxz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 170); 

                auto tg_xxyzz_xxxyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 171); 

                auto tg_xxyzz_xxxyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 172); 

                auto tg_xxyzz_xxxzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 173); 

                auto tg_xxyzz_xxyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 174); 

                auto tg_xxyzz_xxyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 175); 

                auto tg_xxyzz_xxyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 176); 

                auto tg_xxyzz_xxzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 177); 

                auto tg_xxyzz_xyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 178); 

                auto tg_xxyzz_xyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 179); 

                auto tg_xxyzz_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 180); 

                auto tg_xxyzz_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 181); 

                auto tg_xxyzz_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 182); 

                auto tg_xxyzz_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 183); 

                auto tg_xxyzz_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 184); 

                auto tg_xxyzz_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 185); 

                auto tg_xxyzz_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 186); 

                auto tg_xxyzz_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 187); 

                auto tg_xxyzz_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 188); 

                auto tg_xxzzz_xxxxx_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 189); 

                auto tg_xxzzz_xxxxy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 190); 

                auto tg_xxzzz_xxxxz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 191); 

                auto tg_xxzzz_xxxyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 192); 

                auto tg_xxzzz_xxxyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 193); 

                auto tg_xxzzz_xxxzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 194); 

                auto tg_xxzzz_xxyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 195); 

                auto tg_xxzzz_xxyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 196); 

                auto tg_xxzzz_xxyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 197); 

                auto tg_xxzzz_xxzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 198); 

                auto tg_xxzzz_xyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 199); 

                auto tg_xxzzz_xyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 200); 

                auto tg_xxzzz_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 201); 

                auto tg_xxzzz_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 202); 

                auto tg_xxzzz_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 203); 

                auto tg_xxzzz_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 204); 

                auto tg_xxzzz_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 205); 

                auto tg_xxzzz_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 206); 

                auto tg_xxzzz_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 207); 

                auto tg_xxzzz_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 208); 

                auto tg_xxzzz_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 209); 

                auto tg_xyyyy_xxxxx_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 210); 

                auto tg_xyyyy_xxxxy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 211); 

                auto tg_xyyyy_xxxxz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 212); 

                auto tg_xyyyy_xxxyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 213); 

                auto tg_xyyyy_xxxyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 214); 

                auto tg_xyyyy_xxxzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 215); 

                auto tg_xyyyy_xxyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 216); 

                auto tg_xyyyy_xxyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 217); 

                auto tg_xyyyy_xxyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 218); 

                auto tg_xyyyy_xxzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 219); 

                auto tg_xyyyy_xyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 220); 

                auto tg_xyyyy_xyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 221); 

                auto tg_xyyyy_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 222); 

                auto tg_xyyyy_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 223); 

                // set up pointers to integrals

                auto tg_xxxyyz_xxxxxx_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 196); 

                auto tg_xxxyyz_xxxxxy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 197); 

                auto tg_xxxyyz_xxxxxz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 198); 

                auto tg_xxxyyz_xxxxyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 199); 

                auto tg_xxxyyz_xxxxyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 200); 

                auto tg_xxxyyz_xxxxzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 201); 

                auto tg_xxxyyz_xxxyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 202); 

                auto tg_xxxyyz_xxxyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 203); 

                auto tg_xxxyyz_xxxyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 204); 

                auto tg_xxxyyz_xxxzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 205); 

                auto tg_xxxyyz_xxyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 206); 

                auto tg_xxxyyz_xxyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 207); 

                auto tg_xxxyyz_xxyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 208); 

                auto tg_xxxyyz_xxyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 209); 

                auto tg_xxxyyz_xxzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 210); 

                auto tg_xxxyyz_xyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 211); 

                auto tg_xxxyyz_xyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 212); 

                auto tg_xxxyyz_xyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 213); 

                auto tg_xxxyyz_xyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 214); 

                auto tg_xxxyyz_xyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 215); 

                auto tg_xxxyyz_xzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 216); 

                auto tg_xxxyyz_yyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 217); 

                auto tg_xxxyyz_yyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 218); 

                auto tg_xxxyyz_yyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 219); 

                auto tg_xxxyyz_yyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 220); 

                auto tg_xxxyyz_yyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 221); 

                auto tg_xxxyyz_yzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 222); 

                auto tg_xxxyyz_zzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 223); 

                auto tg_xxxyzz_xxxxxx_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 224); 

                auto tg_xxxyzz_xxxxxy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 225); 

                auto tg_xxxyzz_xxxxxz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 226); 

                auto tg_xxxyzz_xxxxyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 227); 

                auto tg_xxxyzz_xxxxyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 228); 

                auto tg_xxxyzz_xxxxzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 229); 

                auto tg_xxxyzz_xxxyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 230); 

                auto tg_xxxyzz_xxxyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 231); 

                auto tg_xxxyzz_xxxyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 232); 

                auto tg_xxxyzz_xxxzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 233); 

                auto tg_xxxyzz_xxyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 234); 

                auto tg_xxxyzz_xxyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 235); 

                auto tg_xxxyzz_xxyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 236); 

                auto tg_xxxyzz_xxyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 237); 

                auto tg_xxxyzz_xxzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 238); 

                auto tg_xxxyzz_xyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 239); 

                auto tg_xxxyzz_xyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 240); 

                auto tg_xxxyzz_xyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 241); 

                auto tg_xxxyzz_xyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 242); 

                auto tg_xxxyzz_xyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 243); 

                auto tg_xxxyzz_xzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 244); 

                auto tg_xxxyzz_yyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 245); 

                auto tg_xxxyzz_yyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 246); 

                auto tg_xxxyzz_yyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 247); 

                auto tg_xxxyzz_yyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 248); 

                auto tg_xxxyzz_yyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 249); 

                auto tg_xxxyzz_yzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 250); 

                auto tg_xxxyzz_zzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 251); 

                auto tg_xxxzzz_xxxxxx_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 252); 

                auto tg_xxxzzz_xxxxxy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 253); 

                auto tg_xxxzzz_xxxxxz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 254); 

                auto tg_xxxzzz_xxxxyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 255); 

                auto tg_xxxzzz_xxxxyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 256); 

                auto tg_xxxzzz_xxxxzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 257); 

                auto tg_xxxzzz_xxxyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 258); 

                auto tg_xxxzzz_xxxyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 259); 

                auto tg_xxxzzz_xxxyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 260); 

                auto tg_xxxzzz_xxxzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 261); 

                auto tg_xxxzzz_xxyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 262); 

                auto tg_xxxzzz_xxyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 263); 

                auto tg_xxxzzz_xxyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 264); 

                auto tg_xxxzzz_xxyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 265); 

                auto tg_xxxzzz_xxzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 266); 

                auto tg_xxxzzz_xyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 267); 

                auto tg_xxxzzz_xyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 268); 

                auto tg_xxxzzz_xyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 269); 

                auto tg_xxxzzz_xyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 270); 

                auto tg_xxxzzz_xyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 271); 

                auto tg_xxxzzz_xzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 272); 

                auto tg_xxxzzz_yyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 273); 

                auto tg_xxxzzz_yyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 274); 

                auto tg_xxxzzz_yyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 275); 

                auto tg_xxxzzz_yyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 276); 

                auto tg_xxxzzz_yyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 277); 

                auto tg_xxxzzz_yzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 278); 

                auto tg_xxxzzz_zzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 279); 

                auto tg_xxyyyy_xxxxxx_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 280); 

                auto tg_xxyyyy_xxxxxy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 281); 

                auto tg_xxyyyy_xxxxxz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 282); 

                auto tg_xxyyyy_xxxxyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 283); 

                auto tg_xxyyyy_xxxxyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 284); 

                auto tg_xxyyyy_xxxxzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 285); 

                auto tg_xxyyyy_xxxyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 286); 

                auto tg_xxyyyy_xxxyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 287); 

                auto tg_xxyyyy_xxxyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 288); 

                auto tg_xxyyyy_xxxzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 289); 

                auto tg_xxyyyy_xxyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 290); 

                auto tg_xxyyyy_xxyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 291); 

                auto tg_xxyyyy_xxyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 292); 

                auto tg_xxyyyy_xxyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 293); 

                // Batch of Integrals (196,294)

                #pragma omp simd aligned(fxn, fza, tg_xxxyyz_xxxxxx_0, tg_xxxyyz_xxxxxy_0, tg_xxxyyz_xxxxxz_0, \
                                         tg_xxxyyz_xxxxyy_0, tg_xxxyyz_xxxxyz_0, tg_xxxyyz_xxxxzz_0, tg_xxxyyz_xxxyyy_0, \
                                         tg_xxxyyz_xxxyyz_0, tg_xxxyyz_xxxyzz_0, tg_xxxyyz_xxxzzz_0, tg_xxxyyz_xxyyyy_0, \
                                         tg_xxxyyz_xxyyyz_0, tg_xxxyyz_xxyyzz_0, tg_xxxyyz_xxyzzz_0, tg_xxxyyz_xxzzzz_0, \
                                         tg_xxxyyz_xyyyyy_0, tg_xxxyyz_xyyyyz_0, tg_xxxyyz_xyyyzz_0, tg_xxxyyz_xyyzzz_0, \
                                         tg_xxxyyz_xyzzzz_0, tg_xxxyyz_xzzzzz_0, tg_xxxyyz_yyyyyy_0, tg_xxxyyz_yyyyyz_0, \
                                         tg_xxxyyz_yyyyzz_0, tg_xxxyyz_yyyzzz_0, tg_xxxyyz_yyzzzz_0, tg_xxxyyz_yzzzzz_0, \
                                         tg_xxxyyz_zzzzzz_0, tg_xxxyzz_xxxxxx_0, tg_xxxyzz_xxxxxy_0, tg_xxxyzz_xxxxxz_0, \
                                         tg_xxxyzz_xxxxyy_0, tg_xxxyzz_xxxxyz_0, tg_xxxyzz_xxxxzz_0, tg_xxxyzz_xxxyyy_0, \
                                         tg_xxxyzz_xxxyyz_0, tg_xxxyzz_xxxyzz_0, tg_xxxyzz_xxxzzz_0, tg_xxxyzz_xxyyyy_0, \
                                         tg_xxxyzz_xxyyyz_0, tg_xxxyzz_xxyyzz_0, tg_xxxyzz_xxyzzz_0, tg_xxxyzz_xxzzzz_0, \
                                         tg_xxxyzz_xyyyyy_0, tg_xxxyzz_xyyyyz_0, tg_xxxyzz_xyyyzz_0, tg_xxxyzz_xyyzzz_0, \
                                         tg_xxxyzz_xyzzzz_0, tg_xxxyzz_xzzzzz_0, tg_xxxyzz_yyyyyy_0, tg_xxxyzz_yyyyyz_0, \
                                         tg_xxxyzz_yyyyzz_0, tg_xxxyzz_yyyzzz_0, tg_xxxyzz_yyzzzz_0, tg_xxxyzz_yzzzzz_0, \
                                         tg_xxxyzz_zzzzzz_0, tg_xxxzzz_xxxxxx_0, tg_xxxzzz_xxxxxy_0, tg_xxxzzz_xxxxxz_0, \
                                         tg_xxxzzz_xxxxyy_0, tg_xxxzzz_xxxxyz_0, tg_xxxzzz_xxxxzz_0, tg_xxxzzz_xxxyyy_0, \
                                         tg_xxxzzz_xxxyyz_0, tg_xxxzzz_xxxyzz_0, tg_xxxzzz_xxxzzz_0, tg_xxxzzz_xxyyyy_0, \
                                         tg_xxxzzz_xxyyyz_0, tg_xxxzzz_xxyyzz_0, tg_xxxzzz_xxyzzz_0, tg_xxxzzz_xxzzzz_0, \
                                         tg_xxxzzz_xyyyyy_0, tg_xxxzzz_xyyyyz_0, tg_xxxzzz_xyyyzz_0, tg_xxxzzz_xyyzzz_0, \
                                         tg_xxxzzz_xyzzzz_0, tg_xxxzzz_xzzzzz_0, tg_xxxzzz_yyyyyy_0, tg_xxxzzz_yyyyyz_0, \
                                         tg_xxxzzz_yyyyzz_0, tg_xxxzzz_yyyzzz_0, tg_xxxzzz_yyzzzz_0, tg_xxxzzz_yzzzzz_0, \
                                         tg_xxxzzz_zzzzzz_0, tg_xxyyyy_xxxxxx_0, tg_xxyyyy_xxxxxy_0, tg_xxyyyy_xxxxxz_0, \
                                         tg_xxyyyy_xxxxyy_0, tg_xxyyyy_xxxxyz_0, tg_xxyyyy_xxxxzz_0, tg_xxyyyy_xxxyyy_0, \
                                         tg_xxyyyy_xxxyyz_0, tg_xxyyyy_xxxyzz_0, tg_xxyyyy_xxxzzz_0, tg_xxyyyy_xxyyyy_0, \
                                         tg_xxyyyy_xxyyyz_0, tg_xxyyyy_xxyyzz_0, tg_xxyyyy_xxyzzz_0, tg_xxyyz_xxxxx_1, \
                                         tg_xxyyz_xxxxxx_0, tg_xxyyz_xxxxxx_1, tg_xxyyz_xxxxxy_0, tg_xxyyz_xxxxxy_1, \
                                         tg_xxyyz_xxxxxz_0, tg_xxyyz_xxxxxz_1, tg_xxyyz_xxxxy_1, tg_xxyyz_xxxxyy_0, \
                                         tg_xxyyz_xxxxyy_1, tg_xxyyz_xxxxyz_0, tg_xxyyz_xxxxyz_1, tg_xxyyz_xxxxz_1, \
                                         tg_xxyyz_xxxxzz_0, tg_xxyyz_xxxxzz_1, tg_xxyyz_xxxyy_1, tg_xxyyz_xxxyyy_0, \
                                         tg_xxyyz_xxxyyy_1, tg_xxyyz_xxxyyz_0, tg_xxyyz_xxxyyz_1, tg_xxyyz_xxxyz_1, \
                                         tg_xxyyz_xxxyzz_0, tg_xxyyz_xxxyzz_1, tg_xxyyz_xxxzz_1, tg_xxyyz_xxxzzz_0, \
                                         tg_xxyyz_xxxzzz_1, tg_xxyyz_xxyyy_1, tg_xxyyz_xxyyyy_0, tg_xxyyz_xxyyyy_1, \
                                         tg_xxyyz_xxyyyz_0, tg_xxyyz_xxyyyz_1, tg_xxyyz_xxyyz_1, tg_xxyyz_xxyyzz_0, \
                                         tg_xxyyz_xxyyzz_1, tg_xxyyz_xxyzz_1, tg_xxyyz_xxyzzz_0, tg_xxyyz_xxyzzz_1, \
                                         tg_xxyyz_xxzzz_1, tg_xxyyz_xxzzzz_0, tg_xxyyz_xxzzzz_1, tg_xxyyz_xyyyy_1, \
                                         tg_xxyyz_xyyyyy_0, tg_xxyyz_xyyyyy_1, tg_xxyyz_xyyyyz_0, tg_xxyyz_xyyyyz_1, \
                                         tg_xxyyz_xyyyz_1, tg_xxyyz_xyyyzz_0, tg_xxyyz_xyyyzz_1, tg_xxyyz_xyyzz_1, \
                                         tg_xxyyz_xyyzzz_0, tg_xxyyz_xyyzzz_1, tg_xxyyz_xyzzz_1, tg_xxyyz_xyzzzz_0, \
                                         tg_xxyyz_xyzzzz_1, tg_xxyyz_xzzzz_1, tg_xxyyz_xzzzzz_0, tg_xxyyz_xzzzzz_1, \
                                         tg_xxyyz_yyyyy_1, tg_xxyyz_yyyyyy_0, tg_xxyyz_yyyyyy_1, tg_xxyyz_yyyyyz_0, \
                                         tg_xxyyz_yyyyyz_1, tg_xxyyz_yyyyz_1, tg_xxyyz_yyyyzz_0, tg_xxyyz_yyyyzz_1, \
                                         tg_xxyyz_yyyzz_1, tg_xxyyz_yyyzzz_0, tg_xxyyz_yyyzzz_1, tg_xxyyz_yyzzz_1, \
                                         tg_xxyyz_yyzzzz_0, tg_xxyyz_yyzzzz_1, tg_xxyyz_yzzzz_1, tg_xxyyz_yzzzzz_0, \
                                         tg_xxyyz_yzzzzz_1, tg_xxyyz_zzzzz_1, tg_xxyyz_zzzzzz_0, tg_xxyyz_zzzzzz_1, \
                                         tg_xxyzz_xxxxx_1, tg_xxyzz_xxxxxx_0, tg_xxyzz_xxxxxx_1, tg_xxyzz_xxxxxy_0, \
                                         tg_xxyzz_xxxxxy_1, tg_xxyzz_xxxxxz_0, tg_xxyzz_xxxxxz_1, tg_xxyzz_xxxxy_1, \
                                         tg_xxyzz_xxxxyy_0, tg_xxyzz_xxxxyy_1, tg_xxyzz_xxxxyz_0, tg_xxyzz_xxxxyz_1, \
                                         tg_xxyzz_xxxxz_1, tg_xxyzz_xxxxzz_0, tg_xxyzz_xxxxzz_1, tg_xxyzz_xxxyy_1, \
                                         tg_xxyzz_xxxyyy_0, tg_xxyzz_xxxyyy_1, tg_xxyzz_xxxyyz_0, tg_xxyzz_xxxyyz_1, \
                                         tg_xxyzz_xxxyz_1, tg_xxyzz_xxxyzz_0, tg_xxyzz_xxxyzz_1, tg_xxyzz_xxxzz_1, \
                                         tg_xxyzz_xxxzzz_0, tg_xxyzz_xxxzzz_1, tg_xxyzz_xxyyy_1, tg_xxyzz_xxyyyy_0, \
                                         tg_xxyzz_xxyyyy_1, tg_xxyzz_xxyyyz_0, tg_xxyzz_xxyyyz_1, tg_xxyzz_xxyyz_1, \
                                         tg_xxyzz_xxyyzz_0, tg_xxyzz_xxyyzz_1, tg_xxyzz_xxyzz_1, tg_xxyzz_xxyzzz_0, \
                                         tg_xxyzz_xxyzzz_1, tg_xxyzz_xxzzz_1, tg_xxyzz_xxzzzz_0, tg_xxyzz_xxzzzz_1, \
                                         tg_xxyzz_xyyyy_1, tg_xxyzz_xyyyyy_0, tg_xxyzz_xyyyyy_1, tg_xxyzz_xyyyyz_0, \
                                         tg_xxyzz_xyyyyz_1, tg_xxyzz_xyyyz_1, tg_xxyzz_xyyyzz_0, tg_xxyzz_xyyyzz_1, \
                                         tg_xxyzz_xyyzz_1, tg_xxyzz_xyyzzz_0, tg_xxyzz_xyyzzz_1, tg_xxyzz_xyzzz_1, \
                                         tg_xxyzz_xyzzzz_0, tg_xxyzz_xyzzzz_1, tg_xxyzz_xzzzz_1, tg_xxyzz_xzzzzz_0, \
                                         tg_xxyzz_xzzzzz_1, tg_xxyzz_yyyyy_1, tg_xxyzz_yyyyyy_0, tg_xxyzz_yyyyyy_1, \
                                         tg_xxyzz_yyyyyz_0, tg_xxyzz_yyyyyz_1, tg_xxyzz_yyyyz_1, tg_xxyzz_yyyyzz_0, \
                                         tg_xxyzz_yyyyzz_1, tg_xxyzz_yyyzz_1, tg_xxyzz_yyyzzz_0, tg_xxyzz_yyyzzz_1, \
                                         tg_xxyzz_yyzzz_1, tg_xxyzz_yyzzzz_0, tg_xxyzz_yyzzzz_1, tg_xxyzz_yzzzz_1, \
                                         tg_xxyzz_yzzzzz_0, tg_xxyzz_yzzzzz_1, tg_xxyzz_zzzzz_1, tg_xxyzz_zzzzzz_0, \
                                         tg_xxyzz_zzzzzz_1, tg_xxzzz_xxxxx_1, tg_xxzzz_xxxxxx_0, tg_xxzzz_xxxxxx_1, \
                                         tg_xxzzz_xxxxxy_0, tg_xxzzz_xxxxxy_1, tg_xxzzz_xxxxxz_0, tg_xxzzz_xxxxxz_1, \
                                         tg_xxzzz_xxxxy_1, tg_xxzzz_xxxxyy_0, tg_xxzzz_xxxxyy_1, tg_xxzzz_xxxxyz_0, \
                                         tg_xxzzz_xxxxyz_1, tg_xxzzz_xxxxz_1, tg_xxzzz_xxxxzz_0, tg_xxzzz_xxxxzz_1, \
                                         tg_xxzzz_xxxyy_1, tg_xxzzz_xxxyyy_0, tg_xxzzz_xxxyyy_1, tg_xxzzz_xxxyyz_0, \
                                         tg_xxzzz_xxxyyz_1, tg_xxzzz_xxxyz_1, tg_xxzzz_xxxyzz_0, tg_xxzzz_xxxyzz_1, \
                                         tg_xxzzz_xxxzz_1, tg_xxzzz_xxxzzz_0, tg_xxzzz_xxxzzz_1, tg_xxzzz_xxyyy_1, \
                                         tg_xxzzz_xxyyyy_0, tg_xxzzz_xxyyyy_1, tg_xxzzz_xxyyyz_0, tg_xxzzz_xxyyyz_1, \
                                         tg_xxzzz_xxyyz_1, tg_xxzzz_xxyyzz_0, tg_xxzzz_xxyyzz_1, tg_xxzzz_xxyzz_1, \
                                         tg_xxzzz_xxyzzz_0, tg_xxzzz_xxyzzz_1, tg_xxzzz_xxzzz_1, tg_xxzzz_xxzzzz_0, \
                                         tg_xxzzz_xxzzzz_1, tg_xxzzz_xyyyy_1, tg_xxzzz_xyyyyy_0, tg_xxzzz_xyyyyy_1, \
                                         tg_xxzzz_xyyyyz_0, tg_xxzzz_xyyyyz_1, tg_xxzzz_xyyyz_1, tg_xxzzz_xyyyzz_0, \
                                         tg_xxzzz_xyyyzz_1, tg_xxzzz_xyyzz_1, tg_xxzzz_xyyzzz_0, tg_xxzzz_xyyzzz_1, \
                                         tg_xxzzz_xyzzz_1, tg_xxzzz_xyzzzz_0, tg_xxzzz_xyzzzz_1, tg_xxzzz_xzzzz_1, \
                                         tg_xxzzz_xzzzzz_0, tg_xxzzz_xzzzzz_1, tg_xxzzz_yyyyy_1, tg_xxzzz_yyyyyy_0, \
                                         tg_xxzzz_yyyyyy_1, tg_xxzzz_yyyyyz_0, tg_xxzzz_yyyyyz_1, tg_xxzzz_yyyyz_1, \
                                         tg_xxzzz_yyyyzz_0, tg_xxzzz_yyyyzz_1, tg_xxzzz_yyyzz_1, tg_xxzzz_yyyzzz_0, \
                                         tg_xxzzz_yyyzzz_1, tg_xxzzz_yyzzz_1, tg_xxzzz_yyzzzz_0, tg_xxzzz_yyzzzz_1, \
                                         tg_xxzzz_yzzzz_1, tg_xxzzz_yzzzzz_0, tg_xxzzz_yzzzzz_1, tg_xxzzz_zzzzz_1, \
                                         tg_xxzzz_zzzzzz_0, tg_xxzzz_zzzzzz_1, tg_xyyyy_xxxxx_1, tg_xyyyy_xxxxxx_0, \
                                         tg_xyyyy_xxxxxx_1, tg_xyyyy_xxxxxy_0, tg_xyyyy_xxxxxy_1, tg_xyyyy_xxxxxz_0, \
                                         tg_xyyyy_xxxxxz_1, tg_xyyyy_xxxxy_1, tg_xyyyy_xxxxyy_0, tg_xyyyy_xxxxyy_1, \
                                         tg_xyyyy_xxxxyz_0, tg_xyyyy_xxxxyz_1, tg_xyyyy_xxxxz_1, tg_xyyyy_xxxxzz_0, \
                                         tg_xyyyy_xxxxzz_1, tg_xyyyy_xxxyy_1, tg_xyyyy_xxxyyy_0, tg_xyyyy_xxxyyy_1, \
                                         tg_xyyyy_xxxyyz_0, tg_xyyyy_xxxyyz_1, tg_xyyyy_xxxyz_1, tg_xyyyy_xxxyzz_0, \
                                         tg_xyyyy_xxxyzz_1, tg_xyyyy_xxxzz_1, tg_xyyyy_xxxzzz_0, tg_xyyyy_xxxzzz_1, \
                                         tg_xyyyy_xxyyy_1, tg_xyyyy_xxyyyy_0, tg_xyyyy_xxyyyy_1, tg_xyyyy_xxyyyz_0, \
                                         tg_xyyyy_xxyyyz_1, tg_xyyyy_xxyyz_1, tg_xyyyy_xxyyzz_0, tg_xyyyy_xxyyzz_1, \
                                         tg_xyyyy_xxyzz_1, tg_xyyyy_xxyzzz_0, tg_xyyyy_xxyzzz_1, tg_xyyyy_xxzzz_1, \
                                         tg_xyyyy_xyyyy_1, tg_xyyyy_xyyyz_1, tg_xyyyy_xyyzz_1, tg_xyyyy_xyzzz_1, \
                                         tg_xyyz_xxxxxx_0, tg_xyyz_xxxxxx_1, tg_xyyz_xxxxxy_0, tg_xyyz_xxxxxy_1, \
                                         tg_xyyz_xxxxxz_0, tg_xyyz_xxxxxz_1, tg_xyyz_xxxxyy_0, tg_xyyz_xxxxyy_1, \
                                         tg_xyyz_xxxxyz_0, tg_xyyz_xxxxyz_1, tg_xyyz_xxxxzz_0, tg_xyyz_xxxxzz_1, \
                                         tg_xyyz_xxxyyy_0, tg_xyyz_xxxyyy_1, tg_xyyz_xxxyyz_0, tg_xyyz_xxxyyz_1, \
                                         tg_xyyz_xxxyzz_0, tg_xyyz_xxxyzz_1, tg_xyyz_xxxzzz_0, tg_xyyz_xxxzzz_1, \
                                         tg_xyyz_xxyyyy_0, tg_xyyz_xxyyyy_1, tg_xyyz_xxyyyz_0, tg_xyyz_xxyyyz_1, \
                                         tg_xyyz_xxyyzz_0, tg_xyyz_xxyyzz_1, tg_xyyz_xxyzzz_0, tg_xyyz_xxyzzz_1, \
                                         tg_xyyz_xxzzzz_0, tg_xyyz_xxzzzz_1, tg_xyyz_xyyyyy_0, tg_xyyz_xyyyyy_1, \
                                         tg_xyyz_xyyyyz_0, tg_xyyz_xyyyyz_1, tg_xyyz_xyyyzz_0, tg_xyyz_xyyyzz_1, \
                                         tg_xyyz_xyyzzz_0, tg_xyyz_xyyzzz_1, tg_xyyz_xyzzzz_0, tg_xyyz_xyzzzz_1, \
                                         tg_xyyz_xzzzzz_0, tg_xyyz_xzzzzz_1, tg_xyyz_yyyyyy_0, tg_xyyz_yyyyyy_1, \
                                         tg_xyyz_yyyyyz_0, tg_xyyz_yyyyyz_1, tg_xyyz_yyyyzz_0, tg_xyyz_yyyyzz_1, \
                                         tg_xyyz_yyyzzz_0, tg_xyyz_yyyzzz_1, tg_xyyz_yyzzzz_0, tg_xyyz_yyzzzz_1, \
                                         tg_xyyz_yzzzzz_0, tg_xyyz_yzzzzz_1, tg_xyyz_zzzzzz_0, tg_xyyz_zzzzzz_1, \
                                         tg_xyzz_xxxxxx_0, tg_xyzz_xxxxxx_1, tg_xyzz_xxxxxy_0, tg_xyzz_xxxxxy_1, \
                                         tg_xyzz_xxxxxz_0, tg_xyzz_xxxxxz_1, tg_xyzz_xxxxyy_0, tg_xyzz_xxxxyy_1, \
                                         tg_xyzz_xxxxyz_0, tg_xyzz_xxxxyz_1, tg_xyzz_xxxxzz_0, tg_xyzz_xxxxzz_1, \
                                         tg_xyzz_xxxyyy_0, tg_xyzz_xxxyyy_1, tg_xyzz_xxxyyz_0, tg_xyzz_xxxyyz_1, \
                                         tg_xyzz_xxxyzz_0, tg_xyzz_xxxyzz_1, tg_xyzz_xxxzzz_0, tg_xyzz_xxxzzz_1, \
                                         tg_xyzz_xxyyyy_0, tg_xyzz_xxyyyy_1, tg_xyzz_xxyyyz_0, tg_xyzz_xxyyyz_1, \
                                         tg_xyzz_xxyyzz_0, tg_xyzz_xxyyzz_1, tg_xyzz_xxyzzz_0, tg_xyzz_xxyzzz_1, \
                                         tg_xyzz_xxzzzz_0, tg_xyzz_xxzzzz_1, tg_xyzz_xyyyyy_0, tg_xyzz_xyyyyy_1, \
                                         tg_xyzz_xyyyyz_0, tg_xyzz_xyyyyz_1, tg_xyzz_xyyyzz_0, tg_xyzz_xyyyzz_1, \
                                         tg_xyzz_xyyzzz_0, tg_xyzz_xyyzzz_1, tg_xyzz_xyzzzz_0, tg_xyzz_xyzzzz_1, \
                                         tg_xyzz_xzzzzz_0, tg_xyzz_xzzzzz_1, tg_xyzz_yyyyyy_0, tg_xyzz_yyyyyy_1, \
                                         tg_xyzz_yyyyyz_0, tg_xyzz_yyyyyz_1, tg_xyzz_yyyyzz_0, tg_xyzz_yyyyzz_1, \
                                         tg_xyzz_yyyzzz_0, tg_xyzz_yyyzzz_1, tg_xyzz_yyzzzz_0, tg_xyzz_yyzzzz_1, \
                                         tg_xyzz_yzzzzz_0, tg_xyzz_yzzzzz_1, tg_xyzz_zzzzzz_0, tg_xyzz_zzzzzz_1, \
                                         tg_xzzz_xxxxxx_0, tg_xzzz_xxxxxx_1, tg_xzzz_xxxxxy_0, tg_xzzz_xxxxxy_1, \
                                         tg_xzzz_xxxxxz_0, tg_xzzz_xxxxxz_1, tg_xzzz_xxxxyy_0, tg_xzzz_xxxxyy_1, \
                                         tg_xzzz_xxxxyz_0, tg_xzzz_xxxxyz_1, tg_xzzz_xxxxzz_0, tg_xzzz_xxxxzz_1, \
                                         tg_xzzz_xxxyyy_0, tg_xzzz_xxxyyy_1, tg_xzzz_xxxyyz_0, tg_xzzz_xxxyyz_1, \
                                         tg_xzzz_xxxyzz_0, tg_xzzz_xxxyzz_1, tg_xzzz_xxxzzz_0, tg_xzzz_xxxzzz_1, \
                                         tg_xzzz_xxyyyy_0, tg_xzzz_xxyyyy_1, tg_xzzz_xxyyyz_0, tg_xzzz_xxyyyz_1, \
                                         tg_xzzz_xxyyzz_0, tg_xzzz_xxyyzz_1, tg_xzzz_xxyzzz_0, tg_xzzz_xxyzzz_1, \
                                         tg_xzzz_xxzzzz_0, tg_xzzz_xxzzzz_1, tg_xzzz_xyyyyy_0, tg_xzzz_xyyyyy_1, \
                                         tg_xzzz_xyyyyz_0, tg_xzzz_xyyyyz_1, tg_xzzz_xyyyzz_0, tg_xzzz_xyyyzz_1, \
                                         tg_xzzz_xyyzzz_0, tg_xzzz_xyyzzz_1, tg_xzzz_xyzzzz_0, tg_xzzz_xyzzzz_1, \
                                         tg_xzzz_xzzzzz_0, tg_xzzz_xzzzzz_1, tg_xzzz_yyyyyy_0, tg_xzzz_yyyyyy_1, \
                                         tg_xzzz_yyyyyz_0, tg_xzzz_yyyyyz_1, tg_xzzz_yyyyzz_0, tg_xzzz_yyyyzz_1, \
                                         tg_xzzz_yyyzzz_0, tg_xzzz_yyyzzz_1, tg_xzzz_yyzzzz_0, tg_xzzz_yyzzzz_1, \
                                         tg_xzzz_yzzzzz_0, tg_xzzz_yzzzzz_1, tg_xzzz_zzzzzz_0, tg_xzzz_zzzzzz_1, \
                                         tg_yyyy_xxxxxx_0, tg_yyyy_xxxxxx_1, tg_yyyy_xxxxxy_0, tg_yyyy_xxxxxy_1, \
                                         tg_yyyy_xxxxxz_0, tg_yyyy_xxxxxz_1, tg_yyyy_xxxxyy_0, tg_yyyy_xxxxyy_1, \
                                         tg_yyyy_xxxxyz_0, tg_yyyy_xxxxyz_1, tg_yyyy_xxxxzz_0, tg_yyyy_xxxxzz_1, \
                                         tg_yyyy_xxxyyy_0, tg_yyyy_xxxyyy_1, tg_yyyy_xxxyyz_0, tg_yyyy_xxxyyz_1, \
                                         tg_yyyy_xxxyzz_0, tg_yyyy_xxxyzz_1, tg_yyyy_xxxzzz_0, tg_yyyy_xxxzzz_1, \
                                         tg_yyyy_xxyyyy_0, tg_yyyy_xxyyyy_1, tg_yyyy_xxyyyz_0, tg_yyyy_xxyyyz_1, \
                                         tg_yyyy_xxyyzz_0, tg_yyyy_xxyyzz_1, tg_yyyy_xxyzzz_0, tg_yyyy_xxyzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxyyz_xxxxxx_0[j] = pb_x * tg_xxyyz_xxxxxx_0[j] + fr * tg_xxyyz_xxxxxx_1[j] + fl1_fx * (tg_xyyz_xxxxxx_0[j] - tg_xyyz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxyyz_xxxxx_1[j];

                    tg_xxxyyz_xxxxxy_0[j] = pb_x * tg_xxyyz_xxxxxy_0[j] + fr * tg_xxyyz_xxxxxy_1[j] + fl1_fx * (tg_xyyz_xxxxxy_0[j] - tg_xyyz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyyz_xxxxy_1[j];

                    tg_xxxyyz_xxxxxz_0[j] = pb_x * tg_xxyyz_xxxxxz_0[j] + fr * tg_xxyyz_xxxxxz_1[j] + fl1_fx * (tg_xyyz_xxxxxz_0[j] - tg_xyyz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyyz_xxxxz_1[j];

                    tg_xxxyyz_xxxxyy_0[j] = pb_x * tg_xxyyz_xxxxyy_0[j] + fr * tg_xxyyz_xxxxyy_1[j] + fl1_fx * (tg_xyyz_xxxxyy_0[j] - tg_xyyz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyz_xxxyy_1[j];

                    tg_xxxyyz_xxxxyz_0[j] = pb_x * tg_xxyyz_xxxxyz_0[j] + fr * tg_xxyyz_xxxxyz_1[j] + fl1_fx * (tg_xyyz_xxxxyz_0[j] - tg_xyyz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyz_xxxyz_1[j];

                    tg_xxxyyz_xxxxzz_0[j] = pb_x * tg_xxyyz_xxxxzz_0[j] + fr * tg_xxyyz_xxxxzz_1[j] + fl1_fx * (tg_xyyz_xxxxzz_0[j] - tg_xyyz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyz_xxxzz_1[j];

                    tg_xxxyyz_xxxyyy_0[j] = pb_x * tg_xxyyz_xxxyyy_0[j] + fr * tg_xxyyz_xxxyyy_1[j] + fl1_fx * (tg_xyyz_xxxyyy_0[j] - tg_xyyz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyz_xxyyy_1[j];

                    tg_xxxyyz_xxxyyz_0[j] = pb_x * tg_xxyyz_xxxyyz_0[j] + fr * tg_xxyyz_xxxyyz_1[j] + fl1_fx * (tg_xyyz_xxxyyz_0[j] - tg_xyyz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyz_xxyyz_1[j];

                    tg_xxxyyz_xxxyzz_0[j] = pb_x * tg_xxyyz_xxxyzz_0[j] + fr * tg_xxyyz_xxxyzz_1[j] + fl1_fx * (tg_xyyz_xxxyzz_0[j] - tg_xyyz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyz_xxyzz_1[j];

                    tg_xxxyyz_xxxzzz_0[j] = pb_x * tg_xxyyz_xxxzzz_0[j] + fr * tg_xxyyz_xxxzzz_1[j] + fl1_fx * (tg_xyyz_xxxzzz_0[j] - tg_xyyz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyz_xxzzz_1[j];

                    tg_xxxyyz_xxyyyy_0[j] = pb_x * tg_xxyyz_xxyyyy_0[j] + fr * tg_xxyyz_xxyyyy_1[j] + fl1_fx * (tg_xyyz_xxyyyy_0[j] - tg_xyyz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxyyz_xyyyy_1[j];

                    tg_xxxyyz_xxyyyz_0[j] = pb_x * tg_xxyyz_xxyyyz_0[j] + fr * tg_xxyyz_xxyyyz_1[j] + fl1_fx * (tg_xyyz_xxyyyz_0[j] - tg_xyyz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyz_xyyyz_1[j];

                    tg_xxxyyz_xxyyzz_0[j] = pb_x * tg_xxyyz_xxyyzz_0[j] + fr * tg_xxyyz_xxyyzz_1[j] + fl1_fx * (tg_xyyz_xxyyzz_0[j] - tg_xyyz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyz_xyyzz_1[j];

                    tg_xxxyyz_xxyzzz_0[j] = pb_x * tg_xxyyz_xxyzzz_0[j] + fr * tg_xxyyz_xxyzzz_1[j] + fl1_fx * (tg_xyyz_xxyzzz_0[j] - tg_xyyz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyz_xyzzz_1[j];

                    tg_xxxyyz_xxzzzz_0[j] = pb_x * tg_xxyyz_xxzzzz_0[j] + fr * tg_xxyyz_xxzzzz_1[j] + fl1_fx * (tg_xyyz_xxzzzz_0[j] - tg_xyyz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyz_xzzzz_1[j];

                    tg_xxxyyz_xyyyyy_0[j] = pb_x * tg_xxyyz_xyyyyy_0[j] + fr * tg_xxyyz_xyyyyy_1[j] + fl1_fx * (tg_xyyz_xyyyyy_0[j] - tg_xyyz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyz_yyyyy_1[j];

                    tg_xxxyyz_xyyyyz_0[j] = pb_x * tg_xxyyz_xyyyyz_0[j] + fr * tg_xxyyz_xyyyyz_1[j] + fl1_fx * (tg_xyyz_xyyyyz_0[j] - tg_xyyz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyz_yyyyz_1[j];

                    tg_xxxyyz_xyyyzz_0[j] = pb_x * tg_xxyyz_xyyyzz_0[j] + fr * tg_xxyyz_xyyyzz_1[j] + fl1_fx * (tg_xyyz_xyyyzz_0[j] - tg_xyyz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyz_yyyzz_1[j];

                    tg_xxxyyz_xyyzzz_0[j] = pb_x * tg_xxyyz_xyyzzz_0[j] + fr * tg_xxyyz_xyyzzz_1[j] + fl1_fx * (tg_xyyz_xyyzzz_0[j] - tg_xyyz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyz_yyzzz_1[j];

                    tg_xxxyyz_xyzzzz_0[j] = pb_x * tg_xxyyz_xyzzzz_0[j] + fr * tg_xxyyz_xyzzzz_1[j] + fl1_fx * (tg_xyyz_xyzzzz_0[j] - tg_xyyz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyz_yzzzz_1[j];

                    tg_xxxyyz_xzzzzz_0[j] = pb_x * tg_xxyyz_xzzzzz_0[j] + fr * tg_xxyyz_xzzzzz_1[j] + fl1_fx * (tg_xyyz_xzzzzz_0[j] - tg_xyyz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyz_zzzzz_1[j];

                    tg_xxxyyz_yyyyyy_0[j] = pb_x * tg_xxyyz_yyyyyy_0[j] + fr * tg_xxyyz_yyyyyy_1[j] + fl1_fx * (tg_xyyz_yyyyyy_0[j] - tg_xyyz_yyyyyy_1[j] * fl1_fza);

                    tg_xxxyyz_yyyyyz_0[j] = pb_x * tg_xxyyz_yyyyyz_0[j] + fr * tg_xxyyz_yyyyyz_1[j] + fl1_fx * (tg_xyyz_yyyyyz_0[j] - tg_xyyz_yyyyyz_1[j] * fl1_fza);

                    tg_xxxyyz_yyyyzz_0[j] = pb_x * tg_xxyyz_yyyyzz_0[j] + fr * tg_xxyyz_yyyyzz_1[j] + fl1_fx * (tg_xyyz_yyyyzz_0[j] - tg_xyyz_yyyyzz_1[j] * fl1_fza);

                    tg_xxxyyz_yyyzzz_0[j] = pb_x * tg_xxyyz_yyyzzz_0[j] + fr * tg_xxyyz_yyyzzz_1[j] + fl1_fx * (tg_xyyz_yyyzzz_0[j] - tg_xyyz_yyyzzz_1[j] * fl1_fza);

                    tg_xxxyyz_yyzzzz_0[j] = pb_x * tg_xxyyz_yyzzzz_0[j] + fr * tg_xxyyz_yyzzzz_1[j] + fl1_fx * (tg_xyyz_yyzzzz_0[j] - tg_xyyz_yyzzzz_1[j] * fl1_fza);

                    tg_xxxyyz_yzzzzz_0[j] = pb_x * tg_xxyyz_yzzzzz_0[j] + fr * tg_xxyyz_yzzzzz_1[j] + fl1_fx * (tg_xyyz_yzzzzz_0[j] - tg_xyyz_yzzzzz_1[j] * fl1_fza);

                    tg_xxxyyz_zzzzzz_0[j] = pb_x * tg_xxyyz_zzzzzz_0[j] + fr * tg_xxyyz_zzzzzz_1[j] + fl1_fx * (tg_xyyz_zzzzzz_0[j] - tg_xyyz_zzzzzz_1[j] * fl1_fza);

                    tg_xxxyzz_xxxxxx_0[j] = pb_x * tg_xxyzz_xxxxxx_0[j] + fr * tg_xxyzz_xxxxxx_1[j] + fl1_fx * (tg_xyzz_xxxxxx_0[j] - tg_xyzz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxyzz_xxxxx_1[j];

                    tg_xxxyzz_xxxxxy_0[j] = pb_x * tg_xxyzz_xxxxxy_0[j] + fr * tg_xxyzz_xxxxxy_1[j] + fl1_fx * (tg_xyzz_xxxxxy_0[j] - tg_xyzz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyzz_xxxxy_1[j];

                    tg_xxxyzz_xxxxxz_0[j] = pb_x * tg_xxyzz_xxxxxz_0[j] + fr * tg_xxyzz_xxxxxz_1[j] + fl1_fx * (tg_xyzz_xxxxxz_0[j] - tg_xyzz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyzz_xxxxz_1[j];

                    tg_xxxyzz_xxxxyy_0[j] = pb_x * tg_xxyzz_xxxxyy_0[j] + fr * tg_xxyzz_xxxxyy_1[j] + fl1_fx * (tg_xyzz_xxxxyy_0[j] - tg_xyzz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyzz_xxxyy_1[j];

                    tg_xxxyzz_xxxxyz_0[j] = pb_x * tg_xxyzz_xxxxyz_0[j] + fr * tg_xxyzz_xxxxyz_1[j] + fl1_fx * (tg_xyzz_xxxxyz_0[j] - tg_xyzz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyzz_xxxyz_1[j];

                    tg_xxxyzz_xxxxzz_0[j] = pb_x * tg_xxyzz_xxxxzz_0[j] + fr * tg_xxyzz_xxxxzz_1[j] + fl1_fx * (tg_xyzz_xxxxzz_0[j] - tg_xyzz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyzz_xxxzz_1[j];

                    tg_xxxyzz_xxxyyy_0[j] = pb_x * tg_xxyzz_xxxyyy_0[j] + fr * tg_xxyzz_xxxyyy_1[j] + fl1_fx * (tg_xyzz_xxxyyy_0[j] - tg_xyzz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyzz_xxyyy_1[j];

                    tg_xxxyzz_xxxyyz_0[j] = pb_x * tg_xxyzz_xxxyyz_0[j] + fr * tg_xxyzz_xxxyyz_1[j] + fl1_fx * (tg_xyzz_xxxyyz_0[j] - tg_xyzz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyzz_xxyyz_1[j];

                    tg_xxxyzz_xxxyzz_0[j] = pb_x * tg_xxyzz_xxxyzz_0[j] + fr * tg_xxyzz_xxxyzz_1[j] + fl1_fx * (tg_xyzz_xxxyzz_0[j] - tg_xyzz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyzz_xxyzz_1[j];

                    tg_xxxyzz_xxxzzz_0[j] = pb_x * tg_xxyzz_xxxzzz_0[j] + fr * tg_xxyzz_xxxzzz_1[j] + fl1_fx * (tg_xyzz_xxxzzz_0[j] - tg_xyzz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyzz_xxzzz_1[j];

                    tg_xxxyzz_xxyyyy_0[j] = pb_x * tg_xxyzz_xxyyyy_0[j] + fr * tg_xxyzz_xxyyyy_1[j] + fl1_fx * (tg_xyzz_xxyyyy_0[j] - tg_xyzz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxyzz_xyyyy_1[j];

                    tg_xxxyzz_xxyyyz_0[j] = pb_x * tg_xxyzz_xxyyyz_0[j] + fr * tg_xxyzz_xxyyyz_1[j] + fl1_fx * (tg_xyzz_xxyyyz_0[j] - tg_xyzz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxyzz_xyyyz_1[j];

                    tg_xxxyzz_xxyyzz_0[j] = pb_x * tg_xxyzz_xxyyzz_0[j] + fr * tg_xxyzz_xxyyzz_1[j] + fl1_fx * (tg_xyzz_xxyyzz_0[j] - tg_xyzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyzz_xyyzz_1[j];

                    tg_xxxyzz_xxyzzz_0[j] = pb_x * tg_xxyzz_xxyzzz_0[j] + fr * tg_xxyzz_xxyzzz_1[j] + fl1_fx * (tg_xyzz_xxyzzz_0[j] - tg_xyzz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyzz_xyzzz_1[j];

                    tg_xxxyzz_xxzzzz_0[j] = pb_x * tg_xxyzz_xxzzzz_0[j] + fr * tg_xxyzz_xxzzzz_1[j] + fl1_fx * (tg_xyzz_xxzzzz_0[j] - tg_xyzz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyzz_xzzzz_1[j];

                    tg_xxxyzz_xyyyyy_0[j] = pb_x * tg_xxyzz_xyyyyy_0[j] + fr * tg_xxyzz_xyyyyy_1[j] + fl1_fx * (tg_xyzz_xyyyyy_0[j] - tg_xyzz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzz_yyyyy_1[j];

                    tg_xxxyzz_xyyyyz_0[j] = pb_x * tg_xxyzz_xyyyyz_0[j] + fr * tg_xxyzz_xyyyyz_1[j] + fl1_fx * (tg_xyzz_xyyyyz_0[j] - tg_xyzz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzz_yyyyz_1[j];

                    tg_xxxyzz_xyyyzz_0[j] = pb_x * tg_xxyzz_xyyyzz_0[j] + fr * tg_xxyzz_xyyyzz_1[j] + fl1_fx * (tg_xyzz_xyyyzz_0[j] - tg_xyzz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzz_yyyzz_1[j];

                    tg_xxxyzz_xyyzzz_0[j] = pb_x * tg_xxyzz_xyyzzz_0[j] + fr * tg_xxyzz_xyyzzz_1[j] + fl1_fx * (tg_xyzz_xyyzzz_0[j] - tg_xyzz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzz_yyzzz_1[j];

                    tg_xxxyzz_xyzzzz_0[j] = pb_x * tg_xxyzz_xyzzzz_0[j] + fr * tg_xxyzz_xyzzzz_1[j] + fl1_fx * (tg_xyzz_xyzzzz_0[j] - tg_xyzz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzz_yzzzz_1[j];

                    tg_xxxyzz_xzzzzz_0[j] = pb_x * tg_xxyzz_xzzzzz_0[j] + fr * tg_xxyzz_xzzzzz_1[j] + fl1_fx * (tg_xyzz_xzzzzz_0[j] - tg_xyzz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzz_zzzzz_1[j];

                    tg_xxxyzz_yyyyyy_0[j] = pb_x * tg_xxyzz_yyyyyy_0[j] + fr * tg_xxyzz_yyyyyy_1[j] + fl1_fx * (tg_xyzz_yyyyyy_0[j] - tg_xyzz_yyyyyy_1[j] * fl1_fza);

                    tg_xxxyzz_yyyyyz_0[j] = pb_x * tg_xxyzz_yyyyyz_0[j] + fr * tg_xxyzz_yyyyyz_1[j] + fl1_fx * (tg_xyzz_yyyyyz_0[j] - tg_xyzz_yyyyyz_1[j] * fl1_fza);

                    tg_xxxyzz_yyyyzz_0[j] = pb_x * tg_xxyzz_yyyyzz_0[j] + fr * tg_xxyzz_yyyyzz_1[j] + fl1_fx * (tg_xyzz_yyyyzz_0[j] - tg_xyzz_yyyyzz_1[j] * fl1_fza);

                    tg_xxxyzz_yyyzzz_0[j] = pb_x * tg_xxyzz_yyyzzz_0[j] + fr * tg_xxyzz_yyyzzz_1[j] + fl1_fx * (tg_xyzz_yyyzzz_0[j] - tg_xyzz_yyyzzz_1[j] * fl1_fza);

                    tg_xxxyzz_yyzzzz_0[j] = pb_x * tg_xxyzz_yyzzzz_0[j] + fr * tg_xxyzz_yyzzzz_1[j] + fl1_fx * (tg_xyzz_yyzzzz_0[j] - tg_xyzz_yyzzzz_1[j] * fl1_fza);

                    tg_xxxyzz_yzzzzz_0[j] = pb_x * tg_xxyzz_yzzzzz_0[j] + fr * tg_xxyzz_yzzzzz_1[j] + fl1_fx * (tg_xyzz_yzzzzz_0[j] - tg_xyzz_yzzzzz_1[j] * fl1_fza);

                    tg_xxxyzz_zzzzzz_0[j] = pb_x * tg_xxyzz_zzzzzz_0[j] + fr * tg_xxyzz_zzzzzz_1[j] + fl1_fx * (tg_xyzz_zzzzzz_0[j] - tg_xyzz_zzzzzz_1[j] * fl1_fza);

                    tg_xxxzzz_xxxxxx_0[j] = pb_x * tg_xxzzz_xxxxxx_0[j] + fr * tg_xxzzz_xxxxxx_1[j] + fl1_fx * (tg_xzzz_xxxxxx_0[j] - tg_xzzz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxzzz_xxxxx_1[j];

                    tg_xxxzzz_xxxxxy_0[j] = pb_x * tg_xxzzz_xxxxxy_0[j] + fr * tg_xxzzz_xxxxxy_1[j] + fl1_fx * (tg_xzzz_xxxxxy_0[j] - tg_xzzz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxzzz_xxxxy_1[j];

                    tg_xxxzzz_xxxxxz_0[j] = pb_x * tg_xxzzz_xxxxxz_0[j] + fr * tg_xxzzz_xxxxxz_1[j] + fl1_fx * (tg_xzzz_xxxxxz_0[j] - tg_xzzz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxzzz_xxxxz_1[j];

                    tg_xxxzzz_xxxxyy_0[j] = pb_x * tg_xxzzz_xxxxyy_0[j] + fr * tg_xxzzz_xxxxyy_1[j] + fl1_fx * (tg_xzzz_xxxxyy_0[j] - tg_xzzz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxzzz_xxxyy_1[j];

                    tg_xxxzzz_xxxxyz_0[j] = pb_x * tg_xxzzz_xxxxyz_0[j] + fr * tg_xxzzz_xxxxyz_1[j] + fl1_fx * (tg_xzzz_xxxxyz_0[j] - tg_xzzz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxzzz_xxxyz_1[j];

                    tg_xxxzzz_xxxxzz_0[j] = pb_x * tg_xxzzz_xxxxzz_0[j] + fr * tg_xxzzz_xxxxzz_1[j] + fl1_fx * (tg_xzzz_xxxxzz_0[j] - tg_xzzz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxzzz_xxxzz_1[j];

                    tg_xxxzzz_xxxyyy_0[j] = pb_x * tg_xxzzz_xxxyyy_0[j] + fr * tg_xxzzz_xxxyyy_1[j] + fl1_fx * (tg_xzzz_xxxyyy_0[j] - tg_xzzz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzzz_xxyyy_1[j];

                    tg_xxxzzz_xxxyyz_0[j] = pb_x * tg_xxzzz_xxxyyz_0[j] + fr * tg_xxzzz_xxxyyz_1[j] + fl1_fx * (tg_xzzz_xxxyyz_0[j] - tg_xzzz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzzz_xxyyz_1[j];

                    tg_xxxzzz_xxxyzz_0[j] = pb_x * tg_xxzzz_xxxyzz_0[j] + fr * tg_xxzzz_xxxyzz_1[j] + fl1_fx * (tg_xzzz_xxxyzz_0[j] - tg_xzzz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzzz_xxyzz_1[j];

                    tg_xxxzzz_xxxzzz_0[j] = pb_x * tg_xxzzz_xxxzzz_0[j] + fr * tg_xxzzz_xxxzzz_1[j] + fl1_fx * (tg_xzzz_xxxzzz_0[j] - tg_xzzz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzzz_xxzzz_1[j];

                    tg_xxxzzz_xxyyyy_0[j] = pb_x * tg_xxzzz_xxyyyy_0[j] + fr * tg_xxzzz_xxyyyy_1[j] + fl1_fx * (tg_xzzz_xxyyyy_0[j] - tg_xzzz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxzzz_xyyyy_1[j];

                    tg_xxxzzz_xxyyyz_0[j] = pb_x * tg_xxzzz_xxyyyz_0[j] + fr * tg_xxzzz_xxyyyz_1[j] + fl1_fx * (tg_xzzz_xxyyyz_0[j] - tg_xzzz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxzzz_xyyyz_1[j];

                    tg_xxxzzz_xxyyzz_0[j] = pb_x * tg_xxzzz_xxyyzz_0[j] + fr * tg_xxzzz_xxyyzz_1[j] + fl1_fx * (tg_xzzz_xxyyzz_0[j] - tg_xzzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxzzz_xyyzz_1[j];

                    tg_xxxzzz_xxyzzz_0[j] = pb_x * tg_xxzzz_xxyzzz_0[j] + fr * tg_xxzzz_xxyzzz_1[j] + fl1_fx * (tg_xzzz_xxyzzz_0[j] - tg_xzzz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxzzz_xyzzz_1[j];

                    tg_xxxzzz_xxzzzz_0[j] = pb_x * tg_xxzzz_xxzzzz_0[j] + fr * tg_xxzzz_xxzzzz_1[j] + fl1_fx * (tg_xzzz_xxzzzz_0[j] - tg_xzzz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxzzz_xzzzz_1[j];

                    tg_xxxzzz_xyyyyy_0[j] = pb_x * tg_xxzzz_xyyyyy_0[j] + fr * tg_xxzzz_xyyyyy_1[j] + fl1_fx * (tg_xzzz_xyyyyy_0[j] - tg_xzzz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzz_yyyyy_1[j];

                    tg_xxxzzz_xyyyyz_0[j] = pb_x * tg_xxzzz_xyyyyz_0[j] + fr * tg_xxzzz_xyyyyz_1[j] + fl1_fx * (tg_xzzz_xyyyyz_0[j] - tg_xzzz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzz_yyyyz_1[j];

                    tg_xxxzzz_xyyyzz_0[j] = pb_x * tg_xxzzz_xyyyzz_0[j] + fr * tg_xxzzz_xyyyzz_1[j] + fl1_fx * (tg_xzzz_xyyyzz_0[j] - tg_xzzz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzz_yyyzz_1[j];

                    tg_xxxzzz_xyyzzz_0[j] = pb_x * tg_xxzzz_xyyzzz_0[j] + fr * tg_xxzzz_xyyzzz_1[j] + fl1_fx * (tg_xzzz_xyyzzz_0[j] - tg_xzzz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzz_yyzzz_1[j];

                    tg_xxxzzz_xyzzzz_0[j] = pb_x * tg_xxzzz_xyzzzz_0[j] + fr * tg_xxzzz_xyzzzz_1[j] + fl1_fx * (tg_xzzz_xyzzzz_0[j] - tg_xzzz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzz_yzzzz_1[j];

                    tg_xxxzzz_xzzzzz_0[j] = pb_x * tg_xxzzz_xzzzzz_0[j] + fr * tg_xxzzz_xzzzzz_1[j] + fl1_fx * (tg_xzzz_xzzzzz_0[j] - tg_xzzz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzz_zzzzz_1[j];

                    tg_xxxzzz_yyyyyy_0[j] = pb_x * tg_xxzzz_yyyyyy_0[j] + fr * tg_xxzzz_yyyyyy_1[j] + fl1_fx * (tg_xzzz_yyyyyy_0[j] - tg_xzzz_yyyyyy_1[j] * fl1_fza);

                    tg_xxxzzz_yyyyyz_0[j] = pb_x * tg_xxzzz_yyyyyz_0[j] + fr * tg_xxzzz_yyyyyz_1[j] + fl1_fx * (tg_xzzz_yyyyyz_0[j] - tg_xzzz_yyyyyz_1[j] * fl1_fza);

                    tg_xxxzzz_yyyyzz_0[j] = pb_x * tg_xxzzz_yyyyzz_0[j] + fr * tg_xxzzz_yyyyzz_1[j] + fl1_fx * (tg_xzzz_yyyyzz_0[j] - tg_xzzz_yyyyzz_1[j] * fl1_fza);

                    tg_xxxzzz_yyyzzz_0[j] = pb_x * tg_xxzzz_yyyzzz_0[j] + fr * tg_xxzzz_yyyzzz_1[j] + fl1_fx * (tg_xzzz_yyyzzz_0[j] - tg_xzzz_yyyzzz_1[j] * fl1_fza);

                    tg_xxxzzz_yyzzzz_0[j] = pb_x * tg_xxzzz_yyzzzz_0[j] + fr * tg_xxzzz_yyzzzz_1[j] + fl1_fx * (tg_xzzz_yyzzzz_0[j] - tg_xzzz_yyzzzz_1[j] * fl1_fza);

                    tg_xxxzzz_yzzzzz_0[j] = pb_x * tg_xxzzz_yzzzzz_0[j] + fr * tg_xxzzz_yzzzzz_1[j] + fl1_fx * (tg_xzzz_yzzzzz_0[j] - tg_xzzz_yzzzzz_1[j] * fl1_fza);

                    tg_xxxzzz_zzzzzz_0[j] = pb_x * tg_xxzzz_zzzzzz_0[j] + fr * tg_xxzzz_zzzzzz_1[j] + fl1_fx * (tg_xzzz_zzzzzz_0[j] - tg_xzzz_zzzzzz_1[j] * fl1_fza);

                    tg_xxyyyy_xxxxxx_0[j] = pb_x * tg_xyyyy_xxxxxx_0[j] + fr * tg_xyyyy_xxxxxx_1[j] + 0.5 * fl1_fx * (tg_yyyy_xxxxxx_0[j] - tg_yyyy_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xyyyy_xxxxx_1[j];

                    tg_xxyyyy_xxxxxy_0[j] = pb_x * tg_xyyyy_xxxxxy_0[j] + fr * tg_xyyyy_xxxxxy_1[j] + 0.5 * fl1_fx * (tg_yyyy_xxxxxy_0[j] - tg_yyyy_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyyy_xxxxy_1[j];

                    tg_xxyyyy_xxxxxz_0[j] = pb_x * tg_xyyyy_xxxxxz_0[j] + fr * tg_xyyyy_xxxxxz_1[j] + 0.5 * fl1_fx * (tg_yyyy_xxxxxz_0[j] - tg_yyyy_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyyy_xxxxz_1[j];

                    tg_xxyyyy_xxxxyy_0[j] = pb_x * tg_xyyyy_xxxxyy_0[j] + fr * tg_xyyyy_xxxxyy_1[j] + 0.5 * fl1_fx * (tg_yyyy_xxxxyy_0[j] - tg_yyyy_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyy_xxxyy_1[j];

                    tg_xxyyyy_xxxxyz_0[j] = pb_x * tg_xyyyy_xxxxyz_0[j] + fr * tg_xyyyy_xxxxyz_1[j] + 0.5 * fl1_fx * (tg_yyyy_xxxxyz_0[j] - tg_yyyy_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyy_xxxyz_1[j];

                    tg_xxyyyy_xxxxzz_0[j] = pb_x * tg_xyyyy_xxxxzz_0[j] + fr * tg_xyyyy_xxxxzz_1[j] + 0.5 * fl1_fx * (tg_yyyy_xxxxzz_0[j] - tg_yyyy_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyy_xxxzz_1[j];

                    tg_xxyyyy_xxxyyy_0[j] = pb_x * tg_xyyyy_xxxyyy_0[j] + fr * tg_xyyyy_xxxyyy_1[j] + 0.5 * fl1_fx * (tg_yyyy_xxxyyy_0[j] - tg_yyyy_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyy_xxyyy_1[j];

                    tg_xxyyyy_xxxyyz_0[j] = pb_x * tg_xyyyy_xxxyyz_0[j] + fr * tg_xyyyy_xxxyyz_1[j] + 0.5 * fl1_fx * (tg_yyyy_xxxyyz_0[j] - tg_yyyy_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyy_xxyyz_1[j];

                    tg_xxyyyy_xxxyzz_0[j] = pb_x * tg_xyyyy_xxxyzz_0[j] + fr * tg_xyyyy_xxxyzz_1[j] + 0.5 * fl1_fx * (tg_yyyy_xxxyzz_0[j] - tg_yyyy_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyy_xxyzz_1[j];

                    tg_xxyyyy_xxxzzz_0[j] = pb_x * tg_xyyyy_xxxzzz_0[j] + fr * tg_xyyyy_xxxzzz_1[j] + 0.5 * fl1_fx * (tg_yyyy_xxxzzz_0[j] - tg_yyyy_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyy_xxzzz_1[j];

                    tg_xxyyyy_xxyyyy_0[j] = pb_x * tg_xyyyy_xxyyyy_0[j] + fr * tg_xyyyy_xxyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyy_xxyyyy_0[j] - tg_yyyy_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyyyy_xyyyy_1[j];

                    tg_xxyyyy_xxyyyz_0[j] = pb_x * tg_xyyyy_xxyyyz_0[j] + fr * tg_xyyyy_xxyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyy_xxyyyz_0[j] - tg_yyyy_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyy_xyyyz_1[j];

                    tg_xxyyyy_xxyyzz_0[j] = pb_x * tg_xyyyy_xxyyzz_0[j] + fr * tg_xyyyy_xxyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyy_xxyyzz_0[j] - tg_yyyy_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyy_xyyzz_1[j];

                    tg_xxyyyy_xxyzzz_0[j] = pb_x * tg_xyyyy_xxyzzz_0[j] + fr * tg_xyyyy_xxyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyy_xxyzzz_0[j] - tg_yyyy_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyy_xyzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSISI_294_392(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (294,392)

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
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_6_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_6_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_5_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
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

                auto tg_xyyyy_xxzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 294); 

                auto tg_xyyyy_xyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 295); 

                auto tg_xyyyy_xyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 296); 

                auto tg_xyyyy_xyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 297); 

                auto tg_xyyyy_xyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 298); 

                auto tg_xyyyy_xyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 299); 

                auto tg_xyyyy_xzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 300); 

                auto tg_xyyyy_yyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 301); 

                auto tg_xyyyy_yyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 302); 

                auto tg_xyyyy_yyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 303); 

                auto tg_xyyyy_yyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 304); 

                auto tg_xyyyy_yyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 305); 

                auto tg_xyyyy_yzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 306); 

                auto tg_xyyyy_zzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 307); 

                auto tg_xyyyz_xxxxxx_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 308); 

                auto tg_xyyyz_xxxxxy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 309); 

                auto tg_xyyyz_xxxxxz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 310); 

                auto tg_xyyyz_xxxxyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 311); 

                auto tg_xyyyz_xxxxyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 312); 

                auto tg_xyyyz_xxxxzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 313); 

                auto tg_xyyyz_xxxyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 314); 

                auto tg_xyyyz_xxxyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 315); 

                auto tg_xyyyz_xxxyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 316); 

                auto tg_xyyyz_xxxzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 317); 

                auto tg_xyyyz_xxyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 318); 

                auto tg_xyyyz_xxyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 319); 

                auto tg_xyyyz_xxyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 320); 

                auto tg_xyyyz_xxyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 321); 

                auto tg_xyyyz_xxzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 322); 

                auto tg_xyyyz_xyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 323); 

                auto tg_xyyyz_xyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 324); 

                auto tg_xyyyz_xyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 325); 

                auto tg_xyyyz_xyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 326); 

                auto tg_xyyyz_xyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 327); 

                auto tg_xyyyz_xzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 328); 

                auto tg_xyyyz_yyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 329); 

                auto tg_xyyyz_yyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 330); 

                auto tg_xyyyz_yyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 331); 

                auto tg_xyyyz_yyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 332); 

                auto tg_xyyyz_yyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 333); 

                auto tg_xyyyz_yzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 334); 

                auto tg_xyyyz_zzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 335); 

                auto tg_xyyzz_xxxxxx_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 336); 

                auto tg_xyyzz_xxxxxy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 337); 

                auto tg_xyyzz_xxxxxz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 338); 

                auto tg_xyyzz_xxxxyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 339); 

                auto tg_xyyzz_xxxxyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 340); 

                auto tg_xyyzz_xxxxzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 341); 

                auto tg_xyyzz_xxxyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 342); 

                auto tg_xyyzz_xxxyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 343); 

                auto tg_xyyzz_xxxyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 344); 

                auto tg_xyyzz_xxxzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 345); 

                auto tg_xyyzz_xxyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 346); 

                auto tg_xyyzz_xxyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 347); 

                auto tg_xyyzz_xxyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 348); 

                auto tg_xyyzz_xxyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 349); 

                auto tg_xyyzz_xxzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 350); 

                auto tg_xyyzz_xyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 351); 

                auto tg_xyyzz_xyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 352); 

                auto tg_xyyzz_xyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 353); 

                auto tg_xyyzz_xyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 354); 

                auto tg_xyyzz_xyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 355); 

                auto tg_xyyzz_xzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 356); 

                auto tg_xyyzz_yyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 357); 

                auto tg_xyyzz_yyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 358); 

                auto tg_xyyzz_yyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 359); 

                auto tg_xyyzz_yyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 360); 

                auto tg_xyyzz_yyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 361); 

                auto tg_xyyzz_yzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 362); 

                auto tg_xyyzz_zzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 363); 

                auto tg_xyzzz_xxxxxx_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 364); 

                auto tg_xyzzz_xxxxxy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 365); 

                auto tg_xyzzz_xxxxxz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 366); 

                auto tg_xyzzz_xxxxyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 367); 

                auto tg_xyzzz_xxxxyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 368); 

                auto tg_xyzzz_xxxxzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 369); 

                auto tg_xyzzz_xxxyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 370); 

                auto tg_xyzzz_xxxyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 371); 

                auto tg_xyzzz_xxxyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 372); 

                auto tg_xyzzz_xxxzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 373); 

                auto tg_xyzzz_xxyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 374); 

                auto tg_xyzzz_xxyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 375); 

                auto tg_xyzzz_xxyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 376); 

                auto tg_xyzzz_xxyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 377); 

                auto tg_xyzzz_xxzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 378); 

                auto tg_xyzzz_xyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 379); 

                auto tg_xyzzz_xyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 380); 

                auto tg_xyzzz_xyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 381); 

                auto tg_xyzzz_xyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 382); 

                auto tg_xyzzz_xyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 383); 

                auto tg_xyzzz_xzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 384); 

                auto tg_xyzzz_yyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 385); 

                auto tg_xyzzz_yyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 386); 

                auto tg_xyzzz_yyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 387); 

                auto tg_xyzzz_yyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 388); 

                auto tg_xyzzz_yyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 389); 

                auto tg_xyzzz_yzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 390); 

                auto tg_xyzzz_zzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 391); 

                auto tg_xyyyy_xxzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 294); 

                auto tg_xyyyy_xyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 295); 

                auto tg_xyyyy_xyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 296); 

                auto tg_xyyyy_xyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 297); 

                auto tg_xyyyy_xyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 298); 

                auto tg_xyyyy_xyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 299); 

                auto tg_xyyyy_xzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 300); 

                auto tg_xyyyy_yyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 301); 

                auto tg_xyyyy_yyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 302); 

                auto tg_xyyyy_yyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 303); 

                auto tg_xyyyy_yyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 304); 

                auto tg_xyyyy_yyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 305); 

                auto tg_xyyyy_yzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 306); 

                auto tg_xyyyy_zzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 307); 

                auto tg_xyyyz_xxxxxx_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 308); 

                auto tg_xyyyz_xxxxxy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 309); 

                auto tg_xyyyz_xxxxxz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 310); 

                auto tg_xyyyz_xxxxyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 311); 

                auto tg_xyyyz_xxxxyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 312); 

                auto tg_xyyyz_xxxxzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 313); 

                auto tg_xyyyz_xxxyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 314); 

                auto tg_xyyyz_xxxyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 315); 

                auto tg_xyyyz_xxxyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 316); 

                auto tg_xyyyz_xxxzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 317); 

                auto tg_xyyyz_xxyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 318); 

                auto tg_xyyyz_xxyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 319); 

                auto tg_xyyyz_xxyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 320); 

                auto tg_xyyyz_xxyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 321); 

                auto tg_xyyyz_xxzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 322); 

                auto tg_xyyyz_xyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 323); 

                auto tg_xyyyz_xyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 324); 

                auto tg_xyyyz_xyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 325); 

                auto tg_xyyyz_xyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 326); 

                auto tg_xyyyz_xyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 327); 

                auto tg_xyyyz_xzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 328); 

                auto tg_xyyyz_yyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 329); 

                auto tg_xyyyz_yyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 330); 

                auto tg_xyyyz_yyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 331); 

                auto tg_xyyyz_yyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 332); 

                auto tg_xyyyz_yyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 333); 

                auto tg_xyyyz_yzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 334); 

                auto tg_xyyyz_zzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 335); 

                auto tg_xyyzz_xxxxxx_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 336); 

                auto tg_xyyzz_xxxxxy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 337); 

                auto tg_xyyzz_xxxxxz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 338); 

                auto tg_xyyzz_xxxxyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 339); 

                auto tg_xyyzz_xxxxyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 340); 

                auto tg_xyyzz_xxxxzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 341); 

                auto tg_xyyzz_xxxyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 342); 

                auto tg_xyyzz_xxxyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 343); 

                auto tg_xyyzz_xxxyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 344); 

                auto tg_xyyzz_xxxzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 345); 

                auto tg_xyyzz_xxyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 346); 

                auto tg_xyyzz_xxyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 347); 

                auto tg_xyyzz_xxyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 348); 

                auto tg_xyyzz_xxyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 349); 

                auto tg_xyyzz_xxzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 350); 

                auto tg_xyyzz_xyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 351); 

                auto tg_xyyzz_xyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 352); 

                auto tg_xyyzz_xyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 353); 

                auto tg_xyyzz_xyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 354); 

                auto tg_xyyzz_xyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 355); 

                auto tg_xyyzz_xzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 356); 

                auto tg_xyyzz_yyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 357); 

                auto tg_xyyzz_yyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 358); 

                auto tg_xyyzz_yyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 359); 

                auto tg_xyyzz_yyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 360); 

                auto tg_xyyzz_yyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 361); 

                auto tg_xyyzz_yzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 362); 

                auto tg_xyyzz_zzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 363); 

                auto tg_xyzzz_xxxxxx_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 364); 

                auto tg_xyzzz_xxxxxy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 365); 

                auto tg_xyzzz_xxxxxz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 366); 

                auto tg_xyzzz_xxxxyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 367); 

                auto tg_xyzzz_xxxxyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 368); 

                auto tg_xyzzz_xxxxzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 369); 

                auto tg_xyzzz_xxxyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 370); 

                auto tg_xyzzz_xxxyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 371); 

                auto tg_xyzzz_xxxyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 372); 

                auto tg_xyzzz_xxxzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 373); 

                auto tg_xyzzz_xxyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 374); 

                auto tg_xyzzz_xxyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 375); 

                auto tg_xyzzz_xxyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 376); 

                auto tg_xyzzz_xxyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 377); 

                auto tg_xyzzz_xxzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 378); 

                auto tg_xyzzz_xyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 379); 

                auto tg_xyzzz_xyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 380); 

                auto tg_xyzzz_xyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 381); 

                auto tg_xyzzz_xyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 382); 

                auto tg_xyzzz_xyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 383); 

                auto tg_xyzzz_xzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 384); 

                auto tg_xyzzz_yyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 385); 

                auto tg_xyzzz_yyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 386); 

                auto tg_xyzzz_yyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 387); 

                auto tg_xyzzz_yyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 388); 

                auto tg_xyzzz_yyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 389); 

                auto tg_xyzzz_yzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 390); 

                auto tg_xyzzz_zzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 391); 

                auto tg_yyyy_xxzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 294); 

                auto tg_yyyy_xyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 295); 

                auto tg_yyyy_xyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 296); 

                auto tg_yyyy_xyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 297); 

                auto tg_yyyy_xyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 298); 

                auto tg_yyyy_xyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 299); 

                auto tg_yyyy_xzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 300); 

                auto tg_yyyy_yyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 301); 

                auto tg_yyyy_yyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 302); 

                auto tg_yyyy_yyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 303); 

                auto tg_yyyy_yyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 304); 

                auto tg_yyyy_yyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 305); 

                auto tg_yyyy_yzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 306); 

                auto tg_yyyy_zzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 307); 

                auto tg_yyyz_xxxxxx_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 308); 

                auto tg_yyyz_xxxxxy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 309); 

                auto tg_yyyz_xxxxxz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 310); 

                auto tg_yyyz_xxxxyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 311); 

                auto tg_yyyz_xxxxyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 312); 

                auto tg_yyyz_xxxxzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 313); 

                auto tg_yyyz_xxxyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 314); 

                auto tg_yyyz_xxxyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 315); 

                auto tg_yyyz_xxxyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 316); 

                auto tg_yyyz_xxxzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 317); 

                auto tg_yyyz_xxyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 318); 

                auto tg_yyyz_xxyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 319); 

                auto tg_yyyz_xxyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 320); 

                auto tg_yyyz_xxyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 321); 

                auto tg_yyyz_xxzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 322); 

                auto tg_yyyz_xyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 323); 

                auto tg_yyyz_xyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 324); 

                auto tg_yyyz_xyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 325); 

                auto tg_yyyz_xyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 326); 

                auto tg_yyyz_xyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 327); 

                auto tg_yyyz_xzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 328); 

                auto tg_yyyz_yyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 329); 

                auto tg_yyyz_yyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 330); 

                auto tg_yyyz_yyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 331); 

                auto tg_yyyz_yyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 332); 

                auto tg_yyyz_yyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 333); 

                auto tg_yyyz_yzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 334); 

                auto tg_yyyz_zzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 335); 

                auto tg_yyzz_xxxxxx_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 336); 

                auto tg_yyzz_xxxxxy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 337); 

                auto tg_yyzz_xxxxxz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 338); 

                auto tg_yyzz_xxxxyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 339); 

                auto tg_yyzz_xxxxyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 340); 

                auto tg_yyzz_xxxxzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 341); 

                auto tg_yyzz_xxxyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 342); 

                auto tg_yyzz_xxxyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 343); 

                auto tg_yyzz_xxxyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 344); 

                auto tg_yyzz_xxxzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 345); 

                auto tg_yyzz_xxyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 346); 

                auto tg_yyzz_xxyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 347); 

                auto tg_yyzz_xxyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 348); 

                auto tg_yyzz_xxyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 349); 

                auto tg_yyzz_xxzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 350); 

                auto tg_yyzz_xyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 351); 

                auto tg_yyzz_xyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 352); 

                auto tg_yyzz_xyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 353); 

                auto tg_yyzz_xyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 354); 

                auto tg_yyzz_xyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 355); 

                auto tg_yyzz_xzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 356); 

                auto tg_yyzz_yyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 357); 

                auto tg_yyzz_yyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 358); 

                auto tg_yyzz_yyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 359); 

                auto tg_yyzz_yyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 360); 

                auto tg_yyzz_yyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 361); 

                auto tg_yyzz_yzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 362); 

                auto tg_yyzz_zzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 363); 

                auto tg_yzzz_xxxxxx_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 364); 

                auto tg_yzzz_xxxxxy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 365); 

                auto tg_yzzz_xxxxxz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 366); 

                auto tg_yzzz_xxxxyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 367); 

                auto tg_yzzz_xxxxyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 368); 

                auto tg_yzzz_xxxxzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 369); 

                auto tg_yzzz_xxxyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 370); 

                auto tg_yzzz_xxxyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 371); 

                auto tg_yzzz_xxxyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 372); 

                auto tg_yzzz_xxxzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 373); 

                auto tg_yzzz_xxyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 374); 

                auto tg_yzzz_xxyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 375); 

                auto tg_yzzz_xxyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 376); 

                auto tg_yzzz_xxyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 377); 

                auto tg_yzzz_xxzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 378); 

                auto tg_yzzz_xyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 379); 

                auto tg_yzzz_xyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 380); 

                auto tg_yzzz_xyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 381); 

                auto tg_yzzz_xyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 382); 

                auto tg_yzzz_xyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 383); 

                auto tg_yzzz_xzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 384); 

                auto tg_yzzz_yyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 385); 

                auto tg_yzzz_yyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 386); 

                auto tg_yzzz_yyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 387); 

                auto tg_yzzz_yyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 388); 

                auto tg_yzzz_yyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 389); 

                auto tg_yzzz_yzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 390); 

                auto tg_yzzz_zzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 391); 

                auto tg_yyyy_xxzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 294); 

                auto tg_yyyy_xyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 295); 

                auto tg_yyyy_xyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 296); 

                auto tg_yyyy_xyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 297); 

                auto tg_yyyy_xyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 298); 

                auto tg_yyyy_xyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 299); 

                auto tg_yyyy_xzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 300); 

                auto tg_yyyy_yyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 301); 

                auto tg_yyyy_yyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 302); 

                auto tg_yyyy_yyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 303); 

                auto tg_yyyy_yyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 304); 

                auto tg_yyyy_yyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 305); 

                auto tg_yyyy_yzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 306); 

                auto tg_yyyy_zzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 307); 

                auto tg_yyyz_xxxxxx_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 308); 

                auto tg_yyyz_xxxxxy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 309); 

                auto tg_yyyz_xxxxxz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 310); 

                auto tg_yyyz_xxxxyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 311); 

                auto tg_yyyz_xxxxyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 312); 

                auto tg_yyyz_xxxxzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 313); 

                auto tg_yyyz_xxxyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 314); 

                auto tg_yyyz_xxxyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 315); 

                auto tg_yyyz_xxxyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 316); 

                auto tg_yyyz_xxxzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 317); 

                auto tg_yyyz_xxyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 318); 

                auto tg_yyyz_xxyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 319); 

                auto tg_yyyz_xxyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 320); 

                auto tg_yyyz_xxyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 321); 

                auto tg_yyyz_xxzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 322); 

                auto tg_yyyz_xyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 323); 

                auto tg_yyyz_xyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 324); 

                auto tg_yyyz_xyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 325); 

                auto tg_yyyz_xyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 326); 

                auto tg_yyyz_xyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 327); 

                auto tg_yyyz_xzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 328); 

                auto tg_yyyz_yyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 329); 

                auto tg_yyyz_yyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 330); 

                auto tg_yyyz_yyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 331); 

                auto tg_yyyz_yyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 332); 

                auto tg_yyyz_yyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 333); 

                auto tg_yyyz_yzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 334); 

                auto tg_yyyz_zzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 335); 

                auto tg_yyzz_xxxxxx_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 336); 

                auto tg_yyzz_xxxxxy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 337); 

                auto tg_yyzz_xxxxxz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 338); 

                auto tg_yyzz_xxxxyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 339); 

                auto tg_yyzz_xxxxyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 340); 

                auto tg_yyzz_xxxxzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 341); 

                auto tg_yyzz_xxxyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 342); 

                auto tg_yyzz_xxxyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 343); 

                auto tg_yyzz_xxxyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 344); 

                auto tg_yyzz_xxxzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 345); 

                auto tg_yyzz_xxyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 346); 

                auto tg_yyzz_xxyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 347); 

                auto tg_yyzz_xxyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 348); 

                auto tg_yyzz_xxyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 349); 

                auto tg_yyzz_xxzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 350); 

                auto tg_yyzz_xyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 351); 

                auto tg_yyzz_xyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 352); 

                auto tg_yyzz_xyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 353); 

                auto tg_yyzz_xyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 354); 

                auto tg_yyzz_xyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 355); 

                auto tg_yyzz_xzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 356); 

                auto tg_yyzz_yyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 357); 

                auto tg_yyzz_yyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 358); 

                auto tg_yyzz_yyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 359); 

                auto tg_yyzz_yyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 360); 

                auto tg_yyzz_yyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 361); 

                auto tg_yyzz_yzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 362); 

                auto tg_yyzz_zzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 363); 

                auto tg_yzzz_xxxxxx_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 364); 

                auto tg_yzzz_xxxxxy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 365); 

                auto tg_yzzz_xxxxxz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 366); 

                auto tg_yzzz_xxxxyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 367); 

                auto tg_yzzz_xxxxyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 368); 

                auto tg_yzzz_xxxxzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 369); 

                auto tg_yzzz_xxxyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 370); 

                auto tg_yzzz_xxxyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 371); 

                auto tg_yzzz_xxxyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 372); 

                auto tg_yzzz_xxxzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 373); 

                auto tg_yzzz_xxyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 374); 

                auto tg_yzzz_xxyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 375); 

                auto tg_yzzz_xxyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 376); 

                auto tg_yzzz_xxyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 377); 

                auto tg_yzzz_xxzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 378); 

                auto tg_yzzz_xyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 379); 

                auto tg_yzzz_xyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 380); 

                auto tg_yzzz_xyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 381); 

                auto tg_yzzz_xyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 382); 

                auto tg_yzzz_xyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 383); 

                auto tg_yzzz_xzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 384); 

                auto tg_yzzz_yyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 385); 

                auto tg_yzzz_yyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 386); 

                auto tg_yzzz_yyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 387); 

                auto tg_yzzz_yyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 388); 

                auto tg_yzzz_yyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 389); 

                auto tg_yzzz_yzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 390); 

                auto tg_yzzz_zzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 391); 

                auto tg_xyyyy_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 224); 

                auto tg_xyyyy_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 225); 

                auto tg_xyyyy_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 226); 

                auto tg_xyyyy_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 227); 

                auto tg_xyyyy_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 228); 

                auto tg_xyyyy_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 229); 

                auto tg_xyyyy_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 230); 

                auto tg_xyyyz_xxxxx_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 231); 

                auto tg_xyyyz_xxxxy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 232); 

                auto tg_xyyyz_xxxxz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 233); 

                auto tg_xyyyz_xxxyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 234); 

                auto tg_xyyyz_xxxyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 235); 

                auto tg_xyyyz_xxxzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 236); 

                auto tg_xyyyz_xxyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 237); 

                auto tg_xyyyz_xxyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 238); 

                auto tg_xyyyz_xxyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 239); 

                auto tg_xyyyz_xxzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 240); 

                auto tg_xyyyz_xyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 241); 

                auto tg_xyyyz_xyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 242); 

                auto tg_xyyyz_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 243); 

                auto tg_xyyyz_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 244); 

                auto tg_xyyyz_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 245); 

                auto tg_xyyyz_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 246); 

                auto tg_xyyyz_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 247); 

                auto tg_xyyyz_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 248); 

                auto tg_xyyyz_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 249); 

                auto tg_xyyyz_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 250); 

                auto tg_xyyyz_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 251); 

                auto tg_xyyzz_xxxxx_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 252); 

                auto tg_xyyzz_xxxxy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 253); 

                auto tg_xyyzz_xxxxz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 254); 

                auto tg_xyyzz_xxxyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 255); 

                auto tg_xyyzz_xxxyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 256); 

                auto tg_xyyzz_xxxzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 257); 

                auto tg_xyyzz_xxyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 258); 

                auto tg_xyyzz_xxyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 259); 

                auto tg_xyyzz_xxyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 260); 

                auto tg_xyyzz_xxzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 261); 

                auto tg_xyyzz_xyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 262); 

                auto tg_xyyzz_xyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 263); 

                auto tg_xyyzz_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 264); 

                auto tg_xyyzz_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 265); 

                auto tg_xyyzz_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 266); 

                auto tg_xyyzz_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 267); 

                auto tg_xyyzz_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 268); 

                auto tg_xyyzz_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 269); 

                auto tg_xyyzz_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 270); 

                auto tg_xyyzz_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 271); 

                auto tg_xyyzz_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 272); 

                auto tg_xyzzz_xxxxx_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 273); 

                auto tg_xyzzz_xxxxy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 274); 

                auto tg_xyzzz_xxxxz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 275); 

                auto tg_xyzzz_xxxyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 276); 

                auto tg_xyzzz_xxxyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 277); 

                auto tg_xyzzz_xxxzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 278); 

                auto tg_xyzzz_xxyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 279); 

                auto tg_xyzzz_xxyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 280); 

                auto tg_xyzzz_xxyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 281); 

                auto tg_xyzzz_xxzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 282); 

                auto tg_xyzzz_xyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 283); 

                auto tg_xyzzz_xyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 284); 

                auto tg_xyzzz_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 285); 

                auto tg_xyzzz_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 286); 

                auto tg_xyzzz_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 287); 

                auto tg_xyzzz_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 288); 

                auto tg_xyzzz_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 289); 

                auto tg_xyzzz_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 290); 

                auto tg_xyzzz_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 291); 

                auto tg_xyzzz_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 292); 

                auto tg_xyzzz_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 293); 

                // set up pointers to integrals

                auto tg_xxyyyy_xxzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 294); 

                auto tg_xxyyyy_xyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 295); 

                auto tg_xxyyyy_xyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 296); 

                auto tg_xxyyyy_xyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 297); 

                auto tg_xxyyyy_xyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 298); 

                auto tg_xxyyyy_xyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 299); 

                auto tg_xxyyyy_xzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 300); 

                auto tg_xxyyyy_yyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 301); 

                auto tg_xxyyyy_yyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 302); 

                auto tg_xxyyyy_yyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 303); 

                auto tg_xxyyyy_yyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 304); 

                auto tg_xxyyyy_yyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 305); 

                auto tg_xxyyyy_yzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 306); 

                auto tg_xxyyyy_zzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 307); 

                auto tg_xxyyyz_xxxxxx_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 308); 

                auto tg_xxyyyz_xxxxxy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 309); 

                auto tg_xxyyyz_xxxxxz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 310); 

                auto tg_xxyyyz_xxxxyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 311); 

                auto tg_xxyyyz_xxxxyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 312); 

                auto tg_xxyyyz_xxxxzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 313); 

                auto tg_xxyyyz_xxxyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 314); 

                auto tg_xxyyyz_xxxyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 315); 

                auto tg_xxyyyz_xxxyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 316); 

                auto tg_xxyyyz_xxxzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 317); 

                auto tg_xxyyyz_xxyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 318); 

                auto tg_xxyyyz_xxyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 319); 

                auto tg_xxyyyz_xxyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 320); 

                auto tg_xxyyyz_xxyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 321); 

                auto tg_xxyyyz_xxzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 322); 

                auto tg_xxyyyz_xyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 323); 

                auto tg_xxyyyz_xyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 324); 

                auto tg_xxyyyz_xyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 325); 

                auto tg_xxyyyz_xyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 326); 

                auto tg_xxyyyz_xyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 327); 

                auto tg_xxyyyz_xzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 328); 

                auto tg_xxyyyz_yyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 329); 

                auto tg_xxyyyz_yyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 330); 

                auto tg_xxyyyz_yyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 331); 

                auto tg_xxyyyz_yyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 332); 

                auto tg_xxyyyz_yyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 333); 

                auto tg_xxyyyz_yzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 334); 

                auto tg_xxyyyz_zzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 335); 

                auto tg_xxyyzz_xxxxxx_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 336); 

                auto tg_xxyyzz_xxxxxy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 337); 

                auto tg_xxyyzz_xxxxxz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 338); 

                auto tg_xxyyzz_xxxxyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 339); 

                auto tg_xxyyzz_xxxxyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 340); 

                auto tg_xxyyzz_xxxxzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 341); 

                auto tg_xxyyzz_xxxyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 342); 

                auto tg_xxyyzz_xxxyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 343); 

                auto tg_xxyyzz_xxxyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 344); 

                auto tg_xxyyzz_xxxzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 345); 

                auto tg_xxyyzz_xxyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 346); 

                auto tg_xxyyzz_xxyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 347); 

                auto tg_xxyyzz_xxyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 348); 

                auto tg_xxyyzz_xxyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 349); 

                auto tg_xxyyzz_xxzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 350); 

                auto tg_xxyyzz_xyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 351); 

                auto tg_xxyyzz_xyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 352); 

                auto tg_xxyyzz_xyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 353); 

                auto tg_xxyyzz_xyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 354); 

                auto tg_xxyyzz_xyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 355); 

                auto tg_xxyyzz_xzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 356); 

                auto tg_xxyyzz_yyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 357); 

                auto tg_xxyyzz_yyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 358); 

                auto tg_xxyyzz_yyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 359); 

                auto tg_xxyyzz_yyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 360); 

                auto tg_xxyyzz_yyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 361); 

                auto tg_xxyyzz_yzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 362); 

                auto tg_xxyyzz_zzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 363); 

                auto tg_xxyzzz_xxxxxx_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 364); 

                auto tg_xxyzzz_xxxxxy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 365); 

                auto tg_xxyzzz_xxxxxz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 366); 

                auto tg_xxyzzz_xxxxyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 367); 

                auto tg_xxyzzz_xxxxyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 368); 

                auto tg_xxyzzz_xxxxzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 369); 

                auto tg_xxyzzz_xxxyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 370); 

                auto tg_xxyzzz_xxxyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 371); 

                auto tg_xxyzzz_xxxyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 372); 

                auto tg_xxyzzz_xxxzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 373); 

                auto tg_xxyzzz_xxyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 374); 

                auto tg_xxyzzz_xxyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 375); 

                auto tg_xxyzzz_xxyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 376); 

                auto tg_xxyzzz_xxyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 377); 

                auto tg_xxyzzz_xxzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 378); 

                auto tg_xxyzzz_xyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 379); 

                auto tg_xxyzzz_xyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 380); 

                auto tg_xxyzzz_xyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 381); 

                auto tg_xxyzzz_xyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 382); 

                auto tg_xxyzzz_xyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 383); 

                auto tg_xxyzzz_xzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 384); 

                auto tg_xxyzzz_yyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 385); 

                auto tg_xxyzzz_yyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 386); 

                auto tg_xxyzzz_yyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 387); 

                auto tg_xxyzzz_yyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 388); 

                auto tg_xxyzzz_yyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 389); 

                auto tg_xxyzzz_yzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 390); 

                auto tg_xxyzzz_zzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 391); 

                // Batch of Integrals (294,392)

                #pragma omp simd aligned(fxn, fza, tg_xxyyyy_xxzzzz_0, tg_xxyyyy_xyyyyy_0, tg_xxyyyy_xyyyyz_0, \
                                         tg_xxyyyy_xyyyzz_0, tg_xxyyyy_xyyzzz_0, tg_xxyyyy_xyzzzz_0, tg_xxyyyy_xzzzzz_0, \
                                         tg_xxyyyy_yyyyyy_0, tg_xxyyyy_yyyyyz_0, tg_xxyyyy_yyyyzz_0, tg_xxyyyy_yyyzzz_0, \
                                         tg_xxyyyy_yyzzzz_0, tg_xxyyyy_yzzzzz_0, tg_xxyyyy_zzzzzz_0, tg_xxyyyz_xxxxxx_0, \
                                         tg_xxyyyz_xxxxxy_0, tg_xxyyyz_xxxxxz_0, tg_xxyyyz_xxxxyy_0, tg_xxyyyz_xxxxyz_0, \
                                         tg_xxyyyz_xxxxzz_0, tg_xxyyyz_xxxyyy_0, tg_xxyyyz_xxxyyz_0, tg_xxyyyz_xxxyzz_0, \
                                         tg_xxyyyz_xxxzzz_0, tg_xxyyyz_xxyyyy_0, tg_xxyyyz_xxyyyz_0, tg_xxyyyz_xxyyzz_0, \
                                         tg_xxyyyz_xxyzzz_0, tg_xxyyyz_xxzzzz_0, tg_xxyyyz_xyyyyy_0, tg_xxyyyz_xyyyyz_0, \
                                         tg_xxyyyz_xyyyzz_0, tg_xxyyyz_xyyzzz_0, tg_xxyyyz_xyzzzz_0, tg_xxyyyz_xzzzzz_0, \
                                         tg_xxyyyz_yyyyyy_0, tg_xxyyyz_yyyyyz_0, tg_xxyyyz_yyyyzz_0, tg_xxyyyz_yyyzzz_0, \
                                         tg_xxyyyz_yyzzzz_0, tg_xxyyyz_yzzzzz_0, tg_xxyyyz_zzzzzz_0, tg_xxyyzz_xxxxxx_0, \
                                         tg_xxyyzz_xxxxxy_0, tg_xxyyzz_xxxxxz_0, tg_xxyyzz_xxxxyy_0, tg_xxyyzz_xxxxyz_0, \
                                         tg_xxyyzz_xxxxzz_0, tg_xxyyzz_xxxyyy_0, tg_xxyyzz_xxxyyz_0, tg_xxyyzz_xxxyzz_0, \
                                         tg_xxyyzz_xxxzzz_0, tg_xxyyzz_xxyyyy_0, tg_xxyyzz_xxyyyz_0, tg_xxyyzz_xxyyzz_0, \
                                         tg_xxyyzz_xxyzzz_0, tg_xxyyzz_xxzzzz_0, tg_xxyyzz_xyyyyy_0, tg_xxyyzz_xyyyyz_0, \
                                         tg_xxyyzz_xyyyzz_0, tg_xxyyzz_xyyzzz_0, tg_xxyyzz_xyzzzz_0, tg_xxyyzz_xzzzzz_0, \
                                         tg_xxyyzz_yyyyyy_0, tg_xxyyzz_yyyyyz_0, tg_xxyyzz_yyyyzz_0, tg_xxyyzz_yyyzzz_0, \
                                         tg_xxyyzz_yyzzzz_0, tg_xxyyzz_yzzzzz_0, tg_xxyyzz_zzzzzz_0, tg_xxyzzz_xxxxxx_0, \
                                         tg_xxyzzz_xxxxxy_0, tg_xxyzzz_xxxxxz_0, tg_xxyzzz_xxxxyy_0, tg_xxyzzz_xxxxyz_0, \
                                         tg_xxyzzz_xxxxzz_0, tg_xxyzzz_xxxyyy_0, tg_xxyzzz_xxxyyz_0, tg_xxyzzz_xxxyzz_0, \
                                         tg_xxyzzz_xxxzzz_0, tg_xxyzzz_xxyyyy_0, tg_xxyzzz_xxyyyz_0, tg_xxyzzz_xxyyzz_0, \
                                         tg_xxyzzz_xxyzzz_0, tg_xxyzzz_xxzzzz_0, tg_xxyzzz_xyyyyy_0, tg_xxyzzz_xyyyyz_0, \
                                         tg_xxyzzz_xyyyzz_0, tg_xxyzzz_xyyzzz_0, tg_xxyzzz_xyzzzz_0, tg_xxyzzz_xzzzzz_0, \
                                         tg_xxyzzz_yyyyyy_0, tg_xxyzzz_yyyyyz_0, tg_xxyzzz_yyyyzz_0, tg_xxyzzz_yyyzzz_0, \
                                         tg_xxyzzz_yyzzzz_0, tg_xxyzzz_yzzzzz_0, tg_xxyzzz_zzzzzz_0, tg_xyyyy_xxzzzz_0, \
                                         tg_xyyyy_xxzzzz_1, tg_xyyyy_xyyyyy_0, tg_xyyyy_xyyyyy_1, tg_xyyyy_xyyyyz_0, \
                                         tg_xyyyy_xyyyyz_1, tg_xyyyy_xyyyzz_0, tg_xyyyy_xyyyzz_1, tg_xyyyy_xyyzzz_0, \
                                         tg_xyyyy_xyyzzz_1, tg_xyyyy_xyzzzz_0, tg_xyyyy_xyzzzz_1, tg_xyyyy_xzzzz_1, \
                                         tg_xyyyy_xzzzzz_0, tg_xyyyy_xzzzzz_1, tg_xyyyy_yyyyy_1, tg_xyyyy_yyyyyy_0, \
                                         tg_xyyyy_yyyyyy_1, tg_xyyyy_yyyyyz_0, tg_xyyyy_yyyyyz_1, tg_xyyyy_yyyyz_1, \
                                         tg_xyyyy_yyyyzz_0, tg_xyyyy_yyyyzz_1, tg_xyyyy_yyyzz_1, tg_xyyyy_yyyzzz_0, \
                                         tg_xyyyy_yyyzzz_1, tg_xyyyy_yyzzz_1, tg_xyyyy_yyzzzz_0, tg_xyyyy_yyzzzz_1, \
                                         tg_xyyyy_yzzzz_1, tg_xyyyy_yzzzzz_0, tg_xyyyy_yzzzzz_1, tg_xyyyy_zzzzz_1, \
                                         tg_xyyyy_zzzzzz_0, tg_xyyyy_zzzzzz_1, tg_xyyyz_xxxxx_1, tg_xyyyz_xxxxxx_0, \
                                         tg_xyyyz_xxxxxx_1, tg_xyyyz_xxxxxy_0, tg_xyyyz_xxxxxy_1, tg_xyyyz_xxxxxz_0, \
                                         tg_xyyyz_xxxxxz_1, tg_xyyyz_xxxxy_1, tg_xyyyz_xxxxyy_0, tg_xyyyz_xxxxyy_1, \
                                         tg_xyyyz_xxxxyz_0, tg_xyyyz_xxxxyz_1, tg_xyyyz_xxxxz_1, tg_xyyyz_xxxxzz_0, \
                                         tg_xyyyz_xxxxzz_1, tg_xyyyz_xxxyy_1, tg_xyyyz_xxxyyy_0, tg_xyyyz_xxxyyy_1, \
                                         tg_xyyyz_xxxyyz_0, tg_xyyyz_xxxyyz_1, tg_xyyyz_xxxyz_1, tg_xyyyz_xxxyzz_0, \
                                         tg_xyyyz_xxxyzz_1, tg_xyyyz_xxxzz_1, tg_xyyyz_xxxzzz_0, tg_xyyyz_xxxzzz_1, \
                                         tg_xyyyz_xxyyy_1, tg_xyyyz_xxyyyy_0, tg_xyyyz_xxyyyy_1, tg_xyyyz_xxyyyz_0, \
                                         tg_xyyyz_xxyyyz_1, tg_xyyyz_xxyyz_1, tg_xyyyz_xxyyzz_0, tg_xyyyz_xxyyzz_1, \
                                         tg_xyyyz_xxyzz_1, tg_xyyyz_xxyzzz_0, tg_xyyyz_xxyzzz_1, tg_xyyyz_xxzzz_1, \
                                         tg_xyyyz_xxzzzz_0, tg_xyyyz_xxzzzz_1, tg_xyyyz_xyyyy_1, tg_xyyyz_xyyyyy_0, \
                                         tg_xyyyz_xyyyyy_1, tg_xyyyz_xyyyyz_0, tg_xyyyz_xyyyyz_1, tg_xyyyz_xyyyz_1, \
                                         tg_xyyyz_xyyyzz_0, tg_xyyyz_xyyyzz_1, tg_xyyyz_xyyzz_1, tg_xyyyz_xyyzzz_0, \
                                         tg_xyyyz_xyyzzz_1, tg_xyyyz_xyzzz_1, tg_xyyyz_xyzzzz_0, tg_xyyyz_xyzzzz_1, \
                                         tg_xyyyz_xzzzz_1, tg_xyyyz_xzzzzz_0, tg_xyyyz_xzzzzz_1, tg_xyyyz_yyyyy_1, \
                                         tg_xyyyz_yyyyyy_0, tg_xyyyz_yyyyyy_1, tg_xyyyz_yyyyyz_0, tg_xyyyz_yyyyyz_1, \
                                         tg_xyyyz_yyyyz_1, tg_xyyyz_yyyyzz_0, tg_xyyyz_yyyyzz_1, tg_xyyyz_yyyzz_1, \
                                         tg_xyyyz_yyyzzz_0, tg_xyyyz_yyyzzz_1, tg_xyyyz_yyzzz_1, tg_xyyyz_yyzzzz_0, \
                                         tg_xyyyz_yyzzzz_1, tg_xyyyz_yzzzz_1, tg_xyyyz_yzzzzz_0, tg_xyyyz_yzzzzz_1, \
                                         tg_xyyyz_zzzzz_1, tg_xyyyz_zzzzzz_0, tg_xyyyz_zzzzzz_1, tg_xyyzz_xxxxx_1, \
                                         tg_xyyzz_xxxxxx_0, tg_xyyzz_xxxxxx_1, tg_xyyzz_xxxxxy_0, tg_xyyzz_xxxxxy_1, \
                                         tg_xyyzz_xxxxxz_0, tg_xyyzz_xxxxxz_1, tg_xyyzz_xxxxy_1, tg_xyyzz_xxxxyy_0, \
                                         tg_xyyzz_xxxxyy_1, tg_xyyzz_xxxxyz_0, tg_xyyzz_xxxxyz_1, tg_xyyzz_xxxxz_1, \
                                         tg_xyyzz_xxxxzz_0, tg_xyyzz_xxxxzz_1, tg_xyyzz_xxxyy_1, tg_xyyzz_xxxyyy_0, \
                                         tg_xyyzz_xxxyyy_1, tg_xyyzz_xxxyyz_0, tg_xyyzz_xxxyyz_1, tg_xyyzz_xxxyz_1, \
                                         tg_xyyzz_xxxyzz_0, tg_xyyzz_xxxyzz_1, tg_xyyzz_xxxzz_1, tg_xyyzz_xxxzzz_0, \
                                         tg_xyyzz_xxxzzz_1, tg_xyyzz_xxyyy_1, tg_xyyzz_xxyyyy_0, tg_xyyzz_xxyyyy_1, \
                                         tg_xyyzz_xxyyyz_0, tg_xyyzz_xxyyyz_1, tg_xyyzz_xxyyz_1, tg_xyyzz_xxyyzz_0, \
                                         tg_xyyzz_xxyyzz_1, tg_xyyzz_xxyzz_1, tg_xyyzz_xxyzzz_0, tg_xyyzz_xxyzzz_1, \
                                         tg_xyyzz_xxzzz_1, tg_xyyzz_xxzzzz_0, tg_xyyzz_xxzzzz_1, tg_xyyzz_xyyyy_1, \
                                         tg_xyyzz_xyyyyy_0, tg_xyyzz_xyyyyy_1, tg_xyyzz_xyyyyz_0, tg_xyyzz_xyyyyz_1, \
                                         tg_xyyzz_xyyyz_1, tg_xyyzz_xyyyzz_0, tg_xyyzz_xyyyzz_1, tg_xyyzz_xyyzz_1, \
                                         tg_xyyzz_xyyzzz_0, tg_xyyzz_xyyzzz_1, tg_xyyzz_xyzzz_1, tg_xyyzz_xyzzzz_0, \
                                         tg_xyyzz_xyzzzz_1, tg_xyyzz_xzzzz_1, tg_xyyzz_xzzzzz_0, tg_xyyzz_xzzzzz_1, \
                                         tg_xyyzz_yyyyy_1, tg_xyyzz_yyyyyy_0, tg_xyyzz_yyyyyy_1, tg_xyyzz_yyyyyz_0, \
                                         tg_xyyzz_yyyyyz_1, tg_xyyzz_yyyyz_1, tg_xyyzz_yyyyzz_0, tg_xyyzz_yyyyzz_1, \
                                         tg_xyyzz_yyyzz_1, tg_xyyzz_yyyzzz_0, tg_xyyzz_yyyzzz_1, tg_xyyzz_yyzzz_1, \
                                         tg_xyyzz_yyzzzz_0, tg_xyyzz_yyzzzz_1, tg_xyyzz_yzzzz_1, tg_xyyzz_yzzzzz_0, \
                                         tg_xyyzz_yzzzzz_1, tg_xyyzz_zzzzz_1, tg_xyyzz_zzzzzz_0, tg_xyyzz_zzzzzz_1, \
                                         tg_xyzzz_xxxxx_1, tg_xyzzz_xxxxxx_0, tg_xyzzz_xxxxxx_1, tg_xyzzz_xxxxxy_0, \
                                         tg_xyzzz_xxxxxy_1, tg_xyzzz_xxxxxz_0, tg_xyzzz_xxxxxz_1, tg_xyzzz_xxxxy_1, \
                                         tg_xyzzz_xxxxyy_0, tg_xyzzz_xxxxyy_1, tg_xyzzz_xxxxyz_0, tg_xyzzz_xxxxyz_1, \
                                         tg_xyzzz_xxxxz_1, tg_xyzzz_xxxxzz_0, tg_xyzzz_xxxxzz_1, tg_xyzzz_xxxyy_1, \
                                         tg_xyzzz_xxxyyy_0, tg_xyzzz_xxxyyy_1, tg_xyzzz_xxxyyz_0, tg_xyzzz_xxxyyz_1, \
                                         tg_xyzzz_xxxyz_1, tg_xyzzz_xxxyzz_0, tg_xyzzz_xxxyzz_1, tg_xyzzz_xxxzz_1, \
                                         tg_xyzzz_xxxzzz_0, tg_xyzzz_xxxzzz_1, tg_xyzzz_xxyyy_1, tg_xyzzz_xxyyyy_0, \
                                         tg_xyzzz_xxyyyy_1, tg_xyzzz_xxyyyz_0, tg_xyzzz_xxyyyz_1, tg_xyzzz_xxyyz_1, \
                                         tg_xyzzz_xxyyzz_0, tg_xyzzz_xxyyzz_1, tg_xyzzz_xxyzz_1, tg_xyzzz_xxyzzz_0, \
                                         tg_xyzzz_xxyzzz_1, tg_xyzzz_xxzzz_1, tg_xyzzz_xxzzzz_0, tg_xyzzz_xxzzzz_1, \
                                         tg_xyzzz_xyyyy_1, tg_xyzzz_xyyyyy_0, tg_xyzzz_xyyyyy_1, tg_xyzzz_xyyyyz_0, \
                                         tg_xyzzz_xyyyyz_1, tg_xyzzz_xyyyz_1, tg_xyzzz_xyyyzz_0, tg_xyzzz_xyyyzz_1, \
                                         tg_xyzzz_xyyzz_1, tg_xyzzz_xyyzzz_0, tg_xyzzz_xyyzzz_1, tg_xyzzz_xyzzz_1, \
                                         tg_xyzzz_xyzzzz_0, tg_xyzzz_xyzzzz_1, tg_xyzzz_xzzzz_1, tg_xyzzz_xzzzzz_0, \
                                         tg_xyzzz_xzzzzz_1, tg_xyzzz_yyyyy_1, tg_xyzzz_yyyyyy_0, tg_xyzzz_yyyyyy_1, \
                                         tg_xyzzz_yyyyyz_0, tg_xyzzz_yyyyyz_1, tg_xyzzz_yyyyz_1, tg_xyzzz_yyyyzz_0, \
                                         tg_xyzzz_yyyyzz_1, tg_xyzzz_yyyzz_1, tg_xyzzz_yyyzzz_0, tg_xyzzz_yyyzzz_1, \
                                         tg_xyzzz_yyzzz_1, tg_xyzzz_yyzzzz_0, tg_xyzzz_yyzzzz_1, tg_xyzzz_yzzzz_1, \
                                         tg_xyzzz_yzzzzz_0, tg_xyzzz_yzzzzz_1, tg_xyzzz_zzzzz_1, tg_xyzzz_zzzzzz_0, \
                                         tg_xyzzz_zzzzzz_1, tg_yyyy_xxzzzz_0, tg_yyyy_xxzzzz_1, tg_yyyy_xyyyyy_0, \
                                         tg_yyyy_xyyyyy_1, tg_yyyy_xyyyyz_0, tg_yyyy_xyyyyz_1, tg_yyyy_xyyyzz_0, \
                                         tg_yyyy_xyyyzz_1, tg_yyyy_xyyzzz_0, tg_yyyy_xyyzzz_1, tg_yyyy_xyzzzz_0, \
                                         tg_yyyy_xyzzzz_1, tg_yyyy_xzzzzz_0, tg_yyyy_xzzzzz_1, tg_yyyy_yyyyyy_0, \
                                         tg_yyyy_yyyyyy_1, tg_yyyy_yyyyyz_0, tg_yyyy_yyyyyz_1, tg_yyyy_yyyyzz_0, \
                                         tg_yyyy_yyyyzz_1, tg_yyyy_yyyzzz_0, tg_yyyy_yyyzzz_1, tg_yyyy_yyzzzz_0, \
                                         tg_yyyy_yyzzzz_1, tg_yyyy_yzzzzz_0, tg_yyyy_yzzzzz_1, tg_yyyy_zzzzzz_0, \
                                         tg_yyyy_zzzzzz_1, tg_yyyz_xxxxxx_0, tg_yyyz_xxxxxx_1, tg_yyyz_xxxxxy_0, \
                                         tg_yyyz_xxxxxy_1, tg_yyyz_xxxxxz_0, tg_yyyz_xxxxxz_1, tg_yyyz_xxxxyy_0, \
                                         tg_yyyz_xxxxyy_1, tg_yyyz_xxxxyz_0, tg_yyyz_xxxxyz_1, tg_yyyz_xxxxzz_0, \
                                         tg_yyyz_xxxxzz_1, tg_yyyz_xxxyyy_0, tg_yyyz_xxxyyy_1, tg_yyyz_xxxyyz_0, \
                                         tg_yyyz_xxxyyz_1, tg_yyyz_xxxyzz_0, tg_yyyz_xxxyzz_1, tg_yyyz_xxxzzz_0, \
                                         tg_yyyz_xxxzzz_1, tg_yyyz_xxyyyy_0, tg_yyyz_xxyyyy_1, tg_yyyz_xxyyyz_0, \
                                         tg_yyyz_xxyyyz_1, tg_yyyz_xxyyzz_0, tg_yyyz_xxyyzz_1, tg_yyyz_xxyzzz_0, \
                                         tg_yyyz_xxyzzz_1, tg_yyyz_xxzzzz_0, tg_yyyz_xxzzzz_1, tg_yyyz_xyyyyy_0, \
                                         tg_yyyz_xyyyyy_1, tg_yyyz_xyyyyz_0, tg_yyyz_xyyyyz_1, tg_yyyz_xyyyzz_0, \
                                         tg_yyyz_xyyyzz_1, tg_yyyz_xyyzzz_0, tg_yyyz_xyyzzz_1, tg_yyyz_xyzzzz_0, \
                                         tg_yyyz_xyzzzz_1, tg_yyyz_xzzzzz_0, tg_yyyz_xzzzzz_1, tg_yyyz_yyyyyy_0, \
                                         tg_yyyz_yyyyyy_1, tg_yyyz_yyyyyz_0, tg_yyyz_yyyyyz_1, tg_yyyz_yyyyzz_0, \
                                         tg_yyyz_yyyyzz_1, tg_yyyz_yyyzzz_0, tg_yyyz_yyyzzz_1, tg_yyyz_yyzzzz_0, \
                                         tg_yyyz_yyzzzz_1, tg_yyyz_yzzzzz_0, tg_yyyz_yzzzzz_1, tg_yyyz_zzzzzz_0, \
                                         tg_yyyz_zzzzzz_1, tg_yyzz_xxxxxx_0, tg_yyzz_xxxxxx_1, tg_yyzz_xxxxxy_0, \
                                         tg_yyzz_xxxxxy_1, tg_yyzz_xxxxxz_0, tg_yyzz_xxxxxz_1, tg_yyzz_xxxxyy_0, \
                                         tg_yyzz_xxxxyy_1, tg_yyzz_xxxxyz_0, tg_yyzz_xxxxyz_1, tg_yyzz_xxxxzz_0, \
                                         tg_yyzz_xxxxzz_1, tg_yyzz_xxxyyy_0, tg_yyzz_xxxyyy_1, tg_yyzz_xxxyyz_0, \
                                         tg_yyzz_xxxyyz_1, tg_yyzz_xxxyzz_0, tg_yyzz_xxxyzz_1, tg_yyzz_xxxzzz_0, \
                                         tg_yyzz_xxxzzz_1, tg_yyzz_xxyyyy_0, tg_yyzz_xxyyyy_1, tg_yyzz_xxyyyz_0, \
                                         tg_yyzz_xxyyyz_1, tg_yyzz_xxyyzz_0, tg_yyzz_xxyyzz_1, tg_yyzz_xxyzzz_0, \
                                         tg_yyzz_xxyzzz_1, tg_yyzz_xxzzzz_0, tg_yyzz_xxzzzz_1, tg_yyzz_xyyyyy_0, \
                                         tg_yyzz_xyyyyy_1, tg_yyzz_xyyyyz_0, tg_yyzz_xyyyyz_1, tg_yyzz_xyyyzz_0, \
                                         tg_yyzz_xyyyzz_1, tg_yyzz_xyyzzz_0, tg_yyzz_xyyzzz_1, tg_yyzz_xyzzzz_0, \
                                         tg_yyzz_xyzzzz_1, tg_yyzz_xzzzzz_0, tg_yyzz_xzzzzz_1, tg_yyzz_yyyyyy_0, \
                                         tg_yyzz_yyyyyy_1, tg_yyzz_yyyyyz_0, tg_yyzz_yyyyyz_1, tg_yyzz_yyyyzz_0, \
                                         tg_yyzz_yyyyzz_1, tg_yyzz_yyyzzz_0, tg_yyzz_yyyzzz_1, tg_yyzz_yyzzzz_0, \
                                         tg_yyzz_yyzzzz_1, tg_yyzz_yzzzzz_0, tg_yyzz_yzzzzz_1, tg_yyzz_zzzzzz_0, \
                                         tg_yyzz_zzzzzz_1, tg_yzzz_xxxxxx_0, tg_yzzz_xxxxxx_1, tg_yzzz_xxxxxy_0, \
                                         tg_yzzz_xxxxxy_1, tg_yzzz_xxxxxz_0, tg_yzzz_xxxxxz_1, tg_yzzz_xxxxyy_0, \
                                         tg_yzzz_xxxxyy_1, tg_yzzz_xxxxyz_0, tg_yzzz_xxxxyz_1, tg_yzzz_xxxxzz_0, \
                                         tg_yzzz_xxxxzz_1, tg_yzzz_xxxyyy_0, tg_yzzz_xxxyyy_1, tg_yzzz_xxxyyz_0, \
                                         tg_yzzz_xxxyyz_1, tg_yzzz_xxxyzz_0, tg_yzzz_xxxyzz_1, tg_yzzz_xxxzzz_0, \
                                         tg_yzzz_xxxzzz_1, tg_yzzz_xxyyyy_0, tg_yzzz_xxyyyy_1, tg_yzzz_xxyyyz_0, \
                                         tg_yzzz_xxyyyz_1, tg_yzzz_xxyyzz_0, tg_yzzz_xxyyzz_1, tg_yzzz_xxyzzz_0, \
                                         tg_yzzz_xxyzzz_1, tg_yzzz_xxzzzz_0, tg_yzzz_xxzzzz_1, tg_yzzz_xyyyyy_0, \
                                         tg_yzzz_xyyyyy_1, tg_yzzz_xyyyyz_0, tg_yzzz_xyyyyz_1, tg_yzzz_xyyyzz_0, \
                                         tg_yzzz_xyyyzz_1, tg_yzzz_xyyzzz_0, tg_yzzz_xyyzzz_1, tg_yzzz_xyzzzz_0, \
                                         tg_yzzz_xyzzzz_1, tg_yzzz_xzzzzz_0, tg_yzzz_xzzzzz_1, tg_yzzz_yyyyyy_0, \
                                         tg_yzzz_yyyyyy_1, tg_yzzz_yyyyyz_0, tg_yzzz_yyyyyz_1, tg_yzzz_yyyyzz_0, \
                                         tg_yzzz_yyyyzz_1, tg_yzzz_yyyzzz_0, tg_yzzz_yyyzzz_1, tg_yzzz_yyzzzz_0, \
                                         tg_yzzz_yyzzzz_1, tg_yzzz_yzzzzz_0, tg_yzzz_yzzzzz_1, tg_yzzz_zzzzzz_0, \
                                         tg_yzzz_zzzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxyyyy_xxzzzz_0[j] = pb_x * tg_xyyyy_xxzzzz_0[j] + fr * tg_xyyyy_xxzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyy_xxzzzz_0[j] - tg_yyyy_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyy_xzzzz_1[j];

                    tg_xxyyyy_xyyyyy_0[j] = pb_x * tg_xyyyy_xyyyyy_0[j] + fr * tg_xyyyy_xyyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyy_xyyyyy_0[j] - tg_yyyy_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyy_yyyyy_1[j];

                    tg_xxyyyy_xyyyyz_0[j] = pb_x * tg_xyyyy_xyyyyz_0[j] + fr * tg_xyyyy_xyyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyy_xyyyyz_0[j] - tg_yyyy_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyy_yyyyz_1[j];

                    tg_xxyyyy_xyyyzz_0[j] = pb_x * tg_xyyyy_xyyyzz_0[j] + fr * tg_xyyyy_xyyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyy_xyyyzz_0[j] - tg_yyyy_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyy_yyyzz_1[j];

                    tg_xxyyyy_xyyzzz_0[j] = pb_x * tg_xyyyy_xyyzzz_0[j] + fr * tg_xyyyy_xyyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyy_xyyzzz_0[j] - tg_yyyy_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyy_yyzzz_1[j];

                    tg_xxyyyy_xyzzzz_0[j] = pb_x * tg_xyyyy_xyzzzz_0[j] + fr * tg_xyyyy_xyzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyy_xyzzzz_0[j] - tg_yyyy_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyy_yzzzz_1[j];

                    tg_xxyyyy_xzzzzz_0[j] = pb_x * tg_xyyyy_xzzzzz_0[j] + fr * tg_xyyyy_xzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyy_xzzzzz_0[j] - tg_yyyy_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyy_zzzzz_1[j];

                    tg_xxyyyy_yyyyyy_0[j] = pb_x * tg_xyyyy_yyyyyy_0[j] + fr * tg_xyyyy_yyyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyy_yyyyyy_0[j] - tg_yyyy_yyyyyy_1[j] * fl1_fza);

                    tg_xxyyyy_yyyyyz_0[j] = pb_x * tg_xyyyy_yyyyyz_0[j] + fr * tg_xyyyy_yyyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyy_yyyyyz_0[j] - tg_yyyy_yyyyyz_1[j] * fl1_fza);

                    tg_xxyyyy_yyyyzz_0[j] = pb_x * tg_xyyyy_yyyyzz_0[j] + fr * tg_xyyyy_yyyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyy_yyyyzz_0[j] - tg_yyyy_yyyyzz_1[j] * fl1_fza);

                    tg_xxyyyy_yyyzzz_0[j] = pb_x * tg_xyyyy_yyyzzz_0[j] + fr * tg_xyyyy_yyyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyy_yyyzzz_0[j] - tg_yyyy_yyyzzz_1[j] * fl1_fza);

                    tg_xxyyyy_yyzzzz_0[j] = pb_x * tg_xyyyy_yyzzzz_0[j] + fr * tg_xyyyy_yyzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyy_yyzzzz_0[j] - tg_yyyy_yyzzzz_1[j] * fl1_fza);

                    tg_xxyyyy_yzzzzz_0[j] = pb_x * tg_xyyyy_yzzzzz_0[j] + fr * tg_xyyyy_yzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyy_yzzzzz_0[j] - tg_yyyy_yzzzzz_1[j] * fl1_fza);

                    tg_xxyyyy_zzzzzz_0[j] = pb_x * tg_xyyyy_zzzzzz_0[j] + fr * tg_xyyyy_zzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyy_zzzzzz_0[j] - tg_yyyy_zzzzzz_1[j] * fl1_fza);

                    tg_xxyyyz_xxxxxx_0[j] = pb_x * tg_xyyyz_xxxxxx_0[j] + fr * tg_xyyyz_xxxxxx_1[j] + 0.5 * fl1_fx * (tg_yyyz_xxxxxx_0[j] - tg_yyyz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xyyyz_xxxxx_1[j];

                    tg_xxyyyz_xxxxxy_0[j] = pb_x * tg_xyyyz_xxxxxy_0[j] + fr * tg_xyyyz_xxxxxy_1[j] + 0.5 * fl1_fx * (tg_yyyz_xxxxxy_0[j] - tg_yyyz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyyz_xxxxy_1[j];

                    tg_xxyyyz_xxxxxz_0[j] = pb_x * tg_xyyyz_xxxxxz_0[j] + fr * tg_xyyyz_xxxxxz_1[j] + 0.5 * fl1_fx * (tg_yyyz_xxxxxz_0[j] - tg_yyyz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyyz_xxxxz_1[j];

                    tg_xxyyyz_xxxxyy_0[j] = pb_x * tg_xyyyz_xxxxyy_0[j] + fr * tg_xyyyz_xxxxyy_1[j] + 0.5 * fl1_fx * (tg_yyyz_xxxxyy_0[j] - tg_yyyz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyz_xxxyy_1[j];

                    tg_xxyyyz_xxxxyz_0[j] = pb_x * tg_xyyyz_xxxxyz_0[j] + fr * tg_xyyyz_xxxxyz_1[j] + 0.5 * fl1_fx * (tg_yyyz_xxxxyz_0[j] - tg_yyyz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyz_xxxyz_1[j];

                    tg_xxyyyz_xxxxzz_0[j] = pb_x * tg_xyyyz_xxxxzz_0[j] + fr * tg_xyyyz_xxxxzz_1[j] + 0.5 * fl1_fx * (tg_yyyz_xxxxzz_0[j] - tg_yyyz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyz_xxxzz_1[j];

                    tg_xxyyyz_xxxyyy_0[j] = pb_x * tg_xyyyz_xxxyyy_0[j] + fr * tg_xyyyz_xxxyyy_1[j] + 0.5 * fl1_fx * (tg_yyyz_xxxyyy_0[j] - tg_yyyz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyz_xxyyy_1[j];

                    tg_xxyyyz_xxxyyz_0[j] = pb_x * tg_xyyyz_xxxyyz_0[j] + fr * tg_xyyyz_xxxyyz_1[j] + 0.5 * fl1_fx * (tg_yyyz_xxxyyz_0[j] - tg_yyyz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyz_xxyyz_1[j];

                    tg_xxyyyz_xxxyzz_0[j] = pb_x * tg_xyyyz_xxxyzz_0[j] + fr * tg_xyyyz_xxxyzz_1[j] + 0.5 * fl1_fx * (tg_yyyz_xxxyzz_0[j] - tg_yyyz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyz_xxyzz_1[j];

                    tg_xxyyyz_xxxzzz_0[j] = pb_x * tg_xyyyz_xxxzzz_0[j] + fr * tg_xyyyz_xxxzzz_1[j] + 0.5 * fl1_fx * (tg_yyyz_xxxzzz_0[j] - tg_yyyz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyz_xxzzz_1[j];

                    tg_xxyyyz_xxyyyy_0[j] = pb_x * tg_xyyyz_xxyyyy_0[j] + fr * tg_xyyyz_xxyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyz_xxyyyy_0[j] - tg_yyyz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyyyz_xyyyy_1[j];

                    tg_xxyyyz_xxyyyz_0[j] = pb_x * tg_xyyyz_xxyyyz_0[j] + fr * tg_xyyyz_xxyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyz_xxyyyz_0[j] - tg_yyyz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyz_xyyyz_1[j];

                    tg_xxyyyz_xxyyzz_0[j] = pb_x * tg_xyyyz_xxyyzz_0[j] + fr * tg_xyyyz_xxyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyz_xxyyzz_0[j] - tg_yyyz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyz_xyyzz_1[j];

                    tg_xxyyyz_xxyzzz_0[j] = pb_x * tg_xyyyz_xxyzzz_0[j] + fr * tg_xyyyz_xxyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyz_xxyzzz_0[j] - tg_yyyz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyz_xyzzz_1[j];

                    tg_xxyyyz_xxzzzz_0[j] = pb_x * tg_xyyyz_xxzzzz_0[j] + fr * tg_xyyyz_xxzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyz_xxzzzz_0[j] - tg_yyyz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyz_xzzzz_1[j];

                    tg_xxyyyz_xyyyyy_0[j] = pb_x * tg_xyyyz_xyyyyy_0[j] + fr * tg_xyyyz_xyyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyz_xyyyyy_0[j] - tg_yyyz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyz_yyyyy_1[j];

                    tg_xxyyyz_xyyyyz_0[j] = pb_x * tg_xyyyz_xyyyyz_0[j] + fr * tg_xyyyz_xyyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyz_xyyyyz_0[j] - tg_yyyz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyz_yyyyz_1[j];

                    tg_xxyyyz_xyyyzz_0[j] = pb_x * tg_xyyyz_xyyyzz_0[j] + fr * tg_xyyyz_xyyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyz_xyyyzz_0[j] - tg_yyyz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyz_yyyzz_1[j];

                    tg_xxyyyz_xyyzzz_0[j] = pb_x * tg_xyyyz_xyyzzz_0[j] + fr * tg_xyyyz_xyyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyz_xyyzzz_0[j] - tg_yyyz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyz_yyzzz_1[j];

                    tg_xxyyyz_xyzzzz_0[j] = pb_x * tg_xyyyz_xyzzzz_0[j] + fr * tg_xyyyz_xyzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyz_xyzzzz_0[j] - tg_yyyz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyz_yzzzz_1[j];

                    tg_xxyyyz_xzzzzz_0[j] = pb_x * tg_xyyyz_xzzzzz_0[j] + fr * tg_xyyyz_xzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyz_xzzzzz_0[j] - tg_yyyz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyz_zzzzz_1[j];

                    tg_xxyyyz_yyyyyy_0[j] = pb_x * tg_xyyyz_yyyyyy_0[j] + fr * tg_xyyyz_yyyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyz_yyyyyy_0[j] - tg_yyyz_yyyyyy_1[j] * fl1_fza);

                    tg_xxyyyz_yyyyyz_0[j] = pb_x * tg_xyyyz_yyyyyz_0[j] + fr * tg_xyyyz_yyyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyz_yyyyyz_0[j] - tg_yyyz_yyyyyz_1[j] * fl1_fza);

                    tg_xxyyyz_yyyyzz_0[j] = pb_x * tg_xyyyz_yyyyzz_0[j] + fr * tg_xyyyz_yyyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyz_yyyyzz_0[j] - tg_yyyz_yyyyzz_1[j] * fl1_fza);

                    tg_xxyyyz_yyyzzz_0[j] = pb_x * tg_xyyyz_yyyzzz_0[j] + fr * tg_xyyyz_yyyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyz_yyyzzz_0[j] - tg_yyyz_yyyzzz_1[j] * fl1_fza);

                    tg_xxyyyz_yyzzzz_0[j] = pb_x * tg_xyyyz_yyzzzz_0[j] + fr * tg_xyyyz_yyzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyz_yyzzzz_0[j] - tg_yyyz_yyzzzz_1[j] * fl1_fza);

                    tg_xxyyyz_yzzzzz_0[j] = pb_x * tg_xyyyz_yzzzzz_0[j] + fr * tg_xyyyz_yzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyz_yzzzzz_0[j] - tg_yyyz_yzzzzz_1[j] * fl1_fza);

                    tg_xxyyyz_zzzzzz_0[j] = pb_x * tg_xyyyz_zzzzzz_0[j] + fr * tg_xyyyz_zzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyz_zzzzzz_0[j] - tg_yyyz_zzzzzz_1[j] * fl1_fza);

                    tg_xxyyzz_xxxxxx_0[j] = pb_x * tg_xyyzz_xxxxxx_0[j] + fr * tg_xyyzz_xxxxxx_1[j] + 0.5 * fl1_fx * (tg_yyzz_xxxxxx_0[j] - tg_yyzz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xyyzz_xxxxx_1[j];

                    tg_xxyyzz_xxxxxy_0[j] = pb_x * tg_xyyzz_xxxxxy_0[j] + fr * tg_xyyzz_xxxxxy_1[j] + 0.5 * fl1_fx * (tg_yyzz_xxxxxy_0[j] - tg_yyzz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyzz_xxxxy_1[j];

                    tg_xxyyzz_xxxxxz_0[j] = pb_x * tg_xyyzz_xxxxxz_0[j] + fr * tg_xyyzz_xxxxxz_1[j] + 0.5 * fl1_fx * (tg_yyzz_xxxxxz_0[j] - tg_yyzz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyzz_xxxxz_1[j];

                    tg_xxyyzz_xxxxyy_0[j] = pb_x * tg_xyyzz_xxxxyy_0[j] + fr * tg_xyyzz_xxxxyy_1[j] + 0.5 * fl1_fx * (tg_yyzz_xxxxyy_0[j] - tg_yyzz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyzz_xxxyy_1[j];

                    tg_xxyyzz_xxxxyz_0[j] = pb_x * tg_xyyzz_xxxxyz_0[j] + fr * tg_xyyzz_xxxxyz_1[j] + 0.5 * fl1_fx * (tg_yyzz_xxxxyz_0[j] - tg_yyzz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyzz_xxxyz_1[j];

                    tg_xxyyzz_xxxxzz_0[j] = pb_x * tg_xyyzz_xxxxzz_0[j] + fr * tg_xyyzz_xxxxzz_1[j] + 0.5 * fl1_fx * (tg_yyzz_xxxxzz_0[j] - tg_yyzz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyzz_xxxzz_1[j];

                    tg_xxyyzz_xxxyyy_0[j] = pb_x * tg_xyyzz_xxxyyy_0[j] + fr * tg_xyyzz_xxxyyy_1[j] + 0.5 * fl1_fx * (tg_yyzz_xxxyyy_0[j] - tg_yyzz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyzz_xxyyy_1[j];

                    tg_xxyyzz_xxxyyz_0[j] = pb_x * tg_xyyzz_xxxyyz_0[j] + fr * tg_xyyzz_xxxyyz_1[j] + 0.5 * fl1_fx * (tg_yyzz_xxxyyz_0[j] - tg_yyzz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyzz_xxyyz_1[j];

                    tg_xxyyzz_xxxyzz_0[j] = pb_x * tg_xyyzz_xxxyzz_0[j] + fr * tg_xyyzz_xxxyzz_1[j] + 0.5 * fl1_fx * (tg_yyzz_xxxyzz_0[j] - tg_yyzz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyzz_xxyzz_1[j];

                    tg_xxyyzz_xxxzzz_0[j] = pb_x * tg_xyyzz_xxxzzz_0[j] + fr * tg_xyyzz_xxxzzz_1[j] + 0.5 * fl1_fx * (tg_yyzz_xxxzzz_0[j] - tg_yyzz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyzz_xxzzz_1[j];

                    tg_xxyyzz_xxyyyy_0[j] = pb_x * tg_xyyzz_xxyyyy_0[j] + fr * tg_xyyzz_xxyyyy_1[j] + 0.5 * fl1_fx * (tg_yyzz_xxyyyy_0[j] - tg_yyzz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyyzz_xyyyy_1[j];

                    tg_xxyyzz_xxyyyz_0[j] = pb_x * tg_xyyzz_xxyyyz_0[j] + fr * tg_xyyzz_xxyyyz_1[j] + 0.5 * fl1_fx * (tg_yyzz_xxyyyz_0[j] - tg_yyzz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyyzz_xyyyz_1[j];

                    tg_xxyyzz_xxyyzz_0[j] = pb_x * tg_xyyzz_xxyyzz_0[j] + fr * tg_xyyzz_xxyyzz_1[j] + 0.5 * fl1_fx * (tg_yyzz_xxyyzz_0[j] - tg_yyzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyzz_xyyzz_1[j];

                    tg_xxyyzz_xxyzzz_0[j] = pb_x * tg_xyyzz_xxyzzz_0[j] + fr * tg_xyyzz_xxyzzz_1[j] + 0.5 * fl1_fx * (tg_yyzz_xxyzzz_0[j] - tg_yyzz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyzz_xyzzz_1[j];

                    tg_xxyyzz_xxzzzz_0[j] = pb_x * tg_xyyzz_xxzzzz_0[j] + fr * tg_xyyzz_xxzzzz_1[j] + 0.5 * fl1_fx * (tg_yyzz_xxzzzz_0[j] - tg_yyzz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyzz_xzzzz_1[j];

                    tg_xxyyzz_xyyyyy_0[j] = pb_x * tg_xyyzz_xyyyyy_0[j] + fr * tg_xyyzz_xyyyyy_1[j] + 0.5 * fl1_fx * (tg_yyzz_xyyyyy_0[j] - tg_yyzz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzz_yyyyy_1[j];

                    tg_xxyyzz_xyyyyz_0[j] = pb_x * tg_xyyzz_xyyyyz_0[j] + fr * tg_xyyzz_xyyyyz_1[j] + 0.5 * fl1_fx * (tg_yyzz_xyyyyz_0[j] - tg_yyzz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzz_yyyyz_1[j];

                    tg_xxyyzz_xyyyzz_0[j] = pb_x * tg_xyyzz_xyyyzz_0[j] + fr * tg_xyyzz_xyyyzz_1[j] + 0.5 * fl1_fx * (tg_yyzz_xyyyzz_0[j] - tg_yyzz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzz_yyyzz_1[j];

                    tg_xxyyzz_xyyzzz_0[j] = pb_x * tg_xyyzz_xyyzzz_0[j] + fr * tg_xyyzz_xyyzzz_1[j] + 0.5 * fl1_fx * (tg_yyzz_xyyzzz_0[j] - tg_yyzz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzz_yyzzz_1[j];

                    tg_xxyyzz_xyzzzz_0[j] = pb_x * tg_xyyzz_xyzzzz_0[j] + fr * tg_xyyzz_xyzzzz_1[j] + 0.5 * fl1_fx * (tg_yyzz_xyzzzz_0[j] - tg_yyzz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzz_yzzzz_1[j];

                    tg_xxyyzz_xzzzzz_0[j] = pb_x * tg_xyyzz_xzzzzz_0[j] + fr * tg_xyyzz_xzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyzz_xzzzzz_0[j] - tg_yyzz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzz_zzzzz_1[j];

                    tg_xxyyzz_yyyyyy_0[j] = pb_x * tg_xyyzz_yyyyyy_0[j] + fr * tg_xyyzz_yyyyyy_1[j] + 0.5 * fl1_fx * (tg_yyzz_yyyyyy_0[j] - tg_yyzz_yyyyyy_1[j] * fl1_fza);

                    tg_xxyyzz_yyyyyz_0[j] = pb_x * tg_xyyzz_yyyyyz_0[j] + fr * tg_xyyzz_yyyyyz_1[j] + 0.5 * fl1_fx * (tg_yyzz_yyyyyz_0[j] - tg_yyzz_yyyyyz_1[j] * fl1_fza);

                    tg_xxyyzz_yyyyzz_0[j] = pb_x * tg_xyyzz_yyyyzz_0[j] + fr * tg_xyyzz_yyyyzz_1[j] + 0.5 * fl1_fx * (tg_yyzz_yyyyzz_0[j] - tg_yyzz_yyyyzz_1[j] * fl1_fza);

                    tg_xxyyzz_yyyzzz_0[j] = pb_x * tg_xyyzz_yyyzzz_0[j] + fr * tg_xyyzz_yyyzzz_1[j] + 0.5 * fl1_fx * (tg_yyzz_yyyzzz_0[j] - tg_yyzz_yyyzzz_1[j] * fl1_fza);

                    tg_xxyyzz_yyzzzz_0[j] = pb_x * tg_xyyzz_yyzzzz_0[j] + fr * tg_xyyzz_yyzzzz_1[j] + 0.5 * fl1_fx * (tg_yyzz_yyzzzz_0[j] - tg_yyzz_yyzzzz_1[j] * fl1_fza);

                    tg_xxyyzz_yzzzzz_0[j] = pb_x * tg_xyyzz_yzzzzz_0[j] + fr * tg_xyyzz_yzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyzz_yzzzzz_0[j] - tg_yyzz_yzzzzz_1[j] * fl1_fza);

                    tg_xxyyzz_zzzzzz_0[j] = pb_x * tg_xyyzz_zzzzzz_0[j] + fr * tg_xyyzz_zzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyzz_zzzzzz_0[j] - tg_yyzz_zzzzzz_1[j] * fl1_fza);

                    tg_xxyzzz_xxxxxx_0[j] = pb_x * tg_xyzzz_xxxxxx_0[j] + fr * tg_xyzzz_xxxxxx_1[j] + 0.5 * fl1_fx * (tg_yzzz_xxxxxx_0[j] - tg_yzzz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xyzzz_xxxxx_1[j];

                    tg_xxyzzz_xxxxxy_0[j] = pb_x * tg_xyzzz_xxxxxy_0[j] + fr * tg_xyzzz_xxxxxy_1[j] + 0.5 * fl1_fx * (tg_yzzz_xxxxxy_0[j] - tg_yzzz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyzzz_xxxxy_1[j];

                    tg_xxyzzz_xxxxxz_0[j] = pb_x * tg_xyzzz_xxxxxz_0[j] + fr * tg_xyzzz_xxxxxz_1[j] + 0.5 * fl1_fx * (tg_yzzz_xxxxxz_0[j] - tg_yzzz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyzzz_xxxxz_1[j];

                    tg_xxyzzz_xxxxyy_0[j] = pb_x * tg_xyzzz_xxxxyy_0[j] + fr * tg_xyzzz_xxxxyy_1[j] + 0.5 * fl1_fx * (tg_yzzz_xxxxyy_0[j] - tg_yzzz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyzzz_xxxyy_1[j];

                    tg_xxyzzz_xxxxyz_0[j] = pb_x * tg_xyzzz_xxxxyz_0[j] + fr * tg_xyzzz_xxxxyz_1[j] + 0.5 * fl1_fx * (tg_yzzz_xxxxyz_0[j] - tg_yzzz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyzzz_xxxyz_1[j];

                    tg_xxyzzz_xxxxzz_0[j] = pb_x * tg_xyzzz_xxxxzz_0[j] + fr * tg_xyzzz_xxxxzz_1[j] + 0.5 * fl1_fx * (tg_yzzz_xxxxzz_0[j] - tg_yzzz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyzzz_xxxzz_1[j];

                    tg_xxyzzz_xxxyyy_0[j] = pb_x * tg_xyzzz_xxxyyy_0[j] + fr * tg_xyzzz_xxxyyy_1[j] + 0.5 * fl1_fx * (tg_yzzz_xxxyyy_0[j] - tg_yzzz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzzz_xxyyy_1[j];

                    tg_xxyzzz_xxxyyz_0[j] = pb_x * tg_xyzzz_xxxyyz_0[j] + fr * tg_xyzzz_xxxyyz_1[j] + 0.5 * fl1_fx * (tg_yzzz_xxxyyz_0[j] - tg_yzzz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzzz_xxyyz_1[j];

                    tg_xxyzzz_xxxyzz_0[j] = pb_x * tg_xyzzz_xxxyzz_0[j] + fr * tg_xyzzz_xxxyzz_1[j] + 0.5 * fl1_fx * (tg_yzzz_xxxyzz_0[j] - tg_yzzz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzzz_xxyzz_1[j];

                    tg_xxyzzz_xxxzzz_0[j] = pb_x * tg_xyzzz_xxxzzz_0[j] + fr * tg_xyzzz_xxxzzz_1[j] + 0.5 * fl1_fx * (tg_yzzz_xxxzzz_0[j] - tg_yzzz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzzz_xxzzz_1[j];

                    tg_xxyzzz_xxyyyy_0[j] = pb_x * tg_xyzzz_xxyyyy_0[j] + fr * tg_xyzzz_xxyyyy_1[j] + 0.5 * fl1_fx * (tg_yzzz_xxyyyy_0[j] - tg_yzzz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyzzz_xyyyy_1[j];

                    tg_xxyzzz_xxyyyz_0[j] = pb_x * tg_xyzzz_xxyyyz_0[j] + fr * tg_xyzzz_xxyyyz_1[j] + 0.5 * fl1_fx * (tg_yzzz_xxyyyz_0[j] - tg_yzzz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyzzz_xyyyz_1[j];

                    tg_xxyzzz_xxyyzz_0[j] = pb_x * tg_xyzzz_xxyyzz_0[j] + fr * tg_xyzzz_xxyyzz_1[j] + 0.5 * fl1_fx * (tg_yzzz_xxyyzz_0[j] - tg_yzzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyzzz_xyyzz_1[j];

                    tg_xxyzzz_xxyzzz_0[j] = pb_x * tg_xyzzz_xxyzzz_0[j] + fr * tg_xyzzz_xxyzzz_1[j] + 0.5 * fl1_fx * (tg_yzzz_xxyzzz_0[j] - tg_yzzz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyzzz_xyzzz_1[j];

                    tg_xxyzzz_xxzzzz_0[j] = pb_x * tg_xyzzz_xxzzzz_0[j] + fr * tg_xyzzz_xxzzzz_1[j] + 0.5 * fl1_fx * (tg_yzzz_xxzzzz_0[j] - tg_yzzz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyzzz_xzzzz_1[j];

                    tg_xxyzzz_xyyyyy_0[j] = pb_x * tg_xyzzz_xyyyyy_0[j] + fr * tg_xyzzz_xyyyyy_1[j] + 0.5 * fl1_fx * (tg_yzzz_xyyyyy_0[j] - tg_yzzz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzz_yyyyy_1[j];

                    tg_xxyzzz_xyyyyz_0[j] = pb_x * tg_xyzzz_xyyyyz_0[j] + fr * tg_xyzzz_xyyyyz_1[j] + 0.5 * fl1_fx * (tg_yzzz_xyyyyz_0[j] - tg_yzzz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzz_yyyyz_1[j];

                    tg_xxyzzz_xyyyzz_0[j] = pb_x * tg_xyzzz_xyyyzz_0[j] + fr * tg_xyzzz_xyyyzz_1[j] + 0.5 * fl1_fx * (tg_yzzz_xyyyzz_0[j] - tg_yzzz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzz_yyyzz_1[j];

                    tg_xxyzzz_xyyzzz_0[j] = pb_x * tg_xyzzz_xyyzzz_0[j] + fr * tg_xyzzz_xyyzzz_1[j] + 0.5 * fl1_fx * (tg_yzzz_xyyzzz_0[j] - tg_yzzz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzz_yyzzz_1[j];

                    tg_xxyzzz_xyzzzz_0[j] = pb_x * tg_xyzzz_xyzzzz_0[j] + fr * tg_xyzzz_xyzzzz_1[j] + 0.5 * fl1_fx * (tg_yzzz_xyzzzz_0[j] - tg_yzzz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzz_yzzzz_1[j];

                    tg_xxyzzz_xzzzzz_0[j] = pb_x * tg_xyzzz_xzzzzz_0[j] + fr * tg_xyzzz_xzzzzz_1[j] + 0.5 * fl1_fx * (tg_yzzz_xzzzzz_0[j] - tg_yzzz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzz_zzzzz_1[j];

                    tg_xxyzzz_yyyyyy_0[j] = pb_x * tg_xyzzz_yyyyyy_0[j] + fr * tg_xyzzz_yyyyyy_1[j] + 0.5 * fl1_fx * (tg_yzzz_yyyyyy_0[j] - tg_yzzz_yyyyyy_1[j] * fl1_fza);

                    tg_xxyzzz_yyyyyz_0[j] = pb_x * tg_xyzzz_yyyyyz_0[j] + fr * tg_xyzzz_yyyyyz_1[j] + 0.5 * fl1_fx * (tg_yzzz_yyyyyz_0[j] - tg_yzzz_yyyyyz_1[j] * fl1_fza);

                    tg_xxyzzz_yyyyzz_0[j] = pb_x * tg_xyzzz_yyyyzz_0[j] + fr * tg_xyzzz_yyyyzz_1[j] + 0.5 * fl1_fx * (tg_yzzz_yyyyzz_0[j] - tg_yzzz_yyyyzz_1[j] * fl1_fza);

                    tg_xxyzzz_yyyzzz_0[j] = pb_x * tg_xyzzz_yyyzzz_0[j] + fr * tg_xyzzz_yyyzzz_1[j] + 0.5 * fl1_fx * (tg_yzzz_yyyzzz_0[j] - tg_yzzz_yyyzzz_1[j] * fl1_fza);

                    tg_xxyzzz_yyzzzz_0[j] = pb_x * tg_xyzzz_yyzzzz_0[j] + fr * tg_xyzzz_yyzzzz_1[j] + 0.5 * fl1_fx * (tg_yzzz_yyzzzz_0[j] - tg_yzzz_yyzzzz_1[j] * fl1_fza);

                    tg_xxyzzz_yzzzzz_0[j] = pb_x * tg_xyzzz_yzzzzz_0[j] + fr * tg_xyzzz_yzzzzz_1[j] + 0.5 * fl1_fx * (tg_yzzz_yzzzzz_0[j] - tg_yzzz_yzzzzz_1[j] * fl1_fza);

                    tg_xxyzzz_zzzzzz_0[j] = pb_x * tg_xyzzz_zzzzzz_0[j] + fr * tg_xyzzz_zzzzzz_1[j] + 0.5 * fl1_fx * (tg_yzzz_zzzzzz_0[j] - tg_yzzz_zzzzzz_1[j] * fl1_fza);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSISI_392_490(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (392,490)

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
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_6_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_6_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_5_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
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

                auto tg_xzzzz_xxxxxx_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 392); 

                auto tg_xzzzz_xxxxxy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 393); 

                auto tg_xzzzz_xxxxxz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 394); 

                auto tg_xzzzz_xxxxyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 395); 

                auto tg_xzzzz_xxxxyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 396); 

                auto tg_xzzzz_xxxxzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 397); 

                auto tg_xzzzz_xxxyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 398); 

                auto tg_xzzzz_xxxyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 399); 

                auto tg_xzzzz_xxxyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 400); 

                auto tg_xzzzz_xxxzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 401); 

                auto tg_xzzzz_xxyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 402); 

                auto tg_xzzzz_xxyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 403); 

                auto tg_xzzzz_xxyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 404); 

                auto tg_xzzzz_xxyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 405); 

                auto tg_xzzzz_xxzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 406); 

                auto tg_xzzzz_xyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 407); 

                auto tg_xzzzz_xyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 408); 

                auto tg_xzzzz_xyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 409); 

                auto tg_xzzzz_xyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 410); 

                auto tg_xzzzz_xyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 411); 

                auto tg_xzzzz_xzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 412); 

                auto tg_xzzzz_yyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 413); 

                auto tg_xzzzz_yyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 414); 

                auto tg_xzzzz_yyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 415); 

                auto tg_xzzzz_yyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 416); 

                auto tg_xzzzz_yyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 417); 

                auto tg_xzzzz_yzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 418); 

                auto tg_xzzzz_zzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 419); 

                auto tg_yyyyy_xxxxxx_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 420); 

                auto tg_yyyyy_xxxxxy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 421); 

                auto tg_yyyyy_xxxxxz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 422); 

                auto tg_yyyyy_xxxxyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 423); 

                auto tg_yyyyy_xxxxyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 424); 

                auto tg_yyyyy_xxxxzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 425); 

                auto tg_yyyyy_xxxyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 426); 

                auto tg_yyyyy_xxxyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 427); 

                auto tg_yyyyy_xxxyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 428); 

                auto tg_yyyyy_xxxzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 429); 

                auto tg_yyyyy_xxyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 430); 

                auto tg_yyyyy_xxyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 431); 

                auto tg_yyyyy_xxyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 432); 

                auto tg_yyyyy_xxyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 433); 

                auto tg_yyyyy_xxzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 434); 

                auto tg_yyyyy_xyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 435); 

                auto tg_yyyyy_xyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 436); 

                auto tg_yyyyy_xyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 437); 

                auto tg_yyyyy_xyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 438); 

                auto tg_yyyyy_xyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 439); 

                auto tg_yyyyy_xzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 440); 

                auto tg_yyyyy_yyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 441); 

                auto tg_yyyyy_yyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 442); 

                auto tg_yyyyy_yyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 443); 

                auto tg_yyyyy_yyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 444); 

                auto tg_yyyyy_yyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 445); 

                auto tg_yyyyy_yzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 446); 

                auto tg_yyyyy_zzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 447); 

                auto tg_yyyyz_xxxxxx_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 448); 

                auto tg_yyyyz_xxxxxy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 449); 

                auto tg_yyyyz_xxxxxz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 450); 

                auto tg_yyyyz_xxxxyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 451); 

                auto tg_yyyyz_xxxxyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 452); 

                auto tg_yyyyz_xxxxzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 453); 

                auto tg_yyyyz_xxxyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 454); 

                auto tg_yyyyz_xxxyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 455); 

                auto tg_yyyyz_xxxyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 456); 

                auto tg_yyyyz_xxxzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 457); 

                auto tg_yyyyz_xxyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 458); 

                auto tg_yyyyz_xxyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 459); 

                auto tg_yyyyz_xxyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 460); 

                auto tg_yyyyz_xxyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 461); 

                auto tg_yyyyz_xxzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 462); 

                auto tg_yyyyz_xyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 463); 

                auto tg_yyyyz_xyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 464); 

                auto tg_yyyyz_xyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 465); 

                auto tg_yyyyz_xyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 466); 

                auto tg_yyyyz_xyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 467); 

                auto tg_yyyyz_xzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 468); 

                auto tg_yyyyz_yyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 469); 

                auto tg_yyyyz_yyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 470); 

                auto tg_yyyyz_yyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 471); 

                auto tg_yyyyz_yyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 472); 

                auto tg_yyyyz_yyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 473); 

                auto tg_yyyyz_yzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 474); 

                auto tg_yyyyz_zzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 475); 

                auto tg_yyyzz_xxxxxx_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 476); 

                auto tg_yyyzz_xxxxxy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 477); 

                auto tg_yyyzz_xxxxxz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 478); 

                auto tg_yyyzz_xxxxyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 479); 

                auto tg_yyyzz_xxxxyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 480); 

                auto tg_yyyzz_xxxxzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 481); 

                auto tg_yyyzz_xxxyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 482); 

                auto tg_yyyzz_xxxyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 483); 

                auto tg_yyyzz_xxxyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 484); 

                auto tg_yyyzz_xxxzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 485); 

                auto tg_yyyzz_xxyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 486); 

                auto tg_yyyzz_xxyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 487); 

                auto tg_yyyzz_xxyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 488); 

                auto tg_yyyzz_xxyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 489); 

                auto tg_xzzzz_xxxxxx_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 392); 

                auto tg_xzzzz_xxxxxy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 393); 

                auto tg_xzzzz_xxxxxz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 394); 

                auto tg_xzzzz_xxxxyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 395); 

                auto tg_xzzzz_xxxxyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 396); 

                auto tg_xzzzz_xxxxzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 397); 

                auto tg_xzzzz_xxxyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 398); 

                auto tg_xzzzz_xxxyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 399); 

                auto tg_xzzzz_xxxyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 400); 

                auto tg_xzzzz_xxxzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 401); 

                auto tg_xzzzz_xxyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 402); 

                auto tg_xzzzz_xxyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 403); 

                auto tg_xzzzz_xxyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 404); 

                auto tg_xzzzz_xxyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 405); 

                auto tg_xzzzz_xxzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 406); 

                auto tg_xzzzz_xyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 407); 

                auto tg_xzzzz_xyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 408); 

                auto tg_xzzzz_xyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 409); 

                auto tg_xzzzz_xyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 410); 

                auto tg_xzzzz_xyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 411); 

                auto tg_xzzzz_xzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 412); 

                auto tg_xzzzz_yyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 413); 

                auto tg_xzzzz_yyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 414); 

                auto tg_xzzzz_yyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 415); 

                auto tg_xzzzz_yyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 416); 

                auto tg_xzzzz_yyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 417); 

                auto tg_xzzzz_yzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 418); 

                auto tg_xzzzz_zzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 419); 

                auto tg_yyyyy_xxxxxx_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 420); 

                auto tg_yyyyy_xxxxxy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 421); 

                auto tg_yyyyy_xxxxxz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 422); 

                auto tg_yyyyy_xxxxyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 423); 

                auto tg_yyyyy_xxxxyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 424); 

                auto tg_yyyyy_xxxxzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 425); 

                auto tg_yyyyy_xxxyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 426); 

                auto tg_yyyyy_xxxyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 427); 

                auto tg_yyyyy_xxxyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 428); 

                auto tg_yyyyy_xxxzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 429); 

                auto tg_yyyyy_xxyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 430); 

                auto tg_yyyyy_xxyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 431); 

                auto tg_yyyyy_xxyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 432); 

                auto tg_yyyyy_xxyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 433); 

                auto tg_yyyyy_xxzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 434); 

                auto tg_yyyyy_xyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 435); 

                auto tg_yyyyy_xyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 436); 

                auto tg_yyyyy_xyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 437); 

                auto tg_yyyyy_xyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 438); 

                auto tg_yyyyy_xyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 439); 

                auto tg_yyyyy_xzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 440); 

                auto tg_yyyyy_yyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 441); 

                auto tg_yyyyy_yyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 442); 

                auto tg_yyyyy_yyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 443); 

                auto tg_yyyyy_yyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 444); 

                auto tg_yyyyy_yyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 445); 

                auto tg_yyyyy_yzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 446); 

                auto tg_yyyyy_zzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 447); 

                auto tg_yyyyz_xxxxxx_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 448); 

                auto tg_yyyyz_xxxxxy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 449); 

                auto tg_yyyyz_xxxxxz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 450); 

                auto tg_yyyyz_xxxxyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 451); 

                auto tg_yyyyz_xxxxyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 452); 

                auto tg_yyyyz_xxxxzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 453); 

                auto tg_yyyyz_xxxyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 454); 

                auto tg_yyyyz_xxxyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 455); 

                auto tg_yyyyz_xxxyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 456); 

                auto tg_yyyyz_xxxzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 457); 

                auto tg_yyyyz_xxyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 458); 

                auto tg_yyyyz_xxyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 459); 

                auto tg_yyyyz_xxyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 460); 

                auto tg_yyyyz_xxyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 461); 

                auto tg_yyyyz_xxzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 462); 

                auto tg_yyyyz_xyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 463); 

                auto tg_yyyyz_xyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 464); 

                auto tg_yyyyz_xyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 465); 

                auto tg_yyyyz_xyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 466); 

                auto tg_yyyyz_xyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 467); 

                auto tg_yyyyz_xzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 468); 

                auto tg_yyyyz_yyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 469); 

                auto tg_yyyyz_yyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 470); 

                auto tg_yyyyz_yyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 471); 

                auto tg_yyyyz_yyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 472); 

                auto tg_yyyyz_yyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 473); 

                auto tg_yyyyz_yzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 474); 

                auto tg_yyyyz_zzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 475); 

                auto tg_yyyzz_xxxxxx_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 476); 

                auto tg_yyyzz_xxxxxy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 477); 

                auto tg_yyyzz_xxxxxz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 478); 

                auto tg_yyyzz_xxxxyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 479); 

                auto tg_yyyzz_xxxxyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 480); 

                auto tg_yyyzz_xxxxzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 481); 

                auto tg_yyyzz_xxxyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 482); 

                auto tg_yyyzz_xxxyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 483); 

                auto tg_yyyzz_xxxyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 484); 

                auto tg_yyyzz_xxxzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 485); 

                auto tg_yyyzz_xxyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 486); 

                auto tg_yyyzz_xxyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 487); 

                auto tg_yyyzz_xxyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 488); 

                auto tg_yyyzz_xxyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 489); 

                auto tg_zzzz_xxxxxx_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 392); 

                auto tg_zzzz_xxxxxy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 393); 

                auto tg_zzzz_xxxxxz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 394); 

                auto tg_zzzz_xxxxyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 395); 

                auto tg_zzzz_xxxxyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 396); 

                auto tg_zzzz_xxxxzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 397); 

                auto tg_zzzz_xxxyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 398); 

                auto tg_zzzz_xxxyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 399); 

                auto tg_zzzz_xxxyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 400); 

                auto tg_zzzz_xxxzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 401); 

                auto tg_zzzz_xxyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 402); 

                auto tg_zzzz_xxyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 403); 

                auto tg_zzzz_xxyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 404); 

                auto tg_zzzz_xxyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 405); 

                auto tg_zzzz_xxzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 406); 

                auto tg_zzzz_xyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 407); 

                auto tg_zzzz_xyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 408); 

                auto tg_zzzz_xyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 409); 

                auto tg_zzzz_xyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 410); 

                auto tg_zzzz_xyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 411); 

                auto tg_zzzz_xzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 412); 

                auto tg_zzzz_yyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 413); 

                auto tg_zzzz_yyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 414); 

                auto tg_zzzz_yyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 415); 

                auto tg_zzzz_yyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 416); 

                auto tg_zzzz_yyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 417); 

                auto tg_zzzz_yzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 418); 

                auto tg_zzzz_zzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 419); 

                auto tg_zzzz_xxxxxx_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 392); 

                auto tg_zzzz_xxxxxy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 393); 

                auto tg_zzzz_xxxxxz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 394); 

                auto tg_zzzz_xxxxyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 395); 

                auto tg_zzzz_xxxxyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 396); 

                auto tg_zzzz_xxxxzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 397); 

                auto tg_zzzz_xxxyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 398); 

                auto tg_zzzz_xxxyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 399); 

                auto tg_zzzz_xxxyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 400); 

                auto tg_zzzz_xxxzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 401); 

                auto tg_zzzz_xxyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 402); 

                auto tg_zzzz_xxyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 403); 

                auto tg_zzzz_xxyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 404); 

                auto tg_zzzz_xxyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 405); 

                auto tg_zzzz_xxzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 406); 

                auto tg_zzzz_xyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 407); 

                auto tg_zzzz_xyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 408); 

                auto tg_zzzz_xyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 409); 

                auto tg_zzzz_xyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 410); 

                auto tg_zzzz_xyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 411); 

                auto tg_zzzz_xzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 412); 

                auto tg_zzzz_yyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 413); 

                auto tg_zzzz_yyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 414); 

                auto tg_zzzz_yyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 415); 

                auto tg_zzzz_yyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 416); 

                auto tg_zzzz_yyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 417); 

                auto tg_zzzz_yzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 418); 

                auto tg_zzzz_zzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 419); 

                auto tg_xzzzz_xxxxx_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 294); 

                auto tg_xzzzz_xxxxy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 295); 

                auto tg_xzzzz_xxxxz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 296); 

                auto tg_xzzzz_xxxyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 297); 

                auto tg_xzzzz_xxxyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 298); 

                auto tg_xzzzz_xxxzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 299); 

                auto tg_xzzzz_xxyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 300); 

                auto tg_xzzzz_xxyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 301); 

                auto tg_xzzzz_xxyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 302); 

                auto tg_xzzzz_xxzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 303); 

                auto tg_xzzzz_xyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 304); 

                auto tg_xzzzz_xyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 305); 

                auto tg_xzzzz_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 306); 

                auto tg_xzzzz_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 307); 

                auto tg_xzzzz_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 308); 

                auto tg_xzzzz_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 309); 

                auto tg_xzzzz_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 310); 

                auto tg_xzzzz_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 311); 

                auto tg_xzzzz_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 312); 

                auto tg_xzzzz_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 313); 

                auto tg_xzzzz_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 314); 

                auto tg_yyyyy_xxxxx_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 315); 

                auto tg_yyyyy_xxxxy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 316); 

                auto tg_yyyyy_xxxxz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 317); 

                auto tg_yyyyy_xxxyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 318); 

                auto tg_yyyyy_xxxyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 319); 

                auto tg_yyyyy_xxxzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 320); 

                auto tg_yyyyy_xxyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 321); 

                auto tg_yyyyy_xxyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 322); 

                auto tg_yyyyy_xxyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 323); 

                auto tg_yyyyy_xxzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 324); 

                auto tg_yyyyy_xyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 325); 

                auto tg_yyyyy_xyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 326); 

                auto tg_yyyyy_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 327); 

                auto tg_yyyyy_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 328); 

                auto tg_yyyyy_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 329); 

                auto tg_yyyyy_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 330); 

                auto tg_yyyyy_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 331); 

                auto tg_yyyyy_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 332); 

                auto tg_yyyyy_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 333); 

                auto tg_yyyyy_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 334); 

                auto tg_yyyyy_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 335); 

                auto tg_yyyyz_xxxxx_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 336); 

                auto tg_yyyyz_xxxxy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 337); 

                auto tg_yyyyz_xxxxz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 338); 

                auto tg_yyyyz_xxxyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 339); 

                auto tg_yyyyz_xxxyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 340); 

                auto tg_yyyyz_xxxzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 341); 

                auto tg_yyyyz_xxyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 342); 

                auto tg_yyyyz_xxyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 343); 

                auto tg_yyyyz_xxyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 344); 

                auto tg_yyyyz_xxzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 345); 

                auto tg_yyyyz_xyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 346); 

                auto tg_yyyyz_xyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 347); 

                auto tg_yyyyz_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 348); 

                auto tg_yyyyz_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 349); 

                auto tg_yyyyz_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 350); 

                auto tg_yyyyz_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 351); 

                auto tg_yyyyz_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 352); 

                auto tg_yyyyz_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 353); 

                auto tg_yyyyz_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 354); 

                auto tg_yyyyz_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 355); 

                auto tg_yyyyz_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 356); 

                auto tg_yyyzz_xxxxx_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 357); 

                auto tg_yyyzz_xxxxy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 358); 

                auto tg_yyyzz_xxxxz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 359); 

                auto tg_yyyzz_xxxyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 360); 

                auto tg_yyyzz_xxxyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 361); 

                auto tg_yyyzz_xxxzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 362); 

                auto tg_yyyzz_xxyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 363); 

                auto tg_yyyzz_xxyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 364); 

                auto tg_yyyzz_xxyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 365); 

                auto tg_yyyzz_xxzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 366); 

                auto tg_yyyzz_xyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 367); 

                auto tg_yyyzz_xyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 368); 

                auto tg_yyyzz_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 369); 

                auto tg_yyyzz_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 370); 

                // set up pointers to integrals

                auto tg_xxzzzz_xxxxxx_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 392); 

                auto tg_xxzzzz_xxxxxy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 393); 

                auto tg_xxzzzz_xxxxxz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 394); 

                auto tg_xxzzzz_xxxxyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 395); 

                auto tg_xxzzzz_xxxxyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 396); 

                auto tg_xxzzzz_xxxxzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 397); 

                auto tg_xxzzzz_xxxyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 398); 

                auto tg_xxzzzz_xxxyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 399); 

                auto tg_xxzzzz_xxxyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 400); 

                auto tg_xxzzzz_xxxzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 401); 

                auto tg_xxzzzz_xxyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 402); 

                auto tg_xxzzzz_xxyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 403); 

                auto tg_xxzzzz_xxyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 404); 

                auto tg_xxzzzz_xxyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 405); 

                auto tg_xxzzzz_xxzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 406); 

                auto tg_xxzzzz_xyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 407); 

                auto tg_xxzzzz_xyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 408); 

                auto tg_xxzzzz_xyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 409); 

                auto tg_xxzzzz_xyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 410); 

                auto tg_xxzzzz_xyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 411); 

                auto tg_xxzzzz_xzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 412); 

                auto tg_xxzzzz_yyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 413); 

                auto tg_xxzzzz_yyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 414); 

                auto tg_xxzzzz_yyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 415); 

                auto tg_xxzzzz_yyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 416); 

                auto tg_xxzzzz_yyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 417); 

                auto tg_xxzzzz_yzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 418); 

                auto tg_xxzzzz_zzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 419); 

                auto tg_xyyyyy_xxxxxx_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 420); 

                auto tg_xyyyyy_xxxxxy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 421); 

                auto tg_xyyyyy_xxxxxz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 422); 

                auto tg_xyyyyy_xxxxyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 423); 

                auto tg_xyyyyy_xxxxyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 424); 

                auto tg_xyyyyy_xxxxzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 425); 

                auto tg_xyyyyy_xxxyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 426); 

                auto tg_xyyyyy_xxxyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 427); 

                auto tg_xyyyyy_xxxyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 428); 

                auto tg_xyyyyy_xxxzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 429); 

                auto tg_xyyyyy_xxyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 430); 

                auto tg_xyyyyy_xxyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 431); 

                auto tg_xyyyyy_xxyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 432); 

                auto tg_xyyyyy_xxyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 433); 

                auto tg_xyyyyy_xxzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 434); 

                auto tg_xyyyyy_xyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 435); 

                auto tg_xyyyyy_xyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 436); 

                auto tg_xyyyyy_xyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 437); 

                auto tg_xyyyyy_xyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 438); 

                auto tg_xyyyyy_xyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 439); 

                auto tg_xyyyyy_xzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 440); 

                auto tg_xyyyyy_yyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 441); 

                auto tg_xyyyyy_yyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 442); 

                auto tg_xyyyyy_yyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 443); 

                auto tg_xyyyyy_yyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 444); 

                auto tg_xyyyyy_yyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 445); 

                auto tg_xyyyyy_yzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 446); 

                auto tg_xyyyyy_zzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 447); 

                auto tg_xyyyyz_xxxxxx_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 448); 

                auto tg_xyyyyz_xxxxxy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 449); 

                auto tg_xyyyyz_xxxxxz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 450); 

                auto tg_xyyyyz_xxxxyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 451); 

                auto tg_xyyyyz_xxxxyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 452); 

                auto tg_xyyyyz_xxxxzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 453); 

                auto tg_xyyyyz_xxxyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 454); 

                auto tg_xyyyyz_xxxyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 455); 

                auto tg_xyyyyz_xxxyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 456); 

                auto tg_xyyyyz_xxxzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 457); 

                auto tg_xyyyyz_xxyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 458); 

                auto tg_xyyyyz_xxyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 459); 

                auto tg_xyyyyz_xxyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 460); 

                auto tg_xyyyyz_xxyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 461); 

                auto tg_xyyyyz_xxzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 462); 

                auto tg_xyyyyz_xyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 463); 

                auto tg_xyyyyz_xyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 464); 

                auto tg_xyyyyz_xyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 465); 

                auto tg_xyyyyz_xyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 466); 

                auto tg_xyyyyz_xyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 467); 

                auto tg_xyyyyz_xzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 468); 

                auto tg_xyyyyz_yyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 469); 

                auto tg_xyyyyz_yyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 470); 

                auto tg_xyyyyz_yyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 471); 

                auto tg_xyyyyz_yyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 472); 

                auto tg_xyyyyz_yyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 473); 

                auto tg_xyyyyz_yzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 474); 

                auto tg_xyyyyz_zzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 475); 

                auto tg_xyyyzz_xxxxxx_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 476); 

                auto tg_xyyyzz_xxxxxy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 477); 

                auto tg_xyyyzz_xxxxxz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 478); 

                auto tg_xyyyzz_xxxxyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 479); 

                auto tg_xyyyzz_xxxxyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 480); 

                auto tg_xyyyzz_xxxxzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 481); 

                auto tg_xyyyzz_xxxyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 482); 

                auto tg_xyyyzz_xxxyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 483); 

                auto tg_xyyyzz_xxxyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 484); 

                auto tg_xyyyzz_xxxzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 485); 

                auto tg_xyyyzz_xxyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 486); 

                auto tg_xyyyzz_xxyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 487); 

                auto tg_xyyyzz_xxyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 488); 

                auto tg_xyyyzz_xxyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 489); 

                // Batch of Integrals (392,490)

                #pragma omp simd aligned(fxn, fza, tg_xxzzzz_xxxxxx_0, tg_xxzzzz_xxxxxy_0, tg_xxzzzz_xxxxxz_0, \
                                         tg_xxzzzz_xxxxyy_0, tg_xxzzzz_xxxxyz_0, tg_xxzzzz_xxxxzz_0, tg_xxzzzz_xxxyyy_0, \
                                         tg_xxzzzz_xxxyyz_0, tg_xxzzzz_xxxyzz_0, tg_xxzzzz_xxxzzz_0, tg_xxzzzz_xxyyyy_0, \
                                         tg_xxzzzz_xxyyyz_0, tg_xxzzzz_xxyyzz_0, tg_xxzzzz_xxyzzz_0, tg_xxzzzz_xxzzzz_0, \
                                         tg_xxzzzz_xyyyyy_0, tg_xxzzzz_xyyyyz_0, tg_xxzzzz_xyyyzz_0, tg_xxzzzz_xyyzzz_0, \
                                         tg_xxzzzz_xyzzzz_0, tg_xxzzzz_xzzzzz_0, tg_xxzzzz_yyyyyy_0, tg_xxzzzz_yyyyyz_0, \
                                         tg_xxzzzz_yyyyzz_0, tg_xxzzzz_yyyzzz_0, tg_xxzzzz_yyzzzz_0, tg_xxzzzz_yzzzzz_0, \
                                         tg_xxzzzz_zzzzzz_0, tg_xyyyyy_xxxxxx_0, tg_xyyyyy_xxxxxy_0, tg_xyyyyy_xxxxxz_0, \
                                         tg_xyyyyy_xxxxyy_0, tg_xyyyyy_xxxxyz_0, tg_xyyyyy_xxxxzz_0, tg_xyyyyy_xxxyyy_0, \
                                         tg_xyyyyy_xxxyyz_0, tg_xyyyyy_xxxyzz_0, tg_xyyyyy_xxxzzz_0, tg_xyyyyy_xxyyyy_0, \
                                         tg_xyyyyy_xxyyyz_0, tg_xyyyyy_xxyyzz_0, tg_xyyyyy_xxyzzz_0, tg_xyyyyy_xxzzzz_0, \
                                         tg_xyyyyy_xyyyyy_0, tg_xyyyyy_xyyyyz_0, tg_xyyyyy_xyyyzz_0, tg_xyyyyy_xyyzzz_0, \
                                         tg_xyyyyy_xyzzzz_0, tg_xyyyyy_xzzzzz_0, tg_xyyyyy_yyyyyy_0, tg_xyyyyy_yyyyyz_0, \
                                         tg_xyyyyy_yyyyzz_0, tg_xyyyyy_yyyzzz_0, tg_xyyyyy_yyzzzz_0, tg_xyyyyy_yzzzzz_0, \
                                         tg_xyyyyy_zzzzzz_0, tg_xyyyyz_xxxxxx_0, tg_xyyyyz_xxxxxy_0, tg_xyyyyz_xxxxxz_0, \
                                         tg_xyyyyz_xxxxyy_0, tg_xyyyyz_xxxxyz_0, tg_xyyyyz_xxxxzz_0, tg_xyyyyz_xxxyyy_0, \
                                         tg_xyyyyz_xxxyyz_0, tg_xyyyyz_xxxyzz_0, tg_xyyyyz_xxxzzz_0, tg_xyyyyz_xxyyyy_0, \
                                         tg_xyyyyz_xxyyyz_0, tg_xyyyyz_xxyyzz_0, tg_xyyyyz_xxyzzz_0, tg_xyyyyz_xxzzzz_0, \
                                         tg_xyyyyz_xyyyyy_0, tg_xyyyyz_xyyyyz_0, tg_xyyyyz_xyyyzz_0, tg_xyyyyz_xyyzzz_0, \
                                         tg_xyyyyz_xyzzzz_0, tg_xyyyyz_xzzzzz_0, tg_xyyyyz_yyyyyy_0, tg_xyyyyz_yyyyyz_0, \
                                         tg_xyyyyz_yyyyzz_0, tg_xyyyyz_yyyzzz_0, tg_xyyyyz_yyzzzz_0, tg_xyyyyz_yzzzzz_0, \
                                         tg_xyyyyz_zzzzzz_0, tg_xyyyzz_xxxxxx_0, tg_xyyyzz_xxxxxy_0, tg_xyyyzz_xxxxxz_0, \
                                         tg_xyyyzz_xxxxyy_0, tg_xyyyzz_xxxxyz_0, tg_xyyyzz_xxxxzz_0, tg_xyyyzz_xxxyyy_0, \
                                         tg_xyyyzz_xxxyyz_0, tg_xyyyzz_xxxyzz_0, tg_xyyyzz_xxxzzz_0, tg_xyyyzz_xxyyyy_0, \
                                         tg_xyyyzz_xxyyyz_0, tg_xyyyzz_xxyyzz_0, tg_xyyyzz_xxyzzz_0, tg_xzzzz_xxxxx_1, \
                                         tg_xzzzz_xxxxxx_0, tg_xzzzz_xxxxxx_1, tg_xzzzz_xxxxxy_0, tg_xzzzz_xxxxxy_1, \
                                         tg_xzzzz_xxxxxz_0, tg_xzzzz_xxxxxz_1, tg_xzzzz_xxxxy_1, tg_xzzzz_xxxxyy_0, \
                                         tg_xzzzz_xxxxyy_1, tg_xzzzz_xxxxyz_0, tg_xzzzz_xxxxyz_1, tg_xzzzz_xxxxz_1, \
                                         tg_xzzzz_xxxxzz_0, tg_xzzzz_xxxxzz_1, tg_xzzzz_xxxyy_1, tg_xzzzz_xxxyyy_0, \
                                         tg_xzzzz_xxxyyy_1, tg_xzzzz_xxxyyz_0, tg_xzzzz_xxxyyz_1, tg_xzzzz_xxxyz_1, \
                                         tg_xzzzz_xxxyzz_0, tg_xzzzz_xxxyzz_1, tg_xzzzz_xxxzz_1, tg_xzzzz_xxxzzz_0, \
                                         tg_xzzzz_xxxzzz_1, tg_xzzzz_xxyyy_1, tg_xzzzz_xxyyyy_0, tg_xzzzz_xxyyyy_1, \
                                         tg_xzzzz_xxyyyz_0, tg_xzzzz_xxyyyz_1, tg_xzzzz_xxyyz_1, tg_xzzzz_xxyyzz_0, \
                                         tg_xzzzz_xxyyzz_1, tg_xzzzz_xxyzz_1, tg_xzzzz_xxyzzz_0, tg_xzzzz_xxyzzz_1, \
                                         tg_xzzzz_xxzzz_1, tg_xzzzz_xxzzzz_0, tg_xzzzz_xxzzzz_1, tg_xzzzz_xyyyy_1, \
                                         tg_xzzzz_xyyyyy_0, tg_xzzzz_xyyyyy_1, tg_xzzzz_xyyyyz_0, tg_xzzzz_xyyyyz_1, \
                                         tg_xzzzz_xyyyz_1, tg_xzzzz_xyyyzz_0, tg_xzzzz_xyyyzz_1, tg_xzzzz_xyyzz_1, \
                                         tg_xzzzz_xyyzzz_0, tg_xzzzz_xyyzzz_1, tg_xzzzz_xyzzz_1, tg_xzzzz_xyzzzz_0, \
                                         tg_xzzzz_xyzzzz_1, tg_xzzzz_xzzzz_1, tg_xzzzz_xzzzzz_0, tg_xzzzz_xzzzzz_1, \
                                         tg_xzzzz_yyyyy_1, tg_xzzzz_yyyyyy_0, tg_xzzzz_yyyyyy_1, tg_xzzzz_yyyyyz_0, \
                                         tg_xzzzz_yyyyyz_1, tg_xzzzz_yyyyz_1, tg_xzzzz_yyyyzz_0, tg_xzzzz_yyyyzz_1, \
                                         tg_xzzzz_yyyzz_1, tg_xzzzz_yyyzzz_0, tg_xzzzz_yyyzzz_1, tg_xzzzz_yyzzz_1, \
                                         tg_xzzzz_yyzzzz_0, tg_xzzzz_yyzzzz_1, tg_xzzzz_yzzzz_1, tg_xzzzz_yzzzzz_0, \
                                         tg_xzzzz_yzzzzz_1, tg_xzzzz_zzzzz_1, tg_xzzzz_zzzzzz_0, tg_xzzzz_zzzzzz_1, \
                                         tg_yyyyy_xxxxx_1, tg_yyyyy_xxxxxx_0, tg_yyyyy_xxxxxx_1, tg_yyyyy_xxxxxy_0, \
                                         tg_yyyyy_xxxxxy_1, tg_yyyyy_xxxxxz_0, tg_yyyyy_xxxxxz_1, tg_yyyyy_xxxxy_1, \
                                         tg_yyyyy_xxxxyy_0, tg_yyyyy_xxxxyy_1, tg_yyyyy_xxxxyz_0, tg_yyyyy_xxxxyz_1, \
                                         tg_yyyyy_xxxxz_1, tg_yyyyy_xxxxzz_0, tg_yyyyy_xxxxzz_1, tg_yyyyy_xxxyy_1, \
                                         tg_yyyyy_xxxyyy_0, tg_yyyyy_xxxyyy_1, tg_yyyyy_xxxyyz_0, tg_yyyyy_xxxyyz_1, \
                                         tg_yyyyy_xxxyz_1, tg_yyyyy_xxxyzz_0, tg_yyyyy_xxxyzz_1, tg_yyyyy_xxxzz_1, \
                                         tg_yyyyy_xxxzzz_0, tg_yyyyy_xxxzzz_1, tg_yyyyy_xxyyy_1, tg_yyyyy_xxyyyy_0, \
                                         tg_yyyyy_xxyyyy_1, tg_yyyyy_xxyyyz_0, tg_yyyyy_xxyyyz_1, tg_yyyyy_xxyyz_1, \
                                         tg_yyyyy_xxyyzz_0, tg_yyyyy_xxyyzz_1, tg_yyyyy_xxyzz_1, tg_yyyyy_xxyzzz_0, \
                                         tg_yyyyy_xxyzzz_1, tg_yyyyy_xxzzz_1, tg_yyyyy_xxzzzz_0, tg_yyyyy_xxzzzz_1, \
                                         tg_yyyyy_xyyyy_1, tg_yyyyy_xyyyyy_0, tg_yyyyy_xyyyyy_1, tg_yyyyy_xyyyyz_0, \
                                         tg_yyyyy_xyyyyz_1, tg_yyyyy_xyyyz_1, tg_yyyyy_xyyyzz_0, tg_yyyyy_xyyyzz_1, \
                                         tg_yyyyy_xyyzz_1, tg_yyyyy_xyyzzz_0, tg_yyyyy_xyyzzz_1, tg_yyyyy_xyzzz_1, \
                                         tg_yyyyy_xyzzzz_0, tg_yyyyy_xyzzzz_1, tg_yyyyy_xzzzz_1, tg_yyyyy_xzzzzz_0, \
                                         tg_yyyyy_xzzzzz_1, tg_yyyyy_yyyyy_1, tg_yyyyy_yyyyyy_0, tg_yyyyy_yyyyyy_1, \
                                         tg_yyyyy_yyyyyz_0, tg_yyyyy_yyyyyz_1, tg_yyyyy_yyyyz_1, tg_yyyyy_yyyyzz_0, \
                                         tg_yyyyy_yyyyzz_1, tg_yyyyy_yyyzz_1, tg_yyyyy_yyyzzz_0, tg_yyyyy_yyyzzz_1, \
                                         tg_yyyyy_yyzzz_1, tg_yyyyy_yyzzzz_0, tg_yyyyy_yyzzzz_1, tg_yyyyy_yzzzz_1, \
                                         tg_yyyyy_yzzzzz_0, tg_yyyyy_yzzzzz_1, tg_yyyyy_zzzzz_1, tg_yyyyy_zzzzzz_0, \
                                         tg_yyyyy_zzzzzz_1, tg_yyyyz_xxxxx_1, tg_yyyyz_xxxxxx_0, tg_yyyyz_xxxxxx_1, \
                                         tg_yyyyz_xxxxxy_0, tg_yyyyz_xxxxxy_1, tg_yyyyz_xxxxxz_0, tg_yyyyz_xxxxxz_1, \
                                         tg_yyyyz_xxxxy_1, tg_yyyyz_xxxxyy_0, tg_yyyyz_xxxxyy_1, tg_yyyyz_xxxxyz_0, \
                                         tg_yyyyz_xxxxyz_1, tg_yyyyz_xxxxz_1, tg_yyyyz_xxxxzz_0, tg_yyyyz_xxxxzz_1, \
                                         tg_yyyyz_xxxyy_1, tg_yyyyz_xxxyyy_0, tg_yyyyz_xxxyyy_1, tg_yyyyz_xxxyyz_0, \
                                         tg_yyyyz_xxxyyz_1, tg_yyyyz_xxxyz_1, tg_yyyyz_xxxyzz_0, tg_yyyyz_xxxyzz_1, \
                                         tg_yyyyz_xxxzz_1, tg_yyyyz_xxxzzz_0, tg_yyyyz_xxxzzz_1, tg_yyyyz_xxyyy_1, \
                                         tg_yyyyz_xxyyyy_0, tg_yyyyz_xxyyyy_1, tg_yyyyz_xxyyyz_0, tg_yyyyz_xxyyyz_1, \
                                         tg_yyyyz_xxyyz_1, tg_yyyyz_xxyyzz_0, tg_yyyyz_xxyyzz_1, tg_yyyyz_xxyzz_1, \
                                         tg_yyyyz_xxyzzz_0, tg_yyyyz_xxyzzz_1, tg_yyyyz_xxzzz_1, tg_yyyyz_xxzzzz_0, \
                                         tg_yyyyz_xxzzzz_1, tg_yyyyz_xyyyy_1, tg_yyyyz_xyyyyy_0, tg_yyyyz_xyyyyy_1, \
                                         tg_yyyyz_xyyyyz_0, tg_yyyyz_xyyyyz_1, tg_yyyyz_xyyyz_1, tg_yyyyz_xyyyzz_0, \
                                         tg_yyyyz_xyyyzz_1, tg_yyyyz_xyyzz_1, tg_yyyyz_xyyzzz_0, tg_yyyyz_xyyzzz_1, \
                                         tg_yyyyz_xyzzz_1, tg_yyyyz_xyzzzz_0, tg_yyyyz_xyzzzz_1, tg_yyyyz_xzzzz_1, \
                                         tg_yyyyz_xzzzzz_0, tg_yyyyz_xzzzzz_1, tg_yyyyz_yyyyy_1, tg_yyyyz_yyyyyy_0, \
                                         tg_yyyyz_yyyyyy_1, tg_yyyyz_yyyyyz_0, tg_yyyyz_yyyyyz_1, tg_yyyyz_yyyyz_1, \
                                         tg_yyyyz_yyyyzz_0, tg_yyyyz_yyyyzz_1, tg_yyyyz_yyyzz_1, tg_yyyyz_yyyzzz_0, \
                                         tg_yyyyz_yyyzzz_1, tg_yyyyz_yyzzz_1, tg_yyyyz_yyzzzz_0, tg_yyyyz_yyzzzz_1, \
                                         tg_yyyyz_yzzzz_1, tg_yyyyz_yzzzzz_0, tg_yyyyz_yzzzzz_1, tg_yyyyz_zzzzz_1, \
                                         tg_yyyyz_zzzzzz_0, tg_yyyyz_zzzzzz_1, tg_yyyzz_xxxxx_1, tg_yyyzz_xxxxxx_0, \
                                         tg_yyyzz_xxxxxx_1, tg_yyyzz_xxxxxy_0, tg_yyyzz_xxxxxy_1, tg_yyyzz_xxxxxz_0, \
                                         tg_yyyzz_xxxxxz_1, tg_yyyzz_xxxxy_1, tg_yyyzz_xxxxyy_0, tg_yyyzz_xxxxyy_1, \
                                         tg_yyyzz_xxxxyz_0, tg_yyyzz_xxxxyz_1, tg_yyyzz_xxxxz_1, tg_yyyzz_xxxxzz_0, \
                                         tg_yyyzz_xxxxzz_1, tg_yyyzz_xxxyy_1, tg_yyyzz_xxxyyy_0, tg_yyyzz_xxxyyy_1, \
                                         tg_yyyzz_xxxyyz_0, tg_yyyzz_xxxyyz_1, tg_yyyzz_xxxyz_1, tg_yyyzz_xxxyzz_0, \
                                         tg_yyyzz_xxxyzz_1, tg_yyyzz_xxxzz_1, tg_yyyzz_xxxzzz_0, tg_yyyzz_xxxzzz_1, \
                                         tg_yyyzz_xxyyy_1, tg_yyyzz_xxyyyy_0, tg_yyyzz_xxyyyy_1, tg_yyyzz_xxyyyz_0, \
                                         tg_yyyzz_xxyyyz_1, tg_yyyzz_xxyyz_1, tg_yyyzz_xxyyzz_0, tg_yyyzz_xxyyzz_1, \
                                         tg_yyyzz_xxyzz_1, tg_yyyzz_xxyzzz_0, tg_yyyzz_xxyzzz_1, tg_yyyzz_xxzzz_1, \
                                         tg_yyyzz_xyyyy_1, tg_yyyzz_xyyyz_1, tg_yyyzz_xyyzz_1, tg_yyyzz_xyzzz_1, \
                                         tg_zzzz_xxxxxx_0, tg_zzzz_xxxxxx_1, tg_zzzz_xxxxxy_0, tg_zzzz_xxxxxy_1, \
                                         tg_zzzz_xxxxxz_0, tg_zzzz_xxxxxz_1, tg_zzzz_xxxxyy_0, tg_zzzz_xxxxyy_1, \
                                         tg_zzzz_xxxxyz_0, tg_zzzz_xxxxyz_1, tg_zzzz_xxxxzz_0, tg_zzzz_xxxxzz_1, \
                                         tg_zzzz_xxxyyy_0, tg_zzzz_xxxyyy_1, tg_zzzz_xxxyyz_0, tg_zzzz_xxxyyz_1, \
                                         tg_zzzz_xxxyzz_0, tg_zzzz_xxxyzz_1, tg_zzzz_xxxzzz_0, tg_zzzz_xxxzzz_1, \
                                         tg_zzzz_xxyyyy_0, tg_zzzz_xxyyyy_1, tg_zzzz_xxyyyz_0, tg_zzzz_xxyyyz_1, \
                                         tg_zzzz_xxyyzz_0, tg_zzzz_xxyyzz_1, tg_zzzz_xxyzzz_0, tg_zzzz_xxyzzz_1, \
                                         tg_zzzz_xxzzzz_0, tg_zzzz_xxzzzz_1, tg_zzzz_xyyyyy_0, tg_zzzz_xyyyyy_1, \
                                         tg_zzzz_xyyyyz_0, tg_zzzz_xyyyyz_1, tg_zzzz_xyyyzz_0, tg_zzzz_xyyyzz_1, \
                                         tg_zzzz_xyyzzz_0, tg_zzzz_xyyzzz_1, tg_zzzz_xyzzzz_0, tg_zzzz_xyzzzz_1, \
                                         tg_zzzz_xzzzzz_0, tg_zzzz_xzzzzz_1, tg_zzzz_yyyyyy_0, tg_zzzz_yyyyyy_1, \
                                         tg_zzzz_yyyyyz_0, tg_zzzz_yyyyyz_1, tg_zzzz_yyyyzz_0, tg_zzzz_yyyyzz_1, \
                                         tg_zzzz_yyyzzz_0, tg_zzzz_yyyzzz_1, tg_zzzz_yyzzzz_0, tg_zzzz_yyzzzz_1, \
                                         tg_zzzz_yzzzzz_0, tg_zzzz_yzzzzz_1, tg_zzzz_zzzzzz_0, tg_zzzz_zzzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxzzzz_xxxxxx_0[j] = pb_x * tg_xzzzz_xxxxxx_0[j] + fr * tg_xzzzz_xxxxxx_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxxxx_0[j] - tg_zzzz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xzzzz_xxxxx_1[j];

                    tg_xxzzzz_xxxxxy_0[j] = pb_x * tg_xzzzz_xxxxxy_0[j] + fr * tg_xzzzz_xxxxxy_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxxxy_0[j] - tg_zzzz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xzzzz_xxxxy_1[j];

                    tg_xxzzzz_xxxxxz_0[j] = pb_x * tg_xzzzz_xxxxxz_0[j] + fr * tg_xzzzz_xxxxxz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxxxz_0[j] - tg_zzzz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xzzzz_xxxxz_1[j];

                    tg_xxzzzz_xxxxyy_0[j] = pb_x * tg_xzzzz_xxxxyy_0[j] + fr * tg_xzzzz_xxxxyy_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxxyy_0[j] - tg_zzzz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzzzz_xxxyy_1[j];

                    tg_xxzzzz_xxxxyz_0[j] = pb_x * tg_xzzzz_xxxxyz_0[j] + fr * tg_xzzzz_xxxxyz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxxyz_0[j] - tg_zzzz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzzzz_xxxyz_1[j];

                    tg_xxzzzz_xxxxzz_0[j] = pb_x * tg_xzzzz_xxxxzz_0[j] + fr * tg_xzzzz_xxxxzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxxzz_0[j] - tg_zzzz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzzzz_xxxzz_1[j];

                    tg_xxzzzz_xxxyyy_0[j] = pb_x * tg_xzzzz_xxxyyy_0[j] + fr * tg_xzzzz_xxxyyy_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxyyy_0[j] - tg_zzzz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzzz_xxyyy_1[j];

                    tg_xxzzzz_xxxyyz_0[j] = pb_x * tg_xzzzz_xxxyyz_0[j] + fr * tg_xzzzz_xxxyyz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxyyz_0[j] - tg_zzzz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzzz_xxyyz_1[j];

                    tg_xxzzzz_xxxyzz_0[j] = pb_x * tg_xzzzz_xxxyzz_0[j] + fr * tg_xzzzz_xxxyzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxyzz_0[j] - tg_zzzz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzzz_xxyzz_1[j];

                    tg_xxzzzz_xxxzzz_0[j] = pb_x * tg_xzzzz_xxxzzz_0[j] + fr * tg_xzzzz_xxxzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxzzz_0[j] - tg_zzzz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzzz_xxzzz_1[j];

                    tg_xxzzzz_xxyyyy_0[j] = pb_x * tg_xzzzz_xxyyyy_0[j] + fr * tg_xzzzz_xxyyyy_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxyyyy_0[j] - tg_zzzz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xzzzz_xyyyy_1[j];

                    tg_xxzzzz_xxyyyz_0[j] = pb_x * tg_xzzzz_xxyyyz_0[j] + fr * tg_xzzzz_xxyyyz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxyyyz_0[j] - tg_zzzz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xzzzz_xyyyz_1[j];

                    tg_xxzzzz_xxyyzz_0[j] = pb_x * tg_xzzzz_xxyyzz_0[j] + fr * tg_xzzzz_xxyyzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxyyzz_0[j] - tg_zzzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xzzzz_xyyzz_1[j];

                    tg_xxzzzz_xxyzzz_0[j] = pb_x * tg_xzzzz_xxyzzz_0[j] + fr * tg_xzzzz_xxyzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxyzzz_0[j] - tg_zzzz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xzzzz_xyzzz_1[j];

                    tg_xxzzzz_xxzzzz_0[j] = pb_x * tg_xzzzz_xxzzzz_0[j] + fr * tg_xzzzz_xxzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxzzzz_0[j] - tg_zzzz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xzzzz_xzzzz_1[j];

                    tg_xxzzzz_xyyyyy_0[j] = pb_x * tg_xzzzz_xyyyyy_0[j] + fr * tg_xzzzz_xyyyyy_1[j] + 0.5 * fl1_fx * (tg_zzzz_xyyyyy_0[j] - tg_zzzz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzz_yyyyy_1[j];

                    tg_xxzzzz_xyyyyz_0[j] = pb_x * tg_xzzzz_xyyyyz_0[j] + fr * tg_xzzzz_xyyyyz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xyyyyz_0[j] - tg_zzzz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzz_yyyyz_1[j];

                    tg_xxzzzz_xyyyzz_0[j] = pb_x * tg_xzzzz_xyyyzz_0[j] + fr * tg_xzzzz_xyyyzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xyyyzz_0[j] - tg_zzzz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzz_yyyzz_1[j];

                    tg_xxzzzz_xyyzzz_0[j] = pb_x * tg_xzzzz_xyyzzz_0[j] + fr * tg_xzzzz_xyyzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xyyzzz_0[j] - tg_zzzz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzz_yyzzz_1[j];

                    tg_xxzzzz_xyzzzz_0[j] = pb_x * tg_xzzzz_xyzzzz_0[j] + fr * tg_xzzzz_xyzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xyzzzz_0[j] - tg_zzzz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzz_yzzzz_1[j];

                    tg_xxzzzz_xzzzzz_0[j] = pb_x * tg_xzzzz_xzzzzz_0[j] + fr * tg_xzzzz_xzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xzzzzz_0[j] - tg_zzzz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzz_zzzzz_1[j];

                    tg_xxzzzz_yyyyyy_0[j] = pb_x * tg_xzzzz_yyyyyy_0[j] + fr * tg_xzzzz_yyyyyy_1[j] + 0.5 * fl1_fx * (tg_zzzz_yyyyyy_0[j] - tg_zzzz_yyyyyy_1[j] * fl1_fza);

                    tg_xxzzzz_yyyyyz_0[j] = pb_x * tg_xzzzz_yyyyyz_0[j] + fr * tg_xzzzz_yyyyyz_1[j] + 0.5 * fl1_fx * (tg_zzzz_yyyyyz_0[j] - tg_zzzz_yyyyyz_1[j] * fl1_fza);

                    tg_xxzzzz_yyyyzz_0[j] = pb_x * tg_xzzzz_yyyyzz_0[j] + fr * tg_xzzzz_yyyyzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_yyyyzz_0[j] - tg_zzzz_yyyyzz_1[j] * fl1_fza);

                    tg_xxzzzz_yyyzzz_0[j] = pb_x * tg_xzzzz_yyyzzz_0[j] + fr * tg_xzzzz_yyyzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_yyyzzz_0[j] - tg_zzzz_yyyzzz_1[j] * fl1_fza);

                    tg_xxzzzz_yyzzzz_0[j] = pb_x * tg_xzzzz_yyzzzz_0[j] + fr * tg_xzzzz_yyzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_yyzzzz_0[j] - tg_zzzz_yyzzzz_1[j] * fl1_fza);

                    tg_xxzzzz_yzzzzz_0[j] = pb_x * tg_xzzzz_yzzzzz_0[j] + fr * tg_xzzzz_yzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_yzzzzz_0[j] - tg_zzzz_yzzzzz_1[j] * fl1_fza);

                    tg_xxzzzz_zzzzzz_0[j] = pb_x * tg_xzzzz_zzzzzz_0[j] + fr * tg_xzzzz_zzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_zzzzzz_0[j] - tg_zzzz_zzzzzz_1[j] * fl1_fza);

                    tg_xyyyyy_xxxxxx_0[j] = pb_x * tg_yyyyy_xxxxxx_0[j] + fr * tg_yyyyy_xxxxxx_1[j] + 3.0 * fl1_fxn * tg_yyyyy_xxxxx_1[j];

                    tg_xyyyyy_xxxxxy_0[j] = pb_x * tg_yyyyy_xxxxxy_0[j] + fr * tg_yyyyy_xxxxxy_1[j] + 2.5 * fl1_fxn * tg_yyyyy_xxxxy_1[j];

                    tg_xyyyyy_xxxxxz_0[j] = pb_x * tg_yyyyy_xxxxxz_0[j] + fr * tg_yyyyy_xxxxxz_1[j] + 2.5 * fl1_fxn * tg_yyyyy_xxxxz_1[j];

                    tg_xyyyyy_xxxxyy_0[j] = pb_x * tg_yyyyy_xxxxyy_0[j] + fr * tg_yyyyy_xxxxyy_1[j] + 2.0 * fl1_fxn * tg_yyyyy_xxxyy_1[j];

                    tg_xyyyyy_xxxxyz_0[j] = pb_x * tg_yyyyy_xxxxyz_0[j] + fr * tg_yyyyy_xxxxyz_1[j] + 2.0 * fl1_fxn * tg_yyyyy_xxxyz_1[j];

                    tg_xyyyyy_xxxxzz_0[j] = pb_x * tg_yyyyy_xxxxzz_0[j] + fr * tg_yyyyy_xxxxzz_1[j] + 2.0 * fl1_fxn * tg_yyyyy_xxxzz_1[j];

                    tg_xyyyyy_xxxyyy_0[j] = pb_x * tg_yyyyy_xxxyyy_0[j] + fr * tg_yyyyy_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_yyyyy_xxyyy_1[j];

                    tg_xyyyyy_xxxyyz_0[j] = pb_x * tg_yyyyy_xxxyyz_0[j] + fr * tg_yyyyy_xxxyyz_1[j] + 1.5 * fl1_fxn * tg_yyyyy_xxyyz_1[j];

                    tg_xyyyyy_xxxyzz_0[j] = pb_x * tg_yyyyy_xxxyzz_0[j] + fr * tg_yyyyy_xxxyzz_1[j] + 1.5 * fl1_fxn * tg_yyyyy_xxyzz_1[j];

                    tg_xyyyyy_xxxzzz_0[j] = pb_x * tg_yyyyy_xxxzzz_0[j] + fr * tg_yyyyy_xxxzzz_1[j] + 1.5 * fl1_fxn * tg_yyyyy_xxzzz_1[j];

                    tg_xyyyyy_xxyyyy_0[j] = pb_x * tg_yyyyy_xxyyyy_0[j] + fr * tg_yyyyy_xxyyyy_1[j] + fl1_fxn * tg_yyyyy_xyyyy_1[j];

                    tg_xyyyyy_xxyyyz_0[j] = pb_x * tg_yyyyy_xxyyyz_0[j] + fr * tg_yyyyy_xxyyyz_1[j] + fl1_fxn * tg_yyyyy_xyyyz_1[j];

                    tg_xyyyyy_xxyyzz_0[j] = pb_x * tg_yyyyy_xxyyzz_0[j] + fr * tg_yyyyy_xxyyzz_1[j] + fl1_fxn * tg_yyyyy_xyyzz_1[j];

                    tg_xyyyyy_xxyzzz_0[j] = pb_x * tg_yyyyy_xxyzzz_0[j] + fr * tg_yyyyy_xxyzzz_1[j] + fl1_fxn * tg_yyyyy_xyzzz_1[j];

                    tg_xyyyyy_xxzzzz_0[j] = pb_x * tg_yyyyy_xxzzzz_0[j] + fr * tg_yyyyy_xxzzzz_1[j] + fl1_fxn * tg_yyyyy_xzzzz_1[j];

                    tg_xyyyyy_xyyyyy_0[j] = pb_x * tg_yyyyy_xyyyyy_0[j] + fr * tg_yyyyy_xyyyyy_1[j] + 0.5 * fl1_fxn * tg_yyyyy_yyyyy_1[j];

                    tg_xyyyyy_xyyyyz_0[j] = pb_x * tg_yyyyy_xyyyyz_0[j] + fr * tg_yyyyy_xyyyyz_1[j] + 0.5 * fl1_fxn * tg_yyyyy_yyyyz_1[j];

                    tg_xyyyyy_xyyyzz_0[j] = pb_x * tg_yyyyy_xyyyzz_0[j] + fr * tg_yyyyy_xyyyzz_1[j] + 0.5 * fl1_fxn * tg_yyyyy_yyyzz_1[j];

                    tg_xyyyyy_xyyzzz_0[j] = pb_x * tg_yyyyy_xyyzzz_0[j] + fr * tg_yyyyy_xyyzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyy_yyzzz_1[j];

                    tg_xyyyyy_xyzzzz_0[j] = pb_x * tg_yyyyy_xyzzzz_0[j] + fr * tg_yyyyy_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyy_yzzzz_1[j];

                    tg_xyyyyy_xzzzzz_0[j] = pb_x * tg_yyyyy_xzzzzz_0[j] + fr * tg_yyyyy_xzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyy_zzzzz_1[j];

                    tg_xyyyyy_yyyyyy_0[j] = pb_x * tg_yyyyy_yyyyyy_0[j] + fr * tg_yyyyy_yyyyyy_1[j];

                    tg_xyyyyy_yyyyyz_0[j] = pb_x * tg_yyyyy_yyyyyz_0[j] + fr * tg_yyyyy_yyyyyz_1[j];

                    tg_xyyyyy_yyyyzz_0[j] = pb_x * tg_yyyyy_yyyyzz_0[j] + fr * tg_yyyyy_yyyyzz_1[j];

                    tg_xyyyyy_yyyzzz_0[j] = pb_x * tg_yyyyy_yyyzzz_0[j] + fr * tg_yyyyy_yyyzzz_1[j];

                    tg_xyyyyy_yyzzzz_0[j] = pb_x * tg_yyyyy_yyzzzz_0[j] + fr * tg_yyyyy_yyzzzz_1[j];

                    tg_xyyyyy_yzzzzz_0[j] = pb_x * tg_yyyyy_yzzzzz_0[j] + fr * tg_yyyyy_yzzzzz_1[j];

                    tg_xyyyyy_zzzzzz_0[j] = pb_x * tg_yyyyy_zzzzzz_0[j] + fr * tg_yyyyy_zzzzzz_1[j];

                    tg_xyyyyz_xxxxxx_0[j] = pb_x * tg_yyyyz_xxxxxx_0[j] + fr * tg_yyyyz_xxxxxx_1[j] + 3.0 * fl1_fxn * tg_yyyyz_xxxxx_1[j];

                    tg_xyyyyz_xxxxxy_0[j] = pb_x * tg_yyyyz_xxxxxy_0[j] + fr * tg_yyyyz_xxxxxy_1[j] + 2.5 * fl1_fxn * tg_yyyyz_xxxxy_1[j];

                    tg_xyyyyz_xxxxxz_0[j] = pb_x * tg_yyyyz_xxxxxz_0[j] + fr * tg_yyyyz_xxxxxz_1[j] + 2.5 * fl1_fxn * tg_yyyyz_xxxxz_1[j];

                    tg_xyyyyz_xxxxyy_0[j] = pb_x * tg_yyyyz_xxxxyy_0[j] + fr * tg_yyyyz_xxxxyy_1[j] + 2.0 * fl1_fxn * tg_yyyyz_xxxyy_1[j];

                    tg_xyyyyz_xxxxyz_0[j] = pb_x * tg_yyyyz_xxxxyz_0[j] + fr * tg_yyyyz_xxxxyz_1[j] + 2.0 * fl1_fxn * tg_yyyyz_xxxyz_1[j];

                    tg_xyyyyz_xxxxzz_0[j] = pb_x * tg_yyyyz_xxxxzz_0[j] + fr * tg_yyyyz_xxxxzz_1[j] + 2.0 * fl1_fxn * tg_yyyyz_xxxzz_1[j];

                    tg_xyyyyz_xxxyyy_0[j] = pb_x * tg_yyyyz_xxxyyy_0[j] + fr * tg_yyyyz_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_yyyyz_xxyyy_1[j];

                    tg_xyyyyz_xxxyyz_0[j] = pb_x * tg_yyyyz_xxxyyz_0[j] + fr * tg_yyyyz_xxxyyz_1[j] + 1.5 * fl1_fxn * tg_yyyyz_xxyyz_1[j];

                    tg_xyyyyz_xxxyzz_0[j] = pb_x * tg_yyyyz_xxxyzz_0[j] + fr * tg_yyyyz_xxxyzz_1[j] + 1.5 * fl1_fxn * tg_yyyyz_xxyzz_1[j];

                    tg_xyyyyz_xxxzzz_0[j] = pb_x * tg_yyyyz_xxxzzz_0[j] + fr * tg_yyyyz_xxxzzz_1[j] + 1.5 * fl1_fxn * tg_yyyyz_xxzzz_1[j];

                    tg_xyyyyz_xxyyyy_0[j] = pb_x * tg_yyyyz_xxyyyy_0[j] + fr * tg_yyyyz_xxyyyy_1[j] + fl1_fxn * tg_yyyyz_xyyyy_1[j];

                    tg_xyyyyz_xxyyyz_0[j] = pb_x * tg_yyyyz_xxyyyz_0[j] + fr * tg_yyyyz_xxyyyz_1[j] + fl1_fxn * tg_yyyyz_xyyyz_1[j];

                    tg_xyyyyz_xxyyzz_0[j] = pb_x * tg_yyyyz_xxyyzz_0[j] + fr * tg_yyyyz_xxyyzz_1[j] + fl1_fxn * tg_yyyyz_xyyzz_1[j];

                    tg_xyyyyz_xxyzzz_0[j] = pb_x * tg_yyyyz_xxyzzz_0[j] + fr * tg_yyyyz_xxyzzz_1[j] + fl1_fxn * tg_yyyyz_xyzzz_1[j];

                    tg_xyyyyz_xxzzzz_0[j] = pb_x * tg_yyyyz_xxzzzz_0[j] + fr * tg_yyyyz_xxzzzz_1[j] + fl1_fxn * tg_yyyyz_xzzzz_1[j];

                    tg_xyyyyz_xyyyyy_0[j] = pb_x * tg_yyyyz_xyyyyy_0[j] + fr * tg_yyyyz_xyyyyy_1[j] + 0.5 * fl1_fxn * tg_yyyyz_yyyyy_1[j];

                    tg_xyyyyz_xyyyyz_0[j] = pb_x * tg_yyyyz_xyyyyz_0[j] + fr * tg_yyyyz_xyyyyz_1[j] + 0.5 * fl1_fxn * tg_yyyyz_yyyyz_1[j];

                    tg_xyyyyz_xyyyzz_0[j] = pb_x * tg_yyyyz_xyyyzz_0[j] + fr * tg_yyyyz_xyyyzz_1[j] + 0.5 * fl1_fxn * tg_yyyyz_yyyzz_1[j];

                    tg_xyyyyz_xyyzzz_0[j] = pb_x * tg_yyyyz_xyyzzz_0[j] + fr * tg_yyyyz_xyyzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyz_yyzzz_1[j];

                    tg_xyyyyz_xyzzzz_0[j] = pb_x * tg_yyyyz_xyzzzz_0[j] + fr * tg_yyyyz_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyz_yzzzz_1[j];

                    tg_xyyyyz_xzzzzz_0[j] = pb_x * tg_yyyyz_xzzzzz_0[j] + fr * tg_yyyyz_xzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyz_zzzzz_1[j];

                    tg_xyyyyz_yyyyyy_0[j] = pb_x * tg_yyyyz_yyyyyy_0[j] + fr * tg_yyyyz_yyyyyy_1[j];

                    tg_xyyyyz_yyyyyz_0[j] = pb_x * tg_yyyyz_yyyyyz_0[j] + fr * tg_yyyyz_yyyyyz_1[j];

                    tg_xyyyyz_yyyyzz_0[j] = pb_x * tg_yyyyz_yyyyzz_0[j] + fr * tg_yyyyz_yyyyzz_1[j];

                    tg_xyyyyz_yyyzzz_0[j] = pb_x * tg_yyyyz_yyyzzz_0[j] + fr * tg_yyyyz_yyyzzz_1[j];

                    tg_xyyyyz_yyzzzz_0[j] = pb_x * tg_yyyyz_yyzzzz_0[j] + fr * tg_yyyyz_yyzzzz_1[j];

                    tg_xyyyyz_yzzzzz_0[j] = pb_x * tg_yyyyz_yzzzzz_0[j] + fr * tg_yyyyz_yzzzzz_1[j];

                    tg_xyyyyz_zzzzzz_0[j] = pb_x * tg_yyyyz_zzzzzz_0[j] + fr * tg_yyyyz_zzzzzz_1[j];

                    tg_xyyyzz_xxxxxx_0[j] = pb_x * tg_yyyzz_xxxxxx_0[j] + fr * tg_yyyzz_xxxxxx_1[j] + 3.0 * fl1_fxn * tg_yyyzz_xxxxx_1[j];

                    tg_xyyyzz_xxxxxy_0[j] = pb_x * tg_yyyzz_xxxxxy_0[j] + fr * tg_yyyzz_xxxxxy_1[j] + 2.5 * fl1_fxn * tg_yyyzz_xxxxy_1[j];

                    tg_xyyyzz_xxxxxz_0[j] = pb_x * tg_yyyzz_xxxxxz_0[j] + fr * tg_yyyzz_xxxxxz_1[j] + 2.5 * fl1_fxn * tg_yyyzz_xxxxz_1[j];

                    tg_xyyyzz_xxxxyy_0[j] = pb_x * tg_yyyzz_xxxxyy_0[j] + fr * tg_yyyzz_xxxxyy_1[j] + 2.0 * fl1_fxn * tg_yyyzz_xxxyy_1[j];

                    tg_xyyyzz_xxxxyz_0[j] = pb_x * tg_yyyzz_xxxxyz_0[j] + fr * tg_yyyzz_xxxxyz_1[j] + 2.0 * fl1_fxn * tg_yyyzz_xxxyz_1[j];

                    tg_xyyyzz_xxxxzz_0[j] = pb_x * tg_yyyzz_xxxxzz_0[j] + fr * tg_yyyzz_xxxxzz_1[j] + 2.0 * fl1_fxn * tg_yyyzz_xxxzz_1[j];

                    tg_xyyyzz_xxxyyy_0[j] = pb_x * tg_yyyzz_xxxyyy_0[j] + fr * tg_yyyzz_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_yyyzz_xxyyy_1[j];

                    tg_xyyyzz_xxxyyz_0[j] = pb_x * tg_yyyzz_xxxyyz_0[j] + fr * tg_yyyzz_xxxyyz_1[j] + 1.5 * fl1_fxn * tg_yyyzz_xxyyz_1[j];

                    tg_xyyyzz_xxxyzz_0[j] = pb_x * tg_yyyzz_xxxyzz_0[j] + fr * tg_yyyzz_xxxyzz_1[j] + 1.5 * fl1_fxn * tg_yyyzz_xxyzz_1[j];

                    tg_xyyyzz_xxxzzz_0[j] = pb_x * tg_yyyzz_xxxzzz_0[j] + fr * tg_yyyzz_xxxzzz_1[j] + 1.5 * fl1_fxn * tg_yyyzz_xxzzz_1[j];

                    tg_xyyyzz_xxyyyy_0[j] = pb_x * tg_yyyzz_xxyyyy_0[j] + fr * tg_yyyzz_xxyyyy_1[j] + fl1_fxn * tg_yyyzz_xyyyy_1[j];

                    tg_xyyyzz_xxyyyz_0[j] = pb_x * tg_yyyzz_xxyyyz_0[j] + fr * tg_yyyzz_xxyyyz_1[j] + fl1_fxn * tg_yyyzz_xyyyz_1[j];

                    tg_xyyyzz_xxyyzz_0[j] = pb_x * tg_yyyzz_xxyyzz_0[j] + fr * tg_yyyzz_xxyyzz_1[j] + fl1_fxn * tg_yyyzz_xyyzz_1[j];

                    tg_xyyyzz_xxyzzz_0[j] = pb_x * tg_yyyzz_xxyzzz_0[j] + fr * tg_yyyzz_xxyzzz_1[j] + fl1_fxn * tg_yyyzz_xyzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSISI_490_588(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (490,588)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {6, -1, -1, -1},
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_6_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_6_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_5_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
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

                auto tg_yyyzz_xxzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 490); 

                auto tg_yyyzz_xyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 491); 

                auto tg_yyyzz_xyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 492); 

                auto tg_yyyzz_xyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 493); 

                auto tg_yyyzz_xyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 494); 

                auto tg_yyyzz_xyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 495); 

                auto tg_yyyzz_xzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 496); 

                auto tg_yyyzz_yyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 497); 

                auto tg_yyyzz_yyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 498); 

                auto tg_yyyzz_yyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 499); 

                auto tg_yyyzz_yyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 500); 

                auto tg_yyyzz_yyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 501); 

                auto tg_yyyzz_yzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 502); 

                auto tg_yyyzz_zzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 503); 

                auto tg_yyzzz_xxxxxx_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 504); 

                auto tg_yyzzz_xxxxxy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 505); 

                auto tg_yyzzz_xxxxxz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 506); 

                auto tg_yyzzz_xxxxyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 507); 

                auto tg_yyzzz_xxxxyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 508); 

                auto tg_yyzzz_xxxxzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 509); 

                auto tg_yyzzz_xxxyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 510); 

                auto tg_yyzzz_xxxyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 511); 

                auto tg_yyzzz_xxxyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 512); 

                auto tg_yyzzz_xxxzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 513); 

                auto tg_yyzzz_xxyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 514); 

                auto tg_yyzzz_xxyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 515); 

                auto tg_yyzzz_xxyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 516); 

                auto tg_yyzzz_xxyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 517); 

                auto tg_yyzzz_xxzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 518); 

                auto tg_yyzzz_xyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 519); 

                auto tg_yyzzz_xyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 520); 

                auto tg_yyzzz_xyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 521); 

                auto tg_yyzzz_xyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 522); 

                auto tg_yyzzz_xyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 523); 

                auto tg_yyzzz_xzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 524); 

                auto tg_yyzzz_yyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 525); 

                auto tg_yyzzz_yyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 526); 

                auto tg_yyzzz_yyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 527); 

                auto tg_yyzzz_yyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 528); 

                auto tg_yyzzz_yyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 529); 

                auto tg_yyzzz_yzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 530); 

                auto tg_yyzzz_zzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 531); 

                auto tg_yzzzz_xxxxxx_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 532); 

                auto tg_yzzzz_xxxxxy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 533); 

                auto tg_yzzzz_xxxxxz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 534); 

                auto tg_yzzzz_xxxxyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 535); 

                auto tg_yzzzz_xxxxyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 536); 

                auto tg_yzzzz_xxxxzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 537); 

                auto tg_yzzzz_xxxyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 538); 

                auto tg_yzzzz_xxxyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 539); 

                auto tg_yzzzz_xxxyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 540); 

                auto tg_yzzzz_xxxzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 541); 

                auto tg_yzzzz_xxyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 542); 

                auto tg_yzzzz_xxyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 543); 

                auto tg_yzzzz_xxyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 544); 

                auto tg_yzzzz_xxyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 545); 

                auto tg_yzzzz_xxzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 546); 

                auto tg_yzzzz_xyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 547); 

                auto tg_yzzzz_xyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 548); 

                auto tg_yzzzz_xyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 549); 

                auto tg_yzzzz_xyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 550); 

                auto tg_yzzzz_xyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 551); 

                auto tg_yzzzz_xzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 552); 

                auto tg_yzzzz_yyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 553); 

                auto tg_yzzzz_yyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 554); 

                auto tg_yzzzz_yyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 555); 

                auto tg_yzzzz_yyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 556); 

                auto tg_yzzzz_yyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 557); 

                auto tg_yzzzz_yzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 558); 

                auto tg_yzzzz_zzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 559); 

                auto tg_zzzzz_xxxxxx_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 560); 

                auto tg_zzzzz_xxxxxy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 561); 

                auto tg_zzzzz_xxxxxz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 562); 

                auto tg_zzzzz_xxxxyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 563); 

                auto tg_zzzzz_xxxxyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 564); 

                auto tg_zzzzz_xxxxzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 565); 

                auto tg_zzzzz_xxxyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 566); 

                auto tg_zzzzz_xxxyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 567); 

                auto tg_zzzzz_xxxyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 568); 

                auto tg_zzzzz_xxxzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 569); 

                auto tg_zzzzz_xxyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 570); 

                auto tg_zzzzz_xxyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 571); 

                auto tg_zzzzz_xxyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 572); 

                auto tg_zzzzz_xxyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 573); 

                auto tg_zzzzz_xxzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 574); 

                auto tg_zzzzz_xyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 575); 

                auto tg_zzzzz_xyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 576); 

                auto tg_zzzzz_xyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 577); 

                auto tg_zzzzz_xyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 578); 

                auto tg_zzzzz_xyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 579); 

                auto tg_zzzzz_xzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 580); 

                auto tg_zzzzz_yyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 581); 

                auto tg_zzzzz_yyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 582); 

                auto tg_zzzzz_yyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 583); 

                auto tg_zzzzz_yyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 584); 

                auto tg_zzzzz_yyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 585); 

                auto tg_zzzzz_yzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 586); 

                auto tg_zzzzz_zzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 587); 

                auto tg_yyyzz_xxzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 490); 

                auto tg_yyyzz_xyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 491); 

                auto tg_yyyzz_xyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 492); 

                auto tg_yyyzz_xyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 493); 

                auto tg_yyyzz_xyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 494); 

                auto tg_yyyzz_xyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 495); 

                auto tg_yyyzz_xzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 496); 

                auto tg_yyyzz_yyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 497); 

                auto tg_yyyzz_yyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 498); 

                auto tg_yyyzz_yyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 499); 

                auto tg_yyyzz_yyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 500); 

                auto tg_yyyzz_yyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 501); 

                auto tg_yyyzz_yzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 502); 

                auto tg_yyyzz_zzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 503); 

                auto tg_yyzzz_xxxxxx_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 504); 

                auto tg_yyzzz_xxxxxy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 505); 

                auto tg_yyzzz_xxxxxz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 506); 

                auto tg_yyzzz_xxxxyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 507); 

                auto tg_yyzzz_xxxxyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 508); 

                auto tg_yyzzz_xxxxzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 509); 

                auto tg_yyzzz_xxxyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 510); 

                auto tg_yyzzz_xxxyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 511); 

                auto tg_yyzzz_xxxyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 512); 

                auto tg_yyzzz_xxxzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 513); 

                auto tg_yyzzz_xxyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 514); 

                auto tg_yyzzz_xxyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 515); 

                auto tg_yyzzz_xxyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 516); 

                auto tg_yyzzz_xxyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 517); 

                auto tg_yyzzz_xxzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 518); 

                auto tg_yyzzz_xyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 519); 

                auto tg_yyzzz_xyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 520); 

                auto tg_yyzzz_xyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 521); 

                auto tg_yyzzz_xyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 522); 

                auto tg_yyzzz_xyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 523); 

                auto tg_yyzzz_xzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 524); 

                auto tg_yyzzz_yyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 525); 

                auto tg_yyzzz_yyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 526); 

                auto tg_yyzzz_yyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 527); 

                auto tg_yyzzz_yyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 528); 

                auto tg_yyzzz_yyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 529); 

                auto tg_yyzzz_yzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 530); 

                auto tg_yyzzz_zzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 531); 

                auto tg_yzzzz_xxxxxx_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 532); 

                auto tg_yzzzz_xxxxxy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 533); 

                auto tg_yzzzz_xxxxxz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 534); 

                auto tg_yzzzz_xxxxyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 535); 

                auto tg_yzzzz_xxxxyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 536); 

                auto tg_yzzzz_xxxxzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 537); 

                auto tg_yzzzz_xxxyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 538); 

                auto tg_yzzzz_xxxyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 539); 

                auto tg_yzzzz_xxxyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 540); 

                auto tg_yzzzz_xxxzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 541); 

                auto tg_yzzzz_xxyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 542); 

                auto tg_yzzzz_xxyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 543); 

                auto tg_yzzzz_xxyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 544); 

                auto tg_yzzzz_xxyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 545); 

                auto tg_yzzzz_xxzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 546); 

                auto tg_yzzzz_xyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 547); 

                auto tg_yzzzz_xyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 548); 

                auto tg_yzzzz_xyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 549); 

                auto tg_yzzzz_xyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 550); 

                auto tg_yzzzz_xyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 551); 

                auto tg_yzzzz_xzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 552); 

                auto tg_yzzzz_yyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 553); 

                auto tg_yzzzz_yyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 554); 

                auto tg_yzzzz_yyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 555); 

                auto tg_yzzzz_yyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 556); 

                auto tg_yzzzz_yyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 557); 

                auto tg_yzzzz_yzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 558); 

                auto tg_yzzzz_zzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 559); 

                auto tg_zzzzz_xxxxxx_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 560); 

                auto tg_zzzzz_xxxxxy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 561); 

                auto tg_zzzzz_xxxxxz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 562); 

                auto tg_zzzzz_xxxxyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 563); 

                auto tg_zzzzz_xxxxyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 564); 

                auto tg_zzzzz_xxxxzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 565); 

                auto tg_zzzzz_xxxyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 566); 

                auto tg_zzzzz_xxxyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 567); 

                auto tg_zzzzz_xxxyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 568); 

                auto tg_zzzzz_xxxzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 569); 

                auto tg_zzzzz_xxyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 570); 

                auto tg_zzzzz_xxyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 571); 

                auto tg_zzzzz_xxyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 572); 

                auto tg_zzzzz_xxyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 573); 

                auto tg_zzzzz_xxzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 574); 

                auto tg_zzzzz_xyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 575); 

                auto tg_zzzzz_xyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 576); 

                auto tg_zzzzz_xyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 577); 

                auto tg_zzzzz_xyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 578); 

                auto tg_zzzzz_xyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 579); 

                auto tg_zzzzz_xzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 580); 

                auto tg_zzzzz_yyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 581); 

                auto tg_zzzzz_yyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 582); 

                auto tg_zzzzz_yyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 583); 

                auto tg_zzzzz_yyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 584); 

                auto tg_zzzzz_yyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 585); 

                auto tg_zzzzz_yzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 586); 

                auto tg_zzzzz_zzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 587); 

                auto tg_yyyzz_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 371); 

                auto tg_yyyzz_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 372); 

                auto tg_yyyzz_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 373); 

                auto tg_yyyzz_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 374); 

                auto tg_yyyzz_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 375); 

                auto tg_yyyzz_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 376); 

                auto tg_yyyzz_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 377); 

                auto tg_yyzzz_xxxxx_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 378); 

                auto tg_yyzzz_xxxxy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 379); 

                auto tg_yyzzz_xxxxz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 380); 

                auto tg_yyzzz_xxxyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 381); 

                auto tg_yyzzz_xxxyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 382); 

                auto tg_yyzzz_xxxzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 383); 

                auto tg_yyzzz_xxyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 384); 

                auto tg_yyzzz_xxyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 385); 

                auto tg_yyzzz_xxyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 386); 

                auto tg_yyzzz_xxzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 387); 

                auto tg_yyzzz_xyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 388); 

                auto tg_yyzzz_xyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 389); 

                auto tg_yyzzz_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 390); 

                auto tg_yyzzz_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 391); 

                auto tg_yyzzz_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 392); 

                auto tg_yyzzz_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 393); 

                auto tg_yyzzz_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 394); 

                auto tg_yyzzz_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 395); 

                auto tg_yyzzz_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 396); 

                auto tg_yyzzz_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 397); 

                auto tg_yyzzz_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 398); 

                auto tg_yzzzz_xxxxx_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 399); 

                auto tg_yzzzz_xxxxy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 400); 

                auto tg_yzzzz_xxxxz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 401); 

                auto tg_yzzzz_xxxyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 402); 

                auto tg_yzzzz_xxxyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 403); 

                auto tg_yzzzz_xxxzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 404); 

                auto tg_yzzzz_xxyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 405); 

                auto tg_yzzzz_xxyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 406); 

                auto tg_yzzzz_xxyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 407); 

                auto tg_yzzzz_xxzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 408); 

                auto tg_yzzzz_xyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 409); 

                auto tg_yzzzz_xyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 410); 

                auto tg_yzzzz_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 411); 

                auto tg_yzzzz_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 412); 

                auto tg_yzzzz_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 413); 

                auto tg_yzzzz_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 414); 

                auto tg_yzzzz_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 415); 

                auto tg_yzzzz_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 416); 

                auto tg_yzzzz_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 417); 

                auto tg_yzzzz_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 418); 

                auto tg_yzzzz_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 419); 

                auto tg_zzzzz_xxxxx_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 420); 

                auto tg_zzzzz_xxxxy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 421); 

                auto tg_zzzzz_xxxxz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 422); 

                auto tg_zzzzz_xxxyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 423); 

                auto tg_zzzzz_xxxyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 424); 

                auto tg_zzzzz_xxxzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 425); 

                auto tg_zzzzz_xxyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 426); 

                auto tg_zzzzz_xxyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 427); 

                auto tg_zzzzz_xxyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 428); 

                auto tg_zzzzz_xxzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 429); 

                auto tg_zzzzz_xyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 430); 

                auto tg_zzzzz_xyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 431); 

                auto tg_zzzzz_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 432); 

                auto tg_zzzzz_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 433); 

                auto tg_zzzzz_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 434); 

                auto tg_zzzzz_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 435); 

                auto tg_zzzzz_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 436); 

                auto tg_zzzzz_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 437); 

                auto tg_zzzzz_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 438); 

                auto tg_zzzzz_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 439); 

                auto tg_zzzzz_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 440); 

                // set up pointers to integrals

                auto tg_xyyyzz_xxzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 490); 

                auto tg_xyyyzz_xyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 491); 

                auto tg_xyyyzz_xyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 492); 

                auto tg_xyyyzz_xyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 493); 

                auto tg_xyyyzz_xyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 494); 

                auto tg_xyyyzz_xyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 495); 

                auto tg_xyyyzz_xzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 496); 

                auto tg_xyyyzz_yyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 497); 

                auto tg_xyyyzz_yyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 498); 

                auto tg_xyyyzz_yyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 499); 

                auto tg_xyyyzz_yyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 500); 

                auto tg_xyyyzz_yyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 501); 

                auto tg_xyyyzz_yzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 502); 

                auto tg_xyyyzz_zzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 503); 

                auto tg_xyyzzz_xxxxxx_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 504); 

                auto tg_xyyzzz_xxxxxy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 505); 

                auto tg_xyyzzz_xxxxxz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 506); 

                auto tg_xyyzzz_xxxxyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 507); 

                auto tg_xyyzzz_xxxxyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 508); 

                auto tg_xyyzzz_xxxxzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 509); 

                auto tg_xyyzzz_xxxyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 510); 

                auto tg_xyyzzz_xxxyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 511); 

                auto tg_xyyzzz_xxxyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 512); 

                auto tg_xyyzzz_xxxzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 513); 

                auto tg_xyyzzz_xxyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 514); 

                auto tg_xyyzzz_xxyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 515); 

                auto tg_xyyzzz_xxyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 516); 

                auto tg_xyyzzz_xxyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 517); 

                auto tg_xyyzzz_xxzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 518); 

                auto tg_xyyzzz_xyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 519); 

                auto tg_xyyzzz_xyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 520); 

                auto tg_xyyzzz_xyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 521); 

                auto tg_xyyzzz_xyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 522); 

                auto tg_xyyzzz_xyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 523); 

                auto tg_xyyzzz_xzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 524); 

                auto tg_xyyzzz_yyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 525); 

                auto tg_xyyzzz_yyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 526); 

                auto tg_xyyzzz_yyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 527); 

                auto tg_xyyzzz_yyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 528); 

                auto tg_xyyzzz_yyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 529); 

                auto tg_xyyzzz_yzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 530); 

                auto tg_xyyzzz_zzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 531); 

                auto tg_xyzzzz_xxxxxx_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 532); 

                auto tg_xyzzzz_xxxxxy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 533); 

                auto tg_xyzzzz_xxxxxz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 534); 

                auto tg_xyzzzz_xxxxyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 535); 

                auto tg_xyzzzz_xxxxyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 536); 

                auto tg_xyzzzz_xxxxzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 537); 

                auto tg_xyzzzz_xxxyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 538); 

                auto tg_xyzzzz_xxxyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 539); 

                auto tg_xyzzzz_xxxyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 540); 

                auto tg_xyzzzz_xxxzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 541); 

                auto tg_xyzzzz_xxyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 542); 

                auto tg_xyzzzz_xxyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 543); 

                auto tg_xyzzzz_xxyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 544); 

                auto tg_xyzzzz_xxyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 545); 

                auto tg_xyzzzz_xxzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 546); 

                auto tg_xyzzzz_xyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 547); 

                auto tg_xyzzzz_xyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 548); 

                auto tg_xyzzzz_xyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 549); 

                auto tg_xyzzzz_xyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 550); 

                auto tg_xyzzzz_xyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 551); 

                auto tg_xyzzzz_xzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 552); 

                auto tg_xyzzzz_yyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 553); 

                auto tg_xyzzzz_yyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 554); 

                auto tg_xyzzzz_yyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 555); 

                auto tg_xyzzzz_yyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 556); 

                auto tg_xyzzzz_yyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 557); 

                auto tg_xyzzzz_yzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 558); 

                auto tg_xyzzzz_zzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 559); 

                auto tg_xzzzzz_xxxxxx_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 560); 

                auto tg_xzzzzz_xxxxxy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 561); 

                auto tg_xzzzzz_xxxxxz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 562); 

                auto tg_xzzzzz_xxxxyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 563); 

                auto tg_xzzzzz_xxxxyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 564); 

                auto tg_xzzzzz_xxxxzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 565); 

                auto tg_xzzzzz_xxxyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 566); 

                auto tg_xzzzzz_xxxyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 567); 

                auto tg_xzzzzz_xxxyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 568); 

                auto tg_xzzzzz_xxxzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 569); 

                auto tg_xzzzzz_xxyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 570); 

                auto tg_xzzzzz_xxyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 571); 

                auto tg_xzzzzz_xxyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 572); 

                auto tg_xzzzzz_xxyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 573); 

                auto tg_xzzzzz_xxzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 574); 

                auto tg_xzzzzz_xyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 575); 

                auto tg_xzzzzz_xyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 576); 

                auto tg_xzzzzz_xyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 577); 

                auto tg_xzzzzz_xyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 578); 

                auto tg_xzzzzz_xyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 579); 

                auto tg_xzzzzz_xzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 580); 

                auto tg_xzzzzz_yyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 581); 

                auto tg_xzzzzz_yyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 582); 

                auto tg_xzzzzz_yyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 583); 

                auto tg_xzzzzz_yyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 584); 

                auto tg_xzzzzz_yyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 585); 

                auto tg_xzzzzz_yzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 586); 

                auto tg_xzzzzz_zzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 587); 

                // Batch of Integrals (490,588)

                #pragma omp simd aligned(fxn, tg_xyyyzz_xxzzzz_0, tg_xyyyzz_xyyyyy_0, tg_xyyyzz_xyyyyz_0, \
                                         tg_xyyyzz_xyyyzz_0, tg_xyyyzz_xyyzzz_0, tg_xyyyzz_xyzzzz_0, tg_xyyyzz_xzzzzz_0, \
                                         tg_xyyyzz_yyyyyy_0, tg_xyyyzz_yyyyyz_0, tg_xyyyzz_yyyyzz_0, tg_xyyyzz_yyyzzz_0, \
                                         tg_xyyyzz_yyzzzz_0, tg_xyyyzz_yzzzzz_0, tg_xyyyzz_zzzzzz_0, tg_xyyzzz_xxxxxx_0, \
                                         tg_xyyzzz_xxxxxy_0, tg_xyyzzz_xxxxxz_0, tg_xyyzzz_xxxxyy_0, tg_xyyzzz_xxxxyz_0, \
                                         tg_xyyzzz_xxxxzz_0, tg_xyyzzz_xxxyyy_0, tg_xyyzzz_xxxyyz_0, tg_xyyzzz_xxxyzz_0, \
                                         tg_xyyzzz_xxxzzz_0, tg_xyyzzz_xxyyyy_0, tg_xyyzzz_xxyyyz_0, tg_xyyzzz_xxyyzz_0, \
                                         tg_xyyzzz_xxyzzz_0, tg_xyyzzz_xxzzzz_0, tg_xyyzzz_xyyyyy_0, tg_xyyzzz_xyyyyz_0, \
                                         tg_xyyzzz_xyyyzz_0, tg_xyyzzz_xyyzzz_0, tg_xyyzzz_xyzzzz_0, tg_xyyzzz_xzzzzz_0, \
                                         tg_xyyzzz_yyyyyy_0, tg_xyyzzz_yyyyyz_0, tg_xyyzzz_yyyyzz_0, tg_xyyzzz_yyyzzz_0, \
                                         tg_xyyzzz_yyzzzz_0, tg_xyyzzz_yzzzzz_0, tg_xyyzzz_zzzzzz_0, tg_xyzzzz_xxxxxx_0, \
                                         tg_xyzzzz_xxxxxy_0, tg_xyzzzz_xxxxxz_0, tg_xyzzzz_xxxxyy_0, tg_xyzzzz_xxxxyz_0, \
                                         tg_xyzzzz_xxxxzz_0, tg_xyzzzz_xxxyyy_0, tg_xyzzzz_xxxyyz_0, tg_xyzzzz_xxxyzz_0, \
                                         tg_xyzzzz_xxxzzz_0, tg_xyzzzz_xxyyyy_0, tg_xyzzzz_xxyyyz_0, tg_xyzzzz_xxyyzz_0, \
                                         tg_xyzzzz_xxyzzz_0, tg_xyzzzz_xxzzzz_0, tg_xyzzzz_xyyyyy_0, tg_xyzzzz_xyyyyz_0, \
                                         tg_xyzzzz_xyyyzz_0, tg_xyzzzz_xyyzzz_0, tg_xyzzzz_xyzzzz_0, tg_xyzzzz_xzzzzz_0, \
                                         tg_xyzzzz_yyyyyy_0, tg_xyzzzz_yyyyyz_0, tg_xyzzzz_yyyyzz_0, tg_xyzzzz_yyyzzz_0, \
                                         tg_xyzzzz_yyzzzz_0, tg_xyzzzz_yzzzzz_0, tg_xyzzzz_zzzzzz_0, tg_xzzzzz_xxxxxx_0, \
                                         tg_xzzzzz_xxxxxy_0, tg_xzzzzz_xxxxxz_0, tg_xzzzzz_xxxxyy_0, tg_xzzzzz_xxxxyz_0, \
                                         tg_xzzzzz_xxxxzz_0, tg_xzzzzz_xxxyyy_0, tg_xzzzzz_xxxyyz_0, tg_xzzzzz_xxxyzz_0, \
                                         tg_xzzzzz_xxxzzz_0, tg_xzzzzz_xxyyyy_0, tg_xzzzzz_xxyyyz_0, tg_xzzzzz_xxyyzz_0, \
                                         tg_xzzzzz_xxyzzz_0, tg_xzzzzz_xxzzzz_0, tg_xzzzzz_xyyyyy_0, tg_xzzzzz_xyyyyz_0, \
                                         tg_xzzzzz_xyyyzz_0, tg_xzzzzz_xyyzzz_0, tg_xzzzzz_xyzzzz_0, tg_xzzzzz_xzzzzz_0, \
                                         tg_xzzzzz_yyyyyy_0, tg_xzzzzz_yyyyyz_0, tg_xzzzzz_yyyyzz_0, tg_xzzzzz_yyyzzz_0, \
                                         tg_xzzzzz_yyzzzz_0, tg_xzzzzz_yzzzzz_0, tg_xzzzzz_zzzzzz_0, tg_yyyzz_xxzzzz_0, \
                                         tg_yyyzz_xxzzzz_1, tg_yyyzz_xyyyyy_0, tg_yyyzz_xyyyyy_1, tg_yyyzz_xyyyyz_0, \
                                         tg_yyyzz_xyyyyz_1, tg_yyyzz_xyyyzz_0, tg_yyyzz_xyyyzz_1, tg_yyyzz_xyyzzz_0, \
                                         tg_yyyzz_xyyzzz_1, tg_yyyzz_xyzzzz_0, tg_yyyzz_xyzzzz_1, tg_yyyzz_xzzzz_1, \
                                         tg_yyyzz_xzzzzz_0, tg_yyyzz_xzzzzz_1, tg_yyyzz_yyyyy_1, tg_yyyzz_yyyyyy_0, \
                                         tg_yyyzz_yyyyyy_1, tg_yyyzz_yyyyyz_0, tg_yyyzz_yyyyyz_1, tg_yyyzz_yyyyz_1, \
                                         tg_yyyzz_yyyyzz_0, tg_yyyzz_yyyyzz_1, tg_yyyzz_yyyzz_1, tg_yyyzz_yyyzzz_0, \
                                         tg_yyyzz_yyyzzz_1, tg_yyyzz_yyzzz_1, tg_yyyzz_yyzzzz_0, tg_yyyzz_yyzzzz_1, \
                                         tg_yyyzz_yzzzz_1, tg_yyyzz_yzzzzz_0, tg_yyyzz_yzzzzz_1, tg_yyyzz_zzzzz_1, \
                                         tg_yyyzz_zzzzzz_0, tg_yyyzz_zzzzzz_1, tg_yyzzz_xxxxx_1, tg_yyzzz_xxxxxx_0, \
                                         tg_yyzzz_xxxxxx_1, tg_yyzzz_xxxxxy_0, tg_yyzzz_xxxxxy_1, tg_yyzzz_xxxxxz_0, \
                                         tg_yyzzz_xxxxxz_1, tg_yyzzz_xxxxy_1, tg_yyzzz_xxxxyy_0, tg_yyzzz_xxxxyy_1, \
                                         tg_yyzzz_xxxxyz_0, tg_yyzzz_xxxxyz_1, tg_yyzzz_xxxxz_1, tg_yyzzz_xxxxzz_0, \
                                         tg_yyzzz_xxxxzz_1, tg_yyzzz_xxxyy_1, tg_yyzzz_xxxyyy_0, tg_yyzzz_xxxyyy_1, \
                                         tg_yyzzz_xxxyyz_0, tg_yyzzz_xxxyyz_1, tg_yyzzz_xxxyz_1, tg_yyzzz_xxxyzz_0, \
                                         tg_yyzzz_xxxyzz_1, tg_yyzzz_xxxzz_1, tg_yyzzz_xxxzzz_0, tg_yyzzz_xxxzzz_1, \
                                         tg_yyzzz_xxyyy_1, tg_yyzzz_xxyyyy_0, tg_yyzzz_xxyyyy_1, tg_yyzzz_xxyyyz_0, \
                                         tg_yyzzz_xxyyyz_1, tg_yyzzz_xxyyz_1, tg_yyzzz_xxyyzz_0, tg_yyzzz_xxyyzz_1, \
                                         tg_yyzzz_xxyzz_1, tg_yyzzz_xxyzzz_0, tg_yyzzz_xxyzzz_1, tg_yyzzz_xxzzz_1, \
                                         tg_yyzzz_xxzzzz_0, tg_yyzzz_xxzzzz_1, tg_yyzzz_xyyyy_1, tg_yyzzz_xyyyyy_0, \
                                         tg_yyzzz_xyyyyy_1, tg_yyzzz_xyyyyz_0, tg_yyzzz_xyyyyz_1, tg_yyzzz_xyyyz_1, \
                                         tg_yyzzz_xyyyzz_0, tg_yyzzz_xyyyzz_1, tg_yyzzz_xyyzz_1, tg_yyzzz_xyyzzz_0, \
                                         tg_yyzzz_xyyzzz_1, tg_yyzzz_xyzzz_1, tg_yyzzz_xyzzzz_0, tg_yyzzz_xyzzzz_1, \
                                         tg_yyzzz_xzzzz_1, tg_yyzzz_xzzzzz_0, tg_yyzzz_xzzzzz_1, tg_yyzzz_yyyyy_1, \
                                         tg_yyzzz_yyyyyy_0, tg_yyzzz_yyyyyy_1, tg_yyzzz_yyyyyz_0, tg_yyzzz_yyyyyz_1, \
                                         tg_yyzzz_yyyyz_1, tg_yyzzz_yyyyzz_0, tg_yyzzz_yyyyzz_1, tg_yyzzz_yyyzz_1, \
                                         tg_yyzzz_yyyzzz_0, tg_yyzzz_yyyzzz_1, tg_yyzzz_yyzzz_1, tg_yyzzz_yyzzzz_0, \
                                         tg_yyzzz_yyzzzz_1, tg_yyzzz_yzzzz_1, tg_yyzzz_yzzzzz_0, tg_yyzzz_yzzzzz_1, \
                                         tg_yyzzz_zzzzz_1, tg_yyzzz_zzzzzz_0, tg_yyzzz_zzzzzz_1, tg_yzzzz_xxxxx_1, \
                                         tg_yzzzz_xxxxxx_0, tg_yzzzz_xxxxxx_1, tg_yzzzz_xxxxxy_0, tg_yzzzz_xxxxxy_1, \
                                         tg_yzzzz_xxxxxz_0, tg_yzzzz_xxxxxz_1, tg_yzzzz_xxxxy_1, tg_yzzzz_xxxxyy_0, \
                                         tg_yzzzz_xxxxyy_1, tg_yzzzz_xxxxyz_0, tg_yzzzz_xxxxyz_1, tg_yzzzz_xxxxz_1, \
                                         tg_yzzzz_xxxxzz_0, tg_yzzzz_xxxxzz_1, tg_yzzzz_xxxyy_1, tg_yzzzz_xxxyyy_0, \
                                         tg_yzzzz_xxxyyy_1, tg_yzzzz_xxxyyz_0, tg_yzzzz_xxxyyz_1, tg_yzzzz_xxxyz_1, \
                                         tg_yzzzz_xxxyzz_0, tg_yzzzz_xxxyzz_1, tg_yzzzz_xxxzz_1, tg_yzzzz_xxxzzz_0, \
                                         tg_yzzzz_xxxzzz_1, tg_yzzzz_xxyyy_1, tg_yzzzz_xxyyyy_0, tg_yzzzz_xxyyyy_1, \
                                         tg_yzzzz_xxyyyz_0, tg_yzzzz_xxyyyz_1, tg_yzzzz_xxyyz_1, tg_yzzzz_xxyyzz_0, \
                                         tg_yzzzz_xxyyzz_1, tg_yzzzz_xxyzz_1, tg_yzzzz_xxyzzz_0, tg_yzzzz_xxyzzz_1, \
                                         tg_yzzzz_xxzzz_1, tg_yzzzz_xxzzzz_0, tg_yzzzz_xxzzzz_1, tg_yzzzz_xyyyy_1, \
                                         tg_yzzzz_xyyyyy_0, tg_yzzzz_xyyyyy_1, tg_yzzzz_xyyyyz_0, tg_yzzzz_xyyyyz_1, \
                                         tg_yzzzz_xyyyz_1, tg_yzzzz_xyyyzz_0, tg_yzzzz_xyyyzz_1, tg_yzzzz_xyyzz_1, \
                                         tg_yzzzz_xyyzzz_0, tg_yzzzz_xyyzzz_1, tg_yzzzz_xyzzz_1, tg_yzzzz_xyzzzz_0, \
                                         tg_yzzzz_xyzzzz_1, tg_yzzzz_xzzzz_1, tg_yzzzz_xzzzzz_0, tg_yzzzz_xzzzzz_1, \
                                         tg_yzzzz_yyyyy_1, tg_yzzzz_yyyyyy_0, tg_yzzzz_yyyyyy_1, tg_yzzzz_yyyyyz_0, \
                                         tg_yzzzz_yyyyyz_1, tg_yzzzz_yyyyz_1, tg_yzzzz_yyyyzz_0, tg_yzzzz_yyyyzz_1, \
                                         tg_yzzzz_yyyzz_1, tg_yzzzz_yyyzzz_0, tg_yzzzz_yyyzzz_1, tg_yzzzz_yyzzz_1, \
                                         tg_yzzzz_yyzzzz_0, tg_yzzzz_yyzzzz_1, tg_yzzzz_yzzzz_1, tg_yzzzz_yzzzzz_0, \
                                         tg_yzzzz_yzzzzz_1, tg_yzzzz_zzzzz_1, tg_yzzzz_zzzzzz_0, tg_yzzzz_zzzzzz_1, \
                                         tg_zzzzz_xxxxx_1, tg_zzzzz_xxxxxx_0, tg_zzzzz_xxxxxx_1, tg_zzzzz_xxxxxy_0, \
                                         tg_zzzzz_xxxxxy_1, tg_zzzzz_xxxxxz_0, tg_zzzzz_xxxxxz_1, tg_zzzzz_xxxxy_1, \
                                         tg_zzzzz_xxxxyy_0, tg_zzzzz_xxxxyy_1, tg_zzzzz_xxxxyz_0, tg_zzzzz_xxxxyz_1, \
                                         tg_zzzzz_xxxxz_1, tg_zzzzz_xxxxzz_0, tg_zzzzz_xxxxzz_1, tg_zzzzz_xxxyy_1, \
                                         tg_zzzzz_xxxyyy_0, tg_zzzzz_xxxyyy_1, tg_zzzzz_xxxyyz_0, tg_zzzzz_xxxyyz_1, \
                                         tg_zzzzz_xxxyz_1, tg_zzzzz_xxxyzz_0, tg_zzzzz_xxxyzz_1, tg_zzzzz_xxxzz_1, \
                                         tg_zzzzz_xxxzzz_0, tg_zzzzz_xxxzzz_1, tg_zzzzz_xxyyy_1, tg_zzzzz_xxyyyy_0, \
                                         tg_zzzzz_xxyyyy_1, tg_zzzzz_xxyyyz_0, tg_zzzzz_xxyyyz_1, tg_zzzzz_xxyyz_1, \
                                         tg_zzzzz_xxyyzz_0, tg_zzzzz_xxyyzz_1, tg_zzzzz_xxyzz_1, tg_zzzzz_xxyzzz_0, \
                                         tg_zzzzz_xxyzzz_1, tg_zzzzz_xxzzz_1, tg_zzzzz_xxzzzz_0, tg_zzzzz_xxzzzz_1, \
                                         tg_zzzzz_xyyyy_1, tg_zzzzz_xyyyyy_0, tg_zzzzz_xyyyyy_1, tg_zzzzz_xyyyyz_0, \
                                         tg_zzzzz_xyyyyz_1, tg_zzzzz_xyyyz_1, tg_zzzzz_xyyyzz_0, tg_zzzzz_xyyyzz_1, \
                                         tg_zzzzz_xyyzz_1, tg_zzzzz_xyyzzz_0, tg_zzzzz_xyyzzz_1, tg_zzzzz_xyzzz_1, \
                                         tg_zzzzz_xyzzzz_0, tg_zzzzz_xyzzzz_1, tg_zzzzz_xzzzz_1, tg_zzzzz_xzzzzz_0, \
                                         tg_zzzzz_xzzzzz_1, tg_zzzzz_yyyyy_1, tg_zzzzz_yyyyyy_0, tg_zzzzz_yyyyyy_1, \
                                         tg_zzzzz_yyyyyz_0, tg_zzzzz_yyyyyz_1, tg_zzzzz_yyyyz_1, tg_zzzzz_yyyyzz_0, \
                                         tg_zzzzz_yyyyzz_1, tg_zzzzz_yyyzz_1, tg_zzzzz_yyyzzz_0, tg_zzzzz_yyyzzz_1, \
                                         tg_zzzzz_yyzzz_1, tg_zzzzz_yyzzzz_0, tg_zzzzz_yyzzzz_1, tg_zzzzz_yzzzz_1, \
                                         tg_zzzzz_yzzzzz_0, tg_zzzzz_yzzzzz_1, tg_zzzzz_zzzzz_1, tg_zzzzz_zzzzzz_0, \
                                         tg_zzzzz_zzzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fxn = fxn[j];

                    double fr = wp_x[j]; 

                    tg_xyyyzz_xxzzzz_0[j] = pb_x * tg_yyyzz_xxzzzz_0[j] + fr * tg_yyyzz_xxzzzz_1[j] + fl1_fxn * tg_yyyzz_xzzzz_1[j];

                    tg_xyyyzz_xyyyyy_0[j] = pb_x * tg_yyyzz_xyyyyy_0[j] + fr * tg_yyyzz_xyyyyy_1[j] + 0.5 * fl1_fxn * tg_yyyzz_yyyyy_1[j];

                    tg_xyyyzz_xyyyyz_0[j] = pb_x * tg_yyyzz_xyyyyz_0[j] + fr * tg_yyyzz_xyyyyz_1[j] + 0.5 * fl1_fxn * tg_yyyzz_yyyyz_1[j];

                    tg_xyyyzz_xyyyzz_0[j] = pb_x * tg_yyyzz_xyyyzz_0[j] + fr * tg_yyyzz_xyyyzz_1[j] + 0.5 * fl1_fxn * tg_yyyzz_yyyzz_1[j];

                    tg_xyyyzz_xyyzzz_0[j] = pb_x * tg_yyyzz_xyyzzz_0[j] + fr * tg_yyyzz_xyyzzz_1[j] + 0.5 * fl1_fxn * tg_yyyzz_yyzzz_1[j];

                    tg_xyyyzz_xyzzzz_0[j] = pb_x * tg_yyyzz_xyzzzz_0[j] + fr * tg_yyyzz_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyzz_yzzzz_1[j];

                    tg_xyyyzz_xzzzzz_0[j] = pb_x * tg_yyyzz_xzzzzz_0[j] + fr * tg_yyyzz_xzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyzz_zzzzz_1[j];

                    tg_xyyyzz_yyyyyy_0[j] = pb_x * tg_yyyzz_yyyyyy_0[j] + fr * tg_yyyzz_yyyyyy_1[j];

                    tg_xyyyzz_yyyyyz_0[j] = pb_x * tg_yyyzz_yyyyyz_0[j] + fr * tg_yyyzz_yyyyyz_1[j];

                    tg_xyyyzz_yyyyzz_0[j] = pb_x * tg_yyyzz_yyyyzz_0[j] + fr * tg_yyyzz_yyyyzz_1[j];

                    tg_xyyyzz_yyyzzz_0[j] = pb_x * tg_yyyzz_yyyzzz_0[j] + fr * tg_yyyzz_yyyzzz_1[j];

                    tg_xyyyzz_yyzzzz_0[j] = pb_x * tg_yyyzz_yyzzzz_0[j] + fr * tg_yyyzz_yyzzzz_1[j];

                    tg_xyyyzz_yzzzzz_0[j] = pb_x * tg_yyyzz_yzzzzz_0[j] + fr * tg_yyyzz_yzzzzz_1[j];

                    tg_xyyyzz_zzzzzz_0[j] = pb_x * tg_yyyzz_zzzzzz_0[j] + fr * tg_yyyzz_zzzzzz_1[j];

                    tg_xyyzzz_xxxxxx_0[j] = pb_x * tg_yyzzz_xxxxxx_0[j] + fr * tg_yyzzz_xxxxxx_1[j] + 3.0 * fl1_fxn * tg_yyzzz_xxxxx_1[j];

                    tg_xyyzzz_xxxxxy_0[j] = pb_x * tg_yyzzz_xxxxxy_0[j] + fr * tg_yyzzz_xxxxxy_1[j] + 2.5 * fl1_fxn * tg_yyzzz_xxxxy_1[j];

                    tg_xyyzzz_xxxxxz_0[j] = pb_x * tg_yyzzz_xxxxxz_0[j] + fr * tg_yyzzz_xxxxxz_1[j] + 2.5 * fl1_fxn * tg_yyzzz_xxxxz_1[j];

                    tg_xyyzzz_xxxxyy_0[j] = pb_x * tg_yyzzz_xxxxyy_0[j] + fr * tg_yyzzz_xxxxyy_1[j] + 2.0 * fl1_fxn * tg_yyzzz_xxxyy_1[j];

                    tg_xyyzzz_xxxxyz_0[j] = pb_x * tg_yyzzz_xxxxyz_0[j] + fr * tg_yyzzz_xxxxyz_1[j] + 2.0 * fl1_fxn * tg_yyzzz_xxxyz_1[j];

                    tg_xyyzzz_xxxxzz_0[j] = pb_x * tg_yyzzz_xxxxzz_0[j] + fr * tg_yyzzz_xxxxzz_1[j] + 2.0 * fl1_fxn * tg_yyzzz_xxxzz_1[j];

                    tg_xyyzzz_xxxyyy_0[j] = pb_x * tg_yyzzz_xxxyyy_0[j] + fr * tg_yyzzz_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_yyzzz_xxyyy_1[j];

                    tg_xyyzzz_xxxyyz_0[j] = pb_x * tg_yyzzz_xxxyyz_0[j] + fr * tg_yyzzz_xxxyyz_1[j] + 1.5 * fl1_fxn * tg_yyzzz_xxyyz_1[j];

                    tg_xyyzzz_xxxyzz_0[j] = pb_x * tg_yyzzz_xxxyzz_0[j] + fr * tg_yyzzz_xxxyzz_1[j] + 1.5 * fl1_fxn * tg_yyzzz_xxyzz_1[j];

                    tg_xyyzzz_xxxzzz_0[j] = pb_x * tg_yyzzz_xxxzzz_0[j] + fr * tg_yyzzz_xxxzzz_1[j] + 1.5 * fl1_fxn * tg_yyzzz_xxzzz_1[j];

                    tg_xyyzzz_xxyyyy_0[j] = pb_x * tg_yyzzz_xxyyyy_0[j] + fr * tg_yyzzz_xxyyyy_1[j] + fl1_fxn * tg_yyzzz_xyyyy_1[j];

                    tg_xyyzzz_xxyyyz_0[j] = pb_x * tg_yyzzz_xxyyyz_0[j] + fr * tg_yyzzz_xxyyyz_1[j] + fl1_fxn * tg_yyzzz_xyyyz_1[j];

                    tg_xyyzzz_xxyyzz_0[j] = pb_x * tg_yyzzz_xxyyzz_0[j] + fr * tg_yyzzz_xxyyzz_1[j] + fl1_fxn * tg_yyzzz_xyyzz_1[j];

                    tg_xyyzzz_xxyzzz_0[j] = pb_x * tg_yyzzz_xxyzzz_0[j] + fr * tg_yyzzz_xxyzzz_1[j] + fl1_fxn * tg_yyzzz_xyzzz_1[j];

                    tg_xyyzzz_xxzzzz_0[j] = pb_x * tg_yyzzz_xxzzzz_0[j] + fr * tg_yyzzz_xxzzzz_1[j] + fl1_fxn * tg_yyzzz_xzzzz_1[j];

                    tg_xyyzzz_xyyyyy_0[j] = pb_x * tg_yyzzz_xyyyyy_0[j] + fr * tg_yyzzz_xyyyyy_1[j] + 0.5 * fl1_fxn * tg_yyzzz_yyyyy_1[j];

                    tg_xyyzzz_xyyyyz_0[j] = pb_x * tg_yyzzz_xyyyyz_0[j] + fr * tg_yyzzz_xyyyyz_1[j] + 0.5 * fl1_fxn * tg_yyzzz_yyyyz_1[j];

                    tg_xyyzzz_xyyyzz_0[j] = pb_x * tg_yyzzz_xyyyzz_0[j] + fr * tg_yyzzz_xyyyzz_1[j] + 0.5 * fl1_fxn * tg_yyzzz_yyyzz_1[j];

                    tg_xyyzzz_xyyzzz_0[j] = pb_x * tg_yyzzz_xyyzzz_0[j] + fr * tg_yyzzz_xyyzzz_1[j] + 0.5 * fl1_fxn * tg_yyzzz_yyzzz_1[j];

                    tg_xyyzzz_xyzzzz_0[j] = pb_x * tg_yyzzz_xyzzzz_0[j] + fr * tg_yyzzz_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_yyzzz_yzzzz_1[j];

                    tg_xyyzzz_xzzzzz_0[j] = pb_x * tg_yyzzz_xzzzzz_0[j] + fr * tg_yyzzz_xzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyzzz_zzzzz_1[j];

                    tg_xyyzzz_yyyyyy_0[j] = pb_x * tg_yyzzz_yyyyyy_0[j] + fr * tg_yyzzz_yyyyyy_1[j];

                    tg_xyyzzz_yyyyyz_0[j] = pb_x * tg_yyzzz_yyyyyz_0[j] + fr * tg_yyzzz_yyyyyz_1[j];

                    tg_xyyzzz_yyyyzz_0[j] = pb_x * tg_yyzzz_yyyyzz_0[j] + fr * tg_yyzzz_yyyyzz_1[j];

                    tg_xyyzzz_yyyzzz_0[j] = pb_x * tg_yyzzz_yyyzzz_0[j] + fr * tg_yyzzz_yyyzzz_1[j];

                    tg_xyyzzz_yyzzzz_0[j] = pb_x * tg_yyzzz_yyzzzz_0[j] + fr * tg_yyzzz_yyzzzz_1[j];

                    tg_xyyzzz_yzzzzz_0[j] = pb_x * tg_yyzzz_yzzzzz_0[j] + fr * tg_yyzzz_yzzzzz_1[j];

                    tg_xyyzzz_zzzzzz_0[j] = pb_x * tg_yyzzz_zzzzzz_0[j] + fr * tg_yyzzz_zzzzzz_1[j];

                    tg_xyzzzz_xxxxxx_0[j] = pb_x * tg_yzzzz_xxxxxx_0[j] + fr * tg_yzzzz_xxxxxx_1[j] + 3.0 * fl1_fxn * tg_yzzzz_xxxxx_1[j];

                    tg_xyzzzz_xxxxxy_0[j] = pb_x * tg_yzzzz_xxxxxy_0[j] + fr * tg_yzzzz_xxxxxy_1[j] + 2.5 * fl1_fxn * tg_yzzzz_xxxxy_1[j];

                    tg_xyzzzz_xxxxxz_0[j] = pb_x * tg_yzzzz_xxxxxz_0[j] + fr * tg_yzzzz_xxxxxz_1[j] + 2.5 * fl1_fxn * tg_yzzzz_xxxxz_1[j];

                    tg_xyzzzz_xxxxyy_0[j] = pb_x * tg_yzzzz_xxxxyy_0[j] + fr * tg_yzzzz_xxxxyy_1[j] + 2.0 * fl1_fxn * tg_yzzzz_xxxyy_1[j];

                    tg_xyzzzz_xxxxyz_0[j] = pb_x * tg_yzzzz_xxxxyz_0[j] + fr * tg_yzzzz_xxxxyz_1[j] + 2.0 * fl1_fxn * tg_yzzzz_xxxyz_1[j];

                    tg_xyzzzz_xxxxzz_0[j] = pb_x * tg_yzzzz_xxxxzz_0[j] + fr * tg_yzzzz_xxxxzz_1[j] + 2.0 * fl1_fxn * tg_yzzzz_xxxzz_1[j];

                    tg_xyzzzz_xxxyyy_0[j] = pb_x * tg_yzzzz_xxxyyy_0[j] + fr * tg_yzzzz_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_yzzzz_xxyyy_1[j];

                    tg_xyzzzz_xxxyyz_0[j] = pb_x * tg_yzzzz_xxxyyz_0[j] + fr * tg_yzzzz_xxxyyz_1[j] + 1.5 * fl1_fxn * tg_yzzzz_xxyyz_1[j];

                    tg_xyzzzz_xxxyzz_0[j] = pb_x * tg_yzzzz_xxxyzz_0[j] + fr * tg_yzzzz_xxxyzz_1[j] + 1.5 * fl1_fxn * tg_yzzzz_xxyzz_1[j];

                    tg_xyzzzz_xxxzzz_0[j] = pb_x * tg_yzzzz_xxxzzz_0[j] + fr * tg_yzzzz_xxxzzz_1[j] + 1.5 * fl1_fxn * tg_yzzzz_xxzzz_1[j];

                    tg_xyzzzz_xxyyyy_0[j] = pb_x * tg_yzzzz_xxyyyy_0[j] + fr * tg_yzzzz_xxyyyy_1[j] + fl1_fxn * tg_yzzzz_xyyyy_1[j];

                    tg_xyzzzz_xxyyyz_0[j] = pb_x * tg_yzzzz_xxyyyz_0[j] + fr * tg_yzzzz_xxyyyz_1[j] + fl1_fxn * tg_yzzzz_xyyyz_1[j];

                    tg_xyzzzz_xxyyzz_0[j] = pb_x * tg_yzzzz_xxyyzz_0[j] + fr * tg_yzzzz_xxyyzz_1[j] + fl1_fxn * tg_yzzzz_xyyzz_1[j];

                    tg_xyzzzz_xxyzzz_0[j] = pb_x * tg_yzzzz_xxyzzz_0[j] + fr * tg_yzzzz_xxyzzz_1[j] + fl1_fxn * tg_yzzzz_xyzzz_1[j];

                    tg_xyzzzz_xxzzzz_0[j] = pb_x * tg_yzzzz_xxzzzz_0[j] + fr * tg_yzzzz_xxzzzz_1[j] + fl1_fxn * tg_yzzzz_xzzzz_1[j];

                    tg_xyzzzz_xyyyyy_0[j] = pb_x * tg_yzzzz_xyyyyy_0[j] + fr * tg_yzzzz_xyyyyy_1[j] + 0.5 * fl1_fxn * tg_yzzzz_yyyyy_1[j];

                    tg_xyzzzz_xyyyyz_0[j] = pb_x * tg_yzzzz_xyyyyz_0[j] + fr * tg_yzzzz_xyyyyz_1[j] + 0.5 * fl1_fxn * tg_yzzzz_yyyyz_1[j];

                    tg_xyzzzz_xyyyzz_0[j] = pb_x * tg_yzzzz_xyyyzz_0[j] + fr * tg_yzzzz_xyyyzz_1[j] + 0.5 * fl1_fxn * tg_yzzzz_yyyzz_1[j];

                    tg_xyzzzz_xyyzzz_0[j] = pb_x * tg_yzzzz_xyyzzz_0[j] + fr * tg_yzzzz_xyyzzz_1[j] + 0.5 * fl1_fxn * tg_yzzzz_yyzzz_1[j];

                    tg_xyzzzz_xyzzzz_0[j] = pb_x * tg_yzzzz_xyzzzz_0[j] + fr * tg_yzzzz_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_yzzzz_yzzzz_1[j];

                    tg_xyzzzz_xzzzzz_0[j] = pb_x * tg_yzzzz_xzzzzz_0[j] + fr * tg_yzzzz_xzzzzz_1[j] + 0.5 * fl1_fxn * tg_yzzzz_zzzzz_1[j];

                    tg_xyzzzz_yyyyyy_0[j] = pb_x * tg_yzzzz_yyyyyy_0[j] + fr * tg_yzzzz_yyyyyy_1[j];

                    tg_xyzzzz_yyyyyz_0[j] = pb_x * tg_yzzzz_yyyyyz_0[j] + fr * tg_yzzzz_yyyyyz_1[j];

                    tg_xyzzzz_yyyyzz_0[j] = pb_x * tg_yzzzz_yyyyzz_0[j] + fr * tg_yzzzz_yyyyzz_1[j];

                    tg_xyzzzz_yyyzzz_0[j] = pb_x * tg_yzzzz_yyyzzz_0[j] + fr * tg_yzzzz_yyyzzz_1[j];

                    tg_xyzzzz_yyzzzz_0[j] = pb_x * tg_yzzzz_yyzzzz_0[j] + fr * tg_yzzzz_yyzzzz_1[j];

                    tg_xyzzzz_yzzzzz_0[j] = pb_x * tg_yzzzz_yzzzzz_0[j] + fr * tg_yzzzz_yzzzzz_1[j];

                    tg_xyzzzz_zzzzzz_0[j] = pb_x * tg_yzzzz_zzzzzz_0[j] + fr * tg_yzzzz_zzzzzz_1[j];

                    tg_xzzzzz_xxxxxx_0[j] = pb_x * tg_zzzzz_xxxxxx_0[j] + fr * tg_zzzzz_xxxxxx_1[j] + 3.0 * fl1_fxn * tg_zzzzz_xxxxx_1[j];

                    tg_xzzzzz_xxxxxy_0[j] = pb_x * tg_zzzzz_xxxxxy_0[j] + fr * tg_zzzzz_xxxxxy_1[j] + 2.5 * fl1_fxn * tg_zzzzz_xxxxy_1[j];

                    tg_xzzzzz_xxxxxz_0[j] = pb_x * tg_zzzzz_xxxxxz_0[j] + fr * tg_zzzzz_xxxxxz_1[j] + 2.5 * fl1_fxn * tg_zzzzz_xxxxz_1[j];

                    tg_xzzzzz_xxxxyy_0[j] = pb_x * tg_zzzzz_xxxxyy_0[j] + fr * tg_zzzzz_xxxxyy_1[j] + 2.0 * fl1_fxn * tg_zzzzz_xxxyy_1[j];

                    tg_xzzzzz_xxxxyz_0[j] = pb_x * tg_zzzzz_xxxxyz_0[j] + fr * tg_zzzzz_xxxxyz_1[j] + 2.0 * fl1_fxn * tg_zzzzz_xxxyz_1[j];

                    tg_xzzzzz_xxxxzz_0[j] = pb_x * tg_zzzzz_xxxxzz_0[j] + fr * tg_zzzzz_xxxxzz_1[j] + 2.0 * fl1_fxn * tg_zzzzz_xxxzz_1[j];

                    tg_xzzzzz_xxxyyy_0[j] = pb_x * tg_zzzzz_xxxyyy_0[j] + fr * tg_zzzzz_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_zzzzz_xxyyy_1[j];

                    tg_xzzzzz_xxxyyz_0[j] = pb_x * tg_zzzzz_xxxyyz_0[j] + fr * tg_zzzzz_xxxyyz_1[j] + 1.5 * fl1_fxn * tg_zzzzz_xxyyz_1[j];

                    tg_xzzzzz_xxxyzz_0[j] = pb_x * tg_zzzzz_xxxyzz_0[j] + fr * tg_zzzzz_xxxyzz_1[j] + 1.5 * fl1_fxn * tg_zzzzz_xxyzz_1[j];

                    tg_xzzzzz_xxxzzz_0[j] = pb_x * tg_zzzzz_xxxzzz_0[j] + fr * tg_zzzzz_xxxzzz_1[j] + 1.5 * fl1_fxn * tg_zzzzz_xxzzz_1[j];

                    tg_xzzzzz_xxyyyy_0[j] = pb_x * tg_zzzzz_xxyyyy_0[j] + fr * tg_zzzzz_xxyyyy_1[j] + fl1_fxn * tg_zzzzz_xyyyy_1[j];

                    tg_xzzzzz_xxyyyz_0[j] = pb_x * tg_zzzzz_xxyyyz_0[j] + fr * tg_zzzzz_xxyyyz_1[j] + fl1_fxn * tg_zzzzz_xyyyz_1[j];

                    tg_xzzzzz_xxyyzz_0[j] = pb_x * tg_zzzzz_xxyyzz_0[j] + fr * tg_zzzzz_xxyyzz_1[j] + fl1_fxn * tg_zzzzz_xyyzz_1[j];

                    tg_xzzzzz_xxyzzz_0[j] = pb_x * tg_zzzzz_xxyzzz_0[j] + fr * tg_zzzzz_xxyzzz_1[j] + fl1_fxn * tg_zzzzz_xyzzz_1[j];

                    tg_xzzzzz_xxzzzz_0[j] = pb_x * tg_zzzzz_xxzzzz_0[j] + fr * tg_zzzzz_xxzzzz_1[j] + fl1_fxn * tg_zzzzz_xzzzz_1[j];

                    tg_xzzzzz_xyyyyy_0[j] = pb_x * tg_zzzzz_xyyyyy_0[j] + fr * tg_zzzzz_xyyyyy_1[j] + 0.5 * fl1_fxn * tg_zzzzz_yyyyy_1[j];

                    tg_xzzzzz_xyyyyz_0[j] = pb_x * tg_zzzzz_xyyyyz_0[j] + fr * tg_zzzzz_xyyyyz_1[j] + 0.5 * fl1_fxn * tg_zzzzz_yyyyz_1[j];

                    tg_xzzzzz_xyyyzz_0[j] = pb_x * tg_zzzzz_xyyyzz_0[j] + fr * tg_zzzzz_xyyyzz_1[j] + 0.5 * fl1_fxn * tg_zzzzz_yyyzz_1[j];

                    tg_xzzzzz_xyyzzz_0[j] = pb_x * tg_zzzzz_xyyzzz_0[j] + fr * tg_zzzzz_xyyzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzz_yyzzz_1[j];

                    tg_xzzzzz_xyzzzz_0[j] = pb_x * tg_zzzzz_xyzzzz_0[j] + fr * tg_zzzzz_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzz_yzzzz_1[j];

                    tg_xzzzzz_xzzzzz_0[j] = pb_x * tg_zzzzz_xzzzzz_0[j] + fr * tg_zzzzz_xzzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzz_zzzzz_1[j];

                    tg_xzzzzz_yyyyyy_0[j] = pb_x * tg_zzzzz_yyyyyy_0[j] + fr * tg_zzzzz_yyyyyy_1[j];

                    tg_xzzzzz_yyyyyz_0[j] = pb_x * tg_zzzzz_yyyyyz_0[j] + fr * tg_zzzzz_yyyyyz_1[j];

                    tg_xzzzzz_yyyyzz_0[j] = pb_x * tg_zzzzz_yyyyzz_0[j] + fr * tg_zzzzz_yyyyzz_1[j];

                    tg_xzzzzz_yyyzzz_0[j] = pb_x * tg_zzzzz_yyyzzz_0[j] + fr * tg_zzzzz_yyyzzz_1[j];

                    tg_xzzzzz_yyzzzz_0[j] = pb_x * tg_zzzzz_yyzzzz_0[j] + fr * tg_zzzzz_yyzzzz_1[j];

                    tg_xzzzzz_yzzzzz_0[j] = pb_x * tg_zzzzz_yzzzzz_0[j] + fr * tg_zzzzz_yzzzzz_1[j];

                    tg_xzzzzz_zzzzzz_0[j] = pb_x * tg_zzzzz_zzzzzz_0[j] + fr * tg_zzzzz_zzzzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSISI_588_686(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (588,686)

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
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_6_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_6_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_5_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
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

                auto tg_yyyyy_xxxxxx_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 420); 

                auto tg_yyyyy_xxxxxy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 421); 

                auto tg_yyyyy_xxxxxz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 422); 

                auto tg_yyyyy_xxxxyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 423); 

                auto tg_yyyyy_xxxxyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 424); 

                auto tg_yyyyy_xxxxzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 425); 

                auto tg_yyyyy_xxxyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 426); 

                auto tg_yyyyy_xxxyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 427); 

                auto tg_yyyyy_xxxyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 428); 

                auto tg_yyyyy_xxxzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 429); 

                auto tg_yyyyy_xxyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 430); 

                auto tg_yyyyy_xxyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 431); 

                auto tg_yyyyy_xxyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 432); 

                auto tg_yyyyy_xxyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 433); 

                auto tg_yyyyy_xxzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 434); 

                auto tg_yyyyy_xyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 435); 

                auto tg_yyyyy_xyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 436); 

                auto tg_yyyyy_xyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 437); 

                auto tg_yyyyy_xyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 438); 

                auto tg_yyyyy_xyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 439); 

                auto tg_yyyyy_xzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 440); 

                auto tg_yyyyy_yyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 441); 

                auto tg_yyyyy_yyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 442); 

                auto tg_yyyyy_yyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 443); 

                auto tg_yyyyy_yyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 444); 

                auto tg_yyyyy_yyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 445); 

                auto tg_yyyyy_yzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 446); 

                auto tg_yyyyy_zzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 447); 

                auto tg_yyyyz_xxxxxx_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 448); 

                auto tg_yyyyz_xxxxxy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 449); 

                auto tg_yyyyz_xxxxxz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 450); 

                auto tg_yyyyz_xxxxyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 451); 

                auto tg_yyyyz_xxxxyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 452); 

                auto tg_yyyyz_xxxxzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 453); 

                auto tg_yyyyz_xxxyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 454); 

                auto tg_yyyyz_xxxyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 455); 

                auto tg_yyyyz_xxxyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 456); 

                auto tg_yyyyz_xxxzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 457); 

                auto tg_yyyyz_xxyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 458); 

                auto tg_yyyyz_xxyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 459); 

                auto tg_yyyyz_xxyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 460); 

                auto tg_yyyyz_xxyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 461); 

                auto tg_yyyyz_xxzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 462); 

                auto tg_yyyyz_xyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 463); 

                auto tg_yyyyz_xyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 464); 

                auto tg_yyyyz_xyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 465); 

                auto tg_yyyyz_xyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 466); 

                auto tg_yyyyz_xyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 467); 

                auto tg_yyyyz_xzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 468); 

                auto tg_yyyyz_yyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 469); 

                auto tg_yyyyz_yyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 470); 

                auto tg_yyyyz_yyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 471); 

                auto tg_yyyyz_yyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 472); 

                auto tg_yyyyz_yyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 473); 

                auto tg_yyyyz_yzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 474); 

                auto tg_yyyyz_zzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 475); 

                auto tg_yyyzz_xxxxxx_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 476); 

                auto tg_yyyzz_xxxxxy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 477); 

                auto tg_yyyzz_xxxxxz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 478); 

                auto tg_yyyzz_xxxxyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 479); 

                auto tg_yyyzz_xxxxyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 480); 

                auto tg_yyyzz_xxxxzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 481); 

                auto tg_yyyzz_xxxyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 482); 

                auto tg_yyyzz_xxxyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 483); 

                auto tg_yyyzz_xxxyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 484); 

                auto tg_yyyzz_xxxzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 485); 

                auto tg_yyyzz_xxyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 486); 

                auto tg_yyyzz_xxyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 487); 

                auto tg_yyyzz_xxyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 488); 

                auto tg_yyyzz_xxyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 489); 

                auto tg_yyyzz_xxzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 490); 

                auto tg_yyyzz_xyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 491); 

                auto tg_yyyzz_xyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 492); 

                auto tg_yyyzz_xyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 493); 

                auto tg_yyyzz_xyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 494); 

                auto tg_yyyzz_xyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 495); 

                auto tg_yyyzz_xzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 496); 

                auto tg_yyyzz_yyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 497); 

                auto tg_yyyzz_yyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 498); 

                auto tg_yyyzz_yyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 499); 

                auto tg_yyyzz_yyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 500); 

                auto tg_yyyzz_yyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 501); 

                auto tg_yyyzz_yzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 502); 

                auto tg_yyyzz_zzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 503); 

                auto tg_yyzzz_xxxxxx_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 504); 

                auto tg_yyzzz_xxxxxy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 505); 

                auto tg_yyzzz_xxxxxz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 506); 

                auto tg_yyzzz_xxxxyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 507); 

                auto tg_yyzzz_xxxxyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 508); 

                auto tg_yyzzz_xxxxzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 509); 

                auto tg_yyzzz_xxxyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 510); 

                auto tg_yyzzz_xxxyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 511); 

                auto tg_yyzzz_xxxyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 512); 

                auto tg_yyzzz_xxxzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 513); 

                auto tg_yyzzz_xxyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 514); 

                auto tg_yyzzz_xxyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 515); 

                auto tg_yyzzz_xxyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 516); 

                auto tg_yyzzz_xxyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 517); 

                auto tg_yyyyy_xxxxxx_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 420); 

                auto tg_yyyyy_xxxxxy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 421); 

                auto tg_yyyyy_xxxxxz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 422); 

                auto tg_yyyyy_xxxxyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 423); 

                auto tg_yyyyy_xxxxyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 424); 

                auto tg_yyyyy_xxxxzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 425); 

                auto tg_yyyyy_xxxyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 426); 

                auto tg_yyyyy_xxxyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 427); 

                auto tg_yyyyy_xxxyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 428); 

                auto tg_yyyyy_xxxzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 429); 

                auto tg_yyyyy_xxyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 430); 

                auto tg_yyyyy_xxyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 431); 

                auto tg_yyyyy_xxyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 432); 

                auto tg_yyyyy_xxyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 433); 

                auto tg_yyyyy_xxzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 434); 

                auto tg_yyyyy_xyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 435); 

                auto tg_yyyyy_xyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 436); 

                auto tg_yyyyy_xyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 437); 

                auto tg_yyyyy_xyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 438); 

                auto tg_yyyyy_xyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 439); 

                auto tg_yyyyy_xzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 440); 

                auto tg_yyyyy_yyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 441); 

                auto tg_yyyyy_yyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 442); 

                auto tg_yyyyy_yyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 443); 

                auto tg_yyyyy_yyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 444); 

                auto tg_yyyyy_yyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 445); 

                auto tg_yyyyy_yzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 446); 

                auto tg_yyyyy_zzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 447); 

                auto tg_yyyyz_xxxxxx_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 448); 

                auto tg_yyyyz_xxxxxy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 449); 

                auto tg_yyyyz_xxxxxz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 450); 

                auto tg_yyyyz_xxxxyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 451); 

                auto tg_yyyyz_xxxxyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 452); 

                auto tg_yyyyz_xxxxzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 453); 

                auto tg_yyyyz_xxxyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 454); 

                auto tg_yyyyz_xxxyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 455); 

                auto tg_yyyyz_xxxyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 456); 

                auto tg_yyyyz_xxxzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 457); 

                auto tg_yyyyz_xxyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 458); 

                auto tg_yyyyz_xxyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 459); 

                auto tg_yyyyz_xxyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 460); 

                auto tg_yyyyz_xxyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 461); 

                auto tg_yyyyz_xxzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 462); 

                auto tg_yyyyz_xyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 463); 

                auto tg_yyyyz_xyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 464); 

                auto tg_yyyyz_xyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 465); 

                auto tg_yyyyz_xyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 466); 

                auto tg_yyyyz_xyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 467); 

                auto tg_yyyyz_xzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 468); 

                auto tg_yyyyz_yyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 469); 

                auto tg_yyyyz_yyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 470); 

                auto tg_yyyyz_yyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 471); 

                auto tg_yyyyz_yyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 472); 

                auto tg_yyyyz_yyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 473); 

                auto tg_yyyyz_yzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 474); 

                auto tg_yyyyz_zzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 475); 

                auto tg_yyyzz_xxxxxx_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 476); 

                auto tg_yyyzz_xxxxxy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 477); 

                auto tg_yyyzz_xxxxxz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 478); 

                auto tg_yyyzz_xxxxyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 479); 

                auto tg_yyyzz_xxxxyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 480); 

                auto tg_yyyzz_xxxxzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 481); 

                auto tg_yyyzz_xxxyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 482); 

                auto tg_yyyzz_xxxyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 483); 

                auto tg_yyyzz_xxxyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 484); 

                auto tg_yyyzz_xxxzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 485); 

                auto tg_yyyzz_xxyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 486); 

                auto tg_yyyzz_xxyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 487); 

                auto tg_yyyzz_xxyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 488); 

                auto tg_yyyzz_xxyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 489); 

                auto tg_yyyzz_xxzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 490); 

                auto tg_yyyzz_xyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 491); 

                auto tg_yyyzz_xyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 492); 

                auto tg_yyyzz_xyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 493); 

                auto tg_yyyzz_xyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 494); 

                auto tg_yyyzz_xyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 495); 

                auto tg_yyyzz_xzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 496); 

                auto tg_yyyzz_yyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 497); 

                auto tg_yyyzz_yyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 498); 

                auto tg_yyyzz_yyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 499); 

                auto tg_yyyzz_yyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 500); 

                auto tg_yyyzz_yyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 501); 

                auto tg_yyyzz_yzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 502); 

                auto tg_yyyzz_zzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 503); 

                auto tg_yyzzz_xxxxxx_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 504); 

                auto tg_yyzzz_xxxxxy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 505); 

                auto tg_yyzzz_xxxxxz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 506); 

                auto tg_yyzzz_xxxxyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 507); 

                auto tg_yyzzz_xxxxyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 508); 

                auto tg_yyzzz_xxxxzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 509); 

                auto tg_yyzzz_xxxyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 510); 

                auto tg_yyzzz_xxxyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 511); 

                auto tg_yyzzz_xxxyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 512); 

                auto tg_yyzzz_xxxzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 513); 

                auto tg_yyzzz_xxyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 514); 

                auto tg_yyzzz_xxyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 515); 

                auto tg_yyzzz_xxyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 516); 

                auto tg_yyzzz_xxyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 517); 

                auto tg_yyyy_xxxxxx_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 280); 

                auto tg_yyyy_xxxxxy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 281); 

                auto tg_yyyy_xxxxxz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 282); 

                auto tg_yyyy_xxxxyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 283); 

                auto tg_yyyy_xxxxyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 284); 

                auto tg_yyyy_xxxxzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 285); 

                auto tg_yyyy_xxxyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 286); 

                auto tg_yyyy_xxxyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 287); 

                auto tg_yyyy_xxxyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 288); 

                auto tg_yyyy_xxxzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 289); 

                auto tg_yyyy_xxyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 290); 

                auto tg_yyyy_xxyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 291); 

                auto tg_yyyy_xxyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 292); 

                auto tg_yyyy_xxyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 293); 

                auto tg_yyyy_xxzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 294); 

                auto tg_yyyy_xyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 295); 

                auto tg_yyyy_xyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 296); 

                auto tg_yyyy_xyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 297); 

                auto tg_yyyy_xyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 298); 

                auto tg_yyyy_xyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 299); 

                auto tg_yyyy_xzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 300); 

                auto tg_yyyy_yyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 301); 

                auto tg_yyyy_yyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 302); 

                auto tg_yyyy_yyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 303); 

                auto tg_yyyy_yyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 304); 

                auto tg_yyyy_yyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 305); 

                auto tg_yyyy_yzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 306); 

                auto tg_yyyy_zzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 307); 

                auto tg_yyyz_xxxxxx_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 308); 

                auto tg_yyyz_xxxxxy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 309); 

                auto tg_yyyz_xxxxxz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 310); 

                auto tg_yyyz_xxxxyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 311); 

                auto tg_yyyz_xxxxyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 312); 

                auto tg_yyyz_xxxxzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 313); 

                auto tg_yyyz_xxxyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 314); 

                auto tg_yyyz_xxxyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 315); 

                auto tg_yyyz_xxxyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 316); 

                auto tg_yyyz_xxxzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 317); 

                auto tg_yyyz_xxyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 318); 

                auto tg_yyyz_xxyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 319); 

                auto tg_yyyz_xxyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 320); 

                auto tg_yyyz_xxyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 321); 

                auto tg_yyyz_xxzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 322); 

                auto tg_yyyz_xyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 323); 

                auto tg_yyyz_xyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 324); 

                auto tg_yyyz_xyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 325); 

                auto tg_yyyz_xyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 326); 

                auto tg_yyyz_xyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 327); 

                auto tg_yyyz_xzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 328); 

                auto tg_yyyz_yyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 329); 

                auto tg_yyyz_yyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 330); 

                auto tg_yyyz_yyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 331); 

                auto tg_yyyz_yyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 332); 

                auto tg_yyyz_yyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 333); 

                auto tg_yyyz_yzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 334); 

                auto tg_yyyz_zzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 335); 

                auto tg_yyzz_xxxxxx_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 336); 

                auto tg_yyzz_xxxxxy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 337); 

                auto tg_yyzz_xxxxxz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 338); 

                auto tg_yyzz_xxxxyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 339); 

                auto tg_yyzz_xxxxyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 340); 

                auto tg_yyzz_xxxxzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 341); 

                auto tg_yyzz_xxxyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 342); 

                auto tg_yyzz_xxxyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 343); 

                auto tg_yyzz_xxxyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 344); 

                auto tg_yyzz_xxxzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 345); 

                auto tg_yyzz_xxyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 346); 

                auto tg_yyzz_xxyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 347); 

                auto tg_yyzz_xxyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 348); 

                auto tg_yyzz_xxyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 349); 

                auto tg_yyzz_xxzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 350); 

                auto tg_yyzz_xyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 351); 

                auto tg_yyzz_xyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 352); 

                auto tg_yyzz_xyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 353); 

                auto tg_yyzz_xyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 354); 

                auto tg_yyzz_xyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 355); 

                auto tg_yyzz_xzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 356); 

                auto tg_yyzz_yyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 357); 

                auto tg_yyzz_yyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 358); 

                auto tg_yyzz_yyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 359); 

                auto tg_yyzz_yyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 360); 

                auto tg_yyzz_yyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 361); 

                auto tg_yyzz_yzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 362); 

                auto tg_yyzz_zzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 363); 

                auto tg_yzzz_xxxxxx_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 364); 

                auto tg_yzzz_xxxxxy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 365); 

                auto tg_yzzz_xxxxxz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 366); 

                auto tg_yzzz_xxxxyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 367); 

                auto tg_yzzz_xxxxyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 368); 

                auto tg_yzzz_xxxxzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 369); 

                auto tg_yzzz_xxxyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 370); 

                auto tg_yzzz_xxxyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 371); 

                auto tg_yzzz_xxxyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 372); 

                auto tg_yzzz_xxxzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 373); 

                auto tg_yzzz_xxyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 374); 

                auto tg_yzzz_xxyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 375); 

                auto tg_yzzz_xxyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 376); 

                auto tg_yzzz_xxyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 377); 

                auto tg_yyyy_xxxxxx_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 280); 

                auto tg_yyyy_xxxxxy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 281); 

                auto tg_yyyy_xxxxxz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 282); 

                auto tg_yyyy_xxxxyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 283); 

                auto tg_yyyy_xxxxyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 284); 

                auto tg_yyyy_xxxxzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 285); 

                auto tg_yyyy_xxxyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 286); 

                auto tg_yyyy_xxxyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 287); 

                auto tg_yyyy_xxxyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 288); 

                auto tg_yyyy_xxxzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 289); 

                auto tg_yyyy_xxyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 290); 

                auto tg_yyyy_xxyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 291); 

                auto tg_yyyy_xxyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 292); 

                auto tg_yyyy_xxyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 293); 

                auto tg_yyyy_xxzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 294); 

                auto tg_yyyy_xyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 295); 

                auto tg_yyyy_xyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 296); 

                auto tg_yyyy_xyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 297); 

                auto tg_yyyy_xyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 298); 

                auto tg_yyyy_xyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 299); 

                auto tg_yyyy_xzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 300); 

                auto tg_yyyy_yyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 301); 

                auto tg_yyyy_yyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 302); 

                auto tg_yyyy_yyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 303); 

                auto tg_yyyy_yyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 304); 

                auto tg_yyyy_yyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 305); 

                auto tg_yyyy_yzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 306); 

                auto tg_yyyy_zzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 307); 

                auto tg_yyyz_xxxxxx_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 308); 

                auto tg_yyyz_xxxxxy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 309); 

                auto tg_yyyz_xxxxxz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 310); 

                auto tg_yyyz_xxxxyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 311); 

                auto tg_yyyz_xxxxyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 312); 

                auto tg_yyyz_xxxxzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 313); 

                auto tg_yyyz_xxxyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 314); 

                auto tg_yyyz_xxxyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 315); 

                auto tg_yyyz_xxxyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 316); 

                auto tg_yyyz_xxxzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 317); 

                auto tg_yyyz_xxyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 318); 

                auto tg_yyyz_xxyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 319); 

                auto tg_yyyz_xxyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 320); 

                auto tg_yyyz_xxyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 321); 

                auto tg_yyyz_xxzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 322); 

                auto tg_yyyz_xyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 323); 

                auto tg_yyyz_xyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 324); 

                auto tg_yyyz_xyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 325); 

                auto tg_yyyz_xyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 326); 

                auto tg_yyyz_xyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 327); 

                auto tg_yyyz_xzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 328); 

                auto tg_yyyz_yyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 329); 

                auto tg_yyyz_yyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 330); 

                auto tg_yyyz_yyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 331); 

                auto tg_yyyz_yyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 332); 

                auto tg_yyyz_yyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 333); 

                auto tg_yyyz_yzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 334); 

                auto tg_yyyz_zzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 335); 

                auto tg_yyzz_xxxxxx_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 336); 

                auto tg_yyzz_xxxxxy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 337); 

                auto tg_yyzz_xxxxxz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 338); 

                auto tg_yyzz_xxxxyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 339); 

                auto tg_yyzz_xxxxyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 340); 

                auto tg_yyzz_xxxxzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 341); 

                auto tg_yyzz_xxxyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 342); 

                auto tg_yyzz_xxxyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 343); 

                auto tg_yyzz_xxxyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 344); 

                auto tg_yyzz_xxxzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 345); 

                auto tg_yyzz_xxyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 346); 

                auto tg_yyzz_xxyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 347); 

                auto tg_yyzz_xxyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 348); 

                auto tg_yyzz_xxyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 349); 

                auto tg_yyzz_xxzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 350); 

                auto tg_yyzz_xyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 351); 

                auto tg_yyzz_xyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 352); 

                auto tg_yyzz_xyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 353); 

                auto tg_yyzz_xyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 354); 

                auto tg_yyzz_xyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 355); 

                auto tg_yyzz_xzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 356); 

                auto tg_yyzz_yyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 357); 

                auto tg_yyzz_yyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 358); 

                auto tg_yyzz_yyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 359); 

                auto tg_yyzz_yyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 360); 

                auto tg_yyzz_yyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 361); 

                auto tg_yyzz_yzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 362); 

                auto tg_yyzz_zzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 363); 

                auto tg_yzzz_xxxxxx_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 364); 

                auto tg_yzzz_xxxxxy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 365); 

                auto tg_yzzz_xxxxxz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 366); 

                auto tg_yzzz_xxxxyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 367); 

                auto tg_yzzz_xxxxyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 368); 

                auto tg_yzzz_xxxxzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 369); 

                auto tg_yzzz_xxxyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 370); 

                auto tg_yzzz_xxxyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 371); 

                auto tg_yzzz_xxxyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 372); 

                auto tg_yzzz_xxxzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 373); 

                auto tg_yzzz_xxyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 374); 

                auto tg_yzzz_xxyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 375); 

                auto tg_yzzz_xxyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 376); 

                auto tg_yzzz_xxyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 377); 

                auto tg_yyyyy_xxxxx_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 315); 

                auto tg_yyyyy_xxxxy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 316); 

                auto tg_yyyyy_xxxxz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 317); 

                auto tg_yyyyy_xxxyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 318); 

                auto tg_yyyyy_xxxyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 319); 

                auto tg_yyyyy_xxxzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 320); 

                auto tg_yyyyy_xxyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 321); 

                auto tg_yyyyy_xxyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 322); 

                auto tg_yyyyy_xxyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 323); 

                auto tg_yyyyy_xxzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 324); 

                auto tg_yyyyy_xyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 325); 

                auto tg_yyyyy_xyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 326); 

                auto tg_yyyyy_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 327); 

                auto tg_yyyyy_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 328); 

                auto tg_yyyyy_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 329); 

                auto tg_yyyyy_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 330); 

                auto tg_yyyyy_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 331); 

                auto tg_yyyyy_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 332); 

                auto tg_yyyyy_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 333); 

                auto tg_yyyyy_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 334); 

                auto tg_yyyyy_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 335); 

                auto tg_yyyyz_xxxxx_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 336); 

                auto tg_yyyyz_xxxxy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 337); 

                auto tg_yyyyz_xxxxz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 338); 

                auto tg_yyyyz_xxxyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 339); 

                auto tg_yyyyz_xxxyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 340); 

                auto tg_yyyyz_xxxzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 341); 

                auto tg_yyyyz_xxyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 342); 

                auto tg_yyyyz_xxyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 343); 

                auto tg_yyyyz_xxyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 344); 

                auto tg_yyyyz_xxzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 345); 

                auto tg_yyyyz_xyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 346); 

                auto tg_yyyyz_xyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 347); 

                auto tg_yyyyz_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 348); 

                auto tg_yyyyz_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 349); 

                auto tg_yyyyz_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 350); 

                auto tg_yyyyz_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 351); 

                auto tg_yyyyz_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 352); 

                auto tg_yyyyz_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 353); 

                auto tg_yyyyz_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 354); 

                auto tg_yyyyz_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 355); 

                auto tg_yyyyz_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 356); 

                auto tg_yyyzz_xxxxx_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 357); 

                auto tg_yyyzz_xxxxy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 358); 

                auto tg_yyyzz_xxxxz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 359); 

                auto tg_yyyzz_xxxyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 360); 

                auto tg_yyyzz_xxxyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 361); 

                auto tg_yyyzz_xxxzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 362); 

                auto tg_yyyzz_xxyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 363); 

                auto tg_yyyzz_xxyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 364); 

                auto tg_yyyzz_xxyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 365); 

                auto tg_yyyzz_xxzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 366); 

                auto tg_yyyzz_xyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 367); 

                auto tg_yyyzz_xyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 368); 

                auto tg_yyyzz_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 369); 

                auto tg_yyyzz_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 370); 

                auto tg_yyyzz_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 371); 

                auto tg_yyyzz_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 372); 

                auto tg_yyyzz_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 373); 

                auto tg_yyyzz_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 374); 

                auto tg_yyyzz_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 375); 

                auto tg_yyyzz_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 376); 

                auto tg_yyyzz_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 377); 

                auto tg_yyzzz_xxxxx_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 378); 

                auto tg_yyzzz_xxxxy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 379); 

                auto tg_yyzzz_xxxxz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 380); 

                auto tg_yyzzz_xxxyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 381); 

                auto tg_yyzzz_xxxyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 382); 

                auto tg_yyzzz_xxxzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 383); 

                auto tg_yyzzz_xxyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 384); 

                auto tg_yyzzz_xxyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 385); 

                auto tg_yyzzz_xxyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 386); 

                auto tg_yyzzz_xxzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 387); 

                // set up pointers to integrals

                auto tg_yyyyyy_xxxxxx_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 588); 

                auto tg_yyyyyy_xxxxxy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 589); 

                auto tg_yyyyyy_xxxxxz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 590); 

                auto tg_yyyyyy_xxxxyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 591); 

                auto tg_yyyyyy_xxxxyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 592); 

                auto tg_yyyyyy_xxxxzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 593); 

                auto tg_yyyyyy_xxxyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 594); 

                auto tg_yyyyyy_xxxyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 595); 

                auto tg_yyyyyy_xxxyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 596); 

                auto tg_yyyyyy_xxxzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 597); 

                auto tg_yyyyyy_xxyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 598); 

                auto tg_yyyyyy_xxyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 599); 

                auto tg_yyyyyy_xxyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 600); 

                auto tg_yyyyyy_xxyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 601); 

                auto tg_yyyyyy_xxzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 602); 

                auto tg_yyyyyy_xyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 603); 

                auto tg_yyyyyy_xyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 604); 

                auto tg_yyyyyy_xyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 605); 

                auto tg_yyyyyy_xyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 606); 

                auto tg_yyyyyy_xyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 607); 

                auto tg_yyyyyy_xzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 608); 

                auto tg_yyyyyy_yyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 609); 

                auto tg_yyyyyy_yyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 610); 

                auto tg_yyyyyy_yyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 611); 

                auto tg_yyyyyy_yyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 612); 

                auto tg_yyyyyy_yyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 613); 

                auto tg_yyyyyy_yzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 614); 

                auto tg_yyyyyy_zzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 615); 

                auto tg_yyyyyz_xxxxxx_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 616); 

                auto tg_yyyyyz_xxxxxy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 617); 

                auto tg_yyyyyz_xxxxxz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 618); 

                auto tg_yyyyyz_xxxxyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 619); 

                auto tg_yyyyyz_xxxxyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 620); 

                auto tg_yyyyyz_xxxxzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 621); 

                auto tg_yyyyyz_xxxyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 622); 

                auto tg_yyyyyz_xxxyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 623); 

                auto tg_yyyyyz_xxxyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 624); 

                auto tg_yyyyyz_xxxzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 625); 

                auto tg_yyyyyz_xxyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 626); 

                auto tg_yyyyyz_xxyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 627); 

                auto tg_yyyyyz_xxyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 628); 

                auto tg_yyyyyz_xxyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 629); 

                auto tg_yyyyyz_xxzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 630); 

                auto tg_yyyyyz_xyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 631); 

                auto tg_yyyyyz_xyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 632); 

                auto tg_yyyyyz_xyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 633); 

                auto tg_yyyyyz_xyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 634); 

                auto tg_yyyyyz_xyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 635); 

                auto tg_yyyyyz_xzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 636); 

                auto tg_yyyyyz_yyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 637); 

                auto tg_yyyyyz_yyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 638); 

                auto tg_yyyyyz_yyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 639); 

                auto tg_yyyyyz_yyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 640); 

                auto tg_yyyyyz_yyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 641); 

                auto tg_yyyyyz_yzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 642); 

                auto tg_yyyyyz_zzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 643); 

                auto tg_yyyyzz_xxxxxx_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 644); 

                auto tg_yyyyzz_xxxxxy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 645); 

                auto tg_yyyyzz_xxxxxz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 646); 

                auto tg_yyyyzz_xxxxyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 647); 

                auto tg_yyyyzz_xxxxyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 648); 

                auto tg_yyyyzz_xxxxzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 649); 

                auto tg_yyyyzz_xxxyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 650); 

                auto tg_yyyyzz_xxxyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 651); 

                auto tg_yyyyzz_xxxyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 652); 

                auto tg_yyyyzz_xxxzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 653); 

                auto tg_yyyyzz_xxyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 654); 

                auto tg_yyyyzz_xxyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 655); 

                auto tg_yyyyzz_xxyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 656); 

                auto tg_yyyyzz_xxyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 657); 

                auto tg_yyyyzz_xxzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 658); 

                auto tg_yyyyzz_xyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 659); 

                auto tg_yyyyzz_xyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 660); 

                auto tg_yyyyzz_xyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 661); 

                auto tg_yyyyzz_xyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 662); 

                auto tg_yyyyzz_xyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 663); 

                auto tg_yyyyzz_xzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 664); 

                auto tg_yyyyzz_yyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 665); 

                auto tg_yyyyzz_yyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 666); 

                auto tg_yyyyzz_yyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 667); 

                auto tg_yyyyzz_yyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 668); 

                auto tg_yyyyzz_yyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 669); 

                auto tg_yyyyzz_yzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 670); 

                auto tg_yyyyzz_zzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 671); 

                auto tg_yyyzzz_xxxxxx_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 672); 

                auto tg_yyyzzz_xxxxxy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 673); 

                auto tg_yyyzzz_xxxxxz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 674); 

                auto tg_yyyzzz_xxxxyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 675); 

                auto tg_yyyzzz_xxxxyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 676); 

                auto tg_yyyzzz_xxxxzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 677); 

                auto tg_yyyzzz_xxxyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 678); 

                auto tg_yyyzzz_xxxyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 679); 

                auto tg_yyyzzz_xxxyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 680); 

                auto tg_yyyzzz_xxxzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 681); 

                auto tg_yyyzzz_xxyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 682); 

                auto tg_yyyzzz_xxyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 683); 

                auto tg_yyyzzz_xxyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 684); 

                auto tg_yyyzzz_xxyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 685); 

                // Batch of Integrals (588,686)

                #pragma omp simd aligned(fxn, fza, tg_yyyy_xxxxxx_0, tg_yyyy_xxxxxx_1, tg_yyyy_xxxxxy_0, \
                                         tg_yyyy_xxxxxy_1, tg_yyyy_xxxxxz_0, tg_yyyy_xxxxxz_1, tg_yyyy_xxxxyy_0, \
                                         tg_yyyy_xxxxyy_1, tg_yyyy_xxxxyz_0, tg_yyyy_xxxxyz_1, tg_yyyy_xxxxzz_0, \
                                         tg_yyyy_xxxxzz_1, tg_yyyy_xxxyyy_0, tg_yyyy_xxxyyy_1, tg_yyyy_xxxyyz_0, \
                                         tg_yyyy_xxxyyz_1, tg_yyyy_xxxyzz_0, tg_yyyy_xxxyzz_1, tg_yyyy_xxxzzz_0, \
                                         tg_yyyy_xxxzzz_1, tg_yyyy_xxyyyy_0, tg_yyyy_xxyyyy_1, tg_yyyy_xxyyyz_0, \
                                         tg_yyyy_xxyyyz_1, tg_yyyy_xxyyzz_0, tg_yyyy_xxyyzz_1, tg_yyyy_xxyzzz_0, \
                                         tg_yyyy_xxyzzz_1, tg_yyyy_xxzzzz_0, tg_yyyy_xxzzzz_1, tg_yyyy_xyyyyy_0, \
                                         tg_yyyy_xyyyyy_1, tg_yyyy_xyyyyz_0, tg_yyyy_xyyyyz_1, tg_yyyy_xyyyzz_0, \
                                         tg_yyyy_xyyyzz_1, tg_yyyy_xyyzzz_0, tg_yyyy_xyyzzz_1, tg_yyyy_xyzzzz_0, \
                                         tg_yyyy_xyzzzz_1, tg_yyyy_xzzzzz_0, tg_yyyy_xzzzzz_1, tg_yyyy_yyyyyy_0, \
                                         tg_yyyy_yyyyyy_1, tg_yyyy_yyyyyz_0, tg_yyyy_yyyyyz_1, tg_yyyy_yyyyzz_0, \
                                         tg_yyyy_yyyyzz_1, tg_yyyy_yyyzzz_0, tg_yyyy_yyyzzz_1, tg_yyyy_yyzzzz_0, \
                                         tg_yyyy_yyzzzz_1, tg_yyyy_yzzzzz_0, tg_yyyy_yzzzzz_1, tg_yyyy_zzzzzz_0, \
                                         tg_yyyy_zzzzzz_1, tg_yyyyy_xxxxx_1, tg_yyyyy_xxxxxx_0, tg_yyyyy_xxxxxx_1, \
                                         tg_yyyyy_xxxxxy_0, tg_yyyyy_xxxxxy_1, tg_yyyyy_xxxxxz_0, tg_yyyyy_xxxxxz_1, \
                                         tg_yyyyy_xxxxy_1, tg_yyyyy_xxxxyy_0, tg_yyyyy_xxxxyy_1, tg_yyyyy_xxxxyz_0, \
                                         tg_yyyyy_xxxxyz_1, tg_yyyyy_xxxxz_1, tg_yyyyy_xxxxzz_0, tg_yyyyy_xxxxzz_1, \
                                         tg_yyyyy_xxxyy_1, tg_yyyyy_xxxyyy_0, tg_yyyyy_xxxyyy_1, tg_yyyyy_xxxyyz_0, \
                                         tg_yyyyy_xxxyyz_1, tg_yyyyy_xxxyz_1, tg_yyyyy_xxxyzz_0, tg_yyyyy_xxxyzz_1, \
                                         tg_yyyyy_xxxzz_1, tg_yyyyy_xxxzzz_0, tg_yyyyy_xxxzzz_1, tg_yyyyy_xxyyy_1, \
                                         tg_yyyyy_xxyyyy_0, tg_yyyyy_xxyyyy_1, tg_yyyyy_xxyyyz_0, tg_yyyyy_xxyyyz_1, \
                                         tg_yyyyy_xxyyz_1, tg_yyyyy_xxyyzz_0, tg_yyyyy_xxyyzz_1, tg_yyyyy_xxyzz_1, \
                                         tg_yyyyy_xxyzzz_0, tg_yyyyy_xxyzzz_1, tg_yyyyy_xxzzz_1, tg_yyyyy_xxzzzz_0, \
                                         tg_yyyyy_xxzzzz_1, tg_yyyyy_xyyyy_1, tg_yyyyy_xyyyyy_0, tg_yyyyy_xyyyyy_1, \
                                         tg_yyyyy_xyyyyz_0, tg_yyyyy_xyyyyz_1, tg_yyyyy_xyyyz_1, tg_yyyyy_xyyyzz_0, \
                                         tg_yyyyy_xyyyzz_1, tg_yyyyy_xyyzz_1, tg_yyyyy_xyyzzz_0, tg_yyyyy_xyyzzz_1, \
                                         tg_yyyyy_xyzzz_1, tg_yyyyy_xyzzzz_0, tg_yyyyy_xyzzzz_1, tg_yyyyy_xzzzz_1, \
                                         tg_yyyyy_xzzzzz_0, tg_yyyyy_xzzzzz_1, tg_yyyyy_yyyyy_1, tg_yyyyy_yyyyyy_0, \
                                         tg_yyyyy_yyyyyy_1, tg_yyyyy_yyyyyz_0, tg_yyyyy_yyyyyz_1, tg_yyyyy_yyyyz_1, \
                                         tg_yyyyy_yyyyzz_0, tg_yyyyy_yyyyzz_1, tg_yyyyy_yyyzz_1, tg_yyyyy_yyyzzz_0, \
                                         tg_yyyyy_yyyzzz_1, tg_yyyyy_yyzzz_1, tg_yyyyy_yyzzzz_0, tg_yyyyy_yyzzzz_1, \
                                         tg_yyyyy_yzzzz_1, tg_yyyyy_yzzzzz_0, tg_yyyyy_yzzzzz_1, tg_yyyyy_zzzzz_1, \
                                         tg_yyyyy_zzzzzz_0, tg_yyyyy_zzzzzz_1, tg_yyyyyy_xxxxxx_0, tg_yyyyyy_xxxxxy_0, \
                                         tg_yyyyyy_xxxxxz_0, tg_yyyyyy_xxxxyy_0, tg_yyyyyy_xxxxyz_0, tg_yyyyyy_xxxxzz_0, \
                                         tg_yyyyyy_xxxyyy_0, tg_yyyyyy_xxxyyz_0, tg_yyyyyy_xxxyzz_0, tg_yyyyyy_xxxzzz_0, \
                                         tg_yyyyyy_xxyyyy_0, tg_yyyyyy_xxyyyz_0, tg_yyyyyy_xxyyzz_0, tg_yyyyyy_xxyzzz_0, \
                                         tg_yyyyyy_xxzzzz_0, tg_yyyyyy_xyyyyy_0, tg_yyyyyy_xyyyyz_0, tg_yyyyyy_xyyyzz_0, \
                                         tg_yyyyyy_xyyzzz_0, tg_yyyyyy_xyzzzz_0, tg_yyyyyy_xzzzzz_0, tg_yyyyyy_yyyyyy_0, \
                                         tg_yyyyyy_yyyyyz_0, tg_yyyyyy_yyyyzz_0, tg_yyyyyy_yyyzzz_0, tg_yyyyyy_yyzzzz_0, \
                                         tg_yyyyyy_yzzzzz_0, tg_yyyyyy_zzzzzz_0, tg_yyyyyz_xxxxxx_0, tg_yyyyyz_xxxxxy_0, \
                                         tg_yyyyyz_xxxxxz_0, tg_yyyyyz_xxxxyy_0, tg_yyyyyz_xxxxyz_0, tg_yyyyyz_xxxxzz_0, \
                                         tg_yyyyyz_xxxyyy_0, tg_yyyyyz_xxxyyz_0, tg_yyyyyz_xxxyzz_0, tg_yyyyyz_xxxzzz_0, \
                                         tg_yyyyyz_xxyyyy_0, tg_yyyyyz_xxyyyz_0, tg_yyyyyz_xxyyzz_0, tg_yyyyyz_xxyzzz_0, \
                                         tg_yyyyyz_xxzzzz_0, tg_yyyyyz_xyyyyy_0, tg_yyyyyz_xyyyyz_0, tg_yyyyyz_xyyyzz_0, \
                                         tg_yyyyyz_xyyzzz_0, tg_yyyyyz_xyzzzz_0, tg_yyyyyz_xzzzzz_0, tg_yyyyyz_yyyyyy_0, \
                                         tg_yyyyyz_yyyyyz_0, tg_yyyyyz_yyyyzz_0, tg_yyyyyz_yyyzzz_0, tg_yyyyyz_yyzzzz_0, \
                                         tg_yyyyyz_yzzzzz_0, tg_yyyyyz_zzzzzz_0, tg_yyyyz_xxxxx_1, tg_yyyyz_xxxxxx_0, \
                                         tg_yyyyz_xxxxxx_1, tg_yyyyz_xxxxxy_0, tg_yyyyz_xxxxxy_1, tg_yyyyz_xxxxxz_0, \
                                         tg_yyyyz_xxxxxz_1, tg_yyyyz_xxxxy_1, tg_yyyyz_xxxxyy_0, tg_yyyyz_xxxxyy_1, \
                                         tg_yyyyz_xxxxyz_0, tg_yyyyz_xxxxyz_1, tg_yyyyz_xxxxz_1, tg_yyyyz_xxxxzz_0, \
                                         tg_yyyyz_xxxxzz_1, tg_yyyyz_xxxyy_1, tg_yyyyz_xxxyyy_0, tg_yyyyz_xxxyyy_1, \
                                         tg_yyyyz_xxxyyz_0, tg_yyyyz_xxxyyz_1, tg_yyyyz_xxxyz_1, tg_yyyyz_xxxyzz_0, \
                                         tg_yyyyz_xxxyzz_1, tg_yyyyz_xxxzz_1, tg_yyyyz_xxxzzz_0, tg_yyyyz_xxxzzz_1, \
                                         tg_yyyyz_xxyyy_1, tg_yyyyz_xxyyyy_0, tg_yyyyz_xxyyyy_1, tg_yyyyz_xxyyyz_0, \
                                         tg_yyyyz_xxyyyz_1, tg_yyyyz_xxyyz_1, tg_yyyyz_xxyyzz_0, tg_yyyyz_xxyyzz_1, \
                                         tg_yyyyz_xxyzz_1, tg_yyyyz_xxyzzz_0, tg_yyyyz_xxyzzz_1, tg_yyyyz_xxzzz_1, \
                                         tg_yyyyz_xxzzzz_0, tg_yyyyz_xxzzzz_1, tg_yyyyz_xyyyy_1, tg_yyyyz_xyyyyy_0, \
                                         tg_yyyyz_xyyyyy_1, tg_yyyyz_xyyyyz_0, tg_yyyyz_xyyyyz_1, tg_yyyyz_xyyyz_1, \
                                         tg_yyyyz_xyyyzz_0, tg_yyyyz_xyyyzz_1, tg_yyyyz_xyyzz_1, tg_yyyyz_xyyzzz_0, \
                                         tg_yyyyz_xyyzzz_1, tg_yyyyz_xyzzz_1, tg_yyyyz_xyzzzz_0, tg_yyyyz_xyzzzz_1, \
                                         tg_yyyyz_xzzzz_1, tg_yyyyz_xzzzzz_0, tg_yyyyz_xzzzzz_1, tg_yyyyz_yyyyy_1, \
                                         tg_yyyyz_yyyyyy_0, tg_yyyyz_yyyyyy_1, tg_yyyyz_yyyyyz_0, tg_yyyyz_yyyyyz_1, \
                                         tg_yyyyz_yyyyz_1, tg_yyyyz_yyyyzz_0, tg_yyyyz_yyyyzz_1, tg_yyyyz_yyyzz_1, \
                                         tg_yyyyz_yyyzzz_0, tg_yyyyz_yyyzzz_1, tg_yyyyz_yyzzz_1, tg_yyyyz_yyzzzz_0, \
                                         tg_yyyyz_yyzzzz_1, tg_yyyyz_yzzzz_1, tg_yyyyz_yzzzzz_0, tg_yyyyz_yzzzzz_1, \
                                         tg_yyyyz_zzzzz_1, tg_yyyyz_zzzzzz_0, tg_yyyyz_zzzzzz_1, tg_yyyyzz_xxxxxx_0, \
                                         tg_yyyyzz_xxxxxy_0, tg_yyyyzz_xxxxxz_0, tg_yyyyzz_xxxxyy_0, tg_yyyyzz_xxxxyz_0, \
                                         tg_yyyyzz_xxxxzz_0, tg_yyyyzz_xxxyyy_0, tg_yyyyzz_xxxyyz_0, tg_yyyyzz_xxxyzz_0, \
                                         tg_yyyyzz_xxxzzz_0, tg_yyyyzz_xxyyyy_0, tg_yyyyzz_xxyyyz_0, tg_yyyyzz_xxyyzz_0, \
                                         tg_yyyyzz_xxyzzz_0, tg_yyyyzz_xxzzzz_0, tg_yyyyzz_xyyyyy_0, tg_yyyyzz_xyyyyz_0, \
                                         tg_yyyyzz_xyyyzz_0, tg_yyyyzz_xyyzzz_0, tg_yyyyzz_xyzzzz_0, tg_yyyyzz_xzzzzz_0, \
                                         tg_yyyyzz_yyyyyy_0, tg_yyyyzz_yyyyyz_0, tg_yyyyzz_yyyyzz_0, tg_yyyyzz_yyyzzz_0, \
                                         tg_yyyyzz_yyzzzz_0, tg_yyyyzz_yzzzzz_0, tg_yyyyzz_zzzzzz_0, tg_yyyz_xxxxxx_0, \
                                         tg_yyyz_xxxxxx_1, tg_yyyz_xxxxxy_0, tg_yyyz_xxxxxy_1, tg_yyyz_xxxxxz_0, \
                                         tg_yyyz_xxxxxz_1, tg_yyyz_xxxxyy_0, tg_yyyz_xxxxyy_1, tg_yyyz_xxxxyz_0, \
                                         tg_yyyz_xxxxyz_1, tg_yyyz_xxxxzz_0, tg_yyyz_xxxxzz_1, tg_yyyz_xxxyyy_0, \
                                         tg_yyyz_xxxyyy_1, tg_yyyz_xxxyyz_0, tg_yyyz_xxxyyz_1, tg_yyyz_xxxyzz_0, \
                                         tg_yyyz_xxxyzz_1, tg_yyyz_xxxzzz_0, tg_yyyz_xxxzzz_1, tg_yyyz_xxyyyy_0, \
                                         tg_yyyz_xxyyyy_1, tg_yyyz_xxyyyz_0, tg_yyyz_xxyyyz_1, tg_yyyz_xxyyzz_0, \
                                         tg_yyyz_xxyyzz_1, tg_yyyz_xxyzzz_0, tg_yyyz_xxyzzz_1, tg_yyyz_xxzzzz_0, \
                                         tg_yyyz_xxzzzz_1, tg_yyyz_xyyyyy_0, tg_yyyz_xyyyyy_1, tg_yyyz_xyyyyz_0, \
                                         tg_yyyz_xyyyyz_1, tg_yyyz_xyyyzz_0, tg_yyyz_xyyyzz_1, tg_yyyz_xyyzzz_0, \
                                         tg_yyyz_xyyzzz_1, tg_yyyz_xyzzzz_0, tg_yyyz_xyzzzz_1, tg_yyyz_xzzzzz_0, \
                                         tg_yyyz_xzzzzz_1, tg_yyyz_yyyyyy_0, tg_yyyz_yyyyyy_1, tg_yyyz_yyyyyz_0, \
                                         tg_yyyz_yyyyyz_1, tg_yyyz_yyyyzz_0, tg_yyyz_yyyyzz_1, tg_yyyz_yyyzzz_0, \
                                         tg_yyyz_yyyzzz_1, tg_yyyz_yyzzzz_0, tg_yyyz_yyzzzz_1, tg_yyyz_yzzzzz_0, \
                                         tg_yyyz_yzzzzz_1, tg_yyyz_zzzzzz_0, tg_yyyz_zzzzzz_1, tg_yyyzz_xxxxx_1, \
                                         tg_yyyzz_xxxxxx_0, tg_yyyzz_xxxxxx_1, tg_yyyzz_xxxxxy_0, tg_yyyzz_xxxxxy_1, \
                                         tg_yyyzz_xxxxxz_0, tg_yyyzz_xxxxxz_1, tg_yyyzz_xxxxy_1, tg_yyyzz_xxxxyy_0, \
                                         tg_yyyzz_xxxxyy_1, tg_yyyzz_xxxxyz_0, tg_yyyzz_xxxxyz_1, tg_yyyzz_xxxxz_1, \
                                         tg_yyyzz_xxxxzz_0, tg_yyyzz_xxxxzz_1, tg_yyyzz_xxxyy_1, tg_yyyzz_xxxyyy_0, \
                                         tg_yyyzz_xxxyyy_1, tg_yyyzz_xxxyyz_0, tg_yyyzz_xxxyyz_1, tg_yyyzz_xxxyz_1, \
                                         tg_yyyzz_xxxyzz_0, tg_yyyzz_xxxyzz_1, tg_yyyzz_xxxzz_1, tg_yyyzz_xxxzzz_0, \
                                         tg_yyyzz_xxxzzz_1, tg_yyyzz_xxyyy_1, tg_yyyzz_xxyyyy_0, tg_yyyzz_xxyyyy_1, \
                                         tg_yyyzz_xxyyyz_0, tg_yyyzz_xxyyyz_1, tg_yyyzz_xxyyz_1, tg_yyyzz_xxyyzz_0, \
                                         tg_yyyzz_xxyyzz_1, tg_yyyzz_xxyzz_1, tg_yyyzz_xxyzzz_0, tg_yyyzz_xxyzzz_1, \
                                         tg_yyyzz_xxzzz_1, tg_yyyzz_xxzzzz_0, tg_yyyzz_xxzzzz_1, tg_yyyzz_xyyyy_1, \
                                         tg_yyyzz_xyyyyy_0, tg_yyyzz_xyyyyy_1, tg_yyyzz_xyyyyz_0, tg_yyyzz_xyyyyz_1, \
                                         tg_yyyzz_xyyyz_1, tg_yyyzz_xyyyzz_0, tg_yyyzz_xyyyzz_1, tg_yyyzz_xyyzz_1, \
                                         tg_yyyzz_xyyzzz_0, tg_yyyzz_xyyzzz_1, tg_yyyzz_xyzzz_1, tg_yyyzz_xyzzzz_0, \
                                         tg_yyyzz_xyzzzz_1, tg_yyyzz_xzzzz_1, tg_yyyzz_xzzzzz_0, tg_yyyzz_xzzzzz_1, \
                                         tg_yyyzz_yyyyy_1, tg_yyyzz_yyyyyy_0, tg_yyyzz_yyyyyy_1, tg_yyyzz_yyyyyz_0, \
                                         tg_yyyzz_yyyyyz_1, tg_yyyzz_yyyyz_1, tg_yyyzz_yyyyzz_0, tg_yyyzz_yyyyzz_1, \
                                         tg_yyyzz_yyyzz_1, tg_yyyzz_yyyzzz_0, tg_yyyzz_yyyzzz_1, tg_yyyzz_yyzzz_1, \
                                         tg_yyyzz_yyzzzz_0, tg_yyyzz_yyzzzz_1, tg_yyyzz_yzzzz_1, tg_yyyzz_yzzzzz_0, \
                                         tg_yyyzz_yzzzzz_1, tg_yyyzz_zzzzz_1, tg_yyyzz_zzzzzz_0, tg_yyyzz_zzzzzz_1, \
                                         tg_yyyzzz_xxxxxx_0, tg_yyyzzz_xxxxxy_0, tg_yyyzzz_xxxxxz_0, tg_yyyzzz_xxxxyy_0, \
                                         tg_yyyzzz_xxxxyz_0, tg_yyyzzz_xxxxzz_0, tg_yyyzzz_xxxyyy_0, tg_yyyzzz_xxxyyz_0, \
                                         tg_yyyzzz_xxxyzz_0, tg_yyyzzz_xxxzzz_0, tg_yyyzzz_xxyyyy_0, tg_yyyzzz_xxyyyz_0, \
                                         tg_yyyzzz_xxyyzz_0, tg_yyyzzz_xxyzzz_0, tg_yyzz_xxxxxx_0, tg_yyzz_xxxxxx_1, \
                                         tg_yyzz_xxxxxy_0, tg_yyzz_xxxxxy_1, tg_yyzz_xxxxxz_0, tg_yyzz_xxxxxz_1, \
                                         tg_yyzz_xxxxyy_0, tg_yyzz_xxxxyy_1, tg_yyzz_xxxxyz_0, tg_yyzz_xxxxyz_1, \
                                         tg_yyzz_xxxxzz_0, tg_yyzz_xxxxzz_1, tg_yyzz_xxxyyy_0, tg_yyzz_xxxyyy_1, \
                                         tg_yyzz_xxxyyz_0, tg_yyzz_xxxyyz_1, tg_yyzz_xxxyzz_0, tg_yyzz_xxxyzz_1, \
                                         tg_yyzz_xxxzzz_0, tg_yyzz_xxxzzz_1, tg_yyzz_xxyyyy_0, tg_yyzz_xxyyyy_1, \
                                         tg_yyzz_xxyyyz_0, tg_yyzz_xxyyyz_1, tg_yyzz_xxyyzz_0, tg_yyzz_xxyyzz_1, \
                                         tg_yyzz_xxyzzz_0, tg_yyzz_xxyzzz_1, tg_yyzz_xxzzzz_0, tg_yyzz_xxzzzz_1, \
                                         tg_yyzz_xyyyyy_0, tg_yyzz_xyyyyy_1, tg_yyzz_xyyyyz_0, tg_yyzz_xyyyyz_1, \
                                         tg_yyzz_xyyyzz_0, tg_yyzz_xyyyzz_1, tg_yyzz_xyyzzz_0, tg_yyzz_xyyzzz_1, \
                                         tg_yyzz_xyzzzz_0, tg_yyzz_xyzzzz_1, tg_yyzz_xzzzzz_0, tg_yyzz_xzzzzz_1, \
                                         tg_yyzz_yyyyyy_0, tg_yyzz_yyyyyy_1, tg_yyzz_yyyyyz_0, tg_yyzz_yyyyyz_1, \
                                         tg_yyzz_yyyyzz_0, tg_yyzz_yyyyzz_1, tg_yyzz_yyyzzz_0, tg_yyzz_yyyzzz_1, \
                                         tg_yyzz_yyzzzz_0, tg_yyzz_yyzzzz_1, tg_yyzz_yzzzzz_0, tg_yyzz_yzzzzz_1, \
                                         tg_yyzz_zzzzzz_0, tg_yyzz_zzzzzz_1, tg_yyzzz_xxxxx_1, tg_yyzzz_xxxxxx_0, \
                                         tg_yyzzz_xxxxxx_1, tg_yyzzz_xxxxxy_0, tg_yyzzz_xxxxxy_1, tg_yyzzz_xxxxxz_0, \
                                         tg_yyzzz_xxxxxz_1, tg_yyzzz_xxxxy_1, tg_yyzzz_xxxxyy_0, tg_yyzzz_xxxxyy_1, \
                                         tg_yyzzz_xxxxyz_0, tg_yyzzz_xxxxyz_1, tg_yyzzz_xxxxz_1, tg_yyzzz_xxxxzz_0, \
                                         tg_yyzzz_xxxxzz_1, tg_yyzzz_xxxyy_1, tg_yyzzz_xxxyyy_0, tg_yyzzz_xxxyyy_1, \
                                         tg_yyzzz_xxxyyz_0, tg_yyzzz_xxxyyz_1, tg_yyzzz_xxxyz_1, tg_yyzzz_xxxyzz_0, \
                                         tg_yyzzz_xxxyzz_1, tg_yyzzz_xxxzz_1, tg_yyzzz_xxxzzz_0, tg_yyzzz_xxxzzz_1, \
                                         tg_yyzzz_xxyyy_1, tg_yyzzz_xxyyyy_0, tg_yyzzz_xxyyyy_1, tg_yyzzz_xxyyyz_0, \
                                         tg_yyzzz_xxyyyz_1, tg_yyzzz_xxyyz_1, tg_yyzzz_xxyyzz_0, tg_yyzzz_xxyyzz_1, \
                                         tg_yyzzz_xxyzz_1, tg_yyzzz_xxyzzz_0, tg_yyzzz_xxyzzz_1, tg_yyzzz_xxzzz_1, \
                                         tg_yzzz_xxxxxx_0, tg_yzzz_xxxxxx_1, tg_yzzz_xxxxxy_0, tg_yzzz_xxxxxy_1, \
                                         tg_yzzz_xxxxxz_0, tg_yzzz_xxxxxz_1, tg_yzzz_xxxxyy_0, tg_yzzz_xxxxyy_1, \
                                         tg_yzzz_xxxxyz_0, tg_yzzz_xxxxyz_1, tg_yzzz_xxxxzz_0, tg_yzzz_xxxxzz_1, \
                                         tg_yzzz_xxxyyy_0, tg_yzzz_xxxyyy_1, tg_yzzz_xxxyyz_0, tg_yzzz_xxxyyz_1, \
                                         tg_yzzz_xxxyzz_0, tg_yzzz_xxxyzz_1, tg_yzzz_xxxzzz_0, tg_yzzz_xxxzzz_1, \
                                         tg_yzzz_xxyyyy_0, tg_yzzz_xxyyyy_1, tg_yzzz_xxyyyz_0, tg_yzzz_xxyyyz_1, \
                                         tg_yzzz_xxyyzz_0, tg_yzzz_xxyyzz_1, tg_yzzz_xxyzzz_0, tg_yzzz_xxyzzz_1, wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_y[j]; 

                    tg_yyyyyy_xxxxxx_0[j] = pb_y * tg_yyyyy_xxxxxx_0[j] + fr * tg_yyyyy_xxxxxx_1[j] + 2.5 * fl1_fx * (tg_yyyy_xxxxxx_0[j] - tg_yyyy_xxxxxx_1[j] * fl1_fza);

                    tg_yyyyyy_xxxxxy_0[j] = pb_y * tg_yyyyy_xxxxxy_0[j] + fr * tg_yyyyy_xxxxxy_1[j] + 2.5 * fl1_fx * (tg_yyyy_xxxxxy_0[j] - tg_yyyy_xxxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyy_xxxxx_1[j];

                    tg_yyyyyy_xxxxxz_0[j] = pb_y * tg_yyyyy_xxxxxz_0[j] + fr * tg_yyyyy_xxxxxz_1[j] + 2.5 * fl1_fx * (tg_yyyy_xxxxxz_0[j] - tg_yyyy_xxxxxz_1[j] * fl1_fza);

                    tg_yyyyyy_xxxxyy_0[j] = pb_y * tg_yyyyy_xxxxyy_0[j] + fr * tg_yyyyy_xxxxyy_1[j] + 2.5 * fl1_fx * (tg_yyyy_xxxxyy_0[j] - tg_yyyy_xxxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyyy_xxxxy_1[j];

                    tg_yyyyyy_xxxxyz_0[j] = pb_y * tg_yyyyy_xxxxyz_0[j] + fr * tg_yyyyy_xxxxyz_1[j] + 2.5 * fl1_fx * (tg_yyyy_xxxxyz_0[j] - tg_yyyy_xxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyy_xxxxz_1[j];

                    tg_yyyyyy_xxxxzz_0[j] = pb_y * tg_yyyyy_xxxxzz_0[j] + fr * tg_yyyyy_xxxxzz_1[j] + 2.5 * fl1_fx * (tg_yyyy_xxxxzz_0[j] - tg_yyyy_xxxxzz_1[j] * fl1_fza);

                    tg_yyyyyy_xxxyyy_0[j] = pb_y * tg_yyyyy_xxxyyy_0[j] + fr * tg_yyyyy_xxxyyy_1[j] + 2.5 * fl1_fx * (tg_yyyy_xxxyyy_0[j] - tg_yyyy_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyy_xxxyy_1[j];

                    tg_yyyyyy_xxxyyz_0[j] = pb_y * tg_yyyyy_xxxyyz_0[j] + fr * tg_yyyyy_xxxyyz_1[j] + 2.5 * fl1_fx * (tg_yyyy_xxxyyz_0[j] - tg_yyyy_xxxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyy_xxxyz_1[j];

                    tg_yyyyyy_xxxyzz_0[j] = pb_y * tg_yyyyy_xxxyzz_0[j] + fr * tg_yyyyy_xxxyzz_1[j] + 2.5 * fl1_fx * (tg_yyyy_xxxyzz_0[j] - tg_yyyy_xxxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyy_xxxzz_1[j];

                    tg_yyyyyy_xxxzzz_0[j] = pb_y * tg_yyyyy_xxxzzz_0[j] + fr * tg_yyyyy_xxxzzz_1[j] + 2.5 * fl1_fx * (tg_yyyy_xxxzzz_0[j] - tg_yyyy_xxxzzz_1[j] * fl1_fza);

                    tg_yyyyyy_xxyyyy_0[j] = pb_y * tg_yyyyy_xxyyyy_0[j] + fr * tg_yyyyy_xxyyyy_1[j] + 2.5 * fl1_fx * (tg_yyyy_xxyyyy_0[j] - tg_yyyy_xxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyy_xxyyy_1[j];

                    tg_yyyyyy_xxyyyz_0[j] = pb_y * tg_yyyyy_xxyyyz_0[j] + fr * tg_yyyyy_xxyyyz_1[j] + 2.5 * fl1_fx * (tg_yyyy_xxyyyz_0[j] - tg_yyyy_xxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyy_xxyyz_1[j];

                    tg_yyyyyy_xxyyzz_0[j] = pb_y * tg_yyyyy_xxyyzz_0[j] + fr * tg_yyyyy_xxyyzz_1[j] + 2.5 * fl1_fx * (tg_yyyy_xxyyzz_0[j] - tg_yyyy_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyy_xxyzz_1[j];

                    tg_yyyyyy_xxyzzz_0[j] = pb_y * tg_yyyyy_xxyzzz_0[j] + fr * tg_yyyyy_xxyzzz_1[j] + 2.5 * fl1_fx * (tg_yyyy_xxyzzz_0[j] - tg_yyyy_xxyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyy_xxzzz_1[j];

                    tg_yyyyyy_xxzzzz_0[j] = pb_y * tg_yyyyy_xxzzzz_0[j] + fr * tg_yyyyy_xxzzzz_1[j] + 2.5 * fl1_fx * (tg_yyyy_xxzzzz_0[j] - tg_yyyy_xxzzzz_1[j] * fl1_fza);

                    tg_yyyyyy_xyyyyy_0[j] = pb_y * tg_yyyyy_xyyyyy_0[j] + fr * tg_yyyyy_xyyyyy_1[j] + 2.5 * fl1_fx * (tg_yyyy_xyyyyy_0[j] - tg_yyyy_xyyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyyy_xyyyy_1[j];

                    tg_yyyyyy_xyyyyz_0[j] = pb_y * tg_yyyyy_xyyyyz_0[j] + fr * tg_yyyyy_xyyyyz_1[j] + 2.5 * fl1_fx * (tg_yyyy_xyyyyz_0[j] - tg_yyyy_xyyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyy_xyyyz_1[j];

                    tg_yyyyyy_xyyyzz_0[j] = pb_y * tg_yyyyy_xyyyzz_0[j] + fr * tg_yyyyy_xyyyzz_1[j] + 2.5 * fl1_fx * (tg_yyyy_xyyyzz_0[j] - tg_yyyy_xyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyy_xyyzz_1[j];

                    tg_yyyyyy_xyyzzz_0[j] = pb_y * tg_yyyyy_xyyzzz_0[j] + fr * tg_yyyyy_xyyzzz_1[j] + 2.5 * fl1_fx * (tg_yyyy_xyyzzz_0[j] - tg_yyyy_xyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyy_xyzzz_1[j];

                    tg_yyyyyy_xyzzzz_0[j] = pb_y * tg_yyyyy_xyzzzz_0[j] + fr * tg_yyyyy_xyzzzz_1[j] + 2.5 * fl1_fx * (tg_yyyy_xyzzzz_0[j] - tg_yyyy_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyy_xzzzz_1[j];

                    tg_yyyyyy_xzzzzz_0[j] = pb_y * tg_yyyyy_xzzzzz_0[j] + fr * tg_yyyyy_xzzzzz_1[j] + 2.5 * fl1_fx * (tg_yyyy_xzzzzz_0[j] - tg_yyyy_xzzzzz_1[j] * fl1_fza);

                    tg_yyyyyy_yyyyyy_0[j] = pb_y * tg_yyyyy_yyyyyy_0[j] + fr * tg_yyyyy_yyyyyy_1[j] + 2.5 * fl1_fx * (tg_yyyy_yyyyyy_0[j] - tg_yyyy_yyyyyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yyyyy_yyyyy_1[j];

                    tg_yyyyyy_yyyyyz_0[j] = pb_y * tg_yyyyy_yyyyyz_0[j] + fr * tg_yyyyy_yyyyyz_1[j] + 2.5 * fl1_fx * (tg_yyyy_yyyyyz_0[j] - tg_yyyy_yyyyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyyy_yyyyz_1[j];

                    tg_yyyyyy_yyyyzz_0[j] = pb_y * tg_yyyyy_yyyyzz_0[j] + fr * tg_yyyyy_yyyyzz_1[j] + 2.5 * fl1_fx * (tg_yyyy_yyyyzz_0[j] - tg_yyyy_yyyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyy_yyyzz_1[j];

                    tg_yyyyyy_yyyzzz_0[j] = pb_y * tg_yyyyy_yyyzzz_0[j] + fr * tg_yyyyy_yyyzzz_1[j] + 2.5 * fl1_fx * (tg_yyyy_yyyzzz_0[j] - tg_yyyy_yyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyy_yyzzz_1[j];

                    tg_yyyyyy_yyzzzz_0[j] = pb_y * tg_yyyyy_yyzzzz_0[j] + fr * tg_yyyyy_yyzzzz_1[j] + 2.5 * fl1_fx * (tg_yyyy_yyzzzz_0[j] - tg_yyyy_yyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyy_yzzzz_1[j];

                    tg_yyyyyy_yzzzzz_0[j] = pb_y * tg_yyyyy_yzzzzz_0[j] + fr * tg_yyyyy_yzzzzz_1[j] + 2.5 * fl1_fx * (tg_yyyy_yzzzzz_0[j] - tg_yyyy_yzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyy_zzzzz_1[j];

                    tg_yyyyyy_zzzzzz_0[j] = pb_y * tg_yyyyy_zzzzzz_0[j] + fr * tg_yyyyy_zzzzzz_1[j] + 2.5 * fl1_fx * (tg_yyyy_zzzzzz_0[j] - tg_yyyy_zzzzzz_1[j] * fl1_fza);

                    tg_yyyyyz_xxxxxx_0[j] = pb_y * tg_yyyyz_xxxxxx_0[j] + fr * tg_yyyyz_xxxxxx_1[j] + 2.0 * fl1_fx * (tg_yyyz_xxxxxx_0[j] - tg_yyyz_xxxxxx_1[j] * fl1_fza);

                    tg_yyyyyz_xxxxxy_0[j] = pb_y * tg_yyyyz_xxxxxy_0[j] + fr * tg_yyyyz_xxxxxy_1[j] + 2.0 * fl1_fx * (tg_yyyz_xxxxxy_0[j] - tg_yyyz_xxxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyz_xxxxx_1[j];

                    tg_yyyyyz_xxxxxz_0[j] = pb_y * tg_yyyyz_xxxxxz_0[j] + fr * tg_yyyyz_xxxxxz_1[j] + 2.0 * fl1_fx * (tg_yyyz_xxxxxz_0[j] - tg_yyyz_xxxxxz_1[j] * fl1_fza);

                    tg_yyyyyz_xxxxyy_0[j] = pb_y * tg_yyyyz_xxxxyy_0[j] + fr * tg_yyyyz_xxxxyy_1[j] + 2.0 * fl1_fx * (tg_yyyz_xxxxyy_0[j] - tg_yyyz_xxxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyyz_xxxxy_1[j];

                    tg_yyyyyz_xxxxyz_0[j] = pb_y * tg_yyyyz_xxxxyz_0[j] + fr * tg_yyyyz_xxxxyz_1[j] + 2.0 * fl1_fx * (tg_yyyz_xxxxyz_0[j] - tg_yyyz_xxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyz_xxxxz_1[j];

                    tg_yyyyyz_xxxxzz_0[j] = pb_y * tg_yyyyz_xxxxzz_0[j] + fr * tg_yyyyz_xxxxzz_1[j] + 2.0 * fl1_fx * (tg_yyyz_xxxxzz_0[j] - tg_yyyz_xxxxzz_1[j] * fl1_fza);

                    tg_yyyyyz_xxxyyy_0[j] = pb_y * tg_yyyyz_xxxyyy_0[j] + fr * tg_yyyyz_xxxyyy_1[j] + 2.0 * fl1_fx * (tg_yyyz_xxxyyy_0[j] - tg_yyyz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyz_xxxyy_1[j];

                    tg_yyyyyz_xxxyyz_0[j] = pb_y * tg_yyyyz_xxxyyz_0[j] + fr * tg_yyyyz_xxxyyz_1[j] + 2.0 * fl1_fx * (tg_yyyz_xxxyyz_0[j] - tg_yyyz_xxxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyz_xxxyz_1[j];

                    tg_yyyyyz_xxxyzz_0[j] = pb_y * tg_yyyyz_xxxyzz_0[j] + fr * tg_yyyyz_xxxyzz_1[j] + 2.0 * fl1_fx * (tg_yyyz_xxxyzz_0[j] - tg_yyyz_xxxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyz_xxxzz_1[j];

                    tg_yyyyyz_xxxzzz_0[j] = pb_y * tg_yyyyz_xxxzzz_0[j] + fr * tg_yyyyz_xxxzzz_1[j] + 2.0 * fl1_fx * (tg_yyyz_xxxzzz_0[j] - tg_yyyz_xxxzzz_1[j] * fl1_fza);

                    tg_yyyyyz_xxyyyy_0[j] = pb_y * tg_yyyyz_xxyyyy_0[j] + fr * tg_yyyyz_xxyyyy_1[j] + 2.0 * fl1_fx * (tg_yyyz_xxyyyy_0[j] - tg_yyyz_xxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyz_xxyyy_1[j];

                    tg_yyyyyz_xxyyyz_0[j] = pb_y * tg_yyyyz_xxyyyz_0[j] + fr * tg_yyyyz_xxyyyz_1[j] + 2.0 * fl1_fx * (tg_yyyz_xxyyyz_0[j] - tg_yyyz_xxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyz_xxyyz_1[j];

                    tg_yyyyyz_xxyyzz_0[j] = pb_y * tg_yyyyz_xxyyzz_0[j] + fr * tg_yyyyz_xxyyzz_1[j] + 2.0 * fl1_fx * (tg_yyyz_xxyyzz_0[j] - tg_yyyz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyz_xxyzz_1[j];

                    tg_yyyyyz_xxyzzz_0[j] = pb_y * tg_yyyyz_xxyzzz_0[j] + fr * tg_yyyyz_xxyzzz_1[j] + 2.0 * fl1_fx * (tg_yyyz_xxyzzz_0[j] - tg_yyyz_xxyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyz_xxzzz_1[j];

                    tg_yyyyyz_xxzzzz_0[j] = pb_y * tg_yyyyz_xxzzzz_0[j] + fr * tg_yyyyz_xxzzzz_1[j] + 2.0 * fl1_fx * (tg_yyyz_xxzzzz_0[j] - tg_yyyz_xxzzzz_1[j] * fl1_fza);

                    tg_yyyyyz_xyyyyy_0[j] = pb_y * tg_yyyyz_xyyyyy_0[j] + fr * tg_yyyyz_xyyyyy_1[j] + 2.0 * fl1_fx * (tg_yyyz_xyyyyy_0[j] - tg_yyyz_xyyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyyz_xyyyy_1[j];

                    tg_yyyyyz_xyyyyz_0[j] = pb_y * tg_yyyyz_xyyyyz_0[j] + fr * tg_yyyyz_xyyyyz_1[j] + 2.0 * fl1_fx * (tg_yyyz_xyyyyz_0[j] - tg_yyyz_xyyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyz_xyyyz_1[j];

                    tg_yyyyyz_xyyyzz_0[j] = pb_y * tg_yyyyz_xyyyzz_0[j] + fr * tg_yyyyz_xyyyzz_1[j] + 2.0 * fl1_fx * (tg_yyyz_xyyyzz_0[j] - tg_yyyz_xyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyz_xyyzz_1[j];

                    tg_yyyyyz_xyyzzz_0[j] = pb_y * tg_yyyyz_xyyzzz_0[j] + fr * tg_yyyyz_xyyzzz_1[j] + 2.0 * fl1_fx * (tg_yyyz_xyyzzz_0[j] - tg_yyyz_xyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyz_xyzzz_1[j];

                    tg_yyyyyz_xyzzzz_0[j] = pb_y * tg_yyyyz_xyzzzz_0[j] + fr * tg_yyyyz_xyzzzz_1[j] + 2.0 * fl1_fx * (tg_yyyz_xyzzzz_0[j] - tg_yyyz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyz_xzzzz_1[j];

                    tg_yyyyyz_xzzzzz_0[j] = pb_y * tg_yyyyz_xzzzzz_0[j] + fr * tg_yyyyz_xzzzzz_1[j] + 2.0 * fl1_fx * (tg_yyyz_xzzzzz_0[j] - tg_yyyz_xzzzzz_1[j] * fl1_fza);

                    tg_yyyyyz_yyyyyy_0[j] = pb_y * tg_yyyyz_yyyyyy_0[j] + fr * tg_yyyyz_yyyyyy_1[j] + 2.0 * fl1_fx * (tg_yyyz_yyyyyy_0[j] - tg_yyyz_yyyyyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yyyyz_yyyyy_1[j];

                    tg_yyyyyz_yyyyyz_0[j] = pb_y * tg_yyyyz_yyyyyz_0[j] + fr * tg_yyyyz_yyyyyz_1[j] + 2.0 * fl1_fx * (tg_yyyz_yyyyyz_0[j] - tg_yyyz_yyyyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyyz_yyyyz_1[j];

                    tg_yyyyyz_yyyyzz_0[j] = pb_y * tg_yyyyz_yyyyzz_0[j] + fr * tg_yyyyz_yyyyzz_1[j] + 2.0 * fl1_fx * (tg_yyyz_yyyyzz_0[j] - tg_yyyz_yyyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyz_yyyzz_1[j];

                    tg_yyyyyz_yyyzzz_0[j] = pb_y * tg_yyyyz_yyyzzz_0[j] + fr * tg_yyyyz_yyyzzz_1[j] + 2.0 * fl1_fx * (tg_yyyz_yyyzzz_0[j] - tg_yyyz_yyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyz_yyzzz_1[j];

                    tg_yyyyyz_yyzzzz_0[j] = pb_y * tg_yyyyz_yyzzzz_0[j] + fr * tg_yyyyz_yyzzzz_1[j] + 2.0 * fl1_fx * (tg_yyyz_yyzzzz_0[j] - tg_yyyz_yyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyz_yzzzz_1[j];

                    tg_yyyyyz_yzzzzz_0[j] = pb_y * tg_yyyyz_yzzzzz_0[j] + fr * tg_yyyyz_yzzzzz_1[j] + 2.0 * fl1_fx * (tg_yyyz_yzzzzz_0[j] - tg_yyyz_yzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyz_zzzzz_1[j];

                    tg_yyyyyz_zzzzzz_0[j] = pb_y * tg_yyyyz_zzzzzz_0[j] + fr * tg_yyyyz_zzzzzz_1[j] + 2.0 * fl1_fx * (tg_yyyz_zzzzzz_0[j] - tg_yyyz_zzzzzz_1[j] * fl1_fza);

                    tg_yyyyzz_xxxxxx_0[j] = pb_y * tg_yyyzz_xxxxxx_0[j] + fr * tg_yyyzz_xxxxxx_1[j] + 1.5 * fl1_fx * (tg_yyzz_xxxxxx_0[j] - tg_yyzz_xxxxxx_1[j] * fl1_fza);

                    tg_yyyyzz_xxxxxy_0[j] = pb_y * tg_yyyzz_xxxxxy_0[j] + fr * tg_yyyzz_xxxxxy_1[j] + 1.5 * fl1_fx * (tg_yyzz_xxxxxy_0[j] - tg_yyzz_xxxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzz_xxxxx_1[j];

                    tg_yyyyzz_xxxxxz_0[j] = pb_y * tg_yyyzz_xxxxxz_0[j] + fr * tg_yyyzz_xxxxxz_1[j] + 1.5 * fl1_fx * (tg_yyzz_xxxxxz_0[j] - tg_yyzz_xxxxxz_1[j] * fl1_fza);

                    tg_yyyyzz_xxxxyy_0[j] = pb_y * tg_yyyzz_xxxxyy_0[j] + fr * tg_yyyzz_xxxxyy_1[j] + 1.5 * fl1_fx * (tg_yyzz_xxxxyy_0[j] - tg_yyzz_xxxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyzz_xxxxy_1[j];

                    tg_yyyyzz_xxxxyz_0[j] = pb_y * tg_yyyzz_xxxxyz_0[j] + fr * tg_yyyzz_xxxxyz_1[j] + 1.5 * fl1_fx * (tg_yyzz_xxxxyz_0[j] - tg_yyzz_xxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzz_xxxxz_1[j];

                    tg_yyyyzz_xxxxzz_0[j] = pb_y * tg_yyyzz_xxxxzz_0[j] + fr * tg_yyyzz_xxxxzz_1[j] + 1.5 * fl1_fx * (tg_yyzz_xxxxzz_0[j] - tg_yyzz_xxxxzz_1[j] * fl1_fza);

                    tg_yyyyzz_xxxyyy_0[j] = pb_y * tg_yyyzz_xxxyyy_0[j] + fr * tg_yyyzz_xxxyyy_1[j] + 1.5 * fl1_fx * (tg_yyzz_xxxyyy_0[j] - tg_yyzz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyzz_xxxyy_1[j];

                    tg_yyyyzz_xxxyyz_0[j] = pb_y * tg_yyyzz_xxxyyz_0[j] + fr * tg_yyyzz_xxxyyz_1[j] + 1.5 * fl1_fx * (tg_yyzz_xxxyyz_0[j] - tg_yyzz_xxxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyzz_xxxyz_1[j];

                    tg_yyyyzz_xxxyzz_0[j] = pb_y * tg_yyyzz_xxxyzz_0[j] + fr * tg_yyyzz_xxxyzz_1[j] + 1.5 * fl1_fx * (tg_yyzz_xxxyzz_0[j] - tg_yyzz_xxxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzz_xxxzz_1[j];

                    tg_yyyyzz_xxxzzz_0[j] = pb_y * tg_yyyzz_xxxzzz_0[j] + fr * tg_yyyzz_xxxzzz_1[j] + 1.5 * fl1_fx * (tg_yyzz_xxxzzz_0[j] - tg_yyzz_xxxzzz_1[j] * fl1_fza);

                    tg_yyyyzz_xxyyyy_0[j] = pb_y * tg_yyyzz_xxyyyy_0[j] + fr * tg_yyyzz_xxyyyy_1[j] + 1.5 * fl1_fx * (tg_yyzz_xxyyyy_0[j] - tg_yyzz_xxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyzz_xxyyy_1[j];

                    tg_yyyyzz_xxyyyz_0[j] = pb_y * tg_yyyzz_xxyyyz_0[j] + fr * tg_yyyzz_xxyyyz_1[j] + 1.5 * fl1_fx * (tg_yyzz_xxyyyz_0[j] - tg_yyzz_xxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyzz_xxyyz_1[j];

                    tg_yyyyzz_xxyyzz_0[j] = pb_y * tg_yyyzz_xxyyzz_0[j] + fr * tg_yyyzz_xxyyzz_1[j] + 1.5 * fl1_fx * (tg_yyzz_xxyyzz_0[j] - tg_yyzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyzz_xxyzz_1[j];

                    tg_yyyyzz_xxyzzz_0[j] = pb_y * tg_yyyzz_xxyzzz_0[j] + fr * tg_yyyzz_xxyzzz_1[j] + 1.5 * fl1_fx * (tg_yyzz_xxyzzz_0[j] - tg_yyzz_xxyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzz_xxzzz_1[j];

                    tg_yyyyzz_xxzzzz_0[j] = pb_y * tg_yyyzz_xxzzzz_0[j] + fr * tg_yyyzz_xxzzzz_1[j] + 1.5 * fl1_fx * (tg_yyzz_xxzzzz_0[j] - tg_yyzz_xxzzzz_1[j] * fl1_fza);

                    tg_yyyyzz_xyyyyy_0[j] = pb_y * tg_yyyzz_xyyyyy_0[j] + fr * tg_yyyzz_xyyyyy_1[j] + 1.5 * fl1_fx * (tg_yyzz_xyyyyy_0[j] - tg_yyzz_xyyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyzz_xyyyy_1[j];

                    tg_yyyyzz_xyyyyz_0[j] = pb_y * tg_yyyzz_xyyyyz_0[j] + fr * tg_yyyzz_xyyyyz_1[j] + 1.5 * fl1_fx * (tg_yyzz_xyyyyz_0[j] - tg_yyzz_xyyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyzz_xyyyz_1[j];

                    tg_yyyyzz_xyyyzz_0[j] = pb_y * tg_yyyzz_xyyyzz_0[j] + fr * tg_yyyzz_xyyyzz_1[j] + 1.5 * fl1_fx * (tg_yyzz_xyyyzz_0[j] - tg_yyzz_xyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyzz_xyyzz_1[j];

                    tg_yyyyzz_xyyzzz_0[j] = pb_y * tg_yyyzz_xyyzzz_0[j] + fr * tg_yyyzz_xyyzzz_1[j] + 1.5 * fl1_fx * (tg_yyzz_xyyzzz_0[j] - tg_yyzz_xyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyzz_xyzzz_1[j];

                    tg_yyyyzz_xyzzzz_0[j] = pb_y * tg_yyyzz_xyzzzz_0[j] + fr * tg_yyyzz_xyzzzz_1[j] + 1.5 * fl1_fx * (tg_yyzz_xyzzzz_0[j] - tg_yyzz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzz_xzzzz_1[j];

                    tg_yyyyzz_xzzzzz_0[j] = pb_y * tg_yyyzz_xzzzzz_0[j] + fr * tg_yyyzz_xzzzzz_1[j] + 1.5 * fl1_fx * (tg_yyzz_xzzzzz_0[j] - tg_yyzz_xzzzzz_1[j] * fl1_fza);

                    tg_yyyyzz_yyyyyy_0[j] = pb_y * tg_yyyzz_yyyyyy_0[j] + fr * tg_yyyzz_yyyyyy_1[j] + 1.5 * fl1_fx * (tg_yyzz_yyyyyy_0[j] - tg_yyzz_yyyyyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yyyzz_yyyyy_1[j];

                    tg_yyyyzz_yyyyyz_0[j] = pb_y * tg_yyyzz_yyyyyz_0[j] + fr * tg_yyyzz_yyyyyz_1[j] + 1.5 * fl1_fx * (tg_yyzz_yyyyyz_0[j] - tg_yyzz_yyyyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyzz_yyyyz_1[j];

                    tg_yyyyzz_yyyyzz_0[j] = pb_y * tg_yyyzz_yyyyzz_0[j] + fr * tg_yyyzz_yyyyzz_1[j] + 1.5 * fl1_fx * (tg_yyzz_yyyyzz_0[j] - tg_yyzz_yyyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyzz_yyyzz_1[j];

                    tg_yyyyzz_yyyzzz_0[j] = pb_y * tg_yyyzz_yyyzzz_0[j] + fr * tg_yyyzz_yyyzzz_1[j] + 1.5 * fl1_fx * (tg_yyzz_yyyzzz_0[j] - tg_yyzz_yyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyzz_yyzzz_1[j];

                    tg_yyyyzz_yyzzzz_0[j] = pb_y * tg_yyyzz_yyzzzz_0[j] + fr * tg_yyyzz_yyzzzz_1[j] + 1.5 * fl1_fx * (tg_yyzz_yyzzzz_0[j] - tg_yyzz_yyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyzz_yzzzz_1[j];

                    tg_yyyyzz_yzzzzz_0[j] = pb_y * tg_yyyzz_yzzzzz_0[j] + fr * tg_yyyzz_yzzzzz_1[j] + 1.5 * fl1_fx * (tg_yyzz_yzzzzz_0[j] - tg_yyzz_yzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzz_zzzzz_1[j];

                    tg_yyyyzz_zzzzzz_0[j] = pb_y * tg_yyyzz_zzzzzz_0[j] + fr * tg_yyyzz_zzzzzz_1[j] + 1.5 * fl1_fx * (tg_yyzz_zzzzzz_0[j] - tg_yyzz_zzzzzz_1[j] * fl1_fza);

                    tg_yyyzzz_xxxxxx_0[j] = pb_y * tg_yyzzz_xxxxxx_0[j] + fr * tg_yyzzz_xxxxxx_1[j] + fl1_fx * (tg_yzzz_xxxxxx_0[j] - tg_yzzz_xxxxxx_1[j] * fl1_fza);

                    tg_yyyzzz_xxxxxy_0[j] = pb_y * tg_yyzzz_xxxxxy_0[j] + fr * tg_yyzzz_xxxxxy_1[j] + fl1_fx * (tg_yzzz_xxxxxy_0[j] - tg_yzzz_xxxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzz_xxxxx_1[j];

                    tg_yyyzzz_xxxxxz_0[j] = pb_y * tg_yyzzz_xxxxxz_0[j] + fr * tg_yyzzz_xxxxxz_1[j] + fl1_fx * (tg_yzzz_xxxxxz_0[j] - tg_yzzz_xxxxxz_1[j] * fl1_fza);

                    tg_yyyzzz_xxxxyy_0[j] = pb_y * tg_yyzzz_xxxxyy_0[j] + fr * tg_yyzzz_xxxxyy_1[j] + fl1_fx * (tg_yzzz_xxxxyy_0[j] - tg_yzzz_xxxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyzzz_xxxxy_1[j];

                    tg_yyyzzz_xxxxyz_0[j] = pb_y * tg_yyzzz_xxxxyz_0[j] + fr * tg_yyzzz_xxxxyz_1[j] + fl1_fx * (tg_yzzz_xxxxyz_0[j] - tg_yzzz_xxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzz_xxxxz_1[j];

                    tg_yyyzzz_xxxxzz_0[j] = pb_y * tg_yyzzz_xxxxzz_0[j] + fr * tg_yyzzz_xxxxzz_1[j] + fl1_fx * (tg_yzzz_xxxxzz_0[j] - tg_yzzz_xxxxzz_1[j] * fl1_fza);

                    tg_yyyzzz_xxxyyy_0[j] = pb_y * tg_yyzzz_xxxyyy_0[j] + fr * tg_yyzzz_xxxyyy_1[j] + fl1_fx * (tg_yzzz_xxxyyy_0[j] - tg_yzzz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzzz_xxxyy_1[j];

                    tg_yyyzzz_xxxyyz_0[j] = pb_y * tg_yyzzz_xxxyyz_0[j] + fr * tg_yyzzz_xxxyyz_1[j] + fl1_fx * (tg_yzzz_xxxyyz_0[j] - tg_yzzz_xxxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyzzz_xxxyz_1[j];

                    tg_yyyzzz_xxxyzz_0[j] = pb_y * tg_yyzzz_xxxyzz_0[j] + fr * tg_yyzzz_xxxyzz_1[j] + fl1_fx * (tg_yzzz_xxxyzz_0[j] - tg_yzzz_xxxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzz_xxxzz_1[j];

                    tg_yyyzzz_xxxzzz_0[j] = pb_y * tg_yyzzz_xxxzzz_0[j] + fr * tg_yyzzz_xxxzzz_1[j] + fl1_fx * (tg_yzzz_xxxzzz_0[j] - tg_yzzz_xxxzzz_1[j] * fl1_fza);

                    tg_yyyzzz_xxyyyy_0[j] = pb_y * tg_yyzzz_xxyyyy_0[j] + fr * tg_yyzzz_xxyyyy_1[j] + fl1_fx * (tg_yzzz_xxyyyy_0[j] - tg_yzzz_xxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyzzz_xxyyy_1[j];

                    tg_yyyzzz_xxyyyz_0[j] = pb_y * tg_yyzzz_xxyyyz_0[j] + fr * tg_yyzzz_xxyyyz_1[j] + fl1_fx * (tg_yzzz_xxyyyz_0[j] - tg_yzzz_xxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzzz_xxyyz_1[j];

                    tg_yyyzzz_xxyyzz_0[j] = pb_y * tg_yyzzz_xxyyzz_0[j] + fr * tg_yyzzz_xxyyzz_1[j] + fl1_fx * (tg_yzzz_xxyyzz_0[j] - tg_yzzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyzzz_xxyzz_1[j];

                    tg_yyyzzz_xxyzzz_0[j] = pb_y * tg_yyzzz_xxyzzz_0[j] + fr * tg_yyzzz_xxyzzz_1[j] + fl1_fx * (tg_yzzz_xxyzzz_0[j] - tg_yzzz_xxyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzz_xxzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSISI_686_784(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (686,784)

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
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_6_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_6_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_5_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
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

                auto tg_yyzzz_xxzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 518); 

                auto tg_yyzzz_xyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 519); 

                auto tg_yyzzz_xyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 520); 

                auto tg_yyzzz_xyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 521); 

                auto tg_yyzzz_xyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 522); 

                auto tg_yyzzz_xyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 523); 

                auto tg_yyzzz_xzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 524); 

                auto tg_yyzzz_yyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 525); 

                auto tg_yyzzz_yyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 526); 

                auto tg_yyzzz_yyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 527); 

                auto tg_yyzzz_yyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 528); 

                auto tg_yyzzz_yyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 529); 

                auto tg_yyzzz_yzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 530); 

                auto tg_yyzzz_zzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 531); 

                auto tg_yzzzz_xxxxxx_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 532); 

                auto tg_yzzzz_xxxxxy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 533); 

                auto tg_yzzzz_xxxxxz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 534); 

                auto tg_yzzzz_xxxxyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 535); 

                auto tg_yzzzz_xxxxyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 536); 

                auto tg_yzzzz_xxxxzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 537); 

                auto tg_yzzzz_xxxyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 538); 

                auto tg_yzzzz_xxxyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 539); 

                auto tg_yzzzz_xxxyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 540); 

                auto tg_yzzzz_xxxzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 541); 

                auto tg_yzzzz_xxyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 542); 

                auto tg_yzzzz_xxyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 543); 

                auto tg_yzzzz_xxyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 544); 

                auto tg_yzzzz_xxyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 545); 

                auto tg_yzzzz_xxzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 546); 

                auto tg_yzzzz_xyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 547); 

                auto tg_yzzzz_xyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 548); 

                auto tg_yzzzz_xyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 549); 

                auto tg_yzzzz_xyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 550); 

                auto tg_yzzzz_xyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 551); 

                auto tg_yzzzz_xzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 552); 

                auto tg_yzzzz_yyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 553); 

                auto tg_yzzzz_yyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 554); 

                auto tg_yzzzz_yyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 555); 

                auto tg_yzzzz_yyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 556); 

                auto tg_yzzzz_yyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 557); 

                auto tg_yzzzz_yzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 558); 

                auto tg_yzzzz_zzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 559); 

                auto tg_zzzzz_xxxxxx_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 560); 

                auto tg_zzzzz_xxxxxy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 561); 

                auto tg_zzzzz_xxxxxz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 562); 

                auto tg_zzzzz_xxxxyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 563); 

                auto tg_zzzzz_xxxxyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 564); 

                auto tg_zzzzz_xxxxzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 565); 

                auto tg_zzzzz_xxxyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 566); 

                auto tg_zzzzz_xxxyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 567); 

                auto tg_zzzzz_xxxyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 568); 

                auto tg_zzzzz_xxxzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 569); 

                auto tg_zzzzz_xxyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 570); 

                auto tg_zzzzz_xxyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 571); 

                auto tg_zzzzz_xxyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 572); 

                auto tg_zzzzz_xxyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 573); 

                auto tg_zzzzz_xxzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 574); 

                auto tg_zzzzz_xyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 575); 

                auto tg_zzzzz_xyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 576); 

                auto tg_zzzzz_xyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 577); 

                auto tg_zzzzz_xyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 578); 

                auto tg_zzzzz_xyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 579); 

                auto tg_zzzzz_xzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 580); 

                auto tg_zzzzz_yyyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 581); 

                auto tg_zzzzz_yyyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 582); 

                auto tg_zzzzz_yyyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 583); 

                auto tg_zzzzz_yyyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 584); 

                auto tg_zzzzz_yyzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 585); 

                auto tg_zzzzz_yzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 586); 

                auto tg_zzzzz_zzzzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 587); 

                auto tg_yyzzz_xxzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 518); 

                auto tg_yyzzz_xyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 519); 

                auto tg_yyzzz_xyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 520); 

                auto tg_yyzzz_xyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 521); 

                auto tg_yyzzz_xyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 522); 

                auto tg_yyzzz_xyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 523); 

                auto tg_yyzzz_xzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 524); 

                auto tg_yyzzz_yyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 525); 

                auto tg_yyzzz_yyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 526); 

                auto tg_yyzzz_yyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 527); 

                auto tg_yyzzz_yyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 528); 

                auto tg_yyzzz_yyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 529); 

                auto tg_yyzzz_yzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 530); 

                auto tg_yyzzz_zzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 531); 

                auto tg_yzzzz_xxxxxx_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 532); 

                auto tg_yzzzz_xxxxxy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 533); 

                auto tg_yzzzz_xxxxxz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 534); 

                auto tg_yzzzz_xxxxyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 535); 

                auto tg_yzzzz_xxxxyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 536); 

                auto tg_yzzzz_xxxxzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 537); 

                auto tg_yzzzz_xxxyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 538); 

                auto tg_yzzzz_xxxyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 539); 

                auto tg_yzzzz_xxxyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 540); 

                auto tg_yzzzz_xxxzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 541); 

                auto tg_yzzzz_xxyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 542); 

                auto tg_yzzzz_xxyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 543); 

                auto tg_yzzzz_xxyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 544); 

                auto tg_yzzzz_xxyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 545); 

                auto tg_yzzzz_xxzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 546); 

                auto tg_yzzzz_xyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 547); 

                auto tg_yzzzz_xyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 548); 

                auto tg_yzzzz_xyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 549); 

                auto tg_yzzzz_xyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 550); 

                auto tg_yzzzz_xyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 551); 

                auto tg_yzzzz_xzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 552); 

                auto tg_yzzzz_yyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 553); 

                auto tg_yzzzz_yyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 554); 

                auto tg_yzzzz_yyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 555); 

                auto tg_yzzzz_yyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 556); 

                auto tg_yzzzz_yyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 557); 

                auto tg_yzzzz_yzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 558); 

                auto tg_yzzzz_zzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 559); 

                auto tg_zzzzz_xxxxxx_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 560); 

                auto tg_zzzzz_xxxxxy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 561); 

                auto tg_zzzzz_xxxxxz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 562); 

                auto tg_zzzzz_xxxxyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 563); 

                auto tg_zzzzz_xxxxyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 564); 

                auto tg_zzzzz_xxxxzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 565); 

                auto tg_zzzzz_xxxyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 566); 

                auto tg_zzzzz_xxxyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 567); 

                auto tg_zzzzz_xxxyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 568); 

                auto tg_zzzzz_xxxzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 569); 

                auto tg_zzzzz_xxyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 570); 

                auto tg_zzzzz_xxyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 571); 

                auto tg_zzzzz_xxyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 572); 

                auto tg_zzzzz_xxyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 573); 

                auto tg_zzzzz_xxzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 574); 

                auto tg_zzzzz_xyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 575); 

                auto tg_zzzzz_xyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 576); 

                auto tg_zzzzz_xyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 577); 

                auto tg_zzzzz_xyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 578); 

                auto tg_zzzzz_xyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 579); 

                auto tg_zzzzz_xzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 580); 

                auto tg_zzzzz_yyyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 581); 

                auto tg_zzzzz_yyyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 582); 

                auto tg_zzzzz_yyyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 583); 

                auto tg_zzzzz_yyyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 584); 

                auto tg_zzzzz_yyzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 585); 

                auto tg_zzzzz_yzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 586); 

                auto tg_zzzzz_zzzzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 587); 

                auto tg_yzzz_xxzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 378); 

                auto tg_yzzz_xyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 379); 

                auto tg_yzzz_xyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 380); 

                auto tg_yzzz_xyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 381); 

                auto tg_yzzz_xyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 382); 

                auto tg_yzzz_xyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 383); 

                auto tg_yzzz_xzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 384); 

                auto tg_yzzz_yyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 385); 

                auto tg_yzzz_yyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 386); 

                auto tg_yzzz_yyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 387); 

                auto tg_yzzz_yyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 388); 

                auto tg_yzzz_yyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 389); 

                auto tg_yzzz_yzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 390); 

                auto tg_yzzz_zzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 391); 

                auto tg_zzzz_xxxxxx_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 392); 

                auto tg_zzzz_xxxxxy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 393); 

                auto tg_zzzz_xxxxxz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 394); 

                auto tg_zzzz_xxxxyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 395); 

                auto tg_zzzz_xxxxyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 396); 

                auto tg_zzzz_xxxxzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 397); 

                auto tg_zzzz_xxxyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 398); 

                auto tg_zzzz_xxxyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 399); 

                auto tg_zzzz_xxxyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 400); 

                auto tg_zzzz_xxxzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 401); 

                auto tg_zzzz_xxyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 402); 

                auto tg_zzzz_xxyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 403); 

                auto tg_zzzz_xxyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 404); 

                auto tg_zzzz_xxyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 405); 

                auto tg_zzzz_xxzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 406); 

                auto tg_zzzz_xyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 407); 

                auto tg_zzzz_xyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 408); 

                auto tg_zzzz_xyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 409); 

                auto tg_zzzz_xyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 410); 

                auto tg_zzzz_xyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 411); 

                auto tg_zzzz_xzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 412); 

                auto tg_zzzz_yyyyyy_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 413); 

                auto tg_zzzz_yyyyyz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 414); 

                auto tg_zzzz_yyyyzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 415); 

                auto tg_zzzz_yyyzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 416); 

                auto tg_zzzz_yyzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 417); 

                auto tg_zzzz_yzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 418); 

                auto tg_zzzz_zzzzzz_0 = primBuffer[pidx_g_4_6_m0].data(420 * idx + 419); 

                auto tg_yzzz_xxzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 378); 

                auto tg_yzzz_xyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 379); 

                auto tg_yzzz_xyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 380); 

                auto tg_yzzz_xyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 381); 

                auto tg_yzzz_xyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 382); 

                auto tg_yzzz_xyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 383); 

                auto tg_yzzz_xzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 384); 

                auto tg_yzzz_yyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 385); 

                auto tg_yzzz_yyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 386); 

                auto tg_yzzz_yyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 387); 

                auto tg_yzzz_yyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 388); 

                auto tg_yzzz_yyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 389); 

                auto tg_yzzz_yzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 390); 

                auto tg_yzzz_zzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 391); 

                auto tg_zzzz_xxxxxx_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 392); 

                auto tg_zzzz_xxxxxy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 393); 

                auto tg_zzzz_xxxxxz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 394); 

                auto tg_zzzz_xxxxyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 395); 

                auto tg_zzzz_xxxxyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 396); 

                auto tg_zzzz_xxxxzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 397); 

                auto tg_zzzz_xxxyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 398); 

                auto tg_zzzz_xxxyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 399); 

                auto tg_zzzz_xxxyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 400); 

                auto tg_zzzz_xxxzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 401); 

                auto tg_zzzz_xxyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 402); 

                auto tg_zzzz_xxyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 403); 

                auto tg_zzzz_xxyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 404); 

                auto tg_zzzz_xxyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 405); 

                auto tg_zzzz_xxzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 406); 

                auto tg_zzzz_xyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 407); 

                auto tg_zzzz_xyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 408); 

                auto tg_zzzz_xyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 409); 

                auto tg_zzzz_xyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 410); 

                auto tg_zzzz_xyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 411); 

                auto tg_zzzz_xzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 412); 

                auto tg_zzzz_yyyyyy_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 413); 

                auto tg_zzzz_yyyyyz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 414); 

                auto tg_zzzz_yyyyzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 415); 

                auto tg_zzzz_yyyzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 416); 

                auto tg_zzzz_yyzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 417); 

                auto tg_zzzz_yzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 418); 

                auto tg_zzzz_zzzzzz_1 = primBuffer[pidx_g_4_6_m1].data(420 * idx + 419); 

                auto tg_yyzzz_xyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 388); 

                auto tg_yyzzz_xyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 389); 

                auto tg_yyzzz_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 390); 

                auto tg_yyzzz_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 391); 

                auto tg_yyzzz_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 392); 

                auto tg_yyzzz_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 393); 

                auto tg_yyzzz_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 394); 

                auto tg_yyzzz_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 395); 

                auto tg_yyzzz_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 396); 

                auto tg_yyzzz_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 397); 

                auto tg_yyzzz_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 398); 

                auto tg_yzzzz_xxxxx_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 399); 

                auto tg_yzzzz_xxxxy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 400); 

                auto tg_yzzzz_xxxxz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 401); 

                auto tg_yzzzz_xxxyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 402); 

                auto tg_yzzzz_xxxyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 403); 

                auto tg_yzzzz_xxxzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 404); 

                auto tg_yzzzz_xxyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 405); 

                auto tg_yzzzz_xxyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 406); 

                auto tg_yzzzz_xxyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 407); 

                auto tg_yzzzz_xxzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 408); 

                auto tg_yzzzz_xyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 409); 

                auto tg_yzzzz_xyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 410); 

                auto tg_yzzzz_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 411); 

                auto tg_yzzzz_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 412); 

                auto tg_yzzzz_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 413); 

                auto tg_yzzzz_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 414); 

                auto tg_yzzzz_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 415); 

                auto tg_yzzzz_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 416); 

                auto tg_yzzzz_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 417); 

                auto tg_yzzzz_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 418); 

                auto tg_yzzzz_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 419); 

                auto tg_zzzzz_xxxxx_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 420); 

                auto tg_zzzzz_xxxxy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 421); 

                auto tg_zzzzz_xxxxz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 422); 

                auto tg_zzzzz_xxxyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 423); 

                auto tg_zzzzz_xxxyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 424); 

                auto tg_zzzzz_xxxzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 425); 

                auto tg_zzzzz_xxyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 426); 

                auto tg_zzzzz_xxyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 427); 

                auto tg_zzzzz_xxyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 428); 

                auto tg_zzzzz_xxzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 429); 

                auto tg_zzzzz_xyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 430); 

                auto tg_zzzzz_xyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 431); 

                auto tg_zzzzz_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 432); 

                auto tg_zzzzz_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 433); 

                auto tg_zzzzz_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 434); 

                auto tg_zzzzz_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 435); 

                auto tg_zzzzz_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 436); 

                auto tg_zzzzz_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 437); 

                auto tg_zzzzz_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 438); 

                auto tg_zzzzz_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 439); 

                auto tg_zzzzz_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 440); 

                // set up pointers to integrals

                auto tg_yyyzzz_xxzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 686); 

                auto tg_yyyzzz_xyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 687); 

                auto tg_yyyzzz_xyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 688); 

                auto tg_yyyzzz_xyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 689); 

                auto tg_yyyzzz_xyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 690); 

                auto tg_yyyzzz_xyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 691); 

                auto tg_yyyzzz_xzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 692); 

                auto tg_yyyzzz_yyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 693); 

                auto tg_yyyzzz_yyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 694); 

                auto tg_yyyzzz_yyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 695); 

                auto tg_yyyzzz_yyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 696); 

                auto tg_yyyzzz_yyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 697); 

                auto tg_yyyzzz_yzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 698); 

                auto tg_yyyzzz_zzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 699); 

                auto tg_yyzzzz_xxxxxx_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 700); 

                auto tg_yyzzzz_xxxxxy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 701); 

                auto tg_yyzzzz_xxxxxz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 702); 

                auto tg_yyzzzz_xxxxyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 703); 

                auto tg_yyzzzz_xxxxyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 704); 

                auto tg_yyzzzz_xxxxzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 705); 

                auto tg_yyzzzz_xxxyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 706); 

                auto tg_yyzzzz_xxxyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 707); 

                auto tg_yyzzzz_xxxyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 708); 

                auto tg_yyzzzz_xxxzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 709); 

                auto tg_yyzzzz_xxyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 710); 

                auto tg_yyzzzz_xxyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 711); 

                auto tg_yyzzzz_xxyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 712); 

                auto tg_yyzzzz_xxyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 713); 

                auto tg_yyzzzz_xxzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 714); 

                auto tg_yyzzzz_xyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 715); 

                auto tg_yyzzzz_xyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 716); 

                auto tg_yyzzzz_xyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 717); 

                auto tg_yyzzzz_xyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 718); 

                auto tg_yyzzzz_xyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 719); 

                auto tg_yyzzzz_xzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 720); 

                auto tg_yyzzzz_yyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 721); 

                auto tg_yyzzzz_yyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 722); 

                auto tg_yyzzzz_yyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 723); 

                auto tg_yyzzzz_yyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 724); 

                auto tg_yyzzzz_yyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 725); 

                auto tg_yyzzzz_yzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 726); 

                auto tg_yyzzzz_zzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 727); 

                auto tg_yzzzzz_xxxxxx_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 728); 

                auto tg_yzzzzz_xxxxxy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 729); 

                auto tg_yzzzzz_xxxxxz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 730); 

                auto tg_yzzzzz_xxxxyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 731); 

                auto tg_yzzzzz_xxxxyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 732); 

                auto tg_yzzzzz_xxxxzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 733); 

                auto tg_yzzzzz_xxxyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 734); 

                auto tg_yzzzzz_xxxyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 735); 

                auto tg_yzzzzz_xxxyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 736); 

                auto tg_yzzzzz_xxxzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 737); 

                auto tg_yzzzzz_xxyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 738); 

                auto tg_yzzzzz_xxyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 739); 

                auto tg_yzzzzz_xxyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 740); 

                auto tg_yzzzzz_xxyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 741); 

                auto tg_yzzzzz_xxzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 742); 

                auto tg_yzzzzz_xyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 743); 

                auto tg_yzzzzz_xyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 744); 

                auto tg_yzzzzz_xyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 745); 

                auto tg_yzzzzz_xyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 746); 

                auto tg_yzzzzz_xyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 747); 

                auto tg_yzzzzz_xzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 748); 

                auto tg_yzzzzz_yyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 749); 

                auto tg_yzzzzz_yyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 750); 

                auto tg_yzzzzz_yyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 751); 

                auto tg_yzzzzz_yyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 752); 

                auto tg_yzzzzz_yyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 753); 

                auto tg_yzzzzz_yzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 754); 

                auto tg_yzzzzz_zzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 755); 

                auto tg_zzzzzz_xxxxxx_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 756); 

                auto tg_zzzzzz_xxxxxy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 757); 

                auto tg_zzzzzz_xxxxxz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 758); 

                auto tg_zzzzzz_xxxxyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 759); 

                auto tg_zzzzzz_xxxxyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 760); 

                auto tg_zzzzzz_xxxxzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 761); 

                auto tg_zzzzzz_xxxyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 762); 

                auto tg_zzzzzz_xxxyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 763); 

                auto tg_zzzzzz_xxxyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 764); 

                auto tg_zzzzzz_xxxzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 765); 

                auto tg_zzzzzz_xxyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 766); 

                auto tg_zzzzzz_xxyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 767); 

                auto tg_zzzzzz_xxyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 768); 

                auto tg_zzzzzz_xxyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 769); 

                auto tg_zzzzzz_xxzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 770); 

                auto tg_zzzzzz_xyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 771); 

                auto tg_zzzzzz_xyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 772); 

                auto tg_zzzzzz_xyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 773); 

                auto tg_zzzzzz_xyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 774); 

                auto tg_zzzzzz_xyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 775); 

                auto tg_zzzzzz_xzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 776); 

                auto tg_zzzzzz_yyyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 777); 

                auto tg_zzzzzz_yyyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 778); 

                auto tg_zzzzzz_yyyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 779); 

                auto tg_zzzzzz_yyyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 780); 

                auto tg_zzzzzz_yyzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 781); 

                auto tg_zzzzzz_yzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 782); 

                auto tg_zzzzzz_zzzzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 783); 

                // Batch of Integrals (686,784)

                #pragma omp simd aligned(fxn, fza, tg_yyyzzz_xxzzzz_0, tg_yyyzzz_xyyyyy_0, tg_yyyzzz_xyyyyz_0, \
                                         tg_yyyzzz_xyyyzz_0, tg_yyyzzz_xyyzzz_0, tg_yyyzzz_xyzzzz_0, tg_yyyzzz_xzzzzz_0, \
                                         tg_yyyzzz_yyyyyy_0, tg_yyyzzz_yyyyyz_0, tg_yyyzzz_yyyyzz_0, tg_yyyzzz_yyyzzz_0, \
                                         tg_yyyzzz_yyzzzz_0, tg_yyyzzz_yzzzzz_0, tg_yyyzzz_zzzzzz_0, tg_yyzzz_xxzzzz_0, \
                                         tg_yyzzz_xxzzzz_1, tg_yyzzz_xyyyy_1, tg_yyzzz_xyyyyy_0, tg_yyzzz_xyyyyy_1, \
                                         tg_yyzzz_xyyyyz_0, tg_yyzzz_xyyyyz_1, tg_yyzzz_xyyyz_1, tg_yyzzz_xyyyzz_0, \
                                         tg_yyzzz_xyyyzz_1, tg_yyzzz_xyyzz_1, tg_yyzzz_xyyzzz_0, tg_yyzzz_xyyzzz_1, \
                                         tg_yyzzz_xyzzz_1, tg_yyzzz_xyzzzz_0, tg_yyzzz_xyzzzz_1, tg_yyzzz_xzzzz_1, \
                                         tg_yyzzz_xzzzzz_0, tg_yyzzz_xzzzzz_1, tg_yyzzz_yyyyy_1, tg_yyzzz_yyyyyy_0, \
                                         tg_yyzzz_yyyyyy_1, tg_yyzzz_yyyyyz_0, tg_yyzzz_yyyyyz_1, tg_yyzzz_yyyyz_1, \
                                         tg_yyzzz_yyyyzz_0, tg_yyzzz_yyyyzz_1, tg_yyzzz_yyyzz_1, tg_yyzzz_yyyzzz_0, \
                                         tg_yyzzz_yyyzzz_1, tg_yyzzz_yyzzz_1, tg_yyzzz_yyzzzz_0, tg_yyzzz_yyzzzz_1, \
                                         tg_yyzzz_yzzzz_1, tg_yyzzz_yzzzzz_0, tg_yyzzz_yzzzzz_1, tg_yyzzz_zzzzz_1, \
                                         tg_yyzzz_zzzzzz_0, tg_yyzzz_zzzzzz_1, tg_yyzzzz_xxxxxx_0, tg_yyzzzz_xxxxxy_0, \
                                         tg_yyzzzz_xxxxxz_0, tg_yyzzzz_xxxxyy_0, tg_yyzzzz_xxxxyz_0, tg_yyzzzz_xxxxzz_0, \
                                         tg_yyzzzz_xxxyyy_0, tg_yyzzzz_xxxyyz_0, tg_yyzzzz_xxxyzz_0, tg_yyzzzz_xxxzzz_0, \
                                         tg_yyzzzz_xxyyyy_0, tg_yyzzzz_xxyyyz_0, tg_yyzzzz_xxyyzz_0, tg_yyzzzz_xxyzzz_0, \
                                         tg_yyzzzz_xxzzzz_0, tg_yyzzzz_xyyyyy_0, tg_yyzzzz_xyyyyz_0, tg_yyzzzz_xyyyzz_0, \
                                         tg_yyzzzz_xyyzzz_0, tg_yyzzzz_xyzzzz_0, tg_yyzzzz_xzzzzz_0, tg_yyzzzz_yyyyyy_0, \
                                         tg_yyzzzz_yyyyyz_0, tg_yyzzzz_yyyyzz_0, tg_yyzzzz_yyyzzz_0, tg_yyzzzz_yyzzzz_0, \
                                         tg_yyzzzz_yzzzzz_0, tg_yyzzzz_zzzzzz_0, tg_yzzz_xxzzzz_0, tg_yzzz_xxzzzz_1, \
                                         tg_yzzz_xyyyyy_0, tg_yzzz_xyyyyy_1, tg_yzzz_xyyyyz_0, tg_yzzz_xyyyyz_1, \
                                         tg_yzzz_xyyyzz_0, tg_yzzz_xyyyzz_1, tg_yzzz_xyyzzz_0, tg_yzzz_xyyzzz_1, \
                                         tg_yzzz_xyzzzz_0, tg_yzzz_xyzzzz_1, tg_yzzz_xzzzzz_0, tg_yzzz_xzzzzz_1, \
                                         tg_yzzz_yyyyyy_0, tg_yzzz_yyyyyy_1, tg_yzzz_yyyyyz_0, tg_yzzz_yyyyyz_1, \
                                         tg_yzzz_yyyyzz_0, tg_yzzz_yyyyzz_1, tg_yzzz_yyyzzz_0, tg_yzzz_yyyzzz_1, \
                                         tg_yzzz_yyzzzz_0, tg_yzzz_yyzzzz_1, tg_yzzz_yzzzzz_0, tg_yzzz_yzzzzz_1, \
                                         tg_yzzz_zzzzzz_0, tg_yzzz_zzzzzz_1, tg_yzzzz_xxxxx_1, tg_yzzzz_xxxxxx_0, \
                                         tg_yzzzz_xxxxxx_1, tg_yzzzz_xxxxxy_0, tg_yzzzz_xxxxxy_1, tg_yzzzz_xxxxxz_0, \
                                         tg_yzzzz_xxxxxz_1, tg_yzzzz_xxxxy_1, tg_yzzzz_xxxxyy_0, tg_yzzzz_xxxxyy_1, \
                                         tg_yzzzz_xxxxyz_0, tg_yzzzz_xxxxyz_1, tg_yzzzz_xxxxz_1, tg_yzzzz_xxxxzz_0, \
                                         tg_yzzzz_xxxxzz_1, tg_yzzzz_xxxyy_1, tg_yzzzz_xxxyyy_0, tg_yzzzz_xxxyyy_1, \
                                         tg_yzzzz_xxxyyz_0, tg_yzzzz_xxxyyz_1, tg_yzzzz_xxxyz_1, tg_yzzzz_xxxyzz_0, \
                                         tg_yzzzz_xxxyzz_1, tg_yzzzz_xxxzz_1, tg_yzzzz_xxxzzz_0, tg_yzzzz_xxxzzz_1, \
                                         tg_yzzzz_xxyyy_1, tg_yzzzz_xxyyyy_0, tg_yzzzz_xxyyyy_1, tg_yzzzz_xxyyyz_0, \
                                         tg_yzzzz_xxyyyz_1, tg_yzzzz_xxyyz_1, tg_yzzzz_xxyyzz_0, tg_yzzzz_xxyyzz_1, \
                                         tg_yzzzz_xxyzz_1, tg_yzzzz_xxyzzz_0, tg_yzzzz_xxyzzz_1, tg_yzzzz_xxzzz_1, \
                                         tg_yzzzz_xxzzzz_0, tg_yzzzz_xxzzzz_1, tg_yzzzz_xyyyy_1, tg_yzzzz_xyyyyy_0, \
                                         tg_yzzzz_xyyyyy_1, tg_yzzzz_xyyyyz_0, tg_yzzzz_xyyyyz_1, tg_yzzzz_xyyyz_1, \
                                         tg_yzzzz_xyyyzz_0, tg_yzzzz_xyyyzz_1, tg_yzzzz_xyyzz_1, tg_yzzzz_xyyzzz_0, \
                                         tg_yzzzz_xyyzzz_1, tg_yzzzz_xyzzz_1, tg_yzzzz_xyzzzz_0, tg_yzzzz_xyzzzz_1, \
                                         tg_yzzzz_xzzzz_1, tg_yzzzz_xzzzzz_0, tg_yzzzz_xzzzzz_1, tg_yzzzz_yyyyy_1, \
                                         tg_yzzzz_yyyyyy_0, tg_yzzzz_yyyyyy_1, tg_yzzzz_yyyyyz_0, tg_yzzzz_yyyyyz_1, \
                                         tg_yzzzz_yyyyz_1, tg_yzzzz_yyyyzz_0, tg_yzzzz_yyyyzz_1, tg_yzzzz_yyyzz_1, \
                                         tg_yzzzz_yyyzzz_0, tg_yzzzz_yyyzzz_1, tg_yzzzz_yyzzz_1, tg_yzzzz_yyzzzz_0, \
                                         tg_yzzzz_yyzzzz_1, tg_yzzzz_yzzzz_1, tg_yzzzz_yzzzzz_0, tg_yzzzz_yzzzzz_1, \
                                         tg_yzzzz_zzzzz_1, tg_yzzzz_zzzzzz_0, tg_yzzzz_zzzzzz_1, tg_yzzzzz_xxxxxx_0, \
                                         tg_yzzzzz_xxxxxy_0, tg_yzzzzz_xxxxxz_0, tg_yzzzzz_xxxxyy_0, tg_yzzzzz_xxxxyz_0, \
                                         tg_yzzzzz_xxxxzz_0, tg_yzzzzz_xxxyyy_0, tg_yzzzzz_xxxyyz_0, tg_yzzzzz_xxxyzz_0, \
                                         tg_yzzzzz_xxxzzz_0, tg_yzzzzz_xxyyyy_0, tg_yzzzzz_xxyyyz_0, tg_yzzzzz_xxyyzz_0, \
                                         tg_yzzzzz_xxyzzz_0, tg_yzzzzz_xxzzzz_0, tg_yzzzzz_xyyyyy_0, tg_yzzzzz_xyyyyz_0, \
                                         tg_yzzzzz_xyyyzz_0, tg_yzzzzz_xyyzzz_0, tg_yzzzzz_xyzzzz_0, tg_yzzzzz_xzzzzz_0, \
                                         tg_yzzzzz_yyyyyy_0, tg_yzzzzz_yyyyyz_0, tg_yzzzzz_yyyyzz_0, tg_yzzzzz_yyyzzz_0, \
                                         tg_yzzzzz_yyzzzz_0, tg_yzzzzz_yzzzzz_0, tg_yzzzzz_zzzzzz_0, tg_zzzz_xxxxxx_0, \
                                         tg_zzzz_xxxxxx_1, tg_zzzz_xxxxxy_0, tg_zzzz_xxxxxy_1, tg_zzzz_xxxxxz_0, \
                                         tg_zzzz_xxxxxz_1, tg_zzzz_xxxxyy_0, tg_zzzz_xxxxyy_1, tg_zzzz_xxxxyz_0, \
                                         tg_zzzz_xxxxyz_1, tg_zzzz_xxxxzz_0, tg_zzzz_xxxxzz_1, tg_zzzz_xxxyyy_0, \
                                         tg_zzzz_xxxyyy_1, tg_zzzz_xxxyyz_0, tg_zzzz_xxxyyz_1, tg_zzzz_xxxyzz_0, \
                                         tg_zzzz_xxxyzz_1, tg_zzzz_xxxzzz_0, tg_zzzz_xxxzzz_1, tg_zzzz_xxyyyy_0, \
                                         tg_zzzz_xxyyyy_1, tg_zzzz_xxyyyz_0, tg_zzzz_xxyyyz_1, tg_zzzz_xxyyzz_0, \
                                         tg_zzzz_xxyyzz_1, tg_zzzz_xxyzzz_0, tg_zzzz_xxyzzz_1, tg_zzzz_xxzzzz_0, \
                                         tg_zzzz_xxzzzz_1, tg_zzzz_xyyyyy_0, tg_zzzz_xyyyyy_1, tg_zzzz_xyyyyz_0, \
                                         tg_zzzz_xyyyyz_1, tg_zzzz_xyyyzz_0, tg_zzzz_xyyyzz_1, tg_zzzz_xyyzzz_0, \
                                         tg_zzzz_xyyzzz_1, tg_zzzz_xyzzzz_0, tg_zzzz_xyzzzz_1, tg_zzzz_xzzzzz_0, \
                                         tg_zzzz_xzzzzz_1, tg_zzzz_yyyyyy_0, tg_zzzz_yyyyyy_1, tg_zzzz_yyyyyz_0, \
                                         tg_zzzz_yyyyyz_1, tg_zzzz_yyyyzz_0, tg_zzzz_yyyyzz_1, tg_zzzz_yyyzzz_0, \
                                         tg_zzzz_yyyzzz_1, tg_zzzz_yyzzzz_0, tg_zzzz_yyzzzz_1, tg_zzzz_yzzzzz_0, \
                                         tg_zzzz_yzzzzz_1, tg_zzzz_zzzzzz_0, tg_zzzz_zzzzzz_1, tg_zzzzz_xxxxx_1, \
                                         tg_zzzzz_xxxxxx_0, tg_zzzzz_xxxxxx_1, tg_zzzzz_xxxxxy_0, tg_zzzzz_xxxxxy_1, \
                                         tg_zzzzz_xxxxxz_0, tg_zzzzz_xxxxxz_1, tg_zzzzz_xxxxy_1, tg_zzzzz_xxxxyy_0, \
                                         tg_zzzzz_xxxxyy_1, tg_zzzzz_xxxxyz_0, tg_zzzzz_xxxxyz_1, tg_zzzzz_xxxxz_1, \
                                         tg_zzzzz_xxxxzz_0, tg_zzzzz_xxxxzz_1, tg_zzzzz_xxxyy_1, tg_zzzzz_xxxyyy_0, \
                                         tg_zzzzz_xxxyyy_1, tg_zzzzz_xxxyyz_0, tg_zzzzz_xxxyyz_1, tg_zzzzz_xxxyz_1, \
                                         tg_zzzzz_xxxyzz_0, tg_zzzzz_xxxyzz_1, tg_zzzzz_xxxzz_1, tg_zzzzz_xxxzzz_0, \
                                         tg_zzzzz_xxxzzz_1, tg_zzzzz_xxyyy_1, tg_zzzzz_xxyyyy_0, tg_zzzzz_xxyyyy_1, \
                                         tg_zzzzz_xxyyyz_0, tg_zzzzz_xxyyyz_1, tg_zzzzz_xxyyz_1, tg_zzzzz_xxyyzz_0, \
                                         tg_zzzzz_xxyyzz_1, tg_zzzzz_xxyzz_1, tg_zzzzz_xxyzzz_0, tg_zzzzz_xxyzzz_1, \
                                         tg_zzzzz_xxzzz_1, tg_zzzzz_xxzzzz_0, tg_zzzzz_xxzzzz_1, tg_zzzzz_xyyyy_1, \
                                         tg_zzzzz_xyyyyy_0, tg_zzzzz_xyyyyy_1, tg_zzzzz_xyyyyz_0, tg_zzzzz_xyyyyz_1, \
                                         tg_zzzzz_xyyyz_1, tg_zzzzz_xyyyzz_0, tg_zzzzz_xyyyzz_1, tg_zzzzz_xyyzz_1, \
                                         tg_zzzzz_xyyzzz_0, tg_zzzzz_xyyzzz_1, tg_zzzzz_xyzzz_1, tg_zzzzz_xyzzzz_0, \
                                         tg_zzzzz_xyzzzz_1, tg_zzzzz_xzzzz_1, tg_zzzzz_xzzzzz_0, tg_zzzzz_xzzzzz_1, \
                                         tg_zzzzz_yyyyy_1, tg_zzzzz_yyyyyy_0, tg_zzzzz_yyyyyy_1, tg_zzzzz_yyyyyz_0, \
                                         tg_zzzzz_yyyyyz_1, tg_zzzzz_yyyyz_1, tg_zzzzz_yyyyzz_0, tg_zzzzz_yyyyzz_1, \
                                         tg_zzzzz_yyyzz_1, tg_zzzzz_yyyzzz_0, tg_zzzzz_yyyzzz_1, tg_zzzzz_yyzzz_1, \
                                         tg_zzzzz_yyzzzz_0, tg_zzzzz_yyzzzz_1, tg_zzzzz_yzzzz_1, tg_zzzzz_yzzzzz_0, \
                                         tg_zzzzz_yzzzzz_1, tg_zzzzz_zzzzz_1, tg_zzzzz_zzzzzz_0, tg_zzzzz_zzzzzz_1, \
                                         tg_zzzzzz_xxxxxx_0, tg_zzzzzz_xxxxxy_0, tg_zzzzzz_xxxxxz_0, tg_zzzzzz_xxxxyy_0, \
                                         tg_zzzzzz_xxxxyz_0, tg_zzzzzz_xxxxzz_0, tg_zzzzzz_xxxyyy_0, tg_zzzzzz_xxxyyz_0, \
                                         tg_zzzzzz_xxxyzz_0, tg_zzzzzz_xxxzzz_0, tg_zzzzzz_xxyyyy_0, tg_zzzzzz_xxyyyz_0, \
                                         tg_zzzzzz_xxyyzz_0, tg_zzzzzz_xxyzzz_0, tg_zzzzzz_xxzzzz_0, tg_zzzzzz_xyyyyy_0, \
                                         tg_zzzzzz_xyyyyz_0, tg_zzzzzz_xyyyzz_0, tg_zzzzzz_xyyzzz_0, tg_zzzzzz_xyzzzz_0, \
                                         tg_zzzzzz_xzzzzz_0, tg_zzzzzz_yyyyyy_0, tg_zzzzzz_yyyyyz_0, tg_zzzzzz_yyyyzz_0, \
                                         tg_zzzzzz_yyyzzz_0, tg_zzzzzz_yyzzzz_0, tg_zzzzzz_yzzzzz_0, tg_zzzzzz_zzzzzz_0, wp_y, \
                                         wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_y[j]; 

                    tg_yyyzzz_xxzzzz_0[j] = pb_y * tg_yyzzz_xxzzzz_0[j] + fr * tg_yyzzz_xxzzzz_1[j] + fl1_fx * (tg_yzzz_xxzzzz_0[j] - tg_yzzz_xxzzzz_1[j] * fl1_fza);

                    tg_yyyzzz_xyyyyy_0[j] = pb_y * tg_yyzzz_xyyyyy_0[j] + fr * tg_yyzzz_xyyyyy_1[j] + fl1_fx * (tg_yzzz_xyyyyy_0[j] - tg_yzzz_xyyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyzzz_xyyyy_1[j];

                    tg_yyyzzz_xyyyyz_0[j] = pb_y * tg_yyzzz_xyyyyz_0[j] + fr * tg_yyzzz_xyyyyz_1[j] + fl1_fx * (tg_yzzz_xyyyyz_0[j] - tg_yzzz_xyyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyzzz_xyyyz_1[j];

                    tg_yyyzzz_xyyyzz_0[j] = pb_y * tg_yyzzz_xyyyzz_0[j] + fr * tg_yyzzz_xyyyzz_1[j] + fl1_fx * (tg_yzzz_xyyyzz_0[j] - tg_yzzz_xyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzzz_xyyzz_1[j];

                    tg_yyyzzz_xyyzzz_0[j] = pb_y * tg_yyzzz_xyyzzz_0[j] + fr * tg_yyzzz_xyyzzz_1[j] + fl1_fx * (tg_yzzz_xyyzzz_0[j] - tg_yzzz_xyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyzzz_xyzzz_1[j];

                    tg_yyyzzz_xyzzzz_0[j] = pb_y * tg_yyzzz_xyzzzz_0[j] + fr * tg_yyzzz_xyzzzz_1[j] + fl1_fx * (tg_yzzz_xyzzzz_0[j] - tg_yzzz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzz_xzzzz_1[j];

                    tg_yyyzzz_xzzzzz_0[j] = pb_y * tg_yyzzz_xzzzzz_0[j] + fr * tg_yyzzz_xzzzzz_1[j] + fl1_fx * (tg_yzzz_xzzzzz_0[j] - tg_yzzz_xzzzzz_1[j] * fl1_fza);

                    tg_yyyzzz_yyyyyy_0[j] = pb_y * tg_yyzzz_yyyyyy_0[j] + fr * tg_yyzzz_yyyyyy_1[j] + fl1_fx * (tg_yzzz_yyyyyy_0[j] - tg_yzzz_yyyyyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yyzzz_yyyyy_1[j];

                    tg_yyyzzz_yyyyyz_0[j] = pb_y * tg_yyzzz_yyyyyz_0[j] + fr * tg_yyzzz_yyyyyz_1[j] + fl1_fx * (tg_yzzz_yyyyyz_0[j] - tg_yzzz_yyyyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyzzz_yyyyz_1[j];

                    tg_yyyzzz_yyyyzz_0[j] = pb_y * tg_yyzzz_yyyyzz_0[j] + fr * tg_yyzzz_yyyyzz_1[j] + fl1_fx * (tg_yzzz_yyyyzz_0[j] - tg_yzzz_yyyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyzzz_yyyzz_1[j];

                    tg_yyyzzz_yyyzzz_0[j] = pb_y * tg_yyzzz_yyyzzz_0[j] + fr * tg_yyzzz_yyyzzz_1[j] + fl1_fx * (tg_yzzz_yyyzzz_0[j] - tg_yzzz_yyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzzz_yyzzz_1[j];

                    tg_yyyzzz_yyzzzz_0[j] = pb_y * tg_yyzzz_yyzzzz_0[j] + fr * tg_yyzzz_yyzzzz_1[j] + fl1_fx * (tg_yzzz_yyzzzz_0[j] - tg_yzzz_yyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyzzz_yzzzz_1[j];

                    tg_yyyzzz_yzzzzz_0[j] = pb_y * tg_yyzzz_yzzzzz_0[j] + fr * tg_yyzzz_yzzzzz_1[j] + fl1_fx * (tg_yzzz_yzzzzz_0[j] - tg_yzzz_yzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzz_zzzzz_1[j];

                    tg_yyyzzz_zzzzzz_0[j] = pb_y * tg_yyzzz_zzzzzz_0[j] + fr * tg_yyzzz_zzzzzz_1[j] + fl1_fx * (tg_yzzz_zzzzzz_0[j] - tg_yzzz_zzzzzz_1[j] * fl1_fza);

                    tg_yyzzzz_xxxxxx_0[j] = pb_y * tg_yzzzz_xxxxxx_0[j] + fr * tg_yzzzz_xxxxxx_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxxxx_0[j] - tg_zzzz_xxxxxx_1[j] * fl1_fza);

                    tg_yyzzzz_xxxxxy_0[j] = pb_y * tg_yzzzz_xxxxxy_0[j] + fr * tg_yzzzz_xxxxxy_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxxxy_0[j] - tg_zzzz_xxxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzz_xxxxx_1[j];

                    tg_yyzzzz_xxxxxz_0[j] = pb_y * tg_yzzzz_xxxxxz_0[j] + fr * tg_yzzzz_xxxxxz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxxxz_0[j] - tg_zzzz_xxxxxz_1[j] * fl1_fza);

                    tg_yyzzzz_xxxxyy_0[j] = pb_y * tg_yzzzz_xxxxyy_0[j] + fr * tg_yzzzz_xxxxyy_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxxyy_0[j] - tg_zzzz_xxxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yzzzz_xxxxy_1[j];

                    tg_yyzzzz_xxxxyz_0[j] = pb_y * tg_yzzzz_xxxxyz_0[j] + fr * tg_yzzzz_xxxxyz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxxyz_0[j] - tg_zzzz_xxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzz_xxxxz_1[j];

                    tg_yyzzzz_xxxxzz_0[j] = pb_y * tg_yzzzz_xxxxzz_0[j] + fr * tg_yzzzz_xxxxzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxxzz_0[j] - tg_zzzz_xxxxzz_1[j] * fl1_fza);

                    tg_yyzzzz_xxxyyy_0[j] = pb_y * tg_yzzzz_xxxyyy_0[j] + fr * tg_yzzzz_xxxyyy_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxyyy_0[j] - tg_zzzz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzzz_xxxyy_1[j];

                    tg_yyzzzz_xxxyyz_0[j] = pb_y * tg_yzzzz_xxxyyz_0[j] + fr * tg_yzzzz_xxxyyz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxyyz_0[j] - tg_zzzz_xxxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yzzzz_xxxyz_1[j];

                    tg_yyzzzz_xxxyzz_0[j] = pb_y * tg_yzzzz_xxxyzz_0[j] + fr * tg_yzzzz_xxxyzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxyzz_0[j] - tg_zzzz_xxxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzz_xxxzz_1[j];

                    tg_yyzzzz_xxxzzz_0[j] = pb_y * tg_yzzzz_xxxzzz_0[j] + fr * tg_yzzzz_xxxzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxzzz_0[j] - tg_zzzz_xxxzzz_1[j] * fl1_fza);

                    tg_yyzzzz_xxyyyy_0[j] = pb_y * tg_yzzzz_xxyyyy_0[j] + fr * tg_yzzzz_xxyyyy_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxyyyy_0[j] - tg_zzzz_xxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzzzz_xxyyy_1[j];

                    tg_yyzzzz_xxyyyz_0[j] = pb_y * tg_yzzzz_xxyyyz_0[j] + fr * tg_yzzzz_xxyyyz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxyyyz_0[j] - tg_zzzz_xxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzzz_xxyyz_1[j];

                    tg_yyzzzz_xxyyzz_0[j] = pb_y * tg_yzzzz_xxyyzz_0[j] + fr * tg_yzzzz_xxyyzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxyyzz_0[j] - tg_zzzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yzzzz_xxyzz_1[j];

                    tg_yyzzzz_xxyzzz_0[j] = pb_y * tg_yzzzz_xxyzzz_0[j] + fr * tg_yzzzz_xxyzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxyzzz_0[j] - tg_zzzz_xxyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzz_xxzzz_1[j];

                    tg_yyzzzz_xxzzzz_0[j] = pb_y * tg_yzzzz_xxzzzz_0[j] + fr * tg_yzzzz_xxzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxzzzz_0[j] - tg_zzzz_xxzzzz_1[j] * fl1_fza);

                    tg_yyzzzz_xyyyyy_0[j] = pb_y * tg_yzzzz_xyyyyy_0[j] + fr * tg_yzzzz_xyyyyy_1[j] + 0.5 * fl1_fx * (tg_zzzz_xyyyyy_0[j] - tg_zzzz_xyyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yzzzz_xyyyy_1[j];

                    tg_yyzzzz_xyyyyz_0[j] = pb_y * tg_yzzzz_xyyyyz_0[j] + fr * tg_yzzzz_xyyyyz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xyyyyz_0[j] - tg_zzzz_xyyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzzzz_xyyyz_1[j];

                    tg_yyzzzz_xyyyzz_0[j] = pb_y * tg_yzzzz_xyyyzz_0[j] + fr * tg_yzzzz_xyyyzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xyyyzz_0[j] - tg_zzzz_xyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzzz_xyyzz_1[j];

                    tg_yyzzzz_xyyzzz_0[j] = pb_y * tg_yzzzz_xyyzzz_0[j] + fr * tg_yzzzz_xyyzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xyyzzz_0[j] - tg_zzzz_xyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yzzzz_xyzzz_1[j];

                    tg_yyzzzz_xyzzzz_0[j] = pb_y * tg_yzzzz_xyzzzz_0[j] + fr * tg_yzzzz_xyzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xyzzzz_0[j] - tg_zzzz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzz_xzzzz_1[j];

                    tg_yyzzzz_xzzzzz_0[j] = pb_y * tg_yzzzz_xzzzzz_0[j] + fr * tg_yzzzz_xzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xzzzzz_0[j] - tg_zzzz_xzzzzz_1[j] * fl1_fza);

                    tg_yyzzzz_yyyyyy_0[j] = pb_y * tg_yzzzz_yyyyyy_0[j] + fr * tg_yzzzz_yyyyyy_1[j] + 0.5 * fl1_fx * (tg_zzzz_yyyyyy_0[j] - tg_zzzz_yyyyyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yzzzz_yyyyy_1[j];

                    tg_yyzzzz_yyyyyz_0[j] = pb_y * tg_yzzzz_yyyyyz_0[j] + fr * tg_yzzzz_yyyyyz_1[j] + 0.5 * fl1_fx * (tg_zzzz_yyyyyz_0[j] - tg_zzzz_yyyyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yzzzz_yyyyz_1[j];

                    tg_yyzzzz_yyyyzz_0[j] = pb_y * tg_yzzzz_yyyyzz_0[j] + fr * tg_yzzzz_yyyyzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_yyyyzz_0[j] - tg_zzzz_yyyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzzzz_yyyzz_1[j];

                    tg_yyzzzz_yyyzzz_0[j] = pb_y * tg_yzzzz_yyyzzz_0[j] + fr * tg_yzzzz_yyyzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_yyyzzz_0[j] - tg_zzzz_yyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzzz_yyzzz_1[j];

                    tg_yyzzzz_yyzzzz_0[j] = pb_y * tg_yzzzz_yyzzzz_0[j] + fr * tg_yzzzz_yyzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_yyzzzz_0[j] - tg_zzzz_yyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yzzzz_yzzzz_1[j];

                    tg_yyzzzz_yzzzzz_0[j] = pb_y * tg_yzzzz_yzzzzz_0[j] + fr * tg_yzzzz_yzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_yzzzzz_0[j] - tg_zzzz_yzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzz_zzzzz_1[j];

                    tg_yyzzzz_zzzzzz_0[j] = pb_y * tg_yzzzz_zzzzzz_0[j] + fr * tg_yzzzz_zzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_zzzzzz_0[j] - tg_zzzz_zzzzzz_1[j] * fl1_fza);

                    tg_yzzzzz_xxxxxx_0[j] = pb_y * tg_zzzzz_xxxxxx_0[j] + fr * tg_zzzzz_xxxxxx_1[j];

                    tg_yzzzzz_xxxxxy_0[j] = pb_y * tg_zzzzz_xxxxxy_0[j] + fr * tg_zzzzz_xxxxxy_1[j] + 0.5 * fl1_fxn * tg_zzzzz_xxxxx_1[j];

                    tg_yzzzzz_xxxxxz_0[j] = pb_y * tg_zzzzz_xxxxxz_0[j] + fr * tg_zzzzz_xxxxxz_1[j];

                    tg_yzzzzz_xxxxyy_0[j] = pb_y * tg_zzzzz_xxxxyy_0[j] + fr * tg_zzzzz_xxxxyy_1[j] + fl1_fxn * tg_zzzzz_xxxxy_1[j];

                    tg_yzzzzz_xxxxyz_0[j] = pb_y * tg_zzzzz_xxxxyz_0[j] + fr * tg_zzzzz_xxxxyz_1[j] + 0.5 * fl1_fxn * tg_zzzzz_xxxxz_1[j];

                    tg_yzzzzz_xxxxzz_0[j] = pb_y * tg_zzzzz_xxxxzz_0[j] + fr * tg_zzzzz_xxxxzz_1[j];

                    tg_yzzzzz_xxxyyy_0[j] = pb_y * tg_zzzzz_xxxyyy_0[j] + fr * tg_zzzzz_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_zzzzz_xxxyy_1[j];

                    tg_yzzzzz_xxxyyz_0[j] = pb_y * tg_zzzzz_xxxyyz_0[j] + fr * tg_zzzzz_xxxyyz_1[j] + fl1_fxn * tg_zzzzz_xxxyz_1[j];

                    tg_yzzzzz_xxxyzz_0[j] = pb_y * tg_zzzzz_xxxyzz_0[j] + fr * tg_zzzzz_xxxyzz_1[j] + 0.5 * fl1_fxn * tg_zzzzz_xxxzz_1[j];

                    tg_yzzzzz_xxxzzz_0[j] = pb_y * tg_zzzzz_xxxzzz_0[j] + fr * tg_zzzzz_xxxzzz_1[j];

                    tg_yzzzzz_xxyyyy_0[j] = pb_y * tg_zzzzz_xxyyyy_0[j] + fr * tg_zzzzz_xxyyyy_1[j] + 2.0 * fl1_fxn * tg_zzzzz_xxyyy_1[j];

                    tg_yzzzzz_xxyyyz_0[j] = pb_y * tg_zzzzz_xxyyyz_0[j] + fr * tg_zzzzz_xxyyyz_1[j] + 1.5 * fl1_fxn * tg_zzzzz_xxyyz_1[j];

                    tg_yzzzzz_xxyyzz_0[j] = pb_y * tg_zzzzz_xxyyzz_0[j] + fr * tg_zzzzz_xxyyzz_1[j] + fl1_fxn * tg_zzzzz_xxyzz_1[j];

                    tg_yzzzzz_xxyzzz_0[j] = pb_y * tg_zzzzz_xxyzzz_0[j] + fr * tg_zzzzz_xxyzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzz_xxzzz_1[j];

                    tg_yzzzzz_xxzzzz_0[j] = pb_y * tg_zzzzz_xxzzzz_0[j] + fr * tg_zzzzz_xxzzzz_1[j];

                    tg_yzzzzz_xyyyyy_0[j] = pb_y * tg_zzzzz_xyyyyy_0[j] + fr * tg_zzzzz_xyyyyy_1[j] + 2.5 * fl1_fxn * tg_zzzzz_xyyyy_1[j];

                    tg_yzzzzz_xyyyyz_0[j] = pb_y * tg_zzzzz_xyyyyz_0[j] + fr * tg_zzzzz_xyyyyz_1[j] + 2.0 * fl1_fxn * tg_zzzzz_xyyyz_1[j];

                    tg_yzzzzz_xyyyzz_0[j] = pb_y * tg_zzzzz_xyyyzz_0[j] + fr * tg_zzzzz_xyyyzz_1[j] + 1.5 * fl1_fxn * tg_zzzzz_xyyzz_1[j];

                    tg_yzzzzz_xyyzzz_0[j] = pb_y * tg_zzzzz_xyyzzz_0[j] + fr * tg_zzzzz_xyyzzz_1[j] + fl1_fxn * tg_zzzzz_xyzzz_1[j];

                    tg_yzzzzz_xyzzzz_0[j] = pb_y * tg_zzzzz_xyzzzz_0[j] + fr * tg_zzzzz_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzz_xzzzz_1[j];

                    tg_yzzzzz_xzzzzz_0[j] = pb_y * tg_zzzzz_xzzzzz_0[j] + fr * tg_zzzzz_xzzzzz_1[j];

                    tg_yzzzzz_yyyyyy_0[j] = pb_y * tg_zzzzz_yyyyyy_0[j] + fr * tg_zzzzz_yyyyyy_1[j] + 3.0 * fl1_fxn * tg_zzzzz_yyyyy_1[j];

                    tg_yzzzzz_yyyyyz_0[j] = pb_y * tg_zzzzz_yyyyyz_0[j] + fr * tg_zzzzz_yyyyyz_1[j] + 2.5 * fl1_fxn * tg_zzzzz_yyyyz_1[j];

                    tg_yzzzzz_yyyyzz_0[j] = pb_y * tg_zzzzz_yyyyzz_0[j] + fr * tg_zzzzz_yyyyzz_1[j] + 2.0 * fl1_fxn * tg_zzzzz_yyyzz_1[j];

                    tg_yzzzzz_yyyzzz_0[j] = pb_y * tg_zzzzz_yyyzzz_0[j] + fr * tg_zzzzz_yyyzzz_1[j] + 1.5 * fl1_fxn * tg_zzzzz_yyzzz_1[j];

                    tg_yzzzzz_yyzzzz_0[j] = pb_y * tg_zzzzz_yyzzzz_0[j] + fr * tg_zzzzz_yyzzzz_1[j] + fl1_fxn * tg_zzzzz_yzzzz_1[j];

                    tg_yzzzzz_yzzzzz_0[j] = pb_y * tg_zzzzz_yzzzzz_0[j] + fr * tg_zzzzz_yzzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzz_zzzzz_1[j];

                    tg_yzzzzz_zzzzzz_0[j] = pb_y * tg_zzzzz_zzzzzz_0[j] + fr * tg_zzzzz_zzzzzz_1[j];

                    fr = wp_z[j]; 

                    tg_zzzzzz_xxxxxx_0[j] = pb_z * tg_zzzzz_xxxxxx_0[j] + fr * tg_zzzzz_xxxxxx_1[j] + 2.5 * fl1_fx * (tg_zzzz_xxxxxx_0[j] - tg_zzzz_xxxxxx_1[j] * fl1_fza);

                    tg_zzzzzz_xxxxxy_0[j] = pb_z * tg_zzzzz_xxxxxy_0[j] + fr * tg_zzzzz_xxxxxy_1[j] + 2.5 * fl1_fx * (tg_zzzz_xxxxxy_0[j] - tg_zzzz_xxxxxy_1[j] * fl1_fza);

                    tg_zzzzzz_xxxxxz_0[j] = pb_z * tg_zzzzz_xxxxxz_0[j] + fr * tg_zzzzz_xxxxxz_1[j] + 2.5 * fl1_fx * (tg_zzzz_xxxxxz_0[j] - tg_zzzz_xxxxxz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzz_xxxxx_1[j];

                    tg_zzzzzz_xxxxyy_0[j] = pb_z * tg_zzzzz_xxxxyy_0[j] + fr * tg_zzzzz_xxxxyy_1[j] + 2.5 * fl1_fx * (tg_zzzz_xxxxyy_0[j] - tg_zzzz_xxxxyy_1[j] * fl1_fza);

                    tg_zzzzzz_xxxxyz_0[j] = pb_z * tg_zzzzz_xxxxyz_0[j] + fr * tg_zzzzz_xxxxyz_1[j] + 2.5 * fl1_fx * (tg_zzzz_xxxxyz_0[j] - tg_zzzz_xxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzz_xxxxy_1[j];

                    tg_zzzzzz_xxxxzz_0[j] = pb_z * tg_zzzzz_xxxxzz_0[j] + fr * tg_zzzzz_xxxxzz_1[j] + 2.5 * fl1_fx * (tg_zzzz_xxxxzz_0[j] - tg_zzzz_xxxxzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzzz_xxxxz_1[j];

                    tg_zzzzzz_xxxyyy_0[j] = pb_z * tg_zzzzz_xxxyyy_0[j] + fr * tg_zzzzz_xxxyyy_1[j] + 2.5 * fl1_fx * (tg_zzzz_xxxyyy_0[j] - tg_zzzz_xxxyyy_1[j] * fl1_fza);

                    tg_zzzzzz_xxxyyz_0[j] = pb_z * tg_zzzzz_xxxyyz_0[j] + fr * tg_zzzzz_xxxyyz_1[j] + 2.5 * fl1_fx * (tg_zzzz_xxxyyz_0[j] - tg_zzzz_xxxyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzz_xxxyy_1[j];

                    tg_zzzzzz_xxxyzz_0[j] = pb_z * tg_zzzzz_xxxyzz_0[j] + fr * tg_zzzzz_xxxyzz_1[j] + 2.5 * fl1_fx * (tg_zzzz_xxxyzz_0[j] - tg_zzzz_xxxyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzzz_xxxyz_1[j];

                    tg_zzzzzz_xxxzzz_0[j] = pb_z * tg_zzzzz_xxxzzz_0[j] + fr * tg_zzzzz_xxxzzz_1[j] + 2.5 * fl1_fx * (tg_zzzz_xxxzzz_0[j] - tg_zzzz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzzz_xxxzz_1[j];

                    tg_zzzzzz_xxyyyy_0[j] = pb_z * tg_zzzzz_xxyyyy_0[j] + fr * tg_zzzzz_xxyyyy_1[j] + 2.5 * fl1_fx * (tg_zzzz_xxyyyy_0[j] - tg_zzzz_xxyyyy_1[j] * fl1_fza);

                    tg_zzzzzz_xxyyyz_0[j] = pb_z * tg_zzzzz_xxyyyz_0[j] + fr * tg_zzzzz_xxyyyz_1[j] + 2.5 * fl1_fx * (tg_zzzz_xxyyyz_0[j] - tg_zzzz_xxyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzz_xxyyy_1[j];

                    tg_zzzzzz_xxyyzz_0[j] = pb_z * tg_zzzzz_xxyyzz_0[j] + fr * tg_zzzzz_xxyyzz_1[j] + 2.5 * fl1_fx * (tg_zzzz_xxyyzz_0[j] - tg_zzzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzzz_xxyyz_1[j];

                    tg_zzzzzz_xxyzzz_0[j] = pb_z * tg_zzzzz_xxyzzz_0[j] + fr * tg_zzzzz_xxyzzz_1[j] + 2.5 * fl1_fx * (tg_zzzz_xxyzzz_0[j] - tg_zzzz_xxyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzzz_xxyzz_1[j];

                    tg_zzzzzz_xxzzzz_0[j] = pb_z * tg_zzzzz_xxzzzz_0[j] + fr * tg_zzzzz_xxzzzz_1[j] + 2.5 * fl1_fx * (tg_zzzz_xxzzzz_0[j] - tg_zzzz_xxzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzzzz_xxzzz_1[j];

                    tg_zzzzzz_xyyyyy_0[j] = pb_z * tg_zzzzz_xyyyyy_0[j] + fr * tg_zzzzz_xyyyyy_1[j] + 2.5 * fl1_fx * (tg_zzzz_xyyyyy_0[j] - tg_zzzz_xyyyyy_1[j] * fl1_fza);

                    tg_zzzzzz_xyyyyz_0[j] = pb_z * tg_zzzzz_xyyyyz_0[j] + fr * tg_zzzzz_xyyyyz_1[j] + 2.5 * fl1_fx * (tg_zzzz_xyyyyz_0[j] - tg_zzzz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzz_xyyyy_1[j];

                    tg_zzzzzz_xyyyzz_0[j] = pb_z * tg_zzzzz_xyyyzz_0[j] + fr * tg_zzzzz_xyyyzz_1[j] + 2.5 * fl1_fx * (tg_zzzz_xyyyzz_0[j] - tg_zzzz_xyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzzz_xyyyz_1[j];

                    tg_zzzzzz_xyyzzz_0[j] = pb_z * tg_zzzzz_xyyzzz_0[j] + fr * tg_zzzzz_xyyzzz_1[j] + 2.5 * fl1_fx * (tg_zzzz_xyyzzz_0[j] - tg_zzzz_xyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzzz_xyyzz_1[j];

                    tg_zzzzzz_xyzzzz_0[j] = pb_z * tg_zzzzz_xyzzzz_0[j] + fr * tg_zzzzz_xyzzzz_1[j] + 2.5 * fl1_fx * (tg_zzzz_xyzzzz_0[j] - tg_zzzz_xyzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzzzz_xyzzz_1[j];

                    tg_zzzzzz_xzzzzz_0[j] = pb_z * tg_zzzzz_xzzzzz_0[j] + fr * tg_zzzzz_xzzzzz_1[j] + 2.5 * fl1_fx * (tg_zzzz_xzzzzz_0[j] - tg_zzzz_xzzzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_zzzzz_xzzzz_1[j];

                    tg_zzzzzz_yyyyyy_0[j] = pb_z * tg_zzzzz_yyyyyy_0[j] + fr * tg_zzzzz_yyyyyy_1[j] + 2.5 * fl1_fx * (tg_zzzz_yyyyyy_0[j] - tg_zzzz_yyyyyy_1[j] * fl1_fza);

                    tg_zzzzzz_yyyyyz_0[j] = pb_z * tg_zzzzz_yyyyyz_0[j] + fr * tg_zzzzz_yyyyyz_1[j] + 2.5 * fl1_fx * (tg_zzzz_yyyyyz_0[j] - tg_zzzz_yyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzz_yyyyy_1[j];

                    tg_zzzzzz_yyyyzz_0[j] = pb_z * tg_zzzzz_yyyyzz_0[j] + fr * tg_zzzzz_yyyyzz_1[j] + 2.5 * fl1_fx * (tg_zzzz_yyyyzz_0[j] - tg_zzzz_yyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzzz_yyyyz_1[j];

                    tg_zzzzzz_yyyzzz_0[j] = pb_z * tg_zzzzz_yyyzzz_0[j] + fr * tg_zzzzz_yyyzzz_1[j] + 2.5 * fl1_fx * (tg_zzzz_yyyzzz_0[j] - tg_zzzz_yyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzzz_yyyzz_1[j];

                    tg_zzzzzz_yyzzzz_0[j] = pb_z * tg_zzzzz_yyzzzz_0[j] + fr * tg_zzzzz_yyzzzz_1[j] + 2.5 * fl1_fx * (tg_zzzz_yyzzzz_0[j] - tg_zzzz_yyzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzzzz_yyzzz_1[j];

                    tg_zzzzzz_yzzzzz_0[j] = pb_z * tg_zzzzz_yzzzzz_0[j] + fr * tg_zzzzz_yzzzzz_1[j] + 2.5 * fl1_fx * (tg_zzzz_yzzzzz_0[j] - tg_zzzz_yzzzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_zzzzz_yzzzz_1[j];

                    tg_zzzzzz_zzzzzz_0[j] = pb_z * tg_zzzzz_zzzzzz_0[j] + fr * tg_zzzzz_zzzzzz_1[j] + 2.5 * fl1_fx * (tg_zzzz_zzzzzz_0[j] - tg_zzzz_zzzzzz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_zzzzz_zzzzz_1[j];
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

