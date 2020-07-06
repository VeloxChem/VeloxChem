//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ElectronRepulsionRecFuncForHL.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSHSL(      CMemBlock2D<double>* primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSHSL_0_95(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSHSL_95_190(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 

        erirecfunc::compElectronRepulsionForSHSL_190_285(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSHSL_285_380(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSHSL_380_475(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSHSL_475_569(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSHSL_569_663(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSHSL_663_757(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSHSL_757_851(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSHSL_851_945(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 
    }

    void
    compElectronRepulsionForSHSL_0_95(      CMemBlock2D<double>* primBuffer,
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
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_8_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_8_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_8_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_8_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_8_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_7_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {7, -1, -1, -1}, 
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

                auto tg_xxxx_xxxxxxxx_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx); 

                auto tg_xxxx_xxxxxxxy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 1); 

                auto tg_xxxx_xxxxxxxz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 2); 

                auto tg_xxxx_xxxxxxyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 3); 

                auto tg_xxxx_xxxxxxyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 4); 

                auto tg_xxxx_xxxxxxzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 5); 

                auto tg_xxxx_xxxxxyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 6); 

                auto tg_xxxx_xxxxxyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 7); 

                auto tg_xxxx_xxxxxyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 8); 

                auto tg_xxxx_xxxxxzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 9); 

                auto tg_xxxx_xxxxyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 10); 

                auto tg_xxxx_xxxxyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 11); 

                auto tg_xxxx_xxxxyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 12); 

                auto tg_xxxx_xxxxyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 13); 

                auto tg_xxxx_xxxxzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 14); 

                auto tg_xxxx_xxxyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 15); 

                auto tg_xxxx_xxxyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 16); 

                auto tg_xxxx_xxxyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 17); 

                auto tg_xxxx_xxxyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 18); 

                auto tg_xxxx_xxxyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 19); 

                auto tg_xxxx_xxxzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 20); 

                auto tg_xxxx_xxyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 21); 

                auto tg_xxxx_xxyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 22); 

                auto tg_xxxx_xxyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 23); 

                auto tg_xxxx_xxyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 24); 

                auto tg_xxxx_xxyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 25); 

                auto tg_xxxx_xxyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 26); 

                auto tg_xxxx_xxzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 27); 

                auto tg_xxxx_xyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 28); 

                auto tg_xxxx_xyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 29); 

                auto tg_xxxx_xyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 30); 

                auto tg_xxxx_xyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 31); 

                auto tg_xxxx_xyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 32); 

                auto tg_xxxx_xyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 33); 

                auto tg_xxxx_xyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 34); 

                auto tg_xxxx_xzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 35); 

                auto tg_xxxx_yyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 36); 

                auto tg_xxxx_yyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 37); 

                auto tg_xxxx_yyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 38); 

                auto tg_xxxx_yyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 39); 

                auto tg_xxxx_yyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 40); 

                auto tg_xxxx_yyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 41); 

                auto tg_xxxx_yyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 42); 

                auto tg_xxxx_yzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 43); 

                auto tg_xxxx_zzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 44); 

                auto tg_xxxy_xxxxxxxx_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 45); 

                auto tg_xxxy_xxxxxxxy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 46); 

                auto tg_xxxy_xxxxxxxz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 47); 

                auto tg_xxxy_xxxxxxyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 48); 

                auto tg_xxxy_xxxxxxyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 49); 

                auto tg_xxxy_xxxxxxzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 50); 

                auto tg_xxxy_xxxxxyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 51); 

                auto tg_xxxy_xxxxxyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 52); 

                auto tg_xxxy_xxxxxyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 53); 

                auto tg_xxxy_xxxxxzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 54); 

                auto tg_xxxy_xxxxyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 55); 

                auto tg_xxxy_xxxxyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 56); 

                auto tg_xxxy_xxxxyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 57); 

                auto tg_xxxy_xxxxyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 58); 

                auto tg_xxxy_xxxxzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 59); 

                auto tg_xxxy_xxxyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 60); 

                auto tg_xxxy_xxxyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 61); 

                auto tg_xxxy_xxxyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 62); 

                auto tg_xxxy_xxxyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 63); 

                auto tg_xxxy_xxxyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 64); 

                auto tg_xxxy_xxxzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 65); 

                auto tg_xxxy_xxyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 66); 

                auto tg_xxxy_xxyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 67); 

                auto tg_xxxy_xxyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 68); 

                auto tg_xxxy_xxyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 69); 

                auto tg_xxxy_xxyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 70); 

                auto tg_xxxy_xxyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 71); 

                auto tg_xxxy_xxzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 72); 

                auto tg_xxxy_xyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 73); 

                auto tg_xxxy_xyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 74); 

                auto tg_xxxy_xyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 75); 

                auto tg_xxxy_xyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 76); 

                auto tg_xxxy_xyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 77); 

                auto tg_xxxy_xyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 78); 

                auto tg_xxxy_xyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 79); 

                auto tg_xxxy_xzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 80); 

                auto tg_xxxy_yyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 81); 

                auto tg_xxxy_yyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 82); 

                auto tg_xxxy_yyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 83); 

                auto tg_xxxy_yyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 84); 

                auto tg_xxxy_yyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 85); 

                auto tg_xxxy_yyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 86); 

                auto tg_xxxy_yyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 87); 

                auto tg_xxxy_yzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 88); 

                auto tg_xxxy_zzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 89); 

                auto tg_xxxz_xxxxxxxx_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 90); 

                auto tg_xxxz_xxxxxxxy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 91); 

                auto tg_xxxz_xxxxxxxz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 92); 

                auto tg_xxxz_xxxxxxyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 93); 

                auto tg_xxxz_xxxxxxyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 94); 

                auto tg_xxxx_xxxxxxxx_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx); 

                auto tg_xxxx_xxxxxxxy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 1); 

                auto tg_xxxx_xxxxxxxz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 2); 

                auto tg_xxxx_xxxxxxyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 3); 

                auto tg_xxxx_xxxxxxyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 4); 

                auto tg_xxxx_xxxxxxzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 5); 

                auto tg_xxxx_xxxxxyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 6); 

                auto tg_xxxx_xxxxxyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 7); 

                auto tg_xxxx_xxxxxyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 8); 

                auto tg_xxxx_xxxxxzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 9); 

                auto tg_xxxx_xxxxyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 10); 

                auto tg_xxxx_xxxxyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 11); 

                auto tg_xxxx_xxxxyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 12); 

                auto tg_xxxx_xxxxyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 13); 

                auto tg_xxxx_xxxxzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 14); 

                auto tg_xxxx_xxxyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 15); 

                auto tg_xxxx_xxxyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 16); 

                auto tg_xxxx_xxxyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 17); 

                auto tg_xxxx_xxxyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 18); 

                auto tg_xxxx_xxxyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 19); 

                auto tg_xxxx_xxxzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 20); 

                auto tg_xxxx_xxyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 21); 

                auto tg_xxxx_xxyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 22); 

                auto tg_xxxx_xxyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 23); 

                auto tg_xxxx_xxyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 24); 

                auto tg_xxxx_xxyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 25); 

                auto tg_xxxx_xxyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 26); 

                auto tg_xxxx_xxzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 27); 

                auto tg_xxxx_xyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 28); 

                auto tg_xxxx_xyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 29); 

                auto tg_xxxx_xyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 30); 

                auto tg_xxxx_xyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 31); 

                auto tg_xxxx_xyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 32); 

                auto tg_xxxx_xyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 33); 

                auto tg_xxxx_xyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 34); 

                auto tg_xxxx_xzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 35); 

                auto tg_xxxx_yyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 36); 

                auto tg_xxxx_yyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 37); 

                auto tg_xxxx_yyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 38); 

                auto tg_xxxx_yyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 39); 

                auto tg_xxxx_yyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 40); 

                auto tg_xxxx_yyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 41); 

                auto tg_xxxx_yyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 42); 

                auto tg_xxxx_yzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 43); 

                auto tg_xxxx_zzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 44); 

                auto tg_xxxy_xxxxxxxx_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 45); 

                auto tg_xxxy_xxxxxxxy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 46); 

                auto tg_xxxy_xxxxxxxz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 47); 

                auto tg_xxxy_xxxxxxyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 48); 

                auto tg_xxxy_xxxxxxyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 49); 

                auto tg_xxxy_xxxxxxzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 50); 

                auto tg_xxxy_xxxxxyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 51); 

                auto tg_xxxy_xxxxxyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 52); 

                auto tg_xxxy_xxxxxyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 53); 

                auto tg_xxxy_xxxxxzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 54); 

                auto tg_xxxy_xxxxyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 55); 

                auto tg_xxxy_xxxxyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 56); 

                auto tg_xxxy_xxxxyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 57); 

                auto tg_xxxy_xxxxyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 58); 

                auto tg_xxxy_xxxxzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 59); 

                auto tg_xxxy_xxxyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 60); 

                auto tg_xxxy_xxxyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 61); 

                auto tg_xxxy_xxxyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 62); 

                auto tg_xxxy_xxxyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 63); 

                auto tg_xxxy_xxxyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 64); 

                auto tg_xxxy_xxxzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 65); 

                auto tg_xxxy_xxyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 66); 

                auto tg_xxxy_xxyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 67); 

                auto tg_xxxy_xxyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 68); 

                auto tg_xxxy_xxyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 69); 

                auto tg_xxxy_xxyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 70); 

                auto tg_xxxy_xxyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 71); 

                auto tg_xxxy_xxzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 72); 

                auto tg_xxxy_xyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 73); 

                auto tg_xxxy_xyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 74); 

                auto tg_xxxy_xyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 75); 

                auto tg_xxxy_xyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 76); 

                auto tg_xxxy_xyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 77); 

                auto tg_xxxy_xyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 78); 

                auto tg_xxxy_xyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 79); 

                auto tg_xxxy_xzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 80); 

                auto tg_xxxy_yyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 81); 

                auto tg_xxxy_yyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 82); 

                auto tg_xxxy_yyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 83); 

                auto tg_xxxy_yyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 84); 

                auto tg_xxxy_yyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 85); 

                auto tg_xxxy_yyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 86); 

                auto tg_xxxy_yyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 87); 

                auto tg_xxxy_yzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 88); 

                auto tg_xxxy_zzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 89); 

                auto tg_xxxz_xxxxxxxx_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 90); 

                auto tg_xxxz_xxxxxxxy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 91); 

                auto tg_xxxz_xxxxxxxz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 92); 

                auto tg_xxxz_xxxxxxyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 93); 

                auto tg_xxxz_xxxxxxyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 94); 

                auto tg_xxx_xxxxxxxx_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx); 

                auto tg_xxx_xxxxxxxy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 1); 

                auto tg_xxx_xxxxxxxz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 2); 

                auto tg_xxx_xxxxxxyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 3); 

                auto tg_xxx_xxxxxxyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 4); 

                auto tg_xxx_xxxxxxzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 5); 

                auto tg_xxx_xxxxxyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 6); 

                auto tg_xxx_xxxxxyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 7); 

                auto tg_xxx_xxxxxyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 8); 

                auto tg_xxx_xxxxxzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 9); 

                auto tg_xxx_xxxxyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 10); 

                auto tg_xxx_xxxxyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 11); 

                auto tg_xxx_xxxxyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 12); 

                auto tg_xxx_xxxxyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 13); 

                auto tg_xxx_xxxxzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 14); 

                auto tg_xxx_xxxyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 15); 

                auto tg_xxx_xxxyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 16); 

                auto tg_xxx_xxxyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 17); 

                auto tg_xxx_xxxyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 18); 

                auto tg_xxx_xxxyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 19); 

                auto tg_xxx_xxxzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 20); 

                auto tg_xxx_xxyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 21); 

                auto tg_xxx_xxyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 22); 

                auto tg_xxx_xxyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 23); 

                auto tg_xxx_xxyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 24); 

                auto tg_xxx_xxyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 25); 

                auto tg_xxx_xxyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 26); 

                auto tg_xxx_xxzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 27); 

                auto tg_xxx_xyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 28); 

                auto tg_xxx_xyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 29); 

                auto tg_xxx_xyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 30); 

                auto tg_xxx_xyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 31); 

                auto tg_xxx_xyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 32); 

                auto tg_xxx_xyyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 33); 

                auto tg_xxx_xyzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 34); 

                auto tg_xxx_xzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 35); 

                auto tg_xxx_yyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 36); 

                auto tg_xxx_yyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 37); 

                auto tg_xxx_yyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 38); 

                auto tg_xxx_yyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 39); 

                auto tg_xxx_yyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 40); 

                auto tg_xxx_yyyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 41); 

                auto tg_xxx_yyzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 42); 

                auto tg_xxx_yzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 43); 

                auto tg_xxx_zzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 44); 

                auto tg_xxy_xxxxxxxx_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 45); 

                auto tg_xxy_xxxxxxxy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 46); 

                auto tg_xxy_xxxxxxxz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 47); 

                auto tg_xxy_xxxxxxyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 48); 

                auto tg_xxy_xxxxxxyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 49); 

                auto tg_xxy_xxxxxxzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 50); 

                auto tg_xxy_xxxxxyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 51); 

                auto tg_xxy_xxxxxyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 52); 

                auto tg_xxy_xxxxxyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 53); 

                auto tg_xxy_xxxxxzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 54); 

                auto tg_xxy_xxxxyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 55); 

                auto tg_xxy_xxxxyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 56); 

                auto tg_xxy_xxxxyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 57); 

                auto tg_xxy_xxxxyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 58); 

                auto tg_xxy_xxxxzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 59); 

                auto tg_xxy_xxxyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 60); 

                auto tg_xxy_xxxyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 61); 

                auto tg_xxy_xxxyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 62); 

                auto tg_xxy_xxxyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 63); 

                auto tg_xxy_xxxyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 64); 

                auto tg_xxy_xxxzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 65); 

                auto tg_xxy_xxyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 66); 

                auto tg_xxy_xxyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 67); 

                auto tg_xxy_xxyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 68); 

                auto tg_xxy_xxyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 69); 

                auto tg_xxy_xxyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 70); 

                auto tg_xxy_xxyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 71); 

                auto tg_xxy_xxzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 72); 

                auto tg_xxy_xyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 73); 

                auto tg_xxy_xyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 74); 

                auto tg_xxy_xyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 75); 

                auto tg_xxy_xyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 76); 

                auto tg_xxy_xyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 77); 

                auto tg_xxy_xyyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 78); 

                auto tg_xxy_xyzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 79); 

                auto tg_xxy_xzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 80); 

                auto tg_xxy_yyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 81); 

                auto tg_xxy_yyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 82); 

                auto tg_xxy_yyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 83); 

                auto tg_xxy_yyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 84); 

                auto tg_xxy_yyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 85); 

                auto tg_xxy_yyyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 86); 

                auto tg_xxy_yyzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 87); 

                auto tg_xxy_yzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 88); 

                auto tg_xxy_zzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 89); 

                auto tg_xxz_xxxxxxxx_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 90); 

                auto tg_xxz_xxxxxxxy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 91); 

                auto tg_xxz_xxxxxxxz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 92); 

                auto tg_xxz_xxxxxxyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 93); 

                auto tg_xxz_xxxxxxyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 94); 

                auto tg_xxx_xxxxxxxx_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx); 

                auto tg_xxx_xxxxxxxy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 1); 

                auto tg_xxx_xxxxxxxz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 2); 

                auto tg_xxx_xxxxxxyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 3); 

                auto tg_xxx_xxxxxxyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 4); 

                auto tg_xxx_xxxxxxzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 5); 

                auto tg_xxx_xxxxxyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 6); 

                auto tg_xxx_xxxxxyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 7); 

                auto tg_xxx_xxxxxyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 8); 

                auto tg_xxx_xxxxxzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 9); 

                auto tg_xxx_xxxxyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 10); 

                auto tg_xxx_xxxxyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 11); 

                auto tg_xxx_xxxxyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 12); 

                auto tg_xxx_xxxxyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 13); 

                auto tg_xxx_xxxxzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 14); 

                auto tg_xxx_xxxyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 15); 

                auto tg_xxx_xxxyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 16); 

                auto tg_xxx_xxxyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 17); 

                auto tg_xxx_xxxyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 18); 

                auto tg_xxx_xxxyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 19); 

                auto tg_xxx_xxxzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 20); 

                auto tg_xxx_xxyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 21); 

                auto tg_xxx_xxyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 22); 

                auto tg_xxx_xxyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 23); 

                auto tg_xxx_xxyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 24); 

                auto tg_xxx_xxyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 25); 

                auto tg_xxx_xxyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 26); 

                auto tg_xxx_xxzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 27); 

                auto tg_xxx_xyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 28); 

                auto tg_xxx_xyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 29); 

                auto tg_xxx_xyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 30); 

                auto tg_xxx_xyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 31); 

                auto tg_xxx_xyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 32); 

                auto tg_xxx_xyyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 33); 

                auto tg_xxx_xyzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 34); 

                auto tg_xxx_xzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 35); 

                auto tg_xxx_yyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 36); 

                auto tg_xxx_yyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 37); 

                auto tg_xxx_yyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 38); 

                auto tg_xxx_yyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 39); 

                auto tg_xxx_yyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 40); 

                auto tg_xxx_yyyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 41); 

                auto tg_xxx_yyzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 42); 

                auto tg_xxx_yzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 43); 

                auto tg_xxx_zzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 44); 

                auto tg_xxy_xxxxxxxx_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 45); 

                auto tg_xxy_xxxxxxxy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 46); 

                auto tg_xxy_xxxxxxxz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 47); 

                auto tg_xxy_xxxxxxyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 48); 

                auto tg_xxy_xxxxxxyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 49); 

                auto tg_xxy_xxxxxxzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 50); 

                auto tg_xxy_xxxxxyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 51); 

                auto tg_xxy_xxxxxyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 52); 

                auto tg_xxy_xxxxxyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 53); 

                auto tg_xxy_xxxxxzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 54); 

                auto tg_xxy_xxxxyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 55); 

                auto tg_xxy_xxxxyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 56); 

                auto tg_xxy_xxxxyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 57); 

                auto tg_xxy_xxxxyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 58); 

                auto tg_xxy_xxxxzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 59); 

                auto tg_xxy_xxxyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 60); 

                auto tg_xxy_xxxyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 61); 

                auto tg_xxy_xxxyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 62); 

                auto tg_xxy_xxxyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 63); 

                auto tg_xxy_xxxyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 64); 

                auto tg_xxy_xxxzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 65); 

                auto tg_xxy_xxyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 66); 

                auto tg_xxy_xxyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 67); 

                auto tg_xxy_xxyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 68); 

                auto tg_xxy_xxyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 69); 

                auto tg_xxy_xxyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 70); 

                auto tg_xxy_xxyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 71); 

                auto tg_xxy_xxzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 72); 

                auto tg_xxy_xyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 73); 

                auto tg_xxy_xyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 74); 

                auto tg_xxy_xyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 75); 

                auto tg_xxy_xyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 76); 

                auto tg_xxy_xyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 77); 

                auto tg_xxy_xyyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 78); 

                auto tg_xxy_xyzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 79); 

                auto tg_xxy_xzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 80); 

                auto tg_xxy_yyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 81); 

                auto tg_xxy_yyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 82); 

                auto tg_xxy_yyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 83); 

                auto tg_xxy_yyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 84); 

                auto tg_xxy_yyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 85); 

                auto tg_xxy_yyyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 86); 

                auto tg_xxy_yyzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 87); 

                auto tg_xxy_yzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 88); 

                auto tg_xxy_zzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 89); 

                auto tg_xxz_xxxxxxxx_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 90); 

                auto tg_xxz_xxxxxxxy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 91); 

                auto tg_xxz_xxxxxxxz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 92); 

                auto tg_xxz_xxxxxxyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 93); 

                auto tg_xxz_xxxxxxyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 94); 

                auto tg_xxxx_xxxxxxx_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx); 

                auto tg_xxxx_xxxxxxy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 1); 

                auto tg_xxxx_xxxxxxz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 2); 

                auto tg_xxxx_xxxxxyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 3); 

                auto tg_xxxx_xxxxxyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 4); 

                auto tg_xxxx_xxxxxzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 5); 

                auto tg_xxxx_xxxxyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 6); 

                auto tg_xxxx_xxxxyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 7); 

                auto tg_xxxx_xxxxyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 8); 

                auto tg_xxxx_xxxxzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 9); 

                auto tg_xxxx_xxxyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 10); 

                auto tg_xxxx_xxxyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 11); 

                auto tg_xxxx_xxxyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 12); 

                auto tg_xxxx_xxxyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 13); 

                auto tg_xxxx_xxxzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 14); 

                auto tg_xxxx_xxyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 15); 

                auto tg_xxxx_xxyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 16); 

                auto tg_xxxx_xxyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 17); 

                auto tg_xxxx_xxyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 18); 

                auto tg_xxxx_xxyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 19); 

                auto tg_xxxx_xxzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 20); 

                auto tg_xxxx_xyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 21); 

                auto tg_xxxx_xyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 22); 

                auto tg_xxxx_xyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 23); 

                auto tg_xxxx_xyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 24); 

                auto tg_xxxx_xyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 25); 

                auto tg_xxxx_xyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 26); 

                auto tg_xxxx_xzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 27); 

                auto tg_xxxx_yyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 28); 

                auto tg_xxxx_yyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 29); 

                auto tg_xxxx_yyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 30); 

                auto tg_xxxx_yyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 31); 

                auto tg_xxxx_yyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 32); 

                auto tg_xxxx_yyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 33); 

                auto tg_xxxx_yzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 34); 

                auto tg_xxxx_zzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 35); 

                auto tg_xxxy_xxxxxxx_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 36); 

                auto tg_xxxy_xxxxxxy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 37); 

                auto tg_xxxy_xxxxxxz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 38); 

                auto tg_xxxy_xxxxxyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 39); 

                auto tg_xxxy_xxxxxyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 40); 

                auto tg_xxxy_xxxxxzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 41); 

                auto tg_xxxy_xxxxyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 42); 

                auto tg_xxxy_xxxxyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 43); 

                auto tg_xxxy_xxxxyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 44); 

                auto tg_xxxy_xxxxzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 45); 

                auto tg_xxxy_xxxyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 46); 

                auto tg_xxxy_xxxyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 47); 

                auto tg_xxxy_xxxyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 48); 

                auto tg_xxxy_xxxyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 49); 

                auto tg_xxxy_xxxzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 50); 

                auto tg_xxxy_xxyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 51); 

                auto tg_xxxy_xxyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 52); 

                auto tg_xxxy_xxyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 53); 

                auto tg_xxxy_xxyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 54); 

                auto tg_xxxy_xxyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 55); 

                auto tg_xxxy_xxzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 56); 

                auto tg_xxxy_xyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 57); 

                auto tg_xxxy_xyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 58); 

                auto tg_xxxy_xyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 59); 

                auto tg_xxxy_xyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 60); 

                auto tg_xxxy_xyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 61); 

                auto tg_xxxy_xyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 62); 

                auto tg_xxxy_xzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 63); 

                auto tg_xxxy_yyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 64); 

                auto tg_xxxy_yyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 65); 

                auto tg_xxxy_yyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 66); 

                auto tg_xxxy_yyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 67); 

                auto tg_xxxy_yyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 68); 

                auto tg_xxxy_yyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 69); 

                auto tg_xxxy_yzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 70); 

                auto tg_xxxy_zzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 71); 

                auto tg_xxxz_xxxxxxx_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 72); 

                auto tg_xxxz_xxxxxxy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 73); 

                auto tg_xxxz_xxxxxxz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 74); 

                auto tg_xxxz_xxxxxyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 75); 

                auto tg_xxxz_xxxxxyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 76); 

                // set up pointers to integrals

                auto tg_xxxxx_xxxxxxxx_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx); 

                auto tg_xxxxx_xxxxxxxy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 1); 

                auto tg_xxxxx_xxxxxxxz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 2); 

                auto tg_xxxxx_xxxxxxyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 3); 

                auto tg_xxxxx_xxxxxxyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 4); 

                auto tg_xxxxx_xxxxxxzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 5); 

                auto tg_xxxxx_xxxxxyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 6); 

                auto tg_xxxxx_xxxxxyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 7); 

                auto tg_xxxxx_xxxxxyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 8); 

                auto tg_xxxxx_xxxxxzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 9); 

                auto tg_xxxxx_xxxxyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 10); 

                auto tg_xxxxx_xxxxyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 11); 

                auto tg_xxxxx_xxxxyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 12); 

                auto tg_xxxxx_xxxxyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 13); 

                auto tg_xxxxx_xxxxzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 14); 

                auto tg_xxxxx_xxxyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 15); 

                auto tg_xxxxx_xxxyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 16); 

                auto tg_xxxxx_xxxyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 17); 

                auto tg_xxxxx_xxxyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 18); 

                auto tg_xxxxx_xxxyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 19); 

                auto tg_xxxxx_xxxzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 20); 

                auto tg_xxxxx_xxyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 21); 

                auto tg_xxxxx_xxyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 22); 

                auto tg_xxxxx_xxyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 23); 

                auto tg_xxxxx_xxyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 24); 

                auto tg_xxxxx_xxyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 25); 

                auto tg_xxxxx_xxyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 26); 

                auto tg_xxxxx_xxzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 27); 

                auto tg_xxxxx_xyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 28); 

                auto tg_xxxxx_xyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 29); 

                auto tg_xxxxx_xyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 30); 

                auto tg_xxxxx_xyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 31); 

                auto tg_xxxxx_xyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 32); 

                auto tg_xxxxx_xyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 33); 

                auto tg_xxxxx_xyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 34); 

                auto tg_xxxxx_xzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 35); 

                auto tg_xxxxx_yyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 36); 

                auto tg_xxxxx_yyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 37); 

                auto tg_xxxxx_yyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 38); 

                auto tg_xxxxx_yyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 39); 

                auto tg_xxxxx_yyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 40); 

                auto tg_xxxxx_yyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 41); 

                auto tg_xxxxx_yyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 42); 

                auto tg_xxxxx_yzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 43); 

                auto tg_xxxxx_zzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 44); 

                auto tg_xxxxy_xxxxxxxx_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 45); 

                auto tg_xxxxy_xxxxxxxy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 46); 

                auto tg_xxxxy_xxxxxxxz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 47); 

                auto tg_xxxxy_xxxxxxyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 48); 

                auto tg_xxxxy_xxxxxxyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 49); 

                auto tg_xxxxy_xxxxxxzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 50); 

                auto tg_xxxxy_xxxxxyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 51); 

                auto tg_xxxxy_xxxxxyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 52); 

                auto tg_xxxxy_xxxxxyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 53); 

                auto tg_xxxxy_xxxxxzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 54); 

                auto tg_xxxxy_xxxxyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 55); 

                auto tg_xxxxy_xxxxyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 56); 

                auto tg_xxxxy_xxxxyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 57); 

                auto tg_xxxxy_xxxxyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 58); 

                auto tg_xxxxy_xxxxzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 59); 

                auto tg_xxxxy_xxxyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 60); 

                auto tg_xxxxy_xxxyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 61); 

                auto tg_xxxxy_xxxyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 62); 

                auto tg_xxxxy_xxxyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 63); 

                auto tg_xxxxy_xxxyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 64); 

                auto tg_xxxxy_xxxzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 65); 

                auto tg_xxxxy_xxyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 66); 

                auto tg_xxxxy_xxyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 67); 

                auto tg_xxxxy_xxyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 68); 

                auto tg_xxxxy_xxyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 69); 

                auto tg_xxxxy_xxyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 70); 

                auto tg_xxxxy_xxyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 71); 

                auto tg_xxxxy_xxzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 72); 

                auto tg_xxxxy_xyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 73); 

                auto tg_xxxxy_xyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 74); 

                auto tg_xxxxy_xyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 75); 

                auto tg_xxxxy_xyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 76); 

                auto tg_xxxxy_xyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 77); 

                auto tg_xxxxy_xyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 78); 

                auto tg_xxxxy_xyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 79); 

                auto tg_xxxxy_xzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 80); 

                auto tg_xxxxy_yyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 81); 

                auto tg_xxxxy_yyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 82); 

                auto tg_xxxxy_yyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 83); 

                auto tg_xxxxy_yyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 84); 

                auto tg_xxxxy_yyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 85); 

                auto tg_xxxxy_yyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 86); 

                auto tg_xxxxy_yyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 87); 

                auto tg_xxxxy_yzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 88); 

                auto tg_xxxxy_zzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 89); 

                auto tg_xxxxz_xxxxxxxx_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 90); 

                auto tg_xxxxz_xxxxxxxy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 91); 

                auto tg_xxxxz_xxxxxxxz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 92); 

                auto tg_xxxxz_xxxxxxyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 93); 

                auto tg_xxxxz_xxxxxxyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 94); 

                // Batch of Integrals (0,95)

                #pragma omp simd aligned(fxn, fza, tg_xxx_xxxxxxxx_0, tg_xxx_xxxxxxxx_1, tg_xxx_xxxxxxxy_0, \
                                         tg_xxx_xxxxxxxy_1, tg_xxx_xxxxxxxz_0, tg_xxx_xxxxxxxz_1, tg_xxx_xxxxxxyy_0, \
                                         tg_xxx_xxxxxxyy_1, tg_xxx_xxxxxxyz_0, tg_xxx_xxxxxxyz_1, tg_xxx_xxxxxxzz_0, \
                                         tg_xxx_xxxxxxzz_1, tg_xxx_xxxxxyyy_0, tg_xxx_xxxxxyyy_1, tg_xxx_xxxxxyyz_0, \
                                         tg_xxx_xxxxxyyz_1, tg_xxx_xxxxxyzz_0, tg_xxx_xxxxxyzz_1, tg_xxx_xxxxxzzz_0, \
                                         tg_xxx_xxxxxzzz_1, tg_xxx_xxxxyyyy_0, tg_xxx_xxxxyyyy_1, tg_xxx_xxxxyyyz_0, \
                                         tg_xxx_xxxxyyyz_1, tg_xxx_xxxxyyzz_0, tg_xxx_xxxxyyzz_1, tg_xxx_xxxxyzzz_0, \
                                         tg_xxx_xxxxyzzz_1, tg_xxx_xxxxzzzz_0, tg_xxx_xxxxzzzz_1, tg_xxx_xxxyyyyy_0, \
                                         tg_xxx_xxxyyyyy_1, tg_xxx_xxxyyyyz_0, tg_xxx_xxxyyyyz_1, tg_xxx_xxxyyyzz_0, \
                                         tg_xxx_xxxyyyzz_1, tg_xxx_xxxyyzzz_0, tg_xxx_xxxyyzzz_1, tg_xxx_xxxyzzzz_0, \
                                         tg_xxx_xxxyzzzz_1, tg_xxx_xxxzzzzz_0, tg_xxx_xxxzzzzz_1, tg_xxx_xxyyyyyy_0, \
                                         tg_xxx_xxyyyyyy_1, tg_xxx_xxyyyyyz_0, tg_xxx_xxyyyyyz_1, tg_xxx_xxyyyyzz_0, \
                                         tg_xxx_xxyyyyzz_1, tg_xxx_xxyyyzzz_0, tg_xxx_xxyyyzzz_1, tg_xxx_xxyyzzzz_0, \
                                         tg_xxx_xxyyzzzz_1, tg_xxx_xxyzzzzz_0, tg_xxx_xxyzzzzz_1, tg_xxx_xxzzzzzz_0, \
                                         tg_xxx_xxzzzzzz_1, tg_xxx_xyyyyyyy_0, tg_xxx_xyyyyyyy_1, tg_xxx_xyyyyyyz_0, \
                                         tg_xxx_xyyyyyyz_1, tg_xxx_xyyyyyzz_0, tg_xxx_xyyyyyzz_1, tg_xxx_xyyyyzzz_0, \
                                         tg_xxx_xyyyyzzz_1, tg_xxx_xyyyzzzz_0, tg_xxx_xyyyzzzz_1, tg_xxx_xyyzzzzz_0, \
                                         tg_xxx_xyyzzzzz_1, tg_xxx_xyzzzzzz_0, tg_xxx_xyzzzzzz_1, tg_xxx_xzzzzzzz_0, \
                                         tg_xxx_xzzzzzzz_1, tg_xxx_yyyyyyyy_0, tg_xxx_yyyyyyyy_1, tg_xxx_yyyyyyyz_0, \
                                         tg_xxx_yyyyyyyz_1, tg_xxx_yyyyyyzz_0, tg_xxx_yyyyyyzz_1, tg_xxx_yyyyyzzz_0, \
                                         tg_xxx_yyyyyzzz_1, tg_xxx_yyyyzzzz_0, tg_xxx_yyyyzzzz_1, tg_xxx_yyyzzzzz_0, \
                                         tg_xxx_yyyzzzzz_1, tg_xxx_yyzzzzzz_0, tg_xxx_yyzzzzzz_1, tg_xxx_yzzzzzzz_0, \
                                         tg_xxx_yzzzzzzz_1, tg_xxx_zzzzzzzz_0, tg_xxx_zzzzzzzz_1, tg_xxxx_xxxxxxx_1, \
                                         tg_xxxx_xxxxxxxx_0, tg_xxxx_xxxxxxxx_1, tg_xxxx_xxxxxxxy_0, tg_xxxx_xxxxxxxy_1, \
                                         tg_xxxx_xxxxxxxz_0, tg_xxxx_xxxxxxxz_1, tg_xxxx_xxxxxxy_1, tg_xxxx_xxxxxxyy_0, \
                                         tg_xxxx_xxxxxxyy_1, tg_xxxx_xxxxxxyz_0, tg_xxxx_xxxxxxyz_1, tg_xxxx_xxxxxxz_1, \
                                         tg_xxxx_xxxxxxzz_0, tg_xxxx_xxxxxxzz_1, tg_xxxx_xxxxxyy_1, tg_xxxx_xxxxxyyy_0, \
                                         tg_xxxx_xxxxxyyy_1, tg_xxxx_xxxxxyyz_0, tg_xxxx_xxxxxyyz_1, tg_xxxx_xxxxxyz_1, \
                                         tg_xxxx_xxxxxyzz_0, tg_xxxx_xxxxxyzz_1, tg_xxxx_xxxxxzz_1, tg_xxxx_xxxxxzzz_0, \
                                         tg_xxxx_xxxxxzzz_1, tg_xxxx_xxxxyyy_1, tg_xxxx_xxxxyyyy_0, tg_xxxx_xxxxyyyy_1, \
                                         tg_xxxx_xxxxyyyz_0, tg_xxxx_xxxxyyyz_1, tg_xxxx_xxxxyyz_1, tg_xxxx_xxxxyyzz_0, \
                                         tg_xxxx_xxxxyyzz_1, tg_xxxx_xxxxyzz_1, tg_xxxx_xxxxyzzz_0, tg_xxxx_xxxxyzzz_1, \
                                         tg_xxxx_xxxxzzz_1, tg_xxxx_xxxxzzzz_0, tg_xxxx_xxxxzzzz_1, tg_xxxx_xxxyyyy_1, \
                                         tg_xxxx_xxxyyyyy_0, tg_xxxx_xxxyyyyy_1, tg_xxxx_xxxyyyyz_0, tg_xxxx_xxxyyyyz_1, \
                                         tg_xxxx_xxxyyyz_1, tg_xxxx_xxxyyyzz_0, tg_xxxx_xxxyyyzz_1, tg_xxxx_xxxyyzz_1, \
                                         tg_xxxx_xxxyyzzz_0, tg_xxxx_xxxyyzzz_1, tg_xxxx_xxxyzzz_1, tg_xxxx_xxxyzzzz_0, \
                                         tg_xxxx_xxxyzzzz_1, tg_xxxx_xxxzzzz_1, tg_xxxx_xxxzzzzz_0, tg_xxxx_xxxzzzzz_1, \
                                         tg_xxxx_xxyyyyy_1, tg_xxxx_xxyyyyyy_0, tg_xxxx_xxyyyyyy_1, tg_xxxx_xxyyyyyz_0, \
                                         tg_xxxx_xxyyyyyz_1, tg_xxxx_xxyyyyz_1, tg_xxxx_xxyyyyzz_0, tg_xxxx_xxyyyyzz_1, \
                                         tg_xxxx_xxyyyzz_1, tg_xxxx_xxyyyzzz_0, tg_xxxx_xxyyyzzz_1, tg_xxxx_xxyyzzz_1, \
                                         tg_xxxx_xxyyzzzz_0, tg_xxxx_xxyyzzzz_1, tg_xxxx_xxyzzzz_1, tg_xxxx_xxyzzzzz_0, \
                                         tg_xxxx_xxyzzzzz_1, tg_xxxx_xxzzzzz_1, tg_xxxx_xxzzzzzz_0, tg_xxxx_xxzzzzzz_1, \
                                         tg_xxxx_xyyyyyy_1, tg_xxxx_xyyyyyyy_0, tg_xxxx_xyyyyyyy_1, tg_xxxx_xyyyyyyz_0, \
                                         tg_xxxx_xyyyyyyz_1, tg_xxxx_xyyyyyz_1, tg_xxxx_xyyyyyzz_0, tg_xxxx_xyyyyyzz_1, \
                                         tg_xxxx_xyyyyzz_1, tg_xxxx_xyyyyzzz_0, tg_xxxx_xyyyyzzz_1, tg_xxxx_xyyyzzz_1, \
                                         tg_xxxx_xyyyzzzz_0, tg_xxxx_xyyyzzzz_1, tg_xxxx_xyyzzzz_1, tg_xxxx_xyyzzzzz_0, \
                                         tg_xxxx_xyyzzzzz_1, tg_xxxx_xyzzzzz_1, tg_xxxx_xyzzzzzz_0, tg_xxxx_xyzzzzzz_1, \
                                         tg_xxxx_xzzzzzz_1, tg_xxxx_xzzzzzzz_0, tg_xxxx_xzzzzzzz_1, tg_xxxx_yyyyyyy_1, \
                                         tg_xxxx_yyyyyyyy_0, tg_xxxx_yyyyyyyy_1, tg_xxxx_yyyyyyyz_0, tg_xxxx_yyyyyyyz_1, \
                                         tg_xxxx_yyyyyyz_1, tg_xxxx_yyyyyyzz_0, tg_xxxx_yyyyyyzz_1, tg_xxxx_yyyyyzz_1, \
                                         tg_xxxx_yyyyyzzz_0, tg_xxxx_yyyyyzzz_1, tg_xxxx_yyyyzzz_1, tg_xxxx_yyyyzzzz_0, \
                                         tg_xxxx_yyyyzzzz_1, tg_xxxx_yyyzzzz_1, tg_xxxx_yyyzzzzz_0, tg_xxxx_yyyzzzzz_1, \
                                         tg_xxxx_yyzzzzz_1, tg_xxxx_yyzzzzzz_0, tg_xxxx_yyzzzzzz_1, tg_xxxx_yzzzzzz_1, \
                                         tg_xxxx_yzzzzzzz_0, tg_xxxx_yzzzzzzz_1, tg_xxxx_zzzzzzz_1, tg_xxxx_zzzzzzzz_0, \
                                         tg_xxxx_zzzzzzzz_1, tg_xxxxx_xxxxxxxx_0, tg_xxxxx_xxxxxxxy_0, tg_xxxxx_xxxxxxxz_0, \
                                         tg_xxxxx_xxxxxxyy_0, tg_xxxxx_xxxxxxyz_0, tg_xxxxx_xxxxxxzz_0, tg_xxxxx_xxxxxyyy_0, \
                                         tg_xxxxx_xxxxxyyz_0, tg_xxxxx_xxxxxyzz_0, tg_xxxxx_xxxxxzzz_0, tg_xxxxx_xxxxyyyy_0, \
                                         tg_xxxxx_xxxxyyyz_0, tg_xxxxx_xxxxyyzz_0, tg_xxxxx_xxxxyzzz_0, tg_xxxxx_xxxxzzzz_0, \
                                         tg_xxxxx_xxxyyyyy_0, tg_xxxxx_xxxyyyyz_0, tg_xxxxx_xxxyyyzz_0, tg_xxxxx_xxxyyzzz_0, \
                                         tg_xxxxx_xxxyzzzz_0, tg_xxxxx_xxxzzzzz_0, tg_xxxxx_xxyyyyyy_0, tg_xxxxx_xxyyyyyz_0, \
                                         tg_xxxxx_xxyyyyzz_0, tg_xxxxx_xxyyyzzz_0, tg_xxxxx_xxyyzzzz_0, tg_xxxxx_xxyzzzzz_0, \
                                         tg_xxxxx_xxzzzzzz_0, tg_xxxxx_xyyyyyyy_0, tg_xxxxx_xyyyyyyz_0, tg_xxxxx_xyyyyyzz_0, \
                                         tg_xxxxx_xyyyyzzz_0, tg_xxxxx_xyyyzzzz_0, tg_xxxxx_xyyzzzzz_0, tg_xxxxx_xyzzzzzz_0, \
                                         tg_xxxxx_xzzzzzzz_0, tg_xxxxx_yyyyyyyy_0, tg_xxxxx_yyyyyyyz_0, tg_xxxxx_yyyyyyzz_0, \
                                         tg_xxxxx_yyyyyzzz_0, tg_xxxxx_yyyyzzzz_0, tg_xxxxx_yyyzzzzz_0, tg_xxxxx_yyzzzzzz_0, \
                                         tg_xxxxx_yzzzzzzz_0, tg_xxxxx_zzzzzzzz_0, tg_xxxxy_xxxxxxxx_0, tg_xxxxy_xxxxxxxy_0, \
                                         tg_xxxxy_xxxxxxxz_0, tg_xxxxy_xxxxxxyy_0, tg_xxxxy_xxxxxxyz_0, tg_xxxxy_xxxxxxzz_0, \
                                         tg_xxxxy_xxxxxyyy_0, tg_xxxxy_xxxxxyyz_0, tg_xxxxy_xxxxxyzz_0, tg_xxxxy_xxxxxzzz_0, \
                                         tg_xxxxy_xxxxyyyy_0, tg_xxxxy_xxxxyyyz_0, tg_xxxxy_xxxxyyzz_0, tg_xxxxy_xxxxyzzz_0, \
                                         tg_xxxxy_xxxxzzzz_0, tg_xxxxy_xxxyyyyy_0, tg_xxxxy_xxxyyyyz_0, tg_xxxxy_xxxyyyzz_0, \
                                         tg_xxxxy_xxxyyzzz_0, tg_xxxxy_xxxyzzzz_0, tg_xxxxy_xxxzzzzz_0, tg_xxxxy_xxyyyyyy_0, \
                                         tg_xxxxy_xxyyyyyz_0, tg_xxxxy_xxyyyyzz_0, tg_xxxxy_xxyyyzzz_0, tg_xxxxy_xxyyzzzz_0, \
                                         tg_xxxxy_xxyzzzzz_0, tg_xxxxy_xxzzzzzz_0, tg_xxxxy_xyyyyyyy_0, tg_xxxxy_xyyyyyyz_0, \
                                         tg_xxxxy_xyyyyyzz_0, tg_xxxxy_xyyyyzzz_0, tg_xxxxy_xyyyzzzz_0, tg_xxxxy_xyyzzzzz_0, \
                                         tg_xxxxy_xyzzzzzz_0, tg_xxxxy_xzzzzzzz_0, tg_xxxxy_yyyyyyyy_0, tg_xxxxy_yyyyyyyz_0, \
                                         tg_xxxxy_yyyyyyzz_0, tg_xxxxy_yyyyyzzz_0, tg_xxxxy_yyyyzzzz_0, tg_xxxxy_yyyzzzzz_0, \
                                         tg_xxxxy_yyzzzzzz_0, tg_xxxxy_yzzzzzzz_0, tg_xxxxy_zzzzzzzz_0, tg_xxxxz_xxxxxxxx_0, \
                                         tg_xxxxz_xxxxxxxy_0, tg_xxxxz_xxxxxxxz_0, tg_xxxxz_xxxxxxyy_0, tg_xxxxz_xxxxxxyz_0, \
                                         tg_xxxy_xxxxxxx_1, tg_xxxy_xxxxxxxx_0, tg_xxxy_xxxxxxxx_1, tg_xxxy_xxxxxxxy_0, \
                                         tg_xxxy_xxxxxxxy_1, tg_xxxy_xxxxxxxz_0, tg_xxxy_xxxxxxxz_1, tg_xxxy_xxxxxxy_1, \
                                         tg_xxxy_xxxxxxyy_0, tg_xxxy_xxxxxxyy_1, tg_xxxy_xxxxxxyz_0, tg_xxxy_xxxxxxyz_1, \
                                         tg_xxxy_xxxxxxz_1, tg_xxxy_xxxxxxzz_0, tg_xxxy_xxxxxxzz_1, tg_xxxy_xxxxxyy_1, \
                                         tg_xxxy_xxxxxyyy_0, tg_xxxy_xxxxxyyy_1, tg_xxxy_xxxxxyyz_0, tg_xxxy_xxxxxyyz_1, \
                                         tg_xxxy_xxxxxyz_1, tg_xxxy_xxxxxyzz_0, tg_xxxy_xxxxxyzz_1, tg_xxxy_xxxxxzz_1, \
                                         tg_xxxy_xxxxxzzz_0, tg_xxxy_xxxxxzzz_1, tg_xxxy_xxxxyyy_1, tg_xxxy_xxxxyyyy_0, \
                                         tg_xxxy_xxxxyyyy_1, tg_xxxy_xxxxyyyz_0, tg_xxxy_xxxxyyyz_1, tg_xxxy_xxxxyyz_1, \
                                         tg_xxxy_xxxxyyzz_0, tg_xxxy_xxxxyyzz_1, tg_xxxy_xxxxyzz_1, tg_xxxy_xxxxyzzz_0, \
                                         tg_xxxy_xxxxyzzz_1, tg_xxxy_xxxxzzz_1, tg_xxxy_xxxxzzzz_0, tg_xxxy_xxxxzzzz_1, \
                                         tg_xxxy_xxxyyyy_1, tg_xxxy_xxxyyyyy_0, tg_xxxy_xxxyyyyy_1, tg_xxxy_xxxyyyyz_0, \
                                         tg_xxxy_xxxyyyyz_1, tg_xxxy_xxxyyyz_1, tg_xxxy_xxxyyyzz_0, tg_xxxy_xxxyyyzz_1, \
                                         tg_xxxy_xxxyyzz_1, tg_xxxy_xxxyyzzz_0, tg_xxxy_xxxyyzzz_1, tg_xxxy_xxxyzzz_1, \
                                         tg_xxxy_xxxyzzzz_0, tg_xxxy_xxxyzzzz_1, tg_xxxy_xxxzzzz_1, tg_xxxy_xxxzzzzz_0, \
                                         tg_xxxy_xxxzzzzz_1, tg_xxxy_xxyyyyy_1, tg_xxxy_xxyyyyyy_0, tg_xxxy_xxyyyyyy_1, \
                                         tg_xxxy_xxyyyyyz_0, tg_xxxy_xxyyyyyz_1, tg_xxxy_xxyyyyz_1, tg_xxxy_xxyyyyzz_0, \
                                         tg_xxxy_xxyyyyzz_1, tg_xxxy_xxyyyzz_1, tg_xxxy_xxyyyzzz_0, tg_xxxy_xxyyyzzz_1, \
                                         tg_xxxy_xxyyzzz_1, tg_xxxy_xxyyzzzz_0, tg_xxxy_xxyyzzzz_1, tg_xxxy_xxyzzzz_1, \
                                         tg_xxxy_xxyzzzzz_0, tg_xxxy_xxyzzzzz_1, tg_xxxy_xxzzzzz_1, tg_xxxy_xxzzzzzz_0, \
                                         tg_xxxy_xxzzzzzz_1, tg_xxxy_xyyyyyy_1, tg_xxxy_xyyyyyyy_0, tg_xxxy_xyyyyyyy_1, \
                                         tg_xxxy_xyyyyyyz_0, tg_xxxy_xyyyyyyz_1, tg_xxxy_xyyyyyz_1, tg_xxxy_xyyyyyzz_0, \
                                         tg_xxxy_xyyyyyzz_1, tg_xxxy_xyyyyzz_1, tg_xxxy_xyyyyzzz_0, tg_xxxy_xyyyyzzz_1, \
                                         tg_xxxy_xyyyzzz_1, tg_xxxy_xyyyzzzz_0, tg_xxxy_xyyyzzzz_1, tg_xxxy_xyyzzzz_1, \
                                         tg_xxxy_xyyzzzzz_0, tg_xxxy_xyyzzzzz_1, tg_xxxy_xyzzzzz_1, tg_xxxy_xyzzzzzz_0, \
                                         tg_xxxy_xyzzzzzz_1, tg_xxxy_xzzzzzz_1, tg_xxxy_xzzzzzzz_0, tg_xxxy_xzzzzzzz_1, \
                                         tg_xxxy_yyyyyyy_1, tg_xxxy_yyyyyyyy_0, tg_xxxy_yyyyyyyy_1, tg_xxxy_yyyyyyyz_0, \
                                         tg_xxxy_yyyyyyyz_1, tg_xxxy_yyyyyyz_1, tg_xxxy_yyyyyyzz_0, tg_xxxy_yyyyyyzz_1, \
                                         tg_xxxy_yyyyyzz_1, tg_xxxy_yyyyyzzz_0, tg_xxxy_yyyyyzzz_1, tg_xxxy_yyyyzzz_1, \
                                         tg_xxxy_yyyyzzzz_0, tg_xxxy_yyyyzzzz_1, tg_xxxy_yyyzzzz_1, tg_xxxy_yyyzzzzz_0, \
                                         tg_xxxy_yyyzzzzz_1, tg_xxxy_yyzzzzz_1, tg_xxxy_yyzzzzzz_0, tg_xxxy_yyzzzzzz_1, \
                                         tg_xxxy_yzzzzzz_1, tg_xxxy_yzzzzzzz_0, tg_xxxy_yzzzzzzz_1, tg_xxxy_zzzzzzz_1, \
                                         tg_xxxy_zzzzzzzz_0, tg_xxxy_zzzzzzzz_1, tg_xxxz_xxxxxxx_1, tg_xxxz_xxxxxxxx_0, \
                                         tg_xxxz_xxxxxxxx_1, tg_xxxz_xxxxxxxy_0, tg_xxxz_xxxxxxxy_1, tg_xxxz_xxxxxxxz_0, \
                                         tg_xxxz_xxxxxxxz_1, tg_xxxz_xxxxxxy_1, tg_xxxz_xxxxxxyy_0, tg_xxxz_xxxxxxyy_1, \
                                         tg_xxxz_xxxxxxyz_0, tg_xxxz_xxxxxxyz_1, tg_xxxz_xxxxxxz_1, tg_xxxz_xxxxxyy_1, \
                                         tg_xxxz_xxxxxyz_1, tg_xxy_xxxxxxxx_0, tg_xxy_xxxxxxxx_1, tg_xxy_xxxxxxxy_0, \
                                         tg_xxy_xxxxxxxy_1, tg_xxy_xxxxxxxz_0, tg_xxy_xxxxxxxz_1, tg_xxy_xxxxxxyy_0, \
                                         tg_xxy_xxxxxxyy_1, tg_xxy_xxxxxxyz_0, tg_xxy_xxxxxxyz_1, tg_xxy_xxxxxxzz_0, \
                                         tg_xxy_xxxxxxzz_1, tg_xxy_xxxxxyyy_0, tg_xxy_xxxxxyyy_1, tg_xxy_xxxxxyyz_0, \
                                         tg_xxy_xxxxxyyz_1, tg_xxy_xxxxxyzz_0, tg_xxy_xxxxxyzz_1, tg_xxy_xxxxxzzz_0, \
                                         tg_xxy_xxxxxzzz_1, tg_xxy_xxxxyyyy_0, tg_xxy_xxxxyyyy_1, tg_xxy_xxxxyyyz_0, \
                                         tg_xxy_xxxxyyyz_1, tg_xxy_xxxxyyzz_0, tg_xxy_xxxxyyzz_1, tg_xxy_xxxxyzzz_0, \
                                         tg_xxy_xxxxyzzz_1, tg_xxy_xxxxzzzz_0, tg_xxy_xxxxzzzz_1, tg_xxy_xxxyyyyy_0, \
                                         tg_xxy_xxxyyyyy_1, tg_xxy_xxxyyyyz_0, tg_xxy_xxxyyyyz_1, tg_xxy_xxxyyyzz_0, \
                                         tg_xxy_xxxyyyzz_1, tg_xxy_xxxyyzzz_0, tg_xxy_xxxyyzzz_1, tg_xxy_xxxyzzzz_0, \
                                         tg_xxy_xxxyzzzz_1, tg_xxy_xxxzzzzz_0, tg_xxy_xxxzzzzz_1, tg_xxy_xxyyyyyy_0, \
                                         tg_xxy_xxyyyyyy_1, tg_xxy_xxyyyyyz_0, tg_xxy_xxyyyyyz_1, tg_xxy_xxyyyyzz_0, \
                                         tg_xxy_xxyyyyzz_1, tg_xxy_xxyyyzzz_0, tg_xxy_xxyyyzzz_1, tg_xxy_xxyyzzzz_0, \
                                         tg_xxy_xxyyzzzz_1, tg_xxy_xxyzzzzz_0, tg_xxy_xxyzzzzz_1, tg_xxy_xxzzzzzz_0, \
                                         tg_xxy_xxzzzzzz_1, tg_xxy_xyyyyyyy_0, tg_xxy_xyyyyyyy_1, tg_xxy_xyyyyyyz_0, \
                                         tg_xxy_xyyyyyyz_1, tg_xxy_xyyyyyzz_0, tg_xxy_xyyyyyzz_1, tg_xxy_xyyyyzzz_0, \
                                         tg_xxy_xyyyyzzz_1, tg_xxy_xyyyzzzz_0, tg_xxy_xyyyzzzz_1, tg_xxy_xyyzzzzz_0, \
                                         tg_xxy_xyyzzzzz_1, tg_xxy_xyzzzzzz_0, tg_xxy_xyzzzzzz_1, tg_xxy_xzzzzzzz_0, \
                                         tg_xxy_xzzzzzzz_1, tg_xxy_yyyyyyyy_0, tg_xxy_yyyyyyyy_1, tg_xxy_yyyyyyyz_0, \
                                         tg_xxy_yyyyyyyz_1, tg_xxy_yyyyyyzz_0, tg_xxy_yyyyyyzz_1, tg_xxy_yyyyyzzz_0, \
                                         tg_xxy_yyyyyzzz_1, tg_xxy_yyyyzzzz_0, tg_xxy_yyyyzzzz_1, tg_xxy_yyyzzzzz_0, \
                                         tg_xxy_yyyzzzzz_1, tg_xxy_yyzzzzzz_0, tg_xxy_yyzzzzzz_1, tg_xxy_yzzzzzzz_0, \
                                         tg_xxy_yzzzzzzz_1, tg_xxy_zzzzzzzz_0, tg_xxy_zzzzzzzz_1, tg_xxz_xxxxxxxx_0, \
                                         tg_xxz_xxxxxxxx_1, tg_xxz_xxxxxxxy_0, tg_xxz_xxxxxxxy_1, tg_xxz_xxxxxxxz_0, \
                                         tg_xxz_xxxxxxxz_1, tg_xxz_xxxxxxyy_0, tg_xxz_xxxxxxyy_1, tg_xxz_xxxxxxyz_0, \
                                         tg_xxz_xxxxxxyz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxxx_xxxxxxxx_0[j] = pb_x * tg_xxxx_xxxxxxxx_0[j] + fr * tg_xxxx_xxxxxxxx_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxxxxxx_0[j] - tg_xxx_xxxxxxxx_1[j] * fl1_fza) + 4.0 * fl1_fxn * tg_xxxx_xxxxxxx_1[j];

                    tg_xxxxx_xxxxxxxy_0[j] = pb_x * tg_xxxx_xxxxxxxy_0[j] + fr * tg_xxxx_xxxxxxxy_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxxxxxy_0[j] - tg_xxx_xxxxxxxy_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_xxxx_xxxxxxy_1[j];

                    tg_xxxxx_xxxxxxxz_0[j] = pb_x * tg_xxxx_xxxxxxxz_0[j] + fr * tg_xxxx_xxxxxxxz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxxxxxz_0[j] - tg_xxx_xxxxxxxz_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_xxxx_xxxxxxz_1[j];

                    tg_xxxxx_xxxxxxyy_0[j] = pb_x * tg_xxxx_xxxxxxyy_0[j] + fr * tg_xxxx_xxxxxxyy_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxxxxyy_0[j] - tg_xxx_xxxxxxyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxxx_xxxxxyy_1[j];

                    tg_xxxxx_xxxxxxyz_0[j] = pb_x * tg_xxxx_xxxxxxyz_0[j] + fr * tg_xxxx_xxxxxxyz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxxxxyz_0[j] - tg_xxx_xxxxxxyz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxxx_xxxxxyz_1[j];

                    tg_xxxxx_xxxxxxzz_0[j] = pb_x * tg_xxxx_xxxxxxzz_0[j] + fr * tg_xxxx_xxxxxxzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxxxxzz_0[j] - tg_xxx_xxxxxxzz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxxx_xxxxxzz_1[j];

                    tg_xxxxx_xxxxxyyy_0[j] = pb_x * tg_xxxx_xxxxxyyy_0[j] + fr * tg_xxxx_xxxxxyyy_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxxxyyy_0[j] - tg_xxx_xxxxxyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxx_xxxxyyy_1[j];

                    tg_xxxxx_xxxxxyyz_0[j] = pb_x * tg_xxxx_xxxxxyyz_0[j] + fr * tg_xxxx_xxxxxyyz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxxxyyz_0[j] - tg_xxx_xxxxxyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxx_xxxxyyz_1[j];

                    tg_xxxxx_xxxxxyzz_0[j] = pb_x * tg_xxxx_xxxxxyzz_0[j] + fr * tg_xxxx_xxxxxyzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxxxyzz_0[j] - tg_xxx_xxxxxyzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxx_xxxxyzz_1[j];

                    tg_xxxxx_xxxxxzzz_0[j] = pb_x * tg_xxxx_xxxxxzzz_0[j] + fr * tg_xxxx_xxxxxzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxxxzzz_0[j] - tg_xxx_xxxxxzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxx_xxxxzzz_1[j];

                    tg_xxxxx_xxxxyyyy_0[j] = pb_x * tg_xxxx_xxxxyyyy_0[j] + fr * tg_xxxx_xxxxyyyy_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxxyyyy_0[j] - tg_xxx_xxxxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxx_xxxyyyy_1[j];

                    tg_xxxxx_xxxxyyyz_0[j] = pb_x * tg_xxxx_xxxxyyyz_0[j] + fr * tg_xxxx_xxxxyyyz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxxyyyz_0[j] - tg_xxx_xxxxyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxx_xxxyyyz_1[j];

                    tg_xxxxx_xxxxyyzz_0[j] = pb_x * tg_xxxx_xxxxyyzz_0[j] + fr * tg_xxxx_xxxxyyzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxxyyzz_0[j] - tg_xxx_xxxxyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxx_xxxyyzz_1[j];

                    tg_xxxxx_xxxxyzzz_0[j] = pb_x * tg_xxxx_xxxxyzzz_0[j] + fr * tg_xxxx_xxxxyzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxxyzzz_0[j] - tg_xxx_xxxxyzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxx_xxxyzzz_1[j];

                    tg_xxxxx_xxxxzzzz_0[j] = pb_x * tg_xxxx_xxxxzzzz_0[j] + fr * tg_xxxx_xxxxzzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxxzzzz_0[j] - tg_xxx_xxxxzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxx_xxxzzzz_1[j];

                    tg_xxxxx_xxxyyyyy_0[j] = pb_x * tg_xxxx_xxxyyyyy_0[j] + fr * tg_xxxx_xxxyyyyy_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxyyyyy_0[j] - tg_xxx_xxxyyyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxx_xxyyyyy_1[j];

                    tg_xxxxx_xxxyyyyz_0[j] = pb_x * tg_xxxx_xxxyyyyz_0[j] + fr * tg_xxxx_xxxyyyyz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxyyyyz_0[j] - tg_xxx_xxxyyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxx_xxyyyyz_1[j];

                    tg_xxxxx_xxxyyyzz_0[j] = pb_x * tg_xxxx_xxxyyyzz_0[j] + fr * tg_xxxx_xxxyyyzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxyyyzz_0[j] - tg_xxx_xxxyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxx_xxyyyzz_1[j];

                    tg_xxxxx_xxxyyzzz_0[j] = pb_x * tg_xxxx_xxxyyzzz_0[j] + fr * tg_xxxx_xxxyyzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxyyzzz_0[j] - tg_xxx_xxxyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxx_xxyyzzz_1[j];

                    tg_xxxxx_xxxyzzzz_0[j] = pb_x * tg_xxxx_xxxyzzzz_0[j] + fr * tg_xxxx_xxxyzzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxyzzzz_0[j] - tg_xxx_xxxyzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxx_xxyzzzz_1[j];

                    tg_xxxxx_xxxzzzzz_0[j] = pb_x * tg_xxxx_xxxzzzzz_0[j] + fr * tg_xxxx_xxxzzzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxzzzzz_0[j] - tg_xxx_xxxzzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxx_xxzzzzz_1[j];

                    tg_xxxxx_xxyyyyyy_0[j] = pb_x * tg_xxxx_xxyyyyyy_0[j] + fr * tg_xxxx_xxyyyyyy_1[j] + 2.0 * fl1_fx * (tg_xxx_xxyyyyyy_0[j] - tg_xxx_xxyyyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxx_xyyyyyy_1[j];

                    tg_xxxxx_xxyyyyyz_0[j] = pb_x * tg_xxxx_xxyyyyyz_0[j] + fr * tg_xxxx_xxyyyyyz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxyyyyyz_0[j] - tg_xxx_xxyyyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxx_xyyyyyz_1[j];

                    tg_xxxxx_xxyyyyzz_0[j] = pb_x * tg_xxxx_xxyyyyzz_0[j] + fr * tg_xxxx_xxyyyyzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxyyyyzz_0[j] - tg_xxx_xxyyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxx_xyyyyzz_1[j];

                    tg_xxxxx_xxyyyzzz_0[j] = pb_x * tg_xxxx_xxyyyzzz_0[j] + fr * tg_xxxx_xxyyyzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxyyyzzz_0[j] - tg_xxx_xxyyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxx_xyyyzzz_1[j];

                    tg_xxxxx_xxyyzzzz_0[j] = pb_x * tg_xxxx_xxyyzzzz_0[j] + fr * tg_xxxx_xxyyzzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxyyzzzz_0[j] - tg_xxx_xxyyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxx_xyyzzzz_1[j];

                    tg_xxxxx_xxyzzzzz_0[j] = pb_x * tg_xxxx_xxyzzzzz_0[j] + fr * tg_xxxx_xxyzzzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxyzzzzz_0[j] - tg_xxx_xxyzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxx_xyzzzzz_1[j];

                    tg_xxxxx_xxzzzzzz_0[j] = pb_x * tg_xxxx_xxzzzzzz_0[j] + fr * tg_xxxx_xxzzzzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxzzzzzz_0[j] - tg_xxx_xxzzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxx_xzzzzzz_1[j];

                    tg_xxxxx_xyyyyyyy_0[j] = pb_x * tg_xxxx_xyyyyyyy_0[j] + fr * tg_xxxx_xyyyyyyy_1[j] + 2.0 * fl1_fx * (tg_xxx_xyyyyyyy_0[j] - tg_xxx_xyyyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxx_yyyyyyy_1[j];

                    tg_xxxxx_xyyyyyyz_0[j] = pb_x * tg_xxxx_xyyyyyyz_0[j] + fr * tg_xxxx_xyyyyyyz_1[j] + 2.0 * fl1_fx * (tg_xxx_xyyyyyyz_0[j] - tg_xxx_xyyyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxx_yyyyyyz_1[j];

                    tg_xxxxx_xyyyyyzz_0[j] = pb_x * tg_xxxx_xyyyyyzz_0[j] + fr * tg_xxxx_xyyyyyzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xyyyyyzz_0[j] - tg_xxx_xyyyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxx_yyyyyzz_1[j];

                    tg_xxxxx_xyyyyzzz_0[j] = pb_x * tg_xxxx_xyyyyzzz_0[j] + fr * tg_xxxx_xyyyyzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xyyyyzzz_0[j] - tg_xxx_xyyyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxx_yyyyzzz_1[j];

                    tg_xxxxx_xyyyzzzz_0[j] = pb_x * tg_xxxx_xyyyzzzz_0[j] + fr * tg_xxxx_xyyyzzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xyyyzzzz_0[j] - tg_xxx_xyyyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxx_yyyzzzz_1[j];

                    tg_xxxxx_xyyzzzzz_0[j] = pb_x * tg_xxxx_xyyzzzzz_0[j] + fr * tg_xxxx_xyyzzzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xyyzzzzz_0[j] - tg_xxx_xyyzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxx_yyzzzzz_1[j];

                    tg_xxxxx_xyzzzzzz_0[j] = pb_x * tg_xxxx_xyzzzzzz_0[j] + fr * tg_xxxx_xyzzzzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xyzzzzzz_0[j] - tg_xxx_xyzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxx_yzzzzzz_1[j];

                    tg_xxxxx_xzzzzzzz_0[j] = pb_x * tg_xxxx_xzzzzzzz_0[j] + fr * tg_xxxx_xzzzzzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xzzzzzzz_0[j] - tg_xxx_xzzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxx_zzzzzzz_1[j];

                    tg_xxxxx_yyyyyyyy_0[j] = pb_x * tg_xxxx_yyyyyyyy_0[j] + fr * tg_xxxx_yyyyyyyy_1[j] + 2.0 * fl1_fx * (tg_xxx_yyyyyyyy_0[j] - tg_xxx_yyyyyyyy_1[j] * fl1_fza);

                    tg_xxxxx_yyyyyyyz_0[j] = pb_x * tg_xxxx_yyyyyyyz_0[j] + fr * tg_xxxx_yyyyyyyz_1[j] + 2.0 * fl1_fx * (tg_xxx_yyyyyyyz_0[j] - tg_xxx_yyyyyyyz_1[j] * fl1_fza);

                    tg_xxxxx_yyyyyyzz_0[j] = pb_x * tg_xxxx_yyyyyyzz_0[j] + fr * tg_xxxx_yyyyyyzz_1[j] + 2.0 * fl1_fx * (tg_xxx_yyyyyyzz_0[j] - tg_xxx_yyyyyyzz_1[j] * fl1_fza);

                    tg_xxxxx_yyyyyzzz_0[j] = pb_x * tg_xxxx_yyyyyzzz_0[j] + fr * tg_xxxx_yyyyyzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_yyyyyzzz_0[j] - tg_xxx_yyyyyzzz_1[j] * fl1_fza);

                    tg_xxxxx_yyyyzzzz_0[j] = pb_x * tg_xxxx_yyyyzzzz_0[j] + fr * tg_xxxx_yyyyzzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_yyyyzzzz_0[j] - tg_xxx_yyyyzzzz_1[j] * fl1_fza);

                    tg_xxxxx_yyyzzzzz_0[j] = pb_x * tg_xxxx_yyyzzzzz_0[j] + fr * tg_xxxx_yyyzzzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_yyyzzzzz_0[j] - tg_xxx_yyyzzzzz_1[j] * fl1_fza);

                    tg_xxxxx_yyzzzzzz_0[j] = pb_x * tg_xxxx_yyzzzzzz_0[j] + fr * tg_xxxx_yyzzzzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_yyzzzzzz_0[j] - tg_xxx_yyzzzzzz_1[j] * fl1_fza);

                    tg_xxxxx_yzzzzzzz_0[j] = pb_x * tg_xxxx_yzzzzzzz_0[j] + fr * tg_xxxx_yzzzzzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_yzzzzzzz_0[j] - tg_xxx_yzzzzzzz_1[j] * fl1_fza);

                    tg_xxxxx_zzzzzzzz_0[j] = pb_x * tg_xxxx_zzzzzzzz_0[j] + fr * tg_xxxx_zzzzzzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_zzzzzzzz_0[j] - tg_xxx_zzzzzzzz_1[j] * fl1_fza);

                    tg_xxxxy_xxxxxxxx_0[j] = pb_x * tg_xxxy_xxxxxxxx_0[j] + fr * tg_xxxy_xxxxxxxx_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxxxxxx_0[j] - tg_xxy_xxxxxxxx_1[j] * fl1_fza) + 4.0 * fl1_fxn * tg_xxxy_xxxxxxx_1[j];

                    tg_xxxxy_xxxxxxxy_0[j] = pb_x * tg_xxxy_xxxxxxxy_0[j] + fr * tg_xxxy_xxxxxxxy_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxxxxxy_0[j] - tg_xxy_xxxxxxxy_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_xxxy_xxxxxxy_1[j];

                    tg_xxxxy_xxxxxxxz_0[j] = pb_x * tg_xxxy_xxxxxxxz_0[j] + fr * tg_xxxy_xxxxxxxz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxxxxxz_0[j] - tg_xxy_xxxxxxxz_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_xxxy_xxxxxxz_1[j];

                    tg_xxxxy_xxxxxxyy_0[j] = pb_x * tg_xxxy_xxxxxxyy_0[j] + fr * tg_xxxy_xxxxxxyy_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxxxxyy_0[j] - tg_xxy_xxxxxxyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxxy_xxxxxyy_1[j];

                    tg_xxxxy_xxxxxxyz_0[j] = pb_x * tg_xxxy_xxxxxxyz_0[j] + fr * tg_xxxy_xxxxxxyz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxxxxyz_0[j] - tg_xxy_xxxxxxyz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxxy_xxxxxyz_1[j];

                    tg_xxxxy_xxxxxxzz_0[j] = pb_x * tg_xxxy_xxxxxxzz_0[j] + fr * tg_xxxy_xxxxxxzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxxxxzz_0[j] - tg_xxy_xxxxxxzz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxxy_xxxxxzz_1[j];

                    tg_xxxxy_xxxxxyyy_0[j] = pb_x * tg_xxxy_xxxxxyyy_0[j] + fr * tg_xxxy_xxxxxyyy_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxxxyyy_0[j] - tg_xxy_xxxxxyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxy_xxxxyyy_1[j];

                    tg_xxxxy_xxxxxyyz_0[j] = pb_x * tg_xxxy_xxxxxyyz_0[j] + fr * tg_xxxy_xxxxxyyz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxxxyyz_0[j] - tg_xxy_xxxxxyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxy_xxxxyyz_1[j];

                    tg_xxxxy_xxxxxyzz_0[j] = pb_x * tg_xxxy_xxxxxyzz_0[j] + fr * tg_xxxy_xxxxxyzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxxxyzz_0[j] - tg_xxy_xxxxxyzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxy_xxxxyzz_1[j];

                    tg_xxxxy_xxxxxzzz_0[j] = pb_x * tg_xxxy_xxxxxzzz_0[j] + fr * tg_xxxy_xxxxxzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxxxzzz_0[j] - tg_xxy_xxxxxzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxy_xxxxzzz_1[j];

                    tg_xxxxy_xxxxyyyy_0[j] = pb_x * tg_xxxy_xxxxyyyy_0[j] + fr * tg_xxxy_xxxxyyyy_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxxyyyy_0[j] - tg_xxy_xxxxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxy_xxxyyyy_1[j];

                    tg_xxxxy_xxxxyyyz_0[j] = pb_x * tg_xxxy_xxxxyyyz_0[j] + fr * tg_xxxy_xxxxyyyz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxxyyyz_0[j] - tg_xxy_xxxxyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxy_xxxyyyz_1[j];

                    tg_xxxxy_xxxxyyzz_0[j] = pb_x * tg_xxxy_xxxxyyzz_0[j] + fr * tg_xxxy_xxxxyyzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxxyyzz_0[j] - tg_xxy_xxxxyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxy_xxxyyzz_1[j];

                    tg_xxxxy_xxxxyzzz_0[j] = pb_x * tg_xxxy_xxxxyzzz_0[j] + fr * tg_xxxy_xxxxyzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxxyzzz_0[j] - tg_xxy_xxxxyzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxy_xxxyzzz_1[j];

                    tg_xxxxy_xxxxzzzz_0[j] = pb_x * tg_xxxy_xxxxzzzz_0[j] + fr * tg_xxxy_xxxxzzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxxzzzz_0[j] - tg_xxy_xxxxzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxy_xxxzzzz_1[j];

                    tg_xxxxy_xxxyyyyy_0[j] = pb_x * tg_xxxy_xxxyyyyy_0[j] + fr * tg_xxxy_xxxyyyyy_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxyyyyy_0[j] - tg_xxy_xxxyyyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxy_xxyyyyy_1[j];

                    tg_xxxxy_xxxyyyyz_0[j] = pb_x * tg_xxxy_xxxyyyyz_0[j] + fr * tg_xxxy_xxxyyyyz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxyyyyz_0[j] - tg_xxy_xxxyyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxy_xxyyyyz_1[j];

                    tg_xxxxy_xxxyyyzz_0[j] = pb_x * tg_xxxy_xxxyyyzz_0[j] + fr * tg_xxxy_xxxyyyzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxyyyzz_0[j] - tg_xxy_xxxyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxy_xxyyyzz_1[j];

                    tg_xxxxy_xxxyyzzz_0[j] = pb_x * tg_xxxy_xxxyyzzz_0[j] + fr * tg_xxxy_xxxyyzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxyyzzz_0[j] - tg_xxy_xxxyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxy_xxyyzzz_1[j];

                    tg_xxxxy_xxxyzzzz_0[j] = pb_x * tg_xxxy_xxxyzzzz_0[j] + fr * tg_xxxy_xxxyzzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxyzzzz_0[j] - tg_xxy_xxxyzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxy_xxyzzzz_1[j];

                    tg_xxxxy_xxxzzzzz_0[j] = pb_x * tg_xxxy_xxxzzzzz_0[j] + fr * tg_xxxy_xxxzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxzzzzz_0[j] - tg_xxy_xxxzzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxy_xxzzzzz_1[j];

                    tg_xxxxy_xxyyyyyy_0[j] = pb_x * tg_xxxy_xxyyyyyy_0[j] + fr * tg_xxxy_xxyyyyyy_1[j] + 1.5 * fl1_fx * (tg_xxy_xxyyyyyy_0[j] - tg_xxy_xxyyyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxy_xyyyyyy_1[j];

                    tg_xxxxy_xxyyyyyz_0[j] = pb_x * tg_xxxy_xxyyyyyz_0[j] + fr * tg_xxxy_xxyyyyyz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxyyyyyz_0[j] - tg_xxy_xxyyyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxy_xyyyyyz_1[j];

                    tg_xxxxy_xxyyyyzz_0[j] = pb_x * tg_xxxy_xxyyyyzz_0[j] + fr * tg_xxxy_xxyyyyzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxyyyyzz_0[j] - tg_xxy_xxyyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxy_xyyyyzz_1[j];

                    tg_xxxxy_xxyyyzzz_0[j] = pb_x * tg_xxxy_xxyyyzzz_0[j] + fr * tg_xxxy_xxyyyzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxyyyzzz_0[j] - tg_xxy_xxyyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxy_xyyyzzz_1[j];

                    tg_xxxxy_xxyyzzzz_0[j] = pb_x * tg_xxxy_xxyyzzzz_0[j] + fr * tg_xxxy_xxyyzzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxyyzzzz_0[j] - tg_xxy_xxyyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxy_xyyzzzz_1[j];

                    tg_xxxxy_xxyzzzzz_0[j] = pb_x * tg_xxxy_xxyzzzzz_0[j] + fr * tg_xxxy_xxyzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxyzzzzz_0[j] - tg_xxy_xxyzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxy_xyzzzzz_1[j];

                    tg_xxxxy_xxzzzzzz_0[j] = pb_x * tg_xxxy_xxzzzzzz_0[j] + fr * tg_xxxy_xxzzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxzzzzzz_0[j] - tg_xxy_xxzzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxy_xzzzzzz_1[j];

                    tg_xxxxy_xyyyyyyy_0[j] = pb_x * tg_xxxy_xyyyyyyy_0[j] + fr * tg_xxxy_xyyyyyyy_1[j] + 1.5 * fl1_fx * (tg_xxy_xyyyyyyy_0[j] - tg_xxy_xyyyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxy_yyyyyyy_1[j];

                    tg_xxxxy_xyyyyyyz_0[j] = pb_x * tg_xxxy_xyyyyyyz_0[j] + fr * tg_xxxy_xyyyyyyz_1[j] + 1.5 * fl1_fx * (tg_xxy_xyyyyyyz_0[j] - tg_xxy_xyyyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxy_yyyyyyz_1[j];

                    tg_xxxxy_xyyyyyzz_0[j] = pb_x * tg_xxxy_xyyyyyzz_0[j] + fr * tg_xxxy_xyyyyyzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xyyyyyzz_0[j] - tg_xxy_xyyyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxy_yyyyyzz_1[j];

                    tg_xxxxy_xyyyyzzz_0[j] = pb_x * tg_xxxy_xyyyyzzz_0[j] + fr * tg_xxxy_xyyyyzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xyyyyzzz_0[j] - tg_xxy_xyyyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxy_yyyyzzz_1[j];

                    tg_xxxxy_xyyyzzzz_0[j] = pb_x * tg_xxxy_xyyyzzzz_0[j] + fr * tg_xxxy_xyyyzzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xyyyzzzz_0[j] - tg_xxy_xyyyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxy_yyyzzzz_1[j];

                    tg_xxxxy_xyyzzzzz_0[j] = pb_x * tg_xxxy_xyyzzzzz_0[j] + fr * tg_xxxy_xyyzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xyyzzzzz_0[j] - tg_xxy_xyyzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxy_yyzzzzz_1[j];

                    tg_xxxxy_xyzzzzzz_0[j] = pb_x * tg_xxxy_xyzzzzzz_0[j] + fr * tg_xxxy_xyzzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xyzzzzzz_0[j] - tg_xxy_xyzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxy_yzzzzzz_1[j];

                    tg_xxxxy_xzzzzzzz_0[j] = pb_x * tg_xxxy_xzzzzzzz_0[j] + fr * tg_xxxy_xzzzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xzzzzzzz_0[j] - tg_xxy_xzzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxy_zzzzzzz_1[j];

                    tg_xxxxy_yyyyyyyy_0[j] = pb_x * tg_xxxy_yyyyyyyy_0[j] + fr * tg_xxxy_yyyyyyyy_1[j] + 1.5 * fl1_fx * (tg_xxy_yyyyyyyy_0[j] - tg_xxy_yyyyyyyy_1[j] * fl1_fza);

                    tg_xxxxy_yyyyyyyz_0[j] = pb_x * tg_xxxy_yyyyyyyz_0[j] + fr * tg_xxxy_yyyyyyyz_1[j] + 1.5 * fl1_fx * (tg_xxy_yyyyyyyz_0[j] - tg_xxy_yyyyyyyz_1[j] * fl1_fza);

                    tg_xxxxy_yyyyyyzz_0[j] = pb_x * tg_xxxy_yyyyyyzz_0[j] + fr * tg_xxxy_yyyyyyzz_1[j] + 1.5 * fl1_fx * (tg_xxy_yyyyyyzz_0[j] - tg_xxy_yyyyyyzz_1[j] * fl1_fza);

                    tg_xxxxy_yyyyyzzz_0[j] = pb_x * tg_xxxy_yyyyyzzz_0[j] + fr * tg_xxxy_yyyyyzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_yyyyyzzz_0[j] - tg_xxy_yyyyyzzz_1[j] * fl1_fza);

                    tg_xxxxy_yyyyzzzz_0[j] = pb_x * tg_xxxy_yyyyzzzz_0[j] + fr * tg_xxxy_yyyyzzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_yyyyzzzz_0[j] - tg_xxy_yyyyzzzz_1[j] * fl1_fza);

                    tg_xxxxy_yyyzzzzz_0[j] = pb_x * tg_xxxy_yyyzzzzz_0[j] + fr * tg_xxxy_yyyzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_yyyzzzzz_0[j] - tg_xxy_yyyzzzzz_1[j] * fl1_fza);

                    tg_xxxxy_yyzzzzzz_0[j] = pb_x * tg_xxxy_yyzzzzzz_0[j] + fr * tg_xxxy_yyzzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_yyzzzzzz_0[j] - tg_xxy_yyzzzzzz_1[j] * fl1_fza);

                    tg_xxxxy_yzzzzzzz_0[j] = pb_x * tg_xxxy_yzzzzzzz_0[j] + fr * tg_xxxy_yzzzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_yzzzzzzz_0[j] - tg_xxy_yzzzzzzz_1[j] * fl1_fza);

                    tg_xxxxy_zzzzzzzz_0[j] = pb_x * tg_xxxy_zzzzzzzz_0[j] + fr * tg_xxxy_zzzzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_zzzzzzzz_0[j] - tg_xxy_zzzzzzzz_1[j] * fl1_fza);

                    tg_xxxxz_xxxxxxxx_0[j] = pb_x * tg_xxxz_xxxxxxxx_0[j] + fr * tg_xxxz_xxxxxxxx_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxxxxxx_0[j] - tg_xxz_xxxxxxxx_1[j] * fl1_fza) + 4.0 * fl1_fxn * tg_xxxz_xxxxxxx_1[j];

                    tg_xxxxz_xxxxxxxy_0[j] = pb_x * tg_xxxz_xxxxxxxy_0[j] + fr * tg_xxxz_xxxxxxxy_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxxxxxy_0[j] - tg_xxz_xxxxxxxy_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_xxxz_xxxxxxy_1[j];

                    tg_xxxxz_xxxxxxxz_0[j] = pb_x * tg_xxxz_xxxxxxxz_0[j] + fr * tg_xxxz_xxxxxxxz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxxxxxz_0[j] - tg_xxz_xxxxxxxz_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_xxxz_xxxxxxz_1[j];

                    tg_xxxxz_xxxxxxyy_0[j] = pb_x * tg_xxxz_xxxxxxyy_0[j] + fr * tg_xxxz_xxxxxxyy_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxxxxyy_0[j] - tg_xxz_xxxxxxyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxxz_xxxxxyy_1[j];

                    tg_xxxxz_xxxxxxyz_0[j] = pb_x * tg_xxxz_xxxxxxyz_0[j] + fr * tg_xxxz_xxxxxxyz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxxxxyz_0[j] - tg_xxz_xxxxxxyz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxxz_xxxxxyz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSL_95_190(      CMemBlock2D<double>* primBuffer,
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
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_8_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_8_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_8_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_8_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_8_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_7_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {7, -1, -1, -1}, 
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

                auto tg_xxxz_xxxxxxzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 95); 

                auto tg_xxxz_xxxxxyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 96); 

                auto tg_xxxz_xxxxxyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 97); 

                auto tg_xxxz_xxxxxyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 98); 

                auto tg_xxxz_xxxxxzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 99); 

                auto tg_xxxz_xxxxyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 100); 

                auto tg_xxxz_xxxxyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 101); 

                auto tg_xxxz_xxxxyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 102); 

                auto tg_xxxz_xxxxyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 103); 

                auto tg_xxxz_xxxxzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 104); 

                auto tg_xxxz_xxxyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 105); 

                auto tg_xxxz_xxxyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 106); 

                auto tg_xxxz_xxxyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 107); 

                auto tg_xxxz_xxxyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 108); 

                auto tg_xxxz_xxxyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 109); 

                auto tg_xxxz_xxxzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 110); 

                auto tg_xxxz_xxyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 111); 

                auto tg_xxxz_xxyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 112); 

                auto tg_xxxz_xxyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 113); 

                auto tg_xxxz_xxyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 114); 

                auto tg_xxxz_xxyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 115); 

                auto tg_xxxz_xxyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 116); 

                auto tg_xxxz_xxzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 117); 

                auto tg_xxxz_xyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 118); 

                auto tg_xxxz_xyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 119); 

                auto tg_xxxz_xyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 120); 

                auto tg_xxxz_xyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 121); 

                auto tg_xxxz_xyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 122); 

                auto tg_xxxz_xyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 123); 

                auto tg_xxxz_xyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 124); 

                auto tg_xxxz_xzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 125); 

                auto tg_xxxz_yyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 126); 

                auto tg_xxxz_yyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 127); 

                auto tg_xxxz_yyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 128); 

                auto tg_xxxz_yyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 129); 

                auto tg_xxxz_yyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 130); 

                auto tg_xxxz_yyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 131); 

                auto tg_xxxz_yyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 132); 

                auto tg_xxxz_yzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 133); 

                auto tg_xxxz_zzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 134); 

                auto tg_xxyy_xxxxxxxx_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 135); 

                auto tg_xxyy_xxxxxxxy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 136); 

                auto tg_xxyy_xxxxxxxz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 137); 

                auto tg_xxyy_xxxxxxyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 138); 

                auto tg_xxyy_xxxxxxyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 139); 

                auto tg_xxyy_xxxxxxzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 140); 

                auto tg_xxyy_xxxxxyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 141); 

                auto tg_xxyy_xxxxxyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 142); 

                auto tg_xxyy_xxxxxyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 143); 

                auto tg_xxyy_xxxxxzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 144); 

                auto tg_xxyy_xxxxyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 145); 

                auto tg_xxyy_xxxxyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 146); 

                auto tg_xxyy_xxxxyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 147); 

                auto tg_xxyy_xxxxyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 148); 

                auto tg_xxyy_xxxxzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 149); 

                auto tg_xxyy_xxxyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 150); 

                auto tg_xxyy_xxxyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 151); 

                auto tg_xxyy_xxxyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 152); 

                auto tg_xxyy_xxxyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 153); 

                auto tg_xxyy_xxxyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 154); 

                auto tg_xxyy_xxxzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 155); 

                auto tg_xxyy_xxyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 156); 

                auto tg_xxyy_xxyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 157); 

                auto tg_xxyy_xxyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 158); 

                auto tg_xxyy_xxyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 159); 

                auto tg_xxyy_xxyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 160); 

                auto tg_xxyy_xxyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 161); 

                auto tg_xxyy_xxzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 162); 

                auto tg_xxyy_xyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 163); 

                auto tg_xxyy_xyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 164); 

                auto tg_xxyy_xyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 165); 

                auto tg_xxyy_xyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 166); 

                auto tg_xxyy_xyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 167); 

                auto tg_xxyy_xyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 168); 

                auto tg_xxyy_xyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 169); 

                auto tg_xxyy_xzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 170); 

                auto tg_xxyy_yyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 171); 

                auto tg_xxyy_yyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 172); 

                auto tg_xxyy_yyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 173); 

                auto tg_xxyy_yyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 174); 

                auto tg_xxyy_yyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 175); 

                auto tg_xxyy_yyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 176); 

                auto tg_xxyy_yyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 177); 

                auto tg_xxyy_yzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 178); 

                auto tg_xxyy_zzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 179); 

                auto tg_xxyz_xxxxxxxx_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 180); 

                auto tg_xxyz_xxxxxxxy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 181); 

                auto tg_xxyz_xxxxxxxz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 182); 

                auto tg_xxyz_xxxxxxyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 183); 

                auto tg_xxyz_xxxxxxyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 184); 

                auto tg_xxyz_xxxxxxzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 185); 

                auto tg_xxyz_xxxxxyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 186); 

                auto tg_xxyz_xxxxxyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 187); 

                auto tg_xxyz_xxxxxyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 188); 

                auto tg_xxyz_xxxxxzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 189); 

                auto tg_xxxz_xxxxxxzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 95); 

                auto tg_xxxz_xxxxxyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 96); 

                auto tg_xxxz_xxxxxyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 97); 

                auto tg_xxxz_xxxxxyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 98); 

                auto tg_xxxz_xxxxxzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 99); 

                auto tg_xxxz_xxxxyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 100); 

                auto tg_xxxz_xxxxyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 101); 

                auto tg_xxxz_xxxxyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 102); 

                auto tg_xxxz_xxxxyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 103); 

                auto tg_xxxz_xxxxzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 104); 

                auto tg_xxxz_xxxyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 105); 

                auto tg_xxxz_xxxyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 106); 

                auto tg_xxxz_xxxyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 107); 

                auto tg_xxxz_xxxyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 108); 

                auto tg_xxxz_xxxyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 109); 

                auto tg_xxxz_xxxzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 110); 

                auto tg_xxxz_xxyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 111); 

                auto tg_xxxz_xxyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 112); 

                auto tg_xxxz_xxyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 113); 

                auto tg_xxxz_xxyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 114); 

                auto tg_xxxz_xxyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 115); 

                auto tg_xxxz_xxyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 116); 

                auto tg_xxxz_xxzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 117); 

                auto tg_xxxz_xyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 118); 

                auto tg_xxxz_xyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 119); 

                auto tg_xxxz_xyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 120); 

                auto tg_xxxz_xyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 121); 

                auto tg_xxxz_xyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 122); 

                auto tg_xxxz_xyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 123); 

                auto tg_xxxz_xyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 124); 

                auto tg_xxxz_xzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 125); 

                auto tg_xxxz_yyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 126); 

                auto tg_xxxz_yyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 127); 

                auto tg_xxxz_yyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 128); 

                auto tg_xxxz_yyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 129); 

                auto tg_xxxz_yyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 130); 

                auto tg_xxxz_yyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 131); 

                auto tg_xxxz_yyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 132); 

                auto tg_xxxz_yzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 133); 

                auto tg_xxxz_zzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 134); 

                auto tg_xxyy_xxxxxxxx_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 135); 

                auto tg_xxyy_xxxxxxxy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 136); 

                auto tg_xxyy_xxxxxxxz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 137); 

                auto tg_xxyy_xxxxxxyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 138); 

                auto tg_xxyy_xxxxxxyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 139); 

                auto tg_xxyy_xxxxxxzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 140); 

                auto tg_xxyy_xxxxxyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 141); 

                auto tg_xxyy_xxxxxyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 142); 

                auto tg_xxyy_xxxxxyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 143); 

                auto tg_xxyy_xxxxxzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 144); 

                auto tg_xxyy_xxxxyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 145); 

                auto tg_xxyy_xxxxyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 146); 

                auto tg_xxyy_xxxxyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 147); 

                auto tg_xxyy_xxxxyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 148); 

                auto tg_xxyy_xxxxzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 149); 

                auto tg_xxyy_xxxyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 150); 

                auto tg_xxyy_xxxyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 151); 

                auto tg_xxyy_xxxyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 152); 

                auto tg_xxyy_xxxyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 153); 

                auto tg_xxyy_xxxyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 154); 

                auto tg_xxyy_xxxzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 155); 

                auto tg_xxyy_xxyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 156); 

                auto tg_xxyy_xxyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 157); 

                auto tg_xxyy_xxyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 158); 

                auto tg_xxyy_xxyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 159); 

                auto tg_xxyy_xxyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 160); 

                auto tg_xxyy_xxyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 161); 

                auto tg_xxyy_xxzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 162); 

                auto tg_xxyy_xyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 163); 

                auto tg_xxyy_xyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 164); 

                auto tg_xxyy_xyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 165); 

                auto tg_xxyy_xyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 166); 

                auto tg_xxyy_xyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 167); 

                auto tg_xxyy_xyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 168); 

                auto tg_xxyy_xyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 169); 

                auto tg_xxyy_xzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 170); 

                auto tg_xxyy_yyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 171); 

                auto tg_xxyy_yyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 172); 

                auto tg_xxyy_yyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 173); 

                auto tg_xxyy_yyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 174); 

                auto tg_xxyy_yyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 175); 

                auto tg_xxyy_yyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 176); 

                auto tg_xxyy_yyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 177); 

                auto tg_xxyy_yzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 178); 

                auto tg_xxyy_zzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 179); 

                auto tg_xxyz_xxxxxxxx_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 180); 

                auto tg_xxyz_xxxxxxxy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 181); 

                auto tg_xxyz_xxxxxxxz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 182); 

                auto tg_xxyz_xxxxxxyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 183); 

                auto tg_xxyz_xxxxxxyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 184); 

                auto tg_xxyz_xxxxxxzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 185); 

                auto tg_xxyz_xxxxxyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 186); 

                auto tg_xxyz_xxxxxyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 187); 

                auto tg_xxyz_xxxxxyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 188); 

                auto tg_xxyz_xxxxxzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 189); 

                auto tg_xxz_xxxxxxzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 95); 

                auto tg_xxz_xxxxxyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 96); 

                auto tg_xxz_xxxxxyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 97); 

                auto tg_xxz_xxxxxyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 98); 

                auto tg_xxz_xxxxxzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 99); 

                auto tg_xxz_xxxxyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 100); 

                auto tg_xxz_xxxxyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 101); 

                auto tg_xxz_xxxxyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 102); 

                auto tg_xxz_xxxxyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 103); 

                auto tg_xxz_xxxxzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 104); 

                auto tg_xxz_xxxyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 105); 

                auto tg_xxz_xxxyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 106); 

                auto tg_xxz_xxxyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 107); 

                auto tg_xxz_xxxyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 108); 

                auto tg_xxz_xxxyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 109); 

                auto tg_xxz_xxxzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 110); 

                auto tg_xxz_xxyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 111); 

                auto tg_xxz_xxyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 112); 

                auto tg_xxz_xxyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 113); 

                auto tg_xxz_xxyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 114); 

                auto tg_xxz_xxyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 115); 

                auto tg_xxz_xxyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 116); 

                auto tg_xxz_xxzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 117); 

                auto tg_xxz_xyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 118); 

                auto tg_xxz_xyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 119); 

                auto tg_xxz_xyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 120); 

                auto tg_xxz_xyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 121); 

                auto tg_xxz_xyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 122); 

                auto tg_xxz_xyyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 123); 

                auto tg_xxz_xyzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 124); 

                auto tg_xxz_xzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 125); 

                auto tg_xxz_yyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 126); 

                auto tg_xxz_yyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 127); 

                auto tg_xxz_yyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 128); 

                auto tg_xxz_yyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 129); 

                auto tg_xxz_yyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 130); 

                auto tg_xxz_yyyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 131); 

                auto tg_xxz_yyzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 132); 

                auto tg_xxz_yzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 133); 

                auto tg_xxz_zzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 134); 

                auto tg_xyy_xxxxxxxx_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 135); 

                auto tg_xyy_xxxxxxxy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 136); 

                auto tg_xyy_xxxxxxxz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 137); 

                auto tg_xyy_xxxxxxyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 138); 

                auto tg_xyy_xxxxxxyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 139); 

                auto tg_xyy_xxxxxxzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 140); 

                auto tg_xyy_xxxxxyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 141); 

                auto tg_xyy_xxxxxyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 142); 

                auto tg_xyy_xxxxxyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 143); 

                auto tg_xyy_xxxxxzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 144); 

                auto tg_xyy_xxxxyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 145); 

                auto tg_xyy_xxxxyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 146); 

                auto tg_xyy_xxxxyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 147); 

                auto tg_xyy_xxxxyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 148); 

                auto tg_xyy_xxxxzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 149); 

                auto tg_xyy_xxxyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 150); 

                auto tg_xyy_xxxyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 151); 

                auto tg_xyy_xxxyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 152); 

                auto tg_xyy_xxxyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 153); 

                auto tg_xyy_xxxyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 154); 

                auto tg_xyy_xxxzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 155); 

                auto tg_xyy_xxyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 156); 

                auto tg_xyy_xxyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 157); 

                auto tg_xyy_xxyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 158); 

                auto tg_xyy_xxyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 159); 

                auto tg_xyy_xxyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 160); 

                auto tg_xyy_xxyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 161); 

                auto tg_xyy_xxzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 162); 

                auto tg_xyy_xyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 163); 

                auto tg_xyy_xyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 164); 

                auto tg_xyy_xyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 165); 

                auto tg_xyy_xyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 166); 

                auto tg_xyy_xyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 167); 

                auto tg_xyy_xyyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 168); 

                auto tg_xyy_xyzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 169); 

                auto tg_xyy_xzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 170); 

                auto tg_xyy_yyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 171); 

                auto tg_xyy_yyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 172); 

                auto tg_xyy_yyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 173); 

                auto tg_xyy_yyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 174); 

                auto tg_xyy_yyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 175); 

                auto tg_xyy_yyyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 176); 

                auto tg_xyy_yyzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 177); 

                auto tg_xyy_yzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 178); 

                auto tg_xyy_zzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 179); 

                auto tg_xyz_xxxxxxxx_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 180); 

                auto tg_xyz_xxxxxxxy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 181); 

                auto tg_xyz_xxxxxxxz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 182); 

                auto tg_xyz_xxxxxxyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 183); 

                auto tg_xyz_xxxxxxyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 184); 

                auto tg_xyz_xxxxxxzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 185); 

                auto tg_xyz_xxxxxyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 186); 

                auto tg_xyz_xxxxxyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 187); 

                auto tg_xyz_xxxxxyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 188); 

                auto tg_xyz_xxxxxzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 189); 

                auto tg_xxz_xxxxxxzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 95); 

                auto tg_xxz_xxxxxyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 96); 

                auto tg_xxz_xxxxxyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 97); 

                auto tg_xxz_xxxxxyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 98); 

                auto tg_xxz_xxxxxzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 99); 

                auto tg_xxz_xxxxyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 100); 

                auto tg_xxz_xxxxyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 101); 

                auto tg_xxz_xxxxyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 102); 

                auto tg_xxz_xxxxyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 103); 

                auto tg_xxz_xxxxzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 104); 

                auto tg_xxz_xxxyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 105); 

                auto tg_xxz_xxxyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 106); 

                auto tg_xxz_xxxyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 107); 

                auto tg_xxz_xxxyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 108); 

                auto tg_xxz_xxxyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 109); 

                auto tg_xxz_xxxzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 110); 

                auto tg_xxz_xxyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 111); 

                auto tg_xxz_xxyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 112); 

                auto tg_xxz_xxyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 113); 

                auto tg_xxz_xxyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 114); 

                auto tg_xxz_xxyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 115); 

                auto tg_xxz_xxyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 116); 

                auto tg_xxz_xxzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 117); 

                auto tg_xxz_xyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 118); 

                auto tg_xxz_xyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 119); 

                auto tg_xxz_xyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 120); 

                auto tg_xxz_xyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 121); 

                auto tg_xxz_xyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 122); 

                auto tg_xxz_xyyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 123); 

                auto tg_xxz_xyzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 124); 

                auto tg_xxz_xzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 125); 

                auto tg_xxz_yyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 126); 

                auto tg_xxz_yyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 127); 

                auto tg_xxz_yyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 128); 

                auto tg_xxz_yyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 129); 

                auto tg_xxz_yyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 130); 

                auto tg_xxz_yyyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 131); 

                auto tg_xxz_yyzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 132); 

                auto tg_xxz_yzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 133); 

                auto tg_xxz_zzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 134); 

                auto tg_xyy_xxxxxxxx_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 135); 

                auto tg_xyy_xxxxxxxy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 136); 

                auto tg_xyy_xxxxxxxz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 137); 

                auto tg_xyy_xxxxxxyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 138); 

                auto tg_xyy_xxxxxxyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 139); 

                auto tg_xyy_xxxxxxzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 140); 

                auto tg_xyy_xxxxxyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 141); 

                auto tg_xyy_xxxxxyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 142); 

                auto tg_xyy_xxxxxyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 143); 

                auto tg_xyy_xxxxxzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 144); 

                auto tg_xyy_xxxxyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 145); 

                auto tg_xyy_xxxxyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 146); 

                auto tg_xyy_xxxxyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 147); 

                auto tg_xyy_xxxxyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 148); 

                auto tg_xyy_xxxxzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 149); 

                auto tg_xyy_xxxyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 150); 

                auto tg_xyy_xxxyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 151); 

                auto tg_xyy_xxxyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 152); 

                auto tg_xyy_xxxyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 153); 

                auto tg_xyy_xxxyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 154); 

                auto tg_xyy_xxxzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 155); 

                auto tg_xyy_xxyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 156); 

                auto tg_xyy_xxyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 157); 

                auto tg_xyy_xxyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 158); 

                auto tg_xyy_xxyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 159); 

                auto tg_xyy_xxyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 160); 

                auto tg_xyy_xxyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 161); 

                auto tg_xyy_xxzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 162); 

                auto tg_xyy_xyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 163); 

                auto tg_xyy_xyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 164); 

                auto tg_xyy_xyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 165); 

                auto tg_xyy_xyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 166); 

                auto tg_xyy_xyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 167); 

                auto tg_xyy_xyyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 168); 

                auto tg_xyy_xyzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 169); 

                auto tg_xyy_xzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 170); 

                auto tg_xyy_yyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 171); 

                auto tg_xyy_yyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 172); 

                auto tg_xyy_yyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 173); 

                auto tg_xyy_yyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 174); 

                auto tg_xyy_yyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 175); 

                auto tg_xyy_yyyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 176); 

                auto tg_xyy_yyzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 177); 

                auto tg_xyy_yzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 178); 

                auto tg_xyy_zzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 179); 

                auto tg_xyz_xxxxxxxx_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 180); 

                auto tg_xyz_xxxxxxxy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 181); 

                auto tg_xyz_xxxxxxxz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 182); 

                auto tg_xyz_xxxxxxyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 183); 

                auto tg_xyz_xxxxxxyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 184); 

                auto tg_xyz_xxxxxxzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 185); 

                auto tg_xyz_xxxxxyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 186); 

                auto tg_xyz_xxxxxyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 187); 

                auto tg_xyz_xxxxxyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 188); 

                auto tg_xyz_xxxxxzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 189); 

                auto tg_xxxz_xxxxxzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 77); 

                auto tg_xxxz_xxxxyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 78); 

                auto tg_xxxz_xxxxyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 79); 

                auto tg_xxxz_xxxxyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 80); 

                auto tg_xxxz_xxxxzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 81); 

                auto tg_xxxz_xxxyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 82); 

                auto tg_xxxz_xxxyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 83); 

                auto tg_xxxz_xxxyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 84); 

                auto tg_xxxz_xxxyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 85); 

                auto tg_xxxz_xxxzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 86); 

                auto tg_xxxz_xxyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 87); 

                auto tg_xxxz_xxyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 88); 

                auto tg_xxxz_xxyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 89); 

                auto tg_xxxz_xxyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 90); 

                auto tg_xxxz_xxyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 91); 

                auto tg_xxxz_xxzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 92); 

                auto tg_xxxz_xyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 93); 

                auto tg_xxxz_xyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 94); 

                auto tg_xxxz_xyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 95); 

                auto tg_xxxz_xyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 96); 

                auto tg_xxxz_xyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 97); 

                auto tg_xxxz_xyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 98); 

                auto tg_xxxz_xzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 99); 

                auto tg_xxxz_yyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 100); 

                auto tg_xxxz_yyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 101); 

                auto tg_xxxz_yyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 102); 

                auto tg_xxxz_yyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 103); 

                auto tg_xxxz_yyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 104); 

                auto tg_xxxz_yyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 105); 

                auto tg_xxxz_yzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 106); 

                auto tg_xxxz_zzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 107); 

                auto tg_xxyy_xxxxxxx_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 108); 

                auto tg_xxyy_xxxxxxy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 109); 

                auto tg_xxyy_xxxxxxz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 110); 

                auto tg_xxyy_xxxxxyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 111); 

                auto tg_xxyy_xxxxxyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 112); 

                auto tg_xxyy_xxxxxzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 113); 

                auto tg_xxyy_xxxxyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 114); 

                auto tg_xxyy_xxxxyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 115); 

                auto tg_xxyy_xxxxyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 116); 

                auto tg_xxyy_xxxxzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 117); 

                auto tg_xxyy_xxxyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 118); 

                auto tg_xxyy_xxxyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 119); 

                auto tg_xxyy_xxxyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 120); 

                auto tg_xxyy_xxxyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 121); 

                auto tg_xxyy_xxxzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 122); 

                auto tg_xxyy_xxyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 123); 

                auto tg_xxyy_xxyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 124); 

                auto tg_xxyy_xxyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 125); 

                auto tg_xxyy_xxyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 126); 

                auto tg_xxyy_xxyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 127); 

                auto tg_xxyy_xxzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 128); 

                auto tg_xxyy_xyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 129); 

                auto tg_xxyy_xyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 130); 

                auto tg_xxyy_xyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 131); 

                auto tg_xxyy_xyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 132); 

                auto tg_xxyy_xyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 133); 

                auto tg_xxyy_xyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 134); 

                auto tg_xxyy_xzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 135); 

                auto tg_xxyy_yyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 136); 

                auto tg_xxyy_yyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 137); 

                auto tg_xxyy_yyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 138); 

                auto tg_xxyy_yyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 139); 

                auto tg_xxyy_yyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 140); 

                auto tg_xxyy_yyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 141); 

                auto tg_xxyy_yzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 142); 

                auto tg_xxyy_zzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 143); 

                auto tg_xxyz_xxxxxxx_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 144); 

                auto tg_xxyz_xxxxxxy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 145); 

                auto tg_xxyz_xxxxxxz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 146); 

                auto tg_xxyz_xxxxxyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 147); 

                auto tg_xxyz_xxxxxyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 148); 

                auto tg_xxyz_xxxxxzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 149); 

                auto tg_xxyz_xxxxyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 150); 

                auto tg_xxyz_xxxxyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 151); 

                auto tg_xxyz_xxxxyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 152); 

                auto tg_xxyz_xxxxzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 153); 

                // set up pointers to integrals

                auto tg_xxxxz_xxxxxxzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 95); 

                auto tg_xxxxz_xxxxxyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 96); 

                auto tg_xxxxz_xxxxxyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 97); 

                auto tg_xxxxz_xxxxxyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 98); 

                auto tg_xxxxz_xxxxxzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 99); 

                auto tg_xxxxz_xxxxyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 100); 

                auto tg_xxxxz_xxxxyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 101); 

                auto tg_xxxxz_xxxxyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 102); 

                auto tg_xxxxz_xxxxyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 103); 

                auto tg_xxxxz_xxxxzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 104); 

                auto tg_xxxxz_xxxyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 105); 

                auto tg_xxxxz_xxxyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 106); 

                auto tg_xxxxz_xxxyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 107); 

                auto tg_xxxxz_xxxyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 108); 

                auto tg_xxxxz_xxxyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 109); 

                auto tg_xxxxz_xxxzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 110); 

                auto tg_xxxxz_xxyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 111); 

                auto tg_xxxxz_xxyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 112); 

                auto tg_xxxxz_xxyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 113); 

                auto tg_xxxxz_xxyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 114); 

                auto tg_xxxxz_xxyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 115); 

                auto tg_xxxxz_xxyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 116); 

                auto tg_xxxxz_xxzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 117); 

                auto tg_xxxxz_xyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 118); 

                auto tg_xxxxz_xyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 119); 

                auto tg_xxxxz_xyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 120); 

                auto tg_xxxxz_xyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 121); 

                auto tg_xxxxz_xyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 122); 

                auto tg_xxxxz_xyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 123); 

                auto tg_xxxxz_xyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 124); 

                auto tg_xxxxz_xzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 125); 

                auto tg_xxxxz_yyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 126); 

                auto tg_xxxxz_yyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 127); 

                auto tg_xxxxz_yyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 128); 

                auto tg_xxxxz_yyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 129); 

                auto tg_xxxxz_yyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 130); 

                auto tg_xxxxz_yyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 131); 

                auto tg_xxxxz_yyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 132); 

                auto tg_xxxxz_yzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 133); 

                auto tg_xxxxz_zzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 134); 

                auto tg_xxxyy_xxxxxxxx_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 135); 

                auto tg_xxxyy_xxxxxxxy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 136); 

                auto tg_xxxyy_xxxxxxxz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 137); 

                auto tg_xxxyy_xxxxxxyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 138); 

                auto tg_xxxyy_xxxxxxyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 139); 

                auto tg_xxxyy_xxxxxxzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 140); 

                auto tg_xxxyy_xxxxxyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 141); 

                auto tg_xxxyy_xxxxxyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 142); 

                auto tg_xxxyy_xxxxxyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 143); 

                auto tg_xxxyy_xxxxxzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 144); 

                auto tg_xxxyy_xxxxyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 145); 

                auto tg_xxxyy_xxxxyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 146); 

                auto tg_xxxyy_xxxxyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 147); 

                auto tg_xxxyy_xxxxyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 148); 

                auto tg_xxxyy_xxxxzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 149); 

                auto tg_xxxyy_xxxyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 150); 

                auto tg_xxxyy_xxxyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 151); 

                auto tg_xxxyy_xxxyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 152); 

                auto tg_xxxyy_xxxyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 153); 

                auto tg_xxxyy_xxxyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 154); 

                auto tg_xxxyy_xxxzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 155); 

                auto tg_xxxyy_xxyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 156); 

                auto tg_xxxyy_xxyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 157); 

                auto tg_xxxyy_xxyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 158); 

                auto tg_xxxyy_xxyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 159); 

                auto tg_xxxyy_xxyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 160); 

                auto tg_xxxyy_xxyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 161); 

                auto tg_xxxyy_xxzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 162); 

                auto tg_xxxyy_xyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 163); 

                auto tg_xxxyy_xyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 164); 

                auto tg_xxxyy_xyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 165); 

                auto tg_xxxyy_xyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 166); 

                auto tg_xxxyy_xyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 167); 

                auto tg_xxxyy_xyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 168); 

                auto tg_xxxyy_xyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 169); 

                auto tg_xxxyy_xzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 170); 

                auto tg_xxxyy_yyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 171); 

                auto tg_xxxyy_yyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 172); 

                auto tg_xxxyy_yyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 173); 

                auto tg_xxxyy_yyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 174); 

                auto tg_xxxyy_yyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 175); 

                auto tg_xxxyy_yyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 176); 

                auto tg_xxxyy_yyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 177); 

                auto tg_xxxyy_yzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 178); 

                auto tg_xxxyy_zzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 179); 

                auto tg_xxxyz_xxxxxxxx_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 180); 

                auto tg_xxxyz_xxxxxxxy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 181); 

                auto tg_xxxyz_xxxxxxxz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 182); 

                auto tg_xxxyz_xxxxxxyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 183); 

                auto tg_xxxyz_xxxxxxyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 184); 

                auto tg_xxxyz_xxxxxxzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 185); 

                auto tg_xxxyz_xxxxxyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 186); 

                auto tg_xxxyz_xxxxxyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 187); 

                auto tg_xxxyz_xxxxxyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 188); 

                auto tg_xxxyz_xxxxxzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 189); 

                // Batch of Integrals (95,190)

                #pragma omp simd aligned(fxn, fza, tg_xxxxz_xxxxxxzz_0, tg_xxxxz_xxxxxyyy_0, \
                                         tg_xxxxz_xxxxxyyz_0, tg_xxxxz_xxxxxyzz_0, tg_xxxxz_xxxxxzzz_0, tg_xxxxz_xxxxyyyy_0, \
                                         tg_xxxxz_xxxxyyyz_0, tg_xxxxz_xxxxyyzz_0, tg_xxxxz_xxxxyzzz_0, tg_xxxxz_xxxxzzzz_0, \
                                         tg_xxxxz_xxxyyyyy_0, tg_xxxxz_xxxyyyyz_0, tg_xxxxz_xxxyyyzz_0, tg_xxxxz_xxxyyzzz_0, \
                                         tg_xxxxz_xxxyzzzz_0, tg_xxxxz_xxxzzzzz_0, tg_xxxxz_xxyyyyyy_0, tg_xxxxz_xxyyyyyz_0, \
                                         tg_xxxxz_xxyyyyzz_0, tg_xxxxz_xxyyyzzz_0, tg_xxxxz_xxyyzzzz_0, tg_xxxxz_xxyzzzzz_0, \
                                         tg_xxxxz_xxzzzzzz_0, tg_xxxxz_xyyyyyyy_0, tg_xxxxz_xyyyyyyz_0, tg_xxxxz_xyyyyyzz_0, \
                                         tg_xxxxz_xyyyyzzz_0, tg_xxxxz_xyyyzzzz_0, tg_xxxxz_xyyzzzzz_0, tg_xxxxz_xyzzzzzz_0, \
                                         tg_xxxxz_xzzzzzzz_0, tg_xxxxz_yyyyyyyy_0, tg_xxxxz_yyyyyyyz_0, tg_xxxxz_yyyyyyzz_0, \
                                         tg_xxxxz_yyyyyzzz_0, tg_xxxxz_yyyyzzzz_0, tg_xxxxz_yyyzzzzz_0, tg_xxxxz_yyzzzzzz_0, \
                                         tg_xxxxz_yzzzzzzz_0, tg_xxxxz_zzzzzzzz_0, tg_xxxyy_xxxxxxxx_0, tg_xxxyy_xxxxxxxy_0, \
                                         tg_xxxyy_xxxxxxxz_0, tg_xxxyy_xxxxxxyy_0, tg_xxxyy_xxxxxxyz_0, tg_xxxyy_xxxxxxzz_0, \
                                         tg_xxxyy_xxxxxyyy_0, tg_xxxyy_xxxxxyyz_0, tg_xxxyy_xxxxxyzz_0, tg_xxxyy_xxxxxzzz_0, \
                                         tg_xxxyy_xxxxyyyy_0, tg_xxxyy_xxxxyyyz_0, tg_xxxyy_xxxxyyzz_0, tg_xxxyy_xxxxyzzz_0, \
                                         tg_xxxyy_xxxxzzzz_0, tg_xxxyy_xxxyyyyy_0, tg_xxxyy_xxxyyyyz_0, tg_xxxyy_xxxyyyzz_0, \
                                         tg_xxxyy_xxxyyzzz_0, tg_xxxyy_xxxyzzzz_0, tg_xxxyy_xxxzzzzz_0, tg_xxxyy_xxyyyyyy_0, \
                                         tg_xxxyy_xxyyyyyz_0, tg_xxxyy_xxyyyyzz_0, tg_xxxyy_xxyyyzzz_0, tg_xxxyy_xxyyzzzz_0, \
                                         tg_xxxyy_xxyzzzzz_0, tg_xxxyy_xxzzzzzz_0, tg_xxxyy_xyyyyyyy_0, tg_xxxyy_xyyyyyyz_0, \
                                         tg_xxxyy_xyyyyyzz_0, tg_xxxyy_xyyyyzzz_0, tg_xxxyy_xyyyzzzz_0, tg_xxxyy_xyyzzzzz_0, \
                                         tg_xxxyy_xyzzzzzz_0, tg_xxxyy_xzzzzzzz_0, tg_xxxyy_yyyyyyyy_0, tg_xxxyy_yyyyyyyz_0, \
                                         tg_xxxyy_yyyyyyzz_0, tg_xxxyy_yyyyyzzz_0, tg_xxxyy_yyyyzzzz_0, tg_xxxyy_yyyzzzzz_0, \
                                         tg_xxxyy_yyzzzzzz_0, tg_xxxyy_yzzzzzzz_0, tg_xxxyy_zzzzzzzz_0, tg_xxxyz_xxxxxxxx_0, \
                                         tg_xxxyz_xxxxxxxy_0, tg_xxxyz_xxxxxxxz_0, tg_xxxyz_xxxxxxyy_0, tg_xxxyz_xxxxxxyz_0, \
                                         tg_xxxyz_xxxxxxzz_0, tg_xxxyz_xxxxxyyy_0, tg_xxxyz_xxxxxyyz_0, tg_xxxyz_xxxxxyzz_0, \
                                         tg_xxxyz_xxxxxzzz_0, tg_xxxz_xxxxxxzz_0, tg_xxxz_xxxxxxzz_1, tg_xxxz_xxxxxyyy_0, \
                                         tg_xxxz_xxxxxyyy_1, tg_xxxz_xxxxxyyz_0, tg_xxxz_xxxxxyyz_1, tg_xxxz_xxxxxyzz_0, \
                                         tg_xxxz_xxxxxyzz_1, tg_xxxz_xxxxxzz_1, tg_xxxz_xxxxxzzz_0, tg_xxxz_xxxxxzzz_1, \
                                         tg_xxxz_xxxxyyy_1, tg_xxxz_xxxxyyyy_0, tg_xxxz_xxxxyyyy_1, tg_xxxz_xxxxyyyz_0, \
                                         tg_xxxz_xxxxyyyz_1, tg_xxxz_xxxxyyz_1, tg_xxxz_xxxxyyzz_0, tg_xxxz_xxxxyyzz_1, \
                                         tg_xxxz_xxxxyzz_1, tg_xxxz_xxxxyzzz_0, tg_xxxz_xxxxyzzz_1, tg_xxxz_xxxxzzz_1, \
                                         tg_xxxz_xxxxzzzz_0, tg_xxxz_xxxxzzzz_1, tg_xxxz_xxxyyyy_1, tg_xxxz_xxxyyyyy_0, \
                                         tg_xxxz_xxxyyyyy_1, tg_xxxz_xxxyyyyz_0, tg_xxxz_xxxyyyyz_1, tg_xxxz_xxxyyyz_1, \
                                         tg_xxxz_xxxyyyzz_0, tg_xxxz_xxxyyyzz_1, tg_xxxz_xxxyyzz_1, tg_xxxz_xxxyyzzz_0, \
                                         tg_xxxz_xxxyyzzz_1, tg_xxxz_xxxyzzz_1, tg_xxxz_xxxyzzzz_0, tg_xxxz_xxxyzzzz_1, \
                                         tg_xxxz_xxxzzzz_1, tg_xxxz_xxxzzzzz_0, tg_xxxz_xxxzzzzz_1, tg_xxxz_xxyyyyy_1, \
                                         tg_xxxz_xxyyyyyy_0, tg_xxxz_xxyyyyyy_1, tg_xxxz_xxyyyyyz_0, tg_xxxz_xxyyyyyz_1, \
                                         tg_xxxz_xxyyyyz_1, tg_xxxz_xxyyyyzz_0, tg_xxxz_xxyyyyzz_1, tg_xxxz_xxyyyzz_1, \
                                         tg_xxxz_xxyyyzzz_0, tg_xxxz_xxyyyzzz_1, tg_xxxz_xxyyzzz_1, tg_xxxz_xxyyzzzz_0, \
                                         tg_xxxz_xxyyzzzz_1, tg_xxxz_xxyzzzz_1, tg_xxxz_xxyzzzzz_0, tg_xxxz_xxyzzzzz_1, \
                                         tg_xxxz_xxzzzzz_1, tg_xxxz_xxzzzzzz_0, tg_xxxz_xxzzzzzz_1, tg_xxxz_xyyyyyy_1, \
                                         tg_xxxz_xyyyyyyy_0, tg_xxxz_xyyyyyyy_1, tg_xxxz_xyyyyyyz_0, tg_xxxz_xyyyyyyz_1, \
                                         tg_xxxz_xyyyyyz_1, tg_xxxz_xyyyyyzz_0, tg_xxxz_xyyyyyzz_1, tg_xxxz_xyyyyzz_1, \
                                         tg_xxxz_xyyyyzzz_0, tg_xxxz_xyyyyzzz_1, tg_xxxz_xyyyzzz_1, tg_xxxz_xyyyzzzz_0, \
                                         tg_xxxz_xyyyzzzz_1, tg_xxxz_xyyzzzz_1, tg_xxxz_xyyzzzzz_0, tg_xxxz_xyyzzzzz_1, \
                                         tg_xxxz_xyzzzzz_1, tg_xxxz_xyzzzzzz_0, tg_xxxz_xyzzzzzz_1, tg_xxxz_xzzzzzz_1, \
                                         tg_xxxz_xzzzzzzz_0, tg_xxxz_xzzzzzzz_1, tg_xxxz_yyyyyyy_1, tg_xxxz_yyyyyyyy_0, \
                                         tg_xxxz_yyyyyyyy_1, tg_xxxz_yyyyyyyz_0, tg_xxxz_yyyyyyyz_1, tg_xxxz_yyyyyyz_1, \
                                         tg_xxxz_yyyyyyzz_0, tg_xxxz_yyyyyyzz_1, tg_xxxz_yyyyyzz_1, tg_xxxz_yyyyyzzz_0, \
                                         tg_xxxz_yyyyyzzz_1, tg_xxxz_yyyyzzz_1, tg_xxxz_yyyyzzzz_0, tg_xxxz_yyyyzzzz_1, \
                                         tg_xxxz_yyyzzzz_1, tg_xxxz_yyyzzzzz_0, tg_xxxz_yyyzzzzz_1, tg_xxxz_yyzzzzz_1, \
                                         tg_xxxz_yyzzzzzz_0, tg_xxxz_yyzzzzzz_1, tg_xxxz_yzzzzzz_1, tg_xxxz_yzzzzzzz_0, \
                                         tg_xxxz_yzzzzzzz_1, tg_xxxz_zzzzzzz_1, tg_xxxz_zzzzzzzz_0, tg_xxxz_zzzzzzzz_1, \
                                         tg_xxyy_xxxxxxx_1, tg_xxyy_xxxxxxxx_0, tg_xxyy_xxxxxxxx_1, tg_xxyy_xxxxxxxy_0, \
                                         tg_xxyy_xxxxxxxy_1, tg_xxyy_xxxxxxxz_0, tg_xxyy_xxxxxxxz_1, tg_xxyy_xxxxxxy_1, \
                                         tg_xxyy_xxxxxxyy_0, tg_xxyy_xxxxxxyy_1, tg_xxyy_xxxxxxyz_0, tg_xxyy_xxxxxxyz_1, \
                                         tg_xxyy_xxxxxxz_1, tg_xxyy_xxxxxxzz_0, tg_xxyy_xxxxxxzz_1, tg_xxyy_xxxxxyy_1, \
                                         tg_xxyy_xxxxxyyy_0, tg_xxyy_xxxxxyyy_1, tg_xxyy_xxxxxyyz_0, tg_xxyy_xxxxxyyz_1, \
                                         tg_xxyy_xxxxxyz_1, tg_xxyy_xxxxxyzz_0, tg_xxyy_xxxxxyzz_1, tg_xxyy_xxxxxzz_1, \
                                         tg_xxyy_xxxxxzzz_0, tg_xxyy_xxxxxzzz_1, tg_xxyy_xxxxyyy_1, tg_xxyy_xxxxyyyy_0, \
                                         tg_xxyy_xxxxyyyy_1, tg_xxyy_xxxxyyyz_0, tg_xxyy_xxxxyyyz_1, tg_xxyy_xxxxyyz_1, \
                                         tg_xxyy_xxxxyyzz_0, tg_xxyy_xxxxyyzz_1, tg_xxyy_xxxxyzz_1, tg_xxyy_xxxxyzzz_0, \
                                         tg_xxyy_xxxxyzzz_1, tg_xxyy_xxxxzzz_1, tg_xxyy_xxxxzzzz_0, tg_xxyy_xxxxzzzz_1, \
                                         tg_xxyy_xxxyyyy_1, tg_xxyy_xxxyyyyy_0, tg_xxyy_xxxyyyyy_1, tg_xxyy_xxxyyyyz_0, \
                                         tg_xxyy_xxxyyyyz_1, tg_xxyy_xxxyyyz_1, tg_xxyy_xxxyyyzz_0, tg_xxyy_xxxyyyzz_1, \
                                         tg_xxyy_xxxyyzz_1, tg_xxyy_xxxyyzzz_0, tg_xxyy_xxxyyzzz_1, tg_xxyy_xxxyzzz_1, \
                                         tg_xxyy_xxxyzzzz_0, tg_xxyy_xxxyzzzz_1, tg_xxyy_xxxzzzz_1, tg_xxyy_xxxzzzzz_0, \
                                         tg_xxyy_xxxzzzzz_1, tg_xxyy_xxyyyyy_1, tg_xxyy_xxyyyyyy_0, tg_xxyy_xxyyyyyy_1, \
                                         tg_xxyy_xxyyyyyz_0, tg_xxyy_xxyyyyyz_1, tg_xxyy_xxyyyyz_1, tg_xxyy_xxyyyyzz_0, \
                                         tg_xxyy_xxyyyyzz_1, tg_xxyy_xxyyyzz_1, tg_xxyy_xxyyyzzz_0, tg_xxyy_xxyyyzzz_1, \
                                         tg_xxyy_xxyyzzz_1, tg_xxyy_xxyyzzzz_0, tg_xxyy_xxyyzzzz_1, tg_xxyy_xxyzzzz_1, \
                                         tg_xxyy_xxyzzzzz_0, tg_xxyy_xxyzzzzz_1, tg_xxyy_xxzzzzz_1, tg_xxyy_xxzzzzzz_0, \
                                         tg_xxyy_xxzzzzzz_1, tg_xxyy_xyyyyyy_1, tg_xxyy_xyyyyyyy_0, tg_xxyy_xyyyyyyy_1, \
                                         tg_xxyy_xyyyyyyz_0, tg_xxyy_xyyyyyyz_1, tg_xxyy_xyyyyyz_1, tg_xxyy_xyyyyyzz_0, \
                                         tg_xxyy_xyyyyyzz_1, tg_xxyy_xyyyyzz_1, tg_xxyy_xyyyyzzz_0, tg_xxyy_xyyyyzzz_1, \
                                         tg_xxyy_xyyyzzz_1, tg_xxyy_xyyyzzzz_0, tg_xxyy_xyyyzzzz_1, tg_xxyy_xyyzzzz_1, \
                                         tg_xxyy_xyyzzzzz_0, tg_xxyy_xyyzzzzz_1, tg_xxyy_xyzzzzz_1, tg_xxyy_xyzzzzzz_0, \
                                         tg_xxyy_xyzzzzzz_1, tg_xxyy_xzzzzzz_1, tg_xxyy_xzzzzzzz_0, tg_xxyy_xzzzzzzz_1, \
                                         tg_xxyy_yyyyyyy_1, tg_xxyy_yyyyyyyy_0, tg_xxyy_yyyyyyyy_1, tg_xxyy_yyyyyyyz_0, \
                                         tg_xxyy_yyyyyyyz_1, tg_xxyy_yyyyyyz_1, tg_xxyy_yyyyyyzz_0, tg_xxyy_yyyyyyzz_1, \
                                         tg_xxyy_yyyyyzz_1, tg_xxyy_yyyyyzzz_0, tg_xxyy_yyyyyzzz_1, tg_xxyy_yyyyzzz_1, \
                                         tg_xxyy_yyyyzzzz_0, tg_xxyy_yyyyzzzz_1, tg_xxyy_yyyzzzz_1, tg_xxyy_yyyzzzzz_0, \
                                         tg_xxyy_yyyzzzzz_1, tg_xxyy_yyzzzzz_1, tg_xxyy_yyzzzzzz_0, tg_xxyy_yyzzzzzz_1, \
                                         tg_xxyy_yzzzzzz_1, tg_xxyy_yzzzzzzz_0, tg_xxyy_yzzzzzzz_1, tg_xxyy_zzzzzzz_1, \
                                         tg_xxyy_zzzzzzzz_0, tg_xxyy_zzzzzzzz_1, tg_xxyz_xxxxxxx_1, tg_xxyz_xxxxxxxx_0, \
                                         tg_xxyz_xxxxxxxx_1, tg_xxyz_xxxxxxxy_0, tg_xxyz_xxxxxxxy_1, tg_xxyz_xxxxxxxz_0, \
                                         tg_xxyz_xxxxxxxz_1, tg_xxyz_xxxxxxy_1, tg_xxyz_xxxxxxyy_0, tg_xxyz_xxxxxxyy_1, \
                                         tg_xxyz_xxxxxxyz_0, tg_xxyz_xxxxxxyz_1, tg_xxyz_xxxxxxz_1, tg_xxyz_xxxxxxzz_0, \
                                         tg_xxyz_xxxxxxzz_1, tg_xxyz_xxxxxyy_1, tg_xxyz_xxxxxyyy_0, tg_xxyz_xxxxxyyy_1, \
                                         tg_xxyz_xxxxxyyz_0, tg_xxyz_xxxxxyyz_1, tg_xxyz_xxxxxyz_1, tg_xxyz_xxxxxyzz_0, \
                                         tg_xxyz_xxxxxyzz_1, tg_xxyz_xxxxxzz_1, tg_xxyz_xxxxxzzz_0, tg_xxyz_xxxxxzzz_1, \
                                         tg_xxyz_xxxxyyy_1, tg_xxyz_xxxxyyz_1, tg_xxyz_xxxxyzz_1, tg_xxyz_xxxxzzz_1, \
                                         tg_xxz_xxxxxxzz_0, tg_xxz_xxxxxxzz_1, tg_xxz_xxxxxyyy_0, tg_xxz_xxxxxyyy_1, \
                                         tg_xxz_xxxxxyyz_0, tg_xxz_xxxxxyyz_1, tg_xxz_xxxxxyzz_0, tg_xxz_xxxxxyzz_1, \
                                         tg_xxz_xxxxxzzz_0, tg_xxz_xxxxxzzz_1, tg_xxz_xxxxyyyy_0, tg_xxz_xxxxyyyy_1, \
                                         tg_xxz_xxxxyyyz_0, tg_xxz_xxxxyyyz_1, tg_xxz_xxxxyyzz_0, tg_xxz_xxxxyyzz_1, \
                                         tg_xxz_xxxxyzzz_0, tg_xxz_xxxxyzzz_1, tg_xxz_xxxxzzzz_0, tg_xxz_xxxxzzzz_1, \
                                         tg_xxz_xxxyyyyy_0, tg_xxz_xxxyyyyy_1, tg_xxz_xxxyyyyz_0, tg_xxz_xxxyyyyz_1, \
                                         tg_xxz_xxxyyyzz_0, tg_xxz_xxxyyyzz_1, tg_xxz_xxxyyzzz_0, tg_xxz_xxxyyzzz_1, \
                                         tg_xxz_xxxyzzzz_0, tg_xxz_xxxyzzzz_1, tg_xxz_xxxzzzzz_0, tg_xxz_xxxzzzzz_1, \
                                         tg_xxz_xxyyyyyy_0, tg_xxz_xxyyyyyy_1, tg_xxz_xxyyyyyz_0, tg_xxz_xxyyyyyz_1, \
                                         tg_xxz_xxyyyyzz_0, tg_xxz_xxyyyyzz_1, tg_xxz_xxyyyzzz_0, tg_xxz_xxyyyzzz_1, \
                                         tg_xxz_xxyyzzzz_0, tg_xxz_xxyyzzzz_1, tg_xxz_xxyzzzzz_0, tg_xxz_xxyzzzzz_1, \
                                         tg_xxz_xxzzzzzz_0, tg_xxz_xxzzzzzz_1, tg_xxz_xyyyyyyy_0, tg_xxz_xyyyyyyy_1, \
                                         tg_xxz_xyyyyyyz_0, tg_xxz_xyyyyyyz_1, tg_xxz_xyyyyyzz_0, tg_xxz_xyyyyyzz_1, \
                                         tg_xxz_xyyyyzzz_0, tg_xxz_xyyyyzzz_1, tg_xxz_xyyyzzzz_0, tg_xxz_xyyyzzzz_1, \
                                         tg_xxz_xyyzzzzz_0, tg_xxz_xyyzzzzz_1, tg_xxz_xyzzzzzz_0, tg_xxz_xyzzzzzz_1, \
                                         tg_xxz_xzzzzzzz_0, tg_xxz_xzzzzzzz_1, tg_xxz_yyyyyyyy_0, tg_xxz_yyyyyyyy_1, \
                                         tg_xxz_yyyyyyyz_0, tg_xxz_yyyyyyyz_1, tg_xxz_yyyyyyzz_0, tg_xxz_yyyyyyzz_1, \
                                         tg_xxz_yyyyyzzz_0, tg_xxz_yyyyyzzz_1, tg_xxz_yyyyzzzz_0, tg_xxz_yyyyzzzz_1, \
                                         tg_xxz_yyyzzzzz_0, tg_xxz_yyyzzzzz_1, tg_xxz_yyzzzzzz_0, tg_xxz_yyzzzzzz_1, \
                                         tg_xxz_yzzzzzzz_0, tg_xxz_yzzzzzzz_1, tg_xxz_zzzzzzzz_0, tg_xxz_zzzzzzzz_1, \
                                         tg_xyy_xxxxxxxx_0, tg_xyy_xxxxxxxx_1, tg_xyy_xxxxxxxy_0, tg_xyy_xxxxxxxy_1, \
                                         tg_xyy_xxxxxxxz_0, tg_xyy_xxxxxxxz_1, tg_xyy_xxxxxxyy_0, tg_xyy_xxxxxxyy_1, \
                                         tg_xyy_xxxxxxyz_0, tg_xyy_xxxxxxyz_1, tg_xyy_xxxxxxzz_0, tg_xyy_xxxxxxzz_1, \
                                         tg_xyy_xxxxxyyy_0, tg_xyy_xxxxxyyy_1, tg_xyy_xxxxxyyz_0, tg_xyy_xxxxxyyz_1, \
                                         tg_xyy_xxxxxyzz_0, tg_xyy_xxxxxyzz_1, tg_xyy_xxxxxzzz_0, tg_xyy_xxxxxzzz_1, \
                                         tg_xyy_xxxxyyyy_0, tg_xyy_xxxxyyyy_1, tg_xyy_xxxxyyyz_0, tg_xyy_xxxxyyyz_1, \
                                         tg_xyy_xxxxyyzz_0, tg_xyy_xxxxyyzz_1, tg_xyy_xxxxyzzz_0, tg_xyy_xxxxyzzz_1, \
                                         tg_xyy_xxxxzzzz_0, tg_xyy_xxxxzzzz_1, tg_xyy_xxxyyyyy_0, tg_xyy_xxxyyyyy_1, \
                                         tg_xyy_xxxyyyyz_0, tg_xyy_xxxyyyyz_1, tg_xyy_xxxyyyzz_0, tg_xyy_xxxyyyzz_1, \
                                         tg_xyy_xxxyyzzz_0, tg_xyy_xxxyyzzz_1, tg_xyy_xxxyzzzz_0, tg_xyy_xxxyzzzz_1, \
                                         tg_xyy_xxxzzzzz_0, tg_xyy_xxxzzzzz_1, tg_xyy_xxyyyyyy_0, tg_xyy_xxyyyyyy_1, \
                                         tg_xyy_xxyyyyyz_0, tg_xyy_xxyyyyyz_1, tg_xyy_xxyyyyzz_0, tg_xyy_xxyyyyzz_1, \
                                         tg_xyy_xxyyyzzz_0, tg_xyy_xxyyyzzz_1, tg_xyy_xxyyzzzz_0, tg_xyy_xxyyzzzz_1, \
                                         tg_xyy_xxyzzzzz_0, tg_xyy_xxyzzzzz_1, tg_xyy_xxzzzzzz_0, tg_xyy_xxzzzzzz_1, \
                                         tg_xyy_xyyyyyyy_0, tg_xyy_xyyyyyyy_1, tg_xyy_xyyyyyyz_0, tg_xyy_xyyyyyyz_1, \
                                         tg_xyy_xyyyyyzz_0, tg_xyy_xyyyyyzz_1, tg_xyy_xyyyyzzz_0, tg_xyy_xyyyyzzz_1, \
                                         tg_xyy_xyyyzzzz_0, tg_xyy_xyyyzzzz_1, tg_xyy_xyyzzzzz_0, tg_xyy_xyyzzzzz_1, \
                                         tg_xyy_xyzzzzzz_0, tg_xyy_xyzzzzzz_1, tg_xyy_xzzzzzzz_0, tg_xyy_xzzzzzzz_1, \
                                         tg_xyy_yyyyyyyy_0, tg_xyy_yyyyyyyy_1, tg_xyy_yyyyyyyz_0, tg_xyy_yyyyyyyz_1, \
                                         tg_xyy_yyyyyyzz_0, tg_xyy_yyyyyyzz_1, tg_xyy_yyyyyzzz_0, tg_xyy_yyyyyzzz_1, \
                                         tg_xyy_yyyyzzzz_0, tg_xyy_yyyyzzzz_1, tg_xyy_yyyzzzzz_0, tg_xyy_yyyzzzzz_1, \
                                         tg_xyy_yyzzzzzz_0, tg_xyy_yyzzzzzz_1, tg_xyy_yzzzzzzz_0, tg_xyy_yzzzzzzz_1, \
                                         tg_xyy_zzzzzzzz_0, tg_xyy_zzzzzzzz_1, tg_xyz_xxxxxxxx_0, tg_xyz_xxxxxxxx_1, \
                                         tg_xyz_xxxxxxxy_0, tg_xyz_xxxxxxxy_1, tg_xyz_xxxxxxxz_0, tg_xyz_xxxxxxxz_1, \
                                         tg_xyz_xxxxxxyy_0, tg_xyz_xxxxxxyy_1, tg_xyz_xxxxxxyz_0, tg_xyz_xxxxxxyz_1, \
                                         tg_xyz_xxxxxxzz_0, tg_xyz_xxxxxxzz_1, tg_xyz_xxxxxyyy_0, tg_xyz_xxxxxyyy_1, \
                                         tg_xyz_xxxxxyyz_0, tg_xyz_xxxxxyyz_1, tg_xyz_xxxxxyzz_0, tg_xyz_xxxxxyzz_1, \
                                         tg_xyz_xxxxxzzz_0, tg_xyz_xxxxxzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxxz_xxxxxxzz_0[j] = pb_x * tg_xxxz_xxxxxxzz_0[j] + fr * tg_xxxz_xxxxxxzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxxxxzz_0[j] - tg_xxz_xxxxxxzz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxxz_xxxxxzz_1[j];

                    tg_xxxxz_xxxxxyyy_0[j] = pb_x * tg_xxxz_xxxxxyyy_0[j] + fr * tg_xxxz_xxxxxyyy_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxxxyyy_0[j] - tg_xxz_xxxxxyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxz_xxxxyyy_1[j];

                    tg_xxxxz_xxxxxyyz_0[j] = pb_x * tg_xxxz_xxxxxyyz_0[j] + fr * tg_xxxz_xxxxxyyz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxxxyyz_0[j] - tg_xxz_xxxxxyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxz_xxxxyyz_1[j];

                    tg_xxxxz_xxxxxyzz_0[j] = pb_x * tg_xxxz_xxxxxyzz_0[j] + fr * tg_xxxz_xxxxxyzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxxxyzz_0[j] - tg_xxz_xxxxxyzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxz_xxxxyzz_1[j];

                    tg_xxxxz_xxxxxzzz_0[j] = pb_x * tg_xxxz_xxxxxzzz_0[j] + fr * tg_xxxz_xxxxxzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxxxzzz_0[j] - tg_xxz_xxxxxzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxz_xxxxzzz_1[j];

                    tg_xxxxz_xxxxyyyy_0[j] = pb_x * tg_xxxz_xxxxyyyy_0[j] + fr * tg_xxxz_xxxxyyyy_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxxyyyy_0[j] - tg_xxz_xxxxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxz_xxxyyyy_1[j];

                    tg_xxxxz_xxxxyyyz_0[j] = pb_x * tg_xxxz_xxxxyyyz_0[j] + fr * tg_xxxz_xxxxyyyz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxxyyyz_0[j] - tg_xxz_xxxxyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxz_xxxyyyz_1[j];

                    tg_xxxxz_xxxxyyzz_0[j] = pb_x * tg_xxxz_xxxxyyzz_0[j] + fr * tg_xxxz_xxxxyyzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxxyyzz_0[j] - tg_xxz_xxxxyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxz_xxxyyzz_1[j];

                    tg_xxxxz_xxxxyzzz_0[j] = pb_x * tg_xxxz_xxxxyzzz_0[j] + fr * tg_xxxz_xxxxyzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxxyzzz_0[j] - tg_xxz_xxxxyzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxz_xxxyzzz_1[j];

                    tg_xxxxz_xxxxzzzz_0[j] = pb_x * tg_xxxz_xxxxzzzz_0[j] + fr * tg_xxxz_xxxxzzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxxzzzz_0[j] - tg_xxz_xxxxzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxz_xxxzzzz_1[j];

                    tg_xxxxz_xxxyyyyy_0[j] = pb_x * tg_xxxz_xxxyyyyy_0[j] + fr * tg_xxxz_xxxyyyyy_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxyyyyy_0[j] - tg_xxz_xxxyyyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxz_xxyyyyy_1[j];

                    tg_xxxxz_xxxyyyyz_0[j] = pb_x * tg_xxxz_xxxyyyyz_0[j] + fr * tg_xxxz_xxxyyyyz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxyyyyz_0[j] - tg_xxz_xxxyyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxz_xxyyyyz_1[j];

                    tg_xxxxz_xxxyyyzz_0[j] = pb_x * tg_xxxz_xxxyyyzz_0[j] + fr * tg_xxxz_xxxyyyzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxyyyzz_0[j] - tg_xxz_xxxyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxz_xxyyyzz_1[j];

                    tg_xxxxz_xxxyyzzz_0[j] = pb_x * tg_xxxz_xxxyyzzz_0[j] + fr * tg_xxxz_xxxyyzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxyyzzz_0[j] - tg_xxz_xxxyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxz_xxyyzzz_1[j];

                    tg_xxxxz_xxxyzzzz_0[j] = pb_x * tg_xxxz_xxxyzzzz_0[j] + fr * tg_xxxz_xxxyzzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxyzzzz_0[j] - tg_xxz_xxxyzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxz_xxyzzzz_1[j];

                    tg_xxxxz_xxxzzzzz_0[j] = pb_x * tg_xxxz_xxxzzzzz_0[j] + fr * tg_xxxz_xxxzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxzzzzz_0[j] - tg_xxz_xxxzzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxz_xxzzzzz_1[j];

                    tg_xxxxz_xxyyyyyy_0[j] = pb_x * tg_xxxz_xxyyyyyy_0[j] + fr * tg_xxxz_xxyyyyyy_1[j] + 1.5 * fl1_fx * (tg_xxz_xxyyyyyy_0[j] - tg_xxz_xxyyyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxz_xyyyyyy_1[j];

                    tg_xxxxz_xxyyyyyz_0[j] = pb_x * tg_xxxz_xxyyyyyz_0[j] + fr * tg_xxxz_xxyyyyyz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxyyyyyz_0[j] - tg_xxz_xxyyyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxz_xyyyyyz_1[j];

                    tg_xxxxz_xxyyyyzz_0[j] = pb_x * tg_xxxz_xxyyyyzz_0[j] + fr * tg_xxxz_xxyyyyzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxyyyyzz_0[j] - tg_xxz_xxyyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxz_xyyyyzz_1[j];

                    tg_xxxxz_xxyyyzzz_0[j] = pb_x * tg_xxxz_xxyyyzzz_0[j] + fr * tg_xxxz_xxyyyzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxyyyzzz_0[j] - tg_xxz_xxyyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxz_xyyyzzz_1[j];

                    tg_xxxxz_xxyyzzzz_0[j] = pb_x * tg_xxxz_xxyyzzzz_0[j] + fr * tg_xxxz_xxyyzzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxyyzzzz_0[j] - tg_xxz_xxyyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxz_xyyzzzz_1[j];

                    tg_xxxxz_xxyzzzzz_0[j] = pb_x * tg_xxxz_xxyzzzzz_0[j] + fr * tg_xxxz_xxyzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxyzzzzz_0[j] - tg_xxz_xxyzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxz_xyzzzzz_1[j];

                    tg_xxxxz_xxzzzzzz_0[j] = pb_x * tg_xxxz_xxzzzzzz_0[j] + fr * tg_xxxz_xxzzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxzzzzzz_0[j] - tg_xxz_xxzzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxz_xzzzzzz_1[j];

                    tg_xxxxz_xyyyyyyy_0[j] = pb_x * tg_xxxz_xyyyyyyy_0[j] + fr * tg_xxxz_xyyyyyyy_1[j] + 1.5 * fl1_fx * (tg_xxz_xyyyyyyy_0[j] - tg_xxz_xyyyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxz_yyyyyyy_1[j];

                    tg_xxxxz_xyyyyyyz_0[j] = pb_x * tg_xxxz_xyyyyyyz_0[j] + fr * tg_xxxz_xyyyyyyz_1[j] + 1.5 * fl1_fx * (tg_xxz_xyyyyyyz_0[j] - tg_xxz_xyyyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxz_yyyyyyz_1[j];

                    tg_xxxxz_xyyyyyzz_0[j] = pb_x * tg_xxxz_xyyyyyzz_0[j] + fr * tg_xxxz_xyyyyyzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xyyyyyzz_0[j] - tg_xxz_xyyyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxz_yyyyyzz_1[j];

                    tg_xxxxz_xyyyyzzz_0[j] = pb_x * tg_xxxz_xyyyyzzz_0[j] + fr * tg_xxxz_xyyyyzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xyyyyzzz_0[j] - tg_xxz_xyyyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxz_yyyyzzz_1[j];

                    tg_xxxxz_xyyyzzzz_0[j] = pb_x * tg_xxxz_xyyyzzzz_0[j] + fr * tg_xxxz_xyyyzzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xyyyzzzz_0[j] - tg_xxz_xyyyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxz_yyyzzzz_1[j];

                    tg_xxxxz_xyyzzzzz_0[j] = pb_x * tg_xxxz_xyyzzzzz_0[j] + fr * tg_xxxz_xyyzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xyyzzzzz_0[j] - tg_xxz_xyyzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxz_yyzzzzz_1[j];

                    tg_xxxxz_xyzzzzzz_0[j] = pb_x * tg_xxxz_xyzzzzzz_0[j] + fr * tg_xxxz_xyzzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xyzzzzzz_0[j] - tg_xxz_xyzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxz_yzzzzzz_1[j];

                    tg_xxxxz_xzzzzzzz_0[j] = pb_x * tg_xxxz_xzzzzzzz_0[j] + fr * tg_xxxz_xzzzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xzzzzzzz_0[j] - tg_xxz_xzzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxz_zzzzzzz_1[j];

                    tg_xxxxz_yyyyyyyy_0[j] = pb_x * tg_xxxz_yyyyyyyy_0[j] + fr * tg_xxxz_yyyyyyyy_1[j] + 1.5 * fl1_fx * (tg_xxz_yyyyyyyy_0[j] - tg_xxz_yyyyyyyy_1[j] * fl1_fza);

                    tg_xxxxz_yyyyyyyz_0[j] = pb_x * tg_xxxz_yyyyyyyz_0[j] + fr * tg_xxxz_yyyyyyyz_1[j] + 1.5 * fl1_fx * (tg_xxz_yyyyyyyz_0[j] - tg_xxz_yyyyyyyz_1[j] * fl1_fza);

                    tg_xxxxz_yyyyyyzz_0[j] = pb_x * tg_xxxz_yyyyyyzz_0[j] + fr * tg_xxxz_yyyyyyzz_1[j] + 1.5 * fl1_fx * (tg_xxz_yyyyyyzz_0[j] - tg_xxz_yyyyyyzz_1[j] * fl1_fza);

                    tg_xxxxz_yyyyyzzz_0[j] = pb_x * tg_xxxz_yyyyyzzz_0[j] + fr * tg_xxxz_yyyyyzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_yyyyyzzz_0[j] - tg_xxz_yyyyyzzz_1[j] * fl1_fza);

                    tg_xxxxz_yyyyzzzz_0[j] = pb_x * tg_xxxz_yyyyzzzz_0[j] + fr * tg_xxxz_yyyyzzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_yyyyzzzz_0[j] - tg_xxz_yyyyzzzz_1[j] * fl1_fza);

                    tg_xxxxz_yyyzzzzz_0[j] = pb_x * tg_xxxz_yyyzzzzz_0[j] + fr * tg_xxxz_yyyzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_yyyzzzzz_0[j] - tg_xxz_yyyzzzzz_1[j] * fl1_fza);

                    tg_xxxxz_yyzzzzzz_0[j] = pb_x * tg_xxxz_yyzzzzzz_0[j] + fr * tg_xxxz_yyzzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_yyzzzzzz_0[j] - tg_xxz_yyzzzzzz_1[j] * fl1_fza);

                    tg_xxxxz_yzzzzzzz_0[j] = pb_x * tg_xxxz_yzzzzzzz_0[j] + fr * tg_xxxz_yzzzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_yzzzzzzz_0[j] - tg_xxz_yzzzzzzz_1[j] * fl1_fza);

                    tg_xxxxz_zzzzzzzz_0[j] = pb_x * tg_xxxz_zzzzzzzz_0[j] + fr * tg_xxxz_zzzzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_zzzzzzzz_0[j] - tg_xxz_zzzzzzzz_1[j] * fl1_fza);

                    tg_xxxyy_xxxxxxxx_0[j] = pb_x * tg_xxyy_xxxxxxxx_0[j] + fr * tg_xxyy_xxxxxxxx_1[j] + fl1_fx * (tg_xyy_xxxxxxxx_0[j] - tg_xyy_xxxxxxxx_1[j] * fl1_fza) + 4.0 * fl1_fxn * tg_xxyy_xxxxxxx_1[j];

                    tg_xxxyy_xxxxxxxy_0[j] = pb_x * tg_xxyy_xxxxxxxy_0[j] + fr * tg_xxyy_xxxxxxxy_1[j] + fl1_fx * (tg_xyy_xxxxxxxy_0[j] - tg_xyy_xxxxxxxy_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_xxyy_xxxxxxy_1[j];

                    tg_xxxyy_xxxxxxxz_0[j] = pb_x * tg_xxyy_xxxxxxxz_0[j] + fr * tg_xxyy_xxxxxxxz_1[j] + fl1_fx * (tg_xyy_xxxxxxxz_0[j] - tg_xyy_xxxxxxxz_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_xxyy_xxxxxxz_1[j];

                    tg_xxxyy_xxxxxxyy_0[j] = pb_x * tg_xxyy_xxxxxxyy_0[j] + fr * tg_xxyy_xxxxxxyy_1[j] + fl1_fx * (tg_xyy_xxxxxxyy_0[j] - tg_xyy_xxxxxxyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxyy_xxxxxyy_1[j];

                    tg_xxxyy_xxxxxxyz_0[j] = pb_x * tg_xxyy_xxxxxxyz_0[j] + fr * tg_xxyy_xxxxxxyz_1[j] + fl1_fx * (tg_xyy_xxxxxxyz_0[j] - tg_xyy_xxxxxxyz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxyy_xxxxxyz_1[j];

                    tg_xxxyy_xxxxxxzz_0[j] = pb_x * tg_xxyy_xxxxxxzz_0[j] + fr * tg_xxyy_xxxxxxzz_1[j] + fl1_fx * (tg_xyy_xxxxxxzz_0[j] - tg_xyy_xxxxxxzz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxyy_xxxxxzz_1[j];

                    tg_xxxyy_xxxxxyyy_0[j] = pb_x * tg_xxyy_xxxxxyyy_0[j] + fr * tg_xxyy_xxxxxyyy_1[j] + fl1_fx * (tg_xyy_xxxxxyyy_0[j] - tg_xyy_xxxxxyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyy_xxxxyyy_1[j];

                    tg_xxxyy_xxxxxyyz_0[j] = pb_x * tg_xxyy_xxxxxyyz_0[j] + fr * tg_xxyy_xxxxxyyz_1[j] + fl1_fx * (tg_xyy_xxxxxyyz_0[j] - tg_xyy_xxxxxyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyy_xxxxyyz_1[j];

                    tg_xxxyy_xxxxxyzz_0[j] = pb_x * tg_xxyy_xxxxxyzz_0[j] + fr * tg_xxyy_xxxxxyzz_1[j] + fl1_fx * (tg_xyy_xxxxxyzz_0[j] - tg_xyy_xxxxxyzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyy_xxxxyzz_1[j];

                    tg_xxxyy_xxxxxzzz_0[j] = pb_x * tg_xxyy_xxxxxzzz_0[j] + fr * tg_xxyy_xxxxxzzz_1[j] + fl1_fx * (tg_xyy_xxxxxzzz_0[j] - tg_xyy_xxxxxzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyy_xxxxzzz_1[j];

                    tg_xxxyy_xxxxyyyy_0[j] = pb_x * tg_xxyy_xxxxyyyy_0[j] + fr * tg_xxyy_xxxxyyyy_1[j] + fl1_fx * (tg_xyy_xxxxyyyy_0[j] - tg_xyy_xxxxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyy_xxxyyyy_1[j];

                    tg_xxxyy_xxxxyyyz_0[j] = pb_x * tg_xxyy_xxxxyyyz_0[j] + fr * tg_xxyy_xxxxyyyz_1[j] + fl1_fx * (tg_xyy_xxxxyyyz_0[j] - tg_xyy_xxxxyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyy_xxxyyyz_1[j];

                    tg_xxxyy_xxxxyyzz_0[j] = pb_x * tg_xxyy_xxxxyyzz_0[j] + fr * tg_xxyy_xxxxyyzz_1[j] + fl1_fx * (tg_xyy_xxxxyyzz_0[j] - tg_xyy_xxxxyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyy_xxxyyzz_1[j];

                    tg_xxxyy_xxxxyzzz_0[j] = pb_x * tg_xxyy_xxxxyzzz_0[j] + fr * tg_xxyy_xxxxyzzz_1[j] + fl1_fx * (tg_xyy_xxxxyzzz_0[j] - tg_xyy_xxxxyzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyy_xxxyzzz_1[j];

                    tg_xxxyy_xxxxzzzz_0[j] = pb_x * tg_xxyy_xxxxzzzz_0[j] + fr * tg_xxyy_xxxxzzzz_1[j] + fl1_fx * (tg_xyy_xxxxzzzz_0[j] - tg_xyy_xxxxzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyy_xxxzzzz_1[j];

                    tg_xxxyy_xxxyyyyy_0[j] = pb_x * tg_xxyy_xxxyyyyy_0[j] + fr * tg_xxyy_xxxyyyyy_1[j] + fl1_fx * (tg_xyy_xxxyyyyy_0[j] - tg_xyy_xxxyyyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyy_xxyyyyy_1[j];

                    tg_xxxyy_xxxyyyyz_0[j] = pb_x * tg_xxyy_xxxyyyyz_0[j] + fr * tg_xxyy_xxxyyyyz_1[j] + fl1_fx * (tg_xyy_xxxyyyyz_0[j] - tg_xyy_xxxyyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyy_xxyyyyz_1[j];

                    tg_xxxyy_xxxyyyzz_0[j] = pb_x * tg_xxyy_xxxyyyzz_0[j] + fr * tg_xxyy_xxxyyyzz_1[j] + fl1_fx * (tg_xyy_xxxyyyzz_0[j] - tg_xyy_xxxyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyy_xxyyyzz_1[j];

                    tg_xxxyy_xxxyyzzz_0[j] = pb_x * tg_xxyy_xxxyyzzz_0[j] + fr * tg_xxyy_xxxyyzzz_1[j] + fl1_fx * (tg_xyy_xxxyyzzz_0[j] - tg_xyy_xxxyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyy_xxyyzzz_1[j];

                    tg_xxxyy_xxxyzzzz_0[j] = pb_x * tg_xxyy_xxxyzzzz_0[j] + fr * tg_xxyy_xxxyzzzz_1[j] + fl1_fx * (tg_xyy_xxxyzzzz_0[j] - tg_xyy_xxxyzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyy_xxyzzzz_1[j];

                    tg_xxxyy_xxxzzzzz_0[j] = pb_x * tg_xxyy_xxxzzzzz_0[j] + fr * tg_xxyy_xxxzzzzz_1[j] + fl1_fx * (tg_xyy_xxxzzzzz_0[j] - tg_xyy_xxxzzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyy_xxzzzzz_1[j];

                    tg_xxxyy_xxyyyyyy_0[j] = pb_x * tg_xxyy_xxyyyyyy_0[j] + fr * tg_xxyy_xxyyyyyy_1[j] + fl1_fx * (tg_xyy_xxyyyyyy_0[j] - tg_xyy_xxyyyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxyy_xyyyyyy_1[j];

                    tg_xxxyy_xxyyyyyz_0[j] = pb_x * tg_xxyy_xxyyyyyz_0[j] + fr * tg_xxyy_xxyyyyyz_1[j] + fl1_fx * (tg_xyy_xxyyyyyz_0[j] - tg_xyy_xxyyyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxyy_xyyyyyz_1[j];

                    tg_xxxyy_xxyyyyzz_0[j] = pb_x * tg_xxyy_xxyyyyzz_0[j] + fr * tg_xxyy_xxyyyyzz_1[j] + fl1_fx * (tg_xyy_xxyyyyzz_0[j] - tg_xyy_xxyyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyy_xyyyyzz_1[j];

                    tg_xxxyy_xxyyyzzz_0[j] = pb_x * tg_xxyy_xxyyyzzz_0[j] + fr * tg_xxyy_xxyyyzzz_1[j] + fl1_fx * (tg_xyy_xxyyyzzz_0[j] - tg_xyy_xxyyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyy_xyyyzzz_1[j];

                    tg_xxxyy_xxyyzzzz_0[j] = pb_x * tg_xxyy_xxyyzzzz_0[j] + fr * tg_xxyy_xxyyzzzz_1[j] + fl1_fx * (tg_xyy_xxyyzzzz_0[j] - tg_xyy_xxyyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyy_xyyzzzz_1[j];

                    tg_xxxyy_xxyzzzzz_0[j] = pb_x * tg_xxyy_xxyzzzzz_0[j] + fr * tg_xxyy_xxyzzzzz_1[j] + fl1_fx * (tg_xyy_xxyzzzzz_0[j] - tg_xyy_xxyzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyy_xyzzzzz_1[j];

                    tg_xxxyy_xxzzzzzz_0[j] = pb_x * tg_xxyy_xxzzzzzz_0[j] + fr * tg_xxyy_xxzzzzzz_1[j] + fl1_fx * (tg_xyy_xxzzzzzz_0[j] - tg_xyy_xxzzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyy_xzzzzzz_1[j];

                    tg_xxxyy_xyyyyyyy_0[j] = pb_x * tg_xxyy_xyyyyyyy_0[j] + fr * tg_xxyy_xyyyyyyy_1[j] + fl1_fx * (tg_xyy_xyyyyyyy_0[j] - tg_xyy_xyyyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyy_yyyyyyy_1[j];

                    tg_xxxyy_xyyyyyyz_0[j] = pb_x * tg_xxyy_xyyyyyyz_0[j] + fr * tg_xxyy_xyyyyyyz_1[j] + fl1_fx * (tg_xyy_xyyyyyyz_0[j] - tg_xyy_xyyyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyy_yyyyyyz_1[j];

                    tg_xxxyy_xyyyyyzz_0[j] = pb_x * tg_xxyy_xyyyyyzz_0[j] + fr * tg_xxyy_xyyyyyzz_1[j] + fl1_fx * (tg_xyy_xyyyyyzz_0[j] - tg_xyy_xyyyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyy_yyyyyzz_1[j];

                    tg_xxxyy_xyyyyzzz_0[j] = pb_x * tg_xxyy_xyyyyzzz_0[j] + fr * tg_xxyy_xyyyyzzz_1[j] + fl1_fx * (tg_xyy_xyyyyzzz_0[j] - tg_xyy_xyyyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyy_yyyyzzz_1[j];

                    tg_xxxyy_xyyyzzzz_0[j] = pb_x * tg_xxyy_xyyyzzzz_0[j] + fr * tg_xxyy_xyyyzzzz_1[j] + fl1_fx * (tg_xyy_xyyyzzzz_0[j] - tg_xyy_xyyyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyy_yyyzzzz_1[j];

                    tg_xxxyy_xyyzzzzz_0[j] = pb_x * tg_xxyy_xyyzzzzz_0[j] + fr * tg_xxyy_xyyzzzzz_1[j] + fl1_fx * (tg_xyy_xyyzzzzz_0[j] - tg_xyy_xyyzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyy_yyzzzzz_1[j];

                    tg_xxxyy_xyzzzzzz_0[j] = pb_x * tg_xxyy_xyzzzzzz_0[j] + fr * tg_xxyy_xyzzzzzz_1[j] + fl1_fx * (tg_xyy_xyzzzzzz_0[j] - tg_xyy_xyzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyy_yzzzzzz_1[j];

                    tg_xxxyy_xzzzzzzz_0[j] = pb_x * tg_xxyy_xzzzzzzz_0[j] + fr * tg_xxyy_xzzzzzzz_1[j] + fl1_fx * (tg_xyy_xzzzzzzz_0[j] - tg_xyy_xzzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyy_zzzzzzz_1[j];

                    tg_xxxyy_yyyyyyyy_0[j] = pb_x * tg_xxyy_yyyyyyyy_0[j] + fr * tg_xxyy_yyyyyyyy_1[j] + fl1_fx * (tg_xyy_yyyyyyyy_0[j] - tg_xyy_yyyyyyyy_1[j] * fl1_fza);

                    tg_xxxyy_yyyyyyyz_0[j] = pb_x * tg_xxyy_yyyyyyyz_0[j] + fr * tg_xxyy_yyyyyyyz_1[j] + fl1_fx * (tg_xyy_yyyyyyyz_0[j] - tg_xyy_yyyyyyyz_1[j] * fl1_fza);

                    tg_xxxyy_yyyyyyzz_0[j] = pb_x * tg_xxyy_yyyyyyzz_0[j] + fr * tg_xxyy_yyyyyyzz_1[j] + fl1_fx * (tg_xyy_yyyyyyzz_0[j] - tg_xyy_yyyyyyzz_1[j] * fl1_fza);

                    tg_xxxyy_yyyyyzzz_0[j] = pb_x * tg_xxyy_yyyyyzzz_0[j] + fr * tg_xxyy_yyyyyzzz_1[j] + fl1_fx * (tg_xyy_yyyyyzzz_0[j] - tg_xyy_yyyyyzzz_1[j] * fl1_fza);

                    tg_xxxyy_yyyyzzzz_0[j] = pb_x * tg_xxyy_yyyyzzzz_0[j] + fr * tg_xxyy_yyyyzzzz_1[j] + fl1_fx * (tg_xyy_yyyyzzzz_0[j] - tg_xyy_yyyyzzzz_1[j] * fl1_fza);

                    tg_xxxyy_yyyzzzzz_0[j] = pb_x * tg_xxyy_yyyzzzzz_0[j] + fr * tg_xxyy_yyyzzzzz_1[j] + fl1_fx * (tg_xyy_yyyzzzzz_0[j] - tg_xyy_yyyzzzzz_1[j] * fl1_fza);

                    tg_xxxyy_yyzzzzzz_0[j] = pb_x * tg_xxyy_yyzzzzzz_0[j] + fr * tg_xxyy_yyzzzzzz_1[j] + fl1_fx * (tg_xyy_yyzzzzzz_0[j] - tg_xyy_yyzzzzzz_1[j] * fl1_fza);

                    tg_xxxyy_yzzzzzzz_0[j] = pb_x * tg_xxyy_yzzzzzzz_0[j] + fr * tg_xxyy_yzzzzzzz_1[j] + fl1_fx * (tg_xyy_yzzzzzzz_0[j] - tg_xyy_yzzzzzzz_1[j] * fl1_fza);

                    tg_xxxyy_zzzzzzzz_0[j] = pb_x * tg_xxyy_zzzzzzzz_0[j] + fr * tg_xxyy_zzzzzzzz_1[j] + fl1_fx * (tg_xyy_zzzzzzzz_0[j] - tg_xyy_zzzzzzzz_1[j] * fl1_fza);

                    tg_xxxyz_xxxxxxxx_0[j] = pb_x * tg_xxyz_xxxxxxxx_0[j] + fr * tg_xxyz_xxxxxxxx_1[j] + fl1_fx * (tg_xyz_xxxxxxxx_0[j] - tg_xyz_xxxxxxxx_1[j] * fl1_fza) + 4.0 * fl1_fxn * tg_xxyz_xxxxxxx_1[j];

                    tg_xxxyz_xxxxxxxy_0[j] = pb_x * tg_xxyz_xxxxxxxy_0[j] + fr * tg_xxyz_xxxxxxxy_1[j] + fl1_fx * (tg_xyz_xxxxxxxy_0[j] - tg_xyz_xxxxxxxy_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_xxyz_xxxxxxy_1[j];

                    tg_xxxyz_xxxxxxxz_0[j] = pb_x * tg_xxyz_xxxxxxxz_0[j] + fr * tg_xxyz_xxxxxxxz_1[j] + fl1_fx * (tg_xyz_xxxxxxxz_0[j] - tg_xyz_xxxxxxxz_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_xxyz_xxxxxxz_1[j];

                    tg_xxxyz_xxxxxxyy_0[j] = pb_x * tg_xxyz_xxxxxxyy_0[j] + fr * tg_xxyz_xxxxxxyy_1[j] + fl1_fx * (tg_xyz_xxxxxxyy_0[j] - tg_xyz_xxxxxxyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxyz_xxxxxyy_1[j];

                    tg_xxxyz_xxxxxxyz_0[j] = pb_x * tg_xxyz_xxxxxxyz_0[j] + fr * tg_xxyz_xxxxxxyz_1[j] + fl1_fx * (tg_xyz_xxxxxxyz_0[j] - tg_xyz_xxxxxxyz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxyz_xxxxxyz_1[j];

                    tg_xxxyz_xxxxxxzz_0[j] = pb_x * tg_xxyz_xxxxxxzz_0[j] + fr * tg_xxyz_xxxxxxzz_1[j] + fl1_fx * (tg_xyz_xxxxxxzz_0[j] - tg_xyz_xxxxxxzz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxyz_xxxxxzz_1[j];

                    tg_xxxyz_xxxxxyyy_0[j] = pb_x * tg_xxyz_xxxxxyyy_0[j] + fr * tg_xxyz_xxxxxyyy_1[j] + fl1_fx * (tg_xyz_xxxxxyyy_0[j] - tg_xyz_xxxxxyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyz_xxxxyyy_1[j];

                    tg_xxxyz_xxxxxyyz_0[j] = pb_x * tg_xxyz_xxxxxyyz_0[j] + fr * tg_xxyz_xxxxxyyz_1[j] + fl1_fx * (tg_xyz_xxxxxyyz_0[j] - tg_xyz_xxxxxyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyz_xxxxyyz_1[j];

                    tg_xxxyz_xxxxxyzz_0[j] = pb_x * tg_xxyz_xxxxxyzz_0[j] + fr * tg_xxyz_xxxxxyzz_1[j] + fl1_fx * (tg_xyz_xxxxxyzz_0[j] - tg_xyz_xxxxxyzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyz_xxxxyzz_1[j];

                    tg_xxxyz_xxxxxzzz_0[j] = pb_x * tg_xxyz_xxxxxzzz_0[j] + fr * tg_xxyz_xxxxxzzz_1[j] + fl1_fx * (tg_xyz_xxxxxzzz_0[j] - tg_xyz_xxxxxzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyz_xxxxzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSL_190_285(      CMemBlock2D<double>* primBuffer,
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
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_8_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_8_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_8_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_8_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_8_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_7_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {7, -1, -1, -1}, 
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

                auto tg_xxyz_xxxxyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 190); 

                auto tg_xxyz_xxxxyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 191); 

                auto tg_xxyz_xxxxyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 192); 

                auto tg_xxyz_xxxxyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 193); 

                auto tg_xxyz_xxxxzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 194); 

                auto tg_xxyz_xxxyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 195); 

                auto tg_xxyz_xxxyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 196); 

                auto tg_xxyz_xxxyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 197); 

                auto tg_xxyz_xxxyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 198); 

                auto tg_xxyz_xxxyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 199); 

                auto tg_xxyz_xxxzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 200); 

                auto tg_xxyz_xxyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 201); 

                auto tg_xxyz_xxyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 202); 

                auto tg_xxyz_xxyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 203); 

                auto tg_xxyz_xxyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 204); 

                auto tg_xxyz_xxyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 205); 

                auto tg_xxyz_xxyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 206); 

                auto tg_xxyz_xxzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 207); 

                auto tg_xxyz_xyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 208); 

                auto tg_xxyz_xyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 209); 

                auto tg_xxyz_xyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 210); 

                auto tg_xxyz_xyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 211); 

                auto tg_xxyz_xyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 212); 

                auto tg_xxyz_xyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 213); 

                auto tg_xxyz_xyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 214); 

                auto tg_xxyz_xzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 215); 

                auto tg_xxyz_yyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 216); 

                auto tg_xxyz_yyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 217); 

                auto tg_xxyz_yyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 218); 

                auto tg_xxyz_yyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 219); 

                auto tg_xxyz_yyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 220); 

                auto tg_xxyz_yyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 221); 

                auto tg_xxyz_yyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 222); 

                auto tg_xxyz_yzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 223); 

                auto tg_xxyz_zzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 224); 

                auto tg_xxzz_xxxxxxxx_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 225); 

                auto tg_xxzz_xxxxxxxy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 226); 

                auto tg_xxzz_xxxxxxxz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 227); 

                auto tg_xxzz_xxxxxxyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 228); 

                auto tg_xxzz_xxxxxxyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 229); 

                auto tg_xxzz_xxxxxxzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 230); 

                auto tg_xxzz_xxxxxyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 231); 

                auto tg_xxzz_xxxxxyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 232); 

                auto tg_xxzz_xxxxxyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 233); 

                auto tg_xxzz_xxxxxzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 234); 

                auto tg_xxzz_xxxxyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 235); 

                auto tg_xxzz_xxxxyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 236); 

                auto tg_xxzz_xxxxyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 237); 

                auto tg_xxzz_xxxxyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 238); 

                auto tg_xxzz_xxxxzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 239); 

                auto tg_xxzz_xxxyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 240); 

                auto tg_xxzz_xxxyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 241); 

                auto tg_xxzz_xxxyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 242); 

                auto tg_xxzz_xxxyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 243); 

                auto tg_xxzz_xxxyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 244); 

                auto tg_xxzz_xxxzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 245); 

                auto tg_xxzz_xxyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 246); 

                auto tg_xxzz_xxyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 247); 

                auto tg_xxzz_xxyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 248); 

                auto tg_xxzz_xxyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 249); 

                auto tg_xxzz_xxyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 250); 

                auto tg_xxzz_xxyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 251); 

                auto tg_xxzz_xxzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 252); 

                auto tg_xxzz_xyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 253); 

                auto tg_xxzz_xyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 254); 

                auto tg_xxzz_xyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 255); 

                auto tg_xxzz_xyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 256); 

                auto tg_xxzz_xyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 257); 

                auto tg_xxzz_xyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 258); 

                auto tg_xxzz_xyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 259); 

                auto tg_xxzz_xzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 260); 

                auto tg_xxzz_yyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 261); 

                auto tg_xxzz_yyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 262); 

                auto tg_xxzz_yyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 263); 

                auto tg_xxzz_yyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 264); 

                auto tg_xxzz_yyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 265); 

                auto tg_xxzz_yyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 266); 

                auto tg_xxzz_yyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 267); 

                auto tg_xxzz_yzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 268); 

                auto tg_xxzz_zzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 269); 

                auto tg_xyyy_xxxxxxxx_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 270); 

                auto tg_xyyy_xxxxxxxy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 271); 

                auto tg_xyyy_xxxxxxxz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 272); 

                auto tg_xyyy_xxxxxxyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 273); 

                auto tg_xyyy_xxxxxxyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 274); 

                auto tg_xyyy_xxxxxxzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 275); 

                auto tg_xyyy_xxxxxyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 276); 

                auto tg_xyyy_xxxxxyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 277); 

                auto tg_xyyy_xxxxxyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 278); 

                auto tg_xyyy_xxxxxzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 279); 

                auto tg_xyyy_xxxxyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 280); 

                auto tg_xyyy_xxxxyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 281); 

                auto tg_xyyy_xxxxyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 282); 

                auto tg_xyyy_xxxxyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 283); 

                auto tg_xyyy_xxxxzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 284); 

                auto tg_xxyz_xxxxyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 190); 

                auto tg_xxyz_xxxxyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 191); 

                auto tg_xxyz_xxxxyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 192); 

                auto tg_xxyz_xxxxyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 193); 

                auto tg_xxyz_xxxxzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 194); 

                auto tg_xxyz_xxxyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 195); 

                auto tg_xxyz_xxxyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 196); 

                auto tg_xxyz_xxxyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 197); 

                auto tg_xxyz_xxxyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 198); 

                auto tg_xxyz_xxxyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 199); 

                auto tg_xxyz_xxxzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 200); 

                auto tg_xxyz_xxyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 201); 

                auto tg_xxyz_xxyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 202); 

                auto tg_xxyz_xxyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 203); 

                auto tg_xxyz_xxyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 204); 

                auto tg_xxyz_xxyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 205); 

                auto tg_xxyz_xxyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 206); 

                auto tg_xxyz_xxzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 207); 

                auto tg_xxyz_xyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 208); 

                auto tg_xxyz_xyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 209); 

                auto tg_xxyz_xyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 210); 

                auto tg_xxyz_xyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 211); 

                auto tg_xxyz_xyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 212); 

                auto tg_xxyz_xyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 213); 

                auto tg_xxyz_xyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 214); 

                auto tg_xxyz_xzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 215); 

                auto tg_xxyz_yyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 216); 

                auto tg_xxyz_yyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 217); 

                auto tg_xxyz_yyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 218); 

                auto tg_xxyz_yyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 219); 

                auto tg_xxyz_yyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 220); 

                auto tg_xxyz_yyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 221); 

                auto tg_xxyz_yyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 222); 

                auto tg_xxyz_yzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 223); 

                auto tg_xxyz_zzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 224); 

                auto tg_xxzz_xxxxxxxx_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 225); 

                auto tg_xxzz_xxxxxxxy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 226); 

                auto tg_xxzz_xxxxxxxz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 227); 

                auto tg_xxzz_xxxxxxyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 228); 

                auto tg_xxzz_xxxxxxyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 229); 

                auto tg_xxzz_xxxxxxzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 230); 

                auto tg_xxzz_xxxxxyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 231); 

                auto tg_xxzz_xxxxxyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 232); 

                auto tg_xxzz_xxxxxyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 233); 

                auto tg_xxzz_xxxxxzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 234); 

                auto tg_xxzz_xxxxyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 235); 

                auto tg_xxzz_xxxxyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 236); 

                auto tg_xxzz_xxxxyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 237); 

                auto tg_xxzz_xxxxyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 238); 

                auto tg_xxzz_xxxxzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 239); 

                auto tg_xxzz_xxxyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 240); 

                auto tg_xxzz_xxxyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 241); 

                auto tg_xxzz_xxxyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 242); 

                auto tg_xxzz_xxxyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 243); 

                auto tg_xxzz_xxxyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 244); 

                auto tg_xxzz_xxxzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 245); 

                auto tg_xxzz_xxyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 246); 

                auto tg_xxzz_xxyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 247); 

                auto tg_xxzz_xxyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 248); 

                auto tg_xxzz_xxyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 249); 

                auto tg_xxzz_xxyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 250); 

                auto tg_xxzz_xxyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 251); 

                auto tg_xxzz_xxzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 252); 

                auto tg_xxzz_xyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 253); 

                auto tg_xxzz_xyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 254); 

                auto tg_xxzz_xyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 255); 

                auto tg_xxzz_xyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 256); 

                auto tg_xxzz_xyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 257); 

                auto tg_xxzz_xyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 258); 

                auto tg_xxzz_xyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 259); 

                auto tg_xxzz_xzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 260); 

                auto tg_xxzz_yyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 261); 

                auto tg_xxzz_yyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 262); 

                auto tg_xxzz_yyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 263); 

                auto tg_xxzz_yyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 264); 

                auto tg_xxzz_yyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 265); 

                auto tg_xxzz_yyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 266); 

                auto tg_xxzz_yyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 267); 

                auto tg_xxzz_yzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 268); 

                auto tg_xxzz_zzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 269); 

                auto tg_xyyy_xxxxxxxx_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 270); 

                auto tg_xyyy_xxxxxxxy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 271); 

                auto tg_xyyy_xxxxxxxz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 272); 

                auto tg_xyyy_xxxxxxyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 273); 

                auto tg_xyyy_xxxxxxyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 274); 

                auto tg_xyyy_xxxxxxzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 275); 

                auto tg_xyyy_xxxxxyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 276); 

                auto tg_xyyy_xxxxxyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 277); 

                auto tg_xyyy_xxxxxyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 278); 

                auto tg_xyyy_xxxxxzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 279); 

                auto tg_xyyy_xxxxyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 280); 

                auto tg_xyyy_xxxxyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 281); 

                auto tg_xyyy_xxxxyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 282); 

                auto tg_xyyy_xxxxyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 283); 

                auto tg_xyyy_xxxxzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 284); 

                auto tg_xyz_xxxxyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 190); 

                auto tg_xyz_xxxxyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 191); 

                auto tg_xyz_xxxxyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 192); 

                auto tg_xyz_xxxxyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 193); 

                auto tg_xyz_xxxxzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 194); 

                auto tg_xyz_xxxyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 195); 

                auto tg_xyz_xxxyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 196); 

                auto tg_xyz_xxxyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 197); 

                auto tg_xyz_xxxyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 198); 

                auto tg_xyz_xxxyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 199); 

                auto tg_xyz_xxxzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 200); 

                auto tg_xyz_xxyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 201); 

                auto tg_xyz_xxyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 202); 

                auto tg_xyz_xxyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 203); 

                auto tg_xyz_xxyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 204); 

                auto tg_xyz_xxyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 205); 

                auto tg_xyz_xxyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 206); 

                auto tg_xyz_xxzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 207); 

                auto tg_xyz_xyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 208); 

                auto tg_xyz_xyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 209); 

                auto tg_xyz_xyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 210); 

                auto tg_xyz_xyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 211); 

                auto tg_xyz_xyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 212); 

                auto tg_xyz_xyyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 213); 

                auto tg_xyz_xyzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 214); 

                auto tg_xyz_xzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 215); 

                auto tg_xyz_yyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 216); 

                auto tg_xyz_yyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 217); 

                auto tg_xyz_yyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 218); 

                auto tg_xyz_yyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 219); 

                auto tg_xyz_yyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 220); 

                auto tg_xyz_yyyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 221); 

                auto tg_xyz_yyzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 222); 

                auto tg_xyz_yzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 223); 

                auto tg_xyz_zzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 224); 

                auto tg_xzz_xxxxxxxx_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 225); 

                auto tg_xzz_xxxxxxxy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 226); 

                auto tg_xzz_xxxxxxxz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 227); 

                auto tg_xzz_xxxxxxyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 228); 

                auto tg_xzz_xxxxxxyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 229); 

                auto tg_xzz_xxxxxxzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 230); 

                auto tg_xzz_xxxxxyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 231); 

                auto tg_xzz_xxxxxyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 232); 

                auto tg_xzz_xxxxxyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 233); 

                auto tg_xzz_xxxxxzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 234); 

                auto tg_xzz_xxxxyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 235); 

                auto tg_xzz_xxxxyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 236); 

                auto tg_xzz_xxxxyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 237); 

                auto tg_xzz_xxxxyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 238); 

                auto tg_xzz_xxxxzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 239); 

                auto tg_xzz_xxxyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 240); 

                auto tg_xzz_xxxyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 241); 

                auto tg_xzz_xxxyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 242); 

                auto tg_xzz_xxxyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 243); 

                auto tg_xzz_xxxyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 244); 

                auto tg_xzz_xxxzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 245); 

                auto tg_xzz_xxyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 246); 

                auto tg_xzz_xxyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 247); 

                auto tg_xzz_xxyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 248); 

                auto tg_xzz_xxyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 249); 

                auto tg_xzz_xxyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 250); 

                auto tg_xzz_xxyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 251); 

                auto tg_xzz_xxzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 252); 

                auto tg_xzz_xyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 253); 

                auto tg_xzz_xyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 254); 

                auto tg_xzz_xyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 255); 

                auto tg_xzz_xyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 256); 

                auto tg_xzz_xyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 257); 

                auto tg_xzz_xyyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 258); 

                auto tg_xzz_xyzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 259); 

                auto tg_xzz_xzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 260); 

                auto tg_xzz_yyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 261); 

                auto tg_xzz_yyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 262); 

                auto tg_xzz_yyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 263); 

                auto tg_xzz_yyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 264); 

                auto tg_xzz_yyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 265); 

                auto tg_xzz_yyyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 266); 

                auto tg_xzz_yyzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 267); 

                auto tg_xzz_yzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 268); 

                auto tg_xzz_zzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 269); 

                auto tg_yyy_xxxxxxxx_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 270); 

                auto tg_yyy_xxxxxxxy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 271); 

                auto tg_yyy_xxxxxxxz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 272); 

                auto tg_yyy_xxxxxxyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 273); 

                auto tg_yyy_xxxxxxyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 274); 

                auto tg_yyy_xxxxxxzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 275); 

                auto tg_yyy_xxxxxyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 276); 

                auto tg_yyy_xxxxxyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 277); 

                auto tg_yyy_xxxxxyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 278); 

                auto tg_yyy_xxxxxzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 279); 

                auto tg_yyy_xxxxyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 280); 

                auto tg_yyy_xxxxyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 281); 

                auto tg_yyy_xxxxyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 282); 

                auto tg_yyy_xxxxyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 283); 

                auto tg_yyy_xxxxzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 284); 

                auto tg_xyz_xxxxyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 190); 

                auto tg_xyz_xxxxyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 191); 

                auto tg_xyz_xxxxyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 192); 

                auto tg_xyz_xxxxyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 193); 

                auto tg_xyz_xxxxzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 194); 

                auto tg_xyz_xxxyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 195); 

                auto tg_xyz_xxxyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 196); 

                auto tg_xyz_xxxyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 197); 

                auto tg_xyz_xxxyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 198); 

                auto tg_xyz_xxxyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 199); 

                auto tg_xyz_xxxzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 200); 

                auto tg_xyz_xxyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 201); 

                auto tg_xyz_xxyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 202); 

                auto tg_xyz_xxyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 203); 

                auto tg_xyz_xxyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 204); 

                auto tg_xyz_xxyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 205); 

                auto tg_xyz_xxyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 206); 

                auto tg_xyz_xxzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 207); 

                auto tg_xyz_xyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 208); 

                auto tg_xyz_xyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 209); 

                auto tg_xyz_xyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 210); 

                auto tg_xyz_xyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 211); 

                auto tg_xyz_xyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 212); 

                auto tg_xyz_xyyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 213); 

                auto tg_xyz_xyzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 214); 

                auto tg_xyz_xzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 215); 

                auto tg_xyz_yyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 216); 

                auto tg_xyz_yyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 217); 

                auto tg_xyz_yyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 218); 

                auto tg_xyz_yyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 219); 

                auto tg_xyz_yyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 220); 

                auto tg_xyz_yyyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 221); 

                auto tg_xyz_yyzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 222); 

                auto tg_xyz_yzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 223); 

                auto tg_xyz_zzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 224); 

                auto tg_xzz_xxxxxxxx_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 225); 

                auto tg_xzz_xxxxxxxy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 226); 

                auto tg_xzz_xxxxxxxz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 227); 

                auto tg_xzz_xxxxxxyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 228); 

                auto tg_xzz_xxxxxxyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 229); 

                auto tg_xzz_xxxxxxzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 230); 

                auto tg_xzz_xxxxxyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 231); 

                auto tg_xzz_xxxxxyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 232); 

                auto tg_xzz_xxxxxyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 233); 

                auto tg_xzz_xxxxxzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 234); 

                auto tg_xzz_xxxxyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 235); 

                auto tg_xzz_xxxxyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 236); 

                auto tg_xzz_xxxxyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 237); 

                auto tg_xzz_xxxxyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 238); 

                auto tg_xzz_xxxxzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 239); 

                auto tg_xzz_xxxyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 240); 

                auto tg_xzz_xxxyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 241); 

                auto tg_xzz_xxxyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 242); 

                auto tg_xzz_xxxyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 243); 

                auto tg_xzz_xxxyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 244); 

                auto tg_xzz_xxxzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 245); 

                auto tg_xzz_xxyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 246); 

                auto tg_xzz_xxyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 247); 

                auto tg_xzz_xxyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 248); 

                auto tg_xzz_xxyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 249); 

                auto tg_xzz_xxyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 250); 

                auto tg_xzz_xxyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 251); 

                auto tg_xzz_xxzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 252); 

                auto tg_xzz_xyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 253); 

                auto tg_xzz_xyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 254); 

                auto tg_xzz_xyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 255); 

                auto tg_xzz_xyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 256); 

                auto tg_xzz_xyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 257); 

                auto tg_xzz_xyyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 258); 

                auto tg_xzz_xyzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 259); 

                auto tg_xzz_xzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 260); 

                auto tg_xzz_yyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 261); 

                auto tg_xzz_yyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 262); 

                auto tg_xzz_yyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 263); 

                auto tg_xzz_yyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 264); 

                auto tg_xzz_yyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 265); 

                auto tg_xzz_yyyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 266); 

                auto tg_xzz_yyzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 267); 

                auto tg_xzz_yzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 268); 

                auto tg_xzz_zzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 269); 

                auto tg_yyy_xxxxxxxx_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 270); 

                auto tg_yyy_xxxxxxxy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 271); 

                auto tg_yyy_xxxxxxxz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 272); 

                auto tg_yyy_xxxxxxyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 273); 

                auto tg_yyy_xxxxxxyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 274); 

                auto tg_yyy_xxxxxxzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 275); 

                auto tg_yyy_xxxxxyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 276); 

                auto tg_yyy_xxxxxyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 277); 

                auto tg_yyy_xxxxxyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 278); 

                auto tg_yyy_xxxxxzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 279); 

                auto tg_yyy_xxxxyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 280); 

                auto tg_yyy_xxxxyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 281); 

                auto tg_yyy_xxxxyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 282); 

                auto tg_yyy_xxxxyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 283); 

                auto tg_yyy_xxxxzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 284); 

                auto tg_xxyz_xxxyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 154); 

                auto tg_xxyz_xxxyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 155); 

                auto tg_xxyz_xxxyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 156); 

                auto tg_xxyz_xxxyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 157); 

                auto tg_xxyz_xxxzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 158); 

                auto tg_xxyz_xxyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 159); 

                auto tg_xxyz_xxyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 160); 

                auto tg_xxyz_xxyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 161); 

                auto tg_xxyz_xxyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 162); 

                auto tg_xxyz_xxyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 163); 

                auto tg_xxyz_xxzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 164); 

                auto tg_xxyz_xyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 165); 

                auto tg_xxyz_xyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 166); 

                auto tg_xxyz_xyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 167); 

                auto tg_xxyz_xyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 168); 

                auto tg_xxyz_xyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 169); 

                auto tg_xxyz_xyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 170); 

                auto tg_xxyz_xzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 171); 

                auto tg_xxyz_yyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 172); 

                auto tg_xxyz_yyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 173); 

                auto tg_xxyz_yyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 174); 

                auto tg_xxyz_yyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 175); 

                auto tg_xxyz_yyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 176); 

                auto tg_xxyz_yyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 177); 

                auto tg_xxyz_yzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 178); 

                auto tg_xxyz_zzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 179); 

                auto tg_xxzz_xxxxxxx_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 180); 

                auto tg_xxzz_xxxxxxy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 181); 

                auto tg_xxzz_xxxxxxz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 182); 

                auto tg_xxzz_xxxxxyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 183); 

                auto tg_xxzz_xxxxxyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 184); 

                auto tg_xxzz_xxxxxzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 185); 

                auto tg_xxzz_xxxxyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 186); 

                auto tg_xxzz_xxxxyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 187); 

                auto tg_xxzz_xxxxyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 188); 

                auto tg_xxzz_xxxxzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 189); 

                auto tg_xxzz_xxxyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 190); 

                auto tg_xxzz_xxxyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 191); 

                auto tg_xxzz_xxxyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 192); 

                auto tg_xxzz_xxxyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 193); 

                auto tg_xxzz_xxxzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 194); 

                auto tg_xxzz_xxyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 195); 

                auto tg_xxzz_xxyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 196); 

                auto tg_xxzz_xxyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 197); 

                auto tg_xxzz_xxyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 198); 

                auto tg_xxzz_xxyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 199); 

                auto tg_xxzz_xxzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 200); 

                auto tg_xxzz_xyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 201); 

                auto tg_xxzz_xyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 202); 

                auto tg_xxzz_xyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 203); 

                auto tg_xxzz_xyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 204); 

                auto tg_xxzz_xyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 205); 

                auto tg_xxzz_xyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 206); 

                auto tg_xxzz_xzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 207); 

                auto tg_xxzz_yyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 208); 

                auto tg_xxzz_yyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 209); 

                auto tg_xxzz_yyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 210); 

                auto tg_xxzz_yyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 211); 

                auto tg_xxzz_yyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 212); 

                auto tg_xxzz_yyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 213); 

                auto tg_xxzz_yzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 214); 

                auto tg_xxzz_zzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 215); 

                auto tg_xyyy_xxxxxxx_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 216); 

                auto tg_xyyy_xxxxxxy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 217); 

                auto tg_xyyy_xxxxxxz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 218); 

                auto tg_xyyy_xxxxxyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 219); 

                auto tg_xyyy_xxxxxyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 220); 

                auto tg_xyyy_xxxxxzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 221); 

                auto tg_xyyy_xxxxyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 222); 

                auto tg_xyyy_xxxxyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 223); 

                auto tg_xyyy_xxxxyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 224); 

                auto tg_xyyy_xxxxzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 225); 

                auto tg_xyyy_xxxyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 226); 

                auto tg_xyyy_xxxyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 227); 

                auto tg_xyyy_xxxyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 228); 

                auto tg_xyyy_xxxyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 229); 

                auto tg_xyyy_xxxzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 230); 

                // set up pointers to integrals

                auto tg_xxxyz_xxxxyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 190); 

                auto tg_xxxyz_xxxxyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 191); 

                auto tg_xxxyz_xxxxyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 192); 

                auto tg_xxxyz_xxxxyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 193); 

                auto tg_xxxyz_xxxxzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 194); 

                auto tg_xxxyz_xxxyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 195); 

                auto tg_xxxyz_xxxyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 196); 

                auto tg_xxxyz_xxxyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 197); 

                auto tg_xxxyz_xxxyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 198); 

                auto tg_xxxyz_xxxyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 199); 

                auto tg_xxxyz_xxxzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 200); 

                auto tg_xxxyz_xxyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 201); 

                auto tg_xxxyz_xxyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 202); 

                auto tg_xxxyz_xxyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 203); 

                auto tg_xxxyz_xxyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 204); 

                auto tg_xxxyz_xxyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 205); 

                auto tg_xxxyz_xxyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 206); 

                auto tg_xxxyz_xxzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 207); 

                auto tg_xxxyz_xyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 208); 

                auto tg_xxxyz_xyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 209); 

                auto tg_xxxyz_xyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 210); 

                auto tg_xxxyz_xyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 211); 

                auto tg_xxxyz_xyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 212); 

                auto tg_xxxyz_xyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 213); 

                auto tg_xxxyz_xyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 214); 

                auto tg_xxxyz_xzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 215); 

                auto tg_xxxyz_yyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 216); 

                auto tg_xxxyz_yyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 217); 

                auto tg_xxxyz_yyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 218); 

                auto tg_xxxyz_yyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 219); 

                auto tg_xxxyz_yyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 220); 

                auto tg_xxxyz_yyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 221); 

                auto tg_xxxyz_yyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 222); 

                auto tg_xxxyz_yzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 223); 

                auto tg_xxxyz_zzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 224); 

                auto tg_xxxzz_xxxxxxxx_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 225); 

                auto tg_xxxzz_xxxxxxxy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 226); 

                auto tg_xxxzz_xxxxxxxz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 227); 

                auto tg_xxxzz_xxxxxxyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 228); 

                auto tg_xxxzz_xxxxxxyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 229); 

                auto tg_xxxzz_xxxxxxzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 230); 

                auto tg_xxxzz_xxxxxyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 231); 

                auto tg_xxxzz_xxxxxyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 232); 

                auto tg_xxxzz_xxxxxyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 233); 

                auto tg_xxxzz_xxxxxzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 234); 

                auto tg_xxxzz_xxxxyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 235); 

                auto tg_xxxzz_xxxxyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 236); 

                auto tg_xxxzz_xxxxyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 237); 

                auto tg_xxxzz_xxxxyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 238); 

                auto tg_xxxzz_xxxxzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 239); 

                auto tg_xxxzz_xxxyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 240); 

                auto tg_xxxzz_xxxyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 241); 

                auto tg_xxxzz_xxxyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 242); 

                auto tg_xxxzz_xxxyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 243); 

                auto tg_xxxzz_xxxyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 244); 

                auto tg_xxxzz_xxxzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 245); 

                auto tg_xxxzz_xxyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 246); 

                auto tg_xxxzz_xxyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 247); 

                auto tg_xxxzz_xxyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 248); 

                auto tg_xxxzz_xxyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 249); 

                auto tg_xxxzz_xxyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 250); 

                auto tg_xxxzz_xxyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 251); 

                auto tg_xxxzz_xxzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 252); 

                auto tg_xxxzz_xyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 253); 

                auto tg_xxxzz_xyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 254); 

                auto tg_xxxzz_xyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 255); 

                auto tg_xxxzz_xyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 256); 

                auto tg_xxxzz_xyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 257); 

                auto tg_xxxzz_xyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 258); 

                auto tg_xxxzz_xyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 259); 

                auto tg_xxxzz_xzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 260); 

                auto tg_xxxzz_yyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 261); 

                auto tg_xxxzz_yyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 262); 

                auto tg_xxxzz_yyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 263); 

                auto tg_xxxzz_yyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 264); 

                auto tg_xxxzz_yyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 265); 

                auto tg_xxxzz_yyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 266); 

                auto tg_xxxzz_yyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 267); 

                auto tg_xxxzz_yzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 268); 

                auto tg_xxxzz_zzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 269); 

                auto tg_xxyyy_xxxxxxxx_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 270); 

                auto tg_xxyyy_xxxxxxxy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 271); 

                auto tg_xxyyy_xxxxxxxz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 272); 

                auto tg_xxyyy_xxxxxxyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 273); 

                auto tg_xxyyy_xxxxxxyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 274); 

                auto tg_xxyyy_xxxxxxzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 275); 

                auto tg_xxyyy_xxxxxyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 276); 

                auto tg_xxyyy_xxxxxyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 277); 

                auto tg_xxyyy_xxxxxyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 278); 

                auto tg_xxyyy_xxxxxzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 279); 

                auto tg_xxyyy_xxxxyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 280); 

                auto tg_xxyyy_xxxxyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 281); 

                auto tg_xxyyy_xxxxyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 282); 

                auto tg_xxyyy_xxxxyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 283); 

                auto tg_xxyyy_xxxxzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 284); 

                // Batch of Integrals (190,285)

                #pragma omp simd aligned(fxn, fza, tg_xxxyz_xxxxyyyy_0, tg_xxxyz_xxxxyyyz_0, \
                                         tg_xxxyz_xxxxyyzz_0, tg_xxxyz_xxxxyzzz_0, tg_xxxyz_xxxxzzzz_0, tg_xxxyz_xxxyyyyy_0, \
                                         tg_xxxyz_xxxyyyyz_0, tg_xxxyz_xxxyyyzz_0, tg_xxxyz_xxxyyzzz_0, tg_xxxyz_xxxyzzzz_0, \
                                         tg_xxxyz_xxxzzzzz_0, tg_xxxyz_xxyyyyyy_0, tg_xxxyz_xxyyyyyz_0, tg_xxxyz_xxyyyyzz_0, \
                                         tg_xxxyz_xxyyyzzz_0, tg_xxxyz_xxyyzzzz_0, tg_xxxyz_xxyzzzzz_0, tg_xxxyz_xxzzzzzz_0, \
                                         tg_xxxyz_xyyyyyyy_0, tg_xxxyz_xyyyyyyz_0, tg_xxxyz_xyyyyyzz_0, tg_xxxyz_xyyyyzzz_0, \
                                         tg_xxxyz_xyyyzzzz_0, tg_xxxyz_xyyzzzzz_0, tg_xxxyz_xyzzzzzz_0, tg_xxxyz_xzzzzzzz_0, \
                                         tg_xxxyz_yyyyyyyy_0, tg_xxxyz_yyyyyyyz_0, tg_xxxyz_yyyyyyzz_0, tg_xxxyz_yyyyyzzz_0, \
                                         tg_xxxyz_yyyyzzzz_0, tg_xxxyz_yyyzzzzz_0, tg_xxxyz_yyzzzzzz_0, tg_xxxyz_yzzzzzzz_0, \
                                         tg_xxxyz_zzzzzzzz_0, tg_xxxzz_xxxxxxxx_0, tg_xxxzz_xxxxxxxy_0, tg_xxxzz_xxxxxxxz_0, \
                                         tg_xxxzz_xxxxxxyy_0, tg_xxxzz_xxxxxxyz_0, tg_xxxzz_xxxxxxzz_0, tg_xxxzz_xxxxxyyy_0, \
                                         tg_xxxzz_xxxxxyyz_0, tg_xxxzz_xxxxxyzz_0, tg_xxxzz_xxxxxzzz_0, tg_xxxzz_xxxxyyyy_0, \
                                         tg_xxxzz_xxxxyyyz_0, tg_xxxzz_xxxxyyzz_0, tg_xxxzz_xxxxyzzz_0, tg_xxxzz_xxxxzzzz_0, \
                                         tg_xxxzz_xxxyyyyy_0, tg_xxxzz_xxxyyyyz_0, tg_xxxzz_xxxyyyzz_0, tg_xxxzz_xxxyyzzz_0, \
                                         tg_xxxzz_xxxyzzzz_0, tg_xxxzz_xxxzzzzz_0, tg_xxxzz_xxyyyyyy_0, tg_xxxzz_xxyyyyyz_0, \
                                         tg_xxxzz_xxyyyyzz_0, tg_xxxzz_xxyyyzzz_0, tg_xxxzz_xxyyzzzz_0, tg_xxxzz_xxyzzzzz_0, \
                                         tg_xxxzz_xxzzzzzz_0, tg_xxxzz_xyyyyyyy_0, tg_xxxzz_xyyyyyyz_0, tg_xxxzz_xyyyyyzz_0, \
                                         tg_xxxzz_xyyyyzzz_0, tg_xxxzz_xyyyzzzz_0, tg_xxxzz_xyyzzzzz_0, tg_xxxzz_xyzzzzzz_0, \
                                         tg_xxxzz_xzzzzzzz_0, tg_xxxzz_yyyyyyyy_0, tg_xxxzz_yyyyyyyz_0, tg_xxxzz_yyyyyyzz_0, \
                                         tg_xxxzz_yyyyyzzz_0, tg_xxxzz_yyyyzzzz_0, tg_xxxzz_yyyzzzzz_0, tg_xxxzz_yyzzzzzz_0, \
                                         tg_xxxzz_yzzzzzzz_0, tg_xxxzz_zzzzzzzz_0, tg_xxyyy_xxxxxxxx_0, tg_xxyyy_xxxxxxxy_0, \
                                         tg_xxyyy_xxxxxxxz_0, tg_xxyyy_xxxxxxyy_0, tg_xxyyy_xxxxxxyz_0, tg_xxyyy_xxxxxxzz_0, \
                                         tg_xxyyy_xxxxxyyy_0, tg_xxyyy_xxxxxyyz_0, tg_xxyyy_xxxxxyzz_0, tg_xxyyy_xxxxxzzz_0, \
                                         tg_xxyyy_xxxxyyyy_0, tg_xxyyy_xxxxyyyz_0, tg_xxyyy_xxxxyyzz_0, tg_xxyyy_xxxxyzzz_0, \
                                         tg_xxyyy_xxxxzzzz_0, tg_xxyz_xxxxyyyy_0, tg_xxyz_xxxxyyyy_1, tg_xxyz_xxxxyyyz_0, \
                                         tg_xxyz_xxxxyyyz_1, tg_xxyz_xxxxyyzz_0, tg_xxyz_xxxxyyzz_1, tg_xxyz_xxxxyzzz_0, \
                                         tg_xxyz_xxxxyzzz_1, tg_xxyz_xxxxzzzz_0, tg_xxyz_xxxxzzzz_1, tg_xxyz_xxxyyyy_1, \
                                         tg_xxyz_xxxyyyyy_0, tg_xxyz_xxxyyyyy_1, tg_xxyz_xxxyyyyz_0, tg_xxyz_xxxyyyyz_1, \
                                         tg_xxyz_xxxyyyz_1, tg_xxyz_xxxyyyzz_0, tg_xxyz_xxxyyyzz_1, tg_xxyz_xxxyyzz_1, \
                                         tg_xxyz_xxxyyzzz_0, tg_xxyz_xxxyyzzz_1, tg_xxyz_xxxyzzz_1, tg_xxyz_xxxyzzzz_0, \
                                         tg_xxyz_xxxyzzzz_1, tg_xxyz_xxxzzzz_1, tg_xxyz_xxxzzzzz_0, tg_xxyz_xxxzzzzz_1, \
                                         tg_xxyz_xxyyyyy_1, tg_xxyz_xxyyyyyy_0, tg_xxyz_xxyyyyyy_1, tg_xxyz_xxyyyyyz_0, \
                                         tg_xxyz_xxyyyyyz_1, tg_xxyz_xxyyyyz_1, tg_xxyz_xxyyyyzz_0, tg_xxyz_xxyyyyzz_1, \
                                         tg_xxyz_xxyyyzz_1, tg_xxyz_xxyyyzzz_0, tg_xxyz_xxyyyzzz_1, tg_xxyz_xxyyzzz_1, \
                                         tg_xxyz_xxyyzzzz_0, tg_xxyz_xxyyzzzz_1, tg_xxyz_xxyzzzz_1, tg_xxyz_xxyzzzzz_0, \
                                         tg_xxyz_xxyzzzzz_1, tg_xxyz_xxzzzzz_1, tg_xxyz_xxzzzzzz_0, tg_xxyz_xxzzzzzz_1, \
                                         tg_xxyz_xyyyyyy_1, tg_xxyz_xyyyyyyy_0, tg_xxyz_xyyyyyyy_1, tg_xxyz_xyyyyyyz_0, \
                                         tg_xxyz_xyyyyyyz_1, tg_xxyz_xyyyyyz_1, tg_xxyz_xyyyyyzz_0, tg_xxyz_xyyyyyzz_1, \
                                         tg_xxyz_xyyyyzz_1, tg_xxyz_xyyyyzzz_0, tg_xxyz_xyyyyzzz_1, tg_xxyz_xyyyzzz_1, \
                                         tg_xxyz_xyyyzzzz_0, tg_xxyz_xyyyzzzz_1, tg_xxyz_xyyzzzz_1, tg_xxyz_xyyzzzzz_0, \
                                         tg_xxyz_xyyzzzzz_1, tg_xxyz_xyzzzzz_1, tg_xxyz_xyzzzzzz_0, tg_xxyz_xyzzzzzz_1, \
                                         tg_xxyz_xzzzzzz_1, tg_xxyz_xzzzzzzz_0, tg_xxyz_xzzzzzzz_1, tg_xxyz_yyyyyyy_1, \
                                         tg_xxyz_yyyyyyyy_0, tg_xxyz_yyyyyyyy_1, tg_xxyz_yyyyyyyz_0, tg_xxyz_yyyyyyyz_1, \
                                         tg_xxyz_yyyyyyz_1, tg_xxyz_yyyyyyzz_0, tg_xxyz_yyyyyyzz_1, tg_xxyz_yyyyyzz_1, \
                                         tg_xxyz_yyyyyzzz_0, tg_xxyz_yyyyyzzz_1, tg_xxyz_yyyyzzz_1, tg_xxyz_yyyyzzzz_0, \
                                         tg_xxyz_yyyyzzzz_1, tg_xxyz_yyyzzzz_1, tg_xxyz_yyyzzzzz_0, tg_xxyz_yyyzzzzz_1, \
                                         tg_xxyz_yyzzzzz_1, tg_xxyz_yyzzzzzz_0, tg_xxyz_yyzzzzzz_1, tg_xxyz_yzzzzzz_1, \
                                         tg_xxyz_yzzzzzzz_0, tg_xxyz_yzzzzzzz_1, tg_xxyz_zzzzzzz_1, tg_xxyz_zzzzzzzz_0, \
                                         tg_xxyz_zzzzzzzz_1, tg_xxzz_xxxxxxx_1, tg_xxzz_xxxxxxxx_0, tg_xxzz_xxxxxxxx_1, \
                                         tg_xxzz_xxxxxxxy_0, tg_xxzz_xxxxxxxy_1, tg_xxzz_xxxxxxxz_0, tg_xxzz_xxxxxxxz_1, \
                                         tg_xxzz_xxxxxxy_1, tg_xxzz_xxxxxxyy_0, tg_xxzz_xxxxxxyy_1, tg_xxzz_xxxxxxyz_0, \
                                         tg_xxzz_xxxxxxyz_1, tg_xxzz_xxxxxxz_1, tg_xxzz_xxxxxxzz_0, tg_xxzz_xxxxxxzz_1, \
                                         tg_xxzz_xxxxxyy_1, tg_xxzz_xxxxxyyy_0, tg_xxzz_xxxxxyyy_1, tg_xxzz_xxxxxyyz_0, \
                                         tg_xxzz_xxxxxyyz_1, tg_xxzz_xxxxxyz_1, tg_xxzz_xxxxxyzz_0, tg_xxzz_xxxxxyzz_1, \
                                         tg_xxzz_xxxxxzz_1, tg_xxzz_xxxxxzzz_0, tg_xxzz_xxxxxzzz_1, tg_xxzz_xxxxyyy_1, \
                                         tg_xxzz_xxxxyyyy_0, tg_xxzz_xxxxyyyy_1, tg_xxzz_xxxxyyyz_0, tg_xxzz_xxxxyyyz_1, \
                                         tg_xxzz_xxxxyyz_1, tg_xxzz_xxxxyyzz_0, tg_xxzz_xxxxyyzz_1, tg_xxzz_xxxxyzz_1, \
                                         tg_xxzz_xxxxyzzz_0, tg_xxzz_xxxxyzzz_1, tg_xxzz_xxxxzzz_1, tg_xxzz_xxxxzzzz_0, \
                                         tg_xxzz_xxxxzzzz_1, tg_xxzz_xxxyyyy_1, tg_xxzz_xxxyyyyy_0, tg_xxzz_xxxyyyyy_1, \
                                         tg_xxzz_xxxyyyyz_0, tg_xxzz_xxxyyyyz_1, tg_xxzz_xxxyyyz_1, tg_xxzz_xxxyyyzz_0, \
                                         tg_xxzz_xxxyyyzz_1, tg_xxzz_xxxyyzz_1, tg_xxzz_xxxyyzzz_0, tg_xxzz_xxxyyzzz_1, \
                                         tg_xxzz_xxxyzzz_1, tg_xxzz_xxxyzzzz_0, tg_xxzz_xxxyzzzz_1, tg_xxzz_xxxzzzz_1, \
                                         tg_xxzz_xxxzzzzz_0, tg_xxzz_xxxzzzzz_1, tg_xxzz_xxyyyyy_1, tg_xxzz_xxyyyyyy_0, \
                                         tg_xxzz_xxyyyyyy_1, tg_xxzz_xxyyyyyz_0, tg_xxzz_xxyyyyyz_1, tg_xxzz_xxyyyyz_1, \
                                         tg_xxzz_xxyyyyzz_0, tg_xxzz_xxyyyyzz_1, tg_xxzz_xxyyyzz_1, tg_xxzz_xxyyyzzz_0, \
                                         tg_xxzz_xxyyyzzz_1, tg_xxzz_xxyyzzz_1, tg_xxzz_xxyyzzzz_0, tg_xxzz_xxyyzzzz_1, \
                                         tg_xxzz_xxyzzzz_1, tg_xxzz_xxyzzzzz_0, tg_xxzz_xxyzzzzz_1, tg_xxzz_xxzzzzz_1, \
                                         tg_xxzz_xxzzzzzz_0, tg_xxzz_xxzzzzzz_1, tg_xxzz_xyyyyyy_1, tg_xxzz_xyyyyyyy_0, \
                                         tg_xxzz_xyyyyyyy_1, tg_xxzz_xyyyyyyz_0, tg_xxzz_xyyyyyyz_1, tg_xxzz_xyyyyyz_1, \
                                         tg_xxzz_xyyyyyzz_0, tg_xxzz_xyyyyyzz_1, tg_xxzz_xyyyyzz_1, tg_xxzz_xyyyyzzz_0, \
                                         tg_xxzz_xyyyyzzz_1, tg_xxzz_xyyyzzz_1, tg_xxzz_xyyyzzzz_0, tg_xxzz_xyyyzzzz_1, \
                                         tg_xxzz_xyyzzzz_1, tg_xxzz_xyyzzzzz_0, tg_xxzz_xyyzzzzz_1, tg_xxzz_xyzzzzz_1, \
                                         tg_xxzz_xyzzzzzz_0, tg_xxzz_xyzzzzzz_1, tg_xxzz_xzzzzzz_1, tg_xxzz_xzzzzzzz_0, \
                                         tg_xxzz_xzzzzzzz_1, tg_xxzz_yyyyyyy_1, tg_xxzz_yyyyyyyy_0, tg_xxzz_yyyyyyyy_1, \
                                         tg_xxzz_yyyyyyyz_0, tg_xxzz_yyyyyyyz_1, tg_xxzz_yyyyyyz_1, tg_xxzz_yyyyyyzz_0, \
                                         tg_xxzz_yyyyyyzz_1, tg_xxzz_yyyyyzz_1, tg_xxzz_yyyyyzzz_0, tg_xxzz_yyyyyzzz_1, \
                                         tg_xxzz_yyyyzzz_1, tg_xxzz_yyyyzzzz_0, tg_xxzz_yyyyzzzz_1, tg_xxzz_yyyzzzz_1, \
                                         tg_xxzz_yyyzzzzz_0, tg_xxzz_yyyzzzzz_1, tg_xxzz_yyzzzzz_1, tg_xxzz_yyzzzzzz_0, \
                                         tg_xxzz_yyzzzzzz_1, tg_xxzz_yzzzzzz_1, tg_xxzz_yzzzzzzz_0, tg_xxzz_yzzzzzzz_1, \
                                         tg_xxzz_zzzzzzz_1, tg_xxzz_zzzzzzzz_0, tg_xxzz_zzzzzzzz_1, tg_xyyy_xxxxxxx_1, \
                                         tg_xyyy_xxxxxxxx_0, tg_xyyy_xxxxxxxx_1, tg_xyyy_xxxxxxxy_0, tg_xyyy_xxxxxxxy_1, \
                                         tg_xyyy_xxxxxxxz_0, tg_xyyy_xxxxxxxz_1, tg_xyyy_xxxxxxy_1, tg_xyyy_xxxxxxyy_0, \
                                         tg_xyyy_xxxxxxyy_1, tg_xyyy_xxxxxxyz_0, tg_xyyy_xxxxxxyz_1, tg_xyyy_xxxxxxz_1, \
                                         tg_xyyy_xxxxxxzz_0, tg_xyyy_xxxxxxzz_1, tg_xyyy_xxxxxyy_1, tg_xyyy_xxxxxyyy_0, \
                                         tg_xyyy_xxxxxyyy_1, tg_xyyy_xxxxxyyz_0, tg_xyyy_xxxxxyyz_1, tg_xyyy_xxxxxyz_1, \
                                         tg_xyyy_xxxxxyzz_0, tg_xyyy_xxxxxyzz_1, tg_xyyy_xxxxxzz_1, tg_xyyy_xxxxxzzz_0, \
                                         tg_xyyy_xxxxxzzz_1, tg_xyyy_xxxxyyy_1, tg_xyyy_xxxxyyyy_0, tg_xyyy_xxxxyyyy_1, \
                                         tg_xyyy_xxxxyyyz_0, tg_xyyy_xxxxyyyz_1, tg_xyyy_xxxxyyz_1, tg_xyyy_xxxxyyzz_0, \
                                         tg_xyyy_xxxxyyzz_1, tg_xyyy_xxxxyzz_1, tg_xyyy_xxxxyzzz_0, tg_xyyy_xxxxyzzz_1, \
                                         tg_xyyy_xxxxzzz_1, tg_xyyy_xxxxzzzz_0, tg_xyyy_xxxxzzzz_1, tg_xyyy_xxxyyyy_1, \
                                         tg_xyyy_xxxyyyz_1, tg_xyyy_xxxyyzz_1, tg_xyyy_xxxyzzz_1, tg_xyyy_xxxzzzz_1, \
                                         tg_xyz_xxxxyyyy_0, tg_xyz_xxxxyyyy_1, tg_xyz_xxxxyyyz_0, tg_xyz_xxxxyyyz_1, \
                                         tg_xyz_xxxxyyzz_0, tg_xyz_xxxxyyzz_1, tg_xyz_xxxxyzzz_0, tg_xyz_xxxxyzzz_1, \
                                         tg_xyz_xxxxzzzz_0, tg_xyz_xxxxzzzz_1, tg_xyz_xxxyyyyy_0, tg_xyz_xxxyyyyy_1, \
                                         tg_xyz_xxxyyyyz_0, tg_xyz_xxxyyyyz_1, tg_xyz_xxxyyyzz_0, tg_xyz_xxxyyyzz_1, \
                                         tg_xyz_xxxyyzzz_0, tg_xyz_xxxyyzzz_1, tg_xyz_xxxyzzzz_0, tg_xyz_xxxyzzzz_1, \
                                         tg_xyz_xxxzzzzz_0, tg_xyz_xxxzzzzz_1, tg_xyz_xxyyyyyy_0, tg_xyz_xxyyyyyy_1, \
                                         tg_xyz_xxyyyyyz_0, tg_xyz_xxyyyyyz_1, tg_xyz_xxyyyyzz_0, tg_xyz_xxyyyyzz_1, \
                                         tg_xyz_xxyyyzzz_0, tg_xyz_xxyyyzzz_1, tg_xyz_xxyyzzzz_0, tg_xyz_xxyyzzzz_1, \
                                         tg_xyz_xxyzzzzz_0, tg_xyz_xxyzzzzz_1, tg_xyz_xxzzzzzz_0, tg_xyz_xxzzzzzz_1, \
                                         tg_xyz_xyyyyyyy_0, tg_xyz_xyyyyyyy_1, tg_xyz_xyyyyyyz_0, tg_xyz_xyyyyyyz_1, \
                                         tg_xyz_xyyyyyzz_0, tg_xyz_xyyyyyzz_1, tg_xyz_xyyyyzzz_0, tg_xyz_xyyyyzzz_1, \
                                         tg_xyz_xyyyzzzz_0, tg_xyz_xyyyzzzz_1, tg_xyz_xyyzzzzz_0, tg_xyz_xyyzzzzz_1, \
                                         tg_xyz_xyzzzzzz_0, tg_xyz_xyzzzzzz_1, tg_xyz_xzzzzzzz_0, tg_xyz_xzzzzzzz_1, \
                                         tg_xyz_yyyyyyyy_0, tg_xyz_yyyyyyyy_1, tg_xyz_yyyyyyyz_0, tg_xyz_yyyyyyyz_1, \
                                         tg_xyz_yyyyyyzz_0, tg_xyz_yyyyyyzz_1, tg_xyz_yyyyyzzz_0, tg_xyz_yyyyyzzz_1, \
                                         tg_xyz_yyyyzzzz_0, tg_xyz_yyyyzzzz_1, tg_xyz_yyyzzzzz_0, tg_xyz_yyyzzzzz_1, \
                                         tg_xyz_yyzzzzzz_0, tg_xyz_yyzzzzzz_1, tg_xyz_yzzzzzzz_0, tg_xyz_yzzzzzzz_1, \
                                         tg_xyz_zzzzzzzz_0, tg_xyz_zzzzzzzz_1, tg_xzz_xxxxxxxx_0, tg_xzz_xxxxxxxx_1, \
                                         tg_xzz_xxxxxxxy_0, tg_xzz_xxxxxxxy_1, tg_xzz_xxxxxxxz_0, tg_xzz_xxxxxxxz_1, \
                                         tg_xzz_xxxxxxyy_0, tg_xzz_xxxxxxyy_1, tg_xzz_xxxxxxyz_0, tg_xzz_xxxxxxyz_1, \
                                         tg_xzz_xxxxxxzz_0, tg_xzz_xxxxxxzz_1, tg_xzz_xxxxxyyy_0, tg_xzz_xxxxxyyy_1, \
                                         tg_xzz_xxxxxyyz_0, tg_xzz_xxxxxyyz_1, tg_xzz_xxxxxyzz_0, tg_xzz_xxxxxyzz_1, \
                                         tg_xzz_xxxxxzzz_0, tg_xzz_xxxxxzzz_1, tg_xzz_xxxxyyyy_0, tg_xzz_xxxxyyyy_1, \
                                         tg_xzz_xxxxyyyz_0, tg_xzz_xxxxyyyz_1, tg_xzz_xxxxyyzz_0, tg_xzz_xxxxyyzz_1, \
                                         tg_xzz_xxxxyzzz_0, tg_xzz_xxxxyzzz_1, tg_xzz_xxxxzzzz_0, tg_xzz_xxxxzzzz_1, \
                                         tg_xzz_xxxyyyyy_0, tg_xzz_xxxyyyyy_1, tg_xzz_xxxyyyyz_0, tg_xzz_xxxyyyyz_1, \
                                         tg_xzz_xxxyyyzz_0, tg_xzz_xxxyyyzz_1, tg_xzz_xxxyyzzz_0, tg_xzz_xxxyyzzz_1, \
                                         tg_xzz_xxxyzzzz_0, tg_xzz_xxxyzzzz_1, tg_xzz_xxxzzzzz_0, tg_xzz_xxxzzzzz_1, \
                                         tg_xzz_xxyyyyyy_0, tg_xzz_xxyyyyyy_1, tg_xzz_xxyyyyyz_0, tg_xzz_xxyyyyyz_1, \
                                         tg_xzz_xxyyyyzz_0, tg_xzz_xxyyyyzz_1, tg_xzz_xxyyyzzz_0, tg_xzz_xxyyyzzz_1, \
                                         tg_xzz_xxyyzzzz_0, tg_xzz_xxyyzzzz_1, tg_xzz_xxyzzzzz_0, tg_xzz_xxyzzzzz_1, \
                                         tg_xzz_xxzzzzzz_0, tg_xzz_xxzzzzzz_1, tg_xzz_xyyyyyyy_0, tg_xzz_xyyyyyyy_1, \
                                         tg_xzz_xyyyyyyz_0, tg_xzz_xyyyyyyz_1, tg_xzz_xyyyyyzz_0, tg_xzz_xyyyyyzz_1, \
                                         tg_xzz_xyyyyzzz_0, tg_xzz_xyyyyzzz_1, tg_xzz_xyyyzzzz_0, tg_xzz_xyyyzzzz_1, \
                                         tg_xzz_xyyzzzzz_0, tg_xzz_xyyzzzzz_1, tg_xzz_xyzzzzzz_0, tg_xzz_xyzzzzzz_1, \
                                         tg_xzz_xzzzzzzz_0, tg_xzz_xzzzzzzz_1, tg_xzz_yyyyyyyy_0, tg_xzz_yyyyyyyy_1, \
                                         tg_xzz_yyyyyyyz_0, tg_xzz_yyyyyyyz_1, tg_xzz_yyyyyyzz_0, tg_xzz_yyyyyyzz_1, \
                                         tg_xzz_yyyyyzzz_0, tg_xzz_yyyyyzzz_1, tg_xzz_yyyyzzzz_0, tg_xzz_yyyyzzzz_1, \
                                         tg_xzz_yyyzzzzz_0, tg_xzz_yyyzzzzz_1, tg_xzz_yyzzzzzz_0, tg_xzz_yyzzzzzz_1, \
                                         tg_xzz_yzzzzzzz_0, tg_xzz_yzzzzzzz_1, tg_xzz_zzzzzzzz_0, tg_xzz_zzzzzzzz_1, \
                                         tg_yyy_xxxxxxxx_0, tg_yyy_xxxxxxxx_1, tg_yyy_xxxxxxxy_0, tg_yyy_xxxxxxxy_1, \
                                         tg_yyy_xxxxxxxz_0, tg_yyy_xxxxxxxz_1, tg_yyy_xxxxxxyy_0, tg_yyy_xxxxxxyy_1, \
                                         tg_yyy_xxxxxxyz_0, tg_yyy_xxxxxxyz_1, tg_yyy_xxxxxxzz_0, tg_yyy_xxxxxxzz_1, \
                                         tg_yyy_xxxxxyyy_0, tg_yyy_xxxxxyyy_1, tg_yyy_xxxxxyyz_0, tg_yyy_xxxxxyyz_1, \
                                         tg_yyy_xxxxxyzz_0, tg_yyy_xxxxxyzz_1, tg_yyy_xxxxxzzz_0, tg_yyy_xxxxxzzz_1, \
                                         tg_yyy_xxxxyyyy_0, tg_yyy_xxxxyyyy_1, tg_yyy_xxxxyyyz_0, tg_yyy_xxxxyyyz_1, \
                                         tg_yyy_xxxxyyzz_0, tg_yyy_xxxxyyzz_1, tg_yyy_xxxxyzzz_0, tg_yyy_xxxxyzzz_1, \
                                         tg_yyy_xxxxzzzz_0, tg_yyy_xxxxzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxyz_xxxxyyyy_0[j] = pb_x * tg_xxyz_xxxxyyyy_0[j] + fr * tg_xxyz_xxxxyyyy_1[j] + fl1_fx * (tg_xyz_xxxxyyyy_0[j] - tg_xyz_xxxxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyz_xxxyyyy_1[j];

                    tg_xxxyz_xxxxyyyz_0[j] = pb_x * tg_xxyz_xxxxyyyz_0[j] + fr * tg_xxyz_xxxxyyyz_1[j] + fl1_fx * (tg_xyz_xxxxyyyz_0[j] - tg_xyz_xxxxyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyz_xxxyyyz_1[j];

                    tg_xxxyz_xxxxyyzz_0[j] = pb_x * tg_xxyz_xxxxyyzz_0[j] + fr * tg_xxyz_xxxxyyzz_1[j] + fl1_fx * (tg_xyz_xxxxyyzz_0[j] - tg_xyz_xxxxyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyz_xxxyyzz_1[j];

                    tg_xxxyz_xxxxyzzz_0[j] = pb_x * tg_xxyz_xxxxyzzz_0[j] + fr * tg_xxyz_xxxxyzzz_1[j] + fl1_fx * (tg_xyz_xxxxyzzz_0[j] - tg_xyz_xxxxyzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyz_xxxyzzz_1[j];

                    tg_xxxyz_xxxxzzzz_0[j] = pb_x * tg_xxyz_xxxxzzzz_0[j] + fr * tg_xxyz_xxxxzzzz_1[j] + fl1_fx * (tg_xyz_xxxxzzzz_0[j] - tg_xyz_xxxxzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyz_xxxzzzz_1[j];

                    tg_xxxyz_xxxyyyyy_0[j] = pb_x * tg_xxyz_xxxyyyyy_0[j] + fr * tg_xxyz_xxxyyyyy_1[j] + fl1_fx * (tg_xyz_xxxyyyyy_0[j] - tg_xyz_xxxyyyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyz_xxyyyyy_1[j];

                    tg_xxxyz_xxxyyyyz_0[j] = pb_x * tg_xxyz_xxxyyyyz_0[j] + fr * tg_xxyz_xxxyyyyz_1[j] + fl1_fx * (tg_xyz_xxxyyyyz_0[j] - tg_xyz_xxxyyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyz_xxyyyyz_1[j];

                    tg_xxxyz_xxxyyyzz_0[j] = pb_x * tg_xxyz_xxxyyyzz_0[j] + fr * tg_xxyz_xxxyyyzz_1[j] + fl1_fx * (tg_xyz_xxxyyyzz_0[j] - tg_xyz_xxxyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyz_xxyyyzz_1[j];

                    tg_xxxyz_xxxyyzzz_0[j] = pb_x * tg_xxyz_xxxyyzzz_0[j] + fr * tg_xxyz_xxxyyzzz_1[j] + fl1_fx * (tg_xyz_xxxyyzzz_0[j] - tg_xyz_xxxyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyz_xxyyzzz_1[j];

                    tg_xxxyz_xxxyzzzz_0[j] = pb_x * tg_xxyz_xxxyzzzz_0[j] + fr * tg_xxyz_xxxyzzzz_1[j] + fl1_fx * (tg_xyz_xxxyzzzz_0[j] - tg_xyz_xxxyzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyz_xxyzzzz_1[j];

                    tg_xxxyz_xxxzzzzz_0[j] = pb_x * tg_xxyz_xxxzzzzz_0[j] + fr * tg_xxyz_xxxzzzzz_1[j] + fl1_fx * (tg_xyz_xxxzzzzz_0[j] - tg_xyz_xxxzzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyz_xxzzzzz_1[j];

                    tg_xxxyz_xxyyyyyy_0[j] = pb_x * tg_xxyz_xxyyyyyy_0[j] + fr * tg_xxyz_xxyyyyyy_1[j] + fl1_fx * (tg_xyz_xxyyyyyy_0[j] - tg_xyz_xxyyyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxyz_xyyyyyy_1[j];

                    tg_xxxyz_xxyyyyyz_0[j] = pb_x * tg_xxyz_xxyyyyyz_0[j] + fr * tg_xxyz_xxyyyyyz_1[j] + fl1_fx * (tg_xyz_xxyyyyyz_0[j] - tg_xyz_xxyyyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxyz_xyyyyyz_1[j];

                    tg_xxxyz_xxyyyyzz_0[j] = pb_x * tg_xxyz_xxyyyyzz_0[j] + fr * tg_xxyz_xxyyyyzz_1[j] + fl1_fx * (tg_xyz_xxyyyyzz_0[j] - tg_xyz_xxyyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyz_xyyyyzz_1[j];

                    tg_xxxyz_xxyyyzzz_0[j] = pb_x * tg_xxyz_xxyyyzzz_0[j] + fr * tg_xxyz_xxyyyzzz_1[j] + fl1_fx * (tg_xyz_xxyyyzzz_0[j] - tg_xyz_xxyyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyz_xyyyzzz_1[j];

                    tg_xxxyz_xxyyzzzz_0[j] = pb_x * tg_xxyz_xxyyzzzz_0[j] + fr * tg_xxyz_xxyyzzzz_1[j] + fl1_fx * (tg_xyz_xxyyzzzz_0[j] - tg_xyz_xxyyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyz_xyyzzzz_1[j];

                    tg_xxxyz_xxyzzzzz_0[j] = pb_x * tg_xxyz_xxyzzzzz_0[j] + fr * tg_xxyz_xxyzzzzz_1[j] + fl1_fx * (tg_xyz_xxyzzzzz_0[j] - tg_xyz_xxyzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyz_xyzzzzz_1[j];

                    tg_xxxyz_xxzzzzzz_0[j] = pb_x * tg_xxyz_xxzzzzzz_0[j] + fr * tg_xxyz_xxzzzzzz_1[j] + fl1_fx * (tg_xyz_xxzzzzzz_0[j] - tg_xyz_xxzzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyz_xzzzzzz_1[j];

                    tg_xxxyz_xyyyyyyy_0[j] = pb_x * tg_xxyz_xyyyyyyy_0[j] + fr * tg_xxyz_xyyyyyyy_1[j] + fl1_fx * (tg_xyz_xyyyyyyy_0[j] - tg_xyz_xyyyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyz_yyyyyyy_1[j];

                    tg_xxxyz_xyyyyyyz_0[j] = pb_x * tg_xxyz_xyyyyyyz_0[j] + fr * tg_xxyz_xyyyyyyz_1[j] + fl1_fx * (tg_xyz_xyyyyyyz_0[j] - tg_xyz_xyyyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyz_yyyyyyz_1[j];

                    tg_xxxyz_xyyyyyzz_0[j] = pb_x * tg_xxyz_xyyyyyzz_0[j] + fr * tg_xxyz_xyyyyyzz_1[j] + fl1_fx * (tg_xyz_xyyyyyzz_0[j] - tg_xyz_xyyyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyz_yyyyyzz_1[j];

                    tg_xxxyz_xyyyyzzz_0[j] = pb_x * tg_xxyz_xyyyyzzz_0[j] + fr * tg_xxyz_xyyyyzzz_1[j] + fl1_fx * (tg_xyz_xyyyyzzz_0[j] - tg_xyz_xyyyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyz_yyyyzzz_1[j];

                    tg_xxxyz_xyyyzzzz_0[j] = pb_x * tg_xxyz_xyyyzzzz_0[j] + fr * tg_xxyz_xyyyzzzz_1[j] + fl1_fx * (tg_xyz_xyyyzzzz_0[j] - tg_xyz_xyyyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyz_yyyzzzz_1[j];

                    tg_xxxyz_xyyzzzzz_0[j] = pb_x * tg_xxyz_xyyzzzzz_0[j] + fr * tg_xxyz_xyyzzzzz_1[j] + fl1_fx * (tg_xyz_xyyzzzzz_0[j] - tg_xyz_xyyzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyz_yyzzzzz_1[j];

                    tg_xxxyz_xyzzzzzz_0[j] = pb_x * tg_xxyz_xyzzzzzz_0[j] + fr * tg_xxyz_xyzzzzzz_1[j] + fl1_fx * (tg_xyz_xyzzzzzz_0[j] - tg_xyz_xyzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyz_yzzzzzz_1[j];

                    tg_xxxyz_xzzzzzzz_0[j] = pb_x * tg_xxyz_xzzzzzzz_0[j] + fr * tg_xxyz_xzzzzzzz_1[j] + fl1_fx * (tg_xyz_xzzzzzzz_0[j] - tg_xyz_xzzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyz_zzzzzzz_1[j];

                    tg_xxxyz_yyyyyyyy_0[j] = pb_x * tg_xxyz_yyyyyyyy_0[j] + fr * tg_xxyz_yyyyyyyy_1[j] + fl1_fx * (tg_xyz_yyyyyyyy_0[j] - tg_xyz_yyyyyyyy_1[j] * fl1_fza);

                    tg_xxxyz_yyyyyyyz_0[j] = pb_x * tg_xxyz_yyyyyyyz_0[j] + fr * tg_xxyz_yyyyyyyz_1[j] + fl1_fx * (tg_xyz_yyyyyyyz_0[j] - tg_xyz_yyyyyyyz_1[j] * fl1_fza);

                    tg_xxxyz_yyyyyyzz_0[j] = pb_x * tg_xxyz_yyyyyyzz_0[j] + fr * tg_xxyz_yyyyyyzz_1[j] + fl1_fx * (tg_xyz_yyyyyyzz_0[j] - tg_xyz_yyyyyyzz_1[j] * fl1_fza);

                    tg_xxxyz_yyyyyzzz_0[j] = pb_x * tg_xxyz_yyyyyzzz_0[j] + fr * tg_xxyz_yyyyyzzz_1[j] + fl1_fx * (tg_xyz_yyyyyzzz_0[j] - tg_xyz_yyyyyzzz_1[j] * fl1_fza);

                    tg_xxxyz_yyyyzzzz_0[j] = pb_x * tg_xxyz_yyyyzzzz_0[j] + fr * tg_xxyz_yyyyzzzz_1[j] + fl1_fx * (tg_xyz_yyyyzzzz_0[j] - tg_xyz_yyyyzzzz_1[j] * fl1_fza);

                    tg_xxxyz_yyyzzzzz_0[j] = pb_x * tg_xxyz_yyyzzzzz_0[j] + fr * tg_xxyz_yyyzzzzz_1[j] + fl1_fx * (tg_xyz_yyyzzzzz_0[j] - tg_xyz_yyyzzzzz_1[j] * fl1_fza);

                    tg_xxxyz_yyzzzzzz_0[j] = pb_x * tg_xxyz_yyzzzzzz_0[j] + fr * tg_xxyz_yyzzzzzz_1[j] + fl1_fx * (tg_xyz_yyzzzzzz_0[j] - tg_xyz_yyzzzzzz_1[j] * fl1_fza);

                    tg_xxxyz_yzzzzzzz_0[j] = pb_x * tg_xxyz_yzzzzzzz_0[j] + fr * tg_xxyz_yzzzzzzz_1[j] + fl1_fx * (tg_xyz_yzzzzzzz_0[j] - tg_xyz_yzzzzzzz_1[j] * fl1_fza);

                    tg_xxxyz_zzzzzzzz_0[j] = pb_x * tg_xxyz_zzzzzzzz_0[j] + fr * tg_xxyz_zzzzzzzz_1[j] + fl1_fx * (tg_xyz_zzzzzzzz_0[j] - tg_xyz_zzzzzzzz_1[j] * fl1_fza);

                    tg_xxxzz_xxxxxxxx_0[j] = pb_x * tg_xxzz_xxxxxxxx_0[j] + fr * tg_xxzz_xxxxxxxx_1[j] + fl1_fx * (tg_xzz_xxxxxxxx_0[j] - tg_xzz_xxxxxxxx_1[j] * fl1_fza) + 4.0 * fl1_fxn * tg_xxzz_xxxxxxx_1[j];

                    tg_xxxzz_xxxxxxxy_0[j] = pb_x * tg_xxzz_xxxxxxxy_0[j] + fr * tg_xxzz_xxxxxxxy_1[j] + fl1_fx * (tg_xzz_xxxxxxxy_0[j] - tg_xzz_xxxxxxxy_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_xxzz_xxxxxxy_1[j];

                    tg_xxxzz_xxxxxxxz_0[j] = pb_x * tg_xxzz_xxxxxxxz_0[j] + fr * tg_xxzz_xxxxxxxz_1[j] + fl1_fx * (tg_xzz_xxxxxxxz_0[j] - tg_xzz_xxxxxxxz_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_xxzz_xxxxxxz_1[j];

                    tg_xxxzz_xxxxxxyy_0[j] = pb_x * tg_xxzz_xxxxxxyy_0[j] + fr * tg_xxzz_xxxxxxyy_1[j] + fl1_fx * (tg_xzz_xxxxxxyy_0[j] - tg_xzz_xxxxxxyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxzz_xxxxxyy_1[j];

                    tg_xxxzz_xxxxxxyz_0[j] = pb_x * tg_xxzz_xxxxxxyz_0[j] + fr * tg_xxzz_xxxxxxyz_1[j] + fl1_fx * (tg_xzz_xxxxxxyz_0[j] - tg_xzz_xxxxxxyz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxzz_xxxxxyz_1[j];

                    tg_xxxzz_xxxxxxzz_0[j] = pb_x * tg_xxzz_xxxxxxzz_0[j] + fr * tg_xxzz_xxxxxxzz_1[j] + fl1_fx * (tg_xzz_xxxxxxzz_0[j] - tg_xzz_xxxxxxzz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxzz_xxxxxzz_1[j];

                    tg_xxxzz_xxxxxyyy_0[j] = pb_x * tg_xxzz_xxxxxyyy_0[j] + fr * tg_xxzz_xxxxxyyy_1[j] + fl1_fx * (tg_xzz_xxxxxyyy_0[j] - tg_xzz_xxxxxyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxzz_xxxxyyy_1[j];

                    tg_xxxzz_xxxxxyyz_0[j] = pb_x * tg_xxzz_xxxxxyyz_0[j] + fr * tg_xxzz_xxxxxyyz_1[j] + fl1_fx * (tg_xzz_xxxxxyyz_0[j] - tg_xzz_xxxxxyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxzz_xxxxyyz_1[j];

                    tg_xxxzz_xxxxxyzz_0[j] = pb_x * tg_xxzz_xxxxxyzz_0[j] + fr * tg_xxzz_xxxxxyzz_1[j] + fl1_fx * (tg_xzz_xxxxxyzz_0[j] - tg_xzz_xxxxxyzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxzz_xxxxyzz_1[j];

                    tg_xxxzz_xxxxxzzz_0[j] = pb_x * tg_xxzz_xxxxxzzz_0[j] + fr * tg_xxzz_xxxxxzzz_1[j] + fl1_fx * (tg_xzz_xxxxxzzz_0[j] - tg_xzz_xxxxxzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxzz_xxxxzzz_1[j];

                    tg_xxxzz_xxxxyyyy_0[j] = pb_x * tg_xxzz_xxxxyyyy_0[j] + fr * tg_xxzz_xxxxyyyy_1[j] + fl1_fx * (tg_xzz_xxxxyyyy_0[j] - tg_xzz_xxxxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxzz_xxxyyyy_1[j];

                    tg_xxxzz_xxxxyyyz_0[j] = pb_x * tg_xxzz_xxxxyyyz_0[j] + fr * tg_xxzz_xxxxyyyz_1[j] + fl1_fx * (tg_xzz_xxxxyyyz_0[j] - tg_xzz_xxxxyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxzz_xxxyyyz_1[j];

                    tg_xxxzz_xxxxyyzz_0[j] = pb_x * tg_xxzz_xxxxyyzz_0[j] + fr * tg_xxzz_xxxxyyzz_1[j] + fl1_fx * (tg_xzz_xxxxyyzz_0[j] - tg_xzz_xxxxyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxzz_xxxyyzz_1[j];

                    tg_xxxzz_xxxxyzzz_0[j] = pb_x * tg_xxzz_xxxxyzzz_0[j] + fr * tg_xxzz_xxxxyzzz_1[j] + fl1_fx * (tg_xzz_xxxxyzzz_0[j] - tg_xzz_xxxxyzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxzz_xxxyzzz_1[j];

                    tg_xxxzz_xxxxzzzz_0[j] = pb_x * tg_xxzz_xxxxzzzz_0[j] + fr * tg_xxzz_xxxxzzzz_1[j] + fl1_fx * (tg_xzz_xxxxzzzz_0[j] - tg_xzz_xxxxzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxzz_xxxzzzz_1[j];

                    tg_xxxzz_xxxyyyyy_0[j] = pb_x * tg_xxzz_xxxyyyyy_0[j] + fr * tg_xxzz_xxxyyyyy_1[j] + fl1_fx * (tg_xzz_xxxyyyyy_0[j] - tg_xzz_xxxyyyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzz_xxyyyyy_1[j];

                    tg_xxxzz_xxxyyyyz_0[j] = pb_x * tg_xxzz_xxxyyyyz_0[j] + fr * tg_xxzz_xxxyyyyz_1[j] + fl1_fx * (tg_xzz_xxxyyyyz_0[j] - tg_xzz_xxxyyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzz_xxyyyyz_1[j];

                    tg_xxxzz_xxxyyyzz_0[j] = pb_x * tg_xxzz_xxxyyyzz_0[j] + fr * tg_xxzz_xxxyyyzz_1[j] + fl1_fx * (tg_xzz_xxxyyyzz_0[j] - tg_xzz_xxxyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzz_xxyyyzz_1[j];

                    tg_xxxzz_xxxyyzzz_0[j] = pb_x * tg_xxzz_xxxyyzzz_0[j] + fr * tg_xxzz_xxxyyzzz_1[j] + fl1_fx * (tg_xzz_xxxyyzzz_0[j] - tg_xzz_xxxyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzz_xxyyzzz_1[j];

                    tg_xxxzz_xxxyzzzz_0[j] = pb_x * tg_xxzz_xxxyzzzz_0[j] + fr * tg_xxzz_xxxyzzzz_1[j] + fl1_fx * (tg_xzz_xxxyzzzz_0[j] - tg_xzz_xxxyzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzz_xxyzzzz_1[j];

                    tg_xxxzz_xxxzzzzz_0[j] = pb_x * tg_xxzz_xxxzzzzz_0[j] + fr * tg_xxzz_xxxzzzzz_1[j] + fl1_fx * (tg_xzz_xxxzzzzz_0[j] - tg_xzz_xxxzzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzz_xxzzzzz_1[j];

                    tg_xxxzz_xxyyyyyy_0[j] = pb_x * tg_xxzz_xxyyyyyy_0[j] + fr * tg_xxzz_xxyyyyyy_1[j] + fl1_fx * (tg_xzz_xxyyyyyy_0[j] - tg_xzz_xxyyyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxzz_xyyyyyy_1[j];

                    tg_xxxzz_xxyyyyyz_0[j] = pb_x * tg_xxzz_xxyyyyyz_0[j] + fr * tg_xxzz_xxyyyyyz_1[j] + fl1_fx * (tg_xzz_xxyyyyyz_0[j] - tg_xzz_xxyyyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxzz_xyyyyyz_1[j];

                    tg_xxxzz_xxyyyyzz_0[j] = pb_x * tg_xxzz_xxyyyyzz_0[j] + fr * tg_xxzz_xxyyyyzz_1[j] + fl1_fx * (tg_xzz_xxyyyyzz_0[j] - tg_xzz_xxyyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxzz_xyyyyzz_1[j];

                    tg_xxxzz_xxyyyzzz_0[j] = pb_x * tg_xxzz_xxyyyzzz_0[j] + fr * tg_xxzz_xxyyyzzz_1[j] + fl1_fx * (tg_xzz_xxyyyzzz_0[j] - tg_xzz_xxyyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxzz_xyyyzzz_1[j];

                    tg_xxxzz_xxyyzzzz_0[j] = pb_x * tg_xxzz_xxyyzzzz_0[j] + fr * tg_xxzz_xxyyzzzz_1[j] + fl1_fx * (tg_xzz_xxyyzzzz_0[j] - tg_xzz_xxyyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxzz_xyyzzzz_1[j];

                    tg_xxxzz_xxyzzzzz_0[j] = pb_x * tg_xxzz_xxyzzzzz_0[j] + fr * tg_xxzz_xxyzzzzz_1[j] + fl1_fx * (tg_xzz_xxyzzzzz_0[j] - tg_xzz_xxyzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxzz_xyzzzzz_1[j];

                    tg_xxxzz_xxzzzzzz_0[j] = pb_x * tg_xxzz_xxzzzzzz_0[j] + fr * tg_xxzz_xxzzzzzz_1[j] + fl1_fx * (tg_xzz_xxzzzzzz_0[j] - tg_xzz_xxzzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxzz_xzzzzzz_1[j];

                    tg_xxxzz_xyyyyyyy_0[j] = pb_x * tg_xxzz_xyyyyyyy_0[j] + fr * tg_xxzz_xyyyyyyy_1[j] + fl1_fx * (tg_xzz_xyyyyyyy_0[j] - tg_xzz_xyyyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzz_yyyyyyy_1[j];

                    tg_xxxzz_xyyyyyyz_0[j] = pb_x * tg_xxzz_xyyyyyyz_0[j] + fr * tg_xxzz_xyyyyyyz_1[j] + fl1_fx * (tg_xzz_xyyyyyyz_0[j] - tg_xzz_xyyyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzz_yyyyyyz_1[j];

                    tg_xxxzz_xyyyyyzz_0[j] = pb_x * tg_xxzz_xyyyyyzz_0[j] + fr * tg_xxzz_xyyyyyzz_1[j] + fl1_fx * (tg_xzz_xyyyyyzz_0[j] - tg_xzz_xyyyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzz_yyyyyzz_1[j];

                    tg_xxxzz_xyyyyzzz_0[j] = pb_x * tg_xxzz_xyyyyzzz_0[j] + fr * tg_xxzz_xyyyyzzz_1[j] + fl1_fx * (tg_xzz_xyyyyzzz_0[j] - tg_xzz_xyyyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzz_yyyyzzz_1[j];

                    tg_xxxzz_xyyyzzzz_0[j] = pb_x * tg_xxzz_xyyyzzzz_0[j] + fr * tg_xxzz_xyyyzzzz_1[j] + fl1_fx * (tg_xzz_xyyyzzzz_0[j] - tg_xzz_xyyyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzz_yyyzzzz_1[j];

                    tg_xxxzz_xyyzzzzz_0[j] = pb_x * tg_xxzz_xyyzzzzz_0[j] + fr * tg_xxzz_xyyzzzzz_1[j] + fl1_fx * (tg_xzz_xyyzzzzz_0[j] - tg_xzz_xyyzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzz_yyzzzzz_1[j];

                    tg_xxxzz_xyzzzzzz_0[j] = pb_x * tg_xxzz_xyzzzzzz_0[j] + fr * tg_xxzz_xyzzzzzz_1[j] + fl1_fx * (tg_xzz_xyzzzzzz_0[j] - tg_xzz_xyzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzz_yzzzzzz_1[j];

                    tg_xxxzz_xzzzzzzz_0[j] = pb_x * tg_xxzz_xzzzzzzz_0[j] + fr * tg_xxzz_xzzzzzzz_1[j] + fl1_fx * (tg_xzz_xzzzzzzz_0[j] - tg_xzz_xzzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzz_zzzzzzz_1[j];

                    tg_xxxzz_yyyyyyyy_0[j] = pb_x * tg_xxzz_yyyyyyyy_0[j] + fr * tg_xxzz_yyyyyyyy_1[j] + fl1_fx * (tg_xzz_yyyyyyyy_0[j] - tg_xzz_yyyyyyyy_1[j] * fl1_fza);

                    tg_xxxzz_yyyyyyyz_0[j] = pb_x * tg_xxzz_yyyyyyyz_0[j] + fr * tg_xxzz_yyyyyyyz_1[j] + fl1_fx * (tg_xzz_yyyyyyyz_0[j] - tg_xzz_yyyyyyyz_1[j] * fl1_fza);

                    tg_xxxzz_yyyyyyzz_0[j] = pb_x * tg_xxzz_yyyyyyzz_0[j] + fr * tg_xxzz_yyyyyyzz_1[j] + fl1_fx * (tg_xzz_yyyyyyzz_0[j] - tg_xzz_yyyyyyzz_1[j] * fl1_fza);

                    tg_xxxzz_yyyyyzzz_0[j] = pb_x * tg_xxzz_yyyyyzzz_0[j] + fr * tg_xxzz_yyyyyzzz_1[j] + fl1_fx * (tg_xzz_yyyyyzzz_0[j] - tg_xzz_yyyyyzzz_1[j] * fl1_fza);

                    tg_xxxzz_yyyyzzzz_0[j] = pb_x * tg_xxzz_yyyyzzzz_0[j] + fr * tg_xxzz_yyyyzzzz_1[j] + fl1_fx * (tg_xzz_yyyyzzzz_0[j] - tg_xzz_yyyyzzzz_1[j] * fl1_fza);

                    tg_xxxzz_yyyzzzzz_0[j] = pb_x * tg_xxzz_yyyzzzzz_0[j] + fr * tg_xxzz_yyyzzzzz_1[j] + fl1_fx * (tg_xzz_yyyzzzzz_0[j] - tg_xzz_yyyzzzzz_1[j] * fl1_fza);

                    tg_xxxzz_yyzzzzzz_0[j] = pb_x * tg_xxzz_yyzzzzzz_0[j] + fr * tg_xxzz_yyzzzzzz_1[j] + fl1_fx * (tg_xzz_yyzzzzzz_0[j] - tg_xzz_yyzzzzzz_1[j] * fl1_fza);

                    tg_xxxzz_yzzzzzzz_0[j] = pb_x * tg_xxzz_yzzzzzzz_0[j] + fr * tg_xxzz_yzzzzzzz_1[j] + fl1_fx * (tg_xzz_yzzzzzzz_0[j] - tg_xzz_yzzzzzzz_1[j] * fl1_fza);

                    tg_xxxzz_zzzzzzzz_0[j] = pb_x * tg_xxzz_zzzzzzzz_0[j] + fr * tg_xxzz_zzzzzzzz_1[j] + fl1_fx * (tg_xzz_zzzzzzzz_0[j] - tg_xzz_zzzzzzzz_1[j] * fl1_fza);

                    tg_xxyyy_xxxxxxxx_0[j] = pb_x * tg_xyyy_xxxxxxxx_0[j] + fr * tg_xyyy_xxxxxxxx_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxxxxxx_0[j] - tg_yyy_xxxxxxxx_1[j] * fl1_fza) + 4.0 * fl1_fxn * tg_xyyy_xxxxxxx_1[j];

                    tg_xxyyy_xxxxxxxy_0[j] = pb_x * tg_xyyy_xxxxxxxy_0[j] + fr * tg_xyyy_xxxxxxxy_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxxxxxy_0[j] - tg_yyy_xxxxxxxy_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_xyyy_xxxxxxy_1[j];

                    tg_xxyyy_xxxxxxxz_0[j] = pb_x * tg_xyyy_xxxxxxxz_0[j] + fr * tg_xyyy_xxxxxxxz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxxxxxz_0[j] - tg_yyy_xxxxxxxz_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_xyyy_xxxxxxz_1[j];

                    tg_xxyyy_xxxxxxyy_0[j] = pb_x * tg_xyyy_xxxxxxyy_0[j] + fr * tg_xyyy_xxxxxxyy_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxxxxyy_0[j] - tg_yyy_xxxxxxyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xyyy_xxxxxyy_1[j];

                    tg_xxyyy_xxxxxxyz_0[j] = pb_x * tg_xyyy_xxxxxxyz_0[j] + fr * tg_xyyy_xxxxxxyz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxxxxyz_0[j] - tg_yyy_xxxxxxyz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xyyy_xxxxxyz_1[j];

                    tg_xxyyy_xxxxxxzz_0[j] = pb_x * tg_xyyy_xxxxxxzz_0[j] + fr * tg_xyyy_xxxxxxzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxxxxzz_0[j] - tg_yyy_xxxxxxzz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xyyy_xxxxxzz_1[j];

                    tg_xxyyy_xxxxxyyy_0[j] = pb_x * tg_xyyy_xxxxxyyy_0[j] + fr * tg_xyyy_xxxxxyyy_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxxxyyy_0[j] - tg_yyy_xxxxxyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyy_xxxxyyy_1[j];

                    tg_xxyyy_xxxxxyyz_0[j] = pb_x * tg_xyyy_xxxxxyyz_0[j] + fr * tg_xyyy_xxxxxyyz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxxxyyz_0[j] - tg_yyy_xxxxxyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyy_xxxxyyz_1[j];

                    tg_xxyyy_xxxxxyzz_0[j] = pb_x * tg_xyyy_xxxxxyzz_0[j] + fr * tg_xyyy_xxxxxyzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxxxyzz_0[j] - tg_yyy_xxxxxyzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyy_xxxxyzz_1[j];

                    tg_xxyyy_xxxxxzzz_0[j] = pb_x * tg_xyyy_xxxxxzzz_0[j] + fr * tg_xyyy_xxxxxzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxxxzzz_0[j] - tg_yyy_xxxxxzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyy_xxxxzzz_1[j];

                    tg_xxyyy_xxxxyyyy_0[j] = pb_x * tg_xyyy_xxxxyyyy_0[j] + fr * tg_xyyy_xxxxyyyy_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxxyyyy_0[j] - tg_yyy_xxxxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyy_xxxyyyy_1[j];

                    tg_xxyyy_xxxxyyyz_0[j] = pb_x * tg_xyyy_xxxxyyyz_0[j] + fr * tg_xyyy_xxxxyyyz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxxyyyz_0[j] - tg_yyy_xxxxyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyy_xxxyyyz_1[j];

                    tg_xxyyy_xxxxyyzz_0[j] = pb_x * tg_xyyy_xxxxyyzz_0[j] + fr * tg_xyyy_xxxxyyzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxxyyzz_0[j] - tg_yyy_xxxxyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyy_xxxyyzz_1[j];

                    tg_xxyyy_xxxxyzzz_0[j] = pb_x * tg_xyyy_xxxxyzzz_0[j] + fr * tg_xyyy_xxxxyzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxxyzzz_0[j] - tg_yyy_xxxxyzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyy_xxxyzzz_1[j];

                    tg_xxyyy_xxxxzzzz_0[j] = pb_x * tg_xyyy_xxxxzzzz_0[j] + fr * tg_xyyy_xxxxzzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxxzzzz_0[j] - tg_yyy_xxxxzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyy_xxxzzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSL_285_380(      CMemBlock2D<double>* primBuffer,
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
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_8_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_8_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_8_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_8_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_8_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_7_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {7, -1, -1, -1}, 
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

                auto tg_xyyy_xxxyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 285); 

                auto tg_xyyy_xxxyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 286); 

                auto tg_xyyy_xxxyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 287); 

                auto tg_xyyy_xxxyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 288); 

                auto tg_xyyy_xxxyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 289); 

                auto tg_xyyy_xxxzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 290); 

                auto tg_xyyy_xxyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 291); 

                auto tg_xyyy_xxyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 292); 

                auto tg_xyyy_xxyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 293); 

                auto tg_xyyy_xxyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 294); 

                auto tg_xyyy_xxyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 295); 

                auto tg_xyyy_xxyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 296); 

                auto tg_xyyy_xxzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 297); 

                auto tg_xyyy_xyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 298); 

                auto tg_xyyy_xyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 299); 

                auto tg_xyyy_xyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 300); 

                auto tg_xyyy_xyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 301); 

                auto tg_xyyy_xyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 302); 

                auto tg_xyyy_xyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 303); 

                auto tg_xyyy_xyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 304); 

                auto tg_xyyy_xzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 305); 

                auto tg_xyyy_yyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 306); 

                auto tg_xyyy_yyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 307); 

                auto tg_xyyy_yyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 308); 

                auto tg_xyyy_yyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 309); 

                auto tg_xyyy_yyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 310); 

                auto tg_xyyy_yyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 311); 

                auto tg_xyyy_yyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 312); 

                auto tg_xyyy_yzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 313); 

                auto tg_xyyy_zzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 314); 

                auto tg_xyyz_xxxxxxxx_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 315); 

                auto tg_xyyz_xxxxxxxy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 316); 

                auto tg_xyyz_xxxxxxxz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 317); 

                auto tg_xyyz_xxxxxxyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 318); 

                auto tg_xyyz_xxxxxxyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 319); 

                auto tg_xyyz_xxxxxxzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 320); 

                auto tg_xyyz_xxxxxyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 321); 

                auto tg_xyyz_xxxxxyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 322); 

                auto tg_xyyz_xxxxxyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 323); 

                auto tg_xyyz_xxxxxzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 324); 

                auto tg_xyyz_xxxxyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 325); 

                auto tg_xyyz_xxxxyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 326); 

                auto tg_xyyz_xxxxyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 327); 

                auto tg_xyyz_xxxxyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 328); 

                auto tg_xyyz_xxxxzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 329); 

                auto tg_xyyz_xxxyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 330); 

                auto tg_xyyz_xxxyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 331); 

                auto tg_xyyz_xxxyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 332); 

                auto tg_xyyz_xxxyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 333); 

                auto tg_xyyz_xxxyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 334); 

                auto tg_xyyz_xxxzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 335); 

                auto tg_xyyz_xxyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 336); 

                auto tg_xyyz_xxyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 337); 

                auto tg_xyyz_xxyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 338); 

                auto tg_xyyz_xxyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 339); 

                auto tg_xyyz_xxyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 340); 

                auto tg_xyyz_xxyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 341); 

                auto tg_xyyz_xxzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 342); 

                auto tg_xyyz_xyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 343); 

                auto tg_xyyz_xyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 344); 

                auto tg_xyyz_xyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 345); 

                auto tg_xyyz_xyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 346); 

                auto tg_xyyz_xyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 347); 

                auto tg_xyyz_xyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 348); 

                auto tg_xyyz_xyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 349); 

                auto tg_xyyz_xzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 350); 

                auto tg_xyyz_yyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 351); 

                auto tg_xyyz_yyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 352); 

                auto tg_xyyz_yyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 353); 

                auto tg_xyyz_yyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 354); 

                auto tg_xyyz_yyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 355); 

                auto tg_xyyz_yyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 356); 

                auto tg_xyyz_yyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 357); 

                auto tg_xyyz_yzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 358); 

                auto tg_xyyz_zzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 359); 

                auto tg_xyzz_xxxxxxxx_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 360); 

                auto tg_xyzz_xxxxxxxy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 361); 

                auto tg_xyzz_xxxxxxxz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 362); 

                auto tg_xyzz_xxxxxxyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 363); 

                auto tg_xyzz_xxxxxxyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 364); 

                auto tg_xyzz_xxxxxxzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 365); 

                auto tg_xyzz_xxxxxyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 366); 

                auto tg_xyzz_xxxxxyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 367); 

                auto tg_xyzz_xxxxxyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 368); 

                auto tg_xyzz_xxxxxzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 369); 

                auto tg_xyzz_xxxxyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 370); 

                auto tg_xyzz_xxxxyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 371); 

                auto tg_xyzz_xxxxyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 372); 

                auto tg_xyzz_xxxxyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 373); 

                auto tg_xyzz_xxxxzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 374); 

                auto tg_xyzz_xxxyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 375); 

                auto tg_xyzz_xxxyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 376); 

                auto tg_xyzz_xxxyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 377); 

                auto tg_xyzz_xxxyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 378); 

                auto tg_xyzz_xxxyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 379); 

                auto tg_xyyy_xxxyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 285); 

                auto tg_xyyy_xxxyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 286); 

                auto tg_xyyy_xxxyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 287); 

                auto tg_xyyy_xxxyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 288); 

                auto tg_xyyy_xxxyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 289); 

                auto tg_xyyy_xxxzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 290); 

                auto tg_xyyy_xxyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 291); 

                auto tg_xyyy_xxyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 292); 

                auto tg_xyyy_xxyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 293); 

                auto tg_xyyy_xxyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 294); 

                auto tg_xyyy_xxyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 295); 

                auto tg_xyyy_xxyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 296); 

                auto tg_xyyy_xxzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 297); 

                auto tg_xyyy_xyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 298); 

                auto tg_xyyy_xyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 299); 

                auto tg_xyyy_xyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 300); 

                auto tg_xyyy_xyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 301); 

                auto tg_xyyy_xyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 302); 

                auto tg_xyyy_xyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 303); 

                auto tg_xyyy_xyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 304); 

                auto tg_xyyy_xzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 305); 

                auto tg_xyyy_yyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 306); 

                auto tg_xyyy_yyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 307); 

                auto tg_xyyy_yyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 308); 

                auto tg_xyyy_yyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 309); 

                auto tg_xyyy_yyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 310); 

                auto tg_xyyy_yyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 311); 

                auto tg_xyyy_yyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 312); 

                auto tg_xyyy_yzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 313); 

                auto tg_xyyy_zzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 314); 

                auto tg_xyyz_xxxxxxxx_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 315); 

                auto tg_xyyz_xxxxxxxy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 316); 

                auto tg_xyyz_xxxxxxxz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 317); 

                auto tg_xyyz_xxxxxxyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 318); 

                auto tg_xyyz_xxxxxxyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 319); 

                auto tg_xyyz_xxxxxxzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 320); 

                auto tg_xyyz_xxxxxyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 321); 

                auto tg_xyyz_xxxxxyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 322); 

                auto tg_xyyz_xxxxxyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 323); 

                auto tg_xyyz_xxxxxzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 324); 

                auto tg_xyyz_xxxxyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 325); 

                auto tg_xyyz_xxxxyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 326); 

                auto tg_xyyz_xxxxyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 327); 

                auto tg_xyyz_xxxxyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 328); 

                auto tg_xyyz_xxxxzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 329); 

                auto tg_xyyz_xxxyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 330); 

                auto tg_xyyz_xxxyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 331); 

                auto tg_xyyz_xxxyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 332); 

                auto tg_xyyz_xxxyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 333); 

                auto tg_xyyz_xxxyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 334); 

                auto tg_xyyz_xxxzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 335); 

                auto tg_xyyz_xxyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 336); 

                auto tg_xyyz_xxyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 337); 

                auto tg_xyyz_xxyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 338); 

                auto tg_xyyz_xxyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 339); 

                auto tg_xyyz_xxyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 340); 

                auto tg_xyyz_xxyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 341); 

                auto tg_xyyz_xxzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 342); 

                auto tg_xyyz_xyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 343); 

                auto tg_xyyz_xyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 344); 

                auto tg_xyyz_xyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 345); 

                auto tg_xyyz_xyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 346); 

                auto tg_xyyz_xyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 347); 

                auto tg_xyyz_xyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 348); 

                auto tg_xyyz_xyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 349); 

                auto tg_xyyz_xzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 350); 

                auto tg_xyyz_yyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 351); 

                auto tg_xyyz_yyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 352); 

                auto tg_xyyz_yyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 353); 

                auto tg_xyyz_yyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 354); 

                auto tg_xyyz_yyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 355); 

                auto tg_xyyz_yyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 356); 

                auto tg_xyyz_yyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 357); 

                auto tg_xyyz_yzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 358); 

                auto tg_xyyz_zzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 359); 

                auto tg_xyzz_xxxxxxxx_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 360); 

                auto tg_xyzz_xxxxxxxy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 361); 

                auto tg_xyzz_xxxxxxxz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 362); 

                auto tg_xyzz_xxxxxxyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 363); 

                auto tg_xyzz_xxxxxxyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 364); 

                auto tg_xyzz_xxxxxxzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 365); 

                auto tg_xyzz_xxxxxyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 366); 

                auto tg_xyzz_xxxxxyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 367); 

                auto tg_xyzz_xxxxxyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 368); 

                auto tg_xyzz_xxxxxzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 369); 

                auto tg_xyzz_xxxxyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 370); 

                auto tg_xyzz_xxxxyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 371); 

                auto tg_xyzz_xxxxyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 372); 

                auto tg_xyzz_xxxxyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 373); 

                auto tg_xyzz_xxxxzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 374); 

                auto tg_xyzz_xxxyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 375); 

                auto tg_xyzz_xxxyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 376); 

                auto tg_xyzz_xxxyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 377); 

                auto tg_xyzz_xxxyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 378); 

                auto tg_xyzz_xxxyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 379); 

                auto tg_yyy_xxxyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 285); 

                auto tg_yyy_xxxyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 286); 

                auto tg_yyy_xxxyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 287); 

                auto tg_yyy_xxxyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 288); 

                auto tg_yyy_xxxyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 289); 

                auto tg_yyy_xxxzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 290); 

                auto tg_yyy_xxyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 291); 

                auto tg_yyy_xxyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 292); 

                auto tg_yyy_xxyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 293); 

                auto tg_yyy_xxyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 294); 

                auto tg_yyy_xxyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 295); 

                auto tg_yyy_xxyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 296); 

                auto tg_yyy_xxzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 297); 

                auto tg_yyy_xyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 298); 

                auto tg_yyy_xyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 299); 

                auto tg_yyy_xyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 300); 

                auto tg_yyy_xyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 301); 

                auto tg_yyy_xyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 302); 

                auto tg_yyy_xyyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 303); 

                auto tg_yyy_xyzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 304); 

                auto tg_yyy_xzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 305); 

                auto tg_yyy_yyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 306); 

                auto tg_yyy_yyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 307); 

                auto tg_yyy_yyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 308); 

                auto tg_yyy_yyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 309); 

                auto tg_yyy_yyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 310); 

                auto tg_yyy_yyyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 311); 

                auto tg_yyy_yyzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 312); 

                auto tg_yyy_yzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 313); 

                auto tg_yyy_zzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 314); 

                auto tg_yyz_xxxxxxxx_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 315); 

                auto tg_yyz_xxxxxxxy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 316); 

                auto tg_yyz_xxxxxxxz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 317); 

                auto tg_yyz_xxxxxxyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 318); 

                auto tg_yyz_xxxxxxyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 319); 

                auto tg_yyz_xxxxxxzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 320); 

                auto tg_yyz_xxxxxyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 321); 

                auto tg_yyz_xxxxxyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 322); 

                auto tg_yyz_xxxxxyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 323); 

                auto tg_yyz_xxxxxzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 324); 

                auto tg_yyz_xxxxyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 325); 

                auto tg_yyz_xxxxyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 326); 

                auto tg_yyz_xxxxyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 327); 

                auto tg_yyz_xxxxyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 328); 

                auto tg_yyz_xxxxzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 329); 

                auto tg_yyz_xxxyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 330); 

                auto tg_yyz_xxxyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 331); 

                auto tg_yyz_xxxyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 332); 

                auto tg_yyz_xxxyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 333); 

                auto tg_yyz_xxxyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 334); 

                auto tg_yyz_xxxzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 335); 

                auto tg_yyz_xxyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 336); 

                auto tg_yyz_xxyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 337); 

                auto tg_yyz_xxyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 338); 

                auto tg_yyz_xxyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 339); 

                auto tg_yyz_xxyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 340); 

                auto tg_yyz_xxyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 341); 

                auto tg_yyz_xxzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 342); 

                auto tg_yyz_xyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 343); 

                auto tg_yyz_xyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 344); 

                auto tg_yyz_xyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 345); 

                auto tg_yyz_xyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 346); 

                auto tg_yyz_xyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 347); 

                auto tg_yyz_xyyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 348); 

                auto tg_yyz_xyzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 349); 

                auto tg_yyz_xzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 350); 

                auto tg_yyz_yyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 351); 

                auto tg_yyz_yyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 352); 

                auto tg_yyz_yyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 353); 

                auto tg_yyz_yyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 354); 

                auto tg_yyz_yyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 355); 

                auto tg_yyz_yyyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 356); 

                auto tg_yyz_yyzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 357); 

                auto tg_yyz_yzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 358); 

                auto tg_yyz_zzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 359); 

                auto tg_yzz_xxxxxxxx_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 360); 

                auto tg_yzz_xxxxxxxy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 361); 

                auto tg_yzz_xxxxxxxz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 362); 

                auto tg_yzz_xxxxxxyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 363); 

                auto tg_yzz_xxxxxxyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 364); 

                auto tg_yzz_xxxxxxzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 365); 

                auto tg_yzz_xxxxxyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 366); 

                auto tg_yzz_xxxxxyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 367); 

                auto tg_yzz_xxxxxyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 368); 

                auto tg_yzz_xxxxxzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 369); 

                auto tg_yzz_xxxxyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 370); 

                auto tg_yzz_xxxxyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 371); 

                auto tg_yzz_xxxxyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 372); 

                auto tg_yzz_xxxxyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 373); 

                auto tg_yzz_xxxxzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 374); 

                auto tg_yzz_xxxyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 375); 

                auto tg_yzz_xxxyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 376); 

                auto tg_yzz_xxxyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 377); 

                auto tg_yzz_xxxyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 378); 

                auto tg_yzz_xxxyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 379); 

                auto tg_yyy_xxxyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 285); 

                auto tg_yyy_xxxyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 286); 

                auto tg_yyy_xxxyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 287); 

                auto tg_yyy_xxxyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 288); 

                auto tg_yyy_xxxyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 289); 

                auto tg_yyy_xxxzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 290); 

                auto tg_yyy_xxyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 291); 

                auto tg_yyy_xxyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 292); 

                auto tg_yyy_xxyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 293); 

                auto tg_yyy_xxyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 294); 

                auto tg_yyy_xxyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 295); 

                auto tg_yyy_xxyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 296); 

                auto tg_yyy_xxzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 297); 

                auto tg_yyy_xyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 298); 

                auto tg_yyy_xyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 299); 

                auto tg_yyy_xyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 300); 

                auto tg_yyy_xyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 301); 

                auto tg_yyy_xyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 302); 

                auto tg_yyy_xyyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 303); 

                auto tg_yyy_xyzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 304); 

                auto tg_yyy_xzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 305); 

                auto tg_yyy_yyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 306); 

                auto tg_yyy_yyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 307); 

                auto tg_yyy_yyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 308); 

                auto tg_yyy_yyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 309); 

                auto tg_yyy_yyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 310); 

                auto tg_yyy_yyyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 311); 

                auto tg_yyy_yyzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 312); 

                auto tg_yyy_yzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 313); 

                auto tg_yyy_zzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 314); 

                auto tg_yyz_xxxxxxxx_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 315); 

                auto tg_yyz_xxxxxxxy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 316); 

                auto tg_yyz_xxxxxxxz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 317); 

                auto tg_yyz_xxxxxxyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 318); 

                auto tg_yyz_xxxxxxyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 319); 

                auto tg_yyz_xxxxxxzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 320); 

                auto tg_yyz_xxxxxyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 321); 

                auto tg_yyz_xxxxxyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 322); 

                auto tg_yyz_xxxxxyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 323); 

                auto tg_yyz_xxxxxzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 324); 

                auto tg_yyz_xxxxyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 325); 

                auto tg_yyz_xxxxyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 326); 

                auto tg_yyz_xxxxyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 327); 

                auto tg_yyz_xxxxyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 328); 

                auto tg_yyz_xxxxzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 329); 

                auto tg_yyz_xxxyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 330); 

                auto tg_yyz_xxxyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 331); 

                auto tg_yyz_xxxyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 332); 

                auto tg_yyz_xxxyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 333); 

                auto tg_yyz_xxxyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 334); 

                auto tg_yyz_xxxzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 335); 

                auto tg_yyz_xxyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 336); 

                auto tg_yyz_xxyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 337); 

                auto tg_yyz_xxyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 338); 

                auto tg_yyz_xxyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 339); 

                auto tg_yyz_xxyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 340); 

                auto tg_yyz_xxyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 341); 

                auto tg_yyz_xxzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 342); 

                auto tg_yyz_xyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 343); 

                auto tg_yyz_xyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 344); 

                auto tg_yyz_xyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 345); 

                auto tg_yyz_xyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 346); 

                auto tg_yyz_xyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 347); 

                auto tg_yyz_xyyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 348); 

                auto tg_yyz_xyzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 349); 

                auto tg_yyz_xzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 350); 

                auto tg_yyz_yyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 351); 

                auto tg_yyz_yyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 352); 

                auto tg_yyz_yyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 353); 

                auto tg_yyz_yyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 354); 

                auto tg_yyz_yyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 355); 

                auto tg_yyz_yyyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 356); 

                auto tg_yyz_yyzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 357); 

                auto tg_yyz_yzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 358); 

                auto tg_yyz_zzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 359); 

                auto tg_yzz_xxxxxxxx_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 360); 

                auto tg_yzz_xxxxxxxy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 361); 

                auto tg_yzz_xxxxxxxz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 362); 

                auto tg_yzz_xxxxxxyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 363); 

                auto tg_yzz_xxxxxxyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 364); 

                auto tg_yzz_xxxxxxzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 365); 

                auto tg_yzz_xxxxxyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 366); 

                auto tg_yzz_xxxxxyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 367); 

                auto tg_yzz_xxxxxyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 368); 

                auto tg_yzz_xxxxxzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 369); 

                auto tg_yzz_xxxxyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 370); 

                auto tg_yzz_xxxxyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 371); 

                auto tg_yzz_xxxxyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 372); 

                auto tg_yzz_xxxxyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 373); 

                auto tg_yzz_xxxxzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 374); 

                auto tg_yzz_xxxyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 375); 

                auto tg_yzz_xxxyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 376); 

                auto tg_yzz_xxxyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 377); 

                auto tg_yzz_xxxyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 378); 

                auto tg_yzz_xxxyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 379); 

                auto tg_xyyy_xxyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 231); 

                auto tg_xyyy_xxyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 232); 

                auto tg_xyyy_xxyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 233); 

                auto tg_xyyy_xxyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 234); 

                auto tg_xyyy_xxyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 235); 

                auto tg_xyyy_xxzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 236); 

                auto tg_xyyy_xyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 237); 

                auto tg_xyyy_xyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 238); 

                auto tg_xyyy_xyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 239); 

                auto tg_xyyy_xyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 240); 

                auto tg_xyyy_xyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 241); 

                auto tg_xyyy_xyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 242); 

                auto tg_xyyy_xzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 243); 

                auto tg_xyyy_yyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 244); 

                auto tg_xyyy_yyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 245); 

                auto tg_xyyy_yyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 246); 

                auto tg_xyyy_yyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 247); 

                auto tg_xyyy_yyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 248); 

                auto tg_xyyy_yyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 249); 

                auto tg_xyyy_yzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 250); 

                auto tg_xyyy_zzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 251); 

                auto tg_xyyz_xxxxxxx_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 252); 

                auto tg_xyyz_xxxxxxy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 253); 

                auto tg_xyyz_xxxxxxz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 254); 

                auto tg_xyyz_xxxxxyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 255); 

                auto tg_xyyz_xxxxxyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 256); 

                auto tg_xyyz_xxxxxzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 257); 

                auto tg_xyyz_xxxxyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 258); 

                auto tg_xyyz_xxxxyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 259); 

                auto tg_xyyz_xxxxyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 260); 

                auto tg_xyyz_xxxxzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 261); 

                auto tg_xyyz_xxxyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 262); 

                auto tg_xyyz_xxxyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 263); 

                auto tg_xyyz_xxxyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 264); 

                auto tg_xyyz_xxxyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 265); 

                auto tg_xyyz_xxxzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 266); 

                auto tg_xyyz_xxyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 267); 

                auto tg_xyyz_xxyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 268); 

                auto tg_xyyz_xxyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 269); 

                auto tg_xyyz_xxyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 270); 

                auto tg_xyyz_xxyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 271); 

                auto tg_xyyz_xxzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 272); 

                auto tg_xyyz_xyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 273); 

                auto tg_xyyz_xyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 274); 

                auto tg_xyyz_xyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 275); 

                auto tg_xyyz_xyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 276); 

                auto tg_xyyz_xyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 277); 

                auto tg_xyyz_xyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 278); 

                auto tg_xyyz_xzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 279); 

                auto tg_xyyz_yyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 280); 

                auto tg_xyyz_yyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 281); 

                auto tg_xyyz_yyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 282); 

                auto tg_xyyz_yyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 283); 

                auto tg_xyyz_yyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 284); 

                auto tg_xyyz_yyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 285); 

                auto tg_xyyz_yzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 286); 

                auto tg_xyyz_zzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 287); 

                auto tg_xyzz_xxxxxxx_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 288); 

                auto tg_xyzz_xxxxxxy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 289); 

                auto tg_xyzz_xxxxxxz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 290); 

                auto tg_xyzz_xxxxxyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 291); 

                auto tg_xyzz_xxxxxyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 292); 

                auto tg_xyzz_xxxxxzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 293); 

                auto tg_xyzz_xxxxyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 294); 

                auto tg_xyzz_xxxxyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 295); 

                auto tg_xyzz_xxxxyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 296); 

                auto tg_xyzz_xxxxzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 297); 

                auto tg_xyzz_xxxyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 298); 

                auto tg_xyzz_xxxyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 299); 

                auto tg_xyzz_xxxyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 300); 

                auto tg_xyzz_xxxyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 301); 

                auto tg_xyzz_xxxzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 302); 

                auto tg_xyzz_xxyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 303); 

                auto tg_xyzz_xxyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 304); 

                auto tg_xyzz_xxyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 305); 

                auto tg_xyzz_xxyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 306); 

                auto tg_xyzz_xxyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 307); 

                // set up pointers to integrals

                auto tg_xxyyy_xxxyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 285); 

                auto tg_xxyyy_xxxyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 286); 

                auto tg_xxyyy_xxxyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 287); 

                auto tg_xxyyy_xxxyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 288); 

                auto tg_xxyyy_xxxyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 289); 

                auto tg_xxyyy_xxxzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 290); 

                auto tg_xxyyy_xxyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 291); 

                auto tg_xxyyy_xxyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 292); 

                auto tg_xxyyy_xxyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 293); 

                auto tg_xxyyy_xxyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 294); 

                auto tg_xxyyy_xxyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 295); 

                auto tg_xxyyy_xxyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 296); 

                auto tg_xxyyy_xxzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 297); 

                auto tg_xxyyy_xyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 298); 

                auto tg_xxyyy_xyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 299); 

                auto tg_xxyyy_xyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 300); 

                auto tg_xxyyy_xyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 301); 

                auto tg_xxyyy_xyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 302); 

                auto tg_xxyyy_xyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 303); 

                auto tg_xxyyy_xyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 304); 

                auto tg_xxyyy_xzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 305); 

                auto tg_xxyyy_yyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 306); 

                auto tg_xxyyy_yyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 307); 

                auto tg_xxyyy_yyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 308); 

                auto tg_xxyyy_yyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 309); 

                auto tg_xxyyy_yyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 310); 

                auto tg_xxyyy_yyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 311); 

                auto tg_xxyyy_yyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 312); 

                auto tg_xxyyy_yzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 313); 

                auto tg_xxyyy_zzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 314); 

                auto tg_xxyyz_xxxxxxxx_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 315); 

                auto tg_xxyyz_xxxxxxxy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 316); 

                auto tg_xxyyz_xxxxxxxz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 317); 

                auto tg_xxyyz_xxxxxxyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 318); 

                auto tg_xxyyz_xxxxxxyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 319); 

                auto tg_xxyyz_xxxxxxzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 320); 

                auto tg_xxyyz_xxxxxyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 321); 

                auto tg_xxyyz_xxxxxyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 322); 

                auto tg_xxyyz_xxxxxyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 323); 

                auto tg_xxyyz_xxxxxzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 324); 

                auto tg_xxyyz_xxxxyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 325); 

                auto tg_xxyyz_xxxxyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 326); 

                auto tg_xxyyz_xxxxyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 327); 

                auto tg_xxyyz_xxxxyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 328); 

                auto tg_xxyyz_xxxxzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 329); 

                auto tg_xxyyz_xxxyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 330); 

                auto tg_xxyyz_xxxyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 331); 

                auto tg_xxyyz_xxxyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 332); 

                auto tg_xxyyz_xxxyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 333); 

                auto tg_xxyyz_xxxyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 334); 

                auto tg_xxyyz_xxxzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 335); 

                auto tg_xxyyz_xxyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 336); 

                auto tg_xxyyz_xxyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 337); 

                auto tg_xxyyz_xxyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 338); 

                auto tg_xxyyz_xxyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 339); 

                auto tg_xxyyz_xxyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 340); 

                auto tg_xxyyz_xxyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 341); 

                auto tg_xxyyz_xxzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 342); 

                auto tg_xxyyz_xyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 343); 

                auto tg_xxyyz_xyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 344); 

                auto tg_xxyyz_xyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 345); 

                auto tg_xxyyz_xyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 346); 

                auto tg_xxyyz_xyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 347); 

                auto tg_xxyyz_xyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 348); 

                auto tg_xxyyz_xyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 349); 

                auto tg_xxyyz_xzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 350); 

                auto tg_xxyyz_yyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 351); 

                auto tg_xxyyz_yyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 352); 

                auto tg_xxyyz_yyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 353); 

                auto tg_xxyyz_yyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 354); 

                auto tg_xxyyz_yyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 355); 

                auto tg_xxyyz_yyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 356); 

                auto tg_xxyyz_yyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 357); 

                auto tg_xxyyz_yzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 358); 

                auto tg_xxyyz_zzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 359); 

                auto tg_xxyzz_xxxxxxxx_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 360); 

                auto tg_xxyzz_xxxxxxxy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 361); 

                auto tg_xxyzz_xxxxxxxz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 362); 

                auto tg_xxyzz_xxxxxxyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 363); 

                auto tg_xxyzz_xxxxxxyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 364); 

                auto tg_xxyzz_xxxxxxzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 365); 

                auto tg_xxyzz_xxxxxyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 366); 

                auto tg_xxyzz_xxxxxyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 367); 

                auto tg_xxyzz_xxxxxyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 368); 

                auto tg_xxyzz_xxxxxzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 369); 

                auto tg_xxyzz_xxxxyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 370); 

                auto tg_xxyzz_xxxxyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 371); 

                auto tg_xxyzz_xxxxyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 372); 

                auto tg_xxyzz_xxxxyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 373); 

                auto tg_xxyzz_xxxxzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 374); 

                auto tg_xxyzz_xxxyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 375); 

                auto tg_xxyzz_xxxyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 376); 

                auto tg_xxyzz_xxxyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 377); 

                auto tg_xxyzz_xxxyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 378); 

                auto tg_xxyzz_xxxyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 379); 

                // Batch of Integrals (285,380)

                #pragma omp simd aligned(fxn, fza, tg_xxyyy_xxxyyyyy_0, tg_xxyyy_xxxyyyyz_0, \
                                         tg_xxyyy_xxxyyyzz_0, tg_xxyyy_xxxyyzzz_0, tg_xxyyy_xxxyzzzz_0, tg_xxyyy_xxxzzzzz_0, \
                                         tg_xxyyy_xxyyyyyy_0, tg_xxyyy_xxyyyyyz_0, tg_xxyyy_xxyyyyzz_0, tg_xxyyy_xxyyyzzz_0, \
                                         tg_xxyyy_xxyyzzzz_0, tg_xxyyy_xxyzzzzz_0, tg_xxyyy_xxzzzzzz_0, tg_xxyyy_xyyyyyyy_0, \
                                         tg_xxyyy_xyyyyyyz_0, tg_xxyyy_xyyyyyzz_0, tg_xxyyy_xyyyyzzz_0, tg_xxyyy_xyyyzzzz_0, \
                                         tg_xxyyy_xyyzzzzz_0, tg_xxyyy_xyzzzzzz_0, tg_xxyyy_xzzzzzzz_0, tg_xxyyy_yyyyyyyy_0, \
                                         tg_xxyyy_yyyyyyyz_0, tg_xxyyy_yyyyyyzz_0, tg_xxyyy_yyyyyzzz_0, tg_xxyyy_yyyyzzzz_0, \
                                         tg_xxyyy_yyyzzzzz_0, tg_xxyyy_yyzzzzzz_0, tg_xxyyy_yzzzzzzz_0, tg_xxyyy_zzzzzzzz_0, \
                                         tg_xxyyz_xxxxxxxx_0, tg_xxyyz_xxxxxxxy_0, tg_xxyyz_xxxxxxxz_0, tg_xxyyz_xxxxxxyy_0, \
                                         tg_xxyyz_xxxxxxyz_0, tg_xxyyz_xxxxxxzz_0, tg_xxyyz_xxxxxyyy_0, tg_xxyyz_xxxxxyyz_0, \
                                         tg_xxyyz_xxxxxyzz_0, tg_xxyyz_xxxxxzzz_0, tg_xxyyz_xxxxyyyy_0, tg_xxyyz_xxxxyyyz_0, \
                                         tg_xxyyz_xxxxyyzz_0, tg_xxyyz_xxxxyzzz_0, tg_xxyyz_xxxxzzzz_0, tg_xxyyz_xxxyyyyy_0, \
                                         tg_xxyyz_xxxyyyyz_0, tg_xxyyz_xxxyyyzz_0, tg_xxyyz_xxxyyzzz_0, tg_xxyyz_xxxyzzzz_0, \
                                         tg_xxyyz_xxxzzzzz_0, tg_xxyyz_xxyyyyyy_0, tg_xxyyz_xxyyyyyz_0, tg_xxyyz_xxyyyyzz_0, \
                                         tg_xxyyz_xxyyyzzz_0, tg_xxyyz_xxyyzzzz_0, tg_xxyyz_xxyzzzzz_0, tg_xxyyz_xxzzzzzz_0, \
                                         tg_xxyyz_xyyyyyyy_0, tg_xxyyz_xyyyyyyz_0, tg_xxyyz_xyyyyyzz_0, tg_xxyyz_xyyyyzzz_0, \
                                         tg_xxyyz_xyyyzzzz_0, tg_xxyyz_xyyzzzzz_0, tg_xxyyz_xyzzzzzz_0, tg_xxyyz_xzzzzzzz_0, \
                                         tg_xxyyz_yyyyyyyy_0, tg_xxyyz_yyyyyyyz_0, tg_xxyyz_yyyyyyzz_0, tg_xxyyz_yyyyyzzz_0, \
                                         tg_xxyyz_yyyyzzzz_0, tg_xxyyz_yyyzzzzz_0, tg_xxyyz_yyzzzzzz_0, tg_xxyyz_yzzzzzzz_0, \
                                         tg_xxyyz_zzzzzzzz_0, tg_xxyzz_xxxxxxxx_0, tg_xxyzz_xxxxxxxy_0, tg_xxyzz_xxxxxxxz_0, \
                                         tg_xxyzz_xxxxxxyy_0, tg_xxyzz_xxxxxxyz_0, tg_xxyzz_xxxxxxzz_0, tg_xxyzz_xxxxxyyy_0, \
                                         tg_xxyzz_xxxxxyyz_0, tg_xxyzz_xxxxxyzz_0, tg_xxyzz_xxxxxzzz_0, tg_xxyzz_xxxxyyyy_0, \
                                         tg_xxyzz_xxxxyyyz_0, tg_xxyzz_xxxxyyzz_0, tg_xxyzz_xxxxyzzz_0, tg_xxyzz_xxxxzzzz_0, \
                                         tg_xxyzz_xxxyyyyy_0, tg_xxyzz_xxxyyyyz_0, tg_xxyzz_xxxyyyzz_0, tg_xxyzz_xxxyyzzz_0, \
                                         tg_xxyzz_xxxyzzzz_0, tg_xyyy_xxxyyyyy_0, tg_xyyy_xxxyyyyy_1, tg_xyyy_xxxyyyyz_0, \
                                         tg_xyyy_xxxyyyyz_1, tg_xyyy_xxxyyyzz_0, tg_xyyy_xxxyyyzz_1, tg_xyyy_xxxyyzzz_0, \
                                         tg_xyyy_xxxyyzzz_1, tg_xyyy_xxxyzzzz_0, tg_xyyy_xxxyzzzz_1, tg_xyyy_xxxzzzzz_0, \
                                         tg_xyyy_xxxzzzzz_1, tg_xyyy_xxyyyyy_1, tg_xyyy_xxyyyyyy_0, tg_xyyy_xxyyyyyy_1, \
                                         tg_xyyy_xxyyyyyz_0, tg_xyyy_xxyyyyyz_1, tg_xyyy_xxyyyyz_1, tg_xyyy_xxyyyyzz_0, \
                                         tg_xyyy_xxyyyyzz_1, tg_xyyy_xxyyyzz_1, tg_xyyy_xxyyyzzz_0, tg_xyyy_xxyyyzzz_1, \
                                         tg_xyyy_xxyyzzz_1, tg_xyyy_xxyyzzzz_0, tg_xyyy_xxyyzzzz_1, tg_xyyy_xxyzzzz_1, \
                                         tg_xyyy_xxyzzzzz_0, tg_xyyy_xxyzzzzz_1, tg_xyyy_xxzzzzz_1, tg_xyyy_xxzzzzzz_0, \
                                         tg_xyyy_xxzzzzzz_1, tg_xyyy_xyyyyyy_1, tg_xyyy_xyyyyyyy_0, tg_xyyy_xyyyyyyy_1, \
                                         tg_xyyy_xyyyyyyz_0, tg_xyyy_xyyyyyyz_1, tg_xyyy_xyyyyyz_1, tg_xyyy_xyyyyyzz_0, \
                                         tg_xyyy_xyyyyyzz_1, tg_xyyy_xyyyyzz_1, tg_xyyy_xyyyyzzz_0, tg_xyyy_xyyyyzzz_1, \
                                         tg_xyyy_xyyyzzz_1, tg_xyyy_xyyyzzzz_0, tg_xyyy_xyyyzzzz_1, tg_xyyy_xyyzzzz_1, \
                                         tg_xyyy_xyyzzzzz_0, tg_xyyy_xyyzzzzz_1, tg_xyyy_xyzzzzz_1, tg_xyyy_xyzzzzzz_0, \
                                         tg_xyyy_xyzzzzzz_1, tg_xyyy_xzzzzzz_1, tg_xyyy_xzzzzzzz_0, tg_xyyy_xzzzzzzz_1, \
                                         tg_xyyy_yyyyyyy_1, tg_xyyy_yyyyyyyy_0, tg_xyyy_yyyyyyyy_1, tg_xyyy_yyyyyyyz_0, \
                                         tg_xyyy_yyyyyyyz_1, tg_xyyy_yyyyyyz_1, tg_xyyy_yyyyyyzz_0, tg_xyyy_yyyyyyzz_1, \
                                         tg_xyyy_yyyyyzz_1, tg_xyyy_yyyyyzzz_0, tg_xyyy_yyyyyzzz_1, tg_xyyy_yyyyzzz_1, \
                                         tg_xyyy_yyyyzzzz_0, tg_xyyy_yyyyzzzz_1, tg_xyyy_yyyzzzz_1, tg_xyyy_yyyzzzzz_0, \
                                         tg_xyyy_yyyzzzzz_1, tg_xyyy_yyzzzzz_1, tg_xyyy_yyzzzzzz_0, tg_xyyy_yyzzzzzz_1, \
                                         tg_xyyy_yzzzzzz_1, tg_xyyy_yzzzzzzz_0, tg_xyyy_yzzzzzzz_1, tg_xyyy_zzzzzzz_1, \
                                         tg_xyyy_zzzzzzzz_0, tg_xyyy_zzzzzzzz_1, tg_xyyz_xxxxxxx_1, tg_xyyz_xxxxxxxx_0, \
                                         tg_xyyz_xxxxxxxx_1, tg_xyyz_xxxxxxxy_0, tg_xyyz_xxxxxxxy_1, tg_xyyz_xxxxxxxz_0, \
                                         tg_xyyz_xxxxxxxz_1, tg_xyyz_xxxxxxy_1, tg_xyyz_xxxxxxyy_0, tg_xyyz_xxxxxxyy_1, \
                                         tg_xyyz_xxxxxxyz_0, tg_xyyz_xxxxxxyz_1, tg_xyyz_xxxxxxz_1, tg_xyyz_xxxxxxzz_0, \
                                         tg_xyyz_xxxxxxzz_1, tg_xyyz_xxxxxyy_1, tg_xyyz_xxxxxyyy_0, tg_xyyz_xxxxxyyy_1, \
                                         tg_xyyz_xxxxxyyz_0, tg_xyyz_xxxxxyyz_1, tg_xyyz_xxxxxyz_1, tg_xyyz_xxxxxyzz_0, \
                                         tg_xyyz_xxxxxyzz_1, tg_xyyz_xxxxxzz_1, tg_xyyz_xxxxxzzz_0, tg_xyyz_xxxxxzzz_1, \
                                         tg_xyyz_xxxxyyy_1, tg_xyyz_xxxxyyyy_0, tg_xyyz_xxxxyyyy_1, tg_xyyz_xxxxyyyz_0, \
                                         tg_xyyz_xxxxyyyz_1, tg_xyyz_xxxxyyz_1, tg_xyyz_xxxxyyzz_0, tg_xyyz_xxxxyyzz_1, \
                                         tg_xyyz_xxxxyzz_1, tg_xyyz_xxxxyzzz_0, tg_xyyz_xxxxyzzz_1, tg_xyyz_xxxxzzz_1, \
                                         tg_xyyz_xxxxzzzz_0, tg_xyyz_xxxxzzzz_1, tg_xyyz_xxxyyyy_1, tg_xyyz_xxxyyyyy_0, \
                                         tg_xyyz_xxxyyyyy_1, tg_xyyz_xxxyyyyz_0, tg_xyyz_xxxyyyyz_1, tg_xyyz_xxxyyyz_1, \
                                         tg_xyyz_xxxyyyzz_0, tg_xyyz_xxxyyyzz_1, tg_xyyz_xxxyyzz_1, tg_xyyz_xxxyyzzz_0, \
                                         tg_xyyz_xxxyyzzz_1, tg_xyyz_xxxyzzz_1, tg_xyyz_xxxyzzzz_0, tg_xyyz_xxxyzzzz_1, \
                                         tg_xyyz_xxxzzzz_1, tg_xyyz_xxxzzzzz_0, tg_xyyz_xxxzzzzz_1, tg_xyyz_xxyyyyy_1, \
                                         tg_xyyz_xxyyyyyy_0, tg_xyyz_xxyyyyyy_1, tg_xyyz_xxyyyyyz_0, tg_xyyz_xxyyyyyz_1, \
                                         tg_xyyz_xxyyyyz_1, tg_xyyz_xxyyyyzz_0, tg_xyyz_xxyyyyzz_1, tg_xyyz_xxyyyzz_1, \
                                         tg_xyyz_xxyyyzzz_0, tg_xyyz_xxyyyzzz_1, tg_xyyz_xxyyzzz_1, tg_xyyz_xxyyzzzz_0, \
                                         tg_xyyz_xxyyzzzz_1, tg_xyyz_xxyzzzz_1, tg_xyyz_xxyzzzzz_0, tg_xyyz_xxyzzzzz_1, \
                                         tg_xyyz_xxzzzzz_1, tg_xyyz_xxzzzzzz_0, tg_xyyz_xxzzzzzz_1, tg_xyyz_xyyyyyy_1, \
                                         tg_xyyz_xyyyyyyy_0, tg_xyyz_xyyyyyyy_1, tg_xyyz_xyyyyyyz_0, tg_xyyz_xyyyyyyz_1, \
                                         tg_xyyz_xyyyyyz_1, tg_xyyz_xyyyyyzz_0, tg_xyyz_xyyyyyzz_1, tg_xyyz_xyyyyzz_1, \
                                         tg_xyyz_xyyyyzzz_0, tg_xyyz_xyyyyzzz_1, tg_xyyz_xyyyzzz_1, tg_xyyz_xyyyzzzz_0, \
                                         tg_xyyz_xyyyzzzz_1, tg_xyyz_xyyzzzz_1, tg_xyyz_xyyzzzzz_0, tg_xyyz_xyyzzzzz_1, \
                                         tg_xyyz_xyzzzzz_1, tg_xyyz_xyzzzzzz_0, tg_xyyz_xyzzzzzz_1, tg_xyyz_xzzzzzz_1, \
                                         tg_xyyz_xzzzzzzz_0, tg_xyyz_xzzzzzzz_1, tg_xyyz_yyyyyyy_1, tg_xyyz_yyyyyyyy_0, \
                                         tg_xyyz_yyyyyyyy_1, tg_xyyz_yyyyyyyz_0, tg_xyyz_yyyyyyyz_1, tg_xyyz_yyyyyyz_1, \
                                         tg_xyyz_yyyyyyzz_0, tg_xyyz_yyyyyyzz_1, tg_xyyz_yyyyyzz_1, tg_xyyz_yyyyyzzz_0, \
                                         tg_xyyz_yyyyyzzz_1, tg_xyyz_yyyyzzz_1, tg_xyyz_yyyyzzzz_0, tg_xyyz_yyyyzzzz_1, \
                                         tg_xyyz_yyyzzzz_1, tg_xyyz_yyyzzzzz_0, tg_xyyz_yyyzzzzz_1, tg_xyyz_yyzzzzz_1, \
                                         tg_xyyz_yyzzzzzz_0, tg_xyyz_yyzzzzzz_1, tg_xyyz_yzzzzzz_1, tg_xyyz_yzzzzzzz_0, \
                                         tg_xyyz_yzzzzzzz_1, tg_xyyz_zzzzzzz_1, tg_xyyz_zzzzzzzz_0, tg_xyyz_zzzzzzzz_1, \
                                         tg_xyzz_xxxxxxx_1, tg_xyzz_xxxxxxxx_0, tg_xyzz_xxxxxxxx_1, tg_xyzz_xxxxxxxy_0, \
                                         tg_xyzz_xxxxxxxy_1, tg_xyzz_xxxxxxxz_0, tg_xyzz_xxxxxxxz_1, tg_xyzz_xxxxxxy_1, \
                                         tg_xyzz_xxxxxxyy_0, tg_xyzz_xxxxxxyy_1, tg_xyzz_xxxxxxyz_0, tg_xyzz_xxxxxxyz_1, \
                                         tg_xyzz_xxxxxxz_1, tg_xyzz_xxxxxxzz_0, tg_xyzz_xxxxxxzz_1, tg_xyzz_xxxxxyy_1, \
                                         tg_xyzz_xxxxxyyy_0, tg_xyzz_xxxxxyyy_1, tg_xyzz_xxxxxyyz_0, tg_xyzz_xxxxxyyz_1, \
                                         tg_xyzz_xxxxxyz_1, tg_xyzz_xxxxxyzz_0, tg_xyzz_xxxxxyzz_1, tg_xyzz_xxxxxzz_1, \
                                         tg_xyzz_xxxxxzzz_0, tg_xyzz_xxxxxzzz_1, tg_xyzz_xxxxyyy_1, tg_xyzz_xxxxyyyy_0, \
                                         tg_xyzz_xxxxyyyy_1, tg_xyzz_xxxxyyyz_0, tg_xyzz_xxxxyyyz_1, tg_xyzz_xxxxyyz_1, \
                                         tg_xyzz_xxxxyyzz_0, tg_xyzz_xxxxyyzz_1, tg_xyzz_xxxxyzz_1, tg_xyzz_xxxxyzzz_0, \
                                         tg_xyzz_xxxxyzzz_1, tg_xyzz_xxxxzzz_1, tg_xyzz_xxxxzzzz_0, tg_xyzz_xxxxzzzz_1, \
                                         tg_xyzz_xxxyyyy_1, tg_xyzz_xxxyyyyy_0, tg_xyzz_xxxyyyyy_1, tg_xyzz_xxxyyyyz_0, \
                                         tg_xyzz_xxxyyyyz_1, tg_xyzz_xxxyyyz_1, tg_xyzz_xxxyyyzz_0, tg_xyzz_xxxyyyzz_1, \
                                         tg_xyzz_xxxyyzz_1, tg_xyzz_xxxyyzzz_0, tg_xyzz_xxxyyzzz_1, tg_xyzz_xxxyzzz_1, \
                                         tg_xyzz_xxxyzzzz_0, tg_xyzz_xxxyzzzz_1, tg_xyzz_xxxzzzz_1, tg_xyzz_xxyyyyy_1, \
                                         tg_xyzz_xxyyyyz_1, tg_xyzz_xxyyyzz_1, tg_xyzz_xxyyzzz_1, tg_xyzz_xxyzzzz_1, \
                                         tg_yyy_xxxyyyyy_0, tg_yyy_xxxyyyyy_1, tg_yyy_xxxyyyyz_0, tg_yyy_xxxyyyyz_1, \
                                         tg_yyy_xxxyyyzz_0, tg_yyy_xxxyyyzz_1, tg_yyy_xxxyyzzz_0, tg_yyy_xxxyyzzz_1, \
                                         tg_yyy_xxxyzzzz_0, tg_yyy_xxxyzzzz_1, tg_yyy_xxxzzzzz_0, tg_yyy_xxxzzzzz_1, \
                                         tg_yyy_xxyyyyyy_0, tg_yyy_xxyyyyyy_1, tg_yyy_xxyyyyyz_0, tg_yyy_xxyyyyyz_1, \
                                         tg_yyy_xxyyyyzz_0, tg_yyy_xxyyyyzz_1, tg_yyy_xxyyyzzz_0, tg_yyy_xxyyyzzz_1, \
                                         tg_yyy_xxyyzzzz_0, tg_yyy_xxyyzzzz_1, tg_yyy_xxyzzzzz_0, tg_yyy_xxyzzzzz_1, \
                                         tg_yyy_xxzzzzzz_0, tg_yyy_xxzzzzzz_1, tg_yyy_xyyyyyyy_0, tg_yyy_xyyyyyyy_1, \
                                         tg_yyy_xyyyyyyz_0, tg_yyy_xyyyyyyz_1, tg_yyy_xyyyyyzz_0, tg_yyy_xyyyyyzz_1, \
                                         tg_yyy_xyyyyzzz_0, tg_yyy_xyyyyzzz_1, tg_yyy_xyyyzzzz_0, tg_yyy_xyyyzzzz_1, \
                                         tg_yyy_xyyzzzzz_0, tg_yyy_xyyzzzzz_1, tg_yyy_xyzzzzzz_0, tg_yyy_xyzzzzzz_1, \
                                         tg_yyy_xzzzzzzz_0, tg_yyy_xzzzzzzz_1, tg_yyy_yyyyyyyy_0, tg_yyy_yyyyyyyy_1, \
                                         tg_yyy_yyyyyyyz_0, tg_yyy_yyyyyyyz_1, tg_yyy_yyyyyyzz_0, tg_yyy_yyyyyyzz_1, \
                                         tg_yyy_yyyyyzzz_0, tg_yyy_yyyyyzzz_1, tg_yyy_yyyyzzzz_0, tg_yyy_yyyyzzzz_1, \
                                         tg_yyy_yyyzzzzz_0, tg_yyy_yyyzzzzz_1, tg_yyy_yyzzzzzz_0, tg_yyy_yyzzzzzz_1, \
                                         tg_yyy_yzzzzzzz_0, tg_yyy_yzzzzzzz_1, tg_yyy_zzzzzzzz_0, tg_yyy_zzzzzzzz_1, \
                                         tg_yyz_xxxxxxxx_0, tg_yyz_xxxxxxxx_1, tg_yyz_xxxxxxxy_0, tg_yyz_xxxxxxxy_1, \
                                         tg_yyz_xxxxxxxz_0, tg_yyz_xxxxxxxz_1, tg_yyz_xxxxxxyy_0, tg_yyz_xxxxxxyy_1, \
                                         tg_yyz_xxxxxxyz_0, tg_yyz_xxxxxxyz_1, tg_yyz_xxxxxxzz_0, tg_yyz_xxxxxxzz_1, \
                                         tg_yyz_xxxxxyyy_0, tg_yyz_xxxxxyyy_1, tg_yyz_xxxxxyyz_0, tg_yyz_xxxxxyyz_1, \
                                         tg_yyz_xxxxxyzz_0, tg_yyz_xxxxxyzz_1, tg_yyz_xxxxxzzz_0, tg_yyz_xxxxxzzz_1, \
                                         tg_yyz_xxxxyyyy_0, tg_yyz_xxxxyyyy_1, tg_yyz_xxxxyyyz_0, tg_yyz_xxxxyyyz_1, \
                                         tg_yyz_xxxxyyzz_0, tg_yyz_xxxxyyzz_1, tg_yyz_xxxxyzzz_0, tg_yyz_xxxxyzzz_1, \
                                         tg_yyz_xxxxzzzz_0, tg_yyz_xxxxzzzz_1, tg_yyz_xxxyyyyy_0, tg_yyz_xxxyyyyy_1, \
                                         tg_yyz_xxxyyyyz_0, tg_yyz_xxxyyyyz_1, tg_yyz_xxxyyyzz_0, tg_yyz_xxxyyyzz_1, \
                                         tg_yyz_xxxyyzzz_0, tg_yyz_xxxyyzzz_1, tg_yyz_xxxyzzzz_0, tg_yyz_xxxyzzzz_1, \
                                         tg_yyz_xxxzzzzz_0, tg_yyz_xxxzzzzz_1, tg_yyz_xxyyyyyy_0, tg_yyz_xxyyyyyy_1, \
                                         tg_yyz_xxyyyyyz_0, tg_yyz_xxyyyyyz_1, tg_yyz_xxyyyyzz_0, tg_yyz_xxyyyyzz_1, \
                                         tg_yyz_xxyyyzzz_0, tg_yyz_xxyyyzzz_1, tg_yyz_xxyyzzzz_0, tg_yyz_xxyyzzzz_1, \
                                         tg_yyz_xxyzzzzz_0, tg_yyz_xxyzzzzz_1, tg_yyz_xxzzzzzz_0, tg_yyz_xxzzzzzz_1, \
                                         tg_yyz_xyyyyyyy_0, tg_yyz_xyyyyyyy_1, tg_yyz_xyyyyyyz_0, tg_yyz_xyyyyyyz_1, \
                                         tg_yyz_xyyyyyzz_0, tg_yyz_xyyyyyzz_1, tg_yyz_xyyyyzzz_0, tg_yyz_xyyyyzzz_1, \
                                         tg_yyz_xyyyzzzz_0, tg_yyz_xyyyzzzz_1, tg_yyz_xyyzzzzz_0, tg_yyz_xyyzzzzz_1, \
                                         tg_yyz_xyzzzzzz_0, tg_yyz_xyzzzzzz_1, tg_yyz_xzzzzzzz_0, tg_yyz_xzzzzzzz_1, \
                                         tg_yyz_yyyyyyyy_0, tg_yyz_yyyyyyyy_1, tg_yyz_yyyyyyyz_0, tg_yyz_yyyyyyyz_1, \
                                         tg_yyz_yyyyyyzz_0, tg_yyz_yyyyyyzz_1, tg_yyz_yyyyyzzz_0, tg_yyz_yyyyyzzz_1, \
                                         tg_yyz_yyyyzzzz_0, tg_yyz_yyyyzzzz_1, tg_yyz_yyyzzzzz_0, tg_yyz_yyyzzzzz_1, \
                                         tg_yyz_yyzzzzzz_0, tg_yyz_yyzzzzzz_1, tg_yyz_yzzzzzzz_0, tg_yyz_yzzzzzzz_1, \
                                         tg_yyz_zzzzzzzz_0, tg_yyz_zzzzzzzz_1, tg_yzz_xxxxxxxx_0, tg_yzz_xxxxxxxx_1, \
                                         tg_yzz_xxxxxxxy_0, tg_yzz_xxxxxxxy_1, tg_yzz_xxxxxxxz_0, tg_yzz_xxxxxxxz_1, \
                                         tg_yzz_xxxxxxyy_0, tg_yzz_xxxxxxyy_1, tg_yzz_xxxxxxyz_0, tg_yzz_xxxxxxyz_1, \
                                         tg_yzz_xxxxxxzz_0, tg_yzz_xxxxxxzz_1, tg_yzz_xxxxxyyy_0, tg_yzz_xxxxxyyy_1, \
                                         tg_yzz_xxxxxyyz_0, tg_yzz_xxxxxyyz_1, tg_yzz_xxxxxyzz_0, tg_yzz_xxxxxyzz_1, \
                                         tg_yzz_xxxxxzzz_0, tg_yzz_xxxxxzzz_1, tg_yzz_xxxxyyyy_0, tg_yzz_xxxxyyyy_1, \
                                         tg_yzz_xxxxyyyz_0, tg_yzz_xxxxyyyz_1, tg_yzz_xxxxyyzz_0, tg_yzz_xxxxyyzz_1, \
                                         tg_yzz_xxxxyzzz_0, tg_yzz_xxxxyzzz_1, tg_yzz_xxxxzzzz_0, tg_yzz_xxxxzzzz_1, \
                                         tg_yzz_xxxyyyyy_0, tg_yzz_xxxyyyyy_1, tg_yzz_xxxyyyyz_0, tg_yzz_xxxyyyyz_1, \
                                         tg_yzz_xxxyyyzz_0, tg_yzz_xxxyyyzz_1, tg_yzz_xxxyyzzz_0, tg_yzz_xxxyyzzz_1, \
                                         tg_yzz_xxxyzzzz_0, tg_yzz_xxxyzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxyyy_xxxyyyyy_0[j] = pb_x * tg_xyyy_xxxyyyyy_0[j] + fr * tg_xyyy_xxxyyyyy_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxyyyyy_0[j] - tg_yyy_xxxyyyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyy_xxyyyyy_1[j];

                    tg_xxyyy_xxxyyyyz_0[j] = pb_x * tg_xyyy_xxxyyyyz_0[j] + fr * tg_xyyy_xxxyyyyz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxyyyyz_0[j] - tg_yyy_xxxyyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyy_xxyyyyz_1[j];

                    tg_xxyyy_xxxyyyzz_0[j] = pb_x * tg_xyyy_xxxyyyzz_0[j] + fr * tg_xyyy_xxxyyyzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxyyyzz_0[j] - tg_yyy_xxxyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyy_xxyyyzz_1[j];

                    tg_xxyyy_xxxyyzzz_0[j] = pb_x * tg_xyyy_xxxyyzzz_0[j] + fr * tg_xyyy_xxxyyzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxyyzzz_0[j] - tg_yyy_xxxyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyy_xxyyzzz_1[j];

                    tg_xxyyy_xxxyzzzz_0[j] = pb_x * tg_xyyy_xxxyzzzz_0[j] + fr * tg_xyyy_xxxyzzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxyzzzz_0[j] - tg_yyy_xxxyzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyy_xxyzzzz_1[j];

                    tg_xxyyy_xxxzzzzz_0[j] = pb_x * tg_xyyy_xxxzzzzz_0[j] + fr * tg_xyyy_xxxzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxzzzzz_0[j] - tg_yyy_xxxzzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyy_xxzzzzz_1[j];

                    tg_xxyyy_xxyyyyyy_0[j] = pb_x * tg_xyyy_xxyyyyyy_0[j] + fr * tg_xyyy_xxyyyyyy_1[j] + 0.5 * fl1_fx * (tg_yyy_xxyyyyyy_0[j] - tg_yyy_xxyyyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyyy_xyyyyyy_1[j];

                    tg_xxyyy_xxyyyyyz_0[j] = pb_x * tg_xyyy_xxyyyyyz_0[j] + fr * tg_xyyy_xxyyyyyz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxyyyyyz_0[j] - tg_yyy_xxyyyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyyy_xyyyyyz_1[j];

                    tg_xxyyy_xxyyyyzz_0[j] = pb_x * tg_xyyy_xxyyyyzz_0[j] + fr * tg_xyyy_xxyyyyzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxyyyyzz_0[j] - tg_yyy_xxyyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyy_xyyyyzz_1[j];

                    tg_xxyyy_xxyyyzzz_0[j] = pb_x * tg_xyyy_xxyyyzzz_0[j] + fr * tg_xyyy_xxyyyzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxyyyzzz_0[j] - tg_yyy_xxyyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyy_xyyyzzz_1[j];

                    tg_xxyyy_xxyyzzzz_0[j] = pb_x * tg_xyyy_xxyyzzzz_0[j] + fr * tg_xyyy_xxyyzzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxyyzzzz_0[j] - tg_yyy_xxyyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyy_xyyzzzz_1[j];

                    tg_xxyyy_xxyzzzzz_0[j] = pb_x * tg_xyyy_xxyzzzzz_0[j] + fr * tg_xyyy_xxyzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxyzzzzz_0[j] - tg_yyy_xxyzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyy_xyzzzzz_1[j];

                    tg_xxyyy_xxzzzzzz_0[j] = pb_x * tg_xyyy_xxzzzzzz_0[j] + fr * tg_xyyy_xxzzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxzzzzzz_0[j] - tg_yyy_xxzzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyy_xzzzzzz_1[j];

                    tg_xxyyy_xyyyyyyy_0[j] = pb_x * tg_xyyy_xyyyyyyy_0[j] + fr * tg_xyyy_xyyyyyyy_1[j] + 0.5 * fl1_fx * (tg_yyy_xyyyyyyy_0[j] - tg_yyy_xyyyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyy_yyyyyyy_1[j];

                    tg_xxyyy_xyyyyyyz_0[j] = pb_x * tg_xyyy_xyyyyyyz_0[j] + fr * tg_xyyy_xyyyyyyz_1[j] + 0.5 * fl1_fx * (tg_yyy_xyyyyyyz_0[j] - tg_yyy_xyyyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyy_yyyyyyz_1[j];

                    tg_xxyyy_xyyyyyzz_0[j] = pb_x * tg_xyyy_xyyyyyzz_0[j] + fr * tg_xyyy_xyyyyyzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xyyyyyzz_0[j] - tg_yyy_xyyyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyy_yyyyyzz_1[j];

                    tg_xxyyy_xyyyyzzz_0[j] = pb_x * tg_xyyy_xyyyyzzz_0[j] + fr * tg_xyyy_xyyyyzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xyyyyzzz_0[j] - tg_yyy_xyyyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyy_yyyyzzz_1[j];

                    tg_xxyyy_xyyyzzzz_0[j] = pb_x * tg_xyyy_xyyyzzzz_0[j] + fr * tg_xyyy_xyyyzzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xyyyzzzz_0[j] - tg_yyy_xyyyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyy_yyyzzzz_1[j];

                    tg_xxyyy_xyyzzzzz_0[j] = pb_x * tg_xyyy_xyyzzzzz_0[j] + fr * tg_xyyy_xyyzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xyyzzzzz_0[j] - tg_yyy_xyyzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyy_yyzzzzz_1[j];

                    tg_xxyyy_xyzzzzzz_0[j] = pb_x * tg_xyyy_xyzzzzzz_0[j] + fr * tg_xyyy_xyzzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xyzzzzzz_0[j] - tg_yyy_xyzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyy_yzzzzzz_1[j];

                    tg_xxyyy_xzzzzzzz_0[j] = pb_x * tg_xyyy_xzzzzzzz_0[j] + fr * tg_xyyy_xzzzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xzzzzzzz_0[j] - tg_yyy_xzzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyy_zzzzzzz_1[j];

                    tg_xxyyy_yyyyyyyy_0[j] = pb_x * tg_xyyy_yyyyyyyy_0[j] + fr * tg_xyyy_yyyyyyyy_1[j] + 0.5 * fl1_fx * (tg_yyy_yyyyyyyy_0[j] - tg_yyy_yyyyyyyy_1[j] * fl1_fza);

                    tg_xxyyy_yyyyyyyz_0[j] = pb_x * tg_xyyy_yyyyyyyz_0[j] + fr * tg_xyyy_yyyyyyyz_1[j] + 0.5 * fl1_fx * (tg_yyy_yyyyyyyz_0[j] - tg_yyy_yyyyyyyz_1[j] * fl1_fza);

                    tg_xxyyy_yyyyyyzz_0[j] = pb_x * tg_xyyy_yyyyyyzz_0[j] + fr * tg_xyyy_yyyyyyzz_1[j] + 0.5 * fl1_fx * (tg_yyy_yyyyyyzz_0[j] - tg_yyy_yyyyyyzz_1[j] * fl1_fza);

                    tg_xxyyy_yyyyyzzz_0[j] = pb_x * tg_xyyy_yyyyyzzz_0[j] + fr * tg_xyyy_yyyyyzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_yyyyyzzz_0[j] - tg_yyy_yyyyyzzz_1[j] * fl1_fza);

                    tg_xxyyy_yyyyzzzz_0[j] = pb_x * tg_xyyy_yyyyzzzz_0[j] + fr * tg_xyyy_yyyyzzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_yyyyzzzz_0[j] - tg_yyy_yyyyzzzz_1[j] * fl1_fza);

                    tg_xxyyy_yyyzzzzz_0[j] = pb_x * tg_xyyy_yyyzzzzz_0[j] + fr * tg_xyyy_yyyzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_yyyzzzzz_0[j] - tg_yyy_yyyzzzzz_1[j] * fl1_fza);

                    tg_xxyyy_yyzzzzzz_0[j] = pb_x * tg_xyyy_yyzzzzzz_0[j] + fr * tg_xyyy_yyzzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_yyzzzzzz_0[j] - tg_yyy_yyzzzzzz_1[j] * fl1_fza);

                    tg_xxyyy_yzzzzzzz_0[j] = pb_x * tg_xyyy_yzzzzzzz_0[j] + fr * tg_xyyy_yzzzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_yzzzzzzz_0[j] - tg_yyy_yzzzzzzz_1[j] * fl1_fza);

                    tg_xxyyy_zzzzzzzz_0[j] = pb_x * tg_xyyy_zzzzzzzz_0[j] + fr * tg_xyyy_zzzzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_zzzzzzzz_0[j] - tg_yyy_zzzzzzzz_1[j] * fl1_fza);

                    tg_xxyyz_xxxxxxxx_0[j] = pb_x * tg_xyyz_xxxxxxxx_0[j] + fr * tg_xyyz_xxxxxxxx_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxxxxxx_0[j] - tg_yyz_xxxxxxxx_1[j] * fl1_fza) + 4.0 * fl1_fxn * tg_xyyz_xxxxxxx_1[j];

                    tg_xxyyz_xxxxxxxy_0[j] = pb_x * tg_xyyz_xxxxxxxy_0[j] + fr * tg_xyyz_xxxxxxxy_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxxxxxy_0[j] - tg_yyz_xxxxxxxy_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_xyyz_xxxxxxy_1[j];

                    tg_xxyyz_xxxxxxxz_0[j] = pb_x * tg_xyyz_xxxxxxxz_0[j] + fr * tg_xyyz_xxxxxxxz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxxxxxz_0[j] - tg_yyz_xxxxxxxz_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_xyyz_xxxxxxz_1[j];

                    tg_xxyyz_xxxxxxyy_0[j] = pb_x * tg_xyyz_xxxxxxyy_0[j] + fr * tg_xyyz_xxxxxxyy_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxxxxyy_0[j] - tg_yyz_xxxxxxyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xyyz_xxxxxyy_1[j];

                    tg_xxyyz_xxxxxxyz_0[j] = pb_x * tg_xyyz_xxxxxxyz_0[j] + fr * tg_xyyz_xxxxxxyz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxxxxyz_0[j] - tg_yyz_xxxxxxyz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xyyz_xxxxxyz_1[j];

                    tg_xxyyz_xxxxxxzz_0[j] = pb_x * tg_xyyz_xxxxxxzz_0[j] + fr * tg_xyyz_xxxxxxzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxxxxzz_0[j] - tg_yyz_xxxxxxzz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xyyz_xxxxxzz_1[j];

                    tg_xxyyz_xxxxxyyy_0[j] = pb_x * tg_xyyz_xxxxxyyy_0[j] + fr * tg_xyyz_xxxxxyyy_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxxxyyy_0[j] - tg_yyz_xxxxxyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyz_xxxxyyy_1[j];

                    tg_xxyyz_xxxxxyyz_0[j] = pb_x * tg_xyyz_xxxxxyyz_0[j] + fr * tg_xyyz_xxxxxyyz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxxxyyz_0[j] - tg_yyz_xxxxxyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyz_xxxxyyz_1[j];

                    tg_xxyyz_xxxxxyzz_0[j] = pb_x * tg_xyyz_xxxxxyzz_0[j] + fr * tg_xyyz_xxxxxyzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxxxyzz_0[j] - tg_yyz_xxxxxyzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyz_xxxxyzz_1[j];

                    tg_xxyyz_xxxxxzzz_0[j] = pb_x * tg_xyyz_xxxxxzzz_0[j] + fr * tg_xyyz_xxxxxzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxxxzzz_0[j] - tg_yyz_xxxxxzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyz_xxxxzzz_1[j];

                    tg_xxyyz_xxxxyyyy_0[j] = pb_x * tg_xyyz_xxxxyyyy_0[j] + fr * tg_xyyz_xxxxyyyy_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxxyyyy_0[j] - tg_yyz_xxxxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyz_xxxyyyy_1[j];

                    tg_xxyyz_xxxxyyyz_0[j] = pb_x * tg_xyyz_xxxxyyyz_0[j] + fr * tg_xyyz_xxxxyyyz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxxyyyz_0[j] - tg_yyz_xxxxyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyz_xxxyyyz_1[j];

                    tg_xxyyz_xxxxyyzz_0[j] = pb_x * tg_xyyz_xxxxyyzz_0[j] + fr * tg_xyyz_xxxxyyzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxxyyzz_0[j] - tg_yyz_xxxxyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyz_xxxyyzz_1[j];

                    tg_xxyyz_xxxxyzzz_0[j] = pb_x * tg_xyyz_xxxxyzzz_0[j] + fr * tg_xyyz_xxxxyzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxxyzzz_0[j] - tg_yyz_xxxxyzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyz_xxxyzzz_1[j];

                    tg_xxyyz_xxxxzzzz_0[j] = pb_x * tg_xyyz_xxxxzzzz_0[j] + fr * tg_xyyz_xxxxzzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxxzzzz_0[j] - tg_yyz_xxxxzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyz_xxxzzzz_1[j];

                    tg_xxyyz_xxxyyyyy_0[j] = pb_x * tg_xyyz_xxxyyyyy_0[j] + fr * tg_xyyz_xxxyyyyy_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxyyyyy_0[j] - tg_yyz_xxxyyyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyz_xxyyyyy_1[j];

                    tg_xxyyz_xxxyyyyz_0[j] = pb_x * tg_xyyz_xxxyyyyz_0[j] + fr * tg_xyyz_xxxyyyyz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxyyyyz_0[j] - tg_yyz_xxxyyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyz_xxyyyyz_1[j];

                    tg_xxyyz_xxxyyyzz_0[j] = pb_x * tg_xyyz_xxxyyyzz_0[j] + fr * tg_xyyz_xxxyyyzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxyyyzz_0[j] - tg_yyz_xxxyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyz_xxyyyzz_1[j];

                    tg_xxyyz_xxxyyzzz_0[j] = pb_x * tg_xyyz_xxxyyzzz_0[j] + fr * tg_xyyz_xxxyyzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxyyzzz_0[j] - tg_yyz_xxxyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyz_xxyyzzz_1[j];

                    tg_xxyyz_xxxyzzzz_0[j] = pb_x * tg_xyyz_xxxyzzzz_0[j] + fr * tg_xyyz_xxxyzzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxyzzzz_0[j] - tg_yyz_xxxyzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyz_xxyzzzz_1[j];

                    tg_xxyyz_xxxzzzzz_0[j] = pb_x * tg_xyyz_xxxzzzzz_0[j] + fr * tg_xyyz_xxxzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxzzzzz_0[j] - tg_yyz_xxxzzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyz_xxzzzzz_1[j];

                    tg_xxyyz_xxyyyyyy_0[j] = pb_x * tg_xyyz_xxyyyyyy_0[j] + fr * tg_xyyz_xxyyyyyy_1[j] + 0.5 * fl1_fx * (tg_yyz_xxyyyyyy_0[j] - tg_yyz_xxyyyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyyz_xyyyyyy_1[j];

                    tg_xxyyz_xxyyyyyz_0[j] = pb_x * tg_xyyz_xxyyyyyz_0[j] + fr * tg_xyyz_xxyyyyyz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxyyyyyz_0[j] - tg_yyz_xxyyyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyyz_xyyyyyz_1[j];

                    tg_xxyyz_xxyyyyzz_0[j] = pb_x * tg_xyyz_xxyyyyzz_0[j] + fr * tg_xyyz_xxyyyyzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxyyyyzz_0[j] - tg_yyz_xxyyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyz_xyyyyzz_1[j];

                    tg_xxyyz_xxyyyzzz_0[j] = pb_x * tg_xyyz_xxyyyzzz_0[j] + fr * tg_xyyz_xxyyyzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxyyyzzz_0[j] - tg_yyz_xxyyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyz_xyyyzzz_1[j];

                    tg_xxyyz_xxyyzzzz_0[j] = pb_x * tg_xyyz_xxyyzzzz_0[j] + fr * tg_xyyz_xxyyzzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxyyzzzz_0[j] - tg_yyz_xxyyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyz_xyyzzzz_1[j];

                    tg_xxyyz_xxyzzzzz_0[j] = pb_x * tg_xyyz_xxyzzzzz_0[j] + fr * tg_xyyz_xxyzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxyzzzzz_0[j] - tg_yyz_xxyzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyz_xyzzzzz_1[j];

                    tg_xxyyz_xxzzzzzz_0[j] = pb_x * tg_xyyz_xxzzzzzz_0[j] + fr * tg_xyyz_xxzzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxzzzzzz_0[j] - tg_yyz_xxzzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyz_xzzzzzz_1[j];

                    tg_xxyyz_xyyyyyyy_0[j] = pb_x * tg_xyyz_xyyyyyyy_0[j] + fr * tg_xyyz_xyyyyyyy_1[j] + 0.5 * fl1_fx * (tg_yyz_xyyyyyyy_0[j] - tg_yyz_xyyyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyz_yyyyyyy_1[j];

                    tg_xxyyz_xyyyyyyz_0[j] = pb_x * tg_xyyz_xyyyyyyz_0[j] + fr * tg_xyyz_xyyyyyyz_1[j] + 0.5 * fl1_fx * (tg_yyz_xyyyyyyz_0[j] - tg_yyz_xyyyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyz_yyyyyyz_1[j];

                    tg_xxyyz_xyyyyyzz_0[j] = pb_x * tg_xyyz_xyyyyyzz_0[j] + fr * tg_xyyz_xyyyyyzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xyyyyyzz_0[j] - tg_yyz_xyyyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyz_yyyyyzz_1[j];

                    tg_xxyyz_xyyyyzzz_0[j] = pb_x * tg_xyyz_xyyyyzzz_0[j] + fr * tg_xyyz_xyyyyzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xyyyyzzz_0[j] - tg_yyz_xyyyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyz_yyyyzzz_1[j];

                    tg_xxyyz_xyyyzzzz_0[j] = pb_x * tg_xyyz_xyyyzzzz_0[j] + fr * tg_xyyz_xyyyzzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xyyyzzzz_0[j] - tg_yyz_xyyyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyz_yyyzzzz_1[j];

                    tg_xxyyz_xyyzzzzz_0[j] = pb_x * tg_xyyz_xyyzzzzz_0[j] + fr * tg_xyyz_xyyzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xyyzzzzz_0[j] - tg_yyz_xyyzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyz_yyzzzzz_1[j];

                    tg_xxyyz_xyzzzzzz_0[j] = pb_x * tg_xyyz_xyzzzzzz_0[j] + fr * tg_xyyz_xyzzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xyzzzzzz_0[j] - tg_yyz_xyzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyz_yzzzzzz_1[j];

                    tg_xxyyz_xzzzzzzz_0[j] = pb_x * tg_xyyz_xzzzzzzz_0[j] + fr * tg_xyyz_xzzzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xzzzzzzz_0[j] - tg_yyz_xzzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyz_zzzzzzz_1[j];

                    tg_xxyyz_yyyyyyyy_0[j] = pb_x * tg_xyyz_yyyyyyyy_0[j] + fr * tg_xyyz_yyyyyyyy_1[j] + 0.5 * fl1_fx * (tg_yyz_yyyyyyyy_0[j] - tg_yyz_yyyyyyyy_1[j] * fl1_fza);

                    tg_xxyyz_yyyyyyyz_0[j] = pb_x * tg_xyyz_yyyyyyyz_0[j] + fr * tg_xyyz_yyyyyyyz_1[j] + 0.5 * fl1_fx * (tg_yyz_yyyyyyyz_0[j] - tg_yyz_yyyyyyyz_1[j] * fl1_fza);

                    tg_xxyyz_yyyyyyzz_0[j] = pb_x * tg_xyyz_yyyyyyzz_0[j] + fr * tg_xyyz_yyyyyyzz_1[j] + 0.5 * fl1_fx * (tg_yyz_yyyyyyzz_0[j] - tg_yyz_yyyyyyzz_1[j] * fl1_fza);

                    tg_xxyyz_yyyyyzzz_0[j] = pb_x * tg_xyyz_yyyyyzzz_0[j] + fr * tg_xyyz_yyyyyzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_yyyyyzzz_0[j] - tg_yyz_yyyyyzzz_1[j] * fl1_fza);

                    tg_xxyyz_yyyyzzzz_0[j] = pb_x * tg_xyyz_yyyyzzzz_0[j] + fr * tg_xyyz_yyyyzzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_yyyyzzzz_0[j] - tg_yyz_yyyyzzzz_1[j] * fl1_fza);

                    tg_xxyyz_yyyzzzzz_0[j] = pb_x * tg_xyyz_yyyzzzzz_0[j] + fr * tg_xyyz_yyyzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_yyyzzzzz_0[j] - tg_yyz_yyyzzzzz_1[j] * fl1_fza);

                    tg_xxyyz_yyzzzzzz_0[j] = pb_x * tg_xyyz_yyzzzzzz_0[j] + fr * tg_xyyz_yyzzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_yyzzzzzz_0[j] - tg_yyz_yyzzzzzz_1[j] * fl1_fza);

                    tg_xxyyz_yzzzzzzz_0[j] = pb_x * tg_xyyz_yzzzzzzz_0[j] + fr * tg_xyyz_yzzzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_yzzzzzzz_0[j] - tg_yyz_yzzzzzzz_1[j] * fl1_fza);

                    tg_xxyyz_zzzzzzzz_0[j] = pb_x * tg_xyyz_zzzzzzzz_0[j] + fr * tg_xyyz_zzzzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_zzzzzzzz_0[j] - tg_yyz_zzzzzzzz_1[j] * fl1_fza);

                    tg_xxyzz_xxxxxxxx_0[j] = pb_x * tg_xyzz_xxxxxxxx_0[j] + fr * tg_xyzz_xxxxxxxx_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxxxxxx_0[j] - tg_yzz_xxxxxxxx_1[j] * fl1_fza) + 4.0 * fl1_fxn * tg_xyzz_xxxxxxx_1[j];

                    tg_xxyzz_xxxxxxxy_0[j] = pb_x * tg_xyzz_xxxxxxxy_0[j] + fr * tg_xyzz_xxxxxxxy_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxxxxxy_0[j] - tg_yzz_xxxxxxxy_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_xyzz_xxxxxxy_1[j];

                    tg_xxyzz_xxxxxxxz_0[j] = pb_x * tg_xyzz_xxxxxxxz_0[j] + fr * tg_xyzz_xxxxxxxz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxxxxxz_0[j] - tg_yzz_xxxxxxxz_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_xyzz_xxxxxxz_1[j];

                    tg_xxyzz_xxxxxxyy_0[j] = pb_x * tg_xyzz_xxxxxxyy_0[j] + fr * tg_xyzz_xxxxxxyy_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxxxxyy_0[j] - tg_yzz_xxxxxxyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xyzz_xxxxxyy_1[j];

                    tg_xxyzz_xxxxxxyz_0[j] = pb_x * tg_xyzz_xxxxxxyz_0[j] + fr * tg_xyzz_xxxxxxyz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxxxxyz_0[j] - tg_yzz_xxxxxxyz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xyzz_xxxxxyz_1[j];

                    tg_xxyzz_xxxxxxzz_0[j] = pb_x * tg_xyzz_xxxxxxzz_0[j] + fr * tg_xyzz_xxxxxxzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxxxxzz_0[j] - tg_yzz_xxxxxxzz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xyzz_xxxxxzz_1[j];

                    tg_xxyzz_xxxxxyyy_0[j] = pb_x * tg_xyzz_xxxxxyyy_0[j] + fr * tg_xyzz_xxxxxyyy_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxxxyyy_0[j] - tg_yzz_xxxxxyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyzz_xxxxyyy_1[j];

                    tg_xxyzz_xxxxxyyz_0[j] = pb_x * tg_xyzz_xxxxxyyz_0[j] + fr * tg_xyzz_xxxxxyyz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxxxyyz_0[j] - tg_yzz_xxxxxyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyzz_xxxxyyz_1[j];

                    tg_xxyzz_xxxxxyzz_0[j] = pb_x * tg_xyzz_xxxxxyzz_0[j] + fr * tg_xyzz_xxxxxyzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxxxyzz_0[j] - tg_yzz_xxxxxyzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyzz_xxxxyzz_1[j];

                    tg_xxyzz_xxxxxzzz_0[j] = pb_x * tg_xyzz_xxxxxzzz_0[j] + fr * tg_xyzz_xxxxxzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxxxzzz_0[j] - tg_yzz_xxxxxzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyzz_xxxxzzz_1[j];

                    tg_xxyzz_xxxxyyyy_0[j] = pb_x * tg_xyzz_xxxxyyyy_0[j] + fr * tg_xyzz_xxxxyyyy_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxxyyyy_0[j] - tg_yzz_xxxxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyzz_xxxyyyy_1[j];

                    tg_xxyzz_xxxxyyyz_0[j] = pb_x * tg_xyzz_xxxxyyyz_0[j] + fr * tg_xyzz_xxxxyyyz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxxyyyz_0[j] - tg_yzz_xxxxyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyzz_xxxyyyz_1[j];

                    tg_xxyzz_xxxxyyzz_0[j] = pb_x * tg_xyzz_xxxxyyzz_0[j] + fr * tg_xyzz_xxxxyyzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxxyyzz_0[j] - tg_yzz_xxxxyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyzz_xxxyyzz_1[j];

                    tg_xxyzz_xxxxyzzz_0[j] = pb_x * tg_xyzz_xxxxyzzz_0[j] + fr * tg_xyzz_xxxxyzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxxyzzz_0[j] - tg_yzz_xxxxyzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyzz_xxxyzzz_1[j];

                    tg_xxyzz_xxxxzzzz_0[j] = pb_x * tg_xyzz_xxxxzzzz_0[j] + fr * tg_xyzz_xxxxzzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxxzzzz_0[j] - tg_yzz_xxxxzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyzz_xxxzzzz_1[j];

                    tg_xxyzz_xxxyyyyy_0[j] = pb_x * tg_xyzz_xxxyyyyy_0[j] + fr * tg_xyzz_xxxyyyyy_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxyyyyy_0[j] - tg_yzz_xxxyyyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzz_xxyyyyy_1[j];

                    tg_xxyzz_xxxyyyyz_0[j] = pb_x * tg_xyzz_xxxyyyyz_0[j] + fr * tg_xyzz_xxxyyyyz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxyyyyz_0[j] - tg_yzz_xxxyyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzz_xxyyyyz_1[j];

                    tg_xxyzz_xxxyyyzz_0[j] = pb_x * tg_xyzz_xxxyyyzz_0[j] + fr * tg_xyzz_xxxyyyzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxyyyzz_0[j] - tg_yzz_xxxyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzz_xxyyyzz_1[j];

                    tg_xxyzz_xxxyyzzz_0[j] = pb_x * tg_xyzz_xxxyyzzz_0[j] + fr * tg_xyzz_xxxyyzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxyyzzz_0[j] - tg_yzz_xxxyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzz_xxyyzzz_1[j];

                    tg_xxyzz_xxxyzzzz_0[j] = pb_x * tg_xyzz_xxxyzzzz_0[j] + fr * tg_xyzz_xxxyzzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxyzzzz_0[j] - tg_yzz_xxxyzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzz_xxyzzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSL_380_475(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (380,475)

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
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_8_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_8_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_8_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_8_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_8_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_7_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {7, -1, -1, -1}, 
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

                auto tg_xyzz_xxxzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 380); 

                auto tg_xyzz_xxyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 381); 

                auto tg_xyzz_xxyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 382); 

                auto tg_xyzz_xxyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 383); 

                auto tg_xyzz_xxyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 384); 

                auto tg_xyzz_xxyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 385); 

                auto tg_xyzz_xxyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 386); 

                auto tg_xyzz_xxzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 387); 

                auto tg_xyzz_xyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 388); 

                auto tg_xyzz_xyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 389); 

                auto tg_xyzz_xyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 390); 

                auto tg_xyzz_xyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 391); 

                auto tg_xyzz_xyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 392); 

                auto tg_xyzz_xyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 393); 

                auto tg_xyzz_xyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 394); 

                auto tg_xyzz_xzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 395); 

                auto tg_xyzz_yyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 396); 

                auto tg_xyzz_yyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 397); 

                auto tg_xyzz_yyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 398); 

                auto tg_xyzz_yyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 399); 

                auto tg_xyzz_yyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 400); 

                auto tg_xyzz_yyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 401); 

                auto tg_xyzz_yyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 402); 

                auto tg_xyzz_yzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 403); 

                auto tg_xyzz_zzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 404); 

                auto tg_xzzz_xxxxxxxx_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 405); 

                auto tg_xzzz_xxxxxxxy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 406); 

                auto tg_xzzz_xxxxxxxz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 407); 

                auto tg_xzzz_xxxxxxyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 408); 

                auto tg_xzzz_xxxxxxyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 409); 

                auto tg_xzzz_xxxxxxzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 410); 

                auto tg_xzzz_xxxxxyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 411); 

                auto tg_xzzz_xxxxxyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 412); 

                auto tg_xzzz_xxxxxyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 413); 

                auto tg_xzzz_xxxxxzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 414); 

                auto tg_xzzz_xxxxyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 415); 

                auto tg_xzzz_xxxxyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 416); 

                auto tg_xzzz_xxxxyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 417); 

                auto tg_xzzz_xxxxyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 418); 

                auto tg_xzzz_xxxxzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 419); 

                auto tg_xzzz_xxxyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 420); 

                auto tg_xzzz_xxxyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 421); 

                auto tg_xzzz_xxxyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 422); 

                auto tg_xzzz_xxxyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 423); 

                auto tg_xzzz_xxxyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 424); 

                auto tg_xzzz_xxxzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 425); 

                auto tg_xzzz_xxyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 426); 

                auto tg_xzzz_xxyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 427); 

                auto tg_xzzz_xxyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 428); 

                auto tg_xzzz_xxyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 429); 

                auto tg_xzzz_xxyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 430); 

                auto tg_xzzz_xxyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 431); 

                auto tg_xzzz_xxzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 432); 

                auto tg_xzzz_xyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 433); 

                auto tg_xzzz_xyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 434); 

                auto tg_xzzz_xyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 435); 

                auto tg_xzzz_xyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 436); 

                auto tg_xzzz_xyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 437); 

                auto tg_xzzz_xyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 438); 

                auto tg_xzzz_xyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 439); 

                auto tg_xzzz_xzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 440); 

                auto tg_xzzz_yyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 441); 

                auto tg_xzzz_yyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 442); 

                auto tg_xzzz_yyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 443); 

                auto tg_xzzz_yyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 444); 

                auto tg_xzzz_yyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 445); 

                auto tg_xzzz_yyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 446); 

                auto tg_xzzz_yyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 447); 

                auto tg_xzzz_yzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 448); 

                auto tg_xzzz_zzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 449); 

                auto tg_yyyy_xxxxxxxx_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 450); 

                auto tg_yyyy_xxxxxxxy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 451); 

                auto tg_yyyy_xxxxxxxz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 452); 

                auto tg_yyyy_xxxxxxyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 453); 

                auto tg_yyyy_xxxxxxyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 454); 

                auto tg_yyyy_xxxxxxzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 455); 

                auto tg_yyyy_xxxxxyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 456); 

                auto tg_yyyy_xxxxxyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 457); 

                auto tg_yyyy_xxxxxyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 458); 

                auto tg_yyyy_xxxxxzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 459); 

                auto tg_yyyy_xxxxyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 460); 

                auto tg_yyyy_xxxxyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 461); 

                auto tg_yyyy_xxxxyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 462); 

                auto tg_yyyy_xxxxyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 463); 

                auto tg_yyyy_xxxxzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 464); 

                auto tg_yyyy_xxxyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 465); 

                auto tg_yyyy_xxxyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 466); 

                auto tg_yyyy_xxxyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 467); 

                auto tg_yyyy_xxxyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 468); 

                auto tg_yyyy_xxxyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 469); 

                auto tg_yyyy_xxxzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 470); 

                auto tg_yyyy_xxyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 471); 

                auto tg_yyyy_xxyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 472); 

                auto tg_yyyy_xxyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 473); 

                auto tg_yyyy_xxyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 474); 

                auto tg_xyzz_xxxzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 380); 

                auto tg_xyzz_xxyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 381); 

                auto tg_xyzz_xxyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 382); 

                auto tg_xyzz_xxyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 383); 

                auto tg_xyzz_xxyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 384); 

                auto tg_xyzz_xxyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 385); 

                auto tg_xyzz_xxyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 386); 

                auto tg_xyzz_xxzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 387); 

                auto tg_xyzz_xyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 388); 

                auto tg_xyzz_xyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 389); 

                auto tg_xyzz_xyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 390); 

                auto tg_xyzz_xyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 391); 

                auto tg_xyzz_xyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 392); 

                auto tg_xyzz_xyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 393); 

                auto tg_xyzz_xyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 394); 

                auto tg_xyzz_xzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 395); 

                auto tg_xyzz_yyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 396); 

                auto tg_xyzz_yyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 397); 

                auto tg_xyzz_yyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 398); 

                auto tg_xyzz_yyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 399); 

                auto tg_xyzz_yyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 400); 

                auto tg_xyzz_yyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 401); 

                auto tg_xyzz_yyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 402); 

                auto tg_xyzz_yzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 403); 

                auto tg_xyzz_zzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 404); 

                auto tg_xzzz_xxxxxxxx_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 405); 

                auto tg_xzzz_xxxxxxxy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 406); 

                auto tg_xzzz_xxxxxxxz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 407); 

                auto tg_xzzz_xxxxxxyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 408); 

                auto tg_xzzz_xxxxxxyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 409); 

                auto tg_xzzz_xxxxxxzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 410); 

                auto tg_xzzz_xxxxxyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 411); 

                auto tg_xzzz_xxxxxyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 412); 

                auto tg_xzzz_xxxxxyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 413); 

                auto tg_xzzz_xxxxxzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 414); 

                auto tg_xzzz_xxxxyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 415); 

                auto tg_xzzz_xxxxyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 416); 

                auto tg_xzzz_xxxxyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 417); 

                auto tg_xzzz_xxxxyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 418); 

                auto tg_xzzz_xxxxzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 419); 

                auto tg_xzzz_xxxyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 420); 

                auto tg_xzzz_xxxyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 421); 

                auto tg_xzzz_xxxyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 422); 

                auto tg_xzzz_xxxyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 423); 

                auto tg_xzzz_xxxyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 424); 

                auto tg_xzzz_xxxzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 425); 

                auto tg_xzzz_xxyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 426); 

                auto tg_xzzz_xxyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 427); 

                auto tg_xzzz_xxyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 428); 

                auto tg_xzzz_xxyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 429); 

                auto tg_xzzz_xxyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 430); 

                auto tg_xzzz_xxyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 431); 

                auto tg_xzzz_xxzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 432); 

                auto tg_xzzz_xyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 433); 

                auto tg_xzzz_xyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 434); 

                auto tg_xzzz_xyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 435); 

                auto tg_xzzz_xyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 436); 

                auto tg_xzzz_xyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 437); 

                auto tg_xzzz_xyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 438); 

                auto tg_xzzz_xyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 439); 

                auto tg_xzzz_xzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 440); 

                auto tg_xzzz_yyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 441); 

                auto tg_xzzz_yyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 442); 

                auto tg_xzzz_yyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 443); 

                auto tg_xzzz_yyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 444); 

                auto tg_xzzz_yyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 445); 

                auto tg_xzzz_yyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 446); 

                auto tg_xzzz_yyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 447); 

                auto tg_xzzz_yzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 448); 

                auto tg_xzzz_zzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 449); 

                auto tg_yyyy_xxxxxxxx_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 450); 

                auto tg_yyyy_xxxxxxxy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 451); 

                auto tg_yyyy_xxxxxxxz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 452); 

                auto tg_yyyy_xxxxxxyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 453); 

                auto tg_yyyy_xxxxxxyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 454); 

                auto tg_yyyy_xxxxxxzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 455); 

                auto tg_yyyy_xxxxxyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 456); 

                auto tg_yyyy_xxxxxyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 457); 

                auto tg_yyyy_xxxxxyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 458); 

                auto tg_yyyy_xxxxxzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 459); 

                auto tg_yyyy_xxxxyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 460); 

                auto tg_yyyy_xxxxyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 461); 

                auto tg_yyyy_xxxxyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 462); 

                auto tg_yyyy_xxxxyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 463); 

                auto tg_yyyy_xxxxzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 464); 

                auto tg_yyyy_xxxyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 465); 

                auto tg_yyyy_xxxyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 466); 

                auto tg_yyyy_xxxyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 467); 

                auto tg_yyyy_xxxyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 468); 

                auto tg_yyyy_xxxyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 469); 

                auto tg_yyyy_xxxzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 470); 

                auto tg_yyyy_xxyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 471); 

                auto tg_yyyy_xxyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 472); 

                auto tg_yyyy_xxyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 473); 

                auto tg_yyyy_xxyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 474); 

                auto tg_yzz_xxxzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 380); 

                auto tg_yzz_xxyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 381); 

                auto tg_yzz_xxyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 382); 

                auto tg_yzz_xxyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 383); 

                auto tg_yzz_xxyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 384); 

                auto tg_yzz_xxyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 385); 

                auto tg_yzz_xxyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 386); 

                auto tg_yzz_xxzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 387); 

                auto tg_yzz_xyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 388); 

                auto tg_yzz_xyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 389); 

                auto tg_yzz_xyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 390); 

                auto tg_yzz_xyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 391); 

                auto tg_yzz_xyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 392); 

                auto tg_yzz_xyyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 393); 

                auto tg_yzz_xyzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 394); 

                auto tg_yzz_xzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 395); 

                auto tg_yzz_yyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 396); 

                auto tg_yzz_yyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 397); 

                auto tg_yzz_yyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 398); 

                auto tg_yzz_yyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 399); 

                auto tg_yzz_yyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 400); 

                auto tg_yzz_yyyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 401); 

                auto tg_yzz_yyzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 402); 

                auto tg_yzz_yzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 403); 

                auto tg_yzz_zzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 404); 

                auto tg_zzz_xxxxxxxx_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 405); 

                auto tg_zzz_xxxxxxxy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 406); 

                auto tg_zzz_xxxxxxxz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 407); 

                auto tg_zzz_xxxxxxyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 408); 

                auto tg_zzz_xxxxxxyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 409); 

                auto tg_zzz_xxxxxxzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 410); 

                auto tg_zzz_xxxxxyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 411); 

                auto tg_zzz_xxxxxyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 412); 

                auto tg_zzz_xxxxxyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 413); 

                auto tg_zzz_xxxxxzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 414); 

                auto tg_zzz_xxxxyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 415); 

                auto tg_zzz_xxxxyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 416); 

                auto tg_zzz_xxxxyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 417); 

                auto tg_zzz_xxxxyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 418); 

                auto tg_zzz_xxxxzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 419); 

                auto tg_zzz_xxxyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 420); 

                auto tg_zzz_xxxyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 421); 

                auto tg_zzz_xxxyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 422); 

                auto tg_zzz_xxxyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 423); 

                auto tg_zzz_xxxyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 424); 

                auto tg_zzz_xxxzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 425); 

                auto tg_zzz_xxyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 426); 

                auto tg_zzz_xxyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 427); 

                auto tg_zzz_xxyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 428); 

                auto tg_zzz_xxyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 429); 

                auto tg_zzz_xxyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 430); 

                auto tg_zzz_xxyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 431); 

                auto tg_zzz_xxzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 432); 

                auto tg_zzz_xyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 433); 

                auto tg_zzz_xyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 434); 

                auto tg_zzz_xyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 435); 

                auto tg_zzz_xyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 436); 

                auto tg_zzz_xyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 437); 

                auto tg_zzz_xyyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 438); 

                auto tg_zzz_xyzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 439); 

                auto tg_zzz_xzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 440); 

                auto tg_zzz_yyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 441); 

                auto tg_zzz_yyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 442); 

                auto tg_zzz_yyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 443); 

                auto tg_zzz_yyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 444); 

                auto tg_zzz_yyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 445); 

                auto tg_zzz_yyyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 446); 

                auto tg_zzz_yyzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 447); 

                auto tg_zzz_yzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 448); 

                auto tg_zzz_zzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 449); 

                auto tg_yzz_xxxzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 380); 

                auto tg_yzz_xxyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 381); 

                auto tg_yzz_xxyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 382); 

                auto tg_yzz_xxyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 383); 

                auto tg_yzz_xxyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 384); 

                auto tg_yzz_xxyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 385); 

                auto tg_yzz_xxyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 386); 

                auto tg_yzz_xxzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 387); 

                auto tg_yzz_xyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 388); 

                auto tg_yzz_xyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 389); 

                auto tg_yzz_xyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 390); 

                auto tg_yzz_xyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 391); 

                auto tg_yzz_xyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 392); 

                auto tg_yzz_xyyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 393); 

                auto tg_yzz_xyzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 394); 

                auto tg_yzz_xzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 395); 

                auto tg_yzz_yyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 396); 

                auto tg_yzz_yyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 397); 

                auto tg_yzz_yyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 398); 

                auto tg_yzz_yyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 399); 

                auto tg_yzz_yyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 400); 

                auto tg_yzz_yyyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 401); 

                auto tg_yzz_yyzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 402); 

                auto tg_yzz_yzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 403); 

                auto tg_yzz_zzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 404); 

                auto tg_zzz_xxxxxxxx_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 405); 

                auto tg_zzz_xxxxxxxy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 406); 

                auto tg_zzz_xxxxxxxz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 407); 

                auto tg_zzz_xxxxxxyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 408); 

                auto tg_zzz_xxxxxxyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 409); 

                auto tg_zzz_xxxxxxzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 410); 

                auto tg_zzz_xxxxxyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 411); 

                auto tg_zzz_xxxxxyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 412); 

                auto tg_zzz_xxxxxyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 413); 

                auto tg_zzz_xxxxxzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 414); 

                auto tg_zzz_xxxxyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 415); 

                auto tg_zzz_xxxxyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 416); 

                auto tg_zzz_xxxxyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 417); 

                auto tg_zzz_xxxxyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 418); 

                auto tg_zzz_xxxxzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 419); 

                auto tg_zzz_xxxyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 420); 

                auto tg_zzz_xxxyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 421); 

                auto tg_zzz_xxxyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 422); 

                auto tg_zzz_xxxyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 423); 

                auto tg_zzz_xxxyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 424); 

                auto tg_zzz_xxxzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 425); 

                auto tg_zzz_xxyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 426); 

                auto tg_zzz_xxyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 427); 

                auto tg_zzz_xxyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 428); 

                auto tg_zzz_xxyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 429); 

                auto tg_zzz_xxyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 430); 

                auto tg_zzz_xxyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 431); 

                auto tg_zzz_xxzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 432); 

                auto tg_zzz_xyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 433); 

                auto tg_zzz_xyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 434); 

                auto tg_zzz_xyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 435); 

                auto tg_zzz_xyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 436); 

                auto tg_zzz_xyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 437); 

                auto tg_zzz_xyyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 438); 

                auto tg_zzz_xyzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 439); 

                auto tg_zzz_xzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 440); 

                auto tg_zzz_yyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 441); 

                auto tg_zzz_yyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 442); 

                auto tg_zzz_yyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 443); 

                auto tg_zzz_yyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 444); 

                auto tg_zzz_yyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 445); 

                auto tg_zzz_yyyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 446); 

                auto tg_zzz_yyzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 447); 

                auto tg_zzz_yzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 448); 

                auto tg_zzz_zzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 449); 

                auto tg_xyzz_xxzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 308); 

                auto tg_xyzz_xyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 309); 

                auto tg_xyzz_xyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 310); 

                auto tg_xyzz_xyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 311); 

                auto tg_xyzz_xyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 312); 

                auto tg_xyzz_xyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 313); 

                auto tg_xyzz_xyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 314); 

                auto tg_xyzz_xzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 315); 

                auto tg_xyzz_yyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 316); 

                auto tg_xyzz_yyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 317); 

                auto tg_xyzz_yyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 318); 

                auto tg_xyzz_yyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 319); 

                auto tg_xyzz_yyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 320); 

                auto tg_xyzz_yyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 321); 

                auto tg_xyzz_yzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 322); 

                auto tg_xyzz_zzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 323); 

                auto tg_xzzz_xxxxxxx_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 324); 

                auto tg_xzzz_xxxxxxy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 325); 

                auto tg_xzzz_xxxxxxz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 326); 

                auto tg_xzzz_xxxxxyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 327); 

                auto tg_xzzz_xxxxxyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 328); 

                auto tg_xzzz_xxxxxzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 329); 

                auto tg_xzzz_xxxxyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 330); 

                auto tg_xzzz_xxxxyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 331); 

                auto tg_xzzz_xxxxyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 332); 

                auto tg_xzzz_xxxxzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 333); 

                auto tg_xzzz_xxxyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 334); 

                auto tg_xzzz_xxxyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 335); 

                auto tg_xzzz_xxxyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 336); 

                auto tg_xzzz_xxxyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 337); 

                auto tg_xzzz_xxxzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 338); 

                auto tg_xzzz_xxyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 339); 

                auto tg_xzzz_xxyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 340); 

                auto tg_xzzz_xxyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 341); 

                auto tg_xzzz_xxyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 342); 

                auto tg_xzzz_xxyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 343); 

                auto tg_xzzz_xxzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 344); 

                auto tg_xzzz_xyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 345); 

                auto tg_xzzz_xyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 346); 

                auto tg_xzzz_xyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 347); 

                auto tg_xzzz_xyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 348); 

                auto tg_xzzz_xyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 349); 

                auto tg_xzzz_xyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 350); 

                auto tg_xzzz_xzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 351); 

                auto tg_xzzz_yyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 352); 

                auto tg_xzzz_yyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 353); 

                auto tg_xzzz_yyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 354); 

                auto tg_xzzz_yyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 355); 

                auto tg_xzzz_yyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 356); 

                auto tg_xzzz_yyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 357); 

                auto tg_xzzz_yzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 358); 

                auto tg_xzzz_zzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 359); 

                auto tg_yyyy_xxxxxxx_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 360); 

                auto tg_yyyy_xxxxxxy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 361); 

                auto tg_yyyy_xxxxxxz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 362); 

                auto tg_yyyy_xxxxxyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 363); 

                auto tg_yyyy_xxxxxyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 364); 

                auto tg_yyyy_xxxxxzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 365); 

                auto tg_yyyy_xxxxyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 366); 

                auto tg_yyyy_xxxxyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 367); 

                auto tg_yyyy_xxxxyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 368); 

                auto tg_yyyy_xxxxzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 369); 

                auto tg_yyyy_xxxyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 370); 

                auto tg_yyyy_xxxyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 371); 

                auto tg_yyyy_xxxyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 372); 

                auto tg_yyyy_xxxyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 373); 

                auto tg_yyyy_xxxzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 374); 

                auto tg_yyyy_xxyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 375); 

                auto tg_yyyy_xxyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 376); 

                auto tg_yyyy_xxyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 377); 

                auto tg_yyyy_xxyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 378); 

                auto tg_yyyy_xxyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 379); 

                auto tg_yyyy_xxzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 380); 

                auto tg_yyyy_xyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 381); 

                auto tg_yyyy_xyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 382); 

                auto tg_yyyy_xyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 383); 

                auto tg_yyyy_xyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 384); 

                // set up pointers to integrals

                auto tg_xxyzz_xxxzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 380); 

                auto tg_xxyzz_xxyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 381); 

                auto tg_xxyzz_xxyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 382); 

                auto tg_xxyzz_xxyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 383); 

                auto tg_xxyzz_xxyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 384); 

                auto tg_xxyzz_xxyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 385); 

                auto tg_xxyzz_xxyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 386); 

                auto tg_xxyzz_xxzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 387); 

                auto tg_xxyzz_xyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 388); 

                auto tg_xxyzz_xyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 389); 

                auto tg_xxyzz_xyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 390); 

                auto tg_xxyzz_xyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 391); 

                auto tg_xxyzz_xyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 392); 

                auto tg_xxyzz_xyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 393); 

                auto tg_xxyzz_xyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 394); 

                auto tg_xxyzz_xzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 395); 

                auto tg_xxyzz_yyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 396); 

                auto tg_xxyzz_yyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 397); 

                auto tg_xxyzz_yyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 398); 

                auto tg_xxyzz_yyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 399); 

                auto tg_xxyzz_yyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 400); 

                auto tg_xxyzz_yyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 401); 

                auto tg_xxyzz_yyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 402); 

                auto tg_xxyzz_yzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 403); 

                auto tg_xxyzz_zzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 404); 

                auto tg_xxzzz_xxxxxxxx_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 405); 

                auto tg_xxzzz_xxxxxxxy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 406); 

                auto tg_xxzzz_xxxxxxxz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 407); 

                auto tg_xxzzz_xxxxxxyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 408); 

                auto tg_xxzzz_xxxxxxyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 409); 

                auto tg_xxzzz_xxxxxxzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 410); 

                auto tg_xxzzz_xxxxxyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 411); 

                auto tg_xxzzz_xxxxxyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 412); 

                auto tg_xxzzz_xxxxxyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 413); 

                auto tg_xxzzz_xxxxxzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 414); 

                auto tg_xxzzz_xxxxyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 415); 

                auto tg_xxzzz_xxxxyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 416); 

                auto tg_xxzzz_xxxxyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 417); 

                auto tg_xxzzz_xxxxyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 418); 

                auto tg_xxzzz_xxxxzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 419); 

                auto tg_xxzzz_xxxyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 420); 

                auto tg_xxzzz_xxxyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 421); 

                auto tg_xxzzz_xxxyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 422); 

                auto tg_xxzzz_xxxyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 423); 

                auto tg_xxzzz_xxxyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 424); 

                auto tg_xxzzz_xxxzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 425); 

                auto tg_xxzzz_xxyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 426); 

                auto tg_xxzzz_xxyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 427); 

                auto tg_xxzzz_xxyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 428); 

                auto tg_xxzzz_xxyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 429); 

                auto tg_xxzzz_xxyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 430); 

                auto tg_xxzzz_xxyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 431); 

                auto tg_xxzzz_xxzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 432); 

                auto tg_xxzzz_xyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 433); 

                auto tg_xxzzz_xyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 434); 

                auto tg_xxzzz_xyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 435); 

                auto tg_xxzzz_xyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 436); 

                auto tg_xxzzz_xyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 437); 

                auto tg_xxzzz_xyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 438); 

                auto tg_xxzzz_xyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 439); 

                auto tg_xxzzz_xzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 440); 

                auto tg_xxzzz_yyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 441); 

                auto tg_xxzzz_yyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 442); 

                auto tg_xxzzz_yyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 443); 

                auto tg_xxzzz_yyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 444); 

                auto tg_xxzzz_yyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 445); 

                auto tg_xxzzz_yyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 446); 

                auto tg_xxzzz_yyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 447); 

                auto tg_xxzzz_yzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 448); 

                auto tg_xxzzz_zzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 449); 

                auto tg_xyyyy_xxxxxxxx_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 450); 

                auto tg_xyyyy_xxxxxxxy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 451); 

                auto tg_xyyyy_xxxxxxxz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 452); 

                auto tg_xyyyy_xxxxxxyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 453); 

                auto tg_xyyyy_xxxxxxyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 454); 

                auto tg_xyyyy_xxxxxxzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 455); 

                auto tg_xyyyy_xxxxxyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 456); 

                auto tg_xyyyy_xxxxxyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 457); 

                auto tg_xyyyy_xxxxxyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 458); 

                auto tg_xyyyy_xxxxxzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 459); 

                auto tg_xyyyy_xxxxyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 460); 

                auto tg_xyyyy_xxxxyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 461); 

                auto tg_xyyyy_xxxxyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 462); 

                auto tg_xyyyy_xxxxyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 463); 

                auto tg_xyyyy_xxxxzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 464); 

                auto tg_xyyyy_xxxyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 465); 

                auto tg_xyyyy_xxxyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 466); 

                auto tg_xyyyy_xxxyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 467); 

                auto tg_xyyyy_xxxyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 468); 

                auto tg_xyyyy_xxxyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 469); 

                auto tg_xyyyy_xxxzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 470); 

                auto tg_xyyyy_xxyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 471); 

                auto tg_xyyyy_xxyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 472); 

                auto tg_xyyyy_xxyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 473); 

                auto tg_xyyyy_xxyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 474); 

                // Batch of Integrals (380,475)

                #pragma omp simd aligned(fxn, fza, tg_xxyzz_xxxzzzzz_0, tg_xxyzz_xxyyyyyy_0, \
                                         tg_xxyzz_xxyyyyyz_0, tg_xxyzz_xxyyyyzz_0, tg_xxyzz_xxyyyzzz_0, tg_xxyzz_xxyyzzzz_0, \
                                         tg_xxyzz_xxyzzzzz_0, tg_xxyzz_xxzzzzzz_0, tg_xxyzz_xyyyyyyy_0, tg_xxyzz_xyyyyyyz_0, \
                                         tg_xxyzz_xyyyyyzz_0, tg_xxyzz_xyyyyzzz_0, tg_xxyzz_xyyyzzzz_0, tg_xxyzz_xyyzzzzz_0, \
                                         tg_xxyzz_xyzzzzzz_0, tg_xxyzz_xzzzzzzz_0, tg_xxyzz_yyyyyyyy_0, tg_xxyzz_yyyyyyyz_0, \
                                         tg_xxyzz_yyyyyyzz_0, tg_xxyzz_yyyyyzzz_0, tg_xxyzz_yyyyzzzz_0, tg_xxyzz_yyyzzzzz_0, \
                                         tg_xxyzz_yyzzzzzz_0, tg_xxyzz_yzzzzzzz_0, tg_xxyzz_zzzzzzzz_0, tg_xxzzz_xxxxxxxx_0, \
                                         tg_xxzzz_xxxxxxxy_0, tg_xxzzz_xxxxxxxz_0, tg_xxzzz_xxxxxxyy_0, tg_xxzzz_xxxxxxyz_0, \
                                         tg_xxzzz_xxxxxxzz_0, tg_xxzzz_xxxxxyyy_0, tg_xxzzz_xxxxxyyz_0, tg_xxzzz_xxxxxyzz_0, \
                                         tg_xxzzz_xxxxxzzz_0, tg_xxzzz_xxxxyyyy_0, tg_xxzzz_xxxxyyyz_0, tg_xxzzz_xxxxyyzz_0, \
                                         tg_xxzzz_xxxxyzzz_0, tg_xxzzz_xxxxzzzz_0, tg_xxzzz_xxxyyyyy_0, tg_xxzzz_xxxyyyyz_0, \
                                         tg_xxzzz_xxxyyyzz_0, tg_xxzzz_xxxyyzzz_0, tg_xxzzz_xxxyzzzz_0, tg_xxzzz_xxxzzzzz_0, \
                                         tg_xxzzz_xxyyyyyy_0, tg_xxzzz_xxyyyyyz_0, tg_xxzzz_xxyyyyzz_0, tg_xxzzz_xxyyyzzz_0, \
                                         tg_xxzzz_xxyyzzzz_0, tg_xxzzz_xxyzzzzz_0, tg_xxzzz_xxzzzzzz_0, tg_xxzzz_xyyyyyyy_0, \
                                         tg_xxzzz_xyyyyyyz_0, tg_xxzzz_xyyyyyzz_0, tg_xxzzz_xyyyyzzz_0, tg_xxzzz_xyyyzzzz_0, \
                                         tg_xxzzz_xyyzzzzz_0, tg_xxzzz_xyzzzzzz_0, tg_xxzzz_xzzzzzzz_0, tg_xxzzz_yyyyyyyy_0, \
                                         tg_xxzzz_yyyyyyyz_0, tg_xxzzz_yyyyyyzz_0, tg_xxzzz_yyyyyzzz_0, tg_xxzzz_yyyyzzzz_0, \
                                         tg_xxzzz_yyyzzzzz_0, tg_xxzzz_yyzzzzzz_0, tg_xxzzz_yzzzzzzz_0, tg_xxzzz_zzzzzzzz_0, \
                                         tg_xyyyy_xxxxxxxx_0, tg_xyyyy_xxxxxxxy_0, tg_xyyyy_xxxxxxxz_0, tg_xyyyy_xxxxxxyy_0, \
                                         tg_xyyyy_xxxxxxyz_0, tg_xyyyy_xxxxxxzz_0, tg_xyyyy_xxxxxyyy_0, tg_xyyyy_xxxxxyyz_0, \
                                         tg_xyyyy_xxxxxyzz_0, tg_xyyyy_xxxxxzzz_0, tg_xyyyy_xxxxyyyy_0, tg_xyyyy_xxxxyyyz_0, \
                                         tg_xyyyy_xxxxyyzz_0, tg_xyyyy_xxxxyzzz_0, tg_xyyyy_xxxxzzzz_0, tg_xyyyy_xxxyyyyy_0, \
                                         tg_xyyyy_xxxyyyyz_0, tg_xyyyy_xxxyyyzz_0, tg_xyyyy_xxxyyzzz_0, tg_xyyyy_xxxyzzzz_0, \
                                         tg_xyyyy_xxxzzzzz_0, tg_xyyyy_xxyyyyyy_0, tg_xyyyy_xxyyyyyz_0, tg_xyyyy_xxyyyyzz_0, \
                                         tg_xyyyy_xxyyyzzz_0, tg_xyzz_xxxzzzzz_0, tg_xyzz_xxxzzzzz_1, tg_xyzz_xxyyyyyy_0, \
                                         tg_xyzz_xxyyyyyy_1, tg_xyzz_xxyyyyyz_0, tg_xyzz_xxyyyyyz_1, tg_xyzz_xxyyyyzz_0, \
                                         tg_xyzz_xxyyyyzz_1, tg_xyzz_xxyyyzzz_0, tg_xyzz_xxyyyzzz_1, tg_xyzz_xxyyzzzz_0, \
                                         tg_xyzz_xxyyzzzz_1, tg_xyzz_xxyzzzzz_0, tg_xyzz_xxyzzzzz_1, tg_xyzz_xxzzzzz_1, \
                                         tg_xyzz_xxzzzzzz_0, tg_xyzz_xxzzzzzz_1, tg_xyzz_xyyyyyy_1, tg_xyzz_xyyyyyyy_0, \
                                         tg_xyzz_xyyyyyyy_1, tg_xyzz_xyyyyyyz_0, tg_xyzz_xyyyyyyz_1, tg_xyzz_xyyyyyz_1, \
                                         tg_xyzz_xyyyyyzz_0, tg_xyzz_xyyyyyzz_1, tg_xyzz_xyyyyzz_1, tg_xyzz_xyyyyzzz_0, \
                                         tg_xyzz_xyyyyzzz_1, tg_xyzz_xyyyzzz_1, tg_xyzz_xyyyzzzz_0, tg_xyzz_xyyyzzzz_1, \
                                         tg_xyzz_xyyzzzz_1, tg_xyzz_xyyzzzzz_0, tg_xyzz_xyyzzzzz_1, tg_xyzz_xyzzzzz_1, \
                                         tg_xyzz_xyzzzzzz_0, tg_xyzz_xyzzzzzz_1, tg_xyzz_xzzzzzz_1, tg_xyzz_xzzzzzzz_0, \
                                         tg_xyzz_xzzzzzzz_1, tg_xyzz_yyyyyyy_1, tg_xyzz_yyyyyyyy_0, tg_xyzz_yyyyyyyy_1, \
                                         tg_xyzz_yyyyyyyz_0, tg_xyzz_yyyyyyyz_1, tg_xyzz_yyyyyyz_1, tg_xyzz_yyyyyyzz_0, \
                                         tg_xyzz_yyyyyyzz_1, tg_xyzz_yyyyyzz_1, tg_xyzz_yyyyyzzz_0, tg_xyzz_yyyyyzzz_1, \
                                         tg_xyzz_yyyyzzz_1, tg_xyzz_yyyyzzzz_0, tg_xyzz_yyyyzzzz_1, tg_xyzz_yyyzzzz_1, \
                                         tg_xyzz_yyyzzzzz_0, tg_xyzz_yyyzzzzz_1, tg_xyzz_yyzzzzz_1, tg_xyzz_yyzzzzzz_0, \
                                         tg_xyzz_yyzzzzzz_1, tg_xyzz_yzzzzzz_1, tg_xyzz_yzzzzzzz_0, tg_xyzz_yzzzzzzz_1, \
                                         tg_xyzz_zzzzzzz_1, tg_xyzz_zzzzzzzz_0, tg_xyzz_zzzzzzzz_1, tg_xzzz_xxxxxxx_1, \
                                         tg_xzzz_xxxxxxxx_0, tg_xzzz_xxxxxxxx_1, tg_xzzz_xxxxxxxy_0, tg_xzzz_xxxxxxxy_1, \
                                         tg_xzzz_xxxxxxxz_0, tg_xzzz_xxxxxxxz_1, tg_xzzz_xxxxxxy_1, tg_xzzz_xxxxxxyy_0, \
                                         tg_xzzz_xxxxxxyy_1, tg_xzzz_xxxxxxyz_0, tg_xzzz_xxxxxxyz_1, tg_xzzz_xxxxxxz_1, \
                                         tg_xzzz_xxxxxxzz_0, tg_xzzz_xxxxxxzz_1, tg_xzzz_xxxxxyy_1, tg_xzzz_xxxxxyyy_0, \
                                         tg_xzzz_xxxxxyyy_1, tg_xzzz_xxxxxyyz_0, tg_xzzz_xxxxxyyz_1, tg_xzzz_xxxxxyz_1, \
                                         tg_xzzz_xxxxxyzz_0, tg_xzzz_xxxxxyzz_1, tg_xzzz_xxxxxzz_1, tg_xzzz_xxxxxzzz_0, \
                                         tg_xzzz_xxxxxzzz_1, tg_xzzz_xxxxyyy_1, tg_xzzz_xxxxyyyy_0, tg_xzzz_xxxxyyyy_1, \
                                         tg_xzzz_xxxxyyyz_0, tg_xzzz_xxxxyyyz_1, tg_xzzz_xxxxyyz_1, tg_xzzz_xxxxyyzz_0, \
                                         tg_xzzz_xxxxyyzz_1, tg_xzzz_xxxxyzz_1, tg_xzzz_xxxxyzzz_0, tg_xzzz_xxxxyzzz_1, \
                                         tg_xzzz_xxxxzzz_1, tg_xzzz_xxxxzzzz_0, tg_xzzz_xxxxzzzz_1, tg_xzzz_xxxyyyy_1, \
                                         tg_xzzz_xxxyyyyy_0, tg_xzzz_xxxyyyyy_1, tg_xzzz_xxxyyyyz_0, tg_xzzz_xxxyyyyz_1, \
                                         tg_xzzz_xxxyyyz_1, tg_xzzz_xxxyyyzz_0, tg_xzzz_xxxyyyzz_1, tg_xzzz_xxxyyzz_1, \
                                         tg_xzzz_xxxyyzzz_0, tg_xzzz_xxxyyzzz_1, tg_xzzz_xxxyzzz_1, tg_xzzz_xxxyzzzz_0, \
                                         tg_xzzz_xxxyzzzz_1, tg_xzzz_xxxzzzz_1, tg_xzzz_xxxzzzzz_0, tg_xzzz_xxxzzzzz_1, \
                                         tg_xzzz_xxyyyyy_1, tg_xzzz_xxyyyyyy_0, tg_xzzz_xxyyyyyy_1, tg_xzzz_xxyyyyyz_0, \
                                         tg_xzzz_xxyyyyyz_1, tg_xzzz_xxyyyyz_1, tg_xzzz_xxyyyyzz_0, tg_xzzz_xxyyyyzz_1, \
                                         tg_xzzz_xxyyyzz_1, tg_xzzz_xxyyyzzz_0, tg_xzzz_xxyyyzzz_1, tg_xzzz_xxyyzzz_1, \
                                         tg_xzzz_xxyyzzzz_0, tg_xzzz_xxyyzzzz_1, tg_xzzz_xxyzzzz_1, tg_xzzz_xxyzzzzz_0, \
                                         tg_xzzz_xxyzzzzz_1, tg_xzzz_xxzzzzz_1, tg_xzzz_xxzzzzzz_0, tg_xzzz_xxzzzzzz_1, \
                                         tg_xzzz_xyyyyyy_1, tg_xzzz_xyyyyyyy_0, tg_xzzz_xyyyyyyy_1, tg_xzzz_xyyyyyyz_0, \
                                         tg_xzzz_xyyyyyyz_1, tg_xzzz_xyyyyyz_1, tg_xzzz_xyyyyyzz_0, tg_xzzz_xyyyyyzz_1, \
                                         tg_xzzz_xyyyyzz_1, tg_xzzz_xyyyyzzz_0, tg_xzzz_xyyyyzzz_1, tg_xzzz_xyyyzzz_1, \
                                         tg_xzzz_xyyyzzzz_0, tg_xzzz_xyyyzzzz_1, tg_xzzz_xyyzzzz_1, tg_xzzz_xyyzzzzz_0, \
                                         tg_xzzz_xyyzzzzz_1, tg_xzzz_xyzzzzz_1, tg_xzzz_xyzzzzzz_0, tg_xzzz_xyzzzzzz_1, \
                                         tg_xzzz_xzzzzzz_1, tg_xzzz_xzzzzzzz_0, tg_xzzz_xzzzzzzz_1, tg_xzzz_yyyyyyy_1, \
                                         tg_xzzz_yyyyyyyy_0, tg_xzzz_yyyyyyyy_1, tg_xzzz_yyyyyyyz_0, tg_xzzz_yyyyyyyz_1, \
                                         tg_xzzz_yyyyyyz_1, tg_xzzz_yyyyyyzz_0, tg_xzzz_yyyyyyzz_1, tg_xzzz_yyyyyzz_1, \
                                         tg_xzzz_yyyyyzzz_0, tg_xzzz_yyyyyzzz_1, tg_xzzz_yyyyzzz_1, tg_xzzz_yyyyzzzz_0, \
                                         tg_xzzz_yyyyzzzz_1, tg_xzzz_yyyzzzz_1, tg_xzzz_yyyzzzzz_0, tg_xzzz_yyyzzzzz_1, \
                                         tg_xzzz_yyzzzzz_1, tg_xzzz_yyzzzzzz_0, tg_xzzz_yyzzzzzz_1, tg_xzzz_yzzzzzz_1, \
                                         tg_xzzz_yzzzzzzz_0, tg_xzzz_yzzzzzzz_1, tg_xzzz_zzzzzzz_1, tg_xzzz_zzzzzzzz_0, \
                                         tg_xzzz_zzzzzzzz_1, tg_yyyy_xxxxxxx_1, tg_yyyy_xxxxxxxx_0, tg_yyyy_xxxxxxxx_1, \
                                         tg_yyyy_xxxxxxxy_0, tg_yyyy_xxxxxxxy_1, tg_yyyy_xxxxxxxz_0, tg_yyyy_xxxxxxxz_1, \
                                         tg_yyyy_xxxxxxy_1, tg_yyyy_xxxxxxyy_0, tg_yyyy_xxxxxxyy_1, tg_yyyy_xxxxxxyz_0, \
                                         tg_yyyy_xxxxxxyz_1, tg_yyyy_xxxxxxz_1, tg_yyyy_xxxxxxzz_0, tg_yyyy_xxxxxxzz_1, \
                                         tg_yyyy_xxxxxyy_1, tg_yyyy_xxxxxyyy_0, tg_yyyy_xxxxxyyy_1, tg_yyyy_xxxxxyyz_0, \
                                         tg_yyyy_xxxxxyyz_1, tg_yyyy_xxxxxyz_1, tg_yyyy_xxxxxyzz_0, tg_yyyy_xxxxxyzz_1, \
                                         tg_yyyy_xxxxxzz_1, tg_yyyy_xxxxxzzz_0, tg_yyyy_xxxxxzzz_1, tg_yyyy_xxxxyyy_1, \
                                         tg_yyyy_xxxxyyyy_0, tg_yyyy_xxxxyyyy_1, tg_yyyy_xxxxyyyz_0, tg_yyyy_xxxxyyyz_1, \
                                         tg_yyyy_xxxxyyz_1, tg_yyyy_xxxxyyzz_0, tg_yyyy_xxxxyyzz_1, tg_yyyy_xxxxyzz_1, \
                                         tg_yyyy_xxxxyzzz_0, tg_yyyy_xxxxyzzz_1, tg_yyyy_xxxxzzz_1, tg_yyyy_xxxxzzzz_0, \
                                         tg_yyyy_xxxxzzzz_1, tg_yyyy_xxxyyyy_1, tg_yyyy_xxxyyyyy_0, tg_yyyy_xxxyyyyy_1, \
                                         tg_yyyy_xxxyyyyz_0, tg_yyyy_xxxyyyyz_1, tg_yyyy_xxxyyyz_1, tg_yyyy_xxxyyyzz_0, \
                                         tg_yyyy_xxxyyyzz_1, tg_yyyy_xxxyyzz_1, tg_yyyy_xxxyyzzz_0, tg_yyyy_xxxyyzzz_1, \
                                         tg_yyyy_xxxyzzz_1, tg_yyyy_xxxyzzzz_0, tg_yyyy_xxxyzzzz_1, tg_yyyy_xxxzzzz_1, \
                                         tg_yyyy_xxxzzzzz_0, tg_yyyy_xxxzzzzz_1, tg_yyyy_xxyyyyy_1, tg_yyyy_xxyyyyyy_0, \
                                         tg_yyyy_xxyyyyyy_1, tg_yyyy_xxyyyyyz_0, tg_yyyy_xxyyyyyz_1, tg_yyyy_xxyyyyz_1, \
                                         tg_yyyy_xxyyyyzz_0, tg_yyyy_xxyyyyzz_1, tg_yyyy_xxyyyzz_1, tg_yyyy_xxyyyzzz_0, \
                                         tg_yyyy_xxyyyzzz_1, tg_yyyy_xxyyzzz_1, tg_yyyy_xxyzzzz_1, tg_yyyy_xxzzzzz_1, \
                                         tg_yyyy_xyyyyyy_1, tg_yyyy_xyyyyyz_1, tg_yyyy_xyyyyzz_1, tg_yyyy_xyyyzzz_1, \
                                         tg_yzz_xxxzzzzz_0, tg_yzz_xxxzzzzz_1, tg_yzz_xxyyyyyy_0, tg_yzz_xxyyyyyy_1, \
                                         tg_yzz_xxyyyyyz_0, tg_yzz_xxyyyyyz_1, tg_yzz_xxyyyyzz_0, tg_yzz_xxyyyyzz_1, \
                                         tg_yzz_xxyyyzzz_0, tg_yzz_xxyyyzzz_1, tg_yzz_xxyyzzzz_0, tg_yzz_xxyyzzzz_1, \
                                         tg_yzz_xxyzzzzz_0, tg_yzz_xxyzzzzz_1, tg_yzz_xxzzzzzz_0, tg_yzz_xxzzzzzz_1, \
                                         tg_yzz_xyyyyyyy_0, tg_yzz_xyyyyyyy_1, tg_yzz_xyyyyyyz_0, tg_yzz_xyyyyyyz_1, \
                                         tg_yzz_xyyyyyzz_0, tg_yzz_xyyyyyzz_1, tg_yzz_xyyyyzzz_0, tg_yzz_xyyyyzzz_1, \
                                         tg_yzz_xyyyzzzz_0, tg_yzz_xyyyzzzz_1, tg_yzz_xyyzzzzz_0, tg_yzz_xyyzzzzz_1, \
                                         tg_yzz_xyzzzzzz_0, tg_yzz_xyzzzzzz_1, tg_yzz_xzzzzzzz_0, tg_yzz_xzzzzzzz_1, \
                                         tg_yzz_yyyyyyyy_0, tg_yzz_yyyyyyyy_1, tg_yzz_yyyyyyyz_0, tg_yzz_yyyyyyyz_1, \
                                         tg_yzz_yyyyyyzz_0, tg_yzz_yyyyyyzz_1, tg_yzz_yyyyyzzz_0, tg_yzz_yyyyyzzz_1, \
                                         tg_yzz_yyyyzzzz_0, tg_yzz_yyyyzzzz_1, tg_yzz_yyyzzzzz_0, tg_yzz_yyyzzzzz_1, \
                                         tg_yzz_yyzzzzzz_0, tg_yzz_yyzzzzzz_1, tg_yzz_yzzzzzzz_0, tg_yzz_yzzzzzzz_1, \
                                         tg_yzz_zzzzzzzz_0, tg_yzz_zzzzzzzz_1, tg_zzz_xxxxxxxx_0, tg_zzz_xxxxxxxx_1, \
                                         tg_zzz_xxxxxxxy_0, tg_zzz_xxxxxxxy_1, tg_zzz_xxxxxxxz_0, tg_zzz_xxxxxxxz_1, \
                                         tg_zzz_xxxxxxyy_0, tg_zzz_xxxxxxyy_1, tg_zzz_xxxxxxyz_0, tg_zzz_xxxxxxyz_1, \
                                         tg_zzz_xxxxxxzz_0, tg_zzz_xxxxxxzz_1, tg_zzz_xxxxxyyy_0, tg_zzz_xxxxxyyy_1, \
                                         tg_zzz_xxxxxyyz_0, tg_zzz_xxxxxyyz_1, tg_zzz_xxxxxyzz_0, tg_zzz_xxxxxyzz_1, \
                                         tg_zzz_xxxxxzzz_0, tg_zzz_xxxxxzzz_1, tg_zzz_xxxxyyyy_0, tg_zzz_xxxxyyyy_1, \
                                         tg_zzz_xxxxyyyz_0, tg_zzz_xxxxyyyz_1, tg_zzz_xxxxyyzz_0, tg_zzz_xxxxyyzz_1, \
                                         tg_zzz_xxxxyzzz_0, tg_zzz_xxxxyzzz_1, tg_zzz_xxxxzzzz_0, tg_zzz_xxxxzzzz_1, \
                                         tg_zzz_xxxyyyyy_0, tg_zzz_xxxyyyyy_1, tg_zzz_xxxyyyyz_0, tg_zzz_xxxyyyyz_1, \
                                         tg_zzz_xxxyyyzz_0, tg_zzz_xxxyyyzz_1, tg_zzz_xxxyyzzz_0, tg_zzz_xxxyyzzz_1, \
                                         tg_zzz_xxxyzzzz_0, tg_zzz_xxxyzzzz_1, tg_zzz_xxxzzzzz_0, tg_zzz_xxxzzzzz_1, \
                                         tg_zzz_xxyyyyyy_0, tg_zzz_xxyyyyyy_1, tg_zzz_xxyyyyyz_0, tg_zzz_xxyyyyyz_1, \
                                         tg_zzz_xxyyyyzz_0, tg_zzz_xxyyyyzz_1, tg_zzz_xxyyyzzz_0, tg_zzz_xxyyyzzz_1, \
                                         tg_zzz_xxyyzzzz_0, tg_zzz_xxyyzzzz_1, tg_zzz_xxyzzzzz_0, tg_zzz_xxyzzzzz_1, \
                                         tg_zzz_xxzzzzzz_0, tg_zzz_xxzzzzzz_1, tg_zzz_xyyyyyyy_0, tg_zzz_xyyyyyyy_1, \
                                         tg_zzz_xyyyyyyz_0, tg_zzz_xyyyyyyz_1, tg_zzz_xyyyyyzz_0, tg_zzz_xyyyyyzz_1, \
                                         tg_zzz_xyyyyzzz_0, tg_zzz_xyyyyzzz_1, tg_zzz_xyyyzzzz_0, tg_zzz_xyyyzzzz_1, \
                                         tg_zzz_xyyzzzzz_0, tg_zzz_xyyzzzzz_1, tg_zzz_xyzzzzzz_0, tg_zzz_xyzzzzzz_1, \
                                         tg_zzz_xzzzzzzz_0, tg_zzz_xzzzzzzz_1, tg_zzz_yyyyyyyy_0, tg_zzz_yyyyyyyy_1, \
                                         tg_zzz_yyyyyyyz_0, tg_zzz_yyyyyyyz_1, tg_zzz_yyyyyyzz_0, tg_zzz_yyyyyyzz_1, \
                                         tg_zzz_yyyyyzzz_0, tg_zzz_yyyyyzzz_1, tg_zzz_yyyyzzzz_0, tg_zzz_yyyyzzzz_1, \
                                         tg_zzz_yyyzzzzz_0, tg_zzz_yyyzzzzz_1, tg_zzz_yyzzzzzz_0, tg_zzz_yyzzzzzz_1, \
                                         tg_zzz_yzzzzzzz_0, tg_zzz_yzzzzzzz_1, tg_zzz_zzzzzzzz_0, tg_zzz_zzzzzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxyzz_xxxzzzzz_0[j] = pb_x * tg_xyzz_xxxzzzzz_0[j] + fr * tg_xyzz_xxxzzzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxzzzzz_0[j] - tg_yzz_xxxzzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzz_xxzzzzz_1[j];

                    tg_xxyzz_xxyyyyyy_0[j] = pb_x * tg_xyzz_xxyyyyyy_0[j] + fr * tg_xyzz_xxyyyyyy_1[j] + 0.5 * fl1_fx * (tg_yzz_xxyyyyyy_0[j] - tg_yzz_xxyyyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyzz_xyyyyyy_1[j];

                    tg_xxyzz_xxyyyyyz_0[j] = pb_x * tg_xyzz_xxyyyyyz_0[j] + fr * tg_xyzz_xxyyyyyz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxyyyyyz_0[j] - tg_yzz_xxyyyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyzz_xyyyyyz_1[j];

                    tg_xxyzz_xxyyyyzz_0[j] = pb_x * tg_xyzz_xxyyyyzz_0[j] + fr * tg_xyzz_xxyyyyzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxyyyyzz_0[j] - tg_yzz_xxyyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyzz_xyyyyzz_1[j];

                    tg_xxyzz_xxyyyzzz_0[j] = pb_x * tg_xyzz_xxyyyzzz_0[j] + fr * tg_xyzz_xxyyyzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxyyyzzz_0[j] - tg_yzz_xxyyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyzz_xyyyzzz_1[j];

                    tg_xxyzz_xxyyzzzz_0[j] = pb_x * tg_xyzz_xxyyzzzz_0[j] + fr * tg_xyzz_xxyyzzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxyyzzzz_0[j] - tg_yzz_xxyyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyzz_xyyzzzz_1[j];

                    tg_xxyzz_xxyzzzzz_0[j] = pb_x * tg_xyzz_xxyzzzzz_0[j] + fr * tg_xyzz_xxyzzzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxyzzzzz_0[j] - tg_yzz_xxyzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyzz_xyzzzzz_1[j];

                    tg_xxyzz_xxzzzzzz_0[j] = pb_x * tg_xyzz_xxzzzzzz_0[j] + fr * tg_xyzz_xxzzzzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxzzzzzz_0[j] - tg_yzz_xxzzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyzz_xzzzzzz_1[j];

                    tg_xxyzz_xyyyyyyy_0[j] = pb_x * tg_xyzz_xyyyyyyy_0[j] + fr * tg_xyzz_xyyyyyyy_1[j] + 0.5 * fl1_fx * (tg_yzz_xyyyyyyy_0[j] - tg_yzz_xyyyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzz_yyyyyyy_1[j];

                    tg_xxyzz_xyyyyyyz_0[j] = pb_x * tg_xyzz_xyyyyyyz_0[j] + fr * tg_xyzz_xyyyyyyz_1[j] + 0.5 * fl1_fx * (tg_yzz_xyyyyyyz_0[j] - tg_yzz_xyyyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzz_yyyyyyz_1[j];

                    tg_xxyzz_xyyyyyzz_0[j] = pb_x * tg_xyzz_xyyyyyzz_0[j] + fr * tg_xyzz_xyyyyyzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xyyyyyzz_0[j] - tg_yzz_xyyyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzz_yyyyyzz_1[j];

                    tg_xxyzz_xyyyyzzz_0[j] = pb_x * tg_xyzz_xyyyyzzz_0[j] + fr * tg_xyzz_xyyyyzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xyyyyzzz_0[j] - tg_yzz_xyyyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzz_yyyyzzz_1[j];

                    tg_xxyzz_xyyyzzzz_0[j] = pb_x * tg_xyzz_xyyyzzzz_0[j] + fr * tg_xyzz_xyyyzzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xyyyzzzz_0[j] - tg_yzz_xyyyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzz_yyyzzzz_1[j];

                    tg_xxyzz_xyyzzzzz_0[j] = pb_x * tg_xyzz_xyyzzzzz_0[j] + fr * tg_xyzz_xyyzzzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xyyzzzzz_0[j] - tg_yzz_xyyzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzz_yyzzzzz_1[j];

                    tg_xxyzz_xyzzzzzz_0[j] = pb_x * tg_xyzz_xyzzzzzz_0[j] + fr * tg_xyzz_xyzzzzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xyzzzzzz_0[j] - tg_yzz_xyzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzz_yzzzzzz_1[j];

                    tg_xxyzz_xzzzzzzz_0[j] = pb_x * tg_xyzz_xzzzzzzz_0[j] + fr * tg_xyzz_xzzzzzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xzzzzzzz_0[j] - tg_yzz_xzzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzz_zzzzzzz_1[j];

                    tg_xxyzz_yyyyyyyy_0[j] = pb_x * tg_xyzz_yyyyyyyy_0[j] + fr * tg_xyzz_yyyyyyyy_1[j] + 0.5 * fl1_fx * (tg_yzz_yyyyyyyy_0[j] - tg_yzz_yyyyyyyy_1[j] * fl1_fza);

                    tg_xxyzz_yyyyyyyz_0[j] = pb_x * tg_xyzz_yyyyyyyz_0[j] + fr * tg_xyzz_yyyyyyyz_1[j] + 0.5 * fl1_fx * (tg_yzz_yyyyyyyz_0[j] - tg_yzz_yyyyyyyz_1[j] * fl1_fza);

                    tg_xxyzz_yyyyyyzz_0[j] = pb_x * tg_xyzz_yyyyyyzz_0[j] + fr * tg_xyzz_yyyyyyzz_1[j] + 0.5 * fl1_fx * (tg_yzz_yyyyyyzz_0[j] - tg_yzz_yyyyyyzz_1[j] * fl1_fza);

                    tg_xxyzz_yyyyyzzz_0[j] = pb_x * tg_xyzz_yyyyyzzz_0[j] + fr * tg_xyzz_yyyyyzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_yyyyyzzz_0[j] - tg_yzz_yyyyyzzz_1[j] * fl1_fza);

                    tg_xxyzz_yyyyzzzz_0[j] = pb_x * tg_xyzz_yyyyzzzz_0[j] + fr * tg_xyzz_yyyyzzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_yyyyzzzz_0[j] - tg_yzz_yyyyzzzz_1[j] * fl1_fza);

                    tg_xxyzz_yyyzzzzz_0[j] = pb_x * tg_xyzz_yyyzzzzz_0[j] + fr * tg_xyzz_yyyzzzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_yyyzzzzz_0[j] - tg_yzz_yyyzzzzz_1[j] * fl1_fza);

                    tg_xxyzz_yyzzzzzz_0[j] = pb_x * tg_xyzz_yyzzzzzz_0[j] + fr * tg_xyzz_yyzzzzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_yyzzzzzz_0[j] - tg_yzz_yyzzzzzz_1[j] * fl1_fza);

                    tg_xxyzz_yzzzzzzz_0[j] = pb_x * tg_xyzz_yzzzzzzz_0[j] + fr * tg_xyzz_yzzzzzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_yzzzzzzz_0[j] - tg_yzz_yzzzzzzz_1[j] * fl1_fza);

                    tg_xxyzz_zzzzzzzz_0[j] = pb_x * tg_xyzz_zzzzzzzz_0[j] + fr * tg_xyzz_zzzzzzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_zzzzzzzz_0[j] - tg_yzz_zzzzzzzz_1[j] * fl1_fza);

                    tg_xxzzz_xxxxxxxx_0[j] = pb_x * tg_xzzz_xxxxxxxx_0[j] + fr * tg_xzzz_xxxxxxxx_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxxxxx_0[j] - tg_zzz_xxxxxxxx_1[j] * fl1_fza) + 4.0 * fl1_fxn * tg_xzzz_xxxxxxx_1[j];

                    tg_xxzzz_xxxxxxxy_0[j] = pb_x * tg_xzzz_xxxxxxxy_0[j] + fr * tg_xzzz_xxxxxxxy_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxxxxy_0[j] - tg_zzz_xxxxxxxy_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_xzzz_xxxxxxy_1[j];

                    tg_xxzzz_xxxxxxxz_0[j] = pb_x * tg_xzzz_xxxxxxxz_0[j] + fr * tg_xzzz_xxxxxxxz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxxxxz_0[j] - tg_zzz_xxxxxxxz_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_xzzz_xxxxxxz_1[j];

                    tg_xxzzz_xxxxxxyy_0[j] = pb_x * tg_xzzz_xxxxxxyy_0[j] + fr * tg_xzzz_xxxxxxyy_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxxxyy_0[j] - tg_zzz_xxxxxxyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xzzz_xxxxxyy_1[j];

                    tg_xxzzz_xxxxxxyz_0[j] = pb_x * tg_xzzz_xxxxxxyz_0[j] + fr * tg_xzzz_xxxxxxyz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxxxyz_0[j] - tg_zzz_xxxxxxyz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xzzz_xxxxxyz_1[j];

                    tg_xxzzz_xxxxxxzz_0[j] = pb_x * tg_xzzz_xxxxxxzz_0[j] + fr * tg_xzzz_xxxxxxzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxxxzz_0[j] - tg_zzz_xxxxxxzz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xzzz_xxxxxzz_1[j];

                    tg_xxzzz_xxxxxyyy_0[j] = pb_x * tg_xzzz_xxxxxyyy_0[j] + fr * tg_xzzz_xxxxxyyy_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxxyyy_0[j] - tg_zzz_xxxxxyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xzzz_xxxxyyy_1[j];

                    tg_xxzzz_xxxxxyyz_0[j] = pb_x * tg_xzzz_xxxxxyyz_0[j] + fr * tg_xzzz_xxxxxyyz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxxyyz_0[j] - tg_zzz_xxxxxyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xzzz_xxxxyyz_1[j];

                    tg_xxzzz_xxxxxyzz_0[j] = pb_x * tg_xzzz_xxxxxyzz_0[j] + fr * tg_xzzz_xxxxxyzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxxyzz_0[j] - tg_zzz_xxxxxyzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xzzz_xxxxyzz_1[j];

                    tg_xxzzz_xxxxxzzz_0[j] = pb_x * tg_xzzz_xxxxxzzz_0[j] + fr * tg_xzzz_xxxxxzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxxzzz_0[j] - tg_zzz_xxxxxzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xzzz_xxxxzzz_1[j];

                    tg_xxzzz_xxxxyyyy_0[j] = pb_x * tg_xzzz_xxxxyyyy_0[j] + fr * tg_xzzz_xxxxyyyy_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxyyyy_0[j] - tg_zzz_xxxxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzzz_xxxyyyy_1[j];

                    tg_xxzzz_xxxxyyyz_0[j] = pb_x * tg_xzzz_xxxxyyyz_0[j] + fr * tg_xzzz_xxxxyyyz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxyyyz_0[j] - tg_zzz_xxxxyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzzz_xxxyyyz_1[j];

                    tg_xxzzz_xxxxyyzz_0[j] = pb_x * tg_xzzz_xxxxyyzz_0[j] + fr * tg_xzzz_xxxxyyzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxyyzz_0[j] - tg_zzz_xxxxyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzzz_xxxyyzz_1[j];

                    tg_xxzzz_xxxxyzzz_0[j] = pb_x * tg_xzzz_xxxxyzzz_0[j] + fr * tg_xzzz_xxxxyzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxyzzz_0[j] - tg_zzz_xxxxyzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzzz_xxxyzzz_1[j];

                    tg_xxzzz_xxxxzzzz_0[j] = pb_x * tg_xzzz_xxxxzzzz_0[j] + fr * tg_xzzz_xxxxzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxzzzz_0[j] - tg_zzz_xxxxzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzzz_xxxzzzz_1[j];

                    tg_xxzzz_xxxyyyyy_0[j] = pb_x * tg_xzzz_xxxyyyyy_0[j] + fr * tg_xzzz_xxxyyyyy_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxyyyyy_0[j] - tg_zzz_xxxyyyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzz_xxyyyyy_1[j];

                    tg_xxzzz_xxxyyyyz_0[j] = pb_x * tg_xzzz_xxxyyyyz_0[j] + fr * tg_xzzz_xxxyyyyz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxyyyyz_0[j] - tg_zzz_xxxyyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzz_xxyyyyz_1[j];

                    tg_xxzzz_xxxyyyzz_0[j] = pb_x * tg_xzzz_xxxyyyzz_0[j] + fr * tg_xzzz_xxxyyyzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxyyyzz_0[j] - tg_zzz_xxxyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzz_xxyyyzz_1[j];

                    tg_xxzzz_xxxyyzzz_0[j] = pb_x * tg_xzzz_xxxyyzzz_0[j] + fr * tg_xzzz_xxxyyzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxyyzzz_0[j] - tg_zzz_xxxyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzz_xxyyzzz_1[j];

                    tg_xxzzz_xxxyzzzz_0[j] = pb_x * tg_xzzz_xxxyzzzz_0[j] + fr * tg_xzzz_xxxyzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxyzzzz_0[j] - tg_zzz_xxxyzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzz_xxyzzzz_1[j];

                    tg_xxzzz_xxxzzzzz_0[j] = pb_x * tg_xzzz_xxxzzzzz_0[j] + fr * tg_xzzz_xxxzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxzzzzz_0[j] - tg_zzz_xxxzzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzz_xxzzzzz_1[j];

                    tg_xxzzz_xxyyyyyy_0[j] = pb_x * tg_xzzz_xxyyyyyy_0[j] + fr * tg_xzzz_xxyyyyyy_1[j] + 0.5 * fl1_fx * (tg_zzz_xxyyyyyy_0[j] - tg_zzz_xxyyyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xzzz_xyyyyyy_1[j];

                    tg_xxzzz_xxyyyyyz_0[j] = pb_x * tg_xzzz_xxyyyyyz_0[j] + fr * tg_xzzz_xxyyyyyz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxyyyyyz_0[j] - tg_zzz_xxyyyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xzzz_xyyyyyz_1[j];

                    tg_xxzzz_xxyyyyzz_0[j] = pb_x * tg_xzzz_xxyyyyzz_0[j] + fr * tg_xzzz_xxyyyyzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxyyyyzz_0[j] - tg_zzz_xxyyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xzzz_xyyyyzz_1[j];

                    tg_xxzzz_xxyyyzzz_0[j] = pb_x * tg_xzzz_xxyyyzzz_0[j] + fr * tg_xzzz_xxyyyzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxyyyzzz_0[j] - tg_zzz_xxyyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xzzz_xyyyzzz_1[j];

                    tg_xxzzz_xxyyzzzz_0[j] = pb_x * tg_xzzz_xxyyzzzz_0[j] + fr * tg_xzzz_xxyyzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxyyzzzz_0[j] - tg_zzz_xxyyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xzzz_xyyzzzz_1[j];

                    tg_xxzzz_xxyzzzzz_0[j] = pb_x * tg_xzzz_xxyzzzzz_0[j] + fr * tg_xzzz_xxyzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxyzzzzz_0[j] - tg_zzz_xxyzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xzzz_xyzzzzz_1[j];

                    tg_xxzzz_xxzzzzzz_0[j] = pb_x * tg_xzzz_xxzzzzzz_0[j] + fr * tg_xzzz_xxzzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxzzzzzz_0[j] - tg_zzz_xxzzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xzzz_xzzzzzz_1[j];

                    tg_xxzzz_xyyyyyyy_0[j] = pb_x * tg_xzzz_xyyyyyyy_0[j] + fr * tg_xzzz_xyyyyyyy_1[j] + 0.5 * fl1_fx * (tg_zzz_xyyyyyyy_0[j] - tg_zzz_xyyyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzz_yyyyyyy_1[j];

                    tg_xxzzz_xyyyyyyz_0[j] = pb_x * tg_xzzz_xyyyyyyz_0[j] + fr * tg_xzzz_xyyyyyyz_1[j] + 0.5 * fl1_fx * (tg_zzz_xyyyyyyz_0[j] - tg_zzz_xyyyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzz_yyyyyyz_1[j];

                    tg_xxzzz_xyyyyyzz_0[j] = pb_x * tg_xzzz_xyyyyyzz_0[j] + fr * tg_xzzz_xyyyyyzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xyyyyyzz_0[j] - tg_zzz_xyyyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzz_yyyyyzz_1[j];

                    tg_xxzzz_xyyyyzzz_0[j] = pb_x * tg_xzzz_xyyyyzzz_0[j] + fr * tg_xzzz_xyyyyzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xyyyyzzz_0[j] - tg_zzz_xyyyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzz_yyyyzzz_1[j];

                    tg_xxzzz_xyyyzzzz_0[j] = pb_x * tg_xzzz_xyyyzzzz_0[j] + fr * tg_xzzz_xyyyzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xyyyzzzz_0[j] - tg_zzz_xyyyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzz_yyyzzzz_1[j];

                    tg_xxzzz_xyyzzzzz_0[j] = pb_x * tg_xzzz_xyyzzzzz_0[j] + fr * tg_xzzz_xyyzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xyyzzzzz_0[j] - tg_zzz_xyyzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzz_yyzzzzz_1[j];

                    tg_xxzzz_xyzzzzzz_0[j] = pb_x * tg_xzzz_xyzzzzzz_0[j] + fr * tg_xzzz_xyzzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xyzzzzzz_0[j] - tg_zzz_xyzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzz_yzzzzzz_1[j];

                    tg_xxzzz_xzzzzzzz_0[j] = pb_x * tg_xzzz_xzzzzzzz_0[j] + fr * tg_xzzz_xzzzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xzzzzzzz_0[j] - tg_zzz_xzzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzz_zzzzzzz_1[j];

                    tg_xxzzz_yyyyyyyy_0[j] = pb_x * tg_xzzz_yyyyyyyy_0[j] + fr * tg_xzzz_yyyyyyyy_1[j] + 0.5 * fl1_fx * (tg_zzz_yyyyyyyy_0[j] - tg_zzz_yyyyyyyy_1[j] * fl1_fza);

                    tg_xxzzz_yyyyyyyz_0[j] = pb_x * tg_xzzz_yyyyyyyz_0[j] + fr * tg_xzzz_yyyyyyyz_1[j] + 0.5 * fl1_fx * (tg_zzz_yyyyyyyz_0[j] - tg_zzz_yyyyyyyz_1[j] * fl1_fza);

                    tg_xxzzz_yyyyyyzz_0[j] = pb_x * tg_xzzz_yyyyyyzz_0[j] + fr * tg_xzzz_yyyyyyzz_1[j] + 0.5 * fl1_fx * (tg_zzz_yyyyyyzz_0[j] - tg_zzz_yyyyyyzz_1[j] * fl1_fza);

                    tg_xxzzz_yyyyyzzz_0[j] = pb_x * tg_xzzz_yyyyyzzz_0[j] + fr * tg_xzzz_yyyyyzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_yyyyyzzz_0[j] - tg_zzz_yyyyyzzz_1[j] * fl1_fza);

                    tg_xxzzz_yyyyzzzz_0[j] = pb_x * tg_xzzz_yyyyzzzz_0[j] + fr * tg_xzzz_yyyyzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_yyyyzzzz_0[j] - tg_zzz_yyyyzzzz_1[j] * fl1_fza);

                    tg_xxzzz_yyyzzzzz_0[j] = pb_x * tg_xzzz_yyyzzzzz_0[j] + fr * tg_xzzz_yyyzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_yyyzzzzz_0[j] - tg_zzz_yyyzzzzz_1[j] * fl1_fza);

                    tg_xxzzz_yyzzzzzz_0[j] = pb_x * tg_xzzz_yyzzzzzz_0[j] + fr * tg_xzzz_yyzzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_yyzzzzzz_0[j] - tg_zzz_yyzzzzzz_1[j] * fl1_fza);

                    tg_xxzzz_yzzzzzzz_0[j] = pb_x * tg_xzzz_yzzzzzzz_0[j] + fr * tg_xzzz_yzzzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_yzzzzzzz_0[j] - tg_zzz_yzzzzzzz_1[j] * fl1_fza);

                    tg_xxzzz_zzzzzzzz_0[j] = pb_x * tg_xzzz_zzzzzzzz_0[j] + fr * tg_xzzz_zzzzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_zzzzzzzz_0[j] - tg_zzz_zzzzzzzz_1[j] * fl1_fza);

                    tg_xyyyy_xxxxxxxx_0[j] = pb_x * tg_yyyy_xxxxxxxx_0[j] + fr * tg_yyyy_xxxxxxxx_1[j] + 4.0 * fl1_fxn * tg_yyyy_xxxxxxx_1[j];

                    tg_xyyyy_xxxxxxxy_0[j] = pb_x * tg_yyyy_xxxxxxxy_0[j] + fr * tg_yyyy_xxxxxxxy_1[j] + 3.5 * fl1_fxn * tg_yyyy_xxxxxxy_1[j];

                    tg_xyyyy_xxxxxxxz_0[j] = pb_x * tg_yyyy_xxxxxxxz_0[j] + fr * tg_yyyy_xxxxxxxz_1[j] + 3.5 * fl1_fxn * tg_yyyy_xxxxxxz_1[j];

                    tg_xyyyy_xxxxxxyy_0[j] = pb_x * tg_yyyy_xxxxxxyy_0[j] + fr * tg_yyyy_xxxxxxyy_1[j] + 3.0 * fl1_fxn * tg_yyyy_xxxxxyy_1[j];

                    tg_xyyyy_xxxxxxyz_0[j] = pb_x * tg_yyyy_xxxxxxyz_0[j] + fr * tg_yyyy_xxxxxxyz_1[j] + 3.0 * fl1_fxn * tg_yyyy_xxxxxyz_1[j];

                    tg_xyyyy_xxxxxxzz_0[j] = pb_x * tg_yyyy_xxxxxxzz_0[j] + fr * tg_yyyy_xxxxxxzz_1[j] + 3.0 * fl1_fxn * tg_yyyy_xxxxxzz_1[j];

                    tg_xyyyy_xxxxxyyy_0[j] = pb_x * tg_yyyy_xxxxxyyy_0[j] + fr * tg_yyyy_xxxxxyyy_1[j] + 2.5 * fl1_fxn * tg_yyyy_xxxxyyy_1[j];

                    tg_xyyyy_xxxxxyyz_0[j] = pb_x * tg_yyyy_xxxxxyyz_0[j] + fr * tg_yyyy_xxxxxyyz_1[j] + 2.5 * fl1_fxn * tg_yyyy_xxxxyyz_1[j];

                    tg_xyyyy_xxxxxyzz_0[j] = pb_x * tg_yyyy_xxxxxyzz_0[j] + fr * tg_yyyy_xxxxxyzz_1[j] + 2.5 * fl1_fxn * tg_yyyy_xxxxyzz_1[j];

                    tg_xyyyy_xxxxxzzz_0[j] = pb_x * tg_yyyy_xxxxxzzz_0[j] + fr * tg_yyyy_xxxxxzzz_1[j] + 2.5 * fl1_fxn * tg_yyyy_xxxxzzz_1[j];

                    tg_xyyyy_xxxxyyyy_0[j] = pb_x * tg_yyyy_xxxxyyyy_0[j] + fr * tg_yyyy_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_yyyy_xxxyyyy_1[j];

                    tg_xyyyy_xxxxyyyz_0[j] = pb_x * tg_yyyy_xxxxyyyz_0[j] + fr * tg_yyyy_xxxxyyyz_1[j] + 2.0 * fl1_fxn * tg_yyyy_xxxyyyz_1[j];

                    tg_xyyyy_xxxxyyzz_0[j] = pb_x * tg_yyyy_xxxxyyzz_0[j] + fr * tg_yyyy_xxxxyyzz_1[j] + 2.0 * fl1_fxn * tg_yyyy_xxxyyzz_1[j];

                    tg_xyyyy_xxxxyzzz_0[j] = pb_x * tg_yyyy_xxxxyzzz_0[j] + fr * tg_yyyy_xxxxyzzz_1[j] + 2.0 * fl1_fxn * tg_yyyy_xxxyzzz_1[j];

                    tg_xyyyy_xxxxzzzz_0[j] = pb_x * tg_yyyy_xxxxzzzz_0[j] + fr * tg_yyyy_xxxxzzzz_1[j] + 2.0 * fl1_fxn * tg_yyyy_xxxzzzz_1[j];

                    tg_xyyyy_xxxyyyyy_0[j] = pb_x * tg_yyyy_xxxyyyyy_0[j] + fr * tg_yyyy_xxxyyyyy_1[j] + 1.5 * fl1_fxn * tg_yyyy_xxyyyyy_1[j];

                    tg_xyyyy_xxxyyyyz_0[j] = pb_x * tg_yyyy_xxxyyyyz_0[j] + fr * tg_yyyy_xxxyyyyz_1[j] + 1.5 * fl1_fxn * tg_yyyy_xxyyyyz_1[j];

                    tg_xyyyy_xxxyyyzz_0[j] = pb_x * tg_yyyy_xxxyyyzz_0[j] + fr * tg_yyyy_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_yyyy_xxyyyzz_1[j];

                    tg_xyyyy_xxxyyzzz_0[j] = pb_x * tg_yyyy_xxxyyzzz_0[j] + fr * tg_yyyy_xxxyyzzz_1[j] + 1.5 * fl1_fxn * tg_yyyy_xxyyzzz_1[j];

                    tg_xyyyy_xxxyzzzz_0[j] = pb_x * tg_yyyy_xxxyzzzz_0[j] + fr * tg_yyyy_xxxyzzzz_1[j] + 1.5 * fl1_fxn * tg_yyyy_xxyzzzz_1[j];

                    tg_xyyyy_xxxzzzzz_0[j] = pb_x * tg_yyyy_xxxzzzzz_0[j] + fr * tg_yyyy_xxxzzzzz_1[j] + 1.5 * fl1_fxn * tg_yyyy_xxzzzzz_1[j];

                    tg_xyyyy_xxyyyyyy_0[j] = pb_x * tg_yyyy_xxyyyyyy_0[j] + fr * tg_yyyy_xxyyyyyy_1[j] + fl1_fxn * tg_yyyy_xyyyyyy_1[j];

                    tg_xyyyy_xxyyyyyz_0[j] = pb_x * tg_yyyy_xxyyyyyz_0[j] + fr * tg_yyyy_xxyyyyyz_1[j] + fl1_fxn * tg_yyyy_xyyyyyz_1[j];

                    tg_xyyyy_xxyyyyzz_0[j] = pb_x * tg_yyyy_xxyyyyzz_0[j] + fr * tg_yyyy_xxyyyyzz_1[j] + fl1_fxn * tg_yyyy_xyyyyzz_1[j];

                    tg_xyyyy_xxyyyzzz_0[j] = pb_x * tg_yyyy_xxyyyzzz_0[j] + fr * tg_yyyy_xxyyyzzz_1[j] + fl1_fxn * tg_yyyy_xyyyzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSL_475_569(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (475,569)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {5, -1, -1, -1},
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_8_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_8_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_8_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_7_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {7, -1, -1, -1}, 
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

                auto tg_yyyy_xxyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 475); 

                auto tg_yyyy_xxyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 476); 

                auto tg_yyyy_xxzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 477); 

                auto tg_yyyy_xyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 478); 

                auto tg_yyyy_xyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 479); 

                auto tg_yyyy_xyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 480); 

                auto tg_yyyy_xyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 481); 

                auto tg_yyyy_xyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 482); 

                auto tg_yyyy_xyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 483); 

                auto tg_yyyy_xyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 484); 

                auto tg_yyyy_xzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 485); 

                auto tg_yyyy_yyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 486); 

                auto tg_yyyy_yyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 487); 

                auto tg_yyyy_yyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 488); 

                auto tg_yyyy_yyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 489); 

                auto tg_yyyy_yyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 490); 

                auto tg_yyyy_yyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 491); 

                auto tg_yyyy_yyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 492); 

                auto tg_yyyy_yzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 493); 

                auto tg_yyyy_zzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 494); 

                auto tg_yyyz_xxxxxxxx_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 495); 

                auto tg_yyyz_xxxxxxxy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 496); 

                auto tg_yyyz_xxxxxxxz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 497); 

                auto tg_yyyz_xxxxxxyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 498); 

                auto tg_yyyz_xxxxxxyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 499); 

                auto tg_yyyz_xxxxxxzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 500); 

                auto tg_yyyz_xxxxxyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 501); 

                auto tg_yyyz_xxxxxyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 502); 

                auto tg_yyyz_xxxxxyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 503); 

                auto tg_yyyz_xxxxxzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 504); 

                auto tg_yyyz_xxxxyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 505); 

                auto tg_yyyz_xxxxyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 506); 

                auto tg_yyyz_xxxxyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 507); 

                auto tg_yyyz_xxxxyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 508); 

                auto tg_yyyz_xxxxzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 509); 

                auto tg_yyyz_xxxyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 510); 

                auto tg_yyyz_xxxyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 511); 

                auto tg_yyyz_xxxyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 512); 

                auto tg_yyyz_xxxyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 513); 

                auto tg_yyyz_xxxyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 514); 

                auto tg_yyyz_xxxzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 515); 

                auto tg_yyyz_xxyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 516); 

                auto tg_yyyz_xxyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 517); 

                auto tg_yyyz_xxyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 518); 

                auto tg_yyyz_xxyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 519); 

                auto tg_yyyz_xxyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 520); 

                auto tg_yyyz_xxyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 521); 

                auto tg_yyyz_xxzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 522); 

                auto tg_yyyz_xyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 523); 

                auto tg_yyyz_xyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 524); 

                auto tg_yyyz_xyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 525); 

                auto tg_yyyz_xyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 526); 

                auto tg_yyyz_xyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 527); 

                auto tg_yyyz_xyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 528); 

                auto tg_yyyz_xyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 529); 

                auto tg_yyyz_xzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 530); 

                auto tg_yyyz_yyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 531); 

                auto tg_yyyz_yyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 532); 

                auto tg_yyyz_yyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 533); 

                auto tg_yyyz_yyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 534); 

                auto tg_yyyz_yyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 535); 

                auto tg_yyyz_yyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 536); 

                auto tg_yyyz_yyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 537); 

                auto tg_yyyz_yzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 538); 

                auto tg_yyyz_zzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 539); 

                auto tg_yyzz_xxxxxxxx_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 540); 

                auto tg_yyzz_xxxxxxxy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 541); 

                auto tg_yyzz_xxxxxxxz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 542); 

                auto tg_yyzz_xxxxxxyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 543); 

                auto tg_yyzz_xxxxxxyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 544); 

                auto tg_yyzz_xxxxxxzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 545); 

                auto tg_yyzz_xxxxxyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 546); 

                auto tg_yyzz_xxxxxyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 547); 

                auto tg_yyzz_xxxxxyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 548); 

                auto tg_yyzz_xxxxxzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 549); 

                auto tg_yyzz_xxxxyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 550); 

                auto tg_yyzz_xxxxyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 551); 

                auto tg_yyzz_xxxxyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 552); 

                auto tg_yyzz_xxxxyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 553); 

                auto tg_yyzz_xxxxzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 554); 

                auto tg_yyzz_xxxyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 555); 

                auto tg_yyzz_xxxyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 556); 

                auto tg_yyzz_xxxyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 557); 

                auto tg_yyzz_xxxyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 558); 

                auto tg_yyzz_xxxyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 559); 

                auto tg_yyzz_xxxzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 560); 

                auto tg_yyzz_xxyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 561); 

                auto tg_yyzz_xxyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 562); 

                auto tg_yyzz_xxyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 563); 

                auto tg_yyzz_xxyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 564); 

                auto tg_yyzz_xxyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 565); 

                auto tg_yyzz_xxyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 566); 

                auto tg_yyzz_xxzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 567); 

                auto tg_yyzz_xyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 568); 

                auto tg_yyyy_xxyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 475); 

                auto tg_yyyy_xxyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 476); 

                auto tg_yyyy_xxzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 477); 

                auto tg_yyyy_xyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 478); 

                auto tg_yyyy_xyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 479); 

                auto tg_yyyy_xyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 480); 

                auto tg_yyyy_xyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 481); 

                auto tg_yyyy_xyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 482); 

                auto tg_yyyy_xyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 483); 

                auto tg_yyyy_xyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 484); 

                auto tg_yyyy_xzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 485); 

                auto tg_yyyy_yyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 486); 

                auto tg_yyyy_yyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 487); 

                auto tg_yyyy_yyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 488); 

                auto tg_yyyy_yyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 489); 

                auto tg_yyyy_yyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 490); 

                auto tg_yyyy_yyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 491); 

                auto tg_yyyy_yyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 492); 

                auto tg_yyyy_yzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 493); 

                auto tg_yyyy_zzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 494); 

                auto tg_yyyz_xxxxxxxx_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 495); 

                auto tg_yyyz_xxxxxxxy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 496); 

                auto tg_yyyz_xxxxxxxz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 497); 

                auto tg_yyyz_xxxxxxyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 498); 

                auto tg_yyyz_xxxxxxyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 499); 

                auto tg_yyyz_xxxxxxzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 500); 

                auto tg_yyyz_xxxxxyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 501); 

                auto tg_yyyz_xxxxxyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 502); 

                auto tg_yyyz_xxxxxyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 503); 

                auto tg_yyyz_xxxxxzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 504); 

                auto tg_yyyz_xxxxyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 505); 

                auto tg_yyyz_xxxxyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 506); 

                auto tg_yyyz_xxxxyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 507); 

                auto tg_yyyz_xxxxyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 508); 

                auto tg_yyyz_xxxxzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 509); 

                auto tg_yyyz_xxxyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 510); 

                auto tg_yyyz_xxxyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 511); 

                auto tg_yyyz_xxxyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 512); 

                auto tg_yyyz_xxxyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 513); 

                auto tg_yyyz_xxxyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 514); 

                auto tg_yyyz_xxxzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 515); 

                auto tg_yyyz_xxyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 516); 

                auto tg_yyyz_xxyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 517); 

                auto tg_yyyz_xxyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 518); 

                auto tg_yyyz_xxyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 519); 

                auto tg_yyyz_xxyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 520); 

                auto tg_yyyz_xxyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 521); 

                auto tg_yyyz_xxzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 522); 

                auto tg_yyyz_xyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 523); 

                auto tg_yyyz_xyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 524); 

                auto tg_yyyz_xyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 525); 

                auto tg_yyyz_xyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 526); 

                auto tg_yyyz_xyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 527); 

                auto tg_yyyz_xyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 528); 

                auto tg_yyyz_xyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 529); 

                auto tg_yyyz_xzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 530); 

                auto tg_yyyz_yyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 531); 

                auto tg_yyyz_yyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 532); 

                auto tg_yyyz_yyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 533); 

                auto tg_yyyz_yyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 534); 

                auto tg_yyyz_yyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 535); 

                auto tg_yyyz_yyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 536); 

                auto tg_yyyz_yyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 537); 

                auto tg_yyyz_yzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 538); 

                auto tg_yyyz_zzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 539); 

                auto tg_yyzz_xxxxxxxx_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 540); 

                auto tg_yyzz_xxxxxxxy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 541); 

                auto tg_yyzz_xxxxxxxz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 542); 

                auto tg_yyzz_xxxxxxyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 543); 

                auto tg_yyzz_xxxxxxyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 544); 

                auto tg_yyzz_xxxxxxzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 545); 

                auto tg_yyzz_xxxxxyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 546); 

                auto tg_yyzz_xxxxxyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 547); 

                auto tg_yyzz_xxxxxyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 548); 

                auto tg_yyzz_xxxxxzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 549); 

                auto tg_yyzz_xxxxyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 550); 

                auto tg_yyzz_xxxxyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 551); 

                auto tg_yyzz_xxxxyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 552); 

                auto tg_yyzz_xxxxyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 553); 

                auto tg_yyzz_xxxxzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 554); 

                auto tg_yyzz_xxxyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 555); 

                auto tg_yyzz_xxxyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 556); 

                auto tg_yyzz_xxxyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 557); 

                auto tg_yyzz_xxxyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 558); 

                auto tg_yyzz_xxxyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 559); 

                auto tg_yyzz_xxxzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 560); 

                auto tg_yyzz_xxyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 561); 

                auto tg_yyzz_xxyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 562); 

                auto tg_yyzz_xxyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 563); 

                auto tg_yyzz_xxyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 564); 

                auto tg_yyzz_xxyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 565); 

                auto tg_yyzz_xxyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 566); 

                auto tg_yyzz_xxzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 567); 

                auto tg_yyzz_xyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 568); 

                auto tg_yyyy_xyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 385); 

                auto tg_yyyy_xyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 386); 

                auto tg_yyyy_xzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 387); 

                auto tg_yyyy_yyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 388); 

                auto tg_yyyy_yyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 389); 

                auto tg_yyyy_yyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 390); 

                auto tg_yyyy_yyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 391); 

                auto tg_yyyy_yyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 392); 

                auto tg_yyyy_yyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 393); 

                auto tg_yyyy_yzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 394); 

                auto tg_yyyy_zzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 395); 

                auto tg_yyyz_xxxxxxx_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 396); 

                auto tg_yyyz_xxxxxxy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 397); 

                auto tg_yyyz_xxxxxxz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 398); 

                auto tg_yyyz_xxxxxyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 399); 

                auto tg_yyyz_xxxxxyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 400); 

                auto tg_yyyz_xxxxxzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 401); 

                auto tg_yyyz_xxxxyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 402); 

                auto tg_yyyz_xxxxyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 403); 

                auto tg_yyyz_xxxxyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 404); 

                auto tg_yyyz_xxxxzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 405); 

                auto tg_yyyz_xxxyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 406); 

                auto tg_yyyz_xxxyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 407); 

                auto tg_yyyz_xxxyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 408); 

                auto tg_yyyz_xxxyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 409); 

                auto tg_yyyz_xxxzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 410); 

                auto tg_yyyz_xxyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 411); 

                auto tg_yyyz_xxyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 412); 

                auto tg_yyyz_xxyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 413); 

                auto tg_yyyz_xxyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 414); 

                auto tg_yyyz_xxyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 415); 

                auto tg_yyyz_xxzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 416); 

                auto tg_yyyz_xyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 417); 

                auto tg_yyyz_xyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 418); 

                auto tg_yyyz_xyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 419); 

                auto tg_yyyz_xyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 420); 

                auto tg_yyyz_xyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 421); 

                auto tg_yyyz_xyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 422); 

                auto tg_yyyz_xzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 423); 

                auto tg_yyyz_yyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 424); 

                auto tg_yyyz_yyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 425); 

                auto tg_yyyz_yyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 426); 

                auto tg_yyyz_yyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 427); 

                auto tg_yyyz_yyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 428); 

                auto tg_yyyz_yyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 429); 

                auto tg_yyyz_yzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 430); 

                auto tg_yyyz_zzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 431); 

                auto tg_yyzz_xxxxxxx_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 432); 

                auto tg_yyzz_xxxxxxy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 433); 

                auto tg_yyzz_xxxxxxz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 434); 

                auto tg_yyzz_xxxxxyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 435); 

                auto tg_yyzz_xxxxxyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 436); 

                auto tg_yyzz_xxxxxzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 437); 

                auto tg_yyzz_xxxxyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 438); 

                auto tg_yyzz_xxxxyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 439); 

                auto tg_yyzz_xxxxyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 440); 

                auto tg_yyzz_xxxxzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 441); 

                auto tg_yyzz_xxxyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 442); 

                auto tg_yyzz_xxxyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 443); 

                auto tg_yyzz_xxxyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 444); 

                auto tg_yyzz_xxxyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 445); 

                auto tg_yyzz_xxxzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 446); 

                auto tg_yyzz_xxyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 447); 

                auto tg_yyzz_xxyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 448); 

                auto tg_yyzz_xxyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 449); 

                auto tg_yyzz_xxyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 450); 

                auto tg_yyzz_xxyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 451); 

                auto tg_yyzz_xxzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 452); 

                auto tg_yyzz_xyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 453); 

                auto tg_yyzz_xyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 454); 

                auto tg_yyzz_xyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 455); 

                auto tg_yyzz_xyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 456); 

                auto tg_yyzz_xyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 457); 

                auto tg_yyzz_xyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 458); 

                auto tg_yyzz_xzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 459); 

                auto tg_yyzz_yyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 460); 

                // set up pointers to integrals

                auto tg_xyyyy_xxyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 475); 

                auto tg_xyyyy_xxyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 476); 

                auto tg_xyyyy_xxzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 477); 

                auto tg_xyyyy_xyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 478); 

                auto tg_xyyyy_xyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 479); 

                auto tg_xyyyy_xyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 480); 

                auto tg_xyyyy_xyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 481); 

                auto tg_xyyyy_xyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 482); 

                auto tg_xyyyy_xyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 483); 

                auto tg_xyyyy_xyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 484); 

                auto tg_xyyyy_xzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 485); 

                auto tg_xyyyy_yyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 486); 

                auto tg_xyyyy_yyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 487); 

                auto tg_xyyyy_yyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 488); 

                auto tg_xyyyy_yyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 489); 

                auto tg_xyyyy_yyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 490); 

                auto tg_xyyyy_yyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 491); 

                auto tg_xyyyy_yyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 492); 

                auto tg_xyyyy_yzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 493); 

                auto tg_xyyyy_zzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 494); 

                auto tg_xyyyz_xxxxxxxx_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 495); 

                auto tg_xyyyz_xxxxxxxy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 496); 

                auto tg_xyyyz_xxxxxxxz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 497); 

                auto tg_xyyyz_xxxxxxyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 498); 

                auto tg_xyyyz_xxxxxxyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 499); 

                auto tg_xyyyz_xxxxxxzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 500); 

                auto tg_xyyyz_xxxxxyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 501); 

                auto tg_xyyyz_xxxxxyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 502); 

                auto tg_xyyyz_xxxxxyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 503); 

                auto tg_xyyyz_xxxxxzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 504); 

                auto tg_xyyyz_xxxxyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 505); 

                auto tg_xyyyz_xxxxyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 506); 

                auto tg_xyyyz_xxxxyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 507); 

                auto tg_xyyyz_xxxxyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 508); 

                auto tg_xyyyz_xxxxzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 509); 

                auto tg_xyyyz_xxxyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 510); 

                auto tg_xyyyz_xxxyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 511); 

                auto tg_xyyyz_xxxyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 512); 

                auto tg_xyyyz_xxxyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 513); 

                auto tg_xyyyz_xxxyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 514); 

                auto tg_xyyyz_xxxzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 515); 

                auto tg_xyyyz_xxyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 516); 

                auto tg_xyyyz_xxyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 517); 

                auto tg_xyyyz_xxyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 518); 

                auto tg_xyyyz_xxyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 519); 

                auto tg_xyyyz_xxyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 520); 

                auto tg_xyyyz_xxyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 521); 

                auto tg_xyyyz_xxzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 522); 

                auto tg_xyyyz_xyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 523); 

                auto tg_xyyyz_xyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 524); 

                auto tg_xyyyz_xyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 525); 

                auto tg_xyyyz_xyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 526); 

                auto tg_xyyyz_xyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 527); 

                auto tg_xyyyz_xyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 528); 

                auto tg_xyyyz_xyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 529); 

                auto tg_xyyyz_xzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 530); 

                auto tg_xyyyz_yyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 531); 

                auto tg_xyyyz_yyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 532); 

                auto tg_xyyyz_yyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 533); 

                auto tg_xyyyz_yyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 534); 

                auto tg_xyyyz_yyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 535); 

                auto tg_xyyyz_yyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 536); 

                auto tg_xyyyz_yyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 537); 

                auto tg_xyyyz_yzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 538); 

                auto tg_xyyyz_zzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 539); 

                auto tg_xyyzz_xxxxxxxx_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 540); 

                auto tg_xyyzz_xxxxxxxy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 541); 

                auto tg_xyyzz_xxxxxxxz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 542); 

                auto tg_xyyzz_xxxxxxyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 543); 

                auto tg_xyyzz_xxxxxxyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 544); 

                auto tg_xyyzz_xxxxxxzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 545); 

                auto tg_xyyzz_xxxxxyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 546); 

                auto tg_xyyzz_xxxxxyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 547); 

                auto tg_xyyzz_xxxxxyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 548); 

                auto tg_xyyzz_xxxxxzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 549); 

                auto tg_xyyzz_xxxxyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 550); 

                auto tg_xyyzz_xxxxyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 551); 

                auto tg_xyyzz_xxxxyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 552); 

                auto tg_xyyzz_xxxxyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 553); 

                auto tg_xyyzz_xxxxzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 554); 

                auto tg_xyyzz_xxxyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 555); 

                auto tg_xyyzz_xxxyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 556); 

                auto tg_xyyzz_xxxyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 557); 

                auto tg_xyyzz_xxxyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 558); 

                auto tg_xyyzz_xxxyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 559); 

                auto tg_xyyzz_xxxzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 560); 

                auto tg_xyyzz_xxyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 561); 

                auto tg_xyyzz_xxyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 562); 

                auto tg_xyyzz_xxyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 563); 

                auto tg_xyyzz_xxyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 564); 

                auto tg_xyyzz_xxyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 565); 

                auto tg_xyyzz_xxyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 566); 

                auto tg_xyyzz_xxzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 567); 

                auto tg_xyyzz_xyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 568); 

                // Batch of Integrals (475,569)

                #pragma omp simd aligned(fxn, tg_xyyyy_xxyyzzzz_0, tg_xyyyy_xxyzzzzz_0, tg_xyyyy_xxzzzzzz_0, \
                                         tg_xyyyy_xyyyyyyy_0, tg_xyyyy_xyyyyyyz_0, tg_xyyyy_xyyyyyzz_0, tg_xyyyy_xyyyyzzz_0, \
                                         tg_xyyyy_xyyyzzzz_0, tg_xyyyy_xyyzzzzz_0, tg_xyyyy_xyzzzzzz_0, tg_xyyyy_xzzzzzzz_0, \
                                         tg_xyyyy_yyyyyyyy_0, tg_xyyyy_yyyyyyyz_0, tg_xyyyy_yyyyyyzz_0, tg_xyyyy_yyyyyzzz_0, \
                                         tg_xyyyy_yyyyzzzz_0, tg_xyyyy_yyyzzzzz_0, tg_xyyyy_yyzzzzzz_0, tg_xyyyy_yzzzzzzz_0, \
                                         tg_xyyyy_zzzzzzzz_0, tg_xyyyz_xxxxxxxx_0, tg_xyyyz_xxxxxxxy_0, tg_xyyyz_xxxxxxxz_0, \
                                         tg_xyyyz_xxxxxxyy_0, tg_xyyyz_xxxxxxyz_0, tg_xyyyz_xxxxxxzz_0, tg_xyyyz_xxxxxyyy_0, \
                                         tg_xyyyz_xxxxxyyz_0, tg_xyyyz_xxxxxyzz_0, tg_xyyyz_xxxxxzzz_0, tg_xyyyz_xxxxyyyy_0, \
                                         tg_xyyyz_xxxxyyyz_0, tg_xyyyz_xxxxyyzz_0, tg_xyyyz_xxxxyzzz_0, tg_xyyyz_xxxxzzzz_0, \
                                         tg_xyyyz_xxxyyyyy_0, tg_xyyyz_xxxyyyyz_0, tg_xyyyz_xxxyyyzz_0, tg_xyyyz_xxxyyzzz_0, \
                                         tg_xyyyz_xxxyzzzz_0, tg_xyyyz_xxxzzzzz_0, tg_xyyyz_xxyyyyyy_0, tg_xyyyz_xxyyyyyz_0, \
                                         tg_xyyyz_xxyyyyzz_0, tg_xyyyz_xxyyyzzz_0, tg_xyyyz_xxyyzzzz_0, tg_xyyyz_xxyzzzzz_0, \
                                         tg_xyyyz_xxzzzzzz_0, tg_xyyyz_xyyyyyyy_0, tg_xyyyz_xyyyyyyz_0, tg_xyyyz_xyyyyyzz_0, \
                                         tg_xyyyz_xyyyyzzz_0, tg_xyyyz_xyyyzzzz_0, tg_xyyyz_xyyzzzzz_0, tg_xyyyz_xyzzzzzz_0, \
                                         tg_xyyyz_xzzzzzzz_0, tg_xyyyz_yyyyyyyy_0, tg_xyyyz_yyyyyyyz_0, tg_xyyyz_yyyyyyzz_0, \
                                         tg_xyyyz_yyyyyzzz_0, tg_xyyyz_yyyyzzzz_0, tg_xyyyz_yyyzzzzz_0, tg_xyyyz_yyzzzzzz_0, \
                                         tg_xyyyz_yzzzzzzz_0, tg_xyyyz_zzzzzzzz_0, tg_xyyzz_xxxxxxxx_0, tg_xyyzz_xxxxxxxy_0, \
                                         tg_xyyzz_xxxxxxxz_0, tg_xyyzz_xxxxxxyy_0, tg_xyyzz_xxxxxxyz_0, tg_xyyzz_xxxxxxzz_0, \
                                         tg_xyyzz_xxxxxyyy_0, tg_xyyzz_xxxxxyyz_0, tg_xyyzz_xxxxxyzz_0, tg_xyyzz_xxxxxzzz_0, \
                                         tg_xyyzz_xxxxyyyy_0, tg_xyyzz_xxxxyyyz_0, tg_xyyzz_xxxxyyzz_0, tg_xyyzz_xxxxyzzz_0, \
                                         tg_xyyzz_xxxxzzzz_0, tg_xyyzz_xxxyyyyy_0, tg_xyyzz_xxxyyyyz_0, tg_xyyzz_xxxyyyzz_0, \
                                         tg_xyyzz_xxxyyzzz_0, tg_xyyzz_xxxyzzzz_0, tg_xyyzz_xxxzzzzz_0, tg_xyyzz_xxyyyyyy_0, \
                                         tg_xyyzz_xxyyyyyz_0, tg_xyyzz_xxyyyyzz_0, tg_xyyzz_xxyyyzzz_0, tg_xyyzz_xxyyzzzz_0, \
                                         tg_xyyzz_xxyzzzzz_0, tg_xyyzz_xxzzzzzz_0, tg_xyyzz_xyyyyyyy_0, tg_yyyy_xxyyzzzz_0, \
                                         tg_yyyy_xxyyzzzz_1, tg_yyyy_xxyzzzzz_0, tg_yyyy_xxyzzzzz_1, tg_yyyy_xxzzzzzz_0, \
                                         tg_yyyy_xxzzzzzz_1, tg_yyyy_xyyyyyyy_0, tg_yyyy_xyyyyyyy_1, tg_yyyy_xyyyyyyz_0, \
                                         tg_yyyy_xyyyyyyz_1, tg_yyyy_xyyyyyzz_0, tg_yyyy_xyyyyyzz_1, tg_yyyy_xyyyyzzz_0, \
                                         tg_yyyy_xyyyyzzz_1, tg_yyyy_xyyyzzzz_0, tg_yyyy_xyyyzzzz_1, tg_yyyy_xyyzzzz_1, \
                                         tg_yyyy_xyyzzzzz_0, tg_yyyy_xyyzzzzz_1, tg_yyyy_xyzzzzz_1, tg_yyyy_xyzzzzzz_0, \
                                         tg_yyyy_xyzzzzzz_1, tg_yyyy_xzzzzzz_1, tg_yyyy_xzzzzzzz_0, tg_yyyy_xzzzzzzz_1, \
                                         tg_yyyy_yyyyyyy_1, tg_yyyy_yyyyyyyy_0, tg_yyyy_yyyyyyyy_1, tg_yyyy_yyyyyyyz_0, \
                                         tg_yyyy_yyyyyyyz_1, tg_yyyy_yyyyyyz_1, tg_yyyy_yyyyyyzz_0, tg_yyyy_yyyyyyzz_1, \
                                         tg_yyyy_yyyyyzz_1, tg_yyyy_yyyyyzzz_0, tg_yyyy_yyyyyzzz_1, tg_yyyy_yyyyzzz_1, \
                                         tg_yyyy_yyyyzzzz_0, tg_yyyy_yyyyzzzz_1, tg_yyyy_yyyzzzz_1, tg_yyyy_yyyzzzzz_0, \
                                         tg_yyyy_yyyzzzzz_1, tg_yyyy_yyzzzzz_1, tg_yyyy_yyzzzzzz_0, tg_yyyy_yyzzzzzz_1, \
                                         tg_yyyy_yzzzzzz_1, tg_yyyy_yzzzzzzz_0, tg_yyyy_yzzzzzzz_1, tg_yyyy_zzzzzzz_1, \
                                         tg_yyyy_zzzzzzzz_0, tg_yyyy_zzzzzzzz_1, tg_yyyz_xxxxxxx_1, tg_yyyz_xxxxxxxx_0, \
                                         tg_yyyz_xxxxxxxx_1, tg_yyyz_xxxxxxxy_0, tg_yyyz_xxxxxxxy_1, tg_yyyz_xxxxxxxz_0, \
                                         tg_yyyz_xxxxxxxz_1, tg_yyyz_xxxxxxy_1, tg_yyyz_xxxxxxyy_0, tg_yyyz_xxxxxxyy_1, \
                                         tg_yyyz_xxxxxxyz_0, tg_yyyz_xxxxxxyz_1, tg_yyyz_xxxxxxz_1, tg_yyyz_xxxxxxzz_0, \
                                         tg_yyyz_xxxxxxzz_1, tg_yyyz_xxxxxyy_1, tg_yyyz_xxxxxyyy_0, tg_yyyz_xxxxxyyy_1, \
                                         tg_yyyz_xxxxxyyz_0, tg_yyyz_xxxxxyyz_1, tg_yyyz_xxxxxyz_1, tg_yyyz_xxxxxyzz_0, \
                                         tg_yyyz_xxxxxyzz_1, tg_yyyz_xxxxxzz_1, tg_yyyz_xxxxxzzz_0, tg_yyyz_xxxxxzzz_1, \
                                         tg_yyyz_xxxxyyy_1, tg_yyyz_xxxxyyyy_0, tg_yyyz_xxxxyyyy_1, tg_yyyz_xxxxyyyz_0, \
                                         tg_yyyz_xxxxyyyz_1, tg_yyyz_xxxxyyz_1, tg_yyyz_xxxxyyzz_0, tg_yyyz_xxxxyyzz_1, \
                                         tg_yyyz_xxxxyzz_1, tg_yyyz_xxxxyzzz_0, tg_yyyz_xxxxyzzz_1, tg_yyyz_xxxxzzz_1, \
                                         tg_yyyz_xxxxzzzz_0, tg_yyyz_xxxxzzzz_1, tg_yyyz_xxxyyyy_1, tg_yyyz_xxxyyyyy_0, \
                                         tg_yyyz_xxxyyyyy_1, tg_yyyz_xxxyyyyz_0, tg_yyyz_xxxyyyyz_1, tg_yyyz_xxxyyyz_1, \
                                         tg_yyyz_xxxyyyzz_0, tg_yyyz_xxxyyyzz_1, tg_yyyz_xxxyyzz_1, tg_yyyz_xxxyyzzz_0, \
                                         tg_yyyz_xxxyyzzz_1, tg_yyyz_xxxyzzz_1, tg_yyyz_xxxyzzzz_0, tg_yyyz_xxxyzzzz_1, \
                                         tg_yyyz_xxxzzzz_1, tg_yyyz_xxxzzzzz_0, tg_yyyz_xxxzzzzz_1, tg_yyyz_xxyyyyy_1, \
                                         tg_yyyz_xxyyyyyy_0, tg_yyyz_xxyyyyyy_1, tg_yyyz_xxyyyyyz_0, tg_yyyz_xxyyyyyz_1, \
                                         tg_yyyz_xxyyyyz_1, tg_yyyz_xxyyyyzz_0, tg_yyyz_xxyyyyzz_1, tg_yyyz_xxyyyzz_1, \
                                         tg_yyyz_xxyyyzzz_0, tg_yyyz_xxyyyzzz_1, tg_yyyz_xxyyzzz_1, tg_yyyz_xxyyzzzz_0, \
                                         tg_yyyz_xxyyzzzz_1, tg_yyyz_xxyzzzz_1, tg_yyyz_xxyzzzzz_0, tg_yyyz_xxyzzzzz_1, \
                                         tg_yyyz_xxzzzzz_1, tg_yyyz_xxzzzzzz_0, tg_yyyz_xxzzzzzz_1, tg_yyyz_xyyyyyy_1, \
                                         tg_yyyz_xyyyyyyy_0, tg_yyyz_xyyyyyyy_1, tg_yyyz_xyyyyyyz_0, tg_yyyz_xyyyyyyz_1, \
                                         tg_yyyz_xyyyyyz_1, tg_yyyz_xyyyyyzz_0, tg_yyyz_xyyyyyzz_1, tg_yyyz_xyyyyzz_1, \
                                         tg_yyyz_xyyyyzzz_0, tg_yyyz_xyyyyzzz_1, tg_yyyz_xyyyzzz_1, tg_yyyz_xyyyzzzz_0, \
                                         tg_yyyz_xyyyzzzz_1, tg_yyyz_xyyzzzz_1, tg_yyyz_xyyzzzzz_0, tg_yyyz_xyyzzzzz_1, \
                                         tg_yyyz_xyzzzzz_1, tg_yyyz_xyzzzzzz_0, tg_yyyz_xyzzzzzz_1, tg_yyyz_xzzzzzz_1, \
                                         tg_yyyz_xzzzzzzz_0, tg_yyyz_xzzzzzzz_1, tg_yyyz_yyyyyyy_1, tg_yyyz_yyyyyyyy_0, \
                                         tg_yyyz_yyyyyyyy_1, tg_yyyz_yyyyyyyz_0, tg_yyyz_yyyyyyyz_1, tg_yyyz_yyyyyyz_1, \
                                         tg_yyyz_yyyyyyzz_0, tg_yyyz_yyyyyyzz_1, tg_yyyz_yyyyyzz_1, tg_yyyz_yyyyyzzz_0, \
                                         tg_yyyz_yyyyyzzz_1, tg_yyyz_yyyyzzz_1, tg_yyyz_yyyyzzzz_0, tg_yyyz_yyyyzzzz_1, \
                                         tg_yyyz_yyyzzzz_1, tg_yyyz_yyyzzzzz_0, tg_yyyz_yyyzzzzz_1, tg_yyyz_yyzzzzz_1, \
                                         tg_yyyz_yyzzzzzz_0, tg_yyyz_yyzzzzzz_1, tg_yyyz_yzzzzzz_1, tg_yyyz_yzzzzzzz_0, \
                                         tg_yyyz_yzzzzzzz_1, tg_yyyz_zzzzzzz_1, tg_yyyz_zzzzzzzz_0, tg_yyyz_zzzzzzzz_1, \
                                         tg_yyzz_xxxxxxx_1, tg_yyzz_xxxxxxxx_0, tg_yyzz_xxxxxxxx_1, tg_yyzz_xxxxxxxy_0, \
                                         tg_yyzz_xxxxxxxy_1, tg_yyzz_xxxxxxxz_0, tg_yyzz_xxxxxxxz_1, tg_yyzz_xxxxxxy_1, \
                                         tg_yyzz_xxxxxxyy_0, tg_yyzz_xxxxxxyy_1, tg_yyzz_xxxxxxyz_0, tg_yyzz_xxxxxxyz_1, \
                                         tg_yyzz_xxxxxxz_1, tg_yyzz_xxxxxxzz_0, tg_yyzz_xxxxxxzz_1, tg_yyzz_xxxxxyy_1, \
                                         tg_yyzz_xxxxxyyy_0, tg_yyzz_xxxxxyyy_1, tg_yyzz_xxxxxyyz_0, tg_yyzz_xxxxxyyz_1, \
                                         tg_yyzz_xxxxxyz_1, tg_yyzz_xxxxxyzz_0, tg_yyzz_xxxxxyzz_1, tg_yyzz_xxxxxzz_1, \
                                         tg_yyzz_xxxxxzzz_0, tg_yyzz_xxxxxzzz_1, tg_yyzz_xxxxyyy_1, tg_yyzz_xxxxyyyy_0, \
                                         tg_yyzz_xxxxyyyy_1, tg_yyzz_xxxxyyyz_0, tg_yyzz_xxxxyyyz_1, tg_yyzz_xxxxyyz_1, \
                                         tg_yyzz_xxxxyyzz_0, tg_yyzz_xxxxyyzz_1, tg_yyzz_xxxxyzz_1, tg_yyzz_xxxxyzzz_0, \
                                         tg_yyzz_xxxxyzzz_1, tg_yyzz_xxxxzzz_1, tg_yyzz_xxxxzzzz_0, tg_yyzz_xxxxzzzz_1, \
                                         tg_yyzz_xxxyyyy_1, tg_yyzz_xxxyyyyy_0, tg_yyzz_xxxyyyyy_1, tg_yyzz_xxxyyyyz_0, \
                                         tg_yyzz_xxxyyyyz_1, tg_yyzz_xxxyyyz_1, tg_yyzz_xxxyyyzz_0, tg_yyzz_xxxyyyzz_1, \
                                         tg_yyzz_xxxyyzz_1, tg_yyzz_xxxyyzzz_0, tg_yyzz_xxxyyzzz_1, tg_yyzz_xxxyzzz_1, \
                                         tg_yyzz_xxxyzzzz_0, tg_yyzz_xxxyzzzz_1, tg_yyzz_xxxzzzz_1, tg_yyzz_xxxzzzzz_0, \
                                         tg_yyzz_xxxzzzzz_1, tg_yyzz_xxyyyyy_1, tg_yyzz_xxyyyyyy_0, tg_yyzz_xxyyyyyy_1, \
                                         tg_yyzz_xxyyyyyz_0, tg_yyzz_xxyyyyyz_1, tg_yyzz_xxyyyyz_1, tg_yyzz_xxyyyyzz_0, \
                                         tg_yyzz_xxyyyyzz_1, tg_yyzz_xxyyyzz_1, tg_yyzz_xxyyyzzz_0, tg_yyzz_xxyyyzzz_1, \
                                         tg_yyzz_xxyyzzz_1, tg_yyzz_xxyyzzzz_0, tg_yyzz_xxyyzzzz_1, tg_yyzz_xxyzzzz_1, \
                                         tg_yyzz_xxyzzzzz_0, tg_yyzz_xxyzzzzz_1, tg_yyzz_xxzzzzz_1, tg_yyzz_xxzzzzzz_0, \
                                         tg_yyzz_xxzzzzzz_1, tg_yyzz_xyyyyyy_1, tg_yyzz_xyyyyyyy_0, tg_yyzz_xyyyyyyy_1, \
                                         tg_yyzz_xyyyyyz_1, tg_yyzz_xyyyyzz_1, tg_yyzz_xyyyzzz_1, tg_yyzz_xyyzzzz_1, \
                                         tg_yyzz_xyzzzzz_1, tg_yyzz_xzzzzzz_1, tg_yyzz_yyyyyyy_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fxn = fxn[j];

                    double fr = wp_x[j]; 

                    tg_xyyyy_xxyyzzzz_0[j] = pb_x * tg_yyyy_xxyyzzzz_0[j] + fr * tg_yyyy_xxyyzzzz_1[j] + fl1_fxn * tg_yyyy_xyyzzzz_1[j];

                    tg_xyyyy_xxyzzzzz_0[j] = pb_x * tg_yyyy_xxyzzzzz_0[j] + fr * tg_yyyy_xxyzzzzz_1[j] + fl1_fxn * tg_yyyy_xyzzzzz_1[j];

                    tg_xyyyy_xxzzzzzz_0[j] = pb_x * tg_yyyy_xxzzzzzz_0[j] + fr * tg_yyyy_xxzzzzzz_1[j] + fl1_fxn * tg_yyyy_xzzzzzz_1[j];

                    tg_xyyyy_xyyyyyyy_0[j] = pb_x * tg_yyyy_xyyyyyyy_0[j] + fr * tg_yyyy_xyyyyyyy_1[j] + 0.5 * fl1_fxn * tg_yyyy_yyyyyyy_1[j];

                    tg_xyyyy_xyyyyyyz_0[j] = pb_x * tg_yyyy_xyyyyyyz_0[j] + fr * tg_yyyy_xyyyyyyz_1[j] + 0.5 * fl1_fxn * tg_yyyy_yyyyyyz_1[j];

                    tg_xyyyy_xyyyyyzz_0[j] = pb_x * tg_yyyy_xyyyyyzz_0[j] + fr * tg_yyyy_xyyyyyzz_1[j] + 0.5 * fl1_fxn * tg_yyyy_yyyyyzz_1[j];

                    tg_xyyyy_xyyyyzzz_0[j] = pb_x * tg_yyyy_xyyyyzzz_0[j] + fr * tg_yyyy_xyyyyzzz_1[j] + 0.5 * fl1_fxn * tg_yyyy_yyyyzzz_1[j];

                    tg_xyyyy_xyyyzzzz_0[j] = pb_x * tg_yyyy_xyyyzzzz_0[j] + fr * tg_yyyy_xyyyzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyy_yyyzzzz_1[j];

                    tg_xyyyy_xyyzzzzz_0[j] = pb_x * tg_yyyy_xyyzzzzz_0[j] + fr * tg_yyyy_xyyzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyy_yyzzzzz_1[j];

                    tg_xyyyy_xyzzzzzz_0[j] = pb_x * tg_yyyy_xyzzzzzz_0[j] + fr * tg_yyyy_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyy_yzzzzzz_1[j];

                    tg_xyyyy_xzzzzzzz_0[j] = pb_x * tg_yyyy_xzzzzzzz_0[j] + fr * tg_yyyy_xzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyy_zzzzzzz_1[j];

                    tg_xyyyy_yyyyyyyy_0[j] = pb_x * tg_yyyy_yyyyyyyy_0[j] + fr * tg_yyyy_yyyyyyyy_1[j];

                    tg_xyyyy_yyyyyyyz_0[j] = pb_x * tg_yyyy_yyyyyyyz_0[j] + fr * tg_yyyy_yyyyyyyz_1[j];

                    tg_xyyyy_yyyyyyzz_0[j] = pb_x * tg_yyyy_yyyyyyzz_0[j] + fr * tg_yyyy_yyyyyyzz_1[j];

                    tg_xyyyy_yyyyyzzz_0[j] = pb_x * tg_yyyy_yyyyyzzz_0[j] + fr * tg_yyyy_yyyyyzzz_1[j];

                    tg_xyyyy_yyyyzzzz_0[j] = pb_x * tg_yyyy_yyyyzzzz_0[j] + fr * tg_yyyy_yyyyzzzz_1[j];

                    tg_xyyyy_yyyzzzzz_0[j] = pb_x * tg_yyyy_yyyzzzzz_0[j] + fr * tg_yyyy_yyyzzzzz_1[j];

                    tg_xyyyy_yyzzzzzz_0[j] = pb_x * tg_yyyy_yyzzzzzz_0[j] + fr * tg_yyyy_yyzzzzzz_1[j];

                    tg_xyyyy_yzzzzzzz_0[j] = pb_x * tg_yyyy_yzzzzzzz_0[j] + fr * tg_yyyy_yzzzzzzz_1[j];

                    tg_xyyyy_zzzzzzzz_0[j] = pb_x * tg_yyyy_zzzzzzzz_0[j] + fr * tg_yyyy_zzzzzzzz_1[j];

                    tg_xyyyz_xxxxxxxx_0[j] = pb_x * tg_yyyz_xxxxxxxx_0[j] + fr * tg_yyyz_xxxxxxxx_1[j] + 4.0 * fl1_fxn * tg_yyyz_xxxxxxx_1[j];

                    tg_xyyyz_xxxxxxxy_0[j] = pb_x * tg_yyyz_xxxxxxxy_0[j] + fr * tg_yyyz_xxxxxxxy_1[j] + 3.5 * fl1_fxn * tg_yyyz_xxxxxxy_1[j];

                    tg_xyyyz_xxxxxxxz_0[j] = pb_x * tg_yyyz_xxxxxxxz_0[j] + fr * tg_yyyz_xxxxxxxz_1[j] + 3.5 * fl1_fxn * tg_yyyz_xxxxxxz_1[j];

                    tg_xyyyz_xxxxxxyy_0[j] = pb_x * tg_yyyz_xxxxxxyy_0[j] + fr * tg_yyyz_xxxxxxyy_1[j] + 3.0 * fl1_fxn * tg_yyyz_xxxxxyy_1[j];

                    tg_xyyyz_xxxxxxyz_0[j] = pb_x * tg_yyyz_xxxxxxyz_0[j] + fr * tg_yyyz_xxxxxxyz_1[j] + 3.0 * fl1_fxn * tg_yyyz_xxxxxyz_1[j];

                    tg_xyyyz_xxxxxxzz_0[j] = pb_x * tg_yyyz_xxxxxxzz_0[j] + fr * tg_yyyz_xxxxxxzz_1[j] + 3.0 * fl1_fxn * tg_yyyz_xxxxxzz_1[j];

                    tg_xyyyz_xxxxxyyy_0[j] = pb_x * tg_yyyz_xxxxxyyy_0[j] + fr * tg_yyyz_xxxxxyyy_1[j] + 2.5 * fl1_fxn * tg_yyyz_xxxxyyy_1[j];

                    tg_xyyyz_xxxxxyyz_0[j] = pb_x * tg_yyyz_xxxxxyyz_0[j] + fr * tg_yyyz_xxxxxyyz_1[j] + 2.5 * fl1_fxn * tg_yyyz_xxxxyyz_1[j];

                    tg_xyyyz_xxxxxyzz_0[j] = pb_x * tg_yyyz_xxxxxyzz_0[j] + fr * tg_yyyz_xxxxxyzz_1[j] + 2.5 * fl1_fxn * tg_yyyz_xxxxyzz_1[j];

                    tg_xyyyz_xxxxxzzz_0[j] = pb_x * tg_yyyz_xxxxxzzz_0[j] + fr * tg_yyyz_xxxxxzzz_1[j] + 2.5 * fl1_fxn * tg_yyyz_xxxxzzz_1[j];

                    tg_xyyyz_xxxxyyyy_0[j] = pb_x * tg_yyyz_xxxxyyyy_0[j] + fr * tg_yyyz_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_yyyz_xxxyyyy_1[j];

                    tg_xyyyz_xxxxyyyz_0[j] = pb_x * tg_yyyz_xxxxyyyz_0[j] + fr * tg_yyyz_xxxxyyyz_1[j] + 2.0 * fl1_fxn * tg_yyyz_xxxyyyz_1[j];

                    tg_xyyyz_xxxxyyzz_0[j] = pb_x * tg_yyyz_xxxxyyzz_0[j] + fr * tg_yyyz_xxxxyyzz_1[j] + 2.0 * fl1_fxn * tg_yyyz_xxxyyzz_1[j];

                    tg_xyyyz_xxxxyzzz_0[j] = pb_x * tg_yyyz_xxxxyzzz_0[j] + fr * tg_yyyz_xxxxyzzz_1[j] + 2.0 * fl1_fxn * tg_yyyz_xxxyzzz_1[j];

                    tg_xyyyz_xxxxzzzz_0[j] = pb_x * tg_yyyz_xxxxzzzz_0[j] + fr * tg_yyyz_xxxxzzzz_1[j] + 2.0 * fl1_fxn * tg_yyyz_xxxzzzz_1[j];

                    tg_xyyyz_xxxyyyyy_0[j] = pb_x * tg_yyyz_xxxyyyyy_0[j] + fr * tg_yyyz_xxxyyyyy_1[j] + 1.5 * fl1_fxn * tg_yyyz_xxyyyyy_1[j];

                    tg_xyyyz_xxxyyyyz_0[j] = pb_x * tg_yyyz_xxxyyyyz_0[j] + fr * tg_yyyz_xxxyyyyz_1[j] + 1.5 * fl1_fxn * tg_yyyz_xxyyyyz_1[j];

                    tg_xyyyz_xxxyyyzz_0[j] = pb_x * tg_yyyz_xxxyyyzz_0[j] + fr * tg_yyyz_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_yyyz_xxyyyzz_1[j];

                    tg_xyyyz_xxxyyzzz_0[j] = pb_x * tg_yyyz_xxxyyzzz_0[j] + fr * tg_yyyz_xxxyyzzz_1[j] + 1.5 * fl1_fxn * tg_yyyz_xxyyzzz_1[j];

                    tg_xyyyz_xxxyzzzz_0[j] = pb_x * tg_yyyz_xxxyzzzz_0[j] + fr * tg_yyyz_xxxyzzzz_1[j] + 1.5 * fl1_fxn * tg_yyyz_xxyzzzz_1[j];

                    tg_xyyyz_xxxzzzzz_0[j] = pb_x * tg_yyyz_xxxzzzzz_0[j] + fr * tg_yyyz_xxxzzzzz_1[j] + 1.5 * fl1_fxn * tg_yyyz_xxzzzzz_1[j];

                    tg_xyyyz_xxyyyyyy_0[j] = pb_x * tg_yyyz_xxyyyyyy_0[j] + fr * tg_yyyz_xxyyyyyy_1[j] + fl1_fxn * tg_yyyz_xyyyyyy_1[j];

                    tg_xyyyz_xxyyyyyz_0[j] = pb_x * tg_yyyz_xxyyyyyz_0[j] + fr * tg_yyyz_xxyyyyyz_1[j] + fl1_fxn * tg_yyyz_xyyyyyz_1[j];

                    tg_xyyyz_xxyyyyzz_0[j] = pb_x * tg_yyyz_xxyyyyzz_0[j] + fr * tg_yyyz_xxyyyyzz_1[j] + fl1_fxn * tg_yyyz_xyyyyzz_1[j];

                    tg_xyyyz_xxyyyzzz_0[j] = pb_x * tg_yyyz_xxyyyzzz_0[j] + fr * tg_yyyz_xxyyyzzz_1[j] + fl1_fxn * tg_yyyz_xyyyzzz_1[j];

                    tg_xyyyz_xxyyzzzz_0[j] = pb_x * tg_yyyz_xxyyzzzz_0[j] + fr * tg_yyyz_xxyyzzzz_1[j] + fl1_fxn * tg_yyyz_xyyzzzz_1[j];

                    tg_xyyyz_xxyzzzzz_0[j] = pb_x * tg_yyyz_xxyzzzzz_0[j] + fr * tg_yyyz_xxyzzzzz_1[j] + fl1_fxn * tg_yyyz_xyzzzzz_1[j];

                    tg_xyyyz_xxzzzzzz_0[j] = pb_x * tg_yyyz_xxzzzzzz_0[j] + fr * tg_yyyz_xxzzzzzz_1[j] + fl1_fxn * tg_yyyz_xzzzzzz_1[j];

                    tg_xyyyz_xyyyyyyy_0[j] = pb_x * tg_yyyz_xyyyyyyy_0[j] + fr * tg_yyyz_xyyyyyyy_1[j] + 0.5 * fl1_fxn * tg_yyyz_yyyyyyy_1[j];

                    tg_xyyyz_xyyyyyyz_0[j] = pb_x * tg_yyyz_xyyyyyyz_0[j] + fr * tg_yyyz_xyyyyyyz_1[j] + 0.5 * fl1_fxn * tg_yyyz_yyyyyyz_1[j];

                    tg_xyyyz_xyyyyyzz_0[j] = pb_x * tg_yyyz_xyyyyyzz_0[j] + fr * tg_yyyz_xyyyyyzz_1[j] + 0.5 * fl1_fxn * tg_yyyz_yyyyyzz_1[j];

                    tg_xyyyz_xyyyyzzz_0[j] = pb_x * tg_yyyz_xyyyyzzz_0[j] + fr * tg_yyyz_xyyyyzzz_1[j] + 0.5 * fl1_fxn * tg_yyyz_yyyyzzz_1[j];

                    tg_xyyyz_xyyyzzzz_0[j] = pb_x * tg_yyyz_xyyyzzzz_0[j] + fr * tg_yyyz_xyyyzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyz_yyyzzzz_1[j];

                    tg_xyyyz_xyyzzzzz_0[j] = pb_x * tg_yyyz_xyyzzzzz_0[j] + fr * tg_yyyz_xyyzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyz_yyzzzzz_1[j];

                    tg_xyyyz_xyzzzzzz_0[j] = pb_x * tg_yyyz_xyzzzzzz_0[j] + fr * tg_yyyz_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyz_yzzzzzz_1[j];

                    tg_xyyyz_xzzzzzzz_0[j] = pb_x * tg_yyyz_xzzzzzzz_0[j] + fr * tg_yyyz_xzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyz_zzzzzzz_1[j];

                    tg_xyyyz_yyyyyyyy_0[j] = pb_x * tg_yyyz_yyyyyyyy_0[j] + fr * tg_yyyz_yyyyyyyy_1[j];

                    tg_xyyyz_yyyyyyyz_0[j] = pb_x * tg_yyyz_yyyyyyyz_0[j] + fr * tg_yyyz_yyyyyyyz_1[j];

                    tg_xyyyz_yyyyyyzz_0[j] = pb_x * tg_yyyz_yyyyyyzz_0[j] + fr * tg_yyyz_yyyyyyzz_1[j];

                    tg_xyyyz_yyyyyzzz_0[j] = pb_x * tg_yyyz_yyyyyzzz_0[j] + fr * tg_yyyz_yyyyyzzz_1[j];

                    tg_xyyyz_yyyyzzzz_0[j] = pb_x * tg_yyyz_yyyyzzzz_0[j] + fr * tg_yyyz_yyyyzzzz_1[j];

                    tg_xyyyz_yyyzzzzz_0[j] = pb_x * tg_yyyz_yyyzzzzz_0[j] + fr * tg_yyyz_yyyzzzzz_1[j];

                    tg_xyyyz_yyzzzzzz_0[j] = pb_x * tg_yyyz_yyzzzzzz_0[j] + fr * tg_yyyz_yyzzzzzz_1[j];

                    tg_xyyyz_yzzzzzzz_0[j] = pb_x * tg_yyyz_yzzzzzzz_0[j] + fr * tg_yyyz_yzzzzzzz_1[j];

                    tg_xyyyz_zzzzzzzz_0[j] = pb_x * tg_yyyz_zzzzzzzz_0[j] + fr * tg_yyyz_zzzzzzzz_1[j];

                    tg_xyyzz_xxxxxxxx_0[j] = pb_x * tg_yyzz_xxxxxxxx_0[j] + fr * tg_yyzz_xxxxxxxx_1[j] + 4.0 * fl1_fxn * tg_yyzz_xxxxxxx_1[j];

                    tg_xyyzz_xxxxxxxy_0[j] = pb_x * tg_yyzz_xxxxxxxy_0[j] + fr * tg_yyzz_xxxxxxxy_1[j] + 3.5 * fl1_fxn * tg_yyzz_xxxxxxy_1[j];

                    tg_xyyzz_xxxxxxxz_0[j] = pb_x * tg_yyzz_xxxxxxxz_0[j] + fr * tg_yyzz_xxxxxxxz_1[j] + 3.5 * fl1_fxn * tg_yyzz_xxxxxxz_1[j];

                    tg_xyyzz_xxxxxxyy_0[j] = pb_x * tg_yyzz_xxxxxxyy_0[j] + fr * tg_yyzz_xxxxxxyy_1[j] + 3.0 * fl1_fxn * tg_yyzz_xxxxxyy_1[j];

                    tg_xyyzz_xxxxxxyz_0[j] = pb_x * tg_yyzz_xxxxxxyz_0[j] + fr * tg_yyzz_xxxxxxyz_1[j] + 3.0 * fl1_fxn * tg_yyzz_xxxxxyz_1[j];

                    tg_xyyzz_xxxxxxzz_0[j] = pb_x * tg_yyzz_xxxxxxzz_0[j] + fr * tg_yyzz_xxxxxxzz_1[j] + 3.0 * fl1_fxn * tg_yyzz_xxxxxzz_1[j];

                    tg_xyyzz_xxxxxyyy_0[j] = pb_x * tg_yyzz_xxxxxyyy_0[j] + fr * tg_yyzz_xxxxxyyy_1[j] + 2.5 * fl1_fxn * tg_yyzz_xxxxyyy_1[j];

                    tg_xyyzz_xxxxxyyz_0[j] = pb_x * tg_yyzz_xxxxxyyz_0[j] + fr * tg_yyzz_xxxxxyyz_1[j] + 2.5 * fl1_fxn * tg_yyzz_xxxxyyz_1[j];

                    tg_xyyzz_xxxxxyzz_0[j] = pb_x * tg_yyzz_xxxxxyzz_0[j] + fr * tg_yyzz_xxxxxyzz_1[j] + 2.5 * fl1_fxn * tg_yyzz_xxxxyzz_1[j];

                    tg_xyyzz_xxxxxzzz_0[j] = pb_x * tg_yyzz_xxxxxzzz_0[j] + fr * tg_yyzz_xxxxxzzz_1[j] + 2.5 * fl1_fxn * tg_yyzz_xxxxzzz_1[j];

                    tg_xyyzz_xxxxyyyy_0[j] = pb_x * tg_yyzz_xxxxyyyy_0[j] + fr * tg_yyzz_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_yyzz_xxxyyyy_1[j];

                    tg_xyyzz_xxxxyyyz_0[j] = pb_x * tg_yyzz_xxxxyyyz_0[j] + fr * tg_yyzz_xxxxyyyz_1[j] + 2.0 * fl1_fxn * tg_yyzz_xxxyyyz_1[j];

                    tg_xyyzz_xxxxyyzz_0[j] = pb_x * tg_yyzz_xxxxyyzz_0[j] + fr * tg_yyzz_xxxxyyzz_1[j] + 2.0 * fl1_fxn * tg_yyzz_xxxyyzz_1[j];

                    tg_xyyzz_xxxxyzzz_0[j] = pb_x * tg_yyzz_xxxxyzzz_0[j] + fr * tg_yyzz_xxxxyzzz_1[j] + 2.0 * fl1_fxn * tg_yyzz_xxxyzzz_1[j];

                    tg_xyyzz_xxxxzzzz_0[j] = pb_x * tg_yyzz_xxxxzzzz_0[j] + fr * tg_yyzz_xxxxzzzz_1[j] + 2.0 * fl1_fxn * tg_yyzz_xxxzzzz_1[j];

                    tg_xyyzz_xxxyyyyy_0[j] = pb_x * tg_yyzz_xxxyyyyy_0[j] + fr * tg_yyzz_xxxyyyyy_1[j] + 1.5 * fl1_fxn * tg_yyzz_xxyyyyy_1[j];

                    tg_xyyzz_xxxyyyyz_0[j] = pb_x * tg_yyzz_xxxyyyyz_0[j] + fr * tg_yyzz_xxxyyyyz_1[j] + 1.5 * fl1_fxn * tg_yyzz_xxyyyyz_1[j];

                    tg_xyyzz_xxxyyyzz_0[j] = pb_x * tg_yyzz_xxxyyyzz_0[j] + fr * tg_yyzz_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_yyzz_xxyyyzz_1[j];

                    tg_xyyzz_xxxyyzzz_0[j] = pb_x * tg_yyzz_xxxyyzzz_0[j] + fr * tg_yyzz_xxxyyzzz_1[j] + 1.5 * fl1_fxn * tg_yyzz_xxyyzzz_1[j];

                    tg_xyyzz_xxxyzzzz_0[j] = pb_x * tg_yyzz_xxxyzzzz_0[j] + fr * tg_yyzz_xxxyzzzz_1[j] + 1.5 * fl1_fxn * tg_yyzz_xxyzzzz_1[j];

                    tg_xyyzz_xxxzzzzz_0[j] = pb_x * tg_yyzz_xxxzzzzz_0[j] + fr * tg_yyzz_xxxzzzzz_1[j] + 1.5 * fl1_fxn * tg_yyzz_xxzzzzz_1[j];

                    tg_xyyzz_xxyyyyyy_0[j] = pb_x * tg_yyzz_xxyyyyyy_0[j] + fr * tg_yyzz_xxyyyyyy_1[j] + fl1_fxn * tg_yyzz_xyyyyyy_1[j];

                    tg_xyyzz_xxyyyyyz_0[j] = pb_x * tg_yyzz_xxyyyyyz_0[j] + fr * tg_yyzz_xxyyyyyz_1[j] + fl1_fxn * tg_yyzz_xyyyyyz_1[j];

                    tg_xyyzz_xxyyyyzz_0[j] = pb_x * tg_yyzz_xxyyyyzz_0[j] + fr * tg_yyzz_xxyyyyzz_1[j] + fl1_fxn * tg_yyzz_xyyyyzz_1[j];

                    tg_xyyzz_xxyyyzzz_0[j] = pb_x * tg_yyzz_xxyyyzzz_0[j] + fr * tg_yyzz_xxyyyzzz_1[j] + fl1_fxn * tg_yyzz_xyyyzzz_1[j];

                    tg_xyyzz_xxyyzzzz_0[j] = pb_x * tg_yyzz_xxyyzzzz_0[j] + fr * tg_yyzz_xxyyzzzz_1[j] + fl1_fxn * tg_yyzz_xyyzzzz_1[j];

                    tg_xyyzz_xxyzzzzz_0[j] = pb_x * tg_yyzz_xxyzzzzz_0[j] + fr * tg_yyzz_xxyzzzzz_1[j] + fl1_fxn * tg_yyzz_xyzzzzz_1[j];

                    tg_xyyzz_xxzzzzzz_0[j] = pb_x * tg_yyzz_xxzzzzzz_0[j] + fr * tg_yyzz_xxzzzzzz_1[j] + fl1_fxn * tg_yyzz_xzzzzzz_1[j];

                    tg_xyyzz_xyyyyyyy_0[j] = pb_x * tg_yyzz_xyyyyyyy_0[j] + fr * tg_yyzz_xyyyyyyy_1[j] + 0.5 * fl1_fxn * tg_yyzz_yyyyyyy_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSL_569_663(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (569,663)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {5, -1, -1, -1},
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_8_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_8_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_8_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_7_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {7, -1, -1, -1}, 
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

                auto tg_yyzz_xyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 569); 

                auto tg_yyzz_xyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 570); 

                auto tg_yyzz_xyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 571); 

                auto tg_yyzz_xyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 572); 

                auto tg_yyzz_xyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 573); 

                auto tg_yyzz_xyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 574); 

                auto tg_yyzz_xzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 575); 

                auto tg_yyzz_yyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 576); 

                auto tg_yyzz_yyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 577); 

                auto tg_yyzz_yyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 578); 

                auto tg_yyzz_yyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 579); 

                auto tg_yyzz_yyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 580); 

                auto tg_yyzz_yyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 581); 

                auto tg_yyzz_yyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 582); 

                auto tg_yyzz_yzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 583); 

                auto tg_yyzz_zzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 584); 

                auto tg_yzzz_xxxxxxxx_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 585); 

                auto tg_yzzz_xxxxxxxy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 586); 

                auto tg_yzzz_xxxxxxxz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 587); 

                auto tg_yzzz_xxxxxxyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 588); 

                auto tg_yzzz_xxxxxxyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 589); 

                auto tg_yzzz_xxxxxxzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 590); 

                auto tg_yzzz_xxxxxyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 591); 

                auto tg_yzzz_xxxxxyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 592); 

                auto tg_yzzz_xxxxxyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 593); 

                auto tg_yzzz_xxxxxzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 594); 

                auto tg_yzzz_xxxxyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 595); 

                auto tg_yzzz_xxxxyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 596); 

                auto tg_yzzz_xxxxyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 597); 

                auto tg_yzzz_xxxxyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 598); 

                auto tg_yzzz_xxxxzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 599); 

                auto tg_yzzz_xxxyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 600); 

                auto tg_yzzz_xxxyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 601); 

                auto tg_yzzz_xxxyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 602); 

                auto tg_yzzz_xxxyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 603); 

                auto tg_yzzz_xxxyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 604); 

                auto tg_yzzz_xxxzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 605); 

                auto tg_yzzz_xxyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 606); 

                auto tg_yzzz_xxyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 607); 

                auto tg_yzzz_xxyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 608); 

                auto tg_yzzz_xxyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 609); 

                auto tg_yzzz_xxyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 610); 

                auto tg_yzzz_xxyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 611); 

                auto tg_yzzz_xxzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 612); 

                auto tg_yzzz_xyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 613); 

                auto tg_yzzz_xyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 614); 

                auto tg_yzzz_xyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 615); 

                auto tg_yzzz_xyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 616); 

                auto tg_yzzz_xyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 617); 

                auto tg_yzzz_xyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 618); 

                auto tg_yzzz_xyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 619); 

                auto tg_yzzz_xzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 620); 

                auto tg_yzzz_yyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 621); 

                auto tg_yzzz_yyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 622); 

                auto tg_yzzz_yyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 623); 

                auto tg_yzzz_yyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 624); 

                auto tg_yzzz_yyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 625); 

                auto tg_yzzz_yyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 626); 

                auto tg_yzzz_yyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 627); 

                auto tg_yzzz_yzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 628); 

                auto tg_yzzz_zzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 629); 

                auto tg_zzzz_xxxxxxxx_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 630); 

                auto tg_zzzz_xxxxxxxy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 631); 

                auto tg_zzzz_xxxxxxxz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 632); 

                auto tg_zzzz_xxxxxxyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 633); 

                auto tg_zzzz_xxxxxxyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 634); 

                auto tg_zzzz_xxxxxxzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 635); 

                auto tg_zzzz_xxxxxyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 636); 

                auto tg_zzzz_xxxxxyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 637); 

                auto tg_zzzz_xxxxxyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 638); 

                auto tg_zzzz_xxxxxzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 639); 

                auto tg_zzzz_xxxxyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 640); 

                auto tg_zzzz_xxxxyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 641); 

                auto tg_zzzz_xxxxyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 642); 

                auto tg_zzzz_xxxxyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 643); 

                auto tg_zzzz_xxxxzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 644); 

                auto tg_zzzz_xxxyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 645); 

                auto tg_zzzz_xxxyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 646); 

                auto tg_zzzz_xxxyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 647); 

                auto tg_zzzz_xxxyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 648); 

                auto tg_zzzz_xxxyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 649); 

                auto tg_zzzz_xxxzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 650); 

                auto tg_zzzz_xxyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 651); 

                auto tg_zzzz_xxyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 652); 

                auto tg_zzzz_xxyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 653); 

                auto tg_zzzz_xxyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 654); 

                auto tg_zzzz_xxyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 655); 

                auto tg_zzzz_xxyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 656); 

                auto tg_zzzz_xxzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 657); 

                auto tg_zzzz_xyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 658); 

                auto tg_zzzz_xyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 659); 

                auto tg_zzzz_xyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 660); 

                auto tg_zzzz_xyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 661); 

                auto tg_zzzz_xyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 662); 

                auto tg_yyzz_xyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 569); 

                auto tg_yyzz_xyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 570); 

                auto tg_yyzz_xyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 571); 

                auto tg_yyzz_xyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 572); 

                auto tg_yyzz_xyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 573); 

                auto tg_yyzz_xyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 574); 

                auto tg_yyzz_xzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 575); 

                auto tg_yyzz_yyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 576); 

                auto tg_yyzz_yyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 577); 

                auto tg_yyzz_yyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 578); 

                auto tg_yyzz_yyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 579); 

                auto tg_yyzz_yyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 580); 

                auto tg_yyzz_yyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 581); 

                auto tg_yyzz_yyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 582); 

                auto tg_yyzz_yzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 583); 

                auto tg_yyzz_zzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 584); 

                auto tg_yzzz_xxxxxxxx_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 585); 

                auto tg_yzzz_xxxxxxxy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 586); 

                auto tg_yzzz_xxxxxxxz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 587); 

                auto tg_yzzz_xxxxxxyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 588); 

                auto tg_yzzz_xxxxxxyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 589); 

                auto tg_yzzz_xxxxxxzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 590); 

                auto tg_yzzz_xxxxxyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 591); 

                auto tg_yzzz_xxxxxyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 592); 

                auto tg_yzzz_xxxxxyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 593); 

                auto tg_yzzz_xxxxxzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 594); 

                auto tg_yzzz_xxxxyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 595); 

                auto tg_yzzz_xxxxyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 596); 

                auto tg_yzzz_xxxxyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 597); 

                auto tg_yzzz_xxxxyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 598); 

                auto tg_yzzz_xxxxzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 599); 

                auto tg_yzzz_xxxyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 600); 

                auto tg_yzzz_xxxyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 601); 

                auto tg_yzzz_xxxyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 602); 

                auto tg_yzzz_xxxyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 603); 

                auto tg_yzzz_xxxyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 604); 

                auto tg_yzzz_xxxzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 605); 

                auto tg_yzzz_xxyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 606); 

                auto tg_yzzz_xxyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 607); 

                auto tg_yzzz_xxyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 608); 

                auto tg_yzzz_xxyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 609); 

                auto tg_yzzz_xxyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 610); 

                auto tg_yzzz_xxyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 611); 

                auto tg_yzzz_xxzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 612); 

                auto tg_yzzz_xyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 613); 

                auto tg_yzzz_xyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 614); 

                auto tg_yzzz_xyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 615); 

                auto tg_yzzz_xyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 616); 

                auto tg_yzzz_xyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 617); 

                auto tg_yzzz_xyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 618); 

                auto tg_yzzz_xyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 619); 

                auto tg_yzzz_xzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 620); 

                auto tg_yzzz_yyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 621); 

                auto tg_yzzz_yyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 622); 

                auto tg_yzzz_yyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 623); 

                auto tg_yzzz_yyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 624); 

                auto tg_yzzz_yyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 625); 

                auto tg_yzzz_yyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 626); 

                auto tg_yzzz_yyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 627); 

                auto tg_yzzz_yzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 628); 

                auto tg_yzzz_zzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 629); 

                auto tg_zzzz_xxxxxxxx_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 630); 

                auto tg_zzzz_xxxxxxxy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 631); 

                auto tg_zzzz_xxxxxxxz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 632); 

                auto tg_zzzz_xxxxxxyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 633); 

                auto tg_zzzz_xxxxxxyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 634); 

                auto tg_zzzz_xxxxxxzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 635); 

                auto tg_zzzz_xxxxxyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 636); 

                auto tg_zzzz_xxxxxyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 637); 

                auto tg_zzzz_xxxxxyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 638); 

                auto tg_zzzz_xxxxxzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 639); 

                auto tg_zzzz_xxxxyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 640); 

                auto tg_zzzz_xxxxyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 641); 

                auto tg_zzzz_xxxxyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 642); 

                auto tg_zzzz_xxxxyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 643); 

                auto tg_zzzz_xxxxzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 644); 

                auto tg_zzzz_xxxyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 645); 

                auto tg_zzzz_xxxyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 646); 

                auto tg_zzzz_xxxyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 647); 

                auto tg_zzzz_xxxyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 648); 

                auto tg_zzzz_xxxyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 649); 

                auto tg_zzzz_xxxzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 650); 

                auto tg_zzzz_xxyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 651); 

                auto tg_zzzz_xxyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 652); 

                auto tg_zzzz_xxyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 653); 

                auto tg_zzzz_xxyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 654); 

                auto tg_zzzz_xxyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 655); 

                auto tg_zzzz_xxyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 656); 

                auto tg_zzzz_xxzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 657); 

                auto tg_zzzz_xyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 658); 

                auto tg_zzzz_xyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 659); 

                auto tg_zzzz_xyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 660); 

                auto tg_zzzz_xyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 661); 

                auto tg_zzzz_xyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 662); 

                auto tg_yyzz_yyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 461); 

                auto tg_yyzz_yyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 462); 

                auto tg_yyzz_yyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 463); 

                auto tg_yyzz_yyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 464); 

                auto tg_yyzz_yyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 465); 

                auto tg_yyzz_yzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 466); 

                auto tg_yyzz_zzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 467); 

                auto tg_yzzz_xxxxxxx_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 468); 

                auto tg_yzzz_xxxxxxy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 469); 

                auto tg_yzzz_xxxxxxz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 470); 

                auto tg_yzzz_xxxxxyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 471); 

                auto tg_yzzz_xxxxxyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 472); 

                auto tg_yzzz_xxxxxzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 473); 

                auto tg_yzzz_xxxxyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 474); 

                auto tg_yzzz_xxxxyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 475); 

                auto tg_yzzz_xxxxyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 476); 

                auto tg_yzzz_xxxxzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 477); 

                auto tg_yzzz_xxxyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 478); 

                auto tg_yzzz_xxxyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 479); 

                auto tg_yzzz_xxxyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 480); 

                auto tg_yzzz_xxxyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 481); 

                auto tg_yzzz_xxxzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 482); 

                auto tg_yzzz_xxyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 483); 

                auto tg_yzzz_xxyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 484); 

                auto tg_yzzz_xxyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 485); 

                auto tg_yzzz_xxyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 486); 

                auto tg_yzzz_xxyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 487); 

                auto tg_yzzz_xxzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 488); 

                auto tg_yzzz_xyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 489); 

                auto tg_yzzz_xyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 490); 

                auto tg_yzzz_xyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 491); 

                auto tg_yzzz_xyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 492); 

                auto tg_yzzz_xyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 493); 

                auto tg_yzzz_xyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 494); 

                auto tg_yzzz_xzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 495); 

                auto tg_yzzz_yyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 496); 

                auto tg_yzzz_yyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 497); 

                auto tg_yzzz_yyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 498); 

                auto tg_yzzz_yyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 499); 

                auto tg_yzzz_yyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 500); 

                auto tg_yzzz_yyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 501); 

                auto tg_yzzz_yzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 502); 

                auto tg_yzzz_zzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 503); 

                auto tg_zzzz_xxxxxxx_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 504); 

                auto tg_zzzz_xxxxxxy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 505); 

                auto tg_zzzz_xxxxxxz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 506); 

                auto tg_zzzz_xxxxxyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 507); 

                auto tg_zzzz_xxxxxyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 508); 

                auto tg_zzzz_xxxxxzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 509); 

                auto tg_zzzz_xxxxyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 510); 

                auto tg_zzzz_xxxxyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 511); 

                auto tg_zzzz_xxxxyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 512); 

                auto tg_zzzz_xxxxzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 513); 

                auto tg_zzzz_xxxyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 514); 

                auto tg_zzzz_xxxyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 515); 

                auto tg_zzzz_xxxyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 516); 

                auto tg_zzzz_xxxyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 517); 

                auto tg_zzzz_xxxzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 518); 

                auto tg_zzzz_xxyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 519); 

                auto tg_zzzz_xxyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 520); 

                auto tg_zzzz_xxyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 521); 

                auto tg_zzzz_xxyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 522); 

                auto tg_zzzz_xxyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 523); 

                auto tg_zzzz_xxzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 524); 

                auto tg_zzzz_xyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 525); 

                auto tg_zzzz_xyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 526); 

                auto tg_zzzz_xyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 527); 

                auto tg_zzzz_xyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 528); 

                auto tg_zzzz_xyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 529); 

                auto tg_zzzz_xyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 530); 

                auto tg_zzzz_xzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 531); 

                auto tg_zzzz_yyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 532); 

                auto tg_zzzz_yyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 533); 

                auto tg_zzzz_yyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 534); 

                auto tg_zzzz_yyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 535); 

                auto tg_zzzz_yyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 536); 

                // set up pointers to integrals

                auto tg_xyyzz_xyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 569); 

                auto tg_xyyzz_xyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 570); 

                auto tg_xyyzz_xyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 571); 

                auto tg_xyyzz_xyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 572); 

                auto tg_xyyzz_xyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 573); 

                auto tg_xyyzz_xyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 574); 

                auto tg_xyyzz_xzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 575); 

                auto tg_xyyzz_yyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 576); 

                auto tg_xyyzz_yyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 577); 

                auto tg_xyyzz_yyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 578); 

                auto tg_xyyzz_yyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 579); 

                auto tg_xyyzz_yyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 580); 

                auto tg_xyyzz_yyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 581); 

                auto tg_xyyzz_yyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 582); 

                auto tg_xyyzz_yzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 583); 

                auto tg_xyyzz_zzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 584); 

                auto tg_xyzzz_xxxxxxxx_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 585); 

                auto tg_xyzzz_xxxxxxxy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 586); 

                auto tg_xyzzz_xxxxxxxz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 587); 

                auto tg_xyzzz_xxxxxxyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 588); 

                auto tg_xyzzz_xxxxxxyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 589); 

                auto tg_xyzzz_xxxxxxzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 590); 

                auto tg_xyzzz_xxxxxyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 591); 

                auto tg_xyzzz_xxxxxyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 592); 

                auto tg_xyzzz_xxxxxyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 593); 

                auto tg_xyzzz_xxxxxzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 594); 

                auto tg_xyzzz_xxxxyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 595); 

                auto tg_xyzzz_xxxxyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 596); 

                auto tg_xyzzz_xxxxyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 597); 

                auto tg_xyzzz_xxxxyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 598); 

                auto tg_xyzzz_xxxxzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 599); 

                auto tg_xyzzz_xxxyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 600); 

                auto tg_xyzzz_xxxyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 601); 

                auto tg_xyzzz_xxxyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 602); 

                auto tg_xyzzz_xxxyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 603); 

                auto tg_xyzzz_xxxyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 604); 

                auto tg_xyzzz_xxxzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 605); 

                auto tg_xyzzz_xxyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 606); 

                auto tg_xyzzz_xxyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 607); 

                auto tg_xyzzz_xxyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 608); 

                auto tg_xyzzz_xxyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 609); 

                auto tg_xyzzz_xxyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 610); 

                auto tg_xyzzz_xxyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 611); 

                auto tg_xyzzz_xxzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 612); 

                auto tg_xyzzz_xyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 613); 

                auto tg_xyzzz_xyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 614); 

                auto tg_xyzzz_xyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 615); 

                auto tg_xyzzz_xyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 616); 

                auto tg_xyzzz_xyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 617); 

                auto tg_xyzzz_xyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 618); 

                auto tg_xyzzz_xyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 619); 

                auto tg_xyzzz_xzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 620); 

                auto tg_xyzzz_yyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 621); 

                auto tg_xyzzz_yyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 622); 

                auto tg_xyzzz_yyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 623); 

                auto tg_xyzzz_yyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 624); 

                auto tg_xyzzz_yyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 625); 

                auto tg_xyzzz_yyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 626); 

                auto tg_xyzzz_yyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 627); 

                auto tg_xyzzz_yzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 628); 

                auto tg_xyzzz_zzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 629); 

                auto tg_xzzzz_xxxxxxxx_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 630); 

                auto tg_xzzzz_xxxxxxxy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 631); 

                auto tg_xzzzz_xxxxxxxz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 632); 

                auto tg_xzzzz_xxxxxxyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 633); 

                auto tg_xzzzz_xxxxxxyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 634); 

                auto tg_xzzzz_xxxxxxzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 635); 

                auto tg_xzzzz_xxxxxyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 636); 

                auto tg_xzzzz_xxxxxyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 637); 

                auto tg_xzzzz_xxxxxyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 638); 

                auto tg_xzzzz_xxxxxzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 639); 

                auto tg_xzzzz_xxxxyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 640); 

                auto tg_xzzzz_xxxxyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 641); 

                auto tg_xzzzz_xxxxyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 642); 

                auto tg_xzzzz_xxxxyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 643); 

                auto tg_xzzzz_xxxxzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 644); 

                auto tg_xzzzz_xxxyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 645); 

                auto tg_xzzzz_xxxyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 646); 

                auto tg_xzzzz_xxxyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 647); 

                auto tg_xzzzz_xxxyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 648); 

                auto tg_xzzzz_xxxyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 649); 

                auto tg_xzzzz_xxxzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 650); 

                auto tg_xzzzz_xxyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 651); 

                auto tg_xzzzz_xxyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 652); 

                auto tg_xzzzz_xxyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 653); 

                auto tg_xzzzz_xxyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 654); 

                auto tg_xzzzz_xxyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 655); 

                auto tg_xzzzz_xxyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 656); 

                auto tg_xzzzz_xxzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 657); 

                auto tg_xzzzz_xyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 658); 

                auto tg_xzzzz_xyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 659); 

                auto tg_xzzzz_xyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 660); 

                auto tg_xzzzz_xyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 661); 

                auto tg_xzzzz_xyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 662); 

                // Batch of Integrals (569,663)

                #pragma omp simd aligned(fxn, tg_xyyzz_xyyyyyyz_0, tg_xyyzz_xyyyyyzz_0, tg_xyyzz_xyyyyzzz_0, \
                                         tg_xyyzz_xyyyzzzz_0, tg_xyyzz_xyyzzzzz_0, tg_xyyzz_xyzzzzzz_0, tg_xyyzz_xzzzzzzz_0, \
                                         tg_xyyzz_yyyyyyyy_0, tg_xyyzz_yyyyyyyz_0, tg_xyyzz_yyyyyyzz_0, tg_xyyzz_yyyyyzzz_0, \
                                         tg_xyyzz_yyyyzzzz_0, tg_xyyzz_yyyzzzzz_0, tg_xyyzz_yyzzzzzz_0, tg_xyyzz_yzzzzzzz_0, \
                                         tg_xyyzz_zzzzzzzz_0, tg_xyzzz_xxxxxxxx_0, tg_xyzzz_xxxxxxxy_0, tg_xyzzz_xxxxxxxz_0, \
                                         tg_xyzzz_xxxxxxyy_0, tg_xyzzz_xxxxxxyz_0, tg_xyzzz_xxxxxxzz_0, tg_xyzzz_xxxxxyyy_0, \
                                         tg_xyzzz_xxxxxyyz_0, tg_xyzzz_xxxxxyzz_0, tg_xyzzz_xxxxxzzz_0, tg_xyzzz_xxxxyyyy_0, \
                                         tg_xyzzz_xxxxyyyz_0, tg_xyzzz_xxxxyyzz_0, tg_xyzzz_xxxxyzzz_0, tg_xyzzz_xxxxzzzz_0, \
                                         tg_xyzzz_xxxyyyyy_0, tg_xyzzz_xxxyyyyz_0, tg_xyzzz_xxxyyyzz_0, tg_xyzzz_xxxyyzzz_0, \
                                         tg_xyzzz_xxxyzzzz_0, tg_xyzzz_xxxzzzzz_0, tg_xyzzz_xxyyyyyy_0, tg_xyzzz_xxyyyyyz_0, \
                                         tg_xyzzz_xxyyyyzz_0, tg_xyzzz_xxyyyzzz_0, tg_xyzzz_xxyyzzzz_0, tg_xyzzz_xxyzzzzz_0, \
                                         tg_xyzzz_xxzzzzzz_0, tg_xyzzz_xyyyyyyy_0, tg_xyzzz_xyyyyyyz_0, tg_xyzzz_xyyyyyzz_0, \
                                         tg_xyzzz_xyyyyzzz_0, tg_xyzzz_xyyyzzzz_0, tg_xyzzz_xyyzzzzz_0, tg_xyzzz_xyzzzzzz_0, \
                                         tg_xyzzz_xzzzzzzz_0, tg_xyzzz_yyyyyyyy_0, tg_xyzzz_yyyyyyyz_0, tg_xyzzz_yyyyyyzz_0, \
                                         tg_xyzzz_yyyyyzzz_0, tg_xyzzz_yyyyzzzz_0, tg_xyzzz_yyyzzzzz_0, tg_xyzzz_yyzzzzzz_0, \
                                         tg_xyzzz_yzzzzzzz_0, tg_xyzzz_zzzzzzzz_0, tg_xzzzz_xxxxxxxx_0, tg_xzzzz_xxxxxxxy_0, \
                                         tg_xzzzz_xxxxxxxz_0, tg_xzzzz_xxxxxxyy_0, tg_xzzzz_xxxxxxyz_0, tg_xzzzz_xxxxxxzz_0, \
                                         tg_xzzzz_xxxxxyyy_0, tg_xzzzz_xxxxxyyz_0, tg_xzzzz_xxxxxyzz_0, tg_xzzzz_xxxxxzzz_0, \
                                         tg_xzzzz_xxxxyyyy_0, tg_xzzzz_xxxxyyyz_0, tg_xzzzz_xxxxyyzz_0, tg_xzzzz_xxxxyzzz_0, \
                                         tg_xzzzz_xxxxzzzz_0, tg_xzzzz_xxxyyyyy_0, tg_xzzzz_xxxyyyyz_0, tg_xzzzz_xxxyyyzz_0, \
                                         tg_xzzzz_xxxyyzzz_0, tg_xzzzz_xxxyzzzz_0, tg_xzzzz_xxxzzzzz_0, tg_xzzzz_xxyyyyyy_0, \
                                         tg_xzzzz_xxyyyyyz_0, tg_xzzzz_xxyyyyzz_0, tg_xzzzz_xxyyyzzz_0, tg_xzzzz_xxyyzzzz_0, \
                                         tg_xzzzz_xxyzzzzz_0, tg_xzzzz_xxzzzzzz_0, tg_xzzzz_xyyyyyyy_0, tg_xzzzz_xyyyyyyz_0, \
                                         tg_xzzzz_xyyyyyzz_0, tg_xzzzz_xyyyyzzz_0, tg_xzzzz_xyyyzzzz_0, tg_yyzz_xyyyyyyz_0, \
                                         tg_yyzz_xyyyyyyz_1, tg_yyzz_xyyyyyzz_0, tg_yyzz_xyyyyyzz_1, tg_yyzz_xyyyyzzz_0, \
                                         tg_yyzz_xyyyyzzz_1, tg_yyzz_xyyyzzzz_0, tg_yyzz_xyyyzzzz_1, tg_yyzz_xyyzzzzz_0, \
                                         tg_yyzz_xyyzzzzz_1, tg_yyzz_xyzzzzzz_0, tg_yyzz_xyzzzzzz_1, tg_yyzz_xzzzzzzz_0, \
                                         tg_yyzz_xzzzzzzz_1, tg_yyzz_yyyyyyyy_0, tg_yyzz_yyyyyyyy_1, tg_yyzz_yyyyyyyz_0, \
                                         tg_yyzz_yyyyyyyz_1, tg_yyzz_yyyyyyz_1, tg_yyzz_yyyyyyzz_0, tg_yyzz_yyyyyyzz_1, \
                                         tg_yyzz_yyyyyzz_1, tg_yyzz_yyyyyzzz_0, tg_yyzz_yyyyyzzz_1, tg_yyzz_yyyyzzz_1, \
                                         tg_yyzz_yyyyzzzz_0, tg_yyzz_yyyyzzzz_1, tg_yyzz_yyyzzzz_1, tg_yyzz_yyyzzzzz_0, \
                                         tg_yyzz_yyyzzzzz_1, tg_yyzz_yyzzzzz_1, tg_yyzz_yyzzzzzz_0, tg_yyzz_yyzzzzzz_1, \
                                         tg_yyzz_yzzzzzz_1, tg_yyzz_yzzzzzzz_0, tg_yyzz_yzzzzzzz_1, tg_yyzz_zzzzzzz_1, \
                                         tg_yyzz_zzzzzzzz_0, tg_yyzz_zzzzzzzz_1, tg_yzzz_xxxxxxx_1, tg_yzzz_xxxxxxxx_0, \
                                         tg_yzzz_xxxxxxxx_1, tg_yzzz_xxxxxxxy_0, tg_yzzz_xxxxxxxy_1, tg_yzzz_xxxxxxxz_0, \
                                         tg_yzzz_xxxxxxxz_1, tg_yzzz_xxxxxxy_1, tg_yzzz_xxxxxxyy_0, tg_yzzz_xxxxxxyy_1, \
                                         tg_yzzz_xxxxxxyz_0, tg_yzzz_xxxxxxyz_1, tg_yzzz_xxxxxxz_1, tg_yzzz_xxxxxxzz_0, \
                                         tg_yzzz_xxxxxxzz_1, tg_yzzz_xxxxxyy_1, tg_yzzz_xxxxxyyy_0, tg_yzzz_xxxxxyyy_1, \
                                         tg_yzzz_xxxxxyyz_0, tg_yzzz_xxxxxyyz_1, tg_yzzz_xxxxxyz_1, tg_yzzz_xxxxxyzz_0, \
                                         tg_yzzz_xxxxxyzz_1, tg_yzzz_xxxxxzz_1, tg_yzzz_xxxxxzzz_0, tg_yzzz_xxxxxzzz_1, \
                                         tg_yzzz_xxxxyyy_1, tg_yzzz_xxxxyyyy_0, tg_yzzz_xxxxyyyy_1, tg_yzzz_xxxxyyyz_0, \
                                         tg_yzzz_xxxxyyyz_1, tg_yzzz_xxxxyyz_1, tg_yzzz_xxxxyyzz_0, tg_yzzz_xxxxyyzz_1, \
                                         tg_yzzz_xxxxyzz_1, tg_yzzz_xxxxyzzz_0, tg_yzzz_xxxxyzzz_1, tg_yzzz_xxxxzzz_1, \
                                         tg_yzzz_xxxxzzzz_0, tg_yzzz_xxxxzzzz_1, tg_yzzz_xxxyyyy_1, tg_yzzz_xxxyyyyy_0, \
                                         tg_yzzz_xxxyyyyy_1, tg_yzzz_xxxyyyyz_0, tg_yzzz_xxxyyyyz_1, tg_yzzz_xxxyyyz_1, \
                                         tg_yzzz_xxxyyyzz_0, tg_yzzz_xxxyyyzz_1, tg_yzzz_xxxyyzz_1, tg_yzzz_xxxyyzzz_0, \
                                         tg_yzzz_xxxyyzzz_1, tg_yzzz_xxxyzzz_1, tg_yzzz_xxxyzzzz_0, tg_yzzz_xxxyzzzz_1, \
                                         tg_yzzz_xxxzzzz_1, tg_yzzz_xxxzzzzz_0, tg_yzzz_xxxzzzzz_1, tg_yzzz_xxyyyyy_1, \
                                         tg_yzzz_xxyyyyyy_0, tg_yzzz_xxyyyyyy_1, tg_yzzz_xxyyyyyz_0, tg_yzzz_xxyyyyyz_1, \
                                         tg_yzzz_xxyyyyz_1, tg_yzzz_xxyyyyzz_0, tg_yzzz_xxyyyyzz_1, tg_yzzz_xxyyyzz_1, \
                                         tg_yzzz_xxyyyzzz_0, tg_yzzz_xxyyyzzz_1, tg_yzzz_xxyyzzz_1, tg_yzzz_xxyyzzzz_0, \
                                         tg_yzzz_xxyyzzzz_1, tg_yzzz_xxyzzzz_1, tg_yzzz_xxyzzzzz_0, tg_yzzz_xxyzzzzz_1, \
                                         tg_yzzz_xxzzzzz_1, tg_yzzz_xxzzzzzz_0, tg_yzzz_xxzzzzzz_1, tg_yzzz_xyyyyyy_1, \
                                         tg_yzzz_xyyyyyyy_0, tg_yzzz_xyyyyyyy_1, tg_yzzz_xyyyyyyz_0, tg_yzzz_xyyyyyyz_1, \
                                         tg_yzzz_xyyyyyz_1, tg_yzzz_xyyyyyzz_0, tg_yzzz_xyyyyyzz_1, tg_yzzz_xyyyyzz_1, \
                                         tg_yzzz_xyyyyzzz_0, tg_yzzz_xyyyyzzz_1, tg_yzzz_xyyyzzz_1, tg_yzzz_xyyyzzzz_0, \
                                         tg_yzzz_xyyyzzzz_1, tg_yzzz_xyyzzzz_1, tg_yzzz_xyyzzzzz_0, tg_yzzz_xyyzzzzz_1, \
                                         tg_yzzz_xyzzzzz_1, tg_yzzz_xyzzzzzz_0, tg_yzzz_xyzzzzzz_1, tg_yzzz_xzzzzzz_1, \
                                         tg_yzzz_xzzzzzzz_0, tg_yzzz_xzzzzzzz_1, tg_yzzz_yyyyyyy_1, tg_yzzz_yyyyyyyy_0, \
                                         tg_yzzz_yyyyyyyy_1, tg_yzzz_yyyyyyyz_0, tg_yzzz_yyyyyyyz_1, tg_yzzz_yyyyyyz_1, \
                                         tg_yzzz_yyyyyyzz_0, tg_yzzz_yyyyyyzz_1, tg_yzzz_yyyyyzz_1, tg_yzzz_yyyyyzzz_0, \
                                         tg_yzzz_yyyyyzzz_1, tg_yzzz_yyyyzzz_1, tg_yzzz_yyyyzzzz_0, tg_yzzz_yyyyzzzz_1, \
                                         tg_yzzz_yyyzzzz_1, tg_yzzz_yyyzzzzz_0, tg_yzzz_yyyzzzzz_1, tg_yzzz_yyzzzzz_1, \
                                         tg_yzzz_yyzzzzzz_0, tg_yzzz_yyzzzzzz_1, tg_yzzz_yzzzzzz_1, tg_yzzz_yzzzzzzz_0, \
                                         tg_yzzz_yzzzzzzz_1, tg_yzzz_zzzzzzz_1, tg_yzzz_zzzzzzzz_0, tg_yzzz_zzzzzzzz_1, \
                                         tg_zzzz_xxxxxxx_1, tg_zzzz_xxxxxxxx_0, tg_zzzz_xxxxxxxx_1, tg_zzzz_xxxxxxxy_0, \
                                         tg_zzzz_xxxxxxxy_1, tg_zzzz_xxxxxxxz_0, tg_zzzz_xxxxxxxz_1, tg_zzzz_xxxxxxy_1, \
                                         tg_zzzz_xxxxxxyy_0, tg_zzzz_xxxxxxyy_1, tg_zzzz_xxxxxxyz_0, tg_zzzz_xxxxxxyz_1, \
                                         tg_zzzz_xxxxxxz_1, tg_zzzz_xxxxxxzz_0, tg_zzzz_xxxxxxzz_1, tg_zzzz_xxxxxyy_1, \
                                         tg_zzzz_xxxxxyyy_0, tg_zzzz_xxxxxyyy_1, tg_zzzz_xxxxxyyz_0, tg_zzzz_xxxxxyyz_1, \
                                         tg_zzzz_xxxxxyz_1, tg_zzzz_xxxxxyzz_0, tg_zzzz_xxxxxyzz_1, tg_zzzz_xxxxxzz_1, \
                                         tg_zzzz_xxxxxzzz_0, tg_zzzz_xxxxxzzz_1, tg_zzzz_xxxxyyy_1, tg_zzzz_xxxxyyyy_0, \
                                         tg_zzzz_xxxxyyyy_1, tg_zzzz_xxxxyyyz_0, tg_zzzz_xxxxyyyz_1, tg_zzzz_xxxxyyz_1, \
                                         tg_zzzz_xxxxyyzz_0, tg_zzzz_xxxxyyzz_1, tg_zzzz_xxxxyzz_1, tg_zzzz_xxxxyzzz_0, \
                                         tg_zzzz_xxxxyzzz_1, tg_zzzz_xxxxzzz_1, tg_zzzz_xxxxzzzz_0, tg_zzzz_xxxxzzzz_1, \
                                         tg_zzzz_xxxyyyy_1, tg_zzzz_xxxyyyyy_0, tg_zzzz_xxxyyyyy_1, tg_zzzz_xxxyyyyz_0, \
                                         tg_zzzz_xxxyyyyz_1, tg_zzzz_xxxyyyz_1, tg_zzzz_xxxyyyzz_0, tg_zzzz_xxxyyyzz_1, \
                                         tg_zzzz_xxxyyzz_1, tg_zzzz_xxxyyzzz_0, tg_zzzz_xxxyyzzz_1, tg_zzzz_xxxyzzz_1, \
                                         tg_zzzz_xxxyzzzz_0, tg_zzzz_xxxyzzzz_1, tg_zzzz_xxxzzzz_1, tg_zzzz_xxxzzzzz_0, \
                                         tg_zzzz_xxxzzzzz_1, tg_zzzz_xxyyyyy_1, tg_zzzz_xxyyyyyy_0, tg_zzzz_xxyyyyyy_1, \
                                         tg_zzzz_xxyyyyyz_0, tg_zzzz_xxyyyyyz_1, tg_zzzz_xxyyyyz_1, tg_zzzz_xxyyyyzz_0, \
                                         tg_zzzz_xxyyyyzz_1, tg_zzzz_xxyyyzz_1, tg_zzzz_xxyyyzzz_0, tg_zzzz_xxyyyzzz_1, \
                                         tg_zzzz_xxyyzzz_1, tg_zzzz_xxyyzzzz_0, tg_zzzz_xxyyzzzz_1, tg_zzzz_xxyzzzz_1, \
                                         tg_zzzz_xxyzzzzz_0, tg_zzzz_xxyzzzzz_1, tg_zzzz_xxzzzzz_1, tg_zzzz_xxzzzzzz_0, \
                                         tg_zzzz_xxzzzzzz_1, tg_zzzz_xyyyyyy_1, tg_zzzz_xyyyyyyy_0, tg_zzzz_xyyyyyyy_1, \
                                         tg_zzzz_xyyyyyyz_0, tg_zzzz_xyyyyyyz_1, tg_zzzz_xyyyyyz_1, tg_zzzz_xyyyyyzz_0, \
                                         tg_zzzz_xyyyyyzz_1, tg_zzzz_xyyyyzz_1, tg_zzzz_xyyyyzzz_0, tg_zzzz_xyyyyzzz_1, \
                                         tg_zzzz_xyyyzzz_1, tg_zzzz_xyyyzzzz_0, tg_zzzz_xyyyzzzz_1, tg_zzzz_xyyzzzz_1, \
                                         tg_zzzz_xyzzzzz_1, tg_zzzz_xzzzzzz_1, tg_zzzz_yyyyyyy_1, tg_zzzz_yyyyyyz_1, \
                                         tg_zzzz_yyyyyzz_1, tg_zzzz_yyyyzzz_1, tg_zzzz_yyyzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fxn = fxn[j];

                    double fr = wp_x[j]; 

                    tg_xyyzz_xyyyyyyz_0[j] = pb_x * tg_yyzz_xyyyyyyz_0[j] + fr * tg_yyzz_xyyyyyyz_1[j] + 0.5 * fl1_fxn * tg_yyzz_yyyyyyz_1[j];

                    tg_xyyzz_xyyyyyzz_0[j] = pb_x * tg_yyzz_xyyyyyzz_0[j] + fr * tg_yyzz_xyyyyyzz_1[j] + 0.5 * fl1_fxn * tg_yyzz_yyyyyzz_1[j];

                    tg_xyyzz_xyyyyzzz_0[j] = pb_x * tg_yyzz_xyyyyzzz_0[j] + fr * tg_yyzz_xyyyyzzz_1[j] + 0.5 * fl1_fxn * tg_yyzz_yyyyzzz_1[j];

                    tg_xyyzz_xyyyzzzz_0[j] = pb_x * tg_yyzz_xyyyzzzz_0[j] + fr * tg_yyzz_xyyyzzzz_1[j] + 0.5 * fl1_fxn * tg_yyzz_yyyzzzz_1[j];

                    tg_xyyzz_xyyzzzzz_0[j] = pb_x * tg_yyzz_xyyzzzzz_0[j] + fr * tg_yyzz_xyyzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyzz_yyzzzzz_1[j];

                    tg_xyyzz_xyzzzzzz_0[j] = pb_x * tg_yyzz_xyzzzzzz_0[j] + fr * tg_yyzz_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyzz_yzzzzzz_1[j];

                    tg_xyyzz_xzzzzzzz_0[j] = pb_x * tg_yyzz_xzzzzzzz_0[j] + fr * tg_yyzz_xzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyzz_zzzzzzz_1[j];

                    tg_xyyzz_yyyyyyyy_0[j] = pb_x * tg_yyzz_yyyyyyyy_0[j] + fr * tg_yyzz_yyyyyyyy_1[j];

                    tg_xyyzz_yyyyyyyz_0[j] = pb_x * tg_yyzz_yyyyyyyz_0[j] + fr * tg_yyzz_yyyyyyyz_1[j];

                    tg_xyyzz_yyyyyyzz_0[j] = pb_x * tg_yyzz_yyyyyyzz_0[j] + fr * tg_yyzz_yyyyyyzz_1[j];

                    tg_xyyzz_yyyyyzzz_0[j] = pb_x * tg_yyzz_yyyyyzzz_0[j] + fr * tg_yyzz_yyyyyzzz_1[j];

                    tg_xyyzz_yyyyzzzz_0[j] = pb_x * tg_yyzz_yyyyzzzz_0[j] + fr * tg_yyzz_yyyyzzzz_1[j];

                    tg_xyyzz_yyyzzzzz_0[j] = pb_x * tg_yyzz_yyyzzzzz_0[j] + fr * tg_yyzz_yyyzzzzz_1[j];

                    tg_xyyzz_yyzzzzzz_0[j] = pb_x * tg_yyzz_yyzzzzzz_0[j] + fr * tg_yyzz_yyzzzzzz_1[j];

                    tg_xyyzz_yzzzzzzz_0[j] = pb_x * tg_yyzz_yzzzzzzz_0[j] + fr * tg_yyzz_yzzzzzzz_1[j];

                    tg_xyyzz_zzzzzzzz_0[j] = pb_x * tg_yyzz_zzzzzzzz_0[j] + fr * tg_yyzz_zzzzzzzz_1[j];

                    tg_xyzzz_xxxxxxxx_0[j] = pb_x * tg_yzzz_xxxxxxxx_0[j] + fr * tg_yzzz_xxxxxxxx_1[j] + 4.0 * fl1_fxn * tg_yzzz_xxxxxxx_1[j];

                    tg_xyzzz_xxxxxxxy_0[j] = pb_x * tg_yzzz_xxxxxxxy_0[j] + fr * tg_yzzz_xxxxxxxy_1[j] + 3.5 * fl1_fxn * tg_yzzz_xxxxxxy_1[j];

                    tg_xyzzz_xxxxxxxz_0[j] = pb_x * tg_yzzz_xxxxxxxz_0[j] + fr * tg_yzzz_xxxxxxxz_1[j] + 3.5 * fl1_fxn * tg_yzzz_xxxxxxz_1[j];

                    tg_xyzzz_xxxxxxyy_0[j] = pb_x * tg_yzzz_xxxxxxyy_0[j] + fr * tg_yzzz_xxxxxxyy_1[j] + 3.0 * fl1_fxn * tg_yzzz_xxxxxyy_1[j];

                    tg_xyzzz_xxxxxxyz_0[j] = pb_x * tg_yzzz_xxxxxxyz_0[j] + fr * tg_yzzz_xxxxxxyz_1[j] + 3.0 * fl1_fxn * tg_yzzz_xxxxxyz_1[j];

                    tg_xyzzz_xxxxxxzz_0[j] = pb_x * tg_yzzz_xxxxxxzz_0[j] + fr * tg_yzzz_xxxxxxzz_1[j] + 3.0 * fl1_fxn * tg_yzzz_xxxxxzz_1[j];

                    tg_xyzzz_xxxxxyyy_0[j] = pb_x * tg_yzzz_xxxxxyyy_0[j] + fr * tg_yzzz_xxxxxyyy_1[j] + 2.5 * fl1_fxn * tg_yzzz_xxxxyyy_1[j];

                    tg_xyzzz_xxxxxyyz_0[j] = pb_x * tg_yzzz_xxxxxyyz_0[j] + fr * tg_yzzz_xxxxxyyz_1[j] + 2.5 * fl1_fxn * tg_yzzz_xxxxyyz_1[j];

                    tg_xyzzz_xxxxxyzz_0[j] = pb_x * tg_yzzz_xxxxxyzz_0[j] + fr * tg_yzzz_xxxxxyzz_1[j] + 2.5 * fl1_fxn * tg_yzzz_xxxxyzz_1[j];

                    tg_xyzzz_xxxxxzzz_0[j] = pb_x * tg_yzzz_xxxxxzzz_0[j] + fr * tg_yzzz_xxxxxzzz_1[j] + 2.5 * fl1_fxn * tg_yzzz_xxxxzzz_1[j];

                    tg_xyzzz_xxxxyyyy_0[j] = pb_x * tg_yzzz_xxxxyyyy_0[j] + fr * tg_yzzz_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_yzzz_xxxyyyy_1[j];

                    tg_xyzzz_xxxxyyyz_0[j] = pb_x * tg_yzzz_xxxxyyyz_0[j] + fr * tg_yzzz_xxxxyyyz_1[j] + 2.0 * fl1_fxn * tg_yzzz_xxxyyyz_1[j];

                    tg_xyzzz_xxxxyyzz_0[j] = pb_x * tg_yzzz_xxxxyyzz_0[j] + fr * tg_yzzz_xxxxyyzz_1[j] + 2.0 * fl1_fxn * tg_yzzz_xxxyyzz_1[j];

                    tg_xyzzz_xxxxyzzz_0[j] = pb_x * tg_yzzz_xxxxyzzz_0[j] + fr * tg_yzzz_xxxxyzzz_1[j] + 2.0 * fl1_fxn * tg_yzzz_xxxyzzz_1[j];

                    tg_xyzzz_xxxxzzzz_0[j] = pb_x * tg_yzzz_xxxxzzzz_0[j] + fr * tg_yzzz_xxxxzzzz_1[j] + 2.0 * fl1_fxn * tg_yzzz_xxxzzzz_1[j];

                    tg_xyzzz_xxxyyyyy_0[j] = pb_x * tg_yzzz_xxxyyyyy_0[j] + fr * tg_yzzz_xxxyyyyy_1[j] + 1.5 * fl1_fxn * tg_yzzz_xxyyyyy_1[j];

                    tg_xyzzz_xxxyyyyz_0[j] = pb_x * tg_yzzz_xxxyyyyz_0[j] + fr * tg_yzzz_xxxyyyyz_1[j] + 1.5 * fl1_fxn * tg_yzzz_xxyyyyz_1[j];

                    tg_xyzzz_xxxyyyzz_0[j] = pb_x * tg_yzzz_xxxyyyzz_0[j] + fr * tg_yzzz_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_yzzz_xxyyyzz_1[j];

                    tg_xyzzz_xxxyyzzz_0[j] = pb_x * tg_yzzz_xxxyyzzz_0[j] + fr * tg_yzzz_xxxyyzzz_1[j] + 1.5 * fl1_fxn * tg_yzzz_xxyyzzz_1[j];

                    tg_xyzzz_xxxyzzzz_0[j] = pb_x * tg_yzzz_xxxyzzzz_0[j] + fr * tg_yzzz_xxxyzzzz_1[j] + 1.5 * fl1_fxn * tg_yzzz_xxyzzzz_1[j];

                    tg_xyzzz_xxxzzzzz_0[j] = pb_x * tg_yzzz_xxxzzzzz_0[j] + fr * tg_yzzz_xxxzzzzz_1[j] + 1.5 * fl1_fxn * tg_yzzz_xxzzzzz_1[j];

                    tg_xyzzz_xxyyyyyy_0[j] = pb_x * tg_yzzz_xxyyyyyy_0[j] + fr * tg_yzzz_xxyyyyyy_1[j] + fl1_fxn * tg_yzzz_xyyyyyy_1[j];

                    tg_xyzzz_xxyyyyyz_0[j] = pb_x * tg_yzzz_xxyyyyyz_0[j] + fr * tg_yzzz_xxyyyyyz_1[j] + fl1_fxn * tg_yzzz_xyyyyyz_1[j];

                    tg_xyzzz_xxyyyyzz_0[j] = pb_x * tg_yzzz_xxyyyyzz_0[j] + fr * tg_yzzz_xxyyyyzz_1[j] + fl1_fxn * tg_yzzz_xyyyyzz_1[j];

                    tg_xyzzz_xxyyyzzz_0[j] = pb_x * tg_yzzz_xxyyyzzz_0[j] + fr * tg_yzzz_xxyyyzzz_1[j] + fl1_fxn * tg_yzzz_xyyyzzz_1[j];

                    tg_xyzzz_xxyyzzzz_0[j] = pb_x * tg_yzzz_xxyyzzzz_0[j] + fr * tg_yzzz_xxyyzzzz_1[j] + fl1_fxn * tg_yzzz_xyyzzzz_1[j];

                    tg_xyzzz_xxyzzzzz_0[j] = pb_x * tg_yzzz_xxyzzzzz_0[j] + fr * tg_yzzz_xxyzzzzz_1[j] + fl1_fxn * tg_yzzz_xyzzzzz_1[j];

                    tg_xyzzz_xxzzzzzz_0[j] = pb_x * tg_yzzz_xxzzzzzz_0[j] + fr * tg_yzzz_xxzzzzzz_1[j] + fl1_fxn * tg_yzzz_xzzzzzz_1[j];

                    tg_xyzzz_xyyyyyyy_0[j] = pb_x * tg_yzzz_xyyyyyyy_0[j] + fr * tg_yzzz_xyyyyyyy_1[j] + 0.5 * fl1_fxn * tg_yzzz_yyyyyyy_1[j];

                    tg_xyzzz_xyyyyyyz_0[j] = pb_x * tg_yzzz_xyyyyyyz_0[j] + fr * tg_yzzz_xyyyyyyz_1[j] + 0.5 * fl1_fxn * tg_yzzz_yyyyyyz_1[j];

                    tg_xyzzz_xyyyyyzz_0[j] = pb_x * tg_yzzz_xyyyyyzz_0[j] + fr * tg_yzzz_xyyyyyzz_1[j] + 0.5 * fl1_fxn * tg_yzzz_yyyyyzz_1[j];

                    tg_xyzzz_xyyyyzzz_0[j] = pb_x * tg_yzzz_xyyyyzzz_0[j] + fr * tg_yzzz_xyyyyzzz_1[j] + 0.5 * fl1_fxn * tg_yzzz_yyyyzzz_1[j];

                    tg_xyzzz_xyyyzzzz_0[j] = pb_x * tg_yzzz_xyyyzzzz_0[j] + fr * tg_yzzz_xyyyzzzz_1[j] + 0.5 * fl1_fxn * tg_yzzz_yyyzzzz_1[j];

                    tg_xyzzz_xyyzzzzz_0[j] = pb_x * tg_yzzz_xyyzzzzz_0[j] + fr * tg_yzzz_xyyzzzzz_1[j] + 0.5 * fl1_fxn * tg_yzzz_yyzzzzz_1[j];

                    tg_xyzzz_xyzzzzzz_0[j] = pb_x * tg_yzzz_xyzzzzzz_0[j] + fr * tg_yzzz_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yzzz_yzzzzzz_1[j];

                    tg_xyzzz_xzzzzzzz_0[j] = pb_x * tg_yzzz_xzzzzzzz_0[j] + fr * tg_yzzz_xzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_yzzz_zzzzzzz_1[j];

                    tg_xyzzz_yyyyyyyy_0[j] = pb_x * tg_yzzz_yyyyyyyy_0[j] + fr * tg_yzzz_yyyyyyyy_1[j];

                    tg_xyzzz_yyyyyyyz_0[j] = pb_x * tg_yzzz_yyyyyyyz_0[j] + fr * tg_yzzz_yyyyyyyz_1[j];

                    tg_xyzzz_yyyyyyzz_0[j] = pb_x * tg_yzzz_yyyyyyzz_0[j] + fr * tg_yzzz_yyyyyyzz_1[j];

                    tg_xyzzz_yyyyyzzz_0[j] = pb_x * tg_yzzz_yyyyyzzz_0[j] + fr * tg_yzzz_yyyyyzzz_1[j];

                    tg_xyzzz_yyyyzzzz_0[j] = pb_x * tg_yzzz_yyyyzzzz_0[j] + fr * tg_yzzz_yyyyzzzz_1[j];

                    tg_xyzzz_yyyzzzzz_0[j] = pb_x * tg_yzzz_yyyzzzzz_0[j] + fr * tg_yzzz_yyyzzzzz_1[j];

                    tg_xyzzz_yyzzzzzz_0[j] = pb_x * tg_yzzz_yyzzzzzz_0[j] + fr * tg_yzzz_yyzzzzzz_1[j];

                    tg_xyzzz_yzzzzzzz_0[j] = pb_x * tg_yzzz_yzzzzzzz_0[j] + fr * tg_yzzz_yzzzzzzz_1[j];

                    tg_xyzzz_zzzzzzzz_0[j] = pb_x * tg_yzzz_zzzzzzzz_0[j] + fr * tg_yzzz_zzzzzzzz_1[j];

                    tg_xzzzz_xxxxxxxx_0[j] = pb_x * tg_zzzz_xxxxxxxx_0[j] + fr * tg_zzzz_xxxxxxxx_1[j] + 4.0 * fl1_fxn * tg_zzzz_xxxxxxx_1[j];

                    tg_xzzzz_xxxxxxxy_0[j] = pb_x * tg_zzzz_xxxxxxxy_0[j] + fr * tg_zzzz_xxxxxxxy_1[j] + 3.5 * fl1_fxn * tg_zzzz_xxxxxxy_1[j];

                    tg_xzzzz_xxxxxxxz_0[j] = pb_x * tg_zzzz_xxxxxxxz_0[j] + fr * tg_zzzz_xxxxxxxz_1[j] + 3.5 * fl1_fxn * tg_zzzz_xxxxxxz_1[j];

                    tg_xzzzz_xxxxxxyy_0[j] = pb_x * tg_zzzz_xxxxxxyy_0[j] + fr * tg_zzzz_xxxxxxyy_1[j] + 3.0 * fl1_fxn * tg_zzzz_xxxxxyy_1[j];

                    tg_xzzzz_xxxxxxyz_0[j] = pb_x * tg_zzzz_xxxxxxyz_0[j] + fr * tg_zzzz_xxxxxxyz_1[j] + 3.0 * fl1_fxn * tg_zzzz_xxxxxyz_1[j];

                    tg_xzzzz_xxxxxxzz_0[j] = pb_x * tg_zzzz_xxxxxxzz_0[j] + fr * tg_zzzz_xxxxxxzz_1[j] + 3.0 * fl1_fxn * tg_zzzz_xxxxxzz_1[j];

                    tg_xzzzz_xxxxxyyy_0[j] = pb_x * tg_zzzz_xxxxxyyy_0[j] + fr * tg_zzzz_xxxxxyyy_1[j] + 2.5 * fl1_fxn * tg_zzzz_xxxxyyy_1[j];

                    tg_xzzzz_xxxxxyyz_0[j] = pb_x * tg_zzzz_xxxxxyyz_0[j] + fr * tg_zzzz_xxxxxyyz_1[j] + 2.5 * fl1_fxn * tg_zzzz_xxxxyyz_1[j];

                    tg_xzzzz_xxxxxyzz_0[j] = pb_x * tg_zzzz_xxxxxyzz_0[j] + fr * tg_zzzz_xxxxxyzz_1[j] + 2.5 * fl1_fxn * tg_zzzz_xxxxyzz_1[j];

                    tg_xzzzz_xxxxxzzz_0[j] = pb_x * tg_zzzz_xxxxxzzz_0[j] + fr * tg_zzzz_xxxxxzzz_1[j] + 2.5 * fl1_fxn * tg_zzzz_xxxxzzz_1[j];

                    tg_xzzzz_xxxxyyyy_0[j] = pb_x * tg_zzzz_xxxxyyyy_0[j] + fr * tg_zzzz_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_zzzz_xxxyyyy_1[j];

                    tg_xzzzz_xxxxyyyz_0[j] = pb_x * tg_zzzz_xxxxyyyz_0[j] + fr * tg_zzzz_xxxxyyyz_1[j] + 2.0 * fl1_fxn * tg_zzzz_xxxyyyz_1[j];

                    tg_xzzzz_xxxxyyzz_0[j] = pb_x * tg_zzzz_xxxxyyzz_0[j] + fr * tg_zzzz_xxxxyyzz_1[j] + 2.0 * fl1_fxn * tg_zzzz_xxxyyzz_1[j];

                    tg_xzzzz_xxxxyzzz_0[j] = pb_x * tg_zzzz_xxxxyzzz_0[j] + fr * tg_zzzz_xxxxyzzz_1[j] + 2.0 * fl1_fxn * tg_zzzz_xxxyzzz_1[j];

                    tg_xzzzz_xxxxzzzz_0[j] = pb_x * tg_zzzz_xxxxzzzz_0[j] + fr * tg_zzzz_xxxxzzzz_1[j] + 2.0 * fl1_fxn * tg_zzzz_xxxzzzz_1[j];

                    tg_xzzzz_xxxyyyyy_0[j] = pb_x * tg_zzzz_xxxyyyyy_0[j] + fr * tg_zzzz_xxxyyyyy_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxyyyyy_1[j];

                    tg_xzzzz_xxxyyyyz_0[j] = pb_x * tg_zzzz_xxxyyyyz_0[j] + fr * tg_zzzz_xxxyyyyz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxyyyyz_1[j];

                    tg_xzzzz_xxxyyyzz_0[j] = pb_x * tg_zzzz_xxxyyyzz_0[j] + fr * tg_zzzz_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxyyyzz_1[j];

                    tg_xzzzz_xxxyyzzz_0[j] = pb_x * tg_zzzz_xxxyyzzz_0[j] + fr * tg_zzzz_xxxyyzzz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxyyzzz_1[j];

                    tg_xzzzz_xxxyzzzz_0[j] = pb_x * tg_zzzz_xxxyzzzz_0[j] + fr * tg_zzzz_xxxyzzzz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxyzzzz_1[j];

                    tg_xzzzz_xxxzzzzz_0[j] = pb_x * tg_zzzz_xxxzzzzz_0[j] + fr * tg_zzzz_xxxzzzzz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxzzzzz_1[j];

                    tg_xzzzz_xxyyyyyy_0[j] = pb_x * tg_zzzz_xxyyyyyy_0[j] + fr * tg_zzzz_xxyyyyyy_1[j] + fl1_fxn * tg_zzzz_xyyyyyy_1[j];

                    tg_xzzzz_xxyyyyyz_0[j] = pb_x * tg_zzzz_xxyyyyyz_0[j] + fr * tg_zzzz_xxyyyyyz_1[j] + fl1_fxn * tg_zzzz_xyyyyyz_1[j];

                    tg_xzzzz_xxyyyyzz_0[j] = pb_x * tg_zzzz_xxyyyyzz_0[j] + fr * tg_zzzz_xxyyyyzz_1[j] + fl1_fxn * tg_zzzz_xyyyyzz_1[j];

                    tg_xzzzz_xxyyyzzz_0[j] = pb_x * tg_zzzz_xxyyyzzz_0[j] + fr * tg_zzzz_xxyyyzzz_1[j] + fl1_fxn * tg_zzzz_xyyyzzz_1[j];

                    tg_xzzzz_xxyyzzzz_0[j] = pb_x * tg_zzzz_xxyyzzzz_0[j] + fr * tg_zzzz_xxyyzzzz_1[j] + fl1_fxn * tg_zzzz_xyyzzzz_1[j];

                    tg_xzzzz_xxyzzzzz_0[j] = pb_x * tg_zzzz_xxyzzzzz_0[j] + fr * tg_zzzz_xxyzzzzz_1[j] + fl1_fxn * tg_zzzz_xyzzzzz_1[j];

                    tg_xzzzz_xxzzzzzz_0[j] = pb_x * tg_zzzz_xxzzzzzz_0[j] + fr * tg_zzzz_xxzzzzzz_1[j] + fl1_fxn * tg_zzzz_xzzzzzz_1[j];

                    tg_xzzzz_xyyyyyyy_0[j] = pb_x * tg_zzzz_xyyyyyyy_0[j] + fr * tg_zzzz_xyyyyyyy_1[j] + 0.5 * fl1_fxn * tg_zzzz_yyyyyyy_1[j];

                    tg_xzzzz_xyyyyyyz_0[j] = pb_x * tg_zzzz_xyyyyyyz_0[j] + fr * tg_zzzz_xyyyyyyz_1[j] + 0.5 * fl1_fxn * tg_zzzz_yyyyyyz_1[j];

                    tg_xzzzz_xyyyyyzz_0[j] = pb_x * tg_zzzz_xyyyyyzz_0[j] + fr * tg_zzzz_xyyyyyzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_yyyyyzz_1[j];

                    tg_xzzzz_xyyyyzzz_0[j] = pb_x * tg_zzzz_xyyyyzzz_0[j] + fr * tg_zzzz_xyyyyzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_yyyyzzz_1[j];

                    tg_xzzzz_xyyyzzzz_0[j] = pb_x * tg_zzzz_xyyyzzzz_0[j] + fr * tg_zzzz_xyyyzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_yyyzzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSL_663_757(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (663,757)

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
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_8_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_8_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_8_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_7_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {7, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_8_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_8_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {8, -1, -1, -1}, 
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

                auto tg_yyyy_xxxxxxxx_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 450); 

                auto tg_yyyy_xxxxxxxy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 451); 

                auto tg_yyyy_xxxxxxxz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 452); 

                auto tg_yyyy_xxxxxxyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 453); 

                auto tg_yyyy_xxxxxxyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 454); 

                auto tg_yyyy_xxxxxxzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 455); 

                auto tg_yyyy_xxxxxyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 456); 

                auto tg_yyyy_xxxxxyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 457); 

                auto tg_yyyy_xxxxxyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 458); 

                auto tg_yyyy_xxxxxzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 459); 

                auto tg_yyyy_xxxxyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 460); 

                auto tg_yyyy_xxxxyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 461); 

                auto tg_yyyy_xxxxyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 462); 

                auto tg_yyyy_xxxxyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 463); 

                auto tg_yyyy_xxxxzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 464); 

                auto tg_yyyy_xxxyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 465); 

                auto tg_yyyy_xxxyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 466); 

                auto tg_yyyy_xxxyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 467); 

                auto tg_yyyy_xxxyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 468); 

                auto tg_yyyy_xxxyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 469); 

                auto tg_yyyy_xxxzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 470); 

                auto tg_yyyy_xxyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 471); 

                auto tg_yyyy_xxyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 472); 

                auto tg_yyyy_xxyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 473); 

                auto tg_yyyy_xxyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 474); 

                auto tg_yyyy_xxyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 475); 

                auto tg_yyyy_xxyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 476); 

                auto tg_yyyy_xxzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 477); 

                auto tg_yyyy_xyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 478); 

                auto tg_yyyy_xyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 479); 

                auto tg_yyyy_xyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 480); 

                auto tg_yyyy_xyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 481); 

                auto tg_yyyy_xyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 482); 

                auto tg_yyyy_xyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 483); 

                auto tg_yyyy_xyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 484); 

                auto tg_yyyy_xzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 485); 

                auto tg_yyyy_yyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 486); 

                auto tg_yyyy_yyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 487); 

                auto tg_yyyy_yyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 488); 

                auto tg_yyyy_yyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 489); 

                auto tg_yyyy_yyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 490); 

                auto tg_yyyy_yyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 491); 

                auto tg_yyyy_yyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 492); 

                auto tg_yyyy_yzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 493); 

                auto tg_yyyy_zzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 494); 

                auto tg_yyyz_xxxxxxxx_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 495); 

                auto tg_yyyz_xxxxxxxy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 496); 

                auto tg_yyyz_xxxxxxxz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 497); 

                auto tg_yyyz_xxxxxxyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 498); 

                auto tg_yyyz_xxxxxxyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 499); 

                auto tg_yyyz_xxxxxxzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 500); 

                auto tg_yyyz_xxxxxyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 501); 

                auto tg_yyyz_xxxxxyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 502); 

                auto tg_yyyz_xxxxxyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 503); 

                auto tg_yyyz_xxxxxzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 504); 

                auto tg_yyyz_xxxxyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 505); 

                auto tg_yyyz_xxxxyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 506); 

                auto tg_yyyz_xxxxyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 507); 

                auto tg_yyyz_xxxxyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 508); 

                auto tg_yyyz_xxxxzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 509); 

                auto tg_yyyz_xxxyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 510); 

                auto tg_yyyz_xxxyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 511); 

                auto tg_yyyz_xxxyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 512); 

                auto tg_yyyz_xxxyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 513); 

                auto tg_yyyz_xxxyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 514); 

                auto tg_yyyz_xxxzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 515); 

                auto tg_yyyz_xxyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 516); 

                auto tg_yyyz_xxyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 517); 

                auto tg_yyyz_xxyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 518); 

                auto tg_yyyz_xxyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 519); 

                auto tg_yyyz_xxyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 520); 

                auto tg_yyyz_xxyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 521); 

                auto tg_yyyz_xxzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 522); 

                auto tg_yyyz_xyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 523); 

                auto tg_yyyz_xyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 524); 

                auto tg_yyyz_xyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 525); 

                auto tg_yyyz_xyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 526); 

                auto tg_yyyz_xyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 527); 

                auto tg_yyyz_xyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 528); 

                auto tg_yyyz_xyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 529); 

                auto tg_yyyz_xzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 530); 

                auto tg_yyyz_yyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 531); 

                auto tg_zzzz_xyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 663); 

                auto tg_zzzz_xyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 664); 

                auto tg_zzzz_xzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 665); 

                auto tg_zzzz_yyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 666); 

                auto tg_zzzz_yyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 667); 

                auto tg_zzzz_yyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 668); 

                auto tg_zzzz_yyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 669); 

                auto tg_zzzz_yyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 670); 

                auto tg_zzzz_yyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 671); 

                auto tg_zzzz_yyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 672); 

                auto tg_zzzz_yzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 673); 

                auto tg_zzzz_zzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 674); 

                auto tg_yyyy_xxxxxxxx_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 450); 

                auto tg_yyyy_xxxxxxxy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 451); 

                auto tg_yyyy_xxxxxxxz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 452); 

                auto tg_yyyy_xxxxxxyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 453); 

                auto tg_yyyy_xxxxxxyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 454); 

                auto tg_yyyy_xxxxxxzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 455); 

                auto tg_yyyy_xxxxxyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 456); 

                auto tg_yyyy_xxxxxyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 457); 

                auto tg_yyyy_xxxxxyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 458); 

                auto tg_yyyy_xxxxxzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 459); 

                auto tg_yyyy_xxxxyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 460); 

                auto tg_yyyy_xxxxyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 461); 

                auto tg_yyyy_xxxxyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 462); 

                auto tg_yyyy_xxxxyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 463); 

                auto tg_yyyy_xxxxzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 464); 

                auto tg_yyyy_xxxyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 465); 

                auto tg_yyyy_xxxyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 466); 

                auto tg_yyyy_xxxyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 467); 

                auto tg_yyyy_xxxyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 468); 

                auto tg_yyyy_xxxyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 469); 

                auto tg_yyyy_xxxzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 470); 

                auto tg_yyyy_xxyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 471); 

                auto tg_yyyy_xxyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 472); 

                auto tg_yyyy_xxyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 473); 

                auto tg_yyyy_xxyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 474); 

                auto tg_yyyy_xxyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 475); 

                auto tg_yyyy_xxyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 476); 

                auto tg_yyyy_xxzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 477); 

                auto tg_yyyy_xyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 478); 

                auto tg_yyyy_xyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 479); 

                auto tg_yyyy_xyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 480); 

                auto tg_yyyy_xyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 481); 

                auto tg_yyyy_xyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 482); 

                auto tg_yyyy_xyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 483); 

                auto tg_yyyy_xyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 484); 

                auto tg_yyyy_xzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 485); 

                auto tg_yyyy_yyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 486); 

                auto tg_yyyy_yyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 487); 

                auto tg_yyyy_yyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 488); 

                auto tg_yyyy_yyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 489); 

                auto tg_yyyy_yyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 490); 

                auto tg_yyyy_yyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 491); 

                auto tg_yyyy_yyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 492); 

                auto tg_yyyy_yzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 493); 

                auto tg_yyyy_zzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 494); 

                auto tg_yyyz_xxxxxxxx_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 495); 

                auto tg_yyyz_xxxxxxxy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 496); 

                auto tg_yyyz_xxxxxxxz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 497); 

                auto tg_yyyz_xxxxxxyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 498); 

                auto tg_yyyz_xxxxxxyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 499); 

                auto tg_yyyz_xxxxxxzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 500); 

                auto tg_yyyz_xxxxxyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 501); 

                auto tg_yyyz_xxxxxyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 502); 

                auto tg_yyyz_xxxxxyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 503); 

                auto tg_yyyz_xxxxxzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 504); 

                auto tg_yyyz_xxxxyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 505); 

                auto tg_yyyz_xxxxyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 506); 

                auto tg_yyyz_xxxxyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 507); 

                auto tg_yyyz_xxxxyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 508); 

                auto tg_yyyz_xxxxzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 509); 

                auto tg_yyyz_xxxyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 510); 

                auto tg_yyyz_xxxyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 511); 

                auto tg_yyyz_xxxyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 512); 

                auto tg_yyyz_xxxyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 513); 

                auto tg_yyyz_xxxyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 514); 

                auto tg_yyyz_xxxzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 515); 

                auto tg_yyyz_xxyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 516); 

                auto tg_yyyz_xxyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 517); 

                auto tg_yyyz_xxyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 518); 

                auto tg_yyyz_xxyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 519); 

                auto tg_yyyz_xxyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 520); 

                auto tg_yyyz_xxyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 521); 

                auto tg_yyyz_xxzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 522); 

                auto tg_yyyz_xyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 523); 

                auto tg_yyyz_xyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 524); 

                auto tg_yyyz_xyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 525); 

                auto tg_yyyz_xyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 526); 

                auto tg_yyyz_xyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 527); 

                auto tg_yyyz_xyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 528); 

                auto tg_yyyz_xyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 529); 

                auto tg_yyyz_xzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 530); 

                auto tg_yyyz_yyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 531); 

                auto tg_zzzz_xyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 663); 

                auto tg_zzzz_xyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 664); 

                auto tg_zzzz_xzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 665); 

                auto tg_zzzz_yyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 666); 

                auto tg_zzzz_yyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 667); 

                auto tg_zzzz_yyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 668); 

                auto tg_zzzz_yyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 669); 

                auto tg_zzzz_yyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 670); 

                auto tg_zzzz_yyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 671); 

                auto tg_zzzz_yyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 672); 

                auto tg_zzzz_yzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 673); 

                auto tg_zzzz_zzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 674); 

                auto tg_yyy_xxxxxxxx_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 270); 

                auto tg_yyy_xxxxxxxy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 271); 

                auto tg_yyy_xxxxxxxz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 272); 

                auto tg_yyy_xxxxxxyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 273); 

                auto tg_yyy_xxxxxxyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 274); 

                auto tg_yyy_xxxxxxzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 275); 

                auto tg_yyy_xxxxxyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 276); 

                auto tg_yyy_xxxxxyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 277); 

                auto tg_yyy_xxxxxyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 278); 

                auto tg_yyy_xxxxxzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 279); 

                auto tg_yyy_xxxxyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 280); 

                auto tg_yyy_xxxxyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 281); 

                auto tg_yyy_xxxxyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 282); 

                auto tg_yyy_xxxxyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 283); 

                auto tg_yyy_xxxxzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 284); 

                auto tg_yyy_xxxyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 285); 

                auto tg_yyy_xxxyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 286); 

                auto tg_yyy_xxxyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 287); 

                auto tg_yyy_xxxyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 288); 

                auto tg_yyy_xxxyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 289); 

                auto tg_yyy_xxxzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 290); 

                auto tg_yyy_xxyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 291); 

                auto tg_yyy_xxyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 292); 

                auto tg_yyy_xxyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 293); 

                auto tg_yyy_xxyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 294); 

                auto tg_yyy_xxyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 295); 

                auto tg_yyy_xxyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 296); 

                auto tg_yyy_xxzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 297); 

                auto tg_yyy_xyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 298); 

                auto tg_yyy_xyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 299); 

                auto tg_yyy_xyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 300); 

                auto tg_yyy_xyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 301); 

                auto tg_yyy_xyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 302); 

                auto tg_yyy_xyyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 303); 

                auto tg_yyy_xyzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 304); 

                auto tg_yyy_xzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 305); 

                auto tg_yyy_yyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 306); 

                auto tg_yyy_yyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 307); 

                auto tg_yyy_yyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 308); 

                auto tg_yyy_yyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 309); 

                auto tg_yyy_yyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 310); 

                auto tg_yyy_yyyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 311); 

                auto tg_yyy_yyzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 312); 

                auto tg_yyy_yzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 313); 

                auto tg_yyy_zzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 314); 

                auto tg_yyz_xxxxxxxx_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 315); 

                auto tg_yyz_xxxxxxxy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 316); 

                auto tg_yyz_xxxxxxxz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 317); 

                auto tg_yyz_xxxxxxyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 318); 

                auto tg_yyz_xxxxxxyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 319); 

                auto tg_yyz_xxxxxxzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 320); 

                auto tg_yyz_xxxxxyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 321); 

                auto tg_yyz_xxxxxyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 322); 

                auto tg_yyz_xxxxxyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 323); 

                auto tg_yyz_xxxxxzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 324); 

                auto tg_yyz_xxxxyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 325); 

                auto tg_yyz_xxxxyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 326); 

                auto tg_yyz_xxxxyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 327); 

                auto tg_yyz_xxxxyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 328); 

                auto tg_yyz_xxxxzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 329); 

                auto tg_yyz_xxxyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 330); 

                auto tg_yyz_xxxyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 331); 

                auto tg_yyz_xxxyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 332); 

                auto tg_yyz_xxxyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 333); 

                auto tg_yyz_xxxyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 334); 

                auto tg_yyz_xxxzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 335); 

                auto tg_yyz_xxyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 336); 

                auto tg_yyz_xxyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 337); 

                auto tg_yyz_xxyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 338); 

                auto tg_yyz_xxyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 339); 

                auto tg_yyz_xxyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 340); 

                auto tg_yyz_xxyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 341); 

                auto tg_yyz_xxzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 342); 

                auto tg_yyz_xyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 343); 

                auto tg_yyz_xyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 344); 

                auto tg_yyz_xyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 345); 

                auto tg_yyz_xyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 346); 

                auto tg_yyz_xyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 347); 

                auto tg_yyz_xyyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 348); 

                auto tg_yyz_xyzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 349); 

                auto tg_yyz_xzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 350); 

                auto tg_yyz_yyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 351); 

                auto tg_yyy_xxxxxxxx_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 270); 

                auto tg_yyy_xxxxxxxy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 271); 

                auto tg_yyy_xxxxxxxz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 272); 

                auto tg_yyy_xxxxxxyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 273); 

                auto tg_yyy_xxxxxxyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 274); 

                auto tg_yyy_xxxxxxzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 275); 

                auto tg_yyy_xxxxxyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 276); 

                auto tg_yyy_xxxxxyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 277); 

                auto tg_yyy_xxxxxyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 278); 

                auto tg_yyy_xxxxxzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 279); 

                auto tg_yyy_xxxxyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 280); 

                auto tg_yyy_xxxxyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 281); 

                auto tg_yyy_xxxxyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 282); 

                auto tg_yyy_xxxxyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 283); 

                auto tg_yyy_xxxxzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 284); 

                auto tg_yyy_xxxyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 285); 

                auto tg_yyy_xxxyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 286); 

                auto tg_yyy_xxxyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 287); 

                auto tg_yyy_xxxyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 288); 

                auto tg_yyy_xxxyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 289); 

                auto tg_yyy_xxxzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 290); 

                auto tg_yyy_xxyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 291); 

                auto tg_yyy_xxyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 292); 

                auto tg_yyy_xxyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 293); 

                auto tg_yyy_xxyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 294); 

                auto tg_yyy_xxyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 295); 

                auto tg_yyy_xxyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 296); 

                auto tg_yyy_xxzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 297); 

                auto tg_yyy_xyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 298); 

                auto tg_yyy_xyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 299); 

                auto tg_yyy_xyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 300); 

                auto tg_yyy_xyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 301); 

                auto tg_yyy_xyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 302); 

                auto tg_yyy_xyyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 303); 

                auto tg_yyy_xyzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 304); 

                auto tg_yyy_xzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 305); 

                auto tg_yyy_yyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 306); 

                auto tg_yyy_yyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 307); 

                auto tg_yyy_yyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 308); 

                auto tg_yyy_yyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 309); 

                auto tg_yyy_yyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 310); 

                auto tg_yyy_yyyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 311); 

                auto tg_yyy_yyzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 312); 

                auto tg_yyy_yzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 313); 

                auto tg_yyy_zzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 314); 

                auto tg_yyz_xxxxxxxx_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 315); 

                auto tg_yyz_xxxxxxxy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 316); 

                auto tg_yyz_xxxxxxxz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 317); 

                auto tg_yyz_xxxxxxyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 318); 

                auto tg_yyz_xxxxxxyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 319); 

                auto tg_yyz_xxxxxxzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 320); 

                auto tg_yyz_xxxxxyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 321); 

                auto tg_yyz_xxxxxyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 322); 

                auto tg_yyz_xxxxxyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 323); 

                auto tg_yyz_xxxxxzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 324); 

                auto tg_yyz_xxxxyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 325); 

                auto tg_yyz_xxxxyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 326); 

                auto tg_yyz_xxxxyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 327); 

                auto tg_yyz_xxxxyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 328); 

                auto tg_yyz_xxxxzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 329); 

                auto tg_yyz_xxxyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 330); 

                auto tg_yyz_xxxyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 331); 

                auto tg_yyz_xxxyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 332); 

                auto tg_yyz_xxxyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 333); 

                auto tg_yyz_xxxyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 334); 

                auto tg_yyz_xxxzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 335); 

                auto tg_yyz_xxyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 336); 

                auto tg_yyz_xxyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 337); 

                auto tg_yyz_xxyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 338); 

                auto tg_yyz_xxyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 339); 

                auto tg_yyz_xxyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 340); 

                auto tg_yyz_xxyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 341); 

                auto tg_yyz_xxzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 342); 

                auto tg_yyz_xyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 343); 

                auto tg_yyz_xyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 344); 

                auto tg_yyz_xyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 345); 

                auto tg_yyz_xyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 346); 

                auto tg_yyz_xyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 347); 

                auto tg_yyz_xyyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 348); 

                auto tg_yyz_xyzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 349); 

                auto tg_yyz_xzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 350); 

                auto tg_yyz_yyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 351); 

                auto tg_yyyy_xxxxxxx_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 360); 

                auto tg_yyyy_xxxxxxy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 361); 

                auto tg_yyyy_xxxxxxz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 362); 

                auto tg_yyyy_xxxxxyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 363); 

                auto tg_yyyy_xxxxxyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 364); 

                auto tg_yyyy_xxxxxzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 365); 

                auto tg_yyyy_xxxxyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 366); 

                auto tg_yyyy_xxxxyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 367); 

                auto tg_yyyy_xxxxyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 368); 

                auto tg_yyyy_xxxxzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 369); 

                auto tg_yyyy_xxxyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 370); 

                auto tg_yyyy_xxxyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 371); 

                auto tg_yyyy_xxxyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 372); 

                auto tg_yyyy_xxxyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 373); 

                auto tg_yyyy_xxxzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 374); 

                auto tg_yyyy_xxyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 375); 

                auto tg_yyyy_xxyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 376); 

                auto tg_yyyy_xxyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 377); 

                auto tg_yyyy_xxyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 378); 

                auto tg_yyyy_xxyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 379); 

                auto tg_yyyy_xxzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 380); 

                auto tg_yyyy_xyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 381); 

                auto tg_yyyy_xyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 382); 

                auto tg_yyyy_xyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 383); 

                auto tg_yyyy_xyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 384); 

                auto tg_yyyy_xyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 385); 

                auto tg_yyyy_xyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 386); 

                auto tg_yyyy_xzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 387); 

                auto tg_yyyy_yyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 388); 

                auto tg_yyyy_yyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 389); 

                auto tg_yyyy_yyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 390); 

                auto tg_yyyy_yyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 391); 

                auto tg_yyyy_yyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 392); 

                auto tg_yyyy_yyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 393); 

                auto tg_yyyy_yzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 394); 

                auto tg_yyyy_zzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 395); 

                auto tg_yyyz_xxxxxxx_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 396); 

                auto tg_yyyz_xxxxxxy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 397); 

                auto tg_yyyz_xxxxxxz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 398); 

                auto tg_yyyz_xxxxxyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 399); 

                auto tg_yyyz_xxxxxyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 400); 

                auto tg_yyyz_xxxxxzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 401); 

                auto tg_yyyz_xxxxyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 402); 

                auto tg_yyyz_xxxxyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 403); 

                auto tg_yyyz_xxxxyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 404); 

                auto tg_yyyz_xxxxzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 405); 

                auto tg_yyyz_xxxyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 406); 

                auto tg_yyyz_xxxyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 407); 

                auto tg_yyyz_xxxyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 408); 

                auto tg_yyyz_xxxyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 409); 

                auto tg_yyyz_xxxzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 410); 

                auto tg_yyyz_xxyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 411); 

                auto tg_yyyz_xxyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 412); 

                auto tg_yyyz_xxyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 413); 

                auto tg_yyyz_xxyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 414); 

                auto tg_yyyz_xxyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 415); 

                auto tg_yyyz_xxzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 416); 

                auto tg_yyyz_xyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 417); 

                auto tg_yyyz_xyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 418); 

                auto tg_yyyz_xyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 419); 

                auto tg_yyyz_xyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 420); 

                auto tg_yyyz_xyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 421); 

                auto tg_yyyz_xyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 422); 

                auto tg_yyyz_xzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 423); 

                auto tg_yyyz_yyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 424); 

                auto tg_zzzz_yyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 537); 

                auto tg_zzzz_yzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 538); 

                auto tg_zzzz_zzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 539); 

                // set up pointers to integrals

                auto tg_xzzzz_xyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 663); 

                auto tg_xzzzz_xyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 664); 

                auto tg_xzzzz_xzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 665); 

                auto tg_xzzzz_yyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 666); 

                auto tg_xzzzz_yyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 667); 

                auto tg_xzzzz_yyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 668); 

                auto tg_xzzzz_yyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 669); 

                auto tg_xzzzz_yyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 670); 

                auto tg_xzzzz_yyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 671); 

                auto tg_xzzzz_yyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 672); 

                auto tg_xzzzz_yzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 673); 

                auto tg_xzzzz_zzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 674); 

                auto tg_yyyyy_xxxxxxxx_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 675); 

                auto tg_yyyyy_xxxxxxxy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 676); 

                auto tg_yyyyy_xxxxxxxz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 677); 

                auto tg_yyyyy_xxxxxxyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 678); 

                auto tg_yyyyy_xxxxxxyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 679); 

                auto tg_yyyyy_xxxxxxzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 680); 

                auto tg_yyyyy_xxxxxyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 681); 

                auto tg_yyyyy_xxxxxyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 682); 

                auto tg_yyyyy_xxxxxyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 683); 

                auto tg_yyyyy_xxxxxzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 684); 

                auto tg_yyyyy_xxxxyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 685); 

                auto tg_yyyyy_xxxxyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 686); 

                auto tg_yyyyy_xxxxyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 687); 

                auto tg_yyyyy_xxxxyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 688); 

                auto tg_yyyyy_xxxxzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 689); 

                auto tg_yyyyy_xxxyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 690); 

                auto tg_yyyyy_xxxyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 691); 

                auto tg_yyyyy_xxxyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 692); 

                auto tg_yyyyy_xxxyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 693); 

                auto tg_yyyyy_xxxyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 694); 

                auto tg_yyyyy_xxxzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 695); 

                auto tg_yyyyy_xxyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 696); 

                auto tg_yyyyy_xxyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 697); 

                auto tg_yyyyy_xxyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 698); 

                auto tg_yyyyy_xxyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 699); 

                auto tg_yyyyy_xxyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 700); 

                auto tg_yyyyy_xxyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 701); 

                auto tg_yyyyy_xxzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 702); 

                auto tg_yyyyy_xyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 703); 

                auto tg_yyyyy_xyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 704); 

                auto tg_yyyyy_xyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 705); 

                auto tg_yyyyy_xyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 706); 

                auto tg_yyyyy_xyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 707); 

                auto tg_yyyyy_xyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 708); 

                auto tg_yyyyy_xyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 709); 

                auto tg_yyyyy_xzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 710); 

                auto tg_yyyyy_yyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 711); 

                auto tg_yyyyy_yyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 712); 

                auto tg_yyyyy_yyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 713); 

                auto tg_yyyyy_yyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 714); 

                auto tg_yyyyy_yyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 715); 

                auto tg_yyyyy_yyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 716); 

                auto tg_yyyyy_yyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 717); 

                auto tg_yyyyy_yzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 718); 

                auto tg_yyyyy_zzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 719); 

                auto tg_yyyyz_xxxxxxxx_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 720); 

                auto tg_yyyyz_xxxxxxxy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 721); 

                auto tg_yyyyz_xxxxxxxz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 722); 

                auto tg_yyyyz_xxxxxxyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 723); 

                auto tg_yyyyz_xxxxxxyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 724); 

                auto tg_yyyyz_xxxxxxzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 725); 

                auto tg_yyyyz_xxxxxyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 726); 

                auto tg_yyyyz_xxxxxyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 727); 

                auto tg_yyyyz_xxxxxyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 728); 

                auto tg_yyyyz_xxxxxzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 729); 

                auto tg_yyyyz_xxxxyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 730); 

                auto tg_yyyyz_xxxxyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 731); 

                auto tg_yyyyz_xxxxyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 732); 

                auto tg_yyyyz_xxxxyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 733); 

                auto tg_yyyyz_xxxxzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 734); 

                auto tg_yyyyz_xxxyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 735); 

                auto tg_yyyyz_xxxyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 736); 

                auto tg_yyyyz_xxxyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 737); 

                auto tg_yyyyz_xxxyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 738); 

                auto tg_yyyyz_xxxyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 739); 

                auto tg_yyyyz_xxxzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 740); 

                auto tg_yyyyz_xxyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 741); 

                auto tg_yyyyz_xxyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 742); 

                auto tg_yyyyz_xxyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 743); 

                auto tg_yyyyz_xxyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 744); 

                auto tg_yyyyz_xxyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 745); 

                auto tg_yyyyz_xxyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 746); 

                auto tg_yyyyz_xxzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 747); 

                auto tg_yyyyz_xyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 748); 

                auto tg_yyyyz_xyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 749); 

                auto tg_yyyyz_xyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 750); 

                auto tg_yyyyz_xyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 751); 

                auto tg_yyyyz_xyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 752); 

                auto tg_yyyyz_xyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 753); 

                auto tg_yyyyz_xyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 754); 

                auto tg_yyyyz_xzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 755); 

                auto tg_yyyyz_yyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 756); 

                // Batch of Integrals (663,757)

                #pragma omp simd aligned(fxn, fza, tg_xzzzz_xyyzzzzz_0, tg_xzzzz_xyzzzzzz_0, \
                                         tg_xzzzz_xzzzzzzz_0, tg_xzzzz_yyyyyyyy_0, tg_xzzzz_yyyyyyyz_0, tg_xzzzz_yyyyyyzz_0, \
                                         tg_xzzzz_yyyyyzzz_0, tg_xzzzz_yyyyzzzz_0, tg_xzzzz_yyyzzzzz_0, tg_xzzzz_yyzzzzzz_0, \
                                         tg_xzzzz_yzzzzzzz_0, tg_xzzzz_zzzzzzzz_0, tg_yyy_xxxxxxxx_0, tg_yyy_xxxxxxxx_1, \
                                         tg_yyy_xxxxxxxy_0, tg_yyy_xxxxxxxy_1, tg_yyy_xxxxxxxz_0, tg_yyy_xxxxxxxz_1, \
                                         tg_yyy_xxxxxxyy_0, tg_yyy_xxxxxxyy_1, tg_yyy_xxxxxxyz_0, tg_yyy_xxxxxxyz_1, \
                                         tg_yyy_xxxxxxzz_0, tg_yyy_xxxxxxzz_1, tg_yyy_xxxxxyyy_0, tg_yyy_xxxxxyyy_1, \
                                         tg_yyy_xxxxxyyz_0, tg_yyy_xxxxxyyz_1, tg_yyy_xxxxxyzz_0, tg_yyy_xxxxxyzz_1, \
                                         tg_yyy_xxxxxzzz_0, tg_yyy_xxxxxzzz_1, tg_yyy_xxxxyyyy_0, tg_yyy_xxxxyyyy_1, \
                                         tg_yyy_xxxxyyyz_0, tg_yyy_xxxxyyyz_1, tg_yyy_xxxxyyzz_0, tg_yyy_xxxxyyzz_1, \
                                         tg_yyy_xxxxyzzz_0, tg_yyy_xxxxyzzz_1, tg_yyy_xxxxzzzz_0, tg_yyy_xxxxzzzz_1, \
                                         tg_yyy_xxxyyyyy_0, tg_yyy_xxxyyyyy_1, tg_yyy_xxxyyyyz_0, tg_yyy_xxxyyyyz_1, \
                                         tg_yyy_xxxyyyzz_0, tg_yyy_xxxyyyzz_1, tg_yyy_xxxyyzzz_0, tg_yyy_xxxyyzzz_1, \
                                         tg_yyy_xxxyzzzz_0, tg_yyy_xxxyzzzz_1, tg_yyy_xxxzzzzz_0, tg_yyy_xxxzzzzz_1, \
                                         tg_yyy_xxyyyyyy_0, tg_yyy_xxyyyyyy_1, tg_yyy_xxyyyyyz_0, tg_yyy_xxyyyyyz_1, \
                                         tg_yyy_xxyyyyzz_0, tg_yyy_xxyyyyzz_1, tg_yyy_xxyyyzzz_0, tg_yyy_xxyyyzzz_1, \
                                         tg_yyy_xxyyzzzz_0, tg_yyy_xxyyzzzz_1, tg_yyy_xxyzzzzz_0, tg_yyy_xxyzzzzz_1, \
                                         tg_yyy_xxzzzzzz_0, tg_yyy_xxzzzzzz_1, tg_yyy_xyyyyyyy_0, tg_yyy_xyyyyyyy_1, \
                                         tg_yyy_xyyyyyyz_0, tg_yyy_xyyyyyyz_1, tg_yyy_xyyyyyzz_0, tg_yyy_xyyyyyzz_1, \
                                         tg_yyy_xyyyyzzz_0, tg_yyy_xyyyyzzz_1, tg_yyy_xyyyzzzz_0, tg_yyy_xyyyzzzz_1, \
                                         tg_yyy_xyyzzzzz_0, tg_yyy_xyyzzzzz_1, tg_yyy_xyzzzzzz_0, tg_yyy_xyzzzzzz_1, \
                                         tg_yyy_xzzzzzzz_0, tg_yyy_xzzzzzzz_1, tg_yyy_yyyyyyyy_0, tg_yyy_yyyyyyyy_1, \
                                         tg_yyy_yyyyyyyz_0, tg_yyy_yyyyyyyz_1, tg_yyy_yyyyyyzz_0, tg_yyy_yyyyyyzz_1, \
                                         tg_yyy_yyyyyzzz_0, tg_yyy_yyyyyzzz_1, tg_yyy_yyyyzzzz_0, tg_yyy_yyyyzzzz_1, \
                                         tg_yyy_yyyzzzzz_0, tg_yyy_yyyzzzzz_1, tg_yyy_yyzzzzzz_0, tg_yyy_yyzzzzzz_1, \
                                         tg_yyy_yzzzzzzz_0, tg_yyy_yzzzzzzz_1, tg_yyy_zzzzzzzz_0, tg_yyy_zzzzzzzz_1, \
                                         tg_yyyy_xxxxxxx_1, tg_yyyy_xxxxxxxx_0, tg_yyyy_xxxxxxxx_1, tg_yyyy_xxxxxxxy_0, \
                                         tg_yyyy_xxxxxxxy_1, tg_yyyy_xxxxxxxz_0, tg_yyyy_xxxxxxxz_1, tg_yyyy_xxxxxxy_1, \
                                         tg_yyyy_xxxxxxyy_0, tg_yyyy_xxxxxxyy_1, tg_yyyy_xxxxxxyz_0, tg_yyyy_xxxxxxyz_1, \
                                         tg_yyyy_xxxxxxz_1, tg_yyyy_xxxxxxzz_0, tg_yyyy_xxxxxxzz_1, tg_yyyy_xxxxxyy_1, \
                                         tg_yyyy_xxxxxyyy_0, tg_yyyy_xxxxxyyy_1, tg_yyyy_xxxxxyyz_0, tg_yyyy_xxxxxyyz_1, \
                                         tg_yyyy_xxxxxyz_1, tg_yyyy_xxxxxyzz_0, tg_yyyy_xxxxxyzz_1, tg_yyyy_xxxxxzz_1, \
                                         tg_yyyy_xxxxxzzz_0, tg_yyyy_xxxxxzzz_1, tg_yyyy_xxxxyyy_1, tg_yyyy_xxxxyyyy_0, \
                                         tg_yyyy_xxxxyyyy_1, tg_yyyy_xxxxyyyz_0, tg_yyyy_xxxxyyyz_1, tg_yyyy_xxxxyyz_1, \
                                         tg_yyyy_xxxxyyzz_0, tg_yyyy_xxxxyyzz_1, tg_yyyy_xxxxyzz_1, tg_yyyy_xxxxyzzz_0, \
                                         tg_yyyy_xxxxyzzz_1, tg_yyyy_xxxxzzz_1, tg_yyyy_xxxxzzzz_0, tg_yyyy_xxxxzzzz_1, \
                                         tg_yyyy_xxxyyyy_1, tg_yyyy_xxxyyyyy_0, tg_yyyy_xxxyyyyy_1, tg_yyyy_xxxyyyyz_0, \
                                         tg_yyyy_xxxyyyyz_1, tg_yyyy_xxxyyyz_1, tg_yyyy_xxxyyyzz_0, tg_yyyy_xxxyyyzz_1, \
                                         tg_yyyy_xxxyyzz_1, tg_yyyy_xxxyyzzz_0, tg_yyyy_xxxyyzzz_1, tg_yyyy_xxxyzzz_1, \
                                         tg_yyyy_xxxyzzzz_0, tg_yyyy_xxxyzzzz_1, tg_yyyy_xxxzzzz_1, tg_yyyy_xxxzzzzz_0, \
                                         tg_yyyy_xxxzzzzz_1, tg_yyyy_xxyyyyy_1, tg_yyyy_xxyyyyyy_0, tg_yyyy_xxyyyyyy_1, \
                                         tg_yyyy_xxyyyyyz_0, tg_yyyy_xxyyyyyz_1, tg_yyyy_xxyyyyz_1, tg_yyyy_xxyyyyzz_0, \
                                         tg_yyyy_xxyyyyzz_1, tg_yyyy_xxyyyzz_1, tg_yyyy_xxyyyzzz_0, tg_yyyy_xxyyyzzz_1, \
                                         tg_yyyy_xxyyzzz_1, tg_yyyy_xxyyzzzz_0, tg_yyyy_xxyyzzzz_1, tg_yyyy_xxyzzzz_1, \
                                         tg_yyyy_xxyzzzzz_0, tg_yyyy_xxyzzzzz_1, tg_yyyy_xxzzzzz_1, tg_yyyy_xxzzzzzz_0, \
                                         tg_yyyy_xxzzzzzz_1, tg_yyyy_xyyyyyy_1, tg_yyyy_xyyyyyyy_0, tg_yyyy_xyyyyyyy_1, \
                                         tg_yyyy_xyyyyyyz_0, tg_yyyy_xyyyyyyz_1, tg_yyyy_xyyyyyz_1, tg_yyyy_xyyyyyzz_0, \
                                         tg_yyyy_xyyyyyzz_1, tg_yyyy_xyyyyzz_1, tg_yyyy_xyyyyzzz_0, tg_yyyy_xyyyyzzz_1, \
                                         tg_yyyy_xyyyzzz_1, tg_yyyy_xyyyzzzz_0, tg_yyyy_xyyyzzzz_1, tg_yyyy_xyyzzzz_1, \
                                         tg_yyyy_xyyzzzzz_0, tg_yyyy_xyyzzzzz_1, tg_yyyy_xyzzzzz_1, tg_yyyy_xyzzzzzz_0, \
                                         tg_yyyy_xyzzzzzz_1, tg_yyyy_xzzzzzz_1, tg_yyyy_xzzzzzzz_0, tg_yyyy_xzzzzzzz_1, \
                                         tg_yyyy_yyyyyyy_1, tg_yyyy_yyyyyyyy_0, tg_yyyy_yyyyyyyy_1, tg_yyyy_yyyyyyyz_0, \
                                         tg_yyyy_yyyyyyyz_1, tg_yyyy_yyyyyyz_1, tg_yyyy_yyyyyyzz_0, tg_yyyy_yyyyyyzz_1, \
                                         tg_yyyy_yyyyyzz_1, tg_yyyy_yyyyyzzz_0, tg_yyyy_yyyyyzzz_1, tg_yyyy_yyyyzzz_1, \
                                         tg_yyyy_yyyyzzzz_0, tg_yyyy_yyyyzzzz_1, tg_yyyy_yyyzzzz_1, tg_yyyy_yyyzzzzz_0, \
                                         tg_yyyy_yyyzzzzz_1, tg_yyyy_yyzzzzz_1, tg_yyyy_yyzzzzzz_0, tg_yyyy_yyzzzzzz_1, \
                                         tg_yyyy_yzzzzzz_1, tg_yyyy_yzzzzzzz_0, tg_yyyy_yzzzzzzz_1, tg_yyyy_zzzzzzz_1, \
                                         tg_yyyy_zzzzzzzz_0, tg_yyyy_zzzzzzzz_1, tg_yyyyy_xxxxxxxx_0, tg_yyyyy_xxxxxxxy_0, \
                                         tg_yyyyy_xxxxxxxz_0, tg_yyyyy_xxxxxxyy_0, tg_yyyyy_xxxxxxyz_0, tg_yyyyy_xxxxxxzz_0, \
                                         tg_yyyyy_xxxxxyyy_0, tg_yyyyy_xxxxxyyz_0, tg_yyyyy_xxxxxyzz_0, tg_yyyyy_xxxxxzzz_0, \
                                         tg_yyyyy_xxxxyyyy_0, tg_yyyyy_xxxxyyyz_0, tg_yyyyy_xxxxyyzz_0, tg_yyyyy_xxxxyzzz_0, \
                                         tg_yyyyy_xxxxzzzz_0, tg_yyyyy_xxxyyyyy_0, tg_yyyyy_xxxyyyyz_0, tg_yyyyy_xxxyyyzz_0, \
                                         tg_yyyyy_xxxyyzzz_0, tg_yyyyy_xxxyzzzz_0, tg_yyyyy_xxxzzzzz_0, tg_yyyyy_xxyyyyyy_0, \
                                         tg_yyyyy_xxyyyyyz_0, tg_yyyyy_xxyyyyzz_0, tg_yyyyy_xxyyyzzz_0, tg_yyyyy_xxyyzzzz_0, \
                                         tg_yyyyy_xxyzzzzz_0, tg_yyyyy_xxzzzzzz_0, tg_yyyyy_xyyyyyyy_0, tg_yyyyy_xyyyyyyz_0, \
                                         tg_yyyyy_xyyyyyzz_0, tg_yyyyy_xyyyyzzz_0, tg_yyyyy_xyyyzzzz_0, tg_yyyyy_xyyzzzzz_0, \
                                         tg_yyyyy_xyzzzzzz_0, tg_yyyyy_xzzzzzzz_0, tg_yyyyy_yyyyyyyy_0, tg_yyyyy_yyyyyyyz_0, \
                                         tg_yyyyy_yyyyyyzz_0, tg_yyyyy_yyyyyzzz_0, tg_yyyyy_yyyyzzzz_0, tg_yyyyy_yyyzzzzz_0, \
                                         tg_yyyyy_yyzzzzzz_0, tg_yyyyy_yzzzzzzz_0, tg_yyyyy_zzzzzzzz_0, tg_yyyyz_xxxxxxxx_0, \
                                         tg_yyyyz_xxxxxxxy_0, tg_yyyyz_xxxxxxxz_0, tg_yyyyz_xxxxxxyy_0, tg_yyyyz_xxxxxxyz_0, \
                                         tg_yyyyz_xxxxxxzz_0, tg_yyyyz_xxxxxyyy_0, tg_yyyyz_xxxxxyyz_0, tg_yyyyz_xxxxxyzz_0, \
                                         tg_yyyyz_xxxxxzzz_0, tg_yyyyz_xxxxyyyy_0, tg_yyyyz_xxxxyyyz_0, tg_yyyyz_xxxxyyzz_0, \
                                         tg_yyyyz_xxxxyzzz_0, tg_yyyyz_xxxxzzzz_0, tg_yyyyz_xxxyyyyy_0, tg_yyyyz_xxxyyyyz_0, \
                                         tg_yyyyz_xxxyyyzz_0, tg_yyyyz_xxxyyzzz_0, tg_yyyyz_xxxyzzzz_0, tg_yyyyz_xxxzzzzz_0, \
                                         tg_yyyyz_xxyyyyyy_0, tg_yyyyz_xxyyyyyz_0, tg_yyyyz_xxyyyyzz_0, tg_yyyyz_xxyyyzzz_0, \
                                         tg_yyyyz_xxyyzzzz_0, tg_yyyyz_xxyzzzzz_0, tg_yyyyz_xxzzzzzz_0, tg_yyyyz_xyyyyyyy_0, \
                                         tg_yyyyz_xyyyyyyz_0, tg_yyyyz_xyyyyyzz_0, tg_yyyyz_xyyyyzzz_0, tg_yyyyz_xyyyzzzz_0, \
                                         tg_yyyyz_xyyzzzzz_0, tg_yyyyz_xyzzzzzz_0, tg_yyyyz_xzzzzzzz_0, tg_yyyyz_yyyyyyyy_0, \
                                         tg_yyyz_xxxxxxx_1, tg_yyyz_xxxxxxxx_0, tg_yyyz_xxxxxxxx_1, tg_yyyz_xxxxxxxy_0, \
                                         tg_yyyz_xxxxxxxy_1, tg_yyyz_xxxxxxxz_0, tg_yyyz_xxxxxxxz_1, tg_yyyz_xxxxxxy_1, \
                                         tg_yyyz_xxxxxxyy_0, tg_yyyz_xxxxxxyy_1, tg_yyyz_xxxxxxyz_0, tg_yyyz_xxxxxxyz_1, \
                                         tg_yyyz_xxxxxxz_1, tg_yyyz_xxxxxxzz_0, tg_yyyz_xxxxxxzz_1, tg_yyyz_xxxxxyy_1, \
                                         tg_yyyz_xxxxxyyy_0, tg_yyyz_xxxxxyyy_1, tg_yyyz_xxxxxyyz_0, tg_yyyz_xxxxxyyz_1, \
                                         tg_yyyz_xxxxxyz_1, tg_yyyz_xxxxxyzz_0, tg_yyyz_xxxxxyzz_1, tg_yyyz_xxxxxzz_1, \
                                         tg_yyyz_xxxxxzzz_0, tg_yyyz_xxxxxzzz_1, tg_yyyz_xxxxyyy_1, tg_yyyz_xxxxyyyy_0, \
                                         tg_yyyz_xxxxyyyy_1, tg_yyyz_xxxxyyyz_0, tg_yyyz_xxxxyyyz_1, tg_yyyz_xxxxyyz_1, \
                                         tg_yyyz_xxxxyyzz_0, tg_yyyz_xxxxyyzz_1, tg_yyyz_xxxxyzz_1, tg_yyyz_xxxxyzzz_0, \
                                         tg_yyyz_xxxxyzzz_1, tg_yyyz_xxxxzzz_1, tg_yyyz_xxxxzzzz_0, tg_yyyz_xxxxzzzz_1, \
                                         tg_yyyz_xxxyyyy_1, tg_yyyz_xxxyyyyy_0, tg_yyyz_xxxyyyyy_1, tg_yyyz_xxxyyyyz_0, \
                                         tg_yyyz_xxxyyyyz_1, tg_yyyz_xxxyyyz_1, tg_yyyz_xxxyyyzz_0, tg_yyyz_xxxyyyzz_1, \
                                         tg_yyyz_xxxyyzz_1, tg_yyyz_xxxyyzzz_0, tg_yyyz_xxxyyzzz_1, tg_yyyz_xxxyzzz_1, \
                                         tg_yyyz_xxxyzzzz_0, tg_yyyz_xxxyzzzz_1, tg_yyyz_xxxzzzz_1, tg_yyyz_xxxzzzzz_0, \
                                         tg_yyyz_xxxzzzzz_1, tg_yyyz_xxyyyyy_1, tg_yyyz_xxyyyyyy_0, tg_yyyz_xxyyyyyy_1, \
                                         tg_yyyz_xxyyyyyz_0, tg_yyyz_xxyyyyyz_1, tg_yyyz_xxyyyyz_1, tg_yyyz_xxyyyyzz_0, \
                                         tg_yyyz_xxyyyyzz_1, tg_yyyz_xxyyyzz_1, tg_yyyz_xxyyyzzz_0, tg_yyyz_xxyyyzzz_1, \
                                         tg_yyyz_xxyyzzz_1, tg_yyyz_xxyyzzzz_0, tg_yyyz_xxyyzzzz_1, tg_yyyz_xxyzzzz_1, \
                                         tg_yyyz_xxyzzzzz_0, tg_yyyz_xxyzzzzz_1, tg_yyyz_xxzzzzz_1, tg_yyyz_xxzzzzzz_0, \
                                         tg_yyyz_xxzzzzzz_1, tg_yyyz_xyyyyyy_1, tg_yyyz_xyyyyyyy_0, tg_yyyz_xyyyyyyy_1, \
                                         tg_yyyz_xyyyyyyz_0, tg_yyyz_xyyyyyyz_1, tg_yyyz_xyyyyyz_1, tg_yyyz_xyyyyyzz_0, \
                                         tg_yyyz_xyyyyyzz_1, tg_yyyz_xyyyyzz_1, tg_yyyz_xyyyyzzz_0, tg_yyyz_xyyyyzzz_1, \
                                         tg_yyyz_xyyyzzz_1, tg_yyyz_xyyyzzzz_0, tg_yyyz_xyyyzzzz_1, tg_yyyz_xyyzzzz_1, \
                                         tg_yyyz_xyyzzzzz_0, tg_yyyz_xyyzzzzz_1, tg_yyyz_xyzzzzz_1, tg_yyyz_xyzzzzzz_0, \
                                         tg_yyyz_xyzzzzzz_1, tg_yyyz_xzzzzzz_1, tg_yyyz_xzzzzzzz_0, tg_yyyz_xzzzzzzz_1, \
                                         tg_yyyz_yyyyyyy_1, tg_yyyz_yyyyyyyy_0, tg_yyyz_yyyyyyyy_1, tg_yyz_xxxxxxxx_0, \
                                         tg_yyz_xxxxxxxx_1, tg_yyz_xxxxxxxy_0, tg_yyz_xxxxxxxy_1, tg_yyz_xxxxxxxz_0, \
                                         tg_yyz_xxxxxxxz_1, tg_yyz_xxxxxxyy_0, tg_yyz_xxxxxxyy_1, tg_yyz_xxxxxxyz_0, \
                                         tg_yyz_xxxxxxyz_1, tg_yyz_xxxxxxzz_0, tg_yyz_xxxxxxzz_1, tg_yyz_xxxxxyyy_0, \
                                         tg_yyz_xxxxxyyy_1, tg_yyz_xxxxxyyz_0, tg_yyz_xxxxxyyz_1, tg_yyz_xxxxxyzz_0, \
                                         tg_yyz_xxxxxyzz_1, tg_yyz_xxxxxzzz_0, tg_yyz_xxxxxzzz_1, tg_yyz_xxxxyyyy_0, \
                                         tg_yyz_xxxxyyyy_1, tg_yyz_xxxxyyyz_0, tg_yyz_xxxxyyyz_1, tg_yyz_xxxxyyzz_0, \
                                         tg_yyz_xxxxyyzz_1, tg_yyz_xxxxyzzz_0, tg_yyz_xxxxyzzz_1, tg_yyz_xxxxzzzz_0, \
                                         tg_yyz_xxxxzzzz_1, tg_yyz_xxxyyyyy_0, tg_yyz_xxxyyyyy_1, tg_yyz_xxxyyyyz_0, \
                                         tg_yyz_xxxyyyyz_1, tg_yyz_xxxyyyzz_0, tg_yyz_xxxyyyzz_1, tg_yyz_xxxyyzzz_0, \
                                         tg_yyz_xxxyyzzz_1, tg_yyz_xxxyzzzz_0, tg_yyz_xxxyzzzz_1, tg_yyz_xxxzzzzz_0, \
                                         tg_yyz_xxxzzzzz_1, tg_yyz_xxyyyyyy_0, tg_yyz_xxyyyyyy_1, tg_yyz_xxyyyyyz_0, \
                                         tg_yyz_xxyyyyyz_1, tg_yyz_xxyyyyzz_0, tg_yyz_xxyyyyzz_1, tg_yyz_xxyyyzzz_0, \
                                         tg_yyz_xxyyyzzz_1, tg_yyz_xxyyzzzz_0, tg_yyz_xxyyzzzz_1, tg_yyz_xxyzzzzz_0, \
                                         tg_yyz_xxyzzzzz_1, tg_yyz_xxzzzzzz_0, tg_yyz_xxzzzzzz_1, tg_yyz_xyyyyyyy_0, \
                                         tg_yyz_xyyyyyyy_1, tg_yyz_xyyyyyyz_0, tg_yyz_xyyyyyyz_1, tg_yyz_xyyyyyzz_0, \
                                         tg_yyz_xyyyyyzz_1, tg_yyz_xyyyyzzz_0, tg_yyz_xyyyyzzz_1, tg_yyz_xyyyzzzz_0, \
                                         tg_yyz_xyyyzzzz_1, tg_yyz_xyyzzzzz_0, tg_yyz_xyyzzzzz_1, tg_yyz_xyzzzzzz_0, \
                                         tg_yyz_xyzzzzzz_1, tg_yyz_xzzzzzzz_0, tg_yyz_xzzzzzzz_1, tg_yyz_yyyyyyyy_0, \
                                         tg_yyz_yyyyyyyy_1, tg_zzzz_xyyzzzzz_0, tg_zzzz_xyyzzzzz_1, tg_zzzz_xyzzzzzz_0, \
                                         tg_zzzz_xyzzzzzz_1, tg_zzzz_xzzzzzzz_0, tg_zzzz_xzzzzzzz_1, tg_zzzz_yyyyyyyy_0, \
                                         tg_zzzz_yyyyyyyy_1, tg_zzzz_yyyyyyyz_0, tg_zzzz_yyyyyyyz_1, tg_zzzz_yyyyyyzz_0, \
                                         tg_zzzz_yyyyyyzz_1, tg_zzzz_yyyyyzzz_0, tg_zzzz_yyyyyzzz_1, tg_zzzz_yyyyzzzz_0, \
                                         tg_zzzz_yyyyzzzz_1, tg_zzzz_yyyzzzzz_0, tg_zzzz_yyyzzzzz_1, tg_zzzz_yyzzzzz_1, \
                                         tg_zzzz_yyzzzzzz_0, tg_zzzz_yyzzzzzz_1, tg_zzzz_yzzzzzz_1, tg_zzzz_yzzzzzzz_0, \
                                         tg_zzzz_yzzzzzzz_1, tg_zzzz_zzzzzzz_1, tg_zzzz_zzzzzzzz_0, tg_zzzz_zzzzzzzz_1, wp_x, \
                                         wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xzzzz_xyyzzzzz_0[j] = pb_x * tg_zzzz_xyyzzzzz_0[j] + fr * tg_zzzz_xyyzzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_yyzzzzz_1[j];

                    tg_xzzzz_xyzzzzzz_0[j] = pb_x * tg_zzzz_xyzzzzzz_0[j] + fr * tg_zzzz_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_yzzzzzz_1[j];

                    tg_xzzzz_xzzzzzzz_0[j] = pb_x * tg_zzzz_xzzzzzzz_0[j] + fr * tg_zzzz_xzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_zzzzzzz_1[j];

                    tg_xzzzz_yyyyyyyy_0[j] = pb_x * tg_zzzz_yyyyyyyy_0[j] + fr * tg_zzzz_yyyyyyyy_1[j];

                    tg_xzzzz_yyyyyyyz_0[j] = pb_x * tg_zzzz_yyyyyyyz_0[j] + fr * tg_zzzz_yyyyyyyz_1[j];

                    tg_xzzzz_yyyyyyzz_0[j] = pb_x * tg_zzzz_yyyyyyzz_0[j] + fr * tg_zzzz_yyyyyyzz_1[j];

                    tg_xzzzz_yyyyyzzz_0[j] = pb_x * tg_zzzz_yyyyyzzz_0[j] + fr * tg_zzzz_yyyyyzzz_1[j];

                    tg_xzzzz_yyyyzzzz_0[j] = pb_x * tg_zzzz_yyyyzzzz_0[j] + fr * tg_zzzz_yyyyzzzz_1[j];

                    tg_xzzzz_yyyzzzzz_0[j] = pb_x * tg_zzzz_yyyzzzzz_0[j] + fr * tg_zzzz_yyyzzzzz_1[j];

                    tg_xzzzz_yyzzzzzz_0[j] = pb_x * tg_zzzz_yyzzzzzz_0[j] + fr * tg_zzzz_yyzzzzzz_1[j];

                    tg_xzzzz_yzzzzzzz_0[j] = pb_x * tg_zzzz_yzzzzzzz_0[j] + fr * tg_zzzz_yzzzzzzz_1[j];

                    tg_xzzzz_zzzzzzzz_0[j] = pb_x * tg_zzzz_zzzzzzzz_0[j] + fr * tg_zzzz_zzzzzzzz_1[j];

                    fr = wp_y[j]; 

                    tg_yyyyy_xxxxxxxx_0[j] = pb_y * tg_yyyy_xxxxxxxx_0[j] + fr * tg_yyyy_xxxxxxxx_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxxxxxx_0[j] - tg_yyy_xxxxxxxx_1[j] * fl1_fza);

                    tg_yyyyy_xxxxxxxy_0[j] = pb_y * tg_yyyy_xxxxxxxy_0[j] + fr * tg_yyyy_xxxxxxxy_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxxxxxy_0[j] - tg_yyy_xxxxxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyy_xxxxxxx_1[j];

                    tg_yyyyy_xxxxxxxz_0[j] = pb_y * tg_yyyy_xxxxxxxz_0[j] + fr * tg_yyyy_xxxxxxxz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxxxxxz_0[j] - tg_yyy_xxxxxxxz_1[j] * fl1_fza);

                    tg_yyyyy_xxxxxxyy_0[j] = pb_y * tg_yyyy_xxxxxxyy_0[j] + fr * tg_yyyy_xxxxxxyy_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxxxxyy_0[j] - tg_yyy_xxxxxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyy_xxxxxxy_1[j];

                    tg_yyyyy_xxxxxxyz_0[j] = pb_y * tg_yyyy_xxxxxxyz_0[j] + fr * tg_yyyy_xxxxxxyz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxxxxyz_0[j] - tg_yyy_xxxxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyy_xxxxxxz_1[j];

                    tg_yyyyy_xxxxxxzz_0[j] = pb_y * tg_yyyy_xxxxxxzz_0[j] + fr * tg_yyyy_xxxxxxzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxxxxzz_0[j] - tg_yyy_xxxxxxzz_1[j] * fl1_fza);

                    tg_yyyyy_xxxxxyyy_0[j] = pb_y * tg_yyyy_xxxxxyyy_0[j] + fr * tg_yyyy_xxxxxyyy_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxxxyyy_0[j] - tg_yyy_xxxxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyy_xxxxxyy_1[j];

                    tg_yyyyy_xxxxxyyz_0[j] = pb_y * tg_yyyy_xxxxxyyz_0[j] + fr * tg_yyyy_xxxxxyyz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxxxyyz_0[j] - tg_yyy_xxxxxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyy_xxxxxyz_1[j];

                    tg_yyyyy_xxxxxyzz_0[j] = pb_y * tg_yyyy_xxxxxyzz_0[j] + fr * tg_yyyy_xxxxxyzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxxxyzz_0[j] - tg_yyy_xxxxxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyy_xxxxxzz_1[j];

                    tg_yyyyy_xxxxxzzz_0[j] = pb_y * tg_yyyy_xxxxxzzz_0[j] + fr * tg_yyyy_xxxxxzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxxxzzz_0[j] - tg_yyy_xxxxxzzz_1[j] * fl1_fza);

                    tg_yyyyy_xxxxyyyy_0[j] = pb_y * tg_yyyy_xxxxyyyy_0[j] + fr * tg_yyyy_xxxxyyyy_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxxyyyy_0[j] - tg_yyy_xxxxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyy_xxxxyyy_1[j];

                    tg_yyyyy_xxxxyyyz_0[j] = pb_y * tg_yyyy_xxxxyyyz_0[j] + fr * tg_yyyy_xxxxyyyz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxxyyyz_0[j] - tg_yyy_xxxxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyy_xxxxyyz_1[j];

                    tg_yyyyy_xxxxyyzz_0[j] = pb_y * tg_yyyy_xxxxyyzz_0[j] + fr * tg_yyyy_xxxxyyzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxxyyzz_0[j] - tg_yyy_xxxxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyy_xxxxyzz_1[j];

                    tg_yyyyy_xxxxyzzz_0[j] = pb_y * tg_yyyy_xxxxyzzz_0[j] + fr * tg_yyyy_xxxxyzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxxyzzz_0[j] - tg_yyy_xxxxyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyy_xxxxzzz_1[j];

                    tg_yyyyy_xxxxzzzz_0[j] = pb_y * tg_yyyy_xxxxzzzz_0[j] + fr * tg_yyyy_xxxxzzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxxzzzz_0[j] - tg_yyy_xxxxzzzz_1[j] * fl1_fza);

                    tg_yyyyy_xxxyyyyy_0[j] = pb_y * tg_yyyy_xxxyyyyy_0[j] + fr * tg_yyyy_xxxyyyyy_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxyyyyy_0[j] - tg_yyy_xxxyyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyy_xxxyyyy_1[j];

                    tg_yyyyy_xxxyyyyz_0[j] = pb_y * tg_yyyy_xxxyyyyz_0[j] + fr * tg_yyyy_xxxyyyyz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxyyyyz_0[j] - tg_yyy_xxxyyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyy_xxxyyyz_1[j];

                    tg_yyyyy_xxxyyyzz_0[j] = pb_y * tg_yyyy_xxxyyyzz_0[j] + fr * tg_yyyy_xxxyyyzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxyyyzz_0[j] - tg_yyy_xxxyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyy_xxxyyzz_1[j];

                    tg_yyyyy_xxxyyzzz_0[j] = pb_y * tg_yyyy_xxxyyzzz_0[j] + fr * tg_yyyy_xxxyyzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxyyzzz_0[j] - tg_yyy_xxxyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyy_xxxyzzz_1[j];

                    tg_yyyyy_xxxyzzzz_0[j] = pb_y * tg_yyyy_xxxyzzzz_0[j] + fr * tg_yyyy_xxxyzzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxyzzzz_0[j] - tg_yyy_xxxyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyy_xxxzzzz_1[j];

                    tg_yyyyy_xxxzzzzz_0[j] = pb_y * tg_yyyy_xxxzzzzz_0[j] + fr * tg_yyyy_xxxzzzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxzzzzz_0[j] - tg_yyy_xxxzzzzz_1[j] * fl1_fza);

                    tg_yyyyy_xxyyyyyy_0[j] = pb_y * tg_yyyy_xxyyyyyy_0[j] + fr * tg_yyyy_xxyyyyyy_1[j] + 2.0 * fl1_fx * (tg_yyy_xxyyyyyy_0[j] - tg_yyy_xxyyyyyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yyyy_xxyyyyy_1[j];

                    tg_yyyyy_xxyyyyyz_0[j] = pb_y * tg_yyyy_xxyyyyyz_0[j] + fr * tg_yyyy_xxyyyyyz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxyyyyyz_0[j] - tg_yyy_xxyyyyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyy_xxyyyyz_1[j];

                    tg_yyyyy_xxyyyyzz_0[j] = pb_y * tg_yyyy_xxyyyyzz_0[j] + fr * tg_yyyy_xxyyyyzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxyyyyzz_0[j] - tg_yyy_xxyyyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyy_xxyyyzz_1[j];

                    tg_yyyyy_xxyyyzzz_0[j] = pb_y * tg_yyyy_xxyyyzzz_0[j] + fr * tg_yyyy_xxyyyzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxyyyzzz_0[j] - tg_yyy_xxyyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyy_xxyyzzz_1[j];

                    tg_yyyyy_xxyyzzzz_0[j] = pb_y * tg_yyyy_xxyyzzzz_0[j] + fr * tg_yyyy_xxyyzzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxyyzzzz_0[j] - tg_yyy_xxyyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyy_xxyzzzz_1[j];

                    tg_yyyyy_xxyzzzzz_0[j] = pb_y * tg_yyyy_xxyzzzzz_0[j] + fr * tg_yyyy_xxyzzzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxyzzzzz_0[j] - tg_yyy_xxyzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyy_xxzzzzz_1[j];

                    tg_yyyyy_xxzzzzzz_0[j] = pb_y * tg_yyyy_xxzzzzzz_0[j] + fr * tg_yyyy_xxzzzzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxzzzzzz_0[j] - tg_yyy_xxzzzzzz_1[j] * fl1_fza);

                    tg_yyyyy_xyyyyyyy_0[j] = pb_y * tg_yyyy_xyyyyyyy_0[j] + fr * tg_yyyy_xyyyyyyy_1[j] + 2.0 * fl1_fx * (tg_yyy_xyyyyyyy_0[j] - tg_yyy_xyyyyyyy_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_yyyy_xyyyyyy_1[j];

                    tg_yyyyy_xyyyyyyz_0[j] = pb_y * tg_yyyy_xyyyyyyz_0[j] + fr * tg_yyyy_xyyyyyyz_1[j] + 2.0 * fl1_fx * (tg_yyy_xyyyyyyz_0[j] - tg_yyy_xyyyyyyz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yyyy_xyyyyyz_1[j];

                    tg_yyyyy_xyyyyyzz_0[j] = pb_y * tg_yyyy_xyyyyyzz_0[j] + fr * tg_yyyy_xyyyyyzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xyyyyyzz_0[j] - tg_yyy_xyyyyyzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyy_xyyyyzz_1[j];

                    tg_yyyyy_xyyyyzzz_0[j] = pb_y * tg_yyyy_xyyyyzzz_0[j] + fr * tg_yyyy_xyyyyzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xyyyyzzz_0[j] - tg_yyy_xyyyyzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyy_xyyyzzz_1[j];

                    tg_yyyyy_xyyyzzzz_0[j] = pb_y * tg_yyyy_xyyyzzzz_0[j] + fr * tg_yyyy_xyyyzzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xyyyzzzz_0[j] - tg_yyy_xyyyzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyy_xyyzzzz_1[j];

                    tg_yyyyy_xyyzzzzz_0[j] = pb_y * tg_yyyy_xyyzzzzz_0[j] + fr * tg_yyyy_xyyzzzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xyyzzzzz_0[j] - tg_yyy_xyyzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyy_xyzzzzz_1[j];

                    tg_yyyyy_xyzzzzzz_0[j] = pb_y * tg_yyyy_xyzzzzzz_0[j] + fr * tg_yyyy_xyzzzzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xyzzzzzz_0[j] - tg_yyy_xyzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyy_xzzzzzz_1[j];

                    tg_yyyyy_xzzzzzzz_0[j] = pb_y * tg_yyyy_xzzzzzzz_0[j] + fr * tg_yyyy_xzzzzzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xzzzzzzz_0[j] - tg_yyy_xzzzzzzz_1[j] * fl1_fza);

                    tg_yyyyy_yyyyyyyy_0[j] = pb_y * tg_yyyy_yyyyyyyy_0[j] + fr * tg_yyyy_yyyyyyyy_1[j] + 2.0 * fl1_fx * (tg_yyy_yyyyyyyy_0[j] - tg_yyy_yyyyyyyy_1[j] * fl1_fza) + 4.0 * fl1_fxn * tg_yyyy_yyyyyyy_1[j];

                    tg_yyyyy_yyyyyyyz_0[j] = pb_y * tg_yyyy_yyyyyyyz_0[j] + fr * tg_yyyy_yyyyyyyz_1[j] + 2.0 * fl1_fx * (tg_yyy_yyyyyyyz_0[j] - tg_yyy_yyyyyyyz_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_yyyy_yyyyyyz_1[j];

                    tg_yyyyy_yyyyyyzz_0[j] = pb_y * tg_yyyy_yyyyyyzz_0[j] + fr * tg_yyyy_yyyyyyzz_1[j] + 2.0 * fl1_fx * (tg_yyy_yyyyyyzz_0[j] - tg_yyy_yyyyyyzz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yyyy_yyyyyzz_1[j];

                    tg_yyyyy_yyyyyzzz_0[j] = pb_y * tg_yyyy_yyyyyzzz_0[j] + fr * tg_yyyy_yyyyyzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_yyyyyzzz_0[j] - tg_yyy_yyyyyzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyy_yyyyzzz_1[j];

                    tg_yyyyy_yyyyzzzz_0[j] = pb_y * tg_yyyy_yyyyzzzz_0[j] + fr * tg_yyyy_yyyyzzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_yyyyzzzz_0[j] - tg_yyy_yyyyzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyy_yyyzzzz_1[j];

                    tg_yyyyy_yyyzzzzz_0[j] = pb_y * tg_yyyy_yyyzzzzz_0[j] + fr * tg_yyyy_yyyzzzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_yyyzzzzz_0[j] - tg_yyy_yyyzzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyy_yyzzzzz_1[j];

                    tg_yyyyy_yyzzzzzz_0[j] = pb_y * tg_yyyy_yyzzzzzz_0[j] + fr * tg_yyyy_yyzzzzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_yyzzzzzz_0[j] - tg_yyy_yyzzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyy_yzzzzzz_1[j];

                    tg_yyyyy_yzzzzzzz_0[j] = pb_y * tg_yyyy_yzzzzzzz_0[j] + fr * tg_yyyy_yzzzzzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_yzzzzzzz_0[j] - tg_yyy_yzzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyy_zzzzzzz_1[j];

                    tg_yyyyy_zzzzzzzz_0[j] = pb_y * tg_yyyy_zzzzzzzz_0[j] + fr * tg_yyyy_zzzzzzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_zzzzzzzz_0[j] - tg_yyy_zzzzzzzz_1[j] * fl1_fza);

                    tg_yyyyz_xxxxxxxx_0[j] = pb_y * tg_yyyz_xxxxxxxx_0[j] + fr * tg_yyyz_xxxxxxxx_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxxxxxx_0[j] - tg_yyz_xxxxxxxx_1[j] * fl1_fza);

                    tg_yyyyz_xxxxxxxy_0[j] = pb_y * tg_yyyz_xxxxxxxy_0[j] + fr * tg_yyyz_xxxxxxxy_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxxxxxy_0[j] - tg_yyz_xxxxxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyz_xxxxxxx_1[j];

                    tg_yyyyz_xxxxxxxz_0[j] = pb_y * tg_yyyz_xxxxxxxz_0[j] + fr * tg_yyyz_xxxxxxxz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxxxxxz_0[j] - tg_yyz_xxxxxxxz_1[j] * fl1_fza);

                    tg_yyyyz_xxxxxxyy_0[j] = pb_y * tg_yyyz_xxxxxxyy_0[j] + fr * tg_yyyz_xxxxxxyy_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxxxxyy_0[j] - tg_yyz_xxxxxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyz_xxxxxxy_1[j];

                    tg_yyyyz_xxxxxxyz_0[j] = pb_y * tg_yyyz_xxxxxxyz_0[j] + fr * tg_yyyz_xxxxxxyz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxxxxyz_0[j] - tg_yyz_xxxxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyz_xxxxxxz_1[j];

                    tg_yyyyz_xxxxxxzz_0[j] = pb_y * tg_yyyz_xxxxxxzz_0[j] + fr * tg_yyyz_xxxxxxzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxxxxzz_0[j] - tg_yyz_xxxxxxzz_1[j] * fl1_fza);

                    tg_yyyyz_xxxxxyyy_0[j] = pb_y * tg_yyyz_xxxxxyyy_0[j] + fr * tg_yyyz_xxxxxyyy_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxxxyyy_0[j] - tg_yyz_xxxxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyz_xxxxxyy_1[j];

                    tg_yyyyz_xxxxxyyz_0[j] = pb_y * tg_yyyz_xxxxxyyz_0[j] + fr * tg_yyyz_xxxxxyyz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxxxyyz_0[j] - tg_yyz_xxxxxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyz_xxxxxyz_1[j];

                    tg_yyyyz_xxxxxyzz_0[j] = pb_y * tg_yyyz_xxxxxyzz_0[j] + fr * tg_yyyz_xxxxxyzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxxxyzz_0[j] - tg_yyz_xxxxxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyz_xxxxxzz_1[j];

                    tg_yyyyz_xxxxxzzz_0[j] = pb_y * tg_yyyz_xxxxxzzz_0[j] + fr * tg_yyyz_xxxxxzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxxxzzz_0[j] - tg_yyz_xxxxxzzz_1[j] * fl1_fza);

                    tg_yyyyz_xxxxyyyy_0[j] = pb_y * tg_yyyz_xxxxyyyy_0[j] + fr * tg_yyyz_xxxxyyyy_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxxyyyy_0[j] - tg_yyz_xxxxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyz_xxxxyyy_1[j];

                    tg_yyyyz_xxxxyyyz_0[j] = pb_y * tg_yyyz_xxxxyyyz_0[j] + fr * tg_yyyz_xxxxyyyz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxxyyyz_0[j] - tg_yyz_xxxxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyz_xxxxyyz_1[j];

                    tg_yyyyz_xxxxyyzz_0[j] = pb_y * tg_yyyz_xxxxyyzz_0[j] + fr * tg_yyyz_xxxxyyzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxxyyzz_0[j] - tg_yyz_xxxxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyz_xxxxyzz_1[j];

                    tg_yyyyz_xxxxyzzz_0[j] = pb_y * tg_yyyz_xxxxyzzz_0[j] + fr * tg_yyyz_xxxxyzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxxyzzz_0[j] - tg_yyz_xxxxyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyz_xxxxzzz_1[j];

                    tg_yyyyz_xxxxzzzz_0[j] = pb_y * tg_yyyz_xxxxzzzz_0[j] + fr * tg_yyyz_xxxxzzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxxzzzz_0[j] - tg_yyz_xxxxzzzz_1[j] * fl1_fza);

                    tg_yyyyz_xxxyyyyy_0[j] = pb_y * tg_yyyz_xxxyyyyy_0[j] + fr * tg_yyyz_xxxyyyyy_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxyyyyy_0[j] - tg_yyz_xxxyyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyz_xxxyyyy_1[j];

                    tg_yyyyz_xxxyyyyz_0[j] = pb_y * tg_yyyz_xxxyyyyz_0[j] + fr * tg_yyyz_xxxyyyyz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxyyyyz_0[j] - tg_yyz_xxxyyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyz_xxxyyyz_1[j];

                    tg_yyyyz_xxxyyyzz_0[j] = pb_y * tg_yyyz_xxxyyyzz_0[j] + fr * tg_yyyz_xxxyyyzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxyyyzz_0[j] - tg_yyz_xxxyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyz_xxxyyzz_1[j];

                    tg_yyyyz_xxxyyzzz_0[j] = pb_y * tg_yyyz_xxxyyzzz_0[j] + fr * tg_yyyz_xxxyyzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxyyzzz_0[j] - tg_yyz_xxxyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyz_xxxyzzz_1[j];

                    tg_yyyyz_xxxyzzzz_0[j] = pb_y * tg_yyyz_xxxyzzzz_0[j] + fr * tg_yyyz_xxxyzzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxyzzzz_0[j] - tg_yyz_xxxyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyz_xxxzzzz_1[j];

                    tg_yyyyz_xxxzzzzz_0[j] = pb_y * tg_yyyz_xxxzzzzz_0[j] + fr * tg_yyyz_xxxzzzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxzzzzz_0[j] - tg_yyz_xxxzzzzz_1[j] * fl1_fza);

                    tg_yyyyz_xxyyyyyy_0[j] = pb_y * tg_yyyz_xxyyyyyy_0[j] + fr * tg_yyyz_xxyyyyyy_1[j] + 1.5 * fl1_fx * (tg_yyz_xxyyyyyy_0[j] - tg_yyz_xxyyyyyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yyyz_xxyyyyy_1[j];

                    tg_yyyyz_xxyyyyyz_0[j] = pb_y * tg_yyyz_xxyyyyyz_0[j] + fr * tg_yyyz_xxyyyyyz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxyyyyyz_0[j] - tg_yyz_xxyyyyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyz_xxyyyyz_1[j];

                    tg_yyyyz_xxyyyyzz_0[j] = pb_y * tg_yyyz_xxyyyyzz_0[j] + fr * tg_yyyz_xxyyyyzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxyyyyzz_0[j] - tg_yyz_xxyyyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyz_xxyyyzz_1[j];

                    tg_yyyyz_xxyyyzzz_0[j] = pb_y * tg_yyyz_xxyyyzzz_0[j] + fr * tg_yyyz_xxyyyzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxyyyzzz_0[j] - tg_yyz_xxyyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyz_xxyyzzz_1[j];

                    tg_yyyyz_xxyyzzzz_0[j] = pb_y * tg_yyyz_xxyyzzzz_0[j] + fr * tg_yyyz_xxyyzzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxyyzzzz_0[j] - tg_yyz_xxyyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyz_xxyzzzz_1[j];

                    tg_yyyyz_xxyzzzzz_0[j] = pb_y * tg_yyyz_xxyzzzzz_0[j] + fr * tg_yyyz_xxyzzzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxyzzzzz_0[j] - tg_yyz_xxyzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyz_xxzzzzz_1[j];

                    tg_yyyyz_xxzzzzzz_0[j] = pb_y * tg_yyyz_xxzzzzzz_0[j] + fr * tg_yyyz_xxzzzzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxzzzzzz_0[j] - tg_yyz_xxzzzzzz_1[j] * fl1_fza);

                    tg_yyyyz_xyyyyyyy_0[j] = pb_y * tg_yyyz_xyyyyyyy_0[j] + fr * tg_yyyz_xyyyyyyy_1[j] + 1.5 * fl1_fx * (tg_yyz_xyyyyyyy_0[j] - tg_yyz_xyyyyyyy_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_yyyz_xyyyyyy_1[j];

                    tg_yyyyz_xyyyyyyz_0[j] = pb_y * tg_yyyz_xyyyyyyz_0[j] + fr * tg_yyyz_xyyyyyyz_1[j] + 1.5 * fl1_fx * (tg_yyz_xyyyyyyz_0[j] - tg_yyz_xyyyyyyz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yyyz_xyyyyyz_1[j];

                    tg_yyyyz_xyyyyyzz_0[j] = pb_y * tg_yyyz_xyyyyyzz_0[j] + fr * tg_yyyz_xyyyyyzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xyyyyyzz_0[j] - tg_yyz_xyyyyyzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyz_xyyyyzz_1[j];

                    tg_yyyyz_xyyyyzzz_0[j] = pb_y * tg_yyyz_xyyyyzzz_0[j] + fr * tg_yyyz_xyyyyzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xyyyyzzz_0[j] - tg_yyz_xyyyyzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyz_xyyyzzz_1[j];

                    tg_yyyyz_xyyyzzzz_0[j] = pb_y * tg_yyyz_xyyyzzzz_0[j] + fr * tg_yyyz_xyyyzzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xyyyzzzz_0[j] - tg_yyz_xyyyzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyz_xyyzzzz_1[j];

                    tg_yyyyz_xyyzzzzz_0[j] = pb_y * tg_yyyz_xyyzzzzz_0[j] + fr * tg_yyyz_xyyzzzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xyyzzzzz_0[j] - tg_yyz_xyyzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyz_xyzzzzz_1[j];

                    tg_yyyyz_xyzzzzzz_0[j] = pb_y * tg_yyyz_xyzzzzzz_0[j] + fr * tg_yyyz_xyzzzzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xyzzzzzz_0[j] - tg_yyz_xyzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyz_xzzzzzz_1[j];

                    tg_yyyyz_xzzzzzzz_0[j] = pb_y * tg_yyyz_xzzzzzzz_0[j] + fr * tg_yyyz_xzzzzzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xzzzzzzz_0[j] - tg_yyz_xzzzzzzz_1[j] * fl1_fza);

                    tg_yyyyz_yyyyyyyy_0[j] = pb_y * tg_yyyz_yyyyyyyy_0[j] + fr * tg_yyyz_yyyyyyyy_1[j] + 1.5 * fl1_fx * (tg_yyz_yyyyyyyy_0[j] - tg_yyz_yyyyyyyy_1[j] * fl1_fza) + 4.0 * fl1_fxn * tg_yyyz_yyyyyyy_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSL_757_851(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (757,851)

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
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_8_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_8_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_8_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_8_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_8_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_7_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {7, -1, -1, -1}, 
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

                auto tg_yyyz_yyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 532); 

                auto tg_yyyz_yyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 533); 

                auto tg_yyyz_yyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 534); 

                auto tg_yyyz_yyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 535); 

                auto tg_yyyz_yyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 536); 

                auto tg_yyyz_yyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 537); 

                auto tg_yyyz_yzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 538); 

                auto tg_yyyz_zzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 539); 

                auto tg_yyzz_xxxxxxxx_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 540); 

                auto tg_yyzz_xxxxxxxy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 541); 

                auto tg_yyzz_xxxxxxxz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 542); 

                auto tg_yyzz_xxxxxxyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 543); 

                auto tg_yyzz_xxxxxxyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 544); 

                auto tg_yyzz_xxxxxxzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 545); 

                auto tg_yyzz_xxxxxyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 546); 

                auto tg_yyzz_xxxxxyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 547); 

                auto tg_yyzz_xxxxxyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 548); 

                auto tg_yyzz_xxxxxzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 549); 

                auto tg_yyzz_xxxxyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 550); 

                auto tg_yyzz_xxxxyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 551); 

                auto tg_yyzz_xxxxyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 552); 

                auto tg_yyzz_xxxxyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 553); 

                auto tg_yyzz_xxxxzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 554); 

                auto tg_yyzz_xxxyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 555); 

                auto tg_yyzz_xxxyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 556); 

                auto tg_yyzz_xxxyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 557); 

                auto tg_yyzz_xxxyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 558); 

                auto tg_yyzz_xxxyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 559); 

                auto tg_yyzz_xxxzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 560); 

                auto tg_yyzz_xxyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 561); 

                auto tg_yyzz_xxyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 562); 

                auto tg_yyzz_xxyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 563); 

                auto tg_yyzz_xxyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 564); 

                auto tg_yyzz_xxyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 565); 

                auto tg_yyzz_xxyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 566); 

                auto tg_yyzz_xxzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 567); 

                auto tg_yyzz_xyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 568); 

                auto tg_yyzz_xyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 569); 

                auto tg_yyzz_xyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 570); 

                auto tg_yyzz_xyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 571); 

                auto tg_yyzz_xyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 572); 

                auto tg_yyzz_xyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 573); 

                auto tg_yyzz_xyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 574); 

                auto tg_yyzz_xzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 575); 

                auto tg_yyzz_yyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 576); 

                auto tg_yyzz_yyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 577); 

                auto tg_yyzz_yyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 578); 

                auto tg_yyzz_yyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 579); 

                auto tg_yyzz_yyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 580); 

                auto tg_yyzz_yyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 581); 

                auto tg_yyzz_yyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 582); 

                auto tg_yyzz_yzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 583); 

                auto tg_yyzz_zzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 584); 

                auto tg_yzzz_xxxxxxxx_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 585); 

                auto tg_yzzz_xxxxxxxy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 586); 

                auto tg_yzzz_xxxxxxxz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 587); 

                auto tg_yzzz_xxxxxxyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 588); 

                auto tg_yzzz_xxxxxxyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 589); 

                auto tg_yzzz_xxxxxxzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 590); 

                auto tg_yzzz_xxxxxyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 591); 

                auto tg_yzzz_xxxxxyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 592); 

                auto tg_yzzz_xxxxxyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 593); 

                auto tg_yzzz_xxxxxzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 594); 

                auto tg_yzzz_xxxxyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 595); 

                auto tg_yzzz_xxxxyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 596); 

                auto tg_yzzz_xxxxyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 597); 

                auto tg_yzzz_xxxxyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 598); 

                auto tg_yzzz_xxxxzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 599); 

                auto tg_yzzz_xxxyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 600); 

                auto tg_yzzz_xxxyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 601); 

                auto tg_yzzz_xxxyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 602); 

                auto tg_yzzz_xxxyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 603); 

                auto tg_yzzz_xxxyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 604); 

                auto tg_yzzz_xxxzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 605); 

                auto tg_yzzz_xxyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 606); 

                auto tg_yzzz_xxyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 607); 

                auto tg_yzzz_xxyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 608); 

                auto tg_yzzz_xxyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 609); 

                auto tg_yzzz_xxyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 610); 

                auto tg_yzzz_xxyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 611); 

                auto tg_yzzz_xxzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 612); 

                auto tg_yzzz_xyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 613); 

                auto tg_yzzz_xyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 614); 

                auto tg_yzzz_xyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 615); 

                auto tg_yzzz_xyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 616); 

                auto tg_yzzz_xyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 617); 

                auto tg_yzzz_xyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 618); 

                auto tg_yzzz_xyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 619); 

                auto tg_yzzz_xzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 620); 

                auto tg_yzzz_yyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 621); 

                auto tg_yzzz_yyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 622); 

                auto tg_yzzz_yyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 623); 

                auto tg_yzzz_yyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 624); 

                auto tg_yzzz_yyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 625); 

                auto tg_yyyz_yyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 532); 

                auto tg_yyyz_yyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 533); 

                auto tg_yyyz_yyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 534); 

                auto tg_yyyz_yyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 535); 

                auto tg_yyyz_yyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 536); 

                auto tg_yyyz_yyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 537); 

                auto tg_yyyz_yzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 538); 

                auto tg_yyyz_zzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 539); 

                auto tg_yyzz_xxxxxxxx_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 540); 

                auto tg_yyzz_xxxxxxxy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 541); 

                auto tg_yyzz_xxxxxxxz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 542); 

                auto tg_yyzz_xxxxxxyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 543); 

                auto tg_yyzz_xxxxxxyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 544); 

                auto tg_yyzz_xxxxxxzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 545); 

                auto tg_yyzz_xxxxxyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 546); 

                auto tg_yyzz_xxxxxyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 547); 

                auto tg_yyzz_xxxxxyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 548); 

                auto tg_yyzz_xxxxxzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 549); 

                auto tg_yyzz_xxxxyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 550); 

                auto tg_yyzz_xxxxyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 551); 

                auto tg_yyzz_xxxxyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 552); 

                auto tg_yyzz_xxxxyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 553); 

                auto tg_yyzz_xxxxzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 554); 

                auto tg_yyzz_xxxyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 555); 

                auto tg_yyzz_xxxyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 556); 

                auto tg_yyzz_xxxyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 557); 

                auto tg_yyzz_xxxyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 558); 

                auto tg_yyzz_xxxyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 559); 

                auto tg_yyzz_xxxzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 560); 

                auto tg_yyzz_xxyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 561); 

                auto tg_yyzz_xxyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 562); 

                auto tg_yyzz_xxyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 563); 

                auto tg_yyzz_xxyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 564); 

                auto tg_yyzz_xxyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 565); 

                auto tg_yyzz_xxyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 566); 

                auto tg_yyzz_xxzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 567); 

                auto tg_yyzz_xyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 568); 

                auto tg_yyzz_xyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 569); 

                auto tg_yyzz_xyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 570); 

                auto tg_yyzz_xyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 571); 

                auto tg_yyzz_xyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 572); 

                auto tg_yyzz_xyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 573); 

                auto tg_yyzz_xyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 574); 

                auto tg_yyzz_xzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 575); 

                auto tg_yyzz_yyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 576); 

                auto tg_yyzz_yyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 577); 

                auto tg_yyzz_yyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 578); 

                auto tg_yyzz_yyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 579); 

                auto tg_yyzz_yyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 580); 

                auto tg_yyzz_yyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 581); 

                auto tg_yyzz_yyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 582); 

                auto tg_yyzz_yzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 583); 

                auto tg_yyzz_zzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 584); 

                auto tg_yzzz_xxxxxxxx_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 585); 

                auto tg_yzzz_xxxxxxxy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 586); 

                auto tg_yzzz_xxxxxxxz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 587); 

                auto tg_yzzz_xxxxxxyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 588); 

                auto tg_yzzz_xxxxxxyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 589); 

                auto tg_yzzz_xxxxxxzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 590); 

                auto tg_yzzz_xxxxxyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 591); 

                auto tg_yzzz_xxxxxyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 592); 

                auto tg_yzzz_xxxxxyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 593); 

                auto tg_yzzz_xxxxxzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 594); 

                auto tg_yzzz_xxxxyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 595); 

                auto tg_yzzz_xxxxyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 596); 

                auto tg_yzzz_xxxxyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 597); 

                auto tg_yzzz_xxxxyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 598); 

                auto tg_yzzz_xxxxzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 599); 

                auto tg_yzzz_xxxyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 600); 

                auto tg_yzzz_xxxyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 601); 

                auto tg_yzzz_xxxyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 602); 

                auto tg_yzzz_xxxyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 603); 

                auto tg_yzzz_xxxyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 604); 

                auto tg_yzzz_xxxzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 605); 

                auto tg_yzzz_xxyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 606); 

                auto tg_yzzz_xxyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 607); 

                auto tg_yzzz_xxyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 608); 

                auto tg_yzzz_xxyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 609); 

                auto tg_yzzz_xxyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 610); 

                auto tg_yzzz_xxyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 611); 

                auto tg_yzzz_xxzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 612); 

                auto tg_yzzz_xyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 613); 

                auto tg_yzzz_xyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 614); 

                auto tg_yzzz_xyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 615); 

                auto tg_yzzz_xyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 616); 

                auto tg_yzzz_xyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 617); 

                auto tg_yzzz_xyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 618); 

                auto tg_yzzz_xyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 619); 

                auto tg_yzzz_xzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 620); 

                auto tg_yzzz_yyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 621); 

                auto tg_yzzz_yyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 622); 

                auto tg_yzzz_yyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 623); 

                auto tg_yzzz_yyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 624); 

                auto tg_yzzz_yyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 625); 

                auto tg_yyz_yyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 352); 

                auto tg_yyz_yyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 353); 

                auto tg_yyz_yyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 354); 

                auto tg_yyz_yyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 355); 

                auto tg_yyz_yyyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 356); 

                auto tg_yyz_yyzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 357); 

                auto tg_yyz_yzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 358); 

                auto tg_yyz_zzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 359); 

                auto tg_yzz_xxxxxxxx_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 360); 

                auto tg_yzz_xxxxxxxy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 361); 

                auto tg_yzz_xxxxxxxz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 362); 

                auto tg_yzz_xxxxxxyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 363); 

                auto tg_yzz_xxxxxxyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 364); 

                auto tg_yzz_xxxxxxzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 365); 

                auto tg_yzz_xxxxxyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 366); 

                auto tg_yzz_xxxxxyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 367); 

                auto tg_yzz_xxxxxyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 368); 

                auto tg_yzz_xxxxxzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 369); 

                auto tg_yzz_xxxxyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 370); 

                auto tg_yzz_xxxxyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 371); 

                auto tg_yzz_xxxxyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 372); 

                auto tg_yzz_xxxxyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 373); 

                auto tg_yzz_xxxxzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 374); 

                auto tg_yzz_xxxyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 375); 

                auto tg_yzz_xxxyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 376); 

                auto tg_yzz_xxxyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 377); 

                auto tg_yzz_xxxyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 378); 

                auto tg_yzz_xxxyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 379); 

                auto tg_yzz_xxxzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 380); 

                auto tg_yzz_xxyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 381); 

                auto tg_yzz_xxyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 382); 

                auto tg_yzz_xxyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 383); 

                auto tg_yzz_xxyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 384); 

                auto tg_yzz_xxyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 385); 

                auto tg_yzz_xxyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 386); 

                auto tg_yzz_xxzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 387); 

                auto tg_yzz_xyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 388); 

                auto tg_yzz_xyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 389); 

                auto tg_yzz_xyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 390); 

                auto tg_yzz_xyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 391); 

                auto tg_yzz_xyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 392); 

                auto tg_yzz_xyyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 393); 

                auto tg_yzz_xyzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 394); 

                auto tg_yzz_xzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 395); 

                auto tg_yzz_yyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 396); 

                auto tg_yzz_yyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 397); 

                auto tg_yzz_yyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 398); 

                auto tg_yzz_yyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 399); 

                auto tg_yzz_yyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 400); 

                auto tg_yzz_yyyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 401); 

                auto tg_yzz_yyzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 402); 

                auto tg_yzz_yzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 403); 

                auto tg_yzz_zzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 404); 

                auto tg_zzz_xxxxxxxx_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 405); 

                auto tg_zzz_xxxxxxxy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 406); 

                auto tg_zzz_xxxxxxxz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 407); 

                auto tg_zzz_xxxxxxyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 408); 

                auto tg_zzz_xxxxxxyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 409); 

                auto tg_zzz_xxxxxxzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 410); 

                auto tg_zzz_xxxxxyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 411); 

                auto tg_zzz_xxxxxyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 412); 

                auto tg_zzz_xxxxxyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 413); 

                auto tg_zzz_xxxxxzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 414); 

                auto tg_zzz_xxxxyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 415); 

                auto tg_zzz_xxxxyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 416); 

                auto tg_zzz_xxxxyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 417); 

                auto tg_zzz_xxxxyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 418); 

                auto tg_zzz_xxxxzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 419); 

                auto tg_zzz_xxxyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 420); 

                auto tg_zzz_xxxyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 421); 

                auto tg_zzz_xxxyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 422); 

                auto tg_zzz_xxxyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 423); 

                auto tg_zzz_xxxyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 424); 

                auto tg_zzz_xxxzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 425); 

                auto tg_zzz_xxyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 426); 

                auto tg_zzz_xxyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 427); 

                auto tg_zzz_xxyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 428); 

                auto tg_zzz_xxyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 429); 

                auto tg_zzz_xxyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 430); 

                auto tg_zzz_xxyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 431); 

                auto tg_zzz_xxzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 432); 

                auto tg_zzz_xyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 433); 

                auto tg_zzz_xyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 434); 

                auto tg_zzz_xyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 435); 

                auto tg_zzz_xyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 436); 

                auto tg_zzz_xyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 437); 

                auto tg_zzz_xyyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 438); 

                auto tg_zzz_xyzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 439); 

                auto tg_zzz_xzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 440); 

                auto tg_zzz_yyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 441); 

                auto tg_zzz_yyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 442); 

                auto tg_zzz_yyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 443); 

                auto tg_zzz_yyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 444); 

                auto tg_zzz_yyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 445); 

                auto tg_yyz_yyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 352); 

                auto tg_yyz_yyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 353); 

                auto tg_yyz_yyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 354); 

                auto tg_yyz_yyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 355); 

                auto tg_yyz_yyyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 356); 

                auto tg_yyz_yyzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 357); 

                auto tg_yyz_yzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 358); 

                auto tg_yyz_zzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 359); 

                auto tg_yzz_xxxxxxxx_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 360); 

                auto tg_yzz_xxxxxxxy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 361); 

                auto tg_yzz_xxxxxxxz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 362); 

                auto tg_yzz_xxxxxxyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 363); 

                auto tg_yzz_xxxxxxyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 364); 

                auto tg_yzz_xxxxxxzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 365); 

                auto tg_yzz_xxxxxyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 366); 

                auto tg_yzz_xxxxxyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 367); 

                auto tg_yzz_xxxxxyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 368); 

                auto tg_yzz_xxxxxzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 369); 

                auto tg_yzz_xxxxyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 370); 

                auto tg_yzz_xxxxyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 371); 

                auto tg_yzz_xxxxyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 372); 

                auto tg_yzz_xxxxyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 373); 

                auto tg_yzz_xxxxzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 374); 

                auto tg_yzz_xxxyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 375); 

                auto tg_yzz_xxxyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 376); 

                auto tg_yzz_xxxyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 377); 

                auto tg_yzz_xxxyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 378); 

                auto tg_yzz_xxxyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 379); 

                auto tg_yzz_xxxzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 380); 

                auto tg_yzz_xxyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 381); 

                auto tg_yzz_xxyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 382); 

                auto tg_yzz_xxyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 383); 

                auto tg_yzz_xxyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 384); 

                auto tg_yzz_xxyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 385); 

                auto tg_yzz_xxyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 386); 

                auto tg_yzz_xxzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 387); 

                auto tg_yzz_xyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 388); 

                auto tg_yzz_xyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 389); 

                auto tg_yzz_xyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 390); 

                auto tg_yzz_xyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 391); 

                auto tg_yzz_xyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 392); 

                auto tg_yzz_xyyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 393); 

                auto tg_yzz_xyzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 394); 

                auto tg_yzz_xzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 395); 

                auto tg_yzz_yyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 396); 

                auto tg_yzz_yyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 397); 

                auto tg_yzz_yyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 398); 

                auto tg_yzz_yyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 399); 

                auto tg_yzz_yyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 400); 

                auto tg_yzz_yyyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 401); 

                auto tg_yzz_yyzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 402); 

                auto tg_yzz_yzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 403); 

                auto tg_yzz_zzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 404); 

                auto tg_zzz_xxxxxxxx_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 405); 

                auto tg_zzz_xxxxxxxy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 406); 

                auto tg_zzz_xxxxxxxz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 407); 

                auto tg_zzz_xxxxxxyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 408); 

                auto tg_zzz_xxxxxxyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 409); 

                auto tg_zzz_xxxxxxzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 410); 

                auto tg_zzz_xxxxxyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 411); 

                auto tg_zzz_xxxxxyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 412); 

                auto tg_zzz_xxxxxyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 413); 

                auto tg_zzz_xxxxxzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 414); 

                auto tg_zzz_xxxxyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 415); 

                auto tg_zzz_xxxxyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 416); 

                auto tg_zzz_xxxxyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 417); 

                auto tg_zzz_xxxxyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 418); 

                auto tg_zzz_xxxxzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 419); 

                auto tg_zzz_xxxyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 420); 

                auto tg_zzz_xxxyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 421); 

                auto tg_zzz_xxxyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 422); 

                auto tg_zzz_xxxyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 423); 

                auto tg_zzz_xxxyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 424); 

                auto tg_zzz_xxxzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 425); 

                auto tg_zzz_xxyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 426); 

                auto tg_zzz_xxyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 427); 

                auto tg_zzz_xxyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 428); 

                auto tg_zzz_xxyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 429); 

                auto tg_zzz_xxyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 430); 

                auto tg_zzz_xxyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 431); 

                auto tg_zzz_xxzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 432); 

                auto tg_zzz_xyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 433); 

                auto tg_zzz_xyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 434); 

                auto tg_zzz_xyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 435); 

                auto tg_zzz_xyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 436); 

                auto tg_zzz_xyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 437); 

                auto tg_zzz_xyyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 438); 

                auto tg_zzz_xyzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 439); 

                auto tg_zzz_xzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 440); 

                auto tg_zzz_yyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 441); 

                auto tg_zzz_yyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 442); 

                auto tg_zzz_yyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 443); 

                auto tg_zzz_yyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 444); 

                auto tg_zzz_yyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 445); 

                auto tg_yyyz_yyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 425); 

                auto tg_yyyz_yyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 426); 

                auto tg_yyyz_yyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 427); 

                auto tg_yyyz_yyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 428); 

                auto tg_yyyz_yyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 429); 

                auto tg_yyyz_yzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 430); 

                auto tg_yyyz_zzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 431); 

                auto tg_yyzz_xxxxxxx_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 432); 

                auto tg_yyzz_xxxxxxy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 433); 

                auto tg_yyzz_xxxxxxz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 434); 

                auto tg_yyzz_xxxxxyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 435); 

                auto tg_yyzz_xxxxxyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 436); 

                auto tg_yyzz_xxxxxzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 437); 

                auto tg_yyzz_xxxxyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 438); 

                auto tg_yyzz_xxxxyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 439); 

                auto tg_yyzz_xxxxyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 440); 

                auto tg_yyzz_xxxxzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 441); 

                auto tg_yyzz_xxxyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 442); 

                auto tg_yyzz_xxxyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 443); 

                auto tg_yyzz_xxxyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 444); 

                auto tg_yyzz_xxxyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 445); 

                auto tg_yyzz_xxxzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 446); 

                auto tg_yyzz_xxyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 447); 

                auto tg_yyzz_xxyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 448); 

                auto tg_yyzz_xxyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 449); 

                auto tg_yyzz_xxyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 450); 

                auto tg_yyzz_xxyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 451); 

                auto tg_yyzz_xxzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 452); 

                auto tg_yyzz_xyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 453); 

                auto tg_yyzz_xyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 454); 

                auto tg_yyzz_xyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 455); 

                auto tg_yyzz_xyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 456); 

                auto tg_yyzz_xyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 457); 

                auto tg_yyzz_xyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 458); 

                auto tg_yyzz_xzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 459); 

                auto tg_yyzz_yyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 460); 

                auto tg_yyzz_yyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 461); 

                auto tg_yyzz_yyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 462); 

                auto tg_yyzz_yyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 463); 

                auto tg_yyzz_yyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 464); 

                auto tg_yyzz_yyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 465); 

                auto tg_yyzz_yzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 466); 

                auto tg_yyzz_zzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 467); 

                auto tg_yzzz_xxxxxxx_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 468); 

                auto tg_yzzz_xxxxxxy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 469); 

                auto tg_yzzz_xxxxxxz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 470); 

                auto tg_yzzz_xxxxxyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 471); 

                auto tg_yzzz_xxxxxyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 472); 

                auto tg_yzzz_xxxxxzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 473); 

                auto tg_yzzz_xxxxyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 474); 

                auto tg_yzzz_xxxxyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 475); 

                auto tg_yzzz_xxxxyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 476); 

                auto tg_yzzz_xxxxzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 477); 

                auto tg_yzzz_xxxyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 478); 

                auto tg_yzzz_xxxyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 479); 

                auto tg_yzzz_xxxyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 480); 

                auto tg_yzzz_xxxyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 481); 

                auto tg_yzzz_xxxzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 482); 

                auto tg_yzzz_xxyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 483); 

                auto tg_yzzz_xxyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 484); 

                auto tg_yzzz_xxyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 485); 

                auto tg_yzzz_xxyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 486); 

                auto tg_yzzz_xxyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 487); 

                auto tg_yzzz_xxzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 488); 

                auto tg_yzzz_xyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 489); 

                auto tg_yzzz_xyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 490); 

                auto tg_yzzz_xyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 491); 

                auto tg_yzzz_xyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 492); 

                auto tg_yzzz_xyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 493); 

                auto tg_yzzz_xyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 494); 

                auto tg_yzzz_xzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 495); 

                auto tg_yzzz_yyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 496); 

                auto tg_yzzz_yyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 497); 

                auto tg_yzzz_yyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 498); 

                auto tg_yzzz_yyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 499); 

                auto tg_yzzz_yyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 500); 

                // set up pointers to integrals

                auto tg_yyyyz_yyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 757); 

                auto tg_yyyyz_yyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 758); 

                auto tg_yyyyz_yyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 759); 

                auto tg_yyyyz_yyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 760); 

                auto tg_yyyyz_yyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 761); 

                auto tg_yyyyz_yyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 762); 

                auto tg_yyyyz_yzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 763); 

                auto tg_yyyyz_zzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 764); 

                auto tg_yyyzz_xxxxxxxx_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 765); 

                auto tg_yyyzz_xxxxxxxy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 766); 

                auto tg_yyyzz_xxxxxxxz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 767); 

                auto tg_yyyzz_xxxxxxyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 768); 

                auto tg_yyyzz_xxxxxxyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 769); 

                auto tg_yyyzz_xxxxxxzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 770); 

                auto tg_yyyzz_xxxxxyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 771); 

                auto tg_yyyzz_xxxxxyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 772); 

                auto tg_yyyzz_xxxxxyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 773); 

                auto tg_yyyzz_xxxxxzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 774); 

                auto tg_yyyzz_xxxxyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 775); 

                auto tg_yyyzz_xxxxyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 776); 

                auto tg_yyyzz_xxxxyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 777); 

                auto tg_yyyzz_xxxxyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 778); 

                auto tg_yyyzz_xxxxzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 779); 

                auto tg_yyyzz_xxxyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 780); 

                auto tg_yyyzz_xxxyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 781); 

                auto tg_yyyzz_xxxyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 782); 

                auto tg_yyyzz_xxxyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 783); 

                auto tg_yyyzz_xxxyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 784); 

                auto tg_yyyzz_xxxzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 785); 

                auto tg_yyyzz_xxyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 786); 

                auto tg_yyyzz_xxyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 787); 

                auto tg_yyyzz_xxyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 788); 

                auto tg_yyyzz_xxyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 789); 

                auto tg_yyyzz_xxyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 790); 

                auto tg_yyyzz_xxyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 791); 

                auto tg_yyyzz_xxzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 792); 

                auto tg_yyyzz_xyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 793); 

                auto tg_yyyzz_xyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 794); 

                auto tg_yyyzz_xyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 795); 

                auto tg_yyyzz_xyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 796); 

                auto tg_yyyzz_xyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 797); 

                auto tg_yyyzz_xyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 798); 

                auto tg_yyyzz_xyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 799); 

                auto tg_yyyzz_xzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 800); 

                auto tg_yyyzz_yyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 801); 

                auto tg_yyyzz_yyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 802); 

                auto tg_yyyzz_yyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 803); 

                auto tg_yyyzz_yyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 804); 

                auto tg_yyyzz_yyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 805); 

                auto tg_yyyzz_yyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 806); 

                auto tg_yyyzz_yyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 807); 

                auto tg_yyyzz_yzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 808); 

                auto tg_yyyzz_zzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 809); 

                auto tg_yyzzz_xxxxxxxx_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 810); 

                auto tg_yyzzz_xxxxxxxy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 811); 

                auto tg_yyzzz_xxxxxxxz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 812); 

                auto tg_yyzzz_xxxxxxyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 813); 

                auto tg_yyzzz_xxxxxxyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 814); 

                auto tg_yyzzz_xxxxxxzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 815); 

                auto tg_yyzzz_xxxxxyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 816); 

                auto tg_yyzzz_xxxxxyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 817); 

                auto tg_yyzzz_xxxxxyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 818); 

                auto tg_yyzzz_xxxxxzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 819); 

                auto tg_yyzzz_xxxxyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 820); 

                auto tg_yyzzz_xxxxyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 821); 

                auto tg_yyzzz_xxxxyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 822); 

                auto tg_yyzzz_xxxxyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 823); 

                auto tg_yyzzz_xxxxzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 824); 

                auto tg_yyzzz_xxxyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 825); 

                auto tg_yyzzz_xxxyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 826); 

                auto tg_yyzzz_xxxyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 827); 

                auto tg_yyzzz_xxxyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 828); 

                auto tg_yyzzz_xxxyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 829); 

                auto tg_yyzzz_xxxzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 830); 

                auto tg_yyzzz_xxyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 831); 

                auto tg_yyzzz_xxyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 832); 

                auto tg_yyzzz_xxyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 833); 

                auto tg_yyzzz_xxyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 834); 

                auto tg_yyzzz_xxyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 835); 

                auto tg_yyzzz_xxyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 836); 

                auto tg_yyzzz_xxzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 837); 

                auto tg_yyzzz_xyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 838); 

                auto tg_yyzzz_xyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 839); 

                auto tg_yyzzz_xyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 840); 

                auto tg_yyzzz_xyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 841); 

                auto tg_yyzzz_xyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 842); 

                auto tg_yyzzz_xyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 843); 

                auto tg_yyzzz_xyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 844); 

                auto tg_yyzzz_xzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 845); 

                auto tg_yyzzz_yyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 846); 

                auto tg_yyzzz_yyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 847); 

                auto tg_yyzzz_yyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 848); 

                auto tg_yyzzz_yyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 849); 

                auto tg_yyzzz_yyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 850); 

                // Batch of Integrals (757,851)

                #pragma omp simd aligned(fxn, fza, tg_yyyyz_yyyyyyyz_0, tg_yyyyz_yyyyyyzz_0, \
                                         tg_yyyyz_yyyyyzzz_0, tg_yyyyz_yyyyzzzz_0, tg_yyyyz_yyyzzzzz_0, tg_yyyyz_yyzzzzzz_0, \
                                         tg_yyyyz_yzzzzzzz_0, tg_yyyyz_zzzzzzzz_0, tg_yyyz_yyyyyyyz_0, tg_yyyz_yyyyyyyz_1, \
                                         tg_yyyz_yyyyyyz_1, tg_yyyz_yyyyyyzz_0, tg_yyyz_yyyyyyzz_1, tg_yyyz_yyyyyzz_1, \
                                         tg_yyyz_yyyyyzzz_0, tg_yyyz_yyyyyzzz_1, tg_yyyz_yyyyzzz_1, tg_yyyz_yyyyzzzz_0, \
                                         tg_yyyz_yyyyzzzz_1, tg_yyyz_yyyzzzz_1, tg_yyyz_yyyzzzzz_0, tg_yyyz_yyyzzzzz_1, \
                                         tg_yyyz_yyzzzzz_1, tg_yyyz_yyzzzzzz_0, tg_yyyz_yyzzzzzz_1, tg_yyyz_yzzzzzz_1, \
                                         tg_yyyz_yzzzzzzz_0, tg_yyyz_yzzzzzzz_1, tg_yyyz_zzzzzzz_1, tg_yyyz_zzzzzzzz_0, \
                                         tg_yyyz_zzzzzzzz_1, tg_yyyzz_xxxxxxxx_0, tg_yyyzz_xxxxxxxy_0, tg_yyyzz_xxxxxxxz_0, \
                                         tg_yyyzz_xxxxxxyy_0, tg_yyyzz_xxxxxxyz_0, tg_yyyzz_xxxxxxzz_0, tg_yyyzz_xxxxxyyy_0, \
                                         tg_yyyzz_xxxxxyyz_0, tg_yyyzz_xxxxxyzz_0, tg_yyyzz_xxxxxzzz_0, tg_yyyzz_xxxxyyyy_0, \
                                         tg_yyyzz_xxxxyyyz_0, tg_yyyzz_xxxxyyzz_0, tg_yyyzz_xxxxyzzz_0, tg_yyyzz_xxxxzzzz_0, \
                                         tg_yyyzz_xxxyyyyy_0, tg_yyyzz_xxxyyyyz_0, tg_yyyzz_xxxyyyzz_0, tg_yyyzz_xxxyyzzz_0, \
                                         tg_yyyzz_xxxyzzzz_0, tg_yyyzz_xxxzzzzz_0, tg_yyyzz_xxyyyyyy_0, tg_yyyzz_xxyyyyyz_0, \
                                         tg_yyyzz_xxyyyyzz_0, tg_yyyzz_xxyyyzzz_0, tg_yyyzz_xxyyzzzz_0, tg_yyyzz_xxyzzzzz_0, \
                                         tg_yyyzz_xxzzzzzz_0, tg_yyyzz_xyyyyyyy_0, tg_yyyzz_xyyyyyyz_0, tg_yyyzz_xyyyyyzz_0, \
                                         tg_yyyzz_xyyyyzzz_0, tg_yyyzz_xyyyzzzz_0, tg_yyyzz_xyyzzzzz_0, tg_yyyzz_xyzzzzzz_0, \
                                         tg_yyyzz_xzzzzzzz_0, tg_yyyzz_yyyyyyyy_0, tg_yyyzz_yyyyyyyz_0, tg_yyyzz_yyyyyyzz_0, \
                                         tg_yyyzz_yyyyyzzz_0, tg_yyyzz_yyyyzzzz_0, tg_yyyzz_yyyzzzzz_0, tg_yyyzz_yyzzzzzz_0, \
                                         tg_yyyzz_yzzzzzzz_0, tg_yyyzz_zzzzzzzz_0, tg_yyz_yyyyyyyz_0, tg_yyz_yyyyyyyz_1, \
                                         tg_yyz_yyyyyyzz_0, tg_yyz_yyyyyyzz_1, tg_yyz_yyyyyzzz_0, tg_yyz_yyyyyzzz_1, \
                                         tg_yyz_yyyyzzzz_0, tg_yyz_yyyyzzzz_1, tg_yyz_yyyzzzzz_0, tg_yyz_yyyzzzzz_1, \
                                         tg_yyz_yyzzzzzz_0, tg_yyz_yyzzzzzz_1, tg_yyz_yzzzzzzz_0, tg_yyz_yzzzzzzz_1, \
                                         tg_yyz_zzzzzzzz_0, tg_yyz_zzzzzzzz_1, tg_yyzz_xxxxxxx_1, tg_yyzz_xxxxxxxx_0, \
                                         tg_yyzz_xxxxxxxx_1, tg_yyzz_xxxxxxxy_0, tg_yyzz_xxxxxxxy_1, tg_yyzz_xxxxxxxz_0, \
                                         tg_yyzz_xxxxxxxz_1, tg_yyzz_xxxxxxy_1, tg_yyzz_xxxxxxyy_0, tg_yyzz_xxxxxxyy_1, \
                                         tg_yyzz_xxxxxxyz_0, tg_yyzz_xxxxxxyz_1, tg_yyzz_xxxxxxz_1, tg_yyzz_xxxxxxzz_0, \
                                         tg_yyzz_xxxxxxzz_1, tg_yyzz_xxxxxyy_1, tg_yyzz_xxxxxyyy_0, tg_yyzz_xxxxxyyy_1, \
                                         tg_yyzz_xxxxxyyz_0, tg_yyzz_xxxxxyyz_1, tg_yyzz_xxxxxyz_1, tg_yyzz_xxxxxyzz_0, \
                                         tg_yyzz_xxxxxyzz_1, tg_yyzz_xxxxxzz_1, tg_yyzz_xxxxxzzz_0, tg_yyzz_xxxxxzzz_1, \
                                         tg_yyzz_xxxxyyy_1, tg_yyzz_xxxxyyyy_0, tg_yyzz_xxxxyyyy_1, tg_yyzz_xxxxyyyz_0, \
                                         tg_yyzz_xxxxyyyz_1, tg_yyzz_xxxxyyz_1, tg_yyzz_xxxxyyzz_0, tg_yyzz_xxxxyyzz_1, \
                                         tg_yyzz_xxxxyzz_1, tg_yyzz_xxxxyzzz_0, tg_yyzz_xxxxyzzz_1, tg_yyzz_xxxxzzz_1, \
                                         tg_yyzz_xxxxzzzz_0, tg_yyzz_xxxxzzzz_1, tg_yyzz_xxxyyyy_1, tg_yyzz_xxxyyyyy_0, \
                                         tg_yyzz_xxxyyyyy_1, tg_yyzz_xxxyyyyz_0, tg_yyzz_xxxyyyyz_1, tg_yyzz_xxxyyyz_1, \
                                         tg_yyzz_xxxyyyzz_0, tg_yyzz_xxxyyyzz_1, tg_yyzz_xxxyyzz_1, tg_yyzz_xxxyyzzz_0, \
                                         tg_yyzz_xxxyyzzz_1, tg_yyzz_xxxyzzz_1, tg_yyzz_xxxyzzzz_0, tg_yyzz_xxxyzzzz_1, \
                                         tg_yyzz_xxxzzzz_1, tg_yyzz_xxxzzzzz_0, tg_yyzz_xxxzzzzz_1, tg_yyzz_xxyyyyy_1, \
                                         tg_yyzz_xxyyyyyy_0, tg_yyzz_xxyyyyyy_1, tg_yyzz_xxyyyyyz_0, tg_yyzz_xxyyyyyz_1, \
                                         tg_yyzz_xxyyyyz_1, tg_yyzz_xxyyyyzz_0, tg_yyzz_xxyyyyzz_1, tg_yyzz_xxyyyzz_1, \
                                         tg_yyzz_xxyyyzzz_0, tg_yyzz_xxyyyzzz_1, tg_yyzz_xxyyzzz_1, tg_yyzz_xxyyzzzz_0, \
                                         tg_yyzz_xxyyzzzz_1, tg_yyzz_xxyzzzz_1, tg_yyzz_xxyzzzzz_0, tg_yyzz_xxyzzzzz_1, \
                                         tg_yyzz_xxzzzzz_1, tg_yyzz_xxzzzzzz_0, tg_yyzz_xxzzzzzz_1, tg_yyzz_xyyyyyy_1, \
                                         tg_yyzz_xyyyyyyy_0, tg_yyzz_xyyyyyyy_1, tg_yyzz_xyyyyyyz_0, tg_yyzz_xyyyyyyz_1, \
                                         tg_yyzz_xyyyyyz_1, tg_yyzz_xyyyyyzz_0, tg_yyzz_xyyyyyzz_1, tg_yyzz_xyyyyzz_1, \
                                         tg_yyzz_xyyyyzzz_0, tg_yyzz_xyyyyzzz_1, tg_yyzz_xyyyzzz_1, tg_yyzz_xyyyzzzz_0, \
                                         tg_yyzz_xyyyzzzz_1, tg_yyzz_xyyzzzz_1, tg_yyzz_xyyzzzzz_0, tg_yyzz_xyyzzzzz_1, \
                                         tg_yyzz_xyzzzzz_1, tg_yyzz_xyzzzzzz_0, tg_yyzz_xyzzzzzz_1, tg_yyzz_xzzzzzz_1, \
                                         tg_yyzz_xzzzzzzz_0, tg_yyzz_xzzzzzzz_1, tg_yyzz_yyyyyyy_1, tg_yyzz_yyyyyyyy_0, \
                                         tg_yyzz_yyyyyyyy_1, tg_yyzz_yyyyyyyz_0, tg_yyzz_yyyyyyyz_1, tg_yyzz_yyyyyyz_1, \
                                         tg_yyzz_yyyyyyzz_0, tg_yyzz_yyyyyyzz_1, tg_yyzz_yyyyyzz_1, tg_yyzz_yyyyyzzz_0, \
                                         tg_yyzz_yyyyyzzz_1, tg_yyzz_yyyyzzz_1, tg_yyzz_yyyyzzzz_0, tg_yyzz_yyyyzzzz_1, \
                                         tg_yyzz_yyyzzzz_1, tg_yyzz_yyyzzzzz_0, tg_yyzz_yyyzzzzz_1, tg_yyzz_yyzzzzz_1, \
                                         tg_yyzz_yyzzzzzz_0, tg_yyzz_yyzzzzzz_1, tg_yyzz_yzzzzzz_1, tg_yyzz_yzzzzzzz_0, \
                                         tg_yyzz_yzzzzzzz_1, tg_yyzz_zzzzzzz_1, tg_yyzz_zzzzzzzz_0, tg_yyzz_zzzzzzzz_1, \
                                         tg_yyzzz_xxxxxxxx_0, tg_yyzzz_xxxxxxxy_0, tg_yyzzz_xxxxxxxz_0, tg_yyzzz_xxxxxxyy_0, \
                                         tg_yyzzz_xxxxxxyz_0, tg_yyzzz_xxxxxxzz_0, tg_yyzzz_xxxxxyyy_0, tg_yyzzz_xxxxxyyz_0, \
                                         tg_yyzzz_xxxxxyzz_0, tg_yyzzz_xxxxxzzz_0, tg_yyzzz_xxxxyyyy_0, tg_yyzzz_xxxxyyyz_0, \
                                         tg_yyzzz_xxxxyyzz_0, tg_yyzzz_xxxxyzzz_0, tg_yyzzz_xxxxzzzz_0, tg_yyzzz_xxxyyyyy_0, \
                                         tg_yyzzz_xxxyyyyz_0, tg_yyzzz_xxxyyyzz_0, tg_yyzzz_xxxyyzzz_0, tg_yyzzz_xxxyzzzz_0, \
                                         tg_yyzzz_xxxzzzzz_0, tg_yyzzz_xxyyyyyy_0, tg_yyzzz_xxyyyyyz_0, tg_yyzzz_xxyyyyzz_0, \
                                         tg_yyzzz_xxyyyzzz_0, tg_yyzzz_xxyyzzzz_0, tg_yyzzz_xxyzzzzz_0, tg_yyzzz_xxzzzzzz_0, \
                                         tg_yyzzz_xyyyyyyy_0, tg_yyzzz_xyyyyyyz_0, tg_yyzzz_xyyyyyzz_0, tg_yyzzz_xyyyyzzz_0, \
                                         tg_yyzzz_xyyyzzzz_0, tg_yyzzz_xyyzzzzz_0, tg_yyzzz_xyzzzzzz_0, tg_yyzzz_xzzzzzzz_0, \
                                         tg_yyzzz_yyyyyyyy_0, tg_yyzzz_yyyyyyyz_0, tg_yyzzz_yyyyyyzz_0, tg_yyzzz_yyyyyzzz_0, \
                                         tg_yyzzz_yyyyzzzz_0, tg_yzz_xxxxxxxx_0, tg_yzz_xxxxxxxx_1, tg_yzz_xxxxxxxy_0, \
                                         tg_yzz_xxxxxxxy_1, tg_yzz_xxxxxxxz_0, tg_yzz_xxxxxxxz_1, tg_yzz_xxxxxxyy_0, \
                                         tg_yzz_xxxxxxyy_1, tg_yzz_xxxxxxyz_0, tg_yzz_xxxxxxyz_1, tg_yzz_xxxxxxzz_0, \
                                         tg_yzz_xxxxxxzz_1, tg_yzz_xxxxxyyy_0, tg_yzz_xxxxxyyy_1, tg_yzz_xxxxxyyz_0, \
                                         tg_yzz_xxxxxyyz_1, tg_yzz_xxxxxyzz_0, tg_yzz_xxxxxyzz_1, tg_yzz_xxxxxzzz_0, \
                                         tg_yzz_xxxxxzzz_1, tg_yzz_xxxxyyyy_0, tg_yzz_xxxxyyyy_1, tg_yzz_xxxxyyyz_0, \
                                         tg_yzz_xxxxyyyz_1, tg_yzz_xxxxyyzz_0, tg_yzz_xxxxyyzz_1, tg_yzz_xxxxyzzz_0, \
                                         tg_yzz_xxxxyzzz_1, tg_yzz_xxxxzzzz_0, tg_yzz_xxxxzzzz_1, tg_yzz_xxxyyyyy_0, \
                                         tg_yzz_xxxyyyyy_1, tg_yzz_xxxyyyyz_0, tg_yzz_xxxyyyyz_1, tg_yzz_xxxyyyzz_0, \
                                         tg_yzz_xxxyyyzz_1, tg_yzz_xxxyyzzz_0, tg_yzz_xxxyyzzz_1, tg_yzz_xxxyzzzz_0, \
                                         tg_yzz_xxxyzzzz_1, tg_yzz_xxxzzzzz_0, tg_yzz_xxxzzzzz_1, tg_yzz_xxyyyyyy_0, \
                                         tg_yzz_xxyyyyyy_1, tg_yzz_xxyyyyyz_0, tg_yzz_xxyyyyyz_1, tg_yzz_xxyyyyzz_0, \
                                         tg_yzz_xxyyyyzz_1, tg_yzz_xxyyyzzz_0, tg_yzz_xxyyyzzz_1, tg_yzz_xxyyzzzz_0, \
                                         tg_yzz_xxyyzzzz_1, tg_yzz_xxyzzzzz_0, tg_yzz_xxyzzzzz_1, tg_yzz_xxzzzzzz_0, \
                                         tg_yzz_xxzzzzzz_1, tg_yzz_xyyyyyyy_0, tg_yzz_xyyyyyyy_1, tg_yzz_xyyyyyyz_0, \
                                         tg_yzz_xyyyyyyz_1, tg_yzz_xyyyyyzz_0, tg_yzz_xyyyyyzz_1, tg_yzz_xyyyyzzz_0, \
                                         tg_yzz_xyyyyzzz_1, tg_yzz_xyyyzzzz_0, tg_yzz_xyyyzzzz_1, tg_yzz_xyyzzzzz_0, \
                                         tg_yzz_xyyzzzzz_1, tg_yzz_xyzzzzzz_0, tg_yzz_xyzzzzzz_1, tg_yzz_xzzzzzzz_0, \
                                         tg_yzz_xzzzzzzz_1, tg_yzz_yyyyyyyy_0, tg_yzz_yyyyyyyy_1, tg_yzz_yyyyyyyz_0, \
                                         tg_yzz_yyyyyyyz_1, tg_yzz_yyyyyyzz_0, tg_yzz_yyyyyyzz_1, tg_yzz_yyyyyzzz_0, \
                                         tg_yzz_yyyyyzzz_1, tg_yzz_yyyyzzzz_0, tg_yzz_yyyyzzzz_1, tg_yzz_yyyzzzzz_0, \
                                         tg_yzz_yyyzzzzz_1, tg_yzz_yyzzzzzz_0, tg_yzz_yyzzzzzz_1, tg_yzz_yzzzzzzz_0, \
                                         tg_yzz_yzzzzzzz_1, tg_yzz_zzzzzzzz_0, tg_yzz_zzzzzzzz_1, tg_yzzz_xxxxxxx_1, \
                                         tg_yzzz_xxxxxxxx_0, tg_yzzz_xxxxxxxx_1, tg_yzzz_xxxxxxxy_0, tg_yzzz_xxxxxxxy_1, \
                                         tg_yzzz_xxxxxxxz_0, tg_yzzz_xxxxxxxz_1, tg_yzzz_xxxxxxy_1, tg_yzzz_xxxxxxyy_0, \
                                         tg_yzzz_xxxxxxyy_1, tg_yzzz_xxxxxxyz_0, tg_yzzz_xxxxxxyz_1, tg_yzzz_xxxxxxz_1, \
                                         tg_yzzz_xxxxxxzz_0, tg_yzzz_xxxxxxzz_1, tg_yzzz_xxxxxyy_1, tg_yzzz_xxxxxyyy_0, \
                                         tg_yzzz_xxxxxyyy_1, tg_yzzz_xxxxxyyz_0, tg_yzzz_xxxxxyyz_1, tg_yzzz_xxxxxyz_1, \
                                         tg_yzzz_xxxxxyzz_0, tg_yzzz_xxxxxyzz_1, tg_yzzz_xxxxxzz_1, tg_yzzz_xxxxxzzz_0, \
                                         tg_yzzz_xxxxxzzz_1, tg_yzzz_xxxxyyy_1, tg_yzzz_xxxxyyyy_0, tg_yzzz_xxxxyyyy_1, \
                                         tg_yzzz_xxxxyyyz_0, tg_yzzz_xxxxyyyz_1, tg_yzzz_xxxxyyz_1, tg_yzzz_xxxxyyzz_0, \
                                         tg_yzzz_xxxxyyzz_1, tg_yzzz_xxxxyzz_1, tg_yzzz_xxxxyzzz_0, tg_yzzz_xxxxyzzz_1, \
                                         tg_yzzz_xxxxzzz_1, tg_yzzz_xxxxzzzz_0, tg_yzzz_xxxxzzzz_1, tg_yzzz_xxxyyyy_1, \
                                         tg_yzzz_xxxyyyyy_0, tg_yzzz_xxxyyyyy_1, tg_yzzz_xxxyyyyz_0, tg_yzzz_xxxyyyyz_1, \
                                         tg_yzzz_xxxyyyz_1, tg_yzzz_xxxyyyzz_0, tg_yzzz_xxxyyyzz_1, tg_yzzz_xxxyyzz_1, \
                                         tg_yzzz_xxxyyzzz_0, tg_yzzz_xxxyyzzz_1, tg_yzzz_xxxyzzz_1, tg_yzzz_xxxyzzzz_0, \
                                         tg_yzzz_xxxyzzzz_1, tg_yzzz_xxxzzzz_1, tg_yzzz_xxxzzzzz_0, tg_yzzz_xxxzzzzz_1, \
                                         tg_yzzz_xxyyyyy_1, tg_yzzz_xxyyyyyy_0, tg_yzzz_xxyyyyyy_1, tg_yzzz_xxyyyyyz_0, \
                                         tg_yzzz_xxyyyyyz_1, tg_yzzz_xxyyyyz_1, tg_yzzz_xxyyyyzz_0, tg_yzzz_xxyyyyzz_1, \
                                         tg_yzzz_xxyyyzz_1, tg_yzzz_xxyyyzzz_0, tg_yzzz_xxyyyzzz_1, tg_yzzz_xxyyzzz_1, \
                                         tg_yzzz_xxyyzzzz_0, tg_yzzz_xxyyzzzz_1, tg_yzzz_xxyzzzz_1, tg_yzzz_xxyzzzzz_0, \
                                         tg_yzzz_xxyzzzzz_1, tg_yzzz_xxzzzzz_1, tg_yzzz_xxzzzzzz_0, tg_yzzz_xxzzzzzz_1, \
                                         tg_yzzz_xyyyyyy_1, tg_yzzz_xyyyyyyy_0, tg_yzzz_xyyyyyyy_1, tg_yzzz_xyyyyyyz_0, \
                                         tg_yzzz_xyyyyyyz_1, tg_yzzz_xyyyyyz_1, tg_yzzz_xyyyyyzz_0, tg_yzzz_xyyyyyzz_1, \
                                         tg_yzzz_xyyyyzz_1, tg_yzzz_xyyyyzzz_0, tg_yzzz_xyyyyzzz_1, tg_yzzz_xyyyzzz_1, \
                                         tg_yzzz_xyyyzzzz_0, tg_yzzz_xyyyzzzz_1, tg_yzzz_xyyzzzz_1, tg_yzzz_xyyzzzzz_0, \
                                         tg_yzzz_xyyzzzzz_1, tg_yzzz_xyzzzzz_1, tg_yzzz_xyzzzzzz_0, tg_yzzz_xyzzzzzz_1, \
                                         tg_yzzz_xzzzzzz_1, tg_yzzz_xzzzzzzz_0, tg_yzzz_xzzzzzzz_1, tg_yzzz_yyyyyyy_1, \
                                         tg_yzzz_yyyyyyyy_0, tg_yzzz_yyyyyyyy_1, tg_yzzz_yyyyyyyz_0, tg_yzzz_yyyyyyyz_1, \
                                         tg_yzzz_yyyyyyz_1, tg_yzzz_yyyyyyzz_0, tg_yzzz_yyyyyyzz_1, tg_yzzz_yyyyyzz_1, \
                                         tg_yzzz_yyyyyzzz_0, tg_yzzz_yyyyyzzz_1, tg_yzzz_yyyyzzz_1, tg_yzzz_yyyyzzzz_0, \
                                         tg_yzzz_yyyyzzzz_1, tg_yzzz_yyyzzzz_1, tg_zzz_xxxxxxxx_0, tg_zzz_xxxxxxxx_1, \
                                         tg_zzz_xxxxxxxy_0, tg_zzz_xxxxxxxy_1, tg_zzz_xxxxxxxz_0, tg_zzz_xxxxxxxz_1, \
                                         tg_zzz_xxxxxxyy_0, tg_zzz_xxxxxxyy_1, tg_zzz_xxxxxxyz_0, tg_zzz_xxxxxxyz_1, \
                                         tg_zzz_xxxxxxzz_0, tg_zzz_xxxxxxzz_1, tg_zzz_xxxxxyyy_0, tg_zzz_xxxxxyyy_1, \
                                         tg_zzz_xxxxxyyz_0, tg_zzz_xxxxxyyz_1, tg_zzz_xxxxxyzz_0, tg_zzz_xxxxxyzz_1, \
                                         tg_zzz_xxxxxzzz_0, tg_zzz_xxxxxzzz_1, tg_zzz_xxxxyyyy_0, tg_zzz_xxxxyyyy_1, \
                                         tg_zzz_xxxxyyyz_0, tg_zzz_xxxxyyyz_1, tg_zzz_xxxxyyzz_0, tg_zzz_xxxxyyzz_1, \
                                         tg_zzz_xxxxyzzz_0, tg_zzz_xxxxyzzz_1, tg_zzz_xxxxzzzz_0, tg_zzz_xxxxzzzz_1, \
                                         tg_zzz_xxxyyyyy_0, tg_zzz_xxxyyyyy_1, tg_zzz_xxxyyyyz_0, tg_zzz_xxxyyyyz_1, \
                                         tg_zzz_xxxyyyzz_0, tg_zzz_xxxyyyzz_1, tg_zzz_xxxyyzzz_0, tg_zzz_xxxyyzzz_1, \
                                         tg_zzz_xxxyzzzz_0, tg_zzz_xxxyzzzz_1, tg_zzz_xxxzzzzz_0, tg_zzz_xxxzzzzz_1, \
                                         tg_zzz_xxyyyyyy_0, tg_zzz_xxyyyyyy_1, tg_zzz_xxyyyyyz_0, tg_zzz_xxyyyyyz_1, \
                                         tg_zzz_xxyyyyzz_0, tg_zzz_xxyyyyzz_1, tg_zzz_xxyyyzzz_0, tg_zzz_xxyyyzzz_1, \
                                         tg_zzz_xxyyzzzz_0, tg_zzz_xxyyzzzz_1, tg_zzz_xxyzzzzz_0, tg_zzz_xxyzzzzz_1, \
                                         tg_zzz_xxzzzzzz_0, tg_zzz_xxzzzzzz_1, tg_zzz_xyyyyyyy_0, tg_zzz_xyyyyyyy_1, \
                                         tg_zzz_xyyyyyyz_0, tg_zzz_xyyyyyyz_1, tg_zzz_xyyyyyzz_0, tg_zzz_xyyyyyzz_1, \
                                         tg_zzz_xyyyyzzz_0, tg_zzz_xyyyyzzz_1, tg_zzz_xyyyzzzz_0, tg_zzz_xyyyzzzz_1, \
                                         tg_zzz_xyyzzzzz_0, tg_zzz_xyyzzzzz_1, tg_zzz_xyzzzzzz_0, tg_zzz_xyzzzzzz_1, \
                                         tg_zzz_xzzzzzzz_0, tg_zzz_xzzzzzzz_1, tg_zzz_yyyyyyyy_0, tg_zzz_yyyyyyyy_1, \
                                         tg_zzz_yyyyyyyz_0, tg_zzz_yyyyyyyz_1, tg_zzz_yyyyyyzz_0, tg_zzz_yyyyyyzz_1, \
                                         tg_zzz_yyyyyzzz_0, tg_zzz_yyyyyzzz_1, tg_zzz_yyyyzzzz_0, tg_zzz_yyyyzzzz_1, wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_y[j]; 

                    tg_yyyyz_yyyyyyyz_0[j] = pb_y * tg_yyyz_yyyyyyyz_0[j] + fr * tg_yyyz_yyyyyyyz_1[j] + 1.5 * fl1_fx * (tg_yyz_yyyyyyyz_0[j] - tg_yyz_yyyyyyyz_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_yyyz_yyyyyyz_1[j];

                    tg_yyyyz_yyyyyyzz_0[j] = pb_y * tg_yyyz_yyyyyyzz_0[j] + fr * tg_yyyz_yyyyyyzz_1[j] + 1.5 * fl1_fx * (tg_yyz_yyyyyyzz_0[j] - tg_yyz_yyyyyyzz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yyyz_yyyyyzz_1[j];

                    tg_yyyyz_yyyyyzzz_0[j] = pb_y * tg_yyyz_yyyyyzzz_0[j] + fr * tg_yyyz_yyyyyzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_yyyyyzzz_0[j] - tg_yyz_yyyyyzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyz_yyyyzzz_1[j];

                    tg_yyyyz_yyyyzzzz_0[j] = pb_y * tg_yyyz_yyyyzzzz_0[j] + fr * tg_yyyz_yyyyzzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_yyyyzzzz_0[j] - tg_yyz_yyyyzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyz_yyyzzzz_1[j];

                    tg_yyyyz_yyyzzzzz_0[j] = pb_y * tg_yyyz_yyyzzzzz_0[j] + fr * tg_yyyz_yyyzzzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_yyyzzzzz_0[j] - tg_yyz_yyyzzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyz_yyzzzzz_1[j];

                    tg_yyyyz_yyzzzzzz_0[j] = pb_y * tg_yyyz_yyzzzzzz_0[j] + fr * tg_yyyz_yyzzzzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_yyzzzzzz_0[j] - tg_yyz_yyzzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyz_yzzzzzz_1[j];

                    tg_yyyyz_yzzzzzzz_0[j] = pb_y * tg_yyyz_yzzzzzzz_0[j] + fr * tg_yyyz_yzzzzzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_yzzzzzzz_0[j] - tg_yyz_yzzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyz_zzzzzzz_1[j];

                    tg_yyyyz_zzzzzzzz_0[j] = pb_y * tg_yyyz_zzzzzzzz_0[j] + fr * tg_yyyz_zzzzzzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_zzzzzzzz_0[j] - tg_yyz_zzzzzzzz_1[j] * fl1_fza);

                    tg_yyyzz_xxxxxxxx_0[j] = pb_y * tg_yyzz_xxxxxxxx_0[j] + fr * tg_yyzz_xxxxxxxx_1[j] + fl1_fx * (tg_yzz_xxxxxxxx_0[j] - tg_yzz_xxxxxxxx_1[j] * fl1_fza);

                    tg_yyyzz_xxxxxxxy_0[j] = pb_y * tg_yyzz_xxxxxxxy_0[j] + fr * tg_yyzz_xxxxxxxy_1[j] + fl1_fx * (tg_yzz_xxxxxxxy_0[j] - tg_yzz_xxxxxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzz_xxxxxxx_1[j];

                    tg_yyyzz_xxxxxxxz_0[j] = pb_y * tg_yyzz_xxxxxxxz_0[j] + fr * tg_yyzz_xxxxxxxz_1[j] + fl1_fx * (tg_yzz_xxxxxxxz_0[j] - tg_yzz_xxxxxxxz_1[j] * fl1_fza);

                    tg_yyyzz_xxxxxxyy_0[j] = pb_y * tg_yyzz_xxxxxxyy_0[j] + fr * tg_yyzz_xxxxxxyy_1[j] + fl1_fx * (tg_yzz_xxxxxxyy_0[j] - tg_yzz_xxxxxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyzz_xxxxxxy_1[j];

                    tg_yyyzz_xxxxxxyz_0[j] = pb_y * tg_yyzz_xxxxxxyz_0[j] + fr * tg_yyzz_xxxxxxyz_1[j] + fl1_fx * (tg_yzz_xxxxxxyz_0[j] - tg_yzz_xxxxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzz_xxxxxxz_1[j];

                    tg_yyyzz_xxxxxxzz_0[j] = pb_y * tg_yyzz_xxxxxxzz_0[j] + fr * tg_yyzz_xxxxxxzz_1[j] + fl1_fx * (tg_yzz_xxxxxxzz_0[j] - tg_yzz_xxxxxxzz_1[j] * fl1_fza);

                    tg_yyyzz_xxxxxyyy_0[j] = pb_y * tg_yyzz_xxxxxyyy_0[j] + fr * tg_yyzz_xxxxxyyy_1[j] + fl1_fx * (tg_yzz_xxxxxyyy_0[j] - tg_yzz_xxxxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzz_xxxxxyy_1[j];

                    tg_yyyzz_xxxxxyyz_0[j] = pb_y * tg_yyzz_xxxxxyyz_0[j] + fr * tg_yyzz_xxxxxyyz_1[j] + fl1_fx * (tg_yzz_xxxxxyyz_0[j] - tg_yzz_xxxxxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyzz_xxxxxyz_1[j];

                    tg_yyyzz_xxxxxyzz_0[j] = pb_y * tg_yyzz_xxxxxyzz_0[j] + fr * tg_yyzz_xxxxxyzz_1[j] + fl1_fx * (tg_yzz_xxxxxyzz_0[j] - tg_yzz_xxxxxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzz_xxxxxzz_1[j];

                    tg_yyyzz_xxxxxzzz_0[j] = pb_y * tg_yyzz_xxxxxzzz_0[j] + fr * tg_yyzz_xxxxxzzz_1[j] + fl1_fx * (tg_yzz_xxxxxzzz_0[j] - tg_yzz_xxxxxzzz_1[j] * fl1_fza);

                    tg_yyyzz_xxxxyyyy_0[j] = pb_y * tg_yyzz_xxxxyyyy_0[j] + fr * tg_yyzz_xxxxyyyy_1[j] + fl1_fx * (tg_yzz_xxxxyyyy_0[j] - tg_yzz_xxxxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyzz_xxxxyyy_1[j];

                    tg_yyyzz_xxxxyyyz_0[j] = pb_y * tg_yyzz_xxxxyyyz_0[j] + fr * tg_yyzz_xxxxyyyz_1[j] + fl1_fx * (tg_yzz_xxxxyyyz_0[j] - tg_yzz_xxxxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzz_xxxxyyz_1[j];

                    tg_yyyzz_xxxxyyzz_0[j] = pb_y * tg_yyzz_xxxxyyzz_0[j] + fr * tg_yyzz_xxxxyyzz_1[j] + fl1_fx * (tg_yzz_xxxxyyzz_0[j] - tg_yzz_xxxxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyzz_xxxxyzz_1[j];

                    tg_yyyzz_xxxxyzzz_0[j] = pb_y * tg_yyzz_xxxxyzzz_0[j] + fr * tg_yyzz_xxxxyzzz_1[j] + fl1_fx * (tg_yzz_xxxxyzzz_0[j] - tg_yzz_xxxxyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzz_xxxxzzz_1[j];

                    tg_yyyzz_xxxxzzzz_0[j] = pb_y * tg_yyzz_xxxxzzzz_0[j] + fr * tg_yyzz_xxxxzzzz_1[j] + fl1_fx * (tg_yzz_xxxxzzzz_0[j] - tg_yzz_xxxxzzzz_1[j] * fl1_fza);

                    tg_yyyzz_xxxyyyyy_0[j] = pb_y * tg_yyzz_xxxyyyyy_0[j] + fr * tg_yyzz_xxxyyyyy_1[j] + fl1_fx * (tg_yzz_xxxyyyyy_0[j] - tg_yzz_xxxyyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyzz_xxxyyyy_1[j];

                    tg_yyyzz_xxxyyyyz_0[j] = pb_y * tg_yyzz_xxxyyyyz_0[j] + fr * tg_yyzz_xxxyyyyz_1[j] + fl1_fx * (tg_yzz_xxxyyyyz_0[j] - tg_yzz_xxxyyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyzz_xxxyyyz_1[j];

                    tg_yyyzz_xxxyyyzz_0[j] = pb_y * tg_yyzz_xxxyyyzz_0[j] + fr * tg_yyzz_xxxyyyzz_1[j] + fl1_fx * (tg_yzz_xxxyyyzz_0[j] - tg_yzz_xxxyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzz_xxxyyzz_1[j];

                    tg_yyyzz_xxxyyzzz_0[j] = pb_y * tg_yyzz_xxxyyzzz_0[j] + fr * tg_yyzz_xxxyyzzz_1[j] + fl1_fx * (tg_yzz_xxxyyzzz_0[j] - tg_yzz_xxxyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyzz_xxxyzzz_1[j];

                    tg_yyyzz_xxxyzzzz_0[j] = pb_y * tg_yyzz_xxxyzzzz_0[j] + fr * tg_yyzz_xxxyzzzz_1[j] + fl1_fx * (tg_yzz_xxxyzzzz_0[j] - tg_yzz_xxxyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzz_xxxzzzz_1[j];

                    tg_yyyzz_xxxzzzzz_0[j] = pb_y * tg_yyzz_xxxzzzzz_0[j] + fr * tg_yyzz_xxxzzzzz_1[j] + fl1_fx * (tg_yzz_xxxzzzzz_0[j] - tg_yzz_xxxzzzzz_1[j] * fl1_fza);

                    tg_yyyzz_xxyyyyyy_0[j] = pb_y * tg_yyzz_xxyyyyyy_0[j] + fr * tg_yyzz_xxyyyyyy_1[j] + fl1_fx * (tg_yzz_xxyyyyyy_0[j] - tg_yzz_xxyyyyyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yyzz_xxyyyyy_1[j];

                    tg_yyyzz_xxyyyyyz_0[j] = pb_y * tg_yyzz_xxyyyyyz_0[j] + fr * tg_yyzz_xxyyyyyz_1[j] + fl1_fx * (tg_yzz_xxyyyyyz_0[j] - tg_yzz_xxyyyyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyzz_xxyyyyz_1[j];

                    tg_yyyzz_xxyyyyzz_0[j] = pb_y * tg_yyzz_xxyyyyzz_0[j] + fr * tg_yyzz_xxyyyyzz_1[j] + fl1_fx * (tg_yzz_xxyyyyzz_0[j] - tg_yzz_xxyyyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyzz_xxyyyzz_1[j];

                    tg_yyyzz_xxyyyzzz_0[j] = pb_y * tg_yyzz_xxyyyzzz_0[j] + fr * tg_yyzz_xxyyyzzz_1[j] + fl1_fx * (tg_yzz_xxyyyzzz_0[j] - tg_yzz_xxyyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzz_xxyyzzz_1[j];

                    tg_yyyzz_xxyyzzzz_0[j] = pb_y * tg_yyzz_xxyyzzzz_0[j] + fr * tg_yyzz_xxyyzzzz_1[j] + fl1_fx * (tg_yzz_xxyyzzzz_0[j] - tg_yzz_xxyyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyzz_xxyzzzz_1[j];

                    tg_yyyzz_xxyzzzzz_0[j] = pb_y * tg_yyzz_xxyzzzzz_0[j] + fr * tg_yyzz_xxyzzzzz_1[j] + fl1_fx * (tg_yzz_xxyzzzzz_0[j] - tg_yzz_xxyzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzz_xxzzzzz_1[j];

                    tg_yyyzz_xxzzzzzz_0[j] = pb_y * tg_yyzz_xxzzzzzz_0[j] + fr * tg_yyzz_xxzzzzzz_1[j] + fl1_fx * (tg_yzz_xxzzzzzz_0[j] - tg_yzz_xxzzzzzz_1[j] * fl1_fza);

                    tg_yyyzz_xyyyyyyy_0[j] = pb_y * tg_yyzz_xyyyyyyy_0[j] + fr * tg_yyzz_xyyyyyyy_1[j] + fl1_fx * (tg_yzz_xyyyyyyy_0[j] - tg_yzz_xyyyyyyy_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_yyzz_xyyyyyy_1[j];

                    tg_yyyzz_xyyyyyyz_0[j] = pb_y * tg_yyzz_xyyyyyyz_0[j] + fr * tg_yyzz_xyyyyyyz_1[j] + fl1_fx * (tg_yzz_xyyyyyyz_0[j] - tg_yzz_xyyyyyyz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yyzz_xyyyyyz_1[j];

                    tg_yyyzz_xyyyyyzz_0[j] = pb_y * tg_yyzz_xyyyyyzz_0[j] + fr * tg_yyzz_xyyyyyzz_1[j] + fl1_fx * (tg_yzz_xyyyyyzz_0[j] - tg_yzz_xyyyyyzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyzz_xyyyyzz_1[j];

                    tg_yyyzz_xyyyyzzz_0[j] = pb_y * tg_yyzz_xyyyyzzz_0[j] + fr * tg_yyzz_xyyyyzzz_1[j] + fl1_fx * (tg_yzz_xyyyyzzz_0[j] - tg_yzz_xyyyyzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyzz_xyyyzzz_1[j];

                    tg_yyyzz_xyyyzzzz_0[j] = pb_y * tg_yyzz_xyyyzzzz_0[j] + fr * tg_yyzz_xyyyzzzz_1[j] + fl1_fx * (tg_yzz_xyyyzzzz_0[j] - tg_yzz_xyyyzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzz_xyyzzzz_1[j];

                    tg_yyyzz_xyyzzzzz_0[j] = pb_y * tg_yyzz_xyyzzzzz_0[j] + fr * tg_yyzz_xyyzzzzz_1[j] + fl1_fx * (tg_yzz_xyyzzzzz_0[j] - tg_yzz_xyyzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyzz_xyzzzzz_1[j];

                    tg_yyyzz_xyzzzzzz_0[j] = pb_y * tg_yyzz_xyzzzzzz_0[j] + fr * tg_yyzz_xyzzzzzz_1[j] + fl1_fx * (tg_yzz_xyzzzzzz_0[j] - tg_yzz_xyzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzz_xzzzzzz_1[j];

                    tg_yyyzz_xzzzzzzz_0[j] = pb_y * tg_yyzz_xzzzzzzz_0[j] + fr * tg_yyzz_xzzzzzzz_1[j] + fl1_fx * (tg_yzz_xzzzzzzz_0[j] - tg_yzz_xzzzzzzz_1[j] * fl1_fza);

                    tg_yyyzz_yyyyyyyy_0[j] = pb_y * tg_yyzz_yyyyyyyy_0[j] + fr * tg_yyzz_yyyyyyyy_1[j] + fl1_fx * (tg_yzz_yyyyyyyy_0[j] - tg_yzz_yyyyyyyy_1[j] * fl1_fza) + 4.0 * fl1_fxn * tg_yyzz_yyyyyyy_1[j];

                    tg_yyyzz_yyyyyyyz_0[j] = pb_y * tg_yyzz_yyyyyyyz_0[j] + fr * tg_yyzz_yyyyyyyz_1[j] + fl1_fx * (tg_yzz_yyyyyyyz_0[j] - tg_yzz_yyyyyyyz_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_yyzz_yyyyyyz_1[j];

                    tg_yyyzz_yyyyyyzz_0[j] = pb_y * tg_yyzz_yyyyyyzz_0[j] + fr * tg_yyzz_yyyyyyzz_1[j] + fl1_fx * (tg_yzz_yyyyyyzz_0[j] - tg_yzz_yyyyyyzz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yyzz_yyyyyzz_1[j];

                    tg_yyyzz_yyyyyzzz_0[j] = pb_y * tg_yyzz_yyyyyzzz_0[j] + fr * tg_yyzz_yyyyyzzz_1[j] + fl1_fx * (tg_yzz_yyyyyzzz_0[j] - tg_yzz_yyyyyzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyzz_yyyyzzz_1[j];

                    tg_yyyzz_yyyyzzzz_0[j] = pb_y * tg_yyzz_yyyyzzzz_0[j] + fr * tg_yyzz_yyyyzzzz_1[j] + fl1_fx * (tg_yzz_yyyyzzzz_0[j] - tg_yzz_yyyyzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyzz_yyyzzzz_1[j];

                    tg_yyyzz_yyyzzzzz_0[j] = pb_y * tg_yyzz_yyyzzzzz_0[j] + fr * tg_yyzz_yyyzzzzz_1[j] + fl1_fx * (tg_yzz_yyyzzzzz_0[j] - tg_yzz_yyyzzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzz_yyzzzzz_1[j];

                    tg_yyyzz_yyzzzzzz_0[j] = pb_y * tg_yyzz_yyzzzzzz_0[j] + fr * tg_yyzz_yyzzzzzz_1[j] + fl1_fx * (tg_yzz_yyzzzzzz_0[j] - tg_yzz_yyzzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyzz_yzzzzzz_1[j];

                    tg_yyyzz_yzzzzzzz_0[j] = pb_y * tg_yyzz_yzzzzzzz_0[j] + fr * tg_yyzz_yzzzzzzz_1[j] + fl1_fx * (tg_yzz_yzzzzzzz_0[j] - tg_yzz_yzzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzz_zzzzzzz_1[j];

                    tg_yyyzz_zzzzzzzz_0[j] = pb_y * tg_yyzz_zzzzzzzz_0[j] + fr * tg_yyzz_zzzzzzzz_1[j] + fl1_fx * (tg_yzz_zzzzzzzz_0[j] - tg_yzz_zzzzzzzz_1[j] * fl1_fza);

                    tg_yyzzz_xxxxxxxx_0[j] = pb_y * tg_yzzz_xxxxxxxx_0[j] + fr * tg_yzzz_xxxxxxxx_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxxxxx_0[j] - tg_zzz_xxxxxxxx_1[j] * fl1_fza);

                    tg_yyzzz_xxxxxxxy_0[j] = pb_y * tg_yzzz_xxxxxxxy_0[j] + fr * tg_yzzz_xxxxxxxy_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxxxxy_0[j] - tg_zzz_xxxxxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzz_xxxxxxx_1[j];

                    tg_yyzzz_xxxxxxxz_0[j] = pb_y * tg_yzzz_xxxxxxxz_0[j] + fr * tg_yzzz_xxxxxxxz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxxxxz_0[j] - tg_zzz_xxxxxxxz_1[j] * fl1_fza);

                    tg_yyzzz_xxxxxxyy_0[j] = pb_y * tg_yzzz_xxxxxxyy_0[j] + fr * tg_yzzz_xxxxxxyy_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxxxyy_0[j] - tg_zzz_xxxxxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yzzz_xxxxxxy_1[j];

                    tg_yyzzz_xxxxxxyz_0[j] = pb_y * tg_yzzz_xxxxxxyz_0[j] + fr * tg_yzzz_xxxxxxyz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxxxyz_0[j] - tg_zzz_xxxxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzz_xxxxxxz_1[j];

                    tg_yyzzz_xxxxxxzz_0[j] = pb_y * tg_yzzz_xxxxxxzz_0[j] + fr * tg_yzzz_xxxxxxzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxxxzz_0[j] - tg_zzz_xxxxxxzz_1[j] * fl1_fza);

                    tg_yyzzz_xxxxxyyy_0[j] = pb_y * tg_yzzz_xxxxxyyy_0[j] + fr * tg_yzzz_xxxxxyyy_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxxyyy_0[j] - tg_zzz_xxxxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzz_xxxxxyy_1[j];

                    tg_yyzzz_xxxxxyyz_0[j] = pb_y * tg_yzzz_xxxxxyyz_0[j] + fr * tg_yzzz_xxxxxyyz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxxyyz_0[j] - tg_zzz_xxxxxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yzzz_xxxxxyz_1[j];

                    tg_yyzzz_xxxxxyzz_0[j] = pb_y * tg_yzzz_xxxxxyzz_0[j] + fr * tg_yzzz_xxxxxyzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxxyzz_0[j] - tg_zzz_xxxxxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzz_xxxxxzz_1[j];

                    tg_yyzzz_xxxxxzzz_0[j] = pb_y * tg_yzzz_xxxxxzzz_0[j] + fr * tg_yzzz_xxxxxzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxxzzz_0[j] - tg_zzz_xxxxxzzz_1[j] * fl1_fza);

                    tg_yyzzz_xxxxyyyy_0[j] = pb_y * tg_yzzz_xxxxyyyy_0[j] + fr * tg_yzzz_xxxxyyyy_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxyyyy_0[j] - tg_zzz_xxxxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzzz_xxxxyyy_1[j];

                    tg_yyzzz_xxxxyyyz_0[j] = pb_y * tg_yzzz_xxxxyyyz_0[j] + fr * tg_yzzz_xxxxyyyz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxyyyz_0[j] - tg_zzz_xxxxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzz_xxxxyyz_1[j];

                    tg_yyzzz_xxxxyyzz_0[j] = pb_y * tg_yzzz_xxxxyyzz_0[j] + fr * tg_yzzz_xxxxyyzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxyyzz_0[j] - tg_zzz_xxxxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yzzz_xxxxyzz_1[j];

                    tg_yyzzz_xxxxyzzz_0[j] = pb_y * tg_yzzz_xxxxyzzz_0[j] + fr * tg_yzzz_xxxxyzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxyzzz_0[j] - tg_zzz_xxxxyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzz_xxxxzzz_1[j];

                    tg_yyzzz_xxxxzzzz_0[j] = pb_y * tg_yzzz_xxxxzzzz_0[j] + fr * tg_yzzz_xxxxzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxzzzz_0[j] - tg_zzz_xxxxzzzz_1[j] * fl1_fza);

                    tg_yyzzz_xxxyyyyy_0[j] = pb_y * tg_yzzz_xxxyyyyy_0[j] + fr * tg_yzzz_xxxyyyyy_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxyyyyy_0[j] - tg_zzz_xxxyyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yzzz_xxxyyyy_1[j];

                    tg_yyzzz_xxxyyyyz_0[j] = pb_y * tg_yzzz_xxxyyyyz_0[j] + fr * tg_yzzz_xxxyyyyz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxyyyyz_0[j] - tg_zzz_xxxyyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzzz_xxxyyyz_1[j];

                    tg_yyzzz_xxxyyyzz_0[j] = pb_y * tg_yzzz_xxxyyyzz_0[j] + fr * tg_yzzz_xxxyyyzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxyyyzz_0[j] - tg_zzz_xxxyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzz_xxxyyzz_1[j];

                    tg_yyzzz_xxxyyzzz_0[j] = pb_y * tg_yzzz_xxxyyzzz_0[j] + fr * tg_yzzz_xxxyyzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxyyzzz_0[j] - tg_zzz_xxxyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yzzz_xxxyzzz_1[j];

                    tg_yyzzz_xxxyzzzz_0[j] = pb_y * tg_yzzz_xxxyzzzz_0[j] + fr * tg_yzzz_xxxyzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxyzzzz_0[j] - tg_zzz_xxxyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzz_xxxzzzz_1[j];

                    tg_yyzzz_xxxzzzzz_0[j] = pb_y * tg_yzzz_xxxzzzzz_0[j] + fr * tg_yzzz_xxxzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxzzzzz_0[j] - tg_zzz_xxxzzzzz_1[j] * fl1_fza);

                    tg_yyzzz_xxyyyyyy_0[j] = pb_y * tg_yzzz_xxyyyyyy_0[j] + fr * tg_yzzz_xxyyyyyy_1[j] + 0.5 * fl1_fx * (tg_zzz_xxyyyyyy_0[j] - tg_zzz_xxyyyyyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yzzz_xxyyyyy_1[j];

                    tg_yyzzz_xxyyyyyz_0[j] = pb_y * tg_yzzz_xxyyyyyz_0[j] + fr * tg_yzzz_xxyyyyyz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxyyyyyz_0[j] - tg_zzz_xxyyyyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yzzz_xxyyyyz_1[j];

                    tg_yyzzz_xxyyyyzz_0[j] = pb_y * tg_yzzz_xxyyyyzz_0[j] + fr * tg_yzzz_xxyyyyzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxyyyyzz_0[j] - tg_zzz_xxyyyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzzz_xxyyyzz_1[j];

                    tg_yyzzz_xxyyyzzz_0[j] = pb_y * tg_yzzz_xxyyyzzz_0[j] + fr * tg_yzzz_xxyyyzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxyyyzzz_0[j] - tg_zzz_xxyyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzz_xxyyzzz_1[j];

                    tg_yyzzz_xxyyzzzz_0[j] = pb_y * tg_yzzz_xxyyzzzz_0[j] + fr * tg_yzzz_xxyyzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxyyzzzz_0[j] - tg_zzz_xxyyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yzzz_xxyzzzz_1[j];

                    tg_yyzzz_xxyzzzzz_0[j] = pb_y * tg_yzzz_xxyzzzzz_0[j] + fr * tg_yzzz_xxyzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxyzzzzz_0[j] - tg_zzz_xxyzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzz_xxzzzzz_1[j];

                    tg_yyzzz_xxzzzzzz_0[j] = pb_y * tg_yzzz_xxzzzzzz_0[j] + fr * tg_yzzz_xxzzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxzzzzzz_0[j] - tg_zzz_xxzzzzzz_1[j] * fl1_fza);

                    tg_yyzzz_xyyyyyyy_0[j] = pb_y * tg_yzzz_xyyyyyyy_0[j] + fr * tg_yzzz_xyyyyyyy_1[j] + 0.5 * fl1_fx * (tg_zzz_xyyyyyyy_0[j] - tg_zzz_xyyyyyyy_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_yzzz_xyyyyyy_1[j];

                    tg_yyzzz_xyyyyyyz_0[j] = pb_y * tg_yzzz_xyyyyyyz_0[j] + fr * tg_yzzz_xyyyyyyz_1[j] + 0.5 * fl1_fx * (tg_zzz_xyyyyyyz_0[j] - tg_zzz_xyyyyyyz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yzzz_xyyyyyz_1[j];

                    tg_yyzzz_xyyyyyzz_0[j] = pb_y * tg_yzzz_xyyyyyzz_0[j] + fr * tg_yzzz_xyyyyyzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xyyyyyzz_0[j] - tg_zzz_xyyyyyzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yzzz_xyyyyzz_1[j];

                    tg_yyzzz_xyyyyzzz_0[j] = pb_y * tg_yzzz_xyyyyzzz_0[j] + fr * tg_yzzz_xyyyyzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xyyyyzzz_0[j] - tg_zzz_xyyyyzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzzz_xyyyzzz_1[j];

                    tg_yyzzz_xyyyzzzz_0[j] = pb_y * tg_yzzz_xyyyzzzz_0[j] + fr * tg_yzzz_xyyyzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xyyyzzzz_0[j] - tg_zzz_xyyyzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzz_xyyzzzz_1[j];

                    tg_yyzzz_xyyzzzzz_0[j] = pb_y * tg_yzzz_xyyzzzzz_0[j] + fr * tg_yzzz_xyyzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xyyzzzzz_0[j] - tg_zzz_xyyzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yzzz_xyzzzzz_1[j];

                    tg_yyzzz_xyzzzzzz_0[j] = pb_y * tg_yzzz_xyzzzzzz_0[j] + fr * tg_yzzz_xyzzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xyzzzzzz_0[j] - tg_zzz_xyzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzz_xzzzzzz_1[j];

                    tg_yyzzz_xzzzzzzz_0[j] = pb_y * tg_yzzz_xzzzzzzz_0[j] + fr * tg_yzzz_xzzzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xzzzzzzz_0[j] - tg_zzz_xzzzzzzz_1[j] * fl1_fza);

                    tg_yyzzz_yyyyyyyy_0[j] = pb_y * tg_yzzz_yyyyyyyy_0[j] + fr * tg_yzzz_yyyyyyyy_1[j] + 0.5 * fl1_fx * (tg_zzz_yyyyyyyy_0[j] - tg_zzz_yyyyyyyy_1[j] * fl1_fza) + 4.0 * fl1_fxn * tg_yzzz_yyyyyyy_1[j];

                    tg_yyzzz_yyyyyyyz_0[j] = pb_y * tg_yzzz_yyyyyyyz_0[j] + fr * tg_yzzz_yyyyyyyz_1[j] + 0.5 * fl1_fx * (tg_zzz_yyyyyyyz_0[j] - tg_zzz_yyyyyyyz_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_yzzz_yyyyyyz_1[j];

                    tg_yyzzz_yyyyyyzz_0[j] = pb_y * tg_yzzz_yyyyyyzz_0[j] + fr * tg_yzzz_yyyyyyzz_1[j] + 0.5 * fl1_fx * (tg_zzz_yyyyyyzz_0[j] - tg_zzz_yyyyyyzz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yzzz_yyyyyzz_1[j];

                    tg_yyzzz_yyyyyzzz_0[j] = pb_y * tg_yzzz_yyyyyzzz_0[j] + fr * tg_yzzz_yyyyyzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_yyyyyzzz_0[j] - tg_zzz_yyyyyzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yzzz_yyyyzzz_1[j];

                    tg_yyzzz_yyyyzzzz_0[j] = pb_y * tg_yzzz_yyyyzzzz_0[j] + fr * tg_yzzz_yyyyzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_yyyyzzzz_0[j] - tg_zzz_yyyyzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzzz_yyyzzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSL_851_945(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (851,945)

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
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_8_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_8_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_8_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_8_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_8_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {8, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_7_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {7, -1, -1, -1}, 
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

                auto tg_yzzz_yyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 626); 

                auto tg_yzzz_yyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 627); 

                auto tg_yzzz_yzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 628); 

                auto tg_yzzz_zzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 629); 

                auto tg_zzzz_xxxxxxxx_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 630); 

                auto tg_zzzz_xxxxxxxy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 631); 

                auto tg_zzzz_xxxxxxxz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 632); 

                auto tg_zzzz_xxxxxxyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 633); 

                auto tg_zzzz_xxxxxxyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 634); 

                auto tg_zzzz_xxxxxxzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 635); 

                auto tg_zzzz_xxxxxyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 636); 

                auto tg_zzzz_xxxxxyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 637); 

                auto tg_zzzz_xxxxxyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 638); 

                auto tg_zzzz_xxxxxzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 639); 

                auto tg_zzzz_xxxxyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 640); 

                auto tg_zzzz_xxxxyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 641); 

                auto tg_zzzz_xxxxyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 642); 

                auto tg_zzzz_xxxxyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 643); 

                auto tg_zzzz_xxxxzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 644); 

                auto tg_zzzz_xxxyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 645); 

                auto tg_zzzz_xxxyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 646); 

                auto tg_zzzz_xxxyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 647); 

                auto tg_zzzz_xxxyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 648); 

                auto tg_zzzz_xxxyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 649); 

                auto tg_zzzz_xxxzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 650); 

                auto tg_zzzz_xxyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 651); 

                auto tg_zzzz_xxyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 652); 

                auto tg_zzzz_xxyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 653); 

                auto tg_zzzz_xxyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 654); 

                auto tg_zzzz_xxyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 655); 

                auto tg_zzzz_xxyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 656); 

                auto tg_zzzz_xxzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 657); 

                auto tg_zzzz_xyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 658); 

                auto tg_zzzz_xyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 659); 

                auto tg_zzzz_xyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 660); 

                auto tg_zzzz_xyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 661); 

                auto tg_zzzz_xyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 662); 

                auto tg_zzzz_xyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 663); 

                auto tg_zzzz_xyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 664); 

                auto tg_zzzz_xzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 665); 

                auto tg_zzzz_yyyyyyyy_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 666); 

                auto tg_zzzz_yyyyyyyz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 667); 

                auto tg_zzzz_yyyyyyzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 668); 

                auto tg_zzzz_yyyyyzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 669); 

                auto tg_zzzz_yyyyzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 670); 

                auto tg_zzzz_yyyzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 671); 

                auto tg_zzzz_yyzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 672); 

                auto tg_zzzz_yzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 673); 

                auto tg_zzzz_zzzzzzzz_0 = primBuffer[pidx_g_4_8_m0].data(675 * idx + 674); 

                auto tg_yzzz_yyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 626); 

                auto tg_yzzz_yyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 627); 

                auto tg_yzzz_yzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 628); 

                auto tg_yzzz_zzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 629); 

                auto tg_zzzz_xxxxxxxx_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 630); 

                auto tg_zzzz_xxxxxxxy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 631); 

                auto tg_zzzz_xxxxxxxz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 632); 

                auto tg_zzzz_xxxxxxyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 633); 

                auto tg_zzzz_xxxxxxyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 634); 

                auto tg_zzzz_xxxxxxzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 635); 

                auto tg_zzzz_xxxxxyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 636); 

                auto tg_zzzz_xxxxxyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 637); 

                auto tg_zzzz_xxxxxyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 638); 

                auto tg_zzzz_xxxxxzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 639); 

                auto tg_zzzz_xxxxyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 640); 

                auto tg_zzzz_xxxxyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 641); 

                auto tg_zzzz_xxxxyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 642); 

                auto tg_zzzz_xxxxyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 643); 

                auto tg_zzzz_xxxxzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 644); 

                auto tg_zzzz_xxxyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 645); 

                auto tg_zzzz_xxxyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 646); 

                auto tg_zzzz_xxxyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 647); 

                auto tg_zzzz_xxxyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 648); 

                auto tg_zzzz_xxxyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 649); 

                auto tg_zzzz_xxxzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 650); 

                auto tg_zzzz_xxyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 651); 

                auto tg_zzzz_xxyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 652); 

                auto tg_zzzz_xxyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 653); 

                auto tg_zzzz_xxyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 654); 

                auto tg_zzzz_xxyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 655); 

                auto tg_zzzz_xxyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 656); 

                auto tg_zzzz_xxzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 657); 

                auto tg_zzzz_xyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 658); 

                auto tg_zzzz_xyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 659); 

                auto tg_zzzz_xyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 660); 

                auto tg_zzzz_xyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 661); 

                auto tg_zzzz_xyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 662); 

                auto tg_zzzz_xyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 663); 

                auto tg_zzzz_xyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 664); 

                auto tg_zzzz_xzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 665); 

                auto tg_zzzz_yyyyyyyy_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 666); 

                auto tg_zzzz_yyyyyyyz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 667); 

                auto tg_zzzz_yyyyyyzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 668); 

                auto tg_zzzz_yyyyyzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 669); 

                auto tg_zzzz_yyyyzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 670); 

                auto tg_zzzz_yyyzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 671); 

                auto tg_zzzz_yyzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 672); 

                auto tg_zzzz_yzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 673); 

                auto tg_zzzz_zzzzzzzz_1 = primBuffer[pidx_g_4_8_m1].data(675 * idx + 674); 

                auto tg_zzz_xxxxxxxx_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 405); 

                auto tg_zzz_xxxxxxxy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 406); 

                auto tg_zzz_xxxxxxxz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 407); 

                auto tg_zzz_xxxxxxyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 408); 

                auto tg_zzz_xxxxxxyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 409); 

                auto tg_zzz_xxxxxxzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 410); 

                auto tg_zzz_xxxxxyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 411); 

                auto tg_zzz_xxxxxyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 412); 

                auto tg_zzz_xxxxxyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 413); 

                auto tg_zzz_xxxxxzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 414); 

                auto tg_zzz_xxxxyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 415); 

                auto tg_zzz_xxxxyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 416); 

                auto tg_zzz_xxxxyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 417); 

                auto tg_zzz_xxxxyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 418); 

                auto tg_zzz_xxxxzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 419); 

                auto tg_zzz_xxxyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 420); 

                auto tg_zzz_xxxyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 421); 

                auto tg_zzz_xxxyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 422); 

                auto tg_zzz_xxxyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 423); 

                auto tg_zzz_xxxyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 424); 

                auto tg_zzz_xxxzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 425); 

                auto tg_zzz_xxyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 426); 

                auto tg_zzz_xxyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 427); 

                auto tg_zzz_xxyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 428); 

                auto tg_zzz_xxyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 429); 

                auto tg_zzz_xxyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 430); 

                auto tg_zzz_xxyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 431); 

                auto tg_zzz_xxzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 432); 

                auto tg_zzz_xyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 433); 

                auto tg_zzz_xyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 434); 

                auto tg_zzz_xyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 435); 

                auto tg_zzz_xyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 436); 

                auto tg_zzz_xyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 437); 

                auto tg_zzz_xyyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 438); 

                auto tg_zzz_xyzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 439); 

                auto tg_zzz_xzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 440); 

                auto tg_zzz_yyyyyyyy_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 441); 

                auto tg_zzz_yyyyyyyz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 442); 

                auto tg_zzz_yyyyyyzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 443); 

                auto tg_zzz_yyyyyzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 444); 

                auto tg_zzz_yyyyzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 445); 

                auto tg_zzz_yyyzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 446); 

                auto tg_zzz_yyzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 447); 

                auto tg_zzz_yzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 448); 

                auto tg_zzz_zzzzzzzz_0 = primBuffer[pidx_g_3_8_m0].data(450 * idx + 449); 

                auto tg_zzz_xxxxxxxx_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 405); 

                auto tg_zzz_xxxxxxxy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 406); 

                auto tg_zzz_xxxxxxxz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 407); 

                auto tg_zzz_xxxxxxyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 408); 

                auto tg_zzz_xxxxxxyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 409); 

                auto tg_zzz_xxxxxxzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 410); 

                auto tg_zzz_xxxxxyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 411); 

                auto tg_zzz_xxxxxyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 412); 

                auto tg_zzz_xxxxxyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 413); 

                auto tg_zzz_xxxxxzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 414); 

                auto tg_zzz_xxxxyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 415); 

                auto tg_zzz_xxxxyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 416); 

                auto tg_zzz_xxxxyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 417); 

                auto tg_zzz_xxxxyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 418); 

                auto tg_zzz_xxxxzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 419); 

                auto tg_zzz_xxxyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 420); 

                auto tg_zzz_xxxyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 421); 

                auto tg_zzz_xxxyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 422); 

                auto tg_zzz_xxxyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 423); 

                auto tg_zzz_xxxyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 424); 

                auto tg_zzz_xxxzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 425); 

                auto tg_zzz_xxyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 426); 

                auto tg_zzz_xxyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 427); 

                auto tg_zzz_xxyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 428); 

                auto tg_zzz_xxyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 429); 

                auto tg_zzz_xxyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 430); 

                auto tg_zzz_xxyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 431); 

                auto tg_zzz_xxzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 432); 

                auto tg_zzz_xyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 433); 

                auto tg_zzz_xyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 434); 

                auto tg_zzz_xyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 435); 

                auto tg_zzz_xyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 436); 

                auto tg_zzz_xyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 437); 

                auto tg_zzz_xyyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 438); 

                auto tg_zzz_xyzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 439); 

                auto tg_zzz_xzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 440); 

                auto tg_zzz_yyyyyyyy_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 441); 

                auto tg_zzz_yyyyyyyz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 442); 

                auto tg_zzz_yyyyyyzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 443); 

                auto tg_zzz_yyyyyzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 444); 

                auto tg_zzz_yyyyzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 445); 

                auto tg_zzz_yyyzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 446); 

                auto tg_zzz_yyzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 447); 

                auto tg_zzz_yzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 448); 

                auto tg_zzz_zzzzzzzz_1 = primBuffer[pidx_g_3_8_m1].data(450 * idx + 449); 

                auto tg_yzzz_yyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 501); 

                auto tg_yzzz_yzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 502); 

                auto tg_yzzz_zzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 503); 

                auto tg_zzzz_xxxxxxx_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 504); 

                auto tg_zzzz_xxxxxxy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 505); 

                auto tg_zzzz_xxxxxxz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 506); 

                auto tg_zzzz_xxxxxyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 507); 

                auto tg_zzzz_xxxxxyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 508); 

                auto tg_zzzz_xxxxxzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 509); 

                auto tg_zzzz_xxxxyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 510); 

                auto tg_zzzz_xxxxyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 511); 

                auto tg_zzzz_xxxxyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 512); 

                auto tg_zzzz_xxxxzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 513); 

                auto tg_zzzz_xxxyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 514); 

                auto tg_zzzz_xxxyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 515); 

                auto tg_zzzz_xxxyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 516); 

                auto tg_zzzz_xxxyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 517); 

                auto tg_zzzz_xxxzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 518); 

                auto tg_zzzz_xxyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 519); 

                auto tg_zzzz_xxyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 520); 

                auto tg_zzzz_xxyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 521); 

                auto tg_zzzz_xxyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 522); 

                auto tg_zzzz_xxyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 523); 

                auto tg_zzzz_xxzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 524); 

                auto tg_zzzz_xyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 525); 

                auto tg_zzzz_xyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 526); 

                auto tg_zzzz_xyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 527); 

                auto tg_zzzz_xyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 528); 

                auto tg_zzzz_xyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 529); 

                auto tg_zzzz_xyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 530); 

                auto tg_zzzz_xzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 531); 

                auto tg_zzzz_yyyyyyy_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 532); 

                auto tg_zzzz_yyyyyyz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 533); 

                auto tg_zzzz_yyyyyzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 534); 

                auto tg_zzzz_yyyyzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 535); 

                auto tg_zzzz_yyyzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 536); 

                auto tg_zzzz_yyzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 537); 

                auto tg_zzzz_yzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 538); 

                auto tg_zzzz_zzzzzzz_1 = primBuffer[pidx_g_4_7_m1].data(540 * idx + 539); 

                // set up pointers to integrals

                auto tg_yyzzz_yyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 851); 

                auto tg_yyzzz_yyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 852); 

                auto tg_yyzzz_yzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 853); 

                auto tg_yyzzz_zzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 854); 

                auto tg_yzzzz_xxxxxxxx_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 855); 

                auto tg_yzzzz_xxxxxxxy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 856); 

                auto tg_yzzzz_xxxxxxxz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 857); 

                auto tg_yzzzz_xxxxxxyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 858); 

                auto tg_yzzzz_xxxxxxyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 859); 

                auto tg_yzzzz_xxxxxxzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 860); 

                auto tg_yzzzz_xxxxxyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 861); 

                auto tg_yzzzz_xxxxxyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 862); 

                auto tg_yzzzz_xxxxxyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 863); 

                auto tg_yzzzz_xxxxxzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 864); 

                auto tg_yzzzz_xxxxyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 865); 

                auto tg_yzzzz_xxxxyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 866); 

                auto tg_yzzzz_xxxxyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 867); 

                auto tg_yzzzz_xxxxyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 868); 

                auto tg_yzzzz_xxxxzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 869); 

                auto tg_yzzzz_xxxyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 870); 

                auto tg_yzzzz_xxxyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 871); 

                auto tg_yzzzz_xxxyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 872); 

                auto tg_yzzzz_xxxyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 873); 

                auto tg_yzzzz_xxxyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 874); 

                auto tg_yzzzz_xxxzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 875); 

                auto tg_yzzzz_xxyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 876); 

                auto tg_yzzzz_xxyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 877); 

                auto tg_yzzzz_xxyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 878); 

                auto tg_yzzzz_xxyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 879); 

                auto tg_yzzzz_xxyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 880); 

                auto tg_yzzzz_xxyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 881); 

                auto tg_yzzzz_xxzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 882); 

                auto tg_yzzzz_xyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 883); 

                auto tg_yzzzz_xyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 884); 

                auto tg_yzzzz_xyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 885); 

                auto tg_yzzzz_xyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 886); 

                auto tg_yzzzz_xyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 887); 

                auto tg_yzzzz_xyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 888); 

                auto tg_yzzzz_xyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 889); 

                auto tg_yzzzz_xzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 890); 

                auto tg_yzzzz_yyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 891); 

                auto tg_yzzzz_yyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 892); 

                auto tg_yzzzz_yyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 893); 

                auto tg_yzzzz_yyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 894); 

                auto tg_yzzzz_yyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 895); 

                auto tg_yzzzz_yyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 896); 

                auto tg_yzzzz_yyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 897); 

                auto tg_yzzzz_yzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 898); 

                auto tg_yzzzz_zzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 899); 

                auto tg_zzzzz_xxxxxxxx_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 900); 

                auto tg_zzzzz_xxxxxxxy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 901); 

                auto tg_zzzzz_xxxxxxxz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 902); 

                auto tg_zzzzz_xxxxxxyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 903); 

                auto tg_zzzzz_xxxxxxyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 904); 

                auto tg_zzzzz_xxxxxxzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 905); 

                auto tg_zzzzz_xxxxxyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 906); 

                auto tg_zzzzz_xxxxxyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 907); 

                auto tg_zzzzz_xxxxxyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 908); 

                auto tg_zzzzz_xxxxxzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 909); 

                auto tg_zzzzz_xxxxyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 910); 

                auto tg_zzzzz_xxxxyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 911); 

                auto tg_zzzzz_xxxxyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 912); 

                auto tg_zzzzz_xxxxyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 913); 

                auto tg_zzzzz_xxxxzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 914); 

                auto tg_zzzzz_xxxyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 915); 

                auto tg_zzzzz_xxxyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 916); 

                auto tg_zzzzz_xxxyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 917); 

                auto tg_zzzzz_xxxyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 918); 

                auto tg_zzzzz_xxxyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 919); 

                auto tg_zzzzz_xxxzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 920); 

                auto tg_zzzzz_xxyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 921); 

                auto tg_zzzzz_xxyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 922); 

                auto tg_zzzzz_xxyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 923); 

                auto tg_zzzzz_xxyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 924); 

                auto tg_zzzzz_xxyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 925); 

                auto tg_zzzzz_xxyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 926); 

                auto tg_zzzzz_xxzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 927); 

                auto tg_zzzzz_xyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 928); 

                auto tg_zzzzz_xyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 929); 

                auto tg_zzzzz_xyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 930); 

                auto tg_zzzzz_xyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 931); 

                auto tg_zzzzz_xyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 932); 

                auto tg_zzzzz_xyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 933); 

                auto tg_zzzzz_xyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 934); 

                auto tg_zzzzz_xzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 935); 

                auto tg_zzzzz_yyyyyyyy_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 936); 

                auto tg_zzzzz_yyyyyyyz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 937); 

                auto tg_zzzzz_yyyyyyzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 938); 

                auto tg_zzzzz_yyyyyzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 939); 

                auto tg_zzzzz_yyyyzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 940); 

                auto tg_zzzzz_yyyzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 941); 

                auto tg_zzzzz_yyzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 942); 

                auto tg_zzzzz_yzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 943); 

                auto tg_zzzzz_zzzzzzzz_0 = primBuffer[pidx_g_5_8_m0].data(945 * idx + 944); 

                // Batch of Integrals (851,945)

                #pragma omp simd aligned(fxn, fza, tg_yyzzz_yyyzzzzz_0, tg_yyzzz_yyzzzzzz_0, \
                                         tg_yyzzz_yzzzzzzz_0, tg_yyzzz_zzzzzzzz_0, tg_yzzz_yyyzzzzz_0, tg_yzzz_yyyzzzzz_1, \
                                         tg_yzzz_yyzzzzz_1, tg_yzzz_yyzzzzzz_0, tg_yzzz_yyzzzzzz_1, tg_yzzz_yzzzzzz_1, \
                                         tg_yzzz_yzzzzzzz_0, tg_yzzz_yzzzzzzz_1, tg_yzzz_zzzzzzz_1, tg_yzzz_zzzzzzzz_0, \
                                         tg_yzzz_zzzzzzzz_1, tg_yzzzz_xxxxxxxx_0, tg_yzzzz_xxxxxxxy_0, tg_yzzzz_xxxxxxxz_0, \
                                         tg_yzzzz_xxxxxxyy_0, tg_yzzzz_xxxxxxyz_0, tg_yzzzz_xxxxxxzz_0, tg_yzzzz_xxxxxyyy_0, \
                                         tg_yzzzz_xxxxxyyz_0, tg_yzzzz_xxxxxyzz_0, tg_yzzzz_xxxxxzzz_0, tg_yzzzz_xxxxyyyy_0, \
                                         tg_yzzzz_xxxxyyyz_0, tg_yzzzz_xxxxyyzz_0, tg_yzzzz_xxxxyzzz_0, tg_yzzzz_xxxxzzzz_0, \
                                         tg_yzzzz_xxxyyyyy_0, tg_yzzzz_xxxyyyyz_0, tg_yzzzz_xxxyyyzz_0, tg_yzzzz_xxxyyzzz_0, \
                                         tg_yzzzz_xxxyzzzz_0, tg_yzzzz_xxxzzzzz_0, tg_yzzzz_xxyyyyyy_0, tg_yzzzz_xxyyyyyz_0, \
                                         tg_yzzzz_xxyyyyzz_0, tg_yzzzz_xxyyyzzz_0, tg_yzzzz_xxyyzzzz_0, tg_yzzzz_xxyzzzzz_0, \
                                         tg_yzzzz_xxzzzzzz_0, tg_yzzzz_xyyyyyyy_0, tg_yzzzz_xyyyyyyz_0, tg_yzzzz_xyyyyyzz_0, \
                                         tg_yzzzz_xyyyyzzz_0, tg_yzzzz_xyyyzzzz_0, tg_yzzzz_xyyzzzzz_0, tg_yzzzz_xyzzzzzz_0, \
                                         tg_yzzzz_xzzzzzzz_0, tg_yzzzz_yyyyyyyy_0, tg_yzzzz_yyyyyyyz_0, tg_yzzzz_yyyyyyzz_0, \
                                         tg_yzzzz_yyyyyzzz_0, tg_yzzzz_yyyyzzzz_0, tg_yzzzz_yyyzzzzz_0, tg_yzzzz_yyzzzzzz_0, \
                                         tg_yzzzz_yzzzzzzz_0, tg_yzzzz_zzzzzzzz_0, tg_zzz_xxxxxxxx_0, tg_zzz_xxxxxxxx_1, \
                                         tg_zzz_xxxxxxxy_0, tg_zzz_xxxxxxxy_1, tg_zzz_xxxxxxxz_0, tg_zzz_xxxxxxxz_1, \
                                         tg_zzz_xxxxxxyy_0, tg_zzz_xxxxxxyy_1, tg_zzz_xxxxxxyz_0, tg_zzz_xxxxxxyz_1, \
                                         tg_zzz_xxxxxxzz_0, tg_zzz_xxxxxxzz_1, tg_zzz_xxxxxyyy_0, tg_zzz_xxxxxyyy_1, \
                                         tg_zzz_xxxxxyyz_0, tg_zzz_xxxxxyyz_1, tg_zzz_xxxxxyzz_0, tg_zzz_xxxxxyzz_1, \
                                         tg_zzz_xxxxxzzz_0, tg_zzz_xxxxxzzz_1, tg_zzz_xxxxyyyy_0, tg_zzz_xxxxyyyy_1, \
                                         tg_zzz_xxxxyyyz_0, tg_zzz_xxxxyyyz_1, tg_zzz_xxxxyyzz_0, tg_zzz_xxxxyyzz_1, \
                                         tg_zzz_xxxxyzzz_0, tg_zzz_xxxxyzzz_1, tg_zzz_xxxxzzzz_0, tg_zzz_xxxxzzzz_1, \
                                         tg_zzz_xxxyyyyy_0, tg_zzz_xxxyyyyy_1, tg_zzz_xxxyyyyz_0, tg_zzz_xxxyyyyz_1, \
                                         tg_zzz_xxxyyyzz_0, tg_zzz_xxxyyyzz_1, tg_zzz_xxxyyzzz_0, tg_zzz_xxxyyzzz_1, \
                                         tg_zzz_xxxyzzzz_0, tg_zzz_xxxyzzzz_1, tg_zzz_xxxzzzzz_0, tg_zzz_xxxzzzzz_1, \
                                         tg_zzz_xxyyyyyy_0, tg_zzz_xxyyyyyy_1, tg_zzz_xxyyyyyz_0, tg_zzz_xxyyyyyz_1, \
                                         tg_zzz_xxyyyyzz_0, tg_zzz_xxyyyyzz_1, tg_zzz_xxyyyzzz_0, tg_zzz_xxyyyzzz_1, \
                                         tg_zzz_xxyyzzzz_0, tg_zzz_xxyyzzzz_1, tg_zzz_xxyzzzzz_0, tg_zzz_xxyzzzzz_1, \
                                         tg_zzz_xxzzzzzz_0, tg_zzz_xxzzzzzz_1, tg_zzz_xyyyyyyy_0, tg_zzz_xyyyyyyy_1, \
                                         tg_zzz_xyyyyyyz_0, tg_zzz_xyyyyyyz_1, tg_zzz_xyyyyyzz_0, tg_zzz_xyyyyyzz_1, \
                                         tg_zzz_xyyyyzzz_0, tg_zzz_xyyyyzzz_1, tg_zzz_xyyyzzzz_0, tg_zzz_xyyyzzzz_1, \
                                         tg_zzz_xyyzzzzz_0, tg_zzz_xyyzzzzz_1, tg_zzz_xyzzzzzz_0, tg_zzz_xyzzzzzz_1, \
                                         tg_zzz_xzzzzzzz_0, tg_zzz_xzzzzzzz_1, tg_zzz_yyyyyyyy_0, tg_zzz_yyyyyyyy_1, \
                                         tg_zzz_yyyyyyyz_0, tg_zzz_yyyyyyyz_1, tg_zzz_yyyyyyzz_0, tg_zzz_yyyyyyzz_1, \
                                         tg_zzz_yyyyyzzz_0, tg_zzz_yyyyyzzz_1, tg_zzz_yyyyzzzz_0, tg_zzz_yyyyzzzz_1, \
                                         tg_zzz_yyyzzzzz_0, tg_zzz_yyyzzzzz_1, tg_zzz_yyzzzzzz_0, tg_zzz_yyzzzzzz_1, \
                                         tg_zzz_yzzzzzzz_0, tg_zzz_yzzzzzzz_1, tg_zzz_zzzzzzzz_0, tg_zzz_zzzzzzzz_1, \
                                         tg_zzzz_xxxxxxx_1, tg_zzzz_xxxxxxxx_0, tg_zzzz_xxxxxxxx_1, tg_zzzz_xxxxxxxy_0, \
                                         tg_zzzz_xxxxxxxy_1, tg_zzzz_xxxxxxxz_0, tg_zzzz_xxxxxxxz_1, tg_zzzz_xxxxxxy_1, \
                                         tg_zzzz_xxxxxxyy_0, tg_zzzz_xxxxxxyy_1, tg_zzzz_xxxxxxyz_0, tg_zzzz_xxxxxxyz_1, \
                                         tg_zzzz_xxxxxxz_1, tg_zzzz_xxxxxxzz_0, tg_zzzz_xxxxxxzz_1, tg_zzzz_xxxxxyy_1, \
                                         tg_zzzz_xxxxxyyy_0, tg_zzzz_xxxxxyyy_1, tg_zzzz_xxxxxyyz_0, tg_zzzz_xxxxxyyz_1, \
                                         tg_zzzz_xxxxxyz_1, tg_zzzz_xxxxxyzz_0, tg_zzzz_xxxxxyzz_1, tg_zzzz_xxxxxzz_1, \
                                         tg_zzzz_xxxxxzzz_0, tg_zzzz_xxxxxzzz_1, tg_zzzz_xxxxyyy_1, tg_zzzz_xxxxyyyy_0, \
                                         tg_zzzz_xxxxyyyy_1, tg_zzzz_xxxxyyyz_0, tg_zzzz_xxxxyyyz_1, tg_zzzz_xxxxyyz_1, \
                                         tg_zzzz_xxxxyyzz_0, tg_zzzz_xxxxyyzz_1, tg_zzzz_xxxxyzz_1, tg_zzzz_xxxxyzzz_0, \
                                         tg_zzzz_xxxxyzzz_1, tg_zzzz_xxxxzzz_1, tg_zzzz_xxxxzzzz_0, tg_zzzz_xxxxzzzz_1, \
                                         tg_zzzz_xxxyyyy_1, tg_zzzz_xxxyyyyy_0, tg_zzzz_xxxyyyyy_1, tg_zzzz_xxxyyyyz_0, \
                                         tg_zzzz_xxxyyyyz_1, tg_zzzz_xxxyyyz_1, tg_zzzz_xxxyyyzz_0, tg_zzzz_xxxyyyzz_1, \
                                         tg_zzzz_xxxyyzz_1, tg_zzzz_xxxyyzzz_0, tg_zzzz_xxxyyzzz_1, tg_zzzz_xxxyzzz_1, \
                                         tg_zzzz_xxxyzzzz_0, tg_zzzz_xxxyzzzz_1, tg_zzzz_xxxzzzz_1, tg_zzzz_xxxzzzzz_0, \
                                         tg_zzzz_xxxzzzzz_1, tg_zzzz_xxyyyyy_1, tg_zzzz_xxyyyyyy_0, tg_zzzz_xxyyyyyy_1, \
                                         tg_zzzz_xxyyyyyz_0, tg_zzzz_xxyyyyyz_1, tg_zzzz_xxyyyyz_1, tg_zzzz_xxyyyyzz_0, \
                                         tg_zzzz_xxyyyyzz_1, tg_zzzz_xxyyyzz_1, tg_zzzz_xxyyyzzz_0, tg_zzzz_xxyyyzzz_1, \
                                         tg_zzzz_xxyyzzz_1, tg_zzzz_xxyyzzzz_0, tg_zzzz_xxyyzzzz_1, tg_zzzz_xxyzzzz_1, \
                                         tg_zzzz_xxyzzzzz_0, tg_zzzz_xxyzzzzz_1, tg_zzzz_xxzzzzz_1, tg_zzzz_xxzzzzzz_0, \
                                         tg_zzzz_xxzzzzzz_1, tg_zzzz_xyyyyyy_1, tg_zzzz_xyyyyyyy_0, tg_zzzz_xyyyyyyy_1, \
                                         tg_zzzz_xyyyyyyz_0, tg_zzzz_xyyyyyyz_1, tg_zzzz_xyyyyyz_1, tg_zzzz_xyyyyyzz_0, \
                                         tg_zzzz_xyyyyyzz_1, tg_zzzz_xyyyyzz_1, tg_zzzz_xyyyyzzz_0, tg_zzzz_xyyyyzzz_1, \
                                         tg_zzzz_xyyyzzz_1, tg_zzzz_xyyyzzzz_0, tg_zzzz_xyyyzzzz_1, tg_zzzz_xyyzzzz_1, \
                                         tg_zzzz_xyyzzzzz_0, tg_zzzz_xyyzzzzz_1, tg_zzzz_xyzzzzz_1, tg_zzzz_xyzzzzzz_0, \
                                         tg_zzzz_xyzzzzzz_1, tg_zzzz_xzzzzzz_1, tg_zzzz_xzzzzzzz_0, tg_zzzz_xzzzzzzz_1, \
                                         tg_zzzz_yyyyyyy_1, tg_zzzz_yyyyyyyy_0, tg_zzzz_yyyyyyyy_1, tg_zzzz_yyyyyyyz_0, \
                                         tg_zzzz_yyyyyyyz_1, tg_zzzz_yyyyyyz_1, tg_zzzz_yyyyyyzz_0, tg_zzzz_yyyyyyzz_1, \
                                         tg_zzzz_yyyyyzz_1, tg_zzzz_yyyyyzzz_0, tg_zzzz_yyyyyzzz_1, tg_zzzz_yyyyzzz_1, \
                                         tg_zzzz_yyyyzzzz_0, tg_zzzz_yyyyzzzz_1, tg_zzzz_yyyzzzz_1, tg_zzzz_yyyzzzzz_0, \
                                         tg_zzzz_yyyzzzzz_1, tg_zzzz_yyzzzzz_1, tg_zzzz_yyzzzzzz_0, tg_zzzz_yyzzzzzz_1, \
                                         tg_zzzz_yzzzzzz_1, tg_zzzz_yzzzzzzz_0, tg_zzzz_yzzzzzzz_1, tg_zzzz_zzzzzzz_1, \
                                         tg_zzzz_zzzzzzzz_0, tg_zzzz_zzzzzzzz_1, tg_zzzzz_xxxxxxxx_0, tg_zzzzz_xxxxxxxy_0, \
                                         tg_zzzzz_xxxxxxxz_0, tg_zzzzz_xxxxxxyy_0, tg_zzzzz_xxxxxxyz_0, tg_zzzzz_xxxxxxzz_0, \
                                         tg_zzzzz_xxxxxyyy_0, tg_zzzzz_xxxxxyyz_0, tg_zzzzz_xxxxxyzz_0, tg_zzzzz_xxxxxzzz_0, \
                                         tg_zzzzz_xxxxyyyy_0, tg_zzzzz_xxxxyyyz_0, tg_zzzzz_xxxxyyzz_0, tg_zzzzz_xxxxyzzz_0, \
                                         tg_zzzzz_xxxxzzzz_0, tg_zzzzz_xxxyyyyy_0, tg_zzzzz_xxxyyyyz_0, tg_zzzzz_xxxyyyzz_0, \
                                         tg_zzzzz_xxxyyzzz_0, tg_zzzzz_xxxyzzzz_0, tg_zzzzz_xxxzzzzz_0, tg_zzzzz_xxyyyyyy_0, \
                                         tg_zzzzz_xxyyyyyz_0, tg_zzzzz_xxyyyyzz_0, tg_zzzzz_xxyyyzzz_0, tg_zzzzz_xxyyzzzz_0, \
                                         tg_zzzzz_xxyzzzzz_0, tg_zzzzz_xxzzzzzz_0, tg_zzzzz_xyyyyyyy_0, tg_zzzzz_xyyyyyyz_0, \
                                         tg_zzzzz_xyyyyyzz_0, tg_zzzzz_xyyyyzzz_0, tg_zzzzz_xyyyzzzz_0, tg_zzzzz_xyyzzzzz_0, \
                                         tg_zzzzz_xyzzzzzz_0, tg_zzzzz_xzzzzzzz_0, tg_zzzzz_yyyyyyyy_0, tg_zzzzz_yyyyyyyz_0, \
                                         tg_zzzzz_yyyyyyzz_0, tg_zzzzz_yyyyyzzz_0, tg_zzzzz_yyyyzzzz_0, tg_zzzzz_yyyzzzzz_0, \
                                         tg_zzzzz_yyzzzzzz_0, tg_zzzzz_yzzzzzzz_0, tg_zzzzz_zzzzzzzz_0, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_y[j]; 

                    tg_yyzzz_yyyzzzzz_0[j] = pb_y * tg_yzzz_yyyzzzzz_0[j] + fr * tg_yzzz_yyyzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_yyyzzzzz_0[j] - tg_zzz_yyyzzzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzz_yyzzzzz_1[j];

                    tg_yyzzz_yyzzzzzz_0[j] = pb_y * tg_yzzz_yyzzzzzz_0[j] + fr * tg_yzzz_yyzzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_yyzzzzzz_0[j] - tg_zzz_yyzzzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yzzz_yzzzzzz_1[j];

                    tg_yyzzz_yzzzzzzz_0[j] = pb_y * tg_yzzz_yzzzzzzz_0[j] + fr * tg_yzzz_yzzzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_yzzzzzzz_0[j] - tg_zzz_yzzzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzz_zzzzzzz_1[j];

                    tg_yyzzz_zzzzzzzz_0[j] = pb_y * tg_yzzz_zzzzzzzz_0[j] + fr * tg_yzzz_zzzzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_zzzzzzzz_0[j] - tg_zzz_zzzzzzzz_1[j] * fl1_fza);

                    tg_yzzzz_xxxxxxxx_0[j] = pb_y * tg_zzzz_xxxxxxxx_0[j] + fr * tg_zzzz_xxxxxxxx_1[j];

                    tg_yzzzz_xxxxxxxy_0[j] = pb_y * tg_zzzz_xxxxxxxy_0[j] + fr * tg_zzzz_xxxxxxxy_1[j] + 0.5 * fl1_fxn * tg_zzzz_xxxxxxx_1[j];

                    tg_yzzzz_xxxxxxxz_0[j] = pb_y * tg_zzzz_xxxxxxxz_0[j] + fr * tg_zzzz_xxxxxxxz_1[j];

                    tg_yzzzz_xxxxxxyy_0[j] = pb_y * tg_zzzz_xxxxxxyy_0[j] + fr * tg_zzzz_xxxxxxyy_1[j] + fl1_fxn * tg_zzzz_xxxxxxy_1[j];

                    tg_yzzzz_xxxxxxyz_0[j] = pb_y * tg_zzzz_xxxxxxyz_0[j] + fr * tg_zzzz_xxxxxxyz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xxxxxxz_1[j];

                    tg_yzzzz_xxxxxxzz_0[j] = pb_y * tg_zzzz_xxxxxxzz_0[j] + fr * tg_zzzz_xxxxxxzz_1[j];

                    tg_yzzzz_xxxxxyyy_0[j] = pb_y * tg_zzzz_xxxxxyyy_0[j] + fr * tg_zzzz_xxxxxyyy_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxxxxyy_1[j];

                    tg_yzzzz_xxxxxyyz_0[j] = pb_y * tg_zzzz_xxxxxyyz_0[j] + fr * tg_zzzz_xxxxxyyz_1[j] + fl1_fxn * tg_zzzz_xxxxxyz_1[j];

                    tg_yzzzz_xxxxxyzz_0[j] = pb_y * tg_zzzz_xxxxxyzz_0[j] + fr * tg_zzzz_xxxxxyzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xxxxxzz_1[j];

                    tg_yzzzz_xxxxxzzz_0[j] = pb_y * tg_zzzz_xxxxxzzz_0[j] + fr * tg_zzzz_xxxxxzzz_1[j];

                    tg_yzzzz_xxxxyyyy_0[j] = pb_y * tg_zzzz_xxxxyyyy_0[j] + fr * tg_zzzz_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_zzzz_xxxxyyy_1[j];

                    tg_yzzzz_xxxxyyyz_0[j] = pb_y * tg_zzzz_xxxxyyyz_0[j] + fr * tg_zzzz_xxxxyyyz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxxxyyz_1[j];

                    tg_yzzzz_xxxxyyzz_0[j] = pb_y * tg_zzzz_xxxxyyzz_0[j] + fr * tg_zzzz_xxxxyyzz_1[j] + fl1_fxn * tg_zzzz_xxxxyzz_1[j];

                    tg_yzzzz_xxxxyzzz_0[j] = pb_y * tg_zzzz_xxxxyzzz_0[j] + fr * tg_zzzz_xxxxyzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xxxxzzz_1[j];

                    tg_yzzzz_xxxxzzzz_0[j] = pb_y * tg_zzzz_xxxxzzzz_0[j] + fr * tg_zzzz_xxxxzzzz_1[j];

                    tg_yzzzz_xxxyyyyy_0[j] = pb_y * tg_zzzz_xxxyyyyy_0[j] + fr * tg_zzzz_xxxyyyyy_1[j] + 2.5 * fl1_fxn * tg_zzzz_xxxyyyy_1[j];

                    tg_yzzzz_xxxyyyyz_0[j] = pb_y * tg_zzzz_xxxyyyyz_0[j] + fr * tg_zzzz_xxxyyyyz_1[j] + 2.0 * fl1_fxn * tg_zzzz_xxxyyyz_1[j];

                    tg_yzzzz_xxxyyyzz_0[j] = pb_y * tg_zzzz_xxxyyyzz_0[j] + fr * tg_zzzz_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxxyyzz_1[j];

                    tg_yzzzz_xxxyyzzz_0[j] = pb_y * tg_zzzz_xxxyyzzz_0[j] + fr * tg_zzzz_xxxyyzzz_1[j] + fl1_fxn * tg_zzzz_xxxyzzz_1[j];

                    tg_yzzzz_xxxyzzzz_0[j] = pb_y * tg_zzzz_xxxyzzzz_0[j] + fr * tg_zzzz_xxxyzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xxxzzzz_1[j];

                    tg_yzzzz_xxxzzzzz_0[j] = pb_y * tg_zzzz_xxxzzzzz_0[j] + fr * tg_zzzz_xxxzzzzz_1[j];

                    tg_yzzzz_xxyyyyyy_0[j] = pb_y * tg_zzzz_xxyyyyyy_0[j] + fr * tg_zzzz_xxyyyyyy_1[j] + 3.0 * fl1_fxn * tg_zzzz_xxyyyyy_1[j];

                    tg_yzzzz_xxyyyyyz_0[j] = pb_y * tg_zzzz_xxyyyyyz_0[j] + fr * tg_zzzz_xxyyyyyz_1[j] + 2.5 * fl1_fxn * tg_zzzz_xxyyyyz_1[j];

                    tg_yzzzz_xxyyyyzz_0[j] = pb_y * tg_zzzz_xxyyyyzz_0[j] + fr * tg_zzzz_xxyyyyzz_1[j] + 2.0 * fl1_fxn * tg_zzzz_xxyyyzz_1[j];

                    tg_yzzzz_xxyyyzzz_0[j] = pb_y * tg_zzzz_xxyyyzzz_0[j] + fr * tg_zzzz_xxyyyzzz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxyyzzz_1[j];

                    tg_yzzzz_xxyyzzzz_0[j] = pb_y * tg_zzzz_xxyyzzzz_0[j] + fr * tg_zzzz_xxyyzzzz_1[j] + fl1_fxn * tg_zzzz_xxyzzzz_1[j];

                    tg_yzzzz_xxyzzzzz_0[j] = pb_y * tg_zzzz_xxyzzzzz_0[j] + fr * tg_zzzz_xxyzzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xxzzzzz_1[j];

                    tg_yzzzz_xxzzzzzz_0[j] = pb_y * tg_zzzz_xxzzzzzz_0[j] + fr * tg_zzzz_xxzzzzzz_1[j];

                    tg_yzzzz_xyyyyyyy_0[j] = pb_y * tg_zzzz_xyyyyyyy_0[j] + fr * tg_zzzz_xyyyyyyy_1[j] + 3.5 * fl1_fxn * tg_zzzz_xyyyyyy_1[j];

                    tg_yzzzz_xyyyyyyz_0[j] = pb_y * tg_zzzz_xyyyyyyz_0[j] + fr * tg_zzzz_xyyyyyyz_1[j] + 3.0 * fl1_fxn * tg_zzzz_xyyyyyz_1[j];

                    tg_yzzzz_xyyyyyzz_0[j] = pb_y * tg_zzzz_xyyyyyzz_0[j] + fr * tg_zzzz_xyyyyyzz_1[j] + 2.5 * fl1_fxn * tg_zzzz_xyyyyzz_1[j];

                    tg_yzzzz_xyyyyzzz_0[j] = pb_y * tg_zzzz_xyyyyzzz_0[j] + fr * tg_zzzz_xyyyyzzz_1[j] + 2.0 * fl1_fxn * tg_zzzz_xyyyzzz_1[j];

                    tg_yzzzz_xyyyzzzz_0[j] = pb_y * tg_zzzz_xyyyzzzz_0[j] + fr * tg_zzzz_xyyyzzzz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xyyzzzz_1[j];

                    tg_yzzzz_xyyzzzzz_0[j] = pb_y * tg_zzzz_xyyzzzzz_0[j] + fr * tg_zzzz_xyyzzzzz_1[j] + fl1_fxn * tg_zzzz_xyzzzzz_1[j];

                    tg_yzzzz_xyzzzzzz_0[j] = pb_y * tg_zzzz_xyzzzzzz_0[j] + fr * tg_zzzz_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xzzzzzz_1[j];

                    tg_yzzzz_xzzzzzzz_0[j] = pb_y * tg_zzzz_xzzzzzzz_0[j] + fr * tg_zzzz_xzzzzzzz_1[j];

                    tg_yzzzz_yyyyyyyy_0[j] = pb_y * tg_zzzz_yyyyyyyy_0[j] + fr * tg_zzzz_yyyyyyyy_1[j] + 4.0 * fl1_fxn * tg_zzzz_yyyyyyy_1[j];

                    tg_yzzzz_yyyyyyyz_0[j] = pb_y * tg_zzzz_yyyyyyyz_0[j] + fr * tg_zzzz_yyyyyyyz_1[j] + 3.5 * fl1_fxn * tg_zzzz_yyyyyyz_1[j];

                    tg_yzzzz_yyyyyyzz_0[j] = pb_y * tg_zzzz_yyyyyyzz_0[j] + fr * tg_zzzz_yyyyyyzz_1[j] + 3.0 * fl1_fxn * tg_zzzz_yyyyyzz_1[j];

                    tg_yzzzz_yyyyyzzz_0[j] = pb_y * tg_zzzz_yyyyyzzz_0[j] + fr * tg_zzzz_yyyyyzzz_1[j] + 2.5 * fl1_fxn * tg_zzzz_yyyyzzz_1[j];

                    tg_yzzzz_yyyyzzzz_0[j] = pb_y * tg_zzzz_yyyyzzzz_0[j] + fr * tg_zzzz_yyyyzzzz_1[j] + 2.0 * fl1_fxn * tg_zzzz_yyyzzzz_1[j];

                    tg_yzzzz_yyyzzzzz_0[j] = pb_y * tg_zzzz_yyyzzzzz_0[j] + fr * tg_zzzz_yyyzzzzz_1[j] + 1.5 * fl1_fxn * tg_zzzz_yyzzzzz_1[j];

                    tg_yzzzz_yyzzzzzz_0[j] = pb_y * tg_zzzz_yyzzzzzz_0[j] + fr * tg_zzzz_yyzzzzzz_1[j] + fl1_fxn * tg_zzzz_yzzzzzz_1[j];

                    tg_yzzzz_yzzzzzzz_0[j] = pb_y * tg_zzzz_yzzzzzzz_0[j] + fr * tg_zzzz_yzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_zzzzzzz_1[j];

                    tg_yzzzz_zzzzzzzz_0[j] = pb_y * tg_zzzz_zzzzzzzz_0[j] + fr * tg_zzzz_zzzzzzzz_1[j];

                    fr = wp_z[j]; 

                    tg_zzzzz_xxxxxxxx_0[j] = pb_z * tg_zzzz_xxxxxxxx_0[j] + fr * tg_zzzz_xxxxxxxx_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxxxxxx_0[j] - tg_zzz_xxxxxxxx_1[j] * fl1_fza);

                    tg_zzzzz_xxxxxxxy_0[j] = pb_z * tg_zzzz_xxxxxxxy_0[j] + fr * tg_zzzz_xxxxxxxy_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxxxxxy_0[j] - tg_zzz_xxxxxxxy_1[j] * fl1_fza);

                    tg_zzzzz_xxxxxxxz_0[j] = pb_z * tg_zzzz_xxxxxxxz_0[j] + fr * tg_zzzz_xxxxxxxz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxxxxxz_0[j] - tg_zzz_xxxxxxxz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzz_xxxxxxx_1[j];

                    tg_zzzzz_xxxxxxyy_0[j] = pb_z * tg_zzzz_xxxxxxyy_0[j] + fr * tg_zzzz_xxxxxxyy_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxxxxyy_0[j] - tg_zzz_xxxxxxyy_1[j] * fl1_fza);

                    tg_zzzzz_xxxxxxyz_0[j] = pb_z * tg_zzzz_xxxxxxyz_0[j] + fr * tg_zzzz_xxxxxxyz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxxxxyz_0[j] - tg_zzz_xxxxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzz_xxxxxxy_1[j];

                    tg_zzzzz_xxxxxxzz_0[j] = pb_z * tg_zzzz_xxxxxxzz_0[j] + fr * tg_zzzz_xxxxxxzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxxxxzz_0[j] - tg_zzz_xxxxxxzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzz_xxxxxxz_1[j];

                    tg_zzzzz_xxxxxyyy_0[j] = pb_z * tg_zzzz_xxxxxyyy_0[j] + fr * tg_zzzz_xxxxxyyy_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxxxyyy_0[j] - tg_zzz_xxxxxyyy_1[j] * fl1_fza);

                    tg_zzzzz_xxxxxyyz_0[j] = pb_z * tg_zzzz_xxxxxyyz_0[j] + fr * tg_zzzz_xxxxxyyz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxxxyyz_0[j] - tg_zzz_xxxxxyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzz_xxxxxyy_1[j];

                    tg_zzzzz_xxxxxyzz_0[j] = pb_z * tg_zzzz_xxxxxyzz_0[j] + fr * tg_zzzz_xxxxxyzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxxxyzz_0[j] - tg_zzz_xxxxxyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzz_xxxxxyz_1[j];

                    tg_zzzzz_xxxxxzzz_0[j] = pb_z * tg_zzzz_xxxxxzzz_0[j] + fr * tg_zzzz_xxxxxzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxxxzzz_0[j] - tg_zzz_xxxxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzz_xxxxxzz_1[j];

                    tg_zzzzz_xxxxyyyy_0[j] = pb_z * tg_zzzz_xxxxyyyy_0[j] + fr * tg_zzzz_xxxxyyyy_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxxyyyy_0[j] - tg_zzz_xxxxyyyy_1[j] * fl1_fza);

                    tg_zzzzz_xxxxyyyz_0[j] = pb_z * tg_zzzz_xxxxyyyz_0[j] + fr * tg_zzzz_xxxxyyyz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxxyyyz_0[j] - tg_zzz_xxxxyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzz_xxxxyyy_1[j];

                    tg_zzzzz_xxxxyyzz_0[j] = pb_z * tg_zzzz_xxxxyyzz_0[j] + fr * tg_zzzz_xxxxyyzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxxyyzz_0[j] - tg_zzz_xxxxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzz_xxxxyyz_1[j];

                    tg_zzzzz_xxxxyzzz_0[j] = pb_z * tg_zzzz_xxxxyzzz_0[j] + fr * tg_zzzz_xxxxyzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxxyzzz_0[j] - tg_zzz_xxxxyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzz_xxxxyzz_1[j];

                    tg_zzzzz_xxxxzzzz_0[j] = pb_z * tg_zzzz_xxxxzzzz_0[j] + fr * tg_zzzz_xxxxzzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxxzzzz_0[j] - tg_zzz_xxxxzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzzz_xxxxzzz_1[j];

                    tg_zzzzz_xxxyyyyy_0[j] = pb_z * tg_zzzz_xxxyyyyy_0[j] + fr * tg_zzzz_xxxyyyyy_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxyyyyy_0[j] - tg_zzz_xxxyyyyy_1[j] * fl1_fza);

                    tg_zzzzz_xxxyyyyz_0[j] = pb_z * tg_zzzz_xxxyyyyz_0[j] + fr * tg_zzzz_xxxyyyyz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxyyyyz_0[j] - tg_zzz_xxxyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzz_xxxyyyy_1[j];

                    tg_zzzzz_xxxyyyzz_0[j] = pb_z * tg_zzzz_xxxyyyzz_0[j] + fr * tg_zzzz_xxxyyyzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxyyyzz_0[j] - tg_zzz_xxxyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzz_xxxyyyz_1[j];

                    tg_zzzzz_xxxyyzzz_0[j] = pb_z * tg_zzzz_xxxyyzzz_0[j] + fr * tg_zzzz_xxxyyzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxyyzzz_0[j] - tg_zzz_xxxyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzz_xxxyyzz_1[j];

                    tg_zzzzz_xxxyzzzz_0[j] = pb_z * tg_zzzz_xxxyzzzz_0[j] + fr * tg_zzzz_xxxyzzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxyzzzz_0[j] - tg_zzz_xxxyzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzzz_xxxyzzz_1[j];

                    tg_zzzzz_xxxzzzzz_0[j] = pb_z * tg_zzzz_xxxzzzzz_0[j] + fr * tg_zzzz_xxxzzzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxzzzzz_0[j] - tg_zzz_xxxzzzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_zzzz_xxxzzzz_1[j];

                    tg_zzzzz_xxyyyyyy_0[j] = pb_z * tg_zzzz_xxyyyyyy_0[j] + fr * tg_zzzz_xxyyyyyy_1[j] + 2.0 * fl1_fx * (tg_zzz_xxyyyyyy_0[j] - tg_zzz_xxyyyyyy_1[j] * fl1_fza);

                    tg_zzzzz_xxyyyyyz_0[j] = pb_z * tg_zzzz_xxyyyyyz_0[j] + fr * tg_zzzz_xxyyyyyz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxyyyyyz_0[j] - tg_zzz_xxyyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzz_xxyyyyy_1[j];

                    tg_zzzzz_xxyyyyzz_0[j] = pb_z * tg_zzzz_xxyyyyzz_0[j] + fr * tg_zzzz_xxyyyyzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxyyyyzz_0[j] - tg_zzz_xxyyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzz_xxyyyyz_1[j];

                    tg_zzzzz_xxyyyzzz_0[j] = pb_z * tg_zzzz_xxyyyzzz_0[j] + fr * tg_zzzz_xxyyyzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxyyyzzz_0[j] - tg_zzz_xxyyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzz_xxyyyzz_1[j];

                    tg_zzzzz_xxyyzzzz_0[j] = pb_z * tg_zzzz_xxyyzzzz_0[j] + fr * tg_zzzz_xxyyzzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxyyzzzz_0[j] - tg_zzz_xxyyzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzzz_xxyyzzz_1[j];

                    tg_zzzzz_xxyzzzzz_0[j] = pb_z * tg_zzzz_xxyzzzzz_0[j] + fr * tg_zzzz_xxyzzzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxyzzzzz_0[j] - tg_zzz_xxyzzzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_zzzz_xxyzzzz_1[j];

                    tg_zzzzz_xxzzzzzz_0[j] = pb_z * tg_zzzz_xxzzzzzz_0[j] + fr * tg_zzzz_xxzzzzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxzzzzzz_0[j] - tg_zzz_xxzzzzzz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_zzzz_xxzzzzz_1[j];

                    tg_zzzzz_xyyyyyyy_0[j] = pb_z * tg_zzzz_xyyyyyyy_0[j] + fr * tg_zzzz_xyyyyyyy_1[j] + 2.0 * fl1_fx * (tg_zzz_xyyyyyyy_0[j] - tg_zzz_xyyyyyyy_1[j] * fl1_fza);

                    tg_zzzzz_xyyyyyyz_0[j] = pb_z * tg_zzzz_xyyyyyyz_0[j] + fr * tg_zzzz_xyyyyyyz_1[j] + 2.0 * fl1_fx * (tg_zzz_xyyyyyyz_0[j] - tg_zzz_xyyyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzz_xyyyyyy_1[j];

                    tg_zzzzz_xyyyyyzz_0[j] = pb_z * tg_zzzz_xyyyyyzz_0[j] + fr * tg_zzzz_xyyyyyzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xyyyyyzz_0[j] - tg_zzz_xyyyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzz_xyyyyyz_1[j];

                    tg_zzzzz_xyyyyzzz_0[j] = pb_z * tg_zzzz_xyyyyzzz_0[j] + fr * tg_zzzz_xyyyyzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xyyyyzzz_0[j] - tg_zzz_xyyyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzz_xyyyyzz_1[j];

                    tg_zzzzz_xyyyzzzz_0[j] = pb_z * tg_zzzz_xyyyzzzz_0[j] + fr * tg_zzzz_xyyyzzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xyyyzzzz_0[j] - tg_zzz_xyyyzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzzz_xyyyzzz_1[j];

                    tg_zzzzz_xyyzzzzz_0[j] = pb_z * tg_zzzz_xyyzzzzz_0[j] + fr * tg_zzzz_xyyzzzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xyyzzzzz_0[j] - tg_zzz_xyyzzzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_zzzz_xyyzzzz_1[j];

                    tg_zzzzz_xyzzzzzz_0[j] = pb_z * tg_zzzz_xyzzzzzz_0[j] + fr * tg_zzzz_xyzzzzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xyzzzzzz_0[j] - tg_zzz_xyzzzzzz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_zzzz_xyzzzzz_1[j];

                    tg_zzzzz_xzzzzzzz_0[j] = pb_z * tg_zzzz_xzzzzzzz_0[j] + fr * tg_zzzz_xzzzzzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xzzzzzzz_0[j] - tg_zzz_xzzzzzzz_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_zzzz_xzzzzzz_1[j];

                    tg_zzzzz_yyyyyyyy_0[j] = pb_z * tg_zzzz_yyyyyyyy_0[j] + fr * tg_zzzz_yyyyyyyy_1[j] + 2.0 * fl1_fx * (tg_zzz_yyyyyyyy_0[j] - tg_zzz_yyyyyyyy_1[j] * fl1_fza);

                    tg_zzzzz_yyyyyyyz_0[j] = pb_z * tg_zzzz_yyyyyyyz_0[j] + fr * tg_zzzz_yyyyyyyz_1[j] + 2.0 * fl1_fx * (tg_zzz_yyyyyyyz_0[j] - tg_zzz_yyyyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzz_yyyyyyy_1[j];

                    tg_zzzzz_yyyyyyzz_0[j] = pb_z * tg_zzzz_yyyyyyzz_0[j] + fr * tg_zzzz_yyyyyyzz_1[j] + 2.0 * fl1_fx * (tg_zzz_yyyyyyzz_0[j] - tg_zzz_yyyyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzz_yyyyyyz_1[j];

                    tg_zzzzz_yyyyyzzz_0[j] = pb_z * tg_zzzz_yyyyyzzz_0[j] + fr * tg_zzzz_yyyyyzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_yyyyyzzz_0[j] - tg_zzz_yyyyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzz_yyyyyzz_1[j];

                    tg_zzzzz_yyyyzzzz_0[j] = pb_z * tg_zzzz_yyyyzzzz_0[j] + fr * tg_zzzz_yyyyzzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_yyyyzzzz_0[j] - tg_zzz_yyyyzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzzz_yyyyzzz_1[j];

                    tg_zzzzz_yyyzzzzz_0[j] = pb_z * tg_zzzz_yyyzzzzz_0[j] + fr * tg_zzzz_yyyzzzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_yyyzzzzz_0[j] - tg_zzz_yyyzzzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_zzzz_yyyzzzz_1[j];

                    tg_zzzzz_yyzzzzzz_0[j] = pb_z * tg_zzzz_yyzzzzzz_0[j] + fr * tg_zzzz_yyzzzzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_yyzzzzzz_0[j] - tg_zzz_yyzzzzzz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_zzzz_yyzzzzz_1[j];

                    tg_zzzzz_yzzzzzzz_0[j] = pb_z * tg_zzzz_yzzzzzzz_0[j] + fr * tg_zzzz_yzzzzzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_yzzzzzzz_0[j] - tg_zzz_yzzzzzzz_1[j] * fl1_fza) + 3.5 * fl1_fxn * tg_zzzz_yzzzzzz_1[j];

                    tg_zzzzz_zzzzzzzz_0[j] = pb_z * tg_zzzz_zzzzzzzz_0[j] + fr * tg_zzzz_zzzzzzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_zzzzzzzz_0[j] - tg_zzz_zzzzzzzz_1[j] * fl1_fza) + 4.0 * fl1_fxn * tg_zzzz_zzzzzzz_1[j];
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

